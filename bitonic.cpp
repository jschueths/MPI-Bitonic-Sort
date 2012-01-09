// Joshua McCarville-Schueths
// bitonic.cpp
// 
// This program is an implementation of the bitonic sort.
//

#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <ctime>
#include <cstdlib>

int compare (const void * a, const void * b);
void compareLow(int bit, unsigned int * list, unsigned int list_size);
void compareHigh(int bit, unsigned int list[], unsigned int list_size);
void localSort(unsigned int list[], unsigned int list_size);

using namespace std;

int rank, size;
int main(int argc, char * argv[])
{
  bool random = false;
  unsigned int data_size, elements;
  unsigned int * data;
  unsigned int * file_buffer = new unsigned int[400000];
  ifstream inFile;
  double sTime, eTime;
  bool skip = false;
  // Just get the initialization of the program going.
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  srand(time(NULL));
  if(argc < 2)
  {
    cout << "Not enough parameters" << endl;
    return 0;
  }
  else if(argc == 2)
  {
    random = true;
    elements = atoi(argv[1]);
  }
  else if(argc == 3)
  {
    elements = atoi(argv[2]);
    size = 1;
  }
  
  data_size = elements / size;
  data = new unsigned int[data_size];
  // Still need to read in the data set, break it up, and then sort the original sublist at each node.
  if(random || rank > 1)
  {
    for(unsigned int i = 0; i < data_size; i++)
    {
      data[i] = rand() % (2 * elements);
    }
  }
  else
  {
    if(rank == 0)
    {
      inFile.open(argv[1]);
      for(unsigned int j = 0; j < 400000; j++)
      {
        inFile >> file_buffer[j];
      }
      for(unsigned int i = 0; i < data_size; i++)
      {
        data[i] = file_buffer[i];
      }
    }
    else
    {
      skip = true;
      for(unsigned int i = 200000; i < 400000; i++)
      {
        data[i - 200000] = file_buffer[i];
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  sTime = MPI_Wtime();
  // Sort the local data set
  localSort(data, data_size);
  // Parallel portion of bitonic sort.
  int d = (int) log2(size);
  for(int i = 0; i < d && !skip; i++)
  {
    int window_id = rank; 
    window_id = window_id>>(i + 1);
    for(int j = i; j >= 0; j--)
    {
      int temp_id = rank;
      temp_id = temp_id>>j;
      if((window_id % 2 == 0 && temp_id % 2 == 0) || (window_id % 2 != 0 && temp_id % 2 != 0))
        compareLow(j, data, data_size);
      else
        compareHigh(j, data, data_size);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  eTime = MPI_Wtime();
  if(rank == 0)
    cout << "Time: " << eTime - sTime << endl;
  // Outputs the results from the input.400K iff 2 processors ran it.
  // I covered the case when all the processes ran, but it's just easier to output this way.
  int count = 0;
  if(count == rank && !random)
  {
    for(unsigned int i = 100000; i < 100010; i++)
    {
      cout << data[i] << " ";
      if(count == 100005)
        cout << endl;
    }
    cout << endl;
    for(unsigned int i = 200000; i < 200010; i++)
    {
      cout << data[i] << " ";
      if(count == 200005)
        cout << endl;
    }
    cout << endl;
  }

  delete [] data;
  delete [] file_buffer;
  MPI_Finalize();
  return 0;
}

void localSort(unsigned int list[], unsigned int list_size)
{
  qsort(list, list_size, sizeof(int), compare);
  return;
}

void compareLow(int bit, unsigned int list[], unsigned int list_size)
{
  // Figure out your communication partner.
  unsigned int index = list_size - 1;
  unsigned int temp;
  unsigned int send_count = 0;
  int recv_count;
  unsigned int min;
  bool swap = true;
  unsigned int * send_buffer = new unsigned int[list_size + 1];
  unsigned int * recv_buffer = new unsigned int[list_size + 1];
  int mask = 1<<bit;
  int comm_partner = rank ^ mask;
  MPI_Send(&list[index], 1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD);
  MPI_Recv(&min, 1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if(min >= list[index])
    return;
  for(unsigned int i = 0; i < list_size; i++)
  {
    if(list[i] > min)
    {
      send_buffer[send_count + 1] = list[i];
      send_count++;
    }
    else
    {
      break;
    }
  }
  send_buffer[0] = send_count;
  MPI_Send(send_buffer, send_count, MPI_INT, comm_partner, 0, MPI_COMM_WORLD);
  MPI_Recv(recv_buffer, list_size, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  recv_count = recv_buffer[0];
  for(unsigned int i = 1; i < recv_count + 1; i++)
  {
    if(list[index] < recv_buffer[i])
    {
      list[index] = recv_buffer[i];
      localSort(list, list_size);
    }
    else
    {
      break;
    }
  }
  delete [] send_buffer;
  delete [] recv_buffer;
  return;
}

void compareHigh(int bit, unsigned int list[], unsigned int list_size)
{
  // Figure out your communication partner.
  unsigned int temp;
  bool swap = true;
  unsigned int * send_buffer = new unsigned int[list_size + 1];
  unsigned int * recv_buffer = new unsigned int[list_size + 1];
  unsigned int send_count = 0;
  int recv_count;
  unsigned int max;
  int mask = 1<<bit;
  int comm_partner = rank ^ mask;
  MPI_Recv(&max, 1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Send(&list[0], 1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD);
  if(list[0] >= max)
    return;
  for(unsigned int i = 0; i < list_size; i++)
  {
    if(list[i] < max)
    {
      send_buffer[send_count + 1] = list[i];
      send_count++;
    }
    else
    {
      break;
    }
  }
  send_buffer[0] = send_count;
  MPI_Recv(recv_buffer, list_size, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  recv_count = recv_buffer[0];
  MPI_Send(send_buffer, send_count, MPI_INT, comm_partner, 0, MPI_COMM_WORLD);
  for(unsigned int i = 1; i < recv_count + 1; i++)
  {
    if(recv_buffer[i] > list[0])
    {
      list[0] = recv_buffer[i];
      localSort(list, list_size);
    }
    else
    {
      break;
    }
  }
  delete [] send_buffer;
  delete [] recv_buffer;
  return;
}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

