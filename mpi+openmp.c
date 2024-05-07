#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "numgen.c"

#define RESULT 1

unsigned long isPrime(unsigned long n)
{
  // Corner case
  if (n <= 1)
    return 0;

  // Check from 2 to n-1
  for (int i = 2; i <= sqrt(n); i++)
    if (n % i == 0)
      return 0;

  return 1;
}

unsigned long CheckHowManyPrimesParallel(unsigned long numbers[], int size)
{
  unsigned long count = 0;
#pragma omp parallel for reduction(+ : count)
  for (int i = 0; i < size; i++)
  {
    if (isPrime(numbers[i]))
    {
      count++;
    }
  }
  return count;
}

int main(int argc, char **argv)
{

  Args ins__args;
  parseArgs(&ins__args, &argc, argv);

  // set number of threads
  omp_set_num_threads(ins__args.n_thr);

  // program input argument
  long inputArgument = ins__args.arg;

  struct timeval ins__tstart, ins__tstop;

  int threadsupport;
  int myrank, nproc;
  unsigned long int *numbers;
  // Initialize MPI with desired support for multithreading -- state your desired support level

  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadsupport);

  if (threadsupport < MPI_THREAD_FUNNELED)
  {
    printf("\nThe implementation does not support MPI_THREAD_FUNNELED, it supports level %d\n", threadsupport);
    MPI_Finalize();
    return -1;
  }

  // obtain my rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // and the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (!myrank)
  {
    gettimeofday(&ins__tstart, NULL);
    numbers = (unsigned long int *)malloc(inputArgument * sizeof(unsigned long int));
    numgen(inputArgument, numbers);
  }
  // run your computations here (including MPI communication and OpenMP stuff)
  if (!myrank)
  {
    long dataChunk = inputArgument / ins__args.n_thr;
    for (int i = 1; i < nproc; i++)
    {
      if (i == nproc - 1)
        MPI_Send(numbers + (i - 1) * dataChunk, inputArgument - (i - 1) * dataChunk, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      else
        MPI_Send(numbers + (i - 1) * dataChunk, dataChunk, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
    }

    unsigned long result = 0;
    unsigned long primes;
    for (int i = 1; i < nproc; i++)
    {
      MPI_Recv(&primes, 1, MPI_UNSIGNED_LONG, i, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      result += primes;
    }
    printf("Number of primes: %lu\n", result);
  }
  else
  {
    MPI_Status status;
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    int count;
    MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
    numbers = (unsigned long int *)malloc(count * sizeof(unsigned long int));
    MPI_Recv(numbers, count, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    unsigned long primes = CheckHowManyPrimesParallel(numbers, count);
    MPI_Send(&primes, 1, MPI_UNSIGNED_LONG, 0, RESULT, MPI_COMM_WORLD);
  }

  // synchronize/finalize your computations

  if (!myrank)
  {
    gettimeofday(&ins__tstop, NULL);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  }

  MPI_Finalize();
  free(numbers);
}
