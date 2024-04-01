#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include "numgen.c"

#define RANGESIZE 1
#define DATA 0
#define RESULT 1
#define FINISH 2

int isPrime(int n)
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

unsigned long CheckHowManyPrimes(unsigned long int numbers[], int size)
{
  unsigned long count = 0;
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

  // program input argument
  long inputArgument = ins__args.arg;

  struct timeval ins__tstart, ins__tstop;

  int myrank, nproc;
  unsigned long int *numbers;

  int a = 0, b = inputArgument;
  int range[2];
  int rangeToSend = 0;
  unsigned long result = 0, resulttemp;
  MPI_Status status;

  MPI_Init(&argc, &argv);

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

  if (((b - a) / RANGESIZE) < 2 * (nproc - 1))
  {
    printf("More subranges needed");
    MPI_Finalize();
    return -1;
  }

  // run your computations here (including MPI communication)
  // now the master will distribute the data and slave processes will perform computations
  if (myrank == 0)
  {                                                                                                        // Master
    unsigned long int *numbersToSend = (unsigned long int *)malloc(RANGESIZE * sizeof(unsigned long int)); // TODO: Might be problematic (BUFFER)
    range[0] = a;
    for (int i = 1; i < nproc; i++)
    {
      range[1] = range[0] + RANGESIZE;

      for (int j = range[0]; j < range[1]; j++)
      {
        numbersToSend[j - range[0]] = numbers[j];
      }

      rangeToSend = range[1] - range[0];
      MPI_Send(numbersToSend, rangeToSend, MPI_UNSIGNED_LONG, i, DATA, MPI_COMM_WORLD);

      range[0] = range[1];
    }

    do
    {

      // distribute remaining subranges to the processes which have completed their parts

      MPI_Recv(&resulttemp, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &status);
      result += resulttemp;

      range[1] = range[0] + RANGESIZE;
      if (range[1] > b)
        range[1] = b;

      for (int j = range[0]; j < range[1]; j++)
      {
        numbersToSend[j - range[0]] = numbers[j];
      }

      rangeToSend = range[1] - range[0];
      MPI_Send(numbersToSend, rangeToSend, MPI_UNSIGNED_LONG, status.MPI_SOURCE, DATA, MPI_COMM_WORLD);
      range[0] = range[1];

    } while (range[1] < b);

    for (int i = 0; i < nproc - 1; i++)
    {

      MPI_Recv(&resulttemp, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &status);

      result += resulttemp;
    }

    for (int i = 1; i < nproc; i++)
    {

      MPI_Send(NULL, 0, MPI_UNSIGNED_LONG, i, FINISH, MPI_COMM_WORLD);
    }

    printf("\nMASTER: there are %ld primes\n", result);
    free(numbersToSend);
  }
  else
  { // SLAVE
    int count;
    do
    {

      MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      if (status.MPI_TAG == DATA)
      {
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
        unsigned long int *receivedNumbers = (unsigned long int *)malloc(count * sizeof(unsigned long int)); // TODO: Might be problematic (BUFFER)

        MPI_Recv(receivedNumbers, count, MPI_UNSIGNED_LONG, 0, DATA, MPI_COMM_WORLD, &status);

        resulttemp = CheckHowManyPrimes(receivedNumbers, count);

        MPI_Send(&resulttemp, 1, MPI_UNSIGNED_LONG, 0, RESULT, MPI_COMM_WORLD);
        free(receivedNumbers);
      }

    } while (status.MPI_TAG != FINISH);
  }

  if (!myrank)
  {
    gettimeofday(&ins__tstop, NULL);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  }

  MPI_Finalize();
  return 0;
}
