#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include "numgen.c"

#define RANGESIZE 5
#define DATA 0
#define RESULT 1
#define FINISH 2

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

unsigned long CheckHowManyPrimes(unsigned long numbers[], int size)
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

  MPI_Request *requests;
  int requestcompleted;
  int myrank, nproc;
  unsigned long int *numbers;

  int a = 0, b = inputArgument;
  int range[2];
  unsigned long result = 0;
  unsigned long *resulttemp;
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
  { // Master
    requests = (MPI_Request *)malloc(3 * (nproc - 1) * sizeof(MPI_Request));
    if (!requests)
    {
      printf("\nNot enough memory");
      MPI_Finalize();
      return -1;
    }

    unsigned long int *numbersToSend = (unsigned long *)malloc(2 * RANGESIZE * (nproc - 1) * sizeof(unsigned long));
    if (!numbersToSend)
    {
      printf("\nNot enough memory");
      MPI_Finalize();
      return -1;
    }

    resulttemp = (unsigned long *)malloc((nproc - 1) * sizeof(unsigned long));
    if (!resulttemp)
    {
      printf("\nNot enough memory");
      MPI_Finalize();
      return -1;
    }

    // first distribute some numbers from ranges to all slaves
    range[0] = a;
    for (int i = 1; i < nproc; i++)
    {
      range[1] = range[0] + RANGESIZE;
      for (int j = range[0]; j < range[1]; j++)
      {
        numbersToSend[j - range[0]] = numbers[j];
      }

      MPI_Send(numbersToSend, RANGESIZE, MPI_UNSIGNED_LONG, i, DATA, MPI_COMM_WORLD);

      range[0] = range[1];
    }

    // the first proccount requests will be for receiving, the latter ones for sending
    for (int i = 0; i < 2 * (nproc - 1); i++)
      requests[i] = MPI_REQUEST_NULL; // none active at this point

    // start receiving for results from the slaves
    for (int i = 1; i < nproc; i++)
      MPI_Irecv(&(resulttemp[i - 1]), 1, MPI_UNSIGNED_LONG, i, RESULT, MPI_COMM_WORLD, &(requests[i - 1]));

    // start sending new data parts to the slaves
    for (int i = 1; i < nproc; i++)
    {
      range[1] = range[0] + RANGESIZE;
      for (int j = range[0]; j < range[1]; j++)
      {
        numbersToSend[RANGESIZE * (i - 1) + (j - range[0])] = numbers[j];
      }

      MPI_Isend(&numbersToSend[RANGESIZE * (i - 1)], RANGESIZE, MPI_UNSIGNED_LONG, i, DATA, MPI_COMM_WORLD, &(requests[nproc - 1 + i - 1]));
      range[0] = range[1];
    }

    while (range[1] < b)
    {
      MPI_Waitany(2 * (nproc - 1), requests, &requestcompleted, MPI_STATUS_IGNORE);

      // if it is a result then send new data to the process and add the result
      if (requestcompleted < (nproc - 1))
      {
        result += resulttemp[requestcompleted];

        // first check if the send has terminated
        MPI_Wait(&(requests[nproc - 1 + requestcompleted]), MPI_STATUS_IGNORE);

        // now send some new data portion to this process
        range[1] = range[0] + RANGESIZE;

        for (int j = range[0]; j < range[1]; j++)
        {
          if (j < b) {
            numbersToSend[RANGESIZE * requestcompleted + (j - range[0])] = numbers[j];
          }
          else 
          {
            numbersToSend[RANGESIZE * requestcompleted + (j - range[0])] = 0;
          }
        }

        
        MPI_Isend(&numbersToSend[RANGESIZE * requestcompleted], RANGESIZE, MPI_UNSIGNED_LONG, requestcompleted + 1, DATA, MPI_COMM_WORLD, &(requests[nproc - 1 + requestcompleted]));
        range[0] = range[1];

        // now issue a corresponding recv
        MPI_Irecv(&(resulttemp[requestcompleted]), 1, MPI_UNSIGNED_LONG, requestcompleted + 1, RESULT, MPI_COMM_WORLD, &(requests[requestcompleted]));
      }
    }

    // now send the FINISHING ranges to the slaves and shut down the slaves
    // TODO: Might need to change this part    
    range[0] = range[1];
    for (int i=0;i<RANGESIZE;i++) {
      numbersToSend[i] = 0;
    }
    // unsigned long ZEROARRAY[RANGESIZE] = {0};
    for (int i = 1; i < nproc; i++)
    {
      for (int j = 0; j < RANGESIZE; j++) {        
        numbersToSend[RANGESIZE*(i-1) + RANGESIZE*(nproc - 1) + j] = 0;
      }
      
      MPI_Isend(&numbersToSend[RANGESIZE*(i-1) + RANGESIZE*(nproc - 1)], RANGESIZE, MPI_UNSIGNED_LONG, i, DATA, MPI_COMM_WORLD, &(requests[2 * nproc - 3 + i])); // nproc - 1 + i - 1
    }

    MPI_Waitall(3 * nproc - 3, requests, MPI_STATUSES_IGNORE);

    // now simply add the results
    for (int i = 0; i < (nproc - 1); i++)
    {
      result += resulttemp[i];
    }

    // now receive results for the initial sends
    for (int i = 0; i < (nproc - 1); i++)
    {
      MPI_Recv(&(resulttemp[i]), 1, MPI_UNSIGNED_LONG, i + 1, RESULT, MPI_COMM_WORLD, &status);
      result += resulttemp[i];
    }

    // display the result
    printf("\nMASTER: there are %ld primes\n", result);
    free(numbersToSend);
    free(resulttemp);
    free(requests);
    free(numbers);
  }
  else
  { // SLAVE
    requests = (MPI_Request *)malloc(2 * sizeof(MPI_Request));
    if (!requests)
    {
      printf("\nNot enough memory");
      MPI_Finalize();
      return -1;
    }

    requests[0] = requests[1] = MPI_REQUEST_NULL;
    resulttemp = (unsigned long *)malloc(2 * sizeof(unsigned long));
    if (!resulttemp)
    {
      printf("\nNot enough memory");
      MPI_Finalize();
      return -1;
    }   

    // first receive the initial data    
    unsigned long *receivedNumbers = (unsigned long *)malloc(RANGESIZE * sizeof(unsigned long));
    unsigned long *numbersToRecieve = (unsigned long *)malloc(RANGESIZE * sizeof(unsigned long));
    MPI_Recv(receivedNumbers, RANGESIZE, MPI_UNSIGNED_LONG, 0, DATA, MPI_COMM_WORLD, &status);
    int isRunning = 1;

    while (isRunning)
    {
      isRunning = 0;
      // if there is some data to process
      // before computing the next part start receiving a new data part
      MPI_Irecv(numbersToRecieve, RANGESIZE, MPI_UNSIGNED_LONG, 0, DATA, MPI_COMM_WORLD, &(requests[0]));

      // Compute
      resulttemp[1] = CheckHowManyPrimes(receivedNumbers, RANGESIZE);

      // now finish receiving the new part and finish sending the previous results back to the master
      MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

      for (int j=0; j<RANGESIZE; j++) {
        receivedNumbers[j] = numbersToRecieve[j];
        if (numbersToRecieve[j] != 0) {
          isRunning = 1;
        }
      }
      resulttemp[0] = resulttemp[1];

      // and start sending the results back
      MPI_Isend(&resulttemp[0], 1, MPI_UNSIGNED_LONG, 0, RESULT, MPI_COMM_WORLD, &(requests[1]));
    }

    // now finish sending the last results to the master
	  MPI_Wait (&(requests[1]), MPI_STATUS_IGNORE);

    free(receivedNumbers);
    free(numbersToRecieve);
    free(resulttemp);
    free(requests);
  }

  if (!myrank)
  {
    gettimeofday(&ins__tstop, NULL);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  }

  MPI_Finalize();
  return 0;
}
