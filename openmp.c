#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include "numgen.c"

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

  // set number of threads
  omp_set_num_threads(ins__args.n_thr);

  // program input argument
  long inputArgument = ins__args.arg;
  unsigned long int *numbers = (unsigned long int *)malloc(inputArgument * sizeof(unsigned long int));
  numgen(inputArgument, numbers);

  struct timeval ins__tstart, ins__tstop;
  gettimeofday(&ins__tstart, NULL);

  // run your computations here (including OpenMP stuff)
  long dataChunk = inputArgument / ins__args.n_thr;
  unsigned long count = 0;

#pragma omp parallel for reduction(+ : count)
  for (int i = 0; i < ins__args.n_thr; i++)
  {
    count += CheckHowManyPrimes(numbers + i * dataChunk, dataChunk);
  }

  printf("Number of primes: %lu\n", count);

  // synchronize/finalize your computations
  gettimeofday(&ins__tstop, NULL);
  ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
}
