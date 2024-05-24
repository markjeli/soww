#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <sys/time.h>
#include "numgen.c"

__device__
unsigned long isPrime(unsigned long n)
{
  // Corner case
  if (n <= 1)
    return 0;

  // Check from 2 to n-1
  for (int i = 2; i <= sqrtf(static_cast<float>(n)); i++)
    if (n % i == 0)
      return 0;

  return 1;
}

__global__
void CheckHowManyPrimes(unsigned long numbers[], unsigned long long *primes, long size)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < size)
  {
    if (isPrime(numbers[index]))
    {
      atomicAdd(primes, 1);
    }
  }
}

int main(int argc, char **argv)
{

  Args ins__args;
  parseArgs(&ins__args, &argc, argv);

  // program input argument
  long inputArgument = ins__args.arg;
  unsigned long int *numbers = (unsigned long int *)malloc(inputArgument * sizeof(unsigned long int));
  numgen(inputArgument, numbers);

  struct timeval ins__tstart, ins__tstop;
  gettimeofday(&ins__tstart, NULL);

  // run your CUDA kernel(s) here
  unsigned long *d_numbers;
  unsigned long long *d_primes;
  unsigned long long primes = 0;

  cudaMalloc(&d_numbers, inputArgument * sizeof(unsigned long));
  cudaMalloc(&d_primes, sizeof(unsigned long long));

  cudaMemcpy(d_numbers, numbers, inputArgument * sizeof(unsigned long), cudaMemcpyHostToDevice);
  cudaMemcpy(d_primes, &primes, sizeof(unsigned long long), cudaMemcpyHostToDevice);

  int blockSize = 256;
  int gridSize = (inputArgument + blockSize - 1) / blockSize;

  CheckHowManyPrimes<<<gridSize, blockSize>>>(d_numbers, d_primes, inputArgument);

  cudaMemcpy(&primes, d_primes, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  printf("Number of primes: %lld\n", primes);

  cudaFree(d_numbers);
  cudaFree(d_primes);

  free(numbers);

  // synchronize/finalize your CUDA computations

  gettimeofday(&ins__tstop, NULL);
  ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
}
