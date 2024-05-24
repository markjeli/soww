#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <sys/time.h>
#include "numgen.c"

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

  for (long i = 0; i < inputArgument; i++)
  {
    printf("%ld\n", numbers[i]);
  }

  // synchronize/finalize your CUDA computations

  gettimeofday(&ins__tstop, NULL);
  ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
}
