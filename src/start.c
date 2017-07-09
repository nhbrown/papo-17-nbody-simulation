#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "plummer.h"
#include "hermite.h"

int main(int argc, const char *argv[])
{
  unsigned long seed; /* seed for Mersenne-Twister. */
  int N; /* amount of particles */
  double dt; /* timestep */
  double end_time;
  
  char *foldername;
  
  switch(argc)
  {
    case 4 :
      seed = (unsigned long) time(NULL);
      N = atoi(argv[1]);
      dt = atof(argv[2]);
      end_time = atof(argv[3]);
      break;
      
    case 5 :
      seed = atol(argv[1]);
      N = atoi(argv[2]);
      dt = atof(argv[3]);
      end_time = atof(argv[4]);
      break;

    default : /* if less than 4 or more than 5 arguments are passed, the execution exits normally */
      printf("Invalid input for start.c!\n");
      exit(0);
  }
  
  foldername = startPlummer(seed, N);
  startHermite(N, dt, end_time, foldername);
  
  return 0;
}
