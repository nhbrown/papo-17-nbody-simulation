#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "plummer.h"
#include "output.h"
#include "hermite.h"

int main(int argc, const char *argv[])
{
  unsigned long seed; /* seed for Mersenne-Twister. */
  
  int N; /* amount of particles */
  int DIM; /* dimensions for space */

  double dt; /* timestep */
  double end_time; /* end of simulation */
  
  /* 
  M defines the total mass of the cluster.
  R defines the dimensions of the cluster.
  G defines the gravitational constant.
  */
  static const double M = 1.0;
  static const double R = 1.0;
  static const double G = 1.0;
  
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
  
  createNames();
  
  double *mass;
  
  mass = malloc(N * sizeof(double));
  if(mass == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
  
  double complex **pos;
  double complex **vel;

  double complex **acc;
  double complex **jerk;
  
  pos = malloc(N * sizeof(double complex *));
  vel = malloc(N * sizeof(double complex *));
  acc = malloc(N * sizeof(double complex *));
  jerk = malloc(N * sizeof(double complex *));
  
  for(int i = 0; i < N; ++i)
  {
    pos[i] = malloc(DIM * sizeof(double complex));
    vel[i] = malloc(DIM * sizeof(double complex));
    acc[i] = malloc(DIM * sizeof(double complex));
    jerk[i] = malloc(DIM * sizeof(double complex));
    
    if(pos[i] == NULL || vel[i] == NULL || acc[i] == NULL || jerk[i] == NULL)
    {
      fprintf(stderr, "Out of memory!\n");
      exit(0);
    }
  }
  
  printLog(seed, N, M, R, G, dt, end_time);
  
  startPlummer(seed, N, mass, pos, vel, M, R);
  
  printInitialConditions(N, mass, pos, vel);
  
  startHermite(N, dt, end_time, mass, pos, vel, acc, jerk);
  
  free(mass);
  
  for(int i = 0; i < N; ++i)
  {
    free(pos[i]);
    free(vel[i]);
    free(acc[i]);
    free(jerk[i]);
  }
  
  free(pos);
  free(vel);
  free(acc);
  free(jerk);
  
  return 0;
}
