#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "plummer.h"
#include "output.h"
/* #include "hermite.h" */

unsigned long seed; /* seed for Mersenne-Twister. */
int N; /* amount of particles */
int DIM; /* dimensions */

double dt; /* timestep */
double end_time; /* end of simulation */

void allocateArrays(double *mass, double complex **pos, double complex **vel, double complex **acc, double complex **jerk)
{
  mass = malloc(N * sizeof(double));
  if(mass == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
  
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
}

int main(int argc, const char *argv[])
{ 
  double *mass;

  double complex **pos;
  double complex **vel;

  double complex **acc;
  double complex **jerk;
  
  /* 
  M defines the total mass of the cluster.
  R defines the dimensions of the cluster.
  G defines the gravitational constant.
  */
  static const double M;
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
  
  allocateArrays(mass, pos, vel, acc, jerk);
  createNames();
  
  startPlummer(seed, N, mass, pos, vel, M, R, G);
  printInitialConditions(seed, N, M, R, G, dt, end_time, mass, pos, vel);
  
  /*
  startHermite(N, dt, end_time);
  */
  
  return 0;
}
