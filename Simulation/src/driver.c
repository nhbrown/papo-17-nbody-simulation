/*    
    The following source code provides an entry point for the NBody routine,
    computing user input, initialising and maintaining all important 
    variables and containers, as well as orchestrating the execution.
    Copyright (C) 2017  Nicholas Lee Hickson-Brown, Michael Eidus
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include <mpi.h>
#include "plummer.h"
#include "output.h"
#include "hermite.h"

/* prototypes */
void mallocArrays(int N, int DIM);
void initializeArrays(int N, int DIM);
void freeArrays(void);

double *mass; /* holds masses of all particles */  
double complex *pos, *vel, *acc, *jerk; /* containers for positions, velocities, acceleration and jerks for all particles */

int main(int argc, const char *argv[])
{
  clock_t start = clock();
  
  int world_size = 0, world_rank = 0;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  int N = 0; /* amount of particles */
  int DIM = 3; /* dimensions */
  unsigned long seed = 0; /* seed for Mersenne-Twister. */  

  double dt = 0.0; /* timestep */
  double end_time = 0.0; /* end of simulation */

  static const double M = 1.0; /* total mass of the cluster */
  static const double R = 1.0; /* radius of cluster */
  static const double G = 1.0; /* gravitational constant */
  
  /* input computation */
  switch(argc)
  {
    case 4 : /* when no seed is specified by user */
      seed = (unsigned long) time(NULL);
      N = atoi(argv[1]);
      dt = atof(argv[2]);
      end_time = atof(argv[3]);
      break;
      
    case 5 : /* when seed is specified by user */
      seed = atol(argv[1]);
      N = atoi(argv[2]);
      dt = atof(argv[3]);
      end_time = atof(argv[4]);
      break;

    default : /* if less than 4 or more than 5 arguments are passed, the execution exits normally */
      printf("Invalid input for start.c!\n");
      exit(0);
  }
  
  createNames(); /* creates folder and names for files */
  
  mallocArrays(N, DIM); /* allocates space for arrays */
  
  initializeArrays(N, DIM); /* zero out all elements */
  
  if(world_rank == 0)
  {
    printLog(seed, N, M, R, G, dt, end_time); /* creates and writes to the log file */
  }
  
  startPlummer(seed, N, mass, pos, vel, M, R); /* starts the Plummer Model routine for initial conditions */
  
  if(world_rank == 0)
  {
    printInitialConditions(N, mass, pos, vel); /* creates and writes to the initial conditions file */
  }

  int proc_elem = (N / world_size) * DIM;
  startHermite(N, DIM, dt, end_time, mass, pos, vel, acc, jerk, proc_elem, world_rank); /* starts the Hermite scheme for further computation */
  
  freeArrays(); /* frees allocated space of arrays after computation has finished */
  
  if(world_rank == 0)
  {
    clock_t end = clock();
    double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  
    printf("CPU time used: %f", cpu_time); 
  }
  
  MPI_Finalize();
  
  return 0;
}

/* allocates neccessary space for all arrays */
void mallocArrays(int N, int DIM)
{
  mass = malloc(N * sizeof(double));
  pos = malloc((N * sizeof(double complex)) * DIM);
  vel = malloc((N * sizeof(double complex)) * DIM);
  acc = malloc((N * sizeof(double complex)) * DIM);
  jerk = malloc((N * sizeof(double complex)) * DIM);
  
  if(mass == NULL || pos == NULL || vel == NULL || acc == NULL || jerk == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
}

/* zeroes out all elements of our arrays to get rid of any garbage values which might have been in memory */
void initializeArrays(int N, int DIM)
{
  for(int i = 0; i < N; ++i)
  {
    mass[i] = 0.0;
    
    for(int j = 0; j < DIM; ++j)
    {
      pos[i + j] = vel[i + j] = acc[i + j] = jerk[i + j] = 0.0;
    }
  }
}

/* frees the allocated space of all arrays */
void freeArrays()
{
  free(mass);
  free(pos);
  free(vel);
  free(acc);
  free(jerk);
}
