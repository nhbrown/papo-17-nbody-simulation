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
#include "plummer.h"
#include "output.h"
#include "hermite.h"
#include "ediag.h"

/* prototypes */
void mallocArrays(void);
void freeArrays(void);

/* for measuring cpu time */
clock_t start;
clock_t end;
double cpu_time;

int N; /* amount of particles */
int DIM = 3; /* dimensions */

double *mass; /* holds masses for all particles */
  
double complex **pos; /* holds positions for all particles */
double complex **vel; /* hold velocities for all particles */

double complex **acc; /* holds acceleration for all particles */
double complex **jerk; /* holds jerk for all particles */

/* allocates neccessary space for all arrays */
void mallocArrays()
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

/* frees the allocated space of all arrays */
void freeArrays()
{
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
}

int main(int argc, const char *argv[])
{
  start = clock();
  
  unsigned long seed; /* seed for Mersenne-Twister. */  

  double dt; /* timestep */
  double end_time; /* end of simulation */

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
  
  mallocArrays(); /* allocates space for arrays */
  
  printLog(seed, N, M, R, G, dt, end_time); /* creates and writes to the log file */
  
  startPlummer(seed, N, mass, pos, vel, M, R); /* starts the Plummer Model routine for initial conditions */
  
  printInitialConditions(N, mass, pos, vel); /* creates and writes to the initial conditions file */
    
  energy_diagnostics(N, DIM, 0, mass, pos, vel); /* calculate kinetic, potential and total energy of the cluster at start */
  
  startHermite(N, DIM, dt, end_time, mass, pos, vel, acc, jerk); /* starts the Hermite scheme for further computation */
    
  energy_diagnostics(N, DIM, 0, mass, pos, vel); /* calculate kinetic, potential and total energy of the cluster at end */
  
  freeArrays(); /* frees allocated space of arrays after computation has finished */
  
  end = clock();
  cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  
  printf("CPU time used: %f", cpu_time);
  
  return 0;
}
