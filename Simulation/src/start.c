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
void mallocArrays(int N);
void initializeArrays(int N);
void freeArrays(int N);

double *mass; /* holds masses for all particles */
  
double complex **pos; /* holds positions for all particles */
double complex **vel; /* hold velocities for all particles */

double complex **acc; /* holds acceleration for all particles */
double complex **jerk; /* holds jerk for all particles */

/* computes user input and starts the simulation */
int main(int argc, const char *argv[])
{
  int N = 0; /* amount of particles */
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
  
  mallocArrays(N); /* allocates space for arrays */
  
  initializeArrays(N); /* zero out all elements */
  
  printLog(seed, N, M, R, G, dt, end_time); /* creates and writes to the log file */
  
  startPlummer(seed, N, mass, pos, vel, M, R); /* starts the Plummer Model routine for initial conditions */
  
  printInitialConditions(N, mass, pos, vel); /* creates and writes to the initial conditions file */
    
  energy_diagnostics(N, 0, mass, pos, vel); /* calculate kinetic, potential and total energy of the cluster at start */
  
  startHermite(N, dt, end_time, mass, pos, vel, acc, jerk); /* starts the Hermite scheme for further computation */
    
  energy_diagnostics(N, 1, mass, pos, vel); /* calculate kinetic, potential and total energy of the cluster at end */
  
  freeArrays(N); /* frees allocated space of arrays after computation has finished */
  
  return 0;
}

/* allocates neccessary space for all arrays */
void mallocArrays(int N)
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
    pos[i] = malloc(3 * sizeof(double complex));
    vel[i] = malloc(3 * sizeof(double complex));
    acc[i] = malloc(3 * sizeof(double complex));
    jerk[i] = malloc(3 * sizeof(double complex));
    
    if(pos[i] == NULL || vel[i] == NULL || acc[i] == NULL || jerk[i] == NULL)
    {
      fprintf(stderr, "Out of memory!\n");
      exit(0);
    }
  }
}

/* zeroes out all elements of our arrays to get rid of any garbage values which might have been in memory */
void initializeArrays(int N)
{
  for(int i = 0; i < N; ++i)
  {
    mass[i] = 0.0;
    
    for(int j = 0; j < 3; ++j)
    {
      pos[i][j] = vel[i][j] = acc[i][j] = jerk[i][j] = 0.0;
    }
  }
}

/* frees the allocated space of all arrays */
void freeArrays(int N)
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
