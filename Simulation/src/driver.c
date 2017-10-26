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

#include <complex.h>
#include "hermite.h"
#include "plummer.h"
#include "output.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* specifies the dimensions of space */
#define DIM   3

/* declaring function prototypes */
void callocArrays(int N, int DIM);
void freeArrays(void);

/* pointer to arrays holding mass, position, velocity, acceleration and jerk for all particles */
double *mass; 
double complex *pos, *vel, *acc, *jerk;

/*
 * Function:  main 
 * ====================
 *  entry point for simulation, computes user input, manages
 *  variables and controls computation
 *
 *  argc: amount of command line arguments
 *  argv: holds command line arguments
 *
 *  returns: zero
 * --------------------
 */
int main(int argc, const char *argv[])
{
  clock_t start = clock();
  
  int N = 0; /* amount of particles */
  unsigned long seed = 0; /* seed for Mersenne-Twister. */  

  double dt = 0.0; /* timestep */
  double end_time = 0.0; /* time where simulation ends */

  static const double M = 1.0; /* total mass cluster */
  static const double R = 1.0; /* radius of cluster */
  static const double G = 1.0; /* gravitational constant */
  
  /* computes command line arguments */
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

    default : /* in case more or less arguments are passed than allowed */
      printf("Invalid input for start.c!\n");
      exit(0);
  }
  
  /* check wether user input is allowed or not */
  if(N <= 0 || dt <= 0 || end_time <= 0)
  {
    fprintf(stderr, "Negative values are not allowed!\n");
    exit(0);
  }
  
  createNames(); /* provided by output.h */
  
  callocArrays(N, DIM);
  
  printLog(seed, N, M, R, G, dt, end_time); /* provided by output.h */
  
  startPlummer(seed, N, DIM, mass, pos, vel, M, R); /* provided by plummer.h */
  
  printInitialConditions(N, DIM, mass, pos, vel); /* provided by output.h */
  
  startHermite(N, DIM, dt, end_time, mass, pos, vel, acc, jerk); /* provided by hermite.h */
  
  freeArrays();
  
  /* calculate total cpu time in seconds and print it to default output */
  clock_t end = clock();
  double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  
  printf("CPU time used: %f", cpu_time);
  
  return 0;
}

/*
 * Function:  callocArrays 
 * ====================
 *  Allocates necessary memory for all arrays declared
 *  above and also initializes each index to zero.
 *
 *  N: amount of particles
 *  DIM: dimensions of space
 *
 *  returns: void
 * --------------------
 */
void callocArrays(int N, int DIM)
{
  mass = calloc(N, sizeof(double));
  pos = calloc((N * DIM), sizeof(double complex));
  vel = calloc((N * DIM), sizeof(double complex));
  acc = calloc((N * DIM), sizeof(double complex));
  jerk = calloc((N * DIM), sizeof(double complex));
  
  if(mass == NULL || pos == NULL || vel == NULL || acc == NULL || jerk == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
}

/*
 * Function:  freeArrays 
 * ====================
 *  Frees all previously allocated memory. 
 *
 *  returns: void
 * --------------------
 */
void freeArrays()
{
  free(mass);
  free(pos);
  free(vel);
  free(acc);
  free(jerk);
}
