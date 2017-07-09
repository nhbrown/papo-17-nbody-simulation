/*    
    The following source code is an implementation of the Plummer three-dimensional
    density profile for generating the inital conditions of a globular cluster,
    also known as an initial conditions generator for N-body simulations.
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
#include <math.h>
#include "mersenne.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

/* Seed for Mersenne-Twister. */
unsigned long seed;
/* Total number of particles to be generated. */
int N;

struct body
{
  double complex xpos;
  double complex ypos;
  double complex zpos;
  
  double mass;  
  
  double complex xvel;
  double complex yvel;
  double complex zvel;
};

/* single particle */
struct body p;

/* 
M defines the total mass of the cluster.
R defines the dimensions of the cluster.
G defines the gravitational constant.
*/
static const double M = 1.0;
static const double R = 1.0;
static const double G = 1.0;

double frand(double low, double high)
{
  return low + genrand_real1() * (high - low);
}

void plummer()
{  
  p.mass = M / N; /* mass equilibrium */
  
  double complex radius = R / csqrt((cpow(genrand_real1(), (-2.0/3.0))) - 1.0); /* inverted cumulative mass distribution */
  double complex theta = cacos(frand(-1.0, 1.0)); /* Polar Angle */
  double complex phi = frand(0.0, (2 * M_PI)); /* Azimuthal Angle */
  
  /* conversion from radial to cartesian coordinates */
  p.xpos = radius * csin(theta) * ccos(phi); 
  p.ypos = radius * csin(theta) * csin(phi);
  p.zpos = radius * ccos(theta);
  
  double x = 0.0;
  double y = 0.1;
  
  while(y > pow((x * x * (1.0 - x * x)), 3.5))
  {
    x = frand(0.0, 1.0);
    y = frand(0.0, 0.1);
  }
  
  double complex velocity = x * csqrt(2.0) * cpow((1.0 + radius * radius), -0.25);
  theta = cacos(frand(-1.0, 1.0));
  phi = frand(0.0, (2 * M_PI));
  
  p.xvel = velocity * csin(theta) * ccos(phi);
  p.yvel = velocity * csin(theta) * csin(phi);
  p.zvel = velocity * ccos(theta);
}

char foldername[40];
char logname[80];
char conditionsname[80];

void createNames()
{
  struct tm *sTm;

  time_t now = time(0);
  sTm = gmtime(&now);

  strftime (foldername, sizeof(foldername), "run_%Y_%m_%d_%H:%M:%S", sTm);
  strftime (logname, sizeof(logname), "run_%Y_%m_%d_%H:%M:%S/log_%Y_%m_%d_%H:%M:%S.txt", sTm);
  strftime (conditionsname, sizeof(conditionsname), "run_%Y_%m_%d_%H:%M:%S/initial_conditions.csv", sTm);

  struct stat st = {0};

  if (stat(foldername, &st) == -1)
  {
      mkdir(foldername, 0700);
  }
}

void generateOutput()
{ 
  createNames();
  
  FILE *log;
  log = fopen(logname, "w"); /* writes to new file log_<currentdate>.txt which holds important parameters */

  fprintf(log, "Seed used: %lu \nNumber of particles: %d \nTotal mass of cluster: %f \nDimensions of cluster: %f \nGravitational constant: %f", 
          seed, N, M, R, G);

  fclose(log);

  FILE *conditions;
  conditions = fopen(conditionsname, "w"); /* writes to new file initial_conditions.csv which holds positions */

  for(int j = 0; j < N; ++j)
  {
    plummer();
    fprintf(conditions, "%f, %f, %f, %f, %f, %f, %f\n", 
            creal(p.xpos), creal(p.ypos), creal(p.zpos), p.mass, creal(p.xvel), creal(p.yvel), creal(p.zvel));
  }

  fclose(conditions);
}

/*
Passing a seed as a parameter is optional, if no seed is passed seed is equal to UNIX-clock.
Specifying the amount of particles to generate is always necessary.
If user wishes to specify the seed, the order of arguments needs to be: <executable> seed amount
*/
int main(int argc, const char *argv[])
{
  switch(argc)
  {
    case 2 : /* if one argument is passed, it is assumed to be amount of particles */
      seed = (unsigned long)time(NULL);
      N = atoi(argv[1]);
      break;

    case 3 : /* if two arguments are passed, first one is assumed to be seed */
      seed = atol(argv[1]);
      N = atoi(argv[2]);
      break;

    default : /* if less than 1 or more than 2 arguments are passed, the executions exits */
      printf("Invalid input for plummer.c!\n");
      exit(0);
  }
  
  init_genrand(seed);
  
  generateOutput();
  
  return 0;
}
