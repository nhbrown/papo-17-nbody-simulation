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
#include "plummer.h"

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

/* factor for scaling to standard units */
static const double scale = 16.0 / (3.0 * M_PI);

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
  p.xpos = (radius * csin(theta) * ccos(phi)) / scale; 
  p.ypos = (radius * csin(theta) * csin(phi)) / scale;
  p.zpos = (radius * ccos(theta)) / scale;
  
  double x = 0.0;
  double y = 0.1;
  
  while(y > (x * x * (pow((1.0 - x * x), 3.5))))
  {
    x = frand(0.0, 1.0);
    y = frand(0.0, 0.1);
  }
  
  double complex velocity = x * csqrt(2.0) * cpow((1.0 + radius * radius), -0.25);
  theta = cacos(frand(-1.0, 1.0));
  phi = frand(0.0, (2 * M_PI));
  
  p.xvel = (velocity * csin(theta) * ccos(phi)) * csqrt(scale);
  p.yvel = (velocity * csin(theta) * csin(phi)) * csqrt(scale);
  p.zvel = (velocity * ccos(theta)) * csqrt(scale);
}

/*
Passing a seed as a parameter is optional, if no seed is passed seed is equal to UNIX-clock.
Specifying the amount of particles to generate is always necessary.
If user wishes to specify the seed, the order of arguments needs to be: <executable> seed amount
*/
char * startPlummer(unsigned long s, int amount, double timestep, double end_time)
{
  seed = s;
  N = amount;
  
  init_genrand(seed);
  
  generateOutput(timestep, end_time);
  
  return foldername;
}
