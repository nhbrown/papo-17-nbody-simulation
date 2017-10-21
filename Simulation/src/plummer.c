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
#include "plummer.h"

/* factor for scaling to standard units (Heggie units) */
static const double scale = 16.0 / (3.0 * 3.14159265359);

/* returns pseudo-random number within specified range */
double rrand(double low, double high)
{
  return low + genrand_real1() * (high - low);
}

/* implementation of the Plummer density profile (Plummer Model) */
void plummer(int N, double *mass, double complex *pos, double complex *vel, int i, int mi, double M, double R)
{  
  mass[mi] = M / N; /* mass equilibrium */
  
  double complex radius = R / csqrt((cpow(genrand_real1(), (-2.0/3.0))) - 1.0); /* inverted cumulative mass distribution */
  double complex theta = cacos(rrand(-1.0, 1.0)); /* Polar Angle */
  double complex phi = rrand(0.0, (2 * 3.14159265359)); /* Azimuthal Angle */
  
  /* conversion from radial to cartesian coordinates */
  pos[i] = (radius * csin(theta) * ccos(phi)) / scale; 
  pos[i + 1] = (radius * csin(theta) * csin(phi)) / scale;
  pos[i + 2] = (radius * ccos(theta)) / scale;
  
  double x = 0.0;
  double y = 0.1;
  
  /* Neumann's rejection technique */
  while(y > (x * x * (pow((1.0 - x * x), 3.5))))
  {
    x = rrand(0.0, 1.0);
    y = rrand(0.0, 0.1);
  }
  
  /* distribution function */
  double complex velocity = x * csqrt(2.0) * cpow((1.0 + radius * radius), -0.25);
  theta = cacos(rrand(-1.0, 1.0));
  phi = rrand(0.0, (2 * 3.14159265359));
  
  /* conversion */
  vel[i] = (velocity * csin(theta) * ccos(phi)) * csqrt(scale);
  vel[i + 1] = (velocity * csin(theta) * csin(phi)) * csqrt(scale);
  vel[i + 2] = (velocity * ccos(theta)) * csqrt(scale);
}

/* adjusting positions and velocities of all particles towards center of mass */
void center_of_mass_adjustment(int N, int DIM, double *mass, double complex *pos, double complex *vel)
{
  double complex pos_center[3] = {0, 0, 0}; /* position of center of mass */
  double complex vel_center[3] = {0, 0, 0}; /* velocity of center of mass */
  
  /* measuring position and velocity of center of mass */
  for(int i = 0, mi = 0; i < (N * DIM); i += 3, ++mi) 
  {
    for(int j = 0; j < DIM; ++j)
    {
      pos_center[j] += pos[i + j] * mass[mi];
      vel_center[j] += vel[i + j] * mass[mi];
    }
  }
  
  /* subtracting position and velocity of center of mass from each particle */
  for(int k = 0; k < (N * DIM); k += 3) 
  {
    for(int l = 0; l < DIM; ++l)
    {
      pos[k + l] -= pos_center[l];
      vel[k + l] -= vel_center[l];
    }
  }
}

/* starts the Plummer Model routine to generate initial conditions for our cluster */
void startPlummer(unsigned long seed, int N, int DIM, double *mass, double complex *pos, double complex *vel, double M, double R)
{
  init_genrand(seed); /* seeding Mersenne Twister */

  for(int i = 0, mi = 0; i < (N * DIM); i += 3, ++mi)
  {
    plummer(N, mass, pos, vel, i, mi, M, R); /* generate masses, positions and velocities for specified amount of particles */
  }
  
  center_of_mass_adjustment(N, DIM, mass, pos, vel); /* adjust center of mass for all particles */
}
