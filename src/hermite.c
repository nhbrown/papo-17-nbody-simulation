/*
    The following source-code is an implementation of the fourth order 
    iterated time-symmetric Hermite integrator as described by Kokubo, 
    Yoshinaga & Makino, 1998.
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
#include <complex.h>
#include <string.h>
#include "hermite.h"
#include "output.h"

int N; /* amount of particles */
int DIM; /* dimensions */
double dt; /* timestep */

void acc_jerk(double *mass, double complex **pos, double complex **vel, 
              double complex **acc, double complex **jerk)
{ 
  for(int i = 0; i < N; ++i)
  {
    for(int k = 0; k < 3; ++k)
    {
      acc[i][k] = 0;
      jerk[i][k] = 0;
    }
  }
  
 for(int i = 0; i < N; ++i)
 {
   for(int j = i + 1; j < N; ++j) /* only loops over half of the particles */
   {
     double complex rji[DIM]; /* position vector from particle i to j */
     double complex vji[DIM]; /* velocity vector from particle i to j */
     
     double complex r2; /* rij^2 */
     double complex v2; /* vij^2 */
     double complex rv; /* rij*vij */
     
     for(int k = 0; k < DIM; ++k)
     {
        rji[k] = pos[j][k] - pos[i][k];
        vji[k] = vel[j][k] - vel[i][k];
       
        r2 += rji[k] * rji[k];
        v2 += vji[k] * vji[k];
        rv += rji[k] * vji[k];
     }
     
     rv /= r2;
     double complex r = csqrt(r2); /* absolute value of rij */
     double complex r3 = r * r2;
     
     double complex da[DIM];
     double complex dj[DIM];
     
     for (int k = 0; k < DIM ; k++)
     {
       da[k] = rji[k] / r3;
       dj[k] = (vji[k] - 3 * rv * rji[k]) / r3;
       
       acc[i][k] += mass[j] * da[k];
       acc[j][k] -= mass[i] * da[k];
       
       jerk[i][k] += mass[j] * dj[k];                
       jerk[j][k] -= mass[i] * dj[k];  
     }
   }
 }
}

void hermite(double *mass, double complex **pos, double complex **vel, 
              double complex **acc, double complex **jerk)
{
  double complex old_pos[N][DIM];
  double complex old_vel[N][DIM];  
  double complex old_acc[N][DIM];
  double complex old_jerk[N][DIM];
  
  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < DIM; ++j)
    {
      old_pos[i][j] = pos[i][j];
      old_vel[i][j] = vel[i][j];
      old_acc[i][j] = acc[i][j];
      old_jerk[i][j] = jerk[i][j];
    }
  }
  
  /* prediction for all particles */
  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < DIM; ++j)
    {
      pos[i][j] += vel[i][j] * dt + acc[i][j] * ((dt * dt)/2) + jerk[i][j] * ((dt * dt * dt)/6);
      vel[i][j] += acc[i][j] * dt + jerk[i][j] * ((dt * dt)/2);
    }
  }
  
  acc_jerk(mass, pos, vel, acc, jerk);
  
  /* correction in reversed order of computation */
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < DIM; ++j)
    {
      vel[i][j] = old_vel[i][j] + (old_acc[i][j] + acc[i][j]) * (dt/2) + (old_jerk[i][j] - jerk[i][j]) * ((dt * dt)/12);       
      pos[i][j] = old_pos[i][j] + (old_vel[i][j] + vel[i][j]) * (dt/2) + (old_acc[i][j] - acc[i][j]) * ((dt * dt)/12);
    }
  }
}

void startHermite(int particles, double timestep, double end, double *mass, double complex **pos, double complex **vel, 
                  double complex **acc, double complex **jerk)
{
  double time = 0.0;
  DIM = 3;
  
  N = particles;
  dt = timestep;
  
  acc_jerk(mass, pos, vel, acc, jerk);
  
  int iterations = 0;
  while(time < end_time)
  {
    hermite(mass, pos, vel, acc, jerk);
    time += dt;
    ++iterations;
    printIteration(mass, pos, vel, iterations);
  }
}
