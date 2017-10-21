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

#include <complex.h>
#include <string.h>
#include "hermite.h"
#include "output.h"
#include "ediag.h"

/* calculates acceleration and jerk for all particles */
void acc_jerk(int N, int DIM, double *mass, double complex *pos, double complex *vel, 
              double complex *acc, double complex *jerk)
{ 
  /* default values for acceleration and jerk */
  for(int i = 0; i < N; i += 3)
  {
    for(int j = 0; j < DIM; ++j)
    {
      acc[i + j] = jerk[i + j] = 0;
    }
  }
  
  for(int i = 0, mi = 0; i < N; i += 3, ++mi) /* loops over all particles */
  { 
    for(int j = i + 3, mj = 0; j < N; j += 3, ++mj) /* only loops over half of the particles because force acts equally on both particles (Newton) */
    {
      double complex rji[DIM], vji[DIM]; /* position vector from particle i to j */
      
      for(int k = 0; k < DIM; ++k)
      {
        rji[k] = vji[k] = 0.0;
      }
     
      double complex r2 = 0.0; /* rij^2 */
      double complex rv = 0.0; /* rij*vij */
     
      /* calculating position and velocity vectors */
      for(int k = 0; k < DIM; ++k)
      {
        rji[k] = pos[j + k] - pos[i + k];
        vji[k] = vel[j + k] - vel[i + k];
       
        r2 += rji[k] * rji[k];
        rv += rji[k] * vji[k];
      }
      
      double complex r3 = csqrt(r2) * r2; /* |rij| * rij^2 */
     
      double complex da[DIM], dj[DIM];
      
      for(int k = 0; k < DIM; ++k)
      {
        da[k] = dj[k] = 0.0;
      }
    
      /* calculates new accceleration and jerk for both particles i and j */
      for (int k = 0; k < DIM ; k++)
      {
        da[k] = rji[k] / r3;
        dj[k] = (vji[k] - 3 * (rv / r2) * rji[k]) / r3;
        
        acc[i + k] += mass[mj] * da[k]; /* add positive acceleration to particle i */
        acc[j + k] -= mass[mi] * da[k]; /* add negative acceleration to particle j */
       
        jerk[i + k] += mass[mj] * dj[k]; /* add positive jerk to particle i */                
        jerk[j + k] -= mass[mi] * dj[k]; /* add negative jerk to particle j */
      }
    }
  }
}

/* 
Implementation of the Hermite scheme, calculates new positions and velocities for all particles. 
For mathematical expressions used, see: http://www.ub.uni-heidelberg.de/archiv/6553 (Section 4.5) 
Further information: Kokubo E., Yoshinaga K., Makino J., 1998, MNRAS 297, 1067
*/
void hermite(int N, int DIM, double dt, double *mass, double complex *pos, 
             double complex *vel, double complex *acc, double complex *jerk)
{
  /* storing positions, velocities, acceleration and jerk from last iteration */
  double complex old_pos[N * DIM];
  double complex old_vel[N * DIM];  
  double complex old_acc[N * DIM];
  double complex old_jerk[N * DIM];
  
  memcpy(old_pos, pos, sizeof(old_pos));
  memcpy(old_vel, vel, sizeof(old_vel));
  memcpy(old_acc, acc, sizeof(old_acc));
  memcpy(old_jerk, jerk, sizeof(old_jerk));
  
  /* prediction for all particles (for mathematical expression please see links provided above) */
  for(int i = 0; i < N; i += 3)
  {
    for(int j = 0; j < DIM; ++j)
    {
      pos[i + j] += vel[i + j] * dt + acc[i + j] * ((dt * dt)/2) + jerk[i + j] * ((dt * dt * dt)/6);
      vel[i + j] += acc[i + j] * dt + jerk[i + j] * ((dt * dt)/2);
    }
  }
  
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk); /* get the new acceleration and jerk for all particles*/
  
  /* correction in reversed order of computation (for mathematical expression please see links provided above) 
     reversed order allows the corrected velocities to be used to correct the positions for better energy behaviour */
  for (int i = 0; i < N; i += 3)
  {
    for (int j = 0; j < DIM; ++j)
    {
      vel[i + j] = old_vel[i + j] + (old_acc[i + j] + acc[i + j]) * (dt/2) + (old_jerk[i + j] - jerk[i + j]) * ((dt * dt)/12);       
      pos[i + j] = old_pos[i + j] + (old_vel[i + j] + vel[i + j]) * (dt/2) + (old_acc[i + j] - acc[i + j]) * ((dt * dt)/12);
    }
  }
}

/* entry point for the Hermite integrator, starts computation and continues until end of simulation is reached */
void startHermite(int N, int DIM, double dt, double end_time, double *mass, double complex *pos, double complex *vel, 
                  double complex *acc, double complex *jerk)
{
  double time = 0.0; /* default time */
  int iterations = 0; /* iteration counter, iteration 0 is equal to initial conditions */
  
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk); /* one time calculation to get inital acceleration and jerk for all particles */
  energy_diagnostics(N, DIM, mass, pos, vel); /* get energy diagnostics for initial conditions */
  
  while(time < end_time) /* until user specified end of simulation is reached */
  {
    ++iterations; /* increment iteration counter from last iteration to current iteration */
    hermite(N, DIM, dt, mass, pos, vel, acc, jerk); /* calculate movement for current iteration */
    printIteration(N, iterations, mass, pos, vel); /* print current iteration */
    energy_diagnostics(N, DIM, mass, pos, vel); /* get energy diagnostics for current iteration */
    time += dt; /* add timestep to current time to advance to next iteration */
  }
}
