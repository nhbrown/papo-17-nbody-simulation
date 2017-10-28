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
#include "ediag.h"
#include "hermite.h"
#include "output.h"
#include <string.h>

/*
 * Function:  startHermite 
 * ====================
 *  Entry point for the Hermite scheme. Controls current computation
 *  and checks wether or not end of simulation has been reached.
 *
 *  N: amount of particles
 *  DIM: dimensions of space
 *  dt: timestep
 *  end_time: end of simulation
 *  mass: masses of all particles
 *  pos: positions of all particles
 *  vel: velocity of all particles
 *  acc: acceleration for all particles
 *  jerk: jerk for all particles
 *
 *  returns: void
 * --------------------
 */
void startHermite(int N, int DIM, double dt, double end_time, double *mass, double complex *pos, double complex *vel, 
                  double complex *acc, double complex *jerk)
{
  double time = 0.0; /* default time */
  int iterations = 0; /* iteration counter, iteration 0 is equal to initial conditions */
  
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk); /* calculate inital acceleration and jerk for all particles */
  energy_diagnostics(N, DIM, mass, pos, vel); /* calculate energy diagnostics for initial conditions */
  
  /* continues until specified end of simulation is reached */
  while(time < end_time)
  {
    ++iterations; /* increment iteration counter from last iteration to current iteration */
    
    hermite(N, DIM, dt, mass, pos, vel, acc, jerk); /* calculate movement for current iteration */
    printIteration(N, DIM, iterations, mass, pos, vel); /* provided by output.h */
    energy_diagnostics(N, DIM, mass, pos, vel); /* provided by ediag.h */
    
    time += dt; /* add timestep to current time to advance to next iteration */
  }
}

/*
 * Function:  acc_jerk 
 * ====================
 *  Calculates acceleration and jerk for all particles by
 *  comparing them pairwise. Comparison is optimized by only
 *  calculating pairwise acceleration and jerk once and adding
 *  them to both particles. Cuts down computation time by 
 *  approximately half.
 *
 *  N: amount of particles
 *  DIM: dimensions of space
 *  mass: masses of all particles
 *  pos: positions of all particles
 *  vel: velocity of all particles
 *  acc: acceleration for all particles
 *  jerk: jerk for all particles
 *
 *  returns: void
 * --------------------
 */
void acc_jerk(int N, int DIM, double *mass, double complex *pos, double complex *vel, 
              double complex *acc, double complex *jerk)
{ 
  /* default values for acceleration and jerk */
  for(int i = 0; i < (N * DIM); ++i)
  {
    acc[i] = jerk[i] = 0;
  }
  
  /* loops over all particles */
  for(int i = 0, mi = 0; i < (N * DIM); i += DIM, ++mi)
  { 
    /* only loops over half of the particles because force acts equally on both particles (Newton) */
    for(int j = i + DIM, mj = 0; j < (N * DIM); j += DIM, ++mj)
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
 * Function:  hermite 
 * ====================
 *  Implementation of the Hermite scheme, calculates new positions 
 *  and velocities for all particles. 
 *  Based on Kokubo E., Yoshinaga K., Makino J., 1998, MNRAS 297, 1067
 *
 *  N: amount of particles
 *  DIM: dimensions of space
 *  dt: timestep
 *  mass: masses of all particles
 *  pos: positions of all particles
 *  vel: velocity of all particles
 *  acc: acceleration for all particles
 *  jerk: jerk for all particles
 *
 *  returns: void
 * --------------------
 */
void hermite(int N, int DIM, double dt, double *mass, double complex *pos, 
             double complex *vel, double complex *acc, double complex *jerk)
{
  /* storing positions, velocities, acceleration and jerk from last iteration */
  double complex *old_pos = calloc((N * DIM), sizeof(double complex));
  double complex *old_vel = calloc((N * DIM), sizeof(double complex));
  double complex *old_acc = calloc((N * DIM), sizeof(double complex));
  double complex *old_jerk = calloc((N * DIM), sizeof(double complex));
  
  /* allocation guard */
  if(old_pos == NULL || old_vel == NULL || old_acc == NULL || old_jerk == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
  
  /* copy data from last iteration */
  memcpy(old_pos, pos, (N * DIM));
  memcpy(old_vel, vel, (N * DIM));
  memcpy(old_acc, acc, (N * DIM));
  memcpy(old_jerk, jerk, (N * DIM));
  
  /* prediction for all particles using old values*/
  for(int i = 0; i < (N * DIM); ++i)
  {
    pos[i] += vel[i] * dt + acc[i] * ((dt * dt)/2) + jerk[i] * ((dt * dt * dt)/6);
    vel[i] += acc[i] * dt + jerk[i] * ((dt * dt)/2);
  }
  
  /* calculate new acceleration and jerk for all particles*/
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk);
  
  /* correction in reversed order of computation, for allows the corrected velocities 
     to be used to correct the positions for better energy behaviour */
  for (int i = 0; i < (N * DIM); ++i)
  {
    vel[i] = old_vel[i] + (old_acc[i] + acc[i]) * (dt/2) + (old_jerk[i] - jerk[i]) * ((dt * dt)/12);       
    pos[i] = old_pos[i] + (old_vel[i] + vel[i]) * (dt/2) + (old_acc[i] - acc[i]) * ((dt * dt)/12);
  }
}
