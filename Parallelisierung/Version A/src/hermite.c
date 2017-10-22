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
#include <mpi.h>
#include "output.h"
#include <stdlib.h>
#include <string.h>

int world_rank, world_size, proc_elem;

/* calculates acceleration and jerk for all particles */
void acc_jerk(int N, int DIM, double *mass, double complex *pos, double complex *vel, 
              double complex *acc, double complex *jerk)
{ 
  /* default values for acceleration and jerk */
  for(int i = 0; i < (N * DIM); ++i)
  {
    acc[i] = jerk[i] = 0;
  }
  
  double *local_mass = calloc((N / world_size), sizeof(double));
  
  double complex *local_pos = calloc(proc_elem, sizeof(double complex));
  double complex *local_vel = calloc(proc_elem, sizeof(double complex));
  
  double complex *local_acc = calloc(proc_elem, sizeof(double complex));
  double complex *local_jerk = calloc(proc_elem, sizeof(double complex));
  
  MPI_Scatter(mass, (N / world_size), MPI_DOUBLE, local_mass, (N / world_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Scatter(pos, proc_elem, MPI_C_DOUBLE_COMPLEX, local_pos, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Scatter(vel, proc_elem, MPI_C_DOUBLE_COMPLEX, local_vel, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  MPI_Scatter(acc, proc_elem, MPI_C_DOUBLE_COMPLEX, local_acc, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Scatter(jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, local_jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  for(int i = 0, mi = 0; i < (N * DIM); i += DIM, ++mi) /* loops over all particles */
  { 
    for(int j = i + DIM, mj = 0; j < (N * DIM); j += DIM, ++mj) /* only loops over half of the particles because force acts equally on both particles (Newton) */
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
  
  MPI_Gather(local_mass, (N / world_size), MPI_DOUBLE, mass, (N / world_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Gather(local_pos, proc_elem, MPI_DOUBLE, pos, proc_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(local_vel, proc_elem, MPI_DOUBLE, vel, proc_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Gather(local_acc, proc_elem, MPI_DOUBLE, acc, proc_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(local_jerk, proc_elem, MPI_DOUBLE, jerk, proc_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
  for(int i = 0; i < (N * DIM); ++i)
  {
    pos[i] += vel[i] * dt + acc[i] * ((dt * dt)/2) + jerk[i] * ((dt * dt * dt)/6);
    vel[i] += acc[i] * dt + jerk[i] * ((dt * dt)/2);
  }
  
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk); /* get the new acceleration and jerk for all particles*/
  
  /* correction in reversed order of computation (for mathematical expression please see links provided above) 
     reversed order allows the corrected velocities to be used to correct the positions for better energy behaviour */
  for (int i = 0; i < (N * DIM); ++i)
  {
    vel[i] = old_vel[i] + (old_acc[i] + acc[i]) * (dt/2) + (old_jerk[i] - jerk[i]) * ((dt * dt)/12);       
    pos[i] = old_pos[i] + (old_vel[i] + vel[i]) * (dt/2) + (old_acc[i] - acc[i]) * ((dt * dt)/12);
  }
}

/* entry point for the Hermite integrator, starts computation and continues until end of simulation is reached */
void startHermite(int N, int DIM, double dt, double end_time, double *mass, double complex *pos, double complex *vel, 
                  double complex *acc, double complex *jerk, int rank, int size, int elements)
{
  double time = 0.0; /* default time */
  int iterations = 0; /* iteration counter, iteration 0 is equal to initial conditions */
  
  world_rank = rank;
  world_size = size;
  proc_elem = elements;
  
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk); /* get inital acceleration and jerk for all particles */
  
  if(world_rank == 0)
  {
    energy_diagnostics(N, DIM, mass, pos, vel); /* get energy diagnostics for initial conditions */
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  while(time < end_time) /* until user specified end of simulation is reached */
  {
    ++iterations; /* increment iteration counter from last iteration to current iteration */  
    hermite(N, DIM, dt, mass, pos, vel, acc, jerk); /* calculate movement for current iteration */
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(world_rank == 0)
    {
      printIteration(N, DIM, iterations, mass, pos, vel); /* print current iteration */
      energy_diagnostics(N, DIM, mass, pos, vel); /* get energy diagnostics for current iteration */
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    time += dt; /* add timestep to current time to advance to next iteration */
  }
}
