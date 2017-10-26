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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* rank of process, amount of processes and elements per process */
int world_rank, world_size, proc_elem;

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
 *  rank: rank of process
 *  size: amount of processes
 *  elements: amount of elements per process
 *
 *  returns: void
 * --------------------
 */
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
  
  /* continues until specified end of simulation is reached */
  while(time < end_time)
  {
    ++iterations; /* increment iteration counter from last iteration to current iteration */ 
    
    hermite(N, DIM, dt, mass, pos, vel, acc, jerk); /* calculate movement for current iteration */
    
    if(world_rank == 0)
    {
      printIteration(N, DIM, iterations, mass, pos, vel); /* provided by output.h */
      energy_diagnostics(N, DIM, mass, pos, vel); /* provided by ediag.h */
    }
    
    time += dt; /* add timestep to current time to advance to next iteration */
  }
}

/*
 * Function:  acc_jerk 
 * ====================
 *  Calculates acceleration and jerk for all particles by
 *  comparing them pairwise. Comparison is not optimized,
 *  all pairwise comparisons are calculated twice.
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
  
  /* allocate space for local particles */
  double complex *local_pos = calloc(proc_elem, sizeof(double complex));
  double complex *local_vel = calloc(proc_elem, sizeof(double complex));
  
  double complex *local_acc = calloc(proc_elem, sizeof(double complex));
  double complex *local_jerk = calloc(proc_elem, sizeof(double complex));
  
  /* allocation guard */
  if(local_pos == NULL || local_vel == NULL || local_acc == NULL || local_jerk == NULL)
  {
    fprintf(stderr, "Out of memory!\n");
    exit(0);
  }
  
  /* scatter particle data to processes */
  MPI_Scatter(pos, proc_elem, MPI_C_DOUBLE_COMPLEX, local_pos, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Scatter(vel, proc_elem, MPI_C_DOUBLE_COMPLEX, local_vel, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  MPI_Scatter(acc, proc_elem, MPI_C_DOUBLE_COMPLEX, local_acc, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Scatter(jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, local_jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  /* broadcast current position and velocity to all processes */
  MPI_Bcast(pos, (N * DIM), MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vel, (N * DIM), MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  for(int i = 0, mi = 0; i < proc_elem; i += DIM, ++mi) /* loops over all local particles */
  { 
    for(int j = 0, mj = 0; j < (proc_elem * world_size); j += DIM, ++mj) /* loops over all particles */
    {
      /* guard - passed if particle i and j are not the same */
      if(local_pos[i] != pos[j] && local_pos[i + 1] != pos[j + 1] && local_pos[i + 2] != pos[j + 2]
        && local_vel[i] != vel[j] && local_vel[i + 1] != vel[j + 1] && local_vel[i + 2] != vel[j + 2])
      {
        double complex rji[DIM], vji[DIM]; /* position vector from particle i to j */
      
        for(int k = 0; k < DIM; ++k)
        {
          rji[k] = vji[k] = 0.0;
        }
     
        double complex r2 = 0.0; /* rij^2 */
        double complex rv = 0.0; /* rij*vij */
     
        /* calculate position and velocity vectors */
        for(int k = 0; k < DIM; ++k)
        {
          rji[k] = pos[j + k] - local_pos[i + k];
          vji[k] = vel[j + k] - local_vel[i + k];
       
          r2 += rji[k] * rji[k];
          rv += rji[k] * vji[k];
        }
      
        double complex r3 = csqrt(r2) * r2; /* |rij| * rij^2 */
     
        double complex da[DIM], dj[DIM];
      
        for(int k = 0; k < DIM; ++k)
        {
          da[k] = dj[k] = 0.0;
        }
    
        /* calculates new accceleration and jerk for particle i */
        for (int k = 0; k < DIM ; k++)
        {
          da[k] = rji[k] / r3;
          dj[k] = (vji[k] - 3 * (rv / r2) * rji[k]) / r3;
        
          local_acc[i + k] += mass[mj] * da[k]; /* add positive acceleration to particle i */
          local_jerk[i + k] += mass[mj] * dj[k]; /* add positive jerk to particle i */                
        }
      }
    }
  }
  
  /* gather local particles back to global arrays on root */
  MPI_Gather(local_pos, proc_elem, MPI_C_DOUBLE_COMPLEX, pos, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Gather(local_vel, proc_elem, MPI_C_DOUBLE_COMPLEX, vel, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  MPI_Gather(local_acc, proc_elem, MPI_C_DOUBLE_COMPLEX, acc, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Gather(local_jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, jerk, proc_elem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
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
  double complex old_pos[N * DIM];
  double complex old_vel[N * DIM];  
  double complex old_acc[N * DIM];
  double complex old_jerk[N * DIM];
  
  memcpy(old_pos, pos, sizeof(old_pos));
  memcpy(old_vel, vel, sizeof(old_vel));
  memcpy(old_acc, acc, sizeof(old_acc));
  memcpy(old_jerk, jerk, sizeof(old_jerk));
  
  if(world_rank == 0)
  { 
    /* prediction for all particles using old values*/
    for(int i = 0; i < (N * DIM); ++i)
    {
      pos[i] += vel[i] * dt + acc[i] * ((dt * dt)/2) + jerk[i] * ((dt * dt * dt)/6);
      vel[i] += acc[i] * dt + jerk[i] * ((dt * dt)/2);
    }
  }
  
  /* calculate new acceleration and jerk for all particles*/
  acc_jerk(N, DIM, mass, pos, vel, acc, jerk);
  
  if(world_rank == 0)
  {
    /* correction in reversed order of computation, for allows the corrected velocities 
       to be used to correct the positions for better energy behaviour */
    for (int i = 0; i < (N * DIM); ++i)
    {
      vel[i] = old_vel[i] + (old_acc[i] + acc[i]) * (dt/2) + (old_jerk[i] - jerk[i]) * ((dt * dt)/12);       
      pos[i] = old_pos[i] + (old_vel[i] + vel[i]) * (dt/2) + (old_acc[i] - acc[i]) * ((dt * dt)/12);
    }
  }
}
