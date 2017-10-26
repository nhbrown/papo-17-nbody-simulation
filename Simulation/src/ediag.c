/*    
    The following source code provides methods for energy diagnostics by
    providing functions to calculate kinetic, potential and total energy
    of the cluster.
    
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
#include <math.h>
#include "output.h"
#include <stdio.h>

/* kinetic, potential and total energy of the cluster */
double e_kinetic, e_potential, e_total;

/*
 * Function:  kinetic_energy 
 * ====================
 *  Entry point for energy diagnostics, calls all other functions,
 *  calulates total energy and calls printEnergyDiagnostic.
 *
 *  N: amout of particles
 *  DIM: dimensions of space
 *  mass: masses of all particles
 *  pos: positions of all particles
 *  vel: velocities of all particles
 *
 *  returns: void
 * --------------------
 */
void energy_diagnostics(int N, int DIM, double *mass, double complex *pos, double complex *vel)
{
  kinetic_energy(N, DIM, mass, vel);
  
  potential_energy(N, DIM, mass, pos);
  
  e_total = e_kinetic + e_potential;
  
  printEnergyDiagnostics(e_kinetic, e_potential, e_total); /* provided by output.h */
}

/*
 * Function:  kinetic_energy 
 * ====================
 *  Calculates kinetic energy of the cluster.
 *
 *  N: amout of particles
 *  DIM: dimensions of space
 *  mass: masses of all particles
 *  vel: velocities of all particles
 *
 *  returns: void
 * --------------------
 */
void kinetic_energy(int N, int DIM, double *mass, double complex *vel)
{
  for(int i = 0, mi = 0; i < (N * DIM); i += DIM, ++mi)
  {
    for(int j = 0; j < DIM; ++j)
    {
      e_kinetic += 0.5 * mass[mi] * vel[i + j] * vel[i + j];
    }
  }
}

/*
 * Function:  potential_energy 
 * ====================
 *  Calculates potential energy of the cluster.
 *
 *  N: amout of particles
 *  DIM: dimensions of space
 *  mass: masses of all particles
 *  pos: positions of all particles
 *
 *  returns: void
 * --------------------
 */
void potential_energy(int N, int DIM, double *mass, double complex *pos)
{
  for (int i = 0, mi = 0; i < (N * DIM) ; i += DIM, ++mi)
  {
    for (int j = i + DIM, mj = 0; j < (N * DIM) ; j += DIM, ++mj)
    {
      double complex rij2 = 0;
    
      for (int k = 0; k < DIM ; ++k)
      {
        rij2 += (pos[j + k] - pos[i + k]) * (pos[j + k] - pos[i + k]);
      }
      
      e_potential -= mass[mi] * mass[mj] / sqrt(rij2);
    }
  }
}
