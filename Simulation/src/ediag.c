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

double e_kinetic, e_potential, e_total; /* kinetic, potential and total energy of the cluster */

/* calculates kinetic energy of the cluster */
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

/* calculates potential energy of the cluster */
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

/* entry point for energy diagnostics, also calculates total energy of the cluster */
void energy_diagnostics(int N, int DIM, double *mass, double complex *pos, double complex *vel)
{
  kinetic_energy(N, DIM, mass, vel);
  
  potential_energy(N, DIM, mass, pos);
  
  e_total = e_kinetic + e_potential;
  
  printEnergyDiagnostics(e_kinetic, e_potential, e_total);
}
