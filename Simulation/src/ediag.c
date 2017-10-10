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
#include <math.h>
#include <stdio.h>
#include "ediag.h"
#include "output.h"

double e_kinetic = 0.0; /* kinetic energy of the cluster */
double e_potential = 0.0; /* potential energy of the cluster */
double e_total = 0.0; /* total energy of the cluster */

/* calculates kinetic energy of the cluster */
void kinetic_energy(int N, double *mass, double complex **vel)
{
  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      e_kinetic += 0.5 * mass[i] * vel[i][j] * vel[i][j];
    }
  }
}

/* calculates potential energy of the cluster */
void potential_energy(int N, double *mass, double complex **pos)
{
  for (int i = 0; i < N ; ++i)
  {
    for (int j = i+1; j < N ; ++j)
    {
      double complex rij2 = 0;
    
      for (int k = 0; k < 3 ; ++k)
      {
        rij2 += (pos[j][k] - pos[i][k]) * (pos[j][k] - pos[i][k]);
      }
      
      e_potential -= mass[i] * mass[j] / sqrt(rij2);
    }
  }
}

/* entry point for energy diagnostics, also calculates total energy of the cluster */
void energy_diagnostics(int N, int marker, double *mass, double complex **pos, double complex **vel)
{
  kinetic_energy(N, mass, vel);
  
  potential_energy(N, mass, pos);
  
  e_total = e_kinetic + e_potential;
  
  printEnergyDiagnostics(marker, e_kinetic, e_potential, e_total);
}
