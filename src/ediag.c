#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "ediag.h"

double e_kinetic = 0.0;
double e_potential = 0.0;
double e_total = 0.0;

void kinetic_energy(int N, int DIM, double *mass, double complex **vel)
{
  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < DIM; ++j)
    {
      e_kinetic += 0.5 * mass[i] * vel[i][j] * vel[i][j];
    }
  }
}

void potential_energy(int N, int DIM, double *mass, double complex **pos)
{
  for (int i = 0; i < N ; ++i)
  {
    for (int j = i+1; j < N ; ++j)
    {
      double complex rji[DIM];
      double complex r2;
    
      for (int k = 0; k < DIM ; ++k)
      {
        rji[k] = pos[j][k] - pos[i][k];
      }
    
      for (int k = 0; k < DIM ; ++k)
      {
        r2 += rji[k] * rji[k];
      }
      
      e_potential -= mass[i] * mass[j] / sqrt(r2);
    }
  }
}

void energy_diagnostics(int N, int DIM, double *mass, double complex **pos, double complex **vel)
{
  kinetic_energy(N, DIM, mass, vel);
  
  potential_energy(N, DIM, mass, pos);
  
  e_total = e_kinetic + e_potential;
  
  printf("Kinetic Energy: %f \nPotential Energy: %f \nTotal Energy: %f \n", e_kinetic, e_potential, e_total);
}
