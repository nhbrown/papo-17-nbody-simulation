#include <complex.h>
#include <stdio.h>
#include "ediag.h"

double e_kinetic = 0.0;
double e_potential = 0.0;
double e_total = 0.0;

void kinetic_energy(int N, double complex *mass, double complex **vel)
{
  for(int i = 0; i < N; ++i)
  {
    double complex vel2 = (cpow (vel[i][0] + vel[i][1] + vel[i][2]) 2.0);
    e_kinetic += 0.5 * mass[i] * creal(vel2);
  }
}

void potential_energy(int N, double complex *mass, double complex **pos)
{
  for(int i = 0; i < N; ++i)
  {
    double complex x = 0;
    for(int j = 0; j < N; ++j)
    {
      if(j != i)
      {
        double complex y = (pos[i][0] + pos[i][1] + pos[i][2]) - (pos[j][0] + pos[j][1] + pos[j][2]);
        double complex z = (-1 * mass[j]) * mass[i];
        x += z / csqrt(y * y);
      }
    }
  }
  e_potential = creal(x) / 2;
}

void energy_diagnostics(int N, double complex *mass, double complex **pos, double complex **vel)
{
  kinetic_energy(N, mass, vel);
  
  potential_energy(N, mass, pos);
  
  e_total = e_kinetic + e_potential;
  
  printf("Kinetic Energy: %f \n Potential Energy: %f \n Total Energy: %f \n", e_kinetic, e_potential, e_total);
}
