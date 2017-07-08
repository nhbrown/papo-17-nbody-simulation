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

int N; /* amount of particles */
int DIM; /* dimensions */
double dt; /* timestep */

void acc_jerk(double *mass, double complex (*pos)[DIM], double complex (*vel)[DIM], 
              double complex (*acc)[DIM], double complex (*jerk)[DIM])
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

void hermite(double *mass, double complex (*pos)[DIM], double complex (*vel)[DIM], 
             double complex (*acc)[DIM], double complex (*jerk)[DIM])
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

int main(int argc, const char *argv[])
{
  N = 3;
  DIM = 3;
  dt = 0.01;
  
  double end_time = 1.0;
  double time = 0.0;
  
  double mass[N]; /* mass for all particles */

  double complex pos[N][DIM]; /* positions for all particles */
  double complex vel[N][DIM]; /* velocities for all particles */

  double complex acc[N][DIM]; /* acceleration for all particles */
  double complex jerk[N][DIM]; /* jerk for all particles */
  
  mass[0] = 0.000007;
  mass[1] = 0.000007;
  mass[2] = 0.000007;
  
  pos[0][0] = 0.126388;
  pos[0][1] = 0.074997;
  pos[0][2] = -0.151221;
  
  pos[1][0] = -1.706910;
  pos[1][1] = 1.302327;
  pos[1][2] = 1.616718;
  
  pos[2][0] = 1.489653;
  pos[2][1] = -0.034517;
  pos[2][2] = -1.582992;
  
  vel[0][0] = 0.353411;
  vel[0][1] = 0.237369;
  vel[0][2] = -0.924605;
  
  vel[1][0] = -0.180430;
  vel[1][1] = 0.497624;
  vel[1][2] = -0.009716;
  
  vel[2][0] = -0.156375;
  vel[2][1] = -0.581931;
  vel[2][2] = -0.005356;
  
  double *pmass = mass;
  
  double complex (*ppos)[DIM] = pos;
  double complex (*pvel)[DIM] = vel;
  
  double complex (*pacc)[DIM] = acc;
  double complex (*pjerk)[DIM] = jerk;
  
  acc_jerk(pmass, ppos, pvel, pacc, pjerk);
  
  while(time < end_time)
  {
    hermite(pmass, ppos, pvel, pacc, pjerk);
    time += dt;
    for(int i = 0; i < N; ++i)
    {
      for(int k = 0; k < 1; ++k)
      {        
        printf("Particle %d: %f %f %f %f %f %f %f \n", i, 
               mass[i], creal(pos[i][k]), creal(pos[i][k + 1]), creal(pos[i][k + 2]), 
               creal(vel[i][k]), creal(vel[i][k + 1]), creal(vel[i][k + 2]));
      }
    }
  }
}
