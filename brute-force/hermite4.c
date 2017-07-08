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
#include <string.h>

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

void readConditions(double *mass, double complex (*pos)[DIM], double complex (*vel)[DIM])
{
  FILE *inp;  
  inp = fopen("./run/initial_conditions.csv", "r");
    
  for(int i = 0; i < N; ++i)
  {
    char buffer[75];
    double values[7];
    int index = 0;
   
    char *delim = ", ";
    char *token = NULL;
    
    fgets(buffer, 75, inp);

    for (token = strtok(buffer, delim); token != NULL; token = strtok(NULL, delim))
    {
      char *ptr;
      double value = strtod(token, &ptr);
      values[index] = value;
      ++index;
    }
    
    pos[i][0] = values[0];
    pos[i][1] = values[1];
    pos[i][2] = values[2];
    
    mass[i] = values[3];
    
    vel[i][0] = values[4];
    vel[i][1] = values[5];
    vel[i][2] = values[6];
  }
}

void printIteration(double *mass, double complex (*pos)[DIM], double complex (*vel)[DIM], int iteration)
{
  char buffer[60];
  snprintf(buffer, sizeof(buffer), "./run/iteration_%d.csv", iteration);
  
  FILE *out;
  out = fopen(buffer, "w");
  
  for(int i = 0; i < N; ++i)
  {
    for(int k = 0; k < 1; ++k)
    {        
      fprintf(out, "%f, %f, %f, %f, %f, %f, %f \n", 
              creal(pos[i][k]), creal(pos[i][k + 1]), creal(pos[i][k + 2]), mass[i],
              creal(vel[i][k]), creal(vel[i][k + 1]), creal(vel[i][k + 2]));
    }
  }
  
  fclose(out);
}

int main(int argc, const char *argv[])
{
  double end_time;
  double time = 0.0;
  
  switch(argc)
  {
    case 5 :
      N = atoi(argv[1]);
      DIM = atoi(argv[2]);
      dt = atof(argv[3]);
      end_time = atof(argv[4]);
      break;

    default : /* if less or more than 5 arguments are passed, the execution exits normally */
      printf("Invalid input for hermite4.c!\n");
      exit(0);
  }
  
  double mass[N]; /* mass for all particles */

  double complex pos[N][DIM]; /* positions for all particles */
  double complex vel[N][DIM]; /* velocities for all particles */

  double complex acc[N][DIM]; /* acceleration for all particles */
  double complex jerk[N][DIM]; /* jerk for all particles */

  double *pmass = mass;
  
  double complex (*ppos)[DIM] = pos;
  double complex (*pvel)[DIM] = vel;
  
  double complex (*pacc)[DIM] = acc;
  double complex (*pjerk)[DIM] = jerk;
    
  readConditions(pmass, ppos, pvel);
  
  acc_jerk(pmass, ppos, pvel, pacc, pjerk);
  
  int iterations = 0;
  while(time < end_time)
  {
    hermite(pmass, ppos, pvel, pacc, pjerk);
    time += dt;
    ++iterations;
    printIteration(pmass, ppos, pvel, iterations);
  }
}
