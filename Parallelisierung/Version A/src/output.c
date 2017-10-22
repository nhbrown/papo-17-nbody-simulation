/*    
    The following source code provides methods for writing initial conditions,
    important information and iterations of the computation to *.txt or *.csv files.
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
#include "output.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

char foldername[40]; /* buffer for name of folder */ 
char logname[80]; /* buffer for name of log file */
char conditionsname[80]; /* buffer for name of intial conditions file */
char ediagname[80]; /* buffer for name of energy diagnostics file */

/* creates names and folder for files of the currently running simulation */
void createNames()
{
  struct tm *sTm;

  time_t now = time(0);
  sTm = gmtime(&now);

  strftime (foldername, sizeof(foldername), "run_%Y_%m_%d_%H:%M:%S", sTm);
  strftime (logname, sizeof(logname), "run_%Y_%m_%d_%H:%M:%S/log_%Y_%m_%d_%H:%M:%S.txt", sTm);
  strftime (conditionsname, sizeof(conditionsname), "run_%Y_%m_%d_%H:%M:%S/initial_conditions.csv", sTm);
  strftime (ediagname, sizeof(ediagname), "run_%Y_%m_%d_%H:%M:%S/energy_diagnostics.csv", sTm);

  struct stat st = {0};

  if (stat(foldername, &st) == -1)
  {
      mkdir(foldername, 0700);
  }
}

/* creates and writes to inital conditions file, which holds intial masses, positions and velocities of all particles */
void printInitialConditions(int N, int DIM, double *mass, double complex *pos, double complex *vel)
{ 
  FILE *conditions;
  conditions = fopen(conditionsname, "w");

  for(int i = 0, mi = 0; i < (N * DIM); i += 3, ++mi)
  {     
    fprintf(conditions, "%f, %f, %f, %f, %f, %f, %f \n",
            creal(pos[i]), creal(pos[i + 1]), creal(pos[i + 2]), mass[mi], 
            creal(vel[i]), creal(vel[i + 1]), creal(vel[i + 2]));
  }

  fclose(conditions);
}

/* creates and writes to log file, which holds important information about the current run */
void printLog(unsigned long seed, int N, double M, double R, double G, double timestep, double end_time)
{
  FILE *log;
  log = fopen(logname, "w"); /* writes to new file log_<currentdate>.txt which holds important parameters */

  fprintf(log, "Seed used: %lu \nNumber of particles: %d \n\nTotal mass of cluster: %f \nDimensions of cluster: %f \nGravitational constant: %f \n\nTimestep: %f \nEndtime: %f \n", 
          seed, N, M, R, G, timestep, end_time);

  fclose(log);
}

/* writes energy diagnostics to file */
void printEnergyDiagnostics(double e_kinetic, double e_potential, double e_total)
{
  FILE *ediag;
  
  if(access(ediagname, F_OK) != -1)
  {
    ediag = fopen(ediagname, "a");
    fprintf(ediag, "%f, %f, %f \n", e_kinetic, e_potential, e_total);
  } 
  else 
  {
    ediag = fopen(ediagname, "w");
    fprintf(ediag, "%f, %f, %f \n", e_kinetic, e_potential, e_total);
  }
  
  fclose(ediag);
}

/* creates and writes to new file for specified iteration of the computation */
void printIteration(int N, int DIM, int iteration, double *mass, double complex *pos, double complex *vel)
{
  char buffer[80];
  snprintf(buffer, sizeof(buffer), "./%s/iteration_%d.csv", foldername, iteration);
  
  FILE *out;
  out = fopen(buffer, "w");
  
  for(int i = 0, mi = 0; i < (N * DIM); i += 3, ++mi)
  {
    fprintf(out, "%f, %f, %f, %f, %f, %f, %f \n", 
            creal(pos[i]), creal(pos[i + 1]), creal(pos[i + 2]), mass[mi], 
            creal(vel[i]), creal(vel[i + 1]), creal(vel[i + 2]));
  }
  
  fclose(out);
}
