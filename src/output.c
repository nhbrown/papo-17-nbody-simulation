#include <stdio.h>
#include <time.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "output.h"

char foldername[40]; /* buffer for foldername */ 
char logname[80]; /* buffer for logname */
char conditionsname[80]; /* buffer for conditionsname */

void createNames()
{
  struct tm *sTm;

  time_t now = time(0);
  sTm = gmtime(&now);

  strftime (foldername, sizeof(foldername), "run_%Y_%m_%d_%H:%M:%S", sTm);
  strftime (logname, sizeof(logname), "run_%Y_%m_%d_%H:%M:%S/log_%Y_%m_%d_%H:%M:%S.txt", sTm);
  strftime (conditionsname, sizeof(conditionsname), "run_%Y_%m_%d_%H:%M:%S/initial_conditions.csv", sTm);

  struct stat st = {0};

  if (stat(foldername, &st) == -1)
  {
      mkdir(foldername, 0700);
  }
}

void printInitialConditions(int N, double *mass, double complex **pos, double complex **vel)
{ 
  FILE *conditions;
  conditions = fopen(conditionsname, "w"); /* writes to new file initial_conditions.csv which holds positions */

  for(int i = 0; i < N; ++i)
  {
    for(int k = 0; k < 1; ++k)
    {        
      fprintf(conditions, "%f, %f, %f, %f, %f, %f, %f \n", 
              creal(pos[i][k]), creal(pos[i][k + 1]), creal(pos[i][k + 2]), mass[i],
              creal(vel[i][k]), creal(vel[i][k + 1]), creal(vel[i][k + 2]));
    }
  }

  fclose(conditions);
}

void printLog(unsigned long seed, int N, double M, double R, double G, double timestep, double end_time)
{
  FILE *log;
  log = fopen(logname, "w"); /* writes to new file log_<currentdate>.txt which holds important parameters */

  fprintf(log, "Seed used: %lu \nNumber of particles: %d \nTotal mass of cluster: %f \nDimensions of cluster: %f \nGravitational constant: %f \nTimestep: %f \nEndtime: %f", 
          seed, N, M, R, G, timestep, end_time);

  fclose(log);
}

void printIteration(double *mass, double complex **pos, double complex **vel, int iteration, int N)
{
  char buffer[80];
  snprintf(buffer, sizeof(buffer), "./%s/iteration_%d.csv", foldername, iteration);
  
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
