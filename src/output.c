#include <stdio.h>
#include <stdlib.h>

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

void generateOutput(double timestep, double end_time)
{ 
  createNames();
  
  FILE *log;
  log = fopen(logname, "w"); /* writes to new file log_<currentdate>.txt which holds important parameters */

  fprintf(log, "Seed used: %lu \nNumber of particles: %d \nTotal mass of cluster: %f \nDimensions of cluster: %f \nGravitational constant: %f \nTimestep: %f \nEndtime: %f", 
          seed, N, M, R, G, timestep, end_time);

  fclose(log);

  FILE *conditions;
  conditions = fopen(conditionsname, "w"); /* writes to new file initial_conditions.csv which holds positions */

  /*
  for(int j = 0; j < N; ++j)
  {
    plummer();
    fprintf(conditions, "%f, %f, %f, %f, %f, %f, %f\n", 
            creal(p.xpos), creal(p.ypos), creal(p.zpos), p.mass, creal(p.xvel), creal(p.yvel), creal(p.zvel));
  }
  */

  fclose(conditions);
}
