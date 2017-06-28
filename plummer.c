#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Seed for srand. */
double seed;
/* Total number of particles to be generated. */
int N;

struct body
{
  double xpos;
  double ypos;
  double zpos;
  
  double mass;  
  
  double xvel;
  double yvel;
  double zvel;
}

/* single particle */
body p;

/* 
M defines the total mass of the cluster.
R defines the dimensions of the cluster.
G defines the gravitational constant.
*/
static const double M = 1.0;
static const double R = 1.0;
static const double G = 1.0;

double frand(double low, double high)
{
  return low + rand * (high - low);
}

void plummer()
{
  srand(seed);
  
  p.mass = M / N;
  
  double radius = R / sqrt((pow(rand(), (-2.0/3.0))) - 1.0);
  double theta = acos(frand(-1.0, 1.0));
  double phi = frand(0.0, (2 * M_PI));
  
  p.xpos = radius * sin(theta) * cos(phi);
  p.ypos = radius * sin(theta) * sin(phi);
  p.zpos = radius * cos(theta);
  
  double x = 0.0;
  double y = 0.1;
  
  while(y > pow((x * x * (1.0 - x * x)), 3.5))
  {
    x = frand(0.0, 1.0);
    y = frand(0.0, 0.1);
  }
  
  double velocity = x * sqrt(2.0) * pow((1.0 + radius * radius), -0.25));
  
  p.xvel = velocity * sin(theta) * cos(phi);
  p.yvel = velocity * sin(theta) * sin(phi);
  p.zvel = velocity * cos(theta);
}

/*
Passing a seed as a parameter is optional, if no seed is passed seed is equal to UNIX-clock.
Specifying the amount of particles to generate is always necessary.
If user wishes to specify the seed, the order of arguments needs to be: <executable> seed amount
*/
int main(int argc, const char *argv[])
{
  switch(argc)
  {
    case 2 : /* if one argument is passed, it is assumed to be amount of particles */
      seed = (double)time(NULL);
      N = atoi(argv[1]);
      break;

    case 3 : /* if two arguments are passed, first one is assumed to be seed */
      seed = atol(argv[1]);
      N = atoi(argv[2]);
      break;

    default : /* if less than 1 or more than 2 arguments are passed, the executions exits */
      printf("Invalid input for plummer.c!\n");
      exit(0);
  }
  
  FILE *log;
  log = fopen("log.txt", "w"); /* writes to new file log.txt which holds important parameters */

  fprintf(log, "Seed used: %d \nNumber of particles: %d \nTotal mass of cluster: %d \nDimensions of cluster: %d \n
          Gravitational constant: %d", seed, N, M, R, G);

  fclose(log);

  FILE *output;
  output = fopen("output.csv", "w"); /* writes to new file output.csv which holds positions */

  for(int j = 0; j < N; ++j)
  {
    plummer();
    fprintf(output, "%d, %d, %d, %d, %d, %d, %d \n\n", p.xpos, p.ypos, p.zpos, p.mass, p.xvel, p.yvel, p.zvel);
  }

  fclose(output);
  
  return 0;
}
