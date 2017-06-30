#include <complex.h>
#include <math.h>
#include "mersenne.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Seed for srand. */
double seed;
/* Total number of particles to be generated. */
int N;

struct body
{
  double complex xpos;
  double complex ypos;
  double complex zpos;
  
  double mass;  
  
  double complex xvel;
  double complex yvel;
  double complex zvel;
};

/* single particle */
struct body p;

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
  return low + genrand_real1() * (high - low);
}

void plummer()
{  
  p.mass = M / N; /* mass equilibrium */
  
  double complex radius = R / csqrt((cpow(genrand_real1(), (-2.0/3.0))) - 1.0); /* inverted cumulative mass distribution */
  double complex theta = cacos(frand(-1.0, 1.0)); /* Polar Angle */
  double complex phi = frand(0.0, (2 * M_PI)); /* Azimuthal Angle */
  
  /* conversion from radial to cartesian coordinates */
  p.xpos = radius * csin(theta) * ccos(phi); 
  p.ypos = radius * csin(theta) * csin(phi);
  p.zpos = radius * ccos(theta);
  
  double x = 0.0;
  double y = 0.1;
  
  while(y > pow((x * x * (1.0 - x * x)), 3.5))
  {
    x = frand(0.0, 1.0);
    y = frand(0.0, 0.1);
  }
  
  double complex velocity = x * csqrt(2.0) * cpow((1.0 + radius * radius), -0.25);
  
  p.xvel = velocity * csin(theta) * ccos(phi);
  p.yvel = velocity * csin(theta) * csin(phi);
  p.zvel = velocity * ccos(theta);
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
      seed = atof(argv[1]);
      N = atoi(argv[2]);
      break;

    default : /* if less than 1 or more than 2 arguments are passed, the executions exits */
      printf("Invalid input for plummer.c!\n");
      exit(0);
  }
  
  srand(seed);
  
  FILE *log;
  log = fopen("log.txt", "w"); /* writes to new file log.txt which holds important parameters */

  fprintf(log, "Seed used: %f \nNumber of particles: %d \nTotal mass of cluster: %f \nDimensions of cluster: %f \nGravitational constant: %f", 
          seed, N, M, R, G);

  fclose(log);

  FILE *output;
  output = fopen("output.csv", "w"); /* writes to new file output.csv which holds positions */

  for(int j = 0; j < N; ++j)
  {
    plummer();
    fprintf(output, "%f, %f, %f, %f, %f, %f, %f\n\n", 
            creal(p.xpos), creal(p.ypos), creal(p.zpos), p.mass, creal(p.xvel), creal(p.yvel), creal(p.zvel));
  }

  fclose(output);
  
  return 0;
}
