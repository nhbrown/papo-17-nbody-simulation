#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Seed for the RNG, equal seeds generate equal sequences. */
long seed;
/* Number of particles passed to main. */
int i;

/* Constant values that are added by the RNG. */
const int m = 2147483647; /* modulos */
const int a = 48271;  /* multiplier */
const int q = 44488; /* (m / a) */
const int r = 3399; /* (m % a) */

/* 
Minimum LCG as defined by Park & Miller.
Returns a pseudo-random long in the range of +1 to +2147483647.
*/
int rand_pos()
{ 
  seed = a * (seed % q) - r * (seed / q);
  
  if(seed <= 0)
  {
    seed += m;
  }
  
  return seed;
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
      seed = (unsigned long)time(NULL);
      i = atoi(argv[1]);
      break;

    case 3 : /* if two arguments are passed, first one is assumed to be seed */
      seed = atol(argv[1]);
      i = atoi(argv[2]);
      break;

    default : /* if less than 1 or more than 2 arguments are passed, the executions exits */
      printf("Invalid input for data-gen.c!\n");
      exit(0);
  }
  
  FILE *log;
  log = fopen("log.txt", "w"); /* writes to new file log.txt which holds important parameters */

  fprintf(log, "Seed used: %d \nNumber of particles: %d \n", seed, i);

  fclose(log);

  FILE *output;
  output = fopen("output.csv", "w"); /* writes to new file output.csv which holds positions */

  for(int j = 0; j < i; ++j)
  {
    fprintf(output, "%d, %d, %d\n\n", rand_pos(), rand_pos(), rand_pos());
  }

  fclose(output);
  
  return 0;
}
