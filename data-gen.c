#include <stdio.h>
#include <stdlib.h>
#include <time.h>

long seed;
int i;

const int m = 2147483647; /* modulos */
const int a = 48271;  /* multiplier */
const int q = 44488; /* (m / a) */
const int r = 3399; /* (m % a) */

int rand_pos()
{
  if(seed <= 0)
  {
    seed += m;
  }
  
  seed = a * (seed % q) - r * (seed / q);
  
  return seed;
}

int main(int argc, const char *argv[])
{
  switch(argc)
  {
    case 2 :
      seed = (unsigned long)time(NULL);
      i = atoi(argv[1]);
      break;

    case 3 :
      seed = atol(argv[1]);
      i = atoi(argv[2]);
      break;

    default :
      printf("Invalid input for data-gen.c!\n");
      exit(0);
  }
  
  printf("The seed is %d and the number of particles is %d.\n", seed, i);
  
  for(int j = 0; j < i; ++j)
  {
    printf("Position #%d is: %d \n", j + 1, rand_pos());
  }
}
