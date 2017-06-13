#include <stdio.h>
#include <stdlib.h>
#include <time.h>

long seed;
int i;

const m = 2147483647; /* modulos */
const a = 48271;  /* multiplier */
const q = 44488; /* (m / a) */
const r = 3399; /* (m % a) */

int rand_pos()
{
  if(seed <= 0)
  {
    seed += m;
  }
  
  seed = a * (seed % q) - r * (seed / q);
  
  return seed;
}

void test()
{
  printf("The seed is %d and the number of particles is %d.\n", seed, i);
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
  
  test();
}
