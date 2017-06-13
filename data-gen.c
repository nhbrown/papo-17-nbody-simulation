#include stdio.h

const long seed;
const int i;

int main(int argc, constant char *argv[])
{
  switch(argc)
  {
    case 2 : 
      seed = (unsigned long)time(NULL);
      i = argv[1];
      break;
      
    case 3 :
      seed = argv[1];
      i = argv[2];
      break;
      
    default : 
      printf("Invalid input for data-gen.c");
  }
}
