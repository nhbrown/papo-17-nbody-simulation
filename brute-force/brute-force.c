#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void update(double dt, struct body x)
{
  x.xvel += dt * x.fx / x.mass;
  x.yvel += dt * x.fy / x.mass;
  x.zvel += dt * x.fz / x.mass;
  
  x.xpos += dt * x.xvel;
  x.ypos += dt * x.yvel;
  x.zpos += dt * x.zvel;
}

double distanceBetween(struct body a, struct body b)
{
  double dx = a.xpos - b.xpos;
  double dy = a.ypos - b.ypos;
  double dz = a.zpos - b.zpos;
  
  return csqrt(dx * dx + dy * dy + dz * dz);
}

void resetForces(struct body x)
{
  x.fx = 0.0;
  x.fy = 0.0;
  x.fz = 0.0;
}

void addForces(struct body a, struct body b)
{
  double distance = distanceBetween(a, b);
  double F = (1.0 * a.mass * b.mass) / (dist * dist);
  
  a.fx = F * dx / distance;
  a.fy = F * dy / distance;
}

int main(int argc, const char *argv[])
{
  for(int i = 0; i < N; ++i)
  {
    struct body a = bodies[i];
    resetForces(a);
    
    for(int j = 0; i < N; ++i)
    {
      if(i != j)
      {
        struct body b = bodies[j];
        addForces(a, b);
      }
    }
    
    for(int i = 0; i < N; ++i)
    {
      struct body p = bodies[i];
      update(0.98, p);
    }
  }
}
