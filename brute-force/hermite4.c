#include <stdio.h>
#include <complex.h>

double dt;

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

/* computes acceleration */
double complex * acc(struct body p)
{
  double complex r2[3] = {(p.xpos * p.xpos), (p.ypos * p.ypos), (p.zpos * p.zpos)};
  double complex r3[3] = {(r2[0] * csqrt(r2[0])), (r2[1] * csqrt(r2[1])), (r2[2] * csqrt(r2[2]))};
  double complex res[3] = {(p.xpos * (-p.mass / r3[0])), (p.ypos * (-p.mass / r3[1])), 
                           (p.zpos * (-p.mass / r3[2]))};
  return res;
}

complex double * jerk(struct body p)
{
  double complex r2[3] = {(p.xpos * p.xpos), (p.ypos * p.ypos), (p.zpos * p.zpos)};
  double complex r3[3] = {(r2[1] * csqrt(r2[1])), (r2[2] * csqrt(r2[2])), (r2[3] * csqrt(r2[3]))};
  double complex res[3] = {((p.xvel + p.xpos * (-3 * (p.xpos * p.xvel) / r2[0])) * (-p.mass / r3[0])),
                          ((p.yvel + p.ypos * (-3 * (p.ypos * p.yvel) / r2[1])) * (-p.mass / r3[1])),
                          ((p.zvel + p.zpos * (-3 * (p.zpos * p.zvel) / r2[2])) * (-p.mass / r3[2]))};
  return res;
}

void hermite(struct body p)
{
  /* get current positions and velocities of particle */
  double complex old_pos[3] = {p.xpos, p.ypos, p.zpos};
  double complex old_vel[3] = {p.xvel, p.yvel, p.zvel};
  
  /* get current acceleration and jerk of particle */
  double complex *old_acc = acc(p);
  double complex *old_jerk = jerk(p);
  
  /* computes prediction of new positions based on acceleration and jerk */
  p.xpos += (p.xvel * dt) + old_acc[0] * (dt * dt / 2.0) + old_jerk[0] * (dt * dt * dt / 6.0);
  p.ypos += (p.yvel * dt) + old_acc[1] * (dt * dt / 2.0) + old_jerk[1] * (dt * dt * dt / 6.0);
  p.zpos += (p.zvel * dt) + old_acc[2] * (dt * dt / 2.0) + old_jerk[2] * (dt * dt * dt / 6.0);
  
  /* computes prediction of new velocities based on acceleration and jerk */
  p.xvel += old_acc[0] * dt + old_jerk[0] * (dt * dt / 2.0);
  p.yvel += old_acc[1] * dt + old_jerk[1] * (dt * dt / 2.0);
  p.zvel += old_acc[2] * dt + old_jerk[2] * (dt * dt / 2.0);
  
  double complex *new_acc = acc(p);
  double complex *new_jerk = jerk(p);
  
  /* correction in reversed order of computation */
  p.xvel = old_vel[0] + (old_acc[0] + new_acc[0]) * (dt / 2.0) + (old_jerk[0] - new_jerk[0]) * (dt * dt / 12.0);
  p.yvel = old_vel[1] + (old_acc[1] + new_acc[1]) * (dt / 2.0) + (old_jerk[1] - new_jerk[1]) * (dt * dt / 12.0);
  p.zvel = old_vel[2] + (old_acc[2] + new_acc[2]) * (dt / 2.0) + (old_jerk[2] - new_jerk[2]) * (dt * dt / 12.0);
  
  p.xpos = old_pos[0] + (old_vel[0] + p.xvel) * (dt / 2.0) + (old_acc[0] + new_acc[0]) * (dt * dt / 12.0);
  p.ypos = old_pos[1] + (old_vel[1] + p.yvel) * (dt / 2.0) + (old_acc[1] + new_acc[1]) * (dt * dt / 12.0);
  p.zpos = old_pos[2] + (old_vel[2] + p.zvel) * (dt / 2.0) + (old_acc[2] + new_acc[2]) * (dt * dt / 12.0);    
}

int main(int argc, const char *argv[])
{
  dt = 0.1;
  
  double N = 10;
  double end_time = 1.0;
  double time = 0.0;
  
  struct body bodies[N];
  
  while(time < end_time)
  {
    for(int i = 0; i < N; ++i)
    {
      hermite(bodies[i]);
    }
  }
}
