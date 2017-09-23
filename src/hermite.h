#ifndef HERMITE_H_
#define HERMITE_H_

void startHermite(int particles, double timestep, double end_time, double *mass, double complex **pos, double complex **vel, 
                  double complex **acc, double complex **jerk);

#endif // HERMITE_H_
