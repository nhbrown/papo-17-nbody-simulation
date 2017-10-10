#ifndef HERMITE_H_
#define HERMITE_H_

void acc_jerk(int N, double *mass, double complex **pos, double complex **vel, 
              double complex **acc, double complex **jerk);

void hermite(int N, double dt, double *mass, double complex **pos, 
             double complex **vel, double complex **acc, double complex **jerk);

void startHermite(int N, double dt, double end_time, double *mass, double complex **pos, double complex **vel, 
                  double complex **acc, double complex **jerk);

#endif // HERMITE_H_
