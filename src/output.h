#ifndef OUTPUT_H_
#define OUTPUT_H_

void createNames();

void printInitialConditions(int N, double *mass, double complex **pos, double complex **vel);

void printLog(unsigned long seed, int N, double M, double R, double G, double timestep, double end_time);

void printIteration(double *mass, double complex **pos, double complex **vel, int iteration, int N);

#endif // OUTPUT_H_
