#ifndef OUTPUT_H_
#define OUTPUT_H_

void createNames(void);

void printInitialConditions(int N, double *mass, double complex **pos, double complex **vel);

void printLog(unsigned long seed, int N, double M, double R, double G, double timestep, double end_time);

void printEnergyDiagnostics(int marker, double e_kinetic, double e_potential, double e_total);

void printIteration(double *mass, double complex **pos, double complex **vel, int iteration, int N);

#endif // OUTPUT_H_
