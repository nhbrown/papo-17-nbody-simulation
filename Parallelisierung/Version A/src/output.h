#ifndef OUTPUT_H_
#define OUTPUT_H_

void createNames(void);

void printInitialConditions(int N, int DIM, double *mass, double complex *pos, double complex *vel);

void printLog(unsigned long seed, int N, double M, double R, double G, double timestep, double end_time);

void printEnergyDiagnostics(double e_kinetic, double e_potential, double e_total);

void printIteration(int N, int DIM, int iteration, double *mass, double complex *pos, double complex *vel);

#endif // OUTPUT_H_
