#ifndef OUTPUT_H_
#define OUTPUT_H_

void printIteration(double *mass, double complex **pos, double complex **vel, int iteration, char *foldername);

void printInitialConditions(unsigned long seed, int N, double M, double R, double G, 
                            double timestep, double end_time, double *mass, double complex **pos, double complex **vel);

void createNames();

#endif // OUTPUT_H_
