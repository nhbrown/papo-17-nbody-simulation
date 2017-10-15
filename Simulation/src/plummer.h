#ifndef PLUMMER_H_
#define PLUMMER_H_

double rrand(double low, double high);

void plummer(int N, double *mass, double complex *pos, double complex *vel, int i, double M, double R);

void center_of_mass_adjustment(int N, double *mass, double complex *pos, double complex *vel);

void startPlummer(unsigned long s, int N, double *mass, double complex *pos, double complex *vel, 
                  double M, double R);

#endif // PLUMMER_H_
