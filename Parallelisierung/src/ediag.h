#ifndef EDIAG_H_
#define EDIAG_H_

void kinetic_energy(int N, int DIM, double *mass, double complex *vel);

void potential_energy(int N, int DIM, double *mass, double complex *pos);

void energy_diagnostics(int N, int DIM, double *mass, double complex *pos, double complex *vel);

#endif // EDIAG_H_
