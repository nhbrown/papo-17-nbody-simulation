#ifndef MERSENNE_H_
#define MERSENNE_H_

void init_genrand(unsigned long s);

double genrand_real1(void);

void init_by_array(unsigned long init_key[], int key_length);

unsigned long genrand_int32(void);

long genrand_int31(void);

double genrand_real2(void);

double genrand_real3(void);

double genrand_res53(void);

#endif // MERSENNE_H_
