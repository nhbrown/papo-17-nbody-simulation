nbody: driver.c plummer.c mersenne.c hermite.c output.c ediag.c
	mpicc -lm -o nbody driver.c plummer.c mersenne.c hermite.c output.c ediag.c -O3 -Wall -Wextra
.PHONY : clean
clean:
	-rm nbody nbody.err nbody.out
  
