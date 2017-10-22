nbody: driver.c plummer.c mersenne.c hermite.c output.c ediag.c
	mpicc -g -o nbody driver.c plummer.c mersenne.c hermite.c output.c ediag.c -lm -m64 -pedantic -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wmaybe-uninitialized

.PHONY : clean
clean:
	-rm nbody nbody.err nbody.out
  
