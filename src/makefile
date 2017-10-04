nbody: plummer.c mersenne.c hermite.c
	gcc -o nbody start.c plummer.c mersenne.c hermite.c output.c ediag.c -lm -m64 -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes

.PHONY : clean
clean:
	-rm nbody nbody.err nbody.out
