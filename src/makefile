nbody: plummer.c mersenne.c hermite.c
	gcc -o nbody start.c plummer.c mersenne.c hermite.c output.c ediag.c -lm

.PHONY : clean
clean:
	-rm nbody nbody.err nbody.out
