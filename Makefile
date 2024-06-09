seq: seq.c
	mpicc -O3 -Wall --std=c99 -o seq seq.c

clean_seq:
	rm -f seq

par: par.c
	mpicc -O3 -Wall --std=c99 -o par par.c

clean_par:
	rm -f par

par_collective: par_collective.c
	mpicc -O3 -Wall --std=c99 -o par_collective par_collective.c

clean_par_collective:
	rm -f par_collective

clean:
	rm -f seq par par_collective