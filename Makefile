seq: seq.c
	mpicc -Wall --std=c99 -o seq seq.c

clean_seq:
	rm -f seq

par: par.c
	mpicc -Wall --std=c99 -o par par.c

clean_par:
	rm -f par

par_collective: par_collective.c
	mpicc -Wall --std=c99 -o par_collective par_collective.c

clean_par_collective:
	rm -f par_collective

mpi_rand: mpi_rand.c
	mpicc -Wall --std=c99  -o $@ $<

clean_mpi_rand:
	rm -f mpi_rand

clean:
	rm -f seq par mpi_rand