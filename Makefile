seq: seq.c
	clang -Wall --std=c99 -o seq seq.c

clean_seq:
	rm -f seq

mpi_rand: mpi_rand.c
	mpicc  -Wall --std=c99  -o $@ $<

clean_mpi_rand:
	rm -f mpi_rand