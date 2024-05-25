#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int n_rows = 10, n_cols = 10, n_generations = 2; // num rows, num cols, num generations
    int seed = 42;
    int verbose = 0;
    int verify = 0; // if set we perform verification with the sequential version
    int weak_scaling = 0; // if set input size is per processor
    int density = 27; // in percent
    int opt;
    int rank, size;

    // define arguments
    static struct option long_options[] = {{"num_rows", required_argument, 0, 'r'},
                                           {"num_cols", required_argument, 0, 'c'},
                                           {"num_generations", required_argument, 0, 'g'},
                                           {"seed", required_argument, 0, 's'},
                                           {"density", required_argument, 0, 'd'},
                                           {"verbose", no_argument, 0, 'v'},
                                           {"verify", no_argument, 0, 'x'},
                                           {"weak_scaling", no_argument, 0, 'w'},
                                           {0, 0, 0, 0}};

    MPI_Init(&argc, &argv);

    // argument parsing
    while (1) {
        int option_index = 0;
        opt = getopt_long(argc, argv, "r:c:g:s:d:vxw", long_options, &option_index);

        if (opt == -1) break;

        switch (opt) {
            case 'r':
                n_rows = atoi(optarg);
                break;
            case 'c':
                n_cols = atoi(optarg);
                break;
            case 'g':
                n_generations = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'd':
                density = atoi(optarg);
                break;
            case 'v':
                verbose = 1;
                break;
            case 'x':
                verify = 1;
                break;
            case 'w':
                weak_scaling = 1;
                break;
            default:
                fprintf(stderr, "Usage: %s -r <r> -c <c> -g <g> -s <s> -v -x\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // density bounds checking
    if( density <= 0 || density > 100) {
        fprintf(stderr, "Density should be between 1 and 100\n");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (verbose) {
        printf("World rank is %d\n", rank);
        printf("World size is %d\n", size);
    }

    srand(seed); // guarantee reproducible results

    if (verbose) {
        printf("Parameters:\n");
        printf("seed: %d\n", seed);
        printf("n_rows: %d\n", n_rows);
        printf("n_cols: %d\n", n_cols);
        printf("n_generations: %d\n", n_generations);
        printf("density: %d \n", density);

        if (verify) printf("verification is on\n");
        if (weak_scaling) printf("weak scaling is on\n");
    }

    // set up the communicator

    // generate input already distributed on different processors
    // then two step process:
    // 1. communication
    // 2. local stencil update

    // how to compare with sequential version?
    // no idea currently, shared memory?

    // how to distribute the input (submatrices or subvectors)?
    // only submatrices make sense probably because otherwise communication cost is a lot

    MPI_Finalize();

    return 0;
}