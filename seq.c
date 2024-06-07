#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>

// prints a matrix for debugging
void print_matrix(int n, uint8_t(*matrix)[n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
    fflush(stdout);
}

// fills a matrix with an initial input
// percentage of alive cells can be configured with the density
void fill_matrix(int n, uint8_t(*matrix)[n], int density) {
    int r;

    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            r = rand() % 100;

            if (r < density) {
                r = 1;
            } else {
                r = 0;
            }

            matrix[i][j] = r;
        }
    }
}

int modulo(int x, int y) {
    return x & (y - 1); // we can use this because y is always a power of 2
}

int stencil_plus_operator(int x, int d, int m) {
    return modulo(x+d, m);
}

int stencil_minus_operator(int x, int d, int m) {
    return modulo(x-d+m, m);
}

// applies stencil
uint8_t get_num_alive_cells_in_neighborhood(int n, uint8_t(*matrix)[n], int i, int j) {
    uint8_t num_alive_cells = 0;

    // precompute indices
    const int i_minus_2_idx = stencil_minus_operator(i, 2, n);
    const int j_minus_2_idx = stencil_minus_operator(j, 2, n);
    const int i_minus_1_idx = stencil_minus_operator(i, 1, n);
    const int j_plus_1_idx = stencil_plus_operator(j, 1, n);
    const int i_plus_2_idx = stencil_plus_operator(i, 2, n);
    const int j_plus_2_idx = stencil_plus_operator(j, 2, n);

    num_alive_cells += matrix[i_minus_2_idx][j_minus_2_idx];
    num_alive_cells += matrix[i][j_minus_2_idx];
    num_alive_cells += matrix[i_plus_2_idx][j_minus_2_idx];
    num_alive_cells += matrix[i_minus_1_idx][j];
    num_alive_cells += matrix[i_plus_2_idx][j];
    num_alive_cells += matrix[i_minus_1_idx][j_plus_1_idx];
    num_alive_cells += matrix[i][j_plus_1_idx];
    num_alive_cells += matrix[i_plus_2_idx][j_plus_2_idx];

    return num_alive_cells;
}

uint8_t state_lookup[2][9] = {
        { 0, 0, 0, 1, 0, 0, 0, 0, 0 },
        { 0, 0, 1, 1, 0, 0, 0, 0, 0 }
};

// runs the gol for one iteration
void run_generation(int n, uint8_t(*current_generation)[n], uint8_t(*next_generation)[n]) {
    int num_alive_neighbors;
    uint8_t new_cell_state;
    uint8_t cell_state;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cell_state = current_generation[i][j];
            num_alive_neighbors = get_num_alive_cells_in_neighborhood(n, current_generation, i, j);
            new_cell_state = state_lookup[cell_state][num_alive_neighbors];
            next_generation[i][j] = new_cell_state;
        }
    }
}

// prints the number of dead and alive cells
void print_summary_output(int n, uint8_t(*matrix)[n], int c_generation) {
    int num_alive_cells = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            num_alive_cells += matrix[i][j];
        }
    }

    printf("\n\nOutput after generation %d:\n", c_generation);
    printf("Number of alive cells: %d\n", num_alive_cells);
    printf("Number of dead cells: %d\n", n * n - num_alive_cells);
}

// copies the values from matrix 2 to matrix 1
void copy_matrix(int n, uint8_t(*matrix1)[n], uint8_t(*matrix2)[n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix1[i][j] = matrix2[i][j];
        }
    }
}

int main(int argc, char *argv[]) {
    int n = 10, n_generations = 1000; // num rows, num cols, num generations
    int seed = 42;
    int verbose = 0;
    int density = 28; // in percent
    int opt;
    int rank, size;
    double total_time = 0.0, total_time_generation, start_time, end_time;

    // define arguments
    static struct option long_options[] = {{"number", required_argument, 0, 'n'},
                                           {"num_generations", required_argument, 0, 'g'},
                                           {"seed", required_argument, 0, 's'},
                                           {"density", required_argument, 0, 'd'},
                                           {"verbose", no_argument, 0, 'v'},
                                           {0, 0, 0, 0}};

    MPI_Init(&argc, &argv);

    // argument parsing
    while (1) {
        int option_index = 0;
        opt = getopt_long(argc, argv, "n:g:s:d:v", long_options, &option_index);

        if (opt == -1) break;

        switch (opt) {
            case 'n':
                n = atoi(optarg);
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
            default:
                fprintf(stderr, "Usage: %s -n <n> -g <g> -s <s> -v\n", argv[0]);
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
        printf("number: %d\n", n);
        printf("n_generations: %d\n", n_generations);
        printf("density: %d \n", density);
    }

    // allocate matrices
    uint8_t(*current_generation)[n];
    current_generation = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));
    uint8_t(*next_generation)[n];
    next_generation = (uint8_t(*)[n])malloc(n * n * sizeof(uint8_t));

    // generate initial input
    fill_matrix(n, current_generation, density);

    if (verbose) {
        print_summary_output(n, current_generation, 0);
        print_matrix(n, current_generation);
    }

    // run gol
    for (int c_generation = 1; c_generation <= n_generations; c_generation++) {
        start_time = MPI_Wtime();
        run_generation(n, current_generation, next_generation);
        copy_matrix(n, current_generation, next_generation);
        end_time = MPI_Wtime();

        total_time_generation = (end_time - start_time) * 1e6; // μs
        total_time += total_time_generation;

        if (verbose) {
            print_summary_output(n, current_generation, c_generation);
            print_matrix(n, current_generation);
            printf("Time needed for generation: %f\n", total_time_generation);
        }
    }

    printf("Total time: %.2f μs\n", total_time);
    printf("Average time/generation: %f\n", total_time / n_generations);

    // resource clean up
    free(current_generation);
    free(next_generation);

    MPI_Finalize();

    return 0;
}