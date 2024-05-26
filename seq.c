#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>

int get_offset(int i, int j, int n) {
    return i * n + j;
}

// prints a matrix for debugging
void print_matrix(uint8_t(*matrix), int n) {
    int offset;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            offset = get_offset(i, j, n);
            printf("%d ", matrix[offset]);
        }
        printf("\n");
    }
    fflush(stdout);
}

// fills a matrix with an initial input
// percentage of alive cells can be configured with the density
void fill_matrix(uint8_t(*matrix), int n, int density) {
    int offset, r;

    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            offset = get_offset(i, j, n);
            r = rand() % 100;

            if (r < density) {
                r = 1;
            } else {
                r = 0;
            }

            matrix[offset] = r;
        }
    }
}

uint8_t get_neighbor_value(uint8_t(*matrix), int i, int j, int n) {
    return matrix[get_offset(i, j, n)];
}

int modulo(int x, int y) {
    return x - (y * (x / y));
}

int stencil_plus_operator(int x, int d, int m) {
    return modulo(x+d, m);
}

int stencil_minus_operator(int x, int d, int m) {
    return modulo(x-d+m, m);
}

// applies stencil
uint8_t get_num_alive_cells_in_neighborhood(uint8_t(*matrix), int i, int j, int n) {
    uint8_t num_alive_cells = 0;

    num_alive_cells += get_neighbor_value(matrix, stencil_minus_operator(i, 2, n), stencil_minus_operator(j, 2, n), n);
    num_alive_cells += get_neighbor_value(matrix, i, stencil_minus_operator(j, 2, n), n);
    num_alive_cells += get_neighbor_value(matrix, stencil_plus_operator(i, 2, n), stencil_minus_operator(j, 2, n), n);
    num_alive_cells += get_neighbor_value(matrix, stencil_minus_operator(i, 1, n), j, n);
    num_alive_cells += get_neighbor_value(matrix, stencil_plus_operator(i, 2, n), j, n);
    num_alive_cells += get_neighbor_value(matrix, stencil_minus_operator(i, 1, n), stencil_plus_operator(j, 1, n), n);
    num_alive_cells += get_neighbor_value(matrix, i, stencil_plus_operator(j, 1, n), n);
    num_alive_cells += get_neighbor_value(matrix, stencil_plus_operator(i, 2, n), stencil_plus_operator(j, 2, n), n);

    return num_alive_cells;
}

uint8_t state_lookup[2][9] = {
        { 0, 0, 0, 1, 0, 0, 0, 0, 0 },
        { 0, 0, 1, 1, 0, 0, 0, 0, 0 }
};

// runs the gol for one iteration
void run_generation(uint8_t(*current_generation), uint8_t(*next_generation), int n) {
    int offset, num_alive_neighbors;
    uint8_t new_cell_state;
    uint8_t cell_state;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            offset = get_offset(i, j, n);;
            cell_state = current_generation[offset];
            num_alive_neighbors = get_num_alive_cells_in_neighborhood(current_generation, i, j, n);
            new_cell_state = state_lookup[cell_state][num_alive_neighbors];
            next_generation[offset] = new_cell_state;
        }
    }
}

// prints the number of dead and alive cells
void print_summary_output(uint8_t(*matrix), int n, int c_generation) {
    int num_alive_cells = 0;

    for (int i = 0; i < n * n; i++) {
        num_alive_cells += matrix[i];
    }

    printf("\n\nOutput after generation %d:\n", c_generation);
    printf("Number of alive cells: %d\n", num_alive_cells);
    printf("Number of dead cells: %d\n", n * n - num_alive_cells);
}

// copies the values from matrix 2 to matrix 1
void copy_matrix(uint8_t(*matrix1), uint8_t(*matrix2), int n) {
    for (int i = 0; i < n * n; i++) {
        matrix1[i] = matrix2[i];
    }
}

int main(int argc, char *argv[]) {
    int n = 10, n_generations = 2; // num rows, num cols, num generations
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
        printf("Parameters:");
        printf("seed: %d\n", seed);
        printf("number: %d\n", n);
        printf("n_generations: %d\n", n_generations);
        printf("density: %d \n", density);
    }

    // allocate matrices
    uint8_t *current_generation = (uint8_t *)malloc(n * n * sizeof(uint8_t));
    uint8_t *next_generation = (uint8_t *)malloc(n * n * sizeof(uint8_t));

    // generate initial input
    fill_matrix(current_generation, n, density);

    if (verbose) {
        print_summary_output(current_generation, n, 0);
        print_matrix(current_generation, n);
    }

    // run gol
    for (int c_generation = 1; c_generation <= n_generations; c_generation++) {
        start_time = MPI_Wtime();
        run_generation(current_generation, next_generation, n);
        copy_matrix(current_generation, next_generation, n);
        end_time = MPI_Wtime();

        total_time_generation = (end_time - start_time) * 1e6; // μs
        total_time += total_time_generation;

        if (verbose) {
            print_summary_output(current_generation, n, c_generation);
            print_matrix(current_generation, n);
            printf("Time needed for generation: %f\n", total_time_generation);
        }
    }

    printf("Total time: %f\n", total_time);
    printf("Average time/generation: %f\n", total_time / n_generations);

    // resource clean up
    free(current_generation);
    free(next_generation);

    MPI_Finalize();

    return 0;
}