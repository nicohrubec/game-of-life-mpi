#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// prints a matrix for debugging
void print_matrix(uint8_t(*matrix), int n_rows, int n_cols) {
    int offset;

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            offset = i * n_cols + j;
            printf("%d ", matrix[offset]);
        }
        printf("\n");
    }
    fflush(stdout);
}

// fills a matrix with an initial input
// percentage of alive cells can be configured with the density
void fill_matrix(uint8_t(*matrix), int n_rows, int n_cols, int density) {
    int offset, r;

    for(int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            offset = i * n_cols + j;
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

uint8_t get_neighbor_value(uint8_t(*matrix), int i, int j, int n_rows, int n_cols) {
    int ni = i % n_rows;
    int nj = (j + n_cols) % n_cols;
    int offset = ni * n_cols + nj;

    return matrix[offset];
}

// applies stencil
uint8_t get_num_alive_cells_in_neighborhood(uint8_t(*matrix), int i, int j, int n_rows, int n_cols) {
    uint8_t num_alive_cells = 0;

    num_alive_cells += get_neighbor_value(matrix, i-2, j-2, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i, j-2, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i+2, j-2, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i-1, j, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i+2, j, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i-1, j+1, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i, j+1, n_rows, n_cols);
    num_alive_cells += get_neighbor_value(matrix, i+2, j+2, n_rows, n_cols);

    return num_alive_cells;
}

uint8_t state_lookup[2][9] = {
        { 0, 0, 0, 1, 0, 0, 0, 0, 0 },
        { 0, 0, 1, 1, 0, 0, 0, 0, 0 }
};

// runs the gol for one iteration
void run_generation(uint8_t(*current_generation), uint8_t(*next_generation), int n_rows, int n_cols) {
    int offset, num_alive_neighbors;
    uint8_t new_cell_state;
    uint8_t cell_state;

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            offset = i * n_cols + j;
            cell_state = current_generation[offset];
            num_alive_neighbors = get_num_alive_cells_in_neighborhood(current_generation, i, j, n_rows, n_cols);
            new_cell_state = state_lookup[cell_state][num_alive_neighbors];
            next_generation[offset] = new_cell_state;
        }
    }
}

// prints the number of dead and alive cells
void print_summary_output(uint8_t(*matrix), int n_rows, int n_cols, int c_generation) {
    int num_alive_cells = 0;

    for (int i = 0; i < n_rows * n_cols; i++) {
        num_alive_cells += matrix[i];
    }

    printf("\n\nOutput after generation %d:\n", c_generation);
    printf("Number of alive cells: %d\n", num_alive_cells);
    printf("Number of dead cells: %d\n", n_rows * n_cols - num_alive_cells);
}

// copies the values from matrix 2 to matrix 1
void copy_matrix(uint8_t(*matrix1), uint8_t(*matrix2), int n_rows, int n_cols) {
    for (int i = 0; i < n_rows * n_cols; i++) {
        matrix1[i] = matrix2[i];
    }
}



int main(int argc, char *argv[]) {
    int n_rows = 10, n_cols = 10, n_generations = 2; // num rows, num cols, num generations
    int seed = 42;
    int verbose = 0;
    int density = 27; // in percent
    int opt;

    // define arguments
    static struct option long_options[] = {{"num_rows", required_argument, 0, 'r'},
                                           {"num_cols", required_argument, 0, 'c'},
                                           {"num_generations", required_argument, 0, 'g'},
                                           {"seed", required_argument, 0, 's'},
                                           {"density", required_argument, 0, 'd'},
                                           {"verbose", no_argument, 0, 'v'},
                                           {0, 0, 0, 0}};

    // argument parsing
    while (1) {
        int option_index = 0;
        opt = getopt_long(argc, argv, "r:c:g:s:d:v", long_options, &option_index);

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
            default:
                fprintf(stderr, "Usage: %s -r <r> -c <c> -g <g> -s <s> -v\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // density bounds checking
    if( density <= 0 || density > 100) {
        fprintf(stderr, "Density should be between 1 and 100\n");
        exit(EXIT_FAILURE);
    }

    srand(seed); // guarantee reproducible results

    if (verbose) {
        printf("Parameters:");
        printf("seed: %d\n", seed);
        printf("n_rows: %d\n", n_rows);
        printf("n_cols: %d\n", n_cols);
        printf("n_generations: %d\n", n_generations);
        printf("density: %d \n", density);
    }

    // allocate matrices
    uint8_t *current_generation = (uint8_t *)malloc(n_rows * n_cols * sizeof(uint8_t));
    uint8_t *next_generation = (uint8_t *)malloc(n_rows * n_cols * sizeof(uint8_t));

    // generate initial input
    fill_matrix(current_generation, n_rows, n_cols, density);

    // run gol
    for (int c_generation = 0; c_generation < n_generations; c_generation++) {
        if (verbose) {
            print_summary_output(current_generation, n_rows, n_cols, c_generation);
            print_matrix(current_generation, n_rows, n_cols);
        }

        // todo: get timings
        run_generation(current_generation, next_generation, n_rows, n_cols);
        copy_matrix(current_generation, next_generation, n_rows, n_cols);
    }

    // final output
    print_summary_output(next_generation, n_rows, n_cols, n_generations);
    if (verbose) {
        print_matrix(current_generation, n_rows, n_cols);
    }

    // resource clean up
    free(current_generation);
    free(next_generation);

    return 0;
}