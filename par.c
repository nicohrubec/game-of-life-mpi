#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#include <string.h>

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
            offset = get_offset(i, j, n);
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

int compare_matrices(int n_loc_r, int n_loc_c, uint8_t(*sequential_matrix), uint8_t(*local_matrix)[n_loc_c], int n, int m_offset_r, int m_offset_c) {
    int offset;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i >= m_offset_r && i < m_offset_r + n_loc_r &&
               j >= m_offset_c && j < m_offset_c + n_loc_c) {
                offset = get_offset(i, j, n);
                if (sequential_matrix[offset] != local_matrix[i-m_offset_r][j-m_offset_c]) {
                    printf("Verification failed at index %d %d\n", i, j);
                    return 0;
                }
            }
        }
    }

    return 1;
}

void copy_full_matrix_to_local_matrix(int n_loc_r, int n_loc_c, uint8_t(*sequential_matrix), uint8_t(*local_matrix)[n_loc_c], int n, int m_offset_r, int m_offset_c) {
    int offset;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i >= m_offset_r && i < m_offset_r + n_loc_r &&
               j >= m_offset_c && j < m_offset_c + n_loc_c) {
                offset = get_offset(i, j, n);
                local_matrix[i-m_offset_r][j-m_offset_c] = sequential_matrix[offset];
            }
        }
    }
}

void fill_matrix_par(int n_loc_r, int n_loc_c, uint8_t(*matrix)[n_loc_c], int n, int density, int m_offset_r, int m_offset_c) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int r = rand() % 100;
            if( r < density) {
                r = 1;
            } else {
                r = 0;
            }

            if(i >= m_offset_r && i < m_offset_r + n_loc_r &&
               j >= m_offset_c && j < m_offset_c + n_loc_c) {
                matrix[i - m_offset_r][j - m_offset_c] = r;
            }
        }
    }
}

void print_matrix_par(int n_loc_r, int n_loc_c, uint8_t(*matrix)[n_loc_c], int rank, int size, int c_generation) {
    for (int i = 0; i < size; i++) {
        if (rank == i) {
            printf("%d: local matrix at generation %d\n", rank, c_generation);
            for (int i = 0; i < n_loc_r; i++) {
                for (int j = 0; j < n_loc_c; j++) {
                    printf("%d ", matrix[i][j]);
                }
                printf("\n");
            }
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void copy_matrix_par(int n_loc_r, int n_loc_c, uint8_t(*matrix1)[n_loc_c], uint8_t(*matrix2)[n_loc_c], int rank, int size) {
    for (int i = 0; i < size; i++) {
        if (rank == i) {
            for (int i = 0; i < n_loc_r; i++) {
                for (int j = 0; j < n_loc_c; j++) {
                    matrix1[i][j] = matrix2[i][j];
                }
            }
        }
    }
}

uint8_t get_num_alive_cells_in_par_neighborhood(int n_loc_c, uint8_t(*matrix)[n_loc_c], int i, int j, int rank) {
    uint8_t num_alive_cells = 0;

    num_alive_cells += matrix[i-2][j-2];
    // if (rank == 0) printf("i - 2, j - 2: %d\n", matrix[i-2][j-2]);
    num_alive_cells += matrix[i][j-2];
    // if (rank == 0) printf("i, j - 2: %d\n", matrix[i][j-2]);
    num_alive_cells += matrix[i+2][j-2];
    // if (rank == 0) printf("i + 2, j - 2: %d\n", matrix[i+2][j-2]);
    num_alive_cells += matrix[i-1][j];
    // if (rank == 0) printf("i - 1, j: %d\n", matrix[i-1][j]);
    num_alive_cells += matrix[i+2][j];
    // if (rank == 0) printf("i + 2, j: %d\n", matrix[i+2][j]);
    num_alive_cells += matrix[i-1][j+1];
    // if (rank == 0) printf("i - 1, j + 1: %d\n", matrix[i-1][j+1]);
    num_alive_cells += matrix[i][j+1];
    // if (rank == 0) printf("i, j + 1: %d\n", matrix[i][j+1]);
    num_alive_cells += matrix[i+2][j+2];
    // if (rank == 0) printf("i + 2, j + 2: %d\n", matrix[i+2][j+2]);

    // if (rank == 0) printf("i: %d, j: %d, num alive cells: %d\n", i-2, j-2, num_alive_cells);
    // if (rank == 0) printf("local matrix at i: %d, j: %d is %d\n", i-2, j-2, matrix[i][j]);
    return num_alive_cells;
}

void run_generation_par(int n_loc_r, int n_loc_c, uint8_t(*current_generation)[n_loc_c+4], uint8_t(*next_generation)[n_loc_c], int rank) {
    int num_alive_neighbors;
    uint8_t new_cell_state;
    uint8_t cell_state;

    for (int i = 0; i < n_loc_r; i++) {
        for (int j = 0; j < n_loc_c; j++) {
            cell_state = current_generation[i+2][j+2];
            num_alive_neighbors = get_num_alive_cells_in_par_neighborhood(n_loc_c+4, current_generation, i+2, j+2, rank);
            new_cell_state = state_lookup[cell_state][num_alive_neighbors];
            next_generation[i][j] = new_cell_state;
        }
    }
}

int main(int argc, char *argv[]) {
    int n = 10, n_generations = 2; // num rows, num cols, num generations
    int seed = 42;
    int verbose = 0;
    int verify = 0; // if set we perform verification with the sequential version
    int weak_scaling = 0; // if set input size is per processor
    int density = 28; // in percent
    int opt;
    int rank, size;

    int n_loc_r, n_loc_c;
    int nprows, npcols;
    int prow_idx, pcol_idx;
    int verification_result;

    // cartesian communicator
    int dims[2] = {0, 0}; // set to 0 to dynamically set dimensions based on number of processes
    int pers[2] = {1, 1}; // grid is not periodic
    int coords[2]; // coordinates of current process in the grid
    MPI_Comm cartcomm;
    MPI_Comm cartcomm_reorder;

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
        opt = getopt_long(argc, argv, "n:g:s:d:vxw", long_options, &option_index);

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
            case 'x':
                verify = 1;
                break;
            case 'w':
                weak_scaling = 1;
                break;
            default:
                fprintf(stderr, "Usage: %s -n <n> -g <g> -s <s> -v -x\n", argv[0]);
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

    srand(seed); // guarantee reproducible results
    MPI_Dims_create(size, 2, dims); // get dimensions of process grid

    if (verbose && rank == 0) {
        printf("Parameters:\n");
        printf("seed: %d\n", seed);
        printf("number: %d\n", n);
        printf("n_generations: %d\n", n_generations);
        printf("density: %d \n", density);
        printf("Dimensions created: [%d, %d]\n", dims[0], dims[1]);

        if (verify) printf("verification is on\n");
        if (weak_scaling) printf("weak scaling is on\n");
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, 0, &cartcomm); // create process grid without reordering
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, 1, &cartcomm_reorder); // process grid with reordering

    // check if process reordering happened
    int result;
    MPI_Comm_compare(cartcomm, cartcomm_reorder, &result);
    if (verbose && rank == 0) {
        if (result == MPI_IDENT || result == MPI_CONGRUENT) {
            printf("No process reordering took place!\n");
        } else {
            printf("Processes were reordered!\n");
        }
    }

    MPI_Cart_coords(cartcomm, rank, 2, coords); // get process coordinates in created grid

    prow_idx = coords[0];
    pcol_idx = coords[1];

    nprows = dims[0];
    npcols = dims[1];

    if (verbose) {
        printf("Process coordinates for rank %d: [%d, %d]\n", rank, prow_idx, pcol_idx);
    }

    if( n % nprows != 0 || n % npcols != 0) {
        if( rank == 0 ) {
            fprintf(stderr, "n should be divisible by nprows and npcols\n");
        }
        exit(EXIT_FAILURE);
    }

    // local matrix size
    n_loc_r = n / nprows;
    n_loc_c = n / npcols;
    if (verbose && rank == 0) {
        printf("n_loc_r: %d n_loc_c: %d\n", n_loc_r, n_loc_c);
    }

    // allocate local matrices for parallel computation
    uint8_t(*current_generation_loc)[n_loc_c];
    current_generation_loc = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));
    uint8_t(*next_generation_loc)[n_loc_c];
    next_generation_loc = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(uint8_t));

    // allocate matrices for sequential verification
    uint8_t *current_generation_seq = NULL;
    uint8_t *next_generation_seq = NULL;

    if (verify) {
        current_generation_seq = (uint8_t *)malloc(n * n * sizeof(uint8_t));
        next_generation_seq = (uint8_t *)malloc(n * n * sizeof(uint8_t));
    }

    // get offset of local matrix in global matrix
    int m_offset_r = prow_idx * n_loc_r;
    int m_offset_c = pcol_idx * n_loc_c;
    if( verbose ) {
        printf("%d: prow_idx: %d pcol_idx: %d m_offset_r: %d m_offset_c: %d\n", rank, prow_idx, pcol_idx, m_offset_r, m_offset_c);
    }

    // input generation
    if (verify) {
        // get sequential input
        fill_matrix(current_generation_seq, n, density);

        if (rank == 0 && verbose) {
            printf("Sequential input matrix: \n");
            print_matrix(current_generation_seq, n);
        }

        // copy values into local matrix
        copy_full_matrix_to_local_matrix(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);
    } else {
        // fill local matrix with initial input
        fill_matrix_par(n_loc_r, n_loc_c, current_generation_loc, n, density, m_offset_r, m_offset_c);
    }

    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // verify initial input
    if (verify && verbose) {
        print_matrix_par(n_loc_r, n_loc_c, current_generation_loc, rank, size, 0);
        verification_result = compare_matrices(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);

        if (verification_result) {
            printf("Input verification on rank %d was successful\n", rank);
        } else {
            printf("Input verification on rank %d failed\n", rank);
        }

        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // communication buffers
    uint8_t *top_rows_recv = malloc(2 * n_loc_c * sizeof(uint8_t));
    uint8_t *bottom_rows_recv = malloc(2 * n_loc_c * sizeof(uint8_t));
    uint8_t *left_cols_recv = malloc(2 * (n_loc_r + 4) * sizeof(uint8_t));
    uint8_t *right_cols_recv = malloc(2 * (n_loc_r + 4) * sizeof(uint8_t));
    uint8_t *left_col = malloc(2 * (n_loc_r + 4) * sizeof(uint8_t));
    uint8_t *right_col = malloc(2 * (n_loc_r + 4) * sizeof(uint8_t));

    // intermediate result buffer after top and bottom rows are communicated
    uint8_t(*loc_matrix_and_neighbor_top_bottom_rows)[n_loc_c];
    loc_matrix_and_neighbor_top_bottom_rows = (uint8_t(*)[n_loc_c])malloc((n_loc_r + 4) * n_loc_c * sizeof(uint8_t));

    // final intermediate result buffer with all neighbors
    uint8_t(*full_current_generation_loc)[n_loc_c + 4];
    full_current_generation_loc = (uint8_t(*)[n_loc_c + 4])malloc((n_loc_r + 4) * (n_loc_c + 4) * sizeof(uint8_t));

    // get process neighbors
    int left_proc_neighbor, right_proc_neighbor, top_proc_neighbor, bottom_proc_neighbor;
    MPI_Cart_shift(cartcomm, 1, 1, &left_proc_neighbor, &right_proc_neighbor);
    MPI_Cart_shift(cartcomm, 0, 1, &top_proc_neighbor, &bottom_proc_neighbor);

    // run gol
    for (int c_generation = 1; c_generation <= n_generations; c_generation++) {
        // communicate top and bottom rows with neighboring processes
        MPI_Sendrecv(current_generation_loc[0], 2 * n_loc_c, MPI_UINT8_T, top_proc_neighbor, 0,
                     bottom_rows_recv, 2 * n_loc_c, MPI_UINT8_T, bottom_proc_neighbor, 0, cartcomm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(current_generation_loc[n_loc_r - 2], 2 * n_loc_c, MPI_UINT8_T, bottom_proc_neighbor, 0,
                     top_rows_recv, 2 * n_loc_c, MPI_UINT8_T, top_proc_neighbor, 0, cartcomm, MPI_STATUS_IGNORE);

        // print_matrix_par(n_loc_r, n_loc_c, current_generation_loc, rank, size, c_generation);

        // fill intermediate result buffer for next communication step
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[0][0], top_rows_recv, 2 * n_loc_c * sizeof(uint8_t));
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[2 + n_loc_r][0], bottom_rows_recv, 2 * n_loc_c * sizeof(uint8_t));
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[2][0], current_generation_loc, n_loc_c * n_loc_r * sizeof(uint8_t));

        // print_matrix_par(n_loc_r + 4, n_loc_c, loc_matrix_and_neighbor_top_bottom_rows, rank, size, c_generation);

        // second communication round to get left right cols from neighboring processes
        int i = 0;
        for (int r = 0; i < (n_loc_r + 4); i++, r++) { // copy first and second to last column
            left_col[i] = loc_matrix_and_neighbor_top_bottom_rows[r][0];
            right_col[i] = loc_matrix_and_neighbor_top_bottom_rows[r][n_loc_c - 2];
        }
        for (int r = 0; i < 2 * (n_loc_r + 4); i++, r++) { // copy second and last column
            left_col[i] = loc_matrix_and_neighbor_top_bottom_rows[r][1];
            right_col[i] = loc_matrix_and_neighbor_top_bottom_rows[r][n_loc_c - 1];
        }

        MPI_Sendrecv(left_col, 2 * (n_loc_r + 4), MPI_UINT8_T, left_proc_neighbor, 0,
                     right_cols_recv, 2 * (n_loc_r + 4), MPI_UINT8_T, right_proc_neighbor, 0, cartcomm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(right_col, 2 * (n_loc_r + 4), MPI_UINT8_T, right_proc_neighbor, 0,
                     left_cols_recv, 2 * (n_loc_r + 4), MPI_UINT8_T, left_proc_neighbor, 0, cartcomm, MPI_STATUS_IGNORE);

        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        // copy left columns
        for (int r = 0; r < (n_loc_r + 4); r++) {
            full_current_generation_loc[r][0] = left_cols_recv[r];
            full_current_generation_loc[r][1] = left_cols_recv[(n_loc_r + 4) + r];
        }

        // copy right columns
        for (int r = 0; r < (n_loc_r + 4); r++) {
            full_current_generation_loc[r][(n_loc_c + 4 - 2)] = right_cols_recv[r];
            full_current_generation_loc[r][(n_loc_c + 4 - 1)] = right_cols_recv[(n_loc_r + 4) + r];
        }

        // fill remaining full local matrix
        for (int r = 0; r < (n_loc_r + 4); r++) {
            for (int c = 0; c < n_loc_c; c++) {
                full_current_generation_loc[r][2 + c] = loc_matrix_and_neighbor_top_bottom_rows[r][c];
            }
        }

        // print_matrix_par(n_loc_r + 4, n_loc_c + 4, full_current_generation_loc, rank, size, 0);
        run_generation_par(n_loc_r, n_loc_c, full_current_generation_loc, next_generation_loc, rank);
        copy_matrix_par(n_loc_r, n_loc_c, current_generation_loc, next_generation_loc, rank, size);

        print_matrix_par(n_loc_r, n_loc_c, current_generation_loc, rank, size, c_generation);

        if (verify) {
            // run iteration of sequential version
            run_generation(current_generation_seq, next_generation_seq, n);
            copy_matrix(current_generation_seq, next_generation_seq, n);

            // verification
            verification_result = compare_matrices(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);

            if (verification_result) {
                printf("Verification after generation %d on rank %d was successful\n", c_generation, rank);
            } else {
                printf("Verification after generation %d on rank %d failed\n", c_generation, rank);
            }

            if (verbose && rank == 0) {
                printf("Sequential matrix after generation %d: \n", c_generation);
                print_matrix(current_generation_seq, n);
            }
            // todo maybe: send all verification results to one rank and output overall result only once
        }
    }

    free(current_generation_loc);
    free(next_generation_loc);
    free(top_rows_recv);
    free(bottom_rows_recv);
    free(left_cols_recv);
    free(right_cols_recv);
    free(loc_matrix_and_neighbor_top_bottom_rows);
    free(full_current_generation_loc);

    if (rank == 0 && verify) {
        free(current_generation_seq);
        free(next_generation_seq);
    }

    MPI_Finalize();

    return 0;
}