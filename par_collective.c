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
    num_alive_cells += matrix[i][j-2];
    num_alive_cells += matrix[i+2][j-2];
    num_alive_cells += matrix[i-1][j];
    num_alive_cells += matrix[i+2][j];
    num_alive_cells += matrix[i-1][j+1];
    num_alive_cells += matrix[i][j+1];
    num_alive_cells += matrix[i+2][j+2];

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
    int n = 10, n_generations = 1000; // num rows, num cols, num generations
    int seed = 42;
    int verbose = 0;
    int verify = 0; // if set we perform verification with the sequential version
    int weak_scaling = 0; // if set input size is per processor
    int density = 28; // in percent
    int opt;
    int rank, size;
    double total_time = 0.0, total_time_process, start_time, end_time;

    int n_loc_r, n_loc_c;
    int nprows, npcols;
    int prow_idx, pcol_idx;
    int verification_result;
    int process_verification_result = 0;

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

    MPI_Cart_coords(cartcomm_reorder, rank, 2, coords); // get process coordinates in created grid

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

        copy_full_matrix_to_local_matrix(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);
    } else {
        // fill local matrix with initial input
        fill_matrix_par(n_loc_r, n_loc_c, current_generation_loc, n, density, m_offset_r, m_offset_c);
    }

    if (verbose) {
        print_matrix_par(n_loc_r, n_loc_c, current_generation_loc, rank, size, 0);
    }

    // verify initial input
    if (verify) {
        process_verification_result = compare_matrices(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);
        MPI_Reduce(&process_verification_result, &verification_result, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            if (verification_result) {
                printf("Input verification was successful\n");
            } else {
                printf("Input verification failed\n");
            }
        }

        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int top_bottom_neighbors[2];
    int left_right_neighbors[2];
    int weights[2] = {1, 1};

    // Determine neighbors
    MPI_Cart_shift(cartcomm, 1, 1, &left_right_neighbors[0], &left_right_neighbors[1]);
    MPI_Cart_shift(cartcomm, 0, 1, &top_bottom_neighbors[0], &top_bottom_neighbors[1]);

    MPI_Comm dist_graph_top_bottom;
    MPI_Comm dist_graph_left_right;

    // create top bottom neighbor graph
    MPI_Dist_graph_create_adjacent(
            MPI_COMM_WORLD,
            2, // number of sources
            top_bottom_neighbors, // sources
            weights,
            2, // number of destinations
            top_bottom_neighbors, // destinations
            weights,
            MPI_INFO_NULL,
            0, // no reordering
            &dist_graph_top_bottom
    );
    // create left right neighbor graph
    MPI_Dist_graph_create_adjacent(
            MPI_COMM_WORLD,
            2, // number of sources
            left_right_neighbors, // sources
            weights,
            2, // number of destinations
            left_right_neighbors, // destinations
            weights,
            MPI_INFO_NULL,
            0, // no reordering
            &dist_graph_left_right
    );

    // communication setup
    int sendcounts_top_bottom[2] = {4 * n_loc_c, 4 * n_loc_c};
    int recvcounts_top_bottom[2] = {4 * n_loc_c, 4 * n_loc_c};
    int sendcounts_left_right[2] = {4 * (n_loc_r + 4), 4 * (n_loc_r + 4)};
    int recvcounts_left_right[2] = {4 * (n_loc_r + 4), 4 * (n_loc_r + 4)};

    int sdispls[2] = {0, 0};
    int rdispls[2] = {0, 0};

    // communication buffers
    uint8_t *sendbuf_top_bottom = malloc(4 * n_loc_c * sizeof(uint8_t));
    uint8_t *recvbuf_top_bottom = malloc(4 * n_loc_c * sizeof(uint8_t));
    uint8_t *sendbuf_left_right = malloc(4 * (n_loc_r + 4) * sizeof(uint8_t));
    uint8_t *recvbuf_left_right = malloc(4 * (n_loc_r + 4) * sizeof(uint8_t));

    // intermediate result buffer after top and bottom rows are communicated
    uint8_t(*loc_matrix_and_neighbor_top_bottom_rows)[n_loc_c];
    loc_matrix_and_neighbor_top_bottom_rows = (uint8_t(*)[n_loc_c])malloc((n_loc_r + 4) * n_loc_c * sizeof(uint8_t));

    // final intermediate result buffer with all neighbors
    uint8_t(*full_current_generation_loc)[n_loc_c + 4];
    full_current_generation_loc = (uint8_t(*)[n_loc_c + 4])malloc((n_loc_r + 4) * (n_loc_c + 4) * sizeof(uint8_t));

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // run gol
    for (int c_generation = 1; c_generation <= n_generations; c_generation++) {
        // fill buffers for communication with top and bottom neighbors
        memcpy(&sendbuf_top_bottom[0], current_generation_loc[0], 2 * n_loc_c * sizeof(uint8_t));
        memcpy(&sendbuf_top_bottom[2 * n_loc_c], current_generation_loc[n_loc_r - 2], 2 * n_loc_c * sizeof(uint8_t));

        // first communication round with top and bottom neighbors
        MPI_Neighbor_alltoallv(
                sendbuf_top_bottom, sendcounts_top_bottom, sdispls, MPI_UINT8_T,
                recvbuf_top_bottom, recvcounts_top_bottom, rdispls, MPI_UINT8_T,
                dist_graph_top_bottom
        );

        // fill intermediate result buffer for next communication step
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[0][0], &recvbuf_top_bottom[2 * n_loc_c], 2 * n_loc_c * sizeof(uint8_t));
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[2 + n_loc_r][0], recvbuf_top_bottom, 2 * n_loc_c * sizeof(uint8_t));
        memcpy(&loc_matrix_and_neighbor_top_bottom_rows[2][0], current_generation_loc, n_loc_c * n_loc_r * sizeof(uint8_t));

        // fill buffers for communication with left and right neighbors
        for (int r = 0; r < (n_loc_r + 4); r++) {
            sendbuf_left_right[r] = loc_matrix_and_neighbor_top_bottom_rows[r][0]; // copy first column
            sendbuf_left_right[r+(n_loc_r + 4)] = loc_matrix_and_neighbor_top_bottom_rows[r][1]; // copy second column
            sendbuf_left_right[r+2*(n_loc_r + 4)] = loc_matrix_and_neighbor_top_bottom_rows[r][n_loc_c - 2]; // copy second to last column
            sendbuf_left_right[r+3*(n_loc_r + 4)] = loc_matrix_and_neighbor_top_bottom_rows[r][n_loc_c - 1]; // copy last column
        }

        // second communication round with left and right neighbors
        MPI_Neighbor_alltoallv(
                sendbuf_left_right, sendcounts_left_right, sdispls, MPI_UINT8_T,
                recvbuf_left_right, recvcounts_left_right, rdispls, MPI_UINT8_T,
                dist_graph_left_right
        );

        // fill final buffer
        for (int r = 0; r < (n_loc_r + 4); r++) { // fill left and right cols
            full_current_generation_loc[r][0] = recvbuf_left_right[r+2*(n_loc_r+4)];; // copy first column
            full_current_generation_loc[r][1] = recvbuf_left_right[r+3*(n_loc_r+4)]; // copy second column
            full_current_generation_loc[r][(n_loc_c + 4 - 2)] = recvbuf_left_right[r]; // copy second to last column
            full_current_generation_loc[r][(n_loc_c + 4 - 1)] = recvbuf_left_right[r+(n_loc_r+4)]; // copy last column
        }
        for (int r = 0; r < (n_loc_r + 4); r++) { // fill remaining full local matrix
            for (int c = 0; c < n_loc_c; c++) {
                full_current_generation_loc[r][2 + c] = loc_matrix_and_neighbor_top_bottom_rows[r][c];
            }
        }

        run_generation_par(n_loc_r, n_loc_c, full_current_generation_loc, next_generation_loc, rank);
        copy_matrix_par(n_loc_r, n_loc_c, current_generation_loc, next_generation_loc, rank, size);

        if (verbose) {
            print_matrix_par(n_loc_r, n_loc_c, current_generation_loc, rank, size, c_generation);
        }

        if (verify) {
            // run iteration of sequential version
            run_generation(current_generation_seq, next_generation_seq, n);
            copy_matrix(current_generation_seq, next_generation_seq, n);

            if (c_generation == n_generations) {
                // gather verification results
                process_verification_result = compare_matrices(n_loc_r, n_loc_c, current_generation_seq, current_generation_loc, n, m_offset_r, m_offset_c);
                MPI_Reduce(&process_verification_result, &verification_result, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

                if (rank == 0) {
                    if (verification_result) {
                        printf("Verification after generation %d was successful\n", c_generation);
                    } else {
                        printf("Verification after generation %d failed\n", c_generation);
                    }
                }

                if (verbose && rank == 0) {
                    printf("Sequential matrix after generation %d: \n", c_generation);
                    print_matrix(current_generation_seq, n);
                }
            }
        }
    }

    end_time = MPI_Wtime();
    total_time_process = (end_time - start_time) * 1e6; // μs
    MPI_Reduce(&total_time_process, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Total time: %.2f μs\n", total_time);
    }

    free(current_generation_loc);
    free(next_generation_loc);
    free(sendbuf_top_bottom);
    free(recvbuf_top_bottom);
    free(sendbuf_left_right);
    free(recvbuf_left_right);
    free(loc_matrix_and_neighbor_top_bottom_rows);
    free(full_current_generation_loc);

    if (rank == 0 && verify) {
        free(current_generation_seq);
        free(next_generation_seq);
    }

    MPI_Finalize();

    return 0;
}
