
# Modified Game of Life MPI Implementation

## Overview

This project implements a modified version of Conway's Game of Life using MPI for parallel processing. It includes three different implementations:
- Sequential version (`seq.c`)
- Parallel version with point-to-point communication (`par.c`)
- Parallel version with collective communication (`par_collective.c`)

This project was done as part of the course "High Performance Computing" at TU Vienna. An exhaustive report on the project including detailed performance results and more implementation details is available [here](Report.pdf).

## Game Rules


The implementation uses a modified version of Conway's Game of Life with the following rules:

1. The grid is toroidal (wraps around at edges)
2. Each cell has 8 neighbors in a specific pattern
3. The state of a cell in the next generation is determined by the state of the previous generation according to the following rules:
   - A cell will survive, if it is alive and has a two or three (2 or 3) alive cells in its neighborhood. 
   - A cell will be reborn, if it is dead and has exactly three (3) alive cells in its neighborhood.
   - A cell that is not surviving or reborn, will die or remain dead.

## Building the Project

### Prerequisites
- MPI implementation (e.g., OpenMPI, MPICH)
- C Compiler with C99 support
- Make

### Compilation


Use the provided Makefile to build the different versions:

```bash
# Build sequential version
make seq

# Build parallel version
make par

# Build parallel collective version
make par_collective

# Clean all builds
make clean
```


## Usage

### Command Line Arguments

All implementations support the following arguments:
- `-n, --number`: Grid size (N x N)
- `-g, --num_generations`: Number of generations to simulate
- `-s, --seed`: Random seed for reproducibility
- `-d, --density`: Initial cell density (percentage, 1-100)
- `-v, --verbose`: Enable verbose output
- `-x, --verify`: Enable verification against sequential version
- `-w, --weak_scaling`: Enable weak scaling mode


### Running the Programs
```bash
# Sequential version
mpirun -np 1 ./seq -n 100 -g 1000 -s 42 -d 28

# Parallel versions (with 4 processes)
mpirun -np 4 ./par -n 100 -g 1000 -s 42 -d 28
mpirun -np 4 ./par_collective -n 100 -g 1000 -s 42 -d 28
```


## Implementation Details


### Sequential Version

- Single-process implementation
- Uses modulo operations for boundary conditions
- Serves as reference for correctness verification

### Parallel Version (Point-to-Point)
- Uses 2D domain decomposition
- Creates a Cartesian process grid
- Implements neighbor communication using MPI_Sendrecv
- Neighbor communication is done in two steps:
  - First step: Communicates the two neighboring top and bottom rows with neighboring processes
  - Second step: Communicates the two neighboring left and right columns with neighboring processes

### Parallel Collective Version
- Similar to parallel version but uses collective communication
- Implements neighbor communication using MPI_Neighbor_alltoallv
- Uses distributed graph topology for optimized communication

## Verification

When running with the -x flag, the parallel versions verify their results against the sequential implementation:
- Compares initial state
- Verifies final state after all generations

## Limitations

- Grid size (N) must be divisible by the number of processes in both dimensions

## Output

The program provides:
- Execution time in microseconds
- Number of alive/dead cells (in verbose mode)
- Grid state visualization (in verbose mode)
- Verification results (when enabled)
