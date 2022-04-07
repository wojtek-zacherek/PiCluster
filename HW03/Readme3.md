# HW-03
## Acquiring files.
If you currently don't have access to these source files (i.e. you're not Dr. Leader), you may clone the parent repositroy at this link: https://github.com/wojtek-zacherek/PiCluster.git. The relevant files are under the `HW03/` directory below the parent level.

To compile, simply run the `make` command. The `Makefile` simply compiles all `*.c` files in the working directory to executable files with the `.bin` extension.

To execute, first allocate a number of processors on the computing node:
```
qsub -I -l nodes=<# of nodes>:ppn=<# of processes per node>
```
Once allocated, simply run the executables with command line arguments:
```
mpirun program.bin #arg1 #arg2 ...
```

## Matrix Column-Vector Multiplication
This programs performs a matrix-vector operation by dividing the multiplcation column-wise: each process acquires columns from the matrix and performs a partial summation for each row. The program utilizes MPI_Reduce to quickly sum all the partial sums from the various processes.

This implementation assumes that the number of columns divided per process is distributed evenly. It also generates it's own random matrix and vector of values.

The `alpha` parameter controls how many columns each process will receive, and directly affects `N`, the dimensions of the matrix and vector. You may set this value manually in the source code, recompile, and run to witness changes.

## Matrix Submatrix Multiplication
This program subdivides a matrix into smaller "submatrix" components, with each submatrix assigned to a processor for processing. Each processor performs a partial sum for the given submatrix and vector. This prorgam utilizes MPI_Reduce to quickly sum all the partial-partial sums.

This implementation assumes that the number of processors, `comm_sz`, is a perfect square. Additionally, it also assumes that each submatrix is a square matrix as well. 

The `alpha` parameter controls the dimensioanilty of the overall matrix and submatricies. To increase the size, simply increase the `alpha` parameter.

Note: please run this program with a square number of processors: 1, 4, 9,...

## MergeSort
This program generates a random array of values across multiple processors, sorts them locally with the `qsort` algorithm, and merges upwards towards process 0.

This implemnetation assumes that each initial array has the same number of elements. However, this should work for any number of processes.

The `alpha` parameter controls how many values are generated for each array. To modify the size, modify the `alpha` value in the source, recompile, and re-run.


