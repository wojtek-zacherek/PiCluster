#!/bin/bash
MAXNUM=10000
MAXNUM=$1

mpiexec -n 4 -hostfile ../../hostFiles/host_file_4 ./MatrixSubmatrix.bin
