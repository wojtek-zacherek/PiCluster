#!/bin/bash
MAXNUM=10000
MAXNUM=$1

mpiexec -n 6 -hostfile ../../hostFiles/host_file ./MatrixSubmatrix.bin
