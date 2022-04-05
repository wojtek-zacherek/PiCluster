#!/bin/bash
MAXNUM=10000
MAXNUM=$1

mpiexec -n 2 -hostfile ../../hostFiles/host_file_2 ./MatrixColumnVector.bin
