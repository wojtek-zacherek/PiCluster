#!/bin/bash
MAXNUM=10000
MAXNUM=$1
mpiexec -n 5 -hostfile ../../hostFiles/host_file_5 ./MergeSort
