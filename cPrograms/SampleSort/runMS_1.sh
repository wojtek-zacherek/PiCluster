#!/bin/bash
mpiexec -n 5 -hostfile ../../hostFiles/host_file_1 ./SampleSortKernelized.bin
