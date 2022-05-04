#!/bin/bash
mpiexec -n 5 -hostfile ../../hostFiles/host_file_2 ./SampleSortKernelized.bin
