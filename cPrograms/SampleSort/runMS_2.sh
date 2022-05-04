#!/bin/bash
mpiexec -n 2 -hostfile ../../hostFiles/host_file_2 ./SampleSortKernelized.bin
