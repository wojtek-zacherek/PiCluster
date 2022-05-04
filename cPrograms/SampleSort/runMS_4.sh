#!/bin/bash
mpiexec -n 5 -hostfile ../../hostFiles/host_file_4 ./SampleSortKernelized.bin
