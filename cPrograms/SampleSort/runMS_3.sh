#!/bin/bash
mpiexec -n 5 -hostfile ../../hostFiles/host_file_3 ./SampleSortKernelized.bin
