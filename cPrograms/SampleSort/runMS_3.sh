#!/bin/bash
mpiexec -n 3 -hostfile ../../hostFiles/host_file_3 ./SampleSortKernelized.bin
