#!/bin/bash
mpiexec -n 4 -hostfile ../../hostFiles/host_file_4 ./SampleSortKernelized.bin
