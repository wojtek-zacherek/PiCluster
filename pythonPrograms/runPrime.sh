#!/bin/bash
k=300000
mpiexec -n 6 -hostfile ../hostFiles/host_file python3 prime.py $k
mpiexec -n 5 -hostfile ../hostFiles/host_file_5 python3 prime.py $k
mpiexec -n 4 -hostfile ../hostFiles/host_file_4 python3 prime.py $k
mpiexec -n 3 -hostfile ../hostFiles/host_file_3 python3 prime.py $k
mpiexec -n 2 -hostfile ../hostFiles/host_file_2 python3 prime.py $k
mpiexec -n 1 -hostfile ../hostFiles/host_file_1 python3 prime.py $k
