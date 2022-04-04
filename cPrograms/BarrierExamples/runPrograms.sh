#!/bin/bash
MAXNUM=10000
MAXNUM=$1
echo '--hello world barrier--'
mpiexec -n 6 -hostfile ../../hostFiles/host_file ./hwBarrOrder.bin

echo '--msg pass barrier--'
mpiexec -n 6 -hostfile ../../hostFiles/host_file ./msgPassBarr.bin

echo '--scatter gather--'
mpiexec -n 6 -hostfile ../../hostFiles/host_file ./scatterGather.bin
