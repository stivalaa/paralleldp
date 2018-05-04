#!/bin/bash

# pbs launching script:
 	
#PBS -N knapsack_httslf

#PBS -l walltime=120:00:0

# Reserve 8 cores on a single node (need SMP shared memory)
#PBS -l nodes=1:ppn=8

#PBS -l pmem=2G

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alexs/phd/paralleldp/streamflow

time ./runall_httslf.sh '-y' tango_httslf_output 8

