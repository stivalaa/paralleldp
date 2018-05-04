#!/bin/bash

# pbs launching script:
 	
#PBS -N knapsack_norandomization

#PBS -l walltime=48:00:0

# Reserve 8 cores on a single node (need SMP shared memory)
#PBS -l nodes=1:ppn=8


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alexs/phd/paralleldp/streamflow

time ./runall.sh '-z -y' tango_norandomization_output 8

