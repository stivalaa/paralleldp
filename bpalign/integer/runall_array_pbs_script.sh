#!/bin/bash

# pbs launching script:
 	
#PBS -N parbpalign

#PBS -l walltime=40:00:0

# Reserve 8 cores on a single node (need SMP shared memory)
#PBS -l nodes=1:ppn=8


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alexs/phd/paralleldp/streamflow

time ./runall.bash '-b -s' tango_bottomup_output 0
time ./runall.bash '-s -a' tango_array_output 8

