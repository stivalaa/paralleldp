#!/bin/bash

# pbs launching script:
 	
#PBS -N parallel_hashtables

#PBS -l walltime=00:59:00
#PBS -q run_1_day

# Reserve 8 cores on a single node (need SMP shared memory)
#PBS -l nodes=1:ppn=8

#PBS -l pmem=2G

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alexs/phd/paralleldp/streamflow
export TBB21_INSTALL_DIR=/home/alexs/tbb21_20080605oss
. $TBB21_INSTALL_DIR/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/bin/tbbvars.sh

./timetests.sh httest oahttslftest 8 > oahttslftest.tango.rtab
./timetests.sh httest httslftest 8 > tango.rtab

./timetests.sh httest ../nbds.0.4.3/output/nbdstest 8 > nbdstest.tango.rtab
./timetests.sh httest tbbhashmaptest 8 > tbbhashmap.tango.rtab


