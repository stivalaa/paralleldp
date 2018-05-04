#!/bin/sh
###############################################################################
#
# timetests.sh - Run test harness programs and make table of times
#
# File:    timetests.sh
# Author:  Alex Stivala
# Created: April 2009
#
#
# Run the test harness wit 1 thread, then 2, etc. up to n threads
# and make table for R read.table(,header=TRUE)
#
# Usage:
#  timetests.sh  baseline program max_threads 
#
#     baseline is the baseline program to run (eg httest), this is "0 threads"
#     program is the test harness program to run  (eg httslftest)
#              (takes num_threads as parameter, and writes a 
#               "elapsed time <x> ms" string to stdout where <x> is time in ms
#     max_threads is maximum number of threads to use
#
#   output to stdout is table for reading in R
#
# $Id: timetests.sh 2316 2009-05-05 08:38:59Z astivala $
# 
###############################################################################

# number of times to repeat each run
NUM_REPETITIONS=10


if [ $# -ne 3 ]; then
    echo "Usage: $0 baseline  program  max_threads" 2>&1
    exit 1
fi

baseline=$1
program=$2
max_threads=$3


echo "# Run as: $0 $*"
echo "# at: `date`"
echo "# by: `whoami`"
echo "# on: `uname -a`"

echo "threads iteration ms"

i=1
while [ $i -le $NUM_REPETITIONS ]; do
    ms=`${baseline} | grep "^elapsed time" | cut -d" " -f3`
    echo 0 ${i} ${ms}
    i=`expr $i + 1`
done

numthreads=1
while [ $numthreads -le $max_threads ]; do
    i=1
    while [ $i -le $NUM_REPETITIONS ]; do
        ms=`${program}  ${numthreads} | grep "^elapsed time" | cut -d" " -f3`
        echo ${numthreads} ${i} ${ms}
        i=`expr $i + 1`
    done
    numthreads=`expr $numthreads + 1`
done

