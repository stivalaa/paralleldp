#!/bin/sh
###############################################################################
#
# runall.sh - Run run_instances.py with different numbers of cores
#
# File:    runall.sh
# Author:  Alex Stivala
# Created: April 2009
#
#
# Run the run_intances.py script with n threads, then n-1, etc. down to 1 thread
#
# Usage:
#  runall.sh  knapsack_options outdir max_threads 
#
#     knapsack_options are options to pass to knapsack_program
#     max_threads is maximum number of threads to use
#
# Creates outdir if it does not exist. WARNING: will overwrite files in
# outdir.
#
# Requires run_instances.py to be in PATH
#
# $Id: runall.sh 3190 2010-01-02 23:58:09Z alexs $
# 
###############################################################################

MAXITERATIONS=10

TIMEGUARD=../utils/timeguard
TIMEOUT=25200

if [ $# -ne 3 ]; then
    echo "Usage: $0  knapsack_options outdir max_threads" >&2
    exit 1
fi

knapsack_options=$1
outdir=$2
max_threads=$3

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi


numthreads=$max_threads
while [ $numthreads -ge 1 ]; do
  iter=1
  while [ $iter -le $MAXITERATIONS ]; do
      /usr/bin/time $TIMEGUARD $TIMEOUT run_instances.py  knapsack_oahttslf "-r${numthreads} ${knapsack_options}" instances > ${outdir}/instances_oahttslf_${numthreads}_${iter}.out 2> ${outdir}/instances_oahttslf_${numthreads}_${iter}.err
      iter=`expr $iter + 1`
  done
  numthreads=`expr $numthreads - 1`
done

/usr/bin/time $TIMEGUARD $TIMEOUT run_instances.py  knapsack_simple "${knapsack_options}" instances > ${outdir}/instances_simple.out 2> ${outdir}/instances_simple.err

