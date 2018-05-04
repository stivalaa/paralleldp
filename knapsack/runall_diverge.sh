#!/bin/sh
###############################################################################
#
# runall_diverge.sh - Run run_instances.py with different numbers of cores
#
# File:    runall_diverge.sh
# Author:  Alex Stivala
# Created: April 2009
#
#
# Run the run_intances.py script with 1 thread, then 2, 4, 8, up to n threads
#
# Usage:
#  runall_diverge.sh  knapsack_options outdir max_threads 
#
#     knapsack_options are options to pass to knapsack_program
#     max_threads is maximum number of threads to use
#
# Creates outdir if it does not exist. WARNING: will overwrite files in
# outdir.
#
# Requires run_instances.py to be in PATH
#
# $Id: runall_diverge.sh 3206 2010-01-05 00:04:49Z alexs $
# 
###############################################################################

#MAXITERATIONS=10
MAXITERATIONS=1

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

#/usr/bin/time $TIMEGUARD $TIMEOUT run_instances.py  knapsack_simple "${knapsack_options}" instances > ${outdir}/instances_simple.out 2> ${outdir}/instances_simple.err

numthreads=1
while [ $numthreads -le $max_threads ]; do
  iter=1
  while [ $iter -le $MAXITERATIONS ]; do
      /usr/bin/time $TIMEGUARD $TIMEOUT run_instances.py  knapsack_diverge_oahttslf "-r${numthreads} ${knapsack_options}" instances > ${outdir}/instances_oahttslf_${numthreads}_${iter}.out 2> ${outdir}/instances_oahttslf_${numthreads}_${iter}.err
      iter=`expr $iter + 1`
  done
  numthreads=`expr $numthreads + 1`
done

