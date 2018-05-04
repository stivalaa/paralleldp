#!/bin/sh
###############################################################################
#
# generate_problems.sh - generate knapsack problem instances with gen2
#
# File:    generate_problems.sh
# Author:  Alex Stivala
# Created: May 2009
#
#
# Run the gen2 program from http://www.diku.dk/~pisinger/gen2.c
# to generate set of knapsack problem instances.
#
# Usage:
#   geneate_problems.sh outdir
#
#     outdir is the name of directory to place problem instances in
#
# Creates outdir if it does not exist. WARNING: will overwrite files in
# outdir. Also overwites test.in in cwd (gen2 always writes to test.in)
#
# Requires gen2 to be in PATH
#
# $Id: generate_problems.sh 2432 2009-05-22 06:11:51Z astivala $
# 
###############################################################################

if [ $# -ne 1 ]; then
    echo "Usage: $0 outdir" >&2
    exit 1
fi

outdir=$1

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

NUM_ITEMS=500
COEF_RANGE=500

NUM_INSTANCES=100  # number of instances for each type

for ktype in 1 2 3 4 5
do
  i=1
  while [ $i -le $NUM_INSTANCES ]; do
      gen2 $NUM_ITEMS $COEF_RANGE $ktype $i $NUM_INSTANCES
      name="gen.${ktype}.${NUM_ITEMS}.${COEF_RANGE}.${i}"
      echo $name > ${outdir}/$name
      cat test.in >> ${outdir}/$name
      i=`expr $i + 1`
  done
done
