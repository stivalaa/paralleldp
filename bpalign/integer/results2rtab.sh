#!/bin/sh
###############################################################################
#
# results2rtab.sh - convert parbpalign output to table for R
#
# File:    results2rtab.sh
# Author:  Alex Stivala
# Created: June 2009
#
# Read the output from runall.bash and
# output data frame R  read.table(,header=TrUE)
# Usage:
#     results2rtab.sh <results_dir>
#
# $Id: results2rtab.sh 3041 2009-12-14 06:39:27Z alexs $
# 
###############################################################################

if [ $# -ne 1 ]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

resultsdir=$1

echo "threads iter ms re hc hn"
for resultsfile in ${resultsdir}/bpalign_*_*.out
do
    threads=`basename ${resultsfile} | cut -d'_'  -f2`
    iter=`basename ${resultsfile} | cut -d'_'  -f3 | cut -d'.' -f1`
    echo -n $threads " " $iter " "
    grep -v '^#' ${resultsfile} | awk '
      { ms += $5;
        re += $10;
        hc += $11;
        hn += $12;
      }
      END {
        printf("%d %d %d %d\n",ms / NR, re / NR, hc / NR, hn / NR);
      }
    '
done
