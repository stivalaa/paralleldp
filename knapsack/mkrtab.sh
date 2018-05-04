#!/bin/sh
###############################################################################
#
# mkrtab.sh - make table for R from summarize.sh output
#
# File:    mkrtab.sh
# Author:  Alex Stivala
# Created: May 2009
#
#
# Make table for R read.table(,header=TRUE) from the output of summarize.sh,
# summaries of the average times for knapsack instances.
#
# Usage:
#  summarize.sh | mkrtab.sh 
#
#    output is to stdout.
#    input on stdin is output from summarize.sh
#
# $Id: mkrtab.sh 3115 2009-12-23 21:28:07Z alexs $
# 
###############################################################################

echo "threads iter ms"

grep -v '^#' | while read line
do
    if [ `expr "${line}" : ".*_simple.*"` -ne 0 ]; then
        iter=1
        thread=0
    else
        name=`echo "$line"|cut -d' ' -f1`
        thread=`echo $name | cut -d '_' -f3`
        iter=`echo $name | cut -d '_' -f4 | cut -d '.' -f1`
    fi
    ms=`echo "$line"|awk '{print $2}'`
    echo $thread $iter $ms
done

