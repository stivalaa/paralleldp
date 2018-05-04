#!/bin/sh
#
# $Id: summarize_httslf.sh 3152 2009-12-28 00:07:23Z alexs $
#
# sum up the mean elapsed times (ms) in each STATS line
#

grep '# on:' instances_simple.out
echo "instances_simple.out: " `grep STATS instances_simple.out | awk '{sum += $2} END {print sum}'`
for i in instances_httslf_*_*.out
do
    echo "${i}: " `grep STATS $i | awk '{sum += $2} END {print sum}'`
done
