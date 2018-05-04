#!/bin/sh
#
# $Id: summarize.sh 3043 2009-12-14 23:50:40Z alexs $
#
# sum up the mean elapsed times (ms) in each STATS line
#

grep '# on:' instances_simple.out
echo "instances_simple.out: " `grep STATS instances_simple.out | awk '{sum += $2} END {print sum}'`
for i in instances_oahttslf_*_*.out
do
    echo "${i}: " `grep STATS $i | awk '{sum += $2} END {print sum}'`
done
