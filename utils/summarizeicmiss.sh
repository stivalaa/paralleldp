#!/bin/sh
###############################################################################
#
# summarizeicmiss.sh - make table for R from cputrack output (IC_miss counter)
#
# File:    summarizeicmiss.sh
# Author:  Alex Stivala
# Created: December 2009
#
#
# Make table for R read.table(,header=TRUE) from the output of Solaris 
# cputrack for UltraSPARC T1, as run by collecticmiss.sh,
# with IC_miss and Instr_cnt counters.
#
# Usage:
#  summarizeicmiss.sh outputdir
#
#    output is to stdout.
#    input is assumed to be in the outputdir, collected by collecticmiss.sh
#
# $Id: summarizeicmiss.sh 3184 2010-01-01 04:10:13Z alexs $
# 
###############################################################################

if [ $# -ne 1 ]; then
    echo usage: $0 outputdir >&2
    exit 1
fi
outputdir=$1
echo threads iter ms lwp icmiss instrcnt
( cd $outputdir >/dev/null ;
n=1;
while [ $n -le 32 ]; do 
    i=1
    while  [ $i -le 10 ]; do 
        etime=`grep 'elapsed time' icmisstest_${n}_${i}.out | awk '{print $3}'` 
        grep 'lwp_exit' icmisstest_${n}_${i}.out | while read line ; do
            lwp=`echo $line | awk '{print $2}'`
            icmiss=`echo $line | awk '{print $4}'`
            instrcnt=`echo $line | awk '{print $5}'`
            echo $n $i $etime $lwp $icmiss $instrcnt
        done
        i=`expr $i + 1`
        done 
    n=`expr $n + 1`
done
)
