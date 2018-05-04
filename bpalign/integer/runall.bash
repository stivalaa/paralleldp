#!/bin/bash
#
# runall.bash - run the paralle bpalign program on all pairs of test data 
#
#
# Usage: runall.bash [bpalign options] outdir max_threads
#
# Creates outdir if it does not exist. WARNING: will overwrite files in
# outdir.
#
# Uses bash not sh in order to use arrays and arithmetic expansion
# Uses the runpairs.bash script to actually run the program on input data
#
# $Id: runall.bash 3173 2009-12-30 00:20:04Z alexs $
#


MAXITERATIONS=10

TIMEGUARD=../../utils/timeguard
TIMEOUT=7200

if [ $# -ne 3 ]; then
    echo "Usage: $0  parbpalign_options outdir max_threads" >&2
    exit 1
fi

bpopts=$1
outdir=$2
max_threads=$3

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi



numthreads=0
while [ $numthreads -le $max_threads ]; do
  iter=1
  while [ $iter -le $MAXITERATIONS ]; do
      outfilebase=${outdir}/bpalign_${numthreads}_${iter}
      echo "# Run as: " $0 $*  > ${outfilebase}.out
      echo "# at: " `date` >> ${outfilebase}.out
      echo "# on: " `uname -a` >> ${outfilebase}.out
      cat /dev/null > ${outfilebase}.err
      /usr/bin/time $TIMEGUARD $TIMEOUT runpairs.bash "${bpopts}" ${numthreads} >> ${outfilebase}.out 2> ${outfilebase}.err
      iter=`expr ${iter} + 1`
   done
  numthreads=`expr ${numthreads} + 1`
done


