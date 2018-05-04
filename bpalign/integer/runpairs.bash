#!/bin/bash
#
# runpairs.bash - run the paralle bpalign program on all pairs of test data 
#
#
# Usage: runpairs.bash bpalign options  numthreads
#
# Uses bash not sh in order to use arrays and arithmetic expansion
#
# $Id: runpairs.bash 3054 2009-12-16 01:01:58Z alexs $
#

# location of input secondary structure .bplist files from rnafold2list.py
#DATADIR=${CLUSTER_HOME}/phd/bpalign/datalist
DATADIR=../datalist
PAIRLIST=${DATADIR}/pairs.list

if [ $# -ne 2 ]; then
    echo "Usage: $0  parbpalign_options numthreads" >&2
    exit 1
fi

bpopts=$1
numthreads=$2

cat ${PAIRLIST} | while read line
do
  seqid1=`echo "${line}" | cut -d' ' -f1`
  seqid2=`echo "${line}" | cut -d' ' -f2`
  infile1=${DATADIR}/${seqid1}.bplist
  infile2=${DATADIR}/${seqid2}.bplist
  echo -n ${seqid1} ${seqid2} ' '
  if [ ${numthreads} -eq 0 ]; then
      ./parbpalign_rand ${bpopts} ${infile1} ${infile2}
  else
      ./parbpalign_rand -t${numthreads} ${bpopts} ${infile1} ${infile2}
  fi
done

