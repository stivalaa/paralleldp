#!/bin/sh
#
# buildinputdata.sh - build input bplist files for parbpalign
#
# The buildpmatrices.py script creates data in ./data with RNAfold
# and data2list.sh converts to basepair probability litss in ./datalist
# 
# $Id: buildinputdata.sh 2612 2009-07-06 00:54:51Z astivala $
#

PYTHONPATH=${PYTHONPATH}:../../bpalign
PATH=${PATH}:../../bpalign

if [ ! -d data ]; then
    mkdir data
fi
cd data
../buildbpmatrices.py
cd ..
data2list.sh 
cp data/pairs.list datalist/pairs.list
