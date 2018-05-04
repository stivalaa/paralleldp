#!/bin/sh
#
# data2list.sh - convert all RNAfold -p dotplot postscript files in the
#                data directory to .bplist files in the datalist directory
#                using rnafol2list.py, for use by C bpalign program.
#
# $Id: data2list.sh 251 2007-06-12 05:56:36Z astivala $
#

DATADIR=data   # input directory with _dp.ps files
DATALISTDIR=datalist # output directory for .bplist files

if [ ! -d ${DATALISTDIR} ]; then
    mkdir ${DATALISTDIR}
fi

for i in ${DATADIR}/*_dp.ps
do
  rnafold2list.py $i > ${DATALISTDIR}/`basename $i _dp.ps`.bplist
done
