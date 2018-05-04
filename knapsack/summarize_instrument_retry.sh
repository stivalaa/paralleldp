#!/bin/sh
###############################################################################
#
# summarize_instrument_retry.sh - make table for R of instrumentation counts
#
# File:    summarize_instrument_retry.sh
# Author:  Alex Stivala
# Created: May 2009
#
#
# sum up the reuse and hash counts and retry counts
# in each INSTRUMENT line, and output
# in format for R read.table(,header=TRUE)
#
# Usage:
#  summarize_instrument_retry.sh
#
#    output is to stdout.
#    run from the output directory (reads files in cwd)
#
# Requires GNU gawk
#
# $Id: summarize_instrument_retry.sh 3162 2009-12-29 02:37:16Z alexs $
# 
###############################################################################

grep '# on:' instances_oahttslf_1_1.out
echo "threads re hc hn or"
for i in instances_oahttslf_*.out
do
  thread=`echo ${i} | cut -d '_' -f3| cut -d '.' -f1`
  echo -n "${thread} "
  grep INSTRUMENT ${i} |  gawk '{
    hc = gensub(/.* hc=([0-9]*).*/,"\\1",1); 
    re = gensub(/.*,re=([0-9]*).*/,"\\1",1);
    hn = gensub(/.*,hn=([0-9]*).*/,"\\1",1);
    orx = gensub(/.*,or=([0-9]*).*/,"\\1",1);
    total_hc += hc; 
    total_re += re;
    total_hn += hn;
    total_orx += orx;
    } 
    END {print total_re,total_hc,total_hn,total_orx}'
done



