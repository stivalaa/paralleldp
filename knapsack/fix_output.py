#!/usr/bin/env python
###############################################################################
#
# fix_output.py - Reparse run_instances.py output and regenerate STATS lines
#
# File:    fix_output.py
# Author:  Alex Stivala
# Created: June 2009
#
# Parse run_instances.py output result lines to get time, and compute stats
# (min, max,median,mean).
#
# Output looks like
# 119612 30723474 30961478 68093 17485  gen.2.500.500.99
#
# where the numberc values are in order:
#
# profit total_resuse total_hashcount cputime elapsedtime
#
# Then followed by flags and then problem description.
# We copy this to stdout, and for each problem file write the stats over
# the instances in that file as
#
# STATS: mean_elapsed median_elapsed max_elapsed flags problem-name
#
# times are in ms
#
# If one of the instances exceeds the limit, the stats are not computed
# and we just print
# STATS: LIMIT EXCEEDED description
#
# Only needed because of run of knapsack instrumentation
# with bug in run_instances.py;
# this script used to regnerazte the STATS lines correctly from the
# knapsack output.
#
# Usage:
#     fix_output.py 
#
#     The run_instances output is read from stdin, output is to stdout.
#
#
# $Id: fix_output.py 2482 2009-06-05 01:19:21Z astivala $
# 
###############################################################################

import sys,os,glob

#from numpy import mean,median
# some systems e.g. mungera don't have numpy so can't use it, have to
# implement mean,median here
def mean(l):
   """
   compute arithmetic mean of a list of numbers
   """
   return float(sum(l)) / float(len(l))

def median(l):
   """ compute median of a list of numbers
   """
   s = sorted(l)
   if len(s) % 2 == 0:
       return (float(s[len(s)/2 - 1] + s[len(s)/2]))/2
   else:
       return s[len(s)/2]


def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + "\n")
    sys.exit(1)


    
def main():
    """
    main for fix_output.py
    """
    if len(sys.argv) != 1:
        usage(os.path.basename(sys.argv[0]))

    ktype = last_ktype = None
    elapsed_list = []
    exceeded = False
    for results in sys.stdin:
        if results[0] == '#':
            sys.stdout.write(results)
            continue
        if results[:5] == 'STATS':
            continue # discard old STATS lines
        if results[:10] == 'INSTRUMENT':
            sys.stdout.write(results)
            continue
        sresults = results.split()
        if sresults[0] == "LIMIT":
            description = ' '.join(sresults[3:])
            exceeded = True
            break
        sys.stdout.write(results)
        ttime = int(sresults[3])
        etime = int(sresults[4])
        flags = sresults[5]
        name = ' '.join(sresults[6:])
        last_ktype = ktype
        ktype  = int(name.split('.')[1])
        elapsed_list.append(etime)

        if last_ktype != None and  ktype != last_ktype:
            new_etime = elapsed_list.pop()
            # STATS: mean_elapsed median_elapsed max_elapsed flags poblem-name
            if exceeded:
                sys.stdout.write("STATS: LIMIT EXCEEDED " + description + '\n')
            else:
                sys.stdout.write('STATS: %f %f %d %s %s\n' %
                    (mean(elapsed_list), median(elapsed_list), max(elapsed_list),
                     flags, str(last_ktype))
                            )
            elapsed_list = [new_etime]
            exceeded = False



    # STATS: mean_elapsed median_elapsed max_elapsed flags poblem-name
    if exceeded:
        sys.stdout.write("STATS: LIMIT EXCEEDED " + description + '\n')
    else:
        sys.stdout.write('STATS: %f %f %d %s %s\n' %
            (mean(elapsed_list), median(elapsed_list), max(elapsed_list),
             flags, str(last_ktype))
                    )

if __name__ == "__main__":
    main()

