#!/usr/bin/env python
###############################################################################
#
# run_instances.py - Run knapsack on the generated test instances
#
# File:    run_instances.py
# Author:  Alex Stivala
# Created: April 2009
#
# Run a knapsack implementation on all generated test instances
# (from generate_problems.sh)
#  http://www.dcs.st-and.ac.uk/~ipg/challenge
#
# Also parses the result line to get time, and compute stats
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
# Each instance is run 10 times and the averages for each problem set are
# are over all the repeitions of all the instances.
#
# Usage:
#     run_instances.py knapsack_program knapsack_options instances_dir
#
#     knapsack_program is the solver to run
#     knapsack_options are options to pass to knapsack_program
#     instances_dir is the directory with probloem instance files
#
#     The problem file is read from stdin.
#
#
# $Id: run_instances.py 2481 2009-06-04 23:50:05Z astivala $
# 
###############################################################################

import sys,os,glob
from time import clock,strftime,localtime
import getpass

# number of times to run each instance
NUM_REPETITIONS = 1

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
    sys.stderr.write("Usage: " + progname + " <knapsack-program> <knapsack-options> <instances-dir>\n")
    sys.exit(1)


def instance_generator(instances_dir, ktype):
    """
    Generator function to yield each problem instance in the directory
    file. Each instance is delimited by a blank line.

    Filenames are in the format gen.<ktype>.500.500.<i>

    where <ktype> is the problem type and <i> is the instance number
    
    Paramaters:
        instances_dir -name of directory continaing instances files
        ktype   - type of problem (1,2,3,4,5,...)

    Return value:
        YIELDs a string which is a problem instance
    """
    for instance_file in glob.glob(os.path.join(instances_dir,
                                              'gen.'+ str(ktype)+'.*.*.*')):
        yield open(instance_file).read()
        
    
def main():
    """
    main for run_instances.py
    """
    if len(sys.argv) != 4:
        usage(os.path.basename(sys.argv[0]))

    knapsack_program  = sys.argv[1]
    knapsack_options = sys.argv[2]
    instances_dir = sys.argv[3]


    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    sys.stdout.write("# Run as: " + " ".join(sys.argv) + "\n")
    sys.stdout.write("# at: " + timestamp + "\n")
    sys.stdout.write("# by: " + getpass.getuser() + "\n")
    sys.stdout.write("# on: " + os.popen("uname -a").read())


    for ktype in [1,2,3,4,5]:
        elapsed_list = []
        exceeded = False
        for problem in instance_generator(instances_dir, ktype):
            for i in xrange(NUM_REPETITIONS):
                (in_fh, out_fh, err_fh) = os.popen3(knapsack_program + ' ' + knapsack_options)
                in_fh.write(problem)
                in_fh.close()
                results = out_fh.read()
                errors = err_fh.read()
                if (len(errors) > 0):
                   sys.stderr.write("ERROR : " +  errors + "\n")
                   break
                nlresults = results.split('\n')
                instrument_line = None
                if (len(nlresults) > 2):
                   # get the INSTRUMENT summary line
                   if nlresults[0].split()[0] == "INSTRUMENT":
                      instrument_line = nlresults[0]
                   results = nlresults[1] + '\n' # normal output is last line
                sys.stdout.write(results)
                sresults = results.split()
                if sresults[0] == "LIMIT":
                    description = ' '.join(sresults[3:])
                    exceeded = True
                    break
                ttime = int(sresults[3])
                etime = int(sresults[4])
                flags = sresults[5]
                name = ' '.join(sresults[6:])

                elapsed_list.append(etime)

                if instrument_line:
                   sys.stdout.write(instrument_line + ' ' + flags + ' ' + name+'\n')

        # STATS: mean_elapsed median_elapsed max_elapsed flags poblem-name
        if exceeded:
            sys.stdout.write("STATS: LIMIT EXCEEDED " + description + '\n')
        else:
            sys.stdout.write('STATS: %f %f %d %s %s\n' %
                (mean(elapsed_list), median(elapsed_list), max(elapsed_list),
                 flags, str(ktype))
                        )


if __name__ == "__main__":
    main()

