#!/usr/bin/env python
"""
  calltree2dot.py - convert -c output from knapsack_simple into call tree
                    in .dot format for GraphViz graph layout
  
  input on stdin, output on stdout

  example usage:

  knapsack_simple -c < tiny.in  2> dbgout
  fgrep '(' dbgout | fgrep -v '99999,99999' | calltree2dot.py > calltree.dot

  or (from ../pjs/src/):

  ./knapsack < tiny.in | grep 'dp:' | fgrep -v '99999,99999' | ../../knapsack/calltree2dot.py | dot -Tps > simple.ps
"""

import sys

def label2nodename(label):
    """
    convert label e.g. '(1,2)' to nodename e.g. 'n1c2'
    """
    lsplit = label.split(',')
    return 'n' + lsplit[0][1:] + 'c' + lsplit[1][:-1]

node_dict = {}
adjlist = []
for line in sys.stdin:
    sline = line.split()
    t1 = sline[1]
    t2 = sline[3]
    if not node_dict.has_key(t1):
        node_dict[t1] = True
    if not node_dict.has_key(t2):
        node_dict[t2] = True
    adjlist.append( (t1, t2) )

print 'graph {'
sys.stdout.write('node [shape=circle]')
for node in node_dict.iterkeys():
    sys.stdout.write('    ' + label2nodename(node) + ' [label="k'+node+'"]\n')
for edge in adjlist:
    sys.stdout.write('    ' + label2nodename(edge[0]) + ' -- ' + label2nodename(edge[1]) + '\n')
print '}'
