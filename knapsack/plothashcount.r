###############################################################################
#
# plothashcount.r - plot hashcount graph from summarize_instrument output
#
# File:    plothashcount.r
# Author:  Alex Stivala
# Created: May 2009
#
#
#
# Plot hashcount graph from summarize_instrument.sh output
#
# Usage:
#       R --vanilla --slave -f plothashcount.r --args rtabfile
#
#    rtabfile is name of 3-column (threads, re, hc ) .instrument.rtab
#      file created by summarize_instrument.sh
#    output is PostScript file nnamed foo.hashcount.eps where 
#    foo.instrument.rtab was the input file
#
# $Id: plothashcount.r 2526 2009-06-15 06:23:41Z astivala $
# 
###############################################################################

rtabfile <- commandArgs(trailingOnly=TRUE)
instrtab <- read.table(rtabfile,header=TRUE)
instrtab <- instrtab[sort.list(instrtab$threads),] # sort by threads ascending


# EPS suitable for inserting into LaTeX
postscript(sub('[.]instrument.rtab$','.hashcount.eps',rtabfile),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

x <- instrtab$threads
y <- instrtab$hc
plot(x, y, xlab="threads", ylab="total computations (h)")
lines(x,y)
dev.off()

