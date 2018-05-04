###############################################################################
#
# plottestspeedup.r - plot speedup graph from timetests.sh output
#
# File:    plottestspeedup.r
# Author:  Alex Stivala
# Created: May 2009
#
# Plot speedup graph from timetests.sh output
#
# Usage:
#       R --vanilla --slave -f plottestspeedup.r --args rtabfile
#
#    rtabfile is name of 3-colmn (threads, iteration, ms) .rtab
#      file created by timetests.sh
#    output is PostScript file named foo.eps where foo.rtab was the input file
#
# 
# Requires the gplots package from CRAN to draw error bars.
#
# $Id: plottestspeedup.r 2492 2009-06-07 06:30:53Z astivala $
# 
###############################################################################

library(gplots)

rtabfile <- commandArgs(trailingOnly=TRUE)
timetab <- read.table(rtabfile,header=TRUE)

maxthreads = max(timetab$threads)



# time for "0 threads", the baseline (NB not the multithread program on 1
# thread, but a version compiled with no thread code at all)
# (speedup is relative to this )
basetime = mean(subset(timetab, threads==0)$ms)


# EPS suitable for inserting into LaTeX
postscript(sub('[.]rtab$','.eps',rtabfile),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

x <- subset(timetab, iteration==1)$threads
means <- c()
stdev <- c()
ciw <- c()
for (i in x) { 
  # FIXME this code is ugly and inefficient, should use sapply or something
  means <-  c(means, mean(basetime/subset(timetab, threads==i)$ms))
  thisstdev <-  sqrt(var(basetime/subset(timetab, threads==i)$ms))
  stdev <-  c(stdev, thisstdev)
  n <- length(subset(timetab, threads==i)$ms)
  ciw <- c(ciw, qt(0.975, n) * thisstdev / sqrt(n))
}

means
stdev
ciw

plotCI(x, y=means, uiw=ciw, xlab="threads", ylab="speedup" )
#     main=paste('httslftest',rtabfile) ) 
lines(x,means)
dev.off()
warnings()
