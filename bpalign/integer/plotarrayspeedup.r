###############################################################################
#
# plotarrayspeedup.r - plot speedup graph from results2rtab.sh output
#
# File:    plotarrayspeedup.r
# Author:  Alex Stivala
# Created: May 2009
#
#
#
# Plot speedup graph from results2rtab.sh output for array (rather than
# hashtable) implementation, with bottom-up time as a horizontal line
# to mark "break-even" point.
#
# Usage:
#       R --vanilla --slave -f plotarrayspeedup.r --args
#                                               rtabfile bottomup_rtabfile
#
#    rtabfile is name of 3-column (threads, iter, ms) .rtab
#      file created by results2rtab.sh
#    bottomup_rtablefile  is name of 3-column (threads, iter ms) .rtab file
#      created by results2rtab.sh for the bottom-up implemenmtation (only
#      0-thread value is relevant)
#    output is PostScript file nnamed foo.eps where foo.rtab was the input file
#
# $Id: plotarrayspeedup.r 3117 2009-12-24 00:28:28Z alexs $
# 
###############################################################################

library(gplots)

args <- commandArgs(trailingOnly=TRUE)

rtabfile <- args[1]
bottomuprtabfile <- args[2]

timetab <- read.table(rtabfile,header=TRUE)
timetab <- timetab[sort.list(timetab$threads),] # sort by threads ascending
bottomuprtab <- read.table(bottomuprtabfile,header=TRUE)
bottomuprtab <- subset(bottomuprtab, threads == 0)

# EPS suitable for inserting into LaTeX
postscript(sub('[.]rtab$','.eps',rtabfile),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

# time for "0 threads", the baseline (NB not the multithread program on 1
# thread, but a version compiled with no thread code at all)
# (speedup is relative to this )
basetime = mean(subset(timetab, threads==0)$ms)

x <- subset(timetab, iter==1)$threads
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

y <- means

bottomupspeedup <- basetime / bottomuprtab$ms[1]
plotCI(x, y=means, uiw=ciw, xlab="threads", ylab="speedup")
lines(x,means)
abline(h = bottomupspeedup, lty = 2)
dev.off()

maxspeedup = max(y)
maxspeedupindex = which(y == maxspeedup)
maxthreads = x[maxspeedupindex]
print(paste(format(bottomupspeedup, width=3, digits=3, nsmall=2),
            format(maxspeedup, width=3, digits=3,nsmall=2),
            format(maxthreads, width = 2),
            sep = ' '))

