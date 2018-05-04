###############################################################################
#
# plotspeeduptabtype2.r - plot speedup graphs for 2 table sizes from mkrtab.sh output
#
# File:    plotspeeduptabtype2.r
# Author:  Alex Stivala
# Created: May 2009
#
# Plot speedup graph from mkrtab.sh output
#
# Usage:
#       R --vanilla --slave -f plotspeeduptabtype2.r --args rtabfile
#
#    rtabfile is name of 3-colmn (threads, iter, ms) .rtab
#      file created by mkrtab.sh
#    output is PostScript file named foo-tabtype2.eps where foo.rtab was the input file
#
# 
# Requires the gplots package from CRAN to draw error bars.
#
# $Id: plotspeeduptabtype2.r 3152 2009-12-28 00:07:23Z alexs $
# 
###############################################################################

library(gplots)

rtabfile <- commandArgs(trailingOnly=TRUE)
timetab <- read.table(rtabfile,header=TRUE)
timetab <- timetab[sort.list(timetab$threads),] # sort by threads ascending

maxthreads = max(timetab$threads)

rtab2file <- sub('[.]rtab$', '.httslf.rtab', rtabfile)
timetab2 <- read.table(rtab2file,header=TRUE)
timetab2 <- timetab2[sort.list(timetab2$threads),] # sort by threads ascending


# time for "0 threads", the baseline (NB not the multithread program on 1
# thread, but a version compiled with no thread code at all)
# (speedup is relative to this )
basetime = mean(subset(timetab, threads==0)$ms)
basetime2 = mean(subset(timetab2, threads==0)$ms)


# EPS suitable for inserting into LaTeX
psfile <- sub('[.]rtab$','-tabtype2.eps',rtabfile)
postscript(psfile,
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

#x <- subset(timetab, iter==1)$threads
#means <- c()
#stdev <- c()
#ciw <- c()
#for (i in x) { 
#  # FIXME this code is ugly and inefficient, should use sapply or something
#  means <-  c(means, mean(basetime/subset(timetab, threads==i)$ms))
#  thisstdev <-  sqrt(var(basetime/subset(timetab, threads==i)$ms))
#  stdev <-  c(stdev, thisstdev)
#  n <- length(subset(timetab, threads==i)$ms)
#  ciw <- c(ciw, qt(0.975, n) * thisstdev / sqrt(n))
#}
#plotCI(x, y=means, uiw=ciw, xlab="threads", ylab="speedup" , lty=1,pch=20,col='red')
#lines(x,means,pch=20,col='red')

x <- timetab$threads
y <- basetime / timetab$ms
plot(x, y, uiw=NULL, xlab="threads", ylab="speedup" , lty=1,pch=20,col='red')
lines(x,y,pch=20,col='red')

x <- subset(timetab2, iter==1)$threads
means <- c()
stdev <- c()
ciw <- c()
for (i in x) { 
  # FIXME this code is ugly and inefficient, should use sapply or something
  means <-  c(means, mean(basetime2/subset(timetab2, threads==i)$ms))
  thisstdev <-  sqrt(var(basetime2/subset(timetab2, threads==i)$ms))
  stdev <-  c(stdev, thisstdev)
  n <- length(subset(timetab2, threads==i)$ms)
  ciw <- c(ciw, qt(0.975, n) * thisstdev / sqrt(n))
}
plotCI(x, y=means, uiw=ciw, add=TRUE, lty=2,pch=21,col='green')
lines(x,means,lty=2,col='green')

legend('topleft', col=c('red','green'), lty=c(1,2), pch=c(20,21), legend=c('open addressing hash table','separate chaining hash table'), bty='n')
dev.off()

