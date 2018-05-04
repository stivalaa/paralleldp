###############################################################################
#
# plottestspeedupmulti2.r - plot speedup graph from timetests.sh output
#
# File:    plottestspeedupmulti2.r
# Author:  Alex Stivala
# Created: May 2009
#
# Plot speedup graphs from timetests.sh output for different meethods on
# one graph for each platform.
#
# This is a modified version for mundara (SPARC) and
# mungera (PPC64) results that
# those platforms don't run most methods (they seem to be Intel only).
#
# Usage:
#       R --vanilla --slave -f plottestspeedupmulti2.r
#
#    output is PostScript file named hostname.testspeedupmulti2.eps
#    where hostmame is name of the test machine.
#
# 
# Requires the gplots package from CRAN to draw error bars.
#
# $Id: plottestspeedupmulti2.r 3126 2009-12-25 02:57:36Z alexs $
# 
###############################################################################

library(gplots)


#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
pchvec=c(20,21,22,23,24,25)
namevec=c('open addressing','separate chaining')
fileprefixvec=c('oahttslftest','')

hostnamevec=c('mundara','mungera')


plotspeedup <- function(hostname)
{

  # EPS suitable for inserting into LaTeX
  postscript(paste(hostname,'testspeedupmulti2','eps',sep='.'),
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)


  for (j in 1:length(fileprefixvec)) {

    if (fileprefixvec[j] == '')
      rtabfile <- paste(hostname,'rtab',sep='.')
    else
      rtabfile <- paste(hostname,fileprefixvec[j],'rtab',sep='.')

    
    timetab <- read.table(rtabfile,header=TRUE)

    maxthreads = max(timetab$threads)
    
    # time for "0 threads", the baseline (NB not the multi2thread program on 1
    # thread, but a version compiled with no thread code at all)
    # (speedup is relative to this )
    basetime = mean(subset(timetab, threads==0)$ms)

    x <- subset(timetab, iteration==1)$threads
    means <- c()
    stdev <- c()
    ciw <- c()

    print(hostname)

    for (i in x) { 
      # FIXME this code is ugly and inefficient, should use sapply or something
      means <-  c(means, mean(basetime/subset(timetab, threads==i)$ms))
      thisstdev <-  sqrt(var(basetime/subset(timetab, threads==i)$ms))
      stdev <-  c(stdev, thisstdev)
      n <- length(subset(timetab, threads==i)$ms)
      ciw <- c(ciw, qt(0.975, n) * thisstdev / sqrt(n))
    }
  
    print(means)

    plotCI(x, y=means, uiw=ciw, xlab="threads", ylab="speedup",
           add = (j > 1),
           col=colorvec[j],pch=pchvec[j])
    lines(x,means,col=colorvec[j],lty=ltyvec[j])
  }
  legend('topleft', col=colorvec, lty=ltyvec, pch=pchvec, legend=namevec, bty='n')  
  dev.off()
}

#
# main
#
  
for (hostname in hostnamevec) {
  plotspeedup(hostname)
}

