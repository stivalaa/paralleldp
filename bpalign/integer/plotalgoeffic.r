###############################################################################
#
# plotalgoeffic.r - plot algorithmic efficiency graphs 
#
# File:    plotalgoeffic.r
# Author:  Alex Stivala
# Created: May 2009
#
# Plot algorithmic efficiency graphs from results2rtab.sh output
#
# Usage:
#       R --vanilla --slave -f plotalgoeffic.r
#
#    output is PostScript files named hostname-hashcountoverreuse.eps and
#    hostname-hashcountmulti.eps  and
#    hostname-hashcountincreasemulti.eps
#    where hostmame is name of the test machine.
#
# 
#
# $Id: plotalgoeffic.r 2649 2009-07-13 05:09:58Z astivala $
# 
###############################################################################



#
# globals
#

colorvec=c('deepskyblue4','brown','red','turquoise','blue','purple','green','cyan','gray20','magenta','darkolivegreen2','midnightblue','magenta3','darkseagreen','violetred3','darkslategray3')
ltyvec=c(1,2,4,5,6,1,2,1,5,6,1,2,4,5,6,1,2)
pchvec=c(20,21,22,23,24,25)
namevec=c('RNA basepair matrix alignment')
directoryvec=c('.')
hostnamevec=c('mundara')
numcoresvec=c( 32      )

#
# utility functions
#

# reutrn every 2nd element of a vector (starting at first element)
every2nd <- function(vec) return(vec[which(vec == vec) %% 2 == 1])

#
# graphing functions
#

plothashcount <- function(hostname, numcores)
{

  # EPS suitable for inserting into LaTeX
  # NB use of - not . as filename speartor for PJS use of pdflatex
  outfile <- paste(hostname,'hashcountmulti',sep='-')
  outfile <- paste(outfile,'eps',sep='.')
  postscript(outfile,
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)

  plot(0,0,xlab="threads", ylab="total computations (h)", type='n',
       xlim=c(1,numcores), ylim=c(500000,4.0e+06) )
       #ylim=c(2e+08, 1.7e+10) )

  k = 1
  for (j in 1:length(directoryvec)) {

#    for (hasrand in c('','_norandomization'))  {
    for (hasrand in c(''))  {
      rtabfile <- paste(directoryvec[j],hostname,sep='/')
      rtabfile <- paste(rtabfile, hasrand,sep='')
      rtabfile <- paste(rtabfile,'rtab',sep='.')

      instrumenttab <- read.table(rtabfile,header=TRUE)
      instrumenttab <- instrumenttab[sort.list(instrumenttab$threads),] # sort by threads ascending
  
      x <- instrumenttab$threads
      y <- instrumenttab$hc
      points(x, y,col=colorvec[k],lty=ltyvec[k],pch=pchvec[k])
      lines(x,y,col=colorvec[k],lty=ltyvec[k])         
      k  = k + 1
    }
  }
  legend('topleft', col=colorvec, lty=ltyvec, pch=pchvec, legend=namevec)  
  dev.off()
}



plothashcountincrease <- function(hostname, numcores)
{

  # EPS suitable for inserting into LaTeX
  # NB use of - not . as filename speartor for PJS use of pdflatex
  outfile <- paste(hostname,'hashcountincreasemulti',sep='-')
  outfile <- paste(outfile,'eps',sep='.')
  postscript(outfile,
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)

  plot(0,0,xlab="threads", ylab="increase in total computations (h/h1)", type='n',
       xlim=c(1,numcores), ylim=c(1.0, 8.0) )

  k = 1
  for (j in 1:length(directoryvec)) {

    for (hasrand in c(''))  {  # only WITH randomization
      rtabfile <- paste(directoryvec[j],hostname,sep='/')
      rtabfile <- paste(rtabfile, hasrand,sep='')
      rtabfile <- paste(rtabfile,'rtab',sep='.')

  
      instrumenttab <- read.table(rtabfile,header=TRUE)
      instrumenttab <- instrumenttab[sort.list(instrumenttab$threads),] # sort by threads ascending

      x <- instrumenttab$threads
      y <- instrumenttab$hc / instrumenttab$hc[1]
      points(x, y,col=colorvec[k],lty=ltyvec[k],pch=pchvec[k])
#      lines(x,y,col=colorvec[k],lty=ltyvec[k])
      model <- lm(y ~ x)
      abline(model, col=colorvec[k],lty=ltyvec[k],pch=pchvec[k])
 #     k  = k + 2  # add 2 not 1 so get seem lty/col/pch as other graphs
       k = k + 1
    }
  }

  colorvec2 <- every2nd(colorvec)
  ltyvec2 <- every2nd(ltyvec)
  pchvec2 <- every2nd(pchvec)
  namevec2 <- every2nd(namevec)
  legend('topleft', col=colorvec2, lty=ltyvec2, pch=pchvec2, legend=namevec2)
  dev.off()
}


plotalgoeffic <- function(hostname, numcores)
{

  # EPS suitable for inserting into LaTeX
  # NB use of - not . as filename speartor for PJS use of pdflatex
  outfile <- paste(hostname,'hashcountoverreusemulti',sep='-')
  outfile <- paste(outfile,'eps',sep='.')
  postscript(outfile,
             onefile=FALSE,paper="special",horizontal=FALSE, 
             width = 9, height = 6)

  plot(0,0,xlab="threads", ylab="(reuse/hashcount)", type='n',
       xlim=c(1,numcores), ylim=c(0,100) )

  k = 1
  for (j in 1:length(directoryvec)) {

#    for (hasrand in c('','_norandomization'))  {
    for (hasrand in c(''))  {
      rtabfile <- paste(directoryvec[j],hostname,sep='/')
      rtabfile <- paste(rtabfile, hasrand,sep='')
      rtabfile <- paste(rtabfile,'rtab',sep='.')

      instrumenttab <- read.table(rtabfile,header=TRUE)
      instrumenttab <- instrumenttab[sort.list(instrumenttab$threads),] # sort by threads ascending
      
      x <- instrumenttab$threads
      y <- instrumenttab$re / instrumenttab$hc
      points(x, y,col=colorvec[k],lty=ltyvec[k],pch=pchvec[k])
      lines(x,y,col=colorvec[k],lty=ltyvec[k])         
      k = k + 1
    }
  }
    
  legend('topleft', col=colorvec, lty=ltyvec, pch=pchvec, legend=namevec)  

  dev.off()
}

#
# main
#
  
for (i in 1:length(hostnamevec)) {
  plothashcount(hostnamevec[i], numcoresvec[i])
  plothashcountincrease(hostnamevec[i], numcoresvec[i])
  plotalgoeffic(hostnamevec[i], numcoresvec[i])
}

warnings()
