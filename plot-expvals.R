# plotexpvals is a function which uses Carsten Urbach's hadron library
# to compute expectation values of the plaquette and the rectangle
# for the set of sample files of the tmLQCD software package

# it produces two plots for each sample file, one for the plaquette, the
# other for the rectangle

# the argument "topdir" must be the path to a directory which contains
# a set of subdirectories
# these subdirectories in turn relate to different runs of the different
# sample input files
# each subdirectory must contain "output.data" for the given run

# the names of the subdirectories should begin with the strings in the 
# "samples" vector below

# example directory structure
# topdir
# |
# -----> hmc0_old_serial
# |  |-> hmc0_new_serial 
# |  |-> hmc0_old_mpi
# |  |-> hmc0_new_mpi
# |
# -----> hmc1_[...]
# [...] 

library(hadron)

source("expvals.R")

plotexpvals <- function(topdir) {
  samples <- c("hmc0","hmc1","hmc2","hmc3","hmc4","hmc-cloverdet","hmc-tmcloverdet","hmc-tmcloverdetratio")
  #samples <- c("hmc0","hmc3")
  norectsamples <- c("hmc0","hmc1","hmc-cloverdet","hmc-tmcloverdet","hmc-tmcloverdetratio")
  
  reference <- read.table("reference.dat", fill=TRUE)
   
  pdf(onefile=TRUE,file="expvals.pdf",width=8,height=8)
  for( i in 1:length(samples) ) {
    skip <- FALSE
    pat <- paste("^",samples[i],sep="")
    dirs <- list.files(topdir,pattern=pat,full.names=TRUE)
    
    evals <- expvals(dirs)

    labels <- row.names(evals)
    labels <- append(labels,"reference")

    for( k in c(1:2) ) {
      if( k == 1 ) {  
        value <- evals$plaq
        dvalue <- evals$dplaq
        tauint <- evals$plaqtauint
        dtauint <- evals$plaqdtauint
        ref <- reference[samples[i],]$plaq
        dref <- reference[samples[i],]$dplaq
        name <- "plaquette"
      } else if ( k == 2 ) {
        value <- evals$rec
        dvalue <- evals$drec
        tauint <- evals$rectauint
        dtauint <- evals$recdtauint
        ref <- reference[samples[i],]$rec
        dref <- reference[samples[i],]$drec
        name <- "rectangle"

        # do not attempt to plot rectangle expectation values for samples
        # which do not contain rectangle gauge terms
        for( j in c(1:length(norectsamples)) ) {
          if( samples[i] == norectsamples[j] ) {
            skip = TRUE
          }
        }
      }

    
      if( !skip ) {
        par(mar=c(8.1,9.1,4.1,1.0))
        mp <- plotwitherror(x = c(1:length(evals$plaq)), xlim=c(1,(length(evals$plaq)+1)),
         y = value, dy=dvalue, 
         main=paste(samples[i],paste(name,"exp. value")), 
         xlab="",ylab=name,xaxt="n",pch=16)
      
        axis(1,labels=FALSE)

        plotwitherror(add=TRUE, x = length(evals$plaq)+1, xlim=c(1,(length(evals$plaq)+1)),
          y = ref, dy=dref, 
          col="dark red",pch=16)

        text(c(1:(length(evals$plaq)+1)),par("usr")[3], labels = labels, srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
      }
    }
  }
  dev.off();
}
 
