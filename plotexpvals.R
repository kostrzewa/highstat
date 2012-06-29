# plotexpvals is a function which uses Carsten Urbach's hadron library
# to compute expectation values of the plaquette and the rectangle
# for the set of sample files of the tmLQCD software package
# under the hood it also computes the optimal autocorrelation time
# accoring to U. Wollf's method

# it produces two to three files for each sample file,
# the first file contains a plot for for the plaquette expectation value, one for
# the rectangle expectation value (unless it contains no rectangle gauge contribution), 
# and one plot each for the first 100 trajectories of the plaquette and rectangle
# histories (as available)
# the second and third files contain complete histories of the plaquette and rectangle part

# the argument "topdir" must be the path to a directory which contains
# a set of subdirectories
# these subdirectories in turn relate to different runs of the different
# sample input files
# each subdirectory must contain "output.data" for the given run

# the names of the subdirectories should begin with the strings in the 
# "samples" vector below, separated from the rest of the filename
# by an underscore ('_')

# "norectsamples" is a subset of samples which do not contain a rectangle
# contribution in their gauge action

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

# because of the way that the subfunctions are called this should be run
# in the directory containing the R files

library(hadron)

source("expvals.R")

subdir <- "plots"

plotexpvals <- function(topdir) {
  samples <- c("hmc0","hmc1","hmc2","hmc3","hmc4","hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
  
  # these samples do not contain a rectangular gauge part
  norectsamples <- c("hmc0","hmc1","hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
 
  # these are reference values for the plaquette and rectangle expectation value 
  reference <- read.table("reference.dat", fill=TRUE)
   
  for( i in 1:length(samples) ) {
    pdf(onefile=TRUE,file=paste(subdir,paste(samples[i],"expvals.pdf",sep="_"),sep="/"),width=8,height=8)
    skip <- FALSE
    pat <- paste(paste("^",samples[i],sep=""),"_",sep="")

    # find the subdirectories which refer to the current sample
    dirs <- list.files(topdir,pattern=pat,full.names=TRUE)
    
    evals <- expvals(dirs)

    # extract the names of the subdirectories for this sample set
    labels <- row.names(evals[[1]])
    labels <- append(labels,"reference")

    # we use the same plotting function for the plaquette and rectangle so
    # we collect the values in temporary variables so we have to write the 
    # plotting routine only once
    for( k in c(1:2) ) {
      if( k == 1 ) {  
        value <- evals[[1]]$plaq
        dvalue <- evals[[1]]$dplaq
        valuemed <- evals[[1]]$plaqmed
        tauint <- evals[[1]]$plaqtauint
        dtauint <- evals[[1]]$plaqdtauint
        ref <- reference[samples[i],]$plaq
        dref <- reference[samples[i],]$dplaq
        name <- "plaquette"
      } else if ( k == 2 ) {
        value <- evals[[1]]$rec
        dvalue <- evals[[1]]$drec
        valuemed <- evals[[1]]$recmed
        tauint <- evals[[1]]$rectauint
        dtauint <- evals[[1]]$recdtauint
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

   
      # plot the plaquette / rectangle expectation value, skip if this sample
      # does not contain a rectangle gauge part (the for loop over k)
      if( !skip ) {
        par(mar=c(8.1,9.1,4.1,1.0))

        plotwitherror(x = c(1:length(value)), xlim=c(1,(length(value)+1)),
         y = value, dy=dvalue, 
         main=paste(samples[i],paste(name,"exp. value")), 
         xlab="",ylab=name,xaxt="n",pch=16)
     
        axis(1,labels=FALSE,tick=FALSE)

        plotwitherror(add=TRUE, x = length(value)+1, xlim=c(1,(length(value)+1)),
          y = ref, dy=dref, 
          col="dark red",pch=16,xaxt="n")
        
        # add a vertical grid to make identifying runs easier 
        abline(v=(seq(1,(length(value)+1),1)), col="lightgray",lty="dashed")

        # display median value of the quantity
        points(c(1:length(valuemed)),valuemed) 

        # print names of data sets
        text(c(1:(length(value)+1)),par("usr")[3], labels = labels, 
          srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

        # display autocorrelation times on the plot
        rounded <- round(tauint,digits=5)
        for( k in 1:length(value) ) {
          tauintlabel <- bquote(tau[int] == .(rounded[k]) )
          text(x=k,y=par("usr")[3],labels=tauintlabel,srt=50,
            adj = c(-0.1,-0.1), xpd = TRUE, cex=.9) 
        }
      }
    }

    # evals contains the statistical results in [[1]] and for each output.data 
    # one plaquette and one rectangle history in [[2*j]] and [[2*j+1]] respectively
    # the total number of histories for this sample is thus:
    number <- ((length(evals)-1)/2)

    # plot initial plaquette histories now
    plot(xlab="",ylab="plaquette",y=evals[[2]][1:100],x=c(1:100), 
      main=paste( samples[i], "plaquette history"),xlim=c(1,50*number),xaxt="n",lty=1,type="l")
    for(j in 2:number ) {
      lines(evals[[2*j]][1:100], x=c(1:100)+(j-1)*50) 
    }
    text(x=seq(1,50*number,50),par("usr")[3], labels = labels[1:number], 
     srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
    # close the expvals file for the current sample
    dev.off();

    # plot full plaquette histories in a separate file
    pdf(onefile=TRUE,file=paste(subdir,paste(samples[i],"pl_history.pdf",sep="_"),sep="/"),width=8,height=8)
    for( j in 1:number ) {
      plot(evals[[2*j]],x=c(1:length(evals[[2*j]])),type="l",lty=1,main=paste(labels[j],"plaquette history"),xlab="trajectory",ylab="plaquette")
    }
    dev.off();

    skip <- FALSE
    for( j in c(1:length(norectsamples)) ) {
      if( samples[i] == norectsamples[j] ) {
        skip = TRUE
      }
    }
    
    if( !skip ) {
      # plot rectangle histories if this sample has a rectangle gauge part
      pdf(onefile=TRUE,file=paste(subdir,paste(samples[i],"rec_history.pdf",sep="_"),sep="/"),width=8,height=8)
      for( j in 1:number ) {
        plot(evals[[2*j+1]],x=c(1:length(evals[[2*j+1]])),type="l",lty=1,main=paste(labels[j],"rectangle history"),xlab="trajectory",ylab="rectangle")
      }
      dev.off();
    }
  }
}
 
