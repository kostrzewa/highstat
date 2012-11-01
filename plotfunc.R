source("collect_data.R")

plotfunc <- function(sample) {
  pdf(onefile=TRUE,file=paste(subdir,paste(sample,"expvals.pdf",sep="_"),sep="/"),width=8,height=8)
  pat <- paste(paste("^",sample,sep=""),"_",sep="")

  # find the subdirectories which refer to the current sample
  dirs <- list.files(topdir,pattern=pat,full.names=TRUE)

  data <- collect_data(dirs)

  # when the length of data is 1, none of the output.data in this sample was suitable 
  # either none of them existed or none had a sufficient number of lines
  skipall <- FALSE
  if( length(data) < 2 ) {
    skipall <- TRUE
  }

  if( !skipall ) {

    # extract the names of the subdirectories for this sample set
    labels <- row.names(data[[1]])
    labels <- append(labels,"reference")

    # extract the "addon" from the labels (ie. serial, mpi_hs etc..)
    shortlabels <- extract_addons(labels,sample)

    # we use the same plotting function for the plaquette and rectangle so
    # we collect the values in temporary variables so we have to write the 
    # plotting routine only once
    for( k in c(1:4) ) {
      skip <- FALSE
      if( k == 1 ) {
        ar <- data[[1]]$ar  
        value <- data[[1]]$plaq
        dvalue <- data[[1]]$dplaq
        valuemed <- data[[1]]$plaqmed
        tauint <- data[[1]]$plaqtauint
        dtauint <- data[[1]]$plaqdtauint
        ref <- reference[sample,]$plaq
        dref <- reference[sample,]$dplaq
        name <- "plaquette"
        postname <- "exp. value"
      } else if ( k == 2 ) {
        ar <- data[[1]]$ar
        value <- data[[1]]$rec
        dvalue <- data[[1]]$drec
        valuemed <- data[[1]]$recmed
        tauint <- data[[1]]$rectauint
        dtauint <- data[[1]]$recdtauint
        ref <- reference[sample,]$rec
        dref <- reference[sample,]$drec
        name <- "rectangle"
        postname <- "exp. value"
      } else if ( k == 3 ) {
        ar <- NA
        value <- data[[1]]$trajtime
        dvalue <- data[[1]]$dtrajtime
        valuemed <- data[[1]]$trajtimemed
        tauint <- NA
        dtauint <- NA
        ref <- NA
        dref <- NA
        name <- "trajectory time"
        postname <- "mean"
      } else if ( k == 4 ) {
        ar <- NA
        value <- data[[1]]$cgitnum
        dvalue <- data[[1]]$dcgitnum
        valuemed <- data[[1]]$cgitnummed
        tauint <- NA
        dtauint <- NA
        ref <- NA
        dref <- NA
        name <- "CG iterations"
        postname <- "mean"
      }

      # do not attempt to plot rectangle expectation values for samples
      # which do not contain rectangle gauge terms
      if( k == 2 && sample %in% norectsamples ) {  
        skip = TRUE
      }
    
      # plot the expectation value, skip if this sample does not contain
      # does not contain the particular value (e.g. rectangle in the test above)
      if( !skip ) {
        par(mar=c(6.1,6.1,4.1,1.0))

        title <- paste(sample,paste(name,postname))
        title <- paste(title,"\n")
        title <- paste(title,topdirname)

        plotwitherror(x = c(1:length(value)), xlim=c(1,(length(value)+1)),
          y = value, dy=dvalue, 
          main=title, 
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
        text(c(1:(length(value)+1)),par("usr")[3], labels = shortlabels, 
          srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)


        if( k %in% c(1,2) ) {
          # display autocorrelation times on the plot
          rounded <- round(tauint,digits=5)
          for( j in 1:length(value) ) {
            tauintlabel <- bquote(tau[int] == .(rounded[j]) )
            text(x=j,y=par("usr")[3],labels=tauintlabel,srt=50,
              adj = c(-0.1,-0.1), xpd = TRUE, cex=.9) 
          }

          # display acceptance rate 
          for( j in 1:length(value) ) {
            ARlabel <- bquote(AR == .(ar[j]) )
            text(x=j,y=par("usr")[3],labels=ARlabel,srt=50,
              adj = c(-0.3,-2.0), xpd = TRUE, cex = .9)
          }
        }
      }

    } # loop over k, plotting plaquette, rectangle, CG time etc..

    # data contains the statistical results in [[1]] and for each output.data 
    # one plaquette and one rectangle history in [[2*j]] and [[2*j+1]] respectively
    # the total number of histories for this sample is thus:
    number <- ((length(data)-1)/2)

    title <- paste(sample,"initial plaquette history\n")
    title <- paste(title,topdirname)

    # plot initial plaquette histories now
    plot(xlab="",ylab="plaquette",y=data[[2]][1:100],x=c(1:100), 
      main=title,xlim=c(1,50*number),xaxt="n",lty=1,type="l")
    for(j in 1:number ) {
      lines(data[[2*j]][1:100], x=c(1:100)+(j-1)*50) 
    }
    # add a vertical grid to make identifying runs easier 
    abline(v=(seq(1,50*length(value),50)), col="lightgray",lty="dashed")

    text(x=seq(1,50*number,50),par("usr")[3], labels = shortlabels[1:number], 
      srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
    # close the expvals file for the current sample
    dev.off();

    # plot full plaquette histories in a separate file
    pdf(onefile=TRUE,file=paste(subdir,paste(sample,"pl_history.pdf",sep="_"),sep="/"),width=8,height=8)
    for( j in 1:number ) {
      title <- paste(labels[j],"plaquette history\n")
      title <- paste(title,topdirname)
      plot(data[[2*j]],x=c(1:length(data[[2*j]])),type="l",lty=1,main=title,xlab="trajectory",ylab="plaquette")
    }
    dev.off();

    skip <- FALSE
    for( j in c(1:length(norectsamples)) ) {
      if( sample == norectsamples[j] ) {
        skip = TRUE
      }
    }
    
    if( !skip ) {
      # plot rectangle histories if this sample has a rectangle gauge part
      pdf(onefile=TRUE,file=paste(subdir,paste(sample,"rec_history.pdf",sep="_"),sep="/"),width=8,height=8)
      for( j in 1:number ) {
        plot(data[[2*j+1]],x=c(1:length(data[[2*j+1]])),type="l",lty=1,main=paste(labels[j],"rectangle history"),xlab="trajectory",ylab="rectangle")
      }
      dev.off();
    }
  }
}
