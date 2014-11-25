# This file is part of the "highstat" R script set
# Copyright (C) 2012  Bartosz Kostrzewa

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


source("collect_data.R")
source("extract_addons.R")
source("plot_timeseries.R")

plotfunc <- function(samp) {
  pat <- paste(paste("^",samp,sep=""),"_",sep="")

  # find the subdirectories which refer to the current sample
  dirs <- list.files(topdir,pattern=pat,full.names=TRUE)

  skipall <- FALSE

  if( length(dirs) == 0 ) {
    print(sprintf("No data found for this sample file. Skipping %s!",samp))
    skipall <- TRUE
  } else {
    data <- collect_data(dirs)
  }

  # when the length of data is 1, none of the output.data in this sample was suitable 
  # either none of them existed or none had a sufficient number of lines
  if( !skipall && length(data) < 2 ) {
    skipall <- TRUE
  }

  reftime <- matrix( rep(NA,length(execs)), ncol=length(execs), byrow=TRUE )
  colnames(reftime) <- execs
  rownames(reftime) <- samp
  
  if( !skipall ) {
    save(data,file=sprintf("%s.Rdata",samp))
    pdf(onefile=TRUE,file=paste(subdir,paste(samp,"expvals.pdf",sep="_"),sep="/"),width=8,height=8)

    # we use the same plotting function for all observables
    # we collect the values in temporary variables so we have to write the 
    # plotting routine only once
    for( k in c(1:4) ) {
      hasref <- FALSE
      if( k == 1 ) {
        ar <- data[[1]]$ar  
        value <- data[[1]]$plaq
        dvalue <- data[[1]]$dplaq
        valuemed <- data[[1]]$plaqmed
        tauint <- data[[1]]$plaqtauint
        dtauint <- data[[1]]$plaqdtauint
        ref <- reference[samp,]$plaq
        dref <- reference[samp,]$dplaq
        name <- "plaquette"
        postname <- "exp. value"
        hasref <- TRUE
      } else if ( k == 2 ) {
        ar <- data[[1]]$ar
        value <- data[[1]]$rec
        dvalue <- data[[1]]$drec
        valuemed <- data[[1]]$recmed
        tauint <- data[[1]]$rectauint
        dtauint <- data[[1]]$recdtauint
        ref <- reference[samp,]$rec
        dref <- reference[samp,]$drec
        name <- "rectangle"
        postname <- "exp. value"
        hasref <- TRUE
      } else if ( k == 3 ) {
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
        hasref <- TRUE
      } else if ( k == 4 ) {
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
        hasref <- FALSE

        # we write runtimes.csv so we know how long a run takes
        # output the trajectory time to use as a timebase
        # multiply by 1.2 as a fudge factor because the zeuthen cluster has
        # unpredictable performance! 
        # note: 1.2 * 1000 * trajtime / 3600 = (1.2 * hours for 1000 trajectories)
        # so instead of dividing by 3.6, we divide by 3 to get hours
        lreftime <- matrix(value/3,ncol=length(value),byrow=TRUE)
        colnames(lreftime) <- shortlabels[1:length(shortlabels)-1]
        rownames(lreftime) <- samp
       
        # merge the table full of NA with the table of actual measurements
        # when completed, this will become runtimes.csv
        reftime <- merge(lreftime,reftime,all.x=TRUE,all.y=FALSE,by=intersect(colnames(reftime),colnames(lreftime)))
        rownames(reftime) <- samp
      }
      
      # extract the names of the subdirectories for this sample set
      # these will become the tick labels on the plot after shortening
      labels <- row.names(data[[1]])
      if(hasref) {
        labels <- append(labels,"reference")
      }
      # extract the "addon" from the labels (ie. serial, mpi_hs etc..)
      # these will be the ticklabels on the plot
      shortlabels <- extract_addons(labels,samp)

      # plot the expectation value, skip if this sample does not contain
      # the particular value (the data collection function would have set it
      # to NA in this case )
      if(!( NA %in% value)) {
        par(mar=c(6.1,6.1,4.1,1.0))

        title <- paste(samp,paste(name,postname))
        title <- paste(title,"\n")
        title <- paste(title,name)
        
        miny <- min(value)-0.22*(max(value)-min(value))
        if( miny > min(value)-max(dvalue) ) {
          miny <- min(value)-1.1*max(dvalue)
        }

        maxy <- max(value)+max(dvalue)

        if( hasref ) {
          # set up plot
          plotwitherror(x = c(1:length(value)), xlim=c(1,(length(value)+1)),
            y = value, dy=dvalue, ylim=c(miny,maxy), 
            main=title, 
            xlab="",ylab=name,xaxt="n",pch=16, type='n')
          # add a vertical grid to make identifying runs easier
          abline(v=(seq(1,(length(value)+1),1)), col="lightgray",lty="dashed")
          # now plot
          plotwitherror(x = c(1:length(value)), y = value, dy=dvalue, pch=16, rep=TRUE)
          plotwitherror(rep=TRUE, x = length(value)+1, xlim=c(1,(length(value)+1)),
            y = ref, dy=dref, 
            col="dark red",pch=16,xaxt="n")
        } else {
          # set up plot
          plotwitherror(x = c(1:length(value)), xlim=c(1,length(value)),
            y = value, dy=dvalue, ylim=c(miny,maxy), 
            main=title, 
            xlab="",ylab=name,xaxt="n",pch=16, type='n')
          # add a vertical grid to make identifying runs easier
          abline(v=(seq(1,length(value),1)), col="lightgray",lty="dashed")
          # and plot
          plotwitherror(x = c(1:length(value)), y = value, dy=dvalue, pch=16, rep=TRUE)
        }
     
        axis(1,labels=FALSE,tick=TRUE,tck=-0.007,at=c(1:(length(value)+1)))

        # display median value of the quantity
        points(c(1:length(valuemed)),valuemed) 

        # print short names of data sets (also called "addons" elsewhere)
        # if there is no reference, do not print the final label which is "reference"
        # these are at the locations of the tick labels and rotated for legibility
        text(c(1:(length(shortlabels))),par("usr")[3], labels = shortlabels, 
          srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

        # for the plaquette and rectangle we would like to display 
        # the acceptance rate and autocorrelation time on the plot
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
              adj = c(-0.3,-1.6), xpd = TRUE, cex = .9)
          }
        }
      }

    } # loop over k, plotting plaquette, rectangle, CG time etc..
    # data contains the statistical results in [[1]] and for each output.data 
    # and starting from [[2]] the histories of the plaquette, rectangle, CG iterations
    # and trajectory times
    # so : data[[2+4*(j-1)+x]], x \in [0,1,2,3]

    # the total number of histories for this sample is thus:
    number <- ((length(data)-1)/4)

    title <- paste(samp,"initial plaquette history\n")
    title <- paste(title,runname)

    # plot initial plaquette histories now
    plot(xlab="",ylab="plaquette",y=data[[2]][1:100],x=c(1:100), 
      main=title,xlim=c(1,50*number),xaxt="n",lty=1,type="l")
    for(j in 1:number ) {
      lines(data[[2+4*(j-1)]][1:100], x=c(1:100)+(j-1)*50) 
    }
    # add a vertical grid to make identifying runs easier 
    abline(v=(seq(1,50*length(value),50)), col="lightgray",lty="dashed")

    text(x=seq(1,50*number,50),par("usr")[3], labels = shortlabels[1:number], 
      srt = 30, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

    axis(1,labels=FALSE,tick=TRUE,tck=-0.007,at=(seq(1,50*length(value),50)))

    # close the expvals file for the current sample
    dev.off();

    # plot full plaquette histories in a separate file
    pdf(onefile=TRUE,file=paste(subdir,paste(samp,"pl_history.pdf",sep="_"),sep="/"),width=8,height=8)
    for( j in 1:number ) {
      title <- paste(labels[j],"plaquette history\n")
      title <- paste(title,runname)
      # we start at trajectory 6 arbitrarily to get a smaller yrange
      plot(y=data[[2+4*(j-1)]][ 6:length(data[[2+4*(j-1)]]) ],x=c(6:length(data[[2+4*(j-1)]])),type="p",pch=".",main=title,xlab="trajectory",ylab="plaquette")
    }
    dev.off();

    # and plot the analyzed parts of the histoies 
    # in seperate plots with a histogram, a periodogram and the UWerr plots
    for( j in 1:number ) {
      max <- length(data[[2+4*(j-1)]])
      if(limit) max <- min+trajs
      title <- paste(labels[j],"\n")
      title <- paste(title,runname)

      plot_timeseries(dat=data[[2+4*(j-1)]][ min:max ], trange=c(min,max),
                      pdf.filename=paste(subdir,paste(labels[j],"plaquette_timeseries.pdf",sep="_"),sep='/'),
                      filelabel=title,name="plaquette",plotsize=5,ylab="<P>",
                      titletext=title,periodogram=TRUE)
    }
       
    if( ! (samp %in% norectsamples) ) {
      # plot rectangle histories if this sample has a rectangle gauge part
      pdf(onefile=TRUE,file=paste(subdir,paste(samp,"rec_history.pdf",sep="_"),sep="/"),width=8,height=8)
      for( j in 1:number ) {
        title <- paste(labels[j],"rectangle history\n")
        title <- paste(title,runname)
      # we start at trajectory 6 arbitrarily to get a smaller yrange
        plot(y=data[[2+4*(j-1)+1]][ 6:length(data[[2+4*(j-1)+1]]) ],x=c(6:length(data[[2+4*(j-1)+1]])),type="p",pch=".",main=title,xlab="trajectory",ylab="rectangle")
      }
      dev.off();

      # and plot the analyzed parts of the histoies 
      # in seperate plots with a histogram, a periodogram and the UWerr plots
      for( j in 1:number ) {
        max <- length(data[[2+4*(j-1)+1]])
        if(limit) max <- min+trajs
        title <- paste(labels[j],"\n")
        title <- paste(title,runname)
  
        plot_timeseries(dat=data[[2+4*(j-1)+1]][ min:max ], trange=c(min,max),
                        pdf.filename=paste(subdir,paste(labels[j],"rectangle_timeseries.pdf",sep="_"),sep='/'),
                        filelabel=title,name="rectangle",plotsize=5,ylab="<P>",
                        titletext=title,periodogram=TRUE)
      }
    } # if norectsamples
  }
  return(reftime)
}
