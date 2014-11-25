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


# given a vector of directories containing output.data, call readoutput
# to get expectation value and uwerr of plaquette and rectangular
# plaquette, further, collect first 100 trajectories of plaquette history
# further, collect mean number of CG iterations in derivative of first monomial
# and mean trajectory time

# res is our return value, it will be a list with the expectation value of the plaquette
# and the rectangle (and their errors) in res[[1]] and the first hundred values of
# the plaquette evolution for each output file in res[[n>1]]

source("readoutput.R")

collect_data <- function(dirs) {
  res <- list()
# create an empty data frame so that res[[1]] exists!
  res[[1]] <- data.frame(c())
  
  # when an output file does not exist we need to keep track of how many files have been "skipped"
  skip <- 0
  for( i in 1:length(dirs) ) {
    norect <- FALSE
    print(dirs[i])
    # create a name for the row from the directory name 
    dirsplit <- strsplit(dirs[i],'/')
    name <- dirsplit[[1]][length(dirsplit[[1]])] 

    # test if the current filename refers to one of the samples without a rectangular part
    for( k in c(1:length(norectsamples)) ) {
      pat <- paste(paste("^",norectsamples[k],sep=""),"_",sep="")
      test <- grep(name,pattern=pat,fixed=FALSE,perl=FALSE)
      if( length(test) > 0 ) {
        norect <- TRUE
      }
    }

    # the old output.data format had no trajectory counter
    # the new format therefore shifted all output to the right (format == 1)
    if( length( grep(dirs[i],pattern="5.1.6",fixed=TRUE) ) > 0 ) {
      format <- 0
    } else {
      format <- 1
    }

    # some samples do not have a CG in the HMC so we don't need to plot this
    if( length( grep(name,pattern="nosplit_ndclover") ) > 0 || 
        length( grep(name,pattern="nosplit_nocsw_ndclover") ) > 0 ) {
      nocg <- TRUE
    } else {
      nocg <- FALSE
    }

    ofile <- paste(dirs[i],"/output.data",sep="")
    parafile <- paste(dirs[i],"/output.para",sep="")

    if( file.exists(ofile) && ( as.integer( strsplit( system(paste("wc -l",ofile),intern=TRUE), " ")[[1]][[1]] ) > min+minlength ) ) {
      
      # this should be replaced by an actual calculation
      if( file.exists(parafile) ) {
        grepcommand <- paste("grep 'Acceptance rate'",parafile)
        grepcommand <- paste(grepcommand,"| awk '{print $5}'")

        # when the job is still running the output.para file will not contain the acceptance rate yet
        # so we check the length explicitly
        TAR <- as.numeric(system(intern=TRUE, grepcommand))
        if( length(TAR) >= 1 ) {
          tar <- round(mean(TAR),digits=2)
        } else {
          tar <- NA
        }
      } else {
        tar <- NA
      }

      outdat <- readoutput(filename=ofile,format=format,norect=norect,nocg=nocg)
      
      # add histories to return value
      res[[2+4*(i-1)-4*skip]] <- outdat$plaq.hist
      res[[2+4*(i-1)+1-4*skip]] <- outdat$rect.hist
      res[[2+4*(i-1)+2-4*skip]] <- outdat$cgitnum.hist
      res[[2+4*(i-1)+3-4*skip]] <- outdat$trajtime.hist 

      # the rectangle and CG iterations may be stored as "NA" in outdat
      # because the construction of the return value requires accessing elements of "uwerr",
      # we are going to define some dummy variables here and assign the actual measurements
      # if they exist

      rec <- NA
      drec <- NA
      recmed <- NA
      rectauint <- NA
      recdtauint <- NA

      cgitnum <- NA
      dcgitnum <- NA
      cgitnummed <- NA

      if(!norect) {
        rec <- outdat$rect.uwerr$value
        drec <- outdat$rect.uwerr$dvalue
        recmed <- median(outdat$rect.uwerr$data)
        rectauint <- outdat$rect.uwerr$tauint
        recdtauint <- outdat$rect.uwerr$dtauint
      }

      if(!nocg) {
        cgitnum <- outdat$cgitnum.uwerr$value
        dcgitnum <- outdat$cgitnum.uwerr$dvalue
        cgitnummed <- median(outdat$cgitnum.uwerr$data)
      }

      tresl <- data.frame(row.names=name,
        list(ar=tar,
          trajtime=outdat$trajtime.uwerr$value, dtrajtime=outdat$trajtime.uwerr$dvalue, trajtimemed=median(outdat$trajtime.uwerr$data),
          cgitnum=cgitnum,dcgitnum=dcgitnum,cgitnummed=cgitnummed,
          plaq=outdat$plaq.uwerr$value,dplaq=outdat$plaq.uwerr$dvalue,plaqmed=median(outdat$plaq.uwerr$data),
          plaqtauint=outdat$plaq.uwerr$tauint,plaqdtauint=outdat$plaq.uwerr$dtauint,
          rec=rec,drec=drec,recmed=recmed,rectauint=rectauint,recdtauint=recdtauint))

      # save the current state of our return value
      tres <- res[[1]]
      # collate our new data frame with the current state
      res[[1]] <- rbind(tres,tresl)
    } else if ( ! file.exists(ofile) )  {
      print(paste(ofile,"does not exist! Skipping."))
      skip <- skip+1
    } else if ( as.integer( strsplit( system(paste("wc -l",ofile),intern=TRUE), " ")[[1]][[1]] ) < min+minlength ) {
      print(paste(paste(ofile,"has less than"),paste(min+minlength,"lines! Skipping.")))
      skip <- skip+1
    }
  }
  return(res)
}
