# given a vector of directories containing output.data, call plaqrect
# to get expectation value and uwerr of plaquette and rectangular
# plaquette

source("~/code/R/highstat/plaqrect.R")

expvals <- function(dirs) {
# these are samples without a rectangular guage action, we keep track of this so
# we don't produce excess data
  norectsamples <- c("hmc0","hmc1","hmc-cloverdet","hmc-tmcloverdet","hmc-tmcloverdetratio")
# create an empty data frame so that res exists!
  res <- data.frame(c())
  for(i in 1:length(dirs) ) {
    norect <- FALSE
    print(dirs[i])
    dirsplit <- strsplit(dirs[i],'/')
    name <- dirsplit[[1]][length(dirsplit[[1]])] 

    # test if the current filename refers to one of the samples without a rectangular part
    for( k in c(1:length(norectsamples)) ) {
      test <- grep(name,pattern=norectsamples[k],fixed=TRUE)
      if( length(test) > 0 ) {
        norect <- TRUE
      }
    }

    ofile <- paste(dirs[i],"/output.data",sep="") 
    if(file.exists(ofile)) {
      tres <- plaqrect(ofile,norect)
      
           
      trec <- NA
      tdrec <- NA
      trectauint <- NA
      trecdtauint <- NA
      
      # if we have a rectangle part, read out values into temporary variables
      if( !norect ) {
        trec <- tres$rec$value
        tdrec <- tres$rec$dvalue
        trectauint <- tres$rec$tauint
        trecdtauint <- tres$rec$dtauint
      }


      # data frame will contain NA if rectangle part is missing, values otherwise
      tresl <- data.frame(row.names=name,
        list(plaq=tres$pl$value,dplaq=tres$pl$dvalue,
          plaqtauint=tres$pl$tauint,plaqdtauint=tres$pl$dtauint,
          rec=trec,drec=tdrec,
          rectauint=trectauint,recdtauint=trecdtauint))

      tres <- res
      res <- rbind(tres,tresl)
    } else {
      print(paste(ofile,"does not exist! Skipping."))
    }
  }
  return(res)
}
