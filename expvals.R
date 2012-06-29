# given a vector of directories containing output.data, call plaqrect
# to get expectation value and uwerr of plaquette and rectangular
# plaquette, further, collect first 100 trajectories of plaquette history
# from plaqrect

# res is our return value, it will be a list with the expectation value of the plaquette
# and the rectangle (and their errors) in res[[1]] and the first hundred values of
# the plaquette evolution for each output file in res[[n>1]]

source("plaqrect.R")

expvals <- function(dirs) {
# these are samples without a rectangular guage action, we keep track of this so
# we don't produce excess data
  norectsamples <- c("hmc0","hmc1","hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
  
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
      test <- grep(name,pattern=norectsamples[k],fixed=TRUE)
      if( length(test) > 0 ) {
        norect <- TRUE
      }
    }

    ofile <- paste(dirs[i],"/output.data",sep="") 
    if(file.exists(ofile)) {
      tres <- plaqrect(ofile,norect)
  
      res[[2*i+1-2*skip]] <- tres$plhist
                
      trec <- NA
      tdrec <- NA
      trecmed <- NA
      trectauint <- NA
      trecdtauint <- NA
      trechist <- NA
      
      # if we have a rectangle part, read out values into temporary variables
      if( !norect ) {
        trec <- tres$rec$value
        tdrec <- tres$rec$dvalue
        trecmed <- median(tres$rec$data[10000:length(tres$rec$data)])
        trectauint <- tres$rec$tauint
        trecdtauint <- tres$rec$dtauint
        trechist <- tres$rechist
      }

      # add plaquette history and (possibly NA) rectangle history to return value
      res[[2*i-2*skip]] <- tres$plhist
      res[[2*i+1-2*skip]] <- trechist

      # data frame will contain NA if rectangle part is missing, values otherwise
      tresl <- data.frame(row.names=name,
        list(plaq=tres$pl$value,dplaq=tres$pl$dvalue,
          plaqmed=median(tres$pl$data[10000:length(tres$pl$data)]),
          plaqtauint=tres$pl$tauint,plaqdtauint=tres$pl$dtauint,
          rec=trec,drec=tdrec,recmed=trecmed,
          rectauint=trectauint,recdtauint=trecdtauint))

      # save the current state of our return value
      tres <- res[[1]]
      # collate our new data frame with the current state
      res[[1]] <- rbind(tres,tresl)
    } else {
      print(paste(ofile,"does not exist! Skipping."))
      skip <- skip+1
    }
  }

  return(res)
}
