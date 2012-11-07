# highstat is a function which uses Carsten Urbach's hadron library
# to compute expectation values of the plaquette and the rectangle
# for the set of sample files of the tmLQCD software package
# under the hood it also computes the optimal autocorrelation time
# accoring to U. Wollf's method

# it produces two to three files for each sample file,
# the first file contains a plot for for the plaquette expectation value, one for
# the rectangle expectation value (unless it contains no rectangle gauge contribution),
# one plot for the cumulative CG iterations of the derivative of the first monomial,
# one plot for the mean time per trajactory 
# and one plot each for the first 100 trajectories of the plaquette and rectangle
# histories (as available)
# the second and third files contain complete histories of the plaquette and rectangle part

# the argument "tdir" must be the path to a directory which contains
# a set of subdirectories
# these subdirectories in turn relate to different runs of the different
# sample input files
# each subdirectory must contain "output.data" for the given run

# the names of the subdirectories should begin with the strings in the 
# "samples" vector below, separated from the rest of the filename
# by an underscore ('_')

# the rest of the filename is referred to as an "addon" in several places
# in the code, this "addon" should differentiate the various runs of the
# same sample file. For instance, this could be different parallelizations.
#    e.g. mpi, mpi_hs, serial, hybrid_hs etc..
# or maybe different versions of the executable used to generate the data

# "norectsamples" is a subset of samples which do not contain a rectangle
# contribution in their gauge action

# example directory structure
# tdir
# |
# -----> hmc0_old_serial
# |  |-> hmc0_new_serial 
# |  |-> hmc0_old_mpi
# |  |-> hmc0_new_mpi
# |
# -----> hmc1_[...]
# [...] 

# because of the way that the subfunctions are called this should be run
# in the directory containing the R files by first calling
#    source("highstat.R")

library(hadron)
library(multicore)

source("plotfunc.R")

# we use global variables to enable parallelization


# the samples vector must be in the same order as in genjobscripts for the reference values produced from here to be usable
samples <- c("hmc0","hmc1","hmc2","hmc3","hmc_ndclover","hmc_nosplit_ndclover","hmc_nocsw_ndclover","hmc_nosplit_nocsw_ndclover","hmc_cloverdet","hmc_tmcloverdet","hmc_check_ndclover_tmcloverdet", "hmc_check_ndclover_nocsw_tmcloverdet","hmc_tmcloverdetratio")

# these samples do not contain a rectangular gauge part
norectsamples <- c("hmc0","hmc1","hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
 
# these are reference values for the plaquette and rectangle expectation value 
reference <- read.table("reference.dat", fill=TRUE)

min <- 100
minlength <- 500

topdir
topdirname 
subdir

highstat <- function(tdir) {
  topdir <<- tdir

  # extract the "name" of the top directory
  topdirsplit <- strsplit(topdir,'/')
  topdirname <<- topdirsplit[[1]][ length(topdirsplit[[1]]) ]
 
  # directory where plots will be created
  subdir <<- paste("plots",topdirname,sep="/")
  if( ! file.exists( subdir ) ) {
    dir.create(subdir)
  }
 
  # process samples in parallel, spawning 8 processes
  # change to lapply in case of errors! 
  timelist <- mclapply( samples, FUN=plotfunc , mc.cores = 8 )
  
  
  # collect timing information in a table
  timetable <- timelist[[1]]
  for (i in seq(2,length(timelist))) {
    timetable <- rbind(timetable,timelist[[i]])
  }
  timetable <- aperm(timetable, perm = NULL, resize = TRUE, keep.class = TRUE)
  write.table(timetable,file=paste(subdir,"runtimes.csv",sep="/"),quote=FALSE,row.names=TRUE,col.names=NA) 
}
 
