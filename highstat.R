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
library(multicore)

source("expvals.R")
source("plotfunc.R")

# we use global variables to enable parallelization

samples <- c("hmc0","hmc1","hmc2","hmc3","hmc_tmcloverdetratio","hmc_cloverdet","hmc_tmcloverdet","hmc_ndclover","hmc_nocsw_ndclover","hmc_nosplit_ndclover","hmc_nosplit_nocsw_ndclover","hmc_check_ndclover_tmcloverdet","hmc_check_ndclover_nocsw_tmcloverdet")

# these samples do not contain a rectangular gauge part
norectsamples <- c("hmc0","hmc1","hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
 
# these are reference values for the plaquette and rectangle expectation value 
reference <- read.table("reference.dat", fill=TRUE)

topdir
topdirname 
subdir

plotexpvals <- function(x) {
  topdir <<- x

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
  mclapply( samples, FUN=plotfunc , mc.cores = 8 )
}
 
