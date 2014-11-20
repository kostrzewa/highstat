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
library(parallel)

source("plotfunc.R")

# we use global variables to enable parallelization

#samples <- c(
  #"mpihmc1","mpihmc2","mpihmc3","mpihmc4","mpihmc5","mpihmc6","mpihmc7","mpihmc8","mpihmc211",
  #"hmc0","hmc_repro_0","hmc1","hmc_repro_1","hmc2","hmc_repro_2","hmc3",
  #"hmc_ndclover","hmc_nosplit_ndclover","hmc_nocsw_ndclover","hmc_nosplit_nocsw_ndclover",
  #"hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")#,
  #"hmc_check_ndclover_tmcloverdet", "hmc_check_ndclover_nocsw_tmcloverdet")

#samples <- c("hmc_ndrat")
#samples <- c("hmc0")

samples <- c("2f_L8T8_ndcloverrat","8f_L8T8_ndcloverrat")

execs <- c("openmp","openmp_hs","serial","serial_hs",
  "1D_hybrid_hs_2","2D_hybrid_hs_4","4D_hybrid_hs_16",
  "BGQ_4D_SPI_hybrid_hs_32","BGQ_4D_SPI_hybrid_hs_128",
  "4D_MPI_hs_16","4D_MPI_hs_32","4D_MPI_hs_64","4D_MPI_hs_128","4D_MPI_hs_256",
  "3D_MPI_hs_8","3D_MPI_hs_16","3D_MPI_hs_24","3D_MPI_hs_32","3D_MPI_hs_64",
  "2D_MPI_hs_16","1D_MPI_hs_16",
  "5.1.6_3D_MPI_hs_8","5.1.6_serial")

# for debugging purposes, we can shuffle the vector
# samples <- sample(samples,size=length(samples))

# these samples do not contain a rectangular gauge part
norectsamples <- c("mpihmc2","mpihmc4","mpihmc6","mpihmc8",
  "hmc0","hmc_repro_0","hmc1","hmc_repro_1",
  "hmc_cloverdet","hmc_tmcloverdet","hmc_tmcloverdetratio")
 
# these are reference values for the plaquette and rectangle expectation value 
reference <- read.table("reference.dat", fill=FALSE,sep=" ", header=TRUE, row.names=1)

min <- 150
minlength <- 500

# the limit parameter specifies whether the trajectory series should be cut
# the trajs parameter specifies where it will be cut ( min + trajs )
limit <- FALSE
trajs <- 6000-min

topdir
topdirname 
subdir

highstat <- function(tdir,name) {
  topdir <<- tdir

  # extract the "name" of the top directory
  topdirsplit <- strsplit(topdir,'/')
  topdirname <<- topdirsplit[[1]][ length(topdirsplit[[1]]) ]

 
  # directory where plots will be created
  if(!missing(name)) {
    subdir <<- paste("plots",name,sep="/")
  } else {
    subdir <<- paste("plots",topdirname,sep="/")
  }
  
  if( ! file.exists( subdir ) ) {
    dir.create(subdir)
  }
 
  # process samples in parallel, spawning 8 processes
  # change to lapply in case of errors! 
  # also makes debugging easier in case of errors not related to multicore
  timelist <- mclapply( samples, FUN=plotfunc , mc.cores = 8 ,mc.preschedule=FALSE)
  #timelist <- lapply( samples, FUN=plotfunc )

  # collect timing information in a table
  for (i in seq(1,length(timelist))) {
    if( i==1 ) {
      timetable <- timelist[[i]]
    } else {
      timetable <- rbind(timetable,timelist[[i]])
    }
  }

  # transpose to get correct format for genjobscripts.sh
  timetable <- t(timetable)
  write.table(timetable,file=paste(subdir,"runtimes.csv",sep="/"),quote=FALSE,row.names=TRUE,col.names=NA) 
}
 
