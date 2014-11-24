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


# given output.data of a high statistics run, plaqrect produces expectation
# values and errors based on uwerrprimary
# norect (boolean) signifies that this output.data does not have a column for the
# rectangle
# format=1,0 specifies whether this is an output file in the new format (with n+1 colums)
# or in the old format without an iteration counter (and hence n columns)


library(hadron)

readoutput <- function(filename,norect,format,nocg) {

  data <- read.table(filename)
 
  pcol <- format+1
  reccol <- ncol(data)

  # the trajectory time is in the last or second to last column depending on whether we have a rectangle
  # in the gauge action, so that's fine
  # the CG iterations however depend on which monomials are defined... as a potentially safe choice
  # we take the 5th or 6th column from the right

  if(norect) {
    trajtimecol <- ncol(data)
    cgitnumcol <- ncol(data)-4
  } else {
    trajtimecol <- ncol(data)-1
    cgitnumcol <- ncol(data)-5
  }

  if( length(data[,pcol]) < min+minlength ) {
    return(NA)
  }

  # the limit and trajs parameters are defined in highstat.R
  if( limit && length(data[,pcol]) > min+trajs ) {
    max <- min+trajs
  } else {
    max <- length(data[,pcol])
  }

  trajtimet <- mean(data[min:max,trajtimecol])
  dtrajtimet <- sd(data[min:max,trajtimecol])
  trajtimemedt <- median(data[min:max,trajtimecol])

  # if the sample has no CG we skip that
  if( nocg ) {
    cgitnumt <- NA
    dcgitnumt <- NA
    cgitnummedt <- NA
  } else {
    cgitnumt <- mean(data[min:max,cgitnumcol])
    dcgitnumt <- sd(data[min:max,cgitnumcol])
    cgitnummedt <- median(data[min:max,cgitnumcol])
  }

  plaq <- data[min:max,pcol]
  plaqres <- uwerrprimary(plaq)
  plaqhist <- data[,pcol][0:max]

  if(!norect) {
    rect <- data[min:max,reccol]
    rectres <- uwerrprimary(rect)
    recthist <- data[,reccol][0:max]
  } else {
    rect <- NA
    rectres <- NA
    recthist <- NA
  }

  return(list(cgitnum=cgitnumt,dcgitnum=dcgitnumt,cgitnummed=cgitnummedt,
    trajtime=trajtimet,dtrajtime=dtrajtimet,trajtimemed=trajtimemedt,
    pl=plaqres,rec=rectres,plhist=plaqhist,rechist=recthist))
}
