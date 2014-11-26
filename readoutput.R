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

readoutput <- function(filename,format,norect,nocg) {

  data <- read.table(filename)
 
  plaqcol <- format+1
  rectcol <- ncol(data)

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

  if( length(data[,plaqcol]) < min+minlength ) {
    return(NA)
  }

  # the limit and trajs parameters are defined in highstat.R
  if( limit && length(data[,plaqcol]) > min+trajs ) {
    max <- min+trajs
  } else {
    max <- length(data[,plaqcol])
  }

  plaq.uwerr <- uwerrprimary(data[min:max,plaqcol])
  plaq.hist <- data[,plaqcol]

  trajtime.uwerr <- uwerrprimary(data[min:max,trajtimecol] )
  trajtime.hist <- data[1:max,trajtimecol]

  if( nocg ) {
    cgitnum.uwerr <- NA
    cgitnum.hist <- rep(NA,times=max)
  } else {
    cgitnum.uwerr <- uwerrprimary(data[min:max,cgitnumcol])
    cgitnum.hist <- data[1:max,cgitnumcol]
  }

  if(norect) {
    rect.uwerr <- NA
    rect.hist <- rep(NA,times=max)
  } else {
    rect.uwerr <- uwerrprimary(data[min:max,rectcol])
    rect.hist <- data[1:max,rectcol]
  }

  return(list(plaq.hist=plaq.hist,plaq.uwerr=plaq.uwerr,
              rect.hist=rect.hist,rect.uwerr=rect.uwerr,
              trajtime.hist=trajtime.hist,trajtime.uwerr=trajtime.uwerr,
              cgitnum.hist=cgitnum.hist,cgitnum.uwerr=cgitnum.uwerr))
}
