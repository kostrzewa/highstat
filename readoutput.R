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

readoutput <- function(filename,norect,format,brokenndclover,nocg) {

  data <- read.table(filename)
 
  pcol <- format+1
  reccol <- length(data)

  # temporary workaround for broken output.data due to CLOVERNDTRLOG
  if( brokenndclover ) {
    cgitnumcol <- format+5+2
  } else {
    cgitnumcol <- format+5
  }

  if(norect) {
    trajtimecol <- length(data)
  } else {
    trajtimecol <- length(data)-1
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
  plaqhist <- data[,pcol]

  if(!norect) {
    rect <- data[min:max,reccol]
    rectres <- uwerrprimary(rect)
    recthist <- data[,reccol]
  } else {
    rect <- NA
    rectres <- NA
    recthist <- NA
  }

  return(list(cgitnum=cgitnumt,dcgitnum=dcgitnumt,cgitnummed=cgitnummedt,
    trajtime=trajtimet,dtrajtime=dtrajtimet,trajtimemed=trajtimemedt,
    pl=plaqres,rec=rectres,plhist=plaqhist,rechist=recthist))
}
