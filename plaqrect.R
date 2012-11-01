# given output.data of a high statistics run, plaqrect produces expectation
# values and errors based on uwerrprimary
# norect (boolean) signifies that this output.data does not have a column for the
# rectangle
# format=1,0 specifies whether this is an output file in the new format (with n+1 colums)
# or in the old format without an iteration counter (and hence n columns)


library(hadron)

plaqrect <- function(filename,norect,format,ndclover) {

  data <- read.table(filename)
 
  pcol <- format+1
  reccol <- length(data)

  # temporary workaround for broken output.data due to CLOVERNDTRLOG
  if( ndclover ) {
    cgitnumcol <- format+5+2
  } else {
    cgitnumcol <- format+5
  }

  if(norect) {
    trajtimecol <- length(data)
  } else {
    trajtimecol <- length(data)-1
  }

  if( length(data[,pcol]) < 11000 ) {
    return(NA)
  }

  min <- 5000
  max <- length(data[,pcol])

  trajtimet <- mean(data[min:max,trajtimecol])
  dtrajtimet <- sd(data[min:max,trajtimecol])
  trajtimemedt <- median(data[min:max,trajtimecol])

  cgitnumt <- mean(data[min:max,cgitnumcol])
  dcgitnumt <- sd(data[min:max,cgitnumcol])
  cgitnummedt <- median(data[min:max,cgitnumcol])

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
