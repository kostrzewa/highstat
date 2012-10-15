# given output.data of a high statistics run, plaqrect produces expectation
# values and errors based on uwerrprimary
# norect signifies that this output.data does not have a column for the
# rectangle

library(hadron)

plaqrect <- function(filename,norect) {
  data <- read.table(filename)
  # old format without iteration counter -> data[,1] is the plaquette!
  if( length(data[,2]) < 11000 ) {
    return(NA)
  }

  min <- 10000

  plaq <- data[min:length(data[,2]),2]
  plaqres <- uwerrprimary(plaq)
  plaqhist <- data[,2]

  if(!norect) {
    rect <- data[min:length(data[,length(data)]),length(data)]
    rectres <- uwerrprimary(rect)
    recthist <- data[,length(data)]
  } else {
    rect <- NA
    rectres <- NA
    recthist <- NA
  }

  return(list(pl=plaqres,rec=rectres,plhist=plaqhist,rechist=recthist))
}
