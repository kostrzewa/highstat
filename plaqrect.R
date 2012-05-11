# given output.data of a high statistics run, plaqrect produces expectation
# values and errors based on uwerrprimary
# norect signifies that this output.data does not have a column for the
# rectangle

library(hadron)

plaqrect <- function(filename,norect) {
  data <- read.table(filename)
  min <- 5000
  if( length(data[,1]) > 2000 ) {
    min <- 2000
  }

  plaq <- data[min:length(data[,1]),1]
  plaqres <- uwerrprimary(plaq)
  #plaqsum <- summary(plaqres)

  if(!norect) {
    rect <- data[min:length(data[,length(data)]),length(data)]
    rectres <- uwerrprimary(rect)
  } else {
    rect <- NA
    rectres <- NA
  }
  #rectsum <- summary(rectres)

  return(list(pl=plaqres,rec=rectres))
}
