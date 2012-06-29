# given output.data of a high statistics run, plaqrect produces expectation
# values and errors based on uwerrprimary
# norect signifies that this output.data does not have a column for the
# rectangle

library(hadron)

plaqrect <- function(filename,norect) {
  data <- read.table(filename)
  min <- 10000

  plaq <- data[min:length(data[,1]),1]
  plaqres <- uwerrprimary(plaq)
  plaqhist <- data[,1]

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
