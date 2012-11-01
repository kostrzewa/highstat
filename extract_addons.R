extract_addons <- function(x,split) {
  addon <- c()
  for ( j in 1:length(x) ) {
    vaddon <- strsplit(x[j],split=paste(split,"_",sep=""),fixed=TRUE)
    addon <- append(addon,vaddon[[1]][length(vaddon[[1]])])
  }
  return(addon) 
}
