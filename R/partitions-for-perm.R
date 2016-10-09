#' Converts tree breaks to a more useful vector of breaks.
#' 
#' @param x The splits out of an object created by the `tree` function: `tree$frame$split.cutleft
#' @return A sorted 1D dataframe of breakpoints.
intervalReturn <- function(x = split.cutleft){ 
intervals <- data.frame(split.cutleft)
intervals <- filter(intervals, split.cutleft != "")
intervals <- select(intervals, split.cutleft)
intervals$split.cutleft <- gsub("<", "", intervals$split.cutleft)
intervals$split.cutleft <- as.numeric(intervals$split.cutleft)
return(arrange(intervals, (split.cutleft)))
}

#' Uses the output from `intervalReturn` and the variable we are permuting with respect to assign bins to each observation
#' 
#' @param i output from `intervalReturn`
#' @param x the variable that we are permuting with respect to
#' @return A 1D dataframe of the binned partitions.
partitionCreator <- function(i = interval, x = data) {
  p <- c(1:(1+nrow(i)))
  bins <- data.frame (bin = rep(1, length(x)), x = x)
  for(j in 1:nrow(i)){
    bins[bins$x > i[j,1],1] <- p[j+ 1]
  }
  return(bins$bin)
}