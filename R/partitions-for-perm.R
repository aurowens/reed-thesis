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