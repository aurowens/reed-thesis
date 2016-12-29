library(tidyverse)
library(randomForest)
library(tree)

#################################################################################################
######################################INF TREES##################################################
#################################################################################################


intervalReturn <- function(split.cutleft){ 
  intervals <- data.frame(split.cutleft)
  intervals <- filter(intervals, split.cutleft != "")
  intervals <- dplyr::select(intervals, split.cutleft)
  intervals$split.cutleft <- gsub("<", "", intervals$split.cutleft)
  intervals$split.cutleft <- as.numeric(intervals$split.cutleft)
  return(arrange(intervals, (split.cutleft)))
}

partitionCreator <- function(i = interval, x, data) {
  p <- c(1:(1+nrow(i)))
  bins <- data.frame (bin = rep(1, nrow(data)), v = data[,x])
  for(j in 1:nrow(i)){
    bins[bins$v > i[j,1],1] <- p[j+ 1]
  }
  return(bins$bin)
}

condPermuter <- function(d2, x) {
  for(i in 1:length(unique(d2$p))){
    set.seed(1)
    d2[d2$p == i, x] <- sample(d2[d2$p== i,x])
  }
  return(d2[,x])
}

permuter <- function(xs, xi) {
  fo1 <- as.formula(paste(xi, "~."))
  ti <- tree(fo1, xs)
  if(length(ti$frame$splits[,1]) == 1) {
    stop()
  } else {
    interval <- intervalReturn(ti$frame$splits[,1])
    xs$p <- partitionCreator(interval, xi, xs)
    return(condPermuter(xs, xi))
  }
}


cp.importance <- function(fo, df){
  fo <- as.formula(fo) 
  y <- as.character(fo)[2]
  xs <-df[ , !(names(df) %in% y)]
  names.xs <- names(xs)
  
  rf0 <- randomForest(fo, data = df, ntree=500, mtry = round(sqrt(ncol(df) - 1)))
  
  vi <- data.frame("variable" = names.xs, 
                   "cp.importance" = rep(0,length(names.xs)), 
                   "base.importance" = rep(0, length(names.xs)))
  
  vi[,3] <- t(as.data.frame.list(importance(rf0)))[,1] 
  
  for(i in 1:nrow(vi)) {
    pdf <- df
    pdf[,names.xs[i]] <- permuter(xs, names.xs[i])
    rfi <-  randomForest(fo, pdf, ntree=500, mtry = round(sqrt(ncol(df) - 1)))
    vi[i, 2] <- vi[i,3] - t(as.data.frame.list(importance(rfi)))[i]
  }
  return(vi)
}

