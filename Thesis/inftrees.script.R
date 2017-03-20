################INFTREES############
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
  d2$p <- as.factor(d2$p)
  for(i in 1:length(levels(d2$p))){
    #set.seed(1)
    d2[d2$p == i, x] <- sample(d2[d2$p== i,x])
  }
  return(d2[,x])
}

VIdevBagged<- function(b, x){
  return(as.numeric(as.data.frame(b$mtrees[[1]]$btree$frame)[,c(1,4)] %>% filter(var == x) %>% dplyr::summarise(sum = sum(dev))))
}

VIdev<- function(t2, x){
  return(as.data.frame(t2$frame)[,c(1,3)] %>% 
           filter(var == x) %>% dplyr::summarise(sum = sum(dev)))
}

variable.importance.dev <- function(t1, x, d){
  if(length(t1$frame$splits[,1]) == 1) {
    return(VIdev(t1, x))
  } else {
    interval <- intervalReturn(t1$frame$splits[,1])
    d$p <- partitionCreator(interval, x, d)
    d[,x] <- condPermuter(d, x)
    t2 <- tree(y ~ ., data = d, y = TRUE)
    v1 <- VIdev(t2, x)
    #rm(t1, t2)
    return(v1)
  }
}

cond.varImp <- function(xs, data) {
  x <- names(xs)
  vi <- rep(0,length(x))
  vo <- rep(0, length(x))
  
  t0 <- tree(y~., data)
  
  set.seed(1)
  for(i in 1:length(x)){
    vo[i] <- VIdev(t0, x[i])
  }
  
  set.seed(1)
  for(i in 1:length(x)){
    form <- as.formula(paste(x[i], "~."))
    ti <- tree(form, xs)
    vi[i] <- variable.importance.dev(ti, x[i], data)
  }
  
  v <- data.frame("variable" = x, "permuted_variable_importance" = as.numeric(vi), "base_variable_importance"= as.numeric(vo))
  return(v)
}

