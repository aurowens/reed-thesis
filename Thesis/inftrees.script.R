################INFTREES############

partition.grabber <- function(vars, sp){
  p <- rep("", 1:dim(splits)[1])
  for (i in 1:dim(sp)[1]) {p[i] <- sp[[i]]}
  p[p == ""] <- NA
  p <- as.double(p)
  return(data.frame("cutleft" = p, "vars" = vars))
}

##start w/tree?
inftree <- function(t, x, y) {
  rss <- predict(t, x, y)
  vi <- rep(0, ncol(x))
  for(i in 1:(ncol(d)-1)) {
    x1 <- x
    ti <- tree(y = x[,i], x= x[,-1], data = x)
    p <- partition.grabber(as.character(ti$frame$var), ti$frame$splits) #p is df with vars and their splits
    for( j in nrow(p)) {
      x1[p[i],i] <- sample(x[p[i],i])
    }
    vi[i] <- rss - predict(ti, x1, y)
  }
  return(vi)
}

data("mtcars")
terr <- tree(mpg~., data = mtcars)
