
############################################################################################
#####################################INFFORESTS#############################################
############################################################################################
############################################################################################

rss.tree <- function(b,a) { 
  a <- rep(mean(a),length(a))
  sum((a-b)^2)
}

max.cor <- function(yy,sxs){
  print("max cor")
  max <- list()
  cors <- cbind(rep(0,ncol(sxs)),seq(1,ncol(sxs)))
  for(i in 1:ncol(sxs)){
    cors[i,1] <- cor(yy,sxs[,i])
  }
  if (sum(cors[,1] == max(cors[,1], na.rm = TRUE), na.rm = TRUE) > 1) {
    max1 <- sample(as.numeric(cors[,1] == max(cors[,1], na.rm = TRUE)), 1)  
  } else {
    max1 <- cors[cors[,1] == max(cors[,1], na.rm = TRUE),2]
  }
  
  max[[length(max)+1]] <- sxs[,max1]
  max[[length(max)+1]] <- names(sxs)[max1]
  return(max)
  
}

split.rf <- function(y,x, min) {
  print("split.rf")
  sp <- list()
  if(length(x) <= min) {
    print("split.rf, leaf")
    return(NULL)
  } else {
  op.partition <- optimise(rss.tree, interval = range(x[c(min:(length(x)-min))]), a = y, maximum = FALSE, tol = .01)
      sp[[length(sp)+1]] <- op.partition[[1]] #partition on x
      sp[[length(sp)+1]]<-  mean(y[x< op.partition[[1]]]) #ypred left daughter
      sp[[length(sp)+1]]<- mean(y[x >= op.partition[[1]]]) #ypred right daughter
      sp[[length(sp)+1]]<- x[x < op.partition[[1]]] #xsleftdaughter
      sp[[length(sp)+1]]<- x[x >= op.partition[[1]]] #xsrightdaughter
      sp[[length(sp)+1]] <- rss.tree(x[x < op.partition[[1]]], y[x < op.partition[[1]]]) #ldrss
      print("split.rf done")
      return(sp)
  }
}

##yrd <- y[max >= spl]
##do first   xs <- xs[max >= spl,]
node1 <- function(y, xs, mtry, min){
  print("splitting on a new interval")
  xssrd <- xs[,sample(ncol(xs),mtry)]
  frame <- data.frame("", 0,0,0,0)
  frame[,1] <- as.character(frame[,1])
  if(length(y) < min) {
    frame[1,] <- c("<leaf>",nrow(xssrd), mean(y), 0, 0)
    print("node1 done, leaf")
    return(frame)
  } else {
  maxxr <- max.cor(yy=y,sxs=xs)
  sprd <- split.rf(y, maxxr[[1]], min)
  if(is.null(sprd)) {
    frame[1,] <- c("<leaf>",nrow(xssrd), mean(y), 0, 0)
    print("node1 done, leaf")
    return(frame)
  }
  frame[1,] <- c(maxxr[[2]],
              nrow(xssrd),
              mean(y),
              rss.tree(xssrd, y),
              sprd[[1]])
   # frame[2:4] <- as.numeric(frame[2:4])
  }
  print("node1 done")
  return(frame)
}


noder <- function(y,xs,mtry, tree, min){
  
  ##what i want this function to do:
  ###given y,xs (like an interval of them):
  ####1. find the split between them 
  print("noder:starting node")
  frame <- node1(y,xs, mtry, min)

  tree <-  rbind(tree,frame)
  ####2. call itself on each side of the split
  split <- as.numeric(frame[5])
  spliton <- frame[1]
  
  if(frame[1] == "<leaf>"){
    return(tree)
  }
 
  yhatl <- y[xs[,spliton[[1]]] < split]
  xsl <- xs[xs[,spliton[[1]]] < split,]
  yhatr <- y[xs[,spliton[[1]]] >= split]
  xsr <- xs[xs[,spliton[[1]]] >= split, ]
  
  if(length(yhatl) < min & length(yhatr) < min){
    return(tree)
  } else {
    noder(y = yhatr,xs = xsr,mtry, tree, min)
    noder(y = yhatl,xs = xsl,mtry, tree, min)
  }
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  y <- y[bootsample]
  
  min <- round(nrow(xs)/5)
  frame <- noder(y,xs,mtry, tree = list(), min)
  return(tree)
}

t <- tree.rf(y = d1$y, xs = d1[,1:12], mtry = 5)
y = d1$y
xs = d1[,1:12]
############################################################################################
#####################################GRAVEYARD##############################################
############################################################################################
############################################################################################

rf <- function(form, d, mtry, ntree) { 
  #imputs: formula to be tested and the dset to test on. outputs: tribble,
  #contains: total oob error, every tree
  rforest <- list() #define empty list for the rf. length = ntree + 2
  #where each tree is a list of 6, and there's a vector of sample indexes for each tree 
  B <- list()
  xsi <- c()
  di <- data.frame()
  sampleset <- c()
  
  yn <- gsub( "~.*", "", form)
  yn <- trimws(yn)
  yi <- as.numeric(names(d) == yn)
  y <- d[,yi]
  xs <- d[,-(yi)]
  ratio <- round(.75*nrow(d))
  
  while (length(rforest) < ntree) {
    sampleset <- sample(nrow(d), ratio, replace = TRUE)
    xsi <- sample((ncol(d)-1),mtry, replace = FALSE)
    B[[length(B)+1]] <- sampleset
    di <- cbind(xs[sampleset,xsi], y[sampleset])
    names(di)[ncol(di)] <- yn 
    rforest[[length(rforest)+1]] <- tree(as.formula(form), d = di)
    }
  rforest[[length(rforest)+1]] <- B
  return(rforest)
}