
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
  op.partition <- optimise(rss.tree, interval = range(x), a = y, maximum = FALSE, tol = .01)
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
  maxxr <- max.cor(yy = y,sxs= xssrd)
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


noder <- function(y,xs,mtry, df, min){
  
  ##what i want this function to do:
  ###given y,xs (like an interval of them):
  ####1. find the split between them 
  print("noder:starting node")
  frame <- node1(y,xs, mtry, min)
  df <-  rbind(df,frame)
  ####2. call itself on each side of the split
  split <- as.numeric(frame[5])
  spliton <- frame[1]
  
  if(frame[1] == "<leaf>"){
    print(df)
    return(df)
  }
  print(df)
  yhatl <- y[xs[,spliton[[1]]] < split]
  xsl <- xs[xs[,spliton[[1]]] < split,]
  yhatr <- y[xs[,spliton[[1]]] >= split]
  xsr <- xs[xs[,spliton[[1]]] >= split, ]
  
  if(length(yhatl) < min & length(yhatr) < min){
    print(df)
    return(df)
  } else {
    #y = yhatr
    #xs = xsr
    noder(y = yhatr,xs = xsr,mtry, df, min)
    noder(y = yhatl,xs = xsl,mtry, df, min)
   
  }
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  xss <- xs[,sample(ncol(xs),mtry)]
  y <- y[bootsample]
  min <- round(nrow(xs)/5)
  
####First split   
  first.split <- node1(y,xss,mtry,min)
  if(first.split[1] == "<leaf>"){
    return(first.split)
  }
  
  
  yhatl <- y[xss[,first.split[[1]]] < first.split[[5]]]
  xsl <- xs[xss[,first.split[[1]]] < first.split[[5]],]
  yhatr <- y[xss[,first.split[[1]]] >= first.split[[5]]]
  xsr <- xs[xss[,first.split[[1]]] >= first.split[[5]],]
  
  frame <- first.split
  
  frameld <- noder(yhatl,xsl,mtry, tree = frame, min)
  framerd <- noder(yhatr, xsr, mtry, tree = frame, min)
  tree <- rbind(frame, frameld, framerd)
  return(tree)
}

t <- tree.rf(y = d1$y, xs = d1[,1:12], mtry = 5)
t <- tree.rf(y = iris$Sepal.Length, xs = iris[,2:4], mtry = 2)
y = iris$Sepal.Length
xs = iris[,2:4]
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

noo <- function(c){
  return(moo(c))
}

foo <- function(a) {
  bee <<- 0
  moo <- function(b){
    b <- b +1
    bee <<- 9
    return(b)
  }
  a <- noo(a)
  return(c(a,bee))
}
