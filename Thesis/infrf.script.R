
############################################################################################
#####################################INFFORESTS#############################################
############################################################################################
############################################################################################

rss.tree <- function(b,a) { 
  a <- rep(mean(a),length(a))
  sum((a-b)^2)
}

max.cor <- function(y,xss){
  max <- list()
  cors <- cbind(rep(0,ncol(xss)),seq(1,ncol(xss)))
  for(i in 1:ncol(xss)){
    cors[i,1] <- cor(y,xss[,i])
    }
  max[[length(max)+1]] <- xss[,cors[cors[,1] == max(cors[,1]),2]]
  max[[length(max)+1]] <- names(xss)[cors[cors[,1] == max(cors[,1]),2]]
  return(max)
}

split.rf <- function(y,x, min) {
  sp <- list()
  op.partition <- optimise(rss.tree, interval = range(x), a = y, maximum = FALSE)
  if(length(x[x < op.partition[[1]]]) < min | length(x[x >= op.partition[[1]]]) < min) {
    return(NULL)
    } else {
      sp[[length(sp)+1]] <- op.partition[[1]] #partition on x
      sp[[length(sp)+1]]<-  mean(y[x< op.partition[[1]]]) #ypred left daughter
      sp[[length(sp)+1]]<- mean(y[x >= op.partition[[1]]]) #ypred right daughter
      sp[[length(sp)+1]]<- x[x < op.partition[[1]]] #xsleftdaughter
      sp[[length(sp)+1]]<- x[x >= op.partition[[1]]] #xsrightdaughter
      sp[[length(sp)+1]] <- rss.tree(x[x < op.partition[[1]]], y[x < op.partition[[1]]]) #ldrss
      return(sp)
    }
}



##yrd <- y[max >= spl]
##do first   xs <- xs[max >= spl,]
node1 <- function(y, xs, mtry, spl){
  xssrd <- xs[,sample(ncol(xs),mtry)]
  maxxr <- max.cor(y,xs)
  sprd <- split.rf(y, maxxr[[1]], min = round(nrow(xs)/20))
  frame <- c()
  if(is.null(sprd)) {
    frame <- c("<leaf>",nrow(xssrd), mean(y), rss.tree(maxxr[[1]], y), 0)
  } else {
    frame <- c(maxxr[[2]],
              nrow(xssrd),
              mean(y),
              rss.tree(xssrd, y),
              sprd[[1]])
  }
  return(frame)
}


noder <- function(y,xs,sp,max,mtry, tree){
  
  l1 <- node1(y[max[[1]] < sp],xs[max[[1]] < sp,], mtry, sp)
  r1 <- node1(y[max[[1]] >= sp],xs[max[[1]] >= sp,], mtry, sp)
  
  var.l1 <- l1[1]
  var.r1 <- r1[1]
  
  tree <-  rbind(tree,l1)
  tree <-  rbind(tree,r1)
  
  if(var.r1 == "<leaf>" & var.l1 == "<leaf>"){
    return(tree)
  } else if(var.r1 == "<leaf>" & var.l1 != "<leaf>") {
    print("right daughter null, following left")
    noder(y,xs,as.numeric(l1[5]),list(xs[,l1[1]], l1[1]),mtry, tree)
  } else if (var.r1 != "<leaf>" & var.l1 == "<leaf>"){
    print("left daughter null, following right")
    noder(y,xs,as.numeric(r1[5]),list(xs[,r1[1]], r1[1]),mtry, tree)
  } else {
    print("neither null, following all")
    noder(y,xs,as.numeric(r1[5]),list(xs[,r1[1]], r1[1]),mtry, tree)
    noder(y,xs,as.numeric(l1[5]),list(xs[,l1[1]], l1[1]),mtry, tree)
  }
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  y <- y[bootsample]
  xss <- xs[,sample(ncol(xs),mtry)]

  maxx <- max.cor(y, xss)
  max.name <- maxx[[2]]
  max <- maxx[[1]]
  spli <- split.rf(y, max, min = 10)
  
  frame <- node1(y,xs, sp = spli[[1]], mtry)
  frame <- noder(y,xs,sp = spli[[1]], list(xs[,frame[1]], frame[1]), mtry, tree = frame)
  tree <- frame
  return(tree)
}

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