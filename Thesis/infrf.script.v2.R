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

tree.rf <- function(y,xs, mtry){
  tree <- list()
  
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
  node1 <- function(y, xs, mtry, min){
    print("splitting on a new interval")
    xssrd <- xs[,sample(ncol(xs),mtry)]
    frame <- data.frame("", 0,0,0,0)
    frame[,1] <- as.character(frame[,1])
    if(length(y) < min) {
      frame[1,] <- c("<leaf>",nrow(xssrd), mean(y), 0, 0)
      print("node1 done, leaf")
      tree <<- rbind(tree, frame)
      return(frame)
    } else {
      maxxr <- max.cor(yy = y,sxs= xssrd)
      sprd <- split.rf(y, maxxr[[1]], min)
      if(is.null(sprd)) {
        frame[1,] <- c("<leaf>",nrow(xssrd), mean(y), 0, 0)
        print("node1 done, leaf")
        tree <<- rbind(tree, frame)
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
    tree <<- rbind(tree, frame)
    return(frame)
  }

  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  y <- y[bootsample]
  min <- round(nrow(xs)/5)
  
  noder(y, xs, mtry, df = list(), min)
  
  return(tree)
}

t <- tree.rf(y = iris$Sepal.Length, xs = iris[,2:4], mtry = 2)
