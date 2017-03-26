
############################################################################################
#####################################INFFORESTS#############################################
############################################################################################
############################################################################################

########################################TESTING#############################################

data("iris")
d <- iris[1:4]
y <- iris$Sepal.Length
xs <- iris[2:4]
mtry <- 2
form <- as.formula("Sepal.Length ~.")

tree.rf(y,xs,mtry)
r <- rforest(form, d, mtry, ntree = 2)
###################################BUILDING A TREE#########################################

rss.tree.to.minimize <- function(i,y,x) { 
  yleftpred <- mean(y[x< i])
  yrightpred <- mean(y[x >= i])
  rssl <- sum((y[x< i]-yleftpred)^2)
  rssr <- sum((y[x >= i] - yrightpred)^2)
  return(rssl + rssr)
}

rss.tree <- function(a){
  a1 <- rep(mean(a), length(a))
  sum((a-a1)^2)
}

max.cor <- function(yy,sxs){
  print("max cor")
  print(length(yy))
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
    op.partition <- optimise(rss.tree.to.minimize, interval = range(x), y=y,x=x, maximum = FALSE, tol = .01)
    
    
    if (length(x[x < op.partition[[1]]]) < min | length(x[x >= op.partition[[1]]]) < min){
      print("split.rf, leaf")
      return(NULL)
    }
    sp[[length(sp)+1]] <- op.partition[[1]] #partition on x
    sp[[length(sp)+1]]<-  mean(y[x< op.partition[[1]]]) #ypred left daughter
    sp[[length(sp)+1]]<- mean(y[x >= op.partition[[1]]]) #ypred right daughter
    sp[[length(sp)+1]]<- x[x < op.partition[[1]]] #xsleftdaughter
    sp[[length(sp)+1]]<- x[x >= op.partition[[1]]] #xsrightdaughter
    sp[[length(sp)+1]] <- rss.tree( y[x < op.partition[[1]]]) #ldrss
    print("split.rf done")
    return(sp)
  }
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  tree.frame <- list()
  
  noder <- function(y,xs,mtry, min){
    
    ##what i want this function to do:
    ###given y,xs (like an interval of them):
    ####1. find the split between them 
    print("noder:starting node")
    frame <- node1(y,xs, mtry, min)
    #df <-  rbind(df,frame)
    ####2. call itself on each side of the split
    split <- as.numeric(frame[5])
    spliton <- frame[1]
    
    if(frame[1] == "<leaf>"){
      return(frame)
    } else {
      #print(df)
        yhatl <- y[xs[,spliton[[1]]] < split]
        xsl <- xs[xs[,spliton[[1]]] < split,]
        yhatr <- y[xs[,spliton[[1]]] >= split]
        xsr <- xs[xs[,spliton[[1]]] >= split, ]
    
#    if(length(yhatl) < min & length(yhatr) < min){
#      print(df)
#      return(df)
#    } else {
      #y = yhatr
      #xs = xsr
        noder(y = yhatr,xs = xsr,mtry, min)
        noder(y = yhatl,xs = xsl,mtry, min)
      }
  #  }
  }
  node1 <- function(y, xs, mtry, min){
    print("splitting on a new interval")
    xssrd <- xs[,sample(ncol(xs),mtry)]
    frame <- data.frame("", 0,0,0,0)
    frame[,1] <- as.character(frame[,1])
    
    if(length(y) < min) {
      frame[1,] <- c("<leaf>",nrow(xssrd), rss(y), mean(y), 0)
      print("node1 done, leaf")
      tree.frame <<- rbind(tree.frame, frame)
      return(frame)
      
    } else {
      
      maxxr <- max.cor(yy = y,sxs= xssrd)
      sprd <- split.rf(y, maxxr[[1]], min)
      if(is.null(sprd)) {
        frame[1,] <- c("<leaf>",nrow(xssrd),rss.tree(y),  mean(y), 0)
        print("node1 done, leaf")
        tree.frame <<- rbind(tree.frame, frame)
        return(frame)
      
        } else if(length(maxxr[[1]] < sprd[[1]]) < min |length(maxxr[[1]] >= sprd[[1]]) < min  ) {
        
        frame[1,] <- c("<leaf>",nrow(xssrd),rss.tree(y),mean(y), 0)
        print("node1 done, leaf, too small")
        tree.frame <<- rbind(tree.frame, frame)
        return(frame)
        }
      
      frame[1,] <- c(maxxr[[2]],
                     nrow(xssrd),
                     rss.tree(y),
                     mean(y),
                     sprd[[1]])
      # frame[2:4] <- as.numeric(frame[2:4])
    }
    print("node1 done")
    tree.frame <<- rbind(tree.frame, frame)
    return(frame)
  }

  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  y <- y[bootsample]
  min <- 5
  
  noder(y, xs, mtry,min)
  names(tree.frame) <- c("var", "n", "ypred", "dev", "splits.cutleft")
  tree[[1]] <- bootsample
  tree[[2]] <- tree.frame
  return(tree)
}

###################################BUILDING A FOREST#########################################

rforest <- function(form, d, mtry, ntree) { 
  #imputs: formula to be tested and the dset to test on. outputs: list,
  #contains: total oob error, trees, and bootsample for the tree  
  rforest <- list() #define empty list for the rf. length = ntree*2+1
  #where each tree is a list, bootsample + frame  
  xsi <- c()
  di <- data.frame()
  form <- as.character(form)
  yn <- form[2]
  yi <- as.numeric(names(d) == yn)
  y <- d[,yi]
  xs <- d[,-(yi)]
  
  while (length(rforest) < ntree) {
    rforest[[length(rforest)+1]] <- tree.rf(y,xs,mtry)
  }
  return(rforest)
}


###############################CONDITIONAL INFERENCE#########################################

