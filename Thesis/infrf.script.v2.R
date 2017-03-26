
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

t1 <- tree.rf(y,xs,mtry)
r <- rforest(form, d, mtry, ntree = 10)

p <- predict.tree.rf(t1, y[-t1[[1]]], xs[-t1[[1]],])

vi <- inftrees(t1, y[-t1[[1]]], xs[-t1[[1]],])

virf <- infforest(r,y,xs)
virfStrobl <- conditional.inf(r,y,xs)
###################################GROWING A TREE#########################################

rss.node <- function(i,y,x) { 
  yleftpred <- mean(y[x< i])
  yrightpred <- mean(y[x >= i])
  rssl <- sum((y[x< i]-yleftpred)^2)
  rssr <- sum((y[x >= i] - yrightpred)^2)
  return(rssl + rssr)
}

rss.leaf <- function(a){
  a1 <- rep(mean(a), length(a))
  sum((a-a1)^2)
}

max.cor <- function(yy,sxs){
 
  max <- list()
  cors <- cbind(rep(0,ncol(sxs)),seq(1,ncol(sxs)))
  for(i in 1:ncol(sxs)){
    cors[i,1] <- ifelse(sd(sxs[,i]) == 0,0, cor(yy, sxs[,i]))
  }

  cors[is.na(cors[,1]), 1] <- 0

  if (sum(cors[,1] == max(abs(cors[,1]), na.rm = TRUE), na.rm = TRUE) > 1) {

    if(sum(abs(cors[,1])) == 0){

      return(NULL)
    }
    max1 <- sample(as.numeric(cors[,1] == max(abs(cors[,1]), na.rm = TRUE)), 1)  

  } else {

    max1 <- cors[abs(cors[,1]) == max(abs(cors[,1]), na.rm = TRUE),2]
  }

  
  max[[length(max)+1]] <- sxs[,max1]
  max[[length(max)+1]] <- names(sxs)[max1]
  return(max)
  
}

split.rf <- function(y,x, min) {

  sp <- list()
  if(length(x) <= min) {

    return(NULL)
  } else {
    op.partition <- optimise(rss.node, interval = range(x), y=y,x=x, maximum = FALSE, tol = .01)
    
    
    if (length(x[x < op.partition[[1]]]) < min | length(x[x >= op.partition[[1]]]) < min){

      return(NULL)
    }
    sp[[length(sp)+1]] <- op.partition[[1]] #partition on x
    sp[[length(sp)+1]]<-  mean(y[x< op.partition[[1]]]) #ypred left daughter
    sp[[length(sp)+1]]<- mean(y[x >= op.partition[[1]]]) #ypred right daughter
    sp[[length(sp)+1]]<- x[x < op.partition[[1]]] #xsleftdaughter
    sp[[length(sp)+1]]<- x[x >= op.partition[[1]]] #xsrightdaughter
    sp[[length(sp)+1]] <- rss.leaf( y[x < op.partition[[1]]]) #ldrss

    return(sp)
  }
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  tree.frame <- list()
  leaf.partitions <- list()
  
  noder <- function(y,xs,mtry, min){
    
    ##what i want this function to do:
    ###given y,xs (like an interval of them):
    ####1. find the split between them 

    frame <- node1(y,xs, mtry, min)
    #df <-  rbind(df,frame)
    ####2. call itself on each side of the split
    split <- as.numeric(frame[5])
    spliton <- frame[1]
    
    if(frame[1] == "<leaf>"){
      return(frame)
    } else {

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

    xssrd <- xs[,sample(ncol(xs),mtry)]
    frame <- data.frame("", 0,0,0,0)
    frame[,1] <- as.character(frame[,1])
    leaf.p <- data.frame(rep(0,2), rep(0,2), rep(0,2))
    
    
    if(length(y) < min) {
      frame[1,] <- c("<leaf>",nrow(xssrd), rss(y), mean(y), 0)
 
      tree.frame <<- rbind(tree.frame, frame)
      leaf.p <- sapply(xs, range)
      leaf.partitions[[length(leaf.partitions)+1]] <<- leaf.p
      return(frame)
      
    } else {
      
      maxxr <- max.cor(yy = y,sxs= xssrd)
      sprd <- split.rf(y, maxxr[[1]], min)
      if(is.null(sprd)) {
        frame[1,] <- c("<leaf>",nrow(xssrd),rss.leaf(y),  mean(y), 0)

        tree.frame <<- rbind(tree.frame, frame)
        leaf.p <- sapply(xs, range)
        leaf.partitions[[length(leaf.partitions)+1]] <<- leaf.p
        return(frame)
      
        } else if(length(maxxr[[1]] < sprd[[1]]) < min |length(maxxr[[1]] >= sprd[[1]]) < min  ) {
        
        frame[1,] <- c("<leaf>",nrow(xssrd),rss.leaf(y), mean(y), 0)

        tree.frame <<- rbind(tree.frame, frame)
        leaf.p <- sapply(xs, range)
        leaf.partitions[[length(leaf.partitions)+1]] <<- leaf.p
        return(frame)
        }
      
      frame[1,] <- c(maxxr[[2]],
                     nrow(xssrd),
                     rss.node(sprd[[1]], y, maxxr[[1]]),
                     mean(y),
                     sprd[[1]])
      # frame[2:4] <- as.numeric(frame[2:4])
    }

    tree.frame <<- rbind(tree.frame, frame)
    return(frame)
  }

  bootsample <- sample(length(y), length(y), replace = TRUE)
  xs <- xs[bootsample,]
  y <- y[bootsample]
  min <- 5
  
  noder(y, xs, mtry,min)
  names(tree.frame) <- c("var", "n", "ypred", "dev", "split.cutleft")
  tree[[1]] <- bootsample
  tree[[2]] <- tree.frame
  return(tree)
}

###################################AUX TREE FUNCTIONS#######################################

predict.tree.rf <- function(t,y,xs) {
  t <- t[[2]] 
  t$n <- as.numeric(t$n)
  t$split.cutleft <- as.numeric(t$split.cutleft)
  rdn <- t$n[2]
  ldn <- t$n[1] - rdn 
  ldname <- as.numeric(row.names(t[t$n == ldn,]))
  right.daughter <- t[1:(ldname[length(ldname)]-1),]
  left.daughter <- t[ldname[length(ldname)]:nrow(t),]
  first.split <- t[1,]
  j <- 0
  predictions <- c(rep(100000, length(y)))
###SKIPPING CONDITION
    #there's no need to map out the whole tree and have it on file, each y just needs the tree to 
    #a. check the split condition, i.e. xs[,var] < split.cutleft 
    #b. if yes, go two ahead
    #b.2. if no, go one ahead
    #c. stop at leaf, pred[i] <- ypred
    
  for(i in 1:(length(y))){
    yh <- y[i]
    xsh <- xs[i,]
    
    if(xsh[,first.split$var] < first.split$split.cutleft){
      #left daughter
      ld <- left.daughter[1,]
      j <- 1
      while (ld$var != "<leaf>"){
        if (xsh[,ld$var] < ld$split.cutleft) {
          j <- j+2
        } else {
          j <- j+1
        }
        ld <- left.daughter[j,]
      }
      predictions[i] <- ld$ypred
    } else {
      #right daughter
      rd <- right.daughter[2,]
      j <- 2
      while (rd$var != "<leaf>"){
        if (xsh[,rd$var] < rd$split.cutleft) {
          j <- j+2
        } else {
          j <- j+1
        }
        rd <- right.daughter[j,]
      }
      predictions[i] <- rd$ypred
    }
  }
  return(predictions)
}

###################################GROWING A FOREST#########################################

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


###################################INFFOREST#################################################

infforest <- function(rf,y,xs) {
  vi.frame <- data.frame()
  for(i in 1:(length(rf))){ ##random forest is a list, pairs of trees + bootsamples
    t <- rf[[i]]
    vi.frame <- rbind(vi.frame,inftrees(t, y[-t[[1]]], xs[-t[[1]],]))
    
    ##send tree off to inftrees with the -bootsampled data
    ##get back frame
  }
  
    ##generate distribution of vi's from frames - pval
  names(vi.frame) <- names(xs)
  return(vi.frame)
}

##############################CONDITIONAL VARIABLE IMPORTANCE################################

conditional.inf <- function(rf, y, xs) {
  vi.frame <- data.frame()
  for(i in 1:(length(rf))){ ##random forest is a list, pairs of trees + bootsamples
    t <- rf[[i]]
    vi.frame <- rbind(vi.frame,strobltrees(t, y[-t[[1]]], xs[-t[[1]],]))
  }
  names(vi.frame) <- names(xs)
  return(vi.frame)
}

###################################INFTREES SUPPORT###########################################

VIdev<- function(t2, x){
  return(sum(as.numeric(t2[t2$var == x, 3])))
}

inftrees <- function(t, y, xs) {

  rss.post <- rep(0, ncol(xs))
  vi <- rep(0, ncol(xs))
  rss.pre <- rep(0, ncol(xs))
  
  ta <- t
  t <- t[[2]]
  
 # set.seed(1)
  #for(i in 1:ncol(xs)){
  #  rss.pre[i] <- VIdev(t, names(xs)[i])
  #}
  
  set.seed(1)
  for(i in 1:ncol(xs)){
    xs.for.permuting.i <- xs
    ti <- tree.rf(xs[,i], xs[,-i], mtry = ncol(xs[,-i])) 
    xi.partitions <- predict.tree.rf(ti, xs[,i],  xs[,-i])
    xs.for.permuting.i$groups <- as.factor(xi.partitions)
    for (j in 1:length(levels(xs.for.permuting.i$groups))) {
      xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],i] <-
        sample(xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],i], 
               replace = TRUE)
    }
    
    predictions.for.y.xi.perm <- predict.tree.rf(ta,y,xs.for.permuting.i)
    rss.post[i] <- sum((y - as.numeric(predictions.for.y.xi.perm))^2) 
  }
  rss.post <- rss.post/sum(rss.post)
  names(rss.post) <- names(xs)
  return(rss.post)
}

strobltrees <- function(t,y,xs) {
  rss.post <- rep(0, ncol(xs))
  vi <- rep(0, ncol(xs))
  rss.pre <- rep(0, ncol(xs))
  
  ta <- t
  t <- t[[2]]

  set.seed(1)
  for(i in 1:ncol(xs)){
    xs.for.permuting.i <- xs
#    ti <- tree.rf(xs[,i], xs[,-i], mtry = ncol(xs[,-i])) 
    xi.partitions <- predict.tree.rf(ta, y, xs)
    xs.for.permuting.i$groups <- as.factor(xi.partitions)
    
    for (j in 1:length(levels(xs.for.permuting.i$groups))) {
      
      xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],i] <-
              sample(xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],i], 
               replace = TRUE)
    }
    
    predictions.for.y.xi.perm <- predict.tree.rf(ta,y,xs.for.permuting.i)
    rss.post[i] <- sum((y - as.numeric(predictions.for.y.xi.perm))^2) 
  }
  rss.post <- rss.post/sum(rss.post)
  names(rss.post) <- names(xs)
  return(rss.post)
}
