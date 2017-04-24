###########################################################################################
###########################################################################################
###################################### Ch4 Sim ############################################
###########################################################################################
###########################################################################################

library(ggplot2)
library(MASS)
library(gridExtra)

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
    cors[i,1] <- ifelse(sd(sxs[,i]) == 0,0, suppressWarnings(cor(yy, sxs[,i])))
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

split.rf <- function(y,x, min, parent_rss) {
  
  sp <- list()
  if(length(x) <= min) {
    
    return(NULL)
  } else {
    
    #gives a warning if it cannot find an optimal solution and instead gives the max point
    op.partition <- suppressWarnings(optimise(rss.node, interval = range(x), 
                                              upper = max(x), y=y,x=x, maximum = FALSE))
    
    # if(parent_rss < rss.leaf( y[x < op.partition[[1]]])){
    #   return(NULL)
    # }
    # 
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

tree.rf <- function(y,xs, mtry, bootsampled, min = 5){
  tree <- list()
  tree.frame <- list()
  tree.exc <- list()
  leaf.partitions <- list()
  bootsample <- seq(1, length(y))
  
  
  #if( length(y) != nrow(xs)){
  # return("error, length y != dim xs")
  # }
  
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
      if(is.null(ncol(xs))){
        yhatl <- y[xs < split]
        xsl <- xs[xs < split]
        yhatr <- y[xs >= split]
        xsr <- xs[xs >= split]
      }else{
        yhatl <- y[xs[,spliton[[1]]] < split]
        xsl <- xs[xs[,spliton[[1]]] < split,]
        yhatr <- y[xs[,spliton[[1]]] >= split]
        xsr <- xs[xs[,spliton[[1]]] >= split, ]
      }
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
    onex <- FALSE
    
    if(is.null(ncol(xs))){
      onex <- TRUE
    }
    if(onex == FALSE){
      xssrd <- xs[,sample(ncol(xs),mtry)]
      In <- names(xs) %in% names(xssrd)
      In <- !In 
      NotIn <- names(xs)[In]
    } else {
      xssrd <- as.data.frame(xs)
      NotIn <- NULL
    }
    
    frame <- data.frame("", 0,0,0,0)
    frame[,1] <- as.character(frame[,1])
    leaf.p <- data.frame(rep(0,2), rep(0,2), rep(0,2))
    
    
    if(length(y) < min) {
      frame[1,] <- c("<leaf>",nrow(xssrd), rss.leaf(y), mean(y), 0)
      
      tree.frame <<- rbind(tree.frame, frame)
      leaf.p <- sapply(xs, range)
      leaf.partitions[[length(leaf.partitions)+1]] <<- leaf.p
      return(frame)
      
    } else {
      if (onex){
        maxxr <- list(xs, "X")
      } else {
        maxxr <- max.cor(yy = y,sxs= xssrd)
      }
      sprd <- split.rf(y, maxxr[[1]], min, parent_rss = rss.leaf(y))
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
    tree.exc <<- rbind(tree.exc, NotIn)
    tree.frame <<- rbind(tree.frame, frame)
    return(frame)
  }
  
  if(bootsampled){
    bootsample <- sample(length(y), length(y), replace = TRUE)
    xs <- xs[bootsample,]
    y <- y[bootsample]
    
  }
  #min <-5
  noder(y, xs, mtry,min)
  names(tree.frame) <- c("var", "n", "dev", "ypred","split.cutleft")
  tree[[1]] <- bootsample
  tree[[2]] <- tree.frame
  tree[[3]] <- tree.exc
  return(tree)
}

###################################AUX TREE FUNCTIONS#######################################
predict.tree.rf <- function(t,xs) {
  t <- t[[2]] 
  first.split <- t[1,]
  onex <- FALSE
  
  if(sum(is.na(t$ypred)) > 0){
    return(NULL)
  }
  if(first.split$var == "<leaf>") {
    return(first.split$ypred)
  }
  
  t$n <- as.numeric(t$n)
  t$split.cutleft <- as.numeric(t$split.cutleft)
  rdn <- t$n[2]
  ldn <- t$n[1] - rdn 
  ldname <- as.numeric(row.names(t[t$n == ldn,]))
  
  if(is.null(ncol(xs))){
    onex <- TRUE
  }
  
  if(length(ldname) > 1) {
    
    right.daughter <- t[1:(ldname[length(ldname)]-1),]
    left.daughter <- t[ldname[length(ldname)]:nrow(t),]
    
  } else {
    
    right.daughter <- t[1:(ldname[length(ldname)]),]
    left.daughter <- t[ldname:nrow(t),]
  }
  
  
  j <- 0
  if(onex){
    predictions <- c(rep(100000, length(xs)))
  } else {
    predictions <- c(rep(100000, nrow(xs)))
  }
  ###SKIPPING CONDITION
  #there's no need to map out the whole tree and have it on file, each y just needs the tree to 
  #a. check the split condition, i.e. xs[,var] < split.cutleft 
  #b. if yes, go two ahead
  #b.2. if no, go one ahead
  #c. stop at leaf, pred[i] <- ypred
  
  if(onex) {
    for (i in 1: length(xs)) {
      xsh <- xs[i]
      if(xsh < first.split$split.cutleft){
        #left daughter
        ld <- left.daughter[1,]
        j <- 1
        while (ld$var != "<leaf>"){
          if (xsh < ld$split.cutleft) {
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
          if (xsh < rd$split.cutleft) {
            j <- j+2
          } else {
            j <- j+1
          }
          rd <- right.daughter[j,]
        }
        predictions[i] <- rd$ypred
      }    
    }
  } else {
    
    for(i in 1:(nrow(xs))){
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
  }
  return(predictions)
}

###################################GROWING A FOREST#########################################

rforest <- function(y, xs,  mtry, ntree) { 
  #imputs: formula to be tested and the dset to test on. outputs: list,
  #contains: total oob error, trees, and bootsample for the tree  
  rforest <- list() #define empty list for the rf. length = ntree*2+1
  #where each tree is a list, bootsample + frame  
  xsi <- c()
  di <- data.frame()
  # form <- as.character(form)
  # yn <- form[2]
  # yi <- as.numeric(names(d) == yn)
  # y <- d[,yi]
  # xs <- d[,(yi + 1 - 2)*-1]
  
  while (length(rforest) < ntree) {
    rforest[[length(rforest)+1]] <- tree.rf(y,xs,mtry, bootsampled = TRUE)
  }
  return(rforest)
}


###################################INFFOREST#################################################

infforest <- function(rf,y,xs) {
  v <- inftrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
  vi.frame <- v
  for(i in 2:(length(rf))){ ##random forest is a list, pairs of trees + bootsamples
    t <- rf[[i]]
    vi.frame <- rbind(vi.frame,inftrees(t, y[-t[[1]]], xs[-t[[1]],]))
    
    ##send tree off to inftrees with the -bootsampled data
    ##get back frame
  }
  
  ##generate distribution of vi's from frames - pval
  return(vi.frame)
}


infforest.test <- function(inf){
  p <- rep(10000, ncol(inf))
  mu <- c()
  sd <- c()
  for(i in 1:ncol(inf)){
    mu[i] <- mean(inf[,i])
    sd[i] <- sd(inf[,i])
    p[i] <- pnorm(0,mean = mu[i], sd = sd[i])
  }
  return(data.frame("var" =  names(as.data.frame(inf)), "P(var > 0)"= p, "mean" = mu, "sd" = sd))
}

##############################Y PERMUTED VARIABLE IMPORTANCE################################

yperm.inf <- function(rf,y,xs) {
  v <- ytrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
  v1 <- v
  for (i in 2:length(rf)){
    v1 <-  ytrees(rf[[i]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
    v <- rbind(v,v1)
  }
  
  return(v)
}

##############################CONDITIONAL VARIABLE IMPORTANCE################################

conditional.inf <- function(rf, y, xs) {
  v <- strobltrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
  v1 <- v
  for (i in 2:length(rf)){
    v1 <-  strobltrees(rf[[i]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
    v <- rbind(v,v1)
  }
  
  return(v)
}

###############################PERMUTED VARIABLE IMPORTANCE################################

permuted.inf <- function(rf, y, xs) {
  v <- briemantrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
  
  for (i in 2:length(rf)){
    v <- rbind(v, briemantrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],]))
  }
  
  return(v)
  
}

###################################INFTREES SUPPORT###########################################

inftrees <- function(t, y, xs) {
  
  rss.post <- rep(0, ncol(xs))
  vi <- rep(0, ncol(xs))
  
  ta <- t
  t <- t[[2]]
  
  # set.seed(1)
  #for(i in 1:ncol(xs)){
  #  rss.pre[i] <- VIdev(t, names(xs)[i])
  #}
  
  rss.pre <- sum((y - as.numeric(predict.tree.rf(ta, xs)))^2)/length(y)
  
  set.seed(1)
  for(i in 1:ncol(xs)){
    xs.for.permuting.i <- xs
    ti <- tree.rf(xs[,i], xs[,-i], mtry = ncol(xs[,-i]), bootsample = FALSE) 
    if(nrow(ti[[2]]) <= 3) {
      xs.for.permuting.i[,i] <- sample(xs.for.permuting.i[,i])
    } else {
      xi.partitions <- predict.tree.rf(ti, xs[,-i])
      xs.for.permuting.i$groups <- as.factor(xi.partitions)
      for (j in 1:length(levels(xs.for.permuting.i$groups))) {
        for(n in 1:ncol(xs)){
          xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],n] <-
            sample(xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],n], 
                   length(xs.for.permuting.i[xs.for.permuting.i$group == levels(xs.for.permuting.i$groups)[j],n]))
        }
      }
    }
    predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.for.permuting.i)
    rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
    #ex <- sum(names(xs)[i] ==  ta[[3]])/nrow(ta[[3]])
    # if(ex > .8){
    #  rss.post[i] <- rss.post[i] * (ex)
    # }
  }
  rss.post <- rss.post/median(abs(rss.post))
  names(rss.post) <- names(xs)
  return(rss.post)
}

ytrees <- function(t,y,xs) {
  d <- rbind(y,xs)
  rss.post <- rep(0, ncol(xs))
  ta <- t
  t <- t[[2]]
  predictions <- as.numeric(predict.tree.rf(ta, xs))
  rss.pre <- sum((y - predictions)^2)/length(y)
  xs.for.permuting.i <- xs
  yperm <- y
  for(i in 1:ncol(xs)){
    corrr <- sapply(xs[,-i],cor, y = xs[,i])
    xs.i <- xs[,-i]
    xs.i <-  (xs.i[,abs(corrr) > .2])
    xs.i <- cbind(xs[,i], xs.i)
    if(ncol(xs.i) == 1){
      yperm <- sample(y)
      rss.post[i] <- (rss.pre-sum((yperm - predictions)^2)/length(y) )
      xs.i <- xs
      yperm <- y
    } else {
      fr <- t[t$var %in% names(xs.i)[-1],]
      if(nrow(fr) == 0) {
        yperm <- sample(y)
        rss.post[i] <- (rss.pre-sum((yperm - predictions)^2)/length(y) )
        xs.i <- xs
        yperm <- y
      } else {
        for (j in 1:length(levels(as.factor(fr$var)))) {
          fri <- fr[fr$var == levels(as.factor(fr$var))[j],]
          for(n in 1:nrow(fri)){
            yperm[xs.i[,fri$var[n]] < fri$split.cutleft[n]] <- sample(y[xs.i[,fri$var[n]] < fri$split.cutleft[n]])
          }
        }
        rss.post[i] <- (rss.pre - sum((yperm - predictions)^2)/length(y))
        xs.i <- xs
        yperm <- y
      }
    }
  }
  if(sum(abs(rss.post)) == 0) {
    return(rss.post)
  } else{
    # rss.post <- rss.post/max(abs(rss.post))
    # names(rss.post) <- names(xs)
    return(rss.post)
  }
  
}


strobltrees <- function(t,y,xs) {
  rss.post <- rep(0, ncol(xs))
  
  ta <- t
  t <- t[[2]]
  
  predictions <- as.numeric(predict.tree.rf(ta, xs))
  rss.pre <- sum((y - predictions)^2)/length(y)
  
  xs.for.permuting.i <- xs
  
  for(i in 1:ncol(xs)){
    corrr <- sapply(xs[,-i],cor, y = xs[,i])
    xs.i <- xs[,-i]
    xs.i <-  (xs.i[,abs(corrr) > .2])
    xs.i <- cbind(xs[,i], xs.i)
    
    if(ncol(xs.i) == 1){
      xs.i <- sample(xs.i)
      xs.ii <- xs
      xs.ii[,i] <- xs.i[,1]
      names(xs.ii) <- names(xs)
      predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.ii)
      rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
      xs.i <- xs
    } else {
      
      fr <- t[t$var %in% names(xs.i)[-1],]
      if(nrow(fr) == 0) {
        xs.i <- sample(xs.i)
        xs.ii <- xs
        xs.ii[,i] <- xs.i[,1]
        names(xs.ii) <- names(xs)
        predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.ii)
        rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
        xs.i <- xs
      } else {
        for (j in 1:length(levels(as.factor(fr$var)))) {
          fri <- fr[fr$var == levels(as.factor(fr$var))[j],]
          for(n in 1:nrow(fri)){
            xs.i[xs.i[,fri$var[n]] < fri$split.cutleft[n] ,1] <- sample(xs.i[xs.i[,fri$var[n]] < fri$split.cutleft[n] ,1])
          }
        }
        xs.ii <- xs
        xs.ii[,i] <- xs.i[,1]
        names(xs.ii) <- names(xs)
        predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.ii)
        rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
        xs.i <- xs
        ex <- sum(names(xs)[i] ==  ta[[3]])/nrow(ta[[3]])
        if(ex > .8){
          rss.post[i] <- rss.post[i] * (ex)
        }
      }
    }
  }
  if(sum(abs(rss.post)) == 0) {
    return(rss.post)  
  } else{
    #rss.post <- rss.post/max(abs(rss.post))
    # names(rss.post) <- names(xs)
    return(rss.post)
  }
  
}

breimantrees <- function (t,y,xs){
  ta <- t
  t <- t[[2]]
  
  rss.pre <- sum((y - as.numeric(predict.tree.rf(ta, y, xs)))^2)
  
  rss.post <- rep(0, ncol(xs))
  xs.for.permuting <- xs
  for (i in 1:ncol(xs)){
    xs.for.permuting[,i] <- sample(xs.for.permuting[,i], replace = TRUE)
    predictions.for.y.xi.perm <- predict.tree.rf(ta,y,xs.for.permuting)
    rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2) - rss.pre)
    xs.for.permuting <- xs
  }
  rss.post <- rss.post/max(abs(rss.post))
  names(rss.post) <- names(xs)
  return(rss.post)
}



############SIM
rep(0, 12) -> mu #mean of each variable

diag(12) -> sigma #creates a 12 by 12 diagonal matrix with 1's down the diagonal

sigma -> sigma1
.9 -> sigma1[1,2]
.9 -> sigma1[2,1] #adding the block correlation between the first four variables
.9 -> sigma1[1,3]
.9 -> sigma1[3,1]
.9 -> sigma1[3,2]
.9 -> sigma1[2,3]
.9 -> sigma1[1,4]
.9 -> sigma1[4,1]
.9 -> sigma1[2,4]
.9 -> sigma1[4,2]
.9 -> sigma1[4,3]
.9 -> sigma1[3,4]

1000 -> n #number of observations
set.seed(1)
mvrnorm(n =n, mu = mu, Sigma = sigma1) -> sim1000 #sampling from the multivariate normal 
c(5,5,2,0,-5,-5,-2,0,0,0,0,0) -> bts #these are the betas 

rnorm(1000, mean = 0, sd = .5) -> e #error terms

rep(0, 1000) -> ys #init a vector of zeros

for( i in 1:1000){ #go row by row and create the ys based on the function of e, bts, and sim1000
  ys[i] <- sim1000[i,1]*bts[1]+
    sim1000[i,2]*bts[2]+ 
    sim1000[i,3]*bts[3]+ 
    sim1000[i,4]*bts[4]+ 
    sim1000[i,5]*bts[5]+ 
    sim1000[i,6]*bts[6]+
    sim1000[i,7]*bts[7]+ 
    sim1000[i,8]*bts[8]+ 
    sim1000[i,9]*bts[9]+ 
    sim1000[i,10]*bts[10]+ 
    sim1000[i,11]*bts[11]+ 
    sim1000[i,12]*bts[12] +
    e[i]
}

d1 <- as.data.frame(sim1000)
d1$y <- ys



n <- 10000
c(5,5,2,0,-5,-5,0,0,0,0,0,0) -> bts #these are the betas 
rnorm(1000, mean = 0, sd = .5) -> e #error terms



x1 <- rnorm(1000)

x2 <- 2*log(abs(x1)) + rnorm(1000)

x3 <- 2*log(abs(x2)) + rnorm(1000)

x4 <- 2*log(abs(x3)) + rnorm(1000)

sigma <- diag(8)

mvrnorm(n =n, mu = mu, Sigma = sigma) -> simx5x12

data.frame(x1,x2,x3,x4,simx5x12) -> d2

c("X1","X2","X3","X4","X5","X6", "X7","X8","X9","X10","X11","X12") -> names(d2)

rep(0, 1000) -> ys #init a vector of zeros

for( i in 1:1000){ #go row by row and create the ys based on the function of e, bts, and d2
  ys[i] <- d2[i,1]^2*bts[1]+
    d2[i,2]^2*bts[2]+ 
    d2[i,3]^2*bts[3]+ 
    d2[i,4]^2*bts[4]+ 
    d2[i,5]*bts[5]+ 
    d2[i,6]*bts[6]+
    d2[i,7]*bts[7]+ 
    d2[i,8]*bts[8]+ 
    d2[i,9]*bts[9]+ 
    d2[i,10]*bts[10]+ 
    d2[i,11]*bts[11]+ 
    d2[i,12]*bts[12] +
    e[i]
}

d2$y <- ys


################################MTRY PLOT 
set.seed(1)

rf1 <- rforest( d1$y[1:100], d1[1:100, 1:12],mtry = 4,ntree = 300)
inf1 <- infforest(rf1,  d1$y[1:100], d1[1:100, 1:12])
inf1[inf1 == "NaN"] <- NA
inf1[inf1 == "Inf"] <- NA
inf11 <- inf1[complete.cases(inf1),]
v <- as.data.frame(inf11)

set.seed(1)
r12 <- rforest( d1$y[1:100], d1[1:100, 1:12],mtry = 6,ntree = 300)
v12 <- as.data.frame(infforest(r12,  d1$y[1:100], d1[1:100, 1:12]))
v12[v12 == "NaN"] <- NA
v12[v12 == "Inf"] <- NA
v12 <- v12[complete.cases(v12),]
set.seed(1)
r13 <- rforest( d1$y[1:100], d1[1:100, 1:12],mtry = 8,ntree = 300)
v13 <- as.data.frame(infforest(r13,  d1$y[1:100], d1[1:100, 1:12]))
v13[v13 == "NaN"]<- NA
v13[v13 == "Inf"]<- NA
v13 <- v13[complete.cases(v13),]


b <- data.frame("variable" = names(as.data.frame(inf11)),
                "median" = sapply(as.data.frame(inf11), median),
                "betas" = bts)


b$variable <- factor(b$variable, levels = b$variable)


c <- data.frame("variable" = names(as.data.frame(v12)),
                "median" = sapply(as.data.frame(v12), median),
                "iqr" = sapply(as.data.frame(v12), IQR),
                "betas" = bts)
d <- data.frame("variable" = names(as.data.frame(v13)),
                "median" = sapply(as.data.frame(v13), median),
                "iqr" = sapply(as.data.frame(v13), IQR),
                "betas" = bts)


c$variable <- factor(c$variable, levels = c$variable)

d$variable <- factor(d$variable, levels = d$variable)


b$q1 <- rep(0,12)
b$q3 <- rep(0,12)
for(i in 1:12){
  b$q1[i] <- quantile(inf11[,i], 1/4)
  b$q3[i] <- quantile(inf11[,i], 3/4)
}

c$q1 <- rep(0,12)
c$q3 <- rep(0,12)
for(i in 1:12){
  c$q1[i] <- quantile(v12[,i], 1/4)
  c$q3[i] <- quantile(v12[,i], 3/4)
}
d$q1 <- rep(0,12)
d$q3 <- rep(0,12)
for(i in 1:12){
  d$q1[i] <- quantile(v13[,i], 1/4)
  d$q3[i] <- quantile(v13[,i], 3/4)
}
p <- ggplot(data = b, aes( x = variable, y = median, group = 1)) +
  geom_point() +
  ylab(" ")+
  xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_line()+
  geom_ribbon(ymin = b$q1 , ymax = b$q3, data = b, alpha = .4)+
  scale_y_continuous(limits = c(0,1))


q<- ggplot(data = c, aes( x = variable, y = median, group = 1)) +
  geom_point() +
  ylab(" ")+
  xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_line()+
  geom_ribbon(aes(ymin = c$q1, ymax = c$q3), alpha = .4)+
  scale_y_continuous(limits = c(0,1))


r <- ggplot(data = d, aes( x = variable, y = median, group = 1)) +
  geom_point()+
  ylab(" ")+
  xlab(" ")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_line()+
  geom_ribbon(aes(ymin = d$q1, ymax = d$q3), alpha = .4)+
  scale_y_continuous(limits = c(0,1))

grid.arrange(p,q,r, ncol=3)