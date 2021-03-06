
############################################################################################
#####################################INFFORESTS#############################################
############################################################################################
############################################################################################

########################################TESTING#############################################

#data("iris")
#d <- iris[1:4]
#y <- iris$Sepal.Length
#xs <- iris[2:4]
#mtry <- 2
# form <- as.formula("Sepal.Length ~.")
# 
# ######from chap2.Rmd
# ##d1
# 
# d <- d2[1:100,]
# xs <- d2[1:100, 1:12]
# y <- d2$y[1:100]
#mtry = 6
# form <- as.formula("y~.")
# 
#t1 <- tree.rf(y,xs,mtry, bootsampled = TRUE)
# 
# start <- sys.time()
# r <- rforest(y,xs, mtry, ntree = 100)
# print(start - Sys.time())
# 
# p <- predict.tree.rf(r[[1]], y[-r[[1]][[1]]], xs[-r[[1]][[1]],])
# 
# vi <- inftrees(r[[2]], y[-r[[2]][[1]]], xs[-r[[2]][[1]],])
# 
# virf <- infforest(r,y,xs)
# virfStrobl <- conditional.inf(r,y,xs)
# virf_permuted <- permuted.inf(r,y,xs)
# 
# 
# 
# strobltrees(r[[3]],y[-r[[3]][[1]]],xs[-r[[3]][[1]],])

r1 <- rforest(d2$y[1:200], d2[1:200,1:12], mtry = 8, ntree = 300)
v1 <- conditional.inf(r1, d2$y[1:200], d2[1:200,1:12])

v22 <- v1[, v1== 0]

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
      if(xsh[,first.split$var] < as.numeric(first.split$split.cutleft)){
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

rssforest <- function(rf, y, xs){
  rss <- sum((y[-rf[[1]][[1]]] - as.numeric(predict.tree.rf(rf[[1]],xs[-rf[[1]][[1]],])))^2)
  rss.frame <- rss
  for(i in 2:(length(rf))){ ##random forest is a list, pairs of trees + bootsamples
    t <- rf[[i]]
    rss.frame <- rbind(rss.frame,sum((y[-rf[[i]][[1]]] - as.numeric(predict.tree.rf(t,xs[-rf[[i]][[1]],])))^2))
    
    ##send tree off to inftrees with the -bootsampled data
    ##get back frame
  }
  
  ##generate distribution of vi's from frames - pval
  return(mean(rss.frame))
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
  v <- breimantrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],])
  
  for (i in 2:length(rf)){
    v <- rbind(v, breimantrees(rf[[1]],y[-rf[[1]][[1]]],xs[-rf[[1]][[1]],]))
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
    # if (is.null(ex)){
    #   rss.post[i] <- rss.post[i] *0
    # } else if(ex > .7) {
    # rss.post[i] <- rss.post[i] * (ex)
    # }
  }
  rss.post <- rss.post/max(abs(rss.post))
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
    xs.i <-  (xs.i[,abs(corrr) >= .2])
    xs.i <- cbind(xs[,i], xs.i)
    
    if(ncol(xs.i) == 1 | is.null(dim(xs.i))){
      xs.i <- sample(xs.i)
      xs.ii <- xs
      xs.ii[,i] <- xs.i
      names(xs.ii) <- names(xs)
      predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.ii)
      rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
      xs.i <- xs
    } else {
      fr <- t[t$var %in% names(xs.i)[-1],]
      if(nrow(fr) < 3) {
        xs.i[,1] <- sample(xs.i[,1])
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
            xs.i[xs.i[,fri$var[n]] < fri$split.cutleft[n] ,1] <- sample(xs.i[xs.i[,fri$var[n]] < fri$split.cutleft[n] ,1], length(xs.i[xs.i[,fri$var[n]] < fri$split.cutleft[n] ,1]))
          }
        }
        xs.ii <- xs
        xs.ii[,i] <- xs.i[,1]
        names(xs.ii) <- names(xs)
        predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.ii)
        rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2)/length(y) - rss.pre)
        xs.i <- xs
        # ex <- sum(names(xs)[i] ==  ta[[3]])/nrow(ta[[3]])
        # if(ex > .8){
        #   rss.post[i] <- rss.post[i] * (ex)
        # }
      }
    }
  }
  if(sum(abs(rss.post)) == 0) {
    return(rss.post)  
  } else{
    rss.post <- rss.post/max(abs(rss.post))
    # names(rss.post) <- names(xs)
    return(rss.post)
  }
  
}

breimantrees <- function (t,y,xs){
  ta <- t
  t <- t[[2]]
  
  rss.pre <- sum((y - as.numeric(predict.tree.rf(ta, xs)))^2)
  
  rss.post <- rep(0, ncol(xs))
  xs.for.permuting <- xs
  for (i in 1:ncol(xs)){
    xs.for.permuting[,i] <- sample(xs.for.permuting[,i], replace = TRUE)
    predictions.for.y.xi.perm <- predict.tree.rf(ta,xs.for.permuting)
    rss.post[i] <- abs(sum((y - as.numeric(predictions.for.y.xi.perm))^2) - rss.pre)
    xs.for.permuting <- xs
  }
  rss.post <- rss.post/max(abs(rss.post))
  names(rss.post) <- names(xs)
  return(rss.post)
}