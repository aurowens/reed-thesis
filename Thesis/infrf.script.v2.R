
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
r <- rforest(form, d, mtry, ntree = 50)
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
  print("y")
  print((yy))
  print("xs")
  print(sxs)
  max <- list()
  cors <- cbind(rep(0,ncol(sxs)),seq(1,ncol(sxs)))
  for(i in 1:ncol(sxs)){
    cors[i,1] <- ifelse(sd(sxs[,i]) == 0, -100000, cor(yy, sxs[,i]))
  }
  print("is na cor?")
  print(cors)
  cors[is.na(cors[,1]), 1] <- 0
  print(cors)
  print("cor finding done, lets take a look and see if we have more than 1 max")
  print(sum(cors[,1] == max(cors[,1], na.rm = TRUE), na.rm = TRUE) > 1)
  if (sum(cors[,1] == max(cors[,1], na.rm = TRUE), na.rm = TRUE) > 1) {
    print("we do!")
    max1 <- sample(as.numeric(cors[,1] == max(cors[,1], na.rm = TRUE)), 1)  
    print("our winning var is")
    print(max1)
  } else {
    print("we don't")
    max1 <- cors[cors[,1] == max(cors[,1], na.rm = TRUE),2]
    print("our winning var is")
    print(max1)
  }
 # print(max1)
  
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


###################################INFFOREST#################################################

infforest <- function(rf,d) {
  for(i in 1:(length(rf)/2)){ ##random forest is a list, pairs of trees + bootsamples
    ##send tree off to inftrees with the -bootsampled data
    ##get back frame
  }
    ##generate distribution of vi's from frames - pval
}

###################################INFTREES SUPPORT###########################################

##INFTREES was created to be used with `tree` objects, will need to be reworked 
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
