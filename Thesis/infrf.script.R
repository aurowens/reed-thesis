
############################################################################################
#####################################INFFORESTS#############################################
############################################################################################
############################################################################################

rss <- function(b,a) { 
  a <- rep(mean(a),length(a))
  sum((a-b)^2)
}


splitter <- function(y,x, min) {
  sp <- list()
  op.partition <- optimise(rss, interval = range(x), a = y, maximum = FALSE)
  if(length(x[x < op.partition[[1]]]) < min | length(x[x >= op.partition[[1]]]) < min) {
    return(NULL)
    } else {
      sp[[length(sp)+1]] <- op.partition[[1]] #partition on x
      sp[[length(sp)+1]]<-  mean(y[x< op.partition[[1]]]) #ypred left daughter
      sp[[length(sp)+1]]<- mean(y[x >= op.partition[[1]]]) #ypred right daughter
      sp[[length(sp)+1]]<- x[x < op.partition[[1]]] #xsleftdaughter
      sp[[length(sp)+1]]<- x[x >= op.partition[[1]]] #xsrightdaughter
      sp[[length(sp)+1]] <- rss(x[x < op.partition[[1]]], y[x < op.partition[[1]]]) #ldrss
      return(sp)
    }
}


node <- function (){
  
}

tree.rf <- function(y,xs, mtry){
  tree <- list()
  
  xss <- xs[,sample(ncol(xs),mtry)]
  cors <- cbind(rep(NA,ncol(xss)),seq(1,ncol(xss)))
 
################ONE SPLIT AND TWO DAUGHTERS
#######FIRST SPLIT
  for(i in 1:ncol(xss)){
    cors[i,1] <- cor(y,xss[,i])}
  maxx <- xss[,cors[cors[,1] == max(cors[,1]),2]]
  sp <- splitter(y, maxx, min = 4)
  
  if(is.null(sp)) {
     return(data.frame("ypred" = mean(y), "rss" = rss(xss,y)))
  } else{
    frame <- data.frame("var" = names(xss)[cors[cors[,1] == max(cors[,1]),2]], 
                        "n"= nrow(xss),
                        "yval" = mean(y),
                        "rss" = rss(y,maxx))
    split <- c(paste("x <", sp[[1]]))
       
    frame$split <- split
    frame$var <- as.character(frame$var) 
  }
######LEFT DAUGHTER  
  spld <- sp
  maxxl <- maxx
  sprd <- sp
  maxxr <- maxx

    xssld <- xs[maxxl < spld[[1]],]
    xssld <- xssld[,sample(ncol(xssld),mtry)]
    yld <- y[maxxl < spld[[1]]]
    for(i in 1:ncol(xssld)){
      cors[i,1] <- cor(yld,xssld[,i])}
    maxxl <- xssld[,cors[cors[,1] == max(cors[,1]),2]]
    spld <- splitter(yld, maxxl, min = 4)
    
    if(is.null(spld)) {
      frame[nrow(frame)+1,] <- c("<leaf>",nrow(xssld), mean(yld), rss(xssld, yld), "")
      
    } else {
      frame[nrow(frame)+1,] <- c(names(xssld)[cors[cors[,1] == max(cors[,1]),2]],
                                 nrow(xssld),
                                 mean(yld), 
                                 rss(xssld, yld),
                                 paste("x <", spld[[1]]))
    }
    
##right daughter
    xssrd <- xs[maxxr >= sprd[[1]],]
    xssrd <- xssrd[,sample(ncol(xssrd),mtry)]
    yrd <- y[maxxr >= sprd[[1]]]
    for(i in 1:ncol(xssrd)){
      cors[i,1] <- cor(yrd,xssrd[,i])}
    maxxr <- xssrd[,cors[cors[,1] == max(cors[,1]),2]]
    sprd <- splitter(yrd, maxxr, min = 4)
    
    if(is.null(sprd)) {
      frame[nrow(frame)+1,] <- c("<leaf>",nrow(xssrd), mean(yrd), rss(maxxr, yrd), "")
    } else {
      frame[nrow(frame)+1,] <- c(names(xssrd)[cors[cors[,1] == max(cors[,1]),2]],
                                 nrow(xssrd),
                                 mean(yrd),
                                 rss(xssrd, yrd),
                                 paste("x <", sprd[[1]]))
      
    }
    
    
  
   tree <- frame
  return(tree)
}

##steps for tree fitting
###each split: choose var, choose split, record var, rss, and split

############################################################################################
#####################################GRAVEYARD##############################################
############################################################################################
############################################################################################


split.rf <- function(y,x,min) {
  
  #right split, x < split
  xsplit <- c()
  xsplit[1] <- sort(x)[min]
  ymean <- mean(y[x <= xsplit])
  rssl0 <- rss(y[x <= xsplit], rep(ymean, length(y[x <= xsplit])))
  rssr0 <- rss(y[x > xsplit], rep(ymean, length(y[x > xsplit])))
  rssl <- rssl0
  rssr <- rssr0
  rss1 <- rssl + rssr
  rssl[2] <- 10000000000
  rssr[2] <- 10000000000
  rss1[2] <- 10000000000
  i <- 2
  split <- c()
  while(length(x[x > xsplit[i-1]]) >= min){ 
    xsplit[i] <- min(sort(x)[sort(x) > xsplit[i-1]])
    ymean[i] <- mean(y[x <= xsplit[i]])
    rssl[i] <- rss(y[x <= xsplit[i]], rep(ymean, length(y[x <= xsplit[i]])))
    rssr[i] <- rss(y[x > xsplit[i]], rep(ymean, length(y[x > xsplit[i]])))
    rss1[i] <- rssl[i] + rssr[i]
    i <- i + 1
  }
  split1 <- c("x < " = xsplit[rss1 == min(rss1)], "ypred" = ymean[rss1 == min(rss1)], "rss" =min(rss1))
  return(split1)
}

for(i in 5:145){
  rsses[i-5] <- split.rf(y = iris$Sepal.Length, x = iris$Petal.Length, min = i) [3]
}


#########################TO DO:
  
#####  Imma have to write my own random forest code, y'all 

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
f <- function(a) {
  c <- 0
  if(length(a) > 1){
    a <- a[c(1:(length(a)-1))] 
    c <- c + sum(a)
    f(a)
  } else { 
    return(c)
  }
}
