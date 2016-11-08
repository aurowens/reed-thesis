#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tree)
library(ggplot2)

x1 <- rnorm(1000)

shinyServer(function(input, output) {
  
  sampled <- reactive({
    sampled <- data.frame(y =  input$coefX1*x1 + input$coefX2*x2 + rnorm(1000), 
                       x1, 
                       x2 =  input$r * x1 + rnorm(mean = 0, sd = sqrt(1 - input$r^2), 1000)) 
  })
  
  output$text1 <- renderText({ 
    paste("Your function is y =", input$coefX1, "* X1+", input$coefX2,"*X2 + E")
  })
  
  output$MSE <- renderPlot({
    e <- rep(0,1000)
    data <- sampled()
    set.seed(444)
    for(i in 1:1000){
      t <- tree(y ~x1 + x2, data)
      yhat <- predict(t)
    e[i] <- mean((data$y - yhat)^2)
    }
    
    ggplot(data = as.data.frame(e)) + 
      geom_histogram(aes(x = e), fill = "#a16c01") + 
      ggtitle("The Distribution of MSE of 1000 Trees") + 
      xlab("MSE")+
      ylab(" ")
  
  })
  output$Splits <- renderPlot({
    s <- rep(0,1000)
    data <- sampled()
    #set.seed(444)
    for(i in 1:1000){
      t <- tree(y ~., data)
      splits.cutleft <- t$frame$splits[,1]
      splits.cutleft <- gsub("<", "", splits.cutleft)
      splits.cutleft[splits.cutleft == ""] <- NA
      s[i]<- length(splits.cutleft[complete.cases(splits.cutleft)])
    }
    
    ggplot(data = as.data.frame(s)) + 
      geom_histogram(aes(x = s), fill = "#99183c" ) +
      ggtitle("Distribution of the Number of Splits in 1000 Trees")
  })
  output$FirstSplit <- renderPlot(({
    f <- rep("none",1000)
    data <- sampled()
    #set.seed(444)
    for(i in 1:1000){
      t <- tree(y ~., data)
      f[i] <- as.character(t$frame$var[1])
    }
    
    ggplot(data = as.data.frame(f)) + 
      geom_bar(aes(x = f), fill = "#800000") +
      ggtitle("The First Split in 1000 Trees")
    
  }))

})
