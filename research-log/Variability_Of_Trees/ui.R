
library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("How Does Correlation Affect the Complexity of Trees?"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       sliderInput("r",
                   "Correlation between X1 and X2 ",
                   min = 0,
                   max = 1,
                   value = 0),
       numericInput("coefX1", "Coefficient of X1", value = 0),
       numericInput("coefX2", "Coefficient of X2", value = 0),
       textOutput("text1")
    ),
    mainPanel(
        tabsetPanel(
          tabPanel("MSE", plotOutput("MSE")),
          tabPanel("Number of Splits", plotOutput("Splits")),
          tabPanel("Variable with the First Split",  plotOutput("FirstSplit"))
        )
      )
    )
  )
)
