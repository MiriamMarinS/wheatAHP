# R Shiny.

library(shiny)
library(ggplot2)
library(ggpattern)
library(grid)
library(gridExtra)

Species <- c("BW", "DW", "HT")
Rye <- c("Y", "N")

# See above for the definitions of ui and server
ui <- fluidPage(
  titlePanel(p("Score AHP", style = "color:black")),
  sidebarLayout(
    sidebarPanel(
      sliderInput("range", "Range", value = c(0, 0.05), min = 0, max = 0.05),
      checkboxGroupInput("genotype", "Genotypes", Species),
      checkboxGroupInput("rye", "Rye", Rye),
      img(
        src = "Logo.png",
        width = "120px", height = "50px"
      )
    ),
    mainPanel(
      plotOutput("ggplot")
    )
  )
)

server <- function(input, output) {
  score <- reactive({
    read.table("./results/table_scores.txt", sep = "\t", header = TRUE, row.names = 1)
  })
  
  rawdata <- reactive({
    read.csv("./Input/Data_matrix_AHP.csv", sep = ";", header = TRUE, row.names = 1)
  })
  
  output$bins1 <- renderUI({
    sliderInput("bins1", 
                h3("Bin width #1 "))
                #min = 1,
                #max = max(x),
                #value = (1 + max(x))/10)
  })
  
  output$bins2 <- renderUI({
    sliderInput("bins2", 
                h3("Bin width #2 "))
                #min = 1,
                #max = max(x),
                #value = (1 + max(x))/10)
  })
  
  output$bins3 <- renderUI({
    sliderInput("bins3", 
                h3("Bin width #3 "))
                #min = 1,
                #max = max(x),
                #value = (1 + max(x))/10)
  })
    
    
  output$ggplot <- renderPlot({
    
    score <- score()
    rawdata <- rawdata()
    
    scoredata <- merge(score[, c("mean", "sd")], rawdata[, c("genus", "rye")], by = "row.names")
    names(scoredata)[names(scoredata) == 'Row.names'] <- 'Genotypes'
    
    scoredata <- subset(scoredata, scoredata$genus %in% input$genotype)
    scoredata <- subset(scoredata, scoredata$rye %in% input$rye)
    scoredata <- scoredata[scoredata$mean >= input$range[1] & scoredata$mean <= input$range[2],]
    
    plot1 <- ggplot(scoredata, aes(reorder(Genotypes, mean), y=mean, fill = genus, pattern = rye)) +
      geom_bar(binwidth = input$bin1, stat="identity", color = "black") +
      scale_fill_manual("legend", values = c("BW" = "#F8766D", "DW" = "#00BA38", "HT" = "#619CFF")) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4,
                position=position_dodge(.9)) +
      geom_bar_pattern(stat = "identity", pattern_color = "black",
                   pattern_fill = "black", pattern_spacing = 0.015, pattern_density = 0.02) +
      scale_pattern_manual(values = c(N = "none", Y = "stripe")) +
      xlab("Genotypes") +
      ylab("score") +
      coord_flip() +
      theme_classic()
  
  
    plot2 <- ggplot(scoredata, aes(x=genus, y=mean, fill=genus)) + 
      geom_boxplot(binwidth = input$bin2)+
      labs(title="Score - Species",x="Species", y = "score") +
      theme_bw()
  
    plot3 <- ggplot(scoredata, aes(x=rye, y=mean, fill=rye)) + 
      geom_boxplot(binwidth = input$bin3)+
      labs(title="Score - Rye",x="Rye", y = "score") +
    theme_bw()
  
    grid.arrange(grid.arrange(plot1), grid.arrange(plot2, plot3, ncol = 1), ncol = 2, widths = c(2, 1))
  })
}

shinyApp(ui = ui, server = server)

