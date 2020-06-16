# An Rshiny app to display local Coronavirus cases
# Doug McNeall dougmcneall@gmail.com @dougmcneall

# Needs fixing with this lesson from the tutorial

#https://shiny.rstudio.com/tutorial/written-tutorial/lesson6/

library(shiny)
library(viridis)


ui <- fluidPage(
  titlePanel("Local area lab-confirmed Coronovirus cases"),
  
  sidebarLayout(
    
    sidebarPanel( textInput("area", "Enter Area. Try 'Devon' or 'Exeter'")),
    
    mainPanel( plotOutput("coolplot")
                  
    )
  )
)


server <- function(input, output) {
  
  getLocalDat <- reactive({
    
    dat <- read.csv('https://coronavirus.data.gov.uk/downloads/csv/coronavirus-cases_latest.csv')
    localDat <- dat[dat[,'Area.name']==input$area, ]
  })
  
  
  output$coolplot <- renderPlot({
      # 'area' should be a string that names an area of interest
    localData <- getLocalDat()
      
      #par(las = 1, mar = c(2.5,4,4,1))
    plot(as.Date(localData[, 'Specimen.date']), localData[, 'Daily.lab.confirmed.cases'],
           
    xlab = '', ylab = 'cases', main = paste0('Daily lab-confirmed Coronavirus cases in ', input$area),
           type = 'h', lwd = 5, col = 'skyblue2', bty = 'n')
      grid()
    
  })
  
  
}

shinyApp(ui = ui, server = server)