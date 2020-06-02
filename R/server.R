library(shiny)
source("helpers.R")
#K <- 24*365*4 # every 15 minutes
ColExtraInputs <- list(K = 24*365, # every hour
                       mu_aux = 2,
                       sd_aux = 0.2)

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
  
  ColInputs <- reactive({ list(nsims = input$nsims,
                               n = input$n,
                               mu_EF = input$mu_EF,
                               sd_EF = input$sd_EF,
                               mu_B = input$mu_B,
                               sd_B = input$sd_B,
                               rho = input$rho,
                               rho_B = input$rhoB) })

  SimDat <- reactive({ SimulateData(Sinput=ColInputs(),SExtraPars=ColExtraInputs) })

  output$Scenarioplot <- renderPlot({
    FitAuxModel(FSimDat=SimDat(),Finput=ColInputs(),doGraph=TRUE,returnAuxLm=FALSE)
    })
  
  SimResults <- reactive({ 
    SimulateResults(Zinput=ColInputs(),ZExtraPars=ColExtraInputs,doBoot=input$doBoot,nboots=input$nboots) 
    })
  
#  output$Resultplot <- renderPlot({
#    plotResults(VSimResults=SimResults(),plotCoverage=TRUE,doBoot=input$doBoot,Vinputs=ColInputs())
#  })
  output$ResultData <- renderDataTable({ 
    TableResults(TSimRes=SimResults(), doBoot=input$doBoot, Tinputs=ColInputs())
    })

})