library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Scenario based comparison of estimators for the relative uncertainty about the estimate of the annual average Emission Factor (EF) of a fuel gas"),
    #p(h4("With or without use of an auxiliary (proxy) variable"))
    # Sidebar with a slider input for number of observations
  sidebarPanel(
    numericInput("mu_EF", label = HTML("&mu; <sub>Y</sub>: mean measured EF (tCO<sub>2</sub>/t)"), value = 2.55),
    numericInput("sd_EF", label = HTML("&sigma; <sub>Y</sub>: std dev of measured EF (tCO<sub>2</sub>/t)"), value = 0.1),
    numericInput("mu_B", label = HTML("&mu; <sub>B</sub>: annual mean fuel gas flow rate (t/day)"), value = 600),
    numericInput("sd_B", label = HTML("&sigma; <sub>B</sub>: standard deviation of Flows (t/day)"), value = 50),
    numericInput("n", label = HTML("n: Number of samples with measured EF, auxiliary variable and Flow"), value = 100),
    numericInput("rho", label = HTML("&rho;<sub>YX</sub>: correlation between EF and auxiliary variable"), value = 0.96, min = -1, max = 1),
    numericInput("rhoB", label = HTML("&rho;<sub>YB</sub>: correlation between EF and flow rate"), value = 0, min = -1, max = 1),
    numericInput("nsims", label = HTML("number of simulations of synthetic data"), value = 500, min = 50, max = 10000),
    checkboxInput("doBoot", label = "do Bootstrap estimates", value = FALSE),
    numericInput("nboots", label = HTML("number of bootstrap samples"), value = 1000, min = 50, max = 10000)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Illustration", 
               plotOutput("Scenarioplot", width = "850px", height = "700px")), 
      #tabPanel("Results", plotOutput("Resultplot")),
      tabPanel("Results", dataTableOutput("ResultData")),
      tabPanel("Readme", tags$iframe(style="height:800px; width:100%; scrolling=yes", 
                                     src="SamplingFuelGas_EF_V2.pdf"))
    )
  )
))