#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("~/MassTesting/buildFigures.R")



# Define UI for application that draws a histogram
ui <- navbarPage("Frequent PCR Testing for Airborne Pathogens",
   

   tabsetPanel(
     tabPanel(
       "Intervention Effectiveness",
       # Sidebar with a slider input for number of bins 
       sidebarLayout(
          sidebarPanel(
             sliderInput("testDelay",
                         "Test Delay [hours]:",
                         min = 0,
                         max = 72,
                         value = 24),
             #sliderInput("contactsPerDay","Contacts Per Day:", min = 0,max = 40,value = 13),
              sliderInput("fracTest","Fraction of population regularly testing:", min = 0,max = 0.995,value = 0.90),
             sliderInput("fracIso","Fraction of transmissions prevented by isolation:", min = 0,max = 0.995,value = 0.90),
             checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(1, 2,3,5,7,10,30), selected= (list(1,3,7)), inline = TRUE),
             sliderInput("maskEffect","Fraction of transmissions prevented by masks:", min = 0,max = 0.995,value = 0.0)
             ),
          
          # Show a plot of the generated distribution
          mainPanel(
             plotOutput("controlRegion", height="550px")
          )
       )
      ),
     tabPanel(
       "Economic Cost",   
       sidebarLayout(
         sidebarPanel(
           sliderInput("variableTestCost","Variable Cost Per Test (USD):", min = 1,max = 100,value = 10),
           sliderInput("isolationCost","Cost of supporting case isolation:", min = 0,max = 50000,value = 5000),
          
           sliderInput("fixedAnnualizedDailyTestCost","Annual Fixed Cost per Daily Test Capability:", min = 0.01,max = 10,value = 0.28)
           
         ),
         mainPanel(
           plotOutput("PrevalenceCost", height="500px")
           
         )
         
       )
       # as function of daily infections
       # as function of import rate (with targeted strategies as an option)
       # fixed cost vs maximum testing frequency
      
      ),     
     tabPanel(
        "Outbreak Model"
        # 2 compartment SIR model -inside and outside region
        # outside not controlled much
        # inside applies border control at a certain time that prevents _% of infected cases entering
        # applies mass testing at another time
        # parameterized by normal number of travelers, viral parameters, type of intervention
        
        # show cumulative cost and cumulative infections over time.
        # show fraction infected over time
        # add ability for temporary lockdown?
        
      ),
     tabPanel(
       "Model Parameters",
       sidebarLayout(
         sidebarPanel(
           sliderInput("relativeDeclineSlope","Relative Slope of Viral Decline:", min = 0.1,max = 3.0,value = 1.0),
           sliderInput("maxDaysAfterPeak","Maximum Number of days after peak \n viral load that infection ends:", min = 0,max = 20,value = 30),
           
          sliderInput("maxProbTransmitPerExposure","Maximum Probability of Transmission Per Exposure:", min = 0.1,max = 0.9,value = 0.3), # 0.3 would be consistent with 95% of measles household contacts infected -> 1 - 0.7^8 = 0.94
           
          plotOutput("Infectiousness", height="130px"),
          plotOutput("TestSensitivity", height="130px"),
          sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
          
         ),
         mainPanel(
           plotOutput("Trajectories", height="500px")
           
         )
         
       )
     )

   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$controlRegion <- renderPlot({
     
     #input$contactsPerDay
     inputParams = c(contactsPerHour = 13/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, 
                     precision = input$simPrecision, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, maskEffect = input$maskEffect)
     
     generateControllabilityFigure(24*as.numeric(input$testPeriods), inputParams)
     
     
    })
   
   output$Trajectories <- renderPlot({
     #reactive({
       #req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
        
     inputParams = c( logPeakLoad = 10, contactsPerHour = 13/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     
      plotTrajectories(inputParams) 
       
    # })
   })
   
   output$Infectiousness <- renderPlot({
     inputParams = c(contactsPerHour = 13/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     plotInfectiousness(inputParams) 
   })
   output$TestSensitivity <- renderPlot({
     inputParams = c(contactsPerHour = 13/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     plotTestSensitivity(inputParams) 
   })

   output$PrevalenceCost <- renderPlot({
     inputParams = c(variableTestCost = input$variableTestCost, isolationCost = input$isolationCost, fixedAnnualizedDailyTestCost = input$fixedAnnualizedDailyTestCost)
     plotPrevalenceCost(as.numeric(input$testPeriods), inputParams) 
   })   
}

# Run the application 
shinyApp(ui = ui, server = server)

