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
              sliderInput("fracTest","Fraction of population regularly testing:", min = 0,max = 0.995,value = 0.95),
             sliderInput("fracIso","Fraction of transmissions prevented by isolation:", min = 0,max = 0.995,value = 0.95),
             checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(1, 2,3,5,7,10,30), selected= (list(1,3,7)), inline = TRUE)
             ),
          
          # Show a plot of the generated distribution
          mainPanel(
             plotOutput("controlRegion")
          )
       )
      ),
     tabPanel(
       "Economic Cost"
      
      ),
     tabPanel(
       "Model Parameters",
       sidebarLayout(
         sidebarPanel(
           sliderInput("relativeDeclineSlope","Relative Slope of Viral Decline:", min = 0.1,max = 3.0,value = 1.0),
           sliderInput("maxDaysAfterPeak","Maximum Number of days after peak \n viral load that infection ends:", min = 0,max = 20,value = 30),
           
           plotOutput("viralTrajectory", height="130px"),
           sliderInput("maxProbTransmitPerExposure","Maximum Probability of Transmission Per Exposure:", min = 0.1,max = 0.9,value = 0.2), # 0.3 would be consistent with 95% of measles household contacts infected -> 1 - 0.7^8 = 0.94
           sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
         ),
         mainPanel(
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
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     
     generateControllabilityFigure(24*as.numeric(input$testPeriods), inputParams)
     
     
    })
   
   output$viralTrajectory <- renderPlot({
     #reactive({
       #req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
        
     inputParams = c(relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     
      plotViralLoads(inputParams) 
       
    # })
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

