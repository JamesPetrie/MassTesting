#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("~/MassTesting/ViralLoad.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Frequent Testing"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("testDelay",
                     "Test Delay [hours]:",
                     min = 0,
                     max = 72,
                     value = 24),
         #sliderInput("contactsPerDay","Contacts Per Day:", min = 0,max = 40,value = 13),
         sliderInput("fracIso","Fraction of transmissions prevented by isolation:", min = 0,max = 0.995,value = 0.95),
         sliderInput("fracTest","Fraction of population regularly testing:", min = 0,max = 0.995,value = 0.95),
         checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(0.5, 1, 2,3,4,5,7,10,30), selected= (list(1,3,10))),
         sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("controlRegion")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$controlRegion <- renderPlot({
     
     #input$contactsPerDay
     inputParams = c(contactsPerHour = 13/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, precision = input$simPrecision)
     
     print(input$testPeriods)
     dt = rbindlist(llply(24*as.numeric(input$testPeriods), function(testPeriod){
       params = copy(inputParams)
       params["testPeriod"] = testPeriod
       evaluateStrategy(params)
     }))
     if(nrow(dt)>0){
       freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24)), by = TestPeriod]
       setkey(freqNames, by = "TestPeriod")
       freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
       dt = merge(dt, freqNames, by = "TestPeriod")
       ggplot(dt, aes(x = TimeToPeak/24, y = MaxR0, colour = FreqLabel)) +
         geom_line(linewidth = 1.4) + scale_y_continuous(breaks = 0:15, limits = c(0,10)) + scale_x_continuous(breaks = 0:10) + 
         guides(colour=guide_legend(title="Test Frequency \n [1/Days]")) + xlab("Time to Peak Viral Load [Days]") + ylab("Maximum Controllable R0") 
       
         #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
     }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

