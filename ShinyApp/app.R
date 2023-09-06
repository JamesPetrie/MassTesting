#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)


folder = "/Users/orayasrim/Documents/MassTest/MassTesting/" # '/Users/jpetrie/MassTesting/' #
source(paste0(folder, "ShinyApp/buildFigures.R"))



# Define UI for application that draws a histogram
ui <- navbarPage("Frequent PCR Testing for Airborne Pathogens. Made by James Petrie (jimpetrie@uwaterloo.ca)",
   

   tabsetPanel(
     tabPanel(
       "Intervention Effectiveness",
       # Sidebar with a slider input for number of bins 
       sidebarLayout(
          sidebarPanel(
            h3("Effectiveness Calculation:"),
            withMathJax(p('The graph is created by finding the largest \\(R_0\\) (modified by varying peak viral load) for each testing strategy such that \\(R_e \\leq 1\\). Pathogens below each line can be controlled by that testing strategy.')),
            withMathJax(p('\\(R_e\\) is the effective reproduction number when using a mass-testing strategy. If \\(R_e < 1\\) then the number of infected cases will decrease.')),
            withMathJax(p('$$R_e = R_0 \\cdot (1 - \\gamma  \\cdot \\beta \\cdot \\sigma) \\cdot (1-\\lambda)$$')),
            
                     #sliderInput("contactsPerDay","Contacts Per Day:", min = 0,max = 40,value = 13),
            sliderInput("fracTest",withMathJax(p('\\(\\gamma\\): Fraction of (homogeneous) population testing regularly:')), min = 0,max = 0.995,value = 0.90),
            sliderInput("fracIso",withMathJax(p('\\(\\beta\\):  Isolation effectiveness (fraction reduction in transmissions for detected positives):')), min = 0,max = 0.995,value = 0.90),
            withMathJax(p('\\(\\sigma\\):   Fraction of counterfactual transmissions occurring after receiving a positive test result. The calculation is shown in the Model tab (depending on viral load trajectory, test frequency, and test delay)')),
            sliderInput("testDelay", "Test Delay [hours]:",min = 0, max = 72,value = 12),
            checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(1, 2,3,5,7,10,30), selected= (list(1,3)), inline = TRUE),
            sliderInput("maskEffect",withMathJax(p('\\(\\lambda\\):  Fraction of transmissions prevented by masks:')), min = 0,max = 0.995,value = 0.0),

  
             ),
          
          # Show a plot of the generated distribution
          mainPanel(
            h2("Motivation"),
            HTML("<b>Scenario</b>: Novel airborne pandemic with no available vaccines or treatments. Global elimination unlikely so expecting imported cases<br/>"),
            HTML("<b>Goal</b>: Reduce the number of infections while waiting for vaccines<br/>"),
            HTML("<b>Challenge</b>: Because of the high cost of existing Non-Pharmaceutical Interventions (NPIs), many countries may be unwilling or unable to control the epidemic<br/>"),
            HTML("<b>Proposed solution</b>: Frequent (inexpensive and convenient) saliva PCR testing for most of the population. Generous financial support for isolation of people who test positive.<br/><br/><br/>"),
             plotOutput("controlRegion", height="550px"),

             
             #Todo: show peak image with peak viral load and time to peak, describe how figure generated
          )
       )
      ),
     tabPanel(
       "Model",
       sidebarLayout(
         sidebarPanel(

           sliderInput("logLimitOfDetection","Minimum viral load (log10 copies / ml) for PCR detection :", min = 0.0,max = 6.0,value = 3.0, step = 0.5), 
           plotOutput("TestSensitivity", height="130px"),
           
           sliderInput("maxProbTransmitPerExposure","Maximum Probability of Transmission Per Exposure:", min = 0.1,max = 0.9,value = 0.3), # 0.3 would be consistent with 95% of measles household contacts infected -> 1 - 0.7^8 = 0.94
           sliderInput("contactsPerDay","Contacts per day:", min = 1,max = 50,value = 13), 
           #test
           sliderInput("probTransmitMid","midpoint test infectiousness :", min = 10e3,max = 10e9,value = 10e4, step = 10e1), 
           plotOutput("Infectiousness", height="140px"),
           sliderInput("relativeDeclineSlope","Relative Slope of Viral Decline:", min = 0.1,max = 3.0,value = 1.0),
           sliderInput("maxDaysAfterPeak","Maximum Number of days after peak \n viral load that infection ends:", min = 0,max = 20,value = 30),
           sliderInput("initialLogLoad","Viral load at time of infection (log10 copies / ml):", min = -4, max = 0,value = -2.5)
           
  # todo: computation of expected transmissions after postive test

           # sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
           
         ),
         mainPanel(
           # todo: figure of characteristic curve
           p("The infectiousness and test sensitivity for a pathogen over the course of infection depend on the viral load trajectory. A viral load trajectory can be characterized by the peak viral load and the time taken to reach the peak."),
          
           h3("Example Viral Load Trajectories with \\(R_0=4.5\\)"),
           plotOutput("Trajectories", height="500px"),
           h3("Fraction of transmissions after a positive test"),
           #Todo: Add figures for fraction of transmissions occuring after positive test for each trajectory
           p("\\(\\sigma\\), the fraction of counterfactual transmissions occuring after receiving a positive test  can be computed for each testing strategy and viral load trajectory."),
           plotOutput("FracAfterPositive", height="300px"),
           h3("Calculation"),
           p("Let \\(E[T(x)|\\pi]\\) be the expected number of transmissions on day x after infection, conditional on parameters \\(\\pi\\) describing the viral load trajectory.
             Let \\(P(Test(x) | \\pi)\\) be the probability that a test taken on day \\(x\\) is positive, also conditional on the viral load trajectory. 
             Define the function ProbAllNegative(x, offset, period) as the probability that all samples collected on or before day x are negative, 
             with the time between sequential tests set by the period variable, and the timing relative to infection set by the offset variable."),
           withMathJax(p('$$ProbAllNegative(x, offset, period) = \\prod_{k=0}^{(x-offset)/period}(1-P(Test(k*period + offset) | \\pi))$$')),
           p('The expected fraction of transmissions after a positive test is computed by averaging over test timing offset, with offset ~ unif(0,period):'),
           withMathJax(p('$$\\sigma(testDelay, testPeriod) = \\frac{E_{offset}[\\sum_{j=0}^{\\infty} E[T(j)|\\pi] \\cdot (1 - ProbAllNegative(j-testDelay, offset, testPeriod))]}{\\sum_{j=0}^{\\infty} E[T(j)|\\pi]}$$')),

           
           
           
         )
         
       )
     ),

    tabPanel(
      "Outbreak Response",
      # Sidebar with a slider input for number of bins
      sidebarLayout(
        sidebarPanel(
          sliderInput("normalTestPeriod", "Normal Test Period [days]:",min = 1, max = 30,value = 4),
          sliderInput("outbreakTestPeriod", "Outbreak Test Period [days]:",min = 1, max = 10,value = 1),
          sliderInput("tracingDelay", "Contact Tracing Delay [hours]:",min = 0, max = 48,value = 8),

          sliderInput("fractionTraced", "Fraction of contacts traced (for detected infections):",min = 0, max = 1,value = 0.3),
          sliderInput("probDetectSymptoms", "Probability of detecting symptoms and getting tested:",min = 0, max = 1,value = 0.5),

          sliderInput("timeFromPeaktoSymptoms", "Time from peak viral load to symptoms [hours]:",min = -72, max = 72,value = 0),
          sliderInput("daysToPeak", "Time from infection to peak viral load [days]:",min = 1, max = 15,value = 6),
          sliderInput("outbreakR0", "R0:",min = 1.0, max = 16,value = 2),


          #
        ),

        mainPanel(
          plotOutput("Outbreak", height="500px")
        )
      )
    ),

     # tabPanel(
     #   "Economic Cost",   
     #   sidebarLayout(
     #     sidebarPanel(
     #       sliderInput("variableTestCost","Variable Cost Per Test (USD):", min = 1,max = 100,value = 10),
     #       sliderInput("isolationCost","Cost of supporting case isolation:", min = 0,max = 50000,value = 5000),
     #      
     #       sliderInput("fixedAnnualizedDailyTestCost","Annual Fixed Cost per Daily Test Capability:", min = 0.01,max = 10,value = 0.28),
     #       wellPanel(
     #         helpText(   
     #                     a("SalivaDirect Protocol with $1.21 per sample in reagents. ",     href="https://doi.org/10.1016/j.medj.2020.12.010", target="_blank"),
     #                     HTML("Reduce cost of logistics and labour to below $1 using unstaffed booths with regular sample collection by a scalable service like Uber or Amazon delivery in combination with highly automated PCR labs.")
     #         )
     #       )
     #     ),
     #     mainPanel(
     #       plotOutput("PrevalenceCost", height="500px")
     # 
     #       # todo: add as function of import rate (with targeted strategies as an option)
     #       # Todo: add sources for fixed and variable costs
     #     )
     #     
     #   )
     # 
     #  )


   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #output$exampleImplementation <- 
  
   output$controlRegion <- renderPlot({
     
     #input$contactsPerDay
     inputParams = c(contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, probTransmitMid = input$probTransmitMid,
                     precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                     maskEffect = input$maskEffect, logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad)
     
     generateControllabilityFigure(24*as.numeric(input$testPeriods), inputParams)
     
     
    })
   
   output$Trajectories <- renderPlot({
     #reactive({
       #req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
        
     inputParams = c( logPeakLoad = 10, contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = input$probTransmitMid,
                      relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                      logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15)
     
      plotTrajectories(inputParams) 
       
    # })
   })
   
   output$FracAfterPositive <- renderPlot({
     inputParams = c(  contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = input$probTransmitMid,
                      relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                      logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15)
     
     plotFracTransmissionsAfterPositive(24*as.numeric(input$testPeriods), inputParams)
   })
   
   # output$Infectiousness <- renderPlot({
   #   inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, 
   #                   relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     
     output$Infectiousness <- renderPlot({
       inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = input$probTransmitMid,
                       relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     plotInfectiousness(inputParams) 
   })
   output$TestSensitivity <- renderPlot({
     inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, logLimitOfDetection = input$logLimitOfDetection, probTransmitMid = input$probTransmitMid)
     plotTestSensitivity(inputParams) 
   })

   output$PrevalenceCost <- renderPlot({
     inputParams = c(variableTestCost = input$variableTestCost, isolationCost = input$isolationCost, fixedAnnualizedDailyTestCost = input$fixedAnnualizedDailyTestCost)
     plotPrevalenceCost(as.numeric(input$testPeriods), inputParams) 
   })   
   
  output$Outbreak <- renderPlot({

     # todo: make some input params
     # todo: take R0 as input and compute peak viral load
     inputParams = c(

                     normalTestPeriod = input$normalTestPeriod*24  ,
                     outbreakTestPeriod = input$outbreakTestPeriod*24 ,
                     ContactTracingDelay = input$tracingDelay ,
                     ProbTracedGivenInfectorDetected = input$fractionTraced,
                     ProbDetectSymptoms = input$probDetectSymptoms,
                     timeFromPeakToSymptoms = input$timeFromPeaktoSymptoms,
                     timeToPeak = input$daysToPeak*24,
                     maxTimeAfterPeak= 24*30,
                     probTransmitMid = input$probTransmitMid,
                     logLimitOfDetection = input$logLimitOfDetection, 

                     timeFromPeakTo0 = 24*5,
       contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest,
                     precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak,
                     maskEffect = input$maskEffect, minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad)


     plotOutbreaks(numOutbreaks = 120, endDay = 120, maxSize = 300, R0 = input$outbreakR0, inputParams)
   })

}

# Run the application 
shinyApp(ui = ui, server = server)

