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
library(shinyWidgets)

source("buildFigures.R")



# Define UI for application that draws a histogram
ui <- navbarPage("Frequent PCR Testing for Airborne Pathogens. Made by James Petrie (jimpetrie@uwaterloo.ca)",
   

   tabsetPanel(
     tabPanel(
       "Intervention Effectiveness",
       # Sidebar with a slider input for number of bins 
       sidebarLayout(
          sidebarPanel(
            h3("Effectiveness Calculation:"),
            withMathJax(p('$$R_e = R_0 \\cdot (1 - \\gamma  \\cdot \\beta \\cdot \\sigma) \\cdot (1-\\lambda)$$')),
            
            withMathJax(p('\\(R_e\\) is the effective reproduction number when using a mass-testing strategy. If \\(R_e < 1\\) then the number of infected cases will decrease.')),
            
            withMathJax(p('The graph is created by finding the largest \\(R_0\\) (modified by varying peak viral load) for each testing strategy such that \\(R_e \\leq 1\\). Pathogens below each line can be controlled by that testing strategy.')),
                  #sliderInput("contactsPerDay","Contacts Per Day:", min = 0,max = 40,value = 13),
            sliderInput("fracTest",withMathJax(p('\\(\\gamma\\): Fraction of (homogeneous) population testing regularly:')), min = 0,max = 0.995,value = 0.90),
            sliderInput("fracIso",withMathJax(p('\\(\\beta\\):  Isolation effectiveness (fraction reduction in transmissions for detected positives, conditional on being a person who adheres to testing):')), min = 0,max = 0.995,value = 0.90),
            withMathJax(p('\\(\\sigma\\):   Fraction of counterfactual transmissions occurring after receiving a positive test result. The calculation is shown in the Model tab (depending on viral load trajectory, test frequency, and test delay)')),
            sliderInput("testDelay", "Test Delay [hours]:",min = 0, max = 72,value = 12),
            checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(1, 2,3,5,7,10,30), selected= (list(1,3,7)), inline = TRUE),
            sliderInput("maskEffect",withMathJax(p('\\(\\lambda\\):  Fraction of transmissions prevented by masks:')), min = 0,max = 0.995,value = 0.0),

  
             ),
          
          # Show a plot of the generated distribution
          mainPanel(
            h2("Motivation"),
            HTML("<b>Scenario</b>: Novel airborne pandemic with no available vaccines or treatments. Global elimination unlikely so expecting imported cases<br/>"),
            HTML("<b>Goal</b>: Reduce the number of infections while waiting for vaccines<br/>"),
            HTML("<b>Challenge</b>: Because of the high cost of existing Non-Pharmaceutical Interventions (NPIs), many countries may be unwilling or unable to control the epidemic<br/>"),
            HTML("<b>Proposed solution</b>: Frequent (inexpensive and convenient) saliva PCR testing for most of the population. Generous financial support for isolation of people who test positive.<br/><br/><br/>"),
             plotOutput("controlRegion", height="550px")

             
             #Todo: show peak image with peak viral load and time to peak, describe how figure generated
          )
       )
      ),
     tabPanel(
       "Model",
       sidebarLayout(
         sidebarPanel(

           sliderInput("minLogPCRViralLoad","Minimum viral load (log10 copies / ml) for PCR detection :", min = 0.0,max = 6.0,value = 3.0, step = 0.5), 
           plotOutput("TestSensitivity", height="130px"),
           
           sliderInput("maxProbTransmitPerExposure","Maximum Probability of Transmission Per Exposure:", min = 0.1,max = 0.9,value = 0.3), # 0.3 would be consistent with 95% of measles household contacts infected -> 1 - 0.7^8 = 0.94
           sliderInput("contactsPerDay","Contacts per day:", min = 1,max = 50,value = 13), 
           
           plotOutput("Infectiousness", height="140px"),
           sliderInput("relativeDeclineSlope","Relative Slope of Viral Decline:", min = 0.1,max = 3.0,value = 1.0),
           sliderInput("maxDaysAfterPeak","Maximum Number of days after peak \n viral load that infection ends:", min = 0,max = 20,value = 30),
           sliderInput("initialLogLoad","Viral load at time of infection (log10 copies / ml):", min = -4, max = 0,value = -2.0)
           
  # todo: computation of expected transmissions after postive test

           # sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
           
         ),
         mainPanel(
           # todo: figure of characteristic curve
           p("The infectiousness and test sensitivity for a pathogen over the course of infection depend on the viral load trajectory. A viral load trajectory can be characterized by the peak viral load and the time taken to reach the peak."),
          
           h3("Example Viral Load Trajectories with \\(R_0=3\\)"),
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
       "Economic Cost",   
       sidebarLayout(
         sidebarPanel(
           sliderTextInput("variableTestCost","Variable Cost Per Test (USD):",
                                         choices=c(0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0),
                                         selected=2.0, grid = T),
           #sliderInput("variableTestCost","Variable Cost Per Test (USD):", min = 1,max = 100,value = 10),
           sliderInput("isolationCost","Cost of supporting case isolation:", min = 0,max = 50000,value = 5000),
          
           sliderInput("fixedAnnualizedDailyTestCost","Annual Fixed Cost per Daily Test Capability:", min = 0.01,max = 10,value = 0.28)
           
         ),
         mainPanel(
           plotOutput("PrevalenceCost", height="500px")
           # todo: add as function of import rate (with targeted strategies as an option)
           # Todo: add sources for fixed and variable costs
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
        #sliderInput("tracingDelay", "Contact Tracing Delay [hours]:",min = 0, max = 48,value = 8),
        
        #sliderInput("fractionTraced", "Fraction of contacts traced (for detected infections):",min = 0, max = 1,value = 0.3),
        #sliderInput("probDetectSymptoms", "Probability of detecting symptoms:",min = 0, max = 1,value = 0.5),
        
        #sliderInput("timeFromPeaktoSymptoms", "Time from peak viral load to symptoms [hours]:",min = -72, max = 72,value = 0),
        sliderInput("daysToPeak", "Time from infection to peak viral load [days]:",min = 1, max = 15,value = 6),
        sliderInput("outbreakR0", "R0:",min = 1.0, max = 16.0,value = 2.0, round = FALSE, step = 0.25),
        
        
        #
      ),

      mainPanel(
        plotOutput("Outbreak", height="500px")
      )
    )
  ),
     tabPanel(
       "Example Implementation",

         includeMarkdown("~/MassTesting/implementation.Rmd")



     ),

     tabPanel(

       "Potential Challenges",
       accordion(
         id = "accordion1",
         accordionItem(
           title = "Testing and isolation adherence",
           collapsed = TRUE,
           HTML("- sol: require proof of a recent negative test for daily activities <br/>
                - sol: generously support isolation and verify adherence <br/>
                - response: for many spillover pathogens, R0 is initially fairly low, so even partial adherence may be enough <br/>
                - response: for higher R0 pathogens, partially effective mass testing substitutes for more expensive social distancing<br/>")
         ),
         accordionItem(
           title = "speed for pathogens with short generation time",
           collapsed = TRUE,
           HTML("- sol: reduce pcr delay time <br/>
                - sol: add household quarantine or some contact tracing")
           ),
         accordionItem(
           title = "early pandemic time to functionality (materials, personnel, equipment, logistics)",
           collapsed = TRUE,
           HTML("sol: more analysis of failure modes<br/>
         - sol: run peacetime simulation excercises")
         ),
         accordionItem(
           title = "accuracy of viral load/infectiosness/pcr model for novel pathogens",
           collapsed = TRUE,
           HTML("argument: many common<br/>
         - argument: PCR is fairly well understood -if the RNA is there at a certain concentration it can be found"
                )
         ),
         accordionItem(
           title = "running cost too high",
           collapsed = TRUE,
           HTML("sol: reduce testing frequency in areas with low chance of infections<br/>
         - sol: automate further, develop own reagents")
         ),
         accordionItem(
           title = "population behaviour in different pandemic scenarios",
           collapsed = TRUE,
           HTML("sol for seemingly low IFR: reduce cost and burden of response<br/>
         - sol for high IFR: keep prevalence very low, have sufficient PPE")
         ),
         accordionItem(
           title = "heterogeneous population behavioiur",
           collapsed = TRUE,
           HTML("sol: feedback mechanisms to find and respond to process failures")
         ),     
         accordionItem(
           title = "perverse incentives",
           collapsed = TRUE,
           HTML("sol: calibrate incentives to not substantially exceed perceived cost of infection and isolation")
         )
         
       )

     ),

     tabPanel("Related Work",
              p("Mass testing deployments in western countries have been limited in time or the fraction
of the population participating (e.g. Slovakia testing the whole population twice, but not on
                an ongoing basis [188]). The expense of massively scaling testing was likely a contributing
                factor preventing more widespread deployment. The limited scope of deployment limited
                the possible success because people testing regularly could still be infected by members of
                the broader community, and cases could rebound after time-limited interventions.

                For the original COVID-19 variant, Elbanna and Goldenfeld estimated that testing twice a week would be sufficient to bring Rt below 1 (assuming perfect isolation of those who test positive)

                Mass testing seems to be an important component of China’s “dynamic Covid-zero”
policy (as of July 2022). The independent effect of mass testing is difficult to estimate
                because they are using it in combination with a few other strategies like venue based trac-
                ing and localized lockdowns. The combined approach seems to be effective at containing
                outbreaks, but it is difficult to determine what the cost of this is

                ")

              ),
     tabPanel("Policy Proposal",
              h3("1 Build enough PCR equipment in advance to be able to test every person every day during a pandemic"),
              p(""),
              h3("2 Incorporate mass testing into national pandemic response plans"),
              h3("3 Regularly perform simulation excercises to test system readiness")
              )




   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #output$exampleImplementation <- 
  
   output$controlRegion <- renderPlot({
     
     #input$contactsPerDay
     inputParams = c(contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, 
                     precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                     maskEffect = input$maskEffect, minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad)
     
     generateControllabilityFigure(24*as.numeric(input$testPeriods), inputParams)
     
     
    })
   
   output$Trajectories <- renderPlot({
     #reactive({
       #req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
        
     inputParams = c( logPeakLoad = 10, contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, 
                      relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                      minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad, precision = 0.15)
     
      plotTrajectories(inputParams) 
       
    # })
   })
   
   output$FracAfterPositive <- renderPlot({
     inputParams = c(  contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, 
                      relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                      minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad, precision = 0.15)
     
     plotFracTransmissionsAfterPositive(24*as.numeric(input$testPeriods), inputParams)
   })
   
   output$Infectiousness <- renderPlot({
     inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, 
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     plotInfectiousness(inputParams) 
   })
   output$TestSensitivity <- renderPlot({
     inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, minLogPCRViralLoad = input$minLogPCRViralLoad)
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
                     ContactTracingDelay = 0, #input$tracingDelay ,
                     ProbTracedGivenInfectorDetected = 0, #input$fractionTraced,
                     ProbDetectSymptoms = 0, #input$probDetectSymptoms,
                     timeFromPeakToSymptoms = 0,#input$timeFromPeaktoSymptoms,
                     timeToPeak = input$daysToPeak*24,
                     maxTimeAfterPeak= 24*30, 
                     
                     
                     timeFromPeakTo0 = 24*5,
       contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracQuar = input$fracIso, fracTest = input$fracTest, 
                     precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                     maskEffect = input$maskEffect, minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad)
     
     
     plotOutbreaks(numOutbreaks = 1000, endDay = 120, maxSize = 100, R0 = input$outbreakR0, inputParams) 
   })   

}

# Run the application 
shinyApp(ui = ui, server = server)

