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
             sliderInput("maskEffect","Fraction of transmissions prevented by masks:", min = 0,max = 0.995,value = 0.0),
             h3("Calculation details:"),
             p("𝛾: Fraction of (homogeneous) population testing regularly
               𝛽: Isolation effectiveness (fraction reduction in transmissions for detected positives)
               𝜎: Fraction of counterfactual transmissions occurring after receiving a positive test result
               "),
             p("The expected daily transmissions and daily PCR sensitivity can be estimated for a viral load trajectory using the functions on the Viral Model page."),
             p("By averaging over test timing offsets and test outcomes, the fraction of counterfactual transmissions occuring after a positive test can be computed for each testing strategy and viral load trajectory."),
             
             ),
          
          # Show a plot of the generated distribution
          mainPanel(
             plotOutput("controlRegion", height="550px"),
             h2("Motivation"),
             p("Novel airborne pandemic with no available vaccines or treatments"),
             p("Global elimination unlikely, so expecting imported cases (even with strong border controls)"),
             p("\nGoal: reduce infections while waiting for vaccines"),
             p("Challenge: because of the high cost of existing NPIs, many countries may be unwilling or unable to control the epidemic"),
             p("Proposed solution: frequent saliva PCR testing for most of the population")
             
             #Todo: show peak image with peak viral load and time to peak, describe how figure generated
          )
       )
      ),
     tabPanel(
       "Viral Model",
       sidebarLayout(
         sidebarPanel(

           sliderInput("minLogPCRViralLoad","Minimum viral load (log10 copies / ml) for PCR detection :", min = 0,max = 6,value = 3.0), 
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
          
           h3("Example Viral Trajectories with R0 = 3"),
           plotOutput("Trajectories", height="500px")
           #Todo: Add figures for fraction of transmissions occuring after positive test for each trajectory
           
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
           # todo: add as function of import rate (with targeted strategies as an option)
           # Todo: add sources for fixed and variable costs
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
           title = "speed for pathogens with short generation time",
           collapsed = TRUE,
           "- sol: reduce pcr delay time - sol: add household quarantine or some contact tracing"
         ),
         accordionItem(
           title = "enforcing isolation and testing",
           collapsed = TRUE,
           "- sol: require testing for daily activities
         - sol: validate isolation being followed, and generously support it
           - response: for many spillover pathogens, R0 is initially fairly low, so even partial adherence may be enough
           - response: for higher R0 pathogens, partially effective mass testing substitutes for more expensive social distancing"
         ),
         accordionItem(
           title = "early pandemic time to functionality (materials, personnel, equipment, logistics)",
           collapsed = TRUE,
           "sol: more analysis of failure modes
         - sol: run peacetime simulation excercises"
         ),
         accordionItem(
           title = "accuracy of viral load/infectiosness/pcr model for novel pathogens",
           collapsed = TRUE,
           "argument: many common
         - argument: PCR is fairly well understood -if the RNA is there at a certain concentration it can be found"
         ),
         accordionItem(
           title = "running cost too high",
           collapsed = TRUE,
           "sol: reduce testing frequency in areas with low chance of infections
         - sol: automate further, develop own reagents"
         ),
         accordionItem(
           title = "population behaviour in different pandemic scenarios",
           collapsed = TRUE,
           "sol for seemingly low IFR: reduce cost and burden of response
         - sol for high IFR: keep prevalence very low, have sufficient PPE"
         ),
         accordionItem(
           title = "heterogeneous population behavioiur",
           collapsed = TRUE,
           "sol: feedback mechanisms to find and respond to process failures"
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
     
     # tabPanel(
     #    "Outbreak Model"
     #    # 2 compartment SIR model -inside and outside region
     #    # outside not controlled much
     #    # inside applies border control at a certain time that prevents _% of infected cases entering
     #    # applies mass testing at another time
     #    # parameterized by normal number of travelers, viral parameters, type of intervention
     #    
     #    # show cumulative cost and cumulative infections over time.
     #    # show fraction infected over time
     #    # add ability for temporary lockdown?
     #    
     #  ),


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
   

}

# Run the application 
shinyApp(ui = ui, server = server)

