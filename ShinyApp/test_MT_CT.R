
#test file for MT & CT
folder =  "~/MassTesting/ShinyApp/" # "/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/" #
Rcpp::sourceCpp(paste0(folder, "ViralLoad.cpp"))
source(paste0(folder, "buildFigures.R"))
library(boot)
#plots (daily testing) another will have testing every 4 days for example ( noraml test period and outbreak test period)
#default 1 day and additional testing with 0 or 2 days ( tracing delay)
#fraction traced along x axis () ( probtraced given )
#prob detect symptoms (depends on disease) -> can set to 0 to ignore -> modify based on disease
#time to peak to symptoms -> modify based on disease
#days to peak -> depends on disease
#contact per day is 13, test delay = 10 hours, fraciso = 100 or 1, frac test = 0.9 to do thihnk about turn off contact tracing for people who dont test, 
#maxprobtranmit per exposure = 0.3, relative decline slope = 1, max days after peak. = 20 , mask effect = 0, minlogpcr = +2.5, initial log load = -2.5

# #base
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0 #adjust based on disease
# timeFromPeaktoSymptoms = 1 #adjust based on disease
# daysToPeak = 6 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 3 # adjust based on disease 

#example diseases
# #influenza
#   probTest = .50 #adjust scenario
#   normalTestPeriod = 4 #adjust scenario
#   outbreakTestPeriod = 4 #adjust scenario
#   tracingDelay = 0 #adjust scenario
#   #fractionTraced = 0.1 #adjust this one
#   probDetectSymptoms = 0.693*probTest #adjust based on disease
#   timeFromPeaktoSymptoms = -2 #adjust based on disease
#   daysToPeak = 3 #adjust based on disease
#   contactsPerDay = 13
#   testDelay = 10
#   fracIso = 1
#   fracTest = 0.90
#   maxProbTransmitPerExposure = 0.30
#   relativeDeclineSlope = 1
#   maxDaysAfterPeak = 20
#   maskEffect = 0
#   minLogPCRViralLoad = 2.5
#   initialLogLoad = -2.5
#   R0 = 2.0 #adjust based on disease

# #SARS-CoV
# probTest = .80 #adjust scenario
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0.925*probTest #adjust based on disease
# timeFromPeaktoSymptoms = -5 #adjust based on disease
# daysToPeak = 7 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 2.4 # adjust based on disease

#add covid-19 

minLogPCRViralLoad = 2.5
initialLogLoad = -2.5
maxProbTransmitPerExposure = 0.30
relativeDeclineSlope = 1
contactsPerDay = 13
precision = 2.0



meanfun <- function(data, i){
  d <- data[i, ]
  return(mean(d))   
}


test_period_list <- c(1,2,4,8,16)

disease_list <- c("1918 Influenza", "SARS-CoV-1","SARS-CoV-2")
#disease_list <- c("test")

fractionTraced_list <- seq(0, 0.9, by = 0.1)

timepeak_list <- seq(1,12,2)
prob_test_list <- c(0.10, 0.50, 0.80, 0.90) 
frac_iso_list <- c(.10,.50,.90) #remove 100% for figure B4

runAndComputeRe = function(numOutbreaks, endDay, maxSize, diseaseName, R0, timeToPeak, timeFromPeaktoSymptoms,  probTestSymptoms, 
                           testPeriod, tracingDelay, fractionTraced, testDelay, fracIso, fracTest ){
  params = c(
    normalTestPeriod = testPeriod,
    outbreakTestPeriod = testPeriod ,
    ContactTracingDelay = tracingDelay ,
    ProbTracedGivenInfectorDetected = fractionTraced,
    ProbTestSymptoms = probTestSymptoms,
    timeFromPeakToSymptoms = timeFromPeaktoSymptoms,
    timeToPeak = timeToPeak,
    timeFromPeakTo0 = timeToPeak, 
    #timeFromPeakTo0 = 24*5,
    maxTimeAfterPeak= 24*30,
    contactsPerHour = contactsPerDay/24, testDelay = testDelay, 
    fracIso = fracIso, fracQuar = fracIso, fracTest = fracTest, 
    maxProbTransmitPerExposure = maxProbTransmitPerExposure,
    #maxTimeAfterPeak = 24*30, 
    maxTimeAfterPeak = 24*20,
    maskEffect = 0, minLogPCRViralLoad = minLogPCRViralLoad, initialLogLoad = initialLogLoad, 
    probTransmitMid = 8.9e6,logLimitOfDetection = minLogPCRViralLoad ,infectHParam = 0.51, 
    #precision = precision
    precision = 0.2
    )
  
    
  caseData = rbindlist(llply(1:numOutbreaks, function(i){
    #browser()
    dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
    dt[,RunNumber := i]
    
    return(dt)
  }))
 
  caseData[,hourNotInfectious := InfectedHour + 5*24*8] # todo: fix
  re_filt <- caseData[hourNotInfectious<SimulationEndHour]
  #for bootstrapping
  test <- data.frame(re_filt$NumInfected)
  bo <- boot(test[, "re_filt.NumInfected", drop = FALSE], statistic=meanfun, R=1000)
  
  result <- boot.ci(bo, conf=0.95, type = c("basic"))
  
  lowerBound <- result$basic[1,4]
  
  upperBound <- result$basic[1,5]
  #end bootstrapping
  
  re <- mean(re_filt$NumInfected)
  
  x = data.table(Re = re, Disease = diseaseName, FractionTraced = fractionTraced, R0 = R0, testPeriod = testPeriod,timePeak = timeToPeak, minBound = lowerBound, maxBound = upperBound,fractionIso = fracIso)
  #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, probTestPositive = prob_test)
  #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, fractionIso = frac_iso)
  
  return(x)
}


runAndComputeRe(numOutbreaks = 30, endDay = 100, maxSize = 300, diseaseName = "Example", R0 = 3, timeToPeak = 24*5, timeFromPeaktoSymptoms = 0,  probTestSymptoms = 1.0, 
                testPeriod = 100000*24, tracingDelay = 0, fractionTraced = 1, testDelay = 500, fracIso = 0.9, fracTest = 0)

