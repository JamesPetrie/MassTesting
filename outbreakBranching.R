
require(data.table)
require(ggplot2)
require(Rcpp)
require(plyr)
require(tidyr)

Rcpp::sourceCpp("~/MassTesting/viralLoad.cpp")



params = c(
           ProbDetectSymptoms = 0.5,
           ProbTracedGivenInfectorDetected = 0.5,
           RelativeTransmissionRisk_Detected = 0.05,
           ContactTracingDelay = 3,
           
           
           # viral load params
           testPeriod = 72 ,
           contactsPerHour = 13/24, 
           testDelay = 8, 
           fracIso = 0.9, 
           fracTest = 0.9, 
           maxProbTransmitPerExposure = 0.3, 
           relativeDeclineSlope = 1.0, 
           maxTimeAfterPeak= 24*30, 
           logPeakLoad = 8, 
           timeToPeak = 24*6,
           timeFromPeakToSymptoms = 0,
           timeFromPeakTo0 = 24*6,
           initialLogLoad = -2, 
           minLogPCRViralLoad = 3)

simulateOutbreaks = function(numOutbreaks = 10, endDay = 60, maxSize = 100, params){
  caseData = rbindlist(llply(1:numOutbreaks, function(i){
    dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
    dt[,RunNumber := i]
    return(dt)
  }))
  
  
  dt = caseData[,list(DailyInfected = .N), by = list(Day = floor(InfectedHour/24),RunNumber ) ]
  dt = data.table(dt %>% complete(nesting(RunNumber), Day = seq(0, endDay, 1), fill = list(DailyInfected = 0)))
  
  print(mean(caseData[DetectedHour < 1000*24,DetectedHour - InfectedHour]))
  #ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
  ggplot(dt) + geom_line(aes(x = Day, y = DailyInfected, group = as.factor(RunNumber))) #  + geom_smooth(aes(x = Day, y = DailyInfected))
  
}

simulateOutbreaks(params = params,numOutbreaks = 30)

