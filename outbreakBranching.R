
require(data.table)
require(ggplot2)
require(Rcpp)
require(plyr)
require(tidyr)

Rcpp::sourceCpp("~/MassTesting/viralLoad.cpp")



params = c(
           ProbDetectSymptoms = 0.0,
           ProbTracedGivenInfectorDetected = 0.3,
           RelativeTransmissionRisk_Detected = 0.05,
           ContactTracingDelay = 0,
           
           
           # viral load params
           normalTestPeriod = 96 ,
           outbreakTestPeriod = 24 ,
           contactsPerHour = 13/24, 
           testDelay = 8, 
           fracIso = 0.9,
           fracQuar = 0.9, 
           fracTest = 0.9, 
           maxProbTransmitPerExposure = 0.3, 
           relativeDeclineSlope = 1.0, 
           maxTimeAfterPeak= 24*30, 
           logPeakLoad = 7.0, 
           timeToPeak = 24*4,
           timeFromPeakToSymptoms = 0,
           timeFromPeakTo0 = 24*4,
           initialLogLoad = -2, 
           minLogPCRViralLoad = 3,
           precision = 0.15)

simulateOutbreaks = function(numOutbreaks = 10, endDay = 90, maxSize = 1000, params){
  caseData = rbindlist(llply(1:numOutbreaks, function(i){
    dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
    dt[,RunNumber := i]
    return(dt)
  }))
  
  
  dt = caseData[,list(DailyInfected = .N), by = list(Day = floor(InfectedHour/24),RunNumber ) ]
  dt = data.table(dt %>% complete(nesting(RunNumber), Day = seq(0, endDay, 1), fill = list(DailyInfected = 0)))
  
  #print(mean(caseData[DetectedHour < 1000*24,DetectedHour - InfectedHour]))
  #ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
  sumData = caseData[, sum(pmin(15*24, DetectedHour - InfectedHour, TracedHour - InfectedHour ))/24, by = RunNumber]; 
  print(mean(sumData$V1))
  ggplot(dt) + geom_line(aes(x = Day, y = DailyInfected, group = as.factor(RunNumber)))  + geom_smooth(aes(x = Day, y = DailyInfected))
  
}

simulateOutbreaks(params = params,numOutbreaks = 100)

