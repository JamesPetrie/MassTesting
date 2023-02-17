
require(data.table)
require(ggplot2)
require(Rcpp)

Rcpp::sourceCpp("~/MassTesting/viralLoad.cpp")



params = c(
           ProbDetectSymptoms = 0.0,
           ProbTracedGivenInfectorDetected = 0.2,
           RelativeTransmissionRisk_Detected = 0.1,
           ContactTracingDelay = 24,
           
           
           # viral load params
           testPeriod = 20,
           contactsPerHour = 13/24, 
           testDelay = 24, 
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


caseData = data.table(branchingModel(endDay = 60, maxSize = 500, params)) 

ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
ggplot(caseData, aes(x = InfectedHour)) + geom_histogram()