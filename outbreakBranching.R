
require(data.table)
require(ggplot2)
require(Rcpp)

Rcpp::sourceCpp("~/MassTesting/viralLoad.cpp")



params = c(
           ProbDetectSymptoms = 0.95,
           ProbTracedGivenInfectorDetected = 0.95,
           RelativeTransmissionRisk_Detected = 0.1,
           ContactTracingDelay = 1,
           TestDelay = 1,
           
           
           # viral load params
           contactsPerHour = 13/24, 
           testDelay = 12, 
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


caseData = data.table(branchingModel(endDay = 60, popSize = 1e4, params)) 

ggplot(caseData[DetectedDay < 1000], aes(x = DetectedDay - InfectedDay)) + geom_histogram()
ggplot(caseData, aes(x = InfectedDay)) + geom_histogram()