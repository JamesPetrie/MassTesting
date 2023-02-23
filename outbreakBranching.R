





params = c(
           ProbDetectSymptoms = 0.5,
           ProbTracedGivenInfectorDetected = 0.5,
           RelativeTransmissionRisk_Detected = 0.05,
           ContactTracingDelay = 12 ,
           
           
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
           timeToPeak = 24*5,
           timeFromPeakToSymptoms = 0,
           timeFromPeakTo0 = 24*5,
           initialLogLoad = -2, 
           minLogPCRViralLoad = 3,
           precision = 0.15)

 

simulateOutbreaks(params = params,numOutbreaks = 600)

