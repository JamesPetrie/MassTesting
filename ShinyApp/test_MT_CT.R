#test file for MT & CT
source("/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/ViralLoad.cpp")

params = c(
  normalTestPeriod = input$normalTestPeriod*24  ,
  outbreakTestPeriod = input$outbreakTestPeriod*24 ,
  ContactTracingDelay = input$tracingDelay ,
  ProbTracedGivenInfectorDetected = input$fractionTraced,
  ProbDetectSymptoms = input$probDetectSymptoms,
  timeFromPeakToSymptoms = input$timeFromPeaktoSymptoms,
  timeToPeak = input$daysToPeak*24,
  maxTimeAfterPeak= 24*30, 
  
  
  timeFromPeakTo0 = 24*5,
  contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, 
  precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
  relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
  maskEffect = input$maskEffect, minLogPCRViralLoad = input$minLogPCRViralLoad, initialLogLoad = input$initialLogLoad)



#output final dt 
params["logPeakLoad"] = computePeakViralLoad(params["timeToPeak"], targetR0 =R0, params)
caseData = rbindlist(llply(1:numOutbreaks, function(i){
  dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
  dt[,RunNumber := i]
  return(dt)
}))

dt

