library(boot)
library(data.table)
library(plyr)

Rcpp::sourceCpp("~/MassTesting/ShinyApp/ViralLoad.cpp")

numOutbreaks = 5
endDay = 150
maxSize = 500

params = c(ProbTrace = 0,ProbIso = 0, R0 = 3, MeanDaysToSymptoms = 8,MeanDaysToTransmission = 7,IsoAndTraceSameTime = TRUE)
  
    
caseData = rbindlist(llply(1:numOutbreaks, function(i){
  print(i)
  dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
  dt[,RunNumber := i]
  return(dt)
}))
 
caseData = caseData[InfectorId != -1 & (InfectedHour + 5*24*8 < endDay*24)]
Re <- mean(caseData$NumInfected)

meanfun <- function(data, i){
  d <- data[i, ]
  return(mean(d))   
}

bo <- boot(data.frame(caseData[,NumInfected]), statistic=meanfun, R=1000)
bootResult <- boot.ci(bo, conf=0.95, type = c("basic"))
lowerBound <- bootResult$basic[1,4]
upperBound <- bootResult$basic[1,5]


result = data.table(Re = Re,  ReLowerBound = lowerBound, ReUpperBound = upperBound,R0 = params["R0"],ProbTrace = params["ProbTrace"], ProbIso = params["ProbIso"], IsoAndTraceSameTime = 1 == params["IsoAndTraceSameTime"])

print(result)



