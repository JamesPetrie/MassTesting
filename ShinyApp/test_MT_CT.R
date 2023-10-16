library(boot)
library(data.table)
library(plyr)

Rcpp::sourceCpp("~/MassTesting/ShinyApp/ViralLoad.cpp")

numOutbreaks = 100
endDay = 300
maxSize = 1000

params = c(ProbTrace = 1,ProbIso = 1, R0 = 2.8, IsoAndTraceSameTime = TRUE)
  
    
caseData = rbindlist(llply(1:numOutbreaks, function(i){
  dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
  dt[,RunNumber := i]
  return(dt)
}))
 
caseData = caseData[InfectorId != -1 & (InfectedHour + 5*24*8 < endDay*24)]
Re <- mean(caseData$NumInfected)

bo <- boot(data.frame(caseData[,NumInfected]), statistic=meanfun, R=1000)
bootResult <- boot.ci(bo, conf=0.95, type = c("basic"))
lowerBound <- bootResult$basic[1,4]
upperBound <- bootResult$basic[1,5]


result = data.table(Re = Re,  ReLowerBound = lowerBound, ReUpperBound = upperBound,R0 = params["R0"],ProbTrace = params["ProbTrace"], ProbIso = params["ProbIso"], IsoAndTraceSameTime = 1 == params["IsoAndTraceSameTime"])

print(result)



