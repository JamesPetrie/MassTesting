
require(data.table)
require(ggplot2)
require(Rcpp)


Rcpp::sourceCpp("~/MassTesting/outbreakBranching.cpp")


params = c(
           ProbDetectSymptoms = 0.15,
           ProbTracedGivenInfectorDetected = 0.5,
           
           DailyTransmissionRate = 0.5,
           TimeToInfectiousRate = 0.4,
           RecoveryRate = 0.2,
           
           RelativeTransmissionRisk_Detected = 0.2,
           
           InfectionDetectionDelay = 2,
           ContactTracingDelay = 2,
)

caseData = data.table(branchingModel(endDay = 30, popSize = 1e5, params)) 