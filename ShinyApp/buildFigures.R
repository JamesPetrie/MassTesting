require(ggforce)
require(ggplot2)
require(cowplot)
require(scales)

require(data.table)
require(Rcpp)
require(plyr)
require(tidyr)

folder =  "~/MassTesting/ShinyApp/" # "/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/"

source(paste0(folder, "ViralLoad.R"))
#source("~/MassTesting/outbreakBranching.R")
Rcpp::sourceCpp(paste0(folder, "ViralLoad.cpp"))

theme_set(theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) + theme(legend.background = element_rect(fill = "white")) + theme_half_open() + background_grid()  + 
            theme(text = element_text(size=20), axis.text = element_text(size=20)))


generateControllabilityFigure = function(testPeriods, inputParams){
  dt = rbindlist(llply(testPeriods, function(testPeriod){
    params = copy(inputParams)
    params["testPeriod"] = testPeriod
    evaluateStrategy(params)
  }))
  if(nrow(dt)>0){
    
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24), PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    freqNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    
    dt[,DaysToPeak := TimeToPeak/24]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wuhan)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
      DaysToPeak = c(8, 5.5, 3.5, 10, 10), R0 = c(2.5, 8, 2.5, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak + 0.5*c(1,0,-1,0), R0 = R0 + R0*0.2*c(0, 1,0,-1)), by = Pathogen ]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = PeriodLabel), linewidth = 1.4) + scale_y_log10(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(0.98,20)) + scale_x_continuous(breaks = 0:12) + 
      guides(colour=guide_legend(title="Test period [Days]")) + xlab("Time to Peak Viral Load [Days]") + ylab("R0")+
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, fill = Pathogen, label = Pathogen),size = 0.01 ,label.fontsize = 14, show.legend = F)+  theme(legend.position="bottom") +
      labs(title = "Effect of Mass Testing",  subtitle = "Maximum controllable R0 for different testing strategies")
      #geom_ellipse(data = data.table(), aes(x0 = 5, y0 = 3, a = 1, b = 1, angle = 0), fill = "orange", alpha = 0.4)
    
    #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
  }
  return(p)
}

replaceParams = function(params, timeToPeak,  logPeakLoad){
  newParams = copy(params)
  newParams["logPeakLoad"] = logPeakLoad
  newParams["timeToPeak"] = timeToPeak
  newParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
  return(newParams)
}

computePeakViralLoad = function(timeToPeak, targetR0,params){
  testParams = copy(params)
  testParams["timeToPeak"] = timeToPeak
  if(! "timeFromPeakTo0" %in% names(testParams)){
    testParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
    
  }
  
  cutoffLoad = bisect(a = 0, b = 20, maxiter = 15, fun = function(logPeakLoad){
    testParams["logPeakLoad"] = logPeakLoad
    
    
    R0 = sumTransmissions(0,timeToPeak + testParams["timeFromPeakTo0"], testParams)
    
    return(targetR0 - R0)
  })$root
  
  return(cutoffLoad)
  
}

plotTrajectories = function(params){
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(3,7,10), Time = 24*seq(0,16, length.out = 200))) 
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =3, params), by = TimeToPeak]
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*probTransmit(ViralLoad,params), by = list( ViralLoad, TimeToPeak) ]
  dt[ ,TestSensitivity :=probPositive(ViralLoad, params), by = list( ViralLoad, TimeToPeak) ]
  
  
  peakNames = dt[, list(PeakLabel = paste("Days to Peak: ", TimeToPeak/24)), by = TimeToPeak]
  setkey(peakNames, by = "TimeToPeak")
  peakNames[, PeakLabel := factor(PeakLabel, levels = PeakLabel)]
  dt = merge(dt, peakNames, by = "TimeToPeak")
  
  
  dt[, LogViralLoad := log10(ViralLoad)]
  dtLong = melt(dt, id.vars = c("Time", "PeakLabel"), measure.vars = c("LogViralLoad", "TestSensitivity","DailyTransmissions"))
  #dtLong[, value:=value/max(value), by = variable]
  dtLong[variable == "LogViralLoad", variable := "Log Viral Load"]
  dtLong[variable == "TestSensitivity", variable := "Test Sensitivity"]
  dtLong[variable == "DailyTransmissions", variable := "Daily Transmissions"]
  
  #ggplot(dtLong, aes(x = Time/24, y = value)) + facet_wrap( ~   PeakLabel + variable , nrow = 3, dir = "v") + geom_line() + 
  #theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.y = element_blank()) + xlab("Day Since Infection" ) +
  # scale_x_continuous(breaks = seq(0,16, by = 2))+
  # theme(strip.background = element_blank()) + ggtitle("Example viral load trajectory, test sensitivity, expected transmissions")
  
  p1 = ggplot(dt, aes(x = Time/24, y = LogViralLoad)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x = element_blank()) +# xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
     theme(strip.background = element_blank())  + ylab("Viral Load\n(log10 copies / ml)")
  
  p2 = ggplot(dt, aes(x = Time/24, y = TestSensitivity)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x = element_blank()) +# xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank())  + ylab("Test\nSensitivity")
  
  p3 = ggplot(dt, aes(x = Time/24, y = DailyTransmissions)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab("Expected\nTransmissions")

  plot_grid(p1,p2,p3, ncol = 1,align = "v",rel_heights = c( 1,1,1.1)) 
}

plotInfectiousness = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 15, length.out = 100))
  dt[ ,ExpectedTransmissions :=params["contactsPerHour"]*24*probTransmit(ViralLoad, params), by = ViralLoad ]
  
  ggplot(dt , aes(x = ViralLoad, y = ExpectedTransmissions)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + xlab("Viral Load [RNA copies/ml]") + ylab("Expected\nTransmissions") +
    theme(text = element_text(size=12), axis.text = element_text(size=12)) + labs(title="Expected Transmissions per Day")
  
}

plotTestSensitivity = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 15, length.out = 100))
  dt[ ,ProbPositive :=probPositive(ViralLoad, params), by = ViralLoad ]
  
  ggplot(dt , aes(x = ViralLoad, y = ProbPositive)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + xlab("Viral Load [RNA copies/ml]") + ylab("Test\n Sensitivity") + 
    theme(text = element_text(size=12), axis.text = element_text(size=12)) + labs(title="PCR Test Sensitivity")
  
}


# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotFracTransmissionsAfterPositive = function(testPeriods,params){
  dt = data.table(expand.grid(TestDelay = seq(0, 72, length.out = 12), TestPeriod  = testPeriods, TimeToPeak = 24*c(3,7,10)))
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =3, params), by = TimeToPeak]
  
  dt = rbindlist(llply(1:nrow(dt), function(i){

    newParams = copy(params)
    newParams["logPeakLoad"] = dt[i, LogPeakLoad]
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  peakNames = dt[, list(PeakLabel = paste("Days to Peak: ", TimeToPeak/24)), by = TimeToPeak]
  setkey(peakNames, by = "TimeToPeak")
  peakNames[, PeakLabel := factor(PeakLabel, levels = PeakLabel)]
  dt = merge(dt, peakNames, by = "TimeToPeak")
  p = ggplot(dt, aes(x = TestDelay, y = FracAfterPositive, colour =  as.factor(TestPeriod/24))) + geom_line(linewidth = 1.4)  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.01)) + guides(colour=guide_legend(title="Test Period [Days]")) +
    xlab("Test Delay [Hours]") + ylab("Fraction Transmissions \n After Positive Test")+ 
    theme(legend.position = c(0.7, 0.2)) + theme(text = element_text(size=16)) + facet_wrap(~PeakLabel, nrow = 1)+
    theme(strip.background = element_blank())+  theme(legend.position="bottom") 
  return(p)
}

plotPrevalenceCost = function(testPeriods, params){
  # need cost per test
  # need cost per isolation
  # need test frequencies
  
  variableTestCost = params["variableTestCost"]
  isoCost = params["isolationCost"]
  fixedAnnualizedDailyTestCost = params["fixedAnnualizedDailyTestCost"]
  
  dt = data.table(expand.grid(DailyFracInfected = 10^seq(-7, -1, length.out = 200), TestPeriod  = testPeriods))
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    dailyInfections = dt[i, DailyFracInfected]
    testPeriod = dt[i,TestPeriod]
   
    
    
    testFreq = 1/testPeriod
    

    testCost = testFreq*variableTestCost

    
     
    dailyTestCostPerPerson = testCost
    dailyIsoCostPerPerson = dailyInfections*isoCost
    
    dailyTotalCost = dailyTestCostPerPerson + dailyIsoCostPerPerson
    
    
    
    return(data.table(TotalCost = dailyTotalCost, TestCost = dailyTestCostPerPerson, IsoCost = dailyIsoCostPerPerson,  TestPeriod = testPeriod, FractionInfectedDaily = dailyInfections))
  } ))
  
  dt[,AnnualizedFixedCostPerPerson := fixedAnnualizedDailyTestCost/TestPeriod]

  dt[,PeriodDescription := paste0(TestPeriod, " ($",round(AnnualizedFixedCostPerPerson, digits = 3), " per person per year fixed cost)" )]
  
  periodNames = dt[, list(PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
  periodNames[, TestPeriod := as.numeric(TestPeriod)]
  setkey(periodNames, by = "TestPeriod")
  periodNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
  dt = merge(dt, periodNames, by = "TestPeriod")
  
  
  dtLong = melt(dt, id.vars = c(
    "FractionInfectedDaily", "PeriodDescription"), measure.vars = c("TotalCost",  "IsoCost") )
  
  dtLong[ variable == "TotalCost", variable := "Total Cost"]
  dtLong[ variable == "IsoCost", variable := "Isolation Cost"]
  
 
  
  
  ggplot() +geom_line(data = dtLong, aes(x= FractionInfectedDaily, y = value, linetype = variable, group = paste(variable,PeriodDescription)),  linewidth = 2) + geom_line(data = dtLong[variable == "Total Cost"], aes(x= FractionInfectedDaily, y = value, colour = PeriodDescription), linewidth = 2) +
    scale_x_log10() +scale_y_log10(limits = c(0.05, 64), n.breaks = 8) + labs(x = "Daily Fraction of Population Infected", y = "Daily Cost per Person (USD)") + 
    guides(colour = guide_legend(nrow = 2)) + theme(legend.position = c(0.1, 0.8)) +guides(colour=guide_legend(title="Test period [Days]"),linetype=guide_legend(title="Cost Type")) +
    ggtitle("Daily Testing and Isolation Cost")
  
}

# plotImportCost = function(){
#   
#   sampleNumGenerations = function(Re, overdispersion=0.1, initialSize = 1){
#     if(Re >= 1.0) return(NA)
#     numInfected = initialSize
#     
#     count = 0
#     totalInfected = 0
#     while(count < 1e4 && numInfected > 0){
#       totalInfected = totalInfected + numInfected
#       count = count + 1
#       numInfected = sum(rnbinom(numInfected, size = overdispersion, mu = Re ))
#     }
#     return(data.table(NumGenerations = count, TotalInfected = totalInfected ))
#   }
#   
#   
#   estimateOutbreakDuration = function(Re, initialSize, tau = 5, safetyMargin = 10){
#     numGen = mean(sapply(1:10000, function(x) {dt = sampleNumGenerations(Re = Re, initialSize = initialSize); return(dt$NumGenerations)}))
#     return(tau*numGen + safetyMargin)
#   }
#   #mean(sapply(1:10000, function(x) {dt = sampleNumGenerations(Re = 0.9, initialSize = 10); return(dt$TotalInfected)}))
#   
#   
#   steadyImportCost = function(importRate, popSize, Re, controlTime, dailyCost, infectionCost, detectionSize,  fracPopTarget, probExportPerInfected){
#     
#     numDailyImports = importRate*popSize
#     
#     if(is.na(controlTime)){
#       return(dailyCost + infectionCost*importRate*1/(1-Re))
#     }else{
#       outbreakSize = detectionSize/(1-Re)
#       outbreakCost = infectionCost*outbreakSize/popSize + fracPopTarget*controlTime*dailyCost
#       
#       if(probExportPerInfected*outbreakSize >= 1){
#         
#         return(dailyCost + infectionCost*importRate*1/(1-Re))
#       }
#       
#       
#       expectedNewOutbreaks = 1/(1 - probExportPerInfected*outbreakSize)
#       
#       
#       
#       
#       return(numDailyImports*outbreakCost*expectedNewOutbreaks) # todo: consider recursive outbreak seeding
#     }
#   }
#   
#   popSize = 60e6
#   fracDt = fread("~/MassTesting/FracTransmissions.csv")
#   R0 = 2
#   incentiveCost = 5
#   logisticCost = 8
#   poolSize = 24
#   pcrCost = 48
#   isoCost = 5e3
#   testFraction = 0.95
#   isoFraction = 0.9
#   testPeriod = 2
#   testDelay = 1
#   
#   dt = rbindlist(llply(10^seq(-11, -2, length.out = 200), function(importRate){
#     
#     transFraction = fracDt[TestPeriod == testPeriod & TestDelay == testDelay, FracAfterPositive]
#     
#     testFreq = 1/testPeriod
#     
#     Re = R0*(1 - testFraction*transFraction*isoFraction)
#     
#     if(Re > 1){
#       cost = NA
#       dailyTestCostPerPerson = NA
#       dailyIsoCostPerPerson = NA
#     }
#     else{
#       newInfectFrac = importRate*1/(1-Re) # fraction new infections (import and transmit per day) relative to total pop
#       
#       if(poolSize == 1){
#         testCost = testFreq*testFraction*(incentiveCost + logisticCost + pcrCost)
#       }else{
#         testCost = testFreq*testFraction*(incentiveCost + logisticCost + 1/poolSize*pcrCost + poolSize*newInfectFrac*testPeriod*pcrCost)
#       }
#       
#       gdpPerCapita = 70e3 # gdp per person in USA
#       
#       dailyTestCostPerPerson = testCost*365/gdpPerCapita
#       dailyIsoCostPerPerson = newInfectFrac*testFraction*isoCost*365/gdpPerCapita
#       
#       cost = dailyTestCostPerPerson + dailyIsoCostPerPerson
#     }
#     
#     
#     return(data.table(Cost = cost, TestCost = dailyTestCostPerPerson, IsoCost = dailyIsoCostPerPerson, Re = Re, TestPeriod = testPeriod, ImportRate = importRate, TestDelay = testDelay, FractionInfectedDaily = importRate/(1-Re)))
#     
#   } ))
#   
#   dt[, Strategy := "Continuous Testing"]
#   
#   
#   
#   
#   dtTimeTarget = copy(dt)
#   dtTimeTarget[, Strategy := "Temporal Targeting"]
#   dtTimeTarget[, FracPopTarget := 1]
#   dtTimeTarget[, ProbExportPerInfected := 0]
#   
#   dtSpatioTempTarget = copy(dt)
#   dtSpatioTempTarget[, Strategy := "Spatio-temporal Targeting (10k person divisions)" ]
#   dtSpatioTempTarget[, FracPopTarget := 10e3/60e6]
#   dtSpatioTempTarget[, ProbExportPerInfected := 0.02]
#   
#   
#   dtTarget = rbind(dtSpatioTempTarget, dtTimeTarget)
#   
#   
#   detectionSize = 10
#   dtTarget[ , LocalOutbreakDuration := estimateOutbreakDuration(Re, detectionSize, 5, 10), by = Re]
#   
#   dtTarget[, Cost := steadyImportCost(ImportRate, popSize, Re, LocalOutbreakDuration, TestCost, 5e3/70e3, detectionSize, FracPopTarget, ProbExportPerInfected), by = list(Re, LocalOutbreakDuration, TestCost, ImportRate, TestPeriod, TestDelay,FracPopTarget, ProbExportPerInfected )]
#   
#   dt = rbind(dt, dtTarget, fill = TRUE)
#   
#   dt[,TestPeriod := as.factor(TestPeriod)]
#   ggplot(dt[ TestPeriod ==2 & TestDelay == 1],aes(x= ImportRate, y = Cost, colour =Strategy)) + geom_line(linewidth = 2) + 
#     geom_vline(xintercept = 1e-5, linetype = "dashed", linewidth = 2) + annotate("text", x=6.5*1e-5, y=0.01, label="1000 cases \nper day (UK)", size = 6 ) +
#     geom_vline(xintercept = 1e-9, linetype= "dashed", linewidth = 2) +annotate("text", x=6.5*1e-9, y=0.01, label="0.1 cases \nper day (UK)", size = 6 )+
#     theme(legend.position = "bottom") +  guides(colour = guide_legend(nrow = 3)) +
#     scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + scale_y_log10(limits = c(0.001, 0.5), labels = trans_format("log10", math_format(10^.x)))+
#     labs(x = "Daily Imported Cases Relative to Population Size", y = "Daily Fraction of GDP") + 
#     theme(text = element_text(size=24), axis.text = element_text(size=24))
#     #+ylim(0,0.3)
#   
#   
#   
# }


plotOutbreaks = function(numOutbreaks = 10, endDay = 90, maxSize = 1000, R0, params){
  #output final dt 
    params["logPeakLoad"] = computePeakViralLoad(params["timeToPeak"], targetR0 =R0, params)
    caseData = rbindlist(llply(1:numOutbreaks, function(i){
      dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
      dt[,RunNumber := i]
      return(dt)
    }))
    
    caseData[, OutbreakDetectedHour := min(DetectedHour), by = RunNumber]
    caseData[, OutbreakSize := max(Id) +1 , by = RunNumber]
    
    outbreakOverDt = caseData[OutbreakSize < maxSize, list(OutbreakOverHour = max(InfectedHour) + 14*24), by = RunNumber]
    days = seq(0, endDay, 1)
    fractionDone = sapply(days, function(day){
      return(sum(outbreakOverDt$OutbreakOverHour <= day*24)/numOutbreaks)
    })
    p1 = ggplot(data.table(FractionDone = fractionDone, Day =days)) + geom_line(aes(x = Day, y = FractionDone)) + ylim(0,1)
    
    
    

    dt = caseData[,list(DailyInfected = .N), by = list(Day = floor(InfectedHour/24),RunNumber ) ]
    dt = data.table(dt %>% complete(nesting(RunNumber), Day = seq(0, endDay, 1), fill = list(DailyInfected = 0)))
    setkeyv(dt, c("Day", "RunNumber"))
    dt[ , CumulativeInfected := cumsum(DailyInfected), by = RunNumber]
    
    
    outbreakDetectedDt = caseData[, list(OutbreakDetectedHour = min(OutbreakDetectedHour)), by = RunNumber]
    
    #print(mean(caseData[DetectedHour < 1000*24,DetectedHour - InfectedHour]))
    #ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
    
    #ggplot(dt) + geom_line(aes(x = Day, y = CumulativeInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = CumulativeInfected)) #+ scale_y_log10()
    #p1 = ggplot(dt) + geom_line(aes(x = Day, y = DailyInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = DailyInfected))# + geom_vline(data = outbreakDetectedDt[OutbreakDetectedHour < 24*1000], aes(xintercept = OutbreakDetectedHour/24), alpha = 0.2, colour = 'green')
    
    # days of undetected infection per outbreak
    sumData = caseData[, list(DaysUndetected = sum(pmin(15*24, DetectedHour - InfectedHour, TracedHour - InfectedHour ))/24), by = RunNumber]; 
    meanUndetected = mean(sumData$DaysUndetected)
    p2 = ggplot(sumData, aes(x= DaysUndetected)) + geom_histogram() + geom_vline(xintercept = meanUndetected, colour = "purple", linewidth =2) + xlab("Number of Person-Days \n Infected but Undetected") +
      scale_x_log10( limits = c(1,1000))
    
    preDetectDt = caseData[, list(NumInfectedBeforeOutbreakDetected = sum(InfectedHour <= OutbreakDetectedHour)-1), by = RunNumber]
    p3 = ggplot(preDetectDt, aes(x= NumInfectedBeforeOutbreakDetected)) + geom_histogram() + geom_vline(xintercept = mean(preDetectDt$NumInfectedBeforeOutbreakDetected), colour = "purple", linewidth =2) +
      xlab("Number Infected \n Before Outbreak Detected") +scale_x_continuous(breaks = seq(0,16, by = 1), limits = c(0,16))
    
    casesAfterOutreakDetected = caseData[InfectedHour > OutbreakDetectedHour & (InfectedHour + params["timeToPeak"] + params["timeFromPeakTo0"] <= SimulationEndHour), list(NumInfected)]
    #re redo using dt ( re after the outbreak is detected)
    bootStrapRe = sapply(1:50, function(i){mean(sample(casesAfterOutreakDetected$NumInfected, size = nrow(casesAfterOutreakDetected), replace = TRUE))})
    p4 = ggplot(data = data.table(ReEstimate = bootStrapRe), aes(x= ReEstimate)) + geom_histogram() + scale_x_continuous(limits = c(0,2)) + xlab("Bootstrapped Re Estimate \n (After Outbreak Detected)") + 
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) + ylab("Fraction of Samples")
    
    
    # todo: days of elevated testing per outbreak
    
    # todo: cost of isolation and quarantine per outbreak
    
    # todo: number of infections per outbreak
    
    # todo: number of infections before outbreak detected
    # todo: remove people from Re estimate when still infected when simulation finished. Also add point
    
    
    plot_grid(p1,p2,p3,p4, ncol = 2) 
  }





inputParams = c(contactsPerHour = 13/24, testDelay = 12, fracIso = 0.9, fracTest = 0.9, precision = 0.15, maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30, logPeakLoad = 10, initialLogLoad = -2, minLogPCRViralLoad = 3)

#generateControllabilityFigure(c(24, 48), inputParams)
#plotTrajectories(inputParams)
#plotImportCost()

#plotFracTransmissionsAfterPositive(c(24, 48), inputParams)



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

#plotOutbreaks(params=params, numOutbreaks = 300)
