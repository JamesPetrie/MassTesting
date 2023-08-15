require(ggforce)
require(ggplot2)
require(cowplot)
require(scales)

require(data.table)
require(Rcpp)
require(plyr)
require(tidyr)
require(metR)


source("ViralLoad.R")
#source("~/MassTesting/outbreakBranching.R")
#Rcpp::sourceCpp("ViralLoad.cpp")

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
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, group = Pathogen, label = Pathogen), fill = "plum3",size = 0.00 ,label.fontsize = 14, show.legend = F,  lty = "blank")+  theme(legend.position="bottom") +
      labs(title = "Effect of Mass Testing",  subtitle = "Maximum controllable R0 for different testing strategies")
      #geom_ellipse(data = data.table(), aes(x0 = 5, y0 = 3, a = 1, b = 1, angle = 0), fill = "orange", alpha = 0.4)
    
    #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
  }
  return(p)
}


generate2TestControllabilityFigure = function(testPeriods, params){
  dt = data.table(expand.grid(TestPeriod = testPeriods, TestType = c("Antigen", "PCR")))
  
  dt[TestType == "Antigen", LOD := 5]
  dt[TestType == "Antigen", TestDelay := 0]
  dt[TestType == "PCR", LOD := 3]
  dt[TestType == "PCR", TestDelay :=8]
  
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["logLimitOfDetection"] = dt[i,LOD]
    
    
    
    
    return(merge(dt[i,], evaluateStrategy(newParams)))
    
  }))
  

  if(nrow(dt)>0){
    
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24), PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    freqNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    
    dt[,DaysToPeak := TimeToPeak/24]
    dt[, TestType := factor(TestType, levels = c("PCR", "Antigen"))]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wuhan)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
                            DaysToPeak = c(8, 5.5, 3.5, 10, 10), R0 = c(2.5, 8, 2.5, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak + 0.5*c(1,0,-1,0), R0 = R0 + R0*0.2*c(0, 1,0,-1)), by = Pathogen ]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = as.factor(TestPeriod/24), linetype = TestType), linewidth = 1.4) +
      scale_y_log10(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(0.98,20)) + scale_x_continuous(breaks = 0:12) + 
      guides(colour=guide_legend(title="Test period [Days]")) + guides(linetype=guide_legend(title="Test Type"))  + xlab("Time to Peak Viral Load [Days]") + ylab("R0")+
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, group = Pathogen, label = Pathogen),fill = "plum3",size = 0.0 ,label.fontsize = 14, show.legend = F, lty = "blank")+  
      theme(legend.position="bottom") + 
      labs( subtitle = "Maximum controllable R0 for different testing strategies")
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
  testParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
  
  cutoffLoad = bisect(a = 0, b = 20, maxiter = 15, fun = function(logPeakLoad){
    testParams["logPeakLoad"] = logPeakLoad
    
    
    R0 = sumTransmissions(0,timeToPeak + testParams["timeFromPeakTo0"], testParams)
    
    return(targetR0 - R0)
  })$root
  
  return(cutoffLoad)
  
}

wrangleFracAfterPositive = function(timeToPeak, testPeriod, lod, testDelay, params){
  newParams = copy(params)
  newParams["timeToPeak"] = timeToPeak
  newParams["timeFromPeakTo0"] = timeToPeak
  newParams["testPeriod"] = testPeriod
  newParams["logLimitOfDetection"] = lod
  newParams["testDelay"] = testDelay
  
  return(fracAfterPositive(newParams))
}

plotFracReduction = function(params, testPeriods = c(24, 72), timesToPeak = 24*c(3,6,9), n = 30, showPooled = FALSE){

  dt = data.table(expand.grid(TimeToPeak = timesToPeak, TestPeriod = testPeriods, LOD = seq(2,6, length.out = n), TestDelay = seq(0,48, length.out = n))) 
  
  dt[, FracAfterPositive :=  wrangleFracAfterPositive(TimeToPeak, TestPeriod, LOD, TestDelay, params), by = list(TimeToPeak, TestPeriod, LOD, TestDelay) ]
  # for each of 2, 6, 10 peak time, peak VL of 8
  # for each of 1, 2, 3 day test period
  # for range of test delays
  # for range of LODs
  # for adherence = 1
  # expand.grid
  # compute fraction reduction for each row
  
 

  breaks = c(1,  0.99, 0.98, 0.96, 0.92, 0.84, 0.68, 0.36, 0)
  p = ggplot(dt, aes(TestDelay, LOD, z=  FracAfterPositive)) +scale_y_continuous(breaks = 1:6, labels = label_math(expr=10^.x) ) +
    geom_contour_filled(breaks = breaks, alpha = 0.6 ) + 
    
    #geom_point(aes(x = 0, y = 5), size = 3.5, colour = "blue")+
      
    #geom_hline(yintercept = 3, linetype = "dashed", colour = "blue", size = 1.0) + 
    #geom_textcontour(breaks = breaks, straight = T , position = "jitter"  ) + 
    geom_contour(aes(z = FracAfterPositive),  breaks = breaks, colour = "grey15", linetype = "dashed") + 
    geom_text_contour(aes(z = FracAfterPositive),  breaks = breaks,rotate = TRUE, nudge_x = 1.5, nudge_y = 0.1, skip = 0, colour = "black", label.placer = label_placer_fraction(frac = 0.99)) + 
    scale_fill_manual(values = terrain.colors(11)) + 
    #theme(legend.position = "none")  +
    guides(fill=guide_legend(title="Fraction Reduction\nin Transmissions"))+
    theme(legend.position = "bottom")+
    xlab("Test Delay (hours)") + ylab("Limit of Detection (copies/ml)")  
  
  if(length(testPeriods) == 1){
    p = p + facet_wrap( ~paste("Days to Peak Viral Load:", TimeToPeak/24 ))+ theme(strip.background = element_blank())
  }else{
    p = p + facet_grid( paste("Test Period (days):",  TestPeriod/24)  ~ paste("Days to Peak Viral Load:", TimeToPeak/24 ))+ theme(strip.background = element_blank())
  }
  

  
  
  p = p +   annotate(geom="rect",  xmin = 2, xmax= max(dt$TestDelay), ymin = 2.7, ymax=3.3, fill="blue", alpha=0.2) + 
    annotate("text", x = max(dt$TestDelay) - 5, y = 3.2, label = "PCR")
  
  if(showPooled){
    p = p +   annotate(geom="rect",  xmin = 4, xmax= max(dt$TestDelay), ymin = 3.7, ymax=4.3, fill="blue", alpha=0.2) + 
      annotate("text", x = max(dt$TestDelay) - 14, y = 4.2, label = "10x Pooled PCR")
  }
 
  p = p +   annotate(geom="rect",  xmin = 0, xmax= 2, ymin = 4.5, ymax=5.5, fill="blue", alpha=0.2) + 
    annotate("text", x = 5, y = 5.3, label = "Rapid\nAntigen")
  

  
  return(p)
}

plot3Trajectories = function(params){
  dt = data.table(expand.grid(TimeToPeak = 24*c(6), Time = 24*seq(0,16, length.out = 200))) 
  dt[, LogPeakLoad := 8] #computePeakViralLoad(TimeToPeak, targetR0 = 4.5, params), by = TimeToPeak]
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*probTransmit(ViralLoad,params), by = list( ViralLoad, TimeToPeak) ]
  dt[ ,TestSensitivity :=probPositive(ViralLoad, params), by = list( ViralLoad, TimeToPeak) ]
  
  
  p1 = ggplot(dt, aes(x = Time/24, y = ViralLoad)) + geom_line(linewidth = 1.4)   + xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+ theme(text = element_text(size=12), axis.text = element_text(size=12))+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    theme(strip.background = element_blank())  + ylab("Viral Load (copies / ml)")+ 
    labs(title="Viral Load Trajectory")
  
  
  p2 = ggplot(dt, aes(x = Time/24, y = TestSensitivity))  + geom_line(linewidth = 1.4) +  xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+theme(text = element_text(size=12), axis.text = element_text(size=12))+
    theme(strip.background = element_blank())  + ylab("Test Sensitivity")+ 
    labs(title="Test Sensitivity Trajectory")
  
  p3 = ggplot(dt, aes(x = Time/24, y = DailyTransmissions)) + geom_line(linewidth = 1.4) + 
    xlab("Day Since Infection" ) +theme(text = element_text(size=12), axis.text = element_text(size=12))+
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab("Expected Daily Transmissions")+
    labs(title="Infectiousness Trajectory")
  
  return(list(p1,p2,p3))
}

plotTrajectories = function(params){
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(3,7,10), Time = 24*seq(0,16, length.out = 200))) 
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = 4.5, params), by = TimeToPeak]
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
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab("Expected Daily\nTransmissions")

  plot_grid(p1,p2,p3, ncol = 1,align = "v",rel_heights = c( 1,1,1.1)) 
}

plotInfectiousness = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 12, length.out = 100))
  dt[ ,ExpectedTransmissions :=params["contactsPerHour"]*24*probTransmit(ViralLoad, params), by = ViralLoad ]
  
  p = ggplot(dt , aes(x = ViralLoad, y = ExpectedTransmissions)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
    xlab("Viral Load (copies/ml)") + ylab("Expected Daily\nTransmissions") +
    theme(text = element_text(size=12), axis.text = element_text(size=12)) 
  #+ labs(title="Infectiousness vs. Viral Load")
  return(p)
}

plotTestSensitivity = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 12, length.out = 100))
  dt[ ,ProbPositive :=probPositive(ViralLoad, params), by = ViralLoad ]
  
  p = ggplot(dt , aes(x = ViralLoad, y = ProbPositive)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
    xlab("Viral Load (copies/ml)") + ylab("Test Sensitivity") + 
    theme(text = element_text(size=12), axis.text = element_text(size=12)) 
  #+ labs(title="Test Sensitivity vs. Viral Load")
  return(p)
}


# todo: modify so that test frequency is varied, either PCR (LOD 3, Delay 12) or Antigen (LOD 5, delay 0) are used, and DaysToPeak in c(3,6,9)
# colour by days to peak, linetype by test type
# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotEffectTestFreq = function(params){
  dt = data.table(expand.grid(TestType = c("Antigen", "PCR"), TestPeriod  = floor(24/seq(0.1,2, length.out = 20)), TimeToPeak = 24*c(3,6,9)))
  dt[TestType == "Antigen", LOD := 5]
  dt[TestType == "Antigen", TestDelay := 0]
  dt[TestType == "PCR", LOD := 3]
  dt[TestType == "PCR", TestDelay :=12]
  
  #dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =4.5, params), by = TimeToPeak]
  dt[, LogPeakLoad := 8]
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["logPeakLoad"] = dt[i, LogPeakLoad]
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["logLimitOfDetection"] = dt[i,LOD]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  
  
  dt[TestType == "Antigen", TestDescription := paste0("10^5, 0")]
  
  dt[TestType == "PCR", TestDescription := paste0("10^3, 12")]
    
  
  
  
  p = ggplot(dt,aes(x = 24/TestPeriod, y = FracAfterPositive, colour = as.factor(TimeToPeak/24), linetype = TestDescription)) + geom_line(size = 0.9) +
    xlab("Tests Per Day") + ylab("Fraction Transmissions Prevented") + 
    guides(colour=guide_legend(title="Days to Peak\nViral Load"), linetype = guide_legend(title = "Limit of Detection,\nTest Delay (hours)")) + 
    scale_x_log10(breaks = c(1/8, 1/4, 0.5, 1,2), labels= c("1/8", "1/4", "1/2", "1", "2")) + theme(text = element_text(size=14), axis.text = element_text(size=14))
  
  
  #+scale_y_continuous(trans=logit_trans(), breaks = c(0.16, 0.5, 0.84, 0.92, 0.96, 0.98, 0.99, 0.995))
  
   return(p)
}


# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotEffectTestDelay = function(params){
  dt = data.table(expand.grid(TestDelay = seq(0, 72, length.out = 24), TestPeriod  = 24,LOD = c(3,5), TimeToPeak = 24*c(3,6,9)))
  
  #dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =4.5, params), by = TimeToPeak]
  dt[, LogPeakLoad := 8]
  dt = rbindlist(llply(1:nrow(dt), function(i){

    newParams = copy(params)
    newParams["logPeakLoad"] = dt[i, LogPeakLoad]
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["logLimitOfDetection"] = dt[i,LOD]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  
  p = ggplot(dt, aes(x = TestDelay, y = FracAfterPositive, colour =  as.factor(TimeToPeak/24), linetype = as.factor(paste0("10^", LOD)))) + geom_line(size = 0.9)+ 
    guides(colour=guide_legend(title="Days to Peak\nViral Load"), linetype=guide_legend(title="Limit of Detection\n(copies/ml)")) + ylab("Fraction Transmissions Prevented") +
    xlab("Test Delay (hours)") + theme(text = element_text(size=14), axis.text = element_text(size=14))
  
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

    
    gdpPerCapita = 70e3 # gdp per person in USA
    
    dailyTestCostPerPerson = testCost*365/gdpPerCapita
    dailyIsoCostPerPerson = dailyInfections*isoCost*365/gdpPerCapita
    
    dailyTotalCost = dailyTestCostPerPerson + dailyIsoCostPerPerson
    
    
    
    return(data.table(TotalCost = dailyTotalCost, TestCost = dailyTestCostPerPerson, IsoCost = dailyIsoCostPerPerson,  TestPeriod = testPeriod, FractionInfectedDaily = dailyInfections))
  } ))
  
  dt[,AnnualizedFixedCostPerPerson := fixedAnnualizedDailyTestCost/TestPeriod]

  dt[,PeriodDescription := paste0(TestPeriod, " ($",round(AnnualizedFixedCostPerPerson, digits = 3), " per person per year fixed cost)" )]
  
  dtLong = melt(dt, id.vars = c(
    "FractionInfectedDaily", "PeriodDescription"), measure.vars = c("TotalCost",  "IsoCost") )
  
  dtLong[ variable == "TotalCost", variable := "Total Cost"]
  dtLong[ variable == "IsoCost", variable := "Isolation Cost"]
  
 
  
  
  ggplot() +geom_line(data = dtLong, aes(x= FractionInfectedDaily, y = value, linetype = variable, group = paste(variable,PeriodDescription)),  linewidth = 2) + geom_line(data = dtLong[variable == "Total Cost"], aes(x= FractionInfectedDaily, y = value, colour = PeriodDescription), linewidth = 2) +
    scale_x_log10() +scale_y_log10(limits = c(0.002, 0.5), n.breaks = 8) + labs(x = "Daily Fraction of Population Infected", y = "Fraction of GDP") + 
    guides(colour = guide_legend(nrow = 2)) + theme(legend.position = c(0.1, 0.8)) +guides(colour=guide_legend(title="Test period [Days]"),linetype=guide_legend(title="Cost Type")) +
    ggtitle("Daily Testing and Isolation Cost During Outbreak")
  
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


plotOutbreaks = function(numOutbreaks = 10, endDay = 90, maxSize = 1000, params){
    caseData = rbindlist(llply(1:numOutbreaks, function(i){
      dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
      dt[,RunNumber := i]
      return(dt)
    }))
    
    
    dt = caseData[,list(DailyInfected = .N), by = list(Day = floor(InfectedHour/24),RunNumber ) ]
    dt = data.table(dt %>% complete(nesting(RunNumber), Day = seq(0, endDay, 1), fill = list(DailyInfected = 0)))
    setkeyv(dt, c("Day", "RunNumber"))
    dt[ , CumulativeInfected := cumsum(DailyInfected), by = RunNumber]
    
    
    
    
    #print(mean(caseData[DetectedHour < 1000*24,DetectedHour - InfectedHour]))
    #ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
    
    #ggplot(dt) + geom_line(aes(x = Day, y = CumulativeInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = CumulativeInfected)) #+ scale_y_log10()
    ggplot(dt) + geom_line(aes(x = Day, y = DailyInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = DailyInfected))
    
    # days of undetected infection per outbreak
    sumData = caseData[, sum(pmin(15*24, DetectedHour - InfectedHour, TracedHour - InfectedHour ))/24, by = RunNumber]; 
    print(mean(sumData$V1))
    ggplot(sumData, aes(x= V1)) + geom_histogram()
    
    # todo: days of elevated testing per outbreak
    
    # todo: cost of isolation and quarantine per outbreak
    
    # todo: number of infections per outbreak
}

plotMultiplePreventedTransmissions = function(params){
  newParams = copy(params)
  newParams["testDelay"] = 24
  newParams["testPeriod"] = 48
  p1 = plotPreventedTransmissions(newParams) + theme(text = element_text(size=14)) #+ ggtitle("Test Every 2 Days with 24 Hour Delay") 
  
  newParams["testDelay"] = 8
  newParams["testPeriod"] = 24
  p2 = plotPreventedTransmissions(newParams)  + theme(text = element_text(size=14)) # + ggtitle("Test Every Day with 8 Hour Delay")
  
  p = plot_grid(p1, p2, nrow= 1, labels = c("A", "B"))
  
}

# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotFracTransmissionsAfterPositive = function(testPeriods,params){
  dt = data.table(expand.grid(TestDelay = seq(0, 72, length.out = 12), TestPeriod  = testPeriods, TimeToPeak = 24*c(3,7,10)))
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =4.5, params), by = TimeToPeak]
  
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

plotPreventedTransmissions = function(params){

  dt = data.table(TimeToPeak = params["timeToPeak"])
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = 4.5, params), by = TimeToPeak]
  dt = rbindlist(llply(1:nrow(dt), function(i){
    fracDiscovered = fracDiscoveredByHour( replaceParams(params, dt[i, TimeToPeak], dt[i,LogPeakLoad])) 
    hours = 0:(length(fracDiscovered) -1)
    return(cbind(data.table(Time = hours, FracDiscovered = fracDiscovered), dt[i,]))
  }))
  
  
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions := 24*params["contactsPerHour"]*probTransmit(ViralLoad,params), by = list( ViralLoad, TimeToPeak) ]
  
  fracAdhere = params["fracTest"]*params["fracIso"]
  
  dt[, NonAdhereTransmissions := DailyTransmissions*(1-fracAdhere)]
  dt[, RemainingTransmissions := NonAdhereTransmissions + DailyTransmissions*(fracAdhere*(1-FracDiscovered))]
  
  dt[, Day:= Time/24]
  
  p =  ggplot(dt) + #geom_ribbon(aes(x = Time , ymin = RemainingTransmissions, ymax = HourlyTransmissions), fill = "grey", alpha = 0.6) + 
    geom_line(aes(x = Day , y = DailyTransmissions ),  linetype = "solid", colour = "black", linewidth = 1) +
    geom_ribbon(fill = "orange" , colour = "orange", alpha = 0.6 ,aes(x = Day ,ymin = 0,  ymax =  NonAdhereTransmissions)) + 
    geom_ribbon(fill = "purple" , colour = "purple", alpha = 0.6, aes(x = Day , ymin =  NonAdhereTransmissions, ymax = RemainingTransmissions)) + 
    xlab("Days Since Infection") + ylab("Expected Transmissions per Day")
  
  return(p)
}

generateReportFigures = function(){
  p = plotInfectiousness(c(contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,relativeDeclineSlope = 1.0, maxTimeAfterPeak = 24*30))
  ggsave("~/MassTesting/figures/infectiousness.pdf", p, width = 2.5, height = 2.5,device = "pdf")
  
  p = plotTestSensitivity(c(logLimitOfDetection = 3))
  ggsave("~/MassTesting/figures/testSensitivity.pdf", p, width = 2.5, height = 2.5,device = "pdf")
  
  
  
  plots = plot3Trajectories(c(  contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3, 
                      relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, 
                      logLimitOfDetection = 3, initialLogLoad = -2, precision = 0.25))
  ggsave("~/MassTesting/figures/viralLoad.pdf", plots[[1]], width = 4, height = 4,device = "pdf")
  ggsave("~/MassTesting/figures/testSensitivityVsTime.pdf", plots[[2]], width = 4, height = 4,device = "pdf")
  ggsave("~/MassTesting/figures/infectiousnessVsTime.pdf", plots[[3]], width = 4, height = 4,device = "pdf")
  
  
  p = plotFracReduction(testPeriods = 24*c(1), timesToPeak = 24*c(3,6,9), n = 20, showPooled = FALSE, params=  c(  contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3, 
                                                          relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logPeakLoad = 8,
                                                        initialLogLoad = -2, precision = 0.25))
  ggsave("~/MassTesting/figures/fracReduction2D.pdf", p, width = 10, height = 6,device = "pdf")
  
  p = plotFracReduction(testPeriods =24*c(1,2,4), timesToPeak = c(3,6,9), n = 10, showPooled = TRUE, params=  c(  contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3, 
                                                      relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logPeakLoad = 8,
                                                      initialLogLoad = -2, precision = 0.25))
  ggsave("~/MassTesting/figures/fracReduction2DSupplement.pdf", p, width = 10, height = 12,device = "pdf")
  
  
  p1 = plotEffectTestFreq(params=  c(  contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3, 
                                      relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logPeakLoad = 8,
                                      initialLogLoad = -2, precision = 2))
  
  p2 = plotEffectTestDelay(c(  contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3, 
                              relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logPeakLoad = 8,
                              initialLogLoad = -2, precision = 2))
  p = plot_grid(p2,p1, labels = c("A", "B"))
  
  ggsave("~/MassTesting/figures/effectTestFreqAndDelay.pdf", p, width = 10, height = 6,device = "pdf")
  
  
  # todo: add dashed line for antigen test and remove weekly testing. Use 8 hour delay for PCR
  # show example for high adherence and low adherence side by side
  # discuss effect of population-wide social distancing or masking as multiplying curves by multiple of 1/(1-lambda) (which for log scale looks like a vertical shift - shown in appendix)
  
  p1 = generate2TestControllabilityFigure(24*c(1,3), c( testDelay = 8, fracIso = 0.95, fracTest = 0.95, 
                                                    maskEffect = 0,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                                    relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, initialLogLoad = -2)) 

  
  p2 =  generate2TestControllabilityFigure(24*c(1,3), c( testDelay = 8, fracIso = 0.8, fracTest = 0.70, 
                                                         maskEffect = 0,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                                         relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, initialLogLoad = -2))
  p = plot_grid(p1+ theme(legend.position = "None"),p2+ theme(legend.position = "None"), labels = c("A", "B"))
  grobs <- ggplotGrob(p1)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  p = plot_grid(p, legend, nrow = 2, rel_heights = c(1, 0.05))
  
  ggsave("~/MassTesting/figures/controllabilityComparison.pdf", p, width = 14, height = 8,device = "pdf")
  
  p = generateControllabilityFigure(24*c(1,3,7), c( testDelay = 8, fracIso = 0.9, fracTest = 0.9, 
                                       maskEffect = 0,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                      relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logLimitOfDetection = 3, initialLogLoad = -2))
  ggsave("~/MassTesting/figures/controllability9090.pdf", p, width = 8, height = 8,device = "pdf")
  
  p = generateControllabilityFigure(24*c(1,3,7), c( testDelay = 8, fracIso = 0.95, fracTest = 0.95, 
                                                    maskEffect = 0.4,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                                    relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logLimitOfDetection = 3, initialLogLoad = -2))
  ggsave("~/MassTesting/figures/controllability9595.pdf", p, width = 8, height = 8,device = "pdf")
  
  p = generateControllabilityFigure(24*c(1,3,7), c( testDelay = 8, fracIso = 0.8, fracTest = 0.8, 
                                                    maskEffect = 0.2,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                                    relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logLimitOfDetection = 3, initialLogLoad = -2))
  ggsave("~/MassTesting/figures/controllability8080.pdf", p, width = 8, height = 8,device = "pdf")
  
  p = generateControllabilityFigure(24*c(1,3,7), c( testDelay = 24, fracIso = 0.7, fracTest = 0.5, 
                                                    maskEffect = 0.0,precision = 0.45, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                                                    relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, logLimitOfDetection = 3, initialLogLoad = -2))
  ggsave("~/MassTesting/figures/controllability5070.pdf", p, width = 8, height = 8,device = "pdf")
  
  # todo: change to daily cost per person (in dollars)
  p = plotPrevalenceCost(c(1,3,7), c(variableTestCost = 2, isolationCost = 5000, fixedAnnualizedDailyTestCost = 1))
  ggsave("~/MassTesting/figures/prevalenceCost.pdf", p, width = 7, height = 7,device = "pdf")
  
  
  p = plotMultiplePreventedTransmissions(c(contactsPerHour = 13/24, fracIso = 0.9, fracTest = 0.9, precision = 0.2,
                                                  maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30, 
                                                  logPeakLoad = 10, initialLogLoad = -2, logLimitOfDetection = 3, timeToPeak = 96, timeFromPeakTo0 = 96))
  ggsave("~/MassTesting/figures/preventedTransmissions.pdf", p, width = 10, height = 6,device = "pdf")
  # plot importation cost with vertical lines for UK and Australia?
  
  # plot sensitivity of result w.r.t distance between 50% infectious viral load and test limit of detection
  
  # ? plotOutbreak()
  
}


#inputParams = c(contactsPerHour = 13/24, testDelay = 24, fracIso = 0.9, fracTest = 0.9, precision = 0.15, maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30, logPeakLoad = 10, initialLogLoad = -2, logLimitOfDetection = 3)

#generateControllabilityFigure(c(24, 48), inputParams)
#plotTrajectories(inputParams)
#plotImportCost()

#plotFracTransmissionsAfterPositive(c(24, 48), inputParams)


#plotPreventedTransmissions(c(inputParams, c(testPeriod = 48, timeToPeak = 100, timeFromPeakTo0 = 96)))

