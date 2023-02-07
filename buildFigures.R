require(ggforce)
require(ggplot2)
require(cowplot)
require(scales)

source("~/MassTesting/ViralLoad.R")

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
      labs(title = "Effect of Mass Testing",  subtitle = "Maximum Controllable R0 for different testing strategies")
      #geom_ellipse(data = data.table(), aes(x0 = 5, y0 = 3, a = 1, b = 1, angle = 0), fill = "orange", alpha = 0.4)
    
    #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
  }
  return(p)
}

replacePeakTime = function(params, timeToPeak){
  newParams = copy(params)
  newParams["timeToPeak"] = timeToPeak
  newParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
  return(newParams)
}

plotTrajectories = function(params){
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(3,7,10), Time = 24*seq(0,16, length.out = 100))) 
  dt[ ,ViralLoad := computeViralLoad(Time, replacePeakTime(params, TimeToPeak)), by = list(Time, TimeToPeak)]
  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*probTransmit(ViralLoad, replacePeakTime(params, TimeToPeak)), by = list( ViralLoad, TimeToPeak) ]
  dt[ ,TestSensitivity :=probPositive(ViralLoad, replacePeakTime(params, TimeToPeak)), by = list( ViralLoad, TimeToPeak) ]
  
  
  peakNames = dt[, list(PeakLabel = paste("Days to Peak: ", TimeToPeak/24)), by = TimeToPeak]
  setkey(peakNames, by = "TimeToPeak")
  peakNames[, PeakLabel := factor(PeakLabel, levels = PeakLabel)]
  dt = merge(dt, peakNames, by = "TimeToPeak")
  
  
  dt[, LogViralLoad := log10(ViralLoad)]
  dtLong = melt(dt, id.vars = c("Time", "PeakLabel"), measure.vars = c("LogViralLoad", "TestSensitivity","DailyTransmissions"))
  dtLong[, value:=value/max(value), by = variable]
  dtLong[variable == "LogViralLoad", variable := "Log Viral Load"]
  dtLong[variable == "TestSensitivity", variable := "Test Sensitivity"]
  dtLong[variable == "DailyTransmissions", variable := "Daily Transmissions"]
  
  ggplot(dtLong, aes(x = Time/24, y = value)) + facet_wrap( ~   PeakLabel + variable , nrow = 3, dir = "v") + geom_line() + 
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.y = element_blank()) + xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank()) + ggtitle("Example viral load trajectory, test sensitivity, expected transmissions")
}

plotInfectiousness = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 15, length.out = 100))
  dt[ ,ProbTransmit :=probTransmit(ViralLoad, params), by = ViralLoad ]
  
  ggplot(dt , aes(x = ViralLoad, y = ProbTransmit)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + xlab("Viral Load [RNA copies/ml]") + ylab("Transmission\n Probability") +
    theme(text = element_text(size=12), axis.text = element_text(size=12)) + labs(title="Probability of Transmission per Exposure")
  
}

plotTestSensitivity = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 15, length.out = 100))
  dt[ ,ProbPositive :=probPositive(ViralLoad, params), by = ViralLoad ]
  
  ggplot(dt , aes(x = ViralLoad, y = ProbPositive)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + xlab("Viral Load [RNA copies/ml]") + ylab("Test\n Sensitivity") + 
    theme(text = element_text(size=12), axis.text = element_text(size=12)) + labs(title="PCR Test Sensitivity")
  
}






inputParams = c(contactsPerHour = 13/24, testDelay = 12, fracIso = 0.9, fracTest = 0.9, precision = 0.15, maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30, logPeakLoad = 10)

#generateControllabilityFigure(c(24, 48), inputParams)
plotTrajectories(inputParams)

