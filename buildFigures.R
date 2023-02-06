require(ggforce)
require(ggplot2)
require(cowplot)

source("~/MassTesting/ViralLoad.R")

theme_set(theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) + theme(legend.background = element_rect(fill = "white")) + theme_half_open() + background_grid()  + theme(text = element_text(size=20), axis.text = element_text(size=20)))


generateControllabilityFigure = function(testPeriods, inputParams){
  dt = rbindlist(llply(testPeriods, function(testPeriod){
    params = copy(inputParams)
    params["testPeriod"] = testPeriod
    evaluateStrategy(params)
  }))
  if(nrow(dt)>0){
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    dt[,DaysToPeak := TimeToPeak/24]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wuhan)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
      DaysToPeak = c(8, 5.5, 4, 10, 10), R0 = c(2.5, 8, 2.5, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak +0.5*c(1,0,-1,0), R0 = R0 + 0.5*c(0, 1,0,-1)), by = Pathogen ]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = FreqLabel), linewidth = 1.4) + scale_y_sqrt(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(1,20)) + scale_x_continuous(breaks = 0:12) + 
      guides(colour=guide_legend(title="Maximum Controllable\nR0 with Test Frequency \n [1/Days]")) + xlab("Time to Peak Viral Load [Days]") + ylab("R0")+
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, fill = Pathogen, label = Pathogen), 
                        linewidth = 0.0, label.fontsize = 8)
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

plotViralLoads = function(params){
  params["logPeakLoad"] = 10
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(3,7,10), Time = 24*seq(0,15, length.out = 100))) 
  dt[ ,ViralLoad := computeViralLoad(Time, replacePeakTime(params, TimeToPeak)), by = list(Time, TimeToPeak)]
  ggplot(dt, aes(x = Time/24, y = ViralLoad, group = TimeToPeak)) + geom_line() + scale_y_log10() + labs(y = "Log Viral Load", x = "Days Since Infection") + 
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
}




inputParams = c(contactsPerHour = 13/24, testDelay = 12, fracIso = 0.9, fracTest = 0.9, precision = 0.15, maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30)

#generateControllabilityFigure(c(24, 48), inputParams)
plotViralLoads(inputParams)

