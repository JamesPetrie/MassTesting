
#test file for MT & CT
folder =  "~/MassTesting/ShinyApp/" # "/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/"
Rcpp::sourceCpp(paste0(folder, "ViralLoad.cpp"))
source(paste0(folder, "buildFigures.R"))
library(boot)
#plots (daily testing) another will have testing every 4 days for example ( noraml test period and outbreak test period)
#default 1 day and additional testing with 0 or 2 days ( tracing delay)
#fraction traced along x axis () ( probtraced given )
#prob detect symptoms (depends on disease) -> can set to 0 to ignore -> modify based on disease
#time to peak to symptoms -> modify based on disease
#days to peak -> depends on disease
#contact per day is 13, test delay = 10 hours, fraciso = 100 or 1, frac test = 0.9 to do thihnk about turn off contact tracing for people who dont test, 
#maxprobtranmit per exposure = 0.3, relative decline slope = 1, max days after peak. = 20 , mask effect = 0, minlogpcr = +2.5, initial log load = -2.5

# #base
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0 #adjust based on disease
# timeFromPeaktoSymptoms = 1 #adjust based on disease
# daysToPeak = 6 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 3 # adjust based on disease 

#example diseases
# #influenza
#   probTest = .50 #adjust scenario
#   normalTestPeriod = 4 #adjust scenario
#   outbreakTestPeriod = 4 #adjust scenario
#   tracingDelay = 0 #adjust scenario
#   #fractionTraced = 0.1 #adjust this one
#   probDetectSymptoms = 0.693*probTest #adjust based on disease
#   timeFromPeaktoSymptoms = -2 #adjust based on disease
#   daysToPeak = 3 #adjust based on disease
#   contactsPerDay = 13
#   testDelay = 10
#   fracIso = 1
#   fracTest = 0.90
#   maxProbTransmitPerExposure = 0.30
#   relativeDeclineSlope = 1
#   maxDaysAfterPeak = 20
#   maskEffect = 0
#   minLogPCRViralLoad = 2.5
#   initialLogLoad = -2.5
#   R0 = 2.0 #adjust based on disease

# #SARS-CoV
# probTest = .80 #adjust scenario
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0.925*probTest #adjust based on disease
# timeFromPeaktoSymptoms = -5 #adjust based on disease
# daysToPeak = 7 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 2.4 # adjust based on disease

#add covid-19 

minLogPCRViralLoad = 2.5
initialLogLoad = -2.5
maxProbTransmitPerExposure = 0.30
relativeDeclineSlope = 1
contactsPerDay = 13
precision = 2.0



meanfun <- function(data, i){
  d <- data[i, ]
  return(mean(d))   
}


test_period_list <- c(1,2,4,8,16)

disease_list <- c("1918 Influenza", "SARS-CoV-1","SARS-CoV-2")
#disease_list <- c("test")

fractionTraced_list <- seq(0, 0.9, by = 0.1)

timepeak_list <- seq(1,12,2)
prob_test_list <- c(0.10, 0.50, 0.80, 0.90) 
frac_iso_list <- c(.10,.50,.90) #remove 100% for figure B4

dt_final = rbindlist(llply(timepeak_list, function(time_peak){
  
  rbindlist(llply(prob_test_list, function(prob_test){

    rbindlist(llply(frac_iso_list, function(frac_iso){
    

    rbindlist(llply(test_period_list, function(test_period){
    
    normalTestPeriod = test_period
    outbreakTestPeriod = test_period
    
    rbindlist(llply(disease_list, function(disease_name){
      #normalTestPeriod = 6 #adjust scenario
      #outbreakTestPeriod = 6 #adjust scenario
      tracingDelay = 24 #adjust scenario
      #fractionTraced = 0.1 #adjust this one
      
      contactsPerDay = 13
      testDelay = 10
      #fracIso = 1
      fracIso = frac_iso
      fracTest = 0.90
      maxDaysAfterPeak = 20
      maskEffect = 0

      
      
      if(disease_name == "1918 Influenza"){
        #influenza
        #probTest = .50 #adjust scenario & disease
        probTest = prob_test #adjust scenario & disease
        probDetectSymptoms = 0.693*probTest #adjust based on disease
        timeFromPeaktoSymptoms = -2 #adjust based on disease
        daysToPeak = 3 #adjust based on disease
        R0 = 2.0 #adjust based on disease
      }
      else if(disease_name == "SARS-CoV-1") {
        #SARS-CoV
        #probTest = .80 #adjust scenario
        probTest = prob_test #adjust scenario
        probDetectSymptoms = 0.925*probTest #adjust based on disease
        timeFromPeaktoSymptoms = -5 #adjust based on disease
        daysToPeak = 7 #adjust based on disease
        R0 = 2.4 # adjust based on disease
      }else{
        #test with other values of time to peak
        probTest = .80 #adjust scenario
        probDetectSymptoms = 0.7*probTest #adjust based on disease
        timeFromPeaktoSymptoms = -1 #adjust based on disease
        daysToPeak = time_peak #adjust based on disease
        R0 = 2 # adjust based on disease
      }
      
      numOutbreaks = 800
      #numOutbreaks = 120
      endDay = 120
      maxSize = 300
      
      y = rbindlist(llply(fractionTraced_list, function(fractionTraced){
        params = c(
          normalTestPeriod = normalTestPeriod*24  ,
          outbreakTestPeriod = outbreakTestPeriod*24 ,
          ContactTracingDelay = tracingDelay ,
          ProbTracedGivenInfectorDetected = fractionTraced,
          ProbTestSymptoms = probDetectSymptoms, # todo: delete
          ProbDetectSymptoms = probDetectSymptoms,
          timeFromPeakToSymptoms = timeFromPeaktoSymptoms,
          timeToPeak = daysToPeak*24,
          maxTimeAfterPeak= 24*30,
          timeFromPeakTo0 = 24*5,
          contactsPerHour = contactsPerDay/24, testDelay = testDelay, fracIso = fracIso, fracTest = fracTest, 
          precision = 0.2, maxProbTransmitPerExposure = maxProbTransmitPerExposure,
          relativeDeclineSlope = relativeDeclineSlope, maxTimeAfterPeak = 24*maxDaysAfterPeak, 
          maskEffect = maskEffect, minLogPCRViralLoad = minLogPCRViralLoad, initialLogLoad = initialLogLoad, probTransmitMid = 8.9e6,logLimitOfDetection = 2.5, testPeriod = 1,fracQuar = 1 
          ,infectHParam = 0.51)
        
      
        #output final dt 
        params["logPeakLoad"] = computePeakViralLoad(params["timeToPeak"], targetR0 =R0, params)
        
        caseData = rbindlist(llply(1:numOutbreaks, function(i){
          dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
          dt[,RunNumber := i]
          
          return(dt)
        }))
        
        re_filt <- caseData[hourNotInfectious<SimulationEndHour]
        #for bootstrapping
        test <- data.frame(re_filt$NumInfected)
        bo <- boot(test[, "re_filt.NumInfected", drop = FALSE], statistic=meanfun, R=1000)
        
        result <- boot.ci(bo, conf=0.95, type = c("basic"))
        
        lowerBound <- result$basic[1,4]
        
        upperBound <- result$basic[1,5]
        #end bootstrapping
        
        re <- mean(re_filt$NumInfected)
        
        x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period,timePeak = time_peak, minBound = lowerBound, maxBound = upperBound,probTestPositive = prob_test,fractionIso = frac_iso)
        #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, probTestPositive = prob_test)
        #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, fractionIso = frac_iso)
        
        return(x)
      }))
    }))
  }))
}))
}))
}))




runAndComputeRe = function(numOutbreaks, endDay, maxSize, diseaseName, R0, timeToPeak, timeFromPeaktoSymptoms,  probTestSymptoms, 
                           testPeriod, tracingDelay, fractionTraced, testDelay, fracIso, fracTest ){
  params = c(
    normalTestPeriod = testPeriod,
    outbreakTestPeriod = testPeriod ,
    ContactTracingDelay = tracingDelay ,
    ProbTracedGivenInfectorDetected = fractionTraced,
    ProbTestSymptoms = probTestSymptoms,
    timeFromPeakToSymptoms = timeFromPeaktoSymptoms,
    timeToPeak = timeToPeak,
    timeFromPeakTo0 = timeToPeak, 
    #timeFromPeakTo0 = 24*5,
    maxTimeAfterPeak= 24*30,
    contactsPerHour = contactsPerDay/24, testDelay = testDelay, 
    fracIso = fracIso, fracQuar = fracIso, fracTest = fracTest, 
    maxProbTransmitPerExposure = maxProbTransmitPerExposure,
    #maxTimeAfterPeak = 24*30, 
    maxTimeAfterPeak = 24*20,
    maskEffect = 0, minLogPCRViralLoad = minLogPCRViralLoad, initialLogLoad = initialLogLoad, 
    probTransmitMid = 8.9e6,logLimitOfDetection = minLogPCRViralLoad ,infectHParam = 0.51, 
    #precision = precision
    precision = 0.2
    )
  

  #output final dt 
  params["logPeakLoad"] = computePeakViralLoad(params["timeToPeak"], targetR0 =R0, params)
  
  caseData = rbindlist(llply(1:numOutbreaks, function(i){
    #browser()
    dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
    dt[,RunNumber := i]
    
    return(dt)
  }))
  
  re_filt <- caseData[hourNotInfectious<SimulationEndHour]
  #for bootstrapping
  test <- data.frame(re_filt$NumInfected)
  bo <- boot(test[, "re_filt.NumInfected", drop = FALSE], statistic=meanfun, R=1000)
  
  result <- boot.ci(bo, conf=0.95, type = c("basic"))
  
  lowerBound <- result$basic[1,4]
  
  upperBound <- result$basic[1,5]
  #end bootstrapping
  
  re <- mean(re_filt$NumInfected)
  
  x = data.table(Re = re, Disease = diseaseName, FractionTraced = fractionTraced, R0 = R0, testPeriod = testPeriod,timePeak = timeToPeak, minBound = lowerBound, maxBound = upperBound,fractionIso = fracIso)
  #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, probTestPositive = prob_test)
  #x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period, fractionIso = frac_iso)
  
  return(x)
}


generateDiseaseDt = function(numOutbreaks = 2000, endDay = 120, maxSize = 300, probTestConditionalSymptomsInfluenza, probTestConditionalSymptomsSars,probTestConditionalSymptomsSars2 ){
  test_period_list <- c(1,2,4,8,16)*24
  fractionTraced_list <- seq(0, 0.9, by = 0.1)
  tracingDelay = 24
  testDelay = 10
  #fracIso = 0.9
  fracIso = 0.8
  fracTest = 0.9
  diseaseDt = rbindlist(llply(fractionTraced_list, function(fractionTraced){
    rbindlist(llply(test_period_list, function(testPeriod){
      rbindlist(llply(disease_list, function(disease_name){
        if(disease_name == "1918 Influenza"){
          probSymptoms = 0.693
          probTestSymptoms =   probSymptoms*probTestConditionalSymptomsInfluenza #adjust based on disease
          timeFromPeaktoSymptoms = -2*24  #-2*24 #adjust based on disease
          timeToPeak = 4*24 #adjust based on disease
          R0 = 2.0 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsInfluenza]
          return(dt)
          
        }
        else if(disease_name == "SARS-CoV-1") {
          probSymptoms = 0.925
          probTestSymptoms = probSymptoms*probTestConditionalSymptomsSars #adjust based on disease
          timeFromPeaktoSymptoms = -6*24 #adjust based on disease
          timeToPeak = 10.5*24 #adjust based on disease
          R0 = 2.4 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsSars]
          return(dt)
        }
        else if(disease_name == "SARS-CoV-2") {
          probSymptoms = 0.80 
          probTestSymptoms = probSymptoms*probTestConditionalSymptomsSars2 #adjust based on disease #find prob test conditional same as influenza 
          timeFromPeaktoSymptoms = -0*24 #adjust based on disease
          timeToPeak = 5*24 #adjust based on disease
          R0 = 2.5 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsSars2]
          return(dt)
        }
      }))
    }))
  }))
  return(diseaseDt)
}

# #figure 2
disease_dt <- generateDiseaseDt(numOutbreaks = 2000, probTestConditionalSymptomsInfluenza = 0.2, probTestConditionalSymptomsSars = 0.8,probTestConditionalSymptomsSars2 = 0.2)
disease_dt[, testPeriodDays := testPeriod/24]

a <- ggplot(disease_dt) + aes(
    x = FractionTraced,
    y = Re,
    colour = factor(testPeriodDays),
    group = testPeriodDays,
    ymin =  minBound,
    ymax = maxBound
    
  ) + geom_line(aes(colour = factor(testPeriodDays)))  +  scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = ""
  ) + geom_ribbon( aes(ymin =  minBound,ymax = maxBound,fill = factor(testPeriodDays) ),alpha = 0.3, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  facet_wrap(vars(Disease)) +
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse=TRUE) )

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_2.pdf", a, width = 7, height = 4,device = "pdf")
write.csv(disease_dt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_2.csv", row.names=TRUE)

# generateInfectiontoPeakDt = function(numOutbreaks = 1000, endDay = 120, maxSize = 300, probTestConditionalSymptoms = 0.5){
#   timepeak_list <- seq(1,12,1)*24
#   test_period_list <-  c(2,4,8)*24
#   #fractionTraced_list <- seq(0, 0.9, by = 0.1)
#   fractionTraced_list <- c(0, 0.5,0.9)
#   tracingDelay = 24
#   testDelay = 10
#   fracIso = 0.9
#   fracTest = 0.9
#   disease_name = "any"
#   diseaseDt = rbindlist(llply(fractionTraced_list, function(fractionTraced){
#     rbindlist(llply(test_period_list, function(testPeriod){
#       rbindlist(llply(timepeak_list, function(time_to_peak){
#         probSymptoms = 0.5 # todo: fix 0.70 #adjust scenario
#         probTestSymptoms = probSymptoms*probTestConditionalSymptoms #adjust based on disease
#         timeFromPeaktoSymptoms = -1*24 #adjust based on disease
#         timeToPeak = time_to_peak #adjust based on disease
#         R0 = 2 # adjust based on disease
#         dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
#                              probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
#         return(dt)
#       }))
#     }))
#   }))
#   return(diseaseDt)
# }
# 
# disease_dt <- generateInfectiontoPeakDt()
# disease_dt_TEST <- disease_dt[ FractionTraced == 0.5 | FractionTraced == 0 | FractionTraced == 0.9 ,]
# disease_dt_TEST[, testPeriod := as.character(testPeriod)][testPeriod == "48", testPeriod := "Test Period: 48 hr"]
# disease_dt_TEST[, testPeriod := as.character(testPeriod)][testPeriod == "96", testPeriod := "Test Period: 96 hr"]
# disease_dt_TEST[, testPeriod := as.character(testPeriod)][testPeriod == "192", testPeriod := "Test Period: 192 hr"]
# 
# #figure 3 
# b <- ggplot(disease_dt_TEST) +
#   aes(
#     x = timePeak/24,
#     y = Re,
#     colour = factor(FractionTraced),
#     group = FractionTraced,
#     fill = factor(FractionTraced)
#   ) +
#   geom_line(aes(colour = factor(FractionTraced))) + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of Contacts Traced") + 
#   scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of Contacts Traced")+
#   labs(
#     x = "Time from Infection to Peak Viral Load (days)",
#     y = "Effective Reproduction Number (Re)",
#     title = "Re vs Time from Infection to Peak Viral Load"
#   ) +
#   theme_minimal() +
#   geom_ribbon( aes(ymin =  minBound,ymax = maxBound,fill = factor(FractionTraced) ),alpha = 0.3, colour = NA )+ theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7)) + 
#   facet_wrap(~factor(testPeriod, levels=c("Test Period: 48 hr", "Test Period: 96 hr", "Test Period: 192 hr")))
# 
# ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_infection_peak_2.pdf", b, width = 7, height = 4,device = "pdf")

#figure 3 edited 
generateTestPeriodDt = function(numOutbreaks = 2000, endDay = 120, maxSize = 300, probTestConditionalSymptoms = 0.80){
  timepeak_list <- c(2,4,6,8,10,12)*24
  test_period_list <- floor(24*2^seq(0, 5, length.out = 12))
  #fractionTraced_list <- seq(0, 0.9, by = 0.1)
  fractionTraced_list <- c(0, 0.5,0.9)
  tracingDelay = 24
  testDelay = 10
  fracIso = 0.9
  fracTest = 0.9
  disease_name = "any"
  diseaseDt = rbindlist(llply(fractionTraced_list, function(fractionTraced){
    rbindlist(llply(test_period_list, function(testPeriod){
      rbindlist(llply(timepeak_list, function(time_to_peak){
        probSymptoms = 0.50 #adjust scenario
        probTestSymptoms = probSymptoms*probTestConditionalSymptoms #adjust based on disease
        timeFromPeaktoSymptoms = -1*24 #adjust based on disease
        timeToPeak = time_to_peak #adjust based on disease
        R0 = 2 # adjust based on disease
        dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                             probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
        return(dt)
      }))
    }))
  }))
  return(diseaseDt)
}


testperiod_dt = generateTestPeriodDt(numOutbreaks = 2000, probTestConditionalSymptoms = 0.5)
c <- ggplot(testperiod_dt) +
  aes(
    x = 24/testPeriod,
    y = Re,
    colour = factor(FractionTraced),
    group = FractionTraced,
    fill = factor(FractionTraced)
  ) +
  geom_line() + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of\nContacts Traced") + 
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of\nContacts Traced") + 
  labs(
    x = "Tests Per Day",
    y = "Effective Reproduction Number (Re)",
    title = ""
  ) + geom_ribbon( aes(ymin =minBound,ymax = maxBound),alpha = 0.09, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7)) + 
  facet_wrap(~factor(paste0("Days to peak viral load: " , timePeak/24), levels=c("Days to peak viral load: 2","Days to peak viral load: 4","Days to peak viral load: 6",
                                                                                 "Days to peak viral load: 8","Days to peak viral load: 10","Days to peak viral load: 12")))+ 
  scale_x_log10(breaks = c(1/32, 1/16, 1/8, 1/4, 1/2, 1), labels= c("1/32","1/16", "1/8", "1/4", "1/2", "1"))

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_infection_peak_3.pdf", c, width = 7, height = 4,device = "pdf")
write.csv(testperiod_dt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_infection_peak_3.csv", row.names=TRUE)

#figure A1: 
disease_dt <- generateDiseaseDt(numOutbreaks = 2000, probTestConditionalSymptomsInfluenza = 0.1, probTestConditionalSymptomsSars = 0.1,probTestConditionalSymptomsSars2 = 0.1)
disease_dt[, testPeriodDays := testPeriod/24]

a1 <- ggplot(disease_dt) + aes(
  x = FractionTraced,
  y = Re,
  colour = factor(testPeriodDays),
  group = testPeriodDays,
  ymin =  minBound,
  ymax = maxBound
  
) + geom_line(aes(colour = factor(testPeriodDays)))  +  scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = ""
  ) + geom_ribbon( aes(ymin =  minBound,ymax = maxBound,fill = factor(testPeriodDays) ),alpha = 0.3, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  facet_wrap(vars(Disease))+
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse=TRUE) )

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_10_2.pdf", a1, width = 7, height = 4,device = "pdf")
write.csv(disease_dt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_10_2.csv", row.names=TRUE)

#figure A2:
disease_dt <- generateDiseaseDt(numOutbreaks = 1000, probTestConditionalSymptomsInfluenza = 0.5, probTestConditionalSymptomsSars = 0.5,probTestConditionalSymptomsSars2 = 0.5)
disease_dt[, testPeriodDays := testPeriod/24]

a2 <- ggplot(disease_dt) + aes(
  x = FractionTraced,
  y = Re,
  colour = factor(testPeriodDays),
  group = testPeriodDays,
  ymin =  minBound,
  ymax = maxBound
  
) + geom_line(aes(colour = factor(testPeriodDays)))  +  scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = ""
  ) + geom_ribbon( aes(ymin =  minBound,ymax = maxBound,fill = factor(testPeriodDays) ),alpha = 0.3, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  facet_wrap(vars(Disease))+
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse=TRUE) )

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_50_2.pdf", a2, width = 7, height = 4,device = "pdf")
write.csv(disease_dt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_50_2.csv", row.names=TRUE)

#figure A3: 
disease_dt <- generateDiseaseDt(numOutbreaks = 1000, probTestConditionalSymptomsInfluenza = 0.9, probTestConditionalSymptomsSars = 0.9,probTestConditionalSymptomsSars2 = 0.9)
disease_dt[, testPeriodDays := testPeriod/24]

a3 <- ggplot(disease_dt) + aes(
  x = FractionTraced,
  y = Re,
  colour = factor(testPeriodDays),
  group = testPeriodDays,
  ymin =  minBound,
  ymax = maxBound
  
) + geom_line(aes(colour = factor(testPeriodDays)))  +  scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (days)") +
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = ""
  ) + geom_ribbon( aes(ymin =  minBound,ymax = maxBound,fill = factor(testPeriodDays) ),alpha = 0.3, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  facet_wrap(vars(Disease))+
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse=TRUE) )

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_90_2.pdf", a3, width = 7, height = 4,device = "pdf")
write.csv(disease_dt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_90_2.csv", row.names=TRUE)

#figure B4
generateEffectiveIsoDt = function(numOutbreaks = 1000, endDay = 120, maxSize = 300, probTestConditionalSymptomsInfluenza = 0.2, probTestConditionalSymptomsSars = 0.8 ,probTestConditionalSymptomsSars2 = 0.2){
  test_period_list <- c(2,8)*24
  fractionTraced_list <- c(0.10, 0.50)#seq(0, 0.9, by = 0.1)
  tracingDelay = 24
  testDelay = 10
  #fracIso = 0.9
  frac_iso_list <- c(.10,.50,.90)
  fracTest = 0.9
  diseaseDt = rbindlist(llply(fractionTraced_list, function(fractionTraced){
    rbindlist(llply(test_period_list, function(testPeriod){
      rbindlist(llply(frac_iso_list, function(fracIso){
      rbindlist(llply(disease_list, function(disease_name){
        if(disease_name == "1918 Influenza"){
          probSymptoms = 0.693
          probTestSymptoms =   probSymptoms*probTestConditionalSymptomsInfluenza #adjust based on disease
          timeFromPeaktoSymptoms = -2*24  #-2*24 #adjust based on disease
          timeToPeak = 4*24 #adjust based on disease
          R0 = 2.0 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsInfluenza]
          return(dt)
          
        }
        else if(disease_name == "SARS-CoV-1") {
          probSymptoms = 0.925
          probTestSymptoms = probSymptoms*probTestConditionalSymptomsSars #adjust based on disease
          timeFromPeaktoSymptoms = -6*24 #adjust based on disease
          timeToPeak = 10.5*24 #adjust based on disease
          R0 = 2.4 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsSars2]
          return(dt)
        }
        else if(disease_name == "SARS-CoV-2") {
          probSymptoms = 0.80 
          probTestSymptoms = probSymptoms*probTestConditionalSymptomsSars2 #adjust based on disease #find prob test conditional same as influenza 
          timeFromPeaktoSymptoms = -0*24 #adjust based on disease
          timeToPeak = 5*24 #adjust based on disease
          R0 = 2.5 #adjust based on disease
          dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                               probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
          dt[, probTestConditionalSymptoms := probTestConditionalSymptomsSars]
          return(dt)
        }
      }))
    }))
    }))
  }))
  return(diseaseDt)
}
testdt <- generateEffectiveIsoDt()
d <- ggplot(data = testdt)+
  geom_line(data = testdt[FractionTraced == 0.10 & testPeriod == 48, ], aes(y = Re, x = fractionIso, color ="#66C2A5" )) +   
  geom_ribbon( data = testdt[FractionTraced == 0.10 & testPeriod == 48, ],aes(y = Re, x=  fractionIso,ymin =  minBound,ymax = maxBound,fill = factor(testPeriod) ),alpha = 0.3, fill = "#66C2A5" )+
  geom_line(data = testdt[FractionTraced == 0.50 & testPeriod == 48 ,], aes(y = Re, x=  fractionIso, color ="#FC8D62"))  +
  geom_ribbon( data = testdt[FractionTraced == 0.50 & testPeriod == 48 ,],aes(y = Re, x=  fractionIso,ymin =  minBound,ymax = maxBound,fill = factor(testPeriod) ),alpha = 0.3, fill = "#FC8D62" )+
  geom_line(data = testdt[FractionTraced == 0.10 & testPeriod == 192 ,], aes(y = Re, x = fractionIso, color ="#8DA0CB"))  + 
  geom_ribbon( data = testdt[FractionTraced == 0.10 & testPeriod == 192 ,],aes(y = Re, x=  fractionIso,ymin =  minBound,ymax = maxBound,fill = factor(testPeriod) ),alpha = 0.3, fill = "#8DA0CB" )+
  geom_line(data = testdt[FractionTraced == 0.50 & testPeriod == 192 ,], aes(y = Re, x=  fractionIso, color ="#E78AC3")) +
  geom_ribbon( data = testdt[FractionTraced == 0.50 & testPeriod == 192 ,],aes(y = Re, x=  fractionIso,ymin =  minBound,ymax = maxBound,fill = factor(testPeriod) ),alpha = 0.3, fill = "#E78AC3" )+
  facet_wrap(~Disease)+theme_minimal() + scale_color_identity(name = "Test Period (days) & Proportion of Contacts Traced",
                                                              breaks = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
                                                              labels = c("Period: 2, Traced: 0.10", "Period: 2, Traced: 0.50", "Period: 8, Traced: 0.10", "Period: 8, Traced: 0.50"),
                                                              guide = "legend") + scale_fill_identity(name = "Test Period (days) & Proportion of Contacts Traced",
                                                                                                      breaks = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
                                                                                                      labels = c("Period: 2, Traced: 0.10", "Period: 2, Traced: 0.50", "Period: 8, Traced: 0.10", "Period: 8, Traced: 0.50"),
                                                                                                      guide = "legend")+
  labs(title = "", y =("Effective Reproduction Number (Re)"), x = ("Effectiveness of Isolation/Quarantine")) + guides(color=guide_legend(title="Test Period (days) & Proportion of Contacts Traced")) + 
  theme(legend.title=element_text(size=6),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=6))

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_prop_contacts_isolated_2.pdf", d, width = 7, height = 4,device = "pdf")
write.csv(testdt, "/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_prop_contacts_isolated_2.csv", row.names=TRUE)






generateTestPeriodDt = function(numOutbreaks = 100, endDay = 120, maxSize = 300, probTestConditionalSymptoms = 0.80){
  timepeak_list <- c(2,4,6,8,10,12)*24
  test_period_list <- floor(24*2^seq(0, 5, length.out = 12))
  #fractionTraced_list <- seq(0, 0.9, by = 0.1)
  fractionTraced_list <- c(0, 0.5,0.9)
  tracingDelay = 24
  testDelay = 10
  fracIso = 0.9
  fracTest = 0.9
  disease_name = "any"
  diseaseDt = rbindlist(llply(fractionTraced_list, function(fractionTraced){
    rbindlist(llply(test_period_list, function(testPeriod){
      rbindlist(llply(timepeak_list, function(time_to_peak){
        probSymptoms = 0.50 #adjust scenario
        probTestSymptoms = probSymptoms*probTestConditionalSymptoms #adjust based on disease
        timeFromPeaktoSymptoms = -1*24 #adjust based on disease
        timeToPeak = time_to_peak #adjust based on disease
        R0 = 2 # adjust based on disease
        dt = runAndComputeRe(numOutbreaks, endDay, maxSize, disease_name, R0, timeToPeak, timeFromPeaktoSymptoms, 
                             probTestSymptoms, testPeriod, tracingDelay,fractionTraced, testDelay, fracIso, fracTest )
        return(dt)
      }))
    }))
  }))
  return(diseaseDt)
}


testperiod_dt = generateTestPeriodDt(numOutbreaks = 20, probTestConditionalSymptoms = 0.5)
c <- ggplot(testperiod_dt) +
  aes(
    x = 24/testPeriod,
    y = Re,
    colour = factor(FractionTraced),
    group = FractionTraced,
    fill = factor(FractionTraced)
  ) +
  geom_line() + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of Contacts Traced") + 
  labs(
    x = "Time from Infection to Peak Viral Load (days)",
    y = "Effective Reproduction Number (Re)",
    title = "Re vs Time from Infection to Peak Viral Load"
  ) + geom_ribbon( aes(ymin =  minBound,ymax = maxBound),alpha = 0.09, colour = NA )+ 
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7)) + 
  facet_wrap(~factor(paste0("Days to peak viral load: " , timePeak/24)))+ 
  scale_x_log10(breaks = c(1/32, 1/16, 1/8, 1/4, 1/2, 1), labels= c("1/32","1/16", "1/8", "1/4", "1/2", "1"))

# p <- ggplot(data = dt_final[Disease == "1918 Influenza"]) + geom_line(aes(x = FractionTraced, y = Re),colour = "#00A2FF", size = 0.5)+  theme_minimal() + labs(y = "Effective Reproduction Number (Re)", x = "Proportion of Contacts Traced", title = "Re vs Proportion of Contacts Traced for 1918 Influenza Pandemic")
# 
# #ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_influenza.pdf", p, width = 7, height = 4,device = "pdf")
# 
# p2 <- ggplot(data = dt_final[Disease == "SARS-CoV-1"]) + geom_line(aes(x = FractionTraced, y = Re),colour = "#00A2FF", size = 0.5)+  theme_minimal() + labs(y = "Effective Reproduction Number (Re)", x = "Proportion of Contacts Traced", title = "Re vs Proportion of Contacts Traced for SARS-CoV-1")
# 
# #ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_SARS.pdf", p2, width = 7, height = 4,device = "pdf")
# 
# 
# a <- ggplot(dt_final) +
#   aes(
#     x = FractionTraced,
#     y = Re,
#     colour = testPeriod,
#     group = testPeriod
#   ) +
#   geom_line() +
#   scale_color_distiller(palette = "Set2", direction = 1) +
#   labs(
#     x = "Proportion of Contacts Traced",
#     y = "Effective Reproduction Number (Re)",
#     title = "Re vs Proportion of Contacts Traced"
#   ) +
#   theme_minimal() +
#   facet_wrap(vars(Disease))

# a <- ggplot(dt_final) +
#   aes(
#     x = FractionTraced,
#     y = Re,
#     colour = factor(testPeriod),
#     group = testPeriod,
#     fill = factor(testPeriod)
#   ) +
#   geom_line() + scale_colour_discrete(name = "testPriod") +
#   labs(
#     x = "x lab title ",
#     y = "y lab title",
#     title = "title main"
#   ) +
#   theme_minimal() +
#   facet_wrap(vars(Disease))

#edit filter for real values probTest = .50 and 0.80 for influenza and sars 

# # #figure 2
# a <- ggplot(dt_final) +
#   aes(
#     x = FractionTraced,
#     y = Re,
#     colour = factor(testPeriod),
#     group = testPeriod,
#     fill = factor(testPeriod)
#   ) +
#   geom_line() + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (Days)") +
#   labs(
#     x = "Proportion of Contacts Traced",
#     y = "Effective Reproduction Number (Re)",
#     title = "Re vs Proportion of Contacts Traced"
#   ) +
#   theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+
#   facet_wrap(vars(Disease))
# 
# ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_2.pdf", a, width = 7, height = 4,device = "pdf")
# 
# dt_filt <- dt_final[ FractionTraced == 0.5 | FractionTraced == 0 | FractionTraced == 0.9 ,]
# dt_filt <- dt_filt[testPeriod == 4,]
# 
# #figure 3
# b <- ggplot(dt_filt) +
#   aes(
#     x = timePeak,
#     y = Re,
#     colour = factor(FractionTraced),
#     group = FractionTraced,
#     fill = factor(FractionTraced)
#   ) +
#   geom_line() + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(11), name = "Proportion of Contacts Traced") + 
#   labs(
#     x = "Time from Infection to Peak Viral Load (days)",
#     y = "Effective Reproduction Number (Re)",
#     title = "Re vs Time from Infection to Peak Viral Load"
#   ) +
#   theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))
# 
# ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_infection_peak_2.pdf", b, width = 7, height = 4,device = "pdf")
# 
# dt_10 <- dt_final[probTestPositive == 0.1,]
# dt_50 <- dt_final[probTestPositive == 0.5,]
# dt_90 <- dt_final[probTestPositive == 0.9,]
# 
# #figure A1 to A3: 
# c <- ggplot(dt_90) +
#   aes(
#     x = FractionTraced,
#     y = Re,
#     colour = factor(testPeriod),
#     group = testPeriod,
#     fill = factor(testPeriod)
#   ) +
#   geom_line() + scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(11), name = "Test Period (Days)") + 
#   labs(
#     x = "Proportion of Contacts Traced",
#     y = "Effective Reproduction Number (Re)",
#     title = "Re vs Proportion of Contacts Traced"
#   ) +
#   theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+ 
#   facet_wrap(vars(Disease)) 
# 
# #figure B4
# ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period_probTest_90.pdf", c, width = 7, height = 4,device = "pdf")
# 
# dt_tst_fraciso <- dt_final[testPeriod == 2 | testPeriod == 8,]
# dt_tst_fraciso <- dt_tst_fraciso[FractionTraced == 0.10 | FractionTraced ==.50,]
# 
# d <- ggplot(data = dt_tst_fraciso)+
#   geom_line(data = dt_tst_fraciso[FractionTraced == 0.10 & testPeriod == 2, ], aes(y = Re, x = fractionIso, color ="#66C2A5" )) +
#   geom_line(data = dt_tst_fraciso[FractionTraced == 0.50 & testPeriod == 2 ,], aes(y = Re, x=  fractionIso, color ="#FC8D62"))  +
#   geom_line(data = dt_tst_fraciso[FractionTraced == 0.10 & testPeriod == 8 ,], aes(y = Re, x = fractionIso, color ="#8DA0CB"))  + 
#   geom_line(data = dt_tst_fraciso[FractionTraced == 0.50 & testPeriod == 8 ,], aes(y = Re, x=  fractionIso, color ="#E78AC3")) + 
#   facet_wrap(~Disease)+theme_minimal() + scale_color_identity(name = "Test Period (days) & Proportion of Contacts Traced",
#                                                               breaks = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
#                                                               labels = c("Period: 2, Traced: 0.10", "Period: 2, Traced: 0.50", "Period: 8, Traced: 0.10", "Period: 8, Traced: 0.50"),
#                                                               guide = "legend") +
#   labs(title = "Re vs Effectiveness of Isolation/Quarantine", y =("Effective Reproduction Number (Re)"), x = ("Effectiveness of Isolation/Quarantine")) + guides(color=guide_legend(title="Test Period (days) & Proportion of Contacts Traced")) + theme(legend.title=element_text(size=6),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=6))
# 
# ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_prop_contacts_isolated.pdf", d, width = 7, height = 4,device = "pdf")
