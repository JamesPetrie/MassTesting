





# computes probability of infection for a typical contact given viral load
probTransmit = function(viralLoad, params = NULL){

  return((viralLoad > 1)*(1 - exp(-0.2*(viralLoad^0.51)/(viralLoad^0.51 + (8.9e6)^0.51))))
 }





# computes probability of a positive PCR result given viral load
probPositive = function(viralLoad, params = NULL){
  if(viralLoad > 1e3){return(1)}
  else{return(0)}
}

# compute viral load at a given time relative to infection
# todo: decide if 10^0 is right intercept
computeViralLoad = function(time, viralParams){
  # if(time < 0){
  #   logLoad = 0
  # }else if(time <= viralParams$timeToPeak){
  #   logLoad = viralParams$logPeakLoad*(time/viralParams$timeToPeak)
  # }else if(time <= viralParams$timeToPeak + viralParams$timeFromPeakTo0){
  #   logLoad = viralParams$logPeakLoad*(1 - ((time-viralParams$timeToPeak)/viralParams$timeFromPeakTo0))
  # }else{
  #   logLoad = 0
  # }
  
  logVals = (time >= 0)*(time <= viralParams$timeToPeak)*viralParams$logPeakLoad*(time/viralParams$timeToPeak) +
    (time > viralParams$timeToPeak)*(time <= viralParams$timeToPeak + viralParams$timeFromPeakTo0)*viralParams$logPeakLoad*(1 - ((time-viralParams$timeToPeak)/viralParams$timeFromPeakTo0))
  return(10^logVals)
}


# computes expected transmissions over specified hour range relative to day of infection
sumTransmissions = function(startTime, endTime, params){
  val = integrate( lower = startTime, upper = endTime, f = function(tVals){
    viralLoads = computeViralLoad(tVals, params)
    return(params$contactsPerHour*probTransmit(viralLoads, params))
  })
  return(val$value)
}


#params = list(timeToPeak = 6*24, timeFromPeakTo0 = 8*24, logPeakLoad = 9, contactsPerHour = 13/24)

#sumTransmissions(0,24*12, params)


# computes the fraction of transmissions that occur after receiving a positive test result
fracAfterPositive = function(params){
  # average over test timing offsets
  
  testPeriod = params$testPeriod
  testDelay = params$testDelay
  stopTime = params$timeToPeak + params$timeFromPeakTo0
  
  if(testPeriod <= 0) return(NA)
  
  # todo: fix so doesn't use same point twice (0 and testPeriod)
  transmissionsAfter = mean(sapply(seq(0, testPeriod, length.out = 30), function(offset){
    remainingProb = 1.0
    expectedTransmissions = 0.0
    testTime = offset
    while(testTime <= stopTime ){
      probPos = probPositive(computeViralLoad(testTime, params), params)
      
      expectedTransmissions = expectedTransmissions + remainingProb*probPos*sumTransmissions(testTime + testDelay, stopTime, params)
      
      remainingProb = remainingProb*(1-probPos)
      
      testTime = testTime + testPeriod
    }
    return(expectedTransmissions)
  }))
  
  return(transmissionsAfter/sumTransmissions(0, stopTime, params))
}

params = list(timeToPeak = 4*24, timeFromPeakTo0 = 6*24, logPeakLoad = 11, contactsPerHour = 13/24, testPeriod = 1*24, testDelay = 12)
fracAfterPositive(params)
