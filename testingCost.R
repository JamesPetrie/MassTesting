R0 = 3.0


computeCost = function(testFraction = 0.95, isoFraction = 0.9, testFreq = 0.5){
  transFraction = 0.8 # todo: make function of test frequency, delay, sensitivity, expected transmission function
  importRate = 1e-4
  
  
  testCost = 10 + 40
  isoCost = 10*(100 + 200) + 0.5*3000
  
  Rt = R0*(1 - testFraction*transFraction*isoFraction)
  
  dailyCostPerPerson = testFraction*testFreq*testCost + importRate*1/(1-Rt)*testFraction*isoCost
  gdpPerCapita = 70e3 # gdp per person in USA
  return(dailyCostPerPerson*365/gdpPerCapita)
}



# todo: compute fraction of transmissions after receiving positive test as function of frequency and delay for Covid


# number of importations essen