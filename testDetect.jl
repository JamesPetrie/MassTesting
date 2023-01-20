
using StatsFuns
using DataFrames
using CSV
using OrderedCollections
using Distributions
using Dates


include("/Users/jpetrie/OptimalQuarantine/VariableQuarantine/SimpleQuarantine.jl")

   


testParamsFile = "/Users/jpetrie/OptimalQuarantine/VariableQuarantine/testSensitivityParams2023Jan13.csv"
testParams = DataFrame(CSV.File(testParamsFile))
testParams = testParams[:, 2:end]  


dist = Distributions.LogNormal(1.58, 0.47)
incubationProb = [Distributions.pdf(dist, x) for x = 0:30]
incubationProb[1] = 0

baseParams = OrderedCollections.OrderedDict{Symbol,Any}(
  :pAsymp => 0.2,
  :Rsymp => 1.0,
  :Rasymp => 1.0,
  :incubationProb => incubationProb/sum(incubationProb),
  :scaleParams => [1.0, 1.0],
  :transParams => [-4, 1.85, 5.85, 5.42],
  :infDay => 0,
)
params = deepcopy(baseParams)

transFunc = createTransmissionFunc(params[:infDay], params[:transParams], params[:incubationProb], params[:Rsymp], params[:Rasymp])
avTestSens =  createAvTestSensFunc(testParams, params)


function integrateTransmissions(startTime, incub, asymp)
  if(startTime > 30)
    return(0)
  end

  transmissions = 0
  for day = startTime:30
    transmissions = transmissions + transFunc(day, incub, asymp)
  end
  
  return(transmissions)
end


# average over timing offset, incubation, asymp, time of first positive
# return expected fraction of transmissions occuring after receiving positive result
function fractionAfterPositive(testPeriod, testDelay)
  if testPeriod <= 0
    return(1)
  end

  fracTrans = 0.0 
  for offset in 0:(testPeriod -1)
    for incub in 0:30
      for asymp in [true, false]
        remainingProb = 1/testPeriod # because sampling over this many timing offsets
        remainingProb = remainingProb*incubationProb[incub+1]
        if asymp
          remainingProb = remainingProb*params[:pAsymp]
        else
          remainingProb = remainingProb*(1- params[:pAsymp])
        end
        testTime = offset 
        while testTime <= 30
          probPos = avTestSens(testTime, 0, incub)

          fracTrans = fracTrans + remainingProb*probPos*integrateTransmissions(testTime + testDelay, incub, asymp) 

          remainingProb = remainingProb*(1-probPos)
          testTime = testTime + testPeriod
        end
      end
    end
  end
  return(fracTrans)
end


# compute fractionAfterPositive for a grid of testDelay, testPeriod values
# write results to csv

df = DataFrame(TestPeriod = Number[], TestDelay = Number[], FracAfterPositive = Number[])

for testPeriod in 1:10 
  for testDelay in 0:1:4 
    fracAfter = fractionAfterPositive(testPeriod, testDelay)
    push!(df, [testPeriod  testDelay fracAfter]) 
  end
end

print(df)
CSV.write("/Users/jpetrie/MassTesting/FracTransmissions.csv", df)
