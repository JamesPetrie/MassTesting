
using StatsFuns
using DataFrames
using CSV
using OrderedCollections
using Distributions
using Dates


include("~/OptimalQuarantine/VariableQuarantine/SimpleQuarantine.jl")

   


  testParamsFile = "~/OptimalQuarantine/VariableQuarantine/testSensitivityParams2023Jan13.csv"
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
    #:testParams => testParams, #c(0.2, 1, 1.7, 4.7), # new form
    :scaleParams => [1.0, 1.0],
    #scaleParams => [0.5, 1.3],
    :transParams => [-4, 1.85, 5.85, 5.42],
    :infDay => 0,
    :negTestEffect => 0.5

  )
  params = deepcopy(baseParams)

  transFunc = createTransmissionFunc(params[:infDay], params[:transParams], params[:incubationProb], params[:Rsymp], params[:Rasymp])
  avTestSens =  createAvTestSensFunc(testParams, params)


