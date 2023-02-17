#include <Rcpp.h>
using namespace Rcpp;

#include <random>


// done - step 1: run model from R, plot cases per day
// step 2: add risk of transmission based on viral load
// step 3: verify contact tracing working
// step 4: add frequent testing with parameter for frequency
// step 5: add parameters for success rate of various steps


// object storing info for a person. Contains:
// day of infection
// v2: list of infected contacts (implicitly keep track of uninfected contacts based on SAR)  
// uses value of INT_MAX to indicate state isn't reached 
struct Case {   
  int dayInfected = INT_MAX;
  int dayInfectionDetected = INT_MAX;
  int dayInfectious = INT_MAX;
  int dayNotInfectious = INT_MAX; 
  int dayTraced = INT_MAX;

  Case(int infectedDay, int dayInfectorDetected, const NumericVector params): dayInfected(infectedDay){ 
    dayInfectious = dayInfected + rexp(1, params["TimeToInfectiousRate"])[0];
    
    // contact tracing
    if(dayInfectorDetected != INT_MAX && (runif(1)[0] < params["ProbTracedGivenInfectorDetected"])){
      dayTraced = std::max(dayInfectorDetected,dayInfected) + params["ContactTracingDelay"];
      dayInfectionDetected = dayTraced + params["InfectionDetectionDelay"];
    }
    
    
    dayNotInfectious = dayInfectious + 6; // store to know which cases to skip in transmission generation loop
  }
};








// input: set of parameters describing disease and population
// return: time series with daily true and estimated values for different states
// [[Rcpp::export]]
Rcpp::DataFrame branchingModel(int endDay, int popSize, const NumericVector params){
  std::vector<Case> cases; 

  // add 10 initial cases
  for(int i =0;i<1;i++){
    Case initialCase = Case(0, INT_MAX, params);
    cases.push_back(initialCase);
  }
  
  // iterate over number of days. For each day generate new cases
  for(int day = 0;day<endDay;day++){
    size_t size = cases.size();
    for (size_t i = 0; i < size; ++i){
      // check if they are infectious
      if(day >= cases[i].dayInfectious && day < cases[i].dayNotInfectious){
        // todo: split mixing by age
        double transmitRate = params["DailyTransmissionRate"];
        if(day >= cases[i].dayInfectionDetected){
          // reduce transmissions from detected cases
          transmitRate *= params["RelativeTransmissionRisk_Detected"]; //  Todo: make parameter
        }
        
        int numTransmissions = rpois(1, transmitRate)[0];
        
        if(numTransmissions > 0){
          for(int j = 0; j < numTransmissions; j++){
            Case newCase = Case(day, cases[i].dayInfectionDetected, params);
            cases.push_back(newCase);// for each transmission, create new case and add to end of vector
          }
        }
        
      }
      
    }
  }
  
  size_t numCases = cases.size();
  std::vector<int> infectedDays(numCases); 
  
  for(int i =0; i< numCases; i++){
    infectedDays[i] = cases[i].dayInfected;
  }
  
  // add all case data to caseLog dataframe
  Rcpp::DataFrame df = DataFrame::create( _["InfectedDay"] = infectedDays);
  
  
  return(df);
}




