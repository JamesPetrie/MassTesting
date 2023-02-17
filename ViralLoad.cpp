#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
#include <random>


// done - step 1: run model from R, plot cases per day
// step 2: add risk of transmission based on viral load
// step 3: verify contact tracing working
// step 4: add frequent testing with parameter for frequency
// step 5: add parameters for success rate of various steps





// computes probability of infection for a typical contact given viral load
// [[Rcpp::export]]
inline double probTransmit(double viralLoad, const NumericVector params) {
  return((viralLoad > 1) * (1 - exp(-params["maxProbTransmitPerExposure"] * pow(viralLoad, 0.51) / (pow(viralLoad, 0.51) + pow(8.9e6, 0.51)))));
}

// computes probability of a positive PCR result given viral load
// [[Rcpp::export]]
inline double probPositive(double viralLoad,const NumericVector params) {
  if (viralLoad > pow(10, params["minLogPCRViralLoad"])) {
    return 1;
  } else {
    return 0;
  }
}

// compute viral load at a given time relative to infection
// [[Rcpp::export]]
inline double computeViralLoad(double time, const NumericVector params) {
  double minLogLoad = params["initialLogLoad"]; //-3.0;
  double logLoad;
  if (time < 0) {
    logLoad = minLogLoad;
  } else if (time <= params["timeToPeak"]) {
    logLoad = (params["logPeakLoad"] - minLogLoad)* (time / params["timeToPeak"]) + minLogLoad;
  } else if (time <= params["timeToPeak"] + std::min((double)params["timeFromPeakTo0"],(double)params["maxTimeAfterPeak"])) {
    logLoad = minLogLoad + (params["logPeakLoad"] - minLogLoad) * (1 - ((time - params["timeToPeak"]) / params["timeFromPeakTo0"]));
  } else {
    logLoad = minLogLoad;
  }
  return pow(10, logLoad);
}



// computes expected transmissions over specified hour range relative to day of infection
// [[Rcpp::export]]
double sumTransmissions(double startTime, double endTime,const NumericVector params){
  int n =   10 + 90*params["precision"];
  double a = startTime;
  double b = endTime;
  double h = (b - a)/n;
  double val = h/2.0*(params["contactsPerHour"]*probTransmit(computeViralLoad(a, params), params));
  
  for(int i=1;i<n;i++){
    val += h*(params["contactsPerHour"]*probTransmit(computeViralLoad(a +h*i, params), params));
    
  }
  val += h/2.0*(params["contactsPerHour"]*probTransmit(computeViralLoad(b, params), params));

  return(val);
}

// average over test timing offsets
// [[Rcpp::export]]
double fracAfterPositive(const NumericVector params){

  
  double testPeriod = params["testPeriod"];
  
  double testDelay = params["testDelay"];
  double stopTime = params["timeToPeak"] + std::min((double)params["timeFromPeakTo0"],(double)params["maxTimeAfterPeak"]) ;
  
  
  if(testPeriod <= 0) return(-1.0);
    
// todo: fix so doesn't use same point twice (0 and testPeriod)
  double transmissionsAfter = 0.0;
  int numOffsets = 5 + 20*params["precision"];
  
  for(int i=0;i<numOffsets;i++){
    double offset = ((double)i)/numOffsets*params["testPeriod"];
    double remainingProb = 1.0;
    double testTime = offset;
    while(testTime <= stopTime ){
      double probPos = probPositive(computeViralLoad(testTime, params), params);
      
      transmissionsAfter += 1.0/numOffsets*remainingProb*probPos*sumTransmissions(testTime + testDelay, stopTime, params);
      
      remainingProb = remainingProb*(1-probPos);
      
      testTime = testTime + testPeriod;
    }
  }
  return(transmissionsAfter/sumTransmissions(0, stopTime, params));
}






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
        // double transmitRate = params["DailyTransmissionRate"];
        double viralLoad = computeViralLoad(24.0*(day - cases[i].dayInfected), params);
        double transmitRate = 24.0 * params["contactsPerHour"] * probTransmit(viralLoad, params);
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
