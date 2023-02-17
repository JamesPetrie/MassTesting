#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
#include <random>
#include <queue> 


// done - step 1: run model from R, plot cases per day
// done - step 2: add risk of transmission based on viral load
// done - step 3: verify contact tracing working
// done - step 4: make hourly instead of daily
// step 5: add checks at each timestep time since last test, symptoms,  contact tracing, queue of test results
// step 6: add parameters for success rate of various steps
// step 7: change strategy based on time and case count
// step 8: return number of undetected infectious days, tests per person, number quarantined, number isolated
// optional: make it faster?



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






enum State {NORMAL, QUARANTINE, ISOLATION};

struct TestResult{
  int hourReady;
  bool result;
  TestResult(int hourReady, bool result): hourReady(hourReady), result(result){
    
  }
};

// object storing info for a person. Contains:
// day of infection
// list of infected contacts (implicitly keep track of uninfected contacts based on SAR)  
// uses value of INT_MAX to indicate state isn't reached 
struct Case {   
  int hourInfected = INT_MAX;
  int hourInfectionDetected = INT_MAX;
  int hourInfectious = INT_MAX;
  int hourNotInfectious = INT_MAX; 
  int hourSymptoms = INT_MAX;
  int hourTraced = INT_MAX;
  int hourLastTested = INT_MIN;
  std::vector<Case*> contacts;
  std::queue<TestResult> testResults;
  State state = NORMAL;
  
  Case(int infectedHour, int testPeriod, const NumericVector params): hourInfected(infectedHour){ 
    hourInfectious = hourInfected ;
    hourNotInfectious = hourInfected + params["timeToPeak"] + params["timeFromPeakTo0"];
    if( runif(1)[0] < params["ProbDetectSymptoms"]){
      hourSymptoms = hourInfected + params["timeToPeak"] + params["timeFromPeakToSymptoms"];
    }
    hourLastTested = infectedHour - (int)(testPeriod*runif(1)[0]); //todo: fix test timing offset (problems with changing test period causing distribution to not be uniform)
  }
  
  void getTested(int hour,const NumericVector params){
    double viralLoad = computeViralLoad(hour - hourInfected, params);
    double probPos = probPositive(viralLoad, params);
    bool result;
    if( runif(1)[0] < probPos){
      result = TRUE;
    }else{
      result = FALSE;
    }
    TestResult testResult = TestResult(hour + params["testDelay"], result);
    testResults.push(testResult);
    hourLastTested = hour;
  }

  // todo: create tests queue
  // create state enum
  int update(int hour, int testPeriod ,const NumericVector params ){
    
    
 
    // if no infection detected yet, check for test results
    if(hourInfectionDetected > hour){
      if(!testResults.empty() && testResults.front().hourReady <= hour){
        std::cout << testResults.front().result;
        TestResult testResult = testResults.front();
        testResults.pop();
        // if positive test result returned, move into isolation and notify contacts
        if(testResult.result == TRUE){
          state = ISOLATION;
          hourInfectionDetected = hour;
          //notify contacts
          for(auto it = contacts.begin(); it != contacts.end(); ++it){
            Case* contact = *it;
            if (runif(1)[0] < params["ProbTracedGivenInfectorDetected"]){
              contact->hourTraced = hourInfectionDetected + params["ContactTracingDelay"];
            }
          }
          
        }
      }

    }
    
    if(state == NORMAL){
      // check if need to test again
      if(hour - hourLastTested > testPeriod ){
        getTested(hour, params);
      }
      // todo: check for symptoms
      // todo: check for quarantine (from tracing)
    }else if(state == QUARANTINE){
      // check if need to test again (at potentially different frequency)
      // check if should move back into NORMAL state 
    }
  
    // check for onward transmissions
    if( hour > hourNotInfectious){ 
      return(0);
    }else{
      double viralLoad = computeViralLoad(hour - hourInfected, params);
      double transmitRate = params["contactsPerHour"] * probTransmit(viralLoad, params);
      if(hour >= hourInfectionDetected){
        // modify transmit rate based on state
        if(state == QUARANTINE){
          transmitRate *= params["RelativeTransmissionRisk_Detected"];
        }else if(state == ISOLATION){
          transmitRate *= params["RelativeTransmissionRisk_Detected"];
        }
      }
      
      int numTransmissions = rpois(1, transmitRate)[0];
      return(numTransmissions);
    }
  }
};








// input: set of parameters describing disease and population
// return: time series with daily true and estimated values for different states
// [[Rcpp::export]]
Rcpp::DataFrame branchingModel(int endDay, int maxSize, const NumericVector params){
  std::vector<Case> cases; 
  int testPeriod = params["testPeriod"];
  
  // add 1 initial cases
  for(int i =0;i<1;i++){
    Case initialCase = Case(0, testPeriod,  params);
    cases.push_back(initialCase);
  }
  
  // iterate over number of days. For each day generate new cases
  for(int hour = 0;hour<endDay*24;hour++){
    size_t size = cases.size();
    for (size_t i = 0; i < size; ++i){
      
      int numTransmissions = cases[i].update(hour, testPeriod, params);

      if(numTransmissions > 0){
        for(int j = 0; j < numTransmissions; j++){
          Case newCase = Case(hour, params["testPeriod"], params);
          cases.push_back(newCase);// for each transmission, create new case and add to end of vector
          cases[i].contacts.push_back(&newCase);
        }
      }

    }
    if(cases.size() > maxSize){
      break;
    }

  }
  
  
  size_t numCases = cases.size();
  std::vector<int> infectedHours(numCases); 
  std::vector<int> detectedHours(numCases); 
  
  for(int i =0; i< numCases; i++){
    infectedHours[i] = cases[i].hourInfected;
    detectedHours[i] = cases[i].hourInfectionDetected;
  }
  
  // add all case data to caseLog dataframe
  Rcpp::DataFrame df = DataFrame::create( _["InfectedHour"] = infectedHours, _["DetectedHour"] = detectedHours);
  
  
  return(df);
}
