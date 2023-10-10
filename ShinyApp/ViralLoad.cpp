#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
#include <random>
#include <queue> 
#include <memory>


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
  // todo: check equation against paper that proposed it

  //test
  double midPoint = params["probTransmitMid"];
  //return((viralLoad > 1) * (1 - exp(-params["maxProbTransmitPerExposure"] * pow(viralLoad, 0.51) / (pow(viralLoad, 0.51) + pow(8.9e6, 0.51)))));
  return((viralLoad > 1) * (1 - exp(-params["maxProbTransmitPerExposure"] * pow(viralLoad, 0.51) / (pow(viralLoad, 0.51) + pow(midPoint, 0.51)))));
}

// computes probability of a positive PCR result given viral load - sensitivity of the test
// [[Rcpp::export]] 
inline double probPositive(double viralLoad,const NumericVector params) {
  double maxSensitivity = 0.995; //highest probability that a test is positive 
  double slope = 6.0; //value of k in the Sensitivity model -> changes also when h is changed 
  double midPoint = params["logLimitOfDetection"];
  
  return((viralLoad > 1) * maxSensitivity / (1+ exp(-slope*(log10(viralLoad) - midPoint))));
  
  //if (viralLoad > pow(10, params["minLogPCRViralLoad"])) {
  //  return 1;
  //} else {
  //  return 0;
  //}
}

// compute viral load at a given time relative to infection - the log viral load 
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

// [[Rcpp::export]] // ignore - didnt use anywhere else
NumericVector fracDiscoveredByHour(const NumericVector params){
  
  int testPeriod = params["testPeriod"];
  
  int testDelay = params["testDelay"];
  int stopTime = params["timeToPeak"] + std::min((double)params["timeFromPeakTo0"],(double)params["maxTimeAfterPeak"]) ;
  
  NumericVector fracDiscovered(std::floor(stopTime));
  
  double numOffsets = testPeriod;
  for(int offset =0;offset<testPeriod;offset++){
    double probStillNegative = 1.0;

    for(int time =0;time + testDelay < stopTime;time++){
      if((time - offset)%testPeriod == 0){
        double probPos = probPositive(computeViralLoad(time, params), params);
        probStillNegative = probStillNegative*(1-probPos);
      }
      
      fracDiscovered[time + testDelay] += (1.0/numOffsets)*(1-probStillNegative);
      
    }
  }
  
  return(fracDiscovered);
}


// computes expected transmissions over specified hour range relative to day of infection - integration trap
// [[Rcpp::export]] // used for mass testing and not the CT -> trap integration on the transmission curve 
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

// average over test timing offsets - MT not CT 
// [[Rcpp::export]]
double fracAfterPositive(const NumericVector params){

  
  double testPeriod = params["testPeriod"];
  
  double testDelay = params["testDelay"];
  double stopTime = params["timeToPeak"] + std::min((double)params["timeFromPeakTo0"],(double)params["maxTimeAfterPeak"]) ;
  
  
  if(testPeriod <= 0) return(-1.0);
    
  double transmissionsAfter = 0.0;
  int numOffsets = 6 + 20*params["precision"];
  
  for(int i=0;i<numOffsets;i++){
    double offset = ((double)i)/(numOffsets+1)*params["testPeriod"];
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







//start of MT functions 
//each person can be normal, quarantine or isolation
enum State {NORMAL, QUARANTINE, ISOLATION};

// contains the result int and bool together - contains the test results and when the patient will learn about the results
struct TestResult{
  int hourReady;
  bool result;
  TestResult(int hourReady, bool result): hourReady(hourReady), result(result){ //creates object

  }
};

// object storing info for a person. Contains:
// day of infection
// list of infected contacts (implicitly keep track of uninfected contacts based on SAR)
// uses value of INT_MAX to indicate state isn't reached
struct Case { // for a case we want to to know when they were infected, earliest time they received a positive test, 
  int hourInfected = INT_MAX;
  int hourInfectionDetected = INT_MAX;
  int hourNotInfectious = INT_MAX; //time someone's viral load reaches 0.  
  int hourSymptoms = INT_MAX; //hour they get symptoms
  int hourTraced = INT_MAX; // time they get traced
  int hourLastTested = INT_MIN; // time they were last tested 
  int infectedBy = -1; // index of the person who infected them
  std::vector<Case*> contacts; // vector of pointers to contact
  std::queue<TestResult> testResults; //queue of test results
 State state = NORMAL; 
  bool adheresToTesting; // see if they have adhered to testing

  //constructor for case - have when they were infected, test period ( how often they test), params from mass testing
  Case(int infectedHour, int testPeriod, int infectorId, const NumericVector params): hourInfected(infectedHour), infectedBy(infectorId){
    hourNotInfectious = hourInfected + params["timeToPeak"] + params["timeFromPeakTo0"]; //time not infectious 
    if( runif(1)[0] < params["ProbTestSymptoms"]){ //randomly sample whether this is a person that both develops symptoms and tests because of them
      hourSymptoms = hourInfected + params["timeToPeak"] + params["timeFromPeakToSymptoms"]; // if they do detect symptoms, we record when they started having noticable symptoms
    }
    if( runif(1)[0] < params["fracTest"]){ // of the fraction of the population who tests, they can choose to adhere to getting test or not. 
      adheresToTesting = TRUE;
    }else{
      adheresToTesting = FALSE;
    }

    hourLastTested = infectedHour - (int)(testPeriod*runif(1)[0]); //todo: fix test timing offset (problems with changing test period causing distribution to not be uniform)
  } // get when u were last tested 

  void addContact(Case* contact){
    contacts.push_back(contact);
  } //people this person infected 

  void getTested(int hour,const NumericVector params){
    if(adheresToTesting == FALSE) return; // people who don't follow testing instructions, return nothing 
    
    double viralLoad = computeViralLoad(hour - hourInfected, params); // this is people who test
    double probPos = probPositive(viralLoad, params); // if the test returns positive 
    bool result;
    if( runif(1)[0] < probPos){ //checking if the test actually returns positive
      result = TRUE;
    }else{
      result = FALSE;
    }
    TestResult testResult = TestResult(hour + params["testDelay"], result); // time when u get test results 
    testResults.push(testResult);
    hourLastTested = hour; // keep track of when u last tested. 
  }

  void getTraced(int hour, const NumericVector params ){
    hourTraced = hour + params["ContactTracingDelay"];
  } //time u were contact traced but we need to add a delay because it's not immediate

  // todo: create tests queue
  // create state enum
  int update(int hour, int testPeriod ,const NumericVector params ){
 
    // processing test and decide what to do with people 
    // if no infection detected yet, check for test results
    if(hourInfectionDetected > hour){ //if you don't know you're infected 
      if(!testResults.empty() && testResults.front().hourReady <= hour){ // if you don't have your test results back yet 
        //std::cout << testResults.front().result;
        TestResult testResult = testResults.front();  //pulls out the test result from the queue and take a look at it
        testResults.pop();
        // if positive test result returned, move into isolation and notify contacts
        if(testResult.result == TRUE){ 
          state = ISOLATION; // you isolate if the test is positive and record when they got infected 
          hourInfectionDetected = hour;
          //notify contacts
          //std::cout << contacts.size() << " ";
          for(auto it = contacts.begin(); it != contacts.end(); ++it){ //for each of these contacts this person infected, we sample if they got traced based on prob they got traced
            Case* contact = *it;
            if (params["ProbTracedGivenInfectorDetected"] > runif(1)[0]){ // < params["ProbTracedGivenInfectorDetected"]){
              contact->getTraced(hour, params);
            }
          }

        }
      }

    }

    //choosing when to get tested and quarantine 
    // check whether they should get tested
    if(state == NORMAL){ // they are infected but not aware, everyone modelled is infected 
      // check if need to test again
      if(hour - hourLastTested > testPeriod ){ //if it has been too long, you test again 
        getTested(hour, params);
      }

      if(hour == hourTraced){ // if they were contact traced, 
        state = QUARANTINE; // they were contact traced but they dont know they were positive yet 
    
        if(hour - hourLastTested > 10 ){ // can set this to something else define later -> to say iof they just got tested or if they tested but hvane tgotten resutls back 
        getTested(hour, params); // you get tested 
        }}

      if(hour == hourSymptoms){ //
        if(hour - hourLastTested > 10 ){ // can set this to something else define later 
          getTested(hour, params); // you get tested 
        }}
      
      // todo: check for symptoms
      // todo: check for quarantine (from tracing)
    }else if(state == QUARANTINE){
      // check if need to test again (at potentially different frequency)
      // todo: make parameter
      if(hour - hourLastTested > 24 ){ // test every day if they are in quarantine 
        getTested(hour, params);
      }
      // todo: check if should move back into NORMAL state

    }

    
    // check for onward transmissions
    if( hour > hourNotInfectious){
      return(0);
    }else{ // if the time is not yet at the point when they are not infectious, they could still infect other people 
      double viralLoad = computeViralLoad(hour - hourInfected, params);
      double transmitRate = params["contactsPerHour"] * probTransmit(viralLoad, params);
      //if(hour >= hourInfectionDetected){
        // modify transmit rate based on state
      if(state == QUARANTINE){ // reduction in transmission based on state, modify based on how these people are behaving 
        //transmitRate *= (1 - params["fracQuar"]); // change name of these params -> if either of these are 1, they prevent all transmission 
        transmitRate *= (1 - params["fracIso"]);
      }else if(state == ISOLATION){
        transmitRate *= (1 - params["fracIso"]);
      }
     // }

      int numTransmissions = rpois(1, transmitRate)[0]; // number of transmissions -> pois takes rate give int 
      return(numTransmissions);
    }
  }
};








// input: set of parameters describing disease and population
// return: time series with daily true and estimated values for different states
// [[Rcpp::export]]
Rcpp::DataFrame branchingModel(int endDay, int maxSize, const NumericVector params){
  std::vector<Case> cases; // vector of cases but empty 
  int testPeriod = params["normalTestPeriod"];

  // add 1 initial cases
  for(int i =0;i<1;i++){ // 0 is the hour they were infected, -1 is the index 
    Case initialCase = Case(0, testPeriod, -1,  params);
    cases.push_back(initialCase); //create one case and add to the list 
  }

  int simulationEndHour = endDay*24;
  // iterate over number of days. For each day generate new cases
  //for(int hour = 0;hour<endDay*24;hour++){ // ingnore for now- checks if there is outbreak and change the test period to outbreak test period 
  for(int hour = 0;hour< simulationEndHour;hour++){
    if(hour % 24 == 0){ //start of every dayif anyone tested positive in the last week 
      size_t size = cases.size(); //check number of cases
      bool outbreak = FALSE;
      for (size_t i = 0; i < size; ++i){
        if(cases[i].hourInfectionDetected - hour < 7*24){
          outbreak = TRUE;
          break;
        }
      }
      if(outbreak) testPeriod = params["outbreakTestPeriod"]; // skip for now 
      else testPeriod = params["normalTestPeriod"];
    }

    size_t size = cases.size(); // get size before through the loop because we are adding to the list ( number of cases) as we go. 
    for (size_t i = 0; i < size; ++i){ // go through each cases
      
      //how many that case transmits to at that hour 
      int numTransmissions = cases[i].update(hour, testPeriod, params); //call update date function to get how many each erson  transmits to and more info
      
      if(numTransmissions > 0){
        //std::cout << "Case:" << i << ", num transmissions: " << numTransmissions << ", hour = " <<hour;
        for(int j = 0; j < numTransmissions; j++){ // we create new cases and keeps track of who infected, when infected 
          Case newCase = Case(hour, testPeriod, i, params);
          cases.push_back(newCase);// for each transmission, create new case and add to end of vector
          cases[i].addContact(&cases.back());
        }
      }

    }
    if(cases.size() > maxSize){
      simulationEndHour = hour;
      break;
    }

  }


  size_t numCases = cases.size();
  std::vector<int> infectedHours(numCases);
  std::vector<int> detectedHours(numCases);
  std::vector<int> tracedHours(numCases);
  std::vector<int> infectorIds(numCases);
  std::vector<int> ids(numCases);
  std::vector<int> transmissionCounts(numCases);
  std::vector<int> endHours(numCases);
  std::vector<int> hourNotInfectious(numCases);

//for each of the cases we are making the final data frame
//we keep track of who each case is infected by but we do keep track of how many people erach case infected NumInfected
  for(int i =0; i< numCases; i++){
    infectedHours[i] = cases[i].hourInfected;
    detectedHours[i] = cases[i].hourInfectionDetected;
    tracedHours[i] = cases[i].hourTraced;
    infectorIds[i] = cases[i].infectedBy;
    ids[i] = i;
    transmissionCounts[i] = cases[i].contacts.size();
    endHours[i] = simulationEndHour;
    hourNotInfectious[i] = cases[i].hourNotInfectious;
  }

  // add all case data to caseLog dataframe
  Rcpp::DataFrame df = DataFrame::create( _["InfectedHour"] = infectedHours, _["DetectedHour"] = detectedHours, _["TracedHour"] = tracedHours,
                                          _["InfectorId"] = infectorIds, _["Id"] = ids, _["NumInfected"] = transmissionCounts, _["SimulationEndHour"] = endHours, _["hourNotInfectious"] = hourNotInfectious);


  return(df);
}
