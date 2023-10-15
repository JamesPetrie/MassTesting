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





//start of MT functions 
//each person can be normal, quarantine or isolation
enum State {NORMAL, QUARANTINE, ISOLATION};



// object storing info for a person. Contains:
// day of infection
// list of infected contacts (implicitly keep track of uninfected contacts based on SAR)
// uses value of INT_MAX to indicate state isn't reached
struct Case { // for a case we want to to know when they were infected, earliest time they received a positive test, 
  int hourInfected = INT_MAX;
  int hourInfectionDetected = INT_MAX;
  int hourSymptoms = INT_MAX; //hour they get symptoms
  int hourTraced = INT_MAX; // time they get traced
  int hourLastTested = INT_MIN; // time they were last tested 
  int infectedBy = -1; // index of the person who infected them
  int totalTransmissions = 0;
  std::vector<int> contacts; // vector of pointers to contact
 State state = NORMAL; 
  bool adheresToTesting; // see if they have adhered to testing

  //constructor for case - have when they were infected, test period ( how often they test), params from mass testing
  Case(int infectedHour,  int infectorId, const NumericVector params): hourInfected(infectedHour), infectedBy(infectorId){
    if( runif(1)[0] < params["ProbTestSymptoms"]){ //randomly sample whether this is a person that both develops symptoms and tests because of them
      hourSymptoms = hourInfected + rexp(1, 1.0/(5.0*24))[0]; // if they do detect symptoms, we record when they started having noticable symptoms

    }

  } 

  void recordTransmission(){
    totalTransmissions++;
  }
  
  void addContact(int contactIndex){
    contacts.push_back(contactIndex);
    totalTransmissions++;
  } //people this person infected 


  void getTraced(int hour, const NumericVector params ){
    std::cout << hour;
    //std::cout << "traced";
    hourTraced = hour;
  } //time u were contact traced but we need to add a delay because it's not immediate

  // todo: create tests queue
  // create state enum
  int update(int hour, std::vector<Case>& cases, const NumericVector params ){
 
    // processing test and decide what to do with people 
    // if no infection detected yet, check for test results
      if(hour == hourSymptoms){ // if you don't have your test results back yet 
        state = ISOLATION; // you isolate if the test is positive and record when they got infected 
        hourInfectionDetected = hour;
        //notify contacts
        //std::cout << contacts.size() << " ";
        for(int contactIndex : contacts){ //for each of these contacts this person infected, we sample if they got traced based on prob they got traced
          if (params["ProbTracedGivenInfectorDetected"] > runif(1)[0]){ // < params["ProbTracedGivenInfectorDetected"]){
            cases[contactIndex].getTraced(hour, params);
          }
        }
        

    }

    //choosing when to get tested and quarantine 
    // check whether they should get tested
    if(state == NORMAL && hour >=  hourTraced){ // they are infected but not aware, everyone modelled is infected 
        state = QUARANTINE; // they were contact traced but they dont know they were positive yet 
    }

    
    // check for onward transmissions
    if(state == QUARANTINE || state == ISOLATION){
      return(0);
    }else{ // if the time is not yet at the point when they are not infectious, they could still infect other people 
      double lambda = 1/(24*5.0);
      double R0 = 3;
      double transmitRate = R0*lambda*exp(-lambda*(hour - hourInfected));

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

  
  // add 10 initial cases
  for(int i =0;i<10;i++){ // 0 is the hour they were infected, -1 is the index 
    Case initialCase = Case(0,  -1,  params);
    cases.push_back(initialCase); //create one case and add to the list 
  }

  int simulationEndHour = endDay*24;
  // iterate over number of days. For each day generate new cases
  //for(int hour = 0;hour<endDay*24;hour++){ // ingnore for now- checks if there is outbreak and change the test period to outbreak test period 
  for(int hour = 0;hour< simulationEndHour;hour++){

    size_t size = cases.size(); // get size before through the loop because we are adding to the list ( number of cases) as we go. 
    for (size_t i = 0; i < size; ++i){ // go through each cases
      
      //how many that case transmits to at that hour 
      int numTransmissions = cases[i].update(hour, cases, params); //call update date function to get how many each erson  transmits to and more info
      
      if(numTransmissions > 0){
        //std::cout << "Case:" << i << ", num transmissions: " << numTransmissions << ", hour = " <<hour;
        for(int j = 0; j < numTransmissions; j++){ // we create new cases and keeps track of who infected, when infected 
          if(cases.size() < maxSize - 1){
            Case newCase = Case(hour, i, params);
            cases.push_back(newCase);// for each transmission, create new case and add to end of vector
            cases[i].addContact(cases.size()-1);
          }else{
            cases[i].recordTransmission();
          }
  
        }
      }

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

//for each of the cases we are making the final data frame
//we keep track of who each case is infected by but we do keep track of how many people erach case infected NumInfected
  for(int i =0; i< numCases; i++){
    infectedHours[i] = cases[i].hourInfected;
    detectedHours[i] = cases[i].hourInfectionDetected;
    tracedHours[i] = cases[i].hourTraced;
    infectorIds[i] = cases[i].infectedBy;
    ids[i] = i;
    transmissionCounts[i] = cases[i].totalTransmissions;
    endHours[i] = simulationEndHour;
  }

  // add all case data to caseLog dataframe
  Rcpp::DataFrame df = DataFrame::create( _["InfectedHour"] = infectedHours, _["DetectedHour"] = detectedHours, _["TracedHour"] = tracedHours,
                                          _["InfectorId"] = infectorIds, _["Id"] = ids, _["NumInfected"] = transmissionCounts, _["SimulationEndHour"] = endHours);


  return(df);
}
