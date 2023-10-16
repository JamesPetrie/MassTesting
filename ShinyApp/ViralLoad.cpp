#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
#include <random>
#include <queue> 
#include <memory>


enum State {NORMAL, QUARANTINE, ISOLATION};

// object storing infection timing info for a person. 
// uses value of INT_MAX to indicate state isn't reached
struct Case { 
  int hourInfected = INT_MAX;
  int hourIso = INT_MAX; 
  int hourTraceContacts = INT_MAX;
  int hourTraced = INT_MAX; 
  int infectedBy = -1; // index of the person who infected them
  int totalTransmissions = 0;
  std::vector<int> contacts; // vector of indices of people infected
  State state = NORMAL; 

  Case(int infectedHour,  int infectorId, const NumericVector params): hourInfected(infectedHour), infectedBy(infectorId){
    if(params["ProbIso"] > runif(1)[0]){
      hourIso = hourInfected + rexp(1, 1.0/(5.0*24))[0]; // if they do detect symptoms, we record when they started having noticable symptoms
    }
    
    if(params["IsoAndTraceSameTime"]){
      hourTraceContacts = hourIso;
    }else{
      if(params["ProbIso"] > runif(1)[0]){
        hourTraceContacts = hourInfected + rexp(1, 1.0/(5.0*24))[0]; // if they do detect symptoms, we record when they started having noticable symptoms
      }
    }
  } 

  void recordTransmissions(int numTransmissions){
    totalTransmissions += numTransmissions;
  }
  
  void addContact(int contactIndex){
    contacts.push_back(contactIndex);
  }

  void getTraced(int hour){
    hourTraced = hour;
    if(state == NORMAL){ 
      state = QUARANTINE;
    }
  } 
  
  int update(int hour, std::vector<Case>& cases, const NumericVector params ){
    if(hour == hourIso){
      state = ISOLATION; 
    }
    if(hour == hourTraceContacts){
      //notify contacts
      for(int contactIndex : contacts){ 
        if (params["ProbTrace"] > runif(1)[0]){ 
          cases[contactIndex].getTraced(hour);
        }
      }
    }
    
    // check for onward transmissions
    if(state == QUARANTINE || state == ISOLATION){
      return(0);
    }else{ 
      double lambda = 1/(24*5.0);
      double transmitRate = params["R0"]*lambda*exp(-lambda*(hour - hourInfected));

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
   for(int hour = 0;hour< simulationEndHour;hour++){

    size_t size = cases.size(); // get size before through the loop because we are adding to the list ( number of cases) as we go. 
    for (size_t i = 0; i < size; ++i){ // go through each cases
      int numTransmissions = cases[i].update(hour, cases, params); //call update date function to get how many each person transmits to
      
      if(numTransmissions > 0){
        cases[i].recordTransmissions(numTransmissions);
        for(int j = 0; j < numTransmissions; j++){ // we create new cases and keeps track of who infected, when infected 
          if(cases.size() < maxSize - 1){
            Case newCase = Case(hour, i, params);
            cases.push_back(newCase);// for each transmission, create new case and add to end of vector
            cases[i].addContact(cases.size()-1);
          }
        }
      }
    }
  }

  // At the end of simulation - record the info for each case
  size_t numCases = cases.size();
  std::vector<int> infectedHours(numCases);
  std::vector<int> tracedHours(numCases);
  std::vector<int> infectorIds(numCases);
  std::vector<int> ids(numCases);
  std::vector<int> transmissionCounts(numCases);

  for(int i =0; i< numCases; i++){
    infectedHours[i] = cases[i].hourInfected;
    tracedHours[i] = cases[i].hourTraced;
    infectorIds[i] = cases[i].infectedBy;
    ids[i] = i;
    transmissionCounts[i] = cases[i].totalTransmissions;
  }

  // add all case data to caseLog dataframe
  Rcpp::DataFrame df = DataFrame::create( _["InfectedHour"] = infectedHours,  _["TracedHour"] = tracedHours,
                                          _["InfectorId"] = infectorIds, _["Id"] = ids, _["NumInfected"] = transmissionCounts);

  return(df);
}
