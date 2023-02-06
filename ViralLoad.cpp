#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;



// computes probability of infection for a typical contact given viral load
inline double probTransmit(double viralLoad, const NumericVector params = R_NilValue) {
  return((viralLoad > 1) * (1 - exp(-params["maxProbTransmitPerExposure"] * pow(viralLoad, 0.51) / (pow(viralLoad, 0.51) + pow(8.9e6, 0.51)))));
}

// computes probability of a positive PCR result given viral load
inline double probPositive(double viralLoad,const NumericVector params = R_NilValue) {
  if (viralLoad > 1e3) {
    return 1;
  } else {
    return 0;
  }
}

// compute viral load at a given time relative to infection
// [[Rcpp::export]]
inline double computeViralLoad(double time, const NumericVector params) {
  double logLoad;
  if (time < 0) {
    logLoad = 0;
  } else if (time <= params["timeToPeak"]) {
    logLoad = params["logPeakLoad"] * (time / params["timeToPeak"]);
  } else if (time <= params["timeToPeak"] + std::min((double)params["timeFromPeakTo0"],(double)params["maxTimeAfterPeak"])) {
    logLoad = params["logPeakLoad"] * (1 - ((time - params["timeToPeak"]) / params["timeFromPeakTo0"]));
  } else {
    logLoad = 0;
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

