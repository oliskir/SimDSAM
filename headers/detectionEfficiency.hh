
#ifndef detectionEfficiency_h
#define detectionEfficiency_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include "Math/Interpolator.h"


const int NA=34;
const int NE=12;
extern double EFFmat[NA][NE][3];


extern double pi;
extern double distanceToCollimator;
extern double diameterOfCollimator;
extern bool zeroDegreeAlphaDetection;
extern bool zeroDegreeGammaDetection;
extern double tuneGammaDetEff;


class detectionEfficiency {

public:
  detectionEfficiency ( int Z );
  double getEff( double q, double e ); 

private:
  double qmin,qmax;
  int index;
  double energy[NE];
  double ymax[NE];
  ROOT::Math::Interpolator* Sinter[NA];
};


#endif
