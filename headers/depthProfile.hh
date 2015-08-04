
#ifndef depthProfile_h
#define depthProfile_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include "Math/Interpolator.h"


/* Global variables */


class depthProfile {

public:
  depthProfile ( double val );
  double GetDensity( double x ); // in units of ions/cm
  double GetSample( int N, double x[] );
  double distribution( double *x, double *par ); 
  void SetStrech( double val) {strech=val;};
 
private:
  ROOT::Math::Interpolator* Sinter;
  double xmin, xmax;
  double strech;

};


#endif
