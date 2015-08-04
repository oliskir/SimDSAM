
#ifndef stoppingDataTables_h
#define stoppingDataTables_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include "Math/Interpolator.h"


/* Global variables */
const int nmat=5;
extern string stoppingPowerSource;


class stoppingDataTables {

public:
  stoppingDataTables ( int Z1 );
  string smaterial (int,int);
  double GetStopping( double EA, int index ); // returns stopping power in keV/nm

  /* 
     Below function returns stopping power in keV/nm for a compound 1+2, 
     where c is the number-ratio of 2 relative to 1, and A is the volume swell
  */
  double GetStopping( double EA, int i1, int i2, double c, double A ); 

private:
  ROOT::Math::Interpolator* Sinter_He3;
  ROOT::Math::Interpolator* Sinter_Au;
  ROOT::Math::Interpolator* Sinter_C;
  ROOT::Math::Interpolator* Sinter_Al;
  ROOT::Math::Interpolator* Sinter_Si;

  int Z2[nmat];
  double m[nmat];
  double stoich[nmat];
  double density[nmat];

};


#endif
