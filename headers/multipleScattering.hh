
#ifndef multipleScattering_h
#define multipleScattering_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;


const int Nbig=1e6;


class multipleScattering {

public:
  multipleScattering ( int Z, int A );
  double GetAngle();

private:
  double theMS[Nbig];

};


#endif
