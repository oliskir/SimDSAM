
#ifndef atomicWeights_h
#define atomicWeights_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;




class atomicWeights {

public:
  atomicWeights ();
  double GetAtomicWeight( int iZ ) { return mnist[iZ-1]; };

private:
  double mnist[118];

};


#endif
