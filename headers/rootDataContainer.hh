
#ifndef rootDataContainer_h
#define rootDataContainer_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;

#include "TTree.h"
#include "TFile.h"



class rootDataContainer {

public:
  rootDataContainer () {};
  rootDataContainer ( TTree * t );
  void fill ( double ep, double eg );
  
private:
  TTree * t0;
  void reset ();
  Float_t EPART,EGAMMA;     
};



#endif
