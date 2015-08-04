
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "../headers/rootDataContainer.hh"


rootDataContainer::rootDataContainer ( TTree * t ) { // constructor
  t0 = t;
  t0->Branch ( "EPART" , &EPART , "EPART/F" );
  t0->Branch ( "EGAMMA" , &EGAMMA , "EGAMMA/F" );
};



void rootDataContainer::fill ( double ep, double eg ) {
  this->reset();
  EPART = ep;
  EGAMMA = eg;
  t0->Fill(); 
}



void rootDataContainer::reset() {  
  EPART=0;
  EGAMMA=0;
}