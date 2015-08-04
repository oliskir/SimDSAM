
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "../headers/atomicWeights.hh"



atomicWeights::atomicWeights() { // constructor

  string sname0 = "../data/atomic_weights_NIST.dat";
  char *cname0=new char[sname0.size()+1];
  cname0[sname0.size()]=0;
  memcpy(cname0,sname0.c_str(),sname0.size());
  ifstream atomicweights(cname0);
  int idummy=0;
  if ( atomicweights.is_open() ) {
    for (int i = 0; i < 118; ++i) {
      atomicweights >> idummy;
      atomicweights >> mnist[i];
    }
    atomicweights.close();
  }
  else {
    cout << "Cannot read file with standard atomic weights" << endl; 
  }

};
