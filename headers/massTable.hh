
#ifndef massTable_h
#define massTable_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;



extern double amu;
extern double me;



class massTable {

public:
  massTable ();
  double GetMass( int Z, int A ) { return massmatrix[Z][A]; };

private:
  double massmatrix [200][300];

};


#endif
