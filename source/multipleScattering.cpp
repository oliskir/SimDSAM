
#include <math.h> //math (e.g. powers and square roots)
using namespace std;

#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include "TROOT.h"

#include "../headers/multipleScattering.hh"


TRandom3 RanNumGen_ms (4);



/* 
   Generates a big sample of multiple-scattering angles 
   distributed according to ~x*exp(-x^2) which is a Taylor
   expansion of the exact formula, sin(theta)*exp(-x^2) valid
   for angles below 10 deg, with x=theta/(sqrt(2)*sigma)
*/

multipleScattering::multipleScattering( int Z, int A )
{
  //  theMS = new double[Nevents];
  double tr, xmax=6./sqrt(2.);
  TF1 *f;
  
  /*
  if (Z==2 && A==4) {
	f = new TF1("f","x*exp(-x*x)+0.0025*exp(-(x-1.5)*(x-1.5)/1.95)",0.,xmax);
  }
  else if (Z==12 && A==23) {
	f = new TF1("f","0.95*x*exp(-x*x)+0.0007*exp(-(x-2.0)*(x-2.0)/1.73)",0.,xmax);	  
  }
  else {
	f = new TF1("f","x*exp(-x*x)",0.,xmax); // upper limit corresponds to 6-sigma
  }  
  */
  
  f = new TF1("f","x*exp(-x*x)",0.,xmax); // upper limit corresponds to 6-sigma

  for (int n=0; n<Nbig; n++) {
    theMS[n] = f->GetRandom ( 0, xmax );
  }
};


double multipleScattering::GetAngle() {
  double angle;
  double tr = RanNumGen_ms.Rndm (0) * (double)Nbig;
  int itr = (int)tr;
  if (itr>=0 && itr<Nbig) {
    angle = theMS[itr];
  }
  else {
    angle = 0;
  }
  return angle;
};

