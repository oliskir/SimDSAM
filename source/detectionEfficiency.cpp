
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include "Math/Interpolator.h"
#include <TMath.h>
#include <TRandom3.h>


#include "../headers/detectionEfficiency.hh"



detectionEfficiency::detectionEfficiency( int Z ) { // constructor

  if (Z==0) { // Ge

    index = 0;

    
    if (zeroDegreeGammaDetection) {
      // if 'zeroDegreeGammaDetection' selected, detected gammas between 0 and 1 degrees
      qmin = 0.0;
      qmax = 1.0;
    }
    else {

      
      qmin = 1.0;
      qmax = 34.0;
    
      // loop over energy bins
      for (int ie=0; ie<NE; ie++) {

	// energy
	energy[ie] = EFFmat[0][ie][1];

	// interpolation object
	Sinter[ie] = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);
      
	// data arrays for interpolation
	double ang[NA], eff[NA];

	// fill arrays
	for (int ia=0; ia<NA; ia++) { // loop over angular bins
	  ang[ia] = EFFmat[ia][ie][0]; // angle
	  eff[ia] = EFFmat[ia][ie][2]; // efficiency
//	  cout << ang[ia] << " " << eff[ia] << " (" << energy[ie] << ")" << endl; cin.get();
	}
      
	// interpolation
	Sinter[ie]->SetData(NA, ang, eff);

	// find maximum
	int in=1e3;
	double q, din=(double)in, y;
	ymax[ie]=0;
	for (int i=0; i<in; i++) {
	  double di = (double)i;
	  q = qmin + (qmax-qmin)/din * (di+0.5);
	  y = Sinter[ie]->Eval(q);
	  if (y>ymax[ie]) {ymax[ie]=y;}
	}
	if (ymax[ie]<=0) {
	  cout << "problem in detectionEfficiency.cpp" << endl; 
	  cin.get();
	}
	
      }
    }
    
  }
  else { // Si  (goes from 1 to 0 very rapidly around 8.07 deg)

    index = 1;
    qmin = 0.0;
    if (zeroDegreeAlphaDetection) {
      qmax = 1.0;
    }
    else {	
      qmax = 180.0/pi * TMath::ATan(diameterOfCollimator/2.0/distanceToCollimator);  // 8.07 deg
    }

  }
      
};



double detectionEfficiency::getEff( double q1, double e ) {

  double q=q1;
  if (index==0) {
    q = TMath::ATan( tuneGammaDetEff * TMath::Tan(q1*pi/180.0) ) * 180.0/pi;
  }

  if (q>qmin && q<qmax) {
    if (index==0) {
      if (zeroDegreeGammaDetection) {
	return 1;
      }
      else {
	int ie_min=0;
	double de_min=10000;
	for (int ie=0; ie<NE; ie++) { // find appropriate energy bin
	  double de = abs(e-energy[ie]);
	  if (de<de_min) {
	    de_min = de;
	    ie_min = ie;
	  }
	}
	return TMath::Min( 1.0, TMath::Max( 0.0, Sinter[ie_min]->Eval(q)/ymax[ie_min] ) ); 	
      }

    }
    else if (index==1) { 
      return 1; 
    }
  }
  else {
    return 0;
  }

};
