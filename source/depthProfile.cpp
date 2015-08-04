
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)

#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include "Math/Interpolator.h"

#include "../headers/depthProfile.hh"



depthProfile::depthProfile( double val ) { // constructor

  // strech factor
  strech = val;

  // create interpolation objects
  Sinter = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);

  // Read in SRIM profile
  string sname = "../data/He3profile/srim.dat";
  char *cname=new char[sname.size()+1];
  cname[sname.size()]=0;
  memcpy(cname,sname.c_str(),sname.size());
  ifstream file(cname);
  double x[999],y[999],dummy,sum=0;
  Int_t n=0,j=0;
  double binwidth = 10e-7; // 10-nm bin width (cm)
  if ( file.is_open() ) {
      while (file >> dummy) {
	if (n%2==0) {
	  x[j] = dummy; // depth (nm)
	}
	if (n%2==1) {
	  y[j] = dummy; 
	  sum = sum + binwidth*y[j]; 
	  j++;
	}
	n=n+1;
      }
      file.close();    
  }
  else {
    cout << "srim.dat not found" << endl; 
  }

  // Data arrays for interpolation
  double xvec[j],yvec[j];
  for ( int i = 0; i < j; ++i ) {   
    xvec[i] = x[i];
    yvec[i] = y[i]/sum;
    //    cout << i << " " << xvec[i] << " " << yvec[i] << " " << endl;
  }
      
  // Interpolation
  Sinter->SetData(j, xvec, yvec);

  // beyond this depth the density is set to zero
  xmin = xvec[0]; //0.0;
  xmax = xvec[j-1]; //300.0 * strech;

};




double depthProfile::GetDensity( double_t x ) { // ions/cm
  double xx = x/strech;
  double s=0;
  if (xx<xmin || xx>xmax) {
    s=0;
  }
  else {
    s = Sinter->Eval(xx);
  }
  return max(0.,s);
}


double depthProfile::distribution( double_t *x, double *par ) {
  double s = this->GetDensity(x[0]);
  return s;
}


double depthProfile::GetSample ( int NMC, double x[] ) {
  double s;
  TF1 *f = new TF1("f", this, &depthProfile::distribution, 0., xmax*strech, 0);
  for (int n=0; n<NMC; n++) {
    x[n] = f->GetRandom ( 0., xmax*strech );
  }
  return s;
}
