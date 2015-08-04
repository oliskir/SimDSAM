
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <TMath.h>
#include "Math/Interpolator.h"

#include "../headers/stoppingDataTables.hh"
#include "../headers/atomicWeights.hh"


const double scorr = 0.0; // relative correction to stopping powers

string convertInt (int);



stoppingDataTables::stoppingDataTables( int Z1 ) { // constructor

  // load atomic weights
  atomicWeights *theAtomicWeights = new atomicWeights();

  // create interpolation objects
  Sinter_He3 = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);
  Sinter_Au = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);
  Sinter_C = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);
  Sinter_Al = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);
  Sinter_Si = new ROOT::Math::Interpolator(ROOT::Math::Interpolation::kCSPLINE);


  // Loop over materials //
  for (int it=0; it<nmat; ++it) {

    // for defition of materials, see also the following places in the code:
    // stoppingDataTables::GetStopping
    // particle::energyStraggling
    // particle::angularStraggling

    // 3He
    if (it==0) {
      Z2[it] = 2;
      m[it] = 3.0160293191; // NIST
      stoich[it] = 1;
      density[it] = 1.0; // g/cm3
    }

    // Au
    else if (it==1) { 
      Z2[it] = 79;
      m[it] = theAtomicWeights->GetAtomicWeight( Z2[it] );
      stoich[it] = 1;
      density[it] = 19.30; // g/cm3
    }

    // carbon
    else if (it==2) { 
      Z2[it] = 6;
      m[it] = theAtomicWeights->GetAtomicWeight( Z2[it] );
      stoich[it] = 1;
      density[it] = 2.267; // g/cm3
    }

    // Al
    else if (it==3) { 
      Z2[it] = 13;
      m[it] = theAtomicWeights->GetAtomicWeight( Z2[it] );
      stoich[it] = 1;
      density[it] = 2.70; // g/cm3
    }

    // Si
    else if (it==4) { 
      Z2[it] = 14;
      m[it] = theAtomicWeights->GetAtomicWeight( Z2[it] );
      stoich[it] = 1;
      density[it] = 2.3290; // g/cm3
    }


    // Read in stopping data
    string cz1 = convertInt( Z1 );
    string sname, sdummy;
    if (Z1<3 || stoppingPowerSource=="PASS") { // for ions with Z>=3 chose between different energy-loss tables:
      sname = "../data/G4EMLOW.6.32/icru73/z";
      sdummy = ".dat";
    }
    else if (Z1>3 && stoppingPowerSource=="GEANT") {
      sname = "../data/GEANT4_extracted/Geant_default_model/z";
      sdummy = "_GEANT4.dat";
    }
    sname = sname+cz1;
    string cz2 = smaterial( Z1, Z2[it] );
    sname = sname+"_"+cz2+sdummy;
    char *cname=new char[sname.size()+1];
    cname[sname.size()]=0;
    memcpy(cname,sname.c_str(),sname.size());
    ifstream stopdata(cname);
    double x,EAmin,EAmax,EAx[999],Sx[999];
    Int_t n=0,j=0,nl=999;
    if ( stopdata.is_open() ) {
      stopdata >> EAmin;
      stopdata >> EAmax;
      stopdata >> nl;
      stopdata >> nl;
      while (stopdata >> x && j < nl) {
	if (Z1<3 || stoppingPowerSource=="PASS") {
	  if (n%2==0) EAx[j]=x;
	  if (n%2==1) {Sx[j]=x*(1+scorr); j=j+1;}
	}
	else if (Z1>=3 && stoppingPowerSource=="GEANT") {
	  if (n%3==0) EAx[j]=x;
	  if (n%3==1) Sx[j]=x*(1+scorr);
	  if (n%3==2) j=j+1;
	}
	n=n+1;
      }
      stopdata.close();

      // Scale stopping power according to atomic weight
      // and convert from keV/(mg/cm2) to keV/nm
      if ( Z2[it]<99 ) {
	Double_t m0 = theAtomicWeights->GetAtomicWeight( Z2[it] );
	for ( Int_t i = 0; i < nl; ++i ) {
	  Sx[i] = (m0/m[it]) * Sx[i]; // scale
	  Sx[i] = Sx[i] * density[it] * 0.1; // convert
	}
      }
      
      
      // Data arrays for interpolation
      double xvec[nl+2],yvec[nl+2];
      xvec[0] = TMath::Log10( 0.1 ); 
      yvec[0] = TMath::Log10( 0.001 );
      xvec[nl+1] = TMath::Log10( 10e6 ); 
      yvec[nl+1] = TMath::Log10( 0.001 );
      for ( int i = 0; i < nl; ++i ) {   
	xvec[i+1] = TMath::Log10( EAx[i]*1000 );
	yvec[i+1] = TMath::Log10( Sx[i] );
      }
      
      // Interpolation
      if (it==0) {
	Sinter_He3->SetData(nl+2, xvec, yvec);
      }
      else if (it==1) {
	Sinter_Au->SetData(nl+2, xvec, yvec);
      }
      else if (it==2) {
	Sinter_C->SetData(nl+2, xvec, yvec);
      }
      else if (it==3) {
	Sinter_Al->SetData(nl+2, xvec, yvec);
      }
      else if (it==4) {
	Sinter_Si->SetData(nl+2, xvec, yvec);
      }
      
    }
    else {
      cout << "Z1: " << Z1 << endl;
      cout << "Z2: " << Z2[it] << endl;
      cout << "No stopping data for chosen ion-target combination" << endl; 
    }
    
  } // end loop over materials
  
};




/* Converts integer to string */
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}



/* character sequence in name of data file identifying the target material */
string stoppingDataTables::smaterial( int Z1, int Z2 ) {
  string s;
  if (Z1<=2) {
    s = convertInt( Z2 );
  }
  else {
    if (Z2<99) { // Stopping in single-atomic material
      s = convertInt( Z2 );
    }
    else { // Stopping in compound
      s = "empty";
      if (Z2==99) s = "G4_A-150_TISSUE";
      if (Z2==103) s = "G4_ADIPOSE_TISSUE_ICRP";
      if (Z2==104) s = "G4_AIR";
      if (Z2==106) s = "G4_ALUMINUM_OXIDE";
      if (Z2==119) s = "G4_BONE_COMPACT_ICRU";
      if (Z2==120) s = "G4_BONE_CORTICAL_ICRP";
      if (Z2==126) s = "G4_C-552";
      if (Z2==130) s = "G4_CALCIUM_FLUORIDE";
      if (Z2==134) s = "G4_CARBON_DIOXIDE";
      if (Z2==141) s = "54"; // substitute with Xenon
      if (Z2==169) s = "G4_Pyrex_Glass";
      if (Z2==179) s = "G4_KAPTON";
      if (Z2==185) s = "G4_LITHIUM_FLUORIDE";
      if (Z2==189) s = "G4_LITHIUM_TETRABORATE";
      if (Z2==197) s = "G4_METHANE";
      if (Z2==202) s = "G4_MUSCLE_STRIATED_ICRU";
      if (Z2==209) s = "G4_NYLON-6-6";
      if (Z2==215) s = "G4_PHOTO_EMULSION";
      if (Z2==216) s = "G4_PLASTIC_SC_VINYLTOLUENE";
      if (Z2==219) s = "G4_POLYCARBONATE";
      if (Z2==221) s = "G4_POLYETHYLENE";
      if (Z2==222) s = "G4_MYLAR";
      if (Z2==223) s = "G4_LUCITE";
      if (Z2==226) s = "G4_POLYSTYRENE";
      if (Z2==227) s = "G4_TEFLON";
      if (Z2==238) s = "G4_PROPANE";
      if (Z2==245) s = "G4_SILICON_DIOXIDE";
      if (Z2==252) s = "G4_SODIUM_IODIDE";
      if (Z2==263) s = "G4_TISSUE-METHANE";
      if (Z2==264) s = "G4_TISSUE-PROPANE";
      if (Z2==276) s = "G4_WATER_WAPOR";
      if (Z2==277) s = "G4_WATER";
    }
  }
  return s;
}



double stoppingDataTables::GetStopping ( double EA, int index ) {
  double ealog, slog, s;  
  ealog = log10(EA);
  if (index==0) {
    slog = Sinter_He3->Eval(ealog);
  }
  else if (index==1) {
    slog = Sinter_Au->Eval(ealog);
  }
  else if (index==2) {
    slog = Sinter_C->Eval(ealog);
  }
  else if (index==3) {
    slog = Sinter_Al->Eval(ealog);
  }
  else if (index==4) {
    slog = Sinter_Si->Eval(ealog);
  }
  s = pow(10,slog); // keV per nm
  return s;
}



double stoppingDataTables::GetStopping ( double EA, int i1, int i2, double n, double A ) {
  double s,s1,s2,rho1,rho2,mass1,mass2,c;
  rho1 = density[i1];
  mass1 = m[i1];
  rho2 = density[i2];
  mass2 = m[i2];
  s1 = this->GetStopping( EA, i1);
  s2 = this->GetStopping( EA, i2);
  double mass1grams = mass1 * 1.661e-24;
  c = n/(rho1/mass1grams);
  //  cout << "c: " << c << endl;
  s = 1/(1+A*c) * ( s1 + c * rho1/rho2 * mass2/mass1 * s2 );
  return s;
}
