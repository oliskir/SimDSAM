// #define MY_EXTERN_CPP

#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "../headers/globals.hh"



globals::globals( string inputFileName ) { // constructor

  /* Convert string to char */
  char *ch=new char[inputFileName.size()+1];
  ch[inputFileName.size()]=0;
  memcpy(ch,inputFileName.c_str(),inputFileName.size());


  /* Constants */
  pi = 3.1415927;
  speedOfLight = 299792458;
  amu = 931494.061; // atomic mass unit in keV
  me = 510.998928; // electron mass in keV


  /* User input */
  double val;
  string text, line;
  ifstream inputfile;
  inputfile.open (ch, ios::in );
  if ( inputfile.is_open() ) {

    int n=0;
    while ( inputfile.good() ) {
      inputfile >> text;
      if (text=="#") { // skip comments
	getline (inputfile,line);
      }
      else {	  
	val = atof(text.c_str()); 
	n = n+1;
	if (n==1) Nevents=(int)val;
	if (n==2) ebeam=val; 
	if (n==3) zbeam=(int)val; 
	if (n==4) abeam=(int)val; 
	if (n==5) ztarget=(int)val; 
	if (n==6) atarget=(int)val; 
	if (n==7) zejectile=(int)val; 
	if (n==8) aejectile=(int)val; 
	if (n==9) excitationInitial=val; 
	if (n==10) excitationFinal=val; 
	if (n==11) lifetimeInitial=val; 
	if (n==12) zTarget=val; 
	if (n==13) {
	  if (text=="ON") {doHe3Implantation=true;} 
	  else if (text=="OFF") {doHe3Implantation=false;}
	}
	if (n==14) he3Dose=val; 
	if (n==15) swellConstant=val; 
	if (n==16) tuneImplantationDepth=val; 
	if (n==17) thicknessCarbon=val; 
	if (n==18) gaussianWidthCarbon=val; 
	if (n==19) radiusBeamSpot=val; 
	if (n==20) displacementBeamSpot=val; 
	if (n==21) fwhmISACII=val; 
	if (n==22) distanceToCollimator=val; 
	if (n==23) distanceToGeDetector=val; 
	if (n==24) diameterOfCollimator=val; 
	if (n==25) thetaMaxGamma=val; 
	if (n==26) thetaMaxAlpha=val; 
	if (n==27) stoppingPowerSource=text; 
	if (n==28) {
	  if (text=="ON") {alphaAngularDistribution=true;} 
	  else if (text=="OFF") {alphaAngularDistribution=false;}
        }
	if (n==29) {
	  if (text=="ON") {gammaAngularCorrelation=true;} 
	  else if (text=="OFF") {gammaAngularCorrelation=false;}
        }
	if (n==30) NMonteCarlo=(int)val;
	if (n==31) {
	  if (text=="ON") {doEnergyStraggling=true;} 
	  else if (text=="OFF") {doEnergyStraggling=false;}
	}
	if (n==32) {
	  if (text=="ON") {doAngularStraggling=true;} 
	  else if (text=="OFF") {doAngularStraggling=false;}
	}
	if (n==33) tuneAngularStraggling=val;
	if (n==34) {
	  if (text=="ON") {zeroDegreeGammaDetection=true;} 
	  else if (text=="OFF") {zeroDegreeGammaDetection=false;}
	}
	if (n==35) {
	  if (text=="ON") {zeroDegreeAlphaDetection=true;} 
	  else if (text=="OFF") {zeroDegreeAlphaDetection=false;}
	}
	if (n==36) {
	  if (text=="ON") {addSiResolution=true;} 
	  else if (text=="OFF") {addSiResolution=false;}
	}
	if (n==37) tuneGammaDetEff=val;
	if (n==38) outputFileName=text;
	if (n==39) rootTreeFileName=text;
	if (n==40) rootFileName=text;
	if (n==41) geBinWidth=val;
	if (n==42) siBinWidth=val;
	if (n==43) peakEmin=val;
	if (n==44) peakEmax=val;
	if (n==45) expBackground=val;
	if (n==46) hisNameGe=text;
	if (n==47) hisNameSi=text;
	if (n==48) alphaGateWidth=val;
      }
    }
    inputfile.close();

  }
  else cout << "Unable to open INPUT.CARD" << endl; 

  tuneGammaDetEff = distanceToGeDetector/0.07826;

};




void globals::print() {

  string text;

  cout << endl;
  cout << " ------ USER INPUT ------ " << endl;
  cout << "Number of simulated detected events: " << Nevents << endl;
  cout << "Beam energy: " << ebeam/1000 << " MeV" << endl;
  cout << "Beam Z: " << zbeam << endl;
  cout << "Beam A: " << abeam << endl;
  cout << "Target Z: " << ztarget << endl;
  cout << "Target A: " << atarget << endl;
  cout << "Ejectile Z: " << zejectile << endl;
  cout << "Ejectile A: " << aejectile << endl;
  cout << "Excitation energy (initial state): " << excitationInitial << " keV" << endl;
  cout << "Excitation energy (final state): " << excitationFinal << " keV" << endl;
  cout << "Lifetime (initial state): " << lifetimeInitial << " fs" << endl;
  cout << "Au foil thickness: " << zTarget*1e6 << " um" << endl;
  if (doHe3Implantation) {text="ON";} else {text="OFF";}
  cout << "3He implantation profile: " << text << endl;
  cout << "3He dose: " << he3Dose << " /cm2" << endl;
  cout << "Swell constant, A=" << swellConstant << endl;
  cout << "3He implantation depth multiplied by a factor: " << tuneImplantationDepth << endl;
  cout << "Carbon coating thickness: " << thicknessCarbon*1e6 << " um" << endl;
  cout << "Carbon coating width: " << gaussianWidthCarbon*1e6 << " um" << endl;
  cout << "Beam radius: " << radiusBeamSpot*1e3 << " mm" << endl;
  cout << "Lateral displacement of beam: " << displacementBeamSpot*1e3 << " mm" << endl;
  cout << "Beam energy resolution (FWHM): " << fwhmISACII << " %" << endl;
  cout << "Distance to collimator: " << distanceToCollimator*1e2 << " cm" << endl;
  cout << "Distance to Ge detector: " << distanceToGeDetector*1e2 << " cm" << endl;
  cout << "Diameter of collimator: " << diameterOfCollimator*1e3 << " mm" << endl;
  cout << "Max gamma angle: " << thetaMaxGamma << " deg" << endl;
  cout << "Max alpha angle: " << thetaMaxAlpha << " deg" << endl;
  cout << "Source of stopping-power data: " << stoppingPowerSource << endl;
  if (alphaAngularDistribution) {text="ON";} else {text="OFF";}
  cout << "Alpha angular distribution: " << text << endl;
  if (gammaAngularCorrelation) {text="ON";} else {text="OFF";}
  cout << "Gamma angular correlation: " << text << endl;
  cout << "Size of Monte Carlo samples: " << NMonteCarlo << endl;
  if (doEnergyStraggling) {text="ON";} else {text="OFF";}
  cout << "Energy straggling: " << text << endl;
  if (doAngularStraggling) {text="ON";} else {text="OFF";}
  cout << "Angular straggling: " << text << endl;
  cout << "Half scattering angle multiplied by a factor: " << tuneAngularStraggling << endl;
  if (zeroDegreeAlphaDetection) {text="ON";} else {text="OFF";}
  cout << "Zero-degree alpha detection: " << text << endl;
  if (zeroDegreeGammaDetection) {text="ON";} else {text="OFF";}
  cout << "Zero-degree gamma detection: " << text << endl;
  if (addSiResolution) {text="ON";} else {text="OFF";}
  cout << "Si-detector exp. resolution: " << text << endl;
  cout << "Angular strech factor for gamma detection efficiency: " << tuneGammaDetEff << endl;
  cout << "Results written to: " << outputFileName << " and " << rootTreeFileName << endl;
  cout << "------------------------ " << endl;
  cout << endl;

  // cout << endl;
  // cout << " <><><> <><><> WARNING <><><> <><><> " << endl;
  // cout << " message goes here ... " << endl;
  // cout << " <><><> <><><> <><><> <><><> <><><> " << endl;
  // cout << endl;
  // cout << endl;

};



void globals::inputCardForAnalysis() {

  string rootfile = "alphas_ncal2.root";
  if (zejectile==2 && aejectile==4) rootfile="alphas_ncal2.root";
  else if (zejectile==2 && aejectile==3) rootfile="He3_ncal2.root";

  string dummy="../";
  ofstream outfile;
  outfile.open ( "../ROOTanal/input.anal", ios::out );
  outfile << Nevents << endl;
  outfile << outputFileName << endl;
  outfile << rootFileName << endl;
  outfile << geBinWidth << endl; // Ge bin width
  outfile << siBinWidth << endl; // Si bin width
  outfile << peakEmin << endl; // lower limit of peak region (keV)
  outfile << peakEmax << endl; // upper limit of peak region (keV)
  outfile << expBackground << endl; // background level per bin
  outfile << hisNameGe << endl; // name of exp gamma histogram
  outfile << hisNameSi << endl; // name of exp alpha histogram
  outfile << alphaGateWidth << endl; // width of alpha gate (keV)
  outfile << rootfile << endl; // root file
  outfile << endl;
  outfile.close();

};
