
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
#include <math.h> //math (e.g. powers and square roots)
using namespace std;


#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"

/* personal header files */
#include "../headers/particle.hh"
#include "../headers/decay.hh"
#include "../headers/sequence.hh"
#include "../headers/massTable.hh"
#include "../headers/atomicWeights.hh"
#include "../headers/stoppingDataTables.hh"
#include "../headers/globals.hh"
#include "../headers/depthProfile.hh"
#include "../headers/rootDataContainer.hh"

#include "atomicWeights.cpp"
#include "decay.cpp"
#include "particle.cpp"
#include "sequence.cpp"
#include "massTable.cpp"
#include "stoppingDataTables.cpp"
#include "globals.cpp"
#include "depthProfile.cpp"
#include "detectionEfficiency.cpp"
#include "multipleScattering.cpp"
#include "rootDataContainer.cpp"


// GRIFFIN efficiency array (see also detectionEfficiency.hh and .cpp)
const int NAMAX=34;
const int NEMAX=12;
double EFFmat[NAMAX][NEMAX][3];

// some histograms
TH1F *h_Ealpha;
TH1F *h_Egamma;
TH1F *h_Qalpha;
TH1F *h_Qgamma;
TH1F *h_info;
TH2F *h2d;

/* object for generating random numbers */
TRandom3 RanNumGen_main(1); // set to 0 for random seed

/* +infinity */
const double infinity = 1./0.;

/* Global variables */

bool doLateralDisplacement = false;

/* 
extern double pi;
extern int Nevents;
extern double ebeam;
extern int zbeam;
extern int abeam;
extern int ztarget;
extern int atarget;
extern int zejectile;
extern int aejectile;
extern double excitationInitial;
extern double excitationFinal;
extern double lifetimeInitial;
extern double thetaMaxGamma;
extern double thetaMaxAlpha;
extern int NMonteCarlo;
extern string outputFileName;
*/


double pi;
double speedOfLight; 
double amu; 
double me; 

int Nevents;
double ebeam;
int zbeam;
int abeam;
int ztarget;
int atarget;
int zejectile;
int aejectile;
double excitationInitial;
double excitationFinal;
double lifetimeInitial;

double zTarget;
bool doHe3Implantation;    
double he3Dose;
double swellConstant;
double tuneImplantationDepth;
double thicknessCarbon;
double gaussianWidthCarbon;
double radiusBeamSpot;
double displacementBeamSpot;
double fwhmISACII;

double distanceToCollimator;
double distanceToGeDetector;
double diameterOfCollimator;
double thetaMaxGamma;
double thetaMaxAlpha;
string stoppingPowerSource;
bool alphaAngularDistribution;
bool gammaAngularCorrelation;
int NMonteCarlo;
bool doEnergyStraggling;
bool doAngularStraggling;
double tuneAngularStraggling;
bool zeroDegreeAlphaDetection;
bool zeroDegreeGammaDetection;
bool addSiResolution;
double tuneGammaDetEff;
string outputFileName;
string rootTreeFileName;
string rootFileName;


// 3He depth distribution
depthProfile he3dist(1);


/* function declarations */
void generate_angles ( int, int, double[], double[], double_t kaR, double_t kbR, double thMin, double thMax );
Double_t adbessel ( double_t *x, double_t *par );
extern string convertInt(int number);
void calckR ( double_t Atar, double_t Apr, double_t Mtar, double_t Mpr, double_t m1, double_t m2, double_t Ex, double_t Ecm, double_t& kaR, double_t& kbR);




int main ( int argc, char* argv[] )
{

  string inputFileName;
  if (atoi(argv[1])>0) inputFileName = "../InputCards/INPUTDEF.CARD";
  else inputFileName = argv[2];  
  cout << endl;
//  cout << "User input read from: " << inputFileName << endl;


  // reads in user input from INPUT.CARD
  globals *globalVariables = new globals( inputFileName );
  
  // overwrite certain entries from the input files
  int runNo=0;
  if (atoi(argv[1])>0) {
	inputFileName = "../InputCards/INPUTDEF.CARD";
	Nevents = atof(argv[1]);
	excitationInitial = atof(argv[2]);
	excitationFinal = atof(argv[3]);
	lifetimeInitial = atof(argv[4]);
	gammaAngularCorrelation = atoi(argv[6]);
	string sRunNo = argv[5];
	runNo = atoi(argv[5]);
	if (runNo<100) sRunNo = "0" + sRunNo;
	outputFileName = "../SimOutput/output" + sRunNo + ".dat";
	rootTreeFileName = "../SimOutput/output" + sRunNo + ".root";
	rootFileName = "RootFiles/i" + sRunNo + ".root";
  }

  // print info
//  globalVariables->print();
  globalVariables->inputCardForAnalysis();


  /* load additional libraries */
  gSystem->Load("libMathMore");
  gSystem->Load("libUnuran");


  // define histograms
  h_Ealpha = new TH1F ( "Ea", "Ea", 600, 0., 60000.);
  h_Egamma = new TH1F ( "Eg", "Eg", 6000, 0., 12000.);
  h_Qalpha = new TH1F ( "Qa", "Qa", 90, 0., 45.);
  h_Qgamma = new TH1F ( "Qg", "Qg", 90, 0., 45.);
  h_info = new TH1F ( "info", "info", 1, 0., 1.);
  h2d = new TH2F ( "EQa", "EQa", 45, 0., 45.,300,0.,30000.);

  
  // read in GRIFFIN efficiency data
    /* User input */
  double val;
  ifstream inputfile;
  inputfile.open ("../data/GRIFFIN_efficiency.dat", ios::in );
  if ( inputfile.is_open() ) {
    for (int ie=0; ie<NEMAX; ie++) {
      for (int ia=0; ia<NAMAX; ia++) {
	inputfile >> EFFmat[ia][ie][1]; // energy
	inputfile >> EFFmat[ia][ie][0]; // angle
	inputfile >> EFFmat[ia][ie][2]; // efficiency
	inputfile >> val;               // relative uncertainty (in percent)
      }
    }
    inputfile.close();
  }
  else cout << "Unable to open INPUT.CARD" << endl; 

  
  
  /* create mass table with AME2012 nuclear masses */
  massTable *theMassTable = new massTable();


  /* Create spherical bessel functions */  // (OBSOLETE!)
  TF1 *jBessel[10];
  string fname, index;
  char *cfname;
  for ( int l=0; l<10; l++ ) { 
    fname = "j_";
    index = convertInt(l);
    fname += index;
    cfname = new char [fname.size()+1];
    strcpy (cfname, fname.c_str());
    jBessel[l] = new TF1(cfname, "ROOT::Math::sph_bessel([0],x)", 0, 10);
    jBessel[l]->SetParameter(0, l); // set value of parameter #0
    //    gROOT->GetListOfFunctions()->Add(jBessel[l]);
  }


  // create beam particle
  double mbeam = theMassTable->GetMass(zbeam,abeam);
  double mtarget = theMassTable->GetMass(ztarget,atarget);
  particle beamParticle ( zbeam, abeam, mbeam, infinity );
  beamParticle.elab = ebeam; // set kinetic energy
  double p[3];
  p[0] = 0;
  p[1] = 0;
  p[2] = sqrt( (ebeam+mbeam)*(ebeam+mbeam) - mbeam*mbeam );
  beamParticle.set_plab ( p );
  
  // create target particle
  particle targetParticle ( ztarget, atarget, mtarget, infinity );

  // create compound nucleus
  particle compoundNucleus = beamParticle.fuseWith( targetParticle );
  


  // BEGIN DEFINITION OF REACTION SEQUENCE

  double thMin=0.0, thMax;
  particle * part[2];

  // Simulate reaction and gamma-decay or just reaction?
  int NOD;
  if (excitationInitial==excitationFinal) NOD=1;
  else NOD=2;


  // Angular correlations
  int alphaAD=0, gammaAC=0;
  if (alphaAngularDistribution) alphaAD=21;
  if (gammaAngularCorrelation) gammaAC=21;


  // 0th generation
  // DECAY
  int z11=zejectile, a11=aejectile;
  particle part11 ( z11, a11, theMassTable->GetMass(z11,a11), infinity );
  double tau;
  if (NOD==2) tau=lifetimeInitial/1e15;
  else tau=infinity;
  int z12=zbeam+ztarget-zejectile, a12=abeam+atarget-aejectile;
  particle part12 ( z12, a12, theMassTable->GetMass(z12,a12), tau );
  part12.excitation = excitationInitial;
  part[0] = &part11;
  part[1] = &part12;
  decay d0 = decay ( 2, part, 0 );
  /* It is important that the particles are passed by reference
   so that the particle itself, rather than just a copy of it, is 
   subject to modifications */

  double * ad0_theta = new double[NMonteCarlo]; // theta angular distribution
  double * ad0_phi = new double[NMonteCarlo]; // phi angular distribution
//  cout << "Generating angular distribution for 0th generation ..." << endl;
  thMax = thetaMaxAlpha; // important that this angle is large enough!
  generate_angles ( alphaAD, 0, ad0_theta, ad0_phi, 0, 0, thMin, thMax ); // isotropic in theta (between thMin and thMax) and phi
//  cout << "Done" << endl;
  d0.ftheta = ad0_theta;
  d0.fphi = ad0_phi;


  // 1st generation
  // DECAY
  particle part21 ( 0, 0, 0, infinity );
  particle part22 ( z12, a12, theMassTable->GetMass(z12,a12), infinity );
  part22.excitation = excitationFinal;
  part22.suppressPropagation = true;
  part[0] = &part21;
  part[1] = &part22;
  decay d1 = decay ( 2, part, 1, &part12 );
  double * ad1_theta = new double[NMonteCarlo]; // theta angular distribution
  double * ad1_phi = new double[NMonteCarlo]; // phi angular distribution
//  cout << "Generating angular distribution for 1th generation ..." << endl;
  thMax = thetaMaxGamma; // again, important that this angle is large enough!
  generate_angles ( gammaAC, 0, ad1_theta, ad1_phi, 0, 0, thMin, thMax ); // isotropic in theta and phi
//  cout << "Done" << endl;
  d1.ftheta = ad1_theta;
  d1.fphi = ad1_phi;


  // Construct sequence
  decay * dec = new decay[2];
  dec[0] = d0;
  dec[1] = d1;
  sequence s0 = sequence ( NOD, dec, beamParticle, targetParticle, compoundNucleus );


  // carbon coating thickness distribution (see sequence::initializeVertex())
  if (gaussianWidthCarbon>0) {
    // define distribution function
//    TF1 *f = new TF1("f","1/(exp((x-[0])/[1])+1)",0.,5*(thicknessCarbon+gaussianWidthCarbon));
    TF1 *f = new TF1("f","exp(-0.5*((x-[0])/[1])**2)",0.,5*(thicknessCarbon+gaussianWidthCarbon));
    f->SetParameter(0,thicknessCarbon);
    f->SetParameter(1,gaussianWidthCarbon);
    // create sample
    double * carbonDist = new double[2*NMonteCarlo];
    for (int n=0; n<2*NMonteCarlo; n++) {
      carbonDist[n] = f->GetRandom ( 0., 3*(thicknessCarbon+gaussianWidthCarbon) );
    }
    s0.fCarbon = carbonDist;
  }


  // 3He implantation distribution (see sequence::initializeVertex())
  he3dist.SetStrech(tuneImplantationDepth);
  double * implantDist = new double[NMonteCarlo];
  he3dist.GetSample ( NMonteCarlo, implantDist );
  s0.fImplant = implantDist;


  // final-state particles
  int Nfs=2;
  particle * fs[Nfs];
  fs[0] = &part11; // alpha
  fs[1] = &part21; // gamma

  /* Open output file */
  ofstream outfile;
  char *ch=new char[outputFileName.size()+1];
  ch[outputFileName.size()]=0;
  memcpy(ch,outputFileName.c_str(),outputFileName.size());
  outfile.open ( ch, ios::out );

  
  // Initialise the ntuple
  TTree dsamTree ( "dsamTree", "Simulated data" );
  // Create root data container 
  rootDataContainer rootCont( &dsamTree );
  

  /* Max angles resulting in detection */
  double maxAlphaAngleDetected = 0;
  double maxGammaAngleDetected = 0;

  /* Simulate kinematics of decay sequence */
  cout << "Run #" << runNo << endl;
  cout << "Lifetime: " << lifetimeInitial << " fs" << endl;
  cout << "Begin simulation of " << Nevents << " events" << endl;
  int ndet=0;
  int i=0,i2=0;
  double f1,f2,f12;
  int if12,if12prev;
  while (ndet<Nevents) { 
    //  while (i<Nevents) { 
    i = i+1;
        
    f1 = (double)ndet;
    f2 = (double)Nevents;
    f12 = f1/f2*100.;
    if12 = (int)f12;
    if ( if12>if12prev ) {
	cout << "\r... " << (int)f12 << " %" << flush;
    }
    if12prev = if12;

    int n = i % NMonteCarlo;
    s0.simulate_kinematics (n);

    // 2d-plot for de-bugging:
    h2d->Fill(fs[0]->get_thetalab(),fs[0]->elab,1);

    // propagate final-state particles to detector (if not stopped in Au foil)
    for ( int k=0; k<Nfs; k++ ) {
      if (fs[k]->elab>0) {
	fs[k]->propagateToDetector();
	fs[k]->detection(); 
      }
      else {
	fs[k]->stoppedInFoil();
      }
    }

    bool alphaDetected = fs[0]->isDetected();
    double e0 = fs[0]->getMeasuredEnergy();
    double q0 = fs[0]->getObservedAngle();
    bool gammaDetected = fs[1]->isDetected();
    double e1 = fs[1]->getMeasuredEnergy();
    double q1 = fs[1]->getObservedAngle();
    if (NOD==1) {
      gammaDetected = true;
      e1 = 0;
      q1 = 0;
    }

    // write alpha-gamma coincidence data to file
    if (alphaDetected) i2++;
    if (alphaDetected && gammaDetected) {
      ndet = ndet+1;
      outfile << " " << ndet << " " << e0 << " " << q0 << " " << e1 << " " << q1 << endl;
      rootCont.fill( e0, e1 ); // fill ROOT Tree
      h_Ealpha->Fill(e0,1);
      h_Egamma->Fill(e1,1);
      h_Qalpha->Fill(q0,1);
      h_Qgamma->Fill(q1,1);
    }

    // max angles resulting in detection of both particles:
    if (alphaDetected && gammaDetected) { // alpha and gamma detected
      if (d0.ftheta[n]*180/pi > maxAlphaAngleDetected) { maxAlphaAngleDetected = d0.ftheta[n]*180/pi; }
      if (d1.ftheta[n]*180/pi > maxGammaAngleDetected) { maxGammaAngleDetected = d1.ftheta[n]*180/pi; }
    }

    // After 20,000 events, reduce max emission angles (if possible) to speed up simulation
    if (i==20000) {
      if (thetaMaxAlpha > maxAlphaAngleDetected) {
	thetaMaxAlpha = maxAlphaAngleDetected + 3.0;
	cout << "\nMax alpha angle reduced to: " <<  thetaMaxAlpha << " deg" << endl;
//	cout << "Generating new angular distribution for alphas ..." << endl;
	thMax = thetaMaxAlpha;
	generate_angles ( alphaAD, 0, ad0_theta, ad0_phi, 0, 0, thMin, thMax );
	d0.ftheta = ad0_theta;
//      cout << "Done" << endl;
      }
      if (thetaMaxGamma > maxGammaAngleDetected) {
	thetaMaxGamma = maxGammaAngleDetected + 5.0;
	cout << "Max gamma angle reduced to: " <<  thetaMaxGamma << " deg" << endl;
//	cout << "Generating new angular distribution for gammas ..." << endl;
	thMax = thetaMaxGamma;
	generate_angles ( gammaAC, 0, ad1_theta, ad1_phi, 0, 0, thMin, thMax );
	d1.ftheta = ad1_theta;
//	cout << "Done" << endl;
      }
    }

  }
//  cout << "NTOT: " << i << " " << i2 << endl;
  
  
  int binmax;
  double x;
  // open text file
//  ofstream fila;
//  fila.open ( "../ExEaConv.dat", ios::app );
  // find most probable alpha energy
  binmax = h_Ealpha->GetMaximumBin();
  x = h_Ealpha->GetXaxis()->GetBinCenter(binmax);
//  cout << "\nEa = " << x << " keV" << endl;
//  fila << x << endl;
  // find most probable gamma energy
  binmax = h_Egamma->GetMaximumBin();
  x = h_Egamma->GetXaxis()->GetBinCenter(binmax);
//  cout << "Eg = " << x << " keV" << endl;
  double Eg0 = excitationInitial-excitationFinal;
//  fila << (x-Eg0)/Eg0 << endl;
  // close file
//  fila.close();


  // Open root file, save the ntuple and close the root file
//  string rootFileName = "../SimOutput/test.root";
  TFile ofile ( rootTreeFileName.c_str() , "RECREATE" );
  dsamTree.Write();
  h_Ealpha->Write();
  h_Egamma->Write();
  h_Qalpha->Write();
  h_Qgamma->Write();
  h2d->Write();
  h_info->Fill(0.5,lifetimeInitial);
  h_info->Write();
  ofile.Close();


  // close files
  outfile.close();
  cout << "\nEnd of simulation" << endl;

  // terminate the program:
  return 0;
}







/* generates angular distributions */
void generate_angles ( int ntheta, int nphi, double x[], double y[], double_t kaR, double_t kbR, double thMin, double thMax ) {
  double_t tr;

  /* THETA */
  if ( ntheta==0 ) { // isotropic theta distribution
    for (int n=0; n<NMonteCarlo; n++) {
      tr = RanNumGen_main.Rndm (0); // random number between 0 and 1
      double costh_min = cos(thMax*pi/180.0);
      double costh_max = cos(thMin*pi/180.0);
      double costh = costh_min + (costh_max-costh_min) * tr;
      x[n] = acos( costh );
    }
  }
  else if ( ntheta>=10 && ntheta<=19 ) { // spherical Bessel with j=0-9
    int l = ntheta-10;
    TF1 *ad = new TF1("ad", *adbessel, 0, 10, 3);
    ad->SetParameter (0,l); // set value of parameter 0 (l)
    ad->SetParameter (1,kaR); // set value of parameter 1 (kaR)
    ad->SetParameter (2,kbR); // set value of parameter 2 (kbR)
    for (int n=0; n<NMonteCarlo; n++) {
      x[n] = ad->GetRandom ( 0, 10 );
    }
  }
  else if ( ntheta==9 ) { // fixed value (7 deg)
    for (int n=0; n<NMonteCarlo; n++) {
      x[n] = 7 * pi / 180;
    }
  }
  else if ( ntheta==21 ) { // a0*P_0(costh)+a2*P_2(costh)+a4*P_4(costh)

    /* 
       NB! Here theta is the angle relative to the beam axis (normal kinematics). 
       However, in my simulation theta means something different, so appropriate
       translation is necessary.
    */

    TF1 *ad = new TF1("ad", "[0]+[1]*1./2.*(3.*x*x-1.)+[2]*1./8.*(35.*x*x*x*x-30.*x*x+3.)", -1.0, 1.0);
    ad->SetParameter (0,1.0); // set value of parameter 0 (a0=1)
    ad->SetParameter (1,-0.66); // set value of parameter 1 (a2)
//    ad->SetParameter (1,-0.36); // set value of parameter 1 (a2)
    ad->SetParameter (2,0.0); // set value of parameter 2 (a4)

//    TF1 *ad = new TF1("ad", "1.0-0.9*TMath::Abs(x)",-1.0, 1.0);

    double costh_min = cos(thMax*pi/180.0);
    double costh_max = cos(thMin*pi/180.0);
    for (int n=0; n<NMonteCarlo; n++) {
      double costh = ad->GetRandom ( costh_min, costh_max );
      x[n] = acos(costh);
    }
  }


  /* PHI */
  if ( nphi==0 ) { // isotropic phi distribution
    for (int n=0; n<NMonteCarlo; n++) {
      tr = RanNumGen_main.Rndm (0);
      double phi = 2 * pi * tr;
      y[n] = phi;
    }
  }
};



/* Angular distribution function based on spherical Bessel function */
Double_t adbessel ( Double_t *x, Double_t *par ) {
  Float_t xx=x[0];

  /* name of function */
  int l=(int)par[0];
  string fname = "j_";
  string index = convertInt(l);
  fname += index;
  char *cfname = new char [fname.size()+1];
  strcpy (cfname, fname.c_str());

  /* argument of Bessel function */
  double_t costh = cos(xx);
  double_t kaR = par[1];
  double_t kbR = par[2];
  Float_t qR = sqrt( pow(kaR,2) + pow(kbR,2) - 2*kaR*kbR*costh );

  /* evaluate Bessel function */
  TF1 *f = (TF1*)gROOT->GetFunction( cfname );
  Double_t y = f->Eval( qR, 0, 0, 0);
  Double_t y2 = pow(y,2);
  return y2;
}



/* Calculates kR */
void calckR ( double_t Atar, double_t Apr, double_t Mtar, double_t Mpr, double_t m1, double_t m2, double_t Ex, double_t Ecm, double_t& kaR, double_t& kbR) {

  /* grazing distance */
  double_t R = 1.1 * ( pow(Atar,0.3333) + pow(Apr,0.3333) );

  /* energies and momenta of projectile and target in c.m. */
  double_t Epr = ( pow(Ecm,2) + pow(Mpr,2) - pow(Mtar,2) )/( 2*Ecm );
  double_t Ppr = sqrt( pow(Epr,2) - pow(Mpr,2) );
  double_t Ptar = -Ppr;

  /* p is in units of keV/c so if we divide by hbar we get
     MeV/hbar*c where hbar*c=197000*keV*fm, hence we have the
     wave vector k in units of 1/fm ... */
    
  double_t k1 = ( Ppr - Ptar )/197000;

  /* masses in outgoing channel */
  double_t Me1 = m1;
  double_t Me2 = m2 + Ex;
  
  /* energies and momenta of ejectiles (e1 and e2) in c.m. */
  double_t Ee1 = ( pow(Ecm,2) + pow(Me1,2) - pow(Me2,2) )/( 2*Ecm );
  double_t Pe1 = sqrt( pow(Ee1,2) - pow(Me1,2) );
  double_t Pe2 = -Pe1;
         
  /* wave vector k in units of 1/fm */
  double_t k2 = ( Pe1 - Pe2 )/197000;

  kaR = k1 * R;
  kbR = k2 * R;

};




