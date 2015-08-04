
#ifndef globals_h
#define globals_h


/*
#ifdef MY_EXTERN_CPP
    #define MY_CONFIGURATION_EXTERN
#else
    #define MY_CONFIGURATION_EXTERN extern
#endif


MY_CONFIGURATION_EXTERN double pi;
MY_CONFIGURATION_EXTERN double speedOfLight; 
MY_CONFIGURATION_EXTERN double amu; 
MY_CONFIGURATION_EXTERN double me; 

MY_CONFIGURATION_EXTERN int Nevents;
MY_CONFIGURATION_EXTERN double ebeam;
MY_CONFIGURATION_EXTERN int zbeam;
MY_CONFIGURATION_EXTERN int abeam;
MY_CONFIGURATION_EXTERN int ztarget;
MY_CONFIGURATION_EXTERN int atarget;
MY_CONFIGURATION_EXTERN int zejectile;
MY_CONFIGURATION_EXTERN int aejectile;
MY_CONFIGURATION_EXTERN double excitationInitial;
MY_CONFIGURATION_EXTERN double excitationFinal;
MY_CONFIGURATION_EXTERN double lifetimeInitial;

MY_CONFIGURATION_EXTERN double zTarget;
MY_CONFIGURATION_EXTERN bool doHe3Implantation;    
MY_CONFIGURATION_EXTERN double he3Dose;
MY_CONFIGURATION_EXTERN double swellConstant;
MY_CONFIGURATION_EXTERN double tuneImplantationDepth;
MY_CONFIGURATION_EXTERN double thicknessCarbon;
MY_CONFIGURATION_EXTERN double gaussianWidthCarbon;
MY_CONFIGURATION_EXTERN double radiusBeamSpot;
MY_CONFIGURATION_EXTERN double fwhmISACII;

MY_CONFIGURATION_EXTERN double distanceToCollimator;
MY_CONFIGURATION_EXTERN double distanceToGeDetector;
MY_CONFIGURATION_EXTERN double diameterOfCollimator;

MY_CONFIGURATION_EXTERN double thetaMaxGamma;
MY_CONFIGURATION_EXTERN double thetaMaxAlpha;
MY_CONFIGURATION_EXTERN string stoppingPowerSource;
MY_CONFIGURATION_EXTERN double dx0;
MY_CONFIGURATION_EXTERN double NoIntStep;
MY_CONFIGURATION_EXTERN int NMonteCarlo;
MY_CONFIGURATION_EXTERN bool doEnergyStraggling;
MY_CONFIGURATION_EXTERN bool doAngularStraggling;
MY_CONFIGURATION_EXTERN double tuneAngularStraggling;
MY_CONFIGURATION_EXTERN bool zeroDegreeAlphaDetection;
MY_CONFIGURATION_EXTERN bool zeroDegreeGammaDetection;
MY_CONFIGURATION_EXTERN bool addSiResolution;
MY_CONFIGURATION_EXTERN string outputFileName;
*/


extern double pi;
extern double speedOfLight; 
extern double amu; 
extern double me; 

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

extern double zTarget;
extern bool doHe3Implantation;    
extern double he3Dose;
extern double swellConstant;
extern double tuneImplantationDepth;
extern double thicknessCarbon;
extern double gaussianWidthCarbon;
extern double radiusBeamSpot;
extern double displacementBeamSpot;
extern double fwhmISACII;

extern double distanceToCollimator;
extern double distanceToGeDetector;
extern double diameterOfCollimator;

extern double thetaMaxGamma;
extern double thetaMaxAlpha;
extern string stoppingPowerSource;
extern bool alphaAngularDistribution;
extern bool gammaAngularCorrelation;
extern int NMonteCarlo;
extern bool doEnergyStraggling;
extern bool doAngularStraggling;
extern double tuneAngularStraggling;
extern bool zeroDegreeAlphaDetection;
extern bool zeroDegreeGammaDetection;
extern bool addSiResolution;
extern double tuneGammaDetEff;
extern string outputFileName;
extern string rootTreeFileName;
extern string rootFileName;




class globals {

public:
  globals( string inputFileName );
  void print();
  void inputCardForAnalysis();
  
private:
  double geBinWidth;
  double siBinWidth;
  double peakEmin;
  double peakEmax;
  double expBackground;
  string hisNameGe;
  string hisNameSi;
  double alphaGateWidth;

};


#endif


