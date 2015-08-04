
#ifndef particle_h
#define particle_h 1

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;

#include "stoppingDataTables.hh"
#include "atomicWeights.hh"
#include "multipleScattering.hh"
#include "detectionEfficiency.hh"
#include "depthProfile.hh"



extern double pi;
extern double speedOfLight;
extern double amu;
extern double zTarget;
extern double thicknessHe3;
extern double radiusBeamSpot;
extern double thicknessCarbon;
extern double distanceToCollimator;
extern double distanceToGeDetector;
extern bool doEnergyStraggling;
extern bool doAngularStraggling;
extern double tuneAngularStraggling;
extern bool addSiResolution;
extern bool doHe3Implantation;  
extern double he3Dose;
extern double swellConstant;

// 3He depth distribution (see sim_dsam.cpp)
extern depthProfile he3dist;

extern bool doLateralDisplacement; // see sim_dsam.cpp


// Class 'particle'

class particle {

public:
  particle () {};
  particle ( int z, int a, double m, double tau=0 );
  int charge, amass; double mass;
  double excitation, kineticenergy, momentum[3]; // in the decay rest frame
  double elab;
  void set_dynamic_data ( double x, double e, double p[3] );
  void set_momentum ( double p[3] );
  void print_info ();
  void print_energy () { cout << kineticenergy << " keV" << endl; };
  void print_elab () { cout << elab << " keV" << endl; };
  void print_momentum ();

  void set_plab ( double p[3] );
  double get_plabx() { return plab[0]; };
  double get_plaby() { return plab[1]; };
  double get_plabz() { return plab[2]; };

  void set_lifetime ( double tau ) { lifetime = tau; } ;
  double get_lifetime () { return lifetime; } ;

  void set_position( double r3[3] ) { position[0]=r3[0]; position[1]=r3[1]; position[2]=r3[2]; }; // sets position vector
  double get_x() { return position[0]; }; // returns x-coordinate (m)
  double get_y() { return position[1]; };
  double get_z() { return position[2]; };

  double get_vlabnorm ();
  double get_thetalab(); // returns theta lab angle (deg)
  double get_philab(); // returns phi lab angle (deg)

  bool suppressPropagation;
  void propagate( int mediumIndex, double zStop ); 
  void energyLoss( double time, double dist, int mediumIndex, double &x );  
  /*
    Tracks particle for a fixed period of time if time>0 (and dist=0) and
    along fixed distance in space if dist>0 (and time=0).
    At the end it updates position, elab and plab.
  */
  void energyStraggling( double energy, double dist, int mediumIndex );  // updates elab and plab
  void angularStraggling( double energy, double dist, int mediumIndex );  // updates elab and plab

  void propagateToDetector();

  particle fuseWith( particle part );
  void fusionKinematics( particle part1, particle part2 );

  void detection(); // checks that particle hits detector and, for gamma only, adds experimental resolution
  double getObservedAngle() { return observedAngle; };
  double getMeasuredEnergy() { return measuredEnergy; };
  bool isDetected() { return detected; };
  void stoppedInFoil() { detected = false; };
  void SiDeadLayers( double &de1, double &de2, double &de3 );

  void initialize();

private:
  double plab[3];
  double lifetime;
  stoppingDataTables* STables;
  double position[3];
  atomicWeights *theAtomicWeights;
  multipleScattering *theMS;
  detectionEfficiency *theDetEff;
  double observedAngle; 
  double measuredEnergy; 
  bool detected;

};


#endif
