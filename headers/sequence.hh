
#ifndef sequence_h
#define sequence_h 1

#include "decay.hh"
#include "particle.hh"


extern double pi;
extern double zTarget;
extern double radiusBeamSpot;  
extern double thicknessCarbon;   
extern double fwhmISACII;
extern double gaussianWidthCarbon;
extern bool doHe3Implantation;    



class sequence {

public:
  decay * list;
  sequence () {};
  sequence ( int a, decay b[], particle beam, particle target, particle compound );
  int decays;
  int generations;
  void simulate_kinematics ( int ieve );
  void initializeVertex( int ieve);
  double * fImplant;
  double * fCarbon;

private:
  particle beamPart;
  particle targetPart;
  particle compoundNucl;
};


#endif
