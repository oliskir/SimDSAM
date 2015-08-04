// my first program in C++

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include <TRandom3.h>
#include <math.h> //math (e.g. powers and square roots)
#include <TF1.h>

#include "../headers/sequence.hh"
#include "../headers/decay.hh"
#include "../headers/particle.hh"



TRandom3 RanNumGen_s (3); 



// Constructor for class 'sequence'
/* specifies how many generations and how many decays the sequence consists of
   and specifies the identity of each decay */
sequence::sequence ( int a, decay b[], particle beam, particle target, particle compound ) {
  beamPart = beam;
  targetPart = target;
  compoundNucl = compound;
  list = new decay[a];
  decays = a;
  int gmax=0;
  for (int n=0; n<a; n++) {
    list[n] = b[n];
    int g=b[n].generation;
    if (g>gmax) { gmax=g; }
  };
  generations = gmax;
};



void sequence::simulate_kinematics ( int ieve ) {

  int sg = this->generations; // number of generations
  int sd = this->decays; // number of decays

  initializeVertex(ieve); // this sets the coordinates (x,y,z) of the reaction vertex

  for ( int g=0; g<sg+1; g++ ) { // loop over generations
    for ( int i=0; i<sd; i++ ) { // loop over decays

      decay * d = new decay;
      *d = this->list[i]; // select decay from list
      int gen = d->generation; // generation that this decay belongs to

      if ( gen==g ) { // if decay belongs to generation of present interest

	if (gen==0) {
	  d->parent = &compoundNucl;
	  d->grandparent = NULL;
	}
	else if (gen==1) {
	  d->grandparent = &compoundNucl;	  
	}

	d->simdecay (ieve); // determines energies and momenta of final-state particle in decay rest-frame 
	d->transform (); // updates contents of LAB 4-momentum
	int mediumIndex=1; // 0:3He 1:Au
	double thick;
	if (gaussianWidthCarbon==0) thick=thicknessCarbon;
	else thick=fCarbon[2*ieve+1]; // sample coating thickness
	double zCoat = zTarget + thick; // layer of coating on the back side of the target 
	d->propagateDaughters( mediumIndex, zTarget, zCoat ); // propagation of daughters in medium before decay (including energy loss and scattering)

	//	d->print_info();
      }

    }
  }

};



void sequence::initializeVertex(int ieve) {

  int mediumIndex;
  
  // this particular beam particle:
  particle thisPart = beamPart;
  double pos[3];
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0; // carbon coating thickness

  // sample carbon coating thickness
  if (gaussianWidthCarbon==0) { 
    pos[2] = -thicknessCarbon;
  }
  else {
    pos[2] = -fCarbon[2*ieve]; // (m)
  }
  thisPart.set_position(pos);
//  cout << "Coating = " << -pos[2]*1e6 << " um" << endl; 

  // add ISAC-II resolution (FWHM = 0.2%)
  double resISACII = fwhmISACII/100 * thisPart.elab;
  double sigma = resISACII / 2.35;
  double Eisac = RanNumGen_s.Gaus( thisPart.elab, sigma );
  // update energy and momentum
  thisPart.elab = Eisac;
  double m=thisPart.mass, p[3];
  p[0] = 0;
  p[1] = 0;
  p[2] = sqrt( (Eisac+m)*(Eisac+m) - m*m );
  thisPart.set_plab ( p );

  // Propagate beam to z=0 (energy loss in carbon coating)
  mediumIndex=2; // 0:3He 1:Au 2:C 3:Al 4:Si
  thisPart.propagate( mediumIndex, 0.0 ); 

  // Sample depth of reaction in target
  double depth=0;
  if (doHe3Implantation) {
    depth = fImplant[ieve] * 1e-9; // (m)
  }
//  cout << "3He depth = " << depth*1e6 << " um" << endl; 

  // Propagate beam to z=depth 
  mediumIndex=1;
  thisPart.propagate( mediumIndex, depth ); 

  // calculate excitation energy, kinetic energy and momentum of compound nucleus
  compoundNucl.fusionKinematics( thisPart, targetPart );

  // Sample beam cross-sectional profile
  double r = RanNumGen_s.Rndm(0) * radiusBeamSpot; 
  double q = RanNumGen_s.Rndm(0) * 2*pi;

  // Position in space of compound nucleus
  double R3[3];
  R3[0] = r*cos(q) + displacementBeamSpot;
  R3[1] = r*sin(q);
  R3[2] = depth;
  compoundNucl.set_position(R3);

  //  cout << depth*1e6 << " um" << endl;

};
