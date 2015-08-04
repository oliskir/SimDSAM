// my first program in C++

#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <TVector3.h>
#include <TLorentzVector.h>

#include "../headers/particle.hh"
#include "../headers/decay.hh"



// Constructor for class 'decay'
/* specifies how many decay products there are
   and specifies the identity of each one */

decay::decay ( int a, particle *b[], int d, particle *c1, particle *c2 ) {
  multiplicity = a;
  for (int n=0; n<a; n++) {
    list[n] = b[n];
  };
  parent = c1;
  grandparent = c2;
  generation = d;
};



void decay::print_info () {
  cout << "\n----- DECAY INFO -----" << endl;
  cout << "parent:" << endl;
  parent->print_info() ;
  int mult = multiplicity;
  cout << "daughters:" << endl;
  for (int n=0; n<mult; n++) {
    list[n]->print_info() ;
  };
  cout << "energies and momenta (in rest frame of parent):" << endl;
  for (int n=0; n<mult; n++) {
    list[n]->print_energy() ;
    list[n]->print_momentum() ;
  };
  cout << "energies (in lab):" << endl;
  for (int n=0; n<mult; n++) {
    list[n]->print_elab() ;
  };
};



void decay::simdecay ( int ieve ) {

  double m0 = parent->mass; // mass of parent nucleus
  double x0 = parent->excitation; // excitation of parent nucleus
  int dm = this->multiplicity; // number of final-state particles
  
  // initialize daughters
  for ( int f=0; f<dm; f++ ) {
    list[f]->initialize();
  }

  // set position vectors to coincide with position of parent
  double r3[3];
  r3[0] = parent->get_x();
  r3[1] = parent->get_y();
  r3[2] = parent->get_z();
  for ( int f=0; f<dm; f++ ) { // loop over final-state particles
    list[f]->set_position(r3);
  }

  double mf=0; // total mass of all final-state particles
  double xf=0; // total excitation of all final-state particles
  for ( int f=0; f<dm; f++ ) { // loop over final-state particles
    mf += list[f]->mass; // mass of particle
    xf += list[f]->excitation; // excitation of particle
  }
  double Q = m0+x0 - (mf+xf); // energy released in decay

  if (dm==2) { // for two-body decay, use E and p conservation
    double etot = mf + xf + Q;
    double m1=list[0]->mass, m2=list[1]->mass; // (rest) masses
    double x1=list[0]->excitation, x2=list[1]->excitation; // excitation energies
    double m1x=m1+x1, m2x=m2+x2; // masses (including excitation energies)
    double e1 = ( pow(etot,2) + pow(m1x,2) - pow(m2x,2) ) / (2*etot);
    double e2 = etot - e1;
    double p = sqrt( pow(e1,2) - pow(m1x,2) );
    double k1=e1-m1x, k2=e2-m2x;
    list[0]->kineticenergy = k1;
    list[1]->kineticenergy = k2;

    double theta=0, phi=0; // angles
    theta = this->ftheta[ieve];
    phi = this->fphi[ieve];

    // define coordinate system (nx,ny,nz)
    TVector3 nx, ny, nz;
    if (grandparent==NULL) {
      nx[0] = 1;
      nx[1] = 0;
      nx[2] = 0;
      ny[0] = 0;
      ny[1] = 1;
      ny[2] = 0;
      nz[0] = 0;
      nz[1] = 0;
      nz[2] = 1;
    }
    else {
      double pg[3]; // momentum of grandparent
      pg[0] = grandparent->get_plabx();
      pg[1] = grandparent->get_plaby();
      pg[2] = grandparent->get_plabz();
      double pgnorm = sqrt( pow(pg[0],2) + pow(pg[1],2) + pow(pg[2],2) );
      double pp[3]; // momentum of parent
      pp[0] = parent->get_plabx();
      pp[1] = parent->get_plaby();
      pp[2] = parent->get_plabz();
      double ppnorm = sqrt( pow(pp[0],2) + pow(pp[1],2) + pow(pp[2],2) );
      if (pgnorm>0 && ppnorm>0) {
	TVector3 n0, n2;
	// define directional vectors n0, n2
	n0[0] = pg[0]/pgnorm;
	n0[1] = pg[1]/pgnorm;
	n0[2] = pg[2]/pgnorm;
	n2[0] = pp[0]/ppnorm;
	n2[1] = pp[1]/ppnorm;
	n2[2] = pp[2]/ppnorm;
	// calculate nx,ny,nz
	nz = n2;
	ny = n2.Cross(n0);
	ny = ny.Unit();
	nx = ny.Cross(nz);
	nx = nx.Unit();
      }
      else {
	nx[0] = 1;
	nx[1] = 0;
	nx[2] = 0;
	ny[0] = 0;
	ny[1] = 1;
	ny[2] = 0;
	nz[0] = 0;
	nz[1] = 0;
	nz[2] = 1;
      }
    }
    
    double p1x = sin(theta) * cos(phi) * p;
    double p1y = sin(theta) * sin(phi) * p;
    double p1z = cos(theta) * p;
    double p2x = -p1x;
    double p2y = -p1y;
    double p2z = -p1z;

    for ( int i=0; i<3; i++ ) { 
      list[0]->momentum[i] = p1x*nx[i] + p1y*ny[i] + p1z*nz[i];
      list[1]->momentum[i] = p2x*nx[i] + p2y*ny[i] + p2z*nz[i];
    }

  }
  else {
    cout << "Procedure for simulating N-body decay with N>2 not yet implemented" << endl;
  }

};



// transforms energies and momenta to lab frame
void decay::transform () {

  double m0 = parent->mass; // mass of parent
  double x0 = parent->excitation; // excitation of parent
  double k0 = parent->elab; // kinetic energy of parent
  double p0[3]; // momentum of parent
  p0[0] = parent->get_plabx();
  p0[1] = parent->get_plaby();
  p0[2] = parent->get_plabz();
  double e0 = m0 + x0 + k0; // total energy of parent
  double g0 = e0 / (m0+x0); // gamma factor of parent
  double v0norm = sqrt( 1 - 1/pow(g0,2) ); // speed of parent
  double p0norm = sqrt( pow(p0[0],2) + pow(p0[1],2) + pow(p0[2],2) );
  TVector3 v0; // velocity of parent
  if (p0norm>0) {
    v0[0] = v0norm * p0[0]/p0norm;
    v0[1] = v0norm * p0[1]/p0norm;
    v0[2] = v0norm * p0[2]/p0norm;
  }
  else {
    v0[0] = 0;
    v0[1] = 0;
    v0[2] = 0;
  }
  
  int dm = multiplicity; // number of final-state particles
  for ( int f=0; f<dm; f++ ) { // loop over final-state particles
    double mf = list[f]->mass; // mass of particle
    double xf = list[f]->excitation; // excitation of particle
    double kf = list[f]->kineticenergy; // kinetic energy of particle
    TLorentzVector p4f; // 4-momentum of particle
    double ef = mf + xf + kf; // total energy of particle
    p4f[0] = list[f]->momentum[0];
    p4f[1] = list[f]->momentum[1];
    p4f[2] = list[f]->momentum[2];
    p4f[3] = ef;
    //    lorentz transformation (content of p4f is updated)
    p4f.Boost(v0); 
    double kf_t = p4f[3] - (mf+xf); // transformed kinetic energy
    double pf_t[3]; // transformed momentum
    pf_t[0] = p4f[0];
    pf_t[1] = p4f[1];
    pf_t[2] = p4f[2];
    list[f]->elab = kf_t;
    list[f]->set_plab ( pf_t );
  }

};




void decay::propagateDaughters ( int index, double zStop, double zCoat ) {

  // loop over daughters
  int dm = multiplicity; 
  for ( int f=0; f<dm; f++ ) {
    bool dont = list[f]->suppressPropagation;
    if (!dont) {
      list[f]->propagate( index, zStop ); // propagate through Au(+3He)
      if (list[f]->charge<=2) { // only for H and He ions
//	double e0=list[f]->elab;
	list[f]->propagate( 2, zCoat ); // propagate through coating on back side of foil
//	cout << list[f]->charge << "||  thickness: " << (zCoat-zStop)*1e6 << " um,   energy loss: " << e0-list[f]->elab << " keV" << endl;
      }
    }
  } 

};
