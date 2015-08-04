
#include <iostream> //print statements
#include <string> //words and sentences
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <TRandom3.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "../headers/particle.hh"
#include "../headers/stoppingDataTables.hh"
#include "../headers/massTable.hh"



TRandom3 RanNumGen_p (2); 


/* function declarations */
void pol( double, double, double, double&, double&, double&);
void wfunc ( double tau, double ws[] );




particle::particle ( int z, int a, double m, double tau ){ // constructor
  charge = z;
  amass = a;
  mass = m;
  lifetime = tau;
  excitation = 0;
  kineticenergy = 0; // in the decay rest frame
  elab = 0;
  for (int i=0; i<3; i++) {
    momentum[i] = 0; // in the decay rest frame
    plab[i] = 0;
    position[i] = 0;
  }
  if (z>0) {
    STables = new stoppingDataTables(charge);
  }
  theAtomicWeights = new atomicWeights();  // load atomic weights
  theMS = new multipleScattering( charge, amass );  // Generate sample of multiple-scattering angles
  theDetEff = new detectionEfficiency(charge);  // detection efficiency 
  observedAngle = -1.0;
  detected = false;
  suppressPropagation = false;
};


void particle::initialize() {
  observedAngle = -1.0;
  detected = false;  
};


void particle::print_momentum () {
  double pnorm = sqrt ( pow(momentum[0],2) + pow(momentum[1],2) + pow(momentum[2],2) );
  cout << pnorm/1000 << " MeV/c" << endl; 
};


void particle::set_momentum ( double p[3] ) {
  momentum[0]=p[0];
  momentum[1]=p[1]; 
  momentum[2]=p[2]; 
};


void particle::set_plab ( double p[3] ) {
  plab[0]=p[0];
  plab[1]=p[1]; 
  plab[2]=p[2]; 
};


void particle::set_dynamic_data (double x, double e, double p[3] ) {
  excitation = x;
  kineticenergy = e;
  momentum[0] = p[0];
  momentum[1] = p[1];
  momentum[2] = p[2];
};


void particle::print_info () {
  cout << "Z=" << charge << ", A=" << amass;
  cout << ", Ex=" << excitation << " keV" << endl;
};


double particle::get_vlabnorm () {
  double pnorm = sqrt( plab[0]*plab[0] + plab[1]*plab[1] + plab[2]*plab[2] );
  double Etot = sqrt( pnorm*pnorm + mass*mass );
  double vlabnorm = pnorm / Etot; // in units of c
  return vlabnorm;
};


double particle::get_thetalab () {
  double th,phi,r;
  double px=plab[0], py=plab[1], pz=plab[2];
  pol ( px,py,pz, th,phi,r );
  return th * 180/pi;
};


double particle::get_philab () {
  double th,phi,r;
  double px=plab[0], py=plab[1], pz=plab[2];
  pol ( px,py,pz, th,phi,r );
  return phi * 180/pi;
};


/* converts from carthesian to polar coordinates */
void pol ( double x, double y, double z, double& theta, double& phi, double& r ) {
  r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );    
  theta = acos(z/r);
  if (x>0 && y>=0) { phi = atan(y/x); }
  else if (x<0) { phi = atan(y/x) + pi; }
  else if (x>0 && y<0) { phi = atan(y/x) + 2*pi; }
  else if (x==0 && y>=0) { phi = pi/2; }
  else if (x==0 && y<0) { phi = 3*pi/2; }
};





/* 
   Propagation of particle in medium before decay (including energy loss and scattering)
   (it is implicitly assumed that the propagation distance is small compared 
   to the dimensions of the medium so that particle does not leave the medium)
*/
void particle::propagate ( int index, double zStop ) { 
  double tau = lifetime; // in units of s
  double x; // distance traveled (ouput parameter of method 'energyLoss')
  double eold = elab; // store lab energy for calculation of straggling

  // only propogate particles with  0.1fs < tau < 10ps  
  if (tau > 0.1e-15 && tau < 10e-12) { 
    double y = RanNumGen_p.Rndm (0); // random number between 0 and 1
    double t = -tau * TMath::Log( y ); // sample lifetime according to exponential decay distribution
    if (charge>0 && t>0) {
      this->energyLoss( t, 0, index, x ); // x is the traveled distance (output value)
      if (elab>25.0 && doEnergyStraggling) {this->energyStraggling( eold, x, index );}
      if (elab>25.0 && doAngularStraggling) {this->angularStraggling( eold, x, index );}      
    }
  }

  // for stable particles propagate fixed distance:
  if (isinf(tau)) { // isinf(x) returns true if x is infinity
    double z = position[2]; // current z coordinate
    double th = this->get_thetalab() * pi/180; // current theta angle
    double costh = TMath::Cos(th);
    double dist = 0; 
    if (costh>0.0) {dist=(zStop-z)/costh;} // distance to travel to reach z=zStop
    if (charge>0 && dist>0) {
      this->energyLoss( 0, dist, index, x ); // x=dist
      if (elab>25.0 && doEnergyStraggling) {this->energyStraggling( eold, x, index );}    
      if (elab>25.0 && doAngularStraggling) {this->angularStraggling( eold, x, index );}      
    }
  }

};




void particle::energyLoss( double time, double dist, int index, double &xOutput ) {  
	
  double tol=5.0; // error tolerance (keV)
  double de0=1.0; // initial step size (keV)
  double ealog, s1=0, s2=0, s=0, s3=0, s13=0, s32=0, dxerr=0, err=0, v=0, gamma=0, etot=0, dx=0, dx132=0;
  double e=elab, x=0;
  double u=0, umax, du, dxcorr;

  double th  = this->get_thetalab() * pi/180;
  double phi = this->get_philab()   * pi/180;
  double r3[3]; // unit vector pointing in direction of particle trajectory
  r3[0] = sin(th)*cos(phi);
  r3[1] = sin(th)*sin(phi);
  r3[2] = cos(th);

  if (time>0 && dist==0) {
    umax = time;
  }
  else if (time==0 && dist>0) {
    umax = dist;
  }
  else { 
    cout << "Problem in particle::energyLoss" << endl;
  }


  // initial stopping-power value 
  double z = position[2];
  double znm = z*1e9;
  int iAu=1, iHe3=0;
  double f = he3dist.GetDensity(znm); 
  double numberDensity = f * he3Dose;
  if (e>25) {
    if (index==1 && doHe3Implantation) {  // 3He+Au mix 
      s2 = STables->GetStopping( e/amass, iAu, iHe3, numberDensity, swellConstant );
    }
    else {
      s2 = STables->GetStopping( e/amass, index );
    }
    de0 = 1.0*s2; // Determine initial step size so that dx=1nm
  }


  // discrete energy-loss process
  double cumerr=0, projcumerr=0;
  int nstep=0;
  if (e>25) {
    while ( u < umax ) {
      nstep = nstep + 1;

      // determine stopping power (keV/nm)
      s1 = s2;
      z = position[2];
      if (index==1 && doHe3Implantation) { // 3He+Au mix
	iAu=1, iHe3=0;
	znm = z*1e9;
	f = he3dist.GetDensity(znm); 
	numberDensity = f * he3Dose;
	s2 = STables->GetStopping( (e-de0)/amass, iAu, iHe3, numberDensity, swellConstant );
	s3 = STables->GetStopping( (e-0.5*de0)/amass, iAu, iHe3, numberDensity, swellConstant ); // for error estimation
      }
      else {
	s2 = STables->GetStopping( (e-de0)/amass, index );
	s3 = STables->GetStopping( (e-0.5*de0)/amass, index ); // for error estimation
      }
      s = 0.5*(s1+s2);

      // new energy
      e = e - de0;

      // increment in x
      dx = de0/s;

      // more precise determination of increment in x
      s13 = 0.5*(s1+s3);
      s32 = 0.5*(s3+s2);
      dx132 = 0.5*de0/s13 + 0.5*de0/s32;

      // error estimate
      dxerr = dx - dx132;
      err = s*dxerr;

      // update x
      x = x + dx132;

      // update particle position 
      for (int i=0; i<3; i++) {
	position[i] = position[i] + dx132*1e-9*r3[i]; // new position (m)
      }

      // temporal or spatial increment
      if (time>0 && dist==0) {
        etot = e-0.5*de0+mass;
	gamma = etot/mass;
	v = sqrt( 1 - 1/(gamma*gamma) ) * speedOfLight * 1e9; // nm/s
	f = gamma/v;
	du = dx132*f; // temporal increment (s)

      }
      else if (time==0 && dist>0) {
	f = 1e-9;
	du = dx132*f; // spatial increment (m)
      }
      u = u+du;

      // correct last integration step
      if (u>umax) {
	double dxcorr = (u-umax)/f;
	e = e + s*dxcorr; // correct energy
	for (int i=0; i<3; i++) {
	  position[i] = position[i] - dxcorr*1e-9*r3[i]; // correct position
	}	
	u = umax;
      }

      // if particle exits Au foil, end loop
      if (index==1 && z>zTarget) {
	u = umax + 1;
      }

      // end loop if energy < 25 keV
      if (e<25) {
	u = umax + 1;
	e = 0;
      }


      // determine new step size based on error estimate
      cumerr = cumerr + err; // cumulative error
      projcumerr = cumerr + (umax-u)/du*err; // projected cumulative error
      projcumerr = abs(projcumerr);
      if (projcumerr<tol) dx=min(200.0,dx*1.2); // max step size ~ 200 nm
      else dx=max(0.5,dx/5.0); // min step size ~ 0.5 nm
      de0 = dx*s;        

/*
      if (amass==4) {
        cout << "z = " << z*1e6 << " um,  E = " << e/1000. << " MeV,  dx = " << dx << " nm,  err = " << projcumerr << " keV,  (" << err << ")" << endl;
      }
*/

    }
  }
  else {
    e=0;
  }

  // distance traveled in m
  xOutput = x*1e-9; 

  // update particle energy and momentum (in lab only)
  elab = e;
  double pnorm_old = sqrt( plab[0]*plab[0] + plab[1]*plab[1] + plab[2]*plab[2] );
  double pnorm_new = sqrt( (elab+mass)*(elab+mass) - mass*mass );
  for (int i=0; i<3; i++) {
    plab[i] = pnorm_new/pnorm_old * plab[i];
  }

  if (cumerr>2*tol) {
    cout << "Estimated error exceeds 2*tolerance:  err = " << cumerr << " keV" << endl;
    cout << "(particle::energyLoss)" << endl;
  }
}




void particle::energyStraggling( double e, double dist, int index ) {

  double Z1 = charge;
  double eaMeV = (e/amass) / 1000; // energy per nucleon in MeV
  double Z2, m2, density;

  // 3He
  if (index==0) {
    Z2 = 2;
    m2 = 3.0160293191; // NIST
    density = 1.0; // g/cm3
  }

  // Au
  else if (index==1) { 
    Z2 = 79;
    m2 = theAtomicWeights->GetAtomicWeight( (int)Z2 );
    density = 19.30; // g/cm3
  }

  // carbon
  else if (index==2) { 
    Z2 = 6;
    m2 = theAtomicWeights->GetAtomicWeight( (int)Z2 );
    density = 2.267; // g/cm3
  }

  else {
    cout << "material requested with index=" << index << " not defined in particle::energyStraggling" << endl;
  }
  
  double x0 = (dist * 1e2) * (density * 1e6); // distance traveled in ug/cm2
  double r = Z2/m2;
  double OmegaB = sqrt(0.157*x0*pow(Z1,2)*r);  // Bohr straggling
  double k = 1.1 + 0.47*TMath::Log10(eaMeV);   // k factor
  double Omega = k * OmegaB;   // Straggling (keV)

  // update energy and momentum (lab)
  double mean=0, std=Omega;
  double deGauss = RanNumGen_p.Gaus ( mean, std ); // sample gaussian;
  while ( elab+deGauss<0. ) {  // Gaussian distribution with cut-off; many wasted calls, could be optimized!         
    deGauss = RanNumGen_p.Gaus ( mean, std ); 
  }
  elab = elab + deGauss;
  double pnorm_old = sqrt( plab[0]*plab[0] + plab[1]*plab[1] + plab[2]*plab[2] );
  double pnorm_new = sqrt( (elab+mass)*(elab+mass) - mass*mass );
  for (int i=0; i<3; i++) {
    plab[i] = pnorm_new/pnorm_old * plab[i];
  }

}





void particle::angularStraggling( double e, double dist, int index ) {

  double amuMeV = amu / 1000;
  double Z1 = charge;
  double A1 = amass;
  double TA = (e/amass) / 1000;  // energy per nucleon in MeV
  double lowtau = 0; 

  /*
    The value of 'lowtau' determines the screening function.
    lowtau=0 :  Thomas-Fermi 
    lowtau=1 :  Lenz-Jensen 
  */
  
  double Z2,A2,density,rx;
  // 3He
  if (index==0) {
    Z2 = 2;
    A2 = 3.0160293191; // NIST
    density = 1.0; // g/cm3
  }
  // Au
  else if (index==1) { 
    Z2 = 79;
    A2 = theAtomicWeights->GetAtomicWeight( (int)Z2 );
    density = 19.30; // g/cm3
  }
  // carbon
  else if (index==2) { 
    Z2 = 6;
    A2 = theAtomicWeights->GetAtomicWeight( (int)Z2 );
    density = 2.267; // g/cm3
  }
  else {
    cout << "material requested with index=" << index << " not defined in particle::angularStraggling" << endl;
  }

  rx = (dist * 1e2) * (density * 1e6); // distance traveled in ug/cm2

  // cout << endl;
  // cout << "INPUT:" << endl;
  // cout << "Z1, A1 = " << Z1 << ", " << A1 << endl; 
  // cout << "Z2, A2 = " << Z2 << ", " << A2 << endl; 
  // cout << "E/A = " << TA << " MeV/u" << endl; 
  // cout << "rho*x = " << rx << " ug/cm2" << endl; 
  // if (lowtau==0) cout << "Screening function: Thomas-Fermi" << endl; 
  // if (lowtau==1) cout << "Screening function: Lenz-Jensen" << endl; 
  // cout << endl;

  // speed in units of c (relativistic formula)
  double v = sqrt(1-pow(1/(1+TA/amuMeV),2));
  double aB = Z1*Z2 / (137*v) ;
  //  cout << "Born parameter (nonrelativistic approximation): " << aB << endl;
  //  cout << endl; 

  // Reduced kinetic energy of beam particle (epsilon_p)
  double Ep = TA*A1;
  double Bpt = 0.961*Z1*Z2/pow(A2,0.3333);
  double epsp = Ep/Bpt;

  // ratio of alpha-tilde to alpha for beam particle (Eq. 12)
  double Zpt = pow(pow(Z1,0.6667)+pow(Z2,0.6667),0.5);
  double aap = 15.63*epsp / (pow(A2,0.3333)*Zpt) ;

  // ratio of alpha-tilde to alpha for compound nucleus (Eq. 14)
  double Zct = pow(pow(Z1+Z2,0.6667)+pow(Z2,0.6667),0.5);
  double aac = 15.63*epsp*A1 / ((A1+A2)*pow(A2,0.3333)*Zct) ;

  // tau for beam particle
  double taup = 41.5*rx / (A2*pow(Zpt,2)) ;

  // tau for compound nucleus
  double tauc = 41.5*rx / (A2*pow(Zct,2)) ;

  /*
  cout << "REDUCED PARAMETERS:" << endl;
  cout << "eps_p = " << epsp << endl; 
  cout << "(a/a)_p = " << aap << endl; 
  cout << "(a/a)_c = " << aac << endl; 
  cout << "tau_p = " << taup << endl; 
  cout << "tau_c = " << tauc << endl; 
  cout << endl;
  */

  // //  if (charge==12 && excitation==2908.) {
  //   cout << "Z=" << charge << ",  tau = " << taup << endl; 
  //   cout << endl;
  // }


  // determine parameters m;Cm and calculate alpha-tilde

  double m,Cm,a12,aa,tau,a12n[4];
  double * ws = new double[4];

  aa = aap;
  tau = taup;
    
  if (lowtau==0) { // Thomas-Fermi
    m = 0.311; 
    Cm = 1.05;
  }
  else if (lowtau==1) { // Lenz-Jensen
    m = 0.191;
    Cm = 3.45;
  }

  a12 = Cm*pow(tau,1/(2*m)); // alpha-tilde
  a12n[0] = a12; // deg
  
  // Sigmund 1<tau<5
  m = 0.5;
  Cm = 0.25;
  a12 = Cm*pow(tau,1/(2*m));
  a12n[1] = a12;

  // Anne global
  m = 0.89;
  Cm = 0.92;
  a12 = Cm*pow(tau,1/(2*m));
  a12n[2] = a12;
  
  // Anne tau>1e3
  m = 0.91;
  Cm = 1.00;
  a12 = Cm*pow(tau,1/(2*m));
  a12n[3] = a12;
  
  a12 = 0;
  wfunc ( tau, ws );
  for ( int n=0; n<4; n++ ) {
    a12 = a12 + ws[n] * a12n[n];
  }
    
  double a12_dim = a12/aa*180/pi/1000; // deg
  //  cout << charge << endl;
  //  cout << "HALF SCATTERING ANGLE = " << a12_dim << " deg" << endl; 
    

  /* Lateral displacement according to Marwick and Sigmund (1975) */
  double Gam0;
  if (lowtau==0) { // Thomas-Fermi
    Gam0 = 2.18; 
  }
  else if (lowtau==1) { // Lenz-Jensen
    Gam0 = 2.34;
  }
  double x = log10(tau);
  double Gam = 1.78 + (Gam0-1.78) / ( exp((x+1)/0.8) + 1 );
  double y12 = tau/Gam * a12; // rho-tilde
  double y12_dim = y12 * 1/tau / aa/1000; // fraction of target thickness
  //  cout << endl; 
  //  cout << "HALF LATERAL DISPLACEMENT = " << y12_dim*rx / (density*1e5) << " mm" << endl; 



  // Update plab (using scattering angle)

  double sigma = 2.*a12_dim / 2.35 * tuneAngularStraggling;
  double alpha = sqrt(2.) * sigma * theMS->GetAngle() * pi/180.; // scattering angle
  double beta = 2*pi*RanNumGen_p.Rndm (0); // azimuthal angle
  TVector3 rnew;
  rnew.SetX( sin(alpha)*cos(beta) );
  rnew.SetY( sin(alpha)*sin(beta) );
  rnew.SetZ( cos(alpha) );

  double th  = this->get_thetalab() * pi/180;
  double phi = this->get_philab()   * pi/180;
  TVector3 r0; // unit vector pointing in direction of particle trajectory
  r0.SetX( sin(th)*cos(phi) );
  r0.SetY( sin(th)*sin(phi) );
  r0.SetZ( cos(th) );
  TVector3 zAxis(0,0,1);
  TVector3 rotAxis = zAxis.Cross(r0);
  double rotAngle = th;
  rnew.Rotate( rotAngle, rotAxis );

  double pnorm = sqrt( plab[0]*plab[0] + plab[1]*plab[1] + plab[2]*plab[2] );
  for (int i=0; i<3; i++) {
    plab[i] = pnorm * rnew(i);
  }

  // Update position (using lateral displacement)
  if (doLateralDisplacement) {
    TVector3 vu = rotAxis.Unit();
    TVector3 wu = r0.Cross(vu);
    wu = wu.Unit();
    double randang = 2*pi*RanNumGen_p.Rndm (0); // random angle
    TVector3 displ = cos(randang)*vu + sin(randang)*wu; // direction of lateral displacement
    double hl0 = y12_dim*rx / (density*1e8); // [m]
    double sig = 2.*hl0 / 2.35;
    double hl = RanNumGen_p.Gaus (0.,sig); // sample from gaussian distribution
    displ = displ*hl;
    position[0] = position[0] + displ.X();
    position[1] = position[1] + displ.Y();
    position[2] = position[2] + displ.Z();
  }

}



/* weight function to connect different tau-domains */
void wfunc ( double tau, double ws[] ) {
  // reset ws
  for ( int n=0; n<4; n++ ) {
    ws[n] = 0;
  }
  double d;
  double x = log10(tau);
  // different tau-domains
  if (tau<=0.1){ // <0.1
    ws[0] = 1;
  }
  if (tau>0.1 && tau<=1) { // 0.1-1
    d = 1.0;
    ws[0] = (log10(1)-x)/d;
    ws[1] = (x-log10(0.1))/d;
  }
  if (tau>1 && tau<=5) { // 1-5
    ws[1] = 1;
  }
  if (tau>5 && tau<=40) { // 5-40
    d = log10(40)-log10(5);
    ws[1] = (log10(40)-x)/d;
    ws[2] = (x-log10(5))/d;
  }
  if (tau>40 && tau<=500) { // 40-500
    ws[2] = 1;
  }
  if (tau>500 && tau<=1e3) { // 500-1e3
    d = log10(1e3)-log10(500);
    ws[2] = (log10(1e3)-x)/d;
    ws[3] = (x-log10(500))/d;
  }
  if (tau>1e3) { // 1e3-1e6
    ws[3] = 1;
  }
}




void particle::propagateToDetector() {

  double zDist; // z distance to detector
  if (charge==0) { zDist = distanceToGeDetector; } // Ge
  if (charge>0) { zDist = distanceToCollimator; }  // Si

  double pnorm = sqrt( plab[0]*plab[0] + plab[1]*plab[1] + plab[2]*plab[2] );
  double dir[3];
  for (int i=0; i<3; i++) {
    if (pnorm>0) {dir[i] = plab[i]/pnorm;}
    else {dir[i] = 0.;}
  }
  double t = 0;
  if (dir[2]>0) {t = (zDist - position[2]) / dir[2];} // "time" it takes particle to reach z=zDist
  for (int i=0; i<3; i++) {
    position[i] = position[i] + t*dir[i]; // new position at "time" t
  }

  double xy = sqrt( position[0]*position[0] + position[1]*position[1] );
  observedAngle = TMath::ATan( xy/zDist ) * 180/pi;

};




particle particle::fuseWith( particle part ) {

  int Z1=charge, Z2=part.charge;
  int A1=amass,  A2=part.amass;
  double m1=mass, m2=part.mass;
  int Z = Z1+Z2;
  int A = A1+A2;
  double m = m1+m2;
  particle compNucl ( Z, A, m );
  return compNucl;

};




void particle::fusionKinematics( particle part1, particle part2 ) {

  double m1 = part1.mass;
  double m2 = part2.mass;
  double e1 = m1 + part1.elab;
  double e2 = m2 + part2.elab;
  TVector3 p1( part1.get_plabx(), part1.get_plaby(), part1.get_plabz() );
  TVector3 p2( part2.get_plabx(), part2.get_plaby(), part2.get_plabz() );
  double m=m1+m2;
  double e=e1+e2;
  TVector3 p=p1+p2;
  TLorentzVector pL( p, e );
  TVector3 vCM = pL.BoostVector(); // spatial component divided by time component
  pL.Boost(-vCM);
  double Eexc = pL(3) - m;
  double Ekin = e - (m+Eexc);

  double mom[3];
  mom[0] = p(0);
  mom[1] = p(1);
  mom[2] = p(2);
  this->set_dynamic_data ( Eexc, Ekin, mom );
  this->set_plab ( mom );
  this->elab = Ekin;

  //  cout << Eexc/1000 << " " << Ekin/1000 << " " << mom[2] << endl;

};



void particle::detection() { 
  
  // determine detection efficiency (0<eff<1)
  double eff;
  if ( observedAngle>=0) {
    eff = theDetEff->getEff( observedAngle, elab );
  }
  else {
    eff = 0;
    cout << "particle does not have observedAngle" << endl;
  }

  // generate random number
  double y = RanNumGen_p.Rndm (0); // random number between 0 and 1
  if (y<eff) {detected = true;}
  else {detected = false;}

  // Discard alphas below 12 MeV because they are stopped in the DE detector
  if (charge>0 && elab<12000.0) {detected = false;}

  measuredEnergy = 0;
  if (detected) {

    double elab0 = elab;

    // Subtract energy loss in dead layers of Si detectors
    double de1=0,de2=0,de3=0; 
    if (charge>0) {
      this->SiDeadLayers(de1,de2,de3);
      elab0 = elab0 - de1 - de2 - de3;
    }
    
    // add experimental resolution
    double sigma=0;
    if (charge==0 && elab0>0.0) { // Ge
      sigma = 0.0387 * pow(elab0,0.5222);
    }
    else { // Si
      if (addSiResolution) {
	double fwhm = sqrt( 28.*28. + 46.*46. );
	sigma = fwhm/2.35;
      }
      else {
	sigma = 0.0;
      }
    }
    double dE = RanNumGen_p.Gaus (0,sigma);
    elab0 = elab0 + dE;
    elab0 = TMath::Max( 0.0, elab0 );
    measuredEnergy = elab0;

  }

}




void particle::SiDeadLayers( double &de1, double &de2, double &de3 ) {

  double dl1 = 1.481e-7; // equivalent to 40 ug/cm2 Al  (front side)
  double dl2 = 2.073e-8; // equivalent to 40 ug/cm2 Au  (back side)
  double dl3 = dl1;
  double SiDE = 100e-6; // 100 um Si DE detector

  double th = this->get_thetalab() * pi/180; // current theta angle
  double costh = TMath::Cos(th);
  int index;
  double dx, dist, edummy, x;
  
  // DE front
  index = 3; // Al
  dist=dl1/costh; 
  edummy = elab;
  this->energyLoss( 0, dist, index, x );
  de1 = edummy-elab;
  
  // DE 
  index = 4; // Si
  dist=SiDE/costh; 
  this->energyLoss( 0, dist, index, x );
  
  // DE back
  index = 1; // Au
  dist=dl2/costh; 
  edummy = elab;
  this->energyLoss( 0, dist, index, x );
  de2 = edummy-elab;
  
  // E front
  index = 3; // Al
  dist=dl3/costh; 
  edummy = elab;
  this->energyLoss( 0, dist, index, x );
  de3 = edummy-elab;

}
