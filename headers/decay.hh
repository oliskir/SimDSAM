
#ifndef decay_h
#define decay_h 1

#include "particle.hh"


class decay {

public:
  particle * parent;
  particle * grandparent;
  particle * list[3];
  decay () {};
  decay ( int a, particle *b[], int d, particle *c1 = NULL, particle *c2 = NULL );
  int multiplicity;
  int generation;
  void simdecay ( int ieve );
  void transform ();
  void print_info ();
  double * ftheta;
  double * fphi;

  void propagateDaughters ( int mediumIndex, double zStop, double zCoat );

};


#endif
