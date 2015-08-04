
#include <iostream> //print statements
#include <string> //words and sentences
#include <fstream> //Stream class to both read and write from/to files
#include <sstream>
using namespace std;
#include <math.h> //math (e.g. powers and square roots)
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../headers/massTable.hh"


/* 
   This reads atomic mass evaluation data from file and
   stores it in the mass matrix
*/

massTable::massTable() { // constructor

  // initialize mass matrix
  for ( int i=0; i<200; i++ ) { 
    for ( int j=0; j<300; j++ ) { 
      massmatrix [i][j] = 0;
    }
  }

  string line;
  ifstream ame2003;
  ame2003.open ("../data/mass.mas12", ios::in ); // new mass evaluation !
  if ( ame2003.is_open() ) {

    int n=0;
    while ( ame2003.good() ) {
      if (n<39) {
    	getline (ame2003,line);
      }
      else {
	getline (ame2003,line);	  

	/* these 4 lines convert string (line) to character array (a) */
	char *a=new char[line.size()+1];
	a[line.size()]=0;
	memcpy(a,line.c_str(),line.size());
	
	/* this selects the section of each line where the charge (Z)
	   is written and converts into an integer */
	char *b=new char[4];
	b[0]=a[11];
	b[1]=a[12];
	b[2]=a[13];
	b[3]=0;
	int *Z=new int;
	*Z = atoi(b);
	
	/* this selects the section of each line where the mass number (A)
	   is written and converts into an integer */
	b[0]=a[16];
	b[1]=a[17];
	b[2]=a[18];
	b[3]=0;
	int *A=new int;
	*A = atoi(b);
	
	/* this selects the section of each line where the mass excess 
	   is written and converts into a double*/
	char *c=new char[12+1];
	c[12]=0;
	for ( int i=0; i<12; i++ ) { 
	  c[i]=a[29+i];
	}
	double *dm=new double;
	*dm = atof(c);

	/* Fill mass matrix */
	massmatrix [*Z][*A] = (*A)*amu + (*dm) - (*Z)*me;

	//	if (n<100) {
	//	  cout << *Z << ", " << *A << ", " << *dm << ", " << massmatrix[*Z][*A] << endl;	
	//       	}

      }
      n++;
    }

    ame2003.close();
  }
  else cout << "Unable to open ame2012 file (masses)" << endl; 

};
