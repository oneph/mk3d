//Header file including all the functions that advance the distribution function due to electric 
//fields, and update of the Electric fields.
#include <complex>
#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;

void eadvance(complex<double>*****f,threevector ***B,int nmax,double dt,int gridx,int gridy,int gridz)//Use if monoenergetic due to enhanced stability of the centred trapezium scheme
{
  int n,m,xpos,ypos,zpos;
  double emrat = 1.0;//The ratio of the electic charge over the electron mass
  complex<double> CompNul(0.0,0.0);
  complex<double> i(0.0,1.0);
  for(n=0;n<nmax;n++)//No need to do n=0,m=0 since df/dt = 0
    {
    	for(xpos=0;xpos<gridx;xpos++)
	    	for(ypos=0;ypos<gridy;ypos++)
	  			for(zpos=0;zpos<gridz;zpos++)
	    			{
	    			
	    			
	    			
	    			}
    }
}