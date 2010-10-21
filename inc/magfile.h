//Header file including all the magnetic field function for 3DKalos
#include <complex>
#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;


void magadvance(complex<double>*****f,threevector ***B,int nmax,double dt,int gridx,int gridy,int gridz)//Use if monoenergetic due to enhanced stability of the centred trapezium scheme
{
  int n,m,xpos,ypos,zpos;
  double omegaz,omegax,omegay;//The larmor frequencies qB/m
  double emrat = 1.0;//The ratio of the electic charge over the electron mass
  complex<double> *diag,*supdiag,*subdiag,*in,*out;
  complex<double> CompNul(0.0,0.0);
  complex<double> i(0.0,1.0);
  complex<double> temp = CompNul;
  for(n=1;n<nmax;n++)//No need to do n=0,m=0 since df/dt = 0
    {
      diag = new complex<double>[2*n+1];
      supdiag = new complex<double>[2*n+1];
      subdiag = new complex<double>[2*n+1];
      in = new complex<double>[2*n+1];
      out = new complex<double>[2*n+1];
      for(m=0;m<=2*n;m++)
	{
	  diag[m]=CompNul;
	  supdiag[m]=CompNul;
	  subdiag[m]=CompNul;
	  in[m]=CompNul;
	  out[m]=CompNul;
	}
      for(xpos=0;xpos<gridx;xpos++)
	for(ypos=0;ypos<gridy;ypos++)
	  for(zpos=0;zpos<gridz;zpos++)
	    {
	      omegax = emrat*B[xpos][ypos][zpos].getx();
	      omegay = emrat*B[xpos][ypos][zpos].gety();
	      omegaz = emrat*B[xpos][ypos][zpos].getz();
	      
	      temp = -0.25*dt*n*(n+1)*(((omegay - omegax*i)*f[xpos][ypos][zpos][n][1]) + (omegay + i*omegax)*conj(f[xpos][ypos][zpos][n][1])); 
	      //cout << n << " " << m << " " << temp << " " << f[xpos][ypos][zpos][n][1] << " " <<  conj(f[xpos][ypos][zpos][n][1]) << endl; 
	      in[n] = f[xpos][ypos][zpos][n][0] + temp;
	      diag[n] = 1.0;
	      subdiag[n] = 0.25*dt*n*(n+1)*(omegay - (omegax*i));
	      supdiag[n] = 0.25*dt*n*(n+1)*(omegay + (omegax*i));
	      for(m=1;m<=n;m++)
		{
		  //The positive m equations
		  temp = -(0.25*dt*(n-m)*(n+m+1)*(omegay - i*omegax)*f[xpos][ypos][zpos][n][m+1]) + (0.25*dt*(omegay + i*omegax)*f[xpos][ypos][zpos][n][m-1]);
		  //cout << n << " " << m << " " << temp << endl; 
		  in[n-m] = (1.0 - 0.5*dt*omegaz*m*i)*f[xpos][ypos][zpos][n][m] + temp;
		  diag[n-m] = (1.0 + 0.5*dt*omegaz*m*i);
		  supdiag[n-m] = -0.25*dt*(omegay+i*omegax);
		  subdiag[n-m] = +0.25*dt*(n-m)*(n+m+1)*(omegay-i*omegax);
		  //The negative m equatios
		  temp = -(0.25*dt*(n-m)*(n+m+1)*(omegay + i*omegax)*conj(f[xpos][ypos][zpos][n][m+1])) + 0.25*dt*(omegay - i*omegax)*conj(f[xpos][ypos][zpos][n][m-1]);
		  // cout << n << " " << m << " " << temp << endl; 
		  in[n+m] = (1.0 + 0.5*dt*omegaz*m*i)*conj(f[xpos][ypos][zpos][n][m]) + temp;
  		  diag[n+m] = (1.0 - 0.5*dt*omegaz*m*i);
		  supdiag[n+m] = 0.25*dt*(n-m)*(n+m+1)*(omegay+i*omegax);
		  subdiag[n+m] = -0.25*dt*(omegay-i*omegax);
		}
	      subdiag[0] = CompNul;
	      supdiag[2*n] = CompNul;
	  
	     /*  cout << "In" << endl; */
/* 	      for(m=0;m<=2*n;m++) */
/* 		{ */
/* 		  cout << m << " " << diag[m] << " " << supdiag[m] << " " << subdiag[m] << " " << in[m] << " " << out[m] << endl;  */
/* 		} */

	      tridiag(subdiag,diag,supdiag,in,out,2*n+1);
	      
	   /*    cout << "Out" << endl; */
/* 	      for(m=0;m<=2*n;m++) */
/* 		{ */
/* 		  cout << m << " " << diag[m] << " " << supdiag[m] << " " << subdiag[m] << " " << in[m] << " " << out[m] << endl;  */
/* 		} */
	      
	      for(m=0;m<=n;m++)
		{
		  if(m==0)
		    {
		      f[xpos][ypos][zpos][n][m] = real(out[n-m]);
		    }else
		    {
		      f[xpos][ypos][zpos][n][m] = out[n-m];
		    }
		}
	      
  	    }
      delete [] diag;
      delete [] supdiag;
      delete [] subdiag;
      delete [] in;
      delete [] out;
    }
}
