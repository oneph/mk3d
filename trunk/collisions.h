#include <complex>
#include <iostream>
#include <fstream>
//#include "functions.h"

using namespace std;

void collscatter(complex<double>*****f,const int nmax,double dt,int gridx,int gridy,int gridz,double nu)
{
  const int nmax1 = nmax;
  int n,m,i,j,k;
  double rate = nu;
  double diag[nmax1];
  double subdiag[nmax1];
  double supdiag[nmax1];
  complex<double> in[nmax1];
  complex<double> out[nmax1];
  for(i=0;i<gridx;i++)
    for(j=0;j<gridy;j++)
      for(k=0;k<gridz;k++)
	for(m=0;m<nmax;m++)
	  {
	    for(n=0;n<nmax;n++)
	      {
		diag[n] = (1.0 + (0.5*dt*rate*n*(n+1.0)));
		subdiag[n] = 0.0;
		supdiag[n] = 0.0;
		in[n] = f[i][j][k][n][m];
	      }
	    tridiagre(subdiag,diag,supdiag,in,out,nmax);
	    for(n=0;n<nmax;n++)
	      {	    
		f[i][j][k][n][m] = out[n];
	      }
	  }
}

void collision(complex<double>*****f,const int nmax,double dt,int gridx,int gridy,int gridz,double nu)
{
  const int nmax1 = nmax;
  int n,m,i,j,k;
  double rate = nu;
  double diag[nmax1];
  double subdiag[nmax1];
  double supdiag[nmax1];
  complex<double> in[nmax1];
  complex<double> out[nmax1];
  for(i=0;i<gridx;i++)
    for(j=0;j<gridy;j++)
      for(k=0;k<gridz;k++)
	for(m=0;m<nmax;m++)
	  {
	    for(n=0;n<nmax;n++)
	      {
		diag[n] = (1.0 + (0.25*dt*rate*n*(n+1.0)));
		subdiag[n] = 0.0;
		supdiag[n] = 0.0;
		in[n] = f[i][j][k][n][m]*(1.0 - (0.25*dt*rate*n*(n+1.0)));
	      }
	    tridiagre(subdiag,diag,supdiag,in,out,nmax);
	    for(n=0;n<nmax;n++)
	      {	    
		f[i][j][k][n][m] = out[n];
	      }
	  }
}

