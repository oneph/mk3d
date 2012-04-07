//Header file including all random functions needed within KALOS
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>

using namespace std;

double gaus(double x,double x0,double sig)
{
  double val = 1.0*exp(-1.0*(x-x0)*(x-x0)/(2.0*sig*sig));
  return val;
}

double plgndr(int l,int m, float x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;

  if((m<0)||(m>l)||fabs(x)>1.0)
    {
	  std::cout << "WRONG ARGUMENT FOR plgndr" << std::endl;
      terminate();
    }
 
  pmm = 1.0;
  if(m>0)
    {
      somx2=sqrt((1.0-x)*(1.0+x));
      fact = 1.0;
      for(i=1;i<=m;i++)
	{
	  pmm *= -fact*somx2;
	  fact += 2.0;
	}
    }
  if(l==m)
    return pmm;
  else
      {
	pmmp1=x*(2*m+1)*pmm;
	if(l==(m+1))
	  return pmmp1;
	else
	  {
	    for(ll=m+2;ll<=l;ll++)
	      {
		pll=(x*(2*ll-1.0)*pmmp1 - (ll+m-1.0)*pmm)/(ll-m);
		pmm=pmmp1;
		pmmp1=pll;
	      }
	    return pll;
	  }
      }
}

void triangularlow(int rank,double *input,double **mat)//Solve a lowertriangular matrix by forward substitution
{
  int i;
  double det = 1.0;
  double *sum;
  sum = new double [rank];
  if(mat[0][0] == 0.0)
    {
      cout << "Top element is zero in triangularlow()" << endl;
      terminate();
    }
  for(i=0;i<rank;i++)
    sum[i] = 0.0;
  
  for(i=1;i<rank;i++)//The determinant is solely the sum of the diagonal elements
    {
      for(int j=0;j<i;j++)
	sum[i] = mat[i][j];
    }
  input[0] = input[0]/mat[0][0];
  for(i=1;i<rank;i++)
    {
      input[i] = (input[i] - sum[i])/mat[i][i];
    } 
}

void tridiag(complex<double> a[],complex<double> b[],complex<double> c[],complex<double> r[],complex<double> u[],int n)//r is input vector, u is output vector, a[] is subdiagonal, b[] is diagonal, c[] is superdiagonal
{
  int j;
  complex<double> bet, *p;
  p = new complex<double> [n];
  bet = b[0];
  u[0] = r[0]/bet;
  p[0]=0.0;
  if(b[0] == 0.0)
    {
      cout << "ERROR IN FUNCTION - tridiag()" << endl;
      cout << "B[0] is equal to zero" << endl;
      terminate();
    }
  
  for(j=1;j<n;j++)
    {
      p[j] = c[j-1]/bet;
      bet = b[j] - a[j]*p[j];
      u[j] = (r[j] - a[j]*u[j-1])/bet;
     
    }

  for(j=(n-2);j>=0;j--)
    {
      u[j] -= p[j+1]*u[j+1];
    }

  delete [] p;

}

void tridiagre(double a[],double b[],double c[],complex<double> r[],complex<double> u[],int n)//r is input vector, u is output vector, a[] is subdiagonal, b[] is diagonal, c[] is superdiagonal -  This version is for non-complex matrices.
{
  int j;
  complex<double> bet, *p;
  p = new complex<double> [n];
  bet = b[0];
  u[0] = r[0]/bet;
  if(b[0] == 0.0)
    {
      cout << "B[0] is equal to zero" << endl;
      terminate();
    }
  
  for(j=1;j<n;j++)
    {
      p[j] = c[j-1]/bet;
      bet = b[j] - a[j]*p[j];
      u[j] = (r[j] - a[j]*u[j-1])/bet;
     
    }

  for(j=(n-2);j>=0;j--)
    {
      u[j] -= p[j+1]*u[j+1];
    }

  delete [] p;

}

double fact(int n)//Calculates the factorial of a number
{
  double sum = 0.0;
  for(int i=1;i<=n;i++)
    {
      sum += log(1.0*i);
     }                                                                          
  return exp(sum);
}

double dubfact(int n)//Calculates the double factorial of a number
{
  double sum = 0.0;
  if((n%2)==1)
    {
        for(int i=1;i<=n;i=i+2)
	  {
	    sum += log(1.0*i);
	  } 
	return exp(sum);
    }if((n%2)==0)
      {
        for(int i=2;i<=n;i=i+2)
	  {
	    sum += log(1.0*i);
	  } 
	return exp(sum);
      }if((n==0)||(n==-1))
	return 1;
}

double LegCoeff(int n,int m)//Calculates the coefficents for a legendre expansion of a number
{
  double coeff;
  double num = pow(2.0,1.0*(n-m)/2)*fact((n-m)/2)*dubfact(n+m+1);
  coeff = (1.0*((2*m)+1)*fact(n))/num;
  // cout << n << " " << m << " " << coeff << endl;
  return coeff;
}

void beam(complex<double>*****f,int gridx,int gridy, int gridz,double dz,int nmax)
{
  //  for(int j=0;j<gridz;j++)
    for(int i=(nmax-1);i>=0;i=i-2)
      {
	f[gridx][gridy][gridz][i][0] = LegCoeff(nmax-1,i)*plgndr(i,0,1.0);
      }
}

double smallest(double a,double b)//Function returns the smallest of the two values
{
  if(a<=b)
    {
	return a;
    }else
    {
	return b;
    }
}
    
void curl3d(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri)//Calculates the curl of the vector A, outputs B;
{
    int i,j,k;
    double xcomp,ycomp,zcomp;
    for(i=0;i<gri.getx();i++)
	for(j=0;j<gri.gety();j++)
	    for(k=0;k<gri.getz();k++)
	    {
		xcomp = ((A[i][j+1][k].getz()-A[i][j][k].getz())/ele.gety()) -  ((A[i][j][k+1].gety()-A[i][j][k].gety())/ele.getz());
		ycomp = ((A[i][j][k+1].getx()-A[i][j][k].getx())/ele.getz()) -  ((A[i+1][j][k].getz()-A[i][j][k].getz())/ele.getx());
		zcomp = ((A[i+1][j][k].gety()-A[i][j][k].gety())/ele.getx()) -  ((A[i][j+1][k].getx()-A[i][j][k].getx())/ele.gety());
		B[i][j][k].setx(xcomp);
		B[i][j][k].sety(ycomp);
		B[i][j][k].setz(zcomp);
	    }
}

void swapele(double *data,int a,int b)
{
    double temp = data[a];
    data[a]=data[b];
    data[b]=temp;
}
    
double div3d(threevector ***B,threevector ele,int i,int j,int k)//Calculates the curl of the vector B;
{
    double xcomp,ycomp,zcomp;
    xcomp = (B[i+1][j][k].getx()-B[i][j][k].getx())/ele.getx();
    ycomp = (B[i][j+1][k].gety()-B[i][j][k].gety())/ele.gety();
    zcomp = (B[i][j][k+1].getz()-B[i][j][k].getz())/ele.getz();

    return (xcomp+ycomp+zcomp);
        
}

void curl2dXY(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri)//Calculates the curl of the vector A, outputs B;
{
    int i,j,k;
    double xcomp,ycomp,zcomp;
    for(i=0;i<gri.getx();i++)
	for(j=0;j<gri.gety();j++)
	    for(k=0;k<gri.getz();k++)
	    {
		xcomp = ((A[i][j+1][k].getz()-A[i][j][k].getz())/ele.gety());
		ycomp = ((A[i+1][j][k].getz()-A[i][j][k].getz())/ele.getx());
		zcomp = ((A[i+1][j][k].gety()-A[i][j][k].gety())/ele.getx()) -  ((A[i][j+1][k].getx()-A[i][j][k].getx())/ele.gety());
		B[i][j][k].setx(xcomp);
		B[i][j][k].sety(ycomp);
		B[i][j][k].setz(zcomp);
	    }
}

void curl2dXZ(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri)//Calculates the curl of the vector A, outputs B;
{
    int i,j,k;
    double xcomp,ycomp,zcomp;
    for(i=0;i<gri.getx();i++)
	for(k=0;k<gri.getz();k++)
	{
	    xcomp = -((A[i][j][k+1].gety()-A[i][j][k].gety())/ele.getz());
	    ycomp = 0.0;
	    zcomp = ((A[i+1][j][k].gety()-A[i][j][k].gety())/ele.getx());
	    B[i][j][k].setx(xcomp);
	    B[i][j][k].sety(ycomp);
	    B[i][j][k].setz(zcomp);
	}
}

//Generate a three-dimensional dynamically allocated structure
threevector ***AllocateThreevector3D(int gridx, int gridy, int gridz)
  {
      threevector ***dynamicArray;

      dynamicArray = new threevector**[gridx];
      for(int i=0;i<gridx;i++)  
      {		
	    dynamicArray[i] = new threevector*[gridy];
	    for(int j=0;j<gridy;j++)
	    {
	      dynamicArray[i][j] = new threevector[gridz]; 
	    }
      }
      return dynamicArray;
  }
  
  void FreeThreevector3D(threevector ***dArray,int gridx,int gridy)
  {
    for(int i=0;i<gridx;i++)
    {
	  for(int j=0;j<gridy;j++)
	  {
	    delete [] dArray[i][j];
	  }
	  delete [] dArray[i];
    }
    delete [] dArray;
  }
