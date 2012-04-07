#include <complex>
#include <iostream>
#include <fstream>
//#include "functions.h"
#include <fftw.h>

#include <rfftw.h>

using namespace std;

complex<double> integrate3d(complex<double>*****f,int gridx,int gridy,int gridz,double dx,double dy,double dz,int n,int m)
{
  //A simple integration routine in 3d utilising the trapezium rule
  int i,j,k;
  complex<double> s(0.0,0.0);
  
  for(i=0;i<(gridx);i++)
    for(j=0;j<(gridy);j++)
      for(k=0;k<(gridz);k++)
	{
	  s+= f[i][j][k][n][m];
	}
  
  s=s*dx*dy*dz;

  return s;
}

complex<double> integrate1d(complex<double>*****f,int grid,int pos1,int pos2,double dx,int n,int m,char dirn)
{//This function integrates in 1D along a strip designated by the two coordinates pos1 and pos2 in the direction designated by the char dirn
  complex<double> s(0.0,0.0);
  int i;
  switch(dirn)
    {
    case 'x'://Integration in the x direction
      for(i=0;i<(grid);i++)
	{
	  s+= f[i][pos1][pos2][n][m];
	}
      return s;
      break;
    case 'y':
      for(i=0;i<(grid);i++)
	{
	  s+= f[pos1][i][pos2][n][m];
	}
      return s;
      break;
    case 'z':
      for(i=0;i<(grid);i++)
	{
	  s+= f[pos1][pos2][i][n][m];
	}
      return s;
      break;
    default:
      cout <<"Wrong direction specified within integrate1d()!" << endl;
      terminate();
    }
      
}

threevector CentreOfMass(complex<double> *****f,int gridx,int gridy,int gridz,double dx,double dy,double dz)
{
	int i,j,k;
	threevector Com;
	double mass = real(integrate3d(f,gridx,gridy,gridz,dx,dy,dz,0,0));
	complex<double> s = 0.0;
  
	for(i=0;i<(gridx);i++)
		for(j=0;j<(gridy);j++)
			for(k=0;k<(gridz);k++)
			{
				s+= 1.0*i*dx*f[i][j][k][0][0];
			}
	s=s*dx*dy*dz/mass;
	Com.setx(real(s));
	s=0.0;
	for(i=0;i<(gridx);i++)
		for(j=0;j<(gridy);j++)
			for(k=0;k<(gridz);k++)
			{
				s+= 1.0*j*dy*f[i][j][k][0][0];
			}
	s=s*dx*dy*dz/mass;
	Com.sety(real(s));	
	s=0.0;
	for(i=0;i<(gridx);i++)
		for(j=0;j<(gridy);j++)
			for(k=0;k<(gridz);k++)
			{
				s+= 1.0*k*dz*f[i][j][k][0][0];
			}
	s=s*dx*dy*dz/mass;
	Com.setz(real(s));
	Com.print();
	return Com;
}

void velplane(complex<double>*****f,complex<double>*ret,int nmax,int gridx,int gridy,int gridz,double dx,double dy, double dz,double angle)//A function that recombines the distribution function and plots it in in a plane specified by the angle theta.....eg vx/vy plane is theta = pi/2
{
  complex<double> ii(0.0,1.0);
  const double pi = acos(-1.0); 
  complex<double> norm;
  int n,m,i;
  complex<double> inte1,inte2;
  complex<double> sum0,summ;
  double temp;
  for(i=0;i<100;i++)
    {
      sum0 = 0.0;
      summ = 0.0;
      for(n=0;n<nmax;n++)
	{
	  temp = plgndr(n,0,cos(angle));
	  // inte1 = integrate3d(f,gridx,gridy,gridz,dx,dy,dz,n,0);
	  sum0 += f[gridx/2][gridy/2][gridz/2][n][0]*temp;
	  for(m=1;m<=n;m++)
	    {
	      temp = plgndr(n,m,cos(angle));
	      // inte1 = integrate3d(f,gridx,gridy,gridz,dx,dy,dz,n,m);
	      summ += temp*((2.0*real(f[gridx/2][gridy/2][gridz/2][n][m])*cos(2.0*m*pi*i/100.0)) - (2.0*imag(f[gridx/2][gridy/2][gridz/2][n][m])*sin(2.0*m*pi*i/100.0)));
	    }
	}
      
      ret[i] = sum0 + summ;
     
    }
  
}

void velplanelong(complex<double>*****f,complex<double>*ret,int nmax,int gridx,int gridy,int gridz,double dx,double dy, double dz,double angle)//A function that recombines the distribution function and plots it in in a plane specified by the angle theta.....eg vx/vy plane is theta = pi/2...THIS FUNCTION IS BASICALLY THE SAME AS THE ONE ABOVE EXCEPT!!!! IT LOOPS OVER !!THETA!! rather than THI!!!!! and the angle determines the plane in THI!
{
  complex<double> ii(0.0,1.0);
  const double pi = acos(-1.0); 
  complex<double> norm;
  int i,n,m;
  complex<double> inte1,inte2;
  complex<double> sum0,summ;
  double temp;
  for(i=0;i<100;i++)
    {
      sum0 = 0.0;
      summ = 0.0;
      for(n=0;n<nmax;n++)
	{
	  temp = plgndr(n,0,cos(1.0*pi*i/100));
	  inte1 = integrate3d(f,gridx,gridy,gridz,dx,dy,dz,n,0);
	  sum0 += inte1*temp;
	  // sum0 += f[gridx/2][gridy/2][gridz/2][n][0]*temp;
	  for(m=1;m<=n;m++)
	    {
	      temp = plgndr(n,m,cos(1.0*pi*i/100));
	      inte1 = integrate3d(f,gridx,gridy,gridz,dx,dy,dz,n,m);
	      summ += temp*(2.0*real(inte1)*cos(1.0*m*angle) - (2.0*imag(inte1)*sin(1.0*m*angle)));
	      // summ += temp*((2.0*real(f[gridx/2][gridy/2][gridz/2][n][m])*cos(1.0*m*angle)) - (2.0*imag(f[gridx/2][gridy/2][gridz/2][n][m])*sin(1.0*m*angle)));
	    }
	}
      
      ret[i] = sum0 + summ;
     
    }
  
}

/*Function that reconstructs the distribution function from the spherical harmonic expansion, currently only in the theta and thi dimensions
nmax = number of terms in the expansion
gridx,y,z the point in coordintate space at which you want the function reformed
t - timestep number, for file numbering purposes*/

void thetathirecon(complex<double>*****f,int nmax,int gridx,int gridy,int gridz,int t)
{
  int i,j,k;
  complex<double> sum;
  const double pi = acos(-1.0);
  char filename[35];
  sprintf(filename,"ThetaThi%i.dat",t);
  ofstream outfile;
  outfile.open(filename);
  for(int j=0;j<=100;j++)
    {
      for(int k=0;k<=100;k++)
	{
	  sum = (0.0,0.0);
	  for(int n=0;n<nmax;n++)
	    {
	      sum += 1.0*plgndr(n,0,cos(k*pi/100))*(real(f[gridx][gridy][gridz][n][0]));
	      for(int m=1;m<=n;m++)
		{
		  sum += 2.0*plgndr(n,m,(cos(k*pi/100)))*((real(f[gridx][gridy][gridz][n][m])*cos(m*2.0*j*pi/100)) - (imag(f[gridx][gridy][gridz][n][m])*sin(m*2.0*j*pi/100)));
		      }
	    }
	  outfile << real(sum) << '\t';
	      }
      outfile << endl;
    }
  outfile.close(); 
}

double fftanalysis(complex<double>*****f,threevector grid,double t)
{
    ofstream datfile;
    datfile.open("FFTanalysis.dat",ios::out | ios::app);
    int Nx = int(grid.getx()); 
    int Ny = int(grid.gety());
    int Nz = int(grid.getz());
    int i,j,k;  
    fftwnd_plan p;
    fftw_complex in[Nx][Ny][Nz];
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	    for(k=0;k<Nz;k++)
	    {
		in[i][j][k].re = real(f[i][j][k][0][0])-1.0;
		in[i][j][k].im = 0.0;
	    }
    p = fftw3d_create_plan(Nx,Ny,Nz,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data 
    fftwnd_one(p,&in[0][0][0],NULL);  
    fftwnd_destroy_plan(p);
    cout << sqrt((in[1][0][0].re*in[1][0][0].re) + (in[0][1][0].im*in[0][1][0].im))/(Nx*Ny) << endl;
    datfile << t << '\t' << (in[1][0][0].re) << '\t' << (in[1][0][0].im) << '\t' << (in[0][1][0].re) << '\t' <<(in[0][1][0].im) << endl;
    datfile.close();
    return (sqrt((in[1][0][0].re*in[1][0][0].re) + (in[0][1][0].im*in[0][1][0].im))/(Nx*Ny));
}
    
double fftanalysis2d(complex<double>*****f,threevector grid,double t)
{
    ofstream datfile;
    datfile.open("FFTanalysis.dat",ios::out | ios::app);
    int Nx = int(grid.getx()); 
    int Ny = int(grid.gety());
    int i,j,k;  
    fftwnd_plan p;
    fftw_complex in[Nx][Ny];
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	{
	    in[i][j].re = real(f[i][j][0][0][0])-1.0;
	    in[i][j].im = 0.0;
	}
    p = fftw2d_create_plan(Nx,Ny,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data 
    fftwnd_one(p,&in[0][0],NULL);  
    fftwnd_destroy_plan(p);
    cout << sqrt((in[1][0].re*in[1][0].re) + (in[1][0].im*in[1][0].im))/(Nx*Ny) << endl;
    datfile << t << '\t' << (in[1][0].re) << '\t' << (in[1][0].im) << '\t' << (in[0][1].re) << '\t' <<(in[0][1].im) << endl;
    datfile.close();
    return (sqrt((in[1][0].re*in[1][0].re) + (in[0][1].im*in[0][1].im))/(Nx*Ny));
}
        
double fftanalysis2dXZ(complex<double>*****f,threevector grid,double t)
{
    ofstream datfile;
    datfile.open("FFTanalysis.dat",ios::out | ios::app);
    int Nx = int(grid.getx()); 
    int Nz = int(grid.getz());
    int i,j,k;  
    fftwnd_plan p;
    fftw_complex in[Nx][Nz];
    for(i=0;i<Nx;i++)
	for(j=0;j<Nz;j++)
	{
	    in[i][j].re = real(f[i][0][j][0][0])-1.0;
	    in[i][j].im = 0.0;
	}
    p = fftw2d_create_plan(Nx,Nz,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data 
    fftwnd_one(p,&in[0][0],NULL);  
    fftwnd_destroy_plan(p);
    cout << sqrt((in[1][0].re*in[1][0].re) + (in[1][0].im*in[1][0].im)) << endl;
    datfile << t << '\t' << (in[1][0].re) << '\t' << (in[1][0].im) << '\t' << (in[0][1].re) << '\t' <<(in[0][1].im) << endl;
    datfile.close();
    return (sqrt((in[1][0].re*in[1][0].re) + (in[0][1].im*in[0][1].im)));
}
