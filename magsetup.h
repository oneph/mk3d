#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <fftw.h>
#include <rfftw.h>

void magfieldsetup(threevector ***B,int gridz,double dz,double Lz,int n,double Bo,double deltaB,char magdirn,char addirn)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation
{
  int i,j,k,h;
  double pi = acos(-1.0);
  threevector Equil(0.0,0.0,0.0);
  switch(magdirn)//Direction of the equilibrium magnetic field
    {
    case 'x':
      Equil.setx(1.0);//This is the equilibrium magnetic field
      cout << "Equilibrium field in the X direction" << endl;
      break;
    case 'y':
      Equil.sety(1.0);
      cout << "Equilibrium field in the Y direction" << endl;
      break;
    case 'z':
      Equil.setz(1.0);
      cout << "Equilibrium field in the Z direction" << endl;
      break;
    default:
      cout << "Direction not specified in magfieldsetup" << endl;
      terminate();
    }
  
  double *pert;
  
  switch(addirn)
    {
    case 'x':
      cout << "Advection in the x-direction" << endl;
      pert = new double [gridz];
      for(k=0;k<gridz;k++)
	{
	  for(h=1;h<=1;h++)
	    {
	      pert[k] += (deltaB)*cos(2.0*pi*dz*(k+0.5)*(n)/Lz);
	    }
	  B[k][0][0] = Equil*(Bo + pert[k]); 
	}
      break;
    case 'y':
      cout << "Advection in the y-direction" << endl;
      pert = new double [gridz];
      for(k=0;k<gridz;k++)
	{
	  for(h=1;h<=1;h++)
	    {
	      pert[k] += (deltaB)*cos(2.0*pi*dz*(k+0.5)*(n)/Lz);
	    }
	  B[0][k][0] = Equil*(Bo + pert[k]); 
	}
      break;
    case 'z':
      cout << "Advection in the z-direction" << endl;
      pert = new double [gridz];
      for(k=0;k<gridz;k++)
	{
	  for(h=1;h<=1;h++)
	    {
	      pert[k] += (deltaB)*cos(2.0*pi*dz*(k+0.5)*(n)/Lz);
	    }
	  B[0][0][k] = Equil*(Bo + pert[k]); 
	}
      break;
    default:
      cout << "Wrong advection direction in magfieldsetup" << endl;
      terminate();
    }
  delete [] pert;
}

void magfieldsetup2d(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double deltaB)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation
{
  int i,j,k,h;
  double pi = acos(-1.0);
  double phase1,phase2,phase3,phase4,phase5,phase6;
  cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
  cout << "Equilibrium field in the Z direction" << endl;	  
  phase1=0.0;//2.0*pi*rand()/RAND_MAX;
  phase2=0.0;//2.0*pi*rand()/RAND_MAX;
  phase3=0.0;//2.0*pi*rand()/RAND_MAX;
  phase4=0.0;//2.0*pi*rand()/RAND_MAX;
  phase5=0.0;//2.0*pi*rand()/RAND_MAX;
  phase6=0.0;//2.0*pi*rand()/RAND_MAX;
  for(j=0;j<=(gri.getx());j++)
      for(k=0;k<=(gri.getz());k++)
      {
	  A[j][0][k].sety((Bo*ele.getx()*j) + ((deltaB)*sin(2.0*pi*j*ele.getx()+phase1)+(deltaB/4.0)*sin(4.0*pi*j*ele.getx()+phase2)+(deltaB/16.0)*sin(8.0*pi*j*ele.getx()+phase3))/(2.0*pi));
      }
  curl3d(B,A,sys,ele,gri);
}

void magfieldsetup1dturb(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double dB1)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation
{ 
    ofstream outfile("out.dat");
    ofstream prefilter("prefil.dat");
    int i,j,k,h;
    double pi = acos(-1.0);
    cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
    cout << "Equilibrium field in the Z direction" << endl;	  
    double ra;
    int Nx =int(gri.getx()) + 1;
    fftw_real in[Nx],out[Nx] ;
    rfftw_plan p;
    for(k=0;k<Nx;k++)
      {
	in[k] = 0.0;
	out[k] = 0.0;
      }
    for(k=1;k<=20;k++)
      {
	double ran=1.0*rand()/RAND_MAX - 0.5;
	in[k+1]=ran;
	in[Nx-1-k]=in[k+1];
      }
    p  = rfftw_create_plan(Nx,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE);
    for(k=0;k<Nx;k++)
      {
	prefilter << in[k] << endl;
      }
    rfftw_one(p,in,out);
    rfftw_destroy_plan(p);
    for(k=0;k<Nx;k++)
      A[k][0][0].sety(Bo*(k)*ele.getx() + dB1*out[k]/Nx);
    curl3d(B,A,sys,ele,gri);
}

void magfieldsetup2dturb(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double dB1)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation
{ 
    ofstream specfileRe("SpectrumRe.dat");
    ofstream specfileIm("SpectrumIm.dat");
    int i,j,k,h;
    double pi = acos(-1.0);
    cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
    cout << "Equilibrium field in the Z direction" << endl;	  
    double ra;
    int Nx =int(gri.getx()) + 1;
    int Nz =int(gri.getz()) + 1;
    fftw_complex in[Nx][Nz];
    fftwnd_plan p;
    for(i=0;i<Nx;i++)//This series of loops initialises the grid with a random number between -1 and 1
      {
	for(j=0;j<Nz;j++)
	  {
	    in[i][j].re=0.0;
	    in[i][j].im=0.0;
	    double ran=2.0*rand()/RAND_MAX - 1.0;
	    in[i][j].re=ran;
	  }
      }

    p = fftw2d_create_plan(Nx,Nz,FFTW_FORWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Forward fourier transform to get spectral information
    fftwnd_one(p,&in[0][0],NULL);
    
    for(i=0;i<=Nx/2;i++)//Filter out all high frequency modes
	for(j=0;j<=Nz/2;j++)
	{ 
	    if((i==0)||(j==0))
	    {
		in[i][j].re=0.0;
		in[i][j].im=0.0;
		in[0][Nz-j].re=0.0;
		in[0][Nz-j].im=0.0;
		in[Nx-i][0].re=0.0;
		in[Nx-i][0].im=0.0;
	    }else
	    {
		if(sqrt((i-1)*(i-1)+(j-1)*(j-1))>10)
		{
		    in[i][j].re = 0.0;
		    in[i][j].im = 0.0;
		    in[Nx-i][Nz-j]=in[i][j];
		    in[Nx-i][j]=in[i][j];
		    in[i][Nz-j]=in[i][j];
		}
	    }
	}
    for(i=0;i<=Nx;i++)
    {
	for(j=0;j<=Nz;j++)
	{
	    specfileRe << in[i][j].re << '\t';
	    specfileIm << in[i][j].im << '\t'; 
	}
	specfileRe << endl;specfileIm << endl;
    }
    
    p = fftw2d_create_plan(Nx,Nz,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data
    fftwnd_one(p,&in[0][0],NULL);
    fftwnd_destroy_plan(p);
    for(i=0;i<Nx;i++)
      for(j=0;j<Nz;j++)
	A[i][0][j].sety(Bo*(i)*ele.getx() + dB1*in[i][j].re/(Nx*Nz));
    curl3d(B,A,sys,ele,gri);//Curl of A is the B field
}

void magfieldsetup2dturb2(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double dB1)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation ... This version of setup2Dturb, uses a flat spectrum and random phase.
{ 
    ofstream specfileRe("SpectrumRe.dat");
    ofstream specfileIm("SpectrumIm.dat");
    int i,j,k,h;
    double pi = acos(-1.0);
    cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
    cout << "Equilibrium field in the Z direction" << endl;	  
    double ra;
    int Nx =int(gri.getx());
    int Nz =int(gri.getz());
    fftw_complex in[Nx][Nz];
    fftwnd_plan p;
    for(i=0;i<Nx;i++)
	for(j=0;j<Nz;j++)
	{
	    in[i][j].re = 0.0;
	    in[i][j].im = 0.0; 
	}
    double randphase,weight;
    for(i=0;i<=Nx/2;i++)//This series of loops initialises the grid with a random number between -1 and 1
      {
	for(j=0;j<=Nz/2;j++)
	{ 
	    double r = sqrt((i*i)+(j*j));
	    if(r==0)
		weight = 1.0;
	    else(weight = 1.0*(pow(1.0*r,11.0/6.0)));
	    { 
		if(sqrt(((i*i)+(j*j)))<Nx/2)
		{
		    randphase=2.0*pi*rand()/RAND_MAX; 		  
		    in[i][j].re = (0.5/weight)*cos(randphase);
		    in[i][j].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[i][Nz-j].re = (0.5/weight)*cos(randphase);
		    in[i][Nz-j].im = (0.5/weight)*sin(randphase);
		}else
		{
		    in[i][j].re=0.0;
		    in[i][j].im=0.0;
		}
	    }
	  }
      }
    in[0][0].re = 0.0;
    in[0][0].im = 0.0;
    for(i=0;i<Nx;i++)
    {
	for(j=0;j<Nz;j++)
	{
	    specfileRe << in[i][j].re << '\t';
	    specfileIm << in[i][j].im << '\t'; 
	}
	specfileRe << endl;specfileIm << endl;
    }
    
    
    p = fftw2d_create_plan(Nx,Nz,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data
    fftwnd_one(p,&in[0][0],NULL);
    fftwnd_destroy_plan(p);
    for(i=0;i<Nx;i++)
    {
	for(j=0;j<Nz;j++)
	{
	    specfileRe << in[i][j].re << '\t';
	    specfileIm << in[i][j].im << '\t'; 
	}
	specfileRe << endl;specfileIm << endl;
    }
     
    for(i=0;i<Nx;i++)
      for(j=0;j<Nz;j++)
      	  A[i][0][j].sety(Bo*(i+0.5)*ele.getx() + dB1*in[i][j].re);
    for(i=0;i<Nx;i++)
	A[i][0][Nz].sety(Bo*(i+0.5)*ele.getx() + dB1*in[i][0].re);
    for(j=0;j<Nz;j++)
	A[Nx][0][j].sety(Bo*(i+0.5)*ele.getx() + dB1*in[0][j].re);//Ensures periodicity of the magnetic field profile
    A[Nx][0][Nz].sety(Bo*(Nx)*ele.getx() + dB1*in[0][0].re);
    
    curl2dXZ(B,A,sys,ele,gri);//Curl of A is the B field
 
}


void magfieldsetup3dIso(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double dB1)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation ... This version of setup2Dturb, uses a flat spectrum and random phase.
{  
    srand((unsigned)time(NULL));
    ofstream specfileRe("SpectrumRe.dat");
    ofstream specfileIm("SpectrumIm.dat");
    int i,j,k,h;
    double pi = acos(-1.0);
    cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
    cout << "Equilibrium field in the Z direction" << endl;	  
    double ra;
    int Nx =int(gri.getx());
    int Ny =int(gri.gety());
    int Nz =int(gri.getz());
    fftw_complex in[Nx][Ny][Nz];
    fftw_complex in2[Nx][Ny][Nz];
    fftw_complex in3[Nx][Ny][Nz];
    fftwnd_plan p;
    double randphase,weight;

    for(i=0;i<Nx/2;i++)//This series of loops initialises the grid with a random number between -1 and 1
	for(j=0;j<Ny/2;j++)
	    for(k=0;k<Nz/2;k++)
	    {
		double r = sqrt((i*i)+(j*j)+(k*k));
		if((r<9*Nx/20))
		{
		    if(r==0)
			weight = 1.0;
		    else(weight = 1.0*(pow(1.0*r,11.0/6.0)));//Here you can change the spectral index.....8/3 corresponds to Kolmogorov
		    //  cout << weight << endl;
		    randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the x component of the vector potential
		    in[i][j][k].re = (0.5/weight)*cos(randphase);
		    in[i][j][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[Nx-i-1][j][k].re = (0.5/weight)*cos(randphase);
		    in[Nx-i-1][j][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[i][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in[i][Ny-j-1][k].im = (0.5/weight)*sin(randphase);	   
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[i][j][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in[i][j][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[Nx-i-1][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in[Nx-i-1][Ny-j-1][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[Nx-i-1][j][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in[Nx-i-1][j][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[i][Ny-j-1][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in[i][Ny-j-1][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in[Nx-i-1][Ny-j-1][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in[Nx-i-1][Ny-j-1][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    
		    randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the y component of the vector potential
		    in2[i][j][k].re = (0.5/weight)*cos(randphase);
		    in2[i][j][k].im = (0.5/weight)*sin(randphase);		    
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[Nx-i-1][j][k].re = (0.5/weight)*cos(randphase);
		    in2[Nx-i-1][j][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[i][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in2[i][Ny-j-1][k].im = (0.5/weight)*sin(randphase);	  
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[i][j][Nz-k].re = (0.5/weight)*cos(randphase);
		    in2[i][j][Nz-k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[Nx-i-1][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in2[Nx-i-1][Ny-j-1][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[Nx-i-1][j][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in2[Nx-i-1][j][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[i][Ny-j-1][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in2[i][Ny-j-1][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in2[Nx-i-1][Ny-j-1][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in2[Nx-i-1][Ny-j-1][Nz-k-1].im = (0.5/weight)*sin(randphase);
				   
		    randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the z component of the vector potential
		    in3[i][j][k].re = (0.5/weight)*cos(randphase);
		    in3[i][j][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[Nx-i-1][j][k].re = (0.5/weight)*cos(randphase);
		    in3[Nx-i-1][j][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[i][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in3[i][Ny-j-1][k].im = (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[i][j][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in3[i][j][Nz-k-1].im =  (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[Nx-i-1][Ny-j-1][k].re = (0.5/weight)*cos(randphase);
		    in3[Nx-i-1][Ny-j-1][k].im =  (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[Nx-i-1][j][Nz-k-1].re =  (0.5/weight)*cos(randphase);
		    in3[Nx-i-1][j][Nz-k-1].im =  (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[i][Ny-j-1][Nz-k-1].re =  (0.5/weight)*cos(randphase);
		    in3[i][Ny-j-1][Nz-k-1].im =  (0.5/weight)*sin(randphase);
		    randphase=2.0*pi*rand()/RAND_MAX;
		    in3[Nx-i-1][Ny-j-1][Nz-k-1].re = (0.5/weight)*cos(randphase);
		    in3[Nx-i-1][Ny-j-1][Nz-k-1].im = (0.5/weight)*sin(randphase);
		    
		}else
		{    		    	   
		    in[i][j][k].re = 0.0;
		    in[i][j][k].im = 0.0;
		    in[Nx-i-1][j][k].re = 0.0;
		    in[Nx-i-1][j][k].im = 0.0;
		    in[i][Ny-j-1][k].re = 0.0;
		    in[i][Ny-j-1][k].im = 0.0;
		    in[i][j][Nz-k-1].re = 0.0;
		    in[i][j][Nz-k-1].im = 0.0;
		    in[Nx-i-1][Ny-j-1][k].re = 0.0;
		    in[Nx-i-1][Ny-j-1][k].im = 0.0;
		    in[Nx-i-1][j][Nz-k-1].re = 0.0;
		    in[Nx-i-1][j][Nz-k-1].im = 0.0;
		    in[i][Ny-j-1][Nz-k-1].re = 0.0;
		    in[i][Ny-j-1][Nz-k-1].im = 0.0;
		    in[Nx-i-1][Ny-j-1][Nz-k-1].re = 0.0;
		    in[Nx-i-1][Ny-j-1][Nz-k-1].im = 0.0;
		    in2[i][j][k].re = 0.0;
		    in2[i][j][k].im = 0.0;
		    in2[Nx-i-1][j][k].re = 0.0;
		    in2[Nx-i-1][j][k].im = 0.0;
		    in2[i][Ny-j-1][k].re = 0.0;
		    in2[i][Ny-j-1][k].im = 0.0;
		    in2[i][j][Nz-k-1].re = 0.0;
		    in2[i][j][Nz-k-1].im = 0.0;
		    in2[Nx-i-1][Ny-j-1][k].re = 0.0;
		    in2[Nx-i-1][Ny-j-1][k].im = 0.0;
		    in2[Nx-i-1][j][Nz-k-1].re = 0.0;
		    in2[Nx-i-1][j][Nz-k-1].im = 0.0;
		    in2[i][Ny-j-1][Nz-k-1].re = 0.0;
		    in2[i][Ny-j-1][Nz-k-1].im = 0.0;
		    in2[Nx-i-1][Ny-j-1][Nz-k-1].re = 0.0;
		    in2[Nx-i-1][Ny-j-1][Nz-k-1].im = 0.0;
		    in3[i][j][k].re = 0.0;
		    in3[i][j][k].im = 0.0;
		    in3[Nx-i-1][j][k].re = 0.0;
		    in3[Nx-i-1][j][k].im = 0.0;
		    in3[i][Ny-j-1][k].re = 0.0;
		    in3[i][Ny-j-1][k].im = 0.0;
		    in3[i][j][Nz-k-1].re = 0.0;
		    in3[i][j][Nz-k-1].im = 0.0;
		    in3[Nx-i-1][Ny-j-1][k].re = 0.0;
		    in3[Nx-i-1][Ny-j-1][k].im = 0.0;
		    in3[Nx-i-1][j][Nz-k-1].re = 0.0;
		    in3[Nx-i-1][j][Nz-k-1].im = 0.0;
		    in3[i][Ny-j-1][Nz-k-1].re = 0.0;
		    in3[i][Ny-j-1][Nz-k-1].im = 0.0;
		    in3[Nx-i-1][Ny-j-1][Nz-k-1].re = 0.0;
		    in3[Nx-i-1][Ny-j-1][Nz-k-1].im = 0.0;
		}
	    } 
    in[0][0][0].re =0.0;
    in[0][0][0].im =0.0;
    in2[0][0][0].re =0.0;
    in2[0][0][0].im =0.0;
    in3[0][0][0].re =0.0;
    in3[0][0][0].im =0.0;
    for(i=0;i<Nx;i++)
    {
	for(j=0;j<Nz;j++)
	{
	    specfileRe << in[i][0][j].re << '\t';
	    specfileIm << in[i][0][j].im << '\t'; 
	}
	specfileRe << endl;specfileIm << endl;
    }
    p = fftw3d_create_plan(Nx,Ny,Nz,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data 
    fftwnd_one(p,&in[0][0][0],NULL);  
    fftwnd_one(p,&in2[0][0][0],NULL);
    fftwnd_one(p,&in3[0][0][0],NULL);
    fftwnd_destroy_plan(p);
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	    for(k=0;k<Nz;k++)
	    {
		A[i][j][k].sety(Bo*(i)*ele.getx() + dB1*in[i][j][k].re);
		A[i][j][k].setx(dB1*in2[i][j][k].re);
		A[i][j][k].setz(dB1*in3[i][j][k].re);
	    }
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	{
	    A[i][j][Nz].sety(Bo*(i)*ele.getx() + dB1*in[i][j][0].re);
	    A[i][j][Nz].setx(dB1*in2[i][j][0].re);
	    A[i][j][Nz].setz(dB1*in3[i][j][0].re);
	}  
    for(j=0;j<Ny;j++)
	for(k=0;k<Nz;k++)
	{
	    A[Nx][j][k].sety(Bo*(Nx)*ele.getx() + dB1*in[0][j][k].re);
	    A[Nx][j][k].setx(dB1*in2[0][j][k].re);
	    A[Nx][j][k].setz(dB1*in3[0][j][k].re);
	}  
    for(i=0;i<Nx;i++)
	for(k=0;k<Nz;k++)
	{
	    A[i][Ny][k].sety(Bo*(i)*ele.getx() + dB1*in[i][0][k].re);
	    A[i][Ny][k].setx(dB1*in2[i][0][k].re);
	    A[i][Ny][k].setz(dB1*in3[i][0][k].re);
	}
    for(i=0;i<Nx;i++)
    {
	A[i][Ny][Nz].sety(Bo*(i)*ele.getx() + dB1*in[i][0][0].re);
	A[i][Ny][Nz].setx(dB1*in2[i][0][0].re);
	A[i][Ny][Nz].setz(dB1*in3[i][0][0].re);
    }
    for(j=0;j<Ny;j++)
    {
	A[Nx][j][Nz].sety(Bo*(Nx)*ele.getx() + dB1*in[0][j][0].re);
	A[Nx][j][Nz].setx(dB1*in2[0][j][0].re);
	A[Nx][j][Nz].setz(dB1*in3[0][j][0].re);
    }
    for(k=0;k<Nz;k++)
    {
	A[Nx][Ny][k].sety(Bo*(Nx)*ele.getx() + dB1*in[0][0][k].re);
	A[Nx][Ny][k].setx(dB1*in2[0][0][k].re);
	A[Nx][Ny][k].setz(dB1*in3[0][0][k].re);
    }
    A[Nx][Ny][Nz].sety(Bo*(Nx)*ele.getx() + dB1*in[0][0][0].re);
    A[Nx][Ny][Nz].setx(dB1*in2[0][0][0].re);
    A[Nx][Ny][Nz].setz(dB1*in3[0][0][0].re);
    curl3d(B,A,sys,ele,gri);//Curl of A is the B field
}


void magfieldsetup3d2Dgeometry(threevector ***B,threevector ***A,threevector sys,threevector ele,threevector gri,double Bo,double dB1)//Delta B is the magnetic perturbation size, while n is the number of harmonics in the perturbation ... This version of setup2Dturb, uses a flat spectrum and random phase.
{    
    srand( (unsigned)time( NULL ) );
    ofstream specfileRe("SpectrumRe.dat");
    ofstream specfileIm("SpectrumIm.dat");
    int i,j,k,h;
    double pi = acos(-1.0);
    cout << gri.getx() << " " << gri.gety() << " " << gri.getz() << endl;
    cout << "Equilibrium field in the Z direction - 2D Geometry" << endl;	  
    double ra;
    int Nx =int(gri.getx());
    int Ny =int(gri.gety());
    int Nz =int(gri.getz());
    fftw_complex in[Nx][Ny];
    fftw_complex in2[Nx][Ny];
    fftw_complex in3[Nx][Ny];
    fftwnd_plan p;
    double randphase,weight;

    for(i=0;i<Nx/2;i++)//This series of loops initialises the grid with a random number between -1 and 1
	for(j=0;j<Ny/2;j++)
	{
	    double r = sqrt((i*i)+(j*j));
	    if((r<9*Nx/20))
	    {
		if(r==0)
		    weight = 1.0;
		else(weight = 1.0*(pow(1.0*r,8.0/3.0)));//Here you can change the spectral index.....8/3 corresponds to Kolmogorov
		//  cout << weight << endl;
		randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the x component of the vector potential
		in[i][j].re = (0.5/weight)*cos(randphase);
		in[i][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in[Nx-i-1][j].re = (0.5/weight)*cos(randphase);
		in[Nx-i-1][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in[i][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in[i][Ny-j-1].im = (0.5/weight)*sin(randphase);	   
		randphase=2.0*pi*rand()/RAND_MAX;
		in[Nx-i-1][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in[Nx-i-1][Ny-j-1].im = (0.5/weight)*sin(randphase);

		randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the x component of the vector potential
		in2[i][j].re = (0.5/weight)*cos(randphase);
		in2[i][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in2[Nx-i-1][j].re = (0.5/weight)*cos(randphase);
		in2[Nx-i-1][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in2[i][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in2[i][Ny-j-1].im = (0.5/weight)*sin(randphase);	   
		randphase=2.0*pi*rand()/RAND_MAX;
		in2[Nx-i-1][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in2[Nx-i-1][Ny-j-1].im = (0.5/weight)*sin(randphase);

		randphase=2.0*pi*rand()/RAND_MAX;//Initialisation of the x component of the vector potential
		in3[i][j].re = (0.5/weight)*cos(randphase);
		in3[i][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in3[Nx-i-1][j].re = (0.5/weight)*cos(randphase);
		in3[Nx-i-1][j].im = (0.5/weight)*sin(randphase);
		randphase=2.0*pi*rand()/RAND_MAX;
		in3[i][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in3[i][Ny-j-1].im = (0.5/weight)*sin(randphase);	   
		randphase=2.0*pi*rand()/RAND_MAX;
		in3[Nx-i-1][Ny-j-1].re = (0.5/weight)*cos(randphase);
		in3[Nx-i-1][Ny-j-1].im = (0.5/weight)*sin(randphase);
			    
	    }else
	    {    		    	   
		in[i][j].re = 0.0;
		in[i][j].im = 0.0;
		in[Nx-i-1][j].re = 0.0;
		in[Nx-i-1][j].im = 0.0;
		in[i][Ny-j-1].re = 0.0;
		in[i][Ny-j-1].im = 0.0;
		in[Nx-i-1][Ny-j-1].re = 0.0;
		in[Nx-i-1][Ny-j-1].im = 0.0;

		in2[i][j].re = 0.0;
		in2[i][j].im = 0.0;
		in2[Nx-i-1][j].re = 0.0;
		in2[Nx-i-1][j].im = 0.0;
		in2[i][Ny-j-1].re = 0.0;
		in2[i][Ny-j-1].im = 0.0;
		in2[Nx-i-1][Ny-j-1].re = 0.0;
		in2[Nx-i-1][Ny-j-1].im = 0.0;
		
		in3[i][j].re = 0.0;
		in3[i][j].im = 0.0;
		in3[Nx-i-1][j].re = 0.0;
		in3[Nx-i-1][j].im = 0.0;
		in3[i][Ny-j-1].re = 0.0;
		in3[i][Ny-j-1].im = 0.0;
		in3[Nx-i-1][Ny-j-1].re = 0.0;
		in3[Nx-i-1][Ny-j-1].im = 0.0;
	    }
	} 
    in[0][0].re =0.0;
    in[0][0].im =0.0;
    in2[0][0].re =0.0;
    in2[0][0].im =0.0;
    in3[0][0].re =0.0;
    in3[0][0].im =0.0;
    for(i=0;i<Nx;i++)
    {
	for(j=0;j<Ny;j++)
	{
	    specfileRe << in[i][j].re << '\t';
	    specfileIm << in[i][j].im << '\t';
	}
	specfileRe << endl;specfileIm << endl;
    }
    p = fftw2d_create_plan(Nx,Ny,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE);//Inverse fourier transform, to calculate real data 
    fftwnd_one(p,&in[0][0],NULL);  
    fftwnd_one(p,&in2[0][0],NULL);
    fftwnd_one(p,&in3[0][0],NULL);
    fftwnd_destroy_plan(p);
 
    
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	{
	    A[i][j][0].sety((Bo*(i)*ele.getx()) + dB1*in[i][j].re);
	    A[i][j][0].setx(dB1*in2[i][j].re);
	    A[i][j][0].setz(dB1*in3[i][j].re); 
	}
    for(i=0;i<Nx;i++)
    {
	A[i][Ny][0].sety((Bo*(i)*ele.getx()) + dB1*in[i][0].re);
	A[i][Ny][0].setx(dB1*in2[i][0].re);
	A[i][Ny][0].setz(dB1*in3[i][0].re); 
    }
    for(i=0;i<Ny;i++)
    {
	A[Nx][i][0].sety((Bo*(Nx)*ele.getx()) + dB1*in[0][i].re);
	A[Nx][i][0].setx(dB1*in2[0][i].re);
	A[Nx][i][0].setz(dB1*in3[0][i].re); 
    }
    A[Nx][Ny][0].sety((Bo*(Nx)*ele.getx()) + dB1*in[0][0].re);
    A[Nx][Ny][0].setx(dB1*in2[0][0].re);
    A[Nx][Ny][0].setz(dB1*in3[0][0].re); 
    curl2dXY(B,A,sys,ele,gri);//Curl of A is the B field
  
}

void uniformB(threevector ***B,threevector gri,double Bo)
{
    int i,j,k;
    int Nx=int(gri.getx());
    int Ny=int(gri.gety());
    int Nz=int(gri.getz());
    for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	    for(k=0;k<Nz;k++)
	    {
		B[i][j][k].setx(0.0);
		B[i][j][k].sety(0.0);
		B[i][j][k].setz(Bo);
	    }
}
