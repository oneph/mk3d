//A program to simualate Tokamak transport using the KALOS algorithm
//Author:William A. Hornsby
//Date:16/1/2005
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>
#include "threevector.h"
#include "advection.h"
#include "magfile.h"   
#include "diagnostics.h"
#include "VTKfiles.h"
#include "collisions.h"
#include "magsetup.h"
#include "datastruct.h"
//#include <mpi.h>

using namespace std;

int main()
{   
    srand( (unsigned)time( NULL ) );
    ofstream Bxfile("BxComp.dat");
    ofstream Byfile("ByComp.dat");
    ofstream Bzfile("BzComp.dat");
    ofstream Bmagfile("BmagXY.dat");  
    ofstream Axfile("AxComp.dat");
    ofstream Ayfile("AyComp.dat");
    ofstream Azfile("AzComp.dat");
    ofstream compfile("BcompXYZ.dat");
    ofstream datfile("RunParameters.txt");
    ofstream harmonicfile("HarmonicAmplitudes.dat");
    ofstream fftfile("FFTanalysis.dat");
    int i,j,k,n,m,p,testing;
    
    char filename[20];
    cout.precision(8);
    compfile.precision(16);
    
    //ANY PHYSICAL CONSTANTS ARE DEFINED HERE!!!!
    const double pi = acos(-1.0);  
    const complex<double> ii(0.0,1.0);
    //END OF PHYSICAL CONSTANTS!
    
    double Lx,Ly,Lz;//The maximum extent of the spatial dimension
    int gridx,gridy,gridz,ign;
    double t=0;//The time coordinate
    double dx,dy,dz,dt;//Time and spatial increments
  
    double Bmag,dB,Bmax=0.0;
    int count = 0;//Timestep count - Integer number useful for file naming
    double v;//Velocity magnitude
    double nu;//The collisional frequency
    threevector Masscentre;
    complex<double> sum;
    double totdens;
    time_t curr; 
    //Following code loads parameters from the file (inputdeck.txt)created in external script 
    char *str;
    str = new char [6];
    double num;
    int num1;
    ifstream finput("inputdeck.txt");
    if(!finput)
    {
	  cout << "Inputdeck.txt error - Program aborted" << endl;
	  return -1;
    }
   
    finput >> str >> num1;
    cout << str << " " << num1  << endl;
    gridx = num1;
    finput >> str >> num1;
    cout << str << " " << num1  << endl;
    gridy = num1;
    finput >> str >> num1;
    cout << str << " " << num1  << endl;
    gridz = num1;  
    
    finput >> str >> num;
    cout << str << " " << num << endl;
    Lx = num;
    finput >> str >> num;
    cout << str << " " << num << endl;
    Ly = num;
    finput >> str >> num;
    cout << str << " " << num << endl;
    Lz = num;
    
    finput >> str >> num;
    cout << str << " " << num << endl;
    nu = num;
    finput >> str >> num;
    cout << str << " " << num << endl;
    v = num;  
    finput >> str >> num1;
    cout << str << " " << num1  << endl;
    const int nmax = num1;//The number of harmonics in the expansion !!!!Expansion goes from 0 to nmax-1!!!!
    
    finput >> str >> num;
    cout << str << " " << num << endl;
    Bmag = num;
    finput >> str >> num;
    cout << str << " " << num << endl;
    dB = num;  
    finput >> str >> num;
    cout << str << " " << num << endl;
    ign = int(num);
    finput >> str >> num;
    cout << str << " " << num << endl;
    testing = int(num);
   
    //File load ends here
    curr = time(NULL);
    cout << "The date is = " << ctime(&curr);
    threevector Systemsize(Lx,Ly,Lz);
    threevector Gridsize(gridx,gridy,gridz);
    
    //Distribution function - First array is the x position element,second is y, third is z, fourth array is the n Legengre ploynomial number, and m is the harmonic
    threevector ***B = AllocateThreevector3D(gridx,gridy,gridz); 
    threevector ***A = AllocateThreevector3D(gridx+1,gridy+1,gridz+1); 
    complex<double> *****f = AllocateDistroF(gridx,gridy,gridz,nmax);

    //End pf definition of distro function - Multi dimensional arrays = PITA
    //Cartesian Grid!!!
    dx = Lx/gridx;
    dy = Ly/gridy;
    dz = Lz/gridz;
    threevector Elementsize(dx,dy,dz);
    double scale = smallest(dx,smallest(dy,dz));
    cout << "SMALLEST LENGTH-SCALE = " << scale << endl;
    dt = scale/(4.0*v);
    cout << "TIMESTEP SIZE =" << dt << endl;    
    double harmamp;
    //!!!!!! INITIAL CONDTIIONS!!!!!!!! 
    
    //MAGNETIC FIELD INITIALISATION
    if(dB==0.0)
    {
	  uniformB(B,Gridsize,Bmag);
    }else
      {
	    magfieldsetup3dIso(B,A,Systemsize,Elementsize,Gridsize,Bmag,dB);//Initialisation of the magnetic field structure
      }
    //ENDS    

    for(i=0;i<gridx;i++) 
      {
	  for(k=0;k<gridy;k++) 
	    {
	      Bxfile << B[i][k][gridz/2].getx() << '\t';
	      Byfile << B[i][k][gridz/2].gety() << '\t';	  
	      Bzfile << B[i][k][gridz/2].getz() << '\t';	    
	      Bmagfile <<  B[i][k][gridz/2].mag() << '\t';
	    }
	    Bxfile << endl;	
	    Byfile << endl;
	    Bzfile << endl;
	    Bmagfile << endl;
      }  
    
    for(i=0;i<=gridx;i++)
      {
	  for(k=0;k<=gridy;k++) 
	    {
	      Axfile << A[i][k][gridz/2].getx() << '\t';
	      Ayfile << A[i][k][gridz/2].gety() << '\t'; 
	      Azfile << A[i][k][gridz/2].getz() << '\t';
	    }
	    Axfile << endl;	
	    Ayfile << endl;	
	    Azfile << endl;
      }
	//Free up the memory from the vector potential
    FreeThreevector3D(A,gridx+1,gridy+1);
    
    
    //END OF FIELD INITIALISATION.  ONLY THE B-FIELD DATA REMAINS.....VECTOR POTENTIAL IS DELETED
 
    //!!!!!! INITIAL CONDTIIONS!!!!!!!!
    for(i=0;i<gridx;i++)
	  for(j=0;j<gridy;j++)
	    for(k=0;k<gridz;k++)
		  for(m=0;m<nmax;m++)  
		    for(n=m;n<nmax;n++)
		    {
			  f[i][j][k][n][m] = 0.0;
		    }
    
    for(i=0;i<gridx;i++) 
	  for(j=0;j<gridy;j++)
	    for(k=0;k<gridz;k++) 
	    {
	      f[i][j][k][0][0] = 1.0 + 0.01*cos(1.0*pi*(((j+0.5)*dy/Ly)+((i+0.5)*dx/Lx)))*cos(1.0*pi*(((j+0.5)*dy/Ly)-((i+0.5)*dx/Lx)));
	      //f[i][j][k][0][0] = 1.0*exp(-(i*dx - 0.5)*(i*dx - 0.5)/0.01)*exp(-(j*dx - 0.5)*(j*dx - 0.5)/0.01)*exp(-(k*dz - 0.5)*(k*dz - 0.5)/0.01);//Harmonic density perturbation
	      compfile << (i+0.5)*dx << '\t' << (j+0.5)*dy << '\t' << (k+0.5)*dz << '\t' << B[i][j][k].getx() << '\t' << B[i][j][k].gety()  << '\t' << B[i][j][k].getz() << endl;
	      if(B[i][j][k].mag()>Bmax)
		  {
		    Bmax = B[i][j][k].mag();
		  }
	    }

    //End pf definition of distro function - Multi dimensional arrays = PITA
    //END OF INITIAL CONDITIONS
    
    //Code that writes to the Run Parameters file
    datfile << "Run Parameters for KALOS run on" << " " << ctime(&curr)  << endl;
    datfile << "Length of system  , Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << endl;
    datfile << "Number of Grid cells, Nx = " << gridx << ", Ny = " << gridy << ", Nz = " << gridz << endl;
    datfile << "Number of Harmonics = " << nmax << endl;
    datfile << "Thermal velocity, vth =" << " " << v << endl;
    datfile << "Collision frequency = " << nu << endl;
    datfile << "Timestep =" << dt << endl;
    //End of parameter file.
    
    t=0.0;
    i=0;

    double magan = fftanalysis(f,Gridsize,0);
    double maganini = magan;
    cout << "magini="<<" "<< maganini << endl;
    for(m=0;m<nmax;m++)
	for(p=0;p<=m;p++)
	{
	    harmamp=0.0;
	    for(n=0;n<gridx;n++) 
		for(j=0;j<gridy;j++)
		    for(k=0;k<gridz;k++) 
			harmamp+= abs(f[n][j][k][m][p])/(gridx*gridy*gridz);
	    
	    harmonicfile << (harmamp) << '\t';
	}
    harmonicfile << endl;	   
    sprintf(filename,"TDdens%i.dat",0);
    ofstream ddfile;
    ddfile.open(filename);
    ddfile.precision(16);
    for(n=0;n<gridx;n++) 
	  for(j=0;j<gridy;j++)
	    for(k=0;k<gridz;k++) 
		  ddfile << (n+0.5)*dx << '\t' << (j+0.5)*dy << '\t' << (k+0.5)*dx << '\t' << real(f[n][j][k][0][0]) << '\t' << real(f[n][j][k][1][0]) << '\t' << real(f[n][j][k][1][1]) << '\t' << -imag(f[n][j][k][1][1]) << endl;
    
    thetathirecon(f,nmax,0,0,0,0);
    cout << "Initialised!!!" << endl;
    int rotate = 1;
    if(testing==1)
    {
      cout << "Testing Mode : Code quits before time integration" << endl;
      return 0;
    }
    
    do{  
	
	  advectionx(f,nmax,gridx,gridy,gridz,dx,v,dt,'p');
	  magadvance(f,B,nmax,dt/3.0,gridx,gridy,gridz);
	  collision(f,nmax,dt/3.0,gridx,gridy,gridz,nu);  
	  advectiony(f,nmax,gridy,gridx,gridz,dy,v,dt,'p');
	  magadvance(f,B,nmax,dt/3.0,gridx,gridy,gridz);
	  collision(f,nmax,dt/3.0,gridx,gridy,gridz,nu); 
	  advectionz(f,nmax,gridz,gridx,gridy,dz,v,dt,'p'); 
	  magadvance(f,B,nmax,dt/3.0,gridx,gridy,gridz);
	  collision(f,nmax,dt/3.0,gridx,gridy,gridz,nu); 
	
	  t+=dt;     
	  count++;	       
      	
	  if((count%ign)==0)
	  {
	    i = count/ign;
	    //	vtkfile(f,B,gridx,gridy,gridz,dx,dy,dz,i);
	    cout << i << endl; 
	    cout <<  real(f[0][0][0][0][0]) <<  endl;
	    cout << "TIMESTEP" << "  " << count << "  " << t << endl;

	    //Data output.....average value of each m=0 harmonic as a test to see how many harmonics are required 
	     
	    magan = fftanalysis(f,Gridsize,t);
	    for(m=0;m<nmax;m++)
	      for(p=0;p<=m;p++)
		  {
		    harmamp=0.0;
		    for(n=0;n<gridx;n++) 
		      for(j=0;j<gridy;j++)
		        for(k=0;k<gridz;k++) 
			      harmamp+= abs(f[n][j][k][m][p])/(gridx*gridy*gridz);
		          harmonicfile << harmamp << '\t';
		  }
	      harmonicfile << endl;
	      if((count%(500*ign))==0)
	      {
		    thetathirecon(f,nmax,0,0,0,count);sprintf(filename,"TDdens%i.dat",count);
		    ofstream ddfile;
		    ddfile.open(filename);
		    ddfile.precision(16);
		    for(n=0;n<gridx;n++) 
		      for(j=0;j<gridy;j++)
		        for(k=0;k<gridz;k++) 
		          ddfile << (n+0.5)*dx << '\t' << (j+0.5)*dy << '\t' << (k+0.5)*dx << '\t' << real(f[n][j][k][0][0]) << '\t' << real(f[n][j][k][1][0]) << '\t' << real(f[n][j][k][1][1]) << '\t' << -imag(f[n][j][k][1][1]) << endl;
	       }	
	     }
     }while(magan>(0.01*maganini));
    //}while(t<=1000.0);
    
    FreeThreevector3D(B,gridx,gridy);
    FreeDistrF(f,gridx,gridy,gridy,nmax);
    delete [] str;
    return 0;
    
}   





