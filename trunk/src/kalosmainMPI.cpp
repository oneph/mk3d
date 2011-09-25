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
#include "MPIfunctions.h"
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{   
    int MPIrank,MPIsize,ierr;
    //Initialise MPI
    ierr =  MPI_Init(&argc, &argv);
    if((ierr!=MPI_SUCCESS))
    {
       cout << "Could not initialise MPI" << endl;
       return 0;
    }

    ierr = MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);    
    cout << MPIrank << " of " << MPIsize << " processes" << endl; 

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
    int i,j,k,n,m,p;
    char filename[20];
    cout.precision(8);
    compfile.precision(16);
    
    //ANY PHYSICAL CONSTANTS ARE DEFINED HERE!!!!
    const double pi = acos(-1.0);  
    const complex<double> ii(0.0,1.0);
    //END OF PHYSICAL CONSTANTS!
    
	//Number of harmonics (global)
	int nmaxG;
	//Number of harmonics (local)
	int nmax;
	//The maximum extent of the spatial dimension globally
    double LxG,LyG,LzG;
	//The local extent of the grid
	double Lx,Ly,Lz;
    //Local grid size
	int gridx,gridy,gridz,ign;
	//Global grid sizes
	int gridxG,gridyG,gridzG;
	//Number of processors in each direction
	int N_procs_x,N_procs_y,N_procs_z;
	
    double t=0;//The time coordinate
    double dx,dy,dz,dt;//Time and spatial increments
    double dxG,dyG,dzG,dtG;
  
  	//Parameters for MPI decomposition
    int mpindim = 1;
    int nprocsdim[1];
    int periods[1];
    int reorder;
    MPI_Comm comm1d;
    //---------------------------------
    
    double Bmag,dB,Bmax=0.0;
    int count = 0;//Timestep count - Integer number useful for file naming
    double v;//Velocity magnitude
    double nu;//The collisional frequency
    threevector Masscentre;
    complex<double> sum;
    double totdens;
    time_t curr; 
    char *str;
	//Following code loads parameters from the file (inputdeck.txt)created in external script 
    
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
      if(MPIrank==0)
      {
        cout << str << " " << num1  << endl;
      }
      gridxG = num1;
      finput >> str >> num1;
      if(MPIrank==0)
      {
        cout << str << " " << num1  << endl;
      }
      gridyG = num1;
      finput >> str >> num1;
      if(MPIrank==0)
      {
        cout << str << " " << num1  << endl;
      }
      gridzG = num1;  
    
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      LxG = num;
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      LyG = num;
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      LzG = num;
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      nu = num;
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      v = num;  
      finput >> str >> num1;
      if(MPIrank==0)
      {
        cout << str << " " << num1  << endl;
      }
      nmaxG = num1;//The number of harmonics in the expansion !!!!Expansion goes from 0 to nmax-1!!!!
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      Bmag = num;
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      dB = num;  
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
      ign = int(num);
	
	  //Processor layout
      finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
	  }
	  N_procs_x = int(num);
	  finput >> str >> num;
      if(MPIrank==0)
      {
        cout << str << " " << num << endl;
      }
	  N_procs_y = int(num);  
	  finput >> str >> num;
	  if(MPIrank==0)
	  {
	    cout << str << " " << num << endl;
	  }
	  N_procs_z = int(num);
		
      curr = time(NULL);
	  if(MPIrank==0)
	  {
	    cout << "The date is = " << ctime(&curr);
	  }
	  
	  //Function called which calculates the optimal parallelisation layout
      if((N_procs_x*N_procs_y*N_procs_z)!=MPIsize)
	  {
		  if(MPIrank==0)
		  {
		    cout << "Mismatch between processor numbers and MPI size" << endl;
		  }
		  ierr = MPI_Finalize();
		  return 0;
	  }
		
      //For the time being local values are set to the global ones
      //apart from the z direction
	  gridx = gridxG/N_procs_x;
	  gridy = gridyG/N_procs_y;
	  gridz = gridzG/N_procs_z;
	  Lx = LxG/N_procs_x;
	  Ly = LyG/N_procs_y;
	  Lz = LzG/N_procs_z;
	  nmax = nmaxG;
	  if(gridx<1 || gridy<1 || gridz<1)
	  {
		cout << "Too many processor in one (or more) directions" << endl;
		cout << gridx << " " << gridy << " " << gridz << endl;
		ierr = MPI_Finalize();
		return 0;
	  }
	
	if(MPIrank==0)
	{
	  cout << "Points x y z on each processor" << endl;
      cout << gridx << " " << gridy << " " << gridz << endl;
	}
	//File load ends here
	
	nprocsdim[0] = 4;
	mpindim = 1;
	periods[0] = 1;
	reorder = 1;

    //Setting up MPI Grid.
    ierr = MPI_Cart_create(MPI_COMM_WORLD,mpindim,nprocsdim,periods,reorder,&comm1d);
     
    //Integer determining which direction 
    int dir = 0;
    //Integer determining the shift
    int shift = 1;
    int pleft, pright;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&MPIrank,&pright);
    cout << "Rank" << MPIrank << " " << pright << endl;
    shift = -1;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&MPIrank,&pleft);    
    cout << "Rank" << MPIrank << " " << pleft << endl;
    //----------------------------------------------------------------------------

    threevector Systemsize(Lx,Ly,Lz);
    threevector SystemsizeGlobal(LxG,LyG,LzG);
    threevector Gridsize(gridx,gridy,gridz);
    threevector GridsizeGlobal(gridxG,gridyG,gridzG);
    
    //Distribution function - First array is the x position element,second is y, third is z, fourth array is the n Legengre ploynomial number, and m is the harmonic
    complex<double> *****f;
    //complex<double> ****buffer;
    //buffer = new complex<double> *[gridx];
    //for(i=0
    
    threevector ***B;//Magnetic field vector
    threevector ***A;//Magnetic thingy
    B = new threevector**[gridx];
    A = new threevector**[gridx+1];
    for(i=0;i<=gridx;i++)  
    {		
	A[i] = new threevector*[gridy+1];
	for(j=0;j<=gridy;j++)
	{
	    A[i][j] = new threevector[gridz+1]; 
	}
    }
    for(i=0;i<gridx;i++)  
    {		
	B[i] = new threevector*[gridy];
	for(j=0;j<gridy;j++)
	{
	    B[i][j] = new threevector[gridz]; 
	}
    }
    
    //End pf definition of distro function - Multi dimensional arrays = PITA
    //Cartesian Grid!!!
    dx = Lx/gridx;
    dy = Ly/gridy;
    dz = Lz/gridz;
    
    threevector Elementsize(dx,dy,dz);
    double scale = smallest(dx,smallest(dy,dz));
    if(MPIrank==0)
    {
      cout << "SMALLEST LENGTH-SCALE = " << scale << endl;
    }
    dt = scale/(4.0*v);
    ierr = MPI_Allreduce(&dt,&dtG,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    if(MPIrank==0)
    {
      cout << "TIMESTEP SIZE =" << dtG << endl;    
    }
    
    //!!!!!! INITIAL CONDTIIONS!!!!!!!! 

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

    for(i=0;i<=gridx;i++)//Deallocation of the memory allocated to the vector potential.
    {
	  for(j=0;j<=gridy;j++)
	  {
	    delete [] A[i][j];
	  }
	  delete [] A[i];
    }
    delete [] A;
    
    //END OF FIELD INITIALISATION.  ONLY THE B-FIELD DATA REMAINS.....VECTOR POTENTIAL IS DELETED

    f = new complex<double>****[gridx];
    for(i=0;i<gridx;i++)
    {		
	  f[i] = new complex<double>***[gridy];
   	  for(j=0;j<gridy;j++)
	  {
	    f[i][j] = new complex<double>**[gridz];
	    for(k=0;k<gridz;k++)
	    {
		  f[i][j][k]  = new complex<double>*[nmax];
		  for(m=0;m<nmax;m++)
		  {
		    f[i][j][k][m] = new complex<double>[nmax];
		  }
	    }
 	  }
    }
 
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
	    
//     cout << "Maximum magnetic field strength is " << Bmax << endl;
//     if((2*pi/(4*Bmax*nmax))<dt)
//     {
// 	dt = 2*pi/(Bmax*nmax);
// 	cout << "Timestep redefined due to high field strength  dt=" << dt << endl;
//     }else
//     {
// 	cout << "Timestep defined by CFL  = " << dt << endl;
//     }

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

    /*
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
    */
    
    
    for(i=0;i<gridx;i++)
      {
	    for(j=0;j<gridy;j++)
	      {
	        delete [] B[i][j];
	      }
	    delete [] B[i];
      }
    delete [] B;

    for(i=0;i<gridx;i++) //Deallocation of the memory used for the distro function    
      {
	    for(j=0;j<gridy;j++)
 	      {
	        for(k=0;k<gridz;k++)
	          {
		        for(n=0;n<nmax;n++)	 
		          {
		            delete [] f[i][j][k][n];
		          }
		        delete [] f[i][j][k];
	          }
	        delete [] f[i][j];
	      }
	    delete [] f[i];
      }
    delete [] f;
    delete [] str;
	
    ierr = MPI_Finalize();   
    return 0;
    
}   





