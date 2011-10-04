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
    ierr =  MPI_Init(&argc, &argv);
    if((ierr!=MPI_SUCCESS))
    {
       cout << "Could not initialise MPI" << endl;
       return 0;
    }
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);    
   
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
	int parX = 0,parZ = 0,parY = 0;
	
    double t=0;//The time coordinate
    double dx,dy,dz,dt;//Time and spatial increments
    double dxG,dyG,dzG,dtG;
  
  	//Parameters for MPI decomposition-At the moment set up for 1dimensional decomposition
    int mpindim;
    int nprocsdim[3];
    int periods[3];
    int reorder;
    int pleft, pright, me;
    MPI_Comm comm1d,comm2d,comm3d;
    //---------------------------------
    
    double Bmag,dB,Bmax=0.0;
    int count = 0;//Timestep count - Integer number useful for file naming
    double v;//Velocity magnitude
    double nu;//The collisional frequency
    threevector Masscentre;
    complex<double> sum;
    double totdens;
    time_t curr; 
    char *str,str2;
    char Bc_x,Bc_y,Bc_z, rin;
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
	
	//Boundary conditions
    finput >> str >> rin;
    if(MPIrank==0)
    {
      cout << str << " " << rin << endl;
	}
	Bc_x = rin;
	finput >> str >> rin;
    if(MPIrank==0)
    {
      cout << str << " " << rin << endl;
    }
	Bc_y = rin;  
	finput >> str >> rin;
	if(MPIrank==0)
	{
	  cout << str << " " << rin << endl;
	}
	Bc_z = rin;
		
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
	nprocsdim[0] = N_procs_z;
	nprocsdim[1] = N_procs_x;
	nprocsdim[2] = N_procs_y;

	nmax = nmaxG;//For the time being no parallelisation over the harmonics.
	
	if(gridx<1 || gridy<1 || gridz<1)
	{
	  cout << "Too many processors in one (or more) directions" << endl;
	  cout << gridx << " " << gridy << " " << gridz << endl;
	  ierr = MPI_Finalize();
  	  return 0;
	}
	
	if(MPIrank==0)
	{
	  cout << "Points x y z on each processor" << endl;
      cout << gridx << " " << gridy << " " << gridz << endl;
      cout << "Bc_x" << " " << "Bc_y" << " " << "Bc_z" << endl;
	  cout << Bc_x << " " << Bc_y << " " << Bc_z << endl;
	}
	//File load ends here
	
	//Code finds out how many directions to parallelise over.
	mpindim = 0;
	if(N_procs_x>1)
	{
	  mpindim += 1;
	  parX = 1;
	}
	if(N_procs_y>1)
	{
	  mpindim += 1;
	  parY = 1;
	}
	if(N_procs_z>1)
	{
	  mpindim += 1;
	  parZ = 1;
	}
	
	if(MPIrank==0)
	{
	  cout << "The code will parallelise over " << mpindim << " directions." << endl;
	}
	if(mpindim==0 && MPIrank==0)
	{
	  cout << "You are trying to run a single procesor job with the multi-processor code" << endl;
	}
	
	//--------------------------------------------------------
	
	//Set parallelisation boundary conditions for communicators.
	if(Bc_z=='p')
	{
	  periods[0] = 1;
	}else
		{
		  periods[0] = 0;
		}
		
	if(Bc_x=='p')
	{
	  periods[1] = 1;
	}else
		{
		  periods[1] = 0;
		}
		
	if(Bc_y=='p')
	{
	  periods[2] = 1;
	}else
		{
		  periods[2] = 0;
		}
	reorder = 1;
    //--------------------------------------------------------

    //Setting up MPI Grid and communicators.
    ierr = MPI_Cart_create(MPI_COMM_WORLD,mpindim,nprocsdim,periods,reorder,&comm1d);
    ierr = MPI_Comm_rank(comm1d, &me); 
    //ierr = MPI_Cart_coords(comm1d, me, mpindim, &coords);
    //get_neighbours1D(comm1d,me,pright,pleft);
      
    //Integer determining which direction -> In this case the only direction
    int dir = 0;
    //Integer determining the shift
    int shift = 1;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&me,&pright);
    shift = -1;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&me,&pleft);  
      
    cout << "RankL" << " " << MPIrank << " " << pleft << " " << "RightR" << " " << pright << endl;
    //----------------------------------------------------------------------------

    threevector Systemsize(Lx,Ly,Lz);
    threevector SystemsizeGlobal(LxG,LyG,LzG);
    threevector Gridsize(gridx,gridy,gridz);
    threevector GridsizeGlobal(gridxG,gridyG,gridzG);
    
    //Distribution function - First array is the x position element,second is y, third is z, fourth array is the n Legengre ploynomial number, and m is the harmonic
    complex<double> *****f;
    complex<double> ****buffer;
    
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
	    Axfile << A[i][k][gridz/2].getx() << '\t';
	    Ayfile << A[i][k][gridz/2].gety() << '\t'; 
	    Azfile << A[i][k][gridz/2].getz() << '\t';
	  }
	  Bxfile << endl;	
	  Byfile << endl;
	  Bzfile << endl;
	  Bmagfile << endl;
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
    
    //Buffer for inter-processor communication-Set up for one way decomposition
    buffer = new complex<double>***[gridx];
    for(i=0;i<gridx;i++)
    {		
	  buffer[i] = new complex<double>**[gridy];
   	  for(j=0;j<gridy;j++)
	  {
	    buffer[i][j] = new complex<double>*[nmax];
		for(m=0;m<nmax;m++)
		{
		  buffer[i][j][m] = new complex<double>[nmax];
		}
 	  }
    }
 
    //!!!!!! INITIAL CONDTIIONS!!!!!!!!
    for(i=0;i<gridx;i++) 
	  for(j=0;j<gridy;j++)
	    for(k=0;k<gridz;k++) 
	    {
	      for(m=0;m<nmax;m++)  
		    for(n=m;n<nmax;n++)
		    {
			  f[i][j][k][n][m] = 0.0;
	 	    }
	        f[i][j][k][0][0] = 1.0 + 0.01*cos(1.0*pi*(((j+0.5)*dy/Ly)+((i+0.5)*dx/Lx)))*cos(1.0*pi*(((j+0.5)*dy/Ly)-((i+0.5)*dx/Lx)));
	        compfile << (i+0.5)*dx << '\t' << (j+0.5)*dy << '\t' << (k+0.5)*dz << '\t' << B[i][j][k].getx() << '\t' << B[i][j][k].gety()  << '\t' << B[i][j][k].getz() << endl;
	        if(B[i][j][k].mag()>Bmax)
		    {
		      Bmax = B[i][j][k].mag();
		    }
	    }

    //End pf definition of distro function - Multi dimensional arrays = PITA
    //END OF INITIAL CONDITIONS
    
    //Code that writes to the Run Parameters file
    if(MPIrank==0)
	{
      datfile << "Run Parameters for KALOS run on" << " " << ctime(&curr)  << endl;
      datfile << "Length of system  , Lx = " << LxG << ", Ly = " << LyG << ", Lz = " << LzG << endl;
      datfile << "Boundary conditions, Bc_x = " << Bc_x << ", Bc_y = " << Bc_y << ", Bc_z = " << Bc_z << endl;
      datfile << "Number of Grid cells, Nx = " << gridxG << ", Ny = " << gridyG << ", Nz = " << gridzG << endl;
      datfile << "Number of Harmonics = " << nmax << endl;
      datfile << "Thermal velocity, vth =" << " " << v << endl;
      datfile << "Collision frequency = " << nu << endl;
      datfile << "Timestep =" << dt << endl;
      //End of parameter file.
    }
    
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
    
	if(parX==1)
	{
	  communicate_buffer_x;
	}
	
	if(parY==1)
	{
	  communicate_buffer_y;
	}
	
	if(parZ==1)
	{
	  communicate_buffer_z;
	}
	
	advectionx(f,nmax,gridx,gridy,gridz,dx,v,dt,Bc_x);
	magadvance(f,B,nmax,dt/3.0,gridx,gridy,gridz);
	collision(f,nmax,dt/3.0,gridx,gridy,gridz,nu);  
	
	advectiony(f,nmax,gridy,gridx,gridz,dy,v,dt,Bc_y);
	magadvance(f,B,nmax,dt/3.0,gridx,gridy,gridz);
	collision(f,nmax,dt/3.0,gridx,gridy,gridz,nu); 
	
	advectionz(f,nmax,gridz,gridx,gridy,dz,v,dt,Bc_y); 
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
    
    for(i=0;i<gridx;i++) //Deallocation of the memory used for the distro function    
      {
	    for(j=0;j<gridy;j++)
 	      {
		    for(n=0;n<nmax;n++)	 
		      {
		        delete [] buffer[i][j][n];
		      }
		      delete [] buffer[i][j];
	       }
	       delete [] buffer[i];
        }
    delete [] f;
    delete [] buffer;
    delete [] str;
	
    ierr = MPI_Finalize();   
    return 0;
    
}   





