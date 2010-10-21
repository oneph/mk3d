#include <complex>
#include <iostream>
#include <fstream>

using namespace std;

void vtkfile(complex<double>*****f,threevector ***B,int gridx,int gridy,int gridz,double dx,double dy,double dz,int t)//Reconstructs the distribution function cross the magnetic field....This uses the structured points system of data
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k;
  char filename[35];
  sprintf(filename,"3datat%i.vtk",t);
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Data for timestep" << " " << t << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS" << " " << gridx << " " << gridy << " " << gridz << endl;
  outfile << "ASPECT_RATIO" << " " << dx << " " << dy << " " << dz << endl;
  outfile << "ORIGIN" << " " << 0 << " " << 0 << " " << 0 << endl;
  outfile << "POINT_DATA" << " " << gridx*gridy*gridz << endl;
  outfile << "SCALARS" << " " << "Density" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  for(k=0;k<gridz;k++)
    {
      for(j=0;j<gridy;j++)
	{
	  for(i=0;i<gridx;i++)
	    outfile << real(f[i][j][k][0][0]) << " ";
	}
      outfile << endl;
    }
  outfile.close();

}

void threedreconvtk(complex<double>*****f,int gridx,int gridy,int gridz,double dx,double dy,double dz,int nmax,int t)
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k,n,m;
  const double pi = acos(-1.0);
  double theta,thi,r;//The spherical coordinates
  char filename[35];
  sprintf(filename,"Veldist%i.vtk",t);
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "A reconstruction of the distribution function in velocity space" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << 1 << " " << 51 << " " << 51 << endl;
  outfile << "POINTS" << " " << 2601 << " " << "float" << endl;
  for(i=0;i<=50;i++)//Loop over theta and thi
    for(j=0;j<=50;j++)
      {
	r=1.0;
	theta = 1.0*i*pi/50.0;
	thi = 2.0*j*pi/50.0;
	outfile << r*sin(theta)*cos(thi) << " " << r*sin(theta)*sin(thi) << " " << r*cos(theta) << endl;
      }
  outfile << "POINT_DATA" << " " << 2601 << endl;
  outfile << "SCALARS" << " " << "VelDist" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;

  for(int k=0;k<=50;k++)
    {
      for(int j=0;j<=50;j++)
	{
	  complex<double> sum(0.0,0.0);
	  for(n=0;n<nmax;n++)
	    {
	      sum += 1.0*plgndr(n,0,cos(k*pi/50))*(real(f[0][0][0][n][0]));
	      for(m=1;m<=n;m++)
		{
		  sum += 2.0*plgndr(n,m,cos(k*pi/50))*((real(f[0][0][0][n][m])*cos(m*2.0*j*pi/50)) - (imag(f[0][0][0][n][m])*sin(m*2.0*j*pi/50)));
		}
	    }
	  outfile << real(sum) << '\t';
	}
      outfile << endl;
    }
  outfile.close();

/*   for(i=0;i<=50;i++) */
/*     { */
/*       for(j=0;j<=50;j++) */
/* 	{ */
/* 	  outfile << real(angularreconst(f,nmax,gridx,gridy,gridz,1.0*pi*i/50.0,2.0*pi*j/50.0)) << " "; */
/* 	} */
/*       outfile << endl; */
/*     } */
/*   outfile.close();  */

}

/*void pxpzphaseplot(complex<double>*****f,int gridx,int gridy,int gridz,double dx,double dy,double dz,int nmax,int t)
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k;
  const double pi = acos(-1.0);
  double theta,thi,r;//The spherical coordinates
  char filename[35];
  sprintf(filename,"PxPzphase%i.vtk",t);
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "A reconstruction of the distribution function in velocity space" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << 3 << " " << 51 << " " << 1 << endl;
  outfile << "POINTS" << " " << 153 << " " << "float" << endl;
  for(i=0;i<51;i++)//Loop over theta and thi  
    for(j=0;j<3;j++)
      { 
	r=0.9+(0.1*j);
	theta = 2.0*i*pi/50.0;
	thi = pi/2.0;
	outfile << r*cos(theta) << " " << r*sin(theta)*sin(thi) << " " << 0.0 << endl;
      }
  outfile << "POINT_DATA" << " " << 153 << endl;
  outfile << "SCALARS" << " " << "VelDist" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  for(i=0;i<51;i++) 
    {
      for(j=0;j<3;j++)
	{
	  if(j==1)
	    outfile << real(angularreconst(f,nmax,gridx,gridy,gridz,2.0*pi*i/50.0,pi/2.0)) << " "; 
	  else
	    outfile << 0.0 << " ";
 	}
      outfile << endl;
    }
  outfile.close();

} 

void pxpyphaseplot(complex<double>*****f,int gridx,int gridy,int gridz,double dx,double dy,double dz,int nmax,int t)
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k;
  const double pi = acos(-1.0);
  double temp;
  double theta,thi,r;//The spherical coordinates
  char filename[35];
  sprintf(filename,"PxPyphase%i.vtk",t);
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "A reconstruction of the distribution function in velocity space" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << 3 << " " << 1 << " " << 101 << endl;
  outfile << "POINTS" << " " << 303 << " " << "float" << endl;
  for(i=0;i<101;i++)//Loop over theta and thi  
    for(j=0;j<3;j++)
      { 
	r=0.9+(0.1*j);
	theta = pi/2.0;
	thi = 2.0*i*pi/100.0;
	outfile << 0.0 << " " << r*sin(theta)*sin(thi) << " " << r*sin(theta)*cos(thi) << endl;
      }
  outfile << "POINT_DATA" << " " << 303 << endl;
  outfile << "SCALARS" << " " << "VelDist" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  for(i=0;i<101;i++) 
    {
      for(j=0;j<3;j++)
	{
	  if(j==1)
	    {
	      complex<double> sum0 = 0.0;
	      complex<double> summ = 0.0;
	      for(int n=0;n<nmax;n++)
		{
		  double temp = plgndr(n,0,cos(pi/2.0));
		  sum0 += f[gridx/2][gridy/2][gridz/2][n][0]*temp;
		  for(int m=1;m<=n;m++)
		    {
		      temp = plgndr(n,m,cos(pi/2.0));
		      summ += temp*((2.0*real(f[gridx/2][gridy/2][gridz/2][n][m])*cos(2.0*m*pi*i/100.0)) - (2.0*imag(f[gridx/2][gridy/2][gridz/2][n][m])*sin(2.0*m*pi*i/100.0)));
		    }
		}
	      outfile << real(sum0+summ) << " "; 
	    }
	  else
	    outfile << 0.0 << " ";
	}
      outfile << endl;
    }
  outfile.close();
  
  } */

void vtkfilemagstrength(threevector ***B,int gridx,int gridy,int gridz,double dx,double dy,double dz)//Reconstructs the distribution function cross the magnetic field....This uses the structured points system of data
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k;
  char filename[35];
  sprintf(filename,"Magstrength.vtk");
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Data for timestep" << " " << 0 << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS" << " " << gridx << " " << gridy << " " << gridz << endl;
  outfile << "ASPECT_RATIO" << " " << dx << " " << dy << " " << dz << endl;
  outfile << "ORIGIN" << " " << 0 << " " << 0 << " " << 0 << endl;
  outfile << "POINT_DATA" << " " << gridx*gridy*gridz << endl;
  outfile << "SCALARS" << " " << "Density" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  for(k=0;k<gridz;k++)
    {
      for(j=0;j<gridy;j++)
	{
	  for(i=0;i<gridx;i++)
	    outfile << B[i][j][k].mag() << " ";
	}
      outfile << endl;
    }
  outfile.close();

}

void vtkfilemagvect(threevector ***B,int gridx,int gridy,int gridz,double dx,double dy,double dz)//Reconstructs the distribution function cross the magnetic field....This uses the structured points system of data
{
  //Code for the production of datafiles that are incremently numbered, density, current and pressure components are logged as well as the structure of the magnetic field
  int i,j,k;
  char filename[35];
  sprintf(filename,"MagVectors.vtk");
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Data for timestep" << " " << 0 << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS" << " " << gridx << " " << gridy << " " << gridz << endl;
  outfile << "ASPECT_RATIO" << " " << dx << " " << dy << " " << dz << endl;
  outfile << "ORIGIN" << " " << 0 << " " << 0 << " " << 0 << endl;
  outfile << "POINT_DATA" << " " << gridx*gridy*gridz << endl;
  outfile << "VECTORS" << " " << "Magfields" << " " << "float" << endl;  
  for(i=0;i<gridx;i++)
    {
	for(j=0;j<gridy;j++)
	{
	    for(k=0;k<gridz;k++)
		outfile << B[i][j][k].getx() << '\t' <<  B[i][j][k].gety() << '\t' <<  B[i][j][k].getz() << endl;;
	}
      outfile << endl;
    }
  outfile.close();

}
