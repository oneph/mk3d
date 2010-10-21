#include <complex>
#include <iostream>
#include <fstream>

using namespace std;

void advectionz(complex<double>*****f,int nmax,int grid,int gridx,int gridy,double dz,double v,double dt,char bcon)//A routine that advects one time step forward IN TOKAMAK GEOMETRY X=Z
{
  
  //Advection routine in z using a second order Runge-Kutta method  
  double A,B;
  int pos;
  int n,m,i,j;
  complex<double> **grad,**fint,**fave;
  grad = new complex<double>*[grid+1];
  fint = new complex<double>*[grid+1];
  fave = new complex<double>*[grid+1];
  for(i=0;i<=grid;i++)
    {
      grad[i] = new complex<double>[nmax+1];
      fint[i] = new complex<double>[nmax+1];
      fave[i] = new complex<double>[nmax+1];
    }

  for(i=0;i<gridx;i++)
    for(j=0;j<gridy;j++)
      {
	for(m=0;m<nmax;m++)
	  {
	    for(pos=0;pos<=grid;pos++)
	      for(n=0;n<=nmax;n++)
		{
		  fave[pos][n] = (0.0,0.0);
		  fint[pos][n] = (0.0,0.0); 	 
		  grad[pos][n] = (0.0,0.0);
		}
	    
	    for(n=m;n<nmax;n++)
	      {
		for(pos=1;pos<grid;pos++)//Average over two spatial grid points
		  {
		    fave[pos][n] = 0.5*(f[i][j][pos][n][m] + f[i][j][pos-1][n][m]);
		  }	
		//BOUNDARY CONDITIONS!
		switch(bcon)
		  {
		  case 'p'://Periodic boundary
		    fave[0][n] = 0.5*(f[i][j][0][n][m] + f[i][j][grid-1][n][m]);
		    fave[grid][n] = fave[0][n];
		    break;
		  case 'r'://reflective boundary	 
		    if(((m+n)%2)==1)
		      {
			fave[0][n]=0.0;
			fave[grid][n]=0.0;
		      }else
			{
			  fave[0][n]=f[i][j][0][n][m];
			  fave[grid][n]=f[i][j][grid-1][n][m];
			}
		    break;
		  default:
		    cout << "No boundary condition defined in advectionz()" << endl;
		    terminate();
		  }
		
		//END OF BC
		for(pos=1;pos<grid;pos++)//Calcualte the gradient across each spatial element
		  {
		    grad[pos][n] = (f[i][j][pos][n][m] - f[i][j][pos-1][n][m])/dz;
		  }
		
		//Boundary conditions
		switch(bcon)
		  {
		  case 'p'://Periodic boundary
		    grad[0][n] = (f[i][j][0][n][m] - f[i][j][grid-1][n][m])/dz;
		    grad[grid][n] = grad[0][n];
		    break;
		  case 'r'://reflective boundary	 
		    if(((m+n)%2)==1)
		      {
			  grad[0][n] = 2.0*f[i][j][0][n][m]/dz;
			  grad[grid][n] = -2.0*f[i][j][grid-1][n][m]/dz;
		      }else
			{
			  grad[0][n] = 0.0;
			  grad[grid][n] = 0.0;
			}
		    break;
		  default:
		    cout << "No boundary condition defined in advectionz()" << endl;
		    terminate();
		  }
		
		//End of BC's
	      }
	    
	    for(n=m;n<nmax;n++)
	      {
		
		A = 1.0*(n-m)/(2*n-1);
		B = 1.0*(n+m+1)/(2*n+3);
		
		for(pos=0;pos<=grid;pos++)
		  {
		    if(n==0)
		      {
			fint[pos][n] = fave[pos][n] - 0.5*v*dt*(B*grad[pos][n+1]);
		      }else
			{
			  fint[pos][n] = fave[pos][n] - 0.5*v*dt*(A*grad[pos][n-1] + B*grad[pos][n+1]);
			}
		  }
	      } 
	    
	    
	    //The first half of the Runge Kutta algorithm complete
	    for(n=m;n<nmax;n++)
	      {	  	                 
		for(pos=0;pos<grid;pos++)//Calcualte the gradient across each spatial element
		  {
		    grad[pos][n] = (fint[pos+1][n] - fint[pos][n])/dz;
		  }
		
	      }
	    
	    for(n=m;n<nmax;n++)
	      {
		A = 1.0*(n-m)/(2*n-1);
		B = 1.0*(n+m+1)/(2*n+3);
		for(pos=0;pos<grid;pos++)
		  {
		    if(n==0)
		      {
			f[i][j][pos][n][m] += real(-v*dt*(B*grad[pos][n+1]));
		      }else
			{
			  f[i][j][pos][n][m]+= real(-v*dt*(A*grad[pos][n-1] + B*grad[pos][n+1]));
			}
		  }
	      }
	    
	  }
      }
  
  for(i=0;i<=grid;i++)
    {
      delete [] fint[i];
      delete [] fave[i];
      delete [] grad[i];
    }
  delete[] fint;
  delete[] fave;
  delete[] grad;
}




void advectionx(complex<double> *****f,int nmax,int grid,int gridy,int gridz,double dx,double v,double dt,char bcon)
{
  //Advection routine in y using a second order Runge-Kutta method
  complex<double> ***fint,***fave,***grad;
  complex<double> CompNul(0.0,0.0);
  fint = new complex<double>**[grid+1];
  fave = new complex<double>**[grid+1];
  grad = new complex<double>**[grid+1];
  for(int i=0;i<=grid;i++)
    {
      fint[i] = new complex<double>*[nmax+1];
      fave[i] = new complex<double>*[nmax+1];
      grad[i] = new complex<double>*[nmax+1];
      for(int j=0;j<=nmax;j++)
	{
	  fint[i][j] = new complex<double>[nmax+1];
	  fave[i][j] = new complex<double>[nmax+1];
	  grad[i][j] = new complex<double>[nmax+1];
	}
    }   
  double c1,c2,c3,c4;
  int pos;
  int n,m,j,k;
  for(j=0;j<gridy;j++)
    for(k=0;k<gridz;k++)
      {
  
  
	  for(m=0;m<=nmax;m++)
	    for(pos=0;pos<=grid;pos++)
	      for(n=0;n<=nmax;n++)
	      {
		fave[pos][n][m] = CompNul;
		fint[pos][n][m] = CompNul;
		grad[pos][n][m] = CompNul;
	      }
	
	
	for(m=0;m<nmax;m++)
	  for(n=m;n<nmax;n++)
	    {
	      for(pos=1;pos<grid;pos++)//Average over two spatial grid points
		{
		  fave[pos][n][m] = 0.5*(f[pos][j][k][n][m] + f[pos-1][j][k][n][m]);
		}	
	      
	      // 	//BOUNDARY CONDITIONS!
	      switch(bcon)
		{
		case 'p'://Periodic boundary
		  fave[0][n][m] = 0.5*(f[0][j][k][n][m] + f[grid-1][j][k][n][m]);
		  fave[grid][n][m] = fave[0][n][m];
		  break;
		case 'r'://reflective boundary
		  if((m%2)==1)
		    {
		      fave[0][n][m]=0.0;
		      fave[grid][n][m]=0.0;
		    }else
		      {
			fave[0][n][m]=conj(f[0][j][k][n][m]);
			fave[grid][n][m]=conj(f[grid-1][j][k][n][m]);
		      }
		  break;
		default:
		  cout << "No boundary condition defined in advectiony()" << endl;
		  terminate();
		}
	      
	      //END OF BC
	      for(pos=1;pos<grid;pos++)//Calcualte the gradient across each spatial element
		{
		  grad[pos][n][m] = (f[pos][j][k][n][m] - f[pos-1][j][k][n][m])/dx;
		}
	      
	      //Boundary conditions
	      switch(bcon)
		{
		case 'p'://Periodic boundary
		  grad[0][n][m] = (f[0][j][k][n][m] - f[grid-1][j][k][n][m])/dx;
		  grad[grid][n][m] = grad[0][n][m];
		  break;
		case 'r'://reflective boundary
		  if((m%2)==1)
		    {
		      grad[0][n][m] = 2.0*f[0][j][k][n][m]/dx;
		      grad[grid][n][m] = -2.0*f[grid-1][j][k][n][m]/dx;
		    }else
		      {
			grad[0][n][m] = 0.0;
			grad[grid][n][m] = 0.0;
		      }
		  break;
		default:
		  cout << "No boundary condition defined in advectiony()" << endl;
		  terminate();
		}
	      
	      //End of BC's
	    }
	
	for(m=0;m<nmax;m++)
	  for(n=m;n<nmax;n++)
	    {
	      c1= 1.0/(2*n-1);
	      c2= 1.0*(n-m)*(n-m-1)/(2*n-1);
	      c3= 1.0/(2*n+3);
	      c4= 1.0*(n+m+1)*(n+m+2)/(2*n+3);
	      for(pos=0;pos<=grid;pos++)
		{
		  if(m==0)
		    {
		      if(n==0)
			{
			  fint[pos][n][0] = fave[pos][n][0] - 0.5*dt*v*real(c4*grad[pos][n+1][1]);
			}else
			  {
			    fint[pos][n][0] = fave[pos][n][0] + 0.5*dt*v*real((c2*grad[pos][n-1][1]) - (c4*grad[pos][n+1][1]));
			  }
		    }else
		      {
			fint[pos][n][m] = fave[pos][n][m] - 0.25*v*dt*((c1*grad[pos][n-1][m-1]) - (c2*grad[pos][n-1][m+1])) + 0.25*v*dt*((c3*grad[pos][n+1][m-1]) - (c4*grad[pos][n+1][m+1]));
		      }
		}
	    }
	
	// The first half of the Runge Kutta algorithm complete
	for(m=0;m<nmax;m++)   
	  for(n=m;n<nmax;n++) 	                 
	    for(pos=0;pos<grid;pos++)//Calcualte the gradient across each spatial element
	      {
		grad[pos][n][m] = (fint[pos+1][n][m] - fint[pos][n][m])/dx;
	      }
	
	for(m=0;m<nmax;m++)    
	  for(n=m;n<nmax;n++)
	    {
	      c1= 1.0/(2*n-1);
	      c2= 1.0*(n-m)*(n-m-1)/(2*n-1);
	      c3= 1.0/(2*n+3);
	      c4= 1.0*(n+m+1)*(n+m+2)/(2*n+3);
	      for(pos=0;pos<grid;pos++)
		{
		  if(m==0)
		    {
		      if(n==0)
			{
			  f[pos][j][k][n][0] += -dt*v*real(c4*grad[pos][n+1][1]);
			}else
			  {
			    f[pos][j][k][n][0] += dt*v*real((c2*grad[pos][n-1][1]) - (c4*grad[pos][n+1][1]));
			  }
		    }else
		      {  
			f[pos][j][k][n][m] += -0.5*v*dt*((c1*grad[pos][n-1][m-1]) - (c2*grad[pos][n-1][m+1])) + 0.5*v*dt*((c3*grad[pos][n+1][m-1]) - (c4*grad[pos][n+1][m+1]));
		      }
		}
	    }
      }   

  for(j=0;j<=grid;j++)
  {
      for(k=0;k<=nmax;k++)
      {
	  delete[] fint[j][k];
	  delete[] fave[j][k];
	  delete[] grad[j][k];
      }
      delete [] fint[j];
      delete [] fave[j];
      delete [] grad[j];
  }
  delete [] fint;
  delete [] fave;
  delete [] grad;
}  
      

void advectiony(complex<double> *****f,int nmax,int grid,int gridx,int gridz,double dy,double v,double dt,char bcon)
{
//Advection routine in y using a second order Runge-Kutta method
 complex<double> tempcomp;
 complex<double> i(0.0,1.0);
 complex<double> CompNul(0.0,0.0);
 double c1,c2,c3,c4;
 int pos;
 int n,m,ii,k;
 
  complex<double> ***fint,***fave,***grad;
  fint = new complex<double>**[grid+1];
  fave = new complex<double>**[grid+1];
  grad = new complex<double>**[grid+1];
  for(ii=0;ii<=grid;ii++)
    {
      fint[ii] = new complex<double>*[nmax+1];
      fave[ii] = new complex<double>*[nmax+1];
      grad[ii] = new complex<double>*[nmax+1];
      for(int j=0;j<=nmax;j++)
	{
	  fint[ii][j] = new complex<double>[nmax+1];
	  fave[ii][j] = new complex<double>[nmax+1];
	  grad[ii][j] = new complex<double>[nmax+1];
	}
    }   
  
  for(ii=0;ii<gridx;ii++)
    for(k=0;k<gridz;k++)
      {
	
	for(m=0;m<=nmax;m++)
	  for(pos=0;pos<=grid;pos++)
	    for(n=0;n<=nmax;n++)
	      {
		fave[pos][n][m] = CompNul;
		fint[pos][n][m] = CompNul;
		grad[pos][n][m] = CompNul;
	      } 
	
	for(m=0;m<nmax;m++)
	  for(n=m;n<nmax;n++)
	    {
	      for(pos=1;pos<grid;pos++)//Average over two spatial grid points
		{
 		  fave[pos][n][m] = 0.5*(f[ii][pos][k][n][m] + f[ii][pos-1][k][n][m]);
 		}
	      
	      
	      //BOUNDARY CONDITIONS!
	      switch(bcon)
		{
		case 'p'://Periodic boundary
		  fave[0][n][m] = 0.5*(f[ii][0][k][n][m] + f[ii][grid-1][k][n][m]);
		  fave[grid][n][m] = fave[0][n][m];
		  break;
		case 'r'://reflective boundary	
		  fave[0][n][m]=real(f[ii][0][k][n][m]);
		  fave[grid][n][m]=real(f[ii][grid-1][k][n][m]);
		  break;
		default:
		  cout << "No boundary condition defined in advectiony()" << endl;
		  terminate();
		}
	      
	      //END OF BC
	      for(pos=1;pos<grid;pos++)//Calcualte the gradient across each spatial element
		{
		  grad[pos][n][m] = (f[ii][pos][k][n][m] - f[ii][pos-1][k][n][m])/dy;
		}
	      
	      //Boundary conditions
	      switch(bcon)
		{
		case 'p'://Periodic boundary	    
		  grad[0][n][m] = (f[ii][0][k][n][m] - f[ii][grid-1][k][n][m])/dy;
		  grad[grid][n][m] = grad[0][n][m];
		  break;
		case 'r'://reflective boundary
		  grad[0][n][m] = 2.0*i*imag(f[ii][0][k][n][m])/dy;
		  grad[grid][n][m] =  -2.0*i*imag(f[ii][grid-1][k][n][m])/dy;
		  break;
		default:
		  cout << "No boundary condition defined in advectiony()" << endl;
		  terminate();
		}
	      
	      //End of BC's
	    }
	
	for(m=0;m<nmax;m++)
	  for(n=m;n<nmax;n++)
	    {
	      c1= 1.0/(2*n-1);
	      c2= 1.0*(n-m)*(n-m-1)/(2*n-1);
	      c3= 1.0/(2*n+3);
	      c4= 1.0*(n+m+1)*(n+m+2)/(2*n+3);
	      
	      for(pos=0;pos<=grid;pos++)
		{
		  if(m==0)
		    {
		      if(n==0)
			{
			  fint[pos][n][m] = fave[pos][n][m] - real(0.5*i*dt*v*(c4*grad[pos][n+1][1]));			
			}else
			  {
			    fint[pos][n][m] = fave[pos][n][m] + real(0.5*i*dt*v*(c2*grad[pos][n-1][1] - c4*grad[pos][n+1][1]));			
			  }
		    }else
		      {
			fint[pos][n][m] = fave[pos][n][m] + 0.25*i*dt*v*(c1*grad[pos][n-1][m-1] + c2*grad[pos][n-1][m+1] - c3*grad[pos][n+1][m-1] - c4*grad[pos][n+1][m+1]);
		      }
		}
	    }
	
	// The first half of the Runge Kutta algorithm complete
	for(m=0;m<nmax;m++)   
	  for(n=m;n<nmax;n++) 	                 
	    for(pos=0;pos<grid;pos++)//Calcualte the gradient across each spatial element
	      {
		grad[pos][n][m] = (fint[pos+1][n][m] - fint[pos][n][m])/dy;
	      }
	
	for(m=0;m<nmax;m++)    
	  for(n=m;n<nmax;n++)
	    {
	      c1= 1.0/(2*n-1);
	      c2= 1.0*(n-m)*(n-m-1)/(2*n-1);
	      c3= 1.0/(2*n+3);
	      c4= 1.0*(n+m+1)*(n+m+2)/(2*n+3);
	      for(pos=0;pos<grid;pos++)
		{
		  if(m==0)
		    {
		      if(n==0)
			{		    
			  f[ii][pos][k][n][m] +=  -real(i*dt*v*(c4*grad[pos][n+1][1]));			
			}else
			  {		      
			    f[ii][pos][k][n][m] += real(i*dt*v*(c2*grad[pos][n-1][1] - c4*grad[pos][n+1][1]));			
			  }
		    }else
		      {  
			f[ii][pos][k][n][m] += 0.5*i*dt*v*(c1*grad[pos][n-1][m-1] + c2*grad[pos][n-1][m+1] - c3*grad[pos][n+1][m-1] - c4*grad[pos][n+1][m+1]);
		      }
		}
	    }
      }
  
  for(n=0;n<=grid;n++)
  {
      for(m=0;m<=nmax;m++)
      {
	  delete[] fint[n][m];
	  delete[] fave[n][m];
	  delete[] grad[n][m];
      }
      delete [] fint[n];
      delete [] fave[n];
      delete [] grad[n];
  }
  delete [] fint;
  delete [] fave;
  delete [] grad;
}
