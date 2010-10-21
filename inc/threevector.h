//A class to produce a three vector
#include <iostream>

class threevector
{
 private:
  double xcoord,ycoord,zcoord;
  
 public:
  threevector()
      {
	  double xcoord = 0.0;
	  double ycoord = 0.0;
	  double zcoord = 0.0;
      }

  threevector(double x,double y,double z)
      {
	  xcoord = x;
	  ycoord = y;
	  zcoord = z;
      }
  
  void print()
      {
	  std::cout << xcoord << '\t' << ycoord << '\t' << zcoord << std::endl;
      }
  
  double getx()
      {
	  return xcoord;
      }
  
  double gety()
      {
	  return ycoord;
      }
  
  double getz()
      {
	  return zcoord;
      }
  
  double mag()
      {
	  double a = xcoord*xcoord + ycoord*ycoord + zcoord*zcoord;
	  return sqrt(a);
      }
  
  void setx(double x)
      {
	  xcoord = x;
      }
  
  void sety(double y)
      {
	  ycoord = y;
      }
  
  void setz(double z)
      {
	  zcoord = z;
      }
  
  threevector operator*(double num)
      {
	  threevector temp(xcoord*num,ycoord*num,zcoord*num);
	  return temp;
      }

  threevector operator+(threevector c)
      {
	  threevector temp(xcoord + c.getx(),ycoord + c.gety(),zcoord + c.getz());
	  return temp;
    }
};
