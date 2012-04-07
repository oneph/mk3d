//Multidimensional dynamical arrays templates are defined here. 
#include <complex>
#include <iostream>
#include <fstream>

using namespace std;

template <typename T> 
T **AllocateDynamicArray2D( int x, int y)
{
      T **dynamicArray;

      dynamicArray = new T*[x];
      for( int i = 0 ; i < x ; i++ )
      dynamicArray[i] = new T [y];

      return dynamicArray;
}

template <typename T>
void FreeDynamicArray2D(T** dArray)
{
      delete [] *dArray;
      delete [] dArray;
}

template <typename T> 
T *AllocateDynamicArray1D( int x)
{
      T *dynamicArray;

      dynamicArray = new T[x];
      return dynamicArray;
}

template <typename T>
void FreeDynamicArray1D(T* dArray)
{
      delete [] dArray;
}

template <typename T> 
T ***AllocateDynamicArray3D( int gridx, int gridy, int gridz)
{
      T ***dynamicArray;

      dynamicArray = new T**[gridx];
      for(int i=0;i<gridx;i++)  
      {		
	    dynamicArray[i] = new T*[gridy];
	    for(int j=0;j<gridy;j++)
	    {
	      dynamicArray[i][j] = new T[gridz]; 
	    }
      }

      return dynamicArray;
}

template <typename T>
void FreeDynamicArray3D(T*** dArray)
{
      delete [] *dArray;
      delete [] dArray;
}

