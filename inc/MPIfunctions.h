//Header file including all random functions needed within KALOS
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <mpi.h>
#include <cmath>

using namespace std;

//Buffer commincation for x-direction advection
void communicate_buffer_x(MPI_Comm comm1d,complex<double>*****f,complex<double>****buffer)
{

}

//Buffer commincation for y-direction advection
void communicate_buffer_y(MPI_Comm comm1d,complex<double>*****f,complex<double>****buffer)
{

}

//Buffer commincation for z-direction advection
void communicate_buffer_z(MPI_Comm comm1d,complex<double>*****f,complex<double>***buffer)
{

}

//Sets up communicators and MPI grid
void MPI_Grid_Setup()
{

}

void get_neighbours1D(MPI_Comm comm1d,int me,int pright,int pleft)
{
    int ierr;
    //Integer determining which direction -> In this case the only direction
    int dir = 0;
    //Integer determining the shift
    int shift = 1;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&me,&pright);
    //cout << "RankR" << " " << MPIrank << " " << pright << endl;
    shift = -1;
    ierr = MPI_Cart_shift(comm1d,dir,shift,&me,&pleft);  
}