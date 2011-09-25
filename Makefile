CXX_IFLAGS = -I/usr/local/include -I./inc
CXX_LFLAGS = -L/usr/local/lib
CXX_LIBS = -lrfftw -lfftw

EXEC = ~/MK3d/mk3d/run/mk3d.out

CXX = mpic++

OPTFLAGS = -g -O3 

CXXSOURCE = ~/MK3d/mk3d/src/kalosmainMPI.cpp
CXXOBJECTS = ~/MK3d/mk3d/obj/kalosmainMPI.o 

$(EXEC) : $(CXXOBJECTS)
	$(CXX) $(OPTFLAGS) $(CXXOBJECTS) -o $(EXEC) $(CXX_LFLAGS) $(CXX_LIBS)

$(CXXOBJECTS) : $(CXXSOURCE)
	$(CXX) $(OPTFLAGS) $(CXX_IFLAGS) -o $(CXXOBJECTS) -c $(CXXSOURCE)   

clean:
	rm -f $(CXXOBJECTS) $(EXEC)



 
