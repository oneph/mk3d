CXX_IFLAGS = -I/usr/local/include -I./inc
CXX_LFLAGS = -L/usr/local/lib
CXX_LIBS = -lrfftw -lfftw

CXX = g++

OPTFLAGS = -g 

CXXSOURCE = ./src/kalosmain3D3.cpp
CXXOBJECTS = kalosmain3D3.o 

EXEC = mk3d.out

$(EXEC) : $(CXXOBJECTS)
	$(CXX) $(OPTFLAGS) $(CXXOBJECTS) -o $(EXEC) $(CXX_LFLAGS) $(CXX_LIBS)

$(CXXOBJECTS) : $(CXXSOURCE)
	$(CXX) $(OPTFLAGS) $(CXX_IFLAGS) -c $(CXXSOURCE)  

clean:
	rm -f $(CXXOBJECTS) $(EXEC)



 
