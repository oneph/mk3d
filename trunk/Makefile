CXX_IFLAGS = -I/home/btpp/btpp05/fftw2/include -I./inc
CXX_LFLAGS = -L/home/btpp/btpp05/fftw2/lib
CXX_LIBS = -lrfftw -lfftw

CXX = mpic++

OPTFLAGS = -g -O3 

CXXSOURCE = ~/MK3d/mk3d/src/kalosmain3D3.cpp
CXXOBJECTS = ~/MK3d/mk3d/obj/kalosmain3D3.o

$(EXEC) : $(CXXOBJECTS)
	$(CXX) $(OPTFLAGS) $(CXXOBJECTS) -o $(EXEC) $(CXX_LFLAGS) $(CXX_LIBS)

$(CXXOBJECTS) : $(CXXSOURCE)
	$(CXX) $(OPTFLAGS) $(CXX_IFLAGS) -o $(CXXOBJECTS) -c $(CXXSOURCE)   

clean:
	rm -f $(CXXOBJECTS) $(EXEC)



 
