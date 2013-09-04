CXX=g++
CC=gcc
SHLIB_EXT=.so
#CXX=C:/MinGW32-xy/bin/g++
#CC=C:/MinGW32-xy/bin/gcc
#EXT=.dll
#WINFUNCS_H=winFuncs.h
#WINFUNCS_C=winFuncs.c
#WINFUNCS_O=winFuncs.o

# GNU compiler flags
CXXFLAGS=-c -O2 -fopenmp -fPIC -Wall
CFLAGS=-c -O2 -fopenmp -fPIC -Wall
LD=$(CXX)
LDFLAGS=-shared
RM=rm

.c.o: ; $(CC) $(CFLAGS) $< -o $@
.cpp.o: ; $(CXX) $(CXXFLAGS) $< -o $@

all: SampleModel$(SHLIB_EXT)

clean:
	$(RM) *.o *$(SHLIB_EXT)

libSphere.o: libSphere.c libSphere.h GaussWeights.h  $(WINFUNCS_H)

GaussWeights.o: GaussWeights.c GaussWeights.h

#winFuncs.o: winFuncs.c winFuncs.h

parameters.o: parameters.cpp parameters.hh

sphere.o: sphere.cpp sphere.h parameters.hh

libfunc.o: libfunc.c

SampleModel.o: SampleModel.cpp ModelInfo.h sphere.h

SampleModel$(SHLIB_EXT): SampleModel.o sphere.o parameters.o libSphere.o GaussWeights.o $(WINFUNCS_O) libfunc.o
	$(LD) $(LDFLAGS) $^ -o $@

