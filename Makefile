
#SASVIEW_ROOT=/Users/pkienzle/Source/sasview
SASVIEW_ROOT=/home/pkienzle/src/sasview
#SASVIEW_ROOT=C:/Source/sasview
LIBIGOR_HEADERS=$(SASVIEW_ROOT)/sansmodels/src/libigor
LIBIGOR_SRC=$(SASVIEW_ROOT)/sansmodels/src/libigor
MODEL_HEADERS=$(SASVIEW_ROOT)/sansmodels/include
MODEL_SRC=$(SASVIEW_ROOT)/sansmodels/src/c_models

#CXX=C:/MinGW32-xy/bin/g++
#CC=C:/MinGW32-xy/bin/gcc
#WINFUNCS_H=winFuncs.h
#WINFUNCS_C=winFuncs.c
#WINFUNCS_O=winFuncs.o
#LIBEXT=.dll

# GNU compiler definition for linux
## use -fopenmp on CXXFLAGS/LDFLAGS for openmp
#CC=gcc
#CXX=g++
#CC=/opt/local/bin/gcc-mp-4.7
#CXX=/opt/local/bin/g++-mp-4.7
LIBEXT=.so

CCFLAGS=-Wall -O2 -fPIC -fopenmp
CXXFLAGS=-Wall -O2 -fPIC -fopenmp
INCLUDE=-I. -I$(LIBIGOR_HEADERS) -I$(MODEL_HEADERS)
LDFLAGS=-shared -fopenmp
LIBS=

%.o: %.cpp ; $(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
%.o: $(LIBIGOR_SRC)/%.c ; $(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@
%.o: $(MODEL_SRC)/%.c ; $(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@
%.so: %.o ; $(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

all: cylinder.so

disperser.o: disperser.cpp ModelInfo.h disperser.h

cylinder.o: cylinder.cpp ModelInfo.h disperser.h

cylinder$(LIBEXT): cylinder.o disperser.o GaussWeights.o \
	libCylinder.o libStructureFactor.o libfunc.o

libSphere.o: libSphere.c libSphere.h GaussWeights.h  $(WINFUNCS_H)

GaussWeights.o: GaussWeights.c GaussWeights.h

#winFuncs.o: winFuncs.c winFuncs.h

parameters.o: parameters.cpp parameters.hh

sphere.o: sphere.cpp sphere.h parameters.hh

libfunc.o: libfunc.c

libCylinder.o: libCylinder.c libCylinder.h

libStructureFactorr.o: libStructureFactor.c libStructureFactor.h

SampleModel.o: SampleModel.cpp ModelInfo.h sphere.h

SampleModel$(SHLIB_EXT): SampleModel.o sphere.o parameters.o libSphere.o GaussWeights.o $(WINFUNCS_O) libfunc.o
	$(LD) $(LDFLAGS) $^ -o $@


clean:
	-rm *.o *.obj *.lib *.exp *.so *.dll *~
	
