#SASVIEW_ROOT=/Users/pkienzle/Source/sasview
SASVIEW_ROOT=/home/pkienzle/src/sasview
LIBIGOR_HEADERS=$(SASVIEW_ROOT)/sansmodels/src/libigor
LIBIGOR_SRC=$(SASVIEW_ROOT)/sansmodels/src/libigor
MODEL_HEADERS=$(SASVIEW_ROOT)/sansmodels/include
MODEL_SRC=$(SASVIEW_ROOT)/sansmodels/src/c_models

# GNU compiler definition for linux
## use -fopenmp on CXXFLAGS/LDFLAGS for openmp
CC=gcc
CCFLAGS=-O2 -fPIC
CXX=g++
CXXFLAGS=-Wall -O2 -fPIC
INCLUDE=-I. -I$(LIBIGOR_HEADERS) -I$(MODEL_HEADERS)
LIBEXT=.so
LDFLAGS=-shared
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



clean:
	rm *.o *.so *.dll *~
