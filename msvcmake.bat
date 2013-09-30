@REM as ctypes library
@REM cl cylinder.cpp disperser.cpp libCylinder.c GaussWeights.c libStructureFactor.c libfunc.c winFuncs.c /I. /LD /Ox /EHsc /Fecylinder.so
@
@REM as python extension
cl cylinder.cpp disperser.cpp libCylinder.c GaussWeights.c libStructureFactor.c libfunc.c winFuncs.c /I. /LD /Ox /EHsc /IC:/Python27/include c:\Python27\libs\python27.lib /Fe_cylinder.pyd

