#! /usr/bin/env python

import sys
import os
import numpy
import time

import cplugin

#import pdb; pdb.set_trace()
q = numpy.arange(1000,dtype='d')
disp = {'npts':10, 'width':0.1666666667, 'nsigmas':2}
sphere = {'scale':1.0, 'radius':60.0, 'sldCyl': 2.0, 'sldSolv': 1.0, 'background': 0.0}

from sans.models.CylinderModel import CylinderModel
from sans.models.dispersion_models import GaussianDispersion
oldsm = CylinderModel()

qlist3 = numpy.asarray([0.001, 0.002])
print "run", oldsm.run(0.001), oldsm.run(0.002)
print "runXY", oldsm.runXY(0.001), oldsm.runXY(0.002)
print "evalDist", oldsm.evalDistribution(qlist3)

for k,v in sphere.items(): oldsm.setParam(k,v)
  
gdisp = GaussianDispersion()
oldsm.set_dispersion('radius', gdisp)
for k,v in disp.items(): oldsm.dispersion['radius'][k] = v
print "G run", oldsm.run(0.001), oldsm.run(0.002)
print "G runXY", oldsm.runXY(0.001), oldsm.runXY(0.002)
print "G evalDist", oldsm.evalDistribution(qlist3)

start=time.time()
print "calc old"; sys.stdout.flush()
Iq=oldsm.evalDistribution(q)
#Iq=oldsm.run(q)
dt = time.time()-start
print "oldsm",dt,Iq[len(q)/2],len(q),disp['npts']
sys.stdout.flush()

#NewSphere = cplugin.cplugin('SampleModel.dll')
NewSphere = cplugin.cplugin('cylinder.so')
newsm = NewSphere()
newsm.params = oldsm.params
#newsm.set_dispersion('radius', cplugin.GaussianDispersion(**disp))

start=time.time()
print "evaldist"
Iq=newsm.evalDistribution(q)
#Iq=newsm.run(q)
dt = time.time()-start
print "oldsm",dt,Iq[len(q)/2],len(q),disp['npts']
