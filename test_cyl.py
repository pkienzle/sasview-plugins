#! /usr/bin/env python

import sys
import os
import numpy
import time

import cplugin

#import pdb; pdb.set_trace()
q = numpy.linspace(0.1,0.2,1000)
disp = {'npts':10, 'width':0.1666666667, 'nsigmas':2}

# test model
sphere = {'scale':1.0, 'radius':60.0, 'sldCyl': 2.0, 'sldSolv': 1.0, 'background': 0.0}
from sans.models.CylinderModel import CylinderModel as OldModel
NewModel = cplugin.cplugin('cylinder.so')


from sans.models.dispersion_models import GaussianDispersion
oldsm = OldModel()

qlist3 = numpy.asarray([0.001, 0.002])
print "run", oldsm.run(0.001), oldsm.run(0.002)
print "runXY", oldsm.runXY(0.001), oldsm.runXY(0.002)
print "evalDist", oldsm.evalDistribution(qlist3)

for k,v in sphere.items(): oldsm.setParam(k,v)
  
rdisp = GaussianDispersion()
oldsm.set_dispersion('radius', rdisp)
for k,v in disp.items(): oldsm.dispersion['radius'][k] = v

ldisp = GaussianDispersion()
oldsm.set_dispersion('length', ldisp)
for k,v in disp.items(): oldsm.dispersion['length'][k] = v

start=time.time()
Iq=oldsm.evalDistribution(q)
dt = time.time()-start
print "oldsm vector",dt,Iq[-1]

start=time.time()
for qi in q: Iq=oldsm.evalDistribution(numpy.array([qi]))
dt = time.time()-start
print "oldsm loop",dt,Iq

newsm = NewModel()
newsm.params = oldsm.params
newsm.set_dispersion('radius', cplugin.GaussianDispersion(**disp))
newsm.set_dispersion('length', cplugin.GaussianDispersion(**disp))

start=time.time()
Iq=newsm.evalDistribution(q)
dt = time.time()-start
print "newsm vector",dt,Iq[-1]

start=time.time()
for qi in q: Iq=newsm.evalDistribution(numpy.array([qi]))
#for qi in q: Iq=newsm.calculate_Iq(numpy.array([qi]))
dt = time.time()-start
print "newsm loop",dt,Iq

