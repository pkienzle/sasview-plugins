#! /usr/bin/env python

import sys
import numpy
import time
import cplugin

#import pdb; pdb.set_trace()
q = numpy.arange(1000)
disp = {'npts':35, 'width':1, 'nsigmas':3}

from sans.models.SphereModel import SphereModel
from sans.models.dispersion_models import GaussianDispersion
oldsm = SphereModel()
#oldsm.set_dispersion('radius', GaussianDispersion())
#oldsm.dispersion['radius']['npts'] = disp['npts']
#oldsm.dispersion['radius']['width'] = disp['width']

print "run", oldsm.run(0.001), model.run(0.002)
print "runXY", oldsm.runXY(0.001), model.runXY(0.002)
qlist3 = numpy.asarray([0.001, 0.002])
print "evalDist", self.model.evalDistribution(qlist3)

tart=time.time()
print "calc old"; sys.stdout.flush()
Iq=oldsm.evalDistribution(q)
#Iq=oldsm.run(q)
dt = time.time()-start
print "oldsm",dt,Iq[len(q)/2],len(q),disp['npts']
sys.stdout.flush()

#NewSphere = cplugin.cplugin('SampleModel.dll')
NewSphere = cplugin.cplugin('SampleModel.so')
newsm = NewSphere
newsm.params = oldsm.params
#newsm.set_dispersion('radius', cplugin.GaussianDispersion(**disp))

start=time.time()
Iq=newsm.evalDistribution(q)
#Iq=newsm.run(q)
dt = time.time()-start
print "oldsm",dt,Iq[len(q)/2],len(q),disp['npts']
