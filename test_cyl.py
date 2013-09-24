#! /usr/bin/env python

import sys
import os
import numpy
import time
#import pdb; pdb.set_trace()

def run(name, model, Dispersion):

    q = numpy.linspace(0.1,0.2,1000)
    disp = {'npts':10, 'width':0.1666666667, 'nsigmas':2}

    model_pars = {'scale':1.0, 'radius':60.0, 'sldCyl': 2.0, 'sldSolv': 1.0, 'background': 0.0}

    for k,v in model_pars.items(): model.setParam(k,v)

    qlist3 = numpy.asarray([0.001, 0.002])
    print name,"run", model.run(0.001), model.run(0.002)
    print name,"runXY", model.runXY(0.001), model.runXY(0.002)
    print name,"evalDist", model.evalDistribution(qlist3)

    start=time.time()
    Iq=model.evalDistribution(q)
    dt = time.time()-start
    print name, "vector",dt,Iq[-1]

    start=time.time()
    for qi in q: Iq=model.evalDistribution(numpy.array([qi]))
    dt = time.time()-start
    print name, "loop",dt,Iq



    rdisp = Dispersion()
    model.set_dispersion('radius', rdisp)
    for k,v in disp.items(): model.dispersion['radius'][k] = v

    ldisp = Dispersion()
    model.set_dispersion('length', ldisp)
    for k,v in disp.items(): model.dispersion['length'][k] = v

    start=time.time()
    Iq=model.evalDistribution(q)
    dt = time.time()-start
    print name, "disp vector",dt,Iq[-1]

    start=time.time()
    for qi in q: Iq=model.evalDistribution(numpy.array([qi]))
    dt = time.time()-start
    print name, "disp loop",dt,Iq


def run_old():
    from sans.models.CylinderModel import CylinderModel as Model
    from sans.models.dispersion_models import GaussianDispersion
    model = Model()
    run("old", model, GaussianDispersion)

def run_new():
    import cplugin
    from cplugin import GaussianDispersion
    Model = cplugin.cplugin('cylinder.so')
    model = Model()
    run("new", model, GaussianDispersion)

run_old()
run_new()
