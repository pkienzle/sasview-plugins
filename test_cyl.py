#! /usr/bin/env python

import sys
import os
import numpy
import time
#import pdb; pdb.set_trace()

def _run_one(model, q, q1, q2):
    print "run   %g=>%.5g, %g=>%.5g"%(q1, model.run(q1), q2, model.run(q2))
    print "runXY %g=>%.5g, %g=>%.5g"%(q1, model.runXY(q1), q2, model.runXY(q2))
    Iq = model.evalDistribution(numpy.array((q1,q2),'d'))
    print "eval  %g=>%.5g, %g=>%.5g"%(q1,Iq[0],q2,Iq[1])

    start=time.time()
    Iq=model.evalDistribution(q)
    dt = time.time()-start
    print "vec   %g=>%.5g, %g=>%.5g"%(q[0],Iq[0],q[-1],Iq[-1]), "time %.1f ms"%(dt*1000)

    start=time.time()
    Iq = []
    for qi in q: Iq.append(model.evalDistribution(numpy.array([qi],'d'))[0])
    dt = time.time()-start
    print "loop  %g=>%.5g, %g=>%.5g"%(q[0],Iq[0],q[-1],Iq[-1]), "time %.1f ms"%(dt*1000)

def run(name, model, Dispersion):

    disp = {'npts':15, 'width':0.1666666667, 'nsigmas':2}
    model_pars = {'scale':1.0, 'radius':60.0, 'sldCyl': 2.0, 'sldSolv': 1.0, 'background': 0.0}
    q1,q2 = 0.1, 0.2
    q = numpy.linspace(q1,q2,1000)    #print "q",q.shape,q.dtype

    for k,v in model_pars.items(): model.setParam(k,v)


    print "===", name, "without dispersion"
    _run_one(model, q, q1, q2)


    # Set up dispersion
    rdisp = Dispersion()
    model.set_dispersion('radius', rdisp)
    for k,v in disp.items(): model.dispersion['radius'][k] = v

    ldisp = Dispersion()
    #model.set_dispersion('length', ldisp)
    #for k,v in disp.items(): model.dispersion['length'][k] = v

    print "===", name, "with dispersion"
    _run_one(model, q, q1, q2)


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

#run_old()
run_new()
run_old()
#run_new()
