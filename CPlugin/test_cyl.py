#! /usr/bin/env python

import sys
import os
import numpy
import time
#import pdb; pdb.set_trace()

def _run_one(model, q, q1, q2):
    print "run     %g=>%.5g, %g=>%.5g"%(q1, model.run(q1), q2, model.run(q2))
    Iq = model.evalDistribution(numpy.array((q1,q2),'d'))
    print "eval    %g=>%.5g, %g=>%.5g"%(q1,Iq[0],q2,Iq[1])

    start=time.time()
    Iq=model.evalDistribution(q)
    dt = time.time()-start
    print "128     %g=>%.5g, %g=>%.5g"%(q[0],Iq[0],q[-1],Iq[-1]), "time %.1f ms"%(dt*1000)

    start=time.time()
    Iq = [model.evalDistribution(numpy.array((qi,),'d'))[0] for qi in q]
    dt = time.time()-start
    print "1:128   %g=>%.5g, %g=>%.5g"%(q[0],Iq[0],q[-1],Iq[-1]), "time %.1f ms"%(dt*1000)

    start=time.time()
    Iq = []
    for _ in range(10): Iq = model.evalDistribution(q)
    dt = time.time()-start
    print "128:10  %g=>%.5g, %g=>%.5g"%(q[0],Iq[0],q[-1],Iq[-1]), "time %.1f ms"%(dt*1000)

def _run_two(model, q, q1, q2):
    print "runXY   %g,%g=>%.5e, %g,%g=>%.5e"%(q1, q1, model.runXY([q1,q1]), q2, q2, model.runXY([q2,q2]))
    Iq = model.evalDistribution([numpy.array((q1,q2),'d'), numpy.array((q1,q2),'d')])
    print "evalXY  %g,%g=>%.5e, %g,%g=>%.5e"%(q1,q1,Iq[0],q2,q2,Iq[1])

    qx,qy = numpy.meshgrid(q,q)
    start=time.time()
    Iq=numpy.ascontiguousarray(model.evalDistribution([qx.flatten(),qy.flatten()]),'d').reshape(qx.shape)
    dt = time.time()-start
    print "128x128 %g,%g=>%.5e, %g,%g=>%.5e"%(qx[0,0],qy[0,0],Iq[0,0],qx[-1,-1],qy[-1,-1],Iq[-1,-1]), "time %.1f ms"%(dt*1000)

def run(name, Model, pars, Dispersion):

    disp = {'npts':15, 'width':0.1666666667, 'nsigmas':2}
    q1,q2 = 0.1, 0.2
    q = numpy.linspace(q1,q2,128)    #print "q",q.shape,q.dtype
    model = Model()
    for k,v in pars.items(): model.setParam(k,v)


    print "===", name, "without polydispersity"
    _run_one(model, q, q1, q2)
    _run_two(model, q, q1, q2)


    # Set up dispersion
    rdisp = Dispersion()
    model.set_dispersion('radius', rdisp)
    for k,v in disp.items(): model.dispersion['radius'][k] = v

    ldisp = Dispersion()
    model.set_dispersion('length', ldisp)
    for k,v in disp.items(): model.dispersion['length'][k] = v

    print "===", name, "with %dx%d polydispersity"%(disp['npts'],disp['npts'])
    _run_one(model, q, q1, q2)
    _run_two(model, q, q1, q2)


cyl_pars = {'scale':1.0, 'radius':60.0, 'sldCyl': 2.0e-6, 'sldSolv': 1.0e-6, 'background': 0.0}
def cyl_old():
    from sans.models.CylinderModel import CylinderModel as Model
    from sans.models.dispersion_models import GaussianDispersion
    run("old", Model, cyl_pars, GaussianDispersion)

def cyl_pk():
    import cplugin
    from cplugin import GaussianDispersion
    Model = cplugin.cplugin('cylinder.so')
    run("pk", Model, cyl_pars, GaussianDispersion)

def cyl_py():
    import cplugin
    from cplugin import GaussianDispersion
    import sys; sys.path.append('.')
    import cylinder
    Model = cplugin.pycplugin(cylinder)
    run("py", Model, cyl_pars, GaussianDispersion)

def cyl_dm():
    import PluginModel
    from PluginModel import PluginModelFactory
    from weights import GaussianDispersion
    factory = PluginModelFactory('cylinder.so')
    Model = factory.create_model
    run('dm', Model, cyl_pars, GaussianDispersion)
    

#cyl_new()
#cyl_old()
#cyl_py()
cyl_dm()
