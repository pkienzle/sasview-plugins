# ====================================================
import os
import numpy
from numpy import inf

"""
Dispersions
"""
class GaussianDispersion(object):
    
    def __init__(self, npts=35, width=0, nsigmas=3):
        self.type = 'gaussian'
        self.npts = npts
        self.width = width
        self.nsigmas = nsigmas
        
    def get_pars(self):
        return self.__dict__
        
    def get_weights(self, center, min=-inf, max=+inf, relative=False):
        """
        *center* is the center of the distribution
        
        *min*,*max* are the min, max allowed values
        
        *relative* is True if the width is relative to the center instead of absolute
        
        For polydispersity use relative.  For orientation parameters use absolute.
        """
        npts, width, nsigmas = self.npts, self.width, self.nsigmas

        sigma = width*center if relative else width
            
        if sigma == 0:
            return numpy.array([center],'d'), numpy.array([1.],'d')
        
        
        x = mean + numpy.linspace(-nsigmas*sigma, +nsigmas*sigma, npts)
        #print 'x', x
        
        val = x[(x>=min) & (x<=max)] - mean
        #print 'val', val
        
        result = numpy.exp((val*val)/(-2.0*sig*sig))
        #print 'result', result
        
        return x, result

DISPERSIONS = {
    'gaussian': GaussianDispersion,
    }

def get_weights(dispersion_pars, value, min, max, relative):
    """
    Given a set of dispersion parameters, return a set of weights
    """
    d = DISPERSIONS[dispersion_pars['type']]()
    d.__dict__ = dispersion_pars
    return d.get_weights(value, min, max, relative)
 

# ====================================================
"""
Basis for defining new types of plugin models.

Models should define ParameterInfo, ModelInfo
"""

from numpy import inf

from sans.models.BaseComponent import BaseComponent

API_VERSION = 1

# Parameter type for C call
PT_End          = 0x0
PT_Simple       = 0xAAAAAAA1
PT_Polydisperse = 0xAAAAAAA2

# Parameter flags
PF_Orientation   = 0x01
PF_Magnetic      = 0x02
PF_Unfittable    = 0x04
PF_Integer       = 0x08
PF_Polydisperse  = 0x10
PF_RepeatCount   = 0x20
PF_Repeated      = 0x40

class ParameterInfo(object):
    def __init__(self, name, description, units="", 
            default=0., dispmin=-inf, dispmax=inf, flags=0):
        self.name = name
        self.description = description
        self.units = units
        self.default = default
        self.dispmin = dispmin
        self.dispmax = dispmax
        self.flags = flags

class ModelInfo(object):
    def __init__(self, name, description, parameter_info, version=API_VERSION):
        self.name = name
        self.description = description
        self.parameter_info = parameter_info
        self.version = version

class PluginBase(BaseComponent):
    def __init__(self, model_info):
        """ Initialization"""
        BaseComponent.__init__(self)
        
        self.parameter_info  = model_info.parameter_info

        # ===== Variable state which needs to be copied/saved =====
        self.params = dict((p.name, p.default) 
                           for p in self.parameter_info)
        self.details = dict((p.name, [p.units,None,None]) 
                            for p in self.parameter_info)
        self.dispersion = dict((p.name, GaussianDispersion().get_pars()) 
                               for p in self.parameter_info
                                if p.flags&(PF_Polydisperse|PF_Orientation))
        #list of parameter that start out fixed by default
        self.fixed = []


        # ===== Fixed state that is not changed by the sasview gui =====
        ## Name of the model
        self.name = model_info.name
        self.description = model_info.description
        self.category = "unknown"
        
        self.non_fittable = [p.name for p in self.parameter_info
                             if p.flags&(PF_Unfittable|PF_RepeatCount)]
        self.orientation_params = [p.name for p in self.parameter_info
                                   if p.flags&PF_Orientation]

        self.magnetic_params = [p.name for p in self.parameter_info
                                if p.flags&PF_Magnetic]
                                    
        ## independent parameter name and unit [string]
        self.input_name = "Q"
        self.input_unit = "A^{-1}"
        ## output name and unit  [string]
        self.output_name = "Intensity"
        self.output_unit = "cm^{-1}"

    def calculate_ER(self): 
        """
        Calculate effective radius
        """
        return NotImplemented  
    
    def calculate_VR(self): 
        """
        Calculate volume fraction ratio
        """
        return NotImplemented 

    def calculate_Iq(self, *q):
        """
        Evaluate Iq for the current parameter set.
        
        Use q or qx,qy or qx,qy,qz.
        """
        return NotImplemented

    def evalDistribution(self, qdist):
        """
        Evaluate a distribution of q-values.
        
        **DEPRECATED**: use calculate_Iq instead
        
        * For 1D, a numpy array is expected as input:
        
            evalDistribution(q)
            
        where q is a numpy array.
        
        
        * For 2D, a list of numpy arrays are expected: [qx_prime,qy_prime],
          where 1D arrays,
        
        qx_prime = [ qx[0], qx[1], qx[2], ....]
        and
        qy_prime = [ qy[0], qy[1], qy[2], ....] 
        
        Then get
        q = numpy.sqrt(qx_prime^2+qy_prime^2)
        
        that is a qr in 1D array;
        q = [q[0], q[1], q[2], ....] 
        
        :Note: Due to 2D speed issue, no anisotropic scattering 
            is supported for python models, thus C-models should have
             their own evalDistribution methods.
        
        The method is then called the following way:
        
        evalDistribution(q)
        where q is a numpy array.
        
        :param qdist: ndarray of scalar q-values or list [qx,qy] 
                    where qx,qy are 1D ndarrays 
        
        """
        if isinstance(qdist, list):
            # Check whether we have a list of ndarrays [qx,qy]
            if (len(qdist)!=2 
                    or any(not isinstance(v, numpy.ndarray) for v in qdist)):
                msg = "evalDistribution expects a list of 2 ndarrays"
                raise RuntimeError, msg
            return self.calculate_Iq(*qdist)
                
        elif isinstance(qdist, numpy.ndarray):
            # We have a simple 1D distribution of q-values
            return self.calculate_Iq(qdist)
            
        else:
            mesg = "evalDistribution is expecting an ndarray of scalar q-values"
            mesg += " or a list [qx,qy] where qx,qy are 2D ndarrays."
            raise RuntimeError, mesg

    def run(self, x): 
        """
        Evaluate a single q or [q,phi]
        
        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, list) and len(x) == 2:
            q,phi = x
            return self.calculate_Iq([q*math.cos(phi)], [q*math.sin(phi)])
        else:
            return self.calculate_Iq([float(x)])
    
    def runXY(self, x): 
        """
        Evaluate a single q or [qx,qy]
        
        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, list) and len(x) == 2:
            return self.calculate_Iq(*x)
        else:
            return self.calculate_Iq([float(x)])

    def clone(self):
        """ Returns a new object identical to the current object """
        return copy.deepcopy(self)
    
    def set_dispersion(self, parameter, dispersion):
        """
        model dispersions
        """
        if parameter.lower() in (s.name.lower() for s in self.parameter_info):
            self.dispersion[parameter] = dispersion.get_pars()
        else:
            raise ValueError("%r is not a dispersity or orientation parameter")

    def _get_param_vector(self):
        pars = []
        for p in self.parameter_info:
            # TODO: need to handle multiplicity
            value = self.params[p.name]
            if p.name in self.dispersion:
                relative = (p.flags&PF_Polydisperse) != 0
                w = get_weights(self.dispersion[p.name], value, 
                                p.dispmin, p.dispmax, relative)
                pars.append(w)
            else:
                pars.append(value)
        return pars
    
# ======================================================================
"""
CPlugin wraps a DLL which exports a SAS model API.
"""
import sys

import ctypes
from ctypes import cdll, c_void_p, c_double, c_size_t, c_char_p, addressof, byref
c_double_p = ctypes.POINTER(c_double)


#from sans.models.plugin import ParameterInfo

if sys.platform == "win32":
    from _ctypes import FreeLibrary as dlclose
else:
    from _ctypes import dlclose

class CParameterInfo(ctypes.Structure):
    _fields_ = (
        ( "name", c_char_p ),
        ( "description", c_char_p ),
        ( "units", c_char_p ),
        ( "default", c_double ),
        ( "dispmin", c_double ),
        ( "dispmax", c_double ),
        ( "flags", c_size_t ),
        )

class CModelInfo(ctypes.Structure):
    _fields_ = (
        ( "version", c_size_t ),
        ( "name", c_char_p ),
        ( "description", c_char_p ),
        ( "nparameters", c_size_t),
        ( "parameters", ctypes.POINTER(CParameterInfo) ),
        )
c_model_p = ctypes.POINTER(CModelInfo)

class CEndParameter(ctypes.Structure):
    _fields_ = (
        ( "type", c_size_t ),
        )
    def __init__(self):
        self.type = PT_End
CEND_PARAMETER = CEndParameter()

class CSimpleParameter(ctypes.Structure):
    _fields_ = (
        ( "type", c_size_t ),
        ( "value", c_double ),
        )
    def __init__(self, value):
        self.type = PT_Simple
        self.value = value
        
class CPolydisperseParameter(ctypes.Structure):
    _fields_ = (
        ( "type", c_size_t ),
        ( "length", c_size_t ),
        ( "values", c_double_p ),
        ( "weights", c_double_p ),
        )
    def __init__(self, values, weights):
        self.type = PT_Polydisperse
        self.length = values.size
        self.values = values.ctypes.data_as(c_double_p)
        self.weights = weights.ctypes.data_as(c_double_p)
        print "setting disp",self.type,self.length


def pyparameters_to_c(parameters):
    """
    Convert parameter list from python to C structure array
    
    Returns a pointer array to pass to C and a list of the underlying objects
    so that they can be collected by 
    """
    # First create an array of parameter structures we can send to C.  We need 
    # to hold on to them until we are finished calling the C function.
    phandles = [(CSimpleParameter(p) if isinstance(p,float) 
                 else CPolydisperseParameter(*p))
                for p in parameters]
    # Add a null terminator to the end of the parameter list
    phandles.append(CEND_PARAMETER)
    # Get references to the parameter structures and put them in an array
    cparameters = (c_void_p * len(phandles))(*[addressof(p) for p in phandles])
    print "cpars %x"%addressof(cparameters)
    print "cpars[1] %x"%cparameters[1]
    # Return the array and the underlying list of structures
    return cparameters, phandles

def cparameter_info_to_py(nparameters, cparameters):
    """
    Convert parameter description list from C structure array to python class
    """
    parameters = []
    for index in xrange(nparameters):
        cp = cparameters[index]
        index += 1
        p = ParameterInfo(
            name=str(cp.name),
            description=str(cp.description),
            units=str(cp.units),
            default=float(cp.default),
            dispmin=float(cp.dispmin),
            dispmax=float(cp.dispmax),
            flags=int(cp.flags),
            )
        parameters.append(p)
    return parameters

def cmodel_info_to_py(cmodel_info_p):
    """
    Convert model description from C structure to python class
    """
    cmi = cmodel_info_p[0]
    return ModelInfo(
        version = int(cmi.version),
        name=str(cmi.name),
        description=str(cmi.description),
        parameter_info=cparameter_info_to_py(cmi.nparameters, cmi.parameters),
        )


create_model_prototype = ctypes.CFUNCTYPE(c_void_p, c_void_p)
def def_create_model(lib):
    """
    Wrap create_model from lib if it exists, or return alternate.
    """
    try: 
        lib_create_model = create_model_prototype(lib.create_model)
        def _create_model(self, data):
            return lib_create_model(data)
    except AttributeError:
        def _create_model(self, data):
            return None
    return _create_model

destroy_model_prototype = ctypes.CFUNCTYPE(None, c_void_p)
def def_destroy_model(lib):
    """
    Wrap destroy_model from lib if it exists, or return alternate.
    """
    try:
        lib_destroy_model = destroy_model_prototype(lib.destroy_model)
        def _destroy_model(self, handle):
            return lib_destroy_model(handle)
    except AttributeError:
        def _destroy_model(self, handle): pass
    return _destroy_model

calculate_q_prototype = ctypes.CFUNCTYPE(None, c_void_p, c_void_p, c_size_t, c_double_p, c_double_p)
def def_calculate_q(lib):
    """
    Wrap calculate_q from lib if it exists, or return alternate.
    """
    try:
        lib_calculate_q = calculate_q_prototype(lib.calculate_q)
        def calculate_q(self, cparameters, q):
            "calc q"
            q = numpy.ascontiguousarray(q,'d')
            iq = numpy.empty_like(q)
            lib_calculate_q(self.handle, byref(cparameters), 
                            iq.size, 
                            iq.ctypes.data_as(c_double_p), 
                            q.ctypes.data_as(c_double_p))
            if numpy.isnan(Iq.flat[0]):
                logging.warn(self.name+" calculate_q returns NaN")
            return iq
    except AttributeError:
        def calculate_q(self, cparameters, q): 
            return NotImplemented
    return calculate_q
        
calculate_qxqy_prototype = ctypes.CFUNCTYPE(None, c_void_p, c_void_p, c_size_t, c_double_p, c_double_p, c_double_p)
def def_calculate_qxqy(lib):
    """
    Wrap calculate_qxqy from lib if it exists, or return alternate.
    """
    try:
        lib_calculate_qxqy = calculate_qxqy_prototype(lib.calculate_qxqy)
        def calculate_qxqy(self, cparameters, qx, qy):
            "calc qxy"
            qx,qy = [numpy.ascontiguousarray(v,'d') for v in qx,qy]
            iq = numpy.empty_like(qx)
            lib_calculate_qxqy(self.handle, cparameters, 
                            iq.size, iq.ctypes, qx.ctypes, qy.ctypes,
                            )
            if numpy.isnan(Iq.flat[0]):
                logging.warn(self.name+" calculate_qxqy returns NaN")
            return iq
    except AttributeError:
        def calculate_qxqy(self, parameters, qx, qy):
            q = numpy.sqrt(qx**2 + qy**2)
            return self._calculate_q(handle, parameters, iq, q)
    return calculate_qxqy
        
calculate_qxqyqz_prototype = ctypes.CFUNCTYPE(None, c_void_p, c_void_p, c_size_t, c_double_p, c_double_p, c_double_p, c_double_p)
def def_calculate_qxqyqz(lib):
    """
    Wrap calculate_qxqyqz from lib if it exists, or return alternate.
    """
    try:
        lib_calculate_qxqyqz = calculate_qxqyqz_prototype(lib.calculate_qxqyqz)
        def calculate_qxqyqz(self, cparameters, qx, qy, qz):
            "calc qxyz"
            qx,qy,qz = [numpy.ascontiguousarray(v,'d') for v in qx,qy,qz]
            iq = numpy.empty_like(qx)
            lib_calculate_qxqyqz(self.handle, cparameters, 
                            iq.size, iq.ctypes, 
                            qx.ctypes, qy.ctypes, qz.ctypes,
                            )
            if numpy.isnan(Iq.flat[0]):
                logging.warn(self.name+" calculate_qxqyqz returns NaN")
            return iq
    except AttributeError:
        def calculate_qxqyqz(self, parameters, qx, qy, qz): 
            return self._calculate_qxqz(parameters, qx, qy)
    return calculate_qxqyqz

calculate_ER_prototype = ctypes.CFUNCTYPE(c_double, c_void_p, c_void_p)
def def_calculate_ER(lib):
    """
    Wrap calculate_ER from lib if it exists, or return alternate.
    """
    try:
        lib_calculate_ER = calculate_ER_prototype(lib.calculate_ER)
        def calculate_ER(self, cparameters):
            result = lib_calculate_ER(self.handle, cparameters)
            if numpy.isnan(result):
                logging.warn(self.name+" calculate_ER returns NaN")
            return result
    except AttributeError:
        def calculate_ER(self, parameters):
            return 0.
    return calculate_ER

calculate_VR_prototype = ctypes.CFUNCTYPE(c_double, c_void_p, c_void_p)
def def_calculate_VR(lib):
    """
    Wrap calculate_VR from lib if it exists, or return alternate.
    """
    try:
        lib_calculate_VR = calculate_VR_prototype(lib.calculate_VR)
        def calculate_VR(self, cparameters):
            result = lib_calculate_VR(self.handle, cparameters)
            if numpy.isnan(result):
                logging.warn(self.name+" calculate_VR returns NaN")
            return result
    except AttributeError:
        def calculate_VR(self, parameters):
            return 1.
    return calculate_VR

get_model_info_prototype = ctypes.CFUNCTYPE(c_model_p)
def get_model_info(lib):
    lib.get_model_info.restype = c_model_p
    #cmodel_info = get_model_info_prototype(lib.get_model_info)()
    cmodel_info = lib.get_model_info()
    model_info = cmodel_info_to_py(cmodel_info)
    return model_info
    
def cplugin(path):
    """
    Generate a SAS model from a dll.
    
    *path* : string
    
        path to the C dll, as used by ctypes
        
    *data* : c_void_p
    
        ctypes void pointer to data required by the C model on initialization
    """
    lib = ctypes.cdll[os.path.abspath(path)]
    model_info = get_model_info(lib)

    class CPlugin(PluginBase):
        handle = None # Default to no handle

        # Attach methods from the c library
        _create_model = def_create_model(lib)
        _destroy_model = def_destroy_model(lib)
        _calculate_q = def_calculate_q(lib)
        _calculate_qxqy = def_calculate_qxqy(lib)
        _calculate_qxqyqz = def_calculate_qxqyqz(lib)
        _calculate_ER = def_calculate_ER(lib)
        _calculate_VR = def_calculate_VR(lib)
        _plugin_path = path

        @staticmethod
        def unload_plugin():
            dlclose(lib._handle)
        
        def __init__(self, data=None):
            PluginBase.__init__(self, model_info)
            self.handle = self._create_model(data)

        def __del__(self):
            self.free()

        def free(self):
            if self.handle is not None: 
                #print "calling model destructor for",self.name
                self._destroy_model(self.handle)
                self.handle = None

        def calculate_Iq(self, *q):
            pars = self._get_param_vector()
            cparameters,phandles = pyparameters_to_c(pars)
            if len(q) == 1:
                print "q"
                Iq = self._calculate_q(cparameters, *q)
            elif len(q) == 2:
                print "qxy"
                Iq = self._calculate_qxqy(cparameters, *q)
            elif len(q) == 3:
                print "qxyz"
                Iq = self._calculate_qxqyqz(cparameters, *q)
            else:
                raise TypeError("calculate_Iq expects q or qx,qy or qx,qy,qz")
            return Iq

        def calculate_ER(self):
            pars = self._get_param_vector()
            cparameters,phandles = pyparameters_to_c(pars)
            return self._calculate_ER(cparameters)

        def calculate_VR(self):
            pars = self._get_param_vector()
            cparameters,phandles = pyparameters_to_c(pars)
            return self._calculate_VR(cparameters)

    CPlugin.__name__ = model_info.name
    CPlugin.__doc__ = model_info.description

    return CPlugin
