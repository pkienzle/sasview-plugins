import os
import struct
from ctypes import *

import numpy

from weights import GaussianDispersion, get_weights

API_VERSION = 1

#################################################################################
## helpers

def enum(*sequential, **named): # is a helper function to create enums
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

class LibraryHandle(object): # is used to open and close external library

    # we need to distinguish between windows and other operating systems
    import _ctypes
    dlopen  = _ctypes.LoadLibrary if os.name in ("nt", "ce") else _ctypes.dlopen
    dlclose = _ctypes.FreeLibrary if os.name in ("nt", "ce") else _ctypes.dlclose

    # instance
    def __init__(self):
        self.handle = None
        
    def __del__(self):
        self.close()
            
    def __nonzero__(self):
        return self.handle is not None

    # library
    def open(self, path):
        self.close()
        self.handle = LibraryHandle.dlopen(path)
        
    def close(self):
        if self.handle is not None:
            handle, self.handle = self.handle, None
            LibraryHandle.dlclose(handle)


#################################################################################
## objects of the following types are passed to python

ParameterFlags = enum(
    Orientation   = 0x01,
    Magnetic      = 0x02,
    Unfittable    = 0x04,
    Integer       = 0x08,
    Polydisperse  = 0x10,
    RepeatCount   = 0x20 | 0x04,
    Repeated      = 0x40)

c_double_p = POINTER(c_double)
c_cmodel_p = c_void_p # pointer to unspecified data which can be used by external library (for each created c-model)

class c_parameter_info(Structure):
    _fields_ = [
        ("name"       , c_char_p),
        ("description", c_char_p),
        ("unit"       , c_char_p),
        ("default"    , c_double),
        ("dispmin"    , c_double),
        ("dispmax"    , c_double),
        ("flags"      , c_size_t)]
c_parameter_info_p = POINTER(c_parameter_info)

class c_model_info(Structure):
    _fields_ = [
        ("version"        , c_size_t),
        ("name"           , c_char_p),
        ("description"    , c_char_p),
        ("parameter_count", c_size_t),
        ("parameters"     , c_parameter_info_p)]    
c_model_info_p     = POINTER(c_model_info)


#################################################################################
## objects of the following types are passed to external library

ParameterType = enum(
    End          = 0xAAAAAAA0,
    Simple       = 0xAAAAAAA1,
    Polydisperse = 0xAAAAAAA2)
        
c_data_p       = c_void_p
c_parameters_p = c_void_p

#################################################################################
## only objects of the following types should be used to access external models

class ModelInfo(object): # describes external model
    # instance
    def __init__(self, name, description, parameters):
        self.name        = name
        self.description = description
        self.parameters  = parameters # list of ParameterInfo

        # the following lists define the type of the parameters 
        self.orientation  = [p.name for p in parameters if p.flags & ParameterFlags.Orientation ]
        self.magnetic     = [p.name for p in parameters if p.flags & ParameterFlags.Magnetic    ]
        self.unfittable   = [p.name for p in parameters if p.flags & ParameterFlags.Unfittable  ]
        self.integer      = [p.name for p in parameters if p.flags & ParameterFlags.Integer     ]
        self.polydisperse = [p.name for p in parameters if p.flags & ParameterFlags.Polydisperse]
        
class ParameterInfo(object): # ModelInfo.parameters contains ParameterInfo for each parameter
    # instance
    def __init__(self, name, description, unit, default, dispmin, dispmax, flags):
        self.name        = name
        self.description = description
        self.unit        = unit
        self.default     = default
        self.dispmin     = dispmin
        self.dispmax     = dispmax
        self.flags       = flags


from sans.models.BaseComponent import BaseComponent
class PluginBase(BaseComponent):
    def __init__(self, model_info, parameter_collection):
        """ Initialization"""
        BaseComponent.__init__(self)

        self.model_info = model_info
        pars = model_info.parameters

        # ===== Variable state which needs to be copied/saved =====
        self.params = dict((p.name, p.default) for p in pars)
        self.details = dict((p.name, [p.unit, None, None]) for p in pars)
        self.dispersion = dict((p.name, GaussianDispersion().get_pars()) for p in pars
                               if p.flags & ParameterFlags.Polydisperse)
        #list of parameter that start out fixed by default
        self.fixed = []

        # ===== Fixed state that is not changed by the sasview gui =====
        ## Name of the model
        self.name = model_info.name
        self.description = model_info.description

        self.non_fittable = [p.name for p in pars
                             if p.flags & (ParameterFlags.Unfittable | ParameterFlags.RepeatCount)]
        self.orientation_params = [p.name for p in pars
                                   if p.flags & ParameterFlags.Orientation]

        self.magnetic_params = [p.name for p in pars
                                if p.flags & ParameterFlags.Magnetic]

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
        if isinstance(qdist, (list,tuple)):
            # Check whether we have a list of ndarrays [qx,qy]
            if (len(qdist) != 2
                or any(not isinstance(v, numpy.ndarray) for v in qdist)):
                msg = "evalDistribution expects a list of 2 ndarrays"
                raise RuntimeError(msg)
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
            q, phi = x
            return self.calculate_Iq(numpy.array([q * math.cos(phi)], 'd'), 
                                     numpy.array([q * math.sin(phi)], 'd'))[0]
        else:
            return self.calculate_Iq(numpy.array([float(x)],'d'))[0]

    def runXY(self, x):
        """
        Evaluate a single q or [qx,qy]
        
        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list,tuple)) and len(x) == 2:
            return self.calculate_Iq(*x)[0]
        else:
            return self.calculate_Iq([float(x)])[0]

    def clone(self):
        """ Returns a new object identical to the current object """
        return copy.deepcopy(self)

    def set_dispersion(self, parameter, dispersion):
        """
        model dispersions
        """
        if parameter.lower() in (s.name.lower() for s in self.model_info.parameters):
            self.dispersion[parameter] = dispersion.get_pars()
        else:
            raise ValueError("%r is not a dispersity or orientation parameter")

    def calculate_Iq(self, *q):
        if len(q) == 1:
            iq = self.calculate_q(*q)
        elif len(q) == 2:
            iq = self.calculate_qxqy(*q)
        elif len(q) == 3:
            iq = self.calculate_qxqyqz(*q)
        else:
            raise TypeError("calculate_Iq expects q or qx,qy or qx,qy,qz")
        return iq


class PluginModel(PluginBase): # represents a concrete model with all its parameters. It's used for simulations.
    # instance
    def __init__(self, factory, id, model_info, parameters):
        self.factory    = factory    # factory object which created this PluginModel
        self.id         = id         # instance id
        PluginBase.__init__(self, model_info, parameters)
        #self.model_info = model_info # should be of type ModelInfo
        #self.parameters = parameters # instance of PluginModelParameterCollection

    def __del__(self):
        self.destroy()
        
    # model information
    def get_model_info(self):
        return self.model_info
        
    # model instantiation
    def destroy(self):
        self.factory.destroy_model(self)

    # calculations
    def calculate_q(self, q):
        return self.factory.calculate_q(self, q)
        
    def calculate_qxqy(self, qx, qy):
        return self.factory.calculate_qxqy(self, qx, qy)
        
    def calculate_qxqyqz(self, qx, qy, qz):
        return self.factory.calculate_qxqyqz(self, qx, qy, qz)
        
    def calculate_ER(self):
        return self.factory.calculate_ER(self)
        
    def calculate_VR(self):
        return self.factory.calculate_VR(self)

class PluginModelParameterCollection(object): # allows access to parameters either as PluginModel.parameters.name or PluginModel.parameters["name"]
    # instance
    def __init__(self, parameters):
        self.__dict__.update(parameters)        
    def __len__(self):
        return len(self.__dict__)
    def __getattr__(self, name):
        return self.__dict__[name]
    def __setattr__(self, name, value):
        if not name in self.__dict__:
            raise AttributeError(name)
        self.__dict__[name] = value
    def __delattr__(self, name):
        raise Exception()
    def __getitem__(self, name):
        return self.__dict__[name]
    def __setitem__(self, name, value):
        if not name in self.__dict__:
            raise AttributeError(name)
        self.__dict__[name] = value
    def __iter__(self):
        return iter(self.__dict__)
    def items(self): return self.__dict__.items()
    def values(self): return self.__dict__.values()
    def keys(self): return self.__dict__.keys()
    def update(self, *args, **kw): return self.__dict__.update(*args, **kw)
    
class PolydisperseParameter(object): # if a parameter is flagged as polydisperse then PluginModel.parameters.x will have "values" and "weights" attributes
    # instance
    def __init__(self, values, weights=None):
        self.values = values
        if weights is not None:
            self.weights = weights
        else:
            w = 1.0 / len(values)
            self.weights = [w for v in values]

class PluginModelFactory(object): # does the hard work

    # instance
    def __init__(self, path=None):
        # library
        self.path      = None               # path to loaded external library
        self._modelLib = LibraryHandle()    # handle to external library
        self._cdll     = None               # helper object which provides access to external methods for loaded library
        # functions
        self._get_model_info   = None
        self._create_model     = None
        self._destroy_model    = None
        self._calculate_q      = None
        self._calculate_qxqy   = None
        self._calculate_qxqyqz = None
        self._calculate_ER     = None
        self._calculate_VR     = None
        # created models
        self._next_model_id  = 1            # every model created will get a new id
        self._created_models = {}           # id -> c-model (used to allow us to unload current library on demand)
        
        # load library
        if path is not None:
            self.load(path)
 
    def __del__(self):
        self.unload()
        
    # load and unload library
    def load(self, path):
        if self._modelLib:
            self.unload()

        # open library
        self._created_models = {}
        self._modelLib.open(path)
        self._cdll = CDLL(None, handle=self._modelLib.handle)
        self.path  = path
        try:
            def loadfunction(cdll, name, restype, argtypes, default=None):
                try:
                    f = cdll[name]
                    f.restype  = restype
                    f.argtypes = argtypes
                    return f
                except:
                    if default:
                        return default
                    raise

            def default_create_model(data):
                return None
            def default_destroy_model(cmodel):
                pass
            def default_calculate(cmodel, cparameter_ptrs, n, iq_data, qx=None, qy=None, qz=None):
                nan = float('nan')
                for k in xrange(n):
                    iq_data[k] = nan
                    
            # load functions
            self._get_model_info   = loadfunction(self._cdll, 'get_model_info'  , c_model_info_p, [])
            # model instantiation
            self._create_model     = loadfunction(self._cdll, 'create_model'    , c_cmodel_p, [c_data_p  ], default=default_create_model )
            self._destroy_model    = loadfunction(self._cdll, 'destroy_model'   , None      , [c_cmodel_p], default=default_destroy_model)
            # I/Q calculations
            self._calculate_q      = loadfunction(self._cdll, 'calculate_q'     , None, [c_cmodel_p, c_parameters_p, c_size_t, c_double_p, c_double_p                        ], default=default_calculate)
            self._calculate_qxqy   = loadfunction(self._cdll, 'calculate_qxqy'  , None, [c_cmodel_p, c_parameters_p, c_size_t, c_double_p, c_double_p, c_double_p            ], default=default_calculate)
            self._calculate_qxqyqz = loadfunction(self._cdll, 'calculate_qxqyqz', None, [c_cmodel_p, c_parameters_p, c_size_t, c_double_p, c_double_p, c_double_p, c_double_p], default=default_calculate)
            # other calculations
            self._calculate_ER     = loadfunction(self._cdll, 'calculate_ER'    , c_double, [c_cmodel_p, c_parameters_p])
            self._calculate_VR     = loadfunction(self._cdll, 'calculate_VR'    , c_double, [c_cmodel_p, c_parameters_p])
        except:
            try:
                self.unload()
            except:
                pass
            raise
        
    def unload(self):
        # destroy existing c-models
        if self._destroy_model is not None:
            for cmodel in self._created_models.itervalues():
                self._destroy_model(cmodel)
        self._created_models = {}
        # reset functions
        self._get_model_info   = None
        self._create_model     = None
        self._destroy_model    = None
        self._calculate_q      = None
        self._calculate_qxqy   = None
        self._calculate_qxqyqz = None
        self._calculate_ER     = None
        self._calculate_VR     = None
        # close library
        self._modelLib.close()
        self._cdll = None
        self.path  = None

    # model information
    def get_model_info(self): # generates an instance of ModelInfo
        # get model info
        cmi = self._get_model_info().contents
        if cmi.version != API_VERSION:
            raise Exception()

        # get parameter info
        parameters = []
        for i in xrange(cmi.parameter_count):
            parameter = cmi.parameters[i]
            parameters.extend([ParameterInfo(
                parameter.name,
                parameter.description,
                parameter.unit,
                parameter.default,
                parameter.dispmin,
                parameter.dispmax,
                parameter.flags)])
                                
        return ModelInfo(
            cmi.name,
            cmi.description,
            parameters)

    # model instantiation
    def create_model(self, data=None): # creates a concrete model (PluginModel) which can have an individual set of parameter values
        if self._create_model is None:
            raise Exception()

        # increment id
        current_id          = self._next_model_id
        self._next_model_id += 1
        
        # create cmodel
        self._created_models[current_id] = self._create_model(data)

        model_info         = self.get_model_info()
        default_parameters = PluginModelParameterCollection({
            p.name : (p.default if not p.flags & ParameterFlags.Polydisperse else PolydisperseParameter([p.default]))
            for p in model_info.parameters})

        return PluginModel(self, current_id, model_info, default_parameters)
        
    def destroy_model(self, model): # destroys a concrete model
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._destroy_model is None:
            raise Exception()
        
        try:
            cmodel = self._created_models[model.id]
            self._destroy_model(cmodel)
        finally:
            self._created_models.pop(model.id)
            model.id         = None
            model.factory    = None
            model.model_info = None
            model.parameters = None
        
    # helper
    def _get_cparameters(self, model): # creates a c-array which holds parameter values (expected by c-model)
        is32bit = sizeof(c_size_t) == 4

        #print model_info
        #print parameters
        # determine size and offsets
        data_size = 0
        offsets   = []
        poly = {}
        for p in model.model_info.parameters:
            offsets.append(data_size)
            if p.name not in model.dispersion:
                data_size += sizeof(c_size_t) # size_t type = ParameterType.Simple;
                data_size += sizeof(c_double) # double value;
            else:
                relative = (p.flags & ~ParameterFlags.Polydisperse) != 0
                w = get_weights(model.dispersion[p.name], model.params[p.name],
                                p.dispmin, p.dispmax, relative)
                poly[p.name] = w
                npoints = len(w)/2
                
                data_size += sizeof(c_size_t)            # size_t type = ParameterType.Polydisperse;
                data_size += sizeof(c_size_t)            # size_t npoints;
                data_size += sizeof(c_double) * npoints  # double values[npoints];
                data_size += sizeof(c_double) * npoints  # double weights[npoints];

        offsets.append(data_size)
        data_size += sizeof(c_size_t) # size_t type = ParameterType.End;

        # determine header size
        header_size  = sizeof(c_size_t)                # size_t count;
        header_size += sizeof(c_size_t) * len(offsets) # size_t offsets[count];

        # create buffer
        buffer = create_string_buffer(header_size + data_size)
        # write header
        if is32bit:
            struct.pack_into('I%iI' % len(offsets), buffer, 0, len(offsets), *offsets)
        else:
            struct.pack_into('Q%iQ' % len(offsets), buffer, 0, len(offsets), *offsets)
        # wrtie data
        offset = header_size
        for p in model.model_info.parameters:
            if p.name not in model.dispersion:
                if is32bit:
                    struct.pack_into('I', buffer, offset, ParameterType.Simple)
                else:
                    struct.pack_into('Q', buffer, offset, ParameterType.Simple)
                offset += sizeof(c_size_t) # size_t type = ParameterType.Simple;

                struct.pack_into('=d', buffer, offset, model.params[p.name])
                offset += sizeof(c_double) # double value;
            else:
                w = poly[p.name]
                npoints = len(w)/2

                if is32bit:
                    struct.pack_into('II', buffer, offset, ParameterType.Polydisperse, npoints)
                else:
                    struct.pack_into('QQ', buffer, offset, ParameterType.Polydisperse, npoints)
                offset += sizeof(c_size_t) # size_t type = ParameterType.Polydisperse;
                offset += sizeof(c_size_t) # size_t npoints;

                format = '=%id' % npoints
                struct.pack_into(format, buffer, offset, *w[0::2])
                offset += sizeof(c_double) * npoints  # double values[npoints];
                struct.pack_into(format, buffer, offset, *w[1::2])
                offset += sizeof(c_double) * npoints  # double weights[npoints];

        # write end
        if is32bit:
            struct.pack_into('I', buffer, offset, ParameterType.End)
        else:
            struct.pack_into('Q', buffer, offset, ParameterType.End)

        return buffer
    
    # I/Q calculations
    def calculate_q(self, model, q):
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._calculate_q is None:
            raise Exception()

        cmodel      = self._created_models[model.id]
        cparameters = self._get_cparameters(model)
        
        if q is None:
            self._calculate_q(cmodel, cparameters, 0, None, None)
            return []

        n       = len(q)
        iq_data = (c_double * n)()
        q_data  = (c_double * n)(*q)
        self._calculate_q(cmodel, cparameters, n, iq_data, q_data)
        return list(iq_data)
        
    def calculate_qxqy(self, model, qx, qy):
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._calculate_qxqy is None:
            raise Exception()

        cmodel      = self._created_models[model.id]
        cparameters = self._get_cparameters(model)

        if (qx is None) or (qy is None):
            self._calculate_qxqy(cmodel, cparameters, 0, None, None, None)
            return []

        try: nx = len(qx)
        except: qx,nx = [qx],1
        try: ny = len(qy)
        except: qy,ny = [qy],1
            
        if nx != ny:
            raise Exception()
        
        n = nx
        iq_data = (c_double * n)()
        qx_data = (c_double * n)(*qx)
        qy_data = (c_double * n)(*qy)
        self._calculate_qxqy(cmodel, cparameters, n, iq_data, qx_data, qy_data)
        return list(iq_data)
        
    def calculate_qxqyqz(self, model, qx, qy, qz):
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._calculate_qxqyqz is None:
            raise Exception()

        cmodel      = self._created_models[model.id]
        cparameters = self._get_cparameters(model)

        if (qx is None) or (qy is None) or (qz is None):
            self._calculate_qxqyqz(cmodel, cparameters, 0, None, None, None, None)
            return []

        nx = len(qx)
        ny = len(qy)
        nz = len(qz)
        if (nx != ny) or (nx != nz):
            raise Exception()

        n = nx
        iq_data = (c_double * n)()
        qx_data = (c_double * n)(*qx)
        qy_data = (c_double * n)(*qy)
        qz_data = (c_double * n)(*qz)
        self._calculate_qxqyqz(cmodel, cparameters, n, iq_data, qx_data, qy_data, qz_data)
        return list(iq_data)

    # other calculations
    def calculate_ER(self, model):
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._calculate_ER is None:
            raise Exception()
        
        cmodel      = self._created_models[model.id]
        cparameters = self._get_cparameters(model)
        
        return self._calculate_ER(cmodel, cparameters)
        
    def calculate_VR(self, model):
        if not model.id in self._created_models:
            raise ValueError('model.id')
        if self._calculate_VR is None:
            raise Exception()
        
        cmodel      = self._created_models[model.id]
        cparameters = self._get_cparameters(model)
        
        return self._calculate_VR(cmodel, cparameters)

