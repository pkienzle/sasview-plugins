import numpy
from numpy import inf

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

        sigma = width * center if relative else width

        if sigma == 0:
            return numpy.array([center, 1.], 'd')

        x = center + numpy.linspace(-nsigmas * sigma, +nsigmas * sigma, npts)
        x = x[(x >= min) & (x <= max)]

        val = x - center

        result = numpy.empty(2 * npts, 'd')
        result[0::2] = x
        result[1::2] = numpy.exp((val * val) / (-2.0 * sigma * sigma))

        return result


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
