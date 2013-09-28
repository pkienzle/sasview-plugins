// Volume parameters are computed first to reduce the total computational 
// effort.  Rearrange the looporder parameters so that they do come first.
#include <cmath>
#include <algorithm>
#include <cstdio>

#include <ModelInfo.h>
#include <disperser.h>

#if 0
#include <float.h>
namespace std {
inline bool isfinite(const double c) { return _finite(c) != 0; }
} ;
#endif

/**
 * Initialize the disperse, putting volume parameters first.
 */
DisperserModel::DisperserModel(const ModelInfo &m) : looporder(m.ParameterCount) {
  npars = m.ParameterCount;
  for (unsigned int i=0; i < npars; i++) looporder[i] = i;

  // Move volume parameters to the beginning
  volume_loops = 0;
  for (unsigned int i=0; i < npars; i++) {
    if (m.Parameters[i].Flags&PF_Volume) {
      std::swap(looporder[volume_loops++],looporder[i]);
    }
  }

  // Move magnetism parameters to the end
  magnetic_loops = npars;
  for (unsigned int i=npars; i > volume_loops; i--) {
    if (m.Parameters[i-1].Flags&PF_Magnetic) {
        std::swap(looporder[--magnetic_loops],looporder[i-1]);
    }
  }

  // Move non-magnetism orientation parameters before the magnetism parameters
  orientation_loops = npars;
  for (unsigned int i=magnetic_loops; i > volume_loops; i--) {
    if (m.Parameters[i-1].Flags&PF_Orientation) {
        std::swap(looporder[--orientation_loops],looporder[i-1]);
    }
  }
}

double
DisperserModel::formQxy(const double dp[], double qx, double qy) const
{
    return formQ(dp, std::sqrt(qx*qx+qy*qy));
}

double
DisperserModel::formQxyz(const double dp[], double qx, double qy, double qz) const
{
    return formQxy(dp, qx, qy);
}

void 
Disperser::calc_Q(int nq, const double q[], double Iq[]) {
  _nq = nq;
  _iq = Iq;
  _qx = q;
  _target = IQ;
  loop_Iq();
}

void 
Disperser::calc_Qxy(int nq, const double qx[], const double qy[], double Iq[]) {
  _nq = nq;
  _iq = Iq;
  _qx = qx;
  _qy = qy;
  _target = IQXY;
  loop_Iq();
}

void 
Disperser::calc_Qxyz(int nq, const double qx[], const double qy[], const double qz[], double Iq[]) {
  _nq = nq;
  _iq = Iq;
  _qx = qx;
  _qy = qy;
  _qz = qz;
  _target = IQXYZ;
  loop_Iq();
}

double
Disperser::calc_VR() {
  _volume = 1;
  _Vnorm = (_model.volume_loops == 0 ? 1. : 0.);
  _target = VR;
  loop_par(0, 1.);
  return _Vnorm;
}

double
Disperser::calc_ER() {
  _volume = 1;
  _Vnorm = (_model.volume_loops == 0 ? 1. : 0.);
  _target = ER;
  loop_par(0, 1.);
  return _Vnorm;
}

/**
 * Process an I(Q) calculation.  This clears the Iq result
 * vector, calls loop_par to process the parameters one by
 * one, then normalizes the result.
 */ 
void
Disperser::loop_Iq(void) {
//std::printf("loopQ enter\n");
  _volume = 1; // in case there is no volume normalization
  _Vnorm = (_model.volume_loops == 0 ? 1. : 0.);

  // clear
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _iq[i] = 0.;

  // loop
  loop_par(0, 1.);

  // scale
  const double weight = (_Vnorm == 0. ? 1./_Wnorm : 1./_Vnorm);
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _iq[i] *= weight;
//std::printf("loopQ exit\n");
}

/**
 * Process one level of the dispersion loop.  The loop
 * parameters are ordered so that those which affect the
 * volume come first.  This way the volume calculation
 * does not need to be repeated for those parameters that
 * do not affect the volume.
 */
void 
Disperser::loop_par(unsigned int loop, double weight)
{
//std::printf("loop %d %g\n",loop,weight);
  int p = _model.looporder[loop];
  int n = (p == 0 ? _endpts[0] : (_endpts[p] - _endpts[p-1]));
  int offset = (p == 0 ? 0 : _endpts[p-1]);
  const Weight *w = _pars + offset;
  for (int i=0; i < n; i++) {
    // Set parameter value and weight for this level
    double wi = w[i].weight * weight;
    _dp[p] = w[i].value;
    // Set volume and accumulate volume normalization
    if (loop+1 == _model.volume_loops) {
      if (_target == VR || _target == IQ) {
        _volume = _model.formV(&_dp[0]);
        _Vnorm += wi * _volume;
      } else if (_target == ER) {
        _volume = _model.formER(&_dp[0]);
        _Vnorm += wi * _volume;
      }
    } 
    // Break if just computing ER/VR
    if (loop == _model.volume_loops) {
        if (_target == ER || _target == VR) break;
    } 
    // Go the the next loop level; when computing Iq, ignore
    // orientation and magnetism parameters.
    if (loop < _model.npars 
        || (_target == IQ && loop < _model.orientation_loops)) {
      loop_par(loop+1, wi);
    } else {
      _Wnorm += wi;

      const double weight = wi*_volume;
      const double* dp = &_dp[0];
      if (_target == IQ) {
        #pragma omp parallel for
        for (int i=0; i < _nq; i++) {
          const double result = _model.formQ(dp, _qx[i]);
          // ignore singular results from kernel
          if (std::isfinite(result)) _iq[i] += result * weight;
        }
      } else if (_target == IQXY) {
        #pragma omp parallel for
        for (int i=0; i < _nq; i++) {
          const double result = _model.formQxy(dp, _qx[i], _qy[i]);
          // if (i==0) std::printf("q:%g,%g result:%g weight:%g\n",_qx[i],_qy[i],result,weight);
          // ignore singular results from kernel
          if (std::isfinite(result)) _iq[i] += result * weight;
        }
      } else if (_target == IQXYZ) {
        #pragma omp parallel for
        for (int i=0; i < _nq; i++) {
          const double result = _model.formQxyz(dp, _qx[i], _qy[i], _qz[i]);
          // ignore singular results from kernel
          if (std::isfinite(result)) _iq[i] += result * weight;
        }
      }
    }
  }
}
