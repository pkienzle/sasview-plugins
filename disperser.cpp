#include <cmath>
#include <algorithm>
#include <ModelInfo.h>
#include <disperser.h>

DisperserModel::DisperserModel(const ModelInfo &m) : looporder(m.ParameterCount) {
  // Volume parameters are computed first to reduce the total computational 
  // effort.  Rearrange the looporder parameters so that they do come first.
  //std::cout << "in disperserModel constructor" << std::endl << std::flush;
  //std::cout << "par count " << m.ParameterCount << std::endl << std::flush;
  for (unsigned int i=0; i < m.ParameterCount; i++) looporder[i] = i;
  //std::cout << "zeroed" << std::endl << std::flush;
  vloops = 0;
  for (unsigned int i=0; i < m.ParameterCount; i++) {
    if (m.Parameters[i].Flags&PF_Volume) {
      //std::cout << "swapping " << vloops << " with " << i << std::endl << std::flush;
      std::swap(looporder[vloops++],looporder[i]);
    }
  }
  //std::cout << "inited" << std::endl << std::flush;
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

class LoopQ : public QLooper  {
private:
  const int _nq;
  const double *_q;
  double *_Iq;
public: 
  LoopQ(int nq, const double q[], double Iq[]) 
    : _nq(nq), _q(q), _Iq(Iq) {}
  void clear(void);
  void scale(double weight);
  void loop(const DisperserModel& m, const double dp[], double weight);
} ;
void LoopQ::clear(void) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] = 0.;
}
void LoopQ::scale(double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] *= weight;
}
void LoopQ::loop(const DisperserModel& m, const double dp[], double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) {
    const double result = m.formQ(dp, _q[i]);
    // ignore singular results from kernel 
    if (std::isfinite(result)) _Iq[i] += result * weight;
  }
}

class LoopQxy : public QLooper  {
private:
  const int _nq;
  const double *_qx;
  const double *_qy;
  double *_Iq;
public: 
  LoopQxy(int nq, const double qx[], const double qy[], double Iq[]) 
    : _nq(nq), _qx(qx), _qy(qy), _Iq(Iq) {}
  void clear(void);
  void scale(double weight);
  void loop(const DisperserModel& m, const double dp[], double weight);
} ;
void LoopQxy::clear(void) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] = 0.;
}
void LoopQxy::scale(double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] *= weight;
}
void LoopQxy::loop(const DisperserModel& m, const double dp[], double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) {
    const double result = m.formQxy(dp, _qx[i], _qy[i]);
    // ignore singular results from kernel 
    if (std::isfinite(result)) _Iq[i] += result * weight;
  }
}

class LoopQxyz : public QLooper  {
private:
  const int _nq;
  const double *_qx;
  const double *_qy;
  const double *_qz;
  double *_Iq;
public: 
  LoopQxyz(int nq, const double qx[], const double qy[], const double qz[], double Iq[]) 
    : _nq(nq), _qx(qx), _qy(qy), _qz(qz), _Iq(Iq) {}
  void clear(void);
  void scale(double weight);
  void loop(const DisperserModel& m, const double dp[], double weight);
} ;
void LoopQxyz::clear(void) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] = 0.;
}
void LoopQxyz::scale(double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) _Iq[i] *= weight;
}
void LoopQxyz::loop(const DisperserModel& m, const double dp[], double weight) {
  #pragma omp parallel for
  for (int i=0; i < _nq; i++) {
    const double result = m.formQxyz(dp, _qx[i], _qy[i], _qz[i]);
    // ignore singular results from kernel 
    if (std::isfinite(result)) _Iq[i] += result * weight;
  }
}

void 
Disperser::calc_Q(int nq, const double q[], double Iq[]) {
  _looper = new LoopQ(nq, q, Iq);
  loop_Iq();
  delete _looper;
}

void 
Disperser::calc_Qxy(int nq, const double qx[], const double qy[], double Iq[]) {
  _looper = new LoopQxy(nq, qx, qy, Iq);
  loop_Iq();
  delete _looper;
}

void 
Disperser::calc_Qxyz(int nq, const double qx[], const double qy[], const double qz[], double Iq[]) {
  _looper = new LoopQxyz(nq, qx, qy, qz, Iq);
  loop_Iq();
  delete _looper;
}

double
Disperser::calc_VR() {
  _volume = 1;
  _Vnorm = (_model.vloops == 0 ? 1. : 0.);
  _target = VR;
  loop_par(0, 1.);
  return _Vnorm;
}

double
Disperser::calc_ER() {
  _volume = 1;
  _Vnorm = (_model.vloops == 0 ? 1. : 0.);
  _target = ER;
  loop_par(0, 1.);
  return _Vnorm;
}

void
Disperser::loop_Iq(void) {
  _volume = 1; // in case there is no volume normalization
  _Vnorm = (_model.vloops == 0 ? 1. : 0.);
  _looper->clear();
  _target = IQ;
  loop_par(0, 1.);
  _looper->scale(_Vnorm/_Wnorm);
}

void 
Disperser::loop_par(int loop, double weight)
{
  int p = _model.looporder[loop];
  int n = (p == 0 ? _endpts[0] : _endpts[p] - _endpts[p-1]);
  int offset = (p == 0 ? 0 : _endpts[p-1]);
  const Weight *w = _pars + offset;
  for (int i=0; i < n; i++) {
    double wi = w[i].weight * weight;
    _dp[p] = w[i].value;
    if (loop == _model.vloops-1) {
      if (_target == VR) {
        _volume = wi * _model.formV(&_dp[0]);
        _Vnorm += _volume;
      } else if (_target == ER) {
        _volume = wi * _model.formER(&_dp[0]);
        _Vnorm += _volume;
      }
    } else if (loop == _model.vloops) {
        if (_target == ER || _target == VR) break;
    } 
    // Go the the next loop level
    if (loop < _model.npars) {
      loop_par(loop+1, wi);
    } else {
      _Wnorm += wi;
      _looper->loop(_model, &_dp[0], wi*_volume);
    }
  }
}
