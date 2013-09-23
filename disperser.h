#include <vector>
#include <ModelInfo.h>


class DisperserModel {
public:
    unsigned int vloops;
    unsigned int npars;
    std::vector<int> looporder;
    DisperserModel(const ModelInfo &m);

    virtual double 
    formV(const double dp[]) const = 0;

    virtual double
    formQ(const double dp[], double q) const = 0;

    // Default implementation of Qxy is Q = sqrt(x^2+y^2)
    virtual double
    formQxy(const double dp[], double qx, double qy) const;

    // Default implementation of Qxyz is Qxy ignoring z.
    virtual double
    formQxyz(const double dp[], double qx, double qy, double qz) const;

    virtual double
    formER(const double dp[]) const = 0;

    virtual double
    formVR(const double dp[]) const = 0;
} ;

class QLooper {
public:
    virtual void clear() = 0;
    virtual void scale(double weight) = 0;
    virtual void loop(const DisperserModel& m, const double dp[], double weight) = 0;
} ;

class Disperser {

private:
    enum CalcTarget { IQ, ER, VR } ;
    // The calculator for the model at Q.  This also contains the parameter
    // looping order and the number of volume parameters so that volume is
    // only calculated once.
    const DisperserModel &_model;

    // The input parameters for the dispersion calculation
    const int *_endpts;
    const Weight *_pars;

    // The point at which the model is being calculated.  The loop_par
    // function fills in the elements of the vector as it traverses the loop.
    std::vector<double> _dp;
    
    // The Q vector to calculate.  These are stored in the disperser object
    // rather than passed through the recursion over each looping parameter.

    QLooper *_looper;
    
    // The normalization factors for the dispersion.  These are accumulated
    // as we loop over the dispersion parameters.  The volume and volume
    // normalization are only calculated when looping over the final volume
    // calculator.
    double _volume;
    double _Vnorm;
    double _Wnorm;
    CalcTarget _target;
    
public:
    Disperser(const DisperserModel &model, int endpts[], Weight pars[])
        : _model(model), _endpts(endpts), _pars(pars) , _dp(model.npars)
        { } 
    void calc_Q(int nq, const double q[], double Iq[]);
    void calc_Qxy(int nq, const double qx[], const double qy[], double Iq[]);
    void calc_Qxyz(int nq, const double qx[], const double qy[], const double qz[], double Iq[]);
    double calc_ER();
    double calc_VR();

private:
    void loop_Iq(void);
    void loop_VR(void);
    void loop_ER(void);
    void loop_par(unsigned int loop, double weight);

} ;
