#include <vector>
#include <ModelInfo.h>


class DisperserModel {
public:
    // The order that the parameters are processed in the dispersion loops may
    // be different from the order in which they appear in the model.
    // (1) Volume parameters are processed first.  Once the final volume
    //     parameter has been set, the volume can be computed and the weighted
    //     volume accumulated.  If for example you are dispersing over the
    //     density, you do not need to recalculate the volume for each density.
    // (2) Untagged parameters are processed next.  Calculating I(q) requires
    //     all volume parameters and untagged parameters, but no orientation or
    //     magnetism parameters.
    // (3) Orientation parameters are processed next. Magnetism measurements
    //     are inherently oriented (e.g., field direction on a sphere), but
    //     orientation parameters are not inherently magnetic.
    // (4) Magnetism parameters are last, since they are only needed when
    //     calculating theory curves for magnetic data.
    //
    // Spin state is a property of the model for now, rather than the data.
    // To calculate all four cross sections for spin in and spin out the entire
    // remainder of the model must also be recalculated.  It may be better to
    // have a calc magnetic method which takes a list of the desired cross sections
    // instead, and returns them in a magnetic cross section Iq_magnetic vector with
    // room for all four.
    //
    // Beam parameters (incident intensity and background) are also part of the
    // model, with each model having its own beam parameters.  It may be better to
    // move these out of the models to slightly simplify model development.
    std::vector<int> looporder; // The order the parameters are processed
    unsigned int volume_loops;  // First non-volume parameter
    unsigned int orientation_loops; // First orientation parameter (>volume_loops)
    unsigned int magnetic_loops; // First magnetism parameter (>orienatation_loops)
    unsigned int npars;
    DisperserModel(const ModelInfo &m);

    // Compute the volume of the model using only the parameters marked as PF_VOLUME.
    virtual double 
    formV(const double dp[]) const = 0;

    // Compute I(q) without using orientation or magnetic parameters.
    virtual double
    formQ(const double dp[], double q) const = 0;

    // Compute I(qx,qy) with orientation and magnetism.
    // Default implementation of Qxy is Q = sqrt(qx^2+qy^2).
    // This function is not needed when fitting rotationally symmetric
    // patterns, so don't worry that it is inefficient to compute the
    // same q for different qx/qy pairs.
    virtual double
    formQxy(const double dp[], double qx, double qy) const;

    // Compute I(qx,qy,qz) with orientation and magnetism.
    // Default implementation of Qxyz is Qxy ignoring z.
    virtual double
    formQxyz(const double dp[], double qx, double qy, double qz) const;

    // Compute the effective radius of the shape by averaging over the volume.
    // Note that this does not take into account density or contrast, so the
    // result may be misleading --- a constrast matched sphere will be invisible
    // but still have an effective radius.
    virtual double
    formER(const double dp[]) const = 0;

    virtual double
    formVR(const double dp[]) const = 0;
} ;

class Disperser {

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
    int _nq;
    double *_iq;
    const double *_qx, *_qy, *_qz;

    // The normalization factors for the dispersion.  These are accumulated
    // as we loop over the dispersion parameters.  The volume and volume
    // normalization are only calculated when looping over the final volume
    // calculator.
    double _volume;
    double _Vnorm;
    double _Wnorm;
    
    // Loop over the parameters, calculating one of the possible targets
    // in the inner loop.
    enum LoopTarget { IQ, IQXY, IQXYZ, ER, VR } ;
    LoopTarget _target;
    void loop_par(unsigned int loop, double weight);

    // Start the loop for I(q), I(qx,qy) or I(qx,qy,qz), and normalize the result when complete.
    void loop_Iq(void);
    // Start the loop for VR(), normalizing the result when complete.
    void loop_VR(void);
    // Start the loop for ER(), normalizing the result when complete.
    void loop_ER(void);

} ;
