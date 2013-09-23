#include <cstdio>
#include <cmath>

#include <disperser.h>

// Model definition
extern "C" {
	#include <libCylinder.h>
	#include <libStructureFactor.h>
	#include <libmultifunc/libfunc.h>
}

typedef struct {
  double scale;
  double radius;
  double length;
  double sldCyl;
  double sldSolv;
  double background;
  double cyl_theta;
  double cyl_phi;
  double M0_sld_cyl;
  double M_theta_cyl;
  double M_phi_cyl;
  double M0_sld_solv;
  double M_theta_solv;
  double M_phi_solv;
  double Up_frac_i;
  double Up_frac_f;
  double Up_theta;
} Parameters;

ParameterInfo param_infos[] = {
   // name, description, units, default, min, max, flags
   { "scale", "Scale factor", "", 1.0, 0, +DBL_INF, PF_None },
   { "radius", "Radius of the cylinder", "A", 20.0, 0, +DBL_INF, PF_Polydisperse|PF_Volume },
   { "length", "Length of the cylinder", "A", 400.0, 0, +DBL_INF, PF_Polydisperse|PF_Volume },
   { "sldCyl", "Contrast", "1/A^(2)", 4e-06, 0, +DBL_INF, PF_None },
   { "sldSolv", "sldCyl", "1/A^(2)", 1e-06, 0, +DBL_INF, PF_None },
   { "background", "Incoherent Background [1/cm] 0.00", "1/cm", 0.0, 0, +DBL_INF, PF_None },
   { "cyl_theta", "Orientation of the cylinder axis w/respect incoming beam", "deg", 60.0, 0, +DBL_INF, PF_Orientation|PF_Polydisperse },
   { "cyl_phi", "Orientation of the cylinder in the plane of the detector", "deg", 60.0, 0, +DBL_INF, PF_Orientation|PF_Polydisperse },
   { "M0_sld_cyl", "M0_sld_cyl", "1/A^(2)", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "M_theta_cyl", "M_theta_cyl", "deg", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "M_phi_cyl", "M_phi_cyl", "deg", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "M0_sld_solv", "M0_sld_solv", "1/A^(2)", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "M_theta_solv", "M_theta_solv", "deg", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "M_phi_solv", "M_phi_solv", "deg", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "Up_frac_i", "Up_frac_i", "u/(u+d)", 0.5, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "Up_frac_f", "Up_frac_f", "u/(u+d)", 0.5, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
   { "Up_theta", "Up_theta", "deg", 0.0, 0, +DBL_INF, PF_Magnetic|PF_Orientation },
};

ModelInfo model_info(
    "Cylinder",
    "\
f(q)= 2*(sldCyl - sldSolv)*V*sin(qLcos(alpha/2))\n\
                /[qLcos(alpha/2)]*J1(qRsin(alpha/2))/[qRsin(alpha)]\n\
\n\
P(q,alpha)= scale/V*f(q)^(2)+bkg\n\
V: Volume of the cylinder\n\
R: Radius of the cylinder\n\
L: Length of the cylinder\n\
J1: The bessel function\n\
alpha: angle betweenthe axis of the\n\
cylinder and the q-vector for 1D\n\
:the ouput is P(q)=scale/V*integral\n\
from pi/2 to zero of...\n\
f(q)^(2)*sin(alpha)*dalpha+ bkg\
",
    GetParameterCount(param_infos),
    param_infos);


/**
 * Function to evaluate 2D scattering function
 * @param pars: parameters of the cylinder
 * @param q: q-value
 * @param q_x: q_x / q
 * @param q_y: q_y / q
 * @return: function value
 */
static double cylinder_analytical_2D_scaled(const Parameters &pars, double q, double q_x, double q_y) {
	double cyl_x, cyl_y;//, cyl_z;
	//double q_z;
	double alpha, vol, cos_val;
	double answer = 0.0;
	double form = 0.0;
	//convert angle degree to radian
	double pi = 4.0*atan(1.0);
	double theta = pars.cyl_theta * pi/180.0;
	double phi = pars.cyl_phi * pi/180.0;
	double sld_solv = pars.sldSolv;
	double sld_cyl = pars.sldCyl;
	double m_max = pars.M0_sld_cyl;
	double m_max_solv = pars.M0_sld_solv;
	double contrast = 0.0;

	// Cylinder orientation
	cyl_x = cos(theta) * cos(phi);
	cyl_y = sin(theta);
	//cyl_z = -cos(theta) * sin(phi);
	// q vector
	//q_z = 0.0;

	// Compute the angle btw vector q and the
	// axis of the cylinder
	cos_val = cyl_x*q_x + cyl_y*q_y;// + cyl_z*q_z;

	// The following test should always pass
	if (fabs(cos_val)>1.0) {
	  std::printf("cyl_ana_2D: Unexpected error: |cos(alpha)=%g|>1\n", cos_val);
	  std::printf("cyl_ana_2D: at theta=%g and phi=%g.", theta, phi);
	  return 1.0;
	}

	// Note: cos(alpha) = 0 and 1 will get an
	// undefined value from CylKernel
  alpha = acos( cos_val );
  if (alpha == 0.0){
	alpha = 1.0e-26;
	}
  // Call the IGOR library function to get the kernel
  //answer = CylKernel(q, pars.radius, pars.length/2.0, alpha) / sin(alpha);

	// Call the IGOR library function to get the kernel
	form = CylKernel(q, pars.radius, pars.length/2.0, alpha) / sin(alpha);

	if (m_max < 1.0e-32 && m_max_solv < 1.0e-32){
		// Multiply by contrast^2
		contrast = (pars.sldCyl - pars.sldSolv);
		answer = contrast * contrast * form;
	}
	else{
		double qx = q_x;
		double qy = q_y;
		double s_theta = pars.Up_theta;
		double m_phi = pars.M_phi_cyl;
		double m_theta = pars.M_theta_cyl;
		double m_phi_solv = pars.M_phi_solv;
		double m_theta_solv = pars.M_theta_solv;
		double in_spin = pars.Up_frac_i;
		double out_spin = pars.Up_frac_f;
		polar_sld p_sld;
		polar_sld p_sld_solv;
		p_sld = cal_msld(1, qx, qy, sld_cyl, m_max, m_theta, m_phi, 
						in_spin, out_spin, s_theta);
		p_sld_solv = cal_msld(1, qx, qy, sld_solv, m_max_solv, m_theta_solv, m_phi_solv, 
						in_spin, out_spin, s_theta);
		//up_up 
		if (in_spin > 0.0 && out_spin > 0.0){			 
			answer += ((p_sld.uu- p_sld_solv.uu) * (p_sld.uu- p_sld_solv.uu) * form);
			}
		//down_down
		if (in_spin < 1.0 && out_spin < 1.0){
			answer += ((p_sld.dd - p_sld_solv.dd) * (p_sld.dd - p_sld_solv.dd) * form);
			}
		//up_down
		if (in_spin > 0.0 && out_spin < 1.0){
			answer += ((p_sld.re_ud - p_sld_solv.re_ud) * (p_sld.re_ud - p_sld_solv.re_ud) * form);
			answer += ((p_sld.im_ud - p_sld_solv.im_ud) * (p_sld.im_ud - p_sld_solv.im_ud) * form);
			}
		//down_up	
		if (in_spin < 1.0 && out_spin > 0.0){
			answer += ((p_sld.re_du - p_sld_solv.re_du) * (p_sld.re_du - p_sld_solv.re_du) * form);
			answer += ((p_sld.im_du - p_sld_solv.im_du) * (p_sld.im_du - p_sld_solv.im_du) * form);
			}
	}

	//normalize by cylinder volume
	//NOTE that for this (Fournet) definition of the integral, one must MULTIPLY by Vcyl
	vol = acos(-1.0) * pars.radius * pars.radius * pars.length;
	answer *= vol;

	//convert to [cm-1]
	answer *= 1.0e8;

	//Scale
	answer *= pars.scale;

	// add in the background
	answer += pars.background;

	return answer;
}

/**
 * Function to evaluate 2D scattering function
 * @param pars: parameters of the cylinder
 * @param q: q-value
 * @return: function value
 */
static double cylinder_analytical_2DXY(const Parameters& pars, double qx, double qy) {
  double q;
  q = sqrt(qx*qx+qy*qy);
  return cylinder_analytical_2D_scaled(pars, q, qx/q, qy/q);
}


class Model: public DisperserModel {

public:
    Model(ModelInfo &model_info) : DisperserModel(model_info) {}
	
    double 
    formV(const double dp[]) const {
        const Parameters &p = *(const Parameters*)dp;
        return p.radius * p.radius * p.length; // * M_PI;
    }

    double
    formQ(const double dp[], double q) const {
        //Parameters *p = (Parameters *)dp;
        return CylinderForm(const_cast<double*>(&dp[0]), q);
    }

    double
    formQxy(const double dp[], double qx, double qy) const {
        const Parameters &p = *(const Parameters *)dp;
        return cylinder_analytical_2DXY(p, qx, qy);
    }

    double
    formQxyz(const double dp[], double qx, double qy, double qz) const {
        // Parameters *p = (Parameters *)dp;
        return formQxy(dp, qx, qy);
    }

    double
    formER(const double dp[]) const {
        const Parameters &p = *(const Parameters *)dp;
        return DiamCyl(p.length,p.radius)/2.0;
    }

    double
    formVR(const double dp[]) const {
        //Parameters *p = (Parameters *)dp;
        return 1.0;
    }

} ;

// ====== Boiler-plate external interface definition ====

// TODO: possible minor memory leak if this isn't destroyed when the module is unloaded
static const Model model(model_info);
CExport void* get_model_info() { return &model_info; }
//CExport void* create_model(void* data) { return NULL; }
//CExport void destroy_model(void* ptr) {}
CExport void calculate_q(void* ptr, int pindex[], Weight pars[], size_t nq, double iq[], double q[]) {
	Disperser(model, pindex, pars).calc_Q(nq, q, iq);
}
CExport void calculate_qxqy(void* ptr, int pindex[], Weight pars[], size_t nq, double iq[], double qx[], double qy[]) {
	Disperser(model, pindex, pars).calc_Qxy(nq, qx, qy, iq);
}
CExport void calculate_qxqyqz(void* ptr, int pindex[], Weight pars[], size_t nq, double iq[], double qx[], double qy[], double qz[]) {
	Disperser(model, pindex, pars).calc_Qxyz(nq, qx, qy, qz, iq);
}
CExport double calculate_ER(void* ptr, int pindex[], Weight pars[]) {
	return Disperser(model, pindex, pars).calc_ER();
}
CExport double calculate_VR(void* ptr, int pindex[], Weight pars[]) {
	return Disperser(model, pindex, pars).calc_VR();
}

