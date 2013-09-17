#include <cstring>
#include "ModelInfo.h"
#include "sphere.h"

ParameterInfo param_infos[] = {
    { "scale"       , NULL                      , "Scale factor", 1.0   , -DBL_INF, +DBL_INF, PF_None},
    { "radius"      , "radius of sphere"        , "A"           , 60.0  , -DBL_INF, +DBL_INF, PF_None},
    { "sldSph"      , "the SLD of the sphere"   , "1/A^(2)"     , 2.0e-6, -DBL_INF, +DBL_INF, PF_None},
    { "sldSolv"     , "the SLD of the solvent"  , "1/A^(2)"     , 1.0e-6, -DBL_INF, +DBL_INF, PF_None},
    { "background"  , NULL                      , "1/cm"        , 0.0   , -DBL_INF, +DBL_INF, PF_None},
    { "M0_sld_sph"  , NULL                      , "1/A^(2)"     , 0.0e-6, -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "M_theta_sph" , NULL                      , "deg"         , 0.0   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "M_phi_sph"   , NULL                      , "deg"         , 0.0   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "M0_sld_solv" , NULL                      , "1/A^(2)"     , 0.0e-6, -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "M_theta_solv", NULL                      , "deg"         , 0.0   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "M_phi_solv"  , NULL                      , "deg"         , 0.0   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "Up_frac_i"   , NULL                      , "u/(u+d)"     , 0.5   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "Up_frac_f"   , NULL                      , "u/(u+d)"     , 0.5   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
    { "Up_theta"    , NULL                      , "deg"         , 0.0   , -DBL_INF, +DBL_INF, PF_Orientation | PF_Magnetic},
};

ModelInfo model_info(
    "SphereModel",
    "P(q)=(scale/V)*[3V(sldSph-sldSolv)*(sin(qR)-qRcos(qR))/(qR)^3]^(2)+bkg",
    GetParameterCount(param_infos),
    param_infos);

// helper
bool update_model(void* ptr, void** p) {
std::cout << "update model " << ptr << " " << p << std::endl;

    if ((ptr == NULL) || (p == NULL))
        return false;

    SphereModel* model        = (SphereModel*)ptr;
    Parameter*   model_params = &(model->scale);

    // update model
    ParameterInfo* param_info = param_infos;
std::cout << "npar " <<  model_info.ParameterCount << std::endl;
    for (unsigned int i=0; i < model_info.ParameterCount; i++) {
std::cout << i << " " <<  p << " " << *p << std::endl;
        if (p == NULL)
            return false;

std::cout << "type " << *(size_t*)*p << std::endl;
        switch (*(size_t*)*p) {
        case PT_End:
std::cout << i << " end" << std::endl;
            return false;
        case PT_Simple:
std::cout << i << " mono" << std::endl;
            model_params->value = AsSimpleParameter(*p);
            break;
        case PT_Polydisperse:
std::cout << i << " poly" << std::endl;
            model_params->dispersion = &AsPolydisperseParameter(*p);
            break;
        default:
std::cout << "unknown type " << *(size_t*)*p << std::endl;
        }

        model_params++;
        param_info++;
        p++;
    }
    return true;
}
void show_error(const char* function, int nq = 0, double iq[] = NULL) {
    std::cerr << "error in " << function << " - update_model failed!" << std::endl;
    while (nq-- != 0)
        *iq = DBL_NAN;
}

// model handling
CExport ModelInfo* get_model_info() {
    // ensure that returned pointer points to static data which doesn’t require to be released
    return &model_info;
}
CExport void* create_model(void* data) { // "data" is only provided if associated python script defines it
    SphereModel* model = new SphereModel();
    Parameter*   model_params = &(model->scale);

    // update description after construction of model
    ParameterInfo* param_info = param_infos;
    for (unsigned int i=0; i < model_info.ParameterCount; i++) {
        param_info->DispMin = model_params->min;
        param_info->DispMax = model_params->max;
        
        if (model_params->has_dispersion)
            param_info->Flags |= PF_Polydisperse;
        
        model_params++;
        param_info++;
    }
    
    return model;
}
CExport void destroy_model(void* ptr) {
    delete (SphereModel*)ptr;
}
// calculations
CExport void calculate_q(void* ptr, void** p, int nq, double iq[], double q[]) {
    if (!update_model(ptr, p))
        show_error("calculate_q", nq, iq);
    else
{
std::cout << "err q" << " " << nq << " " << q << " " << iq << std::endl;
std::cout << "err q" << " " << nq << " " << *(q) << " " << *(iq) << std::endl; 
        while (nq-- != 0)
{
            *iq++ = (*(SphereModel*)ptr)(*q++);
std::cout << "err q" << " " << nq << " " << q << " " << iq << std::endl;
std::cout << "err q" << " " << nq << " " << q[-1] << " " << iq[-1] << std::endl; }
std::cout << "done q" << std::endl; }
}
CExport void calculate_qxqy(void* ptr, void** p, int nq, double iq[], double qx[], double qy[]) {
    if (!update_model(ptr, p))
        show_error("calculate_qxqy", nq, iq);
    else
        while (nq-- != 0)
            *iq++ = (*(SphereModel*)ptr)(*qx++, *qy++);
}
CExport void calculate_qxqyqz(void* ptr, void** p, int nq, double iq[], double qx[], double qy[], double qz[]) {
    if (!update_model(ptr, p))
        show_error("calculate_qxqyqz", nq, iq);
    else
        memset(iq, 0, sizeof(double) * nq);
}
CExport double calculate_ER(void* ptr, void** p) {
    if (!update_model(ptr, p)) {
        show_error("calculate_ER");
        return DBL_NAN;
    }
    return ((SphereModel*)ptr)->calculate_ER();
}
CExport double calculate_VR(void* ptr, void** p) {
    if (!update_model(ptr, p)) {
        show_error("calculate_VR");
        return DBL_NAN;
    }
    return ((SphereModel*)ptr)->calculate_VR();
}
