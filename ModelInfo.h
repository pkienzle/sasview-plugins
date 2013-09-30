#ifndef MODELINFO_H
#define MODELINFO_H

#ifdef _MSC_VER
#define CExport extern "C" __declspec(dllexport)
#else
#define CExport extern "C"
#endif

#include <stdint.h>
#include <iostream>
#include <limits>

enum ParameterFlags {
    PF_None         = 0x00,
    PF_Orientation  = 0x01,
    PF_Magnetic     = 0x02,
    PF_Unfittable   = 0x04,
    PF_Integer      = 0x08,
    PF_Polydisperse = 0x10,
    PF_RepeatCount  = 0x20 | PF_Unfittable, // for multiplicity
    PF_Repeated     = 0x40,
    PF_Volume       = 0x80,
};

struct ParameterInfo {
    const char*   Name;
    const char*   Description;
    const char*   Unit;
    double  Default;
    double  DispMin;
    double  DispMax;
    size_t  Flags;
};

#define GetParameterCount(parameters) (sizeof(parameters) / sizeof(ParameterInfo))

struct ModelInfo {
    size_t          Version;
    const char*           ModelName;
    const char*           ModelDescription;
    size_t          ParameterCount;
    const ParameterInfo*  Parameters;
    
    ModelInfo(const char* modelName, const char* modelDescription, size_t parameterCount, const ParameterInfo* parameters) :
        Version(1),
        ModelName(modelName),
        ModelDescription(modelDescription),
        ParameterCount(parameterCount),
        Parameters(parameters) {
    }
};

struct Weight {
    double value, weight;
} ;


CExport void* get_model_info();
CExport void* create_model(void* data);
CExport void destroy_model(void* ptr);
CExport void calculate_q(void* ptr, const int pindex[], const Weight p[], size_t nq, double iq[], const double q[]);
CExport void calculate_qxqy(void* ptr, const int pindex[], const Weight p[], size_t nq, double iq[], const double qx[], const double qy[]);
CExport void calculate_qxqyqz(void* ptr, int pindex[], Weight p[], size_t nq, double iq[], double qx[], double qy[], double qz[]);
CExport double calculate_ER(void* ptr, int pindex[], Weight p[]);
CExport double calculate_VR(void* ptr, int pindex[], Weight p[]);

#ifndef DBL_INF
#define DBL_INF std::numeric_limits<double>::infinity()
#endif
#ifndef DBL_NAN
#define DBL_NAN std::numeric_limits<double>::quiet_NaN()
#endif

#endif // MODELINFO_H
