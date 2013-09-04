#include "ModelInfo.h"

SimpleParameter* GetSimpleParameter(void* p) {
    if (*(int*)p == PT_Simple)
        return (SimpleParameter*)p;
    else
        return NULL;
}
PolydisperseParameter* GetPolydisperseParameter(void* p) {
    if (*(int*)p == PT_Polydisperse)
        return (PolydisperseParameter*)p;
    else
        return NULL;
}
