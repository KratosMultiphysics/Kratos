#define EXPORT __declspec(dllexport)

#include <stdlib.h>
#include "kratos_wrapper.h"
#include <iostream>

using namespace std;
using namespace CSharpKratosWrapper;

extern "C" {

#if defined(KRATOS_COMPILED_IN_WINDOWS)
EXPORT KratosWrapper *CreateInstance() {
    return new KratosWrapper();
}

EXPORT void __stdcall Init(KratosWrapper *instance, char *mdpaPath, char *parametersJsonPath) {
    instance->init(mdpaPath, parametersJsonPath);
}

EXPORT void __stdcall InitWithMDPA(KratosWrapper *instance, char *mdpaPath) {
    instance->init(mdpaPath);
}

EXPORT void __stdcall InitWithSettings(KratosWrapper *instance, char *parametersJsonPath) {
    instance->initWithSettings(parametersJsonPath);
}

EXPORT float *__stdcall GetXCoordinates(KratosWrapper *instance) {
    return instance->getXCoordinates();
}

EXPORT float *__stdcall GetYCoordinates(KratosWrapper *instance) {
    return instance->getYCoordinates();
}

EXPORT float *__stdcall GetZCoordinates(KratosWrapper *instance) {
    return instance->getZCoordinates();
}

EXPORT int __stdcall GetNodesCount(KratosWrapper *instance) {
    return instance->getNodesCount();
}

EXPORT int *__stdcall GetTriangles(KratosWrapper *instance) {
    return instance->getTriangles();
}

EXPORT int __stdcall GetTrianglesCount(KratosWrapper *instance) {
    return instance->getTrianglesCount();
}

EXPORT void __stdcall UpdateNodePos(KratosWrapper *instance, int nodeId, float x, float y, float z) {
    instance->updateNodePos(nodeId, x, y, z);
}

EXPORT void __stdcall Calculate(KratosWrapper *instance) {
    instance->calculate();
}

EXPORT void __stdcall DisposeInstance(KratosWrapper *instance) {
    free(instance);
}

#endif
}
