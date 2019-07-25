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

EXPORT ModelPartWrapper *__stdcall GetRootModelPart(KratosWrapper *instance) {
    return instance->getRootModelPartWrapper();
}

EXPORT void __stdcall Calculate(KratosWrapper *instance) {
    instance->calculate();
}

EXPORT void __stdcall DisposeKratos(KratosWrapper *instance) {
    delete instance;
}

#endif
}
