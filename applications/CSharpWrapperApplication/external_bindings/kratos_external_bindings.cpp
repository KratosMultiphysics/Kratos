#define EXPORT __declspec(dllexport)

#include <stdlib.h>
#include "custom_includes/kratos_wrapper.h"
#include <iostream>

using namespace std;
using namespace CSharpKratosWrapper;

extern "C" {

#if defined(KRATOS_COMPILED_IN_WINDOWS)
EXPORT KratosWrapper *Kratos_CreateInstance() {
    return new KratosWrapper();
}

EXPORT void __stdcall Kratos_Init(KratosWrapper *instance, char *mdpaPath, char *parametersJsonPath) {
    instance->init(mdpaPath, parametersJsonPath);
}

EXPORT void __stdcall Kratos_InitWithMDPA(KratosWrapper *instance, char *mdpaPath) {
    instance->init(mdpaPath);
}

EXPORT void __stdcall Kratos_InitWithSettings(KratosWrapper *instance, char *parametersJsonPath) {
    instance->initWithSettings(parametersJsonPath);
}

EXPORT ModelPartWrapper *__stdcall Kratos_GetRootModelPart(KratosWrapper *instance) {
    return instance->getRootModelPartWrapper();
}

EXPORT void __stdcall Kratos_Calculate(KratosWrapper *instance) {
    instance->calculate();
}

EXPORT void __stdcall Kratos_DisposeKratos(KratosWrapper *instance) {
    delete instance;
}

#endif
}
