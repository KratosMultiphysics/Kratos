#pragma once
#define EXPORT __declspec(dllexport)

#include "kratos_external_bindings.h"

using namespace std;

extern "C" {

#if defined(KRATOS_COMPILED_IN_WINDOWS)

EXPORT Kratos::KratosExecute *KratosExecute_CreateInstance() {
    return new Kratos::KratosExecute();
}

EXPORT int __stdcall Execute(Kratos::KratosExecute* instance, char* workingDirectory, char* projectFile)
    {
    return instance->cpp_geomechanics(workingDirectory, projectFile);
}

#endif
}
