#pragma once
#define EXPORT __declspec(dllexport)

#include "kratos_external_bindings.h"

using namespace std;

extern "C"
{

#if defined(KRATOS_COMPILED_IN_WINDOWS)

    EXPORT Kratos::KratosExecute *KratosExecute_CreateInstance()
    {
        return new Kratos::KratosExecute();
    }

    EXPORT int __stdcall DGeoFlow(Kratos::KratosExecute *instance, char *workingDirectory, char *projectFile,
                                  double minCriticalHead, double maxCriticalHead, double stepCriticalHead, char *criticalHeadBoundaryModelPartName)
    {
    	return instance->geoflow(workingDirectory, projectFile, minCriticalHead, maxCriticalHead, stepCriticalHead, criticalHeadBoundaryModelPartName);
    }

#endif
}
