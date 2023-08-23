#define EXPORT __declspec(dllexport)
#include "kratos_external_bindings.h"

extern "C"
{

#if defined(KRATOS_COMPILED_IN_WINDOWS)

    EXPORT Kratos::KratosExecute* KratosExecute_CreateInstance()
    {
        return new Kratos::KratosExecute();
    }

    EXPORT int __stdcall execute_flow_analysis(Kratos::KratosExecute* instance,
                                               const char*            workingDirectory,
                                               const char*            projectFile,
                                               double                 minCriticalHead,
                                               double                 maxCriticalHead,
                                               double                 stepCriticalHead,
                                               const char*            criticalHeadBoundaryModelPartName,
                                               void __stdcall         logCallback(const char*),
                                               void __stdcall         reportProgress(double),
                                               void __stdcall         reportTextualProgress(const char*),
                                               bool __stdcall         shouldCancel())
    {
        int errorCode = instance->ExecuteFlowAnalysis(workingDirectory,
                                                      projectFile,
                                                      minCriticalHead,
                                                      maxCriticalHead,
                                                      stepCriticalHead,
                                                      criticalHeadBoundaryModelPartName,
                                                      logCallback,
                                                      reportProgress,
                                                      reportTextualProgress,
                                                      shouldCancel);
        return errorCode;
    }



    EXPORT Kratos::KratosGeoSettlement* KratosGeoSettlement_CreateInstance()
    {
        return new Kratos::KratosGeoSettlement{};
    }

    EXPORT int __stdcall runSettlementStage(Kratos::KratosGeoSettlement* instance,
                                            const char*                  workingDirectory,
                                            const char*                  projectFileName,
                                            void __stdcall               logCallback(const char*),
                                            void __stdcall               reportProgress(double),
                                            void __stdcall               reportTextualProgress(const char*),
                                            bool __stdcall               shouldCancel())
    {
        return instance->RunStage(workingDirectory,
                                  projectFileName,
                                  logCallback,
                                  reportProgress,
                                  reportTextualProgress,
                                  shouldCancel);
    }

#endif
}
