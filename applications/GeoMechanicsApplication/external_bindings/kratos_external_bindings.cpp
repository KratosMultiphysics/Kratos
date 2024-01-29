#define EXPORT __declspec(dllexport)
#include "kratos_external_bindings.h"
#include "custom_workflows/custom_workflow_factory.h"

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
        const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(minCriticalHead, maxCriticalHead, stepCriticalHead);
        const Kratos::KratosExecute::CallBackFunctions call_back_functions(logCallback,
                                                                     reportProgress,
                                                                     reportTextualProgress,
                                                                     shouldCancel);

        return instance->ExecuteFlowAnalysis(workingDirectory,
                                                      projectFile,
                                                      critical_head_info,
                                                      criticalHeadBoundaryModelPartName,
                                                      call_back_functions);
    }

    EXPORT Kratos::KratosGeoSettlement* KratosGeoSettlement_CreateInstance()
    {
        return Kratos::CustomWorkflowFactory::CreateKratosGeoSettlement();
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
