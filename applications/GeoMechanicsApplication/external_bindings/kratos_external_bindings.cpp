#include "kratos_external_bindings.h"
#include "custom_workflows/custom_workflow_factory.h"
#include "includes/define.h"

#if defined(KRATOS_COMPILED_IN_WINDOWS)
#define EXPORT __declspec(dllexport)
#define CALL_CONV __stdcall
#else
#define EXPORT __attribute__((visibility("default")))
#define CALL_CONV
#endif

extern "C" {

EXPORT Kratos::KratosExecute* KratosExecute_CreateInstance() { return new Kratos::KratosExecute(); }

EXPORT int CALL_CONV execute_flow_analysis(Kratos::KratosExecute* instance,
                                           const char*            workingDirectory,
                                           const char*            projectFile,
                                           double                 minCriticalHead,
                                           double                 maxCriticalHead,
                                           double                 stepCriticalHead,
                                           const char*            criticalHeadBoundaryModelPartName,
                                           void(CALL_CONV* logCallback)(const char*),
                                           void(CALL_CONV* reportProgress)(double),
                                           void(CALL_CONV* reportTextualProgress)(const char*),
                                           bool(CALL_CONV* shouldCancel)())
{
    const Kratos::KratosExecute::CriticalHeadInfo critical_head_info(
        minCriticalHead, maxCriticalHead, stepCriticalHead);
    const Kratos::KratosExecute::CallBackFunctions call_back_functions(
        logCallback, reportProgress, reportTextualProgress,
        shouldCancel ? shouldCancel : []() { return false; });

    return instance->ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info,
                                         criticalHeadBoundaryModelPartName, call_back_functions);
}

EXPORT Kratos::KratosGeoSettlement* KratosGeoSettlement_CreateInstance()
{
    return Kratos::CustomWorkflowFactory::CreateKratosGeoSettlement();
}

EXPORT int CALL_CONV runSettlementStage(Kratos::KratosGeoSettlement* instance,
                                        const char*                  workingDirectory,
                                        const char*                  projectFileName,
                                        void(CALL_CONV* logCallback)(const char*),
                                        void(CALL_CONV* reportProgress)(double),
                                        void(CALL_CONV* reportTextualProgress)(const char*),
                                        bool(CALL_CONV* shouldCancel)())
{
    return instance->RunStage(workingDirectory, projectFileName, logCallback, reportProgress,
                              reportTextualProgress, shouldCancel);
}
}
