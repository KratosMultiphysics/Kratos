#include "kratos_external_bindings.h"
#include "custom_workflows/custom_workflow_factory.h"

#include "custom_workflows/dgeoflow.h"
#include "custom_workflows/dgeosettlement.h"


namespace Kratos
{

KratosExecute* ExternalBindings::KratosExecute_CreateInstance() { return new KratosExecute(); }

int ExternalBindings::execute_flow_analysis(KratosExecute* instance,
                          const char*            workingDirectory,
                          const char*            projectFile,
                          double                 minCriticalHead,
                          double                 maxCriticalHead,
                          double                 stepCriticalHead,
                          const char*            criticalHeadBoundaryModelPartName,
                          void                   logCallback(const char*),
                          void                   reportProgress(double),
                          void                   reportTextualProgress(const char*),
                          bool                   shouldCancel())
{
    const KratosExecute::CriticalHeadInfo critical_head_info(
        minCriticalHead, maxCriticalHead, stepCriticalHead);
    const KratosExecute::CallBackFunctions call_back_functions(
        logCallback, reportProgress, reportTextualProgress, shouldCancel);

    return instance->ExecuteFlowAnalysis(workingDirectory, projectFile, critical_head_info,
                                         criticalHeadBoundaryModelPartName, call_back_functions);
}

KratosGeoSettlement* ExternalBindings::KratosGeoSettlement_CreateInstance()
{
    return CustomWorkflowFactory::CreateKratosGeoSettlement();
}

int ExternalBindings::runSettlementStage(KratosGeoSettlement* instance,
                       const char*                  workingDirectory,
                       const char*                  projectFileName,
                       void                         logCallback(const char*),
                       void                         reportProgress(double),
                       void                         reportTextualProgress(const char*),
                       bool                         shouldCancel())
{
    return instance->RunStage(workingDirectory, projectFileName, logCallback, reportProgress,
                              reportTextualProgress, shouldCancel);
}
}