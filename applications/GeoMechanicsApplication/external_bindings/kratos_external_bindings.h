#include "includes/kratos_export_api.h"

namespace Kratos
{

class KratosGeoSettlement;
class KratosExecute;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ExternalBindings
{
    KratosExecute* KratosExecute_CreateInstance();

    int execute_flow_analysis(KratosExecute* instance,
                              const char*    workingDirectory,
                              const char*    projectFile,
                              double         minCriticalHead,
                              double         maxCriticalHead,
                              double         stepCriticalHead,
                              const char*    criticalHeadBoundaryModelPartName,
                              void __stdcall logCallback(const char*),
                              void __stdcall reportProgress(double),
                              void __stdcall reportTextualProgress(const char*),
                              bool __stdcall shouldCancel());

    KratosGeoSettlement* KratosGeoSettlement_CreateInstance();

    int runSettlementStage(KratosGeoSettlement* instance,
                           const char*          workingDirectory,
                           const char*          projectFileName,
                           void __stdcall logCallback(const char*),
                           void __stdcall reportProgress(double),
                           void __stdcall reportTextualProgress(const char*),
                           bool __stdcall shouldCancel());
};

} // namespace Kratos