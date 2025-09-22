#include "includes/kratos_export_api.h"

namespace Kratos
{

class KratosGeoSettlement;
class KratosExecute;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ExternalBindings
{
    static KratosExecute* KratosExecute_CreateInstance();

    static int execute_flow_analysis(KratosExecute* instance,
                                     const char*    workingDirectory,
                                     const char*    projectFile,
                                     double         minCriticalHead,
                                     double         maxCriticalHead,
                                     double         stepCriticalHead,
                                     const char*    criticalHeadBoundaryModelPartName,
                                     void           logCallback(const char*),
                                     void           reportProgress(double),
                                     void           reportTextualProgress(const char*),
                                     bool           shouldCancel());

    static KratosGeoSettlement* KratosGeoSettlement_CreateInstance();

    static int runSettlementStage(KratosGeoSettlement* instance,
                                  const char*          workingDirectory,
                                  const char*          projectFileName,
                                  void                 logCallback(const char*),
                                  void                 reportProgress(double),
                                  void                 reportTextualProgress(const char*),
                                  bool                 shouldCancel());
};

} // namespace Kratos