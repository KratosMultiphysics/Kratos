// KRATOS___
//     //   ) )
//    //         ___      ___
//  //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <filesystem>
#include <functional>
#include <memory>

#include "includes/kratos_export_api.h"

namespace Kratos
{

class InputUtility;
class ProcessInfoParser;
class TimeLoopExecutorInterface;
class KratosGeoSettlement;

class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeneric
{
public:
    KratosGeneric(std::unique_ptr<InputUtility> pInputUtility,
                  std::unique_ptr<ProcessInfoParser> pProcessInfoParser,
                  std::unique_ptr<TimeLoopExecutorInterface> pTimeLoopExecutor);

    ~KratosGeneric();

    int Run(const std::filesystem::path& rWorkingDirectory,
            const std::filesystem::path& rProjectParametersFile,
            const std::function<void(const char*)>& rLogCallback,
            const std::function<void(double)>& rProgressDelegate,
            const std::function<void(const char*)>& rReportTextualProgress,
            const std::function<bool()>& rShouldCancel);

private:
    std::unique_ptr<KratosGeoSettlement> mpWorkflow;
};

} // namespace Kratos
