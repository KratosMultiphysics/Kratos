// KRATOS___
//     //   ) )
//    //         ___      ___
//  //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "generic.h"

#include <utility>

#include "custom_workflows/dgeosettlement.hpp"
#include "custom_utilities/input_utility.h"
#include "custom_utilities/process_info_parser.h"
#include "time_loop_executor_interface.h"

namespace Kratos
{

KratosGeneric::KratosGeneric(std::unique_ptr<InputUtility> pInputUtility,
                             std::unique_ptr<ProcessInfoParser> pProcessInfoParser,
                             std::unique_ptr<TimeLoopExecutorInterface> pTimeLoopExecutor)
    : mpWorkflow{std::make_unique<KratosGeoSettlement>(
          std::move(pInputUtility), std::move(pProcessInfoParser), std::move(pTimeLoopExecutor))}
{
}

KratosGeneric::~KratosGeneric() = default;

int KratosGeneric::Run(const std::filesystem::path& rWorkingDirectory,
                       const std::filesystem::path& rProjectParametersFile,
                       const std::function<void(const char*)>& rLogCallback,
                       const std::function<void(double)>& rProgressDelegate,
                       const std::function<void(const char*)>& rReportTextualProgress,
                       const std::function<bool()>& rShouldCancel)
{
    return mpWorkflow->RunStage(rWorkingDirectory, rProjectParametersFile, rLogCallback,
                                rProgressDelegate, rReportTextualProgress, rShouldCancel);
}

} // namespace Kratos
