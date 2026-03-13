// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_workflow_factory.h"
#include "dgeoflow.h"
#include "dgeosettlement.hpp"
#include "time_loop_executor.h"

#include "custom_utilities/file_input_utility.h"
#include "custom_utilities/json_process_info_parser.h"

#include <memory>

namespace Kratos
{

KratosGeoSettlement* CustomWorkflowFactory::CreateKratosGeoSettlement()
{
    return new KratosGeoSettlement{std::make_unique<FileInputUtility>(),
                                   std::make_unique<JsonProcessInfoParser>(),
                                   std::make_unique<TimeLoopExecutor>()};
}

KratosExecute* CustomWorkflowFactory::CreateKratosGeoFlow() { return new KratosExecute{}; }

} // namespace Kratos
