// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf

#include "custom_processes/apply_seepage_boundary_process.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

using namespace std::string_literals;

namespace Kratos
{

ApplySeepageBoundaryProcess::ApplySeepageBoundaryProcess(Model& rModel, const Parameters& rProcessSettings)
{
    KRATOS_ERROR_IF_NOT(rProcessSettings.Has("model_part_name"))
        << "ApplySeepageBoundaryProcess: \"model_part_name\" is required in the process settings"
        << std::endl;
}

std::string ApplySeepageBoundaryProcess::Info() const { return "ApplySeepageBoundaryProcess"s; }

} // namespace Kratos