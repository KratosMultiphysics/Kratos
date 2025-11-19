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
//

#include "find_neighbours_of_interfaces_process.h"
#include "containers/model.h"
#include "custom_utilities/process_utilities.h"
#include "includes/kratos_parameters.h"

#include <string>

using namespace std::string_literals;

namespace Kratos
{
FindNeighboursOfInterfacesProcess::FindNeighboursOfInterfacesProcess(Model& rModel, const Parameters& rProcessSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(
        rModel, rProcessSettings, FindNeighboursOfInterfacesProcess::Info());
}

FindNeighboursOfInterfacesProcess::~FindNeighboursOfInterfacesProcess() = default;

std::string FindNeighboursOfInterfacesProcess::Info() const
{
    return "FindNeighboursOfInterfacesProcess"s;
}

} // namespace Kratos
