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

#include "process_parameters.h"

namespace Kratos
{
ProcessParameters::ProcessParameters(const std::string& rName, const Parameters& rParameters)
    : name{rName}, parameters{rParameters}
{
}

bool ProcessParameters::operator==(const ProcessParameters& rhs) const
{
    return name == rhs.name && parameters.WriteJsonString() == rhs.parameters.WriteJsonString();
}

} // namespace Kratos