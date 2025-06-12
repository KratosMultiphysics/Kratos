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
#pragma once

#include "includes/kratos_parameters.h"

#include <string>

namespace Kratos
{

struct ProcessParameters {
    std::string name;
    Parameters  parameters;

    ProcessParameters(const std::string& rName, const Parameters& rParameters)
        : name{rName}, parameters{rParameters}
    {
    }

    bool operator==(const ProcessParameters& rhs) const
    {
        return name == rhs.name && parameters.WriteJsonString() == rhs.parameters.WriteJsonString();
    }
};

} // namespace Kratos