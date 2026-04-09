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

struct KRATOS_API(GEO_MECHANICS_APPLICATION) ProcessParameters {
    std::string name;
    Parameters  parameters;

    ProcessParameters(std::string Name, Parameters Settings);

    bool operator==(const ProcessParameters& rhs) const;
};

} // namespace Kratos