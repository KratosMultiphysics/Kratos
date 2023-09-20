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

struct ProcessParameters
{
    Parameters parameters;
    std::string name;

    ProcessParameters(const Parameters& rParameters, const std::string& rName) :
            parameters{rParameters},
            name{rName}
    {}

    bool operator==(const ProcessParameters& rhs) const
    {
        return name == rhs.name && parameters.WriteJsonString() == rhs.parameters.WriteJsonString();
    }

};

}