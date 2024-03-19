// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "variables_utilities.hpp"
#include "includes/kratos_components.h"

namespace Kratos
{

const Variable<double>& VariablesUtilities::GetComponentFromVectorVariable(const Variable<array_1d<double, 3>>& rSource,
                                                                           const std::string& rComponent)
{
    return KratosComponents<Variable<double>>::Get(rSource.Name() + "_" + rComponent);
}

} // namespace Kratos
