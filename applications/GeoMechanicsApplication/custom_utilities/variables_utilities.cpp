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

const Variable<double>& VariablesUtilities::GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                           const std::string& rComponent)
{
    return KratosComponents<Variable<double>>::Get(rSourceVariableName + "_" + rComponent);
}

} // namespace Kratos
