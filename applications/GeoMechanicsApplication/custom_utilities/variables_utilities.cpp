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

#include "variables_utilities.hpp"
#include "includes/kratos_components.h"

namespace Kratos
{

const Variable<double>& VariablesUtilities::GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                           const std::string& rComponent)
{
    const auto component_name = rSourceVariableName + "_" + rComponent;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(component_name))
        << "The component \"" << component_name << "\" is not registered!" << std::endl;

    return KratosComponents<Variable<double>>::Get(component_name);
}

} // namespace Kratos
