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
#include "node_utilities.h"
#include "includes/node.h"
#include "variables_utilities.hpp"
#include <string>
#include <vector>

namespace Kratos
{

void NodeUtilities::ApplyUpdatedVectorVariableToNonFixedComponents(Node& rNode,
                                                                   const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                                   const array_1d<double, 3>& rSourceVector)
{
    const std::vector<std::string> components = {"X", "Y", "Z"};
    int                            counter    = 0;
    for (const auto& component : components) {
        const auto& component_variable =
            VariablesUtilities::GetComponentFromVectorVariable(rDestinationVariable, component);

        if (!rNode.IsFixed(component_variable)) {
            rNode.FastGetSolutionStepValue(component_variable, 0) = rSourceVector[counter];
        }
        counter++;
    }
}

} // namespace Kratos
