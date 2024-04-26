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

#include <boost/range/combine.hpp>
#include <string>
#include <vector>

namespace Kratos
{

void NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(Node& rNode,
                                                                    const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                                    const array_1d<double, 3>& rNewValues)
{
    const std::vector<std::string> components = {"X", "Y", "Z"};
    for (const auto& zipped : boost::combine(rNewValues, components)) {
        double      new_value = 0.0;
        std::string component;
        boost::tie(new_value, component) = zipped;

        if (const auto& component_variable =
                VariablesUtilities::GetComponentFromVectorVariable(rDestinationVariable.Name(), component);
            !rNode.IsFixed(component_variable)) {
            rNode.FastGetSolutionStepValue(component_variable, 0) = new_value;
        }
    }
}

} // namespace Kratos
