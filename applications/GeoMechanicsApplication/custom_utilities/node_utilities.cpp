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
                                                                    const array_1d<double, 3>& rNewValues,
                                                                    IndexType SolutionStepIndex)
{
    const std::vector<std::string> components = {"X", "Y", "Z"};
    for (const auto& [new_value, component] : boost::combine(rNewValues, components)) {
        if (const auto& component_variable =
                VariablesUtilities::GetComponentFromVectorVariable(rDestinationVariable.Name(), component);
            !rNode.IsFixed(component_variable)) {
            rNode.FastGetSolutionStepValue(component_variable, SolutionStepIndex) = new_value;
        }
    }
}

void NodeUtilities::AssignUpdatedVectorVariableToNodes(const ModelPart::NodesContainerType& rNodes,
                                                       const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                       const array_1d<double, 3>& rNewValues,
                                                       IndexType                  SolutionStepIndex)
{
    block_for_each(rNodes, [&rDestinationVariable, &rNewValues, SolutionStepIndex](Node& rNode) {
        rNode.FastGetSolutionStepValue(rDestinationVariable, SolutionStepIndex) = rNewValues;
    });
}

std::map<IndexType, IndexType> NodeUtilities::CreateGlobalToLocalNodeIndexMap(const PointerVector<Node>& rNodes)
{
    std::map<IndexType, IndexType> result;
    int                            counter = 0;
    for (const auto& r_node : rNodes) {
        result[r_node.GetId()] = counter;
        ++counter;
    }

    return result;
}
} // namespace Kratos
