//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/apply_function_to_nodes_utility.h"

namespace Kratos
{

void ApplyFunctionToNodesUtility::ApplyFunction(
    const Variable<double>& rVariable,
    const double t,
    const IndexType Step
    )
{
    // Get function
    auto& r_function = *mpFunction;

    if(!mpFunction->UseLocalSystem()) {
        block_for_each(
            mrNodes,r_function,
            [&rVariable,&t,Step](Node<3>& rNode, GenericFunctionUtility& rFunction) {
                const double value = rFunction.CallFunction(rNode.X(), rNode.Y(), rNode.Z(), t, rNode.X0(), rNode.Y0(), rNode.Z0());
                rNode.FastGetSolutionStepValue(rVariable, Step) = value;
            }
        );
    } else {
        block_for_each(
            mrNodes,r_function,
            [&rVariable,&t,Step](Node<3>& rNode, GenericFunctionUtility& rFunction) {
                const double value = rFunction.RotateAndCallFunction(rNode.X(), rNode.Y(), rNode.Z(), t, rNode.X0(), rNode.Y0(), rNode.Z0());
                rNode.FastGetSolutionStepValue(rVariable, Step) = value;
            }
        );
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<double> ApplyFunctionToNodesUtility::ReturnFunction(const double t)
{
    // The first node iterator
    const auto it_node_begin = mrNodes.begin();

    // The vector containing the values
    std::vector<double> values(mrNodes.size());

    // Get function
    auto& r_function = *mpFunction;

    if(!mpFunction->UseLocalSystem()) {
        IndexPartition<std::size_t>(mrNodes.size()).for_each(r_function,
            [&it_node_begin,&t,&values](std::size_t& k, GenericFunctionUtility& rFunction) {
                auto it_node = it_node_begin + k;
                const double value = rFunction.CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        );
    } else {
        IndexPartition<std::size_t>(mrNodes.size()).for_each(r_function,
            [&it_node_begin,&t,&values](std::size_t& k, GenericFunctionUtility& rFunction) {
                auto it_node = it_node_begin + k;
                const double value = rFunction.RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        );
    }

    return values;
}

} /// namespace Kratos
