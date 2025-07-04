//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//
//

// System includes
#include <limits>

// External includes

// Project includes
#include "find_nodal_h_process_max.h"
#include "utilities/variable_utils.h"
#include "droplet_dynamics_application_variables.h"

namespace Kratos
{

template<bool THistorical>
void FindNodalHProcessMax<THistorical>::Execute()
{
    KRATOS_TRY

    // Check if variables are available
    if (THistorical) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_H_MAX )) << "Variable NODAL_H_MAX not in the model part!" << std::endl;
    }

    // Initialize NODAL_H_MAX values
    const double min = std::numeric_limits<double>::lowest();
    if constexpr(THistorical) {
        VariableUtils().SetVariable(NODAL_H_MAX, min, mrModelPart.Nodes());
    } else {
        VariableUtils().SetNonHistoricalVariable(NODAL_H_MAX, min, mrModelPart.Nodes());
    }

    // Calculate the NODAL_H_MAX values
    for(auto& r_element : mrModelPart.Elements()) {
        auto& r_geom = r_element.GetGeometry();
        const SizeType number_of_nodes = r_geom.size();

      for(IndexType k = 0; k < number_of_nodes-1; ++k) {
        double& r_h1 = GetHValue(r_geom[k]);
        for(IndexType l=k+1; l < number_of_nodes; ++l) {
            const double hedge = r_geom[l].Distance(r_geom[k]);
            double& r_h2 = GetHValue(r_geom[l]);

            // Get maximum between the existent value and the considered edge length
            if constexpr(THistorical) {
                r_geom[k].FastGetSolutionStepValue(NODAL_H_MAX) = std::max(r_h1, hedge);
                r_geom[l].FastGetSolutionStepValue(NODAL_H_MAX) = std::max(r_h2, hedge);
            } else {
                r_geom[k].GetValue(NODAL_H_MAX) = std::max(r_h1, hedge);
                r_geom[l].GetValue(NODAL_H_MAX) = std::max(r_h2, hedge);
                }
            }
        }
    }

    // Synchronize between processes
    if constexpr(THistorical) {
        mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(NODAL_H_MAX);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalDataToMin(NODAL_H_MAX);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& FindNodalHProcessMax<true>::GetHValue(NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(NODAL_H_MAX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& FindNodalHProcessMax<false>::GetHValue(NodeType& rNode)
{
    return rNode.GetValue(NODAL_H_MAX);
}

/***********************************************************************************/
/***********************************************************************************/

template class FindNodalHProcessMax<true>;
template class FindNodalHProcessMax<false>;

} // namespace Kratos
