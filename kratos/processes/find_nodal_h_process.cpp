//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <limits>

// External includes

// Project includes
#include "processes/find_nodal_h_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

template<bool THistorical>
void FindNodalHProcess<THistorical>::Execute()
{
    KRATOS_TRY

    // Check if variables are available
    if (THistorical) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_H )) << "Variable NODAL_H not in the model part!" << std::endl;
    }

    // Initialize NODAL_H values
    const double max = std::numeric_limits<double>::max();
    if constexpr(THistorical) {
        VariableUtils().SetVariable(NODAL_H, max, mrModelPart.Nodes());
    } else {
        VariableUtils().SetNonHistoricalVariable(NODAL_H, max, mrModelPart.Nodes());
    }

    // Calculate the NODAL_H values
    for(auto& r_element : mrModelPart.Elements()) {
        auto& r_geom = r_element.GetGeometry();
        const SizeType number_of_nodes = r_geom.size();

        for(IndexType k = 0; k < number_of_nodes-1; ++k) {
            double& r_h1 = GetHValue(r_geom[k]);
            for(IndexType l=k+1; l < number_of_nodes; ++l) {
                const double hedge = r_geom[l].Distance(r_geom[k]);
                double& r_h2 = GetHValue(r_geom[l]);

                // Get minimum between the existent value and the considered edge length
                if constexpr(THistorical) {
                    r_geom[k].FastGetSolutionStepValue(NODAL_H) = std::min(r_h1, hedge);
                    r_geom[l].FastGetSolutionStepValue(NODAL_H) = std::min(r_h2, hedge);
                } else {
                    r_geom[k].GetValue(NODAL_H) = std::min(r_h1, hedge);
                    r_geom[l].GetValue(NODAL_H) = std::min(r_h2, hedge);
                }
            }
        }
    }

    // Synchronize between processes
    if constexpr(THistorical) {
        mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(NODAL_H);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalDataToMin(NODAL_H);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& FindNodalHProcess<true>::GetHValue(NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(NODAL_H);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& FindNodalHProcess<false>::GetHValue(NodeType& rNode)
{
    return rNode.GetValue(NODAL_H);
}

/***********************************************************************************/
/***********************************************************************************/

template class FindNodalHProcess<true>;
template class FindNodalHProcess<false>;

} // namespace Kratos
