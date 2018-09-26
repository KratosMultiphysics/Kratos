//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
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
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(NODAL_H) = std::numeric_limits<double>::max();
    }
    
    for(IndexType i=0; i < mrModelPart.Elements().size(); ++i) {
        auto it_element = mrModelPart.ElementsBegin() + i;
        auto& r_geom = it_element->GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        
        for(IndexType k = 0; k < number_of_nodes-1; ++k) {
            const double h1 = GetValue(r_geom[k]);
            for(IndexType l=k+1; l < number_of_nodes; ++l) {
                double hedge = norm_2(r_geom[l].Coordinates() - r_geom[k].Coordinates());
                const double h2 = GetValue(r_geom[l]);
                
                // Get minimum between the existent value and the considered edge length 
                SetValue(r_geom[k], std::min(h1, hedge));
                SetValue(r_geom[l], std::min(h2, hedge));
            }
        }
    }
    
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(NODAL_H);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double FindNodalHProcess<true>::GetValue(NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(NODAL_H);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double FindNodalHProcess<false>::GetValue(NodeType& rNode)
{
    return rNode.GetValue(NODAL_H);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindNodalHProcess<true>::SetValue(
    NodeType& rNode,
    const double Value
    )
{
    rNode.FastGetSolutionStepValue(NODAL_H) = Value;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindNodalHProcess<false>::SetValue(
    NodeType& rNode,
    const double Value
    )
{
    rNode.SetValue(NODAL_H, Value);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindNodalHProcess<true>::SetInitialValue(NodeIterator itNode)
{
    itNode->FastGetSolutionStepValue(NODAL_H) = std::numeric_limits<double>::max();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindNodalHProcess<false>::SetInitialValue(NodeIterator itNode)
{
    itNode->SetValue(NODAL_H, std::numeric_limits<double>::max());
}

/***********************************************************************************/
/***********************************************************************************/

template class FindNodalHProcess<true>;
template class FindNodalHProcess<false>;

} // namespace Kratos
