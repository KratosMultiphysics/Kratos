//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/transfer_nodal_forces_to_fem.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

TransferNodalForcesToFem::TransferNodalForcesToFem(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void TransferNodalForcesToFem::Execute() 
{
    const auto it_node_begin = mrModelPart.NodesBegin();
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(IS_DEM)) {
            const int node_id = it_node->Id();
            auto& associated_dem = mrDEMModelPart.GetNode(node_id);
            this->UpdateKinematics(it_node, associated_dem);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos