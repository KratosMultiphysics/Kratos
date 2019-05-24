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

#include "custom_processes/update_dem_kinematics.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

UpdateDemKinematicsProcess::UpdateDemKinematicsProcess(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void UpdateDemKinematicsProcess::Execute() 
{
    const auto it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(IS_DEM)) {
            const int node_id = it_node->Id();
            auto& associated_dem = mrDEMModelPart.GetNode(node_id);
            const array_1d<double,3> coordinates = this->GetNodeCoordinates(it_node);
            const array_1d<double,3> displacement = it_node->GetSolutionStepValue(DISPLACEMENT)
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double,3> UpdateDemKinematicsProcess::GetNodeCoordinates(
    const NodeType& rNode
    )
{
    array_1d<double,3> coordinates;
    coordinates[0] = rNode.X();
    coordinates[1] = rNode.Y();
    coordinates[2] = rNode.Z();
    return coordinates;
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos