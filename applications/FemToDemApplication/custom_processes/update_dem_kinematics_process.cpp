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

#include "custom_processes/update_dem_kinematics_process.h"
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
            this->UpdateKinematics(it_node, associated_dem);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double,3> UpdateDemKinematicsProcess::GetNodeCoordinates(
    const NodeIteratorType& rNode
    )
{
    array_1d<double,3> coordinates;
    coordinates[0] = rNode->X();
    coordinates[1] = rNode->Y();
    coordinates[2] = rNode->Z();
    return coordinates;
}

/***********************************************************************************/
/***********************************************************************************/

void UpdateDemKinematicsProcess::UpdateKinematics(
    const NodeIteratorType& rNode,
    NodeType& rDEMNode
    )
{
    const array_1d<double,3> coordinates = this->GetNodeCoordinates(rNode);
    const array_1d<double,3> displacement = rNode->GetSolutionStepValue(DISPLACEMENT);
    const array_1d<double,3> velocity = rNode->GetSolutionStepValue(VELOCITY);

    auto& displ_dem = rDEMNode.GetSolutionStepValue(DISPLACEMENT);
    displ_dem = displacement;
    auto& vel_dem = rDEMNode.GetSolutionStepValue(VELOCITY);
    vel_dem = velocity;

    double& r_x_dem =  rDEMNode.X();
    double& r_y_dem =  rDEMNode.Y();
    double& r_z_dem =  rDEMNode.Z();

    r_x_dem = coordinates[0];
    r_y_dem = coordinates[1];
    r_z_dem = coordinates[2];
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos