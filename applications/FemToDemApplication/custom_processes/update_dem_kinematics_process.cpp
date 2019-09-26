//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
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

void UpdateDemKinematicsProcess::UpdateKinematics(
    const NodeIteratorType& rNode,
    NodeType& rDEMNode
    )
{
    const array_1d<double,3> coordinates = rNode->Coordinates();
    const array_1d<double,3> displacement = rNode->GetSolutionStepValue(DISPLACEMENT);
    const array_1d<double,3> velocity = rNode->GetSolutionStepValue(VELOCITY);

    auto& r_displ_dem = rDEMNode.GetSolutionStepValue(DISPLACEMENT);
    r_displ_dem = displacement;
    auto& r_vel_dem = rDEMNode.GetSolutionStepValue(VELOCITY);
    r_vel_dem = velocity;

    array_1d<double,3>& r_dem_coordinates = rDEMNode.Coordinates();
    r_dem_coordinates = coordinates;
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos