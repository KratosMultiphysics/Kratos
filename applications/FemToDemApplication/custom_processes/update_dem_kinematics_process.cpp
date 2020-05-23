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
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(IS_DEM)) {
            auto p_associated_dem = it_node->GetValue(DEM_PARTICLE_POINTER);
            this->UpdateKinematics(it_node, p_associated_dem->GetGeometry()[0]);
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
    const array_1d<double,3> displacement = rNode->FastGetSolutionStepValue(DISPLACEMENT);
    const array_1d<double,3> delta_displacement = rNode->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    const array_1d<double,3> velocity = rNode->FastGetSolutionStepValue(VELOCITY);

    auto& r_displ_dem = rDEMNode.FastGetSolutionStepValue(DISPLACEMENT);
    r_displ_dem = displacement;
    auto& r_vel_dem = rDEMNode.FastGetSolutionStepValue(VELOCITY);
    r_vel_dem = velocity;

    array_1d<double,3>& r_dem_coordinates = rDEMNode.Coordinates();
    array_1d<double, 3>& r_dem_delta_displ = rDEMNode.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    r_dem_delta_displ = delta_displacement;
    r_dem_coordinates = coordinates;
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos