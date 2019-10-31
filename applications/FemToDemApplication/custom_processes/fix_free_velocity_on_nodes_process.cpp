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

#include "custom_processes/fix_free_velocity_on_nodes_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

FixFreeVelocityOnNodesProcess::FixFreeVelocityOnNodesProcess(
    ModelPart& rModelPart,
    const int rFreeOrFix)
    : mrModelPart(rModelPart),
      mFreeOrFix(rFreeOrFix)
{
}

/***********************************************************************************/
/***********************************************************************************/

void FixFreeVelocityOnNodesProcess::Execute() 
{
    const auto& it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        if (mFreeOrFix == 0) {
            it_node->Fix(VELOCITY_X);
            it_node->Fix(VELOCITY_Y);
            it_node->Fix(VELOCITY_Z);
        } else {
            it_node->Free(VELOCITY_X);
            it_node->Free(VELOCITY_Y);
            it_node->Free(VELOCITY_Z);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos