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

#include "custom_processes/update_flag_no_remesh_femdem_boundary_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

UpdateFlagNoRemeshFemDemBoundaryProcess::UpdateFlagNoRemeshFemDemBoundaryProcess(
    ModelPart& rModelPart : mrModelPart(rModelPart),
{
}

/***********************************************************************************/
/***********************************************************************************/

void UpdateFlagNoRemeshFemDemBoundaryProcess::Execute() 
{
    const auto it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
    }
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos