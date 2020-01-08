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
    auto &r_process_info = mrModelPart.GetProcessInfo();
    const auto& it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        if (mFreeOrFix == 0) {
            it_node->Fix(VELOCITY_X);
            it_node->Fix(VELOCITY_Y);
            it_node->Fix(VELOCITY_Z);

            // We store previous info in order to reset after PFEM Solve
            if (r_process_info[STEP] > 1) { // We have the previous acc
                const auto current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 0);
                const auto prev_acceleration    = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
                auto &r_current_accel_backup = it_node->FastGetSolutionStepValue(ACCELERATION_BACKUP, 0);
                auto &r_prev_accel_backup    = it_node->FastGetSolutionStepValue(ACCELERATION_BACKUP, 1);
                noalias(r_current_accel_backup) = current_acceleration;
                noalias(r_prev_accel_backup)  = prev_acceleration;

                const auto current_displ = it_node->FastGetSolutionStepValue(DISPLACEMENT, 0);
                const auto prev_displ    = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
                auto &r_current_displ_backup = it_node->FastGetSolutionStepValue(DISPLACEMENT_BACKUP, 0);
                auto &r_prev_displ_backup    = it_node->FastGetSolutionStepValue(DISPLACEMENT_BACKUP, 1);
                noalias(r_current_displ_backup) = current_displ;
                noalias(r_prev_displ_backup) = prev_displ;
            }

        } else {
            it_node->Free(VELOCITY_X);
            it_node->Free(VELOCITY_Y);
            it_node->Free(VELOCITY_Z);

            // We store previous info in order to reset after PFEM Solve
            if (r_process_info[STEP] > 1) { // We have the previous acc
                auto &r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 0);
                auto &r_prev_acceleration    = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
                const auto current_accel_backup = it_node->FastGetSolutionStepValue(ACCELERATION_BACKUP, 0);
                const auto prev_accel_backup    = it_node->FastGetSolutionStepValue(ACCELERATION_BACKUP, 1);
                noalias(r_current_acceleration) = current_accel_backup;
                noalias(r_prev_acceleration) = prev_accel_backup;

                auto &r_current_displ = it_node->FastGetSolutionStepValue(DISPLACEMENT, 0);
                auto &r_prev_displ    = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
                const auto current_displ_backup = it_node->FastGetSolutionStepValue(DISPLACEMENT_BACKUP, 0);
                const auto prev_displ_backup    = it_node->FastGetSolutionStepValue(DISPLACEMENT_BACKUP, 1);
                noalias(r_current_displ) = current_displ_backup;
                noalias(r_prev_displ) = prev_displ_backup;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


}  // namespace Kratos
