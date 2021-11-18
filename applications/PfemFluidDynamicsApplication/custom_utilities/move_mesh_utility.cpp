//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
//                   Alejandro Cornejo Velazquez
//                   Carlos Eulogio Flores
//
// System includes


// External includes


// Project includes
#include "move_mesh_utility.h"
#include "pfem_fluid_dynamics_application_variables.h"
#include "utilities/parallel_utilities.h"


namespace  Kratos
{
    void MoveMeshUtility::ResetPfemKinematicValues(ModelPart& rFluidModelPart)
    {
        block_for_each(rFluidModelPart.Nodes(), [&](NodeType& rNode){
            if (rNode.IsNot(RIGID) && rNode.IsNot(SOLID)) { // We update only the fluid part
                auto &r_current_displ = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
                auto &r_current_vel   = rNode.FastGetSolutionStepValue(VELOCITY, 0);
                auto &r_current_acc   = rNode.FastGetSolutionStepValue(ACCELERATION, 0);

                auto &r_old_displ     = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
                auto &r_old_vel       = rNode.FastGetSolutionStepValue(VELOCITY, 1);
                auto &r_old_acc       = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

                auto copy_old_displ   = r_old_displ;
                auto copy_old_vel     = r_old_vel;
                auto copy_old_acc     = r_old_acc;

                auto& r_coordinates = rNode.Coordinates();
                const auto& r_initial_coordinates = rNode.GetInitialPosition();
                noalias(r_coordinates) = r_initial_coordinates + copy_old_displ;

                noalias(r_current_displ) = copy_old_displ;
                noalias(r_current_vel)   = copy_old_vel;
                noalias(r_current_acc)   = copy_old_acc;
            }
        });

        KRATOS_INFO("MoveMeshUtility::ResetPfemKinematicValues") << " PFEM KINEMATICS RESET " << std::endl;
    }
} //  Kratos
