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


namespace  Kratos
{
    void MoveMeshUtility::MovePfemMesh(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF_NOT(rNodes.begin()->SolutionStepsDataHas(DISPLACEMENT_X))
          << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part"
          << std::endl;

        const int num_nodes = static_cast<int>(rNodes.size());

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;

            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
        }

        KRATOS_INFO("MoveMeshUtility") << " PFEM MESH MOVED " << std::endl;
        KRATOS_CATCH("")
    }
    void MoveMeshUtility::ResetPfemKinematicValues(ModelPart& rFluidModelPart)
    {
        const auto it_node_begin = rFluidModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rFluidModelPart.Nodes().size()); ++i) {
            auto it_node = it_node_begin + i;

            if (it_node->IsNot(RIGID) && it_node->IsNot(SOLID)) { // We update only the fluid part
                auto &r_current_displ = it_node->FastGetSolutionStepValue(DISPLACEMENT, 0);
                auto &r_current_vel   = it_node->FastGetSolutionStepValue(VELOCITY, 0);
                auto &r_current_acc   = it_node->FastGetSolutionStepValue(ACCELERATION, 0);

                auto &r_old_displ     = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
                auto &r_old_vel       = it_node->FastGetSolutionStepValue(VELOCITY, 1);
                auto &r_old_acc       = it_node->FastGetSolutionStepValue(ACCELERATION, 1);

                auto copy_old_displ   = r_old_displ;
                auto copy_old_vel     = r_old_vel;
                auto copy_old_acc     = r_old_acc;

                auto& r_coordinates = it_node->Coordinates();
                const auto& r_initial_coordinates = it_node->GetInitialPosition();
                noalias(r_coordinates) = r_initial_coordinates + copy_old_displ;

                noalias(r_current_displ) = copy_old_displ;
                noalias(r_current_vel)   = copy_old_vel;
                noalias(r_current_acc)   = copy_old_acc;
            }
        }

        KRATOS_INFO("MoveMeshUtility::ResetPfemKinematicValues") << " PFEM KINEMATICS RESET " << std::endl;
    }
} //  Kratos
