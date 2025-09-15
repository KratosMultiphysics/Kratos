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
//
// System includes


// External includes


// Project includes
#include "move_mesh_utility.h"
#include "DEM_application_variables.h"


namespace  Kratos
{
    void MoveMeshUtility::MoveDemMesh(NodesContainerType& rNodes, const bool SetDeltaDisplacement) const
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF_NOT(rNodes.begin()->SolutionStepsDataHas(DISPLACEMENT_X))
          << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part"
          << std::endl;

        const int num_nodes = static_cast<int>(rNodes.size());


        if(!SetDeltaDisplacement)
        {
            #pragma omp parallel for
            for(int i = 0; i < num_nodes; ++i)
            {
                auto it_node = rNodes.begin() + i;

                noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
                noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
            }
        }

        else
        {
            #pragma omp parallel for
            for(int i = 0; i < num_nodes; ++i)
            {
                auto it_node = rNodes.begin() + i;

                BoundedVector<double,3> previous_position = it_node->Coordinates();
                noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
                noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);

                noalias(it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) =
                    it_node->Coordinates()-previous_position;

            }
        }

        KRATOS_INFO("MoveMeshUtility") << " DEM MESH MOVED " << std::endl;
        KRATOS_CATCH("")
    }
} //  Kratos
