//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter (based on Riccardo Rossi's work)
//
// System includes


// External includes


// Project includes
#include "move_mesh_utility.h"
#include "DEM_application_variables.h"


namespace  Kratos
{


    void MoveMeshUtility::CalculateDeltaDispCustom(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;

        const int num_nodes = static_cast<int>(rNodes.size());

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            noalias(it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) =
             it_node->GetInitialPosition().Coordinates()+it_node->FastGetSolutionStepValue(DISPLACEMENT)-it_node->Coordinates();
        }

        KRATOS_INFO("MoveMeshUtility") << " DELTA_DISPLACEMENT successfully set " << std::endl;
        KRATOS_CATCH("")
    }

    void MoveMeshUtility::CalculateDeltaDispCustomFromIntermediatePos(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;

        const int num_nodes = static_cast<int>(rNodes.size());

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            noalias(it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) =
             it_node->GetInitialPosition().Coordinates()+it_node->FastGetSolutionStepValue(DISPLACEMENT)-it_node->GetIntermediatePosition();
        }

        KRATOS_INFO("MoveMeshUtility") << " DELTA_DISPLACEMENT successfully set from intermediate position " << std::endl;
        KRATOS_CATCH("")
    }


    void MoveMeshUtility::MoveDemMesh(NodesContainerType& rNodes, const bool& rSetDeltaDisplacement) const
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF_NOT(rNodes.begin()->SolutionStepsDataHas(DISPLACEMENT_X))
          << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part"
          << std::endl;

        const int num_nodes = static_cast<int>(rNodes.size());


        if(!rSetDeltaDisplacement)
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
                auto old_coordinates = it_node->Coordinates();

                noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
                noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);

                /* noalias(it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) =
                    it_node->FastGetSolutionStepValue(DISPLACEMENT,0)
                    -it_node->FastGetSolutionStepValue(DISPLACEMENT,1); */
                noalias(it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) = it_node->Coordinates()-old_coordinates;
            }
        }

        KRATOS_INFO("MoveMeshUtility") << " DEM MESH MOVED " << std::endl;
        KRATOS_CATCH("")
    }

    const bool MoveMeshUtility::CheckContact(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;
        const int num_nodes = static_cast<int>(rNodes.size());
        bool check_flag = false;

        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            BoundedVector<double,3> contact_forces = it_node->FastGetSolutionStepValue(CONTACT_FORCES);

            for(int j = 0; j < 3 ; ++j)
            {
                if (fabs(contact_forces[j])>0.00)
                {
                    check_flag = true;
                    break;
                }
            }
            if (check_flag) break;
        }
        return check_flag;
        KRATOS_CATCH("")
    }


    const bool MoveMeshUtility::CheckIsNearToWall(NodesContainerType& rNodes,Vector& rMaxDisplacement) const
    {
        KRATOS_TRY;
        //const double numerical_limit = std::numeric_limits<double>::epsilon();
        bool check_flag = false;
        const int num_nodes = static_cast<int>(rNodes.size());

        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            Vector current_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
            for (int j = 0; j < 3; ++j)
            {
                if ((rMaxDisplacement[j]-fabs(current_displacement[j]))<0.00)
                {
                    check_flag = true;
                    break;
                }
            }
            if (check_flag) break;
        }
        return check_flag;
        KRATOS_CATCH("")
    }




    const void MoveMeshUtility::SaveCurrentCoordinates(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;

        const int num_nodes = static_cast<int>(rNodes.size());
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            it_node->SetIntermediatePosition(it_node->Coordinates());
        }

        KRATOS_CATCH("")
    }

    const void MoveMeshUtility::ResetCoordinates(NodesContainerType& rNodes) const
    {
        KRATOS_TRY;

        const int num_nodes = static_cast<int>(rNodes.size());
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = rNodes.begin() + i;
            noalias(it_node->Coordinates()) = it_node->GetIntermediatePosition();
        }

        KRATOS_CATCH("")
    }


    const void MoveMeshUtility::SetInactiveElements(ConditionsContainerType& rElementsDEM, const ElementsContainerType& rElementsFEM)
    {
        KRATOS_TRY;
        const int num_elements_dem = static_cast<int>(rElementsDEM.size());
        const int num_elements_fem = static_cast<int>(rElementsFEM.size());
        KRATOS_ERROR_IF(num_elements_dem != num_elements_fem) << " number of nodes is not equal" << std::endl;

        this->mDEMConditions.resize(num_elements_dem);

        #pragma omp parallel for
        for(int i = 0;i < num_elements_fem; ++i)
        {
            auto it_element = rElementsFEM.begin()+i;
            auto& node_a_coordinates = it_element->GetGeometry()[0].Coordinates();
            auto& node_b_coordinates = it_element->GetGeometry()[1].Coordinates();



            std::cout << node_a_coordinates[0] << node_a_coordinates[1] << node_a_coordinates[2] << std::endl;
            for (int j = 0;j < num_elements_dem; ++j)
            {
                auto jt_element = rElementsDEM.begin()+j;
                auto& node_k_coordinates = jt_element->GetGeometry()[0].Coordinates();
                auto& node_l_coordinates = jt_element->GetGeometry()[1].Coordinates();

                if (this->CompareTwoArrays(node_a_coordinates,node_k_coordinates) && this->CompareTwoArrays(node_b_coordinates,node_l_coordinates))
                {
                    this->mDEMConditions[i] = jt_element;
                    break;
                }
                else if (this->CompareTwoArrays(node_a_coordinates,node_l_coordinates) && this->CompareTwoArrays(node_b_coordinates,node_k_coordinates))
                {
                    this->mDEMConditions[i] = jt_element;
                    break;
                }

            }

        }

        KRATOS_CATCH("")
    }

} //  Kratos
