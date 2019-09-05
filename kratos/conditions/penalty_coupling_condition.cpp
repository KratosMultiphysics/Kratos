//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//

// System includes

// External includes
#include "conditions/penalty_coupling_condition.h"

// Project includes

namespace Kratos
{

    void PenaltyCouplingCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        DofVariablesContainer const& r_dof_variables = this->GetProperties().GetValue(DOF_VARIABLES);

        const double penalty = GetProperties()[INITIAL_PENALTY];

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const int number_of_nodes_master = r_geometry_master.size();
        const int number_of_nodes_slave = r_geometry_slave.size();

        const int mat_size = r_dof_variables.size() * (number_of_nodes_master + number_of_nodes_slave);

        // Memory allocation
        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }
        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
        }

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            //FOR DISPLACEMENTS
            Matrix H = ZeroMatrix(3, mat_size);
            IndexType local_matrix_index = 0;
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                for (IndexType j = 0; j < r_dof_variables.size(); ++j)
                {
                    H(j, local_matrix_index++) = N_master(point_number, i);
                }
            }

            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                for (IndexType j = 0; j < r_dof_variables.size(); ++j)
                {
                    H(j, local_matrix_index++) = -N_slave(point_number, i);
                }
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry_master.DeterminantOfJacobian(point_number);

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                    * integration_weight * determinat_jacobian * penalty;
            }
            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                IndexType local_index = 0;
                for (IndexType i = 0; i < number_of_nodes_master; i++)
                {
                    const auto& r_node = r_geometry_master[i];

                    for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                        itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                        u[local_index++] = r_node.FastGetSolutionStepValue(*itDVar);

                    for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                        itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                        u[local_index++] = r_node.FastGetSolutionStepValue(*itCVar);
                }
                for (IndexType i = 0; i < number_of_nodes_slave; i++)
                {
                    const auto& r_node = r_geometry_slave[i];

                    for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                        itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                        u[local_index++] = r_node.FastGetSolutionStepValue(*itDVar);

                    for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                        itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                        u[local_index++] = r_node.FastGetSolutionStepValue(*itCVar);
                }

                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinat_jacobian * penalty;
            }
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        DofVariablesContainer const& r_dof_variables = this->GetProperties().GetValue(DOF_VARIABLES);

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != r_dof_variables.size() * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(r_dof_variables.size() * (number_of_nodes_master + number_of_nodes_slave), false);

        IndexType local_index = 0;

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master[i];

            for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                rResult[local_index++] = r_node.GetDof(*itDVar).EquationId();

            for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                rResult[local_index++] = r_node.GetDof(*itCVar).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave[i];

            for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                rResult[local_index++] = r_node.GetDof(*itDVar).EquationId();

            for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                rResult[local_index++] = r_node.GetDof(*itCVar).EquationId();
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        DofVariablesContainer const& r_dof_variables = this->GetProperties().GetValue(DOF_VARIABLES);

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(r_dof_variables.size() * (number_of_nodes_master + number_of_nodes_slave));

        IndexType local_index = 0;

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                rElementalDofList[local_index++] = r_node.pGetDof(*itDVar);

            for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                rElementalDofList[local_index++] = r_node.pGetDof(*itCVar);
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            for (auto itDVar = r_dof_variables.DoubleVariablesBegin();
                itDVar != r_dof_variables.DoubleVariablesEnd(); ++itDVar)
                rElementalDofList[local_index++] = r_node.pGetDof(*itDVar);

            for (auto itCVar = r_dof_variables.VariableComponentsBegin();
                itCVar != r_dof_variables.VariableComponentsEnd(); ++itCVar)
                rElementalDofList[local_index++] = r_node.pGetDof(*itCVar);
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


