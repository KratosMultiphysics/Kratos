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
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "custom_conditions/penalty_directional_support_condition.h"

// Project includes

namespace Kratos
{

    void PenaltyDirectionalSupportCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];
        array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

        const auto& r_geometry = GetGeometry();
        const int number_of_nodes = r_geometry.size();

        const int mat_size = 3 * number_of_nodes;

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
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();
            array_1d<double, 3> direction_fixed(0.0);
            r_geometry.Normal(direction_fixed, point_number);
            Vector direction_fixed_vector = -direction_fixed / norm_2(direction_fixed);

            //FOR DISPLACEMENTS
            Matrix H = ZeroMatrix(3, mat_size);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                int index = 3 * i;
                H(0, index) = r_N(point_number, i);
                H(1, index + 1) = r_N(point_number, i);
                H(2, index + 2) = r_N(point_number, i);
            }

            Vector H_rotated = prod(direction_fixed_vector, H);

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = r_geometry.DeterminantOfJacobian(point_number);

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += outer_prod(trans(H_rotated), H_rotated)
                    * integration_weight * determinat_jacobian * penalty;
            }
            if (CalculateResidualVectorFlag) {

                if (displacement[0] > 0)
                {
                    //displacement = direction_fixed_vector * displacement[0];
                }

                Vector u(mat_size);
                for (unsigned int i = 0; i < number_of_nodes; i++)
                {
                    const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                    int index = 3 * i;
                    u[index] = disp[0] - displacement[0];
                    u[index + 1] = disp[1] * direction_fixed[1] - displacement[1];
                    u[index + 2] = disp[2] * direction_fixed[2] - displacement[2];
                }

                noalias(rRightHandSideVector) -= prod(outer_prod(trans(H_rotated), H_rotated), u)
                    * integration_weight * determinat_jacobian * penalty;
            }
        }

        KRATOS_CATCH("")
    }

    void PenaltyDirectionalSupportCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const int number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            const unsigned int index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void PenaltyDirectionalSupportCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const int number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


