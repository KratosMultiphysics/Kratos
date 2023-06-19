//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//                   Tobias Teschemacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_penalty_condition.h"

namespace Kratos
{
    void SupportPenaltyCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry, determinant_jacobian_vector);
        } else {
            r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            const Matrix& N = r_geometry.ShapeFunctionsValues();

            //FOR DISPLACEMENTS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                IndexType index = 3 * i;
                H(0, index)     = N(point_number, i);
                H(1, index + 1) = N(point_number, i);
                H(2, index + 2) = N(point_number, i);
            }

            // Differential area
            const double penalty_integration = penalty * integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += prod(trans(H), H) * penalty_integration;
            }
            if (CalculateResidualVectorFlag) {

                const array_1d<double, 3>& displacement = Has(DISPLACEMENT)
                    ? this->GetValue(DISPLACEMENT)
                    : ZeroVector(3);

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index] = (disp[0] - displacement[0]);
                    u[index + 1] = (disp[1] - displacement[1]);
                    u[index + 2] = (disp[2] - displacement[2]);
                }

                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u) * penalty_integration;
            }
        }
        KRATOS_CATCH("")
    }

    void SupportPenaltyCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType number_of_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType point_number = 0; point_number < nb_integration_points; ++point_number)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[point_number];
            J.clear();
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += r_coordinates[k] * r_DN_De(i, m);
                    }
                }
            }

            //Compute the tangent and  the normal to the boundary vector
            array_1d<double, 3> local_tangent;
            GetGeometry().Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[point_number] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }

    int SupportPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyCondition" << std::endl;
        return 0;
    }

    void SupportPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void SupportPenaltyCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    };

} // Namespace Kratos
