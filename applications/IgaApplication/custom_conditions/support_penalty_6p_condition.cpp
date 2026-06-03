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
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_penalty_6p_condition.h"

namespace Kratos
{
    void SupportPenalty6pCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];
        const double penalty_rotation = GetProperties()[PENALTY_ROTATION_FACTOR];

        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();
        const IndexType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;
        const IndexType disp_size = 3 * number_of_nodes;

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

            // Shape function matrix for displacements and rotations
            Matrix H_disp = ZeroMatrix(3, disp_size);
            Matrix H_rot = ZeroMatrix(3, disp_size);

            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                IndexType index = 3 * i;

                H_disp(0, index)     = N(point_number, i);
                H_disp(1, index + 1) = N(point_number, i);
                H_disp(2, index + 2) = N(point_number, i);

                H_rot(0, index)     = N(point_number, i);
                H_rot(1, index + 1) = N(point_number, i);
                H_rot(2, index + 2) = N(point_number, i);
            }

            // Differential area
            const double penalty_integration = penalty * integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];
            const double penalty_rotation_integration = penalty_rotation * integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            // Matrix multiplication blocks
            const Matrix HtH_disp = prod(trans(H_disp), H_disp);
            const Matrix HtH_rot  = prod(trans(H_rot),  H_rot);

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                for (IndexType i = 0; i < number_of_nodes; ++i) {
                    for (IndexType j = 0; j < number_of_nodes; ++j) {
                        for (IndexType ii = 0; ii < 3; ++ii) {
                            for (IndexType jj = 0; jj < 3; ++jj) {
                                rLeftHandSideMatrix(6 * i + ii, 6 * j + jj) += HtH_disp(3 * i + ii, 3 * j + jj) * penalty_integration;
                                rLeftHandSideMatrix(6 * i + 3 + ii, 6 * j + 3 + jj) += HtH_rot(3 * i + ii, 3 * j + jj) * penalty_rotation_integration;
                            }
                        }
                    }
                }
            }
            if (CalculateResidualVectorFlag) {

                const array_1d<double, 3>& displacement = Has(DISPLACEMENT)
                    ? this->GetValue(DISPLACEMENT)
                    : ZeroVector(3);
                const array_1d<double, 3>& rotation = Has(ROTATION)
                    ? this->GetValue(ROTATION)
                    : ZeroVector(3);

                Vector u_disp(disp_size);
                Vector u_rot(disp_size);

                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                    const array_1d<double, 3> rot = r_geometry[i].FastGetSolutionStepValue(ROTATION);

                    IndexType index = 3 * i;
                    u_disp[index] = (disp[0] - displacement[0]);
                    u_disp[index + 1] = (disp[1] - displacement[1]);
                    u_disp[index + 2] = (disp[2] - displacement[2]);
                    
                    u_rot[index] = (rot[0] - rotation[0]);
                    u_rot[index + 1] = (rot[1] - rotation[1]);
                    u_rot[index + 2] = (rot[2] - rotation[2]);
                }

                const Vector rhs_disp = prod(HtH_disp, u_disp) * penalty_integration;
                const Vector rhs_rot  = prod(HtH_rot,  u_rot)  * penalty_rotation_integration;

                for (IndexType i = 0; i < number_of_nodes; ++i) {
                    for (IndexType ii = 0; ii < 3; ++ii) {
                        rRightHandSideVector(6 * i + ii) -= rhs_disp(3 * i + ii);
                        rRightHandSideVector(6 * i + 3 + ii) -= rhs_rot(3 * i + ii);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportPenalty6pCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const IndexType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const IndexType local_space_dimension = rGeometry.LocalSpaceDimension();
        const IndexType number_of_nodes = rGeometry.PointsNumber();

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

    int SupportPenalty6pCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenalty6pCondition" << std::endl;
        return 0;
    }

    void SupportPenalty6pCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rResult.size() != 6 * number_of_nodes)
            rResult.resize(6 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 6;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
        }
    }

    void SupportPenalty6pCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
        }
    };

} // Namespace Kratos
