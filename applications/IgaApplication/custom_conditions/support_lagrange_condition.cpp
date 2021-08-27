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
#include "custom_conditions/support_lagrange_condition.h"


namespace Kratos
{
    void SupportLagrangeCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();
        const SizeType mat_size = number_of_non_zero_nodes * 6;

        Matrix LHS = ZeroMatrix(mat_size, mat_size);
        Vector u(mat_size);

        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        // Integration
        const typename GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Determinant of jacobian 
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry, determinant_jacobian_vector);
        }
        else {
            r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        // non zero node counter for Lagrange Multipliers.
        IndexType counter_n = 0;
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            // Differential area, being 1 for points.
            const double integration = (r_geometry.Dimension() == 0)
                ? 1
                : integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            // loop over Lagrange Multipliers
            for (IndexType i = 0; i < number_of_nodes; i++) {
                // non zero node counter for for displacements.
                IndexType counter_m = 0;
                if (r_N(point_number, i) > shape_function_tolerance) {
                    // loop over shape functions of displacements
                    for (IndexType j = 0; j < number_of_nodes; j++) {
                        if (r_N(point_number, j) > shape_function_tolerance) {
                            const double NN = r_N(point_number, j) * r_N(point_number, i) * integration;

                            // indices in local stiffness matrix.
                            const IndexType ibase = counter_n * 3 + 3 * (number_of_non_zero_nodes);
                            const IndexType jbase = counter_m * 3;

                            // Matrix in following shape:
                            // |0 H^T|
                            // |H 0  |

                            LHS(ibase, jbase) = NN;
                            LHS(ibase + 1, jbase + 1) = NN;
                            LHS(ibase + 2, jbase + 2) = NN;

                            LHS(jbase, ibase) = NN;
                            LHS(jbase + 1, ibase + 1) = NN;
                            LHS(jbase + 2, ibase + 2) = NN;

                            counter_m++;
                        }
                    }
                    counter_n++;
                }
            }
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += LHS;
            }

            if (CalculateResidualVectorFlag) {
                const array_1d<double, 3>& displacement = (Has(DISPLACEMENT))
                    ? this->GetValue(DISPLACEMENT)
                    : ZeroVector(3);

                IndexType counter = 0;
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const array_1d<double, 3>& r_disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

                            IndexType index = 3 * counter;
                            u[index] = (r_disp[0] - displacement[0]);
                            u[index + 1] = (r_disp[1] - displacement[1]);
                            u[index + 2] = (r_disp[2] - displacement[2]);
                            counter++;
                        }
                    }
                }
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const array_1d<double, 3>& r_l_m = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                            IndexType index = 3 * (counter);
                            u[index] = r_l_m[0];
                            u[index + 1] = r_l_m[1];
                            u[index + 2] = r_l_m[2];
                            counter++;
                        }
                    }
                }

                noalias(rRightHandSideVector) -= prod(LHS, u);
            }
        }
        KRATOS_CATCH("")
    }

    void SupportLagrangeCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType nb_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType pnt = 0; pnt < nb_integration_points; pnt++)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[pnt];
            J.clear();
            for (IndexType i = 0; i < nb_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += value * r_DN_De(i, m);
                    }
                }
            }

            //Compute the tangent and  the normal to the boundary vector
            array_1d<double, 3> local_tangent;
            GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[pnt] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }

    void SupportLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();

        if (rResult.size() != 6 * number_of_non_zero_nodes)
            rResult.resize(6 * number_of_non_zero_nodes, false);

        IndexType counter = 0;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const IndexType index = counter * 3;
                    const auto& r_node = r_geometry[i];
                    rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                    rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                    rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
                    counter++;
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const IndexType index = 3 * (counter);
                    const auto& r_node = r_geometry[i];
                    rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
                    rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
                    rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
                    counter++;
                }
            }
        }
    }

    void SupportLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * GetNumberOfNonZeroNodes());

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
                }
            }
        }
    }

    std::size_t SupportLagrangeCondition::GetNumberOfNonZeroNodes() const {
        const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType n = 0; n < r_N.size1(); ++n) {
            for (IndexType m = 0; m < r_N.size2(); ++m) {
                if (r_N(n, m) > shape_function_tolerance) {
                    counter++;
                }
            }
        }
        return counter;
    }

} // Namespace Kratos
