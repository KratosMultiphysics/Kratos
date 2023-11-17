//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "custom_conditions/sbm_support_lagrange_condition.h"


namespace Kratos
{
    void SBMSupportLagrangeCondition::CalculateAll(
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
        const SizeType mat_size = number_of_non_zero_nodes * 2;

        Matrix LHS = ZeroMatrix(mat_size, mat_size);
        Vector u(mat_size);

        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        // Integration
        const typename GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Determinant of jacobian
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);


        // non zero node counter for Lagrange Multipliers.
        IndexType counter_n = 0;
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            // Differential area, being 1 for points.
            const double integration = integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

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
                            const IndexType ibase = counter_n * 1 + 1 * (number_of_non_zero_nodes);
                            const IndexType jbase = counter_m * 1;

                            // Matrix in following shape:
                            // |0 H^T|
                            // |H 0  |

                            LHS(ibase, jbase) = NN;

                            LHS(jbase, ibase) = NN;

                            counter_m++;
                        }
                    }
                    counter_n++;
                }
            }
            noalias(rLeftHandSideMatrix) += LHS;

            if (CalculateResidualVectorFlag) {
                const array_1d<double, 1>& displacement = ZeroVector(1);

                IndexType counter = 0;
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const double r_disp = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);

                            u[counter] = (r_disp - displacement[0]);
                            counter++;
                        }
                    }
                }
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const double r_l_m = r_geometry[i].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
                            u[counter] = r_l_m;
                            counter++;
                        }
                    }
                }

                noalias(rRightHandSideVector) -= prod(LHS, u);
            }
        }
        KRATOS_CATCH("")
    }

    void SBMSupportLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();

        if (rResult.size() != 2 * number_of_non_zero_nodes)
            rResult.resize(2 * number_of_non_zero_nodes, false);

        IndexType counter = 0;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rResult[counter]     = r_node.GetDof(TEMPERATURE).EquationId();
                    counter++;
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rResult[counter]     = r_node.GetDof(SCALAR_LAGRANGE_MULTIPLIER).EquationId();
                    counter++;
                }
            }
        }
    }

    void SBMSupportLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * GetNumberOfNonZeroNodes());

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(SCALAR_LAGRANGE_MULTIPLIER));
                }
            }
        }
    }

    std::size_t SBMSupportLagrangeCondition::GetNumberOfNonZeroNodes() const {
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
