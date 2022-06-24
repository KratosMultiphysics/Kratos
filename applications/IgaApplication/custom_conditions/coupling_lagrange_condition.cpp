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
//                   Tobias Tescheamacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_lagrange_condition.h"

namespace Kratos
{

    void CouplingLagrangeCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType number_of_non_zero_nodes_master = GetNumberOfNonZeroNodesMaster();
        const SizeType number_of_non_zero_nodes_slave = GetNumberOfNonZeroNodesSlave();

        const SizeType mat_size = 6 * number_of_non_zero_nodes_master + 3 * number_of_non_zero_nodes_slave;

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

        Matrix LHS = ZeroMatrix(mat_size, mat_size);
        Vector u(mat_size);

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector);
        }
        else {
            r_geometry_master.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        // non zero node counter for Lagrange Multipliers.
        IndexType counter_n = 0;
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            const Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Vector H;
            H.resize(number_of_nodes_master + number_of_nodes_slave);

            Vector H_lambda;
            H_lambda.resize(number_of_nodes_master);
            H_lambda = ZeroVector(number_of_nodes_master);

            for (IndexType i = 0; i < number_of_nodes_master; ++i)
            {
                H[i] = N_master(point_number, i);
                H_lambda[i] = N_master(point_number, i);
            }

            for (IndexType i = 0; i < number_of_nodes_slave; ++i)
            {
                H[i + number_of_nodes_master] = -N_slave(point_number, i);
            }

            // Differential area
            const double integration = integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            // loop over Lagrange Multipliers
            for (IndexType i = 0; i < number_of_nodes_master; ++i)
	        {
                // non zero node counter for for displacements.
                IndexType counter_m = 0;
                if (H_lambda[i] > shape_function_tolerance) {
                    // loop over shape functions of displacements
                    for (IndexType j = 0; j < number_of_nodes_master + number_of_nodes_slave; ++j) 
                    {
                        if (std::abs(H[j]) > shape_function_tolerance) {
                            const double HH = H[j] * H_lambda[i];
                                
                            const IndexType ibase = counter_n * 3 + 3 * (number_of_non_zero_nodes_master + number_of_non_zero_nodes_slave);
                            const IndexType jbase = counter_m * 3;

                            // Matrix in following shape:
                            // |0 H^T|
                            // |H 0  |

                            LHS(ibase,     jbase)     = HH;
                            LHS(ibase + 1, jbase + 1) = HH;
                            LHS(ibase + 2, jbase + 2) = HH;

                            LHS(jbase,     ibase)     = HH;
                            LHS(jbase + 1, ibase + 1) = HH;
                            LHS(jbase + 2, ibase + 2) = HH;

                            counter_m++;
                        }
                    }
                counter_n++;
                }
            }

            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += LHS * integration;
            }

            if (CalculateResidualVectorFlag) {

                IndexType counter = 0;
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    for (IndexType n = 0; n < N_master.size1(); ++n) 
                    {
                        if (N_master(n, i) > shape_function_tolerance) 
                        {
                            const array_1d<double, 3> r_disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                            IndexType index = 3 * counter;
                            u[index]     = r_disp[0];
                            u[index + 1] = r_disp[1];
                            u[index + 2] = r_disp[2];
                            counter++;
                        }
                    }
                }
                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    for (IndexType n = 0; n < N_slave.size1(); ++n) 
                    {
                        if (N_slave(n, i) > shape_function_tolerance) 
                        {
                            const array_1d<double, 3> r_disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                            IndexType index = 3 * (counter);
                            u[index]     = r_disp[0];
                            u[index + 1] = r_disp[1];
                            u[index + 2] = r_disp[2];
                            counter++;
                        }
                    }
                }
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    for (IndexType n = 0; n < N_master.size1(); ++n) 
                    {
                        if (N_master(n, i) > shape_function_tolerance) 
                        {
                            const array_1d<double, 3> r_l_m = r_geometry_master[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                            IndexType index = 3 * (counter);
                            u[index]     = r_l_m[0];
                            u[index + 1] = r_l_m[1];
                            u[index + 2] = r_l_m[2];
                            counter++;
                        }
                    }
                }

                noalias(rRightHandSideVector) -= prod(LHS, u) * integration;
            }
        }

        KRATOS_CATCH("")
    }

    void CouplingLagrangeCondition::DeterminantOfJacobianInitial(
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

    void CouplingLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const Matrix N_master = r_geometry_master.ShapeFunctionsValues();
        const Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType number_of_non_zero_nodes_master = GetNumberOfNonZeroNodesMaster();
        const SizeType number_of_non_zero_nodes_slave = GetNumberOfNonZeroNodesSlave();

        if (rResult.size() != 6 * number_of_non_zero_nodes_master + 3* number_of_non_zero_nodes_slave)
            rResult.resize(6 * number_of_non_zero_nodes_master + 3* number_of_non_zero_nodes_slave, false);

        IndexType counter = 0;
        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType n = 0; n < N_master.size1(); ++n) {
                if (N_master(n, i) > shape_function_tolerance) {
                    const IndexType index = counter * 3;
                    const auto& r_node = r_geometry_master[i];
                    rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                    rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                    rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
                    counter++;
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            for (IndexType n = 0; n < N_slave.size1(); ++n) {
                if (N_slave(n, i) > shape_function_tolerance) {
                    const IndexType index = 3 * (counter);
                    const auto& r_node = r_geometry_slave[i];
                    rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                    rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                    rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
                    counter++;
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType n = 0; n < N_master.size1(); ++n) {
                if (N_master(n, i) > shape_function_tolerance) {
                    const IndexType index = 3 * (counter);
                    const auto& r_node = r_geometry_master[i];
                    rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
                    rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
                    rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
                    counter++;
                }
            }
        }

        KRATOS_CATCH("")
    }

    void CouplingLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const Matrix N_master = r_geometry_master.ShapeFunctionsValues();
        const Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * GetNumberOfNonZeroNodesMaster() + 3* GetNumberOfNonZeroNodesSlave());

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType n = 0; n < N_master.size1(); ++n) {
                if (N_master(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry_master.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            for (IndexType n = 0; n < N_slave.size1(); ++n) {
                if (N_slave(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry_slave.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
                    rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            for (IndexType n = 0; n < N_master.size1(); ++n) {
                if (N_master(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry_master.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
                    rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
                }
            }
        }

        KRATOS_CATCH("")
    }

    std::size_t CouplingLagrangeCondition::GetNumberOfNonZeroNodesMaster() const {
        const Matrix N_master = GetGeometry().GetGeometryPart(0).ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType n = 0; n < N_master.size1(); ++n) {
            for (IndexType m = 0; m < N_master.size2(); ++m) {
                if (N_master(n, m) > shape_function_tolerance) {
                    counter++;
                }
            }
        }
        return counter;
    }

    std::size_t CouplingLagrangeCondition::GetNumberOfNonZeroNodesSlave() const {
        const Matrix N_slave = GetGeometry().GetGeometryPart(1).ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType n = 0; n < N_slave.size1(); ++n) {
            for (IndexType m = 0; m < N_slave.size2(); ++m) {
                if (N_slave(n, m) > shape_function_tolerance) {
                    counter++;
                }
            }
        }
        return counter;
    }

} // Namespace Kratos


