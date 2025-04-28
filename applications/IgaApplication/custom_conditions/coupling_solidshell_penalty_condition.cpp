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

// Project includes
#include "custom_conditions/coupling_solidshell_penalty_condition.h"
namespace Kratos
{
    void CouplingSolidShellPenaltyCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

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

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();
                       
            //FOR DISPLACEMENTS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; ++i)
            {

                IndexType index = 3 * i;
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N_master(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N_master(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N_master(point_number, i);
            }

            for (IndexType i = 0; i < number_of_nodes_slave; ++i)
            {
                // Minas: Start: calculate global coordinates of integration point
                GeometryType::CoordinatesArrayType  MasterIPGlobalCoordinates, SlaveIPGlobalCoordinates;
                r_geometry_master.GlobalCoordinates(MasterIPGlobalCoordinates, point_number);
                r_geometry_slave.GlobalCoordinates(SlaveIPGlobalCoordinates, point_number);

                Matrix J;
                r_geometry_slave.Jacobian(J, point_number);

                // Get Gradient of Shape Functions for slave geometry
                const IntegrationMethod integration_method_slave = r_geometry_slave.GetDefaultIntegrationMethod();
                const Matrix& shape_functions_gradients_slave = r_geometry_slave.ShapeFunctionsLocalGradients()[point_number];

                Matrix a3x_g1g3B2_g2g3B1 = zero_matrix(3, 3);
                if (CalculateStiffnessMatrixFlag) {
                    // compute the actual base vectors of slave
                    array_1d<double, 3> g1, g2, g3;

                    g1 = column(J, 0);
                    g2 = column(J, 1);
                    store_g1(g1, i);
                    store_g2(g2, i);               
          
                }

                array_1d<double, 3> g1, g2, g3;
                g1 = _g1[i];
                g2 = _g2[i];

                MathUtils<double>::CrossProduct(g3, g1, g2);
                double norm_g1xg2 = norm_2(g3);
                g3 = g3 / norm_2(g3);

                Matrix g3_cross_product_matrix(3, 3, 0);
                g3_cross_product_matrix(0, 1) = -g3(2);
                g3_cross_product_matrix(1, 0) = g3(2);
                g3_cross_product_matrix(0, 2) = -g3(1);
                g3_cross_product_matrix(2, 0) = g3(1);
                g3_cross_product_matrix(1, 2) = -g3(0);
                g3_cross_product_matrix(2, 1) = g3(0);
               
                Matrix g1g3(3, 3, 0);
                for (int ii = 0; ii < 3; ++ii) {
                    for (int jj = 0; jj < 3; ++jj)
                        g1g3(ii, jj) = g1(ii) * g3(jj);
                }
                //print_matrix_3x3(g1g3, "g1g3 :");
                Matrix g2g3(3, 3, 0);
                for (int ii = 0; ii < 3; ++ii) {
                    for (int jj = 0; jj < 3; ++jj)
                        g2g3(ii, jj) = g2(ii) * g3(jj);
                }

                Matrix B1(3, 3, 0), B2(3, 3, 0);
                B1(0, 0) = shape_functions_gradients_slave(i, 0);
                B1(1, 1) = shape_functions_gradients_slave(i, 0);
                B1(2, 2) = shape_functions_gradients_slave(i, 0);

                B2(0, 0) = shape_functions_gradients_slave(i, 1);
                B2(1, 1) = shape_functions_gradients_slave(i, 1);
                B2(2, 2) = shape_functions_gradients_slave(i, 1);
                
                Matrix g1g3B2 = prod(g1g3, B2);
                Matrix g2g3B1 = prod(g2g3, B1);


                Matrix g1g3B2_g2g3B1 = g1g3B2 - g2g3B1;
                a3x_g1g3B2_g2g3B1 = prod(g3_cross_product_matrix, g1g3B2_g2g3B1) / norm_g1xg2;

                if (CalculateStiffnessMatrixFlag) {
                    // Minas Comppute Theta3
                    double theta3 = (MasterIPGlobalCoordinates(0) - SlaveIPGlobalCoordinates(0)) * g3(0) +
                        (MasterIPGlobalCoordinates(1) - SlaveIPGlobalCoordinates(1)) * g3(1) +
                        (MasterIPGlobalCoordinates(2) - SlaveIPGlobalCoordinates(2)) * g3(2);
                    
                    theta3 = theta3 / norm_2(g3) / norm_2(g3);
                    store_theta3(theta3, i);

                    // Minas that SlavePoint + theta3 * A3 = MasterPoint 
                    double diff71 = SlaveIPGlobalCoordinates(0) + g3(0) * theta3 - MasterIPGlobalCoordinates(0);
                    double diff72 = SlaveIPGlobalCoordinates(1) + g3(1) * theta3 - MasterIPGlobalCoordinates(1);
                    double diff73 = SlaveIPGlobalCoordinates(2) + g3(2) * theta3 - MasterIPGlobalCoordinates(2);
                    bool check71 = abs(diff71) < 1E-10;
                    bool check72 = abs(diff72) < 1E-10;
                    bool check73 = abs(diff73) < 1E-10;

                    // Minas TODO : Import thickness from properties
                    double thickness = 2;
                    KRATOS_ERROR_IF(abs(theta3) > thickness / 2) << "Theta 3 > 1" << std::endl;;
                    KRATOS_ERROR_IF(check71 == false || check72 == false || check73 == false) << "Something happened with theta 3" << std::endl;
                }
                double theta3 = _theta3[i];
                Matrix final_matrix = - a3x_g1g3B2_g2g3B1 * theta3;
                
                /* Following code just for checking at random values
                //if ( (i == 8 || i == 0 || i == 4) && Id() == 550) {
                //    if (CalculateStiffnessMatrixFlag) {
                //        KRATOS_INFO("Condition Id() :  ") << Id() <<  
                //                    " | LHS : Node number =  " << i << 
                //                    " | N_slave(point_number, i) = " << N_slave(point_number, i) << 
                //                    " | shape_functions_gradients_slave(i,0)[point_number] = " << shape_functions_gradients_slave(i,0) << std::endl;
                //                    print_matrix_3x3(final_matrix, "LHS final_matrix :");
                //    }

                //    if (CalculateResidualVectorFlag) {
                //        KRATOS_INFO("Condition Id() :  ") << Id() <<
                //            " | RHS : Node number =  " << i <<
                //            " | N_slave(point_number, i) = " << N_slave(point_number, i) <<
                //            " | shape_functions_gradients_slave(i,0)[point_number] = " << shape_functions_gradients_slave(i, 0) << std::endl;
                //            print_matrix_3x3(final_matrix, "RHS final_matrix :");

                //    }
                //} */

                IndexType index = 3 * (i + number_of_nodes_master);
                if (Is(IgaFlags::FIX_DISPLACEMENT_X)) {
                    H(0, index + 0) = - N_slave(point_number, i) - final_matrix(0, 0);
                    H(0, index + 1) =                            - final_matrix(0, 1);
                    H(0, index + 2) =                            - final_matrix(0, 2);
                } 
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y)) {
                    H(1, index + 0) =                            - final_matrix(1, 0);
                    H(1, index + 1) = - N_slave(point_number, i) - final_matrix(1, 1);
                    H(1, index + 2) =                            - final_matrix(1, 2);
                }
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z)) {
                    H(2, index + 0) =                            - final_matrix(2, 0)   ;
                    H(2, index + 1) =                            - final_matrix(2, 1)   ;
                    H(2, index + 2) = - N_slave(point_number, i) - final_matrix(2, 2)   ;
                }
                
               
            }

            // Differential area
            const double penalty_integration = penalty * (2*integration_points[point_number].Weight()) * determinant_jacobian_vector[point_number];

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += prod(trans(H), H) * penalty_integration;
            }
            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * (i + number_of_nodes_master);
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }

                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u) * penalty_integration;
            }
        }

        KRATOS_CATCH("")
    }

    void CouplingSolidShellPenaltyCondition::store_g1(array_1d<double, 3>& g1, const IndexType& pointnumber_Node) {
        _g1[pointnumber_Node] = g1;
    }

    void CouplingSolidShellPenaltyCondition::store_g2(array_1d<double, 3>& g2, const IndexType& pointnumber_Node) {
        _g2[pointnumber_Node] = g2;
    }
    void CouplingSolidShellPenaltyCondition::store_theta3(double theta3, const IndexType& pointnumber_Node) {
        _theta3[pointnumber_Node] = theta3;
    }

    void CouplingSolidShellPenaltyCondition::print_matrix_3x3(Matrix& rMatrix, std::string nameMatrix) {
        KRATOS_INFO(nameMatrix) << rMatrix(0, 0) << " " << rMatrix(0, 1) << " " << rMatrix(0, 2) << std::endl
            << rMatrix(1, 0) << " " << rMatrix(1, 1) << " " << rMatrix(1, 2) << std::endl
            << rMatrix(2, 0) << " " << rMatrix(2, 1) << " " << rMatrix(2, 2) << std::endl;
    }
    void CouplingSolidShellPenaltyCondition::print_vector_3(array_1d<double, 3>& vector, std::string nameMatrix){
        KRATOS_INFO(nameMatrix) << vector(0) << " | " << vector(1) << " | " << vector(2) << std::endl;
    }

    void CouplingSolidShellPenaltyCondition::DeterminantOfJacobianInitial(
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

    int CouplingSolidShellPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyCondition" << std::endl;
        return 0;
    }

    void CouplingSolidShellPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 3 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(3 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 3 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingSolidShellPenaltyCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


