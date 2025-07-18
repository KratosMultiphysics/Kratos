//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Apostolakis Minas
//                   
//                   
//

// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_solidshell3p_penalty_condition.h"
namespace Kratos
{
    void CouplingSolidShell3pPenaltyCondition::CalculateAll(
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

        const SizeType r_number_of_integration_points_slave = r_geometry_slave.IntegrationPointsNumber();
        if (r_number_of_integration_points_slave != integration_points.size()) {
            std::cout << " # of integration points slave is different than master for a quadrature point" << std::endl;
        }
        if (_g1.size() != r_number_of_integration_points_slave)
            _g1.resize(r_number_of_integration_points_slave);
        if (_g2.size() != r_number_of_integration_points_slave)
            _g2.resize(r_number_of_integration_points_slave);
        if (_g3.size() != r_number_of_integration_points_slave)
            _g3.resize(r_number_of_integration_points_slave);
        if (_theta3.size() != r_number_of_integration_points_slave)
            _theta3.resize(r_number_of_integration_points_slave);

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

            Matrix J;
            r_geometry_slave.Jacobian(J, point_number);
            // Minas: Start: calculate global coordinates of integration point
            GeometryType::CoordinatesArrayType  MasterIPGlobalCoordinates, SlaveIPGlobalCoordinates;
            r_geometry_master.GlobalCoordinates(MasterIPGlobalCoordinates, point_number);
            r_geometry_slave.GlobalCoordinates(SlaveIPGlobalCoordinates, point_number);
            double norm_g1_cross_g2;
            if (CalculateStiffnessMatrixFlag) {
                // compute the actual base vectors of slave
                _g1[point_number] = column(J, 0);
                _g2[point_number] = column(J, 1);

                array_1d<double, 3> g1_cross_g2;
                MathUtils<double>::CrossProduct(g1_cross_g2, _g1[point_number], _g2[point_number]);
                norm_g1_cross_g2 = norm_2(g1_cross_g2);
                _g3[point_number] = g1_cross_g2 / norm_g1_cross_g2;

                // Comppute Theta3
                array_1d<double, 3> distance;
                distance[0] =  MasterIPGlobalCoordinates[0] - SlaveIPGlobalCoordinates[0] - _g3[point_number][0] * _theta3[point_number];
                distance[1] =  MasterIPGlobalCoordinates[1] - SlaveIPGlobalCoordinates[1] - _g3[point_number][1] * _theta3[point_number];
                distance[2] =  MasterIPGlobalCoordinates[2] - SlaveIPGlobalCoordinates[2] - _g3[point_number][2] * _theta3[point_number];
                _theta3[point_number] = inner_prod(distance, _g3[point_number]) ;

                // Check that SlavePoint + theta3 * A3 = MasterPoint 
                array_1d<double, 3> check_distance;
                check_distance = SlaveIPGlobalCoordinates + _g3[point_number] * _theta3[point_number] - MasterIPGlobalCoordinates;
                double check01 = norm_2(distance);
                double check02 = norm_2(check_distance);

                double thickness = GetProperties().GetValue(THICKNESS);
                //if ( abs(_theta3[point_number]) - thickness / 2 > 1E-3) {
                //    std::cout << "Theta3 = " << _theta3[point_number] << "  > thickness / 2" << std::endl;
                //    KRATOS_ERROR;
                //}
                if (norm_2(check_distance) > 1E-6) {
                    array_1d<double, 3> g3 = _g3[point_number];
                    double theta3 = _theta3[point_number];
                    std::cout << "Theta3 caclulation is wrong " << norm_2(distance) <<  std::endl;
                    KRATOS_ERROR;
                }
            }

            Matrix OutOfPlaneDeformationFirstVariationMatrix = zero_matrix(3, 3 * number_of_nodes_slave);
            OutOfPlaneDeformationFirstVariation(
                OutOfPlaneDeformationFirstVariationMatrix,
                3 * number_of_nodes_slave,
                _theta3[point_number],
                _g1[point_number],
                _g2[point_number],
                r_geometry_slave.ShapeFunctionsLocalGradients()[point_number]);

            for (IndexType i = 0; i < number_of_nodes_slave; ++i)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                if (Is(IgaFlags::FIX_DISPLACEMENT_X)) {
                    H(0, index + 0) = -N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i + 0);
                    H(0, index + 1) = -OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i + 1);
                    H(0, index + 2) = -OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i + 2);

                }
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y)) {
                    H(1, index + 0) = -OutOfPlaneDeformationFirstVariationMatrix(1, 3 *  i + 0);
                    H(1, index + 1) = -N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(1, 3 * i + 1);
                    H(1, index + 2) = -OutOfPlaneDeformationFirstVariationMatrix(1, 3 *  i + 2);
                }
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z)) {
                    H(2, index + 0) = -OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i + 0);
                    H(2, index + 1) = -OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i + 1);
                    H(2, index + 2) = -N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i + 2);
                }
            }

            // Differential area
            const double penalty_integration = penalty * (integration_points[point_number].Weight()) * determinant_jacobian_vector[point_number];
            double check_weight = integration_points[point_number].Weight(); 
            double check_det = determinant_jacobian_vector[point_number];

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

    void CouplingSolidShell3pPenaltyCondition::DeterminantOfJacobianInitial(
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

    void CouplingSolidShell3pPenaltyCondition::OutOfPlaneDeformationFirstVariation(
        Matrix& OutOfPlaneDeformationWholeMatrix,
        const size_t& mat_size,
        const double theta3,
        const array_1d<double, 3>& A1,
        const array_1d<double, 3>& A2,
        const Matrix& ShapeFunctionsGradientsValues) {

        for (size_t r = 0; r < mat_size; r++) { // row
            // local node number kr and dof direction dirr
            size_t kr = r / 3;
            size_t dirr = r % 3;

            array_1d<double, 3> N_theta1_r(3, 0), N_theta2_r(3, 0), N_basis_function(3, 0);
            N_theta1_r[dirr] = ShapeFunctionsGradientsValues(kr, 0);
            N_theta2_r[dirr] = ShapeFunctionsGradientsValues(kr, 1);
            array_1d<double, 3> Phi_r_cross_Α3 = Calculate_Phi_r_cross_A3(N_theta1_r, N_theta2_r,A1,A2);

            OutOfPlaneDeformationWholeMatrix(0, r) = theta3 * Phi_r_cross_Α3[0];
            OutOfPlaneDeformationWholeMatrix(1, r) = theta3 * Phi_r_cross_Α3[1];
            OutOfPlaneDeformationWholeMatrix(2, r) = theta3 * Phi_r_cross_Α3[2];
           
        }

    }

    array_1d<double, 3> CouplingSolidShell3pPenaltyCondition::Calculate_Phi_r_cross_A3(
        const array_1d<double, 3>& N_theta1_r, 
        const array_1d<double, 3>& N_theta2_r,
        const array_1d<double, 3>& A1, 
        const array_1d<double, 3>& A2) {
        // A1_cross_A2, A1_cross_A2 norm
        array_1d<double, 3> A1_cross_A2, A3, Phi, Phi_r_cross_A3;
        MathUtils<double>::CrossProduct(A1_cross_A2, A1, A2);
        double norm_A1_cross_A2 = norm_2(A1_cross_A2);
        A3 = A1_cross_A2 / norm_A1_cross_A2;
        
        double phi1 = (1 / norm_A1_cross_A2) * inner_prod(N_theta2_r, A3);
        double phi2 = - (1 / norm_A1_cross_A2) * inner_prod(N_theta1_r, A3);
        Phi = phi1 * A1 + phi2 * A2;

        MathUtils<double>::CrossProduct(Phi_r_cross_A3, Phi, A3);
        return Phi_r_cross_A3;
    }


    int CouplingSolidShell3pPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyCondition" << std::endl;
        return 0;
    }

    void CouplingSolidShell3pPenaltyCondition::EquationIdVector(
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

    void CouplingSolidShell3pPenaltyCondition::GetDofList(
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


