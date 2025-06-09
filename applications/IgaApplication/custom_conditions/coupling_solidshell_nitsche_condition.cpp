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
#include "custom_conditions/coupling_solidshell_nitsche_condition.h"

namespace Kratos
{
    ///@name Initialize Functions
    ///@{
    /// 
    
    void CouplingSolidShellNitscheCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY;

       const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_FACTOR];

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

        // initial determinant of jacobian 
        Vector determinant_jacobian_vector_initial(integration_points.size());
        //DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector_initial); 
        r_geometry_master.DeterminantOfJacobian(determinant_jacobian_vector_initial);

        const IntegrationMethod integration_method_master = r_geometry_master.GetDefaultIntegrationMethod();
        const IntegrationMethod integration_method_slave = r_geometry_slave.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = r_geometry_master.ShapeFunctionsLocalGradients(integration_method_master);
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = r_geometry_slave.ShapeFunctionsLocalGradients(integration_method_slave);

        const SizeType r_number_of_integration_points_master = r_geometry_master.IntegrationPointsNumber();
        const SizeType r_number_of_integration_points_slave = r_geometry_slave.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_ab_covariant_vector_master.size() != r_number_of_integration_points_master)
            m_A_ab_covariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_A_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_A_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_B_ab_covariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_B_ab_covariant_vector_slave.resize(r_number_of_integration_points_slave);
        //if (m_dA_vector_master.size() != r_number_of_integration_points_master)
        //    m_dA_vector_master.resize(r_number_of_integration_points_master);
        if (m_dA_vector_slave.size() != r_number_of_integration_points_slave)
            m_dA_vector_slave.resize(r_number_of_integration_points_slave);
        //if (m_T_vector_master.size() != r_number_of_integration_points_master) // Transformations for master are not needed
        //    m_T_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_vector_slave.resize(r_number_of_integration_points_slave);
        //if (m_T_hat_vector_master.size() != r_number_of_integration_points_master) // Transformations for master are not needed
        //    m_T_hat_vector_master.resize(r_number_of_integration_points_master);
        if (m_T_hat_vector_slave.size() != r_number_of_integration_points_slave)
            m_T_hat_vector_slave.resize(r_number_of_integration_points_slave);
        if (m_reference_contravariant_base_master.size() != r_number_of_integration_points_master)
            m_reference_contravariant_base_master.resize(r_number_of_integration_points_master);
        if (m_reference_contravariant_base_slave.size() != r_number_of_integration_points_slave)
            m_reference_contravariant_base_slave.resize(r_number_of_integration_points_slave);
        if (m_n_contravariant_vector_master.size() != r_number_of_integration_points_master)
            m_n_contravariant_vector_master.resize(r_number_of_integration_points_master);
        if (m_n_contravariant_vector_slave.size() != r_number_of_integration_points_slave)
            m_n_contravariant_vector_slave.resize(r_number_of_integration_points_slave);
        if (_theta3.size() != r_number_of_integration_points_slave)
            _theta3.resize(r_number_of_integration_points_slave);
        if (_A1.size() != r_number_of_integration_points_slave)
            _A1.resize(r_number_of_integration_points_slave);
        if (_A2.size() != r_number_of_integration_points_slave)
            _A2.resize(r_number_of_integration_points_slave);
        if (_A3.size() != r_number_of_integration_points_slave)
            _A3.resize(r_number_of_integration_points_slave);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i_master = r_shape_functions_gradients_master[point_number];
            const Matrix& shape_functions_gradients_i_slave = r_shape_functions_gradients_slave[point_number];

            //Compute Kinematics Reference
            if (CalculateStiffnessMatrixFlag) {
                KinematicVariablesSolid kinematic_variables_reference_master(
                    r_geometry_master.WorkingSpaceDimension());

                KinematicVariablesShell kinematic_variables_reference_slave(
                    r_geometry_slave.WorkingSpaceDimension());

                CalculateKinematicsSolid(
                    point_number,
                    kinematic_variables_reference_master, shape_functions_gradients_i_master, ConfigurationType::Reference); // MasterPart is Solid     

                CalculateKinematicsShell(
                    point_number,
                    kinematic_variables_reference_slave, shape_functions_gradients_i_slave, ConfigurationType::Reference);

                m_A_ab_covariant_vector_master[point_number] = kinematic_variables_reference_master.a_ab_covariant;
                m_A_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.a_ab_covariant;
                m_B_ab_covariant_vector_slave[point_number] = kinematic_variables_reference_slave.b_ab_covariant;
                //m_dA_vector_master[point_number] = kinematic_variables_reference_master.dA; // m_dA_vector_master not used for solid
                m_dA_vector_slave[point_number] = kinematic_variables_reference_slave.dA;

                m_n_contravariant_vector_master[point_number] = kinematic_variables_reference_master.n_contravariant;
                m_n_contravariant_vector_slave[point_number] = kinematic_variables_reference_slave.n_contravariant;

                CalculateTransformationShell(kinematic_variables_reference_slave, m_T_vector_slave[point_number],
                    m_T_hat_vector_slave[point_number], m_reference_contravariant_base_slave[point_number]); // m_reference_contravariant_base_slave is not used anywhere
            
                _A1[point_number] = kinematic_variables_reference_slave.a1;
                _A2[point_number] = kinematic_variables_reference_slave.a2;
                _A3[point_number] = kinematic_variables_reference_slave.a3;

                _theta3[point_number] = ComputeTheta3Shell(point_number, kinematic_variables_reference_slave);
            }

            // Compute Kinematics Current
            KinematicVariablesSolid kinematic_variables_master(
                r_geometry_master.WorkingSpaceDimension());
            KinematicVariablesShell kinematic_variables_slave(
                r_geometry_slave.WorkingSpaceDimension());
            CalculateKinematicsSolid(
                point_number,
                kinematic_variables_master, shape_functions_gradients_i_master, ConfigurationType::Current);
            CalculateKinematicsShell(
                point_number,
                kinematic_variables_slave, shape_functions_gradients_i_slave, ConfigurationType::Current);


            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
                r_geometry_slave, GetProperties().GetSubProperties().back(), rCurrentProcessInfo);
            ConstitutiveLaw::Parameters constitutive_law_parameters_master(
                r_geometry_master, GetProperties().GetSubProperties().front(), rCurrentProcessInfo); 

            ConstitutiveVariables constitutive_variables_solid_master(6);
            ConstitutiveVariables constitutive_variables_slave(3);

            CalculateConstitutiveVariablesSolid(
                point_number,
                kinematic_variables_master,
                constitutive_variables_solid_master,
                constitutive_law_parameters_master,
                ConstitutiveLaw::StressMeasure_PK2);

            CalculateConstitutiveVariablesShell(
                point_number,
                _theta3[point_number],
                kinematic_variables_slave,
                constitutive_variables_slave,
                constitutive_law_parameters_slave,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate traction vectors
            array_1d<double, 3> traction_vector_master;
            array_1d<double, 3> traction_vector_slave;

            // CalculateTractionSolid
            CalculateTractionSolid(point_number, traction_vector_master, kinematic_variables_master, constitutive_variables_solid_master);
            CalculateTractionShell(point_number, traction_vector_slave, kinematic_variables_slave, constitutive_variables_slave);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix first_variations_stress_covariant_master = ZeroMatrix(6, 3 * number_of_nodes_master);
            Matrix first_variations_stress_covariant_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            // CalculateFirstVariationStressCovariantSolid
            CalculateFirstVariationStressCovariantSolid(point_number, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_solid_master); 
            CalculateFirstVariationStressCovariantShell(point_number, _theta3[point_number], first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_slave);

            // calculate first variation of traction vectors
            Matrix first_variations_traction_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix first_variations_traction_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);

            // CalculateFirstVariationTraction
            CalculateFirstVariationTractionSolid(point_number, first_variations_traction_master, first_variations_stress_covariant_master, kinematic_variables_master, constitutive_variables_solid_master);
            CalculateFirstVariationTractionShell(point_number, first_variations_traction_slave, first_variations_stress_covariant_slave, kinematic_variables_slave, constitutive_variables_slave);

            Matrix first_variations_traction = ZeroMatrix(3, mat_size);
            for (SizeType i = 0; i < 3 * number_of_nodes_master; ++i) {
                first_variations_traction(0, i) = first_variations_traction_master(0, i);
                first_variations_traction(1, i) = first_variations_traction_master(1, i);
                first_variations_traction(2, i) = first_variations_traction_master(2, i);
            }
            for (SizeType i = 0; i < 3 * number_of_nodes_slave; ++i) {
                first_variations_traction(0, i + 3 * number_of_nodes_master) = - first_variations_traction_slave(0, i);
                first_variations_traction(1, i + 3 * number_of_nodes_master) = - first_variations_traction_slave(1, i);
                first_variations_traction(2, i + 3 * number_of_nodes_master) = - first_variations_traction_slave(2, i);
            }

            ////Compute the NURBS basis functions
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            Matrix r_N_master = ZeroMatrix(3, 3 * number_of_nodes_master);
            Matrix r_N_slave = ZeroMatrix(3, 3 * number_of_nodes_slave);
            Matrix dN_theta1 = ZeroMatrix(3, 3 * number_of_nodes_slave);
            Matrix dN_theta2 = ZeroMatrix(3, 3 * number_of_nodes_slave);


            for (IndexType r = 0; r < number_of_nodes_master; r++)
            {
                r_N_master(0, 3 * r)     =  N_master(point_number, r);
                r_N_master(1, 3 * r + 1) =  N_master(point_number, r);
                r_N_master(2, 3 * r + 2) =  N_master(point_number, r);
            }

            for (IndexType r = 0; r < number_of_nodes_slave; r++)
            {
                r_N_slave(0, 3 * r)     =  N_slave(point_number, r);
                r_N_slave(1, 3 * r + 1) =  N_slave(point_number, r);
                r_N_slave(2, 3 * r + 2) =  N_slave(point_number, r);

                dN_theta1(0, 3 * r)     = shape_functions_gradients_i_slave(r, 0);
                dN_theta1(1, 3 * r + 1) = shape_functions_gradients_i_slave(r, 0);
                dN_theta1(2, 3 * r + 2) = shape_functions_gradients_i_slave(r, 0);

                dN_theta2(0, 3 * r)     = shape_functions_gradients_i_slave(r, 1);
                dN_theta2(1, 3 * r + 1) = shape_functions_gradients_i_slave(r, 1);
                dN_theta2(2, 3 * r + 2) = shape_functions_gradients_i_slave(r, 1);

            }

            //Get the displacement vectors of the previous iteration step
            Vector current_displacement_total = ZeroVector(mat_size);
            Vector current_displacement_master = ZeroVector(3 * number_of_nodes_master);
            Vector current_displacement_slave = ZeroVector(3 * number_of_nodes_slave);

            GetValuesVector(current_displacement_total);

            for (SizeType i = 0; i < 3 * number_of_nodes_master; ++i) {
                current_displacement_master[i] = current_displacement_total[i];
            }
            for (SizeType i = 0; i < 3 * number_of_nodes_slave; ++i) {
                current_displacement_slave[i] = current_displacement_total[i + 3 * number_of_nodes_master];
            }

            array_1d<double, 3> displacement_vector_master;
            array_1d<double, 3> displacement_vector_slave;

            displacement_vector_master = prod(r_N_master, current_displacement_master);

            // For shell displacements, out of plane deformation must be taken into account            
            array_1d<double, 3> A1_cross_A2;
            MathUtils<double>::CrossProduct(A1_cross_A2, _A1[point_number], _A2[point_number]);
            double norm_A1_cross_A2 = norm_2(A1_cross_A2);

            array_1d<double, 3> v_1 = prod(dN_theta1, current_displacement_slave);
            array_1d<double, 3> v_2 = prod(dN_theta2, current_displacement_slave);
            double phi1 = (1 / norm_A1_cross_A2) * inner_prod(v_2, _A3[point_number]);
            double phi2 = - (1 / norm_A1_cross_A2) * inner_prod(v_1, _A3[point_number]);

            array_1d<double, 3> Phi  = phi1 * _A1[point_number] + phi2 * _A2[point_number];
            array_1d<double, 3> Phi_cross_A3;
            MathUtils<double>::CrossProduct(Phi_cross_A3, Phi, _A3[point_number]);

            displacement_vector_slave = prod(r_N_slave, current_displacement_slave) + _theta3[point_number] * Phi_cross_A3;

            // calculate second variation of traction vectors
            Matrix product_second_variations_traction_displacement_master = ZeroMatrix(3 * number_of_nodes_master, 3 * number_of_nodes_master);
            Matrix second_variations_traction_slave = ZeroMatrix(3 * number_of_nodes_slave, 3 * number_of_nodes_slave);
            
            CalculateSecondVariationTractionSolid(
                point_number,
                product_second_variations_traction_displacement_master,
                first_variations_stress_covariant_master,
                displacement_vector_master,
                displacement_vector_slave,
                kinematic_variables_master,
                constitutive_variables_solid_master);

            CalculateSecondVariationTractionShell(
                point_number,
                second_variations_traction_slave,
                _theta3[point_number],
                first_variations_stress_covariant_slave,
                displacement_vector_master,
                displacement_vector_slave,
                kinematic_variables_slave,
                constitutive_variables_slave);

            //Penalty part & RHS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; i++)
            {
                IndexType index = 3 * i;
                H(0, index) = N_master(point_number, i);
                H(1, index + 1) = N_master(point_number, i);
                H(2, index + 2) = N_master(point_number, i);
            }

            Matrix OutOfPlaneDeformationFirstVariationMatrix = zero_matrix(3, 3 * number_of_nodes_slave);
            OutOfPlaneDeformationFirstVariation(
                OutOfPlaneDeformationFirstVariationMatrix,
                3 * number_of_nodes_slave,
                _theta3[point_number],
                _A1[point_number],
                _A2[point_number],
                r_geometry_slave.ShapeFunctionsLocalGradients()[point_number]);

            for (IndexType i = 0; i < number_of_nodes_slave; i++)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                H(0, index + 0) = - N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i)    ;
                H(0, index + 1) = - OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i + 1);
                H(0, index + 2) = - OutOfPlaneDeformationFirstVariationMatrix(0, 3 * i + 2);
                                                                                                               
                H(1, index + 0) = - OutOfPlaneDeformationFirstVariationMatrix(1, 3 * i)    ;
                H(1, index + 1) = - N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(1, 3 * i + 1);
                H(1, index + 2) = - OutOfPlaneDeformationFirstVariationMatrix(1, 3 * i + 2);
                                                                                                               
                H(2, index + 0) = - OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i)    ;
                H(2, index + 1) = - OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i + 1);
                H(2, index + 2) = - N_slave(point_number, i) - OutOfPlaneDeformationFirstVariationMatrix(2, 3 * i + 2);
            }

            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 0.5;

            // Assembly
            if (CalculateStiffnessMatrixFlag) {

                noalias(rLeftHandSideMatrix) += (prod(trans(first_variations_traction), H) + prod(trans(H), first_variations_traction))
                    * integration_weight * determinant_jacobian * -gammaTilde;

                for (IndexType i = 0; i < 3 * number_of_nodes_master; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_master; j++)
                    {
                        rLeftHandSideMatrix(i, j) += product_second_variations_traction_displacement_master(i, j) 
                            * integration_weight * determinant_jacobian * -gammaTilde;
                    }
                }

                for (IndexType i = 0; i < 3 * number_of_nodes_slave; i++)
                {
                    for (IndexType j = 0; j < 3 * number_of_nodes_slave; j++)
                    {
                        rLeftHandSideMatrix(i + 3 * number_of_nodes_master, j + 3 * number_of_nodes_master) += second_variations_traction_slave(i, j) 
                            * integration_weight * determinant_jacobian * -gammaTilde;
                    }
                }

                noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                    * integration_weight * determinant_jacobian * stabilization_parameter;
            }

            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index] = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (IndexType i = 0; i < number_of_nodes_slave; i++)
                {
                    const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * (i + number_of_nodes_master);
                    u[index] = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }

                noalias(rRightHandSideVector) -= (prod(trans(H), traction_vector_master) - prod(trans(H), traction_vector_slave))
                    * integration_weight * determinant_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= (prod(trans(first_variations_traction), displacement_vector_master) - prod(trans(first_variations_traction), displacement_vector_slave))
                    * integration_weight * determinant_jacobian * -gammaTilde;
                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                    * integration_weight * determinant_jacobian * stabilization_parameter;
            }
        }
        KRATOS_CATCH("")
    }

    void CouplingSolidShellNitscheCondition::CalculateKinematicsShell(
        IndexType IntegrationPointIndex,
        KinematicVariablesShell& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration
    )
    {
        KRATOS_TRY; 
        const auto& r_geometry = GetGeometry().GetGeometryPart(1);

        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType number_of_nodes = r_geometry.size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        Vector current_displacement_total = ZeroVector(dimension*(GetGeometry().GetGeometryPart(0).size()+GetGeometry().GetGeometryPart(1).size()));
        Vector current_displacement = ZeroVector(dimension*number_of_nodes);

        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement_total);

        for (SizeType i=0;i<dimension*number_of_nodes;++i){
            current_displacement[i] = current_displacement_total[i+GetGeometry().GetGeometryPart(0).size()*3];
        }

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 0);
            g1[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
            g1[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

            g2[0] += (r_geometry.GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 1);
            g2[1] += (r_geometry.GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
            g2[2] += (r_geometry.GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
        }

        rKinematicVariables.a1 = g1;
        rKinematicVariables.a2 = g2;

        //not-normalized base vector 3
        MathUtils<double>::CrossProduct(rKinematicVariables.a3_tilde, rKinematicVariables.a1, rKinematicVariables.a2);

        //differential area dA
        rKinematicVariables.dA = norm_2(rKinematicVariables.a3_tilde);

        //base vector 3 normalized
        noalias(rKinematicVariables.a3) = rKinematicVariables.a3_tilde / rKinematicVariables.dA;

        //GetCovariantMetric
        rKinematicVariables.a_ab_covariant[0] = pow(rKinematicVariables.a1[0], 2) + pow(rKinematicVariables.a1[1], 2) + pow(rKinematicVariables.a1[2], 2);
        rKinematicVariables.a_ab_covariant[1] = pow(rKinematicVariables.a2[0], 2) + pow(rKinematicVariables.a2[1], 2) + pow(rKinematicVariables.a2[2], 2);
        rKinematicVariables.a_ab_covariant[2] = rKinematicVariables.a1[0] * rKinematicVariables.a2[0] + rKinematicVariables.a1[1] * rKinematicVariables.a2[1] + rKinematicVariables.a1[2] * rKinematicVariables.a2[2];

        //Compute the tangent and  the normal to the boundary vector
        array_1d<double, 3> local_tangent;
        r_geometry.Calculate(LOCAL_TANGENT, local_tangent);

        rKinematicVariables.t = local_tangent[0]*g1 + local_tangent[1]*g2;
        MathUtils<double>::CrossProduct(rKinematicVariables.n, rKinematicVariables.t/norm_2(rKinematicVariables.t), rKinematicVariables.a3);

        // transform the normal into the contavariant basis
        rKinematicVariables.n_contravariant[0] = rKinematicVariables.a1[0]*rKinematicVariables.n[0] + rKinematicVariables.a1[1]*rKinematicVariables.n[1] + rKinematicVariables.a1[2]*rKinematicVariables.n[2];
        rKinematicVariables.n_contravariant[1] = rKinematicVariables.a2[0]*rKinematicVariables.n[0] + rKinematicVariables.a2[1]*rKinematicVariables.n[1] + rKinematicVariables.a2[2]*rKinematicVariables.n[2];
        
        // Kindl From eq 3.28 calculates dda_dαdβ

        // Hesian contatins [a11 | a22 | a12], a11 = da1/dtheta1, a22 = da2/dtheta2, a12 = da1 / dtheta2
        // For a21 : a21 = a12 because a21 = da2/dtheta1 = dd x / (dtheta1 dtheta2) = dtheta1
        
        CalculateHessian(rKinematicVariables.Hessian, r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex));
        Matrix& H = rKinematicVariables.Hessian;
               
        // Kindl Eq 3.28 : b_αβ 
        rKinematicVariables.b_ab_covariant[0] = H(0, 0) * rKinematicVariables.a3[0] + H(1, 0) * rKinematicVariables.a3[1] + H(2, 0) * rKinematicVariables.a3[2];
        rKinematicVariables.b_ab_covariant[1] = H(0, 1) * rKinematicVariables.a3[0] + H(1, 1) * rKinematicVariables.a3[1] + H(2, 1) * rKinematicVariables.a3[2];
        rKinematicVariables.b_ab_covariant[2] = H(0, 2) * rKinematicVariables.a3[0] + H(1, 2) * rKinematicVariables.a3[1] + H(2, 2) * rKinematicVariables.a3[2];

        // Check Numbers
        if (Id() == 550) {
            if (rConfiguration == ConfigurationType::Reference) {
                std::cout << "CalculateKinematicsShell :: Shell :: Reference Congiguration" << std::endl;
                std::cout << "Contravatiants: n1  = ( " << rKinematicVariables.n_contravariant[0] << " | n2:  " << rKinematicVariables.n_contravariant[1] << std::endl;
            }
            else if (rConfiguration == ConfigurationType::Current)
            {
                std::cout << "CalculateKinematicsShell :: Shell :: Current Congiguration" << std::endl;
            }
            std::cout << "A1 = ( " << g1[0] << " | " << g1[1] << " | " << g1[2] << " )" << std::endl;
            std::cout << "A2 = ( " << g2[0] << " | " << g2[1] << " | " << g2[2] << " )" << std::endl;
        }
        KRATOS_CATCH("")
    }

    void CouplingSolidShellNitscheCondition::CalculateKinematicsSolid( 
        IndexType IntegrationPointIndex,
        KinematicVariablesSolid& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration
    )
    {
        KRATOS_TRY; 

        IndexType GeometryPart =  0 ;
        const auto& r_geometry = GetGeometry().GetGeometryPart(GeometryPart);

        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType number_of_nodes = r_geometry.size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);
        Vector g3 = ZeroVector(dimension);

        Vector current_displacement_total = ZeroVector(dimension * (GetGeometry().GetGeometryPart(0).size() + GetGeometry().GetGeometryPart(1).size()));
        Vector current_displacement = ZeroVector(dimension * number_of_nodes);

        if (rConfiguration == ConfigurationType::Current) GetValuesVector(current_displacement_total);

        for (SizeType i = 0; i < dimension * number_of_nodes; ++i) {
            current_displacement[i] = current_displacement_total[i];
        }
       
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            g1[0] += (r_geometry.GetPoint(i).X0() + current_displacement[i * dimension]) * rShapeFunctionGradientValues(i, 0);
            g1[1] += (r_geometry.GetPoint(i).Y0() + current_displacement[(i * dimension) + 1]) * rShapeFunctionGradientValues(i, 0);
            g1[2] += (r_geometry.GetPoint(i).Z0() + current_displacement[(i * dimension) + 2]) * rShapeFunctionGradientValues(i, 0);

            g2[0] += (r_geometry.GetPoint(i).X0() + current_displacement[i * dimension]) * rShapeFunctionGradientValues(i, 1);
            g2[1] += (r_geometry.GetPoint(i).Y0() + current_displacement[(i * dimension) + 1]) * rShapeFunctionGradientValues(i, 1);
            g2[2] += (r_geometry.GetPoint(i).Z0() + current_displacement[(i * dimension) + 2]) * rShapeFunctionGradientValues(i, 1);

            g3[0] += (r_geometry.GetPoint(i).X0() + current_displacement[i * dimension]) * rShapeFunctionGradientValues(i, 2);
            g3[1] += (r_geometry.GetPoint(i).Y0() + current_displacement[(i * dimension) + 1]) * rShapeFunctionGradientValues(i, 2);
            g3[2] += (r_geometry.GetPoint(i).Z0() + current_displacement[(i * dimension) + 2]) * rShapeFunctionGradientValues(i, 2);
        }

        rKinematicVariables.a1 = g1;
        rKinematicVariables.a2 = g2;
        rKinematicVariables.a3 = g3;

        //GetCovariantMetric
   
        rKinematicVariables.a_ab_covariant[0] = inner_prod(rKinematicVariables.a1, rKinematicVariables.a1); // g11
        rKinematicVariables.a_ab_covariant[1] = inner_prod(rKinematicVariables.a2, rKinematicVariables.a2); // g22
        rKinematicVariables.a_ab_covariant[2] = inner_prod(rKinematicVariables.a3, rKinematicVariables.a3); // g33

        rKinematicVariables.a_ab_covariant[3] = inner_prod(rKinematicVariables.a2, rKinematicVariables.a3); // g23
        rKinematicVariables.a_ab_covariant[4] = inner_prod(rKinematicVariables.a1, rKinematicVariables.a3); // g13
        rKinematicVariables.a_ab_covariant[5] = inner_prod(rKinematicVariables.a1, rKinematicVariables.a2); // g12

        // TODO Compute the normal to the boundary vector
        // It should be imported it from Queso directly
        // For this example the n = (1,0,0) is done manually
        rKinematicVariables.n[0] = 1;
        rKinematicVariables.n[1] = 0;
        rKinematicVariables.n[2] = 0;

        // transform the normal into the contavariant basis
        Matrix beta = ZeroMatrix(3, 3);
        Matrix transpose_inv_beta = ZeroMatrix(3, 3);

        for (size_t ii = 0; ii < 3; ii++) {
            beta(0, ii) = g1[ii];
            beta(1, ii) = g2[ii];
            beta(2, ii) = g3[ii];
        }
       
        CalculateTransposeInverseMatrix3x3(transpose_inv_beta, beta);
        array_1d<double, 3> n_covariant;
        n_covariant[0] = transpose_inv_beta(0, 0) * rKinematicVariables.n[0] + transpose_inv_beta(0, 1) * rKinematicVariables.n[1] + transpose_inv_beta(0, 2) * rKinematicVariables.n[2];
        n_covariant[1] = transpose_inv_beta(1, 0) * rKinematicVariables.n[0] + transpose_inv_beta(1, 1) * rKinematicVariables.n[1] + transpose_inv_beta(1, 2) * rKinematicVariables.n[2];
        n_covariant[2] = transpose_inv_beta(2, 0) * rKinematicVariables.n[0] + transpose_inv_beta(2, 1) * rKinematicVariables.n[1] + transpose_inv_beta(2, 2) * rKinematicVariables.n[2];
        
        // store to beta the [gij]
        beta(0, 0) = inner_prod(g1, g1);
        beta(1, 1) = inner_prod(g2, g2);
        beta(2, 2) = inner_prod(g3, g3);
        
        beta(0, 1) = inner_prod(g1, g2);
        beta(1, 0) = beta(0, 1);
        
        beta(0, 2) = inner_prod(g1, g3);
        beta(2, 0) = beta(0, 2);
        
        beta(1, 2) = inner_prod(g2, g3);
        beta(2, 1) = beta(1, 2);

        rKinematicVariables.n_contravariant[0] = beta(0, 0) * n_covariant[0] + beta(0, 1) * n_covariant[1] + beta(0, 2) * n_covariant[2];
        rKinematicVariables.n_contravariant[1] = beta(1, 0) * n_covariant[0] + beta(1, 1) * n_covariant[1] + beta(1, 2) * n_covariant[2];
        rKinematicVariables.n_contravariant[2] = beta(2, 0) * n_covariant[0] + beta(2, 1) * n_covariant[1] + beta(2, 2) * n_covariant[2];

        // Check Numbers
        if (Id() == 550) {
            if (rConfiguration == ConfigurationType::Reference) {
                std::cout << "CalculateKinematicsSolid :: Solid :: Reference Congiguration" << std::endl;
                std::cout << "Contravatiants: n1:   " << rKinematicVariables.n_contravariant[0] << " | n2:  " << rKinematicVariables.n_contravariant[1] << " | n3:  " << rKinematicVariables.n_contravariant[2] << std::endl;
            }
            else if (rConfiguration == ConfigurationType::Current)
            {
                std::cout << "CalculateKinematicsSolid :: Current Congiguration" << std::endl;
            }
            std::cout << "G1 = ( " << g1[0] << " | " << g1[1] << " | " << g1[2] << " )" << std::endl;
            std::cout << "G2 = ( " << g2[0] << " | " << g2[1] << " | " << g2[2] << " )" << std::endl;
            std::cout << "G3 = ( " << g3[0] << " | " << g3[1] << " | " << g3[2] << " )" << std::endl;
        }

        KRATOS_CATCH("")
    }

    void CouplingSolidShellNitscheCondition::CalculateTransposeInverseMatrix3x3(Matrix& TransposedInverted, const Matrix Input) {
        KRATOS_TRY; 

        double determinant = +  Input(0, 0) * (Input(1, 1) * Input(2, 2) - Input(2, 1) * Input(1, 2))
                             -  Input(0, 1) * (Input(1, 0) * Input(2, 2) - Input(1, 2) * Input(2, 0))
                             +  Input(0, 2) * (Input(1, 0) * Input(2, 1) - Input(1, 1) * Input(2, 0));
        
        double invdet = 1 / determinant;

        KRATOS_ERROR_IF(std::abs(determinant) < 1e-12) << "Inverse of matrix not possible. Dterminant is singular" << std::endl;

        TransposedInverted(0, 0) =  (Input(1, 1) * Input(2, 2) - Input(2, 1) * Input(1, 2) )* invdet;
        TransposedInverted(1, 0) = -(Input(0, 1) * Input(2, 2) - Input(0, 2) * Input(2, 1) )* invdet;
        TransposedInverted(2, 0) =  (Input(0, 1) * Input(1, 2) - Input(0, 2) * Input(1, 1) )* invdet;
        TransposedInverted(0, 1) = -(Input(1, 0) * Input(2, 2) - Input(1, 2) * Input(2, 0) )* invdet;
        TransposedInverted(1, 1) =  (Input(0, 0) * Input(2, 2) - Input(0, 2) * Input(2, 0) )* invdet;
        TransposedInverted(2, 1) = -(Input(0, 0) * Input(1, 2) - Input(1, 0) * Input(0, 2) )* invdet;
        TransposedInverted(0, 2) =  (Input(1, 0) * Input(2, 1) - Input(2, 0) * Input(1, 1) )* invdet;
        TransposedInverted(1, 2) = -(Input(0, 0) * Input(2, 1) - Input(2, 0) * Input(0, 1) )* invdet;
        TransposedInverted(2, 2) =  (Input(0, 0) * Input(1, 1) - Input(1, 0) * Input(0, 1) )* invdet;

        KRATOS_CATCH("")
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void CouplingSolidShellNitscheCondition::CalculateTransformationShell(
        const KinematicVariablesShell& rKinematicVariables,
        Matrix& rT, Matrix& rT_hat, array_1d<array_1d<double, 3>,2>& rReferenceContraVariantBase
    )
    {
        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rKinematicVariables.a_ab_covariant[0] * rKinematicVariables.a_ab_covariant[1]
                - rKinematicVariables.a_ab_covariant[2] * rKinematicVariables.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rKinematicVariables.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rKinematicVariables.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rKinematicVariables.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_1 = rKinematicVariables.a1*a_ab_contravariant[0] + rKinematicVariables.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = rKinematicVariables.a1*a_ab_contravariant[2] + rKinematicVariables.a2*a_ab_contravariant[1];

        rReferenceContraVariantBase[0] = a_contravariant_1;
        rReferenceContraVariantBase[1] = a_contravariant_2;

        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        //Transformation matrix T from contravariant to local cartesian basis
        // e * a_contravariant
        Matrix G = ZeroMatrix(2, 2);
        G(0, 0) = inner_prod(e1, a_contravariant_1);
        G(0, 1) = inner_prod(e1, a_contravariant_2);
        G(1, 0) = inner_prod(e2, a_contravariant_1);
        G(1, 1) = inner_prod(e2, a_contravariant_2);

        //Transformation matrix T
        if (rT.size1() != 3 && rT.size2() != 3)
            rT.resize(3, 3);
        noalias(rT) = ZeroMatrix(3, 3);

        rT(0, 0) = pow(G(0, 0), 2);
        rT(0, 1) = pow(G(0, 1), 2);
        rT(0, 2) = 2 * G(0, 0) * G(0, 1);

        rT(1, 0) = pow(G(1, 0), 2);
        rT(1, 1) = pow(G(1, 1), 2);
        rT(1, 2) = 2 * G(1, 0) * G(1, 1);

        rT(2, 0) = 2 * G(0, 0) * G(1, 0);
        rT(2, 1) = 2 * G(0, 1) * G(1, 1);
        rT(2, 2) = 2 * (G(0, 0) * G(1, 1) + G(0, 1) * G(1, 0));

        //Transformation matrix T from local cartesian basis to covariant basis
        if (rT_hat.size1() != 3 && rT_hat.size2() != 3)
            rT_hat.resize(3, 3);
        noalias(rT_hat) = ZeroMatrix(3, 3);

        rT_hat(0, 0) = pow(G(0, 0), 2);
        rT_hat(0, 1) = pow(G(1, 0), 2);
        rT_hat(0, 2) = 2*G(0, 0) * G(1, 0);

        rT_hat(1, 0) = pow(G(0, 1), 2);
        rT_hat(1, 1) = pow(G(1, 1), 2);
        rT_hat(1, 2) = 2*G(0, 1) * G(1, 1);

        rT_hat(2, 0) = G(0, 0) * G(0, 1);
        rT_hat(2, 1) = G(1, 0) * G(1, 1);
        rT_hat(2, 2) = (G(0, 0) * G(1, 1) + G(1, 0) * G(0, 1));
    }

    void CouplingSolidShellNitscheCondition::CalculateConstitutiveVariablesShell(
        IndexType IntegrationPointIndex,
        const double& theta3,
        KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesShell,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector_slave[IntegrationPointIndex]);
        array_1d<double, 3> curvature_vector = m_B_ab_covariant_vector_slave[IntegrationPointIndex] - rActualKinematic.b_ab_covariant  ;
        
        // Kindl 3.33
        array_1d<double, 3> E_ab;
        E_ab[0] = strain_vector[0] + theta3 * curvature_vector[0];
        E_ab[1] = strain_vector[1] + theta3 * curvature_vector[1];
        E_ab[2] = strain_vector[2] + theta3 * curvature_vector[2];

        noalias(rThisConstitutiveVariablesShell.StrainVector) = prod(m_T_vector_slave[IntegrationPointIndex], E_ab);

        // Constitive Matrices DShell
        rValues.SetStrainVector(rThisConstitutiveVariablesShell.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesShell.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesShell.ConstitutiveMatrix); //this is an ouput parameter

        ConstitutiveLaw::Pointer constitutive_law_slave = GetProperties().GetSubProperties().back()[CONSTITUTIVE_LAW];
        constitutive_law_slave->InitializeMaterial(GetProperties(), GetGeometry().GetGeometryPart(1), row(GetGeometry().GetGeometryPart(1).ShapeFunctionsValues(), IntegrationPointIndex));

        constitutive_law_slave->CalculateMaterialResponse(rValues, ThisStressMeasure);

        //Local Cartesian PK2 Stresses
        noalias(rThisConstitutiveVariablesShell.StressVector) = prod(
            trans(rThisConstitutiveVariablesShell.ConstitutiveMatrix), rThisConstitutiveVariablesShell.StrainVector);
    }

    void CouplingSolidShellNitscheCondition::CalculateConstitutiveVariablesSolid(
        IndexType IntegrationPointIndex,
        KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {

        rThisConstitutiveVariablesSolid.StrainVector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector_master[IntegrationPointIndex]);
        // When using Voigt Notation for solids the strains should be ε = [ε_11, ε_22, ε_33, 2*ε_12, 2*ε13, 2*ε23] 
        rThisConstitutiveVariablesSolid.StrainVector[3] *= 2 ;
        rThisConstitutiveVariablesSolid.StrainVector[4] *= 2 ;
        rThisConstitutiveVariablesSolid.StrainVector[5] *= 2 ;

        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        // Constitive Matrices Solid
        rValues.SetStrainVector(rThisConstitutiveVariablesSolid.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesSolid.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesSolid.ConstitutiveMatrix); //this is an ouput parameter

        ConstitutiveLaw::Pointer constitutive_law_master = GetProperties().GetSubProperties().front()[CONSTITUTIVE_LAW];
        constitutive_law_master->InitializeMaterial(GetProperties(), GetGeometry().GetGeometryPart(0), row(GetGeometry().GetGeometryPart(0).ShapeFunctionsValues(), IntegrationPointIndex));

        constitutive_law_master->CalculateMaterialResponse(rValues, ThisStressMeasure);

        noalias(rThisConstitutiveVariablesSolid.StressVector) = prod(
            trans(rThisConstitutiveVariablesSolid.ConstitutiveMatrix), rThisConstitutiveVariablesSolid.StrainVector);

        if (Id() == 550) {
            auto & C = rThisConstitutiveVariablesSolid.ConstitutiveMatrix;
            std::cout << "Solid :: ConstitutiveMatrix size :" << C.size1() << " | " << C.size2() << std::endl;
            for (int ii = 0; ii < C.size1(); ii++) {
                for (int jj = 0; jj < C.size2(); jj++) {
                    std::cout << C(ii, jj) << " | ";
                    
                }
                 std::cout << std::endl;
                
            }
             auto & E = rThisConstitutiveVariablesSolid.StrainVector;
            std::cout << "Solid :: Strains size :" << E.size() << std::endl;
            for (int ii = 0; ii < E.size(); ii++) {
                std::cout << E(ii) << " | ";
                
            }
            
        }
    }

    void CouplingSolidShellNitscheCondition::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe) const
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (IndexType k = 0; k < number_of_nodes; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian(0, 0) += rDDN_DDe(k, 0) * coords[0];
            Hessian(0, 1) += rDDN_DDe(k, 2) * coords[0];
            Hessian(0, 2) += rDDN_DDe(k, 1) * coords[0];

            Hessian(1, 0) += rDDN_DDe(k, 0) * coords[1];
            Hessian(1, 1) += rDDN_DDe(k, 2) * coords[1];
            Hessian(1, 2) += rDDN_DDe(k, 1) * coords[1];

            Hessian(2, 0) += rDDN_DDe(k, 0) * coords[2];
            Hessian(2, 1) += rDDN_DDe(k, 2) * coords[2];
            Hessian(2, 2) += rDDN_DDe(k, 1) * coords[2];
        }
    }



    void CouplingSolidShellNitscheCondition::DeterminantOfJacobianInitial(
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

     void CouplingSolidShellNitscheCondition::CalculateTractionShell(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        array_1d<double, 2> n_contravariant_vector;
        
        stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        
        // Compute the stress components
        Matrix Palphabeta = ZeroMatrix(2, 2);
        Palphabeta(0,0) = stress_vector_covariant[0];
        Palphabeta(1,1) = stress_vector_covariant[1];
        Palphabeta(0,1) = stress_vector_covariant[2];
        Palphabeta(1,0) = Palphabeta(0,1);
        
        // Compute the traction vectors
        rTraction[0] = rActualKinematic.a1[0]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[0]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);

        rTraction[1] = rActualKinematic.a1[1]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[1]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);

        rTraction[2] = rActualKinematic.a1[2]*(Palphabeta(0,0)*n_contravariant_vector[0]+Palphabeta(0,1)*n_contravariant_vector[1]) 
                     + rActualKinematic.a2[2]*(Palphabeta(1,0)*n_contravariant_vector[0]+Palphabeta(1,1)*n_contravariant_vector[1]);
    }

    void CouplingSolidShellNitscheCondition::CalculateTractionSolid(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid)
    {
        array_1d<double, 3> n_contravariant_vector;
        n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        double n1 = n_contravariant_vector[0];
        double n2 = n_contravariant_vector[1];
        double n3 = n_contravariant_vector[2];

        // StressVector : [11,22,33,12,13,22]
        // n_contravariant_vector : [1,2,3]
        // traction = t^i * g_i
        // t^i = stress^ij n_j 
     
        // t^1 = stress^11 * n_1 + stress^12 * n_2 + stress^13 * n_3
        double t1_covariant = rThisConstitutiveVariablesSolid.StressVector[0] * n1 + rThisConstitutiveVariablesSolid.StressVector[5] * n2  + rThisConstitutiveVariablesSolid.StressVector[4] * n3;

        double t2_covariant = rThisConstitutiveVariablesSolid.StressVector[5] * n1 + rThisConstitutiveVariablesSolid.StressVector[1] * n2  + rThisConstitutiveVariablesSolid.StressVector[3] * n3;

        double t3_covariant = rThisConstitutiveVariablesSolid.StressVector[4] * n1 + rThisConstitutiveVariablesSolid.StressVector[3] * n2 + rThisConstitutiveVariablesSolid.StressVector[2] * n3;

        // Compute the traction vectors
        rTraction = t1_covariant * rActualKinematic.a1 + t2_covariant * rActualKinematic.a2 + t3_covariant * rActualKinematic.a3;
    }

    void CouplingSolidShellNitscheCondition::CalculateFirstVariationStressCovariantShell(
        IndexType IntegrationPointIndex,
        const double& theta3,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid)
    {
        KRATOS_TRY; 

        const auto& r_geometry = GetGeometry().GetGeometryPart(1);
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of the Green-Lagrange strains
        Matrix dE_cartesian = ZeroMatrix(3, mat_size);

        Matrix& T_patch = m_T_vector_slave[IntegrationPointIndex];
        
        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            const Matrix& rDDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());

            // Calculate second derivates of the base vectors wrt d2_theta1, d2_theta2^2, d2_theta1,
            array_1d<double, 3> a1_1 = column(rActualKinematic.Hessian, 0);
            array_1d<double, 3> a2_2 = column(rActualKinematic.Hessian, 1);
            array_1d<double, 3> a1_2 = column(rActualKinematic.Hessian, 2);
            
           
            // Calculate variations of the base vectors wrt to r-th dof
            array_1d<double, 3> a1_r(3, 0), a2_r(3, 0), a1_1_r(3, 0), a2_2_r(3, 0), a1_2_r(3, 0);
            a1_r[dirr] = r_DN_De(kr, 0);
            a2_r[dirr] = r_DN_De(kr, 1);

            // Calculate variations of the first derivate of the base vectors wrt to r-th dof
            a1_1_r[dirr]  = rDDN_DDe(kr, 0); 
            a2_2_r[dirr]  = rDDN_DDe(kr, 2);
            a1_2_r[dirr] =  rDDN_DDe(kr, 1);

            // Kindl 5.24
            array_1d<double, 3> a1_r__cross_a2;
            MathUtils<double>::CrossProduct(a1_r__cross_a2, a1_r, rActualKinematic.a2);

            array_1d<double, 3> a1__cross_a2_r;
            MathUtils<double>::CrossProduct(a1__cross_a2_r, rActualKinematic.a1, a2_r);

            array_1d<double, 3> a3tilde_r = a1_r__cross_a2 + a1__cross_a2_r; 

            // Kindl 5.25
            double a3dash_r = inner_prod(rActualKinematic.a3_tilde, a3tilde_r) / rActualKinematic.dA;

            // Kindl 5.26
            array_1d<double, 3>  a3_r = (a3tilde_r * rActualKinematic.dA) - (rActualKinematic.a3_tilde* a3dash_r) / pow(rActualKinematic.dA,2);

            // Kindl 5.29 & 5.28
            double k11_r = - ( inner_prod(a1_1_r, rActualKinematic.a3) + inner_prod(a1_1, a3_r) );
            double k22_r = - ( inner_prod(a2_2_r, rActualKinematic.a3) + inner_prod(a2_2, a3_r) );
            double k12_r = - ( inner_prod(a1_2_r, rActualKinematic.a3) + inner_prod(a1_2, a3_r) );

            array_1d<double, 3> dE_curvilinear;
            // strain Kindl 5.15
            double e11_r = inner_prod(a1_r, rActualKinematic.a1);
            double e22_r = inner_prod(a2_r, rActualKinematic.a2);
            double e12_r = 0.5 * ( inner_prod(a1_r, rActualKinematic.a2) + inner_prod(rActualKinematic.a1, a2_r));
            
            dE_curvilinear[0] = r_DN_De(kr, 0)*rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(kr, 1)*rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(kr, 1));
           
           
            // Kindl 3.36 for first variation
            dE_curvilinear[0] = e11_r + theta3 * k11_r;
            dE_curvilinear[1] = e22_r + theta3 * k22_r;
            dE_curvilinear[2] = e12_r + theta3 * k12_r;

            // Transform to local Cartesian bases
            dE_cartesian(0, r) = T_patch(0, 0)*dE_curvilinear[0] + T_patch(0, 1)*dE_curvilinear[1] + T_patch(0, 2)*dE_curvilinear[2];
            dE_cartesian(1, r) = T_patch(1, 0)*dE_curvilinear[0] + T_patch(1, 1)*dE_curvilinear[1] + T_patch(1, 2)*dE_curvilinear[2];
            dE_cartesian(2, r) = T_patch(2, 0)*dE_curvilinear[0] + T_patch(2, 1)*dE_curvilinear[1] + T_patch(2, 2)*dE_curvilinear[2];

        }

        //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
        Matrix first_variations_stress_cartesian = ZeroMatrix(3, mat_size);
        first_variations_stress_cartesian = prod(rThisConstitutiveVariablesSolid.ConstitutiveMatrix,dE_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
        rFirstVariationStressCovariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], first_variations_stress_cartesian);

        KRATOS_CATCH("")
    }

    void CouplingSolidShellNitscheCondition::CalculateFirstVariationStressCovariantSolid(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid)
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(0);

        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of the Green-Lagrange strains
        Matrix dE = ZeroMatrix(6, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            // strain
            dE(0,r) = r_DN_De(kr, 0) * rActualKinematic.a1(dirr);
            dE(1,r) = r_DN_De(kr, 1) * rActualKinematic.a2(dirr);
            dE(2,r) = r_DN_De(kr, 2) * rActualKinematic.a3(dirr);
            dE(3,r) = (r_DN_De(kr, 1) * rActualKinematic.a3(dirr) + rActualKinematic.a2(dirr) * r_DN_De(kr, 2)); //23
            dE(4,r) = (r_DN_De(kr, 0) * rActualKinematic.a3(dirr) + rActualKinematic.a1(dirr) * r_DN_De(kr, 2)); //13
            dE(5,r) = (r_DN_De(kr, 1) * rActualKinematic.a1(dirr) + rActualKinematic.a2(dirr) * r_DN_De(kr, 0)); //12
        }

        //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
  
        rFirstVariationStressCovariant = prod(rThisConstitutiveVariablesSolid.ConstitutiveMatrix, dE);
    }
    
    void CouplingSolidShellNitscheCondition::CalculateFirstVariationTractionShell(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesShell)
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(1);
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //get the normal vector
        array_1d<double, 2> n_contravariant_vector; 
        n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];

        //Compute the first variation of the traction vectors:
        //1. normal vector * derivative stress covariant
        
        // normal vector * covariant base vector
        Matrix n_a = ZeroMatrix(3, 3); 

        for (IndexType r = 0; r < 3; r++)
        {
            n_a (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
            n_a (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
            n_a (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
        }

        rFirstVariationTraction = prod(n_a, rFirstVariationStressCovariant);

        //2. derivative normal vector * stress covariant
        
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        stress_vector_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], rThisConstitutiveVariablesShell.StressVector);

        Matrix r_DN_Dxi = ZeroMatrix(3, mat_size);
        Matrix r_DN_Deta = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            r_DN_Dxi(0, 3 * r) = r_DN_De(r, 0);
            r_DN_Dxi(1, 3 * r + 1) = r_DN_De(r, 0);
            r_DN_Dxi(2, 3 * r + 2) = r_DN_De(r, 0);

            r_DN_Deta(0, 3 * r) = r_DN_De(r, 1);
            r_DN_Deta(1, 3 * r + 1) = r_DN_De(r, 1);
            r_DN_Deta(2, 3 * r + 2) = r_DN_De(r, 1);
        }

        rFirstVariationTraction += r_DN_Dxi*(n_contravariant_vector[0]*stress_vector_covariant[0] + n_contravariant_vector[1]*stress_vector_covariant[2])+
                                   r_DN_Deta*(n_contravariant_vector[1]*stress_vector_covariant[1] + n_contravariant_vector[0]*stress_vector_covariant[2]);
    }

    void CouplingSolidShellNitscheCondition::CalculateFirstVariationTractionSolid(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid)
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(0);

        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of contravariant coeeficient * base vectors:

        // normal vector * covariant base vector
        Matrix n_a = ZeroMatrix(3, 6);

        //get the normal vector
        array_1d<double, 3> n_contravariant_vector;
        n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        double n1, n2, n3;
        n1 = n_contravariant_vector[0];
        n2 = n_contravariant_vector[1];
        n3 = n_contravariant_vector[2];

        for (IndexType r = 0; r < 3; r++)
        {
            n_a(r, 0) = rActualKinematic.a1[r] * n1;
            n_a(r, 1) = rActualKinematic.a2[r] * n2;
            n_a(r, 2) = rActualKinematic.a3[r] * n3;
            n_a(r, 3) = rActualKinematic.a2[r] * n3 + rActualKinematic.a3[r] * n2; //23 
            n_a(r, 4) = rActualKinematic.a1[r] * n3 + rActualKinematic.a3[r] * n1; //13
            n_a(r, 5) = rActualKinematic.a1[r] * n2 + rActualKinematic.a2[r] * n1; //12
        }

        rFirstVariationTraction = prod(n_a, rFirstVariationStressCovariant);

        //2. derivative normal vector * stress covariant

        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
 
        Matrix r_DN_theta1 = ZeroMatrix(3, mat_size);
        Matrix r_DN_theta2 = ZeroMatrix(3, mat_size);
        Matrix r_DN_theta3 = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            r_DN_theta1(0, 3 * r) = r_DN_De(r, 0);
            r_DN_theta1(1, 3 * r + 1) = r_DN_De(r, 0);
            r_DN_theta1(2, 3 * r + 2) = r_DN_De(r, 0);

            r_DN_theta2(0, 3 * r) = r_DN_De(r, 1);
            r_DN_theta2(1, 3 * r + 1) = r_DN_De(r, 1);
            r_DN_theta2(2, 3 * r + 2) = r_DN_De(r, 1);

            r_DN_theta3(0, 3 * r) = r_DN_De(r, 2);
            r_DN_theta3(1, 3 * r + 1) = r_DN_De(r, 2);
            r_DN_theta3(2, 3 * r + 2) = r_DN_De(r, 2);
        }

        double traction_covariant_1 = n1 * rThisConstitutiveVariablesSolid.StressVector[0] + n2 * rThisConstitutiveVariablesSolid.StressVector[5] 
            + n3 * rThisConstitutiveVariablesSolid.StressVector[4];
        double traction_covariant_2 = n1 * rThisConstitutiveVariablesSolid.StressVector[5] + n2 * rThisConstitutiveVariablesSolid.StressVector[1]
            + n3 * rThisConstitutiveVariablesSolid.StressVector[3];
        double traction_covariant_3 = n1 * rThisConstitutiveVariablesSolid.StressVector[4] + n2 * rThisConstitutiveVariablesSolid.StressVector[3]
            + n3 * rThisConstitutiveVariablesSolid.StressVector[2];

        rFirstVariationTraction += r_DN_theta1 * traction_covariant_1 + r_DN_theta2 * traction_covariant_2 + r_DN_theta3 * traction_covariant_3;
     }

    void CouplingSolidShellNitscheCondition::CalculateSecondVariationTractionSolid(
        IndexType IntegrationPointIndex,
        Matrix& rProductSecondVariationTraction_Displacement,
        const Matrix& rFirstVariationStressCovariant,
        const array_1d<double, 3>& rDisplacementMaster,
        const array_1d<double, 3>& rDisplacementSlave,
        const KinematicVariablesSolid& rActualKinematic,
        const ConstitutiveVariables& rThisConstitutiveVariablesSolid) {

        const auto& r_geometry = GetGeometry().GetGeometryPart(0);

        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3; //dofs

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        Matrix r_DN_theta1 = ZeroMatrix(3, mat_size);
        Matrix r_DN_theta2 = ZeroMatrix(3, mat_size);
        Matrix r_DN_theta3 = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < number_of_control_points; r++)
        {
            r_DN_theta1(0, 3 * r) = r_DN_De(r, 0);
            r_DN_theta1(1, 3 * r + 1) = r_DN_De(r, 0);
            r_DN_theta1(2, 3 * r + 2) = r_DN_De(r, 0);

            r_DN_theta2(0, 3 * r) = r_DN_De(r, 1);
            r_DN_theta2(1, 3 * r + 1) = r_DN_De(r, 1);
            r_DN_theta2(2, 3 * r + 2) = r_DN_De(r, 1);

            r_DN_theta3(0, 3 * r) = r_DN_De(r, 2);
            r_DN_theta3(1, 3 * r + 1) = r_DN_De(r, 2);
            r_DN_theta3(2, 3 * r + 2) = r_DN_De(r, 2);
        }

        //get the normal vector
        array_1d<double, 3> n_contravariant_vector;
        n_contravariant_vector = m_n_contravariant_vector_master[IntegrationPointIndex];
        double n1, n2, n3;
        n1 = n_contravariant_vector[0];
        n2 = n_contravariant_vector[1];
        n3 = n_contravariant_vector[2];

        array_1d<double, 6> ddE, ddS;
        array_1d<double, 3> second_variation_traction_covariant, first_variation_traction_covariant, rSecondVariationTraction;

        for (size_t dk = 0; dk < mat_size; ++dk) {

            // define variations of base vectors wrt dk
            array_1d<double, 3> r_DN_theta1_dk;
            r_DN_theta1_dk(0) = r_DN_theta1(0, dk);
            r_DN_theta1_dk(1) = r_DN_theta1(1, dk);
            r_DN_theta1_dk(2) = r_DN_theta1(2, dk);

            array_1d<double, 3> r_DN_theta2_dk;
            r_DN_theta2_dk(0) = r_DN_theta2(0, dk);
            r_DN_theta2_dk(1) = r_DN_theta2(1, dk);
            r_DN_theta2_dk(2) = r_DN_theta2(2, dk);

            array_1d<double, 3> r_DN_theta3_dk;
            r_DN_theta3_dk(0) = r_DN_theta3(0, dk);
            r_DN_theta3_dk(1) = r_DN_theta3(1, dk);
            r_DN_theta3_dk(2) = r_DN_theta3(2, dk);

            for (size_t dr = 0; dr < mat_size; ++dr) {

                // define variations of base vectors wrt dr
                array_1d<double, 3> r_DN_theta1_dr;
                r_DN_theta1_dr(0) = r_DN_theta1(0, dr);
                r_DN_theta1_dr(1) = r_DN_theta1(1, dr);
                r_DN_theta1_dr(2) = r_DN_theta1(2, dr);

                array_1d<double, 3> r_DN_theta2_dr;
                r_DN_theta2_dr(0) = r_DN_theta2(0, dr);
                r_DN_theta2_dr(1) = r_DN_theta2(1, dr);
                r_DN_theta2_dr(2) = r_DN_theta2(2, dr);

                array_1d<double, 3> r_DN_theta3_dr;
                r_DN_theta3_dr(0) = r_DN_theta3(0, dr);
                r_DN_theta3_dr(1) = r_DN_theta3(1, dr);
                r_DN_theta3_dr(2) = r_DN_theta3(2, dr);

                ddE(0) = inner_prod(r_DN_theta1_dk, r_DN_theta1_dr);
                ddE(1) = inner_prod(r_DN_theta2_dk, r_DN_theta2_dr);
                ddE(2) = inner_prod(r_DN_theta3_dk, r_DN_theta3_dr);
                       
                ddE(3) = inner_prod(r_DN_theta2_dk, r_DN_theta3_dr) + inner_prod(r_DN_theta3_dk, r_DN_theta2_dr); //23
                ddE(4) = inner_prod(r_DN_theta1_dk, r_DN_theta3_dr) + inner_prod(r_DN_theta3_dk, r_DN_theta1_dr); //13
                ddE(5) = inner_prod(r_DN_theta1_dk, r_DN_theta2_dr) + inner_prod(r_DN_theta2_dk, r_DN_theta1_dr); //12 

                // Second variation of stress covariants
                ddS = prod(rThisConstitutiveVariablesSolid.ConstitutiveMatrix, ddE);

                // Second Traction Variation coefficients for covariant base vectors 
                second_variation_traction_covariant[0] = n1 * ddS[0] + n2 * ddS[5] + n3 * ddS[4]; // 
                second_variation_traction_covariant[1] = n1 * ddS[5] + n2 * ddS[1] + n3 * ddS[3]; // 
                second_variation_traction_covariant[2] = n1 * ddS[4] + n2 * ddS[3] + n3 * ddS[2]; // 

                // term1 : secondvariation of traction covariants wrt dk * base vector
                // Second Traction Variation  expressed for CARTESIAN basis vectors
                rSecondVariationTraction[0] = second_variation_traction_covariant[0] * rActualKinematic.a1[0] + second_variation_traction_covariant[1] * rActualKinematic.a2[0]
                    + second_variation_traction_covariant[2] * rActualKinematic.a3[0]; //

                rSecondVariationTraction[1] = second_variation_traction_covariant[0] * rActualKinematic.a1[1] + second_variation_traction_covariant[1] * rActualKinematic.a2[1]
                   + second_variation_traction_covariant[2] * rActualKinematic.a3[1]; //

                rSecondVariationTraction[2] = second_variation_traction_covariant[0] * rActualKinematic.a1[2] + second_variation_traction_covariant[1] * rActualKinematic.a2[2]
                    + second_variation_traction_covariant[2] * rActualKinematic.a3[2]; //

                // term2 : first variation of traction covariants wrt dk * first variation of base vector wrt dr
                first_variation_traction_covariant[0] = n1 * rFirstVariationStressCovariant(0,dk) + n2 * rFirstVariationStressCovariant(5,dk) + n3 * rFirstVariationStressCovariant(4,dk); // 
                first_variation_traction_covariant[1] = n1 * rFirstVariationStressCovariant(5,dk) + n2 * rFirstVariationStressCovariant(1,dk) + n3 * rFirstVariationStressCovariant(3,dk); // 
                first_variation_traction_covariant[2] = n1 * rFirstVariationStressCovariant(4,dk) + n2 * rFirstVariationStressCovariant(3,dk) + n3 * rFirstVariationStressCovariant(2,dk); // 

                rSecondVariationTraction[0] += first_variation_traction_covariant[0] * r_DN_theta1_dr[0] + first_variation_traction_covariant[1] * r_DN_theta2_dr[0] 
                    + first_variation_traction_covariant[2] * r_DN_theta3_dr[0];

                rSecondVariationTraction[1] += first_variation_traction_covariant[0] * r_DN_theta1_dr[1] + first_variation_traction_covariant[1] * r_DN_theta2_dr[1]
                    + first_variation_traction_covariant[2] * r_DN_theta3_dr[1];

                rSecondVariationTraction[2] += first_variation_traction_covariant[0] * r_DN_theta1_dr[2] + first_variation_traction_covariant[1] * r_DN_theta2_dr[2]
                    + first_variation_traction_covariant[2] * r_DN_theta3_dr[2];

                // term3 : first variation of traction covariants wrt dr * first variation of base vector wrt dk
                first_variation_traction_covariant[0] = n1 * rFirstVariationStressCovariant(0, dr) + n2 * rFirstVariationStressCovariant(5, dr) + n3 * rFirstVariationStressCovariant(4, dr); // 
                first_variation_traction_covariant[1] = n1 * rFirstVariationStressCovariant(5, dr) + n2 * rFirstVariationStressCovariant(1, dr) + n3 * rFirstVariationStressCovariant(3, dr); // 
                first_variation_traction_covariant[2] = n1 * rFirstVariationStressCovariant(4, dr) + n2 * rFirstVariationStressCovariant(3, dr) + n3 * rFirstVariationStressCovariant(2, dr); // 
            
                rSecondVariationTraction[0] += first_variation_traction_covariant[0] * r_DN_theta1_dk[0] + first_variation_traction_covariant[1] * r_DN_theta2_dk[0]
                    + first_variation_traction_covariant[2] * r_DN_theta3_dk[0];

                rSecondVariationTraction[1] += first_variation_traction_covariant[0] * r_DN_theta1_dk[1] + first_variation_traction_covariant[1] * r_DN_theta2_dk[1]
                    + first_variation_traction_covariant[2] * r_DN_theta3_dk[1];

                rSecondVariationTraction[2] += first_variation_traction_covariant[0] * r_DN_theta1_dk[2] + first_variation_traction_covariant[1] * r_DN_theta2_dk[2]
                    + first_variation_traction_covariant[2] * r_DN_theta3_dk[2];

                // Inner product of second variation of traction * displacements
                rProductSecondVariationTraction_Displacement(dk, dr)  = ( rSecondVariationTraction[0] * (rDisplacementMaster(0) - rDisplacementSlave(0)) ) + 
                ( rSecondVariationTraction[1] * (rDisplacementMaster(1) - rDisplacementSlave(1)) ) +
                ( rSecondVariationTraction[2] * (rDisplacementMaster(2) - rDisplacementSlave(2))  );
            }
        }
    }

    void CouplingSolidShellNitscheCondition::CalculateSecondVariationTractionShell(
        IndexType IntegrationPointIndex,
        Matrix& rProductSecondVariationTraction_Displacement,
        const double& rtheta3,
        const Matrix& rFirstVariationStressCovariant,
        const array_1d<double, 3>& rDisplacementMaster,
        const array_1d<double, 3>& rDisplacementSlave,
        const KinematicVariablesShell& rActualKinematic,
        const ConstitutiveVariables& rThisConstitutiveVariablesShell)
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(1);

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
 
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        // Calculate second derivates of the base vectors wrt d2_theta1, d2_theta2^2, d_theta1* d_theta2
        array_1d<double, 3> a1_1 = column(rActualKinematic.Hessian, 0);
        array_1d<double, 3> a2_2 = column(rActualKinematic.Hessian, 1);
        array_1d<double, 3> a1_2 = column(rActualKinematic.Hessian, 2);

        const Matrix& rDDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());
        
        Matrix& T_patch = m_T_vector_slave[IntegrationPointIndex];

        std::vector<array_1d<double, 3>> a1_r(mat_size);

        //get the normal vector
        array_1d<double, 2> n_contravariant_vector = m_n_contravariant_vector_slave[IntegrationPointIndex];
        double n1, n2;
        n1 = n_contravariant_vector[0];
        n2 = n_contravariant_vector[1];

        for (IndexType r = 0; r < mat_size; r++) { // row

            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            // Calculate variations of the base vectors wrt to r-th dof
            array_1d<double, 3> a1_r(3, 0), a2_r(3, 0), a1_1_r(3, 0), a2_2_r(3, 0), a1_2_r(3, 0);
            a1_r[dirr] = r_DN_De(kr, 0);
            a2_r[dirr] = r_DN_De(kr, 1);

            // Calculate variations of the first derivate of the base vectors wrt to r-th dof
            a1_1_r[dirr] = rDDN_DDe(kr, 0);
            a2_2_r[dirr] = rDDN_DDe(kr, 2);
            a1_2_r[dirr] = rDDN_DDe(kr, 1);

            // Kindl 5.24
            array_1d<double, 3> a1_r__cross_a2;
            MathUtils<double>::CrossProduct(a1_r__cross_a2, a1_r, rActualKinematic.a2);

            array_1d<double, 3> a1__cross_a2_r;
            MathUtils<double>::CrossProduct(a1__cross_a2_r, rActualKinematic.a1, a2_r);

            array_1d<double, 3> a3tilde_r = a1_r__cross_a2 + a1__cross_a2_r;

            // Kindl 5.25
            double a3dash_r = inner_prod(rActualKinematic.a3_tilde, a3tilde_r) / rActualKinematic.dA;

            // Kindl 5.26
            array_1d<double, 3>  a3_r = (a3tilde_r * rActualKinematic.dA) - (rActualKinematic.a3_tilde * a3dash_r) / pow(rActualKinematic.dA, 2);

            for (IndexType s = 0; s < mat_size; s++){ // column{

                // local node number ks and dof direction dirs
                IndexType ks = s / 3;
                IndexType dirs = s % 3;

                // Calculate variations of the base vectors wrt to s-th dof
                array_1d<double, 3> a1_s(3, 0), a2_s(3, 0), a1_1_s(3, 0), a2_2_s(3, 0), a1_2_s(3, 0);
                a1_s[dirs] = r_DN_De(ks, 0);
                a2_s[dirs] = r_DN_De(ks, 1);

                // Calculate variations of the first derivate of the base vectors wrt to s-th dof
                a1_1_s[dirs] = rDDN_DDe(ks, 0);
                a2_2_s[dirs] = rDDN_DDe(ks, 2);
                a1_2_s[dirs] = rDDN_DDe(ks, 1);

                // Kindl 5.24, 5.25, 5.26 for s-th dof
                array_1d<double, 3> a1_s__cross_a2;
                MathUtils<double>::CrossProduct(a1_s__cross_a2, a1_s, rActualKinematic.a2);

                array_1d<double, 3> a1__cross_a2_s;
                MathUtils<double>::CrossProduct(a1__cross_a2_s, rActualKinematic.a1, a2_s);

                array_1d<double, 3> a3tilde_s = a1_s__cross_a2 + a1__cross_a2_s;

                double a3dash_s = inner_prod(rActualKinematic.a3_tilde, a3tilde_s) / rActualKinematic.dA;

                array_1d<double, 3>  a3_s = (a3tilde_s * rActualKinematic.dA) - (rActualKinematic.a3_tilde * a3dash_s) / pow(rActualKinematic.dA, 2);

                // Kindl 5.30
                array_1d<double, 3> a1_r__cross_a2_s, a1_s__cross_a2_r;
                MathUtils<double>::CrossProduct(a1_r__cross_a2_s, a1_r, a2_s);
                MathUtils<double>::CrossProduct(a1_s__cross_a2_r, a1_s, a2_r);
                array_1d<double, 3> a3tilde_rs = a1_r__cross_a2_s + a1_s__cross_a2_r;

                // Kindl 5.31
                double a3dash_rs = (inner_prod(a3tilde_rs, rActualKinematic.a3_tilde) + inner_prod(a3tilde_r, a3tilde_s)) / rActualKinematic.dA;
                a3dash_rs -= inner_prod(a3tilde_r, rActualKinematic.a3_tilde) * inner_prod(a3tilde_s, rActualKinematic.a3_tilde) / pow(rActualKinematic.dA, 3);

                // Kindl 5.32
                array_1d<double, 3> a3_rs = a3tilde_rs / rActualKinematic.dA - a3tilde_r * a3dash_r / pow(rActualKinematic.dA, 2) - a3tilde_s * a3dash_r / pow(rActualKinematic.dA, 2)
                    - rActualKinematic.a3_tilde * a3dash_rs + rActualKinematic.a3_tilde * (2 * a3dash_r * a3dash_s) / pow(rActualKinematic.dA, 3);

                // Calculate second variations of strains
                
                // strain Kindl 5.19
                double e11_rs = inner_prod(a1_r, a1_s);
                double e22_rs = inner_prod(a2_r, a2_s);
                double e12_rs = 0.5 * (inner_prod(a1_r, a2_s) + inner_prod(a1_s,a2_r));

                // strain Kindl 5.34, 5.35
                double k11_rs = -(inner_prod(a1_1_r, a3_s) + inner_prod(a1_1_s, a3_r) + inner_prod(a1_1, a3_rs));
                double k22_rs = -(inner_prod(a2_2_r, a3_s) + inner_prod(a2_2_s, a3_r) + inner_prod(a2_2, a3_rs));
                double k12_rs = -(inner_prod(a1_2_r, a3_s) + inner_prod(a1_2_s, a3_r) + inner_prod(a1_2, a3_rs)); 

                // Kindl 3.36 for second variation
                array_1d<double, 3> dE_curvilinear;
                dE_curvilinear[0] = e11_rs + rtheta3 * k11_rs;
                dE_curvilinear[1] = e22_rs + rtheta3 * k22_rs;
                dE_curvilinear[2] = e12_rs + rtheta3 * k12_rs;

                // Transform to local Cartesian bases
                array_1d<double, 3> dE_cartesian;
                dE_cartesian[0] = T_patch(0, 0) * dE_curvilinear[0] + T_patch(0, 1) * dE_curvilinear[1] + T_patch(0, 2) * dE_curvilinear[2];
                dE_cartesian[1] = T_patch(1, 0) * dE_curvilinear[0] + T_patch(1, 1) * dE_curvilinear[1] + T_patch(1, 2) * dE_curvilinear[2];
                dE_cartesian[2] = T_patch(2, 0) * dE_curvilinear[0] + T_patch(2, 1) * dE_curvilinear[1] + T_patch(2, 2) * dE_curvilinear[2];

                //Compute the second variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
                array_1d<double, 3> second_variation_stress_cartesian = prod(rThisConstitutiveVariablesShell.ConstitutiveMatrix, dE_cartesian);

                //Transform the second variations of the 2nd Piola-Kichhoff stresses at the covariant bases
                array_1d<double, 3>  second_variation_stress_covariant = prod(m_T_hat_vector_slave[IntegrationPointIndex], second_variation_stress_cartesian);

                array_1d<double, 2> second_variation_traction_covariant;
                second_variation_traction_covariant[0] = second_variation_stress_covariant[0] * n1 + second_variation_stress_covariant[2] * n2;
                second_variation_traction_covariant[1] = second_variation_stress_covariant[2] * n1 + second_variation_stress_covariant[1] * n2;

                array_1d<double, 3> rSecondVariationTraction;

                // term1 : secondvariation of traction covariants wrt dr ds * base vector
                rSecondVariationTraction[0] = second_variation_traction_covariant[0] * rActualKinematic.a1[0] + second_variation_traction_covariant[1] * rActualKinematic.a2[0];
                rSecondVariationTraction[1] = second_variation_traction_covariant[0] * rActualKinematic.a1[1] + second_variation_traction_covariant[1] * rActualKinematic.a2[1];
                rSecondVariationTraction[2] = second_variation_traction_covariant[0] * rActualKinematic.a1[2] + second_variation_traction_covariant[1] * rActualKinematic.a2[2];

                // term2 : first variation of traction covariants wrt dr * first variation of base vector wrt ds
                array_1d<double, 2> first_variation_traction_covariant;
                first_variation_traction_covariant[0] = rFirstVariationStressCovariant(0, r) * n1 + rFirstVariationStressCovariant(2, r) * n2;
                first_variation_traction_covariant[1] = rFirstVariationStressCovariant(2, r) * n1 + rFirstVariationStressCovariant(1, r) * n2;

                rSecondVariationTraction[0] += first_variation_traction_covariant[0] * a1_s[0] + first_variation_traction_covariant[1] * a2_s[0];
                rSecondVariationTraction[1] += first_variation_traction_covariant[0] * a1_s[1] + first_variation_traction_covariant[1] * a2_s[1];
                rSecondVariationTraction[2] += first_variation_traction_covariant[0] * a1_s[2] + first_variation_traction_covariant[1] * a2_s[2];

                // term3 : first variation of traction covariants wrt ds * first variation of base vector wrt dr
                first_variation_traction_covariant[0] = rFirstVariationStressCovariant(0, s) * n1 + rFirstVariationStressCovariant(2, s) * n2;
                first_variation_traction_covariant[1] = rFirstVariationStressCovariant(2, s) * n1 + rFirstVariationStressCovariant(1, s) * n2;

                rSecondVariationTraction[0] += first_variation_traction_covariant[0] * a1_r[0] + first_variation_traction_covariant[1] * a2_r[0];
                rSecondVariationTraction[1] += first_variation_traction_covariant[0] * a1_r[1] + first_variation_traction_covariant[1] * a2_r[1];
                rSecondVariationTraction[2] += first_variation_traction_covariant[0] * a1_r[2] + first_variation_traction_covariant[1] * a2_r[2];

                rProductSecondVariationTraction_Displacement(r, s) = -  rSecondVariationTraction[0] * (rDisplacementMaster(0) - rDisplacementSlave(0)) 
                                                                     -  rSecondVariationTraction[1] * (rDisplacementMaster(1) - rDisplacementSlave(1)) 
                                                                     -  rSecondVariationTraction[2] * (rDisplacementMaster(2) - rDisplacementSlave(2)) ;
            }
        }

    }

    double CouplingSolidShellNitscheCondition::ComputeTheta3Shell(
        IndexType IntegrationPointIndex,
        const KinematicVariablesShell& rCurrentConfiguraitonKinematicVariablesShell) {
        
        const auto& Solid_geometry = GetGeometry().GetGeometryPart(0);
        const auto& Shell_geometry = GetGeometry().GetGeometryPart(1);

        // Get Global coordinates of integrations point
        GeometryType::CoordinatesArrayType  SolidIPGlobalCoordinates, ShellIPGlobalCoordinates;
        Solid_geometry.GlobalCoordinates(SolidIPGlobalCoordinates, IntegrationPointIndex);
        Shell_geometry.GlobalCoordinates(ShellIPGlobalCoordinates, IntegrationPointIndex);

        const auto& A3 = _A3[IntegrationPointIndex];

        double theta3 = (SolidIPGlobalCoordinates(0) - ShellIPGlobalCoordinates(0)) * A3(0) +
            (SolidIPGlobalCoordinates(1) - ShellIPGlobalCoordinates(1)) * A3(1) +
            (SolidIPGlobalCoordinates(2) - ShellIPGlobalCoordinates(2)) * A3(2);

        theta3 = theta3 / pow(norm_2(A3),2);

        // Check that SlavePoint + theta3 * A3 = MasterPoint
        array_1d<double, 3> distance;
        distance = ShellIPGlobalCoordinates + A3 * theta3 - SolidIPGlobalCoordinates;
        double thickness = GetProperties().GetSubProperties().back()[THICKNESS];
        
        if (abs(theta3) > thickness / 2) {
            std::cout << "Theta3 > thickness/2" << std::endl;
            KRATOS_ERROR;
        }
        if (norm_2(distance) > 1E-10) {
            std::cout << "Theta3 caclulation is wrong" << std::endl;
            KRATOS_ERROR;
        }

        return theta3;

    }

    void CouplingSolidShellNitscheCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_control_points_master = r_geometry_master.size();
        const SizeType number_of_control_points_slave = r_geometry_slave.size();
        const SizeType mat_size = (number_of_control_points_master + number_of_control_points_slave) * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points_master; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }

        for (IndexType i = 0; i < number_of_control_points_slave; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = 3 * (i + number_of_control_points_master) ;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void CouplingSolidShellNitscheCondition::OutOfPlaneDeformationFirstVariation(
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
            array_1d<double, 3> Phi_r_cross_A3(3, 0);
            Phi_r_cross_A3 = Calculate_Phi_r_cross_A3(N_theta1_r, N_theta2_r, A1, A2);

            OutOfPlaneDeformationWholeMatrix(0, r) = theta3 * Phi_r_cross_A3[0];
            OutOfPlaneDeformationWholeMatrix(1, r) = theta3 * Phi_r_cross_A3[1];
            OutOfPlaneDeformationWholeMatrix(2, r) = theta3 * Phi_r_cross_A3[2];

        }
    }

    array_1d<double, 3> CouplingSolidShellNitscheCondition::Calculate_Phi_r_cross_A3(
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
        double phi2 = -(1 / norm_A1_cross_A2) * inner_prod(N_theta1_r, A3);
        Phi = phi1 * A1 + phi2 * A2;

        MathUtils<double>::CrossProduct(Phi_r_cross_A3, Phi, A3);
        return Phi_r_cross_A3;
    }

    void CouplingSolidShellNitscheCondition::EquationIdVector(
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

    void CouplingSolidShellNitscheCondition::GetDofList(
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



