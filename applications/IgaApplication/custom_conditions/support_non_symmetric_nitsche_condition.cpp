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
#include "custom_conditions/support_non_symmetric_nitsche_condition.h"
#include "utilities/atomic_utilities.h"

// Project includes

namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void SupportNonSymmetricNitscheCondition::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration
    )
    {
        const auto& r_geometry = GetGeometry();

        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType number_of_nodes = r_geometry.size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);

        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);

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
        GetGeometry().Calculate(LOCAL_TANGENT, local_tangent);

        rKinematicVariables.t = local_tangent[0]*g1 + local_tangent[1]*g2;
        MathUtils<double>::CrossProduct(rKinematicVariables.n, rKinematicVariables.t/norm_2(rKinematicVariables.t), rKinematicVariables.a3);

        // transform the normal into the contavariant basis
        rKinematicVariables.n_contravariant[0] = rKinematicVariables.a1[0]*rKinematicVariables.n[0] + rKinematicVariables.a1[1]*rKinematicVariables.n[1] + rKinematicVariables.a1[2]*rKinematicVariables.n[2];
        rKinematicVariables.n_contravariant[1] = rKinematicVariables.a2[0]*rKinematicVariables.n[0] + rKinematicVariables.a2[1]*rKinematicVariables.n[1] + rKinematicVariables.a2[2]*rKinematicVariables.n[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void SupportNonSymmetricNitscheCondition::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
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

        //Transformation matrix T from local cartesian basis to the covariant one
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

    void SupportNonSymmetricNitscheCondition::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector[IntegrationPointIndex]);
        noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector[IntegrationPointIndex], strain_vector);

        // Constitive Matrices DMembrane
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        ConstitutiveLaw::Pointer constitutive_law = GetProperties().GetSubProperties().front()[CONSTITUTIVE_LAW];
        constitutive_law->InitializeMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(), IntegrationPointIndex));

        constitutive_law->CalculateMaterialResponse(rValues, ThisStressMeasure);
        rThisConstitutiveVariablesMembrane.ConstitutiveMatrix *= GetProperties().GetSubProperties().front()[THICKNESS];

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
    }

    //Prestress Transformation Matrix
    void SupportNonSymmetricNitscheCondition::CalculateTransformationPrestress(
        Matrix& rTransformationPrestress,
        const KinematicVariables& rActualKinematic
    )
    {
        //define base vector in reference plane
        array_1d<double, 3> t1;
        array_1d<double, 3> t2;
        
        if (Has(LOCAL_PRESTRESS_AXIS_1))
        {
            t1 = GetValue(LOCAL_PRESTRESS_AXIS_1);
            MathUtils<double>::CrossProduct(t2, rActualKinematic.a3, t1);
        }
        else if (Has(LOCAL_PRESTRESS_AXIS_1) && Has(LOCAL_PRESTRESS_AXIS_2))
        {
            t1 = GetValue(LOCAL_PRESTRESS_AXIS_1);
            t2 = GetValue(LOCAL_PRESTRESS_AXIS_2); 
        }

        array_1d<double, 3> t1_n = t1/norm_2(t1);
        array_1d<double, 3> t2_n = t2/norm_2(t2);

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rActualKinematic.a_ab_covariant[0] * rActualKinematic.a_ab_covariant[1]
                - rActualKinematic.a_ab_covariant[2] * rActualKinematic.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rActualKinematic.a_ab_covariant[2];

        array_1d<double, 3> a_contravariant_2 = rActualKinematic.a1*a_ab_contravariant[2] + rActualKinematic.a2*a_ab_contravariant[1]; 

        //local cartesian coordinates oriented along the 1st base vector in the ref. config.
        double l_a1 = norm_2(rActualKinematic.a1);
        array_1d<double, 3> e1 = rActualKinematic.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        //Transformation matrix from the projected basis T to the local cartesian basis
        double eG11 = inner_prod(e1,t1_n);
        double eG12 = inner_prod(e1,t2_n);
        double eG21 = inner_prod(e2,t1_n);
        double eG22 = inner_prod(e2,t2_n);
    
        rTransformationPrestress(0,0) = eG11*eG11;
        rTransformationPrestress(0,1) = eG12*eG12;
        rTransformationPrestress(0,2) = 2.0*eG11*eG12;

        rTransformationPrestress(1,0) = eG21*eG21;
        rTransformationPrestress(1,1) = eG22*eG22;
        rTransformationPrestress(1,2) = 2.0*eG21*eG22;

        rTransformationPrestress(2,0) = eG11*eG21;
        rTransformationPrestress(2,1) = eG12*eG22;
        rTransformationPrestress(2,2) = eG11*eG22+eG12*eG21;      
    }  

    void SupportNonSymmetricNitscheCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_FACTOR];

        KRATOS_WATCH(stabilization_parameter);

        const auto& r_geometry = GetGeometry();

        // Size definitions
        const SizeType number_of_nodes = r_geometry.size();

        const SizeType mat_size = 3 * number_of_nodes;

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

        if (Has(DISPLACEMENT))
        {
            // Integration
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

            // initial determinant of jacobian 
            Vector determinant_jacobian_vector_initial(integration_points.size());
            DeterminantOfJacobianInitial(r_geometry, determinant_jacobian_vector_initial);

            const IntegrationMethod integration_method = r_geometry.GetDefaultIntegrationMethod();
            const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(integration_method);

            const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();
            
            // Prepare memory
            if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
                m_A_ab_covariant_vector.resize(r_number_of_integration_points);
            if (m_dA_vector.size() != r_number_of_integration_points)
                m_dA_vector.resize(r_number_of_integration_points);
            if (m_T_vector.size() != r_number_of_integration_points)
                m_T_vector.resize(r_number_of_integration_points);
            if (m_T_hat_vector.size() != r_number_of_integration_points)
                m_T_hat_vector.resize(r_number_of_integration_points);
            if (m_reference_contravariant_base.size() != r_number_of_integration_points)
                m_reference_contravariant_base.resize(r_number_of_integration_points);
            if (m_n_contravariant_vector.size() != r_number_of_integration_points)
                m_n_contravariant_vector.resize(r_number_of_integration_points);

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
            
                //Compute Kinematics Reference
                KinematicVariables kinematic_variables_reference(
                    r_geometry.WorkingSpaceDimension());

                CalculateKinematics(
                    point_number,
                    kinematic_variables_reference,shape_functions_gradients_i, ConfigurationType::Reference);           
                
                m_A_ab_covariant_vector[point_number] = kinematic_variables_reference.a_ab_covariant;

                m_dA_vector[point_number] = kinematic_variables_reference.dA;

                m_n_contravariant_vector[point_number] = kinematic_variables_reference.n_contravariant;

                CalculateTransformation(kinematic_variables_reference, m_T_vector[point_number], m_T_hat_vector[point_number], m_reference_contravariant_base[point_number]);

                // Compute Kinematics
                KinematicVariables kinematic_variables(
                    r_geometry.WorkingSpaceDimension());
                CalculateKinematics(
                    point_number,
                    kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);
                
                // Create constitutive law parameters:
                ConstitutiveLaw::Parameters constitutive_law_parameters(
                    r_geometry, GetProperties().GetSubProperties().front(), rCurrentProcessInfo);

                ConstitutiveVariables constitutive_variables_membrane(3);

                CalculateConstitutiveVariables(
                    point_number,
                    kinematic_variables,
                    constitutive_variables_membrane,
                    constitutive_law_parameters,
                    ConstitutiveLaw::StressMeasure_PK2);

                //Prestress component
                array_1d<double, 3> prestress = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
                array_1d<double, 3> transformed_prestress;

                Matrix T_pre = ZeroMatrix(3, 3);

                if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
                {
                    CalculateTransformationPrestress(T_pre, kinematic_variables);
                    transformed_prestress = prod(T_pre, prestress);
                }
                else //for isotropic prestress case
                {
                    transformed_prestress = prestress;
                }
                
                constitutive_variables_membrane.StressVector += transformed_prestress;

                // calculate traction vectors
                array_1d<double, 3> traction_vector;

                CalculateTraction(point_number, traction_vector, kinematic_variables, constitutive_variables_membrane);

                // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
                Matrix first_variations_stress_covariant = ZeroMatrix(3, 3*number_of_nodes);

                CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant, kinematic_variables, constitutive_variables_membrane);

                // calculate first variation of traction vectors
                Matrix first_variations_traction = ZeroMatrix(3, 3*number_of_nodes);

                CalculateFirstVariationTraction(point_number, first_variations_traction, first_variations_stress_covariant, kinematic_variables, constitutive_variables_membrane);
                
                //Compute the NURBS basis functions
                Matrix N = r_geometry.ShapeFunctionsValues();
                Matrix r_N = ZeroMatrix(3, 3*number_of_nodes);

                for (IndexType r = 0; r < number_of_nodes; r++)
                {
                    r_N(0, 3 * r) = N(0, r);
                    r_N(1, 3 * r + 1) = N(0, r);
                    r_N(2, 3 * r + 2) = N(0, r);
                }

                //Get the displacement vectors of the previous iteration step
                Vector current_displacement = ZeroVector(mat_size);

                GetValuesVector(current_displacement);

                array_1d<double, 3> displacement_vector;

                displacement_vector = prod(r_N, current_displacement);
                const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);
                displacement_vector -= displacement;

                //Compute the necessary products needed for the second variations of the traction vectors
                Matrix Pi = ZeroMatrix(3, 3);

                CalculateSecondVariationTractionProduct(point_number, Pi, kinematic_variables, constitutive_variables_membrane);

                array_1d<double, 3> second_variations_traction_product_vector = prod(trans(Pi), displacement_vector);

                // calculate second variation of traction vectors
                Matrix second_variations_traction = ZeroMatrix(3 * number_of_nodes, 3 * number_of_nodes);

                CalculateSecondVariationTraction(point_number, second_variations_traction, kinematic_variables, first_variations_stress_covariant, displacement_vector, second_variations_traction_product_vector);

                //Penalty part & RHS
                Matrix H = ZeroMatrix(3, mat_size);
                for (IndexType i = 0; i < number_of_nodes; i++)
                {
                    IndexType index = 3 * i;
                        H(0, index) = N(point_number, i);
                        H(1, index + 1) = N(point_number, i);
                        H(2, index + 2) = N(point_number, i);
                }

                // Differential area
                const double integration_weight = integration_points[point_number].Weight();
                const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
                const double gammaTilde = 1.0;

                // Assembly
                if (CalculateStiffnessMatrixFlag) {

                    noalias(rLeftHandSideMatrix) += (prod(trans(first_variations_traction), H) + prod(trans(H), first_variations_traction))
                        * integration_weight * determinant_jacobian * -gammaTilde;

                    for (IndexType i = 0; i < 3 * number_of_nodes; i++)
                    {
                        for (IndexType j = 0; j < 3 * number_of_nodes; j++)
                        {
                            rLeftHandSideMatrix(i, j) += second_variations_traction(i, j) * integration_weight * determinant_jacobian * -gammaTilde;
                        }
                    }

                    noalias(rLeftHandSideMatrix) += prod(trans(H), H)
                        * integration_weight * determinant_jacobian * stabilization_parameter;
                }

                if (CalculateResidualVectorFlag) {
                    
                    const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);
                    Vector u(mat_size);
                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                        IndexType index = 3 * i;
                        u[index]     = (disp[0] - displacement[0]);
                        u[index + 1] = (disp[1] - displacement[1]);
                        u[index + 2] = (disp[2] - displacement[2]);
                    }
                    
                    noalias(rRightHandSideVector) -= (prod(trans(H), traction_vector)) 
                        * integration_weight * determinant_jacobian * -gammaTilde;
                    noalias(rRightHandSideVector) -= (prod(trans(first_variations_traction), displacement_vector))
                        * integration_weight * determinant_jacobian * -gammaTilde;   
                    noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u)
                        * integration_weight * determinant_jacobian * stabilization_parameter;                  
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportNonSymmetricNitscheCondition::DeterminantOfJacobianInitial(
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

    void SupportNonSymmetricNitscheCondition::CalculateNitscheStabilizationMatrix(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // Size definitions
        const SizeType number_of_nodes = r_geometry.size();

        const SizeType mat_size = 3 * (number_of_nodes);

        // Memory allocation
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
     
        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // initial determinant of jacobian 
        Vector determinant_jacobian_vector_initial(integration_points.size());
        DeterminantOfJacobianInitial(r_geometry, determinant_jacobian_vector_initial);

        const IntegrationMethod integration_method = r_geometry.GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(integration_method);
 
        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();
        
        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);
        if (m_T_hat_vector.size() != r_number_of_integration_points)
            m_T_hat_vector.resize(r_number_of_integration_points);
        if (m_reference_contravariant_base.size() != r_number_of_integration_points)
            m_reference_contravariant_base.resize(r_number_of_integration_points);
        if (m_n_contravariant_vector.size() != r_number_of_integration_points)
            m_n_contravariant_vector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
        
            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference(
                r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference,shape_functions_gradients_i, ConfigurationType::Reference);           

            m_A_ab_covariant_vector[point_number] = kinematic_variables_reference.a_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables_reference.dA;

            m_n_contravariant_vector[point_number] = kinematic_variables_reference.n_contravariant;

            CalculateTransformation(kinematic_variables_reference, m_T_vector[point_number], m_T_hat_vector[point_number], m_reference_contravariant_base[point_number]);

            // Compute Kinematics Actual
            KinematicVariables kinematic_variables(
                r_geometry.WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                r_geometry, GetProperties().GetSubProperties().front(), rCurrentProcessInfo);
            ConstitutiveVariables constitutive_variables_membrane(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            //Prestress component
            array_1d<double, 3> prestress = GetProperties().GetSubProperties().front()[PRESTRESS]*GetProperties().GetSubProperties().front()[THICKNESS];
            array_1d<double, 3> transformed_prestress;

            Matrix T_pre = ZeroMatrix(3, 3);

            if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
            {
                CalculateTransformationPrestress(T_pre, kinematic_variables);
                transformed_prestress = prod(T_pre, prestress);
            }
            else //for isotropic prestress case
            {
                transformed_prestress = prestress;
            }
            
            constitutive_variables_membrane.StressVector += transformed_prestress;

            // calculate traction vectors
            array_1d<double, 3> traction_vector;

            CalculateTraction(point_number, traction_vector, kinematic_variables, constitutive_variables_membrane);

            // calculate the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
            Matrix first_variations_stress_covariant = ZeroMatrix(3, 3*number_of_nodes);

            CalculateFirstVariationStressCovariant(point_number, first_variations_stress_covariant, kinematic_variables, constitutive_variables_membrane);

            // calculate first variation of traction vectors
            Matrix first_variations_traction = ZeroMatrix(3, 3*number_of_nodes);

            CalculateFirstVariationTraction(point_number, first_variations_traction, first_variations_stress_covariant, kinematic_variables, constitutive_variables_membrane);

            //Compute the NURBS basis functions
            Matrix N = r_geometry.ShapeFunctionsValues();
            Matrix r_N = ZeroMatrix(3, 3*number_of_nodes);

            for (IndexType r = 0; r < number_of_nodes; r++)
            {
                r_N(0, 3 * r) = N(0, r);
                r_N(1, 3 * r + 1) = N(0, r);
                r_N(2, 3 * r + 2) = N(0, r);
            }

            //Compute the necessary products needed for the second variations of the traction vectors
            Matrix Pi = ZeroMatrix(3, 3);

            CalculateSecondVariationTractionProduct(point_number, Pi, kinematic_variables, constitutive_variables_membrane);

            array_1d<double, 3> second_variations_traction_product_vector = prod(trans(Pi), traction_vector);

            // calculate second variation of traction vectors
            Matrix second_variations_traction = ZeroMatrix(3 * number_of_nodes, 3 * number_of_nodes);

            CalculateSecondVariationTraction(point_number, second_variations_traction, kinematic_variables, first_variations_stress_covariant, traction_vector, second_variations_traction_product_vector);
           
            // Differential area
            const double integration_weight = integration_points[point_number].Weight();
            const double determinant_jacobian = determinant_jacobian_vector_initial[point_number];
            const double gammaTilde = 1.0;

            // Assembly
            noalias(rLeftHandSideMatrix) += 2 * prod(trans(first_variations_traction), first_variations_traction)
                * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;

            for (IndexType i = 0; i < 3 * number_of_nodes; i++)
            {
                for (IndexType j = 0; j < 3 * number_of_nodes; j++)
                {
                    rLeftHandSideMatrix(i, j) += second_variations_traction(i, j) * integration_weight * determinant_jacobian * gammaTilde * gammaTilde;
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportNonSymmetricNitscheCondition::CalculateTraction(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;
        array_1d<double, 2> n_contravariant_vector;

        stress_vector_covariant = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);
        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

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

    void SupportNonSymmetricNitscheCondition::CalculateFirstVariationStressCovariant(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        const auto& r_geometry = GetGeometry();
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        //Compute the first variation of the Green-Lagrange strains
        Matrix dE_cartesian = ZeroMatrix(3, mat_size);
        Matrix T_patch = ZeroMatrix(3, 3);

        T_patch = m_T_vector[IntegrationPointIndex];

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> dE_curvilinear;
            // strain
            dE_curvilinear[0] = r_DN_De(kr, 0)*rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(kr, 1)*rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(kr, 1));

            dE_cartesian(0, r) = T_patch(0, 0)*dE_curvilinear[0] + T_patch(0, 1)*dE_curvilinear[1] + T_patch(0, 2)*dE_curvilinear[2];
            dE_cartesian(1, r) = T_patch(1, 0)*dE_curvilinear[0] + T_patch(1, 1)*dE_curvilinear[1] + T_patch(1, 2)*dE_curvilinear[2];
            dE_cartesian(2, r) = T_patch(2, 0)*dE_curvilinear[0] + T_patch(2, 1)*dE_curvilinear[1] + T_patch(2, 2)*dE_curvilinear[2];
        }

        //Compute the first variations of the 2nd Piola-Kichhoff stresses in the local Cartesian bases
        Matrix first_variations_stress_cartesian = ZeroMatrix(3, mat_size);
        first_variations_stress_cartesian = prod(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix,dE_cartesian);

        //Transform the first variations of the 2nd Piola-Kichhoff stresses at the covariant bases
        rFirstVariationStressCovariant = prod(m_T_hat_vector[IntegrationPointIndex], first_variations_stress_cartesian);
    }

    void SupportNonSymmetricNitscheCondition::CalculateFirstVariationTraction(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        const auto& r_geometry = GetGeometry();
        
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;
        
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);

        array_1d<double, 2> n_contravariant_vector; 

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

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

        ///2. derivative normal vector * stress covariant
        // Transform the 2nd Piola-kirchhoff stresses in the covariant systems
        array_1d<double, 3> stress_vector_covariant;

        stress_vector_covariant = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.StressVector);

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

    void SupportNonSymmetricNitscheCondition::CalculateSecondVariationTractionProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane)
    {
        array_1d<double, 2> n_contravariant_vector;

        rPi = prod(m_T_hat_vector[IntegrationPointIndex], rThisConstitutiveVariablesMembrane.ConstitutiveMatrix);
        rPi = prod(rPi, m_T_vector[IntegrationPointIndex]);

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        // normal vector * covariant base vector
        Matrix n_a = ZeroMatrix(3, 3);

        for (IndexType r = 0; r < 3; r++)
        {
            n_a (r, 0) = rActualKinematic.a1[r] * n_contravariant_vector[0];
            n_a (r, 1) = rActualKinematic.a2[r] * n_contravariant_vector[1];
            n_a (r, 2) = rActualKinematic.a1[r] * n_contravariant_vector[1] + rActualKinematic.a2[r] * n_contravariant_vector[0];
        }

        rPi = prod(n_a, rPi);
    }

    void SupportNonSymmetricNitscheCondition::CalculateSecondVariationTraction(
        IndexType IntegrationPointIndex,
        Matrix& rSecondVariationTraction,
        const KinematicVariables& rActualKinematic,
        Matrix& rFirstVariationStressCovariant, 
        array_1d<double, 3>& rDisplacement,
        array_1d<double, 3>& rSecondVariationTractionProduct)
    {
        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
 
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

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

        array_1d<double, 2> n_contravariant_vector; 

        n_contravariant_vector = m_n_contravariant_vector[IntegrationPointIndex];

        Vector first_variations_stress_covariant_1 = ZeroVector(mat_size);
        Vector first_variations_stress_covariant_2 = ZeroVector(mat_size);
        Vector first_variations_stress_covariant_3 = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            first_variations_stress_covariant_1(r) = rFirstVariationStressCovariant(0, r);
            first_variations_stress_covariant_2(r) = rFirstVariationStressCovariant(1, r);
            first_variations_stress_covariant_3(r) = rFirstVariationStressCovariant(2, r);
        }

        // displacement * first variation stress covariant
        Matrix displacement_dNCovariant1 = ZeroMatrix(3, mat_size);
        Matrix displacement_dNCovariant2 = ZeroMatrix(3, mat_size);
        Matrix displacement_dNCovariant3 = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            displacement_dNCovariant1(0, r) = rDisplacement(0)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1(1, r) = rDisplacement(1)*first_variations_stress_covariant_1(r);
            displacement_dNCovariant1(2, r) = rDisplacement(2)*first_variations_stress_covariant_1(r);

            displacement_dNCovariant2(0, r) = rDisplacement(0)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2(1, r) = rDisplacement(1)*first_variations_stress_covariant_2(r);
            displacement_dNCovariant2(2, r) = rDisplacement(2)*first_variations_stress_covariant_2(r);

            displacement_dNCovariant3(0, r) = rDisplacement(0)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3(1, r) = rDisplacement(1)*first_variations_stress_covariant_3(r);
            displacement_dNCovariant3(2, r) = rDisplacement(2)*first_variations_stress_covariant_3(r);
        }

        rSecondVariationTraction += prod(trans(r_DN_Dxi), r_DN_Dxi)*rSecondVariationTractionProduct(0);
        rSecondVariationTraction += prod(trans(r_DN_Deta), r_DN_Deta)*rSecondVariationTractionProduct(1);
        rSecondVariationTraction += 0.5*(prod(trans(r_DN_Dxi), r_DN_Deta) + prod(trans(r_DN_Deta), r_DN_Dxi))*rSecondVariationTractionProduct(2);

        rSecondVariationTraction += prod(trans(r_DN_Dxi), displacement_dNCovariant1)*n_contravariant_vector(0);
        rSecondVariationTraction += prod(trans(r_DN_Deta), displacement_dNCovariant2)*n_contravariant_vector(1);
        rSecondVariationTraction += (prod(trans(r_DN_Dxi), displacement_dNCovariant3)*n_contravariant_vector(1) + prod(trans(r_DN_Deta), displacement_dNCovariant3)*n_contravariant_vector(0));

        rSecondVariationTraction += prod(trans(displacement_dNCovariant1), r_DN_Dxi)*n_contravariant_vector(0);
        rSecondVariationTraction += prod(trans(displacement_dNCovariant2), r_DN_Deta)*n_contravariant_vector(1);
        rSecondVariationTraction += (prod(trans(displacement_dNCovariant3), r_DN_Dxi)*n_contravariant_vector(1) + prod(trans(displacement_dNCovariant3), r_DN_Deta)*n_contravariant_vector(0));
    }

    void SupportNonSymmetricNitscheCondition::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const auto& r_geometry = GetGeometry();

        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void SupportNonSymmetricNitscheCondition::AddExplicitContribution(
        const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

        #pragma omp critical
        {
            if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
                // KRATOS_WATCH(rRHS)

                for(SizeType i=0; i< number_of_nodes; ++i) {
                    SizeType index = dimension * i;

                    // KRATOS_WATCH(i)

                    array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                    // KRATOS_WATCH(r_force_residual)

                    for(SizeType j=0; j<dimension; ++j) {
                        AtomicAdd(r_force_residual[j], rRHS[index + j]);
                    }

                    // const array_1d<double, 3> disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

                    // KRATOS_WATCH(rRHS[index + 0])
                    // KRATOS_WATCH(rRHS[index + 1])
                    // KRATOS_WATCH(rRHS[index + 2])

                    // KRATOS_WATCH(disp)

                    // KRATOS_WATCH(r_force_residual)
                }
            }
        }

        KRATOS_CATCH( "" )
    }

    void SupportNonSymmetricNitscheCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        const SizeType number_of_nodes= r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void SupportNonSymmetricNitscheCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry = GetGeometry();

        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos
