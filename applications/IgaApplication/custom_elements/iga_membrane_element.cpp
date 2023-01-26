//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/iga_membrane_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void IgaMembraneElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void IgaMembraneElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Assembly
    ///@{

    void IgaMembraneElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);

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

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

            //Compute Kinematics Reference
            KinematicVariables kinematic_variables_reference(
                GetGeometry().WorkingSpaceDimension());

            CalculateKinematics(
                point_number,
                kinematic_variables_reference,shape_functions_gradients_i, ConfigurationType::Reference);
            
            m_A_ab_covariant_vector[point_number] = kinematic_variables_reference.a_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables_reference.dA;

            CalculateTransformation(kinematic_variables_reference, m_T_vector[point_number], m_T_hat_vector[point_number], m_reference_contravariant_base[point_number]);

            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Matrix BMembrane = ZeroMatrix(3, mat_size);
            CalculateBMembrane(
                point_number,
                BMembrane,
                kinematic_variables);

            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            CalculateSecondVariationStrain(
                point_number,
                second_variations_strain,
                kinematic_variables);

            double integration_weight =
                r_integration_points[point_number].Weight()
                * m_dA_vector[point_number];

            //Prestress component
            array_1d<double, 3> prestress = GetProperties()[PRESTRESS]*GetProperties()[THICKNESS];
            array_1d<double, 3> transformed_prestress;

            Matrix T_pre = ZeroMatrix(3, 3);
            if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
            {
                CalculateTransformationPrestress(T_pre, kinematic_variables_reference);
                transformed_prestress = prod(T_pre, prestress);
            }
            else //for isotropic prestress case
            {
                transformed_prestress = prestress;
            }
            
            constitutive_variables_membrane.StressVector += transformed_prestress;

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                //adding membrane contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BMembrane,
                    constitutive_variables_membrane.ConstitutiveMatrix,
                    integration_weight);
                // adding  non-linear-contribution to Stiffness-Matrix
                CalculateAndAddNonlinearKm(
                    rLeftHandSideMatrix,
                    second_variations_strain,
                    constitutive_variables_membrane.StressVector,
                    integration_weight);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.StressVector);
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Kinematics
    ///@{

    void IgaMembraneElement::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ConfigurationType& rConfiguration
    )
    {
        // pass/call this ShapeFunctionsLocalGradients[pnt]
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType number_of_nodes = GetGeometry().size();
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 0);
            g1[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 0);
            g1[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 0);

            g2[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * rShapeFunctionGradientValues(i, 1);
            g2[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * rShapeFunctionGradientValues(i, 1);
            g2[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * rShapeFunctionGradientValues(i, 1);
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

    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void IgaMembraneElement::CalculateTransformation(
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

        // e * a_covariant
        Matrix G_hat = ZeroMatrix(2, 2);
        G_hat(0, 0) = inner_prod(e1, rKinematicVariables.a1);
        G_hat(0, 1) = inner_prod(e1, rKinematicVariables.a2);
        G_hat(1, 0) = inner_prod(e2, rKinematicVariables.a1);
        G_hat(1, 1) = inner_prod(e2, rKinematicVariables.a2);

        //Transformation matrix T_hat
        if (rT_hat.size1() != 3 && rT_hat.size2() != 3)
            rT_hat.resize(3, 3);
        noalias(rT_hat) = ZeroMatrix(3, 3);

        rT_hat(0, 0) = pow(G_hat(0, 0), 2);
        rT_hat(0, 1) = pow(G_hat(1, 0), 2);
        rT_hat(0, 2) = G_hat(0, 0) * G_hat(1, 0);

        rT_hat(1, 0) = pow(G_hat(0, 1), 2);
        rT_hat(1, 1) = pow(G_hat(1, 1), 2);
        rT_hat(1, 2) = G_hat(0, 1) * G_hat(1, 1);

        rT_hat(2, 0) = G_hat(0, 0) * G_hat(0, 1);
        rT_hat(2, 1) = G_hat(1, 0) * G_hat(1, 1);
        rT_hat(2, 2) = (G_hat(0, 0) * G_hat(1, 1) + G_hat(1, 0) * G_hat(0, 1));
    }

    //Prestress Transformation Matrix
    void IgaMembraneElement::CalculateTransformationPrestress(
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

    void IgaMembraneElement::CalculateConstitutiveVariables(
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

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        rThisConstitutiveVariablesMembrane.ConstitutiveMatrix *= GetProperties()[THICKNESS];

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
   }

    void IgaMembraneElement::CalculateBMembrane(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic)
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

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

            rB(0, r) = m_T_vector[IntegrationPointIndex](0, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](0, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](0, 2)*dE_curvilinear[2];
            rB(1, r) = m_T_vector[IntegrationPointIndex](1, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](1, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](1, 2)*dE_curvilinear[2];
            rB(2, r) = m_T_vector[IntegrationPointIndex](2, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](2, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](2, 2)*dE_curvilinear[2];
        }
    }

    void IgaMembraneElement::CalculateSecondVariationStrain(
        IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        const KinematicVariables& rActualKinematic)
    {
        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
 
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        // second variation of strain w.r.t. dofs
        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            for (IndexType s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                IndexType ks = s / 3;
                IndexType dirs = s % 3;

                // strain
                array_1d<double, 3> ddE_cu = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_cu[0] = r_DN_De(kr, 0) * r_DN_De(ks, 0);
                    ddE_cu[1] = r_DN_De(kr, 1) * r_DN_De(ks, 1);
                    ddE_cu[2] = 0.5 * (r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(kr, 1) * r_DN_De(ks, 0));

                    rSecondVariationsStrain.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddE_cu[2];
                }
            }
        }
    }

    ///@}
    ///@name Implicit
    ///@{

    void IgaMembraneElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        // Rayleigh Damping Matrix: alpha*M + beta*K

        // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if (GetProperties().Has(RAYLEIGH_ALPHA))
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

        double beta = 0.0;
        if (GetProperties().Has(RAYLEIGH_BETA))
            beta = GetProperties()[RAYLEIGH_BETA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        // 2.-Calculate StiffnessMatrix and MassMatrix:
        if (std::abs(alpha) < 1E-12 && std::abs(beta) < 1E-12) {
            // no damping specified, only setting the matrix to zero
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType mat_size = number_of_nodes * 3;
            if (rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size) {
                rDampingMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);
        } else if (std::abs(alpha) > 1E-12 && std::abs(beta) < 1E-12) {
            // damping only required with the mass matrix
            CalculateMassMatrix(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= alpha;
        } else if (std::abs(alpha) < 1E-12 && std::abs(beta) > 1E-12) {
            // damping only required with the stiffness matrix
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;
        } else {
            // damping with both mass matrix and stiffness matrix required
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;

            Matrix mass_matrix;
            CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
            noalias(rDampingMatrix) += alpha  * mass_matrix;
        }

        KRATOS_CATCH("")
    }

    void IgaMembraneElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        // Shape function values for all integration points
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            double integration_weight = r_integration_points[point_number].Weight();

            double thickness = this->GetProperties().GetValue(THICKNESS);
            double density = this->GetProperties().GetValue(DENSITY);
            double mass = thickness * density * m_dA_vector[point_number] * integration_weight;

            if (rMassMatrix.size1() != mat_size)
                rMassMatrix.resize(mat_size, mat_size, false);

            rMassMatrix = ZeroMatrix(mat_size, mat_size);

            for (unsigned int r = 0; r<number_of_nodes; r++)
            {
                for (unsigned int s = 0; s<number_of_nodes; s++)
                {
                    rMassMatrix(3 * s, 3 * r) = r_N(point_number, s)*r_N(point_number, r) * mass;
                    rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                    rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
                }
            }
        }
        KRATOS_CATCH("")
    }

    ///@}
    ///@name Stiffness matrix assembly
    ///@{

    inline void IgaMembraneElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight
    )
    {
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(rB), Matrix(prod(rD, rB)));
    }

    inline void IgaMembraneElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& rSecondVariationsStrain,
        const Vector& rSD,
        const double IntegrationWeight)
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        for (IndexType n = 0; n < mat_size; n++)
        {
            for (IndexType m = 0; m <= n; m++)
            {
                double nm = (rSD[0] * rSecondVariationsStrain.B11(n, m)
                    + rSD[1] * rSecondVariationsStrain.B22(n, m)
                    + rSD[2] * rSecondVariationsStrain.B12(n, m)) * IntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if (n != m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
    }

    ///@}
    ///@name Postprocessing
    ///@{

    void IgaMembraneElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        if (rValues.size() != r_integration_points.size())
        {
            rValues.resize(r_integration_points.size());
        }

        if (rVariable==PRINCIPAL_STRESS_1 || rVariable==PRINCIPAL_STRESS_2)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
            {
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
                KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());

                array_1d<double, 3> pk2_stresses;
                CalculatePK2Stresses(point_number, pk2_stresses, kinematic_variables, shape_functions_gradients_i, rCurrentProcessInfo);

                if (rVariable==PRINCIPAL_STRESS_1)
                {
                    double principal_stresses_1 = 0.5*(pk2_stresses[0] + pk2_stresses[1] + sqrt(pow(pk2_stresses[0] - pk2_stresses[1], 2)+4*pow(pk2_stresses[2], 2)));
                    rValues[point_number] = principal_stresses_1;
                }

                if (rVariable==PRINCIPAL_STRESS_2)
                {
                    double principal_stresses_2 = 0.5*(pk2_stresses[0] + pk2_stresses[1] - sqrt(pow(pk2_stresses[0] - pk2_stresses[1], 2)+4*pow(pk2_stresses[2], 2)));
                    rValues[point_number] = principal_stresses_2;
                }
            }
        }
        else
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            { 
                rValues[point_number] = 0.0;
            }
        }
    }

    void IgaMembraneElement::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        if (rValues.size() != r_integration_points.size())
        {
            rValues.resize(r_integration_points.size());
        }

        if (rVariable==PK2_STRESS_VECTOR || rVariable==CAUCHY_STRESS_VECTOR)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
            {
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
                KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());

                if (rVariable==PK2_STRESS_VECTOR)
                {
                    array_1d<double, 3> pk2_stresses = ZeroVector(3);
                    CalculatePK2Stresses(point_number, pk2_stresses, kinematic_variables, shape_functions_gradients_i, rCurrentProcessInfo);

                    rValues[point_number] = pk2_stresses;
                }

                if (rVariable==CAUCHY_STRESS_VECTOR)
                {
                    array_1d<double, 3> cauchy_stresses = ZeroVector(3);
                    CalculateCauchyStresses(point_number, cauchy_stresses, kinematic_variables, shape_functions_gradients_i, rCurrentProcessInfo);

                    rValues[point_number] = cauchy_stresses;
                } 
            }
        }
        else
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            { 
                rValues[point_number] = ZeroVector(3);
            }
        }
    }
    
    void IgaMembraneElement::CalculatePK2Stresses(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rPK2Stresses,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        // Compute Kinematics and Metric
        CalculateKinematics(
            IntegrationPointIndex,
            rKinematicVariables,rShapeFunctionGradientValues, ConfigurationType::Current);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        ConstitutiveVariables constitutive_variables_membrane(3);
        CalculateConstitutiveVariables(
            IntegrationPointIndex,
            rKinematicVariables,
            constitutive_variables_membrane,
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_PK2);

        //Prestress component
        array_1d<double, 3> prestress = GetProperties()[PRESTRESS]*GetProperties()[THICKNESS];
        array_1d<double, 3> prestress_tensor;

        Matrix T_pre = ZeroMatrix(3, 3);
        if (Has(LOCAL_PRESTRESS_AXIS_1)) //for anisotropic prestress case
        {
            CalculateTransformationPrestress(T_pre, rKinematicVariables);
            prestress_tensor = prod(T_pre, prestress);
        }
        else //for isotropic prestress case
        {
            prestress_tensor = prestress;
        }

        rPK2Stresses = prod(constitutive_variables_membrane.ConstitutiveMatrix, constitutive_variables_membrane.StrainVector) + prestress_tensor;
    }

    void IgaMembraneElement::CalculateCauchyStresses(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rCauchyStresses,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        array_1d<double, 3> pk2_stresses = ZeroVector(3);
        CalculatePK2Stresses(IntegrationPointIndex, pk2_stresses, rKinematicVariables, rShapeFunctionGradientValues, rCurrentProcessInfo);

        //Transform PK2 stress into Cauchy stress
        //Deformation Gradient F
        Matrix deformation_gradient = ZeroMatrix(2);

        array_1d<Vector,2> rCurrentCovariantBase;
        rCurrentCovariantBase[0] = rKinematicVariables.a1;
        rCurrentCovariantBase[1] = rKinematicVariables.a2;

        array_1d< array_1d<double, 3>,2> rReferenceContraVariantBase = m_reference_contravariant_base[IntegrationPointIndex];

        for (SizeType i=0;i<2;++i){
            for (SizeType j=0;j<2;++j){
                for (SizeType k=0;k<2;++k){
                    deformation_gradient(j,k) += rCurrentCovariantBase[i][j]*rReferenceContraVariantBase[i][k];
                }
            }
        }

        double det_deformation_gradient = deformation_gradient(0,0)*deformation_gradient(1,1) - deformation_gradient(0,1)*deformation_gradient(1,0);

        Matrix stress_matrix = MathUtils<double>::StressVectorToTensor(pk2_stresses);
        Matrix temp_stress_matrix = prod(deformation_gradient,stress_matrix);
        Matrix temp_stress_matrix_2 = prod(temp_stress_matrix,trans(deformation_gradient));
        Matrix cauchy_stress_matrix = temp_stress_matrix_2 / det_deformation_gradient;

        rCauchyStresses = MathUtils<double>::StressTensorToVector(cauchy_stress_matrix,3);
    }

    void IgaMembraneElement::Calculate(
        const Variable<Matrix>& rVariable, 
        Matrix& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
            rOutput = ZeroMatrix(3);
            array_1d<Vector,2> base_vectors_current_cov;
            const IntegrationMethod integration_method = GetGeometry().GetDefaultIntegrationMethod();
            const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients = GetGeometry().ShapeFunctionsLocalGradients(integration_method);
            const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(integration_method);

            Vector base_1 = ZeroVector(3);
            Vector base_2 = ZeroVector(3);
            Vector base_3 = ZeroVector(3);

            for (SizeType point_number = 0; point_number < r_integration_points.size(); ++point_number){
                
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
                const double integration_weight_i = r_integration_points[point_number].Weight();

                // Compute Kinematics and Metric
                KinematicVariables kinematic_variables(
                    GetGeometry().WorkingSpaceDimension());
                CalculateKinematics(
                    point_number,
                    kinematic_variables,shape_functions_gradients_i, ConfigurationType::Current);
                
                base_1 += kinematic_variables.a1*integration_weight_i;
                base_2 += kinematic_variables.a2*integration_weight_i;
            }

            MathUtils<double>::CrossProduct(base_3, base_1, base_2);
            base_3 /= MathUtils<double>::Norm(base_3);

            column(rOutput,0) = base_1;
            column(rOutput,1) = base_2;
            column(rOutput,2) = base_3;
        }
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void IgaMembraneElement::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void IgaMembraneElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const IndexType index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void IgaMembraneElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const IndexType index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void IgaMembraneElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void IgaMembraneElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int IgaMembraneElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
                << "THICKNESS not provided for element " << this->Id() << std::endl;

            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;
        }

        return 0;
    }

    ///@}

} // Namespace Kratos