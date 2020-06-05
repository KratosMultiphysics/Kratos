//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//                   Tobias Tescheamacher
//                   Riccardo Rossi
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

    void IgaMembraneElement::Initialize()
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
        ProcessInfo& rCurrentProcessInfo,
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
                * m_dA_vector[point_number]
                * GetProperties()[THICKNESS];

            //Define Prestress
            array_1d<double, 3> prestress = GetProperties()[PRESTRESS]*GetProperties()[THICKNESS];

            PrestresstransVariables prestresstrans_variables(3);
            CalculateTransformationmatrixPrestress(
                kinematic_variables,
                prestresstrans_variables 
            );

            array_1d<double, 3> transformed_prestress = prod(prestresstrans_variables.Tpre, prestress);
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
        int step = 0;

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement,step);

        KRATOS_WATCH(current_displacement)

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

    //Prestress Transformation Matrix
    void IgaMembraneElement::CalculateTransformationmatrixPrestress(
        const KinematicVariables& rActualKinematic,
        PrestresstransVariables& rPrestresstransVariables
    )
    {
        //define base vector in reference plane
        //ATTENTION: in some cases the vector must be modified (e.g. catenoid t3_unten=[0, 0, 1])

        array_1d<double, 3> t3_unten;
        t3_unten[0] = 0;
        t3_unten[1] = 1;
        t3_unten[2] = 0;

        array_1d<double, 3> t1_z;
        MathUtils<double>::CrossProduct(t1_z, t3_unten, rActualKinematic.a3);

        array_1d<double, 3> t2;
        MathUtils<double>::CrossProduct(t2, rActualKinematic.a3, t1_z);

        array_1d<double, 3> t1_n = t1_z/norm_2(t1_z);
        array_1d<double, 3> t2_n = t2/norm_2(t2);
        array_1d<double, 3> t3_n = rActualKinematic.a3/norm_2(rActualKinematic.a3);

        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rActualKinematic.a_ab_covariant[0] * rActualKinematic.a_ab_covariant[1]
                - rActualKinematic.a_ab_covariant[2] * rActualKinematic.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant;
        a_ab_contravariant[0] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[1];
        a_ab_contravariant[1] =  inv_det_g_ab * rActualKinematic.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rActualKinematic.a_ab_covariant[2];

        array_1d<double, 3> a_con_1 = rActualKinematic.a1*a_ab_contravariant[0] + rActualKinematic.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_con_2 = rActualKinematic.a1*a_ab_contravariant[2] + rActualKinematic.a2*a_ab_contravariant[1]; 

        //local cartesian coordinates oriented along the 1st base vector in the ref. config.
        double la1 = norm_2(rActualKinematic.a1);
        array_1d<double, 3> e1 = rActualKinematic.a1 / la1;
        double la_con2 = norm_2(a_con_2);
        array_1d<double, 3> e2 = a_con_2 / la_con2;

        //Transformation matrix from the projected basis T to the local cartesian basis
        double eG11 = inner_prod(e1,t1_n);
        double eG12 = inner_prod(e1,t2_n);
        double eG21 = inner_prod(e2,t1_n);
        double eG22 = inner_prod(e2,t2_n);
    
        rPrestresstransVariables.Tpre = ZeroMatrix(3, 3);
        rPrestresstransVariables.Tpre(0,0) = eG11*eG11;
        rPrestresstransVariables.Tpre(0,1) = eG12*eG12;
        rPrestresstransVariables.Tpre(0,2) = 2.0*eG11*eG12;

        rPrestresstransVariables.Tpre(1,0) = eG21*eG21;
        rPrestresstransVariables.Tpre(1,1) = eG22*eG22;
        rPrestresstransVariables.Tpre(1,2) = 2.0*eG21*eG22;

        rPrestresstransVariables.Tpre(2,0) = eG11*eG21;
        rPrestresstransVariables.Tpre(2,1) = eG12*eG22;
        rPrestresstransVariables.Tpre(2,2) = eG11*eG22+eG12*eG21;       
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
    ///@name Dynamic Functions
    ///@{

    void IgaMembraneElement::GetValuesVector(
        Vector& rValues,
        int Step)
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
        int Step)
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
        int Step)
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
        ProcessInfo& rCurrentProcessInfo
    )
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
        ProcessInfo& rCurrentProcessInfo
    )
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

    int IgaMembraneElement::Check(const ProcessInfo& rCurrentProcessInfo)
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

    void IgaMembraneElement::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe)
    {
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian(0, 0) += rDDN_DDe(k, 0)*coords[0];
            Hessian(0, 1) += rDDN_DDe(k, 2)*coords[0];
            Hessian(0, 2) += rDDN_DDe(k, 1)*coords[0];

            Hessian(1, 0) += rDDN_DDe(k, 0)*coords[1];
            Hessian(1, 1) += rDDN_DDe(k, 2)*coords[1];
            Hessian(1, 2) += rDDN_DDe(k, 1)*coords[1];

            Hessian(2, 0) += rDDN_DDe(k, 0)*coords[2];
            Hessian(2, 1) += rDDN_DDe(k, 2)*coords[2];
            Hessian(2, 2) += rDDN_DDe(k, 1)*coords[2];
        }
    }

    ///@}

    void IgaMembraneElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    void IgaMembraneElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    void IgaMembraneElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rValues.size() != r_integration_points.size())
        {
            rValues.resize(r_integration_points.size());
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

        if (rVariable == STRESSES)
        {
            Vector n = ZeroVector(8);
            Vector sigma_top = ZeroVector(3);
            double thickness = GetProperties()[THICKNESS];

            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
            {
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];
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

                //Contravariant metric g_ab_con
                double inv_det_g_ab = 1.0 /
                    (kinematic_variables.a_ab_covariant[0] * kinematic_variables.a_ab_covariant[1]
                        - kinematic_variables.a_ab_covariant[2] * kinematic_variables.a_ab_covariant[2]);

                array_1d<double, 3> a_ab_contravariant;
                a_ab_contravariant[0] =  inv_det_g_ab * kinematic_variables.a_ab_covariant[1];
                a_ab_contravariant[1] =  inv_det_g_ab * kinematic_variables.a_ab_covariant[0];
                a_ab_contravariant[2] = -inv_det_g_ab * kinematic_variables.a_ab_covariant[2];

                //Calculate actual transformation
                Matrix T_actual_vector = ZeroMatrix(r_integration_points.size(),r_integration_points.size());
                Matrix T_hat_actual_vector = ZeroMatrix(r_integration_points.size(),r_integration_points.size());
                array_1d<array_1d<double, 3>,2> actual_contravariant_base;
                CalculateTransformation(kinematic_variables, T_actual_vector, T_hat_actual_vector, actual_contravariant_base);

                double detF = kinematic_variables.dA / m_dA_vector[point_number];

                //Prestress component
                Vector prestress_tensor = ZeroVector(3);
                CalculatePresstressTensor(prestress_tensor,kinematic_variables);

                Vector n_pk2_ca = prod(constitutive_variables_membrane.ConstitutiveMatrix, constitutive_variables_membrane.StrainVector) + prestress_tensor;
                
                //Transform PK2 stress into Cauchy stress
                //Deformation Gradient F
                Matrix deformation_gradient = ZeroMatrix(2);
                double det_deformation_gradient;

                array_1d<Vector,2> rCurrentCovariantBase;
                rCurrentCovariantBase[0] = kinematic_variables.a1;
                rCurrentCovariantBase[1] = kinematic_variables.a2;

                array_1d< array_1d<double, 3>,2> rReferenceContraVariantBase = m_reference_contravariant_base[point_number];

                for (SizeType i=0;i<2;++i){
                    for (SizeType j=0;j<2;++j){
                        for (SizeType k=0;k<2;++k){
                            deformation_gradient(j,k) += rCurrentCovariantBase[i][j]*rReferenceContraVariantBase[i][k];
                        }
                    }
                }

                Matrix stress_matrix = MathUtils<double>::StressVectorToTensor(n_pk2_ca);
                Matrix temp_stress_matrix = prod(deformation_gradient,stress_matrix);
                Matrix temp_stress_matrix_2 = prod(temp_stress_matrix,trans(deformation_gradient));
                Matrix cauchy_stress_matrix = temp_stress_matrix_2 / detF;
   
                Vector n_cau = MathUtils<double>::StressTensorToVector(cauchy_stress_matrix,3);
                
                //PK2 stress
                // sigma_top(0) = n_pk2_ca[0];
                // sigma_top(1) = n_pk2_ca[1];
                // sigma_top(2) = n_pk2_ca[2];

                //Cauchy stress
                sigma_top(0) = n_cau[0];
                sigma_top(1) = n_cau[1];
                sigma_top(2) = n_cau[2];

                rValues[point_number] = sigma_top;
            }
        }
        else if (rVariable == EXTERNAL_FORCES_VECTOR) {

            // Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
            // Vector external_forces_vector = ZeroVector(3);

            // for (int i = 0; i < number_of_nodes; i++)
            // {

            //     const NodeType & iNode = GetGeometry()[i];
            //     const Vector& forces = iNode.GetValue(EXTERNAL_FORCES_VECTOR);

            //     external_forces_vector[0] += N[i] * forces[0];
            //     external_forces_vector[1] += N[i] * forces[1];
            //     external_forces_vector[2] += N[i] * forces[2];
            // }

            // rValues[0] = external_forces_vector;
        }
        else
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            { 
                rValues[point_number] = ZeroVector(3);
            }
        }
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
                // const double integration_weight_i = r_integration_points[point_number].Weight();
                const Matrix& shape_functions_gradients_i = r_shape_functions_gradients[point_number];

                // CovariantBaseVectors(base_vectors_current_cov,shape_functions_gradients_i,ConfigurationType::Reference);
                // base_1 += base_vectors_current_cov[0]*integration_weight_i;
                // base_2 += base_vectors_current_cov[1]*integration_weight_i;

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

    //Prestress postprocessing
    void IgaMembraneElement::CalculatePresstressTensor(
        Vector& rPrestressTensor,
        KinematicVariables& rActualKinematic
    )
    {
        rPrestressTensor.resize(3);

        array_1d<double, 3> prestress = GetProperties()[PRESTRESS]*GetProperties()[THICKNESS];

        PrestresstransVariables prestresstrans_variables(3);
        CalculateTransformationmatrixPrestress(
                rActualKinematic,
                prestresstrans_variables 
        );

        rPrestressTensor = prod(prestresstrans_variables.Tpre, prestress);
    } 

} // Namespace Kratos