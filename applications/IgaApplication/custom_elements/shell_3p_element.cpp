//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_3p_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void Shell3pElement::Initialize()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_B_ab_covariant_vector.size() != r_number_of_integration_points)
            m_B_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(
            WorkingSpaceDimension());

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateKinematics(
                point_number,
                kinematic_variables);

            m_A_ab_covariant_vector[point_number] = kinematic_variables.a_ab_covariant;
            m_B_ab_covariant_vector[point_number] = kinematic_variables.b_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables.dA;

            CalculateTransformation(kinematic_variables, m_T_vector[point_number]);
        }

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void Shell3pElement::InitializeMaterial()
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

    void Shell3pElement::CalculateAll(
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

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);
            ConstitutiveVariables constitutive_variables_curvature(3);
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Matrix BMembrane = ZeroMatrix(3, mat_size);
            Matrix BCurvature = ZeroMatrix(3, mat_size);
            CalculateBMembrane(
                point_number,
                BMembrane,
                kinematic_variables);
            CalculateBCurvature(
                point_number,
                BCurvature,
                kinematic_variables);

            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrainCurvature(
                point_number,
                second_variations_strain,
                second_variations_curvature,
                kinematic_variables);

            double integration_weight =
                r_integration_points[point_number].Weight()
                * m_dA_vector[point_number]
                * GetProperties()[THICKNESS];

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                //adding membrane contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BMembrane,
                    constitutive_variables_membrane.ConstitutiveMatrix,
                    integration_weight);
                //adding curvature contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BCurvature,
                    constitutive_variables_curvature.ConstitutiveMatrix,
                    integration_weight);

                // adding  non-linear-contribution to Stiffness-Matrix
                CalculateAndAddNonlinearKm(
                    rLeftHandSideMatrix,
                    second_variations_strain,
                    constitutive_variables_membrane.StressVector,
                    integration_weight);

                CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                    second_variations_curvature,
                    constitutive_variables_curvature.StressVector,
                    integration_weight);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.StressVector);
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.StressVector);
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Kinematics
    ///@{

    void Shell3pElement::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    )
    {
        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);

        rKinematicVariables.a1 = column(J, 0);
        rKinematicVariables.a2 = column(J, 1);

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

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        rKinematicVariables.b_ab_covariant[0] = H(0, 0)*rKinematicVariables.a3_tilde[0] + H(1, 0)*rKinematicVariables.a3_tilde[1] + H(2, 0)*rKinematicVariables.a3_tilde[2];
        rKinematicVariables.b_ab_covariant[1] = H(0, 1)*rKinematicVariables.a3_tilde[0] + H(1, 1)*rKinematicVariables.a3_tilde[1] + H(2, 1)*rKinematicVariables.a3_tilde[2];
        rKinematicVariables.b_ab_covariant[2] = H(0, 2)*rKinematicVariables.a3_tilde[0] + H(1, 2)*rKinematicVariables.a3_tilde[1] + H(2, 2)*rKinematicVariables.a3_tilde[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void Shell3pElement::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT
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
    }

    void Shell3pElement::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - m_A_ab_covariant_vector[IntegrationPointIndex]);
        noalias(rThisConstitutiveVariablesMembrane.StrainVector) = prod(m_T_vector[IntegrationPointIndex], strain_vector);

        array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector[IntegrationPointIndex];
        noalias(rThisConstitutiveVariablesCurvature.StrainVector) = prod(m_T_vector[IntegrationPointIndex], curvature_vector);

        // Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
        noalias(rThisConstitutiveVariablesCurvature.StressVector) = prod(
            trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);
    }

    void Shell3pElement::CalculateBMembrane(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic)
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size1() != mat_size || rB.size2() != mat_size)
            rB.resize(mat_size, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> dE_curvilinear;
            // strain
            dE_curvilinear[0] = r_DN_De(0, kr)*rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(1, kr)*rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5*(r_DN_De(0, kr)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(1, kr));

            rB(0, r) = m_T_vector[IntegrationPointIndex](0, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](0, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](0, 2)*dE_curvilinear[2];
            rB(1, r) = m_T_vector[IntegrationPointIndex](1, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](1, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](1, 2)*dE_curvilinear[2];
            rB(2, r) = m_T_vector[IntegrationPointIndex](2, 0)*dE_curvilinear[0] + m_T_vector[IntegrationPointIndex](2, 1)*dE_curvilinear[1] + m_T_vector[IntegrationPointIndex](2, 2)*dE_curvilinear[2];
        }
    }

    void Shell3pElement::CalculateBCurvature(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic)
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        double inv_dA = 1 / rActualKinematic.dA;
        double inv_dA3 = 1 / std::pow(rActualKinematic.dA, 3);

        for (int i = 0; i < number_of_control_points; i++)
        {
            unsigned int index = 3 * i;
            //first line
            da3(0, 0) =  0;
            da3(0, 1) = -r_DN_De(0, i) * rActualKinematic.a2[2] + r_DN_De(1, i) * rActualKinematic.a1[2];
            da3(0, 2) =  r_DN_De(0, i) * rActualKinematic.a2[1] - r_DN_De(1, i) * rActualKinematic.a1[1];

            //second line
            da3(1, 0) =  r_DN_De(0, i) * rActualKinematic.a2[2] - r_DN_De(1, i) * rActualKinematic.a1[2];
            da3(1, 1) =  0;
            da3(1, 2) = -r_DN_De(0, i) * rActualKinematic.a2[0] + r_DN_De(1, i) * rActualKinematic.a1[0];

            //third line
            da3(2, 0) = -r_DN_De(0, i) * rActualKinematic.a2[1] + r_DN_De(1, i) * rActualKinematic.a1[1];
            da3(2, 1) =  r_DN_De(0, i) * rActualKinematic.a2[0] - r_DN_De(1, i) * rActualKinematic.a1[0];
            da3(2, 2) =  0;

            for (unsigned int j = 0; j < 3; j++)
            {
                double a3da3la3 = (rActualKinematic.a3_tilde[0] * da3(j, 0) + rActualKinematic.a3_tilde[1] * da3(j, 1) + rActualKinematic.a3_tilde[2] * da3(j, 2)) * inv_dA3;

                dn(j, 0) = da3(j, 0) * inv_dA - rActualKinematic.a3_tilde[0] * a3da3la3;
                dn(j, 1) = da3(j, 1) * inv_dA - rActualKinematic.a3_tilde[1] * a3da3la3;
                dn(j, 2) = da3(j, 2) * inv_dA - rActualKinematic.a3_tilde[2] * a3da3la3;
            }

            Matrix H = ZeroMatrix(3, 3);
            CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

            // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
            b(0, index)     = 0 - (r_DDN_DDe(0, i) * rActualKinematic.a3[0] + H(0, 0) * dn(0, 0) + H(1, 0) * dn(0, 1) + H(2, 0) * dn(0, 2));
            b(0, index + 1) = 0 - (r_DDN_DDe(0, i) * rActualKinematic.a3[1] + H(0, 0) * dn(1, 0) + H(1, 0) * dn(1, 1) + H(2, 0) * dn(1, 2));
            b(0, index + 2) = 0 - (r_DDN_DDe(0, i) * rActualKinematic.a3[2] + H(0, 0) * dn(2, 0) + H(1, 0) * dn(2, 1) + H(2, 0) * dn(2, 2));

            //second line
            b(1, index)     = 0 - (r_DDN_DDe(1, i) * rActualKinematic.a3[0] + H(0, 1) * dn(0, 0) + H(1, 1) * dn(0, 1) + H(2, 1) * dn(0, 2));
            b(1, index + 1) = 0 - (r_DDN_DDe(1, i) * rActualKinematic.a3[1] + H(0, 1) * dn(1, 0) + H(1, 1) * dn(1, 1) + H(2, 1) * dn(1, 2));
            b(1, index + 2) = 0 - (r_DDN_DDe(1, i) * rActualKinematic.a3[2] + H(0, 1) * dn(2, 0) + H(1, 1) * dn(2, 1) + H(2, 1) * dn(2, 2));

            //third line
            b(2, index)     = 0 - (r_DDN_DDe(2, i) * rActualKinematic.a3[0] + H(0, 2) * dn(0, 0) + H(1, 2) * dn(0, 1) + H(2, 2) * dn(0, 2));
            b(2, index + 1) = 0 - (r_DDN_DDe(2, i) * rActualKinematic.a3[1] + H(0, 2) * dn(1, 0) + H(1, 2) * dn(1, 1) + H(2, 2) * dn(1, 2));
            b(2, index + 2) = 0 - (r_DDN_DDe(2, i) * rActualKinematic.a3[2] + H(0, 2) * dn(2, 0) + H(1, 2) * dn(2, 1) + H(2, 2) * dn(2, 2));
        }

        noalias(rB) = -prod(m_T_vector[IntegrationPointIndex], b);

        KRATOS_CATCH("")
    }

    void Shell3pElement::CalculateSecondVariationStrainCurvature(
        IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const KinematicVariables& rActualKinematic)
    {
        const auto& r_geometry = GetGeometry();

        const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        double l_a3 = norm_2(rActualKinematic.a3_tilde);
        double l_a3_3 = pow(l_a3, 3);
        double l_a3_5 = pow(l_a3, 5);
        double inv_l_a3 = 1 / l_a3;
        double inv_l_a3_3 = 1 / l_a3_3;
        double inv_l_a3_5 = 1 / l_a3_5;

        Matrix S_da3 = ZeroMatrix(3, mat_size);
        Vector S_a3_da3 = ZeroVector(mat_size);
        Vector S_a3_da3_l_a3_3 = ZeroVector(mat_size);
        Matrix S_dn = ZeroMatrix(3, mat_size);

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // first variation of strain and curvature w.r.t. dof
        for (int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            array_1d<double, 3> S_dg_1 = ZeroVector(3);
            array_1d<double, 3> S_dg_2 = ZeroVector(3);
            S_dg_1(dirr) = r_DN_De(0, kr);
            S_dg_2(dirr) = r_DN_De(1, kr);

            // curvature
            S_da3(0, r) = S_dg_1(1) * rActualKinematic.a2(2) - S_dg_1(2) * rActualKinematic.a2(1) + rActualKinematic.a1(1) * S_dg_2(2) - rActualKinematic.a1(2) * S_dg_2(1);
            S_da3(1, r) = S_dg_1(2) * rActualKinematic.a2(0) - S_dg_1(0) * rActualKinematic.a2(2) + rActualKinematic.a1(2) * S_dg_2(0) - rActualKinematic.a1(0) * S_dg_2(2);
            S_da3(2, r) = S_dg_1(0) * rActualKinematic.a2(1) - S_dg_1(1) * rActualKinematic.a2(0) + rActualKinematic.a1(0) * S_dg_2(1) - rActualKinematic.a1(1) * S_dg_2(0);

            S_a3_da3[r] = rActualKinematic.a3_tilde[0] * S_da3(0, r) + rActualKinematic.a3_tilde[1] * S_da3(1, r) + rActualKinematic.a3_tilde[2] * S_da3(2, r);
            S_a3_da3_l_a3_3[r] = S_a3_da3[r] * inv_l_a3_3;

            S_dn(0, r) = S_da3(0, r) * inv_l_a3 - rActualKinematic.a3_tilde[0] * S_a3_da3_l_a3_3[r];
            S_dn(1, r) = S_da3(1, r) * inv_l_a3 - rActualKinematic.a3_tilde[1] * S_a3_da3_l_a3_3[r];
            S_dn(2, r) = S_da3(2, r) * inv_l_a3 - rActualKinematic.a3_tilde[2] * S_a3_da3_l_a3_3[r];
        }

        // second variation of strain and curvature w.r.t. dofs
        for (int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            for (int s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                int ks = s / 3;
                int dirs = s % 3;

                // strain
                array_1d<double, 3> ddE_cu = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_cu[0] = r_DN_De(0, kr) * r_DN_De(0, ks);
                    ddE_cu[1] = r_DN_De(1, kr) * r_DN_De(1, ks);
                    ddE_cu[2] = 0.5 * (r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(1, kr) * r_DN_De(0, ks));

                    rSecondVariationsStrain.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddE_cu[2];
                    rSecondVariationsStrain.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddE_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddE_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddE_cu[2];
                }

                // curvature
                array_1d<double, 3> dda3 = ZeroVector(3);
                int dirt = 4 - dirr - dirs;
                int ddir = dirr - dirs;
                     if (ddir == -1) dda3(dirt - 1) =  r_DN_De(0, kr) * r_DN_De(1, ks) - r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir ==  2) dda3(dirt - 1) =  r_DN_De(0, kr) * r_DN_De(1, ks) - r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir ==  1) dda3(dirt - 1) = -r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(0, ks) * r_DN_De(1, kr);
                else if (ddir == -2) dda3(dirt - 1) = -r_DN_De(0, kr) * r_DN_De(1, ks) + r_DN_De(0, ks) * r_DN_De(1, kr);

                double c = -(dda3[0] * rActualKinematic.a3_tilde[0] + dda3[1] * rActualKinematic.a3_tilde[1] + dda3[2] * rActualKinematic.a3_tilde[2]
                    + S_da3(0, r) * S_da3(0, s) + S_da3(1, r) * S_da3(1, s) + S_da3(2, r) * S_da3(2, s)
                    ) * inv_l_a3_3;

                double d = 3.0 * S_a3_da3[r] * S_a3_da3[s] * inv_l_a3_5;

                array_1d<double, 3> ddn = ZeroVector(3);
                ddn[0] = dda3[0] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(0, r) - S_a3_da3_l_a3_3[r] * S_da3(0, s) + (c + d) * rActualKinematic.a3_tilde[0];
                ddn[1] = dda3[1] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(1, r) - S_a3_da3_l_a3_3[r] * S_da3(1, s) + (c + d) * rActualKinematic.a3_tilde[1];
                ddn[2] = dda3[2] * inv_l_a3 - S_a3_da3_l_a3_3[s] * S_da3(2, r) - S_a3_da3_l_a3_3[r] * S_da3(2, s) + (c + d) * rActualKinematic.a3_tilde[2];

                array_1d<double, 3> ddK_cu = ZeroVector(3);
                ddK_cu[0] = r_DDN_DDe(0, kr) * S_dn(dirr, s) + r_DDN_DDe(0, ks) * S_dn(dirs, r)
                    + H(0, 0) * ddn[0] + H(1, 0) * ddn[1] + H(2, 0) * ddn[2];
                ddK_cu[1] = r_DDN_DDe(1, kr) * S_dn(dirr, s) + r_DDN_DDe(1, ks) * S_dn(dirs, r)
                    + H(0, 1) * ddn[0] + H(1, 1) * ddn[1] + H(2, 1) * ddn[2];
                ddK_cu[2] = r_DDN_DDe(2, kr) * S_dn(dirr, s) + r_DDN_DDe(2, ks) * S_dn(dirs, r)
                    + H(0, 2) * ddn[0] + H(1, 2) * ddn[1] + H(2, 2) * ddn[2];

                rSecondVariationsCurvature.B11(r, s) = m_T_vector[IntegrationPointIndex](0, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](0, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](0, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B11(s, r) = rSecondVariationsCurvature.B11(r, s);
                rSecondVariationsCurvature.B22(r, s) = m_T_vector[IntegrationPointIndex](1, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](1, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](1, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B22(s, r) = rSecondVariationsCurvature.B22(r, s);
                rSecondVariationsCurvature.B12(r, s) = m_T_vector[IntegrationPointIndex](2, 0) * ddK_cu[0] + m_T_vector[IntegrationPointIndex](2, 1) * ddK_cu[1] + m_T_vector[IntegrationPointIndex](2, 2) * ddK_cu[2];
                rSecondVariationsCurvature.B12(s, r) = rSecondVariationsCurvature.B12(r, s);
            }
        }
    }

    ///@}
    ///@name Stiffness matrix assembly
    ///@{

    inline void Shell3pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight
    )
    {
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(rB), Matrix(prod(rD, rB)));
    }

    inline void Shell3pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& rSecondVariationsStrain,
        const Vector& rSD,
        const double IntegrationWeight)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
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

    void Shell3pElement::GetValuesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const int index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void Shell3pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void Shell3pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void Shell3pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void Shell3pElement::GetDofList(
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

    int Shell3pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

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

    void Shell3pElement::CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe)
    {
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (size_t k = 0; k < number_of_points; k++)
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

} // Namespace Kratos


