//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Michael Loibl
//


// System includes

// External includes

// Project includes
#include "custom_elements/iga_shell_3p_element.h"


namespace Kratos
{
    void IgaShell3pElement::Initialize()
    {
        KRATOS_TRY

        InitializeMaterial();

        MetricVariables initial_metric(3);
        CalculateMetrics(initial_metric);
        mA_ab_covariant = initial_metric.a_ab_covariant;
        mB_ab_covariant = initial_metric.b_ab_covariant;

        CalculateTransformation(
            initial_metric,
            mT);

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::InitializeMaterial()
    {
        KRATOS_TRY

            if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
                mConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLaw->InitializeMaterial(GetProperties(),
                    GetGeometry(),
                    GetValue(SHAPE_FUNCTION_VALUES)
                );
            }
            else
                KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH("");
    }

    void IgaShell3pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        // definition of problem size
        const int number_of_nodes = NumberOfNodes();
        const int mat_size = NumberOfDofs();

        //reading in of integration weight, shape function values and shape function derivatives
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT) * mdetJ;

        // Compute actual metric variables
        MetricVariables actual_metric(3);
        CalculateMetrics(actual_metric);

        // Calculate Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_bending(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane,
            constitutive_variables_bending,
            Values, ConstitutiveLaw::StressMeasure_PK2);

        // Calculate B matrices
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        Matrix BCurvature = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);
        CalculateBBending(BCurvature, actual_metric);

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // Calculate second variations of strain and bending strain
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_bending(mat_size);
            CalculateSecondVariationStrainAndBending(
                second_variations_strain,
                second_variations_bending,
                actual_metric);

            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
            //adding curvature contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, constitutive_variables_bending.D, integration_weight);

            // adding  non-linear-contribution to Stiffness-Matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                second_variations_strain,
                constitutive_variables_membrane.S,
                integration_weight);

            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                second_variations_bending,
                constitutive_variables_bending.S,
                integration_weight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            noalias(rRightHandSideVector) -=
                integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            noalias(rRightHandSideVector) -=
                integration_weight * prod(trans(BCurvature), constitutive_variables_bending.S);
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void IgaShell3pElement::CalculateMetrics(
        MetricVariables& rMetric
    )
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(),
            DN_De,
            3, 2, rMetric.J);

        // curvilinear base vectors a1 and a2
        rMetric.a1[0] = rMetric.J(0, 0);
        rMetric.a2[0] = rMetric.J(0, 1);
        rMetric.a1[1] = rMetric.J(1, 0);
        rMetric.a2[1] = rMetric.J(1, 1);
        rMetric.a1[2] = rMetric.J(2, 0);
        rMetric.a2[2] = rMetric.J(2, 1);

        // basis vector a3
        MathUtils<double>::CrossProduct(rMetric.a3_tilde, rMetric.a1, rMetric.a2);
        // differential area dA
        rMetric.dA = norm_2(rMetric.a3_tilde);
        // a3 = normalized vector of a3_tilde
        rMetric.a3 = rMetric.a3_tilde / rMetric.dA;

        //GetCovariantMetric
        rMetric.a_ab_covariant[0] = inner_prod(rMetric.a1, rMetric.a1); // pow(rMetric.a1[0], 2) + pow(rMetric.a1[1], 2) + pow(rMetric.a1[2], 2);
        rMetric.a_ab_covariant[1] = inner_prod(rMetric.a2, rMetric.a2); // pow(rMetric.a2[0], 2) + pow(rMetric.a2[1], 2) + pow(rMetric.a2[2], 2);
        rMetric.a_ab_covariant[2] = inner_prod(rMetric.a1, rMetric.a2); // rMetric.a1[0] * rMetric.a2[0] + rMetric.a1[1] * rMetric.a2[1] + rMetric.a1[2] * rMetric.a2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        //rMetric.b_ab_covariant[0] = rMetric.H(0, 0)*rMetric.a3[0] + rMetric.H(1, 0)*rMetric.a3[1] + rMetric.H(2, 0)*rMetric.a3[2];
        //rMetric.b_ab_covariant[1] = rMetric.H(0, 1)*rMetric.a3[0] + rMetric.H(1, 1)*rMetric.a3[1] + rMetric.H(2, 1)*rMetric.a3[2];
        //rMetric.b_ab_covariant[2] = rMetric.H(0, 2)*rMetric.a3[0] + rMetric.H(1, 2)*rMetric.a3[1] + rMetric.H(2, 2)*rMetric.a3[2];
        rMetric.b_ab_covariant = prod(rMetric.H, rMetric.a3);
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to 
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void IgaShell3pElement::CalculateTransformation(
        const MetricVariables& rMetric,
        Matrix& T
    )
    {
        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rMetric.a_ab_covariant[0] * rMetric.a_ab_covariant[1] - rMetric.a_ab_covariant[2] * rMetric.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant = ZeroVector(3);
        a_ab_contravariant[0] = inv_det_g_ab * rMetric.a_ab_covariant[1];
        a_ab_contravariant[1] = inv_det_g_ab * rMetric.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rMetric.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_1 = rMetric.a1*a_ab_contravariant[0] + rMetric.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = rMetric.a1*a_ab_contravariant[2] + rMetric.a2*a_ab_contravariant[1];


        //Local cartesian coordinates
        double l_a1 = norm_2(rMetric.a1);
        array_1d<double, 3> e1 = rMetric.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        // e * a_contravariant
        Matrix G = ZeroMatrix(2, 2);
        G(0, 0) = inner_prod(e1, a_contravariant_1);
        G(0, 1) = inner_prod(e1, a_contravariant_2);
        G(1, 0) = inner_prod(e2, a_contravariant_1);
        G(1, 1) = inner_prod(e2, a_contravariant_2);

        //Transformation matrix T
        T = ZeroMatrix(3, 3);
        T(0, 0) = pow(G(0, 0), 2);
        T(0, 1) = pow(G(0, 1), 2);
        T(0, 2) = 2 * G(0, 0) * G(0, 1);

        T(1, 0) = pow(G(1, 0), 2);
        T(1, 1) = pow(G(1, 1), 2);
        T(1, 2) = 2 * G(1, 0) * G(1, 1);

        T(2, 0) = 2 * G(0, 0) * G(1, 0);
        T(2, 1) = 2 * G(0, 1) * G(1, 1);
        T(2, 2) = 2 * (G(0, 0) * G(1, 1) + G(0, 1) * G(1, 0));

        KRATOS_WATCH(T)
    }
    //************************************************************************************
    //************************************************************************************
    void IgaShell3pElement::CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rConstitutiveVariablesMembrane,
        ConstitutiveVariables& rConstitutiveVariablesBending,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //Computation of strain and bending strain in local cartesian coordinates
        array_1d<double, 3> membrane_strain_voigt = 0.5 * (rActualMetric.a_ab_covariant - mA_ab_covariant);
        array_1d<double, 3> curvature_change_voigt = mB_ab_covariant - rActualMetric.b_ab_covariant;

        rConstitutiveVariablesMembrane.E = prod(mT, membrane_strain_voigt);
        rConstitutiveVariablesBending.E = prod(mT, curvature_change_voigt);

        //Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rConstitutiveVariablesMembrane.S); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rConstitutiveVariablesMembrane.D); //this is an ouput parameter

        //rValues.CheckAllParameters();
        mConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);
        double thickness = this->GetProperties().GetValue(THICKNESS);

        //Integration through thickness
        // Membrane t * D
        // Bending  t^3/12 * D
        rConstitutiveVariablesMembrane.D *= thickness;
        rConstitutiveVariablesBending.D = rConstitutiveVariablesMembrane.D*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rConstitutiveVariablesMembrane.S = prod(
            trans(rConstitutiveVariablesMembrane.D), rConstitutiveVariablesMembrane.E);
        rConstitutiveVariablesBending.S = prod(
            trans(rConstitutiveVariablesBending.D), rConstitutiveVariablesBending.E);

        KRATOS_WATCH(rConstitutiveVariablesMembrane.D)
        KRATOS_WATCH(rConstitutiveVariablesBending.D)
    }

    //***********************************************************************************
    //***********************************************************************************

    void IgaShell3pElement::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        rB = ZeroMatrix(3, mat_size);

        for (int r = 0; r < static_cast<int>(mat_size); r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 3;
            int dirr = r % 3;

            array_1d<double, 3> dE_curvilinear = ZeroVector(3);
            // strain
            dE_curvilinear[0] = DN_De(kr, 0)*rMetric.a1(dirr);
            dE_curvilinear[1] = DN_De(kr, 1)*rMetric.a2(dirr);
            dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*rMetric.a2(dirr)
                                + rMetric.a1(dirr)*DN_De(kr, 1));

            rB(0, r) = mT(0, 0)*dE_curvilinear[0] + mT(0, 1)*dE_curvilinear[1] + mT(0, 2)*dE_curvilinear[2];
            rB(1, r) = mT(1, 0)*dE_curvilinear[0] + mT(1, 1)*dE_curvilinear[1] + mT(1, 2)*dE_curvilinear[2];
            rB(2, r) = mT(2, 0)*dE_curvilinear[0] + mT(2, 1)*dE_curvilinear[1] + mT(2, 2)*dE_curvilinear[2];
        }
    }

    //***********************************************************************************
    //***********************************************************************************
    void IgaShell3pElement::CalculateSecondVariationStrainMembrane(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        for (int r = 0; r<mat_size; r++)
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
                array_1d<double, 3> ddE_curvilinear = ZeroVector(3);
                if (dirr == dirs)
                {
                    ddE_curvilinear[0] = DN_De(kr, 0)*DN_De(ks, 0);
                    ddE_curvilinear[1] = DN_De(kr, 1)*DN_De(ks, 1);
                    ddE_curvilinear[2] = 0.5*(DN_De(kr, 0)*DN_De(ks, 1) + DN_De(kr, 1)*DN_De(ks, 0));
                }

                rSecondVariationsStrain.B11(r, s) = mT(0, 0)*ddE_curvilinear[0]
                    + mT(0, 1)*ddE_curvilinear[1]
                    + mT(0, 2)*ddE_curvilinear[2];
                rSecondVariationsStrain.B22(r, s) = mT(1, 0)*ddE_curvilinear[0]
                    + mT(1, 1)*ddE_curvilinear[1]
                    + mT(1, 2)*ddE_curvilinear[2];
                rSecondVariationsStrain.B12(r, s) = mT(2, 0)*ddE_curvilinear[0]
                    + mT(2, 1)*ddE_curvilinear[1]
                    + mT(2, 2)*ddE_curvilinear[2];
            }
        }
    }

    void IgaShell3pElement::CalculateBBending(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        KRATOS_TRY

        if (this->Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && this->Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
            const int number_of_control_points = GetGeometry().size();
            const int mat_size = number_of_control_points * 3;

            Matrix dg3 = ZeroMatrix(3, 3);
            Matrix dn = ZeroMatrix(3, 3);
            Matrix b = ZeroMatrix(3, mat_size);

            double invdA = 1 / rMetric.dA;
            double inddA3 = 1 / std::pow(rMetric.dA, 3);

            for (int i = 0; i < number_of_control_points; i++)
            {
                unsigned int index = 3 * i;
                //first line
                dg3(0, 0) = 0;
                dg3(0, 1) = -DN_De(i, 0) * rMetric.a2[2] + DN_De(i, 1)*rMetric.a1[2];
                dg3(0, 2) = DN_De(i, 0) * rMetric.a2[1] - DN_De(i, 1)*rMetric.a1[1];

                //second line
                dg3(1, 0) = DN_De(i, 0) * rMetric.a2[2] - DN_De(i, 1)*rMetric.a1[2];
                dg3(1, 1) = 0;
                dg3(1, 2) = -DN_De(i, 0)*rMetric.a2[0] + DN_De(i, 1)*rMetric.a1[0];

                //third line
                dg3(2, 0) = -DN_De(i, 0) * rMetric.a2[1] + DN_De(i, 1) * rMetric.a1[1];
                dg3(2, 1) = DN_De(i, 0) * rMetric.a2[0] - DN_De(i, 1) * rMetric.a1[0];
                dg3(2, 2) = 0;

                //KRATOS_WATCH(dg3)

                for (unsigned int j = 0; j < 3; j++)
                {
                    double g3dg3lg3 = (rMetric.a3_tilde[0] * dg3(j, 0)
                        + rMetric.a3_tilde[1] * dg3(j, 1)
                        + rMetric.a3_tilde[2] * dg3(j, 2))*inddA3;

                    dn(j, 0) = dg3(j, 0)*invdA - rMetric.a3_tilde[0] * g3dg3lg3;
                    dn(j, 1) = dg3(j, 1)*invdA - rMetric.a3_tilde[1] * g3dg3lg3;
                    dn(j, 2) = dg3(j, 2)*invdA - rMetric.a3_tilde[2] * g3dg3lg3;
                }

                // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
                b(0, index) = 0 - (DDN_DDe(i, 0) * rMetric.a3[0] + rMetric.H(0, 0)*dn(0, 0) + rMetric.H(1, 0)*dn(0, 1) + rMetric.H(2, 0)*dn(0, 2));
                b(0, index + 1) = 0 - (DDN_DDe(i, 0) * rMetric.a3[1] + rMetric.H(0, 0)*dn(1, 0) + rMetric.H(1, 0)*dn(1, 1) + rMetric.H(2, 0)*dn(1, 2));
                b(0, index + 2) = 0 - (DDN_DDe(i, 0) * rMetric.a3[2] + rMetric.H(0, 0)*dn(2, 0) + rMetric.H(1, 0)*dn(2, 1) + rMetric.H(2, 0)*dn(2, 2));

                //second line
                b(1, index) = 0 - (DDN_DDe(i, 1) * rMetric.a3[0] + rMetric.H(0, 1)*dn(0, 0) + rMetric.H(1, 1)*dn(0, 1) + rMetric.H(2, 1)*dn(0, 2));
                b(1, index + 1) = 0 - (DDN_DDe(i, 1) * rMetric.a3[1] + rMetric.H(0, 1)*dn(1, 0) + rMetric.H(1, 1)*dn(1, 1) + rMetric.H(2, 1)*dn(1, 2));
                b(1, index + 2) = 0 - (DDN_DDe(i, 1) * rMetric.a3[2] + rMetric.H(0, 1)*dn(2, 0) + rMetric.H(1, 1)*dn(2, 1) + rMetric.H(2, 1)*dn(2, 2));

                //third line
                b(2, index) = 0 - (DDN_DDe(i, 2) * rMetric.a3[0] + rMetric.H(0, 2)*dn(0, 0) + rMetric.H(1, 2)*dn(0, 1) + rMetric.H(2, 2)*dn(0, 2));
                b(2, index + 1) = 0 - (DDN_DDe(i, 2) * rMetric.a3[1] + rMetric.H(0, 2)*dn(1, 0) + rMetric.H(1, 2)*dn(1, 1) + rMetric.H(2, 2)*dn(1, 2));
                b(2, index + 2) = 0 - (DDN_DDe(i, 2) * rMetric.a3[2] + rMetric.H(0, 2)*dn(2, 0) + rMetric.H(1, 2)*dn(2, 1) + rMetric.H(2, 2)*dn(2, 2));
            }

            rB += prod(mT, b);
        }
        else
        {
            KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
        }

        KRATOS_CATCH("")
    }
    void IgaShell3pElement::CalculateSecondVariationStrainAndBending(
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric)
    {
        if (this->Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && this->Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

            const int number_of_control_points = GetGeometry().size();
            const int mat_size = number_of_control_points * 3;

            double lg3 = norm_2(rMetric.a3_tilde);
            double lg3_3 = pow(lg3, 3);
            double lg3_5 = pow(lg3, 5);
            double inv_lg3 = 1 / lg3;
            double inv_lg3_3 = 1 / lg3_3;
            double inv_lg3_5 = 1 / lg3_5;

            Matrix S_dg3 = ZeroMatrix(3, mat_size);
            Vector S_g3dg3 = ZeroVector(mat_size);
            Vector S_g3dg3lg3_3 = ZeroVector(mat_size);
            Matrix S_dn = ZeroMatrix(3, mat_size);
            // first variation of strain and curvature w.r.t. dof
            for (int r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                int kr = r / 3;
                int dirr = r % 3;

                array_1d<double, 3> S_dg_1 = ZeroVector(3);
                array_1d<double, 3> S_dg_2 = ZeroVector(3);
                S_dg_1(dirr) = DN_De(kr, 0);
                S_dg_2(dirr) = DN_De(kr, 1);

                // curvature
                S_dg3(0, r) = S_dg_1(1)*rMetric.a2(2) - S_dg_1(2)*rMetric.a2(1) + rMetric.a1(1)*S_dg_2(2) - rMetric.a1(2)*S_dg_2(1);
                S_dg3(1, r) = S_dg_1(2)*rMetric.a2(0) - S_dg_1(0)*rMetric.a2(2) + rMetric.a1(2)*S_dg_2(0) - rMetric.a1(0)*S_dg_2(2);
                S_dg3(2, r) = S_dg_1(0)*rMetric.a2(1) - S_dg_1(1)*rMetric.a2(0) + rMetric.a1(0)*S_dg_2(1) - rMetric.a1(1)*S_dg_2(0);

                S_g3dg3[r] = rMetric.a3_tilde[0] * S_dg3(0, r) + rMetric.a3_tilde[1] * S_dg3(1, r) + rMetric.a3_tilde[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.a3_tilde[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.a3_tilde[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.a3_tilde[2] * S_g3dg3lg3_3[r];
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
                        ddE_cu[0] = DN_De(kr, 0)*DN_De(ks, 0);
                        ddE_cu[1] = DN_De(kr, 1)*DN_De(ks, 1);
                        ddE_cu[2] = 0.5*(DN_De(kr, 0)*DN_De(ks, 1) + DN_De(kr, 1)*DN_De(ks, 0));

                        rSecondVariationsStrain.B11(r, s) = mT(0, 0)*ddE_cu[0] + mT(0, 1)*ddE_cu[1] + mT(0, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B22(r, s) = mT(1, 0)*ddE_cu[0] + mT(1, 1)*ddE_cu[1] + mT(1, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B12(r, s) = mT(2, 0)*ddE_cu[0] + mT(2, 1)*ddE_cu[1] + mT(2, 2)*ddE_cu[2];
                    }

                    // curvature
                    array_1d<double, 3> ddg3 = ZeroVector(3);
                    int dirt = 4 - dirr - dirs;
                    int ddir = dirr - dirs;
                    if (ddir == -1)      ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 2) ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 1) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == -2) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);

                    double c = -(ddg3[0] * rMetric.a3_tilde[0] + ddg3[1] * rMetric.a3_tilde[1] + ddg3[2] * rMetric.a3_tilde[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.a3_tilde[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.a3_tilde[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.a3_tilde[2];

                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2];
                    ddK_cu[1] = DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2];
                    ddK_cu[2] = DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2];

                    rSecondVariationsCurvature.B11(r, s) = mT(0, 0)*ddK_cu[0] + mT(0, 1)*ddK_cu[1] + mT(0, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B22(r, s) = mT(1, 0)*ddK_cu[0] + mT(1, 1)*ddK_cu[1] + mT(1, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B12(r, s) = mT(2, 0)*ddK_cu[0] + mT(2, 1)*ddK_cu[1] + mT(2, 2)*ddK_cu[2];
                }
            }
        }
    }


    void IgaShell3pElement::CalculateAndAddKm(
        Matrix& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight
    )
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& S,
        const double& rIntegrationWeight)

    {
        KRATOS_TRY
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
            {
                double nm = (S[0] * SecondVariationsStrain.B11(n, m)
                    + S[1] * SecondVariationsStrain.B22(n, m)
                    + S[2] * SecondVariationsStrain.B12(n, m))*rIntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if (n != m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        double density = this->GetProperties().GetValue(DENSITY);
        double mass = thickness * density * mdetJ * integration_weight;

        const int number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfNodes();

        if (rMassMatrix.size1() != mat_size)
            rMassMatrix.resize(mat_size, mat_size, false);
        rMassMatrix = ZeroMatrix(mat_size, mat_size);

        for (unsigned int r = 0; r<number_of_control_points; r++)
        {
            for (unsigned int s = 0; s<number_of_control_points; s++)
            {
                rMassMatrix(3 * s, 3 * r) = N(s)*N(r) * mass;
                rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
            }
        }
        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void IgaShell3pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 5;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void IgaShell3pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(NumberOfDofs());

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    //************************************************************************************
    //************************************************************************************
    int IgaShell3pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        if (DISPLACEMENT.Key() == 0)
            KRATOS_ERROR << "DISPLACEMENT has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_VALUES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_VALUES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        return 0;
    }


} // Namespace Kratos


