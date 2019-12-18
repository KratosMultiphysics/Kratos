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
#include "custom_elements/iga_shell_5p_element.h"


namespace Kratos
{
    void IgaShell5pElement::Initialize()
    {
        KRATOS_TRY

        IgaShell5pElement::Initialize();

        MetricVariables initial_metric(3);
        CalculateMetric(initial_metric);
        mInitialMetric = initial_metric;

        KRATOS_CATCH("")
    }


    void IgaShell5pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        // definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 5;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }

        //reading in of integration weight, shape function values and shape function derivatives
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        Vector   N     = this->GetValue(SHAPE_FUNCTION_VALUES);
        Matrix  DN_De  = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        Matrix DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        //KRATOS_WATCH(N)
        //KRATOS_WATCH(DN_De)
        //KRATOS_WATCH(DDN_DDe)


        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane, constitutive_variables_curvature,
            Values, ConstitutiveLaw::StressMeasure_PK2);

        // calculate B MATRICES
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        Matrix BCurvature = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);
        CalculateBCurvature(BCurvature, actual_metric);

        // Nonlinear Deformation
        SecondVariations second_variations_strain(mat_size);
        SecondVariations second_variations_curvature(mat_size);
        CalculateSecondVariationStrainCurvature(
            second_variations_strain,
            second_variations_curvature,
            actual_metric);

        integration_weight = this->GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA; // *GetProperties()[THICKNESS];

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
            //adding curvature contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, constitutive_variables_curvature.D, integration_weight);

            // adding  non-linear-contribution to Stiffness-Matrix
            //CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
            //    second_variations_strain,
            //    constitutive_variables_membrane.S,
            //    integration_weight);

            //CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
            //    second_variations_curvature,
            //    constitutive_variables_curvature.S,
            //    integration_weight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // operation performed: rRightHandSideVector -= Weight*IntForce
            //noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            //noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.S);
        }

        //KRATOS_WATCH(rLeftHandSideMatrix)

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void IgaShell5pElement::CalculateMetric(
        MetricVariables& metric
    )
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        Jacobian(DN_De, metric.J);


        metric.g1[0] = metric.J(0, 0);
        metric.g2[0] = metric.J(0, 1);
        metric.g1[1] = metric.J(1, 0);
        metric.g2[1] = metric.J(1, 1);
        metric.g1[2] = metric.J(2, 0);
        metric.g2[2] = metric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(metric.g3, metric.g1, metric.g2);
        //differential area dA
        metric.dA = norm_2(metric.g3);
        //normal vector _n
        Vector n = metric.g3 / metric.dA;


        //GetCovariantMetric
        metric.gab[0] = pow(metric.g1[0], 2) + pow(metric.g1[1], 2) + pow(metric.g1[2], 2);
        metric.gab[1] = pow(metric.g2[0], 2) + pow(metric.g2[1], 2) + pow(metric.g2[2], 2);
        metric.gab[2] = metric.g1[0] * metric.g2[0] + metric.g1[1] * metric.g2[1] + metric.g1[2] * metric.g2[2];

        IgaGeometryUtilitites::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            metric.H);

        metric.curvature[0] = metric.H(0, 0)*n[0] + metric.H(1, 0)*n[1] + metric.H(2, 0)*n[2];
        metric.curvature[1] = metric.H(0, 1)*n[0] + metric.H(1, 1)*n[1] + metric.H(2, 1)*n[2];
        metric.curvature[2] = metric.H(0, 2)*n[0] + metric.H(1, 2)*n[1] + metric.H(2, 2)*n[2];


        //contravariant metric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (metric.gab[0] * metric.gab[1] - metric.gab[2] * metric.gab[2]);
        metric.gab_con[0] = invdetGab*metric.gab[1];
        metric.gab_con[2] = -invdetGab*metric.gab[2];
        metric.gab_con[1] = invdetGab*metric.gab[0];


        array_1d<double, 3> g_con_1 = metric.g1*metric.gab_con[0] + metric.g2*metric.gab_con[2];
        array_1d<double, 3> g_con_2 = metric.g1*metric.gab_con[2] + metric.g2*metric.gab_con[1];


        //local cartesian coordinates
        double lg1 = norm_2(metric.g1);
        array_1d<double, 3> e1 = metric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g_con_1);
        mG(0, 1) = inner_prod(e1, g_con_2);
        mG(1, 0) = inner_prod(e2, g_con_1);
        mG(1, 1) = inner_prod(e2, g_con_2);

        metric.Q = ZeroMatrix(3, 3);
        metric.Q(0, 0) = pow(mG(0, 0), 2);
        metric.Q(0, 1) = pow(mG(0, 1), 2);
        metric.Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

        metric.Q(1, 0) = pow(mG(1, 0), 2);
        metric.Q(1, 1) = pow(mG(1, 1), 2);
        metric.Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

        metric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        metric.Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
        metric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));

        //Matrix T_G_E = ZeroMatrix(3, 3);
        //Transformation matrix T from contravariant to local cartesian basis
        double eG11 = inner_prod(e1, metric.g1);
        double eG12 = inner_prod(e1, metric.g2);
        double eG21 = inner_prod(e2, metric.g1);
        double eG22 = inner_prod(e2, metric.g2);

        //metric.Q = ZeroMatrix(3, 3);
        //metric.Q(0, 0) = eG11*eG11;
        //metric.Q(0, 1) = eG12*eG12;
        //metric.Q(0, 2) = 2.0*eG11*eG12;
        //metric.Q(1, 0) = eG21*eG21;
        //metric.Q(1, 1) = eG22*eG22;
        //metric.Q(1, 2) = 2.0*eG21*eG22;
        //metric.Q(2, 0) = 2.0*eG11*eG21;
        //metric.Q(2, 1) = 2.0*eG12*eG22;
        //metric.Q(2, 2) = 2.0*eG11*eG22 + eG12*eG21;

        metric.T = ZeroMatrix(3, 3);
        metric.T(0, 0) = eG11 * eG11;
        metric.T(0, 1) = eG12 * eG12;
        metric.T(0, 2) = 2.0*eG11*eG12;
        metric.T(1, 0) = eG21 * eG21;
        metric.T(1, 1) = eG22 * eG22;
        metric.T(1, 2) = 2.0*eG21*eG22;
        metric.T(2, 0) = 2.0*eG11*eG21;
        metric.T(2, 1) = 2.0*eG12*eG22;
        metric.T(2, 2) = 2.0*(eG11*eG22 + eG12 * eG21);

        //metric.T = ZeroMatrix(3, 3);
        //metric.T(0, 0) = eG11*eG11;
        //metric.T(0, 1) = eG21*eG21;
        //metric.T(0, 2) = 2.0*eG11*eG21;
        //metric.T(1, 0) = eG12*eG12;
        //metric.T(1, 1) = eG22*eG22;
        //metric.T(1, 2) = 2.0*eG12*eG22;
        //metric.T(2, 0) = eG11*eG12;
        //metric.T(2, 1) = eG21*eG22;
        //metric.T(2, 2) = eG11*eG22 + eG12*eG21;

        //KRATOS_WATCH(metric.Q)
        //KRATOS_WATCH(metric.T)
    }
    //************************************************************************************
    //************************************************************************************
    void IgaShell5pElement::CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        Vector strain_vector = ZeroVector(3);
        Vector curvature_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualMetric.gab, mInitialMetric.gab);
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.T, strain_vector);
        CalculateCurvature(curvature_vector, rActualMetric.curvature, mInitialMetric.curvature);
        rThisConstitutiveVariablesCurvature.E = prod(mInitialMetric.T, curvature_vector);

        //Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter

        //rValues.CheckAllParameters();
        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        double thickness = this->GetProperties().GetValue(THICKNESS);

        rThisConstitutiveVariablesMembrane.D *= thickness;
        rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
            trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.S = prod(
            trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void IgaShell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int number_of_control_points = GetGeometry().size();

        if (rResult.size() != 5 * number_of_control_points)
            rResult.resize(5 * number_of_control_points, false);

        //const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 5;
            rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        }

        KRATOS_CATCH("")
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void IgaShell5pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
        }

        KRATOS_CATCH("")
    };

    //************************************************************************************
    //************************************************************************************
    int IgaShell5pElement::Check(const ProcessInfo& rCurrentProcessInfo)
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


