//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Michael Loibl
//


// System includes
#include "utilities/math_utils.h"

// External includes

// Project includes
#include "custom_elements/iga_shell_5p_element.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"


namespace Kratos
{
    void IgaShell5pElement::Initialize()
    {
        KRATOS_TRY

        // KRATOS_WATCH("InitializeStart");

        // Constitutive Law initialisation
        BaseDiscreteElement::Initialize();

        MetricVariables test_metric(3);
        CalculateMetric(test_metric);
        CalculateMetric(m_initial_metric);
        m_Phi = ZeroVector(3);

        m_Dv_D1 = ZeroVector(3);
        m_Dv_D2 = ZeroVector(3);

        // KRATOS_WATCH("InitializeEnd");

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
        
        KRATOS_WATCH("here: CalculateAllStart");

        // definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 5;        // considers number of parameters

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);        // direct initialisation of the structure 'Values' (name) of the type Parameters

        /** decides where the strains, stresses and the constitutive tensor are calculated
         * the following settings determine: strains from element, stresses and constitutive tensor from constitutive law
         */
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
        double integration_weight = GetValue(INTEGRATION_WEIGHT);
        Vector   N     = GetValue(SHAPE_FUNCTION_VALUES);
        Matrix  DN_De  = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        Matrix DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        
        m_Dv_D1 = actual_metric.g1 - m_initial_metric.g1;
        m_Dv_D2 = actual_metric.g2 - m_initial_metric.g2;
        CalculateRotationVector(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane,
            constitutive_variables_curvature,
            Values, 
            ConstitutiveLaw::StressMeasure_PK2);

        // // calculate B MATRICES
        // Matrix BMembrane = ZeroMatrix(3, mat_size);
        // Matrix BCurvature = ZeroMatrix(3, mat_size);
        // CalculateBMembrane(BMembrane, actual_metric);
        // CalculateBCurvature(BCurvature, actual_metric);

        // // Nonlinear Deformation
        // SecondVariations second_variations_strain(mat_size);
        // SecondVariations second_variations_curvature(mat_size);
        // CalculateSecondVariationStrainCurvature(
        //     second_variations_strain,
        //     second_variations_curvature,
        //     actual_metric);

        // integration_weight = this->GetValue(INTEGRATION_WEIGHT) * m_initial_metric.dA; // *GetProperties()[THICKNESS];

        // // LEFT HAND SIDE MATRIX
        // if (CalculateStiffnessMatrixFlag == true)
        // {
        //     //adding membrane contributions to the stiffness matrix
        //     CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
        //     //adding curvature contributions to the stiffness matrix
        //     CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, constitutive_variables_curvature.D, integration_weight);

        //     // adding  non-linear-contribution to Stiffness-Matrix
        //     //CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
        //     //    second_variations_strain,
        //     //    constitutive_variables_membrane.S,
        //     //    integration_weight);

        //     //CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
        //     //    second_variations_curvature,
        //     //    constitutive_variables_curvature.S,
        //     //    integration_weight);
        // }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // operation performed: rRightHandSideVector -= Weight*IntForce
            //noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            //noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.S);
        }

        //KRATOS_WATCH(rLeftHandSideMatrix)
        
        KRATOS_WATCH("here: CalculateAllEnd");

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateMetric(
        MetricVariables& rMetric)
    {
        KRATOS_WATCH("CalculateMetricStart");

        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        
        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(), DN_De, 3, 3, rMetric.J);

        rMetric.g1[0] = rMetric.J(0, 0);
        rMetric.g2[0] = rMetric.J(0, 1);
        rMetric.g1[1] = rMetric.J(1, 0);
        rMetric.g2[1] = rMetric.J(1, 1);
        rMetric.g1[2] = rMetric.J(2, 0);
        rMetric.g2[2] = rMetric.J(2, 1);

        //basis vector g3
        rMetric.g3 = MathUtils<double>::CrossProduct(rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3);
        //normal vector _n
        Vector n = rMetric.g3 / rMetric.dA;

        //GetCovariantmetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        // not checked yet (ML)
        rMetric.curvature[0] = rMetric.H(0, 0)*n[0] + rMetric.H(1, 0)*n[1] + rMetric.H(2, 0)*n[2];
        rMetric.curvature[1] = rMetric.H(0, 1)*n[0] + rMetric.H(1, 1)*n[1] + rMetric.H(2, 1)*n[2];
        rMetric.curvature[2] = rMetric.H(0, 2)*n[0] + rMetric.H(1, 2)*n[1] + rMetric.H(2, 2)*n[2];

        //contravariant metric gab_con
        double inv_det_gab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = inv_det_gab*rMetric.gab[1];
        rMetric.gab_con[1] = inv_det_gab*rMetric.gab[0];
        rMetric.gab_con[2] = -inv_det_gab*rMetric.gab[2];

        // contravariant base vectors g_con
        array_1d<double, 3> g_con_1 = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];

        //local cartesian base vectors, with the following choice its guaranteed that e1 is orthogonal to e2
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        // calculation of Q not checked yet (ML)
        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g_con_1);
        mG(0, 1) = inner_prod(e1, g_con_2);
        mG(1, 0) = inner_prod(e2, g_con_1);
        mG(1, 1) = inner_prod(e2, g_con_2);

        rMetric.Q(0, 0) = pow(mG(0, 0), 2);
        rMetric.Q(0, 1) = pow(mG(0, 1), 2);
        rMetric.Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

        rMetric.Q(1, 0) = pow(mG(1, 0), 2);
        rMetric.Q(1, 1) = pow(mG(1, 1), 2);
        rMetric.Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

        rMetric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        rMetric.Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
        rMetric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));

        //Matrix T_G_E = ZeroMatrix(3, 3);
        //Transformation matrix T from contravariant to local cartesian basis
        double eG11 = inner_prod(e1, rMetric.g1);
        double eG12 = inner_prod(e1, rMetric.g2);
        double eG21 = inner_prod(e2, rMetric.g1);
        double eG22 = inner_prod(e2, rMetric.g2);

        // does this commented part mean that the above calculation is not rigth for sure? (ML)
        //rMetric.Q = ZeroMatrix(3, 3);
        //rMetric.Q(0, 0) = eG11*eG11;
        //rMetric.Q(0, 1) = eG12*eG12;
        //rMetric.Q(0, 2) = 2.0*eG11*eG12;
        //rMetric.Q(1, 0) = eG21*eG21;
        //rMetric.Q(1, 1) = eG22*eG22;
        //rMetric.Q(1, 2) = 2.0*eG21*eG22;
        //rMetric.Q(2, 0) = 2.0*eG11*eG21;
        //rMetric.Q(2, 1) = 2.0*eG12*eG22;
        //rMetric.Q(2, 2) = 2.0*eG11*eG22 + eG12*eG21;

        // calculation of T not checked yet (ML)
        rMetric.T(0, 0) = eG11 * eG11;
        rMetric.T(0, 1) = eG12 * eG12;
        rMetric.T(0, 2) = 2.0*eG11*eG12;
        rMetric.T(1, 0) = eG21 * eG21;
        rMetric.T(1, 1) = eG22 * eG22;
        rMetric.T(1, 2) = 2.0*eG21*eG22;
        rMetric.T(2, 0) = 2.0*eG11*eG21;
        rMetric.T(2, 1) = 2.0*eG12*eG22;
        rMetric.T(2, 2) = 2.0*(eG11*eG22 + eG12 * eG21);

        // does this commented part mean that the above calculation is rigth for sure? (ML)
        //rMetric.T = ZeroMatrix(3, 3);
        //rMetric.T(0, 0) = eG11*eG11;
        //rMetric.T(0, 1) = eG21*eG21;
        //rMetric.T(0, 2) = 2.0*eG11*eG21;
        //rMetric.T(1, 0) = eG12*eG12;
        //rMetric.T(1, 1) = eG22*eG22;
        //rMetric.T(1, 2) = 2.0*eG12*eG22;
        //rMetric.T(2, 0) = eG11*eG12;
        //rMetric.T(2, 1) = eG21*eG22;
        //rMetric.T(2, 2) = eG11*eG22 + eG12*eG21;

        //KRATOS_WATCH(rMetric.Q)
        //KRATOS_WATCH(rMetric.T)
        KRATOS_WATCH("CalculateMetricEnd");
    }
 
    void IgaShell5pElement::CalculateRotationVector(
        MetricVariables& rActualMetric)
    {   
        m_phi1 = MathUtils<double>::Dot3(m_Dv_D2, m_initial_metric.g3)
        / norm_2(MathUtils<double>::CrossProduct(m_initial_metric.g1,m_initial_metric.g2));
        m_phi2 = MathUtils<double>::Dot3(m_Dv_D1, m_initial_metric.g3)
        / norm_2(MathUtils<double>::CrossProduct(m_initial_metric.g1,m_initial_metric.g2));
        m_Phi = m_phi1 * m_initial_metric.g1 + m_phi2 * m_initial_metric.g2; 
    }
    
    void IgaShell5pElement::CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        // membrane strain
        Vector strain_vector = ZeroVector(3);
        // strain due to bending/curvature
        Vector curvature_vector = ZeroVector(3);

        CalculateStrain(strain_vector);
        
        /** compare in the end with these formulas for the strain vector which are used in the existing KL
         * in a first small check they seemed slightly different (ML)
        Vector StrainVector = ZeroVector(3);
        StrainVector[0] = 0.5 * (rActualMetric.gab[0] - m_initial_metric.gab[0]);
        StrainVector[1] = 0.5 * (rActualMetric.gab[1] - m_initial_metric.gab[1]);
        StrainVector[2] = 0.5 * (rActualMetric.gab[2] - m_initial_metric.gab[2]);
        KRATOS_WATCH(StrainVector);
        */

        // membrane strain w.r.t. local cartesian coordinates through transformation
        rThisConstitutiveVariablesMembrane.E = prod(m_initial_metric.T, strain_vector);
        CalculateCurvature(curvature_vector, rActualMetric.curvature, m_initial_metric.curvature);
        
        Vector CurvatureVector = ZeroVector(3);
        CurvatureVector[0] = (rActualMetric.curvature[0] - m_initial_metric.curvature[0]);
        CurvatureVector[1] = (rActualMetric.curvature[1] - m_initial_metric.curvature[1]);
        CurvatureVector[2] = (rActualMetric.curvature[2] - m_initial_metric.curvature[2]);
        
        // strain due to bending w.r.t. local cartesian coordinates through transformation
        // rThisConstitutiveVariablesCurvature.E = prod(m_initial_metric.T, curvature_vector);

        // //Constitive Matrices DMembrane and DCurvature
        // rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        // rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        // rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter

        // //rValues.CheckAllParameters();
        // mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        // double thickness = this->GetProperties().GetValue(THICKNESS);

        // rThisConstitutiveVariablesMembrane.D *= thickness;
        // rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D*(pow(thickness, 2) / 12);

        // //Local Cartesian Forces and Moments
        // rThisConstitutiveVariablesMembrane.S = prod(
        //     trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        // rThisConstitutiveVariablesCurvature.S = prod(
        //     trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);
    }

    void IgaShell5pElement::CalculateStrain(
        Vector& rStrainVector)
    {
        KRATOS_TRY

        rStrainVector[0] = MathUtils<double>::Dot3(m_Dv_D1, m_initial_metric.g1);
        rStrainVector[1] = MathUtils<double>::Dot3(m_Dv_D2, m_initial_metric.g2);
        rStrainVector[2] = MathUtils<double>::Dot3(m_Dv_D1, m_initial_metric.g2) + MathUtils<double>::Dot3(m_Dv_D2, m_initial_metric.g1);

        KRATOS_WATCH(rStrainVector);

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateCurvature(
        Vector& CurvatureVector,
        Vector& bv,
        Vector& bv_ref)
    {
        KRATOS_TRY



        KRATOS_CATCH("")
    }

    void IgaShell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        KRATOS_WATCH("EquationIdVectorStart");

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
    }

    // /***********************************************************************************/
    // /***********************************************************************************/
    void IgaShell5pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_WATCH("GetDofListStart");

        const unsigned int number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) 
        {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            KRATOS_WATCH("GetDofList1");
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            KRATOS_WATCH("GetDofList2");
        }

        KRATOS_WATCH("GetDofListEnd");

        KRATOS_CATCH("")
    }

    // //************************************************************************************
    // //************************************************************************************
    int IgaShell5pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        
        KRATOS_WATCH("CheckStart");
        
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


