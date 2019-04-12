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
#include "utilities/math_utils.h"

// External includes

// Project includes
#include "custom_elements/iga_shell_5p_element.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell5pElement::Initialize()
    {
        KRATOS_TRY

        //Constitutive Law initialisation
        BaseDiscreteElement::Initialize();

        CalculateMetric(mInitialMetric);

        // KRATOS_WATCH(mInitialMetric.Q)

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
        
        // if (m_phi1 == 1)
        // KRATOS_WATCH("here: CalculateAllStart");

        
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

        //reading in of shape function derivatives
        const Matrix&  DN_De  = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
    
        Vector shear_difference_vector = ZeroVector(3);
        // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * g1 + w_alpha(2) * g2)
        Vector w_alpha = ZeroVector(2);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);

        MetricVariables actual_metric(3, 6);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(shear_difference_vector, w_alpha, Dw_alpha_Dbeta, actual_metric.g1, actual_metric.g2);
        ConstitutiveVariables constitutive_variables_membrane(6);
        ConstitutiveVariables constitutive_variables_curvature(6);
        ConstitutiveVariables constitutive_variables_membrane_hr(6);
        ConstitutiveVariables constitutive_variables_curvature_hr(6);
        CalculateConstitutiveVariables(actual_metric, shear_difference_vector,
            w_alpha, Dw_alpha_Dbeta, constitutive_variables_membrane, constitutive_variables_curvature,
            constitutive_variables_membrane_hr, constitutive_variables_curvature_hr, Values, 
            ConstitutiveLaw::StressMeasure_PK2);

        // calculate B MATRICES
        Matrix BMembrane = ZeroMatrix(6, mat_size);
        Matrix BCurvature = ZeroMatrix(6, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);
        CalculateBCurvature(BCurvature, actual_metric);

        double integration_weight = GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA;

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // KRATOS_WATCH(BMembrane)
            
            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrainCurvature(
                second_variations_strain,
                second_variations_curvature,
                actual_metric);
            
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
            
            // KRATOS_WATCH(rLeftHandSideMatrix)

            //adding curvature contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, constitutive_variables_curvature.D, integration_weight);

            // adding  non-linear-contribution to Stiffness-Matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                second_variations_strain,
                constitutive_variables_membrane.S,
                integration_weight);
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                second_variations_curvature,
                constitutive_variables_curvature.S,
                integration_weight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // KRATOS_WATCH(BMembrane)
            // KRATOS_WATCH(constitutive_variables_membrane.S)
            
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.S);

            // KRATOS_WATCH(rRightHandSideVector)
        }

        //KRATOS_WATCH(rLeftHandSideMatrix)

        KRATOS_CATCH("");
    }

    void IgaShell5pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight
    )
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 5;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
            {
                double nm = (SD[0] * SecondVariationsStrain.B11(n, m)
                    + SD[1] * SecondVariationsStrain.B22(n, m)
                    + SD[2] * SecondVariationsStrain.B33(n, m)
                    + SD[3] * SecondVariationsStrain.B12(n, m)
                    + SD[4] * SecondVariationsStrain.B23(n, m)
                    + SD[5] * SecondVariationsStrain.B13(n, m)) * rIntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if(n!=m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateMetric(
        MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(), DN_De, 3, 2, rMetric.J);

        rMetric.g1[0] = rMetric.J(0, 0);
        rMetric.g2[0] = rMetric.J(0, 1);
        rMetric.g1[1] = rMetric.J(1, 0);
        rMetric.g2[1] = rMetric.J(1, 1);
        rMetric.g1[2] = rMetric.J(2, 0);
        rMetric.g2[2] = rMetric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(rMetric.g3_notnorm, rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3_notnorm);
        //normalized basis vector g3
        rMetric.g3 = rMetric.g3_notnorm / rMetric.dA;

        //GetCovariantMetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = pow(rMetric.g3[0], 2) + pow(rMetric.g3[1], 2) + pow(rMetric.g3[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.g3[0] + rMetric.H(1, 0) * rMetric.g3[1] + rMetric.H(2, 0) * rMetric.g3[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.g3[0] + rMetric.H(1, 1) * rMetric.g3[1] + rMetric.H(2, 1) * rMetric.g3[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.g3[0] + rMetric.H(1, 2) * rMetric.g3[1] + rMetric.H(2, 2) * rMetric.g3[2];

        //contravariant rMetric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.gab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.gab_con[1] = invdetGab*rMetric.gab[0];

        array_1d<double, 3> g_con_1 = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];
        // g_con_3 = g3
        
        //local cartesian coordinates
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;
        // e3 = g3 = g_con_3

        // transformation matrix Q from curvilinear to local cartesian coordinate system
        Matrix mG = ZeroMatrix(3, 3);
        mG(0, 0) = inner_prod(g_con_1, e1);
        mG(0, 1) = inner_prod(g_con_1, e2);
        // if (Id() == 4){
        //     KRATOS_WATCH(g_con_1)
        //     KRATOS_WATCH(rMetric.g3)
        // }
        mG(0, 2) = inner_prod(g_con_1, rMetric.g3); // = 0, should be checked (ML)
        mG(1, 0) = inner_prod(g_con_2, e1); // = 0, should be checked -> simplification possible? (ML)
        mG(1, 1) = inner_prod(g_con_2, e2);
        mG(1, 2) = inner_prod(g_con_2, rMetric.g3); // = 0, should be checked (ML)
        mG(2, 0) = inner_prod(rMetric.g3, e1); // = 0, should be checked (ML)
        mG(2, 1) = inner_prod(rMetric.g3, e2); // = 0, should be checked (ML)
        mG(2, 2) = inner_prod(rMetric.g3, rMetric.g3); // = 1, should be checked (ML)
        // if (Id() == 4)
        //     KRATOS_WATCH(mG(0,2))
        for (unsigned int i = 0; i < 3; i++){
            // if (Id() == 4)
            //     KRATOS_WATCH(rMetric.Q(0,2))
            for (unsigned int j = 0; j < 3; j++)
                rMetric.Q(i, j) = pow(mG(i, j), 2); // entries 11 to 33
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
            for (unsigned int j = 0; j < 2; j++){
                rMetric.Q(i, j+3) = 2.00 * mG(i, j) * mG(i, j+1);   // entries 14 to 35
                rMetric.Q(j+3, i) = 2.00 * mG(j, i) * mG(j+1, i);   // entries 41 to 53
            }
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
            rMetric.Q(i, 5) = 2.00 * mG(i, 0) * mG(i, 2); // entries 16 to 36
            // if (Id() == 4)
                // // KRATOS_WATCH(rMetric.Q)
            rMetric.Q(5, i) = 2.00 * mG(0, i) * mG(2, i); // entries 61 to 63
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
            }
        for (unsigned int i = 0; i < 2; i++){
            for (unsigned int j = 0; j < 2; j++)
                rMetric.Q(i+3, j+3) = 2.00 * (mG(i, j) * mG(i+1, j+1) + mG(i, j+1) * mG(i+1, j)); // entries 44 to 55
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
            rMetric.Q(i+3, 5) = 2.00 * (mG(i, 0) * mG(i+1, 2) + mG(i, 2) * mG(i+1, 0)); // entries 46 and 56
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
            rMetric.Q(5, i+3) = 2.00 * (mG(0, i) * mG(2, i+1) + mG(0, i+1) * mG(2, i)); // entries 64 and 65
            // if (Id() == 4)
                // KRATOS_WATCH(rMetric.Q)
        }
        rMetric.Q(5, 5) = 2.00 * (mG(0, 0) * mG(2, 2) + mG(0, 2) * mG(2, 0));

        // if (Id() == 4){
        //     KRATOS_WATCH(GetValue(LOCAL_COORDINATES))
        //     KRATOS_WATCH(rMetric.Q)
        // }

        // faster transformation matrix taking into account that a lot of entries become zero(ML)
        double mG_00 = inner_prod(g_con_1, e1);
        double mG_01 = inner_prod(g_con_1, e2);
        double mG_11 = inner_prod(g_con_2, e2);
        
        rMetric.Q(0, 0) = pow(mG_00, 2);
        rMetric.Q(0, 1) = pow(mG_01, 2);
        rMetric.Q(0, 3) = 2.00 * mG_00 * mG_01;
        rMetric.Q(1, 1) = pow(mG_11, 2);
        rMetric.Q(2, 2) = 1;
        rMetric.Q(3, 1) = 2.00 * mG_01 * mG_11;
        rMetric.Q(3, 3) = 2.00 * mG_00 * mG_11;
        rMetric.Q(4, 4) = 2.00 * mG_11;
        rMetric.Q(5, 4) = 2.00 * mG_01;
        rMetric.Q(5, 5) = 2.00 * mG_00;

        // if (Id() == 4)
        //     KRATOS_WATCH(rMetric.Q)
    }
    
    void IgaShell5pElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rShearDifferenceVector,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveVariables& rThisConstitutiveVariablesHRMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesHRCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        Vector strain_vector = ZeroVector(6);
        Vector curvature_vector = ZeroVector(6);
        // hr = hierarchic rotations
        Vector hr_strain_vector = ZeroVector(6);
        // hr = hierarchic rotations
        Vector hr_curvature_vector = ZeroVector(6);
        
        CalculateStrain(strain_vector, rActualMetric.gab);
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.Q, strain_vector);
        CalculateCurvature(curvature_vector, rActualMetric.curvature);
        rThisConstitutiveVariablesCurvature.E = prod(mInitialMetric.Q, curvature_vector);
        CalculateStrainHR(hr_strain_vector, rShearDifferenceVector, rActualMetric.g1, rActualMetric.g2);
        rThisConstitutiveVariablesHRMembrane.E = prod(mInitialMetric.Q, hr_strain_vector);
        CalculateCurvatureHR(hr_curvature_vector, rw_alpha, rDw_alpha_Dbeta, rActualMetric);
        rThisConstitutiveVariablesHRCurvature.E = prod(mInitialMetric.Q, hr_curvature_vector);

        // 0ML
        //Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter

        //rValues.CheckAllParameters();
        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        double thickness = GetProperties().GetValue(THICKNESS);

        rThisConstitutiveVariablesMembrane.D *= thickness;
        rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D * (pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
            trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.S = prod(
            trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);
    }

    void IgaShell5pElement::CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab)
    {
        KRATOS_TRY

        rStrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.gab[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.gab[1]);
        rStrainVector[3] = 0.5 * (rgab[2] - mInitialMetric.gab[2]);
        // the other entries are (remain) zero (KL)

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature)
    {
        KRATOS_TRY

        rCurvatureVector[0] = (rCurvature[0] - mInitialMetric.curvature[0]);
        rCurvatureVector[1] = (rCurvature[1] - mInitialMetric.curvature[1]);
        rCurvatureVector[3] = (rCurvature[2] - mInitialMetric.curvature[2]);
        // the other entries are (remain) zero (KL)

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateShearDifferenceVector(
        Vector& rShearDifferenceVector,
        Vector& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const Vector& rg1,
        const Vector& rg2)
    {
        KRATOS_TRY; 
        
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int pos = GetGeometry()[0].GetDofPosition(ROTATION_X);
        double w_1, w_2;

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            w_1 = GetGeometry()[i].GetDof(ROTATION_X, pos).GetVariable();
            w_2 = GetGeometry()[i].GetDof(ROTATION_Y, pos + 1).GetVariable();
            rDw_alpha_Dbeta(0, 0) += DN_De(i, 0) * w_1;
            rDw_alpha_Dbeta(0, 1) += DN_De(i, 1) * w_1;
            rDw_alpha_Dbeta(1, 0) += DN_De(i, 0) * w_2;
            rDw_alpha_Dbeta(1, 1) += DN_De(i, 1) * w_2; 
            rw_alpha(0) += N[i] * w_1;
            rw_alpha(1) += N[i] * w_2;
        }

        rShearDifferenceVector = rw_alpha(0) * rg1 + rw_alpha(1) * rg2;

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateStrainHR(
        Vector& rHRStrainVector,
        const Vector& rShearDifferenceVector,
        const Vector& rg1,
        const Vector& rg2)
    {
        rHRStrainVector[4] = 0.5 * inner_prod(rShearDifferenceVector, rg2);
        rHRStrainVector[5] = 0.5 * inner_prod(rShearDifferenceVector, rg1);
        // the other entries are (remain) zero
    }

    void IgaShell5pElement::CalculateCurvatureHR(
        Vector& rHRCurvatureVector,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric)
    {
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        // derivatives of the shear difference vector
        Vector Dw_D1 = ZeroVector(3);
        Vector Dw_D2 = ZeroVector(3);
        Dw_D1 = rDw_alpha_Dbeta(0, 0) * rActualMetric.g1 + rDw_alpha_Dbeta(1, 0) * rActualMetric.g2;
        Dw_D2 = rDw_alpha_Dbeta(0, 1) * rActualMetric.g1 + rDw_alpha_Dbeta(1, 1) * rActualMetric.g2;
        for (unsigned int i = 0; i < 3; i++)
        {
            Dw_D1[i] += rw_alpha(1) * rActualMetric.H(i, 0) + rw_alpha(2) * rActualMetric.H(i, 2);
            Dw_D2[i] += rw_alpha(2) * rActualMetric.H(i, 2) + rw_alpha(2) * rActualMetric.H(i, 1);
        }

        rHRCurvatureVector[0] = inner_prod(Dw_D1, rActualMetric.g1);
        rHRCurvatureVector[1] = inner_prod(Dw_D2, rActualMetric.g2);
        rHRCurvatureVector[3] = 0.5 * (inner_prod(Dw_D1, rActualMetric.g2) + inner_prod(Dw_D2, rActualMetric.g1));
        // the other entries are (remain) zero
    }

    void IgaShell5pElement::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

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

            Vector dE_curvilinear = ZeroVector(3);
            // strain
            dE_curvilinear[0] = DN_De(kr, 0)*rMetric.g1(dirr);
            dE_curvilinear[1] = DN_De(kr, 1)*rMetric.g2(dirr);
            dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*rMetric.g2(dirr) + rMetric.g1(dirr)*DN_De(kr, 1));

            rB(0, r) = mInitialMetric.Q(0, 0)*dE_curvilinear[0] + mInitialMetric.Q(0, 1)*dE_curvilinear[1] + mInitialMetric.Q(0, 2)*dE_curvilinear[2];
            rB(1, r) = mInitialMetric.Q(1, 0)*dE_curvilinear[0] + mInitialMetric.Q(1, 1)*dE_curvilinear[1] + mInitialMetric.Q(1, 2)*dE_curvilinear[2];
            rB(2, r) = mInitialMetric.Q(2, 0)*dE_curvilinear[0] + mInitialMetric.Q(2, 1)*dE_curvilinear[1] + mInitialMetric.Q(2, 2)*dE_curvilinear[2];
        }
        // KRATOS_WATCH(rB)
    }

    void IgaShell5pElement::CalculateBCurvature(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        KRATOS_TRY

        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
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
                dg3(0, 1) = -DN_De(i, 0) * rMetric.g2[2] + DN_De(i, 1)*rMetric.g1[2];
                dg3(0, 2) = DN_De(i, 0) * rMetric.g2[1] - DN_De(i, 1)*rMetric.g1[1];

                //second line
                dg3(1, 0) = DN_De(i, 0) * rMetric.g2[2] - DN_De(i, 1)*rMetric.g1[2];
                dg3(1, 1) = 0;
                dg3(1, 2) = -DN_De(i, 0)*rMetric.g2[0] + DN_De(i, 1)*rMetric.g1[0];

                //third line
                dg3(2, 0) = -DN_De(i, 0) * rMetric.g2[1] + DN_De(i, 1) * rMetric.g1[1];
                dg3(2, 1) = DN_De(i, 0) * rMetric.g2[0] - DN_De(i, 1) * rMetric.g1[0];
                dg3(2, 2) = 0;

                //KRATOS_WATCH(dg3)

                for (unsigned int j = 0; j < 3; j++)
                {
                    double g3dg3lg3 = (rMetric.g3_notnorm[0] * dg3(j, 0) + rMetric.g3_notnorm[1] * dg3(j, 1) + rMetric.g3_notnorm[2] * dg3(j, 2))*inddA3;

                    dn(j, 0) = dg3(j, 0)*invdA - rMetric.g3_notnorm[0] * g3dg3lg3;
                    dn(j, 1) = dg3(j, 1)*invdA - rMetric.g3_notnorm[1] * g3dg3lg3;
                    dn(j, 2) = dg3(j, 2)*invdA - rMetric.g3_notnorm[2] * g3dg3lg3;
                }

                // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
                b(0, index) = 0 - (DDN_DDe(i, 0) * rMetric.g3[0] + rMetric.H(0, 0)*dn(0, 0) + rMetric.H(1, 0)*dn(0, 1) + rMetric.H(2, 0)*dn(0, 2));
                b(0, index + 1) = 0 - (DDN_DDe(i, 0) * rMetric.g3[1] + rMetric.H(0, 0)*dn(1, 0) + rMetric.H(1, 0)*dn(1, 1) + rMetric.H(2, 0)*dn(1, 2));
                b(0, index + 2) = 0 - (DDN_DDe(i, 0) * rMetric.g3[2] + rMetric.H(0, 0)*dn(2, 0) + rMetric.H(1, 0)*dn(2, 1) + rMetric.H(2, 0)*dn(2, 2));

                //second line
                b(1, index) = 0 - (DDN_DDe(i, 1) * rMetric.g3[0] + rMetric.H(0, 1)*dn(0, 0) + rMetric.H(1, 1)*dn(0, 1) + rMetric.H(2, 1)*dn(0, 2));
                b(1, index + 1) = 0 - (DDN_DDe(i, 1) * rMetric.g3[1] + rMetric.H(0, 1)*dn(1, 0) + rMetric.H(1, 1)*dn(1, 1) + rMetric.H(2, 1)*dn(1, 2));
                b(1, index + 2) = 0 - (DDN_DDe(i, 1) * rMetric.g3[2] + rMetric.H(0, 1)*dn(2, 0) + rMetric.H(1, 1)*dn(2, 1) + rMetric.H(2, 1)*dn(2, 2));

                //third line
                b(2, index) = 0 - (DDN_DDe(i, 2) * rMetric.g3[0] + rMetric.H(0, 2)*dn(0, 0) + rMetric.H(1, 2)*dn(0, 1) + rMetric.H(2, 2)*dn(0, 2));
                b(2, index + 1) = 0 - (DDN_DDe(i, 2) * rMetric.g3[1] + rMetric.H(0, 2)*dn(1, 0) + rMetric.H(1, 2)*dn(1, 1) + rMetric.H(2, 2)*dn(1, 2));
                b(2, index + 2) = 0 - (DDN_DDe(i, 2) * rMetric.g3[2] + rMetric.H(0, 2)*dn(2, 0) + rMetric.H(1, 2)*dn(2, 1) + rMetric.H(2, 2)*dn(2, 2));
            }

            rB = - prod(mInitialMetric.Q, b);
        }
        else
        {
            KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
        }
        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateSecondVariationStrainCurvature(
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric)
    {
        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

            const int number_of_control_points = GetGeometry().size();
            const int mat_size = number_of_control_points * 3;

            double lg3_3 = pow(rMetric.dA, 3);
            double lg3_5 = pow(rMetric.dA, 5);
            double inv_lg3 = 1 / rMetric.dA;
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
                S_dg3(0, r) = S_dg_1(1)*rMetric.g2(2) - S_dg_1(2)*rMetric.g2(1) + rMetric.g1(1)*S_dg_2(2) - rMetric.g1(2)*S_dg_2(1);
                S_dg3(1, r) = S_dg_1(2)*rMetric.g2(0) - S_dg_1(0)*rMetric.g2(2) + rMetric.g1(2)*S_dg_2(0) - rMetric.g1(0)*S_dg_2(2);
                S_dg3(2, r) = S_dg_1(0)*rMetric.g2(1) - S_dg_1(1)*rMetric.g2(0) + rMetric.g1(0)*S_dg_2(1) - rMetric.g1(1)*S_dg_2(0);

                S_g3dg3[r] = rMetric.g3_notnorm[0] * S_dg3(0, r) + rMetric.g3_notnorm[1] * S_dg3(1, r) + rMetric.g3_notnorm[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.g3_notnorm[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.g3_notnorm[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.g3_notnorm[2] * S_g3dg3lg3_3[r];
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

                        rSecondVariationsStrain.B11(r, s) = mInitialMetric.Q(0, 0)*ddE_cu[0] + mInitialMetric.Q(0, 1)*ddE_cu[1] + mInitialMetric.Q(0, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B22(r, s) = mInitialMetric.Q(1, 0)*ddE_cu[0] + mInitialMetric.Q(1, 1)*ddE_cu[1] + mInitialMetric.Q(1, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B12(r, s) = mInitialMetric.Q(2, 0)*ddE_cu[0] + mInitialMetric.Q(2, 1)*ddE_cu[1] + mInitialMetric.Q(2, 2)*ddE_cu[2];
                    }

                    // curvature
                    array_1d<double, 3> ddg3 = ZeroVector(3);
                    int dirt = 4 - dirr - dirs;
                    int ddir = dirr - dirs;
                    if (ddir == -1)      ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 2) ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 1) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == -2) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);

                    double c = -(ddg3[0] * rMetric.g3_notnorm[0] + ddg3[1] * rMetric.g3_notnorm[1] + ddg3[2] * rMetric.g3_notnorm[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.g3_notnorm[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.g3_notnorm[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.g3_notnorm[2];

                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2];
                    ddK_cu[1] = DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2];
                    ddK_cu[2] = DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2];

                    rSecondVariationsCurvature.B11(r, s) = mInitialMetric.Q(0, 0)*ddK_cu[0] + mInitialMetric.Q(0, 1)*ddK_cu[1] + mInitialMetric.Q(0, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B11(s, r) = rSecondVariationsCurvature.B11(r, s);
                    rSecondVariationsCurvature.B22(r, s) = mInitialMetric.Q(1, 0)*ddK_cu[0] + mInitialMetric.Q(1, 1)*ddK_cu[1] + mInitialMetric.Q(1, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B22(s, r) = rSecondVariationsCurvature.B22(r, s);
                    rSecondVariationsCurvature.B12(r, s) = mInitialMetric.Q(2, 0)*ddK_cu[0] + mInitialMetric.Q(2, 1)*ddK_cu[1] + mInitialMetric.Q(2, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B12(s, r) = rSecondVariationsCurvature.B12(r, s);
                }
            }
        }
    }

    void IgaShell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int number_of_control_points = GetGeometry().size();

        if (rResult.size() != 5 * number_of_control_points)
            rResult.resize(5 * number_of_control_points, false);

        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 5;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X, pos + 3).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y, pos + 4).EquationId();
        }

        KRATOS_CATCH("")
    };

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
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
        }

        KRATOS_CATCH("")
    };

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
