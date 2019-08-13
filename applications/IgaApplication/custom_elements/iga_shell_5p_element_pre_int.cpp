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
#include "custom_elements/iga_shell_5p_element_pre_int.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell5pElementPreInt::Initialize()
    {
        KRATOS_TRY

        // KRATOS_WATCH("start: Initialize")
        //Constitutive Law initialisation
        BaseDiscreteElement::Initialize();
        // Check whether ConstitutiveLaw is 3D
        if (mConstitutiveLawVector[0]->GetStrainSize() != 6){
            KRATOS_WATCH("ConstitutiveLaw is not 3D.")
            KRATOS_ERROR << "ConstitutiveLaw is not 3D." << std::endl;
        }

        CalculateMetric(mInitialMetric);

        // KRATOS_WATCH(mInitialMetric.Q)

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        
        // KRATOS_WATCH("start: CalculateAll")
        
        // definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
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

        Vector shear_difference_vector = ZeroVector(3);
        // derivatives of the shear difference vector
        Vector Dw_D1 = ZeroVector(3);
        Vector Dw_D2 = ZeroVector(3);
        // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * g1 + w_alpha(2) * g2)
        Vector w_alpha = ZeroVector(2);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);


        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(shear_difference_vector, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);

        ConstitutiveVariables constitutive_variables_membrane(5);
        ConstitutiveVariables constitutive_variables_curvature(5);
        CalculateConstitutiveVariables(actual_metric, shear_difference_vector, w_alpha, Dw_alpha_Dbeta, 
        Dw_D1, Dw_D2, constitutive_variables_membrane, constitutive_variables_curvature, Values, 
        ConstitutiveLaw::StressMeasure_PK2);

        // calculate B MATRICES
        Matrix B_membrane_KL = ZeroMatrix(5, mat_size);
        Matrix B_curvature_KL = ZeroMatrix(5, mat_size);
        CalculateBMembrane(B_membrane_KL, actual_metric);
        CalculateBCurvature(B_curvature_KL, actual_metric);
        Matrix B_membrane_RM = ZeroMatrix(5, mat_size);
        Matrix B_curvature_RM = ZeroMatrix(5, mat_size);            
        SecondVariations second_variations_membrane_RM(mat_size);
        SecondVariations second_variations_curvature_RM(mat_size);
        CalculateVariationsRM(B_membrane_RM, B_curvature_RM, second_variations_membrane_RM, second_variations_curvature_RM, 
            shear_difference_vector, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric, CalculateStiffnessMatrixFlag);
        Matrix B_membrane = ZeroMatrix(5, mat_size);
        B_membrane = B_membrane_KL + B_membrane_RM;
        Matrix B_curvature = ZeroMatrix(5, mat_size);
        B_curvature = B_curvature_KL + B_curvature_RM;
        
        double integration_weight = GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA;

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // Nonlinear Deformation
            SecondVariations second_variations_membrane_KL(mat_size);
            SecondVariations second_variations_curvature_KL(mat_size);        

            CalculateSecondVariationStrainCurvature(
                second_variations_membrane_KL,
                second_variations_curvature_KL,
                actual_metric);
            
            Matrix LeftHandSideMatrixTest = ZeroMatrix(mat_size, mat_size);
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, B_membrane, constitutive_variables_membrane.D, integration_weight);

            //adding curvature contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, B_curvature, constitutive_variables_curvature.D, integration_weight);
            LeftHandSideMatrixTest = rLeftHandSideMatrix - LeftHandSideMatrixTest;

            // adding  non-linear-contribution to Stiffness-Matrix
            SecondVariations second_variations_membrane(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            second_variations_membrane = second_variations_membrane_KL + second_variations_membrane_RM;
            second_variations_curvature = second_variations_curvature_KL + second_variations_curvature_RM;
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix, second_variations_membrane, constitutive_variables_membrane.S,
                integration_weight);
            LeftHandSideMatrixTest = rLeftHandSideMatrix - LeftHandSideMatrixTest;
            // KRATOS_WATCH(LeftHandSideMatrixTest)
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix, second_variations_curvature, constitutive_variables_curvature.S,
                integration_weight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(B_membrane), constitutive_variables_membrane.S);
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(B_curvature), constitutive_variables_curvature.S);
            if (Id() == 4)
                KRATOS_WATCH(rRightHandSideVector)
        }

        KRATOS_CATCH("");
    }

    void IgaShell5pElementPreInt::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double& rIntegrationWeight
    )
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += rIntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY

        const int number_of_nodes = GetGeometry().size();
        const int mat_size = number_of_nodes * 5;

        if (SD.size() != 5)
            KRATOS_ERROR << "Stress size is wrong." << std::endl;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
            {
                double nm = (SD[0] * SecondVariationsStrain.B11(n, m)
                    + SD[1] * SecondVariationsStrain.B22(n, m)
                    + SD[2] * SecondVariationsStrain.B12(n, m)
                    + SD[3] * SecondVariationsStrain.B23(n, m)
                    + SD[4] * SecondVariationsStrain.B13(n, m)) * rIntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if(n!=m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateMetric(
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
        MathUtils<double>::CrossProduct(rMetric.g3_tilde, rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3_tilde);
        //normalized basis vector g3
        rMetric.g3 = rMetric.g3_tilde / rMetric.dA;

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

        // // transformation matrix Q from curvilinear to local cartesian coordinate system
        // faster computation of transformation matrix Q taking into account that a lot of entries become zero (ML)
        // this matrix Q is referring to a VoigtSize 5 with E11, E22, E12, E23, E13
        double mG_00 = inner_prod(e1, g_con_1);
        double mG_10 = inner_prod(e2, g_con_1);
        double mG_11 = inner_prod(e2, g_con_2);
        
        rMetric.Q(0, 0) = pow(mG_00, 2);
        rMetric.Q(1, 0) = pow(mG_10, 2);
        rMetric.Q(1, 1) = pow(mG_11, 2);
        rMetric.Q(1, 2) = 2.00 * mG_10 * mG_11;
        rMetric.Q(2, 0) = 2.00 * mG_00 * mG_10;
        rMetric.Q(2, 2) = 2.00 * mG_00 * mG_11;
        rMetric.Q(3, 3) = 2.00 * mG_11;
        rMetric.Q(3, 4) = 2.00 * mG_10;
        rMetric.Q(4, 4) = 2.00 * mG_00;

    }

    void IgaShell5pElementPreInt::CalculateShearDifferenceVector(
        Vector& rShearDifferenceVector,
        Vector& rDw_D1,
        Vector& rDw_D2,
        Vector& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric)
    {
        KRATOS_TRY; 
        
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int pos = GetGeometry()[0].GetDofPosition(ROTATION_X);
        double w_1, w_2;

        for (unsigned int i = 0; i < number_of_nodes; ++i) 
        {
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            w_1 = GetGeometry()[i].GetDof(ROTATION_X, pos).GetSolutionStepValue();
            w_2 = GetGeometry()[i].GetDof(ROTATION_Y, pos + 1).GetSolutionStepValue();

            rDw_alpha_Dbeta(0, 0) += DN_De(i, 0) * w_1;
            rDw_alpha_Dbeta(0, 1) += DN_De(i, 1) * w_1;
            rDw_alpha_Dbeta(1, 0) += DN_De(i, 0) * w_2;
            rDw_alpha_Dbeta(1, 1) += DN_De(i, 1) * w_2; 
            rw_alpha(0) += N[i] * w_1;
            rw_alpha(1) += N[i] * w_2;
        }

        // derivatives of the shear difference vector
        rDw_D1 = rDw_alpha_Dbeta(0, 0) * rActualMetric.g1 + rDw_alpha_Dbeta(1, 0) * rActualMetric.g2;
        rDw_D2 = rDw_alpha_Dbeta(0, 1) * rActualMetric.g1 + rDw_alpha_Dbeta(1, 1) * rActualMetric.g2;

        for (unsigned int i = 0; i < 3; i++)
        {
            rDw_D1[i] += rw_alpha(0) * rActualMetric.H(i, 0) + rw_alpha(1) * rActualMetric.H(i, 2);
            rDw_D2[i] += rw_alpha(0) * rActualMetric.H(i, 2) + rw_alpha(1) * rActualMetric.H(i, 1);
        }

        rShearDifferenceVector = rw_alpha(0) * rActualMetric.g1 + rw_alpha(1) * rActualMetric.g2;

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rShearDifferenceVector,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        KRATOS_TRY

        Vector strain_vector = ZeroVector(5);
        Vector curvature_vector = ZeroVector(5);
        Vector strain_vector_RM = ZeroVector(5);
        Vector curvature_vector_RM = ZeroVector(5);
        
        // Strain computation in curvilinear space
        CalculateStrain(strain_vector, rActualMetric.gab);
        CalculateCurvature(curvature_vector, rActualMetric.curvature);
        CalculateStrainRM(strain_vector_RM, rShearDifferenceVector, rActualMetric.g1, rActualMetric.g2);
        CalculateCurvatureRM(curvature_vector_RM, rDw_D1, rDw_D2, rActualMetric.g1, rActualMetric.g2);
        rThisConstitutiveVariablesMembrane.E = strain_vector + strain_vector_RM;
        if (Id() == 4){
            KRATOS_WATCH(strain_vector)
            KRATOS_WATCH(strain_vector_RM)
        }
        rThisConstitutiveVariablesCurvature.E = curvature_vector + curvature_vector_RM;
        // Strain transformation to local Cartesian Space with VoigtSize 6 because ConstitutiveLaw is 3D
        ConstitutiveVariables constitutive_variables_membrane(6);
        TransformationCurvilinearStrainSize5ToCartesianStrainSize6(rThisConstitutiveVariablesMembrane.E, constitutive_variables_membrane.E);

        //Constitutive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(constitutive_variables_membrane.E); //this is the input parameter
        rValues.SetStressVector(constitutive_variables_membrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(constitutive_variables_membrane.D); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        // static condensation of  sigma_33
        unsigned int index_i = 0;
        for (unsigned int i = 0; i < 6; i++)
        {
            if (i != 2){
                unsigned int index_j = 0;
                for (unsigned int j = 0; j < 6; j++ ){
                    if (j != 2){
                        rThisConstitutiveVariablesMembrane.D(index_i, index_j) += constitutive_variables_membrane.D(i, j) - 
                            constitutive_variables_membrane.D(i, 2) * constitutive_variables_membrane.D(2, j) 
                            / constitutive_variables_membrane.D(2, 2);
                        index_j++;
                    }
                }
                index_i++;
            }
        }
        
        double thickness = GetProperties().GetValue(THICKNESS);

        // pre-integration
        rThisConstitutiveVariablesMembrane.D *= thickness;
        rThisConstitutiveVariablesMembrane.D(3, 3) = rThisConstitutiveVariablesMembrane.D(3, 3) * 5 / 6;       // ML
        rThisConstitutiveVariablesMembrane.D(3, 4) = rThisConstitutiveVariablesMembrane.D(3, 4) * 5 / 6;
        rThisConstitutiveVariablesMembrane.D(4, 3) = rThisConstitutiveVariablesMembrane.D(3, 4);
        rThisConstitutiveVariablesMembrane.D(4, 4) = rThisConstitutiveVariablesMembrane.D(4, 4) * 5 / 6;
        rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D * (pow(thickness, 2) / 12);

        // Strain Transformation to local Cartesian space with VoigtSize 5
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.Q, rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.E = prod(mInitialMetric.Q, rThisConstitutiveVariablesCurvature.E);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
            trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.S = prod(
            trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab)
    {
        KRATOS_TRY

        rStrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.gab[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.gab[1]);
        rStrainVector[2] = 0.5 * (rgab[2] - mInitialMetric.gab[2]);
        // the other entries are (remain) zero (KL)

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature)
    {
        KRATOS_TRY

        rCurvatureVector[0] = (mInitialMetric.curvature[0] - rCurvature[0]);
        rCurvatureVector[1] = (mInitialMetric.curvature[1] - rCurvature[1]);
        rCurvatureVector[2] = (mInitialMetric.curvature[2] - rCurvature[2]);
        // the other entries are (remain) zero (KL)

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateStrainRM(
        Vector& rStrainVectorRM,
        const Vector& rShearDifferenceVector,
        const Vector& rg1,
        const Vector& rg2)
    {
        rStrainVectorRM[3] = 0.5 * inner_prod(rShearDifferenceVector, rg2);
        rStrainVectorRM[4] = 0.5 * inner_prod(rShearDifferenceVector, rg1);
        // the other entries are (remain) zero
    }

    void IgaShell5pElementPreInt::CalculateCurvatureRM(
        Vector& rCurvatureVectorRM,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rg1,
        const Vector& rg2)
    {
        rCurvatureVectorRM[0] = inner_prod(rDw_D1, rg1);
        rCurvatureVectorRM[1] = inner_prod(rDw_D2, rg2);
        rCurvatureVectorRM[2] = 0.5 * (inner_prod(rDw_D1, rg2) + inner_prod(rDw_D2, rg1));
        // the other entries are (remain) zero
    }

    void IgaShell5pElementPreInt::TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector rCurvilinearStrain,
        Vector& rCartesianStrain)
    {
        KRATOS_TRY

        if (rCurvilinearStrain.size() != 5 || rCartesianStrain.size() != 6) 
            KRATOS_ERROR << "Wrong strain size in transformation." << std::endl;
        if (mInitialMetric.Q.size1() != 5 || mInitialMetric.Q.size2() != 5)
            KRATOS_ERROR << "Wrong size of transformation matrix Q." << std::endl;

        // transformation with simplified matrix
        rCartesianStrain[0] = mInitialMetric.Q(0, 0) * rCurvilinearStrain[0];
        rCartesianStrain[1] = mInitialMetric.Q(1, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(1, 1) * rCurvilinearStrain[1] 
            + mInitialMetric.Q(1, 2) * rCurvilinearStrain[2];
        rCartesianStrain[2] = 0.0; // RM
        rCartesianStrain[3] = mInitialMetric.Q(2, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(2, 2) * rCurvilinearStrain[2];
        rCartesianStrain[4] = mInitialMetric.Q(3, 3) * rCurvilinearStrain[3] + mInitialMetric.Q(3, 4) * rCurvilinearStrain[4];
        rCartesianStrain[5] = mInitialMetric.Q(4, 4) * rCurvilinearStrain[4];

        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;

        if (rB.size1() != 5 || rB.size2() != mat_size)
            rB.resize(5, mat_size);
        rB = ZeroMatrix(5, mat_size);

        for (int r = 0; r < static_cast<int>(mat_size); r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 5;
            int dirr = r % 5;

            // "if" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            if (dirr == 0 || dirr == 1 || dirr == 2)
            {
                Vector dE_curvilinear = ZeroVector(3);
                // strain corresponding to E11, E22, E12
                dE_curvilinear[0] = DN_De(kr, 0)*rMetric.g1(dirr);
                dE_curvilinear[1] = DN_De(kr, 1)*rMetric.g2(dirr);
                dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*rMetric.g2(dirr) + rMetric.g1(dirr)*DN_De(kr, 1));

                // calculated with simplified Q (ML)
                rB(0, r) = mInitialMetric.Q(0, 0) * dE_curvilinear[0];
                rB(1, r) = mInitialMetric.Q(1, 0) * dE_curvilinear[0] + mInitialMetric.Q(1, 1) * dE_curvilinear[1] 
                    + mInitialMetric.Q(1, 2) * dE_curvilinear[2];
                rB(2, r) = mInitialMetric.Q(2, 0) * dE_curvilinear[0] + mInitialMetric.Q(2, 2) * dE_curvilinear[2];
                // all other entries of rB are (remain) zero
            }
        }
    }

    void IgaShell5pElementPreInt::CalculateBCurvature(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        KRATOS_TRY

        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
            const unsigned int number_of_nodes = GetGeometry().size();
            const unsigned int mat_size_KL = number_of_nodes * 3;
            const unsigned int mat_size = number_of_nodes * 5;
            
            if (rB.size1() != 5 || rB.size2() != mat_size)
                rB.resize(5, mat_size);
            rB = ZeroMatrix(5, mat_size);

            Matrix dg3 = ZeroMatrix(3, 3);
            Matrix dn = ZeroMatrix(3, 3);
            Matrix b = ZeroMatrix(3, mat_size_KL);

            double invdA = 1 / rMetric.dA;
            double inddA3 = 1 / std::pow(rMetric.dA, 3);

            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                unsigned int index_KL = 3 * i;
                unsigned int index = 5 * i;
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

                for (unsigned int j = 0; j < 3; j++)
                {
                    double g3dg3lg3 = (rMetric.g3_tilde[0] * dg3(j, 0) + rMetric.g3_tilde[1] * dg3(j, 1) + rMetric.g3_tilde[2] * dg3(j, 2))*inddA3;

                    dn(j, 0) = dg3(j, 0)*invdA - rMetric.g3_tilde[0] * g3dg3lg3;
                    dn(j, 1) = dg3(j, 1)*invdA - rMetric.g3_tilde[1] * g3dg3lg3;
                    dn(j, 2) = dg3(j, 2)*invdA - rMetric.g3_tilde[2] * g3dg3lg3;
                }

                // b refers to curvilinear and rB to local cartesian coordinate system
                // "index" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
                for (unsigned int j = 0; j < 3; j++)
                {
                    for (unsigned int k = 0; k < 3; k++)
                    {
                        b(k, index_KL + j) = - (DDN_DDe(i, k) * rMetric.g3[j] + rMetric.H(0, k) * dn(j, 0) 
                        + rMetric.H(1, k) * dn(j, 1) + rMetric.H(2, k) * dn(j, 2));
                    }
                    rB(0, index + j) = mInitialMetric.Q(0, 0) * b(0, index_KL + j);
                    rB(1, index + j) = mInitialMetric.Q(1, 0) * b(0, index_KL + j) + mInitialMetric.Q(1, 1) * b(1, index_KL + j)
                        + mInitialMetric.Q(1, 2) * b(2, index_KL + j);
                    rB(2, index + j) = mInitialMetric.Q(2, 0) * b(0, index_KL + j) + mInitialMetric.Q(2, 2) * b(2, index_KL + j);
                }
                
            }

            for (unsigned int r = 0; r < mat_size_KL; r++){
                for (unsigned int i = 0; i < 3; i++){
                        if (b(i, r) < pow(10, -100) && b(i, r) > 0){
                            KRATOS_WATCH(i)
                            KRATOS_WATCH(r)
                            KRATOS_WATCH(b(i, r))
                        }
                }
            }
            for (unsigned int r = 0; r < mat_size; r++){
                for (unsigned int i = 0; i < 5; i++){
                        if (rB(i, r) < pow(10, -100) && rB(i, r) > 0){
                            KRATOS_WATCH(i)
                            KRATOS_WATCH(r)
                            KRATOS_WATCH(rB(i, r))
                        }
                }
            }
        }
        else
        {
            KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
        }
        // KRATOS_WATCH("End: CalculateBCurvature")
        KRATOS_CATCH("")
    }

    void IgaShell5pElementPreInt::CalculateSecondVariationStrainCurvature(
        SecondVariations& rSecondVariationsMembrane,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric)
    {
        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

            const unsigned int number_of_nodes = GetGeometry().size();
            const unsigned int mat_size_KL = number_of_nodes * 3;
            const unsigned int mat_size = number_of_nodes * 5;
           
            if (rSecondVariationsMembrane.B11.size1() != mat_size || rSecondVariationsMembrane.B11.size2() != mat_size)
                rSecondVariationsMembrane.B11.resize(mat_size, mat_size);
            rSecondVariationsMembrane.B11 = ZeroMatrix(mat_size, mat_size);
            if (rSecondVariationsMembrane.B22.size1() != mat_size || rSecondVariationsMembrane.B23.size2() != mat_size)
                rSecondVariationsMembrane.B22.resize(mat_size, mat_size);
            rSecondVariationsMembrane.B22 = ZeroMatrix(mat_size, mat_size);        
            if (rSecondVariationsMembrane.B12.size1() != mat_size || rSecondVariationsMembrane.B13.size2() != mat_size)
                rSecondVariationsMembrane.B12.resize(mat_size, mat_size);
            rSecondVariationsMembrane.B12 = ZeroMatrix(mat_size, mat_size);
            
            if (rSecondVariationsCurvature.B11.size1() != mat_size || rSecondVariationsCurvature.B11.size2() != mat_size)
                rSecondVariationsCurvature.B11.resize(mat_size, mat_size);
            rSecondVariationsCurvature.B11 = ZeroMatrix(mat_size, mat_size);
            if (rSecondVariationsCurvature.B22.size1() != mat_size || rSecondVariationsCurvature.B22.size2() != mat_size)
                rSecondVariationsCurvature.B22.resize(mat_size, mat_size);
            rSecondVariationsCurvature.B22 = ZeroMatrix(mat_size, mat_size);        
            if (rSecondVariationsCurvature.B12.size1() != mat_size || rSecondVariationsCurvature.B12.size2() != mat_size)
                rSecondVariationsCurvature.B12.resize(mat_size, mat_size);
            rSecondVariationsCurvature.B12 = ZeroMatrix(mat_size, mat_size);   

            double lg3_3 = pow(rMetric.dA, 3);
            double lg3_5 = pow(rMetric.dA, 5);
            double inv_lg3 = 1 / rMetric.dA;
            double inv_lg3_3 = 1 / lg3_3;
            double inv_lg3_5 = 1 / lg3_5;

            SecondVariations second_variations_membrane_KL(mat_size_KL);
            SecondVariations second_variations_curvature_KL(mat_size_KL);
            Matrix S_dg3 = ZeroMatrix(3, mat_size_KL);
            Vector S_g3dg3 = ZeroVector(mat_size_KL);
            Vector S_g3dg3lg3_3 = ZeroVector(mat_size_KL);
            Matrix S_dn = ZeroMatrix(3, mat_size_KL);
            // first variation of strain and curvature w.r.t. dof
            for (unsigned int r = 0; r < mat_size_KL; r++)
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

                S_g3dg3[r] = rMetric.g3_tilde[0] * S_dg3(0, r) + rMetric.g3_tilde[1] * S_dg3(1, r) + rMetric.g3_tilde[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.g3_tilde[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.g3_tilde[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.g3_tilde[2] * S_g3dg3lg3_3[r];
            }

            // second variation of strain and curvature w.r.t. dofs
            for (int r = 0; r < mat_size_KL; r++)
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

                        // calculated with simplified Q (ML)
                        second_variations_membrane_KL.B11(r, s) = mInitialMetric.Q(0, 0) * ddE_cu[0];
                        second_variations_membrane_KL.B11(s, r) = second_variations_membrane_KL.B11(r, s);
                        second_variations_membrane_KL.B22(r, s) = mInitialMetric.Q(1, 0) * ddE_cu[0] + mInitialMetric.Q(1, 1) * ddE_cu[1]
                            + mInitialMetric.Q(1, 2) * ddE_cu[2];
                        second_variations_membrane_KL.B22(s, r) = second_variations_membrane_KL.B22(r, s);
                        second_variations_membrane_KL.B12(r, s) = mInitialMetric.Q(2, 0) * ddE_cu[0] + mInitialMetric.Q(2, 2) * ddE_cu[2];
                        second_variations_membrane_KL.B12(s, r) = second_variations_membrane_KL.B12(r, s);
                    }

                    // curvature
                    array_1d<double, 3> ddg3 = ZeroVector(3);
                    int dirt = 4 - dirr - dirs;
                    int ddir = dirr - dirs;
                    if (ddir == -1)      ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 2) ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 1) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == -2) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);

                    double c = -(ddg3[0] * rMetric.g3_tilde[0] + ddg3[1] * rMetric.g3_tilde[1] + ddg3[2] * rMetric.g3_tilde[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.g3_tilde[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.g3_tilde[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.g3_tilde[2];

                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = - (DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2]);
                    ddK_cu[1] = - (DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2]);
                    ddK_cu[2] = - (DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2]);

                    // calculated with simplified Q (ML)
                    second_variations_curvature_KL.B11(r, s) = mInitialMetric.Q(0, 0) * ddK_cu[0];
                    second_variations_curvature_KL.B11(s, r) = second_variations_curvature_KL.B11(r, s);
                    second_variations_curvature_KL.B22(r, s) = mInitialMetric.Q(1, 0) * ddK_cu[0] + mInitialMetric.Q(1, 1) * ddK_cu[1] 
                        + mInitialMetric.Q(1, 2) * ddK_cu[2];
                    second_variations_curvature_KL.B22(s, r) = second_variations_curvature_KL.B22(r, s);
                    second_variations_curvature_KL.B12(r, s) = mInitialMetric.Q(2, 0) * ddK_cu[0] + mInitialMetric.Q(2, 2) * ddK_cu[2];
                    second_variations_curvature_KL.B12(s, r) = second_variations_curvature_KL.B12(r, s);
                }
            }

            // transfer KL-second-variations to RM-second-variations
            for (unsigned int r = 0; r < mat_size; r++) {
                unsigned int kr = r / 5;
                unsigned int dirr = r % 5;
                unsigned int r_KL = kr * 3 + dirr;
                if (dirr != 3 && dirr != 4){
                    for (unsigned int s = 0; s < mat_size; s++){
                        unsigned int ks = s / 5;
                        unsigned int dirs = s % 5;
                        unsigned int s_KL = ks * 3 + dirs;
                        if (dirs != 3 && dirs != 4){
                            rSecondVariationsMembrane.B11(r, s) = second_variations_membrane_KL.B11(r_KL, s_KL);
                            rSecondVariationsMembrane.B22(r, s) = second_variations_membrane_KL.B22(r_KL, s_KL);
                            rSecondVariationsMembrane.B12(r, s) = second_variations_membrane_KL.B12(r_KL, s_KL);
                            rSecondVariationsCurvature.B11(r, s) = second_variations_curvature_KL.B11(r_KL, s_KL);
                            rSecondVariationsCurvature.B22(r, s) = second_variations_curvature_KL.B22(r_KL, s_KL);
                            rSecondVariationsCurvature.B12(r, s) = second_variations_curvature_KL.B12(r_KL, s_KL);
                        }
                    }
                }
            }
        }
    }

    void IgaShell5pElementPreInt::CalculateVariationsRM(        
        Matrix& rBMembraneRM,
        Matrix& rBCurvatureRM,
        SecondVariations& rSecondVariationsMembraneRM,
        SecondVariations& rSecondVariationsCurvatureRM,
        const Vector& rShearDifferenceVector,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag)
    {
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;
        
        if (rBMembraneRM.size1() != 5 || rBMembraneRM.size2() != mat_size)
            rBMembraneRM.resize(5, mat_size);
        rBMembraneRM = ZeroMatrix(5, mat_size);
        if (rBCurvatureRM.size1() != 5 || rBCurvatureRM.size2() != mat_size)
            rBCurvatureRM.resize(5, mat_size);
        rBCurvatureRM = ZeroMatrix(5, mat_size);
        if (rSecondVariationsMembraneRM.B11.size1() != mat_size || rSecondVariationsMembraneRM.B11.size2() != mat_size)
            rSecondVariationsMembraneRM.B11.resize(mat_size, mat_size);
        rSecondVariationsMembraneRM.B11 = ZeroMatrix(mat_size, mat_size);
        if (rSecondVariationsMembraneRM.B23.size1() != mat_size || rSecondVariationsMembraneRM.B23.size2() != mat_size)
            rSecondVariationsMembraneRM.B23.resize(mat_size, mat_size);
        rSecondVariationsMembraneRM.B23 = ZeroMatrix(mat_size, mat_size);        
        if (rSecondVariationsMembraneRM.B13.size1() != mat_size || rSecondVariationsMembraneRM.B13.size2() != mat_size)
            rSecondVariationsMembraneRM.B13.resize(mat_size, mat_size);
        rSecondVariationsMembraneRM.B13 = ZeroMatrix(mat_size, mat_size);
        
        if (rSecondVariationsCurvatureRM.B11.size1() != mat_size || rSecondVariationsCurvatureRM.B11.size2() != mat_size)
            rSecondVariationsCurvatureRM.B11.resize(mat_size, mat_size);
        rSecondVariationsCurvatureRM.B11 = ZeroMatrix(mat_size, mat_size);
        if (rSecondVariationsCurvatureRM.B22.size1() != mat_size || rSecondVariationsCurvatureRM.B22.size2() != mat_size)
            rSecondVariationsCurvatureRM.B22.resize(mat_size, mat_size);
        rSecondVariationsCurvatureRM.B22 = ZeroMatrix(mat_size, mat_size);        
        if (rSecondVariationsCurvatureRM.B12.size1() != mat_size || rSecondVariationsCurvatureRM.B12.size2() != mat_size)
            rSecondVariationsCurvatureRM.B12.resize(mat_size, mat_size);
        rSecondVariationsCurvatureRM.B12 = ZeroMatrix(mat_size, mat_size);          
        // 1. First strain variation
        Matrix Dw_Dr = ZeroMatrix(3, mat_size);
        
        for (unsigned int r = 0; r < mat_size; r++){
            // local node number kr and dof direction dirr
            int kr = r / 5;
            int dirr = r % 5;
            
            // the two entries E23 and E13 w.r.t. the curvilinear coord. sys. are saved in dE_cur
            array_1d<double, 2> dE_cur = ZeroVector(2);
            
            if (dirr == 0 || dirr == 1 || dirr == 2){
                Dw_Dr(dirr, r) = rw_alpha(0) * DN_De(kr, 0) + rw_alpha(1) * DN_De(kr, 1);
                dE_cur[0] = 0.5 * (rShearDifferenceVector(dirr) * DN_De(kr, 1));
                dE_cur[1] = 0.5 * (rShearDifferenceVector(dirr) * DN_De(kr, 0));
            }
            else if(dirr == 3){
                Dw_Dr(0, r) = N(kr) * rActualMetric.g1(0);
                Dw_Dr(1, r) = N(kr) * rActualMetric.g1(1);
                Dw_Dr(2, r) = N(kr) * rActualMetric.g1(2);
            }
            else {
                Dw_Dr(0, r) = N(kr) * rActualMetric.g2(0);
                Dw_Dr(1, r) = N(kr) * rActualMetric.g2(1);
                Dw_Dr(2, r) = N(kr) * rActualMetric.g2(2);                
            }
            dE_cur[0] += 0.5 * (Dw_Dr(0, r) * rActualMetric.g2(0) + Dw_Dr(1, r) * rActualMetric.g2(1) + Dw_Dr(2, r) * rActualMetric.g2(2));
            dE_cur[1] += 0.5 * (Dw_Dr(0, r) * rActualMetric.g1(0) + Dw_Dr(1, r) * rActualMetric.g1(1) + Dw_Dr(2, r) * rActualMetric.g1(2));

            // calculated with the simplified Q (ML)
            rBMembraneRM(3, r) = mInitialMetric.Q(3, 3) * dE_cur[0] + mInitialMetric.Q(3, 4) * dE_cur[1];
            rBMembraneRM(4, r) = mInitialMetric.Q(4, 4) * dE_cur[1];
            // the other entries are (remain) zero
        
            // 2. First curvature variation
            array_1d<double, 3> dK_cu = ZeroVector(3);
            array_1d<double, 3> DDw_DD1r = ZeroVector(3);
            array_1d<double, 3> DDw_DD2r = ZeroVector(3);

            if (dirr == 0 || dirr == 1 || dirr == 2){
                DDw_DD1r[dirr] = rDw_alpha_Dbeta(0, 0) * DN_De(kr, 0)  + rw_alpha[0] * DDN_DDe(kr, 0)
                    + rDw_alpha_Dbeta(1, 0) * DN_De(kr, 1) + rw_alpha[1] * DDN_DDe(kr, 2);
                DDw_DD2r[dirr] = rDw_alpha_Dbeta(0, 1) * DN_De(kr, 0) + rw_alpha[0] * DDN_DDe(kr, 2)
                    + rDw_alpha_Dbeta(1, 1) * DN_De(kr, 1) + rw_alpha[1] * DDN_DDe(kr, 1);
                dK_cu[0] = rDw_D1[dirr] * DN_De(kr, 0);
                dK_cu[1] = rDw_D2[dirr] * DN_De(kr, 1);
                dK_cu[2] = 0.5 * (rDw_D1[dirr] * DN_De(kr, 1) + rDw_D2[dirr] * DN_De(kr, 0));
            }
            else if (dirr == 3){
                DDw_DD1r += DN_De(kr, 0) * rActualMetric.g1;
                DDw_DD1r[0] += N(kr) * rActualMetric.H(0, 0);
                DDw_DD1r[1] += N(kr) * rActualMetric.H(1, 0);
                DDw_DD1r[2] += N(kr) * rActualMetric.H(2, 0);
                DDw_DD2r += DN_De(kr, 1) * rActualMetric.g1;
                DDw_DD2r[0] += N(kr) * rActualMetric.H(0, 2);
                DDw_DD2r[1] += N(kr) * rActualMetric.H(1, 2);
                DDw_DD2r[2] += N(kr) * rActualMetric.H(2, 2);                
            }
            else if (dirr == 4){
                DDw_DD1r += DN_De(kr, 0) * rActualMetric.g2;
                DDw_DD1r[0] += N(kr) * rActualMetric.H(0, 2);
                DDw_DD1r[1] += N(kr) * rActualMetric.H(1, 2);
                DDw_DD1r[2] += N(kr) * rActualMetric.H(2, 2);
                DDw_DD2r += DN_De(kr, 1) * rActualMetric.g2;
                DDw_DD2r[0] += N(kr) * rActualMetric.H(0, 1);
                DDw_DD2r[1] += N(kr) * rActualMetric.H(1, 1);
                DDw_DD2r[2] += N(kr) * rActualMetric.H(2, 1);
            }
            dK_cu[0] += inner_prod(DDw_DD1r, rActualMetric.g1);
            dK_cu[1] += inner_prod(DDw_DD2r, rActualMetric.g2);
            dK_cu[2] += 0.5 * (inner_prod(DDw_DD1r, rActualMetric.g2) + inner_prod(DDw_DD2r, rActualMetric.g1));

            // calculated with simplified Q (ML)
            rBCurvatureRM(0, r) = mInitialMetric.Q(0, 0) * dK_cu[0];
            rBCurvatureRM(1, r) = mInitialMetric.Q(1, 0) * dK_cu[0] + mInitialMetric.Q(1, 1) * dK_cu[1] + mInitialMetric.Q(1, 2) * dK_cu[2];
            rBCurvatureRM(2, r) = mInitialMetric.Q(2, 0) * dK_cu[0] + mInitialMetric.Q(2, 2) * dK_cu[2];
            // all other entries are (remain) zero

            // 3. Second Strain Variation
            if (rCalculateStiffnessMatrixFlag == true){
                for (unsigned int s = 0; s <= r; s++)
                {
                    // local node number ks and dof direction dirs
                    int ks = s / 5;
                    int dirs = s % 5;
                    
                    array_1d <double, 3> DDw_DDrs = ZeroVector(3);
                    array_1d <double, 2> ddE_cu = ZeroVector(2);

                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        ddE_cu[0] = 0.5 * (Dw_Dr(dirr, r) * DN_De(kr, 1));
                        ddE_cu[1] = 0.5 * (Dw_Dr(dirr, r) * DN_De(kr, 0));
                    }
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        ddE_cu[0] += 0.5 * (Dw_Dr(dirs, r) * DN_De(ks, 1));
                        ddE_cu[1] += 0.5 * (Dw_Dr(dirs, r) * DN_De(ks, 0));
                    }
                    if (dirr == 3 && (dirs == 0 || dirs == 1 || dirs == 2))
                        DDw_DDrs(dirs) += N(kr) * DN_De(ks, 0);
                    else if (dirr == 4 && (dirs == 0 || dirs == 1 || dirs == 2))
                        DDw_DDrs(dirs) += N(kr) * DN_De(ks, 1);
                    if (dirs == 3 && (dirr == 0 || dirr == 1 || dirr == 2))
                        DDw_DDrs(dirr) += N(ks) * DN_De(kr, 0);       
                    else if (dirs == 4 && (dirr == 0 || dirr == 1 || dirr == 2))
                        DDw_DDrs(dirr) += N(ks) * DN_De(kr, 1);
                    ddE_cu[0] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.g2);
                    ddE_cu[1] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.g1);

                    // calculated with simplified Q (ML)
                    rSecondVariationsMembraneRM.B23(r, s) = mInitialMetric.Q(3, 3) * ddE_cu[0] + mInitialMetric.Q(3, 4) * ddE_cu[1];
                    rSecondVariationsMembraneRM.B23(s, r) = rSecondVariationsMembraneRM.B23(r, s);
                    rSecondVariationsMembraneRM.B13(r, s) = mInitialMetric.Q(4, 4) * ddE_cu[1];
                    rSecondVariationsMembraneRM.B13(s, r) = rSecondVariationsMembraneRM.B13(r, s);
                    // all other entries are (remain) zero
                    
                    // 4. Second curvature variation
                    array_1d <double, 3> DDDw_DDD1rs = ZeroVector(3);
                    array_1d <double, 3> DDDw_DDD2rs = ZeroVector(3);
                    array_1d <double, 3> ddK_cu = ZeroVector(3);
                    array_1d<double, 3> DDw_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDw_DD2s = ZeroVector(3);

                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        DDw_DD1s[dirs] = rDw_alpha_Dbeta(0, 0) * DN_De(ks, 0)  + rw_alpha[0] * DDN_DDe(ks, 0)
                        + rDw_alpha_Dbeta(1, 0) * DN_De(ks, 1) + rw_alpha[1] * DDN_DDe(ks, 2);
                        DDw_DD2s[dirs] = rDw_alpha_Dbeta(0, 1) * DN_De(ks, 0) + rw_alpha[0] * DDN_DDe(ks, 2)
                        + rDw_alpha_Dbeta(1, 1) * DN_De(ks, 1) + rw_alpha[1] * DDN_DDe(ks, 1);
                    }
                    else if (dirs == 3){
                        DDw_DD1s += DN_De(ks, 0) * rActualMetric.g1;
                        DDw_DD1s[0] += N(ks) * rActualMetric.H(0, 0);
                        DDw_DD1s[1] += N(ks) * rActualMetric.H(1, 0);
                        DDw_DD1s[2] += N(ks) * rActualMetric.H(2, 0);
                        DDw_DD2s += DN_De(ks, 1) * rActualMetric.g1;
                        DDw_DD2s[0] += N(ks) * rActualMetric.H(0, 2);
                        DDw_DD2s[1] += N(ks) * rActualMetric.H(1, 2);
                        DDw_DD2s[2] += N(ks) * rActualMetric.H(2, 2);                
                    }
                    else if (dirs == 4){
                        DDw_DD1s += DN_De(ks, 0) * rActualMetric.g2;
                        DDw_DD1s[0] += N(ks) * rActualMetric.H(0, 2);
                        DDw_DD1s[1] += N(ks) * rActualMetric.H(1, 2);
                        DDw_DD1s[2] += N(ks) * rActualMetric.H(2, 2);
                        DDw_DD2s += DN_De(ks, 1) * rActualMetric.g2;
                        DDw_DD2s[0] += N(ks) * rActualMetric.H(0, 1);
                        DDw_DD2s[1] += N(ks) * rActualMetric.H(1, 1);
                        DDw_DD2s[2] += N(ks) * rActualMetric.H(2, 1);
                    }
                    
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        ddK_cu[0] = DDw_DD1s[dirr] * DN_De(kr, 0);
                        ddK_cu[1] = DDw_DD2s[dirr] * DN_De(kr, 1);
                        ddK_cu[2] = 0.5 * (DDw_DD1s[dirr] * DN_De(kr, 1) + DDw_DD2s[dirr] * DN_De(kr, 0));
                        }
                    else if (dirr == 3 && (dirs == 0 || dirs == 1 || dirs == 2)){
                        DDDw_DDD1rs[dirs] = DN_De(kr, 0) * DN_De(ks, 0) + N(kr) * DDN_DDe(ks, 0);
                        DDDw_DDD2rs[dirs] = DN_De(kr, 1) * DN_De(ks, 0) + N(kr) * DDN_DDe(ks, 2);
                    }
                    else if (dirr == 4 && (dirs == 0 || dirs == 1 || dirs == 2)){
                        DDDw_DDD1rs[dirs] += DN_De(kr, 0) * DN_De(ks, 1) + N(kr) * DDN_DDe(ks, 2);
                        DDDw_DDD2rs[dirs] += DN_De(kr, 1) * DN_De(ks, 1) + N(kr) * DDN_DDe(ks, 1);
                    }
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        ddK_cu[0] += DDw_DD1r[dirs] * DN_De(ks, 0);
                        ddK_cu[1] += DDw_DD2r[dirs] * DN_De(ks, 1);
                        ddK_cu[2] += 0.5 * (DDw_DD1r[dirs] * DN_De(ks, 1) + DDw_DD2r[dirs] * DN_De(ks, 0));
                        }
                    else if (dirs == 3 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += DN_De(ks, 0) * DN_De(kr, 0) + N(ks) * DDN_DDe(kr, 0);
                        DDDw_DDD2rs[dirr] += DN_De(ks, 1) * DN_De(kr, 0) + N(ks) * DDN_DDe(kr, 2);
                    }
                    else if (dirs == 4 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += DN_De(ks, 0) * DN_De(kr, 1) + N(ks) * DDN_DDe(kr, 2);
                        DDDw_DDD2rs[dirr] += DN_De(ks, 1) * DN_De(kr, 1) + N(ks) * DDN_DDe(kr, 1);
                    }
                    ddK_cu[0] += inner_prod(DDDw_DDD1rs, rActualMetric.g1);
                    ddK_cu[1] += inner_prod(DDDw_DDD2rs, rActualMetric.g2);
                    ddK_cu[2] += 0.5 * (inner_prod(DDDw_DDD1rs, rActualMetric.g2) + inner_prod(DDDw_DDD2rs, rActualMetric.g1));

                    // calculated with simplified Q (ML)
                    rSecondVariationsCurvatureRM.B11(r, s) = mInitialMetric.Q(0, 0) * ddK_cu[0];
                    rSecondVariationsCurvatureRM.B11(s, r) = rSecondVariationsCurvatureRM.B11(r, s);
                    rSecondVariationsCurvatureRM.B22(r, s) = mInitialMetric.Q(1, 0) * ddK_cu[0] + mInitialMetric.Q(1, 1) * ddK_cu[1]
                        + mInitialMetric.Q(1, 2) * ddK_cu[2];
                    rSecondVariationsCurvatureRM.B22(s, r) = rSecondVariationsCurvatureRM.B22(r, s);
                    rSecondVariationsCurvatureRM.B12(r, s) = mInitialMetric.Q(2, 0) * ddK_cu[0] + mInitialMetric.Q(2, 2) * ddK_cu[2];
                    rSecondVariationsCurvatureRM.B12(s, r) = rSecondVariationsCurvatureRM.B12(r, s);
                    // all other entries are (remain) zero
                }
            }
        }
        // KRATOS_WATCH("end: CalculateVariationsRM")
    }

    void IgaShell5pElementPreInt::Calculate(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        // if (rValues.size() != 1)
        // {
        //     rValues.resize(1);
        // }

        // if (rVariable == STRESSES)
        // {
        //     // Create constitutive law parameters:
        //     ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        //     // Set constitutive law flags:
        //     Flags& ConstitutiveLawOptions = Values.GetOptions();
        //     ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        //     ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        //     ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        //     MetricVariables actual_metric(3);
        //     CalculateMetric(actual_metric);
        //     ConstitutiveVariables constitutive_variables_membrane(3);
        //     ConstitutiveVariables constitutive_variables_curvature(3);
        //     CalculateConstitutiveVariables(actual_metric,
        //         constitutive_variables_membrane, constitutive_variables_curvature,
        //         Values, ConstitutiveLaw::StressMeasure_PK2);

        //     double detF = actual_metric.dA / mInitialMetric.dA;

        //     Vector n_pk2_ca = prod(constitutive_variables_membrane.D, constitutive_variables_membrane.E);
        //     Vector n_pk2_con = prod(mInitialMetric.T, n_pk2_ca);
        //     Vector n_cau = 1.0 / detF*n_pk2_con;

        //     Vector n = ZeroVector(8);
        //     // Cauchy normal force in normalized g1,g2
        //     n[0] = std::sqrt(actual_metric.gab[0] / actual_metric.gab_con[0])*n_cau[0];
        //     n[1] = std::sqrt(actual_metric.gab[1] / actual_metric.gab_con[1])*n_cau[1];
        //     n[2] = std::sqrt(actual_metric.gab[0] / actual_metric.gab_con[1])*n_cau[2];
        //     // Cauchy normal force in local cartesian e1,e2
        //     array_1d<double, 3> n_e = prod(actual_metric.T, n_cau);
        //     n[3] = n_e[0];
        //     n[4] = n_e[1];
        //     n[5] = n_e[2];
        //     // Principal normal forces
        //     n[6] = 0.5*(n_e[0] + n_e[1] + std::sqrt(pow(n_e[0] - n_e[1], 2) + 4.0*pow(n_e[2], 2)));
        //     n[7] = 0.5*(n_e[0] + n_e[1] - std::sqrt(pow(n_e[0] - n_e[1], 2) + 4.0*pow(n_e[2], 2)));

        //     // -------------------  moments -------------------------
        //     // PK2 moment in local cartesian E1,E2
        //     Vector m_pk2_local_ca = prod(constitutive_variables_curvature.D, constitutive_variables_curvature.E);
        //     // PK2 moment in G1,G2
        //     Vector m_pk2_con = prod(mInitialMetric.T, m_pk2_local_ca);
        //     // Cauchy moment in g1,g2
        //     array_1d<double, 3> m_cau = 1.0 / detF*m_pk2_con;


        //     Vector m = ZeroVector(8);
        //     // Cauchy moment in normalized g1,g2
        //     m[0] = std::sqrt(actual_metric.gab[0] / actual_metric.gab_con[0])*m_cau[0];
        //     m[1] = std::sqrt(actual_metric.gab[1] / actual_metric.gab_con[1])*m_cau[1];
        //     m[2] = std::sqrt(actual_metric.gab[0] / actual_metric.gab_con[1])*m_cau[2];
        //     // Cauchy moment in local cartesian e1,e2
        //     Vector m_e = prod(actual_metric.T, m_cau);
        //     m[3] = m_e[0];
        //     m[4] = m_e[1];
        //     m[5] = m_e[2];
        //     // principal moment
        //     m[6] = 0.5*(m_e[0] + m_e[1] + std::sqrt(pow(m_e[0] - m_e[1], 2) + 4.0*pow(m_e[2], 2)));
        //     m[7] = 0.5*(m_e[0] + m_e[1] - std::sqrt(pow(m_e[0] - m_e[1], 2) + 4.0*pow(m_e[2], 2)));

        //     double thickness = this->GetProperties().GetValue(THICKNESS);

        //     double W = pow(thickness, 2) / 6.0;
        //     Vector sigma_top = ZeroVector(3);
        //     sigma_top(0) = m[3] / W + n[3] / thickness;
        //     sigma_top(1) = m[4] / W + n[4] / thickness;
        //     sigma_top(2) = m[5] / W + n[5] / thickness;

        //     rValues[0] = sigma_top;
        // }
        // else if (rVariable == EXTERNAL_FORCES_VECTOR) {     // there is also a INTERNAL_FORCES_VECTOR (ML)
        //     const int& number_of_nodes = GetGeometry().size();
        //     Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
        //     Vector external_forces_vector = ZeroVector(3);
        //     for (int i = 0; i < number_of_nodes; i++)
        //     {
        //         const NodeType & iNode = GetGeometry()[i];
        //         const Vector& forces = iNode.GetValue(EXTERNAL_FORCES_VECTOR);

        //         external_forces_vector[0] += N[i] * forces[0];
        //         external_forces_vector[1] += N[i] * forces[1];
        //         external_forces_vector[2] += N[i] * forces[2];
        //     }
        //     rValues[0] = external_forces_vector;
        // }
        // else
        // {
        //     rValues[0] = ZeroVector(3);
        // }
    }

    void IgaShell5pElementPreInt::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        // KRATOS_WATCH("EquationIdVector")

        const unsigned int number_of_nodes = GetGeometry().size();

        if (rResult.size() != 5 * number_of_nodes)
            rResult.resize(5 * number_of_nodes, false);

        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
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
    }

    void IgaShell5pElementPreInt::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        // KRATOS_WATCH("GetDofList")

        const unsigned int number_of_nodes = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
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
    }

    int IgaShell5pElementPreInt::Check(const ProcessInfo& rCurrentProcessInfo)
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
