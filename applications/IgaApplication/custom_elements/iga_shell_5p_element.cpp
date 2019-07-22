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

        // KRATOS_WATCH("start: Initialize")
        // Constitutive Law initialisation
        BaseDiscreteElement::Initialize();
        // Check whether ConstitutiveLaw is 3D
        if (mConstitutiveLawVector[0]->GetStrainSize() != 6){
            KRATOS_WATCH("ConstitutiveLaw is not 3D.")
            KRATOS_ERROR << "ConstitutiveLaw is not 3D." << std::endl;
        }

        CalculateMetric(mInitialMetric);
        
        mZeta = 0.0;

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

        // shear difference vector
        array_1d<double, 3> w = ZeroVector(3);
        // derivatives of the shear difference vector
        array_1d<double, 3> Dw_D1 = ZeroVector(3);
        array_1d<double, 3> Dw_D2 = ZeroVector(3);
        // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * a1 + w_alpha(2) * a2)
        array_1d<double, 2> w_alpha = ZeroVector(2);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);

        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);
        double dV = 0.0;
        array_1d<double, 3> G1 = ZeroVector(3);
        array_1d<double, 3> G2 = ZeroVector(3);
        array_1d<double, 3> G1xG2 = ZeroVector(3);
        double thickness = GetProperties().GetValue(THICKNESS);

        for (unsigned int Gauss_index = 0; Gauss_index < mGaussQuadratureThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussQuadratureThickness.zeta(Gauss_index);
            
            // Differential Volume
            CalculateInitialBaseVectorsGLinearized(mZeta, G1, G2);
            MathUtils<double>::CrossProduct(G1xG2, G1, G2);
            dV = inner_prod(G1xG2, mInitialMetric.a3_KL);

            Matrix B = ZeroMatrix(5, mat_size);
            SecondVariations second_variations(mat_size);

            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, w, Dw_D1, Dw_D2, constitutive_variables, Values, 
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            CalculateB(B, actual_metric);
            if(Id()==4)
                KRATOS_WATCH(B)
            CalculateVariationsRM(B, second_variations, w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, 
                actual_metric, CalculateStiffnessMatrixFlag);
            if(Id()==4)
                KRATOS_WATCH(B)
                
            double integration_weight = mGaussQuadratureThickness.integration_weight_thickness(Gauss_index) * 
                GetValue(INTEGRATION_WEIGHT) * dV * thickness / 2.0;

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                // Nonlinear Deformation
                CalculateSecondVariations(
                    second_variations,
                    actual_metric);
                
                // adding linear contributions to the stiffness matrix
                CalculateAndAddKm(rLeftHandSideMatrix, B, constitutive_variables.D, integration_weight);

                // adding  non-linear contribution to stiffness matrix
                CalculateAndAddNonlinearKm(rLeftHandSideMatrix, second_variations, constitutive_variables.S,
                    integration_weight);
            }

            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), constitutive_variables.S);
            }
        }

        // if (Id() == 1){
        //     KRATOS_WATCH(rLeftHandSideMatrix)
        //     KRATOS_WATCH(rRightHandSideVector)
        // }

        KRATOS_CATCH("");
    }

    void IgaShell5pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double& rIntegrationWeight
    )
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += rIntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        // if (Id() == 4){
        //     KRATOS_WATCH(rIntegrationWeight)
        //     KRATOS_WATCH(B)
        //     KRATOS_WATCH(rLeftHandSideMatrix)
        // }
        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateAndAddNonlinearKm(
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

    void IgaShell5pElement::CalculateMetric(
        MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(), DN_De, 3, 2, rMetric.J);

        rMetric.a1[0] = rMetric.J(0, 0);
        rMetric.a2[0] = rMetric.J(0, 1);
        rMetric.a1[1] = rMetric.J(1, 0);
        rMetric.a2[1] = rMetric.J(1, 1);
        rMetric.a1[2] = rMetric.J(2, 0);
        rMetric.a2[2] = rMetric.J(2, 1);

        //basis vector a3_KL
        MathUtils<double>::CrossProduct(rMetric.a3_KL_tilde, rMetric.a1, rMetric.a2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.a3_KL_tilde);
        //normalized basis vector a3_KL
        rMetric.a3_KL = rMetric.a3_KL_tilde / rMetric.dA;
        
        //GetCovariantMetric
        rMetric.a_ab[0] = pow(rMetric.a1[0], 2) + pow(rMetric.a1[1], 2) + pow(rMetric.a1[2], 2);
        rMetric.a_ab[1] = pow(rMetric.a2[0], 2) + pow(rMetric.a2[1], 2) + pow(rMetric.a2[2], 2);
        rMetric.a_ab[2] = rMetric.a1[0] * rMetric.a2[0] + rMetric.a1[1] * rMetric.a2[1] + rMetric.a1[2] * rMetric.a2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);
        
        // base vector derivatives
        for (unsigned int i = 0; i < 3; i++)
        {
            rMetric.Da1_D1[i] = rMetric.H(i, 0);
            rMetric.Da2_D2[i] = rMetric.H(i, 1);
            rMetric.Da1_D2[i] = rMetric.H(i, 2);
        }

        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.a3_KL[0] + rMetric.H(1, 0) * rMetric.a3_KL[1] + rMetric.H(2, 0) * rMetric.a3_KL[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.a3_KL[0] + rMetric.H(1, 1) * rMetric.a3_KL[1] + rMetric.H(2, 1) * rMetric.a3_KL[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.a3_KL[0] + rMetric.H(1, 2) * rMetric.a3_KL[1] + rMetric.H(2, 2) * rMetric.a3_KL[2];

        //contravariant rMetric a_ab_con and base vectors g_con
        double inv_deta_ab = 1.0 / (rMetric.a_ab[0] * rMetric.a_ab[1] - rMetric.a_ab[2] * rMetric.a_ab[2]);
        rMetric.a_ab_con[0] = inv_deta_ab*rMetric.a_ab[1];
        rMetric.a_ab_con[2] = -inv_deta_ab*rMetric.a_ab[2];
        rMetric.a_ab_con[1] = inv_deta_ab*rMetric.a_ab[0];

        rMetric.a1_con = rMetric.a1*rMetric.a_ab_con[0] + rMetric.a2*rMetric.a_ab_con[2];
        rMetric.a2_con = rMetric.a1*rMetric.a_ab_con[2] + rMetric.a2*rMetric.a_ab_con[1];
        // g_con_3 = a3_KL
        
        //local cartesian coordinates
        double lg1 = norm_2(rMetric.a1);
        array_1d<double, 3> e1 = rMetric.a1 / lg1;
        double lg_con2 = norm_2(rMetric.a2_con);
        array_1d<double, 3> e2 = rMetric.a2_con / lg_con2;
        // e3 = a3_KL = g_con_3

        // // transformation matrix Q from curvilinear to local cartesian coordinate system
        // faster computation of transformation matrix Q taking into account that a lot of entries become zero (ML)
        // this matrix Q is referring to a VoigtSize 5 with E11, E22, E12, E23, E13
        double mG_00 = inner_prod(e1, rMetric.a1_con);
        double mG_10 = inner_prod(e2, rMetric.a1_con);
        double mG_11 = inner_prod(e2, rMetric.a2_con);
        
        rMetric.Q(0, 0) = pow(mG_00, 2);
        rMetric.Q(1, 0) = pow(mG_10, 2);
        rMetric.Q(1, 1) = pow(mG_11, 2);
        rMetric.Q(1, 2) = 2.00 * mG_10 * mG_11;
        rMetric.Q(2, 0) = 2.00 * mG_00 * mG_10;
        rMetric.Q(2, 2) = 2.00 * mG_00 * mG_11;
        rMetric.Q(3, 3) = 2.00 * mG_11;
        rMetric.Q(3, 4) = 2.00 * mG_10;
        rMetric.Q(4, 4) = 2.00 * mG_00;

        // transformation matrix TransCartToCov from local Cartesian to covariant basis
        rMetric.TransCartToCov = trans(rMetric.Q);
        // division by 2.0 because not used for strains but for stresses (strains have e.g. entries with 2e12)
        rMetric.TransCartToCov(2, 1) = rMetric.TransCartToCov(2, 1) / 2.0;
        rMetric.TransCartToCov(2, 2) = rMetric.TransCartToCov(2, 2) / 2.0;
        rMetric.TransCartToCov(3, 3) = rMetric.TransCartToCov(3, 3) / 2.0;
        rMetric.TransCartToCov(4, 3) = rMetric.TransCartToCov(4, 3) / 2.0;
        rMetric.TransCartToCov(4, 4) = rMetric.TransCartToCov(4, 4) / 2.0;

        mG_00 = inner_prod(e1, rMetric.a1);
        double mG_01 = inner_prod(e1, rMetric.a2);
        mG_11 = inner_prod(e2, rMetric.a2);
        // transformation matrix TransCovToCart from covariant to local Cartesian basis
        rMetric.TransCovToCart(0, 0) = pow(mG_00, 2);
        rMetric.TransCovToCart(0, 1) = pow(mG_01, 2);
        rMetric.TransCovToCart(0, 2) = 2 * mG_00 * mG_01;
        rMetric.TransCovToCart(1, 1) = pow(mG_11, 2);
        rMetric.TransCovToCart(2, 1) = mG_01 * mG_11;
        rMetric.TransCovToCart(2, 2) = mG_00 * mG_11;
        rMetric.TransCovToCart(3, 3) = mG_11;
        rMetric.TransCovToCart(4, 3) = mG_01;
        rMetric.TransCovToCart(4, 4) = mG_00;
    }

    void IgaShell5pElement::CalculateShearDifferenceVector(
        array_1d<double, 3>& rw,
        array_1d<double, 3>& rDw_D1,
        array_1d<double, 3>& rDw_D2,
        array_1d<double, 2>& rw_alpha,
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
        rDw_D1 = rDw_alpha_Dbeta(0, 0) * rActualMetric.a1 + rDw_alpha_Dbeta(1, 0) * rActualMetric.a2;
        rDw_D2 = rDw_alpha_Dbeta(0, 1) * rActualMetric.a1 + rDw_alpha_Dbeta(1, 1) * rActualMetric.a2;

        for (unsigned int i = 0; i < 3; i++)
        {
            rDw_D1[i] += rw_alpha(0) * rActualMetric.H(i, 0) + rw_alpha(1) * rActualMetric.H(i, 2);
            rDw_D2[i] += rw_alpha(0) * rActualMetric.H(i, 2) + rw_alpha(1) * rActualMetric.H(i, 1);
        }

        rw = rw_alpha(0) * rActualMetric.a1 + rw_alpha(1) * rActualMetric.a2;

        KRATOS_CATCH("")
    }
    
    void IgaShell5pElement::CalculateInitialBaseVectorsGLinearized(
        const double& rZeta,
        array_1d<double, 3>&      rG1,
        array_1d<double, 3>&      rG2)
    {
        Vector DA3_D1 = ZeroVector(3);
        Vector DA3_D2 = ZeroVector(3);
        Vector DA1_D1xA2 = ZeroVector(3);
        Vector A1xDA2_D1 = ZeroVector(3);
        Vector DA1_D2xA2 = ZeroVector(3);
        Vector A1xDA2_D2 = ZeroVector(3);
        Vector G1xG2 = ZeroVector(3);

        MathUtils<double>::CrossProduct(DA1_D1xA2, mInitialMetric.Da1_D1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D1, mInitialMetric.a1, mInitialMetric.Da1_D2); // DA1_D2 = DA2_D1
        MathUtils<double>::CrossProduct(DA1_D2xA2, mInitialMetric.Da1_D2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D2, mInitialMetric.a1, mInitialMetric.Da2_D2);
        DA3_D1 = ((DA1_D1xA2 + A1xDA2_D1) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(DA1_D1xA2 + A1xDA2_D1)) 
            / (mInitialMetric.dA * mInitialMetric.dA);
        DA3_D2 = ((DA1_D2xA2 + A1xDA2_D2) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(DA1_D2xA2 + A1xDA2_D2))
            / (mInitialMetric.dA * mInitialMetric.dA);

        // covariant base vectors of the shell body in the reference configuration
        rG1 = mInitialMetric.a1 + rZeta * DA3_D1;
        rG2 = mInitialMetric.a2 + rZeta * DA3_D2;
        // G3 = A3
    }

    void IgaShell5pElement::CalculateActualBaseVectorsgLinearized(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        array_1d<double, 3>&      rg1,            ///< erster Basisvektor (o)
        array_1d<double, 3>&      rg2,            ///< zweiter Basisvektor (o)
        array_1d<double, 3>&      rg3             ///< dritter Basisvektor (o)
        )
    {
        KRATOS_TRY;
        
        const int number_of_control_points = GetGeometry().size();
        double thickness = GetProperties().GetValue(THICKNESS);

        /* AKTUELLE KONFIGURATION */
        double DdA_D1 = 0.0;
        double DdA_D2 = 0.0;
        array_1d<double, 3> Da1xa2_D1 = ZeroVector(3);
        array_1d<double, 3> Da1xa2_D2 = ZeroVector(3);
        array_1d<double, 3> Da3_KL_D1 = ZeroVector(3);
        array_1d<double, 3> Da3_KL_D2 = ZeroVector(3);

        // shape functions
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        /* Ableitung von a3_KL nach der Richtung alpha */
        /* Ableitung des Zaehlers von a3_KL nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> Da1_D1xa2, a1xDa2_D1;
        MathUtils<double>::CrossProduct(Da1_D1xa2, rActualMetric.Da1_D1, rActualMetric.a2);
        MathUtils<double>::CrossProduct(a1xDa2_D1, rActualMetric.a1, rActualMetric.Da1_D2);
        Da1xa2_D1 = Da1_D1xa2 + a1xDa2_D1;

        /* Ableitung des Zaehlers von a3_KL nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> Da1_D2xa2, a1xDa2_D2;
        MathUtils<double>::CrossProduct(Da1_D2xa2, rActualMetric.Da1_D2, rActualMetric.a2);
        MathUtils<double>::CrossProduct(a1xDa2_D2, rActualMetric.a1, rActualMetric.Da2_D2);
        Da1xa2_D2 = Da1_D2xa2 + a1xDa2_D2;

        /* Ableitung des Nenners von a3_KL nach 1 */
        DdA_D1 = (Da1xa2_D1[0] * rActualMetric.a3_KL_tilde[0] 
            + Da1xa2_D1[1] * rActualMetric.a3_KL_tilde[1] + Da1xa2_D1[2] * rActualMetric.a3_KL_tilde[2]) / rActualMetric.dA;

        /* Ableitung des Nenners von a3_KL nach 2 */
        DdA_D2 = (Da1xa2_D2[0] * rActualMetric.a3_KL_tilde[0] 
            + Da1xa2_D2[1] * rActualMetric.a3_KL_tilde[1] + Da1xa2_D2[2] * rActualMetric.a3_KL_tilde[2]) / rActualMetric.dA;

        /* Ableitung von a3_KL mit Quotientenregel */
        for (unsigned int i = 0; i < 3; i++)
        {
            Da3_KL_D1[i] = (Da1xa2_D1[i] * rActualMetric.dA - rActualMetric.a3_KL_tilde[i] * DdA_D1) / 
                (rActualMetric.dA * rActualMetric.dA);
            Da3_KL_D2[i] = (Da1xa2_D2[i] * rActualMetric.dA - rActualMetric.a3_KL_tilde[i] * DdA_D2) / 
                (rActualMetric.dA * rActualMetric.dA);
        }

        /* Kovariante Basisvektoren */
        rg1 = rActualMetric.a1 + mZeta*(Da3_KL_D1 + rDw_D1);
        rg2 = rActualMetric.a2 + mZeta*(Da3_KL_D2 + rDw_D2);
        rg3 = rActualMetric.a3_KL + rw;     // g3 = a3

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateDeformationGradient(
        const array_1d<double, 3> rG1,
        const array_1d<double, 3> rG2,
        const array_1d<double, 3> rg1,
        const array_1d<double, 3> rg2,
        const array_1d<double, 3> rg3,
        Matrix& rF,
        double& rdetF)
    {
        array_1d<double, 3> G1con(3);                // Kontravariante Basisvektoren
        array_1d<double, 3> G2con(3);                // Kontravariante Basisvektoren
        array_1d<double, 3> G3con(3);                // Kontravariante Basisvektoren

        /* - Initialisieren aller gesuchten Groessen  */
        rF.resize(3, 3);           // Deformationsgradient
        rF = ZeroMatrix(3, 3);
        rdetF = 0.0;                 // Determinante des Deformationsgradienten

        // Covariant Metric
        array_1d<double, 3> G_ab = ZeroVector(3);
        G_ab[0] = pow(mInitialMetric.a1[0], 2) + pow(mInitialMetric.a1[1], 2) + pow(mInitialMetric.a1[2], 2);
        G_ab[1] = pow(mInitialMetric.a2[0], 2) + pow(mInitialMetric.a2[1], 2) + pow(mInitialMetric.a2[2], 2);
        G_ab[2] = mInitialMetric.a1[0] * mInitialMetric.a2[0] + mInitialMetric.a1[1] * mInitialMetric.a2[1] 
            + mInitialMetric.a1[2] * mInitialMetric.a2[2];
        // Contravariant Metric
        array_1d<double, 3> G_ab_con = ZeroVector(3);
        double inv_detG_ab = 1.0 / (G_ab[0] * G_ab[1] - G_ab[2] * G_ab[2]);
        G_ab_con[0] = inv_detG_ab * G_ab[1];
        G_ab_con[2] = -inv_detG_ab * G_ab[2];
        G_ab_con[1] = inv_detG_ab * G_ab[0];
        // Contravariant base vectors
        array_1d<double, 3> G1_con = ZeroVector(3);
        array_1d<double, 3> G2_con = ZeroVector(3);
        G1_con = rG1 * G_ab_con[0] + rG2 * G_ab_con[2];
        G2_con = rG1 * G_ab_con[2] + rG2 * G_ab_con[1];
    
        // deformation gradient and its determinant
        rF = outer_prod(rg1, G1_con) + outer_prod(rg2, G2_con) + outer_prod(rg3, mInitialMetric.a3_KL);
        rdetF = rF(0, 0) * (rF(1, 1) * rF(2, 2) - rF(2, 1) * rF(1, 2)) 
            - rF(1, 0) * (rF(0, 1) * rF(2, 2) - rF(2, 1) * rF(0, 2))
            + rF(2, 0) * (rF(0, 1) * rF(1, 2) - rF(1, 1) * rF(0, 2));
    }

    void IgaShell5pElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        KRATOS_TRY

        Vector strain_vector = ZeroVector(5);
        Vector strain_vector_RM = ZeroVector(5);
        
        // Strain computation in curvilinear space
        CalculateStrain(strain_vector, rActualMetric.a_ab, rActualMetric.curvature);
        CalculateStrainRM(strain_vector_RM, rw, rDw_D1, rDw_D2, rActualMetric.a1, rActualMetric.a2);
        rThisConstitutiveVariables.E = strain_vector + strain_vector_RM;
  
        // Strain transformation to local Cartesian Space with VoigtSize 6 because ConstitutiveLaw is 3D
        ConstitutiveVariables constitutive_variables(6);
        TransformationCurvilinearStrainSize5ToCartesianStrainSize6(rThisConstitutiveVariables.E, constitutive_variables.E);

        //Constitutive Matrix D
        rValues.SetStrainVector(constitutive_variables.E); //this is the input parameter
        rValues.SetStressVector(constitutive_variables.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(constitutive_variables.D); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        // static condensation of  sigma_33
        unsigned int index_i = 0;
        for (unsigned int i = 0; i < 6; i++)
        {
            if (i != 2){
                unsigned int index_j = 0;
                for (unsigned int j = 0; j < 6; j++ ){
                    if (j != 2){
                        rThisConstitutiveVariables.D(index_i, index_j) += constitutive_variables.D(i, j) - 
                            constitutive_variables.D(i, 2) * constitutive_variables.D(2, j) 
                            / constitutive_variables.D(2, 2);
                        index_j++;
                    }
                }
                index_i++;
            }
        }
        
        // Strain Transformation to local Cartesian space with VoigtSize 5
        rThisConstitutiveVariables.E = prod(mInitialMetric.Q, rThisConstitutiveVariables.E);

        //Local Cartesian Stresses
        rThisConstitutiveVariables.S = prod(
            trans(rThisConstitutiveVariables.D), rThisConstitutiveVariables.E);

        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab,
        const Vector& rCurvature)
    {
        KRATOS_TRY
        
        double thickness = GetProperties().GetValue(THICKNESS);

        rStrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.a_ab[0]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[0] - rCurvature[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.a_ab[1]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[1] - rCurvature[1]);
        rStrainVector[2] = 0.5 * (rgab[2] - mInitialMetric.a_ab[2]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[2] - rCurvature[2]);
        // the other entries are (remain) zero (KL)
        
        KRATOS_CATCH("")
    }

    void IgaShell5pElement::CalculateStrainRM(
        Vector& rStrainVectorRM,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rg1,
        const Vector& rg2)
    {
        double thickness = GetProperties().GetValue(THICKNESS);
                
        rStrainVectorRM[0] = mZeta * thickness/2.0 * inner_prod(rDw_D1, rg1);
        rStrainVectorRM[1] = mZeta * thickness/2.0 * inner_prod(rDw_D2, rg2);
        rStrainVectorRM[2] = mZeta * thickness/2.0 * 0.5 * (inner_prod(rDw_D1, rg2) + inner_prod(rDw_D2, rg1));
        rStrainVectorRM[3] = 0.5 * inner_prod(rw, rg2);
        rStrainVectorRM[4] = 0.5 * inner_prod(rw, rg1);
    }

    void IgaShell5pElement::TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector& rCurvilinearStrain,
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

    void IgaShell5pElement::CalculateB(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const double thickness = GetProperties().GetValue(THICKNESS);

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size_KL = number_of_nodes * 3;
        const unsigned int mat_size = number_of_nodes * 5;

        // membrane part
        for (unsigned int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 5;
            int dirr = r % 5;

            Vector dE_curvilinear = ZeroVector(3);
            // "if" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            if (dirr == 0 || dirr == 1 || dirr == 2)
            {
                // strain corresponding to E11, E22, E12
                dE_curvilinear[0] = DN_De(kr, 0)*rMetric.a1(dirr);
                dE_curvilinear[1] = DN_De(kr, 1)*rMetric.a2(dirr);
                dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*rMetric.a2(dirr) + rMetric.a1(dirr)*DN_De(kr, 1));
            }
            // calculated with simplified Q (ML)
            rB(0, r) += mInitialMetric.Q(0, 0) * dE_curvilinear[0] + mInitialMetric.Q(0, 1)*dE_curvilinear[1] 
                + mInitialMetric.Q(0, 2)*dE_curvilinear[2] ;
            rB(1, r) += mInitialMetric.Q(1, 0) * dE_curvilinear[0] + mInitialMetric.Q(1, 1) * dE_curvilinear[1] 
                + mInitialMetric.Q(1, 2) * dE_curvilinear[2];
            rB(2, r) += mInitialMetric.Q(2, 0)*dE_curvilinear[0] + mInitialMetric.Q(2, 1)*dE_curvilinear[1] 
                + mInitialMetric.Q(2, 2)*dE_curvilinear[2];

            // all other entries of rB are (remain) zero
            
        }

        // curvature part
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

            for (unsigned int j = 0; j < 3; j++)
            {
                double g3dg3lg3 = (rMetric.a3_KL_tilde[0] * dg3(j, 0) + rMetric.a3_KL_tilde[1] * dg3(j, 1) + rMetric.a3_KL_tilde[2] * dg3(j, 2))*inddA3;

                dn(j, 0) = dg3(j, 0)*invdA - rMetric.a3_KL_tilde[0] * g3dg3lg3;
                dn(j, 1) = dg3(j, 1)*invdA - rMetric.a3_KL_tilde[1] * g3dg3lg3;
                dn(j, 2) = dg3(j, 2)*invdA - rMetric.a3_KL_tilde[2] * g3dg3lg3;
            }

            // b refers to curvilinear and rB to local cartesian coordinate system
            // "index" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            for (unsigned int j = 0; j < 3; j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    b(k, index_KL + j) = - mZeta * thickness / 2.0 * (DDN_DDe(i, k) * rMetric.a3_KL[j] + rMetric.H(0, k) * dn(j, 0) 
                    + rMetric.H(1, k) * dn(j, 1) + rMetric.H(2, k) * dn(j, 2));
                }
                rB(0, index + j) += mInitialMetric.Q(0, 0) * b(0, index_KL + j);
                rB(1, index + j) += mInitialMetric.Q(1, 0) * b(0, index_KL + j) + mInitialMetric.Q(1, 1) * b(1, index_KL + j)
                    + mInitialMetric.Q(1, 2) * b(2, index_KL + j);
                rB(2, index + j) += mInitialMetric.Q(2, 0) * b(0, index_KL + j) + mInitialMetric.Q(2, 2) * b(2, index_KL + j);
            }
        }
    }

    void IgaShell5pElement::CalculateSecondVariations(
        SecondVariations& rSecondVariations,
        const MetricVariables& rMetric)
    {
        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
            const double thickness = GetProperties().GetValue(THICKNESS);

            const int number_of_nodes = GetGeometry().size();
            const int mat_size_KL = number_of_nodes * 3;
            const int mat_size = number_of_nodes * 5;
           
            double lg3_3 = pow(rMetric.dA, 3);
            double lg3_5 = pow(rMetric.dA, 5);
            double inv_lg3 = 1 / rMetric.dA;
            double inv_lg3_3 = 1 / lg3_3;
            double inv_lg3_5 = 1 / lg3_5;

            SecondVariations second_variations_KL(mat_size_KL);
            Matrix S_dg3 = ZeroMatrix(3, mat_size_KL);
            Vector S_g3dg3 = ZeroVector(mat_size_KL);
            Vector S_g3dg3lg3_3 = ZeroVector(mat_size_KL);
            Matrix S_dn = ZeroMatrix(3, mat_size_KL);
            // first variation of strain and curvature w.r.t. dof
            for (int r = 0; r < mat_size_KL; r++)
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

                S_g3dg3[r] = rMetric.a3_KL_tilde[0] * S_dg3(0, r) + rMetric.a3_KL_tilde[1] * S_dg3(1, r) + rMetric.a3_KL_tilde[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.a3_KL_tilde[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.a3_KL_tilde[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.a3_KL_tilde[2] * S_g3dg3lg3_3[r];
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
                        second_variations_KL.B11(r, s) += mInitialMetric.Q(0, 0) * ddE_cu[0];
                        second_variations_KL.B22(r, s) += mInitialMetric.Q(1, 0) * ddE_cu[0] + mInitialMetric.Q(1, 1) * ddE_cu[1]
                            + mInitialMetric.Q(1, 2) * ddE_cu[2];
                        second_variations_KL.B12(r, s) += mInitialMetric.Q(2, 0) * ddE_cu[0] + mInitialMetric.Q(2, 2) * ddE_cu[2];
                        if (r != s){
                            second_variations_KL.B11(s, r) += second_variations_KL.B11(r, s);
                            second_variations_KL.B22(s, r) += second_variations_KL.B22(r, s);
                            second_variations_KL.B12(s, r) += second_variations_KL.B12(r, s);
                        }
                    }

                    // curvature
                    array_1d<double, 3> ddg3 = ZeroVector(3);
                    int dirt = 4 - dirr - dirs;
                    int ddir = dirr - dirs;
                    if (ddir == -1)      ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 2) ddg3(dirt - 1) = DN_De(kr, 0)*DN_De(ks, 1) - DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == 1) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);
                    else if (ddir == -2) ddg3(dirt - 1) = -DN_De(kr, 0)*DN_De(ks, 1) + DN_De(ks, 0)*DN_De(kr, 1);

                    double c = -(ddg3[0] * rMetric.a3_KL_tilde[0] + ddg3[1] * rMetric.a3_KL_tilde[1] + ddg3[2] * rMetric.a3_KL_tilde[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.a3_KL_tilde[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.a3_KL_tilde[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.a3_KL_tilde[2];

                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = mZeta * thickness / 2.0 * (DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2]);
                    ddK_cu[1] = mZeta * thickness / 2.0 * (DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2]);
                    ddK_cu[2] = mZeta * thickness / 2.0 * (DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2]);

                    // calculated with simplified Q (ML)
                    second_variations_KL.B11(r, s) += mInitialMetric.Q(0, 0) * ddK_cu[0];
                    second_variations_KL.B22(r, s) += mInitialMetric.Q(1, 0) * ddK_cu[0] + mInitialMetric.Q(1, 1) * ddK_cu[1] 
                        + mInitialMetric.Q(1, 2) * ddK_cu[2];
                    second_variations_KL.B12(r, s) += mInitialMetric.Q(2, 0) * ddK_cu[0] + mInitialMetric.Q(2, 2) * ddK_cu[2];
                    if (r != s){
                        second_variations_KL.B11(s, r) += second_variations_KL.B11(r, s);
                        second_variations_KL.B22(s, r) += second_variations_KL.B22(r, s);
                        second_variations_KL.B12(s, r) += second_variations_KL.B12(r, s);
                    }
                }
            }

            // transfer KL-second-variations to RM-second-variations
            for (unsigned int r = 0; r < mat_size; r++) {
                unsigned int kr = r / 5;
                unsigned int dirr = r % 5;
                unsigned int r_KL = kr * 3 + dirr;
                if (dirr != 3 || dirr != 4){
                    for (unsigned int s = 0; s < mat_size; s++){
                        unsigned int ks = s / 5;
                        unsigned int dirs = s % 5;
                        unsigned int s_KL = ks * 3 + dirs;
                        if (dirs != 3 || dirs != 4){
                            rSecondVariations.B11(r, s) += second_variations_KL.B11(r_KL, s_KL);
                            rSecondVariations.B22(r, s) += second_variations_KL.B22(r_KL, s_KL);
                            rSecondVariations.B12(r, s) += second_variations_KL.B12(r_KL, s_KL);               
                        }
                    }
                }
            }
        }
    }

    void IgaShell5pElement::CalculateVariationsRM(        
        Matrix& rB,
        SecondVariations& rSecondVariations,
        const Vector& rw,
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
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;
        
        // Vector testvector = ZeroVector(3); // ML
        // Vector testvector1 = ZeroVector(3); // ML
        // Matrix testmatrix = ZeroMatrix(mat_size, mat_size); // ML
        // double testdouble = 0.0; // ML
        // // ML: Anfang Oesterle kopiert
        // Matrix bop = ZeroMatrix(5, mat_size);
        // Matrix bopVergleich = ZeroMatrix(5, mat_size);
        // int count1 = 0;
        // int count2 = 0;
        // double dw_1_dr;
        // double dw_2_dr;
        // double dw_1dT1_dr;
        // double dw_2dT1_dr;
        // double dw_1dT2_dr;
        // double dw_2dT2_dr;
        // array_1d<double, 3> a1_dr;
        // array_1d<double, 3> a2_dr;
        // array_1d<double, 3> da1_dT1_dr;
        // array_1d<double, 3> da1_dT2_dr;
        // array_1d<double, 3> da2_dT1_dr;
        // array_1d<double, 3> da2_dT2_dr;

        // array_1d<double, 3> dw_dr;
        // array_1d<double, 3> dw_dT1_dr;
        // array_1d<double, 3> dw_dT2_dr;
        // array_1d<double, 3> Dg1_D1;
        // array_1d<double, 3> Dg2_D2;
        // array_1d<double, 3> Dg1_D2;
        // for (unsigned int i = 0; i < 3; i++)
        // {
        //     Dg1_D1[i] = rActualMetric.H(i, 0);
        //     Dg2_D2[i] = rActualMetric.H(i, 1);
        //     Dg1_D2[i] = rActualMetric.H(i, 2);
        // }
        // for (int j = 0; j < number_of_nodes; j++)
        // {
        //     for (int i = 0; i < 5; i++)
        //     {
        //     if (i == 0)
        //     {
        //         a1_dr[0] = DN_De(count2, 0);
        //         a1_dr[1] = 0.0;
        //         a1_dr[2] = 0.0;
        //         a2_dr[0] = DN_De(count2, 1);
        //         a2_dr[1] = 0.0;
        //         a2_dr[2] = 0.0;

        //         da1_dT1_dr[0] = DDN_DDe(count2, 0);
        //         da1_dT1_dr[1] = 0.0;
        //         da1_dT1_dr[2] = 0.0;
        //         da1_dT2_dr[0] = DDN_DDe(count2, 2);
        //         da1_dT2_dr[1] = 0.0;
        //         da1_dT2_dr[2] = 0.0;
        //         da2_dT1_dr[0] = DDN_DDe(count2, 2);
        //         da2_dT1_dr[1] = 0.0;
        //         da2_dT1_dr[2] = 0.0;
        //         da2_dT2_dr[0] = DDN_DDe(count2, 1);
        //         da2_dT2_dr[1] = 0.0;
        //         da2_dT2_dr[2] = 0.0;

        //         dw_1_dr = 0.0;
        //         dw_2_dr = 0.0;
        //         dw_1dT1_dr = 0.0;
        //         dw_2dT1_dr = 0.0;
        //         dw_1dT2_dr = 0.0;
        //         dw_2dT2_dr = 0.0;
        //     }
        //     else if (i == 1)
        //     {
        //         a1_dr[0] = 0;
        //         a1_dr[1] = DN_De(count2, 0);
        //         a1_dr[2] = 0;
        //         a2_dr[0] = 0;
        //         a2_dr[1] = DN_De(count2, 1);
        //         a2_dr[2] = 0;

        //         da1_dT1_dr[0] = 0.0;
        //         da1_dT1_dr[1] = DDN_DDe(count2, 0);
        //         da1_dT1_dr[2] = 0.0;
        //         da1_dT2_dr[0] = 0.0;
        //         da1_dT2_dr[1] = DDN_DDe(count2, 2);
        //         da1_dT2_dr[2] = 0.0;
        //         da2_dT1_dr[0] = 0.0;
        //         da2_dT1_dr[1] = DDN_DDe(count2, 2);
        //         da2_dT1_dr[2] = 0.0;
        //         da2_dT2_dr[0] = 0.0;
        //         da2_dT2_dr[1] = DDN_DDe(count2, 1);
        //         da2_dT2_dr[2] = 0.0;

        //         dw_1_dr = 0.0;
        //         dw_2_dr = 0.0;
        //         dw_1dT1_dr = 0.0;
        //         dw_2dT1_dr = 0.0;
        //         dw_1dT2_dr = 0.0;
        //         dw_2dT2_dr = 0.0;
        //     }
        //     else if (i == 2)
        //     {
        //         a1_dr[0] = 0;
        //         a1_dr[1] = 0;
        //         a1_dr[2] = DN_De(count2, 0);
        //         a2_dr[0] = 0;
        //         a2_dr[1] = 0;
        //         a2_dr[2] = DN_De(count2, 1);

        //         da1_dT1_dr[0] = 0.0;
        //         da1_dT1_dr[1] = 0.0;
        //         da1_dT1_dr[2] = DDN_DDe(count2, 0);
        //         da1_dT2_dr[0] = 0.0;
        //         da1_dT2_dr[1] = 0.0;
        //         da1_dT2_dr[2] = DDN_DDe(count2, 2);
        //         da2_dT1_dr[0] = 0.0;
        //         da2_dT1_dr[1] = 0.0;
        //         da2_dT1_dr[2] = DDN_DDe(count2, 2);
        //         da2_dT2_dr[0] = 0.0;
        //         da2_dT2_dr[1] = 0.0;
        //         da2_dT2_dr[2] = DDN_DDe(count2, 1);

        //         dw_1_dr = 0.0;
        //         dw_2_dr = 0.0;
        //         dw_1dT1_dr = 0.0;
        //         dw_2dT1_dr = 0.0;
        //         dw_1dT2_dr = 0.0;
        //         dw_2dT2_dr = 0.0;
        //     }
        //     else if (i == 3)
        //     {
        //         a1_dr[0] = 0;
        //         a1_dr[1] = 0;
        //         a1_dr[2] = 0;
        //         a2_dr[0] = 0;
        //         a2_dr[1] = 0;
        //         a2_dr[2] = 0;
        //         dw_1_dr = N(count2);
        //         dw_2_dr = 0.0;
        //         dw_1dT1_dr = DN_De(count2, 0);
        //         dw_2dT1_dr = 0.0;
        //         dw_1dT2_dr = DN_De(count2, 1);
        //         dw_2dT2_dr = 0.0;
        //     }
        //     else if (i == 4)
        //     {
        //         a1_dr[0] = 0;
        //         a1_dr[1] = 0;
        //         a1_dr[2] = 0;
        //         a2_dr[0] = 0;
        //         a2_dr[1] = 0;
        //         a2_dr[2] = 0;
        //         dw_1_dr = 0.0;
        //         dw_2_dr = N(count2);
        //         dw_1dT1_dr = 0.0;
        //         dw_2dT1_dr = DN_De(count2, 0);
        //         dw_1dT2_dr = 0.0;
        //         dw_2dT2_dr = DN_De(count2, 1);
        //     }

        //     /* Partielle Ableitung/Variation nach FHG (_dr entspricht ,r) */
        //     dw_dr = dw_1_dr*rActualMetric.a1 + rw_alpha[0]*a1_dr + dw_2_dr*rActualMetric.a2 + rw_alpha[1]*a2_dr;
        //     dw_dT1_dr = dw_1dT1_dr*rActualMetric.a1 + rDw_alpha_Dbeta(0, 0)*a1_dr + dw_1_dr*Dg1_D1 + rw_alpha[0]*da1_dT1_dr +
        //                 dw_2dT1_dr*rActualMetric.a2 + rDw_alpha_Dbeta(1, 0)*a2_dr + dw_2_dr*Dg1_D2 + rw_alpha[1]*da2_dT1_dr;
        //     dw_dT2_dr = dw_1dT2_dr*rActualMetric.a1 + rDw_alpha_Dbeta(0, 1)*a1_dr + dw_1_dr*Dg1_D2 + rw_alpha[0]*da1_dT2_dr +
        //                 dw_2dT2_dr*rActualMetric.a2 + rDw_alpha_Dbeta(1, 1)*a2_dr + dw_2_dr*Dg2_D2 + rw_alpha[1]*da2_dT2_dr;

        //     bop(0,count1) = bop(0,count1) + mZeta*thickness/2.0*(inner_prod(a1_dr, rDw_D1) + inner_prod(rActualMetric.a1, dw_dT1_dr));
        //     bop(2,count1) = bop(2,count1) + 0.5 * mZeta*thickness/2.0*(inner_prod(a1_dr,rDw_D2) + inner_prod(rActualMetric.a1,dw_dT2_dr) + inner_prod(a2_dr,rDw_D1) + inner_prod(rActualMetric.a2,dw_dT1_dr));
        //     bop(1,count1) = bop(1,count1) + mZeta*thickness/2.0*(inner_prod(a2_dr,rDw_D2) + inner_prod(rActualMetric.a2,dw_dT2_dr));
        //     bop(4,count1) = bop(4,count1) + 0.5 * (inner_prod(dw_dr,rActualMetric.a1) + inner_prod(rw,a1_dr));
        //     bop(3,count1) = bop(3,count1) + 0.5 * (inner_prod(dw_dr,rActualMetric.a2) + inner_prod(rw,a2_dr));

        //     count1 ++;
        //     }
        //     count2 ++;
        // }
        // rB += prod(mInitialMetric.Q, bop);
        // // if(Id() == 4){
        // //     KRATOS_WATCH(mInitialMetric.Q)
        // //     KRATOS_WATCH(bop)
        // //     KRATOS_WATCH(rB)
        // // }

        // //second variations Oesterle
        // Matrix dE11_dkl = ZeroMatrix(mat_size, mat_size);
        // Matrix dE12_dkl = ZeroMatrix(mat_size, mat_size);
        // Matrix dE22_dkl = ZeroMatrix(mat_size, mat_size);
        // Matrix dE13_dkl = ZeroMatrix(mat_size, mat_size);
        // Matrix dE23_dkl = ZeroMatrix(mat_size, mat_size);
        // Matrix dE11_dklVergleich = ZeroMatrix(mat_size, mat_size);
        // Matrix dE12_dklVergleich = ZeroMatrix(mat_size, mat_size);
        // Matrix dE22_dklVergleich = ZeroMatrix(mat_size, mat_size);
        // Matrix dE13_dklVergleich = ZeroMatrix(mat_size, mat_size);
        // Matrix dE23_dklVergleich = ZeroMatrix(mat_size, mat_size);
        // array_1d<double, 3> da1xa2k;
        // array_1d<double, 3> da1xa2l;
        // array_1d<double, 3> dda1xa2kl;
        // array_1d<double, 3> da1k;
        // array_1d<double, 3> da2k;
        // array_1d<double, 3> da3k;
        // array_1d<double, 3> da1l;
        // array_1d<double, 3> da2l;
        // array_1d<double, 3> da3l;
        // array_1d<double, 3> da11k;
        // array_1d<double, 3> da12k;
        // array_1d<double, 3> da21k;
        // array_1d<double, 3> da22k;
        // array_1d<double, 3> da11l;
        // array_1d<double, 3> da12l;
        // array_1d<double, 3> da21l;
        // array_1d<double, 3> da22l;
        // array_1d<double, 3> dda3kl;

        // array_1d<double, 3> dw_dk;
        // array_1d<double, 3> dw_dl;
        // array_1d<double, 3> ddw_dkl;
        // array_1d<double, 3> dw_dT1_dk;
        // array_1d<double, 3> dw_dT2_dk;
        // array_1d<double, 3> dw_dT1_dl;
        // array_1d<double, 3> dw_dT2_dl;
        // array_1d<double, 3> ddw_dT1_dkl;
        // array_1d<double, 3> ddw_dT2_dkl;
        // double dw1_dk, dw2_dk, dw1_dl, dw2_dl;
        // double dw1_dT1_dk, dw2_dT1_dk, dw1_dT2_dk, dw2_dT2_dk,
        //         dw1_dT1_dl, dw2_dT1_dl, dw1_dT2_dl, dw2_dT2_dl;
        // int i, j, k, l;
        // double da1xa2k_a1xa2;
        // double da1xa2k_a1xa2_hact;
        // double da1xa2l_a1xa2;
        // double da1xa2l_a1xa2_hact;
        
        // // array_1d<double, 3> Dg1_D1;
        // // array_1d<double, 3> Dg2_D2;
        // // array_1d<double, 3> Dg1_D2;
        // // for (unsigned int i = 0; i < 3; i++)
        // // {
        // //     Dg1_D1[i] = mInitialMetric.H(i, 0);
        // //     Dg2_D2[i] = mInitialMetric.H(i, 1);
        // //     Dg1_D2[i] = mInitialMetric.H(i, 2);
        // // }

        // for (k = 0; k < number_of_nodes; k++)
        // {
        //     for (i = 0; i < 5; i++)
        //     {
        //     for (l = 0; l < number_of_nodes; l++)
        //     {
        //         for (j = 0; j < 5; j++)
        //         {
        //         /* Anteile aus 3p-Formulierung (nicht veraendert) */
        //         if (i < 3 && j < 3)
        //         {
        //         if (i == 0)
        //         {
        //             da1k[0] = DN_De(k, 0);
        //             da1k[1] = 0;
        //             da1k[2] = 0;
        //             da2k[0] = DN_De(k, 1);
        //             da2k[1] = 0;
        //             da2k[2] = 0;
        //         }
        //         else if (i == 1)
        //         {
        //             da1k[0] = 0;
        //             da1k[1] = DN_De(k, 0);
        //             da1k[2] = 0;
        //             da2k[0] = 0;
        //             da2k[1] = DN_De(k, 1);
        //             da2k[2] = 0;
        //         }
        //         else
        //         {
        //             da1k[0] = 0;
        //             da1k[1] = 0;
        //             da1k[2] = DN_De(k, 0);
        //             da2k[0] = 0;
        //             da2k[1] = 0;
        //             da2k[2] = DN_De(k, 1);
        //         }

        //         array_1d<double, 3> da1kxa2, a1xda2k;
        //         MathUtils<double>::CrossProduct(da1kxa2, da1k, rActualMetric.a2);
        //         MathUtils<double>::CrossProduct(a1xda2k, rActualMetric.a1, da2k);
        //         da1xa2k = da1kxa2 + a1xda2k;
        //         da1xa2k_a1xa2 = (da1xa2k[0] * rActualMetric.a3_KL_tilde[0] + da1xa2k[1] * rActualMetric.a3_KL_tilde[1] + da1xa2k[2] * rActualMetric.a3_KL_tilde[2]);
        //         da1xa2k_a1xa2_hact = da1xa2k_a1xa2 / (rActualMetric.dA * rActualMetric.dA * rActualMetric.dA);
        //         da3k[0] = da1xa2k[0] / rActualMetric.dA - rActualMetric.a3_KL_tilde[0] * da1xa2k_a1xa2_hact;
        //         da3k[1] = da1xa2k[1] / rActualMetric.dA - rActualMetric.a3_KL_tilde[1] * da1xa2k_a1xa2_hact;
        //         da3k[2] = da1xa2k[2] / rActualMetric.dA - rActualMetric.a3_KL_tilde[2] * da1xa2k_a1xa2_hact;

        //         if (j == 0)
        //         {
        //             da1l[0] = DN_De(l, 0);
        //             da1l[1] = 0;
        //             da1l[2] = 0;
        //             da2l[0] = DN_De(l, 1);
        //             da2l[1] = 0;
        //             da2l[2] = 0;
        //         }
        //         else if (j == 1)
        //         {
        //             da1l[0] = 0;
        //             da1l[1] = DN_De(l, 0);
        //             da1l[2] = 0;
        //             da2l[0] = 0;
        //             da2l[1] = DN_De(l, 1);
        //             da2l[2] = 0;
        //         }
        //         else
        //         {
        //             da1l[0] = 0;
        //             da1l[1] = 0;
        //             da1l[2] = DN_De(l, 0);
        //             da2l[0] = 0;
        //             da2l[1] = 0;
        //             da2l[2] = DN_De(l, 1);
        //         }

        //         array_1d<double, 3> da1lxa2, a1xda2l;
        //         MathUtils<double>::CrossProduct(da1lxa2, da1l, rActualMetric.a2);
        //         MathUtils<double>::CrossProduct(a1xda2l, rActualMetric.a1, da2l);
        //         da1xa2l = da1lxa2 + a1xda2l;
        //         da1xa2l_a1xa2 = (da1xa2l[0] * rActualMetric.a3_KL_tilde[0] + da1xa2l[1] * rActualMetric.a3_KL_tilde[1] + da1xa2l[2] * rActualMetric.a3_KL_tilde[2]);
        //         da1xa2l_a1xa2_hact = da1xa2l_a1xa2 / (rActualMetric.dA * rActualMetric.dA * rActualMetric.dA);
        //         da3l[0] = da1xa2l[0] / rActualMetric.dA - rActualMetric.a3_KL_tilde[0] * da1xa2l_a1xa2_hact;
        //         da3l[1] = da1xa2l[1] / rActualMetric.dA - rActualMetric.a3_KL_tilde[1] * da1xa2l_a1xa2_hact;
        //         da3l[2] = da1xa2l[2] / rActualMetric.dA - rActualMetric.a3_KL_tilde[2] * da1xa2l_a1xa2_hact;

        //         if (i == 0)
        //         {
        //             if (j == 0)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = 0;
        //             }
        //             else if (j == 1)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = DN_De(k, 0) * DN_De(l, 1) - DN_De(l, 0) * DN_De(k, 1);
        //             }
        //             else if (j == 2)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = -DN_De(k, 0) * DN_De(l, 1) + DN_De(l, 0) * DN_De(k, 1);
        //             dda1xa2kl[2] = 0;
        //             }
        //         }

        //         else if (i == 1)
        //         {
        //             if (j == 0)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = -DN_De(k, 0) * DN_De(l, 1) + DN_De(l, 0) * DN_De(k, 1);
        //             }
        //             else if (j == 1)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = 0;
        //             }
        //             else if (j == 2)
        //             {
        //             dda1xa2kl[0] = DN_De(k, 0) * DN_De(l, 1) - DN_De(l, 0) * DN_De(k, 1);
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = 0;
        //             }
        //         }

        //         else
        //         {
        //             if (j == 0)
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = DN_De(k, 0) * DN_De(l, 1) - DN_De(l, 0) * DN_De(k, 1);
        //             dda1xa2kl[2] = 0;
        //             }
        //             else if (j == 1)
        //             {
        //             dda1xa2kl[0] = -DN_De(k, 0) * DN_De(l, 1) + DN_De(l, 0) * DN_De(k, 1);
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = 0;
        //             }
        //             else
        //             {
        //             dda1xa2kl[0] = 0;
        //             dda1xa2kl[1] = 0;
        //             dda1xa2kl[2] = 0;
        //             }
        //         }

        //         // C = -(dda1xa2kl[0] * rActualMetric.a3_KL_tilde[0] + dda1xa2kl[1] * rActualMetric.a3_KL_tilde[1]
        //         //     + dda1xa2kl[2] * rActualMetric.a3_KL_tilde[2] + da1xa2k[0] * da1xa2l[0]
        //         //     + da1xa2k[1] * da1xa2l[1] + da1xa2k[2] * da1xa2l[2])
        //         //     / (rActualMetric.dA * rActualMetric.dA * rActualMetric.dA);

        //         // D = 3.0 * (da1xa2k_a1xa2) * (da1xa2l_a1xa2) / (rActualMetric.dA * rActualMetric.dA * rActualMetric.dA * rActualMetric.dA * rActualMetric.dA);

        //         // dda3kl[0] = dda1xa2kl[0] / rActualMetric.dA - da1xa2l_a1xa2_hact * da1xa2k[0]
        //         //             - da1xa2k_a1xa2_hact * da1xa2l[0] + C * rActualMetric.a3_KL_tilde[0] + D * rActualMetric.a3_KL_tilde[0];
        //         // dda3kl[1] = dda1xa2kl[1] / rActualMetric.dA - da1xa2l_a1xa2_hact * da1xa2k[1]
        //         //             - da1xa2k_a1xa2_hact * da1xa2l[1] + C * rActualMetric.a3_KL_tilde[1] + D * rActualMetric.a3_KL_tilde[1];
        //         // dda3kl[2] = dda1xa2kl[2] / rActualMetric.dA - da1xa2l_a1xa2_hact * da1xa2k[2]
        //         //             - da1xa2k_a1xa2_hact * da1xa2l[2] + C * rActualMetric.a3_KL_tilde[2] + D * rActualMetric.a3_KL_tilde[2];

        //         // /* Membrananteil */
        //         // if (i == j)
        //         // {
        //         //     dE11_dkl[5 * k + i][5 * l + j] = DN_De(k, 0) * DN_De(l, 0);
        //         //     dE12_dkl[5 * k + i][5 * l + j] = 0.5 * (DN_De(k, 0) * DN_De(l, 1) + DN_De(l, 0) * DN_De(k, 1));
        //         //     dE22_dkl[5 * k + i][5 * l + j] = DN_De(k, 1) * DN_De(l, 1);
        //         // }
        //         // else
        //         // {
        //         //     dE11_dkl[5 * k + i][5 * l + j] = 0;
        //         //     dE12_dkl[5 * k + i][5 * l + j] = 0;
        //         //     dE22_dkl[5 * k + i][5 * l + j] = 0;
        //         // }

        //         // /* Kruemmungsanteil */
        //         // dE11_dkl[5 * k + i][5 * l + j] = dE11_dkl[5 * k + i][5 * l + j]
        //         //     - mZeta * thickness / 2.0
        //         //         * (DDN_DDe(k, 0) * da3l[i] + DDN_DDe(l, 0) * da3k[j]
        //         //             + Dg1_D1[0] * dda3kl[0] + Dg1_D1[1] * dda3kl[1]
        //         //             + Dg1_D1[2] * dda3kl[2]);
        //         // dE12_dkl[5 * k + i][5 * l + j] = dE12_dkl[5 * k + i][5 * l + j]
        //         //     - mZeta * thickness / 2.0
        //         //         * (DDN_DDe(k, 2) * da3l[i] + DDN_DDe(l, 2) * da3k[j]
        //         //             + Dg1_D2[0] * dda3kl[0] + Dg1_D2[1] * dda3kl[1]
        //         //             + Dg1_D2[2] * dda3kl[2]);
        //         // dE22_dkl[5 * k + i][5 * l + j] = dE22_dkl[5 * k + i][5 * l + j]
        //         //     - mZeta * thickness / 2.0
        //         //         * (DDN_DDe(k, 1) * da3l[i] + DDN_DDe(l, 1) * da3k[j]
        //         //             + Dg2_D2[0] * dda3kl[0] + Dg2_D2[1] * dda3kl[1]
        //         //             + Dg2_D2[2] * dda3kl[2]);
        //         }

        //         else
        //         {
        //         /* Anteile aus 5p-Formulierung  */

        //         /* initialisieren */
        //         dw1_dk = 0.0;
        //         dw1_dT1_dk = 0.0;
        //         dw1_dT2_dk = 0.0;

        //         dw2_dk = 0.0;
        //         dw2_dT1_dk = 0.0;
        //         dw2_dT2_dk = 0.0;

        //         /* Ableitung nach dem ersten FHG k */
        //         if (i == 0)
        //         {
        //             da1k[0] = DN_De(k, 0);
        //             da1k[1] = 0.0;
        //             da1k[2] = 0.0;
        //             da2k[0] = DN_De(k, 1);
        //             da2k[1] = 0.0;
        //             da2k[2] = 0.0;
        //             da11k[0] = DDN_DDe(k, 0);
        //             da11k[1] = 0.0;
        //             da11k[2] = 0.0;
        //             da12k[0] = DDN_DDe(k, 2);
        //             da12k[1] = 0.0;
        //             da12k[2] = 0.0;
        //             da22k[0] = DDN_DDe(k, 1);
        //             da22k[1] = 0.0;
        //             da22k[2] = 0.0;
        //             da21k = da12k;

        //             dw1_dk = 0.0;
        //             dw1_dT1_dk = 0.0;
        //             dw1_dT2_dk = 0.0;

        //             dw2_dk = 0.0;
        //             dw2_dT1_dk = 0.0;
        //             dw2_dT2_dk = 0.0;
        //         }
        //         else if (i == 1)
        //         {
        //             da1k[0] = 0.0;
        //             da1k[1] = DN_De(k, 0);
        //             da1k[2] = 0.0;
        //             da2k[0] = 0.0;
        //             da2k[1] = DN_De(k, 1);
        //             da2k[2] = 0.0;
        //             da11k[0] = 0.0;
        //             da11k[1] = DDN_DDe(k, 0);
        //             da11k[2] = 0.0;
        //             da12k[0] = 0.0;
        //             da12k[1] = DDN_DDe(k, 2);
        //             da12k[2] = 0.0;
        //             da22k[0] = 0.0;
        //             da22k[1] = DDN_DDe(k, 1);
        //             da22k[2] = 0.0;
        //             da21k = da12k;

        //             dw1_dk = 0.0;
        //             dw1_dT1_dk = 0.0;
        //             dw1_dT2_dk = 0.0;

        //             dw2_dk = 0.0;
        //             dw2_dT1_dk = 0.0;
        //             dw2_dT2_dk = 0.0;
        //         }
        //         else if (i == 2)
        //         {
        //             da1k[0] = 0.0;
        //             da1k[1] = 0.0;
        //             da1k[2] = DN_De(k, 0);
        //             da2k[0] = 0.0;
        //             da2k[1] = 0.0;
        //             da2k[2] = DN_De(k, 1);
        //             da11k[0] = 0.0;
        //             da11k[1] = 0.0;
        //             da11k[2] = DDN_DDe(k, 0);
        //             da12k[0] = 0.0;
        //             da12k[1] = 0.0;
        //             da12k[2] = DDN_DDe(k, 2);
        //             da22k[0] = 0.0;
        //             da22k[1] = 0.0;
        //             da22k[2] = DDN_DDe(k, 1);
        //             da21k = da12k;

        //             dw1_dk = 0.0;
        //             dw1_dT1_dk = 0.0;
        //             dw1_dT2_dk = 0.0;

        //             dw2_dk = 0.0;
        //             dw2_dT1_dk = 0.0;
        //             dw2_dT2_dk = 0.0;
        //         }
        //         else if (i == 3)
        //         {
        //             da1k[0] = 0.0;
        //             da1k[1] = 0.0;
        //             da1k[2] = 0.0;
        //             da2k[0] = 0.0;
        //             da2k[1] = 0.0;
        //             da2k[2] = 0.0;
        //             da11k[0] = 0.0;
        //             da11k[1] = 0.0;
        //             da11k[2] = 0.0;
        //             da12k[0] = 0.0;
        //             da12k[1] = 0.0;
        //             da12k[2] = 0.0;
        //             da22k[0] = 0.0;
        //             da22k[1] = 0.0;
        //             da22k[2] = 0.0;
        //             da21k = da12k;

        //             dw1_dk = N(k);
        //             dw1_dT1_dk = DN_De(k, 0);
        //             dw1_dT2_dk = DN_De(k, 1);

        //             dw2_dk = 0.0;
        //             dw2_dT1_dk = 0.0;
        //             dw2_dT2_dk = 0.0;
        //         }
        //         else if (i == 4)
        //         {
        //             da1k[0] = 0.0;
        //             da1k[1] = 0.0;
        //             da1k[2] = 0.0;
        //             da2k[0] = 0.0;
        //             da2k[1] = 0.0;
        //             da2k[2] = 0.0;
        //             da11k[0] = 0.0;
        //             da11k[1] = 0.0;
        //             da11k[2] = 0.0;
        //             da12k[0] = 0.0;
        //             da12k[1] = 0.0;
        //             da12k[2] = 0.0;
        //             da22k[0] = 0.0;
        //             da22k[1] = 0.0;
        //             da22k[2] = 0.0;
        //             da21k = da12k;

        //             dw1_dk = 0.0;
        //             dw1_dT1_dk = 0.0;
        //             dw1_dT2_dk = 0.0;

        //             dw2_dk = N(k);
        //             dw2_dT1_dk = DN_De(k, 0);
        //             dw2_dT2_dk = DN_De(k, 1);
        //         }
                

        //         dw_dk = dw1_dk*rActualMetric.a1 + rw_alpha(0)*da1k + dw2_dk*rActualMetric.a2 + rw_alpha(1)*da2k;
        //         dw_dT1_dk = dw1_dT1_dk*rActualMetric.a1 + rDw_alpha_Dbeta(0, 0)*da1k + dw1_dk*Dg1_D1 + rw_alpha(0)*da11k
        //                     + dw2_dT1_dk*rActualMetric.a2 + rDw_alpha_Dbeta(1, 0)*da2k + dw2_dk*Dg1_D2 + rw_alpha(1)*da21k;
        //         dw_dT2_dk = dw1_dT2_dk*rActualMetric.a1 + rDw_alpha_Dbeta(0, 1)*da1k + dw1_dk*Dg1_D2 + rw_alpha(0)*da12k
        //                     + dw2_dT2_dk*rActualMetric.a2 + rDw_alpha_Dbeta(1, 1)*da2k + dw2_dk*Dg2_D2 + rw_alpha(1)*da22k;

        //         /* initialisieren */
        //         dw1_dl = 0.0;
        //         dw1_dT1_dl = 0.0;
        //         dw1_dT2_dl = 0.0;

        //         dw2_dl = 0.0;
        //         dw2_dT1_dl = 0.0;
        //         dw2_dT2_dl = 0.0;

        //         /* Ableitung nach dem zweiten FHG l */
        //         if (j == 0)
        //         {
        //             da1l[0] = DN_De(l, 0);
        //             da1l[1] = 0.0;
        //             da1l[2] = 0.0;
        //             da2l[0] = DN_De(l, 1);
        //             da2l[1] = 0.0;
        //             da2l[2] = 0.0;
        //             da11l[0] = DDN_DDe(l, 0);
        //             da11l[1] = 0.0;
        //             da11l[2] = 0.0;
        //             da12l[0] = DDN_DDe(l, 2);
        //             da12l[1] = 0.0;
        //             da12l[2] = 0.0;
        //             da22l[0] = DDN_DDe(l, 1);
        //             da22l[1] = 0.0;
        //             da22l[2] = 0.0;
        //             da21l = da12l;

        //             dw1_dl = 0.0;
        //             dw1_dT1_dl = 0.0;
        //             dw1_dT2_dl = 0.0;

        //             dw2_dl = 0.0;
        //             dw2_dT1_dl = 0.0;
        //             dw2_dT2_dl = 0.0;
        //         }
        //         else if (j == 1)
        //         {
        //             da1l[0] = 0.0;
        //             da1l[1] = DN_De(l, 0);
        //             da1l[2] = 0.0;
        //             da2l[0] = 0.0;
        //             da2l[1] = DN_De(l, 1);
        //             da2l[2] = 0.0;
        //             da11l[0] = 0.0;
        //             da11l[1] = DDN_DDe(l, 0);
        //             da11l[2] = 0.0;
        //             da12l[0] = 0.0;
        //             da12l[1] = DDN_DDe(l, 2);
        //             da12l[2] = 0.0;
        //             da22l[0] = 0.0;
        //             da22l[1] = DDN_DDe(l, 1);
        //             da22l[2] = 0.0;
        //             da21l = da12l;

        //             dw1_dl = 0.0;
        //             dw1_dT1_dl = 0.0;
        //             dw1_dT2_dl = 0.0;

        //             dw2_dl = 0.0;
        //             dw2_dT1_dl = 0.0;
        //             dw2_dT2_dl = 0.0;
        //         }
        //         else if (j == 2)
        //         {
        //             da1l[0] = 0.0;
        //             da1l[1] = 0.0;
        //             da1l[2] = DN_De(l, 0);
        //             da2l[0] = 0.0;
        //             da2l[1] = 0.0;
        //             da2l[2] = DN_De(l, 1);
        //             da11l[0] = 0.0;
        //             da11l[1] = 0.0;
        //             da11l[2] = DDN_DDe(l, 0);
        //             da12l[0] = 0.0;
        //             da12l[1] = 0.0;
        //             da12l[2] = DDN_DDe(l, 2);
        //             da22l[0] = 0.0;
        //             da22l[1] = 0.0;
        //             da22l[2] = DDN_DDe(l, 1);
        //             da21l = da12l;

        //             dw1_dl = 0.0;
        //             dw1_dT1_dl = 0.0;
        //             dw1_dT2_dl = 0.0;

        //             dw2_dl = 0.0;
        //             dw2_dT1_dl = 0.0;
        //             dw2_dT2_dl = 0.0;
        //         }
        //         else if (j == 3)
        //         {
        //             da1l[0] = 0.0;
        //             da1l[1] = 0.0;
        //             da1l[2] = 0.0;
        //             da2l[0] = 0.0;
        //             da2l[1] = 0.0;
        //             da2l[2] = 0.0;
        //             da11l[0] = 0.0;
        //             da11l[1] = 0.0;
        //             da11l[2] = 0.0;
        //             da12l[0] = 0.0;
        //             da12l[1] = 0.0;
        //             da12l[2] = 0.0;
        //             da22l[0] = 0.0;
        //             da22l[1] = 0.0;
        //             da22l[2] = 0.0;
        //             da21l = da12l;

        //             dw1_dl = N(l);
        //             dw1_dT1_dl = DN_De(l, 0);
        //             dw1_dT2_dl = DN_De(l, 1);

        //             dw2_dl = 0.0;
        //             dw2_dT1_dl = 0.0;
        //             dw2_dT2_dl = 0.0;
        //         }
        //         else if (j == 4)
        //         {
        //             da1l[0] = 0.0;
        //             da1l[1] = 0.0;
        //             da1l[2] = 0.0;
        //             da2l[0] = 0.0;
        //             da2l[1] = 0.0;
        //             da2l[2] = 0.0;
        //             da11l[0] = 0.0;
        //             da11l[1] = 0.0;
        //             da11l[2] = 0.0;
        //             da12l[0] = 0.0;
        //             da12l[1] = 0.0;
        //             da12l[2] = 0.0;
        //             da22l[0] = 0.0;
        //             da22l[1] = 0.0;
        //             da22l[2] = 0.0;
        //             da21l = da12l;

        //             dw1_dl = 0.0;
        //             dw1_dT1_dl = 0.0;
        //             dw1_dT2_dl = 0.0;

        //             dw2_dl = N(l);
        //             dw2_dT1_dl = DN_De(l, 0);
        //             dw2_dT2_dl = DN_De(l, 1);
        //         }

        //         dw_dl = dw1_dl*rActualMetric.a1 + rw_alpha(0)*da1l + dw2_dl*rActualMetric.a2 + rw_alpha(1)*da2l;
        //         dw_dT1_dl = dw1_dT1_dl*rActualMetric.a1 + rDw_alpha_Dbeta(0, 0)*da1l + dw1_dl*Dg1_D1 + rw_alpha(0)*da11l
        //                     + dw2_dT1_dl*rActualMetric.a2 + rDw_alpha_Dbeta(1, 0)*da2l + dw2_dl*Dg1_D2 + rw_alpha(1)*da21l ;
        //         dw_dT2_dl = dw1_dT2_dl*rActualMetric.a1 + rDw_alpha_Dbeta(0, 1)*da1l + dw1_dl*Dg1_D2 + rw_alpha(0)*da12l
        //                     + dw2_dT2_dl*rActualMetric.a2 + rDw_alpha_Dbeta(1, 1)*da2l + dw2_dl*Dg2_D2 + rw_alpha(1)*da22l ;


        //         /* Ableitung nach beiden FHG kl */
        //         ddw_dkl = dw1_dk*da1l + dw1_dl*da1k + dw2_dk*da2l + dw2_dl*da2k;
        //         ddw_dT1_dkl = dw1_dT1_dk*da1l + dw1_dT1_dl*da1k + dw1_dk*da11l + dw1_dl*da11k
        //                     + dw2_dT1_dk*da2l + dw2_dT1_dl*da2k + dw2_dk*da21l + dw2_dl*da21k;
        //         ddw_dT2_dkl = dw1_dT2_dk*da1l + dw1_dT2_dl*da1k + dw1_dk*da12l + dw1_dl*da12k
        //                     + dw2_dT2_dk*da2l + dw2_dT2_dl*da2k + dw2_dk*da22l + dw2_dl*da22k;

        //         /* Anteile fuer Querschubverzerrungen (nur konstant in Dickenrichtung) */
        //         dE13_dkl(5 * k + i, 5 * l + j) = dE13_dkl(5 * k + i,5 * l + j)
        //                                         + 0.5*(inner_prod(da1k, dw_dl) + inner_prod(da1l, dw_dk) + inner_prod(rActualMetric.a1, ddw_dkl));
        //         dE23_dkl(5 * k + i, 5 * l + j) = dE23_dkl(5 * k + i,5 * l + j)
        //                                         + 0.5*(inner_prod(da2k, dw_dl) + inner_prod(da2l, dw_dk) + inner_prod(rActualMetric.a2, ddw_dkl));
        //         /* Anteile fuer Inplane-Verzerrungen (nur linear in Dickenrichtung theta^3) */
        //         if (Id() == 4){
        //             if (k==0 && i==4 && l==0 && j==1){
        //                 KRATOS_WATCH(dE11_dkl(5 * k + i, 5 * l + j))
        //             }
        //         }
        //         // dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i,5 * l + j)
        //         //                                 + mZeta*thickness/2.0*0.5*(inner_prod(da1k,dw_dT1_dl) + inner_prod(da1l,dw_dT1_dk) + inner_prod(rActualMetric.a1,ddw_dT1_dkl)
        //         //                                         + inner_prod(da1k,dw_dT1_dl) + inner_prod(da1l,dw_dT1_dk) + inner_prod(rActualMetric.a1,ddw_dT1_dkl));
        //         // dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i,5 * l + j)
        //         //                                 + mZeta*thickness/2.0*inner_prod(da1k,dw_dT1_dl);
        //         dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i,5 * l + j)
        //                                         + mZeta*thickness/2.0*inner_prod(da1l,dw_dT1_dk);
        //         // dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i,5 * l + j)
        //         //                                 + mZeta*thickness /2.0 * 0.5 * (inner_prod(rActualMetric.a1,ddw_dT1_dkl)
        //         //                                          + inner_prod(rActualMetric.a1,ddw_dT1_dkl));
        //         dE12_dkl(5 * k + i, 5 * l + j) = dE12_dkl(5 * k + i,5 * l + j)
        //                                         + mZeta*thickness/2.0*0.5*(inner_prod(da1k,dw_dT2_dl) + inner_prod(da1l,dw_dT2_dk) + inner_prod(rActualMetric.a1,ddw_dT2_dkl)
        //                                                 + inner_prod(da2k,dw_dT1_dl) + inner_prod(da2l,dw_dT1_dk) + inner_prod(rActualMetric.a2,ddw_dT1_dkl));
        //         dE22_dkl(5 * k + i, 5 * l + j) = dE22_dkl(5 * k + i,5 * l + j)
        //                                         + mZeta*thickness/2.0*0.5*(inner_prod(da2k,dw_dT2_dl) + inner_prod(da2l,dw_dT2_dk) + inner_prod(rActualMetric.a2,ddw_dT2_dkl)
        //                                                 + inner_prod(da2k,dw_dT2_dl) + inner_prod(da2l,dw_dT2_dk) + inner_prod(rActualMetric.a2,ddw_dT2_dkl));
        //         rSecondVariations.B11(5*k+i, 5*l+j) = mInitialMetric.Q(0, 0) * dE11_dkl(5*k+i, 5*l+j);
        //         rSecondVariations.B22(5*k+i, 5*l+j) = mInitialMetric.Q(1, 0) * dE11_dkl(5*k+i, 5*l+j) + mInitialMetric.Q(1, 1) * dE22_dkl(5*k+i, 5*l+j)
        //             + mInitialMetric.Q(1, 2) * dE12_dkl(5*k+i, 5*l+j);
        //         rSecondVariations.B12(5*k+i, 5*l+j) = mInitialMetric.Q(2, 0) * dE11_dkl(5*k+i, 5*l+j) + mInitialMetric.Q(2, 2) * dE12_dkl(5*k+i, 5*l+j);
        //         rSecondVariations.B23(5*k+i, 5*l+j) = mInitialMetric.Q(3, 3) * dE23_dkl(5*k+i, 5*l+j) + mInitialMetric.Q(3, 4) * dE13_dkl(5*k+i, 5*l+j);
        //         rSecondVariations.B13(5*k+i, 5*l+j) = mInitialMetric.Q(4, 4) * dE13_dkl(5*k+i, 5*l+j);
        //         // if (Id() == 4){
        //         //     testdouble = dE11_dkl(5*k+i, 5*l+j);
        //         //     if ((testdouble-dE11_dkl(5*k+i, 5*l+j)) != 0)
        //         //         KRATOS_WATCH(testdouble-dE11_dkl(5*k+i, 5*l+j))
        //         //     if (k==0 && i==4 && l==0 && j==1){
        //         //         KRATOS_WATCH(da1k)
        //         //         KRATOS_WATCH(dw_dT1_dl)
        //         //         KRATOS_WATCH(2*inner_prod(da1k, dw_dT1_dl))
        //         //         KRATOS_WATCH(da1l)
        //         //         KRATOS_WATCH(dw_dT1_dk)
        //         //         KRATOS_WATCH(inner_prod(da1l, dw_dT1_dk))
        //         //         KRATOS_WATCH(2*inner_prod(da1l, dw_dT1_dk))
        //         //         KRATOS_WATCH(dE11_dkl(5 * k + i, 5 * l + j))
        //         //         KRATOS_WATCH(ddw_dT1_dkl)
        //         //         testvector = dw_dT1_dk;
        //         //         testvector1 = da1l;
        //         //     }
        //         // }
                
        //         }
        //         }
        //     }
        //     }
        // }
        // // ML: Ende kopiert Oesterle

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
                dE_cur[0] = 0.5 * (rw(dirr) * DN_De(kr, 1));
                dE_cur[1] = 0.5 * (rw(dirr) * DN_De(kr, 0));
            }
            else if(dirr == 3){
                Dw_Dr(0, r) = N(kr) * rActualMetric.a1(0);
                Dw_Dr(1, r) = N(kr) * rActualMetric.a1(1);
                Dw_Dr(2, r) = N(kr) * rActualMetric.a1(2);
            }
            else {
                Dw_Dr(0, r) = N(kr) * rActualMetric.a2(0);
                Dw_Dr(1, r) = N(kr) * rActualMetric.a2(1);
                Dw_Dr(2, r) = N(kr) * rActualMetric.a2(2);                
            }
            dE_cur[0] += 0.5 * (Dw_Dr(0, r) * rActualMetric.a2(0) + Dw_Dr(1, r) * rActualMetric.a2(1) + Dw_Dr(2, r) * rActualMetric.a2(2));
            dE_cur[1] += 0.5 * (Dw_Dr(0, r) * rActualMetric.a1(0) + Dw_Dr(1, r) * rActualMetric.a1(1) + Dw_Dr(2, r) * rActualMetric.a1(2));

            // calculated with the simplified Q (ML)
            rB(3, r) += mInitialMetric.Q(3, 3) * dE_cur[0] + mInitialMetric.Q(3, 4) * dE_cur[1];
            rB(4, r) += mInitialMetric.Q(4, 4) * dE_cur[1];
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
                DDw_DD1r += DN_De(kr, 0) * rActualMetric.a1;
                DDw_DD1r[0] += N(kr) * rActualMetric.H(0, 0);
                DDw_DD1r[1] += N(kr) * rActualMetric.H(1, 0);
                DDw_DD1r[2] += N(kr) * rActualMetric.H(2, 0);
                DDw_DD2r += DN_De(kr, 1) * rActualMetric.a1;
                DDw_DD2r[0] += N(kr) * rActualMetric.H(0, 2);
                DDw_DD2r[1] += N(kr) * rActualMetric.H(1, 2);
                DDw_DD2r[2] += N(kr) * rActualMetric.H(2, 2);                
            }
            else if (dirr == 4){
                DDw_DD1r += DN_De(kr, 0) * rActualMetric.a2;
                DDw_DD1r[0] += N(kr) * rActualMetric.H(0, 2);
                DDw_DD1r[1] += N(kr) * rActualMetric.H(1, 2);
                DDw_DD1r[2] += N(kr) * rActualMetric.H(2, 2);
                DDw_DD2r += DN_De(kr, 1) * rActualMetric.a2;
                DDw_DD2r[0] += N(kr) * rActualMetric.H(0, 1);
                DDw_DD2r[1] += N(kr) * rActualMetric.H(1, 1);
                DDw_DD2r[2] += N(kr) * rActualMetric.H(2, 1);
            }
            dK_cu[0] += inner_prod(DDw_DD1r, rActualMetric.a1);
            dK_cu[1] += inner_prod(DDw_DD2r, rActualMetric.a2);
            dK_cu[2] += 0.5 * (inner_prod(DDw_DD1r, rActualMetric.a2) + inner_prod(DDw_DD2r, rActualMetric.a1));

            // calculated with simplified Q (ML)
            rB(0, r) += mZeta * thickness / 2.0 * mInitialMetric.Q(0, 0) * dK_cu[0];
            rB(1, r) += mZeta * thickness / 2.0 * (mInitialMetric.Q(1, 0) * dK_cu[0] + mInitialMetric.Q(1, 1) * dK_cu[1] 
                + mInitialMetric.Q(1, 2) * dK_cu[2]);
            rB(2, r) += mZeta * thickness / 2.0 * (mInitialMetric.Q(2, 0) * dK_cu[0] + mInitialMetric.Q(2, 2) * dK_cu[2]);
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
                    ddE_cu[0] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.a2);
                    ddE_cu[1] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.a1);

                    // calculated with simplified Q (ML)
                    rSecondVariations.B23(r, s) += mInitialMetric.Q(3, 3) * ddE_cu[0] + mInitialMetric.Q(3, 4) * ddE_cu[1];
                    rSecondVariations.B13(r, s) += mInitialMetric.Q(4, 4) * ddE_cu[1];
                    if (r != s){
                        rSecondVariations.B23(s, r) += rSecondVariations.B23(r, s);
                        rSecondVariations.B13(s, r) += rSecondVariations.B13(r, s);
                    }
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
                        DDw_DD1s += DN_De(ks, 0) * rActualMetric.a1;
                        DDw_DD1s[0] += N(ks) * rActualMetric.H(0, 0);
                        DDw_DD1s[1] += N(ks) * rActualMetric.H(1, 0);
                        DDw_DD1s[2] += N(ks) * rActualMetric.H(2, 0);
                        DDw_DD2s += DN_De(ks, 1) * rActualMetric.a1;
                        DDw_DD2s[0] += N(ks) * rActualMetric.H(0, 2);
                        DDw_DD2s[1] += N(ks) * rActualMetric.H(1, 2);
                        DDw_DD2s[2] += N(ks) * rActualMetric.H(2, 2);                
                    }
                    else if (dirs == 4){
                        DDw_DD1s += DN_De(ks, 0) * rActualMetric.a2;
                        DDw_DD1s[0] += N(ks) * rActualMetric.H(0, 2);
                        DDw_DD1s[1] += N(ks) * rActualMetric.H(1, 2);
                        DDw_DD1s[2] += N(ks) * rActualMetric.H(2, 2);
                        DDw_DD2s += DN_De(ks, 1) * rActualMetric.a2;
                        DDw_DD2s[0] += N(ks) * rActualMetric.H(0, 1);
                        DDw_DD2s[1] += N(ks) * rActualMetric.H(1, 1);
                        DDw_DD2s[2] += N(ks) * rActualMetric.H(2, 1);
                    }
                    
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        ddK_cu[0] = mZeta * thickness /2.0 * DDw_DD1s[dirr] * DN_De(kr, 0);
                        ddK_cu[1] = mZeta * thickness /2.0 * DDw_DD2s[dirr] * DN_De(kr, 1);
                        ddK_cu[2] = mZeta * thickness /2.0 * 0.5 * (DDw_DD1s[dirr] * DN_De(kr, 1) + DDw_DD2s[dirr] * DN_De(kr, 0));
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
                        ddK_cu[0] += mZeta * thickness /2.0 * DDw_DD1r[dirs] * DN_De(ks, 0);
                        ddK_cu[1] += mZeta * thickness /2.0 * DDw_DD2r[dirs] * DN_De(ks, 1);
                        ddK_cu[2] += mZeta * thickness /2.0 * 0.5 * (DDw_DD1r[dirs] * DN_De(ks, 1) + DDw_DD2r[dirs] * DN_De(ks, 0));
                        }
                    else if (dirs == 3 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += DN_De(ks, 0) * DN_De(kr, 0) + N(ks) * DDN_DDe(kr, 0);
                        DDDw_DDD2rs[dirr] += DN_De(ks, 1) * DN_De(kr, 0) + N(ks) * DDN_DDe(kr, 2);
                    }
                    else if (dirs == 4 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += DN_De(ks, 0) * DN_De(kr, 1) + N(ks) * DDN_DDe(kr, 2);
                        DDDw_DDD2rs[dirr] += DN_De(ks, 1) * DN_De(kr, 1) + N(ks) * DDN_DDe(kr, 1);
                    }
                    ddK_cu[0] += mZeta * thickness / 2.0 * inner_prod(DDDw_DDD1rs, rActualMetric.a1);
                    ddK_cu[1] += mZeta * thickness / 2.0 * inner_prod(DDDw_DDD2rs, rActualMetric.a2);
                    ddK_cu[2] += mZeta * thickness / 2.0 * 0.5 * (inner_prod(DDDw_DDD1rs, rActualMetric.a2) + inner_prod(DDDw_DDD2rs, rActualMetric.a1));

                    // calculated with simplified Q (ML)
                    rSecondVariations.B11(r, s) += mInitialMetric.Q(0, 0) * ddK_cu[0];
                    rSecondVariations.B22(r, s) += mInitialMetric.Q(1, 0) * ddK_cu[0] + mInitialMetric.Q(1, 1) * ddK_cu[1]
                        + mInitialMetric.Q(1, 2) * ddK_cu[2];
                    rSecondVariations.B12(r, s) += mInitialMetric.Q(2, 0) * ddK_cu[0] + mInitialMetric.Q(2, 2) * ddK_cu[2];
                    if (r != s){
                        rSecondVariations.B11(s, r) += rSecondVariations.B11(r, s);
                        rSecondVariations.B22(s, r) += rSecondVariations.B22(r, s);
                        rSecondVariations.B12(s, r) += rSecondVariations.B12(r, s);
                    }
                    // all other entries are (remain) zero
                }
            }
        }
    }

    void IgaShell5pElement::Calculate(
        const Variable<double>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions = Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // shear difference vector
        array_1d<double, 3> w = ZeroVector(3);
        // derivatives of the shear difference vector
        array_1d<double, 3> Dw_D1 = ZeroVector(3);
        array_1d<double, 3> Dw_D2 = ZeroVector(3);
        // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * a1 + w_alpha(2) * a2)
        array_1d<double, 2> w_alpha = ZeroVector(2);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);

        std::vector<array_1d<double, 5>> stress_pk2_cart(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_pk2_cov(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_cau_cov(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_cau_cart(mGaussQuadratureThickness.num_GP_thickness);

        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);
        
        double thickness = GetProperties().GetValue(THICKNESS);

        // the Gauss-Points start from bottom to top
        for (unsigned int Gauss_index = 0; Gauss_index < mGaussQuadratureThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussQuadratureThickness.zeta(Gauss_index);

            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, w, 
                Dw_D1, Dw_D2, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

            array_1d<double, 3> G1 = ZeroVector(3);
            array_1d<double, 3> G2 = ZeroVector(3);
            array_1d<double, 3> g1 = ZeroVector(3);
            array_1d<double, 3> g2 = ZeroVector(3);
            array_1d<double, 3> g3 = ZeroVector(3);
            Matrix F = ZeroMatrix(3, 3);
            double detF = 0.0;
            CalculateInitialBaseVectorsGLinearized(mZeta, G1, G2);
            CalculateActualBaseVectorsgLinearized(actual_metric, w, Dw_D1, Dw_D2, g1, g2, g3);
            CalculateDeformationGradient(G1, G2, g1, g2, g3, F, detF);

            // stresses at GP
            stress_pk2_cart[Gauss_index] = constitutive_variables.S;
            stress_pk2_cov[Gauss_index] = prod(mInitialMetric.TransCartToCov, stress_pk2_cart[Gauss_index]);
            stress_cau_cov[Gauss_index] = stress_pk2_cov[Gauss_index] / detF;
            stress_cau_cart[Gauss_index] = prod(actual_metric.TransCovToCart, stress_cau_cov[Gauss_index]);
        }
    
        // Cauchy stress at midspan
        array_1d<double, 5> stress_cau_cart_mid;
        for (unsigned int i = 0; i < 5; i++)
            stress_cau_cart_mid[i] = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][i] + stress_cau_cart[0][i]) / 2.0;

        // internal forces n11, n22, n12, n23, n13
        array_1d<double, 5> n = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1] + stress_cau_cart[0]) / 2.0 * 
            thickness;

        // internal moments m11, m22, m12
        array_1d<double, 3> m;
        m[0] = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][0] - stress_cau_cart_mid[0]) * thickness * thickness / 
            (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        m[1] = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][1] - stress_cau_cart_mid[1]) * thickness * thickness / 
            (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        m[2] = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][2] - stress_cau_cart_mid[2]) * thickness * thickness / 
            (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        
        // stresses at the top (positive theta_3 direction)
        array_1d<double, 5> stress_cau_cart_top = stress_cau_cart_mid + (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1] 
            - stress_cau_cart_mid) / mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1);
        // stresses at the bottom (negative theta_3 direction)
        array_1d<double, 5> stress_cau_cart_bottom = stress_cau_cart_mid + (stress_cau_cart[0] - stress_cau_cart_mid) / 
            mGaussQuadratureThickness.zeta(0);

        if (rVariable == STRESS_CAUCHY_TOP_11)
            rValues = stress_cau_cart_top[0];
        else if (rVariable == STRESS_CAUCHY_TOP_22)
            rValues = stress_cau_cart_top[1];
        else if (rVariable == STRESS_CAUCHY_TOP_12)
            rValues = stress_cau_cart_top[2];
        else if (rVariable == STRESS_CAUCHY_TOP_23)
            rValues = stress_cau_cart_top[3];
        else if (rVariable == STRESS_CAUCHY_TOP_13)
            rValues = stress_cau_cart_top[4];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_11)
            rValues = stress_cau_cart_bottom[0];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_22)
            rValues = stress_cau_cart_bottom[1];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_12)
            rValues = stress_cau_cart_bottom[2];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_23)
            rValues = stress_cau_cart_bottom[3];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_13)
            rValues = stress_cau_cart_bottom[4];
        else if (rVariable == INTERNAL_FORCE_11)
            rValues = n[0];
        else if (rVariable == INTERNAL_FORCE_22)
            rValues = n[1];
        else if (rVariable == INTERNAL_FORCE_12)
            rValues = n[2];
        else if (rVariable == INTERNAL_MOMENT_11)
            rValues = m[0];
        else if (rVariable == INTERNAL_MOMENT_12)
            rValues = m[2];
        else if (rVariable == INTERNAL_MOMENT_22)
            rValues = m[1];
        else if (rVariable == SHEAR_FORCE_1)
            rValues = n[4];     // q1=n13
        else if (rVariable == SHEAR_FORCE_2)
            rValues = n[3];     // q2=n23
        else{
            KRATOS_WATCH("No results for desired variable available in Calculate of IgaShell5pElement.")
            rValues = 0.0;
        }
    }

    void IgaShell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

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

    void IgaShell5pElement::GetDofList(
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
