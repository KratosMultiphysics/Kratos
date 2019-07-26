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
#include "custom_elements/iga_shell_5p_element_stuttgart.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell5pElementStuttgart::Initialize()
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

    void IgaShell5pElementStuttgart::CalculateAll(
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
        array_1d<double, 3> G1, G2, G1xG2;
        double thickness = GetProperties().GetValue(THICKNESS);

        for (unsigned int Gauss_index = 0; Gauss_index < mGaussQuadratureThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussQuadratureThickness.zeta(Gauss_index);
            
            // Differential Volume
            G1 = ZeroVector(3);
            G2 = ZeroVector(3);
            G1xG2 = ZeroVector(3);
            CalculateInitialBaseVectorsGLinearized(G1, G2);
            MathUtils<double>::CrossProduct(G1xG2, G1, G2);
            dV = inner_prod(G1xG2, mInitialMetric.a3_KL);

            array_1d<double, 5> strain_vector = ZeroVector(5);
            Matrix B = ZeroMatrix(5, mat_size);
            boperator_nln_linearisiert(B, strain_vector, actual_metric, CalculateStiffnessMatrixFlag);

            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, strain_vector, constitutive_variables, Values, 
                ConstitutiveLaw::StressMeasure_PK2);

            double integration_weight = mGaussQuadratureThickness.integration_weight_thickness(Gauss_index) * 
                GetValue(INTEGRATION_WEIGHT) * dV * thickness / 2.0;

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                Matrix IKg = ZeroMatrix(mat_size, mat_size);
                kgeom_linearisiert(IKg, constitutive_variables.S, actual_metric);
                
                // adding linear contributions to the stiffness matrix
                noalias(rLeftHandSideMatrix) += integration_weight * prod(trans(B), Matrix(prod(constitutive_variables.D, B)));

                // adding  non-linear contribution to stiffness matrix
                noalias(rLeftHandSideMatrix) += integration_weight * IKg;
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

    void IgaShell5pElementStuttgart::CalculateAndAddKm(
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

    void IgaShell5pElementStuttgart::CalculateAndAddNonlinearKm(
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

    void IgaShell5pElementStuttgart::CalculateMetric(
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
        MathUtils<double>::CrossProduct(rMetric.a3_tilde, rMetric.a1, rMetric.a2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.a3_tilde);
        //normalized basis vector a3_KL
        rMetric.a3_KL = rMetric.a3_tilde / rMetric.dA;
        
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
        double mG_01 = inner_prod(e1, rMetric.a2_con);
        double mG_02 = inner_prod(e1, rMetric.a3_KL);
        double mG_10 = inner_prod(e2, rMetric.a1_con);
        double mG_11 = inner_prod(e2, rMetric.a2_con);
        double mG_12 = inner_prod(e2, rMetric.a3_KL);
        double mG_20 = inner_prod(rMetric.a3_KL, rMetric.a1_con);
        double mG_21 = inner_prod(rMetric.a3_KL, rMetric.a2_con);
        double mG_22 = inner_prod(rMetric.a3_KL, rMetric.a3_KL);
        
        rMetric.Q(0, 0) = pow(mG_00, 2);
        rMetric.Q(0, 1) = pow(mG_01, 2);
        rMetric.Q(0, 2) = 2.00 * mG_00 * mG_01;
        rMetric.Q(0, 3) = 2.00 * mG_01 * mG_02;
        rMetric.Q(0, 4) = 2.00 * mG_00 * mG_02;

        rMetric.Q(1, 0) = pow(mG_10, 2);
        rMetric.Q(1, 1) = pow(mG_11, 2);
        rMetric.Q(1, 2) = 2.00 * mG_10 * mG_11;
        rMetric.Q(1, 3) = 2.00 * mG_11 * mG_12;
        rMetric.Q(1, 4) = 2.00 * mG_10 * mG_12;

        rMetric.Q(2, 0) = 2.00 * mG_00 * mG_10;
        rMetric.Q(2, 1) = 2.00 * mG_01 * mG_11;
        rMetric.Q(2, 2) = 2.00 * (mG_00 * mG_11 + mG_01 * mG_10);
        rMetric.Q(2, 3) = 2.00 * (mG_01 * mG_12 + mG_02 * mG_11);
        rMetric.Q(2, 4) = 2.00 * (mG_00 * mG_12 + mG_02 * mG_10);

        rMetric.Q(3, 0) = 2.00 * mG_10 * mG_20;
        rMetric.Q(3, 1) = 2.00 * mG_11 * mG_21;
        rMetric.Q(3, 2) = 2.00 * (mG_10 * mG_21 + mG_11 * mG_20);
        rMetric.Q(3, 3) = 2.00 * (mG_11 * mG_22 + mG_12 * mG_21);
        rMetric.Q(3, 4) = 2.00 * (mG_10 * mG_22 + mG_12 * mG_20);

        rMetric.Q(4, 0) = 2.00 * (mG_00 * mG_20);
        rMetric.Q(4, 1) = 2.00 * (mG_01 * mG_21);
        rMetric.Q(4, 2) = 2.00 * (mG_00 * mG_21 + mG_01 * mG_20);
        rMetric.Q(4, 3) = 2.00 * (mG_01 * mG_22 + mG_02 * mG_21);
        rMetric.Q(4, 4) = 2.00 * (mG_00 * mG_22 + mG_02 * mG_20);

        // transformation matrix TransCartToCov from local Cartesian to covariant basis
        rMetric.TransCartToCov = trans(rMetric.Q);
        // division by 2.0 because not used for strains but for stresses (strains have e.g. entries with 2e12)
        rMetric.TransCartToCov(2, 1) = rMetric.TransCartToCov(2, 1) / 2.0;
        rMetric.TransCartToCov(2, 2) = rMetric.TransCartToCov(2, 2) / 2.0;
        rMetric.TransCartToCov(3, 3) = rMetric.TransCartToCov(3, 3) / 2.0;
        rMetric.TransCartToCov(4, 3) = rMetric.TransCartToCov(4, 3) / 2.0;
        rMetric.TransCartToCov(4, 4) = rMetric.TransCartToCov(4, 4) / 2.0;

        mG_00 = inner_prod(e1, rMetric.a1);
        mG_01 = inner_prod(e1, rMetric.a2);
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

        if(Id() == 4)
            KRATOS_WATCH(rMetric.Q)
    }

    void IgaShell5pElementStuttgart::CalculateShearDifferenceVector(
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
    
    void IgaShell5pElementStuttgart::CalculateInitialBaseVectorsGLinearized(
        array_1d<double, 3>&      rG1,
        array_1d<double, 3>&      rG2)
    {
        array_1d<double, 3> DA3_D1 = ZeroVector(3);
        array_1d<double, 3> DA3_D2 = ZeroVector(3);
        array_1d<double, 3> DA1_D1xA2 = ZeroVector(3);
        array_1d<double, 3> A1xDA2_D1 = ZeroVector(3);
        array_1d<double, 3> DA1_D2xA2 = ZeroVector(3);
        array_1d<double, 3> A1xDA2_D2 = ZeroVector(3);

        MathUtils<double>::CrossProduct(DA1_D1xA2, mInitialMetric.Da1_D1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D1, mInitialMetric.a1, mInitialMetric.Da1_D2); // DA1_D2 = DA2_D1
        MathUtils<double>::CrossProduct(DA1_D2xA2, mInitialMetric.Da1_D2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D2, mInitialMetric.a1, mInitialMetric.Da2_D2);
        DA3_D1 = ((DA1_D1xA2 + A1xDA2_D1) * mInitialMetric.dA - mInitialMetric.a3_tilde * norm_2(DA1_D1xA2 + A1xDA2_D1)) 
            / (mInitialMetric.dA * mInitialMetric.dA);
        DA3_D2 = ((DA1_D2xA2 + A1xDA2_D2) * mInitialMetric.dA - mInitialMetric.a3_tilde * norm_2(DA1_D2xA2 + A1xDA2_D2))
            / (mInitialMetric.dA * mInitialMetric.dA);

        // covariant base vectors of the shell body in the reference configuration
        rG1 = mInitialMetric.a1 + mZeta * DA3_D1;
        rG2 = mInitialMetric.a2 + mZeta * DA3_D2;
        // G3 = A3
    }

    void IgaShell5pElementStuttgart::CalculateActualBaseVectorsgLinearized(
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
        DdA_D1 = (Da1xa2_D1[0] * rActualMetric.a3_tilde[0] 
            + Da1xa2_D1[1] * rActualMetric.a3_tilde[1] + Da1xa2_D1[2] * rActualMetric.a3_tilde[2]) / rActualMetric.dA;

        /* Ableitung des Nenners von a3_KL nach 2 */
        DdA_D2 = (Da1xa2_D2[0] * rActualMetric.a3_tilde[0] 
            + Da1xa2_D2[1] * rActualMetric.a3_tilde[1] + Da1xa2_D2[2] * rActualMetric.a3_tilde[2]) / rActualMetric.dA;

        /* Ableitung von a3_KL mit Quotientenregel */
        for (unsigned int i = 0; i < 3; i++)
        {
            Da3_KL_D1[i] = (Da1xa2_D1[i] * rActualMetric.dA - rActualMetric.a3_tilde[i] * DdA_D1) / 
                (rActualMetric.dA * rActualMetric.dA);
            Da3_KL_D2[i] = (Da1xa2_D2[i] * rActualMetric.dA - rActualMetric.a3_tilde[i] * DdA_D2) / 
                (rActualMetric.dA * rActualMetric.dA);
        }

        /* Kovariante Basisvektoren */
        rg1 = rActualMetric.a1 + mZeta*(Da3_KL_D1 + rDw_D1);
        rg2 = rActualMetric.a2 + mZeta*(Da3_KL_D2 + rDw_D2);
        rg3 = rActualMetric.a3_KL + rw;     // g3 = a3

        KRATOS_CATCH("")
    }

    void IgaShell5pElementStuttgart::CalculateDeformationGradient(
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

    void IgaShell5pElementStuttgart::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const array_1d<double, 5>& rStrainVector,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        KRATOS_TRY

        rThisConstitutiveVariables.E = rStrainVector;
  
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

    void IgaShell5pElementStuttgart::CalculateStrain(
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

    void IgaShell5pElementStuttgart::CalculateStrainRM(
        Vector& rStrainVectorRM,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rg1,
        const Vector& rg2)
    {
        double thickness = GetProperties().GetValue(THICKNESS);
                
        rStrainVectorRM[0] = thickness/2.0*mZeta *inner_prod(rDw_D1, rg1);     // MLt
        rStrainVectorRM[1] = thickness/2.0*mZeta *inner_prod(rDw_D2, rg2);     // MLt
        rStrainVectorRM[2] = thickness/2.0*mZeta * 0.5 *(inner_prod(rDw_D1, rg2) + inner_prod(rDw_D2, rg1));     // MLt
        rStrainVectorRM[3] = 0.5 *inner_prod(rw, rg2);
        rStrainVectorRM[4] = 0.5 *inner_prod(rw, rg1);
    }

    void IgaShell5pElementStuttgart::TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector& rCurvilinearStrain,
        Vector& rCartesianStrain)
    {
        KRATOS_TRY

        if (rCurvilinearStrain.size() != 5 || rCartesianStrain.size() != 6) 
            KRATOS_ERROR << "Wrong strain size in transformation." << std::endl;
        if (mInitialMetric.Q.size1() != 5 || mInitialMetric.Q.size2() != 5)
            KRATOS_ERROR << "Wrong size of transformation matrix Q." << std::endl;

        // transformation with simplified matrix
        rCartesianStrain[0] = mInitialMetric.Q(0, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(0, 1) * rCurvilinearStrain[1] +
            mInitialMetric.Q(0, 2) * rCurvilinearStrain[2] + mInitialMetric.Q(0, 3) * rCurvilinearStrain[3] + 
            mInitialMetric.Q(0, 4) * rCurvilinearStrain[4];
        rCartesianStrain[1] = mInitialMetric.Q(1, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(1, 1) * rCurvilinearStrain[1] 
            + mInitialMetric.Q(1, 2) * rCurvilinearStrain[2] + mInitialMetric.Q(1, 3) * rCurvilinearStrain[3] +
            mInitialMetric.Q(1, 4) * rCurvilinearStrain[4];
        rCartesianStrain[2] = 0.0; // RM
        rCartesianStrain[3] = mInitialMetric.Q(2, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(2, 1) * rCurvilinearStrain[1] +
            mInitialMetric.Q(2, 2) * rCurvilinearStrain[2] + mInitialMetric.Q(2, 3) * rCurvilinearStrain[3] +
            mInitialMetric.Q(2, 4) * rCurvilinearStrain[4];
        rCartesianStrain[4] = mInitialMetric.Q(3, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(3, 1) * rCurvilinearStrain[1] +
            mInitialMetric.Q(3, 2) * rCurvilinearStrain[2] + mInitialMetric.Q(3, 3) * rCurvilinearStrain[3] + 
            mInitialMetric.Q(3, 4) * rCurvilinearStrain[4];
        rCartesianStrain[5] = mInitialMetric.Q(4, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(4, 1) * rCurvilinearStrain[1] +
            mInitialMetric.Q(4, 2) * rCurvilinearStrain[2] + mInitialMetric.Q(4, 3) * rCurvilinearStrain[3] + 
            mInitialMetric.Q(4, 4) * rCurvilinearStrain[4];

        KRATOS_CATCH("")
    }

    void IgaShell5pElementStuttgart::CalculateB(
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
                double g3dg3lg3 = (rMetric.a3_tilde[0] * dg3(j, 0) + rMetric.a3_tilde[1] * dg3(j, 1) + rMetric.a3_tilde[2] * dg3(j, 2))*inddA3;

                dn(j, 0) = dg3(j, 0)*invdA - rMetric.a3_tilde[0] * g3dg3lg3;
                dn(j, 1) = dg3(j, 1)*invdA - rMetric.a3_tilde[1] * g3dg3lg3;
                dn(j, 2) = dg3(j, 2)*invdA - rMetric.a3_tilde[2] * g3dg3lg3;
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

    void IgaShell5pElementStuttgart::CalculateSecondVariations(
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

                S_g3dg3[r] = rMetric.a3_tilde[0] * S_dg3(0, r) + rMetric.a3_tilde[1] * S_dg3(1, r) + rMetric.a3_tilde[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.a3_tilde[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.a3_tilde[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.a3_tilde[2] * S_g3dg3lg3_3[r];
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

                    double c = -(ddg3[0] * rMetric.a3_tilde[0] + ddg3[1] * rMetric.a3_tilde[1] + ddg3[2] * rMetric.a3_tilde[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.a3_tilde[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.a3_tilde[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.a3_tilde[2];

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

    void IgaShell5pElementStuttgart::boperator_nln_linearisiert(        
        Matrix& rB,
        array_1d<double, 5>& Egl,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag)
    {
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;

        const Vector& funct = N;
        const Matrix& deriv = DN_De;
        const Matrix& s_deriv = DDN_DDe;
        const unsigned int num_node = number_of_nodes;
        const unsigned int dof_per_node = 5;

        static array_1d<double, 3> v_d1;               // Verschiebungsableitung nach 1
        static array_1d<double, 3> v_d2;               // Verschiebungsableitung nach 2

        /* REFERENZKONFIGURATION */
        int i, j, k;
        double h, deriv_sqrt_normA1crossA2, deriv_radicant_normA1crossA2_dT1,
            deriv_radicant_normA1crossA2_dT2, dnormA1crossA2_dT1, dnormA1crossA2_dT2;
        array_1d<double, 3> A1;
        array_1d<double, 3> A2;
        array_1d<double, 3> A3;
        array_1d<double, 3> A1xA2;
        array_1d<double, 3> dA1_dT1;
        array_1d<double, 3> dA2_dT2;
        array_1d<double, 3> dA1_dT2;
        array_1d<double, 3> dA2_dT1;
        array_1d<double, 3> d_A1crossA2_dT1;
        array_1d<double, 3> d_A1crossA2_dT2;
        array_1d<double, 3> dA3_dT1;
        array_1d<double, 3> dA3_dT2;

        /* Initialisierung aller benoetigter Vektoren */
        A1 = ZeroVector(3);
        A2 = ZeroVector(3);
        A3 = ZeroVector(3);
        A1xA2 = ZeroVector(3);
        dA1_dT1 = ZeroVector(3);
        dA2_dT2 = ZeroVector(3);
        dA1_dT2 = ZeroVector(3);
        dA2_dT1 = ZeroVector(3);
        d_A1crossA2_dT1 = ZeroVector(3);
        d_A1crossA2_dT2 = ZeroVector(3);
        deriv_radicant_normA1crossA2_dT2 = 0.0;
        deriv_radicant_normA1crossA2_dT2 = 0.0;
        dnormA1crossA2_dT1 = 0.0;
        dnormA1crossA2_dT2 = 0.0;
        dA3_dT1 = ZeroVector(3);
        dA3_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der Referenzkonfiguration */
        A1 = mInitialMetric.a1;
        A2 = mInitialMetric.a2;

        /* Dritter kovarianter Basisvektor senkrecht zur Mittelflaeche der Referenzkonfiguration */

        A1xA2 = mInitialMetric.a3_tilde;
        h = sqrt(A1xA2[0] * A1xA2[0] + A1xA2[1] * A1xA2[1] + A1xA2[2] * A1xA2[2]);

        for (k = 0; k < 3; k++)
        {
            A3[k] = A1xA2[k] / h * thickness / 2;
        }

        /* Ableitungen von A1 und A2 nach der Richtung alpha */
        dA1_dT1 = mInitialMetric.Da1_D1;
        dA2_dT2 = mInitialMetric.Da2_D2;
        dA1_dT2 = mInitialMetric.Da1_D2;
        dA2_dT1 = dA1_dT2;

        /* Ableitung von A3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von A3 nach 1: (A1xA2)'= A1'xA2 + A1xA2' */
        d_A1crossA2_dT1 = CrossProduct(dA1_dT1, A2) + CrossProduct(A1, dA2_dT1);

        /* Ableitung des Zaehlers von A3 nach 2: (A1xA2)'= A1'xA2 + A1xA2' */
        d_A1crossA2_dT2 = CrossProduct(dA1_dT2, A2) + CrossProduct(A1, dA2_dT2);

        /* Ableitung des Nenners von A3 */
        deriv_sqrt_normA1crossA2 = 1.0 / (2.0 * sqrt(A1xA2[0] * A1xA2[0] + A1xA2[1] * A1xA2[1] + A1xA2[2] * A1xA2[2]));

        /* Ableitung des Nenners von A3 nach 1 */
        deriv_radicant_normA1crossA2_dT1 = 2.0 * d_A1crossA2_dT1[0] * A1xA2[0] + 2.0 * d_A1crossA2_dT1[1] * A1xA2[1] + 2.0 * d_A1crossA2_dT1[2] * A1xA2[2];
        dnormA1crossA2_dT1 = deriv_sqrt_normA1crossA2 * deriv_radicant_normA1crossA2_dT1;

        /* Ableitung des Nenners von A3 nach 2 */
        deriv_radicant_normA1crossA2_dT2 = 2.0 * d_A1crossA2_dT2[0] * A1xA2[0] + 2.0 * d_A1crossA2_dT2[1] * A1xA2[1] + 2.0 * d_A1crossA2_dT2[2] * A1xA2[2];
        dnormA1crossA2_dT2 = deriv_sqrt_normA1crossA2 * deriv_radicant_normA1crossA2_dT2;

        /* Ableitung von A3 mit Quotientenregel */
        for (k = 0; k < 3; k++)
        {
            dA3_dT1[k] = (d_A1crossA2_dT1[k] * h - A1xA2[k] * dnormA1crossA2_dT1) / (h * h) * thickness / 2.0;
            dA3_dT2[k] = (d_A1crossA2_dT2[k] * h - A1xA2[k] * dnormA1crossA2_dT2) / (h * h) * thickness / 2.0;
        }

        /* AKTUELLE KONFIGURATION */
        int m;
        double hact, deriv_sqrt_norma1crossa2, deriv_radicant_norma1crossa2_dT1,
            deriv_radicant_norma1crossa2_dT2, dnorma1crossa2_dT1, dnorma1crossa2_dT2;
        array_1d<double, 3> a1;
        array_1d<double, 3> a2;
        array_1d<double, 3> a3;
        array_1d<double, 3> a3_KL;
        array_1d<double, 3> a1xa2;
        array_1d<double, 3> da1_dT1;
        array_1d<double, 3> da2_dT2;
        array_1d<double, 3> da1_dT2;
        array_1d<double, 3> da2_dT1;
        array_1d<double, 3> d_a1crossa2_dT1;
        array_1d<double, 3> d_a1crossa2_dT2;
        array_1d<double, 3> da3_dT1;
        array_1d<double, 3> da3_dT2;
        array_1d<double, 3> w;
        array_1d<double, 3> da1norm_dT1;
        array_1d<double, 3> da1norm_dT2;
        array_1d<double, 3> da2norm_dT1;
        array_1d<double, 3> da2norm_dT2;


        /* Initialisierung aller benoetigter Vektoren */
        a1 = ZeroVector(3);
        a2 = ZeroVector(3);
        a3 = ZeroVector(3);
        a3_KL = ZeroVector(3);
        a1xa2 = ZeroVector(3);
        da1_dT1 = ZeroVector(3);
        da2_dT2 = ZeroVector(3);
        da1_dT2 = ZeroVector(3);
        da2_dT1 = ZeroVector(3);
        d_a1crossa2_dT1 = ZeroVector(3);
        d_a1crossa2_dT2 = ZeroVector(3);
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        dnorma1crossa2_dT1 = 0.0;
        dnorma1crossa2_dT2 = 0.0;
        da3_dT1 = ZeroVector(3);
        da3_dT2 = ZeroVector(3);
        w = ZeroVector(3);
        da1norm_dT1 = ZeroVector(3);
        da1norm_dT2 = ZeroVector(3);
        da2norm_dT1 = ZeroVector(3);
        da2norm_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der aktuellen Konfiguration */
        a1 = rActualMetric.a1;
        a2 = rActualMetric.a2;


        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        a1xa2 = rActualMetric.a3_tilde;
        hact = sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]);

        for (m = 0; m < 3; m++)
            a3_KL[m] = a1xa2[m] / hact * thickness / 2.0;


        /*  Normierte Basisvektoren in aktueller Konfiguration */
        double              a1n = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
        double              a2n = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
        array_1d<double, 3>     a1norm = (1.0/a1n)*a1;
        array_1d<double, 3>     a2norm = (1.0/a2n)*a2;
        double w_1,w_2;
        w_1=0.0;
        w_2=0.0;

        /* aktuelle Komponenten vom Schubdifferenzvektor w */
        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
        for (k = 0; k < num_node; k++)
        {
            w_1   += funct[k] * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            w_2   += funct[k] * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }  

        /* linearisierter Schubdifferenzvektor w */
        w = (w_1)*a1 + (w_2)*a2;

        /* a3 nach hier. 5p-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
        a3 = a3_KL + w;

        /* Ableitungen von a1 und a2 nach der Richtung alpha */
        da1_dT1 = rActualMetric.Da1_D1;
        da2_dT2 = rActualMetric.Da2_D2;
        da1_dT2 = rActualMetric.Da1_D2;
        da2_dT1 = da1_dT2;

        /* Ableitung von a3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von a3 nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        d_a1crossa2_dT1 = CrossProduct(da1_dT1, a2) + CrossProduct(a1, da2_dT1);

        /* Ableitung des Zaehlers von a3 nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        d_a1crossa2_dT2 = CrossProduct(da1_dT2, a2) + CrossProduct(a1, da2_dT2);

        /* Ableitung des Nenners von a3 */
        deriv_sqrt_norma1crossa2 = 1.0 / (2.0 * sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]));

        /* Ableitung des Nenners von a3 nach 1 */
        deriv_radicant_norma1crossa2_dT1 = 2.0 * d_a1crossa2_dT1[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT1[1] * a1xa2[1] + 2.0 * d_a1crossa2_dT1[2] * a1xa2[2];
        dnorma1crossa2_dT1 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT1;

        /* Ableitung des Nenners von a3 nach 2 */
        deriv_radicant_norma1crossa2_dT2 = 2.0 * d_a1crossa2_dT2[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT2[1] * a1xa2[1] + 2.0 * d_a1crossa2_dT2[2] * a1xa2[2];
        dnorma1crossa2_dT2 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT2;

        /* Ableitung von A3 mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_dT1[m] = (d_a1crossa2_dT1[m] * hact - a1xa2[m] * dnorma1crossa2_dT1) / (hact * hact) * thickness / 2.0;
            da3_dT2[m] = (d_a1crossa2_dT2[m] * hact - a1xa2[m] * dnorma1crossa2_dT2) / (hact * hact) * thickness / 2.0;
        }

        /* Partielle Ableitungen von Schubdifferenzvektor w nach alpha/beta */
        double dw_1dT1,dw_2dT1,dw_1dT2,dw_2dT2;
        dw_1dT1=0.0;
        dw_2dT1=0.0;
        dw_1dT2=0.0;
        dw_2dT2=0.0;
        
        for (k = 0; k < num_node; k++)
        {
            dw_1dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw_2dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
            dw_1dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw_2dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }

        array_1d<double, 3>  dw_dT1 = dw_1dT1*a1 + w_1*da1_dT1 + dw_2dT1*a2 + w_2*da2_dT1;
        array_1d<double, 3>  dw_dT2 = dw_1dT2*a1 + w_1*da1_dT2 + dw_2dT2*a2 + w_2*da2_dT2;


        /* Nichtlinearer B-Operator */

        array_1d<double, 3> da1xa2;
        array_1d<double, 3> da1;
        array_1d<double, 3> da2;
        array_1d<double, 3> da3;

        da1xa2 = ZeroVector(3);
        da1 = ZeroVector(3);
        da2 = ZeroVector(3);
        da3 = ZeroVector(3);

        double da1xa2_a1xa2;
        double da1xa2_a1xa2_hact;

        int count1 = 0;
        int count2 = 0;


        /* Anteile aus 3P-Formulierung (nicht veraendert) */
        for (j = 0; j < num_node; j++)
        {
            for (i = 0; i <dof_per_node; i++)
            {
            if (i == 0)
            {
                da1[0] = deriv(count2, 0);
                da1[1] = 0;
                da1[2] = 0;
                da2[0] = deriv(count2, 1);
                da2[1] = 0;
                da2[2] = 0;
            }
            else if (i == 1)
            {
                da1[0] = 0;
                da1[1] = deriv(count2, 0);
                da1[2] = 0;
                da2[0] = 0;
                da2[1] = deriv(count2, 1);
                da2[2] = 0;
            }
            else if (i == 2)
            {
                da1[0] = 0;
                da1[1] = 0;
                da1[2] = deriv(count2, 0);
                da2[0] = 0;
                da2[1] = 0;
                da2[2] = deriv(count2, 1);
            }

            da1xa2 = CrossProduct(da1, a2) + CrossProduct(a1, da2);
            da1xa2_a1xa2 = (da1xa2[0] * a1xa2[0] + da1xa2[1] * a1xa2[1]
                + da1xa2[2] * a1xa2[2]);
            da1xa2_a1xa2_hact = da1xa2_a1xa2 / (hact * hact * hact);
            da3[0] = da1xa2[0] / hact - a1xa2[0] * da1xa2_a1xa2_hact;
            da3[1] = da1xa2[1] / hact - a1xa2[1] * da1xa2_a1xa2_hact;
            da3[2] = da1xa2[2] / hact - a1xa2[2] * da1xa2_a1xa2_hact;

            array_1d<double, 5> dE_cur_help = ZeroVector(5);
            if (i<3)
            {
                /* a3_KL anstatt a3! ansonsten wie bei 3P-Schale */
                dE_cur_help[0] = deriv(count2, 0)*a1[i] - mZeta*(s_deriv(count2, 0)*a3_KL[i] + thickness/2.0*(da1_dT1[0]*da3[0] + da1_dT1[1]*da3[1] + da1_dT1[2]*da3[2])) ;		// FML: thickness/2.0 - nicht in 2531??;	ML: warum a1, da1_T1 und da3 anstatt A1, dA1_T1 und A3?? Umformung möglich
                dE_cur_help[2] = 0.5*(deriv(count2, 1)*a1[i] + deriv(count2, 0)*a2[i])																					// ML: Umformung möglich
                                            - mZeta*(s_deriv(count2, 2)*a3_KL[i] + thickness/2.0*(da1_dT2[0]*da3[0] + da1_dT2[1]*da3[1] + da1_dT2[2]*da3[2]));
                dE_cur_help[1] = deriv(count2, 1)*a2[i] - mZeta*(s_deriv(count2, 1)*a3_KL[i] + thickness/2.0*(da2_dT2[0]*da3[0] + da2_dT2[1]*da3[1] + da2_dT2[2]*da3[2]))  ;
            }
            rB(0, j*5 + i) += mInitialMetric.Q(0, 0) * dE_cur_help[0] + mInitialMetric.Q(0, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(0, 2) * dE_cur_help[2] + mInitialMetric.Q(0, 3) * dE_cur_help[3] + mInitialMetric.Q(0, 4) * dE_cur_help[4];
            rB(1, j*5 + i) += mInitialMetric.Q(1, 0) * dE_cur_help[0] + mInitialMetric.Q(1, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(1, 2) * dE_cur_help[2] + mInitialMetric.Q(1, 3) * dE_cur_help[3] + mInitialMetric.Q(1, 4) * dE_cur_help[4];
            rB(2, j*5 + i) += mInitialMetric.Q(2, 0) * dE_cur_help[0] + mInitialMetric.Q(2, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(2, 2) * dE_cur_help[2] + mInitialMetric.Q(2, 3) * dE_cur_help[3] + mInitialMetric.Q(2, 4) * dE_cur_help[4];
            rB(3, j*5 + i) += mInitialMetric.Q(3, 0) * dE_cur_help[0] + mInitialMetric.Q(3, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(3, 2) * dE_cur_help[2] + mInitialMetric.Q(3, 3) * dE_cur_help[3] + mInitialMetric.Q(3, 4) * dE_cur_help[4];
            rB(4, j*5 + i) += mInitialMetric.Q(4, 0) * dE_cur_help[0] + mInitialMetric.Q(4, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(4, 2) * dE_cur_help[2] + mInitialMetric.Q(4, 3) * dE_cur_help[3] + mInitialMetric.Q(4, 4) * dE_cur_help[4];
            }
            count2 ++;
        }


        /* Zusaetzliche Anteile aus linearisiertem Schubvektor w */
        double dw_1_dr;
        double dw_2_dr;
        double dw_1dT1_dr;
        double dw_2dT1_dr;
        double dw_1dT2_dr;
        double dw_2dT2_dr;
        count1 = 0;
        count2 = 0;

        array_1d<double, 3> a1_dr;
        array_1d<double, 3> a2_dr;
        array_1d<double, 3> da1_dT1_dr;
        array_1d<double, 3> da1_dT2_dr;
        array_1d<double, 3> da2_dT1_dr;
        array_1d<double, 3> da2_dT2_dr;

        array_1d<double, 3>  dw_dr;
        array_1d<double, 3>  dw_dT1_dr;
        array_1d<double, 3>  dw_dT2_dr;

        a1_dr = ZeroVector(3);
        a2_dr = ZeroVector(3);
        da1_dT1_dr = ZeroVector(3);
        da1_dT2_dr = ZeroVector(3);
        da2_dT1_dr = ZeroVector(3);
        da2_dT2_dr = ZeroVector(3);
        dw_dr = ZeroVector(3);
        dw_dT1_dr = ZeroVector(3);
        dw_dT2_dr = ZeroVector(3);

        for (j = 0; j < num_node; j++)
        {
            for (i = 0; i < dof_per_node; i++)
            {
            if (i == 0)
            {
                a1_dr[0] = deriv(count2, 0);
                a1_dr[1] = 0.0;
                a1_dr[2] = 0.0;
                a2_dr[0] = deriv(count2, 1);
                a2_dr[1] = 0.0;
                a2_dr[2] = 0.0;

                da1_dT1_dr[0] = s_deriv(count2, 0);
                da1_dT1_dr[1] = 0.0;
                da1_dT1_dr[2] = 0.0;
                da1_dT2_dr[0] = s_deriv(count2, 2);
                da1_dT2_dr[1] = 0.0;
                da1_dT2_dr[2] = 0.0;
                da2_dT1_dr[0] = s_deriv(count2, 2);
                da2_dT1_dr[1] = 0.0;
                da2_dT1_dr[2] = 0.0;
                da2_dT2_dr[0] = s_deriv(count2, 1);
                da2_dT2_dr[1] = 0.0;
                da2_dT2_dr[2] = 0.0;

                dw_2_dr = 0.0;
                dw_1_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 1)
            {
                a1_dr[0] = 0;
                a1_dr[1] = deriv(count2, 0);
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = deriv(count2, 1);
                a2_dr[2] = 0;

                da1_dT1_dr[0] = 0.0;
                da1_dT1_dr[1] = s_deriv(count2, 0);
                da1_dT1_dr[2] = 0.0;
                da1_dT2_dr[0] = 0.0;
                da1_dT2_dr[1] = s_deriv(count2, 2);
                da1_dT2_dr[2] = 0.0;
                da2_dT1_dr[0] = 0.0;
                da2_dT1_dr[1] = s_deriv(count2, 2);
                da2_dT1_dr[2] = 0.0;
                da2_dT2_dr[0] = 0.0;
                da2_dT2_dr[1] = s_deriv(count2, 1);
                da2_dT2_dr[2] = 0.0;

                dw_1_dr = 0.0;
                dw_2_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 2)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = deriv(count2, 0);
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = deriv(count2, 1);

                da1_dT1_dr[0] = 0.0;
                da1_dT1_dr[1] = 0.0;
                da1_dT1_dr[2] = s_deriv(count2, 0);
                da1_dT2_dr[0] = 0.0;
                da1_dT2_dr[1] = 0.0;
                da1_dT2_dr[2] = s_deriv(count2, 2);
                da2_dT1_dr[0] = 0.0;
                da2_dT1_dr[1] = 0.0;
                da2_dT1_dr[2] = s_deriv(count2, 2);
                da2_dT2_dr[0] = 0.0;
                da2_dT2_dr[1] = 0.0;
                da2_dT2_dr[2] = s_deriv(count2, 1);

                dw_1_dr = 0.0;
                dw_2_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 3)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;
                dw_1_dr = funct[count2];
                dw_2_dr = 0.0;
                dw_1dT1_dr = deriv(count2, 0);
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = deriv(count2, 1);
                dw_2dT2_dr = 0.0;
            }
            else if (i == 4)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;
                dw_1_dr = 0.0;
                dw_2_dr = funct[count2];
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = deriv(count2, 0);
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = deriv(count2, 1);
            }

            /* Partielle Ableitung/Variation nach FHG (_dr entspricht ,r) */
            dw_dr = dw_1_dr*a1 + w_1*a1_dr + dw_2_dr*a2 + w_2*a2_dr;
            dw_dT1_dr = dw_1dT1_dr*a1 + dw_1dT1*a1_dr + dw_1_dr*da1_dT1 + w_1*da1_dT1_dr +
                        dw_2dT1_dr*a2 + dw_2dT1*a2_dr + dw_2_dr*da2_dT1 + w_2*da2_dT1_dr;
            dw_dT2_dr = dw_1dT2_dr*a1 + dw_1dT2*a1_dr + dw_1_dr*da1_dT2 + w_1*da1_dT2_dr +
                        dw_2dT2_dr*a2 + dw_2dT2*a2_dr + dw_2_dr*da2_dT2 + w_2*da2_dT2_dr;

            array_1d<double, 5> dE_cur_help = ZeroVector(5);

            dE_cur_help[0] = thickness/2.0*mZeta *(inner_prod(a1_dr, dw_dT1) + inner_prod(a1, dw_dT1_dr));     // MLt
            dE_cur_help[2] = thickness/2.0*0.5 * mZeta*(inner_prod(a1_dr, dw_dT2) + inner_prod(a1, dw_dT2_dr) + inner_prod(a2_dr, dw_dT1) + inner_prod(a2, dw_dT1_dr));     // MLt
            dE_cur_help[1] = thickness/2.0*mZeta*(inner_prod(a2_dr, dw_dT2) + inner_prod(a2, dw_dT2_dr));     // MLt
            dE_cur_help[4] = 0.5 *(inner_prod(dw_dr, a1) + inner_prod(w, a1_dr));
            dE_cur_help[3] = 0.5 *(inner_prod(dw_dr, a2) + inner_prod(w, a2_dr));

            rB(0, count1) += mInitialMetric.Q(0, 0) * dE_cur_help[0] + mInitialMetric.Q(0, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(0, 2) * dE_cur_help[2] + mInitialMetric.Q(0, 3) * dE_cur_help[3] + mInitialMetric.Q(0, 4) * dE_cur_help[4];
            rB(1, count1) += mInitialMetric.Q(1, 0) * dE_cur_help[0] + mInitialMetric.Q(1, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(1, 2) * dE_cur_help[2] + mInitialMetric.Q(1, 3) * dE_cur_help[3] + mInitialMetric.Q(1, 4) * dE_cur_help[4];
            rB(2, count1) += mInitialMetric.Q(2, 0) * dE_cur_help[0] + mInitialMetric.Q(2, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(2, 2) * dE_cur_help[2] + mInitialMetric.Q(2, 3) * dE_cur_help[3] + mInitialMetric.Q(2, 4) * dE_cur_help[4];
            rB(3, count1) += mInitialMetric.Q(3, 0) * dE_cur_help[0] + mInitialMetric.Q(3, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(3, 2) * dE_cur_help[2] + mInitialMetric.Q(3, 3) * dE_cur_help[3] + mInitialMetric.Q(3, 4) * dE_cur_help[4];
            rB(4, count1) += mInitialMetric.Q(4, 0) * dE_cur_help[0] + mInitialMetric.Q(4, 1) * dE_cur_help[1] + 
                mInitialMetric.Q(4, 2) * dE_cur_help[2] + mInitialMetric.Q(4, 3) * dE_cur_help[3] + mInitialMetric.Q(4, 4) * dE_cur_help[4];

        // #if 0
        //     /* zusaetzlicher Term in den Schubverzerrungen, linear in theta^3 */
        //     bop[3][count1] = bop[3][count1] + mZeta*(dw_dT1_dr*w + dw_dT1*dw_dr);
        //     bop[4][count1] = bop[4][count1] + mZeta*(dw_dT2_dr*w + dw_dT2*dw_dr);
        // #endif

            count1 ++;
            }
            count2 ++;
        }

        /* Verschiebungsableitungen */
        v_d1 = ZeroVector(3);
        v_d2 = ZeroVector(3);

        for (k = 0; k < num_node; k++)
        {
            v_d1[0] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_X, pos).GetSolutionStepValue();
            v_d1[1] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_Y, pos + 1).GetSolutionStepValue();
            v_d1[2] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_Z, pos + 2).GetSolutionStepValue();
        }

        for (k = 0; k < num_node; k++)
        {
            v_d2[0] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_X, pos).GetSolutionStepValue();
            v_d2[1] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_Y, pos + 1).GetSolutionStepValue();
            v_d2[2] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_Z, pos + 2).GetSolutionStepValue();
        }

        /* Berechnung der Green-Lagrange Verzerrungen (Voigt-Notation), Anteile von KL(3p)-Schale */
        /* E11 */
        Egl[0] = A1[0]*v_d1[0]+A1[1]*v_d1[1]+A1[2]*v_d1[2] + 0.5*(v_d1[0]*v_d1[0]+v_d1[1]*v_d1[1]+v_d1[2]*v_d1[2])
                + mZeta*(da3_dT1[0]*v_d1[0]+da3_dT1[1]*v_d1[1]+da3_dT1[2]*v_d1[2]+A1[0]*(da3_dT1[0]-dA3_dT1[0])+A1[1]*(da3_dT1[1]-dA3_dT1[1])+A1[2]*(da3_dT1[2]-dA3_dT1[2]));

        /* 2*E12 */
        Egl[2] = 0.5 * (A1[0]*v_d2[0]+A1[1]*v_d2[1]+A1[2]*v_d2[2] + A2[0]*v_d1[0]+A2[1]*v_d1[1]+A2[2]*v_d1[2] + v_d1[0]*v_d2[0]+v_d1[1]*v_d2[1]+v_d1[2]*v_d2[2]
                +mZeta*(A1[0]*(da3_dT2[0]-dA3_dT2[0])+A1[1]*(da3_dT2[1]-dA3_dT2[1])+A1[2]*(da3_dT2[2]-dA3_dT2[2])
                        +A2[0]*(da3_dT1[0]-dA3_dT1[0])+A2[1]*(da3_dT1[1]-dA3_dT1[1])+A2[2]*(da3_dT1[2]-dA3_dT1[2])
                        +v_d1[0]*da3_dT2[0]+v_d1[1]*da3_dT2[1]+v_d1[2]*da3_dT2[2]+v_d2[0]*da3_dT1[0]+v_d2[1]*da3_dT1[1]+v_d2[2]*da3_dT1[2]));

        /* E22 */
        Egl[1] = A2[0]*v_d2[0]+A2[1]*v_d2[1]+A2[2]*v_d2[2] + 0.5*(v_d2[0]*v_d2[0]+v_d2[1]*v_d2[1]+v_d2[2]*v_d2[2])
                + mZeta*(da3_dT2[0]*v_d2[0]+da3_dT2[1]*v_d2[1]+da3_dT2[2]*v_d2[2]+A2[0]*(da3_dT2[0]-dA3_dT2[0])+A2[1]*(da3_dT2[1]-dA3_dT2[1])+A2[2]*(da3_dT2[2]-dA3_dT2[2]));





        /* Zusatzanteile in Green-Lagrange Verzerrungen durch hier. 5-Parameter-Kinematik */
        /* E11 */
        Egl[0] = Egl[0] + thickness/2.0*mZeta*0.5*(inner_prod(a1, dw_dT1) + inner_prod(a1, dw_dT1));    // MLt
        /* 2*E12 */
        Egl[2] = Egl[1] + thickness/2.0*0.5 * mZeta*1.0*(inner_prod(a1, dw_dT2) + inner_prod(a2, dw_dT1));    // MLt
        /* E22 */
        Egl[1] = Egl[2] + thickness/2.0*mZeta*0.5*(inner_prod(a2, dw_dT2) + inner_prod(a2, dw_dT2));    // MLt
        /* 2*E13 */
        Egl[4] =  0.5 * inner_prod(w, a1);
        /* 2*E23 */
        Egl[3] =  0.5 * inner_prod(w, a2);

        // #if 0
        //     /* zusaetzlicher Term in den Schubverzerrungen, linear in theta^3 */
        // /* 2*E13 */
        // Egl[3] = Egl[3] + mZeta*dw_dT1*w;
        // /* 2*E23 */
        // Egl[4] = Egl[4] + mZeta*dw_dT2*w;
        // #endif
    }

    void IgaShell5pElementStuttgart::kgeom_linearisiert(
        Matrix&              IKg,                     ///< Integrand des geometrischen Steifigkeitsmatrix (o)
        const array_1d<double, 5>&               S,                       ///< Zweite Piola-Kirchhoff-Spannungen (i)
        const MetricVariables& rActualMetric
        )
    {
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int mat_size = number_of_nodes * 5;

        const Vector& funct = N;
        const Matrix& deriv = DN_De;
        const Matrix& s_deriv = DDN_DDe;
        const unsigned int num_node = number_of_nodes;
        const unsigned int dof_per_node = 5;

        /* AKTUELLE KONFIGURATION */
        int m, k, l;
        double hact, deriv_sqrt_norma1crossa2, deriv_radicant_norma1crossa2_dT1,
            deriv_radicant_norma1crossa2_dT2, dnorma1crossa2_dT1, dnorma1crossa2_dT2;
        array_1d<double, 3> a1 = ZeroVector(3);
        array_1d<double, 3> a2 = ZeroVector(3);
        array_1d<double, 3> a3 = ZeroVector(3);
        array_1d<double, 3> a3_KL = ZeroVector(3);
        array_1d<double, 3> a1xa2 = ZeroVector(3);
        array_1d<double, 3> da1_dT1 = ZeroVector(3);
        array_1d<double, 3> da1_dT2 = ZeroVector(3);
        array_1d<double, 3> da2_dT1 = ZeroVector(3);
        array_1d<double, 3> da2_dT2 = ZeroVector(3);
        array_1d<double, 3> d_a1crossa2_dT1 = ZeroVector(3);
        array_1d<double, 3> d_a1crossa2_dT2 = ZeroVector(3);
        array_1d<double, 3> da3_KL_dT1 = ZeroVector(3);
        array_1d<double, 3> da3_KL_dT2 = ZeroVector(3);

        /* Initialisierung aller benoetigten Vektoren */
        a1 = ZeroVector(3);
        a2 = ZeroVector(3);
        a3 = ZeroVector(3);
        a3_KL = ZeroVector(3);
        a1xa2 = ZeroVector(3);
        da1_dT1 = ZeroVector(3);
        da1_dT2 = ZeroVector(3);
        da2_dT1 = ZeroVector(3);
        da2_dT2 = ZeroVector(3);
        d_a1crossa2_dT1 = ZeroVector(3);
        d_a1crossa2_dT2 = ZeroVector(3);
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        dnorma1crossa2_dT1 = 0.0;
        dnorma1crossa2_dT2 = 0.0;
        da3_KL_dT1 = ZeroVector(3);
        da3_KL_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der aktuellen Konfiguration */
        a1 = rActualMetric.a1;
        a2 = rActualMetric.a2;

        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        a1xa2 = rActualMetric.a3_tilde;
        hact = sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]);

        for (m = 0; m < 3; m++)
            a3_KL[m] = a1xa2[m] / hact * thickness / 2.0;

        /* Schubdifferenzvektor: w = w^1*a1 + w^2*a2 */
        array_1d<double, 3> w;
        double w1, w2;
        w = ZeroVector(3);
        w1 = 0.0;
        w2 = 0.0;

        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
        for (m = 0; m < num_node; m++)
        {
            w1   += funct[m] * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            w2   += funct[m] * GetGeometry()[m].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }  
        
        w = w1 * a1 + w2 * a2;

        /* a3 nach hier. 5p-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
        a3 = a3_KL + w;

        /* Ableitungen von a1 und a2 nach der Richtung alpha */
        da1_dT1 = rActualMetric.Da1_D1;
        da2_dT2 = rActualMetric.Da2_D2;
        da1_dT2 = rActualMetric.Da1_D2;
        da2_dT1 = da1_dT2;

        /* Ableitung von a3_KL nach der Richtung alpha */

        /* Ableitung des Zaehlers von a3_KL nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        d_a1crossa2_dT1 = CrossProduct(da1_dT1, a2) + CrossProduct(a1, da2_dT1);

        /* Ableitung des Zaehlers von a3_KL nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        d_a1crossa2_dT2 = CrossProduct(da1_dT2, a2) + CrossProduct(a1, da2_dT2);

        /* Ableitung des Nenners von a3_KL */
        deriv_sqrt_norma1crossa2 = 1.0 / (2.0 * sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]));

        /* Ableitung des Nenners von a3_KL nach 1 */
        deriv_radicant_norma1crossa2_dT1 = 2.0 * d_a1crossa2_dT1[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT1[1] * a1xa2[1] + 2.0 * d_a1crossa2_dT1[2] * a1xa2[2];
        dnorma1crossa2_dT1 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT1;

        /* Ableitung des Nenners von a3_KL nach 2 */
        deriv_radicant_norma1crossa2_dT2 = 2.0 * d_a1crossa2_dT2[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT2[1] * a1xa2[1] + 2.0 * d_a1crossa2_dT2[2] * a1xa2[2];
        dnorma1crossa2_dT2 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT2;

        /* Ableitung von a3_KL mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_KL_dT1[m] = (d_a1crossa2_dT1[m] * hact - a1xa2[m] * dnorma1crossa2_dT1) / (hact * hact) * thickness / 2.0;
            da3_KL_dT2[m] = (d_a1crossa2_dT2[m] * hact - a1xa2[m] * dnorma1crossa2_dT2) / (hact * hact) * thickness / 2.0;
        }

        /* Ableitung des Schubdifferenzvektors nach der Richtung alpha */
        array_1d<double, 3> dw_dT1;
        array_1d<double, 3> dw_dT2;
        double dw1_dT1, dw2_dT1, dw1_dT2, dw2_dT2;
        dw_dT1 = ZeroVector(3);
        dw_dT2 = ZeroVector(3);
        dw1_dT1 = 0.0;
        dw1_dT2 = 0.0;
        dw2_dT1 = 0.0;
        dw2_dT2 = 0.0;

        for (m = 0; m < num_node; m++)
        {
            dw1_dT1 += deriv(m, 0) * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw1_dT2 += deriv(m, 1) * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();

            dw2_dT1 += deriv(m, 0) * GetGeometry()[m].GetDof(ROTATION_X, pos + 4).GetSolutionStepValue();
            dw2_dT2 += deriv(m, 1) * GetGeometry()[m].GetDof(ROTATION_X, pos + 4).GetSolutionStepValue();
        }

        dw_dT1 = dw1_dT1 * a1 + w1 * da1_dT1 + dw2_dT1 * a2 + w2 * da2_dT1;
        dw_dT2 = dw1_dT2 * a1 + w1 * da1_dT2 + dw2_dT2 * a2 + w2 * da2_dT2;


        /* - INTEGRAND DER GEOMETRISCHEN STEIFIGKEITSMATRIX - */

        array_1d<double, 3> da1xa2k;
        array_1d<double, 3> da1xa2l;
        array_1d<double, 3> dda1xa2kl;
        array_1d<double, 3> da1k;
        array_1d<double, 3> da2k;
        array_1d<double, 3> da3k;
        array_1d<double, 3> da1l;
        array_1d<double, 3> da2l;
        array_1d<double, 3> da3l;
        array_1d<double, 3> da11k;
        array_1d<double, 3> da12k;
        array_1d<double, 3> da21k;
        array_1d<double, 3> da22k;
        array_1d<double, 3> da11l;
        array_1d<double, 3> da12l;
        array_1d<double, 3> da21l;
        array_1d<double, 3> da22l;
        array_1d<double, 3> dda3kl;

        array_1d<double, 3> dw_dk;
        array_1d<double, 3> dw_dl;
        array_1d<double, 3> ddw_dkl;
        array_1d<double, 3> dw_dT1_dk;
        array_1d<double, 3> dw_dT2_dk;
        array_1d<double, 3> dw_dT1_dl;
        array_1d<double, 3> dw_dT2_dl;
        array_1d<double, 3> ddw_dT1_dkl;
        array_1d<double, 3> ddw_dT2_dkl;
        double dw1_dk, dw2_dk, dw1_dl, dw2_dl;
        double dw1_dT1_dk, dw2_dT1_dk, dw1_dT2_dk, dw2_dT2_dk,
                dw1_dT1_dl, dw2_dT1_dl, dw1_dT2_dl, dw2_dT2_dl;


        dw_dk = ZeroVector(3);
        dw_dl = ZeroVector(3);
        ddw_dkl = ZeroVector(3);
        dw_dT1_dk = ZeroVector(3);
        dw_dT2_dk = ZeroVector(3);
        dw_dT1_dl = ZeroVector(3);
        dw_dT2_dl = ZeroVector(3);
        ddw_dT1_dkl = ZeroVector(3);
        ddw_dT2_dkl = ZeroVector(3);

        da1xa2k = ZeroVector(3);
        da1xa2l = ZeroVector(3);
        dda1xa2kl = ZeroVector(3);
        da3k = ZeroVector(3);
        da3l = ZeroVector(3);
        dda3kl = ZeroVector(3);

        int i, j;
        double da1xa2k_a1xa2;
        double da1xa2k_a1xa2_hact;
        double da1xa2l_a1xa2;
        double da1xa2l_a1xa2_hact;
        double C, D;

        Matrix dE11_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E11 nach FHG abgeleitet
        Matrix dE12_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E12 nach FHG abgeleitet
        Matrix dE22_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E22 nach FHG abgeleitet
        Matrix dE13_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E13 nach FHG abgeleitet
        Matrix dE23_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E23 nach FHG abgeleitet

        for (k = 0; k < num_node; k++)
        {
            for (i = 0; i < dof_per_node; i++)
            {
            for (l = 0; l < num_node; l++)
            {
                for (j = 0; j < dof_per_node; j++)
                {
                /* Anteile aus 3p-Formulierung (nicht veraendert) */
                if (i < 3 && j < 3)
                {
                if (i == 0)
                {
                    da1k[0] = deriv(k, 0);
                    da1k[1] = 0;
                    da1k[2] = 0;
                    da2k[0] =deriv(k, 1);
                    da2k[1] = 0;
                    da2k[2] = 0;
                }
                else if (i == 1)
                {
                    da1k[0] = 0;
                    da1k[1] = deriv(k, 0);
                    da1k[2] = 0;
                    da2k[0] = 0;
                    da2k[1] = deriv(k, 1);
                    da2k[2] = 0;
                }
                else
                {
                    da1k[0] = 0;
                    da1k[1] = 0;
                    da1k[2] = deriv(k, 0);
                    da2k[0] = 0;
                    da2k[1] = 0;
                    da2k[2] = deriv(k, 1);
                }

                da1xa2k = CrossProduct(da1k, a2) + CrossProduct(a1, da2k);
                da1xa2k_a1xa2 = (da1xa2k[0] * a1xa2[0] + da1xa2k[1] * a1xa2[1] + da1xa2k[2] * a1xa2[2]);
                da1xa2k_a1xa2_hact = da1xa2k_a1xa2 / (hact * hact * hact);
                da3k[0] = da1xa2k[0] / hact - a1xa2[0] * da1xa2k_a1xa2_hact;
                da3k[1] = da1xa2k[1] / hact - a1xa2[1] * da1xa2k_a1xa2_hact;
                da3k[2] = da1xa2k[2] / hact - a1xa2[2] * da1xa2k_a1xa2_hact;

                if (j == 0)
                {
                    da1l[0] = deriv(l, 0);
                    da1l[1] = 0;
                    da1l[2] = 0;
                    da2l[0] = deriv(l, 1);
                    da2l[1] = 0;
                    da2l[2] = 0;
                }
                else if (j == 1)
                {
                    da1l[0] = 0;
                    da1l[1] = deriv(l, 0);
                    da1l[2] = 0;
                    da2l[0] = 0;
                    da2l[1] = deriv(l, 1);
                    da2l[2] = 0;
                }
                else
                {
                    da1l[0] = 0;
                    da1l[1] = 0;
                    da1l[2] = deriv(l, 0);
                    da2l[0] = 0;
                    da2l[1] = 0;
                    da2l[2] = deriv(l, 1);
                }

                da1xa2l = CrossProduct(da1l, a2) + CrossProduct(a1, da2l);
                da1xa2l_a1xa2 = (da1xa2l[0] * a1xa2[0] + da1xa2l[1] * a1xa2[1] + da1xa2l[2] * a1xa2[2]);
                da1xa2l_a1xa2_hact = da1xa2l_a1xa2 / (hact * hact * hact);
                da3l[0] = da1xa2l[0] / hact - a1xa2[0] * da1xa2l_a1xa2_hact;
                da3l[1] = da1xa2l[1] / hact - a1xa2[1] * da1xa2l_a1xa2_hact;
                da3l[2] = da1xa2l[2] / hact - a1xa2[2] * da1xa2l_a1xa2_hact;

                if (i == 0)
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) *deriv(k, 1);
                    }
                    else if (j == 2)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) *deriv(k, 1);
                    dda1xa2kl[2] = 0;
                    }
                }

                else if (i == 1)
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) *deriv(k, 1);
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 2)
                    {
                    dda1xa2kl[0] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) *deriv(k, 1);
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                }

                else
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) *deriv(k, 1);
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) *deriv(k, 1);
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                }

                C = -(dda1xa2kl[0] * a1xa2[0] + dda1xa2kl[1] * a1xa2[1]
                    + dda1xa2kl[2] * a1xa2[2] + da1xa2k[0] * da1xa2l[0]
                    + da1xa2k[1] * da1xa2l[1] + da1xa2k[2] * da1xa2l[2])
                    / (hact * hact * hact);

                D = 3.0 * (da1xa2k_a1xa2) * (da1xa2l_a1xa2) / (hact * hact * hact * hact * hact);

                dda3kl[0] = dda1xa2kl[0] / hact - da1xa2l_a1xa2_hact * da1xa2k[0]
                            - da1xa2k_a1xa2_hact * da1xa2l[0] + C * a1xa2[0] + D * a1xa2[0];
                dda3kl[1] = dda1xa2kl[1] / hact - da1xa2l_a1xa2_hact * da1xa2k[1]
                            - da1xa2k_a1xa2_hact * da1xa2l[1] + C * a1xa2[1] + D * a1xa2[1];
                dda3kl[2] = dda1xa2kl[2] / hact - da1xa2l_a1xa2_hact * da1xa2k[2]
                            - da1xa2k_a1xa2_hact * da1xa2l[2] + C * a1xa2[2] + D * a1xa2[2];

                /* Membrananteil */
                if (i == j)
                {
                    dE11_dkl(5 * k + i, 5 * l + j) = deriv(k, 0) * deriv(l, 0);
                    dE12_dkl(5 * k + i, 5 * l + j) = 0.5 * (deriv(k, 0) * deriv(l, 1) + deriv(l, 0) *deriv(k, 1));
                    dE22_dkl(5 * k + i, 5 * l + j) =deriv(k, 1) * deriv(l, 1);
                }
                else
                {
                    dE11_dkl(5 * k + i, 5 * l + j) = 0;
                    dE12_dkl(5 * k + i, 5 * l + j) = 0;
                    dE22_dkl(5 * k + i, 5 * l + j) = 0;
                }

                /* Kruemmungsanteil */
                dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 0) * da3l[i] + s_deriv(l, 0) * da3k[j]
                            + da1_dT1[0] * dda3kl[0] + da1_dT1[1] * dda3kl[1]
                            + da1_dT1[2] * dda3kl[2]);
                dE12_dkl(5 * k + i, 5 * l + j) = dE12_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 2) * da3l[i] + s_deriv(l, 2) * da3k[j]
                            + da1_dT2[0] * dda3kl[0] + da1_dT2[1] * dda3kl[1]
                            + da1_dT2[2] * dda3kl[2]);
                dE22_dkl(5 * k + i, 5 * l + j) = dE22_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 1) * da3l[i] + s_deriv(l, 1) * da3k[j]
                            + da2_dT2[0] * dda3kl[0] + da2_dT2[1] * dda3kl[1]
                            + da2_dT2[2] * dda3kl[2]);
                }

                else
                {
                /* Anteile aus 5p-Formulierung  */

                /* initialisieren */
                dw1_dk = 0.0;
                dw1_dT1_dk = 0.0;
                dw1_dT2_dk = 0.0;

                dw2_dk = 0.0;
                dw2_dT1_dk = 0.0;
                dw2_dT2_dk = 0.0;

                /* Ableitung nach dem ersten FHG k */
                if (i == 0)
                {
                    da1k[0] = deriv(k, 0);
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] =deriv(k, 1);
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = s_deriv(k, 0);
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = s_deriv(k, 2);
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = s_deriv(k, 1);
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 1)
                {
                    da1k[0] = 0.0;
                    da1k[1] = deriv(k, 0);
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] =deriv(k, 1);
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = s_deriv(k, 0);
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = s_deriv(k, 2);
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = s_deriv(k, 1);
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 2)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = deriv(k, 0);
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] =deriv(k, 1);
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = s_deriv(k, 0);
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = s_deriv(k, 2);
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = s_deriv(k, 1);
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 3)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = funct[k];
                    dw1_dT1_dk = deriv(k, 0);
                    dw1_dT2_dk =deriv(k, 1);

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 4)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = funct[k];
                    dw2_dT1_dk = deriv(k, 0);
                    dw2_dT2_dk =deriv(k, 1);
                }
                

                dw_dk = dw1_dk*a1 + w1*da1k + dw2_dk*a2 + w2*da2k;
                dw_dT1_dk = dw1_dT1_dk*a1 + dw1_dT1*da1k + dw1_dk*da1_dT1 + w1*da11k
                            + dw2_dT1_dk*a2 + dw2_dT1*da2k + dw2_dk*da2_dT1 + w2*da21k;
                dw_dT2_dk = dw1_dT2_dk*a1 + dw1_dT2*da1k + dw1_dk*da1_dT2 + w1*da12k
                            + dw2_dT2_dk*a2 + dw2_dT2*da2k + dw2_dk*da2_dT2 + w2*da22k;

                /* initialisieren */
                dw1_dl = 0.0;
                dw1_dT1_dl = 0.0;
                dw1_dT2_dl = 0.0;

                dw2_dl = 0.0;
                dw2_dT1_dl = 0.0;
                dw2_dT2_dl = 0.0;

                /* Ableitung nach dem zweiten FHG l */
                if (j == 0)
                {
                    da1l[0] = deriv(l, 0);
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = deriv(l, 1);
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = s_deriv(l, 0);
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = s_deriv(l, 2);
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = s_deriv(l, 1);
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 1)
                {
                    da1l[0] = 0.0;
                    da1l[1] = deriv(l, 0);
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = deriv(l, 1);
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = s_deriv(l, 0);
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = s_deriv(l, 2);
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = s_deriv(l, 1);
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 2)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = deriv(l, 0);
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = deriv(l, 1);
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = s_deriv(l, 0);
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = s_deriv(l, 2);
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = s_deriv(l, 1);
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 3)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = funct[l];
                    dw1_dT1_dl = deriv(l, 0);
                    dw1_dT2_dl = deriv(l, 1);

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 4)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = funct[l];
                    dw2_dT1_dl = deriv(l, 0);
                    dw2_dT2_dl = deriv(l, 1);
                }

                dw_dl = dw1_dl*a1 + w1*da1l + dw2_dl*a2 + w2*da2l;
                dw_dT1_dl = dw1_dT1_dl*a1 + dw1_dT1*da1l + dw1_dl*da1_dT1 + w1*da11l
                            + dw2_dT1_dl*a2 + dw2_dT1*da2l + dw2_dl*da2_dT1 + w2*da21l ;
                dw_dT2_dl = dw1_dT2_dl*a1 + dw1_dT2*da1l + dw1_dl*da1_dT2 + w1*da12l
                            + dw2_dT2_dl*a2 + dw2_dT2*da2l + dw2_dl*da2_dT2 + w2*da22l ;


                /* Ableitung nach beiden FHG kl */
                ddw_dkl = dw1_dk*da1l + dw1_dl*da1k + dw2_dk*da2l + dw2_dl*da2k;
                ddw_dT1_dkl = dw1_dT1_dk*da1l + dw1_dT1_dl*da1k + dw1_dk*da11l + dw1_dl*da11k
                            + dw2_dT1_dk*da2l + dw2_dT1_dl*da2k + dw2_dk*da21l + dw2_dl*da21k;
                ddw_dT2_dkl = dw1_dT2_dk*da1l + dw1_dT2_dl*da1k + dw1_dk*da12l + dw1_dl*da12k
                            + dw2_dT2_dk*da2l + dw2_dT2_dl*da2k + dw2_dk*da22l + dw2_dl*da22k;

                /* Anteile fuer Querschubverzerrungen (nur konstant in Dickenrichtung) */
                dE13_dkl(5 * k + i, 5 * l + j) = dE13_dkl(5 * k + i, 5 * l + j)
                                                + 0.5*(inner_prod(da1k, dw_dl) + inner_prod(da1l, dw_dk) + inner_prod(a1, ddw_dkl));
                dE23_dkl(5 * k + i, 5 * l + j) = dE23_dkl(5 * k + i, 5 * l + j)
                                                + 0.5*(inner_prod(da2k, dw_dl) + inner_prod(da2l, dw_dk) + inner_prod(a2, ddw_dkl));


                // #if 0
                //         /* Anteile fuer Querschubverzerrungen (LINEARE Zusatzterme in Dickenrichtung theta^3) */
                //             dE13_dkl(5 * k + i, 5 * l + j) = dE13_dkl(5 * k + i, 5 * l + j)
                //                                             + 0.5*mZeta*(ddw_dT1_dkl*w + dw_dT1_dk*dw_dl + dw_dT1_dl*dw_dk + dw_dT1*ddw_dkl);
                //             dE23_dkl(5 * k + i, 5 * l + j) = dE23_dkl(5 * k + i, 5 * l + j)
                //                                             + 0.5*mZeta*(ddw_dT2_dkl*w + dw_dT2_dk*dw_dl + dw_dT2_dl*dw_dk + dw_dT2*ddw_dkl);
                // #endif


                /* Anteile fuer Inplane-Verzerrungen (nur linear in Dickenrichtung theta^3) */
                dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i, 5 * l + j)
                                                + thickness/2.0*mZeta*0.5*(inner_prod(da1k, dw_dT1_dl) + inner_prod(da1l, dw_dT1_dk) + inner_prod(a1, ddw_dT1_dkl)
                                                        + inner_prod(da1k, dw_dT1_dl) + inner_prod(da1l, dw_dT1_dk) + inner_prod(a1, ddw_dT1_dkl));     // MLt
                dE12_dkl(5 * k + i, 5 * l + j) = dE12_dkl(5 * k + i, 5 * l + j)
                                                + thickness/2.0*mZeta*0.5*(inner_prod(da1k, dw_dT2_dl) + inner_prod(da1l, dw_dT2_dk) + inner_prod(a1, ddw_dT2_dkl)
                                                        + inner_prod(da2k, dw_dT1_dl) + inner_prod(da2l, dw_dT1_dk) + inner_prod(a2, ddw_dT1_dkl));     // MLt
                dE22_dkl(5 * k + i, 5 * l + j) = dE22_dkl(5 * k + i, 5 * l + j)
                                                + thickness/2.0*mZeta*0.5*(inner_prod(da2k, dw_dT2_dl) + inner_prod(da2l, dw_dT2_dk) + inner_prod(a2, ddw_dT2_dkl)
                                                        + inner_prod(da2k, dw_dT2_dl) + inner_prod(da2l, dw_dT2_dk) + inner_prod(a2, ddw_dT2_dkl));     // MLt
                }

                double dE11_dkl_car = 0.0;
                double dE22_dkl_car = 0.0;
                double dE12_dkl_car = 0.0;
                double dE23_dkl_car = 0.0;
                double dE13_dkl_car = 0.0;
                
                dE11_dkl_car += mInitialMetric.Q(0, 0) * dE11_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(0, 1) * dE22_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(0, 2) * dE12_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(0, 3) * dE23_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(0, 4) * dE13_dkl(5 * k + i, 5 * l + j);
                dE22_dkl_car += mInitialMetric.Q(1, 0) * dE11_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(1, 1) * dE22_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(1, 2) * dE12_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(1, 3) * dE23_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(1, 4) * dE13_dkl(5 * k + i, 5 * l + j);
                dE12_dkl_car += mInitialMetric.Q(2, 0) * dE11_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(2, 1) * dE22_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(2, 2) * dE12_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(2, 3) * dE23_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(2, 4) * dE13_dkl(5 * k + i, 5 * l + j);
                dE23_dkl_car += mInitialMetric.Q(3, 0) * dE11_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(3, 1) * dE22_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(3, 2) * dE12_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(3, 3) * dE23_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(3, 4) * dE13_dkl(5 * k + i, 5 * l + j);
                dE13_dkl_car += mInitialMetric.Q(4, 0) * dE11_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(4, 1) * dE22_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(4, 2) * dE12_dkl(5 * k + i, 5 * l + j) + mInitialMetric.Q(4, 3) * dE23_dkl(5 * k + i, 5 * l + j) + 
                    mInitialMetric.Q(4, 4) * dE13_dkl(5 * k + i, 5 * l + j);                
                
                /* Integrand von k_g */
                IKg(5 * k + i, 5 * l + j) =   dE11_dkl_car * S[0]
                                            + dE22_dkl_car * S[1]
                                            +   dE12_dkl_car * S[2]
                                            + dE23_dkl_car * S[3]
                                            + dE13_dkl_car * S[4];


                }
            }
        }
     }
    }

    void IgaShell5pElementStuttgart::Calculate(
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

            array_1d<double, 5> strain_vector;
            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, strain_vector, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

            array_1d<double, 3> G1 = ZeroVector(3);
            array_1d<double, 3> G2 = ZeroVector(3);
            array_1d<double, 3> g1 = ZeroVector(3);
            array_1d<double, 3> g2 = ZeroVector(3);
            array_1d<double, 3> g3 = ZeroVector(3);
            Matrix F = ZeroMatrix(3, 3);
            double detF = 0.0;
            CalculateInitialBaseVectorsGLinearized(G1, G2);
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
            KRATOS_WATCH("No results for desired variable available in Calculate of IgaShell5pElementStuttgart.")
            rValues = 0.0;
        }
    }

    void IgaShell5pElementStuttgart::EquationIdVector(
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

    void IgaShell5pElementStuttgart::GetDofList(
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

    int IgaShell5pElementStuttgart::Check(const ProcessInfo& rCurrentProcessInfo)
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
