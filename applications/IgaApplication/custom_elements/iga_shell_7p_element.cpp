//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Michael Loibl (michael.loibl@unibw.de)
//                   based on the work of Ralf Echter (7p linear), Bastian Oesterle (5p non-linear) and 
//                      A. Usta (7p non-linear, master thesis at the University of Stuttgart) on a 
//                      hierarchic shell formulation
/**
    * Description:
    * Isogeometric hierarchic Reissner-Mindlin shell element parameterized by 7 parameters (7p)
    * The 7 parameters: three translations (u,v,w), three components of the shear difference vector 
    *   (w_1,w_2, w_3) and w_bar
    * The 7p element takes shear deformations and thickness changes into account according to the 3D 
    *   shell theory.
    * The hierarchy means that the shell element builds on top of the Reissner-Mindlin 
    *   (iga_shell_5p_element) and the Kirchhoff-Love shell element (iga_shell_3p_element) in the sense 
    *   that the existing equations remain unchanged and additional terms are added. This remains a 
    *   theoretical idea in the case of the 7p shell and is not reflected in the concrete implementation.
    * 
    * Implementation:
    * So far there are no new parameters defined, instead ROTATION_X = w_1 and ROTATION_Y = w_2.
    * The implementation was not revised w.r.t. style, comments and naming.
    * 
    * Attention:
    * In the last version of this element, there were still open questions and problems. For further 
    *   information please do not hesitate to contact the author.
*/


// System includes
#include "utilities/math_utils.h"

// External includes
#include <iostream>     // MLout
#include <fstream>      // MLout

// Project includes
#include "custom_elements/iga_shell_7p_element.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell7pElement::Initialize()
    {
        KRATOS_TRY

        // Constitutive Law initialisation
        BaseDiscreteElement::Initialize();

        CalculateMetric(mInitialMetric);
        
        mZeta = 0.0;
        mInitialTransConToCar = ZeroMatrix(6, 6);

        KRATOS_CATCH("")
    }

    void IgaShell7pElement::CalculateAll(
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
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 7;

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

        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);  // FML: test ab initio initialization of shape function vector N

        // shear difference vector
        array_1d<double, 3> w = ZeroVector(3);
        // derivatives of the shear difference vector
        array_1d<double, 3> Dw_D1 = ZeroVector(3);
        array_1d<double, 3> Dw_D2 = ZeroVector(3);
        // components w_i of the shear difference vector which calculates as (w_i(1) * e1 + w_i(2) * e2 + w_i(3) * e3)
        array_1d<double, 3> w_i = ZeroVector(3);
        // derivatives of the components w_i
        Matrix Dw_i_Dalpha = ZeroMatrix(3, 2);
        // 7th parameter
        double w_bar = 0.0;
        double dV = 0.0;

        MetricVariables actual_metric(3, 6);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_i, Dw_i_Dalpha, actual_metric);
        for (unsigned int i = 0; i < number_of_control_points; ++i){
            w_bar += N(i) * GetGeometry()[i].GetDof(W_BAR).GetSolutionStepValue();
        }

        double thickness = GetProperties().GetValue(THICKNESS);

        for (unsigned int Gauss_index = 0; Gauss_index < mGaussQuadratureThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussQuadratureThickness.zeta(Gauss_index);
            
            // Differential Volume
            array_1d<double, 3> G1(3, 0.0), G2(3, 0.0), G1_con(3, 0.0), G2_con(3, 0.0);
            CalculateInitialBaseVectorsLinearised(G1, G2, G1_con, G2_con);
            dV = inner_prod(MathUtils<double>::CrossProduct(G1, G2), mInitialMetric.a3_KL);   // mInitialMetric.a3_KL = G3
            
            // TransformationMatrix
            CalculateTransformationMatrixInitialTransConToCar(G1_con, G2_con);

            Matrix B = ZeroMatrix(6, mat_size);
            SecondVariations second_variations(mat_size);
            ConstitutiveVariables constitutive_variables(6);
            CalculateConstitutiveVariables(actual_metric, w, Dw_D1, Dw_D2, w_bar, constitutive_variables, Values, 
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            CalculateB(B, actual_metric);
            CalculateVariations7p(B, second_variations, w, Dw_D1, Dw_D2, w_i, Dw_i_Dalpha, w_bar,
                actual_metric, CalculateStiffnessMatrixFlag);
            
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
            if(Id()==1){
                KRATOS_WATCH(B)
                KRATOS_WATCH(constitutive_variables.D)
                KRATOS_WATCH(integration_weight)
                KRATOS_WATCH(constitutive_variables.S)
            }
            
            // MLout
            if(Id()==1 && mcount==0 && mZeta < 0){
                std::ofstream myfile;
                myfile.open("B-7p.csv");
                for(unsigned int i=0; i<6;i++){
                    for(unsigned int j=0; j<mat_size;j++){
                        myfile<<B(i,j)<<",";
                    }
                    myfile << "\n";
                }
                myfile.close();
            }
        } // GP,t loop
        if(Id()==1 && rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
            KRATOS_WATCH(rCurrentProcessInfo.GetValue(NL_ITERATION_NUMBER))

        KRATOS_CATCH("");
    }

    void IgaShell7pElement::CalculateAndAddKm(
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

    void IgaShell7pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 7;

        if (SD.size() != 6){
            KRATOS_WATCH("Stress size is wrong.")
            KRATOS_ERROR << "Stress size is wrong." << std::endl;
        }

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

    void IgaShell7pElement::CalculateMetric(
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
        rMetric.gab[0] = pow(rMetric.a1[0], 2) + pow(rMetric.a1[1], 2) + pow(rMetric.a1[2], 2);
        rMetric.gab[1] = pow(rMetric.a2[0], 2) + pow(rMetric.a2[1], 2) + pow(rMetric.a2[2], 2);
        rMetric.gab[2] = rMetric.a1[0] * rMetric.a2[0] + rMetric.a1[1] * rMetric.a2[1] + rMetric.a1[2] * rMetric.a2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);
        
        // derivatives of a1 and a2 w.r.t. curvilinear coordinates
        for (unsigned int i = 0; i < 3; i++)
        {
            rMetric.Da1_D1[i] = mInitialMetric.H(i, 0);
            rMetric.Da2_D2[i] = mInitialMetric.H(i, 1);
            rMetric.Da1_D2[i] = mInitialMetric.H(i, 2);
        }
        
        // derivatives of a3_KL w.r.t. curvilinear coordinates
        double DdA_D1, DdA_D2;
        array_1d<double, 3> Da1_D1xa2, a1xDa2_D1, Da1_D2xa2, a1xDa2_D2;
        MathUtils<double>::CrossProduct(Da1_D1xa2, rMetric.Da1_D1, rMetric.a2);
        MathUtils<double>::CrossProduct(a1xDa2_D1, rMetric.a1, rMetric.Da1_D2);
        rMetric.Da3_KL_tilde_D1 = Da1_D1xa2 + a1xDa2_D1;
        MathUtils<double>::CrossProduct(Da1_D2xa2, rMetric.Da1_D2, rMetric.a2);
        MathUtils<double>::CrossProduct(a1xDa2_D2, rMetric.a1, rMetric.Da2_D2);
        rMetric.Da3_KL_tilde_D2 = Da1_D2xa2 + a1xDa2_D2;
        DdA_D1 = (rMetric.Da3_KL_tilde_D1[0] * rMetric.a3_KL_tilde[0] 
            + rMetric.Da3_KL_tilde_D1[1] * rMetric.a3_KL_tilde[1] + rMetric.Da3_KL_tilde_D1[2] * rMetric.a3_KL_tilde[2]) / rMetric.dA;
        DdA_D2 = (rMetric.Da3_KL_tilde_D2[0] * rMetric.a3_KL_tilde[0] 
            + rMetric.Da3_KL_tilde_D2[1] * rMetric.a3_KL_tilde[1] + rMetric.Da3_KL_tilde_D2[2] * rMetric.a3_KL_tilde[2]) / rMetric.dA;
        rMetric.Da3_KL_D1 = (rMetric.Da3_KL_tilde_D1 * rMetric.dA - rMetric.a3_KL_tilde * DdA_D1) / 
            (rMetric.dA * rMetric.dA);
        rMetric.Da3_KL_D2 = (rMetric.Da3_KL_tilde_D2 * rMetric.dA - rMetric.a3_KL_tilde * DdA_D2) / 
            (rMetric.dA * rMetric.dA);

        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.a3_KL[0] + rMetric.H(1, 0) * rMetric.a3_KL[1] + rMetric.H(2, 0) * rMetric.a3_KL[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.a3_KL[0] + rMetric.H(1, 1) * rMetric.a3_KL[1] + rMetric.H(2, 1) * rMetric.a3_KL[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.a3_KL[0] + rMetric.H(1, 2) * rMetric.a3_KL[1] + rMetric.H(2, 2) * rMetric.a3_KL[2];

        //contravariant rMetric a_ab_con and base vectors g_con
        //Vector a_ab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.a_ab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.a_ab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.a_ab_con[1] = invdetGab*rMetric.gab[0];
        // contravariant base vectors a_con
        rMetric.a1_con = rMetric.a1*rMetric.a_ab_con[0] + rMetric.a2*rMetric.a_ab_con[2];
        rMetric.a2_con = rMetric.a1*rMetric.a_ab_con[2] + rMetric.a2*rMetric.a_ab_con[1];
        // g_con_3 = a3_KL
        
        //local cartesian coordinates
        double lg1 = norm_2(rMetric.a1);
        array_1d<double, 3> e1 = rMetric.a1 / lg1;
        double lg_con2 = norm_2(rMetric.a2_con);
        array_1d<double, 3> e2 = rMetric.a2_con / lg_con2;
        // e3 = a3_KL = g_con_3

        // // transformation matrix T_con_to_car from curvilinear to local cartesian coordinate system
        // faster computation of transformation matrix T_con_to_car taking into account that a lot of entries become zero (ML)
        // this matrix T_con_to_car is referring to a VoigtSize 6 with E11, E22, E33, E12, E23, E13
        double mG_00 = inner_prod(e1, rMetric.a1_con);
        double mG_10 = inner_prod(e2, rMetric.a1_con);
        double mG_11 = inner_prod(e2, rMetric.a2_con);
        
        rMetric.T_con_to_car(0, 0) = pow(mG_00, 2);
        rMetric.T_con_to_car(1, 0) = pow(mG_10, 2);
        rMetric.T_con_to_car(1, 1) = pow(mG_11, 2);
        rMetric.T_con_to_car(1, 3) = 2.00 * mG_10 * mG_11;
        rMetric.T_con_to_car(2, 2) = 1.00;
        rMetric.T_con_to_car(3, 0) = 2.00 * mG_00 * mG_10;
        rMetric.T_con_to_car(3, 3) = 2.00 * mG_00 * mG_11;
        rMetric.T_con_to_car(4, 4) = 2.00 * mG_11;
        rMetric.T_con_to_car(4, 5) = 2.00 * mG_10;
        rMetric.T_con_to_car(5, 5) = 2.00 * mG_00;

        // transformation matrix T_car_to_cov from local Cartesian to covariant basis
        rMetric.T_car_to_cov(0, 0) = pow(mG_00, 2);
        rMetric.T_car_to_cov(0, 1) = pow(mG_10, 2);
        rMetric.T_car_to_cov(0, 3) = 2.0 * mG_00 * mG_10;
        rMetric.T_car_to_cov(1, 1) = pow(mG_11, 2);
        rMetric.T_car_to_cov(2, 2) = 1.00;
        rMetric.T_car_to_cov(3, 1) = mG_10 * mG_11;
        rMetric.T_car_to_cov(3, 3) = mG_00 * mG_11;
        rMetric.T_car_to_cov(4, 4) = mG_11;
        rMetric.T_car_to_cov(5, 4) = mG_10;
        rMetric.T_car_to_cov(5, 5) = mG_00;

        mG_00 = inner_prod(e1, rMetric.a1);
        double mG_01 = inner_prod(e1, rMetric.a2);
        mG_11 = inner_prod(e2, rMetric.a2);
        // transformation matrix T_cov_to_car from covariant to local Cartesian basis
        rMetric.T_cov_to_car(0, 0) = pow(mG_00, 2);
        rMetric.T_cov_to_car(0, 1) = pow(mG_01, 2);
        rMetric.T_cov_to_car(0, 2) = 2 * mG_00 * mG_01;
        rMetric.T_cov_to_car(1, 1) = pow(mG_11, 2);
        rMetric.T_cov_to_car(2, 2) = 1.0;
        rMetric.T_cov_to_car(3, 1) = mG_01 * mG_11;
        rMetric.T_cov_to_car(3, 3) = mG_00 * mG_11;
        rMetric.T_cov_to_car(4, 4) = mG_11;
        rMetric.T_cov_to_car(5, 4) = mG_01;
        rMetric.T_cov_to_car(5, 5) = mG_00;
    }

    void IgaShell7pElement::CalculateTransformationMatrixInitialTransConToCar(
        const array_1d<double, 3>& rG1_con, 
        const array_1d<double, 3>& rG2_con)
    {
        //local cartesian coordinates
        array_1d<double, 3> e1 = mInitialMetric.a1 / norm_2(mInitialMetric.a1);
        array_1d<double, 3> e2 = mInitialMetric.a2_con / norm_2(mInitialMetric.a2_con);
        // e3 = a3_KL = g_con_3

        double mG_00 = inner_prod(e1, rG1_con);
        double mG_10 = inner_prod(e2, rG1_con);
        double mG_11 = inner_prod(e2, rG2_con);
        
        mInitialTransConToCar = ZeroMatrix(6, 6);
        mInitialTransConToCar(0, 0) = pow(mG_00, 2);
        mInitialTransConToCar(1, 0) = pow(mG_10, 2);
        mInitialTransConToCar(1, 1) = pow(mG_11, 2);
        mInitialTransConToCar(1, 3) = 2.00 * mG_10 * mG_11;
        mInitialTransConToCar(2, 2) = 1.00;
        mInitialTransConToCar(3, 0) = 2.00 * mG_00 * mG_10;
        mInitialTransConToCar(3, 3) = 2.00 * mG_00 * mG_11;
        mInitialTransConToCar(4, 4) = 2.00 * mG_11;
        mInitialTransConToCar(4, 5) = 2.00 * mG_10;
        mInitialTransConToCar(5, 5) = 2.00 * mG_00;
    }

    void IgaShell7pElement::CalculateTransformationMatrixInitialTransCarToCov(
        Matrix& rInitialTransCarToCov,
        const array_1d<double, 3>& rG1_con, 
        const array_1d<double, 3>& rG2_con)
    {
        //local cartesian coordinates
        array_1d<double, 3> e1 = mInitialMetric.a1 / norm_2(mInitialMetric.a1);
        array_1d<double, 3> e2 = mInitialMetric.a2_con / norm_2(mInitialMetric.a2_con);
        // e3 = a3_KL = g_con_3

        double mG_00 = inner_prod(e1, rG1_con);
        double mG_10 = inner_prod(e2, rG1_con);
        double mG_11 = inner_prod(e2, rG2_con);

        // transformation matrix T_car_to_cov from local Cartesian to covariant basis
        rInitialTransCarToCov(0, 0) = pow(mG_00, 2);
        rInitialTransCarToCov(0, 1) = pow(mG_10, 2);
        rInitialTransCarToCov(0, 3) = 2.0 * mG_00 * mG_10;
        rInitialTransCarToCov(1, 1) = pow(mG_11, 2);
        rInitialTransCarToCov(2, 2) = 1.00;
        rInitialTransCarToCov(3, 1) = mG_10 * mG_11;
        rInitialTransCarToCov(3, 3) = mG_00 * mG_11;
        rInitialTransCarToCov(4, 4) = mG_11;
        rInitialTransCarToCov(5, 4) = mG_10;
        rInitialTransCarToCov(5, 5) = mG_00;
    }

    void IgaShell7pElement::CalculateTransformationMatrixActualTransCovToCar(
        Matrix& rActualTransCovToCar,
        const Vector& rg1,
        const Vector& rg2,
        const Vector& rg3,
        const Vector& ra2_con,
        const Vector& ra3_KL)
    {
        // local cartesian coordinates, direction of a2_con and g2_con is the same -> it is possible to use ra2_con instead of g2_con 
        array_1d<double, 3> e1 = rg1 / norm_2(rg1);
        array_1d<double, 3> e2 = ra2_con / norm_2(ra2_con);
        // e3 = a3_KL

        double mG_00 = inner_prod(e1, rg1);
        double mG_01 = inner_prod(e1, rg2);
        double mG_02 = inner_prod(e1, rg3);
        double mG_11 = inner_prod(e2, rg2);
        double mG_12 = inner_prod(e2, rg3);
        double mG_22 = inner_prod(ra3_KL, rg3);
        // transformation matrix T_cov_to_car from covariant to local Cartesian basis
        rActualTransCovToCar(0, 0) = pow(mG_00, 2);
        rActualTransCovToCar(0, 1) = pow(mG_01, 2);
        rActualTransCovToCar(0, 2) = 2 * mG_00 * mG_01;
        rActualTransCovToCar(0, 3) = 2 * mG_00 * mG_01;
        rActualTransCovToCar(0, 4) = 2 * mG_01 * mG_02;
        rActualTransCovToCar(0, 5) = 2 * mG_00 * mG_02;
        rActualTransCovToCar(1, 1) = pow(mG_11, 2);
        rActualTransCovToCar(1, 2) = pow(mG_12, 2);
        rActualTransCovToCar(1, 4) = 2 * mG_11 * mG_12;
        rActualTransCovToCar(2, 2) = pow(mG_22, 2);
        rActualTransCovToCar(3, 1) = mG_01 * mG_11;
        rActualTransCovToCar(3, 2) = mG_02 * mG_12;
        rActualTransCovToCar(3, 3) = mG_00 * mG_11;
        rActualTransCovToCar(3, 4) = mG_01 * mG_12 + mG_02 * mG_11;
        rActualTransCovToCar(3, 5) = mG_00 * mG_12;
        rActualTransCovToCar(4, 2) = mG_12 * mG_22;
        rActualTransCovToCar(4, 4) = mG_11;
        rActualTransCovToCar(5, 2) = mG_02 * mG_22;
        rActualTransCovToCar(5, 4) = mG_01 * mG_22;
        rActualTransCovToCar(5, 5) = mG_00 * mG_22;
    }

    void IgaShell7pElement::CalculateShearDifferenceVector(
        array_1d<double, 3>& rw,
        array_1d<double, 3>& rDw_D1,
        array_1d<double, 3>& rDw_D2,
        array_1d<double, 3>& rw_i,
        Matrix& rDw_i_Dalpha,
        const MetricVariables& rActualMetric)
    {
        KRATOS_TRY; 
        
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const unsigned int number_of_control_points = GetGeometry().size();
        double w_1, w_2, w_3;

        for (unsigned int i = 0; i < number_of_control_points; ++i) 
        {
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            w_1 = GetGeometry()[i].GetDof(ROTATION_X).GetSolutionStepValue();
            w_2 = GetGeometry()[i].GetDof(ROTATION_Y).GetSolutionStepValue();
            w_3 = GetGeometry()[i].GetDof(ROTATION_Z).GetSolutionStepValue();

            rDw_i_Dalpha(0, 0) += DN_De(i, 0) * w_1;
            rDw_i_Dalpha(0, 1) += DN_De(i, 1) * w_1;
            rDw_i_Dalpha(1, 0) += DN_De(i, 0) * w_2;
            rDw_i_Dalpha(1, 1) += DN_De(i, 1) * w_2;
            rDw_i_Dalpha(2, 0) += DN_De(i, 0) * w_3; 
            rDw_i_Dalpha(2, 1) += DN_De(i, 1) * w_3; 

            rw_i(0) += N[i] * w_1;
            rw_i(1) += N[i] * w_2;
            rw_i(2) += N[i] * w_3;
        }

        // derivatives of the shear difference vector
        for (unsigned int i = 0; i < 3; i++){
            rDw_D1[i] = rDw_i_Dalpha(i, 0);
            rDw_D2[i] = rDw_i_Dalpha(i, 1);
            // shear difference vector rw is defined w.r.t. a Cartesian coordinate system
            rw[i] = rw_i[i];
        }

        KRATOS_CATCH("")
    }

    void IgaShell7pElement::CalculateInitialBaseVectorsLinearised(
        array_1d<double, 3>& rG1,
        array_1d<double, 3>& rG2,
        array_1d<double, 3>& rG1_con,
        array_1d<double, 3>& rG2_con)
    {
        double thickness = GetProperties().GetValue(THICKNESS);
        
        array_1d<double, 3> DA3_D1(3, 0.0), DA3_D2(3, 0.0), DA1_D1xA2(3, 0.0), A1xDA2_D1(3, 0.0), DA1_D2xA2(3, 0.0), A1xDA2_D2(3, 0.0), 
            G1xG2(3, 0.0), G_ab(3, 0.0), G_ab_con(3, 0.0);

        MathUtils<double>::CrossProduct(DA1_D1xA2, mInitialMetric.Da1_D1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D1, mInitialMetric.a1, mInitialMetric.Da1_D2); // DA1_D2 = DA2_D1
        MathUtils<double>::CrossProduct(DA1_D2xA2, mInitialMetric.Da1_D2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D2, mInitialMetric.a1, mInitialMetric.Da2_D2);
        DA3_D1 = ((DA1_D1xA2 + A1xDA2_D1) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(DA1_D1xA2 + A1xDA2_D1)) 
            / (mInitialMetric.dA * mInitialMetric.dA);
        DA3_D2 = ((DA1_D2xA2 + A1xDA2_D2) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(DA1_D2xA2 + A1xDA2_D2))
            / (mInitialMetric.dA * mInitialMetric.dA);

        // covariant base vectors of the shell body in the reference configuration
        rG1 = mInitialMetric.a1 + thickness / 2.0 * mZeta * DA3_D1;
        rG2 = mInitialMetric.a2 + thickness / 2.0 * mZeta * DA3_D2;
        // G3 = A3

                // covariant metric
        G_ab[0] = pow(rG1[0], 2) + pow(rG1[1], 2) + pow(rG1[2], 2);
        G_ab[1] = pow(rG2[0], 2) + pow(rG2[1], 2) + pow(rG2[2], 2);
        G_ab[2] = rG1[0] * rG2[0] + rG1[1] * rG2[1] + rG1[2] * rG2[2];

        // contravariant metric
        double inv_detG_ab = 1.0 / (G_ab[0] * G_ab[1] - G_ab[2] * G_ab[2]);
        G_ab_con[0] = inv_detG_ab * G_ab[1];
        G_ab_con[2] = -inv_detG_ab * G_ab[2];
        G_ab_con[1] = inv_detG_ab * G_ab[0];

        // contravariant base vectors a_con
        rG1_con = rG1 * G_ab_con[0] + rG2 * G_ab_con[2];
        rG2_con = rG1 * G_ab_con[2] + rG2 * G_ab_con[1];
        // G3_con = A3
    }

    void IgaShell7pElement::CalculateActualBaseVectorsLinearised(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const double& rw_bar,
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
        double Dw_bar_D1 = 0.0;
        double Dw_bar_D2 = 0.0;

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

        // derivatives of w_bar w.r.t. curvilinear coordinates
        for(IndexType i = 0; i < number_of_control_points; i++){
            Dw_bar_D1 += DN_De(i, 0) * GetGeometry()[i].GetDof(W_BAR).GetSolutionStepValue();
            Dw_bar_D2 += DN_De(i, 1) * GetGeometry()[i].GetDof(W_BAR).GetSolutionStepValue();
        }
        /* Kovariante Basisvektoren */
        rg1 = rActualMetric.a1 + thickness / 2.0 * mZeta * (Da3_KL_D1 + rDw_D1) + 
            pow(thickness / 2.0 * mZeta, 2) * (Dw_bar_D1 * (rActualMetric.a3_KL + rw) + rw_bar * (Da3_KL_D1 + rDw_D1));
        rg2 = rActualMetric.a2 + thickness / 2.0 * mZeta * (Da3_KL_D2 + rDw_D2) + 
            pow(thickness / 2.0 * mZeta, 2) * (Dw_bar_D2 * (rActualMetric.a3_KL + rw) + rw_bar * (Da3_KL_D2 + rDw_D2));
        rg3 = rActualMetric.a3_KL + rw + thickness * mZeta * rw_bar * (rActualMetric.a3_KL + rw);

        KRATOS_CATCH("")
    }

    void IgaShell7pElement::CalculateDeformationGradient(
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

    void IgaShell7pElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const double& rw_bar,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        KRATOS_TRY

        Vector strain_vector = ZeroVector(6);
        Vector strain_vector_7p = ZeroVector(6);
        
        // Strain computation in curvilinear space
        CalculateStrain(strain_vector, rActualMetric.gab, rActualMetric.curvature);
        CalculateStrain7p(strain_vector_7p, rw, rDw_D1, rDw_D2, rw_bar, rActualMetric);
        rThisConstitutiveVariables.E = strain_vector + strain_vector_7p;
        rThisConstitutiveVariables.E = prod(mInitialTransConToCar, rThisConstitutiveVariables.E);
        if(Id()==1){
            KRATOS_WATCH(strain_vector)
            KRATOS_WATCH(strain_vector_7p)
            KRATOS_WATCH(rThisConstitutiveVariables.E)
        }
        //Constitutive Matrix D
        rValues.SetStrainVector(rThisConstitutiveVariables.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        //Local Cartesian Stresses
        rThisConstitutiveVariables.S = prod(
            trans(rThisConstitutiveVariables.D), rThisConstitutiveVariables.E);

        KRATOS_CATCH("")
    }

    void IgaShell7pElement::CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab,
        const Vector& rCurvature)
    {
        KRATOS_TRY
        
        double thickness = GetProperties().GetValue(THICKNESS);

        rStrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.gab[0]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[0] - rCurvature[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.gab[1]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[1] - rCurvature[1]);
        rStrainVector[3] = 0.5 * (rgab[2] - mInitialMetric.gab[2]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[2] - rCurvature[2]);
        // the other entries are (remain) zero (KL)
        
        KRATOS_CATCH("")
    }

    void IgaShell7pElement::CalculateStrain7p(
        Vector& rStrainVector7p,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const double& rw_bar,
        const MetricVariables rActualMetric)
    {
        double thickness = GetProperties().GetValue(THICKNESS);
                
        rStrainVector7p[0] = mZeta * thickness/2.0 * inner_prod(rDw_D1, rActualMetric.a1);
        rStrainVector7p[1] = mZeta * thickness/2.0 * inner_prod(rDw_D2, rActualMetric.a2);
        rStrainVector7p[2] = 0.5 *  inner_prod(rw, rw) + inner_prod(rActualMetric.a3_KL, rw) + mZeta * thickness * rw_bar * 
            (inner_prod(rActualMetric.a3_KL, rActualMetric.a3_KL) + 2.0 * inner_prod(rActualMetric.a3_KL, rw) + inner_prod(rw, rw));
        rStrainVector7p[3] = mZeta * thickness/2.0 * 0.5 * (inner_prod(rDw_D1, rActualMetric.a2) + inner_prod(rDw_D2, rActualMetric.a1));
        rStrainVector7p[4] = 0.5 * inner_prod(rw, rActualMetric.a2) + mZeta * thickness / 2.0 * (inner_prod(rDw_D2, rActualMetric.a3_KL) + 
            inner_prod(rw, rActualMetric.Da3_KL_D2) + inner_prod(rDw_D2, rw) + 2.0 * rw_bar * inner_prod(rActualMetric.a2, rw));
        rStrainVector7p[5] = 0.5 * inner_prod(rw, rActualMetric.a1) + mZeta * thickness / 2.0 * (inner_prod(rDw_D1, rActualMetric.a3_KL) +
            inner_prod(rw, rActualMetric.Da3_KL_D1) + inner_prod(rDw_D1, rw) + 2.0 * rw_bar * inner_prod(rActualMetric.a1, rw));
    }

    void IgaShell7pElement::CalculateB(
        Matrix& rB,
        const MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const double thickness = GetProperties().GetValue(THICKNESS);

        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size_KL = number_of_control_points * 3;
        const unsigned int mat_size = number_of_control_points * 7;

        // membrane part
        for (unsigned int r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 7;
            int dirr = r % 7;

            Vector dE_curvilinear = ZeroVector(3);
            // "if" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            if (dirr == 0 || dirr == 1 || dirr == 2)
            {
                // strain corresponding to E11, E22, E12
                dE_curvilinear[0] = DN_De(kr, 0)*rMetric.a1(dirr);
                dE_curvilinear[1] = DN_De(kr, 1)*rMetric.a2(dirr);
                dE_curvilinear[2] = 0.5*(DN_De(kr, 0)*rMetric.a2(dirr) + rMetric.a1(dirr)*DN_De(kr, 1));
            }
            rB(0, r) += mInitialTransConToCar(0, 0) * dE_curvilinear[0] + mInitialTransConToCar(0, 1)*dE_curvilinear[1] 
                + mInitialTransConToCar(0, 3)*dE_curvilinear[2] ;
            rB(1, r) += mInitialTransConToCar(1, 0) * dE_curvilinear[0] + mInitialTransConToCar(1, 1) * dE_curvilinear[1] 
                + mInitialTransConToCar(1, 3) * dE_curvilinear[2];
            rB(3, r) += mInitialTransConToCar(3, 0)*dE_curvilinear[0] + mInitialTransConToCar(3, 1)*dE_curvilinear[1] 
                + mInitialTransConToCar(3, 3)*dE_curvilinear[2];

            // all other entries of rB are (remain) zero
            
        }

        // curvature part
        Matrix dg3 = ZeroMatrix(3, 3);
        Matrix dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size_KL);

        double invdA = 1 / rMetric.dA;
        double inddA3 = 1 / std::pow(rMetric.dA, 3);

        for (unsigned int i = 0; i < number_of_control_points; i++)
        {
            unsigned int index_KL = 3 * i;
            unsigned int index = 7 * i;
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
            // "index" guarantees that there are zero entries corresponding to the new parameters w_1, w_2, w_3 and w_bar
            for (unsigned int j = 0; j < 3; j++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    b(k, index_KL + j) = - mZeta * thickness / 2.0 * (DDN_DDe(i, k) * rMetric.a3_KL[j] + rMetric.H(0, k) * dn(j, 0) 
                    + rMetric.H(1, k) * dn(j, 1) + rMetric.H(2, k) * dn(j, 2));
                }
                rB(0, index + j) += mInitialTransConToCar(0, 0) * b(0, index_KL + j);
                rB(1, index + j) += mInitialTransConToCar(1, 0) * b(0, index_KL + j) + mInitialTransConToCar(1, 1) * b(1, index_KL + j)
                    + mInitialTransConToCar(1, 3) * b(2, index_KL + j);
                rB(3, index + j) += mInitialTransConToCar(3, 0) * b(0, index_KL + j) + mInitialTransConToCar(3, 3) * b(2, index_KL + j);
            }
        }
    }

    void IgaShell7pElement::CalculateSecondVariations(
        SecondVariations& rSecondVariations,
        const MetricVariables& rMetric)
    {
        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
            const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
            const double thickness = GetProperties().GetValue(THICKNESS);

            const int number_of_control_points = GetGeometry().size();
            const int mat_size_KL = number_of_control_points * 3;
            const int mat_size = number_of_control_points * 7;
           
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

                        // calculated with simplified T_con_to_car (ML)
                        second_variations_KL.B11(r, s) += mInitialTransConToCar(0, 0) * ddE_cu[0];
                        second_variations_KL.B22(r, s) += mInitialTransConToCar(1, 0) * ddE_cu[0] + mInitialTransConToCar(1, 1) * ddE_cu[1]
                            + mInitialTransConToCar(1, 2) * ddE_cu[2];
                        second_variations_KL.B12(r, s) += mInitialTransConToCar(2, 0) * ddE_cu[0] + mInitialTransConToCar(2, 2) * ddE_cu[2];
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
                    ddK_cu[0] = - mZeta * thickness / 2.0 * (DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2]);
                    ddK_cu[1] = - mZeta * thickness / 2.0 * (DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2]);
                    ddK_cu[2] = - mZeta * thickness / 2.0 * (DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2]);

                    // calculated with simplified T_con_to_car (ML)
                    second_variations_KL.B11(r, s) += mInitialTransConToCar(0, 0) * ddK_cu[0];
                    second_variations_KL.B22(r, s) += mInitialTransConToCar(1, 0) * ddK_cu[0] + mInitialTransConToCar(1, 1) * ddK_cu[1] 
                        + mInitialTransConToCar(1, 2) * ddK_cu[2];
                    second_variations_KL.B12(r, s) += mInitialTransConToCar(2, 0) * ddK_cu[0] + mInitialTransConToCar(2, 2) * ddK_cu[2];
                    if (r != s){
                        second_variations_KL.B11(s, r) = second_variations_KL.B11(r, s);
                        second_variations_KL.B22(s, r) = second_variations_KL.B22(r, s);
                        second_variations_KL.B12(s, r) = second_variations_KL.B12(r, s);
                    }
                }
            }

            // transfer KL-second-variations to RM-second-variations
            for (unsigned int r = 0; r < mat_size; r++) {
                unsigned int kr = r / 7;
                unsigned int dirr = r % 7;
                unsigned int r_KL = kr * 3 + dirr;
                if (dirr != 3 && dirr != 4 && dirr != 5 && dirr != 6){
                    for (unsigned int s = 0; s < mat_size; s++){
                        unsigned int ks = s / 7;
                        unsigned int dirs = s % 7;
                        unsigned int s_KL = ks * 3 + dirs;
                        if (dirr != 3 && dirr != 4 && dirr != 5 && dirr != 6){
                            rSecondVariations.B11(r, s) += second_variations_KL.B11(r_KL, s_KL);
                            rSecondVariations.B22(r, s) += second_variations_KL.B22(r_KL, s_KL);
                            rSecondVariations.B12(r, s) += second_variations_KL.B12(r_KL, s_KL);                
                        }
                    }
                }
            }
        }
    }

    void IgaShell7pElement::CalculateVariations7p(        
        Matrix& rB,
        SecondVariations& rSecondVariations,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rw_i,
        const Matrix& rDw_i_Dalpha,
        const double& rw_bar,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag)
    {
        // KRATOS_WATCH("begin:CalculateVariations7p")
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 7;

        const Vector& funct = N;
        const Matrix& deriv = DN_De;
        const Matrix& s_deriv = DDN_DDe;
        const unsigned int num_node = number_of_control_points;
        const unsigned int dof_per_node = 7;

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
        deriv_radicant_normA1crossA2_dT1 = 0.0;
        deriv_radicant_normA1crossA2_dT2 = 0.0;
        dnormA1crossA2_dT1 = 0.0;
        dnormA1crossA2_dT2 = 0.0;
        dA3_dT1 = ZeroVector(3);
        dA3_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der Referenzkonfiguration */
        A1 = mInitialMetric.a1;
        A2 = mInitialMetric.a2;

        /* Dritter kovarianter Basisvektor senkrecht zur Mittelflaeche der Referenzkonfiguration */
        A1xA2 = mInitialMetric.a3_KL_tilde;
        h = mInitialMetric.dA;
        A3 = mInitialMetric.a3_KL; // FML: bei Oesterle mit t/2 multipliziert

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
            dA3_dT1[k] = (d_A1crossA2_dT1[k] * h - A1xA2[k] * dnormA1crossA2_dT1) / (h * h);    // FML: bei Oesterle mit t/2 multipliziert
            dA3_dT2[k] = (d_A1crossA2_dT2[k] * h - A1xA2[k] * dnormA1crossA2_dT2) / (h * h);    // FML: bei Oesterle mit t/2 multipliziert
        }

        /* AKTUELLE KONFIGURATION */
        int m;
        double hact, deriv_sqrt_norma1crossa2, deriv_radicant_norma1crossa2_dT1,
            deriv_radicant_norma1crossa2_dT2, dnorma1crossa2_dT1, dnorma1crossa2_dT2 , wbar;
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
        array_1d<double, 3> da3_KL_dT1;
        array_1d<double, 3> da3_KL_dT2;
        array_1d<double, 3> w;
        array_1d<double, 3> da1norm_dT1;
        array_1d<double, 3> da1norm_dT2;
        array_1d<double, 3> da2norm_dT1;
        array_1d<double, 3> da2norm_dT2;

        array_1d<double, 3> e1, e2, e3;     // GLobal kartesische basisvektoren

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
        deriv_radicant_norma1crossa2_dT1 = 0.0;
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        dnorma1crossa2_dT1 = 0.0;
        dnorma1crossa2_dT2 = 0.0;
        wbar=0.0;
        da3_KL_dT1 = ZeroVector(3);
        da3_KL_dT2 = ZeroVector(3);
        w = ZeroVector(3);

        /* Einfuehrung eines global kartesischen Koordinatensystems*/
        e1 = ZeroVector(3);
        e2 = ZeroVector(3);
        e3 = ZeroVector(3);

        e1[0] = 1;
        e1[1] = 0;
        e1[2] = 0;
        e2[0] = 0;
        e2[1] = 1;
        e2[2] = 0;
        e3[0] = 0;
        e3[1] = 0;
        e3[2] = 1;

	    ///////////////////////////////////
        /* Kovariante Basisvektoren der Mittelflaeche der aktuellen Konfiguration */
        a1 = rActualMetric.a1;
        a2 = rActualMetric.a2;

        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        a1xa2 = rActualMetric.a3_KL_tilde;
        hact = sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]);

        for (m = 0; m < 3; m++)
            a3_KL[m] = a1xa2[m] / hact;     // FML: bei Oesterle multipliziert mit t/2

        double w1,w2,w3;
        w1=0.0;
        w2=0.0;
        w3=0.0;
        /* aktuelle Komponenten vom Schubdifferenzvektor w */
        for (k = 0; k < num_node; k++)
        {
            w1   += funct[k] * GetGeometry()[k].GetDof(ROTATION_X).GetSolutionStepValue();
            w2   += funct[k] * GetGeometry()[k].GetDof(ROTATION_Y).GetSolutionStepValue();
            w3   += funct[k] * GetGeometry()[k].GetDof(ROTATION_Z).GetSolutionStepValue();
            wbar += funct[k] * GetGeometry()[k].GetDof(W_BAR).GetSolutionStepValue();
        }  
        w = w1*e1 + w2*e2 + w3*e3;

        /* a3 nach hier. 7P-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
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

        /* Ableitung von a3 mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_KL_dT1[m] = (d_a1crossa2_dT1[m] * hact - a1xa2[m] * dnorma1crossa2_dT1) / (hact * hact);    // FML: bei Oesterle mit t/2 multipliziert
            da3_KL_dT2[m] = (d_a1crossa2_dT2[m] * hact - a1xa2[m] * dnorma1crossa2_dT2) / (hact * hact);    // FML: bei Oesterle mit t/2 multipliziert
        }

        /* Partielle Ableitungen von Schubdifferenzvektor w nach alpha/beta */
        double dw1dT1,dw2dT1,dw3dT1,dw1dT2,dw2dT2,dw3dT2;
        dw1dT1=0.0;
        dw2dT1=0.0;
        dw3dT1=0.0;
        dw1dT2=0.0;
        dw2dT2=0.0;
        dw3dT2=0.0;

        for (k = 0; k < num_node; k++)
        {
            dw1dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_X).GetSolutionStepValue();
            dw2dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_Y).GetSolutionStepValue();
            dw3dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_Z).GetSolutionStepValue();
            dw1dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_X).GetSolutionStepValue();
            dw2dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_Y).GetSolutionStepValue();
            dw3dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_Z).GetSolutionStepValue();
        }

        array_1d<double, 3>  dw_dT1 = dw1dT1*e1 + dw2dT1*e2 + dw3dT1*e3;
        array_1d<double, 3>  dw_dT2 = dw1dT2*e1 + dw2dT2*e2 + dw3dT2*e3;


        /* Nichtlinearer B-Operator */

        array_1d<double, 3> da1xa2;
        array_1d<double, 3> da1;
        array_1d<double, 3> da2;
        array_1d<double, 3> da3_KL_dr;

        da1xa2 = ZeroVector(3);
        da1 = ZeroVector(3);
        da2 = ZeroVector(3);
        da3_KL_dr = ZeroVector(3);

        double da1xa2_a1xa2;
        double da1xa2_a1xa2_hact;

        int count1 = 0;
        int count2 = 0;

        
            //////////////////// Additional terms for the expression "a3_KL,alpha,r" ///////////////////////
        double dnorma1crossa2_dr, ddnorma1crossa2_dT1dr, ddnorma1crossa2_dT2dr;
        array_1d<double, 3> dda1crossa2_dT1dr;
        array_1d<double, 3> dda1crossa2_dT2dr;
        array_1d<double, 3> dda3_KL_dT1dr;
        array_1d<double, 3> dda3_KL_dT2dr;
	  

        dnorma1crossa2_dr = 0.0;
        ddnorma1crossa2_dT1dr = 0.0;
        ddnorma1crossa2_dT2dr = 0.0;
        dda1crossa2_dT1dr = ZeroVector(3);
        dda1crossa2_dT2dr = ZeroVector(3);
        dda3_KL_dT1dr = ZeroVector(3);
        dda3_KL_dT2dr = ZeroVector(3);

        ////////////// Definition of a1,1,r ; a1,2,r=a,2,1,r and a2,2,r. They are needed to calculate "a3_KL,alpha,r" /////////////////////

        array_1d<double, 3> dda1dT1dr;
        array_1d<double, 3> dda1dT2dr;
        array_1d<double, 3> dda2dT1dr;
        array_1d<double, 3> dda2dT2dr;

        dda1dT1dr = ZeroVector(3);
        dda1dT2dr = ZeroVector(3);
        dda2dT1dr = ZeroVector(3);
        dda2dT2dr = ZeroVector(3);
	  
        ////////////////////////////////////////////////////////////////////////////////////////////////

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
            da3_KL_dr[0] = da1xa2[0] / hact - a1xa2[0] * da1xa2_a1xa2_hact;
            da3_KL_dr[1] = da1xa2[1] / hact - a1xa2[1] * da1xa2_a1xa2_hact;
            da3_KL_dr[2] = da1xa2[2] / hact - a1xa2[2] * da1xa2_a1xa2_hact;


            if (i<3)
            {
                /* a3_KL anstatt a3! ansonsten wie bei 3P-Schale */
                // bop[0][j*7 + i] = bop[0][j*7 + i] + deriv(count2, 0)*a1[i] - mZeta * thickness / 2.0*(s_deriv(count2, 0)*a3_KL[i] + thickness/2.0*(da1_dT1[0]*da3_KL_dr[0] + da1_dT1[1]*da3_KL_dr[1] + da1_dT1[2]*da3_KL_dr[2])) ;
                // bop[1][j*7 + i] = bop[1][j*7 + i] + 2.0*(0.5*(deriv(count2, 1)*a1[i] + deriv(count2, 0)*a2[i]))
                //                             - 2.0*mZeta * thickness / 2.0*(s_deriv(count2, 2)*a3_KL[i] + thickness/2.0*(da1_dT2[0]*da3_KL_dr[0] + da1_dT2[1]*da3_KL_dr[1] + da1_dT2[2]*da3_KL_dr[2]));
                // bop[2][j*7 + i] = bop[2][j*7 + i] + deriv(count2, 1)*a2[i] - mZeta * thickness / 2.0*(s_deriv(count2, 1)*a3_KL[i] + thickness/2.0*(da2_dT2[0]*da3_KL_dr[0] + da2_dT2[1]*da3_KL_dr[1] + da2_dT2[2]*da3_KL_dr[2]))  ;
            }


            }
            count2 ++;
        }

        /* Zusaetzliche Anteile aus Schubvektor "w" und 7. Parameter "wbar" */
        double dw1_dr;
        double dw2_dr;
        double dw3_dr;
        double dw1dT1_dr;
        double dw2dT1_dr;
        double dw3dT1_dr;
        double dw1dT2_dr;
        double dw2dT2_dr;
        double dw3dT2_dr;
        double dwbar_dr;
        count1 = 0;
        count2 = 0;

        array_1d<double, 3> a1_dr;
        array_1d<double, 3> a2_dr;
        array_1d<double, 3> dw_dr;
        array_1d<double, 3> dw_dT1_dr;
        array_1d<double, 3> dw_dT2_dr;

        a1_dr = ZeroVector(3);
        a2_dr = ZeroVector(3);
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

                dw1_dr = 0.0;
                dw2_dr = 0.0;
                dw3_dr = 0.0;
                dw1dT1_dr = 0.0;
                dw2dT1_dr = 0.0;
                dw3dT1_dr = 0.0;
                dw1dT2_dr = 0.0;
                dw2dT2_dr = 0.0;
                dw3dT2_dr = 0.0;
                
                dwbar_dr = 0.0;
            }
            else if (i == 1)
            {
                a1_dr[0] = 0;
                a1_dr[1] = deriv(count2, 0);
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = deriv(count2, 1);
                a2_dr[2] = 0;

                dw1_dr = 0.0;
                dw2_dr = 0.0;
                dw3_dr = 0.0;
                dw1dT1_dr = 0.0;
                dw2dT1_dr = 0.0;
                dw3dT1_dr = 0.0;
                dw1dT2_dr = 0.0;
                dw2dT2_dr = 0.0;
                dw3dT2_dr = 0.0;
                
                dwbar_dr = 0.0;
            }
            else if (i == 2)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = deriv(count2, 0);
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = deriv(count2, 1);

                dw1_dr = 0.0;
                dw2_dr = 0.0;
                dw3_dr = 0.0;
                dw1dT1_dr = 0.0;
                dw2dT1_dr = 0.0;
                dw3dT1_dr = 0.0;
                dw1dT2_dr = 0.0;
                dw2dT2_dr = 0.0;
                dw3dT2_dr = 0.0;
                
                dwbar_dr = 0.0;
            }
            else if (i == 3)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;

                dw1_dr = funct[count2];
                dw2_dr = 0.0;
                dw3_dr = 0.0;
                dw1dT1_dr = deriv(count2, 0);
                dw2dT1_dr = 0.0;
                dw3dT1_dr = 0.0;
                dw1dT2_dr = deriv(count2, 1);
                dw2dT2_dr = 0.0;
                dw3dT2_dr = 0.0;
                
                dwbar_dr = 0.0;
            }
            else if (i == 4)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;

                dw1_dr = 0.0;
                dw2_dr = funct[count2];
                dw3_dr = 0.0;
                dw1dT1_dr = 0.0;
                dw2dT1_dr = deriv(count2, 0);
                dw3dT1_dr = 0.0;
                dw1dT2_dr = 0.0;
                dw2dT2_dr = deriv(count2, 1);
                dw3dT2_dr = 0.0;
                
                dwbar_dr = 0.0;
            }
            else if (i == 5)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;

                dw1_dr = 0.0;
                dw2_dr = 0.0;
                dw3_dr = funct[count2];
                dw1dT1_dr = 0.0;
                dw2dT1_dr = 0.0;
                dw3dT1_dr = deriv(count2, 0);
                dw1dT2_dr = 0.0;
                dw2dT2_dr = 0.0;
                dw3dT2_dr = deriv(count2, 1);
                
                dwbar_dr = 0.0;
            }
            else if (i == 6)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;

                dw1_dr = 0.0;
                dw2_dr = 0.0;
                dw3_dr = 0.0;
                dw1dT1_dr = 0.0;
                dw2dT1_dr = 0.0;
                dw3dT1_dr = 0.0;
                dw1dT2_dr = 0.0;
                dw2dT2_dr = 0.0;
                dw3dT2_dr = 0.0;
                        
                dwbar_dr = funct[count2];
            }

            da1xa2 = CrossProduct(a1_dr, a2) + CrossProduct(a1, a2_dr);
            da1xa2_a1xa2 = (da1xa2[0] * a1xa2[0] + da1xa2[1] * a1xa2[1]
                + da1xa2[2] * a1xa2[2]);
            da1xa2_a1xa2_hact = da1xa2_a1xa2 / (hact * hact * hact);
            da3_KL_dr[0] = da1xa2[0] / hact - a1xa2[0] * da1xa2_a1xa2_hact;
            da3_KL_dr[1] = da1xa2[1] / hact - a1xa2[1] * da1xa2_a1xa2_hact;
            da3_KL_dr[2] = da1xa2[2] / hact - a1xa2[2] * da1xa2_a1xa2_hact;

            /* Partielle Ableitung/Variation nach FHG (_dr entspricht ,r) */
            dw_dr = dw1_dr*e1 + dw2_dr*e2 + dw3_dr*e3;
            dw_dT1_dr = dw1dT1_dr*e1 +  dw2dT1_dr*e2 + dw3dT1_dr*e3;
            dw_dT2_dr = dw1dT2_dr*e1 +  dw2dT2_dr*e2 + dw3dT2_dr*e3;


            /* These terms are later used in "dda1crossa2_dT1dr" and "dda1crossa2_dT2dr" equations below.  */
            if (i == 0)
            {
                dda1dT1dr[0] = s_deriv(count2, 0);
                dda1dT1dr[1] = 0;
                dda1dT1dr[2] = 0;
                dda2dT2dr[0] = s_deriv(count2, 1);
                dda2dT2dr[1] = 0;
                dda2dT2dr[2] = 0;
                dda1dT2dr[0] = s_deriv(count2, 2);
                dda1dT2dr[1] = 0;
                dda1dT2dr[2] = 0;
            }
            else if (i == 1)
            {
                dda1dT1dr[0] = 0;
                dda1dT1dr[1] = s_deriv(count2, 0);
                dda1dT1dr[2] = 0;
                dda2dT2dr[0] = 0;
                dda2dT2dr[1] = s_deriv(count2, 1);
                dda2dT2dr[2] = 0;
                dda1dT2dr[0] = 0;
                dda1dT2dr[1] = s_deriv(count2, 2);
                dda1dT2dr[2] = 0;
            }
            else if (i == 2)
            {
                dda1dT1dr[0] = 0;
                dda1dT1dr[1] = 0;
                dda1dT1dr[2] = s_deriv(count2, 0);
                dda2dT2dr[0] = 0;
                dda2dT2dr[1] = 0;
                dda2dT2dr[2] = s_deriv(count2, 1);
                dda1dT2dr[0] = 0;
                dda1dT2dr[1] = 0;
                dda1dT2dr[2] = s_deriv(count2, 2);
            }
            else if (i == 3 || i == 4 || i == 5 || i == 6 )
            {
                dda1dT1dr[0] = 0;
                dda1dT1dr[1] = 0;
                dda1dT1dr[2] = 0;
                dda2dT2dr[0] = 0;
                dda2dT2dr[1] = 0;
                dda2dT2dr[2] = 0;
                dda1dT2dr[0] = 0;
                dda1dT2dr[1] = 0;
                dda1dT2dr[2] = 0;
            }

            dda2dT1dr = dda1dT2dr;

            /////////////////////////////////////////////////////////////////////////////////////////////////

            //////////////////// Additional terms are now computed for the expression "a3_KL,alpha,r" ///////////////////////
            /* (a1xa2),1,r & (a1xa2),2,r */

            dda1crossa2_dT1dr = CrossProduct(dda1dT1dr, a2) + CrossProduct(da1_dT1, da2) + CrossProduct(da1, da2_dT1) + CrossProduct(a1, dda2dT1dr);
            dda1crossa2_dT2dr = CrossProduct(dda1dT2dr, a2) + CrossProduct(da1_dT2, da2) + CrossProduct(da1, da2_dT2) + CrossProduct(a1, dda2dT2dr);

            /* (||a1xa2||),r */

            dnorma1crossa2_dr = inner_prod(a1xa2, da1xa2) / hact;

            /* (||a1xa2||),1,r & (||a1xa2||),2,r */

            ddnorma1crossa2_dT1dr = - inner_prod(a1xa2, da1xa2)*inner_prod(a1xa2, d_a1crossa2_dT1) / (hact*hact*hact)
                                                + inner_prod(da1xa2, d_a1crossa2_dT1) / hact + inner_prod(a1xa2, dda1crossa2_dT1dr) / hact;

            ddnorma1crossa2_dT2dr = - inner_prod(a1xa2, da1xa2)*inner_prod(a1xa2, d_a1crossa2_dT2) / (hact*hact*hact)
                                                + inner_prod(da1xa2, d_a1crossa2_dT2) / hact + inner_prod(a1xa2, dda1crossa2_dT2dr) / hact;

            /* now write the equations a3_KL,1,r & a3_KL,2,r */
            for (k = 0; k < 3; k++)
                {

                dda3_KL_dT1dr[k] = ((dda1crossa2_dT1dr[k]*hact + d_a1crossa2_dT1[k]*dnorma1crossa2_dr - da1xa2[k]*dnorma1crossa2_dT1 - a1xa2[k]*ddnorma1crossa2_dT1dr)*hact*hact
                                    - 2.0*(d_a1crossa2_dT1[k]*hact - a1xa2[k]*dnorma1crossa2_dT1)*hact*dnorma1crossa2_dr )  * thickness / (2.0*hact*hact*hact*hact) ;

                dda3_KL_dT2dr[k] = ((dda1crossa2_dT2dr[k]*hact + d_a1crossa2_dT2[k]*dnorma1crossa2_dr - da1xa2[k]*dnorma1crossa2_dT2 - a1xa2[k]*ddnorma1crossa2_dT2dr)*hact*hact
                                    - 2.0*(d_a1crossa2_dT2[k]*hact - a1xa2[k]*dnorma1crossa2_dT2)*hact*dnorma1crossa2_dr ) * thickness / (2.0*hact*hact*hact*hact)  ;

                }

            /////////////////////////////////////////////////////////////////////////////////////////////////
            array_1d<double, 6> dE_cur_help = ZeroVector(6);
            if (i<3)
                    {
                    dE_cur_help[5] = 0.5 * mZeta * thickness / 2.0*inner_prod(dda3_KL_dT1dr, w);
                    dE_cur_help[4] = 0.5 * mZeta * thickness / 2.0*inner_prod(dda3_KL_dT2dr, w);
                    }

            //additional terms coming from fully fully exact director are also added!
            dE_cur_help[0] = mZeta * thickness / 2.0*(inner_prod(a1_dr, dw_dT1) + inner_prod(a1, dw_dT1_dr));
            dE_cur_help[3] = 0.5 * mZeta * thickness / 2.0*(inner_prod(a1_dr, dw_dT2) + inner_prod(a1, dw_dT2_dr) + inner_prod(a2_dr, dw_dT1) + inner_prod(a2, dw_dT1_dr));
            dE_cur_help[1] = mZeta * thickness / 2.0*(inner_prod(a2_dr, dw_dT2) + inner_prod(a2, dw_dT2_dr));
            dE_cur_help[5] += 0.5 * (inner_prod(dw_dr, a1) + inner_prod(w, a1_dr) + mZeta * thickness / 2.0*(/*dda3_KL_dT1dr*w  + */ inner_prod(da3_KL_dT1, dw_dr) + inner_prod(da3_KL_dr, dw_dT1)+ inner_prod(a3_KL, dw_dT1_dr) + inner_prod(dw_dT1_dr, w) + inner_prod(dw_dT1, dw_dr) + 2.0*dwbar_dr*inner_prod(a1, w) + 2.0*wbar*inner_prod(a1_dr, w) + 2.0*wbar*inner_prod(a1, dw_dr))); // FML: bei Oesterle zusätzliche t/2
            dE_cur_help[4] += 0.5 * (inner_prod(dw_dr, a2) + inner_prod(w, a2_dr) + mZeta * thickness / 2.0*(/*dda3_KL_dT2dr*w  + */ inner_prod(da3_KL_dT2, dw_dr) + inner_prod(da3_KL_dr, dw_dT2)+ inner_prod(a3_KL, dw_dT2_dr) + inner_prod(dw_dT2_dr, w) + inner_prod(dw_dT2, dw_dr) + 2.0*dwbar_dr*inner_prod(a2, w) + 2.0*wbar*inner_prod(a2_dr, w) + 2.0*wbar*inner_prod(a2, dw_dr))); // FML: bei Oesterle zusätzliche t/2
            dE_cur_help[2] = inner_prod(w, dw_dr) + inner_prod(a3_KL, dw_dr) + inner_prod(da3_KL_dr, w) + 4.0*mZeta * thickness / 2.0*wbar*(inner_prod(a3_KL, dw_dr) + inner_prod(da3_KL_dr, a3_KL) + inner_prod(w, da3_KL_dr) + inner_prod(w, dw_dr)) + 2.0*mZeta * thickness / 2.0*dwbar_dr*(inner_prod(a3_KL, a3_KL) + 2.0*inner_prod(w, a3_KL) + inner_prod(w, w)); // FML: bei Oesterle zusätzliche t/2

            // Commented terms multiplied with 'w' are moved to the other loop where the 'dda3_KL_dT1dr' and 'dda3_KL_dT2dr' are computed above.
            rB(0, count1) += mInitialTransConToCar(0, 0) * dE_cur_help[0];
            rB(1, count1) += mInitialTransConToCar(1, 0) * dE_cur_help[0] + mInitialTransConToCar(1, 1) * dE_cur_help[1] + mInitialTransConToCar(1, 3) * dE_cur_help[3];
            rB(2, count1) += dE_cur_help[2];
            rB(3, count1) += mInitialTransConToCar(3, 0) * dE_cur_help[0] + mInitialTransConToCar(3, 3) * dE_cur_help[3];
            rB(4, count1) += mInitialTransConToCar(4, 4) * dE_cur_help[4] + mInitialTransConToCar(4, 5) * dE_cur_help[5];
            rB(5, count1) += mInitialTransConToCar(5, 5) * dE_cur_help[5];


            count1 ++;
            }
            count2 ++;
        }

        // 1. First strain variation
        
        for (unsigned int r = 0; r < mat_size; r++){
            // local node number kr and dof direction dirr
            int kr = r / 7;
            int dirr = r % 7;
            
            array_1d<double, 3> Da1_Dr = ZeroVector(3);
            array_1d<double, 3> Da2_Dr = ZeroVector(3);
            array_1d<double, 3> Da3_tilde_Dr = ZeroVector(3);
            array_1d<double, 3> Da3_KL_Dr = ZeroVector(3);
            array_1d<double, 3> Dw_Dr = ZeroVector(3);

            // the three entries E33, E23 and E13 w.r.t. the curvilinear coord. sys. are saved in dE_cur
            array_1d<double, 3> dE_cur = ZeroVector(3);

            if (dirr == 0 || dirr == 1 || dirr == 2){
                array_1d<double, 3> Da1_Drxa2, a1xDa2_Dr;
                Da1_Dr[dirr] = DN_De(kr, 0);
                Da2_Dr[dirr] = DN_De(kr, 1);
                MathUtils<double>::CrossProduct(Da1_Drxa2, Da1_Dr, rActualMetric.a2);
                MathUtils<double>::CrossProduct(a1xDa2_Dr, rActualMetric.a1, Da2_Dr);
                Da3_tilde_Dr = Da1_Drxa2 + a1xDa2_Dr;
                double Da1xa2_Dr_a1xa2 = inner_prod(Da3_tilde_Dr, rActualMetric.a3_KL_tilde);
                Da3_KL_Dr = Da3_tilde_Dr / rActualMetric.dA 
                    - rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 3);
                dE_cur[0] = inner_prod(Da3_KL_Dr, rw);
                dE_cur[1] = 0.5 * (rw(dirr) * DN_De(kr, 1));
                dE_cur[2] = 0.5 * (rw(dirr) * DN_De(kr, 0));
            }
            else if(dirr == 3){
                Dw_Dr[0] = N(kr);
                dE_cur[0] = N(kr) * (rw(0) + rActualMetric.a3_KL(0));
            }
            else if(dirr == 4){
                Dw_Dr[1] = N(kr);
                dE_cur[0] = N(kr) * (rw(1) + rActualMetric.a3_KL(1));               
            }
            else if(dirr == 5){
                Dw_Dr[2] = N(kr);
                dE_cur[0] = N(kr) * (rw(2) + rActualMetric.a3_KL(2));
            }
            if(dirr == 3 || dirr == 4 || dirr == 5){
                dE_cur[1] += 0.5 * (Dw_Dr[0] * rActualMetric.a2(0) + Dw_Dr[1] * rActualMetric.a2(1) + Dw_Dr[2] * rActualMetric.a2(2));
                dE_cur[2] += 0.5 * (Dw_Dr[0] * rActualMetric.a1(0) + Dw_Dr[1] * rActualMetric.a1(1) + Dw_Dr[2] * rActualMetric.a1(2));
            }

            // calculated with the simplified T_con_to_car (ML)
            // rB(2, r) += dE_cur[0];
            // rB(4, r) += mInitialTransConToCar(4, 4) * dE_cur[1] + mInitialTransConToCar(4, 5) * dE_cur[2];
            // rB(5, r) += mInitialTransConToCar(5, 5) * dE_cur[2];
            // the other entries are (remain) zero
        
            // 2. First curvature variation
            array_1d<double, 3> DDw_DD1r = ZeroVector(3);
            array_1d<double, 3> DDw_DD2r = ZeroVector(3);
            array_1d<double, 3> DDa1_DD1r = ZeroVector(3);
            array_1d<double, 3> DDa1_DD2r = ZeroVector(3);
            array_1d<double, 3> DDa2_DD2r = ZeroVector(3);
            array_1d<double, 3> DDa3_KL_DD1r = ZeroVector(3);
            array_1d<double, 3> DDa3_KL_DD2r = ZeroVector(3);
            array_1d<double, 3> DDa3_tilde_DD1r = ZeroVector(3);
            array_1d<double, 3> DDa3_tilde_DD2r = ZeroVector(3);
            array_1d<double, 3> Da1_D1xa2, a1xDa2_D1, Da1_D2xa2, a1xDa2_D2;
            array_1d<double, 6> dK_cu = ZeroVector(6);

            MathUtils<double>::CrossProduct(Da1_D1xa2, rActualMetric.Da1_D1, rActualMetric.a2);
            MathUtils<double>::CrossProduct(a1xDa2_D1, rActualMetric.a1, rActualMetric.Da1_D2);
            MathUtils<double>::CrossProduct(Da1_D2xa2, rActualMetric.Da1_D2, rActualMetric.a2);
            MathUtils<double>::CrossProduct(a1xDa2_D2, rActualMetric.a1, rActualMetric.Da2_D2);

            if (dirr == 0 || dirr == 1 || dirr == 2){
                DDa1_DD1r[dirr] = DDN_DDe(kr, 0);
                DDa1_DD2r[dirr] = DDN_DDe(kr, 2);
                DDa2_DD2r[dirr] = DDN_DDe(kr, 1);
                array_1d<double, 3> DDa1_DD1rxa2, Da1_D1xDa2_Dr, Da1_DrxDa2_D1, a1xDDa2_DD1r, DDa1_DD2rxa2, Da1_D2xDa2_Dr, Da1_DrxDa2_D2, 
                    a1xDDa2_DD2r, Da1_Drxa2, a1xDa2_Dr;
                MathUtils<double>::CrossProduct(DDa1_DD1rxa2, DDa1_DD1r, rActualMetric.a2);
                MathUtils<double>::CrossProduct(Da1_D1xDa2_Dr, rActualMetric.Da1_D1, Da2_Dr);
                MathUtils<double>::CrossProduct(Da1_DrxDa2_D1, Da1_Dr, rActualMetric.Da1_D2);
                MathUtils<double>::CrossProduct(a1xDDa2_DD1r, rActualMetric.a1, DDa1_DD2r);
                MathUtils<double>::CrossProduct(DDa1_DD2rxa2, DDa1_DD2r, rActualMetric.a2);
                MathUtils<double>::CrossProduct(Da1_D2xDa2_Dr, rActualMetric.Da1_D2, Da2_Dr);
                MathUtils<double>::CrossProduct(Da1_DrxDa2_D2, Da1_Dr, rActualMetric.Da2_D2);
                MathUtils<double>::CrossProduct(a1xDDa2_DD2r, rActualMetric.a1, DDa2_DD2r);
                MathUtils<double>::CrossProduct(Da1_Drxa2, Da1_Dr, rActualMetric.a2);
                MathUtils<double>::CrossProduct(a1xDa2_Dr, rActualMetric.a1, Da2_Dr);
                DDa3_tilde_DD1r = DDa1_DD1rxa2 + Da1_D1xDa2_Dr + Da1_DrxDa2_D1 + a1xDDa2_DD1r;
                DDa3_tilde_DD2r = DDa1_DD2rxa2 + Da1_D2xDa2_Dr + Da1_DrxDa2_D2 + a1xDDa2_DD2r;
                // FML: compare calculation of DDa3_KL_DD1r with "Schale_NURBS_7P.cpp, line 2405" (same for DDa3_KL_DD1s)
                DDa3_KL_DD1r = DDa3_tilde_DD1r / rActualMetric.dA - 
                    2.0 * (Da1_D1xa2 + a1xDa2_D1) * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 3) -
                    rActualMetric.a3_KL_tilde * (inner_prod(Da3_tilde_Dr, (Da1_D1xa2 + a1xDa2_D1)) 
                    + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD1r)) / 
                    pow(rActualMetric.dA, 3) + 
                    3.0 * rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, (Da1_D1xa2 + a1xDa2_D1)) * 
                    inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 5);
                DDa3_KL_DD2r = DDa3_tilde_DD2r / rActualMetric.dA - 
                    2.0 * (Da1_D2xa2 + a1xDa2_D2) * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 3) -
                    rActualMetric.a3_KL_tilde * (inner_prod(Da3_tilde_Dr, (Da1_D2xa2 + a1xDa2_D2)) 
                    + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD2r)) / 
                    pow(rActualMetric.dA, 3) + 
                    3.0 * rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, (Da1_D2xa2 + a1xDa2_D2)) * 
                    inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 5);
                dK_cu[0] = rDw_D1[dirr] * DN_De(kr, 0);
                dK_cu[1] = rDw_D2[dirr] * DN_De(kr, 1);
                dK_cu[2] = 4.0 * rw_bar * (inner_prod(Da3_KL_Dr, rActualMetric.a3_KL) + inner_prod(rw, Da3_KL_Dr));
                dK_cu[3] = 0.5 * (rDw_D1[dirr] * DN_De(kr, 1) + rDw_D2[dirr] * DN_De(kr, 0));
                dK_cu[4] = 0.5 * (inner_prod(rDw_D2, Da3_KL_Dr) + inner_prod(rw, DDa3_KL_DD2r)) + rw_bar * DN_De(kr, 1) * rw[dirr];
                dK_cu[5] = 0.5 * (inner_prod(rDw_D1, Da3_KL_Dr) + inner_prod(rw, DDa3_KL_DD1r)) + rw_bar * DN_De(kr, 0) * rw[dirr];
            }
            if (dirr == 3 || dirr == 4 || dirr == 5){
                DDw_DD1r[dirr - 3] = DN_De(kr, 0);
                DDw_DD2r[dirr - 3] = DN_De(kr, 0);
                dK_cu[2] = 4.0 * rw_bar * (N[kr] * (rActualMetric.a3_KL[dirr - 3] + rw[dirr - 3]));
                dK_cu[4] = 0.5 * (DN_De(kr, 1) * (rActualMetric.a3_KL[dirr - 3] + rw[dirr - 3]) 
                    + N[kr] * (rActualMetric.Da3_KL_D2[dirr - 3] + rDw_D2[dirr - 3])) + rw_bar * N[kr] * rActualMetric.a2[dirr - 3];
                dK_cu[5] = 0.5 * (DN_De(kr, 0) * (rActualMetric.a3_KL[dirr - 3] + rw[dirr - 3]) 
                    + N[kr] * (rActualMetric.Da3_KL_D1[dirr - 3] + rDw_D1[dirr - 3])) + rw_bar * N[kr] * rActualMetric.a1[dirr - 3];
            }
            else if (dirr == 6){
                    dK_cu[2] = 2.0 * N[kr] * (inner_prod(rActualMetric.a3_KL, rActualMetric.a3_KL) 
                        + 2.0 * inner_prod(rw, rActualMetric.a3_KL) + inner_prod(rw, rw));
                    dK_cu[4] = N[kr] * inner_prod(rActualMetric.a2, rw);
                    dK_cu[5] = N[kr] * inner_prod(rActualMetric.a1, rw);
            }
            dK_cu[0] += inner_prod(DDw_DD1r, rActualMetric.a1);
            dK_cu[1] += inner_prod(DDw_DD2r, rActualMetric.a2);
            dK_cu[3] += 0.5 * (inner_prod(DDw_DD1r, rActualMetric.a2) + inner_prod(DDw_DD2r, rActualMetric.a1));

            // calculated with simplified T_con_to_car (ML)
            // rB(0, r) += mZeta * thickness / 2.0 * mInitialTransConToCar(0, 0) * dK_cu[0];
            // rB(1, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(1, 0) * dK_cu[0] + mInitialTransConToCar(1, 1) * dK_cu[1] 
            //     + mInitialTransConToCar(1, 3) * dK_cu[3]);
            // rB(2, r) += mZeta * thickness / 2.0 * dK_cu[2];
            // rB(3, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(3, 0) * dK_cu[0] + mInitialTransConToCar(3, 3) * dK_cu[3]);
            // rB(4, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(4, 4) * dK_cu[4] + mInitialTransConToCar(4, 5) * dK_cu[5]);
            // rB(5, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(5, 5) * dK_cu[5]);

            // 3. Second Strain Variation
            if (rCalculateStiffnessMatrixFlag == true){
                for (unsigned int s = 0; s <= r; s++)
                {
                    // local node number ks and dof direction dirs
                    int ks = s / 7;
                    int dirs = s % 7;
                    
                    array_1d<double, 3> Da1_Ds = ZeroVector(3);
                    array_1d<double, 3> Da2_Ds = ZeroVector(3);
                    array_1d<double, 3> Da1_Dsxa2, a1xDa2_Ds;
                    array_1d<double, 3> Da3_tilde_Ds = ZeroVector(3);
                    array_1d<double, 3> DDa3_tilde_DDrs = ZeroVector(3);
                    array_1d<double, 3> DDa3_KL_DDrs = ZeroVector(3);
                    array_1d<double, 3> Dw_Ds = ZeroVector(3);
                    array_1d<double, 3> ddE_cu = ZeroVector(3);
                    double DdA_Ds = 0.0;
                    double DdA_Dr = 0.0;

                    if (dirr == 0){
                        if (dirs == 1)
                            DDa3_tilde_DDrs[2] = DN_De(kr, 0) * DN_De(ks, 1) - DN_De(ks, 0) * DN_De(kr, 1);
                        else if (dirs == 2)
                            DDa3_tilde_DDrs[1] = -DN_De(kr, 0) * DN_De(ks, 1) + DN_De(ks, 0) * DN_De(kr, 1);
                    }
                    else if (dirr == 1){
                        if (dirs == 0)
                            DDa3_tilde_DDrs[2] = -DN_De(kr, 0) * DN_De(ks, 1) + DN_De(ks, 0) * DN_De(kr, 1);
                        else if (dirs == 2)
                            DDa3_tilde_DDrs[0] = DN_De(kr, 0) * DN_De(ks, 1) - DN_De(ks, 0) * DN_De(kr, 1);
                    }
                    else if (dirr == 2){
                        if (dirs == 0)
                            DDa3_tilde_DDrs[1] = DN_De(kr, 0) * DN_De(ks, 1) - DN_De(ks, 0) * DN_De(kr, 1);
                        else if (dirs == 1)
                            DDa3_tilde_DDrs[0] = -DN_De(kr, 0) * DN_De(ks, 1) + DN_De(ks, 0) * DN_De(kr, 1);
                    }
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        Da1_Ds[dirs] = DN_De(ks, 0);
                        Da2_Ds[dirs] = DN_De(ks, 1);
                        MathUtils<double>::CrossProduct(Da1_Dsxa2, Da1_Ds, rActualMetric.a2);
                        MathUtils<double>::CrossProduct(a1xDa2_Ds, rActualMetric.a1, Da2_Ds);
                        Da3_tilde_Ds = Da1_Dsxa2 + a1xDa2_Ds;
                        DdA_Ds = inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / rActualMetric.dA;
                        DDa3_KL_DDrs += -DdA_Ds * Da3_tilde_Dr / (rActualMetric.dA * rActualMetric.dA);
                        ddE_cu[1] += 0.5 * (Dw_Dr[dirs] * DN_De(ks, 1));
                        ddE_cu[2] += 0.5 * (Dw_Dr[dirs] * DN_De(ks, 0));
                        if (dirr == 0 || dirr == 1 || dirr == 2)
                            DDa3_tilde_DDrs = CrossProduct(Da1_Dr, Da2_Ds) + CrossProduct(Da1_Ds, Da2_Dr);
                        else if (dirr == 3 || dirr == 4 || dirr == 5)
                            ddE_cu[0] += DdA_Ds * N[kr];
                    }
                    if (dirs == 3 || dirs == 4 || dirs == 5)
                        Dw_Ds[dirs - 3] = N(ks);
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        DdA_Dr = inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / rActualMetric.dA;
                        DDa3_KL_DDrs += -DdA_Dr * Da3_tilde_Ds / (rActualMetric.dA * rActualMetric.dA);
                        if (dirs == 0 || dirs == 1 || dirs == 2){
                            DDa3_KL_DDrs += DDa3_tilde_DDrs / rActualMetric.dA - 
                                (inner_prod(DDa3_tilde_DDrs, rActualMetric.a3_KL_tilde) + 
                                inner_prod(Da3_tilde_Dr, Da3_tilde_Ds)) * rActualMetric.a3_KL_tilde / rActualMetric.dA +
                                3.0 * DdA_Dr * DdA_Ds * rActualMetric.a3_KL_tilde / pow(rActualMetric.dA, 3);
                        }
                        else if (dirs == 3 || dirs == 4 || dirs == 5)
                            ddE_cu[0] += DdA_Dr * N[ks];
                        ddE_cu[1] = 0.5 * (Dw_Ds[dirr] * DN_De(kr, 1));
                        ddE_cu[2] = 0.5 * (Dw_Ds[dirr] * DN_De(kr, 0));
                    }
                    if ((dirr == 3 || dirr == 4 || dirr == 5) && (dirs == 3 || dirs == 4 || dirs ==5))
                        ddE_cu[0] += N[kr] * N[ks];
                    ddE_cu[0] += inner_prod(DDa3_KL_DDrs, rw);

                    // calculated with simplified T_con_to_car (ML)
                    rSecondVariations.B33(r, s) += ddE_cu[0];
                    rSecondVariations.B23(r, s) += mInitialTransConToCar(4, 4) * ddE_cu[1] + mInitialTransConToCar(4, 5) * ddE_cu[2];
                    rSecondVariations.B13(r, s) += mInitialTransConToCar(5, 5) * ddE_cu[2];
                    if (r != s){
                        rSecondVariations.B33(s, r) += rSecondVariations.B33(r, s);
                        rSecondVariations.B23(s, r) += rSecondVariations.B23(r, s);
                        rSecondVariations.B13(s, r) += rSecondVariations.B13(r, s);
                    }
                    // all other entries are (remain) zero
                    
                    // 4. Second curvature variation
                    array_1d<double, 3> DDa1_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDa1_DD2s = ZeroVector(3);
                    array_1d<double, 3> DDa2_DD2s = ZeroVector(3);
                    array_1d<double, 3> Da3_KL_Ds = ZeroVector(3);
                    array_1d<double, 3> DDa3_KL_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDa3_KL_DD2s = ZeroVector(3);
                    array_1d <double, 3> DDDa3_KL_DDD1rs = ZeroVector(3);
                    array_1d <double, 3> DDDa3_KL_DDD2rs = ZeroVector(3);
                    array_1d<double, 3> DDa3_tilde_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDa3_tilde_DD2s = ZeroVector(3);
                    array_1d<double, 3> DDDa3_tilde_DDD1rs = ZeroVector(3);
                    array_1d<double, 3> DDDa3_tilde_DDD2rs = ZeroVector(3);
                    double Da3_bar_Dr = 0.0;
                    double DDa3_bar_DD1r = 0.0;
                    double DDa3_bar_DD2r = 0.0;

                    array_1d<double, 3> DDw_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDw_DD2s = ZeroVector(3);
                    double Dwbar_Dr = 0.0;
                    double Dwbar_Ds = 0.0;
                    double DDdA_DDrs, DDDdA_DDD1rs, DDDdA_DDD2rs, dP1dl, dR1dl, dS1dl, dP2dl, dR2dl, dS2dl;
                    
                    array_1d<double, 3> M1 = ZeroVector(3);
                    array_1d<double, 3> T1 = ZeroVector(3);
                    array_1d<double, 3> M2 = ZeroVector(3);
                    array_1d<double, 3> T2 = ZeroVector(3);
                    array_1d<double, 3> dM1_ds = ZeroVector(3);
                    array_1d<double, 3> dT1_ds = ZeroVector(3);
                    array_1d<double, 3> dM2_ds = ZeroVector(3);
                    array_1d<double, 3> dT2_ds = ZeroVector(3);
                    array_1d <double, 6> ddK_cu = ZeroVector(6);

                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        if(dirr == 0 || dirr == 1 || dirr == 2){
                            DDa1_DD1s[dirr] = DDN_DDe(kr, 0);
                            DDa1_DD2s[dirr] = DDN_DDe(kr, 2);
                            DDa2_DD2s[dirr] = DDN_DDe(kr, 1);
                        }
                        array_1d<double, 3> DDa1_DD1sxa2, Da1_D1xDa2_Ds, Da1_DsxDa2_D1, a1xDDa2_DD1s, DDa1_DD2sxa2, Da1_D2xDa2_Ds, 
                            Da1_DsxDa2_D2, a1xDDa2_DD2s;
                        MathUtils<double>::CrossProduct(DDa1_DD1sxa2, DDa1_DD1s, rActualMetric.a2);
                        MathUtils<double>::CrossProduct(Da1_D1xDa2_Ds, rActualMetric.Da1_D1, Da2_Ds);
                        MathUtils<double>::CrossProduct(Da1_DsxDa2_D1, Da1_Ds, rActualMetric.Da1_D2);
                        MathUtils<double>::CrossProduct(a1xDDa2_DD1s, rActualMetric.a1, DDa1_DD2s);
                        MathUtils<double>::CrossProduct(DDa1_DD2sxa2, DDa1_DD2s, rActualMetric.a2);
                        MathUtils<double>::CrossProduct(Da1_D2xDa2_Ds, rActualMetric.Da1_D2, Da2_Ds);
                        MathUtils<double>::CrossProduct(Da1_DsxDa2_D2, Da1_Ds, rActualMetric.Da2_D2);
                        MathUtils<double>::CrossProduct(a1xDDa2_DD2s, rActualMetric.a1, DDa2_DD2s);
                        DDa3_tilde_DD1s = DDa1_DD1sxa2 + Da1_D1xDa2_Ds + Da1_DsxDa2_D1 + a1xDDa2_DD1s;
                        DDa3_tilde_DD2s = DDa1_DD2sxa2 + Da1_D2xDa2_Ds + Da1_DsxDa2_D2 + a1xDDa2_DD2s;
                        DDa3_KL_DD1s = DDa3_tilde_DD1s / rActualMetric.dA - 
                            2.0 * (Da1_D1xa2 + a1xDa2_D1) * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 3) -
                            rActualMetric.a3_KL_tilde * (inner_prod(Da3_tilde_Ds, (Da1_D1xa2 + a1xDa2_D1)) 
                            + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD1s)) / 
                            pow(rActualMetric.dA, 3) + 
                            3.0 * rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, (Da1_D1xa2 + a1xDa2_D1)) * 
                            inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 5);
                        DDa3_KL_DD2s = DDa3_tilde_DD2s / rActualMetric.dA - 
                            2.0 * (Da1_D2xa2 + a1xDa2_D2) * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 3) -
                            rActualMetric.a3_KL_tilde * (inner_prod(Da3_tilde_Ds, (Da1_D2xa2 + a1xDa2_D2)) 
                            + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD2s)) / 
                            pow(rActualMetric.dA, 3) + 
                            3.0 * rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, (Da1_D2xa2 + a1xDa2_D2)) * 
                            inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 5);
                    }
                    else if (dirs == 3 || dirs == 4 || dirs == 5){
                        DDw_DD1s[dirs - 3] += DN_De(ks, 0);
                        DDw_DD2s[dirs - 3] += DN_De(ks, 0);
                    }
                    else if (dirs == 6){
                        Dwbar_Ds = N(ks);
                    }
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        Da3_KL_Ds = Da3_tilde_Ds / rActualMetric.dA 
                            - rActualMetric.a3_KL_tilde * inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 3);
                        Da3_bar_Dr = inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / rActualMetric.dA;
                        DDa3_bar_DD1r = (inner_prod(Da3_tilde_Dr, rActualMetric.Da3_KL_tilde_D1) + 
                            inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD1r)) / rActualMetric.dA -
                            inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1) * 
                            inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 3);
                        DDa3_bar_DD2r = (inner_prod(Da3_tilde_Dr, rActualMetric.Da3_KL_tilde_D2) + 
                            inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD2r)) / rActualMetric.dA -
                            inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D2) * 
                            inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) / pow(rActualMetric.dA, 3);
                        if (dirs == 0 || dirs == 1 || dirs == 2){
                            double Da3_bar_D1 = inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1) / rActualMetric.dA;
                            double Da3_bar_D2 = inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D2) / rActualMetric.dA;
                            double Da3_bar_Ds = inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / rActualMetric.dA;
                            double DDa3_bar_DD1s = (inner_prod(Da3_tilde_Ds, rActualMetric.Da3_KL_tilde_D1) + 
                                inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD1s)) / rActualMetric.dA -
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1) * 
                                inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 3);
                            double DDa3_bar_DD2s = (inner_prod(Da3_tilde_Ds, rActualMetric.Da3_KL_tilde_D2) + 
                                inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD2s)) / rActualMetric.dA -
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D2) * 
                                inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) / pow(rActualMetric.dA, 3);
                            // Single terms of DDa3_KL_DD1alpha
                            M1 = DDa3_tilde_DD1r*rActualMetric.dA + rActualMetric.Da3_KL_tilde_D1*Da3_bar_Dr - Da3_tilde_Dr*Da3_bar_D1 - 
                                rActualMetric.a3_KL_tilde*DDa3_bar_DD1r;
                            T1 = rActualMetric.Da3_KL_tilde_D1*rActualMetric.dA - rActualMetric.a3_KL_tilde*Da3_bar_D1;
                            M2 = DDa3_tilde_DD2r*rActualMetric.dA + rActualMetric.Da3_KL_tilde_D2*Da3_bar_Dr - Da3_tilde_Dr*Da3_bar_D2 - 
                                rActualMetric.a3_KL_tilde*DDa3_bar_DD2r;
                            T2 = rActualMetric.Da3_KL_tilde_D2*rActualMetric.dA - rActualMetric.a3_KL_tilde*Da3_bar_D2;
                            // expressions a3_KL_tilde,rs and a3_KL_tilde,alpha,rs are needed
                            DDDa3_tilde_DDD1rs = CrossProduct(DDa1_DD1r, Da2_Ds) + CrossProduct(DDa1_DD1s, Da2_Dr) + 
                                CrossProduct(Da1_Dr, DDa1_DD2s) + CrossProduct(Da1_Ds,DDa1_DD2r);
                            DDDa3_tilde_DDD2rs = CrossProduct(DDa1_DD2r, Da2_Ds) + CrossProduct(DDa1_DD2s, Da2_Dr) + 
                                CrossProduct(Da1_Dr, DDa2_DD2s) + CrossProduct(Da1_Ds, DDa2_DD2r);
                            // expressions dA,rs and dA,alpha,rs are needed
                            DDdA_DDrs = -inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Ds) * 
                                inner_prod(rActualMetric.a3_KL_tilde,Da3_tilde_Dr) / pow(rActualMetric.dA, 3) + 
                                (inner_prod(Da3_tilde_Ds, Da3_tilde_Dr) + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DDrs)) /
                                rActualMetric.dA;
                            dP1dl = 3.0*Da3_bar_Ds*(inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) *
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1))/pow(rActualMetric.dA, 4) - 
                                (inner_prod(Da3_tilde_Ds, Da3_tilde_Dr) * inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1) 
                                + inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DDrs) * 
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D1)+
                                inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr)*inner_prod(Da3_tilde_Ds,rActualMetric.Da3_KL_tilde_D1) + 
                                inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr)*inner_prod(rActualMetric.a3_KL_tilde, DDa3_tilde_DD1s)) /
                                pow(rActualMetric.dA, 3);
                            dR1dl = -Da3_bar_Ds*inner_prod(Da3_tilde_Dr,rActualMetric.Da3_KL_tilde_D1)/(rActualMetric.dA*rActualMetric.dA) + 
                                (inner_prod(DDa3_tilde_DDrs, rActualMetric.Da3_KL_tilde_D1) + inner_prod(Da3_tilde_Dr, DDa3_tilde_DD1s)) /
                                rActualMetric.dA;
                            dS1dl = -Da3_bar_Ds*inner_prod(rActualMetric.a3_KL_tilde,DDa3_tilde_DD1r)/(rActualMetric.dA*rActualMetric.dA)+ 
                                (inner_prod(Da3_tilde_Ds,DDa3_tilde_DD1r) + inner_prod(rActualMetric.a3_KL_tilde, DDDa3_tilde_DDD1rs)) /
                                rActualMetric.dA;
                            dP2dl = 3.0*Da3_bar_Ds*(inner_prod(rActualMetric.a3_KL_tilde, Da3_tilde_Dr) *
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D2))/ pow(rActualMetric.dA, 4) - 
                                (inner_prod(Da3_tilde_Ds, Da3_tilde_Dr)*inner_prod(rActualMetric.a3_KL_tilde,rActualMetric.Da3_KL_tilde_D2) + 
                                inner_prod(rActualMetric.a3_KL_tilde,DDa3_tilde_DDrs) * 
                                inner_prod(rActualMetric.a3_KL_tilde, rActualMetric.Da3_KL_tilde_D2) +
                                inner_prod(rActualMetric.a3_KL_tilde,Da3_tilde_Dr)*inner_prod(Da3_tilde_Ds,rActualMetric.Da3_KL_tilde_D2) + 
                                inner_prod(rActualMetric.a3_KL_tilde,Da3_tilde_Dr)*inner_prod(rActualMetric.a3_KL_tilde,DDa3_tilde_DD2s)) / 
                                pow(rActualMetric.dA, 3);
                            dR2dl = -Da3_bar_Ds*inner_prod(Da3_tilde_Dr,rActualMetric.Da3_KL_tilde_D2)/(rActualMetric.dA*rActualMetric.dA) + 
                                (inner_prod(DDa3_tilde_DDrs,rActualMetric.Da3_KL_tilde_D2) + 
                                inner_prod(Da3_tilde_Dr,DDa3_tilde_DD2s)) /
                                rActualMetric.dA;
                            dS2dl = -Da3_bar_Ds*inner_prod(rActualMetric.a3_KL_tilde,DDa3_tilde_DD2r)/(rActualMetric.dA*rActualMetric.dA) + 
                                (inner_prod(Da3_tilde_Ds,DDa3_tilde_DD2r) + inner_prod(rActualMetric.a3_KL_tilde, DDDa3_tilde_DDD2rs)) /
                                rActualMetric.dA;
                            DDDdA_DDD1rs = dP1dl + dR1dl + dS1dl;
                            DDDdA_DDD2rs = dP2dl + dR2dl + dS2dl;
                            // M1,s &  T1,s & M2,s & T2,s
                            dM1_ds = DDDa3_tilde_DDD1rs*rActualMetric.dA + DDa3_tilde_DD1r*Da3_bar_Ds + DDa3_tilde_DD1s*Da3_bar_Dr + 
                                rActualMetric.Da3_KL_tilde_D1*DDdA_DDrs - DDa3_tilde_DDrs*Da3_bar_D1 - Da3_tilde_Dr*DDa3_bar_DD1s - 
                                Da3_tilde_Ds*DDa3_bar_DD1r - rActualMetric.a3_KL_tilde*DDDdA_DDD1rs;
                            dT1_ds = DDa3_tilde_DD1s*rActualMetric.dA + rActualMetric.Da3_KL_tilde_D1*Da3_bar_Ds - Da3_tilde_Ds*Da3_bar_D1 - 
                                rActualMetric.a3_KL_tilde*DDa3_bar_DD1s;
                            dM2_ds = DDDa3_tilde_DDD2rs*rActualMetric.dA + DDa3_tilde_DD2r*Da3_bar_Ds + DDa3_tilde_DD2s*Da3_bar_Dr + 
                                rActualMetric.Da3_KL_tilde_D2*DDdA_DDrs - DDa3_tilde_DDrs*Da3_bar_D2 - Da3_tilde_Dr*DDa3_bar_DD2s - 
                                Da3_tilde_Ds*DDa3_bar_DD2r - rActualMetric.a3_KL_tilde*DDDdA_DDD2rs;
                            dT2_ds = DDa3_tilde_DD2s*rActualMetric.dA + rActualMetric.Da3_KL_tilde_D2*Da3_bar_Ds - Da3_tilde_Ds*Da3_bar_D2 - 
                                rActualMetric.a3_KL_tilde*DDa3_bar_DD2s;
                            // a3_KL,alpha,rs
                            DDDa3_KL_DDD1rs = ((dM1_ds*rActualMetric.dA*rActualMetric.dA + 2.0*M1*rActualMetric.dA*Da3_bar_Ds - 
                                2.0*dT1_ds*rActualMetric.dA*Da3_bar_Dr - 2.0*T1*Da3_bar_Ds*Da3_bar_Dr - 2.0*T1*rActualMetric.dA*DDdA_DDrs) -
                                4.0*Da3_bar_Ds*(M1*rActualMetric.dA - 2.0*T1*Da3_bar_Dr))/ pow(rActualMetric.dA, 4);
                            DDDa3_KL_DDD2rs = ((dM2_ds*rActualMetric.dA*rActualMetric.dA + 2.0*M2*rActualMetric.dA*Da3_bar_Ds - 
                                2.0*dT2_ds*rActualMetric.dA*Da3_bar_Dr - 2.0*T2*Da3_bar_Ds*Da3_bar_Dr - 2.0*T2*rActualMetric.dA*DDdA_DDrs) -
                                4.0*Da3_bar_Ds*(M2*rActualMetric.dA - 2.0*T2*Da3_bar_Dr))/ pow(rActualMetric.dA, 4);
                        }
                        ddK_cu[0] = DDw_DD1s[dirr] * DN_De(kr, 0);
                        ddK_cu[1] = DDw_DD2s[dirr] * DN_De(kr, 1);
                        ddK_cu[3] = 0.5 * (DDw_DD1s[dirr] * DN_De(kr, 1) + DDw_DD2s[dirr] * DN_De(kr, 0));
                    }
                    else if (dirr == 6)
                        Dwbar_Dr = N(kr);
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        ddK_cu[0] += DDw_DD1r[dirs] * DN_De(ks, 0);
                        ddK_cu[1] += DDw_DD2r[dirs] * DN_De(ks, 1);
                        ddK_cu[3] += 0.5 * (DDw_DD1r[dirs] * DN_De(ks, 1) + DDw_DD2r[dirs] * DN_De(ks, 0));
                    }
                    ddK_cu[2] += 4.0 * Dwbar_Dr * (inner_prod(rActualMetric.a3_KL, Da3_KL_Ds) + inner_prod(Da3_KL_Ds, rw) +
                        inner_prod(rActualMetric.a3_KL, Dw_Ds) + inner_prod(rw, Dw_Ds)) +
                        4.0 * Dwbar_Ds * (inner_prod(Da3_KL_Dr, rActualMetric.a3_KL) + inner_prod(Da3_KL_Dr, rw) +
                        inner_prod(rActualMetric.a3_KL, Dw_Dr) + inner_prod(rw, Dw_Dr)) +
                        4.0 * rw_bar * (inner_prod(DDa3_KL_DDrs, rActualMetric.a3_KL) + inner_prod(Da3_KL_Dr, Da3_KL_Ds) +
                        inner_prod(Dw_Dr, Da3_KL_Ds) + inner_prod(Dw_Ds, Da3_KL_Dr) + inner_prod(DDa3_KL_DDrs, rw) + 
                        inner_prod(Dw_Ds, Dw_Dr));
                    ddK_cu[4] += 0.5 * (inner_prod(DDDa3_KL_DDD2rs, rw) + inner_prod(DDa3_KL_DD2r, Dw_Ds) + 
                        inner_prod(DDa3_KL_DD2s, Dw_Dr) + inner_prod(DDa3_KL_DDrs, rDw_D2) + inner_prod(Da3_KL_Dr, DDw_DD2s) + 
                        inner_prod(Da3_KL_Ds, DDw_DD2r) + inner_prod(DDw_DD2r, Dw_Ds) + inner_prod(DDw_DD2s, Dw_Dr)) +
                        Dwbar_Dr * (inner_prod(Da2_Ds, rw) + inner_prod(rActualMetric.a2, Dw_Ds)) +
                        Dwbar_Ds * (inner_prod(Da2_Dr, rw) + inner_prod(rActualMetric.a2, Dw_Dr)) +
                        rw_bar * (inner_prod(Da2_Dr, Dw_Ds) + inner_prod(Da2_Ds, Dw_Dr));
                    ddK_cu[5] += 0.5 * (inner_prod(DDDa3_KL_DDD1rs, rw) + inner_prod(DDa3_KL_DD1r, Dw_Ds) + 
                        inner_prod(DDa3_KL_DD1s, Dw_Dr) + inner_prod(DDa3_KL_DDrs, rDw_D1) + inner_prod(Da3_KL_Dr, DDw_DD1s) + 
                        inner_prod(Da3_KL_Ds, DDw_DD1r) + inner_prod(DDw_DD1r, Dw_Ds) + inner_prod(DDw_DD1s, Dw_Dr)) +
                        Dwbar_Dr * (inner_prod(Da1_Ds, rw) + inner_prod(rActualMetric.a1, Dw_Ds)) +
                        Dwbar_Ds * (inner_prod(Da1_Dr, rw) + inner_prod(rActualMetric.a1, Dw_Dr)) +
                        rw_bar * (inner_prod(Da1_Dr, Dw_Ds) + inner_prod(Da1_Ds, Dw_Dr));

                    // calculated with simplified T_con_to_car (ML)
                    rSecondVariations.B11(r, s) += mZeta * thickness /2.0 * mInitialTransConToCar(0, 0) * ddK_cu[0];
                    rSecondVariations.B22(r, s) += mZeta * thickness /2.0 * (mInitialTransConToCar(1, 0) * ddK_cu[0] + 
                        mInitialTransConToCar(1, 1) * ddK_cu[1] + mInitialTransConToCar(1, 3) * ddK_cu[3]);
                    rSecondVariations.B33(r, s) += ddK_cu[2];
                    rSecondVariations.B12(r, s) += mZeta * thickness /2.0 * (mInitialTransConToCar(3, 0) * ddK_cu[0] + 
                        mInitialTransConToCar(3, 3) * ddK_cu[3]);
                    rSecondVariations.B23(r, s) += mZeta * thickness / 2.0 * (mInitialTransConToCar(4, 4) * ddK_cu[4] +
                        mInitialTransConToCar(4, 5) * ddK_cu[5]);
                    rSecondVariations.B13(r, s) += mZeta * thickness / 2.0 * mInitialTransConToCar(5, 5) * ddK_cu[5];
                    if (r != s){
                        rSecondVariations.B11(s, r) = rSecondVariations.B11(r, s);
                        rSecondVariations.B22(s, r) = rSecondVariations.B22(r, s);
                        rSecondVariations.B33(s, r) = rSecondVariations.B33(r, s);
                        rSecondVariations.B12(s, r) = rSecondVariations.B12(r, s);
                        rSecondVariations.B23(s, r) = rSecondVariations.B23(r, s);
                        rSecondVariations.B13(s, r) = rSecondVariations.B13(r, s);
                    }
                    // all other entries are (remain) zero
                }   // 2nd for-loop s
            }   // if second variations
        }   //1st for-loop r
    }

    void IgaShell7pElement::CalculateDifferentialVolume(
        double& rdV)
    {
        Vector Dg3_D1 = ZeroVector(3);
        Vector Dg3_D2 = ZeroVector(3);
        Vector Dg1_D1xg2 = ZeroVector(3);
        Vector g1xDg2_D1 = ZeroVector(3);
        Vector Dg1_D2xg2 = ZeroVector(3);
        Vector g1xDg2_D2 = ZeroVector(3);
        Vector G1xG2 = ZeroVector(3);

        array_1d<double, 3> Dg1_D1;
        array_1d<double, 3> Dg2_D2;
        array_1d<double, 3> Dg1_D2;
        for (unsigned int i = 0; i < 3; i++)
        {
            Dg1_D1[i] = mInitialMetric.H(i, 0);
            Dg2_D2[i] = mInitialMetric.H(i, 1);
            Dg1_D2[i] = mInitialMetric.H(i, 2);
        }

        MathUtils<double>::CrossProduct(Dg1_D1xg2, Dg1_D1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(g1xDg2_D1, mInitialMetric.a1, Dg1_D2); // Dg1_D2 = Dg2_D1
        MathUtils<double>::CrossProduct(Dg1_D2xg2, Dg1_D2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(g1xDg2_D2, mInitialMetric.a1, Dg2_D2);
        Dg3_D1 = ((Dg1_D1xg2 + g1xDg2_D1) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(Dg1_D1xg2 + g1xDg2_D1)) 
            / (mInitialMetric.dA * mInitialMetric.dA);
        Dg3_D2 = ((Dg1_D2xg2 + g1xDg2_D2) * mInitialMetric.dA - mInitialMetric.a3_KL_tilde * norm_2(Dg1_D2xg2 + g1xDg2_D2))
            / (mInitialMetric.dA * mInitialMetric.dA);

        // covariant base vectors of the shell body in the reference configuration
        array_1d<double, 3> G1 = mInitialMetric.a1 + mZeta * Dg3_D1;
        array_1d<double, 3> G2 = mInitialMetric.a2 + mZeta * Dg3_D2;
        
        MathUtils<double>::CrossProduct(G1xG2, G1, G2);

        rdV = inner_prod(G1xG2, mInitialMetric.a3_KL);
    }

    void IgaShell7pElement::Calculate(
        const Variable<Vector>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int number_of_control_points = GetGeometry().size();        

        // shape functions
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);

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
        array_1d<double, 3> w_alpha = ZeroVector(3);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);
        // 7th parameter
        double w_bar = 0.0;

        std::vector<array_1d<double, 6>> stress_pk2_cart(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_pk2_cov(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_cau_cov(mGaussQuadratureThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_cau_cart(mGaussQuadratureThickness.num_GP_thickness);

        MetricVariables actual_metric(3, 6);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);
        for (unsigned int i = 0; i < number_of_control_points; ++i){
            w_bar += N(i) * GetGeometry()[i].GetDof(W_BAR).GetSolutionStepValue();
        }
        
        double thickness = GetProperties().GetValue(THICKNESS);

        // the Gauss-Points start from bottom to top
        for (unsigned int Gauss_index = 0; Gauss_index < mGaussQuadratureThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussQuadratureThickness.zeta(Gauss_index);

            array_1d<double, 3> G1(3, 0.0), G2(3, 0.0), G1_con(3, 0.0), G2_con(3, 0.0), g1(3, 0.0), g2(3, 0.0), g3(3, 0.0);
            Matrix F = ZeroMatrix(3, 3);
            double detF = 0.0;
            CalculateInitialBaseVectorsLinearised(G1, G2, G1_con, G2_con);
            CalculateActualBaseVectorsLinearised(actual_metric, w, Dw_D1, Dw_D2, w_bar, g1, g2, g3);
            CalculateDeformationGradient(G1, G2, g1, g2, g3, F, detF);

            // TransformationMatrices
            Matrix initial_trans_car_to_cov = ZeroMatrix(5, 5);
            Matrix actual_trans_cov_to_car = ZeroMatrix(5, 5);
            CalculateTransformationMatrixInitialTransConToCar(G1_con, G2_con);
            CalculateTransformationMatrixInitialTransCarToCov(initial_trans_car_to_cov, G1_con, G2_con);
            CalculateTransformationMatrixActualTransCovToCar(actual_trans_cov_to_car, g1, g2, g3, actual_metric.a2_con, actual_metric.a3_KL);

            ConstitutiveVariables constitutive_variables(6);
            CalculateConstitutiveVariables(actual_metric, w, 
                Dw_D1, Dw_D2, w_bar, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

            // stresses at GP
            stress_pk2_cart[Gauss_index] = constitutive_variables.S;
            stress_pk2_cov[Gauss_index] = prod(initial_trans_car_to_cov, stress_pk2_cart[Gauss_index]);
            stress_cau_cov[Gauss_index] = stress_pk2_cov[Gauss_index] / detF;
            stress_cau_cart[Gauss_index] = prod(actual_trans_cov_to_car, stress_cau_cov[Gauss_index]);
        }
    
        // Cauchy stress at midspan
        array_1d<double, 6> stress_cau_cart_mid;
        stress_cau_cart_mid = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1] + stress_cau_cart[0]) / 2.0;
        
        if (rVariable == STRESS_CAUCHY_TOP_11){
            rValues = stress_cau_cart_mid[0] + (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][0] 
                - stress_cau_cart_mid[0]) / mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1);
        }
        else if (rVariable == STRESS_CAUCHY_TOP_22){
            rValues = stress_cau_cart_mid[1] + (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][1] 
            - stress_cau_cart_mid[1]) / mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1);
        }
        else if (rVariable == STRESS_CAUCHY_TOP_12){
            rValues = stress_cau_cart_mid[3] + (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][3] 
            - stress_cau_cart_mid[3]) / mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1);
        }
        else if (rVariable == STRESS_CAUCHY_BOTTOM_11){
            rValues = stress_cau_cart_mid[0] + (stress_cau_cart[0][0] - stress_cau_cart_mid[0]) / 
            mGaussQuadratureThickness.zeta(0);
        }
        else if (rVariable == STRESS_CAUCHY_BOTTOM_22){
            rValues = stress_cau_cart_mid[1] + (stress_cau_cart[0][1] - stress_cau_cart_mid[1]) / 
            mGaussQuadratureThickness.zeta(0);
        }
        else if (rVariable == STRESS_CAUCHY_BOTTOM_12){
            rValues = stress_cau_cart_mid[3] + (stress_cau_cart[0][3] - stress_cau_cart_mid[3]) / 
            mGaussQuadratureThickness.zeta(0);
        }
        else if (rVariable == INTERNAL_FORCE_11)
            rValues = stress_cau_cart_mid[0] * thickness;
        else if (rVariable == INTERNAL_FORCE_22)
            rValues = stress_cau_cart_mid[1] * thickness;
        else if (rVariable == INTERNAL_FORCE_12)
            rValues = stress_cau_cart_mid[3] * thickness;
        else if (rVariable == INTERNAL_MOMENT_11){
            rValues = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][0] - stress_cau_cart_mid[0]) 
                * thickness * thickness / (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        }
        else if (rVariable == INTERNAL_MOMENT_22){
            rValues = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][1] - stress_cau_cart_mid[1]) * thickness * thickness / 
            (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        }
        else if (rVariable == INTERNAL_MOMENT_12){
            rValues = (stress_cau_cart[mGaussQuadratureThickness.num_GP_thickness-1][3] - stress_cau_cart_mid[3]) * thickness * thickness / 
            (mGaussQuadratureThickness.zeta(mGaussQuadratureThickness.num_GP_thickness-1) * 6);
        }
        else if (rVariable == SHEAR_FORCE_1){
            if(mGaussQuadratureThickness.num_GP_thickness != 3){
                KRATOS_WATCH("Calculation formula of SHEAR_FORCE_1 is implemented for 3 Gauss points in thickness direction and can not be used in this case.")
                KRATOS_ERROR << "Calculation formula of SHEAR_FORCE_1 is implemented for 3 Gauss points in thickness direction and can not be used in this case." << std::endl;
            }
            rValues = thickness * stress_cau_cart_mid[5];
        }
        else if (rVariable == SHEAR_FORCE_2){
            if(mGaussQuadratureThickness.num_GP_thickness != 3){
                KRATOS_WATCH("Calculation formula of SHEAR_FORCE_2 is implemented for 3 Gauss points in thickness direction and can not be used in this case.")
                KRATOS_ERROR << "Calculation formula of SHEAR_FORCE_2 is implemented for 3 Gauss points in thickness direction and can not be used in this case." << std::endl;
            }
            rValues = thickness * stress_cau_cart_mid[4];
        }
        else{
            KRATOS_WATCH("No results for desired variable available in Calculate of IgaShell5pElement.")
            rValues = 0.0;
        }
    }

    void IgaShell7pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        // KRATOS_WATCH("EquationIdVector")

        const unsigned int number_of_control_points = GetGeometry().size();

        if (rResult.size() != 7 * number_of_control_points)
            rResult.resize(7 * number_of_control_points, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 7;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
            rResult[index + 6] = GetGeometry()[i].GetDof(W_BAR).EquationId();
        }

        KRATOS_CATCH("")
    }

    void IgaShell7pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        // KRATOS_WATCH("GetDofList")

        const unsigned int number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(7 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(W_BAR));
        }

        KRATOS_CATCH("")
    }

    int IgaShell7pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        if (DISPLACEMENT.Key() == 0)
            KRATOS_ERROR << "DISPLACEMENT has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_VALUES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_VALUES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        // Check whether ConstitutiveLaw is 3D
        if (mConstitutiveLawVector[0]->GetStrainSize() != 6){
            KRATOS_WATCH("ConstitutiveLaw is not 3D.")  // not Kratos style, but output is not shown with KRATOS_ERROR 
            KRATOS_ERROR << "ConstitutiveLaw is not 3D." << std::endl;
        }
        return 0;
    }

} // Namespace Kratos
