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
//                   based on the work of Ralf Echter and Bastian Oesterle on a hierarchic shell formulation
/**
* Description:
* Isogeometric hierarchic Reissner-Mindlin shell element parameterized by 5 parameters (5p)
* The 5 parameters: three translations (u,v,w) and two parameters (w_1,w_2) building the hierarchic shear 
*   difference vector
* The 5p element takes shear deformations into account according to the Reissner-Mindlin shell theory.
* The hierarchy means that the shell element builds on top of the Kirchhoff-Love shell element 
*   (iga_shell_3p_element) in the sense that the existing equations remain unchanged and additional terms are 
*   added.
* 
* Implementation:
* So far there are no new parameters defined, instead ROTATION_X = w_1 and ROTATION_Y = w_2.
*
* Attention:
* In the last version of this element, there were still open questions and problems. For further information please do not hesitate to contact the author.
*/


// System includes

// External includes

// Project includes
#include "custom_elements/shell_5p_hierarchic_element.h"

namespace Kratos
{
    void Shell5pHierarchicElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        //Constitutive Law initialisation
        InitializeMaterial();

        CalculateMetric(mInitialMetric);

        mZeta = 0.0;
        mInitialTransConToCar = ZeroMatrix(5, 5);
    }

    void Shell5pHierarchicElement::InitializeMaterial()
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
            mConstitutiveLawVector[point_number]->InitializeMaterial(
                r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void Shell5pHierarchicElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        // definition of problem size
        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size = number_of_control_points * 5;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        // initialization of shear difference vector
        array_1d<double, 3> w(3, 0.0);  // shear difference vector
        array_1d<double, 3> Dw_D1(3, 0.0), Dw_D2(3, 0.0);   // derivatives of the shear difference vector
        array_1d<double, 2> w_alpha(2, 0.0);    // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * a1 + w_alpha(2) * a2)
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);   // derivatives of the components w_alpha (lines = component w_alpha; columns = coorindates for which the component is derived)

        // computation of actual metric and shear difference vector
        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);

        array_1d<double, 3> G1, G2, G1xG2;
        double thickness = GetProperties().GetValue(THICKNESS);

        // Gauss Integration over thickness
        for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussIntegrationThickness.zeta(Gauss_index);
            
            // Differential Volume
            array_1d<double, 3> G1(3, 0.0), G2(3, 0.0), G1_con(3, 0.0), G2_con(3, 0.0);
            CalculateInitialBaseVectorsLinearised(G1, G2, G1_con, G2_con);
            double dV = inner_prod(MathUtils<double>::CrossProduct(G1, G2), mInitialMetric.a3_kirchhoff_love);   // mInitialMetric.a3_kirchhoff_love = G3

            // Transformation Matrix
            CalculateTransformationMatrixInitialTransConToCar(G1_con, G2_con);

            // calculate strains and stresses
            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, w, Dw_D1, Dw_D2, constitutive_variables, Values, 
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate variations (second variations according to Kirchhoff-Love are later computed, if needed)
            Matrix B = ZeroMatrix(5, mat_size);
            SecondVariations second_variations(mat_size);
            CalculateB(B, actual_metric); // only Kirchhoff-Love contribution
            CalculateVariationsReissnerMindlin(B, second_variations, w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, 
                actual_metric, CalculateStiffnessMatrixFlag);   // only Reissner-Mindlin contribution, the contributions in B are added to that one of the Kirchhoff-Love part

            double integration_weight = mGaussIntegrationThickness.integration_weight_thickness(Gauss_index) * 
                GetGeometry().IntegrationPoints()[0].Weight()  * dV * thickness / 2.0;

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                // Second variations (Kirchhoff-Love)
                CalculateSecondVariations(second_variations, actual_metric);
                
                // adding linear contributions to the stiffness matrix
                CalculateAndAddKm(rLeftHandSideMatrix, B, constitutive_variables.D, integration_weight);
                // adding  non-linear contribution to stiffness matrix
                CalculateAndAddNonlinearKm(rLeftHandSideMatrix, second_variations, constitutive_variables.S, integration_weight);
            }

            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), constitutive_variables.S);
            }
        } // loop GP_thickness

        KRATOS_CATCH("");
    }

    void Shell5pHierarchicElement::CalculateAndAddKm(
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

    void Shell5pHierarchicElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY

        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size = number_of_control_points * 5;

        KRATOS_DEBUG_ERROR_IF(SD.size() != 5) << "Wrong size of stress vector." << std::endl;

        for (IndexType n = 0; n < mat_size; n++)
        {
            for (IndexType m = 0; m <= n; m++)
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

    void Shell5pHierarchicElement::CalculateMetric(
        MetricVariables& rMetric)
    {
        const IndexType IntegrationPointIndex = 0;
        const GeometryType& r_geometry = GetGeometry();

        r_geometry.Jacobian(rMetric.J, IntegrationPointIndex);

        rMetric.a1[0] = rMetric.J(0, 0);
        rMetric.a2[0] = rMetric.J(0, 1);
        rMetric.a1[1] = rMetric.J(1, 0);
        rMetric.a2[1] = rMetric.J(1, 1);
        rMetric.a1[2] = rMetric.J(2, 0);
        rMetric.a2[2] = rMetric.J(2, 1);

        // base vector a3_kirchhoff_love_tilde
        MathUtils<double>::CrossProduct(rMetric.a3_kirchhoff_love_tilde, rMetric.a1, rMetric.a2);
        // differential area dA
        rMetric.dA = norm_2(rMetric.a3_kirchhoff_love_tilde);
        // normalized base vector a3_kirchhoff_love
        rMetric.a3_kirchhoff_love = rMetric.a3_kirchhoff_love_tilde / rMetric.dA;
        
        // covariant metric
        rMetric.a_ab[0] = pow(rMetric.a1[0], 2) + pow(rMetric.a1[1], 2) + pow(rMetric.a1[2], 2);
        rMetric.a_ab[1] = pow(rMetric.a2[0], 2) + pow(rMetric.a2[1], 2) + pow(rMetric.a2[2], 2);
        rMetric.a_ab[2] = rMetric.a1[0] * rMetric.a2[0] + rMetric.a1[1] * rMetric.a2[1] + rMetric.a1[2] * rMetric.a2[2];

        CalculateHessian(
            rMetric.H, GetGeometry().ShapeFunctionDerivatives(
                2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // base vector derivatives
        for (IndexType i = 0; i < 3; i++)
        {
            rMetric.Da1_D1[i] = rMetric.H(i, 0);
            rMetric.Da2_D2[i] = rMetric.H(i, 1);
            rMetric.Da1_D2[i] = rMetric.H(i, 2);
        }

        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.a3_kirchhoff_love[0] + rMetric.H(1, 0) * rMetric.a3_kirchhoff_love[1] + rMetric.H(2, 0) * rMetric.a3_kirchhoff_love[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.a3_kirchhoff_love[0] + rMetric.H(1, 1) * rMetric.a3_kirchhoff_love[1] + rMetric.H(2, 1) * rMetric.a3_kirchhoff_love[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.a3_kirchhoff_love[0] + rMetric.H(1, 2) * rMetric.a3_kirchhoff_love[1] + rMetric.H(2, 2) * rMetric.a3_kirchhoff_love[2];

        // contravariant metric a_ab_con
        double inverse_determinant_a_ab = 1.0 / (rMetric.a_ab[0] * rMetric.a_ab[1] - rMetric.a_ab[2] * rMetric.a_ab[2]);
        rMetric.a_ab_con[0] = inverse_determinant_a_ab*rMetric.a_ab[1];
        rMetric.a_ab_con[2] = -inverse_determinant_a_ab*rMetric.a_ab[2];
        rMetric.a_ab_con[1] = inverse_determinant_a_ab*rMetric.a_ab[0];
        // contravariant base vectors a_con
        rMetric.a1_con = rMetric.a1*rMetric.a_ab_con[0] + rMetric.a2*rMetric.a_ab_con[2];
        rMetric.a2_con = rMetric.a1*rMetric.a_ab_con[2] + rMetric.a2*rMetric.a_ab_con[1];
        // a3_con = a3_kirchhoff_love
    }

    void Shell5pHierarchicElement::CalculateTransformationMatrixInitialTransConToCar(
        const array_1d<double, 3>& rG1_con, 
        const array_1d<double, 3>& rG2_con)
    {
        // initial local cartesian coordinates
        double lg1 = norm_2(mInitialMetric.a1);
        array_1d<double, 3> E1 = mInitialMetric.a1 / lg1;
        double lg_con2 = norm_2(mInitialMetric.a2_con);
        array_1d<double, 3> E2 = mInitialMetric.a2_con / lg_con2;
        // e3 = a3_kirchhoff_love = a3_con

        double G_00 = inner_prod(E1, rG1_con);
        double G_10 = inner_prod(E2, rG1_con);
        double G_11 = inner_prod(E2, rG2_con);
        
        // some entries remain zero due to the particular definition of the local Cartesian coordinate system
        mInitialTransConToCar = ZeroMatrix(5, 5);
        mInitialTransConToCar(0, 0) = pow(G_00, 2);
        mInitialTransConToCar(1, 0) = pow(G_10, 2);
        mInitialTransConToCar(1, 1) = pow(G_11, 2);
        mInitialTransConToCar(1, 2) = 2.00 * G_10 * G_11;
        mInitialTransConToCar(2, 0) = 2.00 * G_00 * G_10;
        mInitialTransConToCar(2, 2) = 2.00 * G_00 * G_11;
        mInitialTransConToCar(3, 3) = 2.00 * G_11;
        mInitialTransConToCar(3, 4) = 2.00 * G_10;
        mInitialTransConToCar(4, 4) = 2.00 * G_00;
    }

    void Shell5pHierarchicElement::CalculateTransformationMatrixInitialTransCarToCov(Matrix& rInitialTransCarToCov)
    {
        rInitialTransCarToCov = trans(mInitialTransConToCar);
        // division by 2.0 because not used for strains but for stresses (strains have e.g. entries with 2e12)
        rInitialTransCarToCov(2, 0) = rInitialTransCarToCov(2, 0) / 2.0;
        rInitialTransCarToCov(2, 1) = rInitialTransCarToCov(2, 1) / 2.0;
        rInitialTransCarToCov(2, 2) = rInitialTransCarToCov(2, 2) / 2.0;
        rInitialTransCarToCov(3, 3) = rInitialTransCarToCov(3, 3) / 2.0;
        rInitialTransCarToCov(4, 3) = rInitialTransCarToCov(4, 3) / 2.0;
        rInitialTransCarToCov(4, 4) = rInitialTransCarToCov(4, 4) / 2.0;
    }

    void Shell5pHierarchicElement::CalculateTransformationMatrixActualTransCovToCar(
        Matrix& rActualTransCovToCar,
        const Vector& rg1,
        const Vector& rg2,
        const Vector& rg3,
        const Vector& ra2_con,
        const Vector& ra3_kirchhoff_love)
    {
        // local cartesian coordinates, direction of a2_con and g2_con is the same -> it is possible to use ra2_con instead of g2_con 
        array_1d<double, 3> e1 = rg1 / norm_2(rg1);
        array_1d<double, 3> e2 = ra2_con / norm_2(ra2_con);
        // e3 = a3_kirchhoff_love

        double G_00 = inner_prod(e1, rg1);
        double G_01 = inner_prod(e1, rg2);
        double G_02 = inner_prod(e1, rg3);
        double G_11 = inner_prod(e2, rg2);
        double G_12 = inner_prod(e2, rg3);
        double G_22 = inner_prod(ra3_kirchhoff_love, rg3);
 
        rActualTransCovToCar(0, 0) = pow(G_00, 2);
        rActualTransCovToCar(0, 1) = pow(G_01, 2);
        rActualTransCovToCar(0, 2) = 2 * G_00 * G_01;
        rActualTransCovToCar(0, 3) = 2 * G_01 * G_02;
        rActualTransCovToCar(0, 4) = 2 * G_00 * G_02;
        rActualTransCovToCar(1, 1) = pow(G_11, 2);
        rActualTransCovToCar(1, 3) = 2 * G_11 * G_12;
        rActualTransCovToCar(2, 1) = G_01 * G_11;
        rActualTransCovToCar(2, 2) = G_00 * G_11;
        rActualTransCovToCar(2, 3) = G_01 * G_12 + G_02 * G_11;
        rActualTransCovToCar(2, 4) = G_00 * G_12;
        rActualTransCovToCar(3, 3) = G_11;
        rActualTransCovToCar(4, 3) = G_01 * G_22;
        rActualTransCovToCar(4, 4) = G_00 * G_22;
    }

    void Shell5pHierarchicElement::CalculateShearDifferenceVector(
        array_1d<double, 3>& rw,
        array_1d<double, 3>& rDw_D1,
        array_1d<double, 3>& rDw_D2,
        array_1d<double, 2>& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex)
    {
        const GeometryType& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const SizeType number_of_control_points = r_geometry.size();

        for (IndexType i = 0; i < number_of_control_points; ++i) 
        {
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector)
            double w_1 = r_geometry[i].GetDof(ROTATION_X).GetSolutionStepValue();
            double w_2 = r_geometry[i].GetDof(ROTATION_Y).GetSolutionStepValue();
            rDw_alpha_Dbeta(0, 0) += r_DN_De(i, 0) * w_1;
            rDw_alpha_Dbeta(0, 1) += r_DN_De(i, 1) * w_1;
            rDw_alpha_Dbeta(1, 0) += r_DN_De(i, 0) * w_2;
            rDw_alpha_Dbeta(1, 1) += r_DN_De(i, 1) * w_2;
            rw_alpha(0) += r_N(IntegrationPointIndex, i) * w_1;
            rw_alpha(1) += r_N(IntegrationPointIndex, i) * w_2;
        }

        // derivatives of the shear difference vector
        rDw_D1 = rDw_alpha_Dbeta(0, 0) * rActualMetric.a1 + rDw_alpha_Dbeta(1, 0) * rActualMetric.a2;
        rDw_D2 = rDw_alpha_Dbeta(0, 1) * rActualMetric.a1 + rDw_alpha_Dbeta(1, 1) * rActualMetric.a2;

        for (IndexType i = 0; i < 3; i++)
        {
            rDw_D1[i] += rw_alpha(0) * rActualMetric.H(i, 0) + rw_alpha(1) * rActualMetric.H(i, 2);
            rDw_D2[i] += rw_alpha(0) * rActualMetric.H(i, 2) + rw_alpha(1) * rActualMetric.H(i, 1);
        }

        rw = rw_alpha(0) * rActualMetric.a1 + rw_alpha(1) * rActualMetric.a2;
    }
    
    void Shell5pHierarchicElement::CalculateInitialBaseVectorsLinearised(
        array_1d<double, 3>& rG1,
        array_1d<double, 3>& rG2,
        array_1d<double, 3>& rG1_con,
        array_1d<double, 3>& rG2_con)
    {
        double thickness = GetProperties()[THICKNESS];
        
        array_1d<double, 3> DA3_D1(3, 0.0), DA3_D2(3, 0.0), DA1_D1xA2(3, 0.0), A1xDA2_D1(3, 0.0), DA1_D2xA2(3, 0.0), A1xDA2_D2(3, 0.0), 
            G1xG2(3, 0.0), G_ab(3, 0.0), G_ab_con(3, 0.0);

        MathUtils<double>::CrossProduct(DA1_D1xA2, mInitialMetric.Da1_D1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D1, mInitialMetric.a1, mInitialMetric.Da1_D2); // DA1_D2 = DA2_D1
        MathUtils<double>::CrossProduct(DA1_D2xA2, mInitialMetric.Da1_D2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xDA2_D2, mInitialMetric.a1, mInitialMetric.Da2_D2);
        DA3_D1 = ((DA1_D1xA2 + A1xDA2_D1) * mInitialMetric.dA - mInitialMetric.a3_kirchhoff_love_tilde * norm_2(DA1_D1xA2 + A1xDA2_D1)) 
            / (mInitialMetric.dA * mInitialMetric.dA);
        DA3_D2 = ((DA1_D2xA2 + A1xDA2_D2) * mInitialMetric.dA - mInitialMetric.a3_kirchhoff_love_tilde * norm_2(DA1_D2xA2 + A1xDA2_D2))
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

        // contravariant base vectors of the shell body in the reference configuration
        rG1_con = rG1 * G_ab_con[0] + rG2 * G_ab_con[2];
        rG2_con = rG1 * G_ab_con[2] + rG2 * G_ab_con[1];
        // G3_con = A3
    }

    void Shell5pHierarchicElement::CalculateActualBaseVectorsLinearised(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        array_1d<double, 3>& rg1,
        array_1d<double, 3>& rg2,
        array_1d<double, 3>& rg3
        )
    {
        double thickness = GetProperties()[THICKNESS];

        double DdA_D1 = 0.0;
        double DdA_D2 = 0.0;
        array_1d<double, 3> Da1xa2_D1(3, 0.0), Da1xa2_D2(3, 0.0);
        array_1d<double, 3> Da3_KL_D1(3, 0.0), Da3_KL_D2(3, 0.0);

        // derivatives of a3_kirchhoff_love w.r.t. theta_alpha
        // derivative of a3_kirchhoff_love w.r.t. theta1: (a1xa2)'= a1'xa2 + a1xa2'
        Da1xa2_D1 = MathUtils<double>::CrossProduct(rActualMetric.Da1_D1, rActualMetric.a2) + MathUtils<double>::CrossProduct(rActualMetric.a1, rActualMetric.Da1_D2);

        // derivative of a3_kirchhoff_love w.r.t. theta2: (a1xa2)'= a1'xa2 + a1xa2'
        array_1d<double, 3> Da1_D2xa2, a1xDa2_D2;
        MathUtils<double>::CrossProduct(Da1_D2xa2, rActualMetric.Da1_D2, rActualMetric.a2);
        MathUtils<double>::CrossProduct(a1xDa2_D2, rActualMetric.a1, rActualMetric.Da2_D2);
        Da1xa2_D2 = Da1_D2xa2 + a1xDa2_D2;

        /* Ableitung des Nenners von a3_kirchhoff_love nach 1 */
        DdA_D1 = (Da1xa2_D1[0] * rActualMetric.a3_kirchhoff_love_tilde[0] 
            + Da1xa2_D1[1] * rActualMetric.a3_kirchhoff_love_tilde[1] + Da1xa2_D1[2] * rActualMetric.a3_kirchhoff_love_tilde[2]) / rActualMetric.dA;

        /* Ableitung des Nenners von a3_kirchhoff_love nach 2 */
        DdA_D2 = (Da1xa2_D2[0] * rActualMetric.a3_kirchhoff_love_tilde[0] 
            + Da1xa2_D2[1] * rActualMetric.a3_kirchhoff_love_tilde[1] + Da1xa2_D2[2] * rActualMetric.a3_kirchhoff_love_tilde[2]) / rActualMetric.dA;

        /* Ableitung von a3_kirchhoff_love mit Quotientenregel */
        for (unsigned int i = 0; i < 3; i++)
        {
            Da3_KL_D1[i] = (Da1xa2_D1[i] * rActualMetric.dA - rActualMetric.a3_kirchhoff_love_tilde[i] * DdA_D1) / 
                (rActualMetric.dA * rActualMetric.dA);
            Da3_KL_D2[i] = (Da1xa2_D2[i] * rActualMetric.dA - rActualMetric.a3_kirchhoff_love_tilde[i] * DdA_D2) / 
                (rActualMetric.dA * rActualMetric.dA);
        }

        /* Kovariante Basisvektoren */
        rg1 = rActualMetric.a1 + thickness / 2.0 * mZeta*(Da3_KL_D1 + rDw_D1);
        rg2 = rActualMetric.a2 + thickness / 2.0 * mZeta*(Da3_KL_D2 + rDw_D2);
        rg3 = rActualMetric.a3_kirchhoff_love + rw;     // g3 = a3
    }

    void Shell5pHierarchicElement::CalculateDeformationGradient(
        const array_1d<double, 3> rG1,
        const array_1d<double, 3> rG2,
        const array_1d<double, 3> rg1,
        const array_1d<double, 3> rg2,
        const array_1d<double, 3> rg3,
        Matrix& rF,
        double& rdetF)
    {
        // Initialization
        rF.resize(3, 3);
        rF = ZeroMatrix(3, 3);
        rdetF = 0.0;

        // Covariant metric
        array_1d<double, 3> G_ab(3, 0.0);
        G_ab[0] = pow(mInitialMetric.a1[0], 2) + pow(mInitialMetric.a1[1], 2) + pow(mInitialMetric.a1[2], 2);
        G_ab[1] = pow(mInitialMetric.a2[0], 2) + pow(mInitialMetric.a2[1], 2) + pow(mInitialMetric.a2[2], 2);
        G_ab[2] = mInitialMetric.a1[0] * mInitialMetric.a2[0] + mInitialMetric.a1[1] * mInitialMetric.a2[1] 
            + mInitialMetric.a1[2] * mInitialMetric.a2[2];
        // Contravariant metric
        array_1d<double, 3> G_ab_con(3, 0.0);
        double inv_detG_ab = 1.0 / (G_ab[0] * G_ab[1] - G_ab[2] * G_ab[2]);
        G_ab_con[0] = inv_detG_ab * G_ab[1];
        G_ab_con[2] = -inv_detG_ab * G_ab[2];
        G_ab_con[1] = inv_detG_ab * G_ab[0];
        // Contravariant base vectors
        array_1d<double, 3> G1_con(3, 0.0), G2_con(3, 0.0);
        G1_con = rG1 * G_ab_con[0] + rG2 * G_ab_con[2];
        G2_con = rG1 * G_ab_con[2] + rG2 * G_ab_con[1];
    
        // deformation gradient and its determinant
        rF = outer_prod(rg1, G1_con) + outer_prod(rg2, G2_con) + outer_prod(rg3, mInitialMetric.a3_kirchhoff_love);
        rdetF = rF(0, 0) * (rF(1, 1) * rF(2, 2) - rF(2, 1) * rF(1, 2)) 
            - rF(1, 0) * (rF(0, 1) * rF(2, 2) - rF(2, 1) * rF(0, 2))
            + rF(2, 0) * (rF(0, 1) * rF(1, 2) - rF(1, 1) * rF(0, 2));
    }

    void Shell5pHierarchicElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        KRATOS_TRY

        array_1d<double, 5> strain_vector(5, 0.0), strain_vector_reissner_mindlin(5, 0.0);
        
        // Strain computation in curvilinear space
        CalculateStrain(strain_vector, rActualMetric.a_ab, rActualMetric.curvature);
        CalculateStrainReissnerMindlin(strain_vector_reissner_mindlin, rw, rDw_D1, rDw_D2, rActualMetric.a1, rActualMetric.a2);
        rThisConstitutiveVariables.E = strain_vector + strain_vector_reissner_mindlin;

        // Strain transformation to local Cartesian Space with VoigtSize 6 because ConstitutiveLaw is 3D
        ConstitutiveVariables constitutive_variables(6);
        TransformationCurvilinearStrainSize5ToCartesianStrainSize6(rThisConstitutiveVariables.E, constitutive_variables.E);

        //Constitutive Matrix D
        rValues.SetStrainVector(constitutive_variables.E); //this is the input parameter
        rValues.SetStressVector(constitutive_variables.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(constitutive_variables.D); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);
        // static condensation of  sigma_33
        IndexType index_i = 0;
        for (IndexType i = 0; i < 6; i++)
        {
            if (i != 2){
                IndexType index_j = 0;
                for (IndexType j = 0; j < 6; j++ ){
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
        rThisConstitutiveVariables.E = prod(mInitialTransConToCar, rThisConstitutiveVariables.E);

        //Local Cartesian Stresses
        rThisConstitutiveVariables.S = prod(
            trans(rThisConstitutiveVariables.D), rThisConstitutiveVariables.E);

        KRATOS_CATCH("")
    }

    void Shell5pHierarchicElement::CalculateStrain(
        array_1d<double, 5>& rStrainVector,
        const Vector& ra_ab,
        const Vector& rCurvature)
    {
        KRATOS_TRY
        
        double thickness = GetProperties().GetValue(THICKNESS);

        rStrainVector[0] = 0.5 * (ra_ab[0] - mInitialMetric.a_ab[0]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[0] - rCurvature[0]);
        rStrainVector[1] = 0.5 * (ra_ab[1] - mInitialMetric.a_ab[1]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[1] - rCurvature[1]);
        rStrainVector[2] = 0.5 * (ra_ab[2] - mInitialMetric.a_ab[2]) + mZeta * thickness / 2.0 * (mInitialMetric.curvature[2] - rCurvature[2]);
        // the other entries are (remain) zero (Kirchhoff-Love)
        
        KRATOS_CATCH("")
    }

    void Shell5pHierarchicElement::CalculateStrainReissnerMindlin(
        array_1d<double, 5>& rStrainVectorReissnerMindlin,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& ra1,
        const Vector& ra2)
    {
        double thickness = GetProperties().GetValue(THICKNESS);
                
        rStrainVectorReissnerMindlin[0] = mZeta * thickness/2.0 * inner_prod(rDw_D1, ra1);
        rStrainVectorReissnerMindlin[1] = mZeta * thickness/2.0 * inner_prod(rDw_D2, ra2);
        rStrainVectorReissnerMindlin[2] = mZeta * thickness/2.0 * 0.5 * (inner_prod(rDw_D1, ra2) + inner_prod(rDw_D2, ra1));
        rStrainVectorReissnerMindlin[3] = 0.5 * inner_prod(rw, ra2);
        rStrainVectorReissnerMindlin[4] = 0.5 * inner_prod(rw, ra1);
    }

    void Shell5pHierarchicElement::TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector& rCurvilinearStrain,
        Vector& rCartesianStrain)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rCurvilinearStrain.size() != 5 || rCartesianStrain.size() != 6) << "Wrong strain size in transformation." << std::endl;
        KRATOS_DEBUG_ERROR_IF(mInitialTransConToCar.size1() != 5 || mInitialTransConToCar.size2() != 5) << "Wrong size of transformation matrix mInitialTransConToCar." << std::endl;

        // transformation with simplified matrix
        rCartesianStrain[0] = mInitialTransConToCar(0, 0) * rCurvilinearStrain[0];
        rCartesianStrain[1] = mInitialTransConToCar(1, 0) * rCurvilinearStrain[0] + mInitialTransConToCar(1, 1) * rCurvilinearStrain[1] 
            + mInitialTransConToCar(1, 2) * rCurvilinearStrain[2];
        rCartesianStrain[2] = 0.0; // RM
        rCartesianStrain[3] = mInitialTransConToCar(2, 0) * rCurvilinearStrain[0] + mInitialTransConToCar(2, 2) * rCurvilinearStrain[2];
        rCartesianStrain[4] = mInitialTransConToCar(3, 3) * rCurvilinearStrain[3] + mInitialTransConToCar(3, 4) * rCurvilinearStrain[4];
        rCartesianStrain[5] = mInitialTransConToCar(4, 4) * rCurvilinearStrain[4];

        KRATOS_CATCH("")
    }

    void Shell5pHierarchicElement::CalculateB(
        Matrix& rB,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex)
    {
        const GeometryType& r_geometry = GetGeometry();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);
        const double thickness = GetProperties().GetValue(THICKNESS);

        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size_KL = number_of_control_points * 3;
        SizeType mat_size = number_of_control_points * 5;

        // membrane part
        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            int kr = r / 5;
            int dirr = r % 5;

            array_1d<double, 3> dE_curvilinear(3, 0.0);
            // "if" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            if (dirr == 0 || dirr == 1 || dirr == 2)
            {
                // strain corresponding to E11, E22, E12
                dE_curvilinear[0] = r_DN_De(kr, 0)*rActualMetric.a1(dirr);
                dE_curvilinear[1] = r_DN_De(kr, 1)*rActualMetric.a2(dirr);
                dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualMetric.a2(dirr) + rActualMetric.a1(dirr)* r_DN_De(kr, 1));
            }
            // calculated with simplified mInitialTransConToCar
            rB(0, r) += mInitialTransConToCar(0, 0) * dE_curvilinear[0] + mInitialTransConToCar(0, 1) * dE_curvilinear[1] 
                + mInitialTransConToCar(0, 2) * dE_curvilinear[2] ;
            rB(1, r) += mInitialTransConToCar(1, 0) * dE_curvilinear[0] + mInitialTransConToCar(1, 1) * dE_curvilinear[1] 
                + mInitialTransConToCar(1, 2) * dE_curvilinear[2];
            rB(2, r) += mInitialTransConToCar(2, 0) * dE_curvilinear[0] + mInitialTransConToCar(2, 1) * dE_curvilinear[1] 
                + mInitialTransConToCar(2, 2) * dE_curvilinear[2];

            // all other entries of rB are (remain) zero
            
        }

        // curvature part
        Matrix dg3 = ZeroMatrix(3, 3);
        Matrix dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size_KL);

        double invdA = 1 / rActualMetric.dA;
        double inddA3 = 1 / std::pow(rActualMetric.dA, 3);

        for (unsigned int i = 0; i < number_of_control_points; i++)
        {
            unsigned int index_KL = 3 * i;
            unsigned int index = 5 * i;
            //first line
            dg3(0, 0) = 0;
            dg3(0, 1) = -r_DN_De(i, 0) * rActualMetric.a2[2] + r_DN_De(i, 1)*rActualMetric.a1[2];
            dg3(0, 2) = r_DN_De(i, 0) * rActualMetric.a2[1] - r_DN_De(i, 1)*rActualMetric.a1[1];

            //second line
            dg3(1, 0) = r_DN_De(i, 0) * rActualMetric.a2[2] - r_DN_De(i, 1)*rActualMetric.a1[2];
            dg3(1, 1) = 0;
            dg3(1, 2) = -r_DN_De(i, 0)*rActualMetric.a2[0] + r_DN_De(i, 1)*rActualMetric.a1[0];

            //third line
            dg3(2, 0) = -r_DN_De(i, 0) * rActualMetric.a2[1] + r_DN_De(i, 1) * rActualMetric.a1[1];
            dg3(2, 1) = r_DN_De(i, 0) * rActualMetric.a2[0] - r_DN_De(i, 1) * rActualMetric.a1[0];
            dg3(2, 2) = 0;

            for (unsigned int j = 0; j < 3; j++)
            {
                double g3dg3lg3 = (rActualMetric.a3_kirchhoff_love_tilde[0] * dg3(j, 0) + rActualMetric.a3_kirchhoff_love_tilde[1] * dg3(j, 1) + rActualMetric.a3_kirchhoff_love_tilde[2] * dg3(j, 2))*inddA3;

                dn(j, 0) = dg3(j, 0)*invdA - rActualMetric.a3_kirchhoff_love_tilde[0] * g3dg3lg3;
                dn(j, 1) = dg3(j, 1)*invdA - rActualMetric.a3_kirchhoff_love_tilde[1] * g3dg3lg3;
                dn(j, 2) = dg3(j, 2)*invdA - rActualMetric.a3_kirchhoff_love_tilde[2] * g3dg3lg3;
            }

            // b refers to curvilinear and rB to local cartesian coordinate system
            // "index" guarantees that there are zero entries corresponding to the new parameters w_1 and w_2
            for (unsigned int j = 0; j < 3; j++)
            {
                b(0, index_KL + j) = -mZeta * thickness / 2.0 * (r_DDN_DDe(i, 0) * rActualMetric.a3_kirchhoff_love[j] + rActualMetric.H(0, 0) * dn(j, 0)
                    + rActualMetric.H(1, 0) * dn(j, 1) + rActualMetric.H(2, 0) * dn(j, 2));
                b(1, index_KL + j) = -mZeta * thickness / 2.0 * (r_DDN_DDe(i, 2) * rActualMetric.a3_kirchhoff_love[j] + rActualMetric.H(0, 1) * dn(j, 0)
                    + rActualMetric.H(1, 1) * dn(j, 1) + rActualMetric.H(2, 1) * dn(j, 2));
                b(2, index_KL + j) = -mZeta * thickness / 2.0 * (r_DDN_DDe(i, 1) * rActualMetric.a3_kirchhoff_love[j] + rActualMetric.H(0, 2) * dn(j, 0)
                    + rActualMetric.H(1, 2) * dn(j, 1) + rActualMetric.H(2, 2) * dn(j, 2));

                rB(0, index + j) += mInitialTransConToCar(0, 0) * b(0, index_KL + j);
                rB(1, index + j) += mInitialTransConToCar(1, 0) * b(0, index_KL + j) + mInitialTransConToCar(1, 1) * b(1, index_KL + j)
                    + mInitialTransConToCar(1, 2) * b(2, index_KL + j);
                rB(2, index + j) += mInitialTransConToCar(2, 0) * b(0, index_KL + j) + mInitialTransConToCar(2, 2) * b(2, index_KL + j);
            }
        }
    }

    void Shell5pHierarchicElement::CalculateSecondVariations(
        SecondVariations& rSecondVariations,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex)
    {
        const GeometryType& r_geometry = GetGeometry();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);
        const double thickness = GetProperties().GetValue(THICKNESS);

        SizeType number_of_control_points = r_geometry.size();
        SizeType mat_size_KL = number_of_control_points * 3;
        SizeType mat_size = number_of_control_points * 5;
           
        double lg3_3 = pow(rActualMetric.dA, 3);
        double lg3_5 = pow(rActualMetric.dA, 5);
        double inv_lg3 = 1 / rActualMetric.dA;
        double inv_lg3_3 = 1 / lg3_3;
        double inv_lg3_5 = 1 / lg3_5;

        SecondVariations second_variations_KL(mat_size_KL);
        Matrix S_dg3 = ZeroMatrix(3, mat_size_KL);
        Vector S_g3dg3 = ZeroVector(mat_size_KL);
        Vector S_g3dg3lg3_3 = ZeroVector(mat_size_KL);
        Matrix S_dn = ZeroMatrix(3, mat_size_KL);
        // first variation of strain and curvature w.r.t. dof
        for (IndexType r = 0; r < mat_size_KL; r++)
        {
            // local node number kr and dof direction dirr
            unsigned int kr = r / 3;
            unsigned int dirr = r % 3;

            array_1d<double, 3> S_dg_1(3, 0.0), S_dg_2(3, 0.0);
            S_dg_1(dirr) = r_DN_De(kr, 0);
            S_dg_2(dirr) = r_DN_De(kr, 1);

            // curvature
            S_dg3(0, r) = S_dg_1(1)*rActualMetric.a2(2) - S_dg_1(2)*rActualMetric.a2(1) + rActualMetric.a1(1)*S_dg_2(2) - rActualMetric.a1(2)*S_dg_2(1);
            S_dg3(1, r) = S_dg_1(2)*rActualMetric.a2(0) - S_dg_1(0)*rActualMetric.a2(2) + rActualMetric.a1(2)*S_dg_2(0) - rActualMetric.a1(0)*S_dg_2(2);
            S_dg3(2, r) = S_dg_1(0)*rActualMetric.a2(1) - S_dg_1(1)*rActualMetric.a2(0) + rActualMetric.a1(0)*S_dg_2(1) - rActualMetric.a1(1)*S_dg_2(0);

            S_g3dg3[r] = rActualMetric.a3_kirchhoff_love_tilde[0] * S_dg3(0, r) + rActualMetric.a3_kirchhoff_love_tilde[1] * S_dg3(1, r) + rActualMetric.a3_kirchhoff_love_tilde[2] * S_dg3(2, r);
            S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

            S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rActualMetric.a3_kirchhoff_love_tilde[0] * S_g3dg3lg3_3[r];
            S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rActualMetric.a3_kirchhoff_love_tilde[1] * S_g3dg3lg3_3[r];
            S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rActualMetric.a3_kirchhoff_love_tilde[2] * S_g3dg3lg3_3[r];
        }

        // second variation of strain and curvature w.r.t. dofs
        for (IndexType r = 0; r < mat_size_KL; r++)
        {
            // local node number kr and dof direction dirr
            unsigned int kr = r / 3;
            unsigned int dirr = r % 3;

            for (IndexType s = 0; s <= r; s++)
            {
                // local node number ks and dof direction dirs
                unsigned int ks = s / 3;
                unsigned int dirs = s % 3;

                // strain
                array_1d<double, 3> ddE_cu(3, 0.0);
                if (dirr == dirs)
                {
                    ddE_cu[0] = r_DN_De(kr, 0)* r_DN_De(ks, 0);
                    ddE_cu[1] = r_DN_De(kr, 1)* r_DN_De(ks, 1);
                    ddE_cu[2] = 0.5*(r_DN_De(kr, 0)* r_DN_De(ks, 1) + r_DN_De(kr, 1)* r_DN_De(ks, 0));

                    // calculated with reduced mInitialTransConToCar
                    second_variations_KL.B11(r, s) += mInitialTransConToCar(0, 0) * ddE_cu[0];
                    second_variations_KL.B22(r, s) += mInitialTransConToCar(1, 0) * ddE_cu[0] + mInitialTransConToCar(1, 1) * ddE_cu[1]
                        + mInitialTransConToCar(1, 2) * ddE_cu[2];
                    second_variations_KL.B12(r, s) += mInitialTransConToCar(2, 0) * ddE_cu[0] + mInitialTransConToCar(2, 2) * ddE_cu[2];
                }

                // curvature
                array_1d<double, 3> ddg3(3, 0.0);
                int dirt = 4 - dirr - dirs;
                int ddir = dirr - dirs;
                if (ddir == -1)      ddg3(dirt - 1) =  r_DN_De(kr, 0) * r_DN_De(ks, 1) - r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == 2)  ddg3(dirt - 1) =  r_DN_De(kr, 0) * r_DN_De(ks, 1) - r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == 1)  ddg3(dirt - 1) = -r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(ks, 0) * r_DN_De(kr, 1);
                else if (ddir == -2) ddg3(dirt - 1) = -r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_DN_De(ks, 0) * r_DN_De(kr, 1);

                double c = -(ddg3[0] * rActualMetric.a3_kirchhoff_love_tilde[0] + ddg3[1] * rActualMetric.a3_kirchhoff_love_tilde[1] + ddg3[2] * rActualMetric.a3_kirchhoff_love_tilde[2]
                    + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                    )*inv_lg3_3;

                double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                array_1d<double, 3> ddn(3, 0.0);
                ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rActualMetric.a3_kirchhoff_love_tilde[0];
                ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rActualMetric.a3_kirchhoff_love_tilde[1];
                ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rActualMetric.a3_kirchhoff_love_tilde[2];

                array_1d<double, 3> ddK_cu(3, 0.0);
                ddK_cu[0] = - mZeta * thickness / 2.0 * (r_DDN_DDe(kr, 0)*S_dn(dirr, s) + r_DDN_DDe(ks, 0)*S_dn(dirs, r)
                    + rActualMetric.H(0, 0)*ddn[0] + rActualMetric.H(1, 0)*ddn[1] + rActualMetric.H(2, 0)*ddn[2]);
                ddK_cu[1] = - mZeta * thickness / 2.0 * (r_DDN_DDe(kr, 2)*S_dn(dirr, s) + r_DDN_DDe(ks, 2)*S_dn(dirs, r)
                    + rActualMetric.H(0, 1)*ddn[0] + rActualMetric.H(1, 1)*ddn[1] + rActualMetric.H(2, 1)*ddn[2]);
                ddK_cu[2] = - mZeta * thickness / 2.0 * (r_DDN_DDe(kr, 1)*S_dn(dirr, s) + r_DDN_DDe(ks, 1)*S_dn(dirs, r)
                    + rActualMetric.H(0, 2)*ddn[0] + rActualMetric.H(1, 2)*ddn[1] + rActualMetric.H(2, 2)*ddn[2]);

                // calculated with reduced mInitialTransConToCar
                second_variations_KL.B11(r, s) += mInitialTransConToCar(0, 0) * ddK_cu[0];
                second_variations_KL.B22(r, s) += mInitialTransConToCar(1, 0) * ddK_cu[0] + mInitialTransConToCar(1, 1) * ddK_cu[1] 
                    + mInitialTransConToCar(1, 2) * ddK_cu[2];
                second_variations_KL.B12(r, s) += mInitialTransConToCar(2, 0) * ddK_cu[0] + mInitialTransConToCar(2, 2) * ddK_cu[2];
                    
                // symmetric B matrices
                if (r != s){
                    second_variations_KL.B11(s, r) = second_variations_KL.B11(r, s);
                    second_variations_KL.B22(s, r) = second_variations_KL.B22(r, s);
                    second_variations_KL.B12(s, r) = second_variations_KL.B12(r, s);
                }
            }
        }

        // transfer Kirchhoff-Love-second-variations to Reissner-Mindlin-second-variations
        for (IndexType r = 0; r < mat_size; r++) {
            unsigned int kr = r / 5;
            unsigned int dirr = r % 5;
            unsigned int r_KL = kr * 3 + dirr;
            if (dirr != 3 && dirr != 4){
                for (unsigned int s = 0; s < mat_size; s++){
                    unsigned int ks = s / 5;
                    unsigned int dirs = s % 5;
                    unsigned int s_KL = ks * 3 + dirs;
                    if (dirs != 3 && dirs != 4){
                        rSecondVariations.B11(r, s) += second_variations_KL.B11(r_KL, s_KL);
                        rSecondVariations.B22(r, s) += second_variations_KL.B22(r_KL, s_KL);
                        rSecondVariations.B12(r, s) += second_variations_KL.B12(r_KL, s_KL);               
                    }
                }
            }
        }
    }

    void Shell5pHierarchicElement::CalculateVariationsReissnerMindlin(
        Matrix& rB,
        SecondVariations& rSecondVariations,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag,
        IndexType IntegrationPointIndex)
    {
        const GeometryType& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int number_of_control_points = r_geometry.size();
        const unsigned int mat_size = number_of_control_points * 5;

        // 1. First strain variation
        Matrix Dw_Dr = ZeroMatrix(3, mat_size);
        
        for (unsigned int r = 0; r < mat_size; r++){
            // local node number kr and dof direction dirr
            int kr = r / 5;
            int dirr = r % 5;
            
            // the two entries E23 and E13 w.r.t. the curvilinear coord. sys. are saved in dE_cur
            array_1d<double, 2> dE_cur = ZeroVector(2);
            
            if (dirr == 0 || dirr == 1 || dirr == 2){
                Dw_Dr(dirr, r) = rw_alpha(0) * r_DN_De(kr, 0) + rw_alpha(1) * r_DN_De(kr, 1);
                dE_cur[0] = 0.5 * (rw(dirr) * r_DN_De(kr, 1));
                dE_cur[1] = 0.5 * (rw(dirr) * r_DN_De(kr, 0));
            }
            else if(dirr == 3){
                Dw_Dr(0, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a1(0);
                Dw_Dr(1, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a1(1);
                Dw_Dr(2, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a1(2);
            }
            else {
                Dw_Dr(0, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a2(0);
                Dw_Dr(1, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a2(1);
                Dw_Dr(2, r) = r_N(IntegrationPointIndex, kr) * rActualMetric.a2(2);
            }
            dE_cur[0] += 0.5 * (Dw_Dr(0, r) * rActualMetric.a2(0) + Dw_Dr(1, r) * rActualMetric.a2(1) + Dw_Dr(2, r) * rActualMetric.a2(2));
            dE_cur[1] += 0.5 * (Dw_Dr(0, r) * rActualMetric.a1(0) + Dw_Dr(1, r) * rActualMetric.a1(1) + Dw_Dr(2, r) * rActualMetric.a1(2));
            
            // calculated with the reduced mInitialTransConToCar
            rB(3, r) += mInitialTransConToCar(3, 3) * dE_cur[0] + mInitialTransConToCar(3, 4) * dE_cur[1];
            rB(4, r) += mInitialTransConToCar(4, 4) * dE_cur[1];
            // the other entries are (remain) zero

            // 2. First curvature variation
            array_1d<double, 3> dK_cu = ZeroVector(3);
            array_1d<double, 3> DDw_DD1r = ZeroVector(3);
            array_1d<double, 3> DDw_DD2r = ZeroVector(3);

            if (dirr == 0 || dirr == 1 || dirr == 2){
                DDw_DD1r[dirr] = rDw_alpha_Dbeta(0, 0) * r_DN_De(kr, 0)  + rw_alpha[0] * r_DDN_DDe(kr, 0)
                    + rDw_alpha_Dbeta(1, 0) * r_DN_De(kr, 1) + rw_alpha[1] * r_DDN_DDe(kr, 1);
                DDw_DD2r[dirr] = rDw_alpha_Dbeta(0, 1) * r_DN_De(kr, 0) + rw_alpha[0] * r_DDN_DDe(kr, 1)
                    + rDw_alpha_Dbeta(1, 1) * r_DN_De(kr, 1) + rw_alpha[1] * r_DDN_DDe(kr, 2);
                dK_cu[0] = rDw_D1[dirr] * r_DN_De(kr, 0);
                dK_cu[1] = rDw_D2[dirr] * r_DN_De(kr, 1);
                dK_cu[2] = 0.5 * (rDw_D1[dirr] * r_DN_De(kr, 1) + rDw_D2[dirr] * r_DN_De(kr, 0));
            }
            else if (dirr == 3){
                DDw_DD1r += r_DN_De(kr, 0) * rActualMetric.a1;
                DDw_DD1r[0] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(0, 0);
                DDw_DD1r[1] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(1, 0);
                DDw_DD1r[2] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(2, 0);
                DDw_DD2r += r_DN_De(kr, 1) * rActualMetric.a1;
                DDw_DD2r[0] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(0, 2);
                DDw_DD2r[1] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(1, 2);
                DDw_DD2r[2] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(2, 2);
            }
            else if (dirr == 4){
                DDw_DD1r += r_DN_De(kr, 0) * rActualMetric.a2;
                DDw_DD1r[0] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(0, 2);
                DDw_DD1r[1] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(1, 2);
                DDw_DD1r[2] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(2, 2);
                DDw_DD2r += r_DN_De(kr, 1) * rActualMetric.a2;
                DDw_DD2r[0] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(0, 1);
                DDw_DD2r[1] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(1, 1);
                DDw_DD2r[2] += r_N(IntegrationPointIndex, kr) * rActualMetric.H(2, 1);
            }
            dK_cu[0] += inner_prod(DDw_DD1r, rActualMetric.a1);
            dK_cu[1] += inner_prod(DDw_DD2r, rActualMetric.a2);
            dK_cu[2] += 0.5 * (inner_prod(DDw_DD1r, rActualMetric.a2) + inner_prod(DDw_DD2r, rActualMetric.a1));

            // calculated with reduced mInitialTransConToCar
            rB(0, r) += mZeta * thickness / 2.0 * mInitialTransConToCar(0, 0) * dK_cu[0];
            rB(1, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(1, 0) * dK_cu[0] + mInitialTransConToCar(1, 1) * dK_cu[1] 
                + mInitialTransConToCar(1, 2) * dK_cu[2]);
            rB(2, r) += mZeta * thickness / 2.0 * (mInitialTransConToCar(2, 0) * dK_cu[0] + mInitialTransConToCar(2, 2) * dK_cu[2]);
            // all other entries are (remain) zero

            // 3. Second Strain Variation
            if (rCalculateStiffnessMatrixFlag == true){
                for (unsigned int s = 0; s <= r; s++)
                {
                    // local node number ks and dof direction dirs
                    int ks = s / 5;
                    int dirs = s % 5;
                    
                    array_1d <double, 3> Dw_Ds = ZeroVector(3);
                    array_1d <double, 3> DDw_DDrs = ZeroVector(3);
                    array_1d <double, 2> ddE_cu = ZeroVector(2);

                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        Dw_Ds(dirs) = rw_alpha(0) * r_DN_De(ks, 0) + rw_alpha(1) * r_DN_De(ks, 1);
                    }
                    else if(dirs == 3){
                        Dw_Ds(0) = r_N(IntegrationPointIndex, ks) * rActualMetric.a1(0);
                        Dw_Ds(1) = r_N(IntegrationPointIndex, ks) * rActualMetric.a1(1);
                        Dw_Ds(2) = r_N(IntegrationPointIndex, ks) * rActualMetric.a1(2);
                    }
                    else {
                        Dw_Ds(0) = r_N(IntegrationPointIndex, ks) * rActualMetric.a2(0);
                        Dw_Ds(1) = r_N(IntegrationPointIndex, ks) * rActualMetric.a2(1);
                        Dw_Ds(2) = r_N(IntegrationPointIndex, ks) * rActualMetric.a2(2);
                    }
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        ddE_cu[0] = 0.5 * Dw_Ds(dirr) * r_DN_De(kr, 1);
                        ddE_cu[1] = 0.5 * Dw_Ds(dirr) * r_DN_De(kr, 0);
                    }
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        ddE_cu[0] += 0.5 * Dw_Dr(dirs, r) * r_DN_De(ks, 1);
                        ddE_cu[1] += 0.5 * Dw_Dr(dirs, r) * r_DN_De(ks, 0);
                    }
                    if (dirr == 3 && (dirs == 0 || dirs == 1 || dirs == 2))
                        DDw_DDrs(dirs) += r_N(IntegrationPointIndex, kr) * r_DN_De(ks, 0);
                    else if (dirr == 4 && (dirs == 0 || dirs == 1 || dirs == 2))
                        DDw_DDrs(dirs) += r_N(IntegrationPointIndex, kr) * r_DN_De(ks, 1);
                    if (dirs == 3 && (dirr == 0 || dirr == 1 || dirr == 2))
                        DDw_DDrs(dirr) += r_N(IntegrationPointIndex, ks) * r_DN_De(kr, 0);
                    else if (dirs == 4 && (dirr == 0 || dirr == 1 || dirr == 2))
                        DDw_DDrs(dirr) += r_N(IntegrationPointIndex, ks) * r_DN_De(kr, 1);
                    ddE_cu[0] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.a2);
                    ddE_cu[1] += 0.5 * inner_prod(DDw_DDrs, rActualMetric.a1);

                    // calculated with reduced mInitialTransConToCar
                    rSecondVariations.B23(r, s) += mInitialTransConToCar(3, 3) * ddE_cu[0] + mInitialTransConToCar(3, 4) * ddE_cu[1];
                    rSecondVariations.B13(r, s) += mInitialTransConToCar(4, 4) * ddE_cu[1];
                    if (r != s){
                        rSecondVariations.B23(s, r) = rSecondVariations.B23(r, s);
                        rSecondVariations.B13(s, r) = rSecondVariations.B13(r, s);
                    }
                    // all other entries are (remain) zero
                    
                    // 4. Second curvature variation
                    array_1d <double, 3> DDDw_DDD1rs = ZeroVector(3);
                    array_1d <double, 3> DDDw_DDD2rs = ZeroVector(3);
                    array_1d <double, 3> ddK_cu = ZeroVector(3);
                    array_1d<double, 3> DDw_DD1s = ZeroVector(3);
                    array_1d<double, 3> DDw_DD2s = ZeroVector(3);

                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        DDw_DD1s[dirs] = rDw_alpha_Dbeta(0, 0) * r_DN_De(ks, 0)  + rw_alpha[0] * r_DDN_DDe(ks, 0)
                        + rDw_alpha_Dbeta(1, 0) * r_DN_De(ks, 1) + rw_alpha[1] * r_DDN_DDe(ks, 1);
                        DDw_DD2s[dirs] = rDw_alpha_Dbeta(0, 1) * r_DN_De(ks, 0) + rw_alpha[0] * r_DDN_DDe(ks, 1)
                        + rDw_alpha_Dbeta(1, 1) * r_DN_De(ks, 1) + rw_alpha[1] * r_DDN_DDe(ks, 2);
                    }
                    else if (dirs == 3){
                        DDw_DD1s += r_DN_De(ks, 0) * rActualMetric.a1;
                        DDw_DD1s[0] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(0, 0);
                        DDw_DD1s[1] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(1, 0);
                        DDw_DD1s[2] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(2, 0);
                        DDw_DD2s += r_DN_De(ks, 1) * rActualMetric.a1;
                        DDw_DD2s[0] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(0, 2);
                        DDw_DD2s[1] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(1, 2);
                        DDw_DD2s[2] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(2, 2);
                    }
                    else if (dirs == 4){
                        DDw_DD1s += r_DN_De(ks, 0) * rActualMetric.a2;
                        DDw_DD1s[0] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(0, 2);
                        DDw_DD1s[1] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(1, 2);
                        DDw_DD1s[2] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(2, 2);
                        DDw_DD2s += r_DN_De(ks, 1) * rActualMetric.a2;
                        DDw_DD2s[0] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(0, 1);
                        DDw_DD2s[1] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(1, 1);
                        DDw_DD2s[2] += r_N(IntegrationPointIndex, ks) * rActualMetric.H(2, 1);
                    }
                    
                    if (dirr == 0 || dirr == 1 || dirr == 2){
                        ddK_cu[0] = mZeta * thickness /2.0 * DDw_DD1s[dirr] * r_DN_De(kr, 0);
                        ddK_cu[1] = mZeta * thickness /2.0 * DDw_DD2s[dirr] * r_DN_De(kr, 1);
                        ddK_cu[2] = mZeta * thickness /2.0 * 0.5 * (DDw_DD1s[dirr] * r_DN_De(kr, 1) + DDw_DD2s[dirr] * r_DN_De(kr, 0));
                        }
                    else if (dirr == 3 && (dirs == 0 || dirs == 1 || dirs == 2)){
                        DDDw_DDD1rs[dirs] = r_DN_De(kr, 0) * r_DN_De(ks, 0) + r_N(IntegrationPointIndex, kr) * r_DDN_DDe(ks, 0);
                        DDDw_DDD2rs[dirs] = r_DN_De(kr, 1) * r_DN_De(ks, 0) + r_N(IntegrationPointIndex, kr) * r_DDN_DDe(ks, 1);
                    }
                    else if (dirr == 4 && (dirs == 0 || dirs == 1 || dirs == 2)){
                        DDDw_DDD1rs[dirs] += r_DN_De(kr, 0) * r_DN_De(ks, 1) + r_N(IntegrationPointIndex, kr) * r_DDN_DDe(ks, 1);
                        DDDw_DDD2rs[dirs] += r_DN_De(kr, 1) * r_DN_De(ks, 1) + r_N(IntegrationPointIndex, kr) * r_DDN_DDe(ks, 2);
                    }
                    if (dirs == 0 || dirs == 1 || dirs == 2){
                        ddK_cu[0] += mZeta * thickness /2.0 * DDw_DD1r[dirs] * r_DN_De(ks, 0);
                        ddK_cu[1] += mZeta * thickness /2.0 * DDw_DD2r[dirs] * r_DN_De(ks, 1);
                        ddK_cu[2] += mZeta * thickness /2.0 * 0.5 * (DDw_DD1r[dirs] * r_DN_De(ks, 1) + DDw_DD2r[dirs] * r_DN_De(ks, 0));
                        }
                    else if (dirs == 3 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += r_DN_De(ks, 0) * r_DN_De(kr, 0) + r_N(IntegrationPointIndex, ks) * r_DDN_DDe(kr, 0);
                        DDDw_DDD2rs[dirr] += r_DN_De(ks, 1) * r_DN_De(kr, 0) + r_N(IntegrationPointIndex, ks) * r_DDN_DDe(kr, 1);
                    }
                    else if (dirs == 4 && (dirr == 0 || dirr == 1 || dirr == 2)){
                        DDDw_DDD1rs[dirr] += r_DN_De(ks, 0) * r_DN_De(kr, 1) + r_N(IntegrationPointIndex, ks) * r_DDN_DDe(kr, 1);
                        DDDw_DDD2rs[dirr] += r_DN_De(ks, 1) * r_DN_De(kr, 1) + r_N(IntegrationPointIndex, ks) * r_DDN_DDe(kr, 2);
                    }
                    ddK_cu[0] += mZeta * thickness / 2.0 * inner_prod(DDDw_DDD1rs, rActualMetric.a1);
                    ddK_cu[1] += mZeta * thickness / 2.0 * inner_prod(DDDw_DDD2rs, rActualMetric.a2);
                    ddK_cu[2] += mZeta * thickness / 2.0 * 0.5 * (inner_prod(DDDw_DDD1rs, rActualMetric.a2) + inner_prod(DDDw_DDD2rs, rActualMetric.a1));

                    // calculated with reduced mInitialTransConToCar
                    rSecondVariations.B11(r, s) += mInitialTransConToCar(0, 0) * ddK_cu[0];
                    rSecondVariations.B22(r, s) += mInitialTransConToCar(1, 0) * ddK_cu[0] + mInitialTransConToCar(1, 1) * ddK_cu[1]
                        + mInitialTransConToCar(1, 2) * ddK_cu[2];
                    rSecondVariations.B12(r, s) += mInitialTransConToCar(2, 0) * ddK_cu[0] + mInitialTransConToCar(2, 2) * ddK_cu[2];
                    if (r != s){
                        rSecondVariations.B11(s, r) = rSecondVariations.B11(r, s);
                        rSecondVariations.B22(s, r) = rSecondVariations.B22(r, s);
                        rSecondVariations.B12(s, r) = rSecondVariations.B12(r, s);
                    }
                    // all other entries are (remain) zero
                    // if(Id()==5 && mcount==1)
                    //     KRATOS_WATCH(rSecondVariations.B11)
                }
            }
        } // r-loop
    }

    void Shell5pHierarchicElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        // Create constitutive law parameters
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        // Set constitutive law flags
        Flags& ConstitutiveLawOptions = Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // shear difference vector
        array_1d<double, 3> w(3, 0.0);
        // derivatives of the shear difference vector
        array_1d<double, 3> Dw_D1(3, 0.0), Dw_D2(3, 0.0);
        // components w_alpha of the shear difference vector which calculates as (w_alpha(1) * a1 + w_alpha(2) * a2)
        array_1d<double, 2> w_alpha(2, 0.0);
        // derivatives of the components w_alpha
        Matrix Dw_alpha_Dbeta = ZeroMatrix(2, 2);

        std::vector<array_1d<double, 5>> stress_pk2_cart(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_pk2_cov(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_cau_cov(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 5>> stress_cau_cart(mGaussIntegrationThickness.num_GP_thickness);

        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);
        CalculateShearDifferenceVector(w, Dw_D1, Dw_D2, w_alpha, Dw_alpha_Dbeta, actual_metric);

        // the Gauss-Points start from bottom to top
        for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
        {
            mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

            array_1d<double, 3> G1(3, 0.0), G2(3, 0.0), G1_con(3, 0.0), G2_con(3, 0.0), g1(3, 0.0), g2(3, 0.0), g3(3, 0.0);
            Matrix F = ZeroMatrix(3, 3);
            double detF = 0.0;
            CalculateInitialBaseVectorsLinearised(G1, G2, G1_con, G2_con);
            CalculateActualBaseVectorsLinearised(actual_metric, w, Dw_D1, Dw_D2, g1, g2, g3);
            CalculateDeformationGradient(G1, G2, g1, g2, g3, F, detF);
            
            // Transformation matrices
            Matrix initial_trans_car_to_cov = ZeroMatrix(5, 5);
            Matrix actual_trans_cov_to_car = ZeroMatrix(5, 5);
            CalculateTransformationMatrixInitialTransConToCar(G1_con, G2_con);
            CalculateTransformationMatrixInitialTransCarToCov(initial_trans_car_to_cov);
            CalculateTransformationMatrixActualTransCovToCar(actual_trans_cov_to_car, g1, g2, g3, actual_metric.a2_con, actual_metric.a3_kirchhoff_love);

            ConstitutiveVariables constitutive_variables(5);
            CalculateConstitutiveVariables(actual_metric, w, 
                Dw_D1, Dw_D2, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

            // stresses at GP
            stress_pk2_cart[Gauss_index] = constitutive_variables.S;
            stress_pk2_cov[Gauss_index] = prod(initial_trans_car_to_cov, stress_pk2_cart[Gauss_index]);
            stress_cau_cov[Gauss_index] = stress_pk2_cov[Gauss_index] / detF;
            stress_cau_cart[Gauss_index] = prod(actual_trans_cov_to_car, stress_cau_cov[Gauss_index]);
        }   // loop GP_thickness
    
        // Cauchy stress at midspan
        array_1d<double, 5> stress_cau_cart_mid;
        stress_cau_cart_mid = (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness-1] + stress_cau_cart[0]) / 2.0;

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
        {

            if (rVariable == CAUCHY_STRESS_TOP_XX) {
                rOutput[point_number] = stress_cau_cart_mid[0] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][0]
                    - stress_cau_cart_mid[0]) / mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_YY) {
                rOutput[point_number] = stress_cau_cart_mid[1] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][1]
                    - stress_cau_cart_mid[1]) / mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_XY) {
                rOutput[point_number] = stress_cau_cart_mid[2] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][2]
                    - stress_cau_cart_mid[2]) / mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XX) {
                rOutput[point_number] = stress_cau_cart_mid[0] + (stress_cau_cart[0][0] - stress_cau_cart_mid[0]) /
                    mGaussIntegrationThickness.zeta(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_YY) {
                rOutput[point_number] = stress_cau_cart_mid[1] + (stress_cau_cart[0][1] - stress_cau_cart_mid[1]) /
                    mGaussIntegrationThickness.zeta(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XY) {
                rOutput[point_number] = stress_cau_cart_mid[2] + (stress_cau_cart[0][2] - stress_cau_cart_mid[2]) /
                    mGaussIntegrationThickness.zeta(0);
            }
            else if (rVariable == MEMBRANE_FORCE_XX) {
                rOutput[point_number] = stress_cau_cart_mid[0] * GetProperties()[THICKNESS];
            }
            else if (rVariable == MEMBRANE_FORCE_YY) {
                rOutput[point_number] = stress_cau_cart_mid[1] * GetProperties()[THICKNESS];
            }
            else if (rVariable == MEMBRANE_FORCE_XY) {
                rOutput[point_number] = stress_cau_cart_mid[2] * GetProperties()[THICKNESS];
            }
            else if (rVariable == INTERNAL_MOMENT_XX) {
                rOutput[point_number] = (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][0] - stress_cau_cart_mid[0])
                    * pow(GetProperties()[THICKNESS], 2) / (mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1) * 6);
            }
            else if (rVariable == INTERNAL_MOMENT_YY) {
                rOutput[point_number] = (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][1] - stress_cau_cart_mid[1]) * pow(GetProperties()[THICKNESS], 2) /
                    (mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1) * 6);
            }
            else if (rVariable == INTERNAL_MOMENT_XY) {
                rOutput[point_number] = (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][2] - stress_cau_cart_mid[2]) * pow(GetProperties()[THICKNESS], 2) /
                    (mGaussIntegrationThickness.zeta(mGaussIntegrationThickness.num_GP_thickness - 1) * 6);
            }
            else if (rVariable == SHEAR_FORCE_1) {
                rOutput[point_number] = GetProperties()[THICKNESS] * stress_cau_cart_mid[4];
            }
            else if (rVariable == SHEAR_FORCE_2) {
                rOutput[point_number] = GetProperties()[THICKNESS] * stress_cau_cart_mid[3];
            }
            else {
                KRATOS_WATCH("No results for desired variable available in Calculate of Shell5pHierarchicElement.")
            }
        }

        KRATOS_CATCH("");
    }

    void Shell5pHierarchicElement::CalculateHessian(
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

            Hessian(0, 0) += rDDN_DDe(k, 0) * coords[0];
            Hessian(0, 1) += rDDN_DDe(k, 2) * coords[0];
            Hessian(0, 2) += rDDN_DDe(k, 1) * coords[0];

            Hessian(1, 0) += rDDN_DDe(k, 0) * coords[1];
            Hessian(1, 1) += rDDN_DDe(k, 2) * coords[1];
            Hessian(1, 2) += rDDN_DDe(k, 1) * coords[1];

            Hessian(2, 0) += rDDN_DDe(k, 0) * coords[2];
            Hessian(2, 1) += rDDN_DDe(k, 2) * coords[2];
            Hessian(2, 2) += rDDN_DDe(k, 1) * coords[2];
        }
    }

    void Shell5pHierarchicElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY

        SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 5 * number_of_control_points)
            rResult.resize(5 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            IndexType index = i * 5;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector)
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        }

        KRATOS_CATCH("")
    }

    void Shell5pHierarchicElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY

        SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector)
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
        }

        KRATOS_CATCH("")
    }

    int Shell5pHierarchicElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {

        // Check whether ConstitutiveLaw is 3D
        KRATOS_ERROR_IF(mConstitutiveLawVector[0]->GetStrainSize() != 6) << "ConstitutiveLaw is not 3D." << std::endl;

        return 0;
    }

} // Namespace Kratos
