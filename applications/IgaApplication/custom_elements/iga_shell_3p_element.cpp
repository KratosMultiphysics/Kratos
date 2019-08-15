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
#include "custom_elements/iga_shell_3p_element.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell3pElement::Initialize()
    {
        KRATOS_TRY

        // KRATOS_WATCH("start: Initialize")
        // Constitutive Law initialisation
        BaseDiscreteElement::Initialize();
        // Check whether ConstitutiveLaw is 2D
        if (mConstitutiveLawVector[0]->GetStrainSize() != 3){
            KRATOS_WATCH("ConstitutiveLaw is not 2D.")
            KRATOS_ERROR << "ConstitutiveLaw is not 2D." << std::endl;
        }

        CalculateMetric(mInitialMetric);

        // KRATOS_WATCH("end: Initialize")

        KRATOS_CATCH("")
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

        // KRATOS_WATCH("here: CalculateAllStart")

        // definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 3;

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

        //reading in of shape function derivatives
        const Matrix&  DN_De  = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane, constitutive_variables_curvature,
            Values, ConstitutiveLaw::StressMeasure_PK2);
        if(Id()==1 && mcount ==1){
            KRATOS_WATCH(constitutive_variables_membrane.S)
            KRATOS_WATCH(constitutive_variables_curvature.S)
        }
        // calculate B MATRICES
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        Matrix BCurvature = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);
        CalculateBCurvature(BCurvature, actual_metric);
        // if (Id() == 4)
        //     KRATOS_WATCH(BMembrane)

        double integration_weight = GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA;
        // if (Id() == 4)
        //     KRATOS_WATCH(mInitialMetric.dA)

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrainCurvature(
                second_variations_strain,
                second_variations_curvature,
                actual_metric);

            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
            // KRATOS_WATCH(BMembrane)
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
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.S);
        }

        KRATOS_CATCH("");
    }

    void IgaShell3pElement::CalculateAndAddKm(
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

    void IgaShell3pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        for (int n = 0; n < mat_size; n++)
        {
            for (int m = 0; m <= n; m++)
            {
                double nm = (SD[0] * SecondVariationsStrain.B11(n, m)
                    + SD[1] * SecondVariationsStrain.B22(n, m)
                    + SD[2] * SecondVariationsStrain.B12(n, m))*rIntegrationWeight;

                rLeftHandSideMatrix(n, m) += nm;
                if(n!=m)
                    rLeftHandSideMatrix(m, n) += nm;
            }
        }
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateMetric(
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

        //basis vector g3_tilde
        MathUtils<double>::CrossProduct(rMetric.g3_tilde, rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3_tilde);
        //normalized basis vector g3
        rMetric.g3 = rMetric.g3_tilde / rMetric.dA;

        //GetCovariantMetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        // curvature
        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.g3[0] + rMetric.H(1, 0) * rMetric.g3[1] + rMetric.H(2, 0) * rMetric.g3[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.g3[0] + rMetric.H(1, 1) * rMetric.g3[1] + rMetric.H(2, 1) * rMetric.g3[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.g3[0] + rMetric.H(1, 2) * rMetric.g3[1] + rMetric.H(2, 2) * rMetric.g3[2];

        if (Has(SHAPE_FUNCTION_LOCAL_THIRD_DERIVATIVES))
        {
            // second derivatives of the base vectors w.r.t. the curvilinear coords.
            const Matrix& DDDN_DDDe = GetValue(SHAPE_FUNCTION_LOCAL_THIRD_DERIVATIVES);
            IgaGeometryUtilities::CalculateSecondDerivativesOfBaseVectors(
                GetGeometry(),
                DDDN_DDDe,
                rMetric.DDg1_DD11,
                rMetric.DDg1_DD12,
                rMetric.DDg2_DD21,
                rMetric.DDg2_DD22);
            
            // derivative of the curvature
            array_1d<double, 3> cross1;
            array_1d<double, 3> cross2;
            array_1d<double, 3> Dg1_D1;
            array_1d<double, 3> Dg2_D2;
            array_1d<double, 3> Dg1_D2;
            for (unsigned int i = 0; i < 3; i++)
            {
                Dg1_D1[i] = rMetric.H(i, 0);
                Dg2_D2[i] = rMetric.H(i, 1);
                Dg1_D2[i] = rMetric.H(i, 2);
            }
            MathUtils<double>::CrossProduct(cross1, Dg1_D1, rMetric.g2);
            MathUtils<double>::CrossProduct(cross2, rMetric.g1, Dg1_D2);
            array_1d<double, 3> Dg3_D1 = (rMetric.dA * (cross1 + cross2) - rMetric.g3_tilde * norm_2(cross1 + cross2)) 
                / pow(rMetric.dA, 2);
            MathUtils<double>::CrossProduct(cross1, Dg1_D2, rMetric.g2);
            MathUtils<double>::CrossProduct(cross2, rMetric.g1, Dg2_D2);
            array_1d<double, 3> Dg3_D2 = (rMetric.dA * (cross1 + cross2) - rMetric.g3_tilde * norm_2(cross1 + cross2))
                / pow(rMetric.dA, 2);
            rMetric.Dcurvature_D1[0] = inner_prod(rMetric.DDg1_DD11, rMetric.g3) + inner_prod(Dg1_D1, Dg3_D1);
            rMetric.Dcurvature_D1[1] = inner_prod(rMetric.DDg2_DD21, rMetric.g3) + inner_prod(Dg2_D2, Dg3_D1);
            rMetric.Dcurvature_D1[2] = inner_prod(rMetric.DDg1_DD12, rMetric.g3) + inner_prod(Dg1_D2, Dg3_D1);
            rMetric.Dcurvature_D2[0] = inner_prod(rMetric.DDg1_DD12, rMetric.g3) + inner_prod(Dg1_D1, Dg3_D2);
            rMetric.Dcurvature_D2[1] = inner_prod(rMetric.DDg2_DD22, rMetric.g3) + inner_prod(Dg2_D2, Dg3_D2);
            rMetric.Dcurvature_D2[2] = inner_prod(rMetric.DDg2_DD21, rMetric.g3) + inner_prod(Dg1_D2, Dg3_D2);
        }

        //contravariant rMetric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.gab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.gab_con[1] = invdetGab*rMetric.gab[0];

        array_1d<double, 3> g1_con = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g2_con = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];

        //local cartesian coordinates
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g2_con);
        array_1d<double, 3> e2 = g2_con / lg_con2;

        // transformation matrix Q from contravariant to local cartesian coordinate system
        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g1_con);
        mG(0, 1) = inner_prod(e1, g2_con);
        mG(1, 0) = inner_prod(e2, g1_con);
        mG(1, 1) = inner_prod(e2, g2_con);

        rMetric.Q(0, 0) = pow(mG(0, 0), 2);
        rMetric.Q(0, 1) = pow(mG(0, 1), 2);
        rMetric.Q(0, 2) = 2.00 * mG(0, 0) * mG(0, 1);

        rMetric.Q(1, 0) = pow(mG(1, 0), 2);
        rMetric.Q(1, 1) = pow(mG(1, 1), 2);
        rMetric.Q(1, 2) = 2.00 * mG(1, 0) * mG(1, 1);

        rMetric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        rMetric.Q(2, 1) = 2.00 * mG(0, 1) * mG(1, 1);
        rMetric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1) * mG(1, 0));

        // transformation matrix TransCartToCov from local Cartesian to covariant basis
        rMetric.TransCartToCov = trans(rMetric.Q);
        for (unsigned int i = 0; i < 3; i++){
            rMetric.TransCartToCov(2, i) = rMetric.TransCartToCov(i, 2) / 2;    // because it is not used for strains
        }

        double mG_00 = inner_prod(e1, rMetric.g1);
        double mG_01 = inner_prod(e1, rMetric.g2);
        double mG_11 = inner_prod(e2, rMetric.g2);
        // transformation matrix TransCovToCart from covariant to local Cartesian basis
        rMetric.TransCovToCart(0, 0) = pow(mG_00, 2);
        rMetric.TransCovToCart(0, 1) = pow(mG_01, 2);
        rMetric.TransCovToCart(0, 2) = 2 * mG_00 * mG_01;
        rMetric.TransCovToCart(1, 1) = pow(mG_11, 2);
        rMetric.TransCovToCart(2, 1) = mG_01 * mG_11;
        rMetric.TransCovToCart(2, 2) = mG_00 * mG_11;
    }

    void IgaShell3pElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        Vector strain_vector = ZeroVector(3);
        Vector curvature_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualMetric.gab);
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.Q, strain_vector);
        CalculateCurvature(curvature_vector, rActualMetric.curvature);
        rThisConstitutiveVariablesCurvature.E = prod(mInitialMetric.Q, curvature_vector);

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

    void IgaShell3pElement::CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab)
    {
        KRATOS_TRY

        rStrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.gab[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.gab[1]);
        rStrainVector[2] = 0.5 * (rgab[2] - mInitialMetric.gab[2]);

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature)
    {
        KRATOS_TRY

        rCurvatureVector[0] = (mInitialMetric.curvature[0] - rCurvature[0]);
        rCurvatureVector[1] = (mInitialMetric.curvature[1] - rCurvature[1]);
        rCurvatureVector[2] = (mInitialMetric.curvature[2] - rCurvature[2]);

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateBMembrane(
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
    }

    void IgaShell3pElement::CalculateBCurvature(
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

                for (unsigned int j = 0; j < 3; j++)
                {
                    double g3dg3lg3 = (rMetric.g3_tilde[0] * dg3(j, 0) + rMetric.g3_tilde[1] * dg3(j, 1) + rMetric.g3_tilde[2] * dg3(j, 2))*inddA3;

                    dn(j, 0) = dg3(j, 0)*invdA - rMetric.g3_tilde[0] * g3dg3lg3;
                    dn(j, 1) = dg3(j, 1)*invdA - rMetric.g3_tilde[1] * g3dg3lg3;
                    dn(j, 2) = dg3(j, 2)*invdA - rMetric.g3_tilde[2] * g3dg3lg3;
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

            rB = prod(mInitialMetric.Q, b);
        }
        else
        {
            KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
        }
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateSecondVariationStrainCurvature(
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

                S_g3dg3[r] = rMetric.g3_tilde[0] * S_dg3(0, r) + rMetric.g3_tilde[1] * S_dg3(1, r) + rMetric.g3_tilde[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.g3_tilde[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.g3_tilde[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.g3_tilde[2] * S_g3dg3lg3_3[r];
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

                    double c = -(ddg3[0] * rMetric.g3_tilde[0] + ddg3[1] * rMetric.g3_tilde[1] + ddg3[2] * rMetric.g3_tilde[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.g3_tilde[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.g3_tilde[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.g3_tilde[2];

                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = - (DDN_DDe(kr, 0)*S_dn(dirr, s) + DDN_DDe(ks, 0)*S_dn(dirs, r)       // MLFÄ
                        + rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2]);
                    ddK_cu[1] = - (DDN_DDe(kr, 1)*S_dn(dirr, s) + DDN_DDe(ks, 1)*S_dn(dirs, r)       // MLFÄ
                        + rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2]);
                    ddK_cu[2] = - (DDN_DDe(kr, 2)*S_dn(dirr, s) + DDN_DDe(ks, 2)*S_dn(dirs, r)       // MLFÄ
                        + rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2]);

                    rSecondVariationsCurvature.B11(r, s) = mInitialMetric.Q(0, 0)*ddK_cu[0] + mInitialMetric.Q(0, 1)*ddK_cu[1] + mInitialMetric.Q(0, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B22(r, s) = mInitialMetric.Q(1, 0)*ddK_cu[0] + mInitialMetric.Q(1, 1)*ddK_cu[1] + mInitialMetric.Q(1, 2)*ddK_cu[2];
                    rSecondVariationsCurvature.B12(r, s) = mInitialMetric.Q(2, 0)*ddK_cu[0] + mInitialMetric.Q(2, 1)*ddK_cu[1] + mInitialMetric.Q(2, 2)*ddK_cu[2];
                }
            }
        }
    }
    
    void IgaShell3pElement::Calculate(
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

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);
        CalculateConstitutiveVariables(actual_metric,
            constitutive_variables_membrane, constitutive_variables_curvature,
            Values, ConstitutiveLaw::StressMeasure_PK2);

        double detF = actual_metric.dA / mInitialMetric.dA; // should be checked (ML)

        // stresses
        double thickness = GetProperties().GetValue(THICKNESS);
        array_1d<double, 3> membrane_stress_pk2_cart = constitutive_variables_membrane.S / thickness;
        array_1d<double, 3> membrane_stress_pk2_cov = prod(mInitialMetric.TransCartToCov, membrane_stress_pk2_cart);
        array_1d<double, 3> membrane_stress_cau_cov = membrane_stress_pk2_cov / detF;
        array_1d<double, 3> membrane_stress_cau_cart = prod(actual_metric.TransCovToCart, membrane_stress_cau_cov);
        array_1d<double, 3> bending_stress_pk2_cart = constitutive_variables_curvature.S / pow(thickness, 3) * 12;
        array_1d<double, 3> bending_stress_pk2_cov = prod(mInitialMetric.TransCartToCov, bending_stress_pk2_cart);
        array_1d<double, 3> bending_stress_cau_cov = bending_stress_pk2_cov / detF;
        array_1d<double, 3> bending_stress_cau_cart = prod(actual_metric.TransCovToCart, bending_stress_cau_cov);

        // internal forces
        array_1d<double, 3> n = membrane_stress_cau_cart * thickness;

        // internal moments
        array_1d<double, 3> m = bending_stress_cau_cart * pow(thickness, 3) / 12;
        
        // stresses at the top (positive theta_3 direction)
        array_1d<double, 3> stress_cau_cart_top = membrane_stress_cau_cart + thickness / 2 * bending_stress_cau_cart;
        // stresses at the bottom (negative theta_3 direction)
        array_1d<double, 3> stress_cau_cart_bottom = membrane_stress_cau_cart - thickness / 2 * bending_stress_cau_cart;

        if (rVariable == STRESS_CAUCHY_11)
            rValues = membrane_stress_cau_cart[0];
        else if (rVariable == STRESS_CAUCHY_22)
            rValues = membrane_stress_cau_cart[1];
        else if (rVariable == STRESS_CAUCHY_12)
            rValues = membrane_stress_cau_cart[2];
        else if (rVariable == STRESS_CAUCHY_TOP_11)
            rValues = stress_cau_cart_top[0];
        else if (rVariable == STRESS_CAUCHY_TOP_22)
            rValues = stress_cau_cart_top[1];
        else if (rVariable == STRESS_CAUCHY_TOP_12)
            rValues = stress_cau_cart_top[2];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_11)
            rValues = stress_cau_cart_bottom[0];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_22)
            rValues = stress_cau_cart_bottom[1];
        else if (rVariable == STRESS_CAUCHY_BOTTOM_12)
            rValues = stress_cau_cart_bottom[2];
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
        //shear force
        else if (Has(SHAPE_FUNCTION_LOCAL_THIRD_DERIVATIVES)){
            array_1d<double, 2> q = ZeroVector(2);
            CalculateShearForce(q, actual_metric, constitutive_variables_curvature);
            if (rVariable == SHEAR_FORCE_1)
                rValues = q[0];
            else if (rVariable == SHEAR_FORCE_2)
                rValues = q[1];
        }
        else
            rValues = 0.0;
    }

    void IgaShell3pElement::CalculateShearForce(
        array_1d<double, 2>& rq, 
        const MetricVariables& rActualMetric,
        const ConstitutiveVariables& rConstitutiveVariablesCurvature)
    {
        std::vector<array_1d<double, 3>> dk_cu(2);
        dk_cu[0] = mInitialMetric.Dcurvature_D1 - rActualMetric.Dcurvature_D1;
        dk_cu[1] = mInitialMetric.Dcurvature_D2 - rActualMetric.Dcurvature_D2;

        Vector k_cu = ZeroVector(3);
        CalculateCurvature(k_cu, rActualMetric.curvature);
        array_1d<double, 3> m_ca = prod(rConstitutiveVariablesCurvature.D, k_cu);
        array_1d<double, 3> m_cu = prod(mInitialMetric.TransCartToCov, m_ca);

        // derivative of the transformation matrix Q (contravariant to local Cartesian basis)
        std::vector<Matrix> DQ_Dalpha_init(2, ZeroMatrix(3, 3));
        std::vector<Matrix> DTransCartToCov_Dalpha_init(2, ZeroMatrix(3, 3));
        CalculateDerivativeTransformationMatrices(DQ_Dalpha_init, DTransCartToCov_Dalpha_init);

        // derivative of the moment
        std::vector<array_1d<double, 3>> dm_ca(2);
        std::vector<array_1d<double, 3>> dm_cu(2);
        array_1d<double, 3> dk_1_ca = prod(mInitialMetric.Q, dk_cu[0]) + prod(DQ_Dalpha_init[0], k_cu);
        array_1d<double, 3> dk_2_ca = prod(mInitialMetric.Q, dk_cu[1]) + prod(DQ_Dalpha_init[1], k_cu);
        dm_ca[0] = prod(rConstitutiveVariablesCurvature.D, dk_1_ca);
        dm_ca[1] = prod(rConstitutiveVariablesCurvature.D, dk_2_ca);

        dm_cu[0] = prod(mInitialMetric.TransCartToCov, dm_ca[0]) + prod(DTransCartToCov_Dalpha_init[0], m_ca);
        dm_cu[1] = prod(mInitialMetric.TransCartToCov, dm_ca[1]) + prod(DTransCartToCov_Dalpha_init[1], m_ca);

        array_1d<double, 2> q_pk2_cov;
        q_pk2_cov[0] = dm_cu[0](0) / sqrt(rActualMetric.gab[0]) + dm_cu[1](2) / sqrt(rActualMetric.gab[1]);
        q_pk2_cov[1] = dm_cu[1](1) / sqrt(rActualMetric.gab[1]) + dm_cu[0](2) / sqrt(rActualMetric.gab[0]);
        double detF = rActualMetric.dA / mInitialMetric.dA; // should be checked (ML)
        array_1d<double, 2> q_cau_cov = q_pk2_cov / detF;

        // transformation actual covariant basis to local Cartesian basis
        array_1d<double, 3> g1_con = rActualMetric.g1 * rActualMetric.gab_con[0] + rActualMetric.g2 * rActualMetric.gab_con[2];
        array_1d<double, 3> g2_con = rActualMetric.g1 * rActualMetric.gab_con[2] + rActualMetric.g2 * rActualMetric.gab_con[1];
        //local cartesian coordinates
        double lg1 = norm_2(rActualMetric.g1);
        array_1d<double, 3> e1 = rActualMetric.g1 / lg1;
        double lg_con2 = norm_2(g2_con);
        array_1d<double, 3> e2 = g2_con / lg_con2;
        rq[0] = inner_prod(e1, g1_con) * q_cau_cov[0];
        rq[1] = inner_prod(e2, g1_con) * q_cau_cov[0] + inner_prod(e2, g2_con) * q_cau_cov[1];
        }

    void IgaShell3pElement::CalculateDerivativeTransformationMatrices(
        std::vector<Matrix>& rDQ_Dalpha_init,
        std::vector<Matrix>& rDTransCartToCov_Dalpha_init)
    {
        double invdetGab = 1.0 / (mInitialMetric.gab[0] * mInitialMetric.gab[1] - mInitialMetric.gab[2] * mInitialMetric.gab[2]);
        array_1d<double, 3> g1_con = mInitialMetric.g1 * mInitialMetric.gab_con[0] + mInitialMetric.g2 * mInitialMetric.gab_con[2];
        array_1d<double, 3> g2_con = mInitialMetric.g1 * mInitialMetric.gab_con[2] + mInitialMetric.g2 * mInitialMetric.gab_con[1];

        //local cartesian coordinates
        double lg1 = norm_2(mInitialMetric.g1);
        array_1d<double, 3> e1 = mInitialMetric.g1 / lg1;
        double lg_con2 = norm_2(g2_con);
        array_1d<double, 3> e2 = g2_con / lg_con2;

        // in 1.direction
        array_1d<double, 3> Dg1_D1_init;
        array_1d<double, 3> Dg1_D2_init;
        array_1d<double, 3> Dg2_D2_init;

        for (unsigned int i = 0; i < 3; i++)
        {
            Dg1_D1_init[i] = mInitialMetric.H(i, 0);
            Dg2_D2_init[i] = mInitialMetric.H(i, 1);
            Dg1_D2_init[i] = mInitialMetric.H(i, 2);
        }
        array_1d<double, 3> De1_D1 = Dg1_D1_init / lg1 - inner_prod(e1, Dg1_D1_init) * e1 / lg1;
        array_1d<double, 3> tilde_t2 = mInitialMetric.g2 - inner_prod(mInitialMetric.g2, e1) * e1;
        double bar_tilde_t2 = norm_2(tilde_t2);
        array_1d<double, 3> tilde_t2_1 = Dg1_D2_init - (inner_prod(Dg1_D2_init, e1) + inner_prod(mInitialMetric.g2,De1_D1)) * e1 
            - inner_prod(mInitialMetric.g2, e1) * De1_D1;
        array_1d<double, 3> De2_D1 = tilde_t2_1 / bar_tilde_t2 + inner_prod(e2, tilde_t2_1) * e2 / bar_tilde_t2;


        //derivative of covariant base vectors
        double A_1 = 2.0 * inner_prod(Dg1_D1_init, mInitialMetric.g1) * mInitialMetric.gab[1] 
            + 2.0 * mInitialMetric.gab[0] * inner_prod(Dg1_D2_init, mInitialMetric.g2) 
            - 2.0 * mInitialMetric.gab[2] * (inner_prod(Dg1_D1_init, mInitialMetric.g2) + inner_prod(mInitialMetric.g1, Dg1_D2_init)); // check (stands in Carat Code)
        array_1d<double, 3> Dg1_con_D1 = invdetGab 
            * (2.0 * inner_prod(Dg1_D2_init, mInitialMetric.g2) * mInitialMetric.g1 + mInitialMetric.gab[1] * Dg1_D1_init 
            - (inner_prod(Dg1_D1_init, mInitialMetric.g2) + inner_prod(mInitialMetric.g1, Dg1_D2_init)) 
            * mInitialMetric.g2 - mInitialMetric.gab[2] * Dg1_D2_init) 
            - pow(invdetGab, 2) * (mInitialMetric.gab[1] * mInitialMetric.g1 - mInitialMetric.gab[2] * mInitialMetric.g2) * A_1;
        array_1d<double, 3> Dg2_con_D1 = invdetGab
            * (-(inner_prod(Dg1_D2_init, mInitialMetric.g1) + inner_prod(mInitialMetric.g2, Dg1_D1_init)) * mInitialMetric.g1 
            + mInitialMetric.gab[2] * Dg1_D1_init + 2.0 * inner_prod(Dg1_D1_init, mInitialMetric.g1) * mInitialMetric.g2 
            - mInitialMetric.gab[0] * Dg1_D2_init) - pow(invdetGab, 2) 
            * (-mInitialMetric.gab[2] * mInitialMetric.g1 + mInitialMetric.gab[0] * mInitialMetric.g2) * A_1;
  
        double eG11 = inner_prod(e1, g1_con);
        double eG12 = inner_prod(e1, g2_con);
        double eG21 = inner_prod(e2, g1_con);
        double eG22 = inner_prod(e2, g2_con);
        double eG11_d1 = inner_prod(De1_D1, g1_con) + inner_prod(e1, Dg1_con_D1);
        double eG12_d1 = inner_prod(De1_D1, g2_con) + inner_prod(e1, Dg2_con_D1);
        double eG21_d1 = inner_prod(De2_D1, g1_con) + inner_prod(e2, Dg1_con_D1);
        double eG22_d1 = inner_prod(De2_D1, g2_con) + inner_prod(e2, Dg2_con_D1);
        
        // derivative of the transformation matrix Q (contravariant to local Cartesian basis) of the initial configuration w.r.t. theta1
        rDQ_Dalpha_init[0](0,0) = eG11_d1 * eG11 + eG11 * eG11_d1;
        rDQ_Dalpha_init[0](0,1) = eG12_d1 * eG12 + eG12 * eG12_d1;
        rDQ_Dalpha_init[0](0,2) = 2.0 * (eG11_d1 * eG12 + eG11 * eG12_d1);
        rDQ_Dalpha_init[0](1,0) = eG21_d1 * eG21 + eG21 * eG21_d1;  // should be always zero (ML)
        rDQ_Dalpha_init[0](1,1) = eG22_d1 * eG22 + eG22 * eG22_d1;
        rDQ_Dalpha_init[0](1,2) = 2.0 * (eG21_d1 * eG22 + eG21 * eG22_d1);  // should be always zero (ML)
        rDQ_Dalpha_init[0](2,0) = 2.0 * (eG11_d1 * eG21 + eG11 * eG21_d1);  // should be always zero (ML)
        rDQ_Dalpha_init[0](2,1) = 2.0 * (eG12_d1 * eG22 + eG12 * eG22_d1);
        rDQ_Dalpha_init[0](2,2) = 2.0 * (eG11_d1 * eG22 + eG12_d1 * eG21 + eG11 * eG22_d1 + eG12 * eG21_d1);
        // derivative of the transformation matrix TransCartToCov (local Cartesian to covariant basis) of the initial configuration 
        // w.r.t. theta1
        rDTransCartToCov_Dalpha_init[0](0,0) = eG11_d1 * eG11 + eG11 * eG11_d1;
        rDTransCartToCov_Dalpha_init[0](0,1) = eG21_d1 * eG21 + eG21 * eG21_d1;
        rDTransCartToCov_Dalpha_init[0](0,2) = 2.0 * (eG11_d1 * eG21 + eG11 * eG21_d1);
        rDTransCartToCov_Dalpha_init[0](1,0) = eG12_d1 * eG12 + eG12 * eG12_d1;
        rDTransCartToCov_Dalpha_init[0](1,1) = eG22_d1 * eG22 + eG22 * eG22_d1;
        rDTransCartToCov_Dalpha_init[0](1,2) = 2.0 * (eG12_d1 * eG22 + eG12 * eG22_d1);
        rDTransCartToCov_Dalpha_init[0](2,0) = eG11_d1 * eG12 + eG11 * eG12_d1;
        rDTransCartToCov_Dalpha_init[0](2,1) = eG21_d1 * eG22 + eG21 * eG22_d1;
        rDTransCartToCov_Dalpha_init[0](2,2) = eG11_d1 * eG22 + eG21_d1 * eG12 + eG11 * eG22_d1 + eG21 * eG12_d1;

        // in 2.direction
        array_1d<double, 3> De1_D2 = Dg1_D2_init / lg1 + inner_prod(e1, Dg1_D2_init) * e1 / lg1;
        array_1d<double, 3> tilde_t2_2 = Dg2_D2_init - (inner_prod(Dg2_D2_init, e1) + inner_prod(mInitialMetric.g2, De1_D2)) * e1 
            - inner_prod(mInitialMetric.g2, e1) * De1_D2;
        array_1d<double, 3> De2_D2 = tilde_t2_2 / bar_tilde_t2 + inner_prod(e2, tilde_t2_2) * e2 / bar_tilde_t2;
  
        //derivative of covariant base vectors
        double A_2 = 2.0 * inner_prod(Dg1_D2_init, mInitialMetric.g1) * mInitialMetric.gab[1] 
            + 2.0 * mInitialMetric.gab[0] * inner_prod(Dg2_D2_init,mInitialMetric.g2) 
            - 2.0 * mInitialMetric.gab[2] * (inner_prod(Dg1_D2_init, mInitialMetric.g2) + inner_prod(mInitialMetric.g1, Dg2_D2_init));
        array_1d<double, 3> Dg1_con_D2 = invdetGab * (2.0 * inner_prod(Dg2_D2_init, mInitialMetric.g2) * mInitialMetric.g1 
            + mInitialMetric.gab[1] * Dg1_D2_init - (inner_prod(Dg1_D2_init, mInitialMetric.g2) 
            + inner_prod(mInitialMetric.g1, Dg2_D2_init)) * mInitialMetric.g2 - mInitialMetric.gab[2] * Dg2_D2_init) 
            - pow(invdetGab, 2) * (mInitialMetric.gab[1] * mInitialMetric.g1 - mInitialMetric.gab[2] * mInitialMetric.g2) * A_2;
        array_1d<double, 3> Dg2_con_D2 = invdetGab * (-(inner_prod(Dg2_D2_init, mInitialMetric.g1) 
            + inner_prod(mInitialMetric.g2, Dg1_D2_init)) * mInitialMetric.g1 + mInitialMetric.gab[2] * Dg1_D2_init 
            + 2.0 * inner_prod(Dg1_D2_init, mInitialMetric.g1) * mInitialMetric.g2 - mInitialMetric.gab[0] * Dg2_D2_init) 
            - pow(invdetGab, 2) * (-mInitialMetric.gab[2] * mInitialMetric.g1 + mInitialMetric.gab[0] * mInitialMetric.g2) * A_2;

        double eG11_d2 = inner_prod(De1_D2, g1_con) + inner_prod(e1, Dg1_con_D2);
        double eG12_d2 = inner_prod(De1_D2, g2_con) + inner_prod(e1, Dg2_con_D2);
        double eG21_d2 = inner_prod(De2_D2, g1_con) + inner_prod(e2, Dg1_con_D2);
        double eG22_d2 = inner_prod(De2_D2, g2_con) + inner_prod(e2, Dg2_con_D2);

        // derivative of the transformation matrix Q (contravariant to local Cartesian basis) of the initial configuration w.r.t. theta2
        rDQ_Dalpha_init[1](0,0) = eG11_d2 * eG11 + eG11 * eG11_d2;
        rDQ_Dalpha_init[1](0,1) = eG12_d2 * eG12 + eG12 * eG12_d2;
        rDQ_Dalpha_init[1](0,2) = 2.0 * (eG11_d2 * eG12 + eG11 * eG12_d2);
        rDQ_Dalpha_init[1](1,0) = eG21_d2 * eG21 + eG21 * eG21_d2;
        rDQ_Dalpha_init[1](1,1) = eG22_d2 * eG22 + eG22 * eG22_d2;
        rDQ_Dalpha_init[1](1,2) = 2.0 * (eG21_d2 * eG22 + eG21 * eG22_d2);
        rDQ_Dalpha_init[1](2,0) = 2.0 * (eG11_d2 * eG21 + eG11 * eG21_d2);
        rDQ_Dalpha_init[1](2,1) = 2.0 * (eG12_d2 * eG22 + eG12 * eG22_d2);
        rDQ_Dalpha_init[1](2,2) = 2.0 * (eG11_d2 * eG22 + eG12_d2 * eG21 + eG11 * eG22_d2 + eG12 * eG21_d2);
        
        // derivative of the transformation matrix TransCartToCov (local Cartesian to covariant basis) of the initial configuration 
        // w.r.t. theta2
        rDTransCartToCov_Dalpha_init[1](0,0) = eG11_d2 * eG11 + eG11 * eG11_d2;
        rDTransCartToCov_Dalpha_init[1](0,1) = eG21_d2*eG21 + eG21*eG21_d2;
        rDTransCartToCov_Dalpha_init[1](0,2) = 2.0*(eG11_d2*eG21 + eG11*eG21_d2);
        rDTransCartToCov_Dalpha_init[1](1,0) = eG12_d2*eG12 + eG12*eG12_d2;
        rDTransCartToCov_Dalpha_init[1](1,1) = eG22_d2*eG22 + eG22*eG22_d2;
        rDTransCartToCov_Dalpha_init[1](1,2) = 2.0*(eG12_d2*eG22 + eG12*eG22_d2);
        rDTransCartToCov_Dalpha_init[1](2,0) = eG11_d2*eG12 + eG11*eG12_d2;
        rDTransCartToCov_Dalpha_init[1](2,1) = eG21_d2*eG22 + eG21*eG22_d2;
        rDTransCartToCov_Dalpha_init[1](2,2) = eG11_d2*eG22 + eG21_d2*eG12 + eG11*eG22_d2 + eG21*eG12_d2;
    }

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
