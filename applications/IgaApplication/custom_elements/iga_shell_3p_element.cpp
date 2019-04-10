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

        //Constitutive Law initialisation
        BaseDiscreteElement::Initialize();

        CalculateMetric(mInitialMetric);

        // KRATOS_WATCH(mInitialMetric.g3)

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
        
        // if (m_phi1 == 1)
        // KRATOS_WATCH("here: CalculateAllStart");

        
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
            GetGeometry(), DN_De, 3, 3, rMetric.J);

        rMetric.g1[0] = rMetric.J(0, 0);
        rMetric.g2[0] = rMetric.J(0, 1);
        rMetric.g1[1] = rMetric.J(1, 0);
        rMetric.g2[1] = rMetric.J(1, 1);
        rMetric.g1[2] = rMetric.J(2, 0);
        rMetric.g2[2] = rMetric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(rMetric.g3, rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3);
        //normalized basis vector g3_norm
        Vector g3_norm = rMetric.g3 / rMetric.dA;

        //GetCovariantMetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        rMetric.curvature[0] = rMetric.H(0, 0) * g3_norm[0] + rMetric.H(1, 0) * g3_norm[1] + rMetric.H(2, 0) * g3_norm[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * g3_norm[0] + rMetric.H(1, 1) * g3_norm[1] + rMetric.H(2, 1) * g3_norm[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * g3_norm[0] + rMetric.H(1, 2) * g3_norm[1] + rMetric.H(2, 2) * g3_norm[2];

        //contravariant rMetric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.gab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.gab_con[1] = invdetGab*rMetric.gab[0];

        array_1d<double, 3> g_con_1 = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];

        //local cartesian coordinates
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        // transformation matrix Q from curvilinear to local cartesian coordinate system
        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g_con_1);
        mG(0, 1) = inner_prod(e1, g_con_2);
        mG(1, 0) = inner_prod(e2, g_con_1);
        mG(1, 1) = inner_prod(e2, g_con_2);

        rMetric.Q(0, 0) = pow(mG(0, 0), 2);
        rMetric.Q(0, 1) = pow(mG(0, 1), 2);
        rMetric.Q(0, 2) = 2.00 * mG(0, 0) * mG(0, 1);

        rMetric.Q(1, 0) = pow(mG(1, 0), 2);
        rMetric.Q(1, 1) = pow(mG(1, 1), 2);
        rMetric.Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

        rMetric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        rMetric.Q(2, 1) = 2.00 * mG(0, 1) * mG(1, 1);
        rMetric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));
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
        
        // if (m_phi1 == 1)
        //    KRATOS_WATCH(mInitialMetric.Q);

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
        rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
            trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.S = prod(
            trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);
    }

    void IgaShell3pElement::CalculateStrain(
        Vector& StrainVector,
        const Vector& rgab)
    {
        KRATOS_TRY

        StrainVector[0] = 0.5 * (rgab[0] - mInitialMetric.gab[0]);
        StrainVector[1] = 0.5 * (rgab[1] - mInitialMetric.gab[1]);
        StrainVector[2] = 0.5 * (rgab[2] - mInitialMetric.gab[2]);

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateCurvature(
        Vector& CurvatureVector,
        const Vector& rCurvature)
    {
        KRATOS_TRY

        CurvatureVector[0] = (rCurvature[0] - mInitialMetric.curvature[0]);
        CurvatureVector[1] = (rCurvature[1] - mInitialMetric.curvature[1]);
        CurvatureVector[2] = (rCurvature[2] - mInitialMetric.curvature[2]);

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
        // KRATOS_WATCH(rB)
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

            //normal vector n
            Vector n = rMetric.g3 / rMetric.dA;

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
                    double g3dg3lg3 = (rMetric.g3[0] * dg3(j, 0) + rMetric.g3[1] * dg3(j, 1) + rMetric.g3[2] * dg3(j, 2))*inddA3;

                    dn(j, 0) = dg3(j, 0)*invdA - rMetric.g3[0] * g3dg3lg3;
                    dn(j, 1) = dg3(j, 1)*invdA - rMetric.g3[1] * g3dg3lg3;
                    dn(j, 2) = dg3(j, 2)*invdA - rMetric.g3[2] * g3dg3lg3;
                }

                // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
                b(0, index) = 0 - (DDN_DDe(i, 0) * n[0] + rMetric.H(0, 0)*dn(0, 0) + rMetric.H(1, 0)*dn(0, 1) + rMetric.H(2, 0)*dn(0, 2));
                b(0, index + 1) = 0 - (DDN_DDe(i, 0) * n[1] + rMetric.H(0, 0)*dn(1, 0) + rMetric.H(1, 0)*dn(1, 1) + rMetric.H(2, 0)*dn(1, 2));
                b(0, index + 2) = 0 - (DDN_DDe(i, 0) * n[2] + rMetric.H(0, 0)*dn(2, 0) + rMetric.H(1, 0)*dn(2, 1) + rMetric.H(2, 0)*dn(2, 2));

                //second line
                b(1, index) = 0 - (DDN_DDe(i, 1) * n[0] + rMetric.H(0, 1)*dn(0, 0) + rMetric.H(1, 1)*dn(0, 1) + rMetric.H(2, 1)*dn(0, 2));
                b(1, index + 1) = 0 - (DDN_DDe(i, 1) * n[1] + rMetric.H(0, 1)*dn(1, 0) + rMetric.H(1, 1)*dn(1, 1) + rMetric.H(2, 1)*dn(1, 2));
                b(1, index + 2) = 0 - (DDN_DDe(i, 1) * n[2] + rMetric.H(0, 1)*dn(2, 0) + rMetric.H(1, 1)*dn(2, 1) + rMetric.H(2, 1)*dn(2, 2));

                //third line
                b(2, index) = 0 - (DDN_DDe(i, 2) * n[0] + rMetric.H(0, 2)*dn(0, 0) + rMetric.H(1, 2)*dn(0, 1) + rMetric.H(2, 2)*dn(0, 2));
                b(2, index + 1) = 0 - (DDN_DDe(i, 2) * n[1] + rMetric.H(0, 2)*dn(1, 0) + rMetric.H(1, 2)*dn(1, 1) + rMetric.H(2, 2)*dn(1, 2));
                b(2, index + 2) = 0 - (DDN_DDe(i, 2) * n[2] + rMetric.H(0, 2)*dn(2, 0) + rMetric.H(1, 2)*dn(2, 1) + rMetric.H(2, 2)*dn(2, 2));
            }

            rB = - prod(mInitialMetric.Q, b);
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

                S_g3dg3[r] = rMetric.g3[0] * S_dg3(0, r) + rMetric.g3[1] * S_dg3(1, r) + rMetric.g3[2] * S_dg3(2, r);
                S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

                S_dn(0, r) = S_dg3(0, r)*inv_lg3 - rMetric.g3[0] * S_g3dg3lg3_3[r];
                S_dn(1, r) = S_dg3(1, r)*inv_lg3 - rMetric.g3[1] * S_g3dg3lg3_3[r];
                S_dn(2, r) = S_dg3(2, r)*inv_lg3 - rMetric.g3[2] * S_g3dg3lg3_3[r];
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

                    double c = -(ddg3[0] * rMetric.g3[0] + ddg3[1] * rMetric.g3[1] + ddg3[2] * rMetric.g3[2]
                        + S_dg3(0, r)*S_dg3(0, s) + S_dg3(1, r)*S_dg3(1, s) + S_dg3(2, r)*S_dg3(2, s)
                        )*inv_lg3_3;

                    double d = 3.0*S_g3dg3[r] * S_g3dg3[s] * inv_lg3_5;

                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d)*rMetric.g3[0];
                    ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d)*rMetric.g3[1];
                    ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d)*rMetric.g3[2];

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
