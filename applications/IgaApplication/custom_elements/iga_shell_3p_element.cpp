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
#include "custom_elements/iga_shell_3p_element.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"


namespace Kratos
{
    void IgaShell3pElement::Initialize()
    {
        KRATOS_TRY

        KRATOS_WATCH("InitializeStart");

        // Constitutive Law initialisation
        BaseDiscreteElement::Initialize();
        
        // KRATOS_WATCH(m_initial_metric.g3);

        CalculateMetric(m_initial_metric);
        m_Phi = ZeroVector(3);

        m_Dv_D1 = ZeroVector(3);
        m_Dv_D2 = ZeroVector(3);

        KRATOS_WATCH("InitializeEnd");

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
        
        KRATOS_WATCH("here: CalculateAllStart");

        // 1. definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 3;        // considers number of parameters 3

        //2. set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);        // direct initialisation of the structure 'Values' (name) of the type Parameters

        /** decides where the strains, stresses and the constitutive tensor are calculated
         * the following settings determine: strains from element, stresses and constitutive tensor from constitutive law
         */
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        // 3. resizing of solving matrices LHS and RHS      
        // resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }
        // resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }

        // 4. reading in shape function derivatives
        // rDN_De(i, m): derivative of form function at node i w.r.t. convective coordinate m
        Matrix  DN_De  = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        Matrix DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        // 5. calculate actual metric
        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        
        // 6. calculate strains
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

        // 7. calculate B MATRICES
        // size of matrix depending on number of strain entries (here: eps11, eps22, eps12) and mat_size
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        Matrix BCurvature = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric, DN_De, number_of_nodes, mat_size);
        CalculateBCurvature(BCurvature, actual_metric, DN_De, DDN_DDe, number_of_nodes, mat_size);

        double integration_weight = GetValue(INTEGRATION_WEIGHT) * m_initial_metric.dA;

        // 8. left hand side matrix
        if (CalculateStiffnessMatrixFlag == true)
        {
            // 8.1 calculate second variations
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrain(
                second_variations_strain,
                second_variations_curvature,
                actual_metric);

            // 8.2 adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(
                rLeftHandSideMatrix, 
                BMembrane, 
                constitutive_variables_membrane.D, 
                integration_weight);
            // 8.3 adding curvature contributions to the stiffness matrix
            CalculateAndAddKm(
                rLeftHandSideMatrix, 
                BCurvature, 
                constitutive_variables_curvature.D, 
                integration_weight);

            // 8.4 adding  non-linear-contribution to Stiffness-Matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
               second_variations_strain,
               constitutive_variables_membrane.S,
               integration_weight);
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
               second_variations_curvature,
               constitutive_variables_curvature.S,
               integration_weight);
        }

        // 9. left hand side matrix (only added internal forces here, external in LoadCurveDiscreteCondition or LoadPointDiscreteCondition)
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.S);
            noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.S);
        }

        //KRATOS_WATCH(rLeftHandSideMatrix)
        
        KRATOS_WATCH("here: CalculateAllEnd");

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double& rIntegrationWeight)
    {
        KRATOS_TRY
        
        // KRATOS_WATCH("CalculateAndAddKm");

        noalias(rLeftHandSideMatrix) += rIntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
        
        // KRATOS_WATCH(rLeftHandSideMatrix);
        
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateAndAddNonlinearKm(
        MatrixType& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight)

    {
        KRATOS_TRY

        KRATOS_WATCH("CalculateAndAddNonlinearKm");

        noalias(rLeftHandSideMatrix) += rIntegrationWeight * (SD[0] * SecondVariationsStrain.B11
                    + SD[1] * SecondVariationsStrain.B22
                    + SD[2] * SecondVariationsStrain.B12);

        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateMetric(
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

        // normalized basis vector g3 and differential area dA
        rMetric.g3 = MathUtils<double>::CrossProduct(rMetric.g1, rMetric.g2);
        rMetric.dA = norm_2(rMetric.g3);
        rMetric.g3 = rMetric.g3 / rMetric.dA;       // normalization

        //GetCovariantmetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        for (int i = 0; i < 3; ++i )
            {
                rMetric.Dg1_D1[i] = rMetric.H(i,1);
                rMetric.Dg1_D2[i] = rMetric.H(i,2);
                rMetric.Dg2_D1[i] = rMetric.H(i,1);
                rMetric.Dg2_D2[i] = rMetric.H(i,2);
                rMetric.Dg3_D1[i] = rMetric.H(i,1);
                rMetric.Dg3_D2[i] = rMetric.H(i,2);
            }

        KRATOS_WATCH("CalculateMetric1");

        // curvature
        rMetric.curvature[0] = rMetric.H(0, 0)*rMetric.g3[0] + rMetric.H(1, 0)*rMetric.g3[1] + rMetric.H(2, 0)*rMetric.g3[2];
        rMetric.curvature[1] = rMetric.H(0, 1)*rMetric.g3[0] + rMetric.H(1, 1)*rMetric.g3[1] + rMetric.H(2, 1)*rMetric.g3[2];
        rMetric.curvature[2] = rMetric.H(0, 2)*rMetric.g3[0] + rMetric.H(1, 2)*rMetric.g3[1] + rMetric.H(2, 2)*rMetric.g3[2];

        //contravariant metric gab_con
        double inv_det_gab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = inv_det_gab*rMetric.gab[1];
        rMetric.gab_con[1] = inv_det_gab*rMetric.gab[0];
        rMetric.gab_con[2] = -inv_det_gab*rMetric.gab[2];
        
        KRATOS_WATCH("CalculateMetric2");

        // contravariant base vectors g_con
        array_1d<double, 3> g_con_1 = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];

        // local cartesian base vectors, with the following choice its guaranteed that e1 is orthogonal to e2
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        // KRATOS_WATCH("CalculateMetric3");

        // transformation matrix T, from contravariant to local cartesian coordinate system
        double Ge11 = inner_prod(rMetric.g1,e1);
        double Ge12 = inner_prod(rMetric.g1,e2);
        double Ge21 = inner_prod(rMetric.g2,e1);
        double Ge22 = inner_prod(rMetric.g2,e2);

        // calculation of T not checked yet (ML)
        rMetric.T(0, 0) = Ge11 * Ge11;
        rMetric.T(0, 1) = Ge12 * Ge12;
        rMetric.T(0, 2) = Ge11 * Ge12;
        rMetric.T(1, 0) = Ge21 * Ge21;
        rMetric.T(1, 1) = Ge22 * Ge22;
        rMetric.T(1, 2) = Ge21 * Ge22;
        rMetric.T(2, 0) = 2.0 * Ge11 * Ge21;
        rMetric.T(2, 1) = 2.0 * Ge12 * Ge22;
        rMetric.T(2, 2) = Ge11 * Ge22 + Ge12 * Ge21;

        // KRATOS_WATCH("CalculateMetricEnd");

    }
 
    void IgaShell3pElement::CalculateRotationVector(
        const MetricVariables& rActualMetric)
    {   
        m_phi1 = MathUtils<double>::Dot3(m_Dv_D2, m_initial_metric.g3)
        / norm_2(MathUtils<double>::CrossProduct(m_initial_metric.g1,m_initial_metric.g2));
        m_phi2 = MathUtils<double>::Dot3(m_Dv_D1, m_initial_metric.g3)
        / norm_2(MathUtils<double>::CrossProduct(m_initial_metric.g1,m_initial_metric.g2));
        m_Phi = m_phi1 * m_initial_metric.g1 + m_phi2 * m_initial_metric.g2;
        m_DPhi_D1 = m_phi1 * m_initial_metric.Dg1_D1 + m_phi2 * m_initial_metric.Dg2_D1;
        m_DPhi_D2 = m_phi1 * m_initial_metric.Dg1_D2 + m_phi2 * m_initial_metric.Dg2_D2;
    }
    
    void IgaShell3pElement::CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        // membrane strain
        Vector strain_vector = ZeroVector(3);
        // strain due to bending/curvature
        Vector curvature_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualMetric.gab);
        // membrane strain w.r.t. local cartesian coordinates through transformation
        rThisConstitutiveVariablesMembrane.E = prod(m_initial_metric.T, strain_vector);
        CalculateCurvature(curvature_vector, rActualMetric.curvature);
        // strain due to bending w.r.t. local cartesian coordinates through transformation
        rThisConstitutiveVariablesCurvature.E = prod(m_initial_metric.T, curvature_vector);

        //Constitutive Matrices DMembrane and DCurvature
        // the following set functions are dereferencing -> local variables can be changed by changing rValues
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter

        //rValues.CheckAllParameters();
        /** mConstitutiveLawVector is a vector of size 1 (defined in base_discrete_element)
         * CalculateMaterialResponse calculates the stress depending on the strain and updates the constitutive matrix
         */
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

        rStrainVector[0] = 0.5 * (rgab[0] - m_initial_metric.gab[0]);
        rStrainVector[1] = 0.5 * (rgab[1] - m_initial_metric.gab[1]);
        rStrainVector[2] = (rgab[2] - m_initial_metric.gab[2]);
        
        KRATOS_WATCH(rStrainVector);
        
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature)
    {
         KRATOS_TRY

        rCurvatureVector[0] = (rCurvature[0] - m_initial_metric.curvature[0]);
        rCurvatureVector[1] = (rCurvature[1] - m_initial_metric.curvature[1]);
        rCurvatureVector[2] = 0.5 * (rCurvature[2] - m_initial_metric.curvature[2]);

        KRATOS_CATCH("")
    }
    
    void IgaShell3pElement::CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& rMetric,
        const Matrix& rDN_De,
        const unsigned int& rNumberOfNodes,
        const unsigned int& rMatSize)
    {
        KRATOS_TRY
        
        if (rB.size1() != rMatSize || rB.size2() != rMatSize)
            rB.resize(rMatSize, rMatSize);
        rB = ZeroMatrix(3, rMatSize);

        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES))
        {
            if (rB.size1() != 3 || rB.size2() != rMatSize)
                rB.resize(3, rMatSize);
            rB = ZeroMatrix(3, rMatSize);

           for (unsigned int i = 0; i < rNumberOfNodes; i++)
           {
               unsigned int index = 3 * i;

               //first line
               rB(0, index) = rDN_De(i, 0) * rMetric.g1[0];
               rB(0, index + 1) = rDN_De(i, 0) * rMetric.g1[1];
               rB(0, index + 2) = rDN_De(i, 0) * rMetric.g1[2];

               //second line
               rB(1, index) = rDN_De(i, 1) * rMetric.g2[0];
               rB(1, index + 1) = rDN_De(i, 1) * rMetric.g2[1];
               rB(1, index + 2) = rDN_De(i, 1) * rMetric.g2[2];

               //third line
               rB(2, index) = rDN_De(i, 1) * rMetric.g1[0] + rDN_De(i, 0) * rMetric.g2[0];
               rB(2, index + 1) = rDN_De(i, 1) * rMetric.g1[1] + rDN_De(i, 0) * rMetric.g2[1];
               rB(2, index + 2) = rDN_De(i, 1) * rMetric.g1[2] + rDN_De(i, 0) * rMetric.g2[2];
           }
           rB = prod(m_initial_metric.T, rB);
        }
        else
        {
           KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES" << std::endl;
        }
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateBCurvature(
        Matrix& rB,
        const MetricVariables& rMetric,
        const Matrix& rDN_De,
        const Matrix& rDDN_DDe,
        const unsigned int& rNumberOfNodes,
        const unsigned int& rMatSize)
    {
        KRATOS_TRY

        if (Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
        {
            if (rB.size1() != 3 || rB.size2() != rMatSize)
                rB.resize(3, rMatSize);
            rB = ZeroMatrix(3, rMatSize);
            
            // not-normalized normal vector g3
            Vector ng3 = ZeroVector(3);
            // crossproduct, predefined functions not possible because rMetric should be const as reference (ML)
            ng3[0] = rMetric.g1[1] * rMetric.g2[2] - rMetric.g1[2] * rMetric.g2[1];
            ng3[1] = rMetric.g1[2] * rMetric.g2[0] - rMetric.g1[0] * rMetric.g2[2];
            ng3[2] = rMetric.g1[0] * rMetric.g2[1] - rMetric.g1[1] * rMetric.g2[0];

            Matrix ndg3 = ZeroMatrix(3, 3);
            Matrix dg3 = ZeroMatrix(3, 3);

            double invdA = 1 / rMetric.dA;
            double invdA3 = 1 / std::pow(rMetric.dA, 3);

            for (int i = 0; i < rNumberOfNodes; i++)
            {
                unsigned int index = 3 * i;
                
                // equ. (5.24) from Kiendl (2011)
                //first line
                ndg3(0, 0) = 0;
                ndg3(0, 1) = -rDN_De(i, 0) * rMetric.g2[2] + rDN_De(i, 1) * rMetric.g1[2];
                ndg3(0, 2) = rDN_De(i, 0) * rMetric.g2[1] - rDN_De(i, 1) * rMetric.g1[1];

                //second line
                ndg3(1, 0) = rDN_De(i, 0) * rMetric.g2[2] - rDN_De(i, 1) * rMetric.g1[2];
                ndg3(1, 1) = 0;
                ndg3(1, 2) = -rDN_De(i, 0)*rMetric.g2[0] + rDN_De(i, 1) * rMetric.g1[0];

                //third line
                ndg3(2, 0) = -rDN_De(i, 0) * rMetric.g2[1] + rDN_De(i, 1) * rMetric.g1[1];
                ndg3(2, 1) = rDN_De(i, 0) * rMetric.g2[0] - rDN_De(i, 1) * rMetric.g1[0];
                ndg3(2, 2) = 0;

                for (unsigned int j = 0; j < 3; j++)
                {
                    // equ. (5.25) from Kiendl (2011), some operations moved with equ. (5.26) 
                    double ng3ndg3lg3 = (ng3[0] * ndg3(j, 0) + ng3[1] * ndg3(j, 1) + ng3[2] * ndg3(j, 2)) * invdA3;

                    // equ. (5.36) from Kiendl (2011), some operations moved with equ. (5.25)
                    dg3(j, 0) = ndg3(j, 0) * invdA - rMetric.g3[0] * ng3ndg3lg3;
                    dg3(j, 1) = ndg3(j, 1) * invdA - rMetric.g3[1] * ng3ndg3lg3;
                    dg3(j, 2) = ndg3(j, 2) * invdA - rMetric.g3[2] * ng3ndg3lg3;
                }

                // equ. (5.20) from Kiendl (2011)
                // curvature vector [K11,K22,K12] referred to curvilinear coordinate system
                rB(0, index) = 0 - (rDDN_DDe(i, 0) * rMetric.g3[0] + rMetric.H(0, 0)*dg3(0, 0) + rMetric.H(1, 0)*dg3(0, 1) + rMetric.H(2, 0)*dg3(0, 2));
                rB(0, index + 1) = 0 - (rDDN_DDe(i, 0) * rMetric.g3[1] + rMetric.H(0, 0)*dg3(1, 0) + rMetric.H(1, 0)*dg3(1, 1) + rMetric.H(2, 0)*dg3(1, 2));
                rB(0, index + 2) = 0 - (rDDN_DDe(i, 0) * rMetric.g3[2] + rMetric.H(0, 0)*dg3(2, 0) + rMetric.H(1, 0)*dg3(2, 1) + rMetric.H(2, 0)*dg3(2, 2));
     
                //second line
                rB(1, index) = 0 - (rDDN_DDe(i, 1) * rMetric.g3[0] + rMetric.H(0, 1)*dg3(0, 0) + rMetric.H(1, 1)*dg3(0, 1) + rMetric.H(2, 1)*dg3(0, 2));
                rB(1, index + 1) = 0 - (rDDN_DDe(i, 1) * rMetric.g3[1] + rMetric.H(0, 1)*dg3(1, 0) + rMetric.H(1, 1)*dg3(1, 1) + rMetric.H(2, 1)*dg3(1, 2));
                rB(1, index + 2) = 0 - (rDDN_DDe(i, 1) * rMetric.g3[2] + rMetric.H(0, 1)*dg3(2, 0) + rMetric.H(1, 1)*dg3(2, 1) + rMetric.H(2, 1)*dg3(2, 2));

                //third line
                rB(2, index) = 0 - (rDDN_DDe(i, 2) * rMetric.g3[0] + rMetric.H(0, 2)*dg3(0, 0) + rMetric.H(1, 2)*dg3(0, 1) + rMetric.H(2, 2)*dg3(0, 2));
                rB(2, index + 1) = 0 - (rDDN_DDe(i, 2) * rMetric.g3[1] + rMetric.H(0, 2)*dg3(1, 0) + rMetric.H(1, 2)*dg3(1, 1) + rMetric.H(2, 2)*dg3(1, 2));
                rB(2, index + 2) = 0 - (rDDN_DDe(i, 2) * rMetric.g3[2] + rMetric.H(0, 2)*dg3(2, 0) + rMetric.H(1, 2)*dg3(2, 1) + rMetric.H(2, 2)*dg3(2, 2));
            }

            rB = prod(m_initial_metric.T, rB);
        }
        else
        {
            KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
        }
        KRATOS_CATCH("")
    }

    void IgaShell3pElement::CalculateSecondVariationStrain(
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

            // not-normalized normal vector g3
            Vector ng3 = ZeroVector(3);
            // crossproduct, predefined functions not possible because rMetric should be const as reference (ML)
            ng3[0] = rMetric.g1[1] * rMetric.g2[2] - rMetric.g1[2] * rMetric.g2[1];
            ng3[1] = rMetric.g1[2] * rMetric.g2[0] - rMetric.g1[0] * rMetric.g2[2];
            ng3[2] = rMetric.g1[0] * rMetric.g2[1] - rMetric.g1[1] * rMetric.g2[0];
            double dA_3 = pow(rMetric.dA, 3);
            double dA_5 = pow(rMetric.dA, 5);
            double inv_dA = 1 / rMetric.dA;
            double inv_dA_3 = 1 / dA_3;
            double inv_dA_5 = 1 / dA_5;

            Matrix ndg3 = ZeroMatrix(3, mat_size);
            Vector ng3ndg3 = ZeroVector(mat_size);
            Vector ng3ndg3inv_dA_3 = ZeroVector(mat_size);
            Matrix dg3 = ZeroMatrix(3, mat_size);
            // first variation of strain and curvature w.r.t. dof
            for (int r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                int kr = r / 3;
                int dirr = r % 3;

                // equ. (5.24) from Kiendl (2011)
                array_1d<double, 3> dg1 = ZeroVector(3);
                array_1d<double, 3> dg2 = ZeroVector(3);
                dg1(dirr) = DN_De(kr, 0);
                dg2(dirr) = DN_De(kr, 1);
                ndg3(0, r) = dg1(1) * rMetric.g2(2) - dg1(2) * rMetric.g2(1) 
                + rMetric.g1(1) * dg2(2) - rMetric.g1(2) * dg2(1);
                ndg3(1, r) = dg1(2) * rMetric.g2(0) - dg1(0) * rMetric.g2(2) 
                + rMetric.g1(2) * dg2(0) - rMetric.g1(0) * dg2(2);
                ndg3(2, r) = dg1(0) * rMetric.g2(1) - dg1(1) * rMetric.g2(0) 
                + rMetric.g1(0) * dg2(1) - rMetric.g1(1) * dg2(0);

                ng3ndg3[r] = ng3[0] * ndg3(0, r) + ng3[1] * ndg3(1, r) + ng3[2] * ndg3(2, r);
                ng3ndg3inv_dA_3[r] = ng3ndg3[r] * inv_dA_3;

                dg3(0, r) = ndg3(0, r)*inv_dA - ng3[0] * ng3ndg3inv_dA_3[r];
                dg3(1, r) = ndg3(1, r)*inv_dA - ng3[1] * ng3ndg3inv_dA_3[r];
                dg3(2, r) = ndg3(2, r)*inv_dA - ng3[2] * ng3ndg3inv_dA_3[r];
            }

            // second variation of strain and curvature w.r.t. dofs u_r, u_s
            for (int r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                int kr = r / 3;
                int dirr = r % 3;

                // calculates only upper diagonal -> symmetric
                for (int s = 0; s <= r; s++)
                {
                    // local node number ks and dof direction dirs
                    int ks = s / 3;
                    int dirs = s % 3;

                    // strain
                    array_1d<double, 3> ddE_cu = ZeroVector(3);
                    if (dirr == dirs)
                    {
                        // equ. (5.19) from Kiendl (2011)
                        ddE_cu[0] = DN_De(kr, 0)*DN_De(ks, 0);
                        ddE_cu[1] = DN_De(kr, 1)*DN_De(ks, 1);
                        ddE_cu[2] = DN_De(kr, 0)*DN_De(ks, 1) + DN_De(kr, 1)*DN_De(ks, 0);

                        // transformation
                        rSecondVariationsStrain.B11(r, s) = m_initial_metric.T(0, 0)*ddE_cu[0] + m_initial_metric.T(0, 1)*ddE_cu[1] + m_initial_metric.T(0, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B22(r, s) = m_initial_metric.T(1, 0)*ddE_cu[0] + m_initial_metric.T(1, 1)*ddE_cu[1] + m_initial_metric.T(1, 2)*ddE_cu[2];
                        rSecondVariationsStrain.B12(r, s) = m_initial_metric.T(2, 0)*ddE_cu[0] + m_initial_metric.T(2, 1)*ddE_cu[1] + m_initial_metric.T(2, 2)*ddE_cu[2];
                    }

                    // curvature
                    // equ. (5.30) from Kiendl (2011)
                    array_1d<double, 3> ddng3 = ZeroVector(3);
                    int dirt = 3 - dirr - dirs;
                    int ddir = dirr - dirs;
                    if (ddir == -1)      ddng3(dirt) = DN_De(kr, 0) * DN_De(ks, 1) - DN_De(ks, 0) * DN_De(kr, 1);
                    else if (ddir == 2) ddng3(dirt) = DN_De(kr, 0) * DN_De(ks, 1) - DN_De(ks, 0) * DN_De(kr, 1);
                    else if (ddir == 1) ddng3(dirt) = -DN_De(kr, 0) * DN_De(ks, 1) + DN_De(ks, 0) * DN_De(kr, 1);
                    else if (ddir == -2) ddng3(dirt) = -DN_De(kr, 0) * DN_De(ks, 1) + DN_De(ks, 0) * DN_De(kr, 1);

                    double c = -(ddng3[0] * ng3[0] + ddng3[1] * ng3[1] + ddng3[2] * ng3[2]
                        + ndg3(0, r)*ndg3(0, s) + ndg3(1, r)*ndg3(1, s) + ndg3(2, r)*ndg3(2, s)
                        )*inv_dA_3;

                    double d = 3.0 * ng3ndg3[r] * ng3ndg3[s] * inv_dA_5;

                    // equ. (5.32) from Kiendl (2011)
                    array_1d<double, 3> ddn = ZeroVector(3);
                    ddn[0] = ddng3[0] * inv_dA - ng3ndg3inv_dA_3[s] * ndg3(0, r) - ng3ndg3inv_dA_3[r] * ndg3(0, s) + (c + d) * ng3[0];
                    ddn[1] = ddng3[1] * inv_dA - ng3ndg3inv_dA_3[s] * ndg3(1, r) - ng3ndg3inv_dA_3[r] * ndg3(1, s) + (c + d) * ng3[1];
                    ddn[2] = ddng3[2] * inv_dA - ng3ndg3inv_dA_3[s] * ndg3(2, r) - ng3ndg3inv_dA_3[r] * ndg3(2, s) + (c + d) * ng3[2];

                    // equ. (5.34) from Kiendl (2011)
                    array_1d<double, 3> ddK_cu = ZeroVector(3);
                    ddK_cu[0] = DDN_DDe(kr, 0) * dg3(dirr, s) + DDN_DDe(ks, 0) * dg3(dirs, r)
                        + rMetric.H(0, 0) * ddn[0] + rMetric.H(1, 0) * ddn[1] + rMetric.H(2, 0) * ddn[2];
                    ddK_cu[1] = DDN_DDe(kr, 1) * dg3(dirr, s) + DDN_DDe(ks, 1) * dg3(dirs, r)
                        + rMetric.H(0, 1) * ddn[0] + rMetric.H(1, 1) * ddn[1] + rMetric.H(2, 1) * ddn[2];
                    ddK_cu[2] = DDN_DDe(kr, 2) * dg3(dirr, s) + DDN_DDe(ks, 2) * dg3(dirs, r)
                        + rMetric.H(0, 2) * ddn[0] + rMetric.H(1, 2) * ddn[1] + rMetric.H(2, 2) * ddn[2];

                    // transformation
                    rSecondVariationsCurvature.B11(r, s) = m_initial_metric.T(0, 0) * ddK_cu[0] + m_initial_metric.T(0, 1) * ddK_cu[1] + m_initial_metric.T(0, 2) * ddK_cu[2];
                    rSecondVariationsCurvature.B11(s, r) = rSecondVariationsCurvature.B11(r, s);
                    rSecondVariationsCurvature.B22(r, s) = m_initial_metric.T(1, 0) * ddK_cu[0] + m_initial_metric.T(1, 1) * ddK_cu[1] + m_initial_metric.T(1, 2) * ddK_cu[2];
                    rSecondVariationsCurvature.B22(s, r) = rSecondVariationsCurvature.B22(r, s);
                    rSecondVariationsCurvature.B12(r, s) = m_initial_metric.T(2, 0) * ddK_cu[0] + m_initial_metric.T(2, 1) * ddK_cu[1] + m_initial_metric.T(2, 2) * ddK_cu[2];
                    rSecondVariationsCurvature.B12(s, r) = rSecondVariationsCurvature.B12(r, s);
                }
            }
        }
    }

    void IgaShell3pElement::EquationIdVector(
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
            // rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            // rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        }

        KRATOS_CATCH("")
    }

    // /***********************************************************************************/
    // /***********************************************************************************/
    void IgaShell3pElement::GetDofList(
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
            // rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            // rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            KRATOS_WATCH("GetDofList2");
        }

        KRATOS_WATCH("GetDofListEnd");

        KRATOS_CATCH("")
    }

    // //************************************************************************************
    // //************************************************************************************
    int IgaShell3pElement::Check(const ProcessInfo& rCurrentProcessInfo)
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


