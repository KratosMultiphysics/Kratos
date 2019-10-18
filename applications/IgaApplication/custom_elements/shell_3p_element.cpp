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

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_3p_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void Shell3pElement::Initialize()
    {
        KRATOS_TRY

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void Shell3pElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_N.size1())
            mConstitutiveLawVector.resize(r_N.size1());


        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Assembly
    ///@{

    void Shell3pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        // definition of problem size
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 3;

        // Compute Kinematics and Metric
        KinematicVariables kinematic_variables(WorkingSpaceDimension());
        CalculateKinematics(kinematic_variables, rCurrentProcessInfo);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables_membrane(3);
        ConstitutiveVariables constitutive_variables_curvature(3);

        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables_membrane,
            constitutive_variables_curvature,
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_PK2);

        // calculate B MATRICES
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        Matrix BCurvature = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);
        CalculateBCurvature(BCurvature, actual_metric);

        // Nonlinear Deformation
        SecondVariations second_variations_strain(mat_size);
        SecondVariations second_variations_curvature(mat_size);
        CalculateSecondVariationStrainCurvature(
            second_variations_strain,
            second_variations_curvature,
            actual_metric);

        integration_weight = this->GetValue(INTEGRATION_WEIGHT) * mInitialMetric.dA * GetProperties()[THICKNESS];

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, constitutive_variables_membrane.D, integration_weight);
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

    ///@}
    ///@name Kinematics
    ///@{

    void Shell3pElement::CalculateKinematics(
        MetricVariables& metric
    )
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        Jacobian(DN_De, metric.J);


        metric.g1[0] = metric.J(0, 0);
        metric.g2[0] = metric.J(0, 1);
        metric.g1[1] = metric.J(1, 0);
        metric.g2[1] = metric.J(1, 1);
        metric.g1[2] = metric.J(2, 0);
        metric.g2[2] = metric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(metric.g3, metric.g1, metric.g2);
        //differential area dA
        metric.dA = norm_2(metric.g3);
        //normal vector _n
        Vector n = metric.g3 / metric.dA;


        //GetCovariantMetric
        metric.gab[0] = pow(metric.g1[0], 2) + pow(metric.g1[1], 2) + pow(metric.g1[2], 2);
        metric.gab[1] = pow(metric.g2[0], 2) + pow(metric.g2[1], 2) + pow(metric.g2[2], 2);
        metric.gab[2] = metric.g1[0] * metric.g2[0] + metric.g1[1] * metric.g2[1] + metric.g1[2] * metric.g2[2];

        CalculateHessian(metric.H, DDN_DDe);

        metric.curvature[0] = metric.H(0, 0)*n[0] + metric.H(1, 0)*n[1] + metric.H(2, 0)*n[2];
        metric.curvature[1] = metric.H(0, 1)*n[0] + metric.H(1, 1)*n[1] + metric.H(2, 1)*n[2];
        metric.curvature[2] = metric.H(0, 2)*n[0] + metric.H(1, 2)*n[1] + metric.H(2, 2)*n[2];


        //contravariant metric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (metric.gab[0] * metric.gab[1] - metric.gab[2] * metric.gab[2]);
        metric.gab_con[0] = invdetGab*metric.gab[1];
        metric.gab_con[2] = -invdetGab*metric.gab[2];
        metric.gab_con[1] = invdetGab*metric.gab[0];


        array_1d<double, 3> g_con_1 = metric.g1*metric.gab_con[0] + metric.g2*metric.gab_con[2];
        array_1d<double, 3> g_con_2 = metric.g1*metric.gab_con[2] + metric.g2*metric.gab_con[1];


        //local cartesian coordinates
        double lg1 = norm_2(metric.g1);
        array_1d<double, 3> e1 = metric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        Matrix mG = ZeroMatrix(2, 2);
        mG(0, 0) = inner_prod(e1, g_con_1);
        mG(0, 1) = inner_prod(e1, g_con_2);
        mG(1, 0) = inner_prod(e2, g_con_1);
        mG(1, 1) = inner_prod(e2, g_con_2);

        metric.Q = ZeroMatrix(3, 3);
        metric.Q(0, 0) = pow(mG(0, 0), 2);
        metric.Q(0, 1) = pow(mG(0, 1), 2);
        metric.Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

        metric.Q(1, 0) = pow(mG(1, 0), 2);
        metric.Q(1, 1) = pow(mG(1, 1), 2);
        metric.Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

        metric.Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
        metric.Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
        metric.Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));

        //Matrix T_G_E = ZeroMatrix(3, 3);
        //Transformation matrix T from contravariant to local cartesian basis
        double eG11 = inner_prod(e1, metric.g1);
        double eG12 = inner_prod(e1, metric.g2);
        double eG21 = inner_prod(e2, metric.g1);
        double eG22 = inner_prod(e2, metric.g2);

        //metric.Q = ZeroMatrix(3, 3);
        //metric.Q(0, 0) = eG11*eG11;
        //metric.Q(0, 1) = eG12*eG12;
        //metric.Q(0, 2) = 2.0*eG11*eG12;
        //metric.Q(1, 0) = eG21*eG21;
        //metric.Q(1, 1) = eG22*eG22;
        //metric.Q(1, 2) = 2.0*eG21*eG22;
        //metric.Q(2, 0) = 2.0*eG11*eG21;
        //metric.Q(2, 1) = 2.0*eG12*eG22;
        //metric.Q(2, 2) = 2.0*eG11*eG22 + eG12*eG21;

        metric.T = ZeroMatrix(3, 3);
        metric.T(0, 0) = eG11*eG11;
        metric.T(0, 1) = eG21*eG21;
        metric.T(0, 2) = 2.0*eG11*eG21;
        metric.T(1, 0) = eG12*eG12;
        metric.T(1, 1) = eG22*eG22;
        metric.T(1, 2) = 2.0*eG12*eG22;
        metric.T(2, 0) = eG11*eG12;
        metric.T(2, 1) = eG21*eG22;
        metric.T(2, 2) = eG11*eG22 + eG12*eG21;
    }
    //************************************************************************************
    //************************************************************************************
    void Shell3pElement::CalculateConstitutiveVariables(
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        Vector strain_vector = ZeroVector(3);
        Vector curvature_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualKinematic.gab, mInitialMetric.gab);
        rThisConstitutiveVariablesMembrane.E = prod(mInitialMetric.Q, strain_vector);

        CalculateCurvature(curvature_vector, rActualMetric.curvature, mInitialMetric.curvature);
        rThisConstitutiveVariablesCurvature.E = prod(mInitialMetric.Q, curvature_vector);

        // Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.D); //this is an ouput parameter
        //rValues.CheckAllParameters();
        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        rThisConstitutiveVariablesCurvature.D = rThisConstitutiveVariablesMembrane.D*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.S = prod(
            trans(rThisConstitutiveVariablesMembrane.D), rThisConstitutiveVariablesMembrane.E);
        rThisConstitutiveVariablesCurvature.S = prod(
            trans(rThisConstitutiveVariablesCurvature.D), rThisConstitutiveVariablesCurvature.E);
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void Shell3pElement::GetValuesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const int index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void Shell3pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void Shell3pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step)
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void Shell3pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void Shell3pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int Shell3pElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            KRATOS_CHECK_VARIABLE_KEY(THICKNESS)
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
                << "THICKNESS not provided for element " << this->Id() << std::endl;

            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;

            // Check constitutive law
            this->GetProperties().GetValue(CONSTITUTIVE_LAW)->Check(this->GetProperties(), r_geometry, rCurrentProcessInfo);
        }

        return 0;
        KRATOS_CATCH("");
    }

    ///@}

} // Namespace Kratos


