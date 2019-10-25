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
        const auto& r_N = r_geometry.ShapeFunctionsValues();

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

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane(3);
            ConstitutiveVariables constitutive_variables_curvature(3);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Matrix BMembrane = ZeroMatrix(3, mat_size);
            Matrix BCurvature = ZeroMatrix(3, mat_size);
            CalculateBMembrane(
                point_number,
                BMembrane,
                kinematic_variables);
            CalculateBCurvature(
                point_number,
                BCurvature,
                kinematic_variables);

            // Nonlinear Deformation
            SecondVariations second_variations_strain(mat_size);
            SecondVariations second_variations_curvature(mat_size);
            CalculateSecondVariationStrainCurvature(
                second_variations_strain,
                second_variations_curvature,
                actual_metric);

            double integration_weight =
                integration_points[point_number].Weight()
                * mdA_vector[point_number]
                * GetProperties()[THICKNESS];

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                //adding membrane contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BMembrane,
                    constitutive_variables_membrane.ConstitutiveMatrix,
                    integration_weight);
                //adding curvature contributions to the stiffness matrix
                CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, constitutive_variables_curvature.ConstitutiveMatrix, integration_weight);

                // adding  non-linear-contribution to Stiffness-Matrix
                CalculateAndAddNonlinearKm(
                    rLeftHandSideMatrix,
                    second_variations_strain,
                    constitutive_variables_membrane.StressVector,
                    integration_weight);

                CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                    second_variations_curvature,
                    constitutive_variables_curvature.StressVector,
                    integration_weight);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.StressVector);
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.StressVector);
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Kinematics
    ///@{

    void Shell3pElement::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    )
    {
        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);

        rKinematicVariables.a1 = column(J, 0);
        rKinematicVariables.a2 = column(J, 1);

        //not-normalized base vector 3
        MathUtils<double>::CrossProduct(rKinematicVariables.a3_tilde, rKinematicVariables.a1, rKinematicVariables.a2);

        //differential area dA
        rKinematicVariables.dA = norm_2(rKinematicVariables.a3_tilde);

        //base vector 3 normalized
        rKinematicVariables.a3 = rKinematicVariables.a3_tilde / rKinematicVariables.dA;


        //GetCovariantMetric
        rKinematicVariables.a_ab_covariant[0] = pow(rKinematicVariables.a1[0], 2) + pow(rKinematicVariables.a1[1], 2) + pow(rKinematicVariables.a1[2], 2);
        rKinematicVariables.a_ab_covariant[1] = pow(rKinematicVariables.a2[0], 2) + pow(rKinematicVariables.a2[1], 2) + pow(rKinematicVariables.a2[2], 2);
        rKinematicVariables.a_ab_covariant[2] = rKinematicVariables.a1[0] * rKinematicVariables.a2[0] + rKinematicVariables.a1[1] * rKinematicVariables.a2[1] + rKinematicVariables.a1[2] * rKinematicVariables.a2[2];

        Matrix H = ZeroMatrix(3, 3);
        CalculateHessian(H, GetGeometry().ShapeFunctionsThirdDerivatives());

        rKinematicVariables.b_ab_covariant[0] = H(0, 0)*n[0] + H(1, 0)*n[1] + H(2, 0)*n[2];
        rKinematicVariables.b_ab_covariant[1] = H(0, 1)*n[0] + H(1, 1)*n[1] + H(2, 1)*n[2];
        rKinematicVariables.b_ab_covariant[2] = H(0, 2)*n[0] + H(1, 2)*n[1] + H(2, 2)*n[2];
    }

    /* Computes the transformation matrix T from the contravariant curvilinear basis to
    *  the local cartesian basis.
    *  ε_curvilinear is defined: [ε_11, ε_22, ε_12]
    *  The transformation matrix T transforms to voigt notation:
    *  ε_local_cartesian = [ε_11, ε_22, 2*ε_12]
    *
    *  The transformation from ε_12_cu to 2*ε_12_ca is included in T.
    */
    void Shell3pElement::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& T
    )
    {
        //Contravariant metric g_ab_con
        double inv_det_g_ab = 1.0 /
            (rKinematicVariables.a_ab_covariant[0] * rKinematicVariables.a_ab_covariant[1]
                - rKinematicVariables.a_ab_covariant[2] * rKinematicVariables.a_ab_covariant[2]);

        array_1d<double, 3> a_ab_contravariant = ZeroVector(3);
        a_ab_contravariant[0] = inv_det_g_ab * rKinematicVariables.a_ab_covariant[1];
        a_ab_contravariant[1] = inv_det_g_ab * rKinematicVariables.a_ab_covariant[0];
        a_ab_contravariant[2] = -inv_det_g_ab * rKinematicVariables.a_ab_covariant[2];

        //Contravariant base vectors
        array_1d<double, 3> a_contravariant_1 = rKinematicVariables.a1*a_ab_contravariant[0] + rKinematicVariables.a2*a_ab_contravariant[2];
        array_1d<double, 3> a_contravariant_2 = rKinematicVariables.a1*a_ab_contravariant[2] + rKinematicVariables.a2*a_ab_contravariant[1];


        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        double l_a_contravariant_2 = norm_2(a_contravariant_2);
        array_1d<double, 3> e2 = a_contravariant_2 / l_a_contravariant_2;

        // e * a_contravariant
        Matrix G = ZeroMatrix(2, 2);
        G(0, 0) = inner_prod(e1, a_contravariant_1);
        G(0, 1) = inner_prod(e1, a_contravariant_2);
        G(1, 0) = inner_prod(e2, a_contravariant_1);
        G(1, 1) = inner_prod(e2, a_contravariant_2);

        //Transformation matrix T
        T = ZeroMatrix(3, 3);
        T(0, 0) = pow(G(0, 0), 2);
        T(0, 1) = pow(G(0, 1), 2);
        T(0, 2) = 2 * G(0, 0) * G(0, 1);

        T(1, 0) = pow(G(1, 0), 2);
        T(1, 1) = pow(G(1, 1), 2);
        T(1, 2) = 2 * G(1, 0) * G(1, 1);

        T(2, 0) = 2 * G(0, 0) * G(1, 0);
        T(2, 1) = 2 * G(0, 1) * G(1, 1);
        T(2, 2) = 2 * (G(0, 0) * G(1, 1) + G(0, 1) * G(1, 0));

        KRATOS_WATCH(T)
    }

    //************************************************************************************
    //************************************************************************************
    void Shell3pElement::CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
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

        array_1d<double, 3> strain_vector = 0.5 * (rActualKinematic.a_ab_covariant - mA_ab_covariant_vector[IntegrationPointIndex]);
        rThisConstitutiveVariablesMembrane.StrainVector = prod(mT_vector[IntegrationPointIndex], strain_vector);

        array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - mB_ab_covariant_vector[IntegrationPointIndex];
        rThisConstitutiveVariablesCurvature.StrainVector = prod(mT_vector[IntegrationPointIndex], curvature_vector);

        // Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        double thickness = this->GetProperties().GetValue(THICKNESS);
        rThisConstitutiveVariablesCurvature.ConstitutiveMatrix = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix*(pow(thickness, 2) / 12);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.StressVector = prod(
            trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
        rThisConstitutiveVariablesCurvature.StressVector = prod(
            trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);
    }

    void Shell3pElement::CalculateBMembrane(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        KinematicVariables& rActualKinematic)
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size1() != mat_size || rB.size2() != mat_size)
            rB.resize(mat_size, mat_size);
        rB = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 3;
            IndexType dirr = r % 3;

            array_1d<double, 3> dE_curvilinear = ZeroVector(3);
            // strain
            dE_curvilinear[0] = r_DN_De(kr, 0)*rActualKinematic.a1(dirr);
            dE_curvilinear[1] = r_DN_De(kr, 1)*rActualKinematic.a2(dirr);
            dE_curvilinear[2] = 0.5*(r_DN_De(kr, 0)*rActualKinematic.a2(dirr) + rActualKinematic.a1(dirr)*r_DN_De(kr, 1));

            rB(0, r) = mT_vector[IntegrationPointIndex](0, 0)*dE_curvilinear[0] + mT_vector[IntegrationPointIndex](0, 1)*dE_curvilinear[1] + mT_vector[IntegrationPointIndex](0, 2)*dE_curvilinear[2];
            rB(1, r) = mT_vector[IntegrationPointIndex](1, 0)*dE_curvilinear[0] + mT_vector[IntegrationPointIndex](1, 1)*dE_curvilinear[1] + mT_vector[IntegrationPointIndex](1, 2)*dE_curvilinear[2];
            rB(2, r) = mT_vector[IntegrationPointIndex](2, 0)*dE_curvilinear[0] + mT_vector[IntegrationPointIndex](2, 1)*dE_curvilinear[1] + mT_vector[IntegrationPointIndex](2, 2)*dE_curvilinear[2];
        }
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

        //KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            //KRATOS_CHECK_VARIABLE_KEY(THICKNESS)
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
                << "THICKNESS not provided for element " << this->Id() << std::endl;

            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;

            // Check constitutive law
            //this->GetProperties().GetValue(CONSTITUTIVE_LAW)->Check(this->GetProperties(), r_geometry, rCurrentProcessInfo);
        }

        return 0;
        KRATOS_CATCH("");
    }

    void Shell3pElement::CalculateHessian(
        Matrix& Hessian,
        const Matrix& DDN_DDe)
    {
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;
        Hessian.resize(working_space_dimension, working_space_dimension);
        Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

        for (size_t k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian(0, 0) += DDN_DDe(k, 0)*coords[0];
            Hessian(0, 1) += DDN_DDe(k, 1)*coords[0];
            Hessian(0, 2) += DDN_DDe(k, 2)*coords[0];

            Hessian(1, 0) += DDN_DDe(k, 0)*coords[1];
            Hessian(1, 1) += DDN_DDe(k, 1)*coords[1];
            Hessian(1, 2) += DDN_DDe(k, 2)*coords[1];

            Hessian(2, 0) += DDN_DDe(k, 0)*coords[2];
            Hessian(2, 1) += DDN_DDe(k, 1)*coords[2];
            Hessian(2, 2) += DDN_DDe(k, 2)*coords[2];
        }
    }

    ///@}

} // Namespace Kratos


