//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/beam_thick_element_2D.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void BeamThickElement2D::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_11_covariant_vector.size() != r_number_of_integration_points)
            m_A_11_covariant_vector.resize(r_number_of_integration_points);
        if (m_B_11_covariant_vector.size() != r_number_of_integration_points)
            m_B_11_covariant_vector.resize(r_number_of_integration_points);
        if (m_dL_vector.size() != r_number_of_integration_points)
            m_dL_vector.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateKinematics(
                point_number,
                kinematic_variables);
            
            m_A_11_covariant_vector[point_number] = kinematic_variables.a_11_covariant;
            m_B_11_covariant_vector[point_number] = kinematic_variables.b_11_covariant;

            m_dL_vector[point_number] = kinematic_variables.dL;
        }

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void BeamThickElement2D::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void BeamThickElement2D::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
        }
    }

    ///@}
    ///@name Assembly
    ///@{

    void BeamThickElement2D::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane;
            ConstitutiveVariables constitutive_variables_curvature;
            ConstitutiveVariables constitutive_variables_shear;
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_variables_shear,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Vector BMembrane = ZeroVector(mat_size);
            Vector BCurvature = ZeroVector(mat_size);
            Vector BShear = ZeroVector(mat_size);
            CalculateBMembrane(
                point_number,
                BMembrane,
                kinematic_variables);
            CalculateBCurvature(
                point_number,
                BCurvature,
                kinematic_variables);
            CalculateBShear(
                point_number,
                BShear,
                kinematic_variables);

            double integration_weight =
                r_integration_points[point_number].Weight()
                * m_dL_vector[point_number];

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                //adding membrane contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BMembrane,
                    constitutive_variables_membrane.ConstitutiveValue,
                    integration_weight);
                //adding curvature contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BCurvature,
                    constitutive_variables_curvature.ConstitutiveValue,
                    integration_weight);
                //adding shear contributions to the stiffness matrix
                CalculateAndAddKm(
                    rLeftHandSideMatrix,
                    BShear,
                    constitutive_variables_shear.ConstitutiveValue,
                    integration_weight);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= Weight*IntForce
                noalias(rRightHandSideVector) -= integration_weight * trans(BMembrane) * constitutive_variables_membrane.StressValue;
                noalias(rRightHandSideVector) -= integration_weight * trans(BCurvature) * constitutive_variables_curvature.StressValue;
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Kinematics
    ///@{

    void BeamThickElement2D::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);

        for (IndexType i = 0; i < r_geometry.size(); i++)
        {
            rKinematicVariables.a1[0] += r_DN_De(i, 0) * r_geometry[i].X();
            rKinematicVariables.a1[1] += r_DN_De(i, 0) * r_geometry[i].Y();
            rKinematicVariables.a1[2] += r_DN_De(i, 0) * r_geometry[i].Z();
            //Compute the cross sectional rotation scalar
            rKinematicVariables.beta += r_N(IntegrationPointIndex, i) * r_geometry[i].FastGetSolutionStepValue(CROSS_SECTIONAL_ROTATION);
            rKinematicVariables.beta_deriv += r_DN_De(i, 0) * r_geometry[i].FastGetSolutionStepValue(CROSS_SECTIONAL_ROTATION);
        }

        //not-normalized base vector 2
        rKinematicVariables.a2_tilde[0] = -rKinematicVariables.a1[1];
        rKinematicVariables.a2_tilde[1] = rKinematicVariables.a1[0];
        rKinematicVariables.a2_tilde[2] = rKinematicVariables.a1[2];

        //differential length dL
        rKinematicVariables.dL = norm_2(rKinematicVariables.a1);

        //base vector 3 normalized
        noalias(rKinematicVariables.a2) = rKinematicVariables.a2_tilde / rKinematicVariables.dL;

        //GetCovariantMetric
        rKinematicVariables.a_11_covariant = pow(rKinematicVariables.a1[0], 2) + pow(rKinematicVariables.a1[1], 2) + pow(rKinematicVariables.a1[2], 2);

        array_1d<double, 3> H = ZeroVector(3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        rKinematicVariables.b_11_covariant = H[0] * rKinematicVariables.a2[0] + H[1] * rKinematicVariables.a2[1] + H[2] * rKinematicVariables.a2[2];
    }

    void BeamThickElement2D::CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveVariables& rThisConstitutiveVariablesShear,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
    {
        double strain = 0.5 * (rActualKinematic.a_11_covariant - m_A_11_covariant_vector[IntegrationPointIndex]);
        double curvature = rActualKinematic.b_11_covariant - m_B_11_covariant_vector[IntegrationPointIndex]; //TO DO
        double shear = 0.0; //TO DO

        double thickness = this->GetProperties().GetValue(THICKNESS);
        double cross_area = this->GetProperties().GetValue(CROSS_AREA);
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);
        double nu = this->GetProperties().GetValue(POISSON_RATIO);

        rThisConstitutiveVariablesMembrane.ConstitutiveValue = E * cross_area / pow(rActualKinematic.dL, 4);
        rThisConstitutiveVariablesCurvature.ConstitutiveValue = E * cross_area * pow(thickness, 2)/ (12.0 * pow(rActualKinematic.dL, 4));
        rThisConstitutiveVariablesShear.ConstitutiveValue = (E / (2 * (1 + nu)))* (5 / 6) * cross_area / (pow(rActualKinematic.dL, 4));

        // //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.StressValue = rThisConstitutiveVariablesMembrane.ConstitutiveValue * rThisConstitutiveVariablesMembrane.StrainValue;
        rThisConstitutiveVariablesCurvature.StressValue = rThisConstitutiveVariablesCurvature.ConstitutiveValue * rThisConstitutiveVariablesCurvature.StrainValue;
        rThisConstitutiveVariablesShear.StressValue = rThisConstitutiveVariablesShear.ConstitutiveValue * rThisConstitutiveVariablesShear.StrainValue;
    }

    void BeamThickElement2D::CalculateBMembrane(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size() != mat_size)
            rB.resize(mat_size);
        noalias(rB) = ZeroVector(mat_size);

        for (IndexType i = 0; i < number_of_control_points; i++)
        {
            IndexType index = 3 * i;
            
            rB[index] = r_DN_De(i, 0) * rActualKinematic.a1[0];
            rB[index + 1] = r_DN_De(i, 0) * rActualKinematic.a1[1];
        }
    }

    void BeamThickElement2D::CalculateBCurvature(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rB.size() != mat_size)
            rB.resize(mat_size);
        noalias(rB) = ZeroVector(mat_size);

        array_1d<double, 3> H = ZeroVector(3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // compute coeficients needed for the curvature
        double d1 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[1];
        double d2 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[1] * rActualKinematic.a1[1] - 1.0;
        double d3 = 1.0 - (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[0];

        double alpha1 = (H[0] * d1 + H[1] * d3)/rActualKinematic.dL; //A_2_1[0]
        double alpha2 = (H[0] * d2 - H[1] * d1)/rActualKinematic.dL; //A_2_1[1]

        double e1 = 3.0 * H[0] * rActualKinematic.a1[0];
        double e2 = H[0] * rActualKinematic.a1[1] + 2.0 * H[1] * rActualKinematic.a1[0];
        double e3 = H[1] * rActualKinematic.a1[0] + 2.0 * H[0] * rActualKinematic.a1[1];
        double e4 = 3.0 * H[1] * rActualKinematic.a1[1];

        double f1 = rActualKinematic.a1[0] * rActualKinematic.a1[0];
        double f2 = rActualKinematic.a1[0] * rActualKinematic.a1[1];
        double f3 = rActualKinematic.a1[1] * rActualKinematic.a1[1];

        double g1 = -3.0 * f1 * (H[0] * rActualKinematic.a1[0] + H[1] * rActualKinematic.a1[1]);
        double g2 = -3.0 * f2 * (H[0] * rActualKinematic.a1[0] + H[1] * rActualKinematic.a1[1]);
        double g3 = -3.0 * f3 * (H[0] * rActualKinematic.a1[0] + H[1] * rActualKinematic.a1[1]);

        for (IndexType i = 0; i < number_of_control_points; i++)
        {
            IndexType index = 3 * i;
            
            //first part + second part subpart 2 + second part subpart 4 (a,b1,b2,c)
            rB[index] = (alpha1 - rActualKinematic.beta * alpha2 - rActualKinematic.beta_deriv * rActualKinematic.a2[1]) * r_DN_De(i, 0)
                      + ((alpha2 * rActualKinematic.beta_deriv) * rActualKinematic.a1[0] * r_DN_De(i, 0)) 
                      + (r_DDN_DDe(i, 0) * (rActualKinematic.a1[0] * rActualKinematic.beta + rActualKinematic.a1[1]) / rActualKinematic.dL)
                      + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (e3 + e1 * rActualKinematic.beta) + rActualKinematic.a1[1] * (-e1 + e3 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3)); 
                      + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (f2 + f1 * rActualKinematic.beta) + rActualKinematic.a1[1] * (-f1 + f2 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3)); 
                      + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (g2 + g1 * rActualKinematic.beta) + rActualKinematic.a1[1] * (-g1 + g2 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3)); 
            rB[index + 1] = (alpha2 + rActualKinematic.beta * alpha1 + rActualKinematic.beta_deriv * rActualKinematic.a2[0]) * r_DN_De(i, 0)
                          + ((alpha1 * rActualKinematic.beta_deriv) * rActualKinematic.a1[1] * r_DN_De(i, 0))
                          + (r_DDN_DDe(i, 0) * (rActualKinematic.a1[1] * rActualKinematic.beta - rActualKinematic.a1[0]) / rActualKinematic.dL)
                          + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (e4 - e2 * rActualKinematic.beta) + rActualKinematic.a1[1] * (e2 + e4 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3))
                          + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (f3 - f2 * rActualKinematic.beta) + rActualKinematic.a1[1] * (f2 + f3 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3))
                          + (r_DN_De(i, 0) * (rActualKinematic.a1[0] * (g3 - g2 * rActualKinematic.beta) + rActualKinematic.a1[1] * (g2 + g3 * rActualKinematic.beta)) / pow(rActualKinematic.dL, 3));
            rB[index + 2] = r_DN_De(i, 0) * rActualKinematic.a1[0] * rActualKinematic.a1[1] * 2.0 / rActualKinematic.dL 
                          + r_N(IntegrationPointIndex, i) * (-rActualKinematic.a1[0] * alpha2 + rActualKinematic.a1[1] * alpha1); //second part subpart 1 + second part subpart 3
        }

        KRATOS_CATCH("")
    }

    void BeamThickElement2D::CalculateBShear(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size() != mat_size)
            rB.resize(mat_size);
        noalias(rB) = ZeroVector(mat_size);

        array_1d<double, 3> H = ZeroVector(3);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // compute coeficients needed 
        double d1 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[1];
        double d2 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[1] * rActualKinematic.a1[1] - 1.0;
        double d3 = 1.0 - (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[0];

        double alpha1 = (H[0] * d1 + H[1] * d3)/rActualKinematic.dL;
        double alpha2 = (H[0] * d2 - H[1] * d1)/rActualKinematic.dL;

        for (IndexType i = 0; i < number_of_control_points; i++)
        {
            IndexType index = 3 * i;
            
            rB[index] = (r_DN_De(i, 0) * rActualKinematic.a2[0]) + ((-alpha1 + alpha2 * rActualKinematic.beta) * rActualKinematic.a1[0] * r_DN_De(i, 0));
                        // first part + second part sub part 2
            rB[index + 1] = (r_DN_De(i, 0) * rActualKinematic.a2[1]) + ((-alpha1 - alpha2 * rActualKinematic.beta) * rActualKinematic.a1[1] * r_DN_De(i, 0));
            rB[index + 2] = r_N(IntegrationPointIndex, i) * rActualKinematic.a1[0] * rActualKinematic.a1[1] * 2.0 / rActualKinematic.dL; //second part sub part 1
        }
    }

    ///@}
    ///@name Stiffness matrix assembly
    ///@{

    inline void BeamThickElement2D::CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Vector& rB,
        const double& rD,
        const double IntegrationWeight
    ) const
    {
        noalias(rLeftHandSideMatrix) += IntegrationWeight * outer_prod(trans(rB), rD * rB);
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void BeamThickElement2D::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void BeamThickElement2D::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const IndexType index = i * 2;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
        }
    }

    void BeamThickElement2D::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const IndexType index = i * 2;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
        }
    }

    void BeamThickElement2D::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);
    
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(CROSS_SECTIONAL_ROTATION).EquationId();
        }

        KRATOS_CATCH("")
    };

    void BeamThickElement2D::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(CROSS_SECTIONAL_ROTATION));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int BeamThickElement2D::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Verify that the constitutive law has the correct dimension
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
                << "THICKNESS not provided for element " << this->Id() << std::endl;
        }

        return 0;
    }

    void BeamThickElement2D::CalculateHessian(
        array_1d<double, 3>& Hessian,
        const Matrix& rDDN_DDe) const
    {
        const SizeType number_of_points = GetGeometry().size();
        const SizeType working_space_dimension = 3;

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian[0] += rDDN_DDe(k, 0)*coords[0];
            Hessian[1] += rDDN_DDe(k, 0)*coords[1];
            Hessian[2] += rDDN_DDe(k, 0)*coords[2];
        }
    }

    ///@}

} // Namespace Kratos


