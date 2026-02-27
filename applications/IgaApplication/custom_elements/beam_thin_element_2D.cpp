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
#include "custom_elements/beam_thin_element_2D.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void BeamThinElement2D::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const std::size_t r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_11_covariant_vector.size() != r_number_of_integration_points)
            m_A_11_covariant_vector.resize(r_number_of_integration_points);
        if (m_B_11_covariant_vector.size() != r_number_of_integration_points)
            m_B_11_covariant_vector.resize(r_number_of_integration_points);
        if (m_dL_vector.size() != r_number_of_integration_points)
            m_dL_vector.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateKinematics(point_number, kinematic_variables);
            
            m_A_11_covariant_vector[point_number] = kinematic_variables.a_11_covariant;
            m_B_11_covariant_vector[point_number] = kinematic_variables.b_11_covariant;

            m_dL_vector[point_number] = kinematic_variables.dL;
        }

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void BeamThinElement2D::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        const std::size_t r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Assembly
    ///@{

    void BeamThinElement2D::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) const
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        // definition of problem size
        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 2;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(point_number, kinematic_variables);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables_membrane;
            ConstitutiveVariables constitutive_variables_curvature;
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables_membrane,
                constitutive_variables_curvature,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            // calculate B MATRICES
            Vector BMembrane = ZeroVector(mat_size);
            Vector BCurvature = ZeroVector(mat_size);

            CalculateBMembrane(point_number, BMembrane, kinematic_variables);
            CalculateBCurvature(point_number, BCurvature, kinematic_variables);

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

    void BeamThinElement2D::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    ) const
    {
        const GeometryType& r_geometry = GetGeometry();
        const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);

        for (IndexType i = 0; i < r_geometry.size(); i++)
        {
            rKinematicVariables.a1[0] += r_DN_De(i, 0) * r_geometry[i].X();
            rKinematicVariables.a1[1] += r_DN_De(i, 0) * r_geometry[i].Y();
        }

        //not-normalized base vector 2
        rKinematicVariables.a2_tilde[0] = -rKinematicVariables.a1[1];
        rKinematicVariables.a2_tilde[1] = rKinematicVariables.a1[0];

        //differential length dL
        rKinematicVariables.dL = norm_2(rKinematicVariables.a1);

        //base vector 3 normalized
        noalias(rKinematicVariables.a2) = rKinematicVariables.a2_tilde / rKinematicVariables.dL;

        //GetCovariantMetric
        rKinematicVariables.a_11_covariant = std::pow(rKinematicVariables.a1[0], 2) + std::pow(rKinematicVariables.a1[1], 2);

        array_1d<double, 2> H = ZeroVector(2);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        rKinematicVariables.b_11_covariant = H[0] * rKinematicVariables.a2[0] + H[1] * rKinematicVariables.a2[1];
    }

    void BeamThinElement2D::CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
    {
        double strain = 0.5 * (rActualKinematic.a_11_covariant - m_A_11_covariant_vector[IntegrationPointIndex]);
        double curvature = rActualKinematic.b_11_covariant - m_B_11_covariant_vector[IntegrationPointIndex];

        rThisConstitutiveVariablesMembrane.StrainValue = strain;
        rThisConstitutiveVariablesCurvature.StrainValue = curvature;

        double thickness = this->GetProperties().GetValue(THICKNESS);
        double cross_area = this->GetProperties().GetValue(CROSS_AREA);
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);

        rThisConstitutiveVariablesMembrane.ConstitutiveValue = E * cross_area / std::pow(rActualKinematic.dL, 4);
        rThisConstitutiveVariablesCurvature.ConstitutiveValue = E * cross_area * std::pow(thickness, 2)/ (12.0 * std::pow(rActualKinematic.dL, 4));

        // //Local Cartesian Forces and Moments
        rThisConstitutiveVariablesMembrane.StressValue = rThisConstitutiveVariablesMembrane.ConstitutiveValue * rThisConstitutiveVariablesMembrane.StrainValue;
        rThisConstitutiveVariablesCurvature.StressValue = rThisConstitutiveVariablesCurvature.ConstitutiveValue * rThisConstitutiveVariablesCurvature.StrainValue;
    }

    void BeamThinElement2D::CalculateBMembrane(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 2;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        if (rB.size() != mat_size)
            rB.resize(mat_size);
        noalias(rB) = ZeroVector(mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            // strain
            rB[r] = r_DN_De(kr, 0) * rActualKinematic.a1(dirr);
        }
    }

    void BeamThinElement2D::CalculateBCurvature(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();

        const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 2;

        if (rB.size() != mat_size)
            rB.resize(mat_size);
        noalias(rB) = ZeroVector(mat_size);

        array_1d<double, 2> H = ZeroVector(2);
        CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // compute coeficients needed for the curvature
        double d1 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[1];
        double d2 = (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[1] * rActualKinematic.a1[1] - 1.0;
        double d3 = 1.0 - (1/(rActualKinematic.dL * rActualKinematic.dL)) * rActualKinematic.a1[0] * rActualKinematic.a1[0];

        double alpha1 = (H[0] * d1 + H[1] * d3)/rActualKinematic.dL;
        double alpha2 = (H[0] * d2 - H[1] * d1)/rActualKinematic.dL;

        for (IndexType i = 0; i < number_of_control_points; i++)
        {
            IndexType index = 2 * i;
            
            rB[index] = - alpha1 * r_DN_De(i, 0) + rActualKinematic.a1[1] * r_DDN_DDe(i, 0) / rActualKinematic.dL;
            rB[index + 1] = - alpha2 * r_DN_De(i, 0) - rActualKinematic.a1[0] * r_DDN_DDe(i, 0) / rActualKinematic.dL;
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Stiffness matrix assembly
    ///@{

    inline void BeamThinElement2D::CalculateAndAddKm(
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

    void BeamThinElement2D::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 3;

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

    void BeamThinElement2D::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const IndexType index = i * 2;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
        }
    }

    void BeamThinElement2D::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const IndexType index = i * 2;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
        }
    }

    void BeamThinElement2D::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const std::size_t number_of_control_points = GetGeometry().size();

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 2;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
        }

        KRATOS_CATCH("")
    };

    void BeamThinElement2D::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const std::size_t number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int BeamThinElement2D::Check(const ProcessInfo& rCurrentProcessInfo) const
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

    void BeamThinElement2D::CalculateHessian(
        array_1d<double, 2>& Hessian,
        const Matrix& rDDN_DDe) const
    {
        const std::size_t number_of_points = GetGeometry().size();

        for (IndexType k = 0; k < number_of_points; k++)
        {
            const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

            Hessian[0] += rDDN_DDe(k, 0)*coords[0];
            Hessian[1] += rDDN_DDe(k, 0)*coords[1];
        }
    }

    ///@}

} // Namespace Kratos


