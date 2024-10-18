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
//                   Pia Halbich
//                   Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/truss_embedded_edge_element.h"



namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void TrussEmbeddedEdgeElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    array_1d<double, 3> TrussEmbeddedEdgeElement::GetActualBaseVector(const Matrix& r_DN_De, const ConfigurationType& rConfiguration)
    {
        const GeometryType& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        array_1d<double, 3> tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);

        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        array_1d<double, 3> actual_base_vector = ZeroVector(dimension);

        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        if (rConfiguration==ConfigurationType::Current) GetValuesVector(current_displacement);

        //basis vectors g1 and g2
        Vector g1 = ZeroVector(dimension);
        Vector g2 = ZeroVector(dimension);

        for (SizeType i=0;i<number_of_nodes;++i){
            g1[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * r_DN_De(i, 0);
            g1[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * r_DN_De(i, 0);
            g1[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * r_DN_De(i, 0);

            g2[0] += (GetGeometry().GetPoint( i ).X0()+current_displacement[i*dimension]) * r_DN_De(i, 1);
            g2[1] += (GetGeometry().GetPoint( i ).Y0()+current_displacement[(i*dimension)+1]) * r_DN_De(i, 1);
            g2[2] += (GetGeometry().GetPoint( i ).Z0()+current_displacement[(i*dimension)+2]) * r_DN_De(i, 1);
        }

        actual_base_vector = g1 * tangents[0] + g2 * tangents[1];

        return actual_base_vector;
    }

    ///@}
    ///@name Assembly
    ///@{

    void TrussEmbeddedEdgeElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (mReferenceBaseVector.size() != r_number_of_integration_points)
            mReferenceBaseVector.resize(r_number_of_integration_points);

        //get properties
        array_1d<double, 3> tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);
        const double E = GetProperties()[YOUNG_MODULUS];
        const double A = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
        {
            // get integration data
            const double& integration_weight = r_integration_points[point_number].Weight();
            const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

            mReferenceBaseVector[point_number] = GetActualBaseVector(r_DN_De, ConfigurationType::Reference);
            const double reference_a = norm_2(mReferenceBaseVector[point_number]);

            // compute base vectors
            const array_1d<double, 3> actual_base_vector = GetActualBaseVector(r_DN_De, ConfigurationType::Current);

            // green-lagrange strain
            const double e11_membrane = 0.5 * (inner_prod(actual_base_vector, actual_base_vector) - inner_prod(mReferenceBaseVector[point_number], mReferenceBaseVector[point_number]));

            // normal force reference_aa
            const double s11_membrane = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

            for (IndexType r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                IndexType kr = r / 3;
                IndexType dirr = r % 3;

                const double epsilon_var_r = actual_base_vector[dirr] *
                    (r_DN_De(kr, 0) * tangents[0]
                    + r_DN_De(kr, 1) * tangents[1]) / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                if (CalculateStiffnessMatrixFlag) {
                    for (IndexType s = 0; s < mat_size; s++)
                    {
                        // local node number ks and dof direction dirs
                        IndexType ks = s / 3;
                        IndexType dirs = s % 3;

                        const double epsilon_var_s =
                            actual_base_vector[dirs] *
                            (r_DN_De(ks, 0) * tangents[0]
                            + r_DN_De(ks, 1) * tangents[1])
                            / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                        rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r * epsilon_var_s * reference_a * integration_weight;

                        if (dirr == dirs) {
                            const double epsilon_var_rs =
                            (r_DN_De(kr, 0) * tangents[0] + r_DN_De(kr, 1) * tangents[1]) *
                            (r_DN_De(ks, 0) * tangents[0] + r_DN_De(ks, 1) * tangents[1]) /inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                            rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs * reference_a * integration_weight;
                        }
                    }
                }
                if (CalculateResidualVectorFlag) {
                    rRightHandSideVector[r] = -s11_membrane * epsilon_var_r * reference_a * integration_weight;
                }
            }
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Implicit
    ///@{

    void TrussEmbeddedEdgeElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        // Rayleigh Damping Matrix: alpha*M + beta*K

        // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if (GetProperties().Has(RAYLEIGH_ALPHA))
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

        double beta = 0.0;
        if (GetProperties().Has(RAYLEIGH_BETA))
            beta = GetProperties()[RAYLEIGH_BETA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        // 2.-Calculate StiffnessMatrix and MassMatrix:
        if (std::abs(alpha) < 1E-12 && std::abs(beta) < 1E-12) {
            // no damping specified, only setting the matrix to zero
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType mat_size = number_of_nodes * 3;
            if (rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size) {
                rDampingMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);
        } else if (std::abs(alpha) > 1E-12 && std::abs(beta) < 1E-12) {
            // damping only required with the mass matrix
            CalculateMassMatrix(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= alpha;
        } else if (std::abs(alpha) < 1E-12 && std::abs(beta) > 1E-12) {
            // damping only required with the stiffness matrix
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;
        } else {
            // damping with both mass matrix and stiffness matrix required
            CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
            rDampingMatrix *= beta;

            Matrix mass_matrix;
            CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
            noalias(rDampingMatrix) += alpha  * mass_matrix;
        }

        KRATOS_CATCH("")
    }

    void TrussEmbeddedEdgeElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = number_of_nodes * 3;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        // Shape function values for all integration points
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            double integration_weight = r_integration_points[point_number].Weight();

            double area = this->GetProperties().GetValue(CROSS_AREA);
            double density = this->GetProperties().GetValue(DENSITY);

            double mass = area * density * norm_2(mReferenceBaseVector[point_number]) * integration_weight;

            if (rMassMatrix.size1() != mat_size)
                rMassMatrix.resize(mat_size, mat_size, false);

            rMassMatrix = ZeroMatrix(mat_size, mat_size);

            for (unsigned int r = 0; r<number_of_nodes; r++)
            {
                for (unsigned int s = 0; s<number_of_nodes; s++)
                {
                    rMassMatrix(3 * s, 3 * r) = r_N(point_number, s)*r_N(point_number, r)*mass;
                    rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                    rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
                }
            }
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Postprocessing
    ///@{

    void TrussEmbeddedEdgeElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rValues.size() != r_integration_points.size())
        {
            rValues.resize(r_integration_points.size());
        }

        //get properties
        const double E = GetProperties()[YOUNG_MODULUS];
        const double A = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        if (rVariable==FORCE_PK2_1D || rVariable==FORCE_CAUCHY_1D)
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                // get integration data
                const Matrix& r_DN_De   = r_geometry.ShapeFunctionLocalGradient(point_number);

                const array_1d<double, 3> actual_base_vector = GetActualBaseVector(r_DN_De, ConfigurationType::Current);
                const double reference_a = norm_2(mReferenceBaseVector[point_number]);
                const double actual_a = norm_2(actual_base_vector);

                const double actual_aa = actual_a * actual_a;
                const double reference_aa = reference_a * reference_a;

                // green-lagrange strain
                const double e11_membrane = 0.5 * (actual_aa - reference_aa);

                // normal force reference_aa
                double normal_force = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector[point_number],mReferenceBaseVector[point_number]);

                if (rVariable==FORCE_PK2_1D)
                {
                    rValues[point_number] = normal_force;
                }

                if (rVariable==FORCE_CAUCHY_1D)
                {
                    rValues[point_number] = normal_force * actual_a / reference_a;
                }
            }
        }
        else
        {
            for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
            {
                rValues[point_number] = 0.0;
            }
        }
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void TrussEmbeddedEdgeElement::GetValuesVector(
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
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void TrussEmbeddedEdgeElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const IndexType index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void TrussEmbeddedEdgeElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const IndexType index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void TrussEmbeddedEdgeElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void TrussEmbeddedEdgeElement::GetDofList(
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
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    ///@}
    ///@name Check
    ///@{

    int TrussEmbeddedEdgeElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        // Verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }
        else
        {
            // Check strain size
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
                << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;
        }

        return 0;
    }

    ///@}

} // Namespace Kratos
