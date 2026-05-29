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
//                   Maram Alkhlaifat
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_6p_element.h"

namespace Kratos
{
  /// Internal variables used for metric transformation
   Shell6pElement::KinematicVariables::KinematicVariables(IndexType Dimension)
    {
        noalias(BaseVector1) = ZeroVector(Dimension); // base vector 1
        noalias(BaseVector2) = ZeroVector(Dimension); // base vector 2
        noalias(NormalVector) = ZeroVector(Dimension);  //base vector 3 normalized
        noalias(NormalVectorTilde) = ZeroVector(Dimension); // not-normalized base vector 3
        DifferentialArea = 1.0; // differential area
    }

    Shell6pElement::ConstitutiveVariables::ConstitutiveVariables(IndexType StrainSize)
    {
        StrainVector       = ZeroVector(StrainSize);
        StressVector       = ZeroVector(StrainSize);
        ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
    }

    ///@name Initialize Functions
    ///@{

    Shell6pElement::Shell6pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    Shell6pElement::Shell6pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Shell6pElement::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const 
    {
        return Kratos::make_intrusive<Shell6pElement>(NewId, pGeom, pProperties);
    }

    Element::Pointer Shell6pElement::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const 
    {
        return Kratos::make_intrusive< Shell6pElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    void Shell6pElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const IndexType number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (mDifferentialAreaVector.size() != number_of_integration_points)
            mDifferentialAreaVector.resize(number_of_integration_points);
        if (mJacobianThicknessDeterminant.size() != number_of_integration_points)
            mJacobianThicknessDeterminant.resize(number_of_integration_points);
        if (mTransformationMatrix.size() != number_of_integration_points)
            mTransformationMatrix.resize(number_of_integration_points);

        KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());

        const double thickness = this->GetProperties().GetValue(THICKNESS);

        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
        {
            CalculateKinematics(point_number, kinematic_variables);

            mDifferentialAreaVector[point_number] = kinematic_variables.DifferentialArea;
            CalculateTransformationFromLocalToGlobalCartesian(kinematic_variables, mTransformationMatrix[point_number]);

            mJacobianThicknessDeterminant[point_number] = ZeroVector(mGaussIntegrationThickness.num_GP_thickness);
            Matrix normal_vector_derivatives = ZeroMatrix(3, 3);
            CalculateNormalVectorDerivatives(point_number, kinematic_variables, normal_vector_derivatives);
            
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                const double zeta = mGaussIntegrationThickness.zeta(Gauss_index);

                // Compute the Jacobian matrix at the initial Gauss point in thickness direction
                Matrix jacobian = ZeroMatrix(3,3);
                double jacobian_det = 0.0;

                for (int i = 0; i < 3; i++) {
                    jacobian(0, i) = kinematic_variables.BaseVector1[i] + (thickness/2) * zeta * normal_vector_derivatives(0, i); 
                    jacobian(1, i) = kinematic_variables.BaseVector2[i] + (thickness/2) * zeta * normal_vector_derivatives(1, i) ; 
                    jacobian(2, i) = kinematic_variables.NormalVector[i] * (thickness/2) ;
                }

                jacobian_det = MathUtils<double>::Det(jacobian);
                mJacobianThicknessDeterminant[point_number][Gauss_index] = jacobian_det;
            }
        }

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    void Shell6pElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const IndexType number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != number_of_integration_points)
            mConstitutiveLawVector.resize(number_of_integration_points);

        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void Shell6pElement::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const IndexType mat_size = GetGeometry().size() * 6;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void Shell6pElement::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const IndexType mat_size = GetGeometry().size() * 6;
        
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        VectorType right_hand_side_vector;
        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void Shell6pElement::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const IndexType number_of_nodes = GetGeometry().size();
        const IndexType mat_size = number_of_nodes * 6;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        MatrixType left_hand_side_matrix;

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void Shell6pElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        // Retrieve geometry, integration points and definition of problem size
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        const IndexType number_of_nodes = r_geometry.size();
        const IndexType mat_size = number_of_nodes * 6;
        const double thickness = this->GetProperties().GetValue(THICKNESS);
        const double density = this->GetProperties().GetValue(DENSITY);
        const double rotational_inertia_factor = thickness * thickness / 12.0;

        if (rMassMatrix.size1() != mat_size)
            rMassMatrix.resize(mat_size, mat_size, false);
        noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

        //// Loop over integration points
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            const double integration_weight = r_integration_points[point_number].Weight();

            //// Loop for zeta (thickness integration points)
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                const double thickness_weight = mGaussIntegrationThickness.integration_weight_thickness(Gauss_index);
                const double jacobian_thickness_det = mJacobianThicknessDeterminant[point_number][Gauss_index];
                const double mass = density * jacobian_thickness_det * thickness_weight * integration_weight;
                
                for (IndexType r = 0; r < number_of_nodes; ++r)
                {
                    for (IndexType s = 0; s < number_of_nodes; ++s)
                    {
                        const double Nrs = r_N(point_number, s) * r_N(point_number, r) * mass;

                        rMassMatrix(6 * s, 6 * r) += Nrs;
                        rMassMatrix(6 * s + 1, 6 * r + 1) += Nrs;
                        rMassMatrix(6 * s + 2, 6 * r + 2) += Nrs;
                        rMassMatrix(6 * s + 3, 6 * r + 3) += Nrs * rotational_inertia_factor;
                        rMassMatrix(6 * s + 4, 6 * r + 4) += Nrs * rotational_inertia_factor;
                        rMassMatrix(6 * s + 5, 6 * r + 5) += Nrs * rotational_inertia_factor;
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void Shell6pElement::CalculateDampingMatrix(
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
            const IndexType number_of_nodes = GetGeometry().size();
            const IndexType mat_size = number_of_nodes * 6;
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

    ///@}
    ///@name Results on Gauss Points
    ///@{

    void Shell6pElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        // Retrieve geometry, integration points and definition of problem size
        const GeometryType& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const IndexType number_of_nodes = r_geometry.size();
        const IndexType mat_size = number_of_nodes * 6;
        const double thickness = this->GetProperties().GetValue(THICKNESS);

        Vector current_displacement = ZeroVector(6*number_of_nodes);
        GetValuesVector(current_displacement,0);

        // Provide a default empty implementation: resize and set zeros
        if (rOutput.size() != r_integration_points.size())
            rOutput.resize(r_integration_points.size());

        std::vector<array_1d<double, 6>> strain_cauchy_cartesian(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_cauchy_cartesian(mGaussIntegrationThickness.num_GP_thickness);

        //// Loop over integration points
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            // Compute Kinematics and Metric components
            KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(point_number, kinematic_variables);

            Matrix normal_vector_derivatives = ZeroMatrix(3, 3);
            CalculateNormalVectorDerivatives(point_number, kinematic_variables, normal_vector_derivatives);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);
            ConstitutiveVariables constitutive_variables(6);

            CalculateConstitutiveVariables(point_number, kinematic_variables, constitutive_variables,
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
            
            // Transform the constitutive matrix to global cartesian coordinates
            constitutive_variables.ConstitutiveMatrix = prod(mTransformationMatrix[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(mTransformationMatrix[point_number]))));
              
            //// Loop for zeta (thickness integration points)
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // Retrieve zeta for the current Gauss point in thickness direction
                const double zeta = mGaussIntegrationThickness.zeta(Gauss_index);

                // Compute the Jacobian matrix at the current Gauss point in thickness direction
                Matrix jacobian = ZeroMatrix(3,3);
                Matrix jacobian_inv = ZeroMatrix(3,3);
                double jacobian_det = 0.0;

                for (int i = 0; i < 3; i++) {
                    jacobian(0, i) = kinematic_variables.BaseVector1[i] + (thickness/2) * zeta * normal_vector_derivatives(0, i); 
                    jacobian(1, i) = kinematic_variables.BaseVector2[i] + (thickness/2) * zeta * normal_vector_derivatives(1, i) ; 
                    jacobian(2, i) = kinematic_variables.NormalVector[i] * (thickness/2) ;
                }

                MathUtils<double>::InvertMatrix(jacobian, jacobian_inv, jacobian_det);

                // Compute the b-operator at the current Gauss point in thickness direction
                Matrix b_operator = ZeroMatrix(6, mat_size);
                CalculateBOperator(point_number, b_operator, zeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);

                strain_cauchy_cartesian[Gauss_index] = prod(b_operator, current_displacement);
                stress_cauchy_cartesian[Gauss_index] = prod(constitutive_variables.ConstitutiveMatrix,strain_cauchy_cartesian[Gauss_index]);
            }
        }
 
        // Extrapolate stress to top/bottom surface using midplane interpolation
        array_1d<double, 6> stress_cauchy_cartesian_midplane;
        stress_cauchy_cartesian_midplane = (stress_cauchy_cartesian[mGaussIntegrationThickness.num_GP_thickness-1] + stress_cauchy_cartesian[0]) / 2.0;

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
        {
            if (rVariable == CAUCHY_STRESS_TOP_XX) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[0] + (stress_cauchy_cartesian[mGaussIntegrationThickness.num_GP_thickness - 1][0]
                    - stress_cauchy_cartesian_midplane[0]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_YY) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[1] + (stress_cauchy_cartesian[mGaussIntegrationThickness.num_GP_thickness - 1][1]
                    - stress_cauchy_cartesian_midplane[1]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_XY) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[3] + (stress_cauchy_cartesian[mGaussIntegrationThickness.num_GP_thickness - 1][3]
                    - stress_cauchy_cartesian_midplane[3]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XX) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[0] + (stress_cauchy_cartesian[0][0] - stress_cauchy_cartesian_midplane[0]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_YY) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[1] + (stress_cauchy_cartesian[0][1] - stress_cauchy_cartesian_midplane[1]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XY) {
                rOutput[point_number] = stress_cauchy_cartesian_midplane[3] + (stress_cauchy_cartesian[0][3] - stress_cauchy_cartesian_midplane[3]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else {
                KRATOS_WATCH("No results for desired variable available in Calculate of Shell6pElement.")
            }
        }   
    }

    void Shell6pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const IndexType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 6 * number_of_control_points)
            rResult.resize(6 * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 6;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();                        
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();                    
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();                    
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X, pos + 3).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y, pos + 4).EquationId(); 
            rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z, pos + 5).EquationId();
        }

        KRATOS_CATCH("")
    };

    void Shell6pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const IndexType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
        }

        KRATOS_CATCH("")
    };

    void Shell6pElement::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const IndexType number_of_control_points = GetGeometry().size();
        const IndexType mat_size = number_of_control_points * 6;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const array_1d<double, 3 >& rotation = GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);
            IndexType index = i * 6;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
            rValues[index + 3] = rotation[0];
            rValues[index + 4] = rotation[1];
            rValues[index + 5] = rotation[2];
        }
    }

    void Shell6pElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
        }
    }

    int Shell6pElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

            // Check whether ConstitutiveLaw is 3D
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 6)
                << "Wrong constitutive law used. Expected strain size is 3 (el id = ) "
                << this->Id() << std::endl;
        }
        return 0;
    }

    ///@}
    ///@name Internal functions
    ///@{

    void Shell6pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        // 0. Retrieve geometry, integration points and definition of problem size
        const GeometryType& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        const IndexType number_of_nodes = r_geometry.size();
        const IndexType mat_size = number_of_nodes * 6;
        const double thickness = this->GetProperties().GetValue(THICKNESS); 

        Vector current_displacement = ZeroVector(6*number_of_nodes);
        GetValuesVector(current_displacement,0);

        //// Loop over integration points
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            // 1. Compute Kinematics and Metric components
            KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(point_number, kinematic_variables);

            Matrix normal_vector_derivatives = ZeroMatrix(3, 3);
            CalculateNormalVectorDerivatives(point_number, kinematic_variables, normal_vector_derivatives);
            
            // 2a. Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);
            ConstitutiveVariables constitutive_variables(6);

            CalculateConstitutiveVariables(point_number, kinematic_variables, constitutive_variables,
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

            // 2b. Transform the constitutive matrix to global cartesian coordinates
            constitutive_variables.ConstitutiveMatrix = prod(mTransformationMatrix[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(mTransformationMatrix[point_number]))));
                
            //// Loop for zeta (thickness integration points)
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // 3.1. Retrieve zeta for the current Gauss point in thickness direction
                const double zeta = mGaussIntegrationThickness.zeta(Gauss_index);

                // 3.2. Compute the Jacobian matrix at the current Gauss point in thickness direction
                Matrix jacobian = ZeroMatrix(3,3);
                Matrix jacobian_inv = ZeroMatrix(3,3);
                double jacobian_det = 0.0;

                for (int i = 0; i < 3; i++) {
                    jacobian(0, i) = kinematic_variables.BaseVector1[i] + (thickness/2) * zeta * normal_vector_derivatives(0, i); 
                    jacobian(1, i) = kinematic_variables.BaseVector2[i] + (thickness/2) * zeta * normal_vector_derivatives(1, i) ; 
                    jacobian(2, i) = kinematic_variables.NormalVector[i] * (thickness/2) ;
                }
                MathUtils<double>::InvertMatrix(jacobian, jacobian_inv, jacobian_det);

                // 3.3. Material stiffness part
                Matrix b_operator = ZeroMatrix(6, mat_size);
                CalculateBOperator(point_number, b_operator, zeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);

                // 3.4. Drilling stiffness part
                Matrix b_drilling = ZeroMatrix(1, mat_size);
                CalculateBDrilling(point_number, b_drilling, jacobian_inv, kinematic_variables);

                // 3.5. Geometric stiffness part
                constitutive_variables.StrainVector = prod(b_operator, current_displacement);
                constitutive_variables.StressVector = prod(constitutive_variables.ConstitutiveMatrix, constitutive_variables.StrainVector);
                Matrix stress_matrix = ZeroMatrix(9,9);
                CalculateStressMatrix(constitutive_variables.StressVector,stress_matrix);

                Matrix b_geometric = ZeroMatrix(9, mat_size);
                CalculateBGeometric(point_number, b_geometric, zeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);
                
                // 3.6. Compute the integration weight 
                double integration_weight = r_integration_points[point_number].Weight() * mJacobianThicknessDeterminant[point_number][Gauss_index]
                                          * mGaussIntegrationThickness.integration_weight_thickness(Gauss_index); 

                // 3.7. Assembly
                if (CalculateStiffnessMatrixFlag == true)
                {
                    // Add the linear stiffness matrix contribution to the element stiffness matrix
                    noalias(rLeftHandSideMatrix) += integration_weight *prod(trans(b_operator), Matrix(prod(constitutive_variables.ConstitutiveMatrix, b_operator)));

                    // Add the drilling stiffness matrix contribution to the element stiffness matrix
                    double drilling_factor = 0.05; //TO DO: user definerd factor to adjust the drilling stiffness contribution
                    noalias(rLeftHandSideMatrix) += drilling_factor * this->GetProperties().GetValue(YOUNG_MODULUS) * integration_weight * prod(trans(b_drilling), Matrix((b_drilling))); 
                    
                    // Add the geometric stiffness matrix contribution to the element stiffness matrix
                    noalias(rLeftHandSideMatrix) += integration_weight *prod(trans(b_geometric), Matrix(prod(stress_matrix, b_geometric)));
                }
                if (CalculateResidualVectorFlag == true)
                {
                    noalias(rRightHandSideVector) -= integration_weight * prod(trans(b_operator), constitutive_variables.StressVector);
                }
            } 
        }
        KRATOS_CATCH("");
    }

    void Shell6pElement::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    ) const
    {
        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);

        rKinematicVariables.BaseVector1 = column(J, 0);
        rKinematicVariables.BaseVector2 = column(J, 1);

        // not-normalized base vector 3
        MathUtils<double>::CrossProduct(rKinematicVariables.NormalVectorTilde, rKinematicVariables.BaseVector1, rKinematicVariables.BaseVector2);

        // differential area DifferentialArea
        rKinematicVariables.DifferentialArea = norm_2(rKinematicVariables.NormalVectorTilde);

        //base vector 3 normalized
        noalias(rKinematicVariables.NormalVector) = rKinematicVariables.NormalVectorTilde / rKinematicVariables.DifferentialArea;
    }

    // Computes transformation matrix from local to global cartesian coordinates
    void Shell6pElement::CalculateTransformationFromLocalToGlobalCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTransformationMatrix
    ) const
    {
        // Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.BaseVector1);
        array_1d<double, 3> local_base_1 = rKinematicVariables.BaseVector1 / l_a1;
        array_1d<double, 3> local_base_3 =  rKinematicVariables.NormalVector;
        array_1d<double, 3> local_base_2;
        MathUtils<double>::CrossProduct(local_base_2, local_base_3, local_base_1);
        
        // Transformation matrix 
        if (rTransformationMatrix.size1() != 6 && rTransformationMatrix.size2() != 6)                                                                 
            rTransformationMatrix.resize(6, 6);
        noalias(rTransformationMatrix) = ZeroMatrix(6, 6);

        for (IndexType i = 0; i < 3; ++i)
        {
            IndexType j = (i + 1) % 3;

            rTransformationMatrix(i, 0) = local_base_1[i] * local_base_1[i];
            rTransformationMatrix(i, 1) = local_base_2[i] * local_base_2[i];
            rTransformationMatrix(i, 2) = local_base_3[i] * local_base_3[i];
            rTransformationMatrix(i, 3) = 2 * local_base_1[i] * local_base_2[i];
            rTransformationMatrix(i, 4) = 2 * local_base_2[i] * local_base_3[i];
            rTransformationMatrix(i, 5) = 2 * local_base_1[i] * local_base_3[i];

            rTransformationMatrix(i + 3, 0) = local_base_1[i] * local_base_1[j];
            rTransformationMatrix(i + 3, 1) = local_base_2[i] * local_base_2[j];
            rTransformationMatrix(i + 3, 2) = local_base_3[i] * local_base_3[j];
            rTransformationMatrix(i + 3, 3) = (local_base_1[i] * local_base_2[j]) + (local_base_2[i] * local_base_1[j]);
            rTransformationMatrix(i + 3, 4) = (local_base_2[i] * local_base_3[j]) + (local_base_3[i] * local_base_2[j]);
            rTransformationMatrix(i + 3, 5) = (local_base_1[i] * local_base_3[j]) + (local_base_3[i] * local_base_1[j]);
        }
    }

    void Shell6pElement::CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        // Constitive Matrices D
        rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //this is an ouput parameter

        // TO DO: for now we support only isotropic linear elastic materials, in the future, this should be computed in the constitutive law.
        // mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        const double poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const double youngs_modulus = this->GetProperties()[YOUNG_MODULUS];
        const double lame_lambda = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
        const double shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
        const double shear_correction_factor = 5.0 / 6.0;
        
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 0) = lame_lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 1) = lame_lambda * poisson_ratio;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 0) = lame_lambda * poisson_ratio;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 1) = lame_lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(3, 3) = lame_lambda * (1.0 - poisson_ratio) / 2.0;
        rThisConstitutiveVariables.ConstitutiveMatrix(4, 4) = shear_modulus * shear_correction_factor;
        rThisConstitutiveVariables.ConstitutiveMatrix(5, 5) = shear_modulus * shear_correction_factor;
    }

    void Shell6pElement::CalculateNormalVectorDerivatives(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        Matrix& rNormalVectorDerivatives) const
    {
        // Get the shape function of second derivatives
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        // Get the area of the element
        double inv_differential_area = 1 / mDifferentialAreaVector[IntegrationPointIndex];
        double inv_differential_area_cube = 1 / std::pow(mDifferentialAreaVector[IntegrationPointIndex], 3);

        // Compute base vector derivatives
        array_1d<double, 3> base_vector1_derivative_11;
        array_1d<double, 3> base_vector1_derivative_12; 
        array_1d<double, 3> base_vector2_derivative_22;

        for (IndexType i = 0; i < GetGeometry().size(); ++i)
        {
            base_vector1_derivative_11[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 0);
            base_vector1_derivative_11[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 0);
            base_vector1_derivative_11[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 0);

            base_vector1_derivative_12[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 1);
            base_vector1_derivative_12[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 1);
            base_vector1_derivative_12[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 1);

            base_vector2_derivative_22[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 2);
            base_vector2_derivative_22[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 2);
            base_vector2_derivative_22[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 2);
        }

        // Compute normal vector derivatives
        array_1d<double, 3> normal_tilde_derivative_1;
        array_1d<double, 3> normal_tilde_derivative_1_term1;
        array_1d<double, 3> normal_tilde_derivative_1_term2;
        array_1d<double, 3> normal_tilde_derivative_2;
        array_1d<double, 3> normal_tilde_derivative_2_term1;
        array_1d<double, 3> normal_tilde_derivative_2_term2;

        MathUtils<double>::CrossProduct(normal_tilde_derivative_1_term1, base_vector1_derivative_11, rKinematicVariables.BaseVector2);
        MathUtils<double>::CrossProduct(normal_tilde_derivative_1_term2, rKinematicVariables.BaseVector1, base_vector1_derivative_12);
        normal_tilde_derivative_1 = normal_tilde_derivative_1_term1 + normal_tilde_derivative_1_term2;

        MathUtils<double>::CrossProduct(normal_tilde_derivative_2_term1, base_vector1_derivative_12, rKinematicVariables.BaseVector2);
        MathUtils<double>::CrossProduct(normal_tilde_derivative_2_term2, rKinematicVariables.BaseVector1, base_vector2_derivative_22);
        normal_tilde_derivative_2 = normal_tilde_derivative_2_term1 + normal_tilde_derivative_2_term2;

        for (IndexType j = 0; j < 3; j++)
        {
            rNormalVectorDerivatives(0, j) = normal_tilde_derivative_1[j] * inv_differential_area - rKinematicVariables.NormalVectorTilde[j] * inner_prod(rKinematicVariables.NormalVectorTilde, normal_tilde_derivative_1) * inv_differential_area_cube;
            rNormalVectorDerivatives(1, j) = normal_tilde_derivative_2[j] * inv_differential_area - rKinematicVariables.NormalVectorTilde[j] * inner_prod(rKinematicVariables.NormalVectorTilde, normal_tilde_derivative_2) * inv_differential_area_cube;
            rNormalVectorDerivatives(2, j) = 0.0;
        }
    }

    void Shell6pElement::CalculateBOperator(
        const IndexType IntegrationPointIndex,
        Matrix& rBOperator,
        const double zeta,
        const Matrix& rJacobianInv,
        const Matrix& rNormalVectorDerivatives,
        const KinematicVariables& rActualKinematic) const
    {
        const IndexType number_of_control_points = GetGeometry().size();
        const IndexType mat_size = number_of_control_points * 6;
        
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const double thickness = this->GetProperties().GetValue(THICKNESS);
        
        // Define shape function derivatives in global coordinates system
        Matrix shape_functions_derivatives_global = ZeroMatrix(number_of_control_points,3);
        Matrix shape_functions_derivatives_local = ZeroMatrix(number_of_control_points,3);
        column(shape_functions_derivatives_local, 0) = column(r_DN_De, 0);
        column(shape_functions_derivatives_local, 1) = column(r_DN_De, 1);
        shape_functions_derivatives_global = trans(prod(rJacobianInv, trans(shape_functions_derivatives_local)));
                                                
        // Normal vector components in global coordinates
        const double normal_x = rActualKinematic.NormalVector[0];
        const double normal_y = rActualKinematic.NormalVector[1];
        const double normal_z = rActualKinematic.NormalVector[2];

        // Normal vector derivatives in global coordinates
        const Matrix normal_vector_derivatives_global = prod(rJacobianInv, rNormalVectorDerivatives);
        const double d_normal_x_dx = normal_vector_derivatives_global(0,0);
        const double d_normal_y_dx = normal_vector_derivatives_global(0,1);
        const double d_normal_z_dx = normal_vector_derivatives_global(0,2);
        const double d_normal_x_dy = normal_vector_derivatives_global(1,0);
        const double d_normal_y_dy = normal_vector_derivatives_global(1,1);
        const double d_normal_z_dy = normal_vector_derivatives_global(1,2);
        const double d_normal_x_dz = normal_vector_derivatives_global(2,0);
        const double d_normal_y_dz = normal_vector_derivatives_global(2,1);
        const double d_normal_z_dz = normal_vector_derivatives_global(2,2);

        // Zeta derivatives in global coordinates  
        const double d_zeta_dx = rJacobianInv(0,2);
        const double d_zeta_dy = rJacobianInv(1,2);
        const double d_zeta_dz = rJacobianInv(2,2);

        if (rBOperator.size1() != 6 || rBOperator.size2() != mat_size)
            rBOperator.resize(6, mat_size);
        noalias(rBOperator) = ZeroMatrix(6, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const IndexType index = i * 6;
            
            //// Displacement DOFs contributions
            rBOperator(0, index)     = shape_functions_derivatives_global(i, 0);
            rBOperator(1, index + 1) = shape_functions_derivatives_global(i, 1);
            rBOperator(2, index + 2) = shape_functions_derivatives_global(i, 2);
            rBOperator(3, index)     = shape_functions_derivatives_global(i, 1);                    
            rBOperator(3, index + 1) = shape_functions_derivatives_global(i, 0);  
            rBOperator(4, index + 1) = shape_functions_derivatives_global(i, 2);
            rBOperator(4, index + 2) = shape_functions_derivatives_global(i, 1);
            rBOperator(5, index)     = shape_functions_derivatives_global(i, 2);                    
            rBOperator(5, index + 2) = shape_functions_derivatives_global(i, 0);

            // Rotation DOFs contributions
            rBOperator(0, index + 4) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (thickness/2);
            rBOperator(0, index + 5) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (thickness/2);

            rBOperator(1, index + 3) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (thickness/2);
            rBOperator(1, index + 5) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (thickness/2);

            rBOperator(2, index + 3) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (thickness/2);
            rBOperator(2, index + 4) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (thickness/2);

            rBOperator(3, index + 3) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (thickness/2);
            rBOperator(3, index + 4) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (thickness/2);
            rBOperator(3, index + 5) = (((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dx + d_zeta_dx * normal_x))) 
                                       -((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))) * (thickness/2);

            rBOperator(4, index + 3) = (((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))  
                                       -((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dz + d_zeta_dz * normal_z)))) * (thickness/2);
            rBOperator(4, index + 4) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (thickness/2);
            rBOperator(4, index + 5) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (thickness/2);

            rBOperator(5, index + 3) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (thickness/2);
            rBOperator(5, index + 4) = (((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dz + d_zeta_dz * normal_z))) 
                                       -((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dx + d_zeta_dx * normal_x)))) * (thickness/2);
            rBOperator(5, index + 5) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (thickness/2);
        }
    }

    void Shell6pElement::CalculateBDrilling(                                                                                         
        const IndexType IntegrationPointIndex,
        Matrix& rBDrilling,
        const Matrix& rJacobianInv,
        const KinematicVariables& rActualKinematic) const
    {
        const IndexType number_of_control_points = GetGeometry().size();
        const IndexType mat_size = number_of_control_points * 6;

        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        // Define shape function derivatives in global coordinates system
        Matrix shape_functions_derivatives_global = ZeroMatrix(number_of_control_points,3);
        Matrix shape_functions_derivatives_local = ZeroMatrix(number_of_control_points,3);
        column(shape_functions_derivatives_local, 0) = column(r_DN_De, 0);
        column(shape_functions_derivatives_local, 1) = column(r_DN_De, 1);
        shape_functions_derivatives_global = trans(prod(rJacobianInv, trans(shape_functions_derivatives_local)));

        if (rBDrilling.size1() != 1|| rBDrilling.size2() != mat_size)                                 
            rBDrilling.resize(1, mat_size);                                                     
        noalias(rBDrilling) = ZeroMatrix(1, mat_size);                                          

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            IndexType index = i * 6;

            rBDrilling(0, index)     = -0.5 * shape_functions_derivatives_global(i, 1);
            rBDrilling(0, index + 1) = 0.5 * shape_functions_derivatives_global(i, 0);
            rBDrilling(0, index + 5) = - r_N(IntegrationPointIndex, i);
        }
    }

    void Shell6pElement::CalculateBGeometric(
        const IndexType IntegrationPointIndex,
        Matrix& rBGeometric,
        const double zeta,
        const Matrix& rJacobianInv,
        const Matrix& rNormalVectorDerivatives,
        const KinematicVariables& rActualKinematic) const
    {
        const IndexType number_of_control_points = GetGeometry().size();
        const IndexType mat_size = number_of_control_points * 6;

        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const double thickness = this->GetProperties().GetValue(THICKNESS); 
        
        // Define shape function derivatives in global coordinates 
        Matrix shape_functions_derivatives_global = ZeroMatrix(number_of_control_points,3);
        Matrix shape_functions_derivatives_local = ZeroMatrix(number_of_control_points,3);
        column(shape_functions_derivatives_local, 0) = column(r_DN_De, 0);
        column(shape_functions_derivatives_local, 1) = column(r_DN_De, 1);
        shape_functions_derivatives_global = trans(prod(rJacobianInv, trans(shape_functions_derivatives_local)));
                                         
        // Normal vector components in global coordinates
        const double normal_x = rActualKinematic.NormalVector[0];
        const double normal_y = rActualKinematic.NormalVector[1];
        const double normal_z = rActualKinematic.NormalVector[2];

        // Normal vector derivatives in global coordinates
        const Matrix normal_vector_derivatives_global = prod(rJacobianInv, rNormalVectorDerivatives);
        const double d_normal_x_dx = normal_vector_derivatives_global(0,0);
        const double d_normal_y_dx = normal_vector_derivatives_global(0,1);
        const double d_normal_z_dx = normal_vector_derivatives_global(0,2);
        const double d_normal_x_dy = normal_vector_derivatives_global(1,0);
        const double d_normal_y_dy = normal_vector_derivatives_global(1,1);
        const double d_normal_z_dy = normal_vector_derivatives_global(1,2);
        const double d_normal_x_dz = normal_vector_derivatives_global(2,0);
        const double d_normal_y_dz = normal_vector_derivatives_global(2,1);
        const double d_normal_z_dz = normal_vector_derivatives_global(2,2);

        // Zeta derivatives in global coordinates  
        const double d_zeta_dx = rJacobianInv(0,2);
        const double d_zeta_dy = rJacobianInv(1,2);
        const double d_zeta_dz = rJacobianInv(2,2);

        if (rBGeometric.size1() != 9 || rBGeometric.size2() != mat_size)
            rBGeometric.resize(9, mat_size);
        noalias(rBGeometric) = ZeroMatrix(9, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const IndexType index = i * 6;

            //// Displacement DOFs contributions
            rBGeometric(0, index)     = shape_functions_derivatives_global(i, 0);
            rBGeometric(1, index)     = shape_functions_derivatives_global(i, 1);
            rBGeometric(2, index)     = shape_functions_derivatives_global(i, 2);
            rBGeometric(3, index + 1) = shape_functions_derivatives_global(i, 0);
            rBGeometric(4, index + 1) = shape_functions_derivatives_global(i, 1);
            rBGeometric(5, index + 1) = shape_functions_derivatives_global(i, 2);
            rBGeometric(6, index + 2) = shape_functions_derivatives_global(i, 0);
            rBGeometric(7, index + 2) = shape_functions_derivatives_global(i, 1);
            rBGeometric(8, index + 2) = shape_functions_derivatives_global(i, 2);

            // Rotation DOFs contributions
            rBGeometric(0, index + 4) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (thickness/2);
            rBGeometric(0, index + 5) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (thickness/2);

            rBGeometric(1, index + 4) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (thickness/2);
            rBGeometric(1, index + 5) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dy + d_zeta_dy * normal_y))) * (thickness/2);

            rBGeometric(2, index + 4) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dz + d_zeta_dz * normal_z))) * (thickness/2);
            rBGeometric(2, index + 5) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (thickness/2);

            rBGeometric(3, index + 3) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (thickness/2); 
            rBGeometric(3, index + 5) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dx + d_zeta_dx * normal_x))) * (thickness/2); 

            rBGeometric(4, index + 3) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (thickness/2); 
            rBGeometric(4, index + 5) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (thickness/2); 

            rBGeometric(5, index + 3) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_z_dz + d_zeta_dz * normal_z))) * (thickness/2); 
            rBGeometric(5, index + 5) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (thickness/2); 

            rBGeometric(6, index + 3) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (thickness/2);  
            rBGeometric(6, index + 4) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dx + d_zeta_dx * normal_x))) * (thickness/2); 

            rBGeometric(7, index + 3) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dy + d_zeta_dy * normal_y))) * (thickness/2);  
            rBGeometric(7, index + 4) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (thickness/2); 

            rBGeometric(8, index + 3) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (thickness/2);  
            rBGeometric(8, index + 4) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (r_N(IntegrationPointIndex, i) * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (thickness/2); 
        }
    }

    void Shell6pElement::CalculateStressMatrix(
        const array_1d<double, 6>& rStressVector,
        Matrix& rStressMatrix
    ) const
    {
        if (rStressMatrix.size1() != 9 || rStressMatrix.size2() != 9)
            rStressMatrix.resize(9, 9);
        noalias(rStressMatrix) = ZeroMatrix(9, 9);

        // Fill the 3x3 symmetric stress block
        rStressMatrix(0, 0) = rStressVector(0);
        rStressMatrix(0, 1) = rStressVector(3);
        rStressMatrix(0, 2) = rStressVector(4);
        rStressMatrix(1, 0) = rStressVector(3);
        rStressMatrix(1, 1) = rStressVector(1);
        rStressMatrix(1, 2) = rStressVector(5);
        rStressMatrix(2, 0) = rStressVector(4);
        rStressMatrix(2, 1) = rStressVector(5);
        rStressMatrix(2, 2) = rStressVector(2);

        // Repeat the same block on the diagonal (blocks at [3,3] and [6,6])
        for (IndexType i = 0; i < 3; ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                rStressMatrix(i + 3, j + 3) = rStressMatrix(i, j);
                rStressMatrix(i + 6, j + 6) = rStressMatrix(i, j);
            }
        }
    }
    ///@}

} // Namespace Kratos


