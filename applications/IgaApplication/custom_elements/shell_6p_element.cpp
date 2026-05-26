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
   Shell6pElement::KinematicVariables::KinematicVariables(std::size_t Dimension)
    {
        // covariant metric
        noalias(MetricCovariant) = ZeroVector(Dimension);
        noalias(CurvatureCovariant) = ZeroVector(Dimension);
        //base vector 1
        noalias(BaseVector1) = ZeroVector(Dimension);
        //base vector 2
        noalias(BaseVector2) = ZeroVector(Dimension);
        //base vector 3 normalized
        noalias(NormalVector) = ZeroVector(Dimension);
        //not-normalized base vector 3
        noalias(NormalVectorTilde) = ZeroVector(Dimension);
        //differential area
        DifferentialArea = 1.0;
    }

    Shell6pElement::ConstitutiveVariables::ConstitutiveVariables(std::size_t StrainSize)
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

        const std::size_t r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (m_A_ab_covariant_vector.size() != r_number_of_integration_points)
            m_A_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_B_ab_covariant_vector.size() != r_number_of_integration_points)
            m_B_ab_covariant_vector.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_T_vector.size() != r_number_of_integration_points)
            m_T_vector.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables(
            GetGeometry().WorkingSpaceDimension());

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            CalculateKinematics(
                point_number,
                kinematic_variables);

            m_A_ab_covariant_vector[point_number] = kinematic_variables.MetricCovariant;
            m_B_ab_covariant_vector[point_number] = kinematic_variables.CurvatureCovariant;

            m_dA_vector[point_number] = kinematic_variables.DifferentialArea;

            CalculateTransformationFromLocalToGlobalCartesian(kinematic_variables, m_T_vector[point_number]);
        }

        InitializeMaterial();
        mZeta = 0.0;

        KRATOS_CATCH("")
    }

    void Shell6pElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

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

    void Shell6pElement::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const std::size_t mat_size = GetGeometry().size() * 6;

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
        const std::size_t mat_size = GetGeometry().size() * 6;
        
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
        const std::size_t number_of_nodes = GetGeometry().size();
        const std::size_t mat_size = number_of_nodes * 6;

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

        const auto& r_geometry = GetGeometry();

        // definition of problem size
        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 6;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        // Shape function values for all integration points
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            double integration_weight = r_integration_points[point_number].Weight();

            double thickness = this->GetProperties().GetValue(THICKNESS);
            double density = this->GetProperties().GetValue(DENSITY);
            double mass = thickness * density * m_dA_vector[point_number] * integration_weight;

            if (rMassMatrix.size1() != mat_size)
                rMassMatrix.resize(mat_size, mat_size, false);

            rMassMatrix = ZeroMatrix(mat_size, mat_size);

            for (unsigned int r = 0; r<number_of_nodes; r++)
            {
                for (unsigned int s = 0; s<number_of_nodes; s++)
                {
                    rMassMatrix(6 * s, 6 * r) = r_N(point_number, s)*r_N(point_number, r) * mass;
                    rMassMatrix(6 * s + 1, 6 * r + 1) = rMassMatrix(6 * s, 6 * r);
                    rMassMatrix(6 * s + 2, 6 * r + 2) = rMassMatrix(6 * s, 6 * r);
                    rMassMatrix(6 * s + 3, 6 * r + 3) = rMassMatrix(6 * s, 6 * r) * thickness * thickness / 12.0;
                    rMassMatrix(6 * s + 4, 6 * r + 4) = rMassMatrix(6 * s, 6 * r) * thickness * thickness / 12.0;
                    rMassMatrix(6 * s + 5, 6 * r + 5) = rMassMatrix(6 * s, 6 * r) * thickness * thickness / 12.0;
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
            const std::size_t number_of_nodes = GetGeometry().size();
            const std::size_t mat_size = number_of_nodes * 6;
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
        // Retrieve geometry and integration points
        const GeometryType& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        // Definition of problem size
        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 6;

        // Provide a default empty implementation: resize and set zeros
        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        double thickness = this->GetProperties().GetValue(THICKNESS); 

        // Initialize components for non-linear analysis
        Vector current_displacement = ZeroVector(6*number_of_nodes);
        GetValuesVector(current_displacement,0);
        std::vector<array_1d<double, 6>> strain_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
     
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
            constitutive_variables.ConstitutiveMatrix = prod(m_T_vector[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(m_T_vector[point_number]))));
              
            // 3. Loop for zeta
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // Retrieve zeta for the current Gauss point in thickness direction
                mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

                // Compute the Jacobian matrix at the current Gauss point in thickness direction
                Matrix jacobian = ZeroMatrix(3,3);
                Matrix jacobian_inv = ZeroMatrix(3,3);
                double jacobian_det = 0.0;

                for (int i = 0; i < 3; i++) {
                    jacobian(0, i) = kinematic_variables.BaseVector1[i] + (thickness/2) * mZeta * normal_vector_derivatives(0, i); 
                    jacobian(1, i) = kinematic_variables.BaseVector2[i] + (thickness/2) * mZeta * normal_vector_derivatives(1, i) ; 
                    jacobian(2, i) = kinematic_variables.NormalVector[i] * (thickness/2) ;
                }

                MathUtils<double>::InvertMatrix(jacobian, jacobian_inv, jacobian_det);

                // Compute the b-operator at the current Gauss point in thickness direction
                Matrix b_operator = ZeroMatrix(6, mat_size);
                CalculateBOperator(point_number, b_operator, mZeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);

                strain_cau_cart[Gauss_index] = prod(b_operator, current_displacement);
                stress_cau_cart[Gauss_index] = prod(constitutive_variables.ConstitutiveMatrix,strain_cau_cart[Gauss_index]);
            }
        }
 
        // Cauchy stress at midspan
        array_1d<double, 6> stress_cau_cart_mid;
        stress_cau_cart_mid = (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness-1] + stress_cau_cart[0]) / 2.0;

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
        {
            if (rVariable == CAUCHY_STRESS_TOP_XX) {
                rOutput[point_number] = stress_cau_cart_mid[0] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][0]
                    - stress_cau_cart_mid[0]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_YY) {
                rOutput[point_number] = stress_cau_cart_mid[1] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][1]
                    - stress_cau_cart_mid[1]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_TOP_XY) {
                rOutput[point_number] = stress_cau_cart_mid[3] + (stress_cau_cart[mGaussIntegrationThickness.num_GP_thickness - 1][3]
                    - stress_cau_cart_mid[3]) / mGaussIntegrationThickness.integration_weight_thickness(mGaussIntegrationThickness.num_GP_thickness - 1);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XX) {
                rOutput[point_number] = stress_cau_cart_mid[0] + (stress_cau_cart[0][0] - stress_cau_cart_mid[0]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_YY) {
                rOutput[point_number] = stress_cau_cart_mid[1] + (stress_cau_cart[0][1] - stress_cau_cart_mid[1]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else if (rVariable == CAUCHY_STRESS_BOTTOM_XY) {
                rOutput[point_number] = stress_cau_cart_mid[3] + (stress_cau_cart[0][3] - stress_cau_cart_mid[3]) /
                    mGaussIntegrationThickness.integration_weight_thickness(0);
            }
            else {
                KRATOS_WATCH("No results for desired variable available in Calculate of Shell6pElement.")
            }
        }   
    }

    void Shell6pElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 >>& rVariable,
        std::vector<array_1d<double, 3 >>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }
    }

    void Shell6pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const std::size_t number_of_control_points = GetGeometry().size();

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

        const std::size_t number_of_control_points = GetGeometry().size();

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
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;

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

        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 6;
        double thickness = this->GetProperties().GetValue(THICKNESS); 

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
            constitutive_variables.ConstitutiveMatrix = prod(m_T_vector[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(m_T_vector[point_number]))));
                
            //// Loop for zeta (thickness integration points)
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // 3.1. Retrieve zeta for the current Gauss point in thickness direction
                mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

                // 3.2. Compute the Jacobian matrix at the current Gauss point in thickness direction
                Matrix jacobian = ZeroMatrix(3,3);
                Matrix jacobian_inv = ZeroMatrix(3,3);
                double jacobian_det = 0.0;

                for (int i = 0; i < 3; i++) {
                    jacobian(0, i) = kinematic_variables.BaseVector1[i] + (thickness/2) * mZeta * normal_vector_derivatives(0, i); 
                    jacobian(1, i) = kinematic_variables.BaseVector2[i] + (thickness/2) * mZeta * normal_vector_derivatives(1, i) ; 
                    jacobian(2, i) = kinematic_variables.NormalVector[i] * (thickness/2) ;
                }
                MathUtils<double>::InvertMatrix(jacobian, jacobian_inv, jacobian_det);

                // 3.3. Material stiffness part
                Matrix b_operator = ZeroMatrix(6, mat_size);
                CalculateBOperator(point_number, b_operator, mZeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);

                // 3.4. Drilling stiffness part
                Matrix b_drilling = ZeroMatrix(1, mat_size);
                CalculateBDrilling(point_number, b_drilling, jacobian_inv, kinematic_variables);

                // 3.5. Geometric stiffness part
                std::vector<array_1d<double, 6>> strain_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
                std::vector<array_1d<double, 6>> stress_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
                Vector current_displacement = ZeroVector(6*number_of_nodes);
                GetValuesVector(current_displacement,0);

                strain_cau_cart[Gauss_index] = prod(b_operator, current_displacement);
                stress_cau_cart[Gauss_index] = prod(constitutive_variables.ConstitutiveMatrix,strain_cau_cart[Gauss_index]);
                Matrix stress_matrix = ZeroMatrix(9,9);
                CalculateStressMatrix(stress_cau_cart[Gauss_index],stress_matrix);

                Matrix b_geometric = ZeroMatrix(9, mat_size);
                CalculateBGeometric(point_number, b_geometric, mZeta, jacobian_inv, normal_vector_derivatives, kinematic_variables);
                
                // 3.6. Compute the integration weight 
                double integration_weight = r_integration_points[point_number].Weight() * jacobian_det
                                          * mGaussIntegrationThickness.integration_weight_thickness(Gauss_index); 

                // 3.7. Assembly
                if (CalculateStiffnessMatrixFlag == true)
                {
                    // Add the linear stiffness matrix contribution to the element stiffness matrix
                    noalias(rLeftHandSideMatrix) += integration_weight *prod(trans(b_operator), Matrix(prod(constitutive_variables.ConstitutiveMatrix, b_operator)));

                    // Add the drilling stiffness matrix contribution to the element stiffness matrix
                    double E = this->GetProperties().GetValue(YOUNG_MODULUS);
                    double drilling_factor = 0.05; 
                    noalias(rLeftHandSideMatrix) += drilling_factor * E  * integration_weight * prod(trans(b_drilling), Matrix((b_drilling))); 
                    
                    // Add the geometric stiffness matrix contribution to the element stiffness matrix
                    noalias(rLeftHandSideMatrix) += integration_weight *prod(trans(b_geometric), Matrix(prod(stress_matrix, b_geometric)));
                }
                if (CalculateResidualVectorFlag == true)
                {
                    noalias(rRightHandSideVector) -= integration_weight * prod(trans(b_operator), stress_cau_cart[Gauss_index]);
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

        //not-normalized base vector 3
        MathUtils<double>::CrossProduct(rKinematicVariables.NormalVectorTilde, rKinematicVariables.BaseVector1, rKinematicVariables.BaseVector2);

        //differential area DifferentialArea
        rKinematicVariables.DifferentialArea = norm_2(rKinematicVariables.NormalVectorTilde);

        //base vector 3 normalized
        noalias(rKinematicVariables.NormalVector) = rKinematicVariables.NormalVectorTilde / rKinematicVariables.DifferentialArea;

        //GetCovariantMetric
        rKinematicVariables.MetricCovariant[0] = pow(rKinematicVariables.BaseVector1[0], 2) + pow(rKinematicVariables.BaseVector1[1], 2) + pow(rKinematicVariables.BaseVector1[2], 2);
        rKinematicVariables.MetricCovariant[1] = pow(rKinematicVariables.BaseVector2[0], 2) + pow(rKinematicVariables.BaseVector2[1], 2) + pow(rKinematicVariables.BaseVector2[2], 2);
        rKinematicVariables.MetricCovariant[2] = rKinematicVariables.BaseVector1[0] * rKinematicVariables.BaseVector2[0] + rKinematicVariables.BaseVector1[1] * rKinematicVariables.BaseVector2[1] + rKinematicVariables.BaseVector1[2] * rKinematicVariables.BaseVector2[2];
    }

    // Computes transformation matrix from local to global cartesian coordinates
    void Shell6pElement::CalculateTransformationFromLocalToGlobalCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTransformationMatrix
    ) const
    {
        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.BaseVector1);
        array_1d<double, 3> local_base_1 = rKinematicVariables.BaseVector1 / l_a1;
        array_1d<double, 3> local_base_3 =  rKinematicVariables.NormalVector;
        array_1d<double, 3> local_base_2;
        MathUtils<double>::CrossProduct(local_base_2, local_base_3, local_base_1);
        
        //Transformation matrix 
        if (rTransformationMatrix.size1() != 6 && rTransformationMatrix.size2() != 6)                                                                 
            rTransformationMatrix.resize(6, 6);
        noalias(rTransformationMatrix) = ZeroMatrix(6, 6);

        for (std::size_t i = 0; i < 3; ++i)
        {
            std::size_t j = (i + 1) % 3;

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

        array_1d<double, 6> strain_vector(6, 0.0);
        //To do: compute strain vector
        noalias(rThisConstitutiveVariables.StrainVector) = prod(m_T_vector[IntegrationPointIndex], strain_vector);

        // Constitive Matrices D
        rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //this is an ouput parameter

        const double nu = this->GetProperties()[POISSON_RATIO];
        const double Emodul = this->GetProperties()[YOUNG_MODULUS];
        double lambda = Emodul / (1.0 - nu * nu);
        double Gmodul = Emodul / (2.0 * (1.0 + nu));
        
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 0) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(0, 1) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 0) = lambda * nu;
        rThisConstitutiveVariables.ConstitutiveMatrix(1, 1) = lambda;
        rThisConstitutiveVariables.ConstitutiveMatrix(3, 3) = lambda * (1 - nu) / 2;
        rThisConstitutiveVariables.ConstitutiveMatrix(4, 4) = Gmodul * 5.0 / 6.0;
        rThisConstitutiveVariables.ConstitutiveMatrix(5, 5) = Gmodul * 5.0 / 6.0;

        //Local Cartesian Stresses
        noalias(rThisConstitutiveVariables.StressVector) = prod(
            trans(rThisConstitutiveVariables.ConstitutiveMatrix), rThisConstitutiveVariables.StrainVector);
    }

    void Shell6pElement::CalculateNormalVectorDerivatives(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        Matrix& rNormalVectorDerivatives) const
    {
        // Get the shape function of second derivatives
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        // Get the area of the element
        double inv_differential_area = 1 / rKinematicVariables.DifferentialArea;
        double inv_differential_area_cube = 1 / std::pow(rKinematicVariables.DifferentialArea, 3);

        // Compute base vector derivatives
        array_1d<double, 3> base_vector1_derivative_11;
        array_1d<double, 3> base_vector1_derivative_12; 
        array_1d<double, 3> base_vector2_derivative_22;

        for (std::size_t i = 0; i < GetGeometry().size(); ++i)
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
        Matrix& rB,
        double zeta,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        double thickness = this->GetProperties().GetValue(THICKNESS);
        
        // Membrane part
        Matrix DN_De_Jn = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De(i, 0) = r_DN_De(i, 0);  // Copy first column
            new_DN_De(i, 1) = r_DN_De(i, 1);  // Copy second column
            new_DN_De(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn = trans(prod(J_inv, trans(new_DN_De)));

        // Bending part
        Matrix DN_De_Jn_bending = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De_bending = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De_bending(i, 0) = r_DN_De(i, 0) * (thickness/2) * zeta;  // Copy first column
            new_DN_De_bending(i, 1) = r_DN_De(i, 1) * (thickness/2) * zeta;  // Copy second column
            new_DN_De_bending(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn_bending = trans(prod(J_inv, trans(new_DN_De_bending)));
                                                
        // y hat 
        const double y1= rActualKinematic.NormalVector[0];
        const double y2= rActualKinematic.NormalVector[1];
        const double y3= rActualKinematic.NormalVector[2];

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix Dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        //compute da3
        array_1d<double, 3> da1_d1;
        array_1d<double, 3> da1_d2; //da1_d2 = da2_d1
        array_1d<double, 3> da2_d2;

        array_1d<double, 3> da3_tilde_d1;
        array_1d<double, 3> da3_tilde_d1_1;
        array_1d<double, 3> da3_tilde_d1_2;
        array_1d<double, 3> da3_tilde_d2;
        array_1d<double, 3> da3_tilde_d2_1;
        array_1d<double, 3> da3_tilde_d2_2;

        if (rB.size1() != 6 || rB.size2() != mat_size)
            rB.resize(6, mat_size);
        noalias(rB) = ZeroMatrix(6, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            // zeta  Derivatives 
            const double dzetadx= J_inv(0,2);
            const double dzetady= J_inv(1,2);
            const double dzetadz= J_inv(2,2); 
          
            Dn = prod(J_inv, dn);

            // y hat Derivative w.r.t x
            const double dy1x= Dn(0,0);
            const double dy2x= Dn(0,1);
            const double dy3x= Dn(0,2);

            // y hat Derivative w.r.t y
            const double dy1y= Dn(1,0);
            const double dy2y= Dn(1,1);
            const double dy3y= Dn(1,2);

            // y hat Derivative w.r.t z
            const double dy1z= Dn(2,0);
            const double dy2z= Dn(2,1);
            const double dy3z= Dn(2,2);

            // Compute rB
            IndexType index = i * 6;
            
            rB(0, index)     = DN_De_Jn(i, 0);
            rB(0, index + 1) = 0;
            rB(0, index + 2) = 0;

            rB(0, index + 3) = 0;
            rB(0, index + 4) = (DN_De_Jn_bending(i, 0) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3x + y3 * dzetadx));
            rB(0, index + 5) = - ((DN_De_Jn_bending(i, 0) * y2) + (r_N(i) * (thickness/2) * (zeta * dy2x + dzetadx * y2)));

            rB(1, index)     = 0;
            rB(1, index + 1) = DN_De_Jn(i, 1);
            rB(1, index + 2) = 0;

            rB(1, index + 3) = - ((DN_De_Jn_bending(i, 1) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3y + dzetady * y3))); 
            rB(1, index + 4) = 0;
            rB(1, index + 5) = (DN_De_Jn_bending(i, 1) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1y + dzetady * y1)); 

            rB(2, index)     = 0;
            rB(2, index + 1) = 0;
            rB(2, index + 2) = DN_De_Jn(i, 2);

            rB(2, index + 3) = (DN_De_Jn_bending(i, 2)   * y2) + (r_N(i)  * (thickness/2) * (zeta * dy2z + dzetadz * y2));  
            rB(2, index + 4) = - ((DN_De_Jn_bending(i, 2) * y1)  + (r_N(i) * (thickness/2) * (zeta  * dy1z + dzetadz *y1))); 
            rB(2, index + 5) = 0;
            
            rB(3, index)     = DN_De_Jn(i, 1);                    
            rB(3, index + 1) = DN_De_Jn(i, 0);  
            rB(3, index + 2) = 0;

            rB(3, index + 3) = - ((DN_De_Jn_bending(i, 0) * y3) +(r_N(i) * (thickness/2) * (zeta * dy3x + dzetadx * y3)));
            rB(3, index + 4) = (DN_De_Jn_bending(i, 1) * y3) +(r_N(i) * (thickness/2) * (zeta * dy3y + dzetady * y3));
            rB(3, index + 5) = ((DN_De_Jn_bending(i, 0) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1x + dzetadx * y1))) - ((DN_De_Jn_bending(i, 1) * y2)+ (r_N(i) * (thickness/2) *  (zeta * dy2y + dzetady * y2))); 

            rB(4, index)     = 0;
            rB(4, index + 1) = DN_De_Jn(i, 2);
            rB(4, index + 2) = DN_De_Jn(i, 1);

            rB(4, index + 3) = ((DN_De_Jn_bending(i, 1) * y2) + (r_N(i) * (thickness/2) * (zeta * dy2y + dzetady * y2)))  - ((DN_De_Jn_bending(i, 2) * y3)+ (r_N (i)  * (thickness/2) * (zeta *dy3z + dzetadz * y3))); 
            rB(4, index + 4) = - ((DN_De_Jn_bending(i, 1) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1y + dzetady * y1))); 
            rB(4, index + 5) = (DN_De_Jn_bending(i, 2) *  y1 ) + (r_N (i) * (thickness/2) * (zeta * dy1z + dzetadz * y1)); 

            rB(5, index)   = DN_De_Jn (i,2);                    
            rB(5, index + 1) = 0;  
            rB(5, index + 2) = DN_De_Jn(i, 0);

            rB(5, index + 3) = ((DN_De_Jn_bending(i, 0)  * y2) + (r_N(i) * (thickness/2) * ( zeta * dy2x +  dzetadx * y2)) ); 
            rB(5, index + 4) = ((DN_De_Jn_bending (i,2)  * y3) + (r_N (i)   * (thickness/2) * ( zeta * dy3z + dzetadz * y3 ))) - ((DN_De_Jn_bending(i, 0) * y1) + (r_N(i) * (thickness/2) * (zeta  * dy1x + dzetadx * y1) )); 
            rB(5, index + 5) = - ((DN_De_Jn_bending (i,2)  * y2 )+(r_N (i) * (thickness/2) * (zeta * dy2z + dzetadz * y2)));  
        }
    }

    void Shell6pElement::CalculateBDrilling(                                                                                         
        const IndexType IntegrationPointIndex,
        Matrix& rBd,
        Matrix& J_inv,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;

        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

        // Membrane part
        Matrix DN_De_Jn = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De(i, 0) = r_DN_De(i, 0);  // Copy first column
            new_DN_De(i, 1) = r_DN_De(i, 1);  // Copy second column
            new_DN_De(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn = trans(prod(J_inv, trans(new_DN_De)));

        if (rBd.size1() != 1|| rBd.size2() != mat_size)                                 
            rBd.resize(1, mat_size);                                                     
        noalias(rBd) = ZeroMatrix(1, mat_size);                                          

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            IndexType index = i * 6;

            rBd(0, index)     = -0.5 * DN_De_Jn(i, 1);
            rBd(0, index + 1) = 0.5 * DN_De_Jn(i, 0);
            rBd(0, index + 2) = 0;
            rBd(0, index + 3) = 0;
            rBd(0, index + 4) = 0;
            rBd(0, index + 5) = - r_N(i);
        }
    }

    void Shell6pElement::CalculateBGeometric(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        double thickness = this->GetProperties().GetValue(THICKNESS); 
        
        // Membrane part
        Matrix DN_De_Jn = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De(i, 0) = r_DN_De(i, 0);  // Copy first column
            new_DN_De(i, 1) = r_DN_De(i, 1);  // Copy second column
            new_DN_De(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn = trans(prod(J_inv, trans(new_DN_De)));

        // Bending part
        Matrix DN_De_Jn_bending = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De_bending = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De_bending(i, 0) = r_DN_De(i, 0) * (thickness/2) * zeta;  // Copy first column
            new_DN_De_bending(i, 1) = r_DN_De(i, 1) * (thickness/2) * zeta;  // Copy second column
            new_DN_De_bending(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn_bending = trans(prod(J_inv, trans(new_DN_De_bending)));
                                                
        // y hat 
        const double y1= rActualKinematic.NormalVector[0];
        const double y2= rActualKinematic.NormalVector[1];
        const double y3= rActualKinematic.NormalVector[2];

        Matrix da3 = ZeroMatrix(3, 3);
        Matrix Dn = ZeroMatrix(3, 3);
        Matrix b = ZeroMatrix(3, mat_size);

        //compute da3
        array_1d<double, 3> da1_d1;
        array_1d<double, 3> da1_d2; //da1_d2 = da2_d1
        array_1d<double, 3> da2_d2;

        array_1d<double, 3> da3_tilde_d1;
        array_1d<double, 3> da3_tilde_d1_1;
        array_1d<double, 3> da3_tilde_d1_2;
        array_1d<double, 3> da3_tilde_d2;
        array_1d<double, 3> da3_tilde_d2_1;
        array_1d<double, 3> da3_tilde_d2_2;

        if (rB.size1() != 9 || rB.size2() != mat_size)
            rB.resize(9, mat_size);
        noalias(rB) = ZeroMatrix(9, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            // zeta  Derivatives 
            const double dzetadx= J_inv(0,2);
            const double dzetady= J_inv(1,2);
            const double dzetadz= J_inv(2,2); 
          
            Dn = prod(J_inv, dn);

            // y hat Derivative w.r.t x
            const double dy1x= Dn(0,0);
            const double dy2x= Dn(0,1);
            const double dy3x= Dn(0,2);

            // y hat Derivative w.r.t y
            const double dy1y= Dn(1,0);
            const double dy2y= Dn(1,1);
            const double dy3y= Dn(1,2);

            // y hat Derivative w.r.t z
            const double dy1z= Dn(2,0);
            const double dy2z= Dn(2,1);
            const double dy3z= Dn(2,2);

            // Compute rB
            IndexType index = i * 6;
            
            rB(0, index)     = DN_De_Jn(i, 0);
            rB(0, index + 1) = 0;
            rB(0, index + 2) = 0;
            rB(0, index + 3) = 0;
            rB(0, index + 4) = (DN_De_Jn_bending(i, 0) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3x + y3 * dzetadx));
            rB(0, index + 5) = - ((DN_De_Jn_bending(i, 0) * y2) + (r_N(i) * (thickness/2) * (zeta * dy2x + dzetadx * y2)));

            rB(1, index)     = DN_De_Jn(i, 1);
            rB(1, index + 1) = 0;
            rB(1, index + 2) = 0;
            rB(1, index + 3) = 0;
            rB(1, index + 4) = (DN_De_Jn_bending(i, 1) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3y + y3 * dzetady));
            rB(1, index + 5) = - ((DN_De_Jn_bending(i, 1) * y2) + (r_N(i) * (thickness/2) * (zeta * dy2y + dzetady * y2)));

            rB(2, index)     = DN_De_Jn(i, 2);
            rB(2, index + 1) = 0;
            rB(2, index + 2) = 0;
            rB(2, index + 3) = 0;
            rB(2, index + 4) = (DN_De_Jn_bending(i, 2) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3z + y3 * dzetadz));
            rB(2, index + 5) = - ((DN_De_Jn_bending(i, 2) * y2) + (r_N(i) * (thickness/2) * (zeta * dy2z + dzetadz * y2)));

            rB(3, index)     = 0;
            rB(3, index + 1) = DN_De_Jn(i, 0);
            rB(3, index + 2) = 0;
            rB(3, index + 3) = - ((DN_De_Jn_bending(i, 0) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3x + dzetadx * y3))); 
            rB(3, index + 4) = 0;
            rB(3, index + 5) = (DN_De_Jn_bending(i, 0) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1x + dzetadx * y1)); 

            rB(4, index)     = 0;
            rB(4, index + 1) = DN_De_Jn(i, 1);
            rB(4, index + 2) = 0;
            rB(4, index + 3) = - ((DN_De_Jn_bending(i, 1) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3y + dzetady * y3))); 
            rB(4, index + 4) = 0;
            rB(4, index + 5) = (DN_De_Jn_bending(i, 1) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1y + dzetady * y1)); 

            rB(5, index)     = 0;
            rB(5, index + 1) = DN_De_Jn(i, 2);
            rB(5, index + 2) = 0;
            rB(5, index + 3) = - ((DN_De_Jn_bending(i, 2) * y3) + (r_N(i) * (thickness/2) * (zeta * dy3z + dzetadz * y3))); 
            rB(5, index + 4) = 0;
            rB(5, index + 5) = (DN_De_Jn_bending(i, 2) * y1) + (r_N(i) * (thickness/2) * (zeta * dy1z + dzetadz * y1)); 

            rB(6, index)     = 0;
            rB(6, index + 1) = 0;
            rB(6, index + 2) = DN_De_Jn(i, 0);
            rB(6, index + 3) = (DN_De_Jn_bending(i, 0)   * y2) + (r_N(i)  * (thickness/2) * (zeta * dy2x + dzetadx * y2));  
            rB(6, index + 4) = - ((DN_De_Jn_bending(i, 0) * y1)  + (r_N(i) * (thickness/2) * (zeta  * dy1x + dzetadx *y1))); 
            rB(6, index + 5) = 0;

            rB(7, index)     = 0;
            rB(7, index + 1) = 0;
            rB(7, index + 2) = DN_De_Jn(i, 1);
            rB(7, index + 3) = (DN_De_Jn_bending(i, 1)   * y2) + (r_N(i)  * (thickness/2) * (zeta * dy2y + dzetady * y2));  
            rB(7, index + 4) = - ((DN_De_Jn_bending(i, 1) * y1)  + (r_N(i) * (thickness/2) * (zeta  * dy1y + dzetady *y1))); 
            rB(7, index + 5) = 0;

            rB(8, index)     = 0;
            rB(8, index + 1) = 0;
            rB(8, index + 2) = DN_De_Jn(i, 2);
            rB(8, index + 3) = (DN_De_Jn_bending(i, 2)   * y2) + (r_N(i)  * (thickness/2) * (zeta * dy2z + dzetadz * y2));  
            rB(8, index + 4) = - ((DN_De_Jn_bending(i, 2) * y1)  + (r_N(i) * (thickness/2) * (zeta  * dy1z + dzetadz *y1))); 
            rB(8, index + 5) = 0;
        }
    }

    void Shell6pElement::CalculateStressMatrix(
        array_1d<double, 6> stress_vector,
        Matrix& stress_matrix
    ) const
    {
        Matrix stress_mat = ZeroMatrix(3,3);
        
        stress_mat(0,0) = stress_vector(0);
        stress_mat(0,1) = stress_vector(3);
        stress_mat(0,2) = stress_vector(4);
        
        stress_mat(1,0) = stress_mat(0,1);
        stress_mat(1,1) = stress_vector(1);
        stress_mat(1,2) = stress_vector(5);
        
        stress_mat(2,0) = stress_mat(0,2);
        stress_mat(2,1) = stress_mat(1,2);
        stress_mat(2,2) = stress_vector(2);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                stress_matrix(i,j) = stress_mat(i,j);
                stress_matrix(i + 3,j + 3) = stress_mat(i,j);
                stress_matrix(i + 6,j + 6) = stress_mat(i,j);
            }
        }
    }
    ///@}

} // Namespace Kratos


