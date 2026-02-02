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
        noalias(a_ab_covariant) = ZeroVector(Dimension);
        noalias(b_ab_covariant) = ZeroVector(Dimension);
        //base vector 1
        noalias(a1) = ZeroVector(Dimension);
        //base vector 2
        noalias(a2) = ZeroVector(Dimension);
        //base vector 3 normalized
        noalias(a3) = ZeroVector(Dimension);
        //not-normalized base vector 3
        noalias(a3_tilde) = ZeroVector(Dimension);
        //differential area
        dA = 1.0;
    }

    Shell6pElement::ConstitutiveVariables::ConstitutiveVariables(std::size_t StrainSize)
    {
        StrainVector       = ZeroVector(StrainSize);
        StressVector       = ZeroVector(StrainSize);
        ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
    }

    Shell6pElement::SecondVariations::SecondVariations(const int& mat_size)
    {
        B11 = ZeroMatrix(mat_size, mat_size);
        B22 = ZeroMatrix(mat_size, mat_size);
        B12 = ZeroMatrix(mat_size, mat_size);
    }


    ///@name Initialize Functions
    ///@{

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

            m_A_ab_covariant_vector[point_number] = kinematic_variables.a_ab_covariant;
            m_B_ab_covariant_vector[point_number] = kinematic_variables.b_ab_covariant;

            m_dA_vector[point_number] = kinematic_variables.dA;

            CalculateTransformation(kinematic_variables, m_T_vector[point_number]);
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

    void Shell6pElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
        }
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
        const GeometryType& r_geometry = GetGeometry();

        // definition of problem size
        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 6;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        Vector current_displacement = ZeroVector(6*number_of_nodes);
        GetValuesVector(current_displacement,0);

        // Initialize strain and stress
        std::vector<array_1d<double, 6>> strain_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
        std::vector<array_1d<double, 6>> stress_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
     
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(point_number, kinematic_variables);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables(6);

            CalculateConstitutiveVariables(point_number, kinematic_variables, constitutive_variables,
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

            constitutive_variables.ConstitutiveMatrix = prod(m_T_vector[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(m_T_vector[point_number]))));
              

            // Loop for zeta
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {
                // Initialization
                Matrix B = ZeroMatrix(6, mat_size);
                Matrix dn = ZeroMatrix(3, 3);
                Matrix DN_De_Jn = ZeroMatrix(number_of_nodes,3);
                Matrix J_inv = ZeroMatrix(3, 3);
                double area = 0.0;

                CalculateJn(point_number, kinematic_variables, mZeta, DN_De_Jn, J_inv, dn, area);

                CalculateB(point_number, B, mZeta, DN_De_Jn, J_inv, dn, kinematic_variables);

                strain_cau_cart[Gauss_index] = prod(B, current_displacement);
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


    ///@}
    ///@name Assembly
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

        const GeometryType& r_geometry = GetGeometry();

        // definition of problem size
        const std::size_t number_of_nodes = r_geometry.size();
        const std::size_t mat_size = number_of_nodes * 6;

        const auto& r_integration_points = r_geometry.IntegrationPoints();

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(point_number, kinematic_variables);
            
            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables(6);

            CalculateConstitutiveVariables(point_number, kinematic_variables, constitutive_variables,
                constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

            constitutive_variables.ConstitutiveMatrix = prod(m_T_vector[point_number], 
                Matrix(prod(constitutive_variables.ConstitutiveMatrix, trans(m_T_vector[point_number]))));
                
            // Loop for zeta
            for (IndexType Gauss_index = 0; Gauss_index < mGaussIntegrationThickness.num_GP_thickness; Gauss_index++)
            {            
                // Initialization
                Matrix B = ZeroMatrix(6, mat_size);
                Matrix B_Drill = ZeroMatrix(1, mat_size);
                Matrix B_Geometric = ZeroMatrix(9, mat_size);
                Matrix dn = ZeroMatrix(3, 3);
                Matrix stress_matrix = ZeroMatrix(9,9);
                Matrix DN_De_Jn = ZeroMatrix(number_of_nodes,3);
                Matrix J_inv = ZeroMatrix(3, 3);
                double area = 0.0;

                mZeta = mGaussIntegrationThickness.zeta(Gauss_index);

                CalculateJn(point_number, kinematic_variables, mZeta, DN_De_Jn, J_inv, dn, area);

                CalculateB(point_number, B, mZeta, DN_De_Jn, J_inv, dn, kinematic_variables);

                CalculateBDrill(point_number, B_Drill, DN_De_Jn, kinematic_variables);

                //Geometric stiffness part (TO DO)
                CalculateBGeometric(point_number, B_Geometric, mZeta, DN_De_Jn, J_inv, dn, kinematic_variables);

                // Initialize strain and stress
                std::vector<array_1d<double, 6>> strain_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
                std::vector<array_1d<double, 6>> stress_cau_cart(mGaussIntegrationThickness.num_GP_thickness);
                Vector current_displacement = ZeroVector(6*number_of_nodes);
                GetValuesVector(current_displacement,0);

                strain_cau_cart[Gauss_index] = prod(B, current_displacement);
                stress_cau_cart[Gauss_index] = prod(constitutive_variables.ConstitutiveMatrix,strain_cau_cart[Gauss_index]);

                CalculateStressMatrix(stress_cau_cart[Gauss_index],stress_matrix);
                
                double integration_weight =
                    r_integration_points[point_number].Weight()
                    * area; // * m_dA_vector[point_number]; 

                Matrix rKm = ZeroMatrix(mat_size, mat_size);
                Matrix rKd = ZeroMatrix(mat_size, mat_size);

                if (CalculateStiffnessMatrixFlag == true)
                {
                    CalculateAndAddKm(rKm, B, constitutive_variables.ConstitutiveMatrix);

                    CalculateAndAddKmBd(rKd, B_Drill);

                    CalculateAndAddK(rLeftHandSideMatrix, rKm, rKd, integration_weight,
                        mGaussIntegrationThickness.integration_weight_thickness(Gauss_index));

                    CalculateAndAddNonlinearKm(rLeftHandSideMatrix, B_Geometric,stress_matrix,
                       integration_weight, mGaussIntegrationThickness.integration_weight_thickness(Gauss_index));
                }
                
                if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                {
                    // operation performed: rRightHandSideVector -= Weight*IntForce
                    noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), stress_cau_cart[Gauss_index]);
                }

            } 
        }
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Implicit
    ///@{
    
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

    ///@}
    ///@name Kinematics
    ///@{

    void Shell6pElement::CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables
    ) const
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
        noalias(rKinematicVariables.a3) = rKinematicVariables.a3_tilde / rKinematicVariables.dA;

        //GetCovariantMetric
        rKinematicVariables.a_ab_covariant[0] = pow(rKinematicVariables.a1[0], 2) + pow(rKinematicVariables.a1[1], 2) + pow(rKinematicVariables.a1[2], 2);
        rKinematicVariables.a_ab_covariant[1] = pow(rKinematicVariables.a2[0], 2) + pow(rKinematicVariables.a2[1], 2) + pow(rKinematicVariables.a2[2], 2);
        rKinematicVariables.a_ab_covariant[2] = rKinematicVariables.a1[0] * rKinematicVariables.a2[0] + rKinematicVariables.a1[1] * rKinematicVariables.a2[1] + rKinematicVariables.a1[2] * rKinematicVariables.a2[2];

        // Matrix H = ZeroMatrix(3, 3);
        // CalculateHessian(H, GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod()));

        // rKinematicVariables.b_ab_covariant[0] = H(0, 0) * rKinematicVariables.a3[0] + H(1, 0) * rKinematicVariables.a3[1] + H(2, 0) * rKinematicVariables.a3[2];
        // rKinematicVariables.b_ab_covariant[1] = H(0, 1) * rKinematicVariables.a3[0] + H(1, 1) * rKinematicVariables.a3[1] + H(2, 1) * rKinematicVariables.a3[2];
        // rKinematicVariables.b_ab_covariant[2] = H(0, 2) * rKinematicVariables.a3[0] + H(1, 2) * rKinematicVariables.a3[1] + H(2, 2) * rKinematicVariables.a3[2];
    }

    // Computes the transformation matrix
    void Shell6pElement::CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT
    ) const
    {
        //Local cartesian coordinates
        double l_a1 = norm_2(rKinematicVariables.a1);
        array_1d<double, 3> e1 = rKinematicVariables.a1 / l_a1;
        array_1d<double, 3> e3 =  rKinematicVariables.a3;
        array_1d<double, 3> e2;
        MathUtils<double>::CrossProduct(e2, e3, e1);
        
        //Transformation matrix 
        if (rT.size1() != 6 && rT.size2() != 6)                                                                 
            rT.resize(6, 6);
        noalias(rT) = ZeroMatrix(6, 6);

        for (std::size_t i = 0; i < 3; ++i)
        {
            std::size_t j = (i + 1) % 3;

            rT(i, 0) = e1[i] * e1[i];
            rT(i, 1) = e2[i] * e2[i];
            rT(i, 2) = e3[i] * e3[i];
            rT(i, 3) = 2 * e1[i] * e2[i];
            rT(i, 4) = 2 * e2[i] * e3[i];
            rT(i, 5) = 2 * e1[i] * e3[i];

            rT(i + 3, 0) = e1[i] * e1[j];
            rT(i + 3, 1) = e2[i] * e2[j];
            rT(i + 3, 2) = e3[i] * e3[j];
            rT(i + 3, 3) = (e1[i] * e2[j]) + (e2[i] * e1[j]);
            rT(i + 3, 4) = (e2[i] * e3[j]) + (e3[i] * e2[j]);
            rT(i + 3, 5) = (e1[i] * e3[j]) + (e3[i] * e1[j]);
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


    void Shell6pElement::CalculateJn(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        double& area) const
    {

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());

        double thickness = this->GetProperties().GetValue(THICKNESS);  

        Matrix J;
        GetGeometry().Jacobian(J, IntegrationPointIndex);
        const std::size_t number_of_control_points = GetGeometry().size();

        rKinematicVariables.a1 = column(J, 0);
        rKinematicVariables.a2 = column(J, 1);
        MathUtils<double>::CrossProduct(rKinematicVariables.a3_tilde, rKinematicVariables.a1, rKinematicVariables.a2);
        rKinematicVariables.dA = norm_2(rKinematicVariables.a3_tilde);
        noalias(rKinematicVariables.a3) = rKinematicVariables.a3_tilde / rKinematicVariables.dA;


        Matrix da3 = ZeroMatrix(3, 3);
        double inv_dA = 1 / rKinematicVariables.dA;
        double inv_dA3 = 1 / std::pow(rKinematicVariables.dA, 3);

        //compute da3
        array_1d<double, 3> da1_d1;
        array_1d<double, 3> da1_d2; //da1_d2 = da2_d1
        array_1d<double, 3> da2_d2;

        for (std::size_t i=0;i<number_of_control_points;++i){
            da1_d1[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 0);
            da1_d1[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 0);
            da1_d1[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 0);

            da1_d2[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 1);
            da1_d2[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 1);
            da1_d2[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 1);

            da2_d2[0] += (GetGeometry().GetPoint( i ).X0()) * r_DDN_DDe(i, 2);
            da2_d2[1] += (GetGeometry().GetPoint( i ).Y0()) * r_DDN_DDe(i, 2);
            da2_d2[2] += (GetGeometry().GetPoint( i ).Z0()) * r_DDN_DDe(i, 2);
        }

        array_1d<double, 3> da3_tilde_d1;
        array_1d<double, 3> da3_tilde_d1_1;
        array_1d<double, 3> da3_tilde_d1_2;
        array_1d<double, 3> da3_tilde_d2;
        array_1d<double, 3> da3_tilde_d2_1;
        array_1d<double, 3> da3_tilde_d2_2;

        MathUtils<double>::CrossProduct(da3_tilde_d1_1, da1_d1, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(da3_tilde_d1_2, rKinematicVariables.a1, da1_d2);
        da3_tilde_d1 = da3_tilde_d1_1 + da3_tilde_d1_2;

        MathUtils<double>::CrossProduct(da3_tilde_d2_1, da1_d2, rKinematicVariables.a2);
        MathUtils<double>::CrossProduct(da3_tilde_d2_2, rKinematicVariables.a1, da2_d2);
        da3_tilde_d2 = da3_tilde_d2_1 + da3_tilde_d2_2;

        for (IndexType j = 0; j < 3; j++)
        {
            dn(0, j) = da3_tilde_d1[j] * inv_dA - rKinematicVariables.a3_tilde[j] * inner_prod(rKinematicVariables.a3_tilde, da3_tilde_d1) * inv_dA3;
            dn(1, j) = da3_tilde_d2[j] * inv_dA - rKinematicVariables.a3_tilde[j] * inner_prod(rKinematicVariables.a3_tilde, da3_tilde_d2) * inv_dA3;
            dn(2, j) = 0.0;
        }

        Matrix Jn = ZeroMatrix(3,3);
        for (int i = 0; i < 3; i++) {
            Jn(0, i) = rKinematicVariables.a1[i] + (thickness/2) * zeta * dn(0, i); 
            Jn(1, i) = rKinematicVariables.a2[i] + (thickness/2) * zeta * dn(1, i) ; 
            Jn(2, i) = rKinematicVariables.a3[i] * (thickness/2) ;
        }

        // Jn(0,0) = 4.0; //Jn must be computed in the reference configuration! Warning: This is just a hard coded test

        Matrix Jn_inv = ZeroMatrix(3,3);
        double det_Jn = 0.0;
        MathUtils<double>::InvertMatrix(Jn, Jn_inv, det_Jn);

        Matrix new_DN_De = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De(i, 0) = r_DN_De(i, 0);  // Copy first column
            new_DN_De(i, 1) = r_DN_De(i, 1);  // Copy second column
            new_DN_De(i, 2) = 0.0;            // Set third column to zero
        }

        DN_De_Jn = trans(prod(Jn_inv, trans(new_DN_De)));  // DN_De_Jn = Jn_inv * r_DN_De that has 3 columns
        J_inv = Jn_inv;
        area = det_Jn;
    }

    void Shell6pElement::CalculateB(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());
        double thickness = this->GetProperties().GetValue(THICKNESS);   

        //Bending part
        Matrix DN_De_Jn_bending = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De_bending = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De_bending(i, 0) = r_DN_De(i, 0) * (thickness/2) * zeta;  // Copy first column
            new_DN_De_bending(i, 1) = r_DN_De(i, 1) * (thickness/2) * zeta;  // Copy second column
            new_DN_De_bending(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn_bending = trans(prod(J_inv, trans(new_DN_De_bending)));
                                                
        // y hat 
        const double y1= rActualKinematic.a3[0];
        const double y2= rActualKinematic.a3[1];
        const double y3= rActualKinematic.a3[2];

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

        double inv_dA = 1 / rActualKinematic.dA;
        double inv_dA3 = 1 / std::pow(rActualKinematic.dA, 3);

        if (rB.size1() != 6 || rB.size2() != mat_size)
            rB.resize(6, mat_size);
        noalias(rB) = ZeroMatrix(6, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
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

    void Shell6pElement::CalculateBDrill(                                                                                         
        const IndexType IntegrationPointIndex,
        Matrix& rBd,
        Matrix& DN_De_Jn,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;

        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        double thickness = this->GetProperties().GetValue(THICKNESS);
        const auto& r_N = GetGeometry().ShapeFunctionsValues();

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
            rBd(0, index + 5) = - r_N (i);
        }
    }

    void Shell6pElement::CalculateBGeometric(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 6;
        const auto& r_N = GetGeometry().ShapeFunctionsValues();
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);
        const Matrix& r_DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, IntegrationPointIndex, GetGeometry().GetDefaultIntegrationMethod());
        double thickness = this->GetProperties().GetValue(THICKNESS);   

        //Bending part
        Matrix DN_De_Jn_bending = ZeroMatrix(number_of_control_points,3);
        Matrix new_DN_De_bending = ZeroMatrix(number_of_control_points,3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            new_DN_De_bending(i, 0) = r_DN_De(i, 0) * (thickness/2) * zeta;  // Copy first column
            new_DN_De_bending(i, 1) = r_DN_De(i, 1) * (thickness/2) * zeta;  // Copy second column
            new_DN_De_bending(i, 2) = 0.0;            // Set third column to zero
        }
        DN_De_Jn_bending = trans(prod(J_inv, trans(new_DN_De_bending)));
                                                
        // y hat 
        const double y1= rActualKinematic.a3[0];
        const double y2= rActualKinematic.a3[1];
        const double y3= rActualKinematic.a3[2];

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

        double inv_dA = 1 / rActualKinematic.dA;
        double inv_dA3 = 1 / std::pow(rActualKinematic.dA, 3);

        if (rB.size1() != 9 || rB.size2() != mat_size)
            rB.resize(9, mat_size);
        noalias(rB) = ZeroMatrix(9, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
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
    ///@name Stiffness matrix assembly
    ///@{

    inline void Shell6pElement::CalculateAndAddK(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rKm,
        const Matrix& rKd,                                                                                                               
        const double IntegrationWeight,
        const double IntegrationWeight_zeta 
    ) const
    {
        double thickness = this->GetProperties().GetValue(THICKNESS);
        noalias(rLeftHandSideMatrix) +=  IntegrationWeight * IntegrationWeight_zeta * (rKm + rKd ); 
        double weight = IntegrationWeight * IntegrationWeight_zeta * thickness/2;                               
    }

 
    inline void Shell6pElement::CalculateAndAddKm(
        MatrixType& rKm,
        const Matrix& rB,
        const Matrix& rD                                                                                                              
    ) const
    {  
        noalias(rKm) += prod(trans(rB), Matrix(prod(rD, rB)));                                              
    }

    
    inline void Shell6pElement::CalculateAndAddKmBd(                                                              
        MatrixType& rKd,
        const Matrix& rBd                                                                                                              
    ) const
    {
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);
        double thickness = this->GetProperties().GetValue(THICKNESS);

        noalias(rKd) += 0.05 * E  * prod(trans(rBd), Matrix((rBd)));                  
    }


    inline void Shell6pElement::CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight,
        const double IntegrationWeight_zeta) const
    {
        noalias(rLeftHandSideMatrix) +=  IntegrationWeight * IntegrationWeight_zeta * prod(trans(rB), Matrix(prod(rD, rB))); 
    }


    ///@}
    ///@name Stress recovery
    ///@{

    ///@}
    ///@name Dynamic Functions
    ///@{

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

    // TO DO
    void Shell6pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 3;

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

    // TO DO
    void Shell6pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const std::size_t number_of_control_points = GetGeometry().size();
        const std::size_t mat_size = number_of_control_points * 3;

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

    ///@}
    ///@name Check
    ///@{

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

    // void Shell6pElement::CalculateHessian(
    //     Matrix& Hessian,
    //     const Matrix& rDDN_DDe) const
    // {
    //     const std::size_t number_of_points = GetGeometry().size();
    //     const std::size_t working_space_dimension = 3;
    //     Hessian.resize(working_space_dimension, working_space_dimension);
    //     Hessian = ZeroMatrix(working_space_dimension, working_space_dimension);

    //     for (IndexType k = 0; k < number_of_points; k++)
    //     {
    //         const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

    //         Hessian(0, 0) += rDDN_DDe(k, 0)*coords[0];
    //         Hessian(0, 1) += rDDN_DDe(k, 2)*coords[0];
    //         Hessian(0, 2) += rDDN_DDe(k, 1)*coords[0];

    //         Hessian(1, 0) += rDDN_DDe(k, 0)*coords[1];
    //         Hessian(1, 1) += rDDN_DDe(k, 2)*coords[1];
    //         Hessian(1, 2) += rDDN_DDe(k, 1)*coords[1];

    //         Hessian(2, 0) += rDDN_DDe(k, 0)*coords[2];
    //         Hessian(2, 1) += rDDN_DDe(k, 2)*coords[2];
    //         Hessian(2, 2) += rDDN_DDe(k, 1)*coords[2];
    //     }
    // }

    // void Shell6pElement::CalculateSecondDerivativesOfBaseVectors(
    //     const Matrix& rDDDN_DDDe,
    //     array_1d<double, 3>& rDDa1_DD11,
    //     array_1d<double, 3>& rDDa1_DD12,
    //     array_1d<double, 3>& rDDa2_DD21,
    //     array_1d<double, 3>& rDDa2_DD22) const
    // {
    //     const std::size_t number_of_points = GetGeometry().size();
    
    //     for (IndexType k = 0; k < number_of_points; k++)
    //     {
    //         const array_1d<double, 3> coords = GetGeometry()[k].Coordinates();

    //         rDDa1_DD11 += rDDDN_DDDe(k, 0) * coords;
    //         rDDa1_DD12 += rDDDN_DDDe(k, 1) * coords;
    //         rDDa2_DD21 += rDDDN_DDDe(k, 2) * coords;
    //         rDDa2_DD22 += rDDDN_DDDe(k, 3) * coords;
    //     }
    // }

    ///@}

} // Namespace Kratos


