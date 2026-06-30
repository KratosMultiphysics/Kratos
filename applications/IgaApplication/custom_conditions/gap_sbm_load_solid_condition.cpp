
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andea Gorgi
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/gap_sbm_load_solid_condition.h"

namespace Kratos
{

void GapSbmLoadSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}


void GapSbmLoadSolidCondition::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetSurrogateGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_control_points = r_geometry.size();
        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, N_sum_vec);
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void GapSbmLoadSolidCondition::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "GapSbmSolidCondition momentarily only supports 2D and 3D conditions, but the current dimension is" << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 3; 
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 2; 
    }

    // Compute the normals
    mNormalParameterSpace = r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    if (mDim == 2)
    {
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        const double integration_weight = r_integration_points[0].Weight()*thickness;
        SetValue(INTEGRATION_WEIGHT, integration_weight);
    }
    else //mDim == 3
    {
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        const double integration_weight = r_integration_points[0].Weight();
        SetValue(INTEGRATION_WEIGHT, integration_weight);        
    }
}

void GapSbmLoadSolidCondition::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();

    const Point&  p_true = r_geometry.Center();            // true boundary
    const Point&  p_sur  = r_surrogate_geometry.Center();  // surrogate

}

void GapSbmLoadSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const std::size_t mat_size = GetSurrogateGeometry().size() * mDim;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void GapSbmLoadSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

void GapSbmLoadSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_true_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    Vector g_N = this->GetValue(FORCE); 

    // Vector g_N = ZeroVector(3); 

    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);

    // 2D
    // const double x = r_true_geometry.Center().X();
    // const double y = r_true_geometry.Center().Y();

    // // // // cosinusoidal
    // g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * mNormalPhysicalSpace[0]; 
    // g_N[1] = E/(1-nu)*(sin(x)*sinh(y)) * mNormalPhysicalSpace[1]; 

    // g_N[0] = E/(1-nu*nu) * mNormalPhysicalSpace[0] +  E/2/(1+nu) * mNormalPhysicalSpace[1]; 
    // g_N[1] = E/2/(1+nu) * mNormalPhysicalSpace[0] + E*nu/(1-nu*nu)* mNormalPhysicalSpace[1]; 

    //3D 
    // const double x = r_true_geometry.Center().X();
    // const double y = r_true_geometry.Center().Y();
    // const double z = r_true_geometry.Center().Z();

    const Vector projection_node_coordinates = r_true_geometry.GetValue(PROJECTION_NODE)->Coordinates();
    const double x = projection_node_coordinates[0];
    const double y = projection_node_coordinates[1];
    const double z = projection_node_coordinates[2];

    Kratos::array_1d<double, 3UL> projection_node_normal = -r_true_geometry.GetValue(PROJECTION_NODE)->GetValue(NORMAL);

    projection_node_normal = mNormalPhysicalSpace;  //DEBUG

    const double c_vol = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    const double c_shear = E / (2.0 * (1.0 + nu));

    const double div_u =
    -std::sin(x) * std::sinh(y) * std::cosh(z)
    -std::sin(y) * std::sinh(z) * std::cosh(x)
    -std::sin(z) * std::sinh(x) * std::cosh(y);

    const double sigma_xx =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(x) * std::sinh(y) * std::cosh(z);

    const double sigma_yy =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(y) * std::sinh(z) * std::cosh(x);

    const double sigma_zz =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(z) * std::sinh(x) * std::cosh(y);

    const double sigma_xy =
        c_shear * (
            std::cos(x) * std::cosh(y) * std::cosh(z)
        + std::cos(y) * std::sinh(z) * std::sinh(x)
        );

    const double sigma_yz =
        c_shear * (
            std::cos(y) * std::cosh(z) * std::cosh(x)
        + std::cos(z) * std::sinh(x) * std::sinh(y)
        );

    const double sigma_xz =
        c_shear * (
            std::cos(x) * std::sinh(y) * std::sinh(z)
        + std::cos(z) * std::cosh(x) * std::cosh(y)
        );

    // lienar case

    // const double sigma_diag = E / (1.0 - 2.0 * nu);
    // const double sigma_off  = E / (1.0 + nu);

    // const double sigma_xx  = sigma_diag; // sigma_xx
    // const double sigma_yy  = sigma_diag; // sigma_yy
    // const double sigma_zz  = sigma_diag; // sigma_zz
    // const double sigma_xy  = sigma_off;  // sigma_xy
    // const double sigma_yz  = sigma_off;  // sigma_yz
    // const double sigma_xz  = sigma_off;  // sigma_xz

    array_1d<double, 3> traction;

    traction[0] =
        sigma_xx   * projection_node_normal[0]
        + sigma_xy * projection_node_normal[1]
        + sigma_xz * projection_node_normal[2];
    
    traction[1] =
        sigma_xy   * projection_node_normal[0]
        + sigma_yy * projection_node_normal[1]
        + sigma_yz * projection_node_normal[2];
    
    traction[2] =
        sigma_xz   * projection_node_normal[0]
        + sigma_yz * projection_node_normal[1]
        + sigma_zz * projection_node_normal[2];

    g_N = traction; //FIXME::



    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType zdim = 0; zdim < mDim; zdim++) {
            
            rRightHandSideVector[mDim*i+zdim] += N_sum_vec(i)*g_N[zdim] * integration_weight;

        }
    }

    KRATOS_CATCH("")
}


    void GapSbmLoadSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != mDim * number_of_control_points)
            rResult.resize(mDim * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * mDim;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            if (mDim == 3) {
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            }
        }
    }

    void GapSbmLoadSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(mDim * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            if (mDim == 3) {
                rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            }
        }
    };


    void GapSbmLoadSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * mDim;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * mDim;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            if (mDim == 3) {
                rValues[index + 2] = displacement[2];
            }
        }
    }

    void GapSbmLoadSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * mDim;
        const SizeType strain_size = (mDim == 2) ? 3 : 6;
    
        if (rB.size1() != strain_size || rB.size2() != mat_size)
            rB.resize(strain_size, mat_size, false);
    
        noalias(rB) = ZeroMatrix(strain_size, mat_size);
    
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * mDim;
    
            if (mDim == 2) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy
    
                rB(2, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(2, index + 1) = r_DN_DX(i, 0);
            }
            else if (mDim == 3) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy
                rB(2, index + 2) = r_DN_DX(i, 2); // ezz
    
                rB(3, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(3, index + 1) = r_DN_DX(i, 0);
    
                rB(4, index + 1) = r_DN_DX(i, 2); // gamma_yz
                rB(4, index + 2) = r_DN_DX(i, 1);
    
                rB(5, index + 0) = r_DN_DX(i, 2); // gamma_xz
                rB(5, index + 2) = r_DN_DX(i, 0);
            }
        }
    }

    void GapSbmLoadSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
    {
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
    }


    void GapSbmLoadSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * mDim;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);
        
        const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B_sum = ZeroMatrix(strain_size,mat_size);
        CalculateB(B_sum, grad_N_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

        const Vector sigma = values_true.GetStressVector();
        Vector sigma_n = ZeroVector(mDim);

        if (mDim == 2) {
            sigma_n[0] = sigma[0] * mNormalPhysicalSpace[0]
                    + sigma[2] * mNormalPhysicalSpace[1];

            sigma_n[1] = sigma[2] * mNormalPhysicalSpace[0]
                    + sigma[1] * mNormalPhysicalSpace[1];
        } else {
            sigma_n[0] = sigma[0] * mNormalPhysicalSpace[0]
                    + sigma[3] * mNormalPhysicalSpace[1]
                    + sigma[5] * mNormalPhysicalSpace[2];

            sigma_n[1] = sigma[3] * mNormalPhysicalSpace[0]
                    + sigma[1] * mNormalPhysicalSpace[1]
                    + sigma[4] * mNormalPhysicalSpace[2];

            sigma_n[2] = sigma[5] * mNormalPhysicalSpace[0]
                    + sigma[4] * mNormalPhysicalSpace[1]
                    + sigma[2] * mNormalPhysicalSpace[2];
        }
        SetValue(NORMAL_STRESS, sigma_n);

        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, sigma[0]);
            SetValue(CAUCHY_STRESS_YY, sigma[1]);
            SetValue(CAUCHY_STRESS_XY, sigma[2]);
        } else {
            SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(sigma));
        }
        // //---------------------

        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);

        array_1d<double, 3> current_displacement = ZeroVector(3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            current_displacement[0] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i];
            current_displacement[1] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i + 1];
            if (mDim == 3) {
                current_displacement[2] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i + 2];
            }
        }
        SetValue(DISPLACEMENT, current_displacement);
    }

void GapSbmLoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmLoadSolidCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_control_points)
    {
        H_sum_vec = ZeroVector(number_of_control_points);
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {
                        
                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i,countDerivativeId); 

                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}

void GapSbmLoadSolidCondition::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points)
    {
        grad_H_sum.resize(3, number_of_control_points);
    }

    // Neumann (Taylor expansion of the gradient)
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each control point
        double H_taylor_term_X = 0.0; 
        double H_taylor_term_Y = 0.0; 
        double H_taylor_term_Z = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                IndexType countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {

                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = shapeFunctionDerivatives(i,countDerivativeId); 
                        
                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x-1, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y-1, mDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z-1);
                        }     
                        countDerivativeId++;
                    }
                }
            }
        }
        grad_H_sum(0,i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1,i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2,i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else 
            grad_H_sum(2,i) = 0;
    }    
}

// Function to compute a single term in the Taylor expansion
double GapSbmLoadSolidCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmLoadSolidCondition::ComputeTaylorTerm3D(
    const double derivative, 
    const double dx, 
    const IndexType k_x, 
    const double dy, 
    const IndexType k_y, 
    const double dz, 
    const IndexType k_z)
{   
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));    
}

} // Namespace Kratos
