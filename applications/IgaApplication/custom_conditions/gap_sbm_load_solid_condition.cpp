
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
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}


void GapSbmLoadSolidCondition::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void GapSbmLoadSolidCondition::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2) << "GapSbmLoadSolidCondition momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    mBasisFunctionsOrder *= 2; 

    // Compute the normals
    mNormalParameterSpace = r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    // const double integration_weight = r_integration_points[0].Weight() * std::abs(detJ0) * thickness;

    const double integration_weight = r_integration_points[0].Weight()*thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void GapSbmLoadSolidCondition::InitializeSbmMemberVariables()
{
    auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];

    mDistanceVectorGap.resize(3);
    noalias(mDistanceVectorGap) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();

    mpSkinProjectionNode = &r_geometry.GetValue(NEIGHBOUR_NODES)[0];

    // Vector temp = r_geometry.GetValue(TEMP_INTERPOLATION_NODES_ID);
    // this->SetValue(TEMP_INTERPOLATION_NODES_ID, temp);

    mTrueNormal = mpSkinProjectionNode->GetValue(NORMAL);
    std::string loopIdentifier = mpSkinProjectionNode->GetValue(IDENTIFIER);

    if (loopIdentifier == "inner")
        mTrueNormal = -mTrueNormal;
    
    mTrueNormal /= norm_2(mTrueNormal);
        
    mDistanceVectorSkin.resize(3);
    noalias(mDistanceVectorSkin) = mpSkinProjectionNode->Coordinates() - r_surrogate_geometry.Center().Coordinates();
    // mDistanceVectorSkin = mDistanceVectorGap; //FIXME:

    this->SetValue(PROJECTION_NODE_COORDINATES, mpSkinProjectionNode->Coordinates());

    // dot product n dot n_tilde
    mTrueDotSurrogateNormal = inner_prod(mNormalPhysicalSpace, mTrueNormal);
    // mTrueDotSurrogateNormal = 1; //FIXME:

    const Point&  p_true = r_geometry.Center();            // true boundary
    const Point&  p_sur  = r_surrogate_geometry.Center();  // surrogate

    // std::ofstream out("centers.txt", std::ios::app);       // append mode
    // out << std::setprecision(15)                           // full precision
    //     << p_true.X() << ' ' << p_true.Y() << ' ' << p_true.Z() << ' '
    //     << p_sur .X() << ' ' << p_sur .Y() << ' ' << p_sur .Z() << '\n';
}

void GapSbmLoadSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetValue(NEIGHBOUR_GEOMETRIES)[0]->size() * 2;

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
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_boundary_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_boundary_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_boundary_sum_vec, mDistanceVectorGap);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_boundary_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_boundary_sum_transposed, mDistanceVectorGap);
    Matrix grad_N_boundary_sum = trans(grad_N_boundary_sum_transposed);

    Matrix B_boundary_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_boundary_sum, grad_N_boundary_sum);

    // compute stress Taylor expansion on skin
    Matrix grad_N_true_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_true_sum_transposed, mDistanceVectorSkin);
    Matrix grad_N_true_sum = trans(grad_N_true_sum_transposed);

    Matrix B_true_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_true_sum, grad_N_true_sum);


    // obtain the tangent constitutive matrix at the boundary position
    
    ConstitutiveLaw::Parameters values_boundary(r_boundary_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_boundary_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_boundary_coefficient_vector);
    Vector old_strain_on_boundary = prod(B_boundary_sum, old_displacement_boundary_coefficient_vector);

    const SizeType strain_size_boundary = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_boundary(strain_size_boundary);
    ApplyConstitutiveLaw(mat_size, old_strain_on_boundary, values_boundary, this_constitutive_variables_boundary);

    const Matrix& r_D_on_boundary = values_boundary.GetConstitutiveMatrix();

    Matrix DB_boundary_sum = prod(r_D_on_boundary, B_boundary_sum);


    // obtain the tangent constitutive matrix at the true position
    // TODO: for damage pass the true value geometry (shape functions evaluate on the exact true location)
    ConstitutiveLaw::Parameters values_true(r_boundary_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_strain_on_true = prod(B_true_sum, old_displacement_boundary_coefficient_vector);

    const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();

    Matrix DB_true_sum = prod(r_D_on_true, B_true_sum);

    // ASSEMBLE
    //-----------------------------------------------------
    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType j = 0; j < number_of_control_points; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    // FLUX 
                    // [sigma(u) \dot n_tilde] * (-w )
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= N_boundary_sum_vec(i)*(DB_boundary_sum(id1, jglob)* mNormalPhysicalSpace[0] + DB_boundary_sum(id2, jglob)* mNormalPhysicalSpace[1]) * integration_weight;
                    
                    // // SBM TERM
                    // // [E(sigma(u)) \dot n] (n*n_tilde) * (-w)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += N_boundary_sum_vec(i)*(DB_true_sum(id1, jglob)* mTrueNormal[0] + DB_true_sum(id2, jglob)* mTrueNormal[1]) 
                                                        * mTrueDotSurrogateNormal * integration_weight;
                }

            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmLoadSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_boundary_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    KRATOS_ERROR_IF(mDim != 2) << "GapSbmLoadSolidCondition momentarily only supports 2D conditions, but the current dimension is " << mDim << std::endl;

    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    Vector N_boundary_sum_vec;
    ComputeTaylorExpansionContribution(N_boundary_sum_vec, mDistanceVectorGap);

    Matrix grad_N_boundary_sum_transposed(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_boundary_sum_transposed, mDistanceVectorGap);
    Matrix grad_N_boundary_sum = trans(grad_N_boundary_sum_transposed);

    Matrix B_boundary_sum = ZeroMatrix(mDim, mat_size);
    CalculateB(B_boundary_sum, grad_N_boundary_sum);

    Matrix grad_N_true_sum_transposed(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_true_sum_transposed, mDistanceVectorSkin);
    Matrix grad_N_true_sum = trans(grad_N_true_sum_transposed);

    Matrix B_true_sum = ZeroMatrix(mDim, mat_size);
    CalculateB(B_true_sum, grad_N_true_sum);

    Vector displacement_coefficients(mat_size);
    GetSolutionCoefficientVector(displacement_coefficients);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();

    // Surrogate (gap) side response
    ConstitutiveLaw::Parameters values_boundary(r_boundary_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_boundary(strain_size);
    Vector strain_boundary = prod(B_boundary_sum, displacement_coefficients);
    ApplyConstitutiveLaw(mat_size, strain_boundary, values_boundary, constitutive_variables_boundary);
    const Vector& stress_boundary = values_boundary.GetStressVector();

    // True boundary response (projection)
    // TODO: damage materials
    ConstitutiveLaw::Parameters values_true(r_boundary_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_true(strain_size);
    Vector strain_true = prod(B_true_sum, displacement_coefficients);
    ApplyConstitutiveLaw(mat_size, strain_true, values_true, constitutive_variables_true);
    const Vector& stress_true = values_true.GetStressVector();

    array_1d<double,3> normal_stress_boundary;
    normal_stress_boundary.clear();
    normal_stress_boundary[0] = stress_boundary[0] * mNormalPhysicalSpace[0] + stress_boundary[2] * mNormalPhysicalSpace[1];
    normal_stress_boundary[1] = stress_boundary[2] * mNormalPhysicalSpace[0] + stress_boundary[1] * mNormalPhysicalSpace[1];

    array_1d<double,3> normal_stress_true;
    normal_stress_true.clear();
    normal_stress_true[0] = stress_true[0] * mTrueNormal[0] + stress_true[2] * mTrueNormal[1];
    normal_stress_true[1] = stress_true[2] * mTrueNormal[0] + stress_true[1] * mTrueNormal[1];

    // retrieve external data

    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);
    Vector g_N = ZeroVector(3);

    const double x = mpSkinProjectionNode->X();
    const double y = mpSkinProjectionNode->Y();
    // // cosinusoidal
    g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[0]; 
    g_N[1] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[1]; 

    // g_N = mpSkinProjectionNode->GetValue(FORCE);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        for (IndexType idim = 0; idim < mDim; ++idim) {
            const IndexType iglob = 2 * i + idim;

            double rhs_contribution = 0.0;

            // Flux term (traction evaluated on surrogate boundary)
            rhs_contribution += N_boundary_sum_vec(i) * normal_stress_boundary[idim] * integration_weight;

            // SBM correction: contribution evaluated on the true boundary
            rhs_contribution -= N_boundary_sum_vec(i) * normal_stress_true[idim] * mTrueDotSurrogateNormal * integration_weight;

            // External load applied on the true boundary
            rhs_contribution += N_boundary_sum_vec(i) * g_N[idim] * mTrueDotSurrogateNormal * integration_weight;

            rRightHandSideVector[iglob] += rhs_contribution;
        }
    }

    KRATOS_CATCH("")
}

    int GapSbmLoadSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void GapSbmLoadSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void GapSbmLoadSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void GapSbmLoadSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void GapSbmLoadSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
            rB(1, r) = r_DN_DX(kr,1) * dirr;
            rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
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
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_boundary_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_boundary_sum_transposed, mDistanceVectorGap);
        Matrix grad_N_boundary_sum = trans(grad_N_boundary_sum_transposed);

        Matrix B_boundary_sum = ZeroMatrix(mDim,mat_size);
        CalculateB(B_boundary_sum, grad_N_boundary_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_boundary(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_boundary_sum, old_displacement_coefficient_vector);

        const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_boundary, this_constitutive_variables_true);

        const Vector sigma = values_boundary.GetStressVector();
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*mNormalPhysicalSpace[0] + sigma[2]*mNormalPhysicalSpace[1];
        sigma_n[1] = sigma[2]*mNormalPhysicalSpace[0] + sigma[1]*mNormalPhysicalSpace[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
    }

void GapSbmLoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmLoadSolidCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec, Vector& rDistanceVector)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType number_of_control_points = r_geometry.PointsNumber();
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
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        
                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i,countDerivativeId); 

                        H_taylor_term += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}

void GapSbmLoadSolidCondition::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum, Vector& rDistanceVector)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType number_of_control_points = r_geometry.PointsNumber();
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
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                IndexType countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {

                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = shapeFunctionDerivatives(i,countDerivativeId); 
                        
                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x-1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y-1, rDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z-1);
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
