
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
#include "custom_conditions/gap_sbm_enhanced_load_solid_condition.h"
#include "geometries/coupling_geometry.h"

namespace Kratos
{

void GapSbmEnhancedLoadSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}


void GapSbmEnhancedLoadSolidCondition::InitializeMaterial()
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

void GapSbmEnhancedLoadSolidCondition::InitializeMemberVariables()
{
    InitializeMemberVariables(GetGeometry());
}

void GapSbmEnhancedLoadSolidCondition::InitializeMemberVariables(
    const GeometryType& r_geometry)
{
    // // Compute class memeber variables

    const auto& r_projected_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3)
        << "GapSbmEnhancedLoadSolidCondition only supports 2D and 3D conditions, but the current dimension is "
        << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
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
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    double integration_weight = r_integration_points[0].Weight();
    if (mDim == 2) {
        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        integration_weight *= thickness;
    }

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void GapSbmEnhancedLoadSolidCondition::InitializeSbmMemberVariables()
{
    auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];

    mDistanceVectorGap.resize(3);
    noalias(mDistanceVectorGap) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();

    mpSkinProjectionNode = r_geometry.GetValue(PROJECTION_NODE).get();

    mTrueNormal = mpSkinProjectionNode->GetValue(NORMAL);
    const std::string loop_identifier =
        mpSkinProjectionNode->GetValue(IDENTIFIER);

    if (loop_identifier == "inner") {
        mTrueNormal = -mTrueNormal;
    }

    const double true_normal_norm = norm_2(mTrueNormal);
    KRATOS_ERROR_IF(
        true_normal_norm <= std::numeric_limits<double>::epsilon())
        << Info() << " has a zero NORMAL on its PROJECTION_NODE.\n";
    mTrueNormal /= true_normal_norm;
        
    mDistanceVectorSkin.resize(3);
    noalias(mDistanceVectorSkin) = mpSkinProjectionNode->Coordinates() - r_surrogate_geometry.Center().Coordinates();
    // mDistanceVectorSkin = mDistanceVectorGap; //FIXME:

    this->SetValue(PROJECTION_NODE_COORDINATES, mpSkinProjectionNode->Coordinates());

    // dot product n dot n_tilde
    mTrueDotSurrogateNormal = inner_prod(mNormalPhysicalSpace, mTrueNormal);
    // mTrueDotSurrogateNormal = 1; //FIXME:

    const Point&  p_true = r_geometry.Center();            // true boundary
    const Point&  p_sur  = r_surrogate_geometry.Center();  // surrogate
}

void GapSbmEnhancedLoadSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetValue(NEIGHBOUR_GEOMETRIES)[0]->size() * mDim;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    const auto& r_surrogate_geometry =
        *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    std::vector<Matrix> taylor_derivatives(mBasisFunctionsOrder);
    for (IndexType order = 1;
         order <= mBasisFunctionsOrder;
         ++order) {
        taylor_derivatives[order - 1] =
            r_surrogate_geometry.ShapeFunctionDerivatives(
                order,
                0,
                this->GetIntegrationMethod());
    }
    mpTaylorDerivativeCache = &taylor_derivatives;
    try {
        CalculateLeftHandSide(
            rLeftHandSideMatrix,
            rCurrentProcessInfo);
        CalculateRightHandSide(
            rRightHandSideVector,
            rCurrentProcessInfo);
    } catch (...) {
        mpTaylorDerivativeCache = nullptr;
        throw;
    }
    mpTaylorDerivativeCache = nullptr;

    KRATOS_CATCH("")
}

void GapSbmEnhancedLoadSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_boundary_geometry = GetBoundaryGeometry();
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

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_boundary_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_boundary_sum, grad_N_boundary_sum);

    // compute stress Taylor expansion on skin
    Matrix grad_N_true_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_true_sum_transposed, mDistanceVectorSkin);
    Matrix grad_N_true_sum = trans(grad_N_true_sum_transposed);

    Matrix B_true_sum = ZeroMatrix(strain_size,mat_size);
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
    Matrix surrogate_traction_tangent(3, mat_size);
    Matrix true_traction_tangent(3, mat_size);
    Vector traction_column;
    for (IndexType column_index = 0;
         column_index < mat_size;
         ++column_index) {
        const Vector boundary_stress_column =
            column(DB_boundary_sum, column_index);
        CalculateTraction(
            boundary_stress_column,
            mNormalPhysicalSpace,
            traction_column);
        for (IndexType component = 0; component < mDim; ++component) {
            surrogate_traction_tangent(component, column_index) =
                traction_column[component];
        }

        const Vector true_stress_column =
            column(DB_true_sum, column_index);
        CalculateTraction(
            true_stress_column,
            mTrueNormal,
            traction_column);
        for (IndexType component = 0; component < mDim; ++component) {
            true_traction_tangent(component, column_index) =
                traction_column[component];
        }
    }

    // ASSEMBLE
    //-----------------------------------------------------
    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType j = 0; j < number_of_control_points; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const IndexType iglob = mDim*i+idim;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const IndexType jglob = mDim*j+jdim;

                    // FLUX 
                    // [sigma(u) \dot n_tilde] * (-w )
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -=
                        N_boundary_sum_vec(i) *
                        surrogate_traction_tangent(idim, jglob) *
                        integration_weight;
                    
                    // // SBM TERM
                    // // [E(sigma(u)) \dot n] (n*n_tilde) * (-w)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) +=
                        N_boundary_sum_vec(i) *
                        true_traction_tangent(idim, jglob) *
                        mTrueDotSurrogateNormal * integration_weight;
                }

            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmEnhancedLoadSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_boundary_geometry = GetBoundaryGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

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

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_boundary_sum = ZeroMatrix(strain_size, mat_size);
    CalculateB(B_boundary_sum, grad_N_boundary_sum);

    Matrix grad_N_true_sum_transposed(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_true_sum_transposed, mDistanceVectorSkin);
    Matrix grad_N_true_sum = trans(grad_N_true_sum_transposed);

    Matrix B_true_sum = ZeroMatrix(strain_size, mat_size);
    CalculateB(B_true_sum, grad_N_true_sum);

    Vector displacement_coefficients(mat_size);
    GetSolutionCoefficientVector(displacement_coefficients);

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

    Vector normal_stress_boundary;
    CalculateTraction(stress_boundary, mNormalPhysicalSpace, normal_stress_boundary);

    Vector normal_stress_true;
    CalculateTraction(stress_true, mTrueNormal, normal_stress_true);

    // retrieve external data
    
    Vector g_N = this->GetValue(FORCE);

    // debug manufactured sol
    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);
    const Vector projection_node_coordinates = mpSkinProjectionNode->Coordinates();
    const double x = projection_node_coordinates[0];
    const double y = projection_node_coordinates[1];
    const double z = projection_node_coordinates[2];

    // Kratos::array_1d<double, 3UL> projection_node_normal = mNormalPhysicalSpace;  //DEBUG

    Kratos::array_1d<double, 3UL> projection_node_normal = mTrueNormal;  

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


    for (IndexType i = 0; i < number_of_control_points; ++i) {
        for (IndexType idim = 0; idim < mDim; ++idim) {
            const IndexType iglob = mDim * i + idim;

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

    int GapSbmEnhancedLoadSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void GapSbmEnhancedLoadSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
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

    void GapSbmEnhancedLoadSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
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


    void GapSbmEnhancedLoadSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
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

    void GapSbmEnhancedLoadSolidCondition::CalculateB(
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
                rB(0, index) = r_DN_DX(i, 0);
                rB(1, index + 1) = r_DN_DX(i, 1);
                rB(2, index) = r_DN_DX(i, 1);
                rB(2, index + 1) = r_DN_DX(i, 0);
            } else {
                rB(0, index) = r_DN_DX(i, 0);
                rB(1, index + 1) = r_DN_DX(i, 1);
                rB(2, index + 2) = r_DN_DX(i, 2);
                rB(3, index) = r_DN_DX(i, 1);
                rB(3, index + 1) = r_DN_DX(i, 0);
                rB(4, index + 1) = r_DN_DX(i, 2);
                rB(4, index + 2) = r_DN_DX(i, 1);
                rB(5, index) = r_DN_DX(i, 2);
                rB(5, index + 2) = r_DN_DX(i, 0);
            }
        }
    }

    void GapSbmEnhancedLoadSolidCondition::CalculateTraction(
        const Vector& rStressVector,
        const array_1d<double, 3>& rNormal,
        Vector& rTraction) const
    {
        if (rTraction.size() != 3) {
            rTraction.resize(3, false);
        }
        noalias(rTraction) = ZeroVector(3);

        if (mDim == 2) {
            rTraction[0] = rStressVector[0] * rNormal[0] + rStressVector[2] * rNormal[1];
            rTraction[1] = rStressVector[2] * rNormal[0] + rStressVector[1] * rNormal[1];
        } else {
            rTraction[0] = rStressVector[0] * rNormal[0] + rStressVector[3] * rNormal[1] + rStressVector[5] * rNormal[2];
            rTraction[1] = rStressVector[3] * rNormal[0] + rStressVector[1] * rNormal[1] + rStressVector[4] * rNormal[2];
            rTraction[2] = rStressVector[5] * rNormal[0] + rStressVector[4] * rNormal[1] + rStressVector[2] * rNormal[2];
        }
    }

    void GapSbmEnhancedLoadSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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


    void GapSbmEnhancedLoadSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * mDim;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_boundary_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_boundary_sum_transposed, mDistanceVectorGap);
        Matrix grad_N_boundary_sum = trans(grad_N_boundary_sum_transposed);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B_boundary_sum = ZeroMatrix(strain_size,mat_size);
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
        Vector sigma_n;
        CalculateTraction(sigma, mNormalPhysicalSpace, sigma_n);
        sigma_n.resize(mDim, true);

        SetValue(NORMAL_STRESS, sigma_n);

        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, sigma[0]);
            SetValue(CAUCHY_STRESS_YY, sigma[1]);
            SetValue(CAUCHY_STRESS_XY, sigma[2]);
        } else {
            SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(sigma));
        }
        // //---------------------
    }

void GapSbmEnhancedLoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmEnhancedLoadSolidCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec, const Vector& rDistanceVector)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType number_of_control_points = r_geometry.PointsNumber();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_control_points)
    {
        H_sum_vec = ZeroVector(number_of_control_points);
    }

    std::vector<Matrix> local_shape_function_derivatives;
    if (mpTaylorDerivativeCache == nullptr) {
        local_shape_function_derivatives.resize(mBasisFunctionsOrder);
        for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
            local_shape_function_derivatives[n-1] =
                r_geometry.ShapeFunctionDerivatives(
                    n, 0, this->GetIntegrationMethod());
        }
    }
    const auto& shape_function_derivatives =
        mpTaylorDerivativeCache != nullptr
            ? *mpTaylorDerivativeCache
            : local_shape_function_derivatives;
    
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                const Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
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
                const Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {
                        
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

void GapSbmEnhancedLoadSolidCondition::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum, const Vector& rDistanceVector)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const SizeType number_of_control_points = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    std::vector<Matrix> local_shape_function_derivatives;
    if (mpTaylorDerivativeCache == nullptr) {
        local_shape_function_derivatives.resize(mBasisFunctionsOrder);
        for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
            local_shape_function_derivatives[n-1] =
                r_geometry.ShapeFunctionDerivatives(
                    n, 0, this->GetIntegrationMethod());
        }
    }
    const auto& shape_function_derivatives =
        mpTaylorDerivativeCache != nullptr
            ? *mpTaylorDerivativeCache
            : local_shape_function_derivatives;

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
                const Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
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
                const Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                IndexType countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {

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
double GapSbmEnhancedLoadSolidCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmEnhancedLoadSolidCondition::ComputeTaylorTerm3D(
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

GapSbmEnhancedLoadSolidConditionBatched::
    GapSbmEnhancedLoadSolidConditionBatched(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
    CompactQuadratureGeometries();
}

GapSbmEnhancedLoadSolidConditionBatched::
    GapSbmEnhancedLoadSolidConditionBatched(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
    CompactQuadratureGeometries();
}

Condition::Pointer GapSbmEnhancedLoadSolidConditionBatched::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<
        GapSbmEnhancedLoadSolidConditionBatched>(
            NewId, pGeom, pProperties);
}

Condition::Pointer GapSbmEnhancedLoadSolidConditionBatched::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    std::vector<GeometryType::Pointer> geometry_parts(
        NumberOfQuadraturePoints(),
        GetGeometry().pGetGeometryPart(0));
    auto p_geometry = Kratos::make_shared<CouplingGeometry<Node>>(
        std::move(geometry_parts));
    auto p_condition = Kratos::make_intrusive<
        GapSbmEnhancedLoadSolidConditionBatched>(
            NewId, p_geometry, pProperties);
    p_condition->mProjectionNodes = mProjectionNodes;
    p_condition->mQuadraturePointReferenceWeights =
        mQuadraturePointReferenceWeights;
    p_condition->mQuadraturePointCoordinates =
        mQuadraturePointCoordinates;
    p_condition->mQuadraturePointNormals = mQuadraturePointNormals;
    p_condition->mQuadraturePointForces = mQuadraturePointForces;
    p_condition->mHasQuadraturePointForce =
        mHasQuadraturePointForce;
    return p_condition;
}

GeometryData::IntegrationMethod
GapSbmEnhancedLoadSolidConditionBatched::GetIntegrationMethod() const
{
    return NumberOfQuadraturePoints() > 0
        ? GetRepresentativeGeometry().GetDefaultIntegrationMethod()
        : BaseType::GetIntegrationMethod();
}

void GapSbmEnhancedLoadSolidConditionBatched::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << Info() << " has no quadrature points.\n";
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_GEOMETRIES) &&
                        GetValue(NEIGHBOUR_GEOMETRIES).size() == 1)
        << Info() << " requires one neighbour geometry.\n";
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CONSTITUTIVE_LAW) &&
                        GetProperties()[CONSTITUTIVE_LAW] != nullptr)
        << Info() << " requires a constitutive law.\n";

    InitializeMemberVariables(GetRepresentativeGeometry());
    const std::size_t number_of_points = NumberOfQuadraturePoints();
    mQuadraturePointWeights.resize(number_of_points, false);
    mConstitutiveLawVector.resize(number_of_points);
    const auto& r_shape_values =
        GetRepresentativeGeometry().ShapeFunctionsValues(
            GetIntegrationMethod());
    double total_weight = 0.0;
    for (std::size_t point_index = 0;
         point_index < number_of_points;
         ++point_index) {
        double weight = mQuadraturePointReferenceWeights[point_index];
        if (mDim == 2) {
            weight *= GetProperties().Has(THICKNESS)
                ? GetProperties()[THICKNESS]
                : 1.0;
        }
        mQuadraturePointWeights[point_index] = weight;
        total_weight += weight;
        mConstitutiveLawVector[point_index] =
            GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[point_index]->InitializeMaterial(
            GetProperties(),
            GetRepresentativeGeometry(),
            row(r_shape_values, 0));
    }
    SetValue(INTEGRATION_WEIGHT, total_weight);
    SetCurrentQuadraturePoint(0);
}

void GapSbmEnhancedLoadSolidConditionBatched::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix,
        &rRightHandSideVector,
        rCurrentProcessInfo);
}

void GapSbmEnhancedLoadSolidConditionBatched::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix, nullptr, rCurrentProcessInfo);
}

void GapSbmEnhancedLoadSolidConditionBatched::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        nullptr, &rRightHandSideVector, rCurrentProcessInfo);
}

void GapSbmEnhancedLoadSolidConditionBatched::
    CalculateAllContributions(
        MatrixType* pLeftHandSideMatrix,
        VectorType* pRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    const std::size_t matrix_size =
        GetValue(NEIGHBOUR_GEOMETRIES)[0]->PointsNumber() * mDim;
    if (pLeftHandSideMatrix) {
        pLeftHandSideMatrix->resize(matrix_size, matrix_size, false);
        noalias(*pLeftHandSideMatrix) =
            ZeroMatrix(matrix_size, matrix_size);
    }
    if (pRightHandSideVector) {
        pRightHandSideVector->resize(matrix_size, false);
        noalias(*pRightHandSideVector) = ZeroVector(matrix_size);
    }

    const auto& r_surrogate_geometry =
        *GetValue(NEIGHBOUR_GEOMETRIES)[0];
    std::vector<Matrix> taylor_derivatives(mBasisFunctionsOrder);
    for (IndexType order = 1;
         order <= mBasisFunctionsOrder;
         ++order) {
        taylor_derivatives[order - 1] =
            r_surrogate_geometry.ShapeFunctionDerivatives(
                order, 0, GetIntegrationMethod());
    }
    mpTaylorDerivativeCache = &taylor_derivatives;
    try {
        Matrix point_lhs;
        Vector point_rhs;
        for (std::size_t point_index = 0;
             point_index < NumberOfQuadraturePoints();
             ++point_index) {
            SetCurrentQuadraturePoint(point_index);
            if (pLeftHandSideMatrix) {
                BaseType::CalculateLeftHandSide(
                    point_lhs, rCurrentProcessInfo);
                noalias(*pLeftHandSideMatrix) += point_lhs;
            }
            if (pRightHandSideVector) {
                BaseType::CalculateRightHandSide(
                    point_rhs, rCurrentProcessInfo);
                noalias(*pRightHandSideVector) += point_rhs;
            }
        }
    } catch (...) {
        mpTaylorDerivativeCache = nullptr;
        throw;
    }
    mpTaylorDerivativeCache = nullptr;
    mpConstitutiveLaw = mConstitutiveLawVector.front();
}

void GapSbmEnhancedLoadSolidConditionBatched::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        SetCurrentQuadraturePoint(point_index);
        BaseType::InitializeSolutionStep(rCurrentProcessInfo);
    }
}

void GapSbmEnhancedLoadSolidConditionBatched::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        SetCurrentQuadraturePoint(point_index);
        BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
    }
}

std::size_t GapSbmEnhancedLoadSolidConditionBatched::
    NumberOfQuadraturePoints() const
{
    return GetGeometry().NumberOfGeometryParts();
}

std::size_t GapSbmEnhancedLoadSolidConditionBatched::
    QuadraturePointsNumber() const
{
    return NumberOfQuadraturePoints();
}

const GapSbmEnhancedLoadSolidConditionBatched::GeometryType&
GapSbmEnhancedLoadSolidConditionBatched::
    GetRepresentativeGeometry() const
{
    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << Info() << " has no representative geometry.\n";
    return GetGeometry().GetGeometryPart(0);
}

const GapSbmEnhancedLoadSolidConditionBatched::GeometryType&
GapSbmEnhancedLoadSolidConditionBatched::GetBoundaryGeometry() const
{
    return GetRepresentativeGeometry();
}

void GapSbmEnhancedLoadSolidConditionBatched::
    CompactQuadratureGeometries()
{
    const std::size_t number_of_points =
        GetGeometry().NumberOfGeometryParts();
    if (number_of_points == 0 ||
        mQuadraturePointCoordinates.size1() == number_of_points) {
        return;
    }
    mProjectionNodes.resize(number_of_points);
    mQuadraturePointReferenceWeights.resize(number_of_points, false);
    mQuadraturePointCoordinates.resize(number_of_points, 3, false);
    mQuadraturePointNormals.resize(number_of_points, 3, false);
    mQuadraturePointForces.resize(number_of_points, 3, false);
    noalias(mQuadraturePointForces) = ZeroMatrix(number_of_points, 3);
    mHasQuadraturePointForce.assign(number_of_points, false);
    auto p_representative = GetGeometry().pGetGeometryPart(0);
    for (std::size_t point_index = 0;
         point_index < number_of_points;
         ++point_index) {
        auto p_geometry = GetGeometry().pGetGeometryPart(point_index);
        KRATOS_ERROR_IF_NOT(p_geometry->Has(PROJECTION_NODE))
            << "Missing PROJECTION_NODE at point " << point_index
            << " of " << Info() << ".\n";
        mProjectionNodes[point_index] =
            p_geometry->GetValue(PROJECTION_NODE);
        const auto method = p_geometry->GetDefaultIntegrationMethod();
        const auto& r_points = p_geometry->IntegrationPoints(method);
        KRATOS_ERROR_IF(r_points.size() != 1)
            << Info() << " requires one integration point per geometry.\n";
        mQuadraturePointReferenceWeights[point_index] =
            r_points.front().Weight();
        const auto center = p_geometry->Center();
        auto normal = p_geometry->Normal(0, method);
        const double normal_norm = norm_2(normal);
        KRATOS_ERROR_IF(normal_norm <= 0.0)
            << "Zero normal at point " << point_index
            << " of " << Info() << ".\n";
        normal /= normal_norm;
        for (std::size_t component = 0; component < 3; ++component) {
            mQuadraturePointCoordinates(point_index, component) =
                center[component];
            mQuadraturePointNormals(point_index, component) =
                normal[component];
        }
        if (p_geometry->Has(FORCE)) {
            const auto& r_force = p_geometry->GetValue(FORCE);
            mHasQuadraturePointForce[point_index] = true;
            for (std::size_t component = 0; component < 3; ++component) {
                mQuadraturePointForces(point_index, component) =
                    r_force[component];
            }
        }
        if (point_index > 0) {
            GetGeometry().SetGeometryPart(point_index, p_representative);
        }
    }
}

void GapSbmEnhancedLoadSolidConditionBatched::
    SetCurrentQuadraturePoint(const std::size_t PointIndex)
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints())
        << "Invalid point index in " << Info() << ".\n";
    mpConstitutiveLaw = mConstitutiveLawVector[PointIndex];
    mpSkinProjectionNode = mProjectionNodes[PointIndex].get();
    mTrueNormal = mpSkinProjectionNode->GetValue(NORMAL);
    const std::string loop_identifier =
        mpSkinProjectionNode->GetValue(IDENTIFIER);
    if (loop_identifier == "inner") {
        mTrueNormal = -mTrueNormal;
    }
    const double true_normal_norm = norm_2(mTrueNormal);
    KRATOS_ERROR_IF(
        true_normal_norm <= std::numeric_limits<double>::epsilon())
        << Info() << " has a zero NORMAL on projection node "
        << mpSkinProjectionNode->Id() << " at quadrature point "
        << PointIndex << ".\n";
    mTrueNormal /= true_normal_norm;

    SetValue(INTEGRATION_WEIGHT, mQuadraturePointWeights[PointIndex]);
    if (mHasQuadraturePointForce[PointIndex]) {
        array_1d<double, 3> force = ZeroVector(3);
        for (std::size_t component = 0; component < 3; ++component) {
            force[component] =
                mQuadraturePointForces(PointIndex, component);
        }
        SetValue(FORCE, force);
    }
    const auto surrogate_center =
        GetValue(NEIGHBOUR_GEOMETRIES)[0]->Center();
    if (mDistanceVectorGap.size() != 3) {
        mDistanceVectorGap.resize(3, false);
        mDistanceVectorSkin.resize(3, false);
    }
    for (std::size_t component = 0; component < 3; ++component) {
        mNormalParameterSpace[component] =
            mQuadraturePointNormals(PointIndex, component);
        mNormalPhysicalSpace[component] =
            mQuadraturePointNormals(PointIndex, component);
        mDistanceVectorGap[component] =
            mQuadraturePointCoordinates(PointIndex, component) -
            surrogate_center[component];
        mDistanceVectorSkin[component] =
            mpSkinProjectionNode->Coordinates()[component] -
            surrogate_center[component];
    }
    mTrueDotSurrogateNormal =
        inner_prod(mNormalPhysicalSpace, mTrueNormal);
}

array_1d<double, 3>
GapSbmEnhancedLoadSolidConditionBatched::
    GetQuadraturePointCoordinates(const std::size_t PointIndex) const
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints())
        << "Invalid point index in " << Info() << ".\n";
    array_1d<double, 3> coordinates = ZeroVector(3);
    for (std::size_t component = 0; component < 3; ++component) {
        coordinates[component] =
            mQuadraturePointCoordinates(PointIndex, component);
    }
    return coordinates;
}

void GapSbmEnhancedLoadSolidConditionBatched::
    SetQuadraturePointForceComponent(
        const std::size_t PointIndex,
        const std::size_t ComponentIndex,
        const double Value)
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints() ||
                    ComponentIndex >= 3)
        << "Invalid force index in " << Info() << ".\n";
    mQuadraturePointForces(PointIndex, ComponentIndex) = Value;
    mHasQuadraturePointForce[PointIndex] = true;
}

void GapSbmEnhancedLoadSolidConditionBatched::save(
    Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(
        rSerializer, GapSbmEnhancedLoadSolidCondition);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.save("ProjectionNodes", mProjectionNodes);
    rSerializer.save("QuadraturePointReferenceWeights",
                     mQuadraturePointReferenceWeights);
    rSerializer.save("QuadraturePointWeights", mQuadraturePointWeights);
    rSerializer.save("QuadraturePointCoordinates",
                     mQuadraturePointCoordinates);
    rSerializer.save("QuadraturePointNormals", mQuadraturePointNormals);
    rSerializer.save("QuadraturePointForces", mQuadraturePointForces);
    rSerializer.save("HasQuadraturePointForce",
                     mHasQuadraturePointForce);
}

void GapSbmEnhancedLoadSolidConditionBatched::load(
    Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(
        rSerializer, GapSbmEnhancedLoadSolidCondition);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.load("ProjectionNodes", mProjectionNodes);
    rSerializer.load("QuadraturePointReferenceWeights",
                     mQuadraturePointReferenceWeights);
    rSerializer.load("QuadraturePointWeights", mQuadraturePointWeights);
    rSerializer.load("QuadraturePointCoordinates",
                     mQuadraturePointCoordinates);
    rSerializer.load("QuadraturePointNormals", mQuadraturePointNormals);
    rSerializer.load("QuadraturePointForces", mQuadraturePointForces);
    rSerializer.load("HasQuadraturePointForce",
                     mHasQuadraturePointForce);
}

} // Namespace Kratos
