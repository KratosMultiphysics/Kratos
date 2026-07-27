//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Andrea Gorgi
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "geometries/coupling_geometry.h"

// Application includes
#include "custom_elements/gap_sbm_solid_element.h"

namespace Kratos
{

GapSbmSolidElement::GapSbmSolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

GapSbmSolidElement::GapSbmSolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer GapSbmSolidElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidElement>(NewId, GetSurrogateGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer GapSbmSolidElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidElement>(NewId, pGeom, pProperties);
}

// Deconstructor

GapSbmSolidElement::~GapSbmSolidElement()
{
}

void GapSbmSolidElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}


void GapSbmSolidElement::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetSurrogateGeometry();
        const Properties& r_properties = GetProperties();        
        const std::size_t number_of_control_points = r_geometry.size();
        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);
    
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, N_sum_vec);

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void GapSbmSolidElement::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "GapSbmSolidElement momentarily only supports 2D and 3D elements, but the current dimension is" << mDim << std::endl;
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 3; 
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 2; 
    }

    // calculate the integration weight
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    if (mDim == 2)
    {
        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        const double integration_weight = r_integration_points[0].Weight()*thickness;
        SetValue(INTEGRATION_WEIGHT, integration_weight);
    }
    if (mDim == 3)
    {
        const double integration_weight = r_integration_points[0].Weight();
        SetValue(INTEGRATION_WEIGHT, integration_weight);
    }
}

void GapSbmSolidElement::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();
}

void GapSbmSolidElement::CalculateLocalSystem(
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

void GapSbmSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const std::size_t mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    ConstitutiveLaw::Parameters values_true(r_surrogate_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);
    const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();

    noalias(rLeftHandSideMatrix) += integration_weight * prod(trans(B_sum), Matrix(prod(r_D_on_true, B_sum))); 

    KRATOS_CATCH("")
}

void GapSbmSolidElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const std::size_t mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    ConstitutiveLaw::Parameters values_true(r_surrogate_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);
    const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Vector& r_stress_vector_on_true = values_true.GetStressVector();

    //-----------------------------------------------------------------------------------
    Vector volume_force_local = this->GetValue(BODY_FORCE);
    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_control_points; ++i ) {
        const std::size_t index = mDim * i;

        for ( IndexType j = 0; j < mDim; ++j )
            rRightHandSideVector[index + j] += integration_weight * N_sum_vec[i] * volume_force_local[j];
    }

    // RHS = ExtForces - K*temp;
    noalias(rRightHandSideVector) -= integration_weight * prod(trans(B_sum), r_stress_vector_on_true); 

    KRATOS_CATCH("")
}

void GapSbmSolidElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_geometry.size();

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

void GapSbmSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(mDim * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        const auto& r_node = r_geometry[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        if (mDim > 2) {
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    }
};




int GapSbmSolidElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_surrogate_geometry.ShapeFunctionsLocalGradients(r_surrogate_geometry.GetDefaultIntegrationMethod());
    const std::size_t dimension = r_DN_De[0].size2();

    KRATOS_ERROR_IF(dimension != 2 && dimension != 3)
        << "GapSbmSolidElement only supports 2D and 3D elements, but the current dimension is "
        << dimension << std::endl;

    // Verify that the constitutive law exists
    if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if (dimension == 2 && this->GetProperties()[CONSTITUTIVE_LAW]->GetStrainSize() != 3)
        {
            KRATOS_ERROR << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
                         << this->Id() << std::endl;
        }
        else if (dimension == 3 && this->GetProperties()[CONSTITUTIVE_LAW]->GetStrainSize() != 6)
        {
            KRATOS_ERROR << "Wrong constitutive law used. This is a 3D element! Expected strain size is 6 (el id = ) "
                         << this->Id() << std::endl;
        }

        if (dimension == 2 && !this->GetProperties().Has(THICKNESS))
        {
            KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
    }
    // Intentionally left blank: Gap-SBM element bypasses default geometry/size checks.
    return 0;
}


Element::IntegrationMethod GapSbmSolidElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void GapSbmSolidElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

    //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * mDim;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);

        const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B_sum = ZeroMatrix(strain_size, mat_size);
        CalculateB(B_sum, grad_N_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

        const Vector sigma = values_true.GetStressVector();
    
        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, sigma[0]);
            SetValue(CAUCHY_STRESS_YY, sigma[1]);
            SetValue(CAUCHY_STRESS_XY, sigma[2]);
        }
        else if (mDim == 3) {
            Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(sigma);
            SetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
        }
        // //---------------------

        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);

        array_1d<double, 3> current_displacement = ZeroVector(3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * mDim;
            for (IndexType j = 0; j < mDim; ++j) {
                current_displacement[j] += N_sum_vec[i] * old_displacement_coefficient_vector[index + j];
            }
        }
        SetValue(DISPLACEMENT, current_displacement);
}

void GapSbmSolidElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}



void GapSbmSolidElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    if (rOutput.size() != r_integration_points.size())
    {
        rOutput.resize(r_integration_points.size());
    }

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
        KRATOS_WATCH(rVariable);
        KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}

void GapSbmSolidElement::CalculateOnIntegrationPoints(
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

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
            KRATOS_WATCH(rVariable);
            KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}

void GapSbmSolidElement::CalculateB(
    Matrix& rB, 
    Matrix& r_DN_DX) const
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
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

void GapSbmSolidElement::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * mDim;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = r_surrogate_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * mDim;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            if (mDim == 3) {
                rValues[index + 2] = displacement[2];
            }
        }
    }

void GapSbmSolidElement::ApplyConstitutiveLaw(std::size_t matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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

void GapSbmSolidElement::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = GetSurrogateGeometry();
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

void GapSbmSolidElement::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = GetSurrogateGeometry();
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
                        IndexType k_z = n - static_cast<IndexType>(k_x) - static_cast<IndexType>(k_y);
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
double GapSbmSolidElement::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmSolidElement::ComputeTaylorTerm3D(
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

GapSbmSolidElement::TaylorDerivativeData
GapSbmSolidElement::CalculateTaylorDerivativeData() const
{
    const auto& r_geometry = GetSurrogateGeometry();
    TaylorDerivativeData derivative_data(mBasisFunctionsOrder);
    for (IndexType order = 1; order <= mBasisFunctionsOrder; ++order) {
        derivative_data[order - 1] = r_geometry.ShapeFunctionDerivatives(
            order,
            0,
            this->GetIntegrationMethod());
    }
    return derivative_data;
}

void GapSbmSolidElement::EvaluateTaylorExpansion(
    const TaylorDerivativeData& rDerivativeData,
    Vector* pShapeFunctions,
    Matrix* pShapeFunctionGradients) const
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points =
        r_geometry.PointsNumber();
    KRATOS_ERROR_IF(rDerivativeData.size() != mBasisFunctionsOrder)
        << "Invalid Taylor derivative data in GapSbmSolidElement #"
        << Id() << ".\n";

    if (pShapeFunctions != nullptr) {
        if (pShapeFunctions->size() != number_of_control_points) {
            pShapeFunctions->resize(number_of_control_points, false);
        }
        noalias(*pShapeFunctions) = row(
            r_geometry.ShapeFunctionsValues(),
            0);
    }

    if (pShapeFunctionGradients != nullptr) {
        if (pShapeFunctionGradients->size1() != 3 ||
            pShapeFunctionGradients->size2() !=
                number_of_control_points) {
            pShapeFunctionGradients->resize(
                3,
                number_of_control_points,
                false);
        }
        noalias(*pShapeFunctionGradients) =
            ZeroMatrix(3, number_of_control_points);
        const auto& r_local_gradients =
            r_geometry.ShapeFunctionsLocalGradients(
                r_geometry.GetDefaultIntegrationMethod())[0];
        for (std::size_t i = 0; i < number_of_control_points; ++i) {
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                (*pShapeFunctionGradients)(component, i) =
                    r_local_gradients(i, component);
            }
        }
    }

    std::array<std::vector<double>, 3> scaled_powers;
    for (std::size_t component = 0; component < 3; ++component) {
        scaled_powers[component].resize(mBasisFunctionsOrder + 1);
        scaled_powers[component][0] = 1.0;
        for (IndexType order = 1;
             order <= mBasisFunctionsOrder;
             ++order) {
            scaled_powers[component][order] =
                scaled_powers[component][order - 1] *
                mDistanceVector[component] /
                static_cast<double>(order);
        }
    }

    for (IndexType order = 1;
         order <= mBasisFunctionsOrder;
         ++order) {
        const Matrix& r_derivatives = rDerivativeData[order - 1];
        if (mDim == 2) {
            for (IndexType y_order = 0;
                 y_order <= order;
                 ++y_order) {
                const IndexType x_order = order - y_order;
                const double value_coefficient =
                    scaled_powers[0][x_order] *
                    scaled_powers[1][y_order];
                const double gradient_x_coefficient = x_order > 0
                    ? scaled_powers[0][x_order - 1] *
                        scaled_powers[1][y_order]
                    : 0.0;
                const double gradient_y_coefficient = y_order > 0
                    ? scaled_powers[0][x_order] *
                        scaled_powers[1][y_order - 1]
                    : 0.0;
                for (std::size_t i = 0;
                     i < number_of_control_points;
                     ++i) {
                    const double derivative =
                        r_derivatives(i, y_order);
                    if (pShapeFunctions != nullptr) {
                        (*pShapeFunctions)[i] +=
                            derivative * value_coefficient;
                    }
                    if (pShapeFunctionGradients != nullptr && order > 1) {
                        (*pShapeFunctionGradients)(0, i) +=
                            derivative * gradient_x_coefficient;
                        (*pShapeFunctionGradients)(1, i) +=
                            derivative * gradient_y_coefficient;
                    }
                }
            }
        } else {
            std::size_t derivative_index = 0;
            for (int x_order = static_cast<int>(order);
                 x_order >= 0;
                 --x_order) {
                for (int y_order =
                         static_cast<int>(order) - x_order;
                     y_order >= 0;
                     --y_order) {
                    const IndexType z_order = order -
                        static_cast<IndexType>(x_order) -
                        static_cast<IndexType>(y_order);
                    const double value_coefficient =
                        scaled_powers[0][x_order] *
                        scaled_powers[1][y_order] *
                        scaled_powers[2][z_order];
                    const double gradient_x_coefficient = x_order > 0
                        ? scaled_powers[0][x_order - 1] *
                            scaled_powers[1][y_order] *
                            scaled_powers[2][z_order]
                        : 0.0;
                    const double gradient_y_coefficient = y_order > 0
                        ? scaled_powers[0][x_order] *
                            scaled_powers[1][y_order - 1] *
                            scaled_powers[2][z_order]
                        : 0.0;
                    const double gradient_z_coefficient = z_order > 0
                        ? scaled_powers[0][x_order] *
                            scaled_powers[1][y_order] *
                            scaled_powers[2][z_order - 1]
                        : 0.0;
                    for (std::size_t i = 0;
                         i < number_of_control_points;
                         ++i) {
                        const double derivative =
                            r_derivatives(i, derivative_index);
                        if (pShapeFunctions != nullptr) {
                            (*pShapeFunctions)[i] +=
                                derivative * value_coefficient;
                        }
                        if (pShapeFunctionGradients != nullptr && order > 1) {
                            (*pShapeFunctionGradients)(0, i) +=
                                derivative * gradient_x_coefficient;
                            (*pShapeFunctionGradients)(1, i) +=
                                derivative * gradient_y_coefficient;
                            (*pShapeFunctionGradients)(2, i) +=
                                derivative * gradient_z_coefficient;
                        }
                    }
                    ++derivative_index;
                }
            }
        }
    }
}

GapSbmSolidElementVolumetric::GapSbmSolidElementVolumetric(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : GapSbmSolidElement(NewId, pGeometry)
{
    CompactQuadratureGeometries();
}

GapSbmSolidElementVolumetric::GapSbmSolidElementVolumetric(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : GapSbmSolidElement(NewId, pGeometry, pProperties)
{
    CompactQuadratureGeometries();
}

Element::Pointer GapSbmSolidElementVolumetric::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << "GapSbmSolidElementVolumetric #" << Id()
        << " has no quadrature-point geometry parts.\n";

    std::vector<GeometryType::Pointer> quadrature_point_geometries;
    quadrature_point_geometries.reserve(NumberOfQuadraturePoints());
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        quadrature_point_geometries.push_back(
            GetGeometry().pGetGeometryPart(point_index));
    }

    auto p_coupling_geometry =
        Kratos::make_shared<CouplingGeometry<Node>>(
            std::move(quadrature_point_geometries));

    return Kratos::make_intrusive<GapSbmSolidElementVolumetric>(
        NewId,
        p_coupling_geometry,
        pProperties);
}

Element::Pointer GapSbmSolidElementVolumetric::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmSolidElementVolumetric>(
        NewId,
        pGeom,
        pProperties);
}

void GapSbmSolidElementVolumetric::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    InitializeMemberVariables();

    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << "GapSbmSolidElementVolumetric #" << Id()
        << " has no quadrature points.\n";
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CONSTITUTIVE_LAW))
        << "A constitutive law needs to be specified for element #"
        << Id() << ".\n";
    KRATOS_ERROR_IF(GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << "The constitutive law is null for element #" << Id() << ".\n";

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_properties = GetProperties();
    const std::size_t number_of_control_points =
        r_surrogate_geometry.PointsNumber();

    const std::size_t number_of_quadrature_points =
        NumberOfQuadraturePoints();
    mConstitutiveLawVector.resize(number_of_quadrature_points);
    KRATOS_ERROR_IF(
        mQuadraturePointReferenceWeights.size() !=
            number_of_quadrature_points ||
        mQuadraturePointCoordinates.size1() !=
            number_of_quadrature_points ||
        mQuadraturePointCoordinates.size2() != 3)
        << "Missing compact quadrature data in "
        << "GapSbmSolidElementVolumetric #" << Id() << ".\n";
    mQuadraturePointWeights.resize(number_of_quadrature_points, false);

    double total_integration_weight = 0.0;
    Vector shape_functions(number_of_control_points);
    for (std::size_t point_index = 0;
         point_index < number_of_quadrature_points;
         ++point_index) {
        double integration_weight =
            mQuadraturePointReferenceWeights[point_index];
        if (mDim == 2) {
            integration_weight *= GetProperties().Has(THICKNESS)
                ? GetProperties()[THICKNESS]
                : 1.0;
        }
        mQuadraturePointWeights[point_index] = integration_weight;
        total_integration_weight += integration_weight;
        SetQuadraturePointDistance(point_index);

        noalias(shape_functions) = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(shape_functions);

        mConstitutiveLawVector[point_index] =
            GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[point_index]->InitializeMaterial(
            r_properties,
            r_surrogate_geometry,
            shape_functions);
    }
    SetValue(INTEGRATION_WEIGHT, total_integration_weight);

    // Keep the inherited pointer valid for generic element infrastructure.
    mpConstitutiveLaw = mConstitutiveLawVector.front();

    KRATOS_CATCH("")
}

void GapSbmSolidElementVolumetric::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix,
        &rRightHandSideVector,
        rCurrentProcessInfo);
}

void GapSbmSolidElementVolumetric::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix,
        nullptr,
        rCurrentProcessInfo);
}

void GapSbmSolidElementVolumetric::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        nullptr,
        &rRightHandSideVector,
        rCurrentProcessInfo);
}

void GapSbmSolidElementVolumetric::CalculateAllContributions(
    MatrixType* pLeftHandSideMatrix,
    VectorType* pRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const std::size_t matrix_size =
        GetSurrogateGeometry().PointsNumber() * mDim;

    if (pLeftHandSideMatrix != nullptr) {
        if (pLeftHandSideMatrix->size1() != matrix_size ||
            pLeftHandSideMatrix->size2() != matrix_size) {
            pLeftHandSideMatrix->resize(
                matrix_size,
                matrix_size,
                false);
        }
        noalias(*pLeftHandSideMatrix) =
            ZeroMatrix(matrix_size, matrix_size);
    }

    if (pRightHandSideVector != nullptr) {
        if (pRightHandSideVector->size() != matrix_size) {
            pRightHandSideVector->resize(matrix_size, false);
        }
        noalias(*pRightHandSideVector) = ZeroVector(matrix_size);
    }

    KRATOS_ERROR_IF(
        mConstitutiveLawVector.size() != NumberOfQuadraturePoints())
        << "GapSbmSolidElementVolumetric #" << Id()
        << " is not initialized consistently: "
        << mConstitutiveLawVector.size() << " constitutive laws for "
        << NumberOfQuadraturePoints() << " quadrature points.\n";
    KRATOS_ERROR_IF(
        mQuadraturePointWeights.size() != NumberOfQuadraturePoints())
        << "GapSbmSolidElementVolumetric #" << Id()
        << " has an inconsistent quadrature cache.\n";

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points =
        r_surrogate_geometry.PointsNumber();
    KRATOS_ERROR_IF_NOT(mConstitutiveLawVector.front())
        << "Null first constitutive law in "
        << "GapSbmSolidElementVolumetric #" << Id() << ".\n";
    const std::size_t strain_size =
        mConstitutiveLawVector.front()->GetStrainSize();

    constexpr std::size_t chunk_size = 64;
    const std::size_t stack_rows = chunk_size * strain_size;

    Matrix B_stack(stack_rows, matrix_size);
    Matrix DB_stack;
    if (pLeftHandSideMatrix != nullptr) {
        DB_stack.resize(stack_rows, matrix_size, false);
    }

    Vector weighted_stress_stack;
    array_1d<double, 3> volume_force = ZeroVector(3);
    if (pRightHandSideVector != nullptr) {
        weighted_stress_stack.resize(stack_rows, false);
        volume_force = GetValue(BODY_FORCE);
    }

    Vector displacement_coefficient_vector(matrix_size);
    GetSolutionCoefficientVector(displacement_coefficient_vector);

    Vector shape_functions(number_of_control_points);
    Matrix shape_function_gradients(
        3,
        number_of_control_points);
    const TaylorDerivativeData taylor_derivative_data =
        CalculateTaylorDerivativeData();

    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters constitutive_parameters(
        r_surrogate_geometry,
        GetProperties(),
        rCurrentProcessInfo);

    Flags& r_options = constitutive_parameters.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    constitutive_parameters.SetStrainVector(
        constitutive_variables.StrainVector);
    constitutive_parameters.SetStressVector(
        constitutive_variables.StressVector);
    constitutive_parameters.SetConstitutiveMatrix(
        constitutive_variables.D);

    for (std::size_t chunk_begin = 0;
         chunk_begin < NumberOfQuadraturePoints();
         chunk_begin += chunk_size) {
        const std::size_t current_chunk_size = std::min(
            chunk_size,
            NumberOfQuadraturePoints() - chunk_begin);

        const std::size_t active_stack_rows =
            current_chunk_size * strain_size;
        for (std::size_t row = 0; row < active_stack_rows; ++row) {
            for (std::size_t column = 0;
                 column < matrix_size;
                 ++column) {
                B_stack(row, column) = 0.0;
            }
        }

        for (std::size_t stack_point_index = 0;
             stack_point_index < current_chunk_size;
             ++stack_point_index) {
            const std::size_t point_index =
                chunk_begin + stack_point_index;
            auto& p_constitutive_law =
                mConstitutiveLawVector[point_index];
            KRATOS_ERROR_IF_NOT(p_constitutive_law)
                << "Null constitutive law at quadrature point "
                << point_index << " of GapSbmSolidElementVolumetric #"
                << Id() << ".\n";
            KRATOS_ERROR_IF(
                p_constitutive_law->GetStrainSize() != strain_size)
                << "Inconsistent strain size at quadrature point "
                << point_index << " of GapSbmSolidElementVolumetric #"
                << Id() << ".\n";

            SetQuadraturePointDistance(point_index);
            EvaluateTaylorExpansion(
                taylor_derivative_data,
                pRightHandSideVector != nullptr
                    ? &shape_functions
                    : nullptr,
                &shape_function_gradients);
            FillBBlock(
                B_stack,
                stack_point_index,
                shape_function_gradients);

            const std::size_t row_offset =
                stack_point_index * strain_size;
            for (std::size_t strain_index = 0;
                 strain_index < strain_size;
                 ++strain_index) {
                double strain_value = 0.0;
                for (std::size_t column = 0;
                     column < matrix_size;
                     ++column) {
                    strain_value +=
                        B_stack(row_offset + strain_index, column) *
                        displacement_coefficient_vector[column];
                }
                constitutive_variables.StrainVector[strain_index] =
                    strain_value;
                constitutive_variables.StressVector[strain_index] = 0.0;
            }
            noalias(constitutive_variables.D) =
                ZeroMatrix(strain_size, strain_size);

            p_constitutive_law->CalculateMaterialResponse(
                constitutive_parameters,
                ConstitutiveLaw::StressMeasure_Cauchy);

            const double integration_weight =
                mQuadraturePointWeights[point_index];

            if (pLeftHandSideMatrix != nullptr) {
                const Matrix& r_constitutive_matrix =
                    constitutive_parameters.GetConstitutiveMatrix();
                for (std::size_t row = 0;
                     row < strain_size;
                     ++row) {
                    for (std::size_t column = 0;
                         column < matrix_size;
                         ++column) {
                        double value = 0.0;
                        for (std::size_t inner = 0;
                             inner < strain_size;
                             ++inner) {
                            value += r_constitutive_matrix(row, inner) *
                                B_stack(
                                    row_offset + inner,
                                    column);
                        }
                        DB_stack(row_offset + row, column) =
                            integration_weight * value;
                    }
                }
            }

            if (pRightHandSideVector != nullptr) {
                array_1d<double, 3> point_volume_force = volume_force;
                if (mHasQuadraturePointBodyForce[point_index]) {
                    for (std::size_t component = 0;
                         component < 3;
                         ++component) {
                        point_volume_force[component] =
                            mQuadraturePointBodyForces(
                                point_index,
                                component);
                    }
                }
                const Vector& r_stress_vector =
                    constitutive_parameters.GetStressVector();
                for (std::size_t strain_index = 0;
                     strain_index < strain_size;
                     ++strain_index) {
                    weighted_stress_stack[
                        row_offset + strain_index] =
                        integration_weight *
                        r_stress_vector[strain_index];
                }

                for (std::size_t i = 0;
                     i < number_of_control_points;
                     ++i) {
                    const std::size_t index = mDim * i;
                    for (std::size_t component = 0;
                         component < mDim;
                         ++component) {
                        (*pRightHandSideVector)[index + component] +=
                            integration_weight *
                            shape_functions[i] *
                            point_volume_force[component];
                    }
                }
            }
        }

        if (pLeftHandSideMatrix != nullptr) {
            for (std::size_t row = 0;
                 row < active_stack_rows;
                 ++row) {
                for (std::size_t i = 0; i < matrix_size; ++i) {
                    const double B_value = B_stack(row, i);
                    if (B_value == 0.0) {
                        continue;
                    }
#pragma omp simd
                    for (std::size_t j = 0; j < matrix_size; ++j) {
                        (*pLeftHandSideMatrix)(i, j) +=
                            B_value * DB_stack(row, j);
                    }
                }
            }
        }

        if (pRightHandSideVector != nullptr) {
            for (std::size_t row = 0;
                 row < active_stack_rows;
                 ++row) {
                const double weighted_stress =
                    weighted_stress_stack[row];
#pragma omp simd
                for (std::size_t i = 0; i < matrix_size; ++i) {
                    (*pRightHandSideVector)[i] -=
                        B_stack(row, i) * weighted_stress;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmSolidElementVolumetric::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    for (auto& p_constitutive_law : mConstitutiveLawVector) {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetSurrogateGeometry(),
            GetProperties(),
            rCurrentProcessInfo);
        p_constitutive_law->InitializeMaterialResponse(
            constitutive_law_parameters,
            ConstitutiveLaw::StressMeasure_Cauchy);
    }
}

void GapSbmSolidElementVolumetric::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points =
        r_surrogate_geometry.PointsNumber();
    const std::size_t matrix_size = number_of_control_points * mDim;

    Vector displacement_coefficient_vector(matrix_size);
    GetSolutionCoefficientVector(displacement_coefficient_vector);

    Vector weighted_stress;
    array_1d<double, 3> weighted_displacement = ZeroVector(3);
    double total_weight = 0.0;
    Vector shape_functions(number_of_control_points);
    Matrix shape_function_gradients(
        3,
        number_of_control_points);
    Matrix B;

    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        const double integration_weight =
            mQuadraturePointWeights[point_index];

        auto& p_constitutive_law =
            mConstitutiveLawVector[point_index];
        const std::size_t strain_size =
            p_constitutive_law->GetStrainSize();
        if (B.size1() != strain_size || B.size2() != matrix_size) {
            B.resize(strain_size, matrix_size, false);
        }
        noalias(B) = ZeroMatrix(strain_size, matrix_size);
        SetQuadraturePointDistance(point_index);
        noalias(shape_function_gradients) =
            ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(
            shape_function_gradients);
        FillBBlock(B, 0, shape_function_gradients);

        Vector strain = prod(B, displacement_coefficient_vector);
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters constitutive_parameters(
            r_surrogate_geometry,
            GetProperties(),
            rCurrentProcessInfo);
        Flags& r_options = constitutive_parameters.GetOptions();
        r_options.Set(
            ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,
            true);
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(
            ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
            true);
        constitutive_parameters.SetStrainVector(strain);
        constitutive_parameters.SetStressVector(
            constitutive_variables.StressVector);
        constitutive_parameters.SetConstitutiveMatrix(
            constitutive_variables.D);
        p_constitutive_law->CalculateMaterialResponse(
            constitutive_parameters,
            ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector& r_stress =
            constitutive_parameters.GetStressVector();
        if (weighted_stress.size() == 0) {
            weighted_stress = ZeroVector(r_stress.size());
        }
        noalias(weighted_stress) += integration_weight * r_stress;

        array_1d<double, 3> point_displacement = ZeroVector(3);
        noalias(shape_functions) =
            ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(shape_functions);
        for (std::size_t i = 0; i < number_of_control_points; ++i) {
            const std::size_t index = i * mDim;
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                point_displacement[component] +=
                    shape_functions[i] *
                    displacement_coefficient_vector[index + component];
            }
        }
        noalias(weighted_displacement) +=
            integration_weight * point_displacement;
        total_weight += integration_weight;

        p_constitutive_law->FinalizeMaterialResponse(
            constitutive_parameters,
            ConstitutiveLaw::StressMeasure_Cauchy);
    }

    if (total_weight > 0.0) {
        weighted_stress /= total_weight;
        weighted_displacement /= total_weight;
        SetValue(DISPLACEMENT, weighted_displacement);

        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, weighted_stress[0]);
            SetValue(CAUCHY_STRESS_YY, weighted_stress[1]);
            SetValue(CAUCHY_STRESS_XY, weighted_stress[2]);
        } else {
            SetValue(
                CAUCHY_STRESS_TENSOR,
                MathUtils<double>::StressVectorToTensor(
                    weighted_stress));
        }
    }
}

void GapSbmSolidElementVolumetric::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(NumberOfQuadraturePoints());
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        if (mConstitutiveLawVector[point_index]->Has(rVariable)) {
            mConstitutiveLawVector[point_index]->GetValue(
                rVariable,
                rOutput[point_index]);
        } else {
            rOutput[point_index] = 0.0;
        }
    }
}

void GapSbmSolidElementVolumetric::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(NumberOfQuadraturePoints());
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        rOutput[point_index] = ZeroVector(3);
        if (mConstitutiveLawVector[point_index]->Has(rVariable)) {
            mConstitutiveLawVector[point_index]->GetValue(
                rVariable,
                rOutput[point_index]);
        }
    }
}

int GapSbmSolidElementVolumetric::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    const int check_result =
        GapSbmSolidElement::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << "GapSbmSolidElementVolumetric #" << Id()
        << " has no quadrature-point geometry parts.\n";
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_GEOMETRIES))
        << "GapSbmSolidElementVolumetric #" << Id()
        << " has no NEIGHBOUR_GEOMETRIES.\n";
    KRATOS_ERROR_IF(GetValue(NEIGHBOUR_GEOMETRIES).size() != 1)
        << "GapSbmSolidElementVolumetric #" << Id()
        << " requires exactly one NEIGHBOUR_GEOMETRIES entry, got "
        << GetValue(NEIGHBOUR_GEOMETRIES).size() << ".\n";

    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        const double weight = GetQuadraturePointWeight(point_index);
        KRATOS_ERROR_IF(!std::isfinite(weight))
            << "Invalid weight " << weight << " at quadrature point "
            << point_index << " of GapSbmSolidElementVolumetric #"
            << Id() << ".\n";
    }

    return check_result;
}

std::size_t GapSbmSolidElementVolumetric::NumberOfQuadraturePoints() const
{
    return GetGeometry().NumberOfGeometryParts();
}

const GapSbmSolidElementVolumetric::GeometryType&
GapSbmSolidElementVolumetric::GetQuadraturePointGeometry(
    const std::size_t PointIndex) const
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints())
        << "Quadrature point index " << PointIndex
        << " is out of range [0, " << NumberOfQuadraturePoints()
        << ") in GapSbmSolidElementVolumetric #" << Id() << ".\n";
    return GetGeometry().GetGeometryPart(PointIndex);
}

double GapSbmSolidElementVolumetric::GetQuadraturePointWeight(
    const std::size_t PointIndex) const
{
    KRATOS_ERROR_IF(PointIndex >= mQuadraturePointReferenceWeights.size())
        << "Quadrature point index " << PointIndex
        << " is out of range in element #" << Id() << ".\n";
    return mQuadraturePointReferenceWeights[PointIndex];
}

void GapSbmSolidElementVolumetric::SetQuadraturePointDistance(
    const std::size_t PointIndex)
{
    if (mDistanceVector.size() != 3) {
        mDistanceVector.resize(3, false);
    }
    const auto surrogate_center =
        GetSurrogateGeometry().Center().Coordinates();
    for (std::size_t component = 0; component < 3; ++component) {
        mDistanceVector[component] =
            mQuadraturePointCoordinates(PointIndex, component) -
            surrogate_center[component];
    }
}

void GapSbmSolidElementVolumetric::CompactQuadratureGeometries()
{
    const std::size_t number_of_points =
        GetGeometry().NumberOfGeometryParts();
    if (number_of_points == 0 ||
        mQuadraturePointCoordinates.size1() == number_of_points) {
        return;
    }

    mQuadraturePointReferenceWeights.resize(number_of_points, false);
    mQuadraturePointCoordinates.resize(number_of_points, 3, false);
    mQuadraturePointBodyForces.resize(number_of_points, 3, false);
    mHasQuadraturePointBodyForce.assign(number_of_points, false);
    noalias(mQuadraturePointBodyForces) = ZeroMatrix(number_of_points, 3);

    auto p_representative_geometry =
        GetGeometry().pGetGeometryPart(0);
    for (std::size_t point_index = 0;
         point_index < number_of_points;
         ++point_index) {
        auto p_quadrature_geometry =
            GetGeometry().pGetGeometryPart(point_index);
        const auto& r_integration_points =
            p_quadrature_geometry->IntegrationPoints(
                p_quadrature_geometry->GetDefaultIntegrationMethod());
        KRATOS_ERROR_IF(r_integration_points.size() != 1)
            << "Each geometry part of GapSbmSolidElementVolumetric #"
            << Id() << " must contain exactly one integration point.\n";
        mQuadraturePointReferenceWeights[point_index] =
            r_integration_points.front().Weight();
        const auto center = p_quadrature_geometry->Center();
        for (std::size_t component = 0; component < 3; ++component) {
            mQuadraturePointCoordinates(point_index, component) =
                center[component];
        }
        if (p_quadrature_geometry->Has(BODY_FORCE)) {
            const auto& r_body_force =
                p_quadrature_geometry->GetValue(BODY_FORCE);
            mHasQuadraturePointBodyForce[point_index] = true;
            for (std::size_t component = 0; component < 3; ++component) {
                mQuadraturePointBodyForces(point_index, component) =
                    r_body_force[component];
            }
        }
        if (point_index > 0) {
            GetGeometry().SetGeometryPart(
                point_index,
                p_representative_geometry);
        }
    }
}

std::size_t
GapSbmSolidElementVolumetric::QuadraturePointsNumber() const
{
    return NumberOfQuadraturePoints();
}

array_1d<double, 3>
GapSbmSolidElementVolumetric::GetQuadraturePointCoordinates(
    const std::size_t PointIndex) const
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints())
        << "Quadrature point index out of range in element #"
        << Id() << ".\n";
    array_1d<double, 3> coordinates = ZeroVector(3);
    for (std::size_t component = 0; component < 3; ++component) {
        coordinates[component] =
            mQuadraturePointCoordinates(PointIndex, component);
    }
    return coordinates;
}

void GapSbmSolidElementVolumetric::
    SetQuadraturePointBodyForceComponent(
        const std::size_t PointIndex,
        const std::size_t ComponentIndex,
        const double Value)
{
    KRATOS_ERROR_IF(
        PointIndex >= NumberOfQuadraturePoints() ||
        ComponentIndex >= 3)
        << "Invalid compact body-force index in element #"
        << Id() << ".\n";
    mQuadraturePointBodyForces(PointIndex, ComponentIndex) = Value;
    mHasQuadraturePointBodyForce[PointIndex] = true;
}

void GapSbmSolidElementVolumetric::FillBBlock(
    Matrix& rBStack,
    const std::size_t StackPointIndex,
    const Matrix& rShapeFunctionGradients) const
{
    const std::size_t number_of_control_points =
        GetSurrogateGeometry().PointsNumber();
    const std::size_t matrix_size = number_of_control_points * mDim;
    const std::size_t strain_size = mDim == 2 ? 3 : 6;
    const std::size_t row_offset = StackPointIndex * strain_size;

    KRATOS_ERROR_IF(
        rBStack.size1() < row_offset + strain_size ||
        rBStack.size2() != matrix_size ||
        rShapeFunctionGradients.size1() != 3 ||
        rShapeFunctionGradients.size2() != number_of_control_points)
        << "Invalid B-stack dimensions in "
        << "GapSbmSolidElementVolumetric #" << Id() << ".\n";

    for (std::size_t i = 0;
         i < number_of_control_points;
         ++i) {
        const std::size_t matrix_index = i * mDim;
        const double derivative_x =
            rShapeFunctionGradients(0, i);
        const double derivative_y =
            rShapeFunctionGradients(1, i);

        if (mDim == 2) {
            rBStack(row_offset, matrix_index) = derivative_x;
            rBStack(row_offset + 1, matrix_index + 1) = derivative_y;
            rBStack(row_offset + 2, matrix_index) = derivative_y;
            rBStack(row_offset + 2, matrix_index + 1) = derivative_x;
        } else {
            const double derivative_z =
                rShapeFunctionGradients(2, i);
            rBStack(row_offset, matrix_index) = derivative_x;
            rBStack(row_offset + 1, matrix_index + 1) = derivative_y;
            rBStack(row_offset + 2, matrix_index + 2) = derivative_z;
            rBStack(row_offset + 3, matrix_index) = derivative_y;
            rBStack(row_offset + 3, matrix_index + 1) = derivative_x;
            rBStack(row_offset + 4, matrix_index + 1) = derivative_z;
            rBStack(row_offset + 4, matrix_index + 2) = derivative_y;
            rBStack(row_offset + 5, matrix_index) = derivative_z;
            rBStack(row_offset + 5, matrix_index + 2) = derivative_x;
        }
    }
}

} // Namespace Kratos
