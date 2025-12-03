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
#include <cmath>

// External includes

// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/cut_sbm_solid_element.h"

namespace Kratos
{

CutSbmSolidElement::CutSbmSolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

CutSbmSolidElement::CutSbmSolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer CutSbmSolidElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<CutSbmSolidElement>(NewId, GetSurrogateGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer CutSbmSolidElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<CutSbmSolidElement>(NewId, pGeom, pProperties);
}

// Deconstructor

CutSbmSolidElement::~CutSbmSolidElement()
{
}

void CutSbmSolidElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}


void CutSbmSolidElement::InitializeMaterial()
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

void CutSbmSolidElement::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2) << "CutSbmSolidElement momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    mBasisFunctionsOrder *= 2; 

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double integration_weight = r_integration_points[0].Weight()*thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void CutSbmSolidElement::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();

    const Point&  p_true = r_geometry.Center();            // true boundary
    const Point&  p_sur  = r_surrogate_geometry.Center();  // surrogate

    std::ofstream out("centers.txt", std::ios::app);       // append mode
    out << std::setprecision(15)                           // full precision
        << p_true.X() << ' ' << p_true.Y() << ' ' << p_true.Z() << ' '
        << p_sur .X() << ' ' << p_sur .Y() << ' ' << p_sur .Z() << '\n';
}

void CutSbmSolidElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetSurrogateGeometry().size() * 2;

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

void CutSbmSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_true_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    
    ConstitutiveLaw::Parameters values_true(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    values_true.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

    const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();

    noalias(rLeftHandSideMatrix) += integration_weight * prod(trans(B_sum), Matrix(prod(r_D_on_true, B_sum))); 

    KRATOS_CATCH("")
}

void CutSbmSolidElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
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

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    
    ConstitutiveLaw::Parameters values_true(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    values_true.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

    const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Vector& r_stress_vector_on_true = values_true.GetStressVector();

    //-----------------------------------------------------------------------------------
    Vector volume_force_local = this->GetValue(BODY_FORCE);
    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_control_points; ++i ) {
        const SizeType index = 2* i;

        for ( IndexType j = 0; j < 2; ++j )
            rRightHandSideVector[index + j] += integration_weight * N_sum_vec[i] * volume_force_local[j];
    }

    // RHS = ExtForces - K*temp;
    noalias(rRightHandSideVector) -= integration_weight * prod(trans(B_sum), r_stress_vector_on_true); 

    
    KRATOS_CATCH("")
}

void CutSbmSolidElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
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

void CutSbmSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    const auto& r_geometry = GetSurrogateGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(2 * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        const auto& r_node = r_geometry[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
    }
};




int CutSbmSolidElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

        // Check strain size
        KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
            << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
            << this->Id() << std::endl;
    }

    // Intentionally left blank: Cut-SBM element bypasses default geometry/size checks.
    // Returning 0 signals success.
    return 0;
    // return Element::Check(rCurrentProcessInfo);
}


Element::IntegrationMethod CutSbmSolidElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void CutSbmSolidElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{

    //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);

        Matrix B_sum = ZeroMatrix(mDim,mat_size);
        CalculateB(B_sum, grad_N_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        values_true.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

        mpConstitutiveLaw->FinalizeMaterialResponse(values_true, ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector sigma = values_true.GetStressVector();
        
        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        const double sigma_xx = sigma[0];
        const double sigma_yy = sigma[1];
        const double sigma_xy = sigma[2];
        const double von_mises_argument = sigma_xx * sigma_xx - sigma_xx * sigma_yy + sigma_yy * sigma_yy + 3.0 * sigma_xy * sigma_xy;
        const double von_mises_stress = von_mises_argument > 0.0 ? std::sqrt(von_mises_argument) : 0.0;
        SetValue(VON_MISES_STRESS, von_mises_stress);
        SetValue(VON_MISES_STRESS_IGA, von_mises_stress);

        // //---------------------


        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);

        Matrix N_sum_matrix = ZeroMatrix(mDim, mDim * number_of_control_points);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            for (IndexType idim = 0; idim < mDim; ++idim) {
                N_sum_matrix(idim, mDim * i + idim) = N_sum_vec(i);
            }
        }

        Vector displacement_true = prod(N_sum_matrix, old_displacement_coefficient_vector);

        SetValue(DISPLACEMENT, displacement_true);
}

void CutSbmSolidElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);
    constitutive_law_parameters.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

    // for (unsigned int i = 0; i < GetSurrogateGeometry().size(); i++) {

    //     std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
    //     outputFile << GetSurrogateGeometry()[i].GetId() << "  " <<GetSurrogateGeometry()[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
    //     outputFile.close();
    // }
}

void CutSbmSolidElement::CalculateOnIntegrationPoints(
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

    bool value_computed = false;

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
        value_computed = true;
    } else {
        double calculated_value = 0.0;
        value_computed = CalculateConstitutiveValue(rVariable, calculated_value, rCurrentProcessInfo);
        if (value_computed) {
            rOutput[0] = calculated_value;
        }
    }

    if (!value_computed) {
        KRATOS_WATCH(rVariable);
        KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED N THE IGA FRAMEWORK");
    }
}
void CutSbmSolidElement::CalculateOnIntegrationPoints(
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

void CutSbmSolidElement::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_surrogate_geometry.size();
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

double CutSbmSolidElement::GetCharacteristicGeometryLengthScalar() const
{
    KRATOS_ERROR_IF_NOT(this->Has(CHARACTERISTIC_GEOMETRY_LENGTH))
        << "CHARACTERISTIC_GEOMETRY_LENGTH not set for element " << this->Id() << std::endl;

    const array_1d<double,3>& characteristic_length_vector = this->GetValue(CHARACTERISTIC_GEOMETRY_LENGTH);
    const double characteristic_length = norm_2(characteristic_length_vector);

    KRATOS_ERROR_IF(characteristic_length <= 0.0)
        << "Non-positive characteristic length computed for element " << this->Id() << std::endl;

    return characteristic_length;
}



void CutSbmSolidElement::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = r_surrogate_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

bool CutSbmSolidElement::CalculateConstitutiveValue(
    const Variable<double>& rVariable,
    double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_true_geometry = GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const SizeType number_of_control_points = r_surrogate_geometry.size();
    const SizeType mat_size = number_of_control_points * mDim;

    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim, mat_size);
    CalculateB(B_sum, grad_N_sum);

    Vector displacement(mat_size);
    GetSolutionCoefficientVector(displacement);
    Vector strain = prod(B_sum, displacement);

    ConstitutiveLaw::Parameters values_true(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Flags& r_constitutive_law_options = values_true.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    values_true.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables constitutive_variables(strain_size);
    values_true.SetStrainVector(strain);
    values_true.SetStressVector(constitutive_variables.StressVector);
    values_true.SetConstitutiveMatrix(constitutive_variables.D);

    mpConstitutiveLaw->CalculateValue(values_true, rVariable, rValue);
    return true;
}

void CutSbmSolidElement::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
{
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    rValues.SetCharacteristicGeometryLength(GetCharacteristicGeometryLengthScalar());
    
    rValues.SetStrainVector(rStrain);
    rValues.SetStressVector(rConstitutiVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
}

void CutSbmSolidElement::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = GetSurrogateGeometry();
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
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
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

                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}


void CutSbmSolidElement::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = GetSurrogateGeometry();
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
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {

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
double CutSbmSolidElement::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double CutSbmSolidElement::ComputeTaylorTerm3D(
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
