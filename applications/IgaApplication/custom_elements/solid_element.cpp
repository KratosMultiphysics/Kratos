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
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/solid_element.h"

#define ANALYTICAL_TANGENT

namespace Kratos
{

SolidElement::SolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

SolidElement::SolidElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer SolidElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer SolidElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SolidElement>(NewId, pGeom, pProperties);
}

// Deconstructor

SolidElement::~SolidElement()
{
}

void SolidElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();

    // calculate the integration weight
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
   
    Matrix InvJ0(2,2);
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    double DetJ0;
    Matrix Jacobian = ZeroMatrix(2,2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
    const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(DetJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
}


void SolidElement::InitializeMaterial()
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

void SolidElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int dim = r_DN_De[0].size2();
    const SizeType mat_size = number_of_control_points * dim;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    KRATOS_ERROR_IF(dim != 2) << "SolidElement momentarily only supports 2D elements, but the current element has dimension " << dim << std::endl;

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,dim);
    Matrix InvJ0(dim,dim);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    double DetJ0;
    Matrix Jacobian = ZeroMatrix(2,2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(r_DN_De[0],InvJ0);

    auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"

    Matrix B = ZeroMatrix(dim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement(mat_size);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);

    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    

    #ifdef ANALYTICAL_TANGENT
    Matrix r_D = ZeroMatrix(3,3);
    AnalyticalConstitutiveMatrix(r_D, old_strain);
    #else
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    #endif


    // KRATOS_WATCH(r_D_analytical)
    // KRATOS_WATCH(r_D)
    // KRATOS_WATCH("-------------")

    noalias(rLeftHandSideMatrix) += int_to_reference_weight * prod(trans(B), Matrix(prod(r_D, B))); 

    KRATOS_CATCH("")
}

void SolidElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

    // reading integration points and local gradients
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int dim = r_DN_De[0].size2();
    const SizeType mat_size = number_of_control_points * dim;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    KRATOS_ERROR_IF(dim != 2) << "SolidElement momentarily only supports 2D elements, but the current element has dimension " << dim << std::endl;
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(2,2);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    double DetJ0;
    Matrix Jacobian = ZeroMatrix(2,2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(r_DN_De[0],InvJ0);

    auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"

    Matrix B = ZeroMatrix(3,mat_size);
    CalculateB(B, DN_DX);

    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement(mat_size);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    //-----------------------------------------------------------------------------------
    Vector volume_force_local = this->GetValue(BODY_FORCE);
    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_control_points; ++i ) {
        const SizeType index = 2* i;

        for ( IndexType j = 0; j < 2; ++j )
            rRightHandSideVector[index + j] += int_to_reference_weight * N[i] * volume_force_local[j];
    }

    // RHS = ExtForces - K*temp;
    

    #ifdef ANALYTICAL_TANGENT
    Vector r_stress_vector = ZeroVector(3);
    AnalyticalStress(r_stress_vector, old_strain);
    #else
    const Vector& r_stress_vector = Values.GetStressVector();
    #endif

    noalias(rRightHandSideVector) -= int_to_reference_weight * prod(trans(B), r_stress_vector); 


    for (unsigned int i = 0; i < GetGeometry().size(); i++) {

        std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
        outputFile << GetGeometry()[i].GetId() << "  " <<GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
        outputFile.close();
    }

    KRATOS_CATCH("")
}

void SolidElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

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

    void SolidElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }

        KRATOS_CATCH("")
    };



int SolidElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

    return Element::Check(rCurrentProcessInfo);
}


Element::IntegrationMethod SolidElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void SolidElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    const SizeType mat_size = number_of_control_points * 2;

    // Initialize DN_DX
    const unsigned int dim = 2;
    Matrix DN_DX(number_of_control_points, 2);
    Matrix InvJ0(dim, dim);

    const GeometryType::ShapeFunctionsGradientsType& r_DN_De =
        r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    double DetJ0;

    Matrix Jacobian = ZeroMatrix(2, 2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian, InvJ0, DetJ0);

    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(r_DN_De[0], InvJ0);

    // calculate the B matrix
    Matrix B = ZeroMatrix(3, mat_size);
    CalculateB(B, DN_DX);

    // Prepare constitutive law parameters
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Flags& r_constitutive_law_options = Values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    ConstitutiveVariables this_constitutive_variables(strain_size);

    Vector old_displacement(mat_size);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(B, old_displacement);
    Values.SetStrainVector(old_strain);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double characteristic_geometry_length = 0.5*sqrt(mesh_size_uv[0]*mesh_size_uv[0] + mesh_size_uv[1]*mesh_size_uv[1]);
    Values.SetCharacteristicGeometryLength(characteristic_geometry_length);

    // compute and set stress solution at the integration point
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

    #ifdef ANALYTICAL_TANGENT
    Vector r_sigma = ZeroVector(3);
    AnalyticalStress(r_sigma, old_strain);
    #else
    const Vector& r_sigma = Values.GetStressVector();
    #endif

    SetValue(CAUCHY_STRESS_XX, r_sigma[0]);
    SetValue(CAUCHY_STRESS_YY, r_sigma[1]);
    SetValue(CAUCHY_STRESS_XY, r_sigma[2]);
    const double sigma_xx = r_sigma[0];
    const double sigma_yy = r_sigma[1];
    const double sigma_xy = r_sigma[2];
    const double von_mises_argument = sigma_xx * sigma_xx - sigma_xx * sigma_yy + sigma_yy * sigma_yy + 3.0 * sigma_xy * sigma_xy;
    const double von_mises_stress = von_mises_argument > 0.0 ? std::sqrt(von_mises_argument) : 0.0;
    SetValue(VON_MISES_STRESS, von_mises_stress);
    SetValue(VON_MISES_STRESS_IGA, von_mises_stress);

    Matrix strain_tensor = ZeroMatrix(2);
    strain_tensor(0, 0) = old_strain[0];
    strain_tensor(0, 1) = old_strain[2]; strain_tensor(1, 0) = old_strain[2];
    strain_tensor(1,1) = old_strain[1];

    SetValue(GREEN_LAGRANGE_STRAIN_TENSOR, strain_tensor);

    // update internal variables before requesting stresses
    mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
}

void SolidElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();
    const SizeType mat_size = number_of_control_points * 2;

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // Initialize DN_DX
    const unsigned int dimension = 2;
    Matrix DN_DX(number_of_control_points, 2);
    Matrix InvJ0(dimension, dimension);
    const GeometryType::ShapeFunctionsGradientsType& DN_De =
        r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    double detJ0;

    Matrix Jacobian = ZeroMatrix(2, 2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian, InvJ0, detJ0);

    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(DN_De[0], InvJ0);

    Matrix B = ZeroMatrix(3, mat_size);
    CalculateB(B, DN_DX);

    // Prepare constitutive law parameters
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Flags& r_constitutive_law_options = Values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    ConstitutiveVariables this_constitutive_variables(strain_size);

    Vector displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(displacement_coefficient_vector);
    Vector old_strain = prod(B, displacement_coefficient_vector);

    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double characteristic_geometry_length = 0.5 * std::sqrt(mesh_size_uv[0] * mesh_size_uv[0] + mesh_size_uv[1] * mesh_size_uv[1]);
    Values.SetCharacteristicGeometryLength(characteristic_geometry_length);

    Values.SetStrainVector(old_strain);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // update internal variables before requesting stresses
    mpConstitutiveLaw->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

    // KRATOS_WATCH("post")
}



void SolidElement::CalculateOnIntegrationPoints(
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
void SolidElement::CalculateOnIntegrationPoints(
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

void SolidElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize(r_integration_points.size());
    }

    if (mpConstitutiveLaw->Has(rVariable)) {
        if (rOutput[0].size() != mpConstitutiveLaw->GetStrainSize()) {
            rOutput[0].resize(mpConstitutiveLaw->GetStrainSize(), false);
        }
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
        KRATOS_WATCH(rVariable);
        KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}


void SolidElement::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
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

void SolidElement::CalculateDeltaPositionMatrix(
        const GeometryType& rGeometry,
        Matrix& rDeltaPosition) const
{
    const SizeType number_of_points = rGeometry.PointsNumber();
    if (rDeltaPosition.size1() != number_of_points || rDeltaPosition.size2() != 3) {
        rDeltaPosition.resize(number_of_points, 3, false);
    }

    for (IndexType i = 0; i < number_of_points; ++i) {
        const auto& r_node = rGeometry[i];
        rDeltaPosition(i, 0) = r_node.X() - r_node.X0();
        rDeltaPosition(i, 1) = r_node.Y() - r_node.Y0();
        rDeltaPosition(i, 2) = r_node.Z() - r_node.Z0();
    }
}



void SolidElement::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

bool SolidElement::CalculateConstitutiveValue(
    const Variable<double>& rVariable,
    double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    const GeometryType::ShapeFunctionsGradientsType& r_DN_De =
        r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int dim = r_DN_De[0].size2();
    const SizeType mat_size = number_of_control_points * dim;

    Matrix DN_DX(number_of_control_points, dim);
    Matrix InvJ0(dim, dim);

    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    double DetJ0;
    Matrix Jacobian = ZeroMatrix(dim, dim);
    for (unsigned int i = 0; i < dim; ++i) {
        for (unsigned int j = 0; j < dim; ++j) {
            Jacobian(i, j) = J0[0](i, j);
        }
    }

    MathUtils<double>::InvertMatrix(Jacobian, InvJ0, DetJ0);
    noalias(DN_DX) = prod(r_DN_De[0], InvJ0);

    Matrix B = ZeroMatrix(3, mat_size);
    CalculateB(B, DN_DX);

    Vector displacement(mat_size);
    GetSolutionCoefficientVector(displacement);
    Vector strain = prod(B, displacement);

    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    Flags& r_constitutive_law_options = Values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double characteristic_geometry_length = 0.5*sqrt(mesh_size_uv[0]*mesh_size_uv[0] + mesh_size_uv[1]*mesh_size_uv[1]);
    Values.SetCharacteristicGeometryLength(characteristic_geometry_length);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables constitutive_variables(strain_size);
    Values.SetStrainVector(strain);
    Values.SetStressVector(constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(constitutive_variables.D);

    mpConstitutiveLaw->CalculateValue(Values, rVariable, rValue);
    return true;
}
void SolidElement::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
{
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double characteristic_geometry_length = 0.5*sqrt(mesh_size_uv[0]*mesh_size_uv[0] + mesh_size_uv[1]*mesh_size_uv[1]);
    rValues.SetCharacteristicGeometryLength(characteristic_geometry_length);
    
    rValues.SetStrainVector(rStrain);
    rValues.SetStressVector(rConstitutiVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
}

} // Namespace Kratos
