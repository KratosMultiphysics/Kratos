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
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/solid_element.h"

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
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

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
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

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

    const Matrix& r_D = Values.GetConstitutiveMatrix();

    noalias(rLeftHandSideMatrix) += int_to_reference_weight * prod(trans(B), Matrix(prod(r_D, B))); 

    for (unsigned int i = 0; i < number_of_control_points; i++) {

        // std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
        // outputFile << r_geometry[i].GetId() << "  " <<r_geometry[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
        // outputFile.close();
    }

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
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

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

    const Vector& r_stress_vector = Values.GetStressVector();
    //-----------------------------------------------------------------------------------
    Vector volume_force_local = this->GetValue(BODY_FORCE);
    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_control_points; ++i ) {
        const SizeType index = 2* i;

        for ( IndexType j = 0; j < 2; ++j )
            rRightHandSideVector[index + j] += int_to_reference_weight * N[i] * volume_force_local[j];
    }

    // RHS = ExtForces - K*temp;
    noalias(rRightHandSideVector) -= int_to_reference_weight * prod(trans(B), r_stress_vector); 

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
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

    // compute and set stress solution at the integration point
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_control_points = r_geometry.size();

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());     

    const SizeType mat_size = number_of_control_points * 2;

    // Initialize DN_DX
    const unsigned int dim = 2;
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(dim,dim);

    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    double DetJ0;

    Matrix Jacobian = ZeroMatrix(2,2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(r_DN_De[0],InvJ0);

    // calculate the B matrix
    Matrix B = ZeroMatrix(3,mat_size);
    CalculateB(B, DN_DX);

    // GET STRESS VECTOR
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    ConstitutiveVariables this_constitutive_variables(strain_size);

    Vector old_displacement(mat_size);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    Values.SetStrainVector(old_strain);

    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

    const Vector sigma = Values.GetStressVector();

    SetValue(CAUCHY_STRESS_XX, sigma[0]);
    SetValue(CAUCHY_STRESS_YY, sigma[1]);
    SetValue(CAUCHY_STRESS_XY, sigma[2]);
    // //---------------------
}

void SolidElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
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

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
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

void SolidElement::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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

} // Namespace Kratos