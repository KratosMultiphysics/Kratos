//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//


// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/mat_variables.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/solid_3D_element.h"

#include "utilities/math_utils.h"
#include "utilities/function_parser_utility.h"


namespace Kratos
{

Solid3DElement::Solid3DElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

Solid3DElement::Solid3DElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer Solid3DElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<Solid3DElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer Solid3DElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<Solid3DElement>(NewId, pGeom, pProperties);
}

// Deconstructor

Solid3DElement::~Solid3DElement()
{
}

void Solid3DElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    
    // const auto& integration_points = this->IntegrationPoints(mThisIntegrationMethod);

    // //Constitutive Law initialisation
    // if ( mConstitutiveLawVector.size() != integration_points.size() )
    //         mConstitutiveLawVector.resize(integration_points.size());

    InitializeMaterial();
}


void Solid3DElement::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        KRATOS_WATCH(N_values)

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
    //KRATOS_WATCH(*mpConstitutiveLaw)

}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure Solid3DElement::GetStressMeasure() const
{
    return ConstitutiveLaw::StressMeasure_PK2;
}

// From classical Laplacian
void Solid3DElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3
    const SizeType mat_size = number_of_points * 3;

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // Compute the normals
    // array_1d<double, 3> normal_parameter_space;


    // r_geometry.Calculate(NORMAL, normal_parameter_space); // Gives the result in the parameter space !!
    // KRATOS_WATCH(normal_parameter_space)

    //-------------------------------------------------------------------------

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,3);
    Matrix InvJ0(3,3);
    

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    Vector GP_parameter_coord(3); 
    GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // check if we have tu put "GetInitialPosition"


    Vector volume_force_local(3);
    /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);

    Vector old_displacement(mat_size);
    GetValuesVector(old_displacement);

    //double factor = E/(1-nu*nu);

    volume_force_local[0] = this->GetValue(BODY_FORCE_X);
    volume_force_local[1] = this->GetValue(BODY_FORCE_Y);
    volume_force_local[2] = this->GetValue(BODY_FORCE_Z);

    //KRATOS_WATCH(volume_force_local)
    // r_geometry.Jacobian(J0, IntegrationPointIndex, this->GetIntegrationMethod());
    
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(dim,dim);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(0,2) = J0[0](0,2);    
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(1,2) = J0[0](1,2);
    Jacobian(2,0) = J0[0](2,0);
    Jacobian(2,1) = J0[0](2,1);
    Jacobian(2,2) = J0[0](2,2);    

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(DN_De[0],InvJ0);

    auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

    // MODIFIED
    Matrix B = ZeroMatrix(6,mat_size);

    CalculateB(B, DN_DX);


    //---------- MODIFIED ----------------------------------------------------------------
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    ConstitutiveVariables this_constitutive_variables(strain_size);

    Vector old_strain = prod(B,old_displacement);
    
    // Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStrainVector(old_strain);

    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2); 

    const Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();

    // KRATOS_WATCH(r_stress_vector)
    // KRATOS_WATCH(r_D)
    // KRATOS_WATCH(old_strain)
    //KRATOS_WATCH(B)
    //-----------------------------------------------------------------------------------
    

    //------------------------------
    noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(trans(B), Matrix(prod(r_D, B))); //

    // KRATOS_WATCH(rLeftHandSideMatrix)

    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_points; ++i ) {
        const SizeType index = 3* i;

        for ( IndexType j = 0; j < 3; ++j )
            rRightHandSideVector[index + j] += IntToReferenceWeight * N[i] * volume_force_local[+j];
    }


    // RHS = ExtForces - K*temp;
    

    // // RHS -= K*temp
    // TO DO 
    // Should be _int{B^T * \sigma}
    //noalias(rRightHandSideVector) -= prod(trans(B),r_stress_vector)*IntToReferenceWeight; 
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, old_displacement);
    //KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH("")
}


// From classical Laplacian
void Solid3DElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
void Solid3DElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);

}





void Solid3DElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();
        const GeometryType& r_geometry = GetGeometry();
        const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3

        if (rResult.size() != dim * number_of_control_points)
            rResult.resize(dim * number_of_control_points, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * dim;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();

        }

        KRATOS_CATCH("")
    };

    void Solid3DElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

        }

        KRATOS_CATCH("")
    };



int Solid3DElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // // Verify that the constitutive law exists
    // if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    // {
    //     KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    // }
    // else
    // {
    //     // Verify that the constitutive law has the correct dimension
    //     KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
    //         << "THICKNESS not provided for element " << this->Id() << std::endl;

    //     // Check strain size
    //     KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
    //         << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
    //         << this->Id() << std::endl;
    // }

    return Element::Check(rCurrentProcessInfo);
}


Element::IntegrationMethod Solid3DElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void Solid3DElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
// {
//     // ConstitutiveLaw::Parameters constitutive_law_parameters(
//     //     GetGeometry(), GetProperties(), rCurrentProcessInfo);

//     // mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

    
// /////////////////////////
//     const auto& r_geometry = GetGeometry();
//     const SizeType nb_nodes = r_geometry.size();
//     const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3

//     // Integration Points
//     const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
//     // Shape function values
//     const Matrix& r_N = r_geometry.ShapeFunctionsValues();

//     GeometryType::JacobiansType J0;
//     r_geometry.Jacobian(J0,this->GetIntegrationMethod());
//     // Get the parameter coordinates
//     Vector GP_parameter_coord(dim); 
//     GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // Only one Integration Points 

   
// }
{
    
    bool required = false;
        if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
            required = true;
        }
    if (required) {
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        const Properties& r_properties = GetProperties();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geom,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points

        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, 0, mThisIntegrationMethod);

        // Compute constitutive law variables
        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, integration_points);

        
        // Call the constitutive law to update material variables
        mpConstitutiveLaw->FinalizeMaterialResponse(Values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mpConstitutiveLaw->FinalizeSolutionStep( r_properties, r_geom, row( N_values, 0 ), rCurrentProcessInfo);
        
    }
}    

void Solid3DElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    
    bool required = false;
        if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
            required = true;
        }
    if (required) {
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        const Properties& r_properties = GetProperties();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geom,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points

        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, 0, mThisIntegrationMethod);

        // Compute constitutive law variables
        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, integration_points);

        
        // Call the constitutive law to update material variables
        mpConstitutiveLaw->InitializeMaterialResponse(Values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mpConstitutiveLaw->InitializeSolutionStep( r_properties, r_geom, row( N_values, 0 ), rCurrentProcessInfo);
        
    }
}    
    
    // ConstitutiveLaw::Parameters constitutive_law_parameters(
    //     GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
    // bool required = false;
    //     if (mConstitutiveLawVector[0]->RequiresInitializeMaterialResponse()) {
    //         required = true;
    //         break;
    //     }
    // }

bool Solid3DElement::UseElementProvidedStrain() const
{
    return false;
}

void Solid3DElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Setting the variables for the CL
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, IntegrationPoints);


    // Actually do the computations in the ConstitutiveLaw in local axes
    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

void Solid3DElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
    
}

// RESULTS ON GAUSS POINTS 11 06 24
void Solid3DElement::CalculateOnIntegrationPoints(
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

    if (rVariable == INTEGRATION_WEIGHT) {
        rOutput[0] = r_integration_points[0].Weight();
    } else if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
        //mpConstitutiveLaw->CalculateValue(rVariable, rOutput[0]);
    }
}

void Solid3DElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    { 
        const auto &r_geometry = GetGeometry();
        const auto &r_integration_points = r_geometry.IntegrationPoints();

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

// void Solid3DElement::CalculateOnIntegrationPoints(
//         const Variable<Vector>& rVariable,
//         std::vector<Vector>& rOutput,
//         const ProcessInfo& rCurrentProcessInfo
//     )
//     { 
//         const auto &r_geometry = GetGeometry();
//         const auto &r_integration_points = r_geometry.IntegrationPoints();

//         if (rOutput.size() != r_integration_points.size())
//         {
//             rOutput.resize(r_integration_points.size());
//         }

//     if (mpConstitutiveLaw->Has(rVariable)) {
//         mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
//     } else {
//         if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
//             rOutput[0].resize(6, false);
//             const auto& r_geometry = GetGeometry();
//             const unsigned int number_of_points = r_geometry.size();
//             const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3
//             const SizeType mat_size = number_of_points * 3;

//             const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//             const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

//             // Initialize DN_DX
//             Matrix DN_DX(number_of_points,3);
//             Matrix InvJ0(3,3);

//             // Initialize Jacobian
//             GeometryType::JacobiansType J0;
//             r_geometry.Jacobian(J0,this->GetIntegrationMethod());
//             Vector GP_parameter_coord(3); 
//             GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // check if we have tu put "GetInitialPosition"

//             double DetJ0;
//             Matrix Jacobian = ZeroMatrix(dim,dim);
//             Jacobian(0,0) = J0[0](0,0);
//             Jacobian(0,1) = J0[0](0,1);
//             Jacobian(0,2) = J0[0](0,2);    
//             Jacobian(1,0) = J0[0](1,0);
//             Jacobian(1,1) = J0[0](1,1);
//             Jacobian(1,2) = J0[0](1,2);
//             Jacobian(2,0) = J0[0](2,0);
//             Jacobian(2,1) = J0[0](2,1);
//             Jacobian(2,2) = J0[0](2,2);    

//             // Calculating inverse jacobian and jacobian determinant
//             MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

//             // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//             noalias(DN_DX) = prod(DN_De[0],InvJ0);


//             Vector old_displacement(mat_size);
//             GetValuesVector(old_displacement);
//             Matrix B = ZeroMatrix(6,mat_size);
//             CalculateB(B, DN_DX); 
//             noalias(rOutput[0]) = prod(B, old_displacement);
//             KRATOS_WATCH(rOutput[0])
//         }
//     }
// }

//------------------------------------------------------------------------------------
// MODIFIED
//------------------------------------------------------------------------------------
// array_1d<double, 3> Solid3DElement::GetBodyForce(
// const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
// const IndexType PointNumber
// ) const
// {


//     // // FUTURE DEVELOPMENTS:
//     // return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
// }

void Solid3DElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{

}

void Solid3DElement::CalculateB(
    Matrix &rB,
    Matrix &r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const GeometryType& r_geometry = GetGeometry();
        const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3
        const SizeType mat_size = number_of_control_points * dim;


        if (rB.size1() != 6 || rB.size2() != mat_size)
            rB.resize(6, mat_size);
        noalias(rB) = ZeroMatrix(6, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / dim;
            IndexType dirr = r % dim;

            // rB(0, r) = (dirr == 0) ? r_DN_DX(kr, 0) : 0.0;  // ε_xx
            // rB(1, r) = (dirr == 1) ? r_DN_DX(kr, 1) : 0.0;  // ε_yy
            // rB(2, r) = (dirr == 2) ? r_DN_DX(kr, 2) : 0.0;  // ε_zz
            // rB(3, r) = (dirr == 0) ? r_DN_DX(kr, 1) : (dirr == 1) ? r_DN_DX(kr, 0) : 0.0;  // ε_xy
            // rB(4, r) = (dirr == 1) ? r_DN_DX(kr, 2) : (dirr == 2) ? r_DN_DX(kr, 1) : 0.0;  // ε_yz
            // rB(5, r) = (dirr == 0) ? r_DN_DX(kr, 2) : (dirr == 2) ? r_DN_DX(kr, 0) : 0.0;  // ε_zx
            if (dirr == 0) {
                    rB(0, r)= r_DN_DX(kr, 0); // dN/dx
                    rB(3, r) = r_DN_DX(kr, 1); // dN/dy
                    rB(4, r) = r_DN_DX(kr, 2); // dN/dz
                }
                else if (dirr == 1) {
                    rB(1, r) = r_DN_DX(kr, 1); // dN/dy
                    rB(3, r) = r_DN_DX(kr, 0); // dN/dx
                    rB(5, r) = r_DN_DX(kr, 2); // dN/dz
                }
                else if (dirr == 2) {
                    rB(2, r) = r_DN_DX(kr, 2); // dN/dz
                    rB(4, r) = r_DN_DX(kr, 0); // dN/dx
                    rB(5, r) = r_DN_DX(kr, 1); // dN/dy
                }
        }
        //KRATOS_WATCH(rB)    

    }



void Solid3DElement::GetValuesVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const GeometryType& r_geometry = GetGeometry();
        const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 3
        const SizeType mat_size = number_of_control_points * dim;


        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];

        }
    }

} // Namespace Kratos