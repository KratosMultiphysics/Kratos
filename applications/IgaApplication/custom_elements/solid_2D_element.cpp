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
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/solid_2D_element.h"

#include "utilities/math_utils.h"
#include "utilities/function_parser_utility.h"


namespace Kratos
{

Solid2DElement::Solid2DElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

Solid2DElement::Solid2DElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer Solid2DElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<Solid2DElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer Solid2DElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<Solid2DElement>(NewId, pGeom, pProperties);
}

// Deconstructor

Solid2DElement::~Solid2DElement()
{
}

void Solid2DElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}


void Solid2DElement::InitializeMaterial()
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


// From classical Laplacian
void Solid2DElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension(); // dim = 2
    const SizeType mat_size = number_of_points * 2;
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS


    //-------------------------------------------------------------------------

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,2);
    Matrix InvJ0(2,2);
    

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    Vector GP_parameter_coord(2); 
    GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // check if we have tu put "GetInitialPosition"


    Vector volume_force_local(2);
    /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Vector old_displacement(mat_size);
    GetValuesVector(old_displacement);

    volume_force_local[0] = this->GetValue(BODY_FORCE_X);// 
    volume_force_local[1] = this->GetValue(BODY_FORCE_Y);//

    // if (norm_2(volume_force_local) > 1e-10)  KRATOS_WATCH(volume_force_local)

    // r_geometry.Jacobian(J0, IntegrationPointIndex, this->GetIntegrationMethod());
    
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(2,2);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = prod(DN_De[0],InvJ0);

    auto N = row(N_gausspoint,0); // these are the N which correspond to the gauss point "i_point"

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);

    // MODIFIED
    Matrix B = ZeroMatrix(3,mat_size);

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

    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy); 

    const Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    //-----------------------------------------------------------------------------------
    

    //------------------------------
    noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(trans(B), Matrix(prod(r_D, B))); //

    // // Calculating the local RHS
    for ( IndexType i = 0; i < number_of_points; ++i ) {
        const SizeType index = 2* i;

        for ( IndexType j = 0; j < 2; ++j )
            rRightHandSideVector[index + j] += IntToReferenceWeight * N[i] * volume_force_local[j];
    }


    // RHS = ExtForces - K*temp;
    

    // // RHS -= K*temp
    // TO DO 
    // Should be _int{B^T * \sigma}
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,old_displacement); 

    KRATOS_CATCH("")
}


// From classical Laplacian
void Solid2DElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
void Solid2DElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}





void Solid2DElement::EquationIdVector(
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

    void Solid2DElement::GetDofList(
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



int Solid2DElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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


Element::IntegrationMethod Solid2DElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void Solid2DElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

    // retrieve integration weight
    double integration_weight = GetValue(INTEGRATION_WEIGHT);

    
/////////////////////////
    const auto& r_geometry = GetGeometry();
    const SizeType nb_nodes = r_geometry.size();

    // Integration Points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    // Shape function values
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    // Get the parameter coordinates
    Vector GP_parameter_coord(2); 
    GP_parameter_coord = r_geometry.Center(); // Only one Integration Points 

    double x_coord_gauss_point = 0;
    double y_coord_gauss_point = 0;
    double rOutput_x = 0;
    double rOutput_y = 0;

    for (IndexType i = 0; i < nb_nodes; ++i)
    {
        // KRATOS_WATCH(r_geometry[i])
        double output_solution_step_value = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
        rOutput_x += r_N(0, i) * output_solution_step_value;

        double output_solution_step_value_y = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_Y);
        rOutput_y += r_N(0, i) * output_solution_step_value_y;

        x_coord_gauss_point += r_N(0, i) * r_geometry[i].X0();
        y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y0();
    }        
    // exit(0);

    std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
    if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput_x << " " << GP_parameter_coord[0] << " " << GP_parameter_coord[1] << " " << integration_weight << std::endl;
        output_file.close();
    } 



    if (x_coord_gauss_point <= 2.0) {
        std::ofstream output_file("txt_files/output_results_GPs_master.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << rOutput_x << " " << rOutput_y << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " << integration_weight<< std::endl;
            output_file.close();
        } 
    } else {
        std::ofstream output_file("txt_files/output_results_GPs_slave.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << rOutput_x << " " << rOutput_y << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " << integration_weight << std::endl;
            output_file.close();
        } 
    }
    


    std::ofstream outputFile("txt_files/Gauss_Point_coordinates.txt", std::ios::app);
    if (!outputFile.is_open())
    {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return;
    }
    outputFile << std::setprecision(14); // Set precision to 10^-14
    outputFile << GP_parameter_coord[0] << "  " << GP_parameter_coord[1] <<"\n";
    outputFile.close();
}

void Solid2DElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
}


// RESULTS ON GAUSS POINTS 11 06 24
void Solid2DElement::CalculateOnIntegrationPoints(
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
void Solid2DElement::CalculateOnIntegrationPoints(
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



//------------------------------------------------------------------------------------
// MODIFIED
//------------------------------------------------------------------------------------
array_1d<double, 3> Solid2DElement::GetBodyForce(
const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
const IndexType PointNumber
) const
{


    // // FUTURE DEVELOPMENTS:
    // return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

void Solid2DElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{

}

void Solid2DElement::CalculateB(
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



void Solid2DElement::GetValuesVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 2 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

} // Namespace Kratos