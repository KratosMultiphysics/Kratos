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
#include "custom_elements/plain_stress_element.h"

#include "utilities/math_utils.h"
#include "utilities/function_parser_utility.h"


namespace Kratos
{

PlainStressElement::PlainStressElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

PlainStressElement::PlainStressElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer PlainStressElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PlainStressElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer PlainStressElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PlainStressElement>(NewId, pGeom, pProperties);
}

// Deconstructor

PlainStressElement::~PlainStressElement()
{
}


// From classical Laplacian
void PlainStressElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,2);
    Matrix InvJ0(2,2);
    Vector temp(mat_size);


    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    Vector GP_parameter_coord(2); 
    GP_parameter_coord = prod(r_geometry.Center(),J0[0]);


    Vector volume_force_local(integration_points.size()*2);
    /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);

    Matrix D = ZeroMatrix(3,3);
    D(0,0) = 1; 
    D(0,1) = nu;
    D(1,0) = nu;
    D(1,1) = 1;
    D(2,2) = (1-nu)/2;
    D *= E/(1-nu*nu);

    double factor = E/(1-nu*nu);
    for(std::size_t i_point = 0; i_point < integration_points.size()*2; i_point = i_point+2)
    {

        // volume_force_local[i_point] = factor*(cos(x)*sinh(y)-nu*cos(y)*sinh(x)) - E/2/(1+nu)*(cos(x)*sinh(y)+cos(y)*sinh(x));
        // volume_force_local[i_point+1] = factor*(sin(y)*cosh(x)-nu*sin(x)*cosh(y)) - E/2/(1+nu)*(-sin(x)*cosh(y) + sin(y)*cosh(x));

        volume_force_local[i_point] = this->GetValue(BODY_FORCE_X);
        volume_force_local[i_point+1] = this->GetValue(BODY_FORCE_Y);
    }

    for(std::size_t i_point = 0; i_point < integration_points.size(); ++i_point)
    {
        const IndexType IntegrationPointIndex = i_point ;
        // r_geometry.Jacobian(J0, IntegrationPointIndex, this->GetIntegrationMethod());
        
        double DetJ0;
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[i_point](0,0);
        Jacobian(0,1) = J0[i_point](0,1);
        Jacobian(1,0) = J0[i_point](1,0);
        Jacobian(1,1) = J0[i_point](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point); // these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0);

        // MODIFIED
        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);
        

        //------------------------------
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(trans(B), Matrix(prod(D, B))); //

        // // Calculating the local RHS
        for ( IndexType i = 0; i < number_of_points; ++i ) {
            const SizeType index = 2* i;

            for ( IndexType j = 0; j < 2; ++j )
                rRightHandSideVector[index + j] += IntToReferenceWeight * N[i] * volume_force_local[i_point+j];
        }
        
    }


    // RHS = ExtForces - K*temp;
    GetValuesVector(temp);

    // // RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")
}


// From classical Laplacian
void PlainStressElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
void PlainStressElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}





void PlainStressElement::EquationIdVector(
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

    void PlainStressElement::GetDofList(
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



int PlainStressElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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


Element::IntegrationMethod PlainStressElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void PlainStressElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // ConstitutiveLaw::Parameters constitutive_law_parameters(
    //     GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
    //     mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(
    //         constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);
    // }

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
    GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // Only one Integration Points 

    double x_coord_gauss_point = 0;
    double y_coord_gauss_point = 0;
    double rOutput = 0;

    for (IndexType i = 0; i < nb_nodes; ++i)
    {
        // KRATOS_WATCH(r_geometry[i])
        double output_solution_step_value = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_Y);
        rOutput += r_N(0, i) * output_solution_step_value;
        x_coord_gauss_point += r_N(0, i) * r_geometry[i].X0();
        y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y0();
    }        
    // exit(0);

    std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
    if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points[0].Weight() << std::endl;
        output_file.close();
    }
}

void PlainStressElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
}

//------------------------------------------------------------------------------------
// MODIFIED
//------------------------------------------------------------------------------------
array_1d<double, 3> PlainStressElement::GetBodyForce(
const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
const IndexType PointNumber
) const
{

    // for(unsigned int node_element = 0; node_element<number_of_nodes; node_element++)
    // {
    //     KRATOS_WATCH(r_geometry[node_element].FastGetSolutionStepValue(BODY_FORCE))
    // }
    //--------------------------------------------------------------------------------------

    // array_1d<double, 3> body_force;
    // for (IndexType i = 0; i < 3; ++i)
    //     body_force[i] = 0.0;

    // const auto& r_properties = GetProperties();
    // double density = 0.0;
    // if (r_properties.Has( DENSITY ))
    //     density = r_properties[DENSITY];

    // if (r_properties.Has( VOLUME_ACCELERATION ))
    //     noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];

    // const auto& r_geometry = GetGeometry();
    // if( r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
    //     Vector N(r_geometry.size());
    //     N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
    //     for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node)
    //         noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    // }

    // return body_force;


    // // FUTURE DEVELOPMENTS:
    // return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

void PlainStressElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    // const auto& r_geometry = GetGeometry();

    // const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // // Shape functions
    // rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    // KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // // Compute B
    // CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // // Compute equivalent F
    // GetValuesVector(rThisKinematicVariables.Displacements);
    // Vector strain_vector(mConstitutiveLawVector[0]->GetStrainSize());
    // noalias(strain_vector) = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    // ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    // rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

void PlainStressElement::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

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



void PlainStressElement::GetValuesVector(
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