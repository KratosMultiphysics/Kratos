// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/stokes_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

StokesElement::StokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

StokesElement::StokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer StokesElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer StokesElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElement>(NewId, pGeom, pProperties);
}

// Deconstructor

StokesElement::~StokesElement()
{
}


// From classical Laplacian
void StokesElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable(); // HeatFlux

    // reading integration points and local gradients
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = DN_De[0].size2();
    
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,dim);
    Vector temp(number_of_points);

    // Initialize Jacobian
    // Matrix InvJ0(dim,dim);
    // GeometryType::JacobiansType J0;
    // r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    const double heat_flux = this->GetValue(r_volume_source_var);
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    // double DetJ0;
    // Matrix Jacobian = ZeroMatrix(dim,dim);
    // Jacobian(0,0) = J0[0](0,0);
    // Jacobian(0,1) = J0[0](0,1);
    // Jacobian(1,0) = J0[0](1,0);
    // Jacobian(1,1) = J0[0](1,1);

    // if (dim > 2) {
    //     Jacobian(0,2) = J0[0](0,2);
    //     Jacobian(1,2) = J0[0](1,2);
    //     Jacobian(2,0) = J0[0](2,0);
    //     Jacobian(2,1) = J0[0](2,1);
    //     Jacobian(2,2) = J0[0](2,2);
    // }
    // // Calculating inverse jacobian and jacobian determinant
    // MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    for(std::size_t i_point = 0; i_point < integration_points.size(); ++i_point)
    {
        noalias(DN_DX) = DN_De[i_point] ; //prod(DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point);
        const double IntToReferenceWeight = integration_points[i_point].Weight(); // * std::abs(DetJ0);

        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity * prod(DN_DX, trans(DN_DX));
        noalias(rRightHandSideVector) += IntToReferenceWeight * heat_flux * N;
    }


    // RHS = ExtForces - K*temp;
    for (unsigned int i = 0; i < number_of_points; i++) {
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);

        // std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
        // outputFile << r_geometry[i].GetId() << "  " <<r_geometry[i].GetDof(r_unknown_var).EquationId() <<"\n";
        // outputFile.close();
    }

    // RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")
}


// From classical Laplacian
void StokesElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
void StokesElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


//************************************************************************************
//************************************************************************************

void StokesElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const unsigned int dim = 2;
        const SizeType number_of_control_points = GetGeometry().size();
        const unsigned int dof_size  = number_of_control_points * (dim+1);
        
        if (rResult.size() != dof_size)
            rResult.resize(dof_size, false);

        const IndexType pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * (dim+1);
            rResult[index]     = GetGeometry()[i].GetDof(VELOCITY_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(VELOCITY_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(PRESSURE, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    void StokesElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(4 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
        }

        KRATOS_CATCH("")
    };



int StokesElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geom = GetGeometry();

    for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);

        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}


Element::IntegrationMethod StokesElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void StokesElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    const auto& r_geometry = GetGeometry();
    const SizeType nb_nodes = r_geometry.size();

    // Integration Points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    // Shape function values
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int dim = DN_De[0].size2();

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;
    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    double rOutput = 0;
    Vector rOutput_gradient = ZeroVector(3);
    for (IndexType i = 0; i < nb_nodes; ++i)
    {
        double output_solution_step_value = r_geometry[i].GetSolutionStepValue(r_unknown_var);
        rOutput += r_N(0, i) * output_solution_step_value;
        for (IndexType idim = 0; idim < dim; idim++) {
            rOutput_gradient(idim) += DN_De[0](i,idim) * output_solution_step_value;
        }
    } 


    #pragma omp critical
    {
        std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() << " " << integration_points[0].Weight() 
                    << " " << rOutput_gradient(0)<< " " << rOutput_gradient(1) << " " << rOutput_gradient(2) << std::endl;
        output_file.close();
        }   
    }  
}

void StokesElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
}



} // Namespace Kratos