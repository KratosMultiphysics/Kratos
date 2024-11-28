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
#include "custom_elements/laplacian_IGA_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

LaplacianIGAElement::LaplacianIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

LaplacianIGAElement::LaplacianIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer LaplacianIGAElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianIGAElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, pGeom, pProperties);
}

// Deconstructor

LaplacianIGAElement::~LaplacianIGAElement()
{
}


// From classical Laplacian
void LaplacianIGAElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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

    const double heat_flux = this->GetValue(r_volume_source_var);
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    // Time-related variables
    const double delta_t = r_process_info[DELTA_TIME];
    const double theta = 1.0; 
    const double heat_flux_old = this->GetValue(FACE_HEAT_FLUX);
    

    noalias(DN_DX) = DN_De[0] ; //prod(DN_De[i_point],InvJ0);
    auto N = row(N_gausspoint,0);

    double previous_temperature = 0.0;
    Vector grad_T_previous = ZeroVector(dim);
    for (unsigned int i = 0; i < number_of_points; ++i)
    {
        previous_temperature += r_geometry[i].FastGetSolutionStepValue(r_unknown_var, 1) * N[i]; // T^(n)
        for (unsigned int d = 0; d < dim; ++d)
        {
            grad_T_previous[d] += DN_DX(i, d) * r_geometry[i].FastGetSolutionStepValue(r_unknown_var, 1); // ∇T^(n)
        }
    }
    
    const double IntToReferenceWeight = integration_points[0].Weight(); // * std::abs(DetJ0);

    noalias(rLeftHandSideMatrix) += theta * IntToReferenceWeight * conductivity * prod(DN_DX, trans(DN_DX));
    noalias(rRightHandSideVector)-= (1.0-theta) * IntToReferenceWeight * conductivity * prod(DN_DX, grad_T_previous);

    noalias(rRightHandSideVector) += theta * IntToReferenceWeight * heat_flux * N;
    noalias(rRightHandSideVector) += (1.0-theta) * IntToReferenceWeight * heat_flux_old * N;

    // Add mass term to LHS
    noalias(rLeftHandSideMatrix) += IntToReferenceWeight * (1.0 / delta_t) * outer_prod(N, N);
    // Add mass term contribution to RHS
    noalias(rRightHandSideVector) += IntToReferenceWeight * N * previous_temperature * 1.0 / delta_t;




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
void LaplacianIGAElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
void LaplacianIGAElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}





void LaplacianIGAElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianIGAElement::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }
}


int LaplacianIGAElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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


Element::IntegrationMethod LaplacianIGAElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void LaplacianIGAElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight() );

    // Set value for the heat_flux_old for next time-step
    const double heat_flux = this->GetValue(HEAT_FLUX);
    SetValue(FACE_HEAT_FLUX, heat_flux);
}

void LaplacianIGAElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
}



} // Namespace Kratos