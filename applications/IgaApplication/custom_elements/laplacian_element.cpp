//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
// System includes
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_elements/laplacian_element.h"


namespace Kratos
{

LaplacianElement::LaplacianElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

LaplacianElement::LaplacianElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer LaplacianElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, pGeom, pProperties);
}

// Deconstructor

LaplacianElement::~LaplacianElement()
{
}


void LaplacianElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    SetValue(INTEGRATION_WEIGHT, r_integration_points[0].Weight());
}

// From classical Laplacian
void LaplacianElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}


// From classical Laplacian
void LaplacianElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity

    // reading integration points and local gradients
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_DN_De[0].size2();
    
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,dim);
    Vector temp(number_of_points);

    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    for(IndexType i_point = 0; i_point < r_integration_points.size(); ++i_point)
    {
        noalias(DN_DX) = r_DN_De[i_point] ; //prod(r_DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point);
        const double int_to_reference_weight = r_integration_points[i_point].Weight(); // * std::abs(DetJ0);

        noalias(rLeftHandSideMatrix) += int_to_reference_weight * conductivity * prod(DN_DX, trans(DN_DX));
    }

    KRATOS_CATCH("")

}

// From modification of lefthandside
void LaplacianIGAElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    KRATOS_TRY
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity

    // reading integration points and local gradients
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_DN_De[0].size2();
    
    //resizing as needed the LHS
    if(rMassMatrix.size1() != number_of_points)
        rMassMatrix.resize(number_of_points,number_of_points,false);
    noalias(rMassMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting M

    // const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);
    const double conductivity = 1.0;

    for(IndexType i_point = 0; i_point < r_integration_points.size(); ++i_point)
    {
        auto N = row(N_gausspoint,i_point);
        const double int_to_reference_weight = r_integration_points[i_point].Weight(); // * std::abs(DetJ0);

        noalias(rMassMatrix) += int_to_reference_weight * conductivity * outer_prod(N, N);
    }

    KRATOS_CATCH("")

}

// From classical Laplacian
void LaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable(); // HeatFlux
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);
    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature

    // reading integration points and local gradients
    const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_DN_De[0].size2();
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,dim);

    const double heat_flux = this->GetValue(r_volume_source_var);

    Vector temperature_old_iteration(number_of_points);
    for (IndexType i = 0; i < number_of_points; i++) {
        temperature_old_iteration[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
    }

    for(IndexType i_point = 0; i_point < r_integration_points.size(); ++i_point)
    {
        noalias(DN_DX) = r_DN_De[i_point] ; //prod(r_DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point);
        const double int_to_reference_weight = r_integration_points[i_point].Weight(); // * std::abs(DetJ0);

        noalias(rRightHandSideVector) -= int_to_reference_weight * conductivity * prod(prod(DN_DX, trans(DN_DX)),temperature_old_iteration);
        noalias(rRightHandSideVector) += int_to_reference_weight * heat_flux * N;
    }

    // for (unsigned int i = 0; i < number_of_points; i++) {
    //     std::ofstream outputFile("txt_files/Id_active_control_points.txt", std::ios::app);
    //     outputFile << r_geometry[i].GetId() << "  " <<r_geometry[i].GetDof(r_unknown_var).EquationId() <<"\n";
    //     outputFile.close();
    // }

    KRATOS_CATCH("")
}

void LaplacianElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_settings.GetGradientVariable();

    auto& r_geometry = const_cast<GeometryType&>(GetGeometry());
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_DN_De[0].size2();

    array_1d<double,3> grad_unknown = ZeroVector(3);
    for (IndexType i = 0; i < number_of_points; ++i) {
        const double value = r_geometry[i].FastGetSolutionStepValue(r_unknown_var);
        for (unsigned int d = 0; d < dim; ++d) {
            grad_unknown[d] += r_DN_De[0](i, d) * value;
        }
    }

    r_geometry.SetValue(r_gradient_var, grad_unknown);

    KRATOS_CATCH("")
}

void LaplacianElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes) {
        rResult.resize(number_of_nodes);
    }
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes) {
        ElementalDofList.resize(number_of_nodes);
    }
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }
}


int LaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const auto& r_volume_source_var = r_settings.GetVolumeSourceVariable();

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


Element::IntegrationMethod LaplacianElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void LaplacianIGAElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight());
}
} // Namespace Kratos
