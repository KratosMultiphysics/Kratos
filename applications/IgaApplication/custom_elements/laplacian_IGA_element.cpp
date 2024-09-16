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
void LaplacianIGAElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable(); // HeatFlux

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
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    // Initialize DN_DX
    Matrix DN_DX(number_of_points,dim);
    Vector temp(number_of_points);

    const double heat_flux = this->GetValue(r_volume_source_var);
    const double conductivity = this->GetProperties().GetValue(r_diffusivity_var);

    for(IndexType i_point = 0; i_point < r_integration_points.size(); ++i_point)
    {
        noalias(DN_DX) = r_DN_De[i_point] ; //prod(r_DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point);
        const double int_to_reference_weight = r_integration_points[i_point].Weight(); // * std::abs(DetJ0);

        noalias(rLeftHandSideMatrix) += int_to_reference_weight * conductivity * prod(DN_DX, trans(DN_DX));
        noalias(rRightHandSideVector) += int_to_reference_weight * heat_flux * N;
    }


    // RHS = ExtForces - K*temp;
    for (IndexType i = 0; i < number_of_points; i++) {
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
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
void LaplacianIGAElement::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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


int LaplacianIGAElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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


Element::IntegrationMethod LaplacianIGAElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

} // Namespace Kratos