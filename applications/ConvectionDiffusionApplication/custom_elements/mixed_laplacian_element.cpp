// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/mixed_laplacian_element.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
MixedLaplacianElement<TDim, TNumNodes>::MixedLaplacianElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
MixedLaplacianElement<TDim, TNumNodes>::MixedLaplacianElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
Element::Pointer MixedLaplacianElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianElement<TDim,TNumNodes>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
Element::Pointer MixedLaplacianElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianElement<TDim,TNumNodes>>(NewId, pGeom, pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
MixedLaplacianElement<TDim, TNumNodes>::~MixedLaplacianElement()
{
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_unknown_grad_var = r_settings.GetGradientVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const auto& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geometry = GetGeometry();

    // LHS initialization
    if(rLeftHandSideMatrix.size1() != ProblemSize || rLeftHandSideMatrix.size2() != ProblemSize) {
        rLeftHandSideMatrix.resize(ProblemSize, ProblemSize, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(ProblemSize, ProblemSize);

    // RHS initialization
    if(rRightHandSideVector.size() != ProblemSize) {
        rRightHandSideVector.resize(ProblemSize,false);
    }
    noalias(rRightHandSideVector) = ZeroVector(ProblemSize);

    // Get nodal data
    array_1d<double,3> temp;
    array_1d<double,TNumNodes> heat_flux_local;
    array_1d<double,TNumNodes> nodal_conductivity;
    BoundedMatrix<double,TNumNodes, TDim> grad_temp;
    for(std::size_t i_node = 0; i_node < TNumNodes; i_node++) {
        temp(i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_var);
        heat_flux_local(i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity(i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_diffusivity_var);
        const auto& r_grad_temp = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_grad_var);
        for (std::size_t d = 0; d < TDim; ++d) {
            grad_temp(i_node, d) = r_grad_temp(d);
        }
    }

    // Reading integration points and local gradients
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const auto& r_N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Allocate Gauss pt. geometry data
    double det_J0;
    array_1d<double,TNumNodes> aux_N;
    Element::GeometryType::JacobiansType J0;
    BoundedMatrix<double, TDim, TDim> inv_J0;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    r_geometry.Jacobian(J0, this->GetIntegrationMethod());

    // Loop Gauss pts.
    const std::size_t n_gauss = r_integration_points.size();
    for(std::size_t i_point = 0; i_point < n_gauss; ++i_point) {
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point], inv_J0, det_J0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(r_DN_De[i_point], inv_J0);

        // Get current Gauss pt. extra data
        noalias(aux_N) = row(r_N,i_point);
        const double w_g = r_integration_points[i_point].Weight() * det_J0;
        const double k_g = inner_prod(aux_N, nodal_conductivity);

    //     // Add mixed Laplacian RHS and LHS contribution
    //     noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //

    //     // Calculating the local RHS
    //     const double qgauss = inner_prod(N, heat_flux_local);

    //     noalias(rRightHandSideVector) += IntToReferenceWeight*qgauss*N;
    }


    // // RHS = ExtForces - K*temp;
    // for (unsigned int i = 0; i < number_of_points; i++)
    //     temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);

    // //axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
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
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::GetDofList(
    DofsVectorType& ElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
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

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
int MixedLaplacianElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedGradientVariable()) << "No Gradient Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_settings.GetGradientVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const auto& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geom = GetGeometry();
    for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); i_node++) {
        const auto& r_node = r_geom[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_gradient_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);

        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
        KRATOS_CHECK_DOF_IN_NODE(KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_X"), r_node);
        KRATOS_CHECK_DOF_IN_NODE(KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Y"), r_node);
        if (TDim == 3) {
            KRATOS_CHECK_DOF_IN_NODE(KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Z"), r_node);
        }
    }

    return Element::Check(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
Element::IntegrationMethod MixedLaplacianElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

template class MixedLaplacianElement<2,3>;
template class MixedLaplacianElement<3,4>;

} // Namespace Kratos
