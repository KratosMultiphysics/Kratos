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
#include "utilities/element_size_calculator.h"
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
    if(rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // RHS initialization
    if(rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize,false);
    }
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Get nodal data
    array_1d<double,TNumNodes> temp;
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
    array_1d<double,TDim> grad_k_g;
    array_1d<double,TDim> grad_temp_j;
    array_1d<double,TDim> aux_grad_N_i;
    array_1d<double,TDim> aux_grad_N_j;
    const std::size_t n_gauss = r_integration_points.size();
    for(std::size_t g = 0; g < n_gauss; ++g) {
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g], inv_J0, det_J0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(r_DN_De[g], inv_J0);

        // Get current Gauss pt. extra data
        noalias(aux_N) = row(r_N, g);
        const double f_g = inner_prod(aux_N, heat_flux_local);
        const double k_g = inner_prod(aux_N, nodal_conductivity);
        noalias(grad_k_g) = prod(trans(DN_DX), nodal_conductivity);
        const double w_g = r_integration_points[g].Weight() * det_J0;

        // Calculate stabilization constants
        const double h = ElementSizeCalculator<TDim, TNumNodes>::AverageElementSize(r_geometry);
        const double tau_temp = 0.1*std::pow(h,2)/k_g;
        const double tau_grad = 0.1;

        // Calculating the local RHS
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
            noalias(aux_grad_N_i) = row(DN_DX, i_node);

            // Volumetric force terms
            rRightHandSideVector(i_node*BlockSize) += w_g * f_g * aux_N(i_node);
            for (std::size_t d = 0; d < TDim; ++d) {
                rRightHandSideVector(i_node*BlockSize + d + 1) -= w_g * tau_temp * aux_grad_N_i(d) * f_g;
            }

            for (std::size_t j_node = 0; j_node < TNumNodes; ++j_node) {
                const double temp_j = temp(j_node);
                noalias(aux_grad_N_j) = row(DN_DX, j_node);
                noalias(grad_temp_j) = row(grad_temp, j_node);

                const double aux_1 = w_g * (1-tau_grad) * k_g * aux_N(j_node);
                const double aux_3 = w_g * (1-tau_grad) * aux_N(i_node) * aux_N(j_node);
                for (std::size_t d = 0; d < TDim; ++d) {

                    rRightHandSideVector(i_node*BlockSize) -= aux_1 * aux_grad_N_i(d) * grad_temp_j(d);
                    rLeftHandSideMatrix(i_node*BlockSize, j_node*BlockSize + d + 1) += aux_1 * aux_grad_N_i(d);

                    const double aux_2 = w_g * tau_grad * k_g * aux_grad_N_i(d) * aux_grad_N_j(d);
                    rRightHandSideVector(i_node*BlockSize) -= aux_2 * temp_j;
                    rLeftHandSideMatrix(i_node*BlockSize, j_node*BlockSize) += aux_2;

                    rRightHandSideVector(i_node*BlockSize + d + 1) -= aux_3 * grad_temp_j(d);
                    rLeftHandSideMatrix(i_node*BlockSize + d + 1, j_node*BlockSize + d + 1) += aux_3;

                    const double aux_4 = w_g * (1-tau_grad) * aux_N(i_node) * aux_grad_N_j(d);
                    rRightHandSideVector(i_node*BlockSize + d + 1) += aux_4 * temp_j;
                    rLeftHandSideMatrix(i_node*BlockSize + d + 1, j_node*BlockSize) -= aux_4;

                    for (std::size_t d2 = 0; d2 < TDim; ++d2) {
                        const double aux_5 = w_g * k_g * tau_temp * aux_grad_N_i(d) * aux_grad_N_j(d2);
                        rRightHandSideVector(i_node*BlockSize + d + 1) -= aux_5 * grad_temp_j(d2);
                        rLeftHandSideMatrix(i_node*BlockSize + d + 1, j_node*BlockSize + d2 + 1) += aux_5;

                        const double aux_6 = w_g * tau_temp *aux_grad_N_i(d) * grad_k_g(d2) * aux_N(j_node);
                        rRightHandSideVector(i_node*BlockSize + d + 1) -= aux_6 * grad_temp_j(d2);
                        rLeftHandSideMatrix(i_node*BlockSize + d + 1, j_node*BlockSize + d2 + 1) += aux_6;
                    }
                }
            }
        }

    }

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
    auto& r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_settings.GetGradientVariable();

    if(rResult.size() != LocalSize) {
        rResult.resize(LocalSize);
    }

    std::size_t local_index = 0;
    const auto& r_geom = GetGeometry();
    const auto& r_gradient_var_x = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_X");
    const auto& r_gradient_var_y = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Y");
    const auto& r_gradient_var_z = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Z");
    const std::size_t unknown_pos = r_geom[0].GetDofPosition(r_unknown_var);
    const std::size_t gradient_x_pos = r_geom[0].GetDofPosition(r_gradient_var_x);
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        rResult[local_index++] = r_geom[i_node].GetDof(r_unknown_var, unknown_pos).EquationId();
        rResult[local_index++] = r_geom[i_node].GetDof(r_gradient_var_x, gradient_x_pos).EquationId();
        rResult[local_index++] = r_geom[i_node].GetDof(r_gradient_var_y, gradient_x_pos + 1).EquationId();
        if (TDim == 3) {
            rResult[local_index++] = r_geom[i_node].GetDof(r_gradient_var_z, gradient_x_pos + 2).EquationId();
        }
    }
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void MixedLaplacianElement<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    auto& r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_settings.GetGradientVariable();

    if(rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    std::size_t local_index = 0;
    const auto& r_geom = GetGeometry();
    const auto& r_gradient_var_x = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_X");
    const auto& r_gradient_var_y = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Y");
    const auto& r_gradient_var_z = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Z");
    const std::size_t unknown_pos = r_geom[0].GetDofPosition(r_unknown_var);
    const std::size_t gradient_x_pos = r_geom[0].GetDofPosition(r_gradient_var_x);
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        rElementalDofList[local_index++] = r_geom[i_node].pGetDof(r_unknown_var, unknown_pos);
        rElementalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_x, gradient_x_pos);
        rElementalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_y, gradient_x_pos + 1);
        if (TDim == 3) {
            rElementalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_z, gradient_x_pos + 2);
        }
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
