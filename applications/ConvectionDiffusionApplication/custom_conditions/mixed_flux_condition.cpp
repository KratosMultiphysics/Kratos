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
#include "includes/convection_diffusion_settings.h"
#include "includes/variables.h"

// Application includes
#include "mixed_flux_condition.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
MixedFluxCondition<TDim,TNumNodes>::MixedFluxCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry)
    : Condition(NewId,pGeometry)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
MixedFluxCondition<TDim,TNumNodes>::MixedFluxCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId,pGeometry,pProperties)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
MixedFluxCondition<TDim,TNumNodes>::~MixedFluxCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
Condition::Pointer MixedFluxCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedFluxCondition<TDim,TNumNodes>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
Condition::Pointer MixedFluxCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedFluxCondition<TDim,TNumNodes>>(NewId, pGeom, pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Zet the zero LHS
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialize RHS vector
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Get convection diffusion settings container data
    ConvectionDiffusionSettings& rSettings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_gradient_var = rSettings.GetGradientVariable();
    const auto& r_surface_source_var = rSettings.GetSurfaceSourceVariable();
    //TODO: Doing this for each condition is super expensive, think about it
    const auto& r_gradient_var_x = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_X");
    const auto& r_gradient_var_y = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Y");
    const auto& r_gradient_var_z = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Z");

    // Get geometry data
    const auto& r_geometry = GetGeometry();
    std::size_t n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    Vector det_J_g_container(n_gauss);
    Matrix N_g_container(n_gauss, NumNodes);
    r_geometry.DeterminantOfJacobian(det_J_g_container, GetIntegrationMethod());
    noalias(N_g_container) = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Get nodal data and check fixity
    // Note that fixity will affect the boundary term to be imposed
    bool is_gradient_fixed = false;
    array_1d<double, NumNodes> external_flux_vect;
    for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        external_flux_vect(i_node) = r_node.FastGetSolutionStepValue(r_surface_source_var);
        if (r_node.IsFixed(r_gradient_var_x) || r_node.IsFixed(r_gradient_var_y) || r_node.IsFixed(r_gradient_var_z)) {
            is_gradient_fixed = true;
        }
    }
    KRATOS_ERROR_IF(is_gradient_fixed) << "Support for gradient Dirichlet BCs is not implemented yet." << std::endl;

    // Add external face heat flux contribution
    array_1d<double, NumNodes> aux_N_g;
    for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate Gauss pt. values
        noalias(aux_N_g) = row(N_g_container, i_gauss);
        const double f_g = inner_prod(aux_N_g, external_flux_vect);
        const double w_g = det_J_g_container(i_gauss) * r_integration_points[i_gauss].Weight();

        const double aux = w_g * f_g;
        for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
            rRightHandSideVector(i_node * BlockSize) += aux * aux_N_g(i_node);
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

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

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    auto& r_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_settings.GetGradientVariable();

    if(rConditionalDofList.size() != LocalSize) {
        rConditionalDofList.resize(LocalSize);
    }

    std::size_t local_index = 0;
    const auto& r_geom = GetGeometry();
    const auto& r_gradient_var_x = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_X");
    const auto& r_gradient_var_y = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Y");
    const auto& r_gradient_var_z = KratosComponents<VariableData>::Get(r_gradient_var.Name()+"_Z");
    const std::size_t unknown_pos = r_geom[0].GetDofPosition(r_unknown_var);
    const std::size_t gradient_x_pos = r_geom[0].GetDofPosition(r_gradient_var_x);
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        rConditionalDofList[local_index++] = r_geom[i_node].pGetDof(r_unknown_var, unknown_pos);
        rConditionalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_x, gradient_x_pos);
        rConditionalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_y, gradient_x_pos + 1);
        if (TDim == 3) {
            rConditionalDofList[local_index++] = r_geom[i_node].pGetDof(r_gradient_var_z, gradient_x_pos + 2);
        }
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
GeometryData::IntegrationMethod MixedFluxCondition<TDim,TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

// Input and Output ///////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
std::string MixedFluxCondition<TDim,TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "MixedFluxCondition #" << Id();
    return buffer.str();
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "MixedFluxCondition #" << Id();
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::PrintData(std::ostream& rOStream) const
{
    rOStream << "MixedFluxCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
MixedFluxCondition<TDim,TNumNodes>::MixedFluxCondition():
    Condition()
{
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template<std::size_t TDim, std::size_t TNumNodes>
void MixedFluxCondition<TDim,TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

template class MixedFluxCondition<2,2>;
template class MixedFluxCondition<3,3>;

}
