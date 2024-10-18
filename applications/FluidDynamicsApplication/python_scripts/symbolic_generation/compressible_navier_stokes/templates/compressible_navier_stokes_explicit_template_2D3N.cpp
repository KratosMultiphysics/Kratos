//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Eduard Gómez
//

//automated_message

// System includes


// External includes


// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"


namespace Kratos {

/**
 * Returns the integration method for computation of midpoint magnitudes.
 * Computation of RHS integration method is chosen in the symbolic generator.
 */
template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicit<2,3>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template <>
void CompressibleNavierStokesExplicit<2,3>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != DofSize) {
        rResult.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<2,3>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (ElementalDofList.size() != DofSize) {
        ElementalDofList.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <>
array_1d<double,3> CompressibleNavierStokesExplicit<2,3>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GetIntegrationMethod());
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= static_cast<double>(NumNodes);
    midpoint_mom /= static_cast<double>(NumNodes);

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / std::pow(midpoint_rho, 2);
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / std::pow(midpoint_rho, 2);
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = 0.0;
    midpoint_rot_v[1] = 0.0;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<2, 3>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GetIntegrationMethod());
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= static_cast<double>(NumNodes);
    midpoint_mom /= static_cast<double>(NumNodes);

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v = ZeroMatrix(3, 3);
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, Dim*NumNodes> mom_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_mom_proj_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    mom_proj *= data.volume / static_cast<double>(NumNodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const IndexType aux = i_node * Dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < Dim; ++d) {
            AtomicAdd(r_mom_proj[d], mom_proj[aux + d]);
        }
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, 3> rho_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_rho_proj_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    rho_proj *= data.volume / static_cast<double>(NumNodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(DENSITY_PROJECTION), rho_proj[i_node]);
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, 3> tot_ener_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_tot_ener_proj_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    tot_ener_proj *= data.volume / static_cast<double>(NumNodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION), tot_ener_proj[i_node]);
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 12> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);
    const auto& DN_DX = data.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 12.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

    if (data.UseOSS)
    {
        //substitute_rhs_2D_OSS
    }
    else
    {
        //substitute_rhs_2D_ASGS
    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(NumNodes);

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2,3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Initialize and fill the mass matrix values
    constexpr double one_six = 1.0 / 6.0;
    constexpr double one_twelve = 1.0 / 12.0;
    constexpr unsigned int size = NumNodes * BlockSize;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six; rMassMatrix(0, 4) = one_twelve; rMassMatrix(0, 8) = one_twelve;
    rMassMatrix(1, 1) = one_six; rMassMatrix(1, 5) = one_twelve; rMassMatrix(1, 9) = one_twelve;
    rMassMatrix(2, 2) = one_six; rMassMatrix(2, 6) = one_twelve; rMassMatrix(2, 10) = one_twelve;
    rMassMatrix(3, 3) = one_six; rMassMatrix(3, 7) = one_twelve; rMassMatrix(3, 11) = one_twelve;
    rMassMatrix(4, 0) = one_twelve; rMassMatrix(4, 4) = one_six; rMassMatrix(4, 8) = one_twelve;
    rMassMatrix(5, 1) = one_twelve; rMassMatrix(5, 5) = one_six; rMassMatrix(5, 9) = one_twelve;
    rMassMatrix(6, 2) = one_twelve; rMassMatrix(6, 6) = one_six; rMassMatrix(6, 10) = one_twelve;
    rMassMatrix(7, 3) = one_twelve; rMassMatrix(7, 7) = one_six; rMassMatrix(7, 11) = one_twelve;
    rMassMatrix(8, 0) = one_twelve; rMassMatrix(8, 4) = one_twelve; rMassMatrix(8, 8) = one_six;
    rMassMatrix(9, 1) = one_twelve; rMassMatrix(9, 5) = one_twelve; rMassMatrix(9, 9) = one_six;
    rMassMatrix(10, 2) = one_twelve; rMassMatrix(10, 6) = one_twelve; rMassMatrix(10, 10) = one_six;
    rMassMatrix(11, 3) = one_twelve; rMassMatrix(11, 7) = one_twelve; rMassMatrix(11, 11) = one_six;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CompressibleNavierStokesExplicit<2,3>;

}
