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
GeometryData::IntegrationMethod CompressibleNavierStokesExplicit<3,4>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

template <>
void CompressibleNavierStokesExplicit<3,4>::EquationIdVector(
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
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Z, mom_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<3,4>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (ElementalDofList.size() != DofSize) {
        ElementalDofList.resize(DofSize);
    }

    unsigned int local_index = 0;
    const auto &r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Z, mom_pos + 2);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <>
array_1d<double,3> CompressibleNavierStokesExplicit<3,4>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GetIntegrationMethod());
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= static_cast<double>(NumNodes);
    midpoint_mom /= static_cast<double>(NumNodes);

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double rho_pow = std::pow(midpoint_rho, 2);
    const double dvz_dy = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy) / rho_pow;
    const double dvy_dz = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz) / rho_pow;
    const double dvx_dz = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz) / rho_pow;
    const double dvz_dx = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx) / rho_pow;
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / rho_pow;
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / rho_pow;
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = dvz_dy - dvy_dz;
    midpoint_rot_v[1] = dvx_dz - dvz_dx;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<3, 4>::CalculateMidPointVelocityGradient() const
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
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_dmz_dz = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
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
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_dmz_dz += r_mom[2] * node_dNdX[2];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= static_cast<double>(NumNodes);
    midpoint_mom /= static_cast<double>(NumNodes);

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v;
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(0,2) = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v(1,2) = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz);
    midpoint_grad_v(2,0) = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx);
    midpoint_grad_v(2,1) = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy);
    midpoint_grad_v(2,2) = (midpoint_dmz_dz * midpoint_rho - midpoint_mom[2] * midpoint_rho_dz);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, Dim*NumNodes> mom_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_mom_proj_3D

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
void CompressibleNavierStokesExplicit<3,4>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, 4> rho_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_rho_proj_3D

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
void CompressibleNavierStokesExplicit<3,4>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, 4> tot_ener_proj;
    const auto& DN_DX = data.DN_DX;

    //substitute_tot_ener_proj_3D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    tot_ener_proj *= data.volume / static_cast<double>(NumNodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION), tot_ener_proj[i_node]);
    }

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokesExplicit<3,4>::CalculateRightHandSideInternal(
    BoundedVector<double, 20> &rRightHandSideBoundedVector,
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
        //substitute_rhs_3D_OSS
    }
    else
    {
        //substitute_rhs_3D_ASGS
    }

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(NumNodes);

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{

    // Initialize and fill the mass matrix values
    constexpr double one_ten = 0.1;
    constexpr double one_twenty = 0.05;
    rMassMatrix = ZeroMatrix(DofSize, DofSize);
    rMassMatrix(0, 0) = one_ten; rMassMatrix(0, 5) = one_twenty; rMassMatrix(0, 10) = one_twenty; rMassMatrix(0,15) = one_twenty;
    rMassMatrix(1, 1) = one_ten; rMassMatrix(1, 6) = one_twenty; rMassMatrix(1, 11) = one_twenty; rMassMatrix(1,16) = one_twenty;
    rMassMatrix(2, 2) = one_ten; rMassMatrix(2, 7) = one_twenty; rMassMatrix(2, 12) = one_twenty; rMassMatrix(2,17) = one_twenty;
    rMassMatrix(3, 3) = one_ten; rMassMatrix(3, 8) = one_twenty; rMassMatrix(3, 13) = one_twenty; rMassMatrix(3,18) = one_twenty;
    rMassMatrix(4, 4) = one_ten; rMassMatrix(4, 9) = one_twenty; rMassMatrix(4, 14) = one_twenty; rMassMatrix(4,19) = one_twenty;
    rMassMatrix(5, 0) = one_twenty; rMassMatrix(5, 5) = one_ten; rMassMatrix(5, 10) = one_twenty; rMassMatrix(5,15) = one_twenty;
    rMassMatrix(6, 1) = one_twenty; rMassMatrix(6, 6) = one_ten; rMassMatrix(6, 11) = one_twenty; rMassMatrix(6,16) = one_twenty;
    rMassMatrix(7, 2) = one_twenty; rMassMatrix(7, 7) = one_ten; rMassMatrix(7, 12) = one_twenty; rMassMatrix(7,17) = one_twenty;
    rMassMatrix(8, 3) = one_twenty; rMassMatrix(8, 8) = one_ten; rMassMatrix(8, 13) = one_twenty; rMassMatrix(8,18) = one_twenty;
    rMassMatrix(9, 4) = one_twenty; rMassMatrix(9, 9) = one_ten; rMassMatrix(9, 14) = one_twenty; rMassMatrix(9,19) = one_twenty;
    rMassMatrix(10, 0) = one_twenty; rMassMatrix(10, 5) = one_twenty; rMassMatrix(10, 10) = one_ten; rMassMatrix(10,15) = one_twenty;
    rMassMatrix(11, 1) = one_twenty; rMassMatrix(11, 6) = one_twenty; rMassMatrix(11, 11) = one_ten; rMassMatrix(11,16) = one_twenty;
    rMassMatrix(12, 2) = one_twenty; rMassMatrix(12, 7) = one_twenty; rMassMatrix(12, 12) = one_ten; rMassMatrix(12,17) = one_twenty;
    rMassMatrix(13, 3) = one_twenty; rMassMatrix(13, 8) = one_twenty; rMassMatrix(13, 13) = one_ten; rMassMatrix(13,18) = one_twenty;
    rMassMatrix(14, 4) = one_twenty; rMassMatrix(14, 9) = one_twenty; rMassMatrix(14, 14) = one_ten; rMassMatrix(14,19) = one_twenty;
    rMassMatrix(15, 0) = one_twenty; rMassMatrix(15, 5) = one_twenty; rMassMatrix(15, 10) = one_twenty; rMassMatrix(15,15) = one_ten;
    rMassMatrix(16, 1) = one_twenty; rMassMatrix(16, 6) = one_twenty; rMassMatrix(16, 11) = one_twenty; rMassMatrix(16,16) = one_ten;
    rMassMatrix(17, 2) = one_twenty; rMassMatrix(17, 7) = one_twenty; rMassMatrix(17, 12) = one_twenty; rMassMatrix(17,17) = one_ten;
    rMassMatrix(18, 3) = one_twenty; rMassMatrix(18, 8) = one_twenty; rMassMatrix(18, 13) = one_twenty; rMassMatrix(18,18) = one_ten;
    rMassMatrix(19, 4) = one_twenty; rMassMatrix(19, 9) = one_twenty; rMassMatrix(19, 14) = one_twenty; rMassMatrix(19,19) = one_ten;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CompressibleNavierStokesExplicit<3,4>;

}