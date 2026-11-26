//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez, Ruben Zorrilla
//

//automated_message

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"


namespace Kratos {

template<>
GeometryData::IntegrationMethod CompressibleNavierStokesExplicit<2,4>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <>
void CompressibleNavierStokesExplicit<2,4>::EquationIdVector(
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
void CompressibleNavierStokesExplicit<2,4>::GetDofList(
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


void ComputeMidpointShapeFunctions(
    const Geometry<Node>& rGeometry,
    Vector& rN,
    Matrix& rDN_DX)
{
    Matrix DN_De;
    Matrix Jinv;

    const auto center = rGeometry.Center();
    rGeometry.ShapeFunctionsValues(rN, center);
    rGeometry.InverseOfJacobian(Jinv, center);

    rGeometry.ShapeFunctionsLocalGradients(DN_De, center);
    GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, rDN_DX);
}


template <>
array_1d<double,3> CompressibleNavierStokesExplicit<2,4>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geometry = GetGeometry();

    Vector N;
    Matrix DN_DX;

    ComputeMidpointShapeFunctions(r_geometry, N, DN_DX);

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double, 3> midpoint_mom = ZeroVector(3);

    for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
        auto& r_node = r_geometry[i_node];

        const auto node_dNdX = row(DN_DX, i_node);
        const auto node_N = N(i_node);

        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += node_N * r_rho;
        midpoint_mom += node_N * r_mom;
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }

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
BoundedMatrix<double, 3, 3> CompressibleNavierStokesExplicit<2, 4>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geometry = GetGeometry();

    Vector N;
    Matrix DN_DX;

    ComputeMidpointShapeFunctions(r_geometry, N, DN_DX);

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
        auto& r_node = r_geometry[i_node];

        const auto node_dNdX = row(DN_DX, i_node);
        const auto node_N = N(i_node);

        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += node_N * r_rho;
        midpoint_mom += node_N * r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }

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
void CompressibleNavierStokesExplicit<2,4>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    constexpr std::size_t vector_size = Dim*NumNodes;
    BoundedVector<double, vector_size> mom_proj = ZeroVector(vector_size);

    Vector N;
    Matrix DN_De;
    Matrix DN_DX;

    Matrix J;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for(const auto& gauss_point: gauss_points)
    {
        const double w = gauss_point.Weight();
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

        double detJ;
        r_geometry.Jacobian(J, gauss_point.Coordinates());
        MathUtils<double>::InvertMatrix(J, Jinv, detJ);

        r_geometry.ShapeFunctionsLocalGradients(DN_De, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, DN_DX);

        BoundedVector<double, vector_size> mom_proj_gauss(DofSize);

        //substitute_mom_proj_2D

        mom_proj += w * detJ * mom_proj_gauss;
    }

    // Assembly the projection contributions
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
void CompressibleNavierStokesExplicit<2,4>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, NumNodes> rho_proj = ZeroVector(NumNodes);

    Vector N;
    Matrix DN_De;
    Matrix DN_DX;

    Matrix J;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for(const auto& gauss_point: gauss_points)
    {
        const double w = gauss_point.Weight();
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

        double detJ;
        r_geometry.Jacobian(J, gauss_point.Coordinates());
        MathUtils<double>::InvertMatrix(J, Jinv, detJ);

        r_geometry.ShapeFunctionsLocalGradients(DN_De, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, DN_DX);

        BoundedVector<double, DofSize> rho_proj_gauss(DofSize);

        //substitute_rho_proj_2D

        rho_proj += w * detJ * rho_proj_gauss;
    }

    // Assembly the projection contributions
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(DENSITY_PROJECTION), rho_proj[i_node]);
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Calculate shock capturing values
    BoundedVector<double, NumNodes> tot_ener_proj = ZeroVector(NumNodes);

    Vector N;
    Matrix DN_De;
    Matrix DN_DX;

    Matrix J;
    Matrix Jinv;

    auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for(const auto& gauss_point: gauss_points)
    {
        const double w = gauss_point.Weight();
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

        double detJ;
        r_geometry.Jacobian(J, gauss_point.Coordinates());
        MathUtils<double>::InvertMatrix(J, Jinv, detJ);

        r_geometry.ShapeFunctionsLocalGradients(DN_De, gauss_point.Coordinates());
        GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, DN_DX);

        BoundedVector<double, DofSize> tot_ener_proj_gauss(DofSize);

        //substitute_tot_ener_proj_2D

        tot_ener_proj += w * detJ * tot_ener_proj_gauss;
    }

    // Assembly the projection contributions
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        AtomicAdd(r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION), tot_ener_proj[i_node]);
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateRightHandSideInternal(
    BoundedVector<double, DofSize> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    std::fill(rRightHandSideBoundedVector.begin(), rRightHandSideBoundedVector.end(), 0.0);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Stabilization parameters
    constexpr double stab_c1 = 12.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

    const auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    Vector N;
    Matrix DN_De;
    Matrix DN_DX;

    Matrix J;
    Matrix Jinv;

    if (data.UseOSS)
    {
        for(const auto& gauss_point: gauss_points)
        {
            const double w = gauss_point.Weight();
            r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

            double detJ;
            r_geometry.Jacobian(J, gauss_point.Coordinates());
            MathUtils<double>::InvertMatrix(J, Jinv, detJ);

            r_geometry.ShapeFunctionsLocalGradients(DN_De, gauss_point.Coordinates());
            GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, DN_DX);

            BoundedVector<double, DofSize> rhs_gauss(DofSize);

            //substitute_rhs_2D_OSS

            rRightHandSideBoundedVector += w * detJ * rhs_gauss;
        }
    }
    else
    {
        for(const auto& gauss_point: gauss_points)
        {

            const double w = gauss_point.Weight();
            r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

            double detJ;
            r_geometry.Jacobian(J, gauss_point.Coordinates());
            MathUtils<double>::InvertMatrix(J, Jinv, detJ);

            r_geometry.ShapeFunctionsLocalGradients(DN_De, gauss_point.Coordinates());
            GeometryUtils::ShapeFunctionsGradients(DN_De, Jinv, DN_DX);

            BoundedVector<double, DofSize> rhs_gauss(DofSize);

            //substitute_rhs_2D_ASGS

            rRightHandSideBoundedVector += w * detJ * rhs_gauss;
        }
    }

    KRATOS_CATCH("")
}


template <>
void CompressibleNavierStokesExplicit<2,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    Vector N;
    Matrix J;
    Matrix Jinv;

    c_matrix<double, NumNodes, NumNodes> M = ZeroMatrix(NumNodes, NumNodes);
    for(const auto& gauss_point: gauss_points)
    {
        const double w = gauss_point.Weight();
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());

        double detJ;
        r_geometry.Jacobian(J, gauss_point.Coordinates());
        MathUtils<double>::InvertMatrix(J, Jinv, detJ);

        noalias(M) += outer_prod(N, N) * detJ * w;
    }

    // Distributing 4x4 matrix to 16x16 matrix
    rMassMatrix = ZeroMatrix(DofSize, DofSize);
    for(std::size_t i_dof=0; i_dof < BlockSize; ++i_dof)
    {
        for(std::size_t i_node=0; i_node < NumNodes; ++i_node)
        {
            const std::size_t I = i_dof*BlockSize + i_node;
            for(std::size_t j_node=0; j_node < NumNodes; ++j_node)
            {
                const std::size_t J = i_dof*BlockSize + j_node;
                rMassMatrix(I, J) = M(i_node, j_node);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2,4>;

}
