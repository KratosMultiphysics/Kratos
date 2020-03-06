//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "shallow_water_2d_3.h"

namespace Kratos
{

int ShallowWater2D3::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(MOMENTUM)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
    KRATOS_CHECK_VARIABLE_KEY(BATHYMETRY)
    KRATOS_CHECK_VARIABLE_KEY(RAIN)
    KRATOS_CHECK_VARIABLE_KEY(MANNING)
    KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME)
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < 3; i++ )
    {
        Node<3>& node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return ierr;

    KRATOS_CATCH("")
}

void ShallowWater2D3::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rResult.size() != 9)
        rResult.resize(9, false); // False says not to preserve existing storage!!

    const GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

void ShallowWater2D3::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rElementalDofList.size() != 9)
        rElementalDofList.resize(9);

    const GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < 3; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}

void ShallowWater2D3::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // Resize and initialize the Left and Right Hand side
    if(rLeftHandSideMatrix.size1() != 9)
        rLeftHandSideMatrix.resize(9, 9, false);

    if(rRightHandSideVector.size() != 9)
        rRightHandSideVector.resize(9, false);

    rLeftHandSideMatrix = ZeroMatrix(9,9);
    rRightHandSideVector = ZeroVector(9);

    // Struct to pass around the data
    ElementData data;
    data.InitializeData(rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    data.GetNodalData(GetGeometry(), DN_DX);

    AddInertiaTerms(rLeftHandSideMatrix, rRightHandSideVector, data);

    AddGradientTerms(rLeftHandSideMatrix, rRightHandSideVector, data, N, DN_DX);

    AddSourceTerms(rLeftHandSideMatrix, rRightHandSideVector, data, N, DN_DX);

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, data.unknown);

    rRightHandSideVector *= area;
    rLeftHandSideMatrix *= area;
}

void ShallowWater2D3::AddInertiaTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData)
{
    BoundedMatrix<double,9,9> mass_matrix = ZeroMatrix(9,9);
    ComputeMassMatrix(mass_matrix, rData);
    rLHS += rData.dt_inv * (mass_matrix);
    rRHS += rData.dt_inv * prod(mass_matrix, rData.prev_unk);
}

void ShallowWater2D3::AddGradientTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    BoundedMatrix<double,9,9> gradient_matrix = ZeroMatrix(9,9);
    ComputeGradientMatrix(gradient_matrix, rData, rN, rDN_DX);
    rLHS += gradient_matrix;
    rRHS += prod(gradient_matrix, rData.depth);
}

void ShallowWater2D3::AddSourceTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    BoundedMatrix<double,9,9> mass_matrix = ZeroMatrix(9,9);
    ComputeMassMatrix(mass_matrix, rData);

    // Friction term
    const double abs_vel = norm_2(rData.velocity);
    const double height4_3 = std::pow(rData.effective_height, 1.33333333333) + rData.irregularity;
    rLHS += rData.gravity * rData.manning2 * abs_vel / height4_3 * mass_matrix;

    // Advection term
    rLHS += rData.velocity_div * mass_matrix;

    // Rain
    rRHS += prod(mass_matrix, rData.rain);
}

void ShallowWater2D3::AddShockCapturingTerm(
    MatrixType& rLHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    BoundedMatrix<double,9,9> diff_matrix = ZeroMatrix(9,9);
    ComputeDiffusionMatrix(diff_matrix, rData, rN, rDN_DX);
    rLHS += diff_matrix;
}

void ShallowWater2D3::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowWater2D3: CalculateLeftHandSide not implemented" << std::endl;
}

void ShallowWater2D3::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowWater2D3: CalculateRightHandSide not implemented" << std::endl;
}

void ShallowWater2D3::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

void ShallowWater2D3::ElementData::InitializeData(const ProcessInfo& rCurrentProcessInfo)
{
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    dt_inv = 1.0 / delta_t;
    lumping_factor = 1.0 / 3.0;;
    stab_factor = rCurrentProcessInfo[DYNAMIC_TAU];
    gravity = rCurrentProcessInfo[GRAVITY_Z];
    irregularity = 1e-2;
}

void ShallowWater2D3::ElementData::GetNodalData(const GeometryType& rGeometry, const BoundedMatrix<double,3,2>& rDN_DX)
{
    height = 0.0;
    flow_rate = ZeroVector(3);
    velocity = ZeroVector(3);
    velocity_div = 0.0;
    manning2 = 0.0;
    wet_fraction = 0.0;
    effective_height = 0.0;

    unsigned int j = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        const auto h = rGeometry[i].FastGetSolutionStepValue(HEIGHT);
        const auto f = rGeometry[i].FastGetSolutionStepValue(MOMENTUM);
        const auto v = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
        const auto hn = rGeometry[i].FastGetSolutionStepValue(HEIGHT, 1);
        const auto fn = rGeometry[i].FastGetSolutionStepValue(MOMENTUM, 1);
        const auto n = rGeometry[i].FastGetSolutionStepValue(EQUIVALENT_MANNING);

        height += h;
        flow_rate += f;
        velocity += v;
        velocity_div += v[0] * rDN_DX(i,0);
        velocity_div += v[1] * rDN_DX(i,1);
        manning2 += n;

        depth[j] = 0;
        rain[j]  = 0;
        unknown[j]  = f[0];
        prev_unk[j]  = fn[0];
        j++;

        depth[j] = 0;
        rain[j]  = 0;
        unknown[j]  = f[1];
        prev_unk[j]  = fn[1];
        j++;

        depth[j] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY);
        rain[j]  = rGeometry[i].FastGetSolutionStepValue(RAIN);
        unknown[j]  = h;
        prev_unk[j]  = hn;
        j++;

        double aux_wet_fraction, aux_effective_height;
        PhaseFunctions(h, aux_wet_fraction, aux_effective_height);
        wet_fraction += aux_wet_fraction;
        effective_height += aux_effective_height;
    }
    height *= lumping_factor;
    flow_rate *= lumping_factor;
    velocity *= lumping_factor;
    manning2 *= lumping_factor;
    manning2 = std::pow(manning2, 2);
    wet_fraction *= lumping_factor;
    effective_height *= lumping_factor;
}

void ShallowWater2D3::ElementData::PhaseFunctions(double Height, double& rWetFraction, double& rEffectiveHeight)
{
    const double unit_height = height / irregularity;
    rWetFraction = 0.5 * (1 + std::erf(2 * unit_height));
    rEffectiveHeight = rWetFraction * Height + irregularity * std::exp(-4 * std::pow(unit_height, 2)) / 4 / std::sqrt(M_PI);
}

void ShallowWater2D3::ComputeMassMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData)
{
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t block = 3 * i;
        rMatrix(block, block) += rData.lumping_factor;
        rMatrix(block+1, block+1) += rData.lumping_factor;
        rMatrix(block+2, block+2) += rData.lumping_factor * rData.wet_fraction;

        // TODO: add stabilization
    }
}

void ShallowWater2D3::ComputeGradientMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++i)
        {
            const size_t j_block = 3 * j;
            /* First component
               A_1 = {{u_1   0   gh},
                      { 0   u_1   0},
                      { 1    0    0}}
             */
            const double g1 = rN[i] * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     += g1 * rData.velocity[0];
            rMatrix(i_block,     j_block + 2) += g1 * rData.gravity * rData.effective_height;
            rMatrix(i_block + 1, j_block + 1) += g1 * rData.velocity[0];
            rMatrix(i_block + 2, j_block)     += g1;

            /* Second component
               A_2 = {{u_2   0    0},
                      { 0   u_2  gh},
                      { 0    1    0}}
             */
            const double g2 = rN[i] * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += g2 * rData.velocity[1];
            rMatrix(i_block + 1, j_block + 1) += g2 * rData.velocity[1];
            rMatrix(i_block + 1, j_block + 2) += g2 * rData.gravity * rData.effective_height;
            rMatrix(i_block + 2, j_block + 1) += g2;

            // TODO: add stabilization
        }
    }
}

void ShallowWater2D3::ComputeDiffusionMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    array_1d<double,2> f_residual;
    double h_residual;
    AlgebraicResidual(f_residual, h_residual, rData, rN, rDN_DX);

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++i)
        {
            const size_t j_block = 3 * j;

            const double d11 = rDN_DX(i,0) * rDN_DX(j,0);
            const double d12 = rDN_DX(i,0) * rDN_DX(j,1);
            const double d21 = rDN_DX(i,1) * rDN_DX(j,0);
            const double d22 = rDN_DX(i,1) * rDN_DX(j,1);

            // TODO: do something
            rMatrix(i_block,     j_block)     += d11;
        }
    }
}

void ShallowWater2D3::AlgebraicResidual(
        array_1d<double,2>& rFlowResidual,
        double& rHeightresidual,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX)
{
    // TODO: compute the residual
}

} // namespace kratos
