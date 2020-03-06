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
#include "utilities/math_utils.h"
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
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_SCALAR1)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_VECTOR1)
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, node)
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
    KRATOS_TRY

    // Resize of the Left and Right Hand side
    if(rLeftHandSideMatrix.size1() != 9)
        rLeftHandSideMatrix.resize(9, 9, false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != 9)
        rRightHandSideVector.resize(9, false);  // False says not to preserve existing storage!!

    // Struct to pass around the data
    ElementData data;
    data.InitializeData(rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    double length = GetGeometry().Length();

    // Mass matrix
    BoundedMatrix<double,9,9> vel_mass_matrix = ZeroMatrix(9,9);
    BoundedMatrix<double,9,9> h_mass_matrix = ZeroMatrix(9,9);
    ComputeMassMatrices(data, vel_mass_matrix, h_mass_matrix);

    // main loop (one Gauss point)

    // Build shape functions and derivatives at the Gauss point
    BoundedMatrix<double,2,9> N_vel = ZeroMatrix(2,9);        // Shape functions matrix (for velocity)
    array_1d<double,9> N_height = ZeroVector(9);              // Shape functions vector (for height)
    array_1d<double,9> DN_DX_vel = ZeroVector(9);             // Gradients vector (for velocity)
    BoundedMatrix<double,2,9> DN_DX_height = ZeroMatrix(2,9); // Gradients matrix (for height)
    for(unsigned int node = 0; node < 3; node++)
    {
        // Velocity divergence
        DN_DX_vel[  node*3] = DN_DX(node,0);
        DN_DX_vel[1+node*3] = DN_DX(node,1);
        // Height gradient
        DN_DX_height(0, 2+node*3) = DN_DX(node,0);
        DN_DX_height(1, 2+node*3) = DN_DX(node,1);
        // Height shape funtions
        N_height[2+node*3] = N[node];
        // Velocity shape functions
        N_vel(0,   node*3) = N[node];
        N_vel(1, 1+node*3) = N[node];
    }

    data.GetNodalData(GetGeometry());

    // Build LHS
    // Wave equation terms
    BoundedMatrix<double,9,9> vel_wave = prod(trans(N_vel),DN_DX_height);
    noalias(rLeftHandSideMatrix)  = data.gravity * vel_wave;           // Add <w*g*grad(h)> to Momentum Eq
    noalias(rLeftHandSideMatrix) += data.height * outer_prod(N_height, DN_DX_vel); // Add <q*h*div(u)> to Mass Eq

    // Inertia terms
    noalias(rLeftHandSideMatrix) += data.dt_inv * vel_mass_matrix; // Velocity lumped mass matrix
    noalias(rLeftHandSideMatrix) += data.dt_inv * h_mass_matrix;   // Height lumped mass matrix

    // Friction term
    const double abs_vel = norm_2(data.velocity);
    const double height4_3 = std::pow(data.height, 1.33333333333) + data.irregularity;
    noalias(rLeftHandSideMatrix) += data.gravity * data.manning2 * abs_vel / height4_3 * vel_mass_matrix;

    // Build RHS
    // Source terms (bathymetry contribution)
    noalias(rRightHandSideVector)  = data.gravity * prod(vel_wave, data.depth);

    // Inertia terms
    noalias(rRightHandSideVector) += data.dt_inv * prod(vel_mass_matrix, data.prev_unk);
    noalias(rRightHandSideVector) += data.dt_inv * prod(h_mass_matrix, data.prev_unk);

    // Computing stabilization terms
    double art = length * data.stab_factor;

    // Mass balance LHS stabilization terms
    BoundedMatrix<double,9,9> diff_h = prod(trans(DN_DX_height), DN_DX_height);
    noalias(rLeftHandSideMatrix) += art * diff_h; // Second order FIC shock capturing
    noalias(rRightHandSideVector) += art * prod(diff_h, data.depth); // Substracting the bottom diffusion

    // Momentum balance stabilization terms
    BoundedMatrix<double,9,9> diff_v = outer_prod(DN_DX_vel, DN_DX_vel);
    noalias(rLeftHandSideMatrix) += art * diff_v; // Second order FIC

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, data.unknown);

    rRightHandSideVector *= area;
    rLeftHandSideMatrix *= area;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowWater2D3::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowWater2D3: CalculateLeftHandSide not implemented" << std::endl;
}

//----------------------------------------------------------------------

void ShallowWater2D3::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowWater2D3: CalculateRightHandSide not implemented" << std::endl;
}

//----------------------------------------------------------------------

void ShallowWater2D3::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VEL_ART_VISC ||
        rVariable == PR_ART_VISC ||
        rVariable == RESIDUAL_NORM ||
        rVariable == MIU)
    {
        for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
            rValues[PointNumber] = double(this->GetValue(rVariable));
    }
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

void ShallowWater2D3::ElementData::GetNodalData(const GeometryType& rGeometry)
{
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

void ShallowWater2D3::ComputeMassMatrices(
    const ElementData& rData,
    BoundedMatrix<double,9,9>& rVelMatrix,
    BoundedMatrix<double,9,9>& rHeightMatrix)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        rVelMatrix(3*i  ,3*i  ) = 1;
        rVelMatrix(3*i+1,3*i+1) = 1;
        rHeightMatrix(3*i+2,3*i+2) = 1;
    }
    rVelMatrix *= rData.lumping_factor;
    rHeightMatrix *= rData.lumping_factor;
}

} // namespace kratos
