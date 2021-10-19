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
#include "utilities/geometry_utilities.h"
#include "includes/mesh_moving_variables.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "custom_friction_laws/friction_laws_factory.h"
#include "shallow_water_2d_3.h"

namespace Kratos
{

int ShallowWater2D3::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ATMOSPHERIC_PRESSURE, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VERTICAL_VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_ACCELERATION, node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return ierr;

    KRATOS_CATCH("")
}

void ShallowWater2D3::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != 9)
        rResult.resize(9, false); // False says not to preserve existing storage!!

    const GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (size_t i = 0; i < 3; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

void ShallowWater2D3::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != 9)
        rElementalDofList.resize(9);

    const GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (size_t i = 0; i < 3; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}

void ShallowWater2D3::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != 9)
        rValues.resize(9, false);

    const GeometryType& r_geom = this->GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < 3; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(HEIGHT, Step);
    }
}

void ShallowWater2D3::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != 9)
        rValues.resize(9, false);

    const GeometryType& r_geom = this->GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < 3; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

void ShallowWater2D3::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "ShallowWater2D3: This method is not supported by the formulation" << std::endl;
}

void ShallowWater2D3::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

void ShallowWater2D3::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
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

    data.boundary_velocity = rCurrentProcessInfo[BOUNDARY_VELOCITY];
    ComputeDampingCoefficient(data.damping, rCurrentProcessInfo[ABSORBING_DISTANCE], rCurrentProcessInfo[DISSIPATION]);

    data.pBottomFriction = FrictionLawsFactory().CreateBottomFrictionLaw(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    data.pSurfaceFriction = FrictionLawsFactory().CreateSurfaceFrictionLaw(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    AddGradientTerms(rLeftHandSideMatrix, rRightHandSideVector, data, N, DN_DX);

    AddSourceTerms(rLeftHandSideMatrix, rRightHandSideVector, data, N, DN_DX);

    AddShockCapturingTerm(rLeftHandSideMatrix, data, DN_DX);

    AddDesingularizationTerm(rLeftHandSideMatrix, data);

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, data.unknown);

    rRightHandSideVector *= area;
    rLeftHandSideMatrix *= area;
}

void ShallowWater2D3::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != 9)
        rMassMatrix.resize(9, 9, false);

    // Struct to pass around the data
    ElementData data;
    data.InitializeData(rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    data.GetNodalData(GetGeometry(), DN_DX);

    BoundedMatrix<double,9,9> mass_matrix = ZeroMatrix(9,9);
    ComputeMassMatrix(mass_matrix, data, N, DN_DX);
    rMassMatrix = area * mass_matrix;
}

void ShallowWater2D3::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 9) {
        rDampingMatrix.resize(9,9,false);
    }

    rDampingMatrix = ZeroMatrix(9,9);

    const double area = this->GetGeometry().Area();
    const double one_third = 1.0 / 3.0;
    const double one_twelve = 1.0 / 12.0;
    const double g = rCurrentProcessInfo[GRAVITY_Z];

    // Computing the diffusion factor (inverse of characteristic time)
    array_1d<double,3> v = ZeroVector(3);
    double h = 0.0;
    for (auto& r_node : this->GetGeometry())
    {
        v += r_node.FastGetSolutionStepValue(VELOCITY);
        h += r_node.FastGetSolutionStepValue(HEIGHT);
    }
    h = std::max(h, 0.0);
    const double lambda = norm_2(v) + std::sqrt(g*h);
    const double c = lambda / this->GetGeometry().Length();

    // Construction of the scalar-based diffusion matrix
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            const double d = (i == j)? one_third - 2*one_twelve : -one_twelve;
            rDampingMatrix(i_block,     j_block    ) = c * area * d;
            rDampingMatrix(i_block + 1, j_block + 1) = c * area * d;
            rDampingMatrix(i_block + 2, j_block + 2) = c * area * d;
        }
    }
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
}

void ShallowWater2D3::AddSourceTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    BoundedMatrix<double,9,9> flow_mass_matrix = ZeroMatrix(9,9);
    BoundedMatrix<double,9,9> height_mass_matrix = ZeroMatrix(9,9);
    ComputeMassMatrix(flow_mass_matrix, height_mass_matrix, rData, rN, rDN_DX);

    // Topography gradient and atmospheric pressure
    array_1d<double,9> source_gradients = ZeroVector(9);
    ComputeGradientVector(source_gradients, rData, rN, rDN_DX);
    rRHS += source_gradients;

    // Friction term
    rLHS += rData.gravity * rData.pBottomFriction->CalculateLHS(rData.height, rData.velocity) * flow_mass_matrix;
    rLHS -= rData.gravity * rData.pSurfaceFriction->CalculateLHS(ZeroVector(3)) * flow_mass_matrix;

    // Newtonian damper for the absorbing boundaries
    rLHS += rData.damping * flow_mass_matrix;
    rRHS += rData.damping * rData.height * prod(flow_mass_matrix, ToNodalVector(rData.boundary_velocity));

    // Rain and mesh acceleration
    const double lumping_factor = 1.0 / 3.0;
    for (size_t i = 0; i < 3; ++i) {
        const size_t block = 3 * i;
        rRHS(block+2) += lumping_factor * rData.rain[i];
    }
    rRHS -= rData.height * prod(flow_mass_matrix, rData.mesh_acc);
}

void ShallowWater2D3::AddShockCapturingTerm(
    MatrixType& rLHS,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    double art_visc;
    double art_diff;
    ShockCapturingParameters(art_visc, art_diff, rData, rDN_DX);

    BoundedMatrix<double,9,9> visc_matrix;
    BoundedMatrix<double,9,9> diff_matrix;
    ShockCapturingViscosityMatrix(visc_matrix, art_visc, rData, rDN_DX);
    ShockCapturingDiffusionMatrix(diff_matrix, art_diff, rData, rDN_DX);
    rLHS += visc_matrix;
    rLHS += diff_matrix;
}

void ShallowWater2D3::AddDesingularizationTerm(
    MatrixType& rLHS,
    const ElementData& rData)
{
    double factor = 1e3 / GetGeometry().Length();
    factor *= 1.0 - ShallowWaterUtilities().WetFraction(rData.height, rData.dry_height);
    for (size_t i = 0; i < 3; ++i) {
        const size_t block = 3 * i;
        rLHS(block, block) += factor;
        rLHS(block + 1, block + 1) += factor;
    }
}

void ShallowWater2D3::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowWater2D3: CalculateLeftHandSide not implemented" << std::endl;
}

void ShallowWater2D3::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType aux_matrix;
    CalculateLocalSystem(aux_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void ShallowWater2D3::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    for (size_t PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

void ShallowWater2D3::ElementData::InitializeData(const ProcessInfo& rCurrentProcessInfo)
{
    stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    shock_stab_factor = rCurrentProcessInfo[SHOCK_STABILIZATION_FACTOR];
    rel_dry_height = rCurrentProcessInfo[RELATIVE_DRY_HEIGHT];
    gravity = rCurrentProcessInfo[GRAVITY_Z];
}

void ShallowWater2D3::ElementData::GetNodalData(const GeometryType& rGeometry, const BoundedMatrix<double,3,2>& rDN_DX)
{
    dry_height = rel_dry_height * rGeometry.Length();
    height = 0.0;
    flow_rate = ZeroVector(3);
    velocity = ZeroVector(3);

    for (size_t i = 0; i < 3; i++)
    {
        const size_t k = 3 * i;

        auto h = rGeometry[i].FastGetSolutionStepValue(HEIGHT);
        const auto f = rGeometry[i].FastGetSolutionStepValue(MOMENTUM);
        h = std::max(0.0, h);

        height += h;
        flow_rate += f;
        velocity += rGeometry[i].FastGetSolutionStepValue(VELOCITY);

        topography[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY);
        rain[i] = rGeometry[i].FastGetSolutionStepValue(RAIN);

        unknown[k]   = f[0];
        unknown[k+1] = f[1];
        unknown[k+2] = h;

        mesh_acc[k]   = rGeometry[i].FastGetSolutionStepValue(MESH_ACCELERATION_X);
        mesh_acc[k+1] = rGeometry[i].FastGetSolutionStepValue(MESH_ACCELERATION_Y);
        mesh_acc[k+2] = 0.0;
    }
    const double lumping_factor = 1.0 / 3.0;
    height = std::max(height, .0);
    height *= lumping_factor;
    flow_rate *= lumping_factor;
    velocity *= lumping_factor;
}

void ShallowWater2D3::ComputeMassMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // The scaling factors
    const double mu_q = 1.0;
    const double mu_h = 1.0;

    // Algebraic factor
    const double one_twelve = 1.0 / 12.0;

    // Stabilization parameters
    using std::pow;
    const double c2 = rData.gravity * rData.height; // c=sqrt(gh)
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double l = StabilizationParameter(rData);
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            // Algebraic mass matrix
            const double n = (i == j)? 2*one_twelve : one_twelve;
            rMatrix(i_block,     j_block)     += mu_q * n;
            rMatrix(i_block + 1, j_block + 1) += mu_q * n;
            rMatrix(i_block + 2, j_block + 2) += mu_h * n;

            /* Stabilization x
             * A1*M
             */
            const double s1_ij = rN[i] * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     += l * mu_q * s1_ij * 2*u_1;
            rMatrix(i_block,     j_block + 2) += l * mu_h * s1_ij * (-pow(u_1,2) + c2);
            rMatrix(i_block + 1, j_block)     += l * mu_q * s1_ij * u_2;
            rMatrix(i_block + 1, j_block + 1) += l * mu_q * s1_ij * u_1;
            rMatrix(i_block + 1, j_block + 2) -= l * mu_h * s1_ij * u_1*u_2;
            rMatrix(i_block + 2, j_block)     += l * mu_q * s1_ij;

             /* Stabilization y
              * A2*M
              */
            const double s2_ij = rN[i] * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += l * mu_q * s2_ij * u_2;
            rMatrix(i_block,     j_block + 1) += l * mu_q * s2_ij * u_1;
            rMatrix(i_block,     j_block + 2) -= l * mu_h * s2_ij * u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += l * mu_q * s2_ij * 2*u_2;
            rMatrix(i_block + 1, j_block + 2) += l * mu_h * s2_ij * (-pow(u_2,2) + c2);
            rMatrix(i_block + 2, j_block + 1) += l * mu_q * s2_ij;
        }
    }
}

void ShallowWater2D3::ComputeMassMatrix(
    BoundedMatrix<double,9,9>& rFlowMatrix,
    BoundedMatrix<double,9,9>& rHeightMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // The scaling factors
    const double mu_q = 1.0;
    const double mu_h = 1.0;

    // Algebraic factor
    const double lumping_factor = 1.0 / 3.0;
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t block = 3 * i;
        rFlowMatrix(block, block) += mu_q * lumping_factor;
        rFlowMatrix(block+1, block+1) += mu_q * lumping_factor;
        rHeightMatrix(block+2, block+2) += mu_h * lumping_factor;
    }

    // Stabilization parameters. NOTE: This stabilization term is not integrated by parts!!
    using std::pow;
    const double c2 = rData.gravity * rData.height; // c=sqrt(gh)
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double l = StabilizationParameter(rData);
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            /**
             * Note: In the source terms the consistent mass matrix
             * is skiped in favour of the lumped mass matrix.
             * It accelerates the convergence rate. The stabilization
             * terms are included to guarantee consistency
             */

            /* Stabilization x
             * A1*M
             */
            const double s1_ij = rN[j] * rDN_DX(i,0);
            rFlowMatrix  (i_block,     j_block)     += l * mu_q * s1_ij * 2*u_1;
            rHeightMatrix(i_block,     j_block + 2) += l * mu_h * s1_ij * (-pow(u_1,2) + c2);
            rFlowMatrix  (i_block + 1, j_block)     += l * mu_q * s1_ij * u_2;
            rFlowMatrix  (i_block + 1, j_block + 1) += l * mu_q * s1_ij * u_1;
            rHeightMatrix(i_block + 1, j_block + 2) -= l * mu_h * s1_ij * u_1*u_2;
            rFlowMatrix  (i_block + 2, j_block)     += l * mu_q * s1_ij;

             /* Stabilization y
              * A2*M
              */
            const double s2_ij = rN[j] * rDN_DX(i,1);
            rFlowMatrix  (i_block,     j_block)     += l * mu_q * s2_ij * u_2;
            rFlowMatrix  (i_block,     j_block + 1) += l * mu_q * s2_ij * u_1;
            rHeightMatrix(i_block,     j_block + 2) -= l * mu_h * s2_ij * u_1*u_2;
            rFlowMatrix  (i_block + 1, j_block + 1) += l * mu_q * s2_ij * 2*u_2;
            rHeightMatrix(i_block + 1, j_block + 2) += l * mu_h * s2_ij * (-pow(u_2,2) + c2);
            rFlowMatrix  (i_block + 2, j_block + 1) += l * mu_q * s2_ij;
        }
    }
}

void ShallowWater2D3::ComputeGradientMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    using std::pow;
    const double c2 = rData.gravity * rData.height; // c=sqrt(gh)
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double l = StabilizationParameter(rData);
    const double w = (rData.height > rData.dry_height)? 1.0 : 0.0;

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            /* First component
             * A_1 = {{2 * u_1   0   -u_1^2 + gh},
             *        {  u_2    u_1   -u_1 * u_2},
             *        {   1      0        0     }}
             */
            const double g1_ij = rN[i] * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     += g1_ij * w * 2*u_1;
            rMatrix(i_block,     j_block + 2) += g1_ij * w * (-u_1*u_1 + c2);
            rMatrix(i_block + 1, j_block)     += g1_ij * w * u_2;
            rMatrix(i_block + 1, j_block + 1) += g1_ij * w * u_1;
            rMatrix(i_block + 1, j_block + 2) -= g1_ij * w * u_1*u_2;
            rMatrix(i_block + 2, j_block)     += g1_ij;

            /* Second component
             * A_2 = {{u_2    u_1      -u_1 * u_2},
             *        { 0   2 * u_2   -u_2^2 + gh},
             *        { 0      1            0    }}
             */
            const double g2_ij = rN[i] * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += g2_ij * w * u_2;
            rMatrix(i_block,     j_block + 1) += g2_ij * w * u_1;
            rMatrix(i_block,     j_block + 2) -= g2_ij * w * u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += g2_ij * w * 2*u_2;
            rMatrix(i_block + 1, j_block + 2) += g2_ij * w * (-u_2*u_2 + c2);
            rMatrix(i_block + 2, j_block + 1) += g2_ij;

            double d_ij;
            /* Stabilization x-x
             * A1*A1
             */
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     += l * d_ij * w * (3*pow(u_1,2) + c2);
            rMatrix(i_block,     j_block + 2) += l * d_ij * w * (-2*pow(u_1,3) + 2*u_1*c2);
            rMatrix(i_block + 1, j_block)     += l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += l * d_ij * w * pow(u_1,2);
            rMatrix(i_block + 1, j_block + 2) += l * d_ij * w * (-2*pow(u_1,2)*u_2 + u_2*c2);
            rMatrix(i_block + 2, j_block)     += l * d_ij * w * 2*u_1;
            rMatrix(i_block + 2, j_block + 2) += l * d_ij * w * (-pow(u_1,2) + c2);

            /* Stabilization y-y
             * A2*A2
             */
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += l * d_ij * w * pow(u_2,2);
            rMatrix(i_block,     j_block + 1) += l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 2) += l * d_ij * w * (-2*u_1*pow(u_2,2) + u_1*c2);
            rMatrix(i_block + 1, j_block + 1) += l * d_ij * w * (3*(pow(u_2,2) + c2));
            rMatrix(i_block + 1, j_block + 2) += l * d_ij * w * (-2*pow(u_2,3) + 2*u_2*c2);
            rMatrix(i_block + 2, j_block + 1) += l * d_ij * w * 2*u_2;
            rMatrix(i_block + 2, j_block + 2) += l * d_ij * w * (-pow(u_2,2) + c2);

            /* Stabilization x-y
             * A1*A2
             */
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     +=  l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 1) +=  l * d_ij * w * (pow(u_1,2)+c2);
            rMatrix(i_block,     j_block + 2) += -l * d_ij * w * 2*pow(u_1,2)*u_2;
            rMatrix(i_block + 1, j_block)     +=  l * d_ij * w * pow(u_2,2);
            rMatrix(i_block + 1, j_block + 1) +=  l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 2) +=  l * d_ij * w * (-2*u_1*pow(u_2,2) + u_1*c2);
            rMatrix(i_block + 2, j_block)     +=  l * d_ij * w * u_2;
            rMatrix(i_block + 2, j_block + 1) +=  l * d_ij * w * u_1;
            rMatrix(i_block + 2, j_block + 2) += -l * d_ij * w * u_1*u_2;

            /* Stabilization y-x
             * A2*A1
             */
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     +=  l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 1) +=  l * d_ij * w * pow(u_1,2);
            rMatrix(i_block,     j_block + 2) +=  l * d_ij * w * (-2*pow(u_1,2)*u_2 + u_2*c2);
            rMatrix(i_block + 1, j_block)     +=  l * d_ij * w * (pow(u_2,2) + c2);
            rMatrix(i_block + 1, j_block + 1) +=  l * d_ij * w * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 2) += -l * d_ij * w * 2*u_1*pow(u_2,2);
            rMatrix(i_block + 2, j_block)     +=  l * d_ij * w * u_2;
            rMatrix(i_block + 2, j_block + 1) +=  l * d_ij * w * u_1;
            rMatrix(i_block + 2, j_block + 2) += -l * d_ij * w * u_1*u_2;
        }
    }
}

void ShallowWater2D3::ComputeDiffusionMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX,
    const BoundedMatrix<double,2,2>& rK1,
    const BoundedMatrix<double,2,2>& rK2,
    const BoundedMatrix<double,2,2>& rKh)
{
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            array_1d<double,2> bi, bj;
            bi[0] = rDN_DX(i,0); bi[1] = rDN_DX(i,1);
            bj[0] = rDN_DX(j,0); bj[1] = rDN_DX(j,1);

            const size_t j_block = 3 * j;

            rMatrix(i_block,     j_block)     += inner_prod(bj, prod(rK1, bi));
            rMatrix(i_block + 1, j_block + 1) += inner_prod(bj, prod(rK2, bi));
            rMatrix(i_block + 2, j_block + 2) += inner_prod(bj, prod(rKh, bi));
        }
    }
}

void ShallowWater2D3::ComputeGradientVector(
    array_1d<double,9>& rVector,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    const double c2 = rData.gravity * rData.height; // c=sqrt(gh)
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double l = StabilizationParameter(rData);
    const auto topography = rData.topography;
    const auto w = (rData.height > rData.dry_height)? 1.0 : 0.0;

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            /* First component */
            rVector[i_block] -= c2 * rN[i] * rDN_DX(j,0) * w * topography[j];

            /* Second component */
            rVector[i_block + 1] -= c2 * rN[i] * rDN_DX(j,1) * w * topography[j];

            /* Stabilization x-x
             * A1*G1
             */
            double d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            rVector[i_block]     -= l * d_ij * w * topography[j] * 2*u_1*c2;
            rVector[i_block + 1] -= l * d_ij * w * topography[j] * u_2*c2;
            rVector[i_block + 2] -= l * d_ij * w * topography[j] * c2;

            /* Stabilization y-y
             * A2*G2
             */
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            rVector[i_block]     -= l * d_ij * w * topography[j] * u_1*c2;
            rVector[i_block + 1] -= l * d_ij * w * topography[j] * 2*u_2*c2;
            rVector[i_block + 2] -= l * d_ij * w * topography[j] * c2;

            /* Stabilization x-y
             * A1*G2
             */
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            rVector[i_block + 1] -= l * d_ij * w * topography[j] * u_1*c2;

            /* Stabilization y-x
             * A2*G1
             */
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            rVector[i_block] -= l * d_ij * w * topography[j] * u_2*c2;
        }
    }
}


void ShallowWater2D3::ShockCapturingParameters(
    double& rArtViscosity,
    double& rArtDiffusion,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // Computation of the residuals and gradients
    array_1d<double,3> flow_residual;
    double height_residual;
    BoundedMatrix<double,3,3> flow_grad;
    array_1d<double,3> height_grad;
    AlgebraicResidual(flow_residual, height_residual, flow_grad, height_grad, rData, rDN_DX);

    // slope limits
    const double min_slope = 0.1;
    const double max_slope = 1.0;
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * rData.height);
    const double min_q_slope = std::max(eigenvalue * 1.0, min_slope);
    const double max_q_slope = std::max(eigenvalue * 10.0, max_slope);

    // Final assembly of the parameters
    const double length = this->GetGeometry().Length();

    const double q_residual_norm = norm_2(flow_residual);
    const double q_grad_frobenius = norm_frobenius(flow_grad);
    const double q_gradient_norm = std::min(std::max(q_grad_frobenius, min_q_slope), max_q_slope);
    rArtViscosity = 0.5 * rData.shock_stab_factor * length * q_residual_norm / q_gradient_norm;

    const double h_residual_norm = std::abs(height_residual);
    const double h_gradient_norm = std::min(std::max(norm_2(height_grad), min_slope), max_slope);
    rArtDiffusion = 0.5 * rData.shock_stab_factor * length * h_residual_norm / h_gradient_norm;
}

void ShallowWater2D3::ShockCapturingViscosityMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const double& rViscosity,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // The derivatives matrix
    BoundedMatrix<double,3,9> b = ZeroMatrix(3,9);
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        b(0, i_block)     = rDN_DX(i,0);
        b(1, i_block + 1) = rDN_DX(i,1);
        b(2, i_block)     = rDN_DX(i,1);
        b(2, i_block + 1) = rDN_DX(i,0);
    }

    // The viscosity previously added by the stabilization
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * rData.height);
    const double stab_viscosity = StabilizationParameter(rData) * std::pow(eigenvalue,2);

    // The orthogonal fourth order tensor
    BoundedMatrix<double,3,3> crosswind_tensor;
    CrossWindTensor(crosswind_tensor, rData.velocity);
    crosswind_tensor *= rViscosity;

    // The streamline fourth order tensor
    BoundedMatrix<double,3,3> streamline_tensor;
    StreamLineTensor(streamline_tensor, rData.velocity);
    streamline_tensor *= std::max(0.0, rViscosity - stab_viscosity);

    // The constitutive tensor
    BoundedMatrix<double,3,3> constitutive_tensor = IdentityMatrix(3,3);
    array_1d<double,3> m;
    m[0] = 1.0; m[1] = 1.0; m[2] = 0.0;
    constitutive_tensor -= outer_prod(m, m) / 3.0;
    constitutive_tensor = prod(constitutive_tensor, crosswind_tensor + streamline_tensor);

    // Assembly of the viscosity matrix
    BoundedMatrix<double,3,9> tmp = prod(constitutive_tensor, b);
    rMatrix = prod(trans(b), tmp);
}

void ShallowWater2D3::ShockCapturingDiffusionMatrix(
    BoundedMatrix<double,9,9>& rMatrix,
    const double& rDiffusivity,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // Output initialization
    rMatrix = ZeroMatrix(9,9);

    // The viscosity previously added by the stabilization
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * rData.height);
    const double stab_diffusivity = StabilizationParameter(rData) * std::pow(eigenvalue,2);

    // The second order crosswind tensor
    BoundedMatrix<double,2,2> crosswind_tensor;
    CrossWindTensor(crosswind_tensor, rData.velocity);
    crosswind_tensor *= rDiffusivity;

    // The second order streamline tensor
    BoundedMatrix<double,2,2> streamline_tensor;
    StreamLineTensor(streamline_tensor, rData.velocity);
    streamline_tensor *= std::max(0.0, rDiffusivity - stab_diffusivity);

    // The constitutive matrix
    BoundedMatrix<double,2,2> constitutive_matrix;
    constitutive_matrix = crosswind_tensor + streamline_tensor;

    // Assembly of the diffusion matrix
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            array_1d<double,2> bi, bj;
            bi[0] = rDN_DX(i,0); bi[1] = rDN_DX(i,1);
            bj[0] = rDN_DX(j,0); bj[1] = rDN_DX(j,1);

            const size_t j_block = 3 * j;

            rMatrix(i_block + 2, j_block + 2) = inner_prod(bj, prod(constitutive_matrix, bi));
        }
    }
}

void ShallowWater2D3::AlgebraicResidual(
    array_1d<double,3>& rFlowResidual,
    double& rHeightresidual,
    BoundedMatrix<double,3,3>& rFlowGrad,
    array_1d<double,3>& rHeightGrad,
    const ElementData& rData,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    array_1d<double,3> flow_acc = ZeroVector(3);
    double height_acc = 0.0;
    double rain = 0.0;
    double vel_div = 0.0;
    double flow_div = 0.0;
    array_1d<double,3> mesh_acc = ZeroVector(3);
    array_1d<double,3> topography_grad = ZeroVector(3);
    rFlowGrad = ZeroMatrix(3);
    rHeightGrad = ZeroVector(3);

    auto& r_geom = GetGeometry();
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t block = 3 * i;

        flow_acc += r_geom[i].FastGetSolutionStepValue(ACCELERATION);
        height_acc += r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY);
        rain += rData.rain[i];
        mesh_acc[0] += rData.mesh_acc[block];
        mesh_acc[1] += rData.mesh_acc[block+1];

        const double f1 = rData.unknown[block];
        const double f2 = rData.unknown[block + 1];
        const double h = rData.unknown[block + 2];
        const auto v = r_geom[i].FastGetSolutionStepValue(VELOCITY);

        flow_div += rDN_DX(i,0) * f1 + rDN_DX(i,1) * f2;
        vel_div += rDN_DX(i,0) * v[0] + rDN_DX(i,1) * v[1];
        rFlowGrad(0,0) += rDN_DX(i,0) * f1;
        rFlowGrad(0,1) += rDN_DX(i,0) * f2;
        rFlowGrad(1,0) += rDN_DX(i,1) * f1;
        rFlowGrad(1,1) += rDN_DX(i,1) * f2;
        rHeightGrad[0] += rDN_DX(i,0) * h;
        rHeightGrad[1] += rDN_DX(i,1) * h;
        topography_grad[0] += rDN_DX(i,0) * rData.topography[i];
        topography_grad[1] += rDN_DX(i,1) * rData.topography[i];
    }
    const double lumping_factor = 1.0 / 3.0;
    flow_acc *= lumping_factor;
    height_acc *= lumping_factor;
    rain *= lumping_factor;
    mesh_acc *= lumping_factor;

    const double c2 = rData.gravity * rData.height;
    const array_1d<double,3> friction = rData.gravity * rData.height * rData.pBottomFriction->CalculateRHS(rData.height, rData.velocity);
    const array_1d<double,3> flux = prod(rData.velocity, trans(rFlowGrad)) + vel_div * rData.flow_rate;
    const array_1d<double,3> damping = rData.damping * (rData.flow_rate - rData.height * rData.boundary_velocity);

    rFlowResidual = flow_acc + flux + c2 * (rHeightGrad + topography_grad) + friction + damping + rData.height * mesh_acc;
    rHeightresidual = height_acc + flow_div + rain;
}

void ShallowWater2D3::StreamLineTensor(BoundedMatrix<double,2,2>& rTensor, const array_1d<double,3>& rVector)
{
    const double e = std::numeric_limits<double>::epsilon(); // small value to avoid division by zero
    array_1d<double,2> aux_vector;
    aux_vector[0] = rVector[0];
    aux_vector[1] = rVector[1];
    rTensor = outer_prod(aux_vector, aux_vector) / (inner_prod(rVector, rVector) + e);
}

void ShallowWater2D3::CrossWindTensor(BoundedMatrix<double,2,2>& rTensor, const array_1d<double,3>& rVector)
{
    StreamLineTensor(rTensor, rVector);
    rTensor = IdentityMatrix(2) - rTensor;
}

void ShallowWater2D3::StreamLineTensor(BoundedMatrix<double,3,3>& rTensor, const array_1d<double,3>& rVector)
{
    BoundedMatrix<double,2,2> stream_line_tensor;
    StreamLineTensor(stream_line_tensor, rVector);
    rTensor(0,0) = stream_line_tensor(0,0);
    rTensor(1,1) = stream_line_tensor(1,1);
    rTensor(0,1) = stream_line_tensor(0,1);
    rTensor(1,0) = stream_line_tensor(0,1);
    rTensor(2,2) = stream_line_tensor(0,1);
    rTensor(0,2) = 0.0;
    rTensor(1,2) = 0.0;
    rTensor(2,0) = 0.0;
    rTensor(2,1) = 0.0;
}

void ShallowWater2D3::CrossWindTensor(BoundedMatrix<double,3,3>& rTensor, const array_1d<double,3>& rVector)
{
    StreamLineTensor(rTensor, rVector);
    rTensor = IdentityMatrix(3) - rTensor;
}

double ShallowWater2D3::StabilizationParameter(const ElementData& rData)
{
    const double e = std::numeric_limits<double>::epsilon(); // small value to avoid division by zero
    const double length = this->GetGeometry().Length();
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * rData.height);
    const double wet_fraction = ShallowWaterUtilities().WetFraction(rData.height, rData.dry_height);

    return length * wet_fraction * rData.stab_factor / (eigenvalue + e);
}

void ShallowWater2D3::ComputeDampingCoefficient(
    double& rDamping,
    const double DistanceThreshold,
    const double MaximumDamping)
{
    if (DistanceThreshold > 0.0) {
        double distance = 0.0;
        for (auto& r_node : GetGeometry()) {
            distance += r_node.FastGetSolutionStepValue(DISTANCE);
        }
        distance /= GetGeometry().size();

        if (distance < DistanceThreshold) {
            const double pow_coeff = 3.0;
            const double smooth_function = std::expm1(std::pow((DistanceThreshold - distance) / DistanceThreshold, pow_coeff)) / std::expm1(1.0);
            rDamping = MaximumDamping * smooth_function;
        } else {
            rDamping = 0.0;
        }
    } else {
        rDamping = 0.0;
    }
}

array_1d<double,9> ShallowWater2D3::ToNodalVector(const array_1d<double,3>& rVector)
{
    array_1d<double,9> nodal_vector;
    for (size_t i = 0; i < 3; ++i) {
        const size_t i_block = 3 * i;
        nodal_vector[i_block]     = rVector[0];
        nodal_vector[i_block + 1] = rVector[1];
        nodal_vector[i_block + 2] = 0.0;
    }
    return nodal_vector;
}

array_1d<double,9> ShallowWater2D3::ToNodalVector(const double& rScalar)
{
    array_1d<double,9> nodal_vector;
    for (size_t i = 0; i < 3; ++i) {
        const size_t i_block = 3 * i;
        nodal_vector[i_block]     = 0.0;
        nodal_vector[i_block + 1] = 0.0;
        nodal_vector[i_block + 2] = rScalar;
    }
    return nodal_vector;
}

} // namespace kratos
