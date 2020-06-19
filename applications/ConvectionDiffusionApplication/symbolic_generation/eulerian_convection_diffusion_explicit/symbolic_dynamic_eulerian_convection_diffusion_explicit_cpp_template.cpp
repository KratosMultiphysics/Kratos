//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#include "symbolic_dynamic_eulerian_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDynamicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDynamicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicDynamicEulerianConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDynamicEulerianConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDynamicEulerianConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != local_size)
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    if (rRightHandSideVector.size() != local_size)
        rRightHandSideVector.resize(local_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Define local variables
    Element::GeometryType::JacobiansType J0;
    Matrix InvJ0(dimension,dimension);
    double DetJ0;

    // Compute Jacobian
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g],InvJ0,DetJ0);

        // Calculate the cartesian derivatives on integration point "g"
        rVariables.DN = prod(DN_De[g],InvJ0);
        // Caluclate N on the gauss point "g"
        rVariables.N = row(N_gausspoint,g);
        // Compute weight
        rVariables.weight = integration_points[g].Weight() * DetJ0;
        // Compute tau
        this->CalculateTau(rVariables);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        if (rCurrentProcessInfo.GetValue(RUNGE_KUTTA_STEP)==1)
        {
            rVariables.unknown_subscale = 0; // because temporal derivative is zero, and delta time is zero
        }
        else
        {
            rVariables.unknown_subscale = mUnknownSubScale(g);
        }

        // Execute standard RHS-LHS build or OSS step
        if (rCurrentProcessInfo.GetValue(OSS_SWITCH) == 1)
        {
            // Update OSS additional term
            this->ComputeOSSGaussPointContribution(rVariables,rLeftHandSideMatrix,rRightHandSideVector);
        }
        else
        {
            // Update rhs and lhs
            this->ComputeGaussPointContribution(rVariables,rLeftHandSideMatrix,rRightHandSideVector);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize output
    if (mUnknownSubScale.size() != 3) // three integration points
        mUnknownSubScale.resize(3, false);
    mUnknownSubScale = ZeroVector(3);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::FinalizeSolutionStep(
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);
    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    // Define local variables
    Element::GeometryType::JacobiansType J0;
    Matrix InvJ0(dimension,dimension);
    double DetJ0;
    // Compute Jacobian
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g],InvJ0,DetJ0);

        // Calculate the cartesian derivatives on integration point "g"
        rVariables.DN = prod(DN_De[g],InvJ0);
        // Caluclate N on the gauss point "g"
        rVariables.N = row(N_gausspoint,g);
        // Compute tau
        this->CalculateTau(rVariables);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rVariables.unknown_subscale = mUnknownSubScale(g);

        // Update unknown belonging to subgrid scale space on gauss integration point g
        this->UpdateUnknownSubgridScaleGaussPoint(rVariables,g);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
int SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::ComputeGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_lhs_2D

    //substitute_rhs_2D

    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_lhs_3D

    //substitute_rhs_3D

    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::ComputeOSSGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_oss_2D

    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeOSSGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_oss_3D

    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/(RK_time_coefficient*delta_time); // mass term
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)); // convective term 1
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
    ElementVariables& rVariables)
{
    // Calculate h
    double h = this->ComputeH(rVariables.DN);
    // Calculate velocity and velocity divergence in the gauss point
    array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
    double div_vel = 0;
    for(unsigned int node_element = 0; node_element<TNumNodes; node_element++)
    {
        for(unsigned int dim = 0; dim < TDim; dim++)
        {
            noalias(vel_gauss) = prod(rVariables.N,rVariables.convective_velocity);
            div_vel += rVariables.DN(node_element,dim)*rVariables.convective_velocity(node_element,dim);
        }
    }
    const double norm_velocity = norm_2(vel_gauss);
    // Estimate tau
    double inv_tau = 0;
    // Dynamic part
    inv_tau += 1.0/rVariables.delta_time;
    // Convection
    inv_tau += 2.0 * norm_velocity / h;
    inv_tau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
    // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
    inv_tau *= rVariables.density * rVariables.specific_heat;
    // Diffusion
    inv_tau += 4.0 * rVariables.diffusivity / (h*h);
    // Limiting
    inv_tau = std::max(inv_tau, 1e-2);
    rVariables.tau = (rVariables.density*rVariables.specific_heat) / inv_tau;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicDynamicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicDynamicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
