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

    // Compute tau
    this->CalculateTau(rVariables);

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
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {
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
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto qstau = rVariables.qstau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double DN_DX_2_1 = rVariables.DN_DX(2, 1);
    // LHS and RHS
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_lhs_2D

    //substitute_rhs_2D
    
    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/3;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/3;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto qstau = rVariables.qstau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double DN_DX_3_2 = rVariables.DN_DX(3,2);
    // LHS and RHS
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_lhs_3D

    //substitute_rhs_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/4;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/4;
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
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double DN_DX_2_1 = rVariables.DN_DX(2, 1);
    // LHS and RHS
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_oss_2D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    noalias(rRightHandSideVector) += rhs * rVariables.volume/3;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeOSSGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    const auto phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double DN_DX_3_2 = rVariables.DN_DX(3,2);
    // LHS and RHS
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    //substitute_oss_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    noalias(rRightHandSideVector) += rhs * rVariables.volume/4;
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
    const auto DN = rVariables.DN_DX;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau[g];
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double DN_DX_2_1 = rVariables.DN_DX(2, 1);

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1)); // convective term 2
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
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau[g];
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double DN_DX_3_2 = rVariables.DN_DX(3,2);
    
    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/(RK_time_coefficient*delta_time); // mass term
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)); // convective term 1
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2)); // convective term 2
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
    double h = this->ComputeH(rVariables.DN_DX);
    // Calculate tau and qstau for each gauss point
    const auto& r_geometry = this->GetGeometry();
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    for(unsigned int g = 0; g<TNumNodes; g++)
    {
	const auto& N = row(N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
        noalias(vel_gauss) = prod(N,rVariables.convective_velocity);
        double div_vel = 0;
        for(unsigned int node_element = 0; node_element<TNumNodes; node_element++)
        {
            for(unsigned int dim = 0; dim < TDim; dim++)
            {
                div_vel += rVariables.DN_DX(node_element,dim)*rVariables.convective_velocity(node_element,dim);
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
        rVariables.tau[g] = (rVariables.density*rVariables.specific_heat) / inv_tau;
        // Estimate quasi-static tau
        double inv_qstau = 0;
        // Convection
        inv_qstau += 2.0 * norm_velocity / h;
        inv_qstau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
        // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
        inv_qstau *= rVariables.density * rVariables.specific_heat;
        // Diffusion
        inv_qstau += 4.0 * rVariables.diffusivity / (h*h);
        // Limiting
        inv_qstau = std::max(inv_qstau, 1e-2);
        rVariables.qstau[g] = (rVariables.density*rVariables.specific_heat) / inv_qstau;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicDynamicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicDynamicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
