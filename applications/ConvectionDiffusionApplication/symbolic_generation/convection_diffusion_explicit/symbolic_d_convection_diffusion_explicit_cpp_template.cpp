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

#include "symbolic_d_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicDConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();

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

    // Execute RHS-LHS build
    this->CalculateLocalSystemInternal(rVariables,rLeftHandSideMatrix,rRightHandSideVector);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    Matrix LeftHandSide;
    this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    VectorType RightHandSide;
    this->CalculateLocalSystem(rLeftHandSideMatrix,RightHandSide,rCurrentProcessInfo);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    VectorType rhs;
    this->CalculateRightHandSide(rhs,rCurrentProcessInfo);
    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable()) += rhs[i_node];
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);
    // Resize and intialize output
    if (mUnknownSubScale.size() != TNumNodes) // number integration points = number nodes
        mUnknownSubScale.resize(TNumNodes, false);
    mUnknownSubScale = ZeroVector(TNumNodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::FinalizeSolutionStep(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);
    // Reading integration points and local gradients
    const auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Iterate over integration points to update subscales on the gauss point
    for (unsigned int g = 0; g < integration_points.size(); g++) {
        // Caluclate N on the gauss point "g"
        rVariables.N = row(rVariables.N_gausspoint,g);
        // Compute tau
        this->CalculateTau(rVariables);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rVariables.unknown_subscale = mUnknownSubScale(g);
        // Update unknown belonging to subgrid scale space on gauss integration point g
        this->UpdateUnknownSubgridScaleGaussPoint(rVariables,g);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    VectorType rhs_oss;
    this->CalculateOrthogonalSubgridScaleSystem(rhs_oss,rCurrentProcessInfo);
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateOrthogonalSubgridScaleSystem(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();

    // Resize and intialize output
    if (rRightHandSideVector.size() != local_size)
        rRightHandSideVector.resize(local_size, false);
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rVariables);

    // Execute OSS step
    this->CalculateOrthogonalSubgridScaleSystemInternal(rVariables,rRightHandSideVector);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
int SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
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
void SymbolicDConvectionDiffusionExplicit<2>::CalculateLocalSystemInternal(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau;
    const auto& prj = rVariables.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double& DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double& DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2, 1);
    // LHS and RHS
    auto& lhs = rVariables.lhs;
    auto& rhs = rVariables.rhs;

    //substitute_lhs_2D

    //substitute_rhs_2D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3>::CalculateLocalSystemInternal(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau;
    const auto& prj = rVariables.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double& DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double& DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double& DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double& DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double& DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double& DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double& DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double& DN_DX_3_2 = rVariables.DN_DX(3,2);
    // LHS and RHS
    auto& lhs = rVariables.lhs;
    auto& rhs = rVariables.rhs;

    //substitute_lhs_3D

    //substitute_rhs_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<2>::CalculateOrthogonalSubgridScaleSystemInternal(
    ElementVariables& rVariables,
    VectorType& rRightHandSideVector)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double& DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double& DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2, 1);
    // RHS
    auto& rhs = rVariables.rhs;

    //substitute_oss_2D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3>::CalculateOrthogonalSubgridScaleSystemInternal(
    ElementVariables& rVariables,
    VectorType& rRightHandSideVector)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double& DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double& DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double& DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double& DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double& DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double& DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double& DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double& DN_DX_3_2 = rVariables.DN_DX(3,2);
    // RHS
    auto& rhs = rVariables.rhs;

    //substitute_oss_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<2>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& N = rVariables.N;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau[g];
    const auto& phi_subscale_gauss = rVariables.unknown_subscale;
    const auto& prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rVariables.DN_DX(0, 0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0, 1);
    const double& DN_DX_1_0 = rVariables.DN_DX(1, 0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1, 1);
    const double& DN_DX_2_0 = rVariables.DN_DX(2, 0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2, 1);

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/(delta_time); // mass term
    // phi_subscale_gauss_new += - (N[0]*phi_acceleration_old[0] + N[1]*phi_acceleration_old[1] + N[2]*phi_acceleration_old[2]); // mass term (use acceleration step n)
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element variables
    const auto& N = rVariables.N;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau[g];
    const auto& phi_subscale_gauss = rVariables.unknown_subscale;
    const auto& prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rVariables.DN_DX(0,0);
    const double& DN_DX_0_1 = rVariables.DN_DX(0,1);
    const double& DN_DX_0_2 = rVariables.DN_DX(0,2);
    const double& DN_DX_1_0 = rVariables.DN_DX(1,0);
    const double& DN_DX_1_1 = rVariables.DN_DX(1,1);
    const double& DN_DX_1_2 = rVariables.DN_DX(1,2);
    const double& DN_DX_2_0 = rVariables.DN_DX(2,0);
    const double& DN_DX_2_1 = rVariables.DN_DX(2,1);
    const double& DN_DX_2_2 = rVariables.DN_DX(2,2);
    const double& DN_DX_3_0 = rVariables.DN_DX(3,0);
    const double& DN_DX_3_1 = rVariables.DN_DX(3,1);
    const double& DN_DX_3_2 = rVariables.DN_DX(3,2);

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)); // convective term 1
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
    ElementVariables& rVariables)
{
    KRATOS_TRY;

    // Calculate h
    double h = this->ComputeH(rVariables.DN_DX);
    // Calculate tau for each gauss point
    for(unsigned int g = 0; g<TNumNodes; g++)
    {
	const auto& N = row(rVariables.N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        array_1d<double, 3 > vel_gauss=ZeroVector(3);
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
        // Diffusion
        inv_tau += 4.0 * rVariables.diffusivity / (h*h);
        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);
        rVariables.tau[g] = (1.0) / inv_tau;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicDConvectionDiffusionExplicit<2>;
template class SymbolicDConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
