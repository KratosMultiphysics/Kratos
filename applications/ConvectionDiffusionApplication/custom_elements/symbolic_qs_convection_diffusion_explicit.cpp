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

#include "symbolic_qs_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicQSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicQSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicQSConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicQSConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicQSConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
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
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
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
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLeftHandSide(
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
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if(rResult.size() != local_size)
        rResult.resize(local_size);
    for (unsigned int i=0; i<local_size; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);
    for (unsigned int i = 0; i < local_size; i++)
    {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = GetGeometry();
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

template <>
void SymbolicQSConvectionDiffusionExplicit<2>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int local_size = 3;
    // Resize, intialize and fill the mass matrix for linear triangular elements
    if (rMassMatrix.size1() != local_size)
        rMassMatrix.resize(local_size, local_size, false);
    noalias(rMassMatrix) = ZeroMatrix(local_size, local_size);
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    rMassMatrix(0,0) = one_six; rMassMatrix(0,1) = one_twelve; rMassMatrix(0,2) = one_twelve;
    rMassMatrix(1,0) = one_twelve; rMassMatrix(1,1) = one_six; rMassMatrix(1,2) = one_twelve;
    rMassMatrix(2,0) = one_twelve; rMassMatrix(2,1) = one_twelve; rMassMatrix(2,2) = one_six;
    // Assumption all the Gauss points have the same weight, so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int local_size = 4;
    // Resize, intialize and fill the mass matrix for linear tetrahedral elements
    if (rMassMatrix.size1() != local_size)
        rMassMatrix.resize(local_size, local_size, false);
    noalias(rMassMatrix) = ZeroMatrix(local_size, local_size);
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    rMassMatrix(0,0) = one_ten; rMassMatrix(0,1) = one_twenty; rMassMatrix(0,2) = one_twenty; rMassMatrix(0,3) = one_twenty;
    rMassMatrix(1,0) = one_twenty; rMassMatrix(1,1) = one_ten; rMassMatrix(1,2) = one_twenty; rMassMatrix(1,3) = one_twenty;
    rMassMatrix(2,0) = one_twenty; rMassMatrix(2,1) = one_twenty; rMassMatrix(2,2) = one_ten; rMassMatrix(2,3) = one_twenty;
    rMassMatrix(3,0) = one_twenty; rMassMatrix(3,1) = one_twenty; rMassMatrix(3,2) = one_twenty; rMassMatrix(3,3) = one_ten;
    // Assumption all the Gauss points have the same weight, so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    auto& r_geometry = GetGeometry();
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
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateOrthogonalSubgridScaleSystem(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
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
int SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
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

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::InitializeEulerianElement(
    ElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    // Getting data for the given geometry and integration method
    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    array_1d<double,TNumNodes> N_aux;
    GeometryUtils::CalculateGeometryData(r_geometry,rVariables.DN_DX,N_aux,rVariables.volume);
    rVariables.N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize some scalar variables
    rVariables.lumping_factor = 1.00 / double(TNumNodes);
    rVariables.diffusivity = 0.0;
    rVariables.dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

    for(unsigned int node_element = 0; node_element<local_size; node_element++)
{
    // Observations
    // * SGS time derivative term approximated as (phi-phi_old)/(RK_time_coefficient*delta_time)
    //   observe that for RK step = 1 ASGS time derivative term = 0 because phi = phi_old (null acceleration wrt step n)
    // * convective velocity and forcing term:
    //   RK step 1: evaluated at previous time step
    //   RK steps 2 and 3: linear interpolation between current and oldprevious time step
    //   RK step 4: evaluated at current time step
    // * convective_velocity = velocity - velocity_mesh
    //   velocity_mesh = 0 in eulerian framework
    if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==1)
    {
        rVariables.RK_time_coefficient = std::numeric_limits<double>::max();
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable(),1);
        rVariables.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[0];
        rVariables.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[1];
        rVariables.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[2];
    }
    else if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==2 || r_process_info.GetValue(RUNGE_KUTTA_STEP)==3)
    {
        rVariables.RK_time_coefficient = 0.5;
        rVariables.forcing[node_element] = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable()) + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable(),1));
        rVariables.convective_velocity(node_element,0) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[0] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[0]);
        rVariables.convective_velocity(node_element,1) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[1] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[1]);
        rVariables.convective_velocity(node_element,2) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[2] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[2]);
    }
    else
    {
        rVariables.RK_time_coefficient = 1.0;
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable());
        rVariables.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[0];
        rVariables.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[1];
        rVariables.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[2];
    }
    rVariables.oss_projection[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetProjectionVariable());
    rVariables.diffusivity += r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetDiffusionVariable());
    rVariables.delta_time = r_process_info[DELTA_TIME];
    rVariables.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable());
    rVariables.unknown_old[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable(),1);
}
    // divide by number of nodes scalar variables
    rVariables.diffusivity *= rVariables.lumping_factor;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<2>::CalculateLocalSystemInternal(
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

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             0.166666666666667*v(0,0);
const double clhs3 =             0.166666666666667*v(2,0);
const double clhs4 =             clhs2 + clhs3 + 0.666666666666667*v(1,0);
const double clhs5 =             DN_DX_0_0*clhs4;
const double clhs6 =             0.166666666666667*v(0,1);
const double clhs7 =             0.166666666666667*v(2,1);
const double clhs8 =             clhs6 + clhs7 + 0.666666666666667*v(1,1);
const double clhs9 =             DN_DX_0_1*clhs8;
const double clhs10 =             clhs5 + clhs9;
const double clhs11 =             clhs0*clhs1*clhs10*tau[1];
const double clhs12 =             0.166666666666667*clhs11;
const double clhs13 =             DN_DX_0_0*v(0,0);
const double clhs14 =             DN_DX_0_1*v(0,1);
const double clhs15 =             DN_DX_1_0*v(1,0);
const double clhs16 =             DN_DX_1_1*v(1,1);
const double clhs17 =             DN_DX_2_0*v(2,0);
const double clhs18 =             DN_DX_2_1*v(2,1);
const double clhs19 =             clhs13 + clhs14 + clhs15 + clhs16 + clhs17 + clhs18;
const double clhs20 =             clhs10*clhs19*tau[1];
const double clhs21 =             0.166666666666667*clhs20;
const double clhs22 =             0.5*clhs13;
const double clhs23 =             0.5*clhs14;
const double clhs24 =             0.5*clhs15;
const double clhs25 =             0.5*clhs16;
const double clhs26 =             0.5*clhs17;
const double clhs27 =             0.5*clhs18;
const double clhs28 =             3*k;
const double clhs29 =             0.166666666666667*v(1,0);
const double clhs30 =             clhs2 + clhs29 + 0.666666666666667*v(2,0);
const double clhs31 =             DN_DX_0_0*clhs30;
const double clhs32 =             0.166666666666667*clhs31;
const double clhs33 =             0.166666666666667*clhs5;
const double clhs34 =             clhs29 + clhs3 + 0.666666666666667*v(0,0);
const double clhs35 =             DN_DX_0_0*clhs34;
const double clhs36 =             0.166666666666667*v(1,1);
const double clhs37 =             clhs36 + clhs6 + 0.666666666666667*v(2,1);
const double clhs38 =             DN_DX_0_1*clhs37;
const double clhs39 =             0.166666666666667*clhs38;
const double clhs40 =             0.166666666666667*clhs9;
const double clhs41 =             clhs36 + clhs7 + 0.666666666666667*v(0,1);
const double clhs42 =             DN_DX_0_1*clhs41;
const double clhs43 =             clhs35 + clhs42;
const double clhs44 =             clhs31 + clhs38;
const double clhs45 =             clhs0*clhs1*clhs43*tau[0];
const double clhs46 =             clhs0*clhs1*clhs44*tau[2];
const double clhs47 =             0.166666666666667*clhs46;
const double clhs48 =             clhs19*clhs43*tau[0];
const double clhs49 =             clhs19*clhs44*tau[2];
const double clhs50 =             0.166666666666667*clhs49;
const double clhs51 =             0.166666666666667*clhs45;
const double clhs52 =             0.166666666666667*clhs48;
const double clhs53 =             0.25*clhs13;
const double clhs54 =             0.25*clhs14;
const double clhs55 =             0.25*clhs15;
const double clhs56 =             0.25*clhs16;
const double clhs57 =             0.25*clhs17;
const double clhs58 =             0.25*clhs18;
const double clhs59 =             3*DN_DX_0_0*k;
const double clhs60 =             DN_DX_1_0*clhs59;
const double clhs61 =             3*DN_DX_0_1*k;
const double clhs62 =             DN_DX_1_1*clhs61;
const double clhs63 =             DN_DX_1_0*clhs30;
const double clhs64 =             0.166666666666667*clhs63;
const double clhs65 =             DN_DX_1_0*clhs4;
const double clhs66 =             0.166666666666667*clhs65;
const double clhs67 =             DN_DX_1_0*clhs34;
const double clhs68 =             DN_DX_1_1*clhs37;
const double clhs69 =             0.166666666666667*clhs68;
const double clhs70 =             DN_DX_1_1*clhs8;
const double clhs71 =             0.166666666666667*clhs70;
const double clhs72 =             DN_DX_1_1*clhs41;
const double clhs73 =             clhs67 + clhs72;
const double clhs74 =             clhs43*tau[0];
const double clhs75 =             clhs73*clhs74;
const double clhs76 =             clhs65 + clhs70;
const double clhs77 =             clhs10*tau[1];
const double clhs78 =             clhs76*clhs77;
const double clhs79 =             clhs63 + clhs68;
const double clhs80 =             clhs44*tau[2];
const double clhs81 =             clhs79*clhs80;
const double clhs82 =             DN_DX_2_0*clhs34;
const double clhs83 =             DN_DX_2_1*clhs41;
const double clhs84 =             clhs82 + clhs83;
const double clhs85 =             DN_DX_2_0*clhs4;
const double clhs86 =             DN_DX_2_1*clhs8;
const double clhs87 =             clhs85 + clhs86;
const double clhs88 =             DN_DX_2_0*clhs30;
const double clhs89 =             DN_DX_2_1*clhs37;
const double clhs90 =             clhs88 + clhs89;
const double clhs91 =             DN_DX_2_0*clhs59 + DN_DX_2_1*clhs61 + clhs53 + clhs54 + clhs55 + clhs56 + clhs57 + clhs58 + clhs74*clhs84 + clhs77*clhs87 + clhs80*clhs90;
const double clhs92 =             0.166666666666667*clhs88;
const double clhs93 =             0.166666666666667*clhs85;
const double clhs94 =             0.166666666666667*clhs89;
const double clhs95 =             0.166666666666667*clhs86;
const double clhs96 =             0.166666666666667*clhs35;
const double clhs97 =             0.166666666666667*clhs42;
const double clhs98 =             clhs0*clhs1*clhs73*tau[0];
const double clhs99 =             clhs0*clhs1*clhs76*tau[1];
const double clhs100 =             0.166666666666667*clhs99;
const double clhs101 =             clhs0*clhs1*clhs79*tau[2];
const double clhs102 =             0.166666666666667*clhs101;
const double clhs103 =             clhs19*clhs73*tau[0];
const double clhs104 =             clhs19*clhs76*tau[1];
const double clhs105 =             0.166666666666667*clhs104;
const double clhs106 =             clhs19*clhs79*tau[2];
const double clhs107 =             0.166666666666667*clhs106;
const double clhs108 =             0.166666666666667*clhs67;
const double clhs109 =             0.166666666666667*clhs72;
const double clhs110 =             0.166666666666667*clhs98;
const double clhs111 =             0.166666666666667*clhs103;
const double clhs112 =             DN_DX_1_0*DN_DX_2_0*clhs28;
const double clhs113 =             DN_DX_1_1*DN_DX_2_1*clhs28;
const double clhs114 =             0.166666666666667*clhs82;
const double clhs115 =             0.166666666666667*clhs83;
const double clhs116 =             clhs73*clhs84*tau[0];
const double clhs117 =             clhs76*clhs87*tau[1];
const double clhs118 =             clhs79*clhs90*tau[2];
const double clhs119 =             clhs0*clhs1*clhs84*tau[0];
const double clhs120 =             clhs0*clhs1*clhs87*tau[1];
const double clhs121 =             0.166666666666667*clhs120;
const double clhs122 =             clhs0*clhs1*clhs90*tau[2];
const double clhs123 =             0.166666666666667*clhs122;
const double clhs124 =             clhs19*clhs84*tau[0];
const double clhs125 =             clhs19*clhs87*tau[1];
const double clhs126 =             0.166666666666667*clhs125;
const double clhs127 =             clhs19*clhs90*tau[2];
const double clhs128 =             0.166666666666667*clhs127;
const double clhs129 =             0.166666666666667*clhs119;
const double clhs130 =             0.166666666666667*clhs124;
            lhs(0,0)=pow(DN_DX_0_0, 2)*clhs28 + pow(DN_DX_0_1, 2)*clhs28 + pow(clhs10, 2)*tau[1] + clhs12 + clhs21 + clhs22 + clhs23 + clhs24 + clhs25 + clhs26 + clhs27 + clhs32 + clhs33 + 0.666666666666667*clhs35 + clhs39 + clhs40 + 0.666666666666667*clhs42 + pow(clhs43, 2)*tau[0] + pow(clhs44, 2)*tau[2] + 0.666666666666667*clhs45 + clhs47 + 0.666666666666667*clhs48 + clhs50;
            lhs(0,1)=0.666666666666667*clhs11 + 0.666666666666667*clhs20 + clhs47 + clhs50 + clhs51 + clhs52 + clhs53 + clhs54 + clhs55 + clhs56 + clhs57 + clhs58 + clhs60 + clhs62 + clhs64 + clhs66 + 0.666666666666667*clhs67 + clhs69 + clhs71 + 0.666666666666667*clhs72 + clhs75 + clhs78 + clhs81;
            lhs(0,2)=clhs12 + clhs21 + 0.666666666666667*clhs46 + 0.666666666666667*clhs49 + clhs51 + clhs52 + 0.666666666666667*clhs82 + 0.666666666666667*clhs83 + clhs91 + clhs92 + clhs93 + clhs94 + clhs95;
            lhs(1,0)=clhs100 + clhs102 + 0.666666666666667*clhs103 + clhs105 + clhs107 + clhs32 + clhs39 + 0.666666666666667*clhs5 + clhs53 + clhs54 + clhs55 + clhs56 + clhs57 + clhs58 + clhs60 + clhs62 + clhs75 + clhs78 + clhs81 + 0.666666666666667*clhs9 + clhs96 + clhs97 + 0.666666666666667*clhs98;
            lhs(1,1)=pow(DN_DX_1_0, 2)*clhs28 + pow(DN_DX_1_1, 2)*clhs28 + clhs102 + 0.666666666666667*clhs104 + clhs107 + clhs108 + clhs109 + clhs110 + clhs111 + clhs22 + clhs23 + clhs24 + clhs25 + clhs26 + clhs27 + clhs64 + 0.666666666666667*clhs65 + clhs69 + 0.666666666666667*clhs70 + pow(clhs73, 2)*tau[0] + pow(clhs76, 2)*tau[1] + pow(clhs79, 2)*tau[2] + 0.666666666666667*clhs99;
            lhs(1,2)=clhs100 + 0.666666666666667*clhs101 + clhs105 + 0.666666666666667*clhs106 + clhs110 + clhs111 + clhs112 + clhs113 + clhs114 + clhs115 + clhs116 + clhs117 + clhs118 + clhs53 + clhs54 + clhs55 + clhs56 + clhs57 + clhs58 + 0.666666666666667*clhs85 + 0.666666666666667*clhs86 + clhs92 + clhs94;
            lhs(2,0)=0.666666666666667*clhs119 + clhs121 + clhs123 + 0.666666666666667*clhs124 + clhs126 + clhs128 + 0.666666666666667*clhs31 + clhs33 + 0.666666666666667*clhs38 + clhs40 + clhs91 + clhs96 + clhs97;
            lhs(2,1)=clhs108 + clhs109 + clhs112 + clhs113 + clhs116 + clhs117 + clhs118 + 0.666666666666667*clhs120 + clhs123 + 0.666666666666667*clhs125 + clhs128 + clhs129 + clhs130 + clhs53 + clhs54 + clhs55 + clhs56 + clhs57 + clhs58 + 0.666666666666667*clhs63 + clhs66 + 0.666666666666667*clhs68 + clhs71;
            lhs(2,2)=pow(DN_DX_2_0, 2)*clhs28 + pow(DN_DX_2_1, 2)*clhs28 + clhs114 + clhs115 + clhs121 + 0.666666666666667*clhs122 + clhs126 + 0.666666666666667*clhs127 + clhs129 + clhs130 + clhs22 + clhs23 + clhs24 + clhs25 + clhs26 + clhs27 + pow(clhs84, 2)*tau[0] + pow(clhs87, 2)*tau[1] + 0.666666666666667*clhs88 + 0.666666666666667*clhs89 + pow(clhs90, 2)*tau[2] + clhs93 + clhs95;


    const double crhs0 =             0.25*f[1];
const double crhs1 =             0.166666666666667*v(0,0);
const double crhs2 =             0.166666666666667*v(2,0);
const double crhs3 =             crhs1 + crhs2 + 0.666666666666667*v(1,0);
const double crhs4 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs5 =             crhs3*crhs4;
const double crhs6 =             -0.166666666666667*crhs5;
const double crhs7 =             0.166666666666667*v(0,1);
const double crhs8 =             0.166666666666667*v(2,1);
const double crhs9 =             crhs7 + crhs8 + 0.666666666666667*v(1,1);
const double crhs10 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs11 =             crhs10*crhs9;
const double crhs12 =             -0.166666666666667*crhs11;
const double crhs13 =             0.166666666666667*phi[0];
const double crhs14 =             0.166666666666667*phi[2];
const double crhs15 =             crhs13 + crhs14 + 0.666666666666667*phi[1];
const double crhs16 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs17 =             crhs15*crhs16;
const double crhs18 =             -0.166666666666667*crhs17;
const double crhs19 =             0.25*f[2];
const double crhs20 =             3*crhs4*k;
const double crhs21 =             3*crhs10*k;
const double crhs22 =             0.166666666666667*v(1,0);
const double crhs23 =             crhs1 + crhs22 + 0.666666666666667*v(2,0);
const double crhs24 =             crhs23*crhs4;
const double crhs25 =             -0.166666666666667*crhs24;
const double crhs26 =             crhs2 + crhs22 + 0.666666666666667*v(0,0);
const double crhs27 =             crhs26*crhs4;
const double crhs28 =             0.166666666666667*v(1,1);
const double crhs29 =             crhs28 + crhs7 + 0.666666666666667*v(2,1);
const double crhs30 =             crhs10*crhs29;
const double crhs31 =             -0.166666666666667*crhs30;
const double crhs32 =             crhs28 + crhs8 + 0.666666666666667*v(0,1);
const double crhs33 =             crhs10*crhs32;
const double crhs34 =             0.166666666666667*phi[1];
const double crhs35 =             crhs13 + crhs34 + 0.666666666666667*phi[2];
const double crhs36 =             crhs16*crhs35;
const double crhs37 =             -0.166666666666667*crhs36;
const double crhs38 =             crhs14 + crhs34 + 0.666666666666667*phi[0];
const double crhs39 =             crhs16*crhs38;
const double crhs40 =             0.166666666666667*f[1];
const double crhs41 =             0.166666666666667*f[2];
const double crhs42 =             crhs40 + crhs41 + 0.666666666666667*f[0];
const double crhs43 =             tau[0]*(DN_DX_0_0*crhs26 + DN_DX_0_1*crhs32);
const double crhs44 =             0.166666666666667*prj[1];
const double crhs45 =             0.166666666666667*prj[2];
const double crhs46 =             crhs44 + crhs45 + 0.666666666666667*prj[0];
const double crhs47 =             0.166666666666667*f[0];
const double crhs48 =             crhs41 + crhs47 + 0.666666666666667*f[1];
const double crhs49 =             tau[1]*(DN_DX_0_0*crhs3 + DN_DX_0_1*crhs9);
const double crhs50 =             0.166666666666667*prj[0];
const double crhs51 =             crhs45 + crhs50 + 0.666666666666667*prj[1];
const double crhs52 =             crhs40 + crhs47 + 0.666666666666667*f[2];
const double crhs53 =             tau[2]*(DN_DX_0_0*crhs23 + DN_DX_0_1*crhs29);
const double crhs54 =             crhs44 + crhs50 + 0.666666666666667*prj[2];
const double crhs55 =             1.0/RK_time_coefficient;
const double crhs56 =             1.0/delta_time;
const double crhs57 =             -0.166666666666667*phi_old[1];
const double crhs58 =             -0.166666666666667*phi_old[2];
const double crhs59 =             crhs55*crhs56*(crhs38 + crhs57 + crhs58 - 0.666666666666667*phi_old[0]);
const double crhs60 =             -0.166666666666667*phi_old[0];
const double crhs61 =             crhs55*crhs56*(crhs15 + crhs58 + crhs60 - 0.666666666666667*phi_old[1]);
const double crhs62 =             crhs55*crhs56*(crhs35 + crhs57 + crhs60 - 0.666666666666667*phi_old[2]);
const double crhs63 =             crhs27 + crhs33;
const double crhs64 =             crhs11 + crhs5;
const double crhs65 =             crhs24 + crhs30;
const double crhs66 =             -0.166666666666667*crhs27 - 0.166666666666667*crhs33 - 0.166666666666667*crhs39 + 0.25*f[0];
const double crhs67 =             tau[0]*(DN_DX_1_0*crhs26 + DN_DX_1_1*crhs32);
const double crhs68 =             tau[1]*(DN_DX_1_0*crhs3 + DN_DX_1_1*crhs9);
const double crhs69 =             tau[2]*(DN_DX_1_0*crhs23 + DN_DX_1_1*crhs29);
const double crhs70 =             tau[0]*(DN_DX_2_0*crhs26 + DN_DX_2_1*crhs32);
const double crhs71 =             tau[1]*(DN_DX_2_0*crhs3 + DN_DX_2_1*crhs9);
const double crhs72 =             tau[2]*(DN_DX_2_0*crhs23 + DN_DX_2_1*crhs29);
            rhs[0]=-DN_DX_0_0*crhs20 - DN_DX_0_1*crhs21 + crhs0 + crhs12 - crhs17*crhs49 + crhs18 + crhs19 + crhs25 - 0.666666666666667*crhs27 + crhs31 - 0.666666666666667*crhs33 - crhs36*crhs53 + crhs37 - crhs39*crhs43 - 0.666666666666667*crhs39 + crhs42*crhs43 + crhs43*crhs46 - crhs43*crhs59 - crhs43*crhs63 + crhs48*crhs49 + crhs49*crhs51 - crhs49*crhs61 - crhs49*crhs64 + crhs52*crhs53 + crhs53*crhs54 - crhs53*crhs62 - crhs53*crhs65 + crhs6 + 0.5*f[0];
            rhs[1]=-DN_DX_1_0*crhs20 - DN_DX_1_1*crhs21 - 0.666666666666667*crhs11 - crhs17*crhs68 - 0.666666666666667*crhs17 + crhs19 + crhs25 + crhs31 - crhs36*crhs69 + crhs37 - crhs39*crhs67 + crhs42*crhs67 + crhs46*crhs67 + crhs48*crhs68 - 0.666666666666667*crhs5 + crhs51*crhs68 + crhs52*crhs69 + crhs54*crhs69 - crhs59*crhs67 - crhs61*crhs68 - crhs62*crhs69 - crhs63*crhs67 - crhs64*crhs68 - crhs65*crhs69 + crhs66 + 0.5*f[1];
            rhs[2]=-DN_DX_2_0*crhs20 - DN_DX_2_1*crhs21 + crhs0 + crhs12 - crhs17*crhs71 + crhs18 - 0.666666666666667*crhs24 - 0.666666666666667*crhs30 - crhs36*crhs72 - 0.666666666666667*crhs36 - crhs39*crhs70 + crhs42*crhs70 + crhs46*crhs70 + crhs48*crhs71 + crhs51*crhs71 + crhs52*crhs72 + crhs54*crhs72 - crhs59*crhs70 + crhs6 - crhs61*crhs71 - crhs62*crhs72 - crhs63*crhs70 - crhs64*crhs71 - crhs65*crhs72 + crhs66 + 0.5*f[2];


    const double local_size = 3;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3>::CalculateLocalSystemInternal(
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

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             0.1381966*v(0,0);
const double clhs3 =             0.1381966*v(2,0);
const double clhs4 =             0.1381966*v(3,0);
const double clhs5 =             clhs3 + clhs4;
const double clhs6 =             clhs2 + clhs5 + 0.5854102*v(1,0);
const double clhs7 =             DN_DX_0_0*clhs6;
const double clhs8 =             0.1381966*v(0,1);
const double clhs9 =             0.1381966*v(2,1);
const double clhs10 =             0.1381966*v(3,1);
const double clhs11 =             clhs10 + clhs9;
const double clhs12 =             clhs11 + clhs8 + 0.5854102*v(1,1);
const double clhs13 =             DN_DX_0_1*clhs12;
const double clhs14 =             0.1381966*v(0,2);
const double clhs15 =             0.1381966*v(2,2);
const double clhs16 =             0.1381966*v(3,2);
const double clhs17 =             clhs15 + clhs16;
const double clhs18 =             clhs14 + clhs17 + 0.5854102*v(1,2);
const double clhs19 =             DN_DX_0_2*clhs18;
const double clhs20 =             clhs13 + clhs19 + clhs7;
const double clhs21 =             clhs0*clhs1*clhs20*tau[1];
const double clhs22 =             0.1381966*clhs21;
const double clhs23 =             DN_DX_0_0*v(0,0);
const double clhs24 =             DN_DX_0_1*v(0,1);
const double clhs25 =             DN_DX_0_2*v(0,2);
const double clhs26 =             DN_DX_1_0*v(1,0);
const double clhs27 =             DN_DX_1_1*v(1,1);
const double clhs28 =             DN_DX_1_2*v(1,2);
const double clhs29 =             DN_DX_2_0*v(2,0);
const double clhs30 =             DN_DX_2_1*v(2,1);
const double clhs31 =             DN_DX_2_2*v(2,2);
const double clhs32 =             DN_DX_3_0*v(3,0);
const double clhs33 =             DN_DX_3_1*v(3,1);
const double clhs34 =             DN_DX_3_2*v(3,2);
const double clhs35 =             clhs23 + clhs24 + clhs25 + clhs26 + clhs27 + clhs28 + clhs29 + clhs30 + clhs31 + clhs32 + clhs33 + clhs34;
const double clhs36 =             clhs20*clhs35*tau[1];
const double clhs37 =             0.1381966*clhs36;
const double clhs38 =             0.40000000301872*clhs23;
const double clhs39 =             0.40000000301872*clhs24;
const double clhs40 =             0.40000000301872*clhs25;
const double clhs41 =             0.40000000301872*clhs26;
const double clhs42 =             0.40000000301872*clhs27;
const double clhs43 =             0.40000000301872*clhs28;
const double clhs44 =             0.40000000301872*clhs29;
const double clhs45 =             0.40000000301872*clhs30;
const double clhs46 =             0.40000000301872*clhs31;
const double clhs47 =             0.40000000301872*clhs32;
const double clhs48 =             0.40000000301872*clhs33;
const double clhs49 =             0.40000000301872*clhs34;
const double clhs50 =             4*k;
const double clhs51 =             0.1381966*v(1,0);
const double clhs52 =             clhs2 + clhs51;
const double clhs53 =             clhs3 + clhs52 + 0.5854102*v(3,0);
const double clhs54 =             DN_DX_0_0*clhs53;
const double clhs55 =             0.1381966*clhs54;
const double clhs56 =             clhs4 + clhs52 + 0.5854102*v(2,0);
const double clhs57 =             DN_DX_0_0*clhs56;
const double clhs58 =             0.1381966*clhs57;
const double clhs59 =             0.1381966*clhs7;
const double clhs60 =             clhs5 + clhs51 + 0.5854102*v(0,0);
const double clhs61 =             DN_DX_0_0*clhs60;
const double clhs62 =             0.1381966*v(1,1);
const double clhs63 =             clhs62 + clhs8;
const double clhs64 =             clhs63 + clhs9 + 0.5854102*v(3,1);
const double clhs65 =             DN_DX_0_1*clhs64;
const double clhs66 =             0.1381966*clhs65;
const double clhs67 =             clhs10 + clhs63 + 0.5854102*v(2,1);
const double clhs68 =             DN_DX_0_1*clhs67;
const double clhs69 =             0.1381966*clhs68;
const double clhs70 =             0.1381966*clhs13;
const double clhs71 =             clhs11 + clhs62 + 0.5854102*v(0,1);
const double clhs72 =             DN_DX_0_1*clhs71;
const double clhs73 =             0.1381966*v(1,2);
const double clhs74 =             clhs14 + clhs73;
const double clhs75 =             clhs15 + clhs74 + 0.5854102*v(3,2);
const double clhs76 =             DN_DX_0_2*clhs75;
const double clhs77 =             0.1381966*clhs76;
const double clhs78 =             clhs16 + clhs74 + 0.5854102*v(2,2);
const double clhs79 =             DN_DX_0_2*clhs78;
const double clhs80 =             0.1381966*clhs79;
const double clhs81 =             0.1381966*clhs19;
const double clhs82 =             clhs17 + clhs73 + 0.5854102*v(0,2);
const double clhs83 =             DN_DX_0_2*clhs82;
const double clhs84 =             clhs61 + clhs72 + clhs83;
const double clhs85 =             clhs57 + clhs68 + clhs79;
const double clhs86 =             clhs54 + clhs65 + clhs76;
const double clhs87 =             clhs0*clhs1*clhs84*tau[0];
const double clhs88 =             clhs0*clhs1*clhs85*tau[2];
const double clhs89 =             0.1381966*clhs88;
const double clhs90 =             clhs0*clhs1*clhs86*tau[3];
const double clhs91 =             0.1381966*clhs90;
const double clhs92 =             clhs35*clhs84*tau[0];
const double clhs93 =             clhs35*clhs85*tau[2];
const double clhs94 =             0.1381966*clhs93;
const double clhs95 =             clhs35*clhs86*tau[3];
const double clhs96 =             0.1381966*clhs95;
const double clhs97 =             0.1381966*clhs87;
const double clhs98 =             0.1381966*clhs92;
const double clhs99 =             0.19999999899376*clhs23;
const double clhs100 =             0.19999999899376*clhs24;
const double clhs101 =             0.19999999899376*clhs25;
const double clhs102 =             0.19999999899376*clhs26;
const double clhs103 =             0.19999999899376*clhs27;
const double clhs104 =             0.19999999899376*clhs28;
const double clhs105 =             0.19999999899376*clhs29;
const double clhs106 =             0.19999999899376*clhs30;
const double clhs107 =             0.19999999899376*clhs31;
const double clhs108 =             0.19999999899376*clhs32;
const double clhs109 =             0.19999999899376*clhs33;
const double clhs110 =             0.19999999899376*clhs34;
const double clhs111 =             4*DN_DX_0_0*k;
const double clhs112 =             DN_DX_1_0*clhs111;
const double clhs113 =             4*DN_DX_0_1*k;
const double clhs114 =             DN_DX_1_1*clhs113;
const double clhs115 =             4*DN_DX_0_2*k;
const double clhs116 =             DN_DX_1_2*clhs115;
const double clhs117 =             DN_DX_1_0*clhs53;
const double clhs118 =             0.1381966*clhs117;
const double clhs119 =             DN_DX_1_0*clhs56;
const double clhs120 =             0.1381966*clhs119;
const double clhs121 =             DN_DX_1_0*clhs6;
const double clhs122 =             0.1381966*clhs121;
const double clhs123 =             DN_DX_1_0*clhs60;
const double clhs124 =             DN_DX_1_1*clhs64;
const double clhs125 =             0.1381966*clhs124;
const double clhs126 =             DN_DX_1_1*clhs67;
const double clhs127 =             0.1381966*clhs126;
const double clhs128 =             DN_DX_1_1*clhs12;
const double clhs129 =             0.1381966*clhs128;
const double clhs130 =             DN_DX_1_1*clhs71;
const double clhs131 =             DN_DX_1_2*clhs75;
const double clhs132 =             0.1381966*clhs131;
const double clhs133 =             DN_DX_1_2*clhs78;
const double clhs134 =             0.1381966*clhs133;
const double clhs135 =             DN_DX_1_2*clhs18;
const double clhs136 =             0.1381966*clhs135;
const double clhs137 =             DN_DX_1_2*clhs82;
const double clhs138 =             clhs123 + clhs130 + clhs137;
const double clhs139 =             clhs84*tau[0];
const double clhs140 =             clhs138*clhs139;
const double clhs141 =             clhs121 + clhs128 + clhs135;
const double clhs142 =             clhs20*tau[1];
const double clhs143 =             clhs141*clhs142;
const double clhs144 =             clhs119 + clhs126 + clhs133;
const double clhs145 =             clhs85*tau[2];
const double clhs146 =             clhs144*clhs145;
const double clhs147 =             clhs117 + clhs124 + clhs131;
const double clhs148 =             clhs86*tau[3];
const double clhs149 =             clhs147*clhs148;
const double clhs150 =             DN_DX_2_0*clhs111;
const double clhs151 =             DN_DX_2_1*clhs113;
const double clhs152 =             DN_DX_2_2*clhs115;
const double clhs153 =             DN_DX_2_0*clhs53;
const double clhs154 =             0.1381966*clhs153;
const double clhs155 =             DN_DX_2_0*clhs56;
const double clhs156 =             0.1381966*clhs155;
const double clhs157 =             DN_DX_2_0*clhs6;
const double clhs158 =             0.1381966*clhs157;
const double clhs159 =             DN_DX_2_0*clhs60;
const double clhs160 =             DN_DX_2_1*clhs64;
const double clhs161 =             0.1381966*clhs160;
const double clhs162 =             DN_DX_2_1*clhs67;
const double clhs163 =             0.1381966*clhs162;
const double clhs164 =             DN_DX_2_1*clhs12;
const double clhs165 =             0.1381966*clhs164;
const double clhs166 =             DN_DX_2_1*clhs71;
const double clhs167 =             DN_DX_2_2*clhs75;
const double clhs168 =             0.1381966*clhs167;
const double clhs169 =             DN_DX_2_2*clhs78;
const double clhs170 =             0.1381966*clhs169;
const double clhs171 =             DN_DX_2_2*clhs18;
const double clhs172 =             0.1381966*clhs171;
const double clhs173 =             DN_DX_2_2*clhs82;
const double clhs174 =             clhs159 + clhs166 + clhs173;
const double clhs175 =             clhs139*clhs174;
const double clhs176 =             clhs157 + clhs164 + clhs171;
const double clhs177 =             clhs142*clhs176;
const double clhs178 =             clhs155 + clhs162 + clhs169;
const double clhs179 =             clhs145*clhs178;
const double clhs180 =             clhs153 + clhs160 + clhs167;
const double clhs181 =             clhs148*clhs180;
const double clhs182 =             DN_DX_3_0*clhs60;
const double clhs183 =             DN_DX_3_1*clhs71;
const double clhs184 =             DN_DX_3_2*clhs82;
const double clhs185 =             clhs182 + clhs183 + clhs184;
const double clhs186 =             DN_DX_3_0*clhs6;
const double clhs187 =             DN_DX_3_1*clhs12;
const double clhs188 =             DN_DX_3_2*clhs18;
const double clhs189 =             clhs186 + clhs187 + clhs188;
const double clhs190 =             DN_DX_3_0*clhs56;
const double clhs191 =             DN_DX_3_1*clhs67;
const double clhs192 =             DN_DX_3_2*clhs78;
const double clhs193 =             clhs190 + clhs191 + clhs192;
const double clhs194 =             DN_DX_3_0*clhs53;
const double clhs195 =             DN_DX_3_1*clhs64;
const double clhs196 =             DN_DX_3_2*clhs75;
const double clhs197 =             clhs194 + clhs195 + clhs196;
const double clhs198 =             DN_DX_3_0*clhs111 + DN_DX_3_1*clhs113 + DN_DX_3_2*clhs115 + clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs139*clhs185 + clhs142*clhs189 + clhs145*clhs193 + clhs148*clhs197 + clhs99;
const double clhs199 =             0.1381966*clhs194;
const double clhs200 =             0.1381966*clhs190;
const double clhs201 =             0.1381966*clhs186;
const double clhs202 =             0.1381966*clhs195;
const double clhs203 =             0.1381966*clhs191;
const double clhs204 =             0.1381966*clhs187;
const double clhs205 =             0.1381966*clhs196;
const double clhs206 =             0.1381966*clhs192;
const double clhs207 =             0.1381966*clhs188;
const double clhs208 =             0.1381966*clhs61;
const double clhs209 =             0.1381966*clhs72;
const double clhs210 =             0.1381966*clhs83;
const double clhs211 =             clhs0*clhs1*clhs138*tau[0];
const double clhs212 =             clhs0*clhs1*clhs141*tau[1];
const double clhs213 =             0.1381966*clhs212;
const double clhs214 =             clhs0*clhs1*clhs144*tau[2];
const double clhs215 =             0.1381966*clhs214;
const double clhs216 =             clhs0*clhs1*clhs147*tau[3];
const double clhs217 =             0.1381966*clhs216;
const double clhs218 =             clhs138*clhs35*tau[0];
const double clhs219 =             clhs141*clhs35*tau[1];
const double clhs220 =             0.1381966*clhs219;
const double clhs221 =             clhs144*clhs35*tau[2];
const double clhs222 =             0.1381966*clhs221;
const double clhs223 =             clhs147*clhs35*tau[3];
const double clhs224 =             0.1381966*clhs223;
const double clhs225 =             0.1381966*clhs123;
const double clhs226 =             0.1381966*clhs130;
const double clhs227 =             0.1381966*clhs137;
const double clhs228 =             0.1381966*clhs211;
const double clhs229 =             0.1381966*clhs218;
const double clhs230 =             4*DN_DX_1_0*k;
const double clhs231 =             DN_DX_2_0*clhs230;
const double clhs232 =             4*DN_DX_1_1*k;
const double clhs233 =             DN_DX_2_1*clhs232;
const double clhs234 =             4*DN_DX_1_2*k;
const double clhs235 =             DN_DX_2_2*clhs234;
const double clhs236 =             0.1381966*clhs159;
const double clhs237 =             0.1381966*clhs166;
const double clhs238 =             0.1381966*clhs173;
const double clhs239 =             clhs138*tau[0];
const double clhs240 =             clhs174*clhs239;
const double clhs241 =             clhs141*tau[1];
const double clhs242 =             clhs176*clhs241;
const double clhs243 =             clhs144*tau[2];
const double clhs244 =             clhs178*clhs243;
const double clhs245 =             clhs147*tau[3];
const double clhs246 =             clhs180*clhs245;
const double clhs247 =             DN_DX_3_0*clhs230;
const double clhs248 =             DN_DX_3_1*clhs232;
const double clhs249 =             DN_DX_3_2*clhs234;
const double clhs250 =             0.1381966*clhs182;
const double clhs251 =             0.1381966*clhs183;
const double clhs252 =             0.1381966*clhs184;
const double clhs253 =             clhs185*clhs239;
const double clhs254 =             clhs189*clhs241;
const double clhs255 =             clhs193*clhs243;
const double clhs256 =             clhs197*clhs245;
const double clhs257 =             clhs0*clhs1*clhs174*tau[0];
const double clhs258 =             clhs0*clhs1*clhs176*tau[1];
const double clhs259 =             0.1381966*clhs258;
const double clhs260 =             clhs0*clhs1*clhs178*tau[2];
const double clhs261 =             0.1381966*clhs260;
const double clhs262 =             clhs0*clhs1*clhs180*tau[3];
const double clhs263 =             0.1381966*clhs262;
const double clhs264 =             clhs174*clhs35*tau[0];
const double clhs265 =             clhs176*clhs35*tau[1];
const double clhs266 =             0.1381966*clhs265;
const double clhs267 =             clhs178*clhs35*tau[2];
const double clhs268 =             0.1381966*clhs267;
const double clhs269 =             clhs180*clhs35*tau[3];
const double clhs270 =             0.1381966*clhs269;
const double clhs271 =             0.1381966*clhs257;
const double clhs272 =             0.1381966*clhs264;
const double clhs273 =             DN_DX_2_0*DN_DX_3_0*clhs50;
const double clhs274 =             DN_DX_2_1*DN_DX_3_1*clhs50;
const double clhs275 =             DN_DX_2_2*DN_DX_3_2*clhs50;
const double clhs276 =             clhs174*clhs185*tau[0];
const double clhs277 =             clhs176*clhs189*tau[1];
const double clhs278 =             clhs178*clhs193*tau[2];
const double clhs279 =             clhs180*clhs197*tau[3];
const double clhs280 =             clhs0*clhs1*clhs185*tau[0];
const double clhs281 =             clhs0*clhs1*clhs189*tau[1];
const double clhs282 =             0.1381966*clhs281;
const double clhs283 =             clhs0*clhs1*clhs193*tau[2];
const double clhs284 =             0.1381966*clhs283;
const double clhs285 =             clhs0*clhs1*clhs197*tau[3];
const double clhs286 =             0.1381966*clhs285;
const double clhs287 =             clhs185*clhs35*tau[0];
const double clhs288 =             clhs189*clhs35*tau[1];
const double clhs289 =             0.1381966*clhs288;
const double clhs290 =             clhs193*clhs35*tau[2];
const double clhs291 =             0.1381966*clhs290;
const double clhs292 =             clhs197*clhs35*tau[3];
const double clhs293 =             0.1381966*clhs292;
const double clhs294 =             0.1381966*clhs280;
const double clhs295 =             0.1381966*clhs287;
            lhs(0,0)=pow(DN_DX_0_0, 2)*clhs50 + pow(DN_DX_0_1, 2)*clhs50 + pow(DN_DX_0_2, 2)*clhs50 + pow(clhs20, 2)*tau[1] + clhs22 + clhs37 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49 + clhs55 + clhs58 + clhs59 + 0.5854102*clhs61 + clhs66 + clhs69 + clhs70 + 0.5854102*clhs72 + clhs77 + clhs80 + clhs81 + 0.5854102*clhs83 + pow(clhs84, 2)*tau[0] + pow(clhs85, 2)*tau[2] + pow(clhs86, 2)*tau[3] + 0.5854102*clhs87 + clhs89 + clhs91 + 0.5854102*clhs92 + clhs94 + clhs96;
            lhs(0,1)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs112 + clhs114 + clhs116 + clhs118 + clhs120 + clhs122 + 0.5854102*clhs123 + clhs125 + clhs127 + clhs129 + 0.5854102*clhs130 + clhs132 + clhs134 + clhs136 + 0.5854102*clhs137 + clhs140 + clhs143 + clhs146 + clhs149 + 0.5854102*clhs21 + 0.5854102*clhs36 + clhs89 + clhs91 + clhs94 + clhs96 + clhs97 + clhs98 + clhs99;
            lhs(0,2)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs150 + clhs151 + clhs152 + clhs154 + clhs156 + clhs158 + 0.5854102*clhs159 + clhs161 + clhs163 + clhs165 + 0.5854102*clhs166 + clhs168 + clhs170 + clhs172 + 0.5854102*clhs173 + clhs175 + clhs177 + clhs179 + clhs181 + clhs22 + clhs37 + 0.5854102*clhs88 + clhs91 + 0.5854102*clhs93 + clhs96 + clhs97 + clhs98 + clhs99;
            lhs(0,3)=0.5854102*clhs182 + 0.5854102*clhs183 + 0.5854102*clhs184 + clhs198 + clhs199 + clhs200 + clhs201 + clhs202 + clhs203 + clhs204 + clhs205 + clhs206 + clhs207 + clhs22 + clhs37 + clhs89 + 0.5854102*clhs90 + clhs94 + 0.5854102*clhs95 + clhs97 + clhs98;
            lhs(1,0)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs112 + clhs114 + clhs116 + 0.5854102*clhs13 + clhs140 + clhs143 + clhs146 + clhs149 + 0.5854102*clhs19 + clhs208 + clhs209 + clhs210 + 0.5854102*clhs211 + clhs213 + clhs215 + clhs217 + 0.5854102*clhs218 + clhs220 + clhs222 + clhs224 + clhs55 + clhs58 + clhs66 + clhs69 + 0.5854102*clhs7 + clhs77 + clhs80 + clhs99;
            lhs(1,1)=pow(DN_DX_1_0, 2)*clhs50 + pow(DN_DX_1_1, 2)*clhs50 + pow(DN_DX_1_2, 2)*clhs50 + clhs118 + clhs120 + 0.5854102*clhs121 + clhs125 + clhs127 + 0.5854102*clhs128 + clhs132 + clhs134 + 0.5854102*clhs135 + pow(clhs138, 2)*tau[0] + pow(clhs141, 2)*tau[1] + pow(clhs144, 2)*tau[2] + pow(clhs147, 2)*tau[3] + 0.5854102*clhs212 + clhs215 + clhs217 + 0.5854102*clhs219 + clhs222 + clhs224 + clhs225 + clhs226 + clhs227 + clhs228 + clhs229 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49;
            lhs(1,2)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs154 + clhs156 + 0.5854102*clhs157 + clhs161 + clhs163 + 0.5854102*clhs164 + clhs168 + clhs170 + 0.5854102*clhs171 + clhs213 + 0.5854102*clhs214 + clhs217 + clhs220 + 0.5854102*clhs221 + clhs224 + clhs228 + clhs229 + clhs231 + clhs233 + clhs235 + clhs236 + clhs237 + clhs238 + clhs240 + clhs242 + clhs244 + clhs246 + clhs99;
            lhs(1,3)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + 0.5854102*clhs186 + 0.5854102*clhs187 + 0.5854102*clhs188 + clhs199 + clhs200 + clhs202 + clhs203 + clhs205 + clhs206 + clhs213 + clhs215 + 0.5854102*clhs216 + clhs220 + clhs222 + 0.5854102*clhs223 + clhs228 + clhs229 + clhs247 + clhs248 + clhs249 + clhs250 + clhs251 + clhs252 + clhs253 + clhs254 + clhs255 + clhs256 + clhs99;
            lhs(2,0)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs150 + clhs151 + clhs152 + clhs175 + clhs177 + clhs179 + clhs181 + clhs208 + clhs209 + clhs210 + 0.5854102*clhs257 + clhs259 + clhs261 + clhs263 + 0.5854102*clhs264 + clhs266 + clhs268 + clhs270 + clhs55 + 0.5854102*clhs57 + clhs59 + clhs66 + 0.5854102*clhs68 + clhs70 + clhs77 + 0.5854102*clhs79 + clhs81 + clhs99;
            lhs(2,1)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + clhs118 + 0.5854102*clhs119 + clhs122 + clhs125 + 0.5854102*clhs126 + clhs129 + clhs132 + 0.5854102*clhs133 + clhs136 + clhs225 + clhs226 + clhs227 + clhs231 + clhs233 + clhs235 + clhs240 + clhs242 + clhs244 + clhs246 + 0.5854102*clhs258 + clhs261 + clhs263 + 0.5854102*clhs265 + clhs268 + clhs270 + clhs271 + clhs272 + clhs99;
            lhs(2,2)=pow(DN_DX_2_0, 2)*clhs50 + pow(DN_DX_2_1, 2)*clhs50 + pow(DN_DX_2_2, 2)*clhs50 + clhs154 + 0.5854102*clhs155 + clhs158 + clhs161 + 0.5854102*clhs162 + clhs165 + clhs168 + 0.5854102*clhs169 + clhs172 + pow(clhs174, 2)*tau[0] + pow(clhs176, 2)*tau[1] + pow(clhs178, 2)*tau[2] + pow(clhs180, 2)*tau[3] + clhs236 + clhs237 + clhs238 + clhs259 + 0.5854102*clhs260 + clhs263 + clhs266 + 0.5854102*clhs267 + clhs270 + clhs271 + clhs272 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49;
            lhs(2,3)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + 0.5854102*clhs190 + 0.5854102*clhs191 + 0.5854102*clhs192 + clhs199 + clhs201 + clhs202 + clhs204 + clhs205 + clhs207 + clhs250 + clhs251 + clhs252 + clhs259 + clhs261 + 0.5854102*clhs262 + clhs266 + clhs268 + 0.5854102*clhs269 + clhs271 + clhs272 + clhs273 + clhs274 + clhs275 + clhs276 + clhs277 + clhs278 + clhs279 + clhs99;
            lhs(3,0)=clhs198 + clhs208 + clhs209 + clhs210 + 0.5854102*clhs280 + clhs282 + clhs284 + clhs286 + 0.5854102*clhs287 + clhs289 + clhs291 + clhs293 + 0.5854102*clhs54 + clhs58 + clhs59 + 0.5854102*clhs65 + clhs69 + clhs70 + 0.5854102*clhs76 + clhs80 + clhs81;
            lhs(3,1)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + 0.5854102*clhs117 + clhs120 + clhs122 + 0.5854102*clhs124 + clhs127 + clhs129 + 0.5854102*clhs131 + clhs134 + clhs136 + clhs225 + clhs226 + clhs227 + clhs247 + clhs248 + clhs249 + clhs253 + clhs254 + clhs255 + clhs256 + 0.5854102*clhs281 + clhs284 + clhs286 + 0.5854102*clhs288 + clhs291 + clhs293 + clhs294 + clhs295 + clhs99;
            lhs(3,2)=clhs100 + clhs101 + clhs102 + clhs103 + clhs104 + clhs105 + clhs106 + clhs107 + clhs108 + clhs109 + clhs110 + 0.5854102*clhs153 + clhs156 + clhs158 + 0.5854102*clhs160 + clhs163 + clhs165 + 0.5854102*clhs167 + clhs170 + clhs172 + clhs236 + clhs237 + clhs238 + clhs273 + clhs274 + clhs275 + clhs276 + clhs277 + clhs278 + clhs279 + clhs282 + 0.5854102*clhs283 + clhs286 + clhs289 + 0.5854102*clhs290 + clhs293 + clhs294 + clhs295 + clhs99;
            lhs(3,3)=pow(DN_DX_3_0, 2)*clhs50 + pow(DN_DX_3_1, 2)*clhs50 + pow(DN_DX_3_2, 2)*clhs50 + pow(clhs185, 2)*tau[0] + pow(clhs189, 2)*tau[1] + pow(clhs193, 2)*tau[2] + 0.5854102*clhs194 + 0.5854102*clhs195 + 0.5854102*clhs196 + pow(clhs197, 2)*tau[3] + clhs200 + clhs201 + clhs203 + clhs204 + clhs206 + clhs207 + clhs250 + clhs251 + clhs252 + clhs282 + clhs284 + 0.5854102*clhs285 + clhs289 + clhs291 + 0.5854102*clhs292 + clhs294 + clhs295 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49;


    const double crhs0 =             0.19999999899376*f[1];
const double crhs1 =             0.1381966*v(0,0);
const double crhs2 =             0.1381966*v(2,0);
const double crhs3 =             0.1381966*v(3,0);
const double crhs4 =             crhs2 + crhs3;
const double crhs5 =             crhs1 + crhs4 + 0.5854102*v(1,0);
const double crhs6 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs7 =             crhs5*crhs6;
const double crhs8 =             -0.1381966*crhs7;
const double crhs9 =             0.1381966*v(0,1);
const double crhs10 =             0.1381966*v(2,1);
const double crhs11 =             0.1381966*v(3,1);
const double crhs12 =             crhs10 + crhs11;
const double crhs13 =             crhs12 + crhs9 + 0.5854102*v(1,1);
const double crhs14 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs15 =             crhs13*crhs14;
const double crhs16 =             -0.1381966*crhs15;
const double crhs17 =             0.1381966*v(0,2);
const double crhs18 =             0.1381966*v(2,2);
const double crhs19 =             0.1381966*v(3,2);
const double crhs20 =             crhs18 + crhs19;
const double crhs21 =             crhs17 + crhs20 + 0.5854102*v(1,2);
const double crhs22 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs23 =             crhs21*crhs22;
const double crhs24 =             -0.1381966*crhs23;
const double crhs25 =             0.1381966*phi[0];
const double crhs26 =             0.5854102*phi[1];
const double crhs27 =             0.1381966*phi[2];
const double crhs28 =             0.1381966*phi[3];
const double crhs29 =             crhs27 + crhs28;
const double crhs30 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs31 =             crhs30*(crhs25 + crhs26 + crhs29);
const double crhs32 =             -0.1381966*crhs31;
const double crhs33 =             0.19999999899376*f[2];
const double crhs34 =             0.19999999899376*f[3];
const double crhs35 =             4*crhs6*k;
const double crhs36 =             4*crhs14*k;
const double crhs37 =             4*crhs22*k;
const double crhs38 =             0.1381966*v(1,0);
const double crhs39 =             crhs1 + crhs38;
const double crhs40 =             crhs2 + crhs39 + 0.5854102*v(3,0);
const double crhs41 =             crhs40*crhs6;
const double crhs42 =             -0.1381966*crhs41;
const double crhs43 =             crhs3 + crhs39 + 0.5854102*v(2,0);
const double crhs44 =             crhs43*crhs6;
const double crhs45 =             -0.1381966*crhs44;
const double crhs46 =             crhs38 + crhs4 + 0.5854102*v(0,0);
const double crhs47 =             crhs46*crhs6;
const double crhs48 =             0.1381966*v(1,1);
const double crhs49 =             crhs48 + crhs9;
const double crhs50 =             crhs10 + crhs49 + 0.5854102*v(3,1);
const double crhs51 =             crhs14*crhs50;
const double crhs52 =             -0.1381966*crhs51;
const double crhs53 =             crhs11 + crhs49 + 0.5854102*v(2,1);
const double crhs54 =             crhs14*crhs53;
const double crhs55 =             -0.1381966*crhs54;
const double crhs56 =             crhs12 + crhs48 + 0.5854102*v(0,1);
const double crhs57 =             crhs14*crhs56;
const double crhs58 =             0.1381966*v(1,2);
const double crhs59 =             crhs17 + crhs58;
const double crhs60 =             crhs18 + crhs59 + 0.5854102*v(3,2);
const double crhs61 =             crhs22*crhs60;
const double crhs62 =             -0.1381966*crhs61;
const double crhs63 =             crhs19 + crhs59 + 0.5854102*v(2,2);
const double crhs64 =             crhs22*crhs63;
const double crhs65 =             -0.1381966*crhs64;
const double crhs66 =             crhs20 + crhs58 + 0.5854102*v(0,2);
const double crhs67 =             crhs22*crhs66;
const double crhs68 =             0.1381966*phi[1];
const double crhs69 =             crhs25 + crhs68;
const double crhs70 =             0.5854102*phi[3];
const double crhs71 =             crhs30*(crhs27 + crhs69 + crhs70);
const double crhs72 =             -0.1381966*crhs71;
const double crhs73 =             0.5854102*phi[2];
const double crhs74 =             crhs30*(crhs28 + crhs69 + crhs73);
const double crhs75 =             -0.1381966*crhs74;
const double crhs76 =             0.5854102*phi[0];
const double crhs77 =             crhs30*(crhs29 + crhs68 + crhs76);
const double crhs78 =             0.1381966*f[1];
const double crhs79 =             0.1381966*f[2];
const double crhs80 =             0.1381966*f[3];
const double crhs81 =             crhs79 + crhs80;
const double crhs82 =             crhs78 + crhs81 + 0.5854102*f[0];
const double crhs83 =             tau[0]*(DN_DX_0_0*crhs46 + DN_DX_0_1*crhs56 + DN_DX_0_2*crhs66);
const double crhs84 =             0.1381966*prj[1];
const double crhs85 =             0.1381966*prj[2];
const double crhs86 =             0.1381966*prj[3];
const double crhs87 =             crhs85 + crhs86;
const double crhs88 =             crhs84 + crhs87 + 0.5854102*prj[0];
const double crhs89 =             0.1381966*f[0];
const double crhs90 =             crhs81 + crhs89 + 0.5854102*f[1];
const double crhs91 =             tau[1]*(DN_DX_0_0*crhs5 + DN_DX_0_1*crhs13 + DN_DX_0_2*crhs21);
const double crhs92 =             0.1381966*prj[0];
const double crhs93 =             crhs87 + crhs92 + 0.5854102*prj[1];
const double crhs94 =             crhs78 + crhs89;
const double crhs95 =             crhs80 + crhs94 + 0.5854102*f[2];
const double crhs96 =             tau[2]*(DN_DX_0_0*crhs43 + DN_DX_0_1*crhs53 + DN_DX_0_2*crhs63);
const double crhs97 =             crhs84 + crhs92;
const double crhs98 =             crhs86 + crhs97 + 0.5854102*prj[2];
const double crhs99 =             crhs79 + crhs94 + 0.5854102*f[3];
const double crhs100 =             tau[3]*(DN_DX_0_0*crhs40 + DN_DX_0_1*crhs50 + DN_DX_0_2*crhs60);
const double crhs101 =             crhs85 + crhs97 + 0.5854102*prj[3];
const double crhs102 =             1.0/RK_time_coefficient;
const double crhs103 =             1.0/delta_time;
const double crhs104 =             -0.1381966*phi_old[2];
const double crhs105 =             -0.1381966*phi_old[3];
const double crhs106 =             crhs104 + crhs105;
const double crhs107 =             -0.1381966*phi_old[1];
const double crhs108 =             crhs102*crhs103*(crhs106 + crhs107 + crhs27 + crhs28 + crhs68 + crhs76 - 0.5854102*phi_old[0]);
const double crhs109 =             -0.1381966*phi_old[0];
const double crhs110 =             crhs102*crhs103*(crhs106 + crhs109 + crhs25 + crhs26 + crhs27 + crhs28 - 0.5854102*phi_old[1]);
const double crhs111 =             crhs107 + crhs109;
const double crhs112 =             crhs102*crhs103*(crhs105 + crhs111 + crhs25 + crhs28 + crhs68 + crhs73 - 0.5854102*phi_old[2]);
const double crhs113 =             crhs102*crhs103*(crhs104 + crhs111 + crhs25 + crhs27 + crhs68 + crhs70 - 0.5854102*phi_old[3]);
const double crhs114 =             crhs47 + crhs57 + crhs67;
const double crhs115 =             crhs15 + crhs23 + crhs7;
const double crhs116 =             crhs44 + crhs54 + crhs64;
const double crhs117 =             crhs41 + crhs51 + crhs61;
const double crhs118 =             0.19999999899376*f[0];
const double crhs119 =             -0.1381966*crhs47;
const double crhs120 =             -0.1381966*crhs57;
const double crhs121 =             -0.1381966*crhs67;
const double crhs122 =             -0.1381966*crhs77;
const double crhs123 =             tau[0]*(DN_DX_1_0*crhs46 + DN_DX_1_1*crhs56 + DN_DX_1_2*crhs66);
const double crhs124 =             tau[1]*(DN_DX_1_0*crhs5 + DN_DX_1_1*crhs13 + DN_DX_1_2*crhs21);
const double crhs125 =             tau[2]*(DN_DX_1_0*crhs43 + DN_DX_1_1*crhs53 + DN_DX_1_2*crhs63);
const double crhs126 =             tau[3]*(DN_DX_1_0*crhs40 + DN_DX_1_1*crhs50 + DN_DX_1_2*crhs60);
const double crhs127 =             crhs0 + crhs118 + crhs119 + crhs120 + crhs121 + crhs122 + crhs16 + crhs24 + crhs32 + crhs8;
const double crhs128 =             tau[0]*(DN_DX_2_0*crhs46 + DN_DX_2_1*crhs56 + DN_DX_2_2*crhs66);
const double crhs129 =             tau[1]*(DN_DX_2_0*crhs5 + DN_DX_2_1*crhs13 + DN_DX_2_2*crhs21);
const double crhs130 =             tau[2]*(DN_DX_2_0*crhs43 + DN_DX_2_1*crhs53 + DN_DX_2_2*crhs63);
const double crhs131 =             tau[3]*(DN_DX_2_0*crhs40 + DN_DX_2_1*crhs50 + DN_DX_2_2*crhs60);
const double crhs132 =             tau[0]*(DN_DX_3_0*crhs46 + DN_DX_3_1*crhs56 + DN_DX_3_2*crhs66);
const double crhs133 =             tau[1]*(DN_DX_3_0*crhs5 + DN_DX_3_1*crhs13 + DN_DX_3_2*crhs21);
const double crhs134 =             tau[2]*(DN_DX_3_0*crhs43 + DN_DX_3_1*crhs53 + DN_DX_3_2*crhs63);
const double crhs135 =             tau[3]*(DN_DX_3_0*crhs40 + DN_DX_3_1*crhs50 + DN_DX_3_2*crhs60);
            rhs[0]=-DN_DX_0_0*crhs35 - DN_DX_0_1*crhs36 - DN_DX_0_2*crhs37 + crhs0 + crhs100*crhs101 - crhs100*crhs113 - crhs100*crhs117 - crhs100*crhs71 + crhs100*crhs99 - crhs108*crhs83 - crhs110*crhs91 - crhs112*crhs96 - crhs114*crhs83 - crhs115*crhs91 - crhs116*crhs96 + crhs16 + crhs24 - crhs31*crhs91 + crhs32 + crhs33 + crhs34 + crhs42 + crhs45 - 0.5854102*crhs47 + crhs52 + crhs55 - 0.5854102*crhs57 + crhs62 + crhs65 - 0.5854102*crhs67 + crhs72 - crhs74*crhs96 + crhs75 - crhs77*crhs83 - 0.5854102*crhs77 + crhs8 + crhs82*crhs83 + crhs83*crhs88 + crhs90*crhs91 + crhs91*crhs93 + crhs95*crhs96 + crhs96*crhs98 + 0.40000000301872*f[0];
            rhs[1]=-DN_DX_1_0*crhs35 - DN_DX_1_1*crhs36 - DN_DX_1_2*crhs37 + crhs101*crhs126 - crhs108*crhs123 - crhs110*crhs124 - crhs112*crhs125 - crhs113*crhs126 - crhs114*crhs123 - crhs115*crhs124 - crhs116*crhs125 - crhs117*crhs126 + crhs118 + crhs119 + crhs120 + crhs121 + crhs122 - crhs123*crhs77 + crhs123*crhs82 + crhs123*crhs88 - crhs124*crhs31 + crhs124*crhs90 + crhs124*crhs93 - crhs125*crhs74 + crhs125*crhs95 + crhs125*crhs98 - crhs126*crhs71 + crhs126*crhs99 - 0.5854102*crhs15 - 0.5854102*crhs23 - 0.5854102*crhs31 + crhs33 + crhs34 + crhs42 + crhs45 + crhs52 + crhs55 + crhs62 + crhs65 - 0.5854102*crhs7 + crhs72 + crhs75 + 0.40000000301872*f[1];
            rhs[2]=-DN_DX_2_0*crhs35 - DN_DX_2_1*crhs36 - DN_DX_2_2*crhs37 + crhs101*crhs131 - crhs108*crhs128 - crhs110*crhs129 - crhs112*crhs130 - crhs113*crhs131 - crhs114*crhs128 - crhs115*crhs129 - crhs116*crhs130 - crhs117*crhs131 + crhs127 - crhs128*crhs77 + crhs128*crhs82 + crhs128*crhs88 - crhs129*crhs31 + crhs129*crhs90 + crhs129*crhs93 - crhs130*crhs74 + crhs130*crhs95 + crhs130*crhs98 - crhs131*crhs71 + crhs131*crhs99 + crhs42 - 0.5854102*crhs44 + crhs52 - 0.5854102*crhs54 + crhs62 - 0.5854102*crhs64 + crhs72 - 0.5854102*crhs74 + 0.40000000301872*f[2] + 0.19999999899376*f[3];
            rhs[3]=-DN_DX_3_0*crhs35 - DN_DX_3_1*crhs36 - DN_DX_3_2*crhs37 + crhs101*crhs135 - crhs108*crhs132 - crhs110*crhs133 - crhs112*crhs134 - crhs113*crhs135 - crhs114*crhs132 - crhs115*crhs133 - crhs116*crhs134 - crhs117*crhs135 + crhs127 - crhs132*crhs77 + crhs132*crhs82 + crhs132*crhs88 - crhs133*crhs31 + crhs133*crhs90 + crhs133*crhs93 - crhs134*crhs74 + crhs134*crhs95 + crhs134*crhs98 - crhs135*crhs71 + crhs135*crhs99 - 0.5854102*crhs41 + crhs45 - 0.5854102*crhs51 + crhs55 - 0.5854102*crhs61 + crhs65 - 0.5854102*crhs71 + crhs75 + 0.19999999899376*f[2] + 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<2>::CalculateOrthogonalSubgridScaleSystemInternal(
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

    const double crhs0 =             -0.25*f[1];
const double crhs1 =             0.166666666666667*v(0,0);
const double crhs2 =             0.166666666666667*v(2,0);
const double crhs3 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs4 =             crhs3*(crhs1 + crhs2 + 0.666666666666667*v(1,0));
const double crhs5 =             0.166666666666667*crhs4;
const double crhs6 =             0.166666666666667*v(0,1);
const double crhs7 =             0.166666666666667*v(2,1);
const double crhs8 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs9 =             crhs8*(crhs6 + crhs7 + 0.666666666666667*v(1,1));
const double crhs10 =             0.166666666666667*crhs9;
const double crhs11 =             1/(RK_time_coefficient*delta_time);
const double crhs12 =             0.111111111111111*phi[1];
const double crhs13 =             -0.111111111111111*phi_old[1];
const double crhs14 =             crhs12 + crhs13;
const double crhs15 =             0.0277777777777778*phi[0];
const double crhs16 =             0.0277777777777778*phi[2];
const double crhs17 =             -0.0277777777777778*phi_old[0];
const double crhs18 =             -0.0277777777777778*phi_old[2];
const double crhs19 =             crhs11*(crhs14 + crhs15 + crhs16 + crhs17 + crhs18);
const double crhs20 =             0.166666666666667*phi[0];
const double crhs21 =             0.166666666666667*phi[2];
const double crhs22 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs23 =             crhs22*(crhs20 + crhs21 + 0.666666666666667*phi[1]);
const double crhs24 =             0.166666666666667*crhs23;
const double crhs25 =             -0.25*f[2];
const double crhs26 =             3*crhs3*k;
const double crhs27 =             3*crhs8*k;
const double crhs28 =             0.166666666666667*v(1,0);
const double crhs29 =             crhs3*(crhs1 + crhs28 + 0.666666666666667*v(2,0));
const double crhs30 =             0.166666666666667*crhs29;
const double crhs31 =             crhs3*(crhs2 + crhs28 + 0.666666666666667*v(0,0));
const double crhs32 =             0.166666666666667*v(1,1);
const double crhs33 =             crhs8*(crhs32 + crhs6 + 0.666666666666667*v(2,1));
const double crhs34 =             0.166666666666667*crhs33;
const double crhs35 =             crhs8*(crhs32 + crhs7 + 0.666666666666667*v(0,1));
const double crhs36 =             0.111111111111111*phi[2];
const double crhs37 =             -0.111111111111111*phi_old[2];
const double crhs38 =             0.0277777777777778*phi[1];
const double crhs39 =             -0.0277777777777778*phi_old[1];
const double crhs40 =             crhs11*(crhs15 + crhs17 + crhs36 + crhs37 + crhs38 + crhs39);
const double crhs41 =             0.166666666666667*phi[1];
const double crhs42 =             crhs22*(crhs20 + crhs41 + 0.666666666666667*phi[2]);
const double crhs43 =             0.166666666666667*crhs42;
const double crhs44 =             crhs22*(crhs21 + crhs41 + 0.666666666666667*phi[0]);
const double crhs45 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs46 =             crhs11*(crhs16 + crhs18 + crhs38 + crhs39 + crhs45) + 0.166666666666667*crhs31 + 0.166666666666667*crhs35 + 0.166666666666667*crhs44 - 0.25*f[0];
            rhs[0]=DN_DX_0_0*crhs26 + DN_DX_0_1*crhs27 + crhs0 + crhs10 + crhs11*(crhs14 + crhs36 + crhs37 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs19 + crhs24 + crhs25 + crhs30 + 0.666666666666667*crhs31 + crhs34 + 0.666666666666667*crhs35 + crhs40 + crhs43 + 0.666666666666667*crhs44 + crhs5 - 0.5*f[0];
            rhs[1]=DN_DX_1_0*crhs26 + DN_DX_1_1*crhs27 + crhs11*(crhs36 + crhs37 + crhs45 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) + 0.666666666666667*crhs23 + crhs25 + crhs30 + crhs34 + 0.666666666666667*crhs4 + crhs40 + crhs43 + crhs46 + 0.666666666666667*crhs9 - 0.5*f[1];
            rhs[2]=DN_DX_2_0*crhs26 + DN_DX_2_1*crhs27 + crhs0 + crhs10 + crhs11*(crhs12 + crhs13 + crhs45 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs19 + crhs24 + 0.666666666666667*crhs29 + 0.666666666666667*crhs33 + 0.666666666666667*crhs42 + crhs46 + crhs5 - 0.5*f[2];


    const double local_size = 3;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3>::CalculateOrthogonalSubgridScaleSystemInternal(
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

    const double crhs0 =             -0.19999999899376*f[1];
const double crhs1 =             0.1381966*v(0,0);
const double crhs2 =             0.1381966*v(2,0);
const double crhs3 =             0.1381966*v(3,0);
const double crhs4 =             crhs2 + crhs3;
const double crhs5 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs6 =             crhs5*(crhs1 + crhs4 + 0.5854102*v(1,0));
const double crhs7 =             0.1381966*crhs6;
const double crhs8 =             0.1381966*v(0,1);
const double crhs9 =             0.1381966*v(2,1);
const double crhs10 =             0.1381966*v(3,1);
const double crhs11 =             crhs10 + crhs9;
const double crhs12 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs13 =             crhs12*(crhs11 + crhs8 + 0.5854102*v(1,1));
const double crhs14 =             0.1381966*crhs13;
const double crhs15 =             0.1381966*v(0,2);
const double crhs16 =             0.1381966*v(2,2);
const double crhs17 =             0.1381966*v(3,2);
const double crhs18 =             crhs16 + crhs17;
const double crhs19 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs20 =             crhs19*(crhs15 + crhs18 + 0.5854102*v(1,2));
const double crhs21 =             0.1381966*crhs20;
const double crhs22 =             1/(RK_time_coefficient*delta_time);
const double crhs23 =             0.08090169924532*phi[1];
const double crhs24 =             -0.08090169924532*phi_old[1];
const double crhs25 =             0.01909830025156*phi[0];
const double crhs26 =             0.01909830025156*phi[2];
const double crhs27 =             0.01909830025156*phi[3];
const double crhs28 =             -0.01909830025156*phi_old[0];
const double crhs29 =             -0.01909830025156*phi_old[2];
const double crhs30 =             -0.01909830025156*phi_old[3];
const double crhs31 =             crhs22*(crhs23 + crhs24 + crhs25 + crhs26 + crhs27 + crhs28 + crhs29 + crhs30);
const double crhs32 =             0.1381966*phi[0];
const double crhs33 =             0.1381966*phi[2];
const double crhs34 =             0.1381966*phi[3];
const double crhs35 =             crhs33 + crhs34;
const double crhs36 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs37 =             crhs36*(crhs32 + crhs35 + 0.5854102*phi[1]);
const double crhs38 =             0.1381966*crhs37;
const double crhs39 =             -0.19999999899376*f[2];
const double crhs40 =             -0.19999999899376*f[3];
const double crhs41 =             4*crhs5*k;
const double crhs42 =             4*crhs12*k;
const double crhs43 =             4*crhs19*k;
const double crhs44 =             0.1381966*v(1,0);
const double crhs45 =             crhs1 + crhs44;
const double crhs46 =             crhs5*(crhs2 + crhs45 + 0.5854102*v(3,0));
const double crhs47 =             0.1381966*crhs46;
const double crhs48 =             crhs5*(crhs3 + crhs45 + 0.5854102*v(2,0));
const double crhs49 =             0.1381966*crhs48;
const double crhs50 =             crhs5*(crhs4 + crhs44 + 0.5854102*v(0,0));
const double crhs51 =             0.1381966*v(1,1);
const double crhs52 =             crhs51 + crhs8;
const double crhs53 =             crhs12*(crhs52 + crhs9 + 0.5854102*v(3,1));
const double crhs54 =             0.1381966*crhs53;
const double crhs55 =             crhs12*(crhs10 + crhs52 + 0.5854102*v(2,1));
const double crhs56 =             0.1381966*crhs55;
const double crhs57 =             crhs12*(crhs11 + crhs51 + 0.5854102*v(0,1));
const double crhs58 =             0.1381966*v(1,2);
const double crhs59 =             crhs15 + crhs58;
const double crhs60 =             crhs19*(crhs16 + crhs59 + 0.5854102*v(3,2));
const double crhs61 =             0.1381966*crhs60;
const double crhs62 =             crhs19*(crhs17 + crhs59 + 0.5854102*v(2,2));
const double crhs63 =             0.1381966*crhs62;
const double crhs64 =             crhs19*(crhs18 + crhs58 + 0.5854102*v(0,2));
const double crhs65 =             0.08090169924532*phi[3];
const double crhs66 =             -0.08090169924532*phi_old[3];
const double crhs67 =             0.01909830025156*phi[1];
const double crhs68 =             -0.01909830025156*phi_old[1];
const double crhs69 =             crhs22*(crhs25 + crhs26 + crhs28 + crhs29 + crhs65 + crhs66 + crhs67 + crhs68);
const double crhs70 =             0.08090169924532*phi[2];
const double crhs71 =             -0.08090169924532*phi_old[2];
const double crhs72 =             crhs22*(crhs25 + crhs27 + crhs28 + crhs30 + crhs67 + crhs68 + crhs70 + crhs71);
const double crhs73 =             crhs65 + crhs66 + crhs70 + crhs71;
const double crhs74 =             0.1381966*phi[1];
const double crhs75 =             crhs32 + crhs74;
const double crhs76 =             crhs36*(crhs33 + crhs75 + 0.5854102*phi[3]);
const double crhs77 =             0.1381966*crhs76;
const double crhs78 =             crhs36*(crhs34 + crhs75 + 0.5854102*phi[2]);
const double crhs79 =             0.1381966*crhs78;
const double crhs80 =             crhs36*(crhs35 + crhs74 + 0.5854102*phi[0]);
const double crhs81 =             -0.19999999899376*f[0];
const double crhs82 =             0.1381966*crhs50;
const double crhs83 =             0.1381966*crhs57;
const double crhs84 =             0.1381966*crhs64;
const double crhs85 =             0.08090169924532*phi[0];
const double crhs86 =             -0.08090169924532*phi_old[0];
const double crhs87 =             crhs22*(crhs26 + crhs27 + crhs29 + crhs30 + crhs67 + crhs68 + crhs85 + crhs86);
const double crhs88 =             0.1381966*crhs80;
const double crhs89 =             crhs0 + crhs14 + crhs21 + crhs31 + crhs38 + crhs7 + crhs81 + crhs82 + crhs83 + crhs84 + crhs87 + crhs88;
const double crhs90 =             crhs23 + crhs24 + crhs85 + crhs86;
            rhs[0]=DN_DX_0_0*crhs41 + DN_DX_0_1*crhs42 + DN_DX_0_2*crhs43 + crhs0 + crhs14 + crhs21 + crhs22*(crhs23 + crhs24 + crhs73 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs31 + crhs38 + crhs39 + crhs40 + crhs47 + crhs49 + 0.5854102*crhs50 + crhs54 + crhs56 + 0.5854102*crhs57 + crhs61 + crhs63 + 0.5854102*crhs64 + crhs69 + crhs7 + crhs72 + crhs77 + crhs79 + 0.5854102*crhs80 - 0.40000000301872*f[0];
            rhs[1]=DN_DX_1_0*crhs41 + DN_DX_1_1*crhs42 + DN_DX_1_2*crhs43 + 0.5854102*crhs13 + 0.5854102*crhs20 + crhs22*(crhs73 + crhs85 + crhs86 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) + 0.5854102*crhs37 + crhs39 + crhs40 + crhs47 + crhs49 + crhs54 + crhs56 + 0.5854102*crhs6 + crhs61 + crhs63 + crhs69 + crhs72 + crhs77 + crhs79 + crhs81 + crhs82 + crhs83 + crhs84 + crhs87 + crhs88 - 0.40000000301872*f[1];
            rhs[2]=DN_DX_2_0*crhs41 + DN_DX_2_1*crhs42 + DN_DX_2_2*crhs43 + crhs22*(crhs65 + crhs66 + crhs90 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) + crhs47 + 0.5854102*crhs48 + crhs54 + 0.5854102*crhs55 + crhs61 + 0.5854102*crhs62 + crhs69 + crhs77 + 0.5854102*crhs78 + crhs89 - 0.40000000301872*f[2] - 0.19999999899376*f[3];
            rhs[3]=DN_DX_3_0*crhs41 + DN_DX_3_1*crhs42 + DN_DX_3_2*crhs43 + crhs22*(crhs70 + crhs71 + crhs90 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) + 0.5854102*crhs46 + crhs49 + 0.5854102*crhs53 + crhs56 + 0.5854102*crhs60 + crhs63 + crhs72 + 0.5854102*crhs76 + crhs79 + crhs89 - 0.19999999899376*f[2] - 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
double SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::ComputeH(
    BoundedMatrix<double,TNumNodes,TDim >& DN_DX)
{
    KRATOS_TRY;

    double h=0.0;
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        double h_inv = 0.0;
        for(unsigned int k=0; k<TDim; k++)
        {
            h_inv += DN_DX(i,k)*DN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    h = sqrt(h)/static_cast<double>(TNumNodes);
    return h;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
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
        inv_tau += rVariables.dynamic_tau/rVariables.delta_time;
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

template< unsigned int TDim, unsigned int TNumNodes >
Element::IntegrationMethod SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicQSConvectionDiffusionExplicit<2>;
template class SymbolicQSConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
