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
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
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
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::CalculateLocalSystemInternal(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau;
    const auto& qstau = rVariables.qstau;
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

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             0.166666666666667*v(0,0);
const double clhs3 =             0.166666666666667*v(2,0);
const double clhs4 =             clhs2 + clhs3 + 0.666666666666667*v(1,0);
const double clhs5 =             0.166666666666667*v(0,1);
const double clhs6 =             0.166666666666667*v(2,1);
const double clhs7 =             clhs5 + clhs6 + 0.666666666666667*v(1,1);
const double clhs8 =             DN_DX_0_0*clhs4 + DN_DX_0_1*clhs7;
const double clhs9 =             clhs0*clhs1*clhs8*tau[1];
const double clhs10 =             0.166666666666667*clhs9;
const double clhs11 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double clhs12 =             clhs11*clhs8*tau[1];
const double clhs13 =             0.166666666666667*clhs12;
const double clhs14 =             clhs0*clhs1;
const double clhs15 =             -0.5*clhs14;
const double clhs16 =             0.166666666666667*v(1,0);
const double clhs17 =             clhs16 + clhs3 + 0.666666666666667*v(0,0);
const double clhs18 =             0.166666666666667*v(1,1);
const double clhs19 =             clhs18 + clhs6 + 0.666666666666667*v(0,1);
const double clhs20 =             DN_DX_0_0*clhs17 + DN_DX_0_1*clhs19;
const double clhs21 =             clhs16 + clhs2 + 0.666666666666667*v(2,0);
const double clhs22 =             clhs18 + clhs5 + 0.666666666666667*v(2,1);
const double clhs23 =             DN_DX_0_0*clhs21 + DN_DX_0_1*clhs22;
const double clhs24 =             clhs0*clhs1*clhs20*tau[0];
const double clhs25 =             clhs0*clhs1*clhs23*tau[2];
const double clhs26 =             0.166666666666667*clhs25;
const double clhs27 =             clhs11*clhs20*tau[0];
const double clhs28 =             clhs11*clhs23*tau[2];
const double clhs29 =             0.166666666666667*clhs28;
const double clhs30 =             0.166666666666667*clhs24;
const double clhs31 =             0.166666666666667*clhs27;
const double clhs32 =             -0.25*clhs14;
const double clhs33 =             DN_DX_1_0*clhs17 + DN_DX_1_1*clhs19;
const double clhs34 =             clhs20*tau[0];
const double clhs35 =             clhs33*clhs34;
const double clhs36 =             DN_DX_1_0*clhs4 + DN_DX_1_1*clhs7;
const double clhs37 =             clhs8*tau[1];
const double clhs38 =             clhs36*clhs37;
const double clhs39 =             DN_DX_1_0*clhs21 + DN_DX_1_1*clhs22;
const double clhs40 =             clhs23*tau[2];
const double clhs41 =             clhs39*clhs40;
const double clhs42 =             DN_DX_2_0*clhs17 + DN_DX_2_1*clhs19;
const double clhs43 =             clhs34*clhs42;
const double clhs44 =             DN_DX_2_0*clhs4 + DN_DX_2_1*clhs7;
const double clhs45 =             clhs37*clhs44;
const double clhs46 =             DN_DX_2_0*clhs21 + DN_DX_2_1*clhs22;
const double clhs47 =             clhs40*clhs46;
const double clhs48 =             clhs0*clhs1*clhs36*tau[1];
const double clhs49 =             0.166666666666667*clhs48;
const double clhs50 =             clhs11*clhs36*tau[1];
const double clhs51 =             0.166666666666667*clhs50;
const double clhs52 =             clhs0*clhs1*clhs33*tau[0];
const double clhs53 =             clhs0*clhs1*clhs39*tau[2];
const double clhs54 =             0.166666666666667*clhs53;
const double clhs55 =             clhs11*clhs33*tau[0];
const double clhs56 =             clhs11*clhs39*tau[2];
const double clhs57 =             0.166666666666667*clhs56;
const double clhs58 =             0.166666666666667*clhs52;
const double clhs59 =             0.166666666666667*clhs55;
const double clhs60 =             clhs32 + clhs33*clhs42*tau[0] + clhs36*clhs44*tau[1] + clhs39*clhs46*tau[2];
const double clhs61 =             clhs0*clhs1*clhs46*tau[2];
const double clhs62 =             0.166666666666667*clhs61;
const double clhs63 =             clhs11*clhs46*tau[2];
const double clhs64 =             0.166666666666667*clhs63;
const double clhs65 =             clhs0*clhs1*clhs42*tau[0];
const double clhs66 =             clhs0*clhs1*clhs44*tau[1];
const double clhs67 =             0.166666666666667*clhs66;
const double clhs68 =             clhs11*clhs42*tau[0];
const double clhs69 =             clhs11*clhs44*tau[1];
const double clhs70 =             0.166666666666667*clhs69;
const double clhs71 =             0.166666666666667*clhs65;
const double clhs72 =             0.166666666666667*clhs68;
            lhs(0,0)=clhs10 + clhs13 + clhs15 + pow(clhs20, 2)*tau[0] + pow(clhs23, 2)*tau[2] + 0.666666666666667*clhs24 + clhs26 + 0.666666666666667*clhs27 + clhs29 + pow(clhs8, 2)*tau[1];
            lhs(0,1)=0.666666666666667*clhs12 + clhs26 + clhs29 + clhs30 + clhs31 + clhs32 + clhs35 + clhs38 + clhs41 + 0.666666666666667*clhs9;
            lhs(0,2)=clhs10 + clhs13 + 0.666666666666667*clhs25 + 0.666666666666667*clhs28 + clhs30 + clhs31 + clhs32 + clhs43 + clhs45 + clhs47;
            lhs(1,0)=clhs32 + clhs35 + clhs38 + clhs41 + clhs49 + clhs51 + 0.666666666666667*clhs52 + clhs54 + 0.666666666666667*clhs55 + clhs57;
            lhs(1,1)=clhs15 + pow(clhs33, 2)*tau[0] + pow(clhs36, 2)*tau[1] + pow(clhs39, 2)*tau[2] + 0.666666666666667*clhs48 + 0.666666666666667*clhs50 + clhs54 + clhs57 + clhs58 + clhs59;
            lhs(1,2)=clhs49 + clhs51 + 0.666666666666667*clhs53 + 0.666666666666667*clhs56 + clhs58 + clhs59 + clhs60;
            lhs(2,0)=clhs32 + clhs43 + clhs45 + clhs47 + clhs62 + clhs64 + 0.666666666666667*clhs65 + clhs67 + 0.666666666666667*clhs68 + clhs70;
            lhs(2,1)=clhs60 + clhs62 + clhs64 + 0.666666666666667*clhs66 + 0.666666666666667*clhs69 + clhs71 + clhs72;
            lhs(2,2)=clhs15 + pow(clhs42, 2)*tau[0] + pow(clhs44, 2)*tau[1] + pow(clhs46, 2)*tau[2] + 0.666666666666667*clhs61 + 0.666666666666667*clhs63 + clhs67 + clhs70 + clhs71 + clhs72;


    const double crhs0 =             phi_subscale_gauss[2]/qstau[2];
const double crhs1 =             1.0/RK_time_coefficient;
const double crhs2 =             1.0/delta_time;
const double crhs3 =             crhs1*crhs2;
const double crhs4 =             0.111111111111111*phi[2];
const double crhs5 =             -0.111111111111111*phi_old[2];
const double crhs6 =             0.0277777777777778*phi[0];
const double crhs7 =             0.0277777777777778*phi[1];
const double crhs8 =             -0.0277777777777778*phi_old[0];
const double crhs9 =             -0.0277777777777778*phi_old[1];
const double crhs10 =             0.166666666666667*crhs0 + crhs3*(crhs4 + crhs5 + crhs6 + crhs7 + crhs8 + crhs9) - 0.25*prj[2];
const double crhs11 =             -0.25*prj[1];
const double crhs12 =             phi_subscale_gauss[0]/qstau[0];
const double crhs13 =             phi_subscale_gauss[1]/qstau[1];
const double crhs14 =             0.166666666666667*crhs13;
const double crhs15 =             0.111111111111111*phi[1];
const double crhs16 =             -0.111111111111111*phi_old[1];
const double crhs17 =             crhs15 + crhs16;
const double crhs18 =             0.0277777777777778*phi[2];
const double crhs19 =             -0.0277777777777778*phi_old[2];
const double crhs20 =             crhs3*(crhs17 + crhs18 + crhs19 + crhs6 + crhs8);
const double crhs21 =             crhs2*phi_subscale_gauss[0]*tau[0];
const double crhs22 =             0.166666666666667*v(1,0);
const double crhs23 =             0.166666666666667*v(2,0);
const double crhs24 =             crhs22 + crhs23 + 0.666666666666667*v(0,0);
const double crhs25 =             0.166666666666667*v(1,1);
const double crhs26 =             0.166666666666667*v(2,1);
const double crhs27 =             crhs25 + crhs26 + 0.666666666666667*v(0,1);
const double crhs28 =             DN_DX_0_0*crhs24 + DN_DX_0_1*crhs27;
const double crhs29 =             crhs2*phi_subscale_gauss[1]*tau[1];
const double crhs30 =             0.166666666666667*v(0,0);
const double crhs31 =             crhs23 + crhs30 + 0.666666666666667*v(1,0);
const double crhs32 =             0.166666666666667*v(0,1);
const double crhs33 =             crhs26 + crhs32 + 0.666666666666667*v(1,1);
const double crhs34 =             DN_DX_0_0*crhs31 + DN_DX_0_1*crhs33;
const double crhs35 =             crhs2*phi_subscale_gauss[2]*tau[2];
const double crhs36 =             crhs22 + crhs30 + 0.666666666666667*v(2,0);
const double crhs37 =             crhs25 + crhs32 + 0.666666666666667*v(2,1);
const double crhs38 =             DN_DX_0_0*crhs36 + DN_DX_0_1*crhs37;
const double crhs39 =             0.166666666666667*f[1];
const double crhs40 =             0.166666666666667*f[2];
const double crhs41 =             crhs39 + crhs40 + 0.666666666666667*f[0];
const double crhs42 =             crhs28*tau[0];
const double crhs43 =             0.166666666666667*prj[1];
const double crhs44 =             0.166666666666667*prj[2];
const double crhs45 =             crhs43 + crhs44 + 0.666666666666667*prj[0];
const double crhs46 =             0.166666666666667*f[0];
const double crhs47 =             crhs40 + crhs46 + 0.666666666666667*f[1];
const double crhs48 =             crhs34*tau[1];
const double crhs49 =             0.166666666666667*prj[0];
const double crhs50 =             crhs44 + crhs49 + 0.666666666666667*prj[1];
const double crhs51 =             crhs39 + crhs46 + 0.666666666666667*f[2];
const double crhs52 =             crhs38*tau[2];
const double crhs53 =             crhs43 + crhs49 + 0.666666666666667*prj[2];
const double crhs54 =             0.166666666666667*phi[1];
const double crhs55 =             0.166666666666667*phi[2];
const double crhs56 =             crhs54 + crhs55 + 0.666666666666667*phi[0];
const double crhs57 =             -0.166666666666667*phi_old[1];
const double crhs58 =             -0.166666666666667*phi_old[2];
const double crhs59 =             crhs1*crhs2*(crhs56 + crhs57 + crhs58 - 0.666666666666667*phi_old[0]);
const double crhs60 =             0.166666666666667*phi[0];
const double crhs61 =             crhs55 + crhs60 + 0.666666666666667*phi[1];
const double crhs62 =             -0.166666666666667*phi_old[0];
const double crhs63 =             crhs1*crhs2*(crhs58 + crhs61 + crhs62 - 0.666666666666667*phi_old[1]);
const double crhs64 =             crhs54 + crhs60 + 0.666666666666667*phi[2];
const double crhs65 =             crhs1*crhs2*(crhs57 + crhs62 + crhs64 - 0.666666666666667*phi_old[2]);
const double crhs66 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs67 =             crhs56*crhs66;
const double crhs68 =             crhs61*crhs66;
const double crhs69 =             crhs64*crhs66;
const double crhs70 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs71 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs72 =             crhs24*crhs70 + crhs27*crhs71;
const double crhs73 =             crhs31*crhs70 + crhs33*crhs71;
const double crhs74 =             crhs36*crhs70 + crhs37*crhs71;
const double crhs75 =             -0.25*prj[0];
const double crhs76 =             0.166666666666667*crhs12;
const double crhs77 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs78 =             crhs3*(crhs18 + crhs19 + crhs7 + crhs77 + crhs9);
const double crhs79 =             DN_DX_1_0*crhs24 + DN_DX_1_1*crhs27;
const double crhs80 =             DN_DX_1_0*crhs31 + DN_DX_1_1*crhs33;
const double crhs81 =             DN_DX_1_0*crhs36 + DN_DX_1_1*crhs37;
const double crhs82 =             crhs79*tau[0];
const double crhs83 =             crhs80*tau[1];
const double crhs84 =             crhs81*tau[2];
const double crhs85 =             DN_DX_2_0*crhs24 + DN_DX_2_1*crhs27;
const double crhs86 =             DN_DX_2_0*crhs31 + DN_DX_2_1*crhs33;
const double crhs87 =             DN_DX_2_0*crhs36 + DN_DX_2_1*crhs37;
const double crhs88 =             crhs85*tau[0];
const double crhs89 =             crhs86*tau[1];
const double crhs90 =             crhs87*tau[2];
            rhs[0]=crhs10 + crhs11 + 0.666666666666667*crhs12 + crhs14 + crhs20 + crhs21*crhs28 + crhs29*crhs34 + crhs3*(crhs17 + crhs4 + crhs5 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs35*crhs38 + crhs41*crhs42 + crhs42*crhs45 - crhs42*crhs59 - crhs42*crhs67 - crhs42*crhs72 + crhs47*crhs48 + crhs48*crhs50 - crhs48*crhs63 - crhs48*crhs68 - crhs48*crhs73 + crhs51*crhs52 + crhs52*crhs53 - crhs52*crhs65 - crhs52*crhs69 - crhs52*crhs74 - 0.5*prj[0];
            rhs[1]=crhs10 + 0.666666666666667*crhs13 + crhs21*crhs79 + crhs29*crhs80 + crhs3*(crhs4 + crhs5 + crhs77 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) + crhs35*crhs81 + crhs41*crhs82 + crhs45*crhs82 + crhs47*crhs83 + crhs50*crhs83 + crhs51*crhs84 + crhs53*crhs84 - crhs59*crhs82 - crhs63*crhs83 - crhs65*crhs84 - crhs67*crhs82 - crhs68*crhs83 - crhs69*crhs84 - crhs72*crhs82 - crhs73*crhs83 - crhs74*crhs84 + crhs75 + crhs76 + crhs78 - 0.5*prj[1];
            rhs[2]=0.666666666666667*crhs0 + crhs11 + crhs14 + crhs20 + crhs21*crhs85 + crhs29*crhs86 + crhs3*(crhs15 + crhs16 + crhs77 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs35*crhs87 + crhs41*crhs88 + crhs45*crhs88 + crhs47*crhs89 + crhs50*crhs89 + crhs51*crhs90 + crhs53*crhs90 - crhs59*crhs88 - crhs63*crhs89 - crhs65*crhs90 - crhs67*crhs88 - crhs68*crhs89 - crhs69*crhs90 - crhs72*crhs88 - crhs73*crhs89 - crhs74*crhs90 + crhs75 + crhs76 + crhs78 - 0.5*prj[2];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::CalculateLocalSystemInternal(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto& v = rVariables.convective_velocity;
    const auto& tau = rVariables.tau;
    const auto& qstau = rVariables.qstau;
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

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             0.1381966*v(0,0);
const double clhs3 =             0.1381966*v(2,0);
const double clhs4 =             0.1381966*v(3,0);
const double clhs5 =             clhs3 + clhs4;
const double clhs6 =             clhs2 + clhs5 + 0.5854102*v(1,0);
const double clhs7 =             0.1381966*v(0,1);
const double clhs8 =             0.1381966*v(2,1);
const double clhs9 =             0.1381966*v(3,1);
const double clhs10 =             clhs8 + clhs9;
const double clhs11 =             clhs10 + clhs7 + 0.5854102*v(1,1);
const double clhs12 =             0.1381966*v(0,2);
const double clhs13 =             0.1381966*v(2,2);
const double clhs14 =             0.1381966*v(3,2);
const double clhs15 =             clhs13 + clhs14;
const double clhs16 =             clhs12 + clhs15 + 0.5854102*v(1,2);
const double clhs17 =             DN_DX_0_0*clhs6 + DN_DX_0_1*clhs11 + DN_DX_0_2*clhs16;
const double clhs18 =             clhs0*clhs1*clhs17*tau[1];
const double clhs19 =             0.1381966*clhs18;
const double clhs20 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double clhs21 =             clhs17*clhs20*tau[1];
const double clhs22 =             0.1381966*clhs21;
const double clhs23 =             clhs0*clhs1;
const double clhs24 =             -0.40000000301872*clhs23;
const double clhs25 =             0.1381966*v(1,0);
const double clhs26 =             clhs25 + clhs5 + 0.5854102*v(0,0);
const double clhs27 =             0.1381966*v(1,1);
const double clhs28 =             clhs10 + clhs27 + 0.5854102*v(0,1);
const double clhs29 =             0.1381966*v(1,2);
const double clhs30 =             clhs15 + clhs29 + 0.5854102*v(0,2);
const double clhs31 =             DN_DX_0_0*clhs26 + DN_DX_0_1*clhs28 + DN_DX_0_2*clhs30;
const double clhs32 =             clhs2 + clhs25;
const double clhs33 =             clhs32 + clhs4 + 0.5854102*v(2,0);
const double clhs34 =             clhs27 + clhs7;
const double clhs35 =             clhs34 + clhs9 + 0.5854102*v(2,1);
const double clhs36 =             clhs12 + clhs29;
const double clhs37 =             clhs14 + clhs36 + 0.5854102*v(2,2);
const double clhs38 =             DN_DX_0_0*clhs33 + DN_DX_0_1*clhs35 + DN_DX_0_2*clhs37;
const double clhs39 =             clhs3 + clhs32 + 0.5854102*v(3,0);
const double clhs40 =             clhs34 + clhs8 + 0.5854102*v(3,1);
const double clhs41 =             clhs13 + clhs36 + 0.5854102*v(3,2);
const double clhs42 =             DN_DX_0_0*clhs39 + DN_DX_0_1*clhs40 + DN_DX_0_2*clhs41;
const double clhs43 =             clhs0*clhs1*clhs31*tau[0];
const double clhs44 =             clhs0*clhs1*clhs38*tau[2];
const double clhs45 =             0.1381966*clhs44;
const double clhs46 =             clhs0*clhs1*clhs42*tau[3];
const double clhs47 =             0.1381966*clhs46;
const double clhs48 =             clhs20*clhs31*tau[0];
const double clhs49 =             clhs20*clhs38*tau[2];
const double clhs50 =             0.1381966*clhs49;
const double clhs51 =             clhs20*clhs42*tau[3];
const double clhs52 =             0.1381966*clhs51;
const double clhs53 =             0.1381966*clhs43;
const double clhs54 =             0.1381966*clhs48;
const double clhs55 =             -0.19999999899376*clhs23;
const double clhs56 =             DN_DX_1_0*clhs26 + DN_DX_1_1*clhs28 + DN_DX_1_2*clhs30;
const double clhs57 =             clhs31*tau[0];
const double clhs58 =             clhs56*clhs57;
const double clhs59 =             DN_DX_1_0*clhs6 + DN_DX_1_1*clhs11 + DN_DX_1_2*clhs16;
const double clhs60 =             clhs17*tau[1];
const double clhs61 =             clhs59*clhs60;
const double clhs62 =             DN_DX_1_0*clhs33 + DN_DX_1_1*clhs35 + DN_DX_1_2*clhs37;
const double clhs63 =             clhs38*tau[2];
const double clhs64 =             clhs62*clhs63;
const double clhs65 =             DN_DX_1_0*clhs39 + DN_DX_1_1*clhs40 + DN_DX_1_2*clhs41;
const double clhs66 =             clhs42*tau[3];
const double clhs67 =             clhs65*clhs66;
const double clhs68 =             DN_DX_2_0*clhs26 + DN_DX_2_1*clhs28 + DN_DX_2_2*clhs30;
const double clhs69 =             clhs57*clhs68;
const double clhs70 =             DN_DX_2_0*clhs6 + DN_DX_2_1*clhs11 + DN_DX_2_2*clhs16;
const double clhs71 =             clhs60*clhs70;
const double clhs72 =             DN_DX_2_0*clhs33 + DN_DX_2_1*clhs35 + DN_DX_2_2*clhs37;
const double clhs73 =             clhs63*clhs72;
const double clhs74 =             DN_DX_2_0*clhs39 + DN_DX_2_1*clhs40 + DN_DX_2_2*clhs41;
const double clhs75 =             clhs66*clhs74;
const double clhs76 =             DN_DX_3_0*clhs26 + DN_DX_3_1*clhs28 + DN_DX_3_2*clhs30;
const double clhs77 =             clhs57*clhs76;
const double clhs78 =             DN_DX_3_0*clhs6 + DN_DX_3_1*clhs11 + DN_DX_3_2*clhs16;
const double clhs79 =             clhs60*clhs78;
const double clhs80 =             DN_DX_3_0*clhs33 + DN_DX_3_1*clhs35 + DN_DX_3_2*clhs37;
const double clhs81 =             clhs63*clhs80;
const double clhs82 =             DN_DX_3_0*clhs39 + DN_DX_3_1*clhs40 + DN_DX_3_2*clhs41;
const double clhs83 =             clhs66*clhs82;
const double clhs84 =             clhs0*clhs1*clhs59*tau[1];
const double clhs85 =             0.1381966*clhs84;
const double clhs86 =             clhs20*clhs59*tau[1];
const double clhs87 =             0.1381966*clhs86;
const double clhs88 =             clhs0*clhs1*clhs56*tau[0];
const double clhs89 =             clhs0*clhs1*clhs62*tau[2];
const double clhs90 =             0.1381966*clhs89;
const double clhs91 =             clhs0*clhs1*clhs65*tau[3];
const double clhs92 =             0.1381966*clhs91;
const double clhs93 =             clhs20*clhs56*tau[0];
const double clhs94 =             clhs20*clhs62*tau[2];
const double clhs95 =             0.1381966*clhs94;
const double clhs96 =             clhs20*clhs65*tau[3];
const double clhs97 =             0.1381966*clhs96;
const double clhs98 =             0.1381966*clhs88;
const double clhs99 =             0.1381966*clhs93;
const double clhs100 =             clhs56*tau[0];
const double clhs101 =             clhs100*clhs68;
const double clhs102 =             clhs59*tau[1];
const double clhs103 =             clhs102*clhs70;
const double clhs104 =             clhs62*tau[2];
const double clhs105 =             clhs104*clhs72;
const double clhs106 =             clhs65*tau[3];
const double clhs107 =             clhs106*clhs74;
const double clhs108 =             clhs100*clhs76;
const double clhs109 =             clhs102*clhs78;
const double clhs110 =             clhs104*clhs80;
const double clhs111 =             clhs106*clhs82;
const double clhs112 =             clhs0*clhs1*clhs72*tau[2];
const double clhs113 =             0.1381966*clhs112;
const double clhs114 =             clhs20*clhs72*tau[2];
const double clhs115 =             0.1381966*clhs114;
const double clhs116 =             clhs0*clhs1*clhs68*tau[0];
const double clhs117 =             clhs0*clhs1*clhs70*tau[1];
const double clhs118 =             0.1381966*clhs117;
const double clhs119 =             clhs0*clhs1*clhs74*tau[3];
const double clhs120 =             0.1381966*clhs119;
const double clhs121 =             clhs20*clhs68*tau[0];
const double clhs122 =             clhs20*clhs70*tau[1];
const double clhs123 =             0.1381966*clhs122;
const double clhs124 =             clhs20*clhs74*tau[3];
const double clhs125 =             0.1381966*clhs124;
const double clhs126 =             0.1381966*clhs116;
const double clhs127 =             0.1381966*clhs121;
const double clhs128 =             clhs55 + clhs68*clhs76*tau[0] + clhs70*clhs78*tau[1] + clhs72*clhs80*tau[2] + clhs74*clhs82*tau[3];
const double clhs129 =             clhs0*clhs1*clhs82*tau[3];
const double clhs130 =             0.1381966*clhs129;
const double clhs131 =             clhs20*clhs82*tau[3];
const double clhs132 =             0.1381966*clhs131;
const double clhs133 =             clhs0*clhs1*clhs76*tau[0];
const double clhs134 =             clhs0*clhs1*clhs78*tau[1];
const double clhs135 =             0.1381966*clhs134;
const double clhs136 =             clhs0*clhs1*clhs80*tau[2];
const double clhs137 =             0.1381966*clhs136;
const double clhs138 =             clhs20*clhs76*tau[0];
const double clhs139 =             clhs20*clhs78*tau[1];
const double clhs140 =             0.1381966*clhs139;
const double clhs141 =             clhs20*clhs80*tau[2];
const double clhs142 =             0.1381966*clhs141;
const double clhs143 =             0.1381966*clhs133;
const double clhs144 =             0.1381966*clhs138;
            lhs(0,0)=pow(clhs17, 2)*tau[1] + clhs19 + clhs22 + clhs24 + pow(clhs31, 2)*tau[0] + pow(clhs38, 2)*tau[2] + pow(clhs42, 2)*tau[3] + 0.5854102*clhs43 + clhs45 + clhs47 + 0.5854102*clhs48 + clhs50 + clhs52;
            lhs(0,1)=0.5854102*clhs18 + 0.5854102*clhs21 + clhs45 + clhs47 + clhs50 + clhs52 + clhs53 + clhs54 + clhs55 + clhs58 + clhs61 + clhs64 + clhs67;
            lhs(0,2)=clhs19 + clhs22 + 0.5854102*clhs44 + clhs47 + 0.5854102*clhs49 + clhs52 + clhs53 + clhs54 + clhs55 + clhs69 + clhs71 + clhs73 + clhs75;
            lhs(0,3)=clhs19 + clhs22 + clhs45 + 0.5854102*clhs46 + clhs50 + 0.5854102*clhs51 + clhs53 + clhs54 + clhs55 + clhs77 + clhs79 + clhs81 + clhs83;
            lhs(1,0)=clhs55 + clhs58 + clhs61 + clhs64 + clhs67 + clhs85 + clhs87 + 0.5854102*clhs88 + clhs90 + clhs92 + 0.5854102*clhs93 + clhs95 + clhs97;
            lhs(1,1)=clhs24 + pow(clhs56, 2)*tau[0] + pow(clhs59, 2)*tau[1] + pow(clhs62, 2)*tau[2] + pow(clhs65, 2)*tau[3] + 0.5854102*clhs84 + 0.5854102*clhs86 + clhs90 + clhs92 + clhs95 + clhs97 + clhs98 + clhs99;
            lhs(1,2)=clhs101 + clhs103 + clhs105 + clhs107 + clhs55 + clhs85 + clhs87 + 0.5854102*clhs89 + clhs92 + 0.5854102*clhs94 + clhs97 + clhs98 + clhs99;
            lhs(1,3)=clhs108 + clhs109 + clhs110 + clhs111 + clhs55 + clhs85 + clhs87 + clhs90 + 0.5854102*clhs91 + clhs95 + 0.5854102*clhs96 + clhs98 + clhs99;
            lhs(2,0)=clhs113 + clhs115 + 0.5854102*clhs116 + clhs118 + clhs120 + 0.5854102*clhs121 + clhs123 + clhs125 + clhs55 + clhs69 + clhs71 + clhs73 + clhs75;
            lhs(2,1)=clhs101 + clhs103 + clhs105 + clhs107 + clhs113 + clhs115 + 0.5854102*clhs117 + clhs120 + 0.5854102*clhs122 + clhs125 + clhs126 + clhs127 + clhs55;
            lhs(2,2)=0.5854102*clhs112 + 0.5854102*clhs114 + clhs118 + clhs120 + clhs123 + clhs125 + clhs126 + clhs127 + clhs24 + pow(clhs68, 2)*tau[0] + pow(clhs70, 2)*tau[1] + pow(clhs72, 2)*tau[2] + pow(clhs74, 2)*tau[3];
            lhs(2,3)=clhs113 + clhs115 + clhs118 + 0.5854102*clhs119 + clhs123 + 0.5854102*clhs124 + clhs126 + clhs127 + clhs128;
            lhs(3,0)=clhs130 + clhs132 + 0.5854102*clhs133 + clhs135 + clhs137 + 0.5854102*clhs138 + clhs140 + clhs142 + clhs55 + clhs77 + clhs79 + clhs81 + clhs83;
            lhs(3,1)=clhs108 + clhs109 + clhs110 + clhs111 + clhs130 + clhs132 + 0.5854102*clhs134 + clhs137 + 0.5854102*clhs139 + clhs142 + clhs143 + clhs144 + clhs55;
            lhs(3,2)=clhs128 + clhs130 + clhs132 + clhs135 + 0.5854102*clhs136 + clhs140 + 0.5854102*clhs141 + clhs143 + clhs144;
            lhs(3,3)=0.5854102*clhs129 + 0.5854102*clhs131 + clhs135 + clhs137 + clhs140 + clhs142 + clhs143 + clhs144 + clhs24 + pow(clhs76, 2)*tau[0] + pow(clhs78, 2)*tau[1] + pow(clhs80, 2)*tau[2] + pow(clhs82, 2)*tau[3];


    const double crhs0 =             phi_subscale_gauss[2]/qstau[2];
const double crhs1 =             0.1381966*crhs0;
const double crhs2 =             phi_subscale_gauss[3]/qstau[3];
const double crhs3 =             0.1381966*crhs2;
const double crhs4 =             1.0/RK_time_coefficient;
const double crhs5 =             1.0/delta_time;
const double crhs6 =             crhs4*crhs5;
const double crhs7 =             0.08090169924532*phi[3];
const double crhs8 =             -0.08090169924532*phi_old[3];
const double crhs9 =             0.01909830025156*phi[0];
const double crhs10 =             0.01909830025156*phi[1];
const double crhs11 =             0.01909830025156*phi[2];
const double crhs12 =             -0.01909830025156*phi_old[0];
const double crhs13 =             -0.01909830025156*phi_old[1];
const double crhs14 =             -0.01909830025156*phi_old[2];
const double crhs15 =             crhs6*(crhs10 + crhs11 + crhs12 + crhs13 + crhs14 + crhs7 + crhs8 + crhs9);
const double crhs16 =             0.08090169924532*phi[2];
const double crhs17 =             -0.08090169924532*phi_old[2];
const double crhs18 =             0.01909830025156*phi[3];
const double crhs19 =             -0.01909830025156*phi_old[3];
const double crhs20 =             crhs6*(crhs10 + crhs12 + crhs13 + crhs16 + crhs17 + crhs18 + crhs19 + crhs9);
const double crhs21 =             crhs1 + crhs15 + crhs20 + crhs3 - 0.19999999899376*prj[2] - 0.19999999899376*prj[3];
const double crhs22 =             -0.19999999899376*prj[1];
const double crhs23 =             phi_subscale_gauss[0]/qstau[0];
const double crhs24 =             phi_subscale_gauss[1]/qstau[1];
const double crhs25 =             0.1381966*crhs24;
const double crhs26 =             0.08090169924532*phi[1];
const double crhs27 =             -0.08090169924532*phi_old[1];
const double crhs28 =             crhs6*(crhs11 + crhs12 + crhs14 + crhs18 + crhs19 + crhs26 + crhs27 + crhs9);
const double crhs29 =             crhs16 + crhs17 + crhs7 + crhs8;
const double crhs30 =             crhs5*phi_subscale_gauss[0]*tau[0];
const double crhs31 =             0.1381966*v(1,0);
const double crhs32 =             0.1381966*v(2,0);
const double crhs33 =             0.1381966*v(3,0);
const double crhs34 =             crhs32 + crhs33;
const double crhs35 =             crhs31 + crhs34 + 0.5854102*v(0,0);
const double crhs36 =             0.1381966*v(1,1);
const double crhs37 =             0.1381966*v(2,1);
const double crhs38 =             0.1381966*v(3,1);
const double crhs39 =             crhs37 + crhs38;
const double crhs40 =             crhs36 + crhs39 + 0.5854102*v(0,1);
const double crhs41 =             0.1381966*v(1,2);
const double crhs42 =             0.1381966*v(2,2);
const double crhs43 =             0.1381966*v(3,2);
const double crhs44 =             crhs42 + crhs43;
const double crhs45 =             crhs41 + crhs44 + 0.5854102*v(0,2);
const double crhs46 =             DN_DX_0_0*crhs35 + DN_DX_0_1*crhs40 + DN_DX_0_2*crhs45;
const double crhs47 =             crhs5*phi_subscale_gauss[1]*tau[1];
const double crhs48 =             0.1381966*v(0,0);
const double crhs49 =             crhs34 + crhs48 + 0.5854102*v(1,0);
const double crhs50 =             0.1381966*v(0,1);
const double crhs51 =             crhs39 + crhs50 + 0.5854102*v(1,1);
const double crhs52 =             0.1381966*v(0,2);
const double crhs53 =             crhs44 + crhs52 + 0.5854102*v(1,2);
const double crhs54 =             DN_DX_0_0*crhs49 + DN_DX_0_1*crhs51 + DN_DX_0_2*crhs53;
const double crhs55 =             crhs5*phi_subscale_gauss[2]*tau[2];
const double crhs56 =             crhs31 + crhs48;
const double crhs57 =             crhs33 + crhs56 + 0.5854102*v(2,0);
const double crhs58 =             crhs36 + crhs50;
const double crhs59 =             crhs38 + crhs58 + 0.5854102*v(2,1);
const double crhs60 =             crhs41 + crhs52;
const double crhs61 =             crhs43 + crhs60 + 0.5854102*v(2,2);
const double crhs62 =             DN_DX_0_0*crhs57 + DN_DX_0_1*crhs59 + DN_DX_0_2*crhs61;
const double crhs63 =             crhs5*phi_subscale_gauss[3]*tau[3];
const double crhs64 =             crhs32 + crhs56 + 0.5854102*v(3,0);
const double crhs65 =             crhs37 + crhs58 + 0.5854102*v(3,1);
const double crhs66 =             crhs42 + crhs60 + 0.5854102*v(3,2);
const double crhs67 =             DN_DX_0_0*crhs64 + DN_DX_0_1*crhs65 + DN_DX_0_2*crhs66;
const double crhs68 =             0.1381966*f[1];
const double crhs69 =             0.1381966*f[2];
const double crhs70 =             0.1381966*f[3];
const double crhs71 =             crhs69 + crhs70;
const double crhs72 =             crhs68 + crhs71 + 0.5854102*f[0];
const double crhs73 =             crhs46*tau[0];
const double crhs74 =             0.1381966*prj[1];
const double crhs75 =             0.1381966*prj[2];
const double crhs76 =             0.1381966*prj[3];
const double crhs77 =             crhs75 + crhs76;
const double crhs78 =             crhs74 + crhs77 + 0.5854102*prj[0];
const double crhs79 =             0.1381966*f[0];
const double crhs80 =             crhs71 + crhs79 + 0.5854102*f[1];
const double crhs81 =             crhs54*tau[1];
const double crhs82 =             0.1381966*prj[0];
const double crhs83 =             crhs77 + crhs82 + 0.5854102*prj[1];
const double crhs84 =             crhs68 + crhs79;
const double crhs85 =             crhs70 + crhs84 + 0.5854102*f[2];
const double crhs86 =             crhs62*tau[2];
const double crhs87 =             crhs74 + crhs82;
const double crhs88 =             crhs76 + crhs87 + 0.5854102*prj[2];
const double crhs89 =             crhs69 + crhs84 + 0.5854102*f[3];
const double crhs90 =             crhs67*tau[3];
const double crhs91 =             crhs75 + crhs87 + 0.5854102*prj[3];
const double crhs92 =             -0.1381966*phi_old[2];
const double crhs93 =             -0.1381966*phi_old[3];
const double crhs94 =             crhs92 + crhs93;
const double crhs95 =             0.5854102*phi[0];
const double crhs96 =             0.1381966*phi[1];
const double crhs97 =             0.1381966*phi[2];
const double crhs98 =             0.1381966*phi[3];
const double crhs99 =             -0.1381966*phi_old[1];
const double crhs100 =             crhs4*crhs5*(crhs94 + crhs95 + crhs96 + crhs97 + crhs98 + crhs99 - 0.5854102*phi_old[0]);
const double crhs101 =             0.1381966*phi[0];
const double crhs102 =             0.5854102*phi[1];
const double crhs103 =             -0.1381966*phi_old[0];
const double crhs104 =             crhs4*crhs5*(crhs101 + crhs102 + crhs103 + crhs94 + crhs97 + crhs98 - 0.5854102*phi_old[1]);
const double crhs105 =             crhs103 + crhs99;
const double crhs106 =             0.5854102*phi[2];
const double crhs107 =             crhs4*crhs5*(crhs101 + crhs105 + crhs106 + crhs93 + crhs96 + crhs98 - 0.5854102*phi_old[2]);
const double crhs108 =             0.5854102*phi[3];
const double crhs109 =             crhs4*crhs5*(crhs101 + crhs105 + crhs108 + crhs92 + crhs96 + crhs97 - 0.5854102*phi_old[3]);
const double crhs110 =             crhs97 + crhs98;
const double crhs111 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs112 =             crhs111*(crhs110 + crhs95 + crhs96);
const double crhs113 =             crhs111*(crhs101 + crhs102 + crhs110);
const double crhs114 =             crhs101 + crhs96;
const double crhs115 =             crhs111*(crhs106 + crhs114 + crhs98);
const double crhs116 =             crhs111*(crhs108 + crhs114 + crhs97);
const double crhs117 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs118 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs119 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs120 =             crhs117*crhs35 + crhs118*crhs40 + crhs119*crhs45;
const double crhs121 =             crhs117*crhs49 + crhs118*crhs51 + crhs119*crhs53;
const double crhs122 =             crhs117*crhs57 + crhs118*crhs59 + crhs119*crhs61;
const double crhs123 =             crhs117*crhs64 + crhs118*crhs65 + crhs119*crhs66;
const double crhs124 =             -0.19999999899376*prj[0];
const double crhs125 =             0.1381966*crhs23;
const double crhs126 =             0.08090169924532*phi[0];
const double crhs127 =             -0.08090169924532*phi_old[0];
const double crhs128 =             crhs6*(crhs10 + crhs11 + crhs126 + crhs127 + crhs13 + crhs14 + crhs18 + crhs19);
const double crhs129 =             DN_DX_1_0*crhs35 + DN_DX_1_1*crhs40 + DN_DX_1_2*crhs45;
const double crhs130 =             DN_DX_1_0*crhs49 + DN_DX_1_1*crhs51 + DN_DX_1_2*crhs53;
const double crhs131 =             DN_DX_1_0*crhs57 + DN_DX_1_1*crhs59 + DN_DX_1_2*crhs61;
const double crhs132 =             DN_DX_1_0*crhs64 + DN_DX_1_1*crhs65 + DN_DX_1_2*crhs66;
const double crhs133 =             crhs129*tau[0];
const double crhs134 =             crhs130*tau[1];
const double crhs135 =             crhs131*tau[2];
const double crhs136 =             crhs132*tau[3];
const double crhs137 =             crhs126 + crhs127 + crhs26 + crhs27;
const double crhs138 =             DN_DX_2_0*crhs35 + DN_DX_2_1*crhs40 + DN_DX_2_2*crhs45;
const double crhs139 =             DN_DX_2_0*crhs49 + DN_DX_2_1*crhs51 + DN_DX_2_2*crhs53;
const double crhs140 =             DN_DX_2_0*crhs57 + DN_DX_2_1*crhs59 + DN_DX_2_2*crhs61;
const double crhs141 =             DN_DX_2_0*crhs64 + DN_DX_2_1*crhs65 + DN_DX_2_2*crhs66;
const double crhs142 =             crhs138*tau[0];
const double crhs143 =             crhs139*tau[1];
const double crhs144 =             crhs140*tau[2];
const double crhs145 =             crhs141*tau[3];
const double crhs146 =             DN_DX_3_0*crhs35 + DN_DX_3_1*crhs40 + DN_DX_3_2*crhs45;
const double crhs147 =             DN_DX_3_0*crhs49 + DN_DX_3_1*crhs51 + DN_DX_3_2*crhs53;
const double crhs148 =             DN_DX_3_0*crhs57 + DN_DX_3_1*crhs59 + DN_DX_3_2*crhs61;
const double crhs149 =             DN_DX_3_0*crhs64 + DN_DX_3_1*crhs65 + DN_DX_3_2*crhs66;
const double crhs150 =             crhs146*tau[0];
const double crhs151 =             crhs147*tau[1];
const double crhs152 =             crhs148*tau[2];
const double crhs153 =             crhs149*tau[3];
            rhs[0]=-crhs100*crhs73 - crhs104*crhs81 - crhs107*crhs86 - crhs109*crhs90 - crhs112*crhs73 - crhs113*crhs81 - crhs115*crhs86 - crhs116*crhs90 - crhs120*crhs73 - crhs121*crhs81 - crhs122*crhs86 - crhs123*crhs90 + crhs21 + crhs22 + 0.5854102*crhs23 + crhs25 + crhs28 + crhs30*crhs46 + crhs47*crhs54 + crhs55*crhs62 + crhs6*(crhs26 + crhs27 + crhs29 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs63*crhs67 + crhs72*crhs73 + crhs73*crhs78 + crhs80*crhs81 + crhs81*crhs83 + crhs85*crhs86 + crhs86*crhs88 + crhs89*crhs90 + crhs90*crhs91 - 0.40000000301872*prj[0];
            rhs[1]=-crhs100*crhs133 - crhs104*crhs134 - crhs107*crhs135 - crhs109*crhs136 - crhs112*crhs133 - crhs113*crhs134 - crhs115*crhs135 - crhs116*crhs136 - crhs120*crhs133 - crhs121*crhs134 - crhs122*crhs135 - crhs123*crhs136 + crhs124 + crhs125 + crhs128 + crhs129*crhs30 + crhs130*crhs47 + crhs131*crhs55 + crhs132*crhs63 + crhs133*crhs72 + crhs133*crhs78 + crhs134*crhs80 + crhs134*crhs83 + crhs135*crhs85 + crhs135*crhs88 + crhs136*crhs89 + crhs136*crhs91 + crhs21 + 0.5854102*crhs24 + crhs6*(crhs126 + crhs127 + crhs29 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) - 0.40000000301872*prj[1];
            rhs[2]=0.5854102*crhs0 - crhs100*crhs142 - crhs104*crhs143 - crhs107*crhs144 - crhs109*crhs145 - crhs112*crhs142 - crhs113*crhs143 - crhs115*crhs144 - crhs116*crhs145 - crhs120*crhs142 - crhs121*crhs143 - crhs122*crhs144 - crhs123*crhs145 + crhs124 + crhs125 + crhs128 + crhs138*crhs30 + crhs139*crhs47 + crhs140*crhs55 + crhs141*crhs63 + crhs142*crhs72 + crhs142*crhs78 + crhs143*crhs80 + crhs143*crhs83 + crhs144*crhs85 + crhs144*crhs88 + crhs145*crhs89 + crhs145*crhs91 + crhs15 + crhs22 + crhs25 + crhs28 + crhs3 + crhs6*(crhs137 + crhs7 + crhs8 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) - 0.40000000301872*prj[2] - 0.19999999899376*prj[3];
            rhs[3]=crhs1 - crhs100*crhs150 - crhs104*crhs151 - crhs107*crhs152 - crhs109*crhs153 - crhs112*crhs150 - crhs113*crhs151 - crhs115*crhs152 - crhs116*crhs153 - crhs120*crhs150 - crhs121*crhs151 - crhs122*crhs152 - crhs123*crhs153 + crhs124 + crhs125 + crhs128 + crhs146*crhs30 + crhs147*crhs47 + crhs148*crhs55 + crhs149*crhs63 + crhs150*crhs72 + crhs150*crhs78 + crhs151*crhs80 + crhs151*crhs83 + crhs152*crhs85 + crhs152*crhs88 + crhs153*crhs89 + crhs153*crhs91 + 0.5854102*crhs2 + crhs20 + crhs22 + crhs25 + crhs28 + crhs6*(crhs137 + crhs16 + crhs17 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) - 0.19999999899376*prj[2] - 0.40000000301872*prj[3];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rLeftHandSideMatrix) += lhs * rVariables.volume/local_size;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::CalculateOrthogonalSubgridScaleSystemInternal(
    ElementVariables& rVariables,
    VectorType& rRightHandSideVector)
{
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

    const double crhs0 =             -0.25*f[1];
const double crhs1 =             1.0/delta_time;
const double crhs2 =             crhs1*phi_subscale_gauss[1];
const double crhs3 =             -0.166666666666667*crhs2;
const double crhs4 =             0.166666666666667*v(0,0);
const double crhs5 =             0.166666666666667*v(2,0);
const double crhs6 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs7 =             crhs6*(crhs4 + crhs5 + 0.666666666666667*v(1,0));
const double crhs8 =             0.166666666666667*crhs7;
const double crhs9 =             0.166666666666667*v(0,1);
const double crhs10 =             0.166666666666667*v(2,1);
const double crhs11 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs12 =             crhs11*(crhs10 + crhs9 + 0.666666666666667*v(1,1));
const double crhs13 =             0.166666666666667*crhs12;
const double crhs14 =             crhs1/RK_time_coefficient;
const double crhs15 =             0.111111111111111*phi[1];
const double crhs16 =             -0.111111111111111*phi_old[1];
const double crhs17 =             crhs15 + crhs16;
const double crhs18 =             0.0277777777777778*phi[0];
const double crhs19 =             0.0277777777777778*phi[2];
const double crhs20 =             -0.0277777777777778*phi_old[0];
const double crhs21 =             -0.0277777777777778*phi_old[2];
const double crhs22 =             crhs14*(crhs17 + crhs18 + crhs19 + crhs20 + crhs21);
const double crhs23 =             0.166666666666667*phi[0];
const double crhs24 =             0.166666666666667*phi[2];
const double crhs25 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs26 =             crhs25*(crhs23 + crhs24 + 0.666666666666667*phi[1]);
const double crhs27 =             0.166666666666667*crhs26;
const double crhs28 =             -0.25*f[2];
const double crhs29 =             crhs1*phi_subscale_gauss[0];
const double crhs30 =             crhs1*phi_subscale_gauss[2];
const double crhs31 =             -0.166666666666667*crhs30;
const double crhs32 =             3*crhs6*k;
const double crhs33 =             3*crhs11*k;
const double crhs34 =             0.166666666666667*v(1,0);
const double crhs35 =             crhs6*(crhs34 + crhs4 + 0.666666666666667*v(2,0));
const double crhs36 =             0.166666666666667*crhs35;
const double crhs37 =             crhs6*(crhs34 + crhs5 + 0.666666666666667*v(0,0));
const double crhs38 =             0.166666666666667*v(1,1);
const double crhs39 =             crhs11*(crhs38 + crhs9 + 0.666666666666667*v(2,1));
const double crhs40 =             0.166666666666667*crhs39;
const double crhs41 =             crhs11*(crhs10 + crhs38 + 0.666666666666667*v(0,1));
const double crhs42 =             0.111111111111111*phi[2];
const double crhs43 =             -0.111111111111111*phi_old[2];
const double crhs44 =             0.0277777777777778*phi[1];
const double crhs45 =             -0.0277777777777778*phi_old[1];
const double crhs46 =             crhs14*(crhs18 + crhs20 + crhs42 + crhs43 + crhs44 + crhs45);
const double crhs47 =             0.166666666666667*phi[1];
const double crhs48 =             crhs25*(crhs23 + crhs47 + 0.666666666666667*phi[2]);
const double crhs49 =             0.166666666666667*crhs48;
const double crhs50 =             crhs25*(crhs24 + crhs47 + 0.666666666666667*phi[0]);
const double crhs51 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs52 =             crhs14*(crhs19 + crhs21 + crhs44 + crhs45 + crhs51) - 0.166666666666667*crhs29 + 0.166666666666667*crhs37 + 0.166666666666667*crhs41 + 0.166666666666667*crhs50 - 0.25*f[0];
            rhs[0]=DN_DX_0_0*crhs32 + DN_DX_0_1*crhs33 + crhs0 + crhs13 + crhs14*(crhs17 + crhs42 + crhs43 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs22 + crhs27 + crhs28 - 0.666666666666667*crhs29 + crhs3 + crhs31 + crhs36 + 0.666666666666667*crhs37 + crhs40 + 0.666666666666667*crhs41 + crhs46 + crhs49 + 0.666666666666667*crhs50 + crhs8 - 0.5*f[0];
            rhs[1]=DN_DX_1_0*crhs32 + DN_DX_1_1*crhs33 + 0.666666666666667*crhs12 + crhs14*(crhs42 + crhs43 + crhs51 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) - 0.666666666666667*crhs2 + 0.666666666666667*crhs26 + crhs28 + crhs31 + crhs36 + crhs40 + crhs46 + crhs49 + crhs52 + 0.666666666666667*crhs7 - 0.5*f[1];
            rhs[2]=DN_DX_2_0*crhs32 + DN_DX_2_1*crhs33 + crhs0 + crhs13 + crhs14*(crhs15 + crhs16 + crhs51 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs22 + crhs27 + crhs3 - 0.666666666666667*crhs30 + 0.666666666666667*crhs35 + 0.666666666666667*crhs39 + 0.666666666666667*crhs48 + crhs52 + crhs8 - 0.5*f[2];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::CalculateOrthogonalSubgridScaleSystemInternal(
    ElementVariables& rVariables,
    VectorType& rRightHandSideVector)
{
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

    const double crhs0 =             -0.19999999899376*f[1];
const double crhs1 =             1.0/delta_time;
const double crhs2 =             crhs1*phi_subscale_gauss[1];
const double crhs3 =             -0.1381966*crhs2;
const double crhs4 =             0.1381966*v(0,0);
const double crhs5 =             0.1381966*v(2,0);
const double crhs6 =             0.1381966*v(3,0);
const double crhs7 =             crhs5 + crhs6;
const double crhs8 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs9 =             crhs8*(crhs4 + crhs7 + 0.5854102*v(1,0));
const double crhs10 =             0.1381966*crhs9;
const double crhs11 =             0.1381966*v(0,1);
const double crhs12 =             0.1381966*v(2,1);
const double crhs13 =             0.1381966*v(3,1);
const double crhs14 =             crhs12 + crhs13;
const double crhs15 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs16 =             crhs15*(crhs11 + crhs14 + 0.5854102*v(1,1));
const double crhs17 =             0.1381966*crhs16;
const double crhs18 =             0.1381966*v(0,2);
const double crhs19 =             0.1381966*v(2,2);
const double crhs20 =             0.1381966*v(3,2);
const double crhs21 =             crhs19 + crhs20;
const double crhs22 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs23 =             crhs22*(crhs18 + crhs21 + 0.5854102*v(1,2));
const double crhs24 =             0.1381966*crhs23;
const double crhs25 =             crhs1/RK_time_coefficient;
const double crhs26 =             0.08090169924532*phi[1];
const double crhs27 =             -0.08090169924532*phi_old[1];
const double crhs28 =             0.01909830025156*phi[0];
const double crhs29 =             0.01909830025156*phi[2];
const double crhs30 =             0.01909830025156*phi[3];
const double crhs31 =             -0.01909830025156*phi_old[0];
const double crhs32 =             -0.01909830025156*phi_old[2];
const double crhs33 =             -0.01909830025156*phi_old[3];
const double crhs34 =             crhs25*(crhs26 + crhs27 + crhs28 + crhs29 + crhs30 + crhs31 + crhs32 + crhs33);
const double crhs35 =             0.1381966*phi[0];
const double crhs36 =             0.1381966*phi[2];
const double crhs37 =             0.1381966*phi[3];
const double crhs38 =             crhs36 + crhs37;
const double crhs39 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs40 =             crhs39*(crhs35 + crhs38 + 0.5854102*phi[1]);
const double crhs41 =             0.1381966*crhs40;
const double crhs42 =             -0.19999999899376*f[2];
const double crhs43 =             -0.19999999899376*f[3];
const double crhs44 =             crhs1*phi_subscale_gauss[0];
const double crhs45 =             crhs1*phi_subscale_gauss[2];
const double crhs46 =             -0.1381966*crhs45;
const double crhs47 =             crhs1*phi_subscale_gauss[3];
const double crhs48 =             -0.1381966*crhs47;
const double crhs49 =             4*crhs8*k;
const double crhs50 =             4*crhs15*k;
const double crhs51 =             4*crhs22*k;
const double crhs52 =             0.1381966*v(1,0);
const double crhs53 =             crhs4 + crhs52;
const double crhs54 =             crhs8*(crhs5 + crhs53 + 0.5854102*v(3,0));
const double crhs55 =             0.1381966*crhs54;
const double crhs56 =             crhs8*(crhs53 + crhs6 + 0.5854102*v(2,0));
const double crhs57 =             0.1381966*crhs56;
const double crhs58 =             crhs8*(crhs52 + crhs7 + 0.5854102*v(0,0));
const double crhs59 =             0.1381966*v(1,1);
const double crhs60 =             crhs11 + crhs59;
const double crhs61 =             crhs15*(crhs12 + crhs60 + 0.5854102*v(3,1));
const double crhs62 =             0.1381966*crhs61;
const double crhs63 =             crhs15*(crhs13 + crhs60 + 0.5854102*v(2,1));
const double crhs64 =             0.1381966*crhs63;
const double crhs65 =             crhs15*(crhs14 + crhs59 + 0.5854102*v(0,1));
const double crhs66 =             0.1381966*v(1,2);
const double crhs67 =             crhs18 + crhs66;
const double crhs68 =             crhs22*(crhs19 + crhs67 + 0.5854102*v(3,2));
const double crhs69 =             0.1381966*crhs68;
const double crhs70 =             crhs22*(crhs20 + crhs67 + 0.5854102*v(2,2));
const double crhs71 =             0.1381966*crhs70;
const double crhs72 =             crhs22*(crhs21 + crhs66 + 0.5854102*v(0,2));
const double crhs73 =             0.08090169924532*phi[3];
const double crhs74 =             -0.08090169924532*phi_old[3];
const double crhs75 =             0.01909830025156*phi[1];
const double crhs76 =             -0.01909830025156*phi_old[1];
const double crhs77 =             crhs25*(crhs28 + crhs29 + crhs31 + crhs32 + crhs73 + crhs74 + crhs75 + crhs76);
const double crhs78 =             0.08090169924532*phi[2];
const double crhs79 =             -0.08090169924532*phi_old[2];
const double crhs80 =             crhs25*(crhs28 + crhs30 + crhs31 + crhs33 + crhs75 + crhs76 + crhs78 + crhs79);
const double crhs81 =             crhs73 + crhs74 + crhs78 + crhs79;
const double crhs82 =             0.1381966*phi[1];
const double crhs83 =             crhs35 + crhs82;
const double crhs84 =             crhs39*(crhs36 + crhs83 + 0.5854102*phi[3]);
const double crhs85 =             0.1381966*crhs84;
const double crhs86 =             crhs39*(crhs37 + crhs83 + 0.5854102*phi[2]);
const double crhs87 =             0.1381966*crhs86;
const double crhs88 =             crhs39*(crhs38 + crhs82 + 0.5854102*phi[0]);
const double crhs89 =             -0.19999999899376*f[0];
const double crhs90 =             -0.1381966*crhs44;
const double crhs91 =             0.1381966*crhs58;
const double crhs92 =             0.1381966*crhs65;
const double crhs93 =             0.1381966*crhs72;
const double crhs94 =             0.08090169924532*phi[0];
const double crhs95 =             -0.08090169924532*phi_old[0];
const double crhs96 =             crhs25*(crhs29 + crhs30 + crhs32 + crhs33 + crhs75 + crhs76 + crhs94 + crhs95);
const double crhs97 =             0.1381966*crhs88;
const double crhs98 =             crhs0 + crhs10 + crhs17 + crhs24 + crhs3 + crhs34 + crhs41 + crhs89 + crhs90 + crhs91 + crhs92 + crhs93 + crhs96 + crhs97;
const double crhs99 =             crhs26 + crhs27 + crhs94 + crhs95;
            rhs[0]=DN_DX_0_0*crhs49 + DN_DX_0_1*crhs50 + DN_DX_0_2*crhs51 + crhs0 + crhs10 + crhs17 + crhs24 + crhs25*(crhs26 + crhs27 + crhs81 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs3 + crhs34 + crhs41 + crhs42 + crhs43 - 0.5854102*crhs44 + crhs46 + crhs48 + crhs55 + crhs57 + 0.5854102*crhs58 + crhs62 + crhs64 + 0.5854102*crhs65 + crhs69 + crhs71 + 0.5854102*crhs72 + crhs77 + crhs80 + crhs85 + crhs87 + 0.5854102*crhs88 - 0.40000000301872*f[0];
            rhs[1]=DN_DX_1_0*crhs49 + DN_DX_1_1*crhs50 + DN_DX_1_2*crhs51 + 0.5854102*crhs16 - 0.5854102*crhs2 + 0.5854102*crhs23 + crhs25*(crhs81 + crhs94 + crhs95 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) + 0.5854102*crhs40 + crhs42 + crhs43 + crhs46 + crhs48 + crhs55 + crhs57 + crhs62 + crhs64 + crhs69 + crhs71 + crhs77 + crhs80 + crhs85 + crhs87 + crhs89 + 0.5854102*crhs9 + crhs90 + crhs91 + crhs92 + crhs93 + crhs96 + crhs97 - 0.40000000301872*f[1];
            rhs[2]=DN_DX_2_0*crhs49 + DN_DX_2_1*crhs50 + DN_DX_2_2*crhs51 + crhs25*(crhs73 + crhs74 + crhs99 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) - 0.5854102*crhs45 + crhs48 + crhs55 + 0.5854102*crhs56 + crhs62 + 0.5854102*crhs63 + crhs69 + 0.5854102*crhs70 + crhs77 + crhs85 + 0.5854102*crhs86 + crhs98 - 0.40000000301872*f[2] - 0.19999999899376*f[3];
            rhs[3]=DN_DX_3_0*crhs49 + DN_DX_3_1*crhs50 + DN_DX_3_2*crhs51 + crhs25*(crhs78 + crhs79 + crhs99 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) + crhs46 - 0.5854102*crhs47 + 0.5854102*crhs54 + crhs57 + 0.5854102*crhs61 + crhs64 + 0.5854102*crhs68 + crhs71 + crhs80 + 0.5854102*crhs84 + crhs87 + crhs98 - 0.19999999899376*f[2] - 0.40000000301872*f[3];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) += rhs * rVariables.volume/local_size;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
    // Retrieve element variables
    const auto& N = rVariables.N;
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
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
    const auto& N = rVariables.N;
    const auto& k = rVariables.diffusivity;
    const auto& f = rVariables.forcing;
    const auto& phi = rVariables.unknown;
    const auto& phi_old = rVariables.unknown_old;
    const auto& delta_time = rVariables.delta_time;
    const auto& RK_time_coefficient = rVariables.RK_time_coefficient;
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
    for(unsigned int g = 0; g<TNumNodes; g++)
    {
	const auto& N = row(rVariables.N_gausspoint,g);
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
