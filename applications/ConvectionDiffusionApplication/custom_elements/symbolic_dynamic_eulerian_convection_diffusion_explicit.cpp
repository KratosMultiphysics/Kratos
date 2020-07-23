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
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    Matrix LeftHandSide;
    this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    VectorType RightHandSide;
    this->CalculateLocalSystem(rLeftHandSideMatrix,RightHandSide,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    ProcessInfo &rCurrentProcessInfo)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
    // Resize and intialize output
    if (mUnknownSubScale.size() != TNumNodes) // number integration points = number nodes
        mUnknownSubScale.resize(TNumNodes, false);
    mUnknownSubScale = ZeroVector(TNumNodes);
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
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    VectorType rhs_oss;
    this->CalculateOrthogonalSubgridScaleSystem(rhs_oss,rCurrentProcessInfo);
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateOrthogonalSubgridScaleSystem(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
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
const double clhs2 =             clhs0*clhs1*tau[1];
const double clhs3 =             -0.0277777777777778*clhs2;
const double clhs4 =             DN_DX_0_0*v(0,0);
const double clhs5 =             DN_DX_0_1*v(0,1);
const double clhs6 =             DN_DX_1_0*v(1,0);
const double clhs7 =             DN_DX_1_1*v(1,1);
const double clhs8 =             DN_DX_2_0*v(2,0);
const double clhs9 =             DN_DX_2_1*v(2,1);
const double clhs10 =             clhs4 + clhs5 + clhs6 + clhs7 + clhs8 + clhs9;
const double clhs11 =             clhs10*tau[1];
const double clhs12 =             -0.0277777777777778*clhs11;
const double clhs13 =             0.166666666666667*v(0,0);
const double clhs14 =             0.166666666666667*v(2,0);
const double clhs15 =             clhs13 + clhs14 + 0.666666666666667*v(1,0);
const double clhs16 =             DN_DX_0_0*clhs15;
const double clhs17 =             0.166666666666667*v(0,1);
const double clhs18 =             0.166666666666667*v(2,1);
const double clhs19 =             clhs17 + clhs18 + 0.666666666666667*v(1,1);
const double clhs20 =             DN_DX_0_1*clhs19;
const double clhs21 =             clhs16 + clhs20;
const double clhs22 =             clhs0*clhs1*clhs21*tau[1];
const double clhs23 =             0.166666666666667*clhs22;
const double clhs24 =             clhs10*clhs21*tau[1];
const double clhs25 =             0.166666666666667*clhs24;
const double clhs26 =             0.5*clhs4;
const double clhs27 =             0.5*clhs5;
const double clhs28 =             0.5*clhs6;
const double clhs29 =             0.5*clhs7;
const double clhs30 =             0.5*clhs8;
const double clhs31 =             0.5*clhs9;
const double clhs32 =             pow(DN_DX_0_0, 2);
const double clhs33 =             3*k;
const double clhs34 =             pow(DN_DX_0_1, 2);
const double clhs35 =             clhs32*k;
const double clhs36 =             clhs34*k;
const double clhs37 =             clhs0*clhs1*tau[0];
const double clhs38 =             clhs0*clhs1*tau[2];
const double clhs39 =             -0.0277777777777778*clhs38;
const double clhs40 =             0.166666666666667*v(1,0);
const double clhs41 =             clhs13 + clhs40 + 0.666666666666667*v(2,0);
const double clhs42 =             DN_DX_0_0*clhs41;
const double clhs43 =             0.166666666666667*clhs42;
const double clhs44 =             0.166666666666667*clhs16;
const double clhs45 =             clhs14 + clhs40 + 0.666666666666667*v(0,0);
const double clhs46 =             DN_DX_0_0*clhs45;
const double clhs47 =             0.166666666666667*v(1,1);
const double clhs48 =             clhs17 + clhs47 + 0.666666666666667*v(2,1);
const double clhs49 =             DN_DX_0_1*clhs48;
const double clhs50 =             0.166666666666667*clhs49;
const double clhs51 =             0.166666666666667*clhs20;
const double clhs52 =             clhs18 + clhs47 + 0.666666666666667*v(0,1);
const double clhs53 =             DN_DX_0_1*clhs52;
const double clhs54 =             clhs10*tau[0];
const double clhs55 =             clhs10*tau[2];
const double clhs56 =             -0.0277777777777778*clhs55;
const double clhs57 =             clhs46 + clhs53;
const double clhs58 =             clhs57*tau[0];
const double clhs59 =             clhs21*tau[1];
const double clhs60 =             -0.166666666666667*clhs59;
const double clhs61 =             clhs42 + clhs49;
const double clhs62 =             clhs61*tau[2];
const double clhs63 =             -0.166666666666667*clhs62;
const double clhs64 =             clhs0*clhs1*clhs57*tau[0];
const double clhs65 =             clhs0*clhs1*clhs61*tau[2];
const double clhs66 =             0.166666666666667*clhs65;
const double clhs67 =             clhs10*clhs57*tau[0];
const double clhs68 =             clhs10*clhs61*tau[2];
const double clhs69 =             0.166666666666667*clhs68;
const double clhs70 =             0.166666666666667*clhs64;
const double clhs71 =             0.166666666666667*clhs67;
const double clhs72 =             0.25*clhs4;
const double clhs73 =             0.25*clhs5;
const double clhs74 =             0.25*clhs6;
const double clhs75 =             0.25*clhs7;
const double clhs76 =             0.25*clhs8;
const double clhs77 =             0.25*clhs9;
const double clhs78 =             3*DN_DX_0_0*k;
const double clhs79 =             DN_DX_1_0*clhs78;
const double clhs80 =             3*DN_DX_0_1*k;
const double clhs81 =             DN_DX_1_1*clhs80;
const double clhs82 =             DN_DX_0_0*DN_DX_1_0*k;
const double clhs83 =             -clhs82*tau[0];
const double clhs84 =             -clhs82*tau[1];
const double clhs85 =             -clhs82*tau[2];
const double clhs86 =             DN_DX_0_1*DN_DX_1_1*k;
const double clhs87 =             -clhs86*tau[0];
const double clhs88 =             -clhs86*tau[1];
const double clhs89 =             -clhs86*tau[2];
const double clhs90 =             -0.111111111111111*clhs37;
const double clhs91 =             -0.111111111111111*clhs2;
const double clhs92 =             DN_DX_1_0*clhs41;
const double clhs93 =             0.166666666666667*clhs92;
const double clhs94 =             DN_DX_1_0*clhs15;
const double clhs95 =             0.166666666666667*clhs94;
const double clhs96 =             DN_DX_1_0*clhs45;
const double clhs97 =             DN_DX_1_1*clhs48;
const double clhs98 =             0.166666666666667*clhs97;
const double clhs99 =             DN_DX_1_1*clhs19;
const double clhs100 =             0.166666666666667*clhs99;
const double clhs101 =             DN_DX_1_1*clhs52;
const double clhs102 =             -0.111111111111111*clhs54;
const double clhs103 =             -0.111111111111111*clhs11;
const double clhs104 =             clhs101 + clhs96;
const double clhs105 =             clhs104*tau[0];
const double clhs106 =             clhs94 + clhs99;
const double clhs107 =             clhs106*tau[1];
const double clhs108 =             -0.166666666666667*clhs107;
const double clhs109 =             clhs92 + clhs97;
const double clhs110 =             clhs109*tau[2];
const double clhs111 =             -0.166666666666667*clhs110;
const double clhs112 =             clhs104*clhs58;
const double clhs113 =             clhs106*clhs59;
const double clhs114 =             clhs109*clhs62;
const double clhs115 =             DN_DX_0_0*DN_DX_2_0*k;
const double clhs116 =             DN_DX_0_1*DN_DX_2_1*k;
const double clhs117 =             -0.111111111111111*clhs38;
const double clhs118 =             -0.111111111111111*clhs55;
const double clhs119 =             DN_DX_2_0*clhs45;
const double clhs120 =             DN_DX_2_1*clhs52;
const double clhs121 =             clhs119 + clhs120;
const double clhs122 =             DN_DX_2_0*clhs15;
const double clhs123 =             DN_DX_2_1*clhs19;
const double clhs124 =             clhs122 + clhs123;
const double clhs125 =             DN_DX_2_0*clhs41;
const double clhs126 =             DN_DX_2_1*clhs48;
const double clhs127 =             clhs125 + clhs126;
const double clhs128 =             DN_DX_2_0*clhs78 + DN_DX_2_1*clhs80 + clhs102 - clhs115*tau[0] - clhs115*tau[1] - clhs115*tau[2] - clhs116*tau[0] - clhs116*tau[1] - clhs116*tau[2] + clhs117 + clhs118 + clhs121*clhs58 + clhs124*clhs59 + clhs127*clhs62 + clhs72 + clhs73 + clhs74 + clhs75 + clhs76 + clhs77 + clhs90;
const double clhs129 =             0.166666666666667*clhs125;
const double clhs130 =             0.166666666666667*clhs122;
const double clhs131 =             0.166666666666667*clhs126;
const double clhs132 =             0.166666666666667*clhs123;
const double clhs133 =             clhs121*tau[0];
const double clhs134 =             clhs124*tau[1];
const double clhs135 =             -0.166666666666667*clhs134;
const double clhs136 =             clhs127*tau[2];
const double clhs137 =             -0.166666666666667*clhs136;
const double clhs138 =             0.166666666666667*clhs46;
const double clhs139 =             0.166666666666667*clhs53;
const double clhs140 =             -0.166666666666667*clhs58;
const double clhs141 =             clhs0*clhs1*clhs104*tau[0];
const double clhs142 =             clhs0*clhs1*clhs106*tau[1];
const double clhs143 =             0.166666666666667*clhs142;
const double clhs144 =             clhs0*clhs1*clhs109*tau[2];
const double clhs145 =             0.166666666666667*clhs144;
const double clhs146 =             clhs10*clhs104*tau[0];
const double clhs147 =             clhs10*clhs106*tau[1];
const double clhs148 =             0.166666666666667*clhs147;
const double clhs149 =             clhs10*clhs109*tau[2];
const double clhs150 =             0.166666666666667*clhs149;
const double clhs151 =             pow(DN_DX_1_0, 2);
const double clhs152 =             pow(DN_DX_1_1, 2);
const double clhs153 =             clhs151*k;
const double clhs154 =             clhs152*k;
const double clhs155 =             -0.0277777777777778*clhs37;
const double clhs156 =             0.166666666666667*clhs96;
const double clhs157 =             0.166666666666667*clhs101;
const double clhs158 =             -0.0277777777777778*clhs54;
const double clhs159 =             -0.166666666666667*clhs105;
const double clhs160 =             0.166666666666667*clhs141;
const double clhs161 =             0.166666666666667*clhs146;
const double clhs162 =             DN_DX_1_0*DN_DX_2_0*clhs33;
const double clhs163 =             DN_DX_1_1*DN_DX_2_1*clhs33;
const double clhs164 =             DN_DX_1_0*DN_DX_2_0*k;
const double clhs165 =             -clhs164*tau[0];
const double clhs166 =             -clhs164*tau[1];
const double clhs167 =             -clhs164*tau[2];
const double clhs168 =             DN_DX_1_1*DN_DX_2_1*k;
const double clhs169 =             -clhs168*tau[0];
const double clhs170 =             -clhs168*tau[1];
const double clhs171 =             -clhs168*tau[2];
const double clhs172 =             0.166666666666667*clhs119;
const double clhs173 =             0.166666666666667*clhs120;
const double clhs174 =             -0.166666666666667*clhs133;
const double clhs175 =             clhs105*clhs121;
const double clhs176 =             clhs107*clhs124;
const double clhs177 =             clhs110*clhs127;
const double clhs178 =             clhs0*clhs1*clhs121*tau[0];
const double clhs179 =             clhs0*clhs1*clhs124*tau[1];
const double clhs180 =             0.166666666666667*clhs179;
const double clhs181 =             clhs0*clhs1*clhs127*tau[2];
const double clhs182 =             0.166666666666667*clhs181;
const double clhs183 =             clhs10*clhs121*tau[0];
const double clhs184 =             clhs10*clhs124*tau[1];
const double clhs185 =             0.166666666666667*clhs184;
const double clhs186 =             clhs10*clhs127*tau[2];
const double clhs187 =             0.166666666666667*clhs186;
const double clhs188 =             0.166666666666667*clhs178;
const double clhs189 =             0.166666666666667*clhs183;
const double clhs190 =             pow(DN_DX_2_0, 2);
const double clhs191 =             pow(DN_DX_2_1, 2);
const double clhs192 =             clhs190*k;
const double clhs193 =             clhs191*k;
            lhs(0,0)=clhs12 + pow(clhs21, 2)*tau[1] + clhs23 + clhs25 + clhs26 + clhs27 + clhs28 + clhs29 + clhs3 + clhs30 + clhs31 + clhs32*clhs33 + clhs33*clhs34 - clhs35*tau[0] - clhs35*tau[1] - clhs35*tau[2] - clhs36*tau[0] - clhs36*tau[1] - clhs36*tau[2] - 0.444444444444444*clhs37 + clhs39 + clhs43 + clhs44 + 0.666666666666667*clhs46 + clhs50 + clhs51 + 0.666666666666667*clhs53 - 0.444444444444444*clhs54 + clhs56 + pow(clhs57, 2)*tau[0] - 0.666666666666667*clhs58 + clhs60 + pow(clhs61, 2)*tau[2] + clhs63 + 0.666666666666667*clhs64 + clhs66 + 0.666666666666667*clhs67 + clhs69;
            lhs(0,1)=clhs100 + 0.666666666666667*clhs101 + clhs102 + clhs103 - 0.666666666666667*clhs105 + clhs108 + clhs111 + clhs112 + clhs113 + clhs114 + 0.666666666666667*clhs22 + 0.666666666666667*clhs24 + clhs39 + clhs56 + clhs66 + clhs69 + clhs70 + clhs71 + clhs72 + clhs73 + clhs74 + clhs75 + clhs76 + clhs77 + clhs79 + clhs81 + clhs83 + clhs84 + clhs85 + clhs87 + clhs88 + clhs89 + clhs90 + clhs91 + clhs93 + clhs95 + 0.666666666666667*clhs96 + clhs98;
            lhs(0,2)=0.666666666666667*clhs119 + clhs12 + 0.666666666666667*clhs120 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 - 0.666666666666667*clhs133 + clhs135 + clhs137 + clhs23 + clhs25 + clhs3 + 0.666666666666667*clhs65 + 0.666666666666667*clhs68 + clhs70 + clhs71;
            lhs(1,0)=clhs102 + clhs103 + clhs112 + clhs113 + clhs114 + clhs138 + clhs139 + clhs140 + 0.666666666666667*clhs141 + clhs143 + clhs145 + 0.666666666666667*clhs146 + clhs148 + clhs150 + 0.666666666666667*clhs16 + 0.666666666666667*clhs20 + clhs39 + clhs43 + clhs50 + clhs56 - 0.666666666666667*clhs59 + clhs63 + clhs72 + clhs73 + clhs74 + clhs75 + clhs76 + clhs77 + clhs79 + clhs81 + clhs83 + clhs84 + clhs85 + clhs87 + clhs88 + clhs89 + clhs90 + clhs91;
            lhs(1,1)=pow(clhs104, 2)*tau[0] + pow(clhs106, 2)*tau[1] - 0.666666666666667*clhs107 + pow(clhs109, 2)*tau[2] - 0.444444444444444*clhs11 + clhs111 + 0.666666666666667*clhs142 + clhs145 + 0.666666666666667*clhs147 + clhs150 + clhs151*clhs33 + clhs152*clhs33 - clhs153*tau[0] - clhs153*tau[1] - clhs153*tau[2] - clhs154*tau[0] - clhs154*tau[1] - clhs154*tau[2] + clhs155 + clhs156 + clhs157 + clhs158 + clhs159 + clhs160 + clhs161 - 0.444444444444444*clhs2 + clhs26 + clhs27 + clhs28 + clhs29 + clhs30 + clhs31 + clhs39 + clhs56 + clhs93 + 0.666666666666667*clhs94 + clhs98 + 0.666666666666667*clhs99;
            lhs(1,2)=clhs103 + clhs117 + clhs118 + 0.666666666666667*clhs122 + 0.666666666666667*clhs123 + clhs129 + clhs131 - 0.666666666666667*clhs134 + clhs137 + clhs143 + 0.666666666666667*clhs144 + clhs148 + 0.666666666666667*clhs149 + clhs155 + clhs158 + clhs160 + clhs161 + clhs162 + clhs163 + clhs165 + clhs166 + clhs167 + clhs169 + clhs170 + clhs171 + clhs172 + clhs173 + clhs174 + clhs175 + clhs176 + clhs177 + clhs72 + clhs73 + clhs74 + clhs75 + clhs76 + clhs77 + clhs91;
            lhs(2,0)=clhs12 + clhs128 + clhs138 + clhs139 + clhs140 + 0.666666666666667*clhs178 + clhs180 + clhs182 + 0.666666666666667*clhs183 + clhs185 + clhs187 + clhs3 + 0.666666666666667*clhs42 + clhs44 + 0.666666666666667*clhs49 + clhs51 + clhs60 - 0.666666666666667*clhs62;
            lhs(2,1)=clhs100 + clhs103 + clhs108 - 0.666666666666667*clhs110 + clhs117 + clhs118 + clhs155 + clhs156 + clhs157 + clhs158 + clhs159 + clhs162 + clhs163 + clhs165 + clhs166 + clhs167 + clhs169 + clhs170 + clhs171 + clhs175 + clhs176 + clhs177 + 0.666666666666667*clhs179 + clhs182 + 0.666666666666667*clhs184 + clhs187 + clhs188 + clhs189 + clhs72 + clhs73 + clhs74 + clhs75 + clhs76 + clhs77 + clhs91 + 0.666666666666667*clhs92 + clhs95 + 0.666666666666667*clhs97;
            lhs(2,2)=clhs12 + pow(clhs121, 2)*tau[0] + pow(clhs124, 2)*tau[1] + 0.666666666666667*clhs125 + 0.666666666666667*clhs126 + pow(clhs127, 2)*tau[2] + clhs130 + clhs132 + clhs135 - 0.666666666666667*clhs136 + clhs155 + clhs158 + clhs172 + clhs173 + clhs174 + clhs180 + 0.666666666666667*clhs181 + clhs185 + 0.666666666666667*clhs186 + clhs188 + clhs189 + clhs190*clhs33 + clhs191*clhs33 - clhs192*tau[0] - clhs192*tau[1] - clhs192*tau[2] - clhs193*tau[0] - clhs193*tau[1] - clhs193*tau[2] + clhs26 + clhs27 + clhs28 + clhs29 + clhs3 + clhs30 + clhs31 - 0.444444444444444*clhs38 - 0.444444444444444*clhs55;


    const double crhs0 =             0.25*f[1];
const double crhs1 =             phi_subscale_gauss[1]*tau[1];
const double crhs2 =             -0.166666666666667*crhs1;
const double crhs3 =             1.0/delta_time;
const double crhs4 =             crhs3*phi_subscale_gauss[1];
const double crhs5 =             0.166666666666667*crhs4;
const double crhs6 =             0.166666666666667*f[0];
const double crhs7 =             0.166666666666667*f[2];
const double crhs8 =             tau[1]*(crhs6 + crhs7 + 0.666666666666667*f[1]);
const double crhs9 =             -0.166666666666667*crhs8;
const double crhs10 =             0.166666666666667*prj[0];
const double crhs11 =             0.166666666666667*prj[2];
const double crhs12 =             tau[1]*(crhs10 + crhs11 + 0.666666666666667*prj[1]);
const double crhs13 =             -0.166666666666667*crhs12;
const double crhs14 =             0.166666666666667*v(0,0);
const double crhs15 =             0.166666666666667*v(2,0);
const double crhs16 =             crhs14 + crhs15 + 0.666666666666667*v(1,0);
const double crhs17 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs18 =             crhs16*crhs17;
const double crhs19 =             -0.166666666666667*crhs18;
const double crhs20 =             0.166666666666667*v(0,1);
const double crhs21 =             0.166666666666667*v(2,1);
const double crhs22 =             crhs20 + crhs21 + 0.666666666666667*v(1,1);
const double crhs23 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs24 =             crhs22*crhs23;
const double crhs25 =             -0.166666666666667*crhs24;
const double crhs26 =             1.0/RK_time_coefficient;
const double crhs27 =             crhs26*crhs3*tau[1];
const double crhs28 =             0.111111111111111*phi[1];
const double crhs29 =             -0.111111111111111*phi_old[1];
const double crhs30 =             crhs28 + crhs29;
const double crhs31 =             0.0277777777777778*phi[0];
const double crhs32 =             0.0277777777777778*phi[2];
const double crhs33 =             -0.0277777777777778*phi_old[0];
const double crhs34 =             -0.0277777777777778*phi_old[2];
const double crhs35 =             crhs27*(crhs30 + crhs31 + crhs32 + crhs33 + crhs34);
const double crhs36 =             0.166666666666667*phi[0];
const double crhs37 =             0.166666666666667*phi[2];
const double crhs38 =             crhs36 + crhs37 + 0.666666666666667*phi[1];
const double crhs39 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs40 =             crhs38*crhs39;
const double crhs41 =             0.166666666666667*crhs40;
const double crhs42 =             -crhs41;
const double crhs43 =             crhs41*tau[1];
const double crhs44 =             tau[1]*(crhs18 + crhs24);
const double crhs45 =             0.166666666666667*crhs44;
const double crhs46 =             0.25*f[2];
const double crhs47 =             phi_subscale_gauss[0]*tau[0];
const double crhs48 =             phi_subscale_gauss[2]*tau[2];
const double crhs49 =             -0.166666666666667*crhs48;
const double crhs50 =             crhs3*phi_subscale_gauss[0];
const double crhs51 =             crhs3*phi_subscale_gauss[2];
const double crhs52 =             0.166666666666667*crhs51;
const double crhs53 =             0.166666666666667*f[1];
const double crhs54 =             tau[0]*(crhs53 + crhs7 + 0.666666666666667*f[0]);
const double crhs55 =             0.166666666666667*prj[1];
const double crhs56 =             tau[0]*(crhs11 + crhs55 + 0.666666666666667*prj[0]);
const double crhs57 =             tau[2]*(crhs53 + crhs6 + 0.666666666666667*f[2]);
const double crhs58 =             -0.166666666666667*crhs57;
const double crhs59 =             tau[2]*(crhs10 + crhs55 + 0.666666666666667*prj[2]);
const double crhs60 =             -0.166666666666667*crhs59;
const double crhs61 =             3*crhs17*k;
const double crhs62 =             3*crhs23*k;
const double crhs63 =             DN_DX_0_0*crhs17*k;
const double crhs64 =             DN_DX_0_1*crhs23*k;
const double crhs65 =             0.166666666666667*v(1,0);
const double crhs66 =             crhs14 + crhs65 + 0.666666666666667*v(2,0);
const double crhs67 =             crhs17*crhs66;
const double crhs68 =             -0.166666666666667*crhs67;
const double crhs69 =             crhs15 + crhs65 + 0.666666666666667*v(0,0);
const double crhs70 =             crhs17*crhs69;
const double crhs71 =             0.166666666666667*v(1,1);
const double crhs72 =             crhs20 + crhs71 + 0.666666666666667*v(2,1);
const double crhs73 =             crhs23*crhs72;
const double crhs74 =             -0.166666666666667*crhs73;
const double crhs75 =             crhs21 + crhs71 + 0.666666666666667*v(0,1);
const double crhs76 =             crhs23*crhs75;
const double crhs77 =             crhs26*crhs3*tau[0];
const double crhs78 =             0.111111111111111*phi[2];
const double crhs79 =             -0.111111111111111*phi_old[2];
const double crhs80 =             crhs26*crhs3*tau[2];
const double crhs81 =             0.0277777777777778*phi[1];
const double crhs82 =             -0.0277777777777778*phi_old[1];
const double crhs83 =             crhs80*(crhs31 + crhs33 + crhs78 + crhs79 + crhs81 + crhs82);
const double crhs84 =             0.166666666666667*phi[1];
const double crhs85 =             crhs36 + crhs84 + 0.666666666666667*phi[2];
const double crhs86 =             crhs39*crhs85;
const double crhs87 =             0.166666666666667*crhs86;
const double crhs88 =             -crhs87;
const double crhs89 =             crhs37 + crhs84 + 0.666666666666667*phi[0];
const double crhs90 =             crhs39*crhs89;
const double crhs91 =             0.666666666666667*crhs90;
const double crhs92 =             crhs3*phi_subscale_gauss[0]*tau[0];
const double crhs93 =             DN_DX_0_0*crhs69 + DN_DX_0_1*crhs75;
const double crhs94 =             crhs3*phi_subscale_gauss[1]*tau[1];
const double crhs95 =             DN_DX_0_0*crhs16 + DN_DX_0_1*crhs22;
const double crhs96 =             crhs3*phi_subscale_gauss[2]*tau[2];
const double crhs97 =             DN_DX_0_0*crhs66 + DN_DX_0_1*crhs72;
const double crhs98 =             crhs87*tau[2];
const double crhs99 =             tau[0]*(crhs70 + crhs76);
const double crhs100 =             tau[2]*(crhs67 + crhs73);
const double crhs101 =             0.166666666666667*crhs100;
const double crhs102 =             -0.166666666666667*phi_old[1];
const double crhs103 =             -0.166666666666667*phi_old[2];
const double crhs104 =             crhs26*crhs3*tau[0]*(crhs102 + crhs103 + crhs89 - 0.666666666666667*phi_old[0]);
const double crhs105 =             -0.166666666666667*phi_old[0];
const double crhs106 =             crhs26*crhs3*tau[1]*(crhs103 + crhs105 + crhs38 - 0.666666666666667*phi_old[1]);
const double crhs107 =             crhs26*crhs3*tau[2]*(crhs102 + crhs105 + crhs85 - 0.666666666666667*phi_old[2]);
const double crhs108 =             crhs39*crhs89*tau[0];
const double crhs109 =             crhs38*crhs39*tau[1];
const double crhs110 =             crhs39*crhs85*tau[2];
const double crhs111 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs112 =             0.166666666666667*crhs90;
const double crhs113 =             crhs112*tau[0] - crhs112 - 0.166666666666667*crhs47 + 0.166666666666667*crhs50 - 0.166666666666667*crhs54 - 0.166666666666667*crhs56 - 0.166666666666667*crhs70 - 0.166666666666667*crhs76 + crhs77*(crhs111 + crhs32 + crhs34 + crhs81 + crhs82) + 0.166666666666667*crhs99 + 0.25*f[0];
const double crhs114 =             DN_DX_1_0*crhs17*k;
const double crhs115 =             DN_DX_1_1*crhs23*k;
const double crhs116 =             0.666666666666667*crhs40;
const double crhs117 =             DN_DX_1_0*crhs69 + DN_DX_1_1*crhs75;
const double crhs118 =             DN_DX_1_0*crhs16 + DN_DX_1_1*crhs22;
const double crhs119 =             DN_DX_1_0*crhs66 + DN_DX_1_1*crhs72;
const double crhs120 =             DN_DX_2_0*crhs17*k;
const double crhs121 =             DN_DX_2_1*crhs23*k;
const double crhs122 =             0.666666666666667*crhs86;
const double crhs123 =             DN_DX_2_0*crhs69 + DN_DX_2_1*crhs75;
const double crhs124 =             DN_DX_2_0*crhs16 + DN_DX_2_1*crhs22;
const double crhs125 =             DN_DX_2_0*crhs66 + DN_DX_2_1*crhs72;
            rhs[0]=-DN_DX_0_0*crhs61 - DN_DX_0_1*crhs62 + crhs0 - crhs100*crhs97 + crhs101 - crhs104*crhs93 - crhs106*crhs95 - crhs107*crhs97 - crhs108*crhs93 - crhs109*crhs95 - crhs110*crhs97 + crhs12*crhs95 + crhs13 + crhs19 + crhs2 + crhs25 + crhs35 + crhs42 + crhs43 - crhs44*crhs95 + crhs45 + crhs46 - 0.666666666666667*crhs47 + crhs49 + crhs5 + 0.666666666666667*crhs50 + crhs52 + crhs54*crhs93 - 0.666666666666667*crhs54 + crhs56*crhs93 - 0.666666666666667*crhs56 + crhs57*crhs97 + crhs58 + crhs59*crhs97 + crhs60 + crhs63*tau[0] + crhs63*tau[1] + crhs63*tau[2] + crhs64*tau[0] + crhs64*tau[1] + crhs64*tau[2] + crhs68 - 0.666666666666667*crhs70 + crhs74 - 0.666666666666667*crhs76 + crhs77*(crhs30 + crhs78 + crhs79 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs8*crhs95 + crhs83 + crhs88 + crhs9 + crhs91*tau[0] - crhs91 + crhs92*crhs93 - crhs93*crhs99 + crhs94*crhs95 + crhs96*crhs97 + crhs98 + 0.666666666666667*crhs99 + 0.5*f[0];
            rhs[1]=-DN_DX_1_0*crhs61 - DN_DX_1_1*crhs62 - 0.666666666666667*crhs1 - crhs100*crhs119 + crhs101 - crhs104*crhs117 - crhs106*crhs118 - crhs107*crhs119 - crhs108*crhs117 - crhs109*crhs118 - crhs110*crhs119 + crhs113 + crhs114*tau[0] + crhs114*tau[1] + crhs114*tau[2] + crhs115*tau[0] + crhs115*tau[1] + crhs115*tau[2] + crhs116*tau[1] - crhs116 + crhs117*crhs54 + crhs117*crhs56 + crhs117*crhs92 - crhs117*crhs99 + crhs118*crhs12 - crhs118*crhs44 + crhs118*crhs8 + crhs118*crhs94 + crhs119*crhs57 + crhs119*crhs59 + crhs119*crhs96 - 0.666666666666667*crhs12 - 0.666666666666667*crhs18 - 0.666666666666667*crhs24 + crhs27*(crhs111 + crhs78 + crhs79 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) + 0.666666666666667*crhs4 + 0.666666666666667*crhs44 + crhs46 + crhs49 + crhs52 + crhs58 + crhs60 + crhs68 + crhs74 - 0.666666666666667*crhs8 + crhs83 + crhs88 + crhs98 + 0.5*f[1];
            rhs[2]=-DN_DX_2_0*crhs61 - DN_DX_2_1*crhs62 + crhs0 - crhs100*crhs125 + 0.666666666666667*crhs100 - crhs104*crhs123 - crhs106*crhs124 - crhs107*crhs125 - crhs108*crhs123 - crhs109*crhs124 - crhs110*crhs125 + crhs113 + crhs12*crhs124 + crhs120*tau[0] + crhs120*tau[1] + crhs120*tau[2] + crhs121*tau[0] + crhs121*tau[1] + crhs121*tau[2] + crhs122*tau[2] - crhs122 + crhs123*crhs54 + crhs123*crhs56 + crhs123*crhs92 - crhs123*crhs99 - crhs124*crhs44 + crhs124*crhs8 + crhs124*crhs94 + crhs125*crhs57 + crhs125*crhs59 + crhs125*crhs96 + crhs13 + crhs19 + crhs2 + crhs25 + crhs35 + crhs42 + crhs43 + crhs45 - 0.666666666666667*crhs48 + crhs5 + 0.666666666666667*crhs51 - 0.666666666666667*crhs57 - 0.666666666666667*crhs59 - 0.666666666666667*crhs67 - 0.666666666666667*crhs73 + crhs80*(crhs111 + crhs28 + crhs29 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs9 + 0.5*f[2];


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
const double clhs50 =             pow(DN_DX_0_0, 2);
const double clhs51 =             4*k;
const double clhs52 =             pow(DN_DX_0_1, 2);
const double clhs53 =             pow(DN_DX_0_2, 2);
const double clhs54 =             clhs50*k;
const double clhs55 =             clhs52*k;
const double clhs56 =             clhs53*k;
const double clhs57 =             clhs0*clhs1*tau[0];
const double clhs58 =             clhs0*clhs1*tau[1];
const double clhs59 =             -0.01909830025156*clhs58;
const double clhs60 =             clhs0*clhs1*tau[2];
const double clhs61 =             -0.01909830025156*clhs60;
const double clhs62 =             clhs0*clhs1*tau[3];
const double clhs63 =             -0.01909830025156*clhs62;
const double clhs64 =             0.1381966*v(1,0);
const double clhs65 =             clhs2 + clhs64;
const double clhs66 =             clhs3 + clhs65 + 0.5854102*v(3,0);
const double clhs67 =             DN_DX_0_0*clhs66;
const double clhs68 =             0.1381966*clhs67;
const double clhs69 =             clhs4 + clhs65 + 0.5854102*v(2,0);
const double clhs70 =             DN_DX_0_0*clhs69;
const double clhs71 =             0.1381966*clhs70;
const double clhs72 =             0.1381966*clhs7;
const double clhs73 =             clhs5 + clhs64 + 0.5854102*v(0,0);
const double clhs74 =             DN_DX_0_0*clhs73;
const double clhs75 =             0.1381966*v(1,1);
const double clhs76 =             clhs75 + clhs8;
const double clhs77 =             clhs76 + clhs9 + 0.5854102*v(3,1);
const double clhs78 =             DN_DX_0_1*clhs77;
const double clhs79 =             0.1381966*clhs78;
const double clhs80 =             clhs10 + clhs76 + 0.5854102*v(2,1);
const double clhs81 =             DN_DX_0_1*clhs80;
const double clhs82 =             0.1381966*clhs81;
const double clhs83 =             0.1381966*clhs13;
const double clhs84 =             clhs11 + clhs75 + 0.5854102*v(0,1);
const double clhs85 =             DN_DX_0_1*clhs84;
const double clhs86 =             0.1381966*v(1,2);
const double clhs87 =             clhs14 + clhs86;
const double clhs88 =             clhs15 + clhs87 + 0.5854102*v(3,2);
const double clhs89 =             DN_DX_0_2*clhs88;
const double clhs90 =             0.1381966*clhs89;
const double clhs91 =             clhs16 + clhs87 + 0.5854102*v(2,2);
const double clhs92 =             DN_DX_0_2*clhs91;
const double clhs93 =             0.1381966*clhs92;
const double clhs94 =             0.1381966*clhs19;
const double clhs95 =             clhs17 + clhs86 + 0.5854102*v(0,2);
const double clhs96 =             DN_DX_0_2*clhs95;
const double clhs97 =             clhs35*tau[0];
const double clhs98 =             clhs35*tau[1];
const double clhs99 =             -0.01909830025156*clhs98;
const double clhs100 =             clhs35*tau[2];
const double clhs101 =             -0.01909830025156*clhs100;
const double clhs102 =             clhs35*tau[3];
const double clhs103 =             -0.01909830025156*clhs102;
const double clhs104 =             clhs74 + clhs85 + clhs96;
const double clhs105 =             clhs104*tau[0];
const double clhs106 =             clhs20*tau[1];
const double clhs107 =             -0.1381966*clhs106;
const double clhs108 =             clhs70 + clhs81 + clhs92;
const double clhs109 =             clhs108*tau[2];
const double clhs110 =             -0.1381966*clhs109;
const double clhs111 =             clhs67 + clhs78 + clhs89;
const double clhs112 =             clhs111*tau[3];
const double clhs113 =             -0.1381966*clhs112;
const double clhs114 =             clhs0*clhs1*clhs104*tau[0];
const double clhs115 =             clhs0*clhs1*clhs108*tau[2];
const double clhs116 =             0.1381966*clhs115;
const double clhs117 =             clhs0*clhs1*clhs111*tau[3];
const double clhs118 =             0.1381966*clhs117;
const double clhs119 =             clhs104*clhs35*tau[0];
const double clhs120 =             clhs108*clhs35*tau[2];
const double clhs121 =             0.1381966*clhs120;
const double clhs122 =             clhs111*clhs35*tau[3];
const double clhs123 =             0.1381966*clhs122;
const double clhs124 =             0.1381966*clhs114;
const double clhs125 =             0.1381966*clhs119;
const double clhs126 =             0.19999999899376*clhs23;
const double clhs127 =             0.19999999899376*clhs24;
const double clhs128 =             0.19999999899376*clhs25;
const double clhs129 =             0.19999999899376*clhs26;
const double clhs130 =             0.19999999899376*clhs27;
const double clhs131 =             0.19999999899376*clhs28;
const double clhs132 =             0.19999999899376*clhs29;
const double clhs133 =             0.19999999899376*clhs30;
const double clhs134 =             0.19999999899376*clhs31;
const double clhs135 =             0.19999999899376*clhs32;
const double clhs136 =             0.19999999899376*clhs33;
const double clhs137 =             0.19999999899376*clhs34;
const double clhs138 =             4*DN_DX_0_0*k;
const double clhs139 =             DN_DX_1_0*clhs138;
const double clhs140 =             4*DN_DX_0_1*k;
const double clhs141 =             DN_DX_1_1*clhs140;
const double clhs142 =             4*DN_DX_0_2*k;
const double clhs143 =             DN_DX_1_2*clhs142;
const double clhs144 =             DN_DX_0_0*DN_DX_1_0*k;
const double clhs145 =             -clhs144*tau[0];
const double clhs146 =             -clhs144*tau[1];
const double clhs147 =             -clhs144*tau[2];
const double clhs148 =             -clhs144*tau[3];
const double clhs149 =             DN_DX_0_1*DN_DX_1_1*k;
const double clhs150 =             -clhs149*tau[0];
const double clhs151 =             -clhs149*tau[1];
const double clhs152 =             -clhs149*tau[2];
const double clhs153 =             -clhs149*tau[3];
const double clhs154 =             DN_DX_0_2*DN_DX_1_2*k;
const double clhs155 =             -clhs154*tau[0];
const double clhs156 =             -clhs154*tau[1];
const double clhs157 =             -clhs154*tau[2];
const double clhs158 =             -clhs154*tau[3];
const double clhs159 =             -0.08090169924532*clhs57;
const double clhs160 =             -0.08090169924532*clhs58;
const double clhs161 =             DN_DX_1_0*clhs66;
const double clhs162 =             0.1381966*clhs161;
const double clhs163 =             DN_DX_1_0*clhs69;
const double clhs164 =             0.1381966*clhs163;
const double clhs165 =             DN_DX_1_0*clhs6;
const double clhs166 =             0.1381966*clhs165;
const double clhs167 =             DN_DX_1_0*clhs73;
const double clhs168 =             DN_DX_1_1*clhs77;
const double clhs169 =             0.1381966*clhs168;
const double clhs170 =             DN_DX_1_1*clhs80;
const double clhs171 =             0.1381966*clhs170;
const double clhs172 =             DN_DX_1_1*clhs12;
const double clhs173 =             0.1381966*clhs172;
const double clhs174 =             DN_DX_1_1*clhs84;
const double clhs175 =             DN_DX_1_2*clhs88;
const double clhs176 =             0.1381966*clhs175;
const double clhs177 =             DN_DX_1_2*clhs91;
const double clhs178 =             0.1381966*clhs177;
const double clhs179 =             DN_DX_1_2*clhs18;
const double clhs180 =             0.1381966*clhs179;
const double clhs181 =             DN_DX_1_2*clhs95;
const double clhs182 =             -0.08090169924532*clhs97;
const double clhs183 =             -0.08090169924532*clhs98;
const double clhs184 =             clhs167 + clhs174 + clhs181;
const double clhs185 =             clhs184*tau[0];
const double clhs186 =             clhs165 + clhs172 + clhs179;
const double clhs187 =             clhs186*tau[1];
const double clhs188 =             -0.1381966*clhs187;
const double clhs189 =             clhs163 + clhs170 + clhs177;
const double clhs190 =             clhs189*tau[2];
const double clhs191 =             -0.1381966*clhs190;
const double clhs192 =             clhs161 + clhs168 + clhs175;
const double clhs193 =             clhs192*tau[3];
const double clhs194 =             -0.1381966*clhs193;
const double clhs195 =             clhs105*clhs184;
const double clhs196 =             clhs106*clhs186;
const double clhs197 =             clhs109*clhs189;
const double clhs198 =             clhs112*clhs192;
const double clhs199 =             DN_DX_2_0*clhs138;
const double clhs200 =             DN_DX_2_1*clhs140;
const double clhs201 =             DN_DX_2_2*clhs142;
const double clhs202 =             DN_DX_0_0*DN_DX_2_0*k;
const double clhs203 =             -clhs202*tau[0];
const double clhs204 =             -clhs202*tau[1];
const double clhs205 =             -clhs202*tau[2];
const double clhs206 =             -clhs202*tau[3];
const double clhs207 =             DN_DX_0_1*DN_DX_2_1*k;
const double clhs208 =             -clhs207*tau[0];
const double clhs209 =             -clhs207*tau[1];
const double clhs210 =             -clhs207*tau[2];
const double clhs211 =             -clhs207*tau[3];
const double clhs212 =             DN_DX_0_2*DN_DX_2_2*k;
const double clhs213 =             -clhs212*tau[0];
const double clhs214 =             -clhs212*tau[1];
const double clhs215 =             -clhs212*tau[2];
const double clhs216 =             -clhs212*tau[3];
const double clhs217 =             -0.08090169924532*clhs60;
const double clhs218 =             DN_DX_2_0*clhs66;
const double clhs219 =             0.1381966*clhs218;
const double clhs220 =             DN_DX_2_0*clhs69;
const double clhs221 =             0.1381966*clhs220;
const double clhs222 =             DN_DX_2_0*clhs6;
const double clhs223 =             0.1381966*clhs222;
const double clhs224 =             DN_DX_2_0*clhs73;
const double clhs225 =             DN_DX_2_1*clhs77;
const double clhs226 =             0.1381966*clhs225;
const double clhs227 =             DN_DX_2_1*clhs80;
const double clhs228 =             0.1381966*clhs227;
const double clhs229 =             DN_DX_2_1*clhs12;
const double clhs230 =             0.1381966*clhs229;
const double clhs231 =             DN_DX_2_1*clhs84;
const double clhs232 =             DN_DX_2_2*clhs88;
const double clhs233 =             0.1381966*clhs232;
const double clhs234 =             DN_DX_2_2*clhs91;
const double clhs235 =             0.1381966*clhs234;
const double clhs236 =             DN_DX_2_2*clhs18;
const double clhs237 =             0.1381966*clhs236;
const double clhs238 =             DN_DX_2_2*clhs95;
const double clhs239 =             -0.08090169924532*clhs100;
const double clhs240 =             clhs224 + clhs231 + clhs238;
const double clhs241 =             clhs240*tau[0];
const double clhs242 =             clhs222 + clhs229 + clhs236;
const double clhs243 =             clhs242*tau[1];
const double clhs244 =             -0.1381966*clhs243;
const double clhs245 =             clhs220 + clhs227 + clhs234;
const double clhs246 =             clhs245*tau[2];
const double clhs247 =             -0.1381966*clhs246;
const double clhs248 =             clhs218 + clhs225 + clhs232;
const double clhs249 =             clhs248*tau[3];
const double clhs250 =             -0.1381966*clhs249;
const double clhs251 =             clhs105*clhs240;
const double clhs252 =             clhs106*clhs242;
const double clhs253 =             clhs109*clhs245;
const double clhs254 =             clhs112*clhs248;
const double clhs255 =             DN_DX_0_0*DN_DX_3_0*k;
const double clhs256 =             DN_DX_0_1*DN_DX_3_1*k;
const double clhs257 =             DN_DX_0_2*DN_DX_3_2*k;
const double clhs258 =             -0.08090169924532*clhs62;
const double clhs259 =             -0.08090169924532*clhs102;
const double clhs260 =             DN_DX_3_0*clhs73;
const double clhs261 =             DN_DX_3_1*clhs84;
const double clhs262 =             DN_DX_3_2*clhs95;
const double clhs263 =             clhs260 + clhs261 + clhs262;
const double clhs264 =             DN_DX_3_0*clhs6;
const double clhs265 =             DN_DX_3_1*clhs12;
const double clhs266 =             DN_DX_3_2*clhs18;
const double clhs267 =             clhs264 + clhs265 + clhs266;
const double clhs268 =             DN_DX_3_0*clhs69;
const double clhs269 =             DN_DX_3_1*clhs80;
const double clhs270 =             DN_DX_3_2*clhs91;
const double clhs271 =             clhs268 + clhs269 + clhs270;
const double clhs272 =             DN_DX_3_0*clhs66;
const double clhs273 =             DN_DX_3_1*clhs77;
const double clhs274 =             DN_DX_3_2*clhs88;
const double clhs275 =             clhs272 + clhs273 + clhs274;
const double clhs276 =             DN_DX_3_0*clhs138 + DN_DX_3_1*clhs140 + DN_DX_3_2*clhs142 + clhs105*clhs263 + clhs106*clhs267 + clhs109*clhs271 + clhs112*clhs275 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs159 + clhs182 - clhs255*tau[0] - clhs255*tau[1] - clhs255*tau[2] - clhs255*tau[3] - clhs256*tau[0] - clhs256*tau[1] - clhs256*tau[2] - clhs256*tau[3] - clhs257*tau[0] - clhs257*tau[1] - clhs257*tau[2] - clhs257*tau[3] + clhs258 + clhs259;
const double clhs277 =             0.1381966*clhs272;
const double clhs278 =             0.1381966*clhs268;
const double clhs279 =             0.1381966*clhs264;
const double clhs280 =             0.1381966*clhs273;
const double clhs281 =             0.1381966*clhs269;
const double clhs282 =             0.1381966*clhs265;
const double clhs283 =             0.1381966*clhs274;
const double clhs284 =             0.1381966*clhs270;
const double clhs285 =             0.1381966*clhs266;
const double clhs286 =             clhs263*tau[0];
const double clhs287 =             clhs267*tau[1];
const double clhs288 =             -0.1381966*clhs287;
const double clhs289 =             clhs271*tau[2];
const double clhs290 =             -0.1381966*clhs289;
const double clhs291 =             clhs275*tau[3];
const double clhs292 =             -0.1381966*clhs291;
const double clhs293 =             0.1381966*clhs74;
const double clhs294 =             0.1381966*clhs85;
const double clhs295 =             0.1381966*clhs96;
const double clhs296 =             -0.1381966*clhs105;
const double clhs297 =             clhs0*clhs1*clhs184*tau[0];
const double clhs298 =             clhs0*clhs1*clhs186*tau[1];
const double clhs299 =             0.1381966*clhs298;
const double clhs300 =             clhs0*clhs1*clhs189*tau[2];
const double clhs301 =             0.1381966*clhs300;
const double clhs302 =             clhs0*clhs1*clhs192*tau[3];
const double clhs303 =             0.1381966*clhs302;
const double clhs304 =             clhs184*clhs35*tau[0];
const double clhs305 =             clhs186*clhs35*tau[1];
const double clhs306 =             0.1381966*clhs305;
const double clhs307 =             clhs189*clhs35*tau[2];
const double clhs308 =             0.1381966*clhs307;
const double clhs309 =             clhs192*clhs35*tau[3];
const double clhs310 =             0.1381966*clhs309;
const double clhs311 =             pow(DN_DX_1_0, 2);
const double clhs312 =             pow(DN_DX_1_1, 2);
const double clhs313 =             pow(DN_DX_1_2, 2);
const double clhs314 =             clhs311*k;
const double clhs315 =             clhs312*k;
const double clhs316 =             clhs313*k;
const double clhs317 =             -0.01909830025156*clhs57;
const double clhs318 =             0.1381966*clhs167;
const double clhs319 =             0.1381966*clhs174;
const double clhs320 =             0.1381966*clhs181;
const double clhs321 =             -0.01909830025156*clhs97;
const double clhs322 =             -0.1381966*clhs185;
const double clhs323 =             0.1381966*clhs297;
const double clhs324 =             0.1381966*clhs304;
const double clhs325 =             4*DN_DX_1_0*k;
const double clhs326 =             DN_DX_2_0*clhs325;
const double clhs327 =             4*DN_DX_1_1*k;
const double clhs328 =             DN_DX_2_1*clhs327;
const double clhs329 =             4*DN_DX_1_2*k;
const double clhs330 =             DN_DX_2_2*clhs329;
const double clhs331 =             DN_DX_1_0*DN_DX_2_0*k;
const double clhs332 =             -clhs331*tau[0];
const double clhs333 =             -clhs331*tau[1];
const double clhs334 =             -clhs331*tau[2];
const double clhs335 =             -clhs331*tau[3];
const double clhs336 =             DN_DX_1_1*DN_DX_2_1*k;
const double clhs337 =             -clhs336*tau[0];
const double clhs338 =             -clhs336*tau[1];
const double clhs339 =             -clhs336*tau[2];
const double clhs340 =             -clhs336*tau[3];
const double clhs341 =             DN_DX_1_2*DN_DX_2_2*k;
const double clhs342 =             -clhs341*tau[0];
const double clhs343 =             -clhs341*tau[1];
const double clhs344 =             -clhs341*tau[2];
const double clhs345 =             -clhs341*tau[3];
const double clhs346 =             0.1381966*clhs224;
const double clhs347 =             0.1381966*clhs231;
const double clhs348 =             0.1381966*clhs238;
const double clhs349 =             -0.1381966*clhs241;
const double clhs350 =             clhs185*clhs240;
const double clhs351 =             clhs187*clhs242;
const double clhs352 =             clhs190*clhs245;
const double clhs353 =             clhs193*clhs248;
const double clhs354 =             DN_DX_3_0*clhs325;
const double clhs355 =             DN_DX_3_1*clhs327;
const double clhs356 =             DN_DX_3_2*clhs329;
const double clhs357 =             DN_DX_1_0*DN_DX_3_0*k;
const double clhs358 =             -clhs357*tau[0];
const double clhs359 =             -clhs357*tau[1];
const double clhs360 =             -clhs357*tau[2];
const double clhs361 =             -clhs357*tau[3];
const double clhs362 =             DN_DX_1_1*DN_DX_3_1*k;
const double clhs363 =             -clhs362*tau[0];
const double clhs364 =             -clhs362*tau[1];
const double clhs365 =             -clhs362*tau[2];
const double clhs366 =             -clhs362*tau[3];
const double clhs367 =             DN_DX_1_2*DN_DX_3_2*k;
const double clhs368 =             -clhs367*tau[0];
const double clhs369 =             -clhs367*tau[1];
const double clhs370 =             -clhs367*tau[2];
const double clhs371 =             -clhs367*tau[3];
const double clhs372 =             0.1381966*clhs260;
const double clhs373 =             0.1381966*clhs261;
const double clhs374 =             0.1381966*clhs262;
const double clhs375 =             -0.1381966*clhs286;
const double clhs376 =             clhs185*clhs263;
const double clhs377 =             clhs187*clhs267;
const double clhs378 =             clhs190*clhs271;
const double clhs379 =             clhs193*clhs275;
const double clhs380 =             clhs0*clhs1*clhs240*tau[0];
const double clhs381 =             clhs0*clhs1*clhs242*tau[1];
const double clhs382 =             0.1381966*clhs381;
const double clhs383 =             clhs0*clhs1*clhs245*tau[2];
const double clhs384 =             0.1381966*clhs383;
const double clhs385 =             clhs0*clhs1*clhs248*tau[3];
const double clhs386 =             0.1381966*clhs385;
const double clhs387 =             clhs240*clhs35*tau[0];
const double clhs388 =             clhs242*clhs35*tau[1];
const double clhs389 =             0.1381966*clhs388;
const double clhs390 =             clhs245*clhs35*tau[2];
const double clhs391 =             0.1381966*clhs390;
const double clhs392 =             clhs248*clhs35*tau[3];
const double clhs393 =             0.1381966*clhs392;
const double clhs394 =             0.1381966*clhs380;
const double clhs395 =             0.1381966*clhs387;
const double clhs396 =             pow(DN_DX_2_0, 2);
const double clhs397 =             pow(DN_DX_2_1, 2);
const double clhs398 =             pow(DN_DX_2_2, 2);
const double clhs399 =             clhs396*k;
const double clhs400 =             clhs397*k;
const double clhs401 =             clhs398*k;
const double clhs402 =             DN_DX_2_0*DN_DX_3_0*clhs51;
const double clhs403 =             DN_DX_2_1*DN_DX_3_1*clhs51;
const double clhs404 =             DN_DX_2_2*DN_DX_3_2*clhs51;
const double clhs405 =             DN_DX_2_0*DN_DX_3_0*k;
const double clhs406 =             -clhs405*tau[0];
const double clhs407 =             -clhs405*tau[1];
const double clhs408 =             -clhs405*tau[2];
const double clhs409 =             -clhs405*tau[3];
const double clhs410 =             DN_DX_2_1*DN_DX_3_1*k;
const double clhs411 =             -clhs410*tau[0];
const double clhs412 =             -clhs410*tau[1];
const double clhs413 =             -clhs410*tau[2];
const double clhs414 =             -clhs410*tau[3];
const double clhs415 =             DN_DX_2_2*DN_DX_3_2*k;
const double clhs416 =             -clhs415*tau[0];
const double clhs417 =             -clhs415*tau[1];
const double clhs418 =             -clhs415*tau[2];
const double clhs419 =             -clhs415*tau[3];
const double clhs420 =             clhs241*clhs263;
const double clhs421 =             clhs243*clhs267;
const double clhs422 =             clhs246*clhs271;
const double clhs423 =             clhs249*clhs275;
const double clhs424 =             clhs0*clhs1*clhs263*tau[0];
const double clhs425 =             clhs0*clhs1*clhs267*tau[1];
const double clhs426 =             0.1381966*clhs425;
const double clhs427 =             clhs0*clhs1*clhs271*tau[2];
const double clhs428 =             0.1381966*clhs427;
const double clhs429 =             clhs0*clhs1*clhs275*tau[3];
const double clhs430 =             0.1381966*clhs429;
const double clhs431 =             clhs263*clhs35*tau[0];
const double clhs432 =             clhs267*clhs35*tau[1];
const double clhs433 =             0.1381966*clhs432;
const double clhs434 =             clhs271*clhs35*tau[2];
const double clhs435 =             0.1381966*clhs434;
const double clhs436 =             clhs275*clhs35*tau[3];
const double clhs437 =             0.1381966*clhs436;
const double clhs438 =             0.1381966*clhs424;
const double clhs439 =             0.1381966*clhs431;
const double clhs440 =             pow(DN_DX_3_0, 2);
const double clhs441 =             pow(DN_DX_3_1, 2);
const double clhs442 =             pow(DN_DX_3_2, 2);
const double clhs443 =             clhs440*k;
const double clhs444 =             clhs441*k;
const double clhs445 =             clhs442*k;
            lhs(0,0)=clhs101 + clhs103 + pow(clhs104, 2)*tau[0] - 0.5854102*clhs105 + clhs107 + pow(clhs108, 2)*tau[2] + clhs110 + pow(clhs111, 2)*tau[3] + clhs113 + 0.5854102*clhs114 + clhs116 + clhs118 + 0.5854102*clhs119 + clhs121 + clhs123 + pow(clhs20, 2)*tau[1] + clhs22 + clhs37 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49 + clhs50*clhs51 + clhs51*clhs52 + clhs51*clhs53 - clhs54*tau[0] - clhs54*tau[1] - clhs54*tau[2] - clhs54*tau[3] - clhs55*tau[0] - clhs55*tau[1] - clhs55*tau[2] - clhs55*tau[3] - clhs56*tau[0] - clhs56*tau[1] - clhs56*tau[2] - clhs56*tau[3] - 0.34270510226404*clhs57 + clhs59 + clhs61 + clhs63 + clhs68 + clhs71 + clhs72 + 0.5854102*clhs74 + clhs79 + clhs82 + clhs83 + 0.5854102*clhs85 + clhs90 + clhs93 + clhs94 + 0.5854102*clhs96 - 0.34270510226404*clhs97 + clhs99;
            lhs(0,1)=clhs101 + clhs103 + clhs116 + clhs118 + clhs121 + clhs123 + clhs124 + clhs125 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs139 + clhs141 + clhs143 + clhs145 + clhs146 + clhs147 + clhs148 + clhs150 + clhs151 + clhs152 + clhs153 + clhs155 + clhs156 + clhs157 + clhs158 + clhs159 + clhs160 + clhs162 + clhs164 + clhs166 + 0.5854102*clhs167 + clhs169 + clhs171 + clhs173 + 0.5854102*clhs174 + clhs176 + clhs178 + clhs180 + 0.5854102*clhs181 + clhs182 + clhs183 - 0.5854102*clhs185 + clhs188 + clhs191 + clhs194 + clhs195 + clhs196 + clhs197 + clhs198 + 0.5854102*clhs21 + 0.5854102*clhs36 + clhs61 + clhs63;
            lhs(0,2)=clhs103 + 0.5854102*clhs115 + clhs118 + 0.5854102*clhs120 + clhs123 + clhs124 + clhs125 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs159 + clhs182 + clhs199 + clhs200 + clhs201 + clhs203 + clhs204 + clhs205 + clhs206 + clhs208 + clhs209 + clhs210 + clhs211 + clhs213 + clhs214 + clhs215 + clhs216 + clhs217 + clhs219 + clhs22 + clhs221 + clhs223 + 0.5854102*clhs224 + clhs226 + clhs228 + clhs230 + 0.5854102*clhs231 + clhs233 + clhs235 + clhs237 + 0.5854102*clhs238 + clhs239 - 0.5854102*clhs241 + clhs244 + clhs247 + clhs250 + clhs251 + clhs252 + clhs253 + clhs254 + clhs37 + clhs59 + clhs63 + clhs99;
            lhs(0,3)=clhs101 + clhs116 + 0.5854102*clhs117 + clhs121 + 0.5854102*clhs122 + clhs124 + clhs125 + clhs22 + 0.5854102*clhs260 + 0.5854102*clhs261 + 0.5854102*clhs262 + clhs276 + clhs277 + clhs278 + clhs279 + clhs280 + clhs281 + clhs282 + clhs283 + clhs284 + clhs285 - 0.5854102*clhs286 + clhs288 + clhs290 + clhs292 + clhs37 + clhs59 + clhs61 + clhs99;
            lhs(1,0)=clhs101 + clhs103 - 0.5854102*clhs106 + clhs110 + clhs113 + clhs126 + clhs127 + clhs128 + clhs129 + 0.5854102*clhs13 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs139 + clhs141 + clhs143 + clhs145 + clhs146 + clhs147 + clhs148 + clhs150 + clhs151 + clhs152 + clhs153 + clhs155 + clhs156 + clhs157 + clhs158 + clhs159 + clhs160 + clhs182 + clhs183 + 0.5854102*clhs19 + clhs195 + clhs196 + clhs197 + clhs198 + clhs293 + clhs294 + clhs295 + clhs296 + 0.5854102*clhs297 + clhs299 + clhs301 + clhs303 + 0.5854102*clhs304 + clhs306 + clhs308 + clhs310 + clhs61 + clhs63 + clhs68 + 0.5854102*clhs7 + clhs71 + clhs79 + clhs82 + clhs90 + clhs93;
            lhs(1,1)=clhs101 + clhs103 + clhs162 + clhs164 + 0.5854102*clhs165 + clhs169 + clhs171 + 0.5854102*clhs172 + clhs176 + clhs178 + 0.5854102*clhs179 + pow(clhs184, 2)*tau[0] + pow(clhs186, 2)*tau[1] - 0.5854102*clhs187 + pow(clhs189, 2)*tau[2] + clhs191 + pow(clhs192, 2)*tau[3] + clhs194 + 0.5854102*clhs298 + clhs301 + clhs303 + 0.5854102*clhs305 + clhs308 + clhs310 + clhs311*clhs51 + clhs312*clhs51 + clhs313*clhs51 - clhs314*tau[0] - clhs314*tau[1] - clhs314*tau[2] - clhs314*tau[3] - clhs315*tau[0] - clhs315*tau[1] - clhs315*tau[2] - clhs315*tau[3] - clhs316*tau[0] - clhs316*tau[1] - clhs316*tau[2] - clhs316*tau[3] + clhs317 + clhs318 + clhs319 + clhs320 + clhs321 + clhs322 + clhs323 + clhs324 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49 - 0.34270510226404*clhs58 + clhs61 + clhs63 - 0.34270510226404*clhs98;
            lhs(1,2)=clhs103 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs160 + clhs183 + clhs217 + clhs219 + clhs221 + 0.5854102*clhs222 + clhs226 + clhs228 + 0.5854102*clhs229 + clhs233 + clhs235 + 0.5854102*clhs236 + clhs239 - 0.5854102*clhs243 + clhs247 + clhs250 + clhs299 + 0.5854102*clhs300 + clhs303 + clhs306 + 0.5854102*clhs307 + clhs310 + clhs317 + clhs321 + clhs323 + clhs324 + clhs326 + clhs328 + clhs330 + clhs332 + clhs333 + clhs334 + clhs335 + clhs337 + clhs338 + clhs339 + clhs340 + clhs342 + clhs343 + clhs344 + clhs345 + clhs346 + clhs347 + clhs348 + clhs349 + clhs350 + clhs351 + clhs352 + clhs353 + clhs63;
            lhs(1,3)=clhs101 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs160 + clhs183 + clhs258 + clhs259 + 0.5854102*clhs264 + 0.5854102*clhs265 + 0.5854102*clhs266 + clhs277 + clhs278 + clhs280 + clhs281 + clhs283 + clhs284 - 0.5854102*clhs287 + clhs290 + clhs292 + clhs299 + clhs301 + 0.5854102*clhs302 + clhs306 + clhs308 + 0.5854102*clhs309 + clhs317 + clhs321 + clhs323 + clhs324 + clhs354 + clhs355 + clhs356 + clhs358 + clhs359 + clhs360 + clhs361 + clhs363 + clhs364 + clhs365 + clhs366 + clhs368 + clhs369 + clhs370 + clhs371 + clhs372 + clhs373 + clhs374 + clhs375 + clhs376 + clhs377 + clhs378 + clhs379 + clhs61;
            lhs(2,0)=clhs103 + clhs107 - 0.5854102*clhs109 + clhs113 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs159 + clhs182 + clhs199 + clhs200 + clhs201 + clhs203 + clhs204 + clhs205 + clhs206 + clhs208 + clhs209 + clhs210 + clhs211 + clhs213 + clhs214 + clhs215 + clhs216 + clhs217 + clhs239 + clhs251 + clhs252 + clhs253 + clhs254 + clhs293 + clhs294 + clhs295 + clhs296 + 0.5854102*clhs380 + clhs382 + clhs384 + clhs386 + 0.5854102*clhs387 + clhs389 + clhs391 + clhs393 + clhs59 + clhs63 + clhs68 + 0.5854102*clhs70 + clhs72 + clhs79 + 0.5854102*clhs81 + clhs83 + clhs90 + 0.5854102*clhs92 + clhs94 + clhs99;
            lhs(2,1)=clhs103 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs160 + clhs162 + 0.5854102*clhs163 + clhs166 + clhs169 + 0.5854102*clhs170 + clhs173 + clhs176 + 0.5854102*clhs177 + clhs180 + clhs183 + clhs188 - 0.5854102*clhs190 + clhs194 + clhs217 + clhs239 + clhs317 + clhs318 + clhs319 + clhs320 + clhs321 + clhs322 + clhs326 + clhs328 + clhs330 + clhs332 + clhs333 + clhs334 + clhs335 + clhs337 + clhs338 + clhs339 + clhs340 + clhs342 + clhs343 + clhs344 + clhs345 + clhs350 + clhs351 + clhs352 + clhs353 + 0.5854102*clhs381 + clhs384 + clhs386 + 0.5854102*clhs388 + clhs391 + clhs393 + clhs394 + clhs395 + clhs63;
            lhs(2,2)=-0.34270510226404*clhs100 + clhs103 + clhs219 + 0.5854102*clhs220 + clhs223 + clhs226 + 0.5854102*clhs227 + clhs230 + clhs233 + 0.5854102*clhs234 + clhs237 + pow(clhs240, 2)*tau[0] + pow(clhs242, 2)*tau[1] + clhs244 + pow(clhs245, 2)*tau[2] - 0.5854102*clhs246 + pow(clhs248, 2)*tau[3] + clhs250 + clhs317 + clhs321 + clhs346 + clhs347 + clhs348 + clhs349 + clhs38 + clhs382 + 0.5854102*clhs383 + clhs386 + clhs389 + clhs39 + 0.5854102*clhs390 + clhs393 + clhs394 + clhs395 + clhs396*clhs51 + clhs397*clhs51 + clhs398*clhs51 - clhs399*tau[0] - clhs399*tau[1] - clhs399*tau[2] - clhs399*tau[3] + clhs40 - clhs400*tau[0] - clhs400*tau[1] - clhs400*tau[2] - clhs400*tau[3] - clhs401*tau[0] - clhs401*tau[1] - clhs401*tau[2] - clhs401*tau[3] + clhs41 + clhs42 + clhs43 + clhs44 + clhs45 + clhs46 + clhs47 + clhs48 + clhs49 + clhs59 - 0.34270510226404*clhs60 + clhs63 + clhs99;
            lhs(2,3)=clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs217 + clhs239 + clhs258 + clhs259 + 0.5854102*clhs268 + 0.5854102*clhs269 + 0.5854102*clhs270 + clhs277 + clhs279 + clhs280 + clhs282 + clhs283 + clhs285 + clhs288 - 0.5854102*clhs289 + clhs292 + clhs317 + clhs321 + clhs372 + clhs373 + clhs374 + clhs375 + clhs382 + clhs384 + 0.5854102*clhs385 + clhs389 + clhs391 + 0.5854102*clhs392 + clhs394 + clhs395 + clhs402 + clhs403 + clhs404 + clhs406 + clhs407 + clhs408 + clhs409 + clhs411 + clhs412 + clhs413 + clhs414 + clhs416 + clhs417 + clhs418 + clhs419 + clhs420 + clhs421 + clhs422 + clhs423 + clhs59 + clhs99;
            lhs(3,0)=clhs101 + clhs107 + clhs110 - 0.5854102*clhs112 + clhs276 + clhs293 + clhs294 + clhs295 + clhs296 + 0.5854102*clhs424 + clhs426 + clhs428 + clhs430 + 0.5854102*clhs431 + clhs433 + clhs435 + clhs437 + clhs59 + clhs61 + 0.5854102*clhs67 + clhs71 + clhs72 + 0.5854102*clhs78 + clhs82 + clhs83 + 0.5854102*clhs89 + clhs93 + clhs94 + clhs99;
            lhs(3,1)=clhs101 + clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs160 + 0.5854102*clhs161 + clhs164 + clhs166 + 0.5854102*clhs168 + clhs171 + clhs173 + 0.5854102*clhs175 + clhs178 + clhs180 + clhs183 + clhs188 + clhs191 - 0.5854102*clhs193 + clhs258 + clhs259 + clhs317 + clhs318 + clhs319 + clhs320 + clhs321 + clhs322 + clhs354 + clhs355 + clhs356 + clhs358 + clhs359 + clhs360 + clhs361 + clhs363 + clhs364 + clhs365 + clhs366 + clhs368 + clhs369 + clhs370 + clhs371 + clhs376 + clhs377 + clhs378 + clhs379 + 0.5854102*clhs425 + clhs428 + clhs430 + 0.5854102*clhs432 + clhs435 + clhs437 + clhs438 + clhs439 + clhs61;
            lhs(3,2)=clhs126 + clhs127 + clhs128 + clhs129 + clhs130 + clhs131 + clhs132 + clhs133 + clhs134 + clhs135 + clhs136 + clhs137 + clhs217 + 0.5854102*clhs218 + clhs221 + clhs223 + 0.5854102*clhs225 + clhs228 + clhs230 + 0.5854102*clhs232 + clhs235 + clhs237 + clhs239 + clhs244 + clhs247 - 0.5854102*clhs249 + clhs258 + clhs259 + clhs317 + clhs321 + clhs346 + clhs347 + clhs348 + clhs349 + clhs402 + clhs403 + clhs404 + clhs406 + clhs407 + clhs408 + clhs409 + clhs411 + clhs412 + clhs413 + clhs414 + clhs416 + clhs417 + clhs418 + clhs419 + clhs420 + clhs421 + clhs422 + clhs423 + clhs426 + 0.5854102*clhs427 + clhs430 + clhs433 + 0.5854102*clhs434 + clhs437 + clhs438 + clhs439 + clhs59 + clhs99;
            lhs(3,3)=clhs101 - 0.34270510226404*clhs102 + pow(clhs263, 2)*tau[0] + pow(clhs267, 2)*tau[1] + pow(clhs271, 2)*tau[2] + 0.5854102*clhs272 + 0.5854102*clhs273 + 0.5854102*clhs274 + pow(clhs275, 2)*tau[3] + clhs278 + clhs279 + clhs281 + clhs282 + clhs284 + clhs285 + clhs288 + clhs290 - 0.5854102*clhs291 + clhs317 + clhs321 + clhs372 + clhs373 + clhs374 + clhs375 + clhs38 + clhs39 + clhs40 + clhs41 + clhs42 + clhs426 + clhs428 + 0.5854102*clhs429 + clhs43 + clhs433 + clhs435 + 0.5854102*clhs436 + clhs438 + clhs439 + clhs44 + clhs440*clhs51 + clhs441*clhs51 + clhs442*clhs51 - clhs443*tau[0] - clhs443*tau[1] - clhs443*tau[2] - clhs443*tau[3] - clhs444*tau[0] - clhs444*tau[1] - clhs444*tau[2] - clhs444*tau[3] - clhs445*tau[0] - clhs445*tau[1] - clhs445*tau[2] - clhs445*tau[3] + clhs45 + clhs46 + clhs47 + clhs48 + clhs49 + clhs59 + clhs61 - 0.34270510226404*clhs62 + clhs99;


    const double crhs0 =             0.19999999899376*f[1];
const double crhs1 =             phi_subscale_gauss[1]*tau[1];
const double crhs2 =             -0.1381966*crhs1;
const double crhs3 =             1.0/delta_time;
const double crhs4 =             crhs3*phi_subscale_gauss[1];
const double crhs5 =             0.1381966*crhs4;
const double crhs6 =             0.1381966*f[0];
const double crhs7 =             0.1381966*f[2];
const double crhs8 =             0.1381966*f[3];
const double crhs9 =             crhs7 + crhs8;
const double crhs10 =             tau[1]*(crhs6 + crhs9 + 0.5854102*f[1]);
const double crhs11 =             -0.1381966*crhs10;
const double crhs12 =             0.1381966*prj[0];
const double crhs13 =             0.1381966*prj[2];
const double crhs14 =             0.1381966*prj[3];
const double crhs15 =             crhs13 + crhs14;
const double crhs16 =             tau[1]*(crhs12 + crhs15 + 0.5854102*prj[1]);
const double crhs17 =             -0.1381966*crhs16;
const double crhs18 =             0.1381966*v(0,0);
const double crhs19 =             0.1381966*v(2,0);
const double crhs20 =             0.1381966*v(3,0);
const double crhs21 =             crhs19 + crhs20;
const double crhs22 =             crhs18 + crhs21 + 0.5854102*v(1,0);
const double crhs23 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs24 =             crhs22*crhs23;
const double crhs25 =             -0.1381966*crhs24;
const double crhs26 =             0.1381966*v(0,1);
const double crhs27 =             0.1381966*v(2,1);
const double crhs28 =             0.1381966*v(3,1);
const double crhs29 =             crhs27 + crhs28;
const double crhs30 =             crhs26 + crhs29 + 0.5854102*v(1,1);
const double crhs31 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs32 =             crhs30*crhs31;
const double crhs33 =             -0.1381966*crhs32;
const double crhs34 =             0.1381966*v(0,2);
const double crhs35 =             0.1381966*v(2,2);
const double crhs36 =             0.1381966*v(3,2);
const double crhs37 =             crhs35 + crhs36;
const double crhs38 =             crhs34 + crhs37 + 0.5854102*v(1,2);
const double crhs39 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs40 =             crhs38*crhs39;
const double crhs41 =             -0.1381966*crhs40;
const double crhs42 =             1.0/RK_time_coefficient;
const double crhs43 =             crhs3*crhs42*tau[1];
const double crhs44 =             0.08090169924532*phi[1];
const double crhs45 =             -0.08090169924532*phi_old[1];
const double crhs46 =             0.01909830025156*phi[0];
const double crhs47 =             0.01909830025156*phi[2];
const double crhs48 =             0.01909830025156*phi[3];
const double crhs49 =             -0.01909830025156*phi_old[0];
const double crhs50 =             -0.01909830025156*phi_old[2];
const double crhs51 =             -0.01909830025156*phi_old[3];
const double crhs52 =             crhs43*(crhs44 + crhs45 + crhs46 + crhs47 + crhs48 + crhs49 + crhs50 + crhs51);
const double crhs53 =             0.1381966*phi[0];
const double crhs54 =             0.5854102*phi[1];
const double crhs55 =             0.1381966*phi[2];
const double crhs56 =             0.1381966*phi[3];
const double crhs57 =             crhs55 + crhs56;
const double crhs58 =             crhs53 + crhs54 + crhs57;
const double crhs59 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs60 =             crhs58*crhs59;
const double crhs61 =             0.1381966*crhs60;
const double crhs62 =             -crhs61;
const double crhs63 =             crhs61*tau[1];
const double crhs64 =             tau[1]*(crhs24 + crhs32 + crhs40);
const double crhs65 =             0.1381966*crhs64;
const double crhs66 =             0.19999999899376*f[2];
const double crhs67 =             0.19999999899376*f[3];
const double crhs68 =             phi_subscale_gauss[0]*tau[0];
const double crhs69 =             phi_subscale_gauss[2]*tau[2];
const double crhs70 =             -0.1381966*crhs69;
const double crhs71 =             phi_subscale_gauss[3]*tau[3];
const double crhs72 =             -0.1381966*crhs71;
const double crhs73 =             crhs3*phi_subscale_gauss[0];
const double crhs74 =             crhs3*phi_subscale_gauss[2];
const double crhs75 =             0.1381966*crhs74;
const double crhs76 =             crhs3*phi_subscale_gauss[3];
const double crhs77 =             0.1381966*crhs76;
const double crhs78 =             0.1381966*f[1];
const double crhs79 =             tau[0]*(crhs78 + crhs9 + 0.5854102*f[0]);
const double crhs80 =             0.1381966*prj[1];
const double crhs81 =             tau[0]*(crhs15 + crhs80 + 0.5854102*prj[0]);
const double crhs82 =             crhs6 + crhs78;
const double crhs83 =             tau[2]*(crhs8 + crhs82 + 0.5854102*f[2]);
const double crhs84 =             -0.1381966*crhs83;
const double crhs85 =             crhs12 + crhs80;
const double crhs86 =             tau[2]*(crhs14 + crhs85 + 0.5854102*prj[2]);
const double crhs87 =             -0.1381966*crhs86;
const double crhs88 =             tau[3]*(crhs7 + crhs82 + 0.5854102*f[3]);
const double crhs89 =             -0.1381966*crhs88;
const double crhs90 =             tau[3]*(crhs13 + crhs85 + 0.5854102*prj[3]);
const double crhs91 =             -0.1381966*crhs90;
const double crhs92 =             4*crhs23*k;
const double crhs93 =             4*crhs31*k;
const double crhs94 =             4*crhs39*k;
const double crhs95 =             DN_DX_0_0*crhs23*k;
const double crhs96 =             DN_DX_0_1*crhs31*k;
const double crhs97 =             DN_DX_0_2*crhs39*k;
const double crhs98 =             0.1381966*v(1,0);
const double crhs99 =             crhs18 + crhs98;
const double crhs100 =             crhs19 + crhs99 + 0.5854102*v(3,0);
const double crhs101 =             crhs100*crhs23;
const double crhs102 =             -0.1381966*crhs101;
const double crhs103 =             crhs20 + crhs99 + 0.5854102*v(2,0);
const double crhs104 =             crhs103*crhs23;
const double crhs105 =             -0.1381966*crhs104;
const double crhs106 =             crhs21 + crhs98 + 0.5854102*v(0,0);
const double crhs107 =             crhs106*crhs23;
const double crhs108 =             0.1381966*v(1,1);
const double crhs109 =             crhs108 + crhs26;
const double crhs110 =             crhs109 + crhs27 + 0.5854102*v(3,1);
const double crhs111 =             crhs110*crhs31;
const double crhs112 =             -0.1381966*crhs111;
const double crhs113 =             crhs109 + crhs28 + 0.5854102*v(2,1);
const double crhs114 =             crhs113*crhs31;
const double crhs115 =             -0.1381966*crhs114;
const double crhs116 =             crhs108 + crhs29 + 0.5854102*v(0,1);
const double crhs117 =             crhs116*crhs31;
const double crhs118 =             0.1381966*v(1,2);
const double crhs119 =             crhs118 + crhs34;
const double crhs120 =             crhs119 + crhs35 + 0.5854102*v(3,2);
const double crhs121 =             crhs120*crhs39;
const double crhs122 =             -0.1381966*crhs121;
const double crhs123 =             crhs119 + crhs36 + 0.5854102*v(2,2);
const double crhs124 =             crhs123*crhs39;
const double crhs125 =             -0.1381966*crhs124;
const double crhs126 =             crhs118 + crhs37 + 0.5854102*v(0,2);
const double crhs127 =             crhs126*crhs39;
const double crhs128 =             crhs3*crhs42*tau[0];
const double crhs129 =             0.08090169924532*phi[2];
const double crhs130 =             0.08090169924532*phi[3];
const double crhs131 =             -0.08090169924532*phi_old[2];
const double crhs132 =             -0.08090169924532*phi_old[3];
const double crhs133 =             crhs129 + crhs130 + crhs131 + crhs132;
const double crhs134 =             crhs3*crhs42*tau[2];
const double crhs135 =             0.01909830025156*phi[1];
const double crhs136 =             -0.01909830025156*phi_old[1];
const double crhs137 =             crhs134*(crhs129 + crhs131 + crhs135 + crhs136 + crhs46 + crhs48 + crhs49 + crhs51);
const double crhs138 =             crhs3*crhs42*tau[3];
const double crhs139 =             crhs138*(crhs130 + crhs132 + crhs135 + crhs136 + crhs46 + crhs47 + crhs49 + crhs50);
const double crhs140 =             0.1381966*phi[1];
const double crhs141 =             crhs140 + crhs53;
const double crhs142 =             0.5854102*phi[3];
const double crhs143 =             crhs141 + crhs142 + crhs55;
const double crhs144 =             crhs143*crhs59;
const double crhs145 =             0.1381966*crhs144;
const double crhs146 =             -crhs145;
const double crhs147 =             0.5854102*phi[2];
const double crhs148 =             crhs141 + crhs147 + crhs56;
const double crhs149 =             crhs148*crhs59;
const double crhs150 =             0.1381966*crhs149;
const double crhs151 =             -crhs150;
const double crhs152 =             0.5854102*phi[0];
const double crhs153 =             crhs140 + crhs152 + crhs57;
const double crhs154 =             crhs153*crhs59;
const double crhs155 =             0.5854102*crhs154;
const double crhs156 =             crhs3*phi_subscale_gauss[0]*tau[0];
const double crhs157 =             DN_DX_0_0*crhs106 + DN_DX_0_1*crhs116 + DN_DX_0_2*crhs126;
const double crhs158 =             crhs3*phi_subscale_gauss[1]*tau[1];
const double crhs159 =             DN_DX_0_0*crhs22 + DN_DX_0_1*crhs30 + DN_DX_0_2*crhs38;
const double crhs160 =             crhs3*phi_subscale_gauss[2]*tau[2];
const double crhs161 =             DN_DX_0_0*crhs103 + DN_DX_0_1*crhs113 + DN_DX_0_2*crhs123;
const double crhs162 =             crhs3*phi_subscale_gauss[3]*tau[3];
const double crhs163 =             DN_DX_0_0*crhs100 + DN_DX_0_1*crhs110 + DN_DX_0_2*crhs120;
const double crhs164 =             crhs150*tau[2];
const double crhs165 =             crhs145*tau[3];
const double crhs166 =             -0.1381966*phi_old[2];
const double crhs167 =             -0.1381966*phi_old[3];
const double crhs168 =             crhs166 + crhs167;
const double crhs169 =             -0.1381966*phi_old[1];
const double crhs170 =             crhs3*crhs42*tau[0]*(crhs140 + crhs152 + crhs168 + crhs169 + crhs55 + crhs56 - 0.5854102*phi_old[0]);
const double crhs171 =             -0.1381966*phi_old[0];
const double crhs172 =             crhs3*crhs42*tau[1]*(crhs168 + crhs171 + crhs53 + crhs54 + crhs55 + crhs56 - 0.5854102*phi_old[1]);
const double crhs173 =             crhs169 + crhs171;
const double crhs174 =             crhs3*crhs42*tau[2]*(crhs140 + crhs147 + crhs167 + crhs173 + crhs53 + crhs56 - 0.5854102*phi_old[2]);
const double crhs175 =             crhs3*crhs42*tau[3]*(crhs140 + crhs142 + crhs166 + crhs173 + crhs53 + crhs55 - 0.5854102*phi_old[3]);
const double crhs176 =             tau[0]*(crhs107 + crhs117 + crhs127);
const double crhs177 =             tau[2]*(crhs104 + crhs114 + crhs124);
const double crhs178 =             0.1381966*crhs177;
const double crhs179 =             tau[3]*(crhs101 + crhs111 + crhs121);
const double crhs180 =             0.1381966*crhs179;
const double crhs181 =             crhs153*crhs59*tau[0];
const double crhs182 =             crhs58*crhs59*tau[1];
const double crhs183 =             crhs148*crhs59*tau[2];
const double crhs184 =             crhs143*crhs59*tau[3];
const double crhs185 =             0.19999999899376*f[0];
const double crhs186 =             -0.1381966*crhs68;
const double crhs187 =             0.1381966*crhs73;
const double crhs188 =             -0.1381966*crhs79;
const double crhs189 =             -0.1381966*crhs81;
const double crhs190 =             -0.1381966*crhs107;
const double crhs191 =             -0.1381966*crhs117;
const double crhs192 =             -0.1381966*crhs127;
const double crhs193 =             0.08090169924532*phi[0];
const double crhs194 =             -0.08090169924532*phi_old[0];
const double crhs195 =             crhs128*(crhs135 + crhs136 + crhs193 + crhs194 + crhs47 + crhs48 + crhs50 + crhs51);
const double crhs196 =             0.1381966*crhs154;
const double crhs197 =             -crhs196;
const double crhs198 =             crhs196*tau[0];
const double crhs199 =             0.1381966*crhs176;
const double crhs200 =             DN_DX_1_0*crhs23*k;
const double crhs201 =             DN_DX_1_1*crhs31*k;
const double crhs202 =             DN_DX_1_2*crhs39*k;
const double crhs203 =             0.5854102*crhs60;
const double crhs204 =             DN_DX_1_0*crhs106 + DN_DX_1_1*crhs116 + DN_DX_1_2*crhs126;
const double crhs205 =             DN_DX_1_0*crhs22 + DN_DX_1_1*crhs30 + DN_DX_1_2*crhs38;
const double crhs206 =             DN_DX_1_0*crhs103 + DN_DX_1_1*crhs113 + DN_DX_1_2*crhs123;
const double crhs207 =             DN_DX_1_0*crhs100 + DN_DX_1_1*crhs110 + DN_DX_1_2*crhs120;
const double crhs208 =             crhs0 + crhs11 + crhs17 + crhs185 + crhs186 + crhs187 + crhs188 + crhs189 + crhs190 + crhs191 + crhs192 + crhs195 + crhs197 + crhs198 + crhs199 + crhs2 + crhs25 + crhs33 + crhs41 + crhs5 + crhs52 + crhs62 + crhs63 + crhs65;
const double crhs209 =             DN_DX_2_0*crhs23*k;
const double crhs210 =             DN_DX_2_1*crhs31*k;
const double crhs211 =             DN_DX_2_2*crhs39*k;
const double crhs212 =             crhs193 + crhs194 + crhs44 + crhs45;
const double crhs213 =             0.5854102*crhs149;
const double crhs214 =             DN_DX_2_0*crhs106 + DN_DX_2_1*crhs116 + DN_DX_2_2*crhs126;
const double crhs215 =             DN_DX_2_0*crhs22 + DN_DX_2_1*crhs30 + DN_DX_2_2*crhs38;
const double crhs216 =             DN_DX_2_0*crhs103 + DN_DX_2_1*crhs113 + DN_DX_2_2*crhs123;
const double crhs217 =             DN_DX_2_0*crhs100 + DN_DX_2_1*crhs110 + DN_DX_2_2*crhs120;
const double crhs218 =             DN_DX_3_0*crhs23*k;
const double crhs219 =             DN_DX_3_1*crhs31*k;
const double crhs220 =             DN_DX_3_2*crhs39*k;
const double crhs221 =             0.5854102*crhs144;
const double crhs222 =             DN_DX_3_0*crhs106 + DN_DX_3_1*crhs116 + DN_DX_3_2*crhs126;
const double crhs223 =             DN_DX_3_0*crhs22 + DN_DX_3_1*crhs30 + DN_DX_3_2*crhs38;
const double crhs224 =             DN_DX_3_0*crhs103 + DN_DX_3_1*crhs113 + DN_DX_3_2*crhs123;
const double crhs225 =             DN_DX_3_0*crhs100 + DN_DX_3_1*crhs110 + DN_DX_3_2*crhs120;
            rhs[0]=-DN_DX_0_0*crhs92 - DN_DX_0_1*crhs93 - DN_DX_0_2*crhs94 + crhs0 + crhs10*crhs159 + crhs102 + crhs105 - 0.5854102*crhs107 + crhs11 + crhs112 + crhs115 - 0.5854102*crhs117 + crhs122 + crhs125 - 0.5854102*crhs127 + crhs128*(crhs133 + crhs44 + crhs45 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs137 + crhs139 + crhs146 + crhs151 + crhs155*tau[0] - crhs155 + crhs156*crhs157 - crhs157*crhs170 - crhs157*crhs176 - crhs157*crhs181 + crhs157*crhs79 + crhs157*crhs81 + crhs158*crhs159 + crhs159*crhs16 - crhs159*crhs172 - crhs159*crhs182 - crhs159*crhs64 + crhs160*crhs161 - crhs161*crhs174 - crhs161*crhs177 - crhs161*crhs183 + crhs161*crhs83 + crhs161*crhs86 + crhs162*crhs163 - crhs163*crhs175 - crhs163*crhs179 - crhs163*crhs184 + crhs163*crhs88 + crhs163*crhs90 + crhs164 + crhs165 + crhs17 + 0.5854102*crhs176 + crhs178 + crhs180 + crhs2 + crhs25 + crhs33 + crhs41 + crhs5 + crhs52 + crhs62 + crhs63 + crhs65 + crhs66 + crhs67 - 0.5854102*crhs68 + crhs70 + crhs72 + 0.5854102*crhs73 + crhs75 + crhs77 - 0.5854102*crhs79 - 0.5854102*crhs81 + crhs84 + crhs87 + crhs89 + crhs91 + crhs95*tau[0] + crhs95*tau[1] + crhs95*tau[2] + crhs95*tau[3] + crhs96*tau[0] + crhs96*tau[1] + crhs96*tau[2] + crhs96*tau[3] + crhs97*tau[0] + crhs97*tau[1] + crhs97*tau[2] + crhs97*tau[3] + 0.40000000301872*f[0];
            rhs[1]=-DN_DX_1_0*crhs92 - DN_DX_1_1*crhs93 - DN_DX_1_2*crhs94 - 0.5854102*crhs1 + crhs10*crhs205 - 0.5854102*crhs10 + crhs102 + crhs105 + crhs112 + crhs115 + crhs122 + crhs125 + crhs137 + crhs139 + crhs146 + crhs151 + crhs156*crhs204 + crhs158*crhs205 + crhs16*crhs205 - 0.5854102*crhs16 + crhs160*crhs206 + crhs162*crhs207 + crhs164 + crhs165 - crhs170*crhs204 - crhs172*crhs205 - crhs174*crhs206 - crhs175*crhs207 - crhs176*crhs204 - crhs177*crhs206 + crhs178 - crhs179*crhs207 + crhs180 - crhs181*crhs204 - crhs182*crhs205 - crhs183*crhs206 - crhs184*crhs207 + crhs185 + crhs186 + crhs187 + crhs188 + crhs189 + crhs190 + crhs191 + crhs192 + crhs195 + crhs197 + crhs198 + crhs199 + crhs200*tau[0] + crhs200*tau[1] + crhs200*tau[2] + crhs200*tau[3] + crhs201*tau[0] + crhs201*tau[1] + crhs201*tau[2] + crhs201*tau[3] + crhs202*tau[0] + crhs202*tau[1] + crhs202*tau[2] + crhs202*tau[3] + crhs203*tau[1] - crhs203 + crhs204*crhs79 + crhs204*crhs81 - crhs205*crhs64 + crhs206*crhs83 + crhs206*crhs86 + crhs207*crhs88 + crhs207*crhs90 - 0.5854102*crhs24 - 0.5854102*crhs32 + 0.5854102*crhs4 - 0.5854102*crhs40 + crhs43*(crhs133 + crhs193 + crhs194 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) + 0.5854102*crhs64 + crhs66 + crhs67 + crhs70 + crhs72 + crhs75 + crhs77 + crhs84 + crhs87 + crhs89 + crhs91 + 0.40000000301872*f[1];
            rhs[2]=-DN_DX_2_0*crhs92 - DN_DX_2_1*crhs93 - DN_DX_2_2*crhs94 + crhs10*crhs215 + crhs102 - 0.5854102*crhs104 + crhs112 - 0.5854102*crhs114 + crhs122 - 0.5854102*crhs124 + crhs134*(crhs130 + crhs132 + crhs212 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) + crhs139 + crhs146 + crhs156*crhs214 + crhs158*crhs215 + crhs16*crhs215 + crhs160*crhs216 + crhs162*crhs217 + crhs165 - crhs170*crhs214 - crhs172*crhs215 - crhs174*crhs216 - crhs175*crhs217 - crhs176*crhs214 - crhs177*crhs216 + 0.5854102*crhs177 - crhs179*crhs217 + crhs180 - crhs181*crhs214 - crhs182*crhs215 - crhs183*crhs216 - crhs184*crhs217 + crhs208 + crhs209*tau[0] + crhs209*tau[1] + crhs209*tau[2] + crhs209*tau[3] + crhs210*tau[0] + crhs210*tau[1] + crhs210*tau[2] + crhs210*tau[3] + crhs211*tau[0] + crhs211*tau[1] + crhs211*tau[2] + crhs211*tau[3] + crhs213*tau[2] - crhs213 + crhs214*crhs79 + crhs214*crhs81 - crhs215*crhs64 + crhs216*crhs83 + crhs216*crhs86 + crhs217*crhs88 + crhs217*crhs90 - 0.5854102*crhs69 + crhs72 + 0.5854102*crhs74 + crhs77 - 0.5854102*crhs83 - 0.5854102*crhs86 + crhs89 + crhs91 + 0.40000000301872*f[2] + 0.19999999899376*f[3];
            rhs[3]=-DN_DX_3_0*crhs92 - DN_DX_3_1*crhs93 - DN_DX_3_2*crhs94 + crhs10*crhs223 - 0.5854102*crhs101 + crhs105 - 0.5854102*crhs111 + crhs115 - 0.5854102*crhs121 + crhs125 + crhs137 + crhs138*(crhs129 + crhs131 + crhs212 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) + crhs151 + crhs156*crhs222 + crhs158*crhs223 + crhs16*crhs223 + crhs160*crhs224 + crhs162*crhs225 + crhs164 - crhs170*crhs222 - crhs172*crhs223 - crhs174*crhs224 - crhs175*crhs225 - crhs176*crhs222 - crhs177*crhs224 + crhs178 - crhs179*crhs225 + 0.5854102*crhs179 - crhs181*crhs222 - crhs182*crhs223 - crhs183*crhs224 - crhs184*crhs225 + crhs208 + crhs218*tau[0] + crhs218*tau[1] + crhs218*tau[2] + crhs218*tau[3] + crhs219*tau[0] + crhs219*tau[1] + crhs219*tau[2] + crhs219*tau[3] + crhs220*tau[0] + crhs220*tau[1] + crhs220*tau[2] + crhs220*tau[3] + crhs221*tau[3] - crhs221 + crhs222*crhs79 + crhs222*crhs81 - crhs223*crhs64 + crhs224*crhs83 + crhs224*crhs86 + crhs225*crhs88 + crhs225*crhs90 + crhs70 - 0.5854102*crhs71 + crhs75 + 0.5854102*crhs76 + crhs84 + crhs87 - 0.5854102*crhs88 - 0.5854102*crhs90 + 0.19999999899376*f[2] + 0.40000000301872*f[3];


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
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
    ElementVariables& rVariables)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicDynamicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicDynamicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
