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
void SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    BoundedVector<double, TNumNodes> rhs;
    this->CalculateRightHandSideInternal(rhs,rCurrentProcessInfo);
    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    const auto& reaction_variable = r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable();
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(reaction_variable) += rhs[i_node];
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
    // Resize and initialize dynamic subscales container
    if (mUnknownSubScale.size() != TNumNodes)
        mUnknownSubScale.resize(TNumNodes, false); // number integration points = number nodes
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
    const auto& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

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

    if (rVariable == SCALAR_PROJECTION) {
        auto& r_geometry = this->GetGeometry();
        const unsigned int local_size = r_geometry.size();
        BoundedVector<double, TNumNodes> rhs_oss;
        this->CalculateOrthogonalSubgridScaleRHSInternal(rhs_oss,rCurrentProcessInfo);
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
        }
    }
    else {
        KRATOS_ERROR << "Variable not implemented to compute OSS projection. Use SCALAR_PROJECTION instead." << std::endl;
    }

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
void SymbolicDConvectionDiffusionExplicit<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rVariables);

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
    // RHS
    auto& rhs = rVariables.rhs;

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
    noalias(rRightHandSideVector) = rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3,4>::CalculateRightHandSideInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rVariables);

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
    // RHS
    auto& rhs = rVariables.rhs;

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
    noalias(rRightHandSideVector) = rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<2,3>::CalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rVariables);

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
    noalias(rRightHandSideVector) = rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3,4>::CalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rVariables);

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
    noalias(rRightHandSideVector) = rhs * rVariables.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<2,3>::UpdateUnknownSubgridScaleGaussPoint(
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
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicDConvectionDiffusionExplicit<3,4>::UpdateUnknownSubgridScaleGaussPoint(
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
    for(unsigned int g = 0; g<TNumNodes; g++) {
	const auto& N = row(rVariables.N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        const array_1d<double,3> vel_gauss = prod(N, rVariables.convective_velocity);
        double div_vel = 0;
        for(unsigned int node_element = 0; node_element<TNumNodes; node_element++) {
            for(unsigned int dim = 0; dim < TDim; dim++) {
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

template class SymbolicDConvectionDiffusionExplicit<2,3>;
template class SymbolicDConvectionDiffusionExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}
