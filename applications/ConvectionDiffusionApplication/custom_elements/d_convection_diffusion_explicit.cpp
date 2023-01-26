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

#include "d_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::DConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : QSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::DConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : QSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::~DConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer DConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer DConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit Convection-Diffusion element.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit Convection-Diffusion element. Call the DCalculateRightHandSideInternal() method instead.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    BoundedVector<double, TNumNodes> rhs;
    this->DCalculateRightHandSideInternal(rhs,rCurrentProcessInfo);
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
void DConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);
    mUnknownSubScale = ZeroVector(TNumNodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::FinalizeSolutionStep(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);
    // Reading integration points and local gradients
    const auto& r_geometry = this->GetGeometry();
    const auto& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Iterate over integration points to update subscales on the gauss point
    for (unsigned int g = 0; g < integration_points.size(); g++) {
        // Caluclate N on the gauss point "g"
        rData.N = row(rData.N_gausspoint,g);
        // Compute tau
        this->DCalculateTau(rData);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rData.unknown_subscale = mUnknownSubScale(g);
        // Update unknown belonging to subgrid scale space on gauss integration point g
        this->UpdateUnknownSubgridScaleGaussPoint(rData,g);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    if (rVariable == p_settings->GetProjectionVariable()) {
        auto& r_geometry = this->GetGeometry();
        const unsigned int local_size = r_geometry.size();
        BoundedVector<double, TNumNodes> rhs_oss;
        this->DCalculateOrthogonalSubgridScaleRHSInternal(rhs_oss,rCurrentProcessInfo);
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
        }
    }
    else {
        BaseType::Calculate(rVariable,Output,rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::DCalculateRightHandSideInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

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
const double crhs26 =             explicit_step_coefficient*tau[1];
const double crhs27 =             0.111111111111111*phi[1];
const double crhs28 =             -0.111111111111111*phi_old[1];
const double crhs29 =             crhs27 + crhs28;
const double crhs30 =             0.0277777777777778*phi[0];
const double crhs31 =             0.0277777777777778*phi[2];
const double crhs32 =             -0.0277777777777778*phi_old[0];
const double crhs33 =             -0.0277777777777778*phi_old[2];
const double crhs34 =             crhs26*(crhs29 + crhs30 + crhs31 + crhs32 + crhs33);
const double crhs35 =             0.166666666666667*phi[0];
const double crhs36 =             0.166666666666667*phi[2];
const double crhs37 =             crhs35 + crhs36 + 0.666666666666667*phi[1];
const double crhs38 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs39 =             crhs37*crhs38;
const double crhs40 =             0.166666666666667*crhs39;
const double crhs41 =             -crhs40;
const double crhs42 =             crhs40*tau[1];
const double crhs43 =             tau[1]*(crhs18 + crhs24);
const double crhs44 =             0.166666666666667*crhs43;
const double crhs45 =             0.25*f[2];
const double crhs46 =             phi_subscale_gauss[0]*tau[0];
const double crhs47 =             phi_subscale_gauss[2]*tau[2];
const double crhs48 =             -0.166666666666667*crhs47;
const double crhs49 =             crhs3*phi_subscale_gauss[0];
const double crhs50 =             crhs3*phi_subscale_gauss[2];
const double crhs51 =             0.166666666666667*crhs50;
const double crhs52 =             0.166666666666667*f[1];
const double crhs53 =             tau[0]*(crhs52 + crhs7 + 0.666666666666667*f[0]);
const double crhs54 =             0.166666666666667*prj[1];
const double crhs55 =             tau[0]*(crhs11 + crhs54 + 0.666666666666667*prj[0]);
const double crhs56 =             tau[2]*(crhs52 + crhs6 + 0.666666666666667*f[2]);
const double crhs57 =             -0.166666666666667*crhs56;
const double crhs58 =             tau[2]*(crhs10 + crhs54 + 0.666666666666667*prj[2]);
const double crhs59 =             -0.166666666666667*crhs58;
const double crhs60 =             3*alpha*crhs17;
const double crhs61 =             3*alpha*crhs23;
const double crhs62 =             DN_DX_0_0*alpha*crhs17;
const double crhs63 =             DN_DX_0_1*alpha*crhs23;
const double crhs64 =             0.166666666666667*v(1,0);
const double crhs65 =             crhs14 + crhs64 + 0.666666666666667*v(2,0);
const double crhs66 =             crhs17*crhs65;
const double crhs67 =             -0.166666666666667*crhs66;
const double crhs68 =             crhs15 + crhs64 + 0.666666666666667*v(0,0);
const double crhs69 =             crhs17*crhs68;
const double crhs70 =             0.166666666666667*v(1,1);
const double crhs71 =             crhs20 + crhs70 + 0.666666666666667*v(2,1);
const double crhs72 =             crhs23*crhs71;
const double crhs73 =             -0.166666666666667*crhs72;
const double crhs74 =             crhs21 + crhs70 + 0.666666666666667*v(0,1);
const double crhs75 =             crhs23*crhs74;
const double crhs76 =             explicit_step_coefficient*tau[0];
const double crhs77 =             0.111111111111111*phi[2];
const double crhs78 =             -0.111111111111111*phi_old[2];
const double crhs79 =             explicit_step_coefficient*tau[2];
const double crhs80 =             0.0277777777777778*phi[1];
const double crhs81 =             -0.0277777777777778*phi_old[1];
const double crhs82 =             crhs79*(crhs30 + crhs32 + crhs77 + crhs78 + crhs80 + crhs81);
const double crhs83 =             0.166666666666667*phi[1];
const double crhs84 =             crhs35 + crhs83 + 0.666666666666667*phi[2];
const double crhs85 =             crhs38*crhs84;
const double crhs86 =             0.166666666666667*crhs85;
const double crhs87 =             -crhs86;
const double crhs88 =             crhs36 + crhs83 + 0.666666666666667*phi[0];
const double crhs89 =             crhs38*crhs88;
const double crhs90 =             0.666666666666667*crhs89;
const double crhs91 =             crhs3*phi_subscale_gauss[0]*tau[0];
const double crhs92 =             DN_DX_0_0*crhs68 + DN_DX_0_1*crhs74;
const double crhs93 =             crhs3*phi_subscale_gauss[1]*tau[1];
const double crhs94 =             DN_DX_0_0*crhs16 + DN_DX_0_1*crhs22;
const double crhs95 =             crhs3*phi_subscale_gauss[2]*tau[2];
const double crhs96 =             DN_DX_0_0*crhs65 + DN_DX_0_1*crhs71;
const double crhs97 =             crhs86*tau[2];
const double crhs98 =             tau[0]*(crhs69 + crhs75);
const double crhs99 =             tau[2]*(crhs66 + crhs72);
const double crhs100 =             0.166666666666667*crhs99;
const double crhs101 =             -0.166666666666667*phi_old[1];
const double crhs102 =             -0.166666666666667*phi_old[2];
const double crhs103 =             explicit_step_coefficient*tau[0]*(crhs101 + crhs102 + crhs88 - 0.666666666666667*phi_old[0]);
const double crhs104 =             -0.166666666666667*phi_old[0];
const double crhs105 =             explicit_step_coefficient*tau[1]*(crhs102 + crhs104 + crhs37 - 0.666666666666667*phi_old[1]);
const double crhs106 =             explicit_step_coefficient*tau[2]*(crhs101 + crhs104 + crhs84 - 0.666666666666667*phi_old[2]);
const double crhs107 =             crhs38*crhs88*tau[0];
const double crhs108 =             crhs37*crhs38*tau[1];
const double crhs109 =             crhs38*crhs84*tau[2];
const double crhs110 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs111 =             0.166666666666667*crhs89;
const double crhs112 =             crhs111*tau[0] - crhs111 - 0.166666666666667*crhs46 + 0.166666666666667*crhs49 - 0.166666666666667*crhs53 - 0.166666666666667*crhs55 - 0.166666666666667*crhs69 - 0.166666666666667*crhs75 + crhs76*(crhs110 + crhs31 + crhs33 + crhs80 + crhs81) + 0.166666666666667*crhs98 + 0.25*f[0];
const double crhs113 =             DN_DX_1_0*alpha*crhs17;
const double crhs114 =             DN_DX_1_1*alpha*crhs23;
const double crhs115 =             0.666666666666667*crhs39;
const double crhs116 =             DN_DX_1_0*crhs68 + DN_DX_1_1*crhs74;
const double crhs117 =             DN_DX_1_0*crhs16 + DN_DX_1_1*crhs22;
const double crhs118 =             DN_DX_1_0*crhs65 + DN_DX_1_1*crhs71;
const double crhs119 =             DN_DX_2_0*alpha*crhs17;
const double crhs120 =             DN_DX_2_1*alpha*crhs23;
const double crhs121 =             0.666666666666667*crhs85;
const double crhs122 =             DN_DX_2_0*crhs68 + DN_DX_2_1*crhs74;
const double crhs123 =             DN_DX_2_0*crhs16 + DN_DX_2_1*crhs22;
const double crhs124 =             DN_DX_2_0*crhs65 + DN_DX_2_1*crhs71;
            rhs[0]=-DN_DX_0_0*crhs60 - DN_DX_0_1*crhs61 + crhs0 + crhs100 - crhs103*crhs92 - crhs105*crhs94 - crhs106*crhs96 - crhs107*crhs92 - crhs108*crhs94 - crhs109*crhs96 + crhs12*crhs94 + crhs13 + crhs19 + crhs2 + crhs25 + crhs34 + crhs41 + crhs42 - crhs43*crhs94 + crhs44 + crhs45 - 0.666666666666667*crhs46 + crhs48 + 0.666666666666667*crhs49 + crhs5 + crhs51 + crhs53*crhs92 - 0.666666666666667*crhs53 + crhs55*crhs92 - 0.666666666666667*crhs55 + crhs56*crhs96 + crhs57 + crhs58*crhs96 + crhs59 + crhs62*tau[0] + crhs62*tau[1] + crhs62*tau[2] + crhs63*tau[0] + crhs63*tau[1] + crhs63*tau[2] + crhs67 - 0.666666666666667*crhs69 + crhs73 - 0.666666666666667*crhs75 + crhs76*(crhs29 + crhs77 + crhs78 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs8*crhs94 + crhs82 + crhs87 + crhs9 + crhs90*tau[0] - crhs90 + crhs91*crhs92 - crhs92*crhs98 + crhs93*crhs94 + crhs95*crhs96 - crhs96*crhs99 + crhs97 + 0.666666666666667*crhs98 + 0.5*f[0];
            rhs[1]=-DN_DX_1_0*crhs60 - DN_DX_1_1*crhs61 - 0.666666666666667*crhs1 + crhs100 - crhs103*crhs116 - crhs105*crhs117 - crhs106*crhs118 - crhs107*crhs116 - crhs108*crhs117 - crhs109*crhs118 + crhs112 + crhs113*tau[0] + crhs113*tau[1] + crhs113*tau[2] + crhs114*tau[0] + crhs114*tau[1] + crhs114*tau[2] + crhs115*tau[1] - crhs115 + crhs116*crhs53 + crhs116*crhs55 + crhs116*crhs91 - crhs116*crhs98 + crhs117*crhs12 - crhs117*crhs43 + crhs117*crhs8 + crhs117*crhs93 + crhs118*crhs56 + crhs118*crhs58 + crhs118*crhs95 - crhs118*crhs99 - 0.666666666666667*crhs12 - 0.666666666666667*crhs18 - 0.666666666666667*crhs24 + crhs26*(crhs110 + crhs77 + crhs78 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) + 0.666666666666667*crhs4 + 0.666666666666667*crhs43 + crhs45 + crhs48 + crhs51 + crhs57 + crhs59 + crhs67 + crhs73 - 0.666666666666667*crhs8 + crhs82 + crhs87 + crhs97 + 0.5*f[1];
            rhs[2]=-DN_DX_2_0*crhs60 - DN_DX_2_1*crhs61 + crhs0 - crhs103*crhs122 - crhs105*crhs123 - crhs106*crhs124 - crhs107*crhs122 - crhs108*crhs123 - crhs109*crhs124 + crhs112 + crhs119*tau[0] + crhs119*tau[1] + crhs119*tau[2] + crhs12*crhs123 + crhs120*tau[0] + crhs120*tau[1] + crhs120*tau[2] + crhs121*tau[2] - crhs121 + crhs122*crhs53 + crhs122*crhs55 + crhs122*crhs91 - crhs122*crhs98 - crhs123*crhs43 + crhs123*crhs8 + crhs123*crhs93 + crhs124*crhs56 + crhs124*crhs58 + crhs124*crhs95 - crhs124*crhs99 + crhs13 + crhs19 + crhs2 + crhs25 + crhs34 + crhs41 + crhs42 + crhs44 - 0.666666666666667*crhs47 + crhs5 + 0.666666666666667*crhs50 - 0.666666666666667*crhs56 - 0.666666666666667*crhs58 - 0.666666666666667*crhs66 - 0.666666666666667*crhs72 + crhs79*(crhs110 + crhs27 + crhs28 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs9 + 0.666666666666667*crhs99 + 0.5*f[2];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<3,4>::DCalculateRightHandSideInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

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
const double crhs42 =             explicit_step_coefficient*tau[1];
const double crhs43 =             0.08090169924532*phi[1];
const double crhs44 =             -0.08090169924532*phi_old[1];
const double crhs45 =             0.01909830025156*phi[0];
const double crhs46 =             0.01909830025156*phi[2];
const double crhs47 =             0.01909830025156*phi[3];
const double crhs48 =             -0.01909830025156*phi_old[0];
const double crhs49 =             -0.01909830025156*phi_old[2];
const double crhs50 =             -0.01909830025156*phi_old[3];
const double crhs51 =             crhs42*(crhs43 + crhs44 + crhs45 + crhs46 + crhs47 + crhs48 + crhs49 + crhs50);
const double crhs52 =             0.1381966*phi[0];
const double crhs53 =             0.5854102*phi[1];
const double crhs54 =             0.1381966*phi[2];
const double crhs55 =             0.1381966*phi[3];
const double crhs56 =             crhs54 + crhs55;
const double crhs57 =             crhs52 + crhs53 + crhs56;
const double crhs58 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs59 =             crhs57*crhs58;
const double crhs60 =             0.1381966*crhs59;
const double crhs61 =             -crhs60;
const double crhs62 =             crhs60*tau[1];
const double crhs63 =             tau[1]*(crhs24 + crhs32 + crhs40);
const double crhs64 =             0.1381966*crhs63;
const double crhs65 =             0.19999999899376*f[2];
const double crhs66 =             0.19999999899376*f[3];
const double crhs67 =             phi_subscale_gauss[0]*tau[0];
const double crhs68 =             phi_subscale_gauss[2]*tau[2];
const double crhs69 =             -0.1381966*crhs68;
const double crhs70 =             phi_subscale_gauss[3]*tau[3];
const double crhs71 =             -0.1381966*crhs70;
const double crhs72 =             crhs3*phi_subscale_gauss[0];
const double crhs73 =             crhs3*phi_subscale_gauss[2];
const double crhs74 =             0.1381966*crhs73;
const double crhs75 =             crhs3*phi_subscale_gauss[3];
const double crhs76 =             0.1381966*crhs75;
const double crhs77 =             0.1381966*f[1];
const double crhs78 =             tau[0]*(crhs77 + crhs9 + 0.5854102*f[0]);
const double crhs79 =             0.1381966*prj[1];
const double crhs80 =             tau[0]*(crhs15 + crhs79 + 0.5854102*prj[0]);
const double crhs81 =             crhs6 + crhs77;
const double crhs82 =             tau[2]*(crhs8 + crhs81 + 0.5854102*f[2]);
const double crhs83 =             -0.1381966*crhs82;
const double crhs84 =             crhs12 + crhs79;
const double crhs85 =             tau[2]*(crhs14 + crhs84 + 0.5854102*prj[2]);
const double crhs86 =             -0.1381966*crhs85;
const double crhs87 =             tau[3]*(crhs7 + crhs81 + 0.5854102*f[3]);
const double crhs88 =             -0.1381966*crhs87;
const double crhs89 =             tau[3]*(crhs13 + crhs84 + 0.5854102*prj[3]);
const double crhs90 =             -0.1381966*crhs89;
const double crhs91 =             4*alpha*crhs23;
const double crhs92 =             4*alpha*crhs31;
const double crhs93 =             4*alpha*crhs39;
const double crhs94 =             DN_DX_0_0*alpha*crhs23;
const double crhs95 =             DN_DX_0_1*alpha*crhs31;
const double crhs96 =             DN_DX_0_2*alpha*crhs39;
const double crhs97 =             0.1381966*v(1,0);
const double crhs98 =             crhs18 + crhs97;
const double crhs99 =             crhs19 + crhs98 + 0.5854102*v(3,0);
const double crhs100 =             crhs23*crhs99;
const double crhs101 =             -0.1381966*crhs100;
const double crhs102 =             crhs20 + crhs98 + 0.5854102*v(2,0);
const double crhs103 =             crhs102*crhs23;
const double crhs104 =             -0.1381966*crhs103;
const double crhs105 =             crhs21 + crhs97 + 0.5854102*v(0,0);
const double crhs106 =             crhs105*crhs23;
const double crhs107 =             0.1381966*v(1,1);
const double crhs108 =             crhs107 + crhs26;
const double crhs109 =             crhs108 + crhs27 + 0.5854102*v(3,1);
const double crhs110 =             crhs109*crhs31;
const double crhs111 =             -0.1381966*crhs110;
const double crhs112 =             crhs108 + crhs28 + 0.5854102*v(2,1);
const double crhs113 =             crhs112*crhs31;
const double crhs114 =             -0.1381966*crhs113;
const double crhs115 =             crhs107 + crhs29 + 0.5854102*v(0,1);
const double crhs116 =             crhs115*crhs31;
const double crhs117 =             0.1381966*v(1,2);
const double crhs118 =             crhs117 + crhs34;
const double crhs119 =             crhs118 + crhs35 + 0.5854102*v(3,2);
const double crhs120 =             crhs119*crhs39;
const double crhs121 =             -0.1381966*crhs120;
const double crhs122 =             crhs118 + crhs36 + 0.5854102*v(2,2);
const double crhs123 =             crhs122*crhs39;
const double crhs124 =             -0.1381966*crhs123;
const double crhs125 =             crhs117 + crhs37 + 0.5854102*v(0,2);
const double crhs126 =             crhs125*crhs39;
const double crhs127 =             explicit_step_coefficient*tau[0];
const double crhs128 =             0.08090169924532*phi[2];
const double crhs129 =             0.08090169924532*phi[3];
const double crhs130 =             -0.08090169924532*phi_old[2];
const double crhs131 =             -0.08090169924532*phi_old[3];
const double crhs132 =             crhs128 + crhs129 + crhs130 + crhs131;
const double crhs133 =             explicit_step_coefficient*tau[2];
const double crhs134 =             0.01909830025156*phi[1];
const double crhs135 =             -0.01909830025156*phi_old[1];
const double crhs136 =             crhs133*(crhs128 + crhs130 + crhs134 + crhs135 + crhs45 + crhs47 + crhs48 + crhs50);
const double crhs137 =             explicit_step_coefficient*tau[3];
const double crhs138 =             crhs137*(crhs129 + crhs131 + crhs134 + crhs135 + crhs45 + crhs46 + crhs48 + crhs49);
const double crhs139 =             0.1381966*phi[1];
const double crhs140 =             crhs139 + crhs52;
const double crhs141 =             0.5854102*phi[3];
const double crhs142 =             crhs140 + crhs141 + crhs54;
const double crhs143 =             crhs142*crhs58;
const double crhs144 =             0.1381966*crhs143;
const double crhs145 =             -crhs144;
const double crhs146 =             0.5854102*phi[2];
const double crhs147 =             crhs140 + crhs146 + crhs55;
const double crhs148 =             crhs147*crhs58;
const double crhs149 =             0.1381966*crhs148;
const double crhs150 =             -crhs149;
const double crhs151 =             0.5854102*phi[0];
const double crhs152 =             crhs139 + crhs151 + crhs56;
const double crhs153 =             crhs152*crhs58;
const double crhs154 =             0.5854102*crhs153;
const double crhs155 =             crhs3*phi_subscale_gauss[0]*tau[0];
const double crhs156 =             DN_DX_0_0*crhs105 + DN_DX_0_1*crhs115 + DN_DX_0_2*crhs125;
const double crhs157 =             crhs3*phi_subscale_gauss[1]*tau[1];
const double crhs158 =             DN_DX_0_0*crhs22 + DN_DX_0_1*crhs30 + DN_DX_0_2*crhs38;
const double crhs159 =             crhs3*phi_subscale_gauss[2]*tau[2];
const double crhs160 =             DN_DX_0_0*crhs102 + DN_DX_0_1*crhs112 + DN_DX_0_2*crhs122;
const double crhs161 =             crhs3*phi_subscale_gauss[3]*tau[3];
const double crhs162 =             DN_DX_0_0*crhs99 + DN_DX_0_1*crhs109 + DN_DX_0_2*crhs119;
const double crhs163 =             crhs149*tau[2];
const double crhs164 =             crhs144*tau[3];
const double crhs165 =             -0.1381966*phi_old[2];
const double crhs166 =             -0.1381966*phi_old[3];
const double crhs167 =             crhs165 + crhs166;
const double crhs168 =             -0.1381966*phi_old[1];
const double crhs169 =             explicit_step_coefficient*tau[0]*(crhs139 + crhs151 + crhs167 + crhs168 + crhs54 + crhs55 - 0.5854102*phi_old[0]);
const double crhs170 =             -0.1381966*phi_old[0];
const double crhs171 =             explicit_step_coefficient*tau[1]*(crhs167 + crhs170 + crhs52 + crhs53 + crhs54 + crhs55 - 0.5854102*phi_old[1]);
const double crhs172 =             crhs168 + crhs170;
const double crhs173 =             explicit_step_coefficient*tau[2]*(crhs139 + crhs146 + crhs166 + crhs172 + crhs52 + crhs55 - 0.5854102*phi_old[2]);
const double crhs174 =             explicit_step_coefficient*tau[3]*(crhs139 + crhs141 + crhs165 + crhs172 + crhs52 + crhs54 - 0.5854102*phi_old[3]);
const double crhs175 =             tau[0]*(crhs106 + crhs116 + crhs126);
const double crhs176 =             tau[2]*(crhs103 + crhs113 + crhs123);
const double crhs177 =             0.1381966*crhs176;
const double crhs178 =             tau[3]*(crhs100 + crhs110 + crhs120);
const double crhs179 =             0.1381966*crhs178;
const double crhs180 =             crhs152*crhs58*tau[0];
const double crhs181 =             crhs57*crhs58*tau[1];
const double crhs182 =             crhs147*crhs58*tau[2];
const double crhs183 =             crhs142*crhs58*tau[3];
const double crhs184 =             0.19999999899376*f[0];
const double crhs185 =             -0.1381966*crhs67;
const double crhs186 =             0.1381966*crhs72;
const double crhs187 =             -0.1381966*crhs78;
const double crhs188 =             -0.1381966*crhs80;
const double crhs189 =             -0.1381966*crhs106;
const double crhs190 =             -0.1381966*crhs116;
const double crhs191 =             -0.1381966*crhs126;
const double crhs192 =             0.08090169924532*phi[0];
const double crhs193 =             -0.08090169924532*phi_old[0];
const double crhs194 =             crhs127*(crhs134 + crhs135 + crhs192 + crhs193 + crhs46 + crhs47 + crhs49 + crhs50);
const double crhs195 =             0.1381966*crhs153;
const double crhs196 =             -crhs195;
const double crhs197 =             crhs195*tau[0];
const double crhs198 =             0.1381966*crhs175;
const double crhs199 =             DN_DX_1_0*alpha*crhs23;
const double crhs200 =             DN_DX_1_1*alpha*crhs31;
const double crhs201 =             DN_DX_1_2*alpha*crhs39;
const double crhs202 =             0.5854102*crhs59;
const double crhs203 =             DN_DX_1_0*crhs105 + DN_DX_1_1*crhs115 + DN_DX_1_2*crhs125;
const double crhs204 =             DN_DX_1_0*crhs22 + DN_DX_1_1*crhs30 + DN_DX_1_2*crhs38;
const double crhs205 =             DN_DX_1_0*crhs102 + DN_DX_1_1*crhs112 + DN_DX_1_2*crhs122;
const double crhs206 =             DN_DX_1_0*crhs99 + DN_DX_1_1*crhs109 + DN_DX_1_2*crhs119;
const double crhs207 =             crhs0 + crhs11 + crhs17 + crhs184 + crhs185 + crhs186 + crhs187 + crhs188 + crhs189 + crhs190 + crhs191 + crhs194 + crhs196 + crhs197 + crhs198 + crhs2 + crhs25 + crhs33 + crhs41 + crhs5 + crhs51 + crhs61 + crhs62 + crhs64;
const double crhs208 =             DN_DX_2_0*alpha*crhs23;
const double crhs209 =             DN_DX_2_1*alpha*crhs31;
const double crhs210 =             DN_DX_2_2*alpha*crhs39;
const double crhs211 =             crhs192 + crhs193 + crhs43 + crhs44;
const double crhs212 =             0.5854102*crhs148;
const double crhs213 =             DN_DX_2_0*crhs105 + DN_DX_2_1*crhs115 + DN_DX_2_2*crhs125;
const double crhs214 =             DN_DX_2_0*crhs22 + DN_DX_2_1*crhs30 + DN_DX_2_2*crhs38;
const double crhs215 =             DN_DX_2_0*crhs102 + DN_DX_2_1*crhs112 + DN_DX_2_2*crhs122;
const double crhs216 =             DN_DX_2_0*crhs99 + DN_DX_2_1*crhs109 + DN_DX_2_2*crhs119;
const double crhs217 =             DN_DX_3_0*alpha*crhs23;
const double crhs218 =             DN_DX_3_1*alpha*crhs31;
const double crhs219 =             DN_DX_3_2*alpha*crhs39;
const double crhs220 =             0.5854102*crhs143;
const double crhs221 =             DN_DX_3_0*crhs105 + DN_DX_3_1*crhs115 + DN_DX_3_2*crhs125;
const double crhs222 =             DN_DX_3_0*crhs22 + DN_DX_3_1*crhs30 + DN_DX_3_2*crhs38;
const double crhs223 =             DN_DX_3_0*crhs102 + DN_DX_3_1*crhs112 + DN_DX_3_2*crhs122;
const double crhs224 =             DN_DX_3_0*crhs99 + DN_DX_3_1*crhs109 + DN_DX_3_2*crhs119;
            rhs[0]=-DN_DX_0_0*crhs91 - DN_DX_0_1*crhs92 - DN_DX_0_2*crhs93 + crhs0 + crhs10*crhs158 + crhs101 + crhs104 - 0.5854102*crhs106 + crhs11 + crhs111 + crhs114 - 0.5854102*crhs116 + crhs121 + crhs124 - 0.5854102*crhs126 + crhs127*(crhs132 + crhs43 + crhs44 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs136 + crhs138 + crhs145 + crhs150 + crhs154*tau[0] - crhs154 + crhs155*crhs156 - crhs156*crhs169 - crhs156*crhs175 - crhs156*crhs180 + crhs156*crhs78 + crhs156*crhs80 + crhs157*crhs158 + crhs158*crhs16 - crhs158*crhs171 - crhs158*crhs181 - crhs158*crhs63 + crhs159*crhs160 - crhs160*crhs173 - crhs160*crhs176 - crhs160*crhs182 + crhs160*crhs82 + crhs160*crhs85 + crhs161*crhs162 - crhs162*crhs174 - crhs162*crhs178 - crhs162*crhs183 + crhs162*crhs87 + crhs162*crhs89 + crhs163 + crhs164 + crhs17 + 0.5854102*crhs175 + crhs177 + crhs179 + crhs2 + crhs25 + crhs33 + crhs41 + crhs5 + crhs51 + crhs61 + crhs62 + crhs64 + crhs65 + crhs66 - 0.5854102*crhs67 + crhs69 + crhs71 + 0.5854102*crhs72 + crhs74 + crhs76 - 0.5854102*crhs78 - 0.5854102*crhs80 + crhs83 + crhs86 + crhs88 + crhs90 + crhs94*tau[0] + crhs94*tau[1] + crhs94*tau[2] + crhs94*tau[3] + crhs95*tau[0] + crhs95*tau[1] + crhs95*tau[2] + crhs95*tau[3] + crhs96*tau[0] + crhs96*tau[1] + crhs96*tau[2] + crhs96*tau[3] + 0.40000000301872*f[0];
            rhs[1]=-DN_DX_1_0*crhs91 - DN_DX_1_1*crhs92 - DN_DX_1_2*crhs93 - 0.5854102*crhs1 + crhs10*crhs204 - 0.5854102*crhs10 + crhs101 + crhs104 + crhs111 + crhs114 + crhs121 + crhs124 + crhs136 + crhs138 + crhs145 + crhs150 + crhs155*crhs203 + crhs157*crhs204 + crhs159*crhs205 + crhs16*crhs204 - 0.5854102*crhs16 + crhs161*crhs206 + crhs163 + crhs164 - crhs169*crhs203 - crhs171*crhs204 - crhs173*crhs205 - crhs174*crhs206 - crhs175*crhs203 - crhs176*crhs205 + crhs177 - crhs178*crhs206 + crhs179 - crhs180*crhs203 - crhs181*crhs204 - crhs182*crhs205 - crhs183*crhs206 + crhs184 + crhs185 + crhs186 + crhs187 + crhs188 + crhs189 + crhs190 + crhs191 + crhs194 + crhs196 + crhs197 + crhs198 + crhs199*tau[0] + crhs199*tau[1] + crhs199*tau[2] + crhs199*tau[3] + crhs200*tau[0] + crhs200*tau[1] + crhs200*tau[2] + crhs200*tau[3] + crhs201*tau[0] + crhs201*tau[1] + crhs201*tau[2] + crhs201*tau[3] + crhs202*tau[1] - crhs202 + crhs203*crhs78 + crhs203*crhs80 - crhs204*crhs63 + crhs205*crhs82 + crhs205*crhs85 + crhs206*crhs87 + crhs206*crhs89 - 0.5854102*crhs24 - 0.5854102*crhs32 + 0.5854102*crhs4 - 0.5854102*crhs40 + crhs42*(crhs132 + crhs192 + crhs193 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) + 0.5854102*crhs63 + crhs65 + crhs66 + crhs69 + crhs71 + crhs74 + crhs76 + crhs83 + crhs86 + crhs88 + crhs90 + 0.40000000301872*f[1];
            rhs[2]=-DN_DX_2_0*crhs91 - DN_DX_2_1*crhs92 - DN_DX_2_2*crhs93 + crhs10*crhs214 + crhs101 - 0.5854102*crhs103 + crhs111 - 0.5854102*crhs113 + crhs121 - 0.5854102*crhs123 + crhs133*(crhs129 + crhs131 + crhs211 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) + crhs138 + crhs145 + crhs155*crhs213 + crhs157*crhs214 + crhs159*crhs215 + crhs16*crhs214 + crhs161*crhs216 + crhs164 - crhs169*crhs213 - crhs171*crhs214 - crhs173*crhs215 - crhs174*crhs216 - crhs175*crhs213 - crhs176*crhs215 + 0.5854102*crhs176 - crhs178*crhs216 + crhs179 - crhs180*crhs213 - crhs181*crhs214 - crhs182*crhs215 - crhs183*crhs216 + crhs207 + crhs208*tau[0] + crhs208*tau[1] + crhs208*tau[2] + crhs208*tau[3] + crhs209*tau[0] + crhs209*tau[1] + crhs209*tau[2] + crhs209*tau[3] + crhs210*tau[0] + crhs210*tau[1] + crhs210*tau[2] + crhs210*tau[3] + crhs212*tau[2] - crhs212 + crhs213*crhs78 + crhs213*crhs80 - crhs214*crhs63 + crhs215*crhs82 + crhs215*crhs85 + crhs216*crhs87 + crhs216*crhs89 - 0.5854102*crhs68 + crhs71 + 0.5854102*crhs73 + crhs76 - 0.5854102*crhs82 - 0.5854102*crhs85 + crhs88 + crhs90 + 0.40000000301872*f[2] + 0.19999999899376*f[3];
            rhs[3]=-DN_DX_3_0*crhs91 - DN_DX_3_1*crhs92 - DN_DX_3_2*crhs93 + crhs10*crhs222 - 0.5854102*crhs100 + crhs104 - 0.5854102*crhs110 + crhs114 - 0.5854102*crhs120 + crhs124 + crhs136 + crhs137*(crhs128 + crhs130 + crhs211 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) + crhs150 + crhs155*crhs221 + crhs157*crhs222 + crhs159*crhs223 + crhs16*crhs222 + crhs161*crhs224 + crhs163 - crhs169*crhs221 - crhs171*crhs222 - crhs173*crhs223 - crhs174*crhs224 - crhs175*crhs221 - crhs176*crhs223 + crhs177 - crhs178*crhs224 + 0.5854102*crhs178 - crhs180*crhs221 - crhs181*crhs222 - crhs182*crhs223 - crhs183*crhs224 + crhs207 + crhs217*tau[0] + crhs217*tau[1] + crhs217*tau[2] + crhs217*tau[3] + crhs218*tau[0] + crhs218*tau[1] + crhs218*tau[2] + crhs218*tau[3] + crhs219*tau[0] + crhs219*tau[1] + crhs219*tau[2] + crhs219*tau[3] + crhs220*tau[3] - crhs220 + crhs221*crhs78 + crhs221*crhs80 - crhs222*crhs63 + crhs223*crhs82 + crhs223*crhs85 + crhs224*crhs87 + crhs224*crhs89 + crhs69 - 0.5854102*crhs70 + crhs74 + 0.5854102*crhs75 + crhs83 + crhs86 - 0.5854102*crhs87 - 0.5854102*crhs89 + 0.19999999899376*f[2] + 0.40000000301872*f[3];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::DCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             -0.25*f[1];
const double crhs1 =             1.0/delta_time;
const double crhs2 =             crhs1*phi_subscale_gauss[1];
const double crhs3 =             -0.166666666666667*crhs2;
const double crhs4 =             0.166666666666667*phi[0];
const double crhs5 =             0.166666666666667*phi[2];
const double crhs6 =             crhs4 + crhs5 + 0.666666666666667*phi[1];
const double crhs7 =             -0.166666666666667*phi_old[0];
const double crhs8 =             -0.166666666666667*phi_old[2];
const double crhs9 =             explicit_step_coefficient*(crhs6 + crhs7 + crhs8 - 0.666666666666667*phi_old[1]);
const double crhs10 =             0.166666666666667*crhs9;
const double crhs11 =             0.166666666666667*v(0,0);
const double crhs12 =             0.166666666666667*v(2,0);
const double crhs13 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs14 =             crhs13*(crhs11 + crhs12 + 0.666666666666667*v(1,0));
const double crhs15 =             0.166666666666667*crhs14;
const double crhs16 =             0.166666666666667*v(0,1);
const double crhs17 =             0.166666666666667*v(2,1);
const double crhs18 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs19 =             crhs18*(crhs16 + crhs17 + 0.666666666666667*v(1,1));
const double crhs20 =             0.166666666666667*crhs19;
const double crhs21 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs22 =             crhs21*crhs6;
const double crhs23 =             0.166666666666667*crhs22;
const double crhs24 =             -0.25*f[2];
const double crhs25 =             crhs1*phi_subscale_gauss[0];
const double crhs26 =             crhs1*phi_subscale_gauss[2];
const double crhs27 =             -0.166666666666667*crhs26;
const double crhs28 =             3*alpha*crhs13;
const double crhs29 =             3*alpha*crhs18;
const double crhs30 =             0.166666666666667*phi[1];
const double crhs31 =             crhs30 + crhs4 + 0.666666666666667*phi[2];
const double crhs32 =             -0.166666666666667*phi_old[1];
const double crhs33 =             explicit_step_coefficient*(crhs31 + crhs32 + crhs7 - 0.666666666666667*phi_old[2]);
const double crhs34 =             0.166666666666667*crhs33;
const double crhs35 =             crhs30 + crhs5 + 0.666666666666667*phi[0];
const double crhs36 =             explicit_step_coefficient*(crhs32 + crhs35 + crhs8 - 0.666666666666667*phi_old[0]);
const double crhs37 =             0.166666666666667*v(1,0);
const double crhs38 =             crhs13*(crhs11 + crhs37 + 0.666666666666667*v(2,0));
const double crhs39 =             0.166666666666667*crhs38;
const double crhs40 =             crhs13*(crhs12 + crhs37 + 0.666666666666667*v(0,0));
const double crhs41 =             0.166666666666667*v(1,1);
const double crhs42 =             crhs18*(crhs16 + crhs41 + 0.666666666666667*v(2,1));
const double crhs43 =             0.166666666666667*crhs42;
const double crhs44 =             crhs18*(crhs17 + crhs41 + 0.666666666666667*v(0,1));
const double crhs45 =             crhs21*crhs31;
const double crhs46 =             0.166666666666667*crhs45;
const double crhs47 =             crhs21*crhs35;
const double crhs48 =             -0.166666666666667*crhs25 + 0.166666666666667*crhs36 + 0.166666666666667*crhs40 + 0.166666666666667*crhs44 + 0.166666666666667*crhs47 - 0.25*f[0];
            rhs[0]=DN_DX_0_0*crhs28 + DN_DX_0_1*crhs29 + crhs0 + crhs10 + crhs15 + crhs20 + crhs23 + crhs24 - 0.666666666666667*crhs25 + crhs27 + crhs3 + crhs34 + 0.666666666666667*crhs36 + crhs39 + 0.666666666666667*crhs40 + crhs43 + 0.666666666666667*crhs44 + crhs46 + 0.666666666666667*crhs47 - 0.5*f[0];
            rhs[1]=DN_DX_1_0*crhs28 + DN_DX_1_1*crhs29 + 0.666666666666667*crhs14 + 0.666666666666667*crhs19 - 0.666666666666667*crhs2 + 0.666666666666667*crhs22 + crhs24 + crhs27 + crhs34 + crhs39 + crhs43 + crhs46 + crhs48 + 0.666666666666667*crhs9 - 0.5*f[1];
            rhs[2]=DN_DX_2_0*crhs28 + DN_DX_2_1*crhs29 + crhs0 + crhs10 + crhs15 + crhs20 + crhs23 - 0.666666666666667*crhs26 + crhs3 + 0.666666666666667*crhs33 + 0.666666666666667*crhs38 + 0.666666666666667*crhs42 + 0.666666666666667*crhs45 + crhs48 - 0.5*f[2];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<3,4>::DCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             -0.19999999899376*f[1];
const double crhs1 =             1.0/delta_time;
const double crhs2 =             crhs1*phi_subscale_gauss[1];
const double crhs3 =             -0.1381966*crhs2;
const double crhs4 =             -0.1381966*phi_old[2];
const double crhs5 =             -0.1381966*phi_old[3];
const double crhs6 =             crhs4 + crhs5;
const double crhs7 =             0.1381966*phi[0];
const double crhs8 =             0.5854102*phi[1];
const double crhs9 =             0.1381966*phi[2];
const double crhs10 =             0.1381966*phi[3];
const double crhs11 =             -0.1381966*phi_old[0];
const double crhs12 =             explicit_step_coefficient*(crhs10 + crhs11 + crhs6 + crhs7 + crhs8 + crhs9 - 0.5854102*phi_old[1]);
const double crhs13 =             0.1381966*crhs12;
const double crhs14 =             0.1381966*v(0,0);
const double crhs15 =             0.1381966*v(2,0);
const double crhs16 =             0.1381966*v(3,0);
const double crhs17 =             crhs15 + crhs16;
const double crhs18 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs19 =             crhs18*(crhs14 + crhs17 + 0.5854102*v(1,0));
const double crhs20 =             0.1381966*crhs19;
const double crhs21 =             0.1381966*v(0,1);
const double crhs22 =             0.1381966*v(2,1);
const double crhs23 =             0.1381966*v(3,1);
const double crhs24 =             crhs22 + crhs23;
const double crhs25 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs26 =             crhs25*(crhs21 + crhs24 + 0.5854102*v(1,1));
const double crhs27 =             0.1381966*crhs26;
const double crhs28 =             0.1381966*v(0,2);
const double crhs29 =             0.1381966*v(2,2);
const double crhs30 =             0.1381966*v(3,2);
const double crhs31 =             crhs29 + crhs30;
const double crhs32 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs33 =             crhs32*(crhs28 + crhs31 + 0.5854102*v(1,2));
const double crhs34 =             0.1381966*crhs33;
const double crhs35 =             crhs10 + crhs9;
const double crhs36 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs37 =             crhs36*(crhs35 + crhs7 + crhs8);
const double crhs38 =             0.1381966*crhs37;
const double crhs39 =             -0.19999999899376*f[2];
const double crhs40 =             -0.19999999899376*f[3];
const double crhs41 =             crhs1*phi_subscale_gauss[0];
const double crhs42 =             crhs1*phi_subscale_gauss[2];
const double crhs43 =             -0.1381966*crhs42;
const double crhs44 =             crhs1*phi_subscale_gauss[3];
const double crhs45 =             -0.1381966*crhs44;
const double crhs46 =             4*alpha*crhs18;
const double crhs47 =             4*alpha*crhs25;
const double crhs48 =             4*alpha*crhs32;
const double crhs49 =             -0.1381966*phi_old[1];
const double crhs50 =             crhs11 + crhs49;
const double crhs51 =             0.1381966*phi[1];
const double crhs52 =             0.5854102*phi[3];
const double crhs53 =             explicit_step_coefficient*(crhs4 + crhs50 + crhs51 + crhs52 + crhs7 + crhs9 - 0.5854102*phi_old[3]);
const double crhs54 =             0.1381966*crhs53;
const double crhs55 =             0.5854102*phi[2];
const double crhs56 =             explicit_step_coefficient*(crhs10 + crhs5 + crhs50 + crhs51 + crhs55 + crhs7 - 0.5854102*phi_old[2]);
const double crhs57 =             0.1381966*crhs56;
const double crhs58 =             0.5854102*phi[0];
const double crhs59 =             explicit_step_coefficient*(crhs10 + crhs49 + crhs51 + crhs58 + crhs6 + crhs9 - 0.5854102*phi_old[0]);
const double crhs60 =             0.1381966*v(1,0);
const double crhs61 =             crhs14 + crhs60;
const double crhs62 =             crhs18*(crhs15 + crhs61 + 0.5854102*v(3,0));
const double crhs63 =             0.1381966*crhs62;
const double crhs64 =             crhs18*(crhs16 + crhs61 + 0.5854102*v(2,0));
const double crhs65 =             0.1381966*crhs64;
const double crhs66 =             crhs18*(crhs17 + crhs60 + 0.5854102*v(0,0));
const double crhs67 =             0.1381966*v(1,1);
const double crhs68 =             crhs21 + crhs67;
const double crhs69 =             crhs25*(crhs22 + crhs68 + 0.5854102*v(3,1));
const double crhs70 =             0.1381966*crhs69;
const double crhs71 =             crhs25*(crhs23 + crhs68 + 0.5854102*v(2,1));
const double crhs72 =             0.1381966*crhs71;
const double crhs73 =             crhs25*(crhs24 + crhs67 + 0.5854102*v(0,1));
const double crhs74 =             0.1381966*v(1,2);
const double crhs75 =             crhs28 + crhs74;
const double crhs76 =             crhs32*(crhs29 + crhs75 + 0.5854102*v(3,2));
const double crhs77 =             0.1381966*crhs76;
const double crhs78 =             crhs32*(crhs30 + crhs75 + 0.5854102*v(2,2));
const double crhs79 =             0.1381966*crhs78;
const double crhs80 =             crhs32*(crhs31 + crhs74 + 0.5854102*v(0,2));
const double crhs81 =             crhs51 + crhs7;
const double crhs82 =             crhs36*(crhs52 + crhs81 + crhs9);
const double crhs83 =             0.1381966*crhs82;
const double crhs84 =             crhs36*(crhs10 + crhs55 + crhs81);
const double crhs85 =             0.1381966*crhs84;
const double crhs86 =             crhs36*(crhs35 + crhs51 + crhs58);
const double crhs87 =             -0.19999999899376*f[0];
const double crhs88 =             -0.1381966*crhs41;
const double crhs89 =             0.1381966*crhs59;
const double crhs90 =             0.1381966*crhs66;
const double crhs91 =             0.1381966*crhs73;
const double crhs92 =             0.1381966*crhs80;
const double crhs93 =             0.1381966*crhs86;
const double crhs94 =             crhs0 + crhs13 + crhs20 + crhs27 + crhs3 + crhs34 + crhs38 + crhs87 + crhs88 + crhs89 + crhs90 + crhs91 + crhs92 + crhs93;
            rhs[0]=DN_DX_0_0*crhs46 + DN_DX_0_1*crhs47 + DN_DX_0_2*crhs48 + crhs0 + crhs13 + crhs20 + crhs27 + crhs3 + crhs34 + crhs38 + crhs39 + crhs40 - 0.5854102*crhs41 + crhs43 + crhs45 + crhs54 + crhs57 + 0.5854102*crhs59 + crhs63 + crhs65 + 0.5854102*crhs66 + crhs70 + crhs72 + 0.5854102*crhs73 + crhs77 + crhs79 + 0.5854102*crhs80 + crhs83 + crhs85 + 0.5854102*crhs86 - 0.40000000301872*f[0];
            rhs[1]=DN_DX_1_0*crhs46 + DN_DX_1_1*crhs47 + DN_DX_1_2*crhs48 + 0.5854102*crhs12 + 0.5854102*crhs19 - 0.5854102*crhs2 + 0.5854102*crhs26 + 0.5854102*crhs33 + 0.5854102*crhs37 + crhs39 + crhs40 + crhs43 + crhs45 + crhs54 + crhs57 + crhs63 + crhs65 + crhs70 + crhs72 + crhs77 + crhs79 + crhs83 + crhs85 + crhs87 + crhs88 + crhs89 + crhs90 + crhs91 + crhs92 + crhs93 - 0.40000000301872*f[1];
            rhs[2]=DN_DX_2_0*crhs46 + DN_DX_2_1*crhs47 + DN_DX_2_2*crhs48 - 0.5854102*crhs42 + crhs45 + crhs54 + 0.5854102*crhs56 + crhs63 + 0.5854102*crhs64 + crhs70 + 0.5854102*crhs71 + crhs77 + 0.5854102*crhs78 + crhs83 + 0.5854102*crhs84 + crhs94 - 0.40000000301872*f[2] - 0.19999999899376*f[3];
            rhs[3]=DN_DX_3_0*crhs46 + DN_DX_3_1*crhs47 + DN_DX_3_2*crhs48 + crhs43 - 0.5854102*crhs44 + 0.5854102*crhs53 + crhs57 + 0.5854102*crhs62 + crhs65 + 0.5854102*crhs69 + crhs72 + 0.5854102*crhs76 + crhs79 + 0.5854102*crhs82 + crhs85 + crhs94 - 0.19999999899376*f[2] - 0.40000000301872*f[3];


    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementData& rData,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element data
    const auto& N = rData.N;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau[g];
    const auto& phi_subscale_gauss = rData.unknown_subscale;
    const auto& prj = rData.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);

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
void DConvectionDiffusionExplicit<3,4>::UpdateUnknownSubgridScaleGaussPoint(
    ElementData& rData,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element data
    const auto& N = rData.N;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau[g];
    const auto& phi_subscale_gauss = rData.unknown_subscale;
    const auto& prj = rData.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);

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
void DConvectionDiffusionExplicit<TDim,TNumNodes>::DCalculateTau(
    ElementData& rData)
{
    KRATOS_TRY;

    // Calculate h
    double h = this->ComputeH(rData.DN_DX);
    // Calculate tau for each gauss point
    for(unsigned int g = 0; g<TNumNodes; g++) {
	const auto& N = row(rData.N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        const array_1d<double,3> vel_gauss = prod(N, rData.convective_velocity);
        double div_vel = 0;
        for(unsigned int node_element = 0; node_element<TNumNodes; node_element++) {
            for(unsigned int dim = 0; dim < TDim; dim++) {
                div_vel += rData.DN_DX(node_element,dim)*rData.convective_velocity(node_element,dim);
            }
        }
        const double norm_velocity = norm_2(vel_gauss);
        // Estimate tau
        // Dynamic part
        double inv_tau = 1.0/rData.delta_time;
        // Convection
        inv_tau += 2.0 * norm_velocity / h;
        inv_tau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
        // Diffusion
        inv_tau += 4.0 * rData.diffusivity / (h*h);
        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);
        rData.tau[g] = (1.0) / inv_tau;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class DConvectionDiffusionExplicit<2,3>;
template class DConvectionDiffusionExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}
