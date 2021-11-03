//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//

#include "two_fluid_navier_stokes_CN.h"
#include "custom_utilities/two_fluid_navier_stokes_CN_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::~TwoFluidNavierStokesCN() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesCN<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesCN>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesCN<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesCN>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (TElementData::ElementManagesTimeIntegration){
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        if (data.IsCut()){
            GeometryType::Pointer p_geom = this->pGetGeometry();
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;
            Matrix int_shape_function_neg;                                       // interface shape functions
            Matrix int_shape_function_enr_neg, int_shape_function_enr_pos;       // interface enriched shape functions
            GeometryType::ShapeFunctionsGradientsType int_shape_derivatives_neg; // interface shape functions derivatives
            Vector int_gauss_pts_weights;                                // interface Gauss points weights
            std::vector<array_1d<double,3>> int_normals_neg;                                 // interface normal vector based on the negative side
            Vector gauss_pts_curvature;                                  // curvatures calculated on interface Gauss points

            ModifiedShapeFunctions::Pointer p_modified_sh_func = pGetModifiedShapeFunctionsUtility(p_geom, data.Distance);

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg,
                p_modified_sh_func);

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side
                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                const unsigned int number_of_gauss_points = gauss_weights.size();
                array_1d<double, NumNodes> Ncenter;
                for (unsigned int i = 0; i < NumNodes; ++i){
                    Ncenter[i] = 1.0/NumNodes;
                }
                for (unsigned int g = 0; g < number_of_gauss_points; ++g){
                    UpdateIntegrationPointData(
                        data,
                        g,
                        gauss_weights[g],
                        row(shape_functions, g),
                        shape_derivatives[g]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                }
            } else {
                MatrixType Vtot = ZeroMatrix(NumNodes * (Dim + 1), NumNodes);
                MatrixType Htot = ZeroMatrix(NumNodes, NumNodes * (Dim + 1));
                MatrixType Kee_tot = ZeroMatrix(NumNodes, NumNodes);
                VectorType rhs_ee_tot = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); ++g_pos){
                    UpdateIntegrationPointData(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);

                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); ++g_neg){
                    UpdateIntegrationPointData(
                        data,
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                if (rCurrentProcessInfo[SURFACE_TENSION]){

                    AddSurfaceTensionContribution(
                        data,
                        p_modified_sh_func,
                        rLeftHandSideMatrix,
                        rRightHandSideVector,
                        Htot,
                        Vtot,
                        Kee_tot,
                        rhs_ee_tot
                    );

                } else{
                    // Without pressure gradient stabilization, volume ratio is checked during condensation
                    // Also, without surface tension, zero pressure difference is penalized
                    CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
                }

            }
        } else {
            //Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            ShapeFunctionDerivativesArrayType shape_derivatives;
            this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
            const unsigned int number_of_gauss_points = gauss_weights.size();
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; ++g){
                UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
    } else{
        KRATOS_ERROR << "TwoFluidNavierStokesCN is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokesCN<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
std::string TwoFluidNavierStokesCN<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokesCN" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr){
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);

    // Calculate current material response
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX,
    const typename TElementData::MatrixRowType& rNenr,
    const typename TElementData::ShapeDerivativesType& rDN_DXenr) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex,Weight,rN,rDN_DX,rNenr,rDN_DXenr);
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesCNData<2, 3>& rData) const
{
    const double theta = rData.theta;
    const BoundedMatrix<double,3,2> mid_step_velocity = theta*rData.Velocity + (1-theta)*rData.Velocity_OldStep1;
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += rDNDX(i,0)*mid_step_velocity(i,0);
        r_strain_rate[1] += rDNDX(i,1)*mid_step_velocity(i,1);
        r_strain_rate[2] += rDNDX(i,0)*mid_step_velocity(i,1) + rDNDX(i,1)*mid_step_velocity(i,0);
    }
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesCNData<3, 4>& rData) const
{
    const double theta = rData.theta;
    const BoundedMatrix<double,4,3> mid_step_velocity = theta*rData.Velocity + (1-theta)*rData.Velocity_OldStep1;;
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(6);
    for (unsigned int i = 0; i < 4; i++) {
        r_strain_rate[0] += rDNDX(i,0)*mid_step_velocity(i,0);
        r_strain_rate[1] += rDNDX(i,1)*mid_step_velocity(i,1);
        r_strain_rate[2] += rDNDX(i,2)*mid_step_velocity(i,2);
        r_strain_rate[3] += rDNDX(i,0)*mid_step_velocity(i,1) + rDNDX(i,1)*mid_step_velocity(i,0);
        r_strain_rate[4] += rDNDX(i,1)*mid_step_velocity(i,2) + rDNDX(i,2)*mid_step_velocity(i,1);
        r_strain_rate[5] += rDNDX(i,0)*mid_step_velocity(i,2) + rDNDX(i,2)*mid_step_velocity(i,0);
    }
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             DN(0,0)*theta;
const double clhs2 =             C(0,2)*DN(0,0);
const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
const double clhs4 =             DN(0,1)*theta;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             mu*stab_c1;
const double clhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs9 =             rho*sqrt(pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs6/stab_c2 + clhs9*h;
const double clhs11 =             clhs10*theta;
const double clhs12 =             1.0/dt;
const double clhs13 =             clhs12*rho;
const double clhs14 =             1.0*clhs13;
const double clhs15 =             DN(0,0)*clhs7 + DN(0,1)*clhs8;
const double clhs16 =             clhs15*theta;
const double clhs17 =             N[0]*rho;
const double clhs18 =             1.0*clhs12;
const double clhs19 =             N[0]*clhs18 + clhs16;
const double clhs20 =             1.0/(clhs13*dyn_tau + clhs6/pow(h, 2) + clhs9*stab_c2/h);
const double clhs21 =             clhs20*pow(rho, 2);
const double clhs22 =             clhs15*clhs21;
const double clhs23 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs24 =             clhs21*clhs23;
const double clhs25 =             N[0]*clhs24;
const double clhs26 =             pow(N[0], 2)*clhs14 + clhs16*clhs17 + clhs19*clhs22 + clhs19*clhs25;
const double clhs27 =             C(0,1)*DN(0,1) + clhs2;
const double clhs28 =             C(1,2)*DN(0,1);
const double clhs29 =             C(2,2)*DN(0,0) + clhs28;
const double clhs30 =             DN(0,0)*clhs10;
const double clhs31 =             DN(0,1)*clhs30;
const double clhs32 =             clhs20*rho;
const double clhs33 =             clhs15*clhs32;
const double clhs34 =             -N[0] + clhs17*clhs20*clhs23 + clhs33;
const double clhs35 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs36 =             C(0,2)*DN(1,0);
const double clhs37 =             C(2,2)*DN(1,1) + clhs36;
const double clhs38 =             DN(0,0)*DN(1,0);
const double clhs39 =             N[0]*clhs14;
const double clhs40 =             N[1]*clhs39;
const double clhs41 =             clhs11*clhs38 + clhs40;
const double clhs42 =             DN(1,0)*clhs7 + DN(1,1)*clhs8;
const double clhs43 =             clhs42*theta;
const double clhs44 =             N[1]*clhs18 + clhs43;
const double clhs45 =             clhs17*clhs43 + clhs22*clhs44 + clhs25*clhs44;
const double clhs46 =             C(0,1)*DN(1,1) + clhs36;
const double clhs47 =             C(1,2)*DN(1,1);
const double clhs48 =             C(2,2)*DN(1,0) + clhs47;
const double clhs49 =             DN(1,1)*clhs30;
const double clhs50 =             DN(0,0)*N[1];
const double clhs51 =             DN(1,0)*N[0];
const double clhs52 =             clhs23*clhs32;
const double clhs53 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs54 =             C(0,2)*DN(2,0);
const double clhs55 =             C(2,2)*DN(2,1) + clhs54;
const double clhs56 =             DN(0,0)*DN(2,0);
const double clhs57 =             N[2]*clhs39;
const double clhs58 =             clhs11*clhs56 + clhs57;
const double clhs59 =             DN(2,0)*clhs7 + DN(2,1)*clhs8;
const double clhs60 =             clhs59*theta;
const double clhs61 =             N[2]*clhs18 + clhs60;
const double clhs62 =             clhs17*clhs60 + clhs22*clhs61 + clhs25*clhs61;
const double clhs63 =             C(0,1)*DN(2,1) + clhs54;
const double clhs64 =             C(1,2)*DN(2,1);
const double clhs65 =             C(2,2)*DN(2,0) + clhs64;
const double clhs66 =             DN(2,1)*clhs30;
const double clhs67 =             DN(0,0)*N[2];
const double clhs68 =             DN(2,0)*N[0];
const double clhs69 =             C(0,1)*DN(0,0) + clhs28;
const double clhs70 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs71 =             pow(DN(0,1), 2);
const double clhs72 =             C(0,1)*DN(1,0) + clhs47;
const double clhs73 =             DN(0,1)*clhs10;
const double clhs74 =             DN(1,0)*clhs73;
const double clhs75 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs76 =             DN(0,1)*DN(1,1);
const double clhs77 =             clhs11*clhs76 + clhs40;
const double clhs78 =             DN(0,1)*N[1];
const double clhs79 =             DN(1,1)*N[0];
const double clhs80 =             C(0,1)*DN(2,0) + clhs64;
const double clhs81 =             DN(2,0)*clhs73;
const double clhs82 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs83 =             DN(0,1)*DN(2,1);
const double clhs84 =             clhs11*clhs83 + clhs57;
const double clhs85 =             DN(0,1)*N[2];
const double clhs86 =             DN(2,1)*N[0];
const double clhs87 =             N[0]*theta;
const double clhs88 =             clhs19*clhs32;
const double clhs89 =             clhs87 + clhs88;
const double clhs90 =             clhs32*clhs44;
const double clhs91 =             clhs20*(clhs38 + clhs76);
const double clhs92 =             clhs32*clhs61;
const double clhs93 =             clhs20*(clhs56 + clhs83);
const double clhs94 =             DN(1,0)*theta;
const double clhs95 =             DN(1,1)*theta;
const double clhs96 =             N[1]*rho;
const double clhs97 =             clhs21*clhs42;
const double clhs98 =             N[1]*clhs24;
const double clhs99 =             clhs16*clhs96 + clhs19*clhs97 + clhs19*clhs98;
const double clhs100 =             clhs32*clhs42;
const double clhs101 =             pow(DN(1,0), 2);
const double clhs102 =             pow(N[1], 2)*clhs14 + clhs43*clhs96 + clhs44*clhs97 + clhs44*clhs98;
const double clhs103 =             DN(1,0)*clhs10;
const double clhs104 =             DN(1,1)*clhs103;
const double clhs105 =             N[1]*clhs52 - N[1] + clhs100;
const double clhs106 =             DN(1,0)*DN(2,0);
const double clhs107 =             N[1]*N[2]*clhs14;
const double clhs108 =             clhs106*clhs11 + clhs107;
const double clhs109 =             clhs60*clhs96 + clhs61*clhs97 + clhs61*clhs98;
const double clhs110 =             DN(2,1)*clhs103;
const double clhs111 =             DN(1,0)*N[2];
const double clhs112 =             DN(2,0)*N[1];
const double clhs113 =             pow(DN(1,1), 2);
const double clhs114 =             DN(2,0)*clhs10;
const double clhs115 =             DN(1,1)*clhs114;
const double clhs116 =             DN(1,1)*DN(2,1);
const double clhs117 =             clhs107 + clhs11*clhs116;
const double clhs118 =             DN(1,1)*N[2];
const double clhs119 =             DN(2,1)*N[1];
const double clhs120 =             N[1]*theta;
const double clhs121 =             clhs120 + clhs90;
const double clhs122 =             clhs20*(clhs106 + clhs116);
const double clhs123 =             DN(2,0)*theta;
const double clhs124 =             DN(2,1)*theta;
const double clhs125 =             N[2]*rho;
const double clhs126 =             clhs21*clhs59;
const double clhs127 =             N[2]*clhs24;
const double clhs128 =             clhs125*clhs16 + clhs126*clhs19 + clhs127*clhs19;
const double clhs129 =             clhs32*clhs59;
const double clhs130 =             clhs125*clhs43 + clhs126*clhs44 + clhs127*clhs44;
const double clhs131 =             pow(DN(2,0), 2);
const double clhs132 =             pow(N[2], 2)*clhs14 + clhs125*clhs60 + clhs126*clhs61 + clhs127*clhs61;
const double clhs133 =             DN(2,1)*clhs114;
const double clhs134 =             N[2]*clhs52 - N[2] + clhs129;
const double clhs135 =             pow(DN(2,1), 2);
const double clhs136 =             N[2]*theta + clhs92;
            lhs(0,0)=clhs0*clhs1 + clhs11*clhs5 + clhs26 + clhs3*clhs4;
            lhs(0,1)=theta*(DN(0,0)*clhs27 + DN(0,1)*clhs29 + clhs31);
            lhs(0,2)=DN(0,0)*clhs34;
            lhs(0,3)=clhs1*clhs35 + clhs37*clhs4 + clhs41 + clhs45;
            lhs(0,4)=theta*(DN(0,0)*clhs46 + DN(0,1)*clhs48 + clhs49);
            lhs(0,5)=DN(1,0)*clhs33 - clhs50 + clhs51*clhs52;
            lhs(0,6)=clhs1*clhs53 + clhs4*clhs55 + clhs58 + clhs62;
            lhs(0,7)=theta*(DN(0,0)*clhs63 + DN(0,1)*clhs65 + clhs66);
            lhs(0,8)=DN(2,0)*clhs33 + clhs52*clhs68 - clhs67;
            lhs(1,0)=theta*(DN(0,0)*clhs3 + DN(0,1)*clhs69 + clhs31);
            lhs(1,1)=clhs1*clhs29 + clhs11*clhs71 + clhs26 + clhs4*clhs70;
            lhs(1,2)=DN(0,1)*clhs34;
            lhs(1,3)=theta*(DN(0,0)*clhs37 + DN(0,1)*clhs72 + clhs74);
            lhs(1,4)=clhs1*clhs48 + clhs4*clhs75 + clhs45 + clhs77;
            lhs(1,5)=DN(1,1)*clhs33 + clhs52*clhs79 - clhs78;
            lhs(1,6)=theta*(DN(0,0)*clhs55 + DN(0,1)*clhs80 + clhs81);
            lhs(1,7)=clhs1*clhs65 + clhs4*clhs82 + clhs62 + clhs84;
            lhs(1,8)=DN(2,1)*clhs33 + clhs52*clhs86 - clhs85;
            lhs(2,0)=DN(0,0)*clhs89;
            lhs(2,1)=DN(0,1)*clhs89;
            lhs(2,2)=clhs20*(clhs5 + clhs71);
            lhs(2,3)=DN(0,0)*clhs90 + DN(1,0)*clhs87;
            lhs(2,4)=DN(0,1)*clhs90 + DN(1,1)*clhs87;
            lhs(2,5)=clhs91;
            lhs(2,6)=DN(0,0)*clhs92 + DN(2,0)*clhs87;
            lhs(2,7)=DN(0,1)*clhs92 + DN(2,1)*clhs87;
            lhs(2,8)=clhs93;
            lhs(3,0)=clhs0*clhs94 + clhs3*clhs95 + clhs41 + clhs99;
            lhs(3,1)=theta*(DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs74);
            lhs(3,2)=DN(0,0)*clhs100 + clhs50*clhs52 - clhs51;
            lhs(3,3)=clhs101*clhs11 + clhs102 + clhs35*clhs94 + clhs37*clhs95;
            lhs(3,4)=theta*(DN(1,0)*clhs46 + DN(1,1)*clhs48 + clhs104);
            lhs(3,5)=DN(1,0)*clhs105;
            lhs(3,6)=clhs108 + clhs109 + clhs53*clhs94 + clhs55*clhs95;
            lhs(3,7)=theta*(DN(1,0)*clhs63 + DN(1,1)*clhs65 + clhs110);
            lhs(3,8)=DN(2,0)*clhs100 - clhs111 + clhs112*clhs52;
            lhs(4,0)=theta*(DN(1,0)*clhs3 + DN(1,1)*clhs69 + clhs49);
            lhs(4,1)=clhs29*clhs94 + clhs70*clhs95 + clhs77 + clhs99;
            lhs(4,2)=DN(0,1)*clhs100 + clhs52*clhs78 - clhs79;
            lhs(4,3)=theta*(DN(1,0)*clhs37 + DN(1,1)*clhs72 + clhs104);
            lhs(4,4)=clhs102 + clhs11*clhs113 + clhs48*clhs94 + clhs75*clhs95;
            lhs(4,5)=DN(1,1)*clhs105;
            lhs(4,6)=theta*(DN(1,0)*clhs55 + DN(1,1)*clhs80 + clhs115);
            lhs(4,7)=clhs109 + clhs117 + clhs65*clhs94 + clhs82*clhs95;
            lhs(4,8)=DN(2,1)*clhs100 - clhs118 + clhs119*clhs52;
            lhs(5,0)=DN(1,0)*clhs88 + clhs50*theta;
            lhs(5,1)=DN(1,1)*clhs88 + clhs78*theta;
            lhs(5,2)=clhs91;
            lhs(5,3)=DN(1,0)*clhs121;
            lhs(5,4)=DN(1,1)*clhs121;
            lhs(5,5)=clhs20*(clhs101 + clhs113);
            lhs(5,6)=DN(1,0)*clhs92 + DN(2,0)*clhs120;
            lhs(5,7)=DN(1,1)*clhs92 + DN(2,1)*clhs120;
            lhs(5,8)=clhs122;
            lhs(6,0)=clhs0*clhs123 + clhs124*clhs3 + clhs128 + clhs58;
            lhs(6,1)=theta*(DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs81);
            lhs(6,2)=DN(0,0)*clhs129 + clhs52*clhs67 - clhs68;
            lhs(6,3)=clhs108 + clhs123*clhs35 + clhs124*clhs37 + clhs130;
            lhs(6,4)=theta*(DN(2,0)*clhs46 + DN(2,1)*clhs48 + clhs115);
            lhs(6,5)=DN(1,0)*clhs129 + clhs111*clhs52 - clhs112;
            lhs(6,6)=clhs11*clhs131 + clhs123*clhs53 + clhs124*clhs55 + clhs132;
            lhs(6,7)=theta*(DN(2,0)*clhs63 + DN(2,1)*clhs65 + clhs133);
            lhs(6,8)=DN(2,0)*clhs134;
            lhs(7,0)=theta*(DN(2,0)*clhs3 + DN(2,1)*clhs69 + clhs66);
            lhs(7,1)=clhs123*clhs29 + clhs124*clhs70 + clhs128 + clhs84;
            lhs(7,2)=DN(0,1)*clhs129 + clhs52*clhs85 - clhs86;
            lhs(7,3)=theta*(DN(2,0)*clhs37 + DN(2,1)*clhs72 + clhs110);
            lhs(7,4)=clhs117 + clhs123*clhs48 + clhs124*clhs75 + clhs130;
            lhs(7,5)=DN(1,1)*clhs129 + clhs118*clhs52 - clhs119;
            lhs(7,6)=theta*(DN(2,0)*clhs55 + DN(2,1)*clhs80 + clhs133);
            lhs(7,7)=clhs11*clhs135 + clhs123*clhs65 + clhs124*clhs82 + clhs132;
            lhs(7,8)=DN(2,1)*clhs134;
            lhs(8,0)=DN(2,0)*clhs88 + clhs67*theta;
            lhs(8,1)=DN(2,1)*clhs88 + clhs85*theta;
            lhs(8,2)=clhs93;
            lhs(8,3)=DN(2,0)*clhs90 + clhs111*theta;
            lhs(8,4)=DN(2,1)*clhs90 + clhs118*theta;
            lhs(8,5)=clhs122;
            lhs(8,6)=DN(2,0)*clhs136;
            lhs(8,7)=DN(2,1)*clhs136;
            lhs(8,8)=clhs20*(clhs131 + clhs135);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             DN(0,0)*theta;
const double clhs2 =             C(0,3)*DN(0,0);
const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
const double clhs4 =             DN(0,1)*theta;
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             DN(0,2)*theta;
const double clhs8 =             pow(DN(0,0), 2);
const double clhs9 =             mu*stab_c1;
const double clhs10 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs11 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs12 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs13 =             rho*sqrt(pow(clhs10, 2) + pow(clhs11, 2) + pow(clhs12, 2));
const double clhs14 =             clhs13*h + clhs9/stab_c2;
const double clhs15 =             clhs14*theta;
const double clhs16 =             1.0/dt;
const double clhs17 =             clhs16*rho;
const double clhs18 =             1.0*clhs17;
const double clhs19 =             DN(0,0)*clhs10 + DN(0,1)*clhs11 + DN(0,2)*clhs12;
const double clhs20 =             clhs19*theta;
const double clhs21 =             N[0]*rho;
const double clhs22 =             1.0*clhs16;
const double clhs23 =             N[0]*clhs22 + clhs20;
const double clhs24 =             1.0/(clhs13*stab_c2/h + clhs17*dyn_tau + clhs9/pow(h, 2));
const double clhs25 =             clhs24*pow(rho, 2);
const double clhs26 =             clhs19*clhs25;
const double clhs27 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs28 =             clhs25*clhs27;
const double clhs29 =             N[0]*clhs28;
const double clhs30 =             pow(N[0], 2)*clhs18 + clhs20*clhs21 + clhs23*clhs26 + clhs23*clhs29;
const double clhs31 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
const double clhs32 =             C(1,3)*DN(0,1);
const double clhs33 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs32;
const double clhs34 =             C(3,5)*DN(0,0);
const double clhs35 =             C(4,5)*DN(0,2);
const double clhs36 =             C(1,5)*DN(0,1) + clhs34 + clhs35;
const double clhs37 =             DN(0,0)*clhs14;
const double clhs38 =             DN(0,1)*clhs37;
const double clhs39 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs40 =             C(3,4)*DN(0,1);
const double clhs41 =             C(2,3)*DN(0,2) + clhs34 + clhs40;
const double clhs42 =             C(2,5)*DN(0,2);
const double clhs43 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs42;
const double clhs44 =             DN(0,2)*clhs37;
const double clhs45 =             clhs24*rho;
const double clhs46 =             clhs19*clhs45;
const double clhs47 =             -N[0] + clhs21*clhs24*clhs27 + clhs46;
const double clhs48 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs49 =             C(0,3)*DN(1,0);
const double clhs50 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs49;
const double clhs51 =             C(0,5)*DN(1,0);
const double clhs52 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs51;
const double clhs53 =             DN(0,0)*DN(1,0);
const double clhs54 =             N[0]*clhs18;
const double clhs55 =             N[1]*clhs54;
const double clhs56 =             clhs15*clhs53 + clhs55;
const double clhs57 =             DN(1,0)*clhs10 + DN(1,1)*clhs11 + DN(1,2)*clhs12;
const double clhs58 =             clhs57*theta;
const double clhs59 =             N[1]*clhs22 + clhs58;
const double clhs60 =             clhs21*clhs58 + clhs26*clhs59 + clhs29*clhs59;
const double clhs61 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs49;
const double clhs62 =             C(1,3)*DN(1,1);
const double clhs63 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs62;
const double clhs64 =             C(3,5)*DN(1,0);
const double clhs65 =             C(4,5)*DN(1,2);
const double clhs66 =             C(1,5)*DN(1,1) + clhs64 + clhs65;
const double clhs67 =             DN(1,1)*clhs37;
const double clhs68 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs51;
const double clhs69 =             C(3,4)*DN(1,1);
const double clhs70 =             C(2,3)*DN(1,2) + clhs64 + clhs69;
const double clhs71 =             C(2,5)*DN(1,2);
const double clhs72 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs71;
const double clhs73 =             DN(1,2)*clhs37;
const double clhs74 =             DN(0,0)*N[1];
const double clhs75 =             DN(1,0)*N[0];
const double clhs76 =             clhs27*clhs45;
const double clhs77 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs78 =             C(0,3)*DN(2,0);
const double clhs79 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs78;
const double clhs80 =             C(0,5)*DN(2,0);
const double clhs81 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs80;
const double clhs82 =             DN(0,0)*DN(2,0);
const double clhs83 =             N[2]*clhs54;
const double clhs84 =             clhs15*clhs82 + clhs83;
const double clhs85 =             DN(2,0)*clhs10 + DN(2,1)*clhs11 + DN(2,2)*clhs12;
const double clhs86 =             clhs85*theta;
const double clhs87 =             N[2]*clhs22 + clhs86;
const double clhs88 =             clhs21*clhs86 + clhs26*clhs87 + clhs29*clhs87;
const double clhs89 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs78;
const double clhs90 =             C(1,3)*DN(2,1);
const double clhs91 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs90;
const double clhs92 =             C(3,5)*DN(2,0);
const double clhs93 =             C(4,5)*DN(2,2);
const double clhs94 =             C(1,5)*DN(2,1) + clhs92 + clhs93;
const double clhs95 =             DN(2,1)*clhs37;
const double clhs96 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs80;
const double clhs97 =             C(3,4)*DN(2,1);
const double clhs98 =             C(2,3)*DN(2,2) + clhs92 + clhs97;
const double clhs99 =             C(2,5)*DN(2,2);
const double clhs100 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs99;
const double clhs101 =             DN(2,2)*clhs37;
const double clhs102 =             DN(0,0)*N[2];
const double clhs103 =             DN(2,0)*N[0];
const double clhs104 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs105 =             C(0,3)*DN(3,0);
const double clhs106 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs105;
const double clhs107 =             C(0,5)*DN(3,0);
const double clhs108 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs107;
const double clhs109 =             DN(0,0)*DN(3,0);
const double clhs110 =             N[3]*clhs54;
const double clhs111 =             clhs109*clhs15 + clhs110;
const double clhs112 =             DN(3,0)*clhs10 + DN(3,1)*clhs11 + DN(3,2)*clhs12;
const double clhs113 =             clhs112*theta;
const double clhs114 =             N[3]*clhs22 + clhs113;
const double clhs115 =             clhs113*clhs21 + clhs114*clhs26 + clhs114*clhs29;
const double clhs116 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs105;
const double clhs117 =             C(1,3)*DN(3,1);
const double clhs118 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs117;
const double clhs119 =             C(3,5)*DN(3,0);
const double clhs120 =             C(4,5)*DN(3,2);
const double clhs121 =             C(1,5)*DN(3,1) + clhs119 + clhs120;
const double clhs122 =             DN(3,1)*clhs37;
const double clhs123 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs107;
const double clhs124 =             C(3,4)*DN(3,1);
const double clhs125 =             C(2,3)*DN(3,2) + clhs119 + clhs124;
const double clhs126 =             C(2,5)*DN(3,2);
const double clhs127 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs126;
const double clhs128 =             DN(3,2)*clhs37;
const double clhs129 =             DN(0,0)*N[3];
const double clhs130 =             DN(3,0)*N[0];
const double clhs131 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs32;
const double clhs132 =             C(0,4)*DN(0,0) + clhs35 + clhs40;
const double clhs133 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs134 =             C(1,4)*DN(0,1);
const double clhs135 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs134;
const double clhs136 =             pow(DN(0,1), 2);
const double clhs137 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs134;
const double clhs138 =             C(2,4)*DN(0,2);
const double clhs139 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs138;
const double clhs140 =             DN(0,1)*clhs14;
const double clhs141 =             DN(0,2)*clhs140;
const double clhs142 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs62;
const double clhs143 =             C(0,4)*DN(1,0) + clhs65 + clhs69;
const double clhs144 =             DN(1,0)*clhs140;
const double clhs145 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs146 =             C(1,4)*DN(1,1);
const double clhs147 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs146;
const double clhs148 =             DN(0,1)*DN(1,1);
const double clhs149 =             clhs148*clhs15;
const double clhs150 =             clhs55 + clhs60;
const double clhs151 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs146;
const double clhs152 =             C(2,4)*DN(1,2);
const double clhs153 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs152;
const double clhs154 =             DN(1,2)*clhs140;
const double clhs155 =             DN(0,1)*N[1];
const double clhs156 =             DN(1,1)*N[0];
const double clhs157 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs90;
const double clhs158 =             C(0,4)*DN(2,0) + clhs93 + clhs97;
const double clhs159 =             DN(2,0)*clhs140;
const double clhs160 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs161 =             C(1,4)*DN(2,1);
const double clhs162 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs161;
const double clhs163 =             DN(0,1)*DN(2,1);
const double clhs164 =             clhs15*clhs163;
const double clhs165 =             clhs83 + clhs88;
const double clhs166 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs161;
const double clhs167 =             C(2,4)*DN(2,2);
const double clhs168 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs167;
const double clhs169 =             DN(2,2)*clhs140;
const double clhs170 =             DN(0,1)*N[2];
const double clhs171 =             DN(2,1)*N[0];
const double clhs172 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs117;
const double clhs173 =             C(0,4)*DN(3,0) + clhs120 + clhs124;
const double clhs174 =             DN(3,0)*clhs140;
const double clhs175 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs176 =             C(1,4)*DN(3,1);
const double clhs177 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs176;
const double clhs178 =             DN(0,1)*DN(3,1);
const double clhs179 =             clhs15*clhs178;
const double clhs180 =             clhs110 + clhs115;
const double clhs181 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs176;
const double clhs182 =             C(2,4)*DN(3,2);
const double clhs183 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs182;
const double clhs184 =             DN(3,2)*clhs140;
const double clhs185 =             DN(0,1)*N[3];
const double clhs186 =             DN(3,1)*N[0];
const double clhs187 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs42;
const double clhs188 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs138;
const double clhs189 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs190 =             pow(DN(0,2), 2);
const double clhs191 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs71;
const double clhs192 =             DN(0,2)*clhs14;
const double clhs193 =             DN(1,0)*clhs192;
const double clhs194 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs152;
const double clhs195 =             DN(1,1)*clhs192;
const double clhs196 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs197 =             DN(0,2)*DN(1,2);
const double clhs198 =             clhs15*clhs197;
const double clhs199 =             DN(0,2)*N[1];
const double clhs200 =             DN(1,2)*N[0];
const double clhs201 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs99;
const double clhs202 =             DN(2,0)*clhs192;
const double clhs203 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs167;
const double clhs204 =             DN(2,1)*clhs192;
const double clhs205 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs206 =             DN(0,2)*DN(2,2);
const double clhs207 =             clhs15*clhs206;
const double clhs208 =             DN(0,2)*N[2];
const double clhs209 =             DN(2,2)*N[0];
const double clhs210 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs126;
const double clhs211 =             DN(3,0)*clhs192;
const double clhs212 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs182;
const double clhs213 =             DN(3,1)*clhs192;
const double clhs214 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs215 =             DN(0,2)*DN(3,2);
const double clhs216 =             clhs15*clhs215;
const double clhs217 =             DN(0,2)*N[3];
const double clhs218 =             DN(3,2)*N[0];
const double clhs219 =             N[0]*theta;
const double clhs220 =             clhs23*clhs45;
const double clhs221 =             clhs219 + clhs220;
const double clhs222 =             clhs45*clhs59;
const double clhs223 =             clhs24*(clhs148 + clhs197 + clhs53);
const double clhs224 =             clhs45*clhs87;
const double clhs225 =             clhs24*(clhs163 + clhs206 + clhs82);
const double clhs226 =             clhs114*clhs45;
const double clhs227 =             clhs24*(clhs109 + clhs178 + clhs215);
const double clhs228 =             DN(1,0)*theta;
const double clhs229 =             DN(1,1)*theta;
const double clhs230 =             DN(1,2)*theta;
const double clhs231 =             N[1]*rho;
const double clhs232 =             clhs25*clhs57;
const double clhs233 =             N[1]*clhs28;
const double clhs234 =             clhs20*clhs231 + clhs23*clhs232 + clhs23*clhs233;
const double clhs235 =             clhs45*clhs57;
const double clhs236 =             pow(DN(1,0), 2);
const double clhs237 =             pow(N[1], 2)*clhs18 + clhs231*clhs58 + clhs232*clhs59 + clhs233*clhs59;
const double clhs238 =             DN(1,0)*clhs14;
const double clhs239 =             DN(1,1)*clhs238;
const double clhs240 =             DN(1,2)*clhs238;
const double clhs241 =             N[1]*clhs76 - N[1] + clhs235;
const double clhs242 =             DN(1,0)*DN(2,0);
const double clhs243 =             N[1]*clhs18;
const double clhs244 =             N[2]*clhs243;
const double clhs245 =             clhs15*clhs242 + clhs244;
const double clhs246 =             clhs231*clhs86 + clhs232*clhs87 + clhs233*clhs87;
const double clhs247 =             DN(2,1)*clhs238;
const double clhs248 =             DN(2,2)*clhs238;
const double clhs249 =             DN(1,0)*N[2];
const double clhs250 =             DN(2,0)*N[1];
const double clhs251 =             DN(1,0)*DN(3,0);
const double clhs252 =             N[3]*clhs243;
const double clhs253 =             clhs15*clhs251 + clhs252;
const double clhs254 =             clhs113*clhs231 + clhs114*clhs232 + clhs114*clhs233;
const double clhs255 =             DN(3,1)*clhs238;
const double clhs256 =             DN(3,2)*clhs238;
const double clhs257 =             DN(1,0)*N[3];
const double clhs258 =             DN(3,0)*N[1];
const double clhs259 =             clhs234 + clhs55;
const double clhs260 =             pow(DN(1,1), 2);
const double clhs261 =             DN(1,1)*clhs14;
const double clhs262 =             DN(1,2)*clhs261;
const double clhs263 =             DN(2,0)*clhs261;
const double clhs264 =             DN(1,1)*DN(2,1);
const double clhs265 =             clhs15*clhs264;
const double clhs266 =             clhs244 + clhs246;
const double clhs267 =             DN(2,2)*clhs261;
const double clhs268 =             DN(1,1)*N[2];
const double clhs269 =             DN(2,1)*N[1];
const double clhs270 =             DN(3,0)*clhs261;
const double clhs271 =             DN(1,1)*DN(3,1);
const double clhs272 =             clhs15*clhs271;
const double clhs273 =             clhs252 + clhs254;
const double clhs274 =             DN(3,2)*clhs261;
const double clhs275 =             DN(1,1)*N[3];
const double clhs276 =             DN(3,1)*N[1];
const double clhs277 =             pow(DN(1,2), 2);
const double clhs278 =             DN(1,2)*clhs14;
const double clhs279 =             DN(2,0)*clhs278;
const double clhs280 =             DN(2,1)*clhs278;
const double clhs281 =             DN(1,2)*DN(2,2);
const double clhs282 =             clhs15*clhs281;
const double clhs283 =             DN(1,2)*N[2];
const double clhs284 =             DN(2,2)*N[1];
const double clhs285 =             DN(3,0)*clhs278;
const double clhs286 =             DN(3,1)*clhs278;
const double clhs287 =             DN(1,2)*DN(3,2);
const double clhs288 =             clhs15*clhs287;
const double clhs289 =             DN(1,2)*N[3];
const double clhs290 =             DN(3,2)*N[1];
const double clhs291 =             N[1]*theta;
const double clhs292 =             clhs222 + clhs291;
const double clhs293 =             clhs24*(clhs242 + clhs264 + clhs281);
const double clhs294 =             clhs24*(clhs251 + clhs271 + clhs287);
const double clhs295 =             DN(2,0)*theta;
const double clhs296 =             DN(2,1)*theta;
const double clhs297 =             DN(2,2)*theta;
const double clhs298 =             N[2]*rho;
const double clhs299 =             clhs25*clhs85;
const double clhs300 =             N[2]*clhs28;
const double clhs301 =             clhs20*clhs298 + clhs23*clhs299 + clhs23*clhs300;
const double clhs302 =             clhs45*clhs85;
const double clhs303 =             clhs298*clhs58 + clhs299*clhs59 + clhs300*clhs59;
const double clhs304 =             pow(DN(2,0), 2);
const double clhs305 =             pow(N[2], 2)*clhs18 + clhs298*clhs86 + clhs299*clhs87 + clhs300*clhs87;
const double clhs306 =             DN(2,0)*clhs14;
const double clhs307 =             DN(2,1)*clhs306;
const double clhs308 =             DN(2,2)*clhs306;
const double clhs309 =             N[2]*clhs76 - N[2] + clhs302;
const double clhs310 =             DN(2,0)*DN(3,0);
const double clhs311 =             N[2]*N[3]*clhs18;
const double clhs312 =             clhs15*clhs310 + clhs311;
const double clhs313 =             clhs113*clhs298 + clhs114*clhs299 + clhs114*clhs300;
const double clhs314 =             DN(3,1)*clhs306;
const double clhs315 =             DN(3,2)*clhs306;
const double clhs316 =             DN(2,0)*N[3];
const double clhs317 =             DN(3,0)*N[2];
const double clhs318 =             clhs301 + clhs83;
const double clhs319 =             clhs244 + clhs303;
const double clhs320 =             pow(DN(2,1), 2);
const double clhs321 =             DN(2,1)*clhs14;
const double clhs322 =             DN(2,2)*clhs321;
const double clhs323 =             DN(3,0)*clhs321;
const double clhs324 =             DN(2,1)*DN(3,1);
const double clhs325 =             clhs15*clhs324;
const double clhs326 =             clhs311 + clhs313;
const double clhs327 =             DN(3,2)*clhs321;
const double clhs328 =             DN(2,1)*N[3];
const double clhs329 =             DN(3,1)*N[2];
const double clhs330 =             pow(DN(2,2), 2);
const double clhs331 =             DN(2,2)*clhs14;
const double clhs332 =             DN(3,0)*clhs331;
const double clhs333 =             DN(3,1)*clhs331;
const double clhs334 =             DN(2,2)*DN(3,2);
const double clhs335 =             clhs15*clhs334;
const double clhs336 =             DN(2,2)*N[3];
const double clhs337 =             DN(3,2)*N[2];
const double clhs338 =             N[2]*theta;
const double clhs339 =             clhs224 + clhs338;
const double clhs340 =             clhs24*(clhs310 + clhs324 + clhs334);
const double clhs341 =             DN(3,0)*theta;
const double clhs342 =             DN(3,1)*theta;
const double clhs343 =             DN(3,2)*theta;
const double clhs344 =             N[3]*rho;
const double clhs345 =             clhs112*clhs25;
const double clhs346 =             N[3]*clhs28;
const double clhs347 =             clhs20*clhs344 + clhs23*clhs345 + clhs23*clhs346;
const double clhs348 =             clhs112*clhs45;
const double clhs349 =             clhs344*clhs58 + clhs345*clhs59 + clhs346*clhs59;
const double clhs350 =             clhs344*clhs86 + clhs345*clhs87 + clhs346*clhs87;
const double clhs351 =             pow(DN(3,0), 2);
const double clhs352 =             pow(N[3], 2)*clhs18 + clhs113*clhs344 + clhs114*clhs345 + clhs114*clhs346;
const double clhs353 =             DN(3,0)*clhs14;
const double clhs354 =             DN(3,1)*clhs353;
const double clhs355 =             DN(3,2)*clhs353;
const double clhs356 =             N[3]*clhs76 - N[3] + clhs348;
const double clhs357 =             clhs110 + clhs347;
const double clhs358 =             clhs252 + clhs349;
const double clhs359 =             clhs311 + clhs350;
const double clhs360 =             pow(DN(3,1), 2);
const double clhs361 =             DN(3,1)*DN(3,2)*clhs14;
const double clhs362 =             pow(DN(3,2), 2);
const double clhs363 =             N[3]*theta + clhs226;
            lhs(0,0)=clhs0*clhs1 + clhs15*clhs8 + clhs3*clhs4 + clhs30 + clhs6*clhs7;
            lhs(0,1)=theta*(DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs36 + clhs38);
            lhs(0,2)=theta*(DN(0,0)*clhs39 + DN(0,1)*clhs41 + DN(0,2)*clhs43 + clhs44);
            lhs(0,3)=DN(0,0)*clhs47;
            lhs(0,4)=clhs1*clhs48 + clhs4*clhs50 + clhs52*clhs7 + clhs56 + clhs60;
            lhs(0,5)=theta*(DN(0,0)*clhs61 + DN(0,1)*clhs63 + DN(0,2)*clhs66 + clhs67);
            lhs(0,6)=theta*(DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs73);
            lhs(0,7)=DN(1,0)*clhs46 - clhs74 + clhs75*clhs76;
            lhs(0,8)=clhs1*clhs77 + clhs4*clhs79 + clhs7*clhs81 + clhs84 + clhs88;
            lhs(0,9)=theta*(DN(0,0)*clhs89 + DN(0,1)*clhs91 + DN(0,2)*clhs94 + clhs95);
            lhs(0,10)=theta*(DN(0,0)*clhs96 + DN(0,1)*clhs98 + DN(0,2)*clhs100 + clhs101);
            lhs(0,11)=DN(2,0)*clhs46 - clhs102 + clhs103*clhs76;
            lhs(0,12)=clhs1*clhs104 + clhs106*clhs4 + clhs108*clhs7 + clhs111 + clhs115;
            lhs(0,13)=theta*(DN(0,0)*clhs116 + DN(0,1)*clhs118 + DN(0,2)*clhs121 + clhs122);
            lhs(0,14)=theta*(DN(0,0)*clhs123 + DN(0,1)*clhs125 + DN(0,2)*clhs127 + clhs128);
            lhs(0,15)=DN(3,0)*clhs46 - clhs129 + clhs130*clhs76;
            lhs(1,0)=theta*(DN(0,0)*clhs3 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs38);
            lhs(1,1)=clhs1*clhs33 + clhs133*clhs4 + clhs135*clhs7 + clhs136*clhs15 + clhs30;
            lhs(1,2)=theta*(DN(0,0)*clhs41 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs141);
            lhs(1,3)=DN(0,1)*clhs47;
            lhs(1,4)=theta*(DN(0,0)*clhs50 + DN(0,1)*clhs142 + DN(0,2)*clhs143 + clhs144);
            lhs(1,5)=clhs1*clhs63 + clhs145*clhs4 + clhs147*clhs7 + clhs149 + clhs150;
            lhs(1,6)=theta*(DN(0,0)*clhs70 + DN(0,1)*clhs151 + DN(0,2)*clhs153 + clhs154);
            lhs(1,7)=DN(1,1)*clhs46 - clhs155 + clhs156*clhs76;
            lhs(1,8)=theta*(DN(0,0)*clhs79 + DN(0,1)*clhs157 + DN(0,2)*clhs158 + clhs159);
            lhs(1,9)=clhs1*clhs91 + clhs160*clhs4 + clhs162*clhs7 + clhs164 + clhs165;
            lhs(1,10)=theta*(DN(0,0)*clhs98 + DN(0,1)*clhs166 + DN(0,2)*clhs168 + clhs169);
            lhs(1,11)=DN(2,1)*clhs46 - clhs170 + clhs171*clhs76;
            lhs(1,12)=theta*(DN(0,0)*clhs106 + DN(0,1)*clhs172 + DN(0,2)*clhs173 + clhs174);
            lhs(1,13)=clhs1*clhs118 + clhs175*clhs4 + clhs177*clhs7 + clhs179 + clhs180;
            lhs(1,14)=theta*(DN(0,0)*clhs125 + DN(0,1)*clhs181 + DN(0,2)*clhs183 + clhs184);
            lhs(1,15)=DN(3,1)*clhs46 - clhs185 + clhs186*clhs76;
            lhs(2,0)=theta*(DN(0,0)*clhs6 + DN(0,1)*clhs132 + DN(0,2)*clhs187 + clhs44);
            lhs(2,1)=theta*(DN(0,0)*clhs36 + DN(0,1)*clhs135 + DN(0,2)*clhs188 + clhs141);
            lhs(2,2)=clhs1*clhs43 + clhs139*clhs4 + clhs15*clhs190 + clhs189*clhs7 + clhs30;
            lhs(2,3)=DN(0,2)*clhs47;
            lhs(2,4)=theta*(DN(0,0)*clhs52 + DN(0,1)*clhs143 + DN(0,2)*clhs191 + clhs193);
            lhs(2,5)=theta*(DN(0,0)*clhs66 + DN(0,1)*clhs147 + DN(0,2)*clhs194 + clhs195);
            lhs(2,6)=clhs1*clhs72 + clhs150 + clhs153*clhs4 + clhs196*clhs7 + clhs198;
            lhs(2,7)=DN(1,2)*clhs46 - clhs199 + clhs200*clhs76;
            lhs(2,8)=theta*(DN(0,0)*clhs81 + DN(0,1)*clhs158 + DN(0,2)*clhs201 + clhs202);
            lhs(2,9)=theta*(DN(0,0)*clhs94 + DN(0,1)*clhs162 + DN(0,2)*clhs203 + clhs204);
            lhs(2,10)=clhs1*clhs100 + clhs165 + clhs168*clhs4 + clhs205*clhs7 + clhs207;
            lhs(2,11)=DN(2,2)*clhs46 - clhs208 + clhs209*clhs76;
            lhs(2,12)=theta*(DN(0,0)*clhs108 + DN(0,1)*clhs173 + DN(0,2)*clhs210 + clhs211);
            lhs(2,13)=theta*(DN(0,0)*clhs121 + DN(0,1)*clhs177 + DN(0,2)*clhs212 + clhs213);
            lhs(2,14)=clhs1*clhs127 + clhs180 + clhs183*clhs4 + clhs214*clhs7 + clhs216;
            lhs(2,15)=DN(3,2)*clhs46 - clhs217 + clhs218*clhs76;
            lhs(3,0)=DN(0,0)*clhs221;
            lhs(3,1)=DN(0,1)*clhs221;
            lhs(3,2)=DN(0,2)*clhs221;
            lhs(3,3)=clhs24*(clhs136 + clhs190 + clhs8);
            lhs(3,4)=DN(0,0)*clhs222 + DN(1,0)*clhs219;
            lhs(3,5)=DN(0,1)*clhs222 + DN(1,1)*clhs219;
            lhs(3,6)=DN(0,2)*clhs222 + DN(1,2)*clhs219;
            lhs(3,7)=clhs223;
            lhs(3,8)=DN(0,0)*clhs224 + DN(2,0)*clhs219;
            lhs(3,9)=DN(0,1)*clhs224 + DN(2,1)*clhs219;
            lhs(3,10)=DN(0,2)*clhs224 + DN(2,2)*clhs219;
            lhs(3,11)=clhs225;
            lhs(3,12)=DN(0,0)*clhs226 + DN(3,0)*clhs219;
            lhs(3,13)=DN(0,1)*clhs226 + DN(3,1)*clhs219;
            lhs(3,14)=DN(0,2)*clhs226 + DN(3,2)*clhs219;
            lhs(3,15)=clhs227;
            lhs(4,0)=clhs0*clhs228 + clhs229*clhs3 + clhs230*clhs6 + clhs234 + clhs56;
            lhs(4,1)=theta*(DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs36 + clhs144);
            lhs(4,2)=theta*(DN(1,0)*clhs39 + DN(1,1)*clhs41 + DN(1,2)*clhs43 + clhs193);
            lhs(4,3)=DN(0,0)*clhs235 + clhs74*clhs76 - clhs75;
            lhs(4,4)=clhs15*clhs236 + clhs228*clhs48 + clhs229*clhs50 + clhs230*clhs52 + clhs237;
            lhs(4,5)=theta*(DN(1,0)*clhs61 + DN(1,1)*clhs63 + DN(1,2)*clhs66 + clhs239);
            lhs(4,6)=theta*(DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs240);
            lhs(4,7)=DN(1,0)*clhs241;
            lhs(4,8)=clhs228*clhs77 + clhs229*clhs79 + clhs230*clhs81 + clhs245 + clhs246;
            lhs(4,9)=theta*(DN(1,0)*clhs89 + DN(1,1)*clhs91 + DN(1,2)*clhs94 + clhs247);
            lhs(4,10)=theta*(DN(1,0)*clhs96 + DN(1,1)*clhs98 + DN(1,2)*clhs100 + clhs248);
            lhs(4,11)=DN(2,0)*clhs235 - clhs249 + clhs250*clhs76;
            lhs(4,12)=clhs104*clhs228 + clhs106*clhs229 + clhs108*clhs230 + clhs253 + clhs254;
            lhs(4,13)=theta*(DN(1,0)*clhs116 + DN(1,1)*clhs118 + DN(1,2)*clhs121 + clhs255);
            lhs(4,14)=theta*(DN(1,0)*clhs123 + DN(1,1)*clhs125 + DN(1,2)*clhs127 + clhs256);
            lhs(4,15)=DN(3,0)*clhs235 - clhs257 + clhs258*clhs76;
            lhs(5,0)=theta*(DN(1,0)*clhs3 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs67);
            lhs(5,1)=clhs133*clhs229 + clhs135*clhs230 + clhs149 + clhs228*clhs33 + clhs259;
            lhs(5,2)=theta*(DN(1,0)*clhs41 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs195);
            lhs(5,3)=DN(0,1)*clhs235 + clhs155*clhs76 - clhs156;
            lhs(5,4)=theta*(DN(1,0)*clhs50 + DN(1,1)*clhs142 + DN(1,2)*clhs143 + clhs239);
            lhs(5,5)=clhs145*clhs229 + clhs147*clhs230 + clhs15*clhs260 + clhs228*clhs63 + clhs237;
            lhs(5,6)=theta*(DN(1,0)*clhs70 + DN(1,1)*clhs151 + DN(1,2)*clhs153 + clhs262);
            lhs(5,7)=DN(1,1)*clhs241;
            lhs(5,8)=theta*(DN(1,0)*clhs79 + DN(1,1)*clhs157 + DN(1,2)*clhs158 + clhs263);
            lhs(5,9)=clhs160*clhs229 + clhs162*clhs230 + clhs228*clhs91 + clhs265 + clhs266;
            lhs(5,10)=theta*(DN(1,0)*clhs98 + DN(1,1)*clhs166 + DN(1,2)*clhs168 + clhs267);
            lhs(5,11)=DN(2,1)*clhs235 - clhs268 + clhs269*clhs76;
            lhs(5,12)=theta*(DN(1,0)*clhs106 + DN(1,1)*clhs172 + DN(1,2)*clhs173 + clhs270);
            lhs(5,13)=clhs118*clhs228 + clhs175*clhs229 + clhs177*clhs230 + clhs272 + clhs273;
            lhs(5,14)=theta*(DN(1,0)*clhs125 + DN(1,1)*clhs181 + DN(1,2)*clhs183 + clhs274);
            lhs(5,15)=DN(3,1)*clhs235 - clhs275 + clhs276*clhs76;
            lhs(6,0)=theta*(DN(1,0)*clhs6 + DN(1,1)*clhs132 + DN(1,2)*clhs187 + clhs73);
            lhs(6,1)=theta*(DN(1,0)*clhs36 + DN(1,1)*clhs135 + DN(1,2)*clhs188 + clhs154);
            lhs(6,2)=clhs139*clhs229 + clhs189*clhs230 + clhs198 + clhs228*clhs43 + clhs259;
            lhs(6,3)=DN(0,2)*clhs235 + clhs199*clhs76 - clhs200;
            lhs(6,4)=theta*(DN(1,0)*clhs52 + DN(1,1)*clhs143 + DN(1,2)*clhs191 + clhs240);
            lhs(6,5)=theta*(DN(1,0)*clhs66 + DN(1,1)*clhs147 + DN(1,2)*clhs194 + clhs262);
            lhs(6,6)=clhs15*clhs277 + clhs153*clhs229 + clhs196*clhs230 + clhs228*clhs72 + clhs237;
            lhs(6,7)=DN(1,2)*clhs241;
            lhs(6,8)=theta*(DN(1,0)*clhs81 + DN(1,1)*clhs158 + DN(1,2)*clhs201 + clhs279);
            lhs(6,9)=theta*(DN(1,0)*clhs94 + DN(1,1)*clhs162 + DN(1,2)*clhs203 + clhs280);
            lhs(6,10)=clhs100*clhs228 + clhs168*clhs229 + clhs205*clhs230 + clhs266 + clhs282;
            lhs(6,11)=DN(2,2)*clhs235 - clhs283 + clhs284*clhs76;
            lhs(6,12)=theta*(DN(1,0)*clhs108 + DN(1,1)*clhs173 + DN(1,2)*clhs210 + clhs285);
            lhs(6,13)=theta*(DN(1,0)*clhs121 + DN(1,1)*clhs177 + DN(1,2)*clhs212 + clhs286);
            lhs(6,14)=clhs127*clhs228 + clhs183*clhs229 + clhs214*clhs230 + clhs273 + clhs288;
            lhs(6,15)=DN(3,2)*clhs235 - clhs289 + clhs290*clhs76;
            lhs(7,0)=DN(1,0)*clhs220 + clhs74*theta;
            lhs(7,1)=DN(1,1)*clhs220 + clhs155*theta;
            lhs(7,2)=DN(1,2)*clhs220 + clhs199*theta;
            lhs(7,3)=clhs223;
            lhs(7,4)=DN(1,0)*clhs292;
            lhs(7,5)=DN(1,1)*clhs292;
            lhs(7,6)=DN(1,2)*clhs292;
            lhs(7,7)=clhs24*(clhs236 + clhs260 + clhs277);
            lhs(7,8)=DN(1,0)*clhs224 + DN(2,0)*clhs291;
            lhs(7,9)=DN(1,1)*clhs224 + DN(2,1)*clhs291;
            lhs(7,10)=DN(1,2)*clhs224 + DN(2,2)*clhs291;
            lhs(7,11)=clhs293;
            lhs(7,12)=DN(1,0)*clhs226 + DN(3,0)*clhs291;
            lhs(7,13)=DN(1,1)*clhs226 + DN(3,1)*clhs291;
            lhs(7,14)=DN(1,2)*clhs226 + DN(3,2)*clhs291;
            lhs(7,15)=clhs294;
            lhs(8,0)=clhs0*clhs295 + clhs296*clhs3 + clhs297*clhs6 + clhs301 + clhs84;
            lhs(8,1)=theta*(DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs36 + clhs159);
            lhs(8,2)=theta*(DN(2,0)*clhs39 + DN(2,1)*clhs41 + DN(2,2)*clhs43 + clhs202);
            lhs(8,3)=DN(0,0)*clhs302 + clhs102*clhs76 - clhs103;
            lhs(8,4)=clhs245 + clhs295*clhs48 + clhs296*clhs50 + clhs297*clhs52 + clhs303;
            lhs(8,5)=theta*(DN(2,0)*clhs61 + DN(2,1)*clhs63 + DN(2,2)*clhs66 + clhs263);
            lhs(8,6)=theta*(DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs279);
            lhs(8,7)=DN(1,0)*clhs302 + clhs249*clhs76 - clhs250;
            lhs(8,8)=clhs15*clhs304 + clhs295*clhs77 + clhs296*clhs79 + clhs297*clhs81 + clhs305;
            lhs(8,9)=theta*(DN(2,0)*clhs89 + DN(2,1)*clhs91 + DN(2,2)*clhs94 + clhs307);
            lhs(8,10)=theta*(DN(2,0)*clhs96 + DN(2,1)*clhs98 + DN(2,2)*clhs100 + clhs308);
            lhs(8,11)=DN(2,0)*clhs309;
            lhs(8,12)=clhs104*clhs295 + clhs106*clhs296 + clhs108*clhs297 + clhs312 + clhs313;
            lhs(8,13)=theta*(DN(2,0)*clhs116 + DN(2,1)*clhs118 + DN(2,2)*clhs121 + clhs314);
            lhs(8,14)=theta*(DN(2,0)*clhs123 + DN(2,1)*clhs125 + DN(2,2)*clhs127 + clhs315);
            lhs(8,15)=DN(3,0)*clhs302 - clhs316 + clhs317*clhs76;
            lhs(9,0)=theta*(DN(2,0)*clhs3 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs95);
            lhs(9,1)=clhs133*clhs296 + clhs135*clhs297 + clhs164 + clhs295*clhs33 + clhs318;
            lhs(9,2)=theta*(DN(2,0)*clhs41 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs204);
            lhs(9,3)=DN(0,1)*clhs302 + clhs170*clhs76 - clhs171;
            lhs(9,4)=theta*(DN(2,0)*clhs50 + DN(2,1)*clhs142 + DN(2,2)*clhs143 + clhs247);
            lhs(9,5)=clhs145*clhs296 + clhs147*clhs297 + clhs265 + clhs295*clhs63 + clhs319;
            lhs(9,6)=theta*(DN(2,0)*clhs70 + DN(2,1)*clhs151 + DN(2,2)*clhs153 + clhs280);
            lhs(9,7)=DN(1,1)*clhs302 + clhs268*clhs76 - clhs269;
            lhs(9,8)=theta*(DN(2,0)*clhs79 + DN(2,1)*clhs157 + DN(2,2)*clhs158 + clhs307);
            lhs(9,9)=clhs15*clhs320 + clhs160*clhs296 + clhs162*clhs297 + clhs295*clhs91 + clhs305;
            lhs(9,10)=theta*(DN(2,0)*clhs98 + DN(2,1)*clhs166 + DN(2,2)*clhs168 + clhs322);
            lhs(9,11)=DN(2,1)*clhs309;
            lhs(9,12)=theta*(DN(2,0)*clhs106 + DN(2,1)*clhs172 + DN(2,2)*clhs173 + clhs323);
            lhs(9,13)=clhs118*clhs295 + clhs175*clhs296 + clhs177*clhs297 + clhs325 + clhs326;
            lhs(9,14)=theta*(DN(2,0)*clhs125 + DN(2,1)*clhs181 + DN(2,2)*clhs183 + clhs327);
            lhs(9,15)=DN(3,1)*clhs302 - clhs328 + clhs329*clhs76;
            lhs(10,0)=theta*(DN(2,0)*clhs6 + DN(2,1)*clhs132 + DN(2,2)*clhs187 + clhs101);
            lhs(10,1)=theta*(DN(2,0)*clhs36 + DN(2,1)*clhs135 + DN(2,2)*clhs188 + clhs169);
            lhs(10,2)=clhs139*clhs296 + clhs189*clhs297 + clhs207 + clhs295*clhs43 + clhs318;
            lhs(10,3)=DN(0,2)*clhs302 + clhs208*clhs76 - clhs209;
            lhs(10,4)=theta*(DN(2,0)*clhs52 + DN(2,1)*clhs143 + DN(2,2)*clhs191 + clhs248);
            lhs(10,5)=theta*(DN(2,0)*clhs66 + DN(2,1)*clhs147 + DN(2,2)*clhs194 + clhs267);
            lhs(10,6)=clhs153*clhs296 + clhs196*clhs297 + clhs282 + clhs295*clhs72 + clhs319;
            lhs(10,7)=DN(1,2)*clhs302 + clhs283*clhs76 - clhs284;
            lhs(10,8)=theta*(DN(2,0)*clhs81 + DN(2,1)*clhs158 + DN(2,2)*clhs201 + clhs308);
            lhs(10,9)=theta*(DN(2,0)*clhs94 + DN(2,1)*clhs162 + DN(2,2)*clhs203 + clhs322);
            lhs(10,10)=clhs100*clhs295 + clhs15*clhs330 + clhs168*clhs296 + clhs205*clhs297 + clhs305;
            lhs(10,11)=DN(2,2)*clhs309;
            lhs(10,12)=theta*(DN(2,0)*clhs108 + DN(2,1)*clhs173 + DN(2,2)*clhs210 + clhs332);
            lhs(10,13)=theta*(DN(2,0)*clhs121 + DN(2,1)*clhs177 + DN(2,2)*clhs212 + clhs333);
            lhs(10,14)=clhs127*clhs295 + clhs183*clhs296 + clhs214*clhs297 + clhs326 + clhs335;
            lhs(10,15)=DN(3,2)*clhs302 - clhs336 + clhs337*clhs76;
            lhs(11,0)=DN(2,0)*clhs220 + clhs102*theta;
            lhs(11,1)=DN(2,1)*clhs220 + clhs170*theta;
            lhs(11,2)=DN(2,2)*clhs220 + clhs208*theta;
            lhs(11,3)=clhs225;
            lhs(11,4)=DN(2,0)*clhs222 + clhs249*theta;
            lhs(11,5)=DN(2,1)*clhs222 + clhs268*theta;
            lhs(11,6)=DN(2,2)*clhs222 + clhs283*theta;
            lhs(11,7)=clhs293;
            lhs(11,8)=DN(2,0)*clhs339;
            lhs(11,9)=DN(2,1)*clhs339;
            lhs(11,10)=DN(2,2)*clhs339;
            lhs(11,11)=clhs24*(clhs304 + clhs320 + clhs330);
            lhs(11,12)=DN(2,0)*clhs226 + DN(3,0)*clhs338;
            lhs(11,13)=DN(2,1)*clhs226 + DN(3,1)*clhs338;
            lhs(11,14)=DN(2,2)*clhs226 + DN(3,2)*clhs338;
            lhs(11,15)=clhs340;
            lhs(12,0)=clhs0*clhs341 + clhs111 + clhs3*clhs342 + clhs343*clhs6 + clhs347;
            lhs(12,1)=theta*(DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs36 + clhs174);
            lhs(12,2)=theta*(DN(3,0)*clhs39 + DN(3,1)*clhs41 + DN(3,2)*clhs43 + clhs211);
            lhs(12,3)=DN(0,0)*clhs348 + clhs129*clhs76 - clhs130;
            lhs(12,4)=clhs253 + clhs341*clhs48 + clhs342*clhs50 + clhs343*clhs52 + clhs349;
            lhs(12,5)=theta*(DN(3,0)*clhs61 + DN(3,1)*clhs63 + DN(3,2)*clhs66 + clhs270);
            lhs(12,6)=theta*(DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs285);
            lhs(12,7)=DN(1,0)*clhs348 + clhs257*clhs76 - clhs258;
            lhs(12,8)=clhs312 + clhs341*clhs77 + clhs342*clhs79 + clhs343*clhs81 + clhs350;
            lhs(12,9)=theta*(DN(3,0)*clhs89 + DN(3,1)*clhs91 + DN(3,2)*clhs94 + clhs323);
            lhs(12,10)=theta*(DN(3,0)*clhs96 + DN(3,1)*clhs98 + DN(3,2)*clhs100 + clhs332);
            lhs(12,11)=DN(2,0)*clhs348 + clhs316*clhs76 - clhs317;
            lhs(12,12)=clhs104*clhs341 + clhs106*clhs342 + clhs108*clhs343 + clhs15*clhs351 + clhs352;
            lhs(12,13)=theta*(DN(3,0)*clhs116 + DN(3,1)*clhs118 + DN(3,2)*clhs121 + clhs354);
            lhs(12,14)=theta*(DN(3,0)*clhs123 + DN(3,1)*clhs125 + DN(3,2)*clhs127 + clhs355);
            lhs(12,15)=DN(3,0)*clhs356;
            lhs(13,0)=theta*(DN(3,0)*clhs3 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs122);
            lhs(13,1)=clhs133*clhs342 + clhs135*clhs343 + clhs179 + clhs33*clhs341 + clhs357;
            lhs(13,2)=theta*(DN(3,0)*clhs41 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs213);
            lhs(13,3)=DN(0,1)*clhs348 + clhs185*clhs76 - clhs186;
            lhs(13,4)=theta*(DN(3,0)*clhs50 + DN(3,1)*clhs142 + DN(3,2)*clhs143 + clhs255);
            lhs(13,5)=clhs145*clhs342 + clhs147*clhs343 + clhs272 + clhs341*clhs63 + clhs358;
            lhs(13,6)=theta*(DN(3,0)*clhs70 + DN(3,1)*clhs151 + DN(3,2)*clhs153 + clhs286);
            lhs(13,7)=DN(1,1)*clhs348 + clhs275*clhs76 - clhs276;
            lhs(13,8)=theta*(DN(3,0)*clhs79 + DN(3,1)*clhs157 + DN(3,2)*clhs158 + clhs314);
            lhs(13,9)=clhs160*clhs342 + clhs162*clhs343 + clhs325 + clhs341*clhs91 + clhs359;
            lhs(13,10)=theta*(DN(3,0)*clhs98 + DN(3,1)*clhs166 + DN(3,2)*clhs168 + clhs333);
            lhs(13,11)=DN(2,1)*clhs348 + clhs328*clhs76 - clhs329;
            lhs(13,12)=theta*(DN(3,0)*clhs106 + DN(3,1)*clhs172 + DN(3,2)*clhs173 + clhs354);
            lhs(13,13)=clhs118*clhs341 + clhs15*clhs360 + clhs175*clhs342 + clhs177*clhs343 + clhs352;
            lhs(13,14)=theta*(DN(3,0)*clhs125 + DN(3,1)*clhs181 + DN(3,2)*clhs183 + clhs361);
            lhs(13,15)=DN(3,1)*clhs356;
            lhs(14,0)=theta*(DN(3,0)*clhs6 + DN(3,1)*clhs132 + DN(3,2)*clhs187 + clhs128);
            lhs(14,1)=theta*(DN(3,0)*clhs36 + DN(3,1)*clhs135 + DN(3,2)*clhs188 + clhs184);
            lhs(14,2)=clhs139*clhs342 + clhs189*clhs343 + clhs216 + clhs341*clhs43 + clhs357;
            lhs(14,3)=DN(0,2)*clhs348 + clhs217*clhs76 - clhs218;
            lhs(14,4)=theta*(DN(3,0)*clhs52 + DN(3,1)*clhs143 + DN(3,2)*clhs191 + clhs256);
            lhs(14,5)=theta*(DN(3,0)*clhs66 + DN(3,1)*clhs147 + DN(3,2)*clhs194 + clhs274);
            lhs(14,6)=clhs153*clhs342 + clhs196*clhs343 + clhs288 + clhs341*clhs72 + clhs358;
            lhs(14,7)=DN(1,2)*clhs348 + clhs289*clhs76 - clhs290;
            lhs(14,8)=theta*(DN(3,0)*clhs81 + DN(3,1)*clhs158 + DN(3,2)*clhs201 + clhs315);
            lhs(14,9)=theta*(DN(3,0)*clhs94 + DN(3,1)*clhs162 + DN(3,2)*clhs203 + clhs327);
            lhs(14,10)=clhs100*clhs341 + clhs168*clhs342 + clhs205*clhs343 + clhs335 + clhs359;
            lhs(14,11)=DN(2,2)*clhs348 + clhs336*clhs76 - clhs337;
            lhs(14,12)=theta*(DN(3,0)*clhs108 + DN(3,1)*clhs173 + DN(3,2)*clhs210 + clhs355);
            lhs(14,13)=theta*(DN(3,0)*clhs121 + DN(3,1)*clhs177 + DN(3,2)*clhs212 + clhs361);
            lhs(14,14)=clhs127*clhs341 + clhs15*clhs362 + clhs183*clhs342 + clhs214*clhs343 + clhs352;
            lhs(14,15)=DN(3,2)*clhs356;
            lhs(15,0)=DN(3,0)*clhs220 + clhs129*theta;
            lhs(15,1)=DN(3,1)*clhs220 + clhs185*theta;
            lhs(15,2)=DN(3,2)*clhs220 + clhs217*theta;
            lhs(15,3)=clhs227;
            lhs(15,4)=DN(3,0)*clhs222 + clhs257*theta;
            lhs(15,5)=DN(3,1)*clhs222 + clhs275*theta;
            lhs(15,6)=DN(3,2)*clhs222 + clhs289*theta;
            lhs(15,7)=clhs294;
            lhs(15,8)=DN(3,0)*clhs224 + clhs316*theta;
            lhs(15,9)=DN(3,1)*clhs224 + clhs328*theta;
            lhs(15,10)=DN(3,2)*clhs224 + clhs336*theta;
            lhs(15,11)=clhs340;
            lhs(15,12)=DN(3,0)*clhs363;
            lhs(15,13)=DN(3,1)*clhs363;
            lhs(15,14)=DN(3,2)*clhs363;
            lhs(15,15)=clhs24*(clhs351 + clhs360 + clhs362);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0));
const double crhs2 =             1.0/dt;
const double crhs3 =             N[0]*rho;
const double crhs4 =             1.0*crhs2*crhs3;
const double crhs5 =             f(0,0)*theta;
const double crhs6 =             theta - 1.0;
const double crhs7 =             f(1,0)*theta;
const double crhs8 =             f(2,0)*theta;
const double crhs9 =             N[0]*(crhs5 - crhs6*fn(0,0)) + N[1]*(-crhs6*fn(1,0) + crhs7) + N[2]*(-crhs6*fn(2,0) + crhs8);
const double crhs10 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs11 =             theta*v(0,0);
const double crhs12 =             crhs11 - crhs6*vn(0,0);
const double crhs13 =             theta*v(1,0);
const double crhs14 =             crhs13 - crhs6*vn(1,0);
const double crhs15 =             theta*v(2,0);
const double crhs16 =             crhs15 - crhs6*vn(2,0);
const double crhs17 =             DN(0,0)*crhs12 + DN(1,0)*crhs14 + DN(2,0)*crhs16;
const double crhs18 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs19 =             crhs10*crhs17 + crhs18*(DN(0,1)*crhs12 + DN(1,1)*crhs14 + DN(2,1)*crhs16);
const double crhs20 =             theta*v(0,1);
const double crhs21 =             crhs20 - crhs6*vn(0,1);
const double crhs22 =             theta*v(1,1);
const double crhs23 =             crhs22 - crhs6*vn(1,1);
const double crhs24 =             theta*v(2,1);
const double crhs25 =             crhs24 - crhs6*vn(2,1);
const double crhs26 =             DN(0,1)*crhs21 + DN(1,1)*crhs23 + DN(2,1)*crhs25;
const double crhs27 =             crhs17 + crhs26 - volume_error_time_ratio;
const double crhs28 =             mu*stab_c1;
const double crhs29 =             rho*sqrt(pow(crhs10, 2) + pow(crhs18, 2));
const double crhs30 =             crhs27*(crhs28/stab_c2 + crhs29*h);
const double crhs31 =             crhs2*rho;
const double crhs32 =             1.0*crhs31;
const double crhs33 =             crhs1*crhs32;
const double crhs34 =             1.0 - theta;
const double crhs35 =             crhs11 + crhs34*vn(0,0);
const double crhs36 =             crhs13 + crhs34*vn(1,0);
const double crhs37 =             crhs15 + crhs34*vn(2,0);
const double crhs38 =             1.0/(crhs28/pow(h, 2) + crhs29*stab_c2/h + crhs31*dyn_tau);
const double crhs39 =             crhs38*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs33 + rho*(crhs10*(DN(0,0)*crhs35 + DN(1,0)*crhs36 + DN(2,0)*crhs37) + crhs18*(DN(0,1)*crhs35 + DN(1,1)*crhs36 + DN(2,1)*crhs37)) - rho*(N[0]*(crhs34*fn(0,0) + crhs5) + N[1]*(crhs34*fn(1,0) + crhs7) + N[2]*(crhs34*fn(2,0) + crhs8)));
const double crhs40 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs41 =             crhs3*crhs40;
const double crhs42 =             rho*(DN(0,0)*crhs10 + DN(0,1)*crhs18);
const double crhs43 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1));
const double crhs44 =             f(0,1)*theta;
const double crhs45 =             f(1,1)*theta;
const double crhs46 =             f(2,1)*theta;
const double crhs47 =             N[0]*(crhs44 - crhs6*fn(0,1)) + N[1]*(crhs45 - crhs6*fn(1,1)) + N[2]*(crhs46 - crhs6*fn(2,1));
const double crhs48 =             crhs10*(DN(0,0)*crhs21 + DN(1,0)*crhs23 + DN(2,0)*crhs25) + crhs18*crhs26;
const double crhs49 =             crhs32*crhs43;
const double crhs50 =             crhs20 + crhs34*vn(0,1);
const double crhs51 =             crhs22 + crhs34*vn(1,1);
const double crhs52 =             crhs24 + crhs34*vn(2,1);
const double crhs53 =             crhs38*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs49 + rho*(crhs10*(DN(0,0)*crhs50 + DN(1,0)*crhs51 + DN(2,0)*crhs52) + crhs18*(DN(0,1)*crhs50 + DN(1,1)*crhs51 + DN(2,1)*crhs52)) - rho*(N[0]*(crhs34*fn(0,1) + crhs44) + N[1]*(crhs34*fn(1,1) + crhs45) + N[2]*(crhs34*fn(2,1) + crhs46)));
const double crhs54 =             N[1]*rho;
const double crhs55 =             crhs40*crhs54;
const double crhs56 =             rho*(DN(1,0)*crhs10 + DN(1,1)*crhs18);
const double crhs57 =             N[2]*rho;
const double crhs58 =             crhs40*crhs57;
const double crhs59 =             rho*(DN(2,0)*crhs10 + DN(2,1)*crhs18);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs30 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - crhs1*crhs4 - crhs19*crhs3 + crhs3*crhs9 - crhs39*crhs41 - crhs39*crhs42;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs30 - DN(0,1)*stress[1] + crhs3*crhs47 - crhs3*crhs48 - crhs4*crhs43 - crhs41*crhs53 - crhs42*crhs53;
            rhs[2]=-DN(0,0)*crhs39 - DN(0,1)*crhs53 - N[0]*crhs27;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs30 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs33 - crhs19*crhs54 - crhs39*crhs55 - crhs39*crhs56 + crhs54*crhs9;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs30 - DN(1,1)*stress[1] - N[1]*crhs49 + crhs47*crhs54 - crhs48*crhs54 - crhs53*crhs55 - crhs53*crhs56;
            rhs[5]=-DN(1,0)*crhs39 - DN(1,1)*crhs53 - N[1]*crhs27;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs30 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs33 - crhs19*crhs57 - crhs39*crhs58 - crhs39*crhs59 + crhs57*crhs9;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs30 - DN(2,1)*stress[1] - N[2]*crhs49 + crhs47*crhs57 - crhs48*crhs57 - crhs53*crhs58 - crhs53*crhs59;
            rhs[8]=-DN(2,0)*crhs39 - DN(2,1)*crhs53 - N[2]*crhs27;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)) + N[3]*(v(3,0) - vn(3,0));
const double crhs2 =             1.0/dt;
const double crhs3 =             N[0]*rho;
const double crhs4 =             1.0*crhs2*crhs3;
const double crhs5 =             f(0,0)*theta;
const double crhs6 =             theta - 1.0;
const double crhs7 =             f(1,0)*theta;
const double crhs8 =             f(2,0)*theta;
const double crhs9 =             f(3,0)*theta;
const double crhs10 =             N[0]*(crhs5 - crhs6*fn(0,0)) + N[1]*(-crhs6*fn(1,0) + crhs7) + N[2]*(-crhs6*fn(2,0) + crhs8) + N[3]*(-crhs6*fn(3,0) + crhs9);
const double crhs11 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs12 =             theta*v(0,0);
const double crhs13 =             crhs12 - crhs6*vn(0,0);
const double crhs14 =             theta*v(1,0);
const double crhs15 =             crhs14 - crhs6*vn(1,0);
const double crhs16 =             theta*v(2,0);
const double crhs17 =             crhs16 - crhs6*vn(2,0);
const double crhs18 =             theta*v(3,0);
const double crhs19 =             crhs18 - crhs6*vn(3,0);
const double crhs20 =             DN(0,0)*crhs13 + DN(1,0)*crhs15 + DN(2,0)*crhs17 + DN(3,0)*crhs19;
const double crhs21 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs22 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs23 =             crhs11*crhs20 + crhs21*(DN(0,1)*crhs13 + DN(1,1)*crhs15 + DN(2,1)*crhs17 + DN(3,1)*crhs19) + crhs22*(DN(0,2)*crhs13 + DN(1,2)*crhs15 + DN(2,2)*crhs17 + DN(3,2)*crhs19);
const double crhs24 =             theta*v(0,1);
const double crhs25 =             crhs24 - crhs6*vn(0,1);
const double crhs26 =             theta*v(1,1);
const double crhs27 =             crhs26 - crhs6*vn(1,1);
const double crhs28 =             theta*v(2,1);
const double crhs29 =             crhs28 - crhs6*vn(2,1);
const double crhs30 =             theta*v(3,1);
const double crhs31 =             crhs30 - crhs6*vn(3,1);
const double crhs32 =             DN(0,1)*crhs25 + DN(1,1)*crhs27 + DN(2,1)*crhs29 + DN(3,1)*crhs31;
const double crhs33 =             theta*v(0,2);
const double crhs34 =             crhs33 - crhs6*vn(0,2);
const double crhs35 =             theta*v(1,2);
const double crhs36 =             crhs35 - crhs6*vn(1,2);
const double crhs37 =             theta*v(2,2);
const double crhs38 =             crhs37 - crhs6*vn(2,2);
const double crhs39 =             theta*v(3,2);
const double crhs40 =             crhs39 - crhs6*vn(3,2);
const double crhs41 =             DN(0,2)*crhs34 + DN(1,2)*crhs36 + DN(2,2)*crhs38 + DN(3,2)*crhs40;
const double crhs42 =             crhs20 + crhs32 + crhs41 - volume_error_time_ratio;
const double crhs43 =             mu*stab_c1;
const double crhs44 =             rho*sqrt(pow(crhs11, 2) + pow(crhs21, 2) + pow(crhs22, 2));
const double crhs45 =             crhs42*(crhs43/stab_c2 + crhs44*h);
const double crhs46 =             crhs2*rho;
const double crhs47 =             1.0*crhs46;
const double crhs48 =             crhs1*crhs47;
const double crhs49 =             1.0 - theta;
const double crhs50 =             crhs12 + crhs49*vn(0,0);
const double crhs51 =             crhs14 + crhs49*vn(1,0);
const double crhs52 =             crhs16 + crhs49*vn(2,0);
const double crhs53 =             crhs18 + crhs49*vn(3,0);
const double crhs54 =             1.0/(crhs43/pow(h, 2) + crhs44*stab_c2/h + crhs46*dyn_tau);
const double crhs55 =             crhs54*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs48 + rho*(crhs11*(DN(0,0)*crhs50 + DN(1,0)*crhs51 + DN(2,0)*crhs52 + DN(3,0)*crhs53) + crhs21*(DN(0,1)*crhs50 + DN(1,1)*crhs51 + DN(2,1)*crhs52 + DN(3,1)*crhs53) + crhs22*(DN(0,2)*crhs50 + DN(1,2)*crhs51 + DN(2,2)*crhs52 + DN(3,2)*crhs53)) - rho*(N[0]*(crhs49*fn(0,0) + crhs5) + N[1]*(crhs49*fn(1,0) + crhs7) + N[2]*(crhs49*fn(2,0) + crhs8) + N[3]*(crhs49*fn(3,0) + crhs9)));
const double crhs56 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs57 =             crhs3*crhs56;
const double crhs58 =             rho*(DN(0,0)*crhs11 + DN(0,1)*crhs21 + DN(0,2)*crhs22);
const double crhs59 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)) + N[3]*(v(3,1) - vn(3,1));
const double crhs60 =             f(0,1)*theta;
const double crhs61 =             f(1,1)*theta;
const double crhs62 =             f(2,1)*theta;
const double crhs63 =             f(3,1)*theta;
const double crhs64 =             N[0]*(-crhs6*fn(0,1) + crhs60) + N[1]*(-crhs6*fn(1,1) + crhs61) + N[2]*(-crhs6*fn(2,1) + crhs62) + N[3]*(-crhs6*fn(3,1) + crhs63);
const double crhs65 =             crhs11*(DN(0,0)*crhs25 + DN(1,0)*crhs27 + DN(2,0)*crhs29 + DN(3,0)*crhs31) + crhs21*crhs32 + crhs22*(DN(0,2)*crhs25 + DN(1,2)*crhs27 + DN(2,2)*crhs29 + DN(3,2)*crhs31);
const double crhs66 =             crhs47*crhs59;
const double crhs67 =             crhs24 + crhs49*vn(0,1);
const double crhs68 =             crhs26 + crhs49*vn(1,1);
const double crhs69 =             crhs28 + crhs49*vn(2,1);
const double crhs70 =             crhs30 + crhs49*vn(3,1);
const double crhs71 =             crhs54*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs66 + rho*(crhs11*(DN(0,0)*crhs67 + DN(1,0)*crhs68 + DN(2,0)*crhs69 + DN(3,0)*crhs70) + crhs21*(DN(0,1)*crhs67 + DN(1,1)*crhs68 + DN(2,1)*crhs69 + DN(3,1)*crhs70) + crhs22*(DN(0,2)*crhs67 + DN(1,2)*crhs68 + DN(2,2)*crhs69 + DN(3,2)*crhs70)) - rho*(N[0]*(crhs49*fn(0,1) + crhs60) + N[1]*(crhs49*fn(1,1) + crhs61) + N[2]*(crhs49*fn(2,1) + crhs62) + N[3]*(crhs49*fn(3,1) + crhs63)));
const double crhs72 =             N[0]*(v(0,2) - vn(0,2)) + N[1]*(v(1,2) - vn(1,2)) + N[2]*(v(2,2) - vn(2,2)) + N[3]*(v(3,2) - vn(3,2));
const double crhs73 =             f(0,2)*theta;
const double crhs74 =             f(1,2)*theta;
const double crhs75 =             f(2,2)*theta;
const double crhs76 =             f(3,2)*theta;
const double crhs77 =             N[0]*(-crhs6*fn(0,2) + crhs73) + N[1]*(-crhs6*fn(1,2) + crhs74) + N[2]*(-crhs6*fn(2,2) + crhs75) + N[3]*(-crhs6*fn(3,2) + crhs76);
const double crhs78 =             crhs11*(DN(0,0)*crhs34 + DN(1,0)*crhs36 + DN(2,0)*crhs38 + DN(3,0)*crhs40) + crhs21*(DN(0,1)*crhs34 + DN(1,1)*crhs36 + DN(2,1)*crhs38 + DN(3,1)*crhs40) + crhs22*crhs41;
const double crhs79 =             crhs47*crhs72;
const double crhs80 =             crhs33 + crhs49*vn(0,2);
const double crhs81 =             crhs35 + crhs49*vn(1,2);
const double crhs82 =             crhs37 + crhs49*vn(2,2);
const double crhs83 =             crhs39 + crhs49*vn(3,2);
const double crhs84 =             crhs54*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + crhs79 + rho*(crhs11*(DN(0,0)*crhs80 + DN(1,0)*crhs81 + DN(2,0)*crhs82 + DN(3,0)*crhs83) + crhs21*(DN(0,1)*crhs80 + DN(1,1)*crhs81 + DN(2,1)*crhs82 + DN(3,1)*crhs83) + crhs22*(DN(0,2)*crhs80 + DN(1,2)*crhs81 + DN(2,2)*crhs82 + DN(3,2)*crhs83)) - rho*(N[0]*(crhs49*fn(0,2) + crhs73) + N[1]*(crhs49*fn(1,2) + crhs74) + N[2]*(crhs49*fn(2,2) + crhs75) + N[3]*(crhs49*fn(3,2) + crhs76)));
const double crhs85 =             N[1]*rho;
const double crhs86 =             crhs56*crhs85;
const double crhs87 =             rho*(DN(1,0)*crhs11 + DN(1,1)*crhs21 + DN(1,2)*crhs22);
const double crhs88 =             N[2]*rho;
const double crhs89 =             crhs56*crhs88;
const double crhs90 =             rho*(DN(2,0)*crhs11 + DN(2,1)*crhs21 + DN(2,2)*crhs22);
const double crhs91 =             N[3]*rho;
const double crhs92 =             crhs56*crhs91;
const double crhs93 =             rho*(DN(3,0)*crhs11 + DN(3,1)*crhs21 + DN(3,2)*crhs22);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs45 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - crhs1*crhs4 + crhs10*crhs3 - crhs23*crhs3 - crhs55*crhs57 - crhs55*crhs58;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs45 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + crhs3*crhs64 - crhs3*crhs65 - crhs4*crhs59 - crhs57*crhs71 - crhs58*crhs71;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs45 - DN(0,2)*stress[2] + crhs3*crhs77 - crhs3*crhs78 - crhs4*crhs72 - crhs57*crhs84 - crhs58*crhs84;
            rhs[3]=-DN(0,0)*crhs55 - DN(0,1)*crhs71 - DN(0,2)*crhs84 - N[0]*crhs42;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs45 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs48 + crhs10*crhs85 - crhs23*crhs85 - crhs55*crhs86 - crhs55*crhs87;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs45 - DN(1,1)*stress[1] - DN(1,2)*stress[4] - N[1]*crhs66 + crhs64*crhs85 - crhs65*crhs85 - crhs71*crhs86 - crhs71*crhs87;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs45 - DN(1,2)*stress[2] - N[1]*crhs79 + crhs77*crhs85 - crhs78*crhs85 - crhs84*crhs86 - crhs84*crhs87;
            rhs[7]=-DN(1,0)*crhs55 - DN(1,1)*crhs71 - DN(1,2)*crhs84 - N[1]*crhs42;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs45 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs48 + crhs10*crhs88 - crhs23*crhs88 - crhs55*crhs89 - crhs55*crhs90;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs45 - DN(2,1)*stress[1] - DN(2,2)*stress[4] - N[2]*crhs66 + crhs64*crhs88 - crhs65*crhs88 - crhs71*crhs89 - crhs71*crhs90;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs45 - DN(2,2)*stress[2] - N[2]*crhs79 + crhs77*crhs88 - crhs78*crhs88 - crhs84*crhs89 - crhs84*crhs90;
            rhs[11]=-DN(2,0)*crhs55 - DN(2,1)*crhs71 - DN(2,2)*crhs84 - N[2]*crhs42;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs45 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs48 + crhs10*crhs91 - crhs23*crhs91 - crhs55*crhs92 - crhs55*crhs93;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs45 - DN(3,1)*stress[1] - DN(3,2)*stress[4] - N[3]*crhs66 + crhs64*crhs91 - crhs65*crhs91 - crhs71*crhs92 - crhs71*crhs93;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs45 - DN(3,2)*stress[2] - N[3]*crhs79 + crhs77*crhs91 - crhs78*crhs91 - crhs84*crhs92 - crhs84*crhs93;
            rhs[15]=-DN(3,0)*crhs55 - DN(3,1)*crhs71 - DN(3,2)*crhs84 - N[3]*crhs42;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p=rData.Pressure;
    const auto &pn=rData.Pressure_OldStep1;

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 =             cV2*rho;
const double cV4 =             cV3*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV5 =             cV3*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV6 =             cV3*(DN(2,0)*cV0 + DN(2,1)*cV1);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV4;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV4;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV4;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV4;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV4;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV4;
            V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV5;
            V(3,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV5;
            V(3,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV5;
            V(4,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV5;
            V(4,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV5;
            V(4,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV5;
            V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV6;
            V(6,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV6;
            V(6,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV6;
            V(7,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV6;
            V(7,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV6;
            V(7,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV6;
            V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             DN(0,0)*theta;
const double cH1 =             1.0/dt;
const double cH2 =             1.0*cH1;
const double cH3 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH4 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH5 =             DN(0,1)*theta;
const double cH6 =             1.0/(cH1*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 =             cH6*rho;
const double cH8 =             cH7*(N[0]*cH2 + cH0*cH3 + cH4*cH5);
const double cH9 =             DN(1,0)*theta;
const double cH10 =             DN(1,1)*theta;
const double cH11 =             cH7*(N[1]*cH2 + cH10*cH4 + cH3*cH9);
const double cH12 =             DN(2,0)*theta;
const double cH13 =             DN(2,1)*theta;
const double cH14 =             cH7*(N[2]*cH2 + cH12*cH3 + cH13*cH4);
            H(0,0)=DNenr(0,0)*cH8 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH8 + Nenr[0]*cH5;
            H(0,2)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DNenr(0,0)*cH11 + Nenr[0]*cH9;
            H(0,4)=DNenr(0,1)*cH11 + Nenr[0]*cH10;
            H(0,5)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DNenr(0,0)*cH14 + Nenr[0]*cH12;
            H(0,7)=DNenr(0,1)*cH14 + Nenr[0]*cH13;
            H(0,8)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DNenr(1,0)*cH8 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH8 + Nenr[1]*cH5;
            H(1,2)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DNenr(1,0)*cH11 + Nenr[1]*cH9;
            H(1,4)=DNenr(1,1)*cH11 + Nenr[1]*cH10;
            H(1,5)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DNenr(1,0)*cH14 + Nenr[1]*cH12;
            H(1,7)=DNenr(1,1)*cH14 + Nenr[1]*cH13;
            H(1,8)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DNenr(2,0)*cH8 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH8 + Nenr[2]*cH5;
            H(2,2)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DNenr(2,0)*cH11 + Nenr[2]*cH9;
            H(2,4)=DNenr(2,1)*cH11 + Nenr[2]*cH10;
            H(2,5)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DNenr(2,0)*cH14 + Nenr[2]*cH12;
            H(2,7)=DNenr(2,1)*cH14 + Nenr[2]*cH13;
            H(2,8)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
            Kee(1,2)=cKee3;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee3;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 =             theta*v(0,0);
const double crhs_ee1 =             theta - 1.0;
const double crhs_ee2 =             theta*v(0,1);
const double crhs_ee3 =             theta*v(1,0);
const double crhs_ee4 =             theta*v(1,1);
const double crhs_ee5 =             theta*v(2,0);
const double crhs_ee6 =             theta*v(2,1);
const double crhs_ee7 =             DN(0,0)*(crhs_ee0 - crhs_ee1*vn(0,0)) + DN(0,1)*(-crhs_ee1*vn(0,1) + crhs_ee2) + DN(1,0)*(-crhs_ee1*vn(1,0) + crhs_ee3) + DN(1,1)*(-crhs_ee1*vn(1,1) + crhs_ee4) + DN(2,0)*(-crhs_ee1*vn(2,0) + crhs_ee5) + DN(2,1)*(-crhs_ee1*vn(2,1) + crhs_ee6) - volume_error_time_ratio;
const double crhs_ee8 =             1.0 - theta;
const double crhs_ee9 =             1.0/dt;
const double crhs_ee10 =             1.0*crhs_ee9;
const double crhs_ee11 =             N[0]*crhs_ee10;
const double crhs_ee12 =             N[1]*crhs_ee10;
const double crhs_ee13 =             N[2]*crhs_ee10;
const double crhs_ee14 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee15 =             crhs_ee0 + crhs_ee8*vn(0,0);
const double crhs_ee16 =             crhs_ee3 + crhs_ee8*vn(1,0);
const double crhs_ee17 =             crhs_ee5 + crhs_ee8*vn(2,0);
const double crhs_ee18 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee19 =             1.0/(crhs_ee9*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee14, 2) + pow(crhs_ee18, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee20 =             crhs_ee19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee8*fn(0,0) + f(0,0)*theta) + N[1]*(crhs_ee8*fn(1,0) + f(1,0)*theta) + N[2]*(crhs_ee8*fn(2,0) + f(2,0)*theta)) + rho*(crhs_ee11*(v(0,0) - vn(0,0)) + crhs_ee12*(v(1,0) - vn(1,0)) + crhs_ee13*(v(2,0) - vn(2,0)) + crhs_ee14*(DN(0,0)*crhs_ee15 + DN(1,0)*crhs_ee16 + DN(2,0)*crhs_ee17) + crhs_ee18*(DN(0,1)*crhs_ee15 + DN(1,1)*crhs_ee16 + DN(2,1)*crhs_ee17)));
const double crhs_ee21 =             crhs_ee2 + crhs_ee8*vn(0,1);
const double crhs_ee22 =             crhs_ee4 + crhs_ee8*vn(1,1);
const double crhs_ee23 =             crhs_ee6 + crhs_ee8*vn(2,1);
const double crhs_ee24 =             crhs_ee19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee8*fn(0,1) + f(0,1)*theta) + N[1]*(crhs_ee8*fn(1,1) + f(1,1)*theta) + N[2]*(crhs_ee8*fn(2,1) + f(2,1)*theta)) + rho*(crhs_ee11*(v(0,1) - vn(0,1)) + crhs_ee12*(v(1,1) - vn(1,1)) + crhs_ee13*(v(2,1) - vn(2,1)) + crhs_ee14*(DN(0,0)*crhs_ee21 + DN(1,0)*crhs_ee22 + DN(2,0)*crhs_ee23) + crhs_ee18*(DN(0,1)*crhs_ee21 + DN(1,1)*crhs_ee22 + DN(2,1)*crhs_ee23)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee20 - DNenr(0,1)*crhs_ee24 - Nenr[0]*crhs_ee7;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee20 - DNenr(1,1)*crhs_ee24 - Nenr[1]*crhs_ee7;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee20 - DNenr(2,1)*crhs_ee24 - Nenr[2]*crhs_ee7;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 =             cV3*rho;
const double cV5 =             cV4*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV6 =             cV4*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV7 =             cV4*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV8 =             cV4*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV5;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV5;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV5;
            V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV5;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV5;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV5;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV5;
            V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV5;
            V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV5;
            V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV5;
            V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV5;
            V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV5;
            V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV6;
            V(4,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV6;
            V(4,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV6;
            V(4,3)=-DN(1,0)*Nenr[3] + DNenr(3,0)*cV6;
            V(5,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV6;
            V(5,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV6;
            V(5,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV6;
            V(5,3)=-DN(1,1)*Nenr[3] + DNenr(3,1)*cV6;
            V(6,0)=-DN(1,2)*Nenr[0] + DNenr(0,2)*cV6;
            V(6,1)=-DN(1,2)*Nenr[1] + DNenr(1,2)*cV6;
            V(6,2)=-DN(1,2)*Nenr[2] + DNenr(2,2)*cV6;
            V(6,3)=-DN(1,2)*Nenr[3] + DNenr(3,2)*cV6;
            V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV7;
            V(8,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV7;
            V(8,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV7;
            V(8,3)=-DN(2,0)*Nenr[3] + DNenr(3,0)*cV7;
            V(9,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV7;
            V(9,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV7;
            V(9,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV7;
            V(9,3)=-DN(2,1)*Nenr[3] + DNenr(3,1)*cV7;
            V(10,0)=-DN(2,2)*Nenr[0] + DNenr(0,2)*cV7;
            V(10,1)=-DN(2,2)*Nenr[1] + DNenr(1,2)*cV7;
            V(10,2)=-DN(2,2)*Nenr[2] + DNenr(2,2)*cV7;
            V(10,3)=-DN(2,2)*Nenr[3] + DNenr(3,2)*cV7;
            V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] + DNenr(0,0)*cV8;
            V(12,1)=-DN(3,0)*Nenr[1] + DNenr(1,0)*cV8;
            V(12,2)=-DN(3,0)*Nenr[2] + DNenr(2,0)*cV8;
            V(12,3)=-DN(3,0)*Nenr[3] + DNenr(3,0)*cV8;
            V(13,0)=-DN(3,1)*Nenr[0] + DNenr(0,1)*cV8;
            V(13,1)=-DN(3,1)*Nenr[1] + DNenr(1,1)*cV8;
            V(13,2)=-DN(3,1)*Nenr[2] + DNenr(2,1)*cV8;
            V(13,3)=-DN(3,1)*Nenr[3] + DNenr(3,1)*cV8;
            V(14,0)=-DN(3,2)*Nenr[0] + DNenr(0,2)*cV8;
            V(14,1)=-DN(3,2)*Nenr[1] + DNenr(1,2)*cV8;
            V(14,2)=-DN(3,2)*Nenr[2] + DNenr(2,2)*cV8;
            V(14,3)=-DN(3,2)*Nenr[3] + DNenr(3,2)*cV8;
            V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             DN(0,0)*theta;
const double cH1 =             1.0/dt;
const double cH2 =             1.0*cH1;
const double cH3 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH4 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH5 =             DN(0,1)*theta;
const double cH6 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH7 =             DN(0,2)*theta;
const double cH8 =             1.0/(cH1*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2) + pow(cH6, 2))/h + mu*stab_c1/pow(h, 2));
const double cH9 =             cH8*rho;
const double cH10 =             cH9*(N[0]*cH2 + cH0*cH3 + cH4*cH5 + cH6*cH7);
const double cH11 =             DN(1,0)*theta;
const double cH12 =             DN(1,1)*theta;
const double cH13 =             DN(1,2)*theta;
const double cH14 =             cH9*(N[1]*cH2 + cH11*cH3 + cH12*cH4 + cH13*cH6);
const double cH15 =             DN(2,0)*theta;
const double cH16 =             DN(2,1)*theta;
const double cH17 =             DN(2,2)*theta;
const double cH18 =             cH9*(N[2]*cH2 + cH15*cH3 + cH16*cH4 + cH17*cH6);
const double cH19 =             DN(3,0)*theta;
const double cH20 =             DN(3,1)*theta;
const double cH21 =             DN(3,2)*theta;
const double cH22 =             cH9*(N[3]*cH2 + cH19*cH3 + cH20*cH4 + cH21*cH6);
            H(0,0)=DNenr(0,0)*cH10 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH10 + Nenr[0]*cH5;
            H(0,2)=DNenr(0,2)*cH10 + Nenr[0]*cH7;
            H(0,3)=cH8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DNenr(0,0)*cH14 + Nenr[0]*cH11;
            H(0,5)=DNenr(0,1)*cH14 + Nenr[0]*cH12;
            H(0,6)=DNenr(0,2)*cH14 + Nenr[0]*cH13;
            H(0,7)=cH8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DNenr(0,0)*cH18 + Nenr[0]*cH15;
            H(0,9)=DNenr(0,1)*cH18 + Nenr[0]*cH16;
            H(0,10)=DNenr(0,2)*cH18 + Nenr[0]*cH17;
            H(0,11)=cH8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DNenr(0,0)*cH22 + Nenr[0]*cH19;
            H(0,13)=DNenr(0,1)*cH22 + Nenr[0]*cH20;
            H(0,14)=DNenr(0,2)*cH22 + Nenr[0]*cH21;
            H(0,15)=cH8*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DNenr(1,0)*cH10 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH10 + Nenr[1]*cH5;
            H(1,2)=DNenr(1,2)*cH10 + Nenr[1]*cH7;
            H(1,3)=cH8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DNenr(1,0)*cH14 + Nenr[1]*cH11;
            H(1,5)=DNenr(1,1)*cH14 + Nenr[1]*cH12;
            H(1,6)=DNenr(1,2)*cH14 + Nenr[1]*cH13;
            H(1,7)=cH8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DNenr(1,0)*cH18 + Nenr[1]*cH15;
            H(1,9)=DNenr(1,1)*cH18 + Nenr[1]*cH16;
            H(1,10)=DNenr(1,2)*cH18 + Nenr[1]*cH17;
            H(1,11)=cH8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DNenr(1,0)*cH22 + Nenr[1]*cH19;
            H(1,13)=DNenr(1,1)*cH22 + Nenr[1]*cH20;
            H(1,14)=DNenr(1,2)*cH22 + Nenr[1]*cH21;
            H(1,15)=cH8*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DNenr(2,0)*cH10 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH10 + Nenr[2]*cH5;
            H(2,2)=DNenr(2,2)*cH10 + Nenr[2]*cH7;
            H(2,3)=cH8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DNenr(2,0)*cH14 + Nenr[2]*cH11;
            H(2,5)=DNenr(2,1)*cH14 + Nenr[2]*cH12;
            H(2,6)=DNenr(2,2)*cH14 + Nenr[2]*cH13;
            H(2,7)=cH8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DNenr(2,0)*cH18 + Nenr[2]*cH15;
            H(2,9)=DNenr(2,1)*cH18 + Nenr[2]*cH16;
            H(2,10)=DNenr(2,2)*cH18 + Nenr[2]*cH17;
            H(2,11)=cH8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DNenr(2,0)*cH22 + Nenr[2]*cH19;
            H(2,13)=DNenr(2,1)*cH22 + Nenr[2]*cH20;
            H(2,14)=DNenr(2,2)*cH22 + Nenr[2]*cH21;
            H(2,15)=cH8*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DNenr(3,0)*cH10 + Nenr[3]*cH0;
            H(3,1)=DNenr(3,1)*cH10 + Nenr[3]*cH5;
            H(3,2)=DNenr(3,2)*cH10 + Nenr[3]*cH7;
            H(3,3)=cH8*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DNenr(3,0)*cH14 + Nenr[3]*cH11;
            H(3,5)=DNenr(3,1)*cH14 + Nenr[3]*cH12;
            H(3,6)=DNenr(3,2)*cH14 + Nenr[3]*cH13;
            H(3,7)=cH8*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DNenr(3,0)*cH18 + Nenr[3]*cH15;
            H(3,9)=DNenr(3,1)*cH18 + Nenr[3]*cH16;
            H(3,10)=DNenr(3,2)*cH18 + Nenr[3]*cH17;
            H(3,11)=cH8*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DNenr(3,0)*cH22 + Nenr[3]*cH19;
            H(3,13)=DNenr(3,1)*cH22 + Nenr[3]*cH20;
            H(3,14)=DNenr(3,2)*cH22 + Nenr[3]*cH21;
            H(3,15)=cH8*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 =             cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 =             cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 =             cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2) + pow(DNenr(0,2), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(0,3)=cKee3;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2) + pow(DNenr(1,2), 2));
            Kee(1,2)=cKee4;
            Kee(1,3)=cKee5;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee4;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2) + pow(DNenr(2,2), 2));
            Kee(2,3)=cKee6;
            Kee(3,0)=cKee3;
            Kee(3,1)=cKee5;
            Kee(3,2)=cKee6;
            Kee(3,3)=cKee0*(pow(DNenr(3,0), 2) + pow(DNenr(3,1), 2) + pow(DNenr(3,2), 2));


    const double crhs_ee0 =             theta*v(0,0);
const double crhs_ee1 =             theta - 1.0;
const double crhs_ee2 =             theta*v(0,1);
const double crhs_ee3 =             theta*v(0,2);
const double crhs_ee4 =             theta*v(1,0);
const double crhs_ee5 =             theta*v(1,1);
const double crhs_ee6 =             theta*v(1,2);
const double crhs_ee7 =             theta*v(2,0);
const double crhs_ee8 =             theta*v(2,1);
const double crhs_ee9 =             theta*v(2,2);
const double crhs_ee10 =             theta*v(3,0);
const double crhs_ee11 =             theta*v(3,1);
const double crhs_ee12 =             theta*v(3,2);
const double crhs_ee13 =             DN(0,0)*(crhs_ee0 - crhs_ee1*vn(0,0)) + DN(0,1)*(-crhs_ee1*vn(0,1) + crhs_ee2) + DN(0,2)*(-crhs_ee1*vn(0,2) + crhs_ee3) + DN(1,0)*(-crhs_ee1*vn(1,0) + crhs_ee4) + DN(1,1)*(-crhs_ee1*vn(1,1) + crhs_ee5) + DN(1,2)*(-crhs_ee1*vn(1,2) + crhs_ee6) + DN(2,0)*(-crhs_ee1*vn(2,0) + crhs_ee7) + DN(2,1)*(-crhs_ee1*vn(2,1) + crhs_ee8) + DN(2,2)*(-crhs_ee1*vn(2,2) + crhs_ee9) + DN(3,0)*(-crhs_ee1*vn(3,0) + crhs_ee10) + DN(3,1)*(-crhs_ee1*vn(3,1) + crhs_ee11) + DN(3,2)*(-crhs_ee1*vn(3,2) + crhs_ee12) - volume_error_time_ratio;
const double crhs_ee14 =             1.0 - theta;
const double crhs_ee15 =             1.0/dt;
const double crhs_ee16 =             1.0*crhs_ee15;
const double crhs_ee17 =             N[0]*crhs_ee16;
const double crhs_ee18 =             N[1]*crhs_ee16;
const double crhs_ee19 =             N[2]*crhs_ee16;
const double crhs_ee20 =             N[3]*crhs_ee16;
const double crhs_ee21 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee22 =             crhs_ee0 + crhs_ee14*vn(0,0);
const double crhs_ee23 =             crhs_ee14*vn(1,0) + crhs_ee4;
const double crhs_ee24 =             crhs_ee14*vn(2,0) + crhs_ee7;
const double crhs_ee25 =             crhs_ee10 + crhs_ee14*vn(3,0);
const double crhs_ee26 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee27 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee28 =             1.0/(crhs_ee15*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee21, 2) + pow(crhs_ee26, 2) + pow(crhs_ee27, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee29 =             crhs_ee28*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee14*fn(0,0) + f(0,0)*theta) + N[1]*(crhs_ee14*fn(1,0) + f(1,0)*theta) + N[2]*(crhs_ee14*fn(2,0) + f(2,0)*theta) + N[3]*(crhs_ee14*fn(3,0) + f(3,0)*theta)) + rho*(crhs_ee17*(v(0,0) - vn(0,0)) + crhs_ee18*(v(1,0) - vn(1,0)) + crhs_ee19*(v(2,0) - vn(2,0)) + crhs_ee20*(v(3,0) - vn(3,0)) + crhs_ee21*(DN(0,0)*crhs_ee22 + DN(1,0)*crhs_ee23 + DN(2,0)*crhs_ee24 + DN(3,0)*crhs_ee25) + crhs_ee26*(DN(0,1)*crhs_ee22 + DN(1,1)*crhs_ee23 + DN(2,1)*crhs_ee24 + DN(3,1)*crhs_ee25) + crhs_ee27*(DN(0,2)*crhs_ee22 + DN(1,2)*crhs_ee23 + DN(2,2)*crhs_ee24 + DN(3,2)*crhs_ee25)));
const double crhs_ee30 =             crhs_ee14*vn(0,1) + crhs_ee2;
const double crhs_ee31 =             crhs_ee14*vn(1,1) + crhs_ee5;
const double crhs_ee32 =             crhs_ee14*vn(2,1) + crhs_ee8;
const double crhs_ee33 =             crhs_ee11 + crhs_ee14*vn(3,1);
const double crhs_ee34 =             crhs_ee28*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee14*fn(0,1) + f(0,1)*theta) + N[1]*(crhs_ee14*fn(1,1) + f(1,1)*theta) + N[2]*(crhs_ee14*fn(2,1) + f(2,1)*theta) + N[3]*(crhs_ee14*fn(3,1) + f(3,1)*theta)) + rho*(crhs_ee17*(v(0,1) - vn(0,1)) + crhs_ee18*(v(1,1) - vn(1,1)) + crhs_ee19*(v(2,1) - vn(2,1)) + crhs_ee20*(v(3,1) - vn(3,1)) + crhs_ee21*(DN(0,0)*crhs_ee30 + DN(1,0)*crhs_ee31 + DN(2,0)*crhs_ee32 + DN(3,0)*crhs_ee33) + crhs_ee26*(DN(0,1)*crhs_ee30 + DN(1,1)*crhs_ee31 + DN(2,1)*crhs_ee32 + DN(3,1)*crhs_ee33) + crhs_ee27*(DN(0,2)*crhs_ee30 + DN(1,2)*crhs_ee31 + DN(2,2)*crhs_ee32 + DN(3,2)*crhs_ee33)));
const double crhs_ee35 =             crhs_ee14*vn(0,2) + crhs_ee3;
const double crhs_ee36 =             crhs_ee14*vn(1,2) + crhs_ee6;
const double crhs_ee37 =             crhs_ee14*vn(2,2) + crhs_ee9;
const double crhs_ee38 =             crhs_ee12 + crhs_ee14*vn(3,2);
const double crhs_ee39 =             crhs_ee28*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee14*fn(0,2) + f(0,2)*theta) + N[1]*(crhs_ee14*fn(1,2) + f(1,2)*theta) + N[2]*(crhs_ee14*fn(2,2) + f(2,2)*theta) + N[3]*(crhs_ee14*fn(3,2) + f(3,2)*theta)) + rho*(crhs_ee17*(v(0,2) - vn(0,2)) + crhs_ee18*(v(1,2) - vn(1,2)) + crhs_ee19*(v(2,2) - vn(2,2)) + crhs_ee20*(v(3,2) - vn(3,2)) + crhs_ee21*(DN(0,0)*crhs_ee35 + DN(1,0)*crhs_ee36 + DN(2,0)*crhs_ee37 + DN(3,0)*crhs_ee38) + crhs_ee26*(DN(0,1)*crhs_ee35 + DN(1,1)*crhs_ee36 + DN(2,1)*crhs_ee37 + DN(3,1)*crhs_ee38) + crhs_ee27*(DN(0,2)*crhs_ee35 + DN(1,2)*crhs_ee36 + DN(2,2)*crhs_ee37 + DN(3,2)*crhs_ee38)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee29 - DNenr(0,1)*crhs_ee34 - DNenr(0,2)*crhs_ee39 - Nenr[0]*crhs_ee13;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee29 - DNenr(1,1)*crhs_ee34 - DNenr(1,2)*crhs_ee39 - Nenr[1]*crhs_ee13;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee29 - DNenr(2,1)*crhs_ee34 - DNenr(2,2)*crhs_ee39 - Nenr[2]*crhs_ee13;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee29 - DNenr(3,1)*crhs_ee34 - DNenr(3,2)*crhs_ee39 - Nenr[3]*crhs_ee13;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::ComputeSplitting(
    TElementData &rData,
    MatrixType &rShapeFunctionsPos,
    MatrixType &rShapeFunctionsNeg,
    MatrixType &rEnrichedShapeFunctionsPos,
    MatrixType &rEnrichedShapeFunctionsNeg,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
{
    // Set the positive and negative enrichment interpolation matrices
    // Note that the enrichment is constructed using the standard shape functions such that:
    // In the negative distance region, the enrichment functions correspondig to the negative
    // distance nodes are null and the positive distance nodes are equal to the standard shape
    // functions. On the contrary, for the positive distance region, the enrichment functions
    // corresponding to the positive distance nodes are null meanwhile the negative distance
    // nodes are equal to the standard. This yields a discontinuous enrichment space.
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the positive side modified shape functions calculator
    pModifiedShapeFunctions->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); ++i){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); ++i){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (pModifiedShapeFunctions->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::ComputeSplitInterface(
    const TElementData &rData,
    MatrixType& rInterfaceShapeFunctionNeg,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Vector& rInterfaceWeightsNeg,
    std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
{
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the Interface negative side shape functions calculator
    pModifiedShapeFunctions->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivativesNeg,
        rInterfaceWeightsNeg,
        GeometryData::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideInterfaceAreaNormals(
        rInterfaceNormalsNeg,
        GeometryData::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); ++gp){
        const double normal_norm = norm_2(rInterfaceNormalsNeg[gp]);
        rInterfaceNormalsNeg[gp] /= normal_norm;
    }

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesCNData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
   return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesCNData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateCurvatureOnInterfaceGaussPoints(
        const Matrix& rInterfaceShapeFunctions,
        Vector& rInterfaceCurvature)
{
    const auto& r_geom = this->GetGeometry();
    const unsigned int n_gpt = rInterfaceShapeFunctions.size1();

    rInterfaceCurvature.resize(n_gpt, false);

    for (unsigned int gpt = 0; gpt < n_gpt; ++gpt){
        double curvature = 0.0;
        for (unsigned int i = 0; i < NumNodes; ++i){
            curvature += rInterfaceShapeFunctions(gpt,i) * r_geom[i].GetValue(CURVATURE);
        }
        rInterfaceCurvature[gpt] = curvature;
    }
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::SurfaceTension(
    const double SurfaceTensionCoefficient,
    const Vector& rCurvature,
    const Vector& rInterfaceWeights,
    const Matrix& rInterfaceShapeFunctions,
    const std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
    VectorType& rRHS)
{
    for (unsigned int intgp = 0; intgp < rInterfaceWeights.size(); ++intgp){
        const double intgp_curv = rCurvature(intgp);
        const double intgp_w = rInterfaceWeights(intgp);
        const auto& intgp_normal = rInterfaceNormalsNeg[intgp];
        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < NumNodes-1; ++dim){
                rRHS[ j*(NumNodes) + dim ] -= SurfaceTensionCoefficient*intgp_normal[dim]
                    *intgp_curv*intgp_w*rInterfaceShapeFunctions(intgp,j);
            }
        }
    }
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::PressureGradientStabilization(
    const TElementData& rData,
    const Vector& rInterfaceWeights,
    const Matrix& rEnrInterfaceShapeFunctionPos,
    const Matrix& rEnrInterfaceShapeFunctionNeg,
    const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivatives,
    MatrixType& rKeeTot,
    VectorType& rRHSeeTot)
{
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);
    VectorType rhs_enr = ZeroVector(NumNodes);

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    double positive_density = 0.0;
    double negative_density = 0.0;
    double positive_viscosity = 0.0;
    double negative_viscosity = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
            positive_density = rData.NodalDensity[i];
            positive_viscosity = rData.NodalDynamicViscosity[i];
        } else{
            enr_pos_interp(i, i) = 1.0;
            negative_density = rData.NodalDensity[i];
            negative_viscosity = rData.NodalDynamicViscosity[i];
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesPos = rInterfaceShapeDerivatives;
    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesNeg = rInterfaceShapeDerivatives;

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesPos[i] = prod(enr_pos_interp, rInterfaceShapeDerivatives[i]);
    }

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesNeg[i] = prod(enr_neg_interp, rInterfaceShapeDerivatives[i]);
    }

    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double element_volume = positive_volume + negative_volume;

    const auto& r_geom = this->GetGeometry();
    const double h_elem = rData.ElementSize;

    double cut_area = 0.0;
    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){
        cut_area += rInterfaceWeights[gp];
    }

    const double density = 1.0/(1.0/positive_density + 1.0/negative_density);
    const double viscosity = 1.0/(1.0/positive_viscosity + 1.0/negative_viscosity);

    // Stabilization parameters
    const double cut_stabilization_coefficient = 1.0;
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    const double dyn_tau = rData.DynamicTau;

    const double dt = rData.DeltaTime;
    const auto &v=rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;

    const double theta = rData.theta;
    const auto vmesh=rData.MeshVelocity;
    const auto vmeshn=rData.MeshVelocityOldStep;
    const BoundedMatrix<double,NumNodes,Dim> v_convection = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){

        Vector vconv = ZeroVector(Dim);
        double positive_weight = 0.0;
        double negative_weight = 0.0;

        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < Dim; ++dim){
                vconv[dim] += (rEnrInterfaceShapeFunctionNeg(gp, j) + rEnrInterfaceShapeFunctionPos(gp, j))
                    *v_convection(j,dim);
            }
            positive_weight += rEnrInterfaceShapeFunctionNeg(gp, j);
            negative_weight += rEnrInterfaceShapeFunctionPos(gp, j);
        }

        const double v_conv_norm = norm_2(vconv);

        const double penalty_coefficient = cut_stabilization_coefficient *
            density * 1.0 / (dyn_tau * density / (0.5*dt) + stab_c1 * viscosity / h_elem / h_elem +
                                stab_c2 * density * v_conv_norm / h_elem) * element_volume / cut_area;

        const auto& r_gp_enriched_interface_shape_derivatives_pos = EnrichedInterfaceShapeDerivativesPos[gp];
        const auto& r_gp_enriched_interface_shape_derivatives_neg = EnrichedInterfaceShapeDerivativesNeg[gp];

        for (unsigned int i = 0; i < NumNodes; ++i){

            for (unsigned int j = 0; j < NumNodes; ++j){

                const auto& r_pressure_gradient_j = r_geom[j].GetValue(PRESSURE_GRADIENT);

                for (unsigned int dim = 0; dim < Dim; ++dim){
                    kee(i, j) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        ( r_gp_enriched_interface_shape_derivatives_pos(j,dim) - r_gp_enriched_interface_shape_derivatives_neg(j,dim) );

                    rhs_enr(i) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        (rEnrInterfaceShapeFunctionNeg(gp, j)/positive_weight - rEnrInterfaceShapeFunctionPos(gp, j)/negative_weight)*
                        r_pressure_gradient_j(dim);
                }
            }
        }
    }

    noalias(rKeeTot) += kee;
    noalias(rRHSeeTot) += rhs_enr;
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CondenseEnrichmentWithContinuity(
    const TElementData &rData,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    const double min_area_ratio = 1e-7;

    // Compute positive side, negative side and total volumes
    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double Vol = positive_volume + negative_volume;

    //We only enrich elements which are not almost empty/full
    if (positive_volume / Vol > min_area_ratio && negative_volume / Vol > min_area_ratio) {

        // Compute the maximum diagonal value in the enrichment stiffness matrix
        double max_diag = 0.0;
        for (unsigned int k = 0; k < NumNodes; ++k){
            if (std::abs(rKeeTot(k, k)) > max_diag){
                max_diag = std::abs(rKeeTot(k, k));
            }
        }
        if (max_diag == 0.0){
            max_diag = 1.0;
        }
        // "weakly" impose continuity
        for (unsigned int i = 0; i < Dim; ++i){
            const double di = std::abs(rData.Distance[i]);
            for (unsigned int j = i + 1; j < NumNodes; ++j){
                const double dj = std::abs(rData.Distance[j]);
                // Check if the edge is cut, if it is, set the penalty constraint
                if (rData.Distance[i] * rData.Distance[j] < 0.0){
                    double sum_d = di + dj;
                    double Ni = dj / sum_d;
                    double Nj = di / sum_d;
                    double penalty_coeff = max_diag * 0.001;
                    rKeeTot(i, i) += penalty_coeff * Ni * Ni;
                    rKeeTot(i, j) -= penalty_coeff * Ni * Nj;
                    rKeeTot(j, i) -= penalty_coeff * Nj * Ni;
                    rKeeTot(j, j) += penalty_coeff * Nj * Nj;
                }
            }
        }

        // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        double det;
        MatrixType inverse_diag(NumNodes, NumNodes);
        MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

        const Matrix tmp = prod(inverse_diag, rHtot);
        noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

        const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
        noalias(rRightHandSideVector) -= prod(rVtot, tmp2);
    }
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CondenseEnrichment(
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    // Enrichment condensation (add to LHS and RHS the enrichment contributions)
    double det;
    MatrixType inverse_diag(NumNodes, NumNodes);
    MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

    const Matrix tmp = prod(inverse_diag, rHtot);
    noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

    const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
    noalias(rRightHandSideVector) -= prod(rVtot, tmp2);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddSurfaceTensionContribution(
    const TElementData& rData,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    VectorType &rRHSeeTot)
{
    // Surface tension coefficient is set in material properties
    const double surface_tension_coefficient = this->GetProperties().GetValue(SURFACE_TENSION_COEFFICIENT);
    Matrix int_shape_function, int_shape_function_enr_neg, int_shape_function_enr_pos;
    GeometryType::ShapeFunctionsGradientsType int_shape_derivatives;
    Vector int_gauss_pts_weights;
    std::vector<array_1d<double,3>> int_normals_neg;
    Vector gauss_pts_curvature;

    ComputeSplitInterface(
        rData,
        int_shape_function,
        int_shape_function_enr_pos,
        int_shape_function_enr_neg,
        int_shape_derivatives,
        int_gauss_pts_weights,
        int_normals_neg,
        pModifiedShapeFunctions);

    CalculateCurvatureOnInterfaceGaussPoints(
        int_shape_function,
        gauss_pts_curvature);

    SurfaceTension(
        surface_tension_coefficient,
        gauss_pts_curvature,
        int_gauss_pts_weights,
        int_shape_function,
        int_normals_neg,
        rRightHandSideVector);

    PressureGradientStabilization(
        rData,
        int_gauss_pts_weights,
        int_shape_function_enr_pos,
        int_shape_function_enr_neg,
        int_shape_derivatives,
        rKeeTot,
        rRHSeeTot);

    CondenseEnrichment(rLeftHandSideMatrix, rRightHandSideVector, rHtot, rVtot, rKeeTot, rRHSeeTot);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo )
{
    if (rVariable == DIVERGENCE){

        const auto& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int num_gauss = IntegrationPoints.size();

        if (rValues.size() != num_gauss){
            rValues.resize(num_gauss);
        }

        Vector gauss_pts_jacobian_determinant = ZeroVector(num_gauss);
        GeometryData::ShapeFunctionsGradientsType DN_DX;
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, gauss_pts_jacobian_determinant, GeometryData::GI_GAUSS_2);

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; ++i_gauss){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; ++nnode){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ++ndim){
                    DVi_DXi += gp_DN_DX(nnode, ndim) * vel[ndim];
                }
            }
            rValues[i_gauss] = DVi_DXi;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>;
template class TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>;

} // namespace Kratos
