//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main author:     Uxue Chasco

#include "two_fluid_navier_stokes_alpha_method.h"
#include "custom_utilities/two_fluid_navier_stokes_alpha_method_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(IndexType NewId)
    : TwoFluidNavierStokes<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : TwoFluidNavierStokes<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : TwoFluidNavierStokes<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : TwoFluidNavierStokes<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::~TwoFluidNavierStokesAlphaMethod() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethod<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethod>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethod<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethod>(NewId, pGeom, pProperties);
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodData<2, 3>& rData) const
{
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const BoundedMatrix<double,3,2> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += rDNDX(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += rDNDX(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += rDNDX(i,0)*velocity_alpha(i,1) + rDNDX(i,1)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodData<3, 4>& rData) const
{
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const BoundedMatrix<double,4,3> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(6);
    for (unsigned int i = 0; i < 4; i++) {
        r_strain_rate[0] += rDNDX(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += rDNDX(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += rDNDX(i,2)*velocity_alpha(i,2);
        r_strain_rate[3] += rDNDX(i,0)*velocity_alpha(i,1) + rDNDX(i,1)*velocity_alpha(i,0);
        r_strain_rate[4] += rDNDX(i,1)*velocity_alpha(i,2) + rDNDX(i,2)*velocity_alpha(i,1);
        r_strain_rate[5] += rDNDX(i,0)*velocity_alpha(i,2) + rDNDX(i,2)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const double dyn_tau = rData.DynamicTau;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const BoundedMatrix<double,3,2> vconv =(vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             1.0/(max_spectral_radius + 1);
const double clhs1 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs2 =             C(0,2)*DN(0,0);
const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
const double clhs4 =             pow(DN(0,0), 2);
const double clhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs7 =             rho*stab_c2*sqrt(pow(clhs5, 2) + pow(clhs6, 2));
const double clhs8 =             clhs7*h/stab_c1 + mu;
const double clhs9 =             DN(0,0)*clhs5 + DN(0,1)*clhs6;
const double clhs10 =             N[0]*rho;
const double clhs11 =             clhs0*(0.5*max_spectral_radius - 1.5);
const double clhs12 =             1.0/dt;
const double clhs13 =             0.5*max_spectral_radius - 1.5;
const double clhs14 =             clhs12*clhs13;
const double clhs15 =             clhs14/(-clhs0 - clhs11 + 0.5);
const double clhs16 =             -N[0]*clhs15 + clhs9;
const double clhs17 =             clhs12*rho;
const double clhs18 =             1.0/(clhs17*dyn_tau + clhs7/h + mu*stab_c1/pow(h, 2));
const double clhs19 =             clhs18*pow(rho, 2);
const double clhs20 =             clhs19*clhs9;
const double clhs21 =             1.0/(clhs0 + clhs11 - 0.5);
const double clhs22 =             clhs13*clhs17*clhs21;
const double clhs23 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs24 =             clhs19*clhs23;
const double clhs25 =             N[0]*clhs24;
const double clhs26 =             pow(N[0], 2)*clhs22 + clhs10*clhs9 + clhs16*clhs20 + clhs16*clhs25;
const double clhs27 =             C(0,1)*DN(0,1) + clhs2;
const double clhs28 =             C(1,2)*DN(0,1);
const double clhs29 =             C(2,2)*DN(0,0) + clhs28;
const double clhs30 =             DN(0,0)*clhs8;
const double clhs31 =             DN(0,1)*clhs30;
const double clhs32 =             clhs18*clhs23;
const double clhs33 =             clhs18*rho;
const double clhs34 =             clhs33*clhs9;
const double clhs35 =             -N[0] + clhs10*clhs32 + clhs34;
const double clhs36 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs37 =             C(0,2)*DN(1,0);
const double clhs38 =             C(2,2)*DN(1,1) + clhs37;
const double clhs39 =             DN(0,0)*DN(1,0);
const double clhs40 =             clhs10*clhs14*clhs21;
const double clhs41 =             N[1]*clhs40;
const double clhs42 =             clhs39*clhs8 + clhs41;
const double clhs43 =             DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs44 =             -N[1]*clhs15 + clhs43;
const double clhs45 =             clhs10*clhs43 + clhs20*clhs44 + clhs25*clhs44;
const double clhs46 =             C(0,1)*DN(1,1) + clhs37;
const double clhs47 =             C(1,2)*DN(1,1);
const double clhs48 =             C(2,2)*DN(1,0) + clhs47;
const double clhs49 =             DN(1,1)*clhs30;
const double clhs50 =             DN(0,0)*N[1];
const double clhs51 =             DN(1,0)*N[0];
const double clhs52 =             clhs23*clhs33;
const double clhs53 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs54 =             C(0,2)*DN(2,0);
const double clhs55 =             C(2,2)*DN(2,1) + clhs54;
const double clhs56 =             DN(0,0)*DN(2,0);
const double clhs57 =             N[2]*clhs40;
const double clhs58 =             clhs56*clhs8 + clhs57;
const double clhs59 =             DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs60 =             -N[2]*clhs15 + clhs59;
const double clhs61 =             clhs10*clhs59 + clhs20*clhs60 + clhs25*clhs60;
const double clhs62 =             C(0,1)*DN(2,1) + clhs54;
const double clhs63 =             C(1,2)*DN(2,1);
const double clhs64 =             C(2,2)*DN(2,0) + clhs63;
const double clhs65 =             DN(2,1)*clhs30;
const double clhs66 =             DN(0,0)*N[2];
const double clhs67 =             DN(2,0)*N[0];
const double clhs68 =             C(0,1)*DN(0,0) + clhs28;
const double clhs69 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs70 =             pow(DN(0,1), 2);
const double clhs71 =             C(0,1)*DN(1,0) + clhs47;
const double clhs72 =             DN(0,1)*clhs8;
const double clhs73 =             DN(1,0)*clhs72;
const double clhs74 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs75 =             DN(0,1)*DN(1,1);
const double clhs76 =             clhs41 + clhs75*clhs8;
const double clhs77 =             DN(0,1)*N[1];
const double clhs78 =             DN(1,1)*N[0];
const double clhs79 =             C(0,1)*DN(2,0) + clhs63;
const double clhs80 =             DN(2,0)*clhs72;
const double clhs81 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs82 =             DN(0,1)*DN(2,1);
const double clhs83 =             clhs57 + clhs8*clhs82;
const double clhs84 =             DN(0,1)*N[2];
const double clhs85 =             DN(2,1)*N[0];
const double clhs86 =             clhs16*clhs33;
const double clhs87 =             clhs0*(N[0] + clhs86);
const double clhs88 =             clhs33*clhs44;
const double clhs89 =             clhs18*(clhs39 + clhs75);
const double clhs90 =             clhs33*clhs60;
const double clhs91 =             clhs18*(clhs56 + clhs82);
const double clhs92 =             N[1]*rho;
const double clhs93 =             clhs19*clhs43;
const double clhs94 =             N[1]*clhs24;
const double clhs95 =             clhs16*clhs93 + clhs16*clhs94 + clhs9*clhs92;
const double clhs96 =             clhs33*clhs43;
const double clhs97 =             pow(DN(1,0), 2);
const double clhs98 =             pow(N[1], 2)*clhs22 + clhs43*clhs92 + clhs44*clhs93 + clhs44*clhs94;
const double clhs99 =             DN(1,0)*clhs8;
const double clhs100 =             DN(1,1)*clhs99;
const double clhs101 =             -N[1] + clhs32*clhs92 + clhs96;
const double clhs102 =             DN(1,0)*DN(2,0);
const double clhs103 =             N[1]*N[2]*clhs22;
const double clhs104 =             clhs102*clhs8 + clhs103;
const double clhs105 =             clhs59*clhs92 + clhs60*clhs93 + clhs60*clhs94;
const double clhs106 =             DN(2,1)*clhs99;
const double clhs107 =             DN(1,0)*N[2];
const double clhs108 =             DN(2,0)*N[1];
const double clhs109 =             pow(DN(1,1), 2);
const double clhs110 =             DN(2,0)*clhs8;
const double clhs111 =             DN(1,1)*clhs110;
const double clhs112 =             DN(1,1)*DN(2,1);
const double clhs113 =             clhs103 + clhs112*clhs8;
const double clhs114 =             DN(1,1)*N[2];
const double clhs115 =             DN(2,1)*N[1];
const double clhs116 =             clhs0*(N[1] + clhs88);
const double clhs117 =             clhs18*(clhs102 + clhs112);
const double clhs118 =             N[2]*rho;
const double clhs119 =             clhs19*clhs59;
const double clhs120 =             N[2]*clhs24;
const double clhs121 =             clhs118*clhs9 + clhs119*clhs16 + clhs120*clhs16;
const double clhs122 =             clhs33*clhs59;
const double clhs123 =             clhs118*clhs43 + clhs119*clhs44 + clhs120*clhs44;
const double clhs124 =             pow(DN(2,0), 2);
const double clhs125 =             pow(N[2], 2)*clhs22 + clhs118*clhs59 + clhs119*clhs60 + clhs120*clhs60;
const double clhs126 =             DN(2,1)*clhs110;
const double clhs127 =             -N[2] + clhs118*clhs32 + clhs122;
const double clhs128 =             pow(DN(2,1), 2);
const double clhs129 =             clhs0*(N[2] + clhs90);
            lhs(0,0)=clhs0*(DN(0,0)*clhs1 + DN(0,1)*clhs3 + clhs26 + clhs4*clhs8);
            lhs(0,1)=clhs0*(DN(0,0)*clhs27 + DN(0,1)*clhs29 + clhs31);
            lhs(0,2)=DN(0,0)*clhs35;
            lhs(0,3)=clhs0*(DN(0,0)*clhs36 + DN(0,1)*clhs38 + clhs42 + clhs45);
            lhs(0,4)=clhs0*(DN(0,0)*clhs46 + DN(0,1)*clhs48 + clhs49);
            lhs(0,5)=DN(1,0)*clhs34 - clhs50 + clhs51*clhs52;
            lhs(0,6)=clhs0*(DN(0,0)*clhs53 + DN(0,1)*clhs55 + clhs58 + clhs61);
            lhs(0,7)=clhs0*(DN(0,0)*clhs62 + DN(0,1)*clhs64 + clhs65);
            lhs(0,8)=DN(2,0)*clhs34 + clhs52*clhs67 - clhs66;
            lhs(1,0)=clhs0*(DN(0,0)*clhs3 + DN(0,1)*clhs68 + clhs31);
            lhs(1,1)=clhs0*(DN(0,0)*clhs29 + DN(0,1)*clhs69 + clhs26 + clhs70*clhs8);
            lhs(1,2)=DN(0,1)*clhs35;
            lhs(1,3)=clhs0*(DN(0,0)*clhs38 + DN(0,1)*clhs71 + clhs73);
            lhs(1,4)=clhs0*(DN(0,0)*clhs48 + DN(0,1)*clhs74 + clhs45 + clhs76);
            lhs(1,5)=DN(1,1)*clhs34 + clhs52*clhs78 - clhs77;
            lhs(1,6)=clhs0*(DN(0,0)*clhs55 + DN(0,1)*clhs79 + clhs80);
            lhs(1,7)=clhs0*(DN(0,0)*clhs64 + DN(0,1)*clhs81 + clhs61 + clhs83);
            lhs(1,8)=DN(2,1)*clhs34 + clhs52*clhs85 - clhs84;
            lhs(2,0)=DN(0,0)*clhs87;
            lhs(2,1)=DN(0,1)*clhs87;
            lhs(2,2)=clhs18*(clhs4 + clhs70);
            lhs(2,3)=clhs0*(DN(0,0)*clhs88 + clhs51);
            lhs(2,4)=clhs0*(DN(0,1)*clhs88 + clhs78);
            lhs(2,5)=clhs89;
            lhs(2,6)=clhs0*(DN(0,0)*clhs90 + clhs67);
            lhs(2,7)=clhs0*(DN(0,1)*clhs90 + clhs85);
            lhs(2,8)=clhs91;
            lhs(3,0)=clhs0*(DN(1,0)*clhs1 + DN(1,1)*clhs3 + clhs42 + clhs95);
            lhs(3,1)=clhs0*(DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs73);
            lhs(3,2)=DN(0,0)*clhs96 + clhs50*clhs52 - clhs51;
            lhs(3,3)=clhs0*(DN(1,0)*clhs36 + DN(1,1)*clhs38 + clhs8*clhs97 + clhs98);
            lhs(3,4)=clhs0*(DN(1,0)*clhs46 + DN(1,1)*clhs48 + clhs100);
            lhs(3,5)=DN(1,0)*clhs101;
            lhs(3,6)=clhs0*(DN(1,0)*clhs53 + DN(1,1)*clhs55 + clhs104 + clhs105);
            lhs(3,7)=clhs0*(DN(1,0)*clhs62 + DN(1,1)*clhs64 + clhs106);
            lhs(3,8)=DN(2,0)*clhs96 - clhs107 + clhs108*clhs52;
            lhs(4,0)=clhs0*(DN(1,0)*clhs3 + DN(1,1)*clhs68 + clhs49);
            lhs(4,1)=clhs0*(DN(1,0)*clhs29 + DN(1,1)*clhs69 + clhs76 + clhs95);
            lhs(4,2)=DN(0,1)*clhs96 + clhs52*clhs77 - clhs78;
            lhs(4,3)=clhs0*(DN(1,0)*clhs38 + DN(1,1)*clhs71 + clhs100);
            lhs(4,4)=clhs0*(DN(1,0)*clhs48 + DN(1,1)*clhs74 + clhs109*clhs8 + clhs98);
            lhs(4,5)=DN(1,1)*clhs101;
            lhs(4,6)=clhs0*(DN(1,0)*clhs55 + DN(1,1)*clhs79 + clhs111);
            lhs(4,7)=clhs0*(DN(1,0)*clhs64 + DN(1,1)*clhs81 + clhs105 + clhs113);
            lhs(4,8)=DN(2,1)*clhs96 - clhs114 + clhs115*clhs52;
            lhs(5,0)=clhs0*(DN(1,0)*clhs86 + clhs50);
            lhs(5,1)=clhs0*(DN(1,1)*clhs86 + clhs77);
            lhs(5,2)=clhs89;
            lhs(5,3)=DN(1,0)*clhs116;
            lhs(5,4)=DN(1,1)*clhs116;
            lhs(5,5)=clhs18*(clhs109 + clhs97);
            lhs(5,6)=clhs0*(DN(1,0)*clhs90 + clhs108);
            lhs(5,7)=clhs0*(DN(1,1)*clhs90 + clhs115);
            lhs(5,8)=clhs117;
            lhs(6,0)=clhs0*(DN(2,0)*clhs1 + DN(2,1)*clhs3 + clhs121 + clhs58);
            lhs(6,1)=clhs0*(DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs80);
            lhs(6,2)=DN(0,0)*clhs122 + clhs52*clhs66 - clhs67;
            lhs(6,3)=clhs0*(DN(2,0)*clhs36 + DN(2,1)*clhs38 + clhs104 + clhs123);
            lhs(6,4)=clhs0*(DN(2,0)*clhs46 + DN(2,1)*clhs48 + clhs111);
            lhs(6,5)=DN(1,0)*clhs122 + clhs107*clhs52 - clhs108;
            lhs(6,6)=clhs0*(DN(2,0)*clhs53 + DN(2,1)*clhs55 + clhs124*clhs8 + clhs125);
            lhs(6,7)=clhs0*(DN(2,0)*clhs62 + DN(2,1)*clhs64 + clhs126);
            lhs(6,8)=DN(2,0)*clhs127;
            lhs(7,0)=clhs0*(DN(2,0)*clhs3 + DN(2,1)*clhs68 + clhs65);
            lhs(7,1)=clhs0*(DN(2,0)*clhs29 + DN(2,1)*clhs69 + clhs121 + clhs83);
            lhs(7,2)=DN(0,1)*clhs122 + clhs52*clhs84 - clhs85;
            lhs(7,3)=clhs0*(DN(2,0)*clhs38 + DN(2,1)*clhs71 + clhs106);
            lhs(7,4)=clhs0*(DN(2,0)*clhs48 + DN(2,1)*clhs74 + clhs113 + clhs123);
            lhs(7,5)=DN(1,1)*clhs122 + clhs114*clhs52 - clhs115;
            lhs(7,6)=clhs0*(DN(2,0)*clhs55 + DN(2,1)*clhs79 + clhs126);
            lhs(7,7)=clhs0*(DN(2,0)*clhs64 + DN(2,1)*clhs81 + clhs125 + clhs128*clhs8);
            lhs(7,8)=DN(2,1)*clhs127;
            lhs(8,0)=clhs0*(DN(2,0)*clhs86 + clhs66);
            lhs(8,1)=clhs0*(DN(2,1)*clhs86 + clhs84);
            lhs(8,2)=clhs91;
            lhs(8,3)=clhs0*(DN(2,0)*clhs88 + clhs107);
            lhs(8,4)=clhs0*(DN(2,1)*clhs88 + clhs114);
            lhs(8,5)=clhs117;
            lhs(8,6)=DN(2,0)*clhs129;
            lhs(8,7)=DN(2,1)*clhs129;
            lhs(8,8)=clhs18*(clhs124 + clhs128);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             1.0/(max_spectral_radius + 1);
const double clhs1 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs2 =             C(0,3)*DN(0,0);
const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
const double clhs4 =             C(0,5)*DN(0,0);
const double clhs5 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs4;
const double clhs6 =             pow(DN(0,0), 2);
const double clhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs10 =             rho*stab_c2*sqrt(pow(clhs7, 2) + pow(clhs8, 2) + pow(clhs9, 2));
const double clhs11 =             clhs10*h/stab_c1 + mu;
const double clhs12 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs13 =             N[0]*rho;
const double clhs14 =             clhs0*(0.5*max_spectral_radius - 1.5);
const double clhs15 =             1.0/dt;
const double clhs16 =             0.5*max_spectral_radius - 1.5;
const double clhs17 =             clhs15*clhs16;
const double clhs18 =             clhs17/(-clhs0 - clhs14 + 0.5);
const double clhs19 =             -N[0]*clhs18 + clhs12;
const double clhs20 =             clhs15*rho;
const double clhs21 =             1.0/(clhs10/h + clhs20*dyn_tau + mu*stab_c1/pow(h, 2));
const double clhs22 =             clhs21*pow(rho, 2);
const double clhs23 =             clhs12*clhs22;
const double clhs24 =             1.0/(clhs0 + clhs14 - 0.5);
const double clhs25 =             clhs16*clhs20*clhs24;
const double clhs26 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs27 =             clhs22*clhs26;
const double clhs28 =             N[0]*clhs27;
const double clhs29 =             pow(N[0], 2)*clhs25 + clhs12*clhs13 + clhs19*clhs23 + clhs19*clhs28;
const double clhs30 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
const double clhs31 =             C(1,3)*DN(0,1);
const double clhs32 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs31;
const double clhs33 =             C(3,5)*DN(0,0);
const double clhs34 =             C(4,5)*DN(0,2);
const double clhs35 =             C(1,5)*DN(0,1) + clhs33 + clhs34;
const double clhs36 =             DN(0,0)*clhs11;
const double clhs37 =             DN(0,1)*clhs36;
const double clhs38 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs4;
const double clhs39 =             C(3,4)*DN(0,1);
const double clhs40 =             C(2,3)*DN(0,2) + clhs33 + clhs39;
const double clhs41 =             C(2,5)*DN(0,2);
const double clhs42 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs41;
const double clhs43 =             DN(0,2)*clhs36;
const double clhs44 =             clhs21*clhs26;
const double clhs45 =             clhs21*rho;
const double clhs46 =             clhs12*clhs45;
const double clhs47 =             -N[0] + clhs13*clhs44 + clhs46;
const double clhs48 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs49 =             C(0,3)*DN(1,0);
const double clhs50 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs49;
const double clhs51 =             C(0,5)*DN(1,0);
const double clhs52 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs51;
const double clhs53 =             DN(0,0)*DN(1,0);
const double clhs54 =             clhs13*clhs17*clhs24;
const double clhs55 =             N[1]*clhs54;
const double clhs56 =             clhs11*clhs53 + clhs55;
const double clhs57 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs58 =             -N[1]*clhs18 + clhs57;
const double clhs59 =             clhs13*clhs57 + clhs23*clhs58 + clhs28*clhs58;
const double clhs60 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs49;
const double clhs61 =             C(1,3)*DN(1,1);
const double clhs62 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs61;
const double clhs63 =             C(3,5)*DN(1,0);
const double clhs64 =             C(4,5)*DN(1,2);
const double clhs65 =             C(1,5)*DN(1,1) + clhs63 + clhs64;
const double clhs66 =             DN(1,1)*clhs36;
const double clhs67 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs51;
const double clhs68 =             C(3,4)*DN(1,1);
const double clhs69 =             C(2,3)*DN(1,2) + clhs63 + clhs68;
const double clhs70 =             C(2,5)*DN(1,2);
const double clhs71 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs70;
const double clhs72 =             DN(1,2)*clhs36;
const double clhs73 =             DN(0,0)*N[1];
const double clhs74 =             DN(1,0)*N[0];
const double clhs75 =             clhs26*clhs45;
const double clhs76 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs77 =             C(0,3)*DN(2,0);
const double clhs78 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs77;
const double clhs79 =             C(0,5)*DN(2,0);
const double clhs80 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs79;
const double clhs81 =             DN(0,0)*DN(2,0);
const double clhs82 =             N[2]*clhs54;
const double clhs83 =             clhs11*clhs81 + clhs82;
const double clhs84 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs85 =             -N[2]*clhs18 + clhs84;
const double clhs86 =             clhs13*clhs84 + clhs23*clhs85 + clhs28*clhs85;
const double clhs87 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs77;
const double clhs88 =             C(1,3)*DN(2,1);
const double clhs89 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs88;
const double clhs90 =             C(3,5)*DN(2,0);
const double clhs91 =             C(4,5)*DN(2,2);
const double clhs92 =             C(1,5)*DN(2,1) + clhs90 + clhs91;
const double clhs93 =             DN(2,1)*clhs36;
const double clhs94 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs79;
const double clhs95 =             C(3,4)*DN(2,1);
const double clhs96 =             C(2,3)*DN(2,2) + clhs90 + clhs95;
const double clhs97 =             C(2,5)*DN(2,2);
const double clhs98 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs97;
const double clhs99 =             DN(2,2)*clhs36;
const double clhs100 =             DN(0,0)*N[2];
const double clhs101 =             DN(2,0)*N[0];
const double clhs102 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs103 =             C(0,3)*DN(3,0);
const double clhs104 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs103;
const double clhs105 =             C(0,5)*DN(3,0);
const double clhs106 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs105;
const double clhs107 =             DN(0,0)*DN(3,0);
const double clhs108 =             N[3]*clhs54;
const double clhs109 =             clhs107*clhs11 + clhs108;
const double clhs110 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs111 =             -N[3]*clhs18 + clhs110;
const double clhs112 =             clhs110*clhs13 + clhs111*clhs23 + clhs111*clhs28;
const double clhs113 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs103;
const double clhs114 =             C(1,3)*DN(3,1);
const double clhs115 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs114;
const double clhs116 =             C(3,5)*DN(3,0);
const double clhs117 =             C(4,5)*DN(3,2);
const double clhs118 =             C(1,5)*DN(3,1) + clhs116 + clhs117;
const double clhs119 =             DN(3,1)*clhs36;
const double clhs120 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs105;
const double clhs121 =             C(3,4)*DN(3,1);
const double clhs122 =             C(2,3)*DN(3,2) + clhs116 + clhs121;
const double clhs123 =             C(2,5)*DN(3,2);
const double clhs124 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs123;
const double clhs125 =             DN(3,2)*clhs36;
const double clhs126 =             DN(0,0)*N[3];
const double clhs127 =             DN(3,0)*N[0];
const double clhs128 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs31;
const double clhs129 =             C(0,4)*DN(0,0) + clhs34 + clhs39;
const double clhs130 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs131 =             C(1,4)*DN(0,1);
const double clhs132 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs131;
const double clhs133 =             pow(DN(0,1), 2);
const double clhs134 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs131;
const double clhs135 =             C(2,4)*DN(0,2);
const double clhs136 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs135;
const double clhs137 =             DN(0,1)*clhs11;
const double clhs138 =             DN(0,2)*clhs137;
const double clhs139 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs61;
const double clhs140 =             C(0,4)*DN(1,0) + clhs64 + clhs68;
const double clhs141 =             DN(1,0)*clhs137;
const double clhs142 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs143 =             C(1,4)*DN(1,1);
const double clhs144 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs143;
const double clhs145 =             DN(0,1)*DN(1,1);
const double clhs146 =             clhs11*clhs145;
const double clhs147 =             clhs55 + clhs59;
const double clhs148 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs143;
const double clhs149 =             C(2,4)*DN(1,2);
const double clhs150 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs149;
const double clhs151 =             DN(1,2)*clhs137;
const double clhs152 =             DN(0,1)*N[1];
const double clhs153 =             DN(1,1)*N[0];
const double clhs154 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs88;
const double clhs155 =             C(0,4)*DN(2,0) + clhs91 + clhs95;
const double clhs156 =             DN(2,0)*clhs137;
const double clhs157 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs158 =             C(1,4)*DN(2,1);
const double clhs159 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs158;
const double clhs160 =             DN(0,1)*DN(2,1);
const double clhs161 =             clhs11*clhs160;
const double clhs162 =             clhs82 + clhs86;
const double clhs163 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs158;
const double clhs164 =             C(2,4)*DN(2,2);
const double clhs165 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs164;
const double clhs166 =             DN(2,2)*clhs137;
const double clhs167 =             DN(0,1)*N[2];
const double clhs168 =             DN(2,1)*N[0];
const double clhs169 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs114;
const double clhs170 =             C(0,4)*DN(3,0) + clhs117 + clhs121;
const double clhs171 =             DN(3,0)*clhs137;
const double clhs172 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs173 =             C(1,4)*DN(3,1);
const double clhs174 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs173;
const double clhs175 =             DN(0,1)*DN(3,1);
const double clhs176 =             clhs11*clhs175;
const double clhs177 =             clhs108 + clhs112;
const double clhs178 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs173;
const double clhs179 =             C(2,4)*DN(3,2);
const double clhs180 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs179;
const double clhs181 =             DN(3,2)*clhs137;
const double clhs182 =             DN(0,1)*N[3];
const double clhs183 =             DN(3,1)*N[0];
const double clhs184 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs41;
const double clhs185 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs135;
const double clhs186 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs187 =             pow(DN(0,2), 2);
const double clhs188 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs70;
const double clhs189 =             DN(0,2)*clhs11;
const double clhs190 =             DN(1,0)*clhs189;
const double clhs191 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs149;
const double clhs192 =             DN(1,1)*clhs189;
const double clhs193 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs194 =             DN(0,2)*DN(1,2);
const double clhs195 =             clhs11*clhs194;
const double clhs196 =             DN(0,2)*N[1];
const double clhs197 =             DN(1,2)*N[0];
const double clhs198 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs97;
const double clhs199 =             DN(2,0)*clhs189;
const double clhs200 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs164;
const double clhs201 =             DN(2,1)*clhs189;
const double clhs202 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs203 =             DN(0,2)*DN(2,2);
const double clhs204 =             clhs11*clhs203;
const double clhs205 =             DN(0,2)*N[2];
const double clhs206 =             DN(2,2)*N[0];
const double clhs207 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs123;
const double clhs208 =             DN(3,0)*clhs189;
const double clhs209 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs179;
const double clhs210 =             DN(3,1)*clhs189;
const double clhs211 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs212 =             DN(0,2)*DN(3,2);
const double clhs213 =             clhs11*clhs212;
const double clhs214 =             DN(0,2)*N[3];
const double clhs215 =             DN(3,2)*N[0];
const double clhs216 =             clhs19*clhs45;
const double clhs217 =             clhs0*(N[0] + clhs216);
const double clhs218 =             clhs45*clhs58;
const double clhs219 =             clhs21*(clhs145 + clhs194 + clhs53);
const double clhs220 =             clhs45*clhs85;
const double clhs221 =             clhs21*(clhs160 + clhs203 + clhs81);
const double clhs222 =             clhs111*clhs45;
const double clhs223 =             clhs21*(clhs107 + clhs175 + clhs212);
const double clhs224 =             N[1]*rho;
const double clhs225 =             clhs22*clhs57;
const double clhs226 =             N[1]*clhs27;
const double clhs227 =             clhs12*clhs224 + clhs19*clhs225 + clhs19*clhs226;
const double clhs228 =             clhs45*clhs57;
const double clhs229 =             pow(DN(1,0), 2);
const double clhs230 =             pow(N[1], 2)*clhs25 + clhs224*clhs57 + clhs225*clhs58 + clhs226*clhs58;
const double clhs231 =             DN(1,0)*clhs11;
const double clhs232 =             DN(1,1)*clhs231;
const double clhs233 =             DN(1,2)*clhs231;
const double clhs234 =             -N[1] + clhs224*clhs44 + clhs228;
const double clhs235 =             DN(1,0)*DN(2,0);
const double clhs236 =             N[1]*clhs25;
const double clhs237 =             N[2]*clhs236;
const double clhs238 =             clhs11*clhs235 + clhs237;
const double clhs239 =             clhs224*clhs84 + clhs225*clhs85 + clhs226*clhs85;
const double clhs240 =             DN(2,1)*clhs231;
const double clhs241 =             DN(2,2)*clhs231;
const double clhs242 =             DN(1,0)*N[2];
const double clhs243 =             DN(2,0)*N[1];
const double clhs244 =             DN(1,0)*DN(3,0);
const double clhs245 =             N[3]*clhs236;
const double clhs246 =             clhs11*clhs244 + clhs245;
const double clhs247 =             clhs110*clhs224 + clhs111*clhs225 + clhs111*clhs226;
const double clhs248 =             DN(3,1)*clhs231;
const double clhs249 =             DN(3,2)*clhs231;
const double clhs250 =             DN(1,0)*N[3];
const double clhs251 =             DN(3,0)*N[1];
const double clhs252 =             clhs227 + clhs55;
const double clhs253 =             pow(DN(1,1), 2);
const double clhs254 =             DN(1,1)*clhs11;
const double clhs255 =             DN(1,2)*clhs254;
const double clhs256 =             DN(2,0)*clhs254;
const double clhs257 =             DN(1,1)*DN(2,1);
const double clhs258 =             clhs11*clhs257;
const double clhs259 =             clhs237 + clhs239;
const double clhs260 =             DN(2,2)*clhs254;
const double clhs261 =             DN(1,1)*N[2];
const double clhs262 =             DN(2,1)*N[1];
const double clhs263 =             DN(3,0)*clhs254;
const double clhs264 =             DN(1,1)*DN(3,1);
const double clhs265 =             clhs11*clhs264;
const double clhs266 =             clhs245 + clhs247;
const double clhs267 =             DN(3,2)*clhs254;
const double clhs268 =             DN(1,1)*N[3];
const double clhs269 =             DN(3,1)*N[1];
const double clhs270 =             pow(DN(1,2), 2);
const double clhs271 =             DN(1,2)*clhs11;
const double clhs272 =             DN(2,0)*clhs271;
const double clhs273 =             DN(2,1)*clhs271;
const double clhs274 =             DN(1,2)*DN(2,2);
const double clhs275 =             clhs11*clhs274;
const double clhs276 =             DN(1,2)*N[2];
const double clhs277 =             DN(2,2)*N[1];
const double clhs278 =             DN(3,0)*clhs271;
const double clhs279 =             DN(3,1)*clhs271;
const double clhs280 =             DN(1,2)*DN(3,2);
const double clhs281 =             clhs11*clhs280;
const double clhs282 =             DN(1,2)*N[3];
const double clhs283 =             DN(3,2)*N[1];
const double clhs284 =             clhs0*(N[1] + clhs218);
const double clhs285 =             clhs21*(clhs235 + clhs257 + clhs274);
const double clhs286 =             clhs21*(clhs244 + clhs264 + clhs280);
const double clhs287 =             N[2]*rho;
const double clhs288 =             clhs22*clhs84;
const double clhs289 =             N[2]*clhs27;
const double clhs290 =             clhs12*clhs287 + clhs19*clhs288 + clhs19*clhs289;
const double clhs291 =             clhs45*clhs84;
const double clhs292 =             clhs287*clhs57 + clhs288*clhs58 + clhs289*clhs58;
const double clhs293 =             pow(DN(2,0), 2);
const double clhs294 =             pow(N[2], 2)*clhs25 + clhs287*clhs84 + clhs288*clhs85 + clhs289*clhs85;
const double clhs295 =             DN(2,0)*clhs11;
const double clhs296 =             DN(2,1)*clhs295;
const double clhs297 =             DN(2,2)*clhs295;
const double clhs298 =             -N[2] + clhs287*clhs44 + clhs291;
const double clhs299 =             DN(2,0)*DN(3,0);
const double clhs300 =             N[2]*N[3]*clhs25;
const double clhs301 =             clhs11*clhs299 + clhs300;
const double clhs302 =             clhs110*clhs287 + clhs111*clhs288 + clhs111*clhs289;
const double clhs303 =             DN(3,1)*clhs295;
const double clhs304 =             DN(3,2)*clhs295;
const double clhs305 =             DN(2,0)*N[3];
const double clhs306 =             DN(3,0)*N[2];
const double clhs307 =             clhs290 + clhs82;
const double clhs308 =             clhs237 + clhs292;
const double clhs309 =             pow(DN(2,1), 2);
const double clhs310 =             DN(2,1)*clhs11;
const double clhs311 =             DN(2,2)*clhs310;
const double clhs312 =             DN(3,0)*clhs310;
const double clhs313 =             DN(2,1)*DN(3,1);
const double clhs314 =             clhs11*clhs313;
const double clhs315 =             clhs300 + clhs302;
const double clhs316 =             DN(3,2)*clhs310;
const double clhs317 =             DN(2,1)*N[3];
const double clhs318 =             DN(3,1)*N[2];
const double clhs319 =             pow(DN(2,2), 2);
const double clhs320 =             DN(2,2)*clhs11;
const double clhs321 =             DN(3,0)*clhs320;
const double clhs322 =             DN(3,1)*clhs320;
const double clhs323 =             DN(2,2)*DN(3,2);
const double clhs324 =             clhs11*clhs323;
const double clhs325 =             DN(2,2)*N[3];
const double clhs326 =             DN(3,2)*N[2];
const double clhs327 =             clhs0*(N[2] + clhs220);
const double clhs328 =             clhs21*(clhs299 + clhs313 + clhs323);
const double clhs329 =             N[3]*rho;
const double clhs330 =             clhs110*clhs22;
const double clhs331 =             N[3]*clhs27;
const double clhs332 =             clhs12*clhs329 + clhs19*clhs330 + clhs19*clhs331;
const double clhs333 =             clhs110*clhs45;
const double clhs334 =             clhs329*clhs57 + clhs330*clhs58 + clhs331*clhs58;
const double clhs335 =             clhs329*clhs84 + clhs330*clhs85 + clhs331*clhs85;
const double clhs336 =             pow(DN(3,0), 2);
const double clhs337 =             pow(N[3], 2)*clhs25 + clhs110*clhs329 + clhs111*clhs330 + clhs111*clhs331;
const double clhs338 =             DN(3,0)*clhs11;
const double clhs339 =             DN(3,1)*clhs338;
const double clhs340 =             DN(3,2)*clhs338;
const double clhs341 =             -N[3] + clhs329*clhs44 + clhs333;
const double clhs342 =             clhs108 + clhs332;
const double clhs343 =             clhs245 + clhs334;
const double clhs344 =             clhs300 + clhs335;
const double clhs345 =             pow(DN(3,1), 2);
const double clhs346 =             DN(3,1)*DN(3,2)*clhs11;
const double clhs347 =             pow(DN(3,2), 2);
const double clhs348 =             clhs0*(N[3] + clhs222);
            lhs(0,0)=clhs0*(DN(0,0)*clhs1 + DN(0,1)*clhs3 + DN(0,2)*clhs5 + clhs11*clhs6 + clhs29);
            lhs(0,1)=clhs0*(DN(0,0)*clhs30 + DN(0,1)*clhs32 + DN(0,2)*clhs35 + clhs37);
            lhs(0,2)=clhs0*(DN(0,0)*clhs38 + DN(0,1)*clhs40 + DN(0,2)*clhs42 + clhs43);
            lhs(0,3)=DN(0,0)*clhs47;
            lhs(0,4)=clhs0*(DN(0,0)*clhs48 + DN(0,1)*clhs50 + DN(0,2)*clhs52 + clhs56 + clhs59);
            lhs(0,5)=clhs0*(DN(0,0)*clhs60 + DN(0,1)*clhs62 + DN(0,2)*clhs65 + clhs66);
            lhs(0,6)=clhs0*(DN(0,0)*clhs67 + DN(0,1)*clhs69 + DN(0,2)*clhs71 + clhs72);
            lhs(0,7)=DN(1,0)*clhs46 - clhs73 + clhs74*clhs75;
            lhs(0,8)=clhs0*(DN(0,0)*clhs76 + DN(0,1)*clhs78 + DN(0,2)*clhs80 + clhs83 + clhs86);
            lhs(0,9)=clhs0*(DN(0,0)*clhs87 + DN(0,1)*clhs89 + DN(0,2)*clhs92 + clhs93);
            lhs(0,10)=clhs0*(DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs99);
            lhs(0,11)=DN(2,0)*clhs46 - clhs100 + clhs101*clhs75;
            lhs(0,12)=clhs0*(DN(0,0)*clhs102 + DN(0,1)*clhs104 + DN(0,2)*clhs106 + clhs109 + clhs112);
            lhs(0,13)=clhs0*(DN(0,0)*clhs113 + DN(0,1)*clhs115 + DN(0,2)*clhs118 + clhs119);
            lhs(0,14)=clhs0*(DN(0,0)*clhs120 + DN(0,1)*clhs122 + DN(0,2)*clhs124 + clhs125);
            lhs(0,15)=DN(3,0)*clhs46 - clhs126 + clhs127*clhs75;
            lhs(1,0)=clhs0*(DN(0,0)*clhs3 + DN(0,1)*clhs128 + DN(0,2)*clhs129 + clhs37);
            lhs(1,1)=clhs0*(DN(0,0)*clhs32 + DN(0,1)*clhs130 + DN(0,2)*clhs132 + clhs11*clhs133 + clhs29);
            lhs(1,2)=clhs0*(DN(0,0)*clhs40 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs138);
            lhs(1,3)=DN(0,1)*clhs47;
            lhs(1,4)=clhs0*(DN(0,0)*clhs50 + DN(0,1)*clhs139 + DN(0,2)*clhs140 + clhs141);
            lhs(1,5)=clhs0*(DN(0,0)*clhs62 + DN(0,1)*clhs142 + DN(0,2)*clhs144 + clhs146 + clhs147);
            lhs(1,6)=clhs0*(DN(0,0)*clhs69 + DN(0,1)*clhs148 + DN(0,2)*clhs150 + clhs151);
            lhs(1,7)=DN(1,1)*clhs46 - clhs152 + clhs153*clhs75;
            lhs(1,8)=clhs0*(DN(0,0)*clhs78 + DN(0,1)*clhs154 + DN(0,2)*clhs155 + clhs156);
            lhs(1,9)=clhs0*(DN(0,0)*clhs89 + DN(0,1)*clhs157 + DN(0,2)*clhs159 + clhs161 + clhs162);
            lhs(1,10)=clhs0*(DN(0,0)*clhs96 + DN(0,1)*clhs163 + DN(0,2)*clhs165 + clhs166);
            lhs(1,11)=DN(2,1)*clhs46 - clhs167 + clhs168*clhs75;
            lhs(1,12)=clhs0*(DN(0,0)*clhs104 + DN(0,1)*clhs169 + DN(0,2)*clhs170 + clhs171);
            lhs(1,13)=clhs0*(DN(0,0)*clhs115 + DN(0,1)*clhs172 + DN(0,2)*clhs174 + clhs176 + clhs177);
            lhs(1,14)=clhs0*(DN(0,0)*clhs122 + DN(0,1)*clhs178 + DN(0,2)*clhs180 + clhs181);
            lhs(1,15)=DN(3,1)*clhs46 - clhs182 + clhs183*clhs75;
            lhs(2,0)=clhs0*(DN(0,0)*clhs5 + DN(0,1)*clhs129 + DN(0,2)*clhs184 + clhs43);
            lhs(2,1)=clhs0*(DN(0,0)*clhs35 + DN(0,1)*clhs132 + DN(0,2)*clhs185 + clhs138);
            lhs(2,2)=clhs0*(DN(0,0)*clhs42 + DN(0,1)*clhs136 + DN(0,2)*clhs186 + clhs11*clhs187 + clhs29);
            lhs(2,3)=DN(0,2)*clhs47;
            lhs(2,4)=clhs0*(DN(0,0)*clhs52 + DN(0,1)*clhs140 + DN(0,2)*clhs188 + clhs190);
            lhs(2,5)=clhs0*(DN(0,0)*clhs65 + DN(0,1)*clhs144 + DN(0,2)*clhs191 + clhs192);
            lhs(2,6)=clhs0*(DN(0,0)*clhs71 + DN(0,1)*clhs150 + DN(0,2)*clhs193 + clhs147 + clhs195);
            lhs(2,7)=DN(1,2)*clhs46 - clhs196 + clhs197*clhs75;
            lhs(2,8)=clhs0*(DN(0,0)*clhs80 + DN(0,1)*clhs155 + DN(0,2)*clhs198 + clhs199);
            lhs(2,9)=clhs0*(DN(0,0)*clhs92 + DN(0,1)*clhs159 + DN(0,2)*clhs200 + clhs201);
            lhs(2,10)=clhs0*(DN(0,0)*clhs98 + DN(0,1)*clhs165 + DN(0,2)*clhs202 + clhs162 + clhs204);
            lhs(2,11)=DN(2,2)*clhs46 - clhs205 + clhs206*clhs75;
            lhs(2,12)=clhs0*(DN(0,0)*clhs106 + DN(0,1)*clhs170 + DN(0,2)*clhs207 + clhs208);
            lhs(2,13)=clhs0*(DN(0,0)*clhs118 + DN(0,1)*clhs174 + DN(0,2)*clhs209 + clhs210);
            lhs(2,14)=clhs0*(DN(0,0)*clhs124 + DN(0,1)*clhs180 + DN(0,2)*clhs211 + clhs177 + clhs213);
            lhs(2,15)=DN(3,2)*clhs46 - clhs214 + clhs215*clhs75;
            lhs(3,0)=DN(0,0)*clhs217;
            lhs(3,1)=DN(0,1)*clhs217;
            lhs(3,2)=DN(0,2)*clhs217;
            lhs(3,3)=clhs21*(clhs133 + clhs187 + clhs6);
            lhs(3,4)=clhs0*(DN(0,0)*clhs218 + clhs74);
            lhs(3,5)=clhs0*(DN(0,1)*clhs218 + clhs153);
            lhs(3,6)=clhs0*(DN(0,2)*clhs218 + clhs197);
            lhs(3,7)=clhs219;
            lhs(3,8)=clhs0*(DN(0,0)*clhs220 + clhs101);
            lhs(3,9)=clhs0*(DN(0,1)*clhs220 + clhs168);
            lhs(3,10)=clhs0*(DN(0,2)*clhs220 + clhs206);
            lhs(3,11)=clhs221;
            lhs(3,12)=clhs0*(DN(0,0)*clhs222 + clhs127);
            lhs(3,13)=clhs0*(DN(0,1)*clhs222 + clhs183);
            lhs(3,14)=clhs0*(DN(0,2)*clhs222 + clhs215);
            lhs(3,15)=clhs223;
            lhs(4,0)=clhs0*(DN(1,0)*clhs1 + DN(1,1)*clhs3 + DN(1,2)*clhs5 + clhs227 + clhs56);
            lhs(4,1)=clhs0*(DN(1,0)*clhs30 + DN(1,1)*clhs32 + DN(1,2)*clhs35 + clhs141);
            lhs(4,2)=clhs0*(DN(1,0)*clhs38 + DN(1,1)*clhs40 + DN(1,2)*clhs42 + clhs190);
            lhs(4,3)=DN(0,0)*clhs228 + clhs73*clhs75 - clhs74;
            lhs(4,4)=clhs0*(DN(1,0)*clhs48 + DN(1,1)*clhs50 + DN(1,2)*clhs52 + clhs11*clhs229 + clhs230);
            lhs(4,5)=clhs0*(DN(1,0)*clhs60 + DN(1,1)*clhs62 + DN(1,2)*clhs65 + clhs232);
            lhs(4,6)=clhs0*(DN(1,0)*clhs67 + DN(1,1)*clhs69 + DN(1,2)*clhs71 + clhs233);
            lhs(4,7)=DN(1,0)*clhs234;
            lhs(4,8)=clhs0*(DN(1,0)*clhs76 + DN(1,1)*clhs78 + DN(1,2)*clhs80 + clhs238 + clhs239);
            lhs(4,9)=clhs0*(DN(1,0)*clhs87 + DN(1,1)*clhs89 + DN(1,2)*clhs92 + clhs240);
            lhs(4,10)=clhs0*(DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs241);
            lhs(4,11)=DN(2,0)*clhs228 - clhs242 + clhs243*clhs75;
            lhs(4,12)=clhs0*(DN(1,0)*clhs102 + DN(1,1)*clhs104 + DN(1,2)*clhs106 + clhs246 + clhs247);
            lhs(4,13)=clhs0*(DN(1,0)*clhs113 + DN(1,1)*clhs115 + DN(1,2)*clhs118 + clhs248);
            lhs(4,14)=clhs0*(DN(1,0)*clhs120 + DN(1,1)*clhs122 + DN(1,2)*clhs124 + clhs249);
            lhs(4,15)=DN(3,0)*clhs228 - clhs250 + clhs251*clhs75;
            lhs(5,0)=clhs0*(DN(1,0)*clhs3 + DN(1,1)*clhs128 + DN(1,2)*clhs129 + clhs66);
            lhs(5,1)=clhs0*(DN(1,0)*clhs32 + DN(1,1)*clhs130 + DN(1,2)*clhs132 + clhs146 + clhs252);
            lhs(5,2)=clhs0*(DN(1,0)*clhs40 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs192);
            lhs(5,3)=DN(0,1)*clhs228 + clhs152*clhs75 - clhs153;
            lhs(5,4)=clhs0*(DN(1,0)*clhs50 + DN(1,1)*clhs139 + DN(1,2)*clhs140 + clhs232);
            lhs(5,5)=clhs0*(DN(1,0)*clhs62 + DN(1,1)*clhs142 + DN(1,2)*clhs144 + clhs11*clhs253 + clhs230);
            lhs(5,6)=clhs0*(DN(1,0)*clhs69 + DN(1,1)*clhs148 + DN(1,2)*clhs150 + clhs255);
            lhs(5,7)=DN(1,1)*clhs234;
            lhs(5,8)=clhs0*(DN(1,0)*clhs78 + DN(1,1)*clhs154 + DN(1,2)*clhs155 + clhs256);
            lhs(5,9)=clhs0*(DN(1,0)*clhs89 + DN(1,1)*clhs157 + DN(1,2)*clhs159 + clhs258 + clhs259);
            lhs(5,10)=clhs0*(DN(1,0)*clhs96 + DN(1,1)*clhs163 + DN(1,2)*clhs165 + clhs260);
            lhs(5,11)=DN(2,1)*clhs228 - clhs261 + clhs262*clhs75;
            lhs(5,12)=clhs0*(DN(1,0)*clhs104 + DN(1,1)*clhs169 + DN(1,2)*clhs170 + clhs263);
            lhs(5,13)=clhs0*(DN(1,0)*clhs115 + DN(1,1)*clhs172 + DN(1,2)*clhs174 + clhs265 + clhs266);
            lhs(5,14)=clhs0*(DN(1,0)*clhs122 + DN(1,1)*clhs178 + DN(1,2)*clhs180 + clhs267);
            lhs(5,15)=DN(3,1)*clhs228 - clhs268 + clhs269*clhs75;
            lhs(6,0)=clhs0*(DN(1,0)*clhs5 + DN(1,1)*clhs129 + DN(1,2)*clhs184 + clhs72);
            lhs(6,1)=clhs0*(DN(1,0)*clhs35 + DN(1,1)*clhs132 + DN(1,2)*clhs185 + clhs151);
            lhs(6,2)=clhs0*(DN(1,0)*clhs42 + DN(1,1)*clhs136 + DN(1,2)*clhs186 + clhs195 + clhs252);
            lhs(6,3)=DN(0,2)*clhs228 + clhs196*clhs75 - clhs197;
            lhs(6,4)=clhs0*(DN(1,0)*clhs52 + DN(1,1)*clhs140 + DN(1,2)*clhs188 + clhs233);
            lhs(6,5)=clhs0*(DN(1,0)*clhs65 + DN(1,1)*clhs144 + DN(1,2)*clhs191 + clhs255);
            lhs(6,6)=clhs0*(DN(1,0)*clhs71 + DN(1,1)*clhs150 + DN(1,2)*clhs193 + clhs11*clhs270 + clhs230);
            lhs(6,7)=DN(1,2)*clhs234;
            lhs(6,8)=clhs0*(DN(1,0)*clhs80 + DN(1,1)*clhs155 + DN(1,2)*clhs198 + clhs272);
            lhs(6,9)=clhs0*(DN(1,0)*clhs92 + DN(1,1)*clhs159 + DN(1,2)*clhs200 + clhs273);
            lhs(6,10)=clhs0*(DN(1,0)*clhs98 + DN(1,1)*clhs165 + DN(1,2)*clhs202 + clhs259 + clhs275);
            lhs(6,11)=DN(2,2)*clhs228 - clhs276 + clhs277*clhs75;
            lhs(6,12)=clhs0*(DN(1,0)*clhs106 + DN(1,1)*clhs170 + DN(1,2)*clhs207 + clhs278);
            lhs(6,13)=clhs0*(DN(1,0)*clhs118 + DN(1,1)*clhs174 + DN(1,2)*clhs209 + clhs279);
            lhs(6,14)=clhs0*(DN(1,0)*clhs124 + DN(1,1)*clhs180 + DN(1,2)*clhs211 + clhs266 + clhs281);
            lhs(6,15)=DN(3,2)*clhs228 - clhs282 + clhs283*clhs75;
            lhs(7,0)=clhs0*(DN(1,0)*clhs216 + clhs73);
            lhs(7,1)=clhs0*(DN(1,1)*clhs216 + clhs152);
            lhs(7,2)=clhs0*(DN(1,2)*clhs216 + clhs196);
            lhs(7,3)=clhs219;
            lhs(7,4)=DN(1,0)*clhs284;
            lhs(7,5)=DN(1,1)*clhs284;
            lhs(7,6)=DN(1,2)*clhs284;
            lhs(7,7)=clhs21*(clhs229 + clhs253 + clhs270);
            lhs(7,8)=clhs0*(DN(1,0)*clhs220 + clhs243);
            lhs(7,9)=clhs0*(DN(1,1)*clhs220 + clhs262);
            lhs(7,10)=clhs0*(DN(1,2)*clhs220 + clhs277);
            lhs(7,11)=clhs285;
            lhs(7,12)=clhs0*(DN(1,0)*clhs222 + clhs251);
            lhs(7,13)=clhs0*(DN(1,1)*clhs222 + clhs269);
            lhs(7,14)=clhs0*(DN(1,2)*clhs222 + clhs283);
            lhs(7,15)=clhs286;
            lhs(8,0)=clhs0*(DN(2,0)*clhs1 + DN(2,1)*clhs3 + DN(2,2)*clhs5 + clhs290 + clhs83);
            lhs(8,1)=clhs0*(DN(2,0)*clhs30 + DN(2,1)*clhs32 + DN(2,2)*clhs35 + clhs156);
            lhs(8,2)=clhs0*(DN(2,0)*clhs38 + DN(2,1)*clhs40 + DN(2,2)*clhs42 + clhs199);
            lhs(8,3)=DN(0,0)*clhs291 + clhs100*clhs75 - clhs101;
            lhs(8,4)=clhs0*(DN(2,0)*clhs48 + DN(2,1)*clhs50 + DN(2,2)*clhs52 + clhs238 + clhs292);
            lhs(8,5)=clhs0*(DN(2,0)*clhs60 + DN(2,1)*clhs62 + DN(2,2)*clhs65 + clhs256);
            lhs(8,6)=clhs0*(DN(2,0)*clhs67 + DN(2,1)*clhs69 + DN(2,2)*clhs71 + clhs272);
            lhs(8,7)=DN(1,0)*clhs291 + clhs242*clhs75 - clhs243;
            lhs(8,8)=clhs0*(DN(2,0)*clhs76 + DN(2,1)*clhs78 + DN(2,2)*clhs80 + clhs11*clhs293 + clhs294);
            lhs(8,9)=clhs0*(DN(2,0)*clhs87 + DN(2,1)*clhs89 + DN(2,2)*clhs92 + clhs296);
            lhs(8,10)=clhs0*(DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs297);
            lhs(8,11)=DN(2,0)*clhs298;
            lhs(8,12)=clhs0*(DN(2,0)*clhs102 + DN(2,1)*clhs104 + DN(2,2)*clhs106 + clhs301 + clhs302);
            lhs(8,13)=clhs0*(DN(2,0)*clhs113 + DN(2,1)*clhs115 + DN(2,2)*clhs118 + clhs303);
            lhs(8,14)=clhs0*(DN(2,0)*clhs120 + DN(2,1)*clhs122 + DN(2,2)*clhs124 + clhs304);
            lhs(8,15)=DN(3,0)*clhs291 - clhs305 + clhs306*clhs75;
            lhs(9,0)=clhs0*(DN(2,0)*clhs3 + DN(2,1)*clhs128 + DN(2,2)*clhs129 + clhs93);
            lhs(9,1)=clhs0*(DN(2,0)*clhs32 + DN(2,1)*clhs130 + DN(2,2)*clhs132 + clhs161 + clhs307);
            lhs(9,2)=clhs0*(DN(2,0)*clhs40 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs201);
            lhs(9,3)=DN(0,1)*clhs291 + clhs167*clhs75 - clhs168;
            lhs(9,4)=clhs0*(DN(2,0)*clhs50 + DN(2,1)*clhs139 + DN(2,2)*clhs140 + clhs240);
            lhs(9,5)=clhs0*(DN(2,0)*clhs62 + DN(2,1)*clhs142 + DN(2,2)*clhs144 + clhs258 + clhs308);
            lhs(9,6)=clhs0*(DN(2,0)*clhs69 + DN(2,1)*clhs148 + DN(2,2)*clhs150 + clhs273);
            lhs(9,7)=DN(1,1)*clhs291 + clhs261*clhs75 - clhs262;
            lhs(9,8)=clhs0*(DN(2,0)*clhs78 + DN(2,1)*clhs154 + DN(2,2)*clhs155 + clhs296);
            lhs(9,9)=clhs0*(DN(2,0)*clhs89 + DN(2,1)*clhs157 + DN(2,2)*clhs159 + clhs11*clhs309 + clhs294);
            lhs(9,10)=clhs0*(DN(2,0)*clhs96 + DN(2,1)*clhs163 + DN(2,2)*clhs165 + clhs311);
            lhs(9,11)=DN(2,1)*clhs298;
            lhs(9,12)=clhs0*(DN(2,0)*clhs104 + DN(2,1)*clhs169 + DN(2,2)*clhs170 + clhs312);
            lhs(9,13)=clhs0*(DN(2,0)*clhs115 + DN(2,1)*clhs172 + DN(2,2)*clhs174 + clhs314 + clhs315);
            lhs(9,14)=clhs0*(DN(2,0)*clhs122 + DN(2,1)*clhs178 + DN(2,2)*clhs180 + clhs316);
            lhs(9,15)=DN(3,1)*clhs291 - clhs317 + clhs318*clhs75;
            lhs(10,0)=clhs0*(DN(2,0)*clhs5 + DN(2,1)*clhs129 + DN(2,2)*clhs184 + clhs99);
            lhs(10,1)=clhs0*(DN(2,0)*clhs35 + DN(2,1)*clhs132 + DN(2,2)*clhs185 + clhs166);
            lhs(10,2)=clhs0*(DN(2,0)*clhs42 + DN(2,1)*clhs136 + DN(2,2)*clhs186 + clhs204 + clhs307);
            lhs(10,3)=DN(0,2)*clhs291 + clhs205*clhs75 - clhs206;
            lhs(10,4)=clhs0*(DN(2,0)*clhs52 + DN(2,1)*clhs140 + DN(2,2)*clhs188 + clhs241);
            lhs(10,5)=clhs0*(DN(2,0)*clhs65 + DN(2,1)*clhs144 + DN(2,2)*clhs191 + clhs260);
            lhs(10,6)=clhs0*(DN(2,0)*clhs71 + DN(2,1)*clhs150 + DN(2,2)*clhs193 + clhs275 + clhs308);
            lhs(10,7)=DN(1,2)*clhs291 + clhs276*clhs75 - clhs277;
            lhs(10,8)=clhs0*(DN(2,0)*clhs80 + DN(2,1)*clhs155 + DN(2,2)*clhs198 + clhs297);
            lhs(10,9)=clhs0*(DN(2,0)*clhs92 + DN(2,1)*clhs159 + DN(2,2)*clhs200 + clhs311);
            lhs(10,10)=clhs0*(DN(2,0)*clhs98 + DN(2,1)*clhs165 + DN(2,2)*clhs202 + clhs11*clhs319 + clhs294);
            lhs(10,11)=DN(2,2)*clhs298;
            lhs(10,12)=clhs0*(DN(2,0)*clhs106 + DN(2,1)*clhs170 + DN(2,2)*clhs207 + clhs321);
            lhs(10,13)=clhs0*(DN(2,0)*clhs118 + DN(2,1)*clhs174 + DN(2,2)*clhs209 + clhs322);
            lhs(10,14)=clhs0*(DN(2,0)*clhs124 + DN(2,1)*clhs180 + DN(2,2)*clhs211 + clhs315 + clhs324);
            lhs(10,15)=DN(3,2)*clhs291 - clhs325 + clhs326*clhs75;
            lhs(11,0)=clhs0*(DN(2,0)*clhs216 + clhs100);
            lhs(11,1)=clhs0*(DN(2,1)*clhs216 + clhs167);
            lhs(11,2)=clhs0*(DN(2,2)*clhs216 + clhs205);
            lhs(11,3)=clhs221;
            lhs(11,4)=clhs0*(DN(2,0)*clhs218 + clhs242);
            lhs(11,5)=clhs0*(DN(2,1)*clhs218 + clhs261);
            lhs(11,6)=clhs0*(DN(2,2)*clhs218 + clhs276);
            lhs(11,7)=clhs285;
            lhs(11,8)=DN(2,0)*clhs327;
            lhs(11,9)=DN(2,1)*clhs327;
            lhs(11,10)=DN(2,2)*clhs327;
            lhs(11,11)=clhs21*(clhs293 + clhs309 + clhs319);
            lhs(11,12)=clhs0*(DN(2,0)*clhs222 + clhs306);
            lhs(11,13)=clhs0*(DN(2,1)*clhs222 + clhs318);
            lhs(11,14)=clhs0*(DN(2,2)*clhs222 + clhs326);
            lhs(11,15)=clhs328;
            lhs(12,0)=clhs0*(DN(3,0)*clhs1 + DN(3,1)*clhs3 + DN(3,2)*clhs5 + clhs109 + clhs332);
            lhs(12,1)=clhs0*(DN(3,0)*clhs30 + DN(3,1)*clhs32 + DN(3,2)*clhs35 + clhs171);
            lhs(12,2)=clhs0*(DN(3,0)*clhs38 + DN(3,1)*clhs40 + DN(3,2)*clhs42 + clhs208);
            lhs(12,3)=DN(0,0)*clhs333 + clhs126*clhs75 - clhs127;
            lhs(12,4)=clhs0*(DN(3,0)*clhs48 + DN(3,1)*clhs50 + DN(3,2)*clhs52 + clhs246 + clhs334);
            lhs(12,5)=clhs0*(DN(3,0)*clhs60 + DN(3,1)*clhs62 + DN(3,2)*clhs65 + clhs263);
            lhs(12,6)=clhs0*(DN(3,0)*clhs67 + DN(3,1)*clhs69 + DN(3,2)*clhs71 + clhs278);
            lhs(12,7)=DN(1,0)*clhs333 + clhs250*clhs75 - clhs251;
            lhs(12,8)=clhs0*(DN(3,0)*clhs76 + DN(3,1)*clhs78 + DN(3,2)*clhs80 + clhs301 + clhs335);
            lhs(12,9)=clhs0*(DN(3,0)*clhs87 + DN(3,1)*clhs89 + DN(3,2)*clhs92 + clhs312);
            lhs(12,10)=clhs0*(DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs321);
            lhs(12,11)=DN(2,0)*clhs333 + clhs305*clhs75 - clhs306;
            lhs(12,12)=clhs0*(DN(3,0)*clhs102 + DN(3,1)*clhs104 + DN(3,2)*clhs106 + clhs11*clhs336 + clhs337);
            lhs(12,13)=clhs0*(DN(3,0)*clhs113 + DN(3,1)*clhs115 + DN(3,2)*clhs118 + clhs339);
            lhs(12,14)=clhs0*(DN(3,0)*clhs120 + DN(3,1)*clhs122 + DN(3,2)*clhs124 + clhs340);
            lhs(12,15)=DN(3,0)*clhs341;
            lhs(13,0)=clhs0*(DN(3,0)*clhs3 + DN(3,1)*clhs128 + DN(3,2)*clhs129 + clhs119);
            lhs(13,1)=clhs0*(DN(3,0)*clhs32 + DN(3,1)*clhs130 + DN(3,2)*clhs132 + clhs176 + clhs342);
            lhs(13,2)=clhs0*(DN(3,0)*clhs40 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs210);
            lhs(13,3)=DN(0,1)*clhs333 + clhs182*clhs75 - clhs183;
            lhs(13,4)=clhs0*(DN(3,0)*clhs50 + DN(3,1)*clhs139 + DN(3,2)*clhs140 + clhs248);
            lhs(13,5)=clhs0*(DN(3,0)*clhs62 + DN(3,1)*clhs142 + DN(3,2)*clhs144 + clhs265 + clhs343);
            lhs(13,6)=clhs0*(DN(3,0)*clhs69 + DN(3,1)*clhs148 + DN(3,2)*clhs150 + clhs279);
            lhs(13,7)=DN(1,1)*clhs333 + clhs268*clhs75 - clhs269;
            lhs(13,8)=clhs0*(DN(3,0)*clhs78 + DN(3,1)*clhs154 + DN(3,2)*clhs155 + clhs303);
            lhs(13,9)=clhs0*(DN(3,0)*clhs89 + DN(3,1)*clhs157 + DN(3,2)*clhs159 + clhs314 + clhs344);
            lhs(13,10)=clhs0*(DN(3,0)*clhs96 + DN(3,1)*clhs163 + DN(3,2)*clhs165 + clhs322);
            lhs(13,11)=DN(2,1)*clhs333 + clhs317*clhs75 - clhs318;
            lhs(13,12)=clhs0*(DN(3,0)*clhs104 + DN(3,1)*clhs169 + DN(3,2)*clhs170 + clhs339);
            lhs(13,13)=clhs0*(DN(3,0)*clhs115 + DN(3,1)*clhs172 + DN(3,2)*clhs174 + clhs11*clhs345 + clhs337);
            lhs(13,14)=clhs0*(DN(3,0)*clhs122 + DN(3,1)*clhs178 + DN(3,2)*clhs180 + clhs346);
            lhs(13,15)=DN(3,1)*clhs341;
            lhs(14,0)=clhs0*(DN(3,0)*clhs5 + DN(3,1)*clhs129 + DN(3,2)*clhs184 + clhs125);
            lhs(14,1)=clhs0*(DN(3,0)*clhs35 + DN(3,1)*clhs132 + DN(3,2)*clhs185 + clhs181);
            lhs(14,2)=clhs0*(DN(3,0)*clhs42 + DN(3,1)*clhs136 + DN(3,2)*clhs186 + clhs213 + clhs342);
            lhs(14,3)=DN(0,2)*clhs333 + clhs214*clhs75 - clhs215;
            lhs(14,4)=clhs0*(DN(3,0)*clhs52 + DN(3,1)*clhs140 + DN(3,2)*clhs188 + clhs249);
            lhs(14,5)=clhs0*(DN(3,0)*clhs65 + DN(3,1)*clhs144 + DN(3,2)*clhs191 + clhs267);
            lhs(14,6)=clhs0*(DN(3,0)*clhs71 + DN(3,1)*clhs150 + DN(3,2)*clhs193 + clhs281 + clhs343);
            lhs(14,7)=DN(1,2)*clhs333 + clhs282*clhs75 - clhs283;
            lhs(14,8)=clhs0*(DN(3,0)*clhs80 + DN(3,1)*clhs155 + DN(3,2)*clhs198 + clhs304);
            lhs(14,9)=clhs0*(DN(3,0)*clhs92 + DN(3,1)*clhs159 + DN(3,2)*clhs200 + clhs316);
            lhs(14,10)=clhs0*(DN(3,0)*clhs98 + DN(3,1)*clhs165 + DN(3,2)*clhs202 + clhs324 + clhs344);
            lhs(14,11)=DN(2,2)*clhs333 + clhs325*clhs75 - clhs326;
            lhs(14,12)=clhs0*(DN(3,0)*clhs106 + DN(3,1)*clhs170 + DN(3,2)*clhs207 + clhs340);
            lhs(14,13)=clhs0*(DN(3,0)*clhs118 + DN(3,1)*clhs174 + DN(3,2)*clhs209 + clhs346);
            lhs(14,14)=clhs0*(DN(3,0)*clhs124 + DN(3,1)*clhs180 + DN(3,2)*clhs211 + clhs11*clhs347 + clhs337);
            lhs(14,15)=DN(3,2)*clhs341;
            lhs(15,0)=clhs0*(DN(3,0)*clhs216 + clhs126);
            lhs(15,1)=clhs0*(DN(3,1)*clhs216 + clhs182);
            lhs(15,2)=clhs0*(DN(3,2)*clhs216 + clhs214);
            lhs(15,3)=clhs223;
            lhs(15,4)=clhs0*(DN(3,0)*clhs218 + clhs250);
            lhs(15,5)=clhs0*(DN(3,1)*clhs218 + clhs268);
            lhs(15,6)=clhs0*(DN(3,2)*clhs218 + clhs282);
            lhs(15,7)=clhs286;
            lhs(15,8)=clhs0*(DN(3,0)*clhs220 + clhs305);
            lhs(15,9)=clhs0*(DN(3,1)*clhs220 + clhs317);
            lhs(15,10)=clhs0*(DN(3,2)*clhs220 + clhs325);
            lhs(15,11)=clhs328;
            lhs(15,12)=DN(3,0)*clhs348;
            lhs(15,13)=DN(3,1)*clhs348;
            lhs(15,14)=DN(3,2)*clhs348;
            lhs(15,15)=clhs21*(clhs336 + clhs345 + clhs347);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,3,2> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             1.0/(max_spectral_radius + 1);
const double crhs2 =             rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)));
const double crhs3 =             -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1)) + volume_error_ratio;
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 =             rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs7 =             crhs3*(crhs6*h/stab_c1 + mu);
const double crhs8 =             v(0,0) - vn(0,0);
const double crhs9 =             crhs1*crhs8 + vn(0,0);
const double crhs10 =             v(1,0) - vn(1,0);
const double crhs11 =             crhs1*crhs10 + vn(1,0);
const double crhs12 =             v(2,0) - vn(2,0);
const double crhs13 =             crhs1*crhs12 + vn(2,0);
const double crhs14 =             rho*(crhs4*(DN(0,0)*crhs9 + DN(1,0)*crhs11 + DN(2,0)*crhs13) + crhs5*(DN(0,1)*crhs9 + DN(1,1)*crhs11 + DN(2,1)*crhs13));
const double crhs15 =             -acceleration_alpha_method(0,0);
const double crhs16 =             1.0/dt;
const double crhs17 =             0.5*max_spectral_radius;
const double crhs18 =             -crhs1;
const double crhs19 =             crhs18 + 0.5;
const double crhs20 =             1.0/(-crhs1*(crhs17 - 1.5) + crhs19);
const double crhs21 =             crhs16*crhs20;
const double crhs22 =             crhs1*(1.5 - crhs17);
const double crhs23 =             crhs20*(crhs1 - crhs22 + 0.5);
const double crhs24 =             0.5*crhs1;
const double crhs25 =             crhs24*(max_spectral_radius - 3);
const double crhs26 =             -acceleration_alpha_method(1,0);
const double crhs27 =             -acceleration_alpha_method(2,0);
const double crhs28 =             N[0]*(acceleration_alpha_method(0,0) - crhs25*(-acceleration_alpha_method(0,0)*crhs23 + crhs15 + crhs21*crhs8)) + N[1]*(acceleration_alpha_method(1,0) - crhs25*(-acceleration_alpha_method(1,0)*crhs23 + crhs10*crhs21 + crhs26)) + N[2]*(acceleration_alpha_method(2,0) - crhs25*(-acceleration_alpha_method(2,0)*crhs23 + crhs12*crhs21 + crhs27));
const double crhs29 =             N[0]*rho;
const double crhs30 =             1.0/(crhs19 + crhs22);
const double crhs31 =             crhs16*crhs30;
const double crhs32 =             crhs30*(crhs18 + crhs22 - 0.5);
const double crhs33 =             crhs24*(3 - max_spectral_radius);
const double crhs34 =             1.0/(crhs16*dyn_tau*rho + crhs6/h + mu*stab_c1/pow(h, 2));
const double crhs35 =             crhs34*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs14 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs33*(acceleration_alpha_method(0,0)*crhs32 + crhs15 + crhs31*crhs8)) + N[1]*(acceleration_alpha_method(1,0) + crhs33*(acceleration_alpha_method(1,0)*crhs32 + crhs10*crhs31 + crhs26)) + N[2]*(acceleration_alpha_method(2,0) + crhs33*(acceleration_alpha_method(2,0)*crhs32 + crhs12*crhs31 + crhs27))));
const double crhs36 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs37 =             crhs29*crhs36;
const double crhs38 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs39 =             rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)));
const double crhs40 =             v(0,1) - vn(0,1);
const double crhs41 =             crhs1*crhs40 + vn(0,1);
const double crhs42 =             v(1,1) - vn(1,1);
const double crhs43 =             crhs1*crhs42 + vn(1,1);
const double crhs44 =             v(2,1) - vn(2,1);
const double crhs45 =             crhs1*crhs44 + vn(2,1);
const double crhs46 =             rho*(crhs4*(DN(0,0)*crhs41 + DN(1,0)*crhs43 + DN(2,0)*crhs45) + crhs5*(DN(0,1)*crhs41 + DN(1,1)*crhs43 + DN(2,1)*crhs45));
const double crhs47 =             -acceleration_alpha_method(0,1);
const double crhs48 =             -acceleration_alpha_method(1,1);
const double crhs49 =             -acceleration_alpha_method(2,1);
const double crhs50 =             N[0]*(acceleration_alpha_method(0,1) - crhs25*(-acceleration_alpha_method(0,1)*crhs23 + crhs21*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) - crhs25*(-acceleration_alpha_method(1,1)*crhs23 + crhs21*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) - crhs25*(-acceleration_alpha_method(2,1)*crhs23 + crhs21*crhs44 + crhs49));
const double crhs51 =             crhs34*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs39 + crhs46 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs33*(acceleration_alpha_method(0,1)*crhs32 + crhs31*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) + crhs33*(acceleration_alpha_method(1,1)*crhs32 + crhs31*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) + crhs33*(acceleration_alpha_method(2,1)*crhs32 + crhs31*crhs44 + crhs49))));
const double crhs52 =             N[1]*rho;
const double crhs53 =             crhs36*crhs52;
const double crhs54 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs55 =             N[2]*rho;
const double crhs56 =             crhs36*crhs55;
const double crhs57 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
            rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs7 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs14 + N[0]*crhs2 - crhs28*crhs29 - crhs35*crhs37 - crhs35*crhs38;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs7 - DN(0,1)*stress[1] + N[0]*crhs39 - N[0]*crhs46 - crhs29*crhs50 - crhs37*crhs51 - crhs38*crhs51;
            rhs[2]=-DN(0,0)*crhs35 - DN(0,1)*crhs51 + N[0]*crhs3;
            rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs7 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs14 + N[1]*crhs2 - crhs28*crhs52 - crhs35*crhs53 - crhs35*crhs54;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs7 - DN(1,1)*stress[1] + N[1]*crhs39 - N[1]*crhs46 - crhs50*crhs52 - crhs51*crhs53 - crhs51*crhs54;
            rhs[5]=-DN(1,0)*crhs35 - DN(1,1)*crhs51 + N[1]*crhs3;
            rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs7 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs14 + N[2]*crhs2 - crhs28*crhs55 - crhs35*crhs56 - crhs35*crhs57;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs7 - DN(2,1)*stress[1] + N[2]*crhs39 - N[2]*crhs46 - crhs50*crhs55 - crhs51*crhs56 - crhs51*crhs57;
            rhs[8]=-DN(2,0)*crhs35 - DN(2,1)*crhs51 + N[2]*crhs3;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_f=1/(1+max_spectral_radius);
    const double dyn_tau = rData.DynamicTau;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             1.0/(max_spectral_radius + 1);
const double crhs2 =             rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs1*(f(3,0) - fn(3,0)) + fn(3,0)));
const double crhs3 =             -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(0,2)*v(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(1,2)*v(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(2,2)*v(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,0)*v(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,1)*v(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + DN(3,2)*v(3,2)) + volume_error_ratio;
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs6 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs7 =             rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2) + pow(crhs6, 2));
const double crhs8 =             crhs3*(crhs7*h/stab_c1 + mu);
const double crhs9 =             v(0,0) - vn(0,0);
const double crhs10 =             crhs1*crhs9 + vn(0,0);
const double crhs11 =             v(1,0) - vn(1,0);
const double crhs12 =             crhs1*crhs11 + vn(1,0);
const double crhs13 =             v(2,0) - vn(2,0);
const double crhs14 =             crhs1*crhs13 + vn(2,0);
const double crhs15 =             v(3,0) - vn(3,0);
const double crhs16 =             crhs1*crhs15 + vn(3,0);
const double crhs17 =             rho*(crhs4*(DN(0,0)*crhs10 + DN(1,0)*crhs12 + DN(2,0)*crhs14 + DN(3,0)*crhs16) + crhs5*(DN(0,1)*crhs10 + DN(1,1)*crhs12 + DN(2,1)*crhs14 + DN(3,1)*crhs16) + crhs6*(DN(0,2)*crhs10 + DN(1,2)*crhs12 + DN(2,2)*crhs14 + DN(3,2)*crhs16));
const double crhs18 =             -acceleration_alpha_method(0,0);
const double crhs19 =             1.0/dt;
const double crhs20 =             0.5*max_spectral_radius;
const double crhs21 =             -crhs1;
const double crhs22 =             crhs21 + 0.5;
const double crhs23 =             1.0/(-crhs1*(crhs20 - 1.5) + crhs22);
const double crhs24 =             crhs19*crhs23;
const double crhs25 =             crhs1*(1.5 - crhs20);
const double crhs26 =             crhs23*(crhs1 - crhs25 + 0.5);
const double crhs27 =             0.5*crhs1;
const double crhs28 =             crhs27*(max_spectral_radius - 3);
const double crhs29 =             -acceleration_alpha_method(1,0);
const double crhs30 =             -acceleration_alpha_method(2,0);
const double crhs31 =             -acceleration_alpha_method(3,0);
const double crhs32 =             N[0]*(acceleration_alpha_method(0,0) - crhs28*(-acceleration_alpha_method(0,0)*crhs26 + crhs18 + crhs24*crhs9)) + N[1]*(acceleration_alpha_method(1,0) - crhs28*(-acceleration_alpha_method(1,0)*crhs26 + crhs11*crhs24 + crhs29)) + N[2]*(acceleration_alpha_method(2,0) - crhs28*(-acceleration_alpha_method(2,0)*crhs26 + crhs13*crhs24 + crhs30)) + N[3]*(acceleration_alpha_method(3,0) - crhs28*(-acceleration_alpha_method(3,0)*crhs26 + crhs15*crhs24 + crhs31));
const double crhs33 =             N[0]*rho;
const double crhs34 =             1.0/(crhs22 + crhs25);
const double crhs35 =             crhs19*crhs34;
const double crhs36 =             crhs34*(crhs21 + crhs25 - 0.5);
const double crhs37 =             crhs27*(3 - max_spectral_radius);
const double crhs38 =             1.0/(crhs19*dyn_tau*rho + crhs7/h + mu*stab_c1/pow(h, 2));
const double crhs39 =             crhs38*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs17 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs37*(acceleration_alpha_method(0,0)*crhs36 + crhs18 + crhs35*crhs9)) + N[1]*(acceleration_alpha_method(1,0) + crhs37*(acceleration_alpha_method(1,0)*crhs36 + crhs11*crhs35 + crhs29)) + N[2]*(acceleration_alpha_method(2,0) + crhs37*(acceleration_alpha_method(2,0)*crhs36 + crhs13*crhs35 + crhs30)) + N[3]*(acceleration_alpha_method(3,0) + crhs37*(acceleration_alpha_method(3,0)*crhs36 + crhs15*crhs35 + crhs31))));
const double crhs40 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs41 =             crhs33*crhs40;
const double crhs42 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6);
const double crhs43 =             rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs1*(f(3,1) - fn(3,1)) + fn(3,1)));
const double crhs44 =             v(0,1) - vn(0,1);
const double crhs45 =             crhs1*crhs44 + vn(0,1);
const double crhs46 =             v(1,1) - vn(1,1);
const double crhs47 =             crhs1*crhs46 + vn(1,1);
const double crhs48 =             v(2,1) - vn(2,1);
const double crhs49 =             crhs1*crhs48 + vn(2,1);
const double crhs50 =             v(3,1) - vn(3,1);
const double crhs51 =             crhs1*crhs50 + vn(3,1);
const double crhs52 =             rho*(crhs4*(DN(0,0)*crhs45 + DN(1,0)*crhs47 + DN(2,0)*crhs49 + DN(3,0)*crhs51) + crhs5*(DN(0,1)*crhs45 + DN(1,1)*crhs47 + DN(2,1)*crhs49 + DN(3,1)*crhs51) + crhs6*(DN(0,2)*crhs45 + DN(1,2)*crhs47 + DN(2,2)*crhs49 + DN(3,2)*crhs51));
const double crhs53 =             -acceleration_alpha_method(0,1);
const double crhs54 =             -acceleration_alpha_method(1,1);
const double crhs55 =             -acceleration_alpha_method(2,1);
const double crhs56 =             -acceleration_alpha_method(3,1);
const double crhs57 =             N[0]*(acceleration_alpha_method(0,1) - crhs28*(-acceleration_alpha_method(0,1)*crhs26 + crhs24*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) - crhs28*(-acceleration_alpha_method(1,1)*crhs26 + crhs24*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) - crhs28*(-acceleration_alpha_method(2,1)*crhs26 + crhs24*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) - crhs28*(-acceleration_alpha_method(3,1)*crhs26 + crhs24*crhs50 + crhs56));
const double crhs58 =             crhs38*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs43 + crhs52 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs37*(acceleration_alpha_method(0,1)*crhs36 + crhs35*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) + crhs37*(acceleration_alpha_method(1,1)*crhs36 + crhs35*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) + crhs37*(acceleration_alpha_method(2,1)*crhs36 + crhs35*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) + crhs37*(acceleration_alpha_method(3,1)*crhs36 + crhs35*crhs50 + crhs56))));
const double crhs59 =             rho*(N[0]*(crhs1*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs1*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs1*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs1*(f(3,2) - fn(3,2)) + fn(3,2)));
const double crhs60 =             v(0,2) - vn(0,2);
const double crhs61 =             crhs1*crhs60 + vn(0,2);
const double crhs62 =             v(1,2) - vn(1,2);
const double crhs63 =             crhs1*crhs62 + vn(1,2);
const double crhs64 =             v(2,2) - vn(2,2);
const double crhs65 =             crhs1*crhs64 + vn(2,2);
const double crhs66 =             v(3,2) - vn(3,2);
const double crhs67 =             crhs1*crhs66 + vn(3,2);
const double crhs68 =             rho*(crhs4*(DN(0,0)*crhs61 + DN(1,0)*crhs63 + DN(2,0)*crhs65 + DN(3,0)*crhs67) + crhs5*(DN(0,1)*crhs61 + DN(1,1)*crhs63 + DN(2,1)*crhs65 + DN(3,1)*crhs67) + crhs6*(DN(0,2)*crhs61 + DN(1,2)*crhs63 + DN(2,2)*crhs65 + DN(3,2)*crhs67));
const double crhs69 =             -acceleration_alpha_method(0,2);
const double crhs70 =             -acceleration_alpha_method(1,2);
const double crhs71 =             -acceleration_alpha_method(2,2);
const double crhs72 =             -acceleration_alpha_method(3,2);
const double crhs73 =             N[0]*(acceleration_alpha_method(0,2) - crhs28*(-acceleration_alpha_method(0,2)*crhs26 + crhs24*crhs60 + crhs69)) + N[1]*(acceleration_alpha_method(1,2) - crhs28*(-acceleration_alpha_method(1,2)*crhs26 + crhs24*crhs62 + crhs70)) + N[2]*(acceleration_alpha_method(2,2) - crhs28*(-acceleration_alpha_method(2,2)*crhs26 + crhs24*crhs64 + crhs71)) + N[3]*(acceleration_alpha_method(3,2) - crhs28*(-acceleration_alpha_method(3,2)*crhs26 + crhs24*crhs66 + crhs72));
const double crhs74 =             crhs38*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs59 + crhs68 + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs37*(acceleration_alpha_method(0,2)*crhs36 + crhs35*crhs60 + crhs69)) + N[1]*(acceleration_alpha_method(1,2) + crhs37*(acceleration_alpha_method(1,2)*crhs36 + crhs35*crhs62 + crhs70)) + N[2]*(acceleration_alpha_method(2,2) + crhs37*(acceleration_alpha_method(2,2)*crhs36 + crhs35*crhs64 + crhs71)) + N[3]*(acceleration_alpha_method(3,2) + crhs37*(acceleration_alpha_method(3,2)*crhs36 + crhs35*crhs66 + crhs72))));
const double crhs75 =             N[1]*rho;
const double crhs76 =             crhs40*crhs75;
const double crhs77 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6);
const double crhs78 =             N[2]*rho;
const double crhs79 =             crhs40*crhs78;
const double crhs80 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6);
const double crhs81 =             N[3]*rho;
const double crhs82 =             crhs40*crhs81;
const double crhs83 =             rho*(DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6);
            rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs17 + N[0]*crhs2 - crhs32*crhs33 - crhs39*crhs41 - crhs39*crhs42;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 + DN(0,1)*crhs8 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs43 - N[0]*crhs52 - crhs33*crhs57 - crhs41*crhs58 - crhs42*crhs58;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 + DN(0,2)*crhs8 - DN(0,2)*stress[2] + N[0]*crhs59 - N[0]*crhs68 - crhs33*crhs73 - crhs41*crhs74 - crhs42*crhs74;
            rhs[3]=-DN(0,0)*crhs39 - DN(0,1)*crhs58 - DN(0,2)*crhs74 + N[0]*crhs3;
            rhs[4]=DN(1,0)*crhs0 + DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs17 + N[1]*crhs2 - crhs32*crhs75 - crhs39*crhs76 - crhs39*crhs77;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 + DN(1,1)*crhs8 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs43 - N[1]*crhs52 - crhs57*crhs75 - crhs58*crhs76 - crhs58*crhs77;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 + DN(1,2)*crhs8 - DN(1,2)*stress[2] + N[1]*crhs59 - N[1]*crhs68 - crhs73*crhs75 - crhs74*crhs76 - crhs74*crhs77;
            rhs[7]=-DN(1,0)*crhs39 - DN(1,1)*crhs58 - DN(1,2)*crhs74 + N[1]*crhs3;
            rhs[8]=DN(2,0)*crhs0 + DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs17 + N[2]*crhs2 - crhs32*crhs78 - crhs39*crhs79 - crhs39*crhs80;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 + DN(2,1)*crhs8 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs43 - N[2]*crhs52 - crhs57*crhs78 - crhs58*crhs79 - crhs58*crhs80;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 + DN(2,2)*crhs8 - DN(2,2)*stress[2] + N[2]*crhs59 - N[2]*crhs68 - crhs73*crhs78 - crhs74*crhs79 - crhs74*crhs80;
            rhs[11]=-DN(2,0)*crhs39 - DN(2,1)*crhs58 - DN(2,2)*crhs74 + N[2]*crhs3;
            rhs[12]=DN(3,0)*crhs0 + DN(3,0)*crhs8 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs17 + N[3]*crhs2 - crhs32*crhs81 - crhs39*crhs82 - crhs39*crhs83;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 + DN(3,1)*crhs8 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs43 - N[3]*crhs52 - crhs57*crhs81 - crhs58*crhs82 - crhs58*crhs83;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 + DN(3,2)*crhs8 - DN(3,2)*stress[2] + N[3]*crhs59 - N[3]*crhs68 - crhs73*crhs81 - crhs74*crhs82 - crhs74*crhs83;
            rhs[15]=-DN(3,0)*crhs39 - DN(3,1)*crhs58 - DN(3,2)*crhs74 + N[3]*crhs3;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p=rData.Pressure;

    const BoundedMatrix<double,3,2> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

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
const double cV4 =             cV3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV5 =             N[0]*cV4;
const double cV6 =             cV3*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV7 =             N[1]*cV4;
const double cV8 =             cV3*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV9 =             N[2]*cV4;
const double cV10 =             cV3*(DN(2,0)*cV0 + DN(2,1)*cV1);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV5 + DNenr(0,0)*cV6;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV5 + DNenr(1,0)*cV6;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV5 + DNenr(2,0)*cV6;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV5 + DNenr(0,1)*cV6;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV5 + DNenr(1,1)*cV6;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV5 + DNenr(2,1)*cV6;
            V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV7 + DNenr(0,0)*cV8;
            V(3,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV7 + DNenr(1,0)*cV8;
            V(3,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV7 + DNenr(2,0)*cV8;
            V(4,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV7 + DNenr(0,1)*cV8;
            V(4,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV7 + DNenr(1,1)*cV8;
            V(4,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV7 + DNenr(2,1)*cV8;
            V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV10 + DNenr(0,0)*cV9;
            V(6,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV9;
            V(6,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV9;
            V(7,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV9;
            V(7,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV9;
            V(7,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV9;
            V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             1.0/(max_spectral_radius + 1);
const double cH1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH3 =             1.0/dt;
const double cH4 =             0.5*cH3*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH5 =             1.0/(cH3*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2));
const double cH6 =             cH5*rho;
const double cH7 =             cH6*(DN(0,0)*cH1 + DN(0,1)*cH2 + N[0]*cH4);
const double cH8 =             cH6*(DN(1,0)*cH1 + DN(1,1)*cH2 + N[1]*cH4);
const double cH9 =             cH6*(DN(2,0)*cH1 + DN(2,1)*cH2 + N[2]*cH4);
            H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH7);
            H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH7);
            H(0,2)=cH5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH8);
            H(0,4)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH8);
            H(0,5)=cH5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH9);
            H(0,7)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH9);
            H(0,8)=cH5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH7);
            H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH7);
            H(1,2)=cH5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH8);
            H(1,4)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH8);
            H(1,5)=cH5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH9);
            H(1,7)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH9);
            H(1,8)=cH5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH7);
            H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH7);
            H(2,2)=cH5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH8);
            H(2,4)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH8);
            H(2,5)=cH5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH9);
            H(2,7)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH9);
            H(2,8)=cH5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


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


    const double crhs_ee0 =             1.0/(max_spectral_radius + 1);
const double crhs_ee1 =             -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1)) + volume_error_ratio;
const double crhs_ee2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee3 =             v(0,0) - vn(0,0);
const double crhs_ee4 =             crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 =             v(1,0) - vn(1,0);
const double crhs_ee6 =             crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 =             v(2,0) - vn(2,0);
const double crhs_ee8 =             crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee10 =             1.0/dt;
const double crhs_ee11 =             crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee12 =             1.0/(crhs_ee11 + 0.5);
const double crhs_ee13 =             crhs_ee10*crhs_ee12;
const double crhs_ee14 =             crhs_ee12*(crhs_ee11 - 0.5);
const double crhs_ee15 =             0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee16 =             1.0/(crhs_ee10*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee2, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee17 =             crhs_ee16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee15*(acceleration_alpha_method(0,0)*crhs_ee14 - acceleration_alpha_method(0,0) + crhs_ee13*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee15*(acceleration_alpha_method(1,0)*crhs_ee14 - acceleration_alpha_method(1,0) + crhs_ee13*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee15*(acceleration_alpha_method(2,0)*crhs_ee14 - acceleration_alpha_method(2,0) + crhs_ee13*crhs_ee7)) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8) + crhs_ee9*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8)));
const double crhs_ee18 =             v(0,1) - vn(0,1);
const double crhs_ee19 =             crhs_ee0*crhs_ee18 + vn(0,1);
const double crhs_ee20 =             v(1,1) - vn(1,1);
const double crhs_ee21 =             crhs_ee0*crhs_ee20 + vn(1,1);
const double crhs_ee22 =             v(2,1) - vn(2,1);
const double crhs_ee23 =             crhs_ee0*crhs_ee22 + vn(2,1);
const double crhs_ee24 =             crhs_ee16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee15*(acceleration_alpha_method(0,1)*crhs_ee14 - acceleration_alpha_method(0,1) + crhs_ee13*crhs_ee18)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee15*(acceleration_alpha_method(1,1)*crhs_ee14 - acceleration_alpha_method(1,1) + crhs_ee13*crhs_ee20)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee15*(acceleration_alpha_method(2,1)*crhs_ee14 - acceleration_alpha_method(2,1) + crhs_ee13*crhs_ee22)) + crhs_ee2*(DN(0,0)*crhs_ee19 + DN(1,0)*crhs_ee21 + DN(2,0)*crhs_ee23) + crhs_ee9*(DN(0,1)*crhs_ee19 + DN(1,1)*crhs_ee21 + DN(2,1)*crhs_ee23)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee17 - DNenr(0,1)*crhs_ee24 + Nenr[0]*crhs_ee1;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee17 - DNenr(1,1)*crhs_ee24 + Nenr[1]*crhs_ee1;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee17 - DNenr(2,1)*crhs_ee24 + Nenr[2]*crhs_ee1;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_f=1/(1-max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeError;

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
const double cV5 =             cV4*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV6 =             N[0]*cV5;
const double cV7 =             cV4*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV8 =             N[1]*cV5;
const double cV9 =             cV4*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV10 =             N[2]*cV5;
const double cV11 =             cV4*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV12 =             N[3]*cV5;
const double cV13 =             cV4*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV6 + DNenr(0,0)*cV7;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV6 + DNenr(1,0)*cV7;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV6 + DNenr(2,0)*cV7;
            V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV6 + DNenr(3,0)*cV7;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV6 + DNenr(0,1)*cV7;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV6 + DNenr(1,1)*cV7;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV6 + DNenr(2,1)*cV7;
            V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV6 + DNenr(3,1)*cV7;
            V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV6 + DNenr(0,2)*cV7;
            V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV6 + DNenr(1,2)*cV7;
            V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV6 + DNenr(2,2)*cV7;
            V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV6 + DNenr(3,2)*cV7;
            V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV8 + DNenr(0,0)*cV9;
            V(4,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV8 + DNenr(1,0)*cV9;
            V(4,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV8 + DNenr(2,0)*cV9;
            V(4,3)=-DN(1,0)*Nenr[3] + DNenr(3,0)*cV8 + DNenr(3,0)*cV9;
            V(5,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV8 + DNenr(0,1)*cV9;
            V(5,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV8 + DNenr(1,1)*cV9;
            V(5,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV8 + DNenr(2,1)*cV9;
            V(5,3)=-DN(1,1)*Nenr[3] + DNenr(3,1)*cV8 + DNenr(3,1)*cV9;
            V(6,0)=-DN(1,2)*Nenr[0] + DNenr(0,2)*cV8 + DNenr(0,2)*cV9;
            V(6,1)=-DN(1,2)*Nenr[1] + DNenr(1,2)*cV8 + DNenr(1,2)*cV9;
            V(6,2)=-DN(1,2)*Nenr[2] + DNenr(2,2)*cV8 + DNenr(2,2)*cV9;
            V(6,3)=-DN(1,2)*Nenr[3] + DNenr(3,2)*cV8 + DNenr(3,2)*cV9;
            V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV10 + DNenr(0,0)*cV11;
            V(8,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV11;
            V(8,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV11;
            V(8,3)=-DN(2,0)*Nenr[3] + DNenr(3,0)*cV10 + DNenr(3,0)*cV11;
            V(9,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV11;
            V(9,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV11;
            V(9,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV11;
            V(9,3)=-DN(2,1)*Nenr[3] + DNenr(3,1)*cV10 + DNenr(3,1)*cV11;
            V(10,0)=-DN(2,2)*Nenr[0] + DNenr(0,2)*cV10 + DNenr(0,2)*cV11;
            V(10,1)=-DN(2,2)*Nenr[1] + DNenr(1,2)*cV10 + DNenr(1,2)*cV11;
            V(10,2)=-DN(2,2)*Nenr[2] + DNenr(2,2)*cV10 + DNenr(2,2)*cV11;
            V(10,3)=-DN(2,2)*Nenr[3] + DNenr(3,2)*cV10 + DNenr(3,2)*cV11;
            V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] + DNenr(0,0)*cV12 + DNenr(0,0)*cV13;
            V(12,1)=-DN(3,0)*Nenr[1] + DNenr(1,0)*cV12 + DNenr(1,0)*cV13;
            V(12,2)=-DN(3,0)*Nenr[2] + DNenr(2,0)*cV12 + DNenr(2,0)*cV13;
            V(12,3)=-DN(3,0)*Nenr[3] + DNenr(3,0)*cV12 + DNenr(3,0)*cV13;
            V(13,0)=-DN(3,1)*Nenr[0] + DNenr(0,1)*cV12 + DNenr(0,1)*cV13;
            V(13,1)=-DN(3,1)*Nenr[1] + DNenr(1,1)*cV12 + DNenr(1,1)*cV13;
            V(13,2)=-DN(3,1)*Nenr[2] + DNenr(2,1)*cV12 + DNenr(2,1)*cV13;
            V(13,3)=-DN(3,1)*Nenr[3] + DNenr(3,1)*cV12 + DNenr(3,1)*cV13;
            V(14,0)=-DN(3,2)*Nenr[0] + DNenr(0,2)*cV12 + DNenr(0,2)*cV13;
            V(14,1)=-DN(3,2)*Nenr[1] + DNenr(1,2)*cV12 + DNenr(1,2)*cV13;
            V(14,2)=-DN(3,2)*Nenr[2] + DNenr(2,2)*cV12 + DNenr(2,2)*cV13;
            V(14,3)=-DN(3,2)*Nenr[3] + DNenr(3,2)*cV12 + DNenr(3,2)*cV13;
            V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             1.0/(max_spectral_radius + 1);
const double cH1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH3 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH4 =             1.0/dt;
const double cH5 =             0.5*cH4*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH6 =             1.0/(cH4*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 =             cH6*rho;
const double cH8 =             cH7*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 + N[0]*cH5);
const double cH9 =             cH7*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 + N[1]*cH5);
const double cH10 =             cH7*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 + N[2]*cH5);
const double cH11 =             cH7*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 + N[3]*cH5);
            H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH8);
            H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH8);
            H(0,2)=cH0*(DN(0,2)*Nenr[0] + DNenr(0,2)*cH8);
            H(0,3)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH9);
            H(0,5)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH9);
            H(0,6)=cH0*(DN(1,2)*Nenr[0] + DNenr(0,2)*cH9);
            H(0,7)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH10);
            H(0,9)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH10);
            H(0,10)=cH0*(DN(2,2)*Nenr[0] + DNenr(0,2)*cH10);
            H(0,11)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=cH0*(DN(3,0)*Nenr[0] + DNenr(0,0)*cH11);
            H(0,13)=cH0*(DN(3,1)*Nenr[0] + DNenr(0,1)*cH11);
            H(0,14)=cH0*(DN(3,2)*Nenr[0] + DNenr(0,2)*cH11);
            H(0,15)=cH6*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH8);
            H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH8);
            H(1,2)=cH0*(DN(0,2)*Nenr[1] + DNenr(1,2)*cH8);
            H(1,3)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH9);
            H(1,5)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH9);
            H(1,6)=cH0*(DN(1,2)*Nenr[1] + DNenr(1,2)*cH9);
            H(1,7)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH10);
            H(1,9)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH10);
            H(1,10)=cH0*(DN(2,2)*Nenr[1] + DNenr(1,2)*cH10);
            H(1,11)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=cH0*(DN(3,0)*Nenr[1] + DNenr(1,0)*cH11);
            H(1,13)=cH0*(DN(3,1)*Nenr[1] + DNenr(1,1)*cH11);
            H(1,14)=cH0*(DN(3,2)*Nenr[1] + DNenr(1,2)*cH11);
            H(1,15)=cH6*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH8);
            H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH8);
            H(2,2)=cH0*(DN(0,2)*Nenr[2] + DNenr(2,2)*cH8);
            H(2,3)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH9);
            H(2,5)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH9);
            H(2,6)=cH0*(DN(1,2)*Nenr[2] + DNenr(2,2)*cH9);
            H(2,7)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH10);
            H(2,9)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH10);
            H(2,10)=cH0*(DN(2,2)*Nenr[2] + DNenr(2,2)*cH10);
            H(2,11)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=cH0*(DN(3,0)*Nenr[2] + DNenr(2,0)*cH11);
            H(2,13)=cH0*(DN(3,1)*Nenr[2] + DNenr(2,1)*cH11);
            H(2,14)=cH0*(DN(3,2)*Nenr[2] + DNenr(2,2)*cH11);
            H(2,15)=cH6*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=cH0*(DN(0,0)*Nenr[3] + DNenr(3,0)*cH8);
            H(3,1)=cH0*(DN(0,1)*Nenr[3] + DNenr(3,1)*cH8);
            H(3,2)=cH0*(DN(0,2)*Nenr[3] + DNenr(3,2)*cH8);
            H(3,3)=cH6*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=cH0*(DN(1,0)*Nenr[3] + DNenr(3,0)*cH9);
            H(3,5)=cH0*(DN(1,1)*Nenr[3] + DNenr(3,1)*cH9);
            H(3,6)=cH0*(DN(1,2)*Nenr[3] + DNenr(3,2)*cH9);
            H(3,7)=cH6*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=cH0*(DN(2,0)*Nenr[3] + DNenr(3,0)*cH10);
            H(3,9)=cH0*(DN(2,1)*Nenr[3] + DNenr(3,1)*cH10);
            H(3,10)=cH0*(DN(2,2)*Nenr[3] + DNenr(3,2)*cH10);
            H(3,11)=cH6*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=cH0*(DN(3,0)*Nenr[3] + DNenr(3,0)*cH11);
            H(3,13)=cH0*(DN(3,1)*Nenr[3] + DNenr(3,1)*cH11);
            H(3,14)=cH0*(DN(3,2)*Nenr[3] + DNenr(3,2)*cH11);
            H(3,15)=cH6*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


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


    const double crhs_ee0 =             1.0/(max_spectral_radius + 1);
const double crhs_ee1 =             -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(0,2)*v(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(1,2)*v(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(2,2)*v(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,0)*v(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,1)*v(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + DN(3,2)*v(3,2)) + volume_error_ratio;
const double crhs_ee2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee3 =             v(0,0) - vn(0,0);
const double crhs_ee4 =             crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 =             v(1,0) - vn(1,0);
const double crhs_ee6 =             crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 =             v(2,0) - vn(2,0);
const double crhs_ee8 =             crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 =             v(3,0) - vn(3,0);
const double crhs_ee10 =             crhs_ee0*crhs_ee9 + vn(3,0);
const double crhs_ee11 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee12 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee13 =             1.0/dt;
const double crhs_ee14 =             crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee15 =             1.0/(crhs_ee14 + 0.5);
const double crhs_ee16 =             crhs_ee13*crhs_ee15;
const double crhs_ee17 =             crhs_ee15*(crhs_ee14 - 0.5);
const double crhs_ee18 =             0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee19 =             1.0/(crhs_ee13*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee11, 2) + pow(crhs_ee12, 2) + pow(crhs_ee2, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee20 =             crhs_ee19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs_ee0*(f(3,0) - fn(3,0)) + fn(3,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee18*(acceleration_alpha_method(0,0)*crhs_ee17 - acceleration_alpha_method(0,0) + crhs_ee16*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee18*(acceleration_alpha_method(1,0)*crhs_ee17 - acceleration_alpha_method(1,0) + crhs_ee16*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee18*(acceleration_alpha_method(2,0)*crhs_ee17 - acceleration_alpha_method(2,0) + crhs_ee16*crhs_ee7)) + N[3]*(acceleration_alpha_method(3,0) + crhs_ee18*(acceleration_alpha_method(3,0)*crhs_ee17 - acceleration_alpha_method(3,0) + crhs_ee16*crhs_ee9)) + crhs_ee11*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8 + DN(3,1)*crhs_ee10) + crhs_ee12*(DN(0,2)*crhs_ee4 + DN(1,2)*crhs_ee6 + DN(2,2)*crhs_ee8 + DN(3,2)*crhs_ee10) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8 + DN(3,0)*crhs_ee10)));
const double crhs_ee21 =             v(0,1) - vn(0,1);
const double crhs_ee22 =             crhs_ee0*crhs_ee21 + vn(0,1);
const double crhs_ee23 =             v(1,1) - vn(1,1);
const double crhs_ee24 =             crhs_ee0*crhs_ee23 + vn(1,1);
const double crhs_ee25 =             v(2,1) - vn(2,1);
const double crhs_ee26 =             crhs_ee0*crhs_ee25 + vn(2,1);
const double crhs_ee27 =             v(3,1) - vn(3,1);
const double crhs_ee28 =             crhs_ee0*crhs_ee27 + vn(3,1);
const double crhs_ee29 =             crhs_ee19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs_ee0*(f(3,1) - fn(3,1)) + fn(3,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee18*(acceleration_alpha_method(0,1)*crhs_ee17 - acceleration_alpha_method(0,1) + crhs_ee16*crhs_ee21)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee18*(acceleration_alpha_method(1,1)*crhs_ee17 - acceleration_alpha_method(1,1) + crhs_ee16*crhs_ee23)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee18*(acceleration_alpha_method(2,1)*crhs_ee17 - acceleration_alpha_method(2,1) + crhs_ee16*crhs_ee25)) + N[3]*(acceleration_alpha_method(3,1) + crhs_ee18*(acceleration_alpha_method(3,1)*crhs_ee17 - acceleration_alpha_method(3,1) + crhs_ee16*crhs_ee27)) + crhs_ee11*(DN(0,1)*crhs_ee22 + DN(1,1)*crhs_ee24 + DN(2,1)*crhs_ee26 + DN(3,1)*crhs_ee28) + crhs_ee12*(DN(0,2)*crhs_ee22 + DN(1,2)*crhs_ee24 + DN(2,2)*crhs_ee26 + DN(3,2)*crhs_ee28) + crhs_ee2*(DN(0,0)*crhs_ee22 + DN(1,0)*crhs_ee24 + DN(2,0)*crhs_ee26 + DN(3,0)*crhs_ee28)));
const double crhs_ee30 =             v(0,2) - vn(0,2);
const double crhs_ee31 =             crhs_ee0*crhs_ee30 + vn(0,2);
const double crhs_ee32 =             v(1,2) - vn(1,2);
const double crhs_ee33 =             crhs_ee0*crhs_ee32 + vn(1,2);
const double crhs_ee34 =             v(2,2) - vn(2,2);
const double crhs_ee35 =             crhs_ee0*crhs_ee34 + vn(2,2);
const double crhs_ee36 =             v(3,2) - vn(3,2);
const double crhs_ee37 =             crhs_ee0*crhs_ee36 + vn(3,2);
const double crhs_ee38 =             crhs_ee19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs_ee0*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs_ee0*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs_ee0*(f(3,2) - fn(3,2)) + fn(3,2))) + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs_ee18*(acceleration_alpha_method(0,2)*crhs_ee17 - acceleration_alpha_method(0,2) + crhs_ee16*crhs_ee30)) + N[1]*(acceleration_alpha_method(1,2) + crhs_ee18*(acceleration_alpha_method(1,2)*crhs_ee17 - acceleration_alpha_method(1,2) + crhs_ee16*crhs_ee32)) + N[2]*(acceleration_alpha_method(2,2) + crhs_ee18*(acceleration_alpha_method(2,2)*crhs_ee17 - acceleration_alpha_method(2,2) + crhs_ee16*crhs_ee34)) + N[3]*(acceleration_alpha_method(3,2) + crhs_ee18*(acceleration_alpha_method(3,2)*crhs_ee17 - acceleration_alpha_method(3,2) + crhs_ee16*crhs_ee36)) + crhs_ee11*(DN(0,1)*crhs_ee31 + DN(1,1)*crhs_ee33 + DN(2,1)*crhs_ee35 + DN(3,1)*crhs_ee37) + crhs_ee12*(DN(0,2)*crhs_ee31 + DN(1,2)*crhs_ee33 + DN(2,2)*crhs_ee35 + DN(3,2)*crhs_ee37) + crhs_ee2*(DN(0,0)*crhs_ee31 + DN(1,0)*crhs_ee33 + DN(2,0)*crhs_ee35 + DN(3,0)*crhs_ee37)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee20 - DNenr(0,1)*crhs_ee29 - DNenr(0,2)*crhs_ee38 + Nenr[0]*crhs_ee1;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee20 - DNenr(1,1)*crhs_ee29 - DNenr(1,2)*crhs_ee38 + Nenr[1]*crhs_ee1;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee20 - DNenr(2,1)*crhs_ee29 - DNenr(2,2)*crhs_ee38 + Nenr[2]*crhs_ee1;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee20 - DNenr(3,1)*crhs_ee29 - DNenr(3,2)*crhs_ee38 + Nenr[3]*crhs_ee1;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::PressureGradientStabilization(
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
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));

    const auto vmesh=rData.MeshVelocity;
    const auto vmeshn=rData.MeshVelocityOldStep;
    const BoundedMatrix<double,NumNodes,Dim> v_convection = (vn-vmeshn)+ alpha_m*((vn-vmeshn)-(v-vmesh));

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
void TwoFluidNavierStokesAlphaMethod<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = TwoFluidNavierStokes<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = TwoFluidNavierStokes<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>;
template class TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>;

} // namespace Kratos
