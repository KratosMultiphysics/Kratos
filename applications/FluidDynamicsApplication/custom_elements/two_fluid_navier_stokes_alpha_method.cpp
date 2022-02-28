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
    const double volume_error_ratio = rData.VolumeError;
    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 = pow(DN(0,0), 2);
const double clhs1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs3 = rho*stab_c2*sqrt(pow(clhs1, 2) + pow(clhs2, 2));
const double clhs4 = clhs3*h/stab_c1 + mu;
const double clhs5 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs6 = 1.0/(max_spectral_radius + 1);
const double clhs7 = DN(0,0)*clhs6;
const double clhs8 = C(0,2)*DN(0,0);
const double clhs9 = C(2,2)*DN(0,1) + clhs8;
const double clhs10 = DN(0,1)*clhs6;
const double clhs11 = DN(0,0)*clhs1 + DN(0,1)*clhs2;
const double clhs12 = clhs6*rho;
const double clhs13 = N[0]*clhs12;
const double clhs14 = pow(N[0], 2);
const double clhs15 = clhs12*volume_error_ratio;
const double clhs16 = N[0]*volume_error_ratio;
const double clhs17 = 0.5*max_spectral_radius - 1.5;
const double clhs18 = N[0]*clhs17;
const double clhs19 = 1.0/dt;
const double clhs20 = clhs6*(0.5*max_spectral_radius - 1.5);
const double clhs21 = clhs19/(-clhs20 - clhs6 + 0.5);
const double clhs22 = clhs11 + clhs16 - clhs18*clhs21;
const double clhs23 = clhs19*rho;
const double clhs24 = 1.0/(clhs23*dyn_tau + clhs3/h + mu*stab_c1/pow(h, 2));
const double clhs25 = clhs24*clhs6*pow(rho, 2);
const double clhs26 = clhs11*clhs25;
const double clhs27 = clhs23*clhs6/(clhs20 + clhs6 - 0.5);
const double clhs28 = clhs17*clhs27;
const double clhs29 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs30 = clhs25*clhs29;
const double clhs31 = N[0]*clhs30;
const double clhs32 = clhs11*clhs13 + clhs14*clhs15 + clhs14*clhs28 + clhs22*clhs26 + clhs22*clhs31;
const double clhs33 = C(0,1)*DN(0,1) + clhs8;
const double clhs34 = C(1,2)*DN(0,1);
const double clhs35 = C(2,2)*DN(0,0) + clhs34;
const double clhs36 = DN(0,0)*clhs4;
const double clhs37 = DN(0,1)*clhs36;
const double clhs38 = clhs24*rho;
const double clhs39 = clhs29*clhs38;
const double clhs40 = clhs11*clhs38;
const double clhs41 = N[0]*clhs39 - N[0] + clhs40;
const double clhs42 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs43 = C(0,2)*DN(1,0);
const double clhs44 = C(2,2)*DN(1,1) + clhs43;
const double clhs45 = DN(0,0)*DN(1,0);
const double clhs46 = N[1]*clhs12;
const double clhs47 = clhs18*clhs27;
const double clhs48 = N[1]*clhs47 + clhs16*clhs46;
const double clhs49 = clhs4*clhs45 + clhs48;
const double clhs50 = DN(1,0)*clhs1 + DN(1,1)*clhs2;
const double clhs51 = N[1]*volume_error_ratio;
const double clhs52 = clhs17*clhs21;
const double clhs53 = -N[1]*clhs52 + clhs50 + clhs51;
const double clhs54 = clhs13*clhs50 + clhs26*clhs53 + clhs31*clhs53;
const double clhs55 = C(0,1)*DN(1,1) + clhs43;
const double clhs56 = C(1,2)*DN(1,1);
const double clhs57 = C(2,2)*DN(1,0) + clhs56;
const double clhs58 = DN(1,1)*clhs36;
const double clhs59 = DN(0,0)*N[1];
const double clhs60 = DN(1,0)*N[0];
const double clhs61 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs62 = C(0,2)*DN(2,0);
const double clhs63 = C(2,2)*DN(2,1) + clhs62;
const double clhs64 = DN(0,0)*DN(2,0);
const double clhs65 = N[2]*clhs12;
const double clhs66 = N[2]*clhs47 + clhs16*clhs65;
const double clhs67 = clhs4*clhs64 + clhs66;
const double clhs68 = DN(2,0)*clhs1 + DN(2,1)*clhs2;
const double clhs69 = -N[2]*clhs52 + N[2]*volume_error_ratio + clhs68;
const double clhs70 = clhs13*clhs68 + clhs26*clhs69 + clhs31*clhs69;
const double clhs71 = C(0,1)*DN(2,1) + clhs62;
const double clhs72 = C(1,2)*DN(2,1);
const double clhs73 = C(2,2)*DN(2,0) + clhs72;
const double clhs74 = DN(2,1)*clhs36;
const double clhs75 = DN(0,0)*N[2];
const double clhs76 = DN(2,0)*N[0];
const double clhs77 = C(0,1)*DN(0,0) + clhs34;
const double clhs78 = pow(DN(0,1), 2);
const double clhs79 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs80 = C(0,1)*DN(1,0) + clhs56;
const double clhs81 = DN(0,1)*clhs4;
const double clhs82 = DN(1,0)*clhs81;
const double clhs83 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs84 = DN(0,1)*DN(1,1);
const double clhs85 = clhs4*clhs84 + clhs48;
const double clhs86 = DN(0,1)*N[1];
const double clhs87 = DN(1,1)*N[0];
const double clhs88 = C(0,1)*DN(2,0) + clhs72;
const double clhs89 = DN(2,0)*clhs81;
const double clhs90 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs91 = DN(0,1)*DN(2,1);
const double clhs92 = clhs4*clhs91 + clhs66;
const double clhs93 = DN(0,1)*N[2];
const double clhs94 = DN(2,1)*N[0];
const double clhs95 = clhs22*clhs38;
const double clhs96 = N[0] + clhs95;
const double clhs97 = clhs38*clhs53;
const double clhs98 = clhs24*(clhs45 + clhs84);
const double clhs99 = clhs38*clhs69;
const double clhs100 = clhs24*(clhs64 + clhs91);
const double clhs101 = DN(1,0)*clhs6;
const double clhs102 = DN(1,1)*clhs6;
const double clhs103 = clhs25*clhs50;
const double clhs104 = N[1]*clhs30;
const double clhs105 = clhs103*clhs22 + clhs104*clhs22 + clhs11*clhs46;
const double clhs106 = clhs38*clhs50;
const double clhs107 = pow(DN(1,0), 2);
const double clhs108 = pow(N[1], 2);
const double clhs109 = clhs103*clhs53 + clhs104*clhs53 + clhs108*clhs15 + clhs108*clhs28 + clhs46*clhs50;
const double clhs110 = DN(1,0)*clhs4;
const double clhs111 = DN(1,1)*clhs110;
const double clhs112 = N[1]*clhs39 - N[1] + clhs106;
const double clhs113 = DN(1,0)*DN(2,0);
const double clhs114 = N[1]*N[2]*clhs28 + clhs51*clhs65;
const double clhs115 = clhs113*clhs4 + clhs114;
const double clhs116 = clhs103*clhs69 + clhs104*clhs69 + clhs46*clhs68;
const double clhs117 = DN(2,1)*clhs110;
const double clhs118 = DN(1,0)*N[2];
const double clhs119 = DN(2,0)*N[1];
const double clhs120 = pow(DN(1,1), 2);
const double clhs121 = DN(2,0)*clhs4;
const double clhs122 = DN(1,1)*clhs121;
const double clhs123 = DN(1,1)*DN(2,1);
const double clhs124 = clhs114 + clhs123*clhs4;
const double clhs125 = DN(1,1)*N[2];
const double clhs126 = DN(2,1)*N[1];
const double clhs127 = N[1] + clhs97;
const double clhs128 = clhs24*(clhs113 + clhs123);
const double clhs129 = DN(2,0)*clhs6;
const double clhs130 = DN(2,1)*clhs6;
const double clhs131 = clhs25*clhs68;
const double clhs132 = N[2]*clhs30;
const double clhs133 = clhs11*clhs65 + clhs131*clhs22 + clhs132*clhs22;
const double clhs134 = clhs38*clhs68;
const double clhs135 = clhs131*clhs53 + clhs132*clhs53 + clhs50*clhs65;
const double clhs136 = pow(DN(2,0), 2);
const double clhs137 = pow(N[2], 2);
const double clhs138 = clhs131*clhs69 + clhs132*clhs69 + clhs137*clhs15 + clhs137*clhs28 + clhs65*clhs68;
const double clhs139 = DN(2,1)*clhs121;
const double clhs140 = N[2]*clhs39 - N[2] + clhs134;
const double clhs141 = pow(DN(2,1), 2);
const double clhs142 = N[2] + clhs99;
lhs(0,0)=clhs0*clhs4 + clhs10*clhs9 + clhs32 + clhs5*clhs7;
lhs(0,1)=clhs10*clhs35 + clhs33*clhs7 + clhs37;
lhs(0,2)=DN(0,0)*clhs41;
lhs(0,3)=clhs10*clhs44 + clhs42*clhs7 + clhs49 + clhs54;
lhs(0,4)=clhs10*clhs57 + clhs55*clhs7 + clhs58;
lhs(0,5)=DN(1,0)*clhs40 + clhs39*clhs60 - clhs59;
lhs(0,6)=clhs10*clhs63 + clhs61*clhs7 + clhs67 + clhs70;
lhs(0,7)=clhs10*clhs73 + clhs7*clhs71 + clhs74;
lhs(0,8)=DN(2,0)*clhs40 + clhs39*clhs76 - clhs75;
lhs(1,0)=clhs10*clhs77 + clhs37 + clhs7*clhs9;
lhs(1,1)=clhs10*clhs79 + clhs32 + clhs35*clhs7 + clhs4*clhs78;
lhs(1,2)=DN(0,1)*clhs41;
lhs(1,3)=clhs10*clhs80 + clhs44*clhs7 + clhs82;
lhs(1,4)=clhs10*clhs83 + clhs54 + clhs57*clhs7 + clhs85;
lhs(1,5)=DN(1,1)*clhs40 + clhs39*clhs87 - clhs86;
lhs(1,6)=clhs10*clhs88 + clhs63*clhs7 + clhs89;
lhs(1,7)=clhs10*clhs90 + clhs7*clhs73 + clhs70 + clhs92;
lhs(1,8)=DN(2,1)*clhs40 + clhs39*clhs94 - clhs93;
lhs(2,0)=clhs7*clhs96;
lhs(2,1)=clhs10*clhs96;
lhs(2,2)=clhs24*(clhs0 + clhs78);
lhs(2,3)=clhs6*(DN(0,0)*clhs97 + clhs60);
lhs(2,4)=clhs6*(DN(0,1)*clhs97 + clhs87);
lhs(2,5)=clhs98;
lhs(2,6)=clhs6*(DN(0,0)*clhs99 + clhs76);
lhs(2,7)=clhs6*(DN(0,1)*clhs99 + clhs94);
lhs(2,8)=clhs100;
lhs(3,0)=clhs101*clhs5 + clhs102*clhs9 + clhs105 + clhs49;
lhs(3,1)=clhs101*clhs33 + clhs102*clhs35 + clhs82;
lhs(3,2)=DN(0,0)*clhs106 + clhs39*clhs59 - clhs60;
lhs(3,3)=clhs101*clhs42 + clhs102*clhs44 + clhs107*clhs4 + clhs109;
lhs(3,4)=clhs101*clhs55 + clhs102*clhs57 + clhs111;
lhs(3,5)=DN(1,0)*clhs112;
lhs(3,6)=clhs101*clhs61 + clhs102*clhs63 + clhs115 + clhs116;
lhs(3,7)=clhs101*clhs71 + clhs102*clhs73 + clhs117;
lhs(3,8)=DN(2,0)*clhs106 - clhs118 + clhs119*clhs39;
lhs(4,0)=clhs101*clhs9 + clhs102*clhs77 + clhs58;
lhs(4,1)=clhs101*clhs35 + clhs102*clhs79 + clhs105 + clhs85;
lhs(4,2)=DN(0,1)*clhs106 + clhs39*clhs86 - clhs87;
lhs(4,3)=clhs101*clhs44 + clhs102*clhs80 + clhs111;
lhs(4,4)=clhs101*clhs57 + clhs102*clhs83 + clhs109 + clhs120*clhs4;
lhs(4,5)=DN(1,1)*clhs112;
lhs(4,6)=clhs101*clhs63 + clhs102*clhs88 + clhs122;
lhs(4,7)=clhs101*clhs73 + clhs102*clhs90 + clhs116 + clhs124;
lhs(4,8)=DN(2,1)*clhs106 - clhs125 + clhs126*clhs39;
lhs(5,0)=clhs6*(DN(1,0)*clhs95 + clhs59);
lhs(5,1)=clhs6*(DN(1,1)*clhs95 + clhs86);
lhs(5,2)=clhs98;
lhs(5,3)=clhs101*clhs127;
lhs(5,4)=clhs102*clhs127;
lhs(5,5)=clhs24*(clhs107 + clhs120);
lhs(5,6)=clhs6*(DN(1,0)*clhs99 + clhs119);
lhs(5,7)=clhs6*(DN(1,1)*clhs99 + clhs126);
lhs(5,8)=clhs128;
lhs(6,0)=clhs129*clhs5 + clhs130*clhs9 + clhs133 + clhs67;
lhs(6,1)=clhs129*clhs33 + clhs130*clhs35 + clhs89;
lhs(6,2)=DN(0,0)*clhs134 + clhs39*clhs75 - clhs76;
lhs(6,3)=clhs115 + clhs129*clhs42 + clhs130*clhs44 + clhs135;
lhs(6,4)=clhs122 + clhs129*clhs55 + clhs130*clhs57;
lhs(6,5)=DN(1,0)*clhs134 + clhs118*clhs39 - clhs119;
lhs(6,6)=clhs129*clhs61 + clhs130*clhs63 + clhs136*clhs4 + clhs138;
lhs(6,7)=clhs129*clhs71 + clhs130*clhs73 + clhs139;
lhs(6,8)=DN(2,0)*clhs140;
lhs(7,0)=clhs129*clhs9 + clhs130*clhs77 + clhs74;
lhs(7,1)=clhs129*clhs35 + clhs130*clhs79 + clhs133 + clhs92;
lhs(7,2)=DN(0,1)*clhs134 + clhs39*clhs93 - clhs94;
lhs(7,3)=clhs117 + clhs129*clhs44 + clhs130*clhs80;
lhs(7,4)=clhs124 + clhs129*clhs57 + clhs130*clhs83 + clhs135;
lhs(7,5)=DN(1,1)*clhs134 + clhs125*clhs39 - clhs126;
lhs(7,6)=clhs129*clhs63 + clhs130*clhs88 + clhs139;
lhs(7,7)=clhs129*clhs73 + clhs130*clhs90 + clhs138 + clhs141*clhs4;
lhs(7,8)=DN(2,1)*clhs140;
lhs(8,0)=clhs6*(DN(2,0)*clhs95 + clhs75);
lhs(8,1)=clhs6*(DN(2,1)*clhs95 + clhs93);
lhs(8,2)=clhs100;
lhs(8,3)=clhs6*(DN(2,0)*clhs97 + clhs118);
lhs(8,4)=clhs6*(DN(2,1)*clhs97 + clhs125);
lhs(8,5)=clhs128;
lhs(8,6)=clhs129*clhs142;
lhs(8,7)=clhs130*clhs142;
lhs(8,8)=clhs24*(clhs136 + clhs141);


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
    const double volume_error_ratio = rData.VolumeError;
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

    const double clhs0 = pow(DN(0,0), 2);
const double clhs1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs3 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs4 = rho*stab_c2*sqrt(pow(clhs1, 2) + pow(clhs2, 2) + pow(clhs3, 2));
const double clhs5 = clhs4*h/stab_c1 + mu;
const double clhs6 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs7 = 1.0/(max_spectral_radius + 1);
const double clhs8 = DN(0,0)*clhs7;
const double clhs9 = C(0,3)*DN(0,0);
const double clhs10 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs9;
const double clhs11 = DN(0,1)*clhs7;
const double clhs12 = C(0,5)*DN(0,0);
const double clhs13 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs12;
const double clhs14 = DN(0,2)*clhs7;
const double clhs15 = DN(0,0)*clhs1 + DN(0,1)*clhs2 + DN(0,2)*clhs3;
const double clhs16 = clhs7*rho;
const double clhs17 = N[0]*clhs16;
const double clhs18 = pow(N[0], 2);
const double clhs19 = clhs16*volume_error_ratio;
const double clhs20 = N[0]*volume_error_ratio;
const double clhs21 = 0.5*max_spectral_radius - 1.5;
const double clhs22 = N[0]*clhs21;
const double clhs23 = 1.0/dt;
const double clhs24 = clhs7*(0.5*max_spectral_radius - 1.5);
const double clhs25 = clhs23/(-clhs24 - clhs7 + 0.5);
const double clhs26 = clhs15 + clhs20 - clhs22*clhs25;
const double clhs27 = clhs23*rho;
const double clhs28 = 1.0/(clhs27*dyn_tau + clhs4/h + mu*stab_c1/pow(h, 2));
const double clhs29 = clhs28*clhs7*pow(rho, 2);
const double clhs30 = clhs15*clhs29;
const double clhs31 = clhs27*clhs7/(clhs24 + clhs7 - 0.5);
const double clhs32 = clhs21*clhs31;
const double clhs33 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs34 = clhs29*clhs33;
const double clhs35 = N[0]*clhs34;
const double clhs36 = clhs15*clhs17 + clhs18*clhs19 + clhs18*clhs32 + clhs26*clhs30 + clhs26*clhs35;
const double clhs37 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs9;
const double clhs38 = C(1,3)*DN(0,1);
const double clhs39 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs38;
const double clhs40 = C(3,5)*DN(0,0);
const double clhs41 = C(4,5)*DN(0,2);
const double clhs42 = C(1,5)*DN(0,1) + clhs40 + clhs41;
const double clhs43 = DN(0,0)*clhs5;
const double clhs44 = DN(0,1)*clhs43;
const double clhs45 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs12;
const double clhs46 = C(3,4)*DN(0,1);
const double clhs47 = C(2,3)*DN(0,2) + clhs40 + clhs46;
const double clhs48 = C(2,5)*DN(0,2);
const double clhs49 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs48;
const double clhs50 = DN(0,2)*clhs43;
const double clhs51 = clhs28*rho;
const double clhs52 = clhs33*clhs51;
const double clhs53 = clhs15*clhs51;
const double clhs54 = N[0]*clhs52 - N[0] + clhs53;
const double clhs55 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs56 = C(0,3)*DN(1,0);
const double clhs57 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs56;
const double clhs58 = C(0,5)*DN(1,0);
const double clhs59 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs58;
const double clhs60 = DN(0,0)*DN(1,0);
const double clhs61 = N[1]*clhs16;
const double clhs62 = clhs22*clhs31;
const double clhs63 = N[1]*clhs62 + clhs20*clhs61;
const double clhs64 = clhs5*clhs60 + clhs63;
const double clhs65 = DN(1,0)*clhs1 + DN(1,1)*clhs2 + DN(1,2)*clhs3;
const double clhs66 = N[1]*volume_error_ratio;
const double clhs67 = clhs21*clhs25;
const double clhs68 = -N[1]*clhs67 + clhs65 + clhs66;
const double clhs69 = clhs17*clhs65 + clhs30*clhs68 + clhs35*clhs68;
const double clhs70 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs56;
const double clhs71 = C(1,3)*DN(1,1);
const double clhs72 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs71;
const double clhs73 = C(3,5)*DN(1,0);
const double clhs74 = C(4,5)*DN(1,2);
const double clhs75 = C(1,5)*DN(1,1) + clhs73 + clhs74;
const double clhs76 = DN(1,1)*clhs43;
const double clhs77 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs58;
const double clhs78 = C(3,4)*DN(1,1);
const double clhs79 = C(2,3)*DN(1,2) + clhs73 + clhs78;
const double clhs80 = C(2,5)*DN(1,2);
const double clhs81 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs80;
const double clhs82 = DN(1,2)*clhs43;
const double clhs83 = DN(0,0)*N[1];
const double clhs84 = DN(1,0)*N[0];
const double clhs85 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs86 = C(0,3)*DN(2,0);
const double clhs87 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs86;
const double clhs88 = C(0,5)*DN(2,0);
const double clhs89 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs88;
const double clhs90 = DN(0,0)*DN(2,0);
const double clhs91 = N[2]*clhs16;
const double clhs92 = N[2]*clhs62 + clhs20*clhs91;
const double clhs93 = clhs5*clhs90 + clhs92;
const double clhs94 = DN(2,0)*clhs1 + DN(2,1)*clhs2 + DN(2,2)*clhs3;
const double clhs95 = N[2]*volume_error_ratio;
const double clhs96 = -N[2]*clhs67 + clhs94 + clhs95;
const double clhs97 = clhs17*clhs94 + clhs30*clhs96 + clhs35*clhs96;
const double clhs98 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs86;
const double clhs99 = C(1,3)*DN(2,1);
const double clhs100 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs99;
const double clhs101 = C(3,5)*DN(2,0);
const double clhs102 = C(4,5)*DN(2,2);
const double clhs103 = C(1,5)*DN(2,1) + clhs101 + clhs102;
const double clhs104 = DN(2,1)*clhs43;
const double clhs105 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs88;
const double clhs106 = C(3,4)*DN(2,1);
const double clhs107 = C(2,3)*DN(2,2) + clhs101 + clhs106;
const double clhs108 = C(2,5)*DN(2,2);
const double clhs109 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs108;
const double clhs110 = DN(2,2)*clhs43;
const double clhs111 = DN(0,0)*N[2];
const double clhs112 = DN(2,0)*N[0];
const double clhs113 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs114 = C(0,3)*DN(3,0);
const double clhs115 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs114;
const double clhs116 = C(0,5)*DN(3,0);
const double clhs117 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs116;
const double clhs118 = DN(0,0)*DN(3,0);
const double clhs119 = N[3]*clhs16;
const double clhs120 = N[3]*clhs62 + clhs119*clhs20;
const double clhs121 = clhs118*clhs5 + clhs120;
const double clhs122 = DN(3,0)*clhs1 + DN(3,1)*clhs2 + DN(3,2)*clhs3;
const double clhs123 = -N[3]*clhs67 + N[3]*volume_error_ratio + clhs122;
const double clhs124 = clhs122*clhs17 + clhs123*clhs30 + clhs123*clhs35;
const double clhs125 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs114;
const double clhs126 = C(1,3)*DN(3,1);
const double clhs127 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs126;
const double clhs128 = C(3,5)*DN(3,0);
const double clhs129 = C(4,5)*DN(3,2);
const double clhs130 = C(1,5)*DN(3,1) + clhs128 + clhs129;
const double clhs131 = DN(3,1)*clhs43;
const double clhs132 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs116;
const double clhs133 = C(3,4)*DN(3,1);
const double clhs134 = C(2,3)*DN(3,2) + clhs128 + clhs133;
const double clhs135 = C(2,5)*DN(3,2);
const double clhs136 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs135;
const double clhs137 = DN(3,2)*clhs43;
const double clhs138 = DN(0,0)*N[3];
const double clhs139 = DN(3,0)*N[0];
const double clhs140 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs38;
const double clhs141 = C(0,4)*DN(0,0) + clhs41 + clhs46;
const double clhs142 = pow(DN(0,1), 2);
const double clhs143 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs144 = C(1,4)*DN(0,1);
const double clhs145 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs144;
const double clhs146 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs144;
const double clhs147 = C(2,4)*DN(0,2);
const double clhs148 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs147;
const double clhs149 = DN(0,1)*clhs5;
const double clhs150 = DN(0,2)*clhs149;
const double clhs151 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs71;
const double clhs152 = C(0,4)*DN(1,0) + clhs74 + clhs78;
const double clhs153 = DN(1,0)*clhs149;
const double clhs154 = DN(0,1)*DN(1,1);
const double clhs155 = clhs154*clhs5;
const double clhs156 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs157 = C(1,4)*DN(1,1);
const double clhs158 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs157;
const double clhs159 = clhs63 + clhs69;
const double clhs160 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs157;
const double clhs161 = C(2,4)*DN(1,2);
const double clhs162 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs161;
const double clhs163 = DN(1,2)*clhs149;
const double clhs164 = DN(0,1)*N[1];
const double clhs165 = DN(1,1)*N[0];
const double clhs166 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs99;
const double clhs167 = C(0,4)*DN(2,0) + clhs102 + clhs106;
const double clhs168 = DN(2,0)*clhs149;
const double clhs169 = DN(0,1)*DN(2,1);
const double clhs170 = clhs169*clhs5;
const double clhs171 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs172 = C(1,4)*DN(2,1);
const double clhs173 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs172;
const double clhs174 = clhs92 + clhs97;
const double clhs175 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs172;
const double clhs176 = C(2,4)*DN(2,2);
const double clhs177 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs176;
const double clhs178 = DN(2,2)*clhs149;
const double clhs179 = DN(0,1)*N[2];
const double clhs180 = DN(2,1)*N[0];
const double clhs181 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs126;
const double clhs182 = C(0,4)*DN(3,0) + clhs129 + clhs133;
const double clhs183 = DN(3,0)*clhs149;
const double clhs184 = DN(0,1)*DN(3,1);
const double clhs185 = clhs184*clhs5;
const double clhs186 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs187 = C(1,4)*DN(3,1);
const double clhs188 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs187;
const double clhs189 = clhs120 + clhs124;
const double clhs190 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs187;
const double clhs191 = C(2,4)*DN(3,2);
const double clhs192 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs191;
const double clhs193 = DN(3,2)*clhs149;
const double clhs194 = DN(0,1)*N[3];
const double clhs195 = DN(3,1)*N[0];
const double clhs196 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs48;
const double clhs197 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs147;
const double clhs198 = pow(DN(0,2), 2);
const double clhs199 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs200 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs80;
const double clhs201 = DN(0,2)*clhs5;
const double clhs202 = DN(1,0)*clhs201;
const double clhs203 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs161;
const double clhs204 = DN(1,1)*clhs201;
const double clhs205 = DN(0,2)*DN(1,2);
const double clhs206 = clhs205*clhs5;
const double clhs207 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs208 = DN(0,2)*N[1];
const double clhs209 = DN(1,2)*N[0];
const double clhs210 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs108;
const double clhs211 = DN(2,0)*clhs201;
const double clhs212 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs176;
const double clhs213 = DN(2,1)*clhs201;
const double clhs214 = DN(0,2)*DN(2,2);
const double clhs215 = clhs214*clhs5;
const double clhs216 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs217 = DN(0,2)*N[2];
const double clhs218 = DN(2,2)*N[0];
const double clhs219 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs135;
const double clhs220 = DN(3,0)*clhs201;
const double clhs221 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs191;
const double clhs222 = DN(3,1)*clhs201;
const double clhs223 = DN(0,2)*DN(3,2);
const double clhs224 = clhs223*clhs5;
const double clhs225 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs226 = DN(0,2)*N[3];
const double clhs227 = DN(3,2)*N[0];
const double clhs228 = clhs26*clhs51;
const double clhs229 = N[0] + clhs228;
const double clhs230 = clhs51*clhs68;
const double clhs231 = clhs28*(clhs154 + clhs205 + clhs60);
const double clhs232 = clhs51*clhs96;
const double clhs233 = clhs28*(clhs169 + clhs214 + clhs90);
const double clhs234 = clhs123*clhs51;
const double clhs235 = clhs28*(clhs118 + clhs184 + clhs223);
const double clhs236 = DN(1,0)*clhs7;
const double clhs237 = DN(1,1)*clhs7;
const double clhs238 = DN(1,2)*clhs7;
const double clhs239 = clhs29*clhs65;
const double clhs240 = N[1]*clhs34;
const double clhs241 = clhs15*clhs61 + clhs239*clhs26 + clhs240*clhs26;
const double clhs242 = clhs51*clhs65;
const double clhs243 = pow(DN(1,0), 2);
const double clhs244 = pow(N[1], 2);
const double clhs245 = clhs19*clhs244 + clhs239*clhs68 + clhs240*clhs68 + clhs244*clhs32 + clhs61*clhs65;
const double clhs246 = DN(1,0)*clhs5;
const double clhs247 = DN(1,1)*clhs246;
const double clhs248 = DN(1,2)*clhs246;
const double clhs249 = N[1]*clhs52 - N[1] + clhs242;
const double clhs250 = DN(1,0)*DN(2,0);
const double clhs251 = N[1]*clhs32;
const double clhs252 = N[2]*clhs251 + clhs66*clhs91;
const double clhs253 = clhs250*clhs5 + clhs252;
const double clhs254 = clhs239*clhs96 + clhs240*clhs96 + clhs61*clhs94;
const double clhs255 = DN(2,1)*clhs246;
const double clhs256 = DN(2,2)*clhs246;
const double clhs257 = DN(1,0)*N[2];
const double clhs258 = DN(2,0)*N[1];
const double clhs259 = DN(1,0)*DN(3,0);
const double clhs260 = N[3]*clhs251 + clhs119*clhs66;
const double clhs261 = clhs259*clhs5 + clhs260;
const double clhs262 = clhs122*clhs61 + clhs123*clhs239 + clhs123*clhs240;
const double clhs263 = DN(3,1)*clhs246;
const double clhs264 = DN(3,2)*clhs246;
const double clhs265 = DN(1,0)*N[3];
const double clhs266 = DN(3,0)*N[1];
const double clhs267 = clhs241 + clhs63;
const double clhs268 = pow(DN(1,1), 2);
const double clhs269 = DN(1,1)*clhs5;
const double clhs270 = DN(1,2)*clhs269;
const double clhs271 = DN(2,0)*clhs269;
const double clhs272 = DN(1,1)*DN(2,1);
const double clhs273 = clhs272*clhs5;
const double clhs274 = clhs252 + clhs254;
const double clhs275 = DN(2,2)*clhs269;
const double clhs276 = DN(1,1)*N[2];
const double clhs277 = DN(2,1)*N[1];
const double clhs278 = DN(3,0)*clhs269;
const double clhs279 = DN(1,1)*DN(3,1);
const double clhs280 = clhs279*clhs5;
const double clhs281 = clhs260 + clhs262;
const double clhs282 = DN(3,2)*clhs269;
const double clhs283 = DN(1,1)*N[3];
const double clhs284 = DN(3,1)*N[1];
const double clhs285 = pow(DN(1,2), 2);
const double clhs286 = DN(1,2)*clhs5;
const double clhs287 = DN(2,0)*clhs286;
const double clhs288 = DN(2,1)*clhs286;
const double clhs289 = DN(1,2)*DN(2,2);
const double clhs290 = clhs289*clhs5;
const double clhs291 = DN(1,2)*N[2];
const double clhs292 = DN(2,2)*N[1];
const double clhs293 = DN(3,0)*clhs286;
const double clhs294 = DN(3,1)*clhs286;
const double clhs295 = DN(1,2)*DN(3,2);
const double clhs296 = clhs295*clhs5;
const double clhs297 = DN(1,2)*N[3];
const double clhs298 = DN(3,2)*N[1];
const double clhs299 = N[1] + clhs230;
const double clhs300 = clhs28*(clhs250 + clhs272 + clhs289);
const double clhs301 = clhs28*(clhs259 + clhs279 + clhs295);
const double clhs302 = DN(2,0)*clhs7;
const double clhs303 = DN(2,1)*clhs7;
const double clhs304 = DN(2,2)*clhs7;
const double clhs305 = clhs29*clhs94;
const double clhs306 = N[2]*clhs34;
const double clhs307 = clhs15*clhs91 + clhs26*clhs305 + clhs26*clhs306;
const double clhs308 = clhs51*clhs94;
const double clhs309 = clhs305*clhs68 + clhs306*clhs68 + clhs65*clhs91;
const double clhs310 = pow(DN(2,0), 2);
const double clhs311 = pow(N[2], 2);
const double clhs312 = clhs19*clhs311 + clhs305*clhs96 + clhs306*clhs96 + clhs311*clhs32 + clhs91*clhs94;
const double clhs313 = DN(2,0)*clhs5;
const double clhs314 = DN(2,1)*clhs313;
const double clhs315 = DN(2,2)*clhs313;
const double clhs316 = N[2]*clhs52 - N[2] + clhs308;
const double clhs317 = DN(2,0)*DN(3,0);
const double clhs318 = N[2]*N[3]*clhs32 + clhs119*clhs95;
const double clhs319 = clhs317*clhs5 + clhs318;
const double clhs320 = clhs122*clhs91 + clhs123*clhs305 + clhs123*clhs306;
const double clhs321 = DN(3,1)*clhs313;
const double clhs322 = DN(3,2)*clhs313;
const double clhs323 = DN(2,0)*N[3];
const double clhs324 = DN(3,0)*N[2];
const double clhs325 = clhs307 + clhs92;
const double clhs326 = clhs252 + clhs309;
const double clhs327 = pow(DN(2,1), 2);
const double clhs328 = DN(2,1)*clhs5;
const double clhs329 = DN(2,2)*clhs328;
const double clhs330 = DN(3,0)*clhs328;
const double clhs331 = DN(2,1)*DN(3,1);
const double clhs332 = clhs331*clhs5;
const double clhs333 = clhs318 + clhs320;
const double clhs334 = DN(3,2)*clhs328;
const double clhs335 = DN(2,1)*N[3];
const double clhs336 = DN(3,1)*N[2];
const double clhs337 = pow(DN(2,2), 2);
const double clhs338 = DN(2,2)*clhs5;
const double clhs339 = DN(3,0)*clhs338;
const double clhs340 = DN(3,1)*clhs338;
const double clhs341 = DN(2,2)*DN(3,2);
const double clhs342 = clhs341*clhs5;
const double clhs343 = DN(2,2)*N[3];
const double clhs344 = DN(3,2)*N[2];
const double clhs345 = N[2] + clhs232;
const double clhs346 = clhs28*(clhs317 + clhs331 + clhs341);
const double clhs347 = DN(3,0)*clhs7;
const double clhs348 = DN(3,1)*clhs7;
const double clhs349 = DN(3,2)*clhs7;
const double clhs350 = clhs122*clhs29;
const double clhs351 = N[3]*clhs34;
const double clhs352 = clhs119*clhs15 + clhs26*clhs350 + clhs26*clhs351;
const double clhs353 = clhs122*clhs51;
const double clhs354 = clhs119*clhs65 + clhs350*clhs68 + clhs351*clhs68;
const double clhs355 = clhs119*clhs94 + clhs350*clhs96 + clhs351*clhs96;
const double clhs356 = pow(DN(3,0), 2);
const double clhs357 = pow(N[3], 2);
const double clhs358 = clhs119*clhs122 + clhs123*clhs350 + clhs123*clhs351 + clhs19*clhs357 + clhs32*clhs357;
const double clhs359 = DN(3,0)*clhs5;
const double clhs360 = DN(3,1)*clhs359;
const double clhs361 = DN(3,2)*clhs359;
const double clhs362 = N[3]*clhs52 - N[3] + clhs353;
const double clhs363 = clhs120 + clhs352;
const double clhs364 = clhs260 + clhs354;
const double clhs365 = clhs318 + clhs355;
const double clhs366 = pow(DN(3,1), 2);
const double clhs367 = DN(3,1)*DN(3,2)*clhs5;
const double clhs368 = pow(DN(3,2), 2);
const double clhs369 = N[3] + clhs234;
lhs(0,0)=clhs0*clhs5 + clhs10*clhs11 + clhs13*clhs14 + clhs36 + clhs6*clhs8;
lhs(0,1)=clhs11*clhs39 + clhs14*clhs42 + clhs37*clhs8 + clhs44;
lhs(0,2)=clhs11*clhs47 + clhs14*clhs49 + clhs45*clhs8 + clhs50;
lhs(0,3)=DN(0,0)*clhs54;
lhs(0,4)=clhs11*clhs57 + clhs14*clhs59 + clhs55*clhs8 + clhs64 + clhs69;
lhs(0,5)=clhs11*clhs72 + clhs14*clhs75 + clhs70*clhs8 + clhs76;
lhs(0,6)=clhs11*clhs79 + clhs14*clhs81 + clhs77*clhs8 + clhs82;
lhs(0,7)=DN(1,0)*clhs53 + clhs52*clhs84 - clhs83;
lhs(0,8)=clhs11*clhs87 + clhs14*clhs89 + clhs8*clhs85 + clhs93 + clhs97;
lhs(0,9)=clhs100*clhs11 + clhs103*clhs14 + clhs104 + clhs8*clhs98;
lhs(0,10)=clhs105*clhs8 + clhs107*clhs11 + clhs109*clhs14 + clhs110;
lhs(0,11)=DN(2,0)*clhs53 - clhs111 + clhs112*clhs52;
lhs(0,12)=clhs11*clhs115 + clhs113*clhs8 + clhs117*clhs14 + clhs121 + clhs124;
lhs(0,13)=clhs11*clhs127 + clhs125*clhs8 + clhs130*clhs14 + clhs131;
lhs(0,14)=clhs11*clhs134 + clhs132*clhs8 + clhs136*clhs14 + clhs137;
lhs(0,15)=DN(3,0)*clhs53 - clhs138 + clhs139*clhs52;
lhs(1,0)=clhs10*clhs8 + clhs11*clhs140 + clhs14*clhs141 + clhs44;
lhs(1,1)=clhs11*clhs143 + clhs14*clhs145 + clhs142*clhs5 + clhs36 + clhs39*clhs8;
lhs(1,2)=clhs11*clhs146 + clhs14*clhs148 + clhs150 + clhs47*clhs8;
lhs(1,3)=DN(0,1)*clhs54;
lhs(1,4)=clhs11*clhs151 + clhs14*clhs152 + clhs153 + clhs57*clhs8;
lhs(1,5)=clhs11*clhs156 + clhs14*clhs158 + clhs155 + clhs159 + clhs72*clhs8;
lhs(1,6)=clhs11*clhs160 + clhs14*clhs162 + clhs163 + clhs79*clhs8;
lhs(1,7)=DN(1,1)*clhs53 - clhs164 + clhs165*clhs52;
lhs(1,8)=clhs11*clhs166 + clhs14*clhs167 + clhs168 + clhs8*clhs87;
lhs(1,9)=clhs100*clhs8 + clhs11*clhs171 + clhs14*clhs173 + clhs170 + clhs174;
lhs(1,10)=clhs107*clhs8 + clhs11*clhs175 + clhs14*clhs177 + clhs178;
lhs(1,11)=DN(2,1)*clhs53 - clhs179 + clhs180*clhs52;
lhs(1,12)=clhs11*clhs181 + clhs115*clhs8 + clhs14*clhs182 + clhs183;
lhs(1,13)=clhs11*clhs186 + clhs127*clhs8 + clhs14*clhs188 + clhs185 + clhs189;
lhs(1,14)=clhs11*clhs190 + clhs134*clhs8 + clhs14*clhs192 + clhs193;
lhs(1,15)=DN(3,1)*clhs53 - clhs194 + clhs195*clhs52;
lhs(2,0)=clhs11*clhs141 + clhs13*clhs8 + clhs14*clhs196 + clhs50;
lhs(2,1)=clhs11*clhs145 + clhs14*clhs197 + clhs150 + clhs42*clhs8;
lhs(2,2)=clhs11*clhs148 + clhs14*clhs199 + clhs198*clhs5 + clhs36 + clhs49*clhs8;
lhs(2,3)=DN(0,2)*clhs54;
lhs(2,4)=clhs11*clhs152 + clhs14*clhs200 + clhs202 + clhs59*clhs8;
lhs(2,5)=clhs11*clhs158 + clhs14*clhs203 + clhs204 + clhs75*clhs8;
lhs(2,6)=clhs11*clhs162 + clhs14*clhs207 + clhs159 + clhs206 + clhs8*clhs81;
lhs(2,7)=DN(1,2)*clhs53 - clhs208 + clhs209*clhs52;
lhs(2,8)=clhs11*clhs167 + clhs14*clhs210 + clhs211 + clhs8*clhs89;
lhs(2,9)=clhs103*clhs8 + clhs11*clhs173 + clhs14*clhs212 + clhs213;
lhs(2,10)=clhs109*clhs8 + clhs11*clhs177 + clhs14*clhs216 + clhs174 + clhs215;
lhs(2,11)=DN(2,2)*clhs53 - clhs217 + clhs218*clhs52;
lhs(2,12)=clhs11*clhs182 + clhs117*clhs8 + clhs14*clhs219 + clhs220;
lhs(2,13)=clhs11*clhs188 + clhs130*clhs8 + clhs14*clhs221 + clhs222;
lhs(2,14)=clhs11*clhs192 + clhs136*clhs8 + clhs14*clhs225 + clhs189 + clhs224;
lhs(2,15)=DN(3,2)*clhs53 - clhs226 + clhs227*clhs52;
lhs(3,0)=clhs229*clhs8;
lhs(3,1)=clhs11*clhs229;
lhs(3,2)=clhs14*clhs229;
lhs(3,3)=clhs28*(clhs0 + clhs142 + clhs198);
lhs(3,4)=clhs7*(DN(0,0)*clhs230 + clhs84);
lhs(3,5)=clhs7*(DN(0,1)*clhs230 + clhs165);
lhs(3,6)=clhs7*(DN(0,2)*clhs230 + clhs209);
lhs(3,7)=clhs231;
lhs(3,8)=clhs7*(DN(0,0)*clhs232 + clhs112);
lhs(3,9)=clhs7*(DN(0,1)*clhs232 + clhs180);
lhs(3,10)=clhs7*(DN(0,2)*clhs232 + clhs218);
lhs(3,11)=clhs233;
lhs(3,12)=clhs7*(DN(0,0)*clhs234 + clhs139);
lhs(3,13)=clhs7*(DN(0,1)*clhs234 + clhs195);
lhs(3,14)=clhs7*(DN(0,2)*clhs234 + clhs227);
lhs(3,15)=clhs235;
lhs(4,0)=clhs10*clhs237 + clhs13*clhs238 + clhs236*clhs6 + clhs241 + clhs64;
lhs(4,1)=clhs153 + clhs236*clhs37 + clhs237*clhs39 + clhs238*clhs42;
lhs(4,2)=clhs202 + clhs236*clhs45 + clhs237*clhs47 + clhs238*clhs49;
lhs(4,3)=DN(0,0)*clhs242 + clhs52*clhs83 - clhs84;
lhs(4,4)=clhs236*clhs55 + clhs237*clhs57 + clhs238*clhs59 + clhs243*clhs5 + clhs245;
lhs(4,5)=clhs236*clhs70 + clhs237*clhs72 + clhs238*clhs75 + clhs247;
lhs(4,6)=clhs236*clhs77 + clhs237*clhs79 + clhs238*clhs81 + clhs248;
lhs(4,7)=DN(1,0)*clhs249;
lhs(4,8)=clhs236*clhs85 + clhs237*clhs87 + clhs238*clhs89 + clhs253 + clhs254;
lhs(4,9)=clhs100*clhs237 + clhs103*clhs238 + clhs236*clhs98 + clhs255;
lhs(4,10)=clhs105*clhs236 + clhs107*clhs237 + clhs109*clhs238 + clhs256;
lhs(4,11)=DN(2,0)*clhs242 - clhs257 + clhs258*clhs52;
lhs(4,12)=clhs113*clhs236 + clhs115*clhs237 + clhs117*clhs238 + clhs261 + clhs262;
lhs(4,13)=clhs125*clhs236 + clhs127*clhs237 + clhs130*clhs238 + clhs263;
lhs(4,14)=clhs132*clhs236 + clhs134*clhs237 + clhs136*clhs238 + clhs264;
lhs(4,15)=DN(3,0)*clhs242 - clhs265 + clhs266*clhs52;
lhs(5,0)=clhs10*clhs236 + clhs140*clhs237 + clhs141*clhs238 + clhs76;
lhs(5,1)=clhs143*clhs237 + clhs145*clhs238 + clhs155 + clhs236*clhs39 + clhs267;
lhs(5,2)=clhs146*clhs237 + clhs148*clhs238 + clhs204 + clhs236*clhs47;
lhs(5,3)=DN(0,1)*clhs242 + clhs164*clhs52 - clhs165;
lhs(5,4)=clhs151*clhs237 + clhs152*clhs238 + clhs236*clhs57 + clhs247;
lhs(5,5)=clhs156*clhs237 + clhs158*clhs238 + clhs236*clhs72 + clhs245 + clhs268*clhs5;
lhs(5,6)=clhs160*clhs237 + clhs162*clhs238 + clhs236*clhs79 + clhs270;
lhs(5,7)=DN(1,1)*clhs249;
lhs(5,8)=clhs166*clhs237 + clhs167*clhs238 + clhs236*clhs87 + clhs271;
lhs(5,9)=clhs100*clhs236 + clhs171*clhs237 + clhs173*clhs238 + clhs273 + clhs274;
lhs(5,10)=clhs107*clhs236 + clhs175*clhs237 + clhs177*clhs238 + clhs275;
lhs(5,11)=DN(2,1)*clhs242 - clhs276 + clhs277*clhs52;
lhs(5,12)=clhs115*clhs236 + clhs181*clhs237 + clhs182*clhs238 + clhs278;
lhs(5,13)=clhs127*clhs236 + clhs186*clhs237 + clhs188*clhs238 + clhs280 + clhs281;
lhs(5,14)=clhs134*clhs236 + clhs190*clhs237 + clhs192*clhs238 + clhs282;
lhs(5,15)=DN(3,1)*clhs242 - clhs283 + clhs284*clhs52;
lhs(6,0)=clhs13*clhs236 + clhs141*clhs237 + clhs196*clhs238 + clhs82;
lhs(6,1)=clhs145*clhs237 + clhs163 + clhs197*clhs238 + clhs236*clhs42;
lhs(6,2)=clhs148*clhs237 + clhs199*clhs238 + clhs206 + clhs236*clhs49 + clhs267;
lhs(6,3)=DN(0,2)*clhs242 + clhs208*clhs52 - clhs209;
lhs(6,4)=clhs152*clhs237 + clhs200*clhs238 + clhs236*clhs59 + clhs248;
lhs(6,5)=clhs158*clhs237 + clhs203*clhs238 + clhs236*clhs75 + clhs270;
lhs(6,6)=clhs162*clhs237 + clhs207*clhs238 + clhs236*clhs81 + clhs245 + clhs285*clhs5;
lhs(6,7)=DN(1,2)*clhs249;
lhs(6,8)=clhs167*clhs237 + clhs210*clhs238 + clhs236*clhs89 + clhs287;
lhs(6,9)=clhs103*clhs236 + clhs173*clhs237 + clhs212*clhs238 + clhs288;
lhs(6,10)=clhs109*clhs236 + clhs177*clhs237 + clhs216*clhs238 + clhs274 + clhs290;
lhs(6,11)=DN(2,2)*clhs242 - clhs291 + clhs292*clhs52;
lhs(6,12)=clhs117*clhs236 + clhs182*clhs237 + clhs219*clhs238 + clhs293;
lhs(6,13)=clhs130*clhs236 + clhs188*clhs237 + clhs221*clhs238 + clhs294;
lhs(6,14)=clhs136*clhs236 + clhs192*clhs237 + clhs225*clhs238 + clhs281 + clhs296;
lhs(6,15)=DN(3,2)*clhs242 - clhs297 + clhs298*clhs52;
lhs(7,0)=clhs7*(DN(1,0)*clhs228 + clhs83);
lhs(7,1)=clhs7*(DN(1,1)*clhs228 + clhs164);
lhs(7,2)=clhs7*(DN(1,2)*clhs228 + clhs208);
lhs(7,3)=clhs231;
lhs(7,4)=clhs236*clhs299;
lhs(7,5)=clhs237*clhs299;
lhs(7,6)=clhs238*clhs299;
lhs(7,7)=clhs28*(clhs243 + clhs268 + clhs285);
lhs(7,8)=clhs7*(DN(1,0)*clhs232 + clhs258);
lhs(7,9)=clhs7*(DN(1,1)*clhs232 + clhs277);
lhs(7,10)=clhs7*(DN(1,2)*clhs232 + clhs292);
lhs(7,11)=clhs300;
lhs(7,12)=clhs7*(DN(1,0)*clhs234 + clhs266);
lhs(7,13)=clhs7*(DN(1,1)*clhs234 + clhs284);
lhs(7,14)=clhs7*(DN(1,2)*clhs234 + clhs298);
lhs(7,15)=clhs301;
lhs(8,0)=clhs10*clhs303 + clhs13*clhs304 + clhs302*clhs6 + clhs307 + clhs93;
lhs(8,1)=clhs168 + clhs302*clhs37 + clhs303*clhs39 + clhs304*clhs42;
lhs(8,2)=clhs211 + clhs302*clhs45 + clhs303*clhs47 + clhs304*clhs49;
lhs(8,3)=DN(0,0)*clhs308 + clhs111*clhs52 - clhs112;
lhs(8,4)=clhs253 + clhs302*clhs55 + clhs303*clhs57 + clhs304*clhs59 + clhs309;
lhs(8,5)=clhs271 + clhs302*clhs70 + clhs303*clhs72 + clhs304*clhs75;
lhs(8,6)=clhs287 + clhs302*clhs77 + clhs303*clhs79 + clhs304*clhs81;
lhs(8,7)=DN(1,0)*clhs308 + clhs257*clhs52 - clhs258;
lhs(8,8)=clhs302*clhs85 + clhs303*clhs87 + clhs304*clhs89 + clhs310*clhs5 + clhs312;
lhs(8,9)=clhs100*clhs303 + clhs103*clhs304 + clhs302*clhs98 + clhs314;
lhs(8,10)=clhs105*clhs302 + clhs107*clhs303 + clhs109*clhs304 + clhs315;
lhs(8,11)=DN(2,0)*clhs316;
lhs(8,12)=clhs113*clhs302 + clhs115*clhs303 + clhs117*clhs304 + clhs319 + clhs320;
lhs(8,13)=clhs125*clhs302 + clhs127*clhs303 + clhs130*clhs304 + clhs321;
lhs(8,14)=clhs132*clhs302 + clhs134*clhs303 + clhs136*clhs304 + clhs322;
lhs(8,15)=DN(3,0)*clhs308 - clhs323 + clhs324*clhs52;
lhs(9,0)=clhs10*clhs302 + clhs104 + clhs140*clhs303 + clhs141*clhs304;
lhs(9,1)=clhs143*clhs303 + clhs145*clhs304 + clhs170 + clhs302*clhs39 + clhs325;
lhs(9,2)=clhs146*clhs303 + clhs148*clhs304 + clhs213 + clhs302*clhs47;
lhs(9,3)=DN(0,1)*clhs308 + clhs179*clhs52 - clhs180;
lhs(9,4)=clhs151*clhs303 + clhs152*clhs304 + clhs255 + clhs302*clhs57;
lhs(9,5)=clhs156*clhs303 + clhs158*clhs304 + clhs273 + clhs302*clhs72 + clhs326;
lhs(9,6)=clhs160*clhs303 + clhs162*clhs304 + clhs288 + clhs302*clhs79;
lhs(9,7)=DN(1,1)*clhs308 + clhs276*clhs52 - clhs277;
lhs(9,8)=clhs166*clhs303 + clhs167*clhs304 + clhs302*clhs87 + clhs314;
lhs(9,9)=clhs100*clhs302 + clhs171*clhs303 + clhs173*clhs304 + clhs312 + clhs327*clhs5;
lhs(9,10)=clhs107*clhs302 + clhs175*clhs303 + clhs177*clhs304 + clhs329;
lhs(9,11)=DN(2,1)*clhs316;
lhs(9,12)=clhs115*clhs302 + clhs181*clhs303 + clhs182*clhs304 + clhs330;
lhs(9,13)=clhs127*clhs302 + clhs186*clhs303 + clhs188*clhs304 + clhs332 + clhs333;
lhs(9,14)=clhs134*clhs302 + clhs190*clhs303 + clhs192*clhs304 + clhs334;
lhs(9,15)=DN(3,1)*clhs308 - clhs335 + clhs336*clhs52;
lhs(10,0)=clhs110 + clhs13*clhs302 + clhs141*clhs303 + clhs196*clhs304;
lhs(10,1)=clhs145*clhs303 + clhs178 + clhs197*clhs304 + clhs302*clhs42;
lhs(10,2)=clhs148*clhs303 + clhs199*clhs304 + clhs215 + clhs302*clhs49 + clhs325;
lhs(10,3)=DN(0,2)*clhs308 + clhs217*clhs52 - clhs218;
lhs(10,4)=clhs152*clhs303 + clhs200*clhs304 + clhs256 + clhs302*clhs59;
lhs(10,5)=clhs158*clhs303 + clhs203*clhs304 + clhs275 + clhs302*clhs75;
lhs(10,6)=clhs162*clhs303 + clhs207*clhs304 + clhs290 + clhs302*clhs81 + clhs326;
lhs(10,7)=DN(1,2)*clhs308 + clhs291*clhs52 - clhs292;
lhs(10,8)=clhs167*clhs303 + clhs210*clhs304 + clhs302*clhs89 + clhs315;
lhs(10,9)=clhs103*clhs302 + clhs173*clhs303 + clhs212*clhs304 + clhs329;
lhs(10,10)=clhs109*clhs302 + clhs177*clhs303 + clhs216*clhs304 + clhs312 + clhs337*clhs5;
lhs(10,11)=DN(2,2)*clhs316;
lhs(10,12)=clhs117*clhs302 + clhs182*clhs303 + clhs219*clhs304 + clhs339;
lhs(10,13)=clhs130*clhs302 + clhs188*clhs303 + clhs221*clhs304 + clhs340;
lhs(10,14)=clhs136*clhs302 + clhs192*clhs303 + clhs225*clhs304 + clhs333 + clhs342;
lhs(10,15)=DN(3,2)*clhs308 - clhs343 + clhs344*clhs52;
lhs(11,0)=clhs7*(DN(2,0)*clhs228 + clhs111);
lhs(11,1)=clhs7*(DN(2,1)*clhs228 + clhs179);
lhs(11,2)=clhs7*(DN(2,2)*clhs228 + clhs217);
lhs(11,3)=clhs233;
lhs(11,4)=clhs7*(DN(2,0)*clhs230 + clhs257);
lhs(11,5)=clhs7*(DN(2,1)*clhs230 + clhs276);
lhs(11,6)=clhs7*(DN(2,2)*clhs230 + clhs291);
lhs(11,7)=clhs300;
lhs(11,8)=clhs302*clhs345;
lhs(11,9)=clhs303*clhs345;
lhs(11,10)=clhs304*clhs345;
lhs(11,11)=clhs28*(clhs310 + clhs327 + clhs337);
lhs(11,12)=clhs7*(DN(2,0)*clhs234 + clhs324);
lhs(11,13)=clhs7*(DN(2,1)*clhs234 + clhs336);
lhs(11,14)=clhs7*(DN(2,2)*clhs234 + clhs344);
lhs(11,15)=clhs346;
lhs(12,0)=clhs10*clhs348 + clhs121 + clhs13*clhs349 + clhs347*clhs6 + clhs352;
lhs(12,1)=clhs183 + clhs347*clhs37 + clhs348*clhs39 + clhs349*clhs42;
lhs(12,2)=clhs220 + clhs347*clhs45 + clhs348*clhs47 + clhs349*clhs49;
lhs(12,3)=DN(0,0)*clhs353 + clhs138*clhs52 - clhs139;
lhs(12,4)=clhs261 + clhs347*clhs55 + clhs348*clhs57 + clhs349*clhs59 + clhs354;
lhs(12,5)=clhs278 + clhs347*clhs70 + clhs348*clhs72 + clhs349*clhs75;
lhs(12,6)=clhs293 + clhs347*clhs77 + clhs348*clhs79 + clhs349*clhs81;
lhs(12,7)=DN(1,0)*clhs353 + clhs265*clhs52 - clhs266;
lhs(12,8)=clhs319 + clhs347*clhs85 + clhs348*clhs87 + clhs349*clhs89 + clhs355;
lhs(12,9)=clhs100*clhs348 + clhs103*clhs349 + clhs330 + clhs347*clhs98;
lhs(12,10)=clhs105*clhs347 + clhs107*clhs348 + clhs109*clhs349 + clhs339;
lhs(12,11)=DN(2,0)*clhs353 + clhs323*clhs52 - clhs324;
lhs(12,12)=clhs113*clhs347 + clhs115*clhs348 + clhs117*clhs349 + clhs356*clhs5 + clhs358;
lhs(12,13)=clhs125*clhs347 + clhs127*clhs348 + clhs130*clhs349 + clhs360;
lhs(12,14)=clhs132*clhs347 + clhs134*clhs348 + clhs136*clhs349 + clhs361;
lhs(12,15)=DN(3,0)*clhs362;
lhs(13,0)=clhs10*clhs347 + clhs131 + clhs140*clhs348 + clhs141*clhs349;
lhs(13,1)=clhs143*clhs348 + clhs145*clhs349 + clhs185 + clhs347*clhs39 + clhs363;
lhs(13,2)=clhs146*clhs348 + clhs148*clhs349 + clhs222 + clhs347*clhs47;
lhs(13,3)=DN(0,1)*clhs353 + clhs194*clhs52 - clhs195;
lhs(13,4)=clhs151*clhs348 + clhs152*clhs349 + clhs263 + clhs347*clhs57;
lhs(13,5)=clhs156*clhs348 + clhs158*clhs349 + clhs280 + clhs347*clhs72 + clhs364;
lhs(13,6)=clhs160*clhs348 + clhs162*clhs349 + clhs294 + clhs347*clhs79;
lhs(13,7)=DN(1,1)*clhs353 + clhs283*clhs52 - clhs284;
lhs(13,8)=clhs166*clhs348 + clhs167*clhs349 + clhs321 + clhs347*clhs87;
lhs(13,9)=clhs100*clhs347 + clhs171*clhs348 + clhs173*clhs349 + clhs332 + clhs365;
lhs(13,10)=clhs107*clhs347 + clhs175*clhs348 + clhs177*clhs349 + clhs340;
lhs(13,11)=DN(2,1)*clhs353 + clhs335*clhs52 - clhs336;
lhs(13,12)=clhs115*clhs347 + clhs181*clhs348 + clhs182*clhs349 + clhs360;
lhs(13,13)=clhs127*clhs347 + clhs186*clhs348 + clhs188*clhs349 + clhs358 + clhs366*clhs5;
lhs(13,14)=clhs134*clhs347 + clhs190*clhs348 + clhs192*clhs349 + clhs367;
lhs(13,15)=DN(3,1)*clhs362;
lhs(14,0)=clhs13*clhs347 + clhs137 + clhs141*clhs348 + clhs196*clhs349;
lhs(14,1)=clhs145*clhs348 + clhs193 + clhs197*clhs349 + clhs347*clhs42;
lhs(14,2)=clhs148*clhs348 + clhs199*clhs349 + clhs224 + clhs347*clhs49 + clhs363;
lhs(14,3)=DN(0,2)*clhs353 + clhs226*clhs52 - clhs227;
lhs(14,4)=clhs152*clhs348 + clhs200*clhs349 + clhs264 + clhs347*clhs59;
lhs(14,5)=clhs158*clhs348 + clhs203*clhs349 + clhs282 + clhs347*clhs75;
lhs(14,6)=clhs162*clhs348 + clhs207*clhs349 + clhs296 + clhs347*clhs81 + clhs364;
lhs(14,7)=DN(1,2)*clhs353 + clhs297*clhs52 - clhs298;
lhs(14,8)=clhs167*clhs348 + clhs210*clhs349 + clhs322 + clhs347*clhs89;
lhs(14,9)=clhs103*clhs347 + clhs173*clhs348 + clhs212*clhs349 + clhs334;
lhs(14,10)=clhs109*clhs347 + clhs177*clhs348 + clhs216*clhs349 + clhs342 + clhs365;
lhs(14,11)=DN(2,2)*clhs353 + clhs343*clhs52 - clhs344;
lhs(14,12)=clhs117*clhs347 + clhs182*clhs348 + clhs219*clhs349 + clhs361;
lhs(14,13)=clhs130*clhs347 + clhs188*clhs348 + clhs221*clhs349 + clhs367;
lhs(14,14)=clhs136*clhs347 + clhs192*clhs348 + clhs225*clhs349 + clhs358 + clhs368*clhs5;
lhs(14,15)=DN(3,2)*clhs362;
lhs(15,0)=clhs7*(DN(3,0)*clhs228 + clhs138);
lhs(15,1)=clhs7*(DN(3,1)*clhs228 + clhs194);
lhs(15,2)=clhs7*(DN(3,2)*clhs228 + clhs226);
lhs(15,3)=clhs235;
lhs(15,4)=clhs7*(DN(3,0)*clhs230 + clhs265);
lhs(15,5)=clhs7*(DN(3,1)*clhs230 + clhs283);
lhs(15,6)=clhs7*(DN(3,2)*clhs230 + clhs297);
lhs(15,7)=clhs301;
lhs(15,8)=clhs7*(DN(3,0)*clhs232 + clhs323);
lhs(15,9)=clhs7*(DN(3,1)*clhs232 + clhs335);
lhs(15,10)=clhs7*(DN(3,2)*clhs232 + clhs343);
lhs(15,11)=clhs346;
lhs(15,12)=clhs347*clhs369;
lhs(15,13)=clhs348*clhs369;
lhs(15,14)=clhs349*clhs369;
lhs(15,15)=clhs28*(clhs356 + clhs366 + clhs368);


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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = 1.0/(max_spectral_radius + 1);
const double crhs2 = rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)));
const double crhs3 = v(0,0) - vn(0,0);
const double crhs4 = crhs1*crhs3 + vn(0,0);
const double crhs5 = v(1,0) - vn(1,0);
const double crhs6 = crhs1*crhs5 + vn(1,0);
const double crhs7 = v(2,0) - vn(2,0);
const double crhs8 = crhs1*crhs7 + vn(2,0);
const double crhs9 = N[0]*crhs4 + N[1]*crhs6 + N[2]*crhs8;
const double crhs10 = N[0]*rho;
const double crhs11 = crhs10*volume_error_ratio;
const double crhs12 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs14 = rho*stab_c2*sqrt(pow(crhs12, 2) + pow(crhs13, 2));
const double crhs15 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs16 = (crhs15 - volume_error_ratio)*(crhs14*h/stab_c1 + mu);
const double crhs17 = rho*(crhs12*(DN(0,0)*crhs4 + DN(1,0)*crhs6 + DN(2,0)*crhs8) + crhs13*(DN(0,1)*crhs4 + DN(1,1)*crhs6 + DN(2,1)*crhs8));
const double crhs18 = -acceleration_alpha_method(0,0);
const double crhs19 = 1.0/dt;
const double crhs20 = 0.5*max_spectral_radius;
const double crhs21 = -crhs1;
const double crhs22 = crhs21 + 0.5;
const double crhs23 = 1.0/(-crhs1*(crhs20 - 1.5) + crhs22);
const double crhs24 = crhs19*crhs23;
const double crhs25 = crhs1*(1.5 - crhs20);
const double crhs26 = crhs23*(crhs1 - crhs25 + 0.5);
const double crhs27 = 0.5*crhs1;
const double crhs28 = crhs27*(max_spectral_radius - 3);
const double crhs29 = -acceleration_alpha_method(1,0);
const double crhs30 = -acceleration_alpha_method(2,0);
const double crhs31 = N[0]*(acceleration_alpha_method(0,0) - crhs28*(-acceleration_alpha_method(0,0)*crhs26 + crhs18 + crhs24*crhs3)) + N[1]*(acceleration_alpha_method(1,0) - crhs28*(-acceleration_alpha_method(1,0)*crhs26 + crhs24*crhs5 + crhs29)) + N[2]*(acceleration_alpha_method(2,0) - crhs28*(-acceleration_alpha_method(2,0)*crhs26 + crhs24*crhs7 + crhs30));
const double crhs32 = rho*volume_error_ratio;
const double crhs33 = crhs32*crhs9;
const double crhs34 = 1.0/(crhs22 + crhs25);
const double crhs35 = crhs19*crhs34;
const double crhs36 = crhs34*(crhs21 + crhs25 - 0.5);
const double crhs37 = crhs27*(3 - max_spectral_radius);
const double crhs38 = 1.0/(crhs14/h + crhs19*dyn_tau*rho + mu*stab_c1/pow(h, 2));
const double crhs39 = crhs38*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs17 - crhs2 + crhs33 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs37*(acceleration_alpha_method(0,0)*crhs36 + crhs18 + crhs3*crhs35)) + N[1]*(acceleration_alpha_method(1,0) + crhs37*(acceleration_alpha_method(1,0)*crhs36 + crhs29 + crhs35*crhs5)) + N[2]*(acceleration_alpha_method(2,0) + crhs37*(acceleration_alpha_method(2,0)*crhs36 + crhs30 + crhs35*crhs7))));
const double crhs40 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs41 = crhs10*crhs40;
const double crhs42 = rho*(DN(0,0)*crhs12 + DN(0,1)*crhs13);
const double crhs43 = rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)));
const double crhs44 = v(0,1) - vn(0,1);
const double crhs45 = crhs1*crhs44 + vn(0,1);
const double crhs46 = v(1,1) - vn(1,1);
const double crhs47 = crhs1*crhs46 + vn(1,1);
const double crhs48 = v(2,1) - vn(2,1);
const double crhs49 = crhs1*crhs48 + vn(2,1);
const double crhs50 = N[0]*crhs45 + N[1]*crhs47 + N[2]*crhs49;
const double crhs51 = rho*(crhs12*(DN(0,0)*crhs45 + DN(1,0)*crhs47 + DN(2,0)*crhs49) + crhs13*(DN(0,1)*crhs45 + DN(1,1)*crhs47 + DN(2,1)*crhs49));
const double crhs52 = -acceleration_alpha_method(0,1);
const double crhs53 = -acceleration_alpha_method(1,1);
const double crhs54 = -acceleration_alpha_method(2,1);
const double crhs55 = N[0]*(acceleration_alpha_method(0,1) - crhs28*(-acceleration_alpha_method(0,1)*crhs26 + crhs24*crhs44 + crhs52)) + N[1]*(acceleration_alpha_method(1,1) - crhs28*(-acceleration_alpha_method(1,1)*crhs26 + crhs24*crhs46 + crhs53)) + N[2]*(acceleration_alpha_method(2,1) - crhs28*(-acceleration_alpha_method(2,1)*crhs26 + crhs24*crhs48 + crhs54));
const double crhs56 = crhs32*crhs50;
const double crhs57 = crhs38*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs43 + crhs51 + crhs56 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs37*(acceleration_alpha_method(0,1)*crhs36 + crhs35*crhs44 + crhs52)) + N[1]*(acceleration_alpha_method(1,1) + crhs37*(acceleration_alpha_method(1,1)*crhs36 + crhs35*crhs46 + crhs53)) + N[2]*(acceleration_alpha_method(2,1) + crhs37*(acceleration_alpha_method(2,1)*crhs36 + crhs35*crhs48 + crhs54))));
const double crhs58 = -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + crhs15) + volume_error_ratio;
const double crhs59 = N[1]*rho;
const double crhs60 = crhs40*crhs59;
const double crhs61 = rho*(DN(1,0)*crhs12 + DN(1,1)*crhs13);
const double crhs62 = N[2]*rho;
const double crhs63 = crhs40*crhs62;
const double crhs64 = rho*(DN(2,0)*crhs12 + DN(2,1)*crhs13);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs16 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs17 + N[0]*crhs2 - crhs10*crhs31 - crhs11*crhs9 - crhs39*crhs41 - crhs39*crhs42;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs16 - DN(0,1)*stress[1] + N[0]*crhs43 - N[0]*crhs51 - crhs10*crhs55 - crhs11*crhs50 - crhs41*crhs57 - crhs42*crhs57;
rhs[2]=-DN(0,0)*crhs39 - DN(0,1)*crhs57 + N[0]*crhs58;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs16 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs17 + N[1]*crhs2 - N[1]*crhs33 - crhs31*crhs59 - crhs39*crhs60 - crhs39*crhs61;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs16 - DN(1,1)*stress[1] + N[1]*crhs43 - N[1]*crhs51 - N[1]*crhs56 - crhs55*crhs59 - crhs57*crhs60 - crhs57*crhs61;
rhs[5]=-DN(1,0)*crhs39 - DN(1,1)*crhs57 + N[1]*crhs58;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs16 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs17 + N[2]*crhs2 - N[2]*crhs33 - crhs31*crhs62 - crhs39*crhs63 - crhs39*crhs64;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs16 - DN(2,1)*stress[1] + N[2]*crhs43 - N[2]*crhs51 - N[2]*crhs56 - crhs55*crhs62 - crhs57*crhs63 - crhs57*crhs64;
rhs[8]=-DN(2,0)*crhs39 - DN(2,1)*crhs57 + N[2]*crhs58;


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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = 1.0/(max_spectral_radius + 1);
const double crhs2 = rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs1*(f(3,0) - fn(3,0)) + fn(3,0)));
const double crhs3 = v(0,0) - vn(0,0);
const double crhs4 = crhs1*crhs3 + vn(0,0);
const double crhs5 = v(1,0) - vn(1,0);
const double crhs6 = crhs1*crhs5 + vn(1,0);
const double crhs7 = v(2,0) - vn(2,0);
const double crhs8 = crhs1*crhs7 + vn(2,0);
const double crhs9 = v(3,0) - vn(3,0);
const double crhs10 = crhs1*crhs9 + vn(3,0);
const double crhs11 = N[0]*crhs4 + N[1]*crhs6 + N[2]*crhs8 + N[3]*crhs10;
const double crhs12 = N[0]*rho;
const double crhs13 = crhs12*volume_error_ratio;
const double crhs14 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs15 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs16 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs17 = rho*stab_c2*sqrt(pow(crhs14, 2) + pow(crhs15, 2) + pow(crhs16, 2));
const double crhs18 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs19 = (crhs18 - volume_error_ratio)*(crhs17*h/stab_c1 + mu);
const double crhs20 = rho*(crhs14*(DN(0,0)*crhs4 + DN(1,0)*crhs6 + DN(2,0)*crhs8 + DN(3,0)*crhs10) + crhs15*(DN(0,1)*crhs4 + DN(1,1)*crhs6 + DN(2,1)*crhs8 + DN(3,1)*crhs10) + crhs16*(DN(0,2)*crhs4 + DN(1,2)*crhs6 + DN(2,2)*crhs8 + DN(3,2)*crhs10));
const double crhs21 = -acceleration_alpha_method(0,0);
const double crhs22 = 1.0/dt;
const double crhs23 = 0.5*max_spectral_radius;
const double crhs24 = -crhs1;
const double crhs25 = crhs24 + 0.5;
const double crhs26 = 1.0/(-crhs1*(crhs23 - 1.5) + crhs25);
const double crhs27 = crhs22*crhs26;
const double crhs28 = crhs1*(1.5 - crhs23);
const double crhs29 = crhs26*(crhs1 - crhs28 + 0.5);
const double crhs30 = 0.5*crhs1;
const double crhs31 = crhs30*(max_spectral_radius - 3);
const double crhs32 = -acceleration_alpha_method(1,0);
const double crhs33 = -acceleration_alpha_method(2,0);
const double crhs34 = -acceleration_alpha_method(3,0);
const double crhs35 = N[0]*(acceleration_alpha_method(0,0) - crhs31*(-acceleration_alpha_method(0,0)*crhs29 + crhs21 + crhs27*crhs3)) + N[1]*(acceleration_alpha_method(1,0) - crhs31*(-acceleration_alpha_method(1,0)*crhs29 + crhs27*crhs5 + crhs32)) + N[2]*(acceleration_alpha_method(2,0) - crhs31*(-acceleration_alpha_method(2,0)*crhs29 + crhs27*crhs7 + crhs33)) + N[3]*(acceleration_alpha_method(3,0) - crhs31*(-acceleration_alpha_method(3,0)*crhs29 + crhs27*crhs9 + crhs34));
const double crhs36 = rho*volume_error_ratio;
const double crhs37 = crhs11*crhs36;
const double crhs38 = 1.0/(crhs25 + crhs28);
const double crhs39 = crhs22*crhs38;
const double crhs40 = crhs38*(crhs24 + crhs28 - 0.5);
const double crhs41 = crhs30*(3 - max_spectral_radius);
const double crhs42 = 1.0/(crhs17/h + crhs22*dyn_tau*rho + mu*stab_c1/pow(h, 2));
const double crhs43 = crhs42*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs2 + crhs20 + crhs37 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs41*(acceleration_alpha_method(0,0)*crhs40 + crhs21 + crhs3*crhs39)) + N[1]*(acceleration_alpha_method(1,0) + crhs41*(acceleration_alpha_method(1,0)*crhs40 + crhs32 + crhs39*crhs5)) + N[2]*(acceleration_alpha_method(2,0) + crhs41*(acceleration_alpha_method(2,0)*crhs40 + crhs33 + crhs39*crhs7)) + N[3]*(acceleration_alpha_method(3,0) + crhs41*(acceleration_alpha_method(3,0)*crhs40 + crhs34 + crhs39*crhs9))));
const double crhs44 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs45 = crhs12*crhs44;
const double crhs46 = rho*(DN(0,0)*crhs14 + DN(0,1)*crhs15 + DN(0,2)*crhs16);
const double crhs47 = rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs1*(f(3,1) - fn(3,1)) + fn(3,1)));
const double crhs48 = v(0,1) - vn(0,1);
const double crhs49 = crhs1*crhs48 + vn(0,1);
const double crhs50 = v(1,1) - vn(1,1);
const double crhs51 = crhs1*crhs50 + vn(1,1);
const double crhs52 = v(2,1) - vn(2,1);
const double crhs53 = crhs1*crhs52 + vn(2,1);
const double crhs54 = v(3,1) - vn(3,1);
const double crhs55 = crhs1*crhs54 + vn(3,1);
const double crhs56 = N[0]*crhs49 + N[1]*crhs51 + N[2]*crhs53 + N[3]*crhs55;
const double crhs57 = rho*(crhs14*(DN(0,0)*crhs49 + DN(1,0)*crhs51 + DN(2,0)*crhs53 + DN(3,0)*crhs55) + crhs15*(DN(0,1)*crhs49 + DN(1,1)*crhs51 + DN(2,1)*crhs53 + DN(3,1)*crhs55) + crhs16*(DN(0,2)*crhs49 + DN(1,2)*crhs51 + DN(2,2)*crhs53 + DN(3,2)*crhs55));
const double crhs58 = -acceleration_alpha_method(0,1);
const double crhs59 = -acceleration_alpha_method(1,1);
const double crhs60 = -acceleration_alpha_method(2,1);
const double crhs61 = -acceleration_alpha_method(3,1);
const double crhs62 = N[0]*(acceleration_alpha_method(0,1) - crhs31*(-acceleration_alpha_method(0,1)*crhs29 + crhs27*crhs48 + crhs58)) + N[1]*(acceleration_alpha_method(1,1) - crhs31*(-acceleration_alpha_method(1,1)*crhs29 + crhs27*crhs50 + crhs59)) + N[2]*(acceleration_alpha_method(2,1) - crhs31*(-acceleration_alpha_method(2,1)*crhs29 + crhs27*crhs52 + crhs60)) + N[3]*(acceleration_alpha_method(3,1) - crhs31*(-acceleration_alpha_method(3,1)*crhs29 + crhs27*crhs54 + crhs61));
const double crhs63 = crhs36*crhs56;
const double crhs64 = crhs42*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs47 + crhs57 + crhs63 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs41*(acceleration_alpha_method(0,1)*crhs40 + crhs39*crhs48 + crhs58)) + N[1]*(acceleration_alpha_method(1,1) + crhs41*(acceleration_alpha_method(1,1)*crhs40 + crhs39*crhs50 + crhs59)) + N[2]*(acceleration_alpha_method(2,1) + crhs41*(acceleration_alpha_method(2,1)*crhs40 + crhs39*crhs52 + crhs60)) + N[3]*(acceleration_alpha_method(3,1) + crhs41*(acceleration_alpha_method(3,1)*crhs40 + crhs39*crhs54 + crhs61))));
const double crhs65 = rho*(N[0]*(crhs1*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs1*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs1*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs1*(f(3,2) - fn(3,2)) + fn(3,2)));
const double crhs66 = v(0,2) - vn(0,2);
const double crhs67 = crhs1*crhs66 + vn(0,2);
const double crhs68 = v(1,2) - vn(1,2);
const double crhs69 = crhs1*crhs68 + vn(1,2);
const double crhs70 = v(2,2) - vn(2,2);
const double crhs71 = crhs1*crhs70 + vn(2,2);
const double crhs72 = v(3,2) - vn(3,2);
const double crhs73 = crhs1*crhs72 + vn(3,2);
const double crhs74 = N[0]*crhs67 + N[1]*crhs69 + N[2]*crhs71 + N[3]*crhs73;
const double crhs75 = rho*(crhs14*(DN(0,0)*crhs67 + DN(1,0)*crhs69 + DN(2,0)*crhs71 + DN(3,0)*crhs73) + crhs15*(DN(0,1)*crhs67 + DN(1,1)*crhs69 + DN(2,1)*crhs71 + DN(3,1)*crhs73) + crhs16*(DN(0,2)*crhs67 + DN(1,2)*crhs69 + DN(2,2)*crhs71 + DN(3,2)*crhs73));
const double crhs76 = -acceleration_alpha_method(0,2);
const double crhs77 = -acceleration_alpha_method(1,2);
const double crhs78 = -acceleration_alpha_method(2,2);
const double crhs79 = -acceleration_alpha_method(3,2);
const double crhs80 = N[0]*(acceleration_alpha_method(0,2) - crhs31*(-acceleration_alpha_method(0,2)*crhs29 + crhs27*crhs66 + crhs76)) + N[1]*(acceleration_alpha_method(1,2) - crhs31*(-acceleration_alpha_method(1,2)*crhs29 + crhs27*crhs68 + crhs77)) + N[2]*(acceleration_alpha_method(2,2) - crhs31*(-acceleration_alpha_method(2,2)*crhs29 + crhs27*crhs70 + crhs78)) + N[3]*(acceleration_alpha_method(3,2) - crhs31*(-acceleration_alpha_method(3,2)*crhs29 + crhs27*crhs72 + crhs79));
const double crhs81 = crhs36*crhs74;
const double crhs82 = crhs42*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs65 + crhs75 + crhs81 + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs41*(acceleration_alpha_method(0,2)*crhs40 + crhs39*crhs66 + crhs76)) + N[1]*(acceleration_alpha_method(1,2) + crhs41*(acceleration_alpha_method(1,2)*crhs40 + crhs39*crhs68 + crhs77)) + N[2]*(acceleration_alpha_method(2,2) + crhs41*(acceleration_alpha_method(2,2)*crhs40 + crhs39*crhs70 + crhs78)) + N[3]*(acceleration_alpha_method(3,2) + crhs41*(acceleration_alpha_method(3,2)*crhs40 + crhs39*crhs72 + crhs79))));
const double crhs83 = -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + crhs18) + volume_error_ratio;
const double crhs84 = N[1]*rho;
const double crhs85 = crhs44*crhs84;
const double crhs86 = rho*(DN(1,0)*crhs14 + DN(1,1)*crhs15 + DN(1,2)*crhs16);
const double crhs87 = N[2]*rho;
const double crhs88 = crhs44*crhs87;
const double crhs89 = rho*(DN(2,0)*crhs14 + DN(2,1)*crhs15 + DN(2,2)*crhs16);
const double crhs90 = N[3]*rho;
const double crhs91 = crhs44*crhs90;
const double crhs92 = rho*(DN(3,0)*crhs14 + DN(3,1)*crhs15 + DN(3,2)*crhs16);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs19 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs2 - N[0]*crhs20 - crhs11*crhs13 - crhs12*crhs35 - crhs43*crhs45 - crhs43*crhs46;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs19 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs47 - N[0]*crhs57 - crhs12*crhs62 - crhs13*crhs56 - crhs45*crhs64 - crhs46*crhs64;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs19 - DN(0,2)*stress[2] + N[0]*crhs65 - N[0]*crhs75 - crhs12*crhs80 - crhs13*crhs74 - crhs45*crhs82 - crhs46*crhs82;
rhs[3]=-DN(0,0)*crhs43 - DN(0,1)*crhs64 - DN(0,2)*crhs82 + N[0]*crhs83;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs19 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs2 - N[1]*crhs20 - N[1]*crhs37 - crhs35*crhs84 - crhs43*crhs85 - crhs43*crhs86;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs19 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs47 - N[1]*crhs57 - N[1]*crhs63 - crhs62*crhs84 - crhs64*crhs85 - crhs64*crhs86;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs19 - DN(1,2)*stress[2] + N[1]*crhs65 - N[1]*crhs75 - N[1]*crhs81 - crhs80*crhs84 - crhs82*crhs85 - crhs82*crhs86;
rhs[7]=-DN(1,0)*crhs43 - DN(1,1)*crhs64 - DN(1,2)*crhs82 + N[1]*crhs83;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs19 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs2 - N[2]*crhs20 - N[2]*crhs37 - crhs35*crhs87 - crhs43*crhs88 - crhs43*crhs89;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs19 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs47 - N[2]*crhs57 - N[2]*crhs63 - crhs62*crhs87 - crhs64*crhs88 - crhs64*crhs89;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs19 - DN(2,2)*stress[2] + N[2]*crhs65 - N[2]*crhs75 - N[2]*crhs81 - crhs80*crhs87 - crhs82*crhs88 - crhs82*crhs89;
rhs[11]=-DN(2,0)*crhs43 - DN(2,1)*crhs64 - DN(2,2)*crhs82 + N[2]*crhs83;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs19 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs2 - N[3]*crhs20 - N[3]*crhs37 - crhs35*crhs90 - crhs43*crhs91 - crhs43*crhs92;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs19 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs47 - N[3]*crhs57 - N[3]*crhs63 - crhs62*crhs90 - crhs64*crhs91 - crhs64*crhs92;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs19 - DN(3,2)*stress[2] + N[3]*crhs65 - N[3]*crhs75 - N[3]*crhs81 - crhs80*crhs90 - crhs82*crhs91 - crhs82*crhs92;
rhs[15]=-DN(3,0)*crhs43 - DN(3,1)*crhs64 - DN(3,2)*crhs82 + N[3]*crhs83;


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

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 = 1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 = cV2*rho;
const double cV4 = cV3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV5 = N[0]*cV4;
const double cV6 = cV3*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV7 = N[1]*cV4;
const double cV8 = cV3*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV9 = N[2]*cV4;
const double cV10 = cV3*(DN(2,0)*cV0 + DN(2,1)*cV1);
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


    const double cH0 = 1.0/(max_spectral_radius + 1);
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH3 = 1.0/dt;
const double cH4 = 0.5*cH3*(max_spectral_radius - 3)/(-cH0*(0.5*max_spectral_radius - 1.5) - cH0 + 0.5);
const double cH5 = 1.0/(cH3*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2));
const double cH6 = cH5*rho;
const double cH7 = cH6*(DN(0,0)*cH1 + DN(0,1)*cH2 - N[0]*cH4 + N[0]*volume_error_ratio);
const double cH8 = cH6*(DN(1,0)*cH1 + DN(1,1)*cH2 - N[1]*cH4 + N[1]*volume_error_ratio);
const double cH9 = cH6*(DN(2,0)*cH1 + DN(2,1)*cH2 - N[2]*cH4 + N[2]*volume_error_ratio);
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


    const double cKee0 = 1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
Kee(1,2)=cKee3;
Kee(2,0)=cKee2;
Kee(2,1)=cKee3;
Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 = 1.0/(max_spectral_radius + 1);
const double crhs_ee1 = -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1)) + volume_error_ratio;
const double crhs_ee2 = v(0,0) - vn(0,0);
const double crhs_ee3 = crhs_ee0*crhs_ee2 + vn(0,0);
const double crhs_ee4 = v(1,0) - vn(1,0);
const double crhs_ee5 = crhs_ee0*crhs_ee4 + vn(1,0);
const double crhs_ee6 = v(2,0) - vn(2,0);
const double crhs_ee7 = crhs_ee0*crhs_ee6 + vn(2,0);
const double crhs_ee8 = rho*volume_error_ratio;
const double crhs_ee9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee11 = 1.0/dt;
const double crhs_ee12 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee13 = 1.0/(crhs_ee12 + 0.5);
const double crhs_ee14 = crhs_ee11*crhs_ee13;
const double crhs_ee15 = crhs_ee13*(crhs_ee12 - 0.5);
const double crhs_ee16 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee17 = 1.0/(crhs_ee11*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee18 = crhs_ee17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + crhs_ee8*(N[0]*crhs_ee3 + N[1]*crhs_ee5 + N[2]*crhs_ee7) - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee16*(acceleration_alpha_method(0,0)*crhs_ee15 - acceleration_alpha_method(0,0) + crhs_ee14*crhs_ee2)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee16*(acceleration_alpha_method(1,0)*crhs_ee15 - acceleration_alpha_method(1,0) + crhs_ee14*crhs_ee4)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee16*(acceleration_alpha_method(2,0)*crhs_ee15 - acceleration_alpha_method(2,0) + crhs_ee14*crhs_ee6)) + crhs_ee10*(DN(0,1)*crhs_ee3 + DN(1,1)*crhs_ee5 + DN(2,1)*crhs_ee7) + crhs_ee9*(DN(0,0)*crhs_ee3 + DN(1,0)*crhs_ee5 + DN(2,0)*crhs_ee7)));
const double crhs_ee19 = v(0,1) - vn(0,1);
const double crhs_ee20 = crhs_ee0*crhs_ee19 + vn(0,1);
const double crhs_ee21 = v(1,1) - vn(1,1);
const double crhs_ee22 = crhs_ee0*crhs_ee21 + vn(1,1);
const double crhs_ee23 = v(2,1) - vn(2,1);
const double crhs_ee24 = crhs_ee0*crhs_ee23 + vn(2,1);
const double crhs_ee25 = crhs_ee17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + crhs_ee8*(N[0]*crhs_ee20 + N[1]*crhs_ee22 + N[2]*crhs_ee24) - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee16*(acceleration_alpha_method(0,1)*crhs_ee15 - acceleration_alpha_method(0,1) + crhs_ee14*crhs_ee19)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee16*(acceleration_alpha_method(1,1)*crhs_ee15 - acceleration_alpha_method(1,1) + crhs_ee14*crhs_ee21)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee16*(acceleration_alpha_method(2,1)*crhs_ee15 - acceleration_alpha_method(2,1) + crhs_ee14*crhs_ee23)) + crhs_ee10*(DN(0,1)*crhs_ee20 + DN(1,1)*crhs_ee22 + DN(2,1)*crhs_ee24) + crhs_ee9*(DN(0,0)*crhs_ee20 + DN(1,0)*crhs_ee22 + DN(2,0)*crhs_ee24)));
rhs_ee[0]=-DNenr(0,0)*crhs_ee18 - DNenr(0,1)*crhs_ee25 + Nenr[0]*crhs_ee1;
rhs_ee[1]=-DNenr(1,0)*crhs_ee18 - DNenr(1,1)*crhs_ee25 + Nenr[1]*crhs_ee1;
rhs_ee[2]=-DNenr(2,0)*crhs_ee18 - DNenr(2,1)*crhs_ee25 + Nenr[2]*crhs_ee1;


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

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 = 1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 = cV3*rho;
const double cV5 = cV4*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV6 = N[0]*cV5;
const double cV7 = cV4*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV8 = N[1]*cV5;
const double cV9 = cV4*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV10 = N[2]*cV5;
const double cV11 = cV4*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV12 = N[3]*cV5;
const double cV13 = cV4*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
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


    const double cH0 = 1.0/(max_spectral_radius + 1);
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH3 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH4 = 1.0/dt;
const double cH5 = 0.5*cH4*(max_spectral_radius - 3)/(-cH0*(0.5*max_spectral_radius - 1.5) - cH0 + 0.5);
const double cH6 = 1.0/(cH4*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 = cH6*rho;
const double cH8 = cH7*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 - N[0]*cH5 + N[0]*volume_error_ratio);
const double cH9 = cH7*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 - N[1]*cH5 + N[1]*volume_error_ratio);
const double cH10 = cH7*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 - N[2]*cH5 + N[2]*volume_error_ratio);
const double cH11 = cH7*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 - N[3]*cH5 + N[3]*volume_error_ratio);
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


    const double cKee0 = 1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 = cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 = cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 = cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
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


    const double crhs_ee0 = 1.0/(max_spectral_radius + 1);
const double crhs_ee1 = -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(0,2)*v(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(1,2)*v(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(2,2)*v(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,0)*v(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,1)*v(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + DN(3,2)*v(3,2)) + volume_error_ratio;
const double crhs_ee2 = v(0,0) - vn(0,0);
const double crhs_ee3 = crhs_ee0*crhs_ee2 + vn(0,0);
const double crhs_ee4 = v(1,0) - vn(1,0);
const double crhs_ee5 = crhs_ee0*crhs_ee4 + vn(1,0);
const double crhs_ee6 = v(2,0) - vn(2,0);
const double crhs_ee7 = crhs_ee0*crhs_ee6 + vn(2,0);
const double crhs_ee8 = v(3,0) - vn(3,0);
const double crhs_ee9 = crhs_ee0*crhs_ee8 + vn(3,0);
const double crhs_ee10 = rho*volume_error_ratio;
const double crhs_ee11 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee12 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee13 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee14 = 1.0/dt;
const double crhs_ee15 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee16 = 1.0/(crhs_ee15 + 0.5);
const double crhs_ee17 = crhs_ee14*crhs_ee16;
const double crhs_ee18 = crhs_ee16*(crhs_ee15 - 0.5);
const double crhs_ee19 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee20 = 1.0/(crhs_ee14*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee11, 2) + pow(crhs_ee12, 2) + pow(crhs_ee13, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee21 = crhs_ee20*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + crhs_ee10*(N[0]*crhs_ee3 + N[1]*crhs_ee5 + N[2]*crhs_ee7 + N[3]*crhs_ee9) - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs_ee0*(f(3,0) - fn(3,0)) + fn(3,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee19*(acceleration_alpha_method(0,0)*crhs_ee18 - acceleration_alpha_method(0,0) + crhs_ee17*crhs_ee2)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee19*(acceleration_alpha_method(1,0)*crhs_ee18 - acceleration_alpha_method(1,0) + crhs_ee17*crhs_ee4)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee19*(acceleration_alpha_method(2,0)*crhs_ee18 - acceleration_alpha_method(2,0) + crhs_ee17*crhs_ee6)) + N[3]*(acceleration_alpha_method(3,0) + crhs_ee19*(acceleration_alpha_method(3,0)*crhs_ee18 - acceleration_alpha_method(3,0) + crhs_ee17*crhs_ee8)) + crhs_ee11*(DN(0,0)*crhs_ee3 + DN(1,0)*crhs_ee5 + DN(2,0)*crhs_ee7 + DN(3,0)*crhs_ee9) + crhs_ee12*(DN(0,1)*crhs_ee3 + DN(1,1)*crhs_ee5 + DN(2,1)*crhs_ee7 + DN(3,1)*crhs_ee9) + crhs_ee13*(DN(0,2)*crhs_ee3 + DN(1,2)*crhs_ee5 + DN(2,2)*crhs_ee7 + DN(3,2)*crhs_ee9)));
const double crhs_ee22 = v(0,1) - vn(0,1);
const double crhs_ee23 = crhs_ee0*crhs_ee22 + vn(0,1);
const double crhs_ee24 = v(1,1) - vn(1,1);
const double crhs_ee25 = crhs_ee0*crhs_ee24 + vn(1,1);
const double crhs_ee26 = v(2,1) - vn(2,1);
const double crhs_ee27 = crhs_ee0*crhs_ee26 + vn(2,1);
const double crhs_ee28 = v(3,1) - vn(3,1);
const double crhs_ee29 = crhs_ee0*crhs_ee28 + vn(3,1);
const double crhs_ee30 = crhs_ee20*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + crhs_ee10*(N[0]*crhs_ee23 + N[1]*crhs_ee25 + N[2]*crhs_ee27 + N[3]*crhs_ee29) - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs_ee0*(f(3,1) - fn(3,1)) + fn(3,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee19*(acceleration_alpha_method(0,1)*crhs_ee18 - acceleration_alpha_method(0,1) + crhs_ee17*crhs_ee22)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee19*(acceleration_alpha_method(1,1)*crhs_ee18 - acceleration_alpha_method(1,1) + crhs_ee17*crhs_ee24)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee19*(acceleration_alpha_method(2,1)*crhs_ee18 - acceleration_alpha_method(2,1) + crhs_ee17*crhs_ee26)) + N[3]*(acceleration_alpha_method(3,1) + crhs_ee19*(acceleration_alpha_method(3,1)*crhs_ee18 - acceleration_alpha_method(3,1) + crhs_ee17*crhs_ee28)) + crhs_ee11*(DN(0,0)*crhs_ee23 + DN(1,0)*crhs_ee25 + DN(2,0)*crhs_ee27 + DN(3,0)*crhs_ee29) + crhs_ee12*(DN(0,1)*crhs_ee23 + DN(1,1)*crhs_ee25 + DN(2,1)*crhs_ee27 + DN(3,1)*crhs_ee29) + crhs_ee13*(DN(0,2)*crhs_ee23 + DN(1,2)*crhs_ee25 + DN(2,2)*crhs_ee27 + DN(3,2)*crhs_ee29)));
const double crhs_ee31 = v(0,2) - vn(0,2);
const double crhs_ee32 = crhs_ee0*crhs_ee31 + vn(0,2);
const double crhs_ee33 = v(1,2) - vn(1,2);
const double crhs_ee34 = crhs_ee0*crhs_ee33 + vn(1,2);
const double crhs_ee35 = v(2,2) - vn(2,2);
const double crhs_ee36 = crhs_ee0*crhs_ee35 + vn(2,2);
const double crhs_ee37 = v(3,2) - vn(3,2);
const double crhs_ee38 = crhs_ee0*crhs_ee37 + vn(3,2);
const double crhs_ee39 = crhs_ee20*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + crhs_ee10*(N[0]*crhs_ee32 + N[1]*crhs_ee34 + N[2]*crhs_ee36 + N[3]*crhs_ee38) - rho*(N[0]*(crhs_ee0*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs_ee0*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs_ee0*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs_ee0*(f(3,2) - fn(3,2)) + fn(3,2))) + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs_ee19*(acceleration_alpha_method(0,2)*crhs_ee18 - acceleration_alpha_method(0,2) + crhs_ee17*crhs_ee31)) + N[1]*(acceleration_alpha_method(1,2) + crhs_ee19*(acceleration_alpha_method(1,2)*crhs_ee18 - acceleration_alpha_method(1,2) + crhs_ee17*crhs_ee33)) + N[2]*(acceleration_alpha_method(2,2) + crhs_ee19*(acceleration_alpha_method(2,2)*crhs_ee18 - acceleration_alpha_method(2,2) + crhs_ee17*crhs_ee35)) + N[3]*(acceleration_alpha_method(3,2) + crhs_ee19*(acceleration_alpha_method(3,2)*crhs_ee18 - acceleration_alpha_method(3,2) + crhs_ee17*crhs_ee37)) + crhs_ee11*(DN(0,0)*crhs_ee32 + DN(1,0)*crhs_ee34 + DN(2,0)*crhs_ee36 + DN(3,0)*crhs_ee38) + crhs_ee12*(DN(0,1)*crhs_ee32 + DN(1,1)*crhs_ee34 + DN(2,1)*crhs_ee36 + DN(3,1)*crhs_ee38) + crhs_ee13*(DN(0,2)*crhs_ee32 + DN(1,2)*crhs_ee34 + DN(2,2)*crhs_ee36 + DN(3,2)*crhs_ee38)));
rhs_ee[0]=-DNenr(0,0)*crhs_ee21 - DNenr(0,1)*crhs_ee30 - DNenr(0,2)*crhs_ee39 + Nenr[0]*crhs_ee1;
rhs_ee[1]=-DNenr(1,0)*crhs_ee21 - DNenr(1,1)*crhs_ee30 - DNenr(1,2)*crhs_ee39 + Nenr[1]*crhs_ee1;
rhs_ee[2]=-DNenr(2,0)*crhs_ee21 - DNenr(2,1)*crhs_ee30 - DNenr(2,2)*crhs_ee39 + Nenr[2]*crhs_ee1;
rhs_ee[3]=-DNenr(3,0)*crhs_ee21 - DNenr(3,1)*crhs_ee30 - DNenr(3,2)*crhs_ee39 + Nenr[3]*crhs_ee1;


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
    const double volume_error_ratio = rData.VolumeError;
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
