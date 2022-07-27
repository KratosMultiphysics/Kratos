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
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSOutletPenaltyContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rLHS)
{
    double c_penalty = 0.0;
    const double mu = rData.EffectiveViscosity;
    const double max_spectral_radius = rData.MaxSprectraRadius;
    double kappa = c_penalty * mu;
    double gauss_weight = rData.Weight;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;

        const double crLHS0 = gauss_weight*kappa/(max_spectral_radius + 1);
    const double crLHS1 = -pow(DN(0,0), 2)*crLHS0;
    const double crLHS2 = DN(0,0)*crLHS0;
    const double crLHS3 = -DN(1,0)*crLHS2;
    const double crLHS4 = -DN(2,0)*crLHS2;
    const double crLHS5 = -pow(DN(1,0), 2)*crLHS0;
    const double crLHS6 = -DN(1,0)*DN(2,0)*crLHS0;
    const double crLHS7 = -pow(DN(2,0), 2)*crLHS0;
    rLHS(0,0)+=crLHS1;
    rLHS(0,1)+=0;
    rLHS(0,2)+=0;
    rLHS(0,3)+=crLHS3;
    rLHS(0,4)+=0;
    rLHS(0,5)+=0;
    rLHS(0,6)+=crLHS4;
    rLHS(0,7)+=0;
    rLHS(0,8)+=0;
    rLHS(1,0)+=0;
    rLHS(1,1)+=crLHS1;
    rLHS(1,2)+=0;
    rLHS(1,3)+=0;
    rLHS(1,4)+=crLHS3;
    rLHS(1,5)+=0;
    rLHS(1,6)+=0;
    rLHS(1,7)+=crLHS4;
    rLHS(1,8)+=0;
    rLHS(2,0)+=0;
    rLHS(2,1)+=0;
    rLHS(2,2)+=0;
    rLHS(2,3)+=0;
    rLHS(2,4)+=0;
    rLHS(2,5)+=0;
    rLHS(2,6)+=0;
    rLHS(2,7)+=0;
    rLHS(2,8)+=0;
    rLHS(3,0)+=crLHS3;
    rLHS(3,1)+=0;
    rLHS(3,2)+=0;
    rLHS(3,3)+=crLHS5;
    rLHS(3,4)+=0;
    rLHS(3,5)+=0;
    rLHS(3,6)+=crLHS6;
    rLHS(3,7)+=0;
    rLHS(3,8)+=0;
    rLHS(4,0)+=0;
    rLHS(4,1)+=crLHS3;
    rLHS(4,2)+=0;
    rLHS(4,3)+=0;
    rLHS(4,4)+=crLHS5;
    rLHS(4,5)+=0;
    rLHS(4,6)+=0;
    rLHS(4,7)+=crLHS6;
    rLHS(4,8)+=0;
    rLHS(5,0)+=0;
    rLHS(5,1)+=0;
    rLHS(5,2)+=0;
    rLHS(5,3)+=0;
    rLHS(5,4)+=0;
    rLHS(5,5)+=0;
    rLHS(5,6)+=0;
    rLHS(5,7)+=0;
    rLHS(5,8)+=0;
    rLHS(6,0)+=crLHS4;
    rLHS(6,1)+=0;
    rLHS(6,2)+=0;
    rLHS(6,3)+=crLHS6;
    rLHS(6,4)+=0;
    rLHS(6,5)+=0;
    rLHS(6,6)+=crLHS7;
    rLHS(6,7)+=0;
    rLHS(6,8)+=0;
    rLHS(7,0)+=0;
    rLHS(7,1)+=crLHS4;
    rLHS(7,2)+=0;
    rLHS(7,3)+=0;
    rLHS(7,4)+=crLHS6;
    rLHS(7,5)+=0;
    rLHS(7,6)+=0;
    rLHS(7,7)+=crLHS7;
    rLHS(7,8)+=0;
    rLHS(8,0)+=0;
    rLHS(8,1)+=0;
    rLHS(8,2)+=0;
    rLHS(8,3)+=0;
    rLHS(8,4)+=0;
    rLHS(8,5)+=0;
    rLHS(8,6)+=0;
    rLHS(8,7)+=0;
    rLHS(8,8)+=0;

}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSOutletPenaltyContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rLHS)
{
    double c_penalty = 10.0;
    const double mu = rData.EffectiveViscosity;
    const double max_spectral_radius = rData.MaxSprectraRadius;
    double kappa = c_penalty * mu;
    double gauss_weight = rData.Weight;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;

        const double crLHS0 = gauss_weight*kappa/(max_spectral_radius + 1);
    const double crLHS1 = -pow(DN(0,0), 2)*crLHS0;
    const double crLHS2 = DN(0,0)*crLHS0;
    const double crLHS3 = -DN(1,0)*crLHS2;
    const double crLHS4 = -DN(2,0)*crLHS2;
    const double crLHS5 = -DN(3,0)*crLHS2;
    const double crLHS6 = -pow(DN(1,0), 2)*crLHS0;
    const double crLHS7 = DN(1,0)*crLHS0;
    const double crLHS8 = -DN(2,0)*crLHS7;
    const double crLHS9 = -DN(3,0)*crLHS7;
    const double crLHS10 = -pow(DN(2,0), 2)*crLHS0;
    const double crLHS11 = -DN(2,0)*DN(3,0)*crLHS0;
    const double crLHS12 = -pow(DN(3,0), 2)*crLHS0;
    rLHS(0,0)+=crLHS1;
    rLHS(0,1)+=0;
    rLHS(0,2)+=0;
    rLHS(0,3)+=0;
    rLHS(0,4)+=crLHS3;
    rLHS(0,5)+=0;
    rLHS(0,6)+=0;
    rLHS(0,7)+=0;
    rLHS(0,8)+=crLHS4;
    rLHS(0,9)+=0;
    rLHS(0,10)+=0;
    rLHS(0,11)+=0;
    rLHS(0,12)+=crLHS5;
    rLHS(0,13)+=0;
    rLHS(0,14)+=0;
    rLHS(0,15)+=0;
    rLHS(1,0)+=0;
    rLHS(1,1)+=crLHS1;
    rLHS(1,2)+=0;
    rLHS(1,3)+=0;
    rLHS(1,4)+=0;
    rLHS(1,5)+=crLHS3;
    rLHS(1,6)+=0;
    rLHS(1,7)+=0;
    rLHS(1,8)+=0;
    rLHS(1,9)+=crLHS4;
    rLHS(1,10)+=0;
    rLHS(1,11)+=0;
    rLHS(1,12)+=0;
    rLHS(1,13)+=crLHS5;
    rLHS(1,14)+=0;
    rLHS(1,15)+=0;
    rLHS(2,0)+=0;
    rLHS(2,1)+=0;
    rLHS(2,2)+=crLHS1;
    rLHS(2,3)+=0;
    rLHS(2,4)+=0;
    rLHS(2,5)+=0;
    rLHS(2,6)+=crLHS3;
    rLHS(2,7)+=0;
    rLHS(2,8)+=0;
    rLHS(2,9)+=0;
    rLHS(2,10)+=crLHS4;
    rLHS(2,11)+=0;
    rLHS(2,12)+=0;
    rLHS(2,13)+=0;
    rLHS(2,14)+=crLHS5;
    rLHS(2,15)+=0;
    rLHS(3,0)+=0;
    rLHS(3,1)+=0;
    rLHS(3,2)+=0;
    rLHS(3,3)+=0;
    rLHS(3,4)+=0;
    rLHS(3,5)+=0;
    rLHS(3,6)+=0;
    rLHS(3,7)+=0;
    rLHS(3,8)+=0;
    rLHS(3,9)+=0;
    rLHS(3,10)+=0;
    rLHS(3,11)+=0;
    rLHS(3,12)+=0;
    rLHS(3,13)+=0;
    rLHS(3,14)+=0;
    rLHS(3,15)+=0;
    rLHS(4,0)+=crLHS3;
    rLHS(4,1)+=0;
    rLHS(4,2)+=0;
    rLHS(4,3)+=0;
    rLHS(4,4)+=crLHS6;
    rLHS(4,5)+=0;
    rLHS(4,6)+=0;
    rLHS(4,7)+=0;
    rLHS(4,8)+=crLHS8;
    rLHS(4,9)+=0;
    rLHS(4,10)+=0;
    rLHS(4,11)+=0;
    rLHS(4,12)+=crLHS9;
    rLHS(4,13)+=0;
    rLHS(4,14)+=0;
    rLHS(4,15)+=0;
    rLHS(5,0)+=0;
    rLHS(5,1)+=crLHS3;
    rLHS(5,2)+=0;
    rLHS(5,3)+=0;
    rLHS(5,4)+=0;
    rLHS(5,5)+=crLHS6;
    rLHS(5,6)+=0;
    rLHS(5,7)+=0;
    rLHS(5,8)+=0;
    rLHS(5,9)+=crLHS8;
    rLHS(5,10)+=0;
    rLHS(5,11)+=0;
    rLHS(5,12)+=0;
    rLHS(5,13)+=crLHS9;
    rLHS(5,14)+=0;
    rLHS(5,15)+=0;
    rLHS(6,0)+=0;
    rLHS(6,1)+=0;
    rLHS(6,2)+=crLHS3;
    rLHS(6,3)+=0;
    rLHS(6,4)+=0;
    rLHS(6,5)+=0;
    rLHS(6,6)+=crLHS6;
    rLHS(6,7)+=0;
    rLHS(6,8)+=0;
    rLHS(6,9)+=0;
    rLHS(6,10)+=crLHS8;
    rLHS(6,11)+=0;
    rLHS(6,12)+=0;
    rLHS(6,13)+=0;
    rLHS(6,14)+=crLHS9;
    rLHS(6,15)+=0;
    rLHS(7,0)+=0;
    rLHS(7,1)+=0;
    rLHS(7,2)+=0;
    rLHS(7,3)+=0;
    rLHS(7,4)+=0;
    rLHS(7,5)+=0;
    rLHS(7,6)+=0;
    rLHS(7,7)+=0;
    rLHS(7,8)+=0;
    rLHS(7,9)+=0;
    rLHS(7,10)+=0;
    rLHS(7,11)+=0;
    rLHS(7,12)+=0;
    rLHS(7,13)+=0;
    rLHS(7,14)+=0;
    rLHS(7,15)+=0;
    rLHS(8,0)+=crLHS4;
    rLHS(8,1)+=0;
    rLHS(8,2)+=0;
    rLHS(8,3)+=0;
    rLHS(8,4)+=crLHS8;
    rLHS(8,5)+=0;
    rLHS(8,6)+=0;
    rLHS(8,7)+=0;
    rLHS(8,8)+=crLHS10;
    rLHS(8,9)+=0;
    rLHS(8,10)+=0;
    rLHS(8,11)+=0;
    rLHS(8,12)+=crLHS11;
    rLHS(8,13)+=0;
    rLHS(8,14)+=0;
    rLHS(8,15)+=0;
    rLHS(9,0)+=0;
    rLHS(9,1)+=crLHS4;
    rLHS(9,2)+=0;
    rLHS(9,3)+=0;
    rLHS(9,4)+=0;
    rLHS(9,5)+=crLHS8;
    rLHS(9,6)+=0;
    rLHS(9,7)+=0;
    rLHS(9,8)+=0;
    rLHS(9,9)+=crLHS10;
    rLHS(9,10)+=0;
    rLHS(9,11)+=0;
    rLHS(9,12)+=0;
    rLHS(9,13)+=crLHS11;
    rLHS(9,14)+=0;
    rLHS(9,15)+=0;
    rLHS(10,0)+=0;
    rLHS(10,1)+=0;
    rLHS(10,2)+=crLHS4;
    rLHS(10,3)+=0;
    rLHS(10,4)+=0;
    rLHS(10,5)+=0;
    rLHS(10,6)+=crLHS8;
    rLHS(10,7)+=0;
    rLHS(10,8)+=0;
    rLHS(10,9)+=0;
    rLHS(10,10)+=crLHS10;
    rLHS(10,11)+=0;
    rLHS(10,12)+=0;
    rLHS(10,13)+=0;
    rLHS(10,14)+=crLHS11;
    rLHS(10,15)+=0;
    rLHS(11,0)+=0;
    rLHS(11,1)+=0;
    rLHS(11,2)+=0;
    rLHS(11,3)+=0;
    rLHS(11,4)+=0;
    rLHS(11,5)+=0;
    rLHS(11,6)+=0;
    rLHS(11,7)+=0;
    rLHS(11,8)+=0;
    rLHS(11,9)+=0;
    rLHS(11,10)+=0;
    rLHS(11,11)+=0;
    rLHS(11,12)+=0;
    rLHS(11,13)+=0;
    rLHS(11,14)+=0;
    rLHS(11,15)+=0;
    rLHS(12,0)+=crLHS5;
    rLHS(12,1)+=0;
    rLHS(12,2)+=0;
    rLHS(12,3)+=0;
    rLHS(12,4)+=crLHS9;
    rLHS(12,5)+=0;
    rLHS(12,6)+=0;
    rLHS(12,7)+=0;
    rLHS(12,8)+=crLHS11;
    rLHS(12,9)+=0;
    rLHS(12,10)+=0;
    rLHS(12,11)+=0;
    rLHS(12,12)+=crLHS12;
    rLHS(12,13)+=0;
    rLHS(12,14)+=0;
    rLHS(12,15)+=0;
    rLHS(13,0)+=0;
    rLHS(13,1)+=crLHS5;
    rLHS(13,2)+=0;
    rLHS(13,3)+=0;
    rLHS(13,4)+=0;
    rLHS(13,5)+=crLHS9;
    rLHS(13,6)+=0;
    rLHS(13,7)+=0;
    rLHS(13,8)+=0;
    rLHS(13,9)+=crLHS11;
    rLHS(13,10)+=0;
    rLHS(13,11)+=0;
    rLHS(13,12)+=0;
    rLHS(13,13)+=crLHS12;
    rLHS(13,14)+=0;
    rLHS(13,15)+=0;
    rLHS(14,0)+=0;
    rLHS(14,1)+=0;
    rLHS(14,2)+=crLHS5;
    rLHS(14,3)+=0;
    rLHS(14,4)+=0;
    rLHS(14,5)+=0;
    rLHS(14,6)+=crLHS9;
    rLHS(14,7)+=0;
    rLHS(14,8)+=0;
    rLHS(14,9)+=0;
    rLHS(14,10)+=crLHS11;
    rLHS(14,11)+=0;
    rLHS(14,12)+=0;
    rLHS(14,13)+=0;
    rLHS(14,14)+=crLHS12;
    rLHS(14,15)+=0;
    rLHS(15,0)+=0;
    rLHS(15,1)+=0;
    rLHS(15,2)+=0;
    rLHS(15,3)+=0;
    rLHS(15,4)+=0;
    rLHS(15,5)+=0;
    rLHS(15,6)+=0;
    rLHS(15,7)+=0;
    rLHS(15,8)+=0;
    rLHS(15,9)+=0;
    rLHS(15,10)+=0;
    rLHS(15,11)+=0;
    rLHS(15,12)+=0;
    rLHS(15,13)+=0;
    rLHS(15,14)+=0;
    rLHS(15,15)+=0;

}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSOutletPenaltyContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    VectorType &rRHS)
{
    double c_penalty = 0.0;

    const double mu = rData.EffectiveViscosity;
    const double max_spectral_radius = rData.MaxSprectraRadius;
    double kappa = c_penalty * mu;
    double gauss_weight = rData.Weight;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;

        const double crRHS0 = 1.0/(max_spectral_radius + 1);
    const double crRHS1 = DN(0,0)*(crRHS0*(v(0,0) - vn(0,0)) + vn(0,0)) + DN(1,0)*(crRHS0*(v(1,0) - vn(1,0)) + vn(1,0)) + DN(2,0)*(crRHS0*(v(2,0) - vn(2,0)) + vn(2,0));
    const double crRHS2 = gauss_weight*kappa;
    const double crRHS3 = DN(0,0)*crRHS2;
    const double crRHS4 = DN(0,0)*(crRHS0*(v(0,1) - vn(0,1)) + vn(0,1)) + DN(1,0)*(crRHS0*(v(1,1) - vn(1,1)) + vn(1,1)) + DN(2,0)*(crRHS0*(v(2,1) - vn(2,1)) + vn(2,1));
    const double crRHS5 = DN(1,0)*crRHS2;
    const double crRHS6 = DN(2,0)*crRHS2;
    rRHS[0]+=crRHS1*crRHS3;
    rRHS[1]+=crRHS3*crRHS4;
    rRHS[2]+=0;
    rRHS[3]+=crRHS1*crRHS5;
    rRHS[4]+=crRHS4*crRHS5;
    rRHS[5]+=0;
    rRHS[6]+=crRHS1*crRHS6;
    rRHS[7]+=crRHS4*crRHS6;
    rRHS[8]+=0;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSOutletPenaltyContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    VectorType &rRHS)
{
    double c_penalty = 10.0;
    const double mu = rData.EffectiveViscosity;
    const double max_spectral_radius = rData.MaxSprectraRadius;
    double kappa = c_penalty * mu;
    double gauss_weight = rData.Weight;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;

        const double crRHS0 = 1.0/(max_spectral_radius + 1);
    const double crRHS1 = DN(0,0)*(crRHS0*(v(0,0) - vn(0,0)) + vn(0,0)) + DN(1,0)*(crRHS0*(v(1,0) - vn(1,0)) + vn(1,0)) + DN(2,0)*(crRHS0*(v(2,0) - vn(2,0)) + vn(2,0)) + DN(3,0)*(crRHS0*(v(3,0) - vn(3,0)) + vn(3,0));
    const double crRHS2 = gauss_weight*kappa;
    const double crRHS3 = DN(0,0)*crRHS2;
    const double crRHS4 = DN(0,0)*(crRHS0*(v(0,1) - vn(0,1)) + vn(0,1)) + DN(1,0)*(crRHS0*(v(1,1) - vn(1,1)) + vn(1,1)) + DN(2,0)*(crRHS0*(v(2,1) - vn(2,1)) + vn(2,1)) + DN(3,0)*(crRHS0*(v(3,1) - vn(3,1)) + vn(3,1));
    const double crRHS5 = DN(0,0)*(crRHS0*(v(0,2) - vn(0,2)) + vn(0,2)) + DN(1,0)*(crRHS0*(v(1,2) - vn(1,2)) + vn(1,2)) + DN(2,0)*(crRHS0*(v(2,2) - vn(2,2)) + vn(2,2)) + DN(3,0)*(crRHS0*(v(3,2) - vn(3,2)) + vn(3,2));
    const double crRHS6 = DN(1,0)*crRHS2;
    const double crRHS7 = DN(2,0)*crRHS2;
    const double crRHS8 = DN(3,0)*crRHS2;
    rRHS[0]+=crRHS1*crRHS3;
    rRHS[1]+=crRHS3*crRHS4;
    rRHS[2]+=crRHS3*crRHS5;
    rRHS[3]+=0;
    rRHS[4]+=crRHS1*crRHS6;
    rRHS[5]+=crRHS4*crRHS6;
    rRHS[6]+=crRHS5*crRHS6;
    rRHS[7]+=0;
    rRHS[8]+=crRHS1*crRHS7;
    rRHS[9]+=crRHS4*crRHS7;
    rRHS[10]+=crRHS5*crRHS7;
    rRHS[11]+=0;
    rRHS[12]+=crRHS1*crRHS8;
    rRHS[13]+=crRHS4*crRHS8;
    rRHS[14]+=crRHS5*crRHS8;
    rRHS[15]+=0;

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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;
    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = 1.0/(max_spectral_radius + 1);
const double clhs2 = DN(0,0)*clhs1;
const double clhs3 = C(0,2)*DN(0,0);
const double clhs4 = C(2,2)*DN(0,1) + clhs3;
const double clhs5 = DN(0,1)*clhs1;
const double clhs6 = pow(DN(0,0), 2);
const double clhs7 = 1.0/dt;
const double clhs8 = clhs7*pow(h, 2)*not_stabilization_cut_elements_mass*rho/(not_stabilization_cut_elements_momentum*stab_c1);
const double clhs9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs11 = DN(0,0)*clhs9 + DN(0,1)*clhs10;
const double clhs12 = clhs1*rho;
const double clhs13 = N[0]*clhs12;
const double clhs14 = clhs1*(0.5*max_spectral_radius - 1.5);
const double clhs15 = 0.5*clhs7*(max_spectral_radius - 3);
const double clhs16 = clhs15/(-clhs1 - clhs14 + 0.5);
const double clhs17 = -N[0]*clhs16 + clhs11;
const double clhs18 = dt*not_stabilization_cut_elements_momentum;
const double clhs19 = clhs11*clhs18;
const double clhs20 = clhs12*clhs19;
const double clhs21 = clhs15/(clhs1 + clhs14 - 0.5);
const double clhs22 = clhs12*clhs21;
const double clhs23 = clhs17*clhs18;
const double clhs24 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs25 = clhs13*clhs24;
const double clhs26 = pow(N[0], 2)*clhs22 + clhs11*clhs13 + clhs17*clhs20 + clhs23*clhs25;
const double clhs27 = C(0,1)*DN(0,1) + clhs3;
const double clhs28 = C(1,2)*DN(0,1);
const double clhs29 = C(2,2)*DN(0,0) + clhs28;
const double clhs30 = DN(0,0)*clhs8;
const double clhs31 = DN(0,1)*clhs30;
const double clhs32 = clhs18*clhs24;
const double clhs33 = N[0]*clhs32 - N[0] + clhs19;
const double clhs34 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs35 = C(0,2)*DN(1,0);
const double clhs36 = C(2,2)*DN(1,1) + clhs35;
const double clhs37 = clhs13*clhs21;
const double clhs38 = N[1]*clhs37;
const double clhs39 = DN(0,0)*DN(1,0);
const double clhs40 = clhs38 + clhs39*clhs8;
const double clhs41 = DN(1,0)*clhs9 + DN(1,1)*clhs10;
const double clhs42 = -N[1]*clhs16 + clhs41;
const double clhs43 = clhs18*clhs42;
const double clhs44 = clhs13*clhs41 + clhs20*clhs42 + clhs25*clhs43;
const double clhs45 = C(0,1)*DN(1,1) + clhs35;
const double clhs46 = C(1,2)*DN(1,1);
const double clhs47 = C(2,2)*DN(1,0) + clhs46;
const double clhs48 = DN(1,1)*clhs30;
const double clhs49 = DN(0,0)*N[1];
const double clhs50 = DN(1,0)*N[0];
const double clhs51 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs52 = C(0,2)*DN(2,0);
const double clhs53 = C(2,2)*DN(2,1) + clhs52;
const double clhs54 = N[2]*clhs37;
const double clhs55 = DN(0,0)*DN(2,0);
const double clhs56 = clhs54 + clhs55*clhs8;
const double clhs57 = DN(2,0)*clhs9 + DN(2,1)*clhs10;
const double clhs58 = -N[2]*clhs16 + clhs57;
const double clhs59 = clhs18*clhs58;
const double clhs60 = clhs13*clhs57 + clhs20*clhs58 + clhs25*clhs59;
const double clhs61 = C(0,1)*DN(2,1) + clhs52;
const double clhs62 = C(1,2)*DN(2,1);
const double clhs63 = C(2,2)*DN(2,0) + clhs62;
const double clhs64 = DN(2,1)*clhs30;
const double clhs65 = DN(0,0)*N[2];
const double clhs66 = DN(2,0)*N[0];
const double clhs67 = C(0,1)*DN(0,0) + clhs28;
const double clhs68 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs69 = pow(DN(0,1), 2);
const double clhs70 = C(0,1)*DN(1,0) + clhs46;
const double clhs71 = DN(0,1)*clhs8;
const double clhs72 = DN(1,0)*clhs71;
const double clhs73 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs74 = DN(0,1)*DN(1,1);
const double clhs75 = clhs38 + clhs74*clhs8;
const double clhs76 = DN(0,1)*N[1];
const double clhs77 = DN(1,1)*N[0];
const double clhs78 = C(0,1)*DN(2,0) + clhs62;
const double clhs79 = DN(2,0)*clhs71;
const double clhs80 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs81 = DN(0,1)*DN(2,1);
const double clhs82 = clhs54 + clhs8*clhs81;
const double clhs83 = DN(0,1)*N[2];
const double clhs84 = DN(2,1)*N[0];
const double clhs85 = N[0] + clhs23;
const double clhs86 = clhs18/rho;
const double clhs87 = clhs86*(clhs39 + clhs74);
const double clhs88 = clhs86*(clhs55 + clhs81);
const double clhs89 = DN(1,0)*clhs1;
const double clhs90 = DN(1,1)*clhs1;
const double clhs91 = N[1]*clhs12;
const double clhs92 = clhs12*clhs23;
const double clhs93 = clhs24*clhs91;
const double clhs94 = clhs11*clhs91 + clhs23*clhs93 + clhs41*clhs92;
const double clhs95 = clhs18*clhs41;
const double clhs96 = pow(DN(1,0), 2);
const double clhs97 = clhs12*clhs95;
const double clhs98 = pow(N[1], 2)*clhs22 + clhs41*clhs91 + clhs42*clhs97 + clhs43*clhs93;
const double clhs99 = DN(1,0)*clhs8;
const double clhs100 = DN(1,1)*clhs99;
const double clhs101 = N[1]*clhs32 - N[1] + clhs95;
const double clhs102 = N[2]*clhs21*clhs91;
const double clhs103 = DN(1,0)*DN(2,0);
const double clhs104 = clhs102 + clhs103*clhs8;
const double clhs105 = clhs57*clhs91 + clhs58*clhs97 + clhs59*clhs93;
const double clhs106 = DN(2,1)*clhs99;
const double clhs107 = DN(1,0)*N[2];
const double clhs108 = DN(2,0)*N[1];
const double clhs109 = pow(DN(1,1), 2);
const double clhs110 = DN(2,0)*clhs8;
const double clhs111 = DN(1,1)*clhs110;
const double clhs112 = DN(1,1)*DN(2,1);
const double clhs113 = clhs102 + clhs112*clhs8;
const double clhs114 = DN(1,1)*N[2];
const double clhs115 = DN(2,1)*N[1];
const double clhs116 = N[1] + clhs43;
const double clhs117 = clhs86*(clhs103 + clhs112);
const double clhs118 = DN(2,0)*clhs1;
const double clhs119 = DN(2,1)*clhs1;
const double clhs120 = N[2]*clhs12;
const double clhs121 = clhs120*clhs24;
const double clhs122 = clhs11*clhs120 + clhs121*clhs23 + clhs57*clhs92;
const double clhs123 = clhs18*clhs57;
const double clhs124 = clhs12*clhs43*clhs57 + clhs120*clhs41 + clhs121*clhs43;
const double clhs125 = pow(DN(2,0), 2);
const double clhs126 = pow(N[2], 2)*clhs22 + clhs12*clhs123*clhs58 + clhs120*clhs57 + clhs121*clhs59;
const double clhs127 = DN(2,1)*clhs110;
const double clhs128 = N[2]*clhs32 - N[2] + clhs123;
const double clhs129 = pow(DN(2,1), 2);
const double clhs130 = N[2] + clhs59;
lhs(0,0)=clhs0*clhs2 + clhs26 + clhs4*clhs5 + clhs6*clhs8;
lhs(0,1)=clhs2*clhs27 + clhs29*clhs5 + clhs31;
lhs(0,2)=DN(0,0)*clhs33;
lhs(0,3)=clhs2*clhs34 + clhs36*clhs5 + clhs40 + clhs44;
lhs(0,4)=clhs2*clhs45 + clhs47*clhs5 + clhs48;
lhs(0,5)=DN(1,0)*clhs19 + clhs32*clhs50 - clhs49;
lhs(0,6)=clhs2*clhs51 + clhs5*clhs53 + clhs56 + clhs60;
lhs(0,7)=clhs2*clhs61 + clhs5*clhs63 + clhs64;
lhs(0,8)=DN(2,0)*clhs19 + clhs32*clhs66 - clhs65;
lhs(1,0)=clhs2*clhs4 + clhs31 + clhs5*clhs67;
lhs(1,1)=clhs2*clhs29 + clhs26 + clhs5*clhs68 + clhs69*clhs8;
lhs(1,2)=DN(0,1)*clhs33;
lhs(1,3)=clhs2*clhs36 + clhs5*clhs70 + clhs72;
lhs(1,4)=clhs2*clhs47 + clhs44 + clhs5*clhs73 + clhs75;
lhs(1,5)=DN(1,1)*clhs19 + clhs32*clhs77 - clhs76;
lhs(1,6)=clhs2*clhs53 + clhs5*clhs78 + clhs79;
lhs(1,7)=clhs2*clhs63 + clhs5*clhs80 + clhs60 + clhs82;
lhs(1,8)=DN(2,1)*clhs19 + clhs32*clhs84 - clhs83;
lhs(2,0)=clhs2*clhs85;
lhs(2,1)=clhs5*clhs85;
lhs(2,2)=clhs86*(clhs6 + clhs69);
lhs(2,3)=clhs1*(DN(0,0)*clhs43 + clhs50);
lhs(2,4)=clhs1*(DN(0,1)*clhs43 + clhs77);
lhs(2,5)=clhs87;
lhs(2,6)=clhs1*(DN(0,0)*clhs59 + clhs66);
lhs(2,7)=clhs1*(DN(0,1)*clhs59 + clhs84);
lhs(2,8)=clhs88;
lhs(3,0)=clhs0*clhs89 + clhs4*clhs90 + clhs40 + clhs94;
lhs(3,1)=clhs27*clhs89 + clhs29*clhs90 + clhs72;
lhs(3,2)=DN(0,0)*clhs95 + clhs32*clhs49 - clhs50;
lhs(3,3)=clhs34*clhs89 + clhs36*clhs90 + clhs8*clhs96 + clhs98;
lhs(3,4)=clhs100 + clhs45*clhs89 + clhs47*clhs90;
lhs(3,5)=DN(1,0)*clhs101;
lhs(3,6)=clhs104 + clhs105 + clhs51*clhs89 + clhs53*clhs90;
lhs(3,7)=clhs106 + clhs61*clhs89 + clhs63*clhs90;
lhs(3,8)=DN(2,0)*clhs95 - clhs107 + clhs108*clhs32;
lhs(4,0)=clhs4*clhs89 + clhs48 + clhs67*clhs90;
lhs(4,1)=clhs29*clhs89 + clhs68*clhs90 + clhs75 + clhs94;
lhs(4,2)=DN(0,1)*clhs95 + clhs32*clhs76 - clhs77;
lhs(4,3)=clhs100 + clhs36*clhs89 + clhs70*clhs90;
lhs(4,4)=clhs109*clhs8 + clhs47*clhs89 + clhs73*clhs90 + clhs98;
lhs(4,5)=DN(1,1)*clhs101;
lhs(4,6)=clhs111 + clhs53*clhs89 + clhs78*clhs90;
lhs(4,7)=clhs105 + clhs113 + clhs63*clhs89 + clhs80*clhs90;
lhs(4,8)=DN(2,1)*clhs95 - clhs114 + clhs115*clhs32;
lhs(5,0)=clhs1*(DN(1,0)*clhs23 + clhs49);
lhs(5,1)=clhs1*(DN(1,1)*clhs23 + clhs76);
lhs(5,2)=clhs87;
lhs(5,3)=clhs116*clhs89;
lhs(5,4)=clhs116*clhs90;
lhs(5,5)=clhs86*(clhs109 + clhs96);
lhs(5,6)=clhs1*(DN(1,0)*clhs59 + clhs108);
lhs(5,7)=clhs1*(DN(1,1)*clhs59 + clhs115);
lhs(5,8)=clhs117;
lhs(6,0)=clhs0*clhs118 + clhs119*clhs4 + clhs122 + clhs56;
lhs(6,1)=clhs118*clhs27 + clhs119*clhs29 + clhs79;
lhs(6,2)=DN(0,0)*clhs123 + clhs32*clhs65 - clhs66;
lhs(6,3)=clhs104 + clhs118*clhs34 + clhs119*clhs36 + clhs124;
lhs(6,4)=clhs111 + clhs118*clhs45 + clhs119*clhs47;
lhs(6,5)=DN(1,0)*clhs123 + clhs107*clhs32 - clhs108;
lhs(6,6)=clhs118*clhs51 + clhs119*clhs53 + clhs125*clhs8 + clhs126;
lhs(6,7)=clhs118*clhs61 + clhs119*clhs63 + clhs127;
lhs(6,8)=DN(2,0)*clhs128;
lhs(7,0)=clhs118*clhs4 + clhs119*clhs67 + clhs64;
lhs(7,1)=clhs118*clhs29 + clhs119*clhs68 + clhs122 + clhs82;
lhs(7,2)=DN(0,1)*clhs123 + clhs32*clhs83 - clhs84;
lhs(7,3)=clhs106 + clhs118*clhs36 + clhs119*clhs70;
lhs(7,4)=clhs113 + clhs118*clhs47 + clhs119*clhs73 + clhs124;
lhs(7,5)=DN(1,1)*clhs123 + clhs114*clhs32 - clhs115;
lhs(7,6)=clhs118*clhs53 + clhs119*clhs78 + clhs127;
lhs(7,7)=clhs118*clhs63 + clhs119*clhs80 + clhs126 + clhs129*clhs8;
lhs(7,8)=DN(2,1)*clhs128;
lhs(8,0)=clhs1*(DN(2,0)*clhs23 + clhs65);
lhs(8,1)=clhs1*(DN(2,1)*clhs23 + clhs83);
lhs(8,2)=clhs88;
lhs(8,3)=clhs1*(DN(2,0)*clhs43 + clhs107);
lhs(8,4)=clhs1*(DN(2,1)*clhs43 + clhs114);
lhs(8,5)=clhs117;
lhs(8,6)=clhs118*clhs130;
lhs(8,7)=clhs119*clhs130;
lhs(8,8)=clhs86*(clhs125 + clhs129);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;

    // Add outlet penalty contribution
    if (this->Is(OUTLET) && !rData.IsCut()) {
        TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSOutletPenaltyContribution(rData, rLHS);
    }
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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;
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

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = 1.0/(max_spectral_radius + 1);
const double clhs2 = DN(0,0)*clhs1;
const double clhs3 = C(0,3)*DN(0,0);
const double clhs4 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs3;
const double clhs5 = DN(0,1)*clhs1;
const double clhs6 = C(0,5)*DN(0,0);
const double clhs7 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs6;
const double clhs8 = DN(0,2)*clhs1;
const double clhs9 = pow(DN(0,0), 2);
const double clhs10 = 1.0/dt;
const double clhs11 = clhs10*pow(h, 2)*not_stabilization_cut_elements_mass*rho/(not_stabilization_cut_elements_momentum*stab_c1);
const double clhs12 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs14 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs15 = DN(0,0)*clhs12 + DN(0,1)*clhs13 + DN(0,2)*clhs14;
const double clhs16 = clhs1*rho;
const double clhs17 = N[0]*clhs16;
const double clhs18 = clhs1*(0.5*max_spectral_radius - 1.5);
const double clhs19 = 0.5*clhs10*(max_spectral_radius - 3);
const double clhs20 = clhs19/(-clhs1 - clhs18 + 0.5);
const double clhs21 = -N[0]*clhs20 + clhs15;
const double clhs22 = dt*not_stabilization_cut_elements_momentum;
const double clhs23 = clhs15*clhs22;
const double clhs24 = clhs16*clhs23;
const double clhs25 = clhs19/(clhs1 + clhs18 - 0.5);
const double clhs26 = clhs16*clhs25;
const double clhs27 = clhs21*clhs22;
const double clhs28 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs29 = clhs17*clhs28;
const double clhs30 = pow(N[0], 2)*clhs26 + clhs15*clhs17 + clhs21*clhs24 + clhs27*clhs29;
const double clhs31 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs3;
const double clhs32 = C(1,3)*DN(0,1);
const double clhs33 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs32;
const double clhs34 = C(3,5)*DN(0,0);
const double clhs35 = C(4,5)*DN(0,2);
const double clhs36 = C(1,5)*DN(0,1) + clhs34 + clhs35;
const double clhs37 = DN(0,0)*clhs11;
const double clhs38 = DN(0,1)*clhs37;
const double clhs39 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs6;
const double clhs40 = C(3,4)*DN(0,1);
const double clhs41 = C(2,3)*DN(0,2) + clhs34 + clhs40;
const double clhs42 = C(2,5)*DN(0,2);
const double clhs43 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs42;
const double clhs44 = DN(0,2)*clhs37;
const double clhs45 = clhs22*clhs28;
const double clhs46 = N[0]*clhs45 - N[0] + clhs23;
const double clhs47 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs48 = C(0,3)*DN(1,0);
const double clhs49 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs48;
const double clhs50 = C(0,5)*DN(1,0);
const double clhs51 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs50;
const double clhs52 = clhs17*clhs25;
const double clhs53 = N[1]*clhs52;
const double clhs54 = DN(0,0)*DN(1,0);
const double clhs55 = clhs11*clhs54 + clhs53;
const double clhs56 = DN(1,0)*clhs12 + DN(1,1)*clhs13 + DN(1,2)*clhs14;
const double clhs57 = -N[1]*clhs20 + clhs56;
const double clhs58 = clhs22*clhs57;
const double clhs59 = clhs17*clhs56 + clhs24*clhs57 + clhs29*clhs58;
const double clhs60 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs48;
const double clhs61 = C(1,3)*DN(1,1);
const double clhs62 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs61;
const double clhs63 = C(3,5)*DN(1,0);
const double clhs64 = C(4,5)*DN(1,2);
const double clhs65 = C(1,5)*DN(1,1) + clhs63 + clhs64;
const double clhs66 = DN(1,1)*clhs37;
const double clhs67 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs50;
const double clhs68 = C(3,4)*DN(1,1);
const double clhs69 = C(2,3)*DN(1,2) + clhs63 + clhs68;
const double clhs70 = C(2,5)*DN(1,2);
const double clhs71 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs70;
const double clhs72 = DN(1,2)*clhs37;
const double clhs73 = DN(0,0)*N[1];
const double clhs74 = DN(1,0)*N[0];
const double clhs75 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs76 = C(0,3)*DN(2,0);
const double clhs77 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs76;
const double clhs78 = C(0,5)*DN(2,0);
const double clhs79 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs78;
const double clhs80 = N[2]*clhs52;
const double clhs81 = DN(0,0)*DN(2,0);
const double clhs82 = clhs11*clhs81 + clhs80;
const double clhs83 = DN(2,0)*clhs12 + DN(2,1)*clhs13 + DN(2,2)*clhs14;
const double clhs84 = -N[2]*clhs20 + clhs83;
const double clhs85 = clhs22*clhs84;
const double clhs86 = clhs17*clhs83 + clhs24*clhs84 + clhs29*clhs85;
const double clhs87 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs76;
const double clhs88 = C(1,3)*DN(2,1);
const double clhs89 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs88;
const double clhs90 = C(3,5)*DN(2,0);
const double clhs91 = C(4,5)*DN(2,2);
const double clhs92 = C(1,5)*DN(2,1) + clhs90 + clhs91;
const double clhs93 = DN(2,1)*clhs37;
const double clhs94 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs78;
const double clhs95 = C(3,4)*DN(2,1);
const double clhs96 = C(2,3)*DN(2,2) + clhs90 + clhs95;
const double clhs97 = C(2,5)*DN(2,2);
const double clhs98 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs97;
const double clhs99 = DN(2,2)*clhs37;
const double clhs100 = DN(0,0)*N[2];
const double clhs101 = DN(2,0)*N[0];
const double clhs102 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs103 = C(0,3)*DN(3,0);
const double clhs104 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs103;
const double clhs105 = C(0,5)*DN(3,0);
const double clhs106 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs105;
const double clhs107 = N[3]*clhs52;
const double clhs108 = DN(0,0)*DN(3,0);
const double clhs109 = clhs107 + clhs108*clhs11;
const double clhs110 = DN(3,0)*clhs12 + DN(3,1)*clhs13 + DN(3,2)*clhs14;
const double clhs111 = -N[3]*clhs20 + clhs110;
const double clhs112 = clhs111*clhs22;
const double clhs113 = clhs110*clhs17 + clhs111*clhs24 + clhs112*clhs29;
const double clhs114 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs103;
const double clhs115 = C(1,3)*DN(3,1);
const double clhs116 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs115;
const double clhs117 = C(3,5)*DN(3,0);
const double clhs118 = C(4,5)*DN(3,2);
const double clhs119 = C(1,5)*DN(3,1) + clhs117 + clhs118;
const double clhs120 = DN(3,1)*clhs37;
const double clhs121 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs105;
const double clhs122 = C(3,4)*DN(3,1);
const double clhs123 = C(2,3)*DN(3,2) + clhs117 + clhs122;
const double clhs124 = C(2,5)*DN(3,2);
const double clhs125 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs124;
const double clhs126 = DN(3,2)*clhs37;
const double clhs127 = DN(0,0)*N[3];
const double clhs128 = DN(3,0)*N[0];
const double clhs129 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs32;
const double clhs130 = C(0,4)*DN(0,0) + clhs35 + clhs40;
const double clhs131 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs132 = C(1,4)*DN(0,1);
const double clhs133 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs132;
const double clhs134 = pow(DN(0,1), 2);
const double clhs135 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs132;
const double clhs136 = C(2,4)*DN(0,2);
const double clhs137 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs136;
const double clhs138 = DN(0,1)*clhs11;
const double clhs139 = DN(0,2)*clhs138;
const double clhs140 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs61;
const double clhs141 = C(0,4)*DN(1,0) + clhs64 + clhs68;
const double clhs142 = DN(1,0)*clhs138;
const double clhs143 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs144 = C(1,4)*DN(1,1);
const double clhs145 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs144;
const double clhs146 = DN(0,1)*DN(1,1);
const double clhs147 = clhs11*clhs146;
const double clhs148 = clhs53 + clhs59;
const double clhs149 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs144;
const double clhs150 = C(2,4)*DN(1,2);
const double clhs151 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs150;
const double clhs152 = DN(1,2)*clhs138;
const double clhs153 = DN(0,1)*N[1];
const double clhs154 = DN(1,1)*N[0];
const double clhs155 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs88;
const double clhs156 = C(0,4)*DN(2,0) + clhs91 + clhs95;
const double clhs157 = DN(2,0)*clhs138;
const double clhs158 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs159 = C(1,4)*DN(2,1);
const double clhs160 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs159;
const double clhs161 = DN(0,1)*DN(2,1);
const double clhs162 = clhs11*clhs161;
const double clhs163 = clhs80 + clhs86;
const double clhs164 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs159;
const double clhs165 = C(2,4)*DN(2,2);
const double clhs166 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs165;
const double clhs167 = DN(2,2)*clhs138;
const double clhs168 = DN(0,1)*N[2];
const double clhs169 = DN(2,1)*N[0];
const double clhs170 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs115;
const double clhs171 = C(0,4)*DN(3,0) + clhs118 + clhs122;
const double clhs172 = DN(3,0)*clhs138;
const double clhs173 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs174 = C(1,4)*DN(3,1);
const double clhs175 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs174;
const double clhs176 = DN(0,1)*DN(3,1);
const double clhs177 = clhs11*clhs176;
const double clhs178 = clhs107 + clhs113;
const double clhs179 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs174;
const double clhs180 = C(2,4)*DN(3,2);
const double clhs181 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs180;
const double clhs182 = DN(3,2)*clhs138;
const double clhs183 = DN(0,1)*N[3];
const double clhs184 = DN(3,1)*N[0];
const double clhs185 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs42;
const double clhs186 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs136;
const double clhs187 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs188 = pow(DN(0,2), 2);
const double clhs189 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs70;
const double clhs190 = DN(0,2)*clhs11;
const double clhs191 = DN(1,0)*clhs190;
const double clhs192 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs150;
const double clhs193 = DN(1,1)*clhs190;
const double clhs194 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs195 = DN(0,2)*DN(1,2);
const double clhs196 = clhs11*clhs195;
const double clhs197 = DN(0,2)*N[1];
const double clhs198 = DN(1,2)*N[0];
const double clhs199 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs97;
const double clhs200 = DN(2,0)*clhs190;
const double clhs201 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs165;
const double clhs202 = DN(2,1)*clhs190;
const double clhs203 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs204 = DN(0,2)*DN(2,2);
const double clhs205 = clhs11*clhs204;
const double clhs206 = DN(0,2)*N[2];
const double clhs207 = DN(2,2)*N[0];
const double clhs208 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs124;
const double clhs209 = DN(3,0)*clhs190;
const double clhs210 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs180;
const double clhs211 = DN(3,1)*clhs190;
const double clhs212 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs213 = DN(0,2)*DN(3,2);
const double clhs214 = clhs11*clhs213;
const double clhs215 = DN(0,2)*N[3];
const double clhs216 = DN(3,2)*N[0];
const double clhs217 = N[0] + clhs27;
const double clhs218 = clhs22/rho;
const double clhs219 = clhs218*(clhs146 + clhs195 + clhs54);
const double clhs220 = clhs218*(clhs161 + clhs204 + clhs81);
const double clhs221 = clhs218*(clhs108 + clhs176 + clhs213);
const double clhs222 = DN(1,0)*clhs1;
const double clhs223 = DN(1,1)*clhs1;
const double clhs224 = DN(1,2)*clhs1;
const double clhs225 = N[1]*clhs16;
const double clhs226 = clhs16*clhs27;
const double clhs227 = clhs225*clhs28;
const double clhs228 = clhs15*clhs225 + clhs226*clhs56 + clhs227*clhs27;
const double clhs229 = clhs22*clhs56;
const double clhs230 = pow(DN(1,0), 2);
const double clhs231 = clhs16*clhs229;
const double clhs232 = pow(N[1], 2)*clhs26 + clhs225*clhs56 + clhs227*clhs58 + clhs231*clhs57;
const double clhs233 = DN(1,0)*clhs11;
const double clhs234 = DN(1,1)*clhs233;
const double clhs235 = DN(1,2)*clhs233;
const double clhs236 = N[1]*clhs45 - N[1] + clhs229;
const double clhs237 = clhs225*clhs25;
const double clhs238 = N[2]*clhs237;
const double clhs239 = DN(1,0)*DN(2,0);
const double clhs240 = clhs11*clhs239 + clhs238;
const double clhs241 = clhs225*clhs83 + clhs227*clhs85 + clhs231*clhs84;
const double clhs242 = DN(2,1)*clhs233;
const double clhs243 = DN(2,2)*clhs233;
const double clhs244 = DN(1,0)*N[2];
const double clhs245 = DN(2,0)*N[1];
const double clhs246 = N[3]*clhs237;
const double clhs247 = DN(1,0)*DN(3,0);
const double clhs248 = clhs11*clhs247 + clhs246;
const double clhs249 = clhs110*clhs225 + clhs111*clhs231 + clhs112*clhs227;
const double clhs250 = DN(3,1)*clhs233;
const double clhs251 = DN(3,2)*clhs233;
const double clhs252 = DN(1,0)*N[3];
const double clhs253 = DN(3,0)*N[1];
const double clhs254 = clhs228 + clhs53;
const double clhs255 = pow(DN(1,1), 2);
const double clhs256 = DN(1,1)*clhs11;
const double clhs257 = DN(1,2)*clhs256;
const double clhs258 = DN(2,0)*clhs256;
const double clhs259 = DN(1,1)*DN(2,1);
const double clhs260 = clhs11*clhs259;
const double clhs261 = clhs238 + clhs241;
const double clhs262 = DN(2,2)*clhs256;
const double clhs263 = DN(1,1)*N[2];
const double clhs264 = DN(2,1)*N[1];
const double clhs265 = DN(3,0)*clhs256;
const double clhs266 = DN(1,1)*DN(3,1);
const double clhs267 = clhs11*clhs266;
const double clhs268 = clhs246 + clhs249;
const double clhs269 = DN(3,2)*clhs256;
const double clhs270 = DN(1,1)*N[3];
const double clhs271 = DN(3,1)*N[1];
const double clhs272 = pow(DN(1,2), 2);
const double clhs273 = DN(1,2)*clhs11;
const double clhs274 = DN(2,0)*clhs273;
const double clhs275 = DN(2,1)*clhs273;
const double clhs276 = DN(1,2)*DN(2,2);
const double clhs277 = clhs11*clhs276;
const double clhs278 = DN(1,2)*N[2];
const double clhs279 = DN(2,2)*N[1];
const double clhs280 = DN(3,0)*clhs273;
const double clhs281 = DN(3,1)*clhs273;
const double clhs282 = DN(1,2)*DN(3,2);
const double clhs283 = clhs11*clhs282;
const double clhs284 = DN(1,2)*N[3];
const double clhs285 = DN(3,2)*N[1];
const double clhs286 = N[1] + clhs58;
const double clhs287 = clhs218*(clhs239 + clhs259 + clhs276);
const double clhs288 = clhs218*(clhs247 + clhs266 + clhs282);
const double clhs289 = DN(2,0)*clhs1;
const double clhs290 = DN(2,1)*clhs1;
const double clhs291 = DN(2,2)*clhs1;
const double clhs292 = N[2]*clhs16;
const double clhs293 = clhs28*clhs292;
const double clhs294 = clhs15*clhs292 + clhs226*clhs83 + clhs27*clhs293;
const double clhs295 = clhs22*clhs83;
const double clhs296 = clhs16*clhs58;
const double clhs297 = clhs292*clhs56 + clhs293*clhs58 + clhs296*clhs83;
const double clhs298 = pow(DN(2,0), 2);
const double clhs299 = clhs16*clhs295;
const double clhs300 = pow(N[2], 2)*clhs26 + clhs292*clhs83 + clhs293*clhs85 + clhs299*clhs84;
const double clhs301 = DN(2,0)*clhs11;
const double clhs302 = DN(2,1)*clhs301;
const double clhs303 = DN(2,2)*clhs301;
const double clhs304 = N[2]*clhs45 - N[2] + clhs295;
const double clhs305 = N[3]*clhs25*clhs292;
const double clhs306 = DN(2,0)*DN(3,0);
const double clhs307 = clhs11*clhs306 + clhs305;
const double clhs308 = clhs110*clhs292 + clhs111*clhs299 + clhs112*clhs293;
const double clhs309 = DN(3,1)*clhs301;
const double clhs310 = DN(3,2)*clhs301;
const double clhs311 = DN(2,0)*N[3];
const double clhs312 = DN(3,0)*N[2];
const double clhs313 = clhs294 + clhs80;
const double clhs314 = clhs238 + clhs297;
const double clhs315 = pow(DN(2,1), 2);
const double clhs316 = DN(2,1)*clhs11;
const double clhs317 = DN(2,2)*clhs316;
const double clhs318 = DN(3,0)*clhs316;
const double clhs319 = DN(2,1)*DN(3,1);
const double clhs320 = clhs11*clhs319;
const double clhs321 = clhs305 + clhs308;
const double clhs322 = DN(3,2)*clhs316;
const double clhs323 = DN(2,1)*N[3];
const double clhs324 = DN(3,1)*N[2];
const double clhs325 = pow(DN(2,2), 2);
const double clhs326 = DN(2,2)*clhs11;
const double clhs327 = DN(3,0)*clhs326;
const double clhs328 = DN(3,1)*clhs326;
const double clhs329 = DN(2,2)*DN(3,2);
const double clhs330 = clhs11*clhs329;
const double clhs331 = DN(2,2)*N[3];
const double clhs332 = DN(3,2)*N[2];
const double clhs333 = N[2] + clhs85;
const double clhs334 = clhs218*(clhs306 + clhs319 + clhs329);
const double clhs335 = DN(3,0)*clhs1;
const double clhs336 = DN(3,1)*clhs1;
const double clhs337 = DN(3,2)*clhs1;
const double clhs338 = N[3]*clhs16;
const double clhs339 = clhs28*clhs338;
const double clhs340 = clhs110*clhs226 + clhs15*clhs338 + clhs27*clhs339;
const double clhs341 = clhs110*clhs22;
const double clhs342 = clhs110*clhs296 + clhs338*clhs56 + clhs339*clhs58;
const double clhs343 = clhs110*clhs16*clhs85 + clhs338*clhs83 + clhs339*clhs85;
const double clhs344 = pow(DN(3,0), 2);
const double clhs345 = pow(N[3], 2)*clhs26 + clhs110*clhs338 + clhs111*clhs16*clhs341 + clhs112*clhs339;
const double clhs346 = DN(3,0)*clhs11;
const double clhs347 = DN(3,1)*clhs346;
const double clhs348 = DN(3,2)*clhs346;
const double clhs349 = N[3]*clhs45 - N[3] + clhs341;
const double clhs350 = clhs107 + clhs340;
const double clhs351 = clhs246 + clhs342;
const double clhs352 = clhs305 + clhs343;
const double clhs353 = pow(DN(3,1), 2);
const double clhs354 = DN(3,1)*DN(3,2)*clhs11;
const double clhs355 = pow(DN(3,2), 2);
const double clhs356 = N[3] + clhs112;
lhs(0,0)=clhs0*clhs2 + clhs11*clhs9 + clhs30 + clhs4*clhs5 + clhs7*clhs8;
lhs(0,1)=clhs2*clhs31 + clhs33*clhs5 + clhs36*clhs8 + clhs38;
lhs(0,2)=clhs2*clhs39 + clhs41*clhs5 + clhs43*clhs8 + clhs44;
lhs(0,3)=DN(0,0)*clhs46;
lhs(0,4)=clhs2*clhs47 + clhs49*clhs5 + clhs51*clhs8 + clhs55 + clhs59;
lhs(0,5)=clhs2*clhs60 + clhs5*clhs62 + clhs65*clhs8 + clhs66;
lhs(0,6)=clhs2*clhs67 + clhs5*clhs69 + clhs71*clhs8 + clhs72;
lhs(0,7)=DN(1,0)*clhs23 + clhs45*clhs74 - clhs73;
lhs(0,8)=clhs2*clhs75 + clhs5*clhs77 + clhs79*clhs8 + clhs82 + clhs86;
lhs(0,9)=clhs2*clhs87 + clhs5*clhs89 + clhs8*clhs92 + clhs93;
lhs(0,10)=clhs2*clhs94 + clhs5*clhs96 + clhs8*clhs98 + clhs99;
lhs(0,11)=DN(2,0)*clhs23 - clhs100 + clhs101*clhs45;
lhs(0,12)=clhs102*clhs2 + clhs104*clhs5 + clhs106*clhs8 + clhs109 + clhs113;
lhs(0,13)=clhs114*clhs2 + clhs116*clhs5 + clhs119*clhs8 + clhs120;
lhs(0,14)=clhs121*clhs2 + clhs123*clhs5 + clhs125*clhs8 + clhs126;
lhs(0,15)=DN(3,0)*clhs23 - clhs127 + clhs128*clhs45;
lhs(1,0)=clhs129*clhs5 + clhs130*clhs8 + clhs2*clhs4 + clhs38;
lhs(1,1)=clhs11*clhs134 + clhs131*clhs5 + clhs133*clhs8 + clhs2*clhs33 + clhs30;
lhs(1,2)=clhs135*clhs5 + clhs137*clhs8 + clhs139 + clhs2*clhs41;
lhs(1,3)=DN(0,1)*clhs46;
lhs(1,4)=clhs140*clhs5 + clhs141*clhs8 + clhs142 + clhs2*clhs49;
lhs(1,5)=clhs143*clhs5 + clhs145*clhs8 + clhs147 + clhs148 + clhs2*clhs62;
lhs(1,6)=clhs149*clhs5 + clhs151*clhs8 + clhs152 + clhs2*clhs69;
lhs(1,7)=DN(1,1)*clhs23 - clhs153 + clhs154*clhs45;
lhs(1,8)=clhs155*clhs5 + clhs156*clhs8 + clhs157 + clhs2*clhs77;
lhs(1,9)=clhs158*clhs5 + clhs160*clhs8 + clhs162 + clhs163 + clhs2*clhs89;
lhs(1,10)=clhs164*clhs5 + clhs166*clhs8 + clhs167 + clhs2*clhs96;
lhs(1,11)=DN(2,1)*clhs23 - clhs168 + clhs169*clhs45;
lhs(1,12)=clhs104*clhs2 + clhs170*clhs5 + clhs171*clhs8 + clhs172;
lhs(1,13)=clhs116*clhs2 + clhs173*clhs5 + clhs175*clhs8 + clhs177 + clhs178;
lhs(1,14)=clhs123*clhs2 + clhs179*clhs5 + clhs181*clhs8 + clhs182;
lhs(1,15)=DN(3,1)*clhs23 - clhs183 + clhs184*clhs45;
lhs(2,0)=clhs130*clhs5 + clhs185*clhs8 + clhs2*clhs7 + clhs44;
lhs(2,1)=clhs133*clhs5 + clhs139 + clhs186*clhs8 + clhs2*clhs36;
lhs(2,2)=clhs11*clhs188 + clhs137*clhs5 + clhs187*clhs8 + clhs2*clhs43 + clhs30;
lhs(2,3)=DN(0,2)*clhs46;
lhs(2,4)=clhs141*clhs5 + clhs189*clhs8 + clhs191 + clhs2*clhs51;
lhs(2,5)=clhs145*clhs5 + clhs192*clhs8 + clhs193 + clhs2*clhs65;
lhs(2,6)=clhs148 + clhs151*clhs5 + clhs194*clhs8 + clhs196 + clhs2*clhs71;
lhs(2,7)=DN(1,2)*clhs23 - clhs197 + clhs198*clhs45;
lhs(2,8)=clhs156*clhs5 + clhs199*clhs8 + clhs2*clhs79 + clhs200;
lhs(2,9)=clhs160*clhs5 + clhs2*clhs92 + clhs201*clhs8 + clhs202;
lhs(2,10)=clhs163 + clhs166*clhs5 + clhs2*clhs98 + clhs203*clhs8 + clhs205;
lhs(2,11)=DN(2,2)*clhs23 - clhs206 + clhs207*clhs45;
lhs(2,12)=clhs106*clhs2 + clhs171*clhs5 + clhs208*clhs8 + clhs209;
lhs(2,13)=clhs119*clhs2 + clhs175*clhs5 + clhs210*clhs8 + clhs211;
lhs(2,14)=clhs125*clhs2 + clhs178 + clhs181*clhs5 + clhs212*clhs8 + clhs214;
lhs(2,15)=DN(3,2)*clhs23 - clhs215 + clhs216*clhs45;
lhs(3,0)=clhs2*clhs217;
lhs(3,1)=clhs217*clhs5;
lhs(3,2)=clhs217*clhs8;
lhs(3,3)=clhs218*(clhs134 + clhs188 + clhs9);
lhs(3,4)=clhs1*(DN(0,0)*clhs58 + clhs74);
lhs(3,5)=clhs1*(DN(0,1)*clhs58 + clhs154);
lhs(3,6)=clhs1*(DN(0,2)*clhs58 + clhs198);
lhs(3,7)=clhs219;
lhs(3,8)=clhs1*(DN(0,0)*clhs85 + clhs101);
lhs(3,9)=clhs1*(DN(0,1)*clhs85 + clhs169);
lhs(3,10)=clhs1*(DN(0,2)*clhs85 + clhs207);
lhs(3,11)=clhs220;
lhs(3,12)=clhs1*(DN(0,0)*clhs112 + clhs128);
lhs(3,13)=clhs1*(DN(0,1)*clhs112 + clhs184);
lhs(3,14)=clhs1*(DN(0,2)*clhs112 + clhs216);
lhs(3,15)=clhs221;
lhs(4,0)=clhs0*clhs222 + clhs223*clhs4 + clhs224*clhs7 + clhs228 + clhs55;
lhs(4,1)=clhs142 + clhs222*clhs31 + clhs223*clhs33 + clhs224*clhs36;
lhs(4,2)=clhs191 + clhs222*clhs39 + clhs223*clhs41 + clhs224*clhs43;
lhs(4,3)=DN(0,0)*clhs229 + clhs45*clhs73 - clhs74;
lhs(4,4)=clhs11*clhs230 + clhs222*clhs47 + clhs223*clhs49 + clhs224*clhs51 + clhs232;
lhs(4,5)=clhs222*clhs60 + clhs223*clhs62 + clhs224*clhs65 + clhs234;
lhs(4,6)=clhs222*clhs67 + clhs223*clhs69 + clhs224*clhs71 + clhs235;
lhs(4,7)=DN(1,0)*clhs236;
lhs(4,8)=clhs222*clhs75 + clhs223*clhs77 + clhs224*clhs79 + clhs240 + clhs241;
lhs(4,9)=clhs222*clhs87 + clhs223*clhs89 + clhs224*clhs92 + clhs242;
lhs(4,10)=clhs222*clhs94 + clhs223*clhs96 + clhs224*clhs98 + clhs243;
lhs(4,11)=DN(2,0)*clhs229 - clhs244 + clhs245*clhs45;
lhs(4,12)=clhs102*clhs222 + clhs104*clhs223 + clhs106*clhs224 + clhs248 + clhs249;
lhs(4,13)=clhs114*clhs222 + clhs116*clhs223 + clhs119*clhs224 + clhs250;
lhs(4,14)=clhs121*clhs222 + clhs123*clhs223 + clhs125*clhs224 + clhs251;
lhs(4,15)=DN(3,0)*clhs229 - clhs252 + clhs253*clhs45;
lhs(5,0)=clhs129*clhs223 + clhs130*clhs224 + clhs222*clhs4 + clhs66;
lhs(5,1)=clhs131*clhs223 + clhs133*clhs224 + clhs147 + clhs222*clhs33 + clhs254;
lhs(5,2)=clhs135*clhs223 + clhs137*clhs224 + clhs193 + clhs222*clhs41;
lhs(5,3)=DN(0,1)*clhs229 + clhs153*clhs45 - clhs154;
lhs(5,4)=clhs140*clhs223 + clhs141*clhs224 + clhs222*clhs49 + clhs234;
lhs(5,5)=clhs11*clhs255 + clhs143*clhs223 + clhs145*clhs224 + clhs222*clhs62 + clhs232;
lhs(5,6)=clhs149*clhs223 + clhs151*clhs224 + clhs222*clhs69 + clhs257;
lhs(5,7)=DN(1,1)*clhs236;
lhs(5,8)=clhs155*clhs223 + clhs156*clhs224 + clhs222*clhs77 + clhs258;
lhs(5,9)=clhs158*clhs223 + clhs160*clhs224 + clhs222*clhs89 + clhs260 + clhs261;
lhs(5,10)=clhs164*clhs223 + clhs166*clhs224 + clhs222*clhs96 + clhs262;
lhs(5,11)=DN(2,1)*clhs229 - clhs263 + clhs264*clhs45;
lhs(5,12)=clhs104*clhs222 + clhs170*clhs223 + clhs171*clhs224 + clhs265;
lhs(5,13)=clhs116*clhs222 + clhs173*clhs223 + clhs175*clhs224 + clhs267 + clhs268;
lhs(5,14)=clhs123*clhs222 + clhs179*clhs223 + clhs181*clhs224 + clhs269;
lhs(5,15)=DN(3,1)*clhs229 - clhs270 + clhs271*clhs45;
lhs(6,0)=clhs130*clhs223 + clhs185*clhs224 + clhs222*clhs7 + clhs72;
lhs(6,1)=clhs133*clhs223 + clhs152 + clhs186*clhs224 + clhs222*clhs36;
lhs(6,2)=clhs137*clhs223 + clhs187*clhs224 + clhs196 + clhs222*clhs43 + clhs254;
lhs(6,3)=DN(0,2)*clhs229 + clhs197*clhs45 - clhs198;
lhs(6,4)=clhs141*clhs223 + clhs189*clhs224 + clhs222*clhs51 + clhs235;
lhs(6,5)=clhs145*clhs223 + clhs192*clhs224 + clhs222*clhs65 + clhs257;
lhs(6,6)=clhs11*clhs272 + clhs151*clhs223 + clhs194*clhs224 + clhs222*clhs71 + clhs232;
lhs(6,7)=DN(1,2)*clhs236;
lhs(6,8)=clhs156*clhs223 + clhs199*clhs224 + clhs222*clhs79 + clhs274;
lhs(6,9)=clhs160*clhs223 + clhs201*clhs224 + clhs222*clhs92 + clhs275;
lhs(6,10)=clhs166*clhs223 + clhs203*clhs224 + clhs222*clhs98 + clhs261 + clhs277;
lhs(6,11)=DN(2,2)*clhs229 - clhs278 + clhs279*clhs45;
lhs(6,12)=clhs106*clhs222 + clhs171*clhs223 + clhs208*clhs224 + clhs280;
lhs(6,13)=clhs119*clhs222 + clhs175*clhs223 + clhs210*clhs224 + clhs281;
lhs(6,14)=clhs125*clhs222 + clhs181*clhs223 + clhs212*clhs224 + clhs268 + clhs283;
lhs(6,15)=DN(3,2)*clhs229 - clhs284 + clhs285*clhs45;
lhs(7,0)=clhs1*(DN(1,0)*clhs27 + clhs73);
lhs(7,1)=clhs1*(DN(1,1)*clhs27 + clhs153);
lhs(7,2)=clhs1*(DN(1,2)*clhs27 + clhs197);
lhs(7,3)=clhs219;
lhs(7,4)=clhs222*clhs286;
lhs(7,5)=clhs223*clhs286;
lhs(7,6)=clhs224*clhs286;
lhs(7,7)=clhs218*(clhs230 + clhs255 + clhs272);
lhs(7,8)=clhs1*(DN(1,0)*clhs85 + clhs245);
lhs(7,9)=clhs1*(DN(1,1)*clhs85 + clhs264);
lhs(7,10)=clhs1*(DN(1,2)*clhs85 + clhs279);
lhs(7,11)=clhs287;
lhs(7,12)=clhs1*(DN(1,0)*clhs112 + clhs253);
lhs(7,13)=clhs1*(DN(1,1)*clhs112 + clhs271);
lhs(7,14)=clhs1*(DN(1,2)*clhs112 + clhs285);
lhs(7,15)=clhs288;
lhs(8,0)=clhs0*clhs289 + clhs290*clhs4 + clhs291*clhs7 + clhs294 + clhs82;
lhs(8,1)=clhs157 + clhs289*clhs31 + clhs290*clhs33 + clhs291*clhs36;
lhs(8,2)=clhs200 + clhs289*clhs39 + clhs290*clhs41 + clhs291*clhs43;
lhs(8,3)=DN(0,0)*clhs295 + clhs100*clhs45 - clhs101;
lhs(8,4)=clhs240 + clhs289*clhs47 + clhs290*clhs49 + clhs291*clhs51 + clhs297;
lhs(8,5)=clhs258 + clhs289*clhs60 + clhs290*clhs62 + clhs291*clhs65;
lhs(8,6)=clhs274 + clhs289*clhs67 + clhs290*clhs69 + clhs291*clhs71;
lhs(8,7)=DN(1,0)*clhs295 + clhs244*clhs45 - clhs245;
lhs(8,8)=clhs11*clhs298 + clhs289*clhs75 + clhs290*clhs77 + clhs291*clhs79 + clhs300;
lhs(8,9)=clhs289*clhs87 + clhs290*clhs89 + clhs291*clhs92 + clhs302;
lhs(8,10)=clhs289*clhs94 + clhs290*clhs96 + clhs291*clhs98 + clhs303;
lhs(8,11)=DN(2,0)*clhs304;
lhs(8,12)=clhs102*clhs289 + clhs104*clhs290 + clhs106*clhs291 + clhs307 + clhs308;
lhs(8,13)=clhs114*clhs289 + clhs116*clhs290 + clhs119*clhs291 + clhs309;
lhs(8,14)=clhs121*clhs289 + clhs123*clhs290 + clhs125*clhs291 + clhs310;
lhs(8,15)=DN(3,0)*clhs295 - clhs311 + clhs312*clhs45;
lhs(9,0)=clhs129*clhs290 + clhs130*clhs291 + clhs289*clhs4 + clhs93;
lhs(9,1)=clhs131*clhs290 + clhs133*clhs291 + clhs162 + clhs289*clhs33 + clhs313;
lhs(9,2)=clhs135*clhs290 + clhs137*clhs291 + clhs202 + clhs289*clhs41;
lhs(9,3)=DN(0,1)*clhs295 + clhs168*clhs45 - clhs169;
lhs(9,4)=clhs140*clhs290 + clhs141*clhs291 + clhs242 + clhs289*clhs49;
lhs(9,5)=clhs143*clhs290 + clhs145*clhs291 + clhs260 + clhs289*clhs62 + clhs314;
lhs(9,6)=clhs149*clhs290 + clhs151*clhs291 + clhs275 + clhs289*clhs69;
lhs(9,7)=DN(1,1)*clhs295 + clhs263*clhs45 - clhs264;
lhs(9,8)=clhs155*clhs290 + clhs156*clhs291 + clhs289*clhs77 + clhs302;
lhs(9,9)=clhs11*clhs315 + clhs158*clhs290 + clhs160*clhs291 + clhs289*clhs89 + clhs300;
lhs(9,10)=clhs164*clhs290 + clhs166*clhs291 + clhs289*clhs96 + clhs317;
lhs(9,11)=DN(2,1)*clhs304;
lhs(9,12)=clhs104*clhs289 + clhs170*clhs290 + clhs171*clhs291 + clhs318;
lhs(9,13)=clhs116*clhs289 + clhs173*clhs290 + clhs175*clhs291 + clhs320 + clhs321;
lhs(9,14)=clhs123*clhs289 + clhs179*clhs290 + clhs181*clhs291 + clhs322;
lhs(9,15)=DN(3,1)*clhs295 - clhs323 + clhs324*clhs45;
lhs(10,0)=clhs130*clhs290 + clhs185*clhs291 + clhs289*clhs7 + clhs99;
lhs(10,1)=clhs133*clhs290 + clhs167 + clhs186*clhs291 + clhs289*clhs36;
lhs(10,2)=clhs137*clhs290 + clhs187*clhs291 + clhs205 + clhs289*clhs43 + clhs313;
lhs(10,3)=DN(0,2)*clhs295 + clhs206*clhs45 - clhs207;
lhs(10,4)=clhs141*clhs290 + clhs189*clhs291 + clhs243 + clhs289*clhs51;
lhs(10,5)=clhs145*clhs290 + clhs192*clhs291 + clhs262 + clhs289*clhs65;
lhs(10,6)=clhs151*clhs290 + clhs194*clhs291 + clhs277 + clhs289*clhs71 + clhs314;
lhs(10,7)=DN(1,2)*clhs295 + clhs278*clhs45 - clhs279;
lhs(10,8)=clhs156*clhs290 + clhs199*clhs291 + clhs289*clhs79 + clhs303;
lhs(10,9)=clhs160*clhs290 + clhs201*clhs291 + clhs289*clhs92 + clhs317;
lhs(10,10)=clhs11*clhs325 + clhs166*clhs290 + clhs203*clhs291 + clhs289*clhs98 + clhs300;
lhs(10,11)=DN(2,2)*clhs304;
lhs(10,12)=clhs106*clhs289 + clhs171*clhs290 + clhs208*clhs291 + clhs327;
lhs(10,13)=clhs119*clhs289 + clhs175*clhs290 + clhs210*clhs291 + clhs328;
lhs(10,14)=clhs125*clhs289 + clhs181*clhs290 + clhs212*clhs291 + clhs321 + clhs330;
lhs(10,15)=DN(3,2)*clhs295 - clhs331 + clhs332*clhs45;
lhs(11,0)=clhs1*(DN(2,0)*clhs27 + clhs100);
lhs(11,1)=clhs1*(DN(2,1)*clhs27 + clhs168);
lhs(11,2)=clhs1*(DN(2,2)*clhs27 + clhs206);
lhs(11,3)=clhs220;
lhs(11,4)=clhs1*(DN(2,0)*clhs58 + clhs244);
lhs(11,5)=clhs1*(DN(2,1)*clhs58 + clhs263);
lhs(11,6)=clhs1*(DN(2,2)*clhs58 + clhs278);
lhs(11,7)=clhs287;
lhs(11,8)=clhs289*clhs333;
lhs(11,9)=clhs290*clhs333;
lhs(11,10)=clhs291*clhs333;
lhs(11,11)=clhs218*(clhs298 + clhs315 + clhs325);
lhs(11,12)=clhs1*(DN(2,0)*clhs112 + clhs312);
lhs(11,13)=clhs1*(DN(2,1)*clhs112 + clhs324);
lhs(11,14)=clhs1*(DN(2,2)*clhs112 + clhs332);
lhs(11,15)=clhs334;
lhs(12,0)=clhs0*clhs335 + clhs109 + clhs336*clhs4 + clhs337*clhs7 + clhs340;
lhs(12,1)=clhs172 + clhs31*clhs335 + clhs33*clhs336 + clhs337*clhs36;
lhs(12,2)=clhs209 + clhs335*clhs39 + clhs336*clhs41 + clhs337*clhs43;
lhs(12,3)=DN(0,0)*clhs341 + clhs127*clhs45 - clhs128;
lhs(12,4)=clhs248 + clhs335*clhs47 + clhs336*clhs49 + clhs337*clhs51 + clhs342;
lhs(12,5)=clhs265 + clhs335*clhs60 + clhs336*clhs62 + clhs337*clhs65;
lhs(12,6)=clhs280 + clhs335*clhs67 + clhs336*clhs69 + clhs337*clhs71;
lhs(12,7)=DN(1,0)*clhs341 + clhs252*clhs45 - clhs253;
lhs(12,8)=clhs307 + clhs335*clhs75 + clhs336*clhs77 + clhs337*clhs79 + clhs343;
lhs(12,9)=clhs318 + clhs335*clhs87 + clhs336*clhs89 + clhs337*clhs92;
lhs(12,10)=clhs327 + clhs335*clhs94 + clhs336*clhs96 + clhs337*clhs98;
lhs(12,11)=DN(2,0)*clhs341 + clhs311*clhs45 - clhs312;
lhs(12,12)=clhs102*clhs335 + clhs104*clhs336 + clhs106*clhs337 + clhs11*clhs344 + clhs345;
lhs(12,13)=clhs114*clhs335 + clhs116*clhs336 + clhs119*clhs337 + clhs347;
lhs(12,14)=clhs121*clhs335 + clhs123*clhs336 + clhs125*clhs337 + clhs348;
lhs(12,15)=DN(3,0)*clhs349;
lhs(13,0)=clhs120 + clhs129*clhs336 + clhs130*clhs337 + clhs335*clhs4;
lhs(13,1)=clhs131*clhs336 + clhs133*clhs337 + clhs177 + clhs33*clhs335 + clhs350;
lhs(13,2)=clhs135*clhs336 + clhs137*clhs337 + clhs211 + clhs335*clhs41;
lhs(13,3)=DN(0,1)*clhs341 + clhs183*clhs45 - clhs184;
lhs(13,4)=clhs140*clhs336 + clhs141*clhs337 + clhs250 + clhs335*clhs49;
lhs(13,5)=clhs143*clhs336 + clhs145*clhs337 + clhs267 + clhs335*clhs62 + clhs351;
lhs(13,6)=clhs149*clhs336 + clhs151*clhs337 + clhs281 + clhs335*clhs69;
lhs(13,7)=DN(1,1)*clhs341 + clhs270*clhs45 - clhs271;
lhs(13,8)=clhs155*clhs336 + clhs156*clhs337 + clhs309 + clhs335*clhs77;
lhs(13,9)=clhs158*clhs336 + clhs160*clhs337 + clhs320 + clhs335*clhs89 + clhs352;
lhs(13,10)=clhs164*clhs336 + clhs166*clhs337 + clhs328 + clhs335*clhs96;
lhs(13,11)=DN(2,1)*clhs341 + clhs323*clhs45 - clhs324;
lhs(13,12)=clhs104*clhs335 + clhs170*clhs336 + clhs171*clhs337 + clhs347;
lhs(13,13)=clhs11*clhs353 + clhs116*clhs335 + clhs173*clhs336 + clhs175*clhs337 + clhs345;
lhs(13,14)=clhs123*clhs335 + clhs179*clhs336 + clhs181*clhs337 + clhs354;
lhs(13,15)=DN(3,1)*clhs349;
lhs(14,0)=clhs126 + clhs130*clhs336 + clhs185*clhs337 + clhs335*clhs7;
lhs(14,1)=clhs133*clhs336 + clhs182 + clhs186*clhs337 + clhs335*clhs36;
lhs(14,2)=clhs137*clhs336 + clhs187*clhs337 + clhs214 + clhs335*clhs43 + clhs350;
lhs(14,3)=DN(0,2)*clhs341 + clhs215*clhs45 - clhs216;
lhs(14,4)=clhs141*clhs336 + clhs189*clhs337 + clhs251 + clhs335*clhs51;
lhs(14,5)=clhs145*clhs336 + clhs192*clhs337 + clhs269 + clhs335*clhs65;
lhs(14,6)=clhs151*clhs336 + clhs194*clhs337 + clhs283 + clhs335*clhs71 + clhs351;
lhs(14,7)=DN(1,2)*clhs341 + clhs284*clhs45 - clhs285;
lhs(14,8)=clhs156*clhs336 + clhs199*clhs337 + clhs310 + clhs335*clhs79;
lhs(14,9)=clhs160*clhs336 + clhs201*clhs337 + clhs322 + clhs335*clhs92;
lhs(14,10)=clhs166*clhs336 + clhs203*clhs337 + clhs330 + clhs335*clhs98 + clhs352;
lhs(14,11)=DN(2,2)*clhs341 + clhs331*clhs45 - clhs332;
lhs(14,12)=clhs106*clhs335 + clhs171*clhs336 + clhs208*clhs337 + clhs348;
lhs(14,13)=clhs119*clhs335 + clhs175*clhs336 + clhs210*clhs337 + clhs354;
lhs(14,14)=clhs11*clhs355 + clhs125*clhs335 + clhs181*clhs336 + clhs212*clhs337 + clhs345;
lhs(14,15)=DN(3,2)*clhs349;
lhs(15,0)=clhs1*(DN(3,0)*clhs27 + clhs127);
lhs(15,1)=clhs1*(DN(3,1)*clhs27 + clhs183);
lhs(15,2)=clhs1*(DN(3,2)*clhs27 + clhs215);
lhs(15,3)=clhs221;
lhs(15,4)=clhs1*(DN(3,0)*clhs58 + clhs252);
lhs(15,5)=clhs1*(DN(3,1)*clhs58 + clhs270);
lhs(15,6)=clhs1*(DN(3,2)*clhs58 + clhs284);
lhs(15,7)=clhs288;
lhs(15,8)=clhs1*(DN(3,0)*clhs85 + clhs311);
lhs(15,9)=clhs1*(DN(3,1)*clhs85 + clhs323);
lhs(15,10)=clhs1*(DN(3,2)*clhs85 + clhs331);
lhs(15,11)=clhs334;
lhs(15,12)=clhs335*clhs356;
lhs(15,13)=clhs336*clhs356;
lhs(15,14)=clhs337*clhs356;
lhs(15,15)=clhs218*(clhs344 + clhs353 + clhs355);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;

    // Add outlet penalty contribution
    if (this->Is(OUTLET)) {
        TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSOutletPenaltyContribution(rData, rLHS);
    }
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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;
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
const double crhs1 = 1.0/dt;
const double crhs2 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs3 = crhs1*pow(h, 2)*not_stabilization_cut_elements_mass*rho*(crhs2 - volume_error_ratio)/(not_stabilization_cut_elements_momentum*stab_c1);
const double crhs4 = 1.0/(max_spectral_radius + 1);
const double crhs5 = rho*(N[0]*(crhs4*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs4*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs4*(f(2,0) - fn(2,0)) + fn(2,0)));
const double crhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs7 = v(0,0) - vn(0,0);
const double crhs8 = crhs4*crhs7 + vn(0,0);
const double crhs9 = v(1,0) - vn(1,0);
const double crhs10 = crhs4*crhs9 + vn(1,0);
const double crhs11 = v(2,0) - vn(2,0);
const double crhs12 = crhs11*crhs4 + vn(2,0);
const double crhs13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs14 = rho*(crhs13*(DN(0,1)*crhs8 + DN(1,1)*crhs10 + DN(2,1)*crhs12) + crhs6*(DN(0,0)*crhs8 + DN(1,0)*crhs10 + DN(2,0)*crhs12));
const double crhs15 = -acceleration_alpha_method(0,0);
const double crhs16 = 0.5*max_spectral_radius;
const double crhs17 = -crhs4;
const double crhs18 = crhs17 + 0.5;
const double crhs19 = 1.0/(crhs18 - crhs4*(crhs16 - 1.5));
const double crhs20 = crhs1*crhs19;
const double crhs21 = crhs4*(1.5 - crhs16);
const double crhs22 = crhs19*(-crhs21 + crhs4 + 0.5);
const double crhs23 = 0.5*crhs4;
const double crhs24 = crhs23*(max_spectral_radius - 3);
const double crhs25 = -acceleration_alpha_method(1,0);
const double crhs26 = -acceleration_alpha_method(2,0);
const double crhs27 = N[0]*(acceleration_alpha_method(0,0) - crhs24*(-acceleration_alpha_method(0,0)*crhs22 + crhs15 + crhs20*crhs7)) + N[1]*(acceleration_alpha_method(1,0) - crhs24*(-acceleration_alpha_method(1,0)*crhs22 + crhs20*crhs9 + crhs25)) + N[2]*(acceleration_alpha_method(2,0) - crhs24*(-acceleration_alpha_method(2,0)*crhs22 + crhs11*crhs20 + crhs26));
const double crhs28 = N[0]*rho;
const double crhs29 = 1.0/(crhs18 + crhs21);
const double crhs30 = crhs1*crhs29;
const double crhs31 = crhs29*(crhs17 + crhs21 - 0.5);
const double crhs32 = crhs23*(3 - max_spectral_radius);
const double crhs33 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs14 - crhs5 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs32*(acceleration_alpha_method(0,0)*crhs31 + crhs15 + crhs30*crhs7)) + N[1]*(acceleration_alpha_method(1,0) + crhs32*(acceleration_alpha_method(1,0)*crhs31 + crhs25 + crhs30*crhs9)) + N[2]*(acceleration_alpha_method(2,0) + crhs32*(acceleration_alpha_method(2,0)*crhs31 + crhs11*crhs30 + crhs26)));
const double crhs34 = dt*not_stabilization_cut_elements_momentum;
const double crhs35 = crhs33*crhs34;
const double crhs36 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs37 = N[0]*crhs36;
const double crhs38 = crhs34*(DN(0,0)*crhs6 + DN(0,1)*crhs13);
const double crhs39 = rho*(N[0]*(crhs4*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs4*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs4*(f(2,1) - fn(2,1)) + fn(2,1)));
const double crhs40 = v(0,1) - vn(0,1);
const double crhs41 = crhs4*crhs40 + vn(0,1);
const double crhs42 = v(1,1) - vn(1,1);
const double crhs43 = crhs4*crhs42 + vn(1,1);
const double crhs44 = v(2,1) - vn(2,1);
const double crhs45 = crhs4*crhs44 + vn(2,1);
const double crhs46 = rho*(crhs13*(DN(0,1)*crhs41 + DN(1,1)*crhs43 + DN(2,1)*crhs45) + crhs6*(DN(0,0)*crhs41 + DN(1,0)*crhs43 + DN(2,0)*crhs45));
const double crhs47 = -acceleration_alpha_method(0,1);
const double crhs48 = -acceleration_alpha_method(1,1);
const double crhs49 = -acceleration_alpha_method(2,1);
const double crhs50 = N[0]*(acceleration_alpha_method(0,1) - crhs24*(-acceleration_alpha_method(0,1)*crhs22 + crhs20*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) - crhs24*(-acceleration_alpha_method(1,1)*crhs22 + crhs20*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) - crhs24*(-acceleration_alpha_method(2,1)*crhs22 + crhs20*crhs44 + crhs49));
const double crhs51 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs39 + crhs46 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs32*(acceleration_alpha_method(0,1)*crhs31 + crhs30*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) + crhs32*(acceleration_alpha_method(1,1)*crhs31 + crhs30*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) + crhs32*(acceleration_alpha_method(2,1)*crhs31 + crhs30*crhs44 + crhs49)));
const double crhs52 = crhs34*crhs51;
const double crhs53 = -crhs4*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + crhs2) + volume_error_ratio;
const double crhs54 = 1.0/rho;
const double crhs55 = crhs35*crhs54;
const double crhs56 = crhs52*crhs54;
const double crhs57 = N[1]*rho;
const double crhs58 = N[1]*crhs36;
const double crhs59 = crhs34*(DN(1,0)*crhs6 + DN(1,1)*crhs13);
const double crhs60 = N[2]*rho;
const double crhs61 = N[2]*crhs36;
const double crhs62 = crhs34*(DN(2,0)*crhs6 + DN(2,1)*crhs13);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs14 + N[0]*crhs5 - crhs27*crhs28 - crhs33*crhs38 - crhs35*crhs37;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] + N[0]*crhs39 - N[0]*crhs46 - crhs28*crhs50 - crhs37*crhs52 - crhs38*crhs51;
rhs[2]=-DN(0,0)*crhs55 - DN(0,1)*crhs56 + N[0]*crhs53;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs14 + N[1]*crhs5 - crhs27*crhs57 - crhs33*crhs59 - crhs35*crhs58;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] + N[1]*crhs39 - N[1]*crhs46 - crhs50*crhs57 - crhs51*crhs59 - crhs52*crhs58;
rhs[5]=-DN(1,0)*crhs55 - DN(1,1)*crhs56 + N[1]*crhs53;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs14 + N[2]*crhs5 - crhs27*crhs60 - crhs33*crhs62 - crhs35*crhs61;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] + N[2]*crhs39 - N[2]*crhs46 - crhs50*crhs60 - crhs51*crhs62 - crhs52*crhs61;
rhs[8]=-DN(2,0)*crhs55 - DN(2,1)*crhs56 + N[2]*crhs53;


    noalias(rRHS) += rData.Weight * rhs;

    // Add outlet penalty contribution
    if (this->Is(OUTLET) && !rData.IsCut())
    {
        TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSOutletPenaltyContribution(rData, rRHS);
    }
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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;
    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = 1.0/dt;
const double crhs2 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs3 = crhs1*pow(h, 2)*not_stabilization_cut_elements_mass*rho*(crhs2 - volume_error_ratio)/(not_stabilization_cut_elements_momentum*stab_c1);
const double crhs4 = 1.0/(max_spectral_radius + 1);
const double crhs5 = rho*(N[0]*(crhs4*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs4*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs4*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs4*(f(3,0) - fn(3,0)) + fn(3,0)));
const double crhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs7 = v(0,0) - vn(0,0);
const double crhs8 = crhs4*crhs7 + vn(0,0);
const double crhs9 = v(1,0) - vn(1,0);
const double crhs10 = crhs4*crhs9 + vn(1,0);
const double crhs11 = v(2,0) - vn(2,0);
const double crhs12 = crhs11*crhs4 + vn(2,0);
const double crhs13 = v(3,0) - vn(3,0);
const double crhs14 = crhs13*crhs4 + vn(3,0);
const double crhs15 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs16 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs17 = rho*(crhs15*(DN(0,1)*crhs8 + DN(1,1)*crhs10 + DN(2,1)*crhs12 + DN(3,1)*crhs14) + crhs16*(DN(0,2)*crhs8 + DN(1,2)*crhs10 + DN(2,2)*crhs12 + DN(3,2)*crhs14) + crhs6*(DN(0,0)*crhs8 + DN(1,0)*crhs10 + DN(2,0)*crhs12 + DN(3,0)*crhs14));
const double crhs18 = -acceleration_alpha_method(0,0);
const double crhs19 = 0.5*max_spectral_radius;
const double crhs20 = -crhs4;
const double crhs21 = crhs20 + 0.5;
const double crhs22 = 1.0/(crhs21 - crhs4*(crhs19 - 1.5));
const double crhs23 = crhs1*crhs22;
const double crhs24 = crhs4*(1.5 - crhs19);
const double crhs25 = crhs22*(-crhs24 + crhs4 + 0.5);
const double crhs26 = 0.5*crhs4;
const double crhs27 = crhs26*(max_spectral_radius - 3);
const double crhs28 = -acceleration_alpha_method(1,0);
const double crhs29 = -acceleration_alpha_method(2,0);
const double crhs30 = -acceleration_alpha_method(3,0);
const double crhs31 = N[0]*(acceleration_alpha_method(0,0) - crhs27*(-acceleration_alpha_method(0,0)*crhs25 + crhs18 + crhs23*crhs7)) + N[1]*(acceleration_alpha_method(1,0) - crhs27*(-acceleration_alpha_method(1,0)*crhs25 + crhs23*crhs9 + crhs28)) + N[2]*(acceleration_alpha_method(2,0) - crhs27*(-acceleration_alpha_method(2,0)*crhs25 + crhs11*crhs23 + crhs29)) + N[3]*(acceleration_alpha_method(3,0) - crhs27*(-acceleration_alpha_method(3,0)*crhs25 + crhs13*crhs23 + crhs30));
const double crhs32 = N[0]*rho;
const double crhs33 = 1.0/(crhs21 + crhs24);
const double crhs34 = crhs1*crhs33;
const double crhs35 = crhs33*(crhs20 + crhs24 - 0.5);
const double crhs36 = crhs26*(3 - max_spectral_radius);
const double crhs37 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs17 - crhs5 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs36*(acceleration_alpha_method(0,0)*crhs35 + crhs18 + crhs34*crhs7)) + N[1]*(acceleration_alpha_method(1,0) + crhs36*(acceleration_alpha_method(1,0)*crhs35 + crhs28 + crhs34*crhs9)) + N[2]*(acceleration_alpha_method(2,0) + crhs36*(acceleration_alpha_method(2,0)*crhs35 + crhs11*crhs34 + crhs29)) + N[3]*(acceleration_alpha_method(3,0) + crhs36*(acceleration_alpha_method(3,0)*crhs35 + crhs13*crhs34 + crhs30)));
const double crhs38 = dt*not_stabilization_cut_elements_momentum;
const double crhs39 = crhs37*crhs38;
const double crhs40 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs41 = N[0]*crhs40;
const double crhs42 = crhs38*(DN(0,0)*crhs6 + DN(0,1)*crhs15 + DN(0,2)*crhs16);
const double crhs43 = rho*(N[0]*(crhs4*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs4*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs4*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs4*(f(3,1) - fn(3,1)) + fn(3,1)));
const double crhs44 = v(0,1) - vn(0,1);
const double crhs45 = crhs4*crhs44 + vn(0,1);
const double crhs46 = v(1,1) - vn(1,1);
const double crhs47 = crhs4*crhs46 + vn(1,1);
const double crhs48 = v(2,1) - vn(2,1);
const double crhs49 = crhs4*crhs48 + vn(2,1);
const double crhs50 = v(3,1) - vn(3,1);
const double crhs51 = crhs4*crhs50 + vn(3,1);
const double crhs52 = rho*(crhs15*(DN(0,1)*crhs45 + DN(1,1)*crhs47 + DN(2,1)*crhs49 + DN(3,1)*crhs51) + crhs16*(DN(0,2)*crhs45 + DN(1,2)*crhs47 + DN(2,2)*crhs49 + DN(3,2)*crhs51) + crhs6*(DN(0,0)*crhs45 + DN(1,0)*crhs47 + DN(2,0)*crhs49 + DN(3,0)*crhs51));
const double crhs53 = -acceleration_alpha_method(0,1);
const double crhs54 = -acceleration_alpha_method(1,1);
const double crhs55 = -acceleration_alpha_method(2,1);
const double crhs56 = -acceleration_alpha_method(3,1);
const double crhs57 = N[0]*(acceleration_alpha_method(0,1) - crhs27*(-acceleration_alpha_method(0,1)*crhs25 + crhs23*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) - crhs27*(-acceleration_alpha_method(1,1)*crhs25 + crhs23*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) - crhs27*(-acceleration_alpha_method(2,1)*crhs25 + crhs23*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) - crhs27*(-acceleration_alpha_method(3,1)*crhs25 + crhs23*crhs50 + crhs56));
const double crhs58 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs43 + crhs52 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs36*(acceleration_alpha_method(0,1)*crhs35 + crhs34*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) + crhs36*(acceleration_alpha_method(1,1)*crhs35 + crhs34*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) + crhs36*(acceleration_alpha_method(2,1)*crhs35 + crhs34*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) + crhs36*(acceleration_alpha_method(3,1)*crhs35 + crhs34*crhs50 + crhs56)));
const double crhs59 = crhs38*crhs41;
const double crhs60 = rho*(N[0]*(crhs4*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs4*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs4*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs4*(f(3,2) - fn(3,2)) + fn(3,2)));
const double crhs61 = v(0,2) - vn(0,2);
const double crhs62 = crhs4*crhs61 + vn(0,2);
const double crhs63 = v(1,2) - vn(1,2);
const double crhs64 = crhs4*crhs63 + vn(1,2);
const double crhs65 = v(2,2) - vn(2,2);
const double crhs66 = crhs4*crhs65 + vn(2,2);
const double crhs67 = v(3,2) - vn(3,2);
const double crhs68 = crhs4*crhs67 + vn(3,2);
const double crhs69 = rho*(crhs15*(DN(0,1)*crhs62 + DN(1,1)*crhs64 + DN(2,1)*crhs66 + DN(3,1)*crhs68) + crhs16*(DN(0,2)*crhs62 + DN(1,2)*crhs64 + DN(2,2)*crhs66 + DN(3,2)*crhs68) + crhs6*(DN(0,0)*crhs62 + DN(1,0)*crhs64 + DN(2,0)*crhs66 + DN(3,0)*crhs68));
const double crhs70 = -acceleration_alpha_method(0,2);
const double crhs71 = -acceleration_alpha_method(1,2);
const double crhs72 = -acceleration_alpha_method(2,2);
const double crhs73 = -acceleration_alpha_method(3,2);
const double crhs74 = N[0]*(acceleration_alpha_method(0,2) - crhs27*(-acceleration_alpha_method(0,2)*crhs25 + crhs23*crhs61 + crhs70)) + N[1]*(acceleration_alpha_method(1,2) - crhs27*(-acceleration_alpha_method(1,2)*crhs25 + crhs23*crhs63 + crhs71)) + N[2]*(acceleration_alpha_method(2,2) - crhs27*(-acceleration_alpha_method(2,2)*crhs25 + crhs23*crhs65 + crhs72)) + N[3]*(acceleration_alpha_method(3,2) - crhs27*(-acceleration_alpha_method(3,2)*crhs25 + crhs23*crhs67 + crhs73));
const double crhs75 = DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs60 + crhs69 + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs36*(acceleration_alpha_method(0,2)*crhs35 + crhs34*crhs61 + crhs70)) + N[1]*(acceleration_alpha_method(1,2) + crhs36*(acceleration_alpha_method(1,2)*crhs35 + crhs34*crhs63 + crhs71)) + N[2]*(acceleration_alpha_method(2,2) + crhs36*(acceleration_alpha_method(2,2)*crhs35 + crhs34*crhs65 + crhs72)) + N[3]*(acceleration_alpha_method(3,2) + crhs36*(acceleration_alpha_method(3,2)*crhs35 + crhs34*crhs67 + crhs73)));
const double crhs76 = -crhs4*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + crhs2) + volume_error_ratio;
const double crhs77 = 1.0/rho;
const double crhs78 = crhs39*crhs77;
const double crhs79 = crhs38*crhs77;
const double crhs80 = crhs58*crhs79;
const double crhs81 = crhs75*crhs79;
const double crhs82 = N[1]*rho;
const double crhs83 = N[1]*crhs40;
const double crhs84 = crhs38*(DN(1,0)*crhs6 + DN(1,1)*crhs15 + DN(1,2)*crhs16);
const double crhs85 = crhs38*crhs83;
const double crhs86 = N[2]*rho;
const double crhs87 = N[2]*crhs40;
const double crhs88 = crhs38*(DN(2,0)*crhs6 + DN(2,1)*crhs15 + DN(2,2)*crhs16);
const double crhs89 = crhs38*crhs87;
const double crhs90 = N[3]*rho;
const double crhs91 = N[3]*crhs40;
const double crhs92 = crhs38*(DN(3,0)*crhs6 + DN(3,1)*crhs15 + DN(3,2)*crhs16);
const double crhs93 = crhs38*crhs91;
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs17 + N[0]*crhs5 - crhs31*crhs32 - crhs37*crhs42 - crhs39*crhs41;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs43 - N[0]*crhs52 - crhs32*crhs57 - crhs42*crhs58 - crhs58*crhs59;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs3 - DN(0,2)*stress[2] + N[0]*crhs60 - N[0]*crhs69 - crhs32*crhs74 - crhs42*crhs75 - crhs59*crhs75;
rhs[3]=-DN(0,0)*crhs78 - DN(0,1)*crhs80 - DN(0,2)*crhs81 + N[0]*crhs76;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs17 + N[1]*crhs5 - crhs31*crhs82 - crhs37*crhs84 - crhs39*crhs83;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs43 - N[1]*crhs52 - crhs57*crhs82 - crhs58*crhs84 - crhs58*crhs85;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs3 - DN(1,2)*stress[2] + N[1]*crhs60 - N[1]*crhs69 - crhs74*crhs82 - crhs75*crhs84 - crhs75*crhs85;
rhs[7]=-DN(1,0)*crhs78 - DN(1,1)*crhs80 - DN(1,2)*crhs81 + N[1]*crhs76;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs17 + N[2]*crhs5 - crhs31*crhs86 - crhs37*crhs88 - crhs39*crhs87;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs43 - N[2]*crhs52 - crhs57*crhs86 - crhs58*crhs88 - crhs58*crhs89;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs3 - DN(2,2)*stress[2] + N[2]*crhs60 - N[2]*crhs69 - crhs74*crhs86 - crhs75*crhs88 - crhs75*crhs89;
rhs[11]=-DN(2,0)*crhs78 - DN(2,1)*crhs80 - DN(2,2)*crhs81 + N[2]*crhs76;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs3 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs17 + N[3]*crhs5 - crhs31*crhs90 - crhs37*crhs92 - crhs39*crhs91;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs3 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs43 - N[3]*crhs52 - crhs57*crhs90 - crhs58*crhs92 - crhs58*crhs93;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs3 - DN(3,2)*stress[2] + N[3]*crhs60 - N[3]*crhs69 - crhs74*crhs90 - crhs75*crhs92 - crhs75*crhs93;
rhs[15]=-DN(3,0)*crhs78 - DN(3,1)*crhs80 - DN(3,2)*crhs81 + N[3]*crhs76;


    noalias(rRHS) += rData.Weight * rhs;

    // Add outlet penalty contribution
    if (this->Is(OUTLET) && !rData.IsCut())
    {
        TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSOutletPenaltyContribution(rData, rRHS);
    }
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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;
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

    const double cV0 = dt*not_stabilization_cut_elements_momentum;
const double cV1 = cV0*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV2 = N[0]*cV1;
const double cV3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV5 = cV0*(DN(0,0)*cV3 + DN(0,1)*cV4);
const double cV6 = cV0/rho;
const double cV7 = N[1]*cV1;
const double cV8 = cV0*(DN(1,0)*cV3 + DN(1,1)*cV4);
const double cV9 = N[2]*cV1;
const double cV10 = cV0*(DN(2,0)*cV3 + DN(2,1)*cV4);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV2 + DNenr(0,0)*cV5;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV2 + DNenr(1,0)*cV5;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV2 + DNenr(2,0)*cV5;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV2 + DNenr(0,1)*cV5;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV2 + DNenr(1,1)*cV5;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV2 + DNenr(2,1)*cV5;
V(2,0)=cV6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV7 + DNenr(0,0)*cV8;
V(3,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV7 + DNenr(1,0)*cV8;
V(3,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV7 + DNenr(2,0)*cV8;
V(4,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV7 + DNenr(0,1)*cV8;
V(4,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV7 + DNenr(1,1)*cV8;
V(4,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV7 + DNenr(2,1)*cV8;
V(5,0)=cV6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV10 + DNenr(0,0)*cV9;
V(6,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV9;
V(6,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV9;
V(7,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV9;
V(7,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV9;
V(7,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV9;
V(8,0)=cV6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
V(8,1)=cV6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
V(8,2)=cV6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 = 1.0/(max_spectral_radius + 1);
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH3 = 0.5*(max_spectral_radius - 3)/(dt*(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5));
const double cH4 = dt*not_stabilization_cut_elements_momentum;
const double cH5 = cH4*(DN(0,0)*cH1 + DN(0,1)*cH2 + N[0]*cH3);
const double cH6 = cH4/rho;
const double cH7 = cH4*(DN(1,0)*cH1 + DN(1,1)*cH2 + N[1]*cH3);
const double cH8 = cH4*(DN(2,0)*cH1 + DN(2,1)*cH2 + N[2]*cH3);
H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH5);
H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH5);
H(0,2)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
H(0,3)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH7);
H(0,4)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH7);
H(0,5)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
H(0,6)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH8);
H(0,7)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH8);
H(0,8)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH5);
H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH5);
H(1,2)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
H(1,3)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH7);
H(1,4)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH7);
H(1,5)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
H(1,6)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH8);
H(1,7)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH8);
H(1,8)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH5);
H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH5);
H(2,2)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
H(2,3)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH7);
H(2,4)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH7);
H(2,5)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
H(2,6)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH8);
H(2,7)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH8);
H(2,8)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 = dt*not_stabilization_cut_elements_momentum/rho;
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
const double crhs_ee2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee3 = v(0,0) - vn(0,0);
const double crhs_ee4 = crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 = v(1,0) - vn(1,0);
const double crhs_ee6 = crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 = v(2,0) - vn(2,0);
const double crhs_ee8 = crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee10 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee11 = 1.0/(crhs_ee10 + 0.5);
const double crhs_ee12 = crhs_ee11/dt;
const double crhs_ee13 = crhs_ee11*(crhs_ee10 - 0.5);
const double crhs_ee14 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee15 = dt*not_stabilization_cut_elements_momentum/rho;
const double crhs_ee16 = crhs_ee15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee14*(acceleration_alpha_method(0,0)*crhs_ee13 - acceleration_alpha_method(0,0) + crhs_ee12*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee14*(acceleration_alpha_method(1,0)*crhs_ee13 - acceleration_alpha_method(1,0) + crhs_ee12*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee14*(acceleration_alpha_method(2,0)*crhs_ee13 - acceleration_alpha_method(2,0) + crhs_ee12*crhs_ee7)) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8) + crhs_ee9*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8)));
const double crhs_ee17 = v(0,1) - vn(0,1);
const double crhs_ee18 = crhs_ee0*crhs_ee17 + vn(0,1);
const double crhs_ee19 = v(1,1) - vn(1,1);
const double crhs_ee20 = crhs_ee0*crhs_ee19 + vn(1,1);
const double crhs_ee21 = v(2,1) - vn(2,1);
const double crhs_ee22 = crhs_ee0*crhs_ee21 + vn(2,1);
const double crhs_ee23 = crhs_ee15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee14*(acceleration_alpha_method(0,1)*crhs_ee13 - acceleration_alpha_method(0,1) + crhs_ee12*crhs_ee17)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee14*(acceleration_alpha_method(1,1)*crhs_ee13 - acceleration_alpha_method(1,1) + crhs_ee12*crhs_ee19)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee14*(acceleration_alpha_method(2,1)*crhs_ee13 - acceleration_alpha_method(2,1) + crhs_ee12*crhs_ee21)) + crhs_ee2*(DN(0,0)*crhs_ee18 + DN(1,0)*crhs_ee20 + DN(2,0)*crhs_ee22) + crhs_ee9*(DN(0,1)*crhs_ee18 + DN(1,1)*crhs_ee20 + DN(2,1)*crhs_ee22)));
rhs_ee[0]=-DNenr(0,0)*crhs_ee16 - DNenr(0,1)*crhs_ee23 + Nenr[0]*crhs_ee1;
rhs_ee[1]=-DNenr(1,0)*crhs_ee16 - DNenr(1,1)*crhs_ee23 + Nenr[1]*crhs_ee1;
rhs_ee[2]=-DNenr(2,0)*crhs_ee16 - DNenr(2,1)*crhs_ee23 + Nenr[2]*crhs_ee1;


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
    const double not_stabilization_cut_elements_momentum=rData.NotStabilizationCutElementsMomentum;
    const double not_stabilization_cut_elements_mass=rData.NotStabilizationCutElementsMass;


    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 = dt*not_stabilization_cut_elements_momentum;
const double cV1 = cV0*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV2 = N[0]*cV1;
const double cV3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV5 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV6 = cV0*(DN(0,0)*cV3 + DN(0,1)*cV4 + DN(0,2)*cV5);
const double cV7 = cV0/rho;
const double cV8 = N[1]*cV1;
const double cV9 = cV0*(DN(1,0)*cV3 + DN(1,1)*cV4 + DN(1,2)*cV5);
const double cV10 = N[2]*cV1;
const double cV11 = cV0*(DN(2,0)*cV3 + DN(2,1)*cV4 + DN(2,2)*cV5);
const double cV12 = N[3]*cV1;
const double cV13 = cV0*(DN(3,0)*cV3 + DN(3,1)*cV4 + DN(3,2)*cV5);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV2 + DNenr(0,0)*cV6;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV2 + DNenr(1,0)*cV6;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV2 + DNenr(2,0)*cV6;
V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV2 + DNenr(3,0)*cV6;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV2 + DNenr(0,1)*cV6;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV2 + DNenr(1,1)*cV6;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV2 + DNenr(2,1)*cV6;
V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV2 + DNenr(3,1)*cV6;
V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV2 + DNenr(0,2)*cV6;
V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV2 + DNenr(1,2)*cV6;
V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV2 + DNenr(2,2)*cV6;
V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV2 + DNenr(3,2)*cV6;
V(3,0)=cV7*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV7*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV7*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV7*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
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
V(7,0)=cV7*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV7*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV7*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV7*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
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
V(11,0)=cV7*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV7*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV7*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV7*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
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
V(15,0)=cV7*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
V(15,1)=cV7*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
V(15,2)=cV7*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
V(15,3)=cV7*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 = 1.0/(max_spectral_radius + 1);
const double cH1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH3 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH4 = 0.5*(max_spectral_radius - 3)/(dt*(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5));
const double cH5 = dt*not_stabilization_cut_elements_momentum;
const double cH6 = cH5*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 + N[0]*cH4);
const double cH7 = cH5/rho;
const double cH8 = cH5*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 + N[1]*cH4);
const double cH9 = cH5*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 + N[2]*cH4);
const double cH10 = cH5*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 + N[3]*cH4);
H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH6);
H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH6);
H(0,2)=cH0*(DN(0,2)*Nenr[0] + DNenr(0,2)*cH6);
H(0,3)=cH7*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
H(0,4)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH8);
H(0,5)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH8);
H(0,6)=cH0*(DN(1,2)*Nenr[0] + DNenr(0,2)*cH8);
H(0,7)=cH7*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
H(0,8)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH9);
H(0,9)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH9);
H(0,10)=cH0*(DN(2,2)*Nenr[0] + DNenr(0,2)*cH9);
H(0,11)=cH7*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
H(0,12)=cH0*(DN(3,0)*Nenr[0] + DNenr(0,0)*cH10);
H(0,13)=cH0*(DN(3,1)*Nenr[0] + DNenr(0,1)*cH10);
H(0,14)=cH0*(DN(3,2)*Nenr[0] + DNenr(0,2)*cH10);
H(0,15)=cH7*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH6);
H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH6);
H(1,2)=cH0*(DN(0,2)*Nenr[1] + DNenr(1,2)*cH6);
H(1,3)=cH7*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
H(1,4)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH8);
H(1,5)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH8);
H(1,6)=cH0*(DN(1,2)*Nenr[1] + DNenr(1,2)*cH8);
H(1,7)=cH7*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
H(1,8)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH9);
H(1,9)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH9);
H(1,10)=cH0*(DN(2,2)*Nenr[1] + DNenr(1,2)*cH9);
H(1,11)=cH7*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
H(1,12)=cH0*(DN(3,0)*Nenr[1] + DNenr(1,0)*cH10);
H(1,13)=cH0*(DN(3,1)*Nenr[1] + DNenr(1,1)*cH10);
H(1,14)=cH0*(DN(3,2)*Nenr[1] + DNenr(1,2)*cH10);
H(1,15)=cH7*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH6);
H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH6);
H(2,2)=cH0*(DN(0,2)*Nenr[2] + DNenr(2,2)*cH6);
H(2,3)=cH7*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
H(2,4)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH8);
H(2,5)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH8);
H(2,6)=cH0*(DN(1,2)*Nenr[2] + DNenr(2,2)*cH8);
H(2,7)=cH7*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
H(2,8)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH9);
H(2,9)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH9);
H(2,10)=cH0*(DN(2,2)*Nenr[2] + DNenr(2,2)*cH9);
H(2,11)=cH7*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
H(2,12)=cH0*(DN(3,0)*Nenr[2] + DNenr(2,0)*cH10);
H(2,13)=cH0*(DN(3,1)*Nenr[2] + DNenr(2,1)*cH10);
H(2,14)=cH0*(DN(3,2)*Nenr[2] + DNenr(2,2)*cH10);
H(2,15)=cH7*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
H(3,0)=cH0*(DN(0,0)*Nenr[3] + DNenr(3,0)*cH6);
H(3,1)=cH0*(DN(0,1)*Nenr[3] + DNenr(3,1)*cH6);
H(3,2)=cH0*(DN(0,2)*Nenr[3] + DNenr(3,2)*cH6);
H(3,3)=cH7*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
H(3,4)=cH0*(DN(1,0)*Nenr[3] + DNenr(3,0)*cH8);
H(3,5)=cH0*(DN(1,1)*Nenr[3] + DNenr(3,1)*cH8);
H(3,6)=cH0*(DN(1,2)*Nenr[3] + DNenr(3,2)*cH8);
H(3,7)=cH7*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
H(3,8)=cH0*(DN(2,0)*Nenr[3] + DNenr(3,0)*cH9);
H(3,9)=cH0*(DN(2,1)*Nenr[3] + DNenr(3,1)*cH9);
H(3,10)=cH0*(DN(2,2)*Nenr[3] + DNenr(3,2)*cH9);
H(3,11)=cH7*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
H(3,12)=cH0*(DN(3,0)*Nenr[3] + DNenr(3,0)*cH10);
H(3,13)=cH0*(DN(3,1)*Nenr[3] + DNenr(3,1)*cH10);
H(3,14)=cH0*(DN(3,2)*Nenr[3] + DNenr(3,2)*cH10);
H(3,15)=cH7*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 = dt*not_stabilization_cut_elements_momentum/rho;
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
const double crhs_ee2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee3 = v(0,0) - vn(0,0);
const double crhs_ee4 = crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 = v(1,0) - vn(1,0);
const double crhs_ee6 = crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 = v(2,0) - vn(2,0);
const double crhs_ee8 = crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 = v(3,0) - vn(3,0);
const double crhs_ee10 = crhs_ee0*crhs_ee9 + vn(3,0);
const double crhs_ee11 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee12 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee13 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee14 = 1.0/(crhs_ee13 + 0.5);
const double crhs_ee15 = crhs_ee14/dt;
const double crhs_ee16 = crhs_ee14*(crhs_ee13 - 0.5);
const double crhs_ee17 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee18 = dt*not_stabilization_cut_elements_momentum/rho;
const double crhs_ee19 = crhs_ee18*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs_ee0*(f(3,0) - fn(3,0)) + fn(3,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee17*(acceleration_alpha_method(0,0)*crhs_ee16 - acceleration_alpha_method(0,0) + crhs_ee15*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee17*(acceleration_alpha_method(1,0)*crhs_ee16 - acceleration_alpha_method(1,0) + crhs_ee15*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee17*(acceleration_alpha_method(2,0)*crhs_ee16 - acceleration_alpha_method(2,0) + crhs_ee15*crhs_ee7)) + N[3]*(acceleration_alpha_method(3,0) + crhs_ee17*(acceleration_alpha_method(3,0)*crhs_ee16 - acceleration_alpha_method(3,0) + crhs_ee15*crhs_ee9)) + crhs_ee11*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8 + DN(3,1)*crhs_ee10) + crhs_ee12*(DN(0,2)*crhs_ee4 + DN(1,2)*crhs_ee6 + DN(2,2)*crhs_ee8 + DN(3,2)*crhs_ee10) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8 + DN(3,0)*crhs_ee10)));
const double crhs_ee20 = v(0,1) - vn(0,1);
const double crhs_ee21 = crhs_ee0*crhs_ee20 + vn(0,1);
const double crhs_ee22 = v(1,1) - vn(1,1);
const double crhs_ee23 = crhs_ee0*crhs_ee22 + vn(1,1);
const double crhs_ee24 = v(2,1) - vn(2,1);
const double crhs_ee25 = crhs_ee0*crhs_ee24 + vn(2,1);
const double crhs_ee26 = v(3,1) - vn(3,1);
const double crhs_ee27 = crhs_ee0*crhs_ee26 + vn(3,1);
const double crhs_ee28 = crhs_ee18*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs_ee0*(f(3,1) - fn(3,1)) + fn(3,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee17*(acceleration_alpha_method(0,1)*crhs_ee16 - acceleration_alpha_method(0,1) + crhs_ee15*crhs_ee20)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee17*(acceleration_alpha_method(1,1)*crhs_ee16 - acceleration_alpha_method(1,1) + crhs_ee15*crhs_ee22)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee17*(acceleration_alpha_method(2,1)*crhs_ee16 - acceleration_alpha_method(2,1) + crhs_ee15*crhs_ee24)) + N[3]*(acceleration_alpha_method(3,1) + crhs_ee17*(acceleration_alpha_method(3,1)*crhs_ee16 - acceleration_alpha_method(3,1) + crhs_ee15*crhs_ee26)) + crhs_ee11*(DN(0,1)*crhs_ee21 + DN(1,1)*crhs_ee23 + DN(2,1)*crhs_ee25 + DN(3,1)*crhs_ee27) + crhs_ee12*(DN(0,2)*crhs_ee21 + DN(1,2)*crhs_ee23 + DN(2,2)*crhs_ee25 + DN(3,2)*crhs_ee27) + crhs_ee2*(DN(0,0)*crhs_ee21 + DN(1,0)*crhs_ee23 + DN(2,0)*crhs_ee25 + DN(3,0)*crhs_ee27)));
const double crhs_ee29 = v(0,2) - vn(0,2);
const double crhs_ee30 = crhs_ee0*crhs_ee29 + vn(0,2);
const double crhs_ee31 = v(1,2) - vn(1,2);
const double crhs_ee32 = crhs_ee0*crhs_ee31 + vn(1,2);
const double crhs_ee33 = v(2,2) - vn(2,2);
const double crhs_ee34 = crhs_ee0*crhs_ee33 + vn(2,2);
const double crhs_ee35 = v(3,2) - vn(3,2);
const double crhs_ee36 = crhs_ee0*crhs_ee35 + vn(3,2);
const double crhs_ee37 = crhs_ee18*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs_ee0*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs_ee0*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs_ee0*(f(3,2) - fn(3,2)) + fn(3,2))) + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs_ee17*(acceleration_alpha_method(0,2)*crhs_ee16 - acceleration_alpha_method(0,2) + crhs_ee15*crhs_ee29)) + N[1]*(acceleration_alpha_method(1,2) + crhs_ee17*(acceleration_alpha_method(1,2)*crhs_ee16 - acceleration_alpha_method(1,2) + crhs_ee15*crhs_ee31)) + N[2]*(acceleration_alpha_method(2,2) + crhs_ee17*(acceleration_alpha_method(2,2)*crhs_ee16 - acceleration_alpha_method(2,2) + crhs_ee15*crhs_ee33)) + N[3]*(acceleration_alpha_method(3,2) + crhs_ee17*(acceleration_alpha_method(3,2)*crhs_ee16 - acceleration_alpha_method(3,2) + crhs_ee15*crhs_ee35)) + crhs_ee11*(DN(0,1)*crhs_ee30 + DN(1,1)*crhs_ee32 + DN(2,1)*crhs_ee34 + DN(3,1)*crhs_ee36) + crhs_ee12*(DN(0,2)*crhs_ee30 + DN(1,2)*crhs_ee32 + DN(2,2)*crhs_ee34 + DN(3,2)*crhs_ee36) + crhs_ee2*(DN(0,0)*crhs_ee30 + DN(1,0)*crhs_ee32 + DN(2,0)*crhs_ee34 + DN(3,0)*crhs_ee36)));
rhs_ee[0]=-DNenr(0,0)*crhs_ee19 - DNenr(0,1)*crhs_ee28 - DNenr(0,2)*crhs_ee37 + Nenr[0]*crhs_ee1;
rhs_ee[1]=-DNenr(1,0)*crhs_ee19 - DNenr(1,1)*crhs_ee28 - DNenr(1,2)*crhs_ee37 + Nenr[1]*crhs_ee1;
rhs_ee[2]=-DNenr(2,0)*crhs_ee19 - DNenr(2,1)*crhs_ee28 - DNenr(2,2)*crhs_ee37 + Nenr[2]*crhs_ee1;
rhs_ee[3]=-DNenr(3,0)*crhs_ee19 - DNenr(3,1)*crhs_ee28 - DNenr(3,2)*crhs_ee37 + Nenr[3]*crhs_ee1;


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
