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
#include "data_containers/two_fluid_navier_stokes_alpha_method/two_fluid_navier_stokes_alpha_method_data.h"

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

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Create new temporary data container
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const unsigned int number_of_gauss_points = gauss_weights.size();

    if (rOutput.size() != number_of_gauss_points){
        rOutput.resize(number_of_gauss_points);
    }

    if (rVariable == ARTIFICIAL_DYNAMIC_VISCOSITY){
        // Iterate over integration points to evaluate the artificial viscosity at each Gauss point
        for (unsigned int g = 0; g < number_of_gauss_points; ++g){
            this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            rOutput[g] = CalculateArtificialDynamicViscositySpecialization(data);
        }
    }
    else{
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <>
double TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::CalculateArtificialDynamicViscositySpecialization(TwoFluidNavierStokesAlphaMethodData<2, 3> &rData) const
{
    // Variables for artificial viscosity calculation
    double artificial_mu = 0.0;
    const double rho = rData.Density;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const auto &acceleration_alpha_method = rData.AccelerationAlphaMethod;
    const double max_spectral_radius = rData.MaxSpectralRadius;
    const double alpha_f = 1 / (1 + max_spectral_radius);
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const BoundedMatrix<double, 3, 2> vconv = (vn - vmeshn) + alpha_f * ((v - vmesh) - (vn - vmeshn));
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const double art_dyn_visc_coeff = 0.8;
    double grad_v_norm = 0.0;

    // Check that velocity gradient norm is non-zero

    grad_v_norm=sqrt(pow(fabs(DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1))), 2));


    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        artificial_mu=0.5*art_dyn_visc_coeff*h*sqrt(pow(fabs(-DN(0,0)*p[0] - DN(1,0)*p[1] - DN(2,0)*p[2] - rho*((DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0)) + (DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1))) - rho*(N[0]*(acceleration_alpha_method(0,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(0,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(0,0) + (v(0,0) - vn(0,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[1]*(acceleration_alpha_method(1,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(1,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(1,0) + (v(1,0) - vn(1,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[2]*(acceleration_alpha_method(2,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(2,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(2,0) + (v(2,0) - vn(2,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1))) + rho*(N[0]*(fn(0,0) + (f(0,0) - fn(0,0))/(max_spectral_radius + 1)) + N[1]*(fn(1,0) + (f(1,0) - fn(1,0))/(max_spectral_radius + 1)) + N[2]*(fn(2,0) + (f(2,0) - fn(2,0))/(max_spectral_radius + 1)))), 2) + pow(fabs(-DN(0,1)*p[0] - DN(1,1)*p[1] - DN(2,1)*p[2] - rho*((DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0)) + (DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1))) - rho*(N[0]*(acceleration_alpha_method(0,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(0,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(0,1) + (v(0,1) - vn(0,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[1]*(acceleration_alpha_method(1,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(1,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(1,1) + (v(1,1) - vn(1,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[2]*(acceleration_alpha_method(2,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(2,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(2,1) + (v(2,1) - vn(2,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1))) + rho*(N[0]*(fn(0,1) + (f(0,1) - fn(0,1))/(max_spectral_radius + 1)) + N[1]*(fn(1,1) + (f(1,1) - fn(1,1))/(max_spectral_radius + 1)) + N[2]*(fn(2,1) + (f(2,1) - fn(2,1))/(max_spectral_radius + 1)))), 2))/sqrt(pow(fabs(DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1))), 2));

    }

    return artificial_mu;
}

template <>
double TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::CalculateArtificialDynamicViscositySpecialization(TwoFluidNavierStokesAlphaMethodData<3, 4> &rData) const
{
    // Variables for artificial viscosity calculation
    double artificial_mu = 0.0;
    const double rho = rData.Density;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const auto &acceleration_alpha_method = rData.AccelerationAlphaMethod;
    const double max_spectral_radius = rData.MaxSpectralRadius;
    const double alpha_f = 1 / (1 + max_spectral_radius);
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const BoundedMatrix<double, 4, 3> vconv = (vn - vmeshn) + alpha_f * ((v - vmesh) - (vn - vmeshn));
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const double art_dyn_visc_coeff = 0.8;
    double grad_v_norm = 0.0;

    // Check that velocity gradient norm is non-zero

    grad_v_norm=sqrt(pow(fabs(DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2));


    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        artificial_mu=0.5*art_dyn_visc_coeff*h*sqrt(pow(fabs(-DN(0,0)*p[0] - DN(1,0)*p[1] - DN(2,0)*p[2] - DN(3,0)*p[3] - rho*((DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1)))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1)))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1)))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*(N[0]*(acceleration_alpha_method(0,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(0,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(0,0) + (v(0,0) - vn(0,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[1]*(acceleration_alpha_method(1,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(1,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(1,0) + (v(1,0) - vn(1,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[2]*(acceleration_alpha_method(2,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(2,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(2,0) + (v(2,0) - vn(2,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[3]*(acceleration_alpha_method(3,0) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(3,0)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(3,0) + (v(3,0) - vn(3,0))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1))) + rho*(N[0]*(fn(0,0) + (f(0,0) - fn(0,0))/(max_spectral_radius + 1)) + N[1]*(fn(1,0) + (f(1,0) - fn(1,0))/(max_spectral_radius + 1)) + N[2]*(fn(2,0) + (f(2,0) - fn(2,0))/(max_spectral_radius + 1)) + N[3]*(fn(3,0) + (f(3,0) - fn(3,0))/(max_spectral_radius + 1)))), 2) + pow(fabs(-DN(0,1)*p[0] - DN(1,1)*p[1] - DN(2,1)*p[2] - DN(3,1)*p[3] - rho*((DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1)))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1)))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1)))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*(N[0]*(acceleration_alpha_method(0,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(0,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(0,1) + (v(0,1) - vn(0,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[1]*(acceleration_alpha_method(1,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(1,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(1,1) + (v(1,1) - vn(1,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[2]*(acceleration_alpha_method(2,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(2,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(2,1) + (v(2,1) - vn(2,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[3]*(acceleration_alpha_method(3,1) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(3,1)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(3,1) + (v(3,1) - vn(3,1))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1))) + rho*(N[0]*(fn(0,1) + (f(0,1) - fn(0,1))/(max_spectral_radius + 1)) + N[1]*(fn(1,1) + (f(1,1) - fn(1,1))/(max_spectral_radius + 1)) + N[2]*(fn(2,1) + (f(2,1) - fn(2,1))/(max_spectral_radius + 1)) + N[3]*(fn(3,1) + (f(3,1) - fn(3,1))/(max_spectral_radius + 1)))), 2) + pow(fabs(-DN(0,2)*p[0] - DN(1,2)*p[1] - DN(2,2)*p[2] - DN(3,2)*p[3] - rho*((DN(0,0)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1)))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1)))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1)))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*(N[0]*(acceleration_alpha_method(0,2) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(0,2)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(0,2) + (v(0,2) - vn(0,2))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[1]*(acceleration_alpha_method(1,2) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(1,2)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(1,2) + (v(1,2) - vn(1,2))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[2]*(acceleration_alpha_method(2,2) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(2,2)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(2,2) + (v(2,2) - vn(2,2))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1)) + N[3]*(acceleration_alpha_method(3,2) + 0.5*(3 - max_spectral_radius)*(acceleration_alpha_method(3,2)*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) - 0.5 - 1/(max_spectral_radius + 1))/(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1)) - acceleration_alpha_method(3,2) + (v(3,2) - vn(3,2))/(dt*(0.5*(3 - max_spectral_radius)/(max_spectral_radius + 1) + 0.5 - 1/(max_spectral_radius + 1))))/(max_spectral_radius + 1))) + rho*(N[0]*(fn(0,2) + (f(0,2) - fn(0,2))/(max_spectral_radius + 1)) + N[1]*(fn(1,2) + (f(1,2) - fn(1,2))/(max_spectral_radius + 1)) + N[2]*(fn(2,2) + (f(2,2) - fn(2,2))/(max_spectral_radius + 1)) + N[3]*(fn(3,2) + (f(3,2) - fn(3,2))/(max_spectral_radius + 1)))), 2))/sqrt(pow(fabs(DN(0,0)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,0)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,0)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,0)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,0)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,1)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,1)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,1)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,1)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,0) + (v(0,0) - vn(0,0))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,0) + (v(1,0) - vn(1,0))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,0) + (v(2,0) - vn(2,0))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,0) + (v(3,0) - vn(3,0))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,1) + (v(0,1) - vn(0,1))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,1) + (v(1,1) - vn(1,1))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,1) + (v(2,1) - vn(2,1))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,1) + (v(3,1) - vn(3,1))/(max_spectral_radius + 1))), 2) + pow(fabs(DN(0,2)*(vn(0,2) + (v(0,2) - vn(0,2))/(max_spectral_radius + 1)) + DN(1,2)*(vn(1,2) + (v(1,2) - vn(1,2))/(max_spectral_radius + 1)) + DN(2,2)*(vn(2,2) + (v(2,2) - vn(2,2))/(max_spectral_radius + 1)) + DN(3,2)*(vn(3,2) + (v(3,2) - vn(3,2))/(max_spectral_radius + 1))), 2));

    }

    return artificial_mu;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodData<2, 3>& rData) const
{
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
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
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
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
    const double max_spectral_radius = rData.MaxSpectralRadius;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
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
const double clhs12 = clhs11*rho;
const double clhs13 = N[0]*clhs6;
const double clhs14 = 1.0/dt;
const double clhs15 = clhs6*(0.5*max_spectral_radius - 1.5);
const double clhs16 = 0.5*max_spectral_radius - 1.5;
const double clhs17 = clhs14*clhs16/(-clhs15 - clhs6 + 0.5);
const double clhs18 = -N[0]*clhs17 + clhs11;
const double clhs19 = clhs14*rho;
const double clhs20 = 1.0/(clhs19*dyn_tau + clhs3/h + mu*stab_c1/pow(h, 2));
const double clhs21 = clhs20*pow(rho, 2);
const double clhs22 = clhs21*clhs6;
const double clhs23 = clhs11*clhs22;
const double clhs24 = clhs16*clhs19/(clhs15 + clhs6 - 0.5);
const double clhs25 = clhs24*clhs6;
const double clhs26 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs27 = clhs21*clhs26;
const double clhs28 = clhs13*clhs27;
const double clhs29 = pow(N[0], 2)*clhs25 + clhs12*clhs13 + clhs18*clhs23 + clhs18*clhs28;
const double clhs30 = C(0,1)*DN(0,1) + clhs8;
const double clhs31 = C(1,2)*DN(0,1);
const double clhs32 = C(2,2)*DN(0,0) + clhs31;
const double clhs33 = DN(0,0)*clhs4;
const double clhs34 = DN(0,1)*clhs33;
const double clhs35 = clhs20*rho;
const double clhs36 = clhs26*clhs35;
const double clhs37 = clhs12*clhs20;
const double clhs38 = N[0]*clhs36 - N[0] + clhs37;
const double clhs39 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs40 = C(0,2)*DN(1,0);
const double clhs41 = C(2,2)*DN(1,1) + clhs40;
const double clhs42 = DN(0,0)*DN(1,0);
const double clhs43 = clhs13*clhs24;
const double clhs44 = N[1]*clhs43;
const double clhs45 = clhs4*clhs42 + clhs44;
const double clhs46 = DN(1,0)*clhs1 + DN(1,1)*clhs2;
const double clhs47 = clhs13*rho;
const double clhs48 = -N[1]*clhs17 + clhs46;
const double clhs49 = clhs23*clhs48 + clhs28*clhs48 + clhs46*clhs47;
const double clhs50 = C(0,1)*DN(1,1) + clhs40;
const double clhs51 = C(1,2)*DN(1,1);
const double clhs52 = C(2,2)*DN(1,0) + clhs51;
const double clhs53 = DN(1,1)*clhs33;
const double clhs54 = DN(0,0)*N[1];
const double clhs55 = DN(1,0)*N[0];
const double clhs56 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs57 = C(0,2)*DN(2,0);
const double clhs58 = C(2,2)*DN(2,1) + clhs57;
const double clhs59 = DN(0,0)*DN(2,0);
const double clhs60 = N[2]*clhs43;
const double clhs61 = clhs4*clhs59 + clhs60;
const double clhs62 = DN(2,0)*clhs1 + DN(2,1)*clhs2;
const double clhs63 = -N[2]*clhs17 + clhs62;
const double clhs64 = clhs23*clhs63 + clhs28*clhs63 + clhs47*clhs62;
const double clhs65 = C(0,1)*DN(2,1) + clhs57;
const double clhs66 = C(1,2)*DN(2,1);
const double clhs67 = C(2,2)*DN(2,0) + clhs66;
const double clhs68 = DN(2,1)*clhs33;
const double clhs69 = DN(0,0)*N[2];
const double clhs70 = DN(2,0)*N[0];
const double clhs71 = C(0,1)*DN(0,0) + clhs31;
const double clhs72 = pow(DN(0,1), 2);
const double clhs73 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs74 = C(0,1)*DN(1,0) + clhs51;
const double clhs75 = DN(0,1)*clhs4;
const double clhs76 = DN(1,0)*clhs75;
const double clhs77 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs78 = DN(0,1)*DN(1,1);
const double clhs79 = clhs4*clhs78 + clhs44;
const double clhs80 = DN(0,1)*N[1];
const double clhs81 = DN(1,1)*N[0];
const double clhs82 = C(0,1)*DN(2,0) + clhs66;
const double clhs83 = DN(2,0)*clhs75;
const double clhs84 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs85 = DN(0,1)*DN(2,1);
const double clhs86 = clhs4*clhs85 + clhs60;
const double clhs87 = DN(0,1)*N[2];
const double clhs88 = DN(2,1)*N[0];
const double clhs89 = clhs18*clhs35;
const double clhs90 = N[0] + clhs89;
const double clhs91 = clhs35*clhs48;
const double clhs92 = clhs20*(clhs42 + clhs78);
const double clhs93 = clhs35*clhs63;
const double clhs94 = clhs20*(clhs59 + clhs85);
const double clhs95 = DN(1,0)*clhs6;
const double clhs96 = DN(1,1)*clhs6;
const double clhs97 = N[1]*clhs6;
const double clhs98 = clhs22*clhs46;
const double clhs99 = clhs27*clhs97;
const double clhs100 = clhs12*clhs97 + clhs18*clhs98 + clhs18*clhs99;
const double clhs101 = clhs35*clhs46;
const double clhs102 = pow(DN(1,0), 2);
const double clhs103 = clhs97*rho;
const double clhs104 = pow(N[1], 2)*clhs25 + clhs103*clhs46 + clhs48*clhs98 + clhs48*clhs99;
const double clhs105 = DN(1,0)*clhs4;
const double clhs106 = DN(1,1)*clhs105;
const double clhs107 = N[1]*clhs36 - N[1] + clhs101;
const double clhs108 = DN(1,0)*DN(2,0);
const double clhs109 = N[2]*clhs24*clhs97;
const double clhs110 = clhs108*clhs4 + clhs109;
const double clhs111 = clhs103*clhs62 + clhs63*clhs98 + clhs63*clhs99;
const double clhs112 = DN(2,1)*clhs105;
const double clhs113 = DN(1,0)*N[2];
const double clhs114 = DN(2,0)*N[1];
const double clhs115 = pow(DN(1,1), 2);
const double clhs116 = DN(2,0)*clhs4;
const double clhs117 = DN(1,1)*clhs116;
const double clhs118 = DN(1,1)*DN(2,1);
const double clhs119 = clhs109 + clhs118*clhs4;
const double clhs120 = DN(1,1)*N[2];
const double clhs121 = DN(2,1)*N[1];
const double clhs122 = N[1] + clhs91;
const double clhs123 = clhs20*(clhs108 + clhs118);
const double clhs124 = DN(2,0)*clhs6;
const double clhs125 = DN(2,1)*clhs6;
const double clhs126 = N[2]*clhs6;
const double clhs127 = clhs22*clhs62;
const double clhs128 = clhs126*clhs27;
const double clhs129 = clhs12*clhs126 + clhs127*clhs18 + clhs128*clhs18;
const double clhs130 = clhs35*clhs62;
const double clhs131 = clhs126*rho;
const double clhs132 = clhs127*clhs48 + clhs128*clhs48 + clhs131*clhs46;
const double clhs133 = pow(DN(2,0), 2);
const double clhs134 = pow(N[2], 2)*clhs25 + clhs127*clhs63 + clhs128*clhs63 + clhs131*clhs62;
const double clhs135 = DN(2,1)*clhs116;
const double clhs136 = N[2]*clhs36 - N[2] + clhs130;
const double clhs137 = pow(DN(2,1), 2);
const double clhs138 = N[2] + clhs93;
lhs(0,0)=clhs0*clhs4 + clhs10*clhs9 + clhs29 + clhs5*clhs7;
lhs(0,1)=clhs10*clhs32 + clhs30*clhs7 + clhs34;
lhs(0,2)=DN(0,0)*clhs38;
lhs(0,3)=clhs10*clhs41 + clhs39*clhs7 + clhs45 + clhs49;
lhs(0,4)=clhs10*clhs52 + clhs50*clhs7 + clhs53;
lhs(0,5)=DN(1,0)*clhs37 + clhs36*clhs55 - clhs54;
lhs(0,6)=clhs10*clhs58 + clhs56*clhs7 + clhs61 + clhs64;
lhs(0,7)=clhs10*clhs67 + clhs65*clhs7 + clhs68;
lhs(0,8)=DN(2,0)*clhs37 + clhs36*clhs70 - clhs69;
lhs(1,0)=clhs10*clhs71 + clhs34 + clhs7*clhs9;
lhs(1,1)=clhs10*clhs73 + clhs29 + clhs32*clhs7 + clhs4*clhs72;
lhs(1,2)=DN(0,1)*clhs38;
lhs(1,3)=clhs10*clhs74 + clhs41*clhs7 + clhs76;
lhs(1,4)=clhs10*clhs77 + clhs49 + clhs52*clhs7 + clhs79;
lhs(1,5)=DN(1,1)*clhs37 + clhs36*clhs81 - clhs80;
lhs(1,6)=clhs10*clhs82 + clhs58*clhs7 + clhs83;
lhs(1,7)=clhs10*clhs84 + clhs64 + clhs67*clhs7 + clhs86;
lhs(1,8)=DN(2,1)*clhs37 + clhs36*clhs88 - clhs87;
lhs(2,0)=clhs7*clhs90;
lhs(2,1)=clhs10*clhs90;
lhs(2,2)=clhs20*(clhs0 + clhs72);
lhs(2,3)=clhs6*(DN(0,0)*clhs91 + clhs55);
lhs(2,4)=clhs6*(DN(0,1)*clhs91 + clhs81);
lhs(2,5)=clhs92;
lhs(2,6)=clhs6*(DN(0,0)*clhs93 + clhs70);
lhs(2,7)=clhs6*(DN(0,1)*clhs93 + clhs88);
lhs(2,8)=clhs94;
lhs(3,0)=clhs100 + clhs45 + clhs5*clhs95 + clhs9*clhs96;
lhs(3,1)=clhs30*clhs95 + clhs32*clhs96 + clhs76;
lhs(3,2)=DN(0,0)*clhs101 + clhs36*clhs54 - clhs55;
lhs(3,3)=clhs102*clhs4 + clhs104 + clhs39*clhs95 + clhs41*clhs96;
lhs(3,4)=clhs106 + clhs50*clhs95 + clhs52*clhs96;
lhs(3,5)=DN(1,0)*clhs107;
lhs(3,6)=clhs110 + clhs111 + clhs56*clhs95 + clhs58*clhs96;
lhs(3,7)=clhs112 + clhs65*clhs95 + clhs67*clhs96;
lhs(3,8)=DN(2,0)*clhs101 - clhs113 + clhs114*clhs36;
lhs(4,0)=clhs53 + clhs71*clhs96 + clhs9*clhs95;
lhs(4,1)=clhs100 + clhs32*clhs95 + clhs73*clhs96 + clhs79;
lhs(4,2)=DN(0,1)*clhs101 + clhs36*clhs80 - clhs81;
lhs(4,3)=clhs106 + clhs41*clhs95 + clhs74*clhs96;
lhs(4,4)=clhs104 + clhs115*clhs4 + clhs52*clhs95 + clhs77*clhs96;
lhs(4,5)=DN(1,1)*clhs107;
lhs(4,6)=clhs117 + clhs58*clhs95 + clhs82*clhs96;
lhs(4,7)=clhs111 + clhs119 + clhs67*clhs95 + clhs84*clhs96;
lhs(4,8)=DN(2,1)*clhs101 - clhs120 + clhs121*clhs36;
lhs(5,0)=clhs6*(DN(1,0)*clhs89 + clhs54);
lhs(5,1)=clhs6*(DN(1,1)*clhs89 + clhs80);
lhs(5,2)=clhs92;
lhs(5,3)=clhs122*clhs95;
lhs(5,4)=clhs122*clhs96;
lhs(5,5)=clhs20*(clhs102 + clhs115);
lhs(5,6)=clhs6*(DN(1,0)*clhs93 + clhs114);
lhs(5,7)=clhs6*(DN(1,1)*clhs93 + clhs121);
lhs(5,8)=clhs123;
lhs(6,0)=clhs124*clhs5 + clhs125*clhs9 + clhs129 + clhs61;
lhs(6,1)=clhs124*clhs30 + clhs125*clhs32 + clhs83;
lhs(6,2)=DN(0,0)*clhs130 + clhs36*clhs69 - clhs70;
lhs(6,3)=clhs110 + clhs124*clhs39 + clhs125*clhs41 + clhs132;
lhs(6,4)=clhs117 + clhs124*clhs50 + clhs125*clhs52;
lhs(6,5)=DN(1,0)*clhs130 + clhs113*clhs36 - clhs114;
lhs(6,6)=clhs124*clhs56 + clhs125*clhs58 + clhs133*clhs4 + clhs134;
lhs(6,7)=clhs124*clhs65 + clhs125*clhs67 + clhs135;
lhs(6,8)=DN(2,0)*clhs136;
lhs(7,0)=clhs124*clhs9 + clhs125*clhs71 + clhs68;
lhs(7,1)=clhs124*clhs32 + clhs125*clhs73 + clhs129 + clhs86;
lhs(7,2)=DN(0,1)*clhs130 + clhs36*clhs87 - clhs88;
lhs(7,3)=clhs112 + clhs124*clhs41 + clhs125*clhs74;
lhs(7,4)=clhs119 + clhs124*clhs52 + clhs125*clhs77 + clhs132;
lhs(7,5)=DN(1,1)*clhs130 + clhs120*clhs36 - clhs121;
lhs(7,6)=clhs124*clhs58 + clhs125*clhs82 + clhs135;
lhs(7,7)=clhs124*clhs67 + clhs125*clhs84 + clhs134 + clhs137*clhs4;
lhs(7,8)=DN(2,1)*clhs136;
lhs(8,0)=clhs6*(DN(2,0)*clhs89 + clhs69);
lhs(8,1)=clhs6*(DN(2,1)*clhs89 + clhs87);
lhs(8,2)=clhs94;
lhs(8,3)=clhs6*(DN(2,0)*clhs91 + clhs113);
lhs(8,4)=clhs6*(DN(2,1)*clhs91 + clhs120);
lhs(8,5)=clhs123;
lhs(8,6)=clhs124*clhs138;
lhs(8,7)=clhs125*clhs138;
lhs(8,8)=clhs20*(clhs133 + clhs137);


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
    const double max_spectral_radius = rData.MaxSpectralRadius;
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
const double clhs16 = clhs15*rho;
const double clhs17 = N[0]*clhs7;
const double clhs18 = 1.0/dt;
const double clhs19 = clhs7*(0.5*max_spectral_radius - 1.5);
const double clhs20 = 0.5*max_spectral_radius - 1.5;
const double clhs21 = clhs18*clhs20/(-clhs19 - clhs7 + 0.5);
const double clhs22 = -N[0]*clhs21 + clhs15;
const double clhs23 = clhs18*rho;
const double clhs24 = 1.0/(clhs23*dyn_tau + clhs4/h + mu*stab_c1/pow(h, 2));
const double clhs25 = clhs24*pow(rho, 2);
const double clhs26 = clhs25*clhs7;
const double clhs27 = clhs15*clhs26;
const double clhs28 = clhs20*clhs23/(clhs19 + clhs7 - 0.5);
const double clhs29 = clhs28*clhs7;
const double clhs30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs31 = clhs25*clhs30;
const double clhs32 = clhs17*clhs31;
const double clhs33 = pow(N[0], 2)*clhs29 + clhs16*clhs17 + clhs22*clhs27 + clhs22*clhs32;
const double clhs34 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs9;
const double clhs35 = C(1,3)*DN(0,1);
const double clhs36 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs35;
const double clhs37 = C(3,5)*DN(0,0);
const double clhs38 = C(4,5)*DN(0,2);
const double clhs39 = C(1,5)*DN(0,1) + clhs37 + clhs38;
const double clhs40 = DN(0,0)*clhs5;
const double clhs41 = DN(0,1)*clhs40;
const double clhs42 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs12;
const double clhs43 = C(3,4)*DN(0,1);
const double clhs44 = C(2,3)*DN(0,2) + clhs37 + clhs43;
const double clhs45 = C(2,5)*DN(0,2);
const double clhs46 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs45;
const double clhs47 = DN(0,2)*clhs40;
const double clhs48 = clhs24*rho;
const double clhs49 = clhs30*clhs48;
const double clhs50 = clhs16*clhs24;
const double clhs51 = N[0]*clhs49 - N[0] + clhs50;
const double clhs52 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs53 = C(0,3)*DN(1,0);
const double clhs54 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs53;
const double clhs55 = C(0,5)*DN(1,0);
const double clhs56 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs55;
const double clhs57 = DN(0,0)*DN(1,0);
const double clhs58 = clhs17*clhs28;
const double clhs59 = N[1]*clhs58;
const double clhs60 = clhs5*clhs57 + clhs59;
const double clhs61 = DN(1,0)*clhs1 + DN(1,1)*clhs2 + DN(1,2)*clhs3;
const double clhs62 = clhs17*rho;
const double clhs63 = -N[1]*clhs21 + clhs61;
const double clhs64 = clhs27*clhs63 + clhs32*clhs63 + clhs61*clhs62;
const double clhs65 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs53;
const double clhs66 = C(1,3)*DN(1,1);
const double clhs67 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs66;
const double clhs68 = C(3,5)*DN(1,0);
const double clhs69 = C(4,5)*DN(1,2);
const double clhs70 = C(1,5)*DN(1,1) + clhs68 + clhs69;
const double clhs71 = DN(1,1)*clhs40;
const double clhs72 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs55;
const double clhs73 = C(3,4)*DN(1,1);
const double clhs74 = C(2,3)*DN(1,2) + clhs68 + clhs73;
const double clhs75 = C(2,5)*DN(1,2);
const double clhs76 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs75;
const double clhs77 = DN(1,2)*clhs40;
const double clhs78 = DN(0,0)*N[1];
const double clhs79 = DN(1,0)*N[0];
const double clhs80 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs81 = C(0,3)*DN(2,0);
const double clhs82 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs81;
const double clhs83 = C(0,5)*DN(2,0);
const double clhs84 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs83;
const double clhs85 = DN(0,0)*DN(2,0);
const double clhs86 = N[2]*clhs58;
const double clhs87 = clhs5*clhs85 + clhs86;
const double clhs88 = DN(2,0)*clhs1 + DN(2,1)*clhs2 + DN(2,2)*clhs3;
const double clhs89 = -N[2]*clhs21 + clhs88;
const double clhs90 = clhs27*clhs89 + clhs32*clhs89 + clhs62*clhs88;
const double clhs91 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs81;
const double clhs92 = C(1,3)*DN(2,1);
const double clhs93 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs92;
const double clhs94 = C(3,5)*DN(2,0);
const double clhs95 = C(4,5)*DN(2,2);
const double clhs96 = C(1,5)*DN(2,1) + clhs94 + clhs95;
const double clhs97 = DN(2,1)*clhs40;
const double clhs98 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs83;
const double clhs99 = C(3,4)*DN(2,1);
const double clhs100 = C(2,3)*DN(2,2) + clhs94 + clhs99;
const double clhs101 = C(2,5)*DN(2,2);
const double clhs102 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs101;
const double clhs103 = DN(2,2)*clhs40;
const double clhs104 = DN(0,0)*N[2];
const double clhs105 = DN(2,0)*N[0];
const double clhs106 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs107 = C(0,3)*DN(3,0);
const double clhs108 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs107;
const double clhs109 = C(0,5)*DN(3,0);
const double clhs110 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs109;
const double clhs111 = DN(0,0)*DN(3,0);
const double clhs112 = N[3]*clhs58;
const double clhs113 = clhs111*clhs5 + clhs112;
const double clhs114 = DN(3,0)*clhs1 + DN(3,1)*clhs2 + DN(3,2)*clhs3;
const double clhs115 = -N[3]*clhs21 + clhs114;
const double clhs116 = clhs114*clhs62 + clhs115*clhs27 + clhs115*clhs32;
const double clhs117 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs107;
const double clhs118 = C(1,3)*DN(3,1);
const double clhs119 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs118;
const double clhs120 = C(3,5)*DN(3,0);
const double clhs121 = C(4,5)*DN(3,2);
const double clhs122 = C(1,5)*DN(3,1) + clhs120 + clhs121;
const double clhs123 = DN(3,1)*clhs40;
const double clhs124 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs109;
const double clhs125 = C(3,4)*DN(3,1);
const double clhs126 = C(2,3)*DN(3,2) + clhs120 + clhs125;
const double clhs127 = C(2,5)*DN(3,2);
const double clhs128 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs127;
const double clhs129 = DN(3,2)*clhs40;
const double clhs130 = DN(0,0)*N[3];
const double clhs131 = DN(3,0)*N[0];
const double clhs132 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs35;
const double clhs133 = C(0,4)*DN(0,0) + clhs38 + clhs43;
const double clhs134 = pow(DN(0,1), 2);
const double clhs135 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs136 = C(1,4)*DN(0,1);
const double clhs137 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs136;
const double clhs138 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs136;
const double clhs139 = C(2,4)*DN(0,2);
const double clhs140 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs139;
const double clhs141 = DN(0,1)*clhs5;
const double clhs142 = DN(0,2)*clhs141;
const double clhs143 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs66;
const double clhs144 = C(0,4)*DN(1,0) + clhs69 + clhs73;
const double clhs145 = DN(1,0)*clhs141;
const double clhs146 = DN(0,1)*DN(1,1);
const double clhs147 = clhs146*clhs5;
const double clhs148 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs149 = C(1,4)*DN(1,1);
const double clhs150 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs149;
const double clhs151 = clhs59 + clhs64;
const double clhs152 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs149;
const double clhs153 = C(2,4)*DN(1,2);
const double clhs154 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs153;
const double clhs155 = DN(1,2)*clhs141;
const double clhs156 = DN(0,1)*N[1];
const double clhs157 = DN(1,1)*N[0];
const double clhs158 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs92;
const double clhs159 = C(0,4)*DN(2,0) + clhs95 + clhs99;
const double clhs160 = DN(2,0)*clhs141;
const double clhs161 = DN(0,1)*DN(2,1);
const double clhs162 = clhs161*clhs5;
const double clhs163 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs164 = C(1,4)*DN(2,1);
const double clhs165 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs164;
const double clhs166 = clhs86 + clhs90;
const double clhs167 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs164;
const double clhs168 = C(2,4)*DN(2,2);
const double clhs169 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs168;
const double clhs170 = DN(2,2)*clhs141;
const double clhs171 = DN(0,1)*N[2];
const double clhs172 = DN(2,1)*N[0];
const double clhs173 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs118;
const double clhs174 = C(0,4)*DN(3,0) + clhs121 + clhs125;
const double clhs175 = DN(3,0)*clhs141;
const double clhs176 = DN(0,1)*DN(3,1);
const double clhs177 = clhs176*clhs5;
const double clhs178 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs179 = C(1,4)*DN(3,1);
const double clhs180 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs179;
const double clhs181 = clhs112 + clhs116;
const double clhs182 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs179;
const double clhs183 = C(2,4)*DN(3,2);
const double clhs184 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs183;
const double clhs185 = DN(3,2)*clhs141;
const double clhs186 = DN(0,1)*N[3];
const double clhs187 = DN(3,1)*N[0];
const double clhs188 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs45;
const double clhs189 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs139;
const double clhs190 = pow(DN(0,2), 2);
const double clhs191 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs192 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs75;
const double clhs193 = DN(0,2)*clhs5;
const double clhs194 = DN(1,0)*clhs193;
const double clhs195 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs153;
const double clhs196 = DN(1,1)*clhs193;
const double clhs197 = DN(0,2)*DN(1,2);
const double clhs198 = clhs197*clhs5;
const double clhs199 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs200 = DN(0,2)*N[1];
const double clhs201 = DN(1,2)*N[0];
const double clhs202 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs101;
const double clhs203 = DN(2,0)*clhs193;
const double clhs204 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs168;
const double clhs205 = DN(2,1)*clhs193;
const double clhs206 = DN(0,2)*DN(2,2);
const double clhs207 = clhs206*clhs5;
const double clhs208 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs209 = DN(0,2)*N[2];
const double clhs210 = DN(2,2)*N[0];
const double clhs211 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs127;
const double clhs212 = DN(3,0)*clhs193;
const double clhs213 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs183;
const double clhs214 = DN(3,1)*clhs193;
const double clhs215 = DN(0,2)*DN(3,2);
const double clhs216 = clhs215*clhs5;
const double clhs217 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs218 = DN(0,2)*N[3];
const double clhs219 = DN(3,2)*N[0];
const double clhs220 = clhs22*clhs48;
const double clhs221 = N[0] + clhs220;
const double clhs222 = clhs48*clhs63;
const double clhs223 = clhs24*(clhs146 + clhs197 + clhs57);
const double clhs224 = clhs48*clhs89;
const double clhs225 = clhs24*(clhs161 + clhs206 + clhs85);
const double clhs226 = clhs115*clhs48;
const double clhs227 = clhs24*(clhs111 + clhs176 + clhs215);
const double clhs228 = DN(1,0)*clhs7;
const double clhs229 = DN(1,1)*clhs7;
const double clhs230 = DN(1,2)*clhs7;
const double clhs231 = N[1]*clhs7;
const double clhs232 = clhs26*clhs61;
const double clhs233 = clhs231*clhs31;
const double clhs234 = clhs16*clhs231 + clhs22*clhs232 + clhs22*clhs233;
const double clhs235 = clhs48*clhs61;
const double clhs236 = pow(DN(1,0), 2);
const double clhs237 = clhs231*rho;
const double clhs238 = pow(N[1], 2)*clhs29 + clhs232*clhs63 + clhs233*clhs63 + clhs237*clhs61;
const double clhs239 = DN(1,0)*clhs5;
const double clhs240 = DN(1,1)*clhs239;
const double clhs241 = DN(1,2)*clhs239;
const double clhs242 = N[1]*clhs49 - N[1] + clhs235;
const double clhs243 = DN(1,0)*DN(2,0);
const double clhs244 = clhs231*clhs28;
const double clhs245 = N[2]*clhs244;
const double clhs246 = clhs243*clhs5 + clhs245;
const double clhs247 = clhs232*clhs89 + clhs233*clhs89 + clhs237*clhs88;
const double clhs248 = DN(2,1)*clhs239;
const double clhs249 = DN(2,2)*clhs239;
const double clhs250 = DN(1,0)*N[2];
const double clhs251 = DN(2,0)*N[1];
const double clhs252 = DN(1,0)*DN(3,0);
const double clhs253 = N[3]*clhs244;
const double clhs254 = clhs252*clhs5 + clhs253;
const double clhs255 = clhs114*clhs237 + clhs115*clhs232 + clhs115*clhs233;
const double clhs256 = DN(3,1)*clhs239;
const double clhs257 = DN(3,2)*clhs239;
const double clhs258 = DN(1,0)*N[3];
const double clhs259 = DN(3,0)*N[1];
const double clhs260 = clhs234 + clhs59;
const double clhs261 = pow(DN(1,1), 2);
const double clhs262 = DN(1,1)*clhs5;
const double clhs263 = DN(1,2)*clhs262;
const double clhs264 = DN(2,0)*clhs262;
const double clhs265 = DN(1,1)*DN(2,1);
const double clhs266 = clhs265*clhs5;
const double clhs267 = clhs245 + clhs247;
const double clhs268 = DN(2,2)*clhs262;
const double clhs269 = DN(1,1)*N[2];
const double clhs270 = DN(2,1)*N[1];
const double clhs271 = DN(3,0)*clhs262;
const double clhs272 = DN(1,1)*DN(3,1);
const double clhs273 = clhs272*clhs5;
const double clhs274 = clhs253 + clhs255;
const double clhs275 = DN(3,2)*clhs262;
const double clhs276 = DN(1,1)*N[3];
const double clhs277 = DN(3,1)*N[1];
const double clhs278 = pow(DN(1,2), 2);
const double clhs279 = DN(1,2)*clhs5;
const double clhs280 = DN(2,0)*clhs279;
const double clhs281 = DN(2,1)*clhs279;
const double clhs282 = DN(1,2)*DN(2,2);
const double clhs283 = clhs282*clhs5;
const double clhs284 = DN(1,2)*N[2];
const double clhs285 = DN(2,2)*N[1];
const double clhs286 = DN(3,0)*clhs279;
const double clhs287 = DN(3,1)*clhs279;
const double clhs288 = DN(1,2)*DN(3,2);
const double clhs289 = clhs288*clhs5;
const double clhs290 = DN(1,2)*N[3];
const double clhs291 = DN(3,2)*N[1];
const double clhs292 = N[1] + clhs222;
const double clhs293 = clhs24*(clhs243 + clhs265 + clhs282);
const double clhs294 = clhs24*(clhs252 + clhs272 + clhs288);
const double clhs295 = DN(2,0)*clhs7;
const double clhs296 = DN(2,1)*clhs7;
const double clhs297 = DN(2,2)*clhs7;
const double clhs298 = N[2]*clhs7;
const double clhs299 = clhs26*clhs88;
const double clhs300 = clhs298*clhs31;
const double clhs301 = clhs16*clhs298 + clhs22*clhs299 + clhs22*clhs300;
const double clhs302 = clhs48*clhs88;
const double clhs303 = clhs298*rho;
const double clhs304 = clhs299*clhs63 + clhs300*clhs63 + clhs303*clhs61;
const double clhs305 = pow(DN(2,0), 2);
const double clhs306 = pow(N[2], 2)*clhs29 + clhs299*clhs89 + clhs300*clhs89 + clhs303*clhs88;
const double clhs307 = DN(2,0)*clhs5;
const double clhs308 = DN(2,1)*clhs307;
const double clhs309 = DN(2,2)*clhs307;
const double clhs310 = N[2]*clhs49 - N[2] + clhs302;
const double clhs311 = DN(2,0)*DN(3,0);
const double clhs312 = N[3]*clhs28*clhs298;
const double clhs313 = clhs311*clhs5 + clhs312;
const double clhs314 = clhs114*clhs303 + clhs115*clhs299 + clhs115*clhs300;
const double clhs315 = DN(3,1)*clhs307;
const double clhs316 = DN(3,2)*clhs307;
const double clhs317 = DN(2,0)*N[3];
const double clhs318 = DN(3,0)*N[2];
const double clhs319 = clhs301 + clhs86;
const double clhs320 = clhs245 + clhs304;
const double clhs321 = pow(DN(2,1), 2);
const double clhs322 = DN(2,1)*clhs5;
const double clhs323 = DN(2,2)*clhs322;
const double clhs324 = DN(3,0)*clhs322;
const double clhs325 = DN(2,1)*DN(3,1);
const double clhs326 = clhs325*clhs5;
const double clhs327 = clhs312 + clhs314;
const double clhs328 = DN(3,2)*clhs322;
const double clhs329 = DN(2,1)*N[3];
const double clhs330 = DN(3,1)*N[2];
const double clhs331 = pow(DN(2,2), 2);
const double clhs332 = DN(2,2)*clhs5;
const double clhs333 = DN(3,0)*clhs332;
const double clhs334 = DN(3,1)*clhs332;
const double clhs335 = DN(2,2)*DN(3,2);
const double clhs336 = clhs335*clhs5;
const double clhs337 = DN(2,2)*N[3];
const double clhs338 = DN(3,2)*N[2];
const double clhs339 = N[2] + clhs224;
const double clhs340 = clhs24*(clhs311 + clhs325 + clhs335);
const double clhs341 = DN(3,0)*clhs7;
const double clhs342 = DN(3,1)*clhs7;
const double clhs343 = DN(3,2)*clhs7;
const double clhs344 = N[3]*clhs7;
const double clhs345 = clhs114*clhs26;
const double clhs346 = clhs31*clhs344;
const double clhs347 = clhs16*clhs344 + clhs22*clhs345 + clhs22*clhs346;
const double clhs348 = clhs114*clhs48;
const double clhs349 = clhs344*rho;
const double clhs350 = clhs345*clhs63 + clhs346*clhs63 + clhs349*clhs61;
const double clhs351 = clhs345*clhs89 + clhs346*clhs89 + clhs349*clhs88;
const double clhs352 = pow(DN(3,0), 2);
const double clhs353 = pow(N[3], 2)*clhs29 + clhs114*clhs349 + clhs115*clhs345 + clhs115*clhs346;
const double clhs354 = DN(3,0)*clhs5;
const double clhs355 = DN(3,1)*clhs354;
const double clhs356 = DN(3,2)*clhs354;
const double clhs357 = N[3]*clhs49 - N[3] + clhs348;
const double clhs358 = clhs112 + clhs347;
const double clhs359 = clhs253 + clhs350;
const double clhs360 = clhs312 + clhs351;
const double clhs361 = pow(DN(3,1), 2);
const double clhs362 = DN(3,1)*DN(3,2)*clhs5;
const double clhs363 = pow(DN(3,2), 2);
const double clhs364 = N[3] + clhs226;
lhs(0,0)=clhs0*clhs5 + clhs10*clhs11 + clhs13*clhs14 + clhs33 + clhs6*clhs8;
lhs(0,1)=clhs11*clhs36 + clhs14*clhs39 + clhs34*clhs8 + clhs41;
lhs(0,2)=clhs11*clhs44 + clhs14*clhs46 + clhs42*clhs8 + clhs47;
lhs(0,3)=DN(0,0)*clhs51;
lhs(0,4)=clhs11*clhs54 + clhs14*clhs56 + clhs52*clhs8 + clhs60 + clhs64;
lhs(0,5)=clhs11*clhs67 + clhs14*clhs70 + clhs65*clhs8 + clhs71;
lhs(0,6)=clhs11*clhs74 + clhs14*clhs76 + clhs72*clhs8 + clhs77;
lhs(0,7)=DN(1,0)*clhs50 + clhs49*clhs79 - clhs78;
lhs(0,8)=clhs11*clhs82 + clhs14*clhs84 + clhs8*clhs80 + clhs87 + clhs90;
lhs(0,9)=clhs11*clhs93 + clhs14*clhs96 + clhs8*clhs91 + clhs97;
lhs(0,10)=clhs100*clhs11 + clhs102*clhs14 + clhs103 + clhs8*clhs98;
lhs(0,11)=DN(2,0)*clhs50 - clhs104 + clhs105*clhs49;
lhs(0,12)=clhs106*clhs8 + clhs108*clhs11 + clhs110*clhs14 + clhs113 + clhs116;
lhs(0,13)=clhs11*clhs119 + clhs117*clhs8 + clhs122*clhs14 + clhs123;
lhs(0,14)=clhs11*clhs126 + clhs124*clhs8 + clhs128*clhs14 + clhs129;
lhs(0,15)=DN(3,0)*clhs50 - clhs130 + clhs131*clhs49;
lhs(1,0)=clhs10*clhs8 + clhs11*clhs132 + clhs133*clhs14 + clhs41;
lhs(1,1)=clhs11*clhs135 + clhs134*clhs5 + clhs137*clhs14 + clhs33 + clhs36*clhs8;
lhs(1,2)=clhs11*clhs138 + clhs14*clhs140 + clhs142 + clhs44*clhs8;
lhs(1,3)=DN(0,1)*clhs51;
lhs(1,4)=clhs11*clhs143 + clhs14*clhs144 + clhs145 + clhs54*clhs8;
lhs(1,5)=clhs11*clhs148 + clhs14*clhs150 + clhs147 + clhs151 + clhs67*clhs8;
lhs(1,6)=clhs11*clhs152 + clhs14*clhs154 + clhs155 + clhs74*clhs8;
lhs(1,7)=DN(1,1)*clhs50 - clhs156 + clhs157*clhs49;
lhs(1,8)=clhs11*clhs158 + clhs14*clhs159 + clhs160 + clhs8*clhs82;
lhs(1,9)=clhs11*clhs163 + clhs14*clhs165 + clhs162 + clhs166 + clhs8*clhs93;
lhs(1,10)=clhs100*clhs8 + clhs11*clhs167 + clhs14*clhs169 + clhs170;
lhs(1,11)=DN(2,1)*clhs50 - clhs171 + clhs172*clhs49;
lhs(1,12)=clhs108*clhs8 + clhs11*clhs173 + clhs14*clhs174 + clhs175;
lhs(1,13)=clhs11*clhs178 + clhs119*clhs8 + clhs14*clhs180 + clhs177 + clhs181;
lhs(1,14)=clhs11*clhs182 + clhs126*clhs8 + clhs14*clhs184 + clhs185;
lhs(1,15)=DN(3,1)*clhs50 - clhs186 + clhs187*clhs49;
lhs(2,0)=clhs11*clhs133 + clhs13*clhs8 + clhs14*clhs188 + clhs47;
lhs(2,1)=clhs11*clhs137 + clhs14*clhs189 + clhs142 + clhs39*clhs8;
lhs(2,2)=clhs11*clhs140 + clhs14*clhs191 + clhs190*clhs5 + clhs33 + clhs46*clhs8;
lhs(2,3)=DN(0,2)*clhs51;
lhs(2,4)=clhs11*clhs144 + clhs14*clhs192 + clhs194 + clhs56*clhs8;
lhs(2,5)=clhs11*clhs150 + clhs14*clhs195 + clhs196 + clhs70*clhs8;
lhs(2,6)=clhs11*clhs154 + clhs14*clhs199 + clhs151 + clhs198 + clhs76*clhs8;
lhs(2,7)=DN(1,2)*clhs50 - clhs200 + clhs201*clhs49;
lhs(2,8)=clhs11*clhs159 + clhs14*clhs202 + clhs203 + clhs8*clhs84;
lhs(2,9)=clhs11*clhs165 + clhs14*clhs204 + clhs205 + clhs8*clhs96;
lhs(2,10)=clhs102*clhs8 + clhs11*clhs169 + clhs14*clhs208 + clhs166 + clhs207;
lhs(2,11)=DN(2,2)*clhs50 - clhs209 + clhs210*clhs49;
lhs(2,12)=clhs11*clhs174 + clhs110*clhs8 + clhs14*clhs211 + clhs212;
lhs(2,13)=clhs11*clhs180 + clhs122*clhs8 + clhs14*clhs213 + clhs214;
lhs(2,14)=clhs11*clhs184 + clhs128*clhs8 + clhs14*clhs217 + clhs181 + clhs216;
lhs(2,15)=DN(3,2)*clhs50 - clhs218 + clhs219*clhs49;
lhs(3,0)=clhs221*clhs8;
lhs(3,1)=clhs11*clhs221;
lhs(3,2)=clhs14*clhs221;
lhs(3,3)=clhs24*(clhs0 + clhs134 + clhs190);
lhs(3,4)=clhs7*(DN(0,0)*clhs222 + clhs79);
lhs(3,5)=clhs7*(DN(0,1)*clhs222 + clhs157);
lhs(3,6)=clhs7*(DN(0,2)*clhs222 + clhs201);
lhs(3,7)=clhs223;
lhs(3,8)=clhs7*(DN(0,0)*clhs224 + clhs105);
lhs(3,9)=clhs7*(DN(0,1)*clhs224 + clhs172);
lhs(3,10)=clhs7*(DN(0,2)*clhs224 + clhs210);
lhs(3,11)=clhs225;
lhs(3,12)=clhs7*(DN(0,0)*clhs226 + clhs131);
lhs(3,13)=clhs7*(DN(0,1)*clhs226 + clhs187);
lhs(3,14)=clhs7*(DN(0,2)*clhs226 + clhs219);
lhs(3,15)=clhs227;
lhs(4,0)=clhs10*clhs229 + clhs13*clhs230 + clhs228*clhs6 + clhs234 + clhs60;
lhs(4,1)=clhs145 + clhs228*clhs34 + clhs229*clhs36 + clhs230*clhs39;
lhs(4,2)=clhs194 + clhs228*clhs42 + clhs229*clhs44 + clhs230*clhs46;
lhs(4,3)=DN(0,0)*clhs235 + clhs49*clhs78 - clhs79;
lhs(4,4)=clhs228*clhs52 + clhs229*clhs54 + clhs230*clhs56 + clhs236*clhs5 + clhs238;
lhs(4,5)=clhs228*clhs65 + clhs229*clhs67 + clhs230*clhs70 + clhs240;
lhs(4,6)=clhs228*clhs72 + clhs229*clhs74 + clhs230*clhs76 + clhs241;
lhs(4,7)=DN(1,0)*clhs242;
lhs(4,8)=clhs228*clhs80 + clhs229*clhs82 + clhs230*clhs84 + clhs246 + clhs247;
lhs(4,9)=clhs228*clhs91 + clhs229*clhs93 + clhs230*clhs96 + clhs248;
lhs(4,10)=clhs100*clhs229 + clhs102*clhs230 + clhs228*clhs98 + clhs249;
lhs(4,11)=DN(2,0)*clhs235 - clhs250 + clhs251*clhs49;
lhs(4,12)=clhs106*clhs228 + clhs108*clhs229 + clhs110*clhs230 + clhs254 + clhs255;
lhs(4,13)=clhs117*clhs228 + clhs119*clhs229 + clhs122*clhs230 + clhs256;
lhs(4,14)=clhs124*clhs228 + clhs126*clhs229 + clhs128*clhs230 + clhs257;
lhs(4,15)=DN(3,0)*clhs235 - clhs258 + clhs259*clhs49;
lhs(5,0)=clhs10*clhs228 + clhs132*clhs229 + clhs133*clhs230 + clhs71;
lhs(5,1)=clhs135*clhs229 + clhs137*clhs230 + clhs147 + clhs228*clhs36 + clhs260;
lhs(5,2)=clhs138*clhs229 + clhs140*clhs230 + clhs196 + clhs228*clhs44;
lhs(5,3)=DN(0,1)*clhs235 + clhs156*clhs49 - clhs157;
lhs(5,4)=clhs143*clhs229 + clhs144*clhs230 + clhs228*clhs54 + clhs240;
lhs(5,5)=clhs148*clhs229 + clhs150*clhs230 + clhs228*clhs67 + clhs238 + clhs261*clhs5;
lhs(5,6)=clhs152*clhs229 + clhs154*clhs230 + clhs228*clhs74 + clhs263;
lhs(5,7)=DN(1,1)*clhs242;
lhs(5,8)=clhs158*clhs229 + clhs159*clhs230 + clhs228*clhs82 + clhs264;
lhs(5,9)=clhs163*clhs229 + clhs165*clhs230 + clhs228*clhs93 + clhs266 + clhs267;
lhs(5,10)=clhs100*clhs228 + clhs167*clhs229 + clhs169*clhs230 + clhs268;
lhs(5,11)=DN(2,1)*clhs235 - clhs269 + clhs270*clhs49;
lhs(5,12)=clhs108*clhs228 + clhs173*clhs229 + clhs174*clhs230 + clhs271;
lhs(5,13)=clhs119*clhs228 + clhs178*clhs229 + clhs180*clhs230 + clhs273 + clhs274;
lhs(5,14)=clhs126*clhs228 + clhs182*clhs229 + clhs184*clhs230 + clhs275;
lhs(5,15)=DN(3,1)*clhs235 - clhs276 + clhs277*clhs49;
lhs(6,0)=clhs13*clhs228 + clhs133*clhs229 + clhs188*clhs230 + clhs77;
lhs(6,1)=clhs137*clhs229 + clhs155 + clhs189*clhs230 + clhs228*clhs39;
lhs(6,2)=clhs140*clhs229 + clhs191*clhs230 + clhs198 + clhs228*clhs46 + clhs260;
lhs(6,3)=DN(0,2)*clhs235 + clhs200*clhs49 - clhs201;
lhs(6,4)=clhs144*clhs229 + clhs192*clhs230 + clhs228*clhs56 + clhs241;
lhs(6,5)=clhs150*clhs229 + clhs195*clhs230 + clhs228*clhs70 + clhs263;
lhs(6,6)=clhs154*clhs229 + clhs199*clhs230 + clhs228*clhs76 + clhs238 + clhs278*clhs5;
lhs(6,7)=DN(1,2)*clhs242;
lhs(6,8)=clhs159*clhs229 + clhs202*clhs230 + clhs228*clhs84 + clhs280;
lhs(6,9)=clhs165*clhs229 + clhs204*clhs230 + clhs228*clhs96 + clhs281;
lhs(6,10)=clhs102*clhs228 + clhs169*clhs229 + clhs208*clhs230 + clhs267 + clhs283;
lhs(6,11)=DN(2,2)*clhs235 - clhs284 + clhs285*clhs49;
lhs(6,12)=clhs110*clhs228 + clhs174*clhs229 + clhs211*clhs230 + clhs286;
lhs(6,13)=clhs122*clhs228 + clhs180*clhs229 + clhs213*clhs230 + clhs287;
lhs(6,14)=clhs128*clhs228 + clhs184*clhs229 + clhs217*clhs230 + clhs274 + clhs289;
lhs(6,15)=DN(3,2)*clhs235 - clhs290 + clhs291*clhs49;
lhs(7,0)=clhs7*(DN(1,0)*clhs220 + clhs78);
lhs(7,1)=clhs7*(DN(1,1)*clhs220 + clhs156);
lhs(7,2)=clhs7*(DN(1,2)*clhs220 + clhs200);
lhs(7,3)=clhs223;
lhs(7,4)=clhs228*clhs292;
lhs(7,5)=clhs229*clhs292;
lhs(7,6)=clhs230*clhs292;
lhs(7,7)=clhs24*(clhs236 + clhs261 + clhs278);
lhs(7,8)=clhs7*(DN(1,0)*clhs224 + clhs251);
lhs(7,9)=clhs7*(DN(1,1)*clhs224 + clhs270);
lhs(7,10)=clhs7*(DN(1,2)*clhs224 + clhs285);
lhs(7,11)=clhs293;
lhs(7,12)=clhs7*(DN(1,0)*clhs226 + clhs259);
lhs(7,13)=clhs7*(DN(1,1)*clhs226 + clhs277);
lhs(7,14)=clhs7*(DN(1,2)*clhs226 + clhs291);
lhs(7,15)=clhs294;
lhs(8,0)=clhs10*clhs296 + clhs13*clhs297 + clhs295*clhs6 + clhs301 + clhs87;
lhs(8,1)=clhs160 + clhs295*clhs34 + clhs296*clhs36 + clhs297*clhs39;
lhs(8,2)=clhs203 + clhs295*clhs42 + clhs296*clhs44 + clhs297*clhs46;
lhs(8,3)=DN(0,0)*clhs302 + clhs104*clhs49 - clhs105;
lhs(8,4)=clhs246 + clhs295*clhs52 + clhs296*clhs54 + clhs297*clhs56 + clhs304;
lhs(8,5)=clhs264 + clhs295*clhs65 + clhs296*clhs67 + clhs297*clhs70;
lhs(8,6)=clhs280 + clhs295*clhs72 + clhs296*clhs74 + clhs297*clhs76;
lhs(8,7)=DN(1,0)*clhs302 + clhs250*clhs49 - clhs251;
lhs(8,8)=clhs295*clhs80 + clhs296*clhs82 + clhs297*clhs84 + clhs305*clhs5 + clhs306;
lhs(8,9)=clhs295*clhs91 + clhs296*clhs93 + clhs297*clhs96 + clhs308;
lhs(8,10)=clhs100*clhs296 + clhs102*clhs297 + clhs295*clhs98 + clhs309;
lhs(8,11)=DN(2,0)*clhs310;
lhs(8,12)=clhs106*clhs295 + clhs108*clhs296 + clhs110*clhs297 + clhs313 + clhs314;
lhs(8,13)=clhs117*clhs295 + clhs119*clhs296 + clhs122*clhs297 + clhs315;
lhs(8,14)=clhs124*clhs295 + clhs126*clhs296 + clhs128*clhs297 + clhs316;
lhs(8,15)=DN(3,0)*clhs302 - clhs317 + clhs318*clhs49;
lhs(9,0)=clhs10*clhs295 + clhs132*clhs296 + clhs133*clhs297 + clhs97;
lhs(9,1)=clhs135*clhs296 + clhs137*clhs297 + clhs162 + clhs295*clhs36 + clhs319;
lhs(9,2)=clhs138*clhs296 + clhs140*clhs297 + clhs205 + clhs295*clhs44;
lhs(9,3)=DN(0,1)*clhs302 + clhs171*clhs49 - clhs172;
lhs(9,4)=clhs143*clhs296 + clhs144*clhs297 + clhs248 + clhs295*clhs54;
lhs(9,5)=clhs148*clhs296 + clhs150*clhs297 + clhs266 + clhs295*clhs67 + clhs320;
lhs(9,6)=clhs152*clhs296 + clhs154*clhs297 + clhs281 + clhs295*clhs74;
lhs(9,7)=DN(1,1)*clhs302 + clhs269*clhs49 - clhs270;
lhs(9,8)=clhs158*clhs296 + clhs159*clhs297 + clhs295*clhs82 + clhs308;
lhs(9,9)=clhs163*clhs296 + clhs165*clhs297 + clhs295*clhs93 + clhs306 + clhs321*clhs5;
lhs(9,10)=clhs100*clhs295 + clhs167*clhs296 + clhs169*clhs297 + clhs323;
lhs(9,11)=DN(2,1)*clhs310;
lhs(9,12)=clhs108*clhs295 + clhs173*clhs296 + clhs174*clhs297 + clhs324;
lhs(9,13)=clhs119*clhs295 + clhs178*clhs296 + clhs180*clhs297 + clhs326 + clhs327;
lhs(9,14)=clhs126*clhs295 + clhs182*clhs296 + clhs184*clhs297 + clhs328;
lhs(9,15)=DN(3,1)*clhs302 - clhs329 + clhs330*clhs49;
lhs(10,0)=clhs103 + clhs13*clhs295 + clhs133*clhs296 + clhs188*clhs297;
lhs(10,1)=clhs137*clhs296 + clhs170 + clhs189*clhs297 + clhs295*clhs39;
lhs(10,2)=clhs140*clhs296 + clhs191*clhs297 + clhs207 + clhs295*clhs46 + clhs319;
lhs(10,3)=DN(0,2)*clhs302 + clhs209*clhs49 - clhs210;
lhs(10,4)=clhs144*clhs296 + clhs192*clhs297 + clhs249 + clhs295*clhs56;
lhs(10,5)=clhs150*clhs296 + clhs195*clhs297 + clhs268 + clhs295*clhs70;
lhs(10,6)=clhs154*clhs296 + clhs199*clhs297 + clhs283 + clhs295*clhs76 + clhs320;
lhs(10,7)=DN(1,2)*clhs302 + clhs284*clhs49 - clhs285;
lhs(10,8)=clhs159*clhs296 + clhs202*clhs297 + clhs295*clhs84 + clhs309;
lhs(10,9)=clhs165*clhs296 + clhs204*clhs297 + clhs295*clhs96 + clhs323;
lhs(10,10)=clhs102*clhs295 + clhs169*clhs296 + clhs208*clhs297 + clhs306 + clhs331*clhs5;
lhs(10,11)=DN(2,2)*clhs310;
lhs(10,12)=clhs110*clhs295 + clhs174*clhs296 + clhs211*clhs297 + clhs333;
lhs(10,13)=clhs122*clhs295 + clhs180*clhs296 + clhs213*clhs297 + clhs334;
lhs(10,14)=clhs128*clhs295 + clhs184*clhs296 + clhs217*clhs297 + clhs327 + clhs336;
lhs(10,15)=DN(3,2)*clhs302 - clhs337 + clhs338*clhs49;
lhs(11,0)=clhs7*(DN(2,0)*clhs220 + clhs104);
lhs(11,1)=clhs7*(DN(2,1)*clhs220 + clhs171);
lhs(11,2)=clhs7*(DN(2,2)*clhs220 + clhs209);
lhs(11,3)=clhs225;
lhs(11,4)=clhs7*(DN(2,0)*clhs222 + clhs250);
lhs(11,5)=clhs7*(DN(2,1)*clhs222 + clhs269);
lhs(11,6)=clhs7*(DN(2,2)*clhs222 + clhs284);
lhs(11,7)=clhs293;
lhs(11,8)=clhs295*clhs339;
lhs(11,9)=clhs296*clhs339;
lhs(11,10)=clhs297*clhs339;
lhs(11,11)=clhs24*(clhs305 + clhs321 + clhs331);
lhs(11,12)=clhs7*(DN(2,0)*clhs226 + clhs318);
lhs(11,13)=clhs7*(DN(2,1)*clhs226 + clhs330);
lhs(11,14)=clhs7*(DN(2,2)*clhs226 + clhs338);
lhs(11,15)=clhs340;
lhs(12,0)=clhs10*clhs342 + clhs113 + clhs13*clhs343 + clhs341*clhs6 + clhs347;
lhs(12,1)=clhs175 + clhs34*clhs341 + clhs342*clhs36 + clhs343*clhs39;
lhs(12,2)=clhs212 + clhs341*clhs42 + clhs342*clhs44 + clhs343*clhs46;
lhs(12,3)=DN(0,0)*clhs348 + clhs130*clhs49 - clhs131;
lhs(12,4)=clhs254 + clhs341*clhs52 + clhs342*clhs54 + clhs343*clhs56 + clhs350;
lhs(12,5)=clhs271 + clhs341*clhs65 + clhs342*clhs67 + clhs343*clhs70;
lhs(12,6)=clhs286 + clhs341*clhs72 + clhs342*clhs74 + clhs343*clhs76;
lhs(12,7)=DN(1,0)*clhs348 + clhs258*clhs49 - clhs259;
lhs(12,8)=clhs313 + clhs341*clhs80 + clhs342*clhs82 + clhs343*clhs84 + clhs351;
lhs(12,9)=clhs324 + clhs341*clhs91 + clhs342*clhs93 + clhs343*clhs96;
lhs(12,10)=clhs100*clhs342 + clhs102*clhs343 + clhs333 + clhs341*clhs98;
lhs(12,11)=DN(2,0)*clhs348 + clhs317*clhs49 - clhs318;
lhs(12,12)=clhs106*clhs341 + clhs108*clhs342 + clhs110*clhs343 + clhs352*clhs5 + clhs353;
lhs(12,13)=clhs117*clhs341 + clhs119*clhs342 + clhs122*clhs343 + clhs355;
lhs(12,14)=clhs124*clhs341 + clhs126*clhs342 + clhs128*clhs343 + clhs356;
lhs(12,15)=DN(3,0)*clhs357;
lhs(13,0)=clhs10*clhs341 + clhs123 + clhs132*clhs342 + clhs133*clhs343;
lhs(13,1)=clhs135*clhs342 + clhs137*clhs343 + clhs177 + clhs341*clhs36 + clhs358;
lhs(13,2)=clhs138*clhs342 + clhs140*clhs343 + clhs214 + clhs341*clhs44;
lhs(13,3)=DN(0,1)*clhs348 + clhs186*clhs49 - clhs187;
lhs(13,4)=clhs143*clhs342 + clhs144*clhs343 + clhs256 + clhs341*clhs54;
lhs(13,5)=clhs148*clhs342 + clhs150*clhs343 + clhs273 + clhs341*clhs67 + clhs359;
lhs(13,6)=clhs152*clhs342 + clhs154*clhs343 + clhs287 + clhs341*clhs74;
lhs(13,7)=DN(1,1)*clhs348 + clhs276*clhs49 - clhs277;
lhs(13,8)=clhs158*clhs342 + clhs159*clhs343 + clhs315 + clhs341*clhs82;
lhs(13,9)=clhs163*clhs342 + clhs165*clhs343 + clhs326 + clhs341*clhs93 + clhs360;
lhs(13,10)=clhs100*clhs341 + clhs167*clhs342 + clhs169*clhs343 + clhs334;
lhs(13,11)=DN(2,1)*clhs348 + clhs329*clhs49 - clhs330;
lhs(13,12)=clhs108*clhs341 + clhs173*clhs342 + clhs174*clhs343 + clhs355;
lhs(13,13)=clhs119*clhs341 + clhs178*clhs342 + clhs180*clhs343 + clhs353 + clhs361*clhs5;
lhs(13,14)=clhs126*clhs341 + clhs182*clhs342 + clhs184*clhs343 + clhs362;
lhs(13,15)=DN(3,1)*clhs357;
lhs(14,0)=clhs129 + clhs13*clhs341 + clhs133*clhs342 + clhs188*clhs343;
lhs(14,1)=clhs137*clhs342 + clhs185 + clhs189*clhs343 + clhs341*clhs39;
lhs(14,2)=clhs140*clhs342 + clhs191*clhs343 + clhs216 + clhs341*clhs46 + clhs358;
lhs(14,3)=DN(0,2)*clhs348 + clhs218*clhs49 - clhs219;
lhs(14,4)=clhs144*clhs342 + clhs192*clhs343 + clhs257 + clhs341*clhs56;
lhs(14,5)=clhs150*clhs342 + clhs195*clhs343 + clhs275 + clhs341*clhs70;
lhs(14,6)=clhs154*clhs342 + clhs199*clhs343 + clhs289 + clhs341*clhs76 + clhs359;
lhs(14,7)=DN(1,2)*clhs348 + clhs290*clhs49 - clhs291;
lhs(14,8)=clhs159*clhs342 + clhs202*clhs343 + clhs316 + clhs341*clhs84;
lhs(14,9)=clhs165*clhs342 + clhs204*clhs343 + clhs328 + clhs341*clhs96;
lhs(14,10)=clhs102*clhs341 + clhs169*clhs342 + clhs208*clhs343 + clhs336 + clhs360;
lhs(14,11)=DN(2,2)*clhs348 + clhs337*clhs49 - clhs338;
lhs(14,12)=clhs110*clhs341 + clhs174*clhs342 + clhs211*clhs343 + clhs356;
lhs(14,13)=clhs122*clhs341 + clhs180*clhs342 + clhs213*clhs343 + clhs362;
lhs(14,14)=clhs128*clhs341 + clhs184*clhs342 + clhs217*clhs343 + clhs353 + clhs363*clhs5;
lhs(14,15)=DN(3,2)*clhs357;
lhs(15,0)=clhs7*(DN(3,0)*clhs220 + clhs130);
lhs(15,1)=clhs7*(DN(3,1)*clhs220 + clhs186);
lhs(15,2)=clhs7*(DN(3,2)*clhs220 + clhs218);
lhs(15,3)=clhs227;
lhs(15,4)=clhs7*(DN(3,0)*clhs222 + clhs258);
lhs(15,5)=clhs7*(DN(3,1)*clhs222 + clhs276);
lhs(15,6)=clhs7*(DN(3,2)*clhs222 + clhs290);
lhs(15,7)=clhs294;
lhs(15,8)=clhs7*(DN(3,0)*clhs224 + clhs317);
lhs(15,9)=clhs7*(DN(3,1)*clhs224 + clhs329);
lhs(15,10)=clhs7*(DN(3,2)*clhs224 + clhs337);
lhs(15,11)=clhs340;
lhs(15,12)=clhs341*clhs364;
lhs(15,13)=clhs342*clhs364;
lhs(15,14)=clhs343*clhs364;
lhs(15,15)=clhs24*(clhs352 + clhs361 + clhs363);


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
    const double max_spectral_radius = rData.MaxSpectralRadius;
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
    const double volume_error_ratio = rData.VolumeErrorRate;

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = 1.0/(max_spectral_radius + 1);
const double crhs2 = rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)));
const double crhs3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs5 = rho*stab_c2*sqrt(pow(crhs3, 2) + pow(crhs4, 2));
const double crhs6 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs7 = (crhs6 - volume_error_ratio)*(crhs5*h/stab_c1 + mu);
const double crhs8 = v(0,0) - vn(0,0);
const double crhs9 = crhs1*crhs8 + vn(0,0);
const double crhs10 = v(1,0) - vn(1,0);
const double crhs11 = crhs1*crhs10 + vn(1,0);
const double crhs12 = v(2,0) - vn(2,0);
const double crhs13 = crhs1*crhs12 + vn(2,0);
const double crhs14 = rho*(crhs3*(DN(0,0)*crhs9 + DN(1,0)*crhs11 + DN(2,0)*crhs13) + crhs4*(DN(0,1)*crhs9 + DN(1,1)*crhs11 + DN(2,1)*crhs13));
const double crhs15 = -acceleration_alpha_method(0,0);
const double crhs16 = 1.0/dt;
const double crhs17 = 0.5*max_spectral_radius;
const double crhs18 = -crhs1;
const double crhs19 = crhs18 + 0.5;
const double crhs20 = 1.0/(-crhs1*(crhs17 - 1.5) + crhs19);
const double crhs21 = crhs16*crhs20;
const double crhs22 = crhs1*(1.5 - crhs17);
const double crhs23 = crhs20*(crhs1 - crhs22 + 0.5);
const double crhs24 = 0.5*crhs1;
const double crhs25 = crhs24*(max_spectral_radius - 3);
const double crhs26 = -acceleration_alpha_method(1,0);
const double crhs27 = -acceleration_alpha_method(2,0);
const double crhs28 = N[0]*(acceleration_alpha_method(0,0) - crhs25*(-acceleration_alpha_method(0,0)*crhs23 + crhs15 + crhs21*crhs8)) + N[1]*(acceleration_alpha_method(1,0) - crhs25*(-acceleration_alpha_method(1,0)*crhs23 + crhs10*crhs21 + crhs26)) + N[2]*(acceleration_alpha_method(2,0) - crhs25*(-acceleration_alpha_method(2,0)*crhs23 + crhs12*crhs21 + crhs27));
const double crhs29 = N[0]*rho;
const double crhs30 = 1.0/(crhs19 + crhs22);
const double crhs31 = crhs16*crhs30;
const double crhs32 = crhs30*(crhs18 + crhs22 - 0.5);
const double crhs33 = crhs24*(3 - max_spectral_radius);
const double crhs34 = 1.0/(crhs16*dyn_tau*rho + crhs5/h + mu*stab_c1/pow(h, 2));
const double crhs35 = crhs34*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs14 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs33*(acceleration_alpha_method(0,0)*crhs32 + crhs15 + crhs31*crhs8)) + N[1]*(acceleration_alpha_method(1,0) + crhs33*(acceleration_alpha_method(1,0)*crhs32 + crhs10*crhs31 + crhs26)) + N[2]*(acceleration_alpha_method(2,0) + crhs33*(acceleration_alpha_method(2,0)*crhs32 + crhs12*crhs31 + crhs27))));
const double crhs36 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs37 = crhs29*crhs36;
const double crhs38 = rho*(DN(0,0)*crhs3 + DN(0,1)*crhs4);
const double crhs39 = rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)));
const double crhs40 = v(0,1) - vn(0,1);
const double crhs41 = crhs1*crhs40 + vn(0,1);
const double crhs42 = v(1,1) - vn(1,1);
const double crhs43 = crhs1*crhs42 + vn(1,1);
const double crhs44 = v(2,1) - vn(2,1);
const double crhs45 = crhs1*crhs44 + vn(2,1);
const double crhs46 = rho*(crhs3*(DN(0,0)*crhs41 + DN(1,0)*crhs43 + DN(2,0)*crhs45) + crhs4*(DN(0,1)*crhs41 + DN(1,1)*crhs43 + DN(2,1)*crhs45));
const double crhs47 = -acceleration_alpha_method(0,1);
const double crhs48 = -acceleration_alpha_method(1,1);
const double crhs49 = -acceleration_alpha_method(2,1);
const double crhs50 = N[0]*(acceleration_alpha_method(0,1) - crhs25*(-acceleration_alpha_method(0,1)*crhs23 + crhs21*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) - crhs25*(-acceleration_alpha_method(1,1)*crhs23 + crhs21*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) - crhs25*(-acceleration_alpha_method(2,1)*crhs23 + crhs21*crhs44 + crhs49));
const double crhs51 = crhs34*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs39 + crhs46 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs33*(acceleration_alpha_method(0,1)*crhs32 + crhs31*crhs40 + crhs47)) + N[1]*(acceleration_alpha_method(1,1) + crhs33*(acceleration_alpha_method(1,1)*crhs32 + crhs31*crhs42 + crhs48)) + N[2]*(acceleration_alpha_method(2,1) + crhs33*(acceleration_alpha_method(2,1)*crhs32 + crhs31*crhs44 + crhs49))));
const double crhs52 = -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + crhs6) + volume_error_ratio;
const double crhs53 = N[1]*rho;
const double crhs54 = crhs36*crhs53;
const double crhs55 = rho*(DN(1,0)*crhs3 + DN(1,1)*crhs4);
const double crhs56 = N[2]*rho;
const double crhs57 = crhs36*crhs56;
const double crhs58 = rho*(DN(2,0)*crhs3 + DN(2,1)*crhs4);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs7 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs14 + N[0]*crhs2 - crhs28*crhs29 - crhs35*crhs37 - crhs35*crhs38;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs7 - DN(0,1)*stress[1] + N[0]*crhs39 - N[0]*crhs46 - crhs29*crhs50 - crhs37*crhs51 - crhs38*crhs51;
rhs[2]=-DN(0,0)*crhs35 - DN(0,1)*crhs51 + N[0]*crhs52;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs7 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs14 + N[1]*crhs2 - crhs28*crhs53 - crhs35*crhs54 - crhs35*crhs55;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs7 - DN(1,1)*stress[1] + N[1]*crhs39 - N[1]*crhs46 - crhs50*crhs53 - crhs51*crhs54 - crhs51*crhs55;
rhs[5]=-DN(1,0)*crhs35 - DN(1,1)*crhs51 + N[1]*crhs52;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs7 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs14 + N[2]*crhs2 - crhs28*crhs56 - crhs35*crhs57 - crhs35*crhs58;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs7 - DN(2,1)*stress[1] + N[2]*crhs39 - N[2]*crhs46 - crhs50*crhs56 - crhs51*crhs57 - crhs51*crhs58;
rhs[8]=-DN(2,0)*crhs35 - DN(2,1)*crhs51 + N[2]*crhs52;


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
    const double max_spectral_radius = rData.MaxSpectralRadius;
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
    const double volume_error_ratio = rData.VolumeErrorRate;

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = 1.0/(max_spectral_radius + 1);
const double crhs2 = rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs1*(f(3,0) - fn(3,0)) + fn(3,0)));
const double crhs3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs5 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs6 = rho*stab_c2*sqrt(pow(crhs3, 2) + pow(crhs4, 2) + pow(crhs5, 2));
const double crhs7 = DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs8 = (crhs7 - volume_error_ratio)*(crhs6*h/stab_c1 + mu);
const double crhs9 = v(0,0) - vn(0,0);
const double crhs10 = crhs1*crhs9 + vn(0,0);
const double crhs11 = v(1,0) - vn(1,0);
const double crhs12 = crhs1*crhs11 + vn(1,0);
const double crhs13 = v(2,0) - vn(2,0);
const double crhs14 = crhs1*crhs13 + vn(2,0);
const double crhs15 = v(3,0) - vn(3,0);
const double crhs16 = crhs1*crhs15 + vn(3,0);
const double crhs17 = rho*(crhs3*(DN(0,0)*crhs10 + DN(1,0)*crhs12 + DN(2,0)*crhs14 + DN(3,0)*crhs16) + crhs4*(DN(0,1)*crhs10 + DN(1,1)*crhs12 + DN(2,1)*crhs14 + DN(3,1)*crhs16) + crhs5*(DN(0,2)*crhs10 + DN(1,2)*crhs12 + DN(2,2)*crhs14 + DN(3,2)*crhs16));
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
const double crhs31 = -acceleration_alpha_method(3,0);
const double crhs32 = N[0]*(acceleration_alpha_method(0,0) - crhs28*(-acceleration_alpha_method(0,0)*crhs26 + crhs18 + crhs24*crhs9)) + N[1]*(acceleration_alpha_method(1,0) - crhs28*(-acceleration_alpha_method(1,0)*crhs26 + crhs11*crhs24 + crhs29)) + N[2]*(acceleration_alpha_method(2,0) - crhs28*(-acceleration_alpha_method(2,0)*crhs26 + crhs13*crhs24 + crhs30)) + N[3]*(acceleration_alpha_method(3,0) - crhs28*(-acceleration_alpha_method(3,0)*crhs26 + crhs15*crhs24 + crhs31));
const double crhs33 = N[0]*rho;
const double crhs34 = 1.0/(crhs22 + crhs25);
const double crhs35 = crhs19*crhs34;
const double crhs36 = crhs34*(crhs21 + crhs25 - 0.5);
const double crhs37 = crhs27*(3 - max_spectral_radius);
const double crhs38 = 1.0/(crhs19*dyn_tau*rho + crhs6/h + mu*stab_c1/pow(h, 2));
const double crhs39 = crhs38*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs17 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs37*(acceleration_alpha_method(0,0)*crhs36 + crhs18 + crhs35*crhs9)) + N[1]*(acceleration_alpha_method(1,0) + crhs37*(acceleration_alpha_method(1,0)*crhs36 + crhs11*crhs35 + crhs29)) + N[2]*(acceleration_alpha_method(2,0) + crhs37*(acceleration_alpha_method(2,0)*crhs36 + crhs13*crhs35 + crhs30)) + N[3]*(acceleration_alpha_method(3,0) + crhs37*(acceleration_alpha_method(3,0)*crhs36 + crhs15*crhs35 + crhs31))));
const double crhs40 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs41 = crhs33*crhs40;
const double crhs42 = rho*(DN(0,0)*crhs3 + DN(0,1)*crhs4 + DN(0,2)*crhs5);
const double crhs43 = rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs1*(f(3,1) - fn(3,1)) + fn(3,1)));
const double crhs44 = v(0,1) - vn(0,1);
const double crhs45 = crhs1*crhs44 + vn(0,1);
const double crhs46 = v(1,1) - vn(1,1);
const double crhs47 = crhs1*crhs46 + vn(1,1);
const double crhs48 = v(2,1) - vn(2,1);
const double crhs49 = crhs1*crhs48 + vn(2,1);
const double crhs50 = v(3,1) - vn(3,1);
const double crhs51 = crhs1*crhs50 + vn(3,1);
const double crhs52 = rho*(crhs3*(DN(0,0)*crhs45 + DN(1,0)*crhs47 + DN(2,0)*crhs49 + DN(3,0)*crhs51) + crhs4*(DN(0,1)*crhs45 + DN(1,1)*crhs47 + DN(2,1)*crhs49 + DN(3,1)*crhs51) + crhs5*(DN(0,2)*crhs45 + DN(1,2)*crhs47 + DN(2,2)*crhs49 + DN(3,2)*crhs51));
const double crhs53 = -acceleration_alpha_method(0,1);
const double crhs54 = -acceleration_alpha_method(1,1);
const double crhs55 = -acceleration_alpha_method(2,1);
const double crhs56 = -acceleration_alpha_method(3,1);
const double crhs57 = N[0]*(acceleration_alpha_method(0,1) - crhs28*(-acceleration_alpha_method(0,1)*crhs26 + crhs24*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) - crhs28*(-acceleration_alpha_method(1,1)*crhs26 + crhs24*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) - crhs28*(-acceleration_alpha_method(2,1)*crhs26 + crhs24*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) - crhs28*(-acceleration_alpha_method(3,1)*crhs26 + crhs24*crhs50 + crhs56));
const double crhs58 = crhs38*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs43 + crhs52 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs37*(acceleration_alpha_method(0,1)*crhs36 + crhs35*crhs44 + crhs53)) + N[1]*(acceleration_alpha_method(1,1) + crhs37*(acceleration_alpha_method(1,1)*crhs36 + crhs35*crhs46 + crhs54)) + N[2]*(acceleration_alpha_method(2,1) + crhs37*(acceleration_alpha_method(2,1)*crhs36 + crhs35*crhs48 + crhs55)) + N[3]*(acceleration_alpha_method(3,1) + crhs37*(acceleration_alpha_method(3,1)*crhs36 + crhs35*crhs50 + crhs56))));
const double crhs59 = rho*(N[0]*(crhs1*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs1*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs1*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs1*(f(3,2) - fn(3,2)) + fn(3,2)));
const double crhs60 = v(0,2) - vn(0,2);
const double crhs61 = crhs1*crhs60 + vn(0,2);
const double crhs62 = v(1,2) - vn(1,2);
const double crhs63 = crhs1*crhs62 + vn(1,2);
const double crhs64 = v(2,2) - vn(2,2);
const double crhs65 = crhs1*crhs64 + vn(2,2);
const double crhs66 = v(3,2) - vn(3,2);
const double crhs67 = crhs1*crhs66 + vn(3,2);
const double crhs68 = rho*(crhs3*(DN(0,0)*crhs61 + DN(1,0)*crhs63 + DN(2,0)*crhs65 + DN(3,0)*crhs67) + crhs4*(DN(0,1)*crhs61 + DN(1,1)*crhs63 + DN(2,1)*crhs65 + DN(3,1)*crhs67) + crhs5*(DN(0,2)*crhs61 + DN(1,2)*crhs63 + DN(2,2)*crhs65 + DN(3,2)*crhs67));
const double crhs69 = -acceleration_alpha_method(0,2);
const double crhs70 = -acceleration_alpha_method(1,2);
const double crhs71 = -acceleration_alpha_method(2,2);
const double crhs72 = -acceleration_alpha_method(3,2);
const double crhs73 = N[0]*(acceleration_alpha_method(0,2) - crhs28*(-acceleration_alpha_method(0,2)*crhs26 + crhs24*crhs60 + crhs69)) + N[1]*(acceleration_alpha_method(1,2) - crhs28*(-acceleration_alpha_method(1,2)*crhs26 + crhs24*crhs62 + crhs70)) + N[2]*(acceleration_alpha_method(2,2) - crhs28*(-acceleration_alpha_method(2,2)*crhs26 + crhs24*crhs64 + crhs71)) + N[3]*(acceleration_alpha_method(3,2) - crhs28*(-acceleration_alpha_method(3,2)*crhs26 + crhs24*crhs66 + crhs72));
const double crhs74 = crhs38*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs59 + crhs68 + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs37*(acceleration_alpha_method(0,2)*crhs36 + crhs35*crhs60 + crhs69)) + N[1]*(acceleration_alpha_method(1,2) + crhs37*(acceleration_alpha_method(1,2)*crhs36 + crhs35*crhs62 + crhs70)) + N[2]*(acceleration_alpha_method(2,2) + crhs37*(acceleration_alpha_method(2,2)*crhs36 + crhs35*crhs64 + crhs71)) + N[3]*(acceleration_alpha_method(3,2) + crhs37*(acceleration_alpha_method(3,2)*crhs36 + crhs35*crhs66 + crhs72))));
const double crhs75 = -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + crhs7) + volume_error_ratio;
const double crhs76 = N[1]*rho;
const double crhs77 = crhs40*crhs76;
const double crhs78 = rho*(DN(1,0)*crhs3 + DN(1,1)*crhs4 + DN(1,2)*crhs5);
const double crhs79 = N[2]*rho;
const double crhs80 = crhs40*crhs79;
const double crhs81 = rho*(DN(2,0)*crhs3 + DN(2,1)*crhs4 + DN(2,2)*crhs5);
const double crhs82 = N[3]*rho;
const double crhs83 = crhs40*crhs82;
const double crhs84 = rho*(DN(3,0)*crhs3 + DN(3,1)*crhs4 + DN(3,2)*crhs5);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs17 + N[0]*crhs2 - crhs32*crhs33 - crhs39*crhs41 - crhs39*crhs42;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs8 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs43 - N[0]*crhs52 - crhs33*crhs57 - crhs41*crhs58 - crhs42*crhs58;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs8 - DN(0,2)*stress[2] + N[0]*crhs59 - N[0]*crhs68 - crhs33*crhs73 - crhs41*crhs74 - crhs42*crhs74;
rhs[3]=-DN(0,0)*crhs39 - DN(0,1)*crhs58 - DN(0,2)*crhs74 + N[0]*crhs75;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs17 + N[1]*crhs2 - crhs32*crhs76 - crhs39*crhs77 - crhs39*crhs78;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs8 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs43 - N[1]*crhs52 - crhs57*crhs76 - crhs58*crhs77 - crhs58*crhs78;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs8 - DN(1,2)*stress[2] + N[1]*crhs59 - N[1]*crhs68 - crhs73*crhs76 - crhs74*crhs77 - crhs74*crhs78;
rhs[7]=-DN(1,0)*crhs39 - DN(1,1)*crhs58 - DN(1,2)*crhs74 + N[1]*crhs75;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs17 + N[2]*crhs2 - crhs32*crhs79 - crhs39*crhs80 - crhs39*crhs81;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs8 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs43 - N[2]*crhs52 - crhs57*crhs79 - crhs58*crhs80 - crhs58*crhs81;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs8 - DN(2,2)*stress[2] + N[2]*crhs59 - N[2]*crhs68 - crhs73*crhs79 - crhs74*crhs80 - crhs74*crhs81;
rhs[11]=-DN(2,0)*crhs39 - DN(2,1)*crhs58 - DN(2,2)*crhs74 + N[2]*crhs75;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs8 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs17 + N[3]*crhs2 - crhs32*crhs82 - crhs39*crhs83 - crhs39*crhs84;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs8 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs43 - N[3]*crhs52 - crhs57*crhs82 - crhs58*crhs83 - crhs58*crhs84;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs8 - DN(3,2)*stress[2] + N[3]*crhs59 - N[3]*crhs68 - crhs73*crhs82 - crhs74*crhs83 - crhs74*crhs84;
rhs[15]=-DN(3,0)*crhs39 - DN(3,1)*crhs58 - DN(3,2)*crhs74 + N[3]*crhs75;


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
    const double max_spectral_radius = rData.MaxSpectralRadius;
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
    const double volume_error_ratio = rData.VolumeErrorRate;

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
const double cH4 = 0.5*cH3*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH5 = 1.0/(cH3*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2));
const double cH6 = cH5*rho;
const double cH7 = cH6*(DN(0,0)*cH1 + DN(0,1)*cH2 + N[0]*cH4);
const double cH8 = cH6*(DN(1,0)*cH1 + DN(1,1)*cH2 + N[1]*cH4);
const double cH9 = cH6*(DN(2,0)*cH1 + DN(2,1)*cH2 + N[2]*cH4);
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
const double crhs_ee2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee3 = v(0,0) - vn(0,0);
const double crhs_ee4 = crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 = v(1,0) - vn(1,0);
const double crhs_ee6 = crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 = v(2,0) - vn(2,0);
const double crhs_ee8 = crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee10 = 1.0/dt;
const double crhs_ee11 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee12 = 1.0/(crhs_ee11 + 0.5);
const double crhs_ee13 = crhs_ee10*crhs_ee12;
const double crhs_ee14 = crhs_ee12*(crhs_ee11 - 0.5);
const double crhs_ee15 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee16 = 1.0/(crhs_ee10*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee2, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee17 = crhs_ee16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee15*(acceleration_alpha_method(0,0)*crhs_ee14 - acceleration_alpha_method(0,0) + crhs_ee13*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee15*(acceleration_alpha_method(1,0)*crhs_ee14 - acceleration_alpha_method(1,0) + crhs_ee13*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee15*(acceleration_alpha_method(2,0)*crhs_ee14 - acceleration_alpha_method(2,0) + crhs_ee13*crhs_ee7)) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8) + crhs_ee9*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8)));
const double crhs_ee18 = v(0,1) - vn(0,1);
const double crhs_ee19 = crhs_ee0*crhs_ee18 + vn(0,1);
const double crhs_ee20 = v(1,1) - vn(1,1);
const double crhs_ee21 = crhs_ee0*crhs_ee20 + vn(1,1);
const double crhs_ee22 = v(2,1) - vn(2,1);
const double crhs_ee23 = crhs_ee0*crhs_ee22 + vn(2,1);
const double crhs_ee24 = crhs_ee16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee15*(acceleration_alpha_method(0,1)*crhs_ee14 - acceleration_alpha_method(0,1) + crhs_ee13*crhs_ee18)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee15*(acceleration_alpha_method(1,1)*crhs_ee14 - acceleration_alpha_method(1,1) + crhs_ee13*crhs_ee20)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee15*(acceleration_alpha_method(2,1)*crhs_ee14 - acceleration_alpha_method(2,1) + crhs_ee13*crhs_ee22)) + crhs_ee2*(DN(0,0)*crhs_ee19 + DN(1,0)*crhs_ee21 + DN(2,0)*crhs_ee23) + crhs_ee9*(DN(0,1)*crhs_ee19 + DN(1,1)*crhs_ee21 + DN(2,1)*crhs_ee23)));
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
    const double max_spectral_radius = rData.MaxSpectralRadius;
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
    const double volume_error_ratio = rData.VolumeErrorRate;

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
const double cH5 = 0.5*cH4*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH6 = 1.0/(cH4*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 = cH6*rho;
const double cH8 = cH7*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 + N[0]*cH5);
const double cH9 = cH7*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 + N[1]*cH5);
const double cH10 = cH7*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 + N[2]*cH5);
const double cH11 = cH7*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 + N[3]*cH5);
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
const double crhs_ee13 = 1.0/dt;
const double crhs_ee14 = crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee15 = 1.0/(crhs_ee14 + 0.5);
const double crhs_ee16 = crhs_ee13*crhs_ee15;
const double crhs_ee17 = crhs_ee15*(crhs_ee14 - 0.5);
const double crhs_ee18 = 0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee19 = 1.0/(crhs_ee13*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee11, 2) + pow(crhs_ee12, 2) + pow(crhs_ee2, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee20 = crhs_ee19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs_ee0*(f(3,0) - fn(3,0)) + fn(3,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee18*(acceleration_alpha_method(0,0)*crhs_ee17 - acceleration_alpha_method(0,0) + crhs_ee16*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee18*(acceleration_alpha_method(1,0)*crhs_ee17 - acceleration_alpha_method(1,0) + crhs_ee16*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee18*(acceleration_alpha_method(2,0)*crhs_ee17 - acceleration_alpha_method(2,0) + crhs_ee16*crhs_ee7)) + N[3]*(acceleration_alpha_method(3,0) + crhs_ee18*(acceleration_alpha_method(3,0)*crhs_ee17 - acceleration_alpha_method(3,0) + crhs_ee16*crhs_ee9)) + crhs_ee11*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8 + DN(3,1)*crhs_ee10) + crhs_ee12*(DN(0,2)*crhs_ee4 + DN(1,2)*crhs_ee6 + DN(2,2)*crhs_ee8 + DN(3,2)*crhs_ee10) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8 + DN(3,0)*crhs_ee10)));
const double crhs_ee21 = v(0,1) - vn(0,1);
const double crhs_ee22 = crhs_ee0*crhs_ee21 + vn(0,1);
const double crhs_ee23 = v(1,1) - vn(1,1);
const double crhs_ee24 = crhs_ee0*crhs_ee23 + vn(1,1);
const double crhs_ee25 = v(2,1) - vn(2,1);
const double crhs_ee26 = crhs_ee0*crhs_ee25 + vn(2,1);
const double crhs_ee27 = v(3,1) - vn(3,1);
const double crhs_ee28 = crhs_ee0*crhs_ee27 + vn(3,1);
const double crhs_ee29 = crhs_ee19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs_ee0*(f(3,1) - fn(3,1)) + fn(3,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee18*(acceleration_alpha_method(0,1)*crhs_ee17 - acceleration_alpha_method(0,1) + crhs_ee16*crhs_ee21)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee18*(acceleration_alpha_method(1,1)*crhs_ee17 - acceleration_alpha_method(1,1) + crhs_ee16*crhs_ee23)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee18*(acceleration_alpha_method(2,1)*crhs_ee17 - acceleration_alpha_method(2,1) + crhs_ee16*crhs_ee25)) + N[3]*(acceleration_alpha_method(3,1) + crhs_ee18*(acceleration_alpha_method(3,1)*crhs_ee17 - acceleration_alpha_method(3,1) + crhs_ee16*crhs_ee27)) + crhs_ee11*(DN(0,1)*crhs_ee22 + DN(1,1)*crhs_ee24 + DN(2,1)*crhs_ee26 + DN(3,1)*crhs_ee28) + crhs_ee12*(DN(0,2)*crhs_ee22 + DN(1,2)*crhs_ee24 + DN(2,2)*crhs_ee26 + DN(3,2)*crhs_ee28) + crhs_ee2*(DN(0,0)*crhs_ee22 + DN(1,0)*crhs_ee24 + DN(2,0)*crhs_ee26 + DN(3,0)*crhs_ee28)));
const double crhs_ee30 = v(0,2) - vn(0,2);
const double crhs_ee31 = crhs_ee0*crhs_ee30 + vn(0,2);
const double crhs_ee32 = v(1,2) - vn(1,2);
const double crhs_ee33 = crhs_ee0*crhs_ee32 + vn(1,2);
const double crhs_ee34 = v(2,2) - vn(2,2);
const double crhs_ee35 = crhs_ee0*crhs_ee34 + vn(2,2);
const double crhs_ee36 = v(3,2) - vn(3,2);
const double crhs_ee37 = crhs_ee0*crhs_ee36 + vn(3,2);
const double crhs_ee38 = crhs_ee19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs_ee0*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs_ee0*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs_ee0*(f(3,2) - fn(3,2)) + fn(3,2))) + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs_ee18*(acceleration_alpha_method(0,2)*crhs_ee17 - acceleration_alpha_method(0,2) + crhs_ee16*crhs_ee30)) + N[1]*(acceleration_alpha_method(1,2) + crhs_ee18*(acceleration_alpha_method(1,2)*crhs_ee17 - acceleration_alpha_method(1,2) + crhs_ee16*crhs_ee32)) + N[2]*(acceleration_alpha_method(2,2) + crhs_ee18*(acceleration_alpha_method(2,2)*crhs_ee17 - acceleration_alpha_method(2,2) + crhs_ee16*crhs_ee34)) + N[3]*(acceleration_alpha_method(3,2) + crhs_ee18*(acceleration_alpha_method(3,2)*crhs_ee17 - acceleration_alpha_method(3,2) + crhs_ee16*crhs_ee36)) + crhs_ee11*(DN(0,1)*crhs_ee31 + DN(1,1)*crhs_ee33 + DN(2,1)*crhs_ee35 + DN(3,1)*crhs_ee37) + crhs_ee12*(DN(0,2)*crhs_ee31 + DN(1,2)*crhs_ee33 + DN(2,2)*crhs_ee35 + DN(3,2)*crhs_ee37) + crhs_ee2*(DN(0,0)*crhs_ee31 + DN(1,0)*crhs_ee33 + DN(2,0)*crhs_ee35 + DN(3,0)*crhs_ee37)));
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
    const double alpha_m = 0.5*((3-rData.MaxSpectralRadius)/(1+rData.MaxSpectralRadius));

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
