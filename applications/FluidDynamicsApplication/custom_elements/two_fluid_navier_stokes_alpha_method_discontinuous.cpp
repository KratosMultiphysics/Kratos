//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

// Application includes
#include "two_fluid_navier_stokes_alpha_method_discontinuous.h"
#include "data_containers/two_fluid_navier_stokes/two_fluid_navier_stokes_alpha_method_discontinuous_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethodDiscontinuous>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethodDiscontinuous>(NewId, pGeom, pProperties);
}


template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize RHS vector
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Resize and intialize LHS matrix
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration){
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        if (data.IsCut()) {
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_ausas_pos, shape_functions_ausas_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_ausas_pos, shape_derivatives_ausas_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_ausas_pos,
                shape_functions_ausas_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_ausas_pos,
                shape_derivatives_ausas_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg);

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side
                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                for (unsigned int g = 0; g < gauss_weights.size(); ++g){
                    this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                }
            } else {
                MatrixType Vtot = ZeroMatrix(LocalSize, NumNodes);
                MatrixType Htot = ZeroMatrix(NumNodes, LocalSize);
                MatrixType Kee_tot = ZeroMatrix(NumNodes, NumNodes);
                VectorType rhs_ee_tot = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); ++g_pos){
                    this->UpdateIntegrationPointDataDiscontinuous(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_ausas_pos, g_pos),
                        shape_derivatives_ausas_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); ++g_neg){
                    this->UpdateIntegrationPointDataDiscontinuous(
                        data,
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_ausas_neg, g_neg),
                        shape_derivatives_ausas_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                // Without pressure gradient stabilization, volume ratio is checked during condensation
                // Also, without surface tension, zero pressure difference is penalized
                this->CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
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
                this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
    } else{
        KRATOS_ERROR << "TwoFluidNavierStokesAlphaMethodDiscontinuous is expected to manage time integration." << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <>
double TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>::CalculateArtificialDynamicViscositySpecialization(TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3> &rData) const
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
    const auto &N_vel = rData.N_vel;
    const auto &DN_vel = rData.DN_DX_vel;
    const double art_dyn_visc_coeff = 0.8;

    double grad_v_norm;
    grad_v_norm=sqrt(pow(fabs(DN_vel(0,0)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,0)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,1)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,1)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0))), 2));


    // Check that velocity gradient norm is non-zero
    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        artificial_mu=0.5*art_dyn_visc_coeff*h*sqrt(pow(fabs(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + rho*((DN_vel(0,0)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0)))*(N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0)) + (DN_vel(0,1)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0)))*(N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1))) + rho*(N_vel[0]*(acceleration_alpha_method(0,0) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(0,0) - acceleration_alpha_method(0,0)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(0,0) - vn(0,0))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0)) + N_vel[1]*(acceleration_alpha_method(1,0) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(1,0) - acceleration_alpha_method(1,0)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(1,0) - vn(1,0))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0)) + N_vel[2]*(acceleration_alpha_method(2,0) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(2,0) - acceleration_alpha_method(2,0)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(2,0) - vn(2,0))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0))) - rho*(N_vel[0]*(fn(0,0) + 1.0*(f(0,0) - fn(0,0))/(max_spectral_radius + 1.0)) + N_vel[1]*(fn(1,0) + 1.0*(f(1,0) - fn(1,0))/(max_spectral_radius + 1.0)) + N_vel[2]*(fn(2,0) + 1.0*(f(2,0) - fn(2,0))/(max_spectral_radius + 1.0)))), 2) + pow(fabs(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + rho*((DN_vel(0,0)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0)))*(N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0)) + (DN_vel(0,1)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0)))*(N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1))) + rho*(N_vel[0]*(acceleration_alpha_method(0,1) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(0,1) - acceleration_alpha_method(0,1)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(0,1) - vn(0,1))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0)) + N_vel[1]*(acceleration_alpha_method(1,1) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(1,1) - acceleration_alpha_method(1,1)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(1,1) - vn(1,1))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0)) + N_vel[2]*(acceleration_alpha_method(2,1) + 0.5*(max_spectral_radius - 3.0)*(acceleration_alpha_method(2,1) - acceleration_alpha_method(2,1)*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) + 0.5 + 1.0/(max_spectral_radius + 1.0))/(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0)) + (v(2,1) - vn(2,1))/(dt*(0.5*(max_spectral_radius - 3.0)/(max_spectral_radius + 1.0) - 0.5 + 1.0/(max_spectral_radius + 1.0))))/(max_spectral_radius + 1.0))) - rho*(N_vel[0]*(fn(0,1) + 1.0*(f(0,1) - fn(0,1))/(max_spectral_radius + 1.0)) + N_vel[1]*(fn(1,1) + 1.0*(f(1,1) - fn(1,1))/(max_spectral_radius + 1.0)) + N_vel[2]*(fn(2,1) + 1.0*(f(2,1) - fn(2,1))/(max_spectral_radius + 1.0)))), 2))/sqrt(pow(fabs(DN_vel(0,0)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,0)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,0)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,0)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,1)*(vn(0,0) + 1.0*(v(0,0) - vn(0,0))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,0) + 1.0*(v(1,0) - vn(1,0))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,0) + 1.0*(v(2,0) - vn(2,0))/(max_spectral_radius + 1.0))), 2) + pow(fabs(DN_vel(0,1)*(vn(0,1) + 1.0*(v(0,1) - vn(0,1))/(max_spectral_radius + 1.0)) + DN_vel(1,1)*(vn(1,1) + 1.0*(v(1,1) - vn(1,1))/(max_spectral_radius + 1.0)) + DN_vel(2,1)*(vn(2,1) + 1.0*(v(2,1) - vn(2,1))/(max_spectral_radius + 1.0))), 2));

    }

    return artificial_mu;
}

template <>
double TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::CalculateArtificialDynamicViscositySpecialization(TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4> &rData) const
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
    const auto &N_vel = rData.N_vel;
    const auto &DN_vel = rData.DN_DX_vel;
    const double art_dyn_visc_coeff = 0.8;

    double grad_v_norm;
    //substitute_artificial_mu_grad_v_norm_3D_4N

    // Check that velocity gradient norm is non-zero
    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        //substitute_artificial_mu_3D_4N
    }

    return artificial_mu;
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>& rData) const
{
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
    const BoundedMatrix<double,3,2> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& r_DN_DX_vel = rData.DN_DX_vel;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += r_DN_DX_vel(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += r_DN_DX_vel(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += r_DN_DX_vel(i,0)*velocity_alpha(i,1) + r_DN_DX_vel(i,1)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>& rData) const
{
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
    const BoundedMatrix<double,4,3> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& r_DN_DX_vel = rData.DN_DX_vel;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(6);
    for (unsigned int i = 0; i < 4; i++) {
        r_strain_rate[0] += r_DN_DX_vel(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += r_DN_DX_vel(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += r_DN_DX_vel(i,2)*velocity_alpha(i,2);
        r_strain_rate[3] += r_DN_DX_vel(i,0)*velocity_alpha(i,1) + r_DN_DX_vel(i,1)*velocity_alpha(i,0);
        r_strain_rate[4] += r_DN_DX_vel(i,1)*velocity_alpha(i,2) + r_DN_DX_vel(i,2)*velocity_alpha(i,1);
        r_strain_rate[5] += r_DN_DX_vel(i,0)*velocity_alpha(i,2) + r_DN_DX_vel(i,2)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Add current Gauss point LHS contribution
    const double gauss_weight = rData.Weight;
            const double crLHS0 = N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0);
        const double crLHS1 = N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1);
        const double crLHS2 = rho*stab_c2*sqrt(pow(crLHS0, 2) + pow(crLHS1, 2));
        const double crLHS3 = crLHS2*h/stab_c1 + mu;
        const double crLHS4 = C(0,0)*DN_vel(0,0) + C(0,2)*DN_vel(0,1);
        const double crLHS5 = 1.0/(max_spectral_radius + 1);
        const double crLHS6 = 1.0*crLHS5;
        const double crLHS7 = DN_vel(0,0)*crLHS6;
        const double crLHS8 = C(0,2)*DN_vel(0,0);
        const double crLHS9 = C(2,2)*DN_vel(0,1) + crLHS8;
        const double crLHS10 = DN_vel(0,1)*crLHS6;
        const double crLHS11 = DN_vel(0,0)*crLHS0;
        const double crLHS12 = DN_vel(0,1)*crLHS1;
        const double crLHS13 = crLHS11 + crLHS12;
        const double crLHS14 = crLHS13*rho;
        const double crLHS15 = 1.0/(max_spectral_radius + 1.0);
        const double crLHS16 = 1.0*crLHS15;
        const double crLHS17 = N_vel[0]*crLHS16;
        const double crLHS18 = 1.0/dt;
        const double crLHS19 = 0.5*crLHS18;
        const double crLHS20 = N_vel[0]*crLHS19;
        const double crLHS21 = 0.5*max_spectral_radius - 1.5;
        const double crLHS22 = crLHS15*crLHS21 + crLHS16 - 0.5;
        const double crLHS23 = max_spectral_radius - 3.0;
        const double crLHS24 = -crLHS23/crLHS22;
        const double crLHS25 = 1.0*crLHS11 + 1.0*crLHS12;
        const double crLHS26 = -crLHS20*crLHS24 + crLHS25;
        const double crLHS27 = crLHS18*rho;
        const double crLHS28 = 1.0/(crLHS2/h + crLHS27*dyn_tau + mu*stab_c1/pow(h, 2));
        const double crLHS29 = crLHS28*pow(rho, 2);
        const double crLHS30 = crLHS16*crLHS29;
        const double crLHS31 = crLHS13*crLHS30;
        const double crLHS32 = 0.5*crLHS15*crLHS23*crLHS27/crLHS22;
        const double crLHS33 = DN_vel(0,0)*vconv(0,0) + DN_vel(0,1)*vconv(0,1) + DN_vel(1,0)*vconv(1,0) + DN_vel(1,1)*vconv(1,1) + DN_vel(2,0)*vconv(2,0) + DN_vel(2,1)*vconv(2,1);
        const double crLHS34 = crLHS29*crLHS33;
        const double crLHS35 = crLHS17*crLHS34;
        const double crLHS36 = pow(N_vel[0], 2)*crLHS32 + crLHS14*crLHS17 + crLHS26*crLHS31 + crLHS26*crLHS35;
        const double crLHS37 = C(0,1)*DN_vel(0,1) + crLHS8;
        const double crLHS38 = C(1,2)*DN_vel(0,1);
        const double crLHS39 = C(2,2)*DN_vel(0,0) + crLHS38;
        const double crLHS40 = DN_vel(0,0)*crLHS3;
        const double crLHS41 = DN_vel(0,1)*crLHS40;
        const double crLHS42 = DN_vel(0,0)*N[0];
        const double crLHS43 = DN(0,0)*crLHS28;
        const double crLHS44 = 1.0*crLHS33*rho;
        const double crLHS45 = N_vel[0]*crLHS44;
        const double crLHS46 = 1.0*crLHS14;
        const double crLHS47 = C(0,0)*DN_vel(1,0) + C(0,2)*DN_vel(1,1);
        const double crLHS48 = C(0,2)*DN_vel(1,0);
        const double crLHS49 = C(2,2)*DN_vel(1,1) + crLHS48;
        const double crLHS50 = N_vel[0]*crLHS32;
        const double crLHS51 = N_vel[1]*crLHS50;
        const double crLHS52 = DN_vel(1,0)*crLHS40 + crLHS51;
        const double crLHS53 = DN_vel(1,0)*crLHS0;
        const double crLHS54 = DN_vel(1,1)*crLHS1;
        const double crLHS55 = crLHS53 + crLHS54;
        const double crLHS56 = crLHS17*rho;
        const double crLHS57 = crLHS19*crLHS24;
        const double crLHS58 = 1.0*crLHS53 + 1.0*crLHS54;
        const double crLHS59 = -N_vel[1]*crLHS57 + crLHS58;
        const double crLHS60 = crLHS31*crLHS59 + crLHS35*crLHS59 + crLHS55*crLHS56;
        const double crLHS61 = C(0,1)*DN_vel(1,1) + crLHS48;
        const double crLHS62 = C(1,2)*DN_vel(1,1);
        const double crLHS63 = C(2,2)*DN_vel(1,0) + crLHS62;
        const double crLHS64 = DN_vel(1,1)*crLHS40;
        const double crLHS65 = DN_vel(0,0)*N[1];
        const double crLHS66 = crLHS28*crLHS45;
        const double crLHS67 = crLHS28*crLHS46;
        const double crLHS68 = C(0,0)*DN_vel(2,0) + C(0,2)*DN_vel(2,1);
        const double crLHS69 = C(0,2)*DN_vel(2,0);
        const double crLHS70 = C(2,2)*DN_vel(2,1) + crLHS69;
        const double crLHS71 = N_vel[2]*crLHS50;
        const double crLHS72 = DN_vel(2,0)*crLHS40 + crLHS71;
        const double crLHS73 = DN_vel(2,0)*crLHS0;
        const double crLHS74 = DN_vel(2,1)*crLHS1;
        const double crLHS75 = crLHS73 + crLHS74;
        const double crLHS76 = 1.0*crLHS73 + 1.0*crLHS74;
        const double crLHS77 = -N_vel[2]*crLHS57 + crLHS76;
        const double crLHS78 = crLHS31*crLHS77 + crLHS35*crLHS77 + crLHS56*crLHS75;
        const double crLHS79 = C(0,1)*DN_vel(2,1) + crLHS69;
        const double crLHS80 = C(1,2)*DN_vel(2,1);
        const double crLHS81 = C(2,2)*DN_vel(2,0) + crLHS80;
        const double crLHS82 = DN_vel(2,1)*crLHS40;
        const double crLHS83 = DN_vel(0,0)*N[2];
        const double crLHS84 = C(0,1)*DN_vel(0,0) + crLHS38;
        const double crLHS85 = C(1,1)*DN_vel(0,1) + C(1,2)*DN_vel(0,0);
        const double crLHS86 = DN_vel(0,1)*N[0];
        const double crLHS87 = DN(0,1)*crLHS28;
        const double crLHS88 = C(0,1)*DN_vel(1,0) + crLHS62;
        const double crLHS89 = DN_vel(0,1)*crLHS3;
        const double crLHS90 = DN_vel(1,0)*crLHS89;
        const double crLHS91 = C(1,1)*DN_vel(1,1) + C(1,2)*DN_vel(1,0);
        const double crLHS92 = DN_vel(1,1)*crLHS89 + crLHS51;
        const double crLHS93 = DN_vel(0,1)*N[1];
        const double crLHS94 = C(0,1)*DN_vel(2,0) + crLHS80;
        const double crLHS95 = DN_vel(2,0)*crLHS89;
        const double crLHS96 = C(1,1)*DN_vel(2,1) + C(1,2)*DN_vel(2,0);
        const double crLHS97 = DN_vel(2,1)*crLHS89 + crLHS71;
        const double crLHS98 = DN_vel(0,1)*N[2];
        const double crLHS99 = -crLHS23/(-crLHS15*crLHS21 - crLHS16 + 0.5);
        const double crLHS100 = crLHS15*rho;
        const double crLHS101 = crLHS100*(crLHS20*crLHS99 + crLHS25);
        const double crLHS102 = 1.0*gauss_weight;
        const double crLHS103 = crLHS102*crLHS28;
        const double crLHS104 = DN_vel(1,0)*N[0];
        const double crLHS105 = crLHS19*crLHS99;
        const double crLHS106 = crLHS100*(N_vel[1]*crLHS105 + crLHS58);
        const double crLHS107 = DN_vel(1,1)*N[0];
        const double crLHS108 = crLHS103*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
        const double crLHS109 = DN_vel(2,0)*N[0];
        const double crLHS110 = crLHS100*(N_vel[2]*crLHS105 + crLHS76);
        const double crLHS111 = DN_vel(2,1)*N[0];
        const double crLHS112 = crLHS103*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
        const double crLHS113 = DN_vel(1,0)*crLHS6;
        const double crLHS114 = DN_vel(1,1)*crLHS6;
        const double crLHS115 = N_vel[1]*crLHS16;
        const double crLHS116 = crLHS30*crLHS55;
        const double crLHS117 = crLHS115*crLHS34;
        const double crLHS118 = crLHS115*crLHS14 + crLHS116*crLHS26 + crLHS117*crLHS26;
        const double crLHS119 = N_vel[1]*crLHS44;
        const double crLHS120 = crLHS55*rho;
        const double crLHS121 = 1.0*crLHS120;
        const double crLHS122 = pow(N_vel[1], 2)*crLHS32 + crLHS115*crLHS120 + crLHS116*crLHS59 + crLHS117*crLHS59;
        const double crLHS123 = DN_vel(1,0)*crLHS3;
        const double crLHS124 = DN_vel(1,1)*crLHS123;
        const double crLHS125 = DN_vel(1,0)*N[1];
        const double crLHS126 = crLHS119*crLHS28;
        const double crLHS127 = crLHS121*crLHS28;
        const double crLHS128 = N_vel[1]*N_vel[2]*crLHS32;
        const double crLHS129 = DN_vel(2,0)*crLHS123 + crLHS128;
        const double crLHS130 = crLHS75*rho;
        const double crLHS131 = crLHS115*crLHS130 + crLHS116*crLHS77 + crLHS117*crLHS77;
        const double crLHS132 = DN_vel(2,1)*crLHS123;
        const double crLHS133 = DN_vel(1,0)*N[2];
        const double crLHS134 = DN_vel(1,1)*N[1];
        const double crLHS135 = DN_vel(1,1)*crLHS3;
        const double crLHS136 = DN_vel(2,0)*crLHS135;
        const double crLHS137 = DN_vel(2,1)*crLHS135 + crLHS128;
        const double crLHS138 = DN_vel(1,1)*N[2];
        const double crLHS139 = crLHS101*crLHS28;
        const double crLHS140 = crLHS106*crLHS28;
        const double crLHS141 = DN_vel(2,0)*N[1];
        const double crLHS142 = crLHS110*crLHS28;
        const double crLHS143 = DN_vel(2,1)*N[1];
        const double crLHS144 = crLHS103*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
        const double crLHS145 = DN_vel(2,0)*crLHS6;
        const double crLHS146 = DN_vel(2,1)*crLHS6;
        const double crLHS147 = N_vel[2]*crLHS16;
        const double crLHS148 = crLHS30*crLHS75;
        const double crLHS149 = crLHS147*crLHS34;
        const double crLHS150 = crLHS14*crLHS147 + crLHS148*crLHS26 + crLHS149*crLHS26;
        const double crLHS151 = N_vel[2]*crLHS44;
        const double crLHS152 = 1.0*crLHS130;
        const double crLHS153 = crLHS120*crLHS147 + crLHS148*crLHS59 + crLHS149*crLHS59;
        const double crLHS154 = crLHS151*crLHS28;
        const double crLHS155 = crLHS152*crLHS28;
        const double crLHS156 = pow(N_vel[2], 2)*crLHS32 + crLHS130*crLHS147 + crLHS148*crLHS77 + crLHS149*crLHS77;
        const double crLHS157 = DN_vel(2,0)*DN_vel(2,1)*crLHS3;
        const double crLHS158 = DN_vel(2,0)*N[2];
        const double crLHS159 = DN_vel(2,1)*N[2];
        rLHS(0,0)+=gauss_weight*(pow(DN_vel(0,0), 2)*crLHS3 + crLHS10*crLHS9 + crLHS36 + crLHS4*crLHS7);
        rLHS(0,1)+=gauss_weight*(crLHS10*crLHS39 + crLHS37*crLHS7 + crLHS41);
        rLHS(0,2)+=gauss_weight*(-crLHS42 + crLHS43*crLHS45 + crLHS43*crLHS46);
        rLHS(0,3)+=gauss_weight*(crLHS10*crLHS49 + crLHS47*crLHS7 + crLHS52 + crLHS60);
        rLHS(0,4)+=gauss_weight*(crLHS10*crLHS63 + crLHS61*crLHS7 + crLHS64);
        rLHS(0,5)+=gauss_weight*(DN(1,0)*crLHS66 + DN(1,0)*crLHS67 - crLHS65);
        rLHS(0,6)+=gauss_weight*(crLHS10*crLHS70 + crLHS68*crLHS7 + crLHS72 + crLHS78);
        rLHS(0,7)+=gauss_weight*(crLHS10*crLHS81 + crLHS7*crLHS79 + crLHS82);
        rLHS(0,8)+=gauss_weight*(DN(2,0)*crLHS66 + DN(2,0)*crLHS67 - crLHS83);
        rLHS(1,0)+=gauss_weight*(crLHS10*crLHS84 + crLHS41 + crLHS7*crLHS9);
        rLHS(1,1)+=gauss_weight*(pow(DN_vel(0,1), 2)*crLHS3 + crLHS10*crLHS85 + crLHS36 + crLHS39*crLHS7);
        rLHS(1,2)+=gauss_weight*(DN(0,1)*crLHS67 + crLHS45*crLHS87 - crLHS86);
        rLHS(1,3)+=gauss_weight*(crLHS10*crLHS88 + crLHS49*crLHS7 + crLHS90);
        rLHS(1,4)+=gauss_weight*(crLHS10*crLHS91 + crLHS60 + crLHS63*crLHS7 + crLHS92);
        rLHS(1,5)+=gauss_weight*(DN(1,1)*crLHS66 + DN(1,1)*crLHS67 - crLHS93);
        rLHS(1,6)+=gauss_weight*(crLHS10*crLHS94 + crLHS7*crLHS70 + crLHS95);
        rLHS(1,7)+=gauss_weight*(crLHS10*crLHS96 + crLHS7*crLHS81 + crLHS78 + crLHS97);
        rLHS(1,8)+=gauss_weight*(DN(2,1)*crLHS66 + DN(2,1)*crLHS67 - crLHS98);
        rLHS(2,0)+=crLHS102*(crLHS101*crLHS43 + crLHS42*crLHS5);
        rLHS(2,1)+=crLHS102*(crLHS101*crLHS87 + crLHS5*crLHS86);
        rLHS(2,2)+=crLHS103*(pow(DN(0,0), 2) + pow(DN(0,1), 2));
        rLHS(2,3)+=crLHS102*(crLHS104*crLHS5 + crLHS106*crLHS43);
        rLHS(2,4)+=crLHS102*(crLHS106*crLHS87 + crLHS107*crLHS5);
        rLHS(2,5)+=crLHS108;
        rLHS(2,6)+=crLHS102*(crLHS109*crLHS5 + crLHS110*crLHS43);
        rLHS(2,7)+=crLHS102*(crLHS110*crLHS87 + crLHS111*crLHS5);
        rLHS(2,8)+=crLHS112;
        rLHS(3,0)+=gauss_weight*(crLHS113*crLHS4 + crLHS114*crLHS9 + crLHS118 + crLHS52);
        rLHS(3,1)+=gauss_weight*(crLHS113*crLHS37 + crLHS114*crLHS39 + crLHS90);
        rLHS(3,2)+=gauss_weight*(-crLHS104 + crLHS119*crLHS43 + crLHS121*crLHS43);
        rLHS(3,3)+=gauss_weight*(pow(DN_vel(1,0), 2)*crLHS3 + crLHS113*crLHS47 + crLHS114*crLHS49 + crLHS122);
        rLHS(3,4)+=gauss_weight*(crLHS113*crLHS61 + crLHS114*crLHS63 + crLHS124);
        rLHS(3,5)+=gauss_weight*(DN(1,0)*crLHS126 + DN(1,0)*crLHS127 - crLHS125);
        rLHS(3,6)+=gauss_weight*(crLHS113*crLHS68 + crLHS114*crLHS70 + crLHS129 + crLHS131);
        rLHS(3,7)+=gauss_weight*(crLHS113*crLHS79 + crLHS114*crLHS81 + crLHS132);
        rLHS(3,8)+=gauss_weight*(DN(2,0)*crLHS126 + DN(2,0)*crLHS127 - crLHS133);
        rLHS(4,0)+=gauss_weight*(crLHS113*crLHS9 + crLHS114*crLHS84 + crLHS64);
        rLHS(4,1)+=gauss_weight*(crLHS113*crLHS39 + crLHS114*crLHS85 + crLHS118 + crLHS92);
        rLHS(4,2)+=gauss_weight*(-crLHS107 + crLHS119*crLHS87 + crLHS121*crLHS87);
        rLHS(4,3)+=gauss_weight*(crLHS113*crLHS49 + crLHS114*crLHS88 + crLHS124);
        rLHS(4,4)+=gauss_weight*(pow(DN_vel(1,1), 2)*crLHS3 + crLHS113*crLHS63 + crLHS114*crLHS91 + crLHS122);
        rLHS(4,5)+=gauss_weight*(DN(1,1)*crLHS126 + DN(1,1)*crLHS127 - crLHS134);
        rLHS(4,6)+=gauss_weight*(crLHS113*crLHS70 + crLHS114*crLHS94 + crLHS136);
        rLHS(4,7)+=gauss_weight*(crLHS113*crLHS81 + crLHS114*crLHS96 + crLHS131 + crLHS137);
        rLHS(4,8)+=gauss_weight*(DN(2,1)*crLHS126 + DN(2,1)*crLHS127 - crLHS138);
        rLHS(5,0)+=crLHS102*(DN(1,0)*crLHS139 + crLHS5*crLHS65);
        rLHS(5,1)+=crLHS102*(DN(1,1)*crLHS139 + crLHS5*crLHS93);
        rLHS(5,2)+=crLHS108;
        rLHS(5,3)+=crLHS102*(DN(1,0)*crLHS140 + crLHS125*crLHS5);
        rLHS(5,4)+=crLHS102*(DN(1,1)*crLHS140 + crLHS134*crLHS5);
        rLHS(5,5)+=crLHS103*(pow(DN(1,0), 2) + pow(DN(1,1), 2));
        rLHS(5,6)+=crLHS102*(DN(1,0)*crLHS142 + crLHS141*crLHS5);
        rLHS(5,7)+=crLHS102*(DN(1,1)*crLHS142 + crLHS143*crLHS5);
        rLHS(5,8)+=crLHS144;
        rLHS(6,0)+=gauss_weight*(crLHS145*crLHS4 + crLHS146*crLHS9 + crLHS150 + crLHS72);
        rLHS(6,1)+=gauss_weight*(crLHS145*crLHS37 + crLHS146*crLHS39 + crLHS95);
        rLHS(6,2)+=gauss_weight*(-crLHS109 + crLHS151*crLHS43 + crLHS152*crLHS43);
        rLHS(6,3)+=gauss_weight*(crLHS129 + crLHS145*crLHS47 + crLHS146*crLHS49 + crLHS153);
        rLHS(6,4)+=gauss_weight*(crLHS136 + crLHS145*crLHS61 + crLHS146*crLHS63);
        rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS154 + DN(1,0)*crLHS155 - crLHS141);
        rLHS(6,6)+=gauss_weight*(pow(DN_vel(2,0), 2)*crLHS3 + crLHS145*crLHS68 + crLHS146*crLHS70 + crLHS156);
        rLHS(6,7)+=gauss_weight*(crLHS145*crLHS79 + crLHS146*crLHS81 + crLHS157);
        rLHS(6,8)+=gauss_weight*(DN(2,0)*crLHS154 + DN(2,0)*crLHS155 - crLHS158);
        rLHS(7,0)+=gauss_weight*(crLHS145*crLHS9 + crLHS146*crLHS84 + crLHS82);
        rLHS(7,1)+=gauss_weight*(crLHS145*crLHS39 + crLHS146*crLHS85 + crLHS150 + crLHS97);
        rLHS(7,2)+=gauss_weight*(-crLHS111 + crLHS151*crLHS87 + crLHS152*crLHS87);
        rLHS(7,3)+=gauss_weight*(crLHS132 + crLHS145*crLHS49 + crLHS146*crLHS88);
        rLHS(7,4)+=gauss_weight*(crLHS137 + crLHS145*crLHS63 + crLHS146*crLHS91 + crLHS153);
        rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS154 + DN(1,1)*crLHS155 - crLHS143);
        rLHS(7,6)+=gauss_weight*(crLHS145*crLHS70 + crLHS146*crLHS94 + crLHS157);
        rLHS(7,7)+=gauss_weight*(pow(DN_vel(2,1), 2)*crLHS3 + crLHS145*crLHS81 + crLHS146*crLHS96 + crLHS156);
        rLHS(7,8)+=gauss_weight*(DN(2,1)*crLHS154 + DN(2,1)*crLHS155 - crLHS159);
        rLHS(8,0)+=crLHS102*(DN(2,0)*crLHS139 + crLHS5*crLHS83);
        rLHS(8,1)+=crLHS102*(DN(2,1)*crLHS139 + crLHS5*crLHS98);
        rLHS(8,2)+=crLHS112;
        rLHS(8,3)+=crLHS102*(DN(2,0)*crLHS140 + crLHS133*crLHS5);
        rLHS(8,4)+=crLHS102*(DN(2,1)*crLHS140 + crLHS138*crLHS5);
        rLHS(8,5)+=crLHS144;
        rLHS(8,6)+=crLHS102*(DN(2,0)*crLHS142 + crLHS158*crLHS5);
        rLHS(8,7)+=crLHS102*(DN(2,1)*crLHS142 + crLHS159*crLHS5);
        rLHS(8,8)+=crLHS103*(pow(DN(2,0), 2) + pow(DN(2,1), 2));

}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Add current Gauss point LHS contribution
    const double gauss_weight = rData.Weight;
    //substitute_lhs_3D
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeErrorRate;

    // Add current Gauss point RHS contribution
    const double gauss_weight = rData.Weight;
            const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
        const double crRHS1 = 1.0/(max_spectral_radius + 1.0);
        const double crRHS2 = rho*(N_vel[0]*(crRHS1*(1.0*f(0,0) - 1.0*fn(0,0)) + fn(0,0)) + N_vel[1]*(crRHS1*(1.0*f(1,0) - 1.0*fn(1,0)) + fn(1,0)) + N_vel[2]*(crRHS1*(1.0*f(2,0) - 1.0*fn(2,0)) + fn(2,0)));
        const double crRHS3 = N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0);
        const double crRHS4 = N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1);
        const double crRHS5 = rho*stab_c2*sqrt(pow(crRHS3, 2) + pow(crRHS4, 2));
        const double crRHS6 = DN_vel(0,0)*v(0,0);
        const double crRHS7 = DN_vel(0,1)*v(0,1);
        const double crRHS8 = DN_vel(1,0)*v(1,0);
        const double crRHS9 = DN_vel(1,1)*v(1,1);
        const double crRHS10 = DN_vel(2,0)*v(2,0);
        const double crRHS11 = DN_vel(2,1)*v(2,1);
        const double crRHS12 = (crRHS5*h/stab_c1 + mu)*(crRHS10 + crRHS11 + crRHS6 + crRHS7 + crRHS8 + crRHS9 - volume_error_ratio);
        const double crRHS13 = 1.0*vn(0,0);
        const double crRHS14 = crRHS1*(-crRHS13 + 1.0*v(0,0)) + vn(0,0);
        const double crRHS15 = 1.0*vn(1,0);
        const double crRHS16 = crRHS1*(-crRHS15 + 1.0*v(1,0)) + vn(1,0);
        const double crRHS17 = 1.0*vn(2,0);
        const double crRHS18 = crRHS1*(-crRHS17 + 1.0*v(2,0)) + vn(2,0);
        const double crRHS19 = rho*(crRHS3*(DN_vel(0,0)*crRHS14 + DN_vel(1,0)*crRHS16 + DN_vel(2,0)*crRHS18) + crRHS4*(DN_vel(0,1)*crRHS14 + DN_vel(1,1)*crRHS16 + DN_vel(2,1)*crRHS18));
        const double crRHS20 = 1.0/dt;
        const double crRHS21 = v(0,0) - vn(0,0);
        const double crRHS22 = 1.0*crRHS1;
        const double crRHS23 = 0.5*max_spectral_radius - 1.5;
        const double crRHS24 = 1.0/(-crRHS1*crRHS23 - crRHS22 + 0.5);
        const double crRHS25 = -crRHS1*crRHS23;
        const double crRHS26 = crRHS22 - crRHS25 + 0.5;
        const double crRHS27 = crRHS24*crRHS26;
        const double crRHS28 = max_spectral_radius - 3.0;
        const double crRHS29 = 0.5*crRHS1;
        const double crRHS30 = crRHS28*crRHS29;
        const double crRHS31 = v(1,0) - vn(1,0);
        const double crRHS32 = v(2,0) - vn(2,0);
        const double crRHS33 = N_vel[0]*(acceleration_alpha_method(0,0) - crRHS30*(-acceleration_alpha_method(0,0)*crRHS27 - acceleration_alpha_method(0,0) + crRHS20*crRHS21*crRHS24)) + N_vel[1]*(acceleration_alpha_method(1,0) - crRHS30*(-acceleration_alpha_method(1,0)*crRHS27 - acceleration_alpha_method(1,0) + crRHS20*crRHS24*crRHS31)) + N_vel[2]*(acceleration_alpha_method(2,0) - crRHS30*(-acceleration_alpha_method(2,0)*crRHS27 - acceleration_alpha_method(2,0) + crRHS20*crRHS24*crRHS32));
        const double crRHS34 = N_vel[0]*rho;
        const double crRHS35 = 1.0/(-crRHS22 + crRHS25 + 0.5);
        const double crRHS36 = crRHS20*crRHS35;
        const double crRHS37 = -crRHS26*crRHS35;
        const double crRHS38 = -crRHS28*crRHS29;
        const double crRHS39 = 1.0/(crRHS20*dyn_tau*rho + crRHS5/h + mu*stab_c1/pow(h, 2));
        const double crRHS40 = crRHS39*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crRHS19 - crRHS2 + rho*(N_vel[0]*(acceleration_alpha_method(0,0) + crRHS38*(acceleration_alpha_method(0,0)*crRHS37 - acceleration_alpha_method(0,0) + crRHS21*crRHS36)) + N_vel[1]*(acceleration_alpha_method(1,0) + crRHS38*(acceleration_alpha_method(1,0)*crRHS37 - acceleration_alpha_method(1,0) + crRHS31*crRHS36)) + N_vel[2]*(acceleration_alpha_method(2,0) + crRHS38*(acceleration_alpha_method(2,0)*crRHS37 - acceleration_alpha_method(2,0) + crRHS32*crRHS36))));
        const double crRHS41 = DN_vel(0,0)*vconv(0,0) + DN_vel(0,1)*vconv(0,1) + DN_vel(1,0)*vconv(1,0) + DN_vel(1,1)*vconv(1,1) + DN_vel(2,0)*vconv(2,0) + DN_vel(2,1)*vconv(2,1);
        const double crRHS42 = crRHS34*crRHS41;
        const double crRHS43 = rho*(DN_vel(0,0)*crRHS3 + DN_vel(0,1)*crRHS4);
        const double crRHS44 = rho*(N_vel[0]*(crRHS1*(1.0*f(0,1) - 1.0*fn(0,1)) + fn(0,1)) + N_vel[1]*(crRHS1*(1.0*f(1,1) - 1.0*fn(1,1)) + fn(1,1)) + N_vel[2]*(crRHS1*(1.0*f(2,1) - 1.0*fn(2,1)) + fn(2,1)));
        const double crRHS45 = 1.0*vn(0,1);
        const double crRHS46 = crRHS1*(-crRHS45 + 1.0*v(0,1)) + vn(0,1);
        const double crRHS47 = 1.0*vn(1,1);
        const double crRHS48 = crRHS1*(-crRHS47 + 1.0*v(1,1)) + vn(1,1);
        const double crRHS49 = 1.0*vn(2,1);
        const double crRHS50 = crRHS1*(-crRHS49 + 1.0*v(2,1)) + vn(2,1);
        const double crRHS51 = rho*(crRHS3*(DN_vel(0,0)*crRHS46 + DN_vel(1,0)*crRHS48 + DN_vel(2,0)*crRHS50) + crRHS4*(DN_vel(0,1)*crRHS46 + DN_vel(1,1)*crRHS48 + DN_vel(2,1)*crRHS50));
        const double crRHS52 = v(0,1) - vn(0,1);
        const double crRHS53 = v(1,1) - vn(1,1);
        const double crRHS54 = v(2,1) - vn(2,1);
        const double crRHS55 = N_vel[0]*(acceleration_alpha_method(0,1) - crRHS30*(-acceleration_alpha_method(0,1)*crRHS27 - acceleration_alpha_method(0,1) + crRHS20*crRHS24*crRHS52)) + N_vel[1]*(acceleration_alpha_method(1,1) - crRHS30*(-acceleration_alpha_method(1,1)*crRHS27 - acceleration_alpha_method(1,1) + crRHS20*crRHS24*crRHS53)) + N_vel[2]*(acceleration_alpha_method(2,1) - crRHS30*(-acceleration_alpha_method(2,1)*crRHS27 - acceleration_alpha_method(2,1) + crRHS20*crRHS24*crRHS54));
        const double crRHS56 = crRHS39*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS44 + crRHS51 + rho*(N_vel[0]*(acceleration_alpha_method(0,1) + crRHS38*(acceleration_alpha_method(0,1)*crRHS37 - acceleration_alpha_method(0,1) + crRHS36*crRHS52)) + N_vel[1]*(acceleration_alpha_method(1,1) + crRHS38*(acceleration_alpha_method(1,1)*crRHS37 - acceleration_alpha_method(1,1) + crRHS36*crRHS53)) + N_vel[2]*(acceleration_alpha_method(2,1) + crRHS38*(acceleration_alpha_method(2,1)*crRHS37 - acceleration_alpha_method(2,1) + crRHS36*crRHS54))));
        const double crRHS57 = volume_error_ratio - (DN_vel(0,0)*crRHS13*max_spectral_radius + DN_vel(0,1)*crRHS45*max_spectral_radius + DN_vel(1,0)*crRHS15*max_spectral_radius + DN_vel(1,1)*crRHS47*max_spectral_radius + DN_vel(2,0)*crRHS17*max_spectral_radius + DN_vel(2,1)*crRHS49*max_spectral_radius + 1.0*crRHS10 + 1.0*crRHS11 + 1.0*crRHS6 + 1.0*crRHS7 + 1.0*crRHS8 + 1.0*crRHS9)/(max_spectral_radius + 1);
        const double crRHS58 = N_vel[1]*rho;
        const double crRHS59 = crRHS41*crRHS58;
        const double crRHS60 = rho*(DN_vel(1,0)*crRHS3 + DN_vel(1,1)*crRHS4);
        const double crRHS61 = N_vel[2]*rho;
        const double crRHS62 = crRHS41*crRHS61;
        const double crRHS63 = rho*(DN_vel(2,0)*crRHS3 + DN_vel(2,1)*crRHS4);
        rRHS[0]+=-gauss_weight*(-DN_vel(0,0)*crRHS0 + DN_vel(0,0)*crRHS12 + DN_vel(0,0)*stress[0] + DN_vel(0,1)*stress[2] + N_vel[0]*crRHS19 - N_vel[0]*crRHS2 + crRHS33*crRHS34 + crRHS40*crRHS42 + crRHS40*crRHS43);
        rRHS[1]+=-gauss_weight*(DN_vel(0,0)*stress[2] - DN_vel(0,1)*crRHS0 + DN_vel(0,1)*crRHS12 + DN_vel(0,1)*stress[1] - N_vel[0]*crRHS44 + N_vel[0]*crRHS51 + crRHS34*crRHS55 + crRHS42*crRHS56 + crRHS43*crRHS56);
        rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS40 + DN(0,1)*crRHS56 - N[0]*crRHS57);
        rRHS[3]+=-gauss_weight*(-DN_vel(1,0)*crRHS0 + DN_vel(1,0)*crRHS12 + DN_vel(1,0)*stress[0] + DN_vel(1,1)*stress[2] + N_vel[1]*crRHS19 - N_vel[1]*crRHS2 + crRHS33*crRHS58 + crRHS40*crRHS59 + crRHS40*crRHS60);
        rRHS[4]+=-gauss_weight*(DN_vel(1,0)*stress[2] - DN_vel(1,1)*crRHS0 + DN_vel(1,1)*crRHS12 + DN_vel(1,1)*stress[1] - N_vel[1]*crRHS44 + N_vel[1]*crRHS51 + crRHS55*crRHS58 + crRHS56*crRHS59 + crRHS56*crRHS60);
        rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS40 + DN(1,1)*crRHS56 - N[1]*crRHS57);
        rRHS[6]+=-gauss_weight*(-DN_vel(2,0)*crRHS0 + DN_vel(2,0)*crRHS12 + DN_vel(2,0)*stress[0] + DN_vel(2,1)*stress[2] + N_vel[2]*crRHS19 - N_vel[2]*crRHS2 + crRHS33*crRHS61 + crRHS40*crRHS62 + crRHS40*crRHS63);
        rRHS[7]+=-gauss_weight*(DN_vel(2,0)*stress[2] - DN_vel(2,1)*crRHS0 + DN_vel(2,1)*crRHS12 + DN_vel(2,1)*stress[1] - N_vel[2]*crRHS44 + N_vel[2]*crRHS51 + crRHS55*crRHS61 + crRHS56*crRHS62 + crRHS56*crRHS63);
        rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS40 + DN(2,1)*crRHS56 - N[2]*crRHS57);

}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeErrorRate;

    // Add current Gauss point RHS contribution
    const double gauss_weight = rData.Weight;
    //substitute_rhs_3D
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Get pressure enrichment shape function values
    const auto &N_enr = rData.Nenr;
    const auto &DN_enr = rData.DN_DXenr;

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

    const double cV0 = N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0);
const double cV1 = N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1);
const double cV2 = 1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 = cV2*rho;
const double cV4 = cV3*(DN_vel(0,0)*vconv(0,0) + DN_vel(0,1)*vconv(0,1) + DN_vel(1,0)*vconv(1,0) + DN_vel(1,1)*vconv(1,1) + DN_vel(2,0)*vconv(2,0) + DN_vel(2,1)*vconv(2,1));
const double cV5 = N_vel[0]*cV4;
const double cV6 = cV3*(DN_vel(0,0)*cV0 + DN_vel(0,1)*cV1);
const double cV7 = N_vel[1]*cV4;
const double cV8 = cV3*(DN_vel(1,0)*cV0 + DN_vel(1,1)*cV1);
const double cV9 = N_vel[2]*cV4;
const double cV10 = cV3*(DN_vel(2,0)*cV0 + DN_vel(2,1)*cV1);
V(0,0)=DN_enr(0,0)*cV5 + DN_enr(0,0)*cV6 - DN_vel(0,0)*N_enr[0];
V(0,1)=DN_enr(1,0)*cV5 + DN_enr(1,0)*cV6 - DN_vel(0,0)*N_enr[1];
V(0,2)=DN_enr(2,0)*cV5 + DN_enr(2,0)*cV6 - DN_vel(0,0)*N_enr[2];
V(1,0)=DN_enr(0,1)*cV5 + DN_enr(0,1)*cV6 - DN_vel(0,1)*N_enr[0];
V(1,1)=DN_enr(1,1)*cV5 + DN_enr(1,1)*cV6 - DN_vel(0,1)*N_enr[1];
V(1,2)=DN_enr(2,1)*cV5 + DN_enr(2,1)*cV6 - DN_vel(0,1)*N_enr[2];
V(2,0)=cV2*(DN(0,0)*DN_enr(0,0) + DN(0,1)*DN_enr(0,1));
V(2,1)=cV2*(DN(0,0)*DN_enr(1,0) + DN(0,1)*DN_enr(1,1));
V(2,2)=cV2*(DN(0,0)*DN_enr(2,0) + DN(0,1)*DN_enr(2,1));
V(3,0)=DN_enr(0,0)*cV7 + DN_enr(0,0)*cV8 - DN_vel(1,0)*N_enr[0];
V(3,1)=DN_enr(1,0)*cV7 + DN_enr(1,0)*cV8 - DN_vel(1,0)*N_enr[1];
V(3,2)=DN_enr(2,0)*cV7 + DN_enr(2,0)*cV8 - DN_vel(1,0)*N_enr[2];
V(4,0)=DN_enr(0,1)*cV7 + DN_enr(0,1)*cV8 - DN_vel(1,1)*N_enr[0];
V(4,1)=DN_enr(1,1)*cV7 + DN_enr(1,1)*cV8 - DN_vel(1,1)*N_enr[1];
V(4,2)=DN_enr(2,1)*cV7 + DN_enr(2,1)*cV8 - DN_vel(1,1)*N_enr[2];
V(5,0)=cV2*(DN(1,0)*DN_enr(0,0) + DN(1,1)*DN_enr(0,1));
V(5,1)=cV2*(DN(1,0)*DN_enr(1,0) + DN(1,1)*DN_enr(1,1));
V(5,2)=cV2*(DN(1,0)*DN_enr(2,0) + DN(1,1)*DN_enr(2,1));
V(6,0)=DN_enr(0,0)*cV10 + DN_enr(0,0)*cV9 - DN_vel(2,0)*N_enr[0];
V(6,1)=DN_enr(1,0)*cV10 + DN_enr(1,0)*cV9 - DN_vel(2,0)*N_enr[1];
V(6,2)=DN_enr(2,0)*cV10 + DN_enr(2,0)*cV9 - DN_vel(2,0)*N_enr[2];
V(7,0)=DN_enr(0,1)*cV10 + DN_enr(0,1)*cV9 - DN_vel(2,1)*N_enr[0];
V(7,1)=DN_enr(1,1)*cV10 + DN_enr(1,1)*cV9 - DN_vel(2,1)*N_enr[1];
V(7,2)=DN_enr(2,1)*cV10 + DN_enr(2,1)*cV9 - DN_vel(2,1)*N_enr[2];
V(8,0)=cV2*(DN(2,0)*DN_enr(0,0) + DN(2,1)*DN_enr(0,1));
V(8,1)=cV2*(DN(2,0)*DN_enr(1,0) + DN(2,1)*DN_enr(1,1));
V(8,2)=cV2*(DN(2,0)*DN_enr(2,0) + DN(2,1)*DN_enr(2,1));


    const double cH0 = 1.0/(max_spectral_radius + 1);
const double cH1 = N_enr[0]*cH0;
const double cH2 = N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0);
const double cH3 = 1.0*cH2;
const double cH4 = N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1);
const double cH5 = 1.0*cH4;
const double cH6 = 1.0/dt;
const double cH7 = 1.0/(max_spectral_radius + 1.0);
const double cH8 = 0.5*cH6*(max_spectral_radius - 3.0)/(cH7*(0.5*max_spectral_radius - 1.5) + 1.0*cH7 - 0.5);
const double cH9 = 1.0/(cH6*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH2, 2) + pow(cH4, 2))/h + mu*stab_c1/pow(h, 2));
const double cH10 = cH7*cH9*rho;
const double cH11 = cH10*(DN_vel(0,0)*cH3 + DN_vel(0,1)*cH5 + N_vel[0]*cH8);
const double cH12 = 1.0*cH9;
const double cH13 = cH10*(DN_vel(1,0)*cH3 + DN_vel(1,1)*cH5 + N_vel[1]*cH8);
const double cH14 = cH10*(DN_vel(2,0)*cH3 + DN_vel(2,1)*cH5 + N_vel[2]*cH8);
const double cH15 = N_enr[1]*cH0;
const double cH16 = N_enr[2]*cH0;
H(0,0)=1.0*DN_enr(0,0)*cH11 + 1.0*DN_vel(0,0)*cH1;
H(0,1)=1.0*DN_enr(0,1)*cH11 + 1.0*DN_vel(0,1)*cH1;
H(0,2)=cH12*(DN(0,0)*DN_enr(0,0) + DN(0,1)*DN_enr(0,1));
H(0,3)=1.0*DN_enr(0,0)*cH13 + 1.0*DN_vel(1,0)*cH1;
H(0,4)=1.0*DN_enr(0,1)*cH13 + 1.0*DN_vel(1,1)*cH1;
H(0,5)=cH12*(DN(1,0)*DN_enr(0,0) + DN(1,1)*DN_enr(0,1));
H(0,6)=1.0*DN_enr(0,0)*cH14 + 1.0*DN_vel(2,0)*cH1;
H(0,7)=1.0*DN_enr(0,1)*cH14 + 1.0*DN_vel(2,1)*cH1;
H(0,8)=cH12*(DN(2,0)*DN_enr(0,0) + DN(2,1)*DN_enr(0,1));
H(1,0)=1.0*DN_enr(1,0)*cH11 + 1.0*DN_vel(0,0)*cH15;
H(1,1)=1.0*DN_enr(1,1)*cH11 + 1.0*DN_vel(0,1)*cH15;
H(1,2)=cH12*(DN(0,0)*DN_enr(1,0) + DN(0,1)*DN_enr(1,1));
H(1,3)=1.0*DN_enr(1,0)*cH13 + 1.0*DN_vel(1,0)*cH15;
H(1,4)=1.0*DN_enr(1,1)*cH13 + 1.0*DN_vel(1,1)*cH15;
H(1,5)=cH12*(DN(1,0)*DN_enr(1,0) + DN(1,1)*DN_enr(1,1));
H(1,6)=1.0*DN_enr(1,0)*cH14 + 1.0*DN_vel(2,0)*cH15;
H(1,7)=1.0*DN_enr(1,1)*cH14 + 1.0*DN_vel(2,1)*cH15;
H(1,8)=cH12*(DN(2,0)*DN_enr(1,0) + DN(2,1)*DN_enr(1,1));
H(2,0)=1.0*DN_enr(2,0)*cH11 + 1.0*DN_vel(0,0)*cH16;
H(2,1)=1.0*DN_enr(2,1)*cH11 + 1.0*DN_vel(0,1)*cH16;
H(2,2)=cH12*(DN(0,0)*DN_enr(2,0) + DN(0,1)*DN_enr(2,1));
H(2,3)=1.0*DN_enr(2,0)*cH13 + 1.0*DN_vel(1,0)*cH16;
H(2,4)=1.0*DN_enr(2,1)*cH13 + 1.0*DN_vel(1,1)*cH16;
H(2,5)=cH12*(DN(1,0)*DN_enr(2,0) + DN(1,1)*DN_enr(2,1));
H(2,6)=1.0*DN_enr(2,0)*cH14 + 1.0*DN_vel(2,0)*cH16;
H(2,7)=1.0*DN_enr(2,1)*cH14 + 1.0*DN_vel(2,1)*cH16;
H(2,8)=cH12*(DN(2,0)*DN_enr(2,0) + DN(2,1)*DN_enr(2,1));


    const double cKee0 = 1.0/(rho*stab_c2*sqrt(pow(N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0), 2) + pow(N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DN_enr(0,0)*DN_enr(1,0) + DN_enr(0,1)*DN_enr(1,1));
const double cKee2 = cKee0*(DN_enr(0,0)*DN_enr(2,0) + DN_enr(0,1)*DN_enr(2,1));
const double cKee3 = cKee0*(DN_enr(1,0)*DN_enr(2,0) + DN_enr(1,1)*DN_enr(2,1));
Kee(0,0)=cKee0*(pow(DN_enr(0,0), 2) + pow(DN_enr(0,1), 2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(pow(DN_enr(1,0), 2) + pow(DN_enr(1,1), 2));
Kee(1,2)=cKee3;
Kee(2,0)=cKee2;
Kee(2,1)=cKee3;
Kee(2,2)=cKee0*(pow(DN_enr(2,0), 2) + pow(DN_enr(2,1), 2));


    const double crhs_ee0 = 1.0*v(0,0);
const double crhs_ee1 = 1.0*v(0,1);
const double crhs_ee2 = 1.0*v(1,0);
const double crhs_ee3 = 1.0*v(1,1);
const double crhs_ee4 = 1.0*v(2,0);
const double crhs_ee5 = 1.0*v(2,1);
const double crhs_ee6 = 1.0*vn(0,0);
const double crhs_ee7 = 1.0*vn(0,1);
const double crhs_ee8 = 1.0*vn(1,0);
const double crhs_ee9 = 1.0*vn(1,1);
const double crhs_ee10 = 1.0*vn(2,0);
const double crhs_ee11 = 1.0*vn(2,1);
const double crhs_ee12 = volume_error_ratio - (DN_vel(0,0)*crhs_ee0 + DN_vel(0,0)*crhs_ee6*max_spectral_radius + DN_vel(0,1)*crhs_ee1 + DN_vel(0,1)*crhs_ee7*max_spectral_radius + DN_vel(1,0)*crhs_ee2 + DN_vel(1,0)*crhs_ee8*max_spectral_radius + DN_vel(1,1)*crhs_ee3 + DN_vel(1,1)*crhs_ee9*max_spectral_radius + DN_vel(2,0)*crhs_ee10*max_spectral_radius + DN_vel(2,0)*crhs_ee4 + DN_vel(2,1)*crhs_ee11*max_spectral_radius + DN_vel(2,1)*crhs_ee5)/(max_spectral_radius + 1);
const double crhs_ee13 = 1.0/(max_spectral_radius + 1.0);
const double crhs_ee14 = N_vel[0]*vconv(0,0) + N_vel[1]*vconv(1,0) + N_vel[2]*vconv(2,0);
const double crhs_ee15 = crhs_ee13*(crhs_ee0 - crhs_ee6) + vn(0,0);
const double crhs_ee16 = crhs_ee13*(crhs_ee2 - crhs_ee8) + vn(1,0);
const double crhs_ee17 = crhs_ee13*(-crhs_ee10 + crhs_ee4) + vn(2,0);
const double crhs_ee18 = N_vel[0]*vconv(0,1) + N_vel[1]*vconv(1,1) + N_vel[2]*vconv(2,1);
const double crhs_ee19 = 1.0/dt;
const double crhs_ee20 = 1.0*crhs_ee13;
const double crhs_ee21 = 1.5 - 0.5*max_spectral_radius;
const double crhs_ee22 = 1.0/(crhs_ee13*crhs_ee21 - crhs_ee20 + 0.5);
const double crhs_ee23 = crhs_ee19*crhs_ee22;
const double crhs_ee24 = crhs_ee22*(crhs_ee13*crhs_ee21 - crhs_ee20 - 0.5);
const double crhs_ee25 = 0.5*crhs_ee13*(3.0 - max_spectral_radius);
const double crhs_ee26 = 1.0/(crhs_ee19*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee14, 2) + pow(crhs_ee18, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee27 = crhs_ee26*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN_enr(0,0)*penr[0] + DN_enr(1,0)*penr[1] + DN_enr(2,0)*penr[2] - rho*(N_vel[0]*(crhs_ee13*(1.0*f(0,0) - 1.0*fn(0,0)) + fn(0,0)) + N_vel[1]*(crhs_ee13*(1.0*f(1,0) - 1.0*fn(1,0)) + fn(1,0)) + N_vel[2]*(crhs_ee13*(1.0*f(2,0) - 1.0*fn(2,0)) + fn(2,0))) + rho*(N_vel[0]*(acceleration_alpha_method(0,0) + crhs_ee25*(acceleration_alpha_method(0,0)*crhs_ee24 - acceleration_alpha_method(0,0) + crhs_ee23*(v(0,0) - vn(0,0)))) + N_vel[1]*(acceleration_alpha_method(1,0) + crhs_ee25*(acceleration_alpha_method(1,0)*crhs_ee24 - acceleration_alpha_method(1,0) + crhs_ee23*(v(1,0) - vn(1,0)))) + N_vel[2]*(acceleration_alpha_method(2,0) + crhs_ee25*(acceleration_alpha_method(2,0)*crhs_ee24 - acceleration_alpha_method(2,0) + crhs_ee23*(v(2,0) - vn(2,0)))) + crhs_ee14*(DN_vel(0,0)*crhs_ee15 + DN_vel(1,0)*crhs_ee16 + DN_vel(2,0)*crhs_ee17) + crhs_ee18*(DN_vel(0,1)*crhs_ee15 + DN_vel(1,1)*crhs_ee16 + DN_vel(2,1)*crhs_ee17)));
const double crhs_ee28 = crhs_ee13*(crhs_ee1 - crhs_ee7) + vn(0,1);
const double crhs_ee29 = crhs_ee13*(crhs_ee3 - crhs_ee9) + vn(1,1);
const double crhs_ee30 = crhs_ee13*(-crhs_ee11 + crhs_ee5) + vn(2,1);
const double crhs_ee31 = crhs_ee26*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN_enr(0,1)*penr[0] + DN_enr(1,1)*penr[1] + DN_enr(2,1)*penr[2] - rho*(N_vel[0]*(crhs_ee13*(1.0*f(0,1) - 1.0*fn(0,1)) + fn(0,1)) + N_vel[1]*(crhs_ee13*(1.0*f(1,1) - 1.0*fn(1,1)) + fn(1,1)) + N_vel[2]*(crhs_ee13*(1.0*f(2,1) - 1.0*fn(2,1)) + fn(2,1))) + rho*(N_vel[0]*(acceleration_alpha_method(0,1) + crhs_ee25*(acceleration_alpha_method(0,1)*crhs_ee24 - acceleration_alpha_method(0,1) + crhs_ee23*(v(0,1) - vn(0,1)))) + N_vel[1]*(acceleration_alpha_method(1,1) + crhs_ee25*(acceleration_alpha_method(1,1)*crhs_ee24 - acceleration_alpha_method(1,1) + crhs_ee23*(v(1,1) - vn(1,1)))) + N_vel[2]*(acceleration_alpha_method(2,1) + crhs_ee25*(acceleration_alpha_method(2,1)*crhs_ee24 - acceleration_alpha_method(2,1) + crhs_ee23*(v(2,1) - vn(2,1)))) + crhs_ee14*(DN_vel(0,0)*crhs_ee28 + DN_vel(1,0)*crhs_ee29 + DN_vel(2,0)*crhs_ee30) + crhs_ee18*(DN_vel(0,1)*crhs_ee28 + DN_vel(1,1)*crhs_ee29 + DN_vel(2,1)*crhs_ee30)));
rhs_ee[0]=-DN_enr(0,0)*crhs_ee27 - DN_enr(0,1)*crhs_ee31 + N_enr[0]*crhs_ee12;
rhs_ee[1]=-DN_enr(1,0)*crhs_ee27 - DN_enr(1,1)*crhs_ee31 + N_enr[1]*crhs_ee12;
rhs_ee[2]=-DN_enr(2,0)*crhs_ee27 - DN_enr(2,1)*crhs_ee31 + N_enr[2]*crhs_ee12;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4> &rData,
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

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get discontinuous shape functions values if the element is intersected
    const auto &N_vel = rData.IsCut() ? rData.N_vel : rData.N;
    const auto &DN_vel = rData.IsCut() ? rData.DN_DX_vel : rData.DN_DX;

    // Get pressure enrichment shape function values
    const auto &N_enr = rData.Nenr;
    const auto &DN_enr = rData.DN_DXenr;

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

    //substitute_enrichment_V_3D

    //substitute_enrichment_H_3D

    //substitute_enrichment_Kee_3D

    //substitute_enrichment_rhs_ee_3D

    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::ComputeSplitting(
    TElementData &rData,
    MatrixType &rStandardShapeFunctionsPos,
    MatrixType &rStandardShapeFunctionsNeg,
    MatrixType &rAusasShapeFunctionsPos,
    MatrixType &rAusasShapeFunctionsNeg,
    MatrixType &rEnrichedShapeFunctionsPos,
    MatrixType &rEnrichedShapeFunctionsNeg,
    GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rAusasShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rAusasShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg)
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

    // Set the standard and Ausas modified shape functions pointers
    auto p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::UniquePointer p_stand_mod_sh_funcs;
    ModifiedShapeFunctions::UniquePointer p_ausas_mod_sh_funcs;
    if constexpr (TElementData::Dim == 2 && TElementData::NumNodes == 3) {
        p_stand_mod_sh_funcs = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
        p_ausas_mod_sh_funcs = Kratos::make_unique<Triangle2D3AusasModifiedShapeFunctions>(p_geom, rData.Distance);
    } else if (TElementData::Dim == 3 && TElementData::NumNodes == 4) {
        p_stand_mod_sh_funcs = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);
        p_ausas_mod_sh_funcs = Kratos::make_unique<Tetrahedra3D4AusasModifiedShapeFunctions>(p_geom, rData.Distance);
    } else {
        KRATOS_ERROR << "This formulation only supports simplicial elements." << std::endl;
    }

    // Call the positive side standard modified shape functions calculator
    p_stand_mod_sh_funcs->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rStandardShapeFunctionsPos,
        rStandardShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the negative side standard modified shape functions calculator
    p_stand_mod_sh_funcs->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rStandardShapeFunctionsNeg,
        rStandardShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the positive side Ausas modified shape functions calculator
    // Note that we overwrite the previously computed weights as these are the same
    p_ausas_mod_sh_funcs->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rAusasShapeFunctionsPos,
        rAusasShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the negative side Ausas modified shape functions calculator
    // Note that we overwrite the previously computed weights as these are the same
    p_ausas_mod_sh_funcs->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rAusasShapeFunctionsNeg,
        rAusasShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rStandardShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rStandardShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rStandardShapeDerivativesPos;
    for (unsigned int i = 0; i < rStandardShapeDerivativesPos.size(); ++i){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rStandardShapeDerivativesPos[i]);
    }

    rEnrichedShapeDerivativesNeg = rStandardShapeDerivativesNeg;
    for (unsigned int i = 0; i < rStandardShapeDerivativesNeg.size(); ++i){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rStandardShapeDerivativesNeg[i]);
    }

    // Get number of splitting divisions from the standard modified shape functions utility
    rData.NumberOfDivisions = (p_stand_mod_sh_funcs->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::UpdateIntegrationPointDataDiscontinuous(
    TElementData& rData,
    IndexType IntegrationPointIndex,
    double IntegrationPointWeight,
    const typename TElementData::MatrixRowType& rStandardShapeFunctions,
    const typename TElementData::ShapeDerivativesType& rStandardShapeFunctionsGradients,
    const typename TElementData::MatrixRowType& rAusasShapeFunctions,
    const typename TElementData::ShapeDerivativesType& rAusasShapeFunctionsGradients,
    const typename TElementData::MatrixRowType& rEnrichedShapeFunctions,
    const typename TElementData::ShapeDerivativesType& rEnrichedShapeFunctionsGradients) const
{
    // Update Gauss point current values
    rData.UpdateGeometryValues(
        IntegrationPointIndex,
        IntegrationPointWeight,
        rStandardShapeFunctions,
        rStandardShapeFunctionsGradients,
        rAusasShapeFunctions,
        rAusasShapeFunctionsGradients,
        rEnrichedShapeFunctions,
        rEnrichedShapeFunctionsGradients);

    // Calculate material response
    const double d_gauss = inner_prod(rData.Distance, rStandardShapeFunctions);
    if (d_gauss > 0.0) {
        rData.CalculateAirMaterialResponse();
    } else {
        this->CalculateMaterialResponse(rData);
    }
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<2, 3>>;
template class TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>;

} // namespace Kratos
