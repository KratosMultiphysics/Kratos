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
#include "custom_utilities/fluid_auxiliary_utilities.h"
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
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // Reset enrichment values to zero
    this->mVelEnrPos = ZeroVector(Dim);
    this->mVelEnrNeg = ZeroVector(Dim);
    this->mPresEnr = ZeroVector(NumNodes);
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
            MatrixType shape_functions_pos, shape_functions_neg;
            MatrixType shape_functions_enr_pos, shape_functions_enr_neg;
            VectorType shape_functions_bubble_pos, shape_functions_bubble_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;
            DenseVector<array_1d<double,Dim>> shape_derivatives_bubble_pos, shape_derivatives_bubble_neg;

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_functions_bubble_pos,
                shape_functions_bubble_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg,
                shape_derivatives_bubble_pos,
                shape_derivatives_bubble_neg);

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
                //TODO: I think we can do all these bounded
                // VectorType rhs_ee_tot = ZeroVector(NumNodes + Dim);
                // MatrixType Vtot = ZeroMatrix(LocalSize, NumNodes + Dim);
                // MatrixType Htot = ZeroMatrix(NumNodes + Dim, LocalSize);
                // MatrixType Kee_tot = ZeroMatrix(NumNodes + Dim, NumNodes + Dim);
                VectorType rhs_ee_tot = ZeroVector(NumNodes + 2*Dim);
                MatrixType Vtot = ZeroMatrix(LocalSize, NumNodes + 2*Dim);
                MatrixType Htot = ZeroMatrix(NumNodes + 2*Dim, LocalSize);
                MatrixType Kee_tot = ZeroMatrix(NumNodes + 2*Dim, NumNodes + 2*Dim);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); ++g_pos){
                    this->UpdateIntegrationPointDataDiscontinuous(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos],
                        shape_functions_bubble_pos[g_pos],
                        shape_derivatives_bubble_pos[g_pos]);
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
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg],
                        shape_functions_bubble_neg[g_neg],
                        shape_derivatives_bubble_neg[g_neg]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                // Without pressure gradient stabilization, volume ratio is checked during condensation
                // Also, without surface tension, zero pressure difference is penalized
                // this->CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
                this->CondenseAndSaveEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
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
    const double art_dyn_visc_coeff = 0.8;

    double grad_v_norm;
    //substitute_artificial_mu_grad_v_norm_2D_3N

    // Check that velocity gradient norm is non-zero
    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        //substitute_artificial_mu_2D_3N
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
    //FIXME: In here I think we should use the enrichment
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
    const BoundedMatrix<double,3,2> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& r_DN_DX_vel = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += r_DN_DX_vel(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += r_DN_DX_vel(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += r_DN_DX_vel(i,0)*velocity_alpha(i,1) + r_DN_DX_vel(i,1)*velocity_alpha(i,0);
    }

    if (FluidAuxiliaryUtilities::IsSplit(rData.Distance)) {
        const auto& DN_enr_vel = rData.DN_DX_enr_vel;
        const bool is_positive = inner_prod(rData.N, rData.Distance) < 0.0 ? false : true;
        if (is_positive) {
            const auto& r_enr_pos = this->mVelEnrPos;
            const array_1d<double,2>& DN_enr_vel_pos = is_positive ? DN_enr_vel : ZeroVector(2);
            r_strain_rate[0] += DN_enr_vel_pos[0]*r_enr_pos[0];
            r_strain_rate[1] += DN_enr_vel_pos[1]*r_enr_pos[1];
            r_strain_rate[2] += DN_enr_vel_pos[0]*r_enr_pos[1] + DN_enr_vel_pos[1]*r_enr_pos[0];
        } else {
            const auto& r_enr_neg = this->mVelEnrNeg;
            const array_1d<double,2>& DN_enr_vel_neg = is_positive ? ZeroVector(2) : DN_enr_vel;
            r_strain_rate[0] += DN_enr_vel_neg[0]*r_enr_neg[0];
            r_strain_rate[1] += DN_enr_vel_neg[1]*r_enr_neg[1];
            r_strain_rate[2] += DN_enr_vel_neg[0]*r_enr_neg[1] + DN_enr_vel_neg[1]*r_enr_neg[0];
        }
    }
}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>& rData) const
{
    //FIXME: In here I think we should use the enrichment
    const double alpha_f = 1/(1+rData.MaxSpectralRadius);
    const BoundedMatrix<double,4,3> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& r_DN_DX_vel = rData.DN_DX;
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

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Add current Gauss point LHS contribution
    const double gauss_weight = rData.Weight;
            const double crLHS0 = pow(DN(0,0), 2);
        const double crLHS1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crLHS2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crLHS3 = rho*stab_c2*sqrt(pow(crLHS1, 2) + pow(crLHS2, 2));
        const double crLHS4 = crLHS3*h/stab_c1 + mu;
        const double crLHS5 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
        const double crLHS6 = 1.0/(max_spectral_radius + 1);
        const double crLHS7 = 1.0*crLHS6;
        const double crLHS8 = DN(0,0)*crLHS7;
        const double crLHS9 = C(0,2)*DN(0,0);
        const double crLHS10 = C(2,2)*DN(0,1) + crLHS9;
        const double crLHS11 = DN(0,1)*crLHS7;
        const double crLHS12 = DN(0,0)*crLHS1;
        const double crLHS13 = DN(0,1)*crLHS2;
        const double crLHS14 = crLHS12 + crLHS13;
        const double crLHS15 = 1.0/(max_spectral_radius + 1.0);
        const double crLHS16 = 1.0*crLHS15;
        const double crLHS17 = crLHS16*rho;
        const double crLHS18 = N[0]*crLHS17;
        const double crLHS19 = max_spectral_radius - 3.0;
        const double crLHS20 = N[0]*crLHS19;
        const double crLHS21 = 0.5*max_spectral_radius - 1.5;
        const double crLHS22 = crLHS15*crLHS21 + crLHS16 - 0.5;
        const double crLHS23 = 1.0/dt;
        const double crLHS24 = 0.5*crLHS23;
        const double crLHS25 = -crLHS24/crLHS22;
        const double crLHS26 = 1.0*crLHS12 + 1.0*crLHS13 - crLHS20*crLHS25;
        const double crLHS27 = crLHS23*rho;
        const double crLHS28 = 1.0/(crLHS27*dyn_tau + crLHS3/h + mu*stab_c1/pow(h, 2));
        const double crLHS29 = crLHS16*crLHS28*pow(rho, 2);
        const double crLHS30 = crLHS14*crLHS29;
        const double crLHS31 = 0.5*crLHS15*crLHS27/crLHS22;
        const double crLHS32 = crLHS19*crLHS31;
        const double crLHS33 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
        const double crLHS34 = crLHS29*crLHS33;
        const double crLHS35 = N[0]*crLHS34;
        const double crLHS36 = pow(N[0], 2)*crLHS32 + crLHS14*crLHS18 + crLHS26*crLHS30 + crLHS26*crLHS35;
        const double crLHS37 = C(0,1)*DN(0,1) + crLHS9;
        const double crLHS38 = C(1,2)*DN(0,1);
        const double crLHS39 = C(2,2)*DN(0,0) + crLHS38;
        const double crLHS40 = DN(0,0)*crLHS4;
        const double crLHS41 = DN(0,1)*crLHS40;
        const double crLHS42 = crLHS28*rho;
        const double crLHS43 = 1.0*crLHS42;
        const double crLHS44 = crLHS33*crLHS43;
        const double crLHS45 = crLHS14*crLHS43;
        const double crLHS46 = gauss_weight*(N[0]*crLHS44 - N[0] + crLHS45);
        const double crLHS47 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
        const double crLHS48 = C(0,2)*DN(1,0);
        const double crLHS49 = C(2,2)*DN(1,1) + crLHS48;
        const double crLHS50 = DN(0,0)*DN(1,0);
        const double crLHS51 = crLHS20*crLHS31;
        const double crLHS52 = N[1]*crLHS51;
        const double crLHS53 = crLHS4*crLHS50 + crLHS52;
        const double crLHS54 = DN(1,0)*crLHS1;
        const double crLHS55 = DN(1,1)*crLHS2;
        const double crLHS56 = crLHS54 + crLHS55;
        const double crLHS57 = crLHS19*crLHS25;
        const double crLHS58 = 1.0*crLHS54 + 1.0*crLHS55;
        const double crLHS59 = -N[1]*crLHS57 + crLHS58;
        const double crLHS60 = crLHS18*crLHS56 + crLHS30*crLHS59 + crLHS35*crLHS59;
        const double crLHS61 = C(0,1)*DN(1,1) + crLHS48;
        const double crLHS62 = C(1,2)*DN(1,1);
        const double crLHS63 = C(2,2)*DN(1,0) + crLHS62;
        const double crLHS64 = DN(1,1)*crLHS40;
        const double crLHS65 = DN(0,0)*N[1];
        const double crLHS66 = DN(1,0)*N[0];
        const double crLHS67 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
        const double crLHS68 = C(0,2)*DN(2,0);
        const double crLHS69 = C(2,2)*DN(2,1) + crLHS68;
        const double crLHS70 = DN(0,0)*DN(2,0);
        const double crLHS71 = N[2]*crLHS51;
        const double crLHS72 = crLHS4*crLHS70 + crLHS71;
        const double crLHS73 = DN(2,0)*crLHS1;
        const double crLHS74 = DN(2,1)*crLHS2;
        const double crLHS75 = crLHS73 + crLHS74;
        const double crLHS76 = 1.0*crLHS73 + 1.0*crLHS74;
        const double crLHS77 = -N[2]*crLHS57 + crLHS76;
        const double crLHS78 = crLHS18*crLHS75 + crLHS30*crLHS77 + crLHS35*crLHS77;
        const double crLHS79 = C(0,1)*DN(2,1) + crLHS68;
        const double crLHS80 = C(1,2)*DN(2,1);
        const double crLHS81 = C(2,2)*DN(2,0) + crLHS80;
        const double crLHS82 = DN(2,1)*crLHS40;
        const double crLHS83 = DN(0,0)*N[2];
        const double crLHS84 = DN(2,0)*N[0];
        const double crLHS85 = C(0,1)*DN(0,0) + crLHS38;
        const double crLHS86 = pow(DN(0,1), 2);
        const double crLHS87 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
        const double crLHS88 = C(0,1)*DN(1,0) + crLHS62;
        const double crLHS89 = DN(0,1)*crLHS4;
        const double crLHS90 = DN(1,0)*crLHS89;
        const double crLHS91 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
        const double crLHS92 = DN(0,1)*DN(1,1);
        const double crLHS93 = crLHS4*crLHS92 + crLHS52;
        const double crLHS94 = DN(0,1)*N[1];
        const double crLHS95 = DN(1,1)*N[0];
        const double crLHS96 = C(0,1)*DN(2,0) + crLHS80;
        const double crLHS97 = DN(2,0)*crLHS89;
        const double crLHS98 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
        const double crLHS99 = DN(0,1)*DN(2,1);
        const double crLHS100 = crLHS4*crLHS99 + crLHS71;
        const double crLHS101 = DN(0,1)*N[2];
        const double crLHS102 = DN(2,1)*N[0];
        const double crLHS103 = N[0]*crLHS6;
        const double crLHS104 = crLHS15*crLHS42;
        const double crLHS105 = crLHS104*crLHS26;
        const double crLHS106 = 1.0*gauss_weight;
        const double crLHS107 = crLHS106*(crLHS103 + crLHS105);
        const double crLHS108 = crLHS106*crLHS28;
        const double crLHS109 = -crLHS19*crLHS24/(-crLHS15*crLHS21 - crLHS16 + 0.5);
        const double crLHS110 = crLHS104*(N[1]*crLHS109 + crLHS58);
        const double crLHS111 = crLHS108*(crLHS50 + crLHS92);
        const double crLHS112 = crLHS104*(N[2]*crLHS109 + crLHS76);
        const double crLHS113 = crLHS108*(crLHS70 + crLHS99);
        const double crLHS114 = DN(1,0)*crLHS7;
        const double crLHS115 = DN(1,1)*crLHS7;
        const double crLHS116 = N[1]*crLHS17;
        const double crLHS117 = crLHS29*crLHS56;
        const double crLHS118 = N[1]*crLHS34;
        const double crLHS119 = crLHS116*crLHS14 + crLHS117*crLHS26 + crLHS118*crLHS26;
        const double crLHS120 = crLHS43*crLHS56;
        const double crLHS121 = pow(DN(1,0), 2);
        const double crLHS122 = pow(N[1], 2)*crLHS32 + crLHS116*crLHS56 + crLHS117*crLHS59 + crLHS118*crLHS59;
        const double crLHS123 = DN(1,0)*crLHS4;
        const double crLHS124 = DN(1,1)*crLHS123;
        const double crLHS125 = gauss_weight*(N[1]*crLHS44 - N[1] + crLHS120);
        const double crLHS126 = DN(1,0)*DN(2,0);
        const double crLHS127 = N[1]*N[2]*crLHS32;
        const double crLHS128 = crLHS126*crLHS4 + crLHS127;
        const double crLHS129 = crLHS116*crLHS75 + crLHS117*crLHS77 + crLHS118*crLHS77;
        const double crLHS130 = DN(2,1)*crLHS123;
        const double crLHS131 = DN(1,0)*N[2];
        const double crLHS132 = DN(2,0)*N[1];
        const double crLHS133 = pow(DN(1,1), 2);
        const double crLHS134 = DN(2,0)*crLHS4;
        const double crLHS135 = DN(1,1)*crLHS134;
        const double crLHS136 = DN(1,1)*DN(2,1);
        const double crLHS137 = crLHS127 + crLHS136*crLHS4;
        const double crLHS138 = DN(1,1)*N[2];
        const double crLHS139 = DN(2,1)*N[1];
        const double crLHS140 = N[1]*crLHS6;
        const double crLHS141 = crLHS104*crLHS59;
        const double crLHS142 = crLHS106*(crLHS140 + crLHS141);
        const double crLHS143 = crLHS108*(crLHS126 + crLHS136);
        const double crLHS144 = DN(2,0)*crLHS7;
        const double crLHS145 = DN(2,1)*crLHS7;
        const double crLHS146 = N[2]*crLHS17;
        const double crLHS147 = crLHS29*crLHS75;
        const double crLHS148 = N[2]*crLHS34;
        const double crLHS149 = crLHS14*crLHS146 + crLHS147*crLHS26 + crLHS148*crLHS26;
        const double crLHS150 = crLHS43*crLHS75;
        const double crLHS151 = crLHS146*crLHS56 + crLHS147*crLHS59 + crLHS148*crLHS59;
        const double crLHS152 = pow(DN(2,0), 2);
        const double crLHS153 = pow(N[2], 2)*crLHS32 + crLHS146*crLHS75 + crLHS147*crLHS77 + crLHS148*crLHS77;
        const double crLHS154 = DN(2,1)*crLHS134;
        const double crLHS155 = gauss_weight*(N[2]*crLHS44 - N[2] + crLHS150);
        const double crLHS156 = pow(DN(2,1), 2);
        const double crLHS157 = crLHS106*(N[2]*crLHS6 + crLHS104*crLHS77);
        rLHS(0,0)+=gauss_weight*(crLHS0*crLHS4 + crLHS10*crLHS11 + crLHS36 + crLHS5*crLHS8);
        rLHS(0,1)+=gauss_weight*(crLHS11*crLHS39 + crLHS37*crLHS8 + crLHS41);
        rLHS(0,2)+=DN(0,0)*crLHS46;
        rLHS(0,3)+=gauss_weight*(crLHS11*crLHS49 + crLHS47*crLHS8 + crLHS53 + crLHS60);
        rLHS(0,4)+=gauss_weight*(crLHS11*crLHS63 + crLHS61*crLHS8 + crLHS64);
        rLHS(0,5)+=gauss_weight*(DN(1,0)*crLHS45 + crLHS44*crLHS66 - crLHS65);
        rLHS(0,6)+=gauss_weight*(crLHS11*crLHS69 + crLHS67*crLHS8 + crLHS72 + crLHS78);
        rLHS(0,7)+=gauss_weight*(crLHS11*crLHS81 + crLHS79*crLHS8 + crLHS82);
        rLHS(0,8)+=gauss_weight*(DN(2,0)*crLHS45 + crLHS44*crLHS84 - crLHS83);
        rLHS(1,0)+=gauss_weight*(crLHS10*crLHS8 + crLHS11*crLHS85 + crLHS41);
        rLHS(1,1)+=gauss_weight*(crLHS11*crLHS87 + crLHS36 + crLHS39*crLHS8 + crLHS4*crLHS86);
        rLHS(1,2)+=DN(0,1)*crLHS46;
        rLHS(1,3)+=gauss_weight*(crLHS11*crLHS88 + crLHS49*crLHS8 + crLHS90);
        rLHS(1,4)+=gauss_weight*(crLHS11*crLHS91 + crLHS60 + crLHS63*crLHS8 + crLHS93);
        rLHS(1,5)+=gauss_weight*(DN(1,1)*crLHS45 + crLHS44*crLHS95 - crLHS94);
        rLHS(1,6)+=gauss_weight*(crLHS11*crLHS96 + crLHS69*crLHS8 + crLHS97);
        rLHS(1,7)+=gauss_weight*(crLHS100 + crLHS11*crLHS98 + crLHS78 + crLHS8*crLHS81);
        rLHS(1,8)+=gauss_weight*(DN(2,1)*crLHS45 - crLHS101 + crLHS102*crLHS44);
        rLHS(2,0)+=DN(0,0)*crLHS107;
        rLHS(2,1)+=DN(0,1)*crLHS107;
        rLHS(2,2)+=crLHS108*(crLHS0 + crLHS86);
        rLHS(2,3)+=crLHS106*(DN(0,0)*crLHS110 + DN(1,0)*crLHS103);
        rLHS(2,4)+=crLHS106*(DN(0,1)*crLHS110 + DN(1,1)*crLHS103);
        rLHS(2,5)+=crLHS111;
        rLHS(2,6)+=crLHS106*(DN(0,0)*crLHS112 + DN(2,0)*crLHS103);
        rLHS(2,7)+=crLHS106*(DN(0,1)*crLHS112 + DN(2,1)*crLHS103);
        rLHS(2,8)+=crLHS113;
        rLHS(3,0)+=gauss_weight*(crLHS10*crLHS115 + crLHS114*crLHS5 + crLHS119 + crLHS53);
        rLHS(3,1)+=gauss_weight*(crLHS114*crLHS37 + crLHS115*crLHS39 + crLHS90);
        rLHS(3,2)+=gauss_weight*(DN(0,0)*crLHS120 + crLHS44*crLHS65 - crLHS66);
        rLHS(3,3)+=gauss_weight*(crLHS114*crLHS47 + crLHS115*crLHS49 + crLHS121*crLHS4 + crLHS122);
        rLHS(3,4)+=gauss_weight*(crLHS114*crLHS61 + crLHS115*crLHS63 + crLHS124);
        rLHS(3,5)+=DN(1,0)*crLHS125;
        rLHS(3,6)+=gauss_weight*(crLHS114*crLHS67 + crLHS115*crLHS69 + crLHS128 + crLHS129);
        rLHS(3,7)+=gauss_weight*(crLHS114*crLHS79 + crLHS115*crLHS81 + crLHS130);
        rLHS(3,8)+=gauss_weight*(DN(2,0)*crLHS120 - crLHS131 + crLHS132*crLHS44);
        rLHS(4,0)+=gauss_weight*(crLHS10*crLHS114 + crLHS115*crLHS85 + crLHS64);
        rLHS(4,1)+=gauss_weight*(crLHS114*crLHS39 + crLHS115*crLHS87 + crLHS119 + crLHS93);
        rLHS(4,2)+=gauss_weight*(DN(0,1)*crLHS120 + crLHS44*crLHS94 - crLHS95);
        rLHS(4,3)+=gauss_weight*(crLHS114*crLHS49 + crLHS115*crLHS88 + crLHS124);
        rLHS(4,4)+=gauss_weight*(crLHS114*crLHS63 + crLHS115*crLHS91 + crLHS122 + crLHS133*crLHS4);
        rLHS(4,5)+=DN(1,1)*crLHS125;
        rLHS(4,6)+=gauss_weight*(crLHS114*crLHS69 + crLHS115*crLHS96 + crLHS135);
        rLHS(4,7)+=gauss_weight*(crLHS114*crLHS81 + crLHS115*crLHS98 + crLHS129 + crLHS137);
        rLHS(4,8)+=gauss_weight*(DN(2,1)*crLHS120 - crLHS138 + crLHS139*crLHS44);
        rLHS(5,0)+=crLHS106*(DN(1,0)*crLHS105 + crLHS6*crLHS65);
        rLHS(5,1)+=crLHS106*(DN(1,1)*crLHS105 + crLHS6*crLHS94);
        rLHS(5,2)+=crLHS111;
        rLHS(5,3)+=DN(1,0)*crLHS142;
        rLHS(5,4)+=DN(1,1)*crLHS142;
        rLHS(5,5)+=crLHS108*(crLHS121 + crLHS133);
        rLHS(5,6)+=crLHS106*(DN(1,0)*crLHS112 + DN(2,0)*crLHS140);
        rLHS(5,7)+=crLHS106*(DN(1,1)*crLHS112 + DN(2,1)*crLHS140);
        rLHS(5,8)+=crLHS143;
        rLHS(6,0)+=gauss_weight*(crLHS10*crLHS145 + crLHS144*crLHS5 + crLHS149 + crLHS72);
        rLHS(6,1)+=gauss_weight*(crLHS144*crLHS37 + crLHS145*crLHS39 + crLHS97);
        rLHS(6,2)+=gauss_weight*(DN(0,0)*crLHS150 + crLHS44*crLHS83 - crLHS84);
        rLHS(6,3)+=gauss_weight*(crLHS128 + crLHS144*crLHS47 + crLHS145*crLHS49 + crLHS151);
        rLHS(6,4)+=gauss_weight*(crLHS135 + crLHS144*crLHS61 + crLHS145*crLHS63);
        rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS150 + crLHS131*crLHS44 - crLHS132);
        rLHS(6,6)+=gauss_weight*(crLHS144*crLHS67 + crLHS145*crLHS69 + crLHS152*crLHS4 + crLHS153);
        rLHS(6,7)+=gauss_weight*(crLHS144*crLHS79 + crLHS145*crLHS81 + crLHS154);
        rLHS(6,8)+=DN(2,0)*crLHS155;
        rLHS(7,0)+=gauss_weight*(crLHS10*crLHS144 + crLHS145*crLHS85 + crLHS82);
        rLHS(7,1)+=gauss_weight*(crLHS100 + crLHS144*crLHS39 + crLHS145*crLHS87 + crLHS149);
        rLHS(7,2)+=gauss_weight*(DN(0,1)*crLHS150 + crLHS101*crLHS44 - crLHS102);
        rLHS(7,3)+=gauss_weight*(crLHS130 + crLHS144*crLHS49 + crLHS145*crLHS88);
        rLHS(7,4)+=gauss_weight*(crLHS137 + crLHS144*crLHS63 + crLHS145*crLHS91 + crLHS151);
        rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS150 + crLHS138*crLHS44 - crLHS139);
        rLHS(7,6)+=gauss_weight*(crLHS144*crLHS69 + crLHS145*crLHS96 + crLHS154);
        rLHS(7,7)+=gauss_weight*(crLHS144*crLHS81 + crLHS145*crLHS98 + crLHS153 + crLHS156*crLHS4);
        rLHS(7,8)+=DN(2,1)*crLHS155;
        rLHS(8,0)+=crLHS106*(DN(2,0)*crLHS105 + crLHS6*crLHS83);
        rLHS(8,1)+=crLHS106*(DN(2,1)*crLHS105 + crLHS101*crLHS6);
        rLHS(8,2)+=crLHS113;
        rLHS(8,3)+=crLHS106*(DN(2,0)*crLHS141 + crLHS131*crLHS6);
        rLHS(8,4)+=crLHS106*(DN(2,1)*crLHS141 + crLHS138*crLHS6);
        rLHS(8,5)+=crLHS143;
        rLHS(8,6)+=DN(2,0)*crLHS157;
        rLHS(8,7)+=DN(2,1)*crLHS157;
        rLHS(8,8)+=crLHS108*(crLHS152 + crLHS156);

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

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeErrorRate;

    // Add current Gauss point RHS contribution
    const double gauss_weight = rData.Weight;
            const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
        const double crRHS1 = 1.0/(max_spectral_radius + 1.0);
        const double crRHS2 = rho*(N[0]*(crRHS1*(1.0*f(0,0) - 1.0*fn(0,0)) + fn(0,0)) + N[1]*(crRHS1*(1.0*f(1,0) - 1.0*fn(1,0)) + fn(1,0)) + N[2]*(crRHS1*(1.0*f(2,0) - 1.0*fn(2,0)) + fn(2,0)));
        const double crRHS3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crRHS4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crRHS5 = rho*stab_c2*sqrt(pow(crRHS3, 2) + pow(crRHS4, 2));
        const double crRHS6 = DN(0,0)*v(0,0);
        const double crRHS7 = DN(0,1)*v(0,1);
        const double crRHS8 = DN(1,0)*v(1,0);
        const double crRHS9 = DN(1,1)*v(1,1);
        const double crRHS10 = DN(2,0)*v(2,0);
        const double crRHS11 = DN(2,1)*v(2,1);
        const double crRHS12 = (crRHS5*h/stab_c1 + mu)*(crRHS10 + crRHS11 + crRHS6 + crRHS7 + crRHS8 + crRHS9 - volume_error_ratio);
        const double crRHS13 = 1.0*vn(0,0);
        const double crRHS14 = crRHS1*(-crRHS13 + 1.0*v(0,0)) + vn(0,0);
        const double crRHS15 = 1.0*vn(1,0);
        const double crRHS16 = crRHS1*(-crRHS15 + 1.0*v(1,0)) + vn(1,0);
        const double crRHS17 = 1.0*vn(2,0);
        const double crRHS18 = crRHS1*(-crRHS17 + 1.0*v(2,0)) + vn(2,0);
        const double crRHS19 = rho*(crRHS3*(DN(0,0)*crRHS14 + DN(1,0)*crRHS16 + DN(2,0)*crRHS18) + crRHS4*(DN(0,1)*crRHS14 + DN(1,1)*crRHS16 + DN(2,1)*crRHS18));
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
        const double crRHS33 = N[0]*(acceleration_alpha_method(0,0) - crRHS30*(-acceleration_alpha_method(0,0)*crRHS27 - acceleration_alpha_method(0,0) + crRHS20*crRHS21*crRHS24)) + N[1]*(acceleration_alpha_method(1,0) - crRHS30*(-acceleration_alpha_method(1,0)*crRHS27 - acceleration_alpha_method(1,0) + crRHS20*crRHS24*crRHS31)) + N[2]*(acceleration_alpha_method(2,0) - crRHS30*(-acceleration_alpha_method(2,0)*crRHS27 - acceleration_alpha_method(2,0) + crRHS20*crRHS24*crRHS32));
        const double crRHS34 = N[0]*rho;
        const double crRHS35 = 1.0/(-crRHS22 + crRHS25 + 0.5);
        const double crRHS36 = crRHS20*crRHS35;
        const double crRHS37 = -crRHS26*crRHS35;
        const double crRHS38 = -crRHS28*crRHS29;
        const double crRHS39 = 1.0/(crRHS20*dyn_tau*rho + crRHS5/h + mu*stab_c1/pow(h, 2));
        const double crRHS40 = crRHS39*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crRHS19 - crRHS2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crRHS38*(acceleration_alpha_method(0,0)*crRHS37 - acceleration_alpha_method(0,0) + crRHS21*crRHS36)) + N[1]*(acceleration_alpha_method(1,0) + crRHS38*(acceleration_alpha_method(1,0)*crRHS37 - acceleration_alpha_method(1,0) + crRHS31*crRHS36)) + N[2]*(acceleration_alpha_method(2,0) + crRHS38*(acceleration_alpha_method(2,0)*crRHS37 - acceleration_alpha_method(2,0) + crRHS32*crRHS36))));
        const double crRHS41 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
        const double crRHS42 = crRHS34*crRHS41;
        const double crRHS43 = rho*(DN(0,0)*crRHS3 + DN(0,1)*crRHS4);
        const double crRHS44 = rho*(N[0]*(crRHS1*(1.0*f(0,1) - 1.0*fn(0,1)) + fn(0,1)) + N[1]*(crRHS1*(1.0*f(1,1) - 1.0*fn(1,1)) + fn(1,1)) + N[2]*(crRHS1*(1.0*f(2,1) - 1.0*fn(2,1)) + fn(2,1)));
        const double crRHS45 = 1.0*vn(0,1);
        const double crRHS46 = crRHS1*(-crRHS45 + 1.0*v(0,1)) + vn(0,1);
        const double crRHS47 = 1.0*vn(1,1);
        const double crRHS48 = crRHS1*(-crRHS47 + 1.0*v(1,1)) + vn(1,1);
        const double crRHS49 = 1.0*vn(2,1);
        const double crRHS50 = crRHS1*(-crRHS49 + 1.0*v(2,1)) + vn(2,1);
        const double crRHS51 = rho*(crRHS3*(DN(0,0)*crRHS46 + DN(1,0)*crRHS48 + DN(2,0)*crRHS50) + crRHS4*(DN(0,1)*crRHS46 + DN(1,1)*crRHS48 + DN(2,1)*crRHS50));
        const double crRHS52 = v(0,1) - vn(0,1);
        const double crRHS53 = v(1,1) - vn(1,1);
        const double crRHS54 = v(2,1) - vn(2,1);
        const double crRHS55 = N[0]*(acceleration_alpha_method(0,1) - crRHS30*(-acceleration_alpha_method(0,1)*crRHS27 - acceleration_alpha_method(0,1) + crRHS20*crRHS24*crRHS52)) + N[1]*(acceleration_alpha_method(1,1) - crRHS30*(-acceleration_alpha_method(1,1)*crRHS27 - acceleration_alpha_method(1,1) + crRHS20*crRHS24*crRHS53)) + N[2]*(acceleration_alpha_method(2,1) - crRHS30*(-acceleration_alpha_method(2,1)*crRHS27 - acceleration_alpha_method(2,1) + crRHS20*crRHS24*crRHS54));
        const double crRHS56 = crRHS39*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS44 + crRHS51 + rho*(N[0]*(acceleration_alpha_method(0,1) + crRHS38*(acceleration_alpha_method(0,1)*crRHS37 - acceleration_alpha_method(0,1) + crRHS36*crRHS52)) + N[1]*(acceleration_alpha_method(1,1) + crRHS38*(acceleration_alpha_method(1,1)*crRHS37 - acceleration_alpha_method(1,1) + crRHS36*crRHS53)) + N[2]*(acceleration_alpha_method(2,1) + crRHS38*(acceleration_alpha_method(2,1)*crRHS37 - acceleration_alpha_method(2,1) + crRHS36*crRHS54))));
        const double crRHS57 = volume_error_ratio - (DN(0,0)*crRHS13*max_spectral_radius + DN(0,1)*crRHS45*max_spectral_radius + DN(1,0)*crRHS15*max_spectral_radius + DN(1,1)*crRHS47*max_spectral_radius + DN(2,0)*crRHS17*max_spectral_radius + DN(2,1)*crRHS49*max_spectral_radius + 1.0*crRHS10 + 1.0*crRHS11 + 1.0*crRHS6 + 1.0*crRHS7 + 1.0*crRHS8 + 1.0*crRHS9)/(max_spectral_radius + 1);
        const double crRHS58 = N[1]*rho;
        const double crRHS59 = crRHS41*crRHS58;
        const double crRHS60 = rho*(DN(1,0)*crRHS3 + DN(1,1)*crRHS4);
        const double crRHS61 = N[2]*rho;
        const double crRHS62 = crRHS41*crRHS61;
        const double crRHS63 = rho*(DN(2,0)*crRHS3 + DN(2,1)*crRHS4);
        rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS12 + DN(0,0)*stress[0] + DN(0,1)*stress[2] + N[0]*crRHS19 - N[0]*crRHS2 + crRHS33*crRHS34 + crRHS40*crRHS42 + crRHS40*crRHS43);
        rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS12 + DN(0,1)*stress[1] - N[0]*crRHS44 + N[0]*crRHS51 + crRHS34*crRHS55 + crRHS42*crRHS56 + crRHS43*crRHS56);
        rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS40 + DN(0,1)*crRHS56 - N[0]*crRHS57);
        rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS12 + DN(1,0)*stress[0] + DN(1,1)*stress[2] + N[1]*crRHS19 - N[1]*crRHS2 + crRHS33*crRHS58 + crRHS40*crRHS59 + crRHS40*crRHS60);
        rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS12 + DN(1,1)*stress[1] - N[1]*crRHS44 + N[1]*crRHS51 + crRHS55*crRHS58 + crRHS56*crRHS59 + crRHS56*crRHS60);
        rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS40 + DN(1,1)*crRHS56 - N[1]*crRHS57);
        rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS12 + DN(2,0)*stress[0] + DN(2,1)*stress[2] + N[2]*crRHS19 - N[2]*crRHS2 + crRHS33*crRHS61 + crRHS40*crRHS62 + crRHS40*crRHS63);
        rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS12 + DN(2,1)*stress[1] - N[2]*crRHS44 + N[2]*crRHS51 + crRHS55*crRHS61 + crRHS56*crRHS62 + crRHS56*crRHS63);
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
    VectorType &rRHSee)
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

    // Get material response data
    const Matrix &C = rData.C;
    const auto &stress = rData.ShearStress;

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get pressure enrichment shape function values
    const auto &N_enr_p = rData.Nenr;
    const auto &DN_enr_p = rData.DN_DXenr;

    // Get velocity enrichment bubble function values
    const auto &N_enr_vel = rData.N_enr_vel;
    const auto &DN_enr_vel = rData.DN_DX_enr_vel;

    const bool is_positive = inner_prod(N, rData.Distance) < 0.0 ? false : true;
    const double N_enr_vel_pos = is_positive ? N_enr_vel : 0.0;
    const double N_enr_vel_neg = is_positive ? 0.0 : N_enr_vel;
    const array_1d<double,2>& DN_enr_vel_pos = is_positive ? DN_enr_vel : ZeroVector(2);
    const array_1d<double,2>& DN_enr_vel_neg = is_positive ? ZeroVector(2) : DN_enr_vel;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeErrorRate;

    // Initialize enrichment DOFs appearing in the enrichment RHS to zero
    // Note that we always initialize them to zero as we do not want to store them
    // array_1d<double, Dim> v_enr = ZeroVector(Dim);
    array_1d<double, 2> v_enr_pos = ZeroVector(2);
    array_1d<double, 2> v_enr_neg = ZeroVector(2);
    array_1d<double, 3> penr = ZeroVector(3);

    // Add current Gauss point enrichment contribution
    const double gauss_weight = rData.Weight;

            const double crV0 = N[0]*rho;
        const double crV1 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crV2 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crV3 = rho*stab_c2*sqrt(pow(crV1, 2) + pow(crV2, 2));
        const double crV4 = 1.0/(crV3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
        const double crV5 = 1.0*crV4;
        const double crV6 = crV5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
        const double crV7 = crV0*crV6;
        const double crV8 = 2.0*crV4;
        const double crV9 = crV8*(DN(0,0)*crV1 + DN(0,1)*crV2);
        const double crV10 = crV9*rho;
        const double crV11 = crV3*h/stab_c1 + mu;
        const double crV12 = DN(0,0)*crV11;
        const double crV13 = DN_enr_vel_pos[0]*crV1 + DN_enr_vel_pos[1]*crV2;
        const double crV14 = pow(rho, 2);
        const double crV15 = crV14*crV9;
        const double crV16 = crV14*crV6;
        const double crV17 = N[0]*crV16;
        const double crV18 = crV0*crV13 + crV13*crV15 + crV13*crV17;
        const double crV19 = crV12*gauss_weight;
        const double crV20 = DN_enr_vel_neg[0]*crV1 + DN_enr_vel_neg[1]*crV2;
        const double crV21 = crV0*crV20 + crV15*crV20 + crV17*crV20;
        const double crV22 = DN(0,1)*crV11;
        const double crV23 = crV22*gauss_weight;
        const double crV24 = crV5*gauss_weight;
        const double crV25 = crV5*rho;
        const double crV26 = crV13*crV25;
        const double crV27 = crV20*crV25;
        const double crV28 = N[1]*rho;
        const double crV29 = crV28*crV6;
        const double crV30 = crV8*(DN(1,0)*crV1 + DN(1,1)*crV2);
        const double crV31 = crV30*rho;
        const double crV32 = DN(1,0)*crV11;
        const double crV33 = crV14*crV30;
        const double crV34 = N[1]*crV16;
        const double crV35 = crV13*crV28 + crV13*crV33 + crV13*crV34;
        const double crV36 = crV32*gauss_weight;
        const double crV37 = crV20*crV28 + crV20*crV33 + crV20*crV34;
        const double crV38 = DN(1,1)*crV11;
        const double crV39 = crV38*gauss_weight;
        const double crV40 = N[2]*rho;
        const double crV41 = crV40*crV6;
        const double crV42 = crV8*(DN(2,0)*crV1 + DN(2,1)*crV2);
        const double crV43 = crV42*rho;
        const double crV44 = DN(2,0)*crV11;
        const double crV45 = crV14*crV42;
        const double crV46 = N[2]*crV16;
        const double crV47 = crV13*crV40 + crV13*crV45 + crV13*crV46;
        const double crV48 = crV44*gauss_weight;
        const double crV49 = crV20*crV40 + crV20*crV45 + crV20*crV46;
        const double crV50 = DN(2,1)*crV11;
        const double crV51 = crV50*gauss_weight;
        rV(0,0)+=gauss_weight*(-DN(0,0)*N_enr_p[0] + DN_enr_p(0,0)*crV10 + DN_enr_p(0,0)*crV7);
        rV(0,1)+=gauss_weight*(-DN(0,0)*N_enr_p[1] + DN_enr_p(1,0)*crV10 + DN_enr_p(1,0)*crV7);
        rV(0,2)+=gauss_weight*(-DN(0,0)*N_enr_p[2] + DN_enr_p(2,0)*crV10 + DN_enr_p(2,0)*crV7);
        rV(0,3)+=gauss_weight*(DN_enr_vel_pos[0]*crV12 + crV18);
        rV(0,4)+=DN_enr_vel_pos[1]*crV19;
        rV(0,5)+=gauss_weight*(DN_enr_vel_neg[0]*crV12 + crV21);
        rV(0,6)+=DN_enr_vel_neg[1]*crV19;
        rV(1,0)+=gauss_weight*(-DN(0,1)*N_enr_p[0] + DN_enr_p(0,1)*crV10 + DN_enr_p(0,1)*crV7);
        rV(1,1)+=gauss_weight*(-DN(0,1)*N_enr_p[1] + DN_enr_p(1,1)*crV10 + DN_enr_p(1,1)*crV7);
        rV(1,2)+=gauss_weight*(-DN(0,1)*N_enr_p[2] + DN_enr_p(2,1)*crV10 + DN_enr_p(2,1)*crV7);
        rV(1,3)+=DN_enr_vel_pos[0]*crV23;
        rV(1,4)+=gauss_weight*(DN_enr_vel_pos[1]*crV22 + crV18);
        rV(1,5)+=DN_enr_vel_neg[0]*crV23;
        rV(1,6)+=gauss_weight*(DN_enr_vel_neg[1]*crV22 + crV21);
        rV(2,0)+=crV24*(DN(0,0)*DN_enr_p(0,0) + DN(0,1)*DN_enr_p(0,1));
        rV(2,1)+=crV24*(DN(0,0)*DN_enr_p(1,0) + DN(0,1)*DN_enr_p(1,1));
        rV(2,2)+=crV24*(DN(0,0)*DN_enr_p(2,0) + DN(0,1)*DN_enr_p(2,1));
        rV(2,3)+=gauss_weight*(DN(0,0)*crV26 + DN_enr_vel_pos[0]*N[0]);
        rV(2,4)+=gauss_weight*(DN(0,1)*crV26 + DN_enr_vel_pos[1]*N[0]);
        rV(2,5)+=gauss_weight*(DN(0,0)*crV27 + DN_enr_vel_neg[0]*N[0]);
        rV(2,6)+=gauss_weight*(DN(0,1)*crV27 + DN_enr_vel_neg[1]*N[0]);
        rV(3,0)+=gauss_weight*(-DN(1,0)*N_enr_p[0] + DN_enr_p(0,0)*crV29 + DN_enr_p(0,0)*crV31);
        rV(3,1)+=gauss_weight*(-DN(1,0)*N_enr_p[1] + DN_enr_p(1,0)*crV29 + DN_enr_p(1,0)*crV31);
        rV(3,2)+=gauss_weight*(-DN(1,0)*N_enr_p[2] + DN_enr_p(2,0)*crV29 + DN_enr_p(2,0)*crV31);
        rV(3,3)+=gauss_weight*(DN_enr_vel_pos[0]*crV32 + crV35);
        rV(3,4)+=DN_enr_vel_pos[1]*crV36;
        rV(3,5)+=gauss_weight*(DN_enr_vel_neg[0]*crV32 + crV37);
        rV(3,6)+=DN_enr_vel_neg[1]*crV36;
        rV(4,0)+=gauss_weight*(-DN(1,1)*N_enr_p[0] + DN_enr_p(0,1)*crV29 + DN_enr_p(0,1)*crV31);
        rV(4,1)+=gauss_weight*(-DN(1,1)*N_enr_p[1] + DN_enr_p(1,1)*crV29 + DN_enr_p(1,1)*crV31);
        rV(4,2)+=gauss_weight*(-DN(1,1)*N_enr_p[2] + DN_enr_p(2,1)*crV29 + DN_enr_p(2,1)*crV31);
        rV(4,3)+=DN_enr_vel_pos[0]*crV39;
        rV(4,4)+=gauss_weight*(DN_enr_vel_pos[1]*crV38 + crV35);
        rV(4,5)+=DN_enr_vel_neg[0]*crV39;
        rV(4,6)+=gauss_weight*(DN_enr_vel_neg[1]*crV38 + crV37);
        rV(5,0)+=crV24*(DN(1,0)*DN_enr_p(0,0) + DN(1,1)*DN_enr_p(0,1));
        rV(5,1)+=crV24*(DN(1,0)*DN_enr_p(1,0) + DN(1,1)*DN_enr_p(1,1));
        rV(5,2)+=crV24*(DN(1,0)*DN_enr_p(2,0) + DN(1,1)*DN_enr_p(2,1));
        rV(5,3)+=gauss_weight*(DN(1,0)*crV26 + DN_enr_vel_pos[0]*N[1]);
        rV(5,4)+=gauss_weight*(DN(1,1)*crV26 + DN_enr_vel_pos[1]*N[1]);
        rV(5,5)+=gauss_weight*(DN(1,0)*crV27 + DN_enr_vel_neg[0]*N[1]);
        rV(5,6)+=gauss_weight*(DN(1,1)*crV27 + DN_enr_vel_neg[1]*N[1]);
        rV(6,0)+=gauss_weight*(-DN(2,0)*N_enr_p[0] + DN_enr_p(0,0)*crV41 + DN_enr_p(0,0)*crV43);
        rV(6,1)+=gauss_weight*(-DN(2,0)*N_enr_p[1] + DN_enr_p(1,0)*crV41 + DN_enr_p(1,0)*crV43);
        rV(6,2)+=gauss_weight*(-DN(2,0)*N_enr_p[2] + DN_enr_p(2,0)*crV41 + DN_enr_p(2,0)*crV43);
        rV(6,3)+=gauss_weight*(DN_enr_vel_pos[0]*crV44 + crV47);
        rV(6,4)+=DN_enr_vel_pos[1]*crV48;
        rV(6,5)+=gauss_weight*(DN_enr_vel_neg[0]*crV44 + crV49);
        rV(6,6)+=DN_enr_vel_neg[1]*crV48;
        rV(7,0)+=gauss_weight*(-DN(2,1)*N_enr_p[0] + DN_enr_p(0,1)*crV41 + DN_enr_p(0,1)*crV43);
        rV(7,1)+=gauss_weight*(-DN(2,1)*N_enr_p[1] + DN_enr_p(1,1)*crV41 + DN_enr_p(1,1)*crV43);
        rV(7,2)+=gauss_weight*(-DN(2,1)*N_enr_p[2] + DN_enr_p(2,1)*crV41 + DN_enr_p(2,1)*crV43);
        rV(7,3)+=DN_enr_vel_pos[0]*crV51;
        rV(7,4)+=gauss_weight*(DN_enr_vel_pos[1]*crV50 + crV47);
        rV(7,5)+=DN_enr_vel_neg[0]*crV51;
        rV(7,6)+=gauss_weight*(DN_enr_vel_neg[1]*crV50 + crV49);
        rV(8,0)+=crV24*(DN(2,0)*DN_enr_p(0,0) + DN(2,1)*DN_enr_p(0,1));
        rV(8,1)+=crV24*(DN(2,0)*DN_enr_p(1,0) + DN(2,1)*DN_enr_p(1,1));
        rV(8,2)+=crV24*(DN(2,0)*DN_enr_p(2,0) + DN(2,1)*DN_enr_p(2,1));
        rV(8,3)+=gauss_weight*(DN(2,0)*crV26 + DN_enr_vel_pos[0]*N[2]);
        rV(8,4)+=gauss_weight*(DN(2,1)*crV26 + DN_enr_vel_pos[1]*N[2]);
        rV(8,5)+=gauss_weight*(DN(2,0)*crV27 + DN_enr_vel_neg[0]*N[2]);
        rV(8,6)+=gauss_weight*(DN(2,1)*crV27 + DN_enr_vel_neg[1]*N[2]);


            const double crH0 = 1.0/(max_spectral_radius + 1);
        const double crH1 = N_enr_p[0]*crH0;
        const double crH2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crH3 = DN(0,0)*crH2;
        const double crH4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crH5 = DN(0,1)*crH4;
        const double crH6 = 0.5*max_spectral_radius - 1.5;
        const double crH7 = N[0]*crH6;
        const double crH8 = 1.0/dt;
        const double crH9 = 1.0/(max_spectral_radius + 1.0);
        const double crH10 = 1.0*crH9;
        const double crH11 = crH10 + crH9*(0.5*max_spectral_radius - 1.5) - 0.5;
        const double crH12 = -crH8/crH11;
        const double crH13 = -crH12*crH7 + 1.0*crH3 + 1.0*crH5;
        const double crH14 = crH8*rho;
        const double crH15 = rho*stab_c2*sqrt(pow(crH2, 2) + pow(crH4, 2));
        const double crH16 = 1.0/(crH14*dyn_tau + crH15/h + mu*stab_c1/pow(h, 2));
        const double crH17 = crH16*rho;
        const double crH18 = crH17*crH9;
        const double crH19 = crH13*crH18;
        const double crH20 = 1.0*gauss_weight;
        const double crH21 = crH16*crH20;
        const double crH22 = DN(1,0)*crH2;
        const double crH23 = DN(1,1)*crH4;
        const double crH24 = crH12*crH6;
        const double crH25 = -N[1]*crH24 + 1.0*crH22 + 1.0*crH23;
        const double crH26 = crH18*crH25;
        const double crH27 = DN(2,0)*crH2;
        const double crH28 = DN(2,1)*crH4;
        const double crH29 = -N[2]*crH24 + 1.0*crH27 + 1.0*crH28;
        const double crH30 = crH18*crH29;
        const double crH31 = N_enr_p[1]*crH0;
        const double crH32 = N_enr_p[2]*crH0;
        const double crH33 = crH15*h/stab_c1 + mu;
        const double crH34 = DN_enr_vel_pos[0]*crH33;
        const double crH35 = crH3 + crH5;
        const double crH36 = crH10*rho;
        const double crH37 = N_enr_vel_pos*crH36;
        const double crH38 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
        const double crH39 = N_enr_vel_pos*crH38;
        const double crH40 = crH10*crH16*pow(rho, 2);
        const double crH41 = crH39*crH40;
        const double crH42 = crH14*crH9/crH11;
        const double crH43 = N_enr_vel_pos*crH42;
        const double crH44 = crH13*crH41 + crH35*crH37 + crH43*crH7;
        const double crH45 = crH34*gauss_weight;
        const double crH46 = 1.0*crH17;
        const double crH47 = crH39*crH46;
        const double crH48 = crH22 + crH23;
        const double crH49 = crH43*crH6;
        const double crH50 = N[1]*crH49 + crH25*crH41 + crH37*crH48;
        const double crH51 = crH27 + crH28;
        const double crH52 = N[2]*crH49 + crH29*crH41 + crH37*crH51;
        const double crH53 = DN_enr_vel_pos[1]*crH33;
        const double crH54 = crH53*gauss_weight;
        const double crH55 = DN_enr_vel_neg[0]*crH33;
        const double crH56 = N_enr_vel_neg*crH36;
        const double crH57 = N_enr_vel_neg*crH38;
        const double crH58 = crH40*crH57;
        const double crH59 = N_enr_vel_neg*crH42;
        const double crH60 = crH13*crH58 + crH35*crH56 + crH59*crH7;
        const double crH61 = crH55*gauss_weight;
        const double crH62 = crH46*crH57;
        const double crH63 = crH59*crH6;
        const double crH64 = N[1]*crH63 + crH25*crH58 + crH48*crH56;
        const double crH65 = N[2]*crH63 + crH29*crH58 + crH51*crH56;
        const double crH66 = DN_enr_vel_neg[1]*crH33;
        const double crH67 = crH66*gauss_weight;
        rH(0,0)+=crH20*(DN(0,0)*crH1 + DN_enr_p(0,0)*crH19);
        rH(0,1)+=crH20*(DN(0,1)*crH1 + DN_enr_p(0,1)*crH19);
        rH(0,2)+=crH21*(DN(0,0)*DN_enr_p(0,0) + DN(0,1)*DN_enr_p(0,1));
        rH(0,3)+=crH20*(DN(1,0)*crH1 + DN_enr_p(0,0)*crH26);
        rH(0,4)+=crH20*(DN(1,1)*crH1 + DN_enr_p(0,1)*crH26);
        rH(0,5)+=crH21*(DN(1,0)*DN_enr_p(0,0) + DN(1,1)*DN_enr_p(0,1));
        rH(0,6)+=crH20*(DN(2,0)*crH1 + DN_enr_p(0,0)*crH30);
        rH(0,7)+=crH20*(DN(2,1)*crH1 + DN_enr_p(0,1)*crH30);
        rH(0,8)+=crH21*(DN(2,0)*DN_enr_p(0,0) + DN(2,1)*DN_enr_p(0,1));
        rH(1,0)+=crH20*(DN(0,0)*crH31 + DN_enr_p(1,0)*crH19);
        rH(1,1)+=crH20*(DN(0,1)*crH31 + DN_enr_p(1,1)*crH19);
        rH(1,2)+=crH21*(DN(0,0)*DN_enr_p(1,0) + DN(0,1)*DN_enr_p(1,1));
        rH(1,3)+=crH20*(DN(1,0)*crH31 + DN_enr_p(1,0)*crH26);
        rH(1,4)+=crH20*(DN(1,1)*crH31 + DN_enr_p(1,1)*crH26);
        rH(1,5)+=crH21*(DN(1,0)*DN_enr_p(1,0) + DN(1,1)*DN_enr_p(1,1));
        rH(1,6)+=crH20*(DN(2,0)*crH31 + DN_enr_p(1,0)*crH30);
        rH(1,7)+=crH20*(DN(2,1)*crH31 + DN_enr_p(1,1)*crH30);
        rH(1,8)+=crH21*(DN(2,0)*DN_enr_p(1,0) + DN(2,1)*DN_enr_p(1,1));
        rH(2,0)+=crH20*(DN(0,0)*crH32 + DN_enr_p(2,0)*crH19);
        rH(2,1)+=crH20*(DN(0,1)*crH32 + DN_enr_p(2,1)*crH19);
        rH(2,2)+=crH21*(DN(0,0)*DN_enr_p(2,0) + DN(0,1)*DN_enr_p(2,1));
        rH(2,3)+=crH20*(DN(1,0)*crH32 + DN_enr_p(2,0)*crH26);
        rH(2,4)+=crH20*(DN(1,1)*crH32 + DN_enr_p(2,1)*crH26);
        rH(2,5)+=crH21*(DN(1,0)*DN_enr_p(2,0) + DN(1,1)*DN_enr_p(2,1));
        rH(2,6)+=crH20*(DN(2,0)*crH32 + DN_enr_p(2,0)*crH30);
        rH(2,7)+=crH20*(DN(2,1)*crH32 + DN_enr_p(2,1)*crH30);
        rH(2,8)+=crH21*(DN(2,0)*DN_enr_p(2,0) + DN(2,1)*DN_enr_p(2,1));
        rH(3,0)+=gauss_weight*(DN(0,0)*crH34 + crH44);
        rH(3,1)+=DN(0,1)*crH45;
        rH(3,2)+=gauss_weight*(DN(0,0)*crH47 - DN_enr_vel_pos[0]*N[0]);
        rH(3,3)+=gauss_weight*(DN(1,0)*crH34 + crH50);
        rH(3,4)+=DN(1,1)*crH45;
        rH(3,5)+=gauss_weight*(DN(1,0)*crH47 - DN_enr_vel_pos[0]*N[1]);
        rH(3,6)+=gauss_weight*(DN(2,0)*crH34 + crH52);
        rH(3,7)+=DN(2,1)*crH45;
        rH(3,8)+=gauss_weight*(DN(2,0)*crH47 - DN_enr_vel_pos[0]*N[2]);
        rH(4,0)+=DN(0,0)*crH54;
        rH(4,1)+=gauss_weight*(DN(0,1)*crH53 + crH44);
        rH(4,2)+=gauss_weight*(DN(0,1)*crH47 - DN_enr_vel_pos[1]*N[0]);
        rH(4,3)+=DN(1,0)*crH54;
        rH(4,4)+=gauss_weight*(DN(1,1)*crH53 + crH50);
        rH(4,5)+=gauss_weight*(DN(1,1)*crH47 - DN_enr_vel_pos[1]*N[1]);
        rH(4,6)+=DN(2,0)*crH54;
        rH(4,7)+=gauss_weight*(DN(2,1)*crH53 + crH52);
        rH(4,8)+=gauss_weight*(DN(2,1)*crH47 - DN_enr_vel_pos[1]*N[2]);
        rH(5,0)+=gauss_weight*(DN(0,0)*crH55 + crH60);
        rH(5,1)+=DN(0,1)*crH61;
        rH(5,2)+=gauss_weight*(DN(0,0)*crH62 - DN_enr_vel_neg[0]*N[0]);
        rH(5,3)+=gauss_weight*(DN(1,0)*crH55 + crH64);
        rH(5,4)+=DN(1,1)*crH61;
        rH(5,5)+=gauss_weight*(DN(1,0)*crH62 - DN_enr_vel_neg[0]*N[1]);
        rH(5,6)+=gauss_weight*(DN(2,0)*crH55 + crH65);
        rH(5,7)+=DN(2,1)*crH61;
        rH(5,8)+=gauss_weight*(DN(2,0)*crH62 - DN_enr_vel_neg[0]*N[2]);
        rH(6,0)+=DN(0,0)*crH67;
        rH(6,1)+=gauss_weight*(DN(0,1)*crH66 + crH60);
        rH(6,2)+=gauss_weight*(DN(0,1)*crH62 - DN_enr_vel_neg[1]*N[0]);
        rH(6,3)+=DN(1,0)*crH67;
        rH(6,4)+=gauss_weight*(DN(1,1)*crH66 + crH64);
        rH(6,5)+=gauss_weight*(DN(1,1)*crH62 - DN_enr_vel_neg[1]*N[1]);
        rH(6,6)+=DN(2,0)*crH67;
        rH(6,7)+=gauss_weight*(DN(2,1)*crH66 + crH65);
        rH(6,8)+=gauss_weight*(DN(2,1)*crH62 - DN_enr_vel_neg[1]*N[2]);


            const double crKee0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crKee1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crKee2 = rho*stab_c2*sqrt(pow(crKee0, 2) + pow(crKee1, 2));
        const double crKee3 = 1.0/(crKee2/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
        const double crKee4 = crKee3*gauss_weight;
        const double crKee5 = crKee4*(DN_enr_p(0,0)*DN_enr_p(1,0) + DN_enr_p(0,1)*DN_enr_p(1,1));
        const double crKee6 = crKee4*(DN_enr_p(0,0)*DN_enr_p(2,0) + DN_enr_p(0,1)*DN_enr_p(2,1));
        const double crKee7 = DN_enr_vel_pos[0]*N_enr_p[0];
        const double crKee8 = DN_enr_vel_pos[0]*crKee0 + DN_enr_vel_pos[1]*crKee1;
        const double crKee9 = crKee3*rho;
        const double crKee10 = crKee8*crKee9;
        const double crKee11 = DN_enr_vel_pos[1]*N_enr_p[0];
        const double crKee12 = DN_enr_vel_neg[0]*N_enr_p[0];
        const double crKee13 = DN_enr_vel_neg[0]*crKee0 + DN_enr_vel_neg[1]*crKee1;
        const double crKee14 = crKee13*crKee9;
        const double crKee15 = DN_enr_vel_neg[1]*N_enr_p[0];
        const double crKee16 = crKee4*(DN_enr_p(1,0)*DN_enr_p(2,0) + DN_enr_p(1,1)*DN_enr_p(2,1));
        const double crKee17 = DN_enr_vel_pos[0]*N_enr_p[1];
        const double crKee18 = DN_enr_vel_pos[1]*N_enr_p[1];
        const double crKee19 = DN_enr_vel_neg[0]*N_enr_p[1];
        const double crKee20 = DN_enr_vel_neg[1]*N_enr_p[1];
        const double crKee21 = DN_enr_vel_pos[0]*N_enr_p[2];
        const double crKee22 = DN_enr_vel_pos[1]*N_enr_p[2];
        const double crKee23 = DN_enr_vel_neg[0]*N_enr_p[2];
        const double crKee24 = DN_enr_vel_neg[1]*N_enr_p[2];
        const double crKee25 = N_enr_vel_pos*rho;
        const double crKee26 = crKee3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
        const double crKee27 = crKee25*crKee26;
        const double crKee28 = crKee2*h/stab_c1 + mu;
        const double crKee29 = crKee26*pow(rho, 2);
        const double crKee30 = N_enr_vel_pos*crKee29;
        const double crKee31 = crKee25*crKee8 + crKee30*crKee8;
        const double crKee32 = DN_enr_vel_pos[0]*crKee28;
        const double crKee33 = crKee32*gauss_weight;
        const double crKee34 = DN_enr_vel_pos[1]*crKee33;
        const double crKee35 = DN_enr_vel_neg[0]*crKee32;
        const double crKee36 = crKee13*crKee25 + crKee13*crKee30;
        const double crKee37 = DN_enr_vel_neg[1]*crKee33;
        const double crKee38 = DN_enr_vel_pos[1]*crKee28;
        const double crKee39 = DN_enr_vel_neg[0]*gauss_weight;
        const double crKee40 = crKee38*crKee39;
        const double crKee41 = DN_enr_vel_neg[1]*crKee38;
        const double crKee42 = N_enr_vel_neg*rho;
        const double crKee43 = crKee26*crKee42;
        const double crKee44 = N_enr_vel_neg*crKee29;
        const double crKee45 = crKee42*crKee8 + crKee44*crKee8;
        const double crKee46 = crKee13*crKee42 + crKee13*crKee44;
        const double crKee47 = DN_enr_vel_neg[1]*crKee28*crKee39;
        rKee(0,0)+=crKee4*(pow(DN_enr_p(0,0), 2) + pow(DN_enr_p(0,1), 2));
        rKee(0,1)+=crKee5;
        rKee(0,2)+=crKee6;
        rKee(0,3)+=gauss_weight*(DN_enr_p(0,0)*crKee10 + crKee7);
        rKee(0,4)+=gauss_weight*(DN_enr_p(0,1)*crKee10 + crKee11);
        rKee(0,5)+=gauss_weight*(DN_enr_p(0,0)*crKee14 + crKee12);
        rKee(0,6)+=gauss_weight*(DN_enr_p(0,1)*crKee14 + crKee15);
        rKee(1,0)+=crKee5;
        rKee(1,1)+=crKee4*(pow(DN_enr_p(1,0), 2) + pow(DN_enr_p(1,1), 2));
        rKee(1,2)+=crKee16;
        rKee(1,3)+=gauss_weight*(DN_enr_p(1,0)*crKee10 + crKee17);
        rKee(1,4)+=gauss_weight*(DN_enr_p(1,1)*crKee10 + crKee18);
        rKee(1,5)+=gauss_weight*(DN_enr_p(1,0)*crKee14 + crKee19);
        rKee(1,6)+=gauss_weight*(DN_enr_p(1,1)*crKee14 + crKee20);
        rKee(2,0)+=crKee6;
        rKee(2,1)+=crKee16;
        rKee(2,2)+=crKee4*(pow(DN_enr_p(2,0), 2) + pow(DN_enr_p(2,1), 2));
        rKee(2,3)+=gauss_weight*(DN_enr_p(2,0)*crKee10 + crKee21);
        rKee(2,4)+=gauss_weight*(DN_enr_p(2,1)*crKee10 + crKee22);
        rKee(2,5)+=gauss_weight*(DN_enr_p(2,0)*crKee14 + crKee23);
        rKee(2,6)+=gauss_weight*(DN_enr_p(2,1)*crKee14 + crKee24);
        rKee(3,0)+=gauss_weight*(DN_enr_p(0,0)*crKee27 - crKee7);
        rKee(3,1)+=gauss_weight*(DN_enr_p(1,0)*crKee27 - crKee17);
        rKee(3,2)+=gauss_weight*(DN_enr_p(2,0)*crKee27 - crKee21);
        rKee(3,3)+=gauss_weight*(pow(DN_enr_vel_pos[0], 2)*crKee28 + crKee31);
        rKee(3,4)+=crKee34;
        rKee(3,5)+=gauss_weight*(crKee35 + crKee36);
        rKee(3,6)+=crKee37;
        rKee(4,0)+=gauss_weight*(DN_enr_p(0,1)*crKee27 - crKee11);
        rKee(4,1)+=gauss_weight*(DN_enr_p(1,1)*crKee27 - crKee18);
        rKee(4,2)+=gauss_weight*(DN_enr_p(2,1)*crKee27 - crKee22);
        rKee(4,3)+=crKee34;
        rKee(4,4)+=gauss_weight*(pow(DN_enr_vel_pos[1], 2)*crKee28 + crKee31);
        rKee(4,5)+=crKee40;
        rKee(4,6)+=gauss_weight*(crKee36 + crKee41);
        rKee(5,0)+=gauss_weight*(DN_enr_p(0,0)*crKee43 - crKee12);
        rKee(5,1)+=gauss_weight*(DN_enr_p(1,0)*crKee43 - crKee19);
        rKee(5,2)+=gauss_weight*(DN_enr_p(2,0)*crKee43 - crKee23);
        rKee(5,3)+=gauss_weight*(crKee35 + crKee45);
        rKee(5,4)+=crKee40;
        rKee(5,5)+=gauss_weight*(pow(DN_enr_vel_neg[0], 2)*crKee28 + crKee46);
        rKee(5,6)+=crKee47;
        rKee(6,0)+=gauss_weight*(DN_enr_p(0,1)*crKee43 - crKee15);
        rKee(6,1)+=gauss_weight*(DN_enr_p(1,1)*crKee43 - crKee20);
        rKee(6,2)+=gauss_weight*(DN_enr_p(2,1)*crKee43 - crKee24);
        rKee(6,3)+=crKee37;
        rKee(6,4)+=gauss_weight*(crKee41 + crKee45);
        rKee(6,5)+=crKee47;
        rKee(6,6)+=gauss_weight*(pow(DN_enr_vel_neg[1], 2)*crKee28 + crKee46);


            const double crRHSee0 = 1.0*v(0,0);
        const double crRHSee1 = 1.0*v(0,1);
        const double crRHSee2 = 1.0*v(1,0);
        const double crRHSee3 = 1.0*v(1,1);
        const double crRHSee4 = 1.0*v(2,0);
        const double crRHSee5 = 1.0*v(2,1);
        const double crRHSee6 = 1.0*vn(0,0);
        const double crRHSee7 = 1.0*vn(0,1);
        const double crRHSee8 = 1.0*vn(1,0);
        const double crRHSee9 = 1.0*vn(1,1);
        const double crRHSee10 = 1.0*vn(2,0);
        const double crRHSee11 = 1.0*vn(2,1);
        const double crRHSee12 = DN_enr_vel_neg[0]*v_enr_neg[0] + DN_enr_vel_pos[0]*v_enr_pos[0];
        const double crRHSee13 = DN_enr_vel_neg[1]*v_enr_neg[1] + DN_enr_vel_pos[1]*v_enr_pos[1];
        const double crRHSee14 = crRHSee12 + crRHSee13 - volume_error_ratio;
        const double crRHSee15 = crRHSee14 + (DN(0,0)*crRHSee0 + DN(0,0)*crRHSee6*max_spectral_radius + DN(0,1)*crRHSee1 + DN(0,1)*crRHSee7*max_spectral_radius + DN(1,0)*crRHSee2 + DN(1,0)*crRHSee8*max_spectral_radius + DN(1,1)*crRHSee3 + DN(1,1)*crRHSee9*max_spectral_radius + DN(2,0)*crRHSee10*max_spectral_radius + DN(2,0)*crRHSee4 + DN(2,1)*crRHSee11*max_spectral_radius + DN(2,1)*crRHSee5)/(max_spectral_radius + 1);
        const double crRHSee16 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crRHSee17 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crRHSee18 = rho*(crRHSee12*crRHSee16 + crRHSee17*(DN_enr_vel_neg[1]*v_enr_neg[0] + DN_enr_vel_pos[1]*v_enr_pos[0]));
        const double crRHSee19 = 1.0/(max_spectral_radius + 1.0);
        const double crRHSee20 = rho*(N[0]*(crRHSee19*(1.0*f(0,0) - 1.0*fn(0,0)) + fn(0,0)) + N[1]*(crRHSee19*(1.0*f(1,0) - 1.0*fn(1,0)) + fn(1,0)) + N[2]*(crRHSee19*(1.0*f(2,0) - 1.0*fn(2,0)) + fn(2,0)));
        const double crRHSee21 = crRHSee19*(crRHSee0 - crRHSee6) + vn(0,0);
        const double crRHSee22 = crRHSee19*(crRHSee2 - crRHSee8) + vn(1,0);
        const double crRHSee23 = crRHSee19*(-crRHSee10 + crRHSee4) + vn(2,0);
        const double crRHSee24 = rho*(crRHSee16*(DN(0,0)*crRHSee21 + DN(1,0)*crRHSee22 + DN(2,0)*crRHSee23) + crRHSee17*(DN(0,1)*crRHSee21 + DN(1,1)*crRHSee22 + DN(2,1)*crRHSee23));
        const double crRHSee25 = v(0,0) - vn(0,0);
        const double crRHSee26 = 1.0/dt;
        const double crRHSee27 = 1.0*crRHSee19;
        const double crRHSee28 = 0.5*max_spectral_radius - 1.5;
        const double crRHSee29 = -crRHSee28;
        const double crRHSee30 = 1.0/(crRHSee19*crRHSee29 - crRHSee27 + 0.5);
        const double crRHSee31 = crRHSee26*crRHSee30;
        const double crRHSee32 = -crRHSee19*crRHSee29 + crRHSee27 + 0.5;
        const double crRHSee33 = -crRHSee30*crRHSee32;
        const double crRHSee34 = max_spectral_radius - 3.0;
        const double crRHSee35 = 0.5*crRHSee19;
        const double crRHSee36 = -crRHSee34*crRHSee35;
        const double crRHSee37 = v(1,0) - vn(1,0);
        const double crRHSee38 = v(2,0) - vn(2,0);
        const double crRHSee39 = rho*stab_c2*sqrt(pow(crRHSee16, 2) + pow(crRHSee17, 2));
        const double crRHSee40 = 1.0/(crRHSee26*dyn_tau*rho + crRHSee39/h + mu*stab_c1/pow(h, 2));
        const double crRHSee41 = crRHSee40*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN_enr_p(0,0)*penr[0] + DN_enr_p(1,0)*penr[1] + DN_enr_p(2,0)*penr[2] + crRHSee18 - crRHSee20 + crRHSee24 + rho*(N[0]*(acceleration_alpha_method(0,0) + crRHSee36*(acceleration_alpha_method(0,0)*crRHSee33 - acceleration_alpha_method(0,0) + crRHSee25*crRHSee31)) + N[1]*(acceleration_alpha_method(1,0) + crRHSee36*(acceleration_alpha_method(1,0)*crRHSee33 - acceleration_alpha_method(1,0) + crRHSee31*crRHSee37)) + N[2]*(acceleration_alpha_method(2,0) + crRHSee36*(acceleration_alpha_method(2,0)*crRHSee33 - acceleration_alpha_method(2,0) + crRHSee31*crRHSee38))));
        const double crRHSee42 = rho*(crRHSee13*crRHSee17 + crRHSee16*(DN_enr_vel_neg[0]*v_enr_neg[1] + DN_enr_vel_pos[0]*v_enr_pos[1]));
        const double crRHSee43 = rho*(N[0]*(crRHSee19*(1.0*f(0,1) - 1.0*fn(0,1)) + fn(0,1)) + N[1]*(crRHSee19*(1.0*f(1,1) - 1.0*fn(1,1)) + fn(1,1)) + N[2]*(crRHSee19*(1.0*f(2,1) - 1.0*fn(2,1)) + fn(2,1)));
        const double crRHSee44 = crRHSee19*(crRHSee1 - crRHSee7) + vn(0,1);
        const double crRHSee45 = crRHSee19*(crRHSee3 - crRHSee9) + vn(1,1);
        const double crRHSee46 = crRHSee19*(-crRHSee11 + crRHSee5) + vn(2,1);
        const double crRHSee47 = rho*(crRHSee16*(DN(0,0)*crRHSee44 + DN(1,0)*crRHSee45 + DN(2,0)*crRHSee46) + crRHSee17*(DN(0,1)*crRHSee44 + DN(1,1)*crRHSee45 + DN(2,1)*crRHSee46));
        const double crRHSee48 = v(0,1) - vn(0,1);
        const double crRHSee49 = v(1,1) - vn(1,1);
        const double crRHSee50 = v(2,1) - vn(2,1);
        const double crRHSee51 = crRHSee40*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN_enr_p(0,1)*penr[0] + DN_enr_p(1,1)*penr[1] + DN_enr_p(2,1)*penr[2] + crRHSee42 - crRHSee43 + crRHSee47 + rho*(N[0]*(acceleration_alpha_method(0,1) + crRHSee36*(acceleration_alpha_method(0,1)*crRHSee33 - acceleration_alpha_method(0,1) + crRHSee31*crRHSee48)) + N[1]*(acceleration_alpha_method(1,1) + crRHSee36*(acceleration_alpha_method(1,1)*crRHSee33 - acceleration_alpha_method(1,1) + crRHSee31*crRHSee49)) + N[2]*(acceleration_alpha_method(2,1) + crRHSee36*(acceleration_alpha_method(2,1)*crRHSee33 - acceleration_alpha_method(2,1) + crRHSee31*crRHSee50))));
        const double crRHSee52 = 1.0*DN_enr_vel_pos[0];
        const double crRHSee53 = 1.0*DN_enr_vel_pos[1];
        const double crRHSee54 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
        const double crRHSee55 = N_enr_p[0]*penr[0] + N_enr_p[1]*penr[1] + N_enr_p[2]*penr[2];
        const double crRHSee56 = (crRHSee39*h/stab_c1 + mu)*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + crRHSee14);
        const double crRHSee57 = 1.0/(-crRHSee19*crRHSee28 - crRHSee27 + 0.5);
        const double crRHSee58 = crRHSee32*crRHSee57;
        const double crRHSee59 = crRHSee34*crRHSee35;
        const double crRHSee60 = N[0]*(acceleration_alpha_method(0,0) - crRHSee59*(-acceleration_alpha_method(0,0)*crRHSee58 - acceleration_alpha_method(0,0) + crRHSee25*crRHSee26*crRHSee57)) + N[1]*(acceleration_alpha_method(1,0) - crRHSee59*(-acceleration_alpha_method(1,0)*crRHSee58 - acceleration_alpha_method(1,0) + crRHSee26*crRHSee37*crRHSee57)) + N[2]*(acceleration_alpha_method(2,0) - crRHSee59*(-acceleration_alpha_method(2,0)*crRHSee58 - acceleration_alpha_method(2,0) + crRHSee26*crRHSee38*crRHSee57));
        const double crRHSee61 = N_enr_vel_pos*rho;
        const double crRHSee62 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
        const double crRHSee63 = crRHSee61*crRHSee62;
        const double crRHSee64 = N[0]*(acceleration_alpha_method(0,1) - crRHSee59*(-acceleration_alpha_method(0,1)*crRHSee58 - acceleration_alpha_method(0,1) + crRHSee26*crRHSee48*crRHSee57)) + N[1]*(acceleration_alpha_method(1,1) - crRHSee59*(-acceleration_alpha_method(1,1)*crRHSee58 - acceleration_alpha_method(1,1) + crRHSee26*crRHSee49*crRHSee57)) + N[2]*(acceleration_alpha_method(2,1) - crRHSee59*(-acceleration_alpha_method(2,1)*crRHSee58 - acceleration_alpha_method(2,1) + crRHSee26*crRHSee50*crRHSee57));
        const double crRHSee65 = 1.0*DN_enr_vel_neg[0];
        const double crRHSee66 = 1.0*DN_enr_vel_neg[1];
        const double crRHSee67 = N_enr_vel_neg*rho;
        const double crRHSee68 = crRHSee62*crRHSee67;
        rRHSee[0]+=-gauss_weight*(DN_enr_p(0,0)*crRHSee41 + DN_enr_p(0,1)*crRHSee51 + N_enr_p[0]*crRHSee15);
        rRHSee[1]+=-gauss_weight*(DN_enr_p(1,0)*crRHSee41 + DN_enr_p(1,1)*crRHSee51 + N_enr_p[1]*crRHSee15);
        rRHSee[2]+=-gauss_weight*(DN_enr_p(2,0)*crRHSee41 + DN_enr_p(2,1)*crRHSee51 + N_enr_p[2]*crRHSee15);
        rRHSee[3]+=-gauss_weight*(-DN_enr_vel_pos[0]*crRHSee54 - DN_enr_vel_pos[0]*crRHSee55 + DN_enr_vel_pos[0]*crRHSee56 + N_enr_vel_pos*crRHSee18 - N_enr_vel_pos*crRHSee20 + N_enr_vel_pos*crRHSee24 + crRHSee41*crRHSee63 + crRHSee52*stress[0] + crRHSee53*stress[2] + crRHSee60*crRHSee61);
        rRHSee[4]+=-gauss_weight*(-DN_enr_vel_pos[1]*crRHSee54 - DN_enr_vel_pos[1]*crRHSee55 + DN_enr_vel_pos[1]*crRHSee56 + N_enr_vel_pos*crRHSee42 - N_enr_vel_pos*crRHSee43 + N_enr_vel_pos*crRHSee47 + crRHSee51*crRHSee63 + crRHSee52*stress[2] + crRHSee53*stress[1] + crRHSee61*crRHSee64);
        rRHSee[5]+=-gauss_weight*(-DN_enr_vel_neg[0]*crRHSee54 - DN_enr_vel_neg[0]*crRHSee55 + DN_enr_vel_neg[0]*crRHSee56 + N_enr_vel_neg*crRHSee18 - N_enr_vel_neg*crRHSee20 + N_enr_vel_neg*crRHSee24 + crRHSee41*crRHSee68 + crRHSee60*crRHSee67 + crRHSee65*stress[0] + crRHSee66*stress[2]);
        rRHSee[6]+=-gauss_weight*(-DN_enr_vel_neg[1]*crRHSee54 - DN_enr_vel_neg[1]*crRHSee55 + DN_enr_vel_neg[1]*crRHSee56 + N_enr_vel_neg*crRHSee42 - N_enr_vel_neg*crRHSee43 + N_enr_vel_neg*crRHSee47 + crRHSee51*crRHSee68 + crRHSee64*crRHSee67 + crRHSee65*stress[2] + crRHSee66*stress[1]);

}

template <>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodDiscontinuousData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHSee)
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

    // Get material response data
    const Matrix &C = rData.C;
    const auto &stress = rData.ShearStress;

    // Get standard shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Get pressure enrichment shape function values
    const auto &N_enr_p = rData.Nenr;
    const auto &DN_enr_p = rData.DN_DXenr;

    // Get velocity enrichment bubble function values
    const auto &N_enr_vel = rData.N_enr_vel;
    const auto &DN_enr_vel = rData.DN_DX_enr_vel;

    const bool is_positive = inner_prod(N, rData.Distance) < 0.0 ? false : true;
    const double N_enr_vel_pos = is_positive ? N_enr_vel : 0.0;
    const double N_enr_vel_neg = is_positive ? 0.0 : N_enr_vel;
    const array_1d<double,3>& DN_enr_vel_pos = is_positive ? DN_enr_vel : ZeroVector(3);
    const array_1d<double,3>& DN_enr_vel_neg = is_positive ? ZeroVector(3) : DN_enr_vel;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_ratio = rData.VolumeErrorRate;

    // Initialize enrichment DOFs appearing in the enrichment RHS to zero
    // Note that we always initialize them to zero as we do not want to store them
    // array_1d<double, Dim> v_enr = ZeroVector(Dim);
    array_1d<double, 3> v_enr_pos = ZeroVector(3);
    array_1d<double, 3> v_enr_neg = ZeroVector(3);
    array_1d<double, 4> penr = ZeroVector(4);

    // Add current Gauss point enrichment contribution
    const double gauss_weight = rData.Weight;

    //substitute_enrichment_V_3D

    //substitute_enrichment_H_3D

    //substitute_enrichment_Kee_3D

    //substitute_enrichment_rhs_ee_3D
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::ComputeSplitting(
    TElementData &rData,
    MatrixType &rStandardShapeFunctionsPos,
    MatrixType &rStandardShapeFunctionsNeg,
    MatrixType &rEnrichedShapeFunctionsPos,
    MatrixType &rEnrichedShapeFunctionsNeg,
    VectorType &rBubbleShapeFunctionsPos,
    VectorType &rBubbleShapeFunctionsNeg,
    GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg,
    DenseVector<array_1d<double,Dim>> & rBubbleShapeDerivativesPos,
    DenseVector<array_1d<double,Dim>> & rBubbleShapeDerivativesNeg)
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
        } else {
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Set the standard modified shape functions pointer
    auto p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::UniquePointer p_stand_mod_sh_funcs;
    if constexpr (TElementData::Dim == 2 && TElementData::NumNodes == 3) {
        p_stand_mod_sh_funcs = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    } else if (TElementData::Dim == 3 && TElementData::NumNodes == 4) {
        p_stand_mod_sh_funcs = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);
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

    // Compute the bubble shape function values
    double dist;
    double abs_dist;
    array_1d<double,Dim> abs_dist_grad;
    const double bubble_scaling = 1.0;

    const SizeType n_pos_gauss = rStandardShapeFunctionsPos.size1();
    rBubbleShapeFunctionsPos.resize(n_pos_gauss);
    rBubbleShapeDerivativesPos.resize(n_pos_gauss);
    for (IndexType g_pos = 0; g_pos < n_pos_gauss; ++g_pos) {
        // Get standard shape functions values at current positive Gauss point
        const auto& r_N_g_pos = row(rStandardShapeFunctionsPos, g_pos);
        const auto& r_DN_g_pos = rStandardShapeDerivativesPos[g_pos];

        // Compute distance interpolations
        dist = 0.0;
        noalias(abs_dist_grad) = ZeroVector(Dim);
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            dist += r_N_g_pos(i_node) * rData.Distance[i_node];
            for (IndexType d = 0; d < Dim; ++d) {
                abs_dist_grad[d] += r_DN_g_pos(i_node, d) * rData.Distance[i_node];
            }
        }
        abs_dist = std::abs(dist);
        abs_dist_grad *= dist > 0.0 ? 1.0 : -1.0;

        // Initialize and calculate the current positive Gauss point bubble function values
        if constexpr (Dim == 2) {
            rBubbleShapeFunctionsPos(g_pos) = bubble_scaling*r_N_g_pos(0)*r_N_g_pos(1)*r_N_g_pos(2)*abs_dist;
            for (IndexType d = 0; d < Dim; ++d) {
                rBubbleShapeDerivativesPos(g_pos)[d] = bubble_scaling * (
                    r_DN_g_pos(0, d)*r_N_g_pos(1)*r_N_g_pos(2)*abs_dist +
                    r_N_g_pos(0)*r_DN_g_pos(1, d)*r_N_g_pos(2)*abs_dist +
                    r_N_g_pos(0)*r_N_g_pos(1)*r_DN_g_pos(2, d)*abs_dist +
                    r_N_g_pos(0)*r_N_g_pos(1)*r_N_g_pos(2)*abs_dist_grad[d]);
            }
        } else {
            rBubbleShapeFunctionsPos(g_pos) = bubble_scaling*r_N_g_pos(0)*r_N_g_pos(1)*r_N_g_pos(2)*r_N_g_pos(3)*abs_dist;
            for (IndexType d = 0; d < Dim; ++d) {
                rBubbleShapeDerivativesPos(g_pos)[d] = bubble_scaling * (
                    r_DN_g_pos(0, d)*r_N_g_pos(1)*r_N_g_pos(2)*r_N_g_pos(3)*abs_dist +
                    r_N_g_pos(0)*r_DN_g_pos(1, d)*r_N_g_pos(2)*r_N_g_pos(3)*abs_dist +
                    r_N_g_pos(0)*r_N_g_pos(1)*r_DN_g_pos(2, d)*r_N_g_pos(3)*abs_dist +
                    r_N_g_pos(0)*r_N_g_pos(1)*r_N_g_pos(2)*r_DN_g_pos(3, d)*abs_dist +
                    r_N_g_pos(0)*r_N_g_pos(1)*r_N_g_pos(2)*r_N_g_pos(3)*abs_dist_grad[d]);
            }
        }
    }

    const SizeType n_neg_gauss = rStandardShapeFunctionsNeg.size1();
    rBubbleShapeFunctionsNeg.resize(n_neg_gauss);
    rBubbleShapeDerivativesNeg.resize(n_neg_gauss);
    for (IndexType g_neg = 0; g_neg < n_neg_gauss; ++g_neg) {
        // Get standard shape functions values at current negative Gauss point
        const auto& r_N_g_neg = row(rStandardShapeFunctionsNeg, g_neg);
        const auto& r_DN_g_neg = rStandardShapeDerivativesNeg[g_neg];

        // Compute distance interpolations
        dist = 0.0;
        noalias(abs_dist_grad) = ZeroVector(Dim);
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            dist += r_N_g_neg(i_node) * rData.Distance[i_node];
            abs_dist += r_N_g_neg(i_node) * std::abs(rData.Distance[i_node]);
            for (IndexType d = 0; d < Dim; ++d) {
                abs_dist_grad[d] += r_DN_g_neg(i_node, d) * rData.Distance[i_node];
            }
        }
        abs_dist = std::abs(dist);
        abs_dist_grad *= dist > 0.0 ? 1.0 : -1.0;

        // Initialize and calculate the current positive Gauss point bubble function values
        if constexpr (Dim == 2) {
            rBubbleShapeFunctionsNeg(g_neg) = bubble_scaling*r_N_g_neg(0)*r_N_g_neg(1)*r_N_g_neg(2)*abs_dist;
            for (IndexType d = 0; d < Dim; ++d) {
                rBubbleShapeDerivativesNeg(g_neg)[d] = bubble_scaling * (
                    r_DN_g_neg(0, d)*r_N_g_neg(1)*r_N_g_neg(2)*abs_dist +
                    r_N_g_neg(0)*r_DN_g_neg(1, d)*r_N_g_neg(2)*abs_dist +
                    r_N_g_neg(0)*r_N_g_neg(1)*r_DN_g_neg(2, d)*abs_dist +
                    r_N_g_neg(0)*r_N_g_neg(1)*r_N_g_neg(2)*abs_dist_grad[d]);
            }
        } else {
            rBubbleShapeFunctionsNeg(g_neg) = bubble_scaling*r_N_g_neg(0)*r_N_g_neg(1)*r_N_g_neg(2)*r_N_g_neg(3)*abs_dist;
            for (IndexType d = 0; d < Dim; ++d) {
                rBubbleShapeDerivativesNeg(g_neg)[d] = bubble_scaling * (
                    r_DN_g_neg(0, d)*r_N_g_neg(1)*r_N_g_neg(2)*r_N_g_neg(3)*abs_dist +
                    r_N_g_neg(0)*r_DN_g_neg(1, d)*r_N_g_neg(2)*r_N_g_neg(3)*abs_dist +
                    r_N_g_neg(0)*r_N_g_neg(1)*r_DN_g_neg(2, d)*r_N_g_neg(3)*abs_dist +
                    r_N_g_neg(0)*r_N_g_neg(1)*r_N_g_neg(2)*r_DN_g_neg(3, d)*abs_dist +
                    r_N_g_neg(0)*r_N_g_neg(1)*r_N_g_neg(2)*r_N_g_neg(3)*abs_dist_grad[d]);
            }
        }
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
    const typename TElementData::MatrixRowType& rEnrichedShapeFunctions,
    const typename TElementData::ShapeDerivativesType& rEnrichedShapeFunctionsGradients,
    const double rBubbleShapeFunction,
    const array_1d<double, Dim>& rBubbleShapeFunctionGradients) const
{
    // Update Gauss point current values
    rData.UpdateGeometryValues(
        IntegrationPointIndex,
        IntegrationPointWeight,
        rStandardShapeFunctions,
        rStandardShapeFunctionsGradients,
        rEnrichedShapeFunctions,
        rEnrichedShapeFunctionsGradients,
        rBubbleShapeFunction,
        rBubbleShapeFunctionGradients);

    // Calculate material response
    const double d_gauss = inner_prod(rData.Distance, rStandardShapeFunctions);
    if (d_gauss > 0.0) {
        rData.CalculateAirMaterialResponse();
    } else {
        this->CalculateMaterialResponse(rData);
    }
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>::CondenseAndSaveEnrichmentWithContinuity(
    const TElementData& rData,
    Matrix& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const MatrixType& rHTot,
    const MatrixType& rVTot,
    MatrixType& rKeeTot,
    const VectorType& rRHSeeTot)
{
    // Compute positive side, negative side and total volumes
    double pos_vol = 0.0;
    double neg_vol = 0.0;
    for (const double w_gauss_pos : rData.w_gauss_pos_side) {
        pos_vol += w_gauss_pos;
    }
    for (const double w_gauss_neg : rData.w_gauss_neg_side) {
        neg_vol += w_gauss_neg;
    }
    const double tot_vol = pos_vol + neg_vol;

    //We only enrich elements which are not almost empty/full
    const double min_area_ratio = 1e-7;
    if (pos_vol / tot_vol > min_area_ratio && neg_vol / tot_vol > min_area_ratio) {

        // Compute the maximum diagonal value in the pressure enrichment stiffness matrix
        double max_diag = 0.0;
        for (IndexType k = 0; k < NumNodes; ++k){
            if (std::abs(rKeeTot(k, k)) > max_diag){
                max_diag = std::abs(rKeeTot(k, k));
            }
        }
        if (max_diag < 1.0e-12){
            max_diag = 1.0;
        }

        // "Weakly" impose continuity of pressure enrichment in the interserced edges
        // Note that this effectively leaves us with a unique pressure enrichment DOF
        for (IndexType i = 0; i < Dim; ++i){
            const double abs_d_i = std::abs(rData.Distance[i]);
            for (IndexType j = i + 1; j < NumNodes; ++j){
                const double abs_d_j = std::abs(rData.Distance[j]);
                // Check if the edge is cut, if it is, set the penalty constraint
                if (rData.Distance[i] * rData.Distance[j] < 0.0){
                    const double sum_d = abs_d_i + abs_d_j;
                    const double N_i = abs_d_j / sum_d;
                    const double N_j = abs_d_i / sum_d;
                    const double penalty_coeff = max_diag * 0.001; // h/BDFVector[0];
                    rKeeTot(i, i) += penalty_coeff * N_i * N_i;
                    rKeeTot(i, j) -= penalty_coeff * N_i * N_j;
                    rKeeTot(j, i) -= penalty_coeff * N_j * N_i;
                    rKeeTot(j, j) += penalty_coeff * N_j * N_j;
                }
            }
        }

        // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        double det;
        BoundedMatrix<double, NumNodes+2*Dim, NumNodes+2*Dim> KeeTotInv;
        MathUtils<double>::InvertMatrix(rKeeTot, KeeTotInv, det);

        const BoundedMatrix<double, NumNodes+2*Dim, LocalSize> tmp = prod(KeeTotInv, rHTot);
        noalias(rLeftHandSideMatrix) -= prod(rVTot, tmp);

        const array_1d<double, NumNodes+2*Dim> tmp2 = prod(KeeTotInv, rRHSeeTot);
        noalias(rRightHandSideVector) -= prod(rVTot, tmp2);

        // Calculate and save the enrichment contributions so we can use them later on
        // These are used for the calculation of the strain for the constitutive law or the postprocess
        const double rho_inf = rData.MaxSpectralRadius;
        const double alpha_f = 1.0 / (rho_inf+1.0);
        array_1d<double, LocalSize> sol_vect;
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            const auto& r_v_i = row(rData.Velocity, i_node);
            const auto& r_v_n_i = row(rData.Velocity_OldStep1, i_node);
            for (IndexType d = 0; d < Dim; ++d) {
                sol_vect[i_node*BlockSize + d] = r_v_n_i[d] + alpha_f*(r_v_i[d]-r_v_n_i[d]);
            }
            sol_vect[i_node*BlockSize + Dim] = rData.Pressure[i_node];
        }

        const Vector H_sol_vect = prod(rHTot, sol_vect);
        const Vector aux_enr_RHS = rRHSeeTot - H_sol_vect;
        const Vector enr_vect = prod(KeeTotInv, aux_enr_RHS);

        // Save the enrichment contributions
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            this->mPresEnr[i_node] = enr_vect[i_node];
        }
        for (IndexType d = 0; d < Dim; ++d) {
            this->mVelEnrPos[d] = enr_vect[NumNodes + d];
            this->mVelEnrNeg[d] = enr_vect[NumNodes + Dim + d];
        }
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
