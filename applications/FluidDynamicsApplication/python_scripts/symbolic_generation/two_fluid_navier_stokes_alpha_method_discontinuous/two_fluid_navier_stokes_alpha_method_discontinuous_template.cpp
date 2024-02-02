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
    //substitute_lhs_2D
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
    //substitute_rhs_2D
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

    //substitute_enrichment_V_2D

    //substitute_enrichment_H_2D

    //substitute_enrichment_Kee_2D

    //substitute_enrichment_rhs_ee_2D
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
