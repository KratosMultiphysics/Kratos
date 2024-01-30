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
            rOutput[g] = this->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        }
    }
    else{
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
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
    rOutput = 0.0; 
    if (rVariable == ARTIFICIAL_DYNAMIC_VISCOSITY)
    {

        // Iterate over integration points to evaluate the artificial viscosity at each Gauss point
        for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        {
            this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            rOutput += CalculateArtificialDynamicViscositySpecialization(data);
        }

        rOutput /= number_of_gauss_points;
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

    //substitute_artificial_mu_grad_v_norm_2D_3N

    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        //substitute_artificial_mu_2D_3N
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

    //substitute_artificial_mu_grad_v_norm_3D_4N

    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        //substitute_artificial_mu_3D_4N
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

    //substitute_lhs_2D

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

    //substitute_lhs_3D

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

    //substitute_rhs_2D

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

    //substitute_rhs_3D

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

    //substitute_enrichment_V_2D

    //substitute_enrichment_H_2D

    //substitute_enrichment_Kee_2D

    //substitute_enrichment_rhs_ee_2D

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
