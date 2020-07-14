//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Mohit Tyagi
//  Co-authors:      Daniel Diez, Pablo Becker, Ruben Zorrilla
//

#include "symbolic_stokes.h"
#include "custom_utilities/symbolic_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
SymbolicStokes<TElementData>::~SymbolicStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
std::string SymbolicStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "SymbolicStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void SymbolicStokes<TElementData>::PrintInfo(
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
void SymbolicStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <>
void SymbolicStokes<SymbolicStokesData<2,3>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<2,3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_lhs_2D3

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,4>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<2,4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_lhs_2D4

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,4> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_lhs_3D4

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,6>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,6> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_lhs_3D6

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,8>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,8> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_lhs_3D8

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,3>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<2,3> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_rhs_2D3

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,4>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<2,4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }
    //substitute_rhs_2D4

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_rhs_3D4

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,6>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,6> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_rhs_3D6

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,8>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,8> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    //constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;
    double dt_inv = 0.0;
    if (dt > 1e-09)
    {
        dt_inv = 1.0/dt;
    }
    if (std::abs(bdf0) < 1e-9)
    {
        dt_inv = 0.0;
    }

    //substitute_rhs_3D8

    noalias(rRHS) += rData.Weight * rhs;
}


template <class TElementData>
void SymbolicStokes<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void SymbolicStokes<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void SymbolicStokes<TElementData>::GetValueOnIntegrationPoints(
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

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; i_gauss++){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; nnode++){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ndim++){
                    DVi_DXi += gp_DN_DX(nnode, ndim) * vel[ndim];
                }
            }
            rValues[i_gauss] = DVi_DXi;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class SymbolicStokes<SymbolicStokesData<2,3>>;
template class SymbolicStokes<SymbolicStokesData<2,4>>;
template class SymbolicStokes<SymbolicStokesData<3,4>>;
template class SymbolicStokes<SymbolicStokesData<3,6>>;
template class SymbolicStokes<SymbolicStokesData<3,8>>;

} // namespace Kratos
