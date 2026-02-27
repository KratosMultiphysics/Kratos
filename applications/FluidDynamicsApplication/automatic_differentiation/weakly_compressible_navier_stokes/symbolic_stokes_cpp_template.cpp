//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// External includes

// System includes

// Project includes

// Application includes
#include "symbolic_stokes.h"
#include "data_containers/stokes/stokes_data.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties)
{}

template <class TElementData>
SymbolicStokes<TElementData>::~SymbolicStokes()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int SymbolicStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
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
std::string SymbolicStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "SymbolicStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void SymbolicStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddBoundaryTraction(
    TElementData& rData,
    const Vector& rUnitNormal,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    // Set the current Gauss pt. Voigt notation normal projection matrix
    BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
    FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal, voigt_normal_projection_matrix);

    // Set the current Gauss pt. strain matrix
    BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX, B_matrix);

    // Compute some Gauss pt. auxiliar matrices
    const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
    const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

    // Fill the pressure to Voigt notation operator matrix
    BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            pres_to_voigt_matrix_op(comp, i*BlockSize+Dim) = rData.N[i];
        }
    }

    // Set the shape functions auxiliar transpose matrix
    BoundedMatrix<double, LocalSize, Dim> N_aux_trans = ZeroMatrix(LocalSize, Dim);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            N_aux_trans(i*BlockSize+comp, comp) = rData.N[i];
        }
    }

    // Contribution coming from the shear stress operator
    noalias(rData.lhs) = prod(N_aux_trans, aux_matrix_ACB);

    // Contribution coming from the pressure terms
    const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
    noalias(rData.lhs) -= prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);

    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    rData.lhs *= rData.Weight;
    noalias(rLHS) -= rData.lhs;
    noalias(rRHS) += prod(rData.lhs,values);
}

template <class TElementData>
void SymbolicStokes<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY
    rOutput = 0.0;

    if (rVariable == HEAT_FLUX) {
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        rOutput = 0.0;
        // Shape functions
        const auto &r_geometry = this->GetGeometry();

        const unsigned int num_nodes = r_geometry.PointsNumber();
        Vector data_N(num_nodes);
        for (unsigned int i = 0; i < num_nodes; i++) {
            data_N[i] = 1.0 / static_cast<double>(num_nodes);
        }
        data.N = data_N;
        // Shape functions gradients
        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);
        data.DN_DX = shape_derivatives[0]; // Note: Valid only for linear elements
        // Compute strain and stress
        this->CalculateMaterialResponse(data);
        // Compute dissipation
        rOutput = inner_prod(data.ShearStress, data.StrainRate);
    }
    else if (rVariable == EQ_STRAIN_RATE || rVariable == EFFECTIVE_VISCOSITY) {
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        rOutput = 0.0;
        // Shape functions
        const auto &r_geometry = this->GetGeometry();
        const unsigned int num_nodes = r_geometry.PointsNumber();
        Vector data_N(num_nodes);
        for (unsigned int i = 0; i < num_nodes; i++) {
            data_N[i] = 1.0 / static_cast<double>(num_nodes);
        }
        data.N = data_N;
        // Shape function gradients
        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);
        data.DN_DX = shape_derivatives[0]; // Note: Valid only for linear elements
        // Compute strain and stress
        this->CalculateMaterialResponse(data);
        const auto &p_constitutive_law = this->GetConstitutiveLaw();
        rOutput = p_constitutive_law->GetValue(rVariable, rOutput);
    }
    else {
        this->Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template <>
void SymbolicStokes< SymbolicStokesData<2,3> >::ComputeGaussPointLHSContribution(
    SymbolicStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_2D3N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes< SymbolicStokesData<2,4> >::ComputeGaussPointLHSContribution(
    SymbolicStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_2D4N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,4>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_3D4N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,6>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,6>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_3D6N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,8>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,8>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_3D8N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,3>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_2D3N

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,4>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<2,4>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_2D4N

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,4>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_3D4N

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,6>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,6>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_3D6N

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,8>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,8>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_3D8N

    noalias(rRHS) += rData.Weight * rhs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void SymbolicStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void SymbolicStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class SymbolicStokes< SymbolicStokesData<2,3> >;
template class SymbolicStokes< SymbolicStokesData<2,4> >;
template class SymbolicStokes< SymbolicStokesData<3,4> >;
template class SymbolicStokes< SymbolicStokesData<3,6> >;
template class SymbolicStokes< SymbolicStokesData<3,8> >;

}
