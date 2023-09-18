//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes


// Application includes
#include "axisymmetric_navier_stokes.h"
#include "data_containers/axisymmetric_navier_stokes/axisymmetric_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::~AxisymmetricNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer AxisymmetricNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer AxisymmetricNavierStokes<TElementData>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricNavierStokes>(NewId, pGeometry, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int AxisymmetricNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

template<class TElementData>
const Parameters AxisymmetricNavierStokes<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE","DENSITY","DYNAMIC_VISCOSITY"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","IS_STRUCTURE","DISPLACEMENT","BODY_FORCE","NODAL_AREA","NODAL_H","ADVPROJ","DIVPROJ","REACTION","REACTION_WATER_PRESSURE","EXTERNAL_PRESSURE","NORMAL","Y_WALL","Q_VALUE"]
        "required_dofs"              : ["VELOCITY_X","VELOCITY_Y","PRESSURE"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Quadrilateral2D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : [],
            "dimension"   : [],
            "strain_size" : [3]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This implements an axisymmetric Navier-Stokes element with quasi-static Variational MultiScales (VMS) stabilization. Viscous behavior is hardcoded to a Newtonian constitutive model. x-direction is assumed to be aligned with the revolution axis meaning that y-direction represents the radial one."
    })");

    return specifications;
}

template <class TElementData>
std::string AxisymmetricNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "AxisymmetricNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
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
void AxisymmetricNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::AddBoundaryTraction(
    TElementData& rData,
    const Vector& rUnitNormal,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    KRATOS_ERROR << "To be implemented" << std::endl;

    // // Set the current Gauss pt. Voigt notation normal projection matrix
    // BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
    // FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal, voigt_normal_projection_matrix);

    // // Set the current Gauss pt. strain matrix
    // BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    // FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX, B_matrix);

    // // Compute some Gauss pt. auxiliar matrices
    // const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
    // const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

    // // Fill the pressure to Voigt notation operator matrix
    // BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
    // for (unsigned int i=0; i<NumNodes; ++i) {
    //     for (unsigned int comp=0; comp<Dim; ++comp) {
    //         pres_to_voigt_matrix_op(comp, i*BlockSize+Dim) = rData.N[i];
    //     }
    // }

    // // Set the shape functions auxiliar transpose matrix
    // BoundedMatrix<double, LocalSize, Dim> N_aux_trans = ZeroMatrix(LocalSize, Dim);
    // for (unsigned int i=0; i<NumNodes; ++i) {
    //     for (unsigned int comp=0; comp<Dim; ++comp) {
    //         N_aux_trans(i*BlockSize+comp, comp) = rData.N[i];
    //     }
    // }

    // // Contribution coming fron the shear stress operator
    // noalias(rData.lhs) = prod(N_aux_trans, aux_matrix_ACB);

    // // Contribution coming from the pressure terms
    // const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
    // noalias(rData.lhs) -= prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);

    // array_1d<double,LocalSize> values;
    // this->GetCurrentValuesVector(rData,values);

    // rData.lhs *= rData.Weight;
    // noalias(rLHS) -= rData.lhs;
    // noalias(rRHS) += prod(rData.lhs,values);
}

template <>
void AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,3> >::ComputeGaussPointLHSContribution(
    AxisymmetricNavierStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3> vconv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    //substitute_lhs_2D3N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    AxisymmetricNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,4> vconv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    //substitute_lhs_2D4N

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    AxisymmetricNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3>& v = rData.Velocity;
    const BoundedMatrix<double,2,3>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,3>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = v - vmesh;
    const BoundedMatrix<double,2,3>& f = rData.BodyForce;
    const array_1d<double,3>& p = rData.Pressure;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    //substitute_rhs_2D3N

    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    AxisymmetricNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,4>& v = rData.Velocity;
    const BoundedMatrix<double,2,4>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,4>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> vconv = v - vmesh;
    const BoundedMatrix<double,2,4>& f = rData.BodyForce;
    const array_1d<double,4>& p = rData.Pressure;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    //substitute_rhs_2D4N

    noalias(rRHS) += rData.Weight * rhs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void AxisymmetricNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void AxisymmetricNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,3> >;
template class AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,4> >;

}