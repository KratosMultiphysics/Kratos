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
#include "low_mach_navier_stokes.h"
#include "data_containers/low_mach_navier_stokes/low_mach_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
LowMachNavierStokes<TElementData>::LowMachNavierStokes(IndexType NewId)
    : BaseType(NewId)
{}

template <class TElementData>
LowMachNavierStokes<TElementData>::LowMachNavierStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
{}

template <class TElementData>
LowMachNavierStokes<TElementData>::LowMachNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{}

template <class TElementData>
LowMachNavierStokes<TElementData>::LowMachNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{}

template <class TElementData>
LowMachNavierStokes<TElementData>::~LowMachNavierStokes()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer LowMachNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LowMachNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer LowMachNavierStokes<TElementData>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LowMachNavierStokes>(NewId, pGeometry, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int LowMachNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }

    // Check variables used by TElementData
    out = TElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    // Check nodes
    const auto& r_geometry = this->GetGeometry();
    for (const auto& r_node : r_geometry) {
        // Check nodal DOFs
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, r_node);
        KRATOS_CHECK_DOF_IN_NODE(TEMPERATURE, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template<class TElementData>
const Parameters LowMachNavierStokes<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE","TEMPERATURE","DENSITY"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","TEMPERATURE","BODY_FORCE","HEAT_FLUX","REACTION","REACTION_WATER_PRESSURE","REACTION_FLUX","EXTERNAL_PRESSURE","NORMAL"]
        "required_dofs"              : ["VELOCITY_X","VELOCITY_Y","PRESSURE","TEMPERATURE"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Quadrilateral2D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
           "type"        : ["Newtonian2DLaw","Newtonian3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This implements a low Mach approximation Navier-Stokes element with quasi-static Variational MultiScales (VMS) stabilization."
    })");

    return specifications;
}

template <class TElementData>
std::string LowMachNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "LowMachNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void LowMachNavierStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public operations

template< class TElementData >
void LowMachNavierStokes<TElementData>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    const unsigned int p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    const unsigned int x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int t_pos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

    unsigned int local_index = 0;
    for (unsigned int i = 0; i < NumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE, p_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X, x_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y, x_pos+1).EquationId();
        if constexpr (Dim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z, x_pos+2).EquationId();
        }
        rResult[local_index++] = r_geometry[i].GetDof(TEMPERATURE, t_pos).EquationId();
    }
}

template< class TElementData >
void LowMachNavierStokes<TElementData>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize) {
         rElementalDofList.resize(LocalSize);
    }

    const unsigned int p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    const unsigned int x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int t_pos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

    unsigned int local_index = 0;
    for (unsigned int i = 0; i < NumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE, p_pos);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X, x_pos);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y, x_pos+1);
        if constexpr (Dim == 3) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z,x_pos+2);
        }
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(TEMPERATURE, t_pos);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void LowMachNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void LowMachNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void LowMachNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void LowMachNavierStokes<TElementData>::AddBoundaryTraction(
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

    // rData.lhs *= 2.0 * Globals::Pi * y * rData.Weight;
    // noalias(rLHS) -= rData.lhs;
    // noalias(rRHS) += prod(rData.lhs,values);
}

template <>
void LowMachNavierStokes< LowMachNavierStokesData<2,3> >::ComputeGaussPointLHSContribution(
    LowMachNavierStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double bdf0 = rData.bdf0;

    // Nodal data
    const auto& r_t_lin = rData.Temperature;
    const BoundedMatrix<double,2,3> u_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    double tau_c;
    double tau_u;
    double tau_t;
    CalculateStabilizationConstants(rData, tau_c, tau_u, tau_t);

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    //substitute_lhs_2D3N
}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double bdf0 = rData.bdf0;

    // Nodal data
    const auto& r_t_lin = rData.Temperature;
    const BoundedMatrix<double,2,4> u_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    double tau_c;
    double tau_u;
    double tau_t;
    CalculateStabilizationConstants(rData, tau_c, tau_u, tau_t);

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    //substitute_lhs_2D4N
}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_stress = rData.ShearStress;

    // Dynamic parameters
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const auto& r_p = rData.Pressure;
    const auto& r_u = rData.Velocity;
    const auto& r_u_n = rData.VelocityOldStep1;
    const auto& r_u_nn = rData.VelocityOldStep2;
    const auto& r_t = rData.Temperature;
    const auto& r_t_n = rData.TemperatureOldStep1;
    const auto& r_t_nn = rData.TemperatureOldStep2;
    const auto& r_g = rData.BodyForce;
    const auto& r_heat_fl = rData.HeatFlux;

    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 3> u_conv = r_u - r_u_mesh;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    double tau_c;
    double tau_u;
    double tau_t;
    CalculateStabilizationConstants(rData, tau_c, tau_u, tau_t);

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    //substitute_rhs_2D3N
}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_stress = rData.ShearStress;

    // Dynamic parameters
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const auto& r_p = rData.Pressure;
    const auto& r_u = rData.Velocity;
    const auto& r_u_n = rData.VelocityOldStep1;
    const auto& r_u_nn = rData.VelocityOldStep2;
    const auto& r_t = rData.Temperature;
    const auto& r_t_n = rData.TemperatureOldStep1;
    const auto& r_t_nn = rData.TemperatureOldStep2;
    const auto& r_g = rData.BodyForce;
    const auto& r_heat_fl = rData.HeatFlux;

    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 4> u_conv = r_u - r_u_mesh;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    double tau_c;
    double tau_u;
    double tau_t;
    CalculateStabilizationConstants(rData, tau_c, tau_u, tau_t);

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    //substitute_rhs_2D4N
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template< class TElementData >
void LowMachNavierStokes<TElementData>::CalculateStabilizationConstants(
    const TElementData& rData,
    double& rTauPressure,
    double& rTauVelocity,
    double& rTauTemperature)
{
    // Get values from data container
    const double h = rData.ElementSize;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    // const double dt = rData.DeltaTime;
    // const double dyn_tau = rData.DynamicTau;
    const double mu = rData.EffectiveViscosity;
    const double gamma = rData.HeatCapacityRatio;
    const double p_th = rData.ThermodynamicPressure;

    // Data (temperature and convective velocity) interpolation to Gauss points
    double t_gauss = 0.0;
    array_1d<double,3> u_conv_gauss = ZeroVector(3);
    for (IndexType i_node = 0; i_node < TElementData::NumNodes; ++i_node) {
        t_gauss += rData.N[i_node] * rData.Temperature[i_node];
        const auto& r_u_i = row(rData.Velocity, i_node);
        const auto& r_u_m_i = row(rData.MeshVelocity, i_node);
        for (IndexType d = 0; d < TElementData::Dim; ++d) {
            u_conv_gauss[d] += rData.N[i_node] * (r_u_i[d] - r_u_m_i[d]);
        }
    }

    // Calculate stabilization constants at Gauss points
    const double norm_u_conv = norm_2(u_conv_gauss); // Convective velocity norm
    const double rho_gauss = (gamma / c_p / (gamma - 1.0)) * p_th / t_gauss; // Density (from equation of state)

    rTauPressure = mu / rho_gauss + stab_c2 * norm_u_conv * h / stab_c1; // Pressure subscale stabilization operator
    rTauVelocity = 1.0 / (stab_c1 * mu / std::pow(h, 2) + stab_c2 * h * norm_u_conv / h); // Velocity subscale stabilization operator
    rTauTemperature = 1.0 / (stab_c1 * kappa / std::pow(h, 2) + stab_c2 * rho_gauss * c_p * norm_u_conv / h); // Temperature subscale stabilization operator;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void LowMachNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}


template< class TElementData >
void LowMachNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class LowMachNavierStokes< LowMachNavierStokesData<2,3> >;
template class LowMachNavierStokes< LowMachNavierStokesData<2,4> >;
// template class LowMachNavierStokes< LowMachNavierStokesData<3,4> >;
// template class LowMachNavierStokes< LowMachNavierStokesData<3,8> >;

}