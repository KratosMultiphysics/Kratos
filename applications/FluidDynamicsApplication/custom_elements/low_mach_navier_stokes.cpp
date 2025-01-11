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
    const double sigma = rData.Resistance;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double bdf0 = rData.bdf0;

    // Nodal data
    const auto& r_u = rData.Velocity;
    const auto& r_t = rData.Temperature;
    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> lin_u_conv = rData.Velocity - rData.MeshVelocity;

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

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = 1.0/(gamma - 1.0);
const double crLHS3 = crLHS2*gamma*p_th;
const double crLHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crLHS5 = 1.0/crLHS4;
const double crLHS6 = 1.0/c_p;
const double crLHS7 = crLHS5*crLHS6;
const double crLHS8 = crLHS3*crLHS7;
const double crLHS9 = crLHS8*gauss_weight;
const double crLHS10 = crLHS9*tau_u;
const double crLHS11 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS12 = r_N[0]*r_N[0];
const double crLHS13 = crLHS12*crLHS5;
const double crLHS14 = r_N[0]*sigma;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS3*(crLHS15 + crLHS18);
const double crLHS20 = tau_u*(crLHS14 + crLHS19*crLHS7);
const double crLHS21 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS22 = bdf0*crLHS3;
const double crLHS23 = 1.0/(crLHS4*crLHS4);
const double crLHS24 = crLHS23*crLHS6;
const double crLHS25 = crLHS22*crLHS24*gauss_weight;
const double crLHS26 = r_DN(0,0)*r_DN(1,0);
const double crLHS27 = r_DN(0,1)*r_DN(1,1);
const double crLHS28 = crLHS10*(crLHS26 + crLHS27);
const double crLHS29 = r_DN(1,0)*r_N[0];
const double crLHS30 = crLHS5*r_N[0];
const double crLHS31 = crLHS30*r_N[1];
const double crLHS32 = -crLHS11*crLHS31;
const double crLHS33 = r_N[1]*sigma;
const double crLHS34 = bdf0*r_N[1];
const double crLHS35 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS36 = crLHS3*(crLHS34 + crLHS35);
const double crLHS37 = tau_u*(crLHS33 + crLHS36*crLHS7);
const double crLHS38 = r_DN(1,1)*r_N[0];
const double crLHS39 = -crLHS21*crLHS31;
const double crLHS40 = crLHS3*gauss_weight;
const double crLHS41 = crLHS40*r_N[1];
const double crLHS42 = crLHS15*crLHS24;
const double crLHS43 = -crLHS41*crLHS42;
const double crLHS44 = r_DN(0,0)*r_DN(2,0);
const double crLHS45 = r_DN(0,1)*r_DN(2,1);
const double crLHS46 = crLHS10*(crLHS44 + crLHS45);
const double crLHS47 = r_DN(2,0)*r_N[0];
const double crLHS48 = crLHS30*r_N[2];
const double crLHS49 = -crLHS11*crLHS48;
const double crLHS50 = r_N[2]*sigma;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS3*(crLHS51 + crLHS52);
const double crLHS54 = tau_u*(crLHS50 + crLHS53*crLHS7);
const double crLHS55 = r_DN(2,1)*r_N[0];
const double crLHS56 = -crLHS21*crLHS48;
const double crLHS57 = crLHS40*r_N[2];
const double crLHS58 = -crLHS42*crLHS57;
const double crLHS59 = crLHS14*tau_u;
const double crLHS60 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crLHS61 = crLHS3*(crLHS11*crLHS16 + crLHS17*crLHS21);
const double crLHS62 = crLHS61*r_N[0];
const double crLHS63 = crLHS24*tau_u;
const double crLHS64 = gauss_weight*(crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_N[0]*tau_u - crLHS59 - crLHS62*crLHS63 - r_N[0]);
const double crLHS65 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS66 = r_C(0,2)*r_DN(0,0);
const double crLHS67 = crLHS66 + r_C(2,2)*r_DN(0,1);
const double crLHS68 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crLHS69 = crLHS30*crLHS68;
const double crLHS70 = -crLHS69 + r_DN(0,0);
const double crLHS71 = crLHS8*tau_c;
const double crLHS72 = crLHS71*r_DN(0,0);
const double crLHS73 = crLHS13*crLHS22;
const double crLHS74 = crLHS3*crLHS30;
const double crLHS75 = crLHS6*crLHS74;
const double crLHS76 = crLHS18*crLHS8;
const double crLHS77 = crLHS60*crLHS75;
const double crLHS78 = crLHS24*crLHS62;
const double crLHS79 = crLHS12*sigma - crLHS14*crLHS20 + crLHS18*crLHS75 + crLHS20*crLHS76 + crLHS20*crLHS77 - crLHS20*crLHS78 + crLHS6*crLHS73;
const double crLHS80 = crLHS66 + r_C(0,1)*r_DN(0,1);
const double crLHS81 = r_C(1,2)*r_DN(0,1);
const double crLHS82 = crLHS81 + r_C(2,2)*r_DN(0,0);
const double crLHS83 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crLHS84 = crLHS30*crLHS83;
const double crLHS85 = -crLHS84 + r_DN(0,1);
const double crLHS86 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0));
const double crLHS87 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1));
const double crLHS88 = crLHS86*r_DN(0,0) + crLHS87*r_DN(0,1);
const double crLHS89 = crLHS15 + crLHS88;
const double crLHS90 = crLHS24*crLHS40*tau_c;
const double crLHS91 = crLHS90*r_DN(0,0);
const double crLHS92 = r_DN(0,0)*r_N[1];
const double crLHS93 = crLHS61*crLHS63;
const double crLHS94 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS95 = r_C(0,2)*r_DN(1,0);
const double crLHS96 = crLHS95 + r_C(2,2)*r_DN(1,1);
const double crLHS97 = crLHS5*r_N[1];
const double crLHS98 = crLHS68*crLHS97;
const double crLHS99 = -crLHS98 + r_DN(1,0);
const double crLHS100 = crLHS3*crLHS97;
const double crLHS101 = crLHS100*crLHS15;
const double crLHS102 = crLHS101*crLHS6 + crLHS14*r_N[1];
const double crLHS103 = crLHS102 - crLHS14*crLHS37 + crLHS35*crLHS75 + crLHS37*crLHS76 + crLHS37*crLHS77 - crLHS37*crLHS78;
const double crLHS104 = crLHS95 + r_C(0,1)*r_DN(1,1);
const double crLHS105 = r_C(1,2)*r_DN(1,1);
const double crLHS106 = crLHS105 + r_C(2,2)*r_DN(1,0);
const double crLHS107 = crLHS83*crLHS97;
const double crLHS108 = -crLHS107 + r_DN(1,1);
const double crLHS109 = crLHS86*r_DN(1,0) + crLHS87*r_DN(1,1);
const double crLHS110 = crLHS109 + crLHS34;
const double crLHS111 = r_DN(0,0)*r_N[2];
const double crLHS112 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS113 = r_C(0,2)*r_DN(2,0);
const double crLHS114 = crLHS113 + r_C(2,2)*r_DN(2,1);
const double crLHS115 = crLHS5*r_N[2];
const double crLHS116 = -crLHS115*crLHS68 + r_DN(2,0);
const double crLHS117 = crLHS115*crLHS3;
const double crLHS118 = crLHS117*crLHS15;
const double crLHS119 = crLHS118*crLHS6 + crLHS14*r_N[2];
const double crLHS120 = crLHS119 - crLHS14*crLHS54 + crLHS52*crLHS75 + crLHS54*crLHS76 + crLHS54*crLHS77 - crLHS54*crLHS78;
const double crLHS121 = crLHS113 + r_C(0,1)*r_DN(2,1);
const double crLHS122 = r_C(1,2)*r_DN(2,1);
const double crLHS123 = crLHS122 + r_C(2,2)*r_DN(2,0);
const double crLHS124 = -crLHS115*crLHS83 + r_DN(2,1);
const double crLHS125 = crLHS86*r_DN(2,0) + crLHS87*r_DN(2,1);
const double crLHS126 = crLHS125 + crLHS51;
const double crLHS127 = crLHS81 + r_C(0,1)*r_DN(0,0);
const double crLHS128 = crLHS71*r_DN(0,1);
const double crLHS129 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS130 = crLHS90*r_DN(0,1);
const double crLHS131 = r_DN(0,1)*r_N[1];
const double crLHS132 = crLHS105 + r_C(0,1)*r_DN(1,0);
const double crLHS133 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS134 = r_DN(0,1)*r_N[2];
const double crLHS135 = crLHS122 + r_C(0,1)*r_DN(2,0);
const double crLHS136 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS137 = crLHS13*crLHS40;
const double crLHS138 = -crLHS19 + dp_th_dt*r_N[0];
const double crLHS139 = tau_t*1.0/(crLHS4*crLHS4*crLHS4);
const double crLHS140 = crLHS139*crLHS62;
const double crLHS141 = crLHS41*crLHS69;
const double crLHS142 = crLHS41*crLHS84;
const double crLHS143 = -crLHS36 + dp_th_dt*r_N[1];
const double crLHS144 = crLHS101 + crLHS26*kappa + crLHS27*kappa - crLHS5*dp_th_dt*r_N[0]*r_N[1];
const double crLHS145 = crLHS57*crLHS69;
const double crLHS146 = crLHS57*crLHS84;
const double crLHS147 = -crLHS53 + dp_th_dt*r_N[2];
const double crLHS148 = crLHS118 + crLHS44*kappa + crLHS45*kappa - crLHS5*dp_th_dt*r_N[0]*r_N[2];
const double crLHS149 = r_DN(1,0)*r_DN(1,0);
const double crLHS150 = r_DN(1,1)*r_DN(1,1);
const double crLHS151 = r_N[1]*r_N[1];
const double crLHS152 = crLHS151*crLHS5;
const double crLHS153 = r_DN(1,0)*r_DN(2,0);
const double crLHS154 = r_DN(1,1)*r_DN(2,1);
const double crLHS155 = crLHS10*(crLHS153 + crLHS154);
const double crLHS156 = r_DN(2,0)*r_N[1];
const double crLHS157 = crLHS97*r_N[2];
const double crLHS158 = -crLHS11*crLHS157;
const double crLHS159 = r_DN(2,1)*r_N[1];
const double crLHS160 = -crLHS157*crLHS21;
const double crLHS161 = -crLHS24*crLHS34*crLHS57;
const double crLHS162 = crLHS33*tau_u;
const double crLHS163 = crLHS71*r_DN(1,0);
const double crLHS164 = crLHS100*crLHS6;
const double crLHS165 = crLHS35*crLHS8;
const double crLHS166 = crLHS164*crLHS60;
const double crLHS167 = crLHS61*r_N[1];
const double crLHS168 = crLHS167*crLHS24;
const double crLHS169 = crLHS102 + crLHS164*crLHS18 + crLHS165*crLHS20 + crLHS166*crLHS20 - crLHS168*crLHS20 - crLHS20*crLHS33;
const double crLHS170 = crLHS90*r_DN(1,0);
const double crLHS171 = gauss_weight*(-crLHS162 - crLHS167*crLHS63 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS172 = crLHS152*crLHS22;
const double crLHS173 = crLHS151*sigma + crLHS164*crLHS35 + crLHS165*crLHS37 + crLHS166*crLHS37 - crLHS168*crLHS37 + crLHS172*crLHS6 - crLHS33*crLHS37;
const double crLHS174 = r_DN(1,0)*r_N[2];
const double crLHS175 = crLHS117*crLHS34;
const double crLHS176 = crLHS175*crLHS6 + crLHS33*r_N[2];
const double crLHS177 = crLHS164*crLHS52 + crLHS165*crLHS54 + crLHS166*crLHS54 - crLHS168*crLHS54 + crLHS176 - crLHS33*crLHS54;
const double crLHS178 = crLHS71*r_DN(1,1);
const double crLHS179 = crLHS90*r_DN(1,1);
const double crLHS180 = r_DN(1,1)*r_N[2];
const double crLHS181 = crLHS139*crLHS167;
const double crLHS182 = crLHS152*crLHS40;
const double crLHS183 = crLHS57*crLHS98;
const double crLHS184 = crLHS107*crLHS57;
const double crLHS185 = crLHS153*kappa + crLHS154*kappa + crLHS175 - crLHS5*dp_th_dt*r_N[1]*r_N[2];
const double crLHS186 = r_DN(2,0)*r_DN(2,0);
const double crLHS187 = r_DN(2,1)*r_DN(2,1);
const double crLHS188 = r_N[2]*r_N[2];
const double crLHS189 = crLHS188*crLHS5;
const double crLHS190 = crLHS50*tau_u;
const double crLHS191 = crLHS71*r_DN(2,0);
const double crLHS192 = crLHS117*crLHS6;
const double crLHS193 = crLHS52*crLHS8;
const double crLHS194 = crLHS192*crLHS60;
const double crLHS195 = crLHS61*r_N[2];
const double crLHS196 = crLHS195*crLHS24;
const double crLHS197 = crLHS119 + crLHS18*crLHS192 + crLHS193*crLHS20 + crLHS194*crLHS20 - crLHS196*crLHS20 - crLHS20*crLHS50;
const double crLHS198 = crLHS90*r_DN(2,0);
const double crLHS199 = crLHS176 + crLHS192*crLHS35 + crLHS193*crLHS37 + crLHS194*crLHS37 - crLHS196*crLHS37 - crLHS37*crLHS50;
const double crLHS200 = gauss_weight*(-crLHS190 - crLHS195*crLHS63 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS201 = crLHS189*crLHS22;
const double crLHS202 = crLHS188*sigma + crLHS192*crLHS52 + crLHS193*crLHS54 + crLHS194*crLHS54 - crLHS196*crLHS54 + crLHS201*crLHS6 - crLHS50*crLHS54;
const double crLHS203 = crLHS71*r_DN(2,1);
const double crLHS204 = crLHS90*r_DN(2,1);
const double crLHS205 = crLHS139*crLHS195;
const double crLHS206 = crLHS189*crLHS40;
rLHS(0,0)+=crLHS10*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS9*(-crLHS11*crLHS13 + crLHS20*r_DN(0,0) + r_DN(0,0)*r_N[0]);
rLHS(0,2)+=crLHS9*(-crLHS13*crLHS21 + crLHS20*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS12*crLHS25;
rLHS(0,4)+=crLHS28;
rLHS(0,5)+=crLHS9*(crLHS29 + crLHS32 + crLHS37*r_DN(0,0));
rLHS(0,6)+=crLHS9*(crLHS37*r_DN(0,1) + crLHS38 + crLHS39);
rLHS(0,7)+=crLHS43;
rLHS(0,8)+=crLHS46;
rLHS(0,9)+=crLHS9*(crLHS47 + crLHS49 + crLHS54*r_DN(0,0));
rLHS(0,10)+=crLHS9*(crLHS54*r_DN(0,1) + crLHS55 + crLHS56);
rLHS(0,11)+=crLHS58;
rLHS(1,0)+=crLHS64*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS65*r_DN(0,0) + crLHS67*r_DN(0,1) + crLHS70*crLHS72 + crLHS79);
rLHS(1,2)+=gauss_weight*(crLHS72*crLHS85 + crLHS80*r_DN(0,0) + crLHS82*r_DN(0,1));
rLHS(1,3)+=-crLHS89*crLHS91;
rLHS(1,4)+=gauss_weight*(crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS29*crLHS93 - crLHS59*r_DN(1,0) - crLHS92);
rLHS(1,5)+=gauss_weight*(crLHS103 + crLHS72*crLHS99 + crLHS94*r_DN(0,0) + crLHS96*r_DN(0,1));
rLHS(1,6)+=gauss_weight*(crLHS104*r_DN(0,0) + crLHS106*r_DN(0,1) + crLHS108*crLHS72);
rLHS(1,7)+=-crLHS110*crLHS91;
rLHS(1,8)+=gauss_weight*(-crLHS111 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS47*crLHS93 - crLHS59*r_DN(2,0));
rLHS(1,9)+=gauss_weight*(crLHS112*r_DN(0,0) + crLHS114*r_DN(0,1) + crLHS116*crLHS72 + crLHS120);
rLHS(1,10)+=gauss_weight*(crLHS121*r_DN(0,0) + crLHS123*r_DN(0,1) + crLHS124*crLHS72);
rLHS(1,11)+=-crLHS126*crLHS91;
rLHS(2,0)+=crLHS64*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS127*r_DN(0,1) + crLHS128*crLHS70 + crLHS67*r_DN(0,0));
rLHS(2,2)+=gauss_weight*(crLHS128*crLHS85 + crLHS129*r_DN(0,1) + crLHS79 + crLHS82*r_DN(0,0));
rLHS(2,3)+=-crLHS130*crLHS89;
rLHS(2,4)+=gauss_weight*(-crLHS131 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS38*crLHS93 - crLHS59*r_DN(1,1));
rLHS(2,5)+=gauss_weight*(crLHS128*crLHS99 + crLHS132*r_DN(0,1) + crLHS96*r_DN(0,0));
rLHS(2,6)+=gauss_weight*(crLHS103 + crLHS106*r_DN(0,0) + crLHS108*crLHS128 + crLHS133*r_DN(0,1));
rLHS(2,7)+=-crLHS110*crLHS130;
rLHS(2,8)+=gauss_weight*(-crLHS134 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS55*crLHS93 - crLHS59*r_DN(2,1));
rLHS(2,9)+=gauss_weight*(crLHS114*r_DN(0,0) + crLHS116*crLHS128 + crLHS135*r_DN(0,1));
rLHS(2,10)+=gauss_weight*(crLHS120 + crLHS123*r_DN(0,0) + crLHS124*crLHS128 + crLHS136*r_DN(0,1));
rLHS(2,11)+=-crLHS126*crLHS130;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS137*crLHS68;
rLHS(3,2)+=crLHS137*crLHS83;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa + crLHS12*crLHS5*dp_th_dt - crLHS138*crLHS140 + crLHS138*crLHS18*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS138*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS138*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS73 - crLHS74*crLHS88);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS141;
rLHS(3,6)+=crLHS142;
rLHS(3,7)+=-gauss_weight*(-crLHS109*crLHS74 - crLHS140*crLHS143 + crLHS143*crLHS18*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS143*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS143*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS144);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS145;
rLHS(3,10)+=crLHS146;
rLHS(3,11)+=-gauss_weight*(-crLHS125*crLHS74 - crLHS140*crLHS147 + crLHS147*crLHS18*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS147*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS147*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS148);
rLHS(4,0)+=crLHS28;
rLHS(4,1)+=crLHS9*(crLHS20*r_DN(1,0) + crLHS32 + crLHS92);
rLHS(4,2)+=crLHS9*(crLHS131 + crLHS20*r_DN(1,1) + crLHS39);
rLHS(4,3)+=crLHS43;
rLHS(4,4)+=crLHS10*(crLHS149 + crLHS150);
rLHS(4,5)+=crLHS9*(-crLHS11*crLHS152 + crLHS37*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS9*(-crLHS152*crLHS21 + crLHS37*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS151*crLHS25;
rLHS(4,8)+=crLHS155;
rLHS(4,9)+=crLHS9*(crLHS156 + crLHS158 + crLHS54*r_DN(1,0));
rLHS(4,10)+=crLHS9*(crLHS159 + crLHS160 + crLHS54*r_DN(1,1));
rLHS(4,11)+=crLHS161;
rLHS(5,0)+=gauss_weight*(-crLHS162*r_DN(0,0) + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS29 - crLHS92*crLHS93);
rLHS(5,1)+=gauss_weight*(crLHS163*crLHS70 + crLHS169 + crLHS65*r_DN(1,0) + crLHS67*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS163*crLHS85 + crLHS80*r_DN(1,0) + crLHS82*r_DN(1,1));
rLHS(5,3)+=-crLHS170*crLHS89;
rLHS(5,4)+=crLHS171*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS163*crLHS99 + crLHS173 + crLHS94*r_DN(1,0) + crLHS96*r_DN(1,1));
rLHS(5,6)+=gauss_weight*(crLHS104*r_DN(1,0) + crLHS106*r_DN(1,1) + crLHS108*crLHS163);
rLHS(5,7)+=-crLHS110*crLHS170;
rLHS(5,8)+=gauss_weight*(-crLHS156*crLHS93 - crLHS162*r_DN(2,0) - crLHS174 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS112*r_DN(1,0) + crLHS114*r_DN(1,1) + crLHS116*crLHS163 + crLHS177);
rLHS(5,10)+=gauss_weight*(crLHS121*r_DN(1,0) + crLHS123*r_DN(1,1) + crLHS124*crLHS163);
rLHS(5,11)+=-crLHS126*crLHS170;
rLHS(6,0)+=gauss_weight*(-crLHS131*crLHS93 - crLHS162*r_DN(0,1) + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS38);
rLHS(6,1)+=gauss_weight*(crLHS127*r_DN(1,1) + crLHS178*crLHS70 + crLHS67*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS129*r_DN(1,1) + crLHS169 + crLHS178*crLHS85 + crLHS82*r_DN(1,0));
rLHS(6,3)+=-crLHS179*crLHS89;
rLHS(6,4)+=crLHS171*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS132*r_DN(1,1) + crLHS178*crLHS99 + crLHS96*r_DN(1,0));
rLHS(6,6)+=gauss_weight*(crLHS106*r_DN(1,0) + crLHS108*crLHS178 + crLHS133*r_DN(1,1) + crLHS173);
rLHS(6,7)+=-crLHS110*crLHS179;
rLHS(6,8)+=gauss_weight*(-crLHS159*crLHS93 - crLHS162*r_DN(2,1) - crLHS180 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS114*r_DN(1,0) + crLHS116*crLHS178 + crLHS135*r_DN(1,1));
rLHS(6,10)+=gauss_weight*(crLHS123*r_DN(1,0) + crLHS124*crLHS178 + crLHS136*r_DN(1,1) + crLHS177);
rLHS(6,11)+=-crLHS126*crLHS179;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS141;
rLHS(7,2)+=crLHS142;
rLHS(7,3)+=-gauss_weight*(-crLHS100*crLHS88 - crLHS138*crLHS181 + crLHS138*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS138*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS138*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS144);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS182*crLHS68;
rLHS(7,6)+=crLHS182*crLHS83;
rLHS(7,7)+=-gauss_weight*(-crLHS100*crLHS109 - crLHS143*crLHS181 + crLHS143*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS143*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS143*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS149*kappa - crLHS150*kappa + crLHS151*crLHS5*dp_th_dt - crLHS172);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS183;
rLHS(7,10)+=crLHS184;
rLHS(7,11)+=-gauss_weight*(-crLHS100*crLHS125 - crLHS147*crLHS181 + crLHS147*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS147*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS147*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS185);
rLHS(8,0)+=crLHS46;
rLHS(8,1)+=crLHS9*(crLHS111 + crLHS20*r_DN(2,0) + crLHS49);
rLHS(8,2)+=crLHS9*(crLHS134 + crLHS20*r_DN(2,1) + crLHS56);
rLHS(8,3)+=crLHS58;
rLHS(8,4)+=crLHS155;
rLHS(8,5)+=crLHS9*(crLHS158 + crLHS174 + crLHS37*r_DN(2,0));
rLHS(8,6)+=crLHS9*(crLHS160 + crLHS180 + crLHS37*r_DN(2,1));
rLHS(8,7)+=crLHS161;
rLHS(8,8)+=crLHS10*(crLHS186 + crLHS187);
rLHS(8,9)+=crLHS9*(-crLHS11*crLHS189 + crLHS54*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS9*(-crLHS189*crLHS21 + crLHS54*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS188*crLHS25;
rLHS(9,0)+=gauss_weight*(-crLHS111*crLHS93 - crLHS190*r_DN(0,0) + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS47);
rLHS(9,1)+=gauss_weight*(crLHS191*crLHS70 + crLHS197 + crLHS65*r_DN(2,0) + crLHS67*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS191*crLHS85 + crLHS80*r_DN(2,0) + crLHS82*r_DN(2,1));
rLHS(9,3)+=-crLHS198*crLHS89;
rLHS(9,4)+=gauss_weight*(-crLHS156 - crLHS174*crLHS93 - crLHS190*r_DN(1,0) + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS191*crLHS99 + crLHS199 + crLHS94*r_DN(2,0) + crLHS96*r_DN(2,1));
rLHS(9,6)+=gauss_weight*(crLHS104*r_DN(2,0) + crLHS106*r_DN(2,1) + crLHS108*crLHS191);
rLHS(9,7)+=-crLHS110*crLHS198;
rLHS(9,8)+=crLHS200*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS112*r_DN(2,0) + crLHS114*r_DN(2,1) + crLHS116*crLHS191 + crLHS202);
rLHS(9,10)+=gauss_weight*(crLHS121*r_DN(2,0) + crLHS123*r_DN(2,1) + crLHS124*crLHS191);
rLHS(9,11)+=-crLHS126*crLHS198;
rLHS(10,0)+=gauss_weight*(-crLHS134*crLHS93 - crLHS190*r_DN(0,1) + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS55);
rLHS(10,1)+=gauss_weight*(crLHS127*r_DN(2,1) + crLHS203*crLHS70 + crLHS67*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS129*r_DN(2,1) + crLHS197 + crLHS203*crLHS85 + crLHS82*r_DN(2,0));
rLHS(10,3)+=-crLHS204*crLHS89;
rLHS(10,4)+=gauss_weight*(-crLHS159 - crLHS180*crLHS93 - crLHS190*r_DN(1,1) + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS60*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS132*r_DN(2,1) + crLHS203*crLHS99 + crLHS96*r_DN(2,0));
rLHS(10,6)+=gauss_weight*(crLHS106*r_DN(2,0) + crLHS108*crLHS203 + crLHS133*r_DN(2,1) + crLHS199);
rLHS(10,7)+=-crLHS110*crLHS204;
rLHS(10,8)+=crLHS200*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS114*r_DN(2,0) + crLHS116*crLHS203 + crLHS135*r_DN(2,1));
rLHS(10,10)+=gauss_weight*(crLHS123*r_DN(2,0) + crLHS124*crLHS203 + crLHS136*r_DN(2,1) + crLHS202);
rLHS(10,11)+=-crLHS126*crLHS204;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS145;
rLHS(11,2)+=crLHS146;
rLHS(11,3)+=-gauss_weight*(-crLHS117*crLHS88 + crLHS138*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS138*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[2]*tau_t - crLHS138*crLHS205 + crLHS138*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS148);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS183;
rLHS(11,6)+=crLHS184;
rLHS(11,7)+=-gauss_weight*(-crLHS109*crLHS117 + crLHS143*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS143*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[2]*tau_t - crLHS143*crLHS205 + crLHS143*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS185);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS206*crLHS68;
rLHS(11,10)+=crLHS206*crLHS83;
rLHS(11,11)+=-gauss_weight*(-crLHS117*crLHS125 + crLHS147*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS147*crLHS2*crLHS23*crLHS60*gamma*p_th*r_N[2]*tau_t - crLHS147*crLHS205 + crLHS147*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS186*kappa - crLHS187*kappa + crLHS188*crLHS5*dp_th_dt - crLHS201);

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double sigma = rData.Resistance;
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
    const auto& r_u = rData.Velocity;
    const auto& r_t = rData.Temperature;
    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> lin_u_conv = rData.Velocity - rData.MeshVelocity;

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

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = 1.0/(gamma - 1.0);
const double crLHS3 = crLHS2*gamma*p_th;
const double crLHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crLHS5 = 1.0/crLHS4;
const double crLHS6 = 1.0/c_p;
const double crLHS7 = crLHS5*crLHS6;
const double crLHS8 = crLHS3*crLHS7;
const double crLHS9 = crLHS8*gauss_weight;
const double crLHS10 = crLHS9*tau_u;
const double crLHS11 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS12 = r_N[0]*r_N[0];
const double crLHS13 = crLHS12*crLHS5;
const double crLHS14 = r_N[0]*sigma;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2] + lin_u_conv(3,0)*r_N[3];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2] + lin_u_conv(3,1)*r_N[3];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS3*(crLHS15 + crLHS18);
const double crLHS20 = tau_u*(crLHS14 + crLHS19*crLHS7);
const double crLHS21 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS22 = bdf0*crLHS3;
const double crLHS23 = 1.0/(crLHS4*crLHS4);
const double crLHS24 = crLHS23*crLHS6;
const double crLHS25 = crLHS22*crLHS24*gauss_weight;
const double crLHS26 = r_DN(0,0)*r_DN(1,0);
const double crLHS27 = r_DN(0,1)*r_DN(1,1);
const double crLHS28 = crLHS10*(crLHS26 + crLHS27);
const double crLHS29 = r_DN(1,0)*r_N[0];
const double crLHS30 = crLHS5*r_N[0];
const double crLHS31 = crLHS30*r_N[1];
const double crLHS32 = -crLHS11*crLHS31;
const double crLHS33 = r_N[1]*sigma;
const double crLHS34 = bdf0*r_N[1];
const double crLHS35 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS36 = crLHS3*(crLHS34 + crLHS35);
const double crLHS37 = tau_u*(crLHS33 + crLHS36*crLHS7);
const double crLHS38 = r_DN(1,1)*r_N[0];
const double crLHS39 = -crLHS21*crLHS31;
const double crLHS40 = crLHS3*gauss_weight;
const double crLHS41 = crLHS40*r_N[1];
const double crLHS42 = crLHS15*crLHS24;
const double crLHS43 = -crLHS41*crLHS42;
const double crLHS44 = r_DN(0,0)*r_DN(2,0);
const double crLHS45 = r_DN(0,1)*r_DN(2,1);
const double crLHS46 = crLHS10*(crLHS44 + crLHS45);
const double crLHS47 = r_DN(2,0)*r_N[0];
const double crLHS48 = crLHS30*r_N[2];
const double crLHS49 = -crLHS11*crLHS48;
const double crLHS50 = r_N[2]*sigma;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS3*(crLHS51 + crLHS52);
const double crLHS54 = tau_u*(crLHS50 + crLHS53*crLHS7);
const double crLHS55 = r_DN(2,1)*r_N[0];
const double crLHS56 = -crLHS21*crLHS48;
const double crLHS57 = crLHS40*r_N[2];
const double crLHS58 = -crLHS42*crLHS57;
const double crLHS59 = r_DN(0,0)*r_DN(3,0);
const double crLHS60 = r_DN(0,1)*r_DN(3,1);
const double crLHS61 = crLHS10*(crLHS59 + crLHS60);
const double crLHS62 = r_DN(3,0)*r_N[0];
const double crLHS63 = crLHS30*r_N[3];
const double crLHS64 = -crLHS11*crLHS63;
const double crLHS65 = r_N[3]*sigma;
const double crLHS66 = bdf0*r_N[3];
const double crLHS67 = crLHS16*r_DN(3,0) + crLHS17*r_DN(3,1);
const double crLHS68 = crLHS3*(crLHS66 + crLHS67);
const double crLHS69 = tau_u*(crLHS65 + crLHS68*crLHS7);
const double crLHS70 = r_DN(3,1)*r_N[0];
const double crLHS71 = -crLHS21*crLHS63;
const double crLHS72 = crLHS40*r_N[3];
const double crLHS73 = -crLHS42*crLHS72;
const double crLHS74 = crLHS14*tau_u;
const double crLHS75 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crLHS76 = crLHS3*(crLHS11*crLHS16 + crLHS17*crLHS21);
const double crLHS77 = crLHS76*r_N[0];
const double crLHS78 = crLHS24*tau_u;
const double crLHS79 = gauss_weight*(crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_N[0]*tau_u - crLHS74 - crLHS77*crLHS78 - r_N[0]);
const double crLHS80 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS81 = r_C(0,2)*r_DN(0,0);
const double crLHS82 = crLHS81 + r_C(2,2)*r_DN(0,1);
const double crLHS83 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crLHS84 = crLHS30*crLHS83;
const double crLHS85 = -crLHS84 + r_DN(0,0);
const double crLHS86 = crLHS8*tau_c;
const double crLHS87 = crLHS86*r_DN(0,0);
const double crLHS88 = crLHS13*crLHS22;
const double crLHS89 = crLHS3*crLHS30;
const double crLHS90 = crLHS6*crLHS89;
const double crLHS91 = crLHS18*crLHS8;
const double crLHS92 = crLHS75*crLHS90;
const double crLHS93 = crLHS24*crLHS77;
const double crLHS94 = crLHS12*sigma - crLHS14*crLHS20 + crLHS18*crLHS90 + crLHS20*crLHS91 + crLHS20*crLHS92 - crLHS20*crLHS93 + crLHS6*crLHS88;
const double crLHS95 = crLHS81 + r_C(0,1)*r_DN(0,1);
const double crLHS96 = r_C(1,2)*r_DN(0,1);
const double crLHS97 = crLHS96 + r_C(2,2)*r_DN(0,0);
const double crLHS98 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crLHS99 = crLHS30*crLHS98;
const double crLHS100 = -crLHS99 + r_DN(0,1);
const double crLHS101 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0)) + r_N[3]*(r_u(3,0) - r_u_mesh(3,0));
const double crLHS102 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1)) + r_N[3]*(r_u(3,1) - r_u_mesh(3,1));
const double crLHS103 = crLHS101*r_DN(0,0) + crLHS102*r_DN(0,1);
const double crLHS104 = crLHS103 + crLHS15;
const double crLHS105 = crLHS24*crLHS40*tau_c;
const double crLHS106 = crLHS105*r_DN(0,0);
const double crLHS107 = r_DN(0,0)*r_N[1];
const double crLHS108 = crLHS76*crLHS78;
const double crLHS109 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS110 = r_C(0,2)*r_DN(1,0);
const double crLHS111 = crLHS110 + r_C(2,2)*r_DN(1,1);
const double crLHS112 = crLHS5*r_N[1];
const double crLHS113 = crLHS112*crLHS83;
const double crLHS114 = -crLHS113 + r_DN(1,0);
const double crLHS115 = crLHS112*crLHS3;
const double crLHS116 = crLHS115*crLHS15;
const double crLHS117 = crLHS116*crLHS6 + crLHS14*r_N[1];
const double crLHS118 = crLHS117 - crLHS14*crLHS37 + crLHS35*crLHS90 + crLHS37*crLHS91 + crLHS37*crLHS92 - crLHS37*crLHS93;
const double crLHS119 = crLHS110 + r_C(0,1)*r_DN(1,1);
const double crLHS120 = r_C(1,2)*r_DN(1,1);
const double crLHS121 = crLHS120 + r_C(2,2)*r_DN(1,0);
const double crLHS122 = crLHS112*crLHS98;
const double crLHS123 = -crLHS122 + r_DN(1,1);
const double crLHS124 = crLHS101*r_DN(1,0) + crLHS102*r_DN(1,1);
const double crLHS125 = crLHS124 + crLHS34;
const double crLHS126 = r_DN(0,0)*r_N[2];
const double crLHS127 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS128 = r_C(0,2)*r_DN(2,0);
const double crLHS129 = crLHS128 + r_C(2,2)*r_DN(2,1);
const double crLHS130 = crLHS5*r_N[2];
const double crLHS131 = crLHS130*crLHS83;
const double crLHS132 = -crLHS131 + r_DN(2,0);
const double crLHS133 = crLHS130*crLHS3;
const double crLHS134 = crLHS133*crLHS15;
const double crLHS135 = crLHS134*crLHS6 + crLHS14*r_N[2];
const double crLHS136 = crLHS135 - crLHS14*crLHS54 + crLHS52*crLHS90 + crLHS54*crLHS91 + crLHS54*crLHS92 - crLHS54*crLHS93;
const double crLHS137 = crLHS128 + r_C(0,1)*r_DN(2,1);
const double crLHS138 = r_C(1,2)*r_DN(2,1);
const double crLHS139 = crLHS138 + r_C(2,2)*r_DN(2,0);
const double crLHS140 = crLHS130*crLHS98;
const double crLHS141 = -crLHS140 + r_DN(2,1);
const double crLHS142 = crLHS101*r_DN(2,0) + crLHS102*r_DN(2,1);
const double crLHS143 = crLHS142 + crLHS51;
const double crLHS144 = r_DN(0,0)*r_N[3];
const double crLHS145 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS146 = r_C(0,2)*r_DN(3,0);
const double crLHS147 = crLHS146 + r_C(2,2)*r_DN(3,1);
const double crLHS148 = crLHS5*r_N[3];
const double crLHS149 = -crLHS148*crLHS83 + r_DN(3,0);
const double crLHS150 = crLHS148*crLHS3;
const double crLHS151 = crLHS15*crLHS150;
const double crLHS152 = crLHS14*r_N[3] + crLHS151*crLHS6;
const double crLHS153 = -crLHS14*crLHS69 + crLHS152 + crLHS67*crLHS90 + crLHS69*crLHS91 + crLHS69*crLHS92 - crLHS69*crLHS93;
const double crLHS154 = crLHS146 + r_C(0,1)*r_DN(3,1);
const double crLHS155 = r_C(1,2)*r_DN(3,1);
const double crLHS156 = crLHS155 + r_C(2,2)*r_DN(3,0);
const double crLHS157 = -crLHS148*crLHS98 + r_DN(3,1);
const double crLHS158 = crLHS101*r_DN(3,0) + crLHS102*r_DN(3,1);
const double crLHS159 = crLHS158 + crLHS66;
const double crLHS160 = crLHS96 + r_C(0,1)*r_DN(0,0);
const double crLHS161 = crLHS86*r_DN(0,1);
const double crLHS162 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS163 = crLHS105*r_DN(0,1);
const double crLHS164 = r_DN(0,1)*r_N[1];
const double crLHS165 = crLHS120 + r_C(0,1)*r_DN(1,0);
const double crLHS166 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS167 = r_DN(0,1)*r_N[2];
const double crLHS168 = crLHS138 + r_C(0,1)*r_DN(2,0);
const double crLHS169 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS170 = r_DN(0,1)*r_N[3];
const double crLHS171 = crLHS155 + r_C(0,1)*r_DN(3,0);
const double crLHS172 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS173 = crLHS13*crLHS40;
const double crLHS174 = -crLHS19 + dp_th_dt*r_N[0];
const double crLHS175 = tau_t*1.0/(crLHS4*crLHS4*crLHS4);
const double crLHS176 = crLHS175*crLHS77;
const double crLHS177 = crLHS41*crLHS84;
const double crLHS178 = crLHS41*crLHS99;
const double crLHS179 = -crLHS36 + dp_th_dt*r_N[1];
const double crLHS180 = crLHS116 + crLHS26*kappa + crLHS27*kappa - crLHS5*dp_th_dt*r_N[0]*r_N[1];
const double crLHS181 = crLHS57*crLHS84;
const double crLHS182 = crLHS57*crLHS99;
const double crLHS183 = -crLHS53 + dp_th_dt*r_N[2];
const double crLHS184 = crLHS134 + crLHS44*kappa + crLHS45*kappa - crLHS5*dp_th_dt*r_N[0]*r_N[2];
const double crLHS185 = crLHS72*crLHS84;
const double crLHS186 = crLHS72*crLHS99;
const double crLHS187 = -crLHS68 + dp_th_dt*r_N[3];
const double crLHS188 = crLHS151 - crLHS5*dp_th_dt*r_N[0]*r_N[3] + crLHS59*kappa + crLHS60*kappa;
const double crLHS189 = r_DN(1,0)*r_DN(1,0);
const double crLHS190 = r_DN(1,1)*r_DN(1,1);
const double crLHS191 = r_N[1]*r_N[1];
const double crLHS192 = crLHS191*crLHS5;
const double crLHS193 = r_DN(1,0)*r_DN(2,0);
const double crLHS194 = r_DN(1,1)*r_DN(2,1);
const double crLHS195 = crLHS10*(crLHS193 + crLHS194);
const double crLHS196 = r_DN(2,0)*r_N[1];
const double crLHS197 = crLHS112*r_N[2];
const double crLHS198 = -crLHS11*crLHS197;
const double crLHS199 = r_DN(2,1)*r_N[1];
const double crLHS200 = -crLHS197*crLHS21;
const double crLHS201 = crLHS24*crLHS34;
const double crLHS202 = -crLHS201*crLHS57;
const double crLHS203 = r_DN(1,0)*r_DN(3,0);
const double crLHS204 = r_DN(1,1)*r_DN(3,1);
const double crLHS205 = crLHS10*(crLHS203 + crLHS204);
const double crLHS206 = r_DN(3,0)*r_N[1];
const double crLHS207 = crLHS112*r_N[3];
const double crLHS208 = -crLHS11*crLHS207;
const double crLHS209 = r_DN(3,1)*r_N[1];
const double crLHS210 = -crLHS207*crLHS21;
const double crLHS211 = -crLHS201*crLHS72;
const double crLHS212 = crLHS33*tau_u;
const double crLHS213 = crLHS86*r_DN(1,0);
const double crLHS214 = crLHS115*crLHS6;
const double crLHS215 = crLHS35*crLHS8;
const double crLHS216 = crLHS214*crLHS75;
const double crLHS217 = crLHS76*r_N[1];
const double crLHS218 = crLHS217*crLHS24;
const double crLHS219 = crLHS117 + crLHS18*crLHS214 + crLHS20*crLHS215 + crLHS20*crLHS216 - crLHS20*crLHS218 - crLHS20*crLHS33;
const double crLHS220 = crLHS105*r_DN(1,0);
const double crLHS221 = gauss_weight*(crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_N[1]*tau_u - crLHS212 - crLHS217*crLHS78 - r_N[1]);
const double crLHS222 = crLHS192*crLHS22;
const double crLHS223 = crLHS191*sigma + crLHS214*crLHS35 + crLHS215*crLHS37 + crLHS216*crLHS37 - crLHS218*crLHS37 + crLHS222*crLHS6 - crLHS33*crLHS37;
const double crLHS224 = r_DN(1,0)*r_N[2];
const double crLHS225 = crLHS133*crLHS34;
const double crLHS226 = crLHS225*crLHS6 + crLHS33*r_N[2];
const double crLHS227 = crLHS214*crLHS52 + crLHS215*crLHS54 + crLHS216*crLHS54 - crLHS218*crLHS54 + crLHS226 - crLHS33*crLHS54;
const double crLHS228 = r_DN(1,0)*r_N[3];
const double crLHS229 = crLHS150*crLHS34;
const double crLHS230 = crLHS229*crLHS6 + crLHS33*r_N[3];
const double crLHS231 = crLHS214*crLHS67 + crLHS215*crLHS69 + crLHS216*crLHS69 - crLHS218*crLHS69 + crLHS230 - crLHS33*crLHS69;
const double crLHS232 = crLHS86*r_DN(1,1);
const double crLHS233 = crLHS105*r_DN(1,1);
const double crLHS234 = r_DN(1,1)*r_N[2];
const double crLHS235 = r_DN(1,1)*r_N[3];
const double crLHS236 = crLHS175*crLHS217;
const double crLHS237 = crLHS192*crLHS40;
const double crLHS238 = crLHS113*crLHS57;
const double crLHS239 = crLHS122*crLHS57;
const double crLHS240 = crLHS193*kappa + crLHS194*kappa + crLHS225 - crLHS5*dp_th_dt*r_N[1]*r_N[2];
const double crLHS241 = crLHS113*crLHS72;
const double crLHS242 = crLHS122*crLHS72;
const double crLHS243 = crLHS203*kappa + crLHS204*kappa + crLHS229 - crLHS5*dp_th_dt*r_N[1]*r_N[3];
const double crLHS244 = r_DN(2,0)*r_DN(2,0);
const double crLHS245 = r_DN(2,1)*r_DN(2,1);
const double crLHS246 = r_N[2]*r_N[2];
const double crLHS247 = crLHS246*crLHS5;
const double crLHS248 = r_DN(2,0)*r_DN(3,0);
const double crLHS249 = r_DN(2,1)*r_DN(3,1);
const double crLHS250 = crLHS10*(crLHS248 + crLHS249);
const double crLHS251 = r_DN(3,0)*r_N[2];
const double crLHS252 = crLHS130*r_N[3];
const double crLHS253 = -crLHS11*crLHS252;
const double crLHS254 = r_DN(3,1)*r_N[2];
const double crLHS255 = -crLHS21*crLHS252;
const double crLHS256 = -crLHS24*crLHS51*crLHS72;
const double crLHS257 = crLHS50*tau_u;
const double crLHS258 = crLHS86*r_DN(2,0);
const double crLHS259 = crLHS133*crLHS6;
const double crLHS260 = crLHS52*crLHS8;
const double crLHS261 = crLHS259*crLHS75;
const double crLHS262 = crLHS76*r_N[2];
const double crLHS263 = crLHS24*crLHS262;
const double crLHS264 = crLHS135 + crLHS18*crLHS259 + crLHS20*crLHS260 + crLHS20*crLHS261 - crLHS20*crLHS263 - crLHS20*crLHS50;
const double crLHS265 = crLHS105*r_DN(2,0);
const double crLHS266 = crLHS226 + crLHS259*crLHS35 + crLHS260*crLHS37 + crLHS261*crLHS37 - crLHS263*crLHS37 - crLHS37*crLHS50;
const double crLHS267 = gauss_weight*(crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_N[2]*tau_u - crLHS257 - crLHS262*crLHS78 - r_N[2]);
const double crLHS268 = crLHS22*crLHS247;
const double crLHS269 = crLHS246*sigma + crLHS259*crLHS52 + crLHS260*crLHS54 + crLHS261*crLHS54 - crLHS263*crLHS54 + crLHS268*crLHS6 - crLHS50*crLHS54;
const double crLHS270 = r_DN(2,0)*r_N[3];
const double crLHS271 = crLHS150*crLHS51;
const double crLHS272 = crLHS271*crLHS6 + crLHS50*r_N[3];
const double crLHS273 = crLHS259*crLHS67 + crLHS260*crLHS69 + crLHS261*crLHS69 - crLHS263*crLHS69 + crLHS272 - crLHS50*crLHS69;
const double crLHS274 = crLHS86*r_DN(2,1);
const double crLHS275 = crLHS105*r_DN(2,1);
const double crLHS276 = r_DN(2,1)*r_N[3];
const double crLHS277 = crLHS175*crLHS262;
const double crLHS278 = crLHS247*crLHS40;
const double crLHS279 = crLHS131*crLHS72;
const double crLHS280 = crLHS140*crLHS72;
const double crLHS281 = crLHS248*kappa + crLHS249*kappa + crLHS271 - crLHS5*dp_th_dt*r_N[2]*r_N[3];
const double crLHS282 = r_DN(3,0)*r_DN(3,0);
const double crLHS283 = r_DN(3,1)*r_DN(3,1);
const double crLHS284 = r_N[3]*r_N[3];
const double crLHS285 = crLHS284*crLHS5;
const double crLHS286 = crLHS65*tau_u;
const double crLHS287 = crLHS86*r_DN(3,0);
const double crLHS288 = crLHS150*crLHS6;
const double crLHS289 = crLHS67*crLHS8;
const double crLHS290 = crLHS288*crLHS75;
const double crLHS291 = crLHS76*r_N[3];
const double crLHS292 = crLHS24*crLHS291;
const double crLHS293 = crLHS152 + crLHS18*crLHS288 + crLHS20*crLHS289 + crLHS20*crLHS290 - crLHS20*crLHS292 - crLHS20*crLHS65;
const double crLHS294 = crLHS105*r_DN(3,0);
const double crLHS295 = crLHS230 + crLHS288*crLHS35 + crLHS289*crLHS37 + crLHS290*crLHS37 - crLHS292*crLHS37 - crLHS37*crLHS65;
const double crLHS296 = crLHS272 + crLHS288*crLHS52 + crLHS289*crLHS54 + crLHS290*crLHS54 - crLHS292*crLHS54 - crLHS54*crLHS65;
const double crLHS297 = gauss_weight*(crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_N[3]*tau_u - crLHS286 - crLHS291*crLHS78 - r_N[3]);
const double crLHS298 = crLHS22*crLHS285;
const double crLHS299 = crLHS284*sigma + crLHS288*crLHS67 + crLHS289*crLHS69 + crLHS290*crLHS69 - crLHS292*crLHS69 + crLHS298*crLHS6 - crLHS65*crLHS69;
const double crLHS300 = crLHS86*r_DN(3,1);
const double crLHS301 = crLHS105*r_DN(3,1);
const double crLHS302 = crLHS175*crLHS291;
const double crLHS303 = crLHS285*crLHS40;
rLHS(0,0)+=crLHS10*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS9*(-crLHS11*crLHS13 + crLHS20*r_DN(0,0) + r_DN(0,0)*r_N[0]);
rLHS(0,2)+=crLHS9*(-crLHS13*crLHS21 + crLHS20*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS12*crLHS25;
rLHS(0,4)+=crLHS28;
rLHS(0,5)+=crLHS9*(crLHS29 + crLHS32 + crLHS37*r_DN(0,0));
rLHS(0,6)+=crLHS9*(crLHS37*r_DN(0,1) + crLHS38 + crLHS39);
rLHS(0,7)+=crLHS43;
rLHS(0,8)+=crLHS46;
rLHS(0,9)+=crLHS9*(crLHS47 + crLHS49 + crLHS54*r_DN(0,0));
rLHS(0,10)+=crLHS9*(crLHS54*r_DN(0,1) + crLHS55 + crLHS56);
rLHS(0,11)+=crLHS58;
rLHS(0,12)+=crLHS61;
rLHS(0,13)+=crLHS9*(crLHS62 + crLHS64 + crLHS69*r_DN(0,0));
rLHS(0,14)+=crLHS9*(crLHS69*r_DN(0,1) + crLHS70 + crLHS71);
rLHS(0,15)+=crLHS73;
rLHS(1,0)+=crLHS79*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS80*r_DN(0,0) + crLHS82*r_DN(0,1) + crLHS85*crLHS87 + crLHS94);
rLHS(1,2)+=gauss_weight*(crLHS100*crLHS87 + crLHS95*r_DN(0,0) + crLHS97*r_DN(0,1));
rLHS(1,3)+=-crLHS104*crLHS106;
rLHS(1,4)+=gauss_weight*(-crLHS107 - crLHS108*crLHS29 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS74*r_DN(1,0));
rLHS(1,5)+=gauss_weight*(crLHS109*r_DN(0,0) + crLHS111*r_DN(0,1) + crLHS114*crLHS87 + crLHS118);
rLHS(1,6)+=gauss_weight*(crLHS119*r_DN(0,0) + crLHS121*r_DN(0,1) + crLHS123*crLHS87);
rLHS(1,7)+=-crLHS106*crLHS125;
rLHS(1,8)+=gauss_weight*(-crLHS108*crLHS47 - crLHS126 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS74*r_DN(2,0));
rLHS(1,9)+=gauss_weight*(crLHS127*r_DN(0,0) + crLHS129*r_DN(0,1) + crLHS132*crLHS87 + crLHS136);
rLHS(1,10)+=gauss_weight*(crLHS137*r_DN(0,0) + crLHS139*r_DN(0,1) + crLHS141*crLHS87);
rLHS(1,11)+=-crLHS106*crLHS143;
rLHS(1,12)+=gauss_weight*(-crLHS108*crLHS62 - crLHS144 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,0)*r_N[0]*tau_u - crLHS74*r_DN(3,0));
rLHS(1,13)+=gauss_weight*(crLHS145*r_DN(0,0) + crLHS147*r_DN(0,1) + crLHS149*crLHS87 + crLHS153);
rLHS(1,14)+=gauss_weight*(crLHS154*r_DN(0,0) + crLHS156*r_DN(0,1) + crLHS157*crLHS87);
rLHS(1,15)+=-crLHS106*crLHS159;
rLHS(2,0)+=crLHS79*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS160*r_DN(0,1) + crLHS161*crLHS85 + crLHS82*r_DN(0,0));
rLHS(2,2)+=gauss_weight*(crLHS100*crLHS161 + crLHS162*r_DN(0,1) + crLHS94 + crLHS97*r_DN(0,0));
rLHS(2,3)+=-crLHS104*crLHS163;
rLHS(2,4)+=gauss_weight*(-crLHS108*crLHS38 - crLHS164 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS74*r_DN(1,1));
rLHS(2,5)+=gauss_weight*(crLHS111*r_DN(0,0) + crLHS114*crLHS161 + crLHS165*r_DN(0,1));
rLHS(2,6)+=gauss_weight*(crLHS118 + crLHS121*r_DN(0,0) + crLHS123*crLHS161 + crLHS166*r_DN(0,1));
rLHS(2,7)+=-crLHS125*crLHS163;
rLHS(2,8)+=gauss_weight*(-crLHS108*crLHS55 - crLHS167 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS74*r_DN(2,1));
rLHS(2,9)+=gauss_weight*(crLHS129*r_DN(0,0) + crLHS132*crLHS161 + crLHS168*r_DN(0,1));
rLHS(2,10)+=gauss_weight*(crLHS136 + crLHS139*r_DN(0,0) + crLHS141*crLHS161 + crLHS169*r_DN(0,1));
rLHS(2,11)+=-crLHS143*crLHS163;
rLHS(2,12)+=gauss_weight*(-crLHS108*crLHS70 - crLHS170 + crLHS18*crLHS2*crLHS5*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,1)*r_N[0]*tau_u - crLHS74*r_DN(3,1));
rLHS(2,13)+=gauss_weight*(crLHS147*r_DN(0,0) + crLHS149*crLHS161 + crLHS171*r_DN(0,1));
rLHS(2,14)+=gauss_weight*(crLHS153 + crLHS156*r_DN(0,0) + crLHS157*crLHS161 + crLHS172*r_DN(0,1));
rLHS(2,15)+=-crLHS159*crLHS163;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS173*crLHS83;
rLHS(3,2)+=crLHS173*crLHS98;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa - crLHS103*crLHS89 + crLHS12*crLHS5*dp_th_dt - crLHS174*crLHS176 + crLHS174*crLHS18*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS174*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS174*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS88);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS177;
rLHS(3,6)+=crLHS178;
rLHS(3,7)+=-gauss_weight*(-crLHS124*crLHS89 - crLHS176*crLHS179 + crLHS179*crLHS18*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS179*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS179*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS180);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS181;
rLHS(3,10)+=crLHS182;
rLHS(3,11)+=-gauss_weight*(-crLHS142*crLHS89 - crLHS176*crLHS183 + crLHS18*crLHS183*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS183*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS183*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS184);
rLHS(3,12)+=0;
rLHS(3,13)+=crLHS185;
rLHS(3,14)+=crLHS186;
rLHS(3,15)+=-gauss_weight*(-crLHS158*crLHS89 - crLHS176*crLHS187 + crLHS18*crLHS187*crLHS2*crLHS23*gamma*p_th*tau_t + crLHS187*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS187*crLHS23*dp_th_dt*r_N[0]*tau_t - crLHS188);
rLHS(4,0)+=crLHS28;
rLHS(4,1)+=crLHS9*(crLHS107 + crLHS20*r_DN(1,0) + crLHS32);
rLHS(4,2)+=crLHS9*(crLHS164 + crLHS20*r_DN(1,1) + crLHS39);
rLHS(4,3)+=crLHS43;
rLHS(4,4)+=crLHS10*(crLHS189 + crLHS190);
rLHS(4,5)+=crLHS9*(-crLHS11*crLHS192 + crLHS37*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS9*(-crLHS192*crLHS21 + crLHS37*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS191*crLHS25;
rLHS(4,8)+=crLHS195;
rLHS(4,9)+=crLHS9*(crLHS196 + crLHS198 + crLHS54*r_DN(1,0));
rLHS(4,10)+=crLHS9*(crLHS199 + crLHS200 + crLHS54*r_DN(1,1));
rLHS(4,11)+=crLHS202;
rLHS(4,12)+=crLHS205;
rLHS(4,13)+=crLHS9*(crLHS206 + crLHS208 + crLHS69*r_DN(1,0));
rLHS(4,14)+=crLHS9*(crLHS209 + crLHS210 + crLHS69*r_DN(1,1));
rLHS(4,15)+=crLHS211;
rLHS(5,0)+=gauss_weight*(-crLHS107*crLHS108 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS212*r_DN(0,0) - crLHS29);
rLHS(5,1)+=gauss_weight*(crLHS213*crLHS85 + crLHS219 + crLHS80*r_DN(1,0) + crLHS82*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS100*crLHS213 + crLHS95*r_DN(1,0) + crLHS97*r_DN(1,1));
rLHS(5,3)+=-crLHS104*crLHS220;
rLHS(5,4)+=crLHS221*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS109*r_DN(1,0) + crLHS111*r_DN(1,1) + crLHS114*crLHS213 + crLHS223);
rLHS(5,6)+=gauss_weight*(crLHS119*r_DN(1,0) + crLHS121*r_DN(1,1) + crLHS123*crLHS213);
rLHS(5,7)+=-crLHS125*crLHS220;
rLHS(5,8)+=gauss_weight*(-crLHS108*crLHS196 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u - crLHS212*r_DN(2,0) - crLHS224);
rLHS(5,9)+=gauss_weight*(crLHS127*r_DN(1,0) + crLHS129*r_DN(1,1) + crLHS132*crLHS213 + crLHS227);
rLHS(5,10)+=gauss_weight*(crLHS137*r_DN(1,0) + crLHS139*r_DN(1,1) + crLHS141*crLHS213);
rLHS(5,11)+=-crLHS143*crLHS220;
rLHS(5,12)+=gauss_weight*(-crLHS108*crLHS206 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,0)*r_N[1]*tau_u - crLHS212*r_DN(3,0) - crLHS228);
rLHS(5,13)+=gauss_weight*(crLHS145*r_DN(1,0) + crLHS147*r_DN(1,1) + crLHS149*crLHS213 + crLHS231);
rLHS(5,14)+=gauss_weight*(crLHS154*r_DN(1,0) + crLHS156*r_DN(1,1) + crLHS157*crLHS213);
rLHS(5,15)+=-crLHS159*crLHS220;
rLHS(6,0)+=gauss_weight*(-crLHS108*crLHS164 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS212*r_DN(0,1) - crLHS38);
rLHS(6,1)+=gauss_weight*(crLHS160*r_DN(1,1) + crLHS232*crLHS85 + crLHS82*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS100*crLHS232 + crLHS162*r_DN(1,1) + crLHS219 + crLHS97*r_DN(1,0));
rLHS(6,3)+=-crLHS104*crLHS233;
rLHS(6,4)+=crLHS221*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS111*r_DN(1,0) + crLHS114*crLHS232 + crLHS165*r_DN(1,1));
rLHS(6,6)+=gauss_weight*(crLHS121*r_DN(1,0) + crLHS123*crLHS232 + crLHS166*r_DN(1,1) + crLHS223);
rLHS(6,7)+=-crLHS125*crLHS233;
rLHS(6,8)+=gauss_weight*(-crLHS108*crLHS199 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u - crLHS212*r_DN(2,1) - crLHS234);
rLHS(6,9)+=gauss_weight*(crLHS129*r_DN(1,0) + crLHS132*crLHS232 + crLHS168*r_DN(1,1));
rLHS(6,10)+=gauss_weight*(crLHS139*r_DN(1,0) + crLHS141*crLHS232 + crLHS169*r_DN(1,1) + crLHS227);
rLHS(6,11)+=-crLHS143*crLHS233;
rLHS(6,12)+=gauss_weight*(-crLHS108*crLHS209 + crLHS2*crLHS35*crLHS5*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,1)*r_N[1]*tau_u - crLHS212*r_DN(3,1) - crLHS235);
rLHS(6,13)+=gauss_weight*(crLHS147*r_DN(1,0) + crLHS149*crLHS232 + crLHS171*r_DN(1,1));
rLHS(6,14)+=gauss_weight*(crLHS156*r_DN(1,0) + crLHS157*crLHS232 + crLHS172*r_DN(1,1) + crLHS231);
rLHS(6,15)+=-crLHS159*crLHS233;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS177;
rLHS(7,2)+=crLHS178;
rLHS(7,3)+=-gauss_weight*(-crLHS103*crLHS115 + crLHS174*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS174*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS174*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS174*crLHS236 - crLHS180);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS237*crLHS83;
rLHS(7,6)+=crLHS237*crLHS98;
rLHS(7,7)+=-gauss_weight*(-crLHS115*crLHS124 + crLHS179*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS179*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS179*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS179*crLHS236 - crLHS189*kappa - crLHS190*kappa + crLHS191*crLHS5*dp_th_dt - crLHS222);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS238;
rLHS(7,10)+=crLHS239;
rLHS(7,11)+=-gauss_weight*(-crLHS115*crLHS142 + crLHS183*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS183*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS183*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS183*crLHS236 - crLHS240);
rLHS(7,12)+=0;
rLHS(7,13)+=crLHS241;
rLHS(7,14)+=crLHS242;
rLHS(7,15)+=-gauss_weight*(-crLHS115*crLHS158 + crLHS187*crLHS2*crLHS23*crLHS35*gamma*p_th*tau_t + crLHS187*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS187*crLHS23*dp_th_dt*r_N[1]*tau_t - crLHS187*crLHS236 - crLHS243);
rLHS(8,0)+=crLHS46;
rLHS(8,1)+=crLHS9*(crLHS126 + crLHS20*r_DN(2,0) + crLHS49);
rLHS(8,2)+=crLHS9*(crLHS167 + crLHS20*r_DN(2,1) + crLHS56);
rLHS(8,3)+=crLHS58;
rLHS(8,4)+=crLHS195;
rLHS(8,5)+=crLHS9*(crLHS198 + crLHS224 + crLHS37*r_DN(2,0));
rLHS(8,6)+=crLHS9*(crLHS200 + crLHS234 + crLHS37*r_DN(2,1));
rLHS(8,7)+=crLHS202;
rLHS(8,8)+=crLHS10*(crLHS244 + crLHS245);
rLHS(8,9)+=crLHS9*(-crLHS11*crLHS247 + crLHS54*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS9*(-crLHS21*crLHS247 + crLHS54*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS246*crLHS25;
rLHS(8,12)+=crLHS250;
rLHS(8,13)+=crLHS9*(crLHS251 + crLHS253 + crLHS69*r_DN(2,0));
rLHS(8,14)+=crLHS9*(crLHS254 + crLHS255 + crLHS69*r_DN(2,1));
rLHS(8,15)+=crLHS256;
rLHS(9,0)+=gauss_weight*(-crLHS108*crLHS126 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS257*r_DN(0,0) - crLHS47);
rLHS(9,1)+=gauss_weight*(crLHS258*crLHS85 + crLHS264 + crLHS80*r_DN(2,0) + crLHS82*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS100*crLHS258 + crLHS95*r_DN(2,0) + crLHS97*r_DN(2,1));
rLHS(9,3)+=-crLHS104*crLHS265;
rLHS(9,4)+=gauss_weight*(-crLHS108*crLHS224 - crLHS196 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u - crLHS257*r_DN(1,0));
rLHS(9,5)+=gauss_weight*(crLHS109*r_DN(2,0) + crLHS111*r_DN(2,1) + crLHS114*crLHS258 + crLHS266);
rLHS(9,6)+=gauss_weight*(crLHS119*r_DN(2,0) + crLHS121*r_DN(2,1) + crLHS123*crLHS258);
rLHS(9,7)+=-crLHS125*crLHS265;
rLHS(9,8)+=crLHS267*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS127*r_DN(2,0) + crLHS129*r_DN(2,1) + crLHS132*crLHS258 + crLHS269);
rLHS(9,10)+=gauss_weight*(crLHS137*r_DN(2,0) + crLHS139*r_DN(2,1) + crLHS141*crLHS258);
rLHS(9,11)+=-crLHS143*crLHS265;
rLHS(9,12)+=gauss_weight*(-crLHS108*crLHS251 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,0)*r_N[2]*tau_u - crLHS257*r_DN(3,0) - crLHS270);
rLHS(9,13)+=gauss_weight*(crLHS145*r_DN(2,0) + crLHS147*r_DN(2,1) + crLHS149*crLHS258 + crLHS273);
rLHS(9,14)+=gauss_weight*(crLHS154*r_DN(2,0) + crLHS156*r_DN(2,1) + crLHS157*crLHS258);
rLHS(9,15)+=-crLHS159*crLHS265;
rLHS(10,0)+=gauss_weight*(-crLHS108*crLHS167 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS257*r_DN(0,1) - crLHS55);
rLHS(10,1)+=gauss_weight*(crLHS160*r_DN(2,1) + crLHS274*crLHS85 + crLHS82*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS100*crLHS274 + crLHS162*r_DN(2,1) + crLHS264 + crLHS97*r_DN(2,0));
rLHS(10,3)+=-crLHS104*crLHS275;
rLHS(10,4)+=gauss_weight*(-crLHS108*crLHS234 - crLHS199 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u - crLHS257*r_DN(1,1));
rLHS(10,5)+=gauss_weight*(crLHS111*r_DN(2,0) + crLHS114*crLHS274 + crLHS165*r_DN(2,1));
rLHS(10,6)+=gauss_weight*(crLHS121*r_DN(2,0) + crLHS123*crLHS274 + crLHS166*r_DN(2,1) + crLHS266);
rLHS(10,7)+=-crLHS125*crLHS275;
rLHS(10,8)+=crLHS267*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS129*r_DN(2,0) + crLHS132*crLHS274 + crLHS168*r_DN(2,1));
rLHS(10,10)+=gauss_weight*(crLHS139*r_DN(2,0) + crLHS141*crLHS274 + crLHS169*r_DN(2,1) + crLHS269);
rLHS(10,11)+=-crLHS143*crLHS275;
rLHS(10,12)+=gauss_weight*(-crLHS108*crLHS254 + crLHS2*crLHS5*crLHS52*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(3,1)*r_N[2]*tau_u - crLHS257*r_DN(3,1) - crLHS276);
rLHS(10,13)+=gauss_weight*(crLHS147*r_DN(2,0) + crLHS149*crLHS274 + crLHS171*r_DN(2,1));
rLHS(10,14)+=gauss_weight*(crLHS156*r_DN(2,0) + crLHS157*crLHS274 + crLHS172*r_DN(2,1) + crLHS273);
rLHS(10,15)+=-crLHS159*crLHS275;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS181;
rLHS(11,2)+=crLHS182;
rLHS(11,3)+=-gauss_weight*(-crLHS103*crLHS133 + crLHS174*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS174*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS174*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS174*crLHS277 - crLHS184);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS238;
rLHS(11,6)+=crLHS239;
rLHS(11,7)+=-gauss_weight*(-crLHS124*crLHS133 + crLHS179*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS179*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS179*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS179*crLHS277 - crLHS240);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS278*crLHS83;
rLHS(11,10)+=crLHS278*crLHS98;
rLHS(11,11)+=-gauss_weight*(-crLHS133*crLHS142 + crLHS183*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS183*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS183*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS183*crLHS277 - crLHS244*kappa - crLHS245*kappa + crLHS246*crLHS5*dp_th_dt - crLHS268);
rLHS(11,12)+=0;
rLHS(11,13)+=crLHS279;
rLHS(11,14)+=crLHS280;
rLHS(11,15)+=-gauss_weight*(-crLHS133*crLHS158 + crLHS187*crLHS2*crLHS23*crLHS52*gamma*p_th*tau_t + crLHS187*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS187*crLHS23*dp_th_dt*r_N[2]*tau_t - crLHS187*crLHS277 - crLHS281);
rLHS(12,0)+=crLHS61;
rLHS(12,1)+=crLHS9*(crLHS144 + crLHS20*r_DN(3,0) + crLHS64);
rLHS(12,2)+=crLHS9*(crLHS170 + crLHS20*r_DN(3,1) + crLHS71);
rLHS(12,3)+=crLHS73;
rLHS(12,4)+=crLHS205;
rLHS(12,5)+=crLHS9*(crLHS208 + crLHS228 + crLHS37*r_DN(3,0));
rLHS(12,6)+=crLHS9*(crLHS210 + crLHS235 + crLHS37*r_DN(3,1));
rLHS(12,7)+=crLHS211;
rLHS(12,8)+=crLHS250;
rLHS(12,9)+=crLHS9*(crLHS253 + crLHS270 + crLHS54*r_DN(3,0));
rLHS(12,10)+=crLHS9*(crLHS255 + crLHS276 + crLHS54*r_DN(3,1));
rLHS(12,11)+=crLHS256;
rLHS(12,12)+=crLHS10*(crLHS282 + crLHS283);
rLHS(12,13)+=crLHS9*(-crLHS11*crLHS285 + crLHS69*r_DN(3,0) + r_DN(3,0)*r_N[3]);
rLHS(12,14)+=crLHS9*(-crLHS21*crLHS285 + crLHS69*r_DN(3,1) + r_DN(3,1)*r_N[3]);
rLHS(12,15)+=-crLHS25*crLHS284;
rLHS(13,0)+=gauss_weight*(-crLHS108*crLHS144 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(0,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,0)*r_N[3]*tau_u - crLHS286*r_DN(0,0) - crLHS62);
rLHS(13,1)+=gauss_weight*(crLHS287*crLHS85 + crLHS293 + crLHS80*r_DN(3,0) + crLHS82*r_DN(3,1));
rLHS(13,2)+=gauss_weight*(crLHS100*crLHS287 + crLHS95*r_DN(3,0) + crLHS97*r_DN(3,1));
rLHS(13,3)+=-crLHS104*crLHS294;
rLHS(13,4)+=gauss_weight*(-crLHS108*crLHS228 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(1,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,0)*r_N[3]*tau_u - crLHS206 - crLHS286*r_DN(1,0));
rLHS(13,5)+=gauss_weight*(crLHS109*r_DN(3,0) + crLHS111*r_DN(3,1) + crLHS114*crLHS287 + crLHS295);
rLHS(13,6)+=gauss_weight*(crLHS119*r_DN(3,0) + crLHS121*r_DN(3,1) + crLHS123*crLHS287);
rLHS(13,7)+=-crLHS125*crLHS294;
rLHS(13,8)+=gauss_weight*(-crLHS108*crLHS270 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(2,0)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,0)*r_N[3]*tau_u - crLHS251 - crLHS286*r_DN(2,0));
rLHS(13,9)+=gauss_weight*(crLHS127*r_DN(3,0) + crLHS129*r_DN(3,1) + crLHS132*crLHS287 + crLHS296);
rLHS(13,10)+=gauss_weight*(crLHS137*r_DN(3,0) + crLHS139*r_DN(3,1) + crLHS141*crLHS287);
rLHS(13,11)+=-crLHS143*crLHS294;
rLHS(13,12)+=crLHS297*r_DN(3,0);
rLHS(13,13)+=gauss_weight*(crLHS145*r_DN(3,0) + crLHS147*r_DN(3,1) + crLHS149*crLHS287 + crLHS299);
rLHS(13,14)+=gauss_weight*(crLHS154*r_DN(3,0) + crLHS156*r_DN(3,1) + crLHS157*crLHS287);
rLHS(13,15)+=-crLHS159*crLHS294;
rLHS(14,0)+=gauss_weight*(-crLHS108*crLHS170 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(0,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(0,1)*r_N[3]*tau_u - crLHS286*r_DN(0,1) - crLHS70);
rLHS(14,1)+=gauss_weight*(crLHS160*r_DN(3,1) + crLHS300*crLHS85 + crLHS82*r_DN(3,0));
rLHS(14,2)+=gauss_weight*(crLHS100*crLHS300 + crLHS162*r_DN(3,1) + crLHS293 + crLHS97*r_DN(3,0));
rLHS(14,3)+=-crLHS104*crLHS301;
rLHS(14,4)+=gauss_weight*(-crLHS108*crLHS235 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(1,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(1,1)*r_N[3]*tau_u - crLHS209 - crLHS286*r_DN(1,1));
rLHS(14,5)+=gauss_weight*(crLHS111*r_DN(3,0) + crLHS114*crLHS300 + crLHS165*r_DN(3,1));
rLHS(14,6)+=gauss_weight*(crLHS121*r_DN(3,0) + crLHS123*crLHS300 + crLHS166*r_DN(3,1) + crLHS295);
rLHS(14,7)+=-crLHS125*crLHS301;
rLHS(14,8)+=gauss_weight*(-crLHS108*crLHS276 + crLHS2*crLHS5*crLHS6*crLHS67*gamma*p_th*r_DN(2,1)*tau_u + crLHS2*crLHS5*crLHS6*crLHS75*gamma*p_th*r_DN(2,1)*r_N[3]*tau_u - crLHS254 - crLHS286*r_DN(2,1));
rLHS(14,9)+=gauss_weight*(crLHS129*r_DN(3,0) + crLHS132*crLHS300 + crLHS168*r_DN(3,1));
rLHS(14,10)+=gauss_weight*(crLHS139*r_DN(3,0) + crLHS141*crLHS300 + crLHS169*r_DN(3,1) + crLHS296);
rLHS(14,11)+=-crLHS143*crLHS301;
rLHS(14,12)+=crLHS297*r_DN(3,1);
rLHS(14,13)+=gauss_weight*(crLHS147*r_DN(3,0) + crLHS149*crLHS300 + crLHS171*r_DN(3,1));
rLHS(14,14)+=gauss_weight*(crLHS156*r_DN(3,0) + crLHS157*crLHS300 + crLHS172*r_DN(3,1) + crLHS299);
rLHS(14,15)+=-crLHS159*crLHS301;
rLHS(15,0)+=0;
rLHS(15,1)+=crLHS185;
rLHS(15,2)+=crLHS186;
rLHS(15,3)+=-gauss_weight*(-crLHS103*crLHS150 + crLHS174*crLHS2*crLHS23*crLHS67*gamma*p_th*tau_t + crLHS174*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS174*crLHS23*dp_th_dt*r_N[3]*tau_t - crLHS174*crLHS302 - crLHS188);
rLHS(15,4)+=0;
rLHS(15,5)+=crLHS241;
rLHS(15,6)+=crLHS242;
rLHS(15,7)+=-gauss_weight*(-crLHS124*crLHS150 + crLHS179*crLHS2*crLHS23*crLHS67*gamma*p_th*tau_t + crLHS179*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS179*crLHS23*dp_th_dt*r_N[3]*tau_t - crLHS179*crLHS302 - crLHS243);
rLHS(15,8)+=0;
rLHS(15,9)+=crLHS279;
rLHS(15,10)+=crLHS280;
rLHS(15,11)+=-gauss_weight*(-crLHS142*crLHS150 + crLHS183*crLHS2*crLHS23*crLHS67*gamma*p_th*tau_t + crLHS183*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS183*crLHS23*dp_th_dt*r_N[3]*tau_t - crLHS183*crLHS302 - crLHS281);
rLHS(15,12)+=0;
rLHS(15,13)+=crLHS303*crLHS83;
rLHS(15,14)+=crLHS303*crLHS98;
rLHS(15,15)+=-gauss_weight*(-crLHS150*crLHS158 + crLHS187*crLHS2*crLHS23*crLHS67*gamma*p_th*tau_t + crLHS187*crLHS2*crLHS23*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS187*crLHS23*dp_th_dt*r_N[3]*tau_t - crLHS187*crLHS302 - crLHS282*kappa - crLHS283*kappa + crLHS284*crLHS5*dp_th_dt - crLHS298);

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double sigma = rData.Resistance;
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
    const auto& r_u_sol_frac = rData.SolidFractionVelocity;

    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 3> lin_u_conv = r_u - r_u_mesh;

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

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]);
const double crRHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crRHS5 = 1.0/crRHS4;
const double crRHS6 = crRHS5*p_th;
const double crRHS7 = crRHS3*crRHS6;
const double crRHS8 = -crRHS7 + dp_th_dt;
const double crRHS9 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crRHS10 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0));
const double crRHS11 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crRHS12 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1));
const double crRHS13 = crRHS6*(crRHS10*crRHS9 + crRHS11*crRHS12);
const double crRHS14 = sigma*(r_N[0]*r_u(0,0) - r_N[0]*r_u_sol_frac(0,0) + r_N[1]*r_u(1,0) - r_N[1]*r_u_sol_frac(1,0) + r_N[2]*r_u(2,0) - r_N[2]*r_u_sol_frac(2,0));
const double crRHS15 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2];
const double crRHS16 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2];
const double crRHS17 = crRHS0*crRHS15 + crRHS16*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS18 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS19 = 1.0/c_p;
const double crRHS20 = gamma/(gamma - 1.0);
const double crRHS21 = crRHS20*crRHS6;
const double crRHS22 = crRHS19*crRHS21;
const double crRHS23 = crRHS14 - crRHS22*(-crRHS17 - crRHS18 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS24 = p_th*tau_u;
const double crRHS25 = crRHS23*crRHS24;
const double crRHS26 = sigma*(r_N[0]*r_u(0,1) - r_N[0]*r_u_sol_frac(0,1) + r_N[1]*r_u(1,1) - r_N[1]*r_u_sol_frac(1,1) + r_N[2]*r_u(2,1) - r_N[2]*r_u_sol_frac(2,1));
const double crRHS27 = crRHS1*crRHS16 + crRHS15*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS28 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS29 = -crRHS22*(-crRHS27 - crRHS28 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1)) + crRHS26 + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS30 = crRHS24*crRHS29;
const double crRHS31 = crRHS19*crRHS20;
const double crRHS32 = crRHS31*crRHS5;
const double crRHS33 = crRHS32*gauss_weight;
const double crRHS34 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS35 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS36 = crRHS22*r_N[0];
const double crRHS37 = sigma*tau_u;
const double crRHS38 = crRHS37*r_N[0];
const double crRHS39 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crRHS40 = crRHS10*crRHS39;
const double crRHS41 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crRHS42 = crRHS12*crRHS41;
const double crRHS43 = crRHS32*tau_c*(-crRHS2 + crRHS40*crRHS6 + crRHS42*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS44 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crRHS45 = crRHS44*tau_u;
const double crRHS46 = crRHS36*crRHS45;
const double crRHS47 = crRHS15*r_DN(0,0) + crRHS16*r_DN(0,1);
const double crRHS48 = crRHS22*tau_u;
const double crRHS49 = crRHS47*crRHS48;
const double crRHS50 = (crRHS11*crRHS16 + crRHS15*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS51 = crRHS50*r_N[0];
const double crRHS52 = crRHS31*crRHS51;
const double crRHS53 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS54 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2];
const double crRHS55 = crRHS39*kappa;
const double crRHS56 = crRHS41*kappa;
const double crRHS57 = crRHS5*dp_th_dt;
const double crRHS58 = crRHS57*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2]);
const double crRHS59 = crRHS20*crRHS7;
const double crRHS60 = crRHS21*(crRHS40 + crRHS42);
const double crRHS61 = tau_t*(-crRHS21*(crRHS15*crRHS39 + crRHS16*crRHS41 + crRHS3) + crRHS54 + crRHS58);
const double crRHS62 = crRHS57*crRHS61;
const double crRHS63 = crRHS21*crRHS61;
const double crRHS64 = crRHS44*crRHS63;
const double crRHS65 = crRHS20*crRHS61*p_th;
const double crRHS66 = crRHS22*r_N[1];
const double crRHS67 = crRHS37*r_N[1];
const double crRHS68 = crRHS45*crRHS66;
const double crRHS69 = crRHS15*r_DN(1,0) + crRHS16*r_DN(1,1);
const double crRHS70 = crRHS48*crRHS69;
const double crRHS71 = crRHS31*crRHS50;
const double crRHS72 = crRHS71*r_N[1];
const double crRHS73 = crRHS50*crRHS65;
const double crRHS74 = crRHS22*r_N[2];
const double crRHS75 = crRHS37*r_N[2];
const double crRHS76 = crRHS45*crRHS74;
const double crRHS77 = crRHS15*r_DN(2,0) + crRHS16*r_DN(2,1);
const double crRHS78 = crRHS48*crRHS77;
const double crRHS79 = crRHS71*r_N[2];
rRHS[0]+=-crRHS33*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS25*r_DN(0,0) + crRHS30*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS14*r_N[0] + crRHS17*crRHS36 + crRHS18*crRHS36 - crRHS23*crRHS38 + crRHS23*crRHS46 + crRHS23*crRHS49 - crRHS25*crRHS52 - crRHS34*r_DN(0,0) - crRHS35*crRHS36 - crRHS43*r_DN(0,0) + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS26*r_N[0] + crRHS27*crRHS36 + crRHS28*crRHS36 - crRHS29*crRHS38 + crRHS29*crRHS46 + crRHS29*crRHS49 - crRHS30*crRHS52 - crRHS34*r_DN(0,1) - crRHS36*crRHS53 - crRHS43*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS47*crRHS63 - crRHS51*crRHS65 + crRHS54*r_N[0] - crRHS55*r_DN(0,0) - crRHS56*r_DN(0,1) + crRHS58*r_N[0] - crRHS59*r_N[0] - crRHS60*r_N[0] + crRHS62*r_N[0] + crRHS64*r_N[0]);
rRHS[4]+=-crRHS33*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS25*r_DN(1,0) + crRHS30*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS14*r_N[1] + crRHS17*crRHS66 + crRHS18*crRHS66 - crRHS23*crRHS67 + crRHS23*crRHS68 + crRHS23*crRHS70 - crRHS25*crRHS72 - crRHS34*r_DN(1,0) - crRHS35*crRHS66 - crRHS43*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS26*r_N[1] + crRHS27*crRHS66 + crRHS28*crRHS66 - crRHS29*crRHS67 + crRHS29*crRHS68 + crRHS29*crRHS70 - crRHS30*crRHS72 - crRHS34*r_DN(1,1) - crRHS43*r_DN(1,1) - crRHS53*crRHS66 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS54*r_N[1] - crRHS55*r_DN(1,0) - crRHS56*r_DN(1,1) + crRHS58*r_N[1] - crRHS59*r_N[1] - crRHS60*r_N[1] + crRHS62*r_N[1] + crRHS63*crRHS69 + crRHS64*r_N[1] - crRHS73*r_N[1]);
rRHS[8]+=-crRHS33*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS25*r_DN(2,0) + crRHS30*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS14*r_N[2] + crRHS17*crRHS74 + crRHS18*crRHS74 - crRHS23*crRHS75 + crRHS23*crRHS76 + crRHS23*crRHS78 - crRHS25*crRHS79 - crRHS34*r_DN(2,0) - crRHS35*crRHS74 - crRHS43*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS26*r_N[2] + crRHS27*crRHS74 + crRHS28*crRHS74 - crRHS29*crRHS75 + crRHS29*crRHS76 + crRHS29*crRHS78 - crRHS30*crRHS79 - crRHS34*r_DN(2,1) - crRHS43*r_DN(2,1) - crRHS53*crRHS74 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS54*r_N[2] - crRHS55*r_DN(2,0) - crRHS56*r_DN(2,1) + crRHS58*r_N[2] - crRHS59*r_N[2] - crRHS60*r_N[2] + crRHS62*r_N[2] + crRHS63*crRHS77 + crRHS64*r_N[2] - crRHS73*r_N[2]);

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const double c_p = rData.SpecificHeat;
    const double sigma = rData.Resistance;
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
    const auto& r_u_sol_frac = rData.SolidFractionVelocity;

    const auto& r_t_lin = rData.Temperature;
    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 4> lin_u_conv = r_u - r_u_mesh;

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

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0) + r_DN(3,0)*r_u(3,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1) + r_DN(3,1)*r_u(3,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]) + r_N[3]*(bdf0*r_t[3] + bdf1*r_t_n[3] + bdf2*r_t_nn[3]);
const double crRHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crRHS5 = 1.0/crRHS4;
const double crRHS6 = crRHS5*p_th;
const double crRHS7 = crRHS3*crRHS6;
const double crRHS8 = -crRHS7 + dp_th_dt;
const double crRHS9 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crRHS10 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0)) + r_N[3]*(r_u(3,0) - r_u_mesh(3,0));
const double crRHS11 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crRHS12 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1)) + r_N[3]*(r_u(3,1) - r_u_mesh(3,1));
const double crRHS13 = crRHS6*(crRHS10*crRHS9 + crRHS11*crRHS12);
const double crRHS14 = sigma*(r_N[0]*r_u(0,0) - r_N[0]*r_u_sol_frac(0,0) + r_N[1]*r_u(1,0) - r_N[1]*r_u_sol_frac(1,0) + r_N[2]*r_u(2,0) - r_N[2]*r_u_sol_frac(2,0) + r_N[3]*r_u(3,0) - r_N[3]*r_u_sol_frac(3,0));
const double crRHS15 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2] + lin_u_conv(3,0)*r_N[3];
const double crRHS16 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2] + lin_u_conv(3,1)*r_N[3];
const double crRHS17 = crRHS0*crRHS15 + crRHS16*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS18 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0)) + r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS19 = 1.0/c_p;
const double crRHS20 = gamma/(gamma - 1.0);
const double crRHS21 = crRHS20*crRHS6;
const double crRHS22 = crRHS19*crRHS21;
const double crRHS23 = crRHS14 - crRHS22*(-crRHS17 - crRHS18 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS24 = p_th*tau_u;
const double crRHS25 = crRHS23*crRHS24;
const double crRHS26 = sigma*(r_N[0]*r_u(0,1) - r_N[0]*r_u_sol_frac(0,1) + r_N[1]*r_u(1,1) - r_N[1]*r_u_sol_frac(1,1) + r_N[2]*r_u(2,1) - r_N[2]*r_u_sol_frac(2,1) + r_N[3]*r_u(3,1) - r_N[3]*r_u_sol_frac(3,1));
const double crRHS27 = crRHS1*crRHS16 + crRHS15*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS28 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1)) + r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS29 = -crRHS22*(-crRHS27 - crRHS28 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1)) + crRHS26 + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS30 = crRHS24*crRHS29;
const double crRHS31 = crRHS19*crRHS20;
const double crRHS32 = crRHS31*crRHS5;
const double crRHS33 = crRHS32*gauss_weight;
const double crRHS34 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS35 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS36 = crRHS22*r_N[0];
const double crRHS37 = sigma*tau_u;
const double crRHS38 = crRHS37*r_N[0];
const double crRHS39 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crRHS40 = crRHS10*crRHS39;
const double crRHS41 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crRHS42 = crRHS12*crRHS41;
const double crRHS43 = crRHS32*tau_c*(-crRHS2 + crRHS40*crRHS6 + crRHS42*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS44 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crRHS45 = crRHS44*tau_u;
const double crRHS46 = crRHS36*crRHS45;
const double crRHS47 = crRHS15*r_DN(0,0) + crRHS16*r_DN(0,1);
const double crRHS48 = crRHS22*tau_u;
const double crRHS49 = crRHS47*crRHS48;
const double crRHS50 = (crRHS11*crRHS16 + crRHS15*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS51 = crRHS50*r_N[0];
const double crRHS52 = crRHS31*crRHS51;
const double crRHS53 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS54 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2] + r_N[3]*r_heat_fl[3];
const double crRHS55 = crRHS39*kappa;
const double crRHS56 = crRHS41*kappa;
const double crRHS57 = crRHS5*dp_th_dt;
const double crRHS58 = crRHS57*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3]);
const double crRHS59 = crRHS20*crRHS7;
const double crRHS60 = crRHS21*(crRHS40 + crRHS42);
const double crRHS61 = tau_t*(-crRHS21*(crRHS15*crRHS39 + crRHS16*crRHS41 + crRHS3) + crRHS54 + crRHS58);
const double crRHS62 = crRHS57*crRHS61;
const double crRHS63 = crRHS21*crRHS61;
const double crRHS64 = crRHS44*crRHS63;
const double crRHS65 = crRHS20*crRHS61*p_th;
const double crRHS66 = crRHS22*r_N[1];
const double crRHS67 = crRHS37*r_N[1];
const double crRHS68 = crRHS45*crRHS66;
const double crRHS69 = crRHS15*r_DN(1,0) + crRHS16*r_DN(1,1);
const double crRHS70 = crRHS48*crRHS69;
const double crRHS71 = crRHS31*crRHS50;
const double crRHS72 = crRHS71*r_N[1];
const double crRHS73 = crRHS50*crRHS65;
const double crRHS74 = crRHS22*r_N[2];
const double crRHS75 = crRHS37*r_N[2];
const double crRHS76 = crRHS45*crRHS74;
const double crRHS77 = crRHS15*r_DN(2,0) + crRHS16*r_DN(2,1);
const double crRHS78 = crRHS48*crRHS77;
const double crRHS79 = crRHS71*r_N[2];
const double crRHS80 = crRHS22*r_N[3];
const double crRHS81 = crRHS37*r_N[3];
const double crRHS82 = crRHS45*crRHS80;
const double crRHS83 = crRHS15*r_DN(3,0) + crRHS16*r_DN(3,1);
const double crRHS84 = crRHS48*crRHS83;
const double crRHS85 = crRHS71*r_N[3];
rRHS[0]+=-crRHS33*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS25*r_DN(0,0) + crRHS30*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS14*r_N[0] + crRHS17*crRHS36 + crRHS18*crRHS36 - crRHS23*crRHS38 + crRHS23*crRHS46 + crRHS23*crRHS49 - crRHS25*crRHS52 - crRHS34*r_DN(0,0) - crRHS35*crRHS36 - crRHS43*r_DN(0,0) + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS26*r_N[0] + crRHS27*crRHS36 + crRHS28*crRHS36 - crRHS29*crRHS38 + crRHS29*crRHS46 + crRHS29*crRHS49 - crRHS30*crRHS52 - crRHS34*r_DN(0,1) - crRHS36*crRHS53 - crRHS43*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS47*crRHS63 - crRHS51*crRHS65 + crRHS54*r_N[0] - crRHS55*r_DN(0,0) - crRHS56*r_DN(0,1) + crRHS58*r_N[0] - crRHS59*r_N[0] - crRHS60*r_N[0] + crRHS62*r_N[0] + crRHS64*r_N[0]);
rRHS[4]+=-crRHS33*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS25*r_DN(1,0) + crRHS30*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS14*r_N[1] + crRHS17*crRHS66 + crRHS18*crRHS66 - crRHS23*crRHS67 + crRHS23*crRHS68 + crRHS23*crRHS70 - crRHS25*crRHS72 - crRHS34*r_DN(1,0) - crRHS35*crRHS66 - crRHS43*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS26*r_N[1] + crRHS27*crRHS66 + crRHS28*crRHS66 - crRHS29*crRHS67 + crRHS29*crRHS68 + crRHS29*crRHS70 - crRHS30*crRHS72 - crRHS34*r_DN(1,1) - crRHS43*r_DN(1,1) - crRHS53*crRHS66 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS54*r_N[1] - crRHS55*r_DN(1,0) - crRHS56*r_DN(1,1) + crRHS58*r_N[1] - crRHS59*r_N[1] - crRHS60*r_N[1] + crRHS62*r_N[1] + crRHS63*crRHS69 + crRHS64*r_N[1] - crRHS73*r_N[1]);
rRHS[8]+=-crRHS33*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS25*r_DN(2,0) + crRHS30*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS14*r_N[2] + crRHS17*crRHS74 + crRHS18*crRHS74 - crRHS23*crRHS75 + crRHS23*crRHS76 + crRHS23*crRHS78 - crRHS25*crRHS79 - crRHS34*r_DN(2,0) - crRHS35*crRHS74 - crRHS43*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS26*r_N[2] + crRHS27*crRHS74 + crRHS28*crRHS74 - crRHS29*crRHS75 + crRHS29*crRHS76 + crRHS29*crRHS78 - crRHS30*crRHS79 - crRHS34*r_DN(2,1) - crRHS43*r_DN(2,1) - crRHS53*crRHS74 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS54*r_N[2] - crRHS55*r_DN(2,0) - crRHS56*r_DN(2,1) + crRHS58*r_N[2] - crRHS59*r_N[2] - crRHS60*r_N[2] + crRHS62*r_N[2] + crRHS63*crRHS77 + crRHS64*r_N[2] - crRHS73*r_N[2]);
rRHS[12]+=-crRHS33*(-crRHS13*r_N[3] + crRHS2*r_N[3] + crRHS25*r_DN(3,0) + crRHS30*r_DN(3,1) + crRHS8*r_N[3]);
rRHS[13]+=-gauss_weight*(crRHS14*r_N[3] + crRHS17*crRHS80 + crRHS18*crRHS80 - crRHS23*crRHS81 + crRHS23*crRHS82 + crRHS23*crRHS84 - crRHS25*crRHS85 - crRHS34*r_DN(3,0) - crRHS35*crRHS80 - crRHS43*r_DN(3,0) + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-gauss_weight*(crRHS26*r_N[3] + crRHS27*crRHS80 + crRHS28*crRHS80 - crRHS29*crRHS81 + crRHS29*crRHS82 + crRHS29*crRHS84 - crRHS30*crRHS85 - crRHS34*r_DN(3,1) - crRHS43*r_DN(3,1) - crRHS53*crRHS80 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=gauss_weight*(crRHS54*r_N[3] - crRHS55*r_DN(3,0) - crRHS56*r_DN(3,1) + crRHS58*r_N[3] - crRHS59*r_N[3] - crRHS60*r_N[3] + crRHS62*r_N[3] + crRHS63*crRHS83 + crRHS64*r_N[3] - crRHS73*r_N[3]);

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
    const double sigma = rData.Resistance;
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

    rTauPressure = mu / rho_gauss + stab_c2 * norm_u_conv * h / stab_c1 + stab_c3 * sigma * h / stab_c1 / rho_gauss; // Pressure subscale stabilization operator
    rTauVelocity = 1.0 / (stab_c1 * mu / std::pow(h, 2) + stab_c2 * rho_gauss * norm_u_conv / h + stab_c3 * sigma / h); // Velocity subscale stabilization operator
    rTauTemperature = 1.0 / (stab_c1 * kappa / std::pow(h, 2) + stab_c2 * rho_gauss * c_p * norm_u_conv / h); // Temperature subscale stabilization operator
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