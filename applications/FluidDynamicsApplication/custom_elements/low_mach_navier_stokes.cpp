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
const double crLHS2 = gamma - 1.0;
const double crLHS3 = 1.0/crLHS2;
const double crLHS4 = crLHS3*gamma*p_th;
const double crLHS5 = crLHS4*gauss_weight;
const double crLHS6 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crLHS7 = 1.0/crLHS6;
const double crLHS8 = 1.0/c_p;
const double crLHS9 = crLHS7*crLHS8;
const double crLHS10 = crLHS5*crLHS9;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS13 = r_N[0]*r_N[0];
const double crLHS14 = crLHS13*crLHS7;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS19*crLHS4;
const double crLHS21 = crLHS9*tau_u;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS24 = bdf0*crLHS4;
const double crLHS25 = 1.0/(crLHS6*crLHS6);
const double crLHS26 = crLHS25*crLHS8;
const double crLHS27 = crLHS24*crLHS26*gauss_weight;
const double crLHS28 = r_DN(0,0)*r_DN(1,0);
const double crLHS29 = r_DN(0,1)*r_DN(1,1);
const double crLHS30 = crLHS11*(crLHS28 + crLHS29);
const double crLHS31 = r_DN(1,0)*r_N[0];
const double crLHS32 = crLHS7*r_N[0];
const double crLHS33 = crLHS32*r_N[1];
const double crLHS34 = -crLHS12*crLHS33;
const double crLHS35 = bdf0*r_N[1];
const double crLHS36 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS37 = crLHS35 + crLHS36;
const double crLHS38 = crLHS37*crLHS4;
const double crLHS39 = crLHS21*crLHS38;
const double crLHS40 = r_DN(1,1)*r_N[0];
const double crLHS41 = -crLHS23*crLHS33;
const double crLHS42 = crLHS5*r_N[1];
const double crLHS43 = crLHS15*crLHS26;
const double crLHS44 = -crLHS42*crLHS43;
const double crLHS45 = r_DN(0,0)*r_DN(2,0);
const double crLHS46 = r_DN(0,1)*r_DN(2,1);
const double crLHS47 = crLHS11*(crLHS45 + crLHS46);
const double crLHS48 = r_DN(2,0)*r_N[0];
const double crLHS49 = crLHS32*r_N[2];
const double crLHS50 = -crLHS12*crLHS49;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS51 + crLHS52;
const double crLHS54 = crLHS4*crLHS53;
const double crLHS55 = crLHS21*crLHS54;
const double crLHS56 = r_DN(2,1)*r_N[0];
const double crLHS57 = -crLHS23*crLHS49;
const double crLHS58 = crLHS5*r_N[2];
const double crLHS59 = -crLHS43*crLHS58;
const double crLHS60 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crLHS61 = crLHS4*r_N[0];
const double crLHS62 = crLHS12*crLHS16 + crLHS17*crLHS23;
const double crLHS63 = crLHS62*tau_u;
const double crLHS64 = crLHS26*crLHS63;
const double crLHS65 = gauss_weight*(crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_N[0]*tau_u - crLHS61*crLHS64 - r_N[0]);
const double crLHS66 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS67 = r_C(0,2)*r_DN(0,0);
const double crLHS68 = crLHS67 + r_C(2,2)*r_DN(0,1);
const double crLHS69 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crLHS70 = crLHS32*crLHS69;
const double crLHS71 = -crLHS70 + r_DN(0,0);
const double crLHS72 = r_DN(0,0)*tau_c;
const double crLHS73 = crLHS4*crLHS9;
const double crLHS74 = crLHS72*crLHS73;
const double crLHS75 = crLHS14*crLHS24;
const double crLHS76 = crLHS32*crLHS4;
const double crLHS77 = crLHS76*crLHS8;
const double crLHS78 = 1.0/(crLHS2*crLHS2)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS79 = crLHS19*crLHS78;
const double crLHS80 = crLHS25*tau_u;
const double crLHS81 = crLHS18*crLHS80;
const double crLHS82 = crLHS79*r_N[0];
const double crLHS83 = crLHS25*crLHS60*tau_u;
const double crLHS84 = 1.0/(crLHS6*crLHS6*crLHS6);
const double crLHS85 = crLHS63*crLHS84;
const double crLHS86 = crLHS18*crLHS77 + crLHS75*crLHS8 + crLHS79*crLHS81 + crLHS82*crLHS83 - crLHS82*crLHS85;
const double crLHS87 = crLHS67 + r_C(0,1)*r_DN(0,1);
const double crLHS88 = r_C(1,2)*r_DN(0,1);
const double crLHS89 = crLHS88 + r_C(2,2)*r_DN(0,0);
const double crLHS90 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crLHS91 = crLHS32*crLHS90;
const double crLHS92 = -crLHS91 + r_DN(0,1);
const double crLHS93 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0));
const double crLHS94 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1));
const double crLHS95 = crLHS93*r_DN(0,0) + crLHS94*r_DN(0,1);
const double crLHS96 = crLHS15 + crLHS95;
const double crLHS97 = crLHS26*crLHS5;
const double crLHS98 = crLHS72*crLHS97;
const double crLHS99 = r_DN(0,0)*r_N[1];
const double crLHS100 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS101 = r_C(0,2)*r_DN(1,0);
const double crLHS102 = crLHS101 + r_C(2,2)*r_DN(1,1);
const double crLHS103 = crLHS7*r_N[1];
const double crLHS104 = crLHS103*crLHS69;
const double crLHS105 = -crLHS104 + r_DN(1,0);
const double crLHS106 = crLHS103*crLHS4;
const double crLHS107 = crLHS106*crLHS15;
const double crLHS108 = crLHS107*crLHS8;
const double crLHS109 = crLHS78*crLHS81;
const double crLHS110 = crLHS78*r_N[0];
const double crLHS111 = crLHS110*crLHS37;
const double crLHS112 = crLHS108 + crLHS109*crLHS37 + crLHS111*crLHS83 - crLHS111*crLHS85 + crLHS36*crLHS77;
const double crLHS113 = crLHS101 + r_C(0,1)*r_DN(1,1);
const double crLHS114 = r_C(1,2)*r_DN(1,1);
const double crLHS115 = crLHS114 + r_C(2,2)*r_DN(1,0);
const double crLHS116 = crLHS103*crLHS90;
const double crLHS117 = -crLHS116 + r_DN(1,1);
const double crLHS118 = crLHS93*r_DN(1,0) + crLHS94*r_DN(1,1);
const double crLHS119 = crLHS118 + crLHS35;
const double crLHS120 = r_DN(0,0)*r_N[2];
const double crLHS121 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS122 = r_C(0,2)*r_DN(2,0);
const double crLHS123 = crLHS122 + r_C(2,2)*r_DN(2,1);
const double crLHS124 = crLHS7*r_N[2];
const double crLHS125 = -crLHS124*crLHS69 + r_DN(2,0);
const double crLHS126 = crLHS124*crLHS4;
const double crLHS127 = crLHS126*crLHS15;
const double crLHS128 = crLHS127*crLHS8;
const double crLHS129 = crLHS110*crLHS53;
const double crLHS130 = crLHS109*crLHS53 + crLHS128 + crLHS129*crLHS83 - crLHS129*crLHS85 + crLHS52*crLHS77;
const double crLHS131 = crLHS122 + r_C(0,1)*r_DN(2,1);
const double crLHS132 = r_C(1,2)*r_DN(2,1);
const double crLHS133 = crLHS132 + r_C(2,2)*r_DN(2,0);
const double crLHS134 = -crLHS124*crLHS90 + r_DN(2,1);
const double crLHS135 = crLHS93*r_DN(2,0) + crLHS94*r_DN(2,1);
const double crLHS136 = crLHS135 + crLHS51;
const double crLHS137 = crLHS88 + r_C(0,1)*r_DN(0,0);
const double crLHS138 = crLHS73*tau_c;
const double crLHS139 = crLHS138*r_DN(0,1);
const double crLHS140 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS141 = crLHS97*tau_c;
const double crLHS142 = crLHS141*r_DN(0,1);
const double crLHS143 = r_DN(0,1)*r_N[1];
const double crLHS144 = crLHS114 + r_C(0,1)*r_DN(1,0);
const double crLHS145 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS146 = r_DN(0,1)*r_N[2];
const double crLHS147 = crLHS132 + r_C(0,1)*r_DN(2,0);
const double crLHS148 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS149 = crLHS14*crLHS5;
const double crLHS150 = -crLHS20 + dp_th_dt*r_N[0];
const double crLHS151 = crLHS62*crLHS84*tau_t;
const double crLHS152 = crLHS42*crLHS70;
const double crLHS153 = crLHS42*crLHS91;
const double crLHS154 = -crLHS38 + dp_th_dt*r_N[1];
const double crLHS155 = crLHS107 + crLHS28*kappa + crLHS29*kappa - crLHS7*dp_th_dt*r_N[0]*r_N[1];
const double crLHS156 = crLHS58*crLHS70;
const double crLHS157 = crLHS58*crLHS91;
const double crLHS158 = -crLHS54 + dp_th_dt*r_N[2];
const double crLHS159 = crLHS127 + crLHS45*kappa + crLHS46*kappa - crLHS7*dp_th_dt*r_N[0]*r_N[2];
const double crLHS160 = r_DN(1,0)*r_DN(1,0);
const double crLHS161 = r_DN(1,1)*r_DN(1,1);
const double crLHS162 = r_N[1]*r_N[1];
const double crLHS163 = crLHS162*crLHS7;
const double crLHS164 = r_DN(1,0)*r_DN(2,0);
const double crLHS165 = r_DN(1,1)*r_DN(2,1);
const double crLHS166 = crLHS11*(crLHS164 + crLHS165);
const double crLHS167 = r_DN(2,0)*r_N[1];
const double crLHS168 = crLHS103*r_N[2];
const double crLHS169 = -crLHS12*crLHS168;
const double crLHS170 = r_DN(2,1)*r_N[1];
const double crLHS171 = -crLHS168*crLHS23;
const double crLHS172 = -crLHS26*crLHS35*crLHS58;
const double crLHS173 = crLHS138*r_DN(1,0);
const double crLHS174 = crLHS106*crLHS8;
const double crLHS175 = crLHS36*crLHS80;
const double crLHS176 = crLHS79*r_N[1];
const double crLHS177 = crLHS108 + crLHS174*crLHS18 + crLHS175*crLHS79 + crLHS176*crLHS83 - crLHS176*crLHS85;
const double crLHS178 = crLHS141*r_DN(1,0);
const double crLHS179 = crLHS4*r_N[1];
const double crLHS180 = gauss_weight*(-crLHS179*crLHS64 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS181 = crLHS163*crLHS24;
const double crLHS182 = crLHS175*crLHS78;
const double crLHS183 = crLHS78*r_N[1];
const double crLHS184 = crLHS183*crLHS37;
const double crLHS185 = crLHS174*crLHS36 + crLHS181*crLHS8 + crLHS182*crLHS37 + crLHS184*crLHS83 - crLHS184*crLHS85;
const double crLHS186 = r_DN(1,0)*r_N[2];
const double crLHS187 = crLHS126*crLHS35;
const double crLHS188 = crLHS187*crLHS8;
const double crLHS189 = crLHS183*crLHS53;
const double crLHS190 = crLHS174*crLHS52 + crLHS182*crLHS53 + crLHS188 + crLHS189*crLHS83 - crLHS189*crLHS85;
const double crLHS191 = crLHS138*r_DN(1,1);
const double crLHS192 = crLHS141*r_DN(1,1);
const double crLHS193 = r_DN(1,1)*r_N[2];
const double crLHS194 = crLHS163*crLHS5;
const double crLHS195 = crLHS104*crLHS58;
const double crLHS196 = crLHS116*crLHS58;
const double crLHS197 = crLHS164*kappa + crLHS165*kappa + crLHS187 - crLHS7*dp_th_dt*r_N[1]*r_N[2];
const double crLHS198 = r_DN(2,0)*r_DN(2,0);
const double crLHS199 = r_DN(2,1)*r_DN(2,1);
const double crLHS200 = r_N[2]*r_N[2];
const double crLHS201 = crLHS200*crLHS7;
const double crLHS202 = crLHS138*r_DN(2,0);
const double crLHS203 = crLHS126*crLHS8;
const double crLHS204 = crLHS52*crLHS80;
const double crLHS205 = crLHS79*r_N[2];
const double crLHS206 = crLHS128 + crLHS18*crLHS203 + crLHS204*crLHS79 + crLHS205*crLHS83 - crLHS205*crLHS85;
const double crLHS207 = crLHS141*r_DN(2,0);
const double crLHS208 = crLHS204*crLHS78;
const double crLHS209 = crLHS78*r_N[2];
const double crLHS210 = crLHS209*crLHS37;
const double crLHS211 = crLHS188 + crLHS203*crLHS36 + crLHS208*crLHS37 + crLHS210*crLHS83 - crLHS210*crLHS85;
const double crLHS212 = crLHS4*r_N[2];
const double crLHS213 = gauss_weight*(-crLHS212*crLHS64 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS214 = crLHS201*crLHS24;
const double crLHS215 = crLHS209*crLHS53;
const double crLHS216 = crLHS203*crLHS52 + crLHS208*crLHS53 + crLHS214*crLHS8 + crLHS215*crLHS83 - crLHS215*crLHS85;
const double crLHS217 = crLHS138*r_DN(2,1);
const double crLHS218 = crLHS141*r_DN(2,1);
const double crLHS219 = crLHS201*crLHS5;
rLHS(0,0)+=crLHS11*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS10*(-crLHS12*crLHS14 + crLHS22*r_DN(0,0) + r_DN(0,0)*r_N[0]);
rLHS(0,2)+=crLHS10*(-crLHS14*crLHS23 + crLHS22*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS13*crLHS27;
rLHS(0,4)+=crLHS30;
rLHS(0,5)+=crLHS10*(crLHS31 + crLHS34 + crLHS39*r_DN(0,0));
rLHS(0,6)+=crLHS10*(crLHS39*r_DN(0,1) + crLHS40 + crLHS41);
rLHS(0,7)+=crLHS44;
rLHS(0,8)+=crLHS47;
rLHS(0,9)+=crLHS10*(crLHS48 + crLHS50 + crLHS55*r_DN(0,0));
rLHS(0,10)+=crLHS10*(crLHS55*r_DN(0,1) + crLHS56 + crLHS57);
rLHS(0,11)+=crLHS59;
rLHS(1,0)+=crLHS65*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS66*r_DN(0,0) + crLHS68*r_DN(0,1) + crLHS71*crLHS74 + crLHS86);
rLHS(1,2)+=gauss_weight*(crLHS74*crLHS92 + crLHS87*r_DN(0,0) + crLHS89*r_DN(0,1));
rLHS(1,3)+=-crLHS96*crLHS98;
rLHS(1,4)+=gauss_weight*(crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS31*crLHS4*crLHS64 - crLHS99);
rLHS(1,5)+=gauss_weight*(crLHS100*r_DN(0,0) + crLHS102*r_DN(0,1) + crLHS105*crLHS74 + crLHS112);
rLHS(1,6)+=gauss_weight*(crLHS113*r_DN(0,0) + crLHS115*r_DN(0,1) + crLHS117*crLHS74);
rLHS(1,7)+=-crLHS119*crLHS98;
rLHS(1,8)+=gauss_weight*(-crLHS120 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS4*crLHS48*crLHS64);
rLHS(1,9)+=gauss_weight*(crLHS121*r_DN(0,0) + crLHS123*r_DN(0,1) + crLHS125*crLHS74 + crLHS130);
rLHS(1,10)+=gauss_weight*(crLHS131*r_DN(0,0) + crLHS133*r_DN(0,1) + crLHS134*crLHS74);
rLHS(1,11)+=-crLHS136*crLHS98;
rLHS(2,0)+=crLHS65*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS137*r_DN(0,1) + crLHS139*crLHS71 + crLHS68*r_DN(0,0));
rLHS(2,2)+=gauss_weight*(crLHS139*crLHS92 + crLHS140*r_DN(0,1) + crLHS86 + crLHS89*r_DN(0,0));
rLHS(2,3)+=-crLHS142*crLHS96;
rLHS(2,4)+=gauss_weight*(-crLHS143 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS4*crLHS40*crLHS64);
rLHS(2,5)+=gauss_weight*(crLHS102*r_DN(0,0) + crLHS105*crLHS139 + crLHS144*r_DN(0,1));
rLHS(2,6)+=gauss_weight*(crLHS112 + crLHS115*r_DN(0,0) + crLHS117*crLHS139 + crLHS145*r_DN(0,1));
rLHS(2,7)+=-crLHS119*crLHS142;
rLHS(2,8)+=gauss_weight*(-crLHS146 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS4*crLHS56*crLHS64);
rLHS(2,9)+=gauss_weight*(crLHS123*r_DN(0,0) + crLHS125*crLHS139 + crLHS147*r_DN(0,1));
rLHS(2,10)+=gauss_weight*(crLHS130 + crLHS133*r_DN(0,0) + crLHS134*crLHS139 + crLHS148*r_DN(0,1));
rLHS(2,11)+=-crLHS136*crLHS142;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS149*crLHS69;
rLHS(3,2)+=crLHS149*crLHS90;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa + crLHS13*crLHS7*dp_th_dt - crLHS150*crLHS151*crLHS61 + crLHS150*crLHS18*crLHS25*crLHS3*gamma*p_th*tau_t + crLHS150*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS150*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS75 - crLHS76*crLHS95);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS152;
rLHS(3,6)+=crLHS153;
rLHS(3,7)+=-gauss_weight*(-crLHS118*crLHS76 - crLHS151*crLHS154*crLHS61 + crLHS154*crLHS18*crLHS25*crLHS3*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS155);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS156;
rLHS(3,10)+=crLHS157;
rLHS(3,11)+=-gauss_weight*(-crLHS135*crLHS76 - crLHS151*crLHS158*crLHS61 + crLHS158*crLHS18*crLHS25*crLHS3*gamma*p_th*tau_t + crLHS158*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS158*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS159);
rLHS(4,0)+=crLHS30;
rLHS(4,1)+=crLHS10*(crLHS22*r_DN(1,0) + crLHS34 + crLHS99);
rLHS(4,2)+=crLHS10*(crLHS143 + crLHS22*r_DN(1,1) + crLHS41);
rLHS(4,3)+=crLHS44;
rLHS(4,4)+=crLHS11*(crLHS160 + crLHS161);
rLHS(4,5)+=crLHS10*(-crLHS12*crLHS163 + crLHS39*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS163*crLHS23 + crLHS39*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS162*crLHS27;
rLHS(4,8)+=crLHS166;
rLHS(4,9)+=crLHS10*(crLHS167 + crLHS169 + crLHS55*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS170 + crLHS171 + crLHS55*r_DN(1,1));
rLHS(4,11)+=crLHS172;
rLHS(5,0)+=gauss_weight*(crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS31 - crLHS4*crLHS64*crLHS99);
rLHS(5,1)+=gauss_weight*(crLHS173*crLHS71 + crLHS177 + crLHS66*r_DN(1,0) + crLHS68*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS173*crLHS92 + crLHS87*r_DN(1,0) + crLHS89*r_DN(1,1));
rLHS(5,3)+=-crLHS178*crLHS96;
rLHS(5,4)+=crLHS180*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS100*r_DN(1,0) + crLHS102*r_DN(1,1) + crLHS105*crLHS173 + crLHS185);
rLHS(5,6)+=gauss_weight*(crLHS113*r_DN(1,0) + crLHS115*r_DN(1,1) + crLHS117*crLHS173);
rLHS(5,7)+=-crLHS119*crLHS178;
rLHS(5,8)+=gauss_weight*(-crLHS167*crLHS4*crLHS64 - crLHS186 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS121*r_DN(1,0) + crLHS123*r_DN(1,1) + crLHS125*crLHS173 + crLHS190);
rLHS(5,10)+=gauss_weight*(crLHS131*r_DN(1,0) + crLHS133*r_DN(1,1) + crLHS134*crLHS173);
rLHS(5,11)+=-crLHS136*crLHS178;
rLHS(6,0)+=gauss_weight*(-crLHS143*crLHS4*crLHS64 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS40);
rLHS(6,1)+=gauss_weight*(crLHS137*r_DN(1,1) + crLHS191*crLHS71 + crLHS68*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS140*r_DN(1,1) + crLHS177 + crLHS191*crLHS92 + crLHS89*r_DN(1,0));
rLHS(6,3)+=-crLHS192*crLHS96;
rLHS(6,4)+=crLHS180*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS102*r_DN(1,0) + crLHS105*crLHS191 + crLHS144*r_DN(1,1));
rLHS(6,6)+=gauss_weight*(crLHS115*r_DN(1,0) + crLHS117*crLHS191 + crLHS145*r_DN(1,1) + crLHS185);
rLHS(6,7)+=-crLHS119*crLHS192;
rLHS(6,8)+=gauss_weight*(-crLHS170*crLHS4*crLHS64 - crLHS193 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS123*r_DN(1,0) + crLHS125*crLHS191 + crLHS147*r_DN(1,1));
rLHS(6,10)+=gauss_weight*(crLHS133*r_DN(1,0) + crLHS134*crLHS191 + crLHS148*r_DN(1,1) + crLHS190);
rLHS(6,11)+=-crLHS136*crLHS192;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS152;
rLHS(7,2)+=crLHS153;
rLHS(7,3)+=-gauss_weight*(-crLHS106*crLHS95 - crLHS150*crLHS151*crLHS179 + crLHS150*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS150*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS150*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS155);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS194*crLHS69;
rLHS(7,6)+=crLHS194*crLHS90;
rLHS(7,7)+=-gauss_weight*(-crLHS106*crLHS118 - crLHS151*crLHS154*crLHS179 + crLHS154*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS160*kappa - crLHS161*kappa + crLHS162*crLHS7*dp_th_dt - crLHS181);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS195;
rLHS(7,10)+=crLHS196;
rLHS(7,11)+=-gauss_weight*(-crLHS106*crLHS135 - crLHS151*crLHS158*crLHS179 + crLHS158*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS158*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS158*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS197);
rLHS(8,0)+=crLHS47;
rLHS(8,1)+=crLHS10*(crLHS120 + crLHS22*r_DN(2,0) + crLHS50);
rLHS(8,2)+=crLHS10*(crLHS146 + crLHS22*r_DN(2,1) + crLHS57);
rLHS(8,3)+=crLHS59;
rLHS(8,4)+=crLHS166;
rLHS(8,5)+=crLHS10*(crLHS169 + crLHS186 + crLHS39*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS171 + crLHS193 + crLHS39*r_DN(2,1));
rLHS(8,7)+=crLHS172;
rLHS(8,8)+=crLHS11*(crLHS198 + crLHS199);
rLHS(8,9)+=crLHS10*(-crLHS12*crLHS201 + crLHS55*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS201*crLHS23 + crLHS55*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS200*crLHS27;
rLHS(9,0)+=gauss_weight*(-crLHS120*crLHS4*crLHS64 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS48);
rLHS(9,1)+=gauss_weight*(crLHS202*crLHS71 + crLHS206 + crLHS66*r_DN(2,0) + crLHS68*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS202*crLHS92 + crLHS87*r_DN(2,0) + crLHS89*r_DN(2,1));
rLHS(9,3)+=-crLHS207*crLHS96;
rLHS(9,4)+=gauss_weight*(-crLHS167 - crLHS186*crLHS4*crLHS64 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS100*r_DN(2,0) + crLHS102*r_DN(2,1) + crLHS105*crLHS202 + crLHS211);
rLHS(9,6)+=gauss_weight*(crLHS113*r_DN(2,0) + crLHS115*r_DN(2,1) + crLHS117*crLHS202);
rLHS(9,7)+=-crLHS119*crLHS207;
rLHS(9,8)+=crLHS213*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS121*r_DN(2,0) + crLHS123*r_DN(2,1) + crLHS125*crLHS202 + crLHS216);
rLHS(9,10)+=gauss_weight*(crLHS131*r_DN(2,0) + crLHS133*r_DN(2,1) + crLHS134*crLHS202);
rLHS(9,11)+=-crLHS136*crLHS207;
rLHS(10,0)+=gauss_weight*(-crLHS146*crLHS4*crLHS64 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS56);
rLHS(10,1)+=gauss_weight*(crLHS137*r_DN(2,1) + crLHS217*crLHS71 + crLHS68*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS140*r_DN(2,1) + crLHS206 + crLHS217*crLHS92 + crLHS89*r_DN(2,0));
rLHS(10,3)+=-crLHS218*crLHS96;
rLHS(10,4)+=gauss_weight*(-crLHS170 - crLHS193*crLHS4*crLHS64 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS60*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS102*r_DN(2,0) + crLHS105*crLHS217 + crLHS144*r_DN(2,1));
rLHS(10,6)+=gauss_weight*(crLHS115*r_DN(2,0) + crLHS117*crLHS217 + crLHS145*r_DN(2,1) + crLHS211);
rLHS(10,7)+=-crLHS119*crLHS218;
rLHS(10,8)+=crLHS213*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS123*r_DN(2,0) + crLHS125*crLHS217 + crLHS147*r_DN(2,1));
rLHS(10,10)+=gauss_weight*(crLHS133*r_DN(2,0) + crLHS134*crLHS217 + crLHS148*r_DN(2,1) + crLHS216);
rLHS(10,11)+=-crLHS136*crLHS218;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS156;
rLHS(11,2)+=crLHS157;
rLHS(11,3)+=-gauss_weight*(-crLHS126*crLHS95 - crLHS150*crLHS151*crLHS212 + crLHS150*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS150*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS150*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS159);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS195;
rLHS(11,6)+=crLHS196;
rLHS(11,7)+=-gauss_weight*(-crLHS118*crLHS126 - crLHS151*crLHS154*crLHS212 + crLHS154*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS197);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS219*crLHS69;
rLHS(11,10)+=crLHS219*crLHS90;
rLHS(11,11)+=-gauss_weight*(-crLHS126*crLHS135 - crLHS151*crLHS158*crLHS212 + crLHS158*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS158*crLHS25*crLHS3*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS158*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS198*kappa - crLHS199*kappa + crLHS200*crLHS7*dp_th_dt - crLHS214);

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
const double crLHS2 = gamma - 1.0;
const double crLHS3 = 1.0/crLHS2;
const double crLHS4 = crLHS3*gamma*p_th;
const double crLHS5 = crLHS4*gauss_weight;
const double crLHS6 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crLHS7 = 1.0/crLHS6;
const double crLHS8 = 1.0/c_p;
const double crLHS9 = crLHS7*crLHS8;
const double crLHS10 = crLHS5*crLHS9;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS13 = r_N[0]*r_N[0];
const double crLHS14 = crLHS13*crLHS7;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2] + lin_u_conv(3,0)*r_N[3];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2] + lin_u_conv(3,1)*r_N[3];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS19*crLHS4;
const double crLHS21 = crLHS9*tau_u;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS24 = bdf0*crLHS4;
const double crLHS25 = 1.0/(crLHS6*crLHS6);
const double crLHS26 = crLHS25*crLHS8;
const double crLHS27 = crLHS24*crLHS26*gauss_weight;
const double crLHS28 = r_DN(0,0)*r_DN(1,0);
const double crLHS29 = r_DN(0,1)*r_DN(1,1);
const double crLHS30 = crLHS11*(crLHS28 + crLHS29);
const double crLHS31 = r_DN(1,0)*r_N[0];
const double crLHS32 = crLHS7*r_N[0];
const double crLHS33 = crLHS32*r_N[1];
const double crLHS34 = -crLHS12*crLHS33;
const double crLHS35 = bdf0*r_N[1];
const double crLHS36 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS37 = crLHS35 + crLHS36;
const double crLHS38 = crLHS37*crLHS4;
const double crLHS39 = crLHS21*crLHS38;
const double crLHS40 = r_DN(1,1)*r_N[0];
const double crLHS41 = -crLHS23*crLHS33;
const double crLHS42 = crLHS5*r_N[1];
const double crLHS43 = crLHS15*crLHS26;
const double crLHS44 = -crLHS42*crLHS43;
const double crLHS45 = r_DN(0,0)*r_DN(2,0);
const double crLHS46 = r_DN(0,1)*r_DN(2,1);
const double crLHS47 = crLHS11*(crLHS45 + crLHS46);
const double crLHS48 = r_DN(2,0)*r_N[0];
const double crLHS49 = crLHS32*r_N[2];
const double crLHS50 = -crLHS12*crLHS49;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS51 + crLHS52;
const double crLHS54 = crLHS4*crLHS53;
const double crLHS55 = crLHS21*crLHS54;
const double crLHS56 = r_DN(2,1)*r_N[0];
const double crLHS57 = -crLHS23*crLHS49;
const double crLHS58 = crLHS5*r_N[2];
const double crLHS59 = -crLHS43*crLHS58;
const double crLHS60 = r_DN(0,0)*r_DN(3,0);
const double crLHS61 = r_DN(0,1)*r_DN(3,1);
const double crLHS62 = crLHS11*(crLHS60 + crLHS61);
const double crLHS63 = r_DN(3,0)*r_N[0];
const double crLHS64 = crLHS32*r_N[3];
const double crLHS65 = -crLHS12*crLHS64;
const double crLHS66 = bdf0*r_N[3];
const double crLHS67 = crLHS16*r_DN(3,0) + crLHS17*r_DN(3,1);
const double crLHS68 = crLHS66 + crLHS67;
const double crLHS69 = crLHS4*crLHS68;
const double crLHS70 = crLHS21*crLHS69;
const double crLHS71 = r_DN(3,1)*r_N[0];
const double crLHS72 = -crLHS23*crLHS64;
const double crLHS73 = crLHS5*r_N[3];
const double crLHS74 = -crLHS43*crLHS73;
const double crLHS75 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crLHS76 = crLHS4*r_N[0];
const double crLHS77 = crLHS12*crLHS16 + crLHS17*crLHS23;
const double crLHS78 = crLHS77*tau_u;
const double crLHS79 = crLHS26*crLHS78;
const double crLHS80 = gauss_weight*(crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_N[0]*tau_u - crLHS76*crLHS79 - r_N[0]);
const double crLHS81 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS82 = r_C(0,2)*r_DN(0,0);
const double crLHS83 = crLHS82 + r_C(2,2)*r_DN(0,1);
const double crLHS84 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crLHS85 = crLHS32*crLHS84;
const double crLHS86 = -crLHS85 + r_DN(0,0);
const double crLHS87 = r_DN(0,0)*tau_c;
const double crLHS88 = crLHS4*crLHS9;
const double crLHS89 = crLHS87*crLHS88;
const double crLHS90 = crLHS14*crLHS24;
const double crLHS91 = crLHS32*crLHS4;
const double crLHS92 = crLHS8*crLHS91;
const double crLHS93 = 1.0/(crLHS2*crLHS2)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS94 = crLHS19*crLHS93;
const double crLHS95 = crLHS25*tau_u;
const double crLHS96 = crLHS18*crLHS95;
const double crLHS97 = crLHS94*r_N[0];
const double crLHS98 = crLHS25*crLHS75*tau_u;
const double crLHS99 = 1.0/(crLHS6*crLHS6*crLHS6);
const double crLHS100 = crLHS78*crLHS99;
const double crLHS101 = -crLHS100*crLHS97 + crLHS18*crLHS92 + crLHS8*crLHS90 + crLHS94*crLHS96 + crLHS97*crLHS98;
const double crLHS102 = crLHS82 + r_C(0,1)*r_DN(0,1);
const double crLHS103 = r_C(1,2)*r_DN(0,1);
const double crLHS104 = crLHS103 + r_C(2,2)*r_DN(0,0);
const double crLHS105 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crLHS106 = crLHS105*crLHS32;
const double crLHS107 = -crLHS106 + r_DN(0,1);
const double crLHS108 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0)) + r_N[3]*(r_u(3,0) - r_u_mesh(3,0));
const double crLHS109 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1)) + r_N[3]*(r_u(3,1) - r_u_mesh(3,1));
const double crLHS110 = crLHS108*r_DN(0,0) + crLHS109*r_DN(0,1);
const double crLHS111 = crLHS110 + crLHS15;
const double crLHS112 = crLHS26*crLHS5;
const double crLHS113 = crLHS112*crLHS87;
const double crLHS114 = r_DN(0,0)*r_N[1];
const double crLHS115 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS116 = r_C(0,2)*r_DN(1,0);
const double crLHS117 = crLHS116 + r_C(2,2)*r_DN(1,1);
const double crLHS118 = crLHS7*r_N[1];
const double crLHS119 = crLHS118*crLHS84;
const double crLHS120 = -crLHS119 + r_DN(1,0);
const double crLHS121 = crLHS118*crLHS4;
const double crLHS122 = crLHS121*crLHS15;
const double crLHS123 = crLHS122*crLHS8;
const double crLHS124 = crLHS93*crLHS96;
const double crLHS125 = crLHS93*r_N[0];
const double crLHS126 = crLHS125*crLHS37;
const double crLHS127 = -crLHS100*crLHS126 + crLHS123 + crLHS124*crLHS37 + crLHS126*crLHS98 + crLHS36*crLHS92;
const double crLHS128 = crLHS116 + r_C(0,1)*r_DN(1,1);
const double crLHS129 = r_C(1,2)*r_DN(1,1);
const double crLHS130 = crLHS129 + r_C(2,2)*r_DN(1,0);
const double crLHS131 = crLHS105*crLHS118;
const double crLHS132 = -crLHS131 + r_DN(1,1);
const double crLHS133 = crLHS108*r_DN(1,0) + crLHS109*r_DN(1,1);
const double crLHS134 = crLHS133 + crLHS35;
const double crLHS135 = r_DN(0,0)*r_N[2];
const double crLHS136 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS137 = r_C(0,2)*r_DN(2,0);
const double crLHS138 = crLHS137 + r_C(2,2)*r_DN(2,1);
const double crLHS139 = crLHS7*r_N[2];
const double crLHS140 = crLHS139*crLHS84;
const double crLHS141 = -crLHS140 + r_DN(2,0);
const double crLHS142 = crLHS139*crLHS4;
const double crLHS143 = crLHS142*crLHS15;
const double crLHS144 = crLHS143*crLHS8;
const double crLHS145 = crLHS125*crLHS53;
const double crLHS146 = -crLHS100*crLHS145 + crLHS124*crLHS53 + crLHS144 + crLHS145*crLHS98 + crLHS52*crLHS92;
const double crLHS147 = crLHS137 + r_C(0,1)*r_DN(2,1);
const double crLHS148 = r_C(1,2)*r_DN(2,1);
const double crLHS149 = crLHS148 + r_C(2,2)*r_DN(2,0);
const double crLHS150 = crLHS105*crLHS139;
const double crLHS151 = -crLHS150 + r_DN(2,1);
const double crLHS152 = crLHS108*r_DN(2,0) + crLHS109*r_DN(2,1);
const double crLHS153 = crLHS152 + crLHS51;
const double crLHS154 = r_DN(0,0)*r_N[3];
const double crLHS155 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS156 = r_C(0,2)*r_DN(3,0);
const double crLHS157 = crLHS156 + r_C(2,2)*r_DN(3,1);
const double crLHS158 = crLHS7*r_N[3];
const double crLHS159 = -crLHS158*crLHS84 + r_DN(3,0);
const double crLHS160 = crLHS158*crLHS4;
const double crLHS161 = crLHS15*crLHS160;
const double crLHS162 = crLHS161*crLHS8;
const double crLHS163 = crLHS125*crLHS68;
const double crLHS164 = -crLHS100*crLHS163 + crLHS124*crLHS68 + crLHS162 + crLHS163*crLHS98 + crLHS67*crLHS92;
const double crLHS165 = crLHS156 + r_C(0,1)*r_DN(3,1);
const double crLHS166 = r_C(1,2)*r_DN(3,1);
const double crLHS167 = crLHS166 + r_C(2,2)*r_DN(3,0);
const double crLHS168 = -crLHS105*crLHS158 + r_DN(3,1);
const double crLHS169 = crLHS108*r_DN(3,0) + crLHS109*r_DN(3,1);
const double crLHS170 = crLHS169 + crLHS66;
const double crLHS171 = crLHS103 + r_C(0,1)*r_DN(0,0);
const double crLHS172 = crLHS88*tau_c;
const double crLHS173 = crLHS172*r_DN(0,1);
const double crLHS174 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS175 = crLHS112*tau_c;
const double crLHS176 = crLHS175*r_DN(0,1);
const double crLHS177 = r_DN(0,1)*r_N[1];
const double crLHS178 = crLHS129 + r_C(0,1)*r_DN(1,0);
const double crLHS179 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS180 = r_DN(0,1)*r_N[2];
const double crLHS181 = crLHS148 + r_C(0,1)*r_DN(2,0);
const double crLHS182 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS183 = r_DN(0,1)*r_N[3];
const double crLHS184 = crLHS166 + r_C(0,1)*r_DN(3,0);
const double crLHS185 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS186 = crLHS14*crLHS5;
const double crLHS187 = -crLHS20 + dp_th_dt*r_N[0];
const double crLHS188 = crLHS77*crLHS99*tau_t;
const double crLHS189 = crLHS42*crLHS85;
const double crLHS190 = crLHS106*crLHS42;
const double crLHS191 = -crLHS38 + dp_th_dt*r_N[1];
const double crLHS192 = crLHS122 + crLHS28*kappa + crLHS29*kappa - crLHS7*dp_th_dt*r_N[0]*r_N[1];
const double crLHS193 = crLHS58*crLHS85;
const double crLHS194 = crLHS106*crLHS58;
const double crLHS195 = -crLHS54 + dp_th_dt*r_N[2];
const double crLHS196 = crLHS143 + crLHS45*kappa + crLHS46*kappa - crLHS7*dp_th_dt*r_N[0]*r_N[2];
const double crLHS197 = crLHS73*crLHS85;
const double crLHS198 = crLHS106*crLHS73;
const double crLHS199 = -crLHS69 + dp_th_dt*r_N[3];
const double crLHS200 = crLHS161 + crLHS60*kappa + crLHS61*kappa - crLHS7*dp_th_dt*r_N[0]*r_N[3];
const double crLHS201 = r_DN(1,0)*r_DN(1,0);
const double crLHS202 = r_DN(1,1)*r_DN(1,1);
const double crLHS203 = r_N[1]*r_N[1];
const double crLHS204 = crLHS203*crLHS7;
const double crLHS205 = r_DN(1,0)*r_DN(2,0);
const double crLHS206 = r_DN(1,1)*r_DN(2,1);
const double crLHS207 = crLHS11*(crLHS205 + crLHS206);
const double crLHS208 = r_DN(2,0)*r_N[1];
const double crLHS209 = crLHS118*r_N[2];
const double crLHS210 = -crLHS12*crLHS209;
const double crLHS211 = r_DN(2,1)*r_N[1];
const double crLHS212 = -crLHS209*crLHS23;
const double crLHS213 = crLHS26*crLHS35;
const double crLHS214 = -crLHS213*crLHS58;
const double crLHS215 = r_DN(1,0)*r_DN(3,0);
const double crLHS216 = r_DN(1,1)*r_DN(3,1);
const double crLHS217 = crLHS11*(crLHS215 + crLHS216);
const double crLHS218 = r_DN(3,0)*r_N[1];
const double crLHS219 = crLHS118*r_N[3];
const double crLHS220 = -crLHS12*crLHS219;
const double crLHS221 = r_DN(3,1)*r_N[1];
const double crLHS222 = -crLHS219*crLHS23;
const double crLHS223 = -crLHS213*crLHS73;
const double crLHS224 = crLHS172*r_DN(1,0);
const double crLHS225 = crLHS121*crLHS8;
const double crLHS226 = crLHS36*crLHS95;
const double crLHS227 = crLHS94*r_N[1];
const double crLHS228 = -crLHS100*crLHS227 + crLHS123 + crLHS18*crLHS225 + crLHS226*crLHS94 + crLHS227*crLHS98;
const double crLHS229 = crLHS175*r_DN(1,0);
const double crLHS230 = crLHS4*r_N[1];
const double crLHS231 = gauss_weight*(-crLHS230*crLHS79 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS232 = crLHS204*crLHS24;
const double crLHS233 = crLHS226*crLHS93;
const double crLHS234 = crLHS93*r_N[1];
const double crLHS235 = crLHS234*crLHS37;
const double crLHS236 = -crLHS100*crLHS235 + crLHS225*crLHS36 + crLHS232*crLHS8 + crLHS233*crLHS37 + crLHS235*crLHS98;
const double crLHS237 = r_DN(1,0)*r_N[2];
const double crLHS238 = crLHS142*crLHS35;
const double crLHS239 = crLHS238*crLHS8;
const double crLHS240 = crLHS234*crLHS53;
const double crLHS241 = -crLHS100*crLHS240 + crLHS225*crLHS52 + crLHS233*crLHS53 + crLHS239 + crLHS240*crLHS98;
const double crLHS242 = r_DN(1,0)*r_N[3];
const double crLHS243 = crLHS160*crLHS35;
const double crLHS244 = crLHS243*crLHS8;
const double crLHS245 = crLHS234*crLHS68;
const double crLHS246 = -crLHS100*crLHS245 + crLHS225*crLHS67 + crLHS233*crLHS68 + crLHS244 + crLHS245*crLHS98;
const double crLHS247 = crLHS172*r_DN(1,1);
const double crLHS248 = crLHS175*r_DN(1,1);
const double crLHS249 = r_DN(1,1)*r_N[2];
const double crLHS250 = r_DN(1,1)*r_N[3];
const double crLHS251 = crLHS204*crLHS5;
const double crLHS252 = crLHS119*crLHS58;
const double crLHS253 = crLHS131*crLHS58;
const double crLHS254 = crLHS205*kappa + crLHS206*kappa + crLHS238 - crLHS7*dp_th_dt*r_N[1]*r_N[2];
const double crLHS255 = crLHS119*crLHS73;
const double crLHS256 = crLHS131*crLHS73;
const double crLHS257 = crLHS215*kappa + crLHS216*kappa + crLHS243 - crLHS7*dp_th_dt*r_N[1]*r_N[3];
const double crLHS258 = r_DN(2,0)*r_DN(2,0);
const double crLHS259 = r_DN(2,1)*r_DN(2,1);
const double crLHS260 = r_N[2]*r_N[2];
const double crLHS261 = crLHS260*crLHS7;
const double crLHS262 = r_DN(2,0)*r_DN(3,0);
const double crLHS263 = r_DN(2,1)*r_DN(3,1);
const double crLHS264 = crLHS11*(crLHS262 + crLHS263);
const double crLHS265 = r_DN(3,0)*r_N[2];
const double crLHS266 = crLHS139*r_N[3];
const double crLHS267 = -crLHS12*crLHS266;
const double crLHS268 = r_DN(3,1)*r_N[2];
const double crLHS269 = -crLHS23*crLHS266;
const double crLHS270 = -crLHS26*crLHS51*crLHS73;
const double crLHS271 = crLHS172*r_DN(2,0);
const double crLHS272 = crLHS142*crLHS8;
const double crLHS273 = crLHS52*crLHS95;
const double crLHS274 = crLHS94*r_N[2];
const double crLHS275 = -crLHS100*crLHS274 + crLHS144 + crLHS18*crLHS272 + crLHS273*crLHS94 + crLHS274*crLHS98;
const double crLHS276 = crLHS175*r_DN(2,0);
const double crLHS277 = crLHS273*crLHS93;
const double crLHS278 = crLHS93*r_N[2];
const double crLHS279 = crLHS278*crLHS37;
const double crLHS280 = -crLHS100*crLHS279 + crLHS239 + crLHS272*crLHS36 + crLHS277*crLHS37 + crLHS279*crLHS98;
const double crLHS281 = crLHS4*r_N[2];
const double crLHS282 = gauss_weight*(-crLHS281*crLHS79 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS283 = crLHS24*crLHS261;
const double crLHS284 = crLHS278*crLHS53;
const double crLHS285 = -crLHS100*crLHS284 + crLHS272*crLHS52 + crLHS277*crLHS53 + crLHS283*crLHS8 + crLHS284*crLHS98;
const double crLHS286 = r_DN(2,0)*r_N[3];
const double crLHS287 = crLHS160*crLHS51;
const double crLHS288 = crLHS287*crLHS8;
const double crLHS289 = crLHS278*crLHS68;
const double crLHS290 = -crLHS100*crLHS289 + crLHS272*crLHS67 + crLHS277*crLHS68 + crLHS288 + crLHS289*crLHS98;
const double crLHS291 = crLHS172*r_DN(2,1);
const double crLHS292 = crLHS175*r_DN(2,1);
const double crLHS293 = r_DN(2,1)*r_N[3];
const double crLHS294 = crLHS261*crLHS5;
const double crLHS295 = crLHS140*crLHS73;
const double crLHS296 = crLHS150*crLHS73;
const double crLHS297 = crLHS262*kappa + crLHS263*kappa + crLHS287 - crLHS7*dp_th_dt*r_N[2]*r_N[3];
const double crLHS298 = r_DN(3,0)*r_DN(3,0);
const double crLHS299 = r_DN(3,1)*r_DN(3,1);
const double crLHS300 = r_N[3]*r_N[3];
const double crLHS301 = crLHS300*crLHS7;
const double crLHS302 = crLHS172*r_DN(3,0);
const double crLHS303 = crLHS160*crLHS8;
const double crLHS304 = crLHS67*crLHS95;
const double crLHS305 = crLHS94*r_N[3];
const double crLHS306 = -crLHS100*crLHS305 + crLHS162 + crLHS18*crLHS303 + crLHS304*crLHS94 + crLHS305*crLHS98;
const double crLHS307 = crLHS175*r_DN(3,0);
const double crLHS308 = crLHS304*crLHS93;
const double crLHS309 = crLHS93*r_N[3];
const double crLHS310 = crLHS309*crLHS37;
const double crLHS311 = -crLHS100*crLHS310 + crLHS244 + crLHS303*crLHS36 + crLHS308*crLHS37 + crLHS310*crLHS98;
const double crLHS312 = crLHS309*crLHS53;
const double crLHS313 = -crLHS100*crLHS312 + crLHS288 + crLHS303*crLHS52 + crLHS308*crLHS53 + crLHS312*crLHS98;
const double crLHS314 = crLHS4*r_N[3];
const double crLHS315 = gauss_weight*(crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_N[3]*tau_u - crLHS314*crLHS79 - r_N[3]);
const double crLHS316 = crLHS24*crLHS301;
const double crLHS317 = crLHS309*crLHS68;
const double crLHS318 = -crLHS100*crLHS317 + crLHS303*crLHS67 + crLHS308*crLHS68 + crLHS316*crLHS8 + crLHS317*crLHS98;
const double crLHS319 = crLHS172*r_DN(3,1);
const double crLHS320 = crLHS175*r_DN(3,1);
const double crLHS321 = crLHS301*crLHS5;
rLHS(0,0)+=crLHS11*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS10*(-crLHS12*crLHS14 + crLHS22*r_DN(0,0) + r_DN(0,0)*r_N[0]);
rLHS(0,2)+=crLHS10*(-crLHS14*crLHS23 + crLHS22*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS13*crLHS27;
rLHS(0,4)+=crLHS30;
rLHS(0,5)+=crLHS10*(crLHS31 + crLHS34 + crLHS39*r_DN(0,0));
rLHS(0,6)+=crLHS10*(crLHS39*r_DN(0,1) + crLHS40 + crLHS41);
rLHS(0,7)+=crLHS44;
rLHS(0,8)+=crLHS47;
rLHS(0,9)+=crLHS10*(crLHS48 + crLHS50 + crLHS55*r_DN(0,0));
rLHS(0,10)+=crLHS10*(crLHS55*r_DN(0,1) + crLHS56 + crLHS57);
rLHS(0,11)+=crLHS59;
rLHS(0,12)+=crLHS62;
rLHS(0,13)+=crLHS10*(crLHS63 + crLHS65 + crLHS70*r_DN(0,0));
rLHS(0,14)+=crLHS10*(crLHS70*r_DN(0,1) + crLHS71 + crLHS72);
rLHS(0,15)+=crLHS74;
rLHS(1,0)+=crLHS80*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS101 + crLHS81*r_DN(0,0) + crLHS83*r_DN(0,1) + crLHS86*crLHS89);
rLHS(1,2)+=gauss_weight*(crLHS102*r_DN(0,0) + crLHS104*r_DN(0,1) + crLHS107*crLHS89);
rLHS(1,3)+=-crLHS111*crLHS113;
rLHS(1,4)+=gauss_weight*(-crLHS114 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS31*crLHS4*crLHS79);
rLHS(1,5)+=gauss_weight*(crLHS115*r_DN(0,0) + crLHS117*r_DN(0,1) + crLHS120*crLHS89 + crLHS127);
rLHS(1,6)+=gauss_weight*(crLHS128*r_DN(0,0) + crLHS130*r_DN(0,1) + crLHS132*crLHS89);
rLHS(1,7)+=-crLHS113*crLHS134;
rLHS(1,8)+=gauss_weight*(-crLHS135 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS4*crLHS48*crLHS79);
rLHS(1,9)+=gauss_weight*(crLHS136*r_DN(0,0) + crLHS138*r_DN(0,1) + crLHS141*crLHS89 + crLHS146);
rLHS(1,10)+=gauss_weight*(crLHS147*r_DN(0,0) + crLHS149*r_DN(0,1) + crLHS151*crLHS89);
rLHS(1,11)+=-crLHS113*crLHS153;
rLHS(1,12)+=gauss_weight*(-crLHS154 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,0)*r_N[0]*tau_u - crLHS4*crLHS63*crLHS79);
rLHS(1,13)+=gauss_weight*(crLHS155*r_DN(0,0) + crLHS157*r_DN(0,1) + crLHS159*crLHS89 + crLHS164);
rLHS(1,14)+=gauss_weight*(crLHS165*r_DN(0,0) + crLHS167*r_DN(0,1) + crLHS168*crLHS89);
rLHS(1,15)+=-crLHS113*crLHS170;
rLHS(2,0)+=crLHS80*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS171*r_DN(0,1) + crLHS173*crLHS86 + crLHS83*r_DN(0,0));
rLHS(2,2)+=gauss_weight*(crLHS101 + crLHS104*r_DN(0,0) + crLHS107*crLHS173 + crLHS174*r_DN(0,1));
rLHS(2,3)+=-crLHS111*crLHS176;
rLHS(2,4)+=gauss_weight*(-crLHS177 + crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS4*crLHS40*crLHS79);
rLHS(2,5)+=gauss_weight*(crLHS117*r_DN(0,0) + crLHS120*crLHS173 + crLHS178*r_DN(0,1));
rLHS(2,6)+=gauss_weight*(crLHS127 + crLHS130*r_DN(0,0) + crLHS132*crLHS173 + crLHS179*r_DN(0,1));
rLHS(2,7)+=-crLHS134*crLHS176;
rLHS(2,8)+=gauss_weight*(crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*tau_u - crLHS180 + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS4*crLHS56*crLHS79);
rLHS(2,9)+=gauss_weight*(crLHS138*r_DN(0,0) + crLHS141*crLHS173 + crLHS181*r_DN(0,1));
rLHS(2,10)+=gauss_weight*(crLHS146 + crLHS149*r_DN(0,0) + crLHS151*crLHS173 + crLHS182*r_DN(0,1));
rLHS(2,11)+=-crLHS153*crLHS176;
rLHS(2,12)+=gauss_weight*(crLHS18*crLHS3*crLHS7*crLHS8*gamma*p_th*r_DN(3,1)*tau_u - crLHS183 + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,1)*r_N[0]*tau_u - crLHS4*crLHS71*crLHS79);
rLHS(2,13)+=gauss_weight*(crLHS157*r_DN(0,0) + crLHS159*crLHS173 + crLHS184*r_DN(0,1));
rLHS(2,14)+=gauss_weight*(crLHS164 + crLHS167*r_DN(0,0) + crLHS168*crLHS173 + crLHS185*r_DN(0,1));
rLHS(2,15)+=-crLHS170*crLHS176;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS186*crLHS84;
rLHS(3,2)+=crLHS105*crLHS186;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa - crLHS110*crLHS91 + crLHS13*crLHS7*dp_th_dt + crLHS18*crLHS187*crLHS25*crLHS3*gamma*p_th*tau_t - crLHS187*crLHS188*crLHS76 + crLHS187*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS187*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS90);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS189;
rLHS(3,6)+=crLHS190;
rLHS(3,7)+=-gauss_weight*(-crLHS133*crLHS91 + crLHS18*crLHS191*crLHS25*crLHS3*gamma*p_th*tau_t - crLHS188*crLHS191*crLHS76 + crLHS191*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS191*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS192);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS193;
rLHS(3,10)+=crLHS194;
rLHS(3,11)+=-gauss_weight*(-crLHS152*crLHS91 + crLHS18*crLHS195*crLHS25*crLHS3*gamma*p_th*tau_t - crLHS188*crLHS195*crLHS76 + crLHS195*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS195*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS196);
rLHS(3,12)+=0;
rLHS(3,13)+=crLHS197;
rLHS(3,14)+=crLHS198;
rLHS(3,15)+=-gauss_weight*(-crLHS169*crLHS91 + crLHS18*crLHS199*crLHS25*crLHS3*gamma*p_th*tau_t - crLHS188*crLHS199*crLHS76 + crLHS199*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[0]*tau_t + crLHS199*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS200);
rLHS(4,0)+=crLHS30;
rLHS(4,1)+=crLHS10*(crLHS114 + crLHS22*r_DN(1,0) + crLHS34);
rLHS(4,2)+=crLHS10*(crLHS177 + crLHS22*r_DN(1,1) + crLHS41);
rLHS(4,3)+=crLHS44;
rLHS(4,4)+=crLHS11*(crLHS201 + crLHS202);
rLHS(4,5)+=crLHS10*(-crLHS12*crLHS204 + crLHS39*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS204*crLHS23 + crLHS39*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS203*crLHS27;
rLHS(4,8)+=crLHS207;
rLHS(4,9)+=crLHS10*(crLHS208 + crLHS210 + crLHS55*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS211 + crLHS212 + crLHS55*r_DN(1,1));
rLHS(4,11)+=crLHS214;
rLHS(4,12)+=crLHS217;
rLHS(4,13)+=crLHS10*(crLHS218 + crLHS220 + crLHS70*r_DN(1,0));
rLHS(4,14)+=crLHS10*(crLHS221 + crLHS222 + crLHS70*r_DN(1,1));
rLHS(4,15)+=crLHS223;
rLHS(5,0)+=gauss_weight*(-crLHS114*crLHS4*crLHS79 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS31);
rLHS(5,1)+=gauss_weight*(crLHS224*crLHS86 + crLHS228 + crLHS81*r_DN(1,0) + crLHS83*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS102*r_DN(1,0) + crLHS104*r_DN(1,1) + crLHS107*crLHS224);
rLHS(5,3)+=-crLHS111*crLHS229;
rLHS(5,4)+=crLHS231*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS115*r_DN(1,0) + crLHS117*r_DN(1,1) + crLHS120*crLHS224 + crLHS236);
rLHS(5,6)+=gauss_weight*(crLHS128*r_DN(1,0) + crLHS130*r_DN(1,1) + crLHS132*crLHS224);
rLHS(5,7)+=-crLHS134*crLHS229;
rLHS(5,8)+=gauss_weight*(-crLHS208*crLHS4*crLHS79 - crLHS237 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS136*r_DN(1,0) + crLHS138*r_DN(1,1) + crLHS141*crLHS224 + crLHS241);
rLHS(5,10)+=gauss_weight*(crLHS147*r_DN(1,0) + crLHS149*r_DN(1,1) + crLHS151*crLHS224);
rLHS(5,11)+=-crLHS153*crLHS229;
rLHS(5,12)+=gauss_weight*(-crLHS218*crLHS4*crLHS79 - crLHS242 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,0)*r_N[1]*tau_u);
rLHS(5,13)+=gauss_weight*(crLHS155*r_DN(1,0) + crLHS157*r_DN(1,1) + crLHS159*crLHS224 + crLHS246);
rLHS(5,14)+=gauss_weight*(crLHS165*r_DN(1,0) + crLHS167*r_DN(1,1) + crLHS168*crLHS224);
rLHS(5,15)+=-crLHS170*crLHS229;
rLHS(6,0)+=gauss_weight*(-crLHS177*crLHS4*crLHS79 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS40);
rLHS(6,1)+=gauss_weight*(crLHS171*r_DN(1,1) + crLHS247*crLHS86 + crLHS83*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS104*r_DN(1,0) + crLHS107*crLHS247 + crLHS174*r_DN(1,1) + crLHS228);
rLHS(6,3)+=-crLHS111*crLHS248;
rLHS(6,4)+=crLHS231*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS117*r_DN(1,0) + crLHS120*crLHS247 + crLHS178*r_DN(1,1));
rLHS(6,6)+=gauss_weight*(crLHS130*r_DN(1,0) + crLHS132*crLHS247 + crLHS179*r_DN(1,1) + crLHS236);
rLHS(6,7)+=-crLHS134*crLHS248;
rLHS(6,8)+=gauss_weight*(-crLHS211*crLHS4*crLHS79 - crLHS249 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS138*r_DN(1,0) + crLHS141*crLHS247 + crLHS181*r_DN(1,1));
rLHS(6,10)+=gauss_weight*(crLHS149*r_DN(1,0) + crLHS151*crLHS247 + crLHS182*r_DN(1,1) + crLHS241);
rLHS(6,11)+=-crLHS153*crLHS248;
rLHS(6,12)+=gauss_weight*(-crLHS221*crLHS4*crLHS79 - crLHS250 + crLHS3*crLHS36*crLHS7*crLHS8*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,1)*r_N[1]*tau_u);
rLHS(6,13)+=gauss_weight*(crLHS157*r_DN(1,0) + crLHS159*crLHS247 + crLHS184*r_DN(1,1));
rLHS(6,14)+=gauss_weight*(crLHS167*r_DN(1,0) + crLHS168*crLHS247 + crLHS185*r_DN(1,1) + crLHS246);
rLHS(6,15)+=-crLHS170*crLHS248;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS189;
rLHS(7,2)+=crLHS190;
rLHS(7,3)+=-gauss_weight*(-crLHS110*crLHS121 - crLHS187*crLHS188*crLHS230 + crLHS187*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS187*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS187*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS192);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS251*crLHS84;
rLHS(7,6)+=crLHS105*crLHS251;
rLHS(7,7)+=-gauss_weight*(-crLHS121*crLHS133 - crLHS188*crLHS191*crLHS230 + crLHS191*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS191*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS191*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS201*kappa - crLHS202*kappa + crLHS203*crLHS7*dp_th_dt - crLHS232);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS252;
rLHS(7,10)+=crLHS253;
rLHS(7,11)+=-gauss_weight*(-crLHS121*crLHS152 - crLHS188*crLHS195*crLHS230 + crLHS195*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS195*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS195*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS254);
rLHS(7,12)+=0;
rLHS(7,13)+=crLHS255;
rLHS(7,14)+=crLHS256;
rLHS(7,15)+=-gauss_weight*(-crLHS121*crLHS169 - crLHS188*crLHS199*crLHS230 + crLHS199*crLHS25*crLHS3*crLHS36*gamma*p_th*tau_t + crLHS199*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[1]*tau_t + crLHS199*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS257);
rLHS(8,0)+=crLHS47;
rLHS(8,1)+=crLHS10*(crLHS135 + crLHS22*r_DN(2,0) + crLHS50);
rLHS(8,2)+=crLHS10*(crLHS180 + crLHS22*r_DN(2,1) + crLHS57);
rLHS(8,3)+=crLHS59;
rLHS(8,4)+=crLHS207;
rLHS(8,5)+=crLHS10*(crLHS210 + crLHS237 + crLHS39*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS212 + crLHS249 + crLHS39*r_DN(2,1));
rLHS(8,7)+=crLHS214;
rLHS(8,8)+=crLHS11*(crLHS258 + crLHS259);
rLHS(8,9)+=crLHS10*(-crLHS12*crLHS261 + crLHS55*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS23*crLHS261 + crLHS55*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS260*crLHS27;
rLHS(8,12)+=crLHS264;
rLHS(8,13)+=crLHS10*(crLHS265 + crLHS267 + crLHS70*r_DN(2,0));
rLHS(8,14)+=crLHS10*(crLHS268 + crLHS269 + crLHS70*r_DN(2,1));
rLHS(8,15)+=crLHS270;
rLHS(9,0)+=gauss_weight*(-crLHS135*crLHS4*crLHS79 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS48);
rLHS(9,1)+=gauss_weight*(crLHS271*crLHS86 + crLHS275 + crLHS81*r_DN(2,0) + crLHS83*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS102*r_DN(2,0) + crLHS104*r_DN(2,1) + crLHS107*crLHS271);
rLHS(9,3)+=-crLHS111*crLHS276;
rLHS(9,4)+=gauss_weight*(-crLHS208 - crLHS237*crLHS4*crLHS79 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS115*r_DN(2,0) + crLHS117*r_DN(2,1) + crLHS120*crLHS271 + crLHS280);
rLHS(9,6)+=gauss_weight*(crLHS128*r_DN(2,0) + crLHS130*r_DN(2,1) + crLHS132*crLHS271);
rLHS(9,7)+=-crLHS134*crLHS276;
rLHS(9,8)+=crLHS282*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS136*r_DN(2,0) + crLHS138*r_DN(2,1) + crLHS141*crLHS271 + crLHS285);
rLHS(9,10)+=gauss_weight*(crLHS147*r_DN(2,0) + crLHS149*r_DN(2,1) + crLHS151*crLHS271);
rLHS(9,11)+=-crLHS153*crLHS276;
rLHS(9,12)+=gauss_weight*(-crLHS265*crLHS4*crLHS79 - crLHS286 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,0)*r_N[2]*tau_u);
rLHS(9,13)+=gauss_weight*(crLHS155*r_DN(2,0) + crLHS157*r_DN(2,1) + crLHS159*crLHS271 + crLHS290);
rLHS(9,14)+=gauss_weight*(crLHS165*r_DN(2,0) + crLHS167*r_DN(2,1) + crLHS168*crLHS271);
rLHS(9,15)+=-crLHS170*crLHS276;
rLHS(10,0)+=gauss_weight*(-crLHS180*crLHS4*crLHS79 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS56);
rLHS(10,1)+=gauss_weight*(crLHS171*r_DN(2,1) + crLHS291*crLHS86 + crLHS83*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS104*r_DN(2,0) + crLHS107*crLHS291 + crLHS174*r_DN(2,1) + crLHS275);
rLHS(10,3)+=-crLHS111*crLHS292;
rLHS(10,4)+=gauss_weight*(-crLHS211 - crLHS249*crLHS4*crLHS79 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS117*r_DN(2,0) + crLHS120*crLHS291 + crLHS178*r_DN(2,1));
rLHS(10,6)+=gauss_weight*(crLHS130*r_DN(2,0) + crLHS132*crLHS291 + crLHS179*r_DN(2,1) + crLHS280);
rLHS(10,7)+=-crLHS134*crLHS292;
rLHS(10,8)+=crLHS282*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS138*r_DN(2,0) + crLHS141*crLHS291 + crLHS181*r_DN(2,1));
rLHS(10,10)+=gauss_weight*(crLHS149*r_DN(2,0) + crLHS151*crLHS291 + crLHS182*r_DN(2,1) + crLHS285);
rLHS(10,11)+=-crLHS153*crLHS292;
rLHS(10,12)+=gauss_weight*(-crLHS268*crLHS4*crLHS79 - crLHS293 + crLHS3*crLHS52*crLHS7*crLHS8*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(3,1)*r_N[2]*tau_u);
rLHS(10,13)+=gauss_weight*(crLHS157*r_DN(2,0) + crLHS159*crLHS291 + crLHS184*r_DN(2,1));
rLHS(10,14)+=gauss_weight*(crLHS167*r_DN(2,0) + crLHS168*crLHS291 + crLHS185*r_DN(2,1) + crLHS290);
rLHS(10,15)+=-crLHS170*crLHS292;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS193;
rLHS(11,2)+=crLHS194;
rLHS(11,3)+=-gauss_weight*(-crLHS110*crLHS142 - crLHS187*crLHS188*crLHS281 + crLHS187*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS187*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS187*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS196);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS252;
rLHS(11,6)+=crLHS253;
rLHS(11,7)+=-gauss_weight*(-crLHS133*crLHS142 - crLHS188*crLHS191*crLHS281 + crLHS191*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS191*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS191*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS254);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS294*crLHS84;
rLHS(11,10)+=crLHS105*crLHS294;
rLHS(11,11)+=-gauss_weight*(-crLHS142*crLHS152 - crLHS188*crLHS195*crLHS281 + crLHS195*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS195*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS195*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS258*kappa - crLHS259*kappa + crLHS260*crLHS7*dp_th_dt - crLHS283);
rLHS(11,12)+=0;
rLHS(11,13)+=crLHS295;
rLHS(11,14)+=crLHS296;
rLHS(11,15)+=-gauss_weight*(-crLHS142*crLHS169 - crLHS188*crLHS199*crLHS281 + crLHS199*crLHS25*crLHS3*crLHS52*gamma*p_th*tau_t + crLHS199*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[2]*tau_t + crLHS199*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS297);
rLHS(12,0)+=crLHS62;
rLHS(12,1)+=crLHS10*(crLHS154 + crLHS22*r_DN(3,0) + crLHS65);
rLHS(12,2)+=crLHS10*(crLHS183 + crLHS22*r_DN(3,1) + crLHS72);
rLHS(12,3)+=crLHS74;
rLHS(12,4)+=crLHS217;
rLHS(12,5)+=crLHS10*(crLHS220 + crLHS242 + crLHS39*r_DN(3,0));
rLHS(12,6)+=crLHS10*(crLHS222 + crLHS250 + crLHS39*r_DN(3,1));
rLHS(12,7)+=crLHS223;
rLHS(12,8)+=crLHS264;
rLHS(12,9)+=crLHS10*(crLHS267 + crLHS286 + crLHS55*r_DN(3,0));
rLHS(12,10)+=crLHS10*(crLHS269 + crLHS293 + crLHS55*r_DN(3,1));
rLHS(12,11)+=crLHS270;
rLHS(12,12)+=crLHS11*(crLHS298 + crLHS299);
rLHS(12,13)+=crLHS10*(-crLHS12*crLHS301 + crLHS70*r_DN(3,0) + r_DN(3,0)*r_N[3]);
rLHS(12,14)+=crLHS10*(-crLHS23*crLHS301 + crLHS70*r_DN(3,1) + r_DN(3,1)*r_N[3]);
rLHS(12,15)+=-crLHS27*crLHS300;
rLHS(13,0)+=gauss_weight*(-crLHS154*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,0)*r_N[3]*tau_u - crLHS63);
rLHS(13,1)+=gauss_weight*(crLHS302*crLHS86 + crLHS306 + crLHS81*r_DN(3,0) + crLHS83*r_DN(3,1));
rLHS(13,2)+=gauss_weight*(crLHS102*r_DN(3,0) + crLHS104*r_DN(3,1) + crLHS107*crLHS302);
rLHS(13,3)+=-crLHS111*crLHS307;
rLHS(13,4)+=gauss_weight*(-crLHS218 - crLHS242*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,0)*r_N[3]*tau_u);
rLHS(13,5)+=gauss_weight*(crLHS115*r_DN(3,0) + crLHS117*r_DN(3,1) + crLHS120*crLHS302 + crLHS311);
rLHS(13,6)+=gauss_weight*(crLHS128*r_DN(3,0) + crLHS130*r_DN(3,1) + crLHS132*crLHS302);
rLHS(13,7)+=-crLHS134*crLHS307;
rLHS(13,8)+=gauss_weight*(-crLHS265 - crLHS286*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,0)*r_N[3]*tau_u);
rLHS(13,9)+=gauss_weight*(crLHS136*r_DN(3,0) + crLHS138*r_DN(3,1) + crLHS141*crLHS302 + crLHS313);
rLHS(13,10)+=gauss_weight*(crLHS147*r_DN(3,0) + crLHS149*r_DN(3,1) + crLHS151*crLHS302);
rLHS(13,11)+=-crLHS153*crLHS307;
rLHS(13,12)+=crLHS315*r_DN(3,0);
rLHS(13,13)+=gauss_weight*(crLHS155*r_DN(3,0) + crLHS157*r_DN(3,1) + crLHS159*crLHS302 + crLHS318);
rLHS(13,14)+=gauss_weight*(crLHS165*r_DN(3,0) + crLHS167*r_DN(3,1) + crLHS168*crLHS302);
rLHS(13,15)+=-crLHS170*crLHS307;
rLHS(14,0)+=gauss_weight*(-crLHS183*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(0,1)*r_N[3]*tau_u - crLHS71);
rLHS(14,1)+=gauss_weight*(crLHS171*r_DN(3,1) + crLHS319*crLHS86 + crLHS83*r_DN(3,0));
rLHS(14,2)+=gauss_weight*(crLHS104*r_DN(3,0) + crLHS107*crLHS319 + crLHS174*r_DN(3,1) + crLHS306);
rLHS(14,3)+=-crLHS111*crLHS320;
rLHS(14,4)+=gauss_weight*(-crLHS221 - crLHS250*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(1,1)*r_N[3]*tau_u);
rLHS(14,5)+=gauss_weight*(crLHS117*r_DN(3,0) + crLHS120*crLHS319 + crLHS178*r_DN(3,1));
rLHS(14,6)+=gauss_weight*(crLHS130*r_DN(3,0) + crLHS132*crLHS319 + crLHS179*r_DN(3,1) + crLHS311);
rLHS(14,7)+=-crLHS134*crLHS320;
rLHS(14,8)+=gauss_weight*(-crLHS268 - crLHS293*crLHS4*crLHS79 + crLHS3*crLHS67*crLHS7*crLHS8*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS7*crLHS75*crLHS8*gamma*p_th*r_DN(2,1)*r_N[3]*tau_u);
rLHS(14,9)+=gauss_weight*(crLHS138*r_DN(3,0) + crLHS141*crLHS319 + crLHS181*r_DN(3,1));
rLHS(14,10)+=gauss_weight*(crLHS149*r_DN(3,0) + crLHS151*crLHS319 + crLHS182*r_DN(3,1) + crLHS313);
rLHS(14,11)+=-crLHS153*crLHS320;
rLHS(14,12)+=crLHS315*r_DN(3,1);
rLHS(14,13)+=gauss_weight*(crLHS157*r_DN(3,0) + crLHS159*crLHS319 + crLHS184*r_DN(3,1));
rLHS(14,14)+=gauss_weight*(crLHS167*r_DN(3,0) + crLHS168*crLHS319 + crLHS185*r_DN(3,1) + crLHS318);
rLHS(14,15)+=-crLHS170*crLHS320;
rLHS(15,0)+=0;
rLHS(15,1)+=crLHS197;
rLHS(15,2)+=crLHS198;
rLHS(15,3)+=-gauss_weight*(-crLHS110*crLHS160 - crLHS187*crLHS188*crLHS314 + crLHS187*crLHS25*crLHS3*crLHS67*gamma*p_th*tau_t + crLHS187*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS187*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS200);
rLHS(15,4)+=0;
rLHS(15,5)+=crLHS255;
rLHS(15,6)+=crLHS256;
rLHS(15,7)+=-gauss_weight*(-crLHS133*crLHS160 - crLHS188*crLHS191*crLHS314 + crLHS191*crLHS25*crLHS3*crLHS67*gamma*p_th*tau_t + crLHS191*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS191*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS257);
rLHS(15,8)+=0;
rLHS(15,9)+=crLHS295;
rLHS(15,10)+=crLHS296;
rLHS(15,11)+=-gauss_weight*(-crLHS152*crLHS160 - crLHS188*crLHS195*crLHS314 + crLHS195*crLHS25*crLHS3*crLHS67*gamma*p_th*tau_t + crLHS195*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS195*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS297);
rLHS(15,12)+=0;
rLHS(15,13)+=crLHS321*crLHS84;
rLHS(15,14)+=crLHS105*crLHS321;
rLHS(15,15)+=-gauss_weight*(-crLHS160*crLHS169 - crLHS188*crLHS199*crLHS314 + crLHS199*crLHS25*crLHS3*crLHS67*gamma*p_th*tau_t + crLHS199*crLHS25*crLHS3*crLHS75*gamma*p_th*r_N[3]*tau_t + crLHS199*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS298*kappa - crLHS299*kappa + crLHS300*crLHS7*dp_th_dt - crLHS316);

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
const double crRHS14 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2];
const double crRHS15 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2];
const double crRHS16 = crRHS0*crRHS14 + crRHS15*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS17 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS18 = 1.0/c_p;
const double crRHS19 = gamma/(gamma - 1.0);
const double crRHS20 = crRHS19*crRHS6;
const double crRHS21 = crRHS18*crRHS20;
const double crRHS22 = -crRHS21*(-crRHS16 - crRHS17 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS23 = p_th*tau_u;
const double crRHS24 = crRHS22*crRHS23;
const double crRHS25 = crRHS1*crRHS15 + crRHS14*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS26 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS27 = -crRHS21*(-crRHS25 - crRHS26 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS28 = crRHS23*crRHS27;
const double crRHS29 = crRHS18*crRHS19;
const double crRHS30 = crRHS29*crRHS5;
const double crRHS31 = crRHS30*gauss_weight;
const double crRHS32 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS33 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS34 = crRHS21*r_N[0];
const double crRHS35 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crRHS36 = crRHS35*tau_u;
const double crRHS37 = crRHS34*crRHS36;
const double crRHS38 = crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1);
const double crRHS39 = crRHS21*tau_u;
const double crRHS40 = crRHS38*crRHS39;
const double crRHS41 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crRHS42 = crRHS10*crRHS41;
const double crRHS43 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crRHS44 = crRHS12*crRHS43;
const double crRHS45 = crRHS30*tau_c*(-crRHS2 + crRHS42*crRHS6 + crRHS44*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS46 = (crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS47 = crRHS46*r_N[0];
const double crRHS48 = crRHS29*crRHS47;
const double crRHS49 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS50 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2];
const double crRHS51 = crRHS41*kappa;
const double crRHS52 = crRHS43*kappa;
const double crRHS53 = crRHS5*dp_th_dt;
const double crRHS54 = crRHS53*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2]);
const double crRHS55 = crRHS19*crRHS7;
const double crRHS56 = crRHS20*(crRHS42 + crRHS44);
const double crRHS57 = tau_t*(-crRHS20*(crRHS14*crRHS41 + crRHS15*crRHS43 + crRHS3) + crRHS50 + crRHS54);
const double crRHS58 = crRHS53*crRHS57;
const double crRHS59 = crRHS20*crRHS57;
const double crRHS60 = crRHS35*crRHS59;
const double crRHS61 = crRHS19*crRHS57*p_th;
const double crRHS62 = crRHS21*r_N[1];
const double crRHS63 = crRHS36*crRHS62;
const double crRHS64 = crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1);
const double crRHS65 = crRHS39*crRHS64;
const double crRHS66 = crRHS29*crRHS46;
const double crRHS67 = crRHS66*r_N[1];
const double crRHS68 = crRHS46*crRHS61;
const double crRHS69 = crRHS21*r_N[2];
const double crRHS70 = crRHS36*crRHS69;
const double crRHS71 = crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1);
const double crRHS72 = crRHS39*crRHS71;
const double crRHS73 = crRHS66*r_N[2];
rRHS[0]+=-crRHS31*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS24*r_DN(0,0) + crRHS28*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS16*crRHS34 + crRHS17*crRHS34 + crRHS22*crRHS37 + crRHS22*crRHS40 - crRHS24*crRHS48 - crRHS32*r_DN(0,0) - crRHS33*crRHS34 - crRHS45*r_DN(0,0) + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS25*crRHS34 + crRHS26*crRHS34 + crRHS27*crRHS37 + crRHS27*crRHS40 - crRHS28*crRHS48 - crRHS32*r_DN(0,1) - crRHS34*crRHS49 - crRHS45*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS38*crRHS59 - crRHS47*crRHS61 + crRHS50*r_N[0] - crRHS51*r_DN(0,0) - crRHS52*r_DN(0,1) + crRHS54*r_N[0] - crRHS55*r_N[0] - crRHS56*r_N[0] + crRHS58*r_N[0] + crRHS60*r_N[0]);
rRHS[4]+=-crRHS31*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS24*r_DN(1,0) + crRHS28*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS16*crRHS62 + crRHS17*crRHS62 + crRHS22*crRHS63 + crRHS22*crRHS65 - crRHS24*crRHS67 - crRHS32*r_DN(1,0) - crRHS33*crRHS62 - crRHS45*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS25*crRHS62 + crRHS26*crRHS62 + crRHS27*crRHS63 + crRHS27*crRHS65 - crRHS28*crRHS67 - crRHS32*r_DN(1,1) - crRHS45*r_DN(1,1) - crRHS49*crRHS62 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS50*r_N[1] - crRHS51*r_DN(1,0) - crRHS52*r_DN(1,1) + crRHS54*r_N[1] - crRHS55*r_N[1] - crRHS56*r_N[1] + crRHS58*r_N[1] + crRHS59*crRHS64 + crRHS60*r_N[1] - crRHS68*r_N[1]);
rRHS[8]+=-crRHS31*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS24*r_DN(2,0) + crRHS28*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS16*crRHS69 + crRHS17*crRHS69 + crRHS22*crRHS70 + crRHS22*crRHS72 - crRHS24*crRHS73 - crRHS32*r_DN(2,0) - crRHS33*crRHS69 - crRHS45*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS25*crRHS69 + crRHS26*crRHS69 + crRHS27*crRHS70 + crRHS27*crRHS72 - crRHS28*crRHS73 - crRHS32*r_DN(2,1) - crRHS45*r_DN(2,1) - crRHS49*crRHS69 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS50*r_N[2] - crRHS51*r_DN(2,0) - crRHS52*r_DN(2,1) + crRHS54*r_N[2] - crRHS55*r_N[2] - crRHS56*r_N[2] + crRHS58*r_N[2] + crRHS59*crRHS71 + crRHS60*r_N[2] - crRHS68*r_N[2]);

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
const double crRHS14 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2] + lin_u_conv(3,0)*r_N[3];
const double crRHS15 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2] + lin_u_conv(3,1)*r_N[3];
const double crRHS16 = crRHS0*crRHS14 + crRHS15*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS17 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0)) + r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS18 = 1.0/c_p;
const double crRHS19 = gamma/(gamma - 1.0);
const double crRHS20 = crRHS19*crRHS6;
const double crRHS21 = crRHS18*crRHS20;
const double crRHS22 = -crRHS21*(-crRHS16 - crRHS17 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS23 = p_th*tau_u;
const double crRHS24 = crRHS22*crRHS23;
const double crRHS25 = crRHS1*crRHS15 + crRHS14*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS26 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1)) + r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS27 = -crRHS21*(-crRHS25 - crRHS26 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS28 = crRHS23*crRHS27;
const double crRHS29 = crRHS18*crRHS19;
const double crRHS30 = crRHS29*crRHS5;
const double crRHS31 = crRHS30*gauss_weight;
const double crRHS32 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS33 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS34 = crRHS21*r_N[0];
const double crRHS35 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crRHS36 = crRHS35*tau_u;
const double crRHS37 = crRHS34*crRHS36;
const double crRHS38 = crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1);
const double crRHS39 = crRHS21*tau_u;
const double crRHS40 = crRHS38*crRHS39;
const double crRHS41 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crRHS42 = crRHS10*crRHS41;
const double crRHS43 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crRHS44 = crRHS12*crRHS43;
const double crRHS45 = crRHS30*tau_c*(-crRHS2 + crRHS42*crRHS6 + crRHS44*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS46 = (crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS47 = crRHS46*r_N[0];
const double crRHS48 = crRHS29*crRHS47;
const double crRHS49 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS50 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2] + r_N[3]*r_heat_fl[3];
const double crRHS51 = crRHS41*kappa;
const double crRHS52 = crRHS43*kappa;
const double crRHS53 = crRHS5*dp_th_dt;
const double crRHS54 = crRHS53*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3]);
const double crRHS55 = crRHS19*crRHS7;
const double crRHS56 = crRHS20*(crRHS42 + crRHS44);
const double crRHS57 = tau_t*(-crRHS20*(crRHS14*crRHS41 + crRHS15*crRHS43 + crRHS3) + crRHS50 + crRHS54);
const double crRHS58 = crRHS53*crRHS57;
const double crRHS59 = crRHS20*crRHS57;
const double crRHS60 = crRHS35*crRHS59;
const double crRHS61 = crRHS19*crRHS57*p_th;
const double crRHS62 = crRHS21*r_N[1];
const double crRHS63 = crRHS36*crRHS62;
const double crRHS64 = crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1);
const double crRHS65 = crRHS39*crRHS64;
const double crRHS66 = crRHS29*crRHS46;
const double crRHS67 = crRHS66*r_N[1];
const double crRHS68 = crRHS46*crRHS61;
const double crRHS69 = crRHS21*r_N[2];
const double crRHS70 = crRHS36*crRHS69;
const double crRHS71 = crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1);
const double crRHS72 = crRHS39*crRHS71;
const double crRHS73 = crRHS66*r_N[2];
const double crRHS74 = crRHS21*r_N[3];
const double crRHS75 = crRHS36*crRHS74;
const double crRHS76 = crRHS14*r_DN(3,0) + crRHS15*r_DN(3,1);
const double crRHS77 = crRHS39*crRHS76;
const double crRHS78 = crRHS66*r_N[3];
rRHS[0]+=-crRHS31*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS24*r_DN(0,0) + crRHS28*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS16*crRHS34 + crRHS17*crRHS34 + crRHS22*crRHS37 + crRHS22*crRHS40 - crRHS24*crRHS48 - crRHS32*r_DN(0,0) - crRHS33*crRHS34 - crRHS45*r_DN(0,0) + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS25*crRHS34 + crRHS26*crRHS34 + crRHS27*crRHS37 + crRHS27*crRHS40 - crRHS28*crRHS48 - crRHS32*r_DN(0,1) - crRHS34*crRHS49 - crRHS45*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS38*crRHS59 - crRHS47*crRHS61 + crRHS50*r_N[0] - crRHS51*r_DN(0,0) - crRHS52*r_DN(0,1) + crRHS54*r_N[0] - crRHS55*r_N[0] - crRHS56*r_N[0] + crRHS58*r_N[0] + crRHS60*r_N[0]);
rRHS[4]+=-crRHS31*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS24*r_DN(1,0) + crRHS28*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS16*crRHS62 + crRHS17*crRHS62 + crRHS22*crRHS63 + crRHS22*crRHS65 - crRHS24*crRHS67 - crRHS32*r_DN(1,0) - crRHS33*crRHS62 - crRHS45*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS25*crRHS62 + crRHS26*crRHS62 + crRHS27*crRHS63 + crRHS27*crRHS65 - crRHS28*crRHS67 - crRHS32*r_DN(1,1) - crRHS45*r_DN(1,1) - crRHS49*crRHS62 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS50*r_N[1] - crRHS51*r_DN(1,0) - crRHS52*r_DN(1,1) + crRHS54*r_N[1] - crRHS55*r_N[1] - crRHS56*r_N[1] + crRHS58*r_N[1] + crRHS59*crRHS64 + crRHS60*r_N[1] - crRHS68*r_N[1]);
rRHS[8]+=-crRHS31*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS24*r_DN(2,0) + crRHS28*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS16*crRHS69 + crRHS17*crRHS69 + crRHS22*crRHS70 + crRHS22*crRHS72 - crRHS24*crRHS73 - crRHS32*r_DN(2,0) - crRHS33*crRHS69 - crRHS45*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS25*crRHS69 + crRHS26*crRHS69 + crRHS27*crRHS70 + crRHS27*crRHS72 - crRHS28*crRHS73 - crRHS32*r_DN(2,1) - crRHS45*r_DN(2,1) - crRHS49*crRHS69 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS50*r_N[2] - crRHS51*r_DN(2,0) - crRHS52*r_DN(2,1) + crRHS54*r_N[2] - crRHS55*r_N[2] - crRHS56*r_N[2] + crRHS58*r_N[2] + crRHS59*crRHS71 + crRHS60*r_N[2] - crRHS68*r_N[2]);
rRHS[12]+=-crRHS31*(-crRHS13*r_N[3] + crRHS2*r_N[3] + crRHS24*r_DN(3,0) + crRHS28*r_DN(3,1) + crRHS8*r_N[3]);
rRHS[13]+=-gauss_weight*(crRHS16*crRHS74 + crRHS17*crRHS74 + crRHS22*crRHS75 + crRHS22*crRHS77 - crRHS24*crRHS78 - crRHS32*r_DN(3,0) - crRHS33*crRHS74 - crRHS45*r_DN(3,0) + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-gauss_weight*(crRHS25*crRHS74 + crRHS26*crRHS74 + crRHS27*crRHS75 + crRHS27*crRHS77 - crRHS28*crRHS78 - crRHS32*r_DN(3,1) - crRHS45*r_DN(3,1) - crRHS49*crRHS74 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=gauss_weight*(crRHS50*r_N[3] - crRHS51*r_DN(3,0) - crRHS52*r_DN(3,1) + crRHS54*r_N[3] - crRHS55*r_N[3] - crRHS56*r_N[3] + crRHS58*r_N[3] + crRHS59*crRHS76 + crRHS60*r_N[3] - crRHS68*r_N[3]);

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