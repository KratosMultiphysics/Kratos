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
const double crLHS2 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crLHS3 = 1.0/crLHS2;
const double crLHS4 = 1.0/c_p;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = 1.0/crLHS5;
const double crLHS7 = crLHS6*gamma*p_th;
const double crLHS8 = crLHS4*crLHS7;
const double crLHS9 = crLHS3*crLHS8;
const double crLHS10 = crLHS9*gauss_weight;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS13 = r_N[0]*r_N[0];
const double crLHS14 = crLHS13*crLHS3;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS19*crLHS7;
const double crLHS21 = crLHS3*crLHS4*tau_u;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS24 = bdf0*crLHS8;
const double crLHS25 = 1.0/(crLHS2*crLHS2);
const double crLHS26 = crLHS25*gauss_weight;
const double crLHS27 = crLHS24*crLHS26;
const double crLHS28 = r_DN(0,0)*r_DN(1,0);
const double crLHS29 = r_DN(0,1)*r_DN(1,1);
const double crLHS30 = crLHS11*(crLHS28 + crLHS29);
const double crLHS31 = r_DN(1,0)*r_N[0];
const double crLHS32 = crLHS3*r_N[0];
const double crLHS33 = crLHS32*r_N[1];
const double crLHS34 = -crLHS12*crLHS33;
const double crLHS35 = bdf0*r_N[1];
const double crLHS36 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS37 = crLHS35 + crLHS36;
const double crLHS38 = crLHS37*crLHS7;
const double crLHS39 = crLHS21*crLHS38;
const double crLHS40 = r_DN(1,1)*r_N[0];
const double crLHS41 = -crLHS23*crLHS33;
const double crLHS42 = crLHS26*crLHS8;
const double crLHS43 = crLHS15*crLHS42;
const double crLHS44 = -crLHS43*r_N[1];
const double crLHS45 = r_DN(0,0)*r_DN(2,0);
const double crLHS46 = r_DN(0,1)*r_DN(2,1);
const double crLHS47 = crLHS11*(crLHS45 + crLHS46);
const double crLHS48 = r_DN(2,0)*r_N[0];
const double crLHS49 = crLHS32*r_N[2];
const double crLHS50 = -crLHS12*crLHS49;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS51 + crLHS52;
const double crLHS54 = crLHS53*crLHS7;
const double crLHS55 = crLHS21*crLHS54;
const double crLHS56 = r_DN(2,1)*r_N[0];
const double crLHS57 = -crLHS23*crLHS49;
const double crLHS58 = -crLHS43*r_N[2];
const double crLHS59 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crLHS60 = crLHS12*crLHS16 + crLHS17*crLHS23;
const double crLHS61 = crLHS60*r_N[0];
const double crLHS62 = crLHS25*tau_u;
const double crLHS63 = crLHS62*crLHS8;
const double crLHS64 = gauss_weight*(crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_N[0]*tau_u - crLHS61*crLHS63 - r_N[0]);
const double crLHS65 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS66 = r_C(0,2)*r_DN(0,0);
const double crLHS67 = crLHS66 + r_C(2,2)*r_DN(0,1);
const double crLHS68 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crLHS69 = crLHS32*crLHS68;
const double crLHS70 = -crLHS69 + r_DN(0,0);
const double crLHS71 = crLHS9*r_DN(0,0);
const double crLHS72 = crLHS71*tau_c;
const double crLHS73 = crLHS32*crLHS8;
const double crLHS74 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS75 = crLHS62*crLHS74;
const double crLHS76 = crLHS18*crLHS75;
const double crLHS77 = crLHS19*crLHS74;
const double crLHS78 = crLHS25*crLHS59*tau_u;
const double crLHS79 = crLHS78*r_N[0];
const double crLHS80 = 1.0/(crLHS2*crLHS2*crLHS2);
const double crLHS81 = crLHS61*crLHS80;
const double crLHS82 = crLHS81*tau_u;
const double crLHS83 = crLHS14*crLHS24 + crLHS18*crLHS73 + crLHS19*crLHS76 + crLHS77*crLHS79 - crLHS77*crLHS82;
const double crLHS84 = crLHS66 + r_C(0,1)*r_DN(0,1);
const double crLHS85 = r_C(1,2)*r_DN(0,1);
const double crLHS86 = crLHS85 + r_C(2,2)*r_DN(0,0);
const double crLHS87 = crLHS71*r_DN(0,1);
const double crLHS88 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crLHS89 = crLHS32*crLHS88;
const double crLHS90 = -crLHS89 + r_DN(0,1);
const double crLHS91 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0));
const double crLHS92 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1));
const double crLHS93 = crLHS91*r_DN(0,0) + crLHS92*r_DN(0,1);
const double crLHS94 = crLHS15 + crLHS93;
const double crLHS95 = crLHS42*tau_c;
const double crLHS96 = crLHS95*r_DN(0,0);
const double crLHS97 = r_DN(0,0)*r_N[1];
const double crLHS98 = crLHS60*crLHS63;
const double crLHS99 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS100 = r_C(0,2)*r_DN(1,0);
const double crLHS101 = crLHS100 + r_C(2,2)*r_DN(1,1);
const double crLHS102 = crLHS3*r_N[1];
const double crLHS103 = crLHS102*crLHS68;
const double crLHS104 = -crLHS103 + r_DN(1,0);
const double crLHS105 = crLHS102*crLHS8;
const double crLHS106 = crLHS105*crLHS15;
const double crLHS107 = crLHS106 + crLHS28*crLHS9;
const double crLHS108 = crLHS37*crLHS74;
const double crLHS109 = crLHS108*crLHS79 - crLHS108*crLHS82 + crLHS36*crLHS73 + crLHS37*crLHS76;
const double crLHS110 = crLHS100 + r_C(0,1)*r_DN(1,1);
const double crLHS111 = r_C(1,2)*r_DN(1,1);
const double crLHS112 = crLHS111 + r_C(2,2)*r_DN(1,0);
const double crLHS113 = crLHS71*r_DN(1,1);
const double crLHS114 = crLHS102*crLHS88;
const double crLHS115 = -crLHS114 + r_DN(1,1);
const double crLHS116 = crLHS91*r_DN(1,0) + crLHS92*r_DN(1,1);
const double crLHS117 = crLHS116 + crLHS35;
const double crLHS118 = r_DN(0,0)*r_N[2];
const double crLHS119 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS120 = r_C(0,2)*r_DN(2,0);
const double crLHS121 = crLHS120 + r_C(2,2)*r_DN(2,1);
const double crLHS122 = crLHS3*r_N[2];
const double crLHS123 = -crLHS122*crLHS68 + r_DN(2,0);
const double crLHS124 = crLHS122*crLHS8;
const double crLHS125 = crLHS124*crLHS15;
const double crLHS126 = crLHS125 + crLHS45*crLHS9;
const double crLHS127 = crLHS53*crLHS74;
const double crLHS128 = crLHS127*crLHS79 - crLHS127*crLHS82 + crLHS52*crLHS73 + crLHS53*crLHS76;
const double crLHS129 = crLHS120 + r_C(0,1)*r_DN(2,1);
const double crLHS130 = r_C(1,2)*r_DN(2,1);
const double crLHS131 = crLHS130 + r_C(2,2)*r_DN(2,0);
const double crLHS132 = crLHS71*r_DN(2,1);
const double crLHS133 = -crLHS122*crLHS88 + r_DN(2,1);
const double crLHS134 = crLHS91*r_DN(2,0) + crLHS92*r_DN(2,1);
const double crLHS135 = crLHS134 + crLHS51;
const double crLHS136 = crLHS85 + r_C(0,1)*r_DN(0,0);
const double crLHS137 = crLHS9*r_DN(0,1);
const double crLHS138 = crLHS137*tau_c;
const double crLHS139 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS140 = crLHS95*r_DN(0,1);
const double crLHS141 = r_DN(0,1)*r_N[1];
const double crLHS142 = crLHS111 + r_C(0,1)*r_DN(1,0);
const double crLHS143 = crLHS137*r_DN(1,0);
const double crLHS144 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS145 = crLHS106 + crLHS29*crLHS9;
const double crLHS146 = r_DN(0,1)*r_N[2];
const double crLHS147 = crLHS130 + r_C(0,1)*r_DN(2,0);
const double crLHS148 = crLHS137*r_DN(2,0);
const double crLHS149 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS150 = crLHS125 + crLHS46*crLHS9;
const double crLHS151 = crLHS7*gauss_weight;
const double crLHS152 = crLHS14*crLHS151;
const double crLHS153 = bdf0*crLHS7;
const double crLHS154 = -crLHS20 + dp_th_dt*r_N[0];
const double crLHS155 = crLHS32*crLHS7;
const double crLHS156 = crLHS7*tau_t;
const double crLHS157 = crLHS156*crLHS81;
const double crLHS158 = crLHS151*r_N[1];
const double crLHS159 = crLHS158*crLHS69;
const double crLHS160 = crLHS158*crLHS89;
const double crLHS161 = -crLHS38 + dp_th_dt*r_N[1];
const double crLHS162 = crLHS102*crLHS7;
const double crLHS163 = crLHS15*crLHS162 + crLHS28*kappa + crLHS29*kappa - crLHS3*dp_th_dt*r_N[0]*r_N[1];
const double crLHS164 = crLHS151*r_N[2];
const double crLHS165 = crLHS164*crLHS69;
const double crLHS166 = crLHS164*crLHS89;
const double crLHS167 = -crLHS54 + dp_th_dt*r_N[2];
const double crLHS168 = crLHS122*crLHS7;
const double crLHS169 = crLHS15*crLHS168 - crLHS3*dp_th_dt*r_N[0]*r_N[2] + crLHS45*kappa + crLHS46*kappa;
const double crLHS170 = r_DN(1,0)*r_DN(1,0);
const double crLHS171 = r_DN(1,1)*r_DN(1,1);
const double crLHS172 = r_N[1]*r_N[1];
const double crLHS173 = crLHS172*crLHS3;
const double crLHS174 = r_DN(1,0)*r_DN(2,0);
const double crLHS175 = r_DN(1,1)*r_DN(2,1);
const double crLHS176 = crLHS11*(crLHS174 + crLHS175);
const double crLHS177 = r_DN(2,0)*r_N[1];
const double crLHS178 = crLHS102*r_N[2];
const double crLHS179 = -crLHS12*crLHS178;
const double crLHS180 = r_DN(2,1)*r_N[1];
const double crLHS181 = -crLHS178*crLHS23;
const double crLHS182 = -crLHS35*crLHS42*r_N[2];
const double crLHS183 = crLHS9*r_DN(1,0);
const double crLHS184 = crLHS183*tau_c;
const double crLHS185 = crLHS36*crLHS75;
const double crLHS186 = crLHS78*r_N[1];
const double crLHS187 = crLHS60*r_N[1];
const double crLHS188 = crLHS80*tau_u;
const double crLHS189 = crLHS187*crLHS188;
const double crLHS190 = crLHS105*crLHS18 + crLHS185*crLHS19 + crLHS186*crLHS77 - crLHS189*crLHS77;
const double crLHS191 = crLHS95*r_DN(1,0);
const double crLHS192 = gauss_weight*(-crLHS187*crLHS63 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS193 = crLHS105*crLHS36 + crLHS108*crLHS186 - crLHS108*crLHS189 + crLHS173*crLHS24 + crLHS185*crLHS37;
const double crLHS194 = crLHS183*r_DN(1,1);
const double crLHS195 = r_DN(1,0)*r_N[2];
const double crLHS196 = crLHS124*crLHS35;
const double crLHS197 = crLHS174*crLHS9 + crLHS196;
const double crLHS198 = crLHS105*crLHS52 + crLHS127*crLHS186 - crLHS127*crLHS189 + crLHS185*crLHS53;
const double crLHS199 = crLHS183*r_DN(2,1);
const double crLHS200 = crLHS70*tau_c;
const double crLHS201 = crLHS9*r_DN(1,1);
const double crLHS202 = crLHS201*tau_c;
const double crLHS203 = crLHS95*r_DN(1,1);
const double crLHS204 = r_DN(1,1)*r_N[2];
const double crLHS205 = crLHS9*r_DN(2,0);
const double crLHS206 = crLHS205*r_DN(1,1);
const double crLHS207 = crLHS175*crLHS9 + crLHS196;
const double crLHS208 = crLHS156*crLHS80;
const double crLHS209 = crLHS187*crLHS208;
const double crLHS210 = crLHS151*crLHS173;
const double crLHS211 = crLHS103*crLHS164;
const double crLHS212 = crLHS114*crLHS164;
const double crLHS213 = crLHS168*crLHS35 + crLHS174*kappa + crLHS175*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[2];
const double crLHS214 = r_DN(2,0)*r_DN(2,0);
const double crLHS215 = r_DN(2,1)*r_DN(2,1);
const double crLHS216 = r_N[2]*r_N[2];
const double crLHS217 = crLHS216*crLHS3;
const double crLHS218 = crLHS52*crLHS75;
const double crLHS219 = crLHS78*r_N[2];
const double crLHS220 = crLHS60*r_N[2];
const double crLHS221 = crLHS188*crLHS220;
const double crLHS222 = crLHS124*crLHS18 + crLHS19*crLHS218 + crLHS219*crLHS77 - crLHS221*crLHS77;
const double crLHS223 = crLHS205*tau_c;
const double crLHS224 = crLHS95*r_DN(2,0);
const double crLHS225 = crLHS108*crLHS219 - crLHS108*crLHS221 + crLHS124*crLHS36 + crLHS218*crLHS37;
const double crLHS226 = gauss_weight*(-crLHS220*crLHS63 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS227 = crLHS124*crLHS52 + crLHS127*crLHS219 - crLHS127*crLHS221 + crLHS217*crLHS24 + crLHS218*crLHS53;
const double crLHS228 = crLHS205*r_DN(2,1);
const double crLHS229 = crLHS9*r_DN(2,1);
const double crLHS230 = crLHS229*tau_c;
const double crLHS231 = crLHS95*r_DN(2,1);
const double crLHS232 = crLHS208*crLHS220;
const double crLHS233 = crLHS151*crLHS217;
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
rLHS(0,11)+=crLHS58;
rLHS(1,0)+=crLHS64*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS0*crLHS9 + crLHS65*r_DN(0,0) + crLHS67*r_DN(0,1) + crLHS70*crLHS72 + crLHS83);
rLHS(1,2)+=gauss_weight*(crLHS72*crLHS90 + crLHS84*r_DN(0,0) + crLHS86*r_DN(0,1) + crLHS87);
rLHS(1,3)+=-crLHS94*crLHS96;
rLHS(1,4)+=gauss_weight*(crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS31*crLHS98 - crLHS97);
rLHS(1,5)+=gauss_weight*(crLHS101*r_DN(0,1) + crLHS104*crLHS72 + crLHS107 + crLHS109 + crLHS99*r_DN(0,0));
rLHS(1,6)+=gauss_weight*(crLHS110*r_DN(0,0) + crLHS112*r_DN(0,1) + crLHS113 + crLHS115*crLHS72);
rLHS(1,7)+=-crLHS117*crLHS96;
rLHS(1,8)+=gauss_weight*(-crLHS118 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS48*crLHS98);
rLHS(1,9)+=gauss_weight*(crLHS119*r_DN(0,0) + crLHS121*r_DN(0,1) + crLHS123*crLHS72 + crLHS126 + crLHS128);
rLHS(1,10)+=gauss_weight*(crLHS129*r_DN(0,0) + crLHS131*r_DN(0,1) + crLHS132 + crLHS133*crLHS72);
rLHS(1,11)+=-crLHS135*crLHS96;
rLHS(2,0)+=crLHS64*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS136*r_DN(0,1) + crLHS138*crLHS70 + crLHS67*r_DN(0,0) + crLHS87);
rLHS(2,2)+=gauss_weight*(crLHS1*crLHS9 + crLHS138*crLHS90 + crLHS139*r_DN(0,1) + crLHS83 + crLHS86*r_DN(0,0));
rLHS(2,3)+=-crLHS140*crLHS94;
rLHS(2,4)+=gauss_weight*(-crLHS141 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS40*crLHS98);
rLHS(2,5)+=gauss_weight*(crLHS101*r_DN(0,0) + crLHS104*crLHS138 + crLHS142*r_DN(0,1) + crLHS143);
rLHS(2,6)+=gauss_weight*(crLHS109 + crLHS112*r_DN(0,0) + crLHS115*crLHS138 + crLHS144*r_DN(0,1) + crLHS145);
rLHS(2,7)+=-crLHS117*crLHS140;
rLHS(2,8)+=gauss_weight*(-crLHS146 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS56*crLHS98);
rLHS(2,9)+=gauss_weight*(crLHS121*r_DN(0,0) + crLHS123*crLHS138 + crLHS147*r_DN(0,1) + crLHS148);
rLHS(2,10)+=gauss_weight*(crLHS128 + crLHS131*r_DN(0,0) + crLHS133*crLHS138 + crLHS149*r_DN(0,1) + crLHS150);
rLHS(2,11)+=-crLHS135*crLHS140;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS152*crLHS68;
rLHS(3,2)+=crLHS152*crLHS88;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa + crLHS13*crLHS3*dp_th_dt - crLHS14*crLHS153 - crLHS154*crLHS157 + crLHS154*crLHS18*crLHS25*crLHS6*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[0]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS155*crLHS93);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS159;
rLHS(3,6)+=crLHS160;
rLHS(3,7)+=-gauss_weight*(-crLHS116*crLHS155 - crLHS157*crLHS161 + crLHS161*crLHS18*crLHS25*crLHS6*gamma*p_th*tau_t + crLHS161*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[0]*tau_t + crLHS161*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS163);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS165;
rLHS(3,10)+=crLHS166;
rLHS(3,11)+=-gauss_weight*(-crLHS134*crLHS155 - crLHS157*crLHS167 + crLHS167*crLHS18*crLHS25*crLHS6*gamma*p_th*tau_t + crLHS167*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[0]*tau_t + crLHS167*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS169);
rLHS(4,0)+=crLHS30;
rLHS(4,1)+=crLHS10*(crLHS22*r_DN(1,0) + crLHS34 + crLHS97);
rLHS(4,2)+=crLHS10*(crLHS141 + crLHS22*r_DN(1,1) + crLHS41);
rLHS(4,3)+=crLHS44;
rLHS(4,4)+=crLHS11*(crLHS170 + crLHS171);
rLHS(4,5)+=crLHS10*(-crLHS12*crLHS173 + crLHS39*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS173*crLHS23 + crLHS39*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS172*crLHS27;
rLHS(4,8)+=crLHS176;
rLHS(4,9)+=crLHS10*(crLHS177 + crLHS179 + crLHS55*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS180 + crLHS181 + crLHS55*r_DN(1,1));
rLHS(4,11)+=crLHS182;
rLHS(5,0)+=gauss_weight*(crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS31 - crLHS97*crLHS98);
rLHS(5,1)+=gauss_weight*(crLHS107 + crLHS184*crLHS70 + crLHS190 + crLHS65*r_DN(1,0) + crLHS67*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS143 + crLHS184*crLHS90 + crLHS84*r_DN(1,0) + crLHS86*r_DN(1,1));
rLHS(5,3)+=-crLHS191*crLHS94;
rLHS(5,4)+=crLHS192*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS101*r_DN(1,1) + crLHS104*crLHS184 + crLHS170*crLHS9 + crLHS193 + crLHS99*r_DN(1,0));
rLHS(5,6)+=gauss_weight*(crLHS110*r_DN(1,0) + crLHS112*r_DN(1,1) + crLHS115*crLHS184 + crLHS194);
rLHS(5,7)+=-crLHS117*crLHS191;
rLHS(5,8)+=gauss_weight*(-crLHS177*crLHS98 - crLHS195 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS119*r_DN(1,0) + crLHS121*r_DN(1,1) + crLHS123*crLHS184 + crLHS197 + crLHS198);
rLHS(5,10)+=gauss_weight*(crLHS129*r_DN(1,0) + crLHS131*r_DN(1,1) + crLHS133*crLHS184 + crLHS199);
rLHS(5,11)+=-crLHS135*crLHS191;
rLHS(6,0)+=gauss_weight*(-crLHS141*crLHS98 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS40);
rLHS(6,1)+=gauss_weight*(crLHS113 + crLHS136*r_DN(1,1) + crLHS200*crLHS201 + crLHS67*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS139*r_DN(1,1) + crLHS145 + crLHS190 + crLHS202*crLHS90 + crLHS86*r_DN(1,0));
rLHS(6,3)+=-crLHS203*crLHS94;
rLHS(6,4)+=crLHS192*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS101*r_DN(1,0) + crLHS104*crLHS202 + crLHS142*r_DN(1,1) + crLHS194);
rLHS(6,6)+=gauss_weight*(crLHS112*r_DN(1,0) + crLHS115*crLHS202 + crLHS144*r_DN(1,1) + crLHS171*crLHS9 + crLHS193);
rLHS(6,7)+=-crLHS117*crLHS203;
rLHS(6,8)+=gauss_weight*(-crLHS180*crLHS98 - crLHS204 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS121*r_DN(1,0) + crLHS123*crLHS202 + crLHS147*r_DN(1,1) + crLHS206);
rLHS(6,10)+=gauss_weight*(crLHS131*r_DN(1,0) + crLHS133*crLHS202 + crLHS149*r_DN(1,1) + crLHS198 + crLHS207);
rLHS(6,11)+=-crLHS135*crLHS203;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS159;
rLHS(7,2)+=crLHS160;
rLHS(7,3)+=-gauss_weight*(-crLHS154*crLHS209 + crLHS154*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[1]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS162*crLHS93 - crLHS163);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS210*crLHS68;
rLHS(7,6)+=crLHS210*crLHS88;
rLHS(7,7)+=-gauss_weight*(-crLHS116*crLHS162 - crLHS153*crLHS173 - crLHS161*crLHS209 + crLHS161*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS161*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[1]*tau_t + crLHS161*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS170*kappa - crLHS171*kappa + crLHS172*crLHS3*dp_th_dt);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS211;
rLHS(7,10)+=crLHS212;
rLHS(7,11)+=-gauss_weight*(-crLHS134*crLHS162 - crLHS167*crLHS209 + crLHS167*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS167*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[1]*tau_t + crLHS167*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS213);
rLHS(8,0)+=crLHS47;
rLHS(8,1)+=crLHS10*(crLHS118 + crLHS22*r_DN(2,0) + crLHS50);
rLHS(8,2)+=crLHS10*(crLHS146 + crLHS22*r_DN(2,1) + crLHS57);
rLHS(8,3)+=crLHS58;
rLHS(8,4)+=crLHS176;
rLHS(8,5)+=crLHS10*(crLHS179 + crLHS195 + crLHS39*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS181 + crLHS204 + crLHS39*r_DN(2,1));
rLHS(8,7)+=crLHS182;
rLHS(8,8)+=crLHS11*(crLHS214 + crLHS215);
rLHS(8,9)+=crLHS10*(-crLHS12*crLHS217 + crLHS55*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS217*crLHS23 + crLHS55*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS216*crLHS27;
rLHS(9,0)+=gauss_weight*(-crLHS118*crLHS98 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS48);
rLHS(9,1)+=gauss_weight*(crLHS126 + crLHS200*crLHS205 + crLHS222 + crLHS65*r_DN(2,0) + crLHS67*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS148 + crLHS223*crLHS90 + crLHS84*r_DN(2,0) + crLHS86*r_DN(2,1));
rLHS(9,3)+=-crLHS224*crLHS94;
rLHS(9,4)+=gauss_weight*(-crLHS177 - crLHS195*crLHS98 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS101*r_DN(2,1) + crLHS104*crLHS223 + crLHS197 + crLHS225 + crLHS99*r_DN(2,0));
rLHS(9,6)+=gauss_weight*(crLHS110*r_DN(2,0) + crLHS112*r_DN(2,1) + crLHS115*crLHS223 + crLHS206);
rLHS(9,7)+=-crLHS117*crLHS224;
rLHS(9,8)+=crLHS226*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS119*r_DN(2,0) + crLHS121*r_DN(2,1) + crLHS123*crLHS223 + crLHS214*crLHS9 + crLHS227);
rLHS(9,10)+=gauss_weight*(crLHS129*r_DN(2,0) + crLHS131*r_DN(2,1) + crLHS133*crLHS223 + crLHS228);
rLHS(9,11)+=-crLHS135*crLHS224;
rLHS(10,0)+=gauss_weight*(-crLHS146*crLHS98 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS56);
rLHS(10,1)+=gauss_weight*(crLHS132 + crLHS136*r_DN(2,1) + crLHS200*crLHS229 + crLHS67*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS139*r_DN(2,1) + crLHS150 + crLHS222 + crLHS230*crLHS90 + crLHS86*r_DN(2,0));
rLHS(10,3)+=-crLHS231*crLHS94;
rLHS(10,4)+=gauss_weight*(-crLHS180 - crLHS204*crLHS98 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS59*crLHS6*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS101*r_DN(2,0) + crLHS104*crLHS230 + crLHS142*r_DN(2,1) + crLHS199);
rLHS(10,6)+=gauss_weight*(crLHS112*r_DN(2,0) + crLHS115*crLHS230 + crLHS144*r_DN(2,1) + crLHS207 + crLHS225);
rLHS(10,7)+=-crLHS117*crLHS231;
rLHS(10,8)+=crLHS226*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS121*r_DN(2,0) + crLHS123*crLHS230 + crLHS147*r_DN(2,1) + crLHS228);
rLHS(10,10)+=gauss_weight*(crLHS131*r_DN(2,0) + crLHS133*crLHS230 + crLHS149*r_DN(2,1) + crLHS215*crLHS9 + crLHS227);
rLHS(10,11)+=-crLHS135*crLHS231;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS165;
rLHS(11,2)+=crLHS166;
rLHS(11,3)+=-gauss_weight*(-crLHS154*crLHS232 + crLHS154*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS154*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[2]*tau_t + crLHS154*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS168*crLHS93 - crLHS169);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS211;
rLHS(11,6)+=crLHS212;
rLHS(11,7)+=-gauss_weight*(-crLHS116*crLHS168 - crLHS161*crLHS232 + crLHS161*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS161*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[2]*tau_t + crLHS161*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS213);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS233*crLHS68;
rLHS(11,10)+=crLHS233*crLHS88;
rLHS(11,11)+=-gauss_weight*(-crLHS134*crLHS168 - crLHS153*crLHS217 - crLHS167*crLHS232 + crLHS167*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS167*crLHS25*crLHS59*crLHS6*gamma*p_th*r_N[2]*tau_t + crLHS167*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS214*kappa - crLHS215*kappa + crLHS216*crLHS3*dp_th_dt);

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
const double crLHS2 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crLHS3 = 1.0/crLHS2;
const double crLHS4 = 1.0/c_p;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = 1.0/crLHS5;
const double crLHS7 = crLHS6*gamma*p_th;
const double crLHS8 = crLHS4*crLHS7;
const double crLHS9 = crLHS3*crLHS8;
const double crLHS10 = crLHS9*gauss_weight;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS13 = r_N[0]*r_N[0];
const double crLHS14 = crLHS13*crLHS3;
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = lin_u_conv(0,0)*r_N[0] + lin_u_conv(1,0)*r_N[1] + lin_u_conv(2,0)*r_N[2] + lin_u_conv(3,0)*r_N[3];
const double crLHS17 = lin_u_conv(0,1)*r_N[0] + lin_u_conv(1,1)*r_N[1] + lin_u_conv(2,1)*r_N[2] + lin_u_conv(3,1)*r_N[3];
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS19*crLHS7;
const double crLHS21 = crLHS3*crLHS4*tau_u;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS24 = bdf0*crLHS8;
const double crLHS25 = 1.0/(crLHS2*crLHS2);
const double crLHS26 = crLHS25*gauss_weight;
const double crLHS27 = crLHS24*crLHS26;
const double crLHS28 = r_DN(0,0)*r_DN(1,0);
const double crLHS29 = r_DN(0,1)*r_DN(1,1);
const double crLHS30 = crLHS11*(crLHS28 + crLHS29);
const double crLHS31 = r_DN(1,0)*r_N[0];
const double crLHS32 = crLHS3*r_N[0];
const double crLHS33 = crLHS32*r_N[1];
const double crLHS34 = -crLHS12*crLHS33;
const double crLHS35 = bdf0*r_N[1];
const double crLHS36 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS37 = crLHS35 + crLHS36;
const double crLHS38 = crLHS37*crLHS7;
const double crLHS39 = crLHS21*crLHS38;
const double crLHS40 = r_DN(1,1)*r_N[0];
const double crLHS41 = -crLHS23*crLHS33;
const double crLHS42 = crLHS26*crLHS8;
const double crLHS43 = crLHS15*crLHS42;
const double crLHS44 = -crLHS43*r_N[1];
const double crLHS45 = r_DN(0,0)*r_DN(2,0);
const double crLHS46 = r_DN(0,1)*r_DN(2,1);
const double crLHS47 = crLHS11*(crLHS45 + crLHS46);
const double crLHS48 = r_DN(2,0)*r_N[0];
const double crLHS49 = crLHS32*r_N[2];
const double crLHS50 = -crLHS12*crLHS49;
const double crLHS51 = bdf0*r_N[2];
const double crLHS52 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS53 = crLHS51 + crLHS52;
const double crLHS54 = crLHS53*crLHS7;
const double crLHS55 = crLHS21*crLHS54;
const double crLHS56 = r_DN(2,1)*r_N[0];
const double crLHS57 = -crLHS23*crLHS49;
const double crLHS58 = -crLHS43*r_N[2];
const double crLHS59 = r_DN(0,0)*r_DN(3,0);
const double crLHS60 = r_DN(0,1)*r_DN(3,1);
const double crLHS61 = crLHS11*(crLHS59 + crLHS60);
const double crLHS62 = r_DN(3,0)*r_N[0];
const double crLHS63 = crLHS32*r_N[3];
const double crLHS64 = -crLHS12*crLHS63;
const double crLHS65 = bdf0*r_N[3];
const double crLHS66 = crLHS16*r_DN(3,0) + crLHS17*r_DN(3,1);
const double crLHS67 = crLHS65 + crLHS66;
const double crLHS68 = crLHS67*crLHS7;
const double crLHS69 = crLHS21*crLHS68;
const double crLHS70 = r_DN(3,1)*r_N[0];
const double crLHS71 = -crLHS23*crLHS63;
const double crLHS72 = -crLHS43*r_N[3];
const double crLHS73 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crLHS74 = crLHS12*crLHS16 + crLHS17*crLHS23;
const double crLHS75 = crLHS74*r_N[0];
const double crLHS76 = crLHS25*tau_u;
const double crLHS77 = crLHS76*crLHS8;
const double crLHS78 = gauss_weight*(crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_u - crLHS75*crLHS77 - r_N[0]);
const double crLHS79 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS80 = r_C(0,2)*r_DN(0,0);
const double crLHS81 = crLHS80 + r_C(2,2)*r_DN(0,1);
const double crLHS82 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crLHS83 = crLHS32*crLHS82;
const double crLHS84 = -crLHS83 + r_DN(0,0);
const double crLHS85 = crLHS9*r_DN(0,0);
const double crLHS86 = crLHS85*tau_c;
const double crLHS87 = crLHS32*crLHS8;
const double crLHS88 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS89 = crLHS76*crLHS88;
const double crLHS90 = crLHS18*crLHS89;
const double crLHS91 = crLHS19*crLHS88;
const double crLHS92 = crLHS25*crLHS73*tau_u;
const double crLHS93 = crLHS92*r_N[0];
const double crLHS94 = 1.0/(crLHS2*crLHS2*crLHS2);
const double crLHS95 = crLHS75*crLHS94;
const double crLHS96 = crLHS95*tau_u;
const double crLHS97 = crLHS14*crLHS24 + crLHS18*crLHS87 + crLHS19*crLHS90 + crLHS91*crLHS93 - crLHS91*crLHS96;
const double crLHS98 = crLHS80 + r_C(0,1)*r_DN(0,1);
const double crLHS99 = r_C(1,2)*r_DN(0,1);
const double crLHS100 = crLHS99 + r_C(2,2)*r_DN(0,0);
const double crLHS101 = crLHS85*r_DN(0,1);
const double crLHS102 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crLHS103 = crLHS102*crLHS32;
const double crLHS104 = -crLHS103 + r_DN(0,1);
const double crLHS105 = r_N[0]*(r_u(0,0) - r_u_mesh(0,0)) + r_N[1]*(r_u(1,0) - r_u_mesh(1,0)) + r_N[2]*(r_u(2,0) - r_u_mesh(2,0)) + r_N[3]*(r_u(3,0) - r_u_mesh(3,0));
const double crLHS106 = r_N[0]*(r_u(0,1) - r_u_mesh(0,1)) + r_N[1]*(r_u(1,1) - r_u_mesh(1,1)) + r_N[2]*(r_u(2,1) - r_u_mesh(2,1)) + r_N[3]*(r_u(3,1) - r_u_mesh(3,1));
const double crLHS107 = crLHS105*r_DN(0,0) + crLHS106*r_DN(0,1);
const double crLHS108 = crLHS107 + crLHS15;
const double crLHS109 = crLHS42*tau_c;
const double crLHS110 = crLHS109*r_DN(0,0);
const double crLHS111 = r_DN(0,0)*r_N[1];
const double crLHS112 = crLHS74*crLHS77;
const double crLHS113 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS114 = r_C(0,2)*r_DN(1,0);
const double crLHS115 = crLHS114 + r_C(2,2)*r_DN(1,1);
const double crLHS116 = crLHS3*r_N[1];
const double crLHS117 = crLHS116*crLHS82;
const double crLHS118 = -crLHS117 + r_DN(1,0);
const double crLHS119 = crLHS116*crLHS8;
const double crLHS120 = crLHS119*crLHS15;
const double crLHS121 = crLHS120 + crLHS28*crLHS9;
const double crLHS122 = crLHS37*crLHS88;
const double crLHS123 = crLHS122*crLHS93 - crLHS122*crLHS96 + crLHS36*crLHS87 + crLHS37*crLHS90;
const double crLHS124 = crLHS114 + r_C(0,1)*r_DN(1,1);
const double crLHS125 = r_C(1,2)*r_DN(1,1);
const double crLHS126 = crLHS125 + r_C(2,2)*r_DN(1,0);
const double crLHS127 = crLHS85*r_DN(1,1);
const double crLHS128 = crLHS102*crLHS116;
const double crLHS129 = -crLHS128 + r_DN(1,1);
const double crLHS130 = crLHS105*r_DN(1,0) + crLHS106*r_DN(1,1);
const double crLHS131 = crLHS130 + crLHS35;
const double crLHS132 = r_DN(0,0)*r_N[2];
const double crLHS133 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS134 = r_C(0,2)*r_DN(2,0);
const double crLHS135 = crLHS134 + r_C(2,2)*r_DN(2,1);
const double crLHS136 = crLHS3*r_N[2];
const double crLHS137 = crLHS136*crLHS82;
const double crLHS138 = -crLHS137 + r_DN(2,0);
const double crLHS139 = crLHS136*crLHS8;
const double crLHS140 = crLHS139*crLHS15;
const double crLHS141 = crLHS140 + crLHS45*crLHS9;
const double crLHS142 = crLHS53*crLHS88;
const double crLHS143 = crLHS142*crLHS93 - crLHS142*crLHS96 + crLHS52*crLHS87 + crLHS53*crLHS90;
const double crLHS144 = crLHS134 + r_C(0,1)*r_DN(2,1);
const double crLHS145 = r_C(1,2)*r_DN(2,1);
const double crLHS146 = crLHS145 + r_C(2,2)*r_DN(2,0);
const double crLHS147 = crLHS85*r_DN(2,1);
const double crLHS148 = crLHS102*crLHS136;
const double crLHS149 = -crLHS148 + r_DN(2,1);
const double crLHS150 = crLHS105*r_DN(2,0) + crLHS106*r_DN(2,1);
const double crLHS151 = crLHS150 + crLHS51;
const double crLHS152 = r_DN(0,0)*r_N[3];
const double crLHS153 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS154 = r_C(0,2)*r_DN(3,0);
const double crLHS155 = crLHS154 + r_C(2,2)*r_DN(3,1);
const double crLHS156 = crLHS3*r_N[3];
const double crLHS157 = -crLHS156*crLHS82 + r_DN(3,0);
const double crLHS158 = crLHS156*crLHS8;
const double crLHS159 = crLHS15*crLHS158;
const double crLHS160 = crLHS159 + crLHS59*crLHS9;
const double crLHS161 = crLHS67*crLHS88;
const double crLHS162 = crLHS161*crLHS93 - crLHS161*crLHS96 + crLHS66*crLHS87 + crLHS67*crLHS90;
const double crLHS163 = crLHS154 + r_C(0,1)*r_DN(3,1);
const double crLHS164 = r_C(1,2)*r_DN(3,1);
const double crLHS165 = crLHS164 + r_C(2,2)*r_DN(3,0);
const double crLHS166 = crLHS85*r_DN(3,1);
const double crLHS167 = -crLHS102*crLHS156 + r_DN(3,1);
const double crLHS168 = crLHS105*r_DN(3,0) + crLHS106*r_DN(3,1);
const double crLHS169 = crLHS168 + crLHS65;
const double crLHS170 = crLHS99 + r_C(0,1)*r_DN(0,0);
const double crLHS171 = crLHS9*r_DN(0,1);
const double crLHS172 = crLHS171*tau_c;
const double crLHS173 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS174 = crLHS109*r_DN(0,1);
const double crLHS175 = r_DN(0,1)*r_N[1];
const double crLHS176 = crLHS125 + r_C(0,1)*r_DN(1,0);
const double crLHS177 = crLHS171*r_DN(1,0);
const double crLHS178 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS179 = crLHS120 + crLHS29*crLHS9;
const double crLHS180 = r_DN(0,1)*r_N[2];
const double crLHS181 = crLHS145 + r_C(0,1)*r_DN(2,0);
const double crLHS182 = crLHS171*r_DN(2,0);
const double crLHS183 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS184 = crLHS140 + crLHS46*crLHS9;
const double crLHS185 = r_DN(0,1)*r_N[3];
const double crLHS186 = crLHS164 + r_C(0,1)*r_DN(3,0);
const double crLHS187 = crLHS171*r_DN(3,0);
const double crLHS188 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS189 = crLHS159 + crLHS60*crLHS9;
const double crLHS190 = crLHS7*gauss_weight;
const double crLHS191 = crLHS14*crLHS190;
const double crLHS192 = bdf0*crLHS7;
const double crLHS193 = -crLHS20 + dp_th_dt*r_N[0];
const double crLHS194 = crLHS32*crLHS7;
const double crLHS195 = crLHS7*tau_t;
const double crLHS196 = crLHS195*crLHS95;
const double crLHS197 = crLHS190*r_N[1];
const double crLHS198 = crLHS197*crLHS83;
const double crLHS199 = crLHS103*crLHS197;
const double crLHS200 = -crLHS38 + dp_th_dt*r_N[1];
const double crLHS201 = crLHS116*crLHS7;
const double crLHS202 = crLHS15*crLHS201 + crLHS28*kappa + crLHS29*kappa - crLHS3*dp_th_dt*r_N[0]*r_N[1];
const double crLHS203 = crLHS190*r_N[2];
const double crLHS204 = crLHS203*crLHS83;
const double crLHS205 = crLHS103*crLHS203;
const double crLHS206 = -crLHS54 + dp_th_dt*r_N[2];
const double crLHS207 = crLHS136*crLHS7;
const double crLHS208 = crLHS15*crLHS207 - crLHS3*dp_th_dt*r_N[0]*r_N[2] + crLHS45*kappa + crLHS46*kappa;
const double crLHS209 = crLHS190*r_N[3];
const double crLHS210 = crLHS209*crLHS83;
const double crLHS211 = crLHS103*crLHS209;
const double crLHS212 = -crLHS68 + dp_th_dt*r_N[3];
const double crLHS213 = crLHS156*crLHS7;
const double crLHS214 = crLHS15*crLHS213 - crLHS3*dp_th_dt*r_N[0]*r_N[3] + crLHS59*kappa + crLHS60*kappa;
const double crLHS215 = r_DN(1,0)*r_DN(1,0);
const double crLHS216 = r_DN(1,1)*r_DN(1,1);
const double crLHS217 = r_N[1]*r_N[1];
const double crLHS218 = crLHS217*crLHS3;
const double crLHS219 = r_DN(1,0)*r_DN(2,0);
const double crLHS220 = r_DN(1,1)*r_DN(2,1);
const double crLHS221 = crLHS11*(crLHS219 + crLHS220);
const double crLHS222 = r_DN(2,0)*r_N[1];
const double crLHS223 = crLHS116*r_N[2];
const double crLHS224 = -crLHS12*crLHS223;
const double crLHS225 = r_DN(2,1)*r_N[1];
const double crLHS226 = -crLHS223*crLHS23;
const double crLHS227 = crLHS35*crLHS42;
const double crLHS228 = -crLHS227*r_N[2];
const double crLHS229 = r_DN(1,0)*r_DN(3,0);
const double crLHS230 = r_DN(1,1)*r_DN(3,1);
const double crLHS231 = crLHS11*(crLHS229 + crLHS230);
const double crLHS232 = r_DN(3,0)*r_N[1];
const double crLHS233 = crLHS116*r_N[3];
const double crLHS234 = -crLHS12*crLHS233;
const double crLHS235 = r_DN(3,1)*r_N[1];
const double crLHS236 = -crLHS23*crLHS233;
const double crLHS237 = -crLHS227*r_N[3];
const double crLHS238 = crLHS9*r_DN(1,0);
const double crLHS239 = crLHS238*tau_c;
const double crLHS240 = crLHS36*crLHS89;
const double crLHS241 = crLHS92*r_N[1];
const double crLHS242 = crLHS74*r_N[1];
const double crLHS243 = crLHS94*tau_u;
const double crLHS244 = crLHS242*crLHS243;
const double crLHS245 = crLHS119*crLHS18 + crLHS19*crLHS240 + crLHS241*crLHS91 - crLHS244*crLHS91;
const double crLHS246 = crLHS109*r_DN(1,0);
const double crLHS247 = gauss_weight*(-crLHS242*crLHS77 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS248 = crLHS119*crLHS36 + crLHS122*crLHS241 - crLHS122*crLHS244 + crLHS218*crLHS24 + crLHS240*crLHS37;
const double crLHS249 = crLHS238*r_DN(1,1);
const double crLHS250 = r_DN(1,0)*r_N[2];
const double crLHS251 = crLHS139*crLHS35;
const double crLHS252 = crLHS219*crLHS9 + crLHS251;
const double crLHS253 = crLHS119*crLHS52 + crLHS142*crLHS241 - crLHS142*crLHS244 + crLHS240*crLHS53;
const double crLHS254 = crLHS238*r_DN(2,1);
const double crLHS255 = r_DN(1,0)*r_N[3];
const double crLHS256 = crLHS158*crLHS35;
const double crLHS257 = crLHS229*crLHS9 + crLHS256;
const double crLHS258 = crLHS119*crLHS66 + crLHS161*crLHS241 - crLHS161*crLHS244 + crLHS240*crLHS67;
const double crLHS259 = crLHS238*r_DN(3,1);
const double crLHS260 = crLHS9*r_DN(1,1);
const double crLHS261 = crLHS260*tau_c;
const double crLHS262 = crLHS109*r_DN(1,1);
const double crLHS263 = r_DN(1,1)*r_N[2];
const double crLHS264 = crLHS260*r_DN(2,0);
const double crLHS265 = crLHS220*crLHS9 + crLHS251;
const double crLHS266 = r_DN(1,1)*r_N[3];
const double crLHS267 = crLHS260*r_DN(3,0);
const double crLHS268 = crLHS230*crLHS9 + crLHS256;
const double crLHS269 = crLHS195*crLHS94;
const double crLHS270 = crLHS242*crLHS269;
const double crLHS271 = crLHS190*crLHS218;
const double crLHS272 = crLHS117*crLHS203;
const double crLHS273 = crLHS128*crLHS203;
const double crLHS274 = crLHS207*crLHS35 + crLHS219*kappa + crLHS220*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[2];
const double crLHS275 = crLHS117*crLHS209;
const double crLHS276 = crLHS128*crLHS209;
const double crLHS277 = crLHS213*crLHS35 + crLHS229*kappa + crLHS230*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[3];
const double crLHS278 = r_DN(2,0)*r_DN(2,0);
const double crLHS279 = r_DN(2,1)*r_DN(2,1);
const double crLHS280 = r_N[2]*r_N[2];
const double crLHS281 = crLHS280*crLHS3;
const double crLHS282 = r_DN(2,0)*r_DN(3,0);
const double crLHS283 = r_DN(2,1)*r_DN(3,1);
const double crLHS284 = crLHS11*(crLHS282 + crLHS283);
const double crLHS285 = r_DN(3,0)*r_N[2];
const double crLHS286 = crLHS136*r_N[3];
const double crLHS287 = -crLHS12*crLHS286;
const double crLHS288 = r_DN(3,1)*r_N[2];
const double crLHS289 = -crLHS23*crLHS286;
const double crLHS290 = -crLHS42*crLHS51*r_N[3];
const double crLHS291 = crLHS9*r_DN(2,0);
const double crLHS292 = crLHS291*tau_c;
const double crLHS293 = crLHS52*crLHS89;
const double crLHS294 = crLHS92*r_N[2];
const double crLHS295 = crLHS74*r_N[2];
const double crLHS296 = crLHS243*crLHS295;
const double crLHS297 = crLHS139*crLHS18 + crLHS19*crLHS293 + crLHS294*crLHS91 - crLHS296*crLHS91;
const double crLHS298 = crLHS109*r_DN(2,0);
const double crLHS299 = crLHS122*crLHS294 - crLHS122*crLHS296 + crLHS139*crLHS36 + crLHS293*crLHS37;
const double crLHS300 = gauss_weight*(-crLHS295*crLHS77 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS301 = crLHS139*crLHS52 + crLHS142*crLHS294 - crLHS142*crLHS296 + crLHS24*crLHS281 + crLHS293*crLHS53;
const double crLHS302 = crLHS291*r_DN(2,1);
const double crLHS303 = r_DN(2,0)*r_N[3];
const double crLHS304 = crLHS158*crLHS51;
const double crLHS305 = crLHS282*crLHS9 + crLHS304;
const double crLHS306 = crLHS139*crLHS66 + crLHS161*crLHS294 - crLHS161*crLHS296 + crLHS293*crLHS67;
const double crLHS307 = crLHS291*r_DN(3,1);
const double crLHS308 = crLHS84*tau_c;
const double crLHS309 = crLHS9*r_DN(2,1);
const double crLHS310 = crLHS309*tau_c;
const double crLHS311 = crLHS109*r_DN(2,1);
const double crLHS312 = r_DN(2,1)*r_N[3];
const double crLHS313 = crLHS9*r_DN(3,0);
const double crLHS314 = crLHS313*r_DN(2,1);
const double crLHS315 = crLHS283*crLHS9 + crLHS304;
const double crLHS316 = crLHS269*crLHS295;
const double crLHS317 = crLHS190*crLHS281;
const double crLHS318 = crLHS137*crLHS209;
const double crLHS319 = crLHS148*crLHS209;
const double crLHS320 = crLHS213*crLHS51 + crLHS282*kappa + crLHS283*kappa - crLHS3*dp_th_dt*r_N[2]*r_N[3];
const double crLHS321 = r_DN(3,0)*r_DN(3,0);
const double crLHS322 = r_DN(3,1)*r_DN(3,1);
const double crLHS323 = r_N[3]*r_N[3];
const double crLHS324 = crLHS3*crLHS323;
const double crLHS325 = crLHS66*crLHS89;
const double crLHS326 = crLHS92*r_N[3];
const double crLHS327 = crLHS74*r_N[3];
const double crLHS328 = crLHS243*crLHS327;
const double crLHS329 = crLHS158*crLHS18 + crLHS19*crLHS325 + crLHS326*crLHS91 - crLHS328*crLHS91;
const double crLHS330 = crLHS313*tau_c;
const double crLHS331 = crLHS109*r_DN(3,0);
const double crLHS332 = crLHS122*crLHS326 - crLHS122*crLHS328 + crLHS158*crLHS36 + crLHS325*crLHS37;
const double crLHS333 = crLHS142*crLHS326 - crLHS142*crLHS328 + crLHS158*crLHS52 + crLHS325*crLHS53;
const double crLHS334 = gauss_weight*(crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_u - crLHS327*crLHS77 - r_N[3]);
const double crLHS335 = crLHS158*crLHS66 + crLHS161*crLHS326 - crLHS161*crLHS328 + crLHS24*crLHS324 + crLHS325*crLHS67;
const double crLHS336 = crLHS313*r_DN(3,1);
const double crLHS337 = crLHS9*r_DN(3,1);
const double crLHS338 = crLHS337*tau_c;
const double crLHS339 = crLHS109*r_DN(3,1);
const double crLHS340 = crLHS269*crLHS327;
const double crLHS341 = crLHS190*crLHS324;
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
rLHS(0,11)+=crLHS58;
rLHS(0,12)+=crLHS61;
rLHS(0,13)+=crLHS10*(crLHS62 + crLHS64 + crLHS69*r_DN(0,0));
rLHS(0,14)+=crLHS10*(crLHS69*r_DN(0,1) + crLHS70 + crLHS71);
rLHS(0,15)+=crLHS72;
rLHS(1,0)+=crLHS78*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS0*crLHS9 + crLHS79*r_DN(0,0) + crLHS81*r_DN(0,1) + crLHS84*crLHS86 + crLHS97);
rLHS(1,2)+=gauss_weight*(crLHS100*r_DN(0,1) + crLHS101 + crLHS104*crLHS86 + crLHS98*r_DN(0,0));
rLHS(1,3)+=-crLHS108*crLHS110;
rLHS(1,4)+=gauss_weight*(-crLHS111 - crLHS112*crLHS31 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u);
rLHS(1,5)+=gauss_weight*(crLHS113*r_DN(0,0) + crLHS115*r_DN(0,1) + crLHS118*crLHS86 + crLHS121 + crLHS123);
rLHS(1,6)+=gauss_weight*(crLHS124*r_DN(0,0) + crLHS126*r_DN(0,1) + crLHS127 + crLHS129*crLHS86);
rLHS(1,7)+=-crLHS110*crLHS131;
rLHS(1,8)+=gauss_weight*(-crLHS112*crLHS48 - crLHS132 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u);
rLHS(1,9)+=gauss_weight*(crLHS133*r_DN(0,0) + crLHS135*r_DN(0,1) + crLHS138*crLHS86 + crLHS141 + crLHS143);
rLHS(1,10)+=gauss_weight*(crLHS144*r_DN(0,0) + crLHS146*r_DN(0,1) + crLHS147 + crLHS149*crLHS86);
rLHS(1,11)+=-crLHS110*crLHS151;
rLHS(1,12)+=gauss_weight*(-crLHS112*crLHS62 - crLHS152 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[0]*tau_u);
rLHS(1,13)+=gauss_weight*(crLHS153*r_DN(0,0) + crLHS155*r_DN(0,1) + crLHS157*crLHS86 + crLHS160 + crLHS162);
rLHS(1,14)+=gauss_weight*(crLHS163*r_DN(0,0) + crLHS165*r_DN(0,1) + crLHS166 + crLHS167*crLHS86);
rLHS(1,15)+=-crLHS110*crLHS169;
rLHS(2,0)+=crLHS78*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS101 + crLHS170*r_DN(0,1) + crLHS172*crLHS84 + crLHS81*r_DN(0,0));
rLHS(2,2)+=gauss_weight*(crLHS1*crLHS9 + crLHS100*r_DN(0,0) + crLHS104*crLHS172 + crLHS173*r_DN(0,1) + crLHS97);
rLHS(2,3)+=-crLHS108*crLHS174;
rLHS(2,4)+=gauss_weight*(-crLHS112*crLHS40 - crLHS175 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u);
rLHS(2,5)+=gauss_weight*(crLHS115*r_DN(0,0) + crLHS118*crLHS172 + crLHS176*r_DN(0,1) + crLHS177);
rLHS(2,6)+=gauss_weight*(crLHS123 + crLHS126*r_DN(0,0) + crLHS129*crLHS172 + crLHS178*r_DN(0,1) + crLHS179);
rLHS(2,7)+=-crLHS131*crLHS174;
rLHS(2,8)+=gauss_weight*(-crLHS112*crLHS56 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u - crLHS180 + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u);
rLHS(2,9)+=gauss_weight*(crLHS135*r_DN(0,0) + crLHS138*crLHS172 + crLHS181*r_DN(0,1) + crLHS182);
rLHS(2,10)+=gauss_weight*(crLHS143 + crLHS146*r_DN(0,0) + crLHS149*crLHS172 + crLHS183*r_DN(0,1) + crLHS184);
rLHS(2,11)+=-crLHS151*crLHS174;
rLHS(2,12)+=gauss_weight*(-crLHS112*crLHS70 + crLHS18*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(3,1)*tau_u - crLHS185 + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[0]*tau_u);
rLHS(2,13)+=gauss_weight*(crLHS155*r_DN(0,0) + crLHS157*crLHS172 + crLHS186*r_DN(0,1) + crLHS187);
rLHS(2,14)+=gauss_weight*(crLHS162 + crLHS165*r_DN(0,0) + crLHS167*crLHS172 + crLHS188*r_DN(0,1) + crLHS189);
rLHS(2,15)+=-crLHS169*crLHS174;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS191*crLHS82;
rLHS(3,2)+=crLHS102*crLHS191;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa - crLHS107*crLHS194 + crLHS13*crLHS3*dp_th_dt - crLHS14*crLHS192 + crLHS18*crLHS193*crLHS25*crLHS6*gamma*p_th*tau_t - crLHS193*crLHS196 + crLHS193*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS193*crLHS25*dp_th_dt*r_N[0]*tau_t);
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS198;
rLHS(3,6)+=crLHS199;
rLHS(3,7)+=-gauss_weight*(-crLHS130*crLHS194 + crLHS18*crLHS200*crLHS25*crLHS6*gamma*p_th*tau_t - crLHS196*crLHS200 + crLHS200*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS200*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS202);
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS204;
rLHS(3,10)+=crLHS205;
rLHS(3,11)+=-gauss_weight*(-crLHS150*crLHS194 + crLHS18*crLHS206*crLHS25*crLHS6*gamma*p_th*tau_t - crLHS196*crLHS206 + crLHS206*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS206*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS208);
rLHS(3,12)+=0;
rLHS(3,13)+=crLHS210;
rLHS(3,14)+=crLHS211;
rLHS(3,15)+=-gauss_weight*(-crLHS168*crLHS194 + crLHS18*crLHS212*crLHS25*crLHS6*gamma*p_th*tau_t - crLHS196*crLHS212 + crLHS212*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS212*crLHS25*dp_th_dt*r_N[0]*tau_t - crLHS214);
rLHS(4,0)+=crLHS30;
rLHS(4,1)+=crLHS10*(crLHS111 + crLHS22*r_DN(1,0) + crLHS34);
rLHS(4,2)+=crLHS10*(crLHS175 + crLHS22*r_DN(1,1) + crLHS41);
rLHS(4,3)+=crLHS44;
rLHS(4,4)+=crLHS11*(crLHS215 + crLHS216);
rLHS(4,5)+=crLHS10*(-crLHS12*crLHS218 + crLHS39*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS218*crLHS23 + crLHS39*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS217*crLHS27;
rLHS(4,8)+=crLHS221;
rLHS(4,9)+=crLHS10*(crLHS222 + crLHS224 + crLHS55*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS225 + crLHS226 + crLHS55*r_DN(1,1));
rLHS(4,11)+=crLHS228;
rLHS(4,12)+=crLHS231;
rLHS(4,13)+=crLHS10*(crLHS232 + crLHS234 + crLHS69*r_DN(1,0));
rLHS(4,14)+=crLHS10*(crLHS235 + crLHS236 + crLHS69*r_DN(1,1));
rLHS(4,15)+=crLHS237;
rLHS(5,0)+=gauss_weight*(-crLHS111*crLHS112 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS31);
rLHS(5,1)+=gauss_weight*(crLHS121 + crLHS239*crLHS84 + crLHS245 + crLHS79*r_DN(1,0) + crLHS81*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS100*r_DN(1,1) + crLHS104*crLHS239 + crLHS177 + crLHS98*r_DN(1,0));
rLHS(5,3)+=-crLHS108*crLHS246;
rLHS(5,4)+=crLHS247*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS113*r_DN(1,0) + crLHS115*r_DN(1,1) + crLHS118*crLHS239 + crLHS215*crLHS9 + crLHS248);
rLHS(5,6)+=gauss_weight*(crLHS124*r_DN(1,0) + crLHS126*r_DN(1,1) + crLHS129*crLHS239 + crLHS249);
rLHS(5,7)+=-crLHS131*crLHS246;
rLHS(5,8)+=gauss_weight*(-crLHS112*crLHS222 - crLHS250 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS133*r_DN(1,0) + crLHS135*r_DN(1,1) + crLHS138*crLHS239 + crLHS252 + crLHS253);
rLHS(5,10)+=gauss_weight*(crLHS144*r_DN(1,0) + crLHS146*r_DN(1,1) + crLHS149*crLHS239 + crLHS254);
rLHS(5,11)+=-crLHS151*crLHS246;
rLHS(5,12)+=gauss_weight*(-crLHS112*crLHS232 - crLHS255 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[1]*tau_u);
rLHS(5,13)+=gauss_weight*(crLHS153*r_DN(1,0) + crLHS155*r_DN(1,1) + crLHS157*crLHS239 + crLHS257 + crLHS258);
rLHS(5,14)+=gauss_weight*(crLHS163*r_DN(1,0) + crLHS165*r_DN(1,1) + crLHS167*crLHS239 + crLHS259);
rLHS(5,15)+=-crLHS169*crLHS246;
rLHS(6,0)+=gauss_weight*(-crLHS112*crLHS175 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS40);
rLHS(6,1)+=gauss_weight*(crLHS127 + crLHS170*r_DN(1,1) + crLHS261*crLHS84 + crLHS81*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS100*r_DN(1,0) + crLHS104*crLHS261 + crLHS173*r_DN(1,1) + crLHS179 + crLHS245);
rLHS(6,3)+=-crLHS108*crLHS262;
rLHS(6,4)+=crLHS247*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS115*r_DN(1,0) + crLHS118*crLHS261 + crLHS176*r_DN(1,1) + crLHS249);
rLHS(6,6)+=gauss_weight*(crLHS126*r_DN(1,0) + crLHS129*crLHS261 + crLHS178*r_DN(1,1) + crLHS216*crLHS9 + crLHS248);
rLHS(6,7)+=-crLHS131*crLHS262;
rLHS(6,8)+=gauss_weight*(-crLHS112*crLHS225 - crLHS263 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS135*r_DN(1,0) + crLHS138*crLHS261 + crLHS181*r_DN(1,1) + crLHS264);
rLHS(6,10)+=gauss_weight*(crLHS146*r_DN(1,0) + crLHS149*crLHS261 + crLHS183*r_DN(1,1) + crLHS253 + crLHS265);
rLHS(6,11)+=-crLHS151*crLHS262;
rLHS(6,12)+=gauss_weight*(-crLHS112*crLHS235 - crLHS266 + crLHS3*crLHS36*crLHS4*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[1]*tau_u);
rLHS(6,13)+=gauss_weight*(crLHS155*r_DN(1,0) + crLHS157*crLHS261 + crLHS186*r_DN(1,1) + crLHS267);
rLHS(6,14)+=gauss_weight*(crLHS165*r_DN(1,0) + crLHS167*crLHS261 + crLHS188*r_DN(1,1) + crLHS258 + crLHS268);
rLHS(6,15)+=-crLHS169*crLHS262;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS198;
rLHS(7,2)+=crLHS199;
rLHS(7,3)+=-gauss_weight*(-crLHS107*crLHS201 + crLHS193*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS193*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS193*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS193*crLHS270 - crLHS202);
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS271*crLHS82;
rLHS(7,6)+=crLHS102*crLHS271;
rLHS(7,7)+=-gauss_weight*(-crLHS130*crLHS201 - crLHS192*crLHS218 + crLHS200*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS200*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS200*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS200*crLHS270 - crLHS215*kappa - crLHS216*kappa + crLHS217*crLHS3*dp_th_dt);
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS272;
rLHS(7,10)+=crLHS273;
rLHS(7,11)+=-gauss_weight*(-crLHS150*crLHS201 + crLHS206*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS206*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS206*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS206*crLHS270 - crLHS274);
rLHS(7,12)+=0;
rLHS(7,13)+=crLHS275;
rLHS(7,14)+=crLHS276;
rLHS(7,15)+=-gauss_weight*(-crLHS168*crLHS201 + crLHS212*crLHS25*crLHS36*crLHS6*gamma*p_th*tau_t + crLHS212*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS212*crLHS25*dp_th_dt*r_N[1]*tau_t - crLHS212*crLHS270 - crLHS277);
rLHS(8,0)+=crLHS47;
rLHS(8,1)+=crLHS10*(crLHS132 + crLHS22*r_DN(2,0) + crLHS50);
rLHS(8,2)+=crLHS10*(crLHS180 + crLHS22*r_DN(2,1) + crLHS57);
rLHS(8,3)+=crLHS58;
rLHS(8,4)+=crLHS221;
rLHS(8,5)+=crLHS10*(crLHS224 + crLHS250 + crLHS39*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS226 + crLHS263 + crLHS39*r_DN(2,1));
rLHS(8,7)+=crLHS228;
rLHS(8,8)+=crLHS11*(crLHS278 + crLHS279);
rLHS(8,9)+=crLHS10*(-crLHS12*crLHS281 + crLHS55*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS23*crLHS281 + crLHS55*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS27*crLHS280;
rLHS(8,12)+=crLHS284;
rLHS(8,13)+=crLHS10*(crLHS285 + crLHS287 + crLHS69*r_DN(2,0));
rLHS(8,14)+=crLHS10*(crLHS288 + crLHS289 + crLHS69*r_DN(2,1));
rLHS(8,15)+=crLHS290;
rLHS(9,0)+=gauss_weight*(-crLHS112*crLHS132 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS48);
rLHS(9,1)+=gauss_weight*(crLHS141 + crLHS292*crLHS84 + crLHS297 + crLHS79*r_DN(2,0) + crLHS81*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS100*r_DN(2,1) + crLHS104*crLHS292 + crLHS182 + crLHS98*r_DN(2,0));
rLHS(9,3)+=-crLHS108*crLHS298;
rLHS(9,4)+=gauss_weight*(-crLHS112*crLHS250 - crLHS222 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS113*r_DN(2,0) + crLHS115*r_DN(2,1) + crLHS118*crLHS292 + crLHS252 + crLHS299);
rLHS(9,6)+=gauss_weight*(crLHS124*r_DN(2,0) + crLHS126*r_DN(2,1) + crLHS129*crLHS292 + crLHS264);
rLHS(9,7)+=-crLHS131*crLHS298;
rLHS(9,8)+=crLHS300*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS133*r_DN(2,0) + crLHS135*r_DN(2,1) + crLHS138*crLHS292 + crLHS278*crLHS9 + crLHS301);
rLHS(9,10)+=gauss_weight*(crLHS144*r_DN(2,0) + crLHS146*r_DN(2,1) + crLHS149*crLHS292 + crLHS302);
rLHS(9,11)+=-crLHS151*crLHS298;
rLHS(9,12)+=gauss_weight*(-crLHS112*crLHS285 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[2]*tau_u - crLHS303);
rLHS(9,13)+=gauss_weight*(crLHS153*r_DN(2,0) + crLHS155*r_DN(2,1) + crLHS157*crLHS292 + crLHS305 + crLHS306);
rLHS(9,14)+=gauss_weight*(crLHS163*r_DN(2,0) + crLHS165*r_DN(2,1) + crLHS167*crLHS292 + crLHS307);
rLHS(9,15)+=-crLHS169*crLHS298;
rLHS(10,0)+=gauss_weight*(-crLHS112*crLHS180 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS56);
rLHS(10,1)+=gauss_weight*(crLHS147 + crLHS170*r_DN(2,1) + crLHS308*crLHS309 + crLHS81*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS100*r_DN(2,0) + crLHS104*crLHS310 + crLHS173*r_DN(2,1) + crLHS184 + crLHS297);
rLHS(10,3)+=-crLHS108*crLHS311;
rLHS(10,4)+=gauss_weight*(-crLHS112*crLHS263 - crLHS225 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS115*r_DN(2,0) + crLHS118*crLHS310 + crLHS176*r_DN(2,1) + crLHS254);
rLHS(10,6)+=gauss_weight*(crLHS126*r_DN(2,0) + crLHS129*crLHS310 + crLHS178*r_DN(2,1) + crLHS265 + crLHS299);
rLHS(10,7)+=-crLHS131*crLHS311;
rLHS(10,8)+=crLHS300*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS135*r_DN(2,0) + crLHS138*crLHS310 + crLHS181*r_DN(2,1) + crLHS302);
rLHS(10,10)+=gauss_weight*(crLHS146*r_DN(2,0) + crLHS149*crLHS310 + crLHS183*r_DN(2,1) + crLHS279*crLHS9 + crLHS301);
rLHS(10,11)+=-crLHS151*crLHS311;
rLHS(10,12)+=gauss_weight*(-crLHS112*crLHS288 + crLHS3*crLHS4*crLHS52*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[2]*tau_u - crLHS312);
rLHS(10,13)+=gauss_weight*(crLHS155*r_DN(2,0) + crLHS157*crLHS310 + crLHS186*r_DN(2,1) + crLHS314);
rLHS(10,14)+=gauss_weight*(crLHS165*r_DN(2,0) + crLHS167*crLHS310 + crLHS188*r_DN(2,1) + crLHS306 + crLHS315);
rLHS(10,15)+=-crLHS169*crLHS311;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS204;
rLHS(11,2)+=crLHS205;
rLHS(11,3)+=-gauss_weight*(-crLHS107*crLHS207 + crLHS193*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS193*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS193*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS193*crLHS316 - crLHS208);
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS272;
rLHS(11,6)+=crLHS273;
rLHS(11,7)+=-gauss_weight*(-crLHS130*crLHS207 + crLHS200*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS200*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS200*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS200*crLHS316 - crLHS274);
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS317*crLHS82;
rLHS(11,10)+=crLHS102*crLHS317;
rLHS(11,11)+=-gauss_weight*(-crLHS150*crLHS207 - crLHS192*crLHS281 + crLHS206*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS206*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS206*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS206*crLHS316 - crLHS278*kappa - crLHS279*kappa + crLHS280*crLHS3*dp_th_dt);
rLHS(11,12)+=0;
rLHS(11,13)+=crLHS318;
rLHS(11,14)+=crLHS319;
rLHS(11,15)+=-gauss_weight*(-crLHS168*crLHS207 + crLHS212*crLHS25*crLHS52*crLHS6*gamma*p_th*tau_t + crLHS212*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS212*crLHS25*dp_th_dt*r_N[2]*tau_t - crLHS212*crLHS316 - crLHS320);
rLHS(12,0)+=crLHS61;
rLHS(12,1)+=crLHS10*(crLHS152 + crLHS22*r_DN(3,0) + crLHS64);
rLHS(12,2)+=crLHS10*(crLHS185 + crLHS22*r_DN(3,1) + crLHS71);
rLHS(12,3)+=crLHS72;
rLHS(12,4)+=crLHS231;
rLHS(12,5)+=crLHS10*(crLHS234 + crLHS255 + crLHS39*r_DN(3,0));
rLHS(12,6)+=crLHS10*(crLHS236 + crLHS266 + crLHS39*r_DN(3,1));
rLHS(12,7)+=crLHS237;
rLHS(12,8)+=crLHS284;
rLHS(12,9)+=crLHS10*(crLHS287 + crLHS303 + crLHS55*r_DN(3,0));
rLHS(12,10)+=crLHS10*(crLHS289 + crLHS312 + crLHS55*r_DN(3,1));
rLHS(12,11)+=crLHS290;
rLHS(12,12)+=crLHS11*(crLHS321 + crLHS322);
rLHS(12,13)+=crLHS10*(-crLHS12*crLHS324 + crLHS69*r_DN(3,0) + r_DN(3,0)*r_N[3]);
rLHS(12,14)+=crLHS10*(-crLHS23*crLHS324 + crLHS69*r_DN(3,1) + r_DN(3,1)*r_N[3]);
rLHS(12,15)+=-crLHS27*crLHS323;
rLHS(13,0)+=gauss_weight*(-crLHS112*crLHS152 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[3]*tau_u - crLHS62);
rLHS(13,1)+=gauss_weight*(crLHS160 + crLHS308*crLHS313 + crLHS329 + crLHS79*r_DN(3,0) + crLHS81*r_DN(3,1));
rLHS(13,2)+=gauss_weight*(crLHS100*r_DN(3,1) + crLHS104*crLHS330 + crLHS187 + crLHS98*r_DN(3,0));
rLHS(13,3)+=-crLHS108*crLHS331;
rLHS(13,4)+=gauss_weight*(-crLHS112*crLHS255 - crLHS232 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[3]*tau_u);
rLHS(13,5)+=gauss_weight*(crLHS113*r_DN(3,0) + crLHS115*r_DN(3,1) + crLHS118*crLHS330 + crLHS257 + crLHS332);
rLHS(13,6)+=gauss_weight*(crLHS124*r_DN(3,0) + crLHS126*r_DN(3,1) + crLHS129*crLHS330 + crLHS267);
rLHS(13,7)+=-crLHS131*crLHS331;
rLHS(13,8)+=gauss_weight*(-crLHS112*crLHS303 - crLHS285 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[3]*tau_u);
rLHS(13,9)+=gauss_weight*(crLHS133*r_DN(3,0) + crLHS135*r_DN(3,1) + crLHS138*crLHS330 + crLHS305 + crLHS333);
rLHS(13,10)+=gauss_weight*(crLHS144*r_DN(3,0) + crLHS146*r_DN(3,1) + crLHS149*crLHS330 + crLHS314);
rLHS(13,11)+=-crLHS151*crLHS331;
rLHS(13,12)+=crLHS334*r_DN(3,0);
rLHS(13,13)+=gauss_weight*(crLHS153*r_DN(3,0) + crLHS155*r_DN(3,1) + crLHS157*crLHS330 + crLHS321*crLHS9 + crLHS335);
rLHS(13,14)+=gauss_weight*(crLHS163*r_DN(3,0) + crLHS165*r_DN(3,1) + crLHS167*crLHS330 + crLHS336);
rLHS(13,15)+=-crLHS169*crLHS331;
rLHS(14,0)+=gauss_weight*(-crLHS112*crLHS185 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[3]*tau_u - crLHS70);
rLHS(14,1)+=gauss_weight*(crLHS166 + crLHS170*r_DN(3,1) + crLHS308*crLHS337 + crLHS81*r_DN(3,0));
rLHS(14,2)+=gauss_weight*(crLHS100*r_DN(3,0) + crLHS104*crLHS338 + crLHS173*r_DN(3,1) + crLHS189 + crLHS329);
rLHS(14,3)+=-crLHS108*crLHS339;
rLHS(14,4)+=gauss_weight*(-crLHS112*crLHS266 - crLHS235 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[3]*tau_u);
rLHS(14,5)+=gauss_weight*(crLHS115*r_DN(3,0) + crLHS118*crLHS338 + crLHS176*r_DN(3,1) + crLHS259);
rLHS(14,6)+=gauss_weight*(crLHS126*r_DN(3,0) + crLHS129*crLHS338 + crLHS178*r_DN(3,1) + crLHS268 + crLHS332);
rLHS(14,7)+=-crLHS131*crLHS339;
rLHS(14,8)+=gauss_weight*(-crLHS112*crLHS312 - crLHS288 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[3]*tau_u);
rLHS(14,9)+=gauss_weight*(crLHS135*r_DN(3,0) + crLHS138*crLHS338 + crLHS181*r_DN(3,1) + crLHS307);
rLHS(14,10)+=gauss_weight*(crLHS146*r_DN(3,0) + crLHS149*crLHS338 + crLHS183*r_DN(3,1) + crLHS315 + crLHS333);
rLHS(14,11)+=-crLHS151*crLHS339;
rLHS(14,12)+=crLHS334*r_DN(3,1);
rLHS(14,13)+=gauss_weight*(crLHS155*r_DN(3,0) + crLHS157*crLHS338 + crLHS186*r_DN(3,1) + crLHS336);
rLHS(14,14)+=gauss_weight*(crLHS165*r_DN(3,0) + crLHS167*crLHS338 + crLHS188*r_DN(3,1) + crLHS322*crLHS9 + crLHS335);
rLHS(14,15)+=-crLHS169*crLHS339;
rLHS(15,0)+=0;
rLHS(15,1)+=crLHS210;
rLHS(15,2)+=crLHS211;
rLHS(15,3)+=-gauss_weight*(-crLHS107*crLHS213 + crLHS193*crLHS25*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS193*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS193*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS193*crLHS340 - crLHS214);
rLHS(15,4)+=0;
rLHS(15,5)+=crLHS275;
rLHS(15,6)+=crLHS276;
rLHS(15,7)+=-gauss_weight*(-crLHS130*crLHS213 + crLHS200*crLHS25*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS200*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS200*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS200*crLHS340 - crLHS277);
rLHS(15,8)+=0;
rLHS(15,9)+=crLHS318;
rLHS(15,10)+=crLHS319;
rLHS(15,11)+=-gauss_weight*(-crLHS150*crLHS213 + crLHS206*crLHS25*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS206*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS206*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS206*crLHS340 - crLHS320);
rLHS(15,12)+=0;
rLHS(15,13)+=crLHS341*crLHS82;
rLHS(15,14)+=crLHS102*crLHS341;
rLHS(15,15)+=-gauss_weight*(-crLHS168*crLHS213 - crLHS192*crLHS324 + crLHS212*crLHS25*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS212*crLHS25*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS212*crLHS25*dp_th_dt*r_N[3]*tau_t - crLHS212*crLHS340 + crLHS3*crLHS323*dp_th_dt - crLHS321*kappa - crLHS322*kappa);

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
const double crRHS35 = crRHS30*r_DN(0,0);
const double crRHS36 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1);
const double crRHS37 = crRHS36*tau_u;
const double crRHS38 = crRHS34*crRHS37;
const double crRHS39 = crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1);
const double crRHS40 = crRHS21*tau_u;
const double crRHS41 = crRHS39*crRHS40;
const double crRHS42 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crRHS43 = crRHS10*crRHS42;
const double crRHS44 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crRHS45 = crRHS12*crRHS44;
const double crRHS46 = tau_c*(-crRHS2 + crRHS43*crRHS6 + crRHS45*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS47 = (crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS48 = crRHS47*r_N[0];
const double crRHS49 = crRHS29*crRHS48;
const double crRHS50 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS51 = crRHS30*r_DN(0,1);
const double crRHS52 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2];
const double crRHS53 = crRHS42*kappa;
const double crRHS54 = crRHS44*kappa;
const double crRHS55 = crRHS5*dp_th_dt;
const double crRHS56 = crRHS55*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2]);
const double crRHS57 = crRHS19*crRHS7;
const double crRHS58 = crRHS20*(crRHS43 + crRHS45);
const double crRHS59 = tau_t*(-crRHS20*(crRHS14*crRHS42 + crRHS15*crRHS44 + crRHS3) + crRHS52 + crRHS56);
const double crRHS60 = crRHS55*crRHS59;
const double crRHS61 = crRHS20*crRHS59;
const double crRHS62 = crRHS36*crRHS61;
const double crRHS63 = crRHS19*crRHS59*p_th;
const double crRHS64 = crRHS21*r_N[1];
const double crRHS65 = crRHS30*r_DN(1,0);
const double crRHS66 = crRHS37*crRHS64;
const double crRHS67 = crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1);
const double crRHS68 = crRHS40*crRHS67;
const double crRHS69 = crRHS29*crRHS47;
const double crRHS70 = crRHS69*r_N[1];
const double crRHS71 = crRHS30*r_DN(1,1);
const double crRHS72 = crRHS47*crRHS63;
const double crRHS73 = crRHS21*r_N[2];
const double crRHS74 = crRHS30*r_DN(2,0);
const double crRHS75 = crRHS37*crRHS73;
const double crRHS76 = crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1);
const double crRHS77 = crRHS40*crRHS76;
const double crRHS78 = crRHS69*r_N[2];
const double crRHS79 = crRHS30*r_DN(2,1);
rRHS[0]+=-crRHS31*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS24*r_DN(0,0) + crRHS28*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS16*crRHS34 + crRHS17*crRHS34 + crRHS2*crRHS35 + crRHS22*crRHS38 + crRHS22*crRHS41 - crRHS24*crRHS49 - crRHS32*r_DN(0,0) - crRHS33*crRHS34 - crRHS35*crRHS46 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS2*crRHS51 + crRHS25*crRHS34 + crRHS26*crRHS34 + crRHS27*crRHS38 + crRHS27*crRHS41 - crRHS28*crRHS49 - crRHS32*r_DN(0,1) - crRHS34*crRHS50 - crRHS46*crRHS51 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS39*crRHS61 - crRHS48*crRHS63 + crRHS52*r_N[0] - crRHS53*r_DN(0,0) - crRHS54*r_DN(0,1) + crRHS56*r_N[0] - crRHS57*r_N[0] - crRHS58*r_N[0] + crRHS60*r_N[0] + crRHS62*r_N[0]);
rRHS[4]+=-crRHS31*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS24*r_DN(1,0) + crRHS28*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS16*crRHS64 + crRHS17*crRHS64 + crRHS2*crRHS65 + crRHS22*crRHS66 + crRHS22*crRHS68 - crRHS24*crRHS70 - crRHS32*r_DN(1,0) - crRHS33*crRHS64 - crRHS46*crRHS65 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS2*crRHS71 + crRHS25*crRHS64 + crRHS26*crRHS64 + crRHS27*crRHS66 + crRHS27*crRHS68 - crRHS28*crRHS70 - crRHS32*r_DN(1,1) - crRHS46*crRHS71 - crRHS50*crRHS64 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS52*r_N[1] - crRHS53*r_DN(1,0) - crRHS54*r_DN(1,1) + crRHS56*r_N[1] - crRHS57*r_N[1] - crRHS58*r_N[1] + crRHS60*r_N[1] + crRHS61*crRHS67 + crRHS62*r_N[1] - crRHS72*r_N[1]);
rRHS[8]+=-crRHS31*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS24*r_DN(2,0) + crRHS28*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS16*crRHS73 + crRHS17*crRHS73 + crRHS2*crRHS74 + crRHS22*crRHS75 + crRHS22*crRHS77 - crRHS24*crRHS78 - crRHS32*r_DN(2,0) - crRHS33*crRHS73 - crRHS46*crRHS74 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS2*crRHS79 + crRHS25*crRHS73 + crRHS26*crRHS73 + crRHS27*crRHS75 + crRHS27*crRHS77 - crRHS28*crRHS78 - crRHS32*r_DN(2,1) - crRHS46*crRHS79 - crRHS50*crRHS73 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS52*r_N[2] - crRHS53*r_DN(2,0) - crRHS54*r_DN(2,1) + crRHS56*r_N[2] - crRHS57*r_N[2] - crRHS58*r_N[2] + crRHS60*r_N[2] + crRHS61*crRHS76 + crRHS62*r_N[2] - crRHS72*r_N[2]);

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
const double crRHS35 = crRHS30*r_DN(0,0);
const double crRHS36 = lin_u_conv(0,0)*r_DN(0,0) + lin_u_conv(0,1)*r_DN(0,1) + lin_u_conv(1,0)*r_DN(1,0) + lin_u_conv(1,1)*r_DN(1,1) + lin_u_conv(2,0)*r_DN(2,0) + lin_u_conv(2,1)*r_DN(2,1) + lin_u_conv(3,0)*r_DN(3,0) + lin_u_conv(3,1)*r_DN(3,1);
const double crRHS37 = crRHS36*tau_u;
const double crRHS38 = crRHS34*crRHS37;
const double crRHS39 = crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1);
const double crRHS40 = crRHS21*tau_u;
const double crRHS41 = crRHS39*crRHS40;
const double crRHS42 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crRHS43 = crRHS10*crRHS42;
const double crRHS44 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crRHS45 = crRHS12*crRHS44;
const double crRHS46 = tau_c*(-crRHS2 + crRHS43*crRHS6 + crRHS45*crRHS6 + crRHS7 - dp_th_dt);
const double crRHS47 = (crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS48 = crRHS47*r_N[0];
const double crRHS49 = crRHS29*crRHS48;
const double crRHS50 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS51 = crRHS30*r_DN(0,1);
const double crRHS52 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2] + r_N[3]*r_heat_fl[3];
const double crRHS53 = crRHS42*kappa;
const double crRHS54 = crRHS44*kappa;
const double crRHS55 = crRHS5*dp_th_dt;
const double crRHS56 = crRHS55*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3]);
const double crRHS57 = crRHS19*crRHS7;
const double crRHS58 = crRHS20*(crRHS43 + crRHS45);
const double crRHS59 = tau_t*(-crRHS20*(crRHS14*crRHS42 + crRHS15*crRHS44 + crRHS3) + crRHS52 + crRHS56);
const double crRHS60 = crRHS55*crRHS59;
const double crRHS61 = crRHS20*crRHS59;
const double crRHS62 = crRHS36*crRHS61;
const double crRHS63 = crRHS19*crRHS59*p_th;
const double crRHS64 = crRHS21*r_N[1];
const double crRHS65 = crRHS30*r_DN(1,0);
const double crRHS66 = crRHS37*crRHS64;
const double crRHS67 = crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1);
const double crRHS68 = crRHS40*crRHS67;
const double crRHS69 = crRHS29*crRHS47;
const double crRHS70 = crRHS69*r_N[1];
const double crRHS71 = crRHS30*r_DN(1,1);
const double crRHS72 = crRHS47*crRHS63;
const double crRHS73 = crRHS21*r_N[2];
const double crRHS74 = crRHS30*r_DN(2,0);
const double crRHS75 = crRHS37*crRHS73;
const double crRHS76 = crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1);
const double crRHS77 = crRHS40*crRHS76;
const double crRHS78 = crRHS69*r_N[2];
const double crRHS79 = crRHS30*r_DN(2,1);
const double crRHS80 = crRHS21*r_N[3];
const double crRHS81 = crRHS30*r_DN(3,0);
const double crRHS82 = crRHS37*crRHS80;
const double crRHS83 = crRHS14*r_DN(3,0) + crRHS15*r_DN(3,1);
const double crRHS84 = crRHS40*crRHS83;
const double crRHS85 = crRHS69*r_N[3];
const double crRHS86 = crRHS30*r_DN(3,1);
rRHS[0]+=-crRHS31*(-crRHS13*r_N[0] + crRHS2*r_N[0] + crRHS24*r_DN(0,0) + crRHS28*r_DN(0,1) + crRHS8*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS16*crRHS34 + crRHS17*crRHS34 + crRHS2*crRHS35 + crRHS22*crRHS38 + crRHS22*crRHS41 - crRHS24*crRHS49 - crRHS32*r_DN(0,0) - crRHS33*crRHS34 - crRHS35*crRHS46 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS2*crRHS51 + crRHS25*crRHS34 + crRHS26*crRHS34 + crRHS27*crRHS38 + crRHS27*crRHS41 - crRHS28*crRHS49 - crRHS32*r_DN(0,1) - crRHS34*crRHS50 - crRHS46*crRHS51 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS39*crRHS61 - crRHS48*crRHS63 + crRHS52*r_N[0] - crRHS53*r_DN(0,0) - crRHS54*r_DN(0,1) + crRHS56*r_N[0] - crRHS57*r_N[0] - crRHS58*r_N[0] + crRHS60*r_N[0] + crRHS62*r_N[0]);
rRHS[4]+=-crRHS31*(-crRHS13*r_N[1] + crRHS2*r_N[1] + crRHS24*r_DN(1,0) + crRHS28*r_DN(1,1) + crRHS8*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS16*crRHS64 + crRHS17*crRHS64 + crRHS2*crRHS65 + crRHS22*crRHS66 + crRHS22*crRHS68 - crRHS24*crRHS70 - crRHS32*r_DN(1,0) - crRHS33*crRHS64 - crRHS46*crRHS65 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS2*crRHS71 + crRHS25*crRHS64 + crRHS26*crRHS64 + crRHS27*crRHS66 + crRHS27*crRHS68 - crRHS28*crRHS70 - crRHS32*r_DN(1,1) - crRHS46*crRHS71 - crRHS50*crRHS64 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS52*r_N[1] - crRHS53*r_DN(1,0) - crRHS54*r_DN(1,1) + crRHS56*r_N[1] - crRHS57*r_N[1] - crRHS58*r_N[1] + crRHS60*r_N[1] + crRHS61*crRHS67 + crRHS62*r_N[1] - crRHS72*r_N[1]);
rRHS[8]+=-crRHS31*(-crRHS13*r_N[2] + crRHS2*r_N[2] + crRHS24*r_DN(2,0) + crRHS28*r_DN(2,1) + crRHS8*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS16*crRHS73 + crRHS17*crRHS73 + crRHS2*crRHS74 + crRHS22*crRHS75 + crRHS22*crRHS77 - crRHS24*crRHS78 - crRHS32*r_DN(2,0) - crRHS33*crRHS73 - crRHS46*crRHS74 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS2*crRHS79 + crRHS25*crRHS73 + crRHS26*crRHS73 + crRHS27*crRHS75 + crRHS27*crRHS77 - crRHS28*crRHS78 - crRHS32*r_DN(2,1) - crRHS46*crRHS79 - crRHS50*crRHS73 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS52*r_N[2] - crRHS53*r_DN(2,0) - crRHS54*r_DN(2,1) + crRHS56*r_N[2] - crRHS57*r_N[2] - crRHS58*r_N[2] + crRHS60*r_N[2] + crRHS61*crRHS76 + crRHS62*r_N[2] - crRHS72*r_N[2]);
rRHS[12]+=-crRHS31*(-crRHS13*r_N[3] + crRHS2*r_N[3] + crRHS24*r_DN(3,0) + crRHS28*r_DN(3,1) + crRHS8*r_N[3]);
rRHS[13]+=-gauss_weight*(crRHS16*crRHS80 + crRHS17*crRHS80 + crRHS2*crRHS81 + crRHS22*crRHS82 + crRHS22*crRHS84 - crRHS24*crRHS85 - crRHS32*r_DN(3,0) - crRHS33*crRHS80 - crRHS46*crRHS81 + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-gauss_weight*(crRHS2*crRHS86 + crRHS25*crRHS80 + crRHS26*crRHS80 + crRHS27*crRHS82 + crRHS27*crRHS84 - crRHS28*crRHS85 - crRHS32*r_DN(3,1) - crRHS46*crRHS86 - crRHS50*crRHS80 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=gauss_weight*(crRHS52*r_N[3] - crRHS53*r_DN(3,0) - crRHS54*r_DN(3,1) + crRHS56*r_N[3] - crRHS57*r_N[3] - crRHS58*r_N[3] + crRHS60*r_N[3] + crRHS61*crRHS83 + crRHS62*r_N[3] - crRHS72*r_N[3]);

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