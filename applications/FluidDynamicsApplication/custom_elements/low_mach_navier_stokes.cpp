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
const double crLHS12 = r_DN(0,0)*r_N[0];
const double crLHS13 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS14 = r_N[0]*r_N[0];
const double crLHS15 = crLHS14*crLHS3;
const double crLHS16 = bdf0*r_N[0];
const double crLHS17 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crLHS18 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crLHS19 = crLHS17*r_DN(0,0) + crLHS18*r_DN(0,1);
const double crLHS20 = crLHS16 + crLHS19;
const double crLHS21 = crLHS20*crLHS7;
const double crLHS22 = crLHS3*crLHS4*tau_u;
const double crLHS23 = crLHS21*crLHS22;
const double crLHS24 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS25 = bdf0*crLHS8;
const double crLHS26 = 1.0/(crLHS2*crLHS2);
const double crLHS27 = crLHS26*gauss_weight;
const double crLHS28 = crLHS25*crLHS27;
const double crLHS29 = r_DN(0,0)*r_DN(1,0);
const double crLHS30 = r_DN(0,1)*r_DN(1,1);
const double crLHS31 = crLHS11*(crLHS29 + crLHS30);
const double crLHS32 = r_DN(1,0)*r_N[0];
const double crLHS33 = crLHS3*r_N[0];
const double crLHS34 = crLHS13*crLHS33;
const double crLHS35 = -crLHS34*r_N[1];
const double crLHS36 = bdf0*r_N[1];
const double crLHS37 = crLHS17*r_DN(1,0) + crLHS18*r_DN(1,1);
const double crLHS38 = crLHS36 + crLHS37;
const double crLHS39 = crLHS38*crLHS7;
const double crLHS40 = crLHS22*crLHS39;
const double crLHS41 = r_DN(1,1)*r_N[0];
const double crLHS42 = crLHS24*crLHS33;
const double crLHS43 = -crLHS42*r_N[1];
const double crLHS44 = crLHS27*crLHS8;
const double crLHS45 = crLHS16*crLHS44;
const double crLHS46 = -crLHS45*r_N[1];
const double crLHS47 = r_DN(0,0)*r_DN(2,0);
const double crLHS48 = r_DN(0,1)*r_DN(2,1);
const double crLHS49 = crLHS11*(crLHS47 + crLHS48);
const double crLHS50 = r_DN(2,0)*r_N[0];
const double crLHS51 = -crLHS34*r_N[2];
const double crLHS52 = bdf0*r_N[2];
const double crLHS53 = crLHS17*r_DN(2,0) + crLHS18*r_DN(2,1);
const double crLHS54 = crLHS52 + crLHS53;
const double crLHS55 = crLHS54*crLHS7;
const double crLHS56 = crLHS22*crLHS55;
const double crLHS57 = r_DN(2,1)*r_N[0];
const double crLHS58 = -crLHS42*r_N[2];
const double crLHS59 = -crLHS45*r_N[2];
const double crLHS60 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1);
const double crLHS61 = crLHS13*crLHS17 + crLHS18*crLHS24;
const double crLHS62 = crLHS61*r_N[0];
const double crLHS63 = crLHS26*tau_u;
const double crLHS64 = crLHS63*crLHS8;
const double crLHS65 = gauss_weight*(crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_N[0]*tau_u - crLHS62*crLHS64 - r_N[0]);
const double crLHS66 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS67 = r_C(0,2)*r_DN(0,0);
const double crLHS68 = crLHS67 + r_C(2,2)*r_DN(0,1);
const double crLHS69 = -crLHS34 + r_DN(0,0);
const double crLHS70 = crLHS9*r_DN(0,0);
const double crLHS71 = crLHS70*tau_c;
const double crLHS72 = crLHS33*crLHS8;
const double crLHS73 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS74 = crLHS63*crLHS73;
const double crLHS75 = crLHS19*crLHS74;
const double crLHS76 = crLHS20*crLHS73;
const double crLHS77 = crLHS26*crLHS60*tau_u;
const double crLHS78 = crLHS77*r_N[0];
const double crLHS79 = 1.0/(crLHS2*crLHS2*crLHS2);
const double crLHS80 = crLHS62*crLHS79;
const double crLHS81 = crLHS80*tau_u;
const double crLHS82 = crLHS15*crLHS25 + crLHS19*crLHS72 + crLHS20*crLHS75 + crLHS76*crLHS78 - crLHS76*crLHS81;
const double crLHS83 = crLHS67 + r_C(0,1)*r_DN(0,1);
const double crLHS84 = r_C(1,2)*r_DN(0,1);
const double crLHS85 = crLHS84 + r_C(2,2)*r_DN(0,0);
const double crLHS86 = crLHS70*r_DN(0,1);
const double crLHS87 = -crLHS42 + r_DN(0,1);
const double crLHS88 = r_DN(0,0)*r_N[1];
const double crLHS89 = crLHS61*crLHS64;
const double crLHS90 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS91 = r_C(0,2)*r_DN(1,0);
const double crLHS92 = crLHS91 + r_C(2,2)*r_DN(1,1);
const double crLHS93 = crLHS3*r_N[1];
const double crLHS94 = crLHS13*crLHS93;
const double crLHS95 = -crLHS94 + r_DN(1,0);
const double crLHS96 = crLHS8*crLHS93;
const double crLHS97 = crLHS16*crLHS96;
const double crLHS98 = crLHS29*crLHS9 + crLHS97;
const double crLHS99 = crLHS38*crLHS73;
const double crLHS100 = crLHS37*crLHS72 + crLHS38*crLHS75 + crLHS78*crLHS99 - crLHS81*crLHS99;
const double crLHS101 = crLHS91 + r_C(0,1)*r_DN(1,1);
const double crLHS102 = r_C(1,2)*r_DN(1,1);
const double crLHS103 = crLHS102 + r_C(2,2)*r_DN(1,0);
const double crLHS104 = crLHS70*r_DN(1,1);
const double crLHS105 = crLHS24*crLHS93;
const double crLHS106 = -crLHS105 + r_DN(1,1);
const double crLHS107 = crLHS36*crLHS44;
const double crLHS108 = r_DN(0,0)*tau_c;
const double crLHS109 = r_DN(0,0)*r_N[2];
const double crLHS110 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS111 = r_C(0,2)*r_DN(2,0);
const double crLHS112 = crLHS111 + r_C(2,2)*r_DN(2,1);
const double crLHS113 = crLHS3*r_N[2];
const double crLHS114 = -crLHS113*crLHS13 + r_DN(2,0);
const double crLHS115 = crLHS113*crLHS8;
const double crLHS116 = crLHS115*crLHS16;
const double crLHS117 = crLHS116 + crLHS47*crLHS9;
const double crLHS118 = crLHS54*crLHS73;
const double crLHS119 = crLHS118*crLHS78 - crLHS118*crLHS81 + crLHS53*crLHS72 + crLHS54*crLHS75;
const double crLHS120 = crLHS111 + r_C(0,1)*r_DN(2,1);
const double crLHS121 = r_C(1,2)*r_DN(2,1);
const double crLHS122 = crLHS121 + r_C(2,2)*r_DN(2,0);
const double crLHS123 = crLHS70*r_DN(2,1);
const double crLHS124 = -crLHS113*crLHS24 + r_DN(2,1);
const double crLHS125 = crLHS44*crLHS52;
const double crLHS126 = crLHS84 + r_C(0,1)*r_DN(0,0);
const double crLHS127 = crLHS9*r_DN(0,1);
const double crLHS128 = crLHS127*tau_c;
const double crLHS129 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS130 = r_DN(0,1)*tau_c;
const double crLHS131 = r_DN(0,1)*r_N[1];
const double crLHS132 = crLHS102 + r_C(0,1)*r_DN(1,0);
const double crLHS133 = crLHS127*r_DN(1,0);
const double crLHS134 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS135 = crLHS30*crLHS9 + crLHS97;
const double crLHS136 = r_DN(0,1)*r_N[2];
const double crLHS137 = crLHS121 + r_C(0,1)*r_DN(2,0);
const double crLHS138 = crLHS127*r_DN(2,0);
const double crLHS139 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS140 = crLHS116 + crLHS48*crLHS9;
const double crLHS141 = bdf0*crLHS7;
const double crLHS142 = crLHS33*crLHS7;
const double crLHS143 = -crLHS21 + dp_th_dt*r_N[0];
const double crLHS144 = crLHS7*tau_t;
const double crLHS145 = crLHS144*crLHS80;
const double crLHS146 = -crLHS39 + dp_th_dt*r_N[1];
const double crLHS147 = crLHS7*crLHS93;
const double crLHS148 = crLHS147*crLHS16 + crLHS29*kappa - crLHS3*dp_th_dt*r_N[0]*r_N[1] + crLHS30*kappa;
const double crLHS149 = -crLHS55 + dp_th_dt*r_N[2];
const double crLHS150 = crLHS113*crLHS7;
const double crLHS151 = crLHS150*crLHS16 - crLHS3*dp_th_dt*r_N[0]*r_N[2] + crLHS47*kappa + crLHS48*kappa;
const double crLHS152 = r_DN(1,0)*r_DN(1,0);
const double crLHS153 = r_DN(1,1)*r_DN(1,1);
const double crLHS154 = r_N[1]*r_N[1];
const double crLHS155 = crLHS154*crLHS3;
const double crLHS156 = r_DN(1,0)*r_DN(2,0);
const double crLHS157 = r_DN(1,1)*r_DN(2,1);
const double crLHS158 = crLHS11*(crLHS156 + crLHS157);
const double crLHS159 = r_DN(2,0)*r_N[1];
const double crLHS160 = -crLHS94*r_N[2];
const double crLHS161 = r_DN(2,1)*r_N[1];
const double crLHS162 = -crLHS105*r_N[2];
const double crLHS163 = -crLHS107*r_N[2];
const double crLHS164 = crLHS9*r_DN(1,0);
const double crLHS165 = crLHS164*tau_c;
const double crLHS166 = crLHS37*crLHS74;
const double crLHS167 = crLHS77*r_N[1];
const double crLHS168 = crLHS61*r_N[1];
const double crLHS169 = crLHS79*tau_u;
const double crLHS170 = crLHS168*crLHS169;
const double crLHS171 = crLHS166*crLHS20 + crLHS167*crLHS76 - crLHS170*crLHS76 + crLHS19*crLHS96;
const double crLHS172 = r_DN(1,0)*tau_c;
const double crLHS173 = gauss_weight*(-crLHS168*crLHS64 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS174 = crLHS155*crLHS25 + crLHS166*crLHS38 + crLHS167*crLHS99 - crLHS170*crLHS99 + crLHS37*crLHS96;
const double crLHS175 = crLHS164*r_DN(1,1);
const double crLHS176 = r_DN(1,0)*r_N[2];
const double crLHS177 = crLHS115*crLHS36;
const double crLHS178 = crLHS156*crLHS9 + crLHS177;
const double crLHS179 = crLHS118*crLHS167 - crLHS118*crLHS170 + crLHS166*crLHS54 + crLHS53*crLHS96;
const double crLHS180 = crLHS164*r_DN(2,1);
const double crLHS181 = crLHS69*tau_c;
const double crLHS182 = crLHS9*r_DN(1,1);
const double crLHS183 = crLHS182*tau_c;
const double crLHS184 = r_DN(1,1)*tau_c;
const double crLHS185 = r_DN(1,1)*r_N[2];
const double crLHS186 = crLHS9*r_DN(2,0);
const double crLHS187 = crLHS186*r_DN(1,1);
const double crLHS188 = crLHS157*crLHS9 + crLHS177;
const double crLHS189 = crLHS144*crLHS79;
const double crLHS190 = crLHS168*crLHS189;
const double crLHS191 = crLHS150*crLHS36 + crLHS156*kappa + crLHS157*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[2];
const double crLHS192 = r_DN(2,0)*r_DN(2,0);
const double crLHS193 = r_DN(2,1)*r_DN(2,1);
const double crLHS194 = r_N[2]*r_N[2];
const double crLHS195 = crLHS194*crLHS3;
const double crLHS196 = crLHS53*crLHS74;
const double crLHS197 = crLHS77*r_N[2];
const double crLHS198 = crLHS61*r_N[2];
const double crLHS199 = crLHS169*crLHS198;
const double crLHS200 = crLHS115*crLHS19 + crLHS196*crLHS20 + crLHS197*crLHS76 - crLHS199*crLHS76;
const double crLHS201 = crLHS186*tau_c;
const double crLHS202 = r_DN(2,0)*tau_c;
const double crLHS203 = crLHS115*crLHS37 + crLHS196*crLHS38 + crLHS197*crLHS99 - crLHS199*crLHS99;
const double crLHS204 = gauss_weight*(-crLHS198*crLHS64 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS205 = crLHS115*crLHS53 + crLHS118*crLHS197 - crLHS118*crLHS199 + crLHS195*crLHS25 + crLHS196*crLHS54;
const double crLHS206 = crLHS186*r_DN(2,1);
const double crLHS207 = crLHS9*r_DN(2,1);
const double crLHS208 = crLHS207*tau_c;
const double crLHS209 = r_DN(2,1)*tau_c;
const double crLHS210 = crLHS189*crLHS198;
rLHS(0,0)+=crLHS11*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS10*(crLHS12 - crLHS13*crLHS15 + crLHS23*r_DN(0,0));
rLHS(0,2)+=crLHS10*(-crLHS15*crLHS24 + crLHS23*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS14*crLHS28;
rLHS(0,4)+=crLHS31;
rLHS(0,5)+=crLHS10*(crLHS32 + crLHS35 + crLHS40*r_DN(0,0));
rLHS(0,6)+=crLHS10*(crLHS40*r_DN(0,1) + crLHS41 + crLHS43);
rLHS(0,7)+=crLHS46;
rLHS(0,8)+=crLHS49;
rLHS(0,9)+=crLHS10*(crLHS50 + crLHS51 + crLHS56*r_DN(0,0));
rLHS(0,10)+=crLHS10*(crLHS56*r_DN(0,1) + crLHS57 + crLHS58);
rLHS(0,11)+=crLHS59;
rLHS(1,0)+=crLHS65*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS0*crLHS9 + crLHS66*r_DN(0,0) + crLHS68*r_DN(0,1) + crLHS69*crLHS71 + crLHS82);
rLHS(1,2)+=gauss_weight*(crLHS71*crLHS87 + crLHS83*r_DN(0,0) + crLHS85*r_DN(0,1) + crLHS86);
rLHS(1,3)+=-crLHS12*crLHS28*tau_c;
rLHS(1,4)+=gauss_weight*(crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS32*crLHS89 - crLHS88);
rLHS(1,5)+=gauss_weight*(crLHS100 + crLHS71*crLHS95 + crLHS90*r_DN(0,0) + crLHS92*r_DN(0,1) + crLHS98);
rLHS(1,6)+=gauss_weight*(crLHS101*r_DN(0,0) + crLHS103*r_DN(0,1) + crLHS104 + crLHS106*crLHS71);
rLHS(1,7)+=-crLHS107*crLHS108;
rLHS(1,8)+=gauss_weight*(-crLHS109 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS50*crLHS89);
rLHS(1,9)+=gauss_weight*(crLHS110*r_DN(0,0) + crLHS112*r_DN(0,1) + crLHS114*crLHS71 + crLHS117 + crLHS119);
rLHS(1,10)+=gauss_weight*(crLHS120*r_DN(0,0) + crLHS122*r_DN(0,1) + crLHS123 + crLHS124*crLHS71);
rLHS(1,11)+=-crLHS108*crLHS125;
rLHS(2,0)+=crLHS65*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS126*r_DN(0,1) + crLHS128*crLHS69 + crLHS68*r_DN(0,0) + crLHS86);
rLHS(2,2)+=gauss_weight*(crLHS1*crLHS9 + crLHS128*crLHS87 + crLHS129*r_DN(0,1) + crLHS82 + crLHS85*r_DN(0,0));
rLHS(2,3)+=-crLHS130*crLHS45;
rLHS(2,4)+=gauss_weight*(-crLHS131 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u - crLHS41*crLHS89);
rLHS(2,5)+=gauss_weight*(crLHS128*crLHS95 + crLHS132*r_DN(0,1) + crLHS133 + crLHS92*r_DN(0,0));
rLHS(2,6)+=gauss_weight*(crLHS100 + crLHS103*r_DN(0,0) + crLHS106*crLHS128 + crLHS134*r_DN(0,1) + crLHS135);
rLHS(2,7)+=-crLHS107*crLHS130;
rLHS(2,8)+=gauss_weight*(-crLHS136 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS57*crLHS89);
rLHS(2,9)+=gauss_weight*(crLHS112*r_DN(0,0) + crLHS114*crLHS128 + crLHS137*r_DN(0,1) + crLHS138);
rLHS(2,10)+=gauss_weight*(crLHS119 + crLHS122*r_DN(0,0) + crLHS124*crLHS128 + crLHS139*r_DN(0,1) + crLHS140);
rLHS(2,11)+=-crLHS125*crLHS130;
rLHS(3,0)+=0;
rLHS(3,1)+=0;
rLHS(3,2)+=0;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa + crLHS14*crLHS3*dp_th_dt - crLHS141*crLHS15 - crLHS142*crLHS19 - crLHS143*crLHS145 + crLHS143*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS143*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS143*crLHS26*dp_th_dt*r_N[0]*tau_t);
rLHS(3,4)+=0;
rLHS(3,5)+=0;
rLHS(3,6)+=0;
rLHS(3,7)+=-gauss_weight*(-crLHS142*crLHS37 - crLHS145*crLHS146 + crLHS146*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS146*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS146*crLHS26*dp_th_dt*r_N[0]*tau_t - crLHS148);
rLHS(3,8)+=0;
rLHS(3,9)+=0;
rLHS(3,10)+=0;
rLHS(3,11)+=-gauss_weight*(-crLHS142*crLHS53 - crLHS145*crLHS149 + crLHS149*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS149*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[0]*tau_t + crLHS149*crLHS26*dp_th_dt*r_N[0]*tau_t - crLHS151);
rLHS(4,0)+=crLHS31;
rLHS(4,1)+=crLHS10*(crLHS23*r_DN(1,0) + crLHS35 + crLHS88);
rLHS(4,2)+=crLHS10*(crLHS131 + crLHS23*r_DN(1,1) + crLHS43);
rLHS(4,3)+=crLHS46;
rLHS(4,4)+=crLHS11*(crLHS152 + crLHS153);
rLHS(4,5)+=crLHS10*(-crLHS13*crLHS155 + crLHS40*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS155*crLHS24 + crLHS40*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS154*crLHS28;
rLHS(4,8)+=crLHS158;
rLHS(4,9)+=crLHS10*(crLHS159 + crLHS160 + crLHS56*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS161 + crLHS162 + crLHS56*r_DN(1,1));
rLHS(4,11)+=crLHS163;
rLHS(5,0)+=gauss_weight*(crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS32 - crLHS88*crLHS89);
rLHS(5,1)+=gauss_weight*(crLHS165*crLHS69 + crLHS171 + crLHS66*r_DN(1,0) + crLHS68*r_DN(1,1) + crLHS98);
rLHS(5,2)+=gauss_weight*(crLHS133 + crLHS165*crLHS87 + crLHS83*r_DN(1,0) + crLHS85*r_DN(1,1));
rLHS(5,3)+=-crLHS172*crLHS45;
rLHS(5,4)+=crLHS173*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS152*crLHS9 + crLHS165*crLHS95 + crLHS174 + crLHS90*r_DN(1,0) + crLHS92*r_DN(1,1));
rLHS(5,6)+=gauss_weight*(crLHS101*r_DN(1,0) + crLHS103*r_DN(1,1) + crLHS106*crLHS165 + crLHS175);
rLHS(5,7)+=-crLHS107*crLHS172;
rLHS(5,8)+=gauss_weight*(-crLHS159*crLHS89 - crLHS176 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS110*r_DN(1,0) + crLHS112*r_DN(1,1) + crLHS114*crLHS165 + crLHS178 + crLHS179);
rLHS(5,10)+=gauss_weight*(crLHS120*r_DN(1,0) + crLHS122*r_DN(1,1) + crLHS124*crLHS165 + crLHS180);
rLHS(5,11)+=-crLHS125*crLHS172;
rLHS(6,0)+=gauss_weight*(-crLHS131*crLHS89 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS41);
rLHS(6,1)+=gauss_weight*(crLHS104 + crLHS126*r_DN(1,1) + crLHS181*crLHS182 + crLHS68*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS129*r_DN(1,1) + crLHS135 + crLHS171 + crLHS183*crLHS87 + crLHS85*r_DN(1,0));
rLHS(6,3)+=-crLHS184*crLHS45;
rLHS(6,4)+=crLHS173*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS132*r_DN(1,1) + crLHS175 + crLHS183*crLHS95 + crLHS92*r_DN(1,0));
rLHS(6,6)+=gauss_weight*(crLHS103*r_DN(1,0) + crLHS106*crLHS183 + crLHS134*r_DN(1,1) + crLHS153*crLHS9 + crLHS174);
rLHS(6,7)+=-crLHS107*crLHS184;
rLHS(6,8)+=gauss_weight*(-crLHS161*crLHS89 - crLHS185 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS112*r_DN(1,0) + crLHS114*crLHS183 + crLHS137*r_DN(1,1) + crLHS187);
rLHS(6,10)+=gauss_weight*(crLHS122*r_DN(1,0) + crLHS124*crLHS183 + crLHS139*r_DN(1,1) + crLHS179 + crLHS188);
rLHS(6,11)+=-crLHS125*crLHS184;
rLHS(7,0)+=0;
rLHS(7,1)+=0;
rLHS(7,2)+=0;
rLHS(7,3)+=-gauss_weight*(-crLHS143*crLHS190 + crLHS143*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS143*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS143*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS147*crLHS19 - crLHS148);
rLHS(7,4)+=0;
rLHS(7,5)+=0;
rLHS(7,6)+=0;
rLHS(7,7)+=-gauss_weight*(-crLHS141*crLHS155 - crLHS146*crLHS190 + crLHS146*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS146*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS146*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS147*crLHS37 - crLHS152*kappa - crLHS153*kappa + crLHS154*crLHS3*dp_th_dt);
rLHS(7,8)+=0;
rLHS(7,9)+=0;
rLHS(7,10)+=0;
rLHS(7,11)+=-gauss_weight*(-crLHS147*crLHS53 - crLHS149*crLHS190 + crLHS149*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS149*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[1]*tau_t + crLHS149*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS191);
rLHS(8,0)+=crLHS49;
rLHS(8,1)+=crLHS10*(crLHS109 + crLHS23*r_DN(2,0) + crLHS51);
rLHS(8,2)+=crLHS10*(crLHS136 + crLHS23*r_DN(2,1) + crLHS58);
rLHS(8,3)+=crLHS59;
rLHS(8,4)+=crLHS158;
rLHS(8,5)+=crLHS10*(crLHS160 + crLHS176 + crLHS40*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS162 + crLHS185 + crLHS40*r_DN(2,1));
rLHS(8,7)+=crLHS163;
rLHS(8,8)+=crLHS11*(crLHS192 + crLHS193);
rLHS(8,9)+=crLHS10*(-crLHS13*crLHS195 + crLHS56*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS195*crLHS24 + crLHS56*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS194*crLHS28;
rLHS(9,0)+=gauss_weight*(-crLHS109*crLHS89 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS50);
rLHS(9,1)+=gauss_weight*(crLHS117 + crLHS181*crLHS186 + crLHS200 + crLHS66*r_DN(2,0) + crLHS68*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS138 + crLHS201*crLHS87 + crLHS83*r_DN(2,0) + crLHS85*r_DN(2,1));
rLHS(9,3)+=-crLHS202*crLHS45;
rLHS(9,4)+=gauss_weight*(-crLHS159 - crLHS176*crLHS89 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS178 + crLHS201*crLHS95 + crLHS203 + crLHS90*r_DN(2,0) + crLHS92*r_DN(2,1));
rLHS(9,6)+=gauss_weight*(crLHS101*r_DN(2,0) + crLHS103*r_DN(2,1) + crLHS106*crLHS201 + crLHS187);
rLHS(9,7)+=-crLHS107*crLHS202;
rLHS(9,8)+=crLHS204*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS110*r_DN(2,0) + crLHS112*r_DN(2,1) + crLHS114*crLHS201 + crLHS192*crLHS9 + crLHS205);
rLHS(9,10)+=gauss_weight*(crLHS120*r_DN(2,0) + crLHS122*r_DN(2,1) + crLHS124*crLHS201 + crLHS206);
rLHS(9,11)+=-crLHS125*crLHS202;
rLHS(10,0)+=gauss_weight*(-crLHS136*crLHS89 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS57);
rLHS(10,1)+=gauss_weight*(crLHS123 + crLHS126*r_DN(2,1) + crLHS181*crLHS207 + crLHS68*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS129*r_DN(2,1) + crLHS140 + crLHS200 + crLHS208*crLHS87 + crLHS85*r_DN(2,0));
rLHS(10,3)+=-crLHS209*crLHS45;
rLHS(10,4)+=gauss_weight*(-crLHS161 - crLHS185*crLHS89 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS60*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS132*r_DN(2,1) + crLHS180 + crLHS208*crLHS95 + crLHS92*r_DN(2,0));
rLHS(10,6)+=gauss_weight*(crLHS103*r_DN(2,0) + crLHS106*crLHS208 + crLHS134*r_DN(2,1) + crLHS188 + crLHS203);
rLHS(10,7)+=-crLHS107*crLHS209;
rLHS(10,8)+=crLHS204*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS112*r_DN(2,0) + crLHS114*crLHS208 + crLHS137*r_DN(2,1) + crLHS206);
rLHS(10,10)+=gauss_weight*(crLHS122*r_DN(2,0) + crLHS124*crLHS208 + crLHS139*r_DN(2,1) + crLHS193*crLHS9 + crLHS205);
rLHS(10,11)+=-crLHS125*crLHS209;
rLHS(11,0)+=0;
rLHS(11,1)+=0;
rLHS(11,2)+=0;
rLHS(11,3)+=-gauss_weight*(-crLHS143*crLHS210 + crLHS143*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS143*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS143*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS150*crLHS19 - crLHS151);
rLHS(11,4)+=0;
rLHS(11,5)+=0;
rLHS(11,6)+=0;
rLHS(11,7)+=-gauss_weight*(-crLHS146*crLHS210 + crLHS146*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS146*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS146*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS150*crLHS37 - crLHS191);
rLHS(11,8)+=0;
rLHS(11,9)+=0;
rLHS(11,10)+=0;
rLHS(11,11)+=-gauss_weight*(-crLHS141*crLHS195 - crLHS149*crLHS210 + crLHS149*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS149*crLHS26*crLHS6*crLHS60*gamma*p_th*r_N[2]*tau_t + crLHS149*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS150*crLHS53 - crLHS192*kappa - crLHS193*kappa + crLHS194*crLHS3*dp_th_dt);

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
const double crLHS12 = r_DN(0,0)*r_N[0];
const double crLHS13 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS14 = r_N[0]*r_N[0];
const double crLHS15 = crLHS14*crLHS3;
const double crLHS16 = bdf0*r_N[0];
const double crLHS17 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crLHS18 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crLHS19 = crLHS17*r_DN(0,0) + crLHS18*r_DN(0,1);
const double crLHS20 = crLHS16 + crLHS19;
const double crLHS21 = crLHS20*crLHS7;
const double crLHS22 = crLHS3*crLHS4*tau_u;
const double crLHS23 = crLHS21*crLHS22;
const double crLHS24 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS25 = bdf0*crLHS8;
const double crLHS26 = 1.0/(crLHS2*crLHS2);
const double crLHS27 = crLHS26*gauss_weight;
const double crLHS28 = crLHS25*crLHS27;
const double crLHS29 = r_DN(0,0)*r_DN(1,0);
const double crLHS30 = r_DN(0,1)*r_DN(1,1);
const double crLHS31 = crLHS11*(crLHS29 + crLHS30);
const double crLHS32 = r_DN(1,0)*r_N[0];
const double crLHS33 = crLHS3*r_N[0];
const double crLHS34 = crLHS13*crLHS33;
const double crLHS35 = -crLHS34*r_N[1];
const double crLHS36 = bdf0*r_N[1];
const double crLHS37 = crLHS17*r_DN(1,0) + crLHS18*r_DN(1,1);
const double crLHS38 = crLHS36 + crLHS37;
const double crLHS39 = crLHS38*crLHS7;
const double crLHS40 = crLHS22*crLHS39;
const double crLHS41 = r_DN(1,1)*r_N[0];
const double crLHS42 = crLHS24*crLHS33;
const double crLHS43 = -crLHS42*r_N[1];
const double crLHS44 = crLHS27*crLHS8;
const double crLHS45 = crLHS16*crLHS44;
const double crLHS46 = -crLHS45*r_N[1];
const double crLHS47 = r_DN(0,0)*r_DN(2,0);
const double crLHS48 = r_DN(0,1)*r_DN(2,1);
const double crLHS49 = crLHS11*(crLHS47 + crLHS48);
const double crLHS50 = r_DN(2,0)*r_N[0];
const double crLHS51 = -crLHS34*r_N[2];
const double crLHS52 = bdf0*r_N[2];
const double crLHS53 = crLHS17*r_DN(2,0) + crLHS18*r_DN(2,1);
const double crLHS54 = crLHS52 + crLHS53;
const double crLHS55 = crLHS54*crLHS7;
const double crLHS56 = crLHS22*crLHS55;
const double crLHS57 = r_DN(2,1)*r_N[0];
const double crLHS58 = -crLHS42*r_N[2];
const double crLHS59 = -crLHS45*r_N[2];
const double crLHS60 = r_DN(0,0)*r_DN(3,0);
const double crLHS61 = r_DN(0,1)*r_DN(3,1);
const double crLHS62 = crLHS11*(crLHS60 + crLHS61);
const double crLHS63 = r_DN(3,0)*r_N[0];
const double crLHS64 = -crLHS34*r_N[3];
const double crLHS65 = bdf0*r_N[3];
const double crLHS66 = crLHS17*r_DN(3,0) + crLHS18*r_DN(3,1);
const double crLHS67 = crLHS65 + crLHS66;
const double crLHS68 = crLHS67*crLHS7;
const double crLHS69 = crLHS22*crLHS68;
const double crLHS70 = r_DN(3,1)*r_N[0];
const double crLHS71 = -crLHS42*r_N[3];
const double crLHS72 = -crLHS45*r_N[3];
const double crLHS73 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1);
const double crLHS74 = crLHS13*crLHS17 + crLHS18*crLHS24;
const double crLHS75 = crLHS74*r_N[0];
const double crLHS76 = crLHS26*tau_u;
const double crLHS77 = crLHS76*crLHS8;
const double crLHS78 = gauss_weight*(crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_u - crLHS75*crLHS77 - r_N[0]);
const double crLHS79 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS80 = r_C(0,2)*r_DN(0,0);
const double crLHS81 = crLHS80 + r_C(2,2)*r_DN(0,1);
const double crLHS82 = -crLHS34 + r_DN(0,0);
const double crLHS83 = crLHS9*r_DN(0,0);
const double crLHS84 = crLHS83*tau_c;
const double crLHS85 = crLHS33*crLHS8;
const double crLHS86 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS87 = crLHS76*crLHS86;
const double crLHS88 = crLHS19*crLHS87;
const double crLHS89 = crLHS20*crLHS86;
const double crLHS90 = crLHS26*crLHS73*tau_u;
const double crLHS91 = crLHS90*r_N[0];
const double crLHS92 = 1.0/(crLHS2*crLHS2*crLHS2);
const double crLHS93 = crLHS75*crLHS92;
const double crLHS94 = crLHS93*tau_u;
const double crLHS95 = crLHS15*crLHS25 + crLHS19*crLHS85 + crLHS20*crLHS88 + crLHS89*crLHS91 - crLHS89*crLHS94;
const double crLHS96 = crLHS80 + r_C(0,1)*r_DN(0,1);
const double crLHS97 = r_C(1,2)*r_DN(0,1);
const double crLHS98 = crLHS97 + r_C(2,2)*r_DN(0,0);
const double crLHS99 = crLHS83*r_DN(0,1);
const double crLHS100 = -crLHS42 + r_DN(0,1);
const double crLHS101 = r_DN(0,0)*r_N[1];
const double crLHS102 = crLHS74*crLHS77;
const double crLHS103 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS104 = r_C(0,2)*r_DN(1,0);
const double crLHS105 = crLHS104 + r_C(2,2)*r_DN(1,1);
const double crLHS106 = crLHS3*r_N[1];
const double crLHS107 = crLHS106*crLHS13;
const double crLHS108 = -crLHS107 + r_DN(1,0);
const double crLHS109 = crLHS106*crLHS8;
const double crLHS110 = crLHS109*crLHS16;
const double crLHS111 = crLHS110 + crLHS29*crLHS9;
const double crLHS112 = crLHS38*crLHS86;
const double crLHS113 = crLHS112*crLHS91 - crLHS112*crLHS94 + crLHS37*crLHS85 + crLHS38*crLHS88;
const double crLHS114 = crLHS104 + r_C(0,1)*r_DN(1,1);
const double crLHS115 = r_C(1,2)*r_DN(1,1);
const double crLHS116 = crLHS115 + r_C(2,2)*r_DN(1,0);
const double crLHS117 = crLHS83*r_DN(1,1);
const double crLHS118 = crLHS106*crLHS24;
const double crLHS119 = -crLHS118 + r_DN(1,1);
const double crLHS120 = crLHS36*crLHS44;
const double crLHS121 = r_DN(0,0)*tau_c;
const double crLHS122 = r_DN(0,0)*r_N[2];
const double crLHS123 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS124 = r_C(0,2)*r_DN(2,0);
const double crLHS125 = crLHS124 + r_C(2,2)*r_DN(2,1);
const double crLHS126 = crLHS3*r_N[2];
const double crLHS127 = crLHS126*crLHS13;
const double crLHS128 = -crLHS127 + r_DN(2,0);
const double crLHS129 = crLHS126*crLHS8;
const double crLHS130 = crLHS129*crLHS16;
const double crLHS131 = crLHS130 + crLHS47*crLHS9;
const double crLHS132 = crLHS54*crLHS86;
const double crLHS133 = crLHS132*crLHS91 - crLHS132*crLHS94 + crLHS53*crLHS85 + crLHS54*crLHS88;
const double crLHS134 = crLHS124 + r_C(0,1)*r_DN(2,1);
const double crLHS135 = r_C(1,2)*r_DN(2,1);
const double crLHS136 = crLHS135 + r_C(2,2)*r_DN(2,0);
const double crLHS137 = crLHS83*r_DN(2,1);
const double crLHS138 = crLHS126*crLHS24;
const double crLHS139 = -crLHS138 + r_DN(2,1);
const double crLHS140 = crLHS44*crLHS52;
const double crLHS141 = r_DN(0,0)*r_N[3];
const double crLHS142 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS143 = r_C(0,2)*r_DN(3,0);
const double crLHS144 = crLHS143 + r_C(2,2)*r_DN(3,1);
const double crLHS145 = crLHS3*r_N[3];
const double crLHS146 = -crLHS13*crLHS145 + r_DN(3,0);
const double crLHS147 = crLHS145*crLHS8;
const double crLHS148 = crLHS147*crLHS16;
const double crLHS149 = crLHS148 + crLHS60*crLHS9;
const double crLHS150 = crLHS67*crLHS86;
const double crLHS151 = crLHS150*crLHS91 - crLHS150*crLHS94 + crLHS66*crLHS85 + crLHS67*crLHS88;
const double crLHS152 = crLHS143 + r_C(0,1)*r_DN(3,1);
const double crLHS153 = r_C(1,2)*r_DN(3,1);
const double crLHS154 = crLHS153 + r_C(2,2)*r_DN(3,0);
const double crLHS155 = crLHS83*r_DN(3,1);
const double crLHS156 = -crLHS145*crLHS24 + r_DN(3,1);
const double crLHS157 = crLHS44*crLHS65;
const double crLHS158 = crLHS97 + r_C(0,1)*r_DN(0,0);
const double crLHS159 = crLHS9*r_DN(0,1);
const double crLHS160 = crLHS159*tau_c;
const double crLHS161 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS162 = r_DN(0,1)*tau_c;
const double crLHS163 = r_DN(0,1)*r_N[1];
const double crLHS164 = crLHS115 + r_C(0,1)*r_DN(1,0);
const double crLHS165 = crLHS159*r_DN(1,0);
const double crLHS166 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS167 = crLHS110 + crLHS30*crLHS9;
const double crLHS168 = r_DN(0,1)*r_N[2];
const double crLHS169 = crLHS135 + r_C(0,1)*r_DN(2,0);
const double crLHS170 = crLHS159*r_DN(2,0);
const double crLHS171 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS172 = crLHS130 + crLHS48*crLHS9;
const double crLHS173 = r_DN(0,1)*r_N[3];
const double crLHS174 = crLHS153 + r_C(0,1)*r_DN(3,0);
const double crLHS175 = crLHS159*r_DN(3,0);
const double crLHS176 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS177 = crLHS148 + crLHS61*crLHS9;
const double crLHS178 = bdf0*crLHS7;
const double crLHS179 = crLHS33*crLHS7;
const double crLHS180 = -crLHS21 + dp_th_dt*r_N[0];
const double crLHS181 = crLHS7*tau_t;
const double crLHS182 = crLHS181*crLHS93;
const double crLHS183 = -crLHS39 + dp_th_dt*r_N[1];
const double crLHS184 = crLHS106*crLHS7;
const double crLHS185 = crLHS16*crLHS184 + crLHS29*kappa - crLHS3*dp_th_dt*r_N[0]*r_N[1] + crLHS30*kappa;
const double crLHS186 = -crLHS55 + dp_th_dt*r_N[2];
const double crLHS187 = crLHS126*crLHS7;
const double crLHS188 = crLHS16*crLHS187 - crLHS3*dp_th_dt*r_N[0]*r_N[2] + crLHS47*kappa + crLHS48*kappa;
const double crLHS189 = -crLHS68 + dp_th_dt*r_N[3];
const double crLHS190 = crLHS145*crLHS7;
const double crLHS191 = crLHS16*crLHS190 - crLHS3*dp_th_dt*r_N[0]*r_N[3] + crLHS60*kappa + crLHS61*kappa;
const double crLHS192 = r_DN(1,0)*r_DN(1,0);
const double crLHS193 = r_DN(1,1)*r_DN(1,1);
const double crLHS194 = r_N[1]*r_N[1];
const double crLHS195 = crLHS194*crLHS3;
const double crLHS196 = r_DN(1,0)*r_DN(2,0);
const double crLHS197 = r_DN(1,1)*r_DN(2,1);
const double crLHS198 = crLHS11*(crLHS196 + crLHS197);
const double crLHS199 = r_DN(2,0)*r_N[1];
const double crLHS200 = -crLHS107*r_N[2];
const double crLHS201 = r_DN(2,1)*r_N[1];
const double crLHS202 = -crLHS118*r_N[2];
const double crLHS203 = -crLHS120*r_N[2];
const double crLHS204 = r_DN(1,0)*r_DN(3,0);
const double crLHS205 = r_DN(1,1)*r_DN(3,1);
const double crLHS206 = crLHS11*(crLHS204 + crLHS205);
const double crLHS207 = r_DN(3,0)*r_N[1];
const double crLHS208 = -crLHS107*r_N[3];
const double crLHS209 = r_DN(3,1)*r_N[1];
const double crLHS210 = -crLHS118*r_N[3];
const double crLHS211 = -crLHS120*r_N[3];
const double crLHS212 = crLHS9*r_DN(1,0);
const double crLHS213 = crLHS212*tau_c;
const double crLHS214 = crLHS37*crLHS87;
const double crLHS215 = crLHS90*r_N[1];
const double crLHS216 = crLHS74*r_N[1];
const double crLHS217 = crLHS92*tau_u;
const double crLHS218 = crLHS216*crLHS217;
const double crLHS219 = crLHS109*crLHS19 + crLHS20*crLHS214 + crLHS215*crLHS89 - crLHS218*crLHS89;
const double crLHS220 = r_DN(1,0)*tau_c;
const double crLHS221 = gauss_weight*(-crLHS216*crLHS77 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_u - r_N[1]);
const double crLHS222 = crLHS109*crLHS37 + crLHS112*crLHS215 - crLHS112*crLHS218 + crLHS195*crLHS25 + crLHS214*crLHS38;
const double crLHS223 = crLHS212*r_DN(1,1);
const double crLHS224 = r_DN(1,0)*r_N[2];
const double crLHS225 = crLHS129*crLHS36;
const double crLHS226 = crLHS196*crLHS9 + crLHS225;
const double crLHS227 = crLHS109*crLHS53 + crLHS132*crLHS215 - crLHS132*crLHS218 + crLHS214*crLHS54;
const double crLHS228 = crLHS212*r_DN(2,1);
const double crLHS229 = r_DN(1,0)*r_N[3];
const double crLHS230 = crLHS147*crLHS36;
const double crLHS231 = crLHS204*crLHS9 + crLHS230;
const double crLHS232 = crLHS109*crLHS66 + crLHS150*crLHS215 - crLHS150*crLHS218 + crLHS214*crLHS67;
const double crLHS233 = crLHS212*r_DN(3,1);
const double crLHS234 = crLHS9*r_DN(1,1);
const double crLHS235 = crLHS234*tau_c;
const double crLHS236 = r_DN(1,1)*tau_c;
const double crLHS237 = r_DN(1,1)*r_N[2];
const double crLHS238 = crLHS234*r_DN(2,0);
const double crLHS239 = crLHS197*crLHS9 + crLHS225;
const double crLHS240 = r_DN(1,1)*r_N[3];
const double crLHS241 = crLHS234*r_DN(3,0);
const double crLHS242 = crLHS205*crLHS9 + crLHS230;
const double crLHS243 = crLHS181*crLHS92;
const double crLHS244 = crLHS216*crLHS243;
const double crLHS245 = crLHS187*crLHS36 + crLHS196*kappa + crLHS197*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[2];
const double crLHS246 = crLHS190*crLHS36 + crLHS204*kappa + crLHS205*kappa - crLHS3*dp_th_dt*r_N[1]*r_N[3];
const double crLHS247 = r_DN(2,0)*r_DN(2,0);
const double crLHS248 = r_DN(2,1)*r_DN(2,1);
const double crLHS249 = r_N[2]*r_N[2];
const double crLHS250 = crLHS249*crLHS3;
const double crLHS251 = r_DN(2,0)*r_DN(3,0);
const double crLHS252 = r_DN(2,1)*r_DN(3,1);
const double crLHS253 = crLHS11*(crLHS251 + crLHS252);
const double crLHS254 = r_DN(3,0)*r_N[2];
const double crLHS255 = -crLHS127*r_N[3];
const double crLHS256 = r_DN(3,1)*r_N[2];
const double crLHS257 = -crLHS138*r_N[3];
const double crLHS258 = -crLHS140*r_N[3];
const double crLHS259 = crLHS9*r_DN(2,0);
const double crLHS260 = crLHS259*tau_c;
const double crLHS261 = crLHS53*crLHS87;
const double crLHS262 = crLHS90*r_N[2];
const double crLHS263 = crLHS74*r_N[2];
const double crLHS264 = crLHS217*crLHS263;
const double crLHS265 = crLHS129*crLHS19 + crLHS20*crLHS261 + crLHS262*crLHS89 - crLHS264*crLHS89;
const double crLHS266 = r_DN(2,0)*tau_c;
const double crLHS267 = crLHS112*crLHS262 - crLHS112*crLHS264 + crLHS129*crLHS37 + crLHS261*crLHS38;
const double crLHS268 = gauss_weight*(-crLHS263*crLHS77 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_u - r_N[2]);
const double crLHS269 = crLHS129*crLHS53 + crLHS132*crLHS262 - crLHS132*crLHS264 + crLHS25*crLHS250 + crLHS261*crLHS54;
const double crLHS270 = crLHS259*r_DN(2,1);
const double crLHS271 = r_DN(2,0)*r_N[3];
const double crLHS272 = crLHS147*crLHS52;
const double crLHS273 = crLHS251*crLHS9 + crLHS272;
const double crLHS274 = crLHS129*crLHS66 + crLHS150*crLHS262 - crLHS150*crLHS264 + crLHS261*crLHS67;
const double crLHS275 = crLHS259*r_DN(3,1);
const double crLHS276 = crLHS82*tau_c;
const double crLHS277 = crLHS9*r_DN(2,1);
const double crLHS278 = crLHS277*tau_c;
const double crLHS279 = r_DN(2,1)*tau_c;
const double crLHS280 = r_DN(2,1)*r_N[3];
const double crLHS281 = crLHS9*r_DN(3,0);
const double crLHS282 = crLHS281*r_DN(2,1);
const double crLHS283 = crLHS252*crLHS9 + crLHS272;
const double crLHS284 = crLHS243*crLHS263;
const double crLHS285 = crLHS190*crLHS52 + crLHS251*kappa + crLHS252*kappa - crLHS3*dp_th_dt*r_N[2]*r_N[3];
const double crLHS286 = r_DN(3,0)*r_DN(3,0);
const double crLHS287 = r_DN(3,1)*r_DN(3,1);
const double crLHS288 = r_N[3]*r_N[3];
const double crLHS289 = crLHS288*crLHS3;
const double crLHS290 = crLHS66*crLHS87;
const double crLHS291 = crLHS90*r_N[3];
const double crLHS292 = crLHS74*r_N[3];
const double crLHS293 = crLHS217*crLHS292;
const double crLHS294 = crLHS147*crLHS19 + crLHS20*crLHS290 + crLHS291*crLHS89 - crLHS293*crLHS89;
const double crLHS295 = crLHS281*tau_c;
const double crLHS296 = r_DN(3,0)*tau_c;
const double crLHS297 = crLHS112*crLHS291 - crLHS112*crLHS293 + crLHS147*crLHS37 + crLHS290*crLHS38;
const double crLHS298 = crLHS132*crLHS291 - crLHS132*crLHS293 + crLHS147*crLHS53 + crLHS290*crLHS54;
const double crLHS299 = gauss_weight*(-crLHS292*crLHS77 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_u - r_N[3]);
const double crLHS300 = crLHS147*crLHS66 + crLHS150*crLHS291 - crLHS150*crLHS293 + crLHS25*crLHS289 + crLHS290*crLHS67;
const double crLHS301 = crLHS281*r_DN(3,1);
const double crLHS302 = crLHS9*r_DN(3,1);
const double crLHS303 = crLHS302*tau_c;
const double crLHS304 = r_DN(3,1)*tau_c;
const double crLHS305 = crLHS243*crLHS292;
rLHS(0,0)+=crLHS11*(crLHS0 + crLHS1);
rLHS(0,1)+=crLHS10*(crLHS12 - crLHS13*crLHS15 + crLHS23*r_DN(0,0));
rLHS(0,2)+=crLHS10*(-crLHS15*crLHS24 + crLHS23*r_DN(0,1) + r_DN(0,1)*r_N[0]);
rLHS(0,3)+=-crLHS14*crLHS28;
rLHS(0,4)+=crLHS31;
rLHS(0,5)+=crLHS10*(crLHS32 + crLHS35 + crLHS40*r_DN(0,0));
rLHS(0,6)+=crLHS10*(crLHS40*r_DN(0,1) + crLHS41 + crLHS43);
rLHS(0,7)+=crLHS46;
rLHS(0,8)+=crLHS49;
rLHS(0,9)+=crLHS10*(crLHS50 + crLHS51 + crLHS56*r_DN(0,0));
rLHS(0,10)+=crLHS10*(crLHS56*r_DN(0,1) + crLHS57 + crLHS58);
rLHS(0,11)+=crLHS59;
rLHS(0,12)+=crLHS62;
rLHS(0,13)+=crLHS10*(crLHS63 + crLHS64 + crLHS69*r_DN(0,0));
rLHS(0,14)+=crLHS10*(crLHS69*r_DN(0,1) + crLHS70 + crLHS71);
rLHS(0,15)+=crLHS72;
rLHS(1,0)+=crLHS78*r_DN(0,0);
rLHS(1,1)+=gauss_weight*(crLHS0*crLHS9 + crLHS79*r_DN(0,0) + crLHS81*r_DN(0,1) + crLHS82*crLHS84 + crLHS95);
rLHS(1,2)+=gauss_weight*(crLHS100*crLHS84 + crLHS96*r_DN(0,0) + crLHS98*r_DN(0,1) + crLHS99);
rLHS(1,3)+=-crLHS12*crLHS28*tau_c;
rLHS(1,4)+=gauss_weight*(-crLHS101 - crLHS102*crLHS32 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u);
rLHS(1,5)+=gauss_weight*(crLHS103*r_DN(0,0) + crLHS105*r_DN(0,1) + crLHS108*crLHS84 + crLHS111 + crLHS113);
rLHS(1,6)+=gauss_weight*(crLHS114*r_DN(0,0) + crLHS116*r_DN(0,1) + crLHS117 + crLHS119*crLHS84);
rLHS(1,7)+=-crLHS120*crLHS121;
rLHS(1,8)+=gauss_weight*(-crLHS102*crLHS50 - crLHS122 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u);
rLHS(1,9)+=gauss_weight*(crLHS123*r_DN(0,0) + crLHS125*r_DN(0,1) + crLHS128*crLHS84 + crLHS131 + crLHS133);
rLHS(1,10)+=gauss_weight*(crLHS134*r_DN(0,0) + crLHS136*r_DN(0,1) + crLHS137 + crLHS139*crLHS84);
rLHS(1,11)+=-crLHS121*crLHS140;
rLHS(1,12)+=gauss_weight*(-crLHS102*crLHS63 - crLHS141 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[0]*tau_u);
rLHS(1,13)+=gauss_weight*(crLHS142*r_DN(0,0) + crLHS144*r_DN(0,1) + crLHS146*crLHS84 + crLHS149 + crLHS151);
rLHS(1,14)+=gauss_weight*(crLHS152*r_DN(0,0) + crLHS154*r_DN(0,1) + crLHS155 + crLHS156*crLHS84);
rLHS(1,15)+=-crLHS121*crLHS157;
rLHS(2,0)+=crLHS78*r_DN(0,1);
rLHS(2,1)+=gauss_weight*(crLHS158*r_DN(0,1) + crLHS160*crLHS82 + crLHS81*r_DN(0,0) + crLHS99);
rLHS(2,2)+=gauss_weight*(crLHS1*crLHS9 + crLHS100*crLHS160 + crLHS161*r_DN(0,1) + crLHS95 + crLHS98*r_DN(0,0));
rLHS(2,3)+=-crLHS162*crLHS45;
rLHS(2,4)+=gauss_weight*(-crLHS102*crLHS41 - crLHS163 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u);
rLHS(2,5)+=gauss_weight*(crLHS105*r_DN(0,0) + crLHS108*crLHS160 + crLHS164*r_DN(0,1) + crLHS165);
rLHS(2,6)+=gauss_weight*(crLHS113 + crLHS116*r_DN(0,0) + crLHS119*crLHS160 + crLHS166*r_DN(0,1) + crLHS167);
rLHS(2,7)+=-crLHS120*crLHS162;
rLHS(2,8)+=gauss_weight*(-crLHS102*crLHS57 - crLHS168 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u);
rLHS(2,9)+=gauss_weight*(crLHS125*r_DN(0,0) + crLHS128*crLHS160 + crLHS169*r_DN(0,1) + crLHS170);
rLHS(2,10)+=gauss_weight*(crLHS133 + crLHS136*r_DN(0,0) + crLHS139*crLHS160 + crLHS171*r_DN(0,1) + crLHS172);
rLHS(2,11)+=-crLHS140*crLHS162;
rLHS(2,12)+=gauss_weight*(-crLHS102*crLHS70 - crLHS173 + crLHS19*crLHS3*crLHS4*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[0]*tau_u);
rLHS(2,13)+=gauss_weight*(crLHS144*r_DN(0,0) + crLHS146*crLHS160 + crLHS174*r_DN(0,1) + crLHS175);
rLHS(2,14)+=gauss_weight*(crLHS151 + crLHS154*r_DN(0,0) + crLHS156*crLHS160 + crLHS176*r_DN(0,1) + crLHS177);
rLHS(2,15)+=-crLHS157*crLHS162;
rLHS(3,0)+=0;
rLHS(3,1)+=0;
rLHS(3,2)+=0;
rLHS(3,3)+=-gauss_weight*(-crLHS0*kappa - crLHS1*kappa + crLHS14*crLHS3*dp_th_dt - crLHS15*crLHS178 - crLHS179*crLHS19 - crLHS180*crLHS182 + crLHS180*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS180*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS180*crLHS26*dp_th_dt*r_N[0]*tau_t);
rLHS(3,4)+=0;
rLHS(3,5)+=0;
rLHS(3,6)+=0;
rLHS(3,7)+=-gauss_weight*(-crLHS179*crLHS37 - crLHS182*crLHS183 + crLHS183*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS183*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS183*crLHS26*dp_th_dt*r_N[0]*tau_t - crLHS185);
rLHS(3,8)+=0;
rLHS(3,9)+=0;
rLHS(3,10)+=0;
rLHS(3,11)+=-gauss_weight*(-crLHS179*crLHS53 - crLHS182*crLHS186 + crLHS186*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS186*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS186*crLHS26*dp_th_dt*r_N[0]*tau_t - crLHS188);
rLHS(3,12)+=0;
rLHS(3,13)+=0;
rLHS(3,14)+=0;
rLHS(3,15)+=-gauss_weight*(-crLHS179*crLHS66 - crLHS182*crLHS189 + crLHS189*crLHS19*crLHS26*crLHS6*gamma*p_th*tau_t + crLHS189*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[0]*tau_t + crLHS189*crLHS26*dp_th_dt*r_N[0]*tau_t - crLHS191);
rLHS(4,0)+=crLHS31;
rLHS(4,1)+=crLHS10*(crLHS101 + crLHS23*r_DN(1,0) + crLHS35);
rLHS(4,2)+=crLHS10*(crLHS163 + crLHS23*r_DN(1,1) + crLHS43);
rLHS(4,3)+=crLHS46;
rLHS(4,4)+=crLHS11*(crLHS192 + crLHS193);
rLHS(4,5)+=crLHS10*(-crLHS13*crLHS195 + crLHS40*r_DN(1,0) + r_DN(1,0)*r_N[1]);
rLHS(4,6)+=crLHS10*(-crLHS195*crLHS24 + crLHS40*r_DN(1,1) + r_DN(1,1)*r_N[1]);
rLHS(4,7)+=-crLHS194*crLHS28;
rLHS(4,8)+=crLHS198;
rLHS(4,9)+=crLHS10*(crLHS199 + crLHS200 + crLHS56*r_DN(1,0));
rLHS(4,10)+=crLHS10*(crLHS201 + crLHS202 + crLHS56*r_DN(1,1));
rLHS(4,11)+=crLHS203;
rLHS(4,12)+=crLHS206;
rLHS(4,13)+=crLHS10*(crLHS207 + crLHS208 + crLHS69*r_DN(1,0));
rLHS(4,14)+=crLHS10*(crLHS209 + crLHS210 + crLHS69*r_DN(1,1));
rLHS(4,15)+=crLHS211;
rLHS(5,0)+=gauss_weight*(-crLHS101*crLHS102 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS32);
rLHS(5,1)+=gauss_weight*(crLHS111 + crLHS213*crLHS82 + crLHS219 + crLHS79*r_DN(1,0) + crLHS81*r_DN(1,1));
rLHS(5,2)+=gauss_weight*(crLHS100*crLHS213 + crLHS165 + crLHS96*r_DN(1,0) + crLHS98*r_DN(1,1));
rLHS(5,3)+=-crLHS220*crLHS45;
rLHS(5,4)+=crLHS221*r_DN(1,0);
rLHS(5,5)+=gauss_weight*(crLHS103*r_DN(1,0) + crLHS105*r_DN(1,1) + crLHS108*crLHS213 + crLHS192*crLHS9 + crLHS222);
rLHS(5,6)+=gauss_weight*(crLHS114*r_DN(1,0) + crLHS116*r_DN(1,1) + crLHS119*crLHS213 + crLHS223);
rLHS(5,7)+=-crLHS120*crLHS220;
rLHS(5,8)+=gauss_weight*(-crLHS102*crLHS199 - crLHS224 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
rLHS(5,9)+=gauss_weight*(crLHS123*r_DN(1,0) + crLHS125*r_DN(1,1) + crLHS128*crLHS213 + crLHS226 + crLHS227);
rLHS(5,10)+=gauss_weight*(crLHS134*r_DN(1,0) + crLHS136*r_DN(1,1) + crLHS139*crLHS213 + crLHS228);
rLHS(5,11)+=-crLHS140*crLHS220;
rLHS(5,12)+=gauss_weight*(-crLHS102*crLHS207 - crLHS229 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[1]*tau_u);
rLHS(5,13)+=gauss_weight*(crLHS142*r_DN(1,0) + crLHS144*r_DN(1,1) + crLHS146*crLHS213 + crLHS231 + crLHS232);
rLHS(5,14)+=gauss_weight*(crLHS152*r_DN(1,0) + crLHS154*r_DN(1,1) + crLHS156*crLHS213 + crLHS233);
rLHS(5,15)+=-crLHS157*crLHS220;
rLHS(6,0)+=gauss_weight*(-crLHS102*crLHS163 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u - crLHS41);
rLHS(6,1)+=gauss_weight*(crLHS117 + crLHS158*r_DN(1,1) + crLHS235*crLHS82 + crLHS81*r_DN(1,0));
rLHS(6,2)+=gauss_weight*(crLHS100*crLHS235 + crLHS161*r_DN(1,1) + crLHS167 + crLHS219 + crLHS98*r_DN(1,0));
rLHS(6,3)+=-crLHS236*crLHS45;
rLHS(6,4)+=crLHS221*r_DN(1,1);
rLHS(6,5)+=gauss_weight*(crLHS105*r_DN(1,0) + crLHS108*crLHS235 + crLHS164*r_DN(1,1) + crLHS223);
rLHS(6,6)+=gauss_weight*(crLHS116*r_DN(1,0) + crLHS119*crLHS235 + crLHS166*r_DN(1,1) + crLHS193*crLHS9 + crLHS222);
rLHS(6,7)+=-crLHS120*crLHS236;
rLHS(6,8)+=gauss_weight*(-crLHS102*crLHS201 - crLHS237 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
rLHS(6,9)+=gauss_weight*(crLHS125*r_DN(1,0) + crLHS128*crLHS235 + crLHS169*r_DN(1,1) + crLHS238);
rLHS(6,10)+=gauss_weight*(crLHS136*r_DN(1,0) + crLHS139*crLHS235 + crLHS171*r_DN(1,1) + crLHS227 + crLHS239);
rLHS(6,11)+=-crLHS140*crLHS236;
rLHS(6,12)+=gauss_weight*(-crLHS102*crLHS209 - crLHS240 + crLHS3*crLHS37*crLHS4*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[1]*tau_u);
rLHS(6,13)+=gauss_weight*(crLHS144*r_DN(1,0) + crLHS146*crLHS235 + crLHS174*r_DN(1,1) + crLHS241);
rLHS(6,14)+=gauss_weight*(crLHS154*r_DN(1,0) + crLHS156*crLHS235 + crLHS176*r_DN(1,1) + crLHS232 + crLHS242);
rLHS(6,15)+=-crLHS157*crLHS236;
rLHS(7,0)+=0;
rLHS(7,1)+=0;
rLHS(7,2)+=0;
rLHS(7,3)+=-gauss_weight*(-crLHS180*crLHS244 + crLHS180*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS180*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS180*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS184*crLHS19 - crLHS185);
rLHS(7,4)+=0;
rLHS(7,5)+=0;
rLHS(7,6)+=0;
rLHS(7,7)+=-gauss_weight*(-crLHS178*crLHS195 - crLHS183*crLHS244 + crLHS183*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS183*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS183*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS184*crLHS37 - crLHS192*kappa - crLHS193*kappa + crLHS194*crLHS3*dp_th_dt);
rLHS(7,8)+=0;
rLHS(7,9)+=0;
rLHS(7,10)+=0;
rLHS(7,11)+=-gauss_weight*(-crLHS184*crLHS53 - crLHS186*crLHS244 + crLHS186*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS186*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS186*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS245);
rLHS(7,12)+=0;
rLHS(7,13)+=0;
rLHS(7,14)+=0;
rLHS(7,15)+=-gauss_weight*(-crLHS184*crLHS66 - crLHS189*crLHS244 + crLHS189*crLHS26*crLHS37*crLHS6*gamma*p_th*tau_t + crLHS189*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[1]*tau_t + crLHS189*crLHS26*dp_th_dt*r_N[1]*tau_t - crLHS246);
rLHS(8,0)+=crLHS49;
rLHS(8,1)+=crLHS10*(crLHS122 + crLHS23*r_DN(2,0) + crLHS51);
rLHS(8,2)+=crLHS10*(crLHS168 + crLHS23*r_DN(2,1) + crLHS58);
rLHS(8,3)+=crLHS59;
rLHS(8,4)+=crLHS198;
rLHS(8,5)+=crLHS10*(crLHS200 + crLHS224 + crLHS40*r_DN(2,0));
rLHS(8,6)+=crLHS10*(crLHS202 + crLHS237 + crLHS40*r_DN(2,1));
rLHS(8,7)+=crLHS203;
rLHS(8,8)+=crLHS11*(crLHS247 + crLHS248);
rLHS(8,9)+=crLHS10*(-crLHS13*crLHS250 + crLHS56*r_DN(2,0) + r_DN(2,0)*r_N[2]);
rLHS(8,10)+=crLHS10*(-crLHS24*crLHS250 + crLHS56*r_DN(2,1) + r_DN(2,1)*r_N[2]);
rLHS(8,11)+=-crLHS249*crLHS28;
rLHS(8,12)+=crLHS253;
rLHS(8,13)+=crLHS10*(crLHS254 + crLHS255 + crLHS69*r_DN(2,0));
rLHS(8,14)+=crLHS10*(crLHS256 + crLHS257 + crLHS69*r_DN(2,1));
rLHS(8,15)+=crLHS258;
rLHS(9,0)+=gauss_weight*(-crLHS102*crLHS122 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS50);
rLHS(9,1)+=gauss_weight*(crLHS131 + crLHS260*crLHS82 + crLHS265 + crLHS79*r_DN(2,0) + crLHS81*r_DN(2,1));
rLHS(9,2)+=gauss_weight*(crLHS100*crLHS260 + crLHS170 + crLHS96*r_DN(2,0) + crLHS98*r_DN(2,1));
rLHS(9,3)+=-crLHS266*crLHS45;
rLHS(9,4)+=gauss_weight*(-crLHS102*crLHS224 - crLHS199 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
rLHS(9,5)+=gauss_weight*(crLHS103*r_DN(2,0) + crLHS105*r_DN(2,1) + crLHS108*crLHS260 + crLHS226 + crLHS267);
rLHS(9,6)+=gauss_weight*(crLHS114*r_DN(2,0) + crLHS116*r_DN(2,1) + crLHS119*crLHS260 + crLHS238);
rLHS(9,7)+=-crLHS120*crLHS266;
rLHS(9,8)+=crLHS268*r_DN(2,0);
rLHS(9,9)+=gauss_weight*(crLHS123*r_DN(2,0) + crLHS125*r_DN(2,1) + crLHS128*crLHS260 + crLHS247*crLHS9 + crLHS269);
rLHS(9,10)+=gauss_weight*(crLHS134*r_DN(2,0) + crLHS136*r_DN(2,1) + crLHS139*crLHS260 + crLHS270);
rLHS(9,11)+=-crLHS140*crLHS266;
rLHS(9,12)+=gauss_weight*(-crLHS102*crLHS254 - crLHS271 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(3,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,0)*r_N[2]*tau_u);
rLHS(9,13)+=gauss_weight*(crLHS142*r_DN(2,0) + crLHS144*r_DN(2,1) + crLHS146*crLHS260 + crLHS273 + crLHS274);
rLHS(9,14)+=gauss_weight*(crLHS152*r_DN(2,0) + crLHS154*r_DN(2,1) + crLHS156*crLHS260 + crLHS275);
rLHS(9,15)+=-crLHS157*crLHS266;
rLHS(10,0)+=gauss_weight*(-crLHS102*crLHS168 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS57);
rLHS(10,1)+=gauss_weight*(crLHS137 + crLHS158*r_DN(2,1) + crLHS276*crLHS277 + crLHS81*r_DN(2,0));
rLHS(10,2)+=gauss_weight*(crLHS100*crLHS278 + crLHS161*r_DN(2,1) + crLHS172 + crLHS265 + crLHS98*r_DN(2,0));
rLHS(10,3)+=-crLHS279*crLHS45;
rLHS(10,4)+=gauss_weight*(-crLHS102*crLHS237 - crLHS201 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
rLHS(10,5)+=gauss_weight*(crLHS105*r_DN(2,0) + crLHS108*crLHS278 + crLHS164*r_DN(2,1) + crLHS228);
rLHS(10,6)+=gauss_weight*(crLHS116*r_DN(2,0) + crLHS119*crLHS278 + crLHS166*r_DN(2,1) + crLHS239 + crLHS267);
rLHS(10,7)+=-crLHS120*crLHS279;
rLHS(10,8)+=crLHS268*r_DN(2,1);
rLHS(10,9)+=gauss_weight*(crLHS125*r_DN(2,0) + crLHS128*crLHS278 + crLHS169*r_DN(2,1) + crLHS270);
rLHS(10,10)+=gauss_weight*(crLHS136*r_DN(2,0) + crLHS139*crLHS278 + crLHS171*r_DN(2,1) + crLHS248*crLHS9 + crLHS269);
rLHS(10,11)+=-crLHS140*crLHS279;
rLHS(10,12)+=gauss_weight*(-crLHS102*crLHS256 - crLHS280 + crLHS3*crLHS4*crLHS53*crLHS6*gamma*p_th*r_DN(3,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(3,1)*r_N[2]*tau_u);
rLHS(10,13)+=gauss_weight*(crLHS144*r_DN(2,0) + crLHS146*crLHS278 + crLHS174*r_DN(2,1) + crLHS282);
rLHS(10,14)+=gauss_weight*(crLHS154*r_DN(2,0) + crLHS156*crLHS278 + crLHS176*r_DN(2,1) + crLHS274 + crLHS283);
rLHS(10,15)+=-crLHS157*crLHS279;
rLHS(11,0)+=0;
rLHS(11,1)+=0;
rLHS(11,2)+=0;
rLHS(11,3)+=-gauss_weight*(crLHS180*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS180*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS180*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS180*crLHS284 - crLHS187*crLHS19 - crLHS188);
rLHS(11,4)+=0;
rLHS(11,5)+=0;
rLHS(11,6)+=0;
rLHS(11,7)+=-gauss_weight*(crLHS183*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS183*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS183*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS183*crLHS284 - crLHS187*crLHS37 - crLHS245);
rLHS(11,8)+=0;
rLHS(11,9)+=0;
rLHS(11,10)+=0;
rLHS(11,11)+=-gauss_weight*(-crLHS178*crLHS250 + crLHS186*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS186*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS186*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS186*crLHS284 - crLHS187*crLHS53 - crLHS247*kappa - crLHS248*kappa + crLHS249*crLHS3*dp_th_dt);
rLHS(11,12)+=0;
rLHS(11,13)+=0;
rLHS(11,14)+=0;
rLHS(11,15)+=-gauss_weight*(-crLHS187*crLHS66 + crLHS189*crLHS26*crLHS53*crLHS6*gamma*p_th*tau_t + crLHS189*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[2]*tau_t + crLHS189*crLHS26*dp_th_dt*r_N[2]*tau_t - crLHS189*crLHS284 - crLHS285);
rLHS(12,0)+=crLHS62;
rLHS(12,1)+=crLHS10*(crLHS141 + crLHS23*r_DN(3,0) + crLHS64);
rLHS(12,2)+=crLHS10*(crLHS173 + crLHS23*r_DN(3,1) + crLHS71);
rLHS(12,3)+=crLHS72;
rLHS(12,4)+=crLHS206;
rLHS(12,5)+=crLHS10*(crLHS208 + crLHS229 + crLHS40*r_DN(3,0));
rLHS(12,6)+=crLHS10*(crLHS210 + crLHS240 + crLHS40*r_DN(3,1));
rLHS(12,7)+=crLHS211;
rLHS(12,8)+=crLHS253;
rLHS(12,9)+=crLHS10*(crLHS255 + crLHS271 + crLHS56*r_DN(3,0));
rLHS(12,10)+=crLHS10*(crLHS257 + crLHS280 + crLHS56*r_DN(3,1));
rLHS(12,11)+=crLHS258;
rLHS(12,12)+=crLHS11*(crLHS286 + crLHS287);
rLHS(12,13)+=crLHS10*(-crLHS13*crLHS289 + crLHS69*r_DN(3,0) + r_DN(3,0)*r_N[3]);
rLHS(12,14)+=crLHS10*(-crLHS24*crLHS289 + crLHS69*r_DN(3,1) + r_DN(3,1)*r_N[3]);
rLHS(12,15)+=-crLHS28*crLHS288;
rLHS(13,0)+=gauss_weight*(-crLHS102*crLHS141 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(0,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,0)*r_N[3]*tau_u - crLHS63);
rLHS(13,1)+=gauss_weight*(crLHS149 + crLHS276*crLHS281 + crLHS294 + crLHS79*r_DN(3,0) + crLHS81*r_DN(3,1));
rLHS(13,2)+=gauss_weight*(crLHS100*crLHS295 + crLHS175 + crLHS96*r_DN(3,0) + crLHS98*r_DN(3,1));
rLHS(13,3)+=-crLHS296*crLHS45;
rLHS(13,4)+=gauss_weight*(-crLHS102*crLHS229 - crLHS207 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(1,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,0)*r_N[3]*tau_u);
rLHS(13,5)+=gauss_weight*(crLHS103*r_DN(3,0) + crLHS105*r_DN(3,1) + crLHS108*crLHS295 + crLHS231 + crLHS297);
rLHS(13,6)+=gauss_weight*(crLHS114*r_DN(3,0) + crLHS116*r_DN(3,1) + crLHS119*crLHS295 + crLHS241);
rLHS(13,7)+=-crLHS120*crLHS296;
rLHS(13,8)+=gauss_weight*(-crLHS102*crLHS271 - crLHS254 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(2,0)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,0)*r_N[3]*tau_u);
rLHS(13,9)+=gauss_weight*(crLHS123*r_DN(3,0) + crLHS125*r_DN(3,1) + crLHS128*crLHS295 + crLHS273 + crLHS298);
rLHS(13,10)+=gauss_weight*(crLHS134*r_DN(3,0) + crLHS136*r_DN(3,1) + crLHS139*crLHS295 + crLHS282);
rLHS(13,11)+=-crLHS140*crLHS296;
rLHS(13,12)+=crLHS299*r_DN(3,0);
rLHS(13,13)+=gauss_weight*(crLHS142*r_DN(3,0) + crLHS144*r_DN(3,1) + crLHS146*crLHS295 + crLHS286*crLHS9 + crLHS300);
rLHS(13,14)+=gauss_weight*(crLHS152*r_DN(3,0) + crLHS154*r_DN(3,1) + crLHS156*crLHS295 + crLHS301);
rLHS(13,15)+=-crLHS157*crLHS296;
rLHS(14,0)+=gauss_weight*(-crLHS102*crLHS173 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(0,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(0,1)*r_N[3]*tau_u - crLHS70);
rLHS(14,1)+=gauss_weight*(crLHS155 + crLHS158*r_DN(3,1) + crLHS276*crLHS302 + crLHS81*r_DN(3,0));
rLHS(14,2)+=gauss_weight*(crLHS100*crLHS303 + crLHS161*r_DN(3,1) + crLHS177 + crLHS294 + crLHS98*r_DN(3,0));
rLHS(14,3)+=-crLHS304*crLHS45;
rLHS(14,4)+=gauss_weight*(-crLHS102*crLHS240 - crLHS209 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(1,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(1,1)*r_N[3]*tau_u);
rLHS(14,5)+=gauss_weight*(crLHS105*r_DN(3,0) + crLHS108*crLHS303 + crLHS164*r_DN(3,1) + crLHS233);
rLHS(14,6)+=gauss_weight*(crLHS116*r_DN(3,0) + crLHS119*crLHS303 + crLHS166*r_DN(3,1) + crLHS242 + crLHS297);
rLHS(14,7)+=-crLHS120*crLHS304;
rLHS(14,8)+=gauss_weight*(-crLHS102*crLHS280 - crLHS256 + crLHS3*crLHS4*crLHS6*crLHS66*gamma*p_th*r_DN(2,1)*tau_u + crLHS3*crLHS4*crLHS6*crLHS73*gamma*p_th*r_DN(2,1)*r_N[3]*tau_u);
rLHS(14,9)+=gauss_weight*(crLHS125*r_DN(3,0) + crLHS128*crLHS303 + crLHS169*r_DN(3,1) + crLHS275);
rLHS(14,10)+=gauss_weight*(crLHS136*r_DN(3,0) + crLHS139*crLHS303 + crLHS171*r_DN(3,1) + crLHS283 + crLHS298);
rLHS(14,11)+=-crLHS140*crLHS304;
rLHS(14,12)+=crLHS299*r_DN(3,1);
rLHS(14,13)+=gauss_weight*(crLHS144*r_DN(3,0) + crLHS146*crLHS303 + crLHS174*r_DN(3,1) + crLHS301);
rLHS(14,14)+=gauss_weight*(crLHS154*r_DN(3,0) + crLHS156*crLHS303 + crLHS176*r_DN(3,1) + crLHS287*crLHS9 + crLHS300);
rLHS(14,15)+=-crLHS157*crLHS304;
rLHS(15,0)+=0;
rLHS(15,1)+=0;
rLHS(15,2)+=0;
rLHS(15,3)+=-gauss_weight*(crLHS180*crLHS26*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS180*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS180*crLHS26*dp_th_dt*r_N[3]*tau_t - crLHS180*crLHS305 - crLHS19*crLHS190 - crLHS191);
rLHS(15,4)+=0;
rLHS(15,5)+=0;
rLHS(15,6)+=0;
rLHS(15,7)+=-gauss_weight*(crLHS183*crLHS26*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS183*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS183*crLHS26*dp_th_dt*r_N[3]*tau_t - crLHS183*crLHS305 - crLHS190*crLHS37 - crLHS246);
rLHS(15,8)+=0;
rLHS(15,9)+=0;
rLHS(15,10)+=0;
rLHS(15,11)+=-gauss_weight*(crLHS186*crLHS26*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS186*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS186*crLHS26*dp_th_dt*r_N[3]*tau_t - crLHS186*crLHS305 - crLHS190*crLHS53 - crLHS285);
rLHS(15,12)+=0;
rLHS(15,13)+=0;
rLHS(15,14)+=0;
rLHS(15,15)+=-gauss_weight*(-crLHS178*crLHS289 + crLHS189*crLHS26*crLHS6*crLHS66*gamma*p_th*tau_t + crLHS189*crLHS26*crLHS6*crLHS73*gamma*p_th*r_N[3]*tau_t + crLHS189*crLHS26*dp_th_dt*r_N[3]*tau_t - crLHS189*crLHS305 - crLHS190*crLHS66 - crLHS286*kappa - crLHS287*kappa + crLHS288*crLHS3*dp_th_dt);

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

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1);
const double crRHS2 = crRHS0 + crRHS1;
const double crRHS3 = crRHS2*p_th;
const double crRHS4 = r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]);
const double crRHS5 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crRHS6 = 1.0/crRHS5;
const double crRHS7 = crRHS6*p_th;
const double crRHS8 = crRHS4*crRHS7;
const double crRHS9 = -crRHS8 + dp_th_dt;
const double crRHS10 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crRHS11 = crRHS10*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0));
const double crRHS12 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crRHS13 = crRHS12*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1));
const double crRHS14 = crRHS7*(crRHS11 + crRHS13);
const double crRHS15 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crRHS16 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crRHS17 = crRHS0*crRHS15 + crRHS16*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS18 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS19 = 1.0/c_p;
const double crRHS20 = gamma/(gamma - 1.0);
const double crRHS21 = crRHS20*crRHS7;
const double crRHS22 = crRHS19*crRHS21;
const double crRHS23 = -crRHS22*(-crRHS17 - crRHS18 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS24 = p_th*tau_u;
const double crRHS25 = crRHS23*crRHS24;
const double crRHS26 = crRHS1*crRHS16 + crRHS15*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS27 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS28 = -crRHS22*(-crRHS26 - crRHS27 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS29 = crRHS24*crRHS28;
const double crRHS30 = crRHS19*crRHS20;
const double crRHS31 = crRHS30*crRHS6;
const double crRHS32 = crRHS31*gauss_weight;
const double crRHS33 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS34 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS35 = crRHS22*r_N[0];
const double crRHS36 = crRHS31*r_DN(0,0);
const double crRHS37 = tau_c*(-crRHS11*crRHS7 - crRHS13*crRHS7 + crRHS2*p_th - crRHS8 + dp_th_dt);
const double crRHS38 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1);
const double crRHS39 = crRHS38*tau_u;
const double crRHS40 = crRHS35*crRHS39;
const double crRHS41 = crRHS15*r_DN(0,0) + crRHS16*r_DN(0,1);
const double crRHS42 = crRHS22*tau_u;
const double crRHS43 = crRHS41*crRHS42;
const double crRHS44 = (crRHS10*crRHS15 + crRHS12*crRHS16)*1.0/(crRHS5*crRHS5);
const double crRHS45 = crRHS44*r_N[0];
const double crRHS46 = crRHS30*crRHS45;
const double crRHS47 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS48 = crRHS31*r_DN(0,1);
const double crRHS49 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2];
const double crRHS50 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crRHS51 = crRHS50*kappa;
const double crRHS52 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crRHS53 = crRHS52*kappa;
const double crRHS54 = crRHS6*dp_th_dt;
const double crRHS55 = crRHS54*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2]);
const double crRHS56 = crRHS20*crRHS8;
const double crRHS57 = crRHS15*crRHS50 + crRHS16*crRHS52;
const double crRHS58 = crRHS21*crRHS57;
const double crRHS59 = tau_t*(-crRHS21*(crRHS4 + crRHS57) + crRHS49 + crRHS55);
const double crRHS60 = crRHS54*crRHS59;
const double crRHS61 = crRHS21*crRHS59;
const double crRHS62 = crRHS38*crRHS61;
const double crRHS63 = crRHS20*crRHS59*p_th;
const double crRHS64 = crRHS22*r_N[1];
const double crRHS65 = crRHS31*r_DN(1,0);
const double crRHS66 = crRHS39*crRHS64;
const double crRHS67 = crRHS15*r_DN(1,0) + crRHS16*r_DN(1,1);
const double crRHS68 = crRHS42*crRHS67;
const double crRHS69 = crRHS30*crRHS44;
const double crRHS70 = crRHS69*r_N[1];
const double crRHS71 = crRHS31*r_DN(1,1);
const double crRHS72 = crRHS44*crRHS63;
const double crRHS73 = crRHS22*r_N[2];
const double crRHS74 = crRHS31*r_DN(2,0);
const double crRHS75 = crRHS39*crRHS73;
const double crRHS76 = crRHS15*r_DN(2,0) + crRHS16*r_DN(2,1);
const double crRHS77 = crRHS42*crRHS76;
const double crRHS78 = crRHS69*r_N[2];
const double crRHS79 = crRHS31*r_DN(2,1);
rRHS[0]+=-crRHS32*(-crRHS14*r_N[0] + crRHS25*r_DN(0,0) + crRHS29*r_DN(0,1) + crRHS3*r_N[0] + crRHS9*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS17*crRHS35 + crRHS18*crRHS35 + crRHS23*crRHS40 + crRHS23*crRHS43 - crRHS25*crRHS46 + crRHS3*crRHS36 - crRHS33*r_DN(0,0) - crRHS34*crRHS35 + crRHS36*crRHS37 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS26*crRHS35 + crRHS27*crRHS35 + crRHS28*crRHS40 + crRHS28*crRHS43 - crRHS29*crRHS46 + crRHS3*crRHS48 - crRHS33*r_DN(0,1) - crRHS35*crRHS47 + crRHS37*crRHS48 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS41*crRHS61 - crRHS45*crRHS63 + crRHS49*r_N[0] - crRHS51*r_DN(0,0) - crRHS53*r_DN(0,1) + crRHS55*r_N[0] - crRHS56*r_N[0] - crRHS58*r_N[0] + crRHS60*r_N[0] + crRHS62*r_N[0]);
rRHS[4]+=-crRHS32*(-crRHS14*r_N[1] + crRHS25*r_DN(1,0) + crRHS29*r_DN(1,1) + crRHS3*r_N[1] + crRHS9*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS17*crRHS64 + crRHS18*crRHS64 + crRHS23*crRHS66 + crRHS23*crRHS68 - crRHS25*crRHS70 + crRHS3*crRHS65 - crRHS33*r_DN(1,0) - crRHS34*crRHS64 + crRHS37*crRHS65 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS26*crRHS64 + crRHS27*crRHS64 + crRHS28*crRHS66 + crRHS28*crRHS68 - crRHS29*crRHS70 + crRHS3*crRHS71 - crRHS33*r_DN(1,1) + crRHS37*crRHS71 - crRHS47*crRHS64 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS49*r_N[1] - crRHS51*r_DN(1,0) - crRHS53*r_DN(1,1) + crRHS55*r_N[1] - crRHS56*r_N[1] - crRHS58*r_N[1] + crRHS60*r_N[1] + crRHS61*crRHS67 + crRHS62*r_N[1] - crRHS72*r_N[1]);
rRHS[8]+=-crRHS32*(-crRHS14*r_N[2] + crRHS25*r_DN(2,0) + crRHS29*r_DN(2,1) + crRHS3*r_N[2] + crRHS9*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS17*crRHS73 + crRHS18*crRHS73 + crRHS23*crRHS75 + crRHS23*crRHS77 - crRHS25*crRHS78 + crRHS3*crRHS74 - crRHS33*r_DN(2,0) - crRHS34*crRHS73 + crRHS37*crRHS74 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS26*crRHS73 + crRHS27*crRHS73 + crRHS28*crRHS75 + crRHS28*crRHS77 - crRHS29*crRHS78 + crRHS3*crRHS79 - crRHS33*r_DN(2,1) + crRHS37*crRHS79 - crRHS47*crRHS73 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS49*r_N[2] - crRHS51*r_DN(2,0) - crRHS53*r_DN(2,1) + crRHS55*r_N[2] - crRHS56*r_N[2] - crRHS58*r_N[2] + crRHS60*r_N[2] + crRHS61*crRHS76 + crRHS62*r_N[2] - crRHS72*r_N[2]);

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

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0) + r_DN(3,0)*r_u(3,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1) + r_DN(3,1)*r_u(3,1);
const double crRHS2 = crRHS0 + crRHS1;
const double crRHS3 = crRHS2*p_th;
const double crRHS4 = r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]) + r_N[3]*(bdf0*r_t[3] + bdf1*r_t_n[3] + bdf2*r_t_nn[3]);
const double crRHS5 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crRHS6 = 1.0/crRHS5;
const double crRHS7 = crRHS6*p_th;
const double crRHS8 = crRHS4*crRHS7;
const double crRHS9 = -crRHS8 + dp_th_dt;
const double crRHS10 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crRHS11 = crRHS10*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0) + r_N[3]*r_u(3,0));
const double crRHS12 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crRHS13 = crRHS12*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1) + r_N[3]*r_u(3,1));
const double crRHS14 = crRHS7*(crRHS11 + crRHS13);
const double crRHS15 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crRHS16 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crRHS17 = crRHS0*crRHS15 + crRHS16*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS18 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0)) + r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS19 = 1.0/c_p;
const double crRHS20 = gamma/(gamma - 1.0);
const double crRHS21 = crRHS20*crRHS7;
const double crRHS22 = crRHS19*crRHS21;
const double crRHS23 = -crRHS22*(-crRHS17 - crRHS18 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS24 = p_th*tau_u;
const double crRHS25 = crRHS23*crRHS24;
const double crRHS26 = crRHS1*crRHS16 + crRHS15*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS27 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1)) + r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS28 = -crRHS22*(-crRHS26 - crRHS27 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS29 = crRHS24*crRHS28;
const double crRHS30 = crRHS19*crRHS20;
const double crRHS31 = crRHS30*crRHS6;
const double crRHS32 = crRHS31*gauss_weight;
const double crRHS33 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS34 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS35 = crRHS22*r_N[0];
const double crRHS36 = crRHS31*r_DN(0,0);
const double crRHS37 = tau_c*(-crRHS11*crRHS7 - crRHS13*crRHS7 + crRHS2*p_th - crRHS8 + dp_th_dt);
const double crRHS38 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1);
const double crRHS39 = crRHS38*tau_u;
const double crRHS40 = crRHS35*crRHS39;
const double crRHS41 = crRHS15*r_DN(0,0) + crRHS16*r_DN(0,1);
const double crRHS42 = crRHS22*tau_u;
const double crRHS43 = crRHS41*crRHS42;
const double crRHS44 = (crRHS10*crRHS15 + crRHS12*crRHS16)*1.0/(crRHS5*crRHS5);
const double crRHS45 = crRHS44*r_N[0];
const double crRHS46 = crRHS30*crRHS45;
const double crRHS47 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS48 = crRHS31*r_DN(0,1);
const double crRHS49 = r_N[0]*r_heat_fl[0] + r_N[1]*r_heat_fl[1] + r_N[2]*r_heat_fl[2] + r_N[3]*r_heat_fl[3];
const double crRHS50 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crRHS51 = crRHS50*kappa;
const double crRHS52 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crRHS53 = crRHS52*kappa;
const double crRHS54 = crRHS6*dp_th_dt;
const double crRHS55 = crRHS54*(r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3]);
const double crRHS56 = crRHS20*crRHS8;
const double crRHS57 = crRHS15*crRHS50 + crRHS16*crRHS52;
const double crRHS58 = crRHS21*crRHS57;
const double crRHS59 = tau_t*(-crRHS21*(crRHS4 + crRHS57) + crRHS49 + crRHS55);
const double crRHS60 = crRHS54*crRHS59;
const double crRHS61 = crRHS21*crRHS59;
const double crRHS62 = crRHS38*crRHS61;
const double crRHS63 = crRHS20*crRHS59*p_th;
const double crRHS64 = crRHS22*r_N[1];
const double crRHS65 = crRHS31*r_DN(1,0);
const double crRHS66 = crRHS39*crRHS64;
const double crRHS67 = crRHS15*r_DN(1,0) + crRHS16*r_DN(1,1);
const double crRHS68 = crRHS42*crRHS67;
const double crRHS69 = crRHS30*crRHS44;
const double crRHS70 = crRHS69*r_N[1];
const double crRHS71 = crRHS31*r_DN(1,1);
const double crRHS72 = crRHS44*crRHS63;
const double crRHS73 = crRHS22*r_N[2];
const double crRHS74 = crRHS31*r_DN(2,0);
const double crRHS75 = crRHS39*crRHS73;
const double crRHS76 = crRHS15*r_DN(2,0) + crRHS16*r_DN(2,1);
const double crRHS77 = crRHS42*crRHS76;
const double crRHS78 = crRHS69*r_N[2];
const double crRHS79 = crRHS31*r_DN(2,1);
const double crRHS80 = crRHS22*r_N[3];
const double crRHS81 = crRHS31*r_DN(3,0);
const double crRHS82 = crRHS39*crRHS80;
const double crRHS83 = crRHS15*r_DN(3,0) + crRHS16*r_DN(3,1);
const double crRHS84 = crRHS42*crRHS83;
const double crRHS85 = crRHS69*r_N[3];
const double crRHS86 = crRHS31*r_DN(3,1);
rRHS[0]+=-crRHS32*(-crRHS14*r_N[0] + crRHS25*r_DN(0,0) + crRHS29*r_DN(0,1) + crRHS3*r_N[0] + crRHS9*r_N[0]);
rRHS[1]+=-gauss_weight*(crRHS17*crRHS35 + crRHS18*crRHS35 + crRHS23*crRHS40 + crRHS23*crRHS43 - crRHS25*crRHS46 + crRHS3*crRHS36 - crRHS33*r_DN(0,0) - crRHS34*crRHS35 + crRHS36*crRHS37 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-gauss_weight*(crRHS26*crRHS35 + crRHS27*crRHS35 + crRHS28*crRHS40 + crRHS28*crRHS43 - crRHS29*crRHS46 + crRHS3*crRHS48 - crRHS33*r_DN(0,1) - crRHS35*crRHS47 + crRHS37*crRHS48 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=gauss_weight*(crRHS41*crRHS61 - crRHS45*crRHS63 + crRHS49*r_N[0] - crRHS51*r_DN(0,0) - crRHS53*r_DN(0,1) + crRHS55*r_N[0] - crRHS56*r_N[0] - crRHS58*r_N[0] + crRHS60*r_N[0] + crRHS62*r_N[0]);
rRHS[4]+=-crRHS32*(-crRHS14*r_N[1] + crRHS25*r_DN(1,0) + crRHS29*r_DN(1,1) + crRHS3*r_N[1] + crRHS9*r_N[1]);
rRHS[5]+=-gauss_weight*(crRHS17*crRHS64 + crRHS18*crRHS64 + crRHS23*crRHS66 + crRHS23*crRHS68 - crRHS25*crRHS70 + crRHS3*crRHS65 - crRHS33*r_DN(1,0) - crRHS34*crRHS64 + crRHS37*crRHS65 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-gauss_weight*(crRHS26*crRHS64 + crRHS27*crRHS64 + crRHS28*crRHS66 + crRHS28*crRHS68 - crRHS29*crRHS70 + crRHS3*crRHS71 - crRHS33*r_DN(1,1) + crRHS37*crRHS71 - crRHS47*crRHS64 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=gauss_weight*(crRHS49*r_N[1] - crRHS51*r_DN(1,0) - crRHS53*r_DN(1,1) + crRHS55*r_N[1] - crRHS56*r_N[1] - crRHS58*r_N[1] + crRHS60*r_N[1] + crRHS61*crRHS67 + crRHS62*r_N[1] - crRHS72*r_N[1]);
rRHS[8]+=-crRHS32*(-crRHS14*r_N[2] + crRHS25*r_DN(2,0) + crRHS29*r_DN(2,1) + crRHS3*r_N[2] + crRHS9*r_N[2]);
rRHS[9]+=-gauss_weight*(crRHS17*crRHS73 + crRHS18*crRHS73 + crRHS23*crRHS75 + crRHS23*crRHS77 - crRHS25*crRHS78 + crRHS3*crRHS74 - crRHS33*r_DN(2,0) - crRHS34*crRHS73 + crRHS37*crRHS74 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-gauss_weight*(crRHS26*crRHS73 + crRHS27*crRHS73 + crRHS28*crRHS75 + crRHS28*crRHS77 - crRHS29*crRHS78 + crRHS3*crRHS79 - crRHS33*r_DN(2,1) + crRHS37*crRHS79 - crRHS47*crRHS73 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=gauss_weight*(crRHS49*r_N[2] - crRHS51*r_DN(2,0) - crRHS53*r_DN(2,1) + crRHS55*r_N[2] - crRHS56*r_N[2] - crRHS58*r_N[2] + crRHS60*r_N[2] + crRHS61*crRHS76 + crRHS62*r_N[2] - crRHS72*r_N[2]);
rRHS[12]+=-crRHS32*(-crRHS14*r_N[3] + crRHS25*r_DN(3,0) + crRHS29*r_DN(3,1) + crRHS3*r_N[3] + crRHS9*r_N[3]);
rRHS[13]+=-gauss_weight*(crRHS17*crRHS80 + crRHS18*crRHS80 + crRHS23*crRHS82 + crRHS23*crRHS84 - crRHS25*crRHS85 + crRHS3*crRHS81 - crRHS33*r_DN(3,0) - crRHS34*crRHS80 + crRHS37*crRHS81 + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-gauss_weight*(crRHS26*crRHS80 + crRHS27*crRHS80 + crRHS28*crRHS82 + crRHS28*crRHS84 - crRHS29*crRHS85 + crRHS3*crRHS86 - crRHS33*r_DN(3,1) + crRHS37*crRHS86 - crRHS47*crRHS80 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=gauss_weight*(crRHS49*r_N[3] - crRHS51*r_DN(3,0) - crRHS53*r_DN(3,1) + crRHS55*r_N[3] - crRHS56*r_N[3] - crRHS58*r_N[3] + crRHS60*r_N[3] + crRHS61*crRHS83 + crRHS62*r_N[3] - crRHS72*r_N[3]);

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