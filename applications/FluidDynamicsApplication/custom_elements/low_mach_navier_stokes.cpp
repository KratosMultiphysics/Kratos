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
    const double check_tolerance = 1.0e-8;
    const auto& r_geometry = this->GetGeometry();
    for (const auto& r_node : r_geometry) {
        // Check nodal DOFs
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE,r_node);
        // Axisymmetry: nodes are in XY plane
        KRATOS_ERROR_IF(std::abs(r_node.Z()) > check_tolerance ) << "Node " << r_node.Id() << "has non-zero Z coordinate." << std::endl;
        // Axisymmetry: check that there are no negative y-coordinates (radius is always positive)
        KRATOS_ERROR_IF(r_node.Y() < 0.0) << "Negative y-coordinate found in node " << r_node.Id() << ". Axisymmetric radius must be positive." << std::endl;
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
    const auto& r_rho = rData.Density;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

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
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = 1.0/c_p;
const double crLHS6 = gamma - 1.0;
const double crLHS7 = 1.0/crLHS6;
const double crLHS8 = crLHS5*crLHS7*gamma*p_th;
const double crLHS9 = crLHS4*crLHS8;
const double crLHS10 = crLHS2*crLHS9;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = crLHS11*(crLHS0 + crLHS1);
const double crLHS13 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS14 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crLHS17 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS9*tau_u;
const double crLHS21 = crLHS19*crLHS20;
const double crLHS22 = crLHS10*(-crLHS13*crLHS14 + crLHS21*r_DN(0,0) + r_DN(0,0)*r_N[0]);
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS24 = crLHS10*(-crLHS14*crLHS23 + crLHS21*r_DN(0,1) + r_DN(0,1)*r_N[0]);
const double crLHS25 = r_DN(0,0)*r_DN(1,0);
const double crLHS26 = r_DN(0,1)*r_DN(1,1);
const double crLHS27 = crLHS11*(crLHS25 + crLHS26);
const double crLHS28 = r_DN(1,0)*r_N[0];
const double crLHS29 = crLHS4*r_N[0];
const double crLHS30 = crLHS13*crLHS29;
const double crLHS31 = -crLHS30*r_N[1];
const double crLHS32 = bdf0*r_N[1];
const double crLHS33 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS34 = crLHS32 + crLHS33;
const double crLHS35 = crLHS20*crLHS34;
const double crLHS36 = crLHS10*(crLHS28 + crLHS31 + crLHS35*r_DN(0,0));
const double crLHS37 = r_DN(1,1)*r_N[0];
const double crLHS38 = crLHS23*crLHS29;
const double crLHS39 = -crLHS38*r_N[1];
const double crLHS40 = crLHS10*(crLHS35*r_DN(0,1) + crLHS37 + crLHS39);
const double crLHS41 = r_DN(0,0)*r_DN(2,0);
const double crLHS42 = r_DN(0,1)*r_DN(2,1);
const double crLHS43 = crLHS11*(crLHS41 + crLHS42);
const double crLHS44 = r_DN(2,0)*r_N[0];
const double crLHS45 = -crLHS30*r_N[2];
const double crLHS46 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS47 = bdf0*r_N[2] + crLHS46;
const double crLHS48 = crLHS20*crLHS47;
const double crLHS49 = crLHS10*(crLHS44 + crLHS45 + crLHS48*r_DN(0,0));
const double crLHS50 = r_DN(2,1)*r_N[0];
const double crLHS51 = -crLHS38*r_N[2];
const double crLHS52 = crLHS10*(crLHS48*r_DN(0,1) + crLHS50 + crLHS51);
const double crLHS53 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1);
const double crLHS54 = crLHS13*crLHS16 + crLHS17*crLHS23;
const double crLHS55 = crLHS54*r_N[0];
const double crLHS56 = 1.0/(crLHS3*crLHS3);
const double crLHS57 = crLHS56*tau_u;
const double crLHS58 = crLHS57*crLHS8;
const double crLHS59 = crLHS2*(crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_N[0]*tau_u - crLHS55*crLHS58 - r_N[0]);
const double crLHS60 = crLHS59*r_DN(0,0);
const double crLHS61 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS62 = r_C(0,2)*r_DN(0,0);
const double crLHS63 = crLHS62 + r_C(2,2)*r_DN(0,1);
const double crLHS64 = -crLHS30 + r_DN(0,0);
const double crLHS65 = crLHS9*r_DN(0,0);
const double crLHS66 = crLHS65*tau_c;
const double crLHS67 = bdf0*crLHS8;
const double crLHS68 = crLHS29*crLHS8;
const double crLHS69 = 1.0/(crLHS6*crLHS6)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS70 = crLHS57*crLHS69;
const double crLHS71 = crLHS18*crLHS70;
const double crLHS72 = crLHS19*crLHS69;
const double crLHS73 = crLHS53*crLHS56*tau_u;
const double crLHS74 = crLHS73*r_N[0];
const double crLHS75 = tau_u*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS76 = crLHS55*crLHS75;
const double crLHS77 = crLHS14*crLHS67 + crLHS18*crLHS68 + crLHS19*crLHS71 + crLHS72*crLHS74 - crLHS72*crLHS76;
const double crLHS78 = crLHS62 + r_C(0,1)*r_DN(0,1);
const double crLHS79 = r_C(1,2)*r_DN(0,1);
const double crLHS80 = crLHS79 + r_C(2,2)*r_DN(0,0);
const double crLHS81 = crLHS65*r_DN(0,1);
const double crLHS82 = -crLHS38 + r_DN(0,1);
const double crLHS83 = r_DN(0,0)*r_N[1];
const double crLHS84 = crLHS54*crLHS58;
const double crLHS85 = crLHS2*(crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(1,0)*tau_u - crLHS28*crLHS84 + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS83);
const double crLHS86 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS87 = r_C(0,2)*r_DN(1,0);
const double crLHS88 = crLHS87 + r_C(2,2)*r_DN(1,1);
const double crLHS89 = crLHS4*r_N[1];
const double crLHS90 = crLHS13*crLHS89;
const double crLHS91 = -crLHS90 + r_DN(1,0);
const double crLHS92 = crLHS8*crLHS89;
const double crLHS93 = crLHS15*crLHS92;
const double crLHS94 = crLHS25*crLHS9 + crLHS93;
const double crLHS95 = crLHS34*crLHS69;
const double crLHS96 = crLHS33*crLHS68 + crLHS34*crLHS71 + crLHS74*crLHS95 - crLHS76*crLHS95;
const double crLHS97 = crLHS87 + r_C(0,1)*r_DN(1,1);
const double crLHS98 = r_C(1,2)*r_DN(1,1);
const double crLHS99 = crLHS98 + r_C(2,2)*r_DN(1,0);
const double crLHS100 = crLHS65*r_DN(1,1);
const double crLHS101 = crLHS23*crLHS89;
const double crLHS102 = -crLHS101 + r_DN(1,1);
const double crLHS103 = r_DN(0,0)*r_N[2];
const double crLHS104 = crLHS2*(-crLHS103 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,0)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS44*crLHS84);
const double crLHS105 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS106 = r_C(0,2)*r_DN(2,0);
const double crLHS107 = crLHS106 + r_C(2,2)*r_DN(2,1);
const double crLHS108 = crLHS4*r_N[2];
const double crLHS109 = -crLHS108*crLHS13 + r_DN(2,0);
const double crLHS110 = crLHS108*crLHS8;
const double crLHS111 = crLHS110*crLHS15;
const double crLHS112 = crLHS111 + crLHS41*crLHS9;
const double crLHS113 = crLHS47*crLHS69;
const double crLHS114 = crLHS113*crLHS74 - crLHS113*crLHS76 + crLHS46*crLHS68 + crLHS47*crLHS71;
const double crLHS115 = crLHS106 + r_C(0,1)*r_DN(2,1);
const double crLHS116 = r_C(1,2)*r_DN(2,1);
const double crLHS117 = crLHS116 + r_C(2,2)*r_DN(2,0);
const double crLHS118 = crLHS65*r_DN(2,1);
const double crLHS119 = -crLHS108*crLHS23 + r_DN(2,1);
const double crLHS120 = crLHS59*r_DN(0,1);
const double crLHS121 = crLHS79 + r_C(0,1)*r_DN(0,0);
const double crLHS122 = crLHS9*r_DN(0,1);
const double crLHS123 = crLHS122*tau_c;
const double crLHS124 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS125 = r_DN(0,1)*r_N[1];
const double crLHS126 = crLHS2*(-crLHS125 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(1,1)*tau_u - crLHS37*crLHS84 + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u);
const double crLHS127 = crLHS98 + r_C(0,1)*r_DN(1,0);
const double crLHS128 = crLHS122*r_DN(1,0);
const double crLHS129 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS130 = crLHS26*crLHS9 + crLHS93;
const double crLHS131 = r_DN(0,1)*r_N[2];
const double crLHS132 = crLHS2*(-crLHS131 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,1)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS50*crLHS84);
const double crLHS133 = crLHS116 + r_C(0,1)*r_DN(2,0);
const double crLHS134 = crLHS122*r_DN(2,0);
const double crLHS135 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS136 = crLHS111 + crLHS42*crLHS9;
const double crLHS137 = crLHS10*(crLHS21*r_DN(1,0) + crLHS31 + crLHS83);
const double crLHS138 = crLHS10*(crLHS125 + crLHS21*r_DN(1,1) + crLHS39);
const double crLHS139 = r_DN(1,0)*r_DN(1,0);
const double crLHS140 = r_DN(1,1)*r_DN(1,1);
const double crLHS141 = crLHS11*(crLHS139 + crLHS140);
const double crLHS142 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS143 = crLHS10*(-crLHS13*crLHS142 + crLHS35*r_DN(1,0) + r_DN(1,0)*r_N[1]);
const double crLHS144 = crLHS10*(-crLHS142*crLHS23 + crLHS35*r_DN(1,1) + r_DN(1,1)*r_N[1]);
const double crLHS145 = r_DN(1,0)*r_DN(2,0);
const double crLHS146 = r_DN(1,1)*r_DN(2,1);
const double crLHS147 = crLHS11*(crLHS145 + crLHS146);
const double crLHS148 = r_DN(2,0)*r_N[1];
const double crLHS149 = -crLHS90*r_N[2];
const double crLHS150 = crLHS10*(crLHS148 + crLHS149 + crLHS48*r_DN(1,0));
const double crLHS151 = r_DN(2,1)*r_N[1];
const double crLHS152 = -crLHS101*r_N[2];
const double crLHS153 = crLHS10*(crLHS151 + crLHS152 + crLHS48*r_DN(1,1));
const double crLHS154 = crLHS2*(-crLHS28 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(0,0)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS83*crLHS84);
const double crLHS155 = crLHS9*r_DN(1,0);
const double crLHS156 = crLHS155*tau_c;
const double crLHS157 = crLHS33*crLHS70;
const double crLHS158 = crLHS72*r_N[1];
const double crLHS159 = crLHS54*crLHS75;
const double crLHS160 = crLHS157*crLHS19 - crLHS158*crLHS159 + crLHS158*crLHS73 + crLHS18*crLHS92;
const double crLHS161 = crLHS2*(crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_N[1]*tau_u - crLHS84*r_N[1] - r_N[1]);
const double crLHS162 = crLHS161*r_DN(1,0);
const double crLHS163 = crLHS95*r_N[1];
const double crLHS164 = crLHS142*crLHS67 + crLHS157*crLHS34 - crLHS159*crLHS163 + crLHS163*crLHS73 + crLHS33*crLHS92;
const double crLHS165 = crLHS155*r_DN(1,1);
const double crLHS166 = r_DN(1,0)*r_N[2];
const double crLHS167 = crLHS2*(-crLHS148*crLHS84 - crLHS166 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,0)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
const double crLHS168 = crLHS110*crLHS32;
const double crLHS169 = crLHS145*crLHS9 + crLHS168;
const double crLHS170 = crLHS113*r_N[1];
const double crLHS171 = crLHS157*crLHS47 - crLHS159*crLHS170 + crLHS170*crLHS73 + crLHS46*crLHS92;
const double crLHS172 = crLHS155*r_DN(2,1);
const double crLHS173 = crLHS2*(-crLHS125*crLHS84 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(0,1)*tau_u - crLHS37 + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u);
const double crLHS174 = crLHS64*tau_c;
const double crLHS175 = crLHS9*r_DN(1,1);
const double crLHS176 = crLHS175*tau_c;
const double crLHS177 = crLHS161*r_DN(1,1);
const double crLHS178 = r_DN(1,1)*r_N[2];
const double crLHS179 = crLHS2*(-crLHS151*crLHS84 - crLHS178 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,1)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
const double crLHS180 = crLHS9*r_DN(2,0);
const double crLHS181 = crLHS180*r_DN(1,1);
const double crLHS182 = crLHS146*crLHS9 + crLHS168;
const double crLHS183 = crLHS10*(crLHS103 + crLHS21*r_DN(2,0) + crLHS45);
const double crLHS184 = crLHS10*(crLHS131 + crLHS21*r_DN(2,1) + crLHS51);
const double crLHS185 = crLHS10*(crLHS149 + crLHS166 + crLHS35*r_DN(2,0));
const double crLHS186 = crLHS10*(crLHS152 + crLHS178 + crLHS35*r_DN(2,1));
const double crLHS187 = r_DN(2,0)*r_DN(2,0);
const double crLHS188 = r_DN(2,1)*r_DN(2,1);
const double crLHS189 = crLHS11*(crLHS187 + crLHS188);
const double crLHS190 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS191 = crLHS10*(-crLHS13*crLHS190 + crLHS48*r_DN(2,0) + r_DN(2,0)*r_N[2]);
const double crLHS192 = crLHS10*(-crLHS190*crLHS23 + crLHS48*r_DN(2,1) + r_DN(2,1)*r_N[2]);
const double crLHS193 = crLHS2*(-crLHS103*crLHS84 + crLHS4*crLHS46*crLHS5*crLHS7*gamma*p_th*r_DN(0,0)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS44);
const double crLHS194 = crLHS46*crLHS70;
const double crLHS195 = crLHS72*r_N[2];
const double crLHS196 = crLHS110*crLHS18 - crLHS159*crLHS195 + crLHS19*crLHS194 + crLHS195*crLHS73;
const double crLHS197 = crLHS180*tau_c;
const double crLHS198 = crLHS2*(-crLHS148 - crLHS166*crLHS84 + crLHS4*crLHS46*crLHS5*crLHS7*gamma*p_th*r_DN(1,0)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
const double crLHS199 = crLHS95*r_N[2];
const double crLHS200 = crLHS110*crLHS33 - crLHS159*crLHS199 + crLHS194*crLHS34 + crLHS199*crLHS73;
const double crLHS201 = crLHS2*(crLHS4*crLHS46*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_N[2]*tau_u - crLHS84*r_N[2] - r_N[2]);
const double crLHS202 = crLHS201*r_DN(2,0);
const double crLHS203 = crLHS113*r_N[2];
const double crLHS204 = crLHS110*crLHS46 - crLHS159*crLHS203 + crLHS190*crLHS67 + crLHS194*crLHS47 + crLHS203*crLHS73;
const double crLHS205 = crLHS180*r_DN(2,1);
const double crLHS206 = crLHS2*(-crLHS131*crLHS84 + crLHS4*crLHS46*crLHS5*crLHS7*gamma*p_th*r_DN(0,1)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS50);
const double crLHS207 = crLHS9*r_DN(2,1);
const double crLHS208 = crLHS207*tau_c;
const double crLHS209 = crLHS2*(-crLHS151 - crLHS178*crLHS84 + crLHS4*crLHS46*crLHS5*crLHS7*gamma*p_th*r_DN(1,1)*tau_u + crLHS4*crLHS5*crLHS53*crLHS7*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
const double crLHS210 = crLHS201*r_DN(2,1);
rLHS(0,0)+=crLHS12;
rLHS(0,1)+=crLHS22;
rLHS(0,2)+=crLHS24;
rLHS(0,3)+=crLHS12;
rLHS(0,4)+=crLHS27;
rLHS(0,5)+=crLHS36;
rLHS(0,6)+=crLHS40;
rLHS(0,7)+=crLHS27;
rLHS(0,8)+=crLHS43;
rLHS(0,9)+=crLHS49;
rLHS(0,10)+=crLHS52;
rLHS(0,11)+=crLHS43;
rLHS(1,0)+=crLHS60;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS9 + crLHS61*r_DN(0,0) + crLHS63*r_DN(0,1) + crLHS64*crLHS66 + crLHS77);
rLHS(1,2)+=crLHS2*(crLHS66*crLHS82 + crLHS78*r_DN(0,0) + crLHS80*r_DN(0,1) + crLHS81);
rLHS(1,3)+=crLHS60;
rLHS(1,4)+=crLHS85;
rLHS(1,5)+=crLHS2*(crLHS66*crLHS91 + crLHS86*r_DN(0,0) + crLHS88*r_DN(0,1) + crLHS94 + crLHS96);
rLHS(1,6)+=crLHS2*(crLHS100 + crLHS102*crLHS66 + crLHS97*r_DN(0,0) + crLHS99*r_DN(0,1));
rLHS(1,7)+=crLHS85;
rLHS(1,8)+=crLHS104;
rLHS(1,9)+=crLHS2*(crLHS105*r_DN(0,0) + crLHS107*r_DN(0,1) + crLHS109*crLHS66 + crLHS112 + crLHS114);
rLHS(1,10)+=crLHS2*(crLHS115*r_DN(0,0) + crLHS117*r_DN(0,1) + crLHS118 + crLHS119*crLHS66);
rLHS(1,11)+=crLHS104;
rLHS(2,0)+=crLHS120;
rLHS(2,1)+=crLHS2*(crLHS121*r_DN(0,1) + crLHS123*crLHS64 + crLHS63*r_DN(0,0) + crLHS81);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS9 + crLHS123*crLHS82 + crLHS124*r_DN(0,1) + crLHS77 + crLHS80*r_DN(0,0));
rLHS(2,3)+=crLHS120;
rLHS(2,4)+=crLHS126;
rLHS(2,5)+=crLHS2*(crLHS123*crLHS91 + crLHS127*r_DN(0,1) + crLHS128 + crLHS88*r_DN(0,0));
rLHS(2,6)+=crLHS2*(crLHS102*crLHS123 + crLHS129*r_DN(0,1) + crLHS130 + crLHS96 + crLHS99*r_DN(0,0));
rLHS(2,7)+=crLHS126;
rLHS(2,8)+=crLHS132;
rLHS(2,9)+=crLHS2*(crLHS107*r_DN(0,0) + crLHS109*crLHS123 + crLHS133*r_DN(0,1) + crLHS134);
rLHS(2,10)+=crLHS2*(crLHS114 + crLHS117*r_DN(0,0) + crLHS119*crLHS123 + crLHS135*r_DN(0,1) + crLHS136);
rLHS(2,11)+=crLHS132;
rLHS(3,0)+=crLHS12;
rLHS(3,1)+=crLHS22;
rLHS(3,2)+=crLHS24;
rLHS(3,3)+=crLHS12;
rLHS(3,4)+=crLHS27;
rLHS(3,5)+=crLHS36;
rLHS(3,6)+=crLHS40;
rLHS(3,7)+=crLHS27;
rLHS(3,8)+=crLHS43;
rLHS(3,9)+=crLHS49;
rLHS(3,10)+=crLHS52;
rLHS(3,11)+=crLHS43;
rLHS(4,0)+=crLHS27;
rLHS(4,1)+=crLHS137;
rLHS(4,2)+=crLHS138;
rLHS(4,3)+=crLHS27;
rLHS(4,4)+=crLHS141;
rLHS(4,5)+=crLHS143;
rLHS(4,6)+=crLHS144;
rLHS(4,7)+=crLHS141;
rLHS(4,8)+=crLHS147;
rLHS(4,9)+=crLHS150;
rLHS(4,10)+=crLHS153;
rLHS(4,11)+=crLHS147;
rLHS(5,0)+=crLHS154;
rLHS(5,1)+=crLHS2*(crLHS156*crLHS64 + crLHS160 + crLHS61*r_DN(1,0) + crLHS63*r_DN(1,1) + crLHS94);
rLHS(5,2)+=crLHS2*(crLHS128 + crLHS156*crLHS82 + crLHS78*r_DN(1,0) + crLHS80*r_DN(1,1));
rLHS(5,3)+=crLHS154;
rLHS(5,4)+=crLHS162;
rLHS(5,5)+=crLHS2*(crLHS139*crLHS9 + crLHS156*crLHS91 + crLHS164 + crLHS86*r_DN(1,0) + crLHS88*r_DN(1,1));
rLHS(5,6)+=crLHS2*(crLHS102*crLHS156 + crLHS165 + crLHS97*r_DN(1,0) + crLHS99*r_DN(1,1));
rLHS(5,7)+=crLHS162;
rLHS(5,8)+=crLHS167;
rLHS(5,9)+=crLHS2*(crLHS105*r_DN(1,0) + crLHS107*r_DN(1,1) + crLHS109*crLHS156 + crLHS169 + crLHS171);
rLHS(5,10)+=crLHS2*(crLHS115*r_DN(1,0) + crLHS117*r_DN(1,1) + crLHS119*crLHS156 + crLHS172);
rLHS(5,11)+=crLHS167;
rLHS(6,0)+=crLHS173;
rLHS(6,1)+=crLHS2*(crLHS100 + crLHS121*r_DN(1,1) + crLHS174*crLHS175 + crLHS63*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS124*r_DN(1,1) + crLHS130 + crLHS160 + crLHS176*crLHS82 + crLHS80*r_DN(1,0));
rLHS(6,3)+=crLHS173;
rLHS(6,4)+=crLHS177;
rLHS(6,5)+=crLHS2*(crLHS127*r_DN(1,1) + crLHS165 + crLHS176*crLHS91 + crLHS88*r_DN(1,0));
rLHS(6,6)+=crLHS2*(crLHS102*crLHS176 + crLHS129*r_DN(1,1) + crLHS140*crLHS9 + crLHS164 + crLHS99*r_DN(1,0));
rLHS(6,7)+=crLHS177;
rLHS(6,8)+=crLHS179;
rLHS(6,9)+=crLHS2*(crLHS107*r_DN(1,0) + crLHS109*crLHS176 + crLHS133*r_DN(1,1) + crLHS181);
rLHS(6,10)+=crLHS2*(crLHS117*r_DN(1,0) + crLHS119*crLHS176 + crLHS135*r_DN(1,1) + crLHS171 + crLHS182);
rLHS(6,11)+=crLHS179;
rLHS(7,0)+=crLHS27;
rLHS(7,1)+=crLHS137;
rLHS(7,2)+=crLHS138;
rLHS(7,3)+=crLHS27;
rLHS(7,4)+=crLHS141;
rLHS(7,5)+=crLHS143;
rLHS(7,6)+=crLHS144;
rLHS(7,7)+=crLHS141;
rLHS(7,8)+=crLHS147;
rLHS(7,9)+=crLHS150;
rLHS(7,10)+=crLHS153;
rLHS(7,11)+=crLHS147;
rLHS(8,0)+=crLHS43;
rLHS(8,1)+=crLHS183;
rLHS(8,2)+=crLHS184;
rLHS(8,3)+=crLHS43;
rLHS(8,4)+=crLHS147;
rLHS(8,5)+=crLHS185;
rLHS(8,6)+=crLHS186;
rLHS(8,7)+=crLHS147;
rLHS(8,8)+=crLHS189;
rLHS(8,9)+=crLHS191;
rLHS(8,10)+=crLHS192;
rLHS(8,11)+=crLHS189;
rLHS(9,0)+=crLHS193;
rLHS(9,1)+=crLHS2*(crLHS112 + crLHS174*crLHS180 + crLHS196 + crLHS61*r_DN(2,0) + crLHS63*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS134 + crLHS197*crLHS82 + crLHS78*r_DN(2,0) + crLHS80*r_DN(2,1));
rLHS(9,3)+=crLHS193;
rLHS(9,4)+=crLHS198;
rLHS(9,5)+=crLHS2*(crLHS169 + crLHS197*crLHS91 + crLHS200 + crLHS86*r_DN(2,0) + crLHS88*r_DN(2,1));
rLHS(9,6)+=crLHS2*(crLHS102*crLHS197 + crLHS181 + crLHS97*r_DN(2,0) + crLHS99*r_DN(2,1));
rLHS(9,7)+=crLHS198;
rLHS(9,8)+=crLHS202;
rLHS(9,9)+=crLHS2*(crLHS105*r_DN(2,0) + crLHS107*r_DN(2,1) + crLHS109*crLHS197 + crLHS187*crLHS9 + crLHS204);
rLHS(9,10)+=crLHS2*(crLHS115*r_DN(2,0) + crLHS117*r_DN(2,1) + crLHS119*crLHS197 + crLHS205);
rLHS(9,11)+=crLHS202;
rLHS(10,0)+=crLHS206;
rLHS(10,1)+=crLHS2*(crLHS118 + crLHS121*r_DN(2,1) + crLHS174*crLHS207 + crLHS63*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS124*r_DN(2,1) + crLHS136 + crLHS196 + crLHS208*crLHS82 + crLHS80*r_DN(2,0));
rLHS(10,3)+=crLHS206;
rLHS(10,4)+=crLHS209;
rLHS(10,5)+=crLHS2*(crLHS127*r_DN(2,1) + crLHS172 + crLHS208*crLHS91 + crLHS88*r_DN(2,0));
rLHS(10,6)+=crLHS2*(crLHS102*crLHS208 + crLHS129*r_DN(2,1) + crLHS182 + crLHS200 + crLHS99*r_DN(2,0));
rLHS(10,7)+=crLHS209;
rLHS(10,8)+=crLHS210;
rLHS(10,9)+=crLHS2*(crLHS107*r_DN(2,0) + crLHS109*crLHS208 + crLHS133*r_DN(2,1) + crLHS205);
rLHS(10,10)+=crLHS2*(crLHS117*r_DN(2,0) + crLHS119*crLHS208 + crLHS135*r_DN(2,1) + crLHS188*crLHS9 + crLHS204);
rLHS(10,11)+=crLHS210;
rLHS(11,0)+=crLHS43;
rLHS(11,1)+=crLHS183;
rLHS(11,2)+=crLHS184;
rLHS(11,3)+=crLHS43;
rLHS(11,4)+=crLHS147;
rLHS(11,5)+=crLHS185;
rLHS(11,6)+=crLHS186;
rLHS(11,7)+=crLHS147;
rLHS(11,8)+=crLHS189;
rLHS(11,9)+=crLHS191;
rLHS(11,10)+=crLHS192;
rLHS(11,11)+=crLHS189;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

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
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = 1.0/c_p;
const double crLHS6 = gamma - 1.0;
const double crLHS7 = 1.0/crLHS6;
const double crLHS8 = crLHS5*crLHS7*gamma*p_th;
const double crLHS9 = crLHS4*crLHS8;
const double crLHS10 = crLHS2*crLHS9;
const double crLHS11 = crLHS10*tau_u;
const double crLHS12 = crLHS11*(crLHS0 + crLHS1);
const double crLHS13 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS14 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS15 = bdf0*r_N[0];
const double crLHS16 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crLHS17 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crLHS18 = crLHS16*r_DN(0,0) + crLHS17*r_DN(0,1);
const double crLHS19 = crLHS15 + crLHS18;
const double crLHS20 = crLHS9*tau_u;
const double crLHS21 = crLHS19*crLHS20;
const double crLHS22 = crLHS10*(-crLHS13*crLHS14 + crLHS21*r_DN(0,0) + r_DN(0,0)*r_N[0]);
const double crLHS23 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS24 = crLHS10*(-crLHS14*crLHS23 + crLHS21*r_DN(0,1) + r_DN(0,1)*r_N[0]);
const double crLHS25 = r_DN(0,0)*r_DN(1,0);
const double crLHS26 = r_DN(0,1)*r_DN(1,1);
const double crLHS27 = crLHS11*(crLHS25 + crLHS26);
const double crLHS28 = r_DN(1,0)*r_N[0];
const double crLHS29 = crLHS4*r_N[0];
const double crLHS30 = crLHS13*crLHS29;
const double crLHS31 = -crLHS30*r_N[1];
const double crLHS32 = bdf0*r_N[1];
const double crLHS33 = crLHS16*r_DN(1,0) + crLHS17*r_DN(1,1);
const double crLHS34 = crLHS32 + crLHS33;
const double crLHS35 = crLHS20*crLHS34;
const double crLHS36 = crLHS10*(crLHS28 + crLHS31 + crLHS35*r_DN(0,0));
const double crLHS37 = r_DN(1,1)*r_N[0];
const double crLHS38 = crLHS23*crLHS29;
const double crLHS39 = -crLHS38*r_N[1];
const double crLHS40 = crLHS10*(crLHS35*r_DN(0,1) + crLHS37 + crLHS39);
const double crLHS41 = r_DN(0,0)*r_DN(2,0);
const double crLHS42 = r_DN(0,1)*r_DN(2,1);
const double crLHS43 = crLHS11*(crLHS41 + crLHS42);
const double crLHS44 = r_DN(2,0)*r_N[0];
const double crLHS45 = -crLHS30*r_N[2];
const double crLHS46 = bdf0*r_N[2];
const double crLHS47 = crLHS16*r_DN(2,0) + crLHS17*r_DN(2,1);
const double crLHS48 = crLHS46 + crLHS47;
const double crLHS49 = crLHS20*crLHS48;
const double crLHS50 = crLHS10*(crLHS44 + crLHS45 + crLHS49*r_DN(0,0));
const double crLHS51 = r_DN(2,1)*r_N[0];
const double crLHS52 = -crLHS38*r_N[2];
const double crLHS53 = crLHS10*(crLHS49*r_DN(0,1) + crLHS51 + crLHS52);
const double crLHS54 = r_DN(0,0)*r_DN(3,0);
const double crLHS55 = r_DN(0,1)*r_DN(3,1);
const double crLHS56 = crLHS11*(crLHS54 + crLHS55);
const double crLHS57 = r_DN(3,0)*r_N[0];
const double crLHS58 = -crLHS30*r_N[3];
const double crLHS59 = crLHS16*r_DN(3,0) + crLHS17*r_DN(3,1);
const double crLHS60 = bdf0*r_N[3] + crLHS59;
const double crLHS61 = crLHS20*crLHS60;
const double crLHS62 = crLHS10*(crLHS57 + crLHS58 + crLHS61*r_DN(0,0));
const double crLHS63 = r_DN(3,1)*r_N[0];
const double crLHS64 = -crLHS38*r_N[3];
const double crLHS65 = crLHS10*(crLHS61*r_DN(0,1) + crLHS63 + crLHS64);
const double crLHS66 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1);
const double crLHS67 = crLHS13*crLHS16 + crLHS17*crLHS23;
const double crLHS68 = crLHS67*r_N[0];
const double crLHS69 = 1.0/(crLHS3*crLHS3);
const double crLHS70 = crLHS69*tau_u;
const double crLHS71 = crLHS70*crLHS8;
const double crLHS72 = crLHS2*(crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_N[0]*tau_u - crLHS68*crLHS71 - r_N[0]);
const double crLHS73 = crLHS72*r_DN(0,0);
const double crLHS74 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS75 = r_C(0,2)*r_DN(0,0);
const double crLHS76 = crLHS75 + r_C(2,2)*r_DN(0,1);
const double crLHS77 = -crLHS30 + r_DN(0,0);
const double crLHS78 = crLHS9*r_DN(0,0);
const double crLHS79 = crLHS78*tau_c;
const double crLHS80 = bdf0*crLHS8;
const double crLHS81 = crLHS29*crLHS8;
const double crLHS82 = 1.0/(crLHS6*crLHS6)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS83 = crLHS70*crLHS82;
const double crLHS84 = crLHS18*crLHS83;
const double crLHS85 = crLHS19*crLHS82;
const double crLHS86 = crLHS66*crLHS69*tau_u;
const double crLHS87 = crLHS86*r_N[0];
const double crLHS88 = tau_u*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS89 = crLHS68*crLHS88;
const double crLHS90 = crLHS14*crLHS80 + crLHS18*crLHS81 + crLHS19*crLHS84 + crLHS85*crLHS87 - crLHS85*crLHS89;
const double crLHS91 = crLHS75 + r_C(0,1)*r_DN(0,1);
const double crLHS92 = r_C(1,2)*r_DN(0,1);
const double crLHS93 = crLHS92 + r_C(2,2)*r_DN(0,0);
const double crLHS94 = crLHS78*r_DN(0,1);
const double crLHS95 = -crLHS38 + r_DN(0,1);
const double crLHS96 = r_DN(0,0)*r_N[1];
const double crLHS97 = crLHS67*crLHS71;
const double crLHS98 = crLHS2*(crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(1,0)*tau_u - crLHS28*crLHS97 + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,0)*r_N[0]*tau_u - crLHS96);
const double crLHS99 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS100 = r_C(0,2)*r_DN(1,0);
const double crLHS101 = crLHS100 + r_C(2,2)*r_DN(1,1);
const double crLHS102 = crLHS4*r_N[1];
const double crLHS103 = crLHS102*crLHS13;
const double crLHS104 = -crLHS103 + r_DN(1,0);
const double crLHS105 = crLHS102*crLHS8;
const double crLHS106 = crLHS105*crLHS15;
const double crLHS107 = crLHS106 + crLHS25*crLHS9;
const double crLHS108 = crLHS34*crLHS82;
const double crLHS109 = crLHS108*crLHS87 - crLHS108*crLHS89 + crLHS33*crLHS81 + crLHS34*crLHS84;
const double crLHS110 = crLHS100 + r_C(0,1)*r_DN(1,1);
const double crLHS111 = r_C(1,2)*r_DN(1,1);
const double crLHS112 = crLHS111 + r_C(2,2)*r_DN(1,0);
const double crLHS113 = crLHS78*r_DN(1,1);
const double crLHS114 = crLHS102*crLHS23;
const double crLHS115 = -crLHS114 + r_DN(1,1);
const double crLHS116 = r_DN(0,0)*r_N[2];
const double crLHS117 = crLHS2*(-crLHS116 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,0)*r_N[0]*tau_u - crLHS44*crLHS97);
const double crLHS118 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS119 = r_C(0,2)*r_DN(2,0);
const double crLHS120 = crLHS119 + r_C(2,2)*r_DN(2,1);
const double crLHS121 = crLHS4*r_N[2];
const double crLHS122 = crLHS121*crLHS13;
const double crLHS123 = -crLHS122 + r_DN(2,0);
const double crLHS124 = crLHS121*crLHS8;
const double crLHS125 = crLHS124*crLHS15;
const double crLHS126 = crLHS125 + crLHS41*crLHS9;
const double crLHS127 = crLHS48*crLHS82;
const double crLHS128 = crLHS127*crLHS87 - crLHS127*crLHS89 + crLHS47*crLHS81 + crLHS48*crLHS84;
const double crLHS129 = crLHS119 + r_C(0,1)*r_DN(2,1);
const double crLHS130 = r_C(1,2)*r_DN(2,1);
const double crLHS131 = crLHS130 + r_C(2,2)*r_DN(2,0);
const double crLHS132 = crLHS78*r_DN(2,1);
const double crLHS133 = crLHS121*crLHS23;
const double crLHS134 = -crLHS133 + r_DN(2,1);
const double crLHS135 = r_DN(0,0)*r_N[3];
const double crLHS136 = crLHS2*(-crLHS135 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(3,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,0)*r_N[0]*tau_u - crLHS57*crLHS97);
const double crLHS137 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS138 = r_C(0,2)*r_DN(3,0);
const double crLHS139 = crLHS138 + r_C(2,2)*r_DN(3,1);
const double crLHS140 = crLHS4*r_N[3];
const double crLHS141 = -crLHS13*crLHS140 + r_DN(3,0);
const double crLHS142 = crLHS140*crLHS8;
const double crLHS143 = crLHS142*crLHS15;
const double crLHS144 = crLHS143 + crLHS54*crLHS9;
const double crLHS145 = crLHS60*crLHS82;
const double crLHS146 = crLHS145*crLHS87 - crLHS145*crLHS89 + crLHS59*crLHS81 + crLHS60*crLHS84;
const double crLHS147 = crLHS138 + r_C(0,1)*r_DN(3,1);
const double crLHS148 = r_C(1,2)*r_DN(3,1);
const double crLHS149 = crLHS148 + r_C(2,2)*r_DN(3,0);
const double crLHS150 = crLHS78*r_DN(3,1);
const double crLHS151 = -crLHS140*crLHS23 + r_DN(3,1);
const double crLHS152 = crLHS72*r_DN(0,1);
const double crLHS153 = crLHS92 + r_C(0,1)*r_DN(0,0);
const double crLHS154 = crLHS9*r_DN(0,1);
const double crLHS155 = crLHS154*tau_c;
const double crLHS156 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS157 = r_DN(0,1)*r_N[1];
const double crLHS158 = crLHS2*(-crLHS157 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(1,1)*tau_u - crLHS37*crLHS97 + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,1)*r_N[0]*tau_u);
const double crLHS159 = crLHS111 + r_C(0,1)*r_DN(1,0);
const double crLHS160 = crLHS154*r_DN(1,0);
const double crLHS161 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS162 = crLHS106 + crLHS26*crLHS9;
const double crLHS163 = r_DN(0,1)*r_N[2];
const double crLHS164 = crLHS2*(-crLHS163 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,1)*r_N[0]*tau_u - crLHS51*crLHS97);
const double crLHS165 = crLHS130 + r_C(0,1)*r_DN(2,0);
const double crLHS166 = crLHS154*r_DN(2,0);
const double crLHS167 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS168 = crLHS125 + crLHS42*crLHS9;
const double crLHS169 = r_DN(0,1)*r_N[3];
const double crLHS170 = crLHS2*(-crLHS169 + crLHS18*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(3,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,1)*r_N[0]*tau_u - crLHS63*crLHS97);
const double crLHS171 = crLHS148 + r_C(0,1)*r_DN(3,0);
const double crLHS172 = crLHS154*r_DN(3,0);
const double crLHS173 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS174 = crLHS143 + crLHS55*crLHS9;
const double crLHS175 = crLHS10*(crLHS21*r_DN(1,0) + crLHS31 + crLHS96);
const double crLHS176 = crLHS10*(crLHS157 + crLHS21*r_DN(1,1) + crLHS39);
const double crLHS177 = r_DN(1,0)*r_DN(1,0);
const double crLHS178 = r_DN(1,1)*r_DN(1,1);
const double crLHS179 = crLHS11*(crLHS177 + crLHS178);
const double crLHS180 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS181 = crLHS10*(-crLHS13*crLHS180 + crLHS35*r_DN(1,0) + r_DN(1,0)*r_N[1]);
const double crLHS182 = crLHS10*(-crLHS180*crLHS23 + crLHS35*r_DN(1,1) + r_DN(1,1)*r_N[1]);
const double crLHS183 = r_DN(1,0)*r_DN(2,0);
const double crLHS184 = r_DN(1,1)*r_DN(2,1);
const double crLHS185 = crLHS11*(crLHS183 + crLHS184);
const double crLHS186 = r_DN(2,0)*r_N[1];
const double crLHS187 = -crLHS103*r_N[2];
const double crLHS188 = crLHS10*(crLHS186 + crLHS187 + crLHS49*r_DN(1,0));
const double crLHS189 = r_DN(2,1)*r_N[1];
const double crLHS190 = -crLHS114*r_N[2];
const double crLHS191 = crLHS10*(crLHS189 + crLHS190 + crLHS49*r_DN(1,1));
const double crLHS192 = r_DN(1,0)*r_DN(3,0);
const double crLHS193 = r_DN(1,1)*r_DN(3,1);
const double crLHS194 = crLHS11*(crLHS192 + crLHS193);
const double crLHS195 = r_DN(3,0)*r_N[1];
const double crLHS196 = -crLHS103*r_N[3];
const double crLHS197 = crLHS10*(crLHS195 + crLHS196 + crLHS61*r_DN(1,0));
const double crLHS198 = r_DN(3,1)*r_N[1];
const double crLHS199 = -crLHS114*r_N[3];
const double crLHS200 = crLHS10*(crLHS198 + crLHS199 + crLHS61*r_DN(1,1));
const double crLHS201 = crLHS2*(-crLHS28 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(0,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,0)*r_N[1]*tau_u - crLHS96*crLHS97);
const double crLHS202 = crLHS9*r_DN(1,0);
const double crLHS203 = crLHS202*tau_c;
const double crLHS204 = crLHS33*crLHS83;
const double crLHS205 = crLHS85*r_N[1];
const double crLHS206 = crLHS67*crLHS88;
const double crLHS207 = crLHS105*crLHS18 + crLHS19*crLHS204 - crLHS205*crLHS206 + crLHS205*crLHS86;
const double crLHS208 = crLHS2*(crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_N[1]*tau_u - crLHS97*r_N[1] - r_N[1]);
const double crLHS209 = crLHS208*r_DN(1,0);
const double crLHS210 = crLHS108*r_N[1];
const double crLHS211 = crLHS105*crLHS33 + crLHS180*crLHS80 + crLHS204*crLHS34 - crLHS206*crLHS210 + crLHS210*crLHS86;
const double crLHS212 = crLHS202*r_DN(1,1);
const double crLHS213 = r_DN(1,0)*r_N[2];
const double crLHS214 = crLHS2*(-crLHS186*crLHS97 - crLHS213 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,0)*r_N[1]*tau_u);
const double crLHS215 = crLHS124*crLHS32;
const double crLHS216 = crLHS183*crLHS9 + crLHS215;
const double crLHS217 = crLHS127*r_N[1];
const double crLHS218 = crLHS105*crLHS47 + crLHS204*crLHS48 - crLHS206*crLHS217 + crLHS217*crLHS86;
const double crLHS219 = crLHS202*r_DN(2,1);
const double crLHS220 = r_DN(1,0)*r_N[3];
const double crLHS221 = crLHS2*(-crLHS195*crLHS97 - crLHS220 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(3,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,0)*r_N[1]*tau_u);
const double crLHS222 = crLHS142*crLHS32;
const double crLHS223 = crLHS192*crLHS9 + crLHS222;
const double crLHS224 = crLHS145*r_N[1];
const double crLHS225 = crLHS105*crLHS59 + crLHS204*crLHS60 - crLHS206*crLHS224 + crLHS224*crLHS86;
const double crLHS226 = crLHS202*r_DN(3,1);
const double crLHS227 = crLHS2*(-crLHS157*crLHS97 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(0,1)*tau_u - crLHS37 + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,1)*r_N[1]*tau_u);
const double crLHS228 = crLHS9*r_DN(1,1);
const double crLHS229 = crLHS228*tau_c;
const double crLHS230 = crLHS208*r_DN(1,1);
const double crLHS231 = r_DN(1,1)*r_N[2];
const double crLHS232 = crLHS2*(-crLHS189*crLHS97 - crLHS231 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(2,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,1)*r_N[1]*tau_u);
const double crLHS233 = crLHS228*r_DN(2,0);
const double crLHS234 = crLHS184*crLHS9 + crLHS215;
const double crLHS235 = r_DN(1,1)*r_N[3];
const double crLHS236 = crLHS2*(-crLHS198*crLHS97 - crLHS235 + crLHS33*crLHS4*crLHS5*crLHS7*gamma*p_th*r_DN(3,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,1)*r_N[1]*tau_u);
const double crLHS237 = crLHS228*r_DN(3,0);
const double crLHS238 = crLHS193*crLHS9 + crLHS222;
const double crLHS239 = crLHS10*(crLHS116 + crLHS21*r_DN(2,0) + crLHS45);
const double crLHS240 = crLHS10*(crLHS163 + crLHS21*r_DN(2,1) + crLHS52);
const double crLHS241 = crLHS10*(crLHS187 + crLHS213 + crLHS35*r_DN(2,0));
const double crLHS242 = crLHS10*(crLHS190 + crLHS231 + crLHS35*r_DN(2,1));
const double crLHS243 = r_DN(2,0)*r_DN(2,0);
const double crLHS244 = r_DN(2,1)*r_DN(2,1);
const double crLHS245 = crLHS11*(crLHS243 + crLHS244);
const double crLHS246 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS247 = crLHS10*(-crLHS13*crLHS246 + crLHS49*r_DN(2,0) + r_DN(2,0)*r_N[2]);
const double crLHS248 = crLHS10*(-crLHS23*crLHS246 + crLHS49*r_DN(2,1) + r_DN(2,1)*r_N[2]);
const double crLHS249 = r_DN(2,0)*r_DN(3,0);
const double crLHS250 = r_DN(2,1)*r_DN(3,1);
const double crLHS251 = crLHS11*(crLHS249 + crLHS250);
const double crLHS252 = r_DN(3,0)*r_N[2];
const double crLHS253 = -crLHS122*r_N[3];
const double crLHS254 = crLHS10*(crLHS252 + crLHS253 + crLHS61*r_DN(2,0));
const double crLHS255 = r_DN(3,1)*r_N[2];
const double crLHS256 = -crLHS133*r_N[3];
const double crLHS257 = crLHS10*(crLHS255 + crLHS256 + crLHS61*r_DN(2,1));
const double crLHS258 = crLHS2*(-crLHS116*crLHS97 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(0,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,0)*r_N[2]*tau_u - crLHS44);
const double crLHS259 = crLHS9*r_DN(2,0);
const double crLHS260 = crLHS259*tau_c;
const double crLHS261 = crLHS47*crLHS83;
const double crLHS262 = crLHS85*r_N[2];
const double crLHS263 = crLHS124*crLHS18 + crLHS19*crLHS261 - crLHS206*crLHS262 + crLHS262*crLHS86;
const double crLHS264 = crLHS2*(-crLHS186 - crLHS213*crLHS97 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(1,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,0)*r_N[2]*tau_u);
const double crLHS265 = crLHS108*r_N[2];
const double crLHS266 = crLHS124*crLHS33 - crLHS206*crLHS265 + crLHS261*crLHS34 + crLHS265*crLHS86;
const double crLHS267 = crLHS2*(crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_N[2]*tau_u - crLHS97*r_N[2] - r_N[2]);
const double crLHS268 = crLHS267*r_DN(2,0);
const double crLHS269 = crLHS127*r_N[2];
const double crLHS270 = crLHS124*crLHS47 - crLHS206*crLHS269 + crLHS246*crLHS80 + crLHS261*crLHS48 + crLHS269*crLHS86;
const double crLHS271 = crLHS259*r_DN(2,1);
const double crLHS272 = r_DN(2,0)*r_N[3];
const double crLHS273 = crLHS2*(-crLHS252*crLHS97 - crLHS272 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(3,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,0)*r_N[2]*tau_u);
const double crLHS274 = crLHS142*crLHS46;
const double crLHS275 = crLHS249*crLHS9 + crLHS274;
const double crLHS276 = crLHS145*r_N[2];
const double crLHS277 = crLHS124*crLHS59 - crLHS206*crLHS276 + crLHS261*crLHS60 + crLHS276*crLHS86;
const double crLHS278 = crLHS259*r_DN(3,1);
const double crLHS279 = crLHS2*(-crLHS163*crLHS97 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(0,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,1)*r_N[2]*tau_u - crLHS51);
const double crLHS280 = crLHS77*tau_c;
const double crLHS281 = crLHS9*r_DN(2,1);
const double crLHS282 = crLHS281*tau_c;
const double crLHS283 = crLHS2*(-crLHS189 - crLHS231*crLHS97 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(1,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,1)*r_N[2]*tau_u);
const double crLHS284 = crLHS267*r_DN(2,1);
const double crLHS285 = r_DN(2,1)*r_N[3];
const double crLHS286 = crLHS2*(-crLHS255*crLHS97 - crLHS285 + crLHS4*crLHS47*crLHS5*crLHS7*gamma*p_th*r_DN(3,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(3,1)*r_N[2]*tau_u);
const double crLHS287 = crLHS9*r_DN(3,0);
const double crLHS288 = crLHS287*r_DN(2,1);
const double crLHS289 = crLHS250*crLHS9 + crLHS274;
const double crLHS290 = crLHS10*(crLHS135 + crLHS21*r_DN(3,0) + crLHS58);
const double crLHS291 = crLHS10*(crLHS169 + crLHS21*r_DN(3,1) + crLHS64);
const double crLHS292 = crLHS10*(crLHS196 + crLHS220 + crLHS35*r_DN(3,0));
const double crLHS293 = crLHS10*(crLHS199 + crLHS235 + crLHS35*r_DN(3,1));
const double crLHS294 = crLHS10*(crLHS253 + crLHS272 + crLHS49*r_DN(3,0));
const double crLHS295 = crLHS10*(crLHS256 + crLHS285 + crLHS49*r_DN(3,1));
const double crLHS296 = r_DN(3,0)*r_DN(3,0);
const double crLHS297 = r_DN(3,1)*r_DN(3,1);
const double crLHS298 = crLHS11*(crLHS296 + crLHS297);
const double crLHS299 = crLHS4*(r_N[3]*r_N[3]);
const double crLHS300 = crLHS10*(-crLHS13*crLHS299 + crLHS61*r_DN(3,0) + r_DN(3,0)*r_N[3]);
const double crLHS301 = crLHS10*(-crLHS23*crLHS299 + crLHS61*r_DN(3,1) + r_DN(3,1)*r_N[3]);
const double crLHS302 = crLHS2*(-crLHS135*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(0,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,0)*r_N[3]*tau_u - crLHS57);
const double crLHS303 = crLHS59*crLHS83;
const double crLHS304 = crLHS85*r_N[3];
const double crLHS305 = crLHS142*crLHS18 + crLHS19*crLHS303 - crLHS206*crLHS304 + crLHS304*crLHS86;
const double crLHS306 = crLHS287*tau_c;
const double crLHS307 = crLHS2*(-crLHS195 - crLHS220*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(1,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,0)*r_N[3]*tau_u);
const double crLHS308 = crLHS108*r_N[3];
const double crLHS309 = crLHS142*crLHS33 - crLHS206*crLHS308 + crLHS303*crLHS34 + crLHS308*crLHS86;
const double crLHS310 = crLHS2*(-crLHS252 - crLHS272*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(2,0)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,0)*r_N[3]*tau_u);
const double crLHS311 = crLHS127*r_N[3];
const double crLHS312 = crLHS142*crLHS47 - crLHS206*crLHS311 + crLHS303*crLHS48 + crLHS311*crLHS86;
const double crLHS313 = crLHS2*(crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_N[3]*tau_u - crLHS97*r_N[3] - r_N[3]);
const double crLHS314 = crLHS313*r_DN(3,0);
const double crLHS315 = crLHS145*r_N[3];
const double crLHS316 = crLHS142*crLHS59 - crLHS206*crLHS315 + crLHS299*crLHS80 + crLHS303*crLHS60 + crLHS315*crLHS86;
const double crLHS317 = crLHS287*r_DN(3,1);
const double crLHS318 = crLHS2*(-crLHS169*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(0,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(0,1)*r_N[3]*tau_u - crLHS63);
const double crLHS319 = crLHS9*r_DN(3,1);
const double crLHS320 = crLHS319*tau_c;
const double crLHS321 = crLHS2*(-crLHS198 - crLHS235*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(1,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(1,1)*r_N[3]*tau_u);
const double crLHS322 = crLHS2*(-crLHS255 - crLHS285*crLHS97 + crLHS4*crLHS5*crLHS59*crLHS7*gamma*p_th*r_DN(2,1)*tau_u + crLHS4*crLHS5*crLHS66*crLHS7*gamma*p_th*r_DN(2,1)*r_N[3]*tau_u);
const double crLHS323 = crLHS313*r_DN(3,1);
rLHS(0,0)+=crLHS12;
rLHS(0,1)+=crLHS22;
rLHS(0,2)+=crLHS24;
rLHS(0,3)+=crLHS12;
rLHS(0,4)+=crLHS27;
rLHS(0,5)+=crLHS36;
rLHS(0,6)+=crLHS40;
rLHS(0,7)+=crLHS27;
rLHS(0,8)+=crLHS43;
rLHS(0,9)+=crLHS50;
rLHS(0,10)+=crLHS53;
rLHS(0,11)+=crLHS43;
rLHS(0,12)+=crLHS56;
rLHS(0,13)+=crLHS62;
rLHS(0,14)+=crLHS65;
rLHS(0,15)+=crLHS56;
rLHS(1,0)+=crLHS73;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS9 + crLHS74*r_DN(0,0) + crLHS76*r_DN(0,1) + crLHS77*crLHS79 + crLHS90);
rLHS(1,2)+=crLHS2*(crLHS79*crLHS95 + crLHS91*r_DN(0,0) + crLHS93*r_DN(0,1) + crLHS94);
rLHS(1,3)+=crLHS73;
rLHS(1,4)+=crLHS98;
rLHS(1,5)+=crLHS2*(crLHS101*r_DN(0,1) + crLHS104*crLHS79 + crLHS107 + crLHS109 + crLHS99*r_DN(0,0));
rLHS(1,6)+=crLHS2*(crLHS110*r_DN(0,0) + crLHS112*r_DN(0,1) + crLHS113 + crLHS115*crLHS79);
rLHS(1,7)+=crLHS98;
rLHS(1,8)+=crLHS117;
rLHS(1,9)+=crLHS2*(crLHS118*r_DN(0,0) + crLHS120*r_DN(0,1) + crLHS123*crLHS79 + crLHS126 + crLHS128);
rLHS(1,10)+=crLHS2*(crLHS129*r_DN(0,0) + crLHS131*r_DN(0,1) + crLHS132 + crLHS134*crLHS79);
rLHS(1,11)+=crLHS117;
rLHS(1,12)+=crLHS136;
rLHS(1,13)+=crLHS2*(crLHS137*r_DN(0,0) + crLHS139*r_DN(0,1) + crLHS141*crLHS79 + crLHS144 + crLHS146);
rLHS(1,14)+=crLHS2*(crLHS147*r_DN(0,0) + crLHS149*r_DN(0,1) + crLHS150 + crLHS151*crLHS79);
rLHS(1,15)+=crLHS136;
rLHS(2,0)+=crLHS152;
rLHS(2,1)+=crLHS2*(crLHS153*r_DN(0,1) + crLHS155*crLHS77 + crLHS76*r_DN(0,0) + crLHS94);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS9 + crLHS155*crLHS95 + crLHS156*r_DN(0,1) + crLHS90 + crLHS93*r_DN(0,0));
rLHS(2,3)+=crLHS152;
rLHS(2,4)+=crLHS158;
rLHS(2,5)+=crLHS2*(crLHS101*r_DN(0,0) + crLHS104*crLHS155 + crLHS159*r_DN(0,1) + crLHS160);
rLHS(2,6)+=crLHS2*(crLHS109 + crLHS112*r_DN(0,0) + crLHS115*crLHS155 + crLHS161*r_DN(0,1) + crLHS162);
rLHS(2,7)+=crLHS158;
rLHS(2,8)+=crLHS164;
rLHS(2,9)+=crLHS2*(crLHS120*r_DN(0,0) + crLHS123*crLHS155 + crLHS165*r_DN(0,1) + crLHS166);
rLHS(2,10)+=crLHS2*(crLHS128 + crLHS131*r_DN(0,0) + crLHS134*crLHS155 + crLHS167*r_DN(0,1) + crLHS168);
rLHS(2,11)+=crLHS164;
rLHS(2,12)+=crLHS170;
rLHS(2,13)+=crLHS2*(crLHS139*r_DN(0,0) + crLHS141*crLHS155 + crLHS171*r_DN(0,1) + crLHS172);
rLHS(2,14)+=crLHS2*(crLHS146 + crLHS149*r_DN(0,0) + crLHS151*crLHS155 + crLHS173*r_DN(0,1) + crLHS174);
rLHS(2,15)+=crLHS170;
rLHS(3,0)+=crLHS12;
rLHS(3,1)+=crLHS22;
rLHS(3,2)+=crLHS24;
rLHS(3,3)+=crLHS12;
rLHS(3,4)+=crLHS27;
rLHS(3,5)+=crLHS36;
rLHS(3,6)+=crLHS40;
rLHS(3,7)+=crLHS27;
rLHS(3,8)+=crLHS43;
rLHS(3,9)+=crLHS50;
rLHS(3,10)+=crLHS53;
rLHS(3,11)+=crLHS43;
rLHS(3,12)+=crLHS56;
rLHS(3,13)+=crLHS62;
rLHS(3,14)+=crLHS65;
rLHS(3,15)+=crLHS56;
rLHS(4,0)+=crLHS27;
rLHS(4,1)+=crLHS175;
rLHS(4,2)+=crLHS176;
rLHS(4,3)+=crLHS27;
rLHS(4,4)+=crLHS179;
rLHS(4,5)+=crLHS181;
rLHS(4,6)+=crLHS182;
rLHS(4,7)+=crLHS179;
rLHS(4,8)+=crLHS185;
rLHS(4,9)+=crLHS188;
rLHS(4,10)+=crLHS191;
rLHS(4,11)+=crLHS185;
rLHS(4,12)+=crLHS194;
rLHS(4,13)+=crLHS197;
rLHS(4,14)+=crLHS200;
rLHS(4,15)+=crLHS194;
rLHS(5,0)+=crLHS201;
rLHS(5,1)+=crLHS2*(crLHS107 + crLHS203*crLHS77 + crLHS207 + crLHS74*r_DN(1,0) + crLHS76*r_DN(1,1));
rLHS(5,2)+=crLHS2*(crLHS160 + crLHS203*crLHS95 + crLHS91*r_DN(1,0) + crLHS93*r_DN(1,1));
rLHS(5,3)+=crLHS201;
rLHS(5,4)+=crLHS209;
rLHS(5,5)+=crLHS2*(crLHS101*r_DN(1,1) + crLHS104*crLHS203 + crLHS177*crLHS9 + crLHS211 + crLHS99*r_DN(1,0));
rLHS(5,6)+=crLHS2*(crLHS110*r_DN(1,0) + crLHS112*r_DN(1,1) + crLHS115*crLHS203 + crLHS212);
rLHS(5,7)+=crLHS209;
rLHS(5,8)+=crLHS214;
rLHS(5,9)+=crLHS2*(crLHS118*r_DN(1,0) + crLHS120*r_DN(1,1) + crLHS123*crLHS203 + crLHS216 + crLHS218);
rLHS(5,10)+=crLHS2*(crLHS129*r_DN(1,0) + crLHS131*r_DN(1,1) + crLHS134*crLHS203 + crLHS219);
rLHS(5,11)+=crLHS214;
rLHS(5,12)+=crLHS221;
rLHS(5,13)+=crLHS2*(crLHS137*r_DN(1,0) + crLHS139*r_DN(1,1) + crLHS141*crLHS203 + crLHS223 + crLHS225);
rLHS(5,14)+=crLHS2*(crLHS147*r_DN(1,0) + crLHS149*r_DN(1,1) + crLHS151*crLHS203 + crLHS226);
rLHS(5,15)+=crLHS221;
rLHS(6,0)+=crLHS227;
rLHS(6,1)+=crLHS2*(crLHS113 + crLHS153*r_DN(1,1) + crLHS229*crLHS77 + crLHS76*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS156*r_DN(1,1) + crLHS162 + crLHS207 + crLHS229*crLHS95 + crLHS93*r_DN(1,0));
rLHS(6,3)+=crLHS227;
rLHS(6,4)+=crLHS230;
rLHS(6,5)+=crLHS2*(crLHS101*r_DN(1,0) + crLHS104*crLHS229 + crLHS159*r_DN(1,1) + crLHS212);
rLHS(6,6)+=crLHS2*(crLHS112*r_DN(1,0) + crLHS115*crLHS229 + crLHS161*r_DN(1,1) + crLHS178*crLHS9 + crLHS211);
rLHS(6,7)+=crLHS230;
rLHS(6,8)+=crLHS232;
rLHS(6,9)+=crLHS2*(crLHS120*r_DN(1,0) + crLHS123*crLHS229 + crLHS165*r_DN(1,1) + crLHS233);
rLHS(6,10)+=crLHS2*(crLHS131*r_DN(1,0) + crLHS134*crLHS229 + crLHS167*r_DN(1,1) + crLHS218 + crLHS234);
rLHS(6,11)+=crLHS232;
rLHS(6,12)+=crLHS236;
rLHS(6,13)+=crLHS2*(crLHS139*r_DN(1,0) + crLHS141*crLHS229 + crLHS171*r_DN(1,1) + crLHS237);
rLHS(6,14)+=crLHS2*(crLHS149*r_DN(1,0) + crLHS151*crLHS229 + crLHS173*r_DN(1,1) + crLHS225 + crLHS238);
rLHS(6,15)+=crLHS236;
rLHS(7,0)+=crLHS27;
rLHS(7,1)+=crLHS175;
rLHS(7,2)+=crLHS176;
rLHS(7,3)+=crLHS27;
rLHS(7,4)+=crLHS179;
rLHS(7,5)+=crLHS181;
rLHS(7,6)+=crLHS182;
rLHS(7,7)+=crLHS179;
rLHS(7,8)+=crLHS185;
rLHS(7,9)+=crLHS188;
rLHS(7,10)+=crLHS191;
rLHS(7,11)+=crLHS185;
rLHS(7,12)+=crLHS194;
rLHS(7,13)+=crLHS197;
rLHS(7,14)+=crLHS200;
rLHS(7,15)+=crLHS194;
rLHS(8,0)+=crLHS43;
rLHS(8,1)+=crLHS239;
rLHS(8,2)+=crLHS240;
rLHS(8,3)+=crLHS43;
rLHS(8,4)+=crLHS185;
rLHS(8,5)+=crLHS241;
rLHS(8,6)+=crLHS242;
rLHS(8,7)+=crLHS185;
rLHS(8,8)+=crLHS245;
rLHS(8,9)+=crLHS247;
rLHS(8,10)+=crLHS248;
rLHS(8,11)+=crLHS245;
rLHS(8,12)+=crLHS251;
rLHS(8,13)+=crLHS254;
rLHS(8,14)+=crLHS257;
rLHS(8,15)+=crLHS251;
rLHS(9,0)+=crLHS258;
rLHS(9,1)+=crLHS2*(crLHS126 + crLHS260*crLHS77 + crLHS263 + crLHS74*r_DN(2,0) + crLHS76*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS166 + crLHS260*crLHS95 + crLHS91*r_DN(2,0) + crLHS93*r_DN(2,1));
rLHS(9,3)+=crLHS258;
rLHS(9,4)+=crLHS264;
rLHS(9,5)+=crLHS2*(crLHS101*r_DN(2,1) + crLHS104*crLHS260 + crLHS216 + crLHS266 + crLHS99*r_DN(2,0));
rLHS(9,6)+=crLHS2*(crLHS110*r_DN(2,0) + crLHS112*r_DN(2,1) + crLHS115*crLHS260 + crLHS233);
rLHS(9,7)+=crLHS264;
rLHS(9,8)+=crLHS268;
rLHS(9,9)+=crLHS2*(crLHS118*r_DN(2,0) + crLHS120*r_DN(2,1) + crLHS123*crLHS260 + crLHS243*crLHS9 + crLHS270);
rLHS(9,10)+=crLHS2*(crLHS129*r_DN(2,0) + crLHS131*r_DN(2,1) + crLHS134*crLHS260 + crLHS271);
rLHS(9,11)+=crLHS268;
rLHS(9,12)+=crLHS273;
rLHS(9,13)+=crLHS2*(crLHS137*r_DN(2,0) + crLHS139*r_DN(2,1) + crLHS141*crLHS260 + crLHS275 + crLHS277);
rLHS(9,14)+=crLHS2*(crLHS147*r_DN(2,0) + crLHS149*r_DN(2,1) + crLHS151*crLHS260 + crLHS278);
rLHS(9,15)+=crLHS273;
rLHS(10,0)+=crLHS279;
rLHS(10,1)+=crLHS2*(crLHS132 + crLHS153*r_DN(2,1) + crLHS280*crLHS281 + crLHS76*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS156*r_DN(2,1) + crLHS168 + crLHS263 + crLHS282*crLHS95 + crLHS93*r_DN(2,0));
rLHS(10,3)+=crLHS279;
rLHS(10,4)+=crLHS283;
rLHS(10,5)+=crLHS2*(crLHS101*r_DN(2,0) + crLHS104*crLHS282 + crLHS159*r_DN(2,1) + crLHS219);
rLHS(10,6)+=crLHS2*(crLHS112*r_DN(2,0) + crLHS115*crLHS282 + crLHS161*r_DN(2,1) + crLHS234 + crLHS266);
rLHS(10,7)+=crLHS283;
rLHS(10,8)+=crLHS284;
rLHS(10,9)+=crLHS2*(crLHS120*r_DN(2,0) + crLHS123*crLHS282 + crLHS165*r_DN(2,1) + crLHS271);
rLHS(10,10)+=crLHS2*(crLHS131*r_DN(2,0) + crLHS134*crLHS282 + crLHS167*r_DN(2,1) + crLHS244*crLHS9 + crLHS270);
rLHS(10,11)+=crLHS284;
rLHS(10,12)+=crLHS286;
rLHS(10,13)+=crLHS2*(crLHS139*r_DN(2,0) + crLHS141*crLHS282 + crLHS171*r_DN(2,1) + crLHS288);
rLHS(10,14)+=crLHS2*(crLHS149*r_DN(2,0) + crLHS151*crLHS282 + crLHS173*r_DN(2,1) + crLHS277 + crLHS289);
rLHS(10,15)+=crLHS286;
rLHS(11,0)+=crLHS43;
rLHS(11,1)+=crLHS239;
rLHS(11,2)+=crLHS240;
rLHS(11,3)+=crLHS43;
rLHS(11,4)+=crLHS185;
rLHS(11,5)+=crLHS241;
rLHS(11,6)+=crLHS242;
rLHS(11,7)+=crLHS185;
rLHS(11,8)+=crLHS245;
rLHS(11,9)+=crLHS247;
rLHS(11,10)+=crLHS248;
rLHS(11,11)+=crLHS245;
rLHS(11,12)+=crLHS251;
rLHS(11,13)+=crLHS254;
rLHS(11,14)+=crLHS257;
rLHS(11,15)+=crLHS251;
rLHS(12,0)+=crLHS56;
rLHS(12,1)+=crLHS290;
rLHS(12,2)+=crLHS291;
rLHS(12,3)+=crLHS56;
rLHS(12,4)+=crLHS194;
rLHS(12,5)+=crLHS292;
rLHS(12,6)+=crLHS293;
rLHS(12,7)+=crLHS194;
rLHS(12,8)+=crLHS251;
rLHS(12,9)+=crLHS294;
rLHS(12,10)+=crLHS295;
rLHS(12,11)+=crLHS251;
rLHS(12,12)+=crLHS298;
rLHS(12,13)+=crLHS300;
rLHS(12,14)+=crLHS301;
rLHS(12,15)+=crLHS298;
rLHS(13,0)+=crLHS302;
rLHS(13,1)+=crLHS2*(crLHS144 + crLHS280*crLHS287 + crLHS305 + crLHS74*r_DN(3,0) + crLHS76*r_DN(3,1));
rLHS(13,2)+=crLHS2*(crLHS172 + crLHS306*crLHS95 + crLHS91*r_DN(3,0) + crLHS93*r_DN(3,1));
rLHS(13,3)+=crLHS302;
rLHS(13,4)+=crLHS307;
rLHS(13,5)+=crLHS2*(crLHS101*r_DN(3,1) + crLHS104*crLHS306 + crLHS223 + crLHS309 + crLHS99*r_DN(3,0));
rLHS(13,6)+=crLHS2*(crLHS110*r_DN(3,0) + crLHS112*r_DN(3,1) + crLHS115*crLHS306 + crLHS237);
rLHS(13,7)+=crLHS307;
rLHS(13,8)+=crLHS310;
rLHS(13,9)+=crLHS2*(crLHS118*r_DN(3,0) + crLHS120*r_DN(3,1) + crLHS123*crLHS306 + crLHS275 + crLHS312);
rLHS(13,10)+=crLHS2*(crLHS129*r_DN(3,0) + crLHS131*r_DN(3,1) + crLHS134*crLHS306 + crLHS288);
rLHS(13,11)+=crLHS310;
rLHS(13,12)+=crLHS314;
rLHS(13,13)+=crLHS2*(crLHS137*r_DN(3,0) + crLHS139*r_DN(3,1) + crLHS141*crLHS306 + crLHS296*crLHS9 + crLHS316);
rLHS(13,14)+=crLHS2*(crLHS147*r_DN(3,0) + crLHS149*r_DN(3,1) + crLHS151*crLHS306 + crLHS317);
rLHS(13,15)+=crLHS314;
rLHS(14,0)+=crLHS318;
rLHS(14,1)+=crLHS2*(crLHS150 + crLHS153*r_DN(3,1) + crLHS280*crLHS319 + crLHS76*r_DN(3,0));
rLHS(14,2)+=crLHS2*(crLHS156*r_DN(3,1) + crLHS174 + crLHS305 + crLHS320*crLHS95 + crLHS93*r_DN(3,0));
rLHS(14,3)+=crLHS318;
rLHS(14,4)+=crLHS321;
rLHS(14,5)+=crLHS2*(crLHS101*r_DN(3,0) + crLHS104*crLHS320 + crLHS159*r_DN(3,1) + crLHS226);
rLHS(14,6)+=crLHS2*(crLHS112*r_DN(3,0) + crLHS115*crLHS320 + crLHS161*r_DN(3,1) + crLHS238 + crLHS309);
rLHS(14,7)+=crLHS321;
rLHS(14,8)+=crLHS322;
rLHS(14,9)+=crLHS2*(crLHS120*r_DN(3,0) + crLHS123*crLHS320 + crLHS165*r_DN(3,1) + crLHS278);
rLHS(14,10)+=crLHS2*(crLHS131*r_DN(3,0) + crLHS134*crLHS320 + crLHS167*r_DN(3,1) + crLHS289 + crLHS312);
rLHS(14,11)+=crLHS322;
rLHS(14,12)+=crLHS323;
rLHS(14,13)+=crLHS2*(crLHS139*r_DN(3,0) + crLHS141*crLHS320 + crLHS171*r_DN(3,1) + crLHS317);
rLHS(14,14)+=crLHS2*(crLHS149*r_DN(3,0) + crLHS151*crLHS320 + crLHS173*r_DN(3,1) + crLHS297*crLHS9 + crLHS316);
rLHS(14,15)+=crLHS323;
rLHS(15,0)+=crLHS56;
rLHS(15,1)+=crLHS290;
rLHS(15,2)+=crLHS291;
rLHS(15,3)+=crLHS56;
rLHS(15,4)+=crLHS194;
rLHS(15,5)+=crLHS292;
rLHS(15,6)+=crLHS293;
rLHS(15,7)+=crLHS194;
rLHS(15,8)+=crLHS251;
rLHS(15,9)+=crLHS294;
rLHS(15,10)+=crLHS295;
rLHS(15,11)+=crLHS251;
rLHS(15,12)+=crLHS298;
rLHS(15,13)+=crLHS300;
rLHS(15,14)+=crLHS301;
rLHS(15,15)+=crLHS298;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_stress = rData.ShearStress;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
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
const double crRHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crRHS5 = 1.0/crRHS4;
const double crRHS6 = crRHS5*p_th;
const double crRHS7 = crRHS6*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]));
const double crRHS8 = -crRHS7 + dp_th_dt;
const double crRHS9 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crRHS10 = crRHS9*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0));
const double crRHS11 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crRHS12 = crRHS11*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1));
const double crRHS13 = crRHS6*(crRHS10 + crRHS12);
const double crRHS14 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crRHS15 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crRHS16 = crRHS0*crRHS14 + crRHS15*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS17 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS18 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS19 = crRHS18*crRHS6;
const double crRHS20 = -crRHS19*(-crRHS16 - crRHS17 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS21 = p_th*tau_u;
const double crRHS22 = crRHS20*crRHS21;
const double crRHS23 = crRHS1*crRHS15 + crRHS14*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS24 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS25 = -crRHS19*(-crRHS23 - crRHS24 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS26 = crRHS21*crRHS25;
const double crRHS27 = gauss_weight*gauss_weight;
const double crRHS28 = crRHS18*crRHS5;
const double crRHS29 = crRHS27*crRHS28;
const double crRHS30 = -crRHS29*(-crRHS13*r_N[0] + crRHS22*r_DN(0,0) + crRHS26*r_DN(0,1) + crRHS3*r_N[0] + crRHS8*r_N[0]);
const double crRHS31 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS32 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS33 = crRHS19*r_N[0];
const double crRHS34 = crRHS28*r_DN(0,0);
const double crRHS35 = tau_c*(-crRHS10*crRHS6 - crRHS12*crRHS6 + crRHS2*p_th - crRHS7 + dp_th_dt);
const double crRHS36 = tau_u*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1));
const double crRHS37 = crRHS33*crRHS36;
const double crRHS38 = crRHS19*tau_u;
const double crRHS39 = crRHS38*(crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1));
const double crRHS40 = crRHS18*(crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS41 = crRHS40*r_N[0];
const double crRHS42 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS43 = crRHS28*r_DN(0,1);
const double crRHS44 = -crRHS29*(-crRHS13*r_N[1] + crRHS22*r_DN(1,0) + crRHS26*r_DN(1,1) + crRHS3*r_N[1] + crRHS8*r_N[1]);
const double crRHS45 = crRHS19*r_N[1];
const double crRHS46 = crRHS28*r_DN(1,0);
const double crRHS47 = crRHS36*crRHS45;
const double crRHS48 = crRHS38*(crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1));
const double crRHS49 = crRHS40*r_N[1];
const double crRHS50 = crRHS28*r_DN(1,1);
const double crRHS51 = -crRHS29*(-crRHS13*r_N[2] + crRHS22*r_DN(2,0) + crRHS26*r_DN(2,1) + crRHS3*r_N[2] + crRHS8*r_N[2]);
const double crRHS52 = crRHS19*r_N[2];
const double crRHS53 = crRHS28*r_DN(2,0);
const double crRHS54 = crRHS36*crRHS52;
const double crRHS55 = crRHS38*(crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1));
const double crRHS56 = crRHS40*r_N[2];
const double crRHS57 = crRHS28*r_DN(2,1);
rRHS[0]+=crRHS30;
rRHS[1]+=-crRHS27*(crRHS16*crRHS33 + crRHS17*crRHS33 + crRHS20*crRHS37 + crRHS20*crRHS39 - crRHS22*crRHS41 + crRHS3*crRHS34 - crRHS31*r_DN(0,0) - crRHS32*crRHS33 + crRHS34*crRHS35 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS27*(crRHS23*crRHS33 + crRHS24*crRHS33 + crRHS25*crRHS37 + crRHS25*crRHS39 - crRHS26*crRHS41 + crRHS3*crRHS43 - crRHS31*r_DN(0,1) - crRHS33*crRHS42 + crRHS35*crRHS43 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS30;
rRHS[4]+=crRHS44;
rRHS[5]+=-crRHS27*(crRHS16*crRHS45 + crRHS17*crRHS45 + crRHS20*crRHS47 + crRHS20*crRHS48 - crRHS22*crRHS49 + crRHS3*crRHS46 - crRHS31*r_DN(1,0) - crRHS32*crRHS45 + crRHS35*crRHS46 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS27*(crRHS23*crRHS45 + crRHS24*crRHS45 + crRHS25*crRHS47 + crRHS25*crRHS48 - crRHS26*crRHS49 + crRHS3*crRHS50 - crRHS31*r_DN(1,1) + crRHS35*crRHS50 - crRHS42*crRHS45 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS44;
rRHS[8]+=crRHS51;
rRHS[9]+=-crRHS27*(crRHS16*crRHS52 + crRHS17*crRHS52 + crRHS20*crRHS54 + crRHS20*crRHS55 - crRHS22*crRHS56 + crRHS3*crRHS53 - crRHS31*r_DN(2,0) - crRHS32*crRHS52 + crRHS35*crRHS53 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS27*(crRHS23*crRHS52 + crRHS24*crRHS52 + crRHS25*crRHS54 + crRHS25*crRHS55 - crRHS26*crRHS56 + crRHS3*crRHS57 - crRHS31*r_DN(2,1) + crRHS35*crRHS57 - crRHS42*crRHS52 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS51;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double gamma = rData.HeatCapacityRatio;

    // Thermodynamic pressure
    const double p_th = rData.ThermodynamicPressure;
    const double dp_th_dt = rData.ThermodynamicPressureDerivative;

    // Material response parameters
    const auto& r_stress = rData.ShearStress;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
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
const double crRHS4 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crRHS5 = 1.0/crRHS4;
const double crRHS6 = crRHS5*p_th;
const double crRHS7 = crRHS6*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]) + r_N[3]*(bdf0*r_t[3] + bdf1*r_t_n[3] + bdf2*r_t_nn[3]));
const double crRHS8 = -crRHS7 + dp_th_dt;
const double crRHS9 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crRHS10 = crRHS9*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0) + r_N[3]*r_u(3,0));
const double crRHS11 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crRHS12 = crRHS11*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1) + r_N[3]*r_u(3,1));
const double crRHS13 = crRHS6*(crRHS10 + crRHS12);
const double crRHS14 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crRHS15 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crRHS16 = crRHS0*crRHS14 + crRHS15*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS17 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0)) + r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS18 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS19 = crRHS18*crRHS6;
const double crRHS20 = -crRHS19*(-crRHS16 - crRHS17 + r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0)) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS21 = p_th*tau_u;
const double crRHS22 = crRHS20*crRHS21;
const double crRHS23 = crRHS1*crRHS15 + crRHS14*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS24 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1)) + r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS25 = -crRHS19*(-crRHS23 - crRHS24 + r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1)) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS26 = crRHS21*crRHS25;
const double crRHS27 = gauss_weight*gauss_weight;
const double crRHS28 = crRHS18*crRHS5;
const double crRHS29 = crRHS27*crRHS28;
const double crRHS30 = -crRHS29*(-crRHS13*r_N[0] + crRHS22*r_DN(0,0) + crRHS26*r_DN(0,1) + crRHS3*r_N[0] + crRHS8*r_N[0]);
const double crRHS31 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS32 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS33 = crRHS19*r_N[0];
const double crRHS34 = crRHS28*r_DN(0,0);
const double crRHS35 = tau_c*(-crRHS10*crRHS6 - crRHS12*crRHS6 + crRHS2*p_th - crRHS7 + dp_th_dt);
const double crRHS36 = tau_u*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1));
const double crRHS37 = crRHS33*crRHS36;
const double crRHS38 = crRHS19*tau_u;
const double crRHS39 = crRHS38*(crRHS14*r_DN(0,0) + crRHS15*r_DN(0,1));
const double crRHS40 = crRHS18*(crRHS11*crRHS15 + crRHS14*crRHS9)*1.0/(crRHS4*crRHS4);
const double crRHS41 = crRHS40*r_N[0];
const double crRHS42 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS43 = crRHS28*r_DN(0,1);
const double crRHS44 = -crRHS29*(-crRHS13*r_N[1] + crRHS22*r_DN(1,0) + crRHS26*r_DN(1,1) + crRHS3*r_N[1] + crRHS8*r_N[1]);
const double crRHS45 = crRHS19*r_N[1];
const double crRHS46 = crRHS28*r_DN(1,0);
const double crRHS47 = crRHS36*crRHS45;
const double crRHS48 = crRHS38*(crRHS14*r_DN(1,0) + crRHS15*r_DN(1,1));
const double crRHS49 = crRHS40*r_N[1];
const double crRHS50 = crRHS28*r_DN(1,1);
const double crRHS51 = -crRHS29*(-crRHS13*r_N[2] + crRHS22*r_DN(2,0) + crRHS26*r_DN(2,1) + crRHS3*r_N[2] + crRHS8*r_N[2]);
const double crRHS52 = crRHS19*r_N[2];
const double crRHS53 = crRHS28*r_DN(2,0);
const double crRHS54 = crRHS36*crRHS52;
const double crRHS55 = crRHS38*(crRHS14*r_DN(2,0) + crRHS15*r_DN(2,1));
const double crRHS56 = crRHS40*r_N[2];
const double crRHS57 = crRHS28*r_DN(2,1);
const double crRHS58 = -crRHS29*(-crRHS13*r_N[3] + crRHS22*r_DN(3,0) + crRHS26*r_DN(3,1) + crRHS3*r_N[3] + crRHS8*r_N[3]);
const double crRHS59 = crRHS19*r_N[3];
const double crRHS60 = crRHS28*r_DN(3,0);
const double crRHS61 = crRHS36*crRHS59;
const double crRHS62 = crRHS38*(crRHS14*r_DN(3,0) + crRHS15*r_DN(3,1));
const double crRHS63 = crRHS40*r_N[3];
const double crRHS64 = crRHS28*r_DN(3,1);
rRHS[0]+=crRHS30;
rRHS[1]+=-crRHS27*(crRHS16*crRHS33 + crRHS17*crRHS33 + crRHS20*crRHS37 + crRHS20*crRHS39 - crRHS22*crRHS41 + crRHS3*crRHS34 - crRHS31*r_DN(0,0) - crRHS32*crRHS33 + crRHS34*crRHS35 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS27*(crRHS23*crRHS33 + crRHS24*crRHS33 + crRHS25*crRHS37 + crRHS25*crRHS39 - crRHS26*crRHS41 + crRHS3*crRHS43 - crRHS31*r_DN(0,1) - crRHS33*crRHS42 + crRHS35*crRHS43 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS30;
rRHS[4]+=crRHS44;
rRHS[5]+=-crRHS27*(crRHS16*crRHS45 + crRHS17*crRHS45 + crRHS20*crRHS47 + crRHS20*crRHS48 - crRHS22*crRHS49 + crRHS3*crRHS46 - crRHS31*r_DN(1,0) - crRHS32*crRHS45 + crRHS35*crRHS46 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS27*(crRHS23*crRHS45 + crRHS24*crRHS45 + crRHS25*crRHS47 + crRHS25*crRHS48 - crRHS26*crRHS49 + crRHS3*crRHS50 - crRHS31*r_DN(1,1) + crRHS35*crRHS50 - crRHS42*crRHS45 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS44;
rRHS[8]+=crRHS51;
rRHS[9]+=-crRHS27*(crRHS16*crRHS52 + crRHS17*crRHS52 + crRHS20*crRHS54 + crRHS20*crRHS55 - crRHS22*crRHS56 + crRHS3*crRHS53 - crRHS31*r_DN(2,0) - crRHS32*crRHS52 + crRHS35*crRHS53 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS27*(crRHS23*crRHS52 + crRHS24*crRHS52 + crRHS25*crRHS54 + crRHS25*crRHS55 - crRHS26*crRHS56 + crRHS3*crRHS57 - crRHS31*r_DN(2,1) + crRHS35*crRHS57 - crRHS42*crRHS52 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS51;
rRHS[12]+=crRHS58;
rRHS[13]+=-crRHS27*(crRHS16*crRHS59 + crRHS17*crRHS59 + crRHS20*crRHS61 + crRHS20*crRHS62 - crRHS22*crRHS63 + crRHS3*crRHS60 - crRHS31*r_DN(3,0) - crRHS32*crRHS59 + crRHS35*crRHS60 + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-crRHS27*(crRHS23*crRHS59 + crRHS24*crRHS59 + crRHS25*crRHS61 + crRHS25*crRHS62 - crRHS26*crRHS63 + crRHS3*crRHS64 - crRHS31*r_DN(3,1) + crRHS35*crRHS64 - crRHS42*crRHS59 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=crRHS58;

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
    const double h = rData.ElementSize;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double dyn_tau = rData.DynamicTau;
    const double mu = rData.EffectiveViscosity;

    double rho_gauss = 0.0;
    array_1d<double,3> u_conv_gauss = ZeroVector(3);
    for (IndexType i_node = 0; i_node < TElementData::NumNodes; ++i_node) {
        rho_gauss += rData.N[i_node] * rData.Density[i_node];
        const auto& r_u_i = row(rData.Velocity, i_node);
        const auto& r_u_m_i = row(rData.MeshVelocity, i_node);
        for (IndexType d = 0; d < TElementData::Dim; ++d) {
            u_conv_gauss[d] += rData.N[i_node] * (r_u_i[d] - r_u_m_i[d]);
        }
    }
    const double norm_u_conv = norm_2(u_conv_gauss);

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