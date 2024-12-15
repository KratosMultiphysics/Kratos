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

template <class TElementData>
void LowMachNavierStokes<TElementData>::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    //TODO: We want to use the constitutive law in the low Mach

    // // Set the constitutive law pointer to null as the axisymmetric element hardcodes a Newtonian fluid viscous behavior
    // // Note that to use a constitutive law the gradient in cylindrical coordinates would require the corresponding stress
    // // implementation in cylindrical coordinates within the mechanical response calculation
    // this->GetConstitutiveLaw() = nullptr;

    KRATOS_CATCH("");
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
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const double tau_c = CalculateTauPressure(rData);
    const double tau_u = CalculateTauVelocity(rData);
    const double tau_t = CalculateTauTemperature(rData);

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = gamma*p_th*1.0/crLHS5*1.0/c_p;
const double crLHS7 = crLHS4*crLHS6;
const double crLHS8 = crLHS2*crLHS7;
const double crLHS9 = crLHS8*tau_u;
const double crLHS10 = crLHS9*(crLHS0 + crLHS1);
const double crLHS11 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crLHS12 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS13 = bdf0*r_N[0];
const double crLHS14 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crLHS15 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crLHS16 = crLHS14*r_DN(0,0) + crLHS15*r_DN(0,1);
const double crLHS17 = crLHS13 + crLHS16;
const double crLHS18 = crLHS7*tau_u;
const double crLHS19 = crLHS17*crLHS18;
const double crLHS20 = crLHS8*(-crLHS11*crLHS12 + crLHS19*r_DN(0,0) + r_DN(0,0)*r_N[0]);
const double crLHS21 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crLHS22 = crLHS8*(-crLHS12*crLHS21 + crLHS19*r_DN(0,1) + r_DN(0,1)*r_N[0]);
const double crLHS23 = r_DN(0,0)*r_DN(1,0);
const double crLHS24 = r_DN(0,1)*r_DN(1,1);
const double crLHS25 = crLHS9*(crLHS23 + crLHS24);
const double crLHS26 = r_DN(1,0)*r_N[0];
const double crLHS27 = crLHS4*r_N[0];
const double crLHS28 = crLHS11*crLHS27;
const double crLHS29 = -crLHS28*r_N[1];
const double crLHS30 = bdf0*r_N[1];
const double crLHS31 = crLHS14*r_DN(1,0) + crLHS15*r_DN(1,1);
const double crLHS32 = crLHS30 + crLHS31;
const double crLHS33 = crLHS18*crLHS32;
const double crLHS34 = crLHS8*(crLHS26 + crLHS29 + crLHS33*r_DN(0,0));
const double crLHS35 = r_DN(1,1)*r_N[0];
const double crLHS36 = crLHS21*crLHS27;
const double crLHS37 = -crLHS36*r_N[1];
const double crLHS38 = crLHS8*(crLHS33*r_DN(0,1) + crLHS35 + crLHS37);
const double crLHS39 = r_DN(0,0)*r_DN(2,0);
const double crLHS40 = r_DN(0,1)*r_DN(2,1);
const double crLHS41 = crLHS9*(crLHS39 + crLHS40);
const double crLHS42 = r_DN(2,0)*r_N[0];
const double crLHS43 = -crLHS28*r_N[2];
const double crLHS44 = crLHS14*r_DN(2,0) + crLHS15*r_DN(2,1);
const double crLHS45 = bdf0*r_N[2] + crLHS44;
const double crLHS46 = crLHS18*crLHS45;
const double crLHS47 = crLHS8*(crLHS42 + crLHS43 + crLHS46*r_DN(0,0));
const double crLHS48 = r_DN(2,1)*r_N[0];
const double crLHS49 = -crLHS36*r_N[2];
const double crLHS50 = crLHS8*(crLHS46*r_DN(0,1) + crLHS48 + crLHS49);
const double crLHS51 = crLHS27*crLHS6;
const double crLHS52 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1);
const double crLHS53 = crLHS52*tau_u;
const double crLHS54 = crLHS16*crLHS18;
const double crLHS55 = crLHS11*crLHS14 + crLHS15*crLHS21;
const double crLHS56 = crLHS55*r_N[0];
const double crLHS57 = 1.0/(crLHS3*crLHS3);
const double crLHS58 = crLHS57*tau_u;
const double crLHS59 = crLHS58*crLHS6;
const double crLHS60 = crLHS2*(crLHS51*crLHS53 + crLHS54 - crLHS56*crLHS59 - r_N[0]);
const double crLHS61 = crLHS60*r_DN(0,0);
const double crLHS62 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS63 = r_C(0,2)*r_DN(0,0);
const double crLHS64 = crLHS63 + r_C(2,2)*r_DN(0,1);
const double crLHS65 = -crLHS28 + r_DN(0,0);
const double crLHS66 = crLHS7*r_DN(0,0);
const double crLHS67 = crLHS66*tau_c;
const double crLHS68 = bdf0*crLHS6;
const double crLHS69 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS70 = crLHS58*crLHS69;
const double crLHS71 = crLHS16*crLHS70;
const double crLHS72 = crLHS17*crLHS69;
const double crLHS73 = crLHS53*crLHS57;
const double crLHS74 = crLHS73*r_N[0];
const double crLHS75 = tau_u*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS76 = crLHS56*crLHS75;
const double crLHS77 = crLHS12*crLHS68 + crLHS16*crLHS51 + crLHS17*crLHS71 + crLHS72*crLHS74 - crLHS72*crLHS76;
const double crLHS78 = crLHS63 + r_C(0,1)*r_DN(0,1);
const double crLHS79 = r_C(1,2)*r_DN(0,1);
const double crLHS80 = crLHS79 + r_C(2,2)*r_DN(0,0);
const double crLHS81 = crLHS66*r_DN(0,1);
const double crLHS82 = -crLHS36 + r_DN(0,1);
const double crLHS83 = r_DN(0,0)*r_N[1];
const double crLHS84 = crLHS18*crLHS52;
const double crLHS85 = crLHS55*crLHS59;
const double crLHS86 = -crLHS2*(-crLHS26*crLHS84 + crLHS26*crLHS85 - crLHS54*r_DN(1,0) + crLHS83);
const double crLHS87 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS88 = r_C(0,2)*r_DN(1,0);
const double crLHS89 = crLHS88 + r_C(2,2)*r_DN(1,1);
const double crLHS90 = crLHS4*r_N[1];
const double crLHS91 = crLHS11*crLHS90;
const double crLHS92 = -crLHS91 + r_DN(1,0);
const double crLHS93 = crLHS6*crLHS90;
const double crLHS94 = crLHS13*crLHS93;
const double crLHS95 = crLHS23*crLHS7 + crLHS94;
const double crLHS96 = crLHS32*crLHS69;
const double crLHS97 = crLHS31*crLHS51 + crLHS32*crLHS71 + crLHS74*crLHS96 - crLHS76*crLHS96;
const double crLHS98 = crLHS88 + r_C(0,1)*r_DN(1,1);
const double crLHS99 = r_C(1,2)*r_DN(1,1);
const double crLHS100 = crLHS99 + r_C(2,2)*r_DN(1,0);
const double crLHS101 = crLHS66*r_DN(1,1);
const double crLHS102 = crLHS21*crLHS90;
const double crLHS103 = -crLHS102 + r_DN(1,1);
const double crLHS104 = r_DN(0,0)*r_N[2];
const double crLHS105 = -crLHS2*(crLHS104 - crLHS42*crLHS84 + crLHS42*crLHS85 - crLHS54*r_DN(2,0));
const double crLHS106 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS107 = r_C(0,2)*r_DN(2,0);
const double crLHS108 = crLHS107 + r_C(2,2)*r_DN(2,1);
const double crLHS109 = crLHS4*r_N[2];
const double crLHS110 = -crLHS109*crLHS11 + r_DN(2,0);
const double crLHS111 = crLHS109*crLHS6;
const double crLHS112 = crLHS111*crLHS13;
const double crLHS113 = crLHS112 + crLHS39*crLHS7;
const double crLHS114 = crLHS45*crLHS69;
const double crLHS115 = crLHS114*crLHS74 - crLHS114*crLHS76 + crLHS44*crLHS51 + crLHS45*crLHS71;
const double crLHS116 = crLHS107 + r_C(0,1)*r_DN(2,1);
const double crLHS117 = r_C(1,2)*r_DN(2,1);
const double crLHS118 = crLHS117 + r_C(2,2)*r_DN(2,0);
const double crLHS119 = crLHS66*r_DN(2,1);
const double crLHS120 = -crLHS109*crLHS21 + r_DN(2,1);
const double crLHS121 = crLHS60*r_DN(0,1);
const double crLHS122 = crLHS79 + r_C(0,1)*r_DN(0,0);
const double crLHS123 = crLHS7*r_DN(0,1);
const double crLHS124 = crLHS123*tau_c;
const double crLHS125 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS126 = r_DN(0,1)*r_N[1];
const double crLHS127 = -crLHS2*(crLHS126 - crLHS35*crLHS84 + crLHS35*crLHS85 - crLHS54*r_DN(1,1));
const double crLHS128 = crLHS99 + r_C(0,1)*r_DN(1,0);
const double crLHS129 = crLHS123*r_DN(1,0);
const double crLHS130 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS131 = crLHS24*crLHS7 + crLHS94;
const double crLHS132 = r_DN(0,1)*r_N[2];
const double crLHS133 = -crLHS2*(crLHS132 - crLHS48*crLHS84 + crLHS48*crLHS85 - crLHS54*r_DN(2,1));
const double crLHS134 = crLHS117 + r_C(0,1)*r_DN(2,0);
const double crLHS135 = crLHS123*r_DN(2,0);
const double crLHS136 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS137 = crLHS112 + crLHS40*crLHS7;
const double crLHS138 = crLHS8*(crLHS19*r_DN(1,0) + crLHS29 + crLHS83);
const double crLHS139 = crLHS8*(crLHS126 + crLHS19*r_DN(1,1) + crLHS37);
const double crLHS140 = r_DN(1,0)*r_DN(1,0);
const double crLHS141 = r_DN(1,1)*r_DN(1,1);
const double crLHS142 = crLHS9*(crLHS140 + crLHS141);
const double crLHS143 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS144 = crLHS8*(-crLHS11*crLHS143 + crLHS33*r_DN(1,0) + r_DN(1,0)*r_N[1]);
const double crLHS145 = crLHS8*(-crLHS143*crLHS21 + crLHS33*r_DN(1,1) + r_DN(1,1)*r_N[1]);
const double crLHS146 = r_DN(1,0)*r_DN(2,0);
const double crLHS147 = r_DN(1,1)*r_DN(2,1);
const double crLHS148 = crLHS9*(crLHS146 + crLHS147);
const double crLHS149 = r_DN(2,0)*r_N[1];
const double crLHS150 = -crLHS91*r_N[2];
const double crLHS151 = crLHS8*(crLHS149 + crLHS150 + crLHS46*r_DN(1,0));
const double crLHS152 = r_DN(2,1)*r_N[1];
const double crLHS153 = -crLHS102*r_N[2];
const double crLHS154 = crLHS8*(crLHS152 + crLHS153 + crLHS46*r_DN(1,1));
const double crLHS155 = crLHS18*crLHS31;
const double crLHS156 = crLHS2*(crLHS155*r_DN(0,0) - crLHS26 + crLHS83*crLHS84 - crLHS83*crLHS85);
const double crLHS157 = crLHS7*r_DN(1,0);
const double crLHS158 = crLHS157*tau_c;
const double crLHS159 = crLHS31*crLHS70;
const double crLHS160 = crLHS72*r_N[1];
const double crLHS161 = crLHS55*crLHS75;
const double crLHS162 = crLHS159*crLHS17 + crLHS16*crLHS93 - crLHS160*crLHS161 + crLHS160*crLHS73;
const double crLHS163 = crLHS2*(crLHS155 + crLHS53*crLHS93 - crLHS85*r_N[1] - r_N[1]);
const double crLHS164 = crLHS163*r_DN(1,0);
const double crLHS165 = crLHS96*r_N[1];
const double crLHS166 = crLHS143*crLHS68 + crLHS159*crLHS32 - crLHS161*crLHS165 + crLHS165*crLHS73 + crLHS31*crLHS93;
const double crLHS167 = crLHS157*r_DN(1,1);
const double crLHS168 = r_DN(1,0)*r_N[2];
const double crLHS169 = -crLHS2*(-crLHS149*crLHS84 + crLHS149*crLHS85 - crLHS155*r_DN(2,0) + crLHS168);
const double crLHS170 = crLHS111*crLHS30;
const double crLHS171 = crLHS146*crLHS7 + crLHS170;
const double crLHS172 = crLHS114*r_N[1];
const double crLHS173 = crLHS159*crLHS45 - crLHS161*crLHS172 + crLHS172*crLHS73 + crLHS44*crLHS93;
const double crLHS174 = crLHS157*r_DN(2,1);
const double crLHS175 = crLHS2*(crLHS126*crLHS84 - crLHS126*crLHS85 + crLHS155*r_DN(0,1) - crLHS35);
const double crLHS176 = crLHS65*tau_c;
const double crLHS177 = crLHS7*r_DN(1,1);
const double crLHS178 = crLHS177*tau_c;
const double crLHS179 = crLHS163*r_DN(1,1);
const double crLHS180 = r_DN(1,1)*r_N[2];
const double crLHS181 = -crLHS2*(-crLHS152*crLHS84 + crLHS152*crLHS85 - crLHS155*r_DN(2,1) + crLHS180);
const double crLHS182 = crLHS7*r_DN(2,0);
const double crLHS183 = crLHS182*r_DN(1,1);
const double crLHS184 = crLHS147*crLHS7 + crLHS170;
const double crLHS185 = crLHS8*(crLHS104 + crLHS19*r_DN(2,0) + crLHS43);
const double crLHS186 = crLHS8*(crLHS132 + crLHS19*r_DN(2,1) + crLHS49);
const double crLHS187 = crLHS8*(crLHS150 + crLHS168 + crLHS33*r_DN(2,0));
const double crLHS188 = crLHS8*(crLHS153 + crLHS180 + crLHS33*r_DN(2,1));
const double crLHS189 = r_DN(2,0)*r_DN(2,0);
const double crLHS190 = r_DN(2,1)*r_DN(2,1);
const double crLHS191 = crLHS9*(crLHS189 + crLHS190);
const double crLHS192 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS193 = crLHS8*(-crLHS11*crLHS192 + crLHS46*r_DN(2,0) + r_DN(2,0)*r_N[2]);
const double crLHS194 = crLHS8*(-crLHS192*crLHS21 + crLHS46*r_DN(2,1) + r_DN(2,1)*r_N[2]);
const double crLHS195 = crLHS18*crLHS44;
const double crLHS196 = crLHS2*(crLHS104*crLHS84 - crLHS104*crLHS85 + crLHS195*r_DN(0,0) - crLHS42);
const double crLHS197 = crLHS44*crLHS70;
const double crLHS198 = crLHS72*r_N[2];
const double crLHS199 = crLHS111*crLHS16 - crLHS161*crLHS198 + crLHS17*crLHS197 + crLHS198*crLHS73;
const double crLHS200 = crLHS182*tau_c;
const double crLHS201 = crLHS2*(-crLHS149 + crLHS168*crLHS84 - crLHS168*crLHS85 + crLHS195*r_DN(1,0));
const double crLHS202 = crLHS96*r_N[2];
const double crLHS203 = crLHS111*crLHS31 - crLHS161*crLHS202 + crLHS197*crLHS32 + crLHS202*crLHS73;
const double crLHS204 = crLHS2*(crLHS111*crLHS53 + crLHS195 - crLHS85*r_N[2] - r_N[2]);
const double crLHS205 = crLHS204*r_DN(2,0);
const double crLHS206 = crLHS114*r_N[2];
const double crLHS207 = crLHS111*crLHS44 - crLHS161*crLHS206 + crLHS192*crLHS68 + crLHS197*crLHS45 + crLHS206*crLHS73;
const double crLHS208 = crLHS182*r_DN(2,1);
const double crLHS209 = crLHS2*(crLHS132*crLHS84 - crLHS132*crLHS85 + crLHS195*r_DN(0,1) - crLHS48);
const double crLHS210 = crLHS7*r_DN(2,1);
const double crLHS211 = crLHS210*tau_c;
const double crLHS212 = crLHS2*(-crLHS152 + crLHS180*crLHS84 - crLHS180*crLHS85 + crLHS195*r_DN(1,1));
const double crLHS213 = crLHS204*r_DN(2,1);
rLHS(0,0)+=crLHS10;
rLHS(0,1)+=crLHS20;
rLHS(0,2)+=crLHS22;
rLHS(0,3)+=crLHS10;
rLHS(0,4)+=crLHS25;
rLHS(0,5)+=crLHS34;
rLHS(0,6)+=crLHS38;
rLHS(0,7)+=crLHS25;
rLHS(0,8)+=crLHS41;
rLHS(0,9)+=crLHS47;
rLHS(0,10)+=crLHS50;
rLHS(0,11)+=crLHS41;
rLHS(1,0)+=crLHS61;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS7 + crLHS62*r_DN(0,0) + crLHS64*r_DN(0,1) + crLHS65*crLHS67 + crLHS77);
rLHS(1,2)+=crLHS2*(crLHS67*crLHS82 + crLHS78*r_DN(0,0) + crLHS80*r_DN(0,1) + crLHS81);
rLHS(1,3)+=crLHS61;
rLHS(1,4)+=crLHS86;
rLHS(1,5)+=crLHS2*(crLHS67*crLHS92 + crLHS87*r_DN(0,0) + crLHS89*r_DN(0,1) + crLHS95 + crLHS97);
rLHS(1,6)+=crLHS2*(crLHS100*r_DN(0,1) + crLHS101 + crLHS103*crLHS67 + crLHS98*r_DN(0,0));
rLHS(1,7)+=crLHS86;
rLHS(1,8)+=crLHS105;
rLHS(1,9)+=crLHS2*(crLHS106*r_DN(0,0) + crLHS108*r_DN(0,1) + crLHS110*crLHS67 + crLHS113 + crLHS115);
rLHS(1,10)+=crLHS2*(crLHS116*r_DN(0,0) + crLHS118*r_DN(0,1) + crLHS119 + crLHS120*crLHS67);
rLHS(1,11)+=crLHS105;
rLHS(2,0)+=crLHS121;
rLHS(2,1)+=crLHS2*(crLHS122*r_DN(0,1) + crLHS124*crLHS65 + crLHS64*r_DN(0,0) + crLHS81);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS7 + crLHS124*crLHS82 + crLHS125*r_DN(0,1) + crLHS77 + crLHS80*r_DN(0,0));
rLHS(2,3)+=crLHS121;
rLHS(2,4)+=crLHS127;
rLHS(2,5)+=crLHS2*(crLHS124*crLHS92 + crLHS128*r_DN(0,1) + crLHS129 + crLHS89*r_DN(0,0));
rLHS(2,6)+=crLHS2*(crLHS100*r_DN(0,0) + crLHS103*crLHS124 + crLHS130*r_DN(0,1) + crLHS131 + crLHS97);
rLHS(2,7)+=crLHS127;
rLHS(2,8)+=crLHS133;
rLHS(2,9)+=crLHS2*(crLHS108*r_DN(0,0) + crLHS110*crLHS124 + crLHS134*r_DN(0,1) + crLHS135);
rLHS(2,10)+=crLHS2*(crLHS115 + crLHS118*r_DN(0,0) + crLHS120*crLHS124 + crLHS136*r_DN(0,1) + crLHS137);
rLHS(2,11)+=crLHS133;
rLHS(3,0)+=crLHS10;
rLHS(3,1)+=crLHS20;
rLHS(3,2)+=crLHS22;
rLHS(3,3)+=crLHS10;
rLHS(3,4)+=crLHS25;
rLHS(3,5)+=crLHS34;
rLHS(3,6)+=crLHS38;
rLHS(3,7)+=crLHS25;
rLHS(3,8)+=crLHS41;
rLHS(3,9)+=crLHS47;
rLHS(3,10)+=crLHS50;
rLHS(3,11)+=crLHS41;
rLHS(4,0)+=crLHS25;
rLHS(4,1)+=crLHS138;
rLHS(4,2)+=crLHS139;
rLHS(4,3)+=crLHS25;
rLHS(4,4)+=crLHS142;
rLHS(4,5)+=crLHS144;
rLHS(4,6)+=crLHS145;
rLHS(4,7)+=crLHS142;
rLHS(4,8)+=crLHS148;
rLHS(4,9)+=crLHS151;
rLHS(4,10)+=crLHS154;
rLHS(4,11)+=crLHS148;
rLHS(5,0)+=crLHS156;
rLHS(5,1)+=crLHS2*(crLHS158*crLHS65 + crLHS162 + crLHS62*r_DN(1,0) + crLHS64*r_DN(1,1) + crLHS95);
rLHS(5,2)+=crLHS2*(crLHS129 + crLHS158*crLHS82 + crLHS78*r_DN(1,0) + crLHS80*r_DN(1,1));
rLHS(5,3)+=crLHS156;
rLHS(5,4)+=crLHS164;
rLHS(5,5)+=crLHS2*(crLHS140*crLHS7 + crLHS158*crLHS92 + crLHS166 + crLHS87*r_DN(1,0) + crLHS89*r_DN(1,1));
rLHS(5,6)+=crLHS2*(crLHS100*r_DN(1,1) + crLHS103*crLHS158 + crLHS167 + crLHS98*r_DN(1,0));
rLHS(5,7)+=crLHS164;
rLHS(5,8)+=crLHS169;
rLHS(5,9)+=crLHS2*(crLHS106*r_DN(1,0) + crLHS108*r_DN(1,1) + crLHS110*crLHS158 + crLHS171 + crLHS173);
rLHS(5,10)+=crLHS2*(crLHS116*r_DN(1,0) + crLHS118*r_DN(1,1) + crLHS120*crLHS158 + crLHS174);
rLHS(5,11)+=crLHS169;
rLHS(6,0)+=crLHS175;
rLHS(6,1)+=crLHS2*(crLHS101 + crLHS122*r_DN(1,1) + crLHS176*crLHS177 + crLHS64*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS125*r_DN(1,1) + crLHS131 + crLHS162 + crLHS178*crLHS82 + crLHS80*r_DN(1,0));
rLHS(6,3)+=crLHS175;
rLHS(6,4)+=crLHS179;
rLHS(6,5)+=crLHS2*(crLHS128*r_DN(1,1) + crLHS167 + crLHS178*crLHS92 + crLHS89*r_DN(1,0));
rLHS(6,6)+=crLHS2*(crLHS100*r_DN(1,0) + crLHS103*crLHS178 + crLHS130*r_DN(1,1) + crLHS141*crLHS7 + crLHS166);
rLHS(6,7)+=crLHS179;
rLHS(6,8)+=crLHS181;
rLHS(6,9)+=crLHS2*(crLHS108*r_DN(1,0) + crLHS110*crLHS178 + crLHS134*r_DN(1,1) + crLHS183);
rLHS(6,10)+=crLHS2*(crLHS118*r_DN(1,0) + crLHS120*crLHS178 + crLHS136*r_DN(1,1) + crLHS173 + crLHS184);
rLHS(6,11)+=crLHS181;
rLHS(7,0)+=crLHS25;
rLHS(7,1)+=crLHS138;
rLHS(7,2)+=crLHS139;
rLHS(7,3)+=crLHS25;
rLHS(7,4)+=crLHS142;
rLHS(7,5)+=crLHS144;
rLHS(7,6)+=crLHS145;
rLHS(7,7)+=crLHS142;
rLHS(7,8)+=crLHS148;
rLHS(7,9)+=crLHS151;
rLHS(7,10)+=crLHS154;
rLHS(7,11)+=crLHS148;
rLHS(8,0)+=crLHS41;
rLHS(8,1)+=crLHS185;
rLHS(8,2)+=crLHS186;
rLHS(8,3)+=crLHS41;
rLHS(8,4)+=crLHS148;
rLHS(8,5)+=crLHS187;
rLHS(8,6)+=crLHS188;
rLHS(8,7)+=crLHS148;
rLHS(8,8)+=crLHS191;
rLHS(8,9)+=crLHS193;
rLHS(8,10)+=crLHS194;
rLHS(8,11)+=crLHS191;
rLHS(9,0)+=crLHS196;
rLHS(9,1)+=crLHS2*(crLHS113 + crLHS176*crLHS182 + crLHS199 + crLHS62*r_DN(2,0) + crLHS64*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS135 + crLHS200*crLHS82 + crLHS78*r_DN(2,0) + crLHS80*r_DN(2,1));
rLHS(9,3)+=crLHS196;
rLHS(9,4)+=crLHS201;
rLHS(9,5)+=crLHS2*(crLHS171 + crLHS200*crLHS92 + crLHS203 + crLHS87*r_DN(2,0) + crLHS89*r_DN(2,1));
rLHS(9,6)+=crLHS2*(crLHS100*r_DN(2,1) + crLHS103*crLHS200 + crLHS183 + crLHS98*r_DN(2,0));
rLHS(9,7)+=crLHS201;
rLHS(9,8)+=crLHS205;
rLHS(9,9)+=crLHS2*(crLHS106*r_DN(2,0) + crLHS108*r_DN(2,1) + crLHS110*crLHS200 + crLHS189*crLHS7 + crLHS207);
rLHS(9,10)+=crLHS2*(crLHS116*r_DN(2,0) + crLHS118*r_DN(2,1) + crLHS120*crLHS200 + crLHS208);
rLHS(9,11)+=crLHS205;
rLHS(10,0)+=crLHS209;
rLHS(10,1)+=crLHS2*(crLHS119 + crLHS122*r_DN(2,1) + crLHS176*crLHS210 + crLHS64*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS125*r_DN(2,1) + crLHS137 + crLHS199 + crLHS211*crLHS82 + crLHS80*r_DN(2,0));
rLHS(10,3)+=crLHS209;
rLHS(10,4)+=crLHS212;
rLHS(10,5)+=crLHS2*(crLHS128*r_DN(2,1) + crLHS174 + crLHS211*crLHS92 + crLHS89*r_DN(2,0));
rLHS(10,6)+=crLHS2*(crLHS100*r_DN(2,0) + crLHS103*crLHS211 + crLHS130*r_DN(2,1) + crLHS184 + crLHS203);
rLHS(10,7)+=crLHS212;
rLHS(10,8)+=crLHS213;
rLHS(10,9)+=crLHS2*(crLHS108*r_DN(2,0) + crLHS110*crLHS211 + crLHS134*r_DN(2,1) + crLHS208);
rLHS(10,10)+=crLHS2*(crLHS118*r_DN(2,0) + crLHS120*crLHS211 + crLHS136*r_DN(2,1) + crLHS190*crLHS7 + crLHS207);
rLHS(10,11)+=crLHS213;
rLHS(11,0)+=crLHS41;
rLHS(11,1)+=crLHS185;
rLHS(11,2)+=crLHS186;
rLHS(11,3)+=crLHS41;
rLHS(11,4)+=crLHS148;
rLHS(11,5)+=crLHS187;
rLHS(11,6)+=crLHS188;
rLHS(11,7)+=crLHS148;
rLHS(11,8)+=crLHS191;
rLHS(11,9)+=crLHS193;
rLHS(11,10)+=crLHS194;
rLHS(11,11)+=crLHS191;

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
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const double tau_c = CalculateTauPressure(rData);
    const double tau_u = CalculateTauVelocity(rData);
    const double tau_t = CalculateTauTemperature(rData);

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = gamma*p_th*1.0/crLHS5*1.0/c_p;
const double crLHS7 = crLHS4*crLHS6;
const double crLHS8 = crLHS2*crLHS7;
const double crLHS9 = crLHS8*tau_u;
const double crLHS10 = crLHS9*(crLHS0 + crLHS1);
const double crLHS11 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crLHS12 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS13 = bdf0*r_N[0];
const double crLHS14 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crLHS15 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crLHS16 = crLHS14*r_DN(0,0) + crLHS15*r_DN(0,1);
const double crLHS17 = crLHS13 + crLHS16;
const double crLHS18 = crLHS7*tau_u;
const double crLHS19 = crLHS17*crLHS18;
const double crLHS20 = crLHS8*(-crLHS11*crLHS12 + crLHS19*r_DN(0,0) + r_DN(0,0)*r_N[0]);
const double crLHS21 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crLHS22 = crLHS8*(-crLHS12*crLHS21 + crLHS19*r_DN(0,1) + r_DN(0,1)*r_N[0]);
const double crLHS23 = r_DN(0,0)*r_DN(1,0);
const double crLHS24 = r_DN(0,1)*r_DN(1,1);
const double crLHS25 = crLHS9*(crLHS23 + crLHS24);
const double crLHS26 = r_DN(1,0)*r_N[0];
const double crLHS27 = crLHS4*r_N[0];
const double crLHS28 = crLHS11*crLHS27;
const double crLHS29 = -crLHS28*r_N[1];
const double crLHS30 = bdf0*r_N[1];
const double crLHS31 = crLHS14*r_DN(1,0) + crLHS15*r_DN(1,1);
const double crLHS32 = crLHS30 + crLHS31;
const double crLHS33 = crLHS18*crLHS32;
const double crLHS34 = crLHS8*(crLHS26 + crLHS29 + crLHS33*r_DN(0,0));
const double crLHS35 = r_DN(1,1)*r_N[0];
const double crLHS36 = crLHS21*crLHS27;
const double crLHS37 = -crLHS36*r_N[1];
const double crLHS38 = crLHS8*(crLHS33*r_DN(0,1) + crLHS35 + crLHS37);
const double crLHS39 = r_DN(0,0)*r_DN(2,0);
const double crLHS40 = r_DN(0,1)*r_DN(2,1);
const double crLHS41 = crLHS9*(crLHS39 + crLHS40);
const double crLHS42 = r_DN(2,0)*r_N[0];
const double crLHS43 = -crLHS28*r_N[2];
const double crLHS44 = bdf0*r_N[2];
const double crLHS45 = crLHS14*r_DN(2,0) + crLHS15*r_DN(2,1);
const double crLHS46 = crLHS44 + crLHS45;
const double crLHS47 = crLHS18*crLHS46;
const double crLHS48 = crLHS8*(crLHS42 + crLHS43 + crLHS47*r_DN(0,0));
const double crLHS49 = r_DN(2,1)*r_N[0];
const double crLHS50 = -crLHS36*r_N[2];
const double crLHS51 = crLHS8*(crLHS47*r_DN(0,1) + crLHS49 + crLHS50);
const double crLHS52 = r_DN(0,0)*r_DN(3,0);
const double crLHS53 = r_DN(0,1)*r_DN(3,1);
const double crLHS54 = crLHS9*(crLHS52 + crLHS53);
const double crLHS55 = r_DN(3,0)*r_N[0];
const double crLHS56 = -crLHS28*r_N[3];
const double crLHS57 = crLHS14*r_DN(3,0) + crLHS15*r_DN(3,1);
const double crLHS58 = bdf0*r_N[3] + crLHS57;
const double crLHS59 = crLHS18*crLHS58;
const double crLHS60 = crLHS8*(crLHS55 + crLHS56 + crLHS59*r_DN(0,0));
const double crLHS61 = r_DN(3,1)*r_N[0];
const double crLHS62 = -crLHS36*r_N[3];
const double crLHS63 = crLHS8*(crLHS59*r_DN(0,1) + crLHS61 + crLHS62);
const double crLHS64 = crLHS27*crLHS6;
const double crLHS65 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1);
const double crLHS66 = crLHS65*tau_u;
const double crLHS67 = crLHS16*crLHS18;
const double crLHS68 = crLHS11*crLHS14 + crLHS15*crLHS21;
const double crLHS69 = crLHS68*r_N[0];
const double crLHS70 = 1.0/(crLHS3*crLHS3);
const double crLHS71 = crLHS70*tau_u;
const double crLHS72 = crLHS6*crLHS71;
const double crLHS73 = crLHS2*(crLHS64*crLHS66 + crLHS67 - crLHS69*crLHS72 - r_N[0]);
const double crLHS74 = crLHS73*r_DN(0,0);
const double crLHS75 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS76 = r_C(0,2)*r_DN(0,0);
const double crLHS77 = crLHS76 + r_C(2,2)*r_DN(0,1);
const double crLHS78 = -crLHS28 + r_DN(0,0);
const double crLHS79 = crLHS7*r_DN(0,0);
const double crLHS80 = crLHS79*tau_c;
const double crLHS81 = bdf0*crLHS6;
const double crLHS82 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS83 = crLHS71*crLHS82;
const double crLHS84 = crLHS16*crLHS83;
const double crLHS85 = crLHS17*crLHS82;
const double crLHS86 = crLHS66*crLHS70;
const double crLHS87 = crLHS86*r_N[0];
const double crLHS88 = tau_u*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS89 = crLHS69*crLHS88;
const double crLHS90 = crLHS12*crLHS81 + crLHS16*crLHS64 + crLHS17*crLHS84 + crLHS85*crLHS87 - crLHS85*crLHS89;
const double crLHS91 = crLHS76 + r_C(0,1)*r_DN(0,1);
const double crLHS92 = r_C(1,2)*r_DN(0,1);
const double crLHS93 = crLHS92 + r_C(2,2)*r_DN(0,0);
const double crLHS94 = crLHS79*r_DN(0,1);
const double crLHS95 = -crLHS36 + r_DN(0,1);
const double crLHS96 = r_DN(0,0)*r_N[1];
const double crLHS97 = crLHS18*crLHS65;
const double crLHS98 = crLHS68*crLHS72;
const double crLHS99 = -crLHS2*(-crLHS26*crLHS97 + crLHS26*crLHS98 - crLHS67*r_DN(1,0) + crLHS96);
const double crLHS100 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS101 = r_C(0,2)*r_DN(1,0);
const double crLHS102 = crLHS101 + r_C(2,2)*r_DN(1,1);
const double crLHS103 = crLHS4*r_N[1];
const double crLHS104 = crLHS103*crLHS11;
const double crLHS105 = -crLHS104 + r_DN(1,0);
const double crLHS106 = crLHS103*crLHS6;
const double crLHS107 = crLHS106*crLHS13;
const double crLHS108 = crLHS107 + crLHS23*crLHS7;
const double crLHS109 = crLHS32*crLHS82;
const double crLHS110 = crLHS109*crLHS87 - crLHS109*crLHS89 + crLHS31*crLHS64 + crLHS32*crLHS84;
const double crLHS111 = crLHS101 + r_C(0,1)*r_DN(1,1);
const double crLHS112 = r_C(1,2)*r_DN(1,1);
const double crLHS113 = crLHS112 + r_C(2,2)*r_DN(1,0);
const double crLHS114 = crLHS79*r_DN(1,1);
const double crLHS115 = crLHS103*crLHS21;
const double crLHS116 = -crLHS115 + r_DN(1,1);
const double crLHS117 = r_DN(0,0)*r_N[2];
const double crLHS118 = -crLHS2*(crLHS117 - crLHS42*crLHS97 + crLHS42*crLHS98 - crLHS67*r_DN(2,0));
const double crLHS119 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS120 = r_C(0,2)*r_DN(2,0);
const double crLHS121 = crLHS120 + r_C(2,2)*r_DN(2,1);
const double crLHS122 = crLHS4*r_N[2];
const double crLHS123 = crLHS11*crLHS122;
const double crLHS124 = -crLHS123 + r_DN(2,0);
const double crLHS125 = crLHS122*crLHS6;
const double crLHS126 = crLHS125*crLHS13;
const double crLHS127 = crLHS126 + crLHS39*crLHS7;
const double crLHS128 = crLHS46*crLHS82;
const double crLHS129 = crLHS128*crLHS87 - crLHS128*crLHS89 + crLHS45*crLHS64 + crLHS46*crLHS84;
const double crLHS130 = crLHS120 + r_C(0,1)*r_DN(2,1);
const double crLHS131 = r_C(1,2)*r_DN(2,1);
const double crLHS132 = crLHS131 + r_C(2,2)*r_DN(2,0);
const double crLHS133 = crLHS79*r_DN(2,1);
const double crLHS134 = crLHS122*crLHS21;
const double crLHS135 = -crLHS134 + r_DN(2,1);
const double crLHS136 = r_DN(0,0)*r_N[3];
const double crLHS137 = -crLHS2*(crLHS136 - crLHS55*crLHS97 + crLHS55*crLHS98 - crLHS67*r_DN(3,0));
const double crLHS138 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS139 = r_C(0,2)*r_DN(3,0);
const double crLHS140 = crLHS139 + r_C(2,2)*r_DN(3,1);
const double crLHS141 = crLHS4*r_N[3];
const double crLHS142 = -crLHS11*crLHS141 + r_DN(3,0);
const double crLHS143 = crLHS141*crLHS6;
const double crLHS144 = crLHS13*crLHS143;
const double crLHS145 = crLHS144 + crLHS52*crLHS7;
const double crLHS146 = crLHS58*crLHS82;
const double crLHS147 = crLHS146*crLHS87 - crLHS146*crLHS89 + crLHS57*crLHS64 + crLHS58*crLHS84;
const double crLHS148 = crLHS139 + r_C(0,1)*r_DN(3,1);
const double crLHS149 = r_C(1,2)*r_DN(3,1);
const double crLHS150 = crLHS149 + r_C(2,2)*r_DN(3,0);
const double crLHS151 = crLHS79*r_DN(3,1);
const double crLHS152 = -crLHS141*crLHS21 + r_DN(3,1);
const double crLHS153 = crLHS73*r_DN(0,1);
const double crLHS154 = crLHS92 + r_C(0,1)*r_DN(0,0);
const double crLHS155 = crLHS7*r_DN(0,1);
const double crLHS156 = crLHS155*tau_c;
const double crLHS157 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS158 = r_DN(0,1)*r_N[1];
const double crLHS159 = -crLHS2*(crLHS158 - crLHS35*crLHS97 + crLHS35*crLHS98 - crLHS67*r_DN(1,1));
const double crLHS160 = crLHS112 + r_C(0,1)*r_DN(1,0);
const double crLHS161 = crLHS155*r_DN(1,0);
const double crLHS162 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS163 = crLHS107 + crLHS24*crLHS7;
const double crLHS164 = r_DN(0,1)*r_N[2];
const double crLHS165 = -crLHS2*(crLHS164 - crLHS49*crLHS97 + crLHS49*crLHS98 - crLHS67*r_DN(2,1));
const double crLHS166 = crLHS131 + r_C(0,1)*r_DN(2,0);
const double crLHS167 = crLHS155*r_DN(2,0);
const double crLHS168 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS169 = crLHS126 + crLHS40*crLHS7;
const double crLHS170 = r_DN(0,1)*r_N[3];
const double crLHS171 = -crLHS2*(crLHS170 - crLHS61*crLHS97 + crLHS61*crLHS98 - crLHS67*r_DN(3,1));
const double crLHS172 = crLHS149 + r_C(0,1)*r_DN(3,0);
const double crLHS173 = crLHS155*r_DN(3,0);
const double crLHS174 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS175 = crLHS144 + crLHS53*crLHS7;
const double crLHS176 = crLHS8*(crLHS19*r_DN(1,0) + crLHS29 + crLHS96);
const double crLHS177 = crLHS8*(crLHS158 + crLHS19*r_DN(1,1) + crLHS37);
const double crLHS178 = r_DN(1,0)*r_DN(1,0);
const double crLHS179 = r_DN(1,1)*r_DN(1,1);
const double crLHS180 = crLHS9*(crLHS178 + crLHS179);
const double crLHS181 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS182 = crLHS8*(-crLHS11*crLHS181 + crLHS33*r_DN(1,0) + r_DN(1,0)*r_N[1]);
const double crLHS183 = crLHS8*(-crLHS181*crLHS21 + crLHS33*r_DN(1,1) + r_DN(1,1)*r_N[1]);
const double crLHS184 = r_DN(1,0)*r_DN(2,0);
const double crLHS185 = r_DN(1,1)*r_DN(2,1);
const double crLHS186 = crLHS9*(crLHS184 + crLHS185);
const double crLHS187 = r_DN(2,0)*r_N[1];
const double crLHS188 = -crLHS104*r_N[2];
const double crLHS189 = crLHS8*(crLHS187 + crLHS188 + crLHS47*r_DN(1,0));
const double crLHS190 = r_DN(2,1)*r_N[1];
const double crLHS191 = -crLHS115*r_N[2];
const double crLHS192 = crLHS8*(crLHS190 + crLHS191 + crLHS47*r_DN(1,1));
const double crLHS193 = r_DN(1,0)*r_DN(3,0);
const double crLHS194 = r_DN(1,1)*r_DN(3,1);
const double crLHS195 = crLHS9*(crLHS193 + crLHS194);
const double crLHS196 = r_DN(3,0)*r_N[1];
const double crLHS197 = -crLHS104*r_N[3];
const double crLHS198 = crLHS8*(crLHS196 + crLHS197 + crLHS59*r_DN(1,0));
const double crLHS199 = r_DN(3,1)*r_N[1];
const double crLHS200 = -crLHS115*r_N[3];
const double crLHS201 = crLHS8*(crLHS199 + crLHS200 + crLHS59*r_DN(1,1));
const double crLHS202 = crLHS18*crLHS31;
const double crLHS203 = crLHS2*(crLHS202*r_DN(0,0) - crLHS26 + crLHS96*crLHS97 - crLHS96*crLHS98);
const double crLHS204 = crLHS7*r_DN(1,0);
const double crLHS205 = crLHS204*tau_c;
const double crLHS206 = crLHS31*crLHS83;
const double crLHS207 = crLHS85*r_N[1];
const double crLHS208 = crLHS68*crLHS88;
const double crLHS209 = crLHS106*crLHS16 + crLHS17*crLHS206 - crLHS207*crLHS208 + crLHS207*crLHS86;
const double crLHS210 = crLHS2*(crLHS106*crLHS66 + crLHS202 - crLHS98*r_N[1] - r_N[1]);
const double crLHS211 = crLHS210*r_DN(1,0);
const double crLHS212 = crLHS109*r_N[1];
const double crLHS213 = crLHS106*crLHS31 + crLHS181*crLHS81 + crLHS206*crLHS32 - crLHS208*crLHS212 + crLHS212*crLHS86;
const double crLHS214 = crLHS204*r_DN(1,1);
const double crLHS215 = r_DN(1,0)*r_N[2];
const double crLHS216 = -crLHS2*(-crLHS187*crLHS97 + crLHS187*crLHS98 - crLHS202*r_DN(2,0) + crLHS215);
const double crLHS217 = crLHS125*crLHS30;
const double crLHS218 = crLHS184*crLHS7 + crLHS217;
const double crLHS219 = crLHS128*r_N[1];
const double crLHS220 = crLHS106*crLHS45 + crLHS206*crLHS46 - crLHS208*crLHS219 + crLHS219*crLHS86;
const double crLHS221 = crLHS204*r_DN(2,1);
const double crLHS222 = r_DN(1,0)*r_N[3];
const double crLHS223 = -crLHS2*(-crLHS196*crLHS97 + crLHS196*crLHS98 - crLHS202*r_DN(3,0) + crLHS222);
const double crLHS224 = crLHS143*crLHS30;
const double crLHS225 = crLHS193*crLHS7 + crLHS224;
const double crLHS226 = crLHS146*r_N[1];
const double crLHS227 = crLHS106*crLHS57 + crLHS206*crLHS58 - crLHS208*crLHS226 + crLHS226*crLHS86;
const double crLHS228 = crLHS204*r_DN(3,1);
const double crLHS229 = crLHS2*(crLHS158*crLHS97 - crLHS158*crLHS98 + crLHS202*r_DN(0,1) - crLHS35);
const double crLHS230 = crLHS7*r_DN(1,1);
const double crLHS231 = crLHS230*tau_c;
const double crLHS232 = crLHS210*r_DN(1,1);
const double crLHS233 = r_DN(1,1)*r_N[2];
const double crLHS234 = -crLHS2*(-crLHS190*crLHS97 + crLHS190*crLHS98 - crLHS202*r_DN(2,1) + crLHS233);
const double crLHS235 = crLHS230*r_DN(2,0);
const double crLHS236 = crLHS185*crLHS7 + crLHS217;
const double crLHS237 = r_DN(1,1)*r_N[3];
const double crLHS238 = -crLHS2*(-crLHS199*crLHS97 + crLHS199*crLHS98 - crLHS202*r_DN(3,1) + crLHS237);
const double crLHS239 = crLHS230*r_DN(3,0);
const double crLHS240 = crLHS194*crLHS7 + crLHS224;
const double crLHS241 = crLHS8*(crLHS117 + crLHS19*r_DN(2,0) + crLHS43);
const double crLHS242 = crLHS8*(crLHS164 + crLHS19*r_DN(2,1) + crLHS50);
const double crLHS243 = crLHS8*(crLHS188 + crLHS215 + crLHS33*r_DN(2,0));
const double crLHS244 = crLHS8*(crLHS191 + crLHS233 + crLHS33*r_DN(2,1));
const double crLHS245 = r_DN(2,0)*r_DN(2,0);
const double crLHS246 = r_DN(2,1)*r_DN(2,1);
const double crLHS247 = crLHS9*(crLHS245 + crLHS246);
const double crLHS248 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS249 = crLHS8*(-crLHS11*crLHS248 + crLHS47*r_DN(2,0) + r_DN(2,0)*r_N[2]);
const double crLHS250 = crLHS8*(-crLHS21*crLHS248 + crLHS47*r_DN(2,1) + r_DN(2,1)*r_N[2]);
const double crLHS251 = r_DN(2,0)*r_DN(3,0);
const double crLHS252 = r_DN(2,1)*r_DN(3,1);
const double crLHS253 = crLHS9*(crLHS251 + crLHS252);
const double crLHS254 = r_DN(3,0)*r_N[2];
const double crLHS255 = -crLHS123*r_N[3];
const double crLHS256 = crLHS8*(crLHS254 + crLHS255 + crLHS59*r_DN(2,0));
const double crLHS257 = r_DN(3,1)*r_N[2];
const double crLHS258 = -crLHS134*r_N[3];
const double crLHS259 = crLHS8*(crLHS257 + crLHS258 + crLHS59*r_DN(2,1));
const double crLHS260 = crLHS18*crLHS45;
const double crLHS261 = crLHS2*(crLHS117*crLHS97 - crLHS117*crLHS98 + crLHS260*r_DN(0,0) - crLHS42);
const double crLHS262 = crLHS7*r_DN(2,0);
const double crLHS263 = crLHS262*tau_c;
const double crLHS264 = crLHS45*crLHS83;
const double crLHS265 = crLHS85*r_N[2];
const double crLHS266 = crLHS125*crLHS16 + crLHS17*crLHS264 - crLHS208*crLHS265 + crLHS265*crLHS86;
const double crLHS267 = crLHS2*(-crLHS187 + crLHS215*crLHS97 - crLHS215*crLHS98 + crLHS260*r_DN(1,0));
const double crLHS268 = crLHS109*r_N[2];
const double crLHS269 = crLHS125*crLHS31 - crLHS208*crLHS268 + crLHS264*crLHS32 + crLHS268*crLHS86;
const double crLHS270 = crLHS2*(crLHS125*crLHS66 + crLHS260 - crLHS98*r_N[2] - r_N[2]);
const double crLHS271 = crLHS270*r_DN(2,0);
const double crLHS272 = crLHS128*r_N[2];
const double crLHS273 = crLHS125*crLHS45 - crLHS208*crLHS272 + crLHS248*crLHS81 + crLHS264*crLHS46 + crLHS272*crLHS86;
const double crLHS274 = crLHS262*r_DN(2,1);
const double crLHS275 = r_DN(2,0)*r_N[3];
const double crLHS276 = -crLHS2*(-crLHS254*crLHS97 + crLHS254*crLHS98 - crLHS260*r_DN(3,0) + crLHS275);
const double crLHS277 = crLHS143*crLHS44;
const double crLHS278 = crLHS251*crLHS7 + crLHS277;
const double crLHS279 = crLHS146*r_N[2];
const double crLHS280 = crLHS125*crLHS57 - crLHS208*crLHS279 + crLHS264*crLHS58 + crLHS279*crLHS86;
const double crLHS281 = crLHS262*r_DN(3,1);
const double crLHS282 = crLHS2*(crLHS164*crLHS97 - crLHS164*crLHS98 + crLHS260*r_DN(0,1) - crLHS49);
const double crLHS283 = crLHS78*tau_c;
const double crLHS284 = crLHS7*r_DN(2,1);
const double crLHS285 = crLHS284*tau_c;
const double crLHS286 = crLHS2*(-crLHS190 + crLHS233*crLHS97 - crLHS233*crLHS98 + crLHS260*r_DN(1,1));
const double crLHS287 = crLHS270*r_DN(2,1);
const double crLHS288 = r_DN(2,1)*r_N[3];
const double crLHS289 = -crLHS2*(-crLHS257*crLHS97 + crLHS257*crLHS98 - crLHS260*r_DN(3,1) + crLHS288);
const double crLHS290 = crLHS7*r_DN(3,0);
const double crLHS291 = crLHS290*r_DN(2,1);
const double crLHS292 = crLHS252*crLHS7 + crLHS277;
const double crLHS293 = crLHS8*(crLHS136 + crLHS19*r_DN(3,0) + crLHS56);
const double crLHS294 = crLHS8*(crLHS170 + crLHS19*r_DN(3,1) + crLHS62);
const double crLHS295 = crLHS8*(crLHS197 + crLHS222 + crLHS33*r_DN(3,0));
const double crLHS296 = crLHS8*(crLHS200 + crLHS237 + crLHS33*r_DN(3,1));
const double crLHS297 = crLHS8*(crLHS255 + crLHS275 + crLHS47*r_DN(3,0));
const double crLHS298 = crLHS8*(crLHS258 + crLHS288 + crLHS47*r_DN(3,1));
const double crLHS299 = r_DN(3,0)*r_DN(3,0);
const double crLHS300 = r_DN(3,1)*r_DN(3,1);
const double crLHS301 = crLHS9*(crLHS299 + crLHS300);
const double crLHS302 = crLHS4*(r_N[3]*r_N[3]);
const double crLHS303 = crLHS8*(-crLHS11*crLHS302 + crLHS59*r_DN(3,0) + r_DN(3,0)*r_N[3]);
const double crLHS304 = crLHS8*(-crLHS21*crLHS302 + crLHS59*r_DN(3,1) + r_DN(3,1)*r_N[3]);
const double crLHS305 = crLHS18*crLHS57;
const double crLHS306 = crLHS2*(crLHS136*crLHS97 - crLHS136*crLHS98 + crLHS305*r_DN(0,0) - crLHS55);
const double crLHS307 = crLHS57*crLHS83;
const double crLHS308 = crLHS85*r_N[3];
const double crLHS309 = crLHS143*crLHS16 + crLHS17*crLHS307 - crLHS208*crLHS308 + crLHS308*crLHS86;
const double crLHS310 = crLHS290*tau_c;
const double crLHS311 = crLHS2*(-crLHS196 + crLHS222*crLHS97 - crLHS222*crLHS98 + crLHS305*r_DN(1,0));
const double crLHS312 = crLHS109*r_N[3];
const double crLHS313 = crLHS143*crLHS31 - crLHS208*crLHS312 + crLHS307*crLHS32 + crLHS312*crLHS86;
const double crLHS314 = crLHS2*(-crLHS254 + crLHS275*crLHS97 - crLHS275*crLHS98 + crLHS305*r_DN(2,0));
const double crLHS315 = crLHS128*r_N[3];
const double crLHS316 = crLHS143*crLHS45 - crLHS208*crLHS315 + crLHS307*crLHS46 + crLHS315*crLHS86;
const double crLHS317 = crLHS2*(crLHS143*crLHS66 + crLHS305 - crLHS98*r_N[3] - r_N[3]);
const double crLHS318 = crLHS317*r_DN(3,0);
const double crLHS319 = crLHS146*r_N[3];
const double crLHS320 = crLHS143*crLHS57 - crLHS208*crLHS319 + crLHS302*crLHS81 + crLHS307*crLHS58 + crLHS319*crLHS86;
const double crLHS321 = crLHS290*r_DN(3,1);
const double crLHS322 = crLHS2*(crLHS170*crLHS97 - crLHS170*crLHS98 + crLHS305*r_DN(0,1) - crLHS61);
const double crLHS323 = crLHS7*r_DN(3,1);
const double crLHS324 = crLHS323*tau_c;
const double crLHS325 = crLHS2*(-crLHS199 + crLHS237*crLHS97 - crLHS237*crLHS98 + crLHS305*r_DN(1,1));
const double crLHS326 = crLHS2*(-crLHS257 + crLHS288*crLHS97 - crLHS288*crLHS98 + crLHS305*r_DN(2,1));
const double crLHS327 = crLHS317*r_DN(3,1);
rLHS(0,0)+=crLHS10;
rLHS(0,1)+=crLHS20;
rLHS(0,2)+=crLHS22;
rLHS(0,3)+=crLHS10;
rLHS(0,4)+=crLHS25;
rLHS(0,5)+=crLHS34;
rLHS(0,6)+=crLHS38;
rLHS(0,7)+=crLHS25;
rLHS(0,8)+=crLHS41;
rLHS(0,9)+=crLHS48;
rLHS(0,10)+=crLHS51;
rLHS(0,11)+=crLHS41;
rLHS(0,12)+=crLHS54;
rLHS(0,13)+=crLHS60;
rLHS(0,14)+=crLHS63;
rLHS(0,15)+=crLHS54;
rLHS(1,0)+=crLHS74;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS7 + crLHS75*r_DN(0,0) + crLHS77*r_DN(0,1) + crLHS78*crLHS80 + crLHS90);
rLHS(1,2)+=crLHS2*(crLHS80*crLHS95 + crLHS91*r_DN(0,0) + crLHS93*r_DN(0,1) + crLHS94);
rLHS(1,3)+=crLHS74;
rLHS(1,4)+=crLHS99;
rLHS(1,5)+=crLHS2*(crLHS100*r_DN(0,0) + crLHS102*r_DN(0,1) + crLHS105*crLHS80 + crLHS108 + crLHS110);
rLHS(1,6)+=crLHS2*(crLHS111*r_DN(0,0) + crLHS113*r_DN(0,1) + crLHS114 + crLHS116*crLHS80);
rLHS(1,7)+=crLHS99;
rLHS(1,8)+=crLHS118;
rLHS(1,9)+=crLHS2*(crLHS119*r_DN(0,0) + crLHS121*r_DN(0,1) + crLHS124*crLHS80 + crLHS127 + crLHS129);
rLHS(1,10)+=crLHS2*(crLHS130*r_DN(0,0) + crLHS132*r_DN(0,1) + crLHS133 + crLHS135*crLHS80);
rLHS(1,11)+=crLHS118;
rLHS(1,12)+=crLHS137;
rLHS(1,13)+=crLHS2*(crLHS138*r_DN(0,0) + crLHS140*r_DN(0,1) + crLHS142*crLHS80 + crLHS145 + crLHS147);
rLHS(1,14)+=crLHS2*(crLHS148*r_DN(0,0) + crLHS150*r_DN(0,1) + crLHS151 + crLHS152*crLHS80);
rLHS(1,15)+=crLHS137;
rLHS(2,0)+=crLHS153;
rLHS(2,1)+=crLHS2*(crLHS154*r_DN(0,1) + crLHS156*crLHS78 + crLHS77*r_DN(0,0) + crLHS94);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS7 + crLHS156*crLHS95 + crLHS157*r_DN(0,1) + crLHS90 + crLHS93*r_DN(0,0));
rLHS(2,3)+=crLHS153;
rLHS(2,4)+=crLHS159;
rLHS(2,5)+=crLHS2*(crLHS102*r_DN(0,0) + crLHS105*crLHS156 + crLHS160*r_DN(0,1) + crLHS161);
rLHS(2,6)+=crLHS2*(crLHS110 + crLHS113*r_DN(0,0) + crLHS116*crLHS156 + crLHS162*r_DN(0,1) + crLHS163);
rLHS(2,7)+=crLHS159;
rLHS(2,8)+=crLHS165;
rLHS(2,9)+=crLHS2*(crLHS121*r_DN(0,0) + crLHS124*crLHS156 + crLHS166*r_DN(0,1) + crLHS167);
rLHS(2,10)+=crLHS2*(crLHS129 + crLHS132*r_DN(0,0) + crLHS135*crLHS156 + crLHS168*r_DN(0,1) + crLHS169);
rLHS(2,11)+=crLHS165;
rLHS(2,12)+=crLHS171;
rLHS(2,13)+=crLHS2*(crLHS140*r_DN(0,0) + crLHS142*crLHS156 + crLHS172*r_DN(0,1) + crLHS173);
rLHS(2,14)+=crLHS2*(crLHS147 + crLHS150*r_DN(0,0) + crLHS152*crLHS156 + crLHS174*r_DN(0,1) + crLHS175);
rLHS(2,15)+=crLHS171;
rLHS(3,0)+=crLHS10;
rLHS(3,1)+=crLHS20;
rLHS(3,2)+=crLHS22;
rLHS(3,3)+=crLHS10;
rLHS(3,4)+=crLHS25;
rLHS(3,5)+=crLHS34;
rLHS(3,6)+=crLHS38;
rLHS(3,7)+=crLHS25;
rLHS(3,8)+=crLHS41;
rLHS(3,9)+=crLHS48;
rLHS(3,10)+=crLHS51;
rLHS(3,11)+=crLHS41;
rLHS(3,12)+=crLHS54;
rLHS(3,13)+=crLHS60;
rLHS(3,14)+=crLHS63;
rLHS(3,15)+=crLHS54;
rLHS(4,0)+=crLHS25;
rLHS(4,1)+=crLHS176;
rLHS(4,2)+=crLHS177;
rLHS(4,3)+=crLHS25;
rLHS(4,4)+=crLHS180;
rLHS(4,5)+=crLHS182;
rLHS(4,6)+=crLHS183;
rLHS(4,7)+=crLHS180;
rLHS(4,8)+=crLHS186;
rLHS(4,9)+=crLHS189;
rLHS(4,10)+=crLHS192;
rLHS(4,11)+=crLHS186;
rLHS(4,12)+=crLHS195;
rLHS(4,13)+=crLHS198;
rLHS(4,14)+=crLHS201;
rLHS(4,15)+=crLHS195;
rLHS(5,0)+=crLHS203;
rLHS(5,1)+=crLHS2*(crLHS108 + crLHS205*crLHS78 + crLHS209 + crLHS75*r_DN(1,0) + crLHS77*r_DN(1,1));
rLHS(5,2)+=crLHS2*(crLHS161 + crLHS205*crLHS95 + crLHS91*r_DN(1,0) + crLHS93*r_DN(1,1));
rLHS(5,3)+=crLHS203;
rLHS(5,4)+=crLHS211;
rLHS(5,5)+=crLHS2*(crLHS100*r_DN(1,0) + crLHS102*r_DN(1,1) + crLHS105*crLHS205 + crLHS178*crLHS7 + crLHS213);
rLHS(5,6)+=crLHS2*(crLHS111*r_DN(1,0) + crLHS113*r_DN(1,1) + crLHS116*crLHS205 + crLHS214);
rLHS(5,7)+=crLHS211;
rLHS(5,8)+=crLHS216;
rLHS(5,9)+=crLHS2*(crLHS119*r_DN(1,0) + crLHS121*r_DN(1,1) + crLHS124*crLHS205 + crLHS218 + crLHS220);
rLHS(5,10)+=crLHS2*(crLHS130*r_DN(1,0) + crLHS132*r_DN(1,1) + crLHS135*crLHS205 + crLHS221);
rLHS(5,11)+=crLHS216;
rLHS(5,12)+=crLHS223;
rLHS(5,13)+=crLHS2*(crLHS138*r_DN(1,0) + crLHS140*r_DN(1,1) + crLHS142*crLHS205 + crLHS225 + crLHS227);
rLHS(5,14)+=crLHS2*(crLHS148*r_DN(1,0) + crLHS150*r_DN(1,1) + crLHS152*crLHS205 + crLHS228);
rLHS(5,15)+=crLHS223;
rLHS(6,0)+=crLHS229;
rLHS(6,1)+=crLHS2*(crLHS114 + crLHS154*r_DN(1,1) + crLHS231*crLHS78 + crLHS77*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS157*r_DN(1,1) + crLHS163 + crLHS209 + crLHS231*crLHS95 + crLHS93*r_DN(1,0));
rLHS(6,3)+=crLHS229;
rLHS(6,4)+=crLHS232;
rLHS(6,5)+=crLHS2*(crLHS102*r_DN(1,0) + crLHS105*crLHS231 + crLHS160*r_DN(1,1) + crLHS214);
rLHS(6,6)+=crLHS2*(crLHS113*r_DN(1,0) + crLHS116*crLHS231 + crLHS162*r_DN(1,1) + crLHS179*crLHS7 + crLHS213);
rLHS(6,7)+=crLHS232;
rLHS(6,8)+=crLHS234;
rLHS(6,9)+=crLHS2*(crLHS121*r_DN(1,0) + crLHS124*crLHS231 + crLHS166*r_DN(1,1) + crLHS235);
rLHS(6,10)+=crLHS2*(crLHS132*r_DN(1,0) + crLHS135*crLHS231 + crLHS168*r_DN(1,1) + crLHS220 + crLHS236);
rLHS(6,11)+=crLHS234;
rLHS(6,12)+=crLHS238;
rLHS(6,13)+=crLHS2*(crLHS140*r_DN(1,0) + crLHS142*crLHS231 + crLHS172*r_DN(1,1) + crLHS239);
rLHS(6,14)+=crLHS2*(crLHS150*r_DN(1,0) + crLHS152*crLHS231 + crLHS174*r_DN(1,1) + crLHS227 + crLHS240);
rLHS(6,15)+=crLHS238;
rLHS(7,0)+=crLHS25;
rLHS(7,1)+=crLHS176;
rLHS(7,2)+=crLHS177;
rLHS(7,3)+=crLHS25;
rLHS(7,4)+=crLHS180;
rLHS(7,5)+=crLHS182;
rLHS(7,6)+=crLHS183;
rLHS(7,7)+=crLHS180;
rLHS(7,8)+=crLHS186;
rLHS(7,9)+=crLHS189;
rLHS(7,10)+=crLHS192;
rLHS(7,11)+=crLHS186;
rLHS(7,12)+=crLHS195;
rLHS(7,13)+=crLHS198;
rLHS(7,14)+=crLHS201;
rLHS(7,15)+=crLHS195;
rLHS(8,0)+=crLHS41;
rLHS(8,1)+=crLHS241;
rLHS(8,2)+=crLHS242;
rLHS(8,3)+=crLHS41;
rLHS(8,4)+=crLHS186;
rLHS(8,5)+=crLHS243;
rLHS(8,6)+=crLHS244;
rLHS(8,7)+=crLHS186;
rLHS(8,8)+=crLHS247;
rLHS(8,9)+=crLHS249;
rLHS(8,10)+=crLHS250;
rLHS(8,11)+=crLHS247;
rLHS(8,12)+=crLHS253;
rLHS(8,13)+=crLHS256;
rLHS(8,14)+=crLHS259;
rLHS(8,15)+=crLHS253;
rLHS(9,0)+=crLHS261;
rLHS(9,1)+=crLHS2*(crLHS127 + crLHS263*crLHS78 + crLHS266 + crLHS75*r_DN(2,0) + crLHS77*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS167 + crLHS263*crLHS95 + crLHS91*r_DN(2,0) + crLHS93*r_DN(2,1));
rLHS(9,3)+=crLHS261;
rLHS(9,4)+=crLHS267;
rLHS(9,5)+=crLHS2*(crLHS100*r_DN(2,0) + crLHS102*r_DN(2,1) + crLHS105*crLHS263 + crLHS218 + crLHS269);
rLHS(9,6)+=crLHS2*(crLHS111*r_DN(2,0) + crLHS113*r_DN(2,1) + crLHS116*crLHS263 + crLHS235);
rLHS(9,7)+=crLHS267;
rLHS(9,8)+=crLHS271;
rLHS(9,9)+=crLHS2*(crLHS119*r_DN(2,0) + crLHS121*r_DN(2,1) + crLHS124*crLHS263 + crLHS245*crLHS7 + crLHS273);
rLHS(9,10)+=crLHS2*(crLHS130*r_DN(2,0) + crLHS132*r_DN(2,1) + crLHS135*crLHS263 + crLHS274);
rLHS(9,11)+=crLHS271;
rLHS(9,12)+=crLHS276;
rLHS(9,13)+=crLHS2*(crLHS138*r_DN(2,0) + crLHS140*r_DN(2,1) + crLHS142*crLHS263 + crLHS278 + crLHS280);
rLHS(9,14)+=crLHS2*(crLHS148*r_DN(2,0) + crLHS150*r_DN(2,1) + crLHS152*crLHS263 + crLHS281);
rLHS(9,15)+=crLHS276;
rLHS(10,0)+=crLHS282;
rLHS(10,1)+=crLHS2*(crLHS133 + crLHS154*r_DN(2,1) + crLHS283*crLHS284 + crLHS77*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS157*r_DN(2,1) + crLHS169 + crLHS266 + crLHS285*crLHS95 + crLHS93*r_DN(2,0));
rLHS(10,3)+=crLHS282;
rLHS(10,4)+=crLHS286;
rLHS(10,5)+=crLHS2*(crLHS102*r_DN(2,0) + crLHS105*crLHS285 + crLHS160*r_DN(2,1) + crLHS221);
rLHS(10,6)+=crLHS2*(crLHS113*r_DN(2,0) + crLHS116*crLHS285 + crLHS162*r_DN(2,1) + crLHS236 + crLHS269);
rLHS(10,7)+=crLHS286;
rLHS(10,8)+=crLHS287;
rLHS(10,9)+=crLHS2*(crLHS121*r_DN(2,0) + crLHS124*crLHS285 + crLHS166*r_DN(2,1) + crLHS274);
rLHS(10,10)+=crLHS2*(crLHS132*r_DN(2,0) + crLHS135*crLHS285 + crLHS168*r_DN(2,1) + crLHS246*crLHS7 + crLHS273);
rLHS(10,11)+=crLHS287;
rLHS(10,12)+=crLHS289;
rLHS(10,13)+=crLHS2*(crLHS140*r_DN(2,0) + crLHS142*crLHS285 + crLHS172*r_DN(2,1) + crLHS291);
rLHS(10,14)+=crLHS2*(crLHS150*r_DN(2,0) + crLHS152*crLHS285 + crLHS174*r_DN(2,1) + crLHS280 + crLHS292);
rLHS(10,15)+=crLHS289;
rLHS(11,0)+=crLHS41;
rLHS(11,1)+=crLHS241;
rLHS(11,2)+=crLHS242;
rLHS(11,3)+=crLHS41;
rLHS(11,4)+=crLHS186;
rLHS(11,5)+=crLHS243;
rLHS(11,6)+=crLHS244;
rLHS(11,7)+=crLHS186;
rLHS(11,8)+=crLHS247;
rLHS(11,9)+=crLHS249;
rLHS(11,10)+=crLHS250;
rLHS(11,11)+=crLHS247;
rLHS(11,12)+=crLHS253;
rLHS(11,13)+=crLHS256;
rLHS(11,14)+=crLHS259;
rLHS(11,15)+=crLHS253;
rLHS(12,0)+=crLHS54;
rLHS(12,1)+=crLHS293;
rLHS(12,2)+=crLHS294;
rLHS(12,3)+=crLHS54;
rLHS(12,4)+=crLHS195;
rLHS(12,5)+=crLHS295;
rLHS(12,6)+=crLHS296;
rLHS(12,7)+=crLHS195;
rLHS(12,8)+=crLHS253;
rLHS(12,9)+=crLHS297;
rLHS(12,10)+=crLHS298;
rLHS(12,11)+=crLHS253;
rLHS(12,12)+=crLHS301;
rLHS(12,13)+=crLHS303;
rLHS(12,14)+=crLHS304;
rLHS(12,15)+=crLHS301;
rLHS(13,0)+=crLHS306;
rLHS(13,1)+=crLHS2*(crLHS145 + crLHS283*crLHS290 + crLHS309 + crLHS75*r_DN(3,0) + crLHS77*r_DN(3,1));
rLHS(13,2)+=crLHS2*(crLHS173 + crLHS310*crLHS95 + crLHS91*r_DN(3,0) + crLHS93*r_DN(3,1));
rLHS(13,3)+=crLHS306;
rLHS(13,4)+=crLHS311;
rLHS(13,5)+=crLHS2*(crLHS100*r_DN(3,0) + crLHS102*r_DN(3,1) + crLHS105*crLHS310 + crLHS225 + crLHS313);
rLHS(13,6)+=crLHS2*(crLHS111*r_DN(3,0) + crLHS113*r_DN(3,1) + crLHS116*crLHS310 + crLHS239);
rLHS(13,7)+=crLHS311;
rLHS(13,8)+=crLHS314;
rLHS(13,9)+=crLHS2*(crLHS119*r_DN(3,0) + crLHS121*r_DN(3,1) + crLHS124*crLHS310 + crLHS278 + crLHS316);
rLHS(13,10)+=crLHS2*(crLHS130*r_DN(3,0) + crLHS132*r_DN(3,1) + crLHS135*crLHS310 + crLHS291);
rLHS(13,11)+=crLHS314;
rLHS(13,12)+=crLHS318;
rLHS(13,13)+=crLHS2*(crLHS138*r_DN(3,0) + crLHS140*r_DN(3,1) + crLHS142*crLHS310 + crLHS299*crLHS7 + crLHS320);
rLHS(13,14)+=crLHS2*(crLHS148*r_DN(3,0) + crLHS150*r_DN(3,1) + crLHS152*crLHS310 + crLHS321);
rLHS(13,15)+=crLHS318;
rLHS(14,0)+=crLHS322;
rLHS(14,1)+=crLHS2*(crLHS151 + crLHS154*r_DN(3,1) + crLHS283*crLHS323 + crLHS77*r_DN(3,0));
rLHS(14,2)+=crLHS2*(crLHS157*r_DN(3,1) + crLHS175 + crLHS309 + crLHS324*crLHS95 + crLHS93*r_DN(3,0));
rLHS(14,3)+=crLHS322;
rLHS(14,4)+=crLHS325;
rLHS(14,5)+=crLHS2*(crLHS102*r_DN(3,0) + crLHS105*crLHS324 + crLHS160*r_DN(3,1) + crLHS228);
rLHS(14,6)+=crLHS2*(crLHS113*r_DN(3,0) + crLHS116*crLHS324 + crLHS162*r_DN(3,1) + crLHS240 + crLHS313);
rLHS(14,7)+=crLHS325;
rLHS(14,8)+=crLHS326;
rLHS(14,9)+=crLHS2*(crLHS121*r_DN(3,0) + crLHS124*crLHS324 + crLHS166*r_DN(3,1) + crLHS281);
rLHS(14,10)+=crLHS2*(crLHS132*r_DN(3,0) + crLHS135*crLHS324 + crLHS168*r_DN(3,1) + crLHS292 + crLHS316);
rLHS(14,11)+=crLHS326;
rLHS(14,12)+=crLHS327;
rLHS(14,13)+=crLHS2*(crLHS140*r_DN(3,0) + crLHS142*crLHS324 + crLHS172*r_DN(3,1) + crLHS321);
rLHS(14,14)+=crLHS2*(crLHS150*r_DN(3,0) + crLHS152*crLHS324 + crLHS174*r_DN(3,1) + crLHS300*crLHS7 + crLHS320);
rLHS(14,15)+=crLHS327;
rLHS(15,0)+=crLHS54;
rLHS(15,1)+=crLHS293;
rLHS(15,2)+=crLHS294;
rLHS(15,3)+=crLHS54;
rLHS(15,4)+=crLHS195;
rLHS(15,5)+=crLHS295;
rLHS(15,6)+=crLHS296;
rLHS(15,7)+=crLHS195;
rLHS(15,8)+=crLHS253;
rLHS(15,9)+=crLHS297;
rLHS(15,10)+=crLHS298;
rLHS(15,11)+=crLHS253;
rLHS(15,12)+=crLHS301;
rLHS(15,13)+=crLHS303;
rLHS(15,14)+=crLHS304;
rLHS(15,15)+=crLHS301;

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
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const double tau_c = CalculateTauPressure(rData);
    const double tau_u = CalculateTauVelocity(rData);
    const double tau_t = CalculateTauTemperature(rData);

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2];
const double crRHS4 = 1.0/crRHS3;
const double crRHS5 = crRHS4*p_th;
const double crRHS6 = -crRHS5*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2])) + dp_th_dt;
const double crRHS7 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2];
const double crRHS8 = crRHS7*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0));
const double crRHS9 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2];
const double crRHS10 = crRHS9*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1));
const double crRHS11 = crRHS5*(crRHS10 + crRHS8);
const double crRHS12 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0));
const double crRHS13 = r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0));
const double crRHS14 = r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS15 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crRHS16 = crRHS0*crRHS15;
const double crRHS17 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crRHS18 = crRHS17*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS19 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS20 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS21 = crRHS20*crRHS5;
const double crRHS22 = -crRHS21*(-crRHS12 - crRHS13 - crRHS14 - crRHS16 - crRHS18 + crRHS19) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS23 = p_th*tau_u;
const double crRHS24 = crRHS22*crRHS23;
const double crRHS25 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1));
const double crRHS26 = r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1));
const double crRHS27 = r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS28 = crRHS15*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS29 = crRHS1*crRHS17;
const double crRHS30 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS31 = -crRHS21*(-crRHS25 - crRHS26 - crRHS27 - crRHS28 - crRHS29 + crRHS30) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS32 = crRHS23*crRHS31;
const double crRHS33 = gauss_weight*gauss_weight;
const double crRHS34 = crRHS20*crRHS4;
const double crRHS35 = crRHS33*crRHS34;
const double crRHS36 = -crRHS35*(-crRHS11*r_N[0] + crRHS2*r_N[0] + crRHS24*r_DN(0,0) + crRHS32*r_DN(0,1) + crRHS6*r_N[0]);
const double crRHS37 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS38 = crRHS21*r_N[0];
const double crRHS39 = crRHS34*r_DN(0,0);
const double crRHS40 = crRHS12 + crRHS13 + crRHS14;
const double crRHS41 = crRHS16 + crRHS18;
const double crRHS42 = tau_c*(-crRHS10*crRHS5 + crRHS2 - crRHS5*crRHS8 + crRHS6);
const double crRHS43 = tau_u*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1));
const double crRHS44 = crRHS38*crRHS43;
const double crRHS45 = crRHS21*tau_u;
const double crRHS46 = crRHS45*(crRHS15*r_DN(0,0) + crRHS17*r_DN(0,1));
const double crRHS47 = crRHS20*(crRHS15*crRHS7 + crRHS17*crRHS9)*1.0/(crRHS3*crRHS3);
const double crRHS48 = crRHS47*r_N[0];
const double crRHS49 = crRHS34*r_DN(0,1);
const double crRHS50 = crRHS25 + crRHS26 + crRHS27;
const double crRHS51 = crRHS28 + crRHS29;
const double crRHS52 = -crRHS35*(-crRHS11*r_N[1] + crRHS2*r_N[1] + crRHS24*r_DN(1,0) + crRHS32*r_DN(1,1) + crRHS6*r_N[1]);
const double crRHS53 = crRHS21*r_N[1];
const double crRHS54 = crRHS34*r_DN(1,0);
const double crRHS55 = crRHS43*crRHS53;
const double crRHS56 = crRHS45*(crRHS15*r_DN(1,0) + crRHS17*r_DN(1,1));
const double crRHS57 = crRHS47*r_N[1];
const double crRHS58 = crRHS34*r_DN(1,1);
const double crRHS59 = -crRHS35*(-crRHS11*r_N[2] + crRHS2*r_N[2] + crRHS24*r_DN(2,0) + crRHS32*r_DN(2,1) + crRHS6*r_N[2]);
const double crRHS60 = crRHS21*r_N[2];
const double crRHS61 = crRHS34*r_DN(2,0);
const double crRHS62 = crRHS43*crRHS60;
const double crRHS63 = crRHS45*(crRHS15*r_DN(2,0) + crRHS17*r_DN(2,1));
const double crRHS64 = crRHS47*r_N[2];
const double crRHS65 = crRHS34*r_DN(2,1);
rRHS[0]+=crRHS36;
rRHS[1]+=-crRHS33*(-crRHS19*crRHS38 + crRHS2*crRHS39 + crRHS22*crRHS44 + crRHS22*crRHS46 - crRHS24*crRHS48 - crRHS37*r_DN(0,0) + crRHS38*crRHS40 + crRHS38*crRHS41 + crRHS39*crRHS42 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS33*(crRHS2*crRHS49 - crRHS30*crRHS38 + crRHS31*crRHS44 + crRHS31*crRHS46 - crRHS32*crRHS48 - crRHS37*r_DN(0,1) + crRHS38*crRHS50 + crRHS38*crRHS51 + crRHS42*crRHS49 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS36;
rRHS[4]+=crRHS52;
rRHS[5]+=-crRHS33*(-crRHS19*crRHS53 + crRHS2*crRHS54 + crRHS22*crRHS55 + crRHS22*crRHS56 - crRHS24*crRHS57 - crRHS37*r_DN(1,0) + crRHS40*crRHS53 + crRHS41*crRHS53 + crRHS42*crRHS54 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS33*(crRHS2*crRHS58 - crRHS30*crRHS53 + crRHS31*crRHS55 + crRHS31*crRHS56 - crRHS32*crRHS57 - crRHS37*r_DN(1,1) + crRHS42*crRHS58 + crRHS50*crRHS53 + crRHS51*crRHS53 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS52;
rRHS[8]+=crRHS59;
rRHS[9]+=-crRHS33*(-crRHS19*crRHS60 + crRHS2*crRHS61 + crRHS22*crRHS62 + crRHS22*crRHS63 - crRHS24*crRHS64 - crRHS37*r_DN(2,0) + crRHS40*crRHS60 + crRHS41*crRHS60 + crRHS42*crRHS61 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS33*(crRHS2*crRHS65 - crRHS30*crRHS60 + crRHS31*crRHS62 + crRHS31*crRHS63 - crRHS32*crRHS64 - crRHS37*r_DN(2,1) + crRHS42*crRHS65 + crRHS50*crRHS60 + crRHS51*crRHS60 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS59;

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
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const double tau_c = CalculateTauPressure(rData);
    const double tau_u = CalculateTauVelocity(rData);
    const double tau_t = CalculateTauTemperature(rData);

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0) + r_DN(3,0)*r_u(3,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1) + r_DN(3,1)*r_u(3,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*r_t_lin[0] + r_N[1]*r_t_lin[1] + r_N[2]*r_t_lin[2] + r_N[3]*r_t_lin[3];
const double crRHS4 = 1.0/crRHS3;
const double crRHS5 = crRHS4*p_th;
const double crRHS6 = -crRHS5*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]) + r_N[3]*(bdf0*r_t[3] + bdf1*r_t_n[3] + bdf2*r_t_nn[3])) + dp_th_dt;
const double crRHS7 = r_DN(0,0)*r_t_lin[0] + r_DN(1,0)*r_t_lin[1] + r_DN(2,0)*r_t_lin[2] + r_DN(3,0)*r_t_lin[3];
const double crRHS8 = crRHS7*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0) + r_N[3]*r_u(3,0));
const double crRHS9 = r_DN(0,1)*r_t_lin[0] + r_DN(1,1)*r_t_lin[1] + r_DN(2,1)*r_t_lin[2] + r_DN(3,1)*r_t_lin[3];
const double crRHS10 = crRHS9*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1) + r_N[3]*r_u(3,1));
const double crRHS11 = crRHS5*(crRHS10 + crRHS8);
const double crRHS12 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0));
const double crRHS13 = r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0));
const double crRHS14 = r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS15 = r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS16 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crRHS17 = crRHS0*crRHS16;
const double crRHS18 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crRHS19 = crRHS18*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS20 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS21 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS22 = crRHS21*crRHS5;
const double crRHS23 = -crRHS22*(-crRHS12 - crRHS13 - crRHS14 - crRHS15 - crRHS17 - crRHS19 + crRHS20) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS24 = p_th*tau_u;
const double crRHS25 = crRHS23*crRHS24;
const double crRHS26 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1));
const double crRHS27 = r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1));
const double crRHS28 = r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS29 = r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS30 = crRHS16*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS31 = crRHS1*crRHS18;
const double crRHS32 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS33 = -crRHS22*(-crRHS26 - crRHS27 - crRHS28 - crRHS29 - crRHS30 - crRHS31 + crRHS32) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS34 = crRHS24*crRHS33;
const double crRHS35 = gauss_weight*gauss_weight;
const double crRHS36 = crRHS21*crRHS4;
const double crRHS37 = crRHS35*crRHS36;
const double crRHS38 = -crRHS37*(-crRHS11*r_N[0] + crRHS2*r_N[0] + crRHS25*r_DN(0,0) + crRHS34*r_DN(0,1) + crRHS6*r_N[0]);
const double crRHS39 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS40 = crRHS22*r_N[0];
const double crRHS41 = crRHS36*r_DN(0,0);
const double crRHS42 = crRHS12 + crRHS13 + crRHS14 + crRHS15;
const double crRHS43 = crRHS17 + crRHS19;
const double crRHS44 = tau_c*(-crRHS10*crRHS5 + crRHS2 - crRHS5*crRHS8 + crRHS6);
const double crRHS45 = tau_u*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1));
const double crRHS46 = crRHS40*crRHS45;
const double crRHS47 = crRHS22*tau_u;
const double crRHS48 = crRHS47*(crRHS16*r_DN(0,0) + crRHS18*r_DN(0,1));
const double crRHS49 = crRHS21*(crRHS16*crRHS7 + crRHS18*crRHS9)*1.0/(crRHS3*crRHS3);
const double crRHS50 = crRHS49*r_N[0];
const double crRHS51 = crRHS36*r_DN(0,1);
const double crRHS52 = crRHS26 + crRHS27 + crRHS28 + crRHS29;
const double crRHS53 = crRHS30 + crRHS31;
const double crRHS54 = -crRHS37*(-crRHS11*r_N[1] + crRHS2*r_N[1] + crRHS25*r_DN(1,0) + crRHS34*r_DN(1,1) + crRHS6*r_N[1]);
const double crRHS55 = crRHS22*r_N[1];
const double crRHS56 = crRHS36*r_DN(1,0);
const double crRHS57 = crRHS45*crRHS55;
const double crRHS58 = crRHS47*(crRHS16*r_DN(1,0) + crRHS18*r_DN(1,1));
const double crRHS59 = crRHS49*r_N[1];
const double crRHS60 = crRHS36*r_DN(1,1);
const double crRHS61 = -crRHS37*(-crRHS11*r_N[2] + crRHS2*r_N[2] + crRHS25*r_DN(2,0) + crRHS34*r_DN(2,1) + crRHS6*r_N[2]);
const double crRHS62 = crRHS22*r_N[2];
const double crRHS63 = crRHS36*r_DN(2,0);
const double crRHS64 = crRHS45*crRHS62;
const double crRHS65 = crRHS47*(crRHS16*r_DN(2,0) + crRHS18*r_DN(2,1));
const double crRHS66 = crRHS49*r_N[2];
const double crRHS67 = crRHS36*r_DN(2,1);
const double crRHS68 = -crRHS37*(-crRHS11*r_N[3] + crRHS2*r_N[3] + crRHS25*r_DN(3,0) + crRHS34*r_DN(3,1) + crRHS6*r_N[3]);
const double crRHS69 = crRHS22*r_N[3];
const double crRHS70 = crRHS36*r_DN(3,0);
const double crRHS71 = crRHS45*crRHS69;
const double crRHS72 = crRHS47*(crRHS16*r_DN(3,0) + crRHS18*r_DN(3,1));
const double crRHS73 = crRHS49*r_N[3];
const double crRHS74 = crRHS36*r_DN(3,1);
rRHS[0]+=crRHS38;
rRHS[1]+=-crRHS35*(crRHS2*crRHS41 - crRHS20*crRHS40 + crRHS23*crRHS46 + crRHS23*crRHS48 - crRHS25*crRHS50 - crRHS39*r_DN(0,0) + crRHS40*crRHS42 + crRHS40*crRHS43 + crRHS41*crRHS44 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS35*(crRHS2*crRHS51 - crRHS32*crRHS40 + crRHS33*crRHS46 + crRHS33*crRHS48 - crRHS34*crRHS50 - crRHS39*r_DN(0,1) + crRHS40*crRHS52 + crRHS40*crRHS53 + crRHS44*crRHS51 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS38;
rRHS[4]+=crRHS54;
rRHS[5]+=-crRHS35*(crRHS2*crRHS56 - crRHS20*crRHS55 + crRHS23*crRHS57 + crRHS23*crRHS58 - crRHS25*crRHS59 - crRHS39*r_DN(1,0) + crRHS42*crRHS55 + crRHS43*crRHS55 + crRHS44*crRHS56 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS35*(crRHS2*crRHS60 - crRHS32*crRHS55 + crRHS33*crRHS57 + crRHS33*crRHS58 - crRHS34*crRHS59 - crRHS39*r_DN(1,1) + crRHS44*crRHS60 + crRHS52*crRHS55 + crRHS53*crRHS55 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS54;
rRHS[8]+=crRHS61;
rRHS[9]+=-crRHS35*(crRHS2*crRHS63 - crRHS20*crRHS62 + crRHS23*crRHS64 + crRHS23*crRHS65 - crRHS25*crRHS66 - crRHS39*r_DN(2,0) + crRHS42*crRHS62 + crRHS43*crRHS62 + crRHS44*crRHS63 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS35*(crRHS2*crRHS67 - crRHS32*crRHS62 + crRHS33*crRHS64 + crRHS33*crRHS65 - crRHS34*crRHS66 - crRHS39*r_DN(2,1) + crRHS44*crRHS67 + crRHS52*crRHS62 + crRHS53*crRHS62 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS61;
rRHS[12]+=crRHS68;
rRHS[13]+=-crRHS35*(crRHS2*crRHS70 - crRHS20*crRHS69 + crRHS23*crRHS71 + crRHS23*crRHS72 - crRHS25*crRHS73 - crRHS39*r_DN(3,0) + crRHS42*crRHS69 + crRHS43*crRHS69 + crRHS44*crRHS70 + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-crRHS35*(crRHS2*crRHS74 - crRHS32*crRHS69 + crRHS33*crRHS71 + crRHS33*crRHS72 - crRHS34*crRHS73 - crRHS39*r_DN(3,1) + crRHS44*crRHS74 + crRHS52*crRHS69 + crRHS53*crRHS69 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=crRHS68;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations


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