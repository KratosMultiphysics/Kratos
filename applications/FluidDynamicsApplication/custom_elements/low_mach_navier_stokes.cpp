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
    const auto& r_rho_n = rData.DensityOldStep1;
    const auto& r_rho_nn = rData.DensityOldStep2;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double alpha = rData.ThermalExpansionCoefficient;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const BoundedMatrix<double,2,3> u_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = gamma*p_th*1.0/crLHS5*1.0/c_p;
const double crLHS7 = crLHS4*crLHS6;
const double crLHS8 = crLHS2*crLHS7;
const double crLHS9 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crLHS10 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crLHS11 = stab_c2*sqrt(crLHS10*crLHS10 + crLHS9*crLHS9);
const double crLHS12 = 1.0*1.0/(crLHS11 + mu*stab_c1*1.0/(h*h));
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = crLHS13*(crLHS0 + crLHS1);
const double crLHS15 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crLHS16 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS17 = crLHS7*r_DN(0,0);
const double crLHS18 = bdf0*r_N[0];
const double crLHS19 = crLHS10*r_DN(0,1) + crLHS9*r_DN(0,0);
const double crLHS20 = crLHS18 + crLHS19;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = crLHS8*(-crLHS15*crLHS16 + crLHS17*crLHS21 + r_DN(0,0)*r_N[0]);
const double crLHS23 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crLHS24 = crLHS7*r_DN(0,1);
const double crLHS25 = crLHS8*(-crLHS16*crLHS23 + crLHS21*crLHS24 + r_DN(0,1)*r_N[0]);
const double crLHS26 = r_DN(0,0)*r_DN(1,0);
const double crLHS27 = r_DN(0,1)*r_DN(1,1);
const double crLHS28 = crLHS13*(crLHS26 + crLHS27);
const double crLHS29 = r_DN(1,0)*r_N[0];
const double crLHS30 = crLHS4*r_N[0];
const double crLHS31 = crLHS15*crLHS30;
const double crLHS32 = -crLHS31*r_N[1];
const double crLHS33 = bdf0*r_N[1];
const double crLHS34 = crLHS10*r_DN(1,1) + crLHS9*r_DN(1,0);
const double crLHS35 = crLHS33 + crLHS34;
const double crLHS36 = crLHS12*crLHS35;
const double crLHS37 = crLHS8*(crLHS17*crLHS36 + crLHS29 + crLHS32);
const double crLHS38 = r_DN(1,1)*r_N[0];
const double crLHS39 = crLHS23*crLHS30;
const double crLHS40 = -crLHS39*r_N[1];
const double crLHS41 = crLHS8*(crLHS24*crLHS36 + crLHS38 + crLHS40);
const double crLHS42 = r_DN(0,0)*r_DN(2,0);
const double crLHS43 = r_DN(0,1)*r_DN(2,1);
const double crLHS44 = crLHS13*(crLHS42 + crLHS43);
const double crLHS45 = r_DN(2,0)*r_N[0];
const double crLHS46 = -crLHS31*r_N[2];
const double crLHS47 = crLHS10*r_DN(2,1) + crLHS9*r_DN(2,0);
const double crLHS48 = bdf0*r_N[2] + crLHS47;
const double crLHS49 = crLHS12*crLHS48;
const double crLHS50 = crLHS8*(crLHS17*crLHS49 + crLHS45 + crLHS46);
const double crLHS51 = r_DN(2,1)*r_N[0];
const double crLHS52 = -crLHS39*r_N[2];
const double crLHS53 = crLHS8*(crLHS24*crLHS49 + crLHS51 + crLHS52);
const double crLHS54 = crLHS30*crLHS6;
const double crLHS55 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1);
const double crLHS56 = crLHS12*crLHS55;
const double crLHS57 = crLHS12*crLHS7;
const double crLHS58 = crLHS19*crLHS57;
const double crLHS59 = 1.0/(crLHS3*crLHS3);
const double crLHS60 = crLHS59*r_N[0];
const double crLHS61 = crLHS10*crLHS23 + crLHS15*crLHS9;
const double crLHS62 = crLHS12*crLHS6*crLHS61;
const double crLHS63 = crLHS2*(-crLHS54*crLHS56 - crLHS58 + crLHS60*crLHS62 + r_N[0]);
const double crLHS64 = -crLHS63*r_DN(0,0);
const double crLHS65 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS66 = r_C(0,2)*r_DN(0,0);
const double crLHS67 = crLHS66 + r_C(2,2)*r_DN(0,1);
const double crLHS68 = -crLHS31 + r_DN(0,0);
const double crLHS69 = crLHS11*h*1.0/stab_c1 + mu*1.0/rho_lin;
const double crLHS70 = crLHS17*crLHS69;
const double crLHS71 = bdf0*crLHS6;
const double crLHS72 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS73 = crLHS59*crLHS72;
const double crLHS74 = crLHS19*crLHS73;
const double crLHS75 = crLHS21*crLHS72;
const double crLHS76 = crLHS55*crLHS60;
const double crLHS77 = crLHS61*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS78 = crLHS77*r_N[0];
const double crLHS79 = crLHS16*crLHS71 + crLHS19*crLHS54 + crLHS21*crLHS74 + crLHS75*crLHS76 - crLHS75*crLHS78;
const double crLHS80 = crLHS66 + r_C(0,1)*r_DN(0,1);
const double crLHS81 = r_C(1,2)*r_DN(0,1);
const double crLHS82 = crLHS81 + r_C(2,2)*r_DN(0,0);
const double crLHS83 = crLHS17*r_DN(0,1);
const double crLHS84 = -crLHS39 + r_DN(0,1);
const double crLHS85 = r_DN(0,0)*r_N[1];
const double crLHS86 = crLHS55*crLHS57;
const double crLHS87 = crLHS7*r_DN(1,0);
const double crLHS88 = crLHS12*crLHS19;
const double crLHS89 = crLHS59*crLHS62;
const double crLHS90 = -crLHS2*(-crLHS29*crLHS86 + crLHS29*crLHS89 + crLHS85 - crLHS87*crLHS88);
const double crLHS91 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS92 = r_C(0,2)*r_DN(1,0);
const double crLHS93 = crLHS92 + r_C(2,2)*r_DN(1,1);
const double crLHS94 = crLHS4*r_N[1];
const double crLHS95 = crLHS15*crLHS94;
const double crLHS96 = -crLHS95 + r_DN(1,0);
const double crLHS97 = crLHS6*crLHS94;
const double crLHS98 = crLHS18*crLHS97;
const double crLHS99 = crLHS26*crLHS7 + crLHS98;
const double crLHS100 = crLHS36*crLHS72;
const double crLHS101 = crLHS100*crLHS76 - crLHS100*crLHS78 + crLHS34*crLHS54 + crLHS36*crLHS74;
const double crLHS102 = crLHS92 + r_C(0,1)*r_DN(1,1);
const double crLHS103 = r_C(1,2)*r_DN(1,1);
const double crLHS104 = crLHS103 + r_C(2,2)*r_DN(1,0);
const double crLHS105 = crLHS17*r_DN(1,1);
const double crLHS106 = crLHS23*crLHS94;
const double crLHS107 = -crLHS106 + r_DN(1,1);
const double crLHS108 = r_DN(0,0)*r_N[2];
const double crLHS109 = crLHS7*r_DN(2,0);
const double crLHS110 = -crLHS2*(crLHS108 - crLHS109*crLHS88 - crLHS45*crLHS86 + crLHS45*crLHS89);
const double crLHS111 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS112 = r_C(0,2)*r_DN(2,0);
const double crLHS113 = crLHS112 + r_C(2,2)*r_DN(2,1);
const double crLHS114 = crLHS4*r_N[2];
const double crLHS115 = -crLHS114*crLHS15 + r_DN(2,0);
const double crLHS116 = crLHS114*crLHS6;
const double crLHS117 = crLHS116*crLHS18;
const double crLHS118 = crLHS117 + crLHS42*crLHS7;
const double crLHS119 = crLHS49*crLHS72;
const double crLHS120 = crLHS119*crLHS76 - crLHS119*crLHS78 + crLHS47*crLHS54 + crLHS49*crLHS74;
const double crLHS121 = crLHS112 + r_C(0,1)*r_DN(2,1);
const double crLHS122 = r_C(1,2)*r_DN(2,1);
const double crLHS123 = crLHS122 + r_C(2,2)*r_DN(2,0);
const double crLHS124 = crLHS17*r_DN(2,1);
const double crLHS125 = -crLHS114*crLHS23 + r_DN(2,1);
const double crLHS126 = -crLHS63*r_DN(0,1);
const double crLHS127 = crLHS81 + r_C(0,1)*r_DN(0,0);
const double crLHS128 = crLHS24*crLHS69;
const double crLHS129 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS130 = r_DN(0,1)*r_N[1];
const double crLHS131 = -crLHS2*(crLHS130 - crLHS38*crLHS86 + crLHS38*crLHS89 - crLHS58*r_DN(1,1));
const double crLHS132 = crLHS103 + r_C(0,1)*r_DN(1,0);
const double crLHS133 = crLHS24*r_DN(1,0);
const double crLHS134 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS135 = crLHS27*crLHS7 + crLHS98;
const double crLHS136 = r_DN(0,1)*r_N[2];
const double crLHS137 = -crLHS2*(crLHS136 - crLHS51*crLHS86 + crLHS51*crLHS89 - crLHS58*r_DN(2,1));
const double crLHS138 = crLHS122 + r_C(0,1)*r_DN(2,0);
const double crLHS139 = crLHS24*r_DN(2,0);
const double crLHS140 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS141 = crLHS117 + crLHS43*crLHS7;
const double crLHS142 = crLHS8*(crLHS21*crLHS87 + crLHS32 + crLHS85);
const double crLHS143 = crLHS57*r_DN(1,1);
const double crLHS144 = crLHS8*(crLHS130 + crLHS143*crLHS20 + crLHS40);
const double crLHS145 = r_DN(1,0)*r_DN(1,0);
const double crLHS146 = r_DN(1,1)*r_DN(1,1);
const double crLHS147 = crLHS13*(crLHS145 + crLHS146);
const double crLHS148 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS149 = crLHS8*(-crLHS148*crLHS15 + crLHS36*crLHS87 + r_DN(1,0)*r_N[1]);
const double crLHS150 = crLHS8*(crLHS143*crLHS35 - crLHS148*crLHS23 + r_DN(1,1)*r_N[1]);
const double crLHS151 = r_DN(1,0)*r_DN(2,0);
const double crLHS152 = r_DN(1,1)*r_DN(2,1);
const double crLHS153 = crLHS13*(crLHS151 + crLHS152);
const double crLHS154 = r_DN(2,0)*r_N[1];
const double crLHS155 = -crLHS95*r_N[2];
const double crLHS156 = crLHS8*(crLHS154 + crLHS155 + crLHS49*crLHS87);
const double crLHS157 = r_DN(2,1)*r_N[1];
const double crLHS158 = -crLHS106*r_N[2];
const double crLHS159 = crLHS8*(crLHS143*crLHS48 + crLHS157 + crLHS158);
const double crLHS160 = crLHS12*crLHS34;
const double crLHS161 = crLHS2*(crLHS160*crLHS17 - crLHS29 + crLHS85*crLHS86 - crLHS85*crLHS89);
const double crLHS162 = crLHS69*crLHS87;
const double crLHS163 = crLHS34*crLHS73;
const double crLHS164 = crLHS55*crLHS73;
const double crLHS165 = crLHS164*r_N[1];
const double crLHS166 = crLHS77*r_N[1];
const double crLHS167 = crLHS163*crLHS21 + crLHS165*crLHS21 - crLHS166*crLHS75 + crLHS19*crLHS97;
const double crLHS168 = crLHS34*crLHS57;
const double crLHS169 = crLHS2*(-crLHS168 - crLHS56*crLHS97 + crLHS89*r_N[1] + r_N[1]);
const double crLHS170 = -crLHS169*r_DN(1,0);
const double crLHS171 = -crLHS100*crLHS166 + crLHS148*crLHS71 + crLHS163*crLHS36 + crLHS165*crLHS36 + crLHS34*crLHS97;
const double crLHS172 = crLHS87*r_DN(1,1);
const double crLHS173 = r_DN(1,0)*r_N[2];
const double crLHS174 = -crLHS2*(-crLHS109*crLHS160 - crLHS154*crLHS86 + crLHS154*crLHS89 + crLHS173);
const double crLHS175 = crLHS116*crLHS33;
const double crLHS176 = crLHS151*crLHS7 + crLHS175;
const double crLHS177 = -crLHS119*crLHS166 + crLHS163*crLHS49 + crLHS165*crLHS49 + crLHS47*crLHS97;
const double crLHS178 = crLHS87*r_DN(2,1);
const double crLHS179 = crLHS2*(crLHS130*crLHS86 - crLHS130*crLHS89 + crLHS160*crLHS24 - crLHS38);
const double crLHS180 = crLHS68*crLHS69;
const double crLHS181 = crLHS7*r_DN(1,1);
const double crLHS182 = crLHS181*crLHS69;
const double crLHS183 = -crLHS169*r_DN(1,1);
const double crLHS184 = r_DN(1,1)*r_N[2];
const double crLHS185 = -crLHS2*(-crLHS157*crLHS86 + crLHS157*crLHS89 - crLHS168*r_DN(2,1) + crLHS184);
const double crLHS186 = crLHS109*r_DN(1,1);
const double crLHS187 = crLHS152*crLHS7 + crLHS175;
const double crLHS188 = crLHS8*(crLHS108 + crLHS109*crLHS21 + crLHS46);
const double crLHS189 = crLHS57*r_DN(2,1);
const double crLHS190 = crLHS8*(crLHS136 + crLHS189*crLHS20 + crLHS52);
const double crLHS191 = crLHS8*(crLHS109*crLHS36 + crLHS155 + crLHS173);
const double crLHS192 = crLHS8*(crLHS158 + crLHS184 + crLHS189*crLHS35);
const double crLHS193 = r_DN(2,0)*r_DN(2,0);
const double crLHS194 = r_DN(2,1)*r_DN(2,1);
const double crLHS195 = crLHS13*(crLHS193 + crLHS194);
const double crLHS196 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS197 = crLHS8*(crLHS109*crLHS49 - crLHS15*crLHS196 + r_DN(2,0)*r_N[2]);
const double crLHS198 = crLHS8*(crLHS189*crLHS48 - crLHS196*crLHS23 + r_DN(2,1)*r_N[2]);
const double crLHS199 = crLHS12*crLHS47;
const double crLHS200 = crLHS2*(crLHS108*crLHS86 - crLHS108*crLHS89 + crLHS17*crLHS199 - crLHS45);
const double crLHS201 = crLHS47*crLHS73;
const double crLHS202 = crLHS164*r_N[2];
const double crLHS203 = crLHS77*r_N[2];
const double crLHS204 = crLHS116*crLHS19 + crLHS201*crLHS21 + crLHS202*crLHS21 - crLHS203*crLHS75;
const double crLHS205 = crLHS109*crLHS69;
const double crLHS206 = crLHS2*(-crLHS154 + crLHS173*crLHS86 - crLHS173*crLHS89 + crLHS199*crLHS87);
const double crLHS207 = -crLHS100*crLHS203 + crLHS116*crLHS34 + crLHS201*crLHS36 + crLHS202*crLHS36;
const double crLHS208 = crLHS47*crLHS57;
const double crLHS209 = crLHS2*(-crLHS116*crLHS56 - crLHS208 + crLHS89*r_N[2] + r_N[2]);
const double crLHS210 = -crLHS209*r_DN(2,0);
const double crLHS211 = crLHS116*crLHS47 - crLHS119*crLHS203 + crLHS196*crLHS71 + crLHS201*crLHS49 + crLHS202*crLHS49;
const double crLHS212 = crLHS109*r_DN(2,1);
const double crLHS213 = crLHS2*(crLHS136*crLHS86 - crLHS136*crLHS89 + crLHS199*crLHS24 - crLHS51);
const double crLHS214 = crLHS7*r_DN(2,1);
const double crLHS215 = crLHS214*crLHS69;
const double crLHS216 = crLHS2*(-crLHS157 + crLHS184*crLHS86 - crLHS184*crLHS89 + crLHS208*r_DN(1,1));
const double crLHS217 = -crLHS209*r_DN(2,1);
rLHS(0,0)+=crLHS14;
rLHS(0,1)+=crLHS22;
rLHS(0,2)+=crLHS25;
rLHS(0,3)+=crLHS14;
rLHS(0,4)+=crLHS28;
rLHS(0,5)+=crLHS37;
rLHS(0,6)+=crLHS41;
rLHS(0,7)+=crLHS28;
rLHS(0,8)+=crLHS44;
rLHS(0,9)+=crLHS50;
rLHS(0,10)+=crLHS53;
rLHS(0,11)+=crLHS44;
rLHS(1,0)+=crLHS64;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS7 + crLHS65*r_DN(0,0) + crLHS67*r_DN(0,1) + crLHS68*crLHS70 + crLHS79);
rLHS(1,2)+=crLHS2*(crLHS70*crLHS84 + crLHS80*r_DN(0,0) + crLHS82*r_DN(0,1) + crLHS83);
rLHS(1,3)+=crLHS64;
rLHS(1,4)+=crLHS90;
rLHS(1,5)+=crLHS2*(crLHS101 + crLHS70*crLHS96 + crLHS91*r_DN(0,0) + crLHS93*r_DN(0,1) + crLHS99);
rLHS(1,6)+=crLHS2*(crLHS102*r_DN(0,0) + crLHS104*r_DN(0,1) + crLHS105 + crLHS107*crLHS70);
rLHS(1,7)+=crLHS90;
rLHS(1,8)+=crLHS110;
rLHS(1,9)+=crLHS2*(crLHS111*r_DN(0,0) + crLHS113*r_DN(0,1) + crLHS115*crLHS70 + crLHS118 + crLHS120);
rLHS(1,10)+=crLHS2*(crLHS121*r_DN(0,0) + crLHS123*r_DN(0,1) + crLHS124 + crLHS125*crLHS70);
rLHS(1,11)+=crLHS110;
rLHS(2,0)+=crLHS126;
rLHS(2,1)+=crLHS2*(crLHS127*r_DN(0,1) + crLHS128*crLHS68 + crLHS67*r_DN(0,0) + crLHS83);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS7 + crLHS128*crLHS84 + crLHS129*r_DN(0,1) + crLHS79 + crLHS82*r_DN(0,0));
rLHS(2,3)+=crLHS126;
rLHS(2,4)+=crLHS131;
rLHS(2,5)+=crLHS2*(crLHS128*crLHS96 + crLHS132*r_DN(0,1) + crLHS133 + crLHS93*r_DN(0,0));
rLHS(2,6)+=crLHS2*(crLHS101 + crLHS104*r_DN(0,0) + crLHS107*crLHS128 + crLHS134*r_DN(0,1) + crLHS135);
rLHS(2,7)+=crLHS131;
rLHS(2,8)+=crLHS137;
rLHS(2,9)+=crLHS2*(crLHS113*r_DN(0,0) + crLHS115*crLHS128 + crLHS138*r_DN(0,1) + crLHS139);
rLHS(2,10)+=crLHS2*(crLHS120 + crLHS123*r_DN(0,0) + crLHS125*crLHS128 + crLHS140*r_DN(0,1) + crLHS141);
rLHS(2,11)+=crLHS137;
rLHS(3,0)+=crLHS14;
rLHS(3,1)+=crLHS22;
rLHS(3,2)+=crLHS25;
rLHS(3,3)+=crLHS14;
rLHS(3,4)+=crLHS28;
rLHS(3,5)+=crLHS37;
rLHS(3,6)+=crLHS41;
rLHS(3,7)+=crLHS28;
rLHS(3,8)+=crLHS44;
rLHS(3,9)+=crLHS50;
rLHS(3,10)+=crLHS53;
rLHS(3,11)+=crLHS44;
rLHS(4,0)+=crLHS28;
rLHS(4,1)+=crLHS142;
rLHS(4,2)+=crLHS144;
rLHS(4,3)+=crLHS28;
rLHS(4,4)+=crLHS147;
rLHS(4,5)+=crLHS149;
rLHS(4,6)+=crLHS150;
rLHS(4,7)+=crLHS147;
rLHS(4,8)+=crLHS153;
rLHS(4,9)+=crLHS156;
rLHS(4,10)+=crLHS159;
rLHS(4,11)+=crLHS153;
rLHS(5,0)+=crLHS161;
rLHS(5,1)+=crLHS2*(crLHS162*crLHS68 + crLHS167 + crLHS65*r_DN(1,0) + crLHS67*r_DN(1,1) + crLHS99);
rLHS(5,2)+=crLHS2*(crLHS133 + crLHS162*crLHS84 + crLHS80*r_DN(1,0) + crLHS82*r_DN(1,1));
rLHS(5,3)+=crLHS161;
rLHS(5,4)+=crLHS170;
rLHS(5,5)+=crLHS2*(crLHS145*crLHS7 + crLHS162*crLHS96 + crLHS171 + crLHS91*r_DN(1,0) + crLHS93*r_DN(1,1));
rLHS(5,6)+=crLHS2*(crLHS102*r_DN(1,0) + crLHS104*r_DN(1,1) + crLHS107*crLHS162 + crLHS172);
rLHS(5,7)+=crLHS170;
rLHS(5,8)+=crLHS174;
rLHS(5,9)+=crLHS2*(crLHS111*r_DN(1,0) + crLHS113*r_DN(1,1) + crLHS115*crLHS162 + crLHS176 + crLHS177);
rLHS(5,10)+=crLHS2*(crLHS121*r_DN(1,0) + crLHS123*r_DN(1,1) + crLHS125*crLHS162 + crLHS178);
rLHS(5,11)+=crLHS174;
rLHS(6,0)+=crLHS179;
rLHS(6,1)+=crLHS2*(crLHS105 + crLHS127*r_DN(1,1) + crLHS180*crLHS181 + crLHS67*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS129*r_DN(1,1) + crLHS135 + crLHS167 + crLHS182*crLHS84 + crLHS82*r_DN(1,0));
rLHS(6,3)+=crLHS179;
rLHS(6,4)+=crLHS183;
rLHS(6,5)+=crLHS2*(crLHS132*r_DN(1,1) + crLHS172 + crLHS182*crLHS96 + crLHS93*r_DN(1,0));
rLHS(6,6)+=crLHS2*(crLHS104*r_DN(1,0) + crLHS107*crLHS182 + crLHS134*r_DN(1,1) + crLHS146*crLHS7 + crLHS171);
rLHS(6,7)+=crLHS183;
rLHS(6,8)+=crLHS185;
rLHS(6,9)+=crLHS2*(crLHS113*r_DN(1,0) + crLHS115*crLHS182 + crLHS138*r_DN(1,1) + crLHS186);
rLHS(6,10)+=crLHS2*(crLHS123*r_DN(1,0) + crLHS125*crLHS182 + crLHS140*r_DN(1,1) + crLHS177 + crLHS187);
rLHS(6,11)+=crLHS185;
rLHS(7,0)+=crLHS28;
rLHS(7,1)+=crLHS142;
rLHS(7,2)+=crLHS144;
rLHS(7,3)+=crLHS28;
rLHS(7,4)+=crLHS147;
rLHS(7,5)+=crLHS149;
rLHS(7,6)+=crLHS150;
rLHS(7,7)+=crLHS147;
rLHS(7,8)+=crLHS153;
rLHS(7,9)+=crLHS156;
rLHS(7,10)+=crLHS159;
rLHS(7,11)+=crLHS153;
rLHS(8,0)+=crLHS44;
rLHS(8,1)+=crLHS188;
rLHS(8,2)+=crLHS190;
rLHS(8,3)+=crLHS44;
rLHS(8,4)+=crLHS153;
rLHS(8,5)+=crLHS191;
rLHS(8,6)+=crLHS192;
rLHS(8,7)+=crLHS153;
rLHS(8,8)+=crLHS195;
rLHS(8,9)+=crLHS197;
rLHS(8,10)+=crLHS198;
rLHS(8,11)+=crLHS195;
rLHS(9,0)+=crLHS200;
rLHS(9,1)+=crLHS2*(crLHS109*crLHS180 + crLHS118 + crLHS204 + crLHS65*r_DN(2,0) + crLHS67*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS139 + crLHS205*crLHS84 + crLHS80*r_DN(2,0) + crLHS82*r_DN(2,1));
rLHS(9,3)+=crLHS200;
rLHS(9,4)+=crLHS206;
rLHS(9,5)+=crLHS2*(crLHS176 + crLHS205*crLHS96 + crLHS207 + crLHS91*r_DN(2,0) + crLHS93*r_DN(2,1));
rLHS(9,6)+=crLHS2*(crLHS102*r_DN(2,0) + crLHS104*r_DN(2,1) + crLHS107*crLHS205 + crLHS186);
rLHS(9,7)+=crLHS206;
rLHS(9,8)+=crLHS210;
rLHS(9,9)+=crLHS2*(crLHS111*r_DN(2,0) + crLHS113*r_DN(2,1) + crLHS115*crLHS205 + crLHS193*crLHS7 + crLHS211);
rLHS(9,10)+=crLHS2*(crLHS121*r_DN(2,0) + crLHS123*r_DN(2,1) + crLHS125*crLHS205 + crLHS212);
rLHS(9,11)+=crLHS210;
rLHS(10,0)+=crLHS213;
rLHS(10,1)+=crLHS2*(crLHS124 + crLHS127*r_DN(2,1) + crLHS180*crLHS214 + crLHS67*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS129*r_DN(2,1) + crLHS141 + crLHS204 + crLHS215*crLHS84 + crLHS82*r_DN(2,0));
rLHS(10,3)+=crLHS213;
rLHS(10,4)+=crLHS216;
rLHS(10,5)+=crLHS2*(crLHS132*r_DN(2,1) + crLHS178 + crLHS215*crLHS96 + crLHS93*r_DN(2,0));
rLHS(10,6)+=crLHS2*(crLHS104*r_DN(2,0) + crLHS107*crLHS215 + crLHS134*r_DN(2,1) + crLHS187 + crLHS207);
rLHS(10,7)+=crLHS216;
rLHS(10,8)+=crLHS217;
rLHS(10,9)+=crLHS2*(crLHS113*r_DN(2,0) + crLHS115*crLHS215 + crLHS138*r_DN(2,1) + crLHS212);
rLHS(10,10)+=crLHS2*(crLHS123*r_DN(2,0) + crLHS125*crLHS215 + crLHS140*r_DN(2,1) + crLHS194*crLHS7 + crLHS211);
rLHS(10,11)+=crLHS217;
rLHS(11,0)+=crLHS44;
rLHS(11,1)+=crLHS188;
rLHS(11,2)+=crLHS190;
rLHS(11,3)+=crLHS44;
rLHS(11,4)+=crLHS153;
rLHS(11,5)+=crLHS191;
rLHS(11,6)+=crLHS192;
rLHS(11,7)+=crLHS153;
rLHS(11,8)+=crLHS195;
rLHS(11,9)+=crLHS197;
rLHS(11,10)+=crLHS198;
rLHS(11,11)+=crLHS195;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const auto& r_rho_n = rData.DensityOldStep1;
    const auto& r_rho_nn = rData.DensityOldStep2;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double alpha = rData.ThermalExpansionCoefficient;

    // Material response parameters
    const auto& r_C = rData.C;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const BoundedMatrix<double,2,4> u_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add LHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crLHS0 = r_DN(0,0)*r_DN(0,0);
const double crLHS1 = r_DN(0,1)*r_DN(0,1);
const double crLHS2 = gauss_weight*gauss_weight;
const double crLHS3 = r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3];
const double crLHS4 = 1.0/crLHS3;
const double crLHS5 = gamma - 1.0;
const double crLHS6 = gamma*p_th*1.0/crLHS5*1.0/c_p;
const double crLHS7 = crLHS4*crLHS6;
const double crLHS8 = crLHS2*crLHS7;
const double crLHS9 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crLHS10 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crLHS11 = stab_c2*sqrt(crLHS10*crLHS10 + crLHS9*crLHS9);
const double crLHS12 = 1.0*1.0/(crLHS11 + mu*stab_c1*1.0/(h*h));
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = crLHS13*(crLHS0 + crLHS1);
const double crLHS15 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crLHS16 = crLHS4*(r_N[0]*r_N[0]);
const double crLHS17 = crLHS7*r_DN(0,0);
const double crLHS18 = bdf0*r_N[0];
const double crLHS19 = crLHS10*r_DN(0,1) + crLHS9*r_DN(0,0);
const double crLHS20 = crLHS18 + crLHS19;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = crLHS8*(-crLHS15*crLHS16 + crLHS17*crLHS21 + r_DN(0,0)*r_N[0]);
const double crLHS23 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crLHS24 = crLHS7*r_DN(0,1);
const double crLHS25 = crLHS8*(-crLHS16*crLHS23 + crLHS21*crLHS24 + r_DN(0,1)*r_N[0]);
const double crLHS26 = r_DN(0,0)*r_DN(1,0);
const double crLHS27 = r_DN(0,1)*r_DN(1,1);
const double crLHS28 = crLHS13*(crLHS26 + crLHS27);
const double crLHS29 = r_DN(1,0)*r_N[0];
const double crLHS30 = crLHS4*r_N[0];
const double crLHS31 = crLHS15*crLHS30;
const double crLHS32 = -crLHS31*r_N[1];
const double crLHS33 = bdf0*r_N[1];
const double crLHS34 = crLHS10*r_DN(1,1) + crLHS9*r_DN(1,0);
const double crLHS35 = crLHS33 + crLHS34;
const double crLHS36 = crLHS12*crLHS35;
const double crLHS37 = crLHS8*(crLHS17*crLHS36 + crLHS29 + crLHS32);
const double crLHS38 = r_DN(1,1)*r_N[0];
const double crLHS39 = crLHS23*crLHS30;
const double crLHS40 = -crLHS39*r_N[1];
const double crLHS41 = crLHS8*(crLHS24*crLHS36 + crLHS38 + crLHS40);
const double crLHS42 = r_DN(0,0)*r_DN(2,0);
const double crLHS43 = r_DN(0,1)*r_DN(2,1);
const double crLHS44 = crLHS13*(crLHS42 + crLHS43);
const double crLHS45 = r_DN(2,0)*r_N[0];
const double crLHS46 = -crLHS31*r_N[2];
const double crLHS47 = bdf0*r_N[2];
const double crLHS48 = crLHS10*r_DN(2,1) + crLHS9*r_DN(2,0);
const double crLHS49 = crLHS47 + crLHS48;
const double crLHS50 = crLHS12*crLHS49;
const double crLHS51 = crLHS8*(crLHS17*crLHS50 + crLHS45 + crLHS46);
const double crLHS52 = r_DN(2,1)*r_N[0];
const double crLHS53 = -crLHS39*r_N[2];
const double crLHS54 = crLHS8*(crLHS24*crLHS50 + crLHS52 + crLHS53);
const double crLHS55 = r_DN(0,0)*r_DN(3,0);
const double crLHS56 = r_DN(0,1)*r_DN(3,1);
const double crLHS57 = crLHS13*(crLHS55 + crLHS56);
const double crLHS58 = r_DN(3,0)*r_N[0];
const double crLHS59 = -crLHS31*r_N[3];
const double crLHS60 = crLHS10*r_DN(3,1) + crLHS9*r_DN(3,0);
const double crLHS61 = bdf0*r_N[3] + crLHS60;
const double crLHS62 = crLHS12*crLHS61;
const double crLHS63 = crLHS8*(crLHS17*crLHS62 + crLHS58 + crLHS59);
const double crLHS64 = r_DN(3,1)*r_N[0];
const double crLHS65 = -crLHS39*r_N[3];
const double crLHS66 = crLHS8*(crLHS24*crLHS62 + crLHS64 + crLHS65);
const double crLHS67 = crLHS30*crLHS6;
const double crLHS68 = r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1);
const double crLHS69 = crLHS12*crLHS68;
const double crLHS70 = crLHS12*crLHS7;
const double crLHS71 = crLHS19*crLHS70;
const double crLHS72 = 1.0/(crLHS3*crLHS3);
const double crLHS73 = crLHS72*r_N[0];
const double crLHS74 = crLHS10*crLHS23 + crLHS15*crLHS9;
const double crLHS75 = crLHS12*crLHS6*crLHS74;
const double crLHS76 = crLHS2*(-crLHS67*crLHS69 - crLHS71 + crLHS73*crLHS75 + r_N[0]);
const double crLHS77 = -crLHS76*r_DN(0,0);
const double crLHS78 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS79 = r_C(0,2)*r_DN(0,0);
const double crLHS80 = crLHS79 + r_C(2,2)*r_DN(0,1);
const double crLHS81 = -crLHS31 + r_DN(0,0);
const double crLHS82 = crLHS11*h*1.0/stab_c1 + mu*1.0/rho_lin;
const double crLHS83 = crLHS17*crLHS82;
const double crLHS84 = bdf0*crLHS6;
const double crLHS85 = 1.0/(crLHS5*crLHS5)*1.0/(c_p*c_p)*(gamma*gamma)*(p_th*p_th);
const double crLHS86 = crLHS72*crLHS85;
const double crLHS87 = crLHS19*crLHS86;
const double crLHS88 = crLHS21*crLHS85;
const double crLHS89 = crLHS68*crLHS73;
const double crLHS90 = crLHS74*1.0/(crLHS3*crLHS3*crLHS3);
const double crLHS91 = crLHS90*r_N[0];
const double crLHS92 = crLHS16*crLHS84 + crLHS19*crLHS67 + crLHS21*crLHS87 + crLHS88*crLHS89 - crLHS88*crLHS91;
const double crLHS93 = crLHS79 + r_C(0,1)*r_DN(0,1);
const double crLHS94 = r_C(1,2)*r_DN(0,1);
const double crLHS95 = crLHS94 + r_C(2,2)*r_DN(0,0);
const double crLHS96 = crLHS17*r_DN(0,1);
const double crLHS97 = -crLHS39 + r_DN(0,1);
const double crLHS98 = r_DN(0,0)*r_N[1];
const double crLHS99 = crLHS68*crLHS70;
const double crLHS100 = crLHS7*r_DN(1,0);
const double crLHS101 = crLHS12*crLHS19;
const double crLHS102 = crLHS72*crLHS75;
const double crLHS103 = -crLHS2*(-crLHS100*crLHS101 + crLHS102*crLHS29 - crLHS29*crLHS99 + crLHS98);
const double crLHS104 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS105 = r_C(0,2)*r_DN(1,0);
const double crLHS106 = crLHS105 + r_C(2,2)*r_DN(1,1);
const double crLHS107 = crLHS4*r_N[1];
const double crLHS108 = crLHS107*crLHS15;
const double crLHS109 = -crLHS108 + r_DN(1,0);
const double crLHS110 = crLHS107*crLHS6;
const double crLHS111 = crLHS110*crLHS18;
const double crLHS112 = crLHS111 + crLHS26*crLHS7;
const double crLHS113 = crLHS36*crLHS85;
const double crLHS114 = crLHS113*crLHS89 - crLHS113*crLHS91 + crLHS34*crLHS67 + crLHS36*crLHS87;
const double crLHS115 = crLHS105 + r_C(0,1)*r_DN(1,1);
const double crLHS116 = r_C(1,2)*r_DN(1,1);
const double crLHS117 = crLHS116 + r_C(2,2)*r_DN(1,0);
const double crLHS118 = crLHS17*r_DN(1,1);
const double crLHS119 = crLHS107*crLHS23;
const double crLHS120 = -crLHS119 + r_DN(1,1);
const double crLHS121 = r_DN(0,0)*r_N[2];
const double crLHS122 = crLHS7*r_DN(2,0);
const double crLHS123 = -crLHS2*(-crLHS101*crLHS122 + crLHS102*crLHS45 + crLHS121 - crLHS45*crLHS99);
const double crLHS124 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS125 = r_C(0,2)*r_DN(2,0);
const double crLHS126 = crLHS125 + r_C(2,2)*r_DN(2,1);
const double crLHS127 = crLHS4*r_N[2];
const double crLHS128 = crLHS127*crLHS15;
const double crLHS129 = -crLHS128 + r_DN(2,0);
const double crLHS130 = crLHS127*crLHS6;
const double crLHS131 = crLHS130*crLHS18;
const double crLHS132 = crLHS131 + crLHS42*crLHS7;
const double crLHS133 = crLHS50*crLHS85;
const double crLHS134 = crLHS133*crLHS89 - crLHS133*crLHS91 + crLHS48*crLHS67 + crLHS50*crLHS87;
const double crLHS135 = crLHS125 + r_C(0,1)*r_DN(2,1);
const double crLHS136 = r_C(1,2)*r_DN(2,1);
const double crLHS137 = crLHS136 + r_C(2,2)*r_DN(2,0);
const double crLHS138 = crLHS17*r_DN(2,1);
const double crLHS139 = crLHS127*crLHS23;
const double crLHS140 = -crLHS139 + r_DN(2,1);
const double crLHS141 = r_DN(0,0)*r_N[3];
const double crLHS142 = crLHS7*r_DN(3,0);
const double crLHS143 = -crLHS2*(-crLHS101*crLHS142 + crLHS102*crLHS58 + crLHS141 - crLHS58*crLHS99);
const double crLHS144 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS145 = r_C(0,2)*r_DN(3,0);
const double crLHS146 = crLHS145 + r_C(2,2)*r_DN(3,1);
const double crLHS147 = crLHS4*r_N[3];
const double crLHS148 = -crLHS147*crLHS15 + r_DN(3,0);
const double crLHS149 = crLHS147*crLHS6;
const double crLHS150 = crLHS149*crLHS18;
const double crLHS151 = crLHS150 + crLHS55*crLHS7;
const double crLHS152 = crLHS62*crLHS85;
const double crLHS153 = crLHS152*crLHS89 - crLHS152*crLHS91 + crLHS60*crLHS67 + crLHS62*crLHS87;
const double crLHS154 = crLHS145 + r_C(0,1)*r_DN(3,1);
const double crLHS155 = r_C(1,2)*r_DN(3,1);
const double crLHS156 = crLHS155 + r_C(2,2)*r_DN(3,0);
const double crLHS157 = crLHS17*r_DN(3,1);
const double crLHS158 = -crLHS147*crLHS23 + r_DN(3,1);
const double crLHS159 = -crLHS76*r_DN(0,1);
const double crLHS160 = crLHS94 + r_C(0,1)*r_DN(0,0);
const double crLHS161 = crLHS24*crLHS82;
const double crLHS162 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS163 = r_DN(0,1)*r_N[1];
const double crLHS164 = crLHS7*r_DN(1,1);
const double crLHS165 = -crLHS2*(-crLHS101*crLHS164 + crLHS102*crLHS38 + crLHS163 - crLHS38*crLHS99);
const double crLHS166 = crLHS116 + r_C(0,1)*r_DN(1,0);
const double crLHS167 = crLHS24*r_DN(1,0);
const double crLHS168 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS169 = crLHS111 + crLHS27*crLHS7;
const double crLHS170 = r_DN(0,1)*r_N[2];
const double crLHS171 = -crLHS2*(crLHS102*crLHS52 + crLHS170 - crLHS52*crLHS99 - crLHS71*r_DN(2,1));
const double crLHS172 = crLHS136 + r_C(0,1)*r_DN(2,0);
const double crLHS173 = crLHS24*r_DN(2,0);
const double crLHS174 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS175 = crLHS131 + crLHS43*crLHS7;
const double crLHS176 = r_DN(0,1)*r_N[3];
const double crLHS177 = -crLHS2*(crLHS102*crLHS64 + crLHS176 - crLHS64*crLHS99 - crLHS71*r_DN(3,1));
const double crLHS178 = crLHS155 + r_C(0,1)*r_DN(3,0);
const double crLHS179 = crLHS24*r_DN(3,0);
const double crLHS180 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS181 = crLHS150 + crLHS56*crLHS7;
const double crLHS182 = crLHS8*(crLHS100*crLHS21 + crLHS32 + crLHS98);
const double crLHS183 = crLHS8*(crLHS163 + crLHS164*crLHS21 + crLHS40);
const double crLHS184 = r_DN(1,0)*r_DN(1,0);
const double crLHS185 = r_DN(1,1)*r_DN(1,1);
const double crLHS186 = crLHS13*(crLHS184 + crLHS185);
const double crLHS187 = crLHS4*(r_N[1]*r_N[1]);
const double crLHS188 = crLHS8*(crLHS100*crLHS36 - crLHS15*crLHS187 + r_DN(1,0)*r_N[1]);
const double crLHS189 = crLHS8*(crLHS164*crLHS36 - crLHS187*crLHS23 + r_DN(1,1)*r_N[1]);
const double crLHS190 = r_DN(1,0)*r_DN(2,0);
const double crLHS191 = r_DN(1,1)*r_DN(2,1);
const double crLHS192 = crLHS13*(crLHS190 + crLHS191);
const double crLHS193 = r_DN(2,0)*r_N[1];
const double crLHS194 = -crLHS108*r_N[2];
const double crLHS195 = crLHS8*(crLHS100*crLHS50 + crLHS193 + crLHS194);
const double crLHS196 = r_DN(2,1)*r_N[1];
const double crLHS197 = -crLHS119*r_N[2];
const double crLHS198 = crLHS8*(crLHS164*crLHS50 + crLHS196 + crLHS197);
const double crLHS199 = r_DN(1,0)*r_DN(3,0);
const double crLHS200 = r_DN(1,1)*r_DN(3,1);
const double crLHS201 = crLHS13*(crLHS199 + crLHS200);
const double crLHS202 = r_DN(3,0)*r_N[1];
const double crLHS203 = -crLHS108*r_N[3];
const double crLHS204 = crLHS8*(crLHS100*crLHS62 + crLHS202 + crLHS203);
const double crLHS205 = r_DN(3,1)*r_N[1];
const double crLHS206 = -crLHS119*r_N[3];
const double crLHS207 = crLHS8*(crLHS164*crLHS62 + crLHS205 + crLHS206);
const double crLHS208 = crLHS12*crLHS34;
const double crLHS209 = crLHS2*(-crLHS102*crLHS98 + crLHS17*crLHS208 - crLHS29 + crLHS98*crLHS99);
const double crLHS210 = crLHS100*crLHS82;
const double crLHS211 = crLHS34*crLHS86;
const double crLHS212 = crLHS68*crLHS86;
const double crLHS213 = crLHS212*r_N[1];
const double crLHS214 = crLHS90*r_N[1];
const double crLHS215 = crLHS110*crLHS19 + crLHS21*crLHS211 + crLHS21*crLHS213 - crLHS214*crLHS88;
const double crLHS216 = crLHS34*crLHS70;
const double crLHS217 = crLHS2*(crLHS102*r_N[1] - crLHS110*crLHS69 - crLHS216 + r_N[1]);
const double crLHS218 = -crLHS217*r_DN(1,0);
const double crLHS219 = crLHS110*crLHS34 - crLHS113*crLHS214 + crLHS187*crLHS84 + crLHS211*crLHS36 + crLHS213*crLHS36;
const double crLHS220 = crLHS100*r_DN(1,1);
const double crLHS221 = r_DN(1,0)*r_N[2];
const double crLHS222 = -crLHS2*(crLHS102*crLHS193 - crLHS122*crLHS208 - crLHS193*crLHS99 + crLHS221);
const double crLHS223 = crLHS130*crLHS33;
const double crLHS224 = crLHS190*crLHS7 + crLHS223;
const double crLHS225 = crLHS110*crLHS48 - crLHS133*crLHS214 + crLHS211*crLHS50 + crLHS213*crLHS50;
const double crLHS226 = crLHS100*r_DN(2,1);
const double crLHS227 = r_DN(1,0)*r_N[3];
const double crLHS228 = -crLHS2*(crLHS102*crLHS202 - crLHS142*crLHS208 - crLHS202*crLHS99 + crLHS227);
const double crLHS229 = crLHS149*crLHS33;
const double crLHS230 = crLHS199*crLHS7 + crLHS229;
const double crLHS231 = crLHS110*crLHS60 - crLHS152*crLHS214 + crLHS211*crLHS62 + crLHS213*crLHS62;
const double crLHS232 = crLHS100*r_DN(3,1);
const double crLHS233 = crLHS2*(-crLHS102*crLHS163 + crLHS163*crLHS99 + crLHS208*crLHS24 - crLHS38);
const double crLHS234 = crLHS164*crLHS82;
const double crLHS235 = -crLHS217*r_DN(1,1);
const double crLHS236 = r_DN(1,1)*r_N[2];
const double crLHS237 = -crLHS2*(crLHS102*crLHS196 - crLHS196*crLHS99 - crLHS216*r_DN(2,1) + crLHS236);
const double crLHS238 = crLHS164*r_DN(2,0);
const double crLHS239 = crLHS191*crLHS7 + crLHS223;
const double crLHS240 = r_DN(1,1)*r_N[3];
const double crLHS241 = -crLHS2*(crLHS102*crLHS205 - crLHS205*crLHS99 - crLHS216*r_DN(3,1) + crLHS240);
const double crLHS242 = crLHS164*r_DN(3,0);
const double crLHS243 = crLHS200*crLHS7 + crLHS229;
const double crLHS244 = crLHS8*(crLHS121 + crLHS122*crLHS21 + crLHS46);
const double crLHS245 = crLHS70*r_DN(2,1);
const double crLHS246 = crLHS8*(crLHS170 + crLHS20*crLHS245 + crLHS53);
const double crLHS247 = crLHS8*(crLHS122*crLHS36 + crLHS194 + crLHS221);
const double crLHS248 = crLHS8*(crLHS197 + crLHS236 + crLHS245*crLHS35);
const double crLHS249 = r_DN(2,0)*r_DN(2,0);
const double crLHS250 = r_DN(2,1)*r_DN(2,1);
const double crLHS251 = crLHS13*(crLHS249 + crLHS250);
const double crLHS252 = crLHS4*(r_N[2]*r_N[2]);
const double crLHS253 = crLHS8*(crLHS122*crLHS50 - crLHS15*crLHS252 + r_DN(2,0)*r_N[2]);
const double crLHS254 = crLHS8*(-crLHS23*crLHS252 + crLHS245*crLHS49 + r_DN(2,1)*r_N[2]);
const double crLHS255 = r_DN(2,0)*r_DN(3,0);
const double crLHS256 = r_DN(2,1)*r_DN(3,1);
const double crLHS257 = crLHS13*(crLHS255 + crLHS256);
const double crLHS258 = r_DN(3,0)*r_N[2];
const double crLHS259 = -crLHS128*r_N[3];
const double crLHS260 = crLHS8*(crLHS122*crLHS62 + crLHS258 + crLHS259);
const double crLHS261 = r_DN(3,1)*r_N[2];
const double crLHS262 = -crLHS139*r_N[3];
const double crLHS263 = crLHS8*(crLHS245*crLHS61 + crLHS261 + crLHS262);
const double crLHS264 = crLHS12*crLHS48;
const double crLHS265 = crLHS2*(-crLHS102*crLHS121 + crLHS121*crLHS99 + crLHS17*crLHS264 - crLHS45);
const double crLHS266 = crLHS122*crLHS82;
const double crLHS267 = crLHS48*crLHS86;
const double crLHS268 = crLHS212*r_N[2];
const double crLHS269 = crLHS90*r_N[2];
const double crLHS270 = crLHS130*crLHS19 + crLHS21*crLHS267 + crLHS21*crLHS268 - crLHS269*crLHS88;
const double crLHS271 = crLHS2*(crLHS100*crLHS264 - crLHS102*crLHS221 - crLHS193 + crLHS221*crLHS99);
const double crLHS272 = -crLHS113*crLHS269 + crLHS130*crLHS34 + crLHS267*crLHS36 + crLHS268*crLHS36;
const double crLHS273 = crLHS48*crLHS70;
const double crLHS274 = crLHS2*(crLHS102*r_N[2] - crLHS130*crLHS69 - crLHS273 + r_N[2]);
const double crLHS275 = -crLHS274*r_DN(2,0);
const double crLHS276 = crLHS130*crLHS48 - crLHS133*crLHS269 + crLHS252*crLHS84 + crLHS267*crLHS50 + crLHS268*crLHS50;
const double crLHS277 = crLHS122*r_DN(2,1);
const double crLHS278 = r_DN(2,0)*r_N[3];
const double crLHS279 = -crLHS2*(crLHS102*crLHS258 - crLHS142*crLHS264 - crLHS258*crLHS99 + crLHS278);
const double crLHS280 = crLHS149*crLHS47;
const double crLHS281 = crLHS255*crLHS7 + crLHS280;
const double crLHS282 = crLHS130*crLHS60 - crLHS152*crLHS269 + crLHS267*crLHS62 + crLHS268*crLHS62;
const double crLHS283 = crLHS122*r_DN(3,1);
const double crLHS284 = crLHS2*(-crLHS102*crLHS170 + crLHS170*crLHS99 + crLHS24*crLHS264 - crLHS52);
const double crLHS285 = crLHS81*crLHS82;
const double crLHS286 = crLHS7*r_DN(2,1);
const double crLHS287 = crLHS286*crLHS82;
const double crLHS288 = crLHS2*(-crLHS102*crLHS236 + crLHS164*crLHS264 - crLHS196 + crLHS236*crLHS99);
const double crLHS289 = -crLHS274*r_DN(2,1);
const double crLHS290 = r_DN(2,1)*r_N[3];
const double crLHS291 = -crLHS2*(crLHS102*crLHS261 - crLHS261*crLHS99 - crLHS273*r_DN(3,1) + crLHS290);
const double crLHS292 = crLHS142*r_DN(2,1);
const double crLHS293 = crLHS256*crLHS7 + crLHS280;
const double crLHS294 = crLHS8*(crLHS141 + crLHS142*crLHS21 + crLHS59);
const double crLHS295 = crLHS70*r_DN(3,1);
const double crLHS296 = crLHS8*(crLHS176 + crLHS20*crLHS295 + crLHS65);
const double crLHS297 = crLHS8*(crLHS142*crLHS36 + crLHS203 + crLHS227);
const double crLHS298 = crLHS8*(crLHS206 + crLHS240 + crLHS295*crLHS35);
const double crLHS299 = crLHS8*(crLHS142*crLHS50 + crLHS259 + crLHS278);
const double crLHS300 = crLHS8*(crLHS262 + crLHS290 + crLHS295*crLHS49);
const double crLHS301 = r_DN(3,0)*r_DN(3,0);
const double crLHS302 = r_DN(3,1)*r_DN(3,1);
const double crLHS303 = crLHS13*(crLHS301 + crLHS302);
const double crLHS304 = crLHS4*(r_N[3]*r_N[3]);
const double crLHS305 = crLHS8*(crLHS142*crLHS62 - crLHS15*crLHS304 + r_DN(3,0)*r_N[3]);
const double crLHS306 = crLHS8*(-crLHS23*crLHS304 + crLHS295*crLHS61 + r_DN(3,1)*r_N[3]);
const double crLHS307 = crLHS12*crLHS60;
const double crLHS308 = crLHS2*(-crLHS102*crLHS141 + crLHS141*crLHS99 + crLHS17*crLHS307 - crLHS58);
const double crLHS309 = crLHS60*crLHS86;
const double crLHS310 = crLHS212*r_N[3];
const double crLHS311 = crLHS90*r_N[3];
const double crLHS312 = crLHS149*crLHS19 + crLHS21*crLHS309 + crLHS21*crLHS310 - crLHS311*crLHS88;
const double crLHS313 = crLHS142*crLHS82;
const double crLHS314 = crLHS2*(crLHS100*crLHS307 - crLHS102*crLHS227 - crLHS202 + crLHS227*crLHS99);
const double crLHS315 = -crLHS113*crLHS311 + crLHS149*crLHS34 + crLHS309*crLHS36 + crLHS310*crLHS36;
const double crLHS316 = crLHS2*(-crLHS102*crLHS278 + crLHS122*crLHS307 - crLHS258 + crLHS278*crLHS99);
const double crLHS317 = -crLHS133*crLHS311 + crLHS149*crLHS48 + crLHS309*crLHS50 + crLHS310*crLHS50;
const double crLHS318 = crLHS60*crLHS70;
const double crLHS319 = crLHS2*(crLHS102*r_N[3] - crLHS149*crLHS69 - crLHS318 + r_N[3]);
const double crLHS320 = -crLHS319*r_DN(3,0);
const double crLHS321 = crLHS149*crLHS60 - crLHS152*crLHS311 + crLHS304*crLHS84 + crLHS309*crLHS62 + crLHS310*crLHS62;
const double crLHS322 = crLHS142*r_DN(3,1);
const double crLHS323 = crLHS2*(-crLHS102*crLHS176 + crLHS176*crLHS99 + crLHS24*crLHS307 - crLHS64);
const double crLHS324 = crLHS7*r_DN(3,1);
const double crLHS325 = crLHS324*crLHS82;
const double crLHS326 = crLHS2*(-crLHS102*crLHS240 + crLHS164*crLHS307 - crLHS205 + crLHS240*crLHS99);
const double crLHS327 = crLHS2*(-crLHS102*crLHS290 - crLHS261 + crLHS290*crLHS99 + crLHS318*r_DN(2,1));
const double crLHS328 = -crLHS319*r_DN(3,1);
rLHS(0,0)+=crLHS14;
rLHS(0,1)+=crLHS22;
rLHS(0,2)+=crLHS25;
rLHS(0,3)+=crLHS14;
rLHS(0,4)+=crLHS28;
rLHS(0,5)+=crLHS37;
rLHS(0,6)+=crLHS41;
rLHS(0,7)+=crLHS28;
rLHS(0,8)+=crLHS44;
rLHS(0,9)+=crLHS51;
rLHS(0,10)+=crLHS54;
rLHS(0,11)+=crLHS44;
rLHS(0,12)+=crLHS57;
rLHS(0,13)+=crLHS63;
rLHS(0,14)+=crLHS66;
rLHS(0,15)+=crLHS57;
rLHS(1,0)+=crLHS77;
rLHS(1,1)+=crLHS2*(crLHS0*crLHS7 + crLHS78*r_DN(0,0) + crLHS80*r_DN(0,1) + crLHS81*crLHS83 + crLHS92);
rLHS(1,2)+=crLHS2*(crLHS83*crLHS97 + crLHS93*r_DN(0,0) + crLHS95*r_DN(0,1) + crLHS96);
rLHS(1,3)+=crLHS77;
rLHS(1,4)+=crLHS103;
rLHS(1,5)+=crLHS2*(crLHS104*r_DN(0,0) + crLHS106*r_DN(0,1) + crLHS109*crLHS83 + crLHS112 + crLHS114);
rLHS(1,6)+=crLHS2*(crLHS115*r_DN(0,0) + crLHS117*r_DN(0,1) + crLHS118 + crLHS120*crLHS83);
rLHS(1,7)+=crLHS103;
rLHS(1,8)+=crLHS123;
rLHS(1,9)+=crLHS2*(crLHS124*r_DN(0,0) + crLHS126*r_DN(0,1) + crLHS129*crLHS83 + crLHS132 + crLHS134);
rLHS(1,10)+=crLHS2*(crLHS135*r_DN(0,0) + crLHS137*r_DN(0,1) + crLHS138 + crLHS140*crLHS83);
rLHS(1,11)+=crLHS123;
rLHS(1,12)+=crLHS143;
rLHS(1,13)+=crLHS2*(crLHS144*r_DN(0,0) + crLHS146*r_DN(0,1) + crLHS148*crLHS83 + crLHS151 + crLHS153);
rLHS(1,14)+=crLHS2*(crLHS154*r_DN(0,0) + crLHS156*r_DN(0,1) + crLHS157 + crLHS158*crLHS83);
rLHS(1,15)+=crLHS143;
rLHS(2,0)+=crLHS159;
rLHS(2,1)+=crLHS2*(crLHS160*r_DN(0,1) + crLHS161*crLHS81 + crLHS80*r_DN(0,0) + crLHS96);
rLHS(2,2)+=crLHS2*(crLHS1*crLHS7 + crLHS161*crLHS97 + crLHS162*r_DN(0,1) + crLHS92 + crLHS95*r_DN(0,0));
rLHS(2,3)+=crLHS159;
rLHS(2,4)+=crLHS165;
rLHS(2,5)+=crLHS2*(crLHS106*r_DN(0,0) + crLHS109*crLHS161 + crLHS166*r_DN(0,1) + crLHS167);
rLHS(2,6)+=crLHS2*(crLHS114 + crLHS117*r_DN(0,0) + crLHS120*crLHS161 + crLHS168*r_DN(0,1) + crLHS169);
rLHS(2,7)+=crLHS165;
rLHS(2,8)+=crLHS171;
rLHS(2,9)+=crLHS2*(crLHS126*r_DN(0,0) + crLHS129*crLHS161 + crLHS172*r_DN(0,1) + crLHS173);
rLHS(2,10)+=crLHS2*(crLHS134 + crLHS137*r_DN(0,0) + crLHS140*crLHS161 + crLHS174*r_DN(0,1) + crLHS175);
rLHS(2,11)+=crLHS171;
rLHS(2,12)+=crLHS177;
rLHS(2,13)+=crLHS2*(crLHS146*r_DN(0,0) + crLHS148*crLHS161 + crLHS178*r_DN(0,1) + crLHS179);
rLHS(2,14)+=crLHS2*(crLHS153 + crLHS156*r_DN(0,0) + crLHS158*crLHS161 + crLHS180*r_DN(0,1) + crLHS181);
rLHS(2,15)+=crLHS177;
rLHS(3,0)+=crLHS14;
rLHS(3,1)+=crLHS22;
rLHS(3,2)+=crLHS25;
rLHS(3,3)+=crLHS14;
rLHS(3,4)+=crLHS28;
rLHS(3,5)+=crLHS37;
rLHS(3,6)+=crLHS41;
rLHS(3,7)+=crLHS28;
rLHS(3,8)+=crLHS44;
rLHS(3,9)+=crLHS51;
rLHS(3,10)+=crLHS54;
rLHS(3,11)+=crLHS44;
rLHS(3,12)+=crLHS57;
rLHS(3,13)+=crLHS63;
rLHS(3,14)+=crLHS66;
rLHS(3,15)+=crLHS57;
rLHS(4,0)+=crLHS28;
rLHS(4,1)+=crLHS182;
rLHS(4,2)+=crLHS183;
rLHS(4,3)+=crLHS28;
rLHS(4,4)+=crLHS186;
rLHS(4,5)+=crLHS188;
rLHS(4,6)+=crLHS189;
rLHS(4,7)+=crLHS186;
rLHS(4,8)+=crLHS192;
rLHS(4,9)+=crLHS195;
rLHS(4,10)+=crLHS198;
rLHS(4,11)+=crLHS192;
rLHS(4,12)+=crLHS201;
rLHS(4,13)+=crLHS204;
rLHS(4,14)+=crLHS207;
rLHS(4,15)+=crLHS201;
rLHS(5,0)+=crLHS209;
rLHS(5,1)+=crLHS2*(crLHS112 + crLHS210*crLHS81 + crLHS215 + crLHS78*r_DN(1,0) + crLHS80*r_DN(1,1));
rLHS(5,2)+=crLHS2*(crLHS167 + crLHS210*crLHS97 + crLHS93*r_DN(1,0) + crLHS95*r_DN(1,1));
rLHS(5,3)+=crLHS209;
rLHS(5,4)+=crLHS218;
rLHS(5,5)+=crLHS2*(crLHS104*r_DN(1,0) + crLHS106*r_DN(1,1) + crLHS109*crLHS210 + crLHS184*crLHS7 + crLHS219);
rLHS(5,6)+=crLHS2*(crLHS115*r_DN(1,0) + crLHS117*r_DN(1,1) + crLHS120*crLHS210 + crLHS220);
rLHS(5,7)+=crLHS218;
rLHS(5,8)+=crLHS222;
rLHS(5,9)+=crLHS2*(crLHS124*r_DN(1,0) + crLHS126*r_DN(1,1) + crLHS129*crLHS210 + crLHS224 + crLHS225);
rLHS(5,10)+=crLHS2*(crLHS135*r_DN(1,0) + crLHS137*r_DN(1,1) + crLHS140*crLHS210 + crLHS226);
rLHS(5,11)+=crLHS222;
rLHS(5,12)+=crLHS228;
rLHS(5,13)+=crLHS2*(crLHS144*r_DN(1,0) + crLHS146*r_DN(1,1) + crLHS148*crLHS210 + crLHS230 + crLHS231);
rLHS(5,14)+=crLHS2*(crLHS154*r_DN(1,0) + crLHS156*r_DN(1,1) + crLHS158*crLHS210 + crLHS232);
rLHS(5,15)+=crLHS228;
rLHS(6,0)+=crLHS233;
rLHS(6,1)+=crLHS2*(crLHS118 + crLHS160*r_DN(1,1) + crLHS234*crLHS81 + crLHS80*r_DN(1,0));
rLHS(6,2)+=crLHS2*(crLHS162*r_DN(1,1) + crLHS169 + crLHS215 + crLHS234*crLHS97 + crLHS95*r_DN(1,0));
rLHS(6,3)+=crLHS233;
rLHS(6,4)+=crLHS235;
rLHS(6,5)+=crLHS2*(crLHS106*r_DN(1,0) + crLHS109*crLHS234 + crLHS166*r_DN(1,1) + crLHS220);
rLHS(6,6)+=crLHS2*(crLHS117*r_DN(1,0) + crLHS120*crLHS234 + crLHS168*r_DN(1,1) + crLHS185*crLHS7 + crLHS219);
rLHS(6,7)+=crLHS235;
rLHS(6,8)+=crLHS237;
rLHS(6,9)+=crLHS2*(crLHS126*r_DN(1,0) + crLHS129*crLHS234 + crLHS172*r_DN(1,1) + crLHS238);
rLHS(6,10)+=crLHS2*(crLHS137*r_DN(1,0) + crLHS140*crLHS234 + crLHS174*r_DN(1,1) + crLHS225 + crLHS239);
rLHS(6,11)+=crLHS237;
rLHS(6,12)+=crLHS241;
rLHS(6,13)+=crLHS2*(crLHS146*r_DN(1,0) + crLHS148*crLHS234 + crLHS178*r_DN(1,1) + crLHS242);
rLHS(6,14)+=crLHS2*(crLHS156*r_DN(1,0) + crLHS158*crLHS234 + crLHS180*r_DN(1,1) + crLHS231 + crLHS243);
rLHS(6,15)+=crLHS241;
rLHS(7,0)+=crLHS28;
rLHS(7,1)+=crLHS182;
rLHS(7,2)+=crLHS183;
rLHS(7,3)+=crLHS28;
rLHS(7,4)+=crLHS186;
rLHS(7,5)+=crLHS188;
rLHS(7,6)+=crLHS189;
rLHS(7,7)+=crLHS186;
rLHS(7,8)+=crLHS192;
rLHS(7,9)+=crLHS195;
rLHS(7,10)+=crLHS198;
rLHS(7,11)+=crLHS192;
rLHS(7,12)+=crLHS201;
rLHS(7,13)+=crLHS204;
rLHS(7,14)+=crLHS207;
rLHS(7,15)+=crLHS201;
rLHS(8,0)+=crLHS44;
rLHS(8,1)+=crLHS244;
rLHS(8,2)+=crLHS246;
rLHS(8,3)+=crLHS44;
rLHS(8,4)+=crLHS192;
rLHS(8,5)+=crLHS247;
rLHS(8,6)+=crLHS248;
rLHS(8,7)+=crLHS192;
rLHS(8,8)+=crLHS251;
rLHS(8,9)+=crLHS253;
rLHS(8,10)+=crLHS254;
rLHS(8,11)+=crLHS251;
rLHS(8,12)+=crLHS257;
rLHS(8,13)+=crLHS260;
rLHS(8,14)+=crLHS263;
rLHS(8,15)+=crLHS257;
rLHS(9,0)+=crLHS265;
rLHS(9,1)+=crLHS2*(crLHS132 + crLHS266*crLHS81 + crLHS270 + crLHS78*r_DN(2,0) + crLHS80*r_DN(2,1));
rLHS(9,2)+=crLHS2*(crLHS173 + crLHS266*crLHS97 + crLHS93*r_DN(2,0) + crLHS95*r_DN(2,1));
rLHS(9,3)+=crLHS265;
rLHS(9,4)+=crLHS271;
rLHS(9,5)+=crLHS2*(crLHS104*r_DN(2,0) + crLHS106*r_DN(2,1) + crLHS109*crLHS266 + crLHS224 + crLHS272);
rLHS(9,6)+=crLHS2*(crLHS115*r_DN(2,0) + crLHS117*r_DN(2,1) + crLHS120*crLHS266 + crLHS238);
rLHS(9,7)+=crLHS271;
rLHS(9,8)+=crLHS275;
rLHS(9,9)+=crLHS2*(crLHS124*r_DN(2,0) + crLHS126*r_DN(2,1) + crLHS129*crLHS266 + crLHS249*crLHS7 + crLHS276);
rLHS(9,10)+=crLHS2*(crLHS135*r_DN(2,0) + crLHS137*r_DN(2,1) + crLHS140*crLHS266 + crLHS277);
rLHS(9,11)+=crLHS275;
rLHS(9,12)+=crLHS279;
rLHS(9,13)+=crLHS2*(crLHS144*r_DN(2,0) + crLHS146*r_DN(2,1) + crLHS148*crLHS266 + crLHS281 + crLHS282);
rLHS(9,14)+=crLHS2*(crLHS154*r_DN(2,0) + crLHS156*r_DN(2,1) + crLHS158*crLHS266 + crLHS283);
rLHS(9,15)+=crLHS279;
rLHS(10,0)+=crLHS284;
rLHS(10,1)+=crLHS2*(crLHS138 + crLHS160*r_DN(2,1) + crLHS285*crLHS286 + crLHS80*r_DN(2,0));
rLHS(10,2)+=crLHS2*(crLHS162*r_DN(2,1) + crLHS175 + crLHS270 + crLHS287*crLHS97 + crLHS95*r_DN(2,0));
rLHS(10,3)+=crLHS284;
rLHS(10,4)+=crLHS288;
rLHS(10,5)+=crLHS2*(crLHS106*r_DN(2,0) + crLHS109*crLHS287 + crLHS166*r_DN(2,1) + crLHS226);
rLHS(10,6)+=crLHS2*(crLHS117*r_DN(2,0) + crLHS120*crLHS287 + crLHS168*r_DN(2,1) + crLHS239 + crLHS272);
rLHS(10,7)+=crLHS288;
rLHS(10,8)+=crLHS289;
rLHS(10,9)+=crLHS2*(crLHS126*r_DN(2,0) + crLHS129*crLHS287 + crLHS172*r_DN(2,1) + crLHS277);
rLHS(10,10)+=crLHS2*(crLHS137*r_DN(2,0) + crLHS140*crLHS287 + crLHS174*r_DN(2,1) + crLHS250*crLHS7 + crLHS276);
rLHS(10,11)+=crLHS289;
rLHS(10,12)+=crLHS291;
rLHS(10,13)+=crLHS2*(crLHS146*r_DN(2,0) + crLHS148*crLHS287 + crLHS178*r_DN(2,1) + crLHS292);
rLHS(10,14)+=crLHS2*(crLHS156*r_DN(2,0) + crLHS158*crLHS287 + crLHS180*r_DN(2,1) + crLHS282 + crLHS293);
rLHS(10,15)+=crLHS291;
rLHS(11,0)+=crLHS44;
rLHS(11,1)+=crLHS244;
rLHS(11,2)+=crLHS246;
rLHS(11,3)+=crLHS44;
rLHS(11,4)+=crLHS192;
rLHS(11,5)+=crLHS247;
rLHS(11,6)+=crLHS248;
rLHS(11,7)+=crLHS192;
rLHS(11,8)+=crLHS251;
rLHS(11,9)+=crLHS253;
rLHS(11,10)+=crLHS254;
rLHS(11,11)+=crLHS251;
rLHS(11,12)+=crLHS257;
rLHS(11,13)+=crLHS260;
rLHS(11,14)+=crLHS263;
rLHS(11,15)+=crLHS257;
rLHS(12,0)+=crLHS57;
rLHS(12,1)+=crLHS294;
rLHS(12,2)+=crLHS296;
rLHS(12,3)+=crLHS57;
rLHS(12,4)+=crLHS201;
rLHS(12,5)+=crLHS297;
rLHS(12,6)+=crLHS298;
rLHS(12,7)+=crLHS201;
rLHS(12,8)+=crLHS257;
rLHS(12,9)+=crLHS299;
rLHS(12,10)+=crLHS300;
rLHS(12,11)+=crLHS257;
rLHS(12,12)+=crLHS303;
rLHS(12,13)+=crLHS305;
rLHS(12,14)+=crLHS306;
rLHS(12,15)+=crLHS303;
rLHS(13,0)+=crLHS308;
rLHS(13,1)+=crLHS2*(crLHS142*crLHS285 + crLHS151 + crLHS312 + crLHS78*r_DN(3,0) + crLHS80*r_DN(3,1));
rLHS(13,2)+=crLHS2*(crLHS179 + crLHS313*crLHS97 + crLHS93*r_DN(3,0) + crLHS95*r_DN(3,1));
rLHS(13,3)+=crLHS308;
rLHS(13,4)+=crLHS314;
rLHS(13,5)+=crLHS2*(crLHS104*r_DN(3,0) + crLHS106*r_DN(3,1) + crLHS109*crLHS313 + crLHS230 + crLHS315);
rLHS(13,6)+=crLHS2*(crLHS115*r_DN(3,0) + crLHS117*r_DN(3,1) + crLHS120*crLHS313 + crLHS242);
rLHS(13,7)+=crLHS314;
rLHS(13,8)+=crLHS316;
rLHS(13,9)+=crLHS2*(crLHS124*r_DN(3,0) + crLHS126*r_DN(3,1) + crLHS129*crLHS313 + crLHS281 + crLHS317);
rLHS(13,10)+=crLHS2*(crLHS135*r_DN(3,0) + crLHS137*r_DN(3,1) + crLHS140*crLHS313 + crLHS292);
rLHS(13,11)+=crLHS316;
rLHS(13,12)+=crLHS320;
rLHS(13,13)+=crLHS2*(crLHS144*r_DN(3,0) + crLHS146*r_DN(3,1) + crLHS148*crLHS313 + crLHS301*crLHS7 + crLHS321);
rLHS(13,14)+=crLHS2*(crLHS154*r_DN(3,0) + crLHS156*r_DN(3,1) + crLHS158*crLHS313 + crLHS322);
rLHS(13,15)+=crLHS320;
rLHS(14,0)+=crLHS323;
rLHS(14,1)+=crLHS2*(crLHS157 + crLHS160*r_DN(3,1) + crLHS285*crLHS324 + crLHS80*r_DN(3,0));
rLHS(14,2)+=crLHS2*(crLHS162*r_DN(3,1) + crLHS181 + crLHS312 + crLHS325*crLHS97 + crLHS95*r_DN(3,0));
rLHS(14,3)+=crLHS323;
rLHS(14,4)+=crLHS326;
rLHS(14,5)+=crLHS2*(crLHS106*r_DN(3,0) + crLHS109*crLHS325 + crLHS166*r_DN(3,1) + crLHS232);
rLHS(14,6)+=crLHS2*(crLHS117*r_DN(3,0) + crLHS120*crLHS325 + crLHS168*r_DN(3,1) + crLHS243 + crLHS315);
rLHS(14,7)+=crLHS326;
rLHS(14,8)+=crLHS327;
rLHS(14,9)+=crLHS2*(crLHS126*r_DN(3,0) + crLHS129*crLHS325 + crLHS172*r_DN(3,1) + crLHS283);
rLHS(14,10)+=crLHS2*(crLHS137*r_DN(3,0) + crLHS140*crLHS325 + crLHS174*r_DN(3,1) + crLHS293 + crLHS317);
rLHS(14,11)+=crLHS327;
rLHS(14,12)+=crLHS328;
rLHS(14,13)+=crLHS2*(crLHS146*r_DN(3,0) + crLHS148*crLHS325 + crLHS178*r_DN(3,1) + crLHS322);
rLHS(14,14)+=crLHS2*(crLHS156*r_DN(3,0) + crLHS158*crLHS325 + crLHS180*r_DN(3,1) + crLHS302*crLHS7 + crLHS321);
rLHS(14,15)+=crLHS328;
rLHS(15,0)+=crLHS57;
rLHS(15,1)+=crLHS294;
rLHS(15,2)+=crLHS296;
rLHS(15,3)+=crLHS57;
rLHS(15,4)+=crLHS201;
rLHS(15,5)+=crLHS297;
rLHS(15,6)+=crLHS298;
rLHS(15,7)+=crLHS201;
rLHS(15,8)+=crLHS257;
rLHS(15,9)+=crLHS299;
rLHS(15,10)+=crLHS300;
rLHS(15,11)+=crLHS257;
rLHS(15,12)+=crLHS303;
rLHS(15,13)+=crLHS305;
rLHS(15,14)+=crLHS306;
rLHS(15,15)+=crLHS303;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const auto& r_rho_n = rData.DensityOldStep1;
    const auto& r_rho_nn = rData.DensityOldStep2;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const auto& r_g = rData.BodyForce;
    const auto& r_heat_fl = rData.HeatFlux;

    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 3> u_conv = r_u - r_u_mesh;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2];
const double crRHS4 = 1.0/crRHS3;
const double crRHS5 = crRHS4*p_th;
const double crRHS6 = crRHS5*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]));
const double crRHS7 = -crRHS6 + dp_th_dt;
const double crRHS8 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2];
const double crRHS9 = crRHS8*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0));
const double crRHS10 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2];
const double crRHS11 = crRHS10*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1));
const double crRHS12 = crRHS5*(crRHS11 + crRHS9);
const double crRHS13 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0));
const double crRHS14 = r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0));
const double crRHS15 = r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS16 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crRHS17 = crRHS0*crRHS16;
const double crRHS18 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crRHS19 = crRHS18*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS20 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS21 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS22 = crRHS21*crRHS5;
const double crRHS23 = -crRHS22*(-crRHS13 - crRHS14 - crRHS15 - crRHS17 - crRHS19 + crRHS20) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2];
const double crRHS24 = stab_c2*sqrt(crRHS16*crRHS16 + crRHS18*crRHS18);
const double crRHS25 = 1.0*1.0/(crRHS24 + mu*stab_c1*1.0/(h*h));
const double crRHS26 = crRHS25*p_th;
const double crRHS27 = crRHS23*crRHS26;
const double crRHS28 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1));
const double crRHS29 = r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1));
const double crRHS30 = r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS31 = crRHS16*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1));
const double crRHS32 = crRHS1*crRHS18;
const double crRHS33 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS34 = -crRHS22*(-crRHS28 - crRHS29 - crRHS30 - crRHS31 - crRHS32 + crRHS33) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2];
const double crRHS35 = crRHS26*crRHS34;
const double crRHS36 = gauss_weight*gauss_weight;
const double crRHS37 = crRHS21*crRHS4;
const double crRHS38 = crRHS36*crRHS37;
const double crRHS39 = -crRHS38*(-crRHS12*r_N[0] + crRHS2*r_N[0] + crRHS27*r_DN(0,0) + crRHS35*r_DN(0,1) + crRHS7*r_N[0]);
const double crRHS40 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS41 = crRHS22*r_N[0];
const double crRHS42 = crRHS37*r_DN(0,0);
const double crRHS43 = crRHS13 + crRHS14 + crRHS15;
const double crRHS44 = crRHS17 + crRHS19;
const double crRHS45 = (crRHS11*crRHS5 - crRHS2 + crRHS5*crRHS9 + crRHS6 - dp_th_dt)*(crRHS24*h*1.0/stab_c1 + mu*1.0/rho_lin);
const double crRHS46 = crRHS25*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1));
const double crRHS47 = crRHS41*crRHS46;
const double crRHS48 = crRHS22*crRHS25;
const double crRHS49 = crRHS48*(crRHS16*r_DN(0,0) + crRHS18*r_DN(0,1));
const double crRHS50 = crRHS21*(crRHS10*crRHS18 + crRHS16*crRHS8)*1.0/(crRHS3*crRHS3);
const double crRHS51 = crRHS50*r_N[0];
const double crRHS52 = crRHS37*r_DN(0,1);
const double crRHS53 = crRHS28 + crRHS29 + crRHS30;
const double crRHS54 = crRHS31 + crRHS32;
const double crRHS55 = -crRHS38*(-crRHS12*r_N[1] + crRHS2*r_N[1] + crRHS27*r_DN(1,0) + crRHS35*r_DN(1,1) + crRHS7*r_N[1]);
const double crRHS56 = crRHS22*r_N[1];
const double crRHS57 = crRHS37*r_DN(1,0);
const double crRHS58 = crRHS46*crRHS56;
const double crRHS59 = crRHS48*(crRHS16*r_DN(1,0) + crRHS18*r_DN(1,1));
const double crRHS60 = crRHS50*r_N[1];
const double crRHS61 = crRHS37*r_DN(1,1);
const double crRHS62 = -crRHS38*(-crRHS12*r_N[2] + crRHS2*r_N[2] + crRHS27*r_DN(2,0) + crRHS35*r_DN(2,1) + crRHS7*r_N[2]);
const double crRHS63 = crRHS22*r_N[2];
const double crRHS64 = crRHS37*r_DN(2,0);
const double crRHS65 = crRHS46*crRHS63;
const double crRHS66 = crRHS48*(crRHS16*r_DN(2,0) + crRHS18*r_DN(2,1));
const double crRHS67 = crRHS50*r_N[2];
const double crRHS68 = crRHS37*r_DN(2,1);
rRHS[0]+=crRHS39;
rRHS[1]+=-crRHS36*(crRHS2*crRHS42 - crRHS20*crRHS41 + crRHS23*crRHS47 + crRHS23*crRHS49 - crRHS27*crRHS51 - crRHS40*r_DN(0,0) + crRHS41*crRHS43 + crRHS41*crRHS44 - crRHS42*crRHS45 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS36*(crRHS2*crRHS52 - crRHS33*crRHS41 + crRHS34*crRHS47 + crRHS34*crRHS49 - crRHS35*crRHS51 - crRHS40*r_DN(0,1) + crRHS41*crRHS53 + crRHS41*crRHS54 - crRHS45*crRHS52 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS39;
rRHS[4]+=crRHS55;
rRHS[5]+=-crRHS36*(crRHS2*crRHS57 - crRHS20*crRHS56 + crRHS23*crRHS58 + crRHS23*crRHS59 - crRHS27*crRHS60 - crRHS40*r_DN(1,0) + crRHS43*crRHS56 + crRHS44*crRHS56 - crRHS45*crRHS57 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS36*(crRHS2*crRHS61 - crRHS33*crRHS56 + crRHS34*crRHS58 + crRHS34*crRHS59 - crRHS35*crRHS60 - crRHS40*r_DN(1,1) - crRHS45*crRHS61 + crRHS53*crRHS56 + crRHS54*crRHS56 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS55;
rRHS[8]+=crRHS62;
rRHS[9]+=-crRHS36*(crRHS2*crRHS64 - crRHS20*crRHS63 + crRHS23*crRHS65 + crRHS23*crRHS66 - crRHS27*crRHS67 - crRHS40*r_DN(2,0) + crRHS43*crRHS63 + crRHS44*crRHS63 - crRHS45*crRHS64 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS36*(crRHS2*crRHS68 - crRHS33*crRHS63 + crRHS34*crRHS65 + crRHS34*crRHS66 - crRHS35*crRHS67 - crRHS40*r_DN(2,1) - crRHS45*crRHS68 + crRHS53*crRHS63 + crRHS54*crRHS63 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS62;

}

template <>
void LowMachNavierStokes<LowMachNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    LowMachNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    // Material parameters
    const auto& r_rho = rData.Density;
    const auto& r_rho_n = rData.DensityOldStep1;
    const auto& r_rho_nn = rData.DensityOldStep2;
    const double c_p = rData.SpecificHeat;
    const double kappa = rData.Conductivity;
    const double mu = rData.EffectiveViscosity;
    const double alpha = rData.ThermalExpansionCoefficient;

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
    const auto& r_g = rData.BodyForce;
    const auto& r_heat_fl = rData.HeatFlux;

    const auto& r_u_mesh = rData.MeshVelocity;
    const BoundedMatrix<double, 2, 4> u_conv = r_u - r_u_mesh;

    // Get shape function values
    const auto& r_N = rData.N;
    const auto& r_DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add RHS Gauss point contribution
    const double gauss_weight = rData.Weight;

    const double crRHS0 = r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0) + r_DN(3,0)*r_u(3,0);
const double crRHS1 = r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1) + r_DN(3,1)*r_u(3,1);
const double crRHS2 = p_th*(crRHS0 + crRHS1);
const double crRHS3 = r_N[0]*r_t[0] + r_N[1]*r_t[1] + r_N[2]*r_t[2] + r_N[3]*r_t[3];
const double crRHS4 = 1.0/crRHS3;
const double crRHS5 = crRHS4*p_th;
const double crRHS6 = crRHS5*(r_N[0]*(bdf0*r_t[0] + bdf1*r_t_n[0] + bdf2*r_t_nn[0]) + r_N[1]*(bdf0*r_t[1] + bdf1*r_t_n[1] + bdf2*r_t_nn[1]) + r_N[2]*(bdf0*r_t[2] + bdf1*r_t_n[2] + bdf2*r_t_nn[2]) + r_N[3]*(bdf0*r_t[3] + bdf1*r_t_n[3] + bdf2*r_t_nn[3]));
const double crRHS7 = -crRHS6 + dp_th_dt;
const double crRHS8 = r_DN(0,0)*r_t[0] + r_DN(1,0)*r_t[1] + r_DN(2,0)*r_t[2] + r_DN(3,0)*r_t[3];
const double crRHS9 = crRHS8*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0) + r_N[3]*r_u(3,0));
const double crRHS10 = r_DN(0,1)*r_t[0] + r_DN(1,1)*r_t[1] + r_DN(2,1)*r_t[2] + r_DN(3,1)*r_t[3];
const double crRHS11 = crRHS10*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1) + r_N[3]*r_u(3,1));
const double crRHS12 = crRHS5*(crRHS11 + crRHS9);
const double crRHS13 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0));
const double crRHS14 = r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0));
const double crRHS15 = r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS16 = r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS17 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crRHS18 = crRHS0*crRHS17;
const double crRHS19 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crRHS20 = crRHS19*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS21 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS22 = gamma*1.0/c_p/(gamma - 1.0);
const double crRHS23 = crRHS22*crRHS5;
const double crRHS24 = -crRHS23*(-crRHS13 - crRHS14 - crRHS15 - crRHS16 - crRHS18 - crRHS20 + crRHS21) + r_DN(0,0)*r_p[0] + r_DN(1,0)*r_p[1] + r_DN(2,0)*r_p[2] + r_DN(3,0)*r_p[3];
const double crRHS25 = stab_c2*sqrt(crRHS17*crRHS17 + crRHS19*crRHS19);
const double crRHS26 = 1.0*1.0/(crRHS25 + mu*stab_c1*1.0/(h*h));
const double crRHS27 = crRHS26*p_th;
const double crRHS28 = crRHS24*crRHS27;
const double crRHS29 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1));
const double crRHS30 = r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1));
const double crRHS31 = r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS32 = r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS33 = crRHS17*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1));
const double crRHS34 = crRHS1*crRHS19;
const double crRHS35 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS36 = -crRHS23*(-crRHS29 - crRHS30 - crRHS31 - crRHS32 - crRHS33 - crRHS34 + crRHS35) + r_DN(0,1)*r_p[0] + r_DN(1,1)*r_p[1] + r_DN(2,1)*r_p[2] + r_DN(3,1)*r_p[3];
const double crRHS37 = crRHS27*crRHS36;
const double crRHS38 = gauss_weight*gauss_weight;
const double crRHS39 = crRHS22*crRHS4;
const double crRHS40 = crRHS38*crRHS39;
const double crRHS41 = -crRHS40*(-crRHS12*r_N[0] + crRHS2*r_N[0] + crRHS28*r_DN(0,0) + crRHS37*r_DN(0,1) + crRHS7*r_N[0]);
const double crRHS42 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS43 = crRHS23*r_N[0];
const double crRHS44 = crRHS39*r_DN(0,0);
const double crRHS45 = crRHS13 + crRHS14 + crRHS15 + crRHS16;
const double crRHS46 = crRHS18 + crRHS20;
const double crRHS47 = (crRHS11*crRHS5 - crRHS2 + crRHS5*crRHS9 + crRHS6 - dp_th_dt)*(crRHS25*h*1.0/stab_c1 + mu*1.0/rho_lin);
const double crRHS48 = crRHS26*(r_DN(0,0)*u_conv(0,0) + r_DN(0,1)*u_conv(0,1) + r_DN(1,0)*u_conv(1,0) + r_DN(1,1)*u_conv(1,1) + r_DN(2,0)*u_conv(2,0) + r_DN(2,1)*u_conv(2,1) + r_DN(3,0)*u_conv(3,0) + r_DN(3,1)*u_conv(3,1));
const double crRHS49 = crRHS43*crRHS48;
const double crRHS50 = crRHS23*crRHS26;
const double crRHS51 = crRHS50*(crRHS17*r_DN(0,0) + crRHS19*r_DN(0,1));
const double crRHS52 = crRHS22*(crRHS10*crRHS19 + crRHS17*crRHS8)*1.0/(crRHS3*crRHS3);
const double crRHS53 = crRHS52*r_N[0];
const double crRHS54 = crRHS39*r_DN(0,1);
const double crRHS55 = crRHS29 + crRHS30 + crRHS31 + crRHS32;
const double crRHS56 = crRHS33 + crRHS34;
const double crRHS57 = -crRHS40*(-crRHS12*r_N[1] + crRHS2*r_N[1] + crRHS28*r_DN(1,0) + crRHS37*r_DN(1,1) + crRHS7*r_N[1]);
const double crRHS58 = crRHS23*r_N[1];
const double crRHS59 = crRHS39*r_DN(1,0);
const double crRHS60 = crRHS48*crRHS58;
const double crRHS61 = crRHS50*(crRHS17*r_DN(1,0) + crRHS19*r_DN(1,1));
const double crRHS62 = crRHS52*r_N[1];
const double crRHS63 = crRHS39*r_DN(1,1);
const double crRHS64 = -crRHS40*(-crRHS12*r_N[2] + crRHS2*r_N[2] + crRHS28*r_DN(2,0) + crRHS37*r_DN(2,1) + crRHS7*r_N[2]);
const double crRHS65 = crRHS23*r_N[2];
const double crRHS66 = crRHS39*r_DN(2,0);
const double crRHS67 = crRHS48*crRHS65;
const double crRHS68 = crRHS50*(crRHS17*r_DN(2,0) + crRHS19*r_DN(2,1));
const double crRHS69 = crRHS52*r_N[2];
const double crRHS70 = crRHS39*r_DN(2,1);
const double crRHS71 = -crRHS40*(-crRHS12*r_N[3] + crRHS2*r_N[3] + crRHS28*r_DN(3,0) + crRHS37*r_DN(3,1) + crRHS7*r_N[3]);
const double crRHS72 = crRHS23*r_N[3];
const double crRHS73 = crRHS39*r_DN(3,0);
const double crRHS74 = crRHS48*crRHS72;
const double crRHS75 = crRHS50*(crRHS17*r_DN(3,0) + crRHS19*r_DN(3,1));
const double crRHS76 = crRHS52*r_N[3];
const double crRHS77 = crRHS39*r_DN(3,1);
rRHS[0]+=crRHS41;
rRHS[1]+=-crRHS38*(crRHS2*crRHS44 - crRHS21*crRHS43 + crRHS24*crRHS49 + crRHS24*crRHS51 - crRHS28*crRHS53 - crRHS42*r_DN(0,0) + crRHS43*crRHS45 + crRHS43*crRHS46 - crRHS44*crRHS47 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS38*(crRHS2*crRHS54 - crRHS35*crRHS43 + crRHS36*crRHS49 + crRHS36*crRHS51 - crRHS37*crRHS53 - crRHS42*r_DN(0,1) + crRHS43*crRHS55 + crRHS43*crRHS56 - crRHS47*crRHS54 + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS41;
rRHS[4]+=crRHS57;
rRHS[5]+=-crRHS38*(crRHS2*crRHS59 - crRHS21*crRHS58 + crRHS24*crRHS60 + crRHS24*crRHS61 - crRHS28*crRHS62 - crRHS42*r_DN(1,0) + crRHS45*crRHS58 + crRHS46*crRHS58 - crRHS47*crRHS59 + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS38*(crRHS2*crRHS63 - crRHS35*crRHS58 + crRHS36*crRHS60 + crRHS36*crRHS61 - crRHS37*crRHS62 - crRHS42*r_DN(1,1) - crRHS47*crRHS63 + crRHS55*crRHS58 + crRHS56*crRHS58 + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS57;
rRHS[8]+=crRHS64;
rRHS[9]+=-crRHS38*(crRHS2*crRHS66 - crRHS21*crRHS65 + crRHS24*crRHS67 + crRHS24*crRHS68 - crRHS28*crRHS69 - crRHS42*r_DN(2,0) + crRHS45*crRHS65 + crRHS46*crRHS65 - crRHS47*crRHS66 + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS38*(crRHS2*crRHS70 - crRHS35*crRHS65 + crRHS36*crRHS67 + crRHS36*crRHS68 - crRHS37*crRHS69 - crRHS42*r_DN(2,1) - crRHS47*crRHS70 + crRHS55*crRHS65 + crRHS56*crRHS65 + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS64;
rRHS[12]+=crRHS71;
rRHS[13]+=-crRHS38*(crRHS2*crRHS73 - crRHS21*crRHS72 + crRHS24*crRHS74 + crRHS24*crRHS75 - crRHS28*crRHS76 - crRHS42*r_DN(3,0) + crRHS45*crRHS72 + crRHS46*crRHS72 - crRHS47*crRHS73 + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-crRHS38*(crRHS2*crRHS77 - crRHS35*crRHS72 + crRHS36*crRHS74 + crRHS36*crRHS75 - crRHS37*crRHS76 - crRHS42*r_DN(3,1) - crRHS47*crRHS77 + crRHS55*crRHS72 + crRHS56*crRHS72 + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=crRHS71;

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