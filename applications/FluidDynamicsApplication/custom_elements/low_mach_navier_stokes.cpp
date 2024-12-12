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

    const double crLHS0 = r_N[0]*r_rho[0] + r_N[1]*r_rho[1] + r_N[2]*r_rho[2];
const double crLHS1 = gauss_weight*gauss_weight;
const double crLHS2 = crLHS1*r_DN(0,0);
const double crLHS3 = crLHS2*r_N[0];
const double crLHS4 = -crLHS0*crLHS3;
const double crLHS5 = crLHS0*r_N[0];
const double crLHS6 = crLHS1*r_DN(0,1);
const double crLHS7 = -crLHS5*crLHS6;
const double crLHS8 = crLHS2*r_N[1];
const double crLHS9 = -crLHS0*crLHS8;
const double crLHS10 = crLHS6*r_N[1];
const double crLHS11 = -crLHS0*crLHS10;
const double crLHS12 = crLHS2*r_N[2];
const double crLHS13 = -crLHS0*crLHS12;
const double crLHS14 = crLHS6*r_N[2];
const double crLHS15 = -crLHS0*crLHS14;
const double crLHS16 = -crLHS3;
const double crLHS17 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS18 = r_C(0,2)*r_DN(0,0);
const double crLHS19 = crLHS18 + r_C(2,2)*r_DN(0,1);
const double crLHS20 = bdf0*crLHS0;
const double crLHS21 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crLHS22 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crLHS23 = crLHS21*r_DN(0,0) + crLHS22*r_DN(0,1);
const double crLHS24 = crLHS20*(r_N[0]*r_N[0]) + crLHS23*crLHS5;
const double crLHS25 = crLHS18 + r_C(0,1)*r_DN(0,1);
const double crLHS26 = r_C(1,2)*r_DN(0,1);
const double crLHS27 = crLHS26 + r_C(2,2)*r_DN(0,0);
const double crLHS28 = -crLHS8;
const double crLHS29 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS30 = r_C(0,2)*r_DN(1,0);
const double crLHS31 = crLHS30 + r_C(2,2)*r_DN(1,1);
const double crLHS32 = crLHS21*r_DN(1,0) + crLHS22*r_DN(1,1);
const double crLHS33 = crLHS20*r_N[0];
const double crLHS34 = crLHS33*r_N[1];
const double crLHS35 = crLHS32*crLHS5 + crLHS34;
const double crLHS36 = crLHS30 + r_C(0,1)*r_DN(1,1);
const double crLHS37 = r_C(1,2)*r_DN(1,1);
const double crLHS38 = crLHS37 + r_C(2,2)*r_DN(1,0);
const double crLHS39 = -crLHS12;
const double crLHS40 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS41 = r_C(0,2)*r_DN(2,0);
const double crLHS42 = crLHS41 + r_C(2,2)*r_DN(2,1);
const double crLHS43 = crLHS21*r_DN(2,0) + crLHS22*r_DN(2,1);
const double crLHS44 = crLHS33*r_N[2];
const double crLHS45 = crLHS43*crLHS5 + crLHS44;
const double crLHS46 = crLHS41 + r_C(0,1)*r_DN(2,1);
const double crLHS47 = r_C(1,2)*r_DN(2,1);
const double crLHS48 = crLHS47 + r_C(2,2)*r_DN(2,0);
const double crLHS49 = -crLHS6*r_N[0];
const double crLHS50 = crLHS26 + r_C(0,1)*r_DN(0,0);
const double crLHS51 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS52 = -crLHS10;
const double crLHS53 = crLHS37 + r_C(0,1)*r_DN(1,0);
const double crLHS54 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS55 = -crLHS14;
const double crLHS56 = crLHS47 + r_C(0,1)*r_DN(2,0);
const double crLHS57 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS58 = crLHS1*r_DN(1,0);
const double crLHS59 = -crLHS5*crLHS58;
const double crLHS60 = crLHS1*r_DN(1,1);
const double crLHS61 = -crLHS5*crLHS60;
const double crLHS62 = crLHS0*r_N[1];
const double crLHS63 = -crLHS58*crLHS62;
const double crLHS64 = -crLHS60*crLHS62;
const double crLHS65 = crLHS58*r_N[2];
const double crLHS66 = -crLHS0*crLHS65;
const double crLHS67 = crLHS60*r_N[2];
const double crLHS68 = -crLHS0*crLHS67;
const double crLHS69 = -crLHS58*r_N[0];
const double crLHS70 = crLHS23*crLHS62 + crLHS34;
const double crLHS71 = -crLHS58*r_N[1];
const double crLHS72 = crLHS20*(r_N[1]*r_N[1]) + crLHS32*crLHS62;
const double crLHS73 = -crLHS65;
const double crLHS74 = crLHS20*r_N[1]*r_N[2];
const double crLHS75 = crLHS43*crLHS62 + crLHS74;
const double crLHS76 = -crLHS60*r_N[0];
const double crLHS77 = -crLHS60*r_N[1];
const double crLHS78 = -crLHS67;
const double crLHS79 = crLHS1*r_DN(2,0);
const double crLHS80 = -crLHS5*crLHS79;
const double crLHS81 = crLHS1*r_DN(2,1);
const double crLHS82 = -crLHS5*crLHS81;
const double crLHS83 = -crLHS62*crLHS79;
const double crLHS84 = -crLHS62*crLHS81;
const double crLHS85 = crLHS0*r_N[2];
const double crLHS86 = -crLHS79*crLHS85;
const double crLHS87 = -crLHS81*crLHS85;
const double crLHS88 = -crLHS79*r_N[0];
const double crLHS89 = crLHS23*crLHS85 + crLHS44;
const double crLHS90 = -crLHS79*r_N[1];
const double crLHS91 = crLHS32*crLHS85 + crLHS74;
const double crLHS92 = -crLHS79*r_N[2];
const double crLHS93 = crLHS20*(r_N[2]*r_N[2]) + crLHS43*crLHS85;
const double crLHS94 = -crLHS81*r_N[0];
const double crLHS95 = -crLHS81*r_N[1];
const double crLHS96 = -crLHS81*r_N[2];
rLHS(0,0)+=0;
rLHS(0,1)+=crLHS4;
rLHS(0,2)+=crLHS7;
rLHS(0,3)+=0;
rLHS(0,4)+=0;
rLHS(0,5)+=crLHS9;
rLHS(0,6)+=crLHS11;
rLHS(0,7)+=0;
rLHS(0,8)+=0;
rLHS(0,9)+=crLHS13;
rLHS(0,10)+=crLHS15;
rLHS(0,11)+=0;
rLHS(1,0)+=crLHS16;
rLHS(1,1)+=crLHS1*(crLHS17*r_DN(0,0) + crLHS19*r_DN(0,1) + crLHS24);
rLHS(1,2)+=crLHS1*(crLHS25*r_DN(0,0) + crLHS27*r_DN(0,1));
rLHS(1,3)+=crLHS16;
rLHS(1,4)+=crLHS28;
rLHS(1,5)+=crLHS1*(crLHS29*r_DN(0,0) + crLHS31*r_DN(0,1) + crLHS35);
rLHS(1,6)+=crLHS1*(crLHS36*r_DN(0,0) + crLHS38*r_DN(0,1));
rLHS(1,7)+=crLHS28;
rLHS(1,8)+=crLHS39;
rLHS(1,9)+=crLHS1*(crLHS40*r_DN(0,0) + crLHS42*r_DN(0,1) + crLHS45);
rLHS(1,10)+=crLHS1*(crLHS46*r_DN(0,0) + crLHS48*r_DN(0,1));
rLHS(1,11)+=crLHS39;
rLHS(2,0)+=crLHS49;
rLHS(2,1)+=crLHS1*(crLHS19*r_DN(0,0) + crLHS50*r_DN(0,1));
rLHS(2,2)+=crLHS1*(crLHS24 + crLHS27*r_DN(0,0) + crLHS51*r_DN(0,1));
rLHS(2,3)+=crLHS49;
rLHS(2,4)+=crLHS52;
rLHS(2,5)+=crLHS1*(crLHS31*r_DN(0,0) + crLHS53*r_DN(0,1));
rLHS(2,6)+=crLHS1*(crLHS35 + crLHS38*r_DN(0,0) + crLHS54*r_DN(0,1));
rLHS(2,7)+=crLHS52;
rLHS(2,8)+=crLHS55;
rLHS(2,9)+=crLHS1*(crLHS42*r_DN(0,0) + crLHS56*r_DN(0,1));
rLHS(2,10)+=crLHS1*(crLHS45 + crLHS48*r_DN(0,0) + crLHS57*r_DN(0,1));
rLHS(2,11)+=crLHS55;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS4;
rLHS(3,2)+=crLHS7;
rLHS(3,3)+=0;
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS9;
rLHS(3,6)+=crLHS11;
rLHS(3,7)+=0;
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS13;
rLHS(3,10)+=crLHS15;
rLHS(3,11)+=0;
rLHS(4,0)+=0;
rLHS(4,1)+=crLHS59;
rLHS(4,2)+=crLHS61;
rLHS(4,3)+=0;
rLHS(4,4)+=0;
rLHS(4,5)+=crLHS63;
rLHS(4,6)+=crLHS64;
rLHS(4,7)+=0;
rLHS(4,8)+=0;
rLHS(4,9)+=crLHS66;
rLHS(4,10)+=crLHS68;
rLHS(4,11)+=0;
rLHS(5,0)+=crLHS69;
rLHS(5,1)+=crLHS1*(crLHS17*r_DN(1,0) + crLHS19*r_DN(1,1) + crLHS70);
rLHS(5,2)+=crLHS1*(crLHS25*r_DN(1,0) + crLHS27*r_DN(1,1));
rLHS(5,3)+=crLHS69;
rLHS(5,4)+=crLHS71;
rLHS(5,5)+=crLHS1*(crLHS29*r_DN(1,0) + crLHS31*r_DN(1,1) + crLHS72);
rLHS(5,6)+=crLHS1*(crLHS36*r_DN(1,0) + crLHS38*r_DN(1,1));
rLHS(5,7)+=crLHS71;
rLHS(5,8)+=crLHS73;
rLHS(5,9)+=crLHS1*(crLHS40*r_DN(1,0) + crLHS42*r_DN(1,1) + crLHS75);
rLHS(5,10)+=crLHS1*(crLHS46*r_DN(1,0) + crLHS48*r_DN(1,1));
rLHS(5,11)+=crLHS73;
rLHS(6,0)+=crLHS76;
rLHS(6,1)+=crLHS1*(crLHS19*r_DN(1,0) + crLHS50*r_DN(1,1));
rLHS(6,2)+=crLHS1*(crLHS27*r_DN(1,0) + crLHS51*r_DN(1,1) + crLHS70);
rLHS(6,3)+=crLHS76;
rLHS(6,4)+=crLHS77;
rLHS(6,5)+=crLHS1*(crLHS31*r_DN(1,0) + crLHS53*r_DN(1,1));
rLHS(6,6)+=crLHS1*(crLHS38*r_DN(1,0) + crLHS54*r_DN(1,1) + crLHS72);
rLHS(6,7)+=crLHS77;
rLHS(6,8)+=crLHS78;
rLHS(6,9)+=crLHS1*(crLHS42*r_DN(1,0) + crLHS56*r_DN(1,1));
rLHS(6,10)+=crLHS1*(crLHS48*r_DN(1,0) + crLHS57*r_DN(1,1) + crLHS75);
rLHS(6,11)+=crLHS78;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS59;
rLHS(7,2)+=crLHS61;
rLHS(7,3)+=0;
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS63;
rLHS(7,6)+=crLHS64;
rLHS(7,7)+=0;
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS66;
rLHS(7,10)+=crLHS68;
rLHS(7,11)+=0;
rLHS(8,0)+=0;
rLHS(8,1)+=crLHS80;
rLHS(8,2)+=crLHS82;
rLHS(8,3)+=0;
rLHS(8,4)+=0;
rLHS(8,5)+=crLHS83;
rLHS(8,6)+=crLHS84;
rLHS(8,7)+=0;
rLHS(8,8)+=0;
rLHS(8,9)+=crLHS86;
rLHS(8,10)+=crLHS87;
rLHS(8,11)+=0;
rLHS(9,0)+=crLHS88;
rLHS(9,1)+=crLHS1*(crLHS17*r_DN(2,0) + crLHS19*r_DN(2,1) + crLHS89);
rLHS(9,2)+=crLHS1*(crLHS25*r_DN(2,0) + crLHS27*r_DN(2,1));
rLHS(9,3)+=crLHS88;
rLHS(9,4)+=crLHS90;
rLHS(9,5)+=crLHS1*(crLHS29*r_DN(2,0) + crLHS31*r_DN(2,1) + crLHS91);
rLHS(9,6)+=crLHS1*(crLHS36*r_DN(2,0) + crLHS38*r_DN(2,1));
rLHS(9,7)+=crLHS90;
rLHS(9,8)+=crLHS92;
rLHS(9,9)+=crLHS1*(crLHS40*r_DN(2,0) + crLHS42*r_DN(2,1) + crLHS93);
rLHS(9,10)+=crLHS1*(crLHS46*r_DN(2,0) + crLHS48*r_DN(2,1));
rLHS(9,11)+=crLHS92;
rLHS(10,0)+=crLHS94;
rLHS(10,1)+=crLHS1*(crLHS19*r_DN(2,0) + crLHS50*r_DN(2,1));
rLHS(10,2)+=crLHS1*(crLHS27*r_DN(2,0) + crLHS51*r_DN(2,1) + crLHS89);
rLHS(10,3)+=crLHS94;
rLHS(10,4)+=crLHS95;
rLHS(10,5)+=crLHS1*(crLHS31*r_DN(2,0) + crLHS53*r_DN(2,1));
rLHS(10,6)+=crLHS1*(crLHS38*r_DN(2,0) + crLHS54*r_DN(2,1) + crLHS91);
rLHS(10,7)+=crLHS95;
rLHS(10,8)+=crLHS96;
rLHS(10,9)+=crLHS1*(crLHS42*r_DN(2,0) + crLHS56*r_DN(2,1));
rLHS(10,10)+=crLHS1*(crLHS48*r_DN(2,0) + crLHS57*r_DN(2,1) + crLHS93);
rLHS(10,11)+=crLHS96;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS80;
rLHS(11,2)+=crLHS82;
rLHS(11,3)+=0;
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS83;
rLHS(11,6)+=crLHS84;
rLHS(11,7)+=0;
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS86;
rLHS(11,10)+=crLHS87;
rLHS(11,11)+=0;

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

    const double crLHS0 = r_N[0]*r_rho[0] + r_N[1]*r_rho[1] + r_N[2]*r_rho[2] + r_N[3]*r_rho[3];
const double crLHS1 = gauss_weight*gauss_weight;
const double crLHS2 = crLHS1*r_DN(0,0);
const double crLHS3 = crLHS2*r_N[0];
const double crLHS4 = -crLHS0*crLHS3;
const double crLHS5 = crLHS0*r_N[0];
const double crLHS6 = crLHS1*r_DN(0,1);
const double crLHS7 = -crLHS5*crLHS6;
const double crLHS8 = crLHS2*r_N[1];
const double crLHS9 = -crLHS0*crLHS8;
const double crLHS10 = crLHS6*r_N[1];
const double crLHS11 = -crLHS0*crLHS10;
const double crLHS12 = crLHS2*r_N[2];
const double crLHS13 = -crLHS0*crLHS12;
const double crLHS14 = crLHS6*r_N[2];
const double crLHS15 = -crLHS0*crLHS14;
const double crLHS16 = crLHS2*r_N[3];
const double crLHS17 = -crLHS0*crLHS16;
const double crLHS18 = crLHS6*r_N[3];
const double crLHS19 = -crLHS0*crLHS18;
const double crLHS20 = -crLHS3;
const double crLHS21 = r_C(0,0)*r_DN(0,0) + r_C(0,2)*r_DN(0,1);
const double crLHS22 = r_C(0,2)*r_DN(0,0);
const double crLHS23 = crLHS22 + r_C(2,2)*r_DN(0,1);
const double crLHS24 = bdf0*crLHS0;
const double crLHS25 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crLHS26 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crLHS27 = crLHS25*r_DN(0,0) + crLHS26*r_DN(0,1);
const double crLHS28 = crLHS24*(r_N[0]*r_N[0]) + crLHS27*crLHS5;
const double crLHS29 = crLHS22 + r_C(0,1)*r_DN(0,1);
const double crLHS30 = r_C(1,2)*r_DN(0,1);
const double crLHS31 = crLHS30 + r_C(2,2)*r_DN(0,0);
const double crLHS32 = -crLHS8;
const double crLHS33 = r_C(0,0)*r_DN(1,0) + r_C(0,2)*r_DN(1,1);
const double crLHS34 = r_C(0,2)*r_DN(1,0);
const double crLHS35 = crLHS34 + r_C(2,2)*r_DN(1,1);
const double crLHS36 = crLHS25*r_DN(1,0) + crLHS26*r_DN(1,1);
const double crLHS37 = crLHS24*r_N[0];
const double crLHS38 = crLHS37*r_N[1];
const double crLHS39 = crLHS36*crLHS5 + crLHS38;
const double crLHS40 = crLHS34 + r_C(0,1)*r_DN(1,1);
const double crLHS41 = r_C(1,2)*r_DN(1,1);
const double crLHS42 = crLHS41 + r_C(2,2)*r_DN(1,0);
const double crLHS43 = -crLHS12;
const double crLHS44 = r_C(0,0)*r_DN(2,0) + r_C(0,2)*r_DN(2,1);
const double crLHS45 = r_C(0,2)*r_DN(2,0);
const double crLHS46 = crLHS45 + r_C(2,2)*r_DN(2,1);
const double crLHS47 = crLHS25*r_DN(2,0) + crLHS26*r_DN(2,1);
const double crLHS48 = crLHS37*r_N[2];
const double crLHS49 = crLHS47*crLHS5 + crLHS48;
const double crLHS50 = crLHS45 + r_C(0,1)*r_DN(2,1);
const double crLHS51 = r_C(1,2)*r_DN(2,1);
const double crLHS52 = crLHS51 + r_C(2,2)*r_DN(2,0);
const double crLHS53 = -crLHS16;
const double crLHS54 = r_C(0,0)*r_DN(3,0) + r_C(0,2)*r_DN(3,1);
const double crLHS55 = r_C(0,2)*r_DN(3,0);
const double crLHS56 = crLHS55 + r_C(2,2)*r_DN(3,1);
const double crLHS57 = crLHS25*r_DN(3,0) + crLHS26*r_DN(3,1);
const double crLHS58 = crLHS37*r_N[3];
const double crLHS59 = crLHS5*crLHS57 + crLHS58;
const double crLHS60 = crLHS55 + r_C(0,1)*r_DN(3,1);
const double crLHS61 = r_C(1,2)*r_DN(3,1);
const double crLHS62 = crLHS61 + r_C(2,2)*r_DN(3,0);
const double crLHS63 = -crLHS6*r_N[0];
const double crLHS64 = crLHS30 + r_C(0,1)*r_DN(0,0);
const double crLHS65 = r_C(1,1)*r_DN(0,1) + r_C(1,2)*r_DN(0,0);
const double crLHS66 = -crLHS10;
const double crLHS67 = crLHS41 + r_C(0,1)*r_DN(1,0);
const double crLHS68 = r_C(1,1)*r_DN(1,1) + r_C(1,2)*r_DN(1,0);
const double crLHS69 = -crLHS14;
const double crLHS70 = crLHS51 + r_C(0,1)*r_DN(2,0);
const double crLHS71 = r_C(1,1)*r_DN(2,1) + r_C(1,2)*r_DN(2,0);
const double crLHS72 = -crLHS18;
const double crLHS73 = crLHS61 + r_C(0,1)*r_DN(3,0);
const double crLHS74 = r_C(1,1)*r_DN(3,1) + r_C(1,2)*r_DN(3,0);
const double crLHS75 = crLHS1*r_DN(1,0);
const double crLHS76 = -crLHS5*crLHS75;
const double crLHS77 = crLHS1*r_DN(1,1);
const double crLHS78 = -crLHS5*crLHS77;
const double crLHS79 = crLHS0*r_N[1];
const double crLHS80 = -crLHS75*crLHS79;
const double crLHS81 = -crLHS77*crLHS79;
const double crLHS82 = crLHS75*r_N[2];
const double crLHS83 = -crLHS0*crLHS82;
const double crLHS84 = crLHS77*r_N[2];
const double crLHS85 = -crLHS0*crLHS84;
const double crLHS86 = crLHS75*r_N[3];
const double crLHS87 = -crLHS0*crLHS86;
const double crLHS88 = crLHS77*r_N[3];
const double crLHS89 = -crLHS0*crLHS88;
const double crLHS90 = -crLHS75*r_N[0];
const double crLHS91 = crLHS27*crLHS79 + crLHS38;
const double crLHS92 = -crLHS75*r_N[1];
const double crLHS93 = crLHS24*(r_N[1]*r_N[1]) + crLHS36*crLHS79;
const double crLHS94 = -crLHS82;
const double crLHS95 = crLHS24*r_N[1];
const double crLHS96 = crLHS95*r_N[2];
const double crLHS97 = crLHS47*crLHS79 + crLHS96;
const double crLHS98 = -crLHS86;
const double crLHS99 = crLHS95*r_N[3];
const double crLHS100 = crLHS57*crLHS79 + crLHS99;
const double crLHS101 = -crLHS77*r_N[0];
const double crLHS102 = -crLHS77*r_N[1];
const double crLHS103 = -crLHS84;
const double crLHS104 = -crLHS88;
const double crLHS105 = crLHS1*r_DN(2,0);
const double crLHS106 = -crLHS105*crLHS5;
const double crLHS107 = crLHS1*r_DN(2,1);
const double crLHS108 = -crLHS107*crLHS5;
const double crLHS109 = -crLHS105*crLHS79;
const double crLHS110 = -crLHS107*crLHS79;
const double crLHS111 = crLHS0*r_N[2];
const double crLHS112 = -crLHS105*crLHS111;
const double crLHS113 = -crLHS107*crLHS111;
const double crLHS114 = crLHS105*r_N[3];
const double crLHS115 = -crLHS0*crLHS114;
const double crLHS116 = crLHS107*r_N[3];
const double crLHS117 = -crLHS0*crLHS116;
const double crLHS118 = -crLHS105*r_N[0];
const double crLHS119 = crLHS111*crLHS27 + crLHS48;
const double crLHS120 = -crLHS105*r_N[1];
const double crLHS121 = crLHS111*crLHS36 + crLHS96;
const double crLHS122 = -crLHS105*r_N[2];
const double crLHS123 = crLHS111*crLHS47 + crLHS24*(r_N[2]*r_N[2]);
const double crLHS124 = -crLHS114;
const double crLHS125 = crLHS24*r_N[2]*r_N[3];
const double crLHS126 = crLHS111*crLHS57 + crLHS125;
const double crLHS127 = -crLHS107*r_N[0];
const double crLHS128 = -crLHS107*r_N[1];
const double crLHS129 = -crLHS107*r_N[2];
const double crLHS130 = -crLHS116;
const double crLHS131 = crLHS1*r_DN(3,0);
const double crLHS132 = -crLHS131*crLHS5;
const double crLHS133 = crLHS1*r_DN(3,1);
const double crLHS134 = -crLHS133*crLHS5;
const double crLHS135 = -crLHS131*crLHS79;
const double crLHS136 = -crLHS133*crLHS79;
const double crLHS137 = -crLHS111*crLHS131;
const double crLHS138 = -crLHS111*crLHS133;
const double crLHS139 = crLHS0*r_N[3];
const double crLHS140 = -crLHS131*crLHS139;
const double crLHS141 = -crLHS133*crLHS139;
const double crLHS142 = -crLHS131*r_N[0];
const double crLHS143 = crLHS139*crLHS27 + crLHS58;
const double crLHS144 = -crLHS131*r_N[1];
const double crLHS145 = crLHS139*crLHS36 + crLHS99;
const double crLHS146 = -crLHS131*r_N[2];
const double crLHS147 = crLHS125 + crLHS139*crLHS47;
const double crLHS148 = -crLHS131*r_N[3];
const double crLHS149 = crLHS139*crLHS57 + crLHS24*(r_N[3]*r_N[3]);
const double crLHS150 = -crLHS133*r_N[0];
const double crLHS151 = -crLHS133*r_N[1];
const double crLHS152 = -crLHS133*r_N[2];
const double crLHS153 = -crLHS133*r_N[3];
rLHS(0,0)+=0;
rLHS(0,1)+=crLHS4;
rLHS(0,2)+=crLHS7;
rLHS(0,3)+=0;
rLHS(0,4)+=0;
rLHS(0,5)+=crLHS9;
rLHS(0,6)+=crLHS11;
rLHS(0,7)+=0;
rLHS(0,8)+=0;
rLHS(0,9)+=crLHS13;
rLHS(0,10)+=crLHS15;
rLHS(0,11)+=0;
rLHS(0,12)+=0;
rLHS(0,13)+=crLHS17;
rLHS(0,14)+=crLHS19;
rLHS(0,15)+=0;
rLHS(1,0)+=crLHS20;
rLHS(1,1)+=crLHS1*(crLHS21*r_DN(0,0) + crLHS23*r_DN(0,1) + crLHS28);
rLHS(1,2)+=crLHS1*(crLHS29*r_DN(0,0) + crLHS31*r_DN(0,1));
rLHS(1,3)+=crLHS20;
rLHS(1,4)+=crLHS32;
rLHS(1,5)+=crLHS1*(crLHS33*r_DN(0,0) + crLHS35*r_DN(0,1) + crLHS39);
rLHS(1,6)+=crLHS1*(crLHS40*r_DN(0,0) + crLHS42*r_DN(0,1));
rLHS(1,7)+=crLHS32;
rLHS(1,8)+=crLHS43;
rLHS(1,9)+=crLHS1*(crLHS44*r_DN(0,0) + crLHS46*r_DN(0,1) + crLHS49);
rLHS(1,10)+=crLHS1*(crLHS50*r_DN(0,0) + crLHS52*r_DN(0,1));
rLHS(1,11)+=crLHS43;
rLHS(1,12)+=crLHS53;
rLHS(1,13)+=crLHS1*(crLHS54*r_DN(0,0) + crLHS56*r_DN(0,1) + crLHS59);
rLHS(1,14)+=crLHS1*(crLHS60*r_DN(0,0) + crLHS62*r_DN(0,1));
rLHS(1,15)+=crLHS53;
rLHS(2,0)+=crLHS63;
rLHS(2,1)+=crLHS1*(crLHS23*r_DN(0,0) + crLHS64*r_DN(0,1));
rLHS(2,2)+=crLHS1*(crLHS28 + crLHS31*r_DN(0,0) + crLHS65*r_DN(0,1));
rLHS(2,3)+=crLHS63;
rLHS(2,4)+=crLHS66;
rLHS(2,5)+=crLHS1*(crLHS35*r_DN(0,0) + crLHS67*r_DN(0,1));
rLHS(2,6)+=crLHS1*(crLHS39 + crLHS42*r_DN(0,0) + crLHS68*r_DN(0,1));
rLHS(2,7)+=crLHS66;
rLHS(2,8)+=crLHS69;
rLHS(2,9)+=crLHS1*(crLHS46*r_DN(0,0) + crLHS70*r_DN(0,1));
rLHS(2,10)+=crLHS1*(crLHS49 + crLHS52*r_DN(0,0) + crLHS71*r_DN(0,1));
rLHS(2,11)+=crLHS69;
rLHS(2,12)+=crLHS72;
rLHS(2,13)+=crLHS1*(crLHS56*r_DN(0,0) + crLHS73*r_DN(0,1));
rLHS(2,14)+=crLHS1*(crLHS59 + crLHS62*r_DN(0,0) + crLHS74*r_DN(0,1));
rLHS(2,15)+=crLHS72;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS4;
rLHS(3,2)+=crLHS7;
rLHS(3,3)+=0;
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS9;
rLHS(3,6)+=crLHS11;
rLHS(3,7)+=0;
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS13;
rLHS(3,10)+=crLHS15;
rLHS(3,11)+=0;
rLHS(3,12)+=0;
rLHS(3,13)+=crLHS17;
rLHS(3,14)+=crLHS19;
rLHS(3,15)+=0;
rLHS(4,0)+=0;
rLHS(4,1)+=crLHS76;
rLHS(4,2)+=crLHS78;
rLHS(4,3)+=0;
rLHS(4,4)+=0;
rLHS(4,5)+=crLHS80;
rLHS(4,6)+=crLHS81;
rLHS(4,7)+=0;
rLHS(4,8)+=0;
rLHS(4,9)+=crLHS83;
rLHS(4,10)+=crLHS85;
rLHS(4,11)+=0;
rLHS(4,12)+=0;
rLHS(4,13)+=crLHS87;
rLHS(4,14)+=crLHS89;
rLHS(4,15)+=0;
rLHS(5,0)+=crLHS90;
rLHS(5,1)+=crLHS1*(crLHS21*r_DN(1,0) + crLHS23*r_DN(1,1) + crLHS91);
rLHS(5,2)+=crLHS1*(crLHS29*r_DN(1,0) + crLHS31*r_DN(1,1));
rLHS(5,3)+=crLHS90;
rLHS(5,4)+=crLHS92;
rLHS(5,5)+=crLHS1*(crLHS33*r_DN(1,0) + crLHS35*r_DN(1,1) + crLHS93);
rLHS(5,6)+=crLHS1*(crLHS40*r_DN(1,0) + crLHS42*r_DN(1,1));
rLHS(5,7)+=crLHS92;
rLHS(5,8)+=crLHS94;
rLHS(5,9)+=crLHS1*(crLHS44*r_DN(1,0) + crLHS46*r_DN(1,1) + crLHS97);
rLHS(5,10)+=crLHS1*(crLHS50*r_DN(1,0) + crLHS52*r_DN(1,1));
rLHS(5,11)+=crLHS94;
rLHS(5,12)+=crLHS98;
rLHS(5,13)+=crLHS1*(crLHS100 + crLHS54*r_DN(1,0) + crLHS56*r_DN(1,1));
rLHS(5,14)+=crLHS1*(crLHS60*r_DN(1,0) + crLHS62*r_DN(1,1));
rLHS(5,15)+=crLHS98;
rLHS(6,0)+=crLHS101;
rLHS(6,1)+=crLHS1*(crLHS23*r_DN(1,0) + crLHS64*r_DN(1,1));
rLHS(6,2)+=crLHS1*(crLHS31*r_DN(1,0) + crLHS65*r_DN(1,1) + crLHS91);
rLHS(6,3)+=crLHS101;
rLHS(6,4)+=crLHS102;
rLHS(6,5)+=crLHS1*(crLHS35*r_DN(1,0) + crLHS67*r_DN(1,1));
rLHS(6,6)+=crLHS1*(crLHS42*r_DN(1,0) + crLHS68*r_DN(1,1) + crLHS93);
rLHS(6,7)+=crLHS102;
rLHS(6,8)+=crLHS103;
rLHS(6,9)+=crLHS1*(crLHS46*r_DN(1,0) + crLHS70*r_DN(1,1));
rLHS(6,10)+=crLHS1*(crLHS52*r_DN(1,0) + crLHS71*r_DN(1,1) + crLHS97);
rLHS(6,11)+=crLHS103;
rLHS(6,12)+=crLHS104;
rLHS(6,13)+=crLHS1*(crLHS56*r_DN(1,0) + crLHS73*r_DN(1,1));
rLHS(6,14)+=crLHS1*(crLHS100 + crLHS62*r_DN(1,0) + crLHS74*r_DN(1,1));
rLHS(6,15)+=crLHS104;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS76;
rLHS(7,2)+=crLHS78;
rLHS(7,3)+=0;
rLHS(7,4)+=0;
rLHS(7,5)+=crLHS80;
rLHS(7,6)+=crLHS81;
rLHS(7,7)+=0;
rLHS(7,8)+=0;
rLHS(7,9)+=crLHS83;
rLHS(7,10)+=crLHS85;
rLHS(7,11)+=0;
rLHS(7,12)+=0;
rLHS(7,13)+=crLHS87;
rLHS(7,14)+=crLHS89;
rLHS(7,15)+=0;
rLHS(8,0)+=0;
rLHS(8,1)+=crLHS106;
rLHS(8,2)+=crLHS108;
rLHS(8,3)+=0;
rLHS(8,4)+=0;
rLHS(8,5)+=crLHS109;
rLHS(8,6)+=crLHS110;
rLHS(8,7)+=0;
rLHS(8,8)+=0;
rLHS(8,9)+=crLHS112;
rLHS(8,10)+=crLHS113;
rLHS(8,11)+=0;
rLHS(8,12)+=0;
rLHS(8,13)+=crLHS115;
rLHS(8,14)+=crLHS117;
rLHS(8,15)+=0;
rLHS(9,0)+=crLHS118;
rLHS(9,1)+=crLHS1*(crLHS119 + crLHS21*r_DN(2,0) + crLHS23*r_DN(2,1));
rLHS(9,2)+=crLHS1*(crLHS29*r_DN(2,0) + crLHS31*r_DN(2,1));
rLHS(9,3)+=crLHS118;
rLHS(9,4)+=crLHS120;
rLHS(9,5)+=crLHS1*(crLHS121 + crLHS33*r_DN(2,0) + crLHS35*r_DN(2,1));
rLHS(9,6)+=crLHS1*(crLHS40*r_DN(2,0) + crLHS42*r_DN(2,1));
rLHS(9,7)+=crLHS120;
rLHS(9,8)+=crLHS122;
rLHS(9,9)+=crLHS1*(crLHS123 + crLHS44*r_DN(2,0) + crLHS46*r_DN(2,1));
rLHS(9,10)+=crLHS1*(crLHS50*r_DN(2,0) + crLHS52*r_DN(2,1));
rLHS(9,11)+=crLHS122;
rLHS(9,12)+=crLHS124;
rLHS(9,13)+=crLHS1*(crLHS126 + crLHS54*r_DN(2,0) + crLHS56*r_DN(2,1));
rLHS(9,14)+=crLHS1*(crLHS60*r_DN(2,0) + crLHS62*r_DN(2,1));
rLHS(9,15)+=crLHS124;
rLHS(10,0)+=crLHS127;
rLHS(10,1)+=crLHS1*(crLHS23*r_DN(2,0) + crLHS64*r_DN(2,1));
rLHS(10,2)+=crLHS1*(crLHS119 + crLHS31*r_DN(2,0) + crLHS65*r_DN(2,1));
rLHS(10,3)+=crLHS127;
rLHS(10,4)+=crLHS128;
rLHS(10,5)+=crLHS1*(crLHS35*r_DN(2,0) + crLHS67*r_DN(2,1));
rLHS(10,6)+=crLHS1*(crLHS121 + crLHS42*r_DN(2,0) + crLHS68*r_DN(2,1));
rLHS(10,7)+=crLHS128;
rLHS(10,8)+=crLHS129;
rLHS(10,9)+=crLHS1*(crLHS46*r_DN(2,0) + crLHS70*r_DN(2,1));
rLHS(10,10)+=crLHS1*(crLHS123 + crLHS52*r_DN(2,0) + crLHS71*r_DN(2,1));
rLHS(10,11)+=crLHS129;
rLHS(10,12)+=crLHS130;
rLHS(10,13)+=crLHS1*(crLHS56*r_DN(2,0) + crLHS73*r_DN(2,1));
rLHS(10,14)+=crLHS1*(crLHS126 + crLHS62*r_DN(2,0) + crLHS74*r_DN(2,1));
rLHS(10,15)+=crLHS130;
rLHS(11,0)+=0;
rLHS(11,1)+=crLHS106;
rLHS(11,2)+=crLHS108;
rLHS(11,3)+=0;
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS109;
rLHS(11,6)+=crLHS110;
rLHS(11,7)+=0;
rLHS(11,8)+=0;
rLHS(11,9)+=crLHS112;
rLHS(11,10)+=crLHS113;
rLHS(11,11)+=0;
rLHS(11,12)+=0;
rLHS(11,13)+=crLHS115;
rLHS(11,14)+=crLHS117;
rLHS(11,15)+=0;
rLHS(12,0)+=0;
rLHS(12,1)+=crLHS132;
rLHS(12,2)+=crLHS134;
rLHS(12,3)+=0;
rLHS(12,4)+=0;
rLHS(12,5)+=crLHS135;
rLHS(12,6)+=crLHS136;
rLHS(12,7)+=0;
rLHS(12,8)+=0;
rLHS(12,9)+=crLHS137;
rLHS(12,10)+=crLHS138;
rLHS(12,11)+=0;
rLHS(12,12)+=0;
rLHS(12,13)+=crLHS140;
rLHS(12,14)+=crLHS141;
rLHS(12,15)+=0;
rLHS(13,0)+=crLHS142;
rLHS(13,1)+=crLHS1*(crLHS143 + crLHS21*r_DN(3,0) + crLHS23*r_DN(3,1));
rLHS(13,2)+=crLHS1*(crLHS29*r_DN(3,0) + crLHS31*r_DN(3,1));
rLHS(13,3)+=crLHS142;
rLHS(13,4)+=crLHS144;
rLHS(13,5)+=crLHS1*(crLHS145 + crLHS33*r_DN(3,0) + crLHS35*r_DN(3,1));
rLHS(13,6)+=crLHS1*(crLHS40*r_DN(3,0) + crLHS42*r_DN(3,1));
rLHS(13,7)+=crLHS144;
rLHS(13,8)+=crLHS146;
rLHS(13,9)+=crLHS1*(crLHS147 + crLHS44*r_DN(3,0) + crLHS46*r_DN(3,1));
rLHS(13,10)+=crLHS1*(crLHS50*r_DN(3,0) + crLHS52*r_DN(3,1));
rLHS(13,11)+=crLHS146;
rLHS(13,12)+=crLHS148;
rLHS(13,13)+=crLHS1*(crLHS149 + crLHS54*r_DN(3,0) + crLHS56*r_DN(3,1));
rLHS(13,14)+=crLHS1*(crLHS60*r_DN(3,0) + crLHS62*r_DN(3,1));
rLHS(13,15)+=crLHS148;
rLHS(14,0)+=crLHS150;
rLHS(14,1)+=crLHS1*(crLHS23*r_DN(3,0) + crLHS64*r_DN(3,1));
rLHS(14,2)+=crLHS1*(crLHS143 + crLHS31*r_DN(3,0) + crLHS65*r_DN(3,1));
rLHS(14,3)+=crLHS150;
rLHS(14,4)+=crLHS151;
rLHS(14,5)+=crLHS1*(crLHS35*r_DN(3,0) + crLHS67*r_DN(3,1));
rLHS(14,6)+=crLHS1*(crLHS145 + crLHS42*r_DN(3,0) + crLHS68*r_DN(3,1));
rLHS(14,7)+=crLHS151;
rLHS(14,8)+=crLHS152;
rLHS(14,9)+=crLHS1*(crLHS46*r_DN(3,0) + crLHS70*r_DN(3,1));
rLHS(14,10)+=crLHS1*(crLHS147 + crLHS52*r_DN(3,0) + crLHS71*r_DN(3,1));
rLHS(14,11)+=crLHS152;
rLHS(14,12)+=crLHS153;
rLHS(14,13)+=crLHS1*(crLHS56*r_DN(3,0) + crLHS73*r_DN(3,1));
rLHS(14,14)+=crLHS1*(crLHS149 + crLHS62*r_DN(3,0) + crLHS74*r_DN(3,1));
rLHS(14,15)+=crLHS153;
rLHS(15,0)+=0;
rLHS(15,1)+=crLHS132;
rLHS(15,2)+=crLHS134;
rLHS(15,3)+=0;
rLHS(15,4)+=0;
rLHS(15,5)+=crLHS135;
rLHS(15,6)+=crLHS136;
rLHS(15,7)+=0;
rLHS(15,8)+=0;
rLHS(15,9)+=crLHS137;
rLHS(15,10)+=crLHS138;
rLHS(15,11)+=0;
rLHS(15,12)+=0;
rLHS(15,13)+=crLHS140;
rLHS(15,14)+=crLHS141;
rLHS(15,15)+=0;

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

    const double crRHS0 = gauss_weight*gauss_weight;
const double crRHS1 = r_N[0]*r_rho[0] + r_N[1]*r_rho[1] + r_N[2]*r_rho[2];
const double crRHS2 = crRHS1*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0));
const double crRHS3 = crRHS1*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1));
const double crRHS4 = r_N[0]*(bdf0*r_rho[0] + bdf1*r_rho_n[0] + bdf2*r_rho_nn[0]) + r_N[1]*(bdf0*r_rho[1] + bdf1*r_rho_n[1] + bdf2*r_rho_nn[1]) + r_N[2]*(bdf0*r_rho[2] + bdf1*r_rho_n[2] + bdf2*r_rho_nn[2]);
const double crRHS5 = crRHS0*(crRHS2*r_DN(0,0) + crRHS3*r_DN(0,1) - crRHS4*r_N[0]);
const double crRHS6 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2];
const double crRHS7 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0);
const double crRHS8 = crRHS1*r_N[0];
const double crRHS9 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0));
const double crRHS10 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0);
const double crRHS11 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1);
const double crRHS12 = crRHS10*(r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0)) + crRHS11*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0));
const double crRHS13 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1);
const double crRHS14 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1));
const double crRHS15 = crRHS10*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1)) + crRHS11*(r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1));
const double crRHS16 = crRHS0*(crRHS2*r_DN(1,0) + crRHS3*r_DN(1,1) - crRHS4*r_N[1]);
const double crRHS17 = crRHS1*r_N[1];
const double crRHS18 = crRHS0*(crRHS2*r_DN(2,0) + crRHS3*r_DN(2,1) - crRHS4*r_N[2]);
const double crRHS19 = crRHS1*r_N[2];
rRHS[0]+=crRHS5;
rRHS[1]+=-crRHS0*(crRHS12*crRHS8 - crRHS6*r_DN(0,0) - crRHS7*crRHS8 + crRHS8*crRHS9 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS0*(-crRHS13*crRHS8 + crRHS14*crRHS8 + crRHS15*crRHS8 - crRHS6*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS5;
rRHS[4]+=crRHS16;
rRHS[5]+=-crRHS0*(crRHS12*crRHS17 - crRHS17*crRHS7 + crRHS17*crRHS9 - crRHS6*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS0*(-crRHS13*crRHS17 + crRHS14*crRHS17 + crRHS15*crRHS17 - crRHS6*r_DN(1,1) + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS16;
rRHS[8]+=crRHS18;
rRHS[9]+=-crRHS0*(crRHS12*crRHS19 - crRHS19*crRHS7 + crRHS19*crRHS9 - crRHS6*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS0*(-crRHS13*crRHS19 + crRHS14*crRHS19 + crRHS15*crRHS19 - crRHS6*r_DN(2,1) + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS18;

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

    const double crRHS0 = gauss_weight*gauss_weight;
const double crRHS1 = r_N[0]*r_rho[0] + r_N[1]*r_rho[1] + r_N[2]*r_rho[2] + r_N[3]*r_rho[3];
const double crRHS2 = crRHS1*(r_N[0]*r_u(0,0) + r_N[1]*r_u(1,0) + r_N[2]*r_u(2,0) + r_N[3]*r_u(3,0));
const double crRHS3 = crRHS1*(r_N[0]*r_u(0,1) + r_N[1]*r_u(1,1) + r_N[2]*r_u(2,1) + r_N[3]*r_u(3,1));
const double crRHS4 = r_N[0]*(bdf0*r_rho[0] + bdf1*r_rho_n[0] + bdf2*r_rho_nn[0]) + r_N[1]*(bdf0*r_rho[1] + bdf1*r_rho_n[1] + bdf2*r_rho_nn[1]) + r_N[2]*(bdf0*r_rho[2] + bdf1*r_rho_n[2] + bdf2*r_rho_nn[2]) + r_N[3]*(bdf0*r_rho[3] + bdf1*r_rho_n[3] + bdf2*r_rho_nn[3]);
const double crRHS5 = crRHS0*(crRHS2*r_DN(0,0) + crRHS3*r_DN(0,1) - crRHS4*r_N[0]);
const double crRHS6 = r_N[0]*r_p[0] + r_N[1]*r_p[1] + r_N[2]*r_p[2] + r_N[3]*r_p[3];
const double crRHS7 = r_N[0]*r_g(0,0) + r_N[1]*r_g(1,0) + r_N[2]*r_g(2,0) + r_N[3]*r_g(3,0);
const double crRHS8 = crRHS1*r_N[0];
const double crRHS9 = r_N[0]*(bdf0*r_u(0,0) + bdf1*r_u_n(0,0) + bdf2*r_u_nn(0,0)) + r_N[1]*(bdf0*r_u(1,0) + bdf1*r_u_n(1,0) + bdf2*r_u_nn(1,0)) + r_N[2]*(bdf0*r_u(2,0) + bdf1*r_u_n(2,0) + bdf2*r_u_nn(2,0)) + r_N[3]*(bdf0*r_u(3,0) + bdf1*r_u_n(3,0) + bdf2*r_u_nn(3,0));
const double crRHS10 = r_N[0]*u_conv(0,0) + r_N[1]*u_conv(1,0) + r_N[2]*u_conv(2,0) + r_N[3]*u_conv(3,0);
const double crRHS11 = r_N[0]*u_conv(0,1) + r_N[1]*u_conv(1,1) + r_N[2]*u_conv(2,1) + r_N[3]*u_conv(3,1);
const double crRHS12 = crRHS10*(r_DN(0,0)*r_u(0,0) + r_DN(1,0)*r_u(1,0) + r_DN(2,0)*r_u(2,0) + r_DN(3,0)*r_u(3,0)) + crRHS11*(r_DN(0,1)*r_u(0,0) + r_DN(1,1)*r_u(1,0) + r_DN(2,1)*r_u(2,0) + r_DN(3,1)*r_u(3,0));
const double crRHS13 = r_N[0]*r_g(0,1) + r_N[1]*r_g(1,1) + r_N[2]*r_g(2,1) + r_N[3]*r_g(3,1);
const double crRHS14 = r_N[0]*(bdf0*r_u(0,1) + bdf1*r_u_n(0,1) + bdf2*r_u_nn(0,1)) + r_N[1]*(bdf0*r_u(1,1) + bdf1*r_u_n(1,1) + bdf2*r_u_nn(1,1)) + r_N[2]*(bdf0*r_u(2,1) + bdf1*r_u_n(2,1) + bdf2*r_u_nn(2,1)) + r_N[3]*(bdf0*r_u(3,1) + bdf1*r_u_n(3,1) + bdf2*r_u_nn(3,1));
const double crRHS15 = crRHS10*(r_DN(0,0)*r_u(0,1) + r_DN(1,0)*r_u(1,1) + r_DN(2,0)*r_u(2,1) + r_DN(3,0)*r_u(3,1)) + crRHS11*(r_DN(0,1)*r_u(0,1) + r_DN(1,1)*r_u(1,1) + r_DN(2,1)*r_u(2,1) + r_DN(3,1)*r_u(3,1));
const double crRHS16 = crRHS0*(crRHS2*r_DN(1,0) + crRHS3*r_DN(1,1) - crRHS4*r_N[1]);
const double crRHS17 = crRHS1*r_N[1];
const double crRHS18 = crRHS0*(crRHS2*r_DN(2,0) + crRHS3*r_DN(2,1) - crRHS4*r_N[2]);
const double crRHS19 = crRHS1*r_N[2];
const double crRHS20 = crRHS0*(crRHS2*r_DN(3,0) + crRHS3*r_DN(3,1) - crRHS4*r_N[3]);
const double crRHS21 = crRHS1*r_N[3];
rRHS[0]+=crRHS5;
rRHS[1]+=-crRHS0*(crRHS12*crRHS8 - crRHS6*r_DN(0,0) - crRHS7*crRHS8 + crRHS8*crRHS9 + r_DN(0,0)*r_stress[0] + r_DN(0,1)*r_stress[2]);
rRHS[2]+=-crRHS0*(-crRHS13*crRHS8 + crRHS14*crRHS8 + crRHS15*crRHS8 - crRHS6*r_DN(0,1) + r_DN(0,0)*r_stress[2] + r_DN(0,1)*r_stress[1]);
rRHS[3]+=crRHS5;
rRHS[4]+=crRHS16;
rRHS[5]+=-crRHS0*(crRHS12*crRHS17 - crRHS17*crRHS7 + crRHS17*crRHS9 - crRHS6*r_DN(1,0) + r_DN(1,0)*r_stress[0] + r_DN(1,1)*r_stress[2]);
rRHS[6]+=-crRHS0*(-crRHS13*crRHS17 + crRHS14*crRHS17 + crRHS15*crRHS17 - crRHS6*r_DN(1,1) + r_DN(1,0)*r_stress[2] + r_DN(1,1)*r_stress[1]);
rRHS[7]+=crRHS16;
rRHS[8]+=crRHS18;
rRHS[9]+=-crRHS0*(crRHS12*crRHS19 - crRHS19*crRHS7 + crRHS19*crRHS9 - crRHS6*r_DN(2,0) + r_DN(2,0)*r_stress[0] + r_DN(2,1)*r_stress[2]);
rRHS[10]+=-crRHS0*(-crRHS13*crRHS19 + crRHS14*crRHS19 + crRHS15*crRHS19 - crRHS6*r_DN(2,1) + r_DN(2,0)*r_stress[2] + r_DN(2,1)*r_stress[1]);
rRHS[11]+=crRHS18;
rRHS[12]+=crRHS20;
rRHS[13]+=-crRHS0*(crRHS12*crRHS21 - crRHS21*crRHS7 + crRHS21*crRHS9 - crRHS6*r_DN(3,0) + r_DN(3,0)*r_stress[0] + r_DN(3,1)*r_stress[2]);
rRHS[14]+=-crRHS0*(-crRHS13*crRHS21 + crRHS14*crRHS21 + crRHS15*crRHS21 - crRHS6*r_DN(3,1) + r_DN(3,0)*r_stress[2] + r_DN(3,1)*r_stress[1]);
rRHS[15]+=crRHS20;

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