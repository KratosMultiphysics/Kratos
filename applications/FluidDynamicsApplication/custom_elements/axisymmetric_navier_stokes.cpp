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
    : BaseType(NewId) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry) {}

template <class TElementData>
AxisymmetricNavierStokes<TElementData>::AxisymmetricNavierStokes(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties) {}

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

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Set the constitutive law pointer to null as the axisymmetric element hardcodes a Newtonian fluid viscous behavior
    // Note that to use a constitutive law the gradient in cylindrical coordinates would require the corresponding stress
    // implementation in cylindrical coordinates within the mechanical response calculation
    this->GetConstitutiveLaw() = nullptr;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int AxisymmetricNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::CalculateMaterialResponse(TElementData &rData) const
{
    rData.EffectiveViscosity = this->GetProperties()[DYNAMIC_VISCOSITY];
}

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

    // rData.lhs *= 2.0 * Globals::Pi * y * rData.Weight;
    // noalias(rLHS) -= rData.lhs;
    // noalias(rRHS) += prod(rData.lhs,values);
}

template <>
void AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,3> >::ComputeGaussPointLHSContribution(
    AxisymmetricNavierStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    // Axisymmetric formulation radius
    const double y = this->CalculateGaussPointRadius(rData);

    // Material parameters
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    // Nodal data
    const BoundedMatrix<double,2,3> v_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add LHS Gauss point contribution
    const double w_gauss = 2.0 * Globals::Pi * y * rData.Weight;

    const double crLHS0 = pow(DN(0,0), 2);
    const double crLHS1 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
    const double crLHS2 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
    const double crLHS3 = rho*stab_c2*sqrt(pow(crLHS1, 2) + pow(crLHS2, 2));
    const double crLHS4 = crLHS3*h/stab_c1 + mu;
    const double crLHS5 = bdf0*rho;
    const double crLHS6 = N[0]*crLHS5;
    const double crLHS7 = 1.0/y;
    const double crLHS8 = crLHS7*mu;
    const double crLHS9 = DN(0,1)*crLHS8;
    const double crLHS10 = DN(0,0)*crLHS1;
    const double crLHS11 = crLHS10*rho;
    const double crLHS12 = DN(0,1)*crLHS2;
    const double crLHS13 = crLHS12*rho;
    const double crLHS14 = crLHS11 + crLHS13 + crLHS6 - crLHS9;
    const double crLHS15 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
    const double crLHS16 = 1.0/(crLHS3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crLHS17 = 1.0*crLHS16;
    const double crLHS18 = crLHS17*rho;
    const double crLHS19 = crLHS18*(N[0]*crLHS15 + crLHS10);
    const double crLHS20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
    const double crLHS21 = crLHS18*(N[0]*crLHS20 + crLHS12);
    const double crLHS22 = pow(DN(0,1), 2);
    const double crLHS23 = pow(N[0], 2);
    const double crLHS24 = N[0]*crLHS11 + N[0]*crLHS13 + crLHS0*mu + crLHS22*mu + crLHS23*crLHS5;
    const double crLHS25 = N[0]*crLHS7;
    const double crLHS26 = DN(0,1) + crLHS25;
    const double crLHS27 = DN(0,0)*w_gauss;
    const double crLHS28 = crLHS27*crLHS4;
    const double crLHS29 = crLHS26*crLHS28;
    const double crLHS30 = N[1]*crLHS5;
    const double crLHS31 = DN(1,1)*crLHS8;
    const double crLHS32 = DN(1,0)*crLHS1;
    const double crLHS33 = crLHS32*rho;
    const double crLHS34 = DN(1,1)*crLHS2;
    const double crLHS35 = crLHS34*rho;
    const double crLHS36 = crLHS30 - crLHS31 + crLHS33 + crLHS35;
    const double crLHS37 = DN(0,0)*DN(1,0);
    const double crLHS38 = DN(0,1)*DN(1,1);
    const double crLHS39 = N[1]*crLHS6 + crLHS37*mu + crLHS38*mu;
    const double crLHS40 = crLHS37*crLHS4 + crLHS39;
    const double crLHS41 = DN(1,0)*N[0];
    const double crLHS42 = crLHS1*rho;
    const double crLHS43 = DN(1,1)*N[0];
    const double crLHS44 = crLHS2*rho;
    const double crLHS45 = crLHS41*crLHS42 + crLHS43*crLHS44;
    const double crLHS46 = N[1]*crLHS7;
    const double crLHS47 = DN(1,1) + crLHS46;
    const double crLHS48 = crLHS28*crLHS47;
    const double crLHS49 = DN(0,0)*N[1];
    const double crLHS50 = N[2]*crLHS5;
    const double crLHS51 = DN(2,1)*crLHS8;
    const double crLHS52 = DN(2,0)*crLHS1;
    const double crLHS53 = crLHS52*rho;
    const double crLHS54 = DN(2,1)*crLHS2;
    const double crLHS55 = crLHS54*rho;
    const double crLHS56 = crLHS50 - crLHS51 + crLHS53 + crLHS55;
    const double crLHS57 = DN(0,0)*DN(2,0);
    const double crLHS58 = DN(0,1)*DN(2,1);
    const double crLHS59 = N[2]*crLHS6 + crLHS57*mu + crLHS58*mu;
    const double crLHS60 = crLHS4*crLHS57 + crLHS59;
    const double crLHS61 = DN(2,0)*N[0];
    const double crLHS62 = DN(2,1)*N[0];
    const double crLHS63 = crLHS42*crLHS61 + crLHS44*crLHS62;
    const double crLHS64 = N[2]*crLHS7;
    const double crLHS65 = DN(2,1) + crLHS64;
    const double crLHS66 = crLHS28*crLHS65;
    const double crLHS67 = DN(0,0)*N[2];
    const double crLHS68 = DN(0,1)*crLHS4;
    const double crLHS69 = crLHS26*crLHS4;
    const double crLHS70 = mu/pow(y, 2);
    const double crLHS71 = N[0]*crLHS70 + crLHS14;
    const double crLHS72 = DN(1,0)*w_gauss;
    const double crLHS73 = crLHS69*crLHS72;
    const double crLHS74 = crLHS25*crLHS4;
    const double crLHS75 = N[1]*crLHS70 + crLHS36;
    const double crLHS76 = DN(2,0)*w_gauss;
    const double crLHS77 = crLHS69*crLHS76;
    const double crLHS78 = N[2]*crLHS70 + crLHS56;
    const double crLHS79 = DN(0,1)*crLHS17;
    const double crLHS80 = crLHS17*crLHS25;
    const double crLHS81 = crLHS17*w_gauss;
    const double crLHS82 = DN(0,0)*crLHS17;
    const double crLHS83 = N[1]*crLHS25;
    const double crLHS84 = crLHS37 + crLHS38;
    const double crLHS85 = N[2]*crLHS25;
    const double crLHS86 = crLHS57 + crLHS58;
    const double crLHS87 = crLHS18*(N[1]*crLHS15 + crLHS32);
    const double crLHS88 = crLHS18*(N[1]*crLHS20 + crLHS34);
    const double crLHS89 = N[1]*crLHS11 + N[1]*crLHS13;
    const double crLHS90 = pow(DN(1,0), 2);
    const double crLHS91 = pow(DN(1,1), 2);
    const double crLHS92 = pow(N[1], 2);
    const double crLHS93 = N[1]*crLHS33 + N[1]*crLHS35 + crLHS5*crLHS92 + crLHS90*mu + crLHS91*mu;
    const double crLHS94 = crLHS4*crLHS72;
    const double crLHS95 = crLHS47*crLHS94;
    const double crLHS96 = DN(1,0)*DN(2,0);
    const double crLHS97 = DN(1,1)*DN(2,1);
    const double crLHS98 = N[2]*crLHS30 + crLHS96*mu + crLHS97*mu;
    const double crLHS99 = crLHS4*crLHS96 + crLHS98;
    const double crLHS100 = DN(2,0)*N[1];
    const double crLHS101 = DN(2,1)*N[1];
    const double crLHS102 = crLHS100*crLHS42 + crLHS101*crLHS44;
    const double crLHS103 = crLHS65*crLHS94;
    const double crLHS104 = DN(1,0)*N[2];
    const double crLHS105 = DN(1,1)*crLHS4;
    const double crLHS106 = crLHS4*crLHS47;
    const double crLHS107 = crLHS106*crLHS76;
    const double crLHS108 = crLHS4*crLHS65;
    const double crLHS109 = DN(1,0)*crLHS17;
    const double crLHS110 = DN(1,1)*crLHS17;
    const double crLHS111 = crLHS17*crLHS46;
    const double crLHS112 = N[2]*crLHS46;
    const double crLHS113 = crLHS96 + crLHS97;
    const double crLHS114 = crLHS18*(N[2]*crLHS15 + crLHS52);
    const double crLHS115 = crLHS18*(N[2]*crLHS20 + crLHS54);
    const double crLHS116 = N[2]*crLHS11 + N[2]*crLHS13;
    const double crLHS117 = N[2]*crLHS33 + N[2]*crLHS35;
    const double crLHS118 = pow(DN(2,0), 2);
    const double crLHS119 = pow(DN(2,1), 2);
    const double crLHS120 = pow(N[2], 2);
    const double crLHS121 = N[2]*crLHS53 + N[2]*crLHS55 + crLHS118*mu + crLHS119*mu + crLHS120*crLHS5;
    const double crLHS122 = crLHS108*crLHS76;
    const double crLHS123 = DN(2,1)*crLHS4;
    const double crLHS124 = DN(2,0)*crLHS17;
    const double crLHS125 = DN(2,1)*crLHS17;
    const double crLHS126 = crLHS17*crLHS64;
    rLHS(0,0)+=w_gauss*(crLHS0*crLHS4 + crLHS14*crLHS19 + crLHS14*crLHS21 + crLHS24);
    rLHS(0,1)+=crLHS29;
    rLHS(0,2)+=crLHS27*(-N[0] + crLHS19 + crLHS21);
    rLHS(0,3)+=w_gauss*(crLHS19*crLHS36 + crLHS21*crLHS36 + crLHS40 + crLHS45);
    rLHS(0,4)+=crLHS48;
    rLHS(0,5)+=w_gauss*(DN(1,0)*crLHS19 + DN(1,0)*crLHS21 - crLHS49);
    rLHS(0,6)+=w_gauss*(crLHS19*crLHS56 + crLHS21*crLHS56 + crLHS60 + crLHS63);
    rLHS(0,7)+=crLHS66;
    rLHS(0,8)+=w_gauss*(DN(2,0)*crLHS19 + DN(2,0)*crLHS21 - crLHS67);
    rLHS(1,0)+=crLHS29;
    rLHS(1,1)+=w_gauss*(crLHS19*crLHS71 + crLHS21*crLHS71 + crLHS24 + crLHS25*crLHS69 + crLHS26*crLHS68);
    rLHS(1,2)+=w_gauss*(DN(0,1)*crLHS19 + DN(0,1)*crLHS21 - N[0]*crLHS26);
    rLHS(1,3)+=crLHS73;
    rLHS(1,4)+=w_gauss*(crLHS19*crLHS75 + crLHS21*crLHS75 + crLHS39 + crLHS45 + crLHS47*crLHS68 + crLHS47*crLHS74);
    rLHS(1,5)+=w_gauss*(DN(1,1)*crLHS19 + DN(1,1)*crLHS21 - N[1]*crLHS26);
    rLHS(1,6)+=crLHS77;
    rLHS(1,7)+=w_gauss*(crLHS19*crLHS78 + crLHS21*crLHS78 + crLHS59 + crLHS63 + crLHS65*crLHS68 + crLHS65*crLHS74);
    rLHS(1,8)+=w_gauss*(DN(2,1)*crLHS19 + DN(2,1)*crLHS21 - N[2]*crLHS26);
    rLHS(2,0)+=crLHS27*(N[0] + crLHS16*(1.0*crLHS11 + 1.0*crLHS13 + 1.0*crLHS6 - 1.0*crLHS9));
    rLHS(2,1)+=w_gauss*(DN(0,1)*N[0] + crLHS23*crLHS7 + crLHS71*crLHS79 - crLHS71*crLHS80);
    rLHS(2,2)+=crLHS81*(-DN(0,1)*crLHS25 + crLHS0 + crLHS22);
    rLHS(2,3)+=w_gauss*(crLHS36*crLHS82 + crLHS41);
    rLHS(2,4)+=w_gauss*(crLHS43 + crLHS75*crLHS79 - crLHS75*crLHS80 + crLHS83);
    rLHS(2,5)+=crLHS81*(-DN(1,1)*crLHS25 + crLHS84);
    rLHS(2,6)+=w_gauss*(crLHS56*crLHS82 + crLHS61);
    rLHS(2,7)+=w_gauss*(crLHS62 + crLHS78*crLHS79 - crLHS78*crLHS80 + crLHS85);
    rLHS(2,8)+=crLHS81*(-DN(2,1)*crLHS25 + crLHS86);
    rLHS(3,0)+=w_gauss*(crLHS14*crLHS87 + crLHS14*crLHS88 + crLHS40 + crLHS89);
    rLHS(3,1)+=crLHS73;
    rLHS(3,2)+=w_gauss*(DN(0,0)*crLHS87 + DN(0,0)*crLHS88 - crLHS41);
    rLHS(3,3)+=w_gauss*(crLHS36*crLHS87 + crLHS36*crLHS88 + crLHS4*crLHS90 + crLHS93);
    rLHS(3,4)+=crLHS95;
    rLHS(3,5)+=crLHS72*(-N[1] + crLHS87 + crLHS88);
    rLHS(3,6)+=w_gauss*(crLHS102 + crLHS56*crLHS87 + crLHS56*crLHS88 + crLHS99);
    rLHS(3,7)+=crLHS103;
    rLHS(3,8)+=w_gauss*(DN(2,0)*crLHS87 + DN(2,0)*crLHS88 - crLHS104);
    rLHS(4,0)+=crLHS48;
    rLHS(4,1)+=w_gauss*(crLHS105*crLHS26 + crLHS39 + crLHS46*crLHS69 + crLHS71*crLHS87 + crLHS71*crLHS88 + crLHS89);
    rLHS(4,2)+=w_gauss*(DN(0,1)*crLHS87 + DN(0,1)*crLHS88 - N[0]*crLHS47);
    rLHS(4,3)+=crLHS95;
    rLHS(4,4)+=w_gauss*(crLHS105*crLHS47 + crLHS106*crLHS46 + crLHS75*crLHS87 + crLHS75*crLHS88 + crLHS93);
    rLHS(4,5)+=w_gauss*(DN(1,1)*crLHS87 + DN(1,1)*crLHS88 - N[1]*crLHS47);
    rLHS(4,6)+=crLHS107;
    rLHS(4,7)+=w_gauss*(crLHS102 + crLHS105*crLHS65 + crLHS108*crLHS46 + crLHS78*crLHS87 + crLHS78*crLHS88 + crLHS98);
    rLHS(4,8)+=w_gauss*(DN(2,1)*crLHS87 + DN(2,1)*crLHS88 - N[2]*crLHS47);
    rLHS(5,0)+=w_gauss*(crLHS109*crLHS14 + crLHS49);
    rLHS(5,1)+=w_gauss*(DN(0,1)*N[1] + crLHS110*crLHS71 - crLHS111*crLHS71 + crLHS83);
    rLHS(5,2)+=crLHS81*(-DN(0,1)*crLHS46 + crLHS84);
    rLHS(5,3)+=crLHS72*(N[1] + crLHS16*(1.0*crLHS30 - 1.0*crLHS31 + 1.0*crLHS33 + 1.0*crLHS35));
    rLHS(5,4)+=w_gauss*(DN(1,1)*N[1] + crLHS110*crLHS75 - crLHS111*crLHS75 + crLHS7*crLHS92);
    rLHS(5,5)+=crLHS81*(-DN(1,1)*crLHS46 + crLHS90 + crLHS91);
    rLHS(5,6)+=w_gauss*(crLHS100 + crLHS109*crLHS56);
    rLHS(5,7)+=w_gauss*(crLHS101 + crLHS110*crLHS78 - crLHS111*crLHS78 + crLHS112);
    rLHS(5,8)+=crLHS81*(-DN(2,1)*crLHS46 + crLHS113);
    rLHS(6,0)+=w_gauss*(crLHS114*crLHS14 + crLHS115*crLHS14 + crLHS116 + crLHS60);
    rLHS(6,1)+=crLHS77;
    rLHS(6,2)+=w_gauss*(DN(0,0)*crLHS114 + DN(0,0)*crLHS115 - crLHS61);
    rLHS(6,3)+=w_gauss*(crLHS114*crLHS36 + crLHS115*crLHS36 + crLHS117 + crLHS99);
    rLHS(6,4)+=crLHS107;
    rLHS(6,5)+=w_gauss*(DN(1,0)*crLHS114 + DN(1,0)*crLHS115 - crLHS100);
    rLHS(6,6)+=w_gauss*(crLHS114*crLHS56 + crLHS115*crLHS56 + crLHS118*crLHS4 + crLHS121);
    rLHS(6,7)+=crLHS122;
    rLHS(6,8)+=crLHS76*(-N[2] + crLHS114 + crLHS115);
    rLHS(7,0)+=crLHS66;
    rLHS(7,1)+=w_gauss*(crLHS114*crLHS71 + crLHS115*crLHS71 + crLHS116 + crLHS123*crLHS26 + crLHS59 + crLHS64*crLHS69);
    rLHS(7,2)+=w_gauss*(DN(0,1)*crLHS114 + DN(0,1)*crLHS115 - N[0]*crLHS65);
    rLHS(7,3)+=crLHS103;
    rLHS(7,4)+=w_gauss*(crLHS106*crLHS64 + crLHS114*crLHS75 + crLHS115*crLHS75 + crLHS117 + crLHS123*crLHS47 + crLHS98);
    rLHS(7,5)+=w_gauss*(DN(1,1)*crLHS114 + DN(1,1)*crLHS115 - N[1]*crLHS65);
    rLHS(7,6)+=crLHS122;
    rLHS(7,7)+=w_gauss*(crLHS108*crLHS64 + crLHS114*crLHS78 + crLHS115*crLHS78 + crLHS121 + crLHS123*crLHS65);
    rLHS(7,8)+=w_gauss*(DN(2,1)*crLHS114 + DN(2,1)*crLHS115 - N[2]*crLHS65);
    rLHS(8,0)+=w_gauss*(crLHS124*crLHS14 + crLHS67);
    rLHS(8,1)+=w_gauss*(DN(0,1)*N[2] + crLHS125*crLHS71 - crLHS126*crLHS71 + crLHS85);
    rLHS(8,2)+=crLHS81*(-DN(0,1)*crLHS64 + crLHS86);
    rLHS(8,3)+=w_gauss*(crLHS104 + crLHS124*crLHS36);
    rLHS(8,4)+=w_gauss*(DN(1,1)*N[2] + crLHS112 + crLHS125*crLHS75 - crLHS126*crLHS75);
    rLHS(8,5)+=crLHS81*(-DN(1,1)*crLHS64 + crLHS113);
    rLHS(8,6)+=crLHS76*(N[2] + crLHS16*(1.0*crLHS50 - 1.0*crLHS51 + 1.0*crLHS53 + 1.0*crLHS55));
    rLHS(8,7)+=w_gauss*(DN(2,1)*N[2] + crLHS120*crLHS7 + crLHS125*crLHS78 - crLHS126*crLHS78);
    rLHS(8,8)+=crLHS81*(-DN(2,1)*crLHS64 + crLHS118 + crLHS119);

}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    AxisymmetricNavierStokesData<2,4>& rData,
    MatrixType& rLHS)
{
    // Axisymmetric formulation radius
    const double y = this->CalculateGaussPointRadius(rData);

    // Material parameters
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    // Nodal data
    const BoundedMatrix<double,2,4> v_conv = rData.Velocity - rData.MeshVelocity;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add LHS Gauss point contribution
    const double w_gauss = 2.0 * Globals::Pi * y * rData.Weight;

    const double crLHS0 = pow(DN(0,0), 2);
    const double crLHS1 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
    const double crLHS2 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
    const double crLHS3 = rho*stab_c2*sqrt(pow(crLHS1, 2) + pow(crLHS2, 2));
    const double crLHS4 = crLHS3*h/stab_c1 + mu;
    const double crLHS5 = bdf0*rho;
    const double crLHS6 = N[0]*crLHS5;
    const double crLHS7 = 1.0/y;
    const double crLHS8 = crLHS7*mu;
    const double crLHS9 = DN(0,1)*crLHS8;
    const double crLHS10 = DN(0,0)*crLHS1;
    const double crLHS11 = crLHS10*rho;
    const double crLHS12 = DN(0,1)*crLHS2;
    const double crLHS13 = crLHS12*rho;
    const double crLHS14 = crLHS11 + crLHS13 + crLHS6 - crLHS9;
    const double crLHS15 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
    const double crLHS16 = 1.0/(crLHS3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crLHS17 = 1.0*crLHS16;
    const double crLHS18 = crLHS17*rho;
    const double crLHS19 = crLHS18*(N[0]*crLHS15 + crLHS10);
    const double crLHS20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
    const double crLHS21 = crLHS18*(N[0]*crLHS20 + crLHS12);
    const double crLHS22 = pow(DN(0,1), 2);
    const double crLHS23 = pow(N[0], 2);
    const double crLHS24 = N[0]*crLHS11 + N[0]*crLHS13 + crLHS0*mu + crLHS22*mu + crLHS23*crLHS5;
    const double crLHS25 = N[0]*crLHS7;
    const double crLHS26 = DN(0,1) + crLHS25;
    const double crLHS27 = DN(0,0)*w_gauss;
    const double crLHS28 = crLHS27*crLHS4;
    const double crLHS29 = crLHS26*crLHS28;
    const double crLHS30 = N[1]*crLHS5;
    const double crLHS31 = DN(1,1)*crLHS8;
    const double crLHS32 = DN(1,0)*crLHS1;
    const double crLHS33 = crLHS32*rho;
    const double crLHS34 = DN(1,1)*crLHS2;
    const double crLHS35 = crLHS34*rho;
    const double crLHS36 = crLHS30 - crLHS31 + crLHS33 + crLHS35;
    const double crLHS37 = DN(0,0)*DN(1,0);
    const double crLHS38 = DN(0,1)*DN(1,1);
    const double crLHS39 = N[1]*crLHS6 + crLHS37*mu + crLHS38*mu;
    const double crLHS40 = crLHS37*crLHS4 + crLHS39;
    const double crLHS41 = DN(1,0)*N[0];
    const double crLHS42 = crLHS1*rho;
    const double crLHS43 = DN(1,1)*N[0];
    const double crLHS44 = crLHS2*rho;
    const double crLHS45 = crLHS41*crLHS42 + crLHS43*crLHS44;
    const double crLHS46 = N[1]*crLHS7;
    const double crLHS47 = DN(1,1) + crLHS46;
    const double crLHS48 = crLHS28*crLHS47;
    const double crLHS49 = DN(0,0)*N[1];
    const double crLHS50 = N[2]*crLHS5;
    const double crLHS51 = DN(2,1)*crLHS8;
    const double crLHS52 = DN(2,0)*crLHS1;
    const double crLHS53 = crLHS52*rho;
    const double crLHS54 = DN(2,1)*crLHS2;
    const double crLHS55 = crLHS54*rho;
    const double crLHS56 = crLHS50 - crLHS51 + crLHS53 + crLHS55;
    const double crLHS57 = DN(0,0)*DN(2,0);
    const double crLHS58 = DN(0,1)*DN(2,1);
    const double crLHS59 = N[2]*crLHS6 + crLHS57*mu + crLHS58*mu;
    const double crLHS60 = crLHS4*crLHS57 + crLHS59;
    const double crLHS61 = DN(2,0)*N[0];
    const double crLHS62 = DN(2,1)*N[0];
    const double crLHS63 = crLHS42*crLHS61 + crLHS44*crLHS62;
    const double crLHS64 = N[2]*crLHS7;
    const double crLHS65 = DN(2,1) + crLHS64;
    const double crLHS66 = crLHS28*crLHS65;
    const double crLHS67 = DN(0,0)*N[2];
    const double crLHS68 = N[3]*crLHS5;
    const double crLHS69 = DN(3,1)*crLHS8;
    const double crLHS70 = DN(3,0)*crLHS1;
    const double crLHS71 = crLHS70*rho;
    const double crLHS72 = DN(3,1)*crLHS2;
    const double crLHS73 = crLHS72*rho;
    const double crLHS74 = crLHS68 - crLHS69 + crLHS71 + crLHS73;
    const double crLHS75 = DN(0,0)*DN(3,0);
    const double crLHS76 = DN(0,1)*DN(3,1);
    const double crLHS77 = N[3]*crLHS6 + crLHS75*mu + crLHS76*mu;
    const double crLHS78 = crLHS4*crLHS75 + crLHS77;
    const double crLHS79 = DN(3,0)*N[0];
    const double crLHS80 = DN(3,1)*N[0];
    const double crLHS81 = crLHS42*crLHS79 + crLHS44*crLHS80;
    const double crLHS82 = N[3]*crLHS7;
    const double crLHS83 = DN(3,1) + crLHS82;
    const double crLHS84 = crLHS28*crLHS83;
    const double crLHS85 = DN(0,0)*N[3];
    const double crLHS86 = DN(0,1)*crLHS4;
    const double crLHS87 = crLHS26*crLHS4;
    const double crLHS88 = mu/pow(y, 2);
    const double crLHS89 = N[0]*crLHS88 + crLHS14;
    const double crLHS90 = DN(1,0)*w_gauss;
    const double crLHS91 = crLHS87*crLHS90;
    const double crLHS92 = crLHS25*crLHS4;
    const double crLHS93 = N[1]*crLHS88 + crLHS36;
    const double crLHS94 = DN(2,0)*w_gauss;
    const double crLHS95 = crLHS87*crLHS94;
    const double crLHS96 = N[2]*crLHS88 + crLHS56;
    const double crLHS97 = DN(3,0)*w_gauss;
    const double crLHS98 = crLHS87*crLHS97;
    const double crLHS99 = N[3]*crLHS88 + crLHS74;
    const double crLHS100 = DN(0,1)*crLHS17;
    const double crLHS101 = crLHS17*crLHS25;
    const double crLHS102 = crLHS17*w_gauss;
    const double crLHS103 = DN(0,0)*crLHS17;
    const double crLHS104 = N[1]*crLHS25;
    const double crLHS105 = crLHS37 + crLHS38;
    const double crLHS106 = N[2]*crLHS25;
    const double crLHS107 = crLHS57 + crLHS58;
    const double crLHS108 = N[3]*crLHS25;
    const double crLHS109 = crLHS75 + crLHS76;
    const double crLHS110 = crLHS18*(N[1]*crLHS15 + crLHS32);
    const double crLHS111 = crLHS18*(N[1]*crLHS20 + crLHS34);
    const double crLHS112 = N[1]*crLHS11 + N[1]*crLHS13;
    const double crLHS113 = pow(DN(1,0), 2);
    const double crLHS114 = pow(DN(1,1), 2);
    const double crLHS115 = pow(N[1], 2);
    const double crLHS116 = N[1]*crLHS33 + N[1]*crLHS35 + crLHS113*mu + crLHS114*mu + crLHS115*crLHS5;
    const double crLHS117 = crLHS4*crLHS90;
    const double crLHS118 = crLHS117*crLHS47;
    const double crLHS119 = DN(1,0)*DN(2,0);
    const double crLHS120 = DN(1,1)*DN(2,1);
    const double crLHS121 = N[2]*crLHS30 + crLHS119*mu + crLHS120*mu;
    const double crLHS122 = crLHS119*crLHS4 + crLHS121;
    const double crLHS123 = DN(2,0)*N[1];
    const double crLHS124 = DN(2,1)*N[1];
    const double crLHS125 = crLHS123*crLHS42 + crLHS124*crLHS44;
    const double crLHS126 = crLHS117*crLHS65;
    const double crLHS127 = DN(1,0)*N[2];
    const double crLHS128 = DN(1,0)*DN(3,0);
    const double crLHS129 = DN(1,1)*DN(3,1);
    const double crLHS130 = N[3]*crLHS30 + crLHS128*mu + crLHS129*mu;
    const double crLHS131 = crLHS128*crLHS4 + crLHS130;
    const double crLHS132 = DN(3,0)*N[1];
    const double crLHS133 = DN(3,1)*N[1];
    const double crLHS134 = crLHS132*crLHS42 + crLHS133*crLHS44;
    const double crLHS135 = crLHS117*crLHS83;
    const double crLHS136 = DN(1,0)*N[3];
    const double crLHS137 = DN(1,1)*crLHS4;
    const double crLHS138 = crLHS4*crLHS47;
    const double crLHS139 = crLHS138*crLHS94;
    const double crLHS140 = crLHS4*crLHS46;
    const double crLHS141 = crLHS138*crLHS97;
    const double crLHS142 = DN(1,0)*crLHS17;
    const double crLHS143 = DN(1,1)*crLHS17;
    const double crLHS144 = crLHS17*crLHS46;
    const double crLHS145 = N[2]*crLHS46;
    const double crLHS146 = crLHS119 + crLHS120;
    const double crLHS147 = N[3]*crLHS46;
    const double crLHS148 = crLHS128 + crLHS129;
    const double crLHS149 = crLHS18*(N[2]*crLHS15 + crLHS52);
    const double crLHS150 = crLHS18*(N[2]*crLHS20 + crLHS54);
    const double crLHS151 = N[2]*crLHS11 + N[2]*crLHS13;
    const double crLHS152 = N[2]*crLHS33 + N[2]*crLHS35;
    const double crLHS153 = pow(DN(2,0), 2);
    const double crLHS154 = pow(DN(2,1), 2);
    const double crLHS155 = pow(N[2], 2);
    const double crLHS156 = N[2]*crLHS53 + N[2]*crLHS55 + crLHS153*mu + crLHS154*mu + crLHS155*crLHS5;
    const double crLHS157 = crLHS4*crLHS94;
    const double crLHS158 = crLHS157*crLHS65;
    const double crLHS159 = DN(2,0)*DN(3,0);
    const double crLHS160 = DN(2,1)*DN(3,1);
    const double crLHS161 = N[3]*crLHS50 + crLHS159*mu + crLHS160*mu;
    const double crLHS162 = crLHS159*crLHS4 + crLHS161;
    const double crLHS163 = DN(3,0)*N[2];
    const double crLHS164 = DN(3,1)*N[2];
    const double crLHS165 = crLHS163*crLHS42 + crLHS164*crLHS44;
    const double crLHS166 = crLHS157*crLHS83;
    const double crLHS167 = DN(2,0)*N[3];
    const double crLHS168 = DN(2,1)*crLHS4;
    const double crLHS169 = crLHS4*crLHS65;
    const double crLHS170 = crLHS169*crLHS97;
    const double crLHS171 = crLHS4*crLHS83;
    const double crLHS172 = DN(2,0)*crLHS17;
    const double crLHS173 = DN(2,1)*crLHS17;
    const double crLHS174 = crLHS17*crLHS64;
    const double crLHS175 = N[3]*crLHS64;
    const double crLHS176 = crLHS159 + crLHS160;
    const double crLHS177 = crLHS18*(N[3]*crLHS15 + crLHS70);
    const double crLHS178 = crLHS18*(N[3]*crLHS20 + crLHS72);
    const double crLHS179 = N[3]*crLHS11 + N[3]*crLHS13;
    const double crLHS180 = N[3]*crLHS33 + N[3]*crLHS35;
    const double crLHS181 = N[3]*crLHS53 + N[3]*crLHS55;
    const double crLHS182 = pow(DN(3,0), 2);
    const double crLHS183 = pow(DN(3,1), 2);
    const double crLHS184 = pow(N[3], 2);
    const double crLHS185 = N[3]*crLHS71 + N[3]*crLHS73 + crLHS182*mu + crLHS183*mu + crLHS184*crLHS5;
    const double crLHS186 = crLHS171*crLHS97;
    const double crLHS187 = DN(3,1)*crLHS4;
    const double crLHS188 = DN(3,0)*crLHS17;
    const double crLHS189 = DN(3,1)*crLHS17;
    const double crLHS190 = crLHS17*crLHS82;
    rLHS(0,0)+=w_gauss*(crLHS0*crLHS4 + crLHS14*crLHS19 + crLHS14*crLHS21 + crLHS24);
    rLHS(0,1)+=crLHS29;
    rLHS(0,2)+=crLHS27*(-N[0] + crLHS19 + crLHS21);
    rLHS(0,3)+=w_gauss*(crLHS19*crLHS36 + crLHS21*crLHS36 + crLHS40 + crLHS45);
    rLHS(0,4)+=crLHS48;
    rLHS(0,5)+=w_gauss*(DN(1,0)*crLHS19 + DN(1,0)*crLHS21 - crLHS49);
    rLHS(0,6)+=w_gauss*(crLHS19*crLHS56 + crLHS21*crLHS56 + crLHS60 + crLHS63);
    rLHS(0,7)+=crLHS66;
    rLHS(0,8)+=w_gauss*(DN(2,0)*crLHS19 + DN(2,0)*crLHS21 - crLHS67);
    rLHS(0,9)+=w_gauss*(crLHS19*crLHS74 + crLHS21*crLHS74 + crLHS78 + crLHS81);
    rLHS(0,10)+=crLHS84;
    rLHS(0,11)+=w_gauss*(DN(3,0)*crLHS19 + DN(3,0)*crLHS21 - crLHS85);
    rLHS(1,0)+=crLHS29;
    rLHS(1,1)+=w_gauss*(crLHS19*crLHS89 + crLHS21*crLHS89 + crLHS24 + crLHS25*crLHS87 + crLHS26*crLHS86);
    rLHS(1,2)+=w_gauss*(DN(0,1)*crLHS19 + DN(0,1)*crLHS21 - N[0]*crLHS26);
    rLHS(1,3)+=crLHS91;
    rLHS(1,4)+=w_gauss*(crLHS19*crLHS93 + crLHS21*crLHS93 + crLHS39 + crLHS45 + crLHS47*crLHS86 + crLHS47*crLHS92);
    rLHS(1,5)+=w_gauss*(DN(1,1)*crLHS19 + DN(1,1)*crLHS21 - N[1]*crLHS26);
    rLHS(1,6)+=crLHS95;
    rLHS(1,7)+=w_gauss*(crLHS19*crLHS96 + crLHS21*crLHS96 + crLHS59 + crLHS63 + crLHS65*crLHS86 + crLHS65*crLHS92);
    rLHS(1,8)+=w_gauss*(DN(2,1)*crLHS19 + DN(2,1)*crLHS21 - N[2]*crLHS26);
    rLHS(1,9)+=crLHS98;
    rLHS(1,10)+=w_gauss*(crLHS19*crLHS99 + crLHS21*crLHS99 + crLHS77 + crLHS81 + crLHS83*crLHS86 + crLHS83*crLHS92);
    rLHS(1,11)+=w_gauss*(DN(3,1)*crLHS19 + DN(3,1)*crLHS21 - N[3]*crLHS26);
    rLHS(2,0)+=crLHS27*(N[0] + crLHS16*(1.0*crLHS11 + 1.0*crLHS13 + 1.0*crLHS6 - 1.0*crLHS9));
    rLHS(2,1)+=w_gauss*(DN(0,1)*N[0] + crLHS100*crLHS89 - crLHS101*crLHS89 + crLHS23*crLHS7);
    rLHS(2,2)+=crLHS102*(-DN(0,1)*crLHS25 + crLHS0 + crLHS22);
    rLHS(2,3)+=w_gauss*(crLHS103*crLHS36 + crLHS41);
    rLHS(2,4)+=w_gauss*(crLHS100*crLHS93 - crLHS101*crLHS93 + crLHS104 + crLHS43);
    rLHS(2,5)+=crLHS102*(-DN(1,1)*crLHS25 + crLHS105);
    rLHS(2,6)+=w_gauss*(crLHS103*crLHS56 + crLHS61);
    rLHS(2,7)+=w_gauss*(crLHS100*crLHS96 - crLHS101*crLHS96 + crLHS106 + crLHS62);
    rLHS(2,8)+=crLHS102*(-DN(2,1)*crLHS25 + crLHS107);
    rLHS(2,9)+=w_gauss*(crLHS103*crLHS74 + crLHS79);
    rLHS(2,10)+=w_gauss*(crLHS100*crLHS99 - crLHS101*crLHS99 + crLHS108 + crLHS80);
    rLHS(2,11)+=crLHS102*(-DN(3,1)*crLHS25 + crLHS109);
    rLHS(3,0)+=w_gauss*(crLHS110*crLHS14 + crLHS111*crLHS14 + crLHS112 + crLHS40);
    rLHS(3,1)+=crLHS91;
    rLHS(3,2)+=w_gauss*(DN(0,0)*crLHS110 + DN(0,0)*crLHS111 - crLHS41);
    rLHS(3,3)+=w_gauss*(crLHS110*crLHS36 + crLHS111*crLHS36 + crLHS113*crLHS4 + crLHS116);
    rLHS(3,4)+=crLHS118;
    rLHS(3,5)+=crLHS90*(-N[1] + crLHS110 + crLHS111);
    rLHS(3,6)+=w_gauss*(crLHS110*crLHS56 + crLHS111*crLHS56 + crLHS122 + crLHS125);
    rLHS(3,7)+=crLHS126;
    rLHS(3,8)+=w_gauss*(DN(2,0)*crLHS110 + DN(2,0)*crLHS111 - crLHS127);
    rLHS(3,9)+=w_gauss*(crLHS110*crLHS74 + crLHS111*crLHS74 + crLHS131 + crLHS134);
    rLHS(3,10)+=crLHS135;
    rLHS(3,11)+=w_gauss*(DN(3,0)*crLHS110 + DN(3,0)*crLHS111 - crLHS136);
    rLHS(4,0)+=crLHS48;
    rLHS(4,1)+=w_gauss*(crLHS110*crLHS89 + crLHS111*crLHS89 + crLHS112 + crLHS137*crLHS26 + crLHS39 + crLHS46*crLHS87);
    rLHS(4,2)+=w_gauss*(DN(0,1)*crLHS110 + DN(0,1)*crLHS111 - N[0]*crLHS47);
    rLHS(4,3)+=crLHS118;
    rLHS(4,4)+=w_gauss*(crLHS110*crLHS93 + crLHS111*crLHS93 + crLHS116 + crLHS137*crLHS47 + crLHS138*crLHS46);
    rLHS(4,5)+=w_gauss*(DN(1,1)*crLHS110 + DN(1,1)*crLHS111 - N[1]*crLHS47);
    rLHS(4,6)+=crLHS139;
    rLHS(4,7)+=w_gauss*(crLHS110*crLHS96 + crLHS111*crLHS96 + crLHS121 + crLHS125 + crLHS137*crLHS65 + crLHS140*crLHS65);
    rLHS(4,8)+=w_gauss*(DN(2,1)*crLHS110 + DN(2,1)*crLHS111 - N[2]*crLHS47);
    rLHS(4,9)+=crLHS141;
    rLHS(4,10)+=w_gauss*(crLHS110*crLHS99 + crLHS111*crLHS99 + crLHS130 + crLHS134 + crLHS137*crLHS83 + crLHS140*crLHS83);
    rLHS(4,11)+=w_gauss*(DN(3,1)*crLHS110 + DN(3,1)*crLHS111 - N[3]*crLHS47);
    rLHS(5,0)+=w_gauss*(crLHS14*crLHS142 + crLHS49);
    rLHS(5,1)+=w_gauss*(DN(0,1)*N[1] + crLHS104 + crLHS143*crLHS89 - crLHS144*crLHS89);
    rLHS(5,2)+=crLHS102*(-DN(0,1)*crLHS46 + crLHS105);
    rLHS(5,3)+=crLHS90*(N[1] + crLHS16*(1.0*crLHS30 - 1.0*crLHS31 + 1.0*crLHS33 + 1.0*crLHS35));
    rLHS(5,4)+=w_gauss*(DN(1,1)*N[1] + crLHS115*crLHS7 + crLHS143*crLHS93 - crLHS144*crLHS93);
    rLHS(5,5)+=crLHS102*(-DN(1,1)*crLHS46 + crLHS113 + crLHS114);
    rLHS(5,6)+=w_gauss*(crLHS123 + crLHS142*crLHS56);
    rLHS(5,7)+=w_gauss*(crLHS124 + crLHS143*crLHS96 - crLHS144*crLHS96 + crLHS145);
    rLHS(5,8)+=crLHS102*(-DN(2,1)*crLHS46 + crLHS146);
    rLHS(5,9)+=w_gauss*(crLHS132 + crLHS142*crLHS74);
    rLHS(5,10)+=w_gauss*(crLHS133 + crLHS143*crLHS99 - crLHS144*crLHS99 + crLHS147);
    rLHS(5,11)+=crLHS102*(-DN(3,1)*crLHS46 + crLHS148);
    rLHS(6,0)+=w_gauss*(crLHS14*crLHS149 + crLHS14*crLHS150 + crLHS151 + crLHS60);
    rLHS(6,1)+=crLHS95;
    rLHS(6,2)+=w_gauss*(DN(0,0)*crLHS149 + DN(0,0)*crLHS150 - crLHS61);
    rLHS(6,3)+=w_gauss*(crLHS122 + crLHS149*crLHS36 + crLHS150*crLHS36 + crLHS152);
    rLHS(6,4)+=crLHS139;
    rLHS(6,5)+=w_gauss*(DN(1,0)*crLHS149 + DN(1,0)*crLHS150 - crLHS123);
    rLHS(6,6)+=w_gauss*(crLHS149*crLHS56 + crLHS150*crLHS56 + crLHS153*crLHS4 + crLHS156);
    rLHS(6,7)+=crLHS158;
    rLHS(6,8)+=crLHS94*(-N[2] + crLHS149 + crLHS150);
    rLHS(6,9)+=w_gauss*(crLHS149*crLHS74 + crLHS150*crLHS74 + crLHS162 + crLHS165);
    rLHS(6,10)+=crLHS166;
    rLHS(6,11)+=w_gauss*(DN(3,0)*crLHS149 + DN(3,0)*crLHS150 - crLHS167);
    rLHS(7,0)+=crLHS66;
    rLHS(7,1)+=w_gauss*(crLHS149*crLHS89 + crLHS150*crLHS89 + crLHS151 + crLHS168*crLHS26 + crLHS59 + crLHS64*crLHS87);
    rLHS(7,2)+=w_gauss*(DN(0,1)*crLHS149 + DN(0,1)*crLHS150 - N[0]*crLHS65);
    rLHS(7,3)+=crLHS126;
    rLHS(7,4)+=w_gauss*(crLHS121 + crLHS138*crLHS64 + crLHS149*crLHS93 + crLHS150*crLHS93 + crLHS152 + crLHS168*crLHS47);
    rLHS(7,5)+=w_gauss*(DN(1,1)*crLHS149 + DN(1,1)*crLHS150 - N[1]*crLHS65);
    rLHS(7,6)+=crLHS158;
    rLHS(7,7)+=w_gauss*(crLHS149*crLHS96 + crLHS150*crLHS96 + crLHS156 + crLHS168*crLHS65 + crLHS169*crLHS64);
    rLHS(7,8)+=w_gauss*(DN(2,1)*crLHS149 + DN(2,1)*crLHS150 - N[2]*crLHS65);
    rLHS(7,9)+=crLHS170;
    rLHS(7,10)+=w_gauss*(crLHS149*crLHS99 + crLHS150*crLHS99 + crLHS161 + crLHS165 + crLHS168*crLHS83 + crLHS171*crLHS64);
    rLHS(7,11)+=w_gauss*(DN(3,1)*crLHS149 + DN(3,1)*crLHS150 - N[3]*crLHS65);
    rLHS(8,0)+=w_gauss*(crLHS14*crLHS172 + crLHS67);
    rLHS(8,1)+=w_gauss*(DN(0,1)*N[2] + crLHS106 + crLHS173*crLHS89 - crLHS174*crLHS89);
    rLHS(8,2)+=crLHS102*(-DN(0,1)*crLHS64 + crLHS107);
    rLHS(8,3)+=w_gauss*(crLHS127 + crLHS172*crLHS36);
    rLHS(8,4)+=w_gauss*(DN(1,1)*N[2] + crLHS145 + crLHS173*crLHS93 - crLHS174*crLHS93);
    rLHS(8,5)+=crLHS102*(-DN(1,1)*crLHS64 + crLHS146);
    rLHS(8,6)+=crLHS94*(N[2] + crLHS16*(1.0*crLHS50 - 1.0*crLHS51 + 1.0*crLHS53 + 1.0*crLHS55));
    rLHS(8,7)+=w_gauss*(DN(2,1)*N[2] + crLHS155*crLHS7 + crLHS173*crLHS96 - crLHS174*crLHS96);
    rLHS(8,8)+=crLHS102*(-DN(2,1)*crLHS64 + crLHS153 + crLHS154);
    rLHS(8,9)+=w_gauss*(crLHS163 + crLHS172*crLHS74);
    rLHS(8,10)+=w_gauss*(crLHS164 + crLHS173*crLHS99 - crLHS174*crLHS99 + crLHS175);
    rLHS(8,11)+=crLHS102*(-DN(3,1)*crLHS64 + crLHS176);
    rLHS(9,0)+=w_gauss*(crLHS14*crLHS177 + crLHS14*crLHS178 + crLHS179 + crLHS78);
    rLHS(9,1)+=crLHS98;
    rLHS(9,2)+=w_gauss*(DN(0,0)*crLHS177 + DN(0,0)*crLHS178 - crLHS79);
    rLHS(9,3)+=w_gauss*(crLHS131 + crLHS177*crLHS36 + crLHS178*crLHS36 + crLHS180);
    rLHS(9,4)+=crLHS141;
    rLHS(9,5)+=w_gauss*(DN(1,0)*crLHS177 + DN(1,0)*crLHS178 - crLHS132);
    rLHS(9,6)+=w_gauss*(crLHS162 + crLHS177*crLHS56 + crLHS178*crLHS56 + crLHS181);
    rLHS(9,7)+=crLHS170;
    rLHS(9,8)+=w_gauss*(DN(2,0)*crLHS177 + DN(2,0)*crLHS178 - crLHS163);
    rLHS(9,9)+=w_gauss*(crLHS177*crLHS74 + crLHS178*crLHS74 + crLHS182*crLHS4 + crLHS185);
    rLHS(9,10)+=crLHS186;
    rLHS(9,11)+=crLHS97*(-N[3] + crLHS177 + crLHS178);
    rLHS(10,0)+=crLHS84;
    rLHS(10,1)+=w_gauss*(crLHS177*crLHS89 + crLHS178*crLHS89 + crLHS179 + crLHS187*crLHS26 + crLHS77 + crLHS82*crLHS87);
    rLHS(10,2)+=w_gauss*(DN(0,1)*crLHS177 + DN(0,1)*crLHS178 - N[0]*crLHS83);
    rLHS(10,3)+=crLHS135;
    rLHS(10,4)+=w_gauss*(crLHS130 + crLHS138*crLHS82 + crLHS177*crLHS93 + crLHS178*crLHS93 + crLHS180 + crLHS187*crLHS47);
    rLHS(10,5)+=w_gauss*(DN(1,1)*crLHS177 + DN(1,1)*crLHS178 - N[1]*crLHS83);
    rLHS(10,6)+=crLHS166;
    rLHS(10,7)+=w_gauss*(crLHS161 + crLHS169*crLHS82 + crLHS177*crLHS96 + crLHS178*crLHS96 + crLHS181 + crLHS187*crLHS65);
    rLHS(10,8)+=w_gauss*(DN(2,1)*crLHS177 + DN(2,1)*crLHS178 - N[2]*crLHS83);
    rLHS(10,9)+=crLHS186;
    rLHS(10,10)+=w_gauss*(crLHS171*crLHS82 + crLHS177*crLHS99 + crLHS178*crLHS99 + crLHS185 + crLHS187*crLHS83);
    rLHS(10,11)+=w_gauss*(DN(3,1)*crLHS177 + DN(3,1)*crLHS178 - N[3]*crLHS83);
    rLHS(11,0)+=w_gauss*(crLHS14*crLHS188 + crLHS85);
    rLHS(11,1)+=w_gauss*(DN(0,1)*N[3] + crLHS108 + crLHS189*crLHS89 - crLHS190*crLHS89);
    rLHS(11,2)+=crLHS102*(-DN(0,1)*crLHS82 + crLHS109);
    rLHS(11,3)+=w_gauss*(crLHS136 + crLHS188*crLHS36);
    rLHS(11,4)+=w_gauss*(DN(1,1)*N[3] + crLHS147 + crLHS189*crLHS93 - crLHS190*crLHS93);
    rLHS(11,5)+=crLHS102*(-DN(1,1)*crLHS82 + crLHS148);
    rLHS(11,6)+=w_gauss*(crLHS167 + crLHS188*crLHS56);
    rLHS(11,7)+=w_gauss*(DN(2,1)*N[3] + crLHS175 + crLHS189*crLHS96 - crLHS190*crLHS96);
    rLHS(11,8)+=crLHS102*(-DN(2,1)*crLHS82 + crLHS176);
    rLHS(11,9)+=crLHS97*(N[3] + crLHS16*(1.0*crLHS68 - 1.0*crLHS69 + 1.0*crLHS71 + 1.0*crLHS73));
    rLHS(11,10)+=w_gauss*(DN(3,1)*N[3] + crLHS184*crLHS7 + crLHS189*crLHS99 - crLHS190*crLHS99);
    rLHS(11,11)+=crLHS102*(-DN(3,1)*crLHS82 + crLHS182 + crLHS183);

}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    AxisymmetricNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    // Axisymmetric formulation radius
    const double y = this->CalculateGaussPointRadius(rData);

    // Material parameters
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const auto& p = rData.Pressure;
    const auto& v = rData.Velocity;
    const auto& v_n = rData.Velocity_OldStep1;
    const auto& v_nn = rData.Velocity_OldStep2;
    const auto& v_mesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> v_conv = v - v_mesh;
    const auto& f = rData.BodyForce;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add RHS Gauss point contribution
    const double w_gauss = 2.0 * Globals::Pi * y * rData.Weight;

    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
    const double crRHS1 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
    const double crRHS2 = DN(0,0)*mu;
    const double crRHS3 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
    const double crRHS4 = crRHS3*mu;
    const double crRHS5 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
    const double crRHS6 = N[0]*crRHS1;
    const double crRHS7 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
    const double crRHS8 = crRHS7*rho;
    const double crRHS9 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
    const double crRHS10 = crRHS9*rho;
    const double crRHS11 = crRHS10*crRHS3;
    const double crRHS12 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
    const double crRHS13 = rho*stab_c2*sqrt(pow(crRHS7, 2) + pow(crRHS9, 2));
    const double crRHS14 = 1.0/y;
    const double crRHS15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
    const double crRHS16 = crRHS14*crRHS15;
    const double crRHS17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double crRHS18 = (crRHS13*h/stab_c1 + mu)*(crRHS1 + crRHS16 + crRHS17);
    const double crRHS19 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
    const double crRHS20 = DN(0,0)*crRHS7 + N[0]*crRHS19;
    const double crRHS21 = 1.0/(crRHS13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crRHS22 = crRHS21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crRHS1*crRHS8 + crRHS11 + crRHS12 - crRHS14*crRHS4 - crRHS5);
    const double crRHS23 = crRHS22*rho;
    const double crRHS24 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
    const double crRHS25 = DN(0,1)*crRHS9 + N[0]*crRHS24;
    const double crRHS26 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
    const double crRHS27 = crRHS17*mu;
    const double crRHS28 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
    const double crRHS29 = N[0]*crRHS14;
    const double crRHS30 = crRHS26*crRHS8;
    const double crRHS31 = N[0]*crRHS17;
    const double crRHS32 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
    const double crRHS33 = crRHS21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crRHS10*crRHS17 - crRHS14*crRHS27 + crRHS15*mu/pow(y, 2) - crRHS28 + crRHS30 + crRHS32);
    const double crRHS34 = crRHS33*rho;
    const double crRHS35 = DN(1,0)*mu;
    const double crRHS36 = N[1]*crRHS1;
    const double crRHS37 = DN(1,0)*crRHS7 + N[1]*crRHS19;
    const double crRHS38 = DN(1,1)*crRHS9 + N[1]*crRHS24;
    const double crRHS39 = N[1]*crRHS14;
    const double crRHS40 = N[1]*crRHS17;
    const double crRHS41 = DN(2,0)*mu;
    const double crRHS42 = N[2]*crRHS1;
    const double crRHS43 = DN(2,0)*crRHS7 + N[2]*crRHS19;
    const double crRHS44 = DN(2,1)*crRHS9 + N[2]*crRHS24;
    const double crRHS45 = N[2]*crRHS14;
    const double crRHS46 = N[2]*crRHS17;
    rRHS[0]+=-w_gauss*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS18 + DN(0,1)*crRHS4 + N[0]*crRHS11 + N[0]*crRHS12 - N[0]*crRHS5 + crRHS1*crRHS2 + crRHS20*crRHS23 + crRHS23*crRHS25 + crRHS6*crRHS8);
    rRHS[1]+=-w_gauss*(DN(0,1)*crRHS18 + DN(0,1)*crRHS27 - N[0]*crRHS28 + N[0]*crRHS30 + N[0]*crRHS32 - crRHS0*(DN(0,1) + crRHS29) + crRHS10*crRHS31 + crRHS18*crRHS29 + crRHS2*crRHS26 + crRHS20*crRHS34 + crRHS25*crRHS34);
    rRHS[2]+=-w_gauss*(DN(0,0)*crRHS22 + DN(0,1)*crRHS33 + N[0]*crRHS16 - crRHS29*crRHS33 + crRHS31 + crRHS6);
    rRHS[3]+=-w_gauss*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS18 + DN(1,1)*crRHS4 + N[1]*crRHS11 + N[1]*crRHS12 - N[1]*crRHS5 + crRHS1*crRHS35 + crRHS23*crRHS37 + crRHS23*crRHS38 + crRHS36*crRHS8);
    rRHS[4]+=-w_gauss*(DN(1,1)*crRHS18 + DN(1,1)*crRHS27 - N[1]*crRHS28 + N[1]*crRHS30 + N[1]*crRHS32 - crRHS0*(DN(1,1) + crRHS39) + crRHS10*crRHS40 + crRHS18*crRHS39 + crRHS26*crRHS35 + crRHS34*crRHS37 + crRHS34*crRHS38);
    rRHS[5]+=-w_gauss*(DN(1,0)*crRHS22 + DN(1,1)*crRHS33 + N[1]*crRHS16 - crRHS33*crRHS39 + crRHS36 + crRHS40);
    rRHS[6]+=-w_gauss*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS18 + DN(2,1)*crRHS4 + N[2]*crRHS11 + N[2]*crRHS12 - N[2]*crRHS5 + crRHS1*crRHS41 + crRHS23*crRHS43 + crRHS23*crRHS44 + crRHS42*crRHS8);
    rRHS[7]+=-w_gauss*(DN(2,1)*crRHS18 + DN(2,1)*crRHS27 - N[2]*crRHS28 + N[2]*crRHS30 + N[2]*crRHS32 - crRHS0*(DN(2,1) + crRHS45) + crRHS10*crRHS46 + crRHS18*crRHS45 + crRHS26*crRHS41 + crRHS34*crRHS43 + crRHS34*crRHS44);
    rRHS[8]+=-w_gauss*(DN(2,0)*crRHS22 + DN(2,1)*crRHS33 + N[2]*crRHS16 - crRHS33*crRHS45 + crRHS42 + crRHS46);

}

template <>
void AxisymmetricNavierStokes<AxisymmetricNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    AxisymmetricNavierStokesData<2,4>& rData,
    VectorType& rRHS)
{
    // Axisymmetric formulation radius
    const double y = this->CalculateGaussPointRadius(rData);

    // Material parameters
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Dynamic parameters
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    // Nodal data
    const auto& p = rData.Pressure;
    const auto& v = rData.Velocity;
    const auto& v_n = rData.Velocity_OldStep1;
    const auto& v_nn = rData.Velocity_OldStep2;
    const auto& v_mesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,4> v_conv = v - v_mesh;
    const auto& f = rData.BodyForce;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const double h = rData.ElementSize;
    const double dyn_tau = rData.DynamicTau;

    // Add RHS Gauss point contribution
    const double w_gauss = 2.0 * Globals::Pi * y * rData.Weight;

    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
    const double crRHS1 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
    const double crRHS2 = DN(0,0)*mu;
    const double crRHS3 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
    const double crRHS4 = crRHS3*mu;
    const double crRHS5 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
    const double crRHS6 = N[0]*crRHS1;
    const double crRHS7 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
    const double crRHS8 = crRHS7*rho;
    const double crRHS9 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
    const double crRHS10 = crRHS9*rho;
    const double crRHS11 = crRHS10*crRHS3;
    const double crRHS12 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
    const double crRHS13 = rho*stab_c2*sqrt(pow(crRHS7, 2) + pow(crRHS9, 2));
    const double crRHS14 = 1.0/y;
    const double crRHS15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
    const double crRHS16 = crRHS14*crRHS15;
    const double crRHS17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
    const double crRHS18 = (crRHS13*h/stab_c1 + mu)*(crRHS1 + crRHS16 + crRHS17);
    const double crRHS19 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
    const double crRHS20 = DN(0,0)*crRHS7 + N[0]*crRHS19;
    const double crRHS21 = 1.0/(crRHS13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crRHS22 = crRHS21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crRHS1*crRHS8 + crRHS11 + crRHS12 - crRHS14*crRHS4 - crRHS5);
    const double crRHS23 = crRHS22*rho;
    const double crRHS24 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
    const double crRHS25 = DN(0,1)*crRHS9 + N[0]*crRHS24;
    const double crRHS26 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
    const double crRHS27 = crRHS17*mu;
    const double crRHS28 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
    const double crRHS29 = N[0]*crRHS14;
    const double crRHS30 = crRHS26*crRHS8;
    const double crRHS31 = N[0]*crRHS17;
    const double crRHS32 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
    const double crRHS33 = crRHS21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crRHS10*crRHS17 - crRHS14*crRHS27 + crRHS15*mu/pow(y, 2) - crRHS28 + crRHS30 + crRHS32);
    const double crRHS34 = crRHS33*rho;
    const double crRHS35 = DN(1,0)*mu;
    const double crRHS36 = N[1]*crRHS1;
    const double crRHS37 = DN(1,0)*crRHS7 + N[1]*crRHS19;
    const double crRHS38 = DN(1,1)*crRHS9 + N[1]*crRHS24;
    const double crRHS39 = N[1]*crRHS14;
    const double crRHS40 = N[1]*crRHS17;
    const double crRHS41 = DN(2,0)*mu;
    const double crRHS42 = N[2]*crRHS1;
    const double crRHS43 = DN(2,0)*crRHS7 + N[2]*crRHS19;
    const double crRHS44 = DN(2,1)*crRHS9 + N[2]*crRHS24;
    const double crRHS45 = N[2]*crRHS14;
    const double crRHS46 = N[2]*crRHS17;
    const double crRHS47 = DN(3,0)*mu;
    const double crRHS48 = N[3]*crRHS1;
    const double crRHS49 = DN(3,0)*crRHS7 + N[3]*crRHS19;
    const double crRHS50 = DN(3,1)*crRHS9 + N[3]*crRHS24;
    const double crRHS51 = N[3]*crRHS14;
    const double crRHS52 = N[3]*crRHS17;
    rRHS[0]+=-w_gauss*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS18 + DN(0,1)*crRHS4 + N[0]*crRHS11 + N[0]*crRHS12 - N[0]*crRHS5 + crRHS1*crRHS2 + crRHS20*crRHS23 + crRHS23*crRHS25 + crRHS6*crRHS8);
    rRHS[1]+=-w_gauss*(DN(0,1)*crRHS18 + DN(0,1)*crRHS27 - N[0]*crRHS28 + N[0]*crRHS30 + N[0]*crRHS32 - crRHS0*(DN(0,1) + crRHS29) + crRHS10*crRHS31 + crRHS18*crRHS29 + crRHS2*crRHS26 + crRHS20*crRHS34 + crRHS25*crRHS34);
    rRHS[2]+=-w_gauss*(DN(0,0)*crRHS22 + DN(0,1)*crRHS33 + N[0]*crRHS16 - crRHS29*crRHS33 + crRHS31 + crRHS6);
    rRHS[3]+=-w_gauss*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS18 + DN(1,1)*crRHS4 + N[1]*crRHS11 + N[1]*crRHS12 - N[1]*crRHS5 + crRHS1*crRHS35 + crRHS23*crRHS37 + crRHS23*crRHS38 + crRHS36*crRHS8);
    rRHS[4]+=-w_gauss*(DN(1,1)*crRHS18 + DN(1,1)*crRHS27 - N[1]*crRHS28 + N[1]*crRHS30 + N[1]*crRHS32 - crRHS0*(DN(1,1) + crRHS39) + crRHS10*crRHS40 + crRHS18*crRHS39 + crRHS26*crRHS35 + crRHS34*crRHS37 + crRHS34*crRHS38);
    rRHS[5]+=-w_gauss*(DN(1,0)*crRHS22 + DN(1,1)*crRHS33 + N[1]*crRHS16 - crRHS33*crRHS39 + crRHS36 + crRHS40);
    rRHS[6]+=-w_gauss*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS18 + DN(2,1)*crRHS4 + N[2]*crRHS11 + N[2]*crRHS12 - N[2]*crRHS5 + crRHS1*crRHS41 + crRHS23*crRHS43 + crRHS23*crRHS44 + crRHS42*crRHS8);
    rRHS[7]+=-w_gauss*(DN(2,1)*crRHS18 + DN(2,1)*crRHS27 - N[2]*crRHS28 + N[2]*crRHS30 + N[2]*crRHS32 - crRHS0*(DN(2,1) + crRHS45) + crRHS10*crRHS46 + crRHS18*crRHS45 + crRHS26*crRHS41 + crRHS34*crRHS43 + crRHS34*crRHS44);
    rRHS[8]+=-w_gauss*(DN(2,0)*crRHS22 + DN(2,1)*crRHS33 + N[2]*crRHS16 - crRHS33*crRHS45 + crRHS42 + crRHS46);
    rRHS[9]+=-w_gauss*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS18 + DN(3,1)*crRHS4 + N[3]*crRHS11 + N[3]*crRHS12 - N[3]*crRHS5 + crRHS1*crRHS47 + crRHS23*crRHS49 + crRHS23*crRHS50 + crRHS48*crRHS8);
    rRHS[10]+=-w_gauss*(DN(3,1)*crRHS18 + DN(3,1)*crRHS27 - N[3]*crRHS28 + N[3]*crRHS30 + N[3]*crRHS32 - crRHS0*(DN(3,1) + crRHS51) + crRHS10*crRHS52 + crRHS18*crRHS51 + crRHS26*crRHS47 + crRHS34*crRHS49 + crRHS34*crRHS50);
    rRHS[11]+=-w_gauss*(DN(3,0)*crRHS22 + DN(3,1)*crRHS33 + N[3]*crRHS16 - crRHS33*crRHS51 + crRHS48 + crRHS52);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template< class TElementData >
double AxisymmetricNavierStokes<TElementData>::CalculateGaussPointRadius(const TElementData& rData) const
{
    double radius = 0.0;
    const auto& r_N = rData.N;
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < NumNodes; ++i) {
        radius += r_N[i] * r_geom[i].Y();
    }
    return radius;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void AxisymmetricNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}


template< class TElementData >
void AxisymmetricNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,3> >;
template class AxisymmetricNavierStokes< AxisymmetricNavierStokesData<2,4> >;

}