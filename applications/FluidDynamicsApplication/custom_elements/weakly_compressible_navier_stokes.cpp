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
#include "weakly_compressible_navier_stokes.h"
#include "custom_utilities/weakly_compressible_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::~WeaklyCompressibleNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokes>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int WeaklyCompressibleNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template<class TElementData>
const Parameters WeaklyCompressibleNavierStokes<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["VORTICITY","Q_VALUE","VORTICITY_MAGNITUDE"],
            "nodal_historical"       : ["VELOCITY","PRESSURE","DENSITY","DYNAMIC_VISCOSITY"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","IS_STRUCTURE","DISPLACEMENT","BODY_FORCE","NODAL_AREA","NODAL_H","ADVPROJ","DIVPROJ","REACTION","REACTION_WATER_PRESSURE","EXTERNAL_PRESSURE","NORMAL","Y_WALL","Q_VALUE"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["Newtonian2DLaw","Newtonian3DLaw","NewtonianTemperatureDependent2DLaw","NewtonianTemperatureDependent3DLaw","Euler2DLaw","Euler3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This implements a weakly compressible Navier-Stokes element with quasi-static Variational MultiScales (VMS) stabilization. Note that this formulation allows a weak coupling with a custom Equation of State that updates the nodal DENSITY from the obtained PRESSURE values. Also note that no viscous behavior is hardcoded, meaning that any fluid constitutive model can be used through a constitutive law."
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TElementData>
std::string WeaklyCompressibleNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "WeaklyCompressibleNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddBoundaryTraction(
    TElementData& rData,
    const Vector& rUnitNormal,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    // Set the current Gauss pt. Voigt notation normal projection matrix
    BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
    FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal, voigt_normal_projection_matrix);

    // Set the current Gauss pt. strain matrix
    BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX, B_matrix);

    // Compute some Gauss pt. auxiliar matrices
    const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
    const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

    // Fill the pressure to Voigt notation operator matrix
    BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            pres_to_voigt_matrix_op(comp, i*BlockSize+Dim) = rData.N[i];
        }
    }

    // Set the shape functions auxiliar transpose matrix
    BoundedMatrix<double, LocalSize, Dim> N_aux_trans = ZeroMatrix(LocalSize, Dim);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            N_aux_trans(i*BlockSize+comp, comp) = rData.N[i];
        }
    }

    // Contribution coming fron the shear stress operator
    noalias(rData.lhs) = prod(N_aux_trans, aux_matrix_ACB);

    // Contribution coming from the pressure terms
    const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
    noalias(rData.lhs) -= prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);

    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    rData.lhs *= rData.Weight;
    noalias(rLHS) -= rData.lhs;
    noalias(rRHS) += prod(rData.lhs,values);
}

template <>
void WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> >::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double sigma = rData.Resistance;

    const double h = rData.ElementSize;
    const array_1d<double,3>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,3,3>& C = rData.C;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = pow(DN(0,0), 2);
const double crLHS4 = sigma*stab_c3;
const double crLHS5 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crLHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS8 = crLHS5*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2));
const double crLHS9 = mu + (crLHS4 + crLHS8*h)/stab_c1;
const double crLHS10 = pow(N[0], 2);
const double crLHS11 = DN(0,0)*crLHS6 + DN(0,1)*crLHS7;
const double crLHS12 = N[0]*crLHS5;
const double crLHS13 = bdf0*crLHS5;
const double crLHS14 = N[0]*bdf0;
const double crLHS15 = crLHS11 + crLHS14;
const double crLHS16 = pow(h, -2);
const double crLHS17 = 1.0/(crLHS16*crLHS4 + crLHS16*mu*stab_c1 + crLHS5*dyn_tau/dt + crLHS8/h);
const double crLHS18 = crLHS17*pow(crLHS5, 2);
const double crLHS19 = crLHS11*crLHS18;
const double crLHS20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS21 = crLHS18*crLHS20;
const double crLHS22 = N[0]*crLHS21;
const double crLHS23 = crLHS17*sigma;
const double crLHS24 = crLHS12*crLHS23;
const double crLHS25 = crLHS10*crLHS13 + crLHS10*sigma + crLHS11*crLHS12 + crLHS15*crLHS19 + crLHS15*crLHS22 - crLHS15*crLHS24;
const double crLHS26 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS27 = C(1,2)*DN(0,1);
const double crLHS28 = C(2,2)*DN(0,0) + crLHS27;
const double crLHS29 = DN(0,0)*crLHS9;
const double crLHS30 = DN(0,1)*crLHS29;
const double crLHS31 = N[0]*sigma;
const double crLHS32 = 1/(crLHS5*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2));
const double crLHS33 = crLHS14*crLHS32;
const double crLHS34 = crLHS33*crLHS9;
const double crLHS35 = crLHS17*crLHS20;
const double crLHS36 = crLHS17*crLHS5;
const double crLHS37 = crLHS11*crLHS36;
const double crLHS38 = gauss_weight*(-N[0] + crLHS12*crLHS35 - crLHS17*crLHS31 + crLHS34 + crLHS37);
const double crLHS39 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS40 = C(0,2)*DN(1,0);
const double crLHS41 = C(2,2)*DN(1,1) + crLHS40;
const double crLHS42 = N[1]*crLHS5;
const double crLHS43 = N[1]*crLHS31 + crLHS14*crLHS42;
const double crLHS44 = DN(1,0)*crLHS29 + crLHS43;
const double crLHS45 = DN(1,0)*crLHS6 + DN(1,1)*crLHS7;
const double crLHS46 = N[1]*bdf0;
const double crLHS47 = crLHS45 + crLHS46;
const double crLHS48 = crLHS12*crLHS45 + crLHS19*crLHS47 + crLHS22*crLHS47 - crLHS24*crLHS47;
const double crLHS49 = C(0,1)*DN(1,1) + crLHS40;
const double crLHS50 = C(1,2)*DN(1,1);
const double crLHS51 = C(2,2)*DN(1,0) + crLHS50;
const double crLHS52 = DN(1,1)*crLHS29;
const double crLHS53 = DN(0,0)*N[1];
const double crLHS54 = DN(1,0)*N[0];
const double crLHS55 = crLHS32*crLHS46;
const double crLHS56 = crLHS20*crLHS36;
const double crLHS57 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS58 = C(0,2)*DN(2,0);
const double crLHS59 = C(2,2)*DN(2,1) + crLHS58;
const double crLHS60 = N[2]*crLHS5;
const double crLHS61 = N[2]*crLHS31 + crLHS14*crLHS60;
const double crLHS62 = DN(2,0)*crLHS29 + crLHS61;
const double crLHS63 = DN(2,0)*crLHS6 + DN(2,1)*crLHS7;
const double crLHS64 = N[2]*bdf0;
const double crLHS65 = crLHS63 + crLHS64;
const double crLHS66 = crLHS12*crLHS63 + crLHS19*crLHS65 + crLHS22*crLHS65 - crLHS24*crLHS65;
const double crLHS67 = C(0,1)*DN(2,1) + crLHS58;
const double crLHS68 = C(1,2)*DN(2,1);
const double crLHS69 = C(2,2)*DN(2,0) + crLHS68;
const double crLHS70 = DN(2,1)*crLHS29;
const double crLHS71 = DN(0,0)*N[2];
const double crLHS72 = DN(2,0)*N[0];
const double crLHS73 = crLHS32*crLHS64;
const double crLHS74 = C(0,1)*DN(0,0) + crLHS27;
const double crLHS75 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS76 = pow(DN(0,1), 2);
const double crLHS77 = C(0,1)*DN(1,0) + crLHS50;
const double crLHS78 = DN(0,1)*crLHS9;
const double crLHS79 = DN(1,0)*crLHS78;
const double crLHS80 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS81 = DN(1,1)*crLHS78 + crLHS43;
const double crLHS82 = DN(0,1)*N[1];
const double crLHS83 = DN(1,1)*N[0];
const double crLHS84 = C(0,1)*DN(2,0) + crLHS68;
const double crLHS85 = DN(2,0)*crLHS78;
const double crLHS86 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS87 = DN(2,1)*crLHS78 + crLHS61;
const double crLHS88 = DN(0,1)*N[2];
const double crLHS89 = DN(2,1)*N[0];
const double crLHS90 = crLHS15*crLHS36;
const double crLHS91 = gauss_weight*(N[0] + crLHS90);
const double crLHS92 = bdf0*crLHS32;
const double crLHS93 = crLHS36*crLHS47;
const double crLHS94 = DN(0,0)*crLHS17;
const double crLHS95 = DN(0,1)*crLHS17;
const double crLHS96 = gauss_weight*(DN(1,0)*crLHS94 + DN(1,1)*crLHS95 + N[1]*crLHS33);
const double crLHS97 = crLHS36*crLHS65;
const double crLHS98 = gauss_weight*(DN(2,0)*crLHS94 + DN(2,1)*crLHS95 + N[2]*crLHS33);
const double crLHS99 = crLHS18*crLHS45;
const double crLHS100 = N[1]*crLHS21;
const double crLHS101 = crLHS23*crLHS42;
const double crLHS102 = crLHS100*crLHS15 - crLHS101*crLHS15 + crLHS11*crLHS42 + crLHS15*crLHS99;
const double crLHS103 = DN(1,0)*crLHS9;
const double crLHS104 = crLHS36*crLHS45;
const double crLHS105 = pow(DN(1,0), 2);
const double crLHS106 = pow(N[1], 2);
const double crLHS107 = crLHS100*crLHS47 - crLHS101*crLHS47 + crLHS106*crLHS13 + crLHS106*sigma + crLHS42*crLHS45 + crLHS47*crLHS99;
const double crLHS108 = DN(1,1)*crLHS103;
const double crLHS109 = N[1]*sigma;
const double crLHS110 = crLHS55*crLHS9;
const double crLHS111 = gauss_weight*(-N[1] + crLHS104 - crLHS109*crLHS17 + crLHS110 + crLHS35*crLHS42);
const double crLHS112 = N[2]*crLHS109 + crLHS46*crLHS60;
const double crLHS113 = DN(2,0)*crLHS103 + crLHS112;
const double crLHS114 = crLHS100*crLHS65 - crLHS101*crLHS65 + crLHS42*crLHS63 + crLHS65*crLHS99;
const double crLHS115 = DN(2,1)*crLHS103;
const double crLHS116 = DN(1,0)*N[2];
const double crLHS117 = DN(2,0)*N[1];
const double crLHS118 = DN(1,1)*crLHS9;
const double crLHS119 = pow(DN(1,1), 2);
const double crLHS120 = DN(2,0)*crLHS118;
const double crLHS121 = DN(2,1)*crLHS118 + crLHS112;
const double crLHS122 = DN(1,1)*N[2];
const double crLHS123 = DN(2,1)*N[1];
const double crLHS124 = gauss_weight*(N[1] + crLHS93);
const double crLHS125 = gauss_weight*(DN(1,0)*DN(2,0)*crLHS17 + DN(1,1)*DN(2,1)*crLHS17 + N[2]*crLHS55);
const double crLHS126 = crLHS18*crLHS63;
const double crLHS127 = N[2]*crLHS21;
const double crLHS128 = crLHS23*crLHS60;
const double crLHS129 = crLHS11*crLHS60 + crLHS126*crLHS15 + crLHS127*crLHS15 - crLHS128*crLHS15;
const double crLHS130 = DN(2,0)*crLHS9;
const double crLHS131 = crLHS36*crLHS63;
const double crLHS132 = crLHS126*crLHS47 + crLHS127*crLHS47 - crLHS128*crLHS47 + crLHS45*crLHS60;
const double crLHS133 = pow(DN(2,0), 2);
const double crLHS134 = pow(N[2], 2);
const double crLHS135 = crLHS126*crLHS65 + crLHS127*crLHS65 - crLHS128*crLHS65 + crLHS13*crLHS134 + crLHS134*sigma + crLHS60*crLHS63;
const double crLHS136 = DN(2,1)*crLHS130;
const double crLHS137 = gauss_weight*(-N[2]*crLHS23 - N[2] + crLHS131 + crLHS35*crLHS60 + crLHS73*crLHS9);
const double crLHS138 = pow(DN(2,1), 2);
const double crLHS139 = gauss_weight*(N[2] + crLHS97);
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS25 + crLHS3*crLHS9);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS26 + DN(0,1)*crLHS28 + crLHS30);
rLHS(0,2)+=DN(0,0)*crLHS38;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS39 + DN(0,1)*crLHS41 + crLHS44 + crLHS48);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS49 + DN(0,1)*crLHS51 + crLHS52);
rLHS(0,5)+=gauss_weight*(DN(1,0)*crLHS37 - crLHS23*crLHS54 + crLHS29*crLHS55 - crLHS53 + crLHS54*crLHS56);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS57 + DN(0,1)*crLHS59 + crLHS62 + crLHS66);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS67 + DN(0,1)*crLHS69 + crLHS70);
rLHS(0,8)+=gauss_weight*(DN(2,0)*crLHS37 - crLHS23*crLHS72 + crLHS29*crLHS73 + crLHS56*crLHS72 - crLHS71);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS74 + crLHS30);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS28 + DN(0,1)*crLHS75 + crLHS25 + crLHS76*crLHS9);
rLHS(1,2)+=DN(0,1)*crLHS38;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS41 + DN(0,1)*crLHS77 + crLHS79);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS51 + DN(0,1)*crLHS80 + crLHS48 + crLHS81);
rLHS(1,5)+=gauss_weight*(DN(1,1)*crLHS37 - crLHS23*crLHS83 + crLHS55*crLHS78 + crLHS56*crLHS83 - crLHS82);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS59 + DN(0,1)*crLHS84 + crLHS85);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS69 + DN(0,1)*crLHS86 + crLHS66 + crLHS87);
rLHS(1,8)+=gauss_weight*(DN(2,1)*crLHS37 - crLHS23*crLHS89 + crLHS56*crLHS89 + crLHS73*crLHS78 - crLHS88);
rLHS(2,0)+=DN(0,0)*crLHS91;
rLHS(2,1)+=DN(0,1)*crLHS91;
rLHS(2,2)+=gauss_weight*(crLHS10*crLHS92 + crLHS17*crLHS3 + crLHS17*crLHS76);
rLHS(2,3)+=gauss_weight*(DN(0,0)*crLHS93 + crLHS54);
rLHS(2,4)+=gauss_weight*(DN(0,1)*crLHS93 + crLHS83);
rLHS(2,5)+=crLHS96;
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS97 + crLHS72);
rLHS(2,7)+=gauss_weight*(DN(0,1)*crLHS97 + crLHS89);
rLHS(2,8)+=crLHS98;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS102 + crLHS44);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS26 + DN(1,1)*crLHS28 + crLHS79);
rLHS(3,2)+=gauss_weight*(DN(0,0)*crLHS104 + crLHS103*crLHS33 - crLHS23*crLHS53 + crLHS53*crLHS56 - crLHS54);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS39 + DN(1,1)*crLHS41 + crLHS105*crLHS9 + crLHS107);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS49 + DN(1,1)*crLHS51 + crLHS108);
rLHS(3,5)+=DN(1,0)*crLHS111;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS57 + DN(1,1)*crLHS59 + crLHS113 + crLHS114);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS67 + DN(1,1)*crLHS69 + crLHS115);
rLHS(3,8)+=gauss_weight*(DN(2,0)*crLHS104 + crLHS103*crLHS73 - crLHS116 - crLHS117*crLHS23 + crLHS117*crLHS56);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS74 + crLHS52);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS28 + DN(1,1)*crLHS75 + crLHS102 + crLHS81);
rLHS(4,2)+=gauss_weight*(DN(0,1)*crLHS104 + crLHS118*crLHS33 - crLHS23*crLHS82 + crLHS56*crLHS82 - crLHS83);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS41 + DN(1,1)*crLHS77 + crLHS108);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS51 + DN(1,1)*crLHS80 + crLHS107 + crLHS119*crLHS9);
rLHS(4,5)+=DN(1,1)*crLHS111;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS59 + DN(1,1)*crLHS84 + crLHS120);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS69 + DN(1,1)*crLHS86 + crLHS114 + crLHS121);
rLHS(4,8)+=gauss_weight*(DN(2,1)*crLHS104 + crLHS118*crLHS73 - crLHS122 - crLHS123*crLHS23 + crLHS123*crLHS56);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS90 + crLHS53);
rLHS(5,1)+=gauss_weight*(DN(1,1)*crLHS90 + crLHS82);
rLHS(5,2)+=crLHS96;
rLHS(5,3)+=DN(1,0)*crLHS124;
rLHS(5,4)+=DN(1,1)*crLHS124;
rLHS(5,5)+=gauss_weight*(crLHS105*crLHS17 + crLHS106*crLHS92 + crLHS119*crLHS17);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS97 + crLHS117);
rLHS(5,7)+=gauss_weight*(DN(1,1)*crLHS97 + crLHS123);
rLHS(5,8)+=crLHS125;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS129 + crLHS62);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS26 + DN(2,1)*crLHS28 + crLHS85);
rLHS(6,2)+=gauss_weight*(DN(0,0)*crLHS131 + crLHS130*crLHS33 - crLHS23*crLHS71 + crLHS56*crLHS71 - crLHS72);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS39 + DN(2,1)*crLHS41 + crLHS113 + crLHS132);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS49 + DN(2,1)*crLHS51 + crLHS120);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS131 - crLHS116*crLHS23 + crLHS116*crLHS56 - crLHS117 + crLHS130*crLHS55);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS57 + DN(2,1)*crLHS59 + crLHS133*crLHS9 + crLHS135);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS67 + DN(2,1)*crLHS69 + crLHS136);
rLHS(6,8)+=DN(2,0)*crLHS137;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS74 + crLHS70);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS28 + DN(2,1)*crLHS75 + crLHS129 + crLHS87);
rLHS(7,2)+=gauss_weight*(DN(0,1)*crLHS131 + DN(2,1)*crLHS34 - crLHS23*crLHS88 + crLHS56*crLHS88 - crLHS89);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS41 + DN(2,1)*crLHS77 + crLHS115);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS51 + DN(2,1)*crLHS80 + crLHS121 + crLHS132);
rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS131 + DN(2,1)*crLHS110 - crLHS122*crLHS23 + crLHS122*crLHS56 - crLHS123);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS59 + DN(2,1)*crLHS84 + crLHS136);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS69 + DN(2,1)*crLHS86 + crLHS135 + crLHS138*crLHS9);
rLHS(7,8)+=DN(2,1)*crLHS137;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS90 + crLHS71);
rLHS(8,1)+=gauss_weight*(DN(2,1)*crLHS90 + crLHS88);
rLHS(8,2)+=crLHS98;
rLHS(8,3)+=gauss_weight*(DN(2,0)*crLHS93 + crLHS116);
rLHS(8,4)+=gauss_weight*(DN(2,1)*crLHS93 + crLHS122);
rLHS(8,5)+=crLHS125;
rLHS(8,6)+=DN(2,0)*crLHS139;
rLHS(8,7)+=DN(2,1)*crLHS139;
rLHS(8,8)+=gauss_weight*(crLHS133*crLHS17 + crLHS134*crLHS92 + crLHS138*crLHS17);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double sigma = rData.Resistance;

    const double h = rData.ElementSize;
    const array_1d<double,4>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,3,4> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,6,6>& C = rData.C;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS1 = C(0,3)*DN(0,0);
const double crLHS2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS1;
const double crLHS3 = C(0,5)*DN(0,0);
const double crLHS4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS3;
const double crLHS5 = pow(DN(0,0), 2);
const double crLHS6 = sigma*stab_c3;
const double crLHS7 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS11 = crLHS7*stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2));
const double crLHS12 = mu + (crLHS11*h + crLHS6)/stab_c1;
const double crLHS13 = pow(N[0], 2);
const double crLHS14 = DN(0,0)*crLHS8 + DN(0,1)*crLHS9 + DN(0,2)*crLHS10;
const double crLHS15 = N[0]*crLHS7;
const double crLHS16 = bdf0*crLHS7;
const double crLHS17 = N[0]*bdf0;
const double crLHS18 = crLHS14 + crLHS17;
const double crLHS19 = pow(h, -2);
const double crLHS20 = 1.0/(crLHS11/h + crLHS19*crLHS6 + crLHS19*mu*stab_c1 + crLHS7*dyn_tau/dt);
const double crLHS21 = crLHS20*pow(crLHS7, 2);
const double crLHS22 = crLHS14*crLHS21;
const double crLHS23 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS24 = crLHS21*crLHS23;
const double crLHS25 = N[0]*crLHS24;
const double crLHS26 = crLHS20*sigma;
const double crLHS27 = crLHS15*crLHS26;
const double crLHS28 = crLHS13*crLHS16 + crLHS13*sigma + crLHS14*crLHS15 + crLHS18*crLHS22 + crLHS18*crLHS25 - crLHS18*crLHS27;
const double crLHS29 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS30 = C(1,3)*DN(0,1);
const double crLHS31 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS30;
const double crLHS32 = C(3,5)*DN(0,0);
const double crLHS33 = C(4,5)*DN(0,2);
const double crLHS34 = C(1,5)*DN(0,1) + crLHS32 + crLHS33;
const double crLHS35 = DN(0,0)*crLHS12;
const double crLHS36 = DN(0,1)*crLHS35;
const double crLHS37 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS38 = C(3,4)*DN(0,1);
const double crLHS39 = C(2,3)*DN(0,2) + crLHS32 + crLHS38;
const double crLHS40 = C(2,5)*DN(0,2);
const double crLHS41 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS40;
const double crLHS42 = DN(0,2)*crLHS35;
const double crLHS43 = N[0]*sigma;
const double crLHS44 = 1/(crLHS7*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2));
const double crLHS45 = crLHS17*crLHS44;
const double crLHS46 = crLHS12*crLHS45;
const double crLHS47 = crLHS20*crLHS23;
const double crLHS48 = crLHS20*crLHS7;
const double crLHS49 = crLHS14*crLHS48;
const double crLHS50 = gauss_weight*(-N[0] + crLHS15*crLHS47 - crLHS20*crLHS43 + crLHS46 + crLHS49);
const double crLHS51 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS52 = C(0,3)*DN(1,0);
const double crLHS53 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS52;
const double crLHS54 = C(0,5)*DN(1,0);
const double crLHS55 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS54;
const double crLHS56 = N[1]*crLHS7;
const double crLHS57 = N[1]*crLHS43 + crLHS17*crLHS56;
const double crLHS58 = DN(1,0)*crLHS35 + crLHS57;
const double crLHS59 = DN(1,0)*crLHS8 + DN(1,1)*crLHS9 + DN(1,2)*crLHS10;
const double crLHS60 = N[1]*bdf0;
const double crLHS61 = crLHS59 + crLHS60;
const double crLHS62 = crLHS15*crLHS59 + crLHS22*crLHS61 + crLHS25*crLHS61 - crLHS27*crLHS61;
const double crLHS63 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS52;
const double crLHS64 = C(1,3)*DN(1,1);
const double crLHS65 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS64;
const double crLHS66 = C(3,5)*DN(1,0);
const double crLHS67 = C(4,5)*DN(1,2);
const double crLHS68 = C(1,5)*DN(1,1) + crLHS66 + crLHS67;
const double crLHS69 = DN(1,1)*crLHS35;
const double crLHS70 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS54;
const double crLHS71 = C(3,4)*DN(1,1);
const double crLHS72 = C(2,3)*DN(1,2) + crLHS66 + crLHS71;
const double crLHS73 = C(2,5)*DN(1,2);
const double crLHS74 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS73;
const double crLHS75 = DN(1,2)*crLHS35;
const double crLHS76 = DN(0,0)*N[1];
const double crLHS77 = DN(1,0)*N[0];
const double crLHS78 = crLHS44*crLHS60;
const double crLHS79 = crLHS23*crLHS48;
const double crLHS80 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS81 = C(0,3)*DN(2,0);
const double crLHS82 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS81;
const double crLHS83 = C(0,5)*DN(2,0);
const double crLHS84 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS83;
const double crLHS85 = N[2]*crLHS7;
const double crLHS86 = N[2]*crLHS43 + crLHS17*crLHS85;
const double crLHS87 = DN(2,0)*crLHS35 + crLHS86;
const double crLHS88 = DN(2,0)*crLHS8 + DN(2,1)*crLHS9 + DN(2,2)*crLHS10;
const double crLHS89 = N[2]*bdf0;
const double crLHS90 = crLHS88 + crLHS89;
const double crLHS91 = crLHS15*crLHS88 + crLHS22*crLHS90 + crLHS25*crLHS90 - crLHS27*crLHS90;
const double crLHS92 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS81;
const double crLHS93 = C(1,3)*DN(2,1);
const double crLHS94 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS93;
const double crLHS95 = C(3,5)*DN(2,0);
const double crLHS96 = C(4,5)*DN(2,2);
const double crLHS97 = C(1,5)*DN(2,1) + crLHS95 + crLHS96;
const double crLHS98 = DN(2,1)*crLHS35;
const double crLHS99 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS83;
const double crLHS100 = C(3,4)*DN(2,1);
const double crLHS101 = C(2,3)*DN(2,2) + crLHS100 + crLHS95;
const double crLHS102 = C(2,5)*DN(2,2);
const double crLHS103 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS102;
const double crLHS104 = DN(2,2)*crLHS35;
const double crLHS105 = DN(0,0)*N[2];
const double crLHS106 = DN(2,0)*N[0];
const double crLHS107 = crLHS44*crLHS89;
const double crLHS108 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS109 = C(0,3)*DN(3,0);
const double crLHS110 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS109;
const double crLHS111 = C(0,5)*DN(3,0);
const double crLHS112 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS111;
const double crLHS113 = N[3]*crLHS7;
const double crLHS114 = N[3]*crLHS43 + crLHS113*crLHS17;
const double crLHS115 = DN(3,0)*crLHS35 + crLHS114;
const double crLHS116 = DN(3,0)*crLHS8 + DN(3,1)*crLHS9 + DN(3,2)*crLHS10;
const double crLHS117 = N[3]*bdf0;
const double crLHS118 = crLHS116 + crLHS117;
const double crLHS119 = crLHS116*crLHS15 + crLHS118*crLHS22 + crLHS118*crLHS25 - crLHS118*crLHS27;
const double crLHS120 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS109;
const double crLHS121 = C(1,3)*DN(3,1);
const double crLHS122 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS121;
const double crLHS123 = C(3,5)*DN(3,0);
const double crLHS124 = C(4,5)*DN(3,2);
const double crLHS125 = C(1,5)*DN(3,1) + crLHS123 + crLHS124;
const double crLHS126 = DN(3,1)*crLHS35;
const double crLHS127 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS111;
const double crLHS128 = C(3,4)*DN(3,1);
const double crLHS129 = C(2,3)*DN(3,2) + crLHS123 + crLHS128;
const double crLHS130 = C(2,5)*DN(3,2);
const double crLHS131 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS130;
const double crLHS132 = DN(3,2)*crLHS35;
const double crLHS133 = DN(0,0)*N[3];
const double crLHS134 = DN(3,0)*N[0];
const double crLHS135 = crLHS117*crLHS44;
const double crLHS136 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS30;
const double crLHS137 = C(0,4)*DN(0,0) + crLHS33 + crLHS38;
const double crLHS138 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS139 = C(1,4)*DN(0,1);
const double crLHS140 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS139;
const double crLHS141 = pow(DN(0,1), 2);
const double crLHS142 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS139;
const double crLHS143 = C(2,4)*DN(0,2);
const double crLHS144 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS143;
const double crLHS145 = DN(0,1)*crLHS12;
const double crLHS146 = DN(0,2)*crLHS145;
const double crLHS147 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS64;
const double crLHS148 = C(0,4)*DN(1,0) + crLHS67 + crLHS71;
const double crLHS149 = DN(1,0)*crLHS145;
const double crLHS150 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS151 = C(1,4)*DN(1,1);
const double crLHS152 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS151;
const double crLHS153 = DN(1,1)*crLHS145;
const double crLHS154 = crLHS57 + crLHS62;
const double crLHS155 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS151;
const double crLHS156 = C(2,4)*DN(1,2);
const double crLHS157 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS156;
const double crLHS158 = DN(1,2)*crLHS145;
const double crLHS159 = DN(0,1)*N[1];
const double crLHS160 = DN(1,1)*N[0];
const double crLHS161 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS93;
const double crLHS162 = C(0,4)*DN(2,0) + crLHS100 + crLHS96;
const double crLHS163 = DN(2,0)*crLHS145;
const double crLHS164 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS165 = C(1,4)*DN(2,1);
const double crLHS166 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS165;
const double crLHS167 = DN(2,1)*crLHS145;
const double crLHS168 = crLHS86 + crLHS91;
const double crLHS169 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS165;
const double crLHS170 = C(2,4)*DN(2,2);
const double crLHS171 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS170;
const double crLHS172 = DN(2,2)*crLHS145;
const double crLHS173 = DN(0,1)*N[2];
const double crLHS174 = DN(2,1)*N[0];
const double crLHS175 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS121;
const double crLHS176 = C(0,4)*DN(3,0) + crLHS124 + crLHS128;
const double crLHS177 = DN(3,0)*crLHS145;
const double crLHS178 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS179 = C(1,4)*DN(3,1);
const double crLHS180 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS179;
const double crLHS181 = DN(3,1)*crLHS145;
const double crLHS182 = crLHS114 + crLHS119;
const double crLHS183 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS179;
const double crLHS184 = C(2,4)*DN(3,2);
const double crLHS185 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS184;
const double crLHS186 = DN(3,2)*crLHS145;
const double crLHS187 = DN(0,1)*N[3];
const double crLHS188 = DN(3,1)*N[0];
const double crLHS189 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS40;
const double crLHS190 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS143;
const double crLHS191 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS192 = pow(DN(0,2), 2);
const double crLHS193 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS73;
const double crLHS194 = DN(0,2)*crLHS12;
const double crLHS195 = DN(1,0)*crLHS194;
const double crLHS196 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS156;
const double crLHS197 = DN(1,1)*crLHS194;
const double crLHS198 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS199 = DN(1,2)*crLHS194;
const double crLHS200 = DN(0,2)*N[1];
const double crLHS201 = DN(1,2)*N[0];
const double crLHS202 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS102;
const double crLHS203 = DN(2,0)*crLHS194;
const double crLHS204 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS170;
const double crLHS205 = DN(2,1)*crLHS194;
const double crLHS206 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS207 = DN(2,2)*crLHS194;
const double crLHS208 = DN(0,2)*N[2];
const double crLHS209 = DN(2,2)*N[0];
const double crLHS210 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS130;
const double crLHS211 = DN(3,0)*crLHS194;
const double crLHS212 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS184;
const double crLHS213 = DN(3,1)*crLHS194;
const double crLHS214 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS215 = DN(3,2)*crLHS194;
const double crLHS216 = DN(0,2)*N[3];
const double crLHS217 = DN(3,2)*N[0];
const double crLHS218 = crLHS18*crLHS48;
const double crLHS219 = gauss_weight*(N[0] + crLHS218);
const double crLHS220 = bdf0*crLHS44;
const double crLHS221 = crLHS48*crLHS61;
const double crLHS222 = DN(0,0)*crLHS20;
const double crLHS223 = DN(0,1)*crLHS20;
const double crLHS224 = DN(0,2)*crLHS20;
const double crLHS225 = gauss_weight*(DN(1,0)*crLHS222 + DN(1,1)*crLHS223 + DN(1,2)*crLHS224 + N[1]*crLHS45);
const double crLHS226 = crLHS48*crLHS90;
const double crLHS227 = gauss_weight*(DN(2,0)*crLHS222 + DN(2,1)*crLHS223 + DN(2,2)*crLHS224 + N[2]*crLHS45);
const double crLHS228 = crLHS118*crLHS48;
const double crLHS229 = gauss_weight*(DN(3,0)*crLHS222 + DN(3,1)*crLHS223 + DN(3,2)*crLHS224 + N[3]*crLHS45);
const double crLHS230 = crLHS21*crLHS59;
const double crLHS231 = N[1]*crLHS24;
const double crLHS232 = crLHS26*crLHS56;
const double crLHS233 = crLHS14*crLHS56 + crLHS18*crLHS230 + crLHS18*crLHS231 - crLHS18*crLHS232;
const double crLHS234 = DN(1,0)*crLHS12;
const double crLHS235 = crLHS48*crLHS59;
const double crLHS236 = pow(DN(1,0), 2);
const double crLHS237 = pow(N[1], 2);
const double crLHS238 = crLHS16*crLHS237 + crLHS230*crLHS61 + crLHS231*crLHS61 - crLHS232*crLHS61 + crLHS237*sigma + crLHS56*crLHS59;
const double crLHS239 = DN(1,1)*crLHS234;
const double crLHS240 = DN(1,2)*crLHS234;
const double crLHS241 = N[1]*sigma;
const double crLHS242 = crLHS12*crLHS78;
const double crLHS243 = gauss_weight*(-N[1] - crLHS20*crLHS241 + crLHS235 + crLHS242 + crLHS47*crLHS56);
const double crLHS244 = N[2]*crLHS241 + crLHS60*crLHS85;
const double crLHS245 = DN(2,0)*crLHS234 + crLHS244;
const double crLHS246 = crLHS230*crLHS90 + crLHS231*crLHS90 - crLHS232*crLHS90 + crLHS56*crLHS88;
const double crLHS247 = DN(2,1)*crLHS234;
const double crLHS248 = DN(2,2)*crLHS234;
const double crLHS249 = DN(1,0)*N[2];
const double crLHS250 = DN(2,0)*N[1];
const double crLHS251 = N[3]*crLHS241 + crLHS113*crLHS60;
const double crLHS252 = DN(3,0)*crLHS234 + crLHS251;
const double crLHS253 = crLHS116*crLHS56 + crLHS118*crLHS230 + crLHS118*crLHS231 - crLHS118*crLHS232;
const double crLHS254 = DN(3,1)*crLHS234;
const double crLHS255 = DN(3,2)*crLHS234;
const double crLHS256 = DN(1,0)*N[3];
const double crLHS257 = DN(3,0)*N[1];
const double crLHS258 = crLHS233 + crLHS57;
const double crLHS259 = DN(1,1)*crLHS12;
const double crLHS260 = pow(DN(1,1), 2);
const double crLHS261 = DN(1,2)*crLHS259;
const double crLHS262 = DN(2,0)*crLHS259;
const double crLHS263 = DN(2,1)*crLHS259;
const double crLHS264 = crLHS244 + crLHS246;
const double crLHS265 = DN(2,2)*crLHS259;
const double crLHS266 = DN(1,1)*N[2];
const double crLHS267 = DN(2,1)*N[1];
const double crLHS268 = DN(3,0)*crLHS259;
const double crLHS269 = DN(3,1)*crLHS259;
const double crLHS270 = crLHS251 + crLHS253;
const double crLHS271 = DN(3,2)*crLHS259;
const double crLHS272 = DN(1,1)*N[3];
const double crLHS273 = DN(3,1)*N[1];
const double crLHS274 = DN(1,2)*crLHS12;
const double crLHS275 = pow(DN(1,2), 2);
const double crLHS276 = DN(2,0)*crLHS274;
const double crLHS277 = DN(2,1)*crLHS274;
const double crLHS278 = DN(2,2)*crLHS274;
const double crLHS279 = DN(1,2)*N[2];
const double crLHS280 = DN(2,2)*N[1];
const double crLHS281 = DN(3,0)*crLHS274;
const double crLHS282 = DN(3,1)*crLHS274;
const double crLHS283 = DN(3,2)*crLHS274;
const double crLHS284 = DN(1,2)*N[3];
const double crLHS285 = DN(3,2)*N[1];
const double crLHS286 = gauss_weight*(N[1] + crLHS221);
const double crLHS287 = DN(1,0)*crLHS20;
const double crLHS288 = DN(1,1)*crLHS20;
const double crLHS289 = DN(1,2)*crLHS20;
const double crLHS290 = gauss_weight*(DN(2,0)*crLHS287 + DN(2,1)*crLHS288 + DN(2,2)*crLHS289 + N[2]*crLHS78);
const double crLHS291 = gauss_weight*(DN(3,0)*crLHS287 + DN(3,1)*crLHS288 + DN(3,2)*crLHS289 + N[3]*crLHS78);
const double crLHS292 = crLHS21*crLHS88;
const double crLHS293 = N[2]*crLHS24;
const double crLHS294 = crLHS26*crLHS85;
const double crLHS295 = crLHS14*crLHS85 + crLHS18*crLHS292 + crLHS18*crLHS293 - crLHS18*crLHS294;
const double crLHS296 = DN(2,0)*crLHS12;
const double crLHS297 = crLHS48*crLHS88;
const double crLHS298 = crLHS292*crLHS61 + crLHS293*crLHS61 - crLHS294*crLHS61 + crLHS59*crLHS85;
const double crLHS299 = pow(DN(2,0), 2);
const double crLHS300 = pow(N[2], 2);
const double crLHS301 = crLHS16*crLHS300 + crLHS292*crLHS90 + crLHS293*crLHS90 - crLHS294*crLHS90 + crLHS300*sigma + crLHS85*crLHS88;
const double crLHS302 = DN(2,1)*crLHS296;
const double crLHS303 = DN(2,2)*crLHS296;
const double crLHS304 = N[2]*sigma;
const double crLHS305 = crLHS107*crLHS12;
const double crLHS306 = gauss_weight*(-N[2] - crLHS20*crLHS304 + crLHS297 + crLHS305 + crLHS47*crLHS85);
const double crLHS307 = N[3]*crLHS304 + crLHS113*crLHS89;
const double crLHS308 = DN(3,0)*crLHS296 + crLHS307;
const double crLHS309 = crLHS116*crLHS85 + crLHS118*crLHS292 + crLHS118*crLHS293 - crLHS118*crLHS294;
const double crLHS310 = DN(3,1)*crLHS296;
const double crLHS311 = DN(3,2)*crLHS296;
const double crLHS312 = DN(2,0)*N[3];
const double crLHS313 = DN(3,0)*N[2];
const double crLHS314 = crLHS295 + crLHS86;
const double crLHS315 = DN(2,1)*crLHS12;
const double crLHS316 = crLHS244 + crLHS298;
const double crLHS317 = pow(DN(2,1), 2);
const double crLHS318 = DN(2,2)*crLHS315;
const double crLHS319 = DN(3,0)*crLHS315;
const double crLHS320 = DN(3,1)*crLHS315;
const double crLHS321 = crLHS307 + crLHS309;
const double crLHS322 = DN(3,2)*crLHS315;
const double crLHS323 = DN(2,1)*N[3];
const double crLHS324 = DN(3,1)*N[2];
const double crLHS325 = DN(2,2)*crLHS12;
const double crLHS326 = pow(DN(2,2), 2);
const double crLHS327 = DN(3,0)*crLHS325;
const double crLHS328 = DN(3,1)*crLHS325;
const double crLHS329 = DN(3,2)*crLHS325;
const double crLHS330 = DN(2,2)*N[3];
const double crLHS331 = DN(3,2)*N[2];
const double crLHS332 = gauss_weight*(N[2] + crLHS226);
const double crLHS333 = gauss_weight*(DN(2,0)*DN(3,0)*crLHS20 + DN(2,1)*DN(3,1)*crLHS20 + DN(2,2)*DN(3,2)*crLHS20 + N[3]*crLHS107);
const double crLHS334 = crLHS116*crLHS21;
const double crLHS335 = N[3]*crLHS24;
const double crLHS336 = crLHS113*crLHS26;
const double crLHS337 = crLHS113*crLHS14 + crLHS18*crLHS334 + crLHS18*crLHS335 - crLHS18*crLHS336;
const double crLHS338 = DN(3,0)*crLHS12;
const double crLHS339 = crLHS116*crLHS48;
const double crLHS340 = crLHS113*crLHS59 + crLHS334*crLHS61 + crLHS335*crLHS61 - crLHS336*crLHS61;
const double crLHS341 = crLHS113*crLHS88 + crLHS334*crLHS90 + crLHS335*crLHS90 - crLHS336*crLHS90;
const double crLHS342 = pow(DN(3,0), 2);
const double crLHS343 = pow(N[3], 2);
const double crLHS344 = crLHS113*crLHS116 + crLHS118*crLHS334 + crLHS118*crLHS335 - crLHS118*crLHS336 + crLHS16*crLHS343 + crLHS343*sigma;
const double crLHS345 = DN(3,1)*crLHS338;
const double crLHS346 = DN(3,2)*crLHS338;
const double crLHS347 = gauss_weight*(-N[3]*crLHS26 - N[3] + crLHS113*crLHS47 + crLHS12*crLHS135 + crLHS339);
const double crLHS348 = crLHS114 + crLHS337;
const double crLHS349 = DN(3,1)*crLHS12;
const double crLHS350 = crLHS251 + crLHS340;
const double crLHS351 = crLHS307 + crLHS341;
const double crLHS352 = pow(DN(3,1), 2);
const double crLHS353 = DN(3,2)*crLHS349;
const double crLHS354 = pow(DN(3,2), 2);
const double crLHS355 = gauss_weight*(N[3] + crLHS228);
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS28);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS29 + DN(0,1)*crLHS31 + DN(0,2)*crLHS34 + crLHS36);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS37 + DN(0,1)*crLHS39 + DN(0,2)*crLHS41 + crLHS42);
rLHS(0,3)+=DN(0,0)*crLHS50;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS51 + DN(0,1)*crLHS53 + DN(0,2)*crLHS55 + crLHS58 + crLHS62);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS65 + DN(0,2)*crLHS68 + crLHS69);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS70 + DN(0,1)*crLHS72 + DN(0,2)*crLHS74 + crLHS75);
rLHS(0,7)+=gauss_weight*(DN(1,0)*crLHS49 - crLHS26*crLHS77 + crLHS35*crLHS78 - crLHS76 + crLHS77*crLHS79);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS80 + DN(0,1)*crLHS82 + DN(0,2)*crLHS84 + crLHS87 + crLHS91);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS92 + DN(0,1)*crLHS94 + DN(0,2)*crLHS97 + crLHS98);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS99 + DN(0,1)*crLHS101 + DN(0,2)*crLHS103 + crLHS104);
rLHS(0,11)+=gauss_weight*(DN(2,0)*crLHS49 - crLHS105 - crLHS106*crLHS26 + crLHS106*crLHS79 + crLHS107*crLHS35);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS108 + DN(0,1)*crLHS110 + DN(0,2)*crLHS112 + crLHS115 + crLHS119);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS120 + DN(0,1)*crLHS122 + DN(0,2)*crLHS125 + crLHS126);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS129 + DN(0,2)*crLHS131 + crLHS132);
rLHS(0,15)+=gauss_weight*(DN(3,0)*crLHS49 - crLHS133 - crLHS134*crLHS26 + crLHS134*crLHS79 + crLHS135*crLHS35);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS136 + DN(0,2)*crLHS137 + crLHS36);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS31 + DN(0,1)*crLHS138 + DN(0,2)*crLHS140 + crLHS12*crLHS141 + crLHS28);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS39 + DN(0,1)*crLHS142 + DN(0,2)*crLHS144 + crLHS146);
rLHS(1,3)+=DN(0,1)*crLHS50;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS53 + DN(0,1)*crLHS147 + DN(0,2)*crLHS148 + crLHS149);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS65 + DN(0,1)*crLHS150 + DN(0,2)*crLHS152 + crLHS153 + crLHS154);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS72 + DN(0,1)*crLHS155 + DN(0,2)*crLHS157 + crLHS158);
rLHS(1,7)+=gauss_weight*(DN(1,1)*crLHS49 + crLHS145*crLHS78 - crLHS159 - crLHS160*crLHS26 + crLHS160*crLHS79);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS82 + DN(0,1)*crLHS161 + DN(0,2)*crLHS162 + crLHS163);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS94 + DN(0,1)*crLHS164 + DN(0,2)*crLHS166 + crLHS167 + crLHS168);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS101 + DN(0,1)*crLHS169 + DN(0,2)*crLHS171 + crLHS172);
rLHS(1,11)+=gauss_weight*(DN(2,1)*crLHS49 + crLHS107*crLHS145 - crLHS173 - crLHS174*crLHS26 + crLHS174*crLHS79);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS175 + DN(0,2)*crLHS176 + crLHS177);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS122 + DN(0,1)*crLHS178 + DN(0,2)*crLHS180 + crLHS181 + crLHS182);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS129 + DN(0,1)*crLHS183 + DN(0,2)*crLHS185 + crLHS186);
rLHS(1,15)+=gauss_weight*(DN(3,1)*crLHS49 + crLHS135*crLHS145 - crLHS187 - crLHS188*crLHS26 + crLHS188*crLHS79);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS137 + DN(0,2)*crLHS189 + crLHS42);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS34 + DN(0,1)*crLHS140 + DN(0,2)*crLHS190 + crLHS146);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS41 + DN(0,1)*crLHS144 + DN(0,2)*crLHS191 + crLHS12*crLHS192 + crLHS28);
rLHS(2,3)+=DN(0,2)*crLHS50;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS55 + DN(0,1)*crLHS148 + DN(0,2)*crLHS193 + crLHS195);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS68 + DN(0,1)*crLHS152 + DN(0,2)*crLHS196 + crLHS197);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS157 + DN(0,2)*crLHS198 + crLHS154 + crLHS199);
rLHS(2,7)+=gauss_weight*(DN(1,2)*crLHS49 + crLHS194*crLHS78 - crLHS200 - crLHS201*crLHS26 + crLHS201*crLHS79);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS84 + DN(0,1)*crLHS162 + DN(0,2)*crLHS202 + crLHS203);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS97 + DN(0,1)*crLHS166 + DN(0,2)*crLHS204 + crLHS205);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS103 + DN(0,1)*crLHS171 + DN(0,2)*crLHS206 + crLHS168 + crLHS207);
rLHS(2,11)+=gauss_weight*(DN(2,2)*crLHS49 + crLHS107*crLHS194 - crLHS208 - crLHS209*crLHS26 + crLHS209*crLHS79);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS176 + DN(0,2)*crLHS210 + crLHS211);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS125 + DN(0,1)*crLHS180 + DN(0,2)*crLHS212 + crLHS213);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS131 + DN(0,1)*crLHS185 + DN(0,2)*crLHS214 + crLHS182 + crLHS215);
rLHS(2,15)+=gauss_weight*(DN(3,2)*crLHS49 + crLHS135*crLHS194 - crLHS216 - crLHS217*crLHS26 + crLHS217*crLHS79);
rLHS(3,0)+=DN(0,0)*crLHS219;
rLHS(3,1)+=DN(0,1)*crLHS219;
rLHS(3,2)+=DN(0,2)*crLHS219;
rLHS(3,3)+=gauss_weight*(crLHS13*crLHS220 + crLHS141*crLHS20 + crLHS192*crLHS20 + crLHS20*crLHS5);
rLHS(3,4)+=gauss_weight*(DN(0,0)*crLHS221 + crLHS77);
rLHS(3,5)+=gauss_weight*(DN(0,1)*crLHS221 + crLHS160);
rLHS(3,6)+=gauss_weight*(DN(0,2)*crLHS221 + crLHS201);
rLHS(3,7)+=crLHS225;
rLHS(3,8)+=gauss_weight*(DN(0,0)*crLHS226 + crLHS106);
rLHS(3,9)+=gauss_weight*(DN(0,1)*crLHS226 + crLHS174);
rLHS(3,10)+=gauss_weight*(DN(0,2)*crLHS226 + crLHS209);
rLHS(3,11)+=crLHS227;
rLHS(3,12)+=gauss_weight*(DN(0,0)*crLHS228 + crLHS134);
rLHS(3,13)+=gauss_weight*(DN(0,1)*crLHS228 + crLHS188);
rLHS(3,14)+=gauss_weight*(DN(0,2)*crLHS228 + crLHS217);
rLHS(3,15)+=crLHS229;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS233 + crLHS58);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS29 + DN(1,1)*crLHS31 + DN(1,2)*crLHS34 + crLHS149);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS37 + DN(1,1)*crLHS39 + DN(1,2)*crLHS41 + crLHS195);
rLHS(4,3)+=gauss_weight*(DN(0,0)*crLHS235 + crLHS234*crLHS45 - crLHS26*crLHS76 + crLHS76*crLHS79 - crLHS77);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS51 + DN(1,1)*crLHS53 + DN(1,2)*crLHS55 + crLHS12*crLHS236 + crLHS238);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS65 + DN(1,2)*crLHS68 + crLHS239);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS70 + DN(1,1)*crLHS72 + DN(1,2)*crLHS74 + crLHS240);
rLHS(4,7)+=DN(1,0)*crLHS243;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS80 + DN(1,1)*crLHS82 + DN(1,2)*crLHS84 + crLHS245 + crLHS246);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS92 + DN(1,1)*crLHS94 + DN(1,2)*crLHS97 + crLHS247);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS99 + DN(1,1)*crLHS101 + DN(1,2)*crLHS103 + crLHS248);
rLHS(4,11)+=gauss_weight*(DN(2,0)*crLHS235 + crLHS107*crLHS234 - crLHS249 - crLHS250*crLHS26 + crLHS250*crLHS79);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS108 + DN(1,1)*crLHS110 + DN(1,2)*crLHS112 + crLHS252 + crLHS253);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS120 + DN(1,1)*crLHS122 + DN(1,2)*crLHS125 + crLHS254);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS129 + DN(1,2)*crLHS131 + crLHS255);
rLHS(4,15)+=gauss_weight*(DN(3,0)*crLHS235 + crLHS135*crLHS234 - crLHS256 - crLHS257*crLHS26 + crLHS257*crLHS79);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS136 + DN(1,2)*crLHS137 + crLHS69);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS31 + DN(1,1)*crLHS138 + DN(1,2)*crLHS140 + crLHS153 + crLHS258);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS39 + DN(1,1)*crLHS142 + DN(1,2)*crLHS144 + crLHS197);
rLHS(5,3)+=gauss_weight*(DN(0,1)*crLHS235 - crLHS159*crLHS26 + crLHS159*crLHS79 - crLHS160 + crLHS259*crLHS45);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS53 + DN(1,1)*crLHS147 + DN(1,2)*crLHS148 + crLHS239);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS65 + DN(1,1)*crLHS150 + DN(1,2)*crLHS152 + crLHS12*crLHS260 + crLHS238);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS72 + DN(1,1)*crLHS155 + DN(1,2)*crLHS157 + crLHS261);
rLHS(5,7)+=DN(1,1)*crLHS243;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS82 + DN(1,1)*crLHS161 + DN(1,2)*crLHS162 + crLHS262);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS94 + DN(1,1)*crLHS164 + DN(1,2)*crLHS166 + crLHS263 + crLHS264);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS101 + DN(1,1)*crLHS169 + DN(1,2)*crLHS171 + crLHS265);
rLHS(5,11)+=gauss_weight*(DN(2,1)*crLHS235 + crLHS107*crLHS259 - crLHS26*crLHS267 - crLHS266 + crLHS267*crLHS79);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS175 + DN(1,2)*crLHS176 + crLHS268);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS122 + DN(1,1)*crLHS178 + DN(1,2)*crLHS180 + crLHS269 + crLHS270);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS129 + DN(1,1)*crLHS183 + DN(1,2)*crLHS185 + crLHS271);
rLHS(5,15)+=gauss_weight*(DN(3,1)*crLHS235 + crLHS135*crLHS259 - crLHS26*crLHS273 - crLHS272 + crLHS273*crLHS79);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS137 + DN(1,2)*crLHS189 + crLHS75);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS34 + DN(1,1)*crLHS140 + DN(1,2)*crLHS190 + crLHS158);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS41 + DN(1,1)*crLHS144 + DN(1,2)*crLHS191 + crLHS199 + crLHS258);
rLHS(6,3)+=gauss_weight*(DN(0,2)*crLHS235 - crLHS200*crLHS26 + crLHS200*crLHS79 - crLHS201 + crLHS274*crLHS45);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS55 + DN(1,1)*crLHS148 + DN(1,2)*crLHS193 + crLHS240);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS68 + DN(1,1)*crLHS152 + DN(1,2)*crLHS196 + crLHS261);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS157 + DN(1,2)*crLHS198 + crLHS12*crLHS275 + crLHS238);
rLHS(6,7)+=DN(1,2)*crLHS243;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS84 + DN(1,1)*crLHS162 + DN(1,2)*crLHS202 + crLHS276);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS97 + DN(1,1)*crLHS166 + DN(1,2)*crLHS204 + crLHS277);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS103 + DN(1,1)*crLHS171 + DN(1,2)*crLHS206 + crLHS264 + crLHS278);
rLHS(6,11)+=gauss_weight*(DN(2,2)*crLHS235 + crLHS107*crLHS274 - crLHS26*crLHS280 - crLHS279 + crLHS280*crLHS79);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS176 + DN(1,2)*crLHS210 + crLHS281);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS125 + DN(1,1)*crLHS180 + DN(1,2)*crLHS212 + crLHS282);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS131 + DN(1,1)*crLHS185 + DN(1,2)*crLHS214 + crLHS270 + crLHS283);
rLHS(6,15)+=gauss_weight*(DN(3,2)*crLHS235 + crLHS135*crLHS274 - crLHS26*crLHS285 - crLHS284 + crLHS285*crLHS79);
rLHS(7,0)+=gauss_weight*(DN(1,0)*crLHS218 + crLHS76);
rLHS(7,1)+=gauss_weight*(DN(1,1)*crLHS218 + crLHS159);
rLHS(7,2)+=gauss_weight*(DN(1,2)*crLHS218 + crLHS200);
rLHS(7,3)+=crLHS225;
rLHS(7,4)+=DN(1,0)*crLHS286;
rLHS(7,5)+=DN(1,1)*crLHS286;
rLHS(7,6)+=DN(1,2)*crLHS286;
rLHS(7,7)+=gauss_weight*(crLHS20*crLHS236 + crLHS20*crLHS260 + crLHS20*crLHS275 + crLHS220*crLHS237);
rLHS(7,8)+=gauss_weight*(DN(1,0)*crLHS226 + crLHS250);
rLHS(7,9)+=gauss_weight*(DN(1,1)*crLHS226 + crLHS267);
rLHS(7,10)+=gauss_weight*(DN(1,2)*crLHS226 + crLHS280);
rLHS(7,11)+=crLHS290;
rLHS(7,12)+=gauss_weight*(DN(1,0)*crLHS228 + crLHS257);
rLHS(7,13)+=gauss_weight*(DN(1,1)*crLHS228 + crLHS273);
rLHS(7,14)+=gauss_weight*(DN(1,2)*crLHS228 + crLHS285);
rLHS(7,15)+=crLHS291;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS295 + crLHS87);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS29 + DN(2,1)*crLHS31 + DN(2,2)*crLHS34 + crLHS163);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS37 + DN(2,1)*crLHS39 + DN(2,2)*crLHS41 + crLHS203);
rLHS(8,3)+=gauss_weight*(DN(0,0)*crLHS297 - crLHS105*crLHS26 + crLHS105*crLHS79 - crLHS106 + crLHS296*crLHS45);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS51 + DN(2,1)*crLHS53 + DN(2,2)*crLHS55 + crLHS245 + crLHS298);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS65 + DN(2,2)*crLHS68 + crLHS262);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS70 + DN(2,1)*crLHS72 + DN(2,2)*crLHS74 + crLHS276);
rLHS(8,7)+=gauss_weight*(DN(1,0)*crLHS297 - crLHS249*crLHS26 + crLHS249*crLHS79 - crLHS250 + crLHS296*crLHS78);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS80 + DN(2,1)*crLHS82 + DN(2,2)*crLHS84 + crLHS12*crLHS299 + crLHS301);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS92 + DN(2,1)*crLHS94 + DN(2,2)*crLHS97 + crLHS302);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS99 + DN(2,1)*crLHS101 + DN(2,2)*crLHS103 + crLHS303);
rLHS(8,11)+=DN(2,0)*crLHS306;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS108 + DN(2,1)*crLHS110 + DN(2,2)*crLHS112 + crLHS308 + crLHS309);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS120 + DN(2,1)*crLHS122 + DN(2,2)*crLHS125 + crLHS310);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS129 + DN(2,2)*crLHS131 + crLHS311);
rLHS(8,15)+=gauss_weight*(DN(3,0)*crLHS297 + crLHS135*crLHS296 - crLHS26*crLHS313 - crLHS312 + crLHS313*crLHS79);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS136 + DN(2,2)*crLHS137 + crLHS98);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS31 + DN(2,1)*crLHS138 + DN(2,2)*crLHS140 + crLHS167 + crLHS314);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS39 + DN(2,1)*crLHS142 + DN(2,2)*crLHS144 + crLHS205);
rLHS(9,3)+=gauss_weight*(DN(0,1)*crLHS297 - crLHS173*crLHS26 + crLHS173*crLHS79 - crLHS174 + crLHS315*crLHS45);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS53 + DN(2,1)*crLHS147 + DN(2,2)*crLHS148 + crLHS247);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS65 + DN(2,1)*crLHS150 + DN(2,2)*crLHS152 + crLHS263 + crLHS316);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS72 + DN(2,1)*crLHS155 + DN(2,2)*crLHS157 + crLHS277);
rLHS(9,7)+=gauss_weight*(DN(1,1)*crLHS297 - crLHS26*crLHS266 + crLHS266*crLHS79 - crLHS267 + crLHS315*crLHS78);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS82 + DN(2,1)*crLHS161 + DN(2,2)*crLHS162 + crLHS302);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS94 + DN(2,1)*crLHS164 + DN(2,2)*crLHS166 + crLHS12*crLHS317 + crLHS301);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS101 + DN(2,1)*crLHS169 + DN(2,2)*crLHS171 + crLHS318);
rLHS(9,11)+=DN(2,1)*crLHS306;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS175 + DN(2,2)*crLHS176 + crLHS319);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS122 + DN(2,1)*crLHS178 + DN(2,2)*crLHS180 + crLHS320 + crLHS321);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS129 + DN(2,1)*crLHS183 + DN(2,2)*crLHS185 + crLHS322);
rLHS(9,15)+=gauss_weight*(DN(3,1)*crLHS297 + crLHS135*crLHS315 - crLHS26*crLHS324 - crLHS323 + crLHS324*crLHS79);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS137 + DN(2,2)*crLHS189 + crLHS104);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS34 + DN(2,1)*crLHS140 + DN(2,2)*crLHS190 + crLHS172);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS41 + DN(2,1)*crLHS144 + DN(2,2)*crLHS191 + crLHS207 + crLHS314);
rLHS(10,3)+=gauss_weight*(DN(0,2)*crLHS297 - crLHS208*crLHS26 + crLHS208*crLHS79 - crLHS209 + crLHS325*crLHS45);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS55 + DN(2,1)*crLHS148 + DN(2,2)*crLHS193 + crLHS248);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS68 + DN(2,1)*crLHS152 + DN(2,2)*crLHS196 + crLHS265);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS157 + DN(2,2)*crLHS198 + crLHS278 + crLHS316);
rLHS(10,7)+=gauss_weight*(DN(1,2)*crLHS297 - crLHS26*crLHS279 + crLHS279*crLHS79 - crLHS280 + crLHS325*crLHS78);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS84 + DN(2,1)*crLHS162 + DN(2,2)*crLHS202 + crLHS303);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS97 + DN(2,1)*crLHS166 + DN(2,2)*crLHS204 + crLHS318);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS103 + DN(2,1)*crLHS171 + DN(2,2)*crLHS206 + crLHS12*crLHS326 + crLHS301);
rLHS(10,11)+=DN(2,2)*crLHS306;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS176 + DN(2,2)*crLHS210 + crLHS327);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS125 + DN(2,1)*crLHS180 + DN(2,2)*crLHS212 + crLHS328);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS131 + DN(2,1)*crLHS185 + DN(2,2)*crLHS214 + crLHS321 + crLHS329);
rLHS(10,15)+=gauss_weight*(DN(3,2)*crLHS297 + crLHS135*crLHS325 - crLHS26*crLHS331 - crLHS330 + crLHS331*crLHS79);
rLHS(11,0)+=gauss_weight*(DN(2,0)*crLHS218 + crLHS105);
rLHS(11,1)+=gauss_weight*(DN(2,1)*crLHS218 + crLHS173);
rLHS(11,2)+=gauss_weight*(DN(2,2)*crLHS218 + crLHS208);
rLHS(11,3)+=crLHS227;
rLHS(11,4)+=gauss_weight*(DN(2,0)*crLHS221 + crLHS249);
rLHS(11,5)+=gauss_weight*(DN(2,1)*crLHS221 + crLHS266);
rLHS(11,6)+=gauss_weight*(DN(2,2)*crLHS221 + crLHS279);
rLHS(11,7)+=crLHS290;
rLHS(11,8)+=DN(2,0)*crLHS332;
rLHS(11,9)+=DN(2,1)*crLHS332;
rLHS(11,10)+=DN(2,2)*crLHS332;
rLHS(11,11)+=gauss_weight*(crLHS20*crLHS299 + crLHS20*crLHS317 + crLHS20*crLHS326 + crLHS220*crLHS300);
rLHS(11,12)+=gauss_weight*(DN(2,0)*crLHS228 + crLHS313);
rLHS(11,13)+=gauss_weight*(DN(2,1)*crLHS228 + crLHS324);
rLHS(11,14)+=gauss_weight*(DN(2,2)*crLHS228 + crLHS331);
rLHS(11,15)+=crLHS333;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS115 + crLHS337);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS29 + DN(3,1)*crLHS31 + DN(3,2)*crLHS34 + crLHS177);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS37 + DN(3,1)*crLHS39 + DN(3,2)*crLHS41 + crLHS211);
rLHS(12,3)+=gauss_weight*(DN(0,0)*crLHS339 - crLHS133*crLHS26 + crLHS133*crLHS79 - crLHS134 + crLHS338*crLHS45);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS51 + DN(3,1)*crLHS53 + DN(3,2)*crLHS55 + crLHS252 + crLHS340);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS63 + DN(3,1)*crLHS65 + DN(3,2)*crLHS68 + crLHS268);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS70 + DN(3,1)*crLHS72 + DN(3,2)*crLHS74 + crLHS281);
rLHS(12,7)+=gauss_weight*(DN(1,0)*crLHS339 - crLHS256*crLHS26 + crLHS256*crLHS79 - crLHS257 + crLHS338*crLHS78);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS80 + DN(3,1)*crLHS82 + DN(3,2)*crLHS84 + crLHS308 + crLHS341);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS92 + DN(3,1)*crLHS94 + DN(3,2)*crLHS97 + crLHS319);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS99 + DN(3,1)*crLHS101 + DN(3,2)*crLHS103 + crLHS327);
rLHS(12,11)+=gauss_weight*(DN(2,0)*crLHS339 + crLHS107*crLHS338 - crLHS26*crLHS312 + crLHS312*crLHS79 - crLHS313);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS108 + DN(3,1)*crLHS110 + DN(3,2)*crLHS112 + crLHS12*crLHS342 + crLHS344);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS120 + DN(3,1)*crLHS122 + DN(3,2)*crLHS125 + crLHS345);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS129 + DN(3,2)*crLHS131 + crLHS346);
rLHS(12,15)+=DN(3,0)*crLHS347;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS136 + DN(3,2)*crLHS137 + crLHS126);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS31 + DN(3,1)*crLHS138 + DN(3,2)*crLHS140 + crLHS181 + crLHS348);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS39 + DN(3,1)*crLHS142 + DN(3,2)*crLHS144 + crLHS213);
rLHS(13,3)+=gauss_weight*(DN(0,1)*crLHS339 - crLHS187*crLHS26 + crLHS187*crLHS79 - crLHS188 + crLHS349*crLHS45);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS53 + DN(3,1)*crLHS147 + DN(3,2)*crLHS148 + crLHS254);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS65 + DN(3,1)*crLHS150 + DN(3,2)*crLHS152 + crLHS269 + crLHS350);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS72 + DN(3,1)*crLHS155 + DN(3,2)*crLHS157 + crLHS282);
rLHS(13,7)+=gauss_weight*(DN(1,1)*crLHS339 - crLHS26*crLHS272 + crLHS272*crLHS79 - crLHS273 + crLHS349*crLHS78);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS82 + DN(3,1)*crLHS161 + DN(3,2)*crLHS162 + crLHS310);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS94 + DN(3,1)*crLHS164 + DN(3,2)*crLHS166 + crLHS320 + crLHS351);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS101 + DN(3,1)*crLHS169 + DN(3,2)*crLHS171 + crLHS328);
rLHS(13,11)+=gauss_weight*(DN(2,1)*crLHS339 + crLHS107*crLHS349 - crLHS26*crLHS323 + crLHS323*crLHS79 - crLHS324);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS110 + DN(3,1)*crLHS175 + DN(3,2)*crLHS176 + crLHS345);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS122 + DN(3,1)*crLHS178 + DN(3,2)*crLHS180 + crLHS12*crLHS352 + crLHS344);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS129 + DN(3,1)*crLHS183 + DN(3,2)*crLHS185 + crLHS353);
rLHS(13,15)+=DN(3,1)*crLHS347;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS137 + DN(3,2)*crLHS189 + crLHS132);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS34 + DN(3,1)*crLHS140 + DN(3,2)*crLHS190 + crLHS186);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS41 + DN(3,1)*crLHS144 + DN(3,2)*crLHS191 + crLHS215 + crLHS348);
rLHS(14,3)+=gauss_weight*(DN(0,2)*crLHS339 + DN(3,2)*crLHS46 - crLHS216*crLHS26 + crLHS216*crLHS79 - crLHS217);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS55 + DN(3,1)*crLHS148 + DN(3,2)*crLHS193 + crLHS255);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS68 + DN(3,1)*crLHS152 + DN(3,2)*crLHS196 + crLHS271);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS74 + DN(3,1)*crLHS157 + DN(3,2)*crLHS198 + crLHS283 + crLHS350);
rLHS(14,7)+=gauss_weight*(DN(1,2)*crLHS339 + DN(3,2)*crLHS242 - crLHS26*crLHS284 + crLHS284*crLHS79 - crLHS285);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS84 + DN(3,1)*crLHS162 + DN(3,2)*crLHS202 + crLHS311);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS97 + DN(3,1)*crLHS166 + DN(3,2)*crLHS204 + crLHS322);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS103 + DN(3,1)*crLHS171 + DN(3,2)*crLHS206 + crLHS329 + crLHS351);
rLHS(14,11)+=gauss_weight*(DN(2,2)*crLHS339 + DN(3,2)*crLHS305 - crLHS26*crLHS330 + crLHS330*crLHS79 - crLHS331);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS112 + DN(3,1)*crLHS176 + DN(3,2)*crLHS210 + crLHS346);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS125 + DN(3,1)*crLHS180 + DN(3,2)*crLHS212 + crLHS353);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS131 + DN(3,1)*crLHS185 + DN(3,2)*crLHS214 + crLHS12*crLHS354 + crLHS344);
rLHS(14,15)+=DN(3,2)*crLHS347;
rLHS(15,0)+=gauss_weight*(DN(3,0)*crLHS218 + crLHS133);
rLHS(15,1)+=gauss_weight*(DN(3,1)*crLHS218 + crLHS187);
rLHS(15,2)+=gauss_weight*(DN(3,2)*crLHS218 + crLHS216);
rLHS(15,3)+=crLHS229;
rLHS(15,4)+=gauss_weight*(DN(3,0)*crLHS221 + crLHS256);
rLHS(15,5)+=gauss_weight*(DN(3,1)*crLHS221 + crLHS272);
rLHS(15,6)+=gauss_weight*(DN(3,2)*crLHS221 + crLHS284);
rLHS(15,7)+=crLHS291;
rLHS(15,8)+=gauss_weight*(DN(3,0)*crLHS226 + crLHS312);
rLHS(15,9)+=gauss_weight*(DN(3,1)*crLHS226 + crLHS323);
rLHS(15,10)+=gauss_weight*(DN(3,2)*crLHS226 + crLHS330);
rLHS(15,11)+=crLHS333;
rLHS(15,12)+=DN(3,0)*crLHS355;
rLHS(15,13)+=DN(3,1)*crLHS355;
rLHS(15,14)+=DN(3,2)*crLHS355;
rLHS(15,15)+=gauss_weight*(crLHS20*crLHS342 + crLHS20*crLHS352 + crLHS20*crLHS354 + crLHS220*crLHS343);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double sigma = rData.Resistance;

    const double h = rData.ElementSize;
    const array_1d<double,3>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3>& v = rData.Velocity;
    const BoundedMatrix<double,2,3>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,3>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = v - vmesh;
    const BoundedMatrix<double,2,3>& f = rData.BodyForce;
    const array_1d<double,3>& p = rData.Pressure;
    const array_1d<double,3>& pn = rData.Pressure_OldStep1;
    const array_1d<double,3>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,3>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;
    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crRHS1 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crRHS2 = N[0]*sigma;
const double crRHS3 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crRHS4 = crRHS3*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crRHS5 = crRHS3*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crRHS6 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS7 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS9 = crRHS3*(crRHS6*crRHS7 + crRHS8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crRHS10 = sigma*stab_c3;
const double crRHS11 = crRHS3*stab_c2*sqrt(pow(crRHS7, 2) + pow(crRHS8, 2));
const double crRHS12 = 1.0/crRHS3;
const double crRHS13 = crRHS12*(crRHS7*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crRHS8*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crRHS14 = crRHS12*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crRHS15 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS16 = crRHS15 + crRHS6;
const double crRHS17 = (mu + (crRHS10 + crRHS11*h)/stab_c1)*(crRHS13 + crRHS14 + crRHS16);
const double crRHS18 = pow(h, -2);
const double crRHS19 = 1.0/(crRHS10*crRHS18 + crRHS11/h + crRHS18*mu*stab_c1 + crRHS3*dyn_tau/dt);
const double crRHS20 = crRHS19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS4 + crRHS5 + crRHS9);
const double crRHS21 = crRHS3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crRHS22 = N[0]*crRHS21;
const double crRHS23 = crRHS3*(DN(0,0)*crRHS7 + DN(0,1)*crRHS8);
const double crRHS24 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crRHS25 = crRHS3*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS26 = crRHS3*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crRHS27 = crRHS3*(crRHS15*crRHS8 + crRHS7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crRHS28 = crRHS19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS25 + crRHS26 + crRHS27);
const double crRHS29 = N[1]*sigma;
const double crRHS30 = N[1]*crRHS21;
const double crRHS31 = crRHS3*(DN(1,0)*crRHS7 + DN(1,1)*crRHS8);
const double crRHS32 = N[2]*sigma;
const double crRHS33 = N[2]*crRHS21;
const double crRHS34 = crRHS3*(DN(2,0)*crRHS7 + DN(2,1)*crRHS8);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS17 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS4 + N[0]*crRHS5 + N[0]*crRHS9 + crRHS1*crRHS2 - crRHS2*crRHS20 + crRHS20*crRHS22 + crRHS20*crRHS23);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS17 + DN(0,1)*stress[1] - N[0]*crRHS25 + N[0]*crRHS26 + N[0]*crRHS27 + crRHS2*crRHS24 - crRHS2*crRHS28 + crRHS22*crRHS28 + crRHS23*crRHS28);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS20 + DN(0,1)*crRHS28 + N[0]*crRHS13 + N[0]*crRHS14 + N[0]*crRHS16);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS17 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS4 + N[1]*crRHS5 + N[1]*crRHS9 + crRHS1*crRHS29 - crRHS20*crRHS29 + crRHS20*crRHS30 + crRHS20*crRHS31);
rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS17 + DN(1,1)*stress[1] - N[1]*crRHS25 + N[1]*crRHS26 + N[1]*crRHS27 + crRHS24*crRHS29 - crRHS28*crRHS29 + crRHS28*crRHS30 + crRHS28*crRHS31);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS20 + DN(1,1)*crRHS28 + N[1]*crRHS13 + N[1]*crRHS14 + N[1]*crRHS16);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS17 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS4 + N[2]*crRHS5 + N[2]*crRHS9 + crRHS1*crRHS32 - crRHS20*crRHS32 + crRHS20*crRHS33 + crRHS20*crRHS34);
rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS17 + DN(2,1)*stress[1] - N[2]*crRHS25 + N[2]*crRHS26 + N[2]*crRHS27 + crRHS24*crRHS32 - crRHS28*crRHS32 + crRHS28*crRHS33 + crRHS28*crRHS34);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS20 + DN(2,1)*crRHS28 + N[2]*crRHS13 + N[2]*crRHS14 + N[2]*crRHS16);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double sigma = rData.Resistance;

    const double h = rData.ElementSize;
    const array_1d<double,4>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,3,4>& v = rData.Velocity;
    const BoundedMatrix<double,3,4>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,3,4>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = v - vmesh;
    const BoundedMatrix<double,3,4>& f = rData.BodyForce;
    const array_1d<double,4>& p = rData.Pressure;
    const array_1d<double,4>& pn = rData.Pressure_OldStep1;
    const array_1d<double,4>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,6>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;
    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crRHS1 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crRHS2 = N[0]*sigma;
const double crRHS3 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crRHS4 = crRHS3*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS5 = crRHS3*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crRHS6 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS7 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS9 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS10 = crRHS3*(crRHS6*crRHS7 + crRHS8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crRHS9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crRHS11 = sigma*stab_c3;
const double crRHS12 = crRHS3*stab_c2*sqrt(pow(crRHS7, 2) + pow(crRHS8, 2) + pow(crRHS9, 2));
const double crRHS13 = 1.0/crRHS3;
const double crRHS14 = crRHS13*(crRHS7*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crRHS8*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]) + crRHS9*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]));
const double crRHS15 = crRHS13*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crRHS16 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS17 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS18 = crRHS16 + crRHS17 + crRHS6;
const double crRHS19 = (mu + (crRHS11 + crRHS12*h)/stab_c1)*(crRHS14 + crRHS15 + crRHS18);
const double crRHS20 = pow(h, -2);
const double crRHS21 = 1.0/(crRHS11*crRHS20 + crRHS12/h + crRHS20*mu*stab_c1 + crRHS3*dyn_tau/dt);
const double crRHS22 = crRHS21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crRHS10 - crRHS4 + crRHS5);
const double crRHS23 = crRHS3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crRHS24 = N[0]*crRHS23;
const double crRHS25 = crRHS3*(DN(0,0)*crRHS7 + DN(0,1)*crRHS8 + DN(0,2)*crRHS9);
const double crRHS26 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crRHS27 = crRHS3*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS28 = crRHS3*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crRHS29 = crRHS3*(crRHS16*crRHS8 + crRHS7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crRHS9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crRHS30 = crRHS21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS27 + crRHS28 + crRHS29);
const double crRHS31 = N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crRHS32 = crRHS3*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS33 = crRHS3*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crRHS34 = crRHS3*(crRHS17*crRHS9 + crRHS7*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crRHS8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crRHS35 = crRHS21*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS32 + crRHS33 + crRHS34);
const double crRHS36 = N[1]*sigma;
const double crRHS37 = N[1]*crRHS23;
const double crRHS38 = crRHS3*(DN(1,0)*crRHS7 + DN(1,1)*crRHS8 + DN(1,2)*crRHS9);
const double crRHS39 = N[2]*sigma;
const double crRHS40 = N[2]*crRHS23;
const double crRHS41 = crRHS3*(DN(2,0)*crRHS7 + DN(2,1)*crRHS8 + DN(2,2)*crRHS9);
const double crRHS42 = N[3]*sigma;
const double crRHS43 = N[3]*crRHS23;
const double crRHS44 = crRHS3*(DN(3,0)*crRHS7 + DN(3,1)*crRHS8 + DN(3,2)*crRHS9);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS19 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] + N[0]*crRHS10 - N[0]*crRHS4 + N[0]*crRHS5 + crRHS1*crRHS2 - crRHS2*crRHS22 + crRHS22*crRHS24 + crRHS22*crRHS25);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS19 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS27 + N[0]*crRHS28 + N[0]*crRHS29 + crRHS2*crRHS26 - crRHS2*crRHS30 + crRHS24*crRHS30 + crRHS25*crRHS30);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS19 + DN(0,2)*stress[2] - N[0]*crRHS32 + N[0]*crRHS33 + N[0]*crRHS34 + crRHS2*crRHS31 - crRHS2*crRHS35 + crRHS24*crRHS35 + crRHS25*crRHS35);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS22 + DN(0,1)*crRHS30 + DN(0,2)*crRHS35 + N[0]*crRHS14 + N[0]*crRHS15 + N[0]*crRHS18);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS19 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] + N[1]*crRHS10 - N[1]*crRHS4 + N[1]*crRHS5 + crRHS1*crRHS36 - crRHS22*crRHS36 + crRHS22*crRHS37 + crRHS22*crRHS38);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS19 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS27 + N[1]*crRHS28 + N[1]*crRHS29 + crRHS26*crRHS36 - crRHS30*crRHS36 + crRHS30*crRHS37 + crRHS30*crRHS38);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS19 + DN(1,2)*stress[2] - N[1]*crRHS32 + N[1]*crRHS33 + N[1]*crRHS34 + crRHS31*crRHS36 - crRHS35*crRHS36 + crRHS35*crRHS37 + crRHS35*crRHS38);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS22 + DN(1,1)*crRHS30 + DN(1,2)*crRHS35 + N[1]*crRHS14 + N[1]*crRHS15 + N[1]*crRHS18);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS19 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] + N[2]*crRHS10 - N[2]*crRHS4 + N[2]*crRHS5 + crRHS1*crRHS39 - crRHS22*crRHS39 + crRHS22*crRHS40 + crRHS22*crRHS41);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS19 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS27 + N[2]*crRHS28 + N[2]*crRHS29 + crRHS26*crRHS39 - crRHS30*crRHS39 + crRHS30*crRHS40 + crRHS30*crRHS41);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS19 + DN(2,2)*stress[2] - N[2]*crRHS32 + N[2]*crRHS33 + N[2]*crRHS34 + crRHS31*crRHS39 - crRHS35*crRHS39 + crRHS35*crRHS40 + crRHS35*crRHS41);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS22 + DN(2,1)*crRHS30 + DN(2,2)*crRHS35 + N[2]*crRHS14 + N[2]*crRHS15 + N[2]*crRHS18);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS19 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] + N[3]*crRHS10 - N[3]*crRHS4 + N[3]*crRHS5 + crRHS1*crRHS42 - crRHS22*crRHS42 + crRHS22*crRHS43 + crRHS22*crRHS44);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS19 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS27 + N[3]*crRHS28 + N[3]*crRHS29 + crRHS26*crRHS42 - crRHS30*crRHS42 + crRHS30*crRHS43 + crRHS30*crRHS44);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS19 + DN(3,2)*stress[2] - N[3]*crRHS32 + N[3]*crRHS33 + N[3]*crRHS34 + crRHS31*crRHS42 - crRHS35*crRHS42 + crRHS35*crRHS43 + crRHS35*crRHS44);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS22 + DN(3,1)*crRHS30 + DN(3,2)*crRHS35 + N[3]*crRHS14 + N[3]*crRHS15 + N[3]*crRHS18);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void WeaklyCompressibleNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void WeaklyCompressibleNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> >;
template class WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> >;

}