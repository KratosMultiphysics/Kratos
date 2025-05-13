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
#include "data_containers/weakly_compressible_navier_stokes/weakly_compressible_navier_stokes_data.h"

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

    // Contribution coming from the shear stress operator
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
const double crLHS3 = DN(0,0)*DN(0,0);
const double crLHS4 = sigma*stab_c3;
const double crLHS5 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crLHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS8 = crLHS5*stab_c2*sqrt(crLHS6*crLHS6 + crLHS7*crLHS7);
const double crLHS9 = mu + (crLHS4 + crLHS8*h)*1.0/stab_c1;
const double crLHS10 = N[0]*N[0];
const double crLHS11 = crLHS5*(DN(0,0)*crLHS6 + DN(0,1)*crLHS7);
const double crLHS12 = bdf0*crLHS5;
const double crLHS13 = N[0]*sigma;
const double crLHS14 = N[0]*crLHS12;
const double crLHS15 = crLHS11 + crLHS13 + crLHS14;
const double crLHS16 = 1.0/h;
const double crLHS17 = 1.0/(crLHS16*crLHS4 + crLHS16*crLHS8 + crLHS5*dyn_tau*1.0/dt + mu*stab_c1*1.0/(h*h));
const double crLHS18 = 1.0*crLHS11;
const double crLHS19 = crLHS17*crLHS18;
const double crLHS20 = 1.0*crLHS13;
const double crLHS21 = crLHS17*crLHS20;
const double crLHS22 = 1.0*crLHS17;
const double crLHS23 = crLHS15*crLHS22;
const double crLHS24 = crLHS5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crLHS25 = N[0]*crLHS24;
const double crLHS26 = N[0]*crLHS11 + crLHS10*crLHS12 + crLHS10*sigma + crLHS15*crLHS19 - crLHS15*crLHS21 + crLHS23*crLHS25;
const double crLHS27 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS28 = C(1,2)*DN(0,1);
const double crLHS29 = C(2,2)*DN(0,0) + crLHS28;
const double crLHS30 = DN(0,0)*crLHS9;
const double crLHS31 = DN(0,1)*crLHS30;
const double crLHS32 = bdf0*1.0/crLHS5/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crLHS33 = N[0]*crLHS32;
const double crLHS34 = crLHS22*crLHS24;
const double crLHS35 = gauss_weight*(N[0]*crLHS34 - N[0] + crLHS19 - crLHS21 + crLHS33*crLHS9);
const double crLHS36 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS37 = C(0,2)*DN(1,0);
const double crLHS38 = C(2,2)*DN(1,1) + crLHS37;
const double crLHS39 = N[1]*crLHS13 + N[1]*crLHS14;
const double crLHS40 = DN(1,0)*crLHS30 + crLHS39;
const double crLHS41 = crLHS5*(DN(1,0)*crLHS6 + DN(1,1)*crLHS7);
const double crLHS42 = N[1]*sigma;
const double crLHS43 = N[1]*crLHS12;
const double crLHS44 = crLHS41 + crLHS42 + crLHS43;
const double crLHS45 = crLHS22*crLHS44;
const double crLHS46 = N[0]*crLHS41 + crLHS19*crLHS44 - crLHS21*crLHS44 + crLHS25*crLHS45;
const double crLHS47 = C(0,1)*DN(1,1) + crLHS37;
const double crLHS48 = C(1,2)*DN(1,1);
const double crLHS49 = C(2,2)*DN(1,0) + crLHS48;
const double crLHS50 = DN(1,1)*crLHS30;
const double crLHS51 = DN(0,0)*N[1];
const double crLHS52 = crLHS32*crLHS9;
const double crLHS53 = DN(1,0)*N[0];
const double crLHS54 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS55 = C(0,2)*DN(2,0);
const double crLHS56 = C(2,2)*DN(2,1) + crLHS55;
const double crLHS57 = N[2]*crLHS13 + N[2]*crLHS14;
const double crLHS58 = DN(2,0)*crLHS30 + crLHS57;
const double crLHS59 = crLHS5*(DN(2,0)*crLHS6 + DN(2,1)*crLHS7);
const double crLHS60 = N[2]*sigma;
const double crLHS61 = N[2]*crLHS12;
const double crLHS62 = crLHS59 + crLHS60 + crLHS61;
const double crLHS63 = crLHS22*crLHS62;
const double crLHS64 = N[0]*crLHS59 + crLHS19*crLHS62 - crLHS21*crLHS62 + crLHS25*crLHS63;
const double crLHS65 = C(0,1)*DN(2,1) + crLHS55;
const double crLHS66 = C(1,2)*DN(2,1);
const double crLHS67 = C(2,2)*DN(2,0) + crLHS66;
const double crLHS68 = DN(2,1)*crLHS30;
const double crLHS69 = DN(0,0)*N[2];
const double crLHS70 = DN(2,0)*N[0];
const double crLHS71 = C(0,1)*DN(0,0) + crLHS28;
const double crLHS72 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS73 = DN(0,1)*DN(0,1);
const double crLHS74 = C(0,1)*DN(1,0) + crLHS48;
const double crLHS75 = DN(0,1)*crLHS9;
const double crLHS76 = DN(1,0)*crLHS75;
const double crLHS77 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS78 = DN(1,1)*crLHS75 + crLHS39;
const double crLHS79 = DN(0,1)*N[1];
const double crLHS80 = DN(1,1)*N[0];
const double crLHS81 = C(0,1)*DN(2,0) + crLHS66;
const double crLHS82 = DN(2,0)*crLHS75;
const double crLHS83 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS84 = DN(2,1)*crLHS75 + crLHS57;
const double crLHS85 = DN(0,1)*N[2];
const double crLHS86 = DN(2,1)*N[0];
const double crLHS87 = gauss_weight*(N[0] + crLHS17*(1.0*crLHS14 + crLHS18 + crLHS20));
const double crLHS88 = DN(0,0)*crLHS22;
const double crLHS89 = DN(0,1)*crLHS22;
const double crLHS90 = gauss_weight*(DN(1,0)*crLHS88 + DN(1,1)*crLHS89 + N[1]*crLHS33);
const double crLHS91 = gauss_weight*(DN(2,0)*crLHS88 + DN(2,1)*crLHS89 + N[2]*crLHS33);
const double crLHS92 = crLHS22*crLHS41;
const double crLHS93 = crLHS22*crLHS42;
const double crLHS94 = N[1]*crLHS24;
const double crLHS95 = N[1]*crLHS11 + crLHS15*crLHS92 - crLHS15*crLHS93 + crLHS23*crLHS94;
const double crLHS96 = DN(1,0)*DN(1,0);
const double crLHS97 = N[1]*N[1];
const double crLHS98 = N[1]*crLHS41 + crLHS12*crLHS97 + crLHS41*crLHS45 - crLHS42*crLHS45 + crLHS45*crLHS94 + crLHS97*sigma;
const double crLHS99 = DN(1,0)*crLHS9;
const double crLHS100 = DN(1,1)*crLHS99;
const double crLHS101 = N[1]*crLHS32;
const double crLHS102 = gauss_weight*(N[1]*crLHS34 - N[1] + crLHS101*crLHS9 + crLHS92 - crLHS93);
const double crLHS103 = N[2]*crLHS42 + N[2]*crLHS43;
const double crLHS104 = DN(2,0)*crLHS99 + crLHS103;
const double crLHS105 = N[1]*crLHS59 + crLHS62*crLHS92 - crLHS62*crLHS93 + crLHS63*crLHS94;
const double crLHS106 = DN(2,1)*crLHS99;
const double crLHS107 = DN(1,0)*N[2];
const double crLHS108 = DN(2,0)*N[1];
const double crLHS109 = DN(1,1)*DN(1,1);
const double crLHS110 = DN(1,1)*crLHS9;
const double crLHS111 = DN(2,0)*crLHS110;
const double crLHS112 = DN(2,1)*crLHS110 + crLHS103;
const double crLHS113 = DN(1,1)*N[2];
const double crLHS114 = DN(2,1)*N[1];
const double crLHS115 = gauss_weight*(N[1] + crLHS17*(1.0*crLHS41 + 1.0*crLHS42 + 1.0*crLHS43));
const double crLHS116 = DN(1,0)*crLHS22;
const double crLHS117 = DN(1,1)*crLHS22;
const double crLHS118 = gauss_weight*(DN(2,0)*crLHS116 + DN(2,1)*crLHS117 + N[2]*crLHS101);
const double crLHS119 = N[2]*crLHS24;
const double crLHS120 = N[2]*crLHS11 + crLHS119*crLHS23 + crLHS23*crLHS59 - crLHS23*crLHS60;
const double crLHS121 = N[2]*crLHS41 + crLHS119*crLHS45 + crLHS45*crLHS59 - crLHS45*crLHS60;
const double crLHS122 = DN(2,0)*DN(2,0);
const double crLHS123 = N[2]*N[2];
const double crLHS124 = N[2]*crLHS59 + crLHS119*crLHS63 + crLHS12*crLHS123 + crLHS123*sigma + crLHS59*crLHS63 - crLHS60*crLHS63;
const double crLHS125 = DN(2,0)*DN(2,1)*crLHS9;
const double crLHS126 = gauss_weight*(N[2]*crLHS34 + N[2]*crLHS52 - N[2] + crLHS22*crLHS59 - crLHS22*crLHS60);
const double crLHS127 = DN(2,1)*DN(2,1);
const double crLHS128 = gauss_weight*(N[2] + crLHS17*(1.0*crLHS59 + 1.0*crLHS60 + 1.0*crLHS61));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS26 + crLHS3*crLHS9);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS27 + DN(0,1)*crLHS29 + crLHS31);
rLHS(0,2)+=DN(0,0)*crLHS35;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS36 + DN(0,1)*crLHS38 + crLHS40 + crLHS46);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS47 + DN(0,1)*crLHS49 + crLHS50);
rLHS(0,5)+=gauss_weight*(DN(1,0)*crLHS19 - DN(1,0)*crLHS21 + crLHS34*crLHS53 + crLHS51*crLHS52 - crLHS51);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS56 + crLHS58 + crLHS64);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS65 + DN(0,1)*crLHS67 + crLHS68);
rLHS(0,8)+=gauss_weight*(DN(2,0)*crLHS19 - DN(2,0)*crLHS21 + crLHS34*crLHS70 + crLHS52*crLHS69 - crLHS69);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS71 + crLHS31);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS29 + DN(0,1)*crLHS72 + crLHS26 + crLHS73*crLHS9);
rLHS(1,2)+=DN(0,1)*crLHS35;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS74 + crLHS76);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS49 + DN(0,1)*crLHS77 + crLHS46 + crLHS78);
rLHS(1,5)+=gauss_weight*(DN(1,1)*crLHS19 - DN(1,1)*crLHS21 + crLHS34*crLHS80 + crLHS52*crLHS79 - crLHS79);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS56 + DN(0,1)*crLHS81 + crLHS82);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS67 + DN(0,1)*crLHS83 + crLHS64 + crLHS84);
rLHS(1,8)+=gauss_weight*(DN(2,1)*crLHS19 - DN(2,1)*crLHS21 + crLHS34*crLHS86 + crLHS52*crLHS85 - crLHS85);
rLHS(2,0)+=DN(0,0)*crLHS87;
rLHS(2,1)+=DN(0,1)*crLHS87;
rLHS(2,2)+=gauss_weight*(crLHS10*crLHS32 + crLHS22*crLHS3 + crLHS22*crLHS73);
rLHS(2,3)+=gauss_weight*(DN(0,0)*crLHS45 + crLHS53);
rLHS(2,4)+=gauss_weight*(DN(0,1)*crLHS45 + crLHS80);
rLHS(2,5)+=crLHS90;
rLHS(2,6)+=gauss_weight*(crLHS62*crLHS88 + crLHS70);
rLHS(2,7)+=gauss_weight*(crLHS62*crLHS89 + crLHS86);
rLHS(2,8)+=crLHS91;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS40 + crLHS95);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS27 + DN(1,1)*crLHS29 + crLHS76);
rLHS(3,2)+=gauss_weight*(crLHS34*crLHS51 + crLHS41*crLHS88 - crLHS42*crLHS88 + crLHS52*crLHS53 - crLHS53);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS36 + DN(1,1)*crLHS38 + crLHS9*crLHS96 + crLHS98);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS47 + DN(1,1)*crLHS49 + crLHS100);
rLHS(3,5)+=DN(1,0)*crLHS102;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS56 + crLHS104 + crLHS105);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS65 + DN(1,1)*crLHS67 + crLHS106);
rLHS(3,8)+=gauss_weight*(DN(2,0)*crLHS92 - DN(2,0)*crLHS93 + crLHS107*crLHS52 - crLHS107 + crLHS108*crLHS34);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS71 + crLHS50);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS29 + DN(1,1)*crLHS72 + crLHS78 + crLHS95);
rLHS(4,2)+=gauss_weight*(crLHS34*crLHS79 + crLHS41*crLHS89 - crLHS42*crLHS89 + crLHS52*crLHS80 - crLHS80);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS74 + crLHS100);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS49 + DN(1,1)*crLHS77 + crLHS109*crLHS9 + crLHS98);
rLHS(4,5)+=DN(1,1)*crLHS102;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS56 + DN(1,1)*crLHS81 + crLHS111);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS67 + DN(1,1)*crLHS83 + crLHS105 + crLHS112);
rLHS(4,8)+=gauss_weight*(DN(2,1)*crLHS92 - DN(2,1)*crLHS93 + crLHS113*crLHS52 - crLHS113 + crLHS114*crLHS34);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS23 + crLHS51);
rLHS(5,1)+=gauss_weight*(DN(1,1)*crLHS23 + crLHS79);
rLHS(5,2)+=crLHS90;
rLHS(5,3)+=DN(1,0)*crLHS115;
rLHS(5,4)+=DN(1,1)*crLHS115;
rLHS(5,5)+=gauss_weight*(crLHS109*crLHS22 + crLHS22*crLHS96 + crLHS32*crLHS97);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS63 + crLHS108);
rLHS(5,7)+=gauss_weight*(DN(1,1)*crLHS63 + crLHS114);
rLHS(5,8)+=crLHS118;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS120 + crLHS58);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS27 + DN(2,1)*crLHS29 + crLHS82);
rLHS(6,2)+=gauss_weight*(crLHS34*crLHS69 + crLHS52*crLHS70 + crLHS59*crLHS88 - crLHS60*crLHS88 - crLHS70);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS36 + DN(2,1)*crLHS38 + crLHS104 + crLHS121);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS47 + DN(2,1)*crLHS49 + crLHS111);
rLHS(6,5)+=gauss_weight*(crLHS107*crLHS34 + crLHS108*crLHS52 - crLHS108 + crLHS116*crLHS59 - crLHS116*crLHS60);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS56 + crLHS122*crLHS9 + crLHS124);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS65 + DN(2,1)*crLHS67 + crLHS125);
rLHS(6,8)+=DN(2,0)*crLHS126;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS71 + crLHS68);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS29 + DN(2,1)*crLHS72 + crLHS120 + crLHS84);
rLHS(7,2)+=gauss_weight*(crLHS34*crLHS85 + crLHS52*crLHS86 + crLHS59*crLHS89 - crLHS60*crLHS89 - crLHS86);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS74 + crLHS106);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS49 + DN(2,1)*crLHS77 + crLHS112 + crLHS121);
rLHS(7,5)+=gauss_weight*(crLHS113*crLHS34 + crLHS114*crLHS52 - crLHS114 + crLHS117*crLHS59 - crLHS117*crLHS60);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS56 + DN(2,1)*crLHS81 + crLHS125);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS67 + DN(2,1)*crLHS83 + crLHS124 + crLHS127*crLHS9);
rLHS(7,8)+=DN(2,1)*crLHS126;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS23 + crLHS69);
rLHS(8,1)+=gauss_weight*(DN(2,1)*crLHS23 + crLHS85);
rLHS(8,2)+=crLHS91;
rLHS(8,3)+=gauss_weight*(DN(2,0)*crLHS45 + crLHS107);
rLHS(8,4)+=gauss_weight*(DN(2,1)*crLHS45 + crLHS113);
rLHS(8,5)+=crLHS118;
rLHS(8,6)+=DN(2,0)*crLHS128;
rLHS(8,7)+=DN(2,1)*crLHS128;
rLHS(8,8)+=gauss_weight*(crLHS122*crLHS22 + crLHS123*crLHS32 + crLHS127*crLHS22);

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
const double crLHS5 = DN(0,0)*DN(0,0);
const double crLHS6 = sigma*stab_c3;
const double crLHS7 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS11 = crLHS7*stab_c2*sqrt(crLHS10*crLHS10 + crLHS8*crLHS8 + crLHS9*crLHS9);
const double crLHS12 = mu + (crLHS11*h + crLHS6)*1.0/stab_c1;
const double crLHS13 = N[0]*N[0];
const double crLHS14 = crLHS7*(DN(0,0)*crLHS8 + DN(0,1)*crLHS9 + DN(0,2)*crLHS10);
const double crLHS15 = bdf0*crLHS7;
const double crLHS16 = N[0]*sigma;
const double crLHS17 = N[0]*crLHS15;
const double crLHS18 = crLHS14 + crLHS16 + crLHS17;
const double crLHS19 = 1.0/h;
const double crLHS20 = 1.0/(crLHS11*crLHS19 + crLHS19*crLHS6 + crLHS7*dyn_tau*1.0/dt + mu*stab_c1*1.0/(h*h));
const double crLHS21 = 1.0*crLHS14;
const double crLHS22 = crLHS20*crLHS21;
const double crLHS23 = 1.0*crLHS16;
const double crLHS24 = crLHS20*crLHS23;
const double crLHS25 = 1.0*crLHS20;
const double crLHS26 = crLHS18*crLHS25;
const double crLHS27 = crLHS7*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crLHS28 = N[0]*crLHS27;
const double crLHS29 = N[0]*crLHS14 + crLHS13*crLHS15 + crLHS13*sigma + crLHS18*crLHS22 - crLHS18*crLHS24 + crLHS26*crLHS28;
const double crLHS30 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS31 = C(1,3)*DN(0,1);
const double crLHS32 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS31;
const double crLHS33 = C(3,5)*DN(0,0);
const double crLHS34 = C(4,5)*DN(0,2);
const double crLHS35 = C(1,5)*DN(0,1) + crLHS33 + crLHS34;
const double crLHS36 = DN(0,0)*crLHS12;
const double crLHS37 = DN(0,1)*crLHS36;
const double crLHS38 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS39 = C(3,4)*DN(0,1);
const double crLHS40 = C(2,3)*DN(0,2) + crLHS33 + crLHS39;
const double crLHS41 = C(2,5)*DN(0,2);
const double crLHS42 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS41;
const double crLHS43 = DN(0,2)*crLHS36;
const double crLHS44 = bdf0*1.0/crLHS7/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crLHS45 = N[0]*crLHS44;
const double crLHS46 = crLHS25*crLHS27;
const double crLHS47 = gauss_weight*(N[0]*crLHS46 - N[0] + crLHS12*crLHS45 + crLHS22 - crLHS24);
const double crLHS48 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS49 = C(0,3)*DN(1,0);
const double crLHS50 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS49;
const double crLHS51 = C(0,5)*DN(1,0);
const double crLHS52 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS51;
const double crLHS53 = N[1]*crLHS16 + N[1]*crLHS17;
const double crLHS54 = DN(1,0)*crLHS36 + crLHS53;
const double crLHS55 = crLHS7*(DN(1,0)*crLHS8 + DN(1,1)*crLHS9 + DN(1,2)*crLHS10);
const double crLHS56 = N[1]*sigma;
const double crLHS57 = N[1]*crLHS15;
const double crLHS58 = crLHS55 + crLHS56 + crLHS57;
const double crLHS59 = crLHS25*crLHS58;
const double crLHS60 = N[0]*crLHS55 + crLHS22*crLHS58 - crLHS24*crLHS58 + crLHS28*crLHS59;
const double crLHS61 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS49;
const double crLHS62 = C(1,3)*DN(1,1);
const double crLHS63 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS62;
const double crLHS64 = C(3,5)*DN(1,0);
const double crLHS65 = C(4,5)*DN(1,2);
const double crLHS66 = C(1,5)*DN(1,1) + crLHS64 + crLHS65;
const double crLHS67 = DN(1,1)*crLHS36;
const double crLHS68 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS51;
const double crLHS69 = C(3,4)*DN(1,1);
const double crLHS70 = C(2,3)*DN(1,2) + crLHS64 + crLHS69;
const double crLHS71 = C(2,5)*DN(1,2);
const double crLHS72 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS71;
const double crLHS73 = DN(1,2)*crLHS36;
const double crLHS74 = DN(0,0)*N[1];
const double crLHS75 = crLHS12*crLHS44;
const double crLHS76 = DN(1,0)*N[0];
const double crLHS77 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS78 = C(0,3)*DN(2,0);
const double crLHS79 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS78;
const double crLHS80 = C(0,5)*DN(2,0);
const double crLHS81 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS80;
const double crLHS82 = N[2]*crLHS16 + N[2]*crLHS17;
const double crLHS83 = DN(2,0)*crLHS36 + crLHS82;
const double crLHS84 = crLHS7*(DN(2,0)*crLHS8 + DN(2,1)*crLHS9 + DN(2,2)*crLHS10);
const double crLHS85 = N[2]*sigma;
const double crLHS86 = N[2]*crLHS15;
const double crLHS87 = crLHS84 + crLHS85 + crLHS86;
const double crLHS88 = crLHS25*crLHS87;
const double crLHS89 = N[0]*crLHS84 + crLHS22*crLHS87 - crLHS24*crLHS87 + crLHS28*crLHS88;
const double crLHS90 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS78;
const double crLHS91 = C(1,3)*DN(2,1);
const double crLHS92 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS91;
const double crLHS93 = C(3,5)*DN(2,0);
const double crLHS94 = C(4,5)*DN(2,2);
const double crLHS95 = C(1,5)*DN(2,1) + crLHS93 + crLHS94;
const double crLHS96 = DN(2,1)*crLHS36;
const double crLHS97 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS80;
const double crLHS98 = C(3,4)*DN(2,1);
const double crLHS99 = C(2,3)*DN(2,2) + crLHS93 + crLHS98;
const double crLHS100 = C(2,5)*DN(2,2);
const double crLHS101 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS100;
const double crLHS102 = DN(2,2)*crLHS36;
const double crLHS103 = DN(0,0)*N[2];
const double crLHS104 = DN(2,0)*N[0];
const double crLHS105 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS106 = C(0,3)*DN(3,0);
const double crLHS107 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS106;
const double crLHS108 = C(0,5)*DN(3,0);
const double crLHS109 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS108;
const double crLHS110 = N[3]*crLHS16 + N[3]*crLHS17;
const double crLHS111 = DN(3,0)*crLHS36 + crLHS110;
const double crLHS112 = crLHS7*(DN(3,0)*crLHS8 + DN(3,1)*crLHS9 + DN(3,2)*crLHS10);
const double crLHS113 = N[3]*sigma;
const double crLHS114 = N[3]*crLHS15;
const double crLHS115 = crLHS112 + crLHS113 + crLHS114;
const double crLHS116 = crLHS115*crLHS25;
const double crLHS117 = N[0]*crLHS112 + crLHS115*crLHS22 - crLHS115*crLHS24 + crLHS116*crLHS28;
const double crLHS118 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS106;
const double crLHS119 = C(1,3)*DN(3,1);
const double crLHS120 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS119;
const double crLHS121 = C(3,5)*DN(3,0);
const double crLHS122 = C(4,5)*DN(3,2);
const double crLHS123 = C(1,5)*DN(3,1) + crLHS121 + crLHS122;
const double crLHS124 = DN(3,1)*crLHS36;
const double crLHS125 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS108;
const double crLHS126 = C(3,4)*DN(3,1);
const double crLHS127 = C(2,3)*DN(3,2) + crLHS121 + crLHS126;
const double crLHS128 = C(2,5)*DN(3,2);
const double crLHS129 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS128;
const double crLHS130 = DN(3,2)*crLHS36;
const double crLHS131 = DN(0,0)*N[3];
const double crLHS132 = DN(3,0)*N[0];
const double crLHS133 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS31;
const double crLHS134 = C(0,4)*DN(0,0) + crLHS34 + crLHS39;
const double crLHS135 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS136 = C(1,4)*DN(0,1);
const double crLHS137 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS136;
const double crLHS138 = DN(0,1)*DN(0,1);
const double crLHS139 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS136;
const double crLHS140 = C(2,4)*DN(0,2);
const double crLHS141 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS140;
const double crLHS142 = DN(0,1)*crLHS12;
const double crLHS143 = DN(0,2)*crLHS142;
const double crLHS144 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS62;
const double crLHS145 = C(0,4)*DN(1,0) + crLHS65 + crLHS69;
const double crLHS146 = DN(1,0)*crLHS142;
const double crLHS147 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS148 = C(1,4)*DN(1,1);
const double crLHS149 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS148;
const double crLHS150 = DN(1,1)*crLHS142;
const double crLHS151 = crLHS53 + crLHS60;
const double crLHS152 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS148;
const double crLHS153 = C(2,4)*DN(1,2);
const double crLHS154 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS153;
const double crLHS155 = DN(1,2)*crLHS142;
const double crLHS156 = DN(0,1)*N[1];
const double crLHS157 = DN(1,1)*N[0];
const double crLHS158 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS91;
const double crLHS159 = C(0,4)*DN(2,0) + crLHS94 + crLHS98;
const double crLHS160 = DN(2,0)*crLHS142;
const double crLHS161 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS162 = C(1,4)*DN(2,1);
const double crLHS163 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS162;
const double crLHS164 = DN(2,1)*crLHS142;
const double crLHS165 = crLHS82 + crLHS89;
const double crLHS166 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS162;
const double crLHS167 = C(2,4)*DN(2,2);
const double crLHS168 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS167;
const double crLHS169 = DN(2,2)*crLHS142;
const double crLHS170 = DN(0,1)*N[2];
const double crLHS171 = DN(2,1)*N[0];
const double crLHS172 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS119;
const double crLHS173 = C(0,4)*DN(3,0) + crLHS122 + crLHS126;
const double crLHS174 = DN(3,0)*crLHS142;
const double crLHS175 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS176 = C(1,4)*DN(3,1);
const double crLHS177 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS176;
const double crLHS178 = DN(3,1)*crLHS142;
const double crLHS179 = crLHS110 + crLHS117;
const double crLHS180 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS176;
const double crLHS181 = C(2,4)*DN(3,2);
const double crLHS182 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS181;
const double crLHS183 = DN(3,2)*crLHS142;
const double crLHS184 = DN(0,1)*N[3];
const double crLHS185 = DN(3,1)*N[0];
const double crLHS186 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS41;
const double crLHS187 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS140;
const double crLHS188 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS189 = DN(0,2)*DN(0,2);
const double crLHS190 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS71;
const double crLHS191 = DN(0,2)*crLHS12;
const double crLHS192 = DN(1,0)*crLHS191;
const double crLHS193 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS153;
const double crLHS194 = DN(1,1)*crLHS191;
const double crLHS195 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS196 = DN(1,2)*crLHS191;
const double crLHS197 = DN(0,2)*N[1];
const double crLHS198 = DN(1,2)*N[0];
const double crLHS199 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS100;
const double crLHS200 = DN(2,0)*crLHS191;
const double crLHS201 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS167;
const double crLHS202 = DN(2,1)*crLHS191;
const double crLHS203 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS204 = DN(2,2)*crLHS191;
const double crLHS205 = DN(0,2)*N[2];
const double crLHS206 = DN(2,2)*N[0];
const double crLHS207 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS128;
const double crLHS208 = DN(3,0)*crLHS191;
const double crLHS209 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS181;
const double crLHS210 = DN(3,1)*crLHS191;
const double crLHS211 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS212 = DN(3,2)*crLHS191;
const double crLHS213 = DN(0,2)*N[3];
const double crLHS214 = DN(3,2)*N[0];
const double crLHS215 = gauss_weight*(N[0] + crLHS20*(1.0*crLHS17 + crLHS21 + crLHS23));
const double crLHS216 = DN(0,0)*crLHS25;
const double crLHS217 = DN(0,1)*crLHS25;
const double crLHS218 = DN(0,2)*crLHS25;
const double crLHS219 = gauss_weight*(DN(1,0)*crLHS216 + DN(1,1)*crLHS217 + DN(1,2)*crLHS218 + N[1]*crLHS45);
const double crLHS220 = gauss_weight*(DN(2,0)*crLHS216 + DN(2,1)*crLHS217 + DN(2,2)*crLHS218 + N[2]*crLHS45);
const double crLHS221 = gauss_weight*(DN(3,0)*crLHS216 + DN(3,1)*crLHS217 + DN(3,2)*crLHS218 + N[3]*crLHS45);
const double crLHS222 = crLHS25*crLHS55;
const double crLHS223 = crLHS25*crLHS56;
const double crLHS224 = N[1]*crLHS27;
const double crLHS225 = N[1]*crLHS14 + crLHS18*crLHS222 - crLHS18*crLHS223 + crLHS224*crLHS26;
const double crLHS226 = DN(1,0)*DN(1,0);
const double crLHS227 = N[1]*N[1];
const double crLHS228 = N[1]*crLHS55 + crLHS15*crLHS227 + crLHS224*crLHS59 + crLHS227*sigma + crLHS55*crLHS59 - crLHS56*crLHS59;
const double crLHS229 = DN(1,0)*crLHS12;
const double crLHS230 = DN(1,1)*crLHS229;
const double crLHS231 = DN(1,2)*crLHS229;
const double crLHS232 = N[1]*crLHS44;
const double crLHS233 = gauss_weight*(N[1]*crLHS46 - N[1] + crLHS12*crLHS232 + crLHS222 - crLHS223);
const double crLHS234 = N[2]*crLHS56 + N[2]*crLHS57;
const double crLHS235 = DN(2,0)*crLHS229 + crLHS234;
const double crLHS236 = N[1]*crLHS84 + crLHS222*crLHS87 - crLHS223*crLHS87 + crLHS224*crLHS88;
const double crLHS237 = DN(2,1)*crLHS229;
const double crLHS238 = DN(2,2)*crLHS229;
const double crLHS239 = DN(1,0)*N[2];
const double crLHS240 = DN(2,0)*N[1];
const double crLHS241 = N[3]*crLHS56 + N[3]*crLHS57;
const double crLHS242 = DN(3,0)*crLHS229 + crLHS241;
const double crLHS243 = N[1]*crLHS112 + crLHS115*crLHS222 - crLHS115*crLHS223 + crLHS116*crLHS224;
const double crLHS244 = DN(3,1)*crLHS229;
const double crLHS245 = DN(3,2)*crLHS229;
const double crLHS246 = DN(1,0)*N[3];
const double crLHS247 = DN(3,0)*N[1];
const double crLHS248 = crLHS225 + crLHS53;
const double crLHS249 = DN(1,1)*DN(1,1);
const double crLHS250 = DN(1,1)*crLHS12;
const double crLHS251 = DN(1,2)*crLHS250;
const double crLHS252 = DN(2,0)*crLHS250;
const double crLHS253 = DN(2,1)*crLHS250;
const double crLHS254 = crLHS234 + crLHS236;
const double crLHS255 = DN(2,2)*crLHS250;
const double crLHS256 = DN(1,1)*N[2];
const double crLHS257 = DN(2,1)*N[1];
const double crLHS258 = DN(3,0)*crLHS250;
const double crLHS259 = DN(3,1)*crLHS250;
const double crLHS260 = crLHS241 + crLHS243;
const double crLHS261 = DN(3,2)*crLHS250;
const double crLHS262 = DN(1,1)*N[3];
const double crLHS263 = DN(3,1)*N[1];
const double crLHS264 = DN(1,2)*DN(1,2);
const double crLHS265 = DN(1,2)*crLHS12;
const double crLHS266 = DN(2,0)*crLHS265;
const double crLHS267 = DN(2,1)*crLHS265;
const double crLHS268 = DN(2,2)*crLHS265;
const double crLHS269 = DN(1,2)*N[2];
const double crLHS270 = DN(2,2)*N[1];
const double crLHS271 = DN(3,0)*crLHS265;
const double crLHS272 = DN(3,1)*crLHS265;
const double crLHS273 = DN(3,2)*crLHS265;
const double crLHS274 = DN(1,2)*N[3];
const double crLHS275 = DN(3,2)*N[1];
const double crLHS276 = gauss_weight*(N[1] + crLHS20*(1.0*crLHS55 + 1.0*crLHS56 + 1.0*crLHS57));
const double crLHS277 = DN(1,0)*crLHS25;
const double crLHS278 = DN(1,1)*crLHS25;
const double crLHS279 = DN(1,2)*crLHS25;
const double crLHS280 = gauss_weight*(DN(2,0)*crLHS277 + DN(2,1)*crLHS278 + DN(2,2)*crLHS279 + N[2]*crLHS232);
const double crLHS281 = gauss_weight*(DN(3,0)*crLHS277 + DN(3,1)*crLHS278 + DN(3,2)*crLHS279 + N[3]*crLHS232);
const double crLHS282 = N[2]*crLHS27;
const double crLHS283 = N[2]*crLHS14 + crLHS26*crLHS282 + crLHS26*crLHS84 - crLHS26*crLHS85;
const double crLHS284 = N[2]*crLHS55 + crLHS282*crLHS59 + crLHS59*crLHS84 - crLHS59*crLHS85;
const double crLHS285 = DN(2,0)*DN(2,0);
const double crLHS286 = N[2]*N[2];
const double crLHS287 = N[2]*crLHS84 + crLHS15*crLHS286 + crLHS282*crLHS88 + crLHS286*sigma + crLHS84*crLHS88 - crLHS85*crLHS88;
const double crLHS288 = DN(2,0)*crLHS12;
const double crLHS289 = DN(2,1)*crLHS288;
const double crLHS290 = DN(2,2)*crLHS288;
const double crLHS291 = crLHS25*crLHS85;
const double crLHS292 = N[2]*crLHS44;
const double crLHS293 = crLHS25*crLHS84;
const double crLHS294 = gauss_weight*(N[2]*crLHS46 - N[2] + crLHS12*crLHS292 - crLHS291 + crLHS293);
const double crLHS295 = N[3]*crLHS85 + N[3]*crLHS86;
const double crLHS296 = DN(3,0)*crLHS288 + crLHS295;
const double crLHS297 = N[2]*crLHS112 - crLHS115*crLHS291 + crLHS115*crLHS293 + crLHS116*crLHS282;
const double crLHS298 = DN(3,1)*crLHS288;
const double crLHS299 = DN(3,2)*crLHS288;
const double crLHS300 = DN(2,0)*N[3];
const double crLHS301 = DN(3,0)*N[2];
const double crLHS302 = crLHS283 + crLHS82;
const double crLHS303 = crLHS234 + crLHS284;
const double crLHS304 = DN(2,1)*DN(2,1);
const double crLHS305 = DN(2,1)*crLHS12;
const double crLHS306 = DN(2,2)*crLHS305;
const double crLHS307 = DN(3,0)*crLHS305;
const double crLHS308 = DN(3,1)*crLHS305;
const double crLHS309 = crLHS295 + crLHS297;
const double crLHS310 = DN(3,2)*crLHS305;
const double crLHS311 = DN(2,1)*N[3];
const double crLHS312 = DN(3,1)*N[2];
const double crLHS313 = DN(2,2)*DN(2,2);
const double crLHS314 = DN(2,2)*crLHS12;
const double crLHS315 = DN(3,0)*crLHS314;
const double crLHS316 = DN(3,1)*crLHS314;
const double crLHS317 = DN(3,2)*crLHS314;
const double crLHS318 = DN(2,2)*N[3];
const double crLHS319 = DN(3,2)*N[2];
const double crLHS320 = gauss_weight*(N[2] + crLHS20*(1.0*crLHS84 + 1.0*crLHS85 + 1.0*crLHS86));
const double crLHS321 = DN(2,0)*crLHS25;
const double crLHS322 = DN(2,1)*crLHS25;
const double crLHS323 = DN(2,2)*crLHS25;
const double crLHS324 = gauss_weight*(DN(3,0)*crLHS321 + DN(3,1)*crLHS322 + DN(3,2)*crLHS323 + N[3]*crLHS292);
const double crLHS325 = N[3]*crLHS27;
const double crLHS326 = N[3]*crLHS14 + crLHS112*crLHS26 - crLHS113*crLHS26 + crLHS26*crLHS325;
const double crLHS327 = N[3]*crLHS55 + crLHS112*crLHS59 - crLHS113*crLHS59 + crLHS325*crLHS59;
const double crLHS328 = N[3]*crLHS84 + crLHS112*crLHS88 - crLHS113*crLHS88 + crLHS325*crLHS88;
const double crLHS329 = DN(3,0)*DN(3,0);
const double crLHS330 = N[3]*N[3];
const double crLHS331 = N[3]*crLHS112 + crLHS112*crLHS116 - crLHS113*crLHS116 + crLHS116*crLHS325 + crLHS15*crLHS330 + crLHS330*sigma;
const double crLHS332 = DN(3,0)*crLHS12;
const double crLHS333 = DN(3,1)*crLHS332;
const double crLHS334 = DN(3,2)*crLHS332;
const double crLHS335 = gauss_weight*(N[3]*crLHS46 + N[3]*crLHS75 - N[3] + crLHS112*crLHS25 - crLHS113*crLHS25);
const double crLHS336 = crLHS110 + crLHS326;
const double crLHS337 = crLHS241 + crLHS327;
const double crLHS338 = crLHS295 + crLHS328;
const double crLHS339 = DN(3,1)*DN(3,1);
const double crLHS340 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS341 = DN(3,2)*DN(3,2);
const double crLHS342 = gauss_weight*(N[3] + crLHS20*(1.0*crLHS112 + 1.0*crLHS113 + 1.0*crLHS114));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS29);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS30 + DN(0,1)*crLHS32 + DN(0,2)*crLHS35 + crLHS37);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS40 + DN(0,2)*crLHS42 + crLHS43);
rLHS(0,3)+=DN(0,0)*crLHS47;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS48 + DN(0,1)*crLHS50 + DN(0,2)*crLHS52 + crLHS54 + crLHS60);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + DN(0,2)*crLHS66 + crLHS67);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS68 + DN(0,1)*crLHS70 + DN(0,2)*crLHS72 + crLHS73);
rLHS(0,7)+=gauss_weight*(DN(1,0)*crLHS22 - DN(1,0)*crLHS24 + crLHS46*crLHS76 + crLHS74*crLHS75 - crLHS74);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS79 + DN(0,2)*crLHS81 + crLHS83 + crLHS89);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS90 + DN(0,1)*crLHS92 + DN(0,2)*crLHS95 + crLHS96);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS97 + DN(0,1)*crLHS99 + DN(0,2)*crLHS101 + crLHS102);
rLHS(0,11)+=gauss_weight*(DN(2,0)*crLHS22 - DN(2,0)*crLHS24 + crLHS103*crLHS75 - crLHS103 + crLHS104*crLHS46);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS105 + DN(0,1)*crLHS107 + DN(0,2)*crLHS109 + crLHS111 + crLHS117);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS118 + DN(0,1)*crLHS120 + DN(0,2)*crLHS123 + crLHS124);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS125 + DN(0,1)*crLHS127 + DN(0,2)*crLHS129 + crLHS130);
rLHS(0,15)+=gauss_weight*(DN(3,0)*crLHS22 - DN(3,0)*crLHS24 + crLHS131*crLHS75 - crLHS131 + crLHS132*crLHS46);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS133 + DN(0,2)*crLHS134 + crLHS37);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS32 + DN(0,1)*crLHS135 + DN(0,2)*crLHS137 + crLHS12*crLHS138 + crLHS29);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS139 + DN(0,2)*crLHS141 + crLHS143);
rLHS(1,3)+=DN(0,1)*crLHS47;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS50 + DN(0,1)*crLHS144 + DN(0,2)*crLHS145 + crLHS146);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS147 + DN(0,2)*crLHS149 + crLHS150 + crLHS151);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS70 + DN(0,1)*crLHS152 + DN(0,2)*crLHS154 + crLHS155);
rLHS(1,7)+=gauss_weight*(DN(1,1)*crLHS22 - DN(1,1)*crLHS24 + crLHS156*crLHS75 - crLHS156 + crLHS157*crLHS46);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS158 + DN(0,2)*crLHS159 + crLHS160);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS92 + DN(0,1)*crLHS161 + DN(0,2)*crLHS163 + crLHS164 + crLHS165);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS99 + DN(0,1)*crLHS166 + DN(0,2)*crLHS168 + crLHS169);
rLHS(1,11)+=gauss_weight*(DN(2,1)*crLHS22 - DN(2,1)*crLHS24 + crLHS170*crLHS75 - crLHS170 + crLHS171*crLHS46);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS107 + DN(0,1)*crLHS172 + DN(0,2)*crLHS173 + crLHS174);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS120 + DN(0,1)*crLHS175 + DN(0,2)*crLHS177 + crLHS178 + crLHS179);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS180 + DN(0,2)*crLHS182 + crLHS183);
rLHS(1,15)+=gauss_weight*(DN(3,1)*crLHS22 - DN(3,1)*crLHS24 + crLHS184*crLHS75 - crLHS184 + crLHS185*crLHS46);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS134 + DN(0,2)*crLHS186 + crLHS43);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS35 + DN(0,1)*crLHS137 + DN(0,2)*crLHS187 + crLHS143);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS42 + DN(0,1)*crLHS141 + DN(0,2)*crLHS188 + crLHS12*crLHS189 + crLHS29);
rLHS(2,3)+=DN(0,2)*crLHS47;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS52 + DN(0,1)*crLHS145 + DN(0,2)*crLHS190 + crLHS192);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS66 + DN(0,1)*crLHS149 + DN(0,2)*crLHS193 + crLHS194);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS72 + DN(0,1)*crLHS154 + DN(0,2)*crLHS195 + crLHS151 + crLHS196);
rLHS(2,7)+=gauss_weight*(DN(1,2)*crLHS22 - DN(1,2)*crLHS24 + crLHS197*crLHS75 - crLHS197 + crLHS198*crLHS46);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS81 + DN(0,1)*crLHS159 + DN(0,2)*crLHS199 + crLHS200);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS95 + DN(0,1)*crLHS163 + DN(0,2)*crLHS201 + crLHS202);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS101 + DN(0,1)*crLHS168 + DN(0,2)*crLHS203 + crLHS165 + crLHS204);
rLHS(2,11)+=gauss_weight*(DN(2,2)*crLHS22 - DN(2,2)*crLHS24 + crLHS205*crLHS75 - crLHS205 + crLHS206*crLHS46);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS109 + DN(0,1)*crLHS173 + DN(0,2)*crLHS207 + crLHS208);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS123 + DN(0,1)*crLHS177 + DN(0,2)*crLHS209 + crLHS210);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS129 + DN(0,1)*crLHS182 + DN(0,2)*crLHS211 + crLHS179 + crLHS212);
rLHS(2,15)+=gauss_weight*(DN(3,2)*crLHS22 - DN(3,2)*crLHS24 + crLHS213*crLHS75 - crLHS213 + crLHS214*crLHS46);
rLHS(3,0)+=DN(0,0)*crLHS215;
rLHS(3,1)+=DN(0,1)*crLHS215;
rLHS(3,2)+=DN(0,2)*crLHS215;
rLHS(3,3)+=gauss_weight*(crLHS13*crLHS44 + crLHS138*crLHS25 + crLHS189*crLHS25 + crLHS25*crLHS5);
rLHS(3,4)+=gauss_weight*(DN(0,0)*crLHS59 + crLHS76);
rLHS(3,5)+=gauss_weight*(DN(0,1)*crLHS59 + crLHS157);
rLHS(3,6)+=gauss_weight*(DN(0,2)*crLHS59 + crLHS198);
rLHS(3,7)+=crLHS219;
rLHS(3,8)+=gauss_weight*(crLHS104 + crLHS216*crLHS87);
rLHS(3,9)+=gauss_weight*(crLHS171 + crLHS217*crLHS87);
rLHS(3,10)+=gauss_weight*(crLHS206 + crLHS218*crLHS87);
rLHS(3,11)+=crLHS220;
rLHS(3,12)+=gauss_weight*(crLHS115*crLHS216 + crLHS132);
rLHS(3,13)+=gauss_weight*(crLHS115*crLHS217 + crLHS185);
rLHS(3,14)+=gauss_weight*(crLHS115*crLHS218 + crLHS214);
rLHS(3,15)+=crLHS221;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS225 + crLHS54);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS30 + DN(1,1)*crLHS32 + DN(1,2)*crLHS35 + crLHS146);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS40 + DN(1,2)*crLHS42 + crLHS192);
rLHS(4,3)+=gauss_weight*(crLHS216*crLHS55 - crLHS216*crLHS56 + crLHS46*crLHS74 + crLHS75*crLHS76 - crLHS76);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS48 + DN(1,1)*crLHS50 + DN(1,2)*crLHS52 + crLHS12*crLHS226 + crLHS228);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + DN(1,2)*crLHS66 + crLHS230);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS68 + DN(1,1)*crLHS70 + DN(1,2)*crLHS72 + crLHS231);
rLHS(4,7)+=DN(1,0)*crLHS233;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS79 + DN(1,2)*crLHS81 + crLHS235 + crLHS236);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS90 + DN(1,1)*crLHS92 + DN(1,2)*crLHS95 + crLHS237);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS97 + DN(1,1)*crLHS99 + DN(1,2)*crLHS101 + crLHS238);
rLHS(4,11)+=gauss_weight*(DN(2,0)*crLHS222 - DN(2,0)*crLHS223 + crLHS239*crLHS75 - crLHS239 + crLHS240*crLHS46);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS105 + DN(1,1)*crLHS107 + DN(1,2)*crLHS109 + crLHS242 + crLHS243);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS118 + DN(1,1)*crLHS120 + DN(1,2)*crLHS123 + crLHS244);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS125 + DN(1,1)*crLHS127 + DN(1,2)*crLHS129 + crLHS245);
rLHS(4,15)+=gauss_weight*(DN(3,0)*crLHS222 - DN(3,0)*crLHS223 + crLHS246*crLHS75 - crLHS246 + crLHS247*crLHS46);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS133 + DN(1,2)*crLHS134 + crLHS67);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS32 + DN(1,1)*crLHS135 + DN(1,2)*crLHS137 + crLHS150 + crLHS248);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS139 + DN(1,2)*crLHS141 + crLHS194);
rLHS(5,3)+=gauss_weight*(crLHS156*crLHS46 + crLHS157*crLHS75 - crLHS157 + crLHS217*crLHS55 - crLHS217*crLHS56);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS50 + DN(1,1)*crLHS144 + DN(1,2)*crLHS145 + crLHS230);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS147 + DN(1,2)*crLHS149 + crLHS12*crLHS249 + crLHS228);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS70 + DN(1,1)*crLHS152 + DN(1,2)*crLHS154 + crLHS251);
rLHS(5,7)+=DN(1,1)*crLHS233;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS158 + DN(1,2)*crLHS159 + crLHS252);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS92 + DN(1,1)*crLHS161 + DN(1,2)*crLHS163 + crLHS253 + crLHS254);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS99 + DN(1,1)*crLHS166 + DN(1,2)*crLHS168 + crLHS255);
rLHS(5,11)+=gauss_weight*(DN(2,1)*crLHS222 - DN(2,1)*crLHS223 + crLHS256*crLHS75 - crLHS256 + crLHS257*crLHS46);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS107 + DN(1,1)*crLHS172 + DN(1,2)*crLHS173 + crLHS258);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS120 + DN(1,1)*crLHS175 + DN(1,2)*crLHS177 + crLHS259 + crLHS260);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS180 + DN(1,2)*crLHS182 + crLHS261);
rLHS(5,15)+=gauss_weight*(DN(3,1)*crLHS222 - DN(3,1)*crLHS223 + crLHS262*crLHS75 - crLHS262 + crLHS263*crLHS46);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS134 + DN(1,2)*crLHS186 + crLHS73);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS35 + DN(1,1)*crLHS137 + DN(1,2)*crLHS187 + crLHS155);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS42 + DN(1,1)*crLHS141 + DN(1,2)*crLHS188 + crLHS196 + crLHS248);
rLHS(6,3)+=gauss_weight*(crLHS197*crLHS46 + crLHS198*crLHS75 - crLHS198 + crLHS218*crLHS55 - crLHS218*crLHS56);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS52 + DN(1,1)*crLHS145 + DN(1,2)*crLHS190 + crLHS231);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS66 + DN(1,1)*crLHS149 + DN(1,2)*crLHS193 + crLHS251);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS72 + DN(1,1)*crLHS154 + DN(1,2)*crLHS195 + crLHS12*crLHS264 + crLHS228);
rLHS(6,7)+=DN(1,2)*crLHS233;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS81 + DN(1,1)*crLHS159 + DN(1,2)*crLHS199 + crLHS266);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS95 + DN(1,1)*crLHS163 + DN(1,2)*crLHS201 + crLHS267);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS101 + DN(1,1)*crLHS168 + DN(1,2)*crLHS203 + crLHS254 + crLHS268);
rLHS(6,11)+=gauss_weight*(DN(2,2)*crLHS222 - DN(2,2)*crLHS223 + crLHS269*crLHS75 - crLHS269 + crLHS270*crLHS46);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS109 + DN(1,1)*crLHS173 + DN(1,2)*crLHS207 + crLHS271);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS123 + DN(1,1)*crLHS177 + DN(1,2)*crLHS209 + crLHS272);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS129 + DN(1,1)*crLHS182 + DN(1,2)*crLHS211 + crLHS260 + crLHS273);
rLHS(6,15)+=gauss_weight*(DN(3,2)*crLHS222 - DN(3,2)*crLHS223 + crLHS274*crLHS75 - crLHS274 + crLHS275*crLHS46);
rLHS(7,0)+=gauss_weight*(DN(1,0)*crLHS26 + crLHS74);
rLHS(7,1)+=gauss_weight*(DN(1,1)*crLHS26 + crLHS156);
rLHS(7,2)+=gauss_weight*(DN(1,2)*crLHS26 + crLHS197);
rLHS(7,3)+=crLHS219;
rLHS(7,4)+=DN(1,0)*crLHS276;
rLHS(7,5)+=DN(1,1)*crLHS276;
rLHS(7,6)+=DN(1,2)*crLHS276;
rLHS(7,7)+=gauss_weight*(crLHS226*crLHS25 + crLHS227*crLHS44 + crLHS249*crLHS25 + crLHS25*crLHS264);
rLHS(7,8)+=gauss_weight*(DN(1,0)*crLHS88 + crLHS240);
rLHS(7,9)+=gauss_weight*(DN(1,1)*crLHS88 + crLHS257);
rLHS(7,10)+=gauss_weight*(DN(1,2)*crLHS88 + crLHS270);
rLHS(7,11)+=crLHS280;
rLHS(7,12)+=gauss_weight*(crLHS115*crLHS277 + crLHS247);
rLHS(7,13)+=gauss_weight*(crLHS115*crLHS278 + crLHS263);
rLHS(7,14)+=gauss_weight*(crLHS115*crLHS279 + crLHS275);
rLHS(7,15)+=crLHS281;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS283 + crLHS83);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS30 + DN(2,1)*crLHS32 + DN(2,2)*crLHS35 + crLHS160);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS40 + DN(2,2)*crLHS42 + crLHS200);
rLHS(8,3)+=gauss_weight*(crLHS103*crLHS46 + crLHS104*crLHS75 - crLHS104 + crLHS216*crLHS84 - crLHS216*crLHS85);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS48 + DN(2,1)*crLHS50 + DN(2,2)*crLHS52 + crLHS235 + crLHS284);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + DN(2,2)*crLHS66 + crLHS252);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS68 + DN(2,1)*crLHS70 + DN(2,2)*crLHS72 + crLHS266);
rLHS(8,7)+=gauss_weight*(crLHS239*crLHS46 + crLHS240*crLHS75 - crLHS240 + crLHS277*crLHS84 - crLHS277*crLHS85);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS79 + DN(2,2)*crLHS81 + crLHS12*crLHS285 + crLHS287);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS90 + DN(2,1)*crLHS92 + DN(2,2)*crLHS95 + crLHS289);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS97 + DN(2,1)*crLHS99 + DN(2,2)*crLHS101 + crLHS290);
rLHS(8,11)+=DN(2,0)*crLHS294;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS105 + DN(2,1)*crLHS107 + DN(2,2)*crLHS109 + crLHS296 + crLHS297);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS118 + DN(2,1)*crLHS120 + DN(2,2)*crLHS123 + crLHS298);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS125 + DN(2,1)*crLHS127 + DN(2,2)*crLHS129 + crLHS299);
rLHS(8,15)+=gauss_weight*(-DN(3,0)*crLHS291 + DN(3,0)*crLHS293 + crLHS300*crLHS75 - crLHS300 + crLHS301*crLHS46);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS133 + DN(2,2)*crLHS134 + crLHS96);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS32 + DN(2,1)*crLHS135 + DN(2,2)*crLHS137 + crLHS164 + crLHS302);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS139 + DN(2,2)*crLHS141 + crLHS202);
rLHS(9,3)+=gauss_weight*(crLHS170*crLHS46 + crLHS171*crLHS75 - crLHS171 + crLHS217*crLHS84 - crLHS217*crLHS85);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS50 + DN(2,1)*crLHS144 + DN(2,2)*crLHS145 + crLHS237);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS147 + DN(2,2)*crLHS149 + crLHS253 + crLHS303);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS70 + DN(2,1)*crLHS152 + DN(2,2)*crLHS154 + crLHS267);
rLHS(9,7)+=gauss_weight*(crLHS256*crLHS46 + crLHS257*crLHS75 - crLHS257 + crLHS278*crLHS84 - crLHS278*crLHS85);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS158 + DN(2,2)*crLHS159 + crLHS289);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS92 + DN(2,1)*crLHS161 + DN(2,2)*crLHS163 + crLHS12*crLHS304 + crLHS287);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS99 + DN(2,1)*crLHS166 + DN(2,2)*crLHS168 + crLHS306);
rLHS(9,11)+=DN(2,1)*crLHS294;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS107 + DN(2,1)*crLHS172 + DN(2,2)*crLHS173 + crLHS307);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS120 + DN(2,1)*crLHS175 + DN(2,2)*crLHS177 + crLHS308 + crLHS309);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS180 + DN(2,2)*crLHS182 + crLHS310);
rLHS(9,15)+=gauss_weight*(-DN(3,1)*crLHS291 + DN(3,1)*crLHS293 + crLHS311*crLHS75 - crLHS311 + crLHS312*crLHS46);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS134 + DN(2,2)*crLHS186 + crLHS102);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS35 + DN(2,1)*crLHS137 + DN(2,2)*crLHS187 + crLHS169);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS42 + DN(2,1)*crLHS141 + DN(2,2)*crLHS188 + crLHS204 + crLHS302);
rLHS(10,3)+=gauss_weight*(crLHS205*crLHS46 + crLHS206*crLHS75 - crLHS206 + crLHS218*crLHS84 - crLHS218*crLHS85);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS52 + DN(2,1)*crLHS145 + DN(2,2)*crLHS190 + crLHS238);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS66 + DN(2,1)*crLHS149 + DN(2,2)*crLHS193 + crLHS255);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS72 + DN(2,1)*crLHS154 + DN(2,2)*crLHS195 + crLHS268 + crLHS303);
rLHS(10,7)+=gauss_weight*(crLHS269*crLHS46 + crLHS270*crLHS75 - crLHS270 + crLHS279*crLHS84 - crLHS279*crLHS85);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS81 + DN(2,1)*crLHS159 + DN(2,2)*crLHS199 + crLHS290);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS95 + DN(2,1)*crLHS163 + DN(2,2)*crLHS201 + crLHS306);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS101 + DN(2,1)*crLHS168 + DN(2,2)*crLHS203 + crLHS12*crLHS313 + crLHS287);
rLHS(10,11)+=DN(2,2)*crLHS294;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS109 + DN(2,1)*crLHS173 + DN(2,2)*crLHS207 + crLHS315);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS123 + DN(2,1)*crLHS177 + DN(2,2)*crLHS209 + crLHS316);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS129 + DN(2,1)*crLHS182 + DN(2,2)*crLHS211 + crLHS309 + crLHS317);
rLHS(10,15)+=gauss_weight*(-DN(3,2)*crLHS291 + DN(3,2)*crLHS293 + crLHS318*crLHS75 - crLHS318 + crLHS319*crLHS46);
rLHS(11,0)+=gauss_weight*(DN(2,0)*crLHS26 + crLHS103);
rLHS(11,1)+=gauss_weight*(DN(2,1)*crLHS26 + crLHS170);
rLHS(11,2)+=gauss_weight*(DN(2,2)*crLHS26 + crLHS205);
rLHS(11,3)+=crLHS220;
rLHS(11,4)+=gauss_weight*(DN(2,0)*crLHS59 + crLHS239);
rLHS(11,5)+=gauss_weight*(DN(2,1)*crLHS59 + crLHS256);
rLHS(11,6)+=gauss_weight*(DN(2,2)*crLHS59 + crLHS269);
rLHS(11,7)+=crLHS280;
rLHS(11,8)+=DN(2,0)*crLHS320;
rLHS(11,9)+=DN(2,1)*crLHS320;
rLHS(11,10)+=DN(2,2)*crLHS320;
rLHS(11,11)+=gauss_weight*(crLHS25*crLHS285 + crLHS25*crLHS304 + crLHS25*crLHS313 + crLHS286*crLHS44);
rLHS(11,12)+=gauss_weight*(DN(2,0)*crLHS116 + crLHS301);
rLHS(11,13)+=gauss_weight*(DN(2,1)*crLHS116 + crLHS312);
rLHS(11,14)+=gauss_weight*(DN(2,2)*crLHS116 + crLHS319);
rLHS(11,15)+=crLHS324;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS111 + crLHS326);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS30 + DN(3,1)*crLHS32 + DN(3,2)*crLHS35 + crLHS174);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS38 + DN(3,1)*crLHS40 + DN(3,2)*crLHS42 + crLHS208);
rLHS(12,3)+=gauss_weight*(crLHS112*crLHS216 - crLHS113*crLHS216 + crLHS131*crLHS46 + crLHS132*crLHS75 - crLHS132);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS48 + DN(3,1)*crLHS50 + DN(3,2)*crLHS52 + crLHS242 + crLHS327);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS61 + DN(3,1)*crLHS63 + DN(3,2)*crLHS66 + crLHS258);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS68 + DN(3,1)*crLHS70 + DN(3,2)*crLHS72 + crLHS271);
rLHS(12,7)+=gauss_weight*(crLHS112*crLHS277 - crLHS113*crLHS277 + crLHS246*crLHS46 + crLHS247*crLHS75 - crLHS247);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS79 + DN(3,2)*crLHS81 + crLHS296 + crLHS328);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS90 + DN(3,1)*crLHS92 + DN(3,2)*crLHS95 + crLHS307);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS97 + DN(3,1)*crLHS99 + DN(3,2)*crLHS101 + crLHS315);
rLHS(12,11)+=gauss_weight*(crLHS112*crLHS321 - crLHS113*crLHS321 + crLHS300*crLHS46 + crLHS301*crLHS75 - crLHS301);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS105 + DN(3,1)*crLHS107 + DN(3,2)*crLHS109 + crLHS12*crLHS329 + crLHS331);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS118 + DN(3,1)*crLHS120 + DN(3,2)*crLHS123 + crLHS333);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS125 + DN(3,1)*crLHS127 + DN(3,2)*crLHS129 + crLHS334);
rLHS(12,15)+=DN(3,0)*crLHS335;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS133 + DN(3,2)*crLHS134 + crLHS124);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS32 + DN(3,1)*crLHS135 + DN(3,2)*crLHS137 + crLHS178 + crLHS336);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS40 + DN(3,1)*crLHS139 + DN(3,2)*crLHS141 + crLHS210);
rLHS(13,3)+=gauss_weight*(crLHS112*crLHS217 - crLHS113*crLHS217 + crLHS184*crLHS46 + crLHS185*crLHS75 - crLHS185);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS50 + DN(3,1)*crLHS144 + DN(3,2)*crLHS145 + crLHS244);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS63 + DN(3,1)*crLHS147 + DN(3,2)*crLHS149 + crLHS259 + crLHS337);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS70 + DN(3,1)*crLHS152 + DN(3,2)*crLHS154 + crLHS272);
rLHS(13,7)+=gauss_weight*(crLHS112*crLHS278 - crLHS113*crLHS278 + crLHS262*crLHS46 + crLHS263*crLHS75 - crLHS263);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS79 + DN(3,1)*crLHS158 + DN(3,2)*crLHS159 + crLHS298);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS92 + DN(3,1)*crLHS161 + DN(3,2)*crLHS163 + crLHS308 + crLHS338);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS99 + DN(3,1)*crLHS166 + DN(3,2)*crLHS168 + crLHS316);
rLHS(13,11)+=gauss_weight*(crLHS112*crLHS322 - crLHS113*crLHS322 + crLHS311*crLHS46 + crLHS312*crLHS75 - crLHS312);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS107 + DN(3,1)*crLHS172 + DN(3,2)*crLHS173 + crLHS333);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS120 + DN(3,1)*crLHS175 + DN(3,2)*crLHS177 + crLHS12*crLHS339 + crLHS331);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS180 + DN(3,2)*crLHS182 + crLHS340);
rLHS(13,15)+=DN(3,1)*crLHS335;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS134 + DN(3,2)*crLHS186 + crLHS130);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS35 + DN(3,1)*crLHS137 + DN(3,2)*crLHS187 + crLHS183);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS42 + DN(3,1)*crLHS141 + DN(3,2)*crLHS188 + crLHS212 + crLHS336);
rLHS(14,3)+=gauss_weight*(crLHS112*crLHS218 - crLHS113*crLHS218 + crLHS213*crLHS46 + crLHS214*crLHS75 - crLHS214);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS52 + DN(3,1)*crLHS145 + DN(3,2)*crLHS190 + crLHS245);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS66 + DN(3,1)*crLHS149 + DN(3,2)*crLHS193 + crLHS261);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS72 + DN(3,1)*crLHS154 + DN(3,2)*crLHS195 + crLHS273 + crLHS337);
rLHS(14,7)+=gauss_weight*(crLHS112*crLHS279 - crLHS113*crLHS279 + crLHS274*crLHS46 + crLHS275*crLHS75 - crLHS275);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS81 + DN(3,1)*crLHS159 + DN(3,2)*crLHS199 + crLHS299);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS95 + DN(3,1)*crLHS163 + DN(3,2)*crLHS201 + crLHS310);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS101 + DN(3,1)*crLHS168 + DN(3,2)*crLHS203 + crLHS317 + crLHS338);
rLHS(14,11)+=gauss_weight*(crLHS112*crLHS323 - crLHS113*crLHS323 + crLHS318*crLHS46 + crLHS319*crLHS75 - crLHS319);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS109 + DN(3,1)*crLHS173 + DN(3,2)*crLHS207 + crLHS334);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS123 + DN(3,1)*crLHS177 + DN(3,2)*crLHS209 + crLHS340);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS129 + DN(3,1)*crLHS182 + DN(3,2)*crLHS211 + crLHS12*crLHS341 + crLHS331);
rLHS(14,15)+=DN(3,2)*crLHS335;
rLHS(15,0)+=gauss_weight*(DN(3,0)*crLHS26 + crLHS131);
rLHS(15,1)+=gauss_weight*(DN(3,1)*crLHS26 + crLHS184);
rLHS(15,2)+=gauss_weight*(DN(3,2)*crLHS26 + crLHS213);
rLHS(15,3)+=crLHS221;
rLHS(15,4)+=gauss_weight*(DN(3,0)*crLHS59 + crLHS246);
rLHS(15,5)+=gauss_weight*(DN(3,1)*crLHS59 + crLHS262);
rLHS(15,6)+=gauss_weight*(DN(3,2)*crLHS59 + crLHS274);
rLHS(15,7)+=crLHS281;
rLHS(15,8)+=gauss_weight*(DN(3,0)*crLHS88 + crLHS300);
rLHS(15,9)+=gauss_weight*(DN(3,1)*crLHS88 + crLHS311);
rLHS(15,10)+=gauss_weight*(DN(3,2)*crLHS88 + crLHS318);
rLHS(15,11)+=crLHS324;
rLHS(15,12)+=DN(3,0)*crLHS342;
rLHS(15,13)+=DN(3,1)*crLHS342;
rLHS(15,14)+=DN(3,2)*crLHS342;
rLHS(15,15)+=gauss_weight*(crLHS25*crLHS329 + crLHS25*crLHS339 + crLHS25*crLHS341 + crLHS330*crLHS44);

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
    const BoundedMatrix<double,2,3>& r_v_sol_frac = rData.SolidFractionVelocity;
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
const double crRHS1 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crRHS2 = crRHS1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crRHS3 = N[0]*r_v_sol_frac(0,0);
const double crRHS4 = N[1]*r_v_sol_frac(1,0);
const double crRHS5 = N[2]*r_v_sol_frac(2,0);
const double crRHS6 = N[0]*v(0,0);
const double crRHS7 = N[1]*v(1,0);
const double crRHS8 = N[2]*v(2,0);
const double crRHS9 = crRHS3 + crRHS4 + crRHS5 - crRHS6 - crRHS7 - crRHS8;
const double crRHS10 = N[0]*sigma;
const double crRHS11 = crRHS1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crRHS12 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS13 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS14 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS15 = crRHS1*(crRHS12*crRHS13 + crRHS14*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crRHS16 = sigma*stab_c3;
const double crRHS17 = crRHS1*stab_c2*sqrt(crRHS13*crRHS13 + crRHS14*crRHS14);
const double crRHS18 = 1.0/crRHS1;
const double crRHS19 = crRHS18*(crRHS13*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crRHS14*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crRHS20 = crRHS18*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crRHS21 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS22 = crRHS12 + crRHS21;
const double crRHS23 = (crRHS19 + crRHS20 + crRHS22)*(mu + (crRHS16 + crRHS17*h)*1.0/stab_c1);
const double crRHS24 = 1.0/h;
const double crRHS25 = 1.0*1.0/(crRHS1*dyn_tau*1.0/dt + crRHS16*crRHS24 + crRHS17*crRHS24 + mu*stab_c1*1.0/(h*h));
const double crRHS26 = crRHS25*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crRHS11 + crRHS15 - crRHS2 + sigma*(-crRHS3 - crRHS4 - crRHS5 + crRHS6 + crRHS7 + crRHS8));
const double crRHS27 = crRHS1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crRHS28 = N[0]*crRHS27;
const double crRHS29 = crRHS1*(DN(0,0)*crRHS13 + DN(0,1)*crRHS14);
const double crRHS30 = crRHS1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS31 = N[0]*r_v_sol_frac(0,1);
const double crRHS32 = N[1]*r_v_sol_frac(1,1);
const double crRHS33 = N[2]*r_v_sol_frac(2,1);
const double crRHS34 = N[0]*v(0,1);
const double crRHS35 = N[1]*v(1,1);
const double crRHS36 = N[2]*v(2,1);
const double crRHS37 = crRHS31 + crRHS32 + crRHS33 - crRHS34 - crRHS35 - crRHS36;
const double crRHS38 = crRHS1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crRHS39 = crRHS1*(crRHS13*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crRHS14*crRHS21);
const double crRHS40 = crRHS25*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS30 + crRHS38 + crRHS39 + sigma*(-crRHS31 - crRHS32 - crRHS33 + crRHS34 + crRHS35 + crRHS36));
const double crRHS41 = N[1]*sigma;
const double crRHS42 = N[1]*crRHS27;
const double crRHS43 = crRHS1*(DN(1,0)*crRHS13 + DN(1,1)*crRHS14);
const double crRHS44 = N[2]*sigma;
const double crRHS45 = N[2]*crRHS27;
const double crRHS46 = crRHS1*(DN(2,0)*crRHS13 + DN(2,1)*crRHS14);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS23 + DN(0,0)*stress[0] + DN(0,1)*stress[2] + N[0]*crRHS11 + N[0]*crRHS15 - N[0]*crRHS2 - crRHS10*crRHS26 - crRHS10*crRHS9 + crRHS26*crRHS28 + crRHS26*crRHS29);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS23 + DN(0,1)*stress[1] - N[0]*crRHS30 + N[0]*crRHS38 + N[0]*crRHS39 - crRHS10*crRHS37 - crRHS10*crRHS40 + crRHS28*crRHS40 + crRHS29*crRHS40);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS26 + DN(0,1)*crRHS40 + N[0]*crRHS19 + N[0]*crRHS20 + N[0]*crRHS22);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS23 + DN(1,0)*stress[0] + DN(1,1)*stress[2] + N[1]*crRHS11 + N[1]*crRHS15 - N[1]*crRHS2 - crRHS26*crRHS41 + crRHS26*crRHS42 + crRHS26*crRHS43 - crRHS41*crRHS9);
rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS23 + DN(1,1)*stress[1] - N[1]*crRHS30 + N[1]*crRHS38 + N[1]*crRHS39 - crRHS37*crRHS41 - crRHS40*crRHS41 + crRHS40*crRHS42 + crRHS40*crRHS43);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS26 + DN(1,1)*crRHS40 + N[1]*crRHS19 + N[1]*crRHS20 + N[1]*crRHS22);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS23 + DN(2,0)*stress[0] + DN(2,1)*stress[2] + N[2]*crRHS11 + N[2]*crRHS15 - N[2]*crRHS2 - crRHS26*crRHS44 + crRHS26*crRHS45 + crRHS26*crRHS46 - crRHS44*crRHS9);
rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS23 + DN(2,1)*stress[1] - N[2]*crRHS30 + N[2]*crRHS38 + N[2]*crRHS39 - crRHS37*crRHS44 - crRHS40*crRHS44 + crRHS40*crRHS45 + crRHS40*crRHS46);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS26 + DN(2,1)*crRHS40 + N[2]*crRHS19 + N[2]*crRHS20 + N[2]*crRHS22);

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
    const BoundedMatrix<double,3,4>& r_v_sol_frac = rData.SolidFractionVelocity;
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
const double crRHS1 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crRHS2 = crRHS1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS3 = N[0]*r_v_sol_frac(0,0);
const double crRHS4 = N[1]*r_v_sol_frac(1,0);
const double crRHS5 = N[2]*r_v_sol_frac(2,0);
const double crRHS6 = N[3]*r_v_sol_frac(3,0);
const double crRHS7 = N[0]*v(0,0);
const double crRHS8 = N[1]*v(1,0);
const double crRHS9 = N[2]*v(2,0);
const double crRHS10 = N[3]*v(3,0);
const double crRHS11 = -crRHS10 + crRHS3 + crRHS4 + crRHS5 + crRHS6 - crRHS7 - crRHS8 - crRHS9;
const double crRHS12 = N[0]*sigma;
const double crRHS13 = crRHS1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crRHS14 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS15 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS16 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS17 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS18 = crRHS1*(crRHS14*crRHS15 + crRHS16*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crRHS17*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crRHS19 = sigma*stab_c3;
const double crRHS20 = crRHS1*stab_c2*sqrt(crRHS15*crRHS15 + crRHS16*crRHS16 + crRHS17*crRHS17);
const double crRHS21 = 1.0/crRHS1;
const double crRHS22 = crRHS21*(crRHS15*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crRHS16*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]) + crRHS17*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]));
const double crRHS23 = crRHS21*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crRHS24 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS25 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS26 = crRHS14 + crRHS24 + crRHS25;
const double crRHS27 = (crRHS22 + crRHS23 + crRHS26)*(mu + (crRHS19 + crRHS20*h)*1.0/stab_c1);
const double crRHS28 = 1.0/h;
const double crRHS29 = 1.0*1.0/(crRHS1*dyn_tau*1.0/dt + crRHS19*crRHS28 + crRHS20*crRHS28 + mu*stab_c1*1.0/(h*h));
const double crRHS30 = crRHS29*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crRHS13 + crRHS18 - crRHS2 + sigma*(crRHS10 - crRHS3 - crRHS4 - crRHS5 - crRHS6 + crRHS7 + crRHS8 + crRHS9));
const double crRHS31 = crRHS1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crRHS32 = N[0]*crRHS31;
const double crRHS33 = crRHS1*(DN(0,0)*crRHS15 + DN(0,1)*crRHS16 + DN(0,2)*crRHS17);
const double crRHS34 = crRHS1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS35 = N[0]*r_v_sol_frac(0,1);
const double crRHS36 = N[1]*r_v_sol_frac(1,1);
const double crRHS37 = N[2]*r_v_sol_frac(2,1);
const double crRHS38 = N[3]*r_v_sol_frac(3,1);
const double crRHS39 = N[0]*v(0,1);
const double crRHS40 = N[1]*v(1,1);
const double crRHS41 = N[2]*v(2,1);
const double crRHS42 = N[3]*v(3,1);
const double crRHS43 = crRHS35 + crRHS36 + crRHS37 + crRHS38 - crRHS39 - crRHS40 - crRHS41 - crRHS42;
const double crRHS44 = crRHS1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crRHS45 = crRHS1*(crRHS15*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crRHS16*crRHS24 + crRHS17*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crRHS46 = crRHS29*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS34 + crRHS44 + crRHS45 + sigma*(-crRHS35 - crRHS36 - crRHS37 - crRHS38 + crRHS39 + crRHS40 + crRHS41 + crRHS42));
const double crRHS47 = crRHS1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS48 = N[0]*r_v_sol_frac(0,2);
const double crRHS49 = N[1]*r_v_sol_frac(1,2);
const double crRHS50 = N[2]*r_v_sol_frac(2,2);
const double crRHS51 = N[3]*r_v_sol_frac(3,2);
const double crRHS52 = N[0]*v(0,2);
const double crRHS53 = N[1]*v(1,2);
const double crRHS54 = N[2]*v(2,2);
const double crRHS55 = N[3]*v(3,2);
const double crRHS56 = crRHS48 + crRHS49 + crRHS50 + crRHS51 - crRHS52 - crRHS53 - crRHS54 - crRHS55;
const double crRHS57 = crRHS1*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crRHS58 = crRHS1*(crRHS15*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crRHS16*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crRHS17*crRHS25);
const double crRHS59 = crRHS29*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS47 + crRHS57 + crRHS58 + sigma*(-crRHS48 - crRHS49 - crRHS50 - crRHS51 + crRHS52 + crRHS53 + crRHS54 + crRHS55));
const double crRHS60 = N[1]*sigma;
const double crRHS61 = N[1]*crRHS31;
const double crRHS62 = crRHS1*(DN(1,0)*crRHS15 + DN(1,1)*crRHS16 + DN(1,2)*crRHS17);
const double crRHS63 = N[2]*sigma;
const double crRHS64 = N[2]*crRHS31;
const double crRHS65 = crRHS1*(DN(2,0)*crRHS15 + DN(2,1)*crRHS16 + DN(2,2)*crRHS17);
const double crRHS66 = N[3]*sigma;
const double crRHS67 = N[3]*crRHS31;
const double crRHS68 = crRHS1*(DN(3,0)*crRHS15 + DN(3,1)*crRHS16 + DN(3,2)*crRHS17);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS27 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] + N[0]*crRHS13 + N[0]*crRHS18 - N[0]*crRHS2 - crRHS11*crRHS12 - crRHS12*crRHS30 + crRHS30*crRHS32 + crRHS30*crRHS33);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS27 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS34 + N[0]*crRHS44 + N[0]*crRHS45 - crRHS12*crRHS43 - crRHS12*crRHS46 + crRHS32*crRHS46 + crRHS33*crRHS46);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS27 + DN(0,2)*stress[2] - N[0]*crRHS47 + N[0]*crRHS57 + N[0]*crRHS58 - crRHS12*crRHS56 - crRHS12*crRHS59 + crRHS32*crRHS59 + crRHS33*crRHS59);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS30 + DN(0,1)*crRHS46 + DN(0,2)*crRHS59 + N[0]*crRHS22 + N[0]*crRHS23 + N[0]*crRHS26);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS27 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] + N[1]*crRHS13 + N[1]*crRHS18 - N[1]*crRHS2 - crRHS11*crRHS60 - crRHS30*crRHS60 + crRHS30*crRHS61 + crRHS30*crRHS62);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS27 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS34 + N[1]*crRHS44 + N[1]*crRHS45 - crRHS43*crRHS60 - crRHS46*crRHS60 + crRHS46*crRHS61 + crRHS46*crRHS62);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS27 + DN(1,2)*stress[2] - N[1]*crRHS47 + N[1]*crRHS57 + N[1]*crRHS58 - crRHS56*crRHS60 - crRHS59*crRHS60 + crRHS59*crRHS61 + crRHS59*crRHS62);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS30 + DN(1,1)*crRHS46 + DN(1,2)*crRHS59 + N[1]*crRHS22 + N[1]*crRHS23 + N[1]*crRHS26);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS27 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] + N[2]*crRHS13 + N[2]*crRHS18 - N[2]*crRHS2 - crRHS11*crRHS63 - crRHS30*crRHS63 + crRHS30*crRHS64 + crRHS30*crRHS65);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS27 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS34 + N[2]*crRHS44 + N[2]*crRHS45 - crRHS43*crRHS63 - crRHS46*crRHS63 + crRHS46*crRHS64 + crRHS46*crRHS65);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS27 + DN(2,2)*stress[2] - N[2]*crRHS47 + N[2]*crRHS57 + N[2]*crRHS58 - crRHS56*crRHS63 - crRHS59*crRHS63 + crRHS59*crRHS64 + crRHS59*crRHS65);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS30 + DN(2,1)*crRHS46 + DN(2,2)*crRHS59 + N[2]*crRHS22 + N[2]*crRHS23 + N[2]*crRHS26);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS27 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] + N[3]*crRHS13 + N[3]*crRHS18 - N[3]*crRHS2 - crRHS11*crRHS66 - crRHS30*crRHS66 + crRHS30*crRHS67 + crRHS30*crRHS68);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS27 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS34 + N[3]*crRHS44 + N[3]*crRHS45 - crRHS43*crRHS66 - crRHS46*crRHS66 + crRHS46*crRHS67 + crRHS46*crRHS68);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS27 + DN(3,2)*stress[2] - N[3]*crRHS47 + N[3]*crRHS57 + N[3]*crRHS58 - crRHS56*crRHS66 - crRHS59*crRHS66 + crRHS59*crRHS67 + crRHS59*crRHS68);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS30 + DN(3,1)*crRHS46 + DN(3,2)*crRHS59 + N[3]*crRHS22 + N[3]*crRHS23 + N[3]*crRHS26);

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
