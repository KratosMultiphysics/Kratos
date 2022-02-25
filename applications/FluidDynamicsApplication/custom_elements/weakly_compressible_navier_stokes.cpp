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

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = pow(DN(0,0), 2);
const double clhs4 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double clhs5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs7 = clhs4*stab_c2*sqrt(pow(clhs5, 2) + pow(clhs6, 2));
const double clhs8 = clhs7*h/stab_c1 + mu;
const double clhs9 = DN(0,0)*clhs5 + DN(0,1)*clhs6;
const double clhs10 = N[0]*clhs4;
const double clhs11 = pow(N[0], 2);
const double clhs12 = bdf0*clhs4;
const double clhs13 = N[0]*bdf0;
const double clhs14 = clhs13 + clhs9;
const double clhs15 = 1.0/(clhs4*dyn_tau/dt + clhs7/h + mu*stab_c1/pow(h, 2));
const double clhs16 = clhs15*pow(clhs4, 2);
const double clhs17 = clhs16*clhs9;
const double clhs18 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs19 = clhs16*clhs18;
const double clhs20 = N[0]*clhs19;
const double clhs21 = clhs10*clhs9 + clhs11*clhs12 + clhs14*clhs17 + clhs14*clhs20;
const double clhs22 = C(0,1)*DN(0,1) + clhs1;
const double clhs23 = C(1,2)*DN(0,1);
const double clhs24 = C(2,2)*DN(0,0) + clhs23;
const double clhs25 = DN(0,0)*clhs8;
const double clhs26 = DN(0,1)*clhs25;
const double clhs27 = 1/(clhs4*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2));
const double clhs28 = clhs13*clhs27;
const double clhs29 = clhs28*clhs8;
const double clhs30 = clhs15*clhs18;
const double clhs31 = clhs15*clhs4;
const double clhs32 = clhs31*clhs9;
const double clhs33 = -N[0] + clhs10*clhs30 + clhs29 + clhs32;
const double clhs34 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs35 = C(0,2)*DN(1,0);
const double clhs36 = C(2,2)*DN(1,1) + clhs35;
const double clhs37 = N[1]*clhs4;
const double clhs38 = clhs13*clhs37;
const double clhs39 = DN(1,0)*clhs25 + clhs38;
const double clhs40 = DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs41 = N[1]*bdf0;
const double clhs42 = clhs40 + clhs41;
const double clhs43 = clhs10*clhs40 + clhs17*clhs42 + clhs20*clhs42;
const double clhs44 = C(0,1)*DN(1,1) + clhs35;
const double clhs45 = C(1,2)*DN(1,1);
const double clhs46 = C(2,2)*DN(1,0) + clhs45;
const double clhs47 = DN(1,1)*clhs25;
const double clhs48 = DN(0,0)*N[1];
const double clhs49 = clhs27*clhs41;
const double clhs50 = DN(1,0)*N[0];
const double clhs51 = clhs18*clhs31;
const double clhs52 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs53 = C(0,2)*DN(2,0);
const double clhs54 = C(2,2)*DN(2,1) + clhs53;
const double clhs55 = N[2]*clhs4;
const double clhs56 = clhs13*clhs55;
const double clhs57 = DN(2,0)*clhs25 + clhs56;
const double clhs58 = DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs59 = N[2]*bdf0;
const double clhs60 = clhs58 + clhs59;
const double clhs61 = clhs10*clhs58 + clhs17*clhs60 + clhs20*clhs60;
const double clhs62 = C(0,1)*DN(2,1) + clhs53;
const double clhs63 = C(1,2)*DN(2,1);
const double clhs64 = C(2,2)*DN(2,0) + clhs63;
const double clhs65 = DN(2,1)*clhs25;
const double clhs66 = DN(0,0)*N[2];
const double clhs67 = clhs27*clhs59;
const double clhs68 = DN(2,0)*N[0];
const double clhs69 = C(0,1)*DN(0,0) + clhs23;
const double clhs70 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs71 = pow(DN(0,1), 2);
const double clhs72 = C(0,1)*DN(1,0) + clhs45;
const double clhs73 = DN(0,1)*clhs8;
const double clhs74 = DN(1,0)*clhs73;
const double clhs75 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs76 = DN(1,1)*clhs73 + clhs38;
const double clhs77 = DN(0,1)*N[1];
const double clhs78 = DN(1,1)*N[0];
const double clhs79 = C(0,1)*DN(2,0) + clhs63;
const double clhs80 = DN(2,0)*clhs73;
const double clhs81 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs82 = DN(2,1)*clhs73 + clhs56;
const double clhs83 = DN(0,1)*N[2];
const double clhs84 = DN(2,1)*N[0];
const double clhs85 = clhs14*clhs31;
const double clhs86 = N[0] + clhs85;
const double clhs87 = bdf0*clhs27;
const double clhs88 = clhs31*clhs42;
const double clhs89 = DN(0,0)*clhs15;
const double clhs90 = DN(0,1)*clhs15;
const double clhs91 = DN(1,0)*clhs89 + DN(1,1)*clhs90 + N[1]*clhs28;
const double clhs92 = clhs31*clhs60;
const double clhs93 = DN(2,0)*clhs89 + DN(2,1)*clhs90 + N[2]*clhs28;
const double clhs94 = clhs16*clhs40;
const double clhs95 = N[1]*clhs19;
const double clhs96 = clhs14*clhs94 + clhs14*clhs95 + clhs37*clhs9;
const double clhs97 = DN(1,0)*clhs8;
const double clhs98 = clhs31*clhs40;
const double clhs99 = pow(DN(1,0), 2);
const double clhs100 = pow(N[1], 2);
const double clhs101 = clhs100*clhs12 + clhs37*clhs40 + clhs42*clhs94 + clhs42*clhs95;
const double clhs102 = DN(1,1)*clhs97;
const double clhs103 = clhs49*clhs8;
const double clhs104 = -N[1] + clhs103 + clhs30*clhs37 + clhs98;
const double clhs105 = clhs41*clhs55;
const double clhs106 = DN(2,0)*clhs97 + clhs105;
const double clhs107 = clhs37*clhs58 + clhs60*clhs94 + clhs60*clhs95;
const double clhs108 = DN(2,1)*clhs97;
const double clhs109 = DN(1,0)*N[2];
const double clhs110 = DN(2,0)*N[1];
const double clhs111 = DN(1,1)*clhs8;
const double clhs112 = pow(DN(1,1), 2);
const double clhs113 = DN(2,0)*clhs111;
const double clhs114 = DN(2,1)*clhs111 + clhs105;
const double clhs115 = DN(1,1)*N[2];
const double clhs116 = DN(2,1)*N[1];
const double clhs117 = N[1] + clhs88;
const double clhs118 = DN(1,0)*DN(2,0)*clhs15 + DN(1,1)*DN(2,1)*clhs15 + N[2]*clhs49;
const double clhs119 = clhs16*clhs58;
const double clhs120 = N[2]*clhs19;
const double clhs121 = clhs119*clhs14 + clhs120*clhs14 + clhs55*clhs9;
const double clhs122 = DN(2,0)*clhs8;
const double clhs123 = clhs31*clhs58;
const double clhs124 = clhs119*clhs42 + clhs120*clhs42 + clhs40*clhs55;
const double clhs125 = pow(DN(2,0), 2);
const double clhs126 = pow(N[2], 2);
const double clhs127 = clhs119*clhs60 + clhs12*clhs126 + clhs120*clhs60 + clhs55*clhs58;
const double clhs128 = DN(2,1)*clhs122;
const double clhs129 = -N[2] + clhs123 + clhs30*clhs55 + clhs67*clhs8;
const double clhs130 = pow(DN(2,1), 2);
const double clhs131 = N[2] + clhs92;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs21 + clhs3*clhs8;
lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + clhs26;
lhs(0,2)=DN(0,0)*clhs33;
lhs(0,3)=DN(0,0)*clhs34 + DN(0,1)*clhs36 + clhs39 + clhs43;
lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + clhs47;
lhs(0,5)=DN(1,0)*clhs32 + clhs25*clhs49 - clhs48 + clhs50*clhs51;
lhs(0,6)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + clhs57 + clhs61;
lhs(0,7)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + clhs65;
lhs(0,8)=DN(2,0)*clhs32 + clhs25*clhs67 + clhs51*clhs68 - clhs66;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs69 + clhs26;
lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs70 + clhs21 + clhs71*clhs8;
lhs(1,2)=DN(0,1)*clhs33;
lhs(1,3)=DN(0,0)*clhs36 + DN(0,1)*clhs72 + clhs74;
lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs75 + clhs43 + clhs76;
lhs(1,5)=DN(1,1)*clhs32 + clhs49*clhs73 + clhs51*clhs78 - clhs77;
lhs(1,6)=DN(0,0)*clhs54 + DN(0,1)*clhs79 + clhs80;
lhs(1,7)=DN(0,0)*clhs64 + DN(0,1)*clhs81 + clhs61 + clhs82;
lhs(1,8)=DN(2,1)*clhs32 + clhs51*clhs84 + clhs67*clhs73 - clhs83;
lhs(2,0)=DN(0,0)*clhs86;
lhs(2,1)=DN(0,1)*clhs86;
lhs(2,2)=clhs11*clhs87 + clhs15*clhs3 + clhs15*clhs71;
lhs(2,3)=DN(0,0)*clhs88 + clhs50;
lhs(2,4)=DN(0,1)*clhs88 + clhs78;
lhs(2,5)=clhs91;
lhs(2,6)=DN(0,0)*clhs92 + clhs68;
lhs(2,7)=DN(0,1)*clhs92 + clhs84;
lhs(2,8)=clhs93;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs39 + clhs96;
lhs(3,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + clhs74;
lhs(3,2)=DN(0,0)*clhs98 + clhs28*clhs97 + clhs48*clhs51 - clhs50;
lhs(3,3)=DN(1,0)*clhs34 + DN(1,1)*clhs36 + clhs101 + clhs8*clhs99;
lhs(3,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + clhs102;
lhs(3,5)=DN(1,0)*clhs104;
lhs(3,6)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + clhs106 + clhs107;
lhs(3,7)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + clhs108;
lhs(3,8)=DN(2,0)*clhs98 - clhs109 + clhs110*clhs51 + clhs67*clhs97;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs69 + clhs47;
lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs70 + clhs76 + clhs96;
lhs(4,2)=DN(0,1)*clhs98 + clhs111*clhs28 + clhs51*clhs77 - clhs78;
lhs(4,3)=DN(1,0)*clhs36 + DN(1,1)*clhs72 + clhs102;
lhs(4,4)=DN(1,0)*clhs46 + DN(1,1)*clhs75 + clhs101 + clhs112*clhs8;
lhs(4,5)=DN(1,1)*clhs104;
lhs(4,6)=DN(1,0)*clhs54 + DN(1,1)*clhs79 + clhs113;
lhs(4,7)=DN(1,0)*clhs64 + DN(1,1)*clhs81 + clhs107 + clhs114;
lhs(4,8)=DN(2,1)*clhs98 + clhs111*clhs67 - clhs115 + clhs116*clhs51;
lhs(5,0)=DN(1,0)*clhs85 + clhs48;
lhs(5,1)=DN(1,1)*clhs85 + clhs77;
lhs(5,2)=clhs91;
lhs(5,3)=DN(1,0)*clhs117;
lhs(5,4)=DN(1,1)*clhs117;
lhs(5,5)=clhs100*clhs87 + clhs112*clhs15 + clhs15*clhs99;
lhs(5,6)=DN(1,0)*clhs92 + clhs110;
lhs(5,7)=DN(1,1)*clhs92 + clhs116;
lhs(5,8)=clhs118;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs121 + clhs57;
lhs(6,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + clhs80;
lhs(6,2)=DN(0,0)*clhs123 + clhs122*clhs28 + clhs51*clhs66 - clhs68;
lhs(6,3)=DN(2,0)*clhs34 + DN(2,1)*clhs36 + clhs106 + clhs124;
lhs(6,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + clhs113;
lhs(6,5)=DN(1,0)*clhs123 + clhs109*clhs51 - clhs110 + clhs122*clhs49;
lhs(6,6)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + clhs125*clhs8 + clhs127;
lhs(6,7)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + clhs128;
lhs(6,8)=DN(2,0)*clhs129;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs69 + clhs65;
lhs(7,1)=DN(2,0)*clhs24 + DN(2,1)*clhs70 + clhs121 + clhs82;
lhs(7,2)=DN(0,1)*clhs123 + DN(2,1)*clhs29 + clhs51*clhs83 - clhs84;
lhs(7,3)=DN(2,0)*clhs36 + DN(2,1)*clhs72 + clhs108;
lhs(7,4)=DN(2,0)*clhs46 + DN(2,1)*clhs75 + clhs114 + clhs124;
lhs(7,5)=DN(1,1)*clhs123 + DN(2,1)*clhs103 + clhs115*clhs51 - clhs116;
lhs(7,6)=DN(2,0)*clhs54 + DN(2,1)*clhs79 + clhs128;
lhs(7,7)=DN(2,0)*clhs64 + DN(2,1)*clhs81 + clhs127 + clhs130*clhs8;
lhs(7,8)=DN(2,1)*clhs129;
lhs(8,0)=DN(2,0)*clhs85 + clhs66;
lhs(8,1)=DN(2,1)*clhs85 + clhs83;
lhs(8,2)=clhs93;
lhs(8,3)=DN(2,0)*clhs88 + clhs109;
lhs(8,4)=DN(2,1)*clhs88 + clhs115;
lhs(8,5)=clhs118;
lhs(8,6)=DN(2,0)*clhs131;
lhs(8,7)=DN(2,1)*clhs131;
lhs(8,8)=clhs125*clhs15 + clhs126*clhs87 + clhs130*clhs15;


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double clhs7 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs9 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs10 = clhs6*stab_c2*sqrt(pow(clhs7, 2) + pow(clhs8, 2) + pow(clhs9, 2));
const double clhs11 = clhs10*h/stab_c1 + mu;
const double clhs12 = DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs13 = N[0]*clhs6;
const double clhs14 = pow(N[0], 2);
const double clhs15 = bdf0*clhs6;
const double clhs16 = N[0]*bdf0;
const double clhs17 = clhs12 + clhs16;
const double clhs18 = 1.0/(clhs10/h + clhs6*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double clhs19 = clhs18*pow(clhs6, 2);
const double clhs20 = clhs12*clhs19;
const double clhs21 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs22 = clhs19*clhs21;
const double clhs23 = N[0]*clhs22;
const double clhs24 = clhs12*clhs13 + clhs14*clhs15 + clhs17*clhs20 + clhs17*clhs23;
const double clhs25 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs26 = C(1,3)*DN(0,1);
const double clhs27 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs26;
const double clhs28 = C(3,5)*DN(0,0);
const double clhs29 = C(4,5)*DN(0,2);
const double clhs30 = C(1,5)*DN(0,1) + clhs28 + clhs29;
const double clhs31 = DN(0,0)*clhs11;
const double clhs32 = DN(0,1)*clhs31;
const double clhs33 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs34 = C(3,4)*DN(0,1);
const double clhs35 = C(2,3)*DN(0,2) + clhs28 + clhs34;
const double clhs36 = C(2,5)*DN(0,2);
const double clhs37 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs36;
const double clhs38 = DN(0,2)*clhs31;
const double clhs39 = 1/(clhs6*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2));
const double clhs40 = clhs16*clhs39;
const double clhs41 = clhs11*clhs40;
const double clhs42 = clhs18*clhs21;
const double clhs43 = clhs18*clhs6;
const double clhs44 = clhs12*clhs43;
const double clhs45 = -N[0] + clhs13*clhs42 + clhs41 + clhs44;
const double clhs46 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs47 = C(0,3)*DN(1,0);
const double clhs48 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs47;
const double clhs49 = C(0,5)*DN(1,0);
const double clhs50 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs49;
const double clhs51 = N[1]*clhs6;
const double clhs52 = clhs16*clhs51;
const double clhs53 = DN(1,0)*clhs31 + clhs52;
const double clhs54 = DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs55 = N[1]*bdf0;
const double clhs56 = clhs54 + clhs55;
const double clhs57 = clhs13*clhs54 + clhs20*clhs56 + clhs23*clhs56;
const double clhs58 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs47;
const double clhs59 = C(1,3)*DN(1,1);
const double clhs60 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs59;
const double clhs61 = C(3,5)*DN(1,0);
const double clhs62 = C(4,5)*DN(1,2);
const double clhs63 = C(1,5)*DN(1,1) + clhs61 + clhs62;
const double clhs64 = DN(1,1)*clhs31;
const double clhs65 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs49;
const double clhs66 = C(3,4)*DN(1,1);
const double clhs67 = C(2,3)*DN(1,2) + clhs61 + clhs66;
const double clhs68 = C(2,5)*DN(1,2);
const double clhs69 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs68;
const double clhs70 = DN(1,2)*clhs31;
const double clhs71 = DN(0,0)*N[1];
const double clhs72 = clhs39*clhs55;
const double clhs73 = DN(1,0)*N[0];
const double clhs74 = clhs21*clhs43;
const double clhs75 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs76 = C(0,3)*DN(2,0);
const double clhs77 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs76;
const double clhs78 = C(0,5)*DN(2,0);
const double clhs79 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs78;
const double clhs80 = N[2]*clhs6;
const double clhs81 = clhs16*clhs80;
const double clhs82 = DN(2,0)*clhs31 + clhs81;
const double clhs83 = DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs84 = N[2]*bdf0;
const double clhs85 = clhs83 + clhs84;
const double clhs86 = clhs13*clhs83 + clhs20*clhs85 + clhs23*clhs85;
const double clhs87 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs76;
const double clhs88 = C(1,3)*DN(2,1);
const double clhs89 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs88;
const double clhs90 = C(3,5)*DN(2,0);
const double clhs91 = C(4,5)*DN(2,2);
const double clhs92 = C(1,5)*DN(2,1) + clhs90 + clhs91;
const double clhs93 = DN(2,1)*clhs31;
const double clhs94 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs78;
const double clhs95 = C(3,4)*DN(2,1);
const double clhs96 = C(2,3)*DN(2,2) + clhs90 + clhs95;
const double clhs97 = C(2,5)*DN(2,2);
const double clhs98 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs97;
const double clhs99 = DN(2,2)*clhs31;
const double clhs100 = DN(0,0)*N[2];
const double clhs101 = clhs39*clhs84;
const double clhs102 = DN(2,0)*N[0];
const double clhs103 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs104 = C(0,3)*DN(3,0);
const double clhs105 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs104;
const double clhs106 = C(0,5)*DN(3,0);
const double clhs107 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs106;
const double clhs108 = N[3]*clhs6;
const double clhs109 = clhs108*clhs16;
const double clhs110 = DN(3,0)*clhs31 + clhs109;
const double clhs111 = DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs112 = N[3]*bdf0;
const double clhs113 = clhs111 + clhs112;
const double clhs114 = clhs111*clhs13 + clhs113*clhs20 + clhs113*clhs23;
const double clhs115 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs104;
const double clhs116 = C(1,3)*DN(3,1);
const double clhs117 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs116;
const double clhs118 = C(3,5)*DN(3,0);
const double clhs119 = C(4,5)*DN(3,2);
const double clhs120 = C(1,5)*DN(3,1) + clhs118 + clhs119;
const double clhs121 = DN(3,1)*clhs31;
const double clhs122 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs106;
const double clhs123 = C(3,4)*DN(3,1);
const double clhs124 = C(2,3)*DN(3,2) + clhs118 + clhs123;
const double clhs125 = C(2,5)*DN(3,2);
const double clhs126 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs125;
const double clhs127 = DN(3,2)*clhs31;
const double clhs128 = DN(0,0)*N[3];
const double clhs129 = clhs112*clhs39;
const double clhs130 = DN(3,0)*N[0];
const double clhs131 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs26;
const double clhs132 = C(0,4)*DN(0,0) + clhs29 + clhs34;
const double clhs133 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs134 = C(1,4)*DN(0,1);
const double clhs135 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs134;
const double clhs136 = pow(DN(0,1), 2);
const double clhs137 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs134;
const double clhs138 = C(2,4)*DN(0,2);
const double clhs139 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs138;
const double clhs140 = DN(0,1)*clhs11;
const double clhs141 = DN(0,2)*clhs140;
const double clhs142 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs59;
const double clhs143 = C(0,4)*DN(1,0) + clhs62 + clhs66;
const double clhs144 = DN(1,0)*clhs140;
const double clhs145 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs146 = C(1,4)*DN(1,1);
const double clhs147 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs146;
const double clhs148 = DN(1,1)*clhs140;
const double clhs149 = clhs52 + clhs57;
const double clhs150 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs146;
const double clhs151 = C(2,4)*DN(1,2);
const double clhs152 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs151;
const double clhs153 = DN(1,2)*clhs140;
const double clhs154 = DN(0,1)*N[1];
const double clhs155 = DN(1,1)*N[0];
const double clhs156 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs88;
const double clhs157 = C(0,4)*DN(2,0) + clhs91 + clhs95;
const double clhs158 = DN(2,0)*clhs140;
const double clhs159 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs160 = C(1,4)*DN(2,1);
const double clhs161 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs160;
const double clhs162 = DN(2,1)*clhs140;
const double clhs163 = clhs81 + clhs86;
const double clhs164 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs160;
const double clhs165 = C(2,4)*DN(2,2);
const double clhs166 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs165;
const double clhs167 = DN(2,2)*clhs140;
const double clhs168 = DN(0,1)*N[2];
const double clhs169 = DN(2,1)*N[0];
const double clhs170 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs116;
const double clhs171 = C(0,4)*DN(3,0) + clhs119 + clhs123;
const double clhs172 = DN(3,0)*clhs140;
const double clhs173 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs174 = C(1,4)*DN(3,1);
const double clhs175 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs174;
const double clhs176 = DN(3,1)*clhs140;
const double clhs177 = clhs109 + clhs114;
const double clhs178 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs174;
const double clhs179 = C(2,4)*DN(3,2);
const double clhs180 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs179;
const double clhs181 = DN(3,2)*clhs140;
const double clhs182 = DN(0,1)*N[3];
const double clhs183 = DN(3,1)*N[0];
const double clhs184 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs36;
const double clhs185 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs138;
const double clhs186 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs187 = pow(DN(0,2), 2);
const double clhs188 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs68;
const double clhs189 = DN(0,2)*clhs11;
const double clhs190 = DN(1,0)*clhs189;
const double clhs191 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs151;
const double clhs192 = DN(1,1)*clhs189;
const double clhs193 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs194 = DN(1,2)*clhs189;
const double clhs195 = DN(0,2)*N[1];
const double clhs196 = DN(1,2)*N[0];
const double clhs197 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs97;
const double clhs198 = DN(2,0)*clhs189;
const double clhs199 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs165;
const double clhs200 = DN(2,1)*clhs189;
const double clhs201 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs202 = DN(2,2)*clhs189;
const double clhs203 = DN(0,2)*N[2];
const double clhs204 = DN(2,2)*N[0];
const double clhs205 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs125;
const double clhs206 = DN(3,0)*clhs189;
const double clhs207 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs179;
const double clhs208 = DN(3,1)*clhs189;
const double clhs209 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs210 = DN(3,2)*clhs189;
const double clhs211 = DN(0,2)*N[3];
const double clhs212 = DN(3,2)*N[0];
const double clhs213 = clhs17*clhs43;
const double clhs214 = N[0] + clhs213;
const double clhs215 = bdf0*clhs39;
const double clhs216 = clhs43*clhs56;
const double clhs217 = DN(0,0)*clhs18;
const double clhs218 = DN(0,1)*clhs18;
const double clhs219 = DN(0,2)*clhs18;
const double clhs220 = DN(1,0)*clhs217 + DN(1,1)*clhs218 + DN(1,2)*clhs219 + N[1]*clhs40;
const double clhs221 = clhs43*clhs85;
const double clhs222 = DN(2,0)*clhs217 + DN(2,1)*clhs218 + DN(2,2)*clhs219 + N[2]*clhs40;
const double clhs223 = clhs113*clhs43;
const double clhs224 = DN(3,0)*clhs217 + DN(3,1)*clhs218 + DN(3,2)*clhs219 + N[3]*clhs40;
const double clhs225 = clhs19*clhs54;
const double clhs226 = N[1]*clhs22;
const double clhs227 = clhs12*clhs51 + clhs17*clhs225 + clhs17*clhs226;
const double clhs228 = DN(1,0)*clhs11;
const double clhs229 = clhs43*clhs54;
const double clhs230 = pow(DN(1,0), 2);
const double clhs231 = pow(N[1], 2);
const double clhs232 = clhs15*clhs231 + clhs225*clhs56 + clhs226*clhs56 + clhs51*clhs54;
const double clhs233 = DN(1,1)*clhs228;
const double clhs234 = DN(1,2)*clhs228;
const double clhs235 = clhs11*clhs72;
const double clhs236 = -N[1] + clhs229 + clhs235 + clhs42*clhs51;
const double clhs237 = clhs55*clhs80;
const double clhs238 = DN(2,0)*clhs228 + clhs237;
const double clhs239 = clhs225*clhs85 + clhs226*clhs85 + clhs51*clhs83;
const double clhs240 = DN(2,1)*clhs228;
const double clhs241 = DN(2,2)*clhs228;
const double clhs242 = DN(1,0)*N[2];
const double clhs243 = DN(2,0)*N[1];
const double clhs244 = clhs108*clhs55;
const double clhs245 = DN(3,0)*clhs228 + clhs244;
const double clhs246 = clhs111*clhs51 + clhs113*clhs225 + clhs113*clhs226;
const double clhs247 = DN(3,1)*clhs228;
const double clhs248 = DN(3,2)*clhs228;
const double clhs249 = DN(1,0)*N[3];
const double clhs250 = DN(3,0)*N[1];
const double clhs251 = clhs227 + clhs52;
const double clhs252 = DN(1,1)*clhs11;
const double clhs253 = pow(DN(1,1), 2);
const double clhs254 = DN(1,2)*clhs252;
const double clhs255 = DN(2,0)*clhs252;
const double clhs256 = DN(2,1)*clhs252;
const double clhs257 = clhs237 + clhs239;
const double clhs258 = DN(2,2)*clhs252;
const double clhs259 = DN(1,1)*N[2];
const double clhs260 = DN(2,1)*N[1];
const double clhs261 = DN(3,0)*clhs252;
const double clhs262 = DN(3,1)*clhs252;
const double clhs263 = clhs244 + clhs246;
const double clhs264 = DN(3,2)*clhs252;
const double clhs265 = DN(1,1)*N[3];
const double clhs266 = DN(3,1)*N[1];
const double clhs267 = DN(1,2)*clhs11;
const double clhs268 = pow(DN(1,2), 2);
const double clhs269 = DN(2,0)*clhs267;
const double clhs270 = DN(2,1)*clhs267;
const double clhs271 = DN(2,2)*clhs267;
const double clhs272 = DN(1,2)*N[2];
const double clhs273 = DN(2,2)*N[1];
const double clhs274 = DN(3,0)*clhs267;
const double clhs275 = DN(3,1)*clhs267;
const double clhs276 = DN(3,2)*clhs267;
const double clhs277 = DN(1,2)*N[3];
const double clhs278 = DN(3,2)*N[1];
const double clhs279 = N[1] + clhs216;
const double clhs280 = DN(1,0)*clhs18;
const double clhs281 = DN(1,1)*clhs18;
const double clhs282 = DN(1,2)*clhs18;
const double clhs283 = DN(2,0)*clhs280 + DN(2,1)*clhs281 + DN(2,2)*clhs282 + N[2]*clhs72;
const double clhs284 = DN(3,0)*clhs280 + DN(3,1)*clhs281 + DN(3,2)*clhs282 + N[3]*clhs72;
const double clhs285 = clhs19*clhs83;
const double clhs286 = N[2]*clhs22;
const double clhs287 = clhs12*clhs80 + clhs17*clhs285 + clhs17*clhs286;
const double clhs288 = DN(2,0)*clhs11;
const double clhs289 = clhs43*clhs83;
const double clhs290 = clhs285*clhs56 + clhs286*clhs56 + clhs54*clhs80;
const double clhs291 = pow(DN(2,0), 2);
const double clhs292 = pow(N[2], 2);
const double clhs293 = clhs15*clhs292 + clhs285*clhs85 + clhs286*clhs85 + clhs80*clhs83;
const double clhs294 = DN(2,1)*clhs288;
const double clhs295 = DN(2,2)*clhs288;
const double clhs296 = clhs101*clhs11;
const double clhs297 = -N[2] + clhs289 + clhs296 + clhs42*clhs80;
const double clhs298 = clhs108*clhs84;
const double clhs299 = DN(3,0)*clhs288 + clhs298;
const double clhs300 = clhs111*clhs80 + clhs113*clhs285 + clhs113*clhs286;
const double clhs301 = DN(3,1)*clhs288;
const double clhs302 = DN(3,2)*clhs288;
const double clhs303 = DN(2,0)*N[3];
const double clhs304 = DN(3,0)*N[2];
const double clhs305 = clhs287 + clhs81;
const double clhs306 = DN(2,1)*clhs11;
const double clhs307 = clhs237 + clhs290;
const double clhs308 = pow(DN(2,1), 2);
const double clhs309 = DN(2,2)*clhs306;
const double clhs310 = DN(3,0)*clhs306;
const double clhs311 = DN(3,1)*clhs306;
const double clhs312 = clhs298 + clhs300;
const double clhs313 = DN(3,2)*clhs306;
const double clhs314 = DN(2,1)*N[3];
const double clhs315 = DN(3,1)*N[2];
const double clhs316 = DN(2,2)*clhs11;
const double clhs317 = pow(DN(2,2), 2);
const double clhs318 = DN(3,0)*clhs316;
const double clhs319 = DN(3,1)*clhs316;
const double clhs320 = DN(3,2)*clhs316;
const double clhs321 = DN(2,2)*N[3];
const double clhs322 = DN(3,2)*N[2];
const double clhs323 = N[2] + clhs221;
const double clhs324 = DN(2,0)*DN(3,0)*clhs18 + DN(2,1)*DN(3,1)*clhs18 + DN(2,2)*DN(3,2)*clhs18 + N[3]*clhs101;
const double clhs325 = clhs111*clhs19;
const double clhs326 = N[3]*clhs22;
const double clhs327 = clhs108*clhs12 + clhs17*clhs325 + clhs17*clhs326;
const double clhs328 = DN(3,0)*clhs11;
const double clhs329 = clhs111*clhs43;
const double clhs330 = clhs108*clhs54 + clhs325*clhs56 + clhs326*clhs56;
const double clhs331 = clhs108*clhs83 + clhs325*clhs85 + clhs326*clhs85;
const double clhs332 = pow(DN(3,0), 2);
const double clhs333 = pow(N[3], 2);
const double clhs334 = clhs108*clhs111 + clhs113*clhs325 + clhs113*clhs326 + clhs15*clhs333;
const double clhs335 = DN(3,1)*clhs328;
const double clhs336 = DN(3,2)*clhs328;
const double clhs337 = -N[3] + clhs108*clhs42 + clhs11*clhs129 + clhs329;
const double clhs338 = clhs109 + clhs327;
const double clhs339 = DN(3,1)*clhs11;
const double clhs340 = clhs244 + clhs330;
const double clhs341 = clhs298 + clhs331;
const double clhs342 = pow(DN(3,1), 2);
const double clhs343 = DN(3,2)*clhs339;
const double clhs344 = pow(DN(3,2), 2);
const double clhs345 = N[3] + clhs223;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs24;
lhs(0,1)=DN(0,0)*clhs25 + DN(0,1)*clhs27 + DN(0,2)*clhs30 + clhs32;
lhs(0,2)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + DN(0,2)*clhs37 + clhs38;
lhs(0,3)=DN(0,0)*clhs45;
lhs(0,4)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs50 + clhs53 + clhs57;
lhs(0,5)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs63 + clhs64;
lhs(0,6)=DN(0,0)*clhs65 + DN(0,1)*clhs67 + DN(0,2)*clhs69 + clhs70;
lhs(0,7)=DN(1,0)*clhs44 + clhs31*clhs72 - clhs71 + clhs73*clhs74;
lhs(0,8)=DN(0,0)*clhs75 + DN(0,1)*clhs77 + DN(0,2)*clhs79 + clhs82 + clhs86;
lhs(0,9)=DN(0,0)*clhs87 + DN(0,1)*clhs89 + DN(0,2)*clhs92 + clhs93;
lhs(0,10)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs99;
lhs(0,11)=DN(2,0)*clhs44 - clhs100 + clhs101*clhs31 + clhs102*clhs74;
lhs(0,12)=DN(0,0)*clhs103 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + clhs110 + clhs114;
lhs(0,13)=DN(0,0)*clhs115 + DN(0,1)*clhs117 + DN(0,2)*clhs120 + clhs121;
lhs(0,14)=DN(0,0)*clhs122 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs127;
lhs(0,15)=DN(3,0)*clhs44 - clhs128 + clhs129*clhs31 + clhs130*clhs74;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs32;
lhs(1,1)=DN(0,0)*clhs27 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs11*clhs136 + clhs24;
lhs(1,2)=DN(0,0)*clhs35 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs141;
lhs(1,3)=DN(0,1)*clhs45;
lhs(1,4)=DN(0,0)*clhs48 + DN(0,1)*clhs142 + DN(0,2)*clhs143 + clhs144;
lhs(1,5)=DN(0,0)*clhs60 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148 + clhs149;
lhs(1,6)=DN(0,0)*clhs67 + DN(0,1)*clhs150 + DN(0,2)*clhs152 + clhs153;
lhs(1,7)=DN(1,1)*clhs44 + clhs140*clhs72 - clhs154 + clhs155*clhs74;
lhs(1,8)=DN(0,0)*clhs77 + DN(0,1)*clhs156 + DN(0,2)*clhs157 + clhs158;
lhs(1,9)=DN(0,0)*clhs89 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162 + clhs163;
lhs(1,10)=DN(0,0)*clhs96 + DN(0,1)*clhs164 + DN(0,2)*clhs166 + clhs167;
lhs(1,11)=DN(2,1)*clhs44 + clhs101*clhs140 - clhs168 + clhs169*clhs74;
lhs(1,12)=DN(0,0)*clhs105 + DN(0,1)*clhs170 + DN(0,2)*clhs171 + clhs172;
lhs(1,13)=DN(0,0)*clhs117 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176 + clhs177;
lhs(1,14)=DN(0,0)*clhs124 + DN(0,1)*clhs178 + DN(0,2)*clhs180 + clhs181;
lhs(1,15)=DN(3,1)*clhs44 + clhs129*clhs140 - clhs182 + clhs183*clhs74;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs132 + DN(0,2)*clhs184 + clhs38;
lhs(2,1)=DN(0,0)*clhs30 + DN(0,1)*clhs135 + DN(0,2)*clhs185 + clhs141;
lhs(2,2)=DN(0,0)*clhs37 + DN(0,1)*clhs139 + DN(0,2)*clhs186 + clhs11*clhs187 + clhs24;
lhs(2,3)=DN(0,2)*clhs45;
lhs(2,4)=DN(0,0)*clhs50 + DN(0,1)*clhs143 + DN(0,2)*clhs188 + clhs190;
lhs(2,5)=DN(0,0)*clhs63 + DN(0,1)*clhs147 + DN(0,2)*clhs191 + clhs192;
lhs(2,6)=DN(0,0)*clhs69 + DN(0,1)*clhs152 + DN(0,2)*clhs193 + clhs149 + clhs194;
lhs(2,7)=DN(1,2)*clhs44 + clhs189*clhs72 - clhs195 + clhs196*clhs74;
lhs(2,8)=DN(0,0)*clhs79 + DN(0,1)*clhs157 + DN(0,2)*clhs197 + clhs198;
lhs(2,9)=DN(0,0)*clhs92 + DN(0,1)*clhs161 + DN(0,2)*clhs199 + clhs200;
lhs(2,10)=DN(0,0)*clhs98 + DN(0,1)*clhs166 + DN(0,2)*clhs201 + clhs163 + clhs202;
lhs(2,11)=DN(2,2)*clhs44 + clhs101*clhs189 - clhs203 + clhs204*clhs74;
lhs(2,12)=DN(0,0)*clhs107 + DN(0,1)*clhs171 + DN(0,2)*clhs205 + clhs206;
lhs(2,13)=DN(0,0)*clhs120 + DN(0,1)*clhs175 + DN(0,2)*clhs207 + clhs208;
lhs(2,14)=DN(0,0)*clhs126 + DN(0,1)*clhs180 + DN(0,2)*clhs209 + clhs177 + clhs210;
lhs(2,15)=DN(3,2)*clhs44 + clhs129*clhs189 - clhs211 + clhs212*clhs74;
lhs(3,0)=DN(0,0)*clhs214;
lhs(3,1)=DN(0,1)*clhs214;
lhs(3,2)=DN(0,2)*clhs214;
lhs(3,3)=clhs136*clhs18 + clhs14*clhs215 + clhs18*clhs187 + clhs18*clhs5;
lhs(3,4)=DN(0,0)*clhs216 + clhs73;
lhs(3,5)=DN(0,1)*clhs216 + clhs155;
lhs(3,6)=DN(0,2)*clhs216 + clhs196;
lhs(3,7)=clhs220;
lhs(3,8)=DN(0,0)*clhs221 + clhs102;
lhs(3,9)=DN(0,1)*clhs221 + clhs169;
lhs(3,10)=DN(0,2)*clhs221 + clhs204;
lhs(3,11)=clhs222;
lhs(3,12)=DN(0,0)*clhs223 + clhs130;
lhs(3,13)=DN(0,1)*clhs223 + clhs183;
lhs(3,14)=DN(0,2)*clhs223 + clhs212;
lhs(3,15)=clhs224;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs227 + clhs53;
lhs(4,1)=DN(1,0)*clhs25 + DN(1,1)*clhs27 + DN(1,2)*clhs30 + clhs144;
lhs(4,2)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + DN(1,2)*clhs37 + clhs190;
lhs(4,3)=DN(0,0)*clhs229 + clhs228*clhs40 + clhs71*clhs74 - clhs73;
lhs(4,4)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs50 + clhs11*clhs230 + clhs232;
lhs(4,5)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs63 + clhs233;
lhs(4,6)=DN(1,0)*clhs65 + DN(1,1)*clhs67 + DN(1,2)*clhs69 + clhs234;
lhs(4,7)=DN(1,0)*clhs236;
lhs(4,8)=DN(1,0)*clhs75 + DN(1,1)*clhs77 + DN(1,2)*clhs79 + clhs238 + clhs239;
lhs(4,9)=DN(1,0)*clhs87 + DN(1,1)*clhs89 + DN(1,2)*clhs92 + clhs240;
lhs(4,10)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs241;
lhs(4,11)=DN(2,0)*clhs229 + clhs101*clhs228 - clhs242 + clhs243*clhs74;
lhs(4,12)=DN(1,0)*clhs103 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + clhs245 + clhs246;
lhs(4,13)=DN(1,0)*clhs115 + DN(1,1)*clhs117 + DN(1,2)*clhs120 + clhs247;
lhs(4,14)=DN(1,0)*clhs122 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs248;
lhs(4,15)=DN(3,0)*clhs229 + clhs129*clhs228 - clhs249 + clhs250*clhs74;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs64;
lhs(5,1)=DN(1,0)*clhs27 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs148 + clhs251;
lhs(5,2)=DN(1,0)*clhs35 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs192;
lhs(5,3)=DN(0,1)*clhs229 + clhs154*clhs74 - clhs155 + clhs252*clhs40;
lhs(5,4)=DN(1,0)*clhs48 + DN(1,1)*clhs142 + DN(1,2)*clhs143 + clhs233;
lhs(5,5)=DN(1,0)*clhs60 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs11*clhs253 + clhs232;
lhs(5,6)=DN(1,0)*clhs67 + DN(1,1)*clhs150 + DN(1,2)*clhs152 + clhs254;
lhs(5,7)=DN(1,1)*clhs236;
lhs(5,8)=DN(1,0)*clhs77 + DN(1,1)*clhs156 + DN(1,2)*clhs157 + clhs255;
lhs(5,9)=DN(1,0)*clhs89 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs256 + clhs257;
lhs(5,10)=DN(1,0)*clhs96 + DN(1,1)*clhs164 + DN(1,2)*clhs166 + clhs258;
lhs(5,11)=DN(2,1)*clhs229 + clhs101*clhs252 - clhs259 + clhs260*clhs74;
lhs(5,12)=DN(1,0)*clhs105 + DN(1,1)*clhs170 + DN(1,2)*clhs171 + clhs261;
lhs(5,13)=DN(1,0)*clhs117 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs262 + clhs263;
lhs(5,14)=DN(1,0)*clhs124 + DN(1,1)*clhs178 + DN(1,2)*clhs180 + clhs264;
lhs(5,15)=DN(3,1)*clhs229 + clhs129*clhs252 - clhs265 + clhs266*clhs74;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs132 + DN(1,2)*clhs184 + clhs70;
lhs(6,1)=DN(1,0)*clhs30 + DN(1,1)*clhs135 + DN(1,2)*clhs185 + clhs153;
lhs(6,2)=DN(1,0)*clhs37 + DN(1,1)*clhs139 + DN(1,2)*clhs186 + clhs194 + clhs251;
lhs(6,3)=DN(0,2)*clhs229 + clhs195*clhs74 - clhs196 + clhs267*clhs40;
lhs(6,4)=DN(1,0)*clhs50 + DN(1,1)*clhs143 + DN(1,2)*clhs188 + clhs234;
lhs(6,5)=DN(1,0)*clhs63 + DN(1,1)*clhs147 + DN(1,2)*clhs191 + clhs254;
lhs(6,6)=DN(1,0)*clhs69 + DN(1,1)*clhs152 + DN(1,2)*clhs193 + clhs11*clhs268 + clhs232;
lhs(6,7)=DN(1,2)*clhs236;
lhs(6,8)=DN(1,0)*clhs79 + DN(1,1)*clhs157 + DN(1,2)*clhs197 + clhs269;
lhs(6,9)=DN(1,0)*clhs92 + DN(1,1)*clhs161 + DN(1,2)*clhs199 + clhs270;
lhs(6,10)=DN(1,0)*clhs98 + DN(1,1)*clhs166 + DN(1,2)*clhs201 + clhs257 + clhs271;
lhs(6,11)=DN(2,2)*clhs229 + clhs101*clhs267 - clhs272 + clhs273*clhs74;
lhs(6,12)=DN(1,0)*clhs107 + DN(1,1)*clhs171 + DN(1,2)*clhs205 + clhs274;
lhs(6,13)=DN(1,0)*clhs120 + DN(1,1)*clhs175 + DN(1,2)*clhs207 + clhs275;
lhs(6,14)=DN(1,0)*clhs126 + DN(1,1)*clhs180 + DN(1,2)*clhs209 + clhs263 + clhs276;
lhs(6,15)=DN(3,2)*clhs229 + clhs129*clhs267 - clhs277 + clhs278*clhs74;
lhs(7,0)=DN(1,0)*clhs213 + clhs71;
lhs(7,1)=DN(1,1)*clhs213 + clhs154;
lhs(7,2)=DN(1,2)*clhs213 + clhs195;
lhs(7,3)=clhs220;
lhs(7,4)=DN(1,0)*clhs279;
lhs(7,5)=DN(1,1)*clhs279;
lhs(7,6)=DN(1,2)*clhs279;
lhs(7,7)=clhs18*clhs230 + clhs18*clhs253 + clhs18*clhs268 + clhs215*clhs231;
lhs(7,8)=DN(1,0)*clhs221 + clhs243;
lhs(7,9)=DN(1,1)*clhs221 + clhs260;
lhs(7,10)=DN(1,2)*clhs221 + clhs273;
lhs(7,11)=clhs283;
lhs(7,12)=DN(1,0)*clhs223 + clhs250;
lhs(7,13)=DN(1,1)*clhs223 + clhs266;
lhs(7,14)=DN(1,2)*clhs223 + clhs278;
lhs(7,15)=clhs284;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs287 + clhs82;
lhs(8,1)=DN(2,0)*clhs25 + DN(2,1)*clhs27 + DN(2,2)*clhs30 + clhs158;
lhs(8,2)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + DN(2,2)*clhs37 + clhs198;
lhs(8,3)=DN(0,0)*clhs289 + clhs100*clhs74 - clhs102 + clhs288*clhs40;
lhs(8,4)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs50 + clhs238 + clhs290;
lhs(8,5)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs63 + clhs255;
lhs(8,6)=DN(2,0)*clhs65 + DN(2,1)*clhs67 + DN(2,2)*clhs69 + clhs269;
lhs(8,7)=DN(1,0)*clhs289 + clhs242*clhs74 - clhs243 + clhs288*clhs72;
lhs(8,8)=DN(2,0)*clhs75 + DN(2,1)*clhs77 + DN(2,2)*clhs79 + clhs11*clhs291 + clhs293;
lhs(8,9)=DN(2,0)*clhs87 + DN(2,1)*clhs89 + DN(2,2)*clhs92 + clhs294;
lhs(8,10)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs295;
lhs(8,11)=DN(2,0)*clhs297;
lhs(8,12)=DN(2,0)*clhs103 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs299 + clhs300;
lhs(8,13)=DN(2,0)*clhs115 + DN(2,1)*clhs117 + DN(2,2)*clhs120 + clhs301;
lhs(8,14)=DN(2,0)*clhs122 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs302;
lhs(8,15)=DN(3,0)*clhs289 + clhs129*clhs288 - clhs303 + clhs304*clhs74;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs93;
lhs(9,1)=DN(2,0)*clhs27 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs162 + clhs305;
lhs(9,2)=DN(2,0)*clhs35 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs200;
lhs(9,3)=DN(0,1)*clhs289 + clhs168*clhs74 - clhs169 + clhs306*clhs40;
lhs(9,4)=DN(2,0)*clhs48 + DN(2,1)*clhs142 + DN(2,2)*clhs143 + clhs240;
lhs(9,5)=DN(2,0)*clhs60 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs256 + clhs307;
lhs(9,6)=DN(2,0)*clhs67 + DN(2,1)*clhs150 + DN(2,2)*clhs152 + clhs270;
lhs(9,7)=DN(1,1)*clhs289 + clhs259*clhs74 - clhs260 + clhs306*clhs72;
lhs(9,8)=DN(2,0)*clhs77 + DN(2,1)*clhs156 + DN(2,2)*clhs157 + clhs294;
lhs(9,9)=DN(2,0)*clhs89 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs11*clhs308 + clhs293;
lhs(9,10)=DN(2,0)*clhs96 + DN(2,1)*clhs164 + DN(2,2)*clhs166 + clhs309;
lhs(9,11)=DN(2,1)*clhs297;
lhs(9,12)=DN(2,0)*clhs105 + DN(2,1)*clhs170 + DN(2,2)*clhs171 + clhs310;
lhs(9,13)=DN(2,0)*clhs117 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs311 + clhs312;
lhs(9,14)=DN(2,0)*clhs124 + DN(2,1)*clhs178 + DN(2,2)*clhs180 + clhs313;
lhs(9,15)=DN(3,1)*clhs289 + clhs129*clhs306 - clhs314 + clhs315*clhs74;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs132 + DN(2,2)*clhs184 + clhs99;
lhs(10,1)=DN(2,0)*clhs30 + DN(2,1)*clhs135 + DN(2,2)*clhs185 + clhs167;
lhs(10,2)=DN(2,0)*clhs37 + DN(2,1)*clhs139 + DN(2,2)*clhs186 + clhs202 + clhs305;
lhs(10,3)=DN(0,2)*clhs289 + clhs203*clhs74 - clhs204 + clhs316*clhs40;
lhs(10,4)=DN(2,0)*clhs50 + DN(2,1)*clhs143 + DN(2,2)*clhs188 + clhs241;
lhs(10,5)=DN(2,0)*clhs63 + DN(2,1)*clhs147 + DN(2,2)*clhs191 + clhs258;
lhs(10,6)=DN(2,0)*clhs69 + DN(2,1)*clhs152 + DN(2,2)*clhs193 + clhs271 + clhs307;
lhs(10,7)=DN(1,2)*clhs289 + clhs272*clhs74 - clhs273 + clhs316*clhs72;
lhs(10,8)=DN(2,0)*clhs79 + DN(2,1)*clhs157 + DN(2,2)*clhs197 + clhs295;
lhs(10,9)=DN(2,0)*clhs92 + DN(2,1)*clhs161 + DN(2,2)*clhs199 + clhs309;
lhs(10,10)=DN(2,0)*clhs98 + DN(2,1)*clhs166 + DN(2,2)*clhs201 + clhs11*clhs317 + clhs293;
lhs(10,11)=DN(2,2)*clhs297;
lhs(10,12)=DN(2,0)*clhs107 + DN(2,1)*clhs171 + DN(2,2)*clhs205 + clhs318;
lhs(10,13)=DN(2,0)*clhs120 + DN(2,1)*clhs175 + DN(2,2)*clhs207 + clhs319;
lhs(10,14)=DN(2,0)*clhs126 + DN(2,1)*clhs180 + DN(2,2)*clhs209 + clhs312 + clhs320;
lhs(10,15)=DN(3,2)*clhs289 + clhs129*clhs316 - clhs321 + clhs322*clhs74;
lhs(11,0)=DN(2,0)*clhs213 + clhs100;
lhs(11,1)=DN(2,1)*clhs213 + clhs168;
lhs(11,2)=DN(2,2)*clhs213 + clhs203;
lhs(11,3)=clhs222;
lhs(11,4)=DN(2,0)*clhs216 + clhs242;
lhs(11,5)=DN(2,1)*clhs216 + clhs259;
lhs(11,6)=DN(2,2)*clhs216 + clhs272;
lhs(11,7)=clhs283;
lhs(11,8)=DN(2,0)*clhs323;
lhs(11,9)=DN(2,1)*clhs323;
lhs(11,10)=DN(2,2)*clhs323;
lhs(11,11)=clhs18*clhs291 + clhs18*clhs308 + clhs18*clhs317 + clhs215*clhs292;
lhs(11,12)=DN(2,0)*clhs223 + clhs304;
lhs(11,13)=DN(2,1)*clhs223 + clhs315;
lhs(11,14)=DN(2,2)*clhs223 + clhs322;
lhs(11,15)=clhs324;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs110 + clhs327;
lhs(12,1)=DN(3,0)*clhs25 + DN(3,1)*clhs27 + DN(3,2)*clhs30 + clhs172;
lhs(12,2)=DN(3,0)*clhs33 + DN(3,1)*clhs35 + DN(3,2)*clhs37 + clhs206;
lhs(12,3)=DN(0,0)*clhs329 + clhs128*clhs74 - clhs130 + clhs328*clhs40;
lhs(12,4)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs50 + clhs245 + clhs330;
lhs(12,5)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs63 + clhs261;
lhs(12,6)=DN(3,0)*clhs65 + DN(3,1)*clhs67 + DN(3,2)*clhs69 + clhs274;
lhs(12,7)=DN(1,0)*clhs329 + clhs249*clhs74 - clhs250 + clhs328*clhs72;
lhs(12,8)=DN(3,0)*clhs75 + DN(3,1)*clhs77 + DN(3,2)*clhs79 + clhs299 + clhs331;
lhs(12,9)=DN(3,0)*clhs87 + DN(3,1)*clhs89 + DN(3,2)*clhs92 + clhs310;
lhs(12,10)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs318;
lhs(12,11)=DN(2,0)*clhs329 + clhs101*clhs328 + clhs303*clhs74 - clhs304;
lhs(12,12)=DN(3,0)*clhs103 + DN(3,1)*clhs105 + DN(3,2)*clhs107 + clhs11*clhs332 + clhs334;
lhs(12,13)=DN(3,0)*clhs115 + DN(3,1)*clhs117 + DN(3,2)*clhs120 + clhs335;
lhs(12,14)=DN(3,0)*clhs122 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs336;
lhs(12,15)=DN(3,0)*clhs337;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs121;
lhs(13,1)=DN(3,0)*clhs27 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs176 + clhs338;
lhs(13,2)=DN(3,0)*clhs35 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs208;
lhs(13,3)=DN(0,1)*clhs329 + clhs182*clhs74 - clhs183 + clhs339*clhs40;
lhs(13,4)=DN(3,0)*clhs48 + DN(3,1)*clhs142 + DN(3,2)*clhs143 + clhs247;
lhs(13,5)=DN(3,0)*clhs60 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs262 + clhs340;
lhs(13,6)=DN(3,0)*clhs67 + DN(3,1)*clhs150 + DN(3,2)*clhs152 + clhs275;
lhs(13,7)=DN(1,1)*clhs329 + clhs265*clhs74 - clhs266 + clhs339*clhs72;
lhs(13,8)=DN(3,0)*clhs77 + DN(3,1)*clhs156 + DN(3,2)*clhs157 + clhs301;
lhs(13,9)=DN(3,0)*clhs89 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs311 + clhs341;
lhs(13,10)=DN(3,0)*clhs96 + DN(3,1)*clhs164 + DN(3,2)*clhs166 + clhs319;
lhs(13,11)=DN(2,1)*clhs329 + clhs101*clhs339 + clhs314*clhs74 - clhs315;
lhs(13,12)=DN(3,0)*clhs105 + DN(3,1)*clhs170 + DN(3,2)*clhs171 + clhs335;
lhs(13,13)=DN(3,0)*clhs117 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs11*clhs342 + clhs334;
lhs(13,14)=DN(3,0)*clhs124 + DN(3,1)*clhs178 + DN(3,2)*clhs180 + clhs343;
lhs(13,15)=DN(3,1)*clhs337;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs132 + DN(3,2)*clhs184 + clhs127;
lhs(14,1)=DN(3,0)*clhs30 + DN(3,1)*clhs135 + DN(3,2)*clhs185 + clhs181;
lhs(14,2)=DN(3,0)*clhs37 + DN(3,1)*clhs139 + DN(3,2)*clhs186 + clhs210 + clhs338;
lhs(14,3)=DN(0,2)*clhs329 + DN(3,2)*clhs41 + clhs211*clhs74 - clhs212;
lhs(14,4)=DN(3,0)*clhs50 + DN(3,1)*clhs143 + DN(3,2)*clhs188 + clhs248;
lhs(14,5)=DN(3,0)*clhs63 + DN(3,1)*clhs147 + DN(3,2)*clhs191 + clhs264;
lhs(14,6)=DN(3,0)*clhs69 + DN(3,1)*clhs152 + DN(3,2)*clhs193 + clhs276 + clhs340;
lhs(14,7)=DN(1,2)*clhs329 + DN(3,2)*clhs235 + clhs277*clhs74 - clhs278;
lhs(14,8)=DN(3,0)*clhs79 + DN(3,1)*clhs157 + DN(3,2)*clhs197 + clhs302;
lhs(14,9)=DN(3,0)*clhs92 + DN(3,1)*clhs161 + DN(3,2)*clhs199 + clhs313;
lhs(14,10)=DN(3,0)*clhs98 + DN(3,1)*clhs166 + DN(3,2)*clhs201 + clhs320 + clhs341;
lhs(14,11)=DN(2,2)*clhs329 + DN(3,2)*clhs296 + clhs321*clhs74 - clhs322;
lhs(14,12)=DN(3,0)*clhs107 + DN(3,1)*clhs171 + DN(3,2)*clhs205 + clhs336;
lhs(14,13)=DN(3,0)*clhs120 + DN(3,1)*clhs175 + DN(3,2)*clhs207 + clhs343;
lhs(14,14)=DN(3,0)*clhs126 + DN(3,1)*clhs180 + DN(3,2)*clhs209 + clhs11*clhs344 + clhs334;
lhs(14,15)=DN(3,2)*clhs337;
lhs(15,0)=DN(3,0)*clhs213 + clhs128;
lhs(15,1)=DN(3,1)*clhs213 + clhs182;
lhs(15,2)=DN(3,2)*clhs213 + clhs211;
lhs(15,3)=clhs224;
lhs(15,4)=DN(3,0)*clhs216 + clhs249;
lhs(15,5)=DN(3,1)*clhs216 + clhs265;
lhs(15,6)=DN(3,2)*clhs216 + clhs277;
lhs(15,7)=clhs284;
lhs(15,8)=DN(3,0)*clhs221 + clhs303;
lhs(15,9)=DN(3,1)*clhs221 + clhs314;
lhs(15,10)=DN(3,2)*clhs221 + clhs321;
lhs(15,11)=clhs324;
lhs(15,12)=DN(3,0)*clhs345;
lhs(15,13)=DN(3,1)*clhs345;
lhs(15,14)=DN(3,2)*clhs345;
lhs(15,15)=clhs18*clhs332 + clhs18*clhs342 + clhs18*clhs344 + clhs215*clhs333;


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crhs2 = crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs3 = crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs4 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs7 = crhs1*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs8 = crhs1*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
const double crhs9 = 1.0/crhs1;
const double crhs10 = crhs9*(crhs5*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crhs6*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crhs11 = crhs9*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crhs12 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs13 = crhs12 + crhs4;
const double crhs14 = (crhs8*h/stab_c1 + mu)*(crhs10 + crhs11 + crhs13);
const double crhs15 = 1.0/(crhs1*dyn_tau/dt + crhs8/h + mu*stab_c1/pow(h, 2));
const double crhs16 = crhs15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs2 + crhs3 + crhs7);
const double crhs17 = crhs1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs18 = N[0]*crhs17;
const double crhs19 = crhs1*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
const double crhs20 = crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs21 = crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs22 = crhs1*(crhs12*crhs6 + crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs23 = crhs15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs20 + crhs21 + crhs22);
const double crhs24 = N[1]*crhs17;
const double crhs25 = crhs1*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
const double crhs26 = N[2]*crhs17;
const double crhs27 = crhs1*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs7 - crhs16*crhs18 - crhs16*crhs19;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs14 - DN(0,1)*stress[1] + N[0]*crhs20 - N[0]*crhs21 - N[0]*crhs22 - crhs18*crhs23 - crhs19*crhs23;
rhs[2]=-DN(0,0)*crhs16 - DN(0,1)*crhs23 - N[0]*crhs10 - N[0]*crhs11 - N[0]*crhs13;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs7 - crhs16*crhs24 - crhs16*crhs25;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs14 - DN(1,1)*stress[1] + N[1]*crhs20 - N[1]*crhs21 - N[1]*crhs22 - crhs23*crhs24 - crhs23*crhs25;
rhs[5]=-DN(1,0)*crhs16 - DN(1,1)*crhs23 - N[1]*crhs10 - N[1]*crhs11 - N[1]*crhs13;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs7 - crhs16*crhs26 - crhs16*crhs27;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs14 - DN(2,1)*stress[1] + N[2]*crhs20 - N[2]*crhs21 - N[2]*crhs22 - crhs23*crhs26 - crhs23*crhs27;
rhs[8]=-DN(2,0)*crhs16 - DN(2,1)*crhs23 - N[2]*crhs10 - N[2]*crhs11 - N[2]*crhs13;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crhs2 = crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs3 = crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs7 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs8 = crhs1*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs7*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs9 = crhs1*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2) + pow(crhs7, 2));
const double crhs10 = 1.0/crhs1;
const double crhs11 = crhs10*(crhs5*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crhs6*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]) + crhs7*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]));
const double crhs12 = crhs10*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crhs13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs14 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs15 = crhs13 + crhs14 + crhs4;
const double crhs16 = (crhs9*h/stab_c1 + mu)*(crhs11 + crhs12 + crhs15);
const double crhs17 = 1.0/(crhs1*dyn_tau/dt + crhs9/h + mu*stab_c1/pow(h, 2));
const double crhs18 = crhs17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs2 + crhs3 + crhs8);
const double crhs19 = crhs1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs20 = N[0]*crhs19;
const double crhs21 = crhs1*(DN(0,0)*crhs5 + DN(0,1)*crhs6 + DN(0,2)*crhs7);
const double crhs22 = crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs23 = crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs24 = crhs1*(crhs13*crhs6 + crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs7*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs25 = crhs17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs22 + crhs23 + crhs24);
const double crhs26 = crhs1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs27 = crhs1*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs28 = crhs1*(crhs14*crhs7 + crhs5*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs6*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs29 = crhs17*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs26 + crhs27 + crhs28);
const double crhs30 = N[1]*crhs19;
const double crhs31 = crhs1*(DN(1,0)*crhs5 + DN(1,1)*crhs6 + DN(1,2)*crhs7);
const double crhs32 = N[2]*crhs19;
const double crhs33 = crhs1*(DN(2,0)*crhs5 + DN(2,1)*crhs6 + DN(2,2)*crhs7);
const double crhs34 = N[3]*crhs19;
const double crhs35 = crhs1*(DN(3,0)*crhs5 + DN(3,1)*crhs6 + DN(3,2)*crhs7);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs16 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs8 - crhs18*crhs20 - crhs18*crhs21;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs16 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs22 - N[0]*crhs23 - N[0]*crhs24 - crhs20*crhs25 - crhs21*crhs25;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs16 - DN(0,2)*stress[2] + N[0]*crhs26 - N[0]*crhs27 - N[0]*crhs28 - crhs20*crhs29 - crhs21*crhs29;
rhs[3]=-DN(0,0)*crhs18 - DN(0,1)*crhs25 - DN(0,2)*crhs29 - N[0]*crhs11 - N[0]*crhs12 - N[0]*crhs15;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs16 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs8 - crhs18*crhs30 - crhs18*crhs31;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs16 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs22 - N[1]*crhs23 - N[1]*crhs24 - crhs25*crhs30 - crhs25*crhs31;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs16 - DN(1,2)*stress[2] + N[1]*crhs26 - N[1]*crhs27 - N[1]*crhs28 - crhs29*crhs30 - crhs29*crhs31;
rhs[7]=-DN(1,0)*crhs18 - DN(1,1)*crhs25 - DN(1,2)*crhs29 - N[1]*crhs11 - N[1]*crhs12 - N[1]*crhs15;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs16 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs8 - crhs18*crhs32 - crhs18*crhs33;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs16 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs22 - N[2]*crhs23 - N[2]*crhs24 - crhs25*crhs32 - crhs25*crhs33;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs16 - DN(2,2)*stress[2] + N[2]*crhs26 - N[2]*crhs27 - N[2]*crhs28 - crhs29*crhs32 - crhs29*crhs33;
rhs[11]=-DN(2,0)*crhs18 - DN(2,1)*crhs25 - DN(2,2)*crhs29 - N[2]*crhs11 - N[2]*crhs12 - N[2]*crhs15;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs16 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs2 - N[3]*crhs3 - N[3]*crhs8 - crhs18*crhs34 - crhs18*crhs35;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs16 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs22 - N[3]*crhs23 - N[3]*crhs24 - crhs25*crhs34 - crhs25*crhs35;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs16 - DN(3,2)*stress[2] + N[3]*crhs26 - N[3]*crhs27 - N[3]*crhs28 - crhs29*crhs34 - crhs29*crhs35;
rhs[15]=-DN(3,0)*crhs18 - DN(3,1)*crhs25 - DN(3,2)*crhs29 - N[3]*crhs11 - N[3]*crhs12 - N[3]*crhs15;


    noalias(rRHS) += rData.Weight * rhs;
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