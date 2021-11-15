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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             C(0,2)*DN(0,0);
const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
const double clhs3 =             pow(DN(0,0), 2);
const double clhs4 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double clhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs7 =             clhs4*stab_c2*sqrt(pow(clhs5, 2) + pow(clhs6, 2));
const double clhs8 =             clhs7*h/stab_c1 + mu;
const double clhs9 =             pow(N[0], 2);
const double clhs10 =             bdf0*clhs4;
const double clhs11 =             N[0]*clhs4;
const double clhs12 =             DN(0,0)*clhs5 + DN(0,1)*clhs6;
const double clhs13 =             N[0]*bdf0 + clhs12;
const double clhs14 =             pow(clhs4, 2);
const double clhs15 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs16 =             1.0/(clhs4*dyn_tau/dt + clhs7/h + mu*stab_c1/pow(h, 2));
const double clhs17 =             1.0*N[0]*clhs14*clhs15*clhs16;
const double clhs18 =             1.0*clhs12*clhs14*clhs16;
const double clhs19 =             clhs10*clhs9 + clhs11*clhs12 + clhs13*clhs17 + clhs13*clhs18;
const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
const double clhs21 =             C(1,2)*DN(0,1);
const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,0)*clhs8;
const double clhs24 =             DN(0,1)*clhs23;
const double clhs25 =             pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], -2);
const double clhs26 =             1.0/clhs4;
const double clhs27 =             N[0]*bdf0*clhs25*clhs26;
const double clhs28 =             1.0*clhs15*clhs16;
const double clhs29 =             1.0*clhs16*clhs4;
const double clhs30 =             clhs12*clhs29;
const double clhs31 =             -N[0] + clhs11*clhs28 + clhs27*clhs8 + clhs30;
const double clhs32 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs33 =             C(0,2)*DN(1,0);
const double clhs34 =             C(2,2)*DN(1,1) + clhs33;
const double clhs35 =             DN(1,0)*clhs23;
const double clhs36 =             N[0]*bdf0*clhs4;
const double clhs37 =             N[1]*clhs36;
const double clhs38 =             DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs39 =             N[1]*bdf0;
const double clhs40 =             clhs38 + clhs39;
const double clhs41 =             clhs11*clhs38 + clhs17*clhs40 + clhs18*clhs40 + clhs37;
const double clhs42 =             C(0,1)*DN(1,1) + clhs33;
const double clhs43 =             C(1,2)*DN(1,1);
const double clhs44 =             C(2,2)*DN(1,0) + clhs43;
const double clhs45 =             DN(1,1)*clhs23;
const double clhs46 =             DN(0,0)*N[1];
const double clhs47 =             bdf0*clhs25*clhs26*clhs8;
const double clhs48 =             DN(1,0)*N[0];
const double clhs49 =             1.0*clhs15*clhs16*clhs4;
const double clhs50 =             1.0*DN(1,0)*clhs16*clhs4;
const double clhs51 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs52 =             C(0,2)*DN(2,0);
const double clhs53 =             C(2,2)*DN(2,1) + clhs52;
const double clhs54 =             DN(2,0)*clhs23;
const double clhs55 =             N[2]*clhs36;
const double clhs56 =             DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs57 =             N[2]*bdf0;
const double clhs58 =             clhs56 + clhs57;
const double clhs59 =             clhs11*clhs56 + clhs17*clhs58 + clhs18*clhs58 + clhs55;
const double clhs60 =             C(0,1)*DN(2,1) + clhs52;
const double clhs61 =             C(1,2)*DN(2,1);
const double clhs62 =             C(2,2)*DN(2,0) + clhs61;
const double clhs63 =             DN(2,1)*clhs23;
const double clhs64 =             DN(0,0)*N[2];
const double clhs65 =             DN(2,0)*N[0];
const double clhs66 =             C(0,1)*DN(0,0) + clhs21;
const double clhs67 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs68 =             pow(DN(0,1), 2);
const double clhs69 =             C(0,1)*DN(1,0) + clhs43;
const double clhs70 =             DN(0,1)*clhs8;
const double clhs71 =             DN(1,0)*clhs70;
const double clhs72 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs73 =             DN(1,1)*clhs70;
const double clhs74 =             DN(0,1)*N[1];
const double clhs75 =             DN(1,1)*N[0];
const double clhs76 =             1.0*DN(1,1)*clhs16*clhs4;
const double clhs77 =             C(0,1)*DN(2,0) + clhs61;
const double clhs78 =             DN(2,0)*clhs70;
const double clhs79 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs80 =             DN(2,1)*clhs70;
const double clhs81 =             DN(0,1)*N[2];
const double clhs82 =             DN(2,1)*N[0];
const double clhs83 =             clhs13*clhs29;
const double clhs84 =             N[0] + clhs83;
const double clhs85 =             bdf0*clhs25*clhs26;
const double clhs86 =             1.0*clhs16;
const double clhs87 =             1.0*DN(0,0)*clhs16*clhs4;
const double clhs88 =             1.0*DN(0,1)*clhs16*clhs4;
const double clhs89 =             1.0*DN(0,0)*clhs16;
const double clhs90 =             1.0*DN(0,1)*clhs16;
const double clhs91 =             DN(1,0)*clhs89 + DN(1,1)*clhs90 + N[1]*clhs27;
const double clhs92 =             DN(2,0)*clhs89 + DN(2,1)*clhs90 + N[2]*clhs27;
const double clhs93 =             N[1]*clhs4;
const double clhs94 =             1.0*N[1]*clhs14*clhs15*clhs16;
const double clhs95 =             1.0*clhs14*clhs16*clhs38;
const double clhs96 =             clhs12*clhs93 + clhs13*clhs94 + clhs13*clhs95 + clhs37;
const double clhs97 =             pow(DN(1,0), 2);
const double clhs98 =             pow(N[1], 2);
const double clhs99 =             clhs10*clhs98 + clhs38*clhs93 + clhs40*clhs94 + clhs40*clhs95;
const double clhs100 =             DN(1,0)*clhs8;
const double clhs101 =             DN(1,1)*clhs100;
const double clhs102 =             clhs25*clhs26*clhs8;
const double clhs103 =             clhs29*clhs38;
const double clhs104 =             -N[1] + clhs102*clhs39 + clhs103 + clhs28*clhs93;
const double clhs105 =             DN(2,0)*clhs100;
const double clhs106 =             N[1]*N[2]*bdf0;
const double clhs107 =             clhs106*clhs4;
const double clhs108 =             clhs107 + clhs56*clhs93 + clhs58*clhs94 + clhs58*clhs95;
const double clhs109 =             DN(2,1)*clhs100;
const double clhs110 =             DN(1,0)*N[2];
const double clhs111 =             DN(2,0)*N[1];
const double clhs112 =             pow(DN(1,1), 2);
const double clhs113 =             DN(1,1)*clhs8;
const double clhs114 =             DN(2,0)*clhs113;
const double clhs115 =             DN(2,1)*clhs113;
const double clhs116 =             DN(1,1)*N[2];
const double clhs117 =             DN(2,1)*N[1];
const double clhs118 =             clhs29*clhs40;
const double clhs119 =             N[1] + clhs118;
const double clhs120 =             1.0*DN(1,0)*DN(2,0)*clhs16 + 1.0*DN(1,1)*DN(2,1)*clhs16 + clhs106*clhs25*clhs26;
const double clhs121 =             N[2]*clhs4;
const double clhs122 =             1.0*N[2]*clhs14*clhs15*clhs16;
const double clhs123 =             1.0*clhs14*clhs16*clhs56;
const double clhs124 =             clhs12*clhs121 + clhs122*clhs13 + clhs123*clhs13 + clhs55;
const double clhs125 =             clhs107 + clhs121*clhs38 + clhs122*clhs40 + clhs123*clhs40;
const double clhs126 =             pow(DN(2,0), 2);
const double clhs127 =             pow(N[2], 2);
const double clhs128 =             clhs10*clhs127 + clhs121*clhs56 + clhs122*clhs58 + clhs123*clhs58;
const double clhs129 =             DN(2,0)*DN(2,1)*clhs8;
const double clhs130 =             -N[2] + clhs102*clhs57 + clhs121*clhs28 + clhs29*clhs56;
const double clhs131 =             pow(DN(2,1), 2);
const double clhs132 =             N[2] + clhs29*clhs58;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs8;
            lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
            lhs(0,2)=DN(0,0)*clhs31;
            lhs(0,3)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + clhs35 + clhs41;
            lhs(0,4)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + clhs45;
            lhs(0,5)=clhs12*clhs50 + clhs46*clhs47 - clhs46 + clhs48*clhs49;
            lhs(0,6)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + clhs54 + clhs59;
            lhs(0,7)=DN(0,0)*clhs60 + DN(0,1)*clhs62 + clhs63;
            lhs(0,8)=DN(2,0)*clhs30 + clhs47*clhs64 + clhs49*clhs65 - clhs64;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs66 + clhs24;
            lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs67 + clhs19 + clhs68*clhs8;
            lhs(1,2)=DN(0,1)*clhs31;
            lhs(1,3)=DN(0,0)*clhs34 + DN(0,1)*clhs69 + clhs71;
            lhs(1,4)=DN(0,0)*clhs44 + DN(0,1)*clhs72 + clhs41 + clhs73;
            lhs(1,5)=clhs12*clhs76 + clhs47*clhs74 + clhs49*clhs75 - clhs74;
            lhs(1,6)=DN(0,0)*clhs53 + DN(0,1)*clhs77 + clhs78;
            lhs(1,7)=DN(0,0)*clhs62 + DN(0,1)*clhs79 + clhs59 + clhs80;
            lhs(1,8)=DN(2,1)*clhs30 + clhs47*clhs81 + clhs49*clhs82 - clhs81;
            lhs(2,0)=DN(0,0)*clhs84;
            lhs(2,1)=DN(0,1)*clhs84;
            lhs(2,2)=clhs3*clhs86 + clhs68*clhs86 + clhs85*clhs9;
            lhs(2,3)=clhs40*clhs87 + clhs48;
            lhs(2,4)=clhs40*clhs88 + clhs75;
            lhs(2,5)=clhs91;
            lhs(2,6)=clhs58*clhs87 + clhs65;
            lhs(2,7)=clhs58*clhs88 + clhs82;
            lhs(2,8)=clhs92;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs35 + clhs96;
            lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs71;
            lhs(3,2)=clhs38*clhs87 + clhs46*clhs49 + clhs47*clhs48 - clhs48;
            lhs(3,3)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + clhs8*clhs97 + clhs99;
            lhs(3,4)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + clhs101;
            lhs(3,5)=DN(1,0)*clhs104;
            lhs(3,6)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + clhs105 + clhs108;
            lhs(3,7)=DN(1,0)*clhs60 + DN(1,1)*clhs62 + clhs109;
            lhs(3,8)=DN(2,0)*clhs103 + clhs110*clhs47 - clhs110 + clhs111*clhs49;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs66 + clhs45;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs67 + clhs73 + clhs96;
            lhs(4,2)=clhs38*clhs88 + clhs47*clhs75 + clhs49*clhs74 - clhs75;
            lhs(4,3)=DN(1,0)*clhs34 + DN(1,1)*clhs69 + clhs101;
            lhs(4,4)=DN(1,0)*clhs44 + DN(1,1)*clhs72 + clhs112*clhs8 + clhs99;
            lhs(4,5)=DN(1,1)*clhs104;
            lhs(4,6)=DN(1,0)*clhs53 + DN(1,1)*clhs77 + clhs114;
            lhs(4,7)=DN(1,0)*clhs62 + DN(1,1)*clhs79 + clhs108 + clhs115;
            lhs(4,8)=DN(2,1)*clhs103 + clhs116*clhs47 - clhs116 + clhs117*clhs49;
            lhs(5,0)=clhs13*clhs50 + clhs46;
            lhs(5,1)=clhs13*clhs76 + clhs74;
            lhs(5,2)=clhs91;
            lhs(5,3)=DN(1,0)*clhs119;
            lhs(5,4)=DN(1,1)*clhs119;
            lhs(5,5)=clhs112*clhs86 + clhs85*clhs98 + clhs86*clhs97;
            lhs(5,6)=clhs111 + clhs50*clhs58;
            lhs(5,7)=clhs117 + clhs58*clhs76;
            lhs(5,8)=clhs120;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs124 + clhs54;
            lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs78;
            lhs(6,2)=clhs47*clhs65 + clhs49*clhs64 + clhs56*clhs87 - clhs65;
            lhs(6,3)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + clhs105 + clhs125;
            lhs(6,4)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + clhs114;
            lhs(6,5)=clhs110*clhs49 + clhs111*clhs47 - clhs111 + clhs50*clhs56;
            lhs(6,6)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + clhs126*clhs8 + clhs128;
            lhs(6,7)=DN(2,0)*clhs60 + DN(2,1)*clhs62 + clhs129;
            lhs(6,8)=DN(2,0)*clhs130;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs66 + clhs63;
            lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs67 + clhs124 + clhs80;
            lhs(7,2)=clhs47*clhs82 + clhs49*clhs81 + clhs56*clhs88 - clhs82;
            lhs(7,3)=DN(2,0)*clhs34 + DN(2,1)*clhs69 + clhs109;
            lhs(7,4)=DN(2,0)*clhs44 + DN(2,1)*clhs72 + clhs115 + clhs125;
            lhs(7,5)=clhs116*clhs49 + clhs117*clhs47 - clhs117 + clhs56*clhs76;
            lhs(7,6)=DN(2,0)*clhs53 + DN(2,1)*clhs77 + clhs129;
            lhs(7,7)=DN(2,0)*clhs62 + DN(2,1)*clhs79 + clhs128 + clhs131*clhs8;
            lhs(7,8)=DN(2,1)*clhs130;
            lhs(8,0)=DN(2,0)*clhs83 + clhs64;
            lhs(8,1)=DN(2,1)*clhs83 + clhs81;
            lhs(8,2)=clhs92;
            lhs(8,3)=DN(2,0)*clhs118 + clhs110;
            lhs(8,4)=DN(2,1)*clhs118 + clhs116;
            lhs(8,5)=clhs120;
            lhs(8,6)=DN(2,0)*clhs132;
            lhs(8,7)=DN(2,1)*clhs132;
            lhs(8,8)=clhs126*clhs86 + clhs127*clhs85 + clhs131*clhs86;


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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             C(0,3)*DN(0,0);
const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 =             C(0,5)*DN(0,0);
const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double clhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs10 =             clhs6*stab_c2*sqrt(pow(clhs7, 2) + pow(clhs8, 2) + pow(clhs9, 2));
const double clhs11 =             clhs10*h/stab_c1 + mu;
const double clhs12 =             pow(N[0], 2);
const double clhs13 =             bdf0*clhs6;
const double clhs14 =             N[0]*clhs6;
const double clhs15 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs16 =             N[0]*bdf0 + clhs15;
const double clhs17 =             pow(clhs6, 2);
const double clhs18 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs19 =             1.0/(clhs10/h + clhs6*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double clhs20 =             1.0*N[0]*clhs17*clhs18*clhs19;
const double clhs21 =             1.0*clhs15*clhs17*clhs19;
const double clhs22 =             clhs12*clhs13 + clhs14*clhs15 + clhs16*clhs20 + clhs16*clhs21;
const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs24 =             C(1,3)*DN(0,1);
const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
const double clhs26 =             C(3,5)*DN(0,0);
const double clhs27 =             C(4,5)*DN(0,2);
const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
const double clhs29 =             DN(0,0)*clhs11;
const double clhs30 =             DN(0,1)*clhs29;
const double clhs31 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs32 =             C(3,4)*DN(0,1);
const double clhs33 =             C(2,3)*DN(0,2) + clhs26 + clhs32;
const double clhs34 =             C(2,5)*DN(0,2);
const double clhs35 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
const double clhs36 =             DN(0,2)*clhs29;
const double clhs37 =             pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], -2);
const double clhs38 =             1.0/clhs6;
const double clhs39 =             N[0]*bdf0*clhs37*clhs38;
const double clhs40 =             1.0*clhs18*clhs19;
const double clhs41 =             1.0*clhs19*clhs6;
const double clhs42 =             clhs15*clhs41;
const double clhs43 =             -N[0] + clhs11*clhs39 + clhs14*clhs40 + clhs42;
const double clhs44 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs45 =             C(0,3)*DN(1,0);
const double clhs46 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs45;
const double clhs47 =             C(0,5)*DN(1,0);
const double clhs48 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs47;
const double clhs49 =             DN(1,0)*clhs29;
const double clhs50 =             N[0]*bdf0*clhs6;
const double clhs51 =             N[1]*clhs50;
const double clhs52 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs53 =             N[1]*bdf0 + clhs52;
const double clhs54 =             clhs14*clhs52 + clhs20*clhs53 + clhs21*clhs53 + clhs51;
const double clhs55 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs45;
const double clhs56 =             C(1,3)*DN(1,1);
const double clhs57 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs56;
const double clhs58 =             C(3,5)*DN(1,0);
const double clhs59 =             C(4,5)*DN(1,2);
const double clhs60 =             C(1,5)*DN(1,1) + clhs58 + clhs59;
const double clhs61 =             DN(1,1)*clhs29;
const double clhs62 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs47;
const double clhs63 =             C(3,4)*DN(1,1);
const double clhs64 =             C(2,3)*DN(1,2) + clhs58 + clhs63;
const double clhs65 =             C(2,5)*DN(1,2);
const double clhs66 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs65;
const double clhs67 =             DN(1,2)*clhs29;
const double clhs68 =             DN(0,0)*N[1];
const double clhs69 =             bdf0*clhs11*clhs37*clhs38;
const double clhs70 =             DN(1,0)*N[0];
const double clhs71 =             1.0*clhs18*clhs19*clhs6;
const double clhs72 =             1.0*DN(1,0)*clhs19*clhs6;
const double clhs73 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs74 =             C(0,3)*DN(2,0);
const double clhs75 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs74;
const double clhs76 =             C(0,5)*DN(2,0);
const double clhs77 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs76;
const double clhs78 =             DN(2,0)*clhs29;
const double clhs79 =             N[2]*clhs50;
const double clhs80 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs81 =             N[2]*bdf0;
const double clhs82 =             clhs80 + clhs81;
const double clhs83 =             clhs14*clhs80 + clhs20*clhs82 + clhs21*clhs82 + clhs79;
const double clhs84 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs74;
const double clhs85 =             C(1,3)*DN(2,1);
const double clhs86 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs85;
const double clhs87 =             C(3,5)*DN(2,0);
const double clhs88 =             C(4,5)*DN(2,2);
const double clhs89 =             C(1,5)*DN(2,1) + clhs87 + clhs88;
const double clhs90 =             DN(2,1)*clhs29;
const double clhs91 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs76;
const double clhs92 =             C(3,4)*DN(2,1);
const double clhs93 =             C(2,3)*DN(2,2) + clhs87 + clhs92;
const double clhs94 =             C(2,5)*DN(2,2);
const double clhs95 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs94;
const double clhs96 =             DN(2,2)*clhs29;
const double clhs97 =             DN(0,0)*N[2];
const double clhs98 =             DN(2,0)*N[0];
const double clhs99 =             1.0*DN(2,0)*clhs19*clhs6;
const double clhs100 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs101 =             C(0,3)*DN(3,0);
const double clhs102 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs101;
const double clhs103 =             C(0,5)*DN(3,0);
const double clhs104 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs103;
const double clhs105 =             DN(3,0)*clhs29;
const double clhs106 =             N[3]*clhs50;
const double clhs107 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs108 =             N[3]*bdf0;
const double clhs109 =             clhs107 + clhs108;
const double clhs110 =             clhs106 + clhs107*clhs14 + clhs109*clhs20 + clhs109*clhs21;
const double clhs111 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs101;
const double clhs112 =             C(1,3)*DN(3,1);
const double clhs113 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs112;
const double clhs114 =             C(3,5)*DN(3,0);
const double clhs115 =             C(4,5)*DN(3,2);
const double clhs116 =             C(1,5)*DN(3,1) + clhs114 + clhs115;
const double clhs117 =             DN(3,1)*clhs29;
const double clhs118 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs103;
const double clhs119 =             C(3,4)*DN(3,1);
const double clhs120 =             C(2,3)*DN(3,2) + clhs114 + clhs119;
const double clhs121 =             C(2,5)*DN(3,2);
const double clhs122 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs121;
const double clhs123 =             DN(3,2)*clhs29;
const double clhs124 =             DN(0,0)*N[3];
const double clhs125 =             DN(3,0)*N[0];
const double clhs126 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs127 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs128 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs129 =             C(1,4)*DN(0,1);
const double clhs130 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs129;
const double clhs131 =             pow(DN(0,1), 2);
const double clhs132 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs129;
const double clhs133 =             C(2,4)*DN(0,2);
const double clhs134 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs133;
const double clhs135 =             DN(0,1)*clhs11;
const double clhs136 =             DN(0,2)*clhs135;
const double clhs137 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs56;
const double clhs138 =             C(0,4)*DN(1,0) + clhs59 + clhs63;
const double clhs139 =             DN(1,0)*clhs135;
const double clhs140 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs141 =             C(1,4)*DN(1,1);
const double clhs142 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs141;
const double clhs143 =             DN(1,1)*clhs135;
const double clhs144 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs141;
const double clhs145 =             C(2,4)*DN(1,2);
const double clhs146 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs145;
const double clhs147 =             DN(1,2)*clhs135;
const double clhs148 =             DN(0,1)*N[1];
const double clhs149 =             DN(1,1)*N[0];
const double clhs150 =             1.0*DN(1,1)*clhs19*clhs6;
const double clhs151 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs85;
const double clhs152 =             C(0,4)*DN(2,0) + clhs88 + clhs92;
const double clhs153 =             DN(2,0)*clhs135;
const double clhs154 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs155 =             C(1,4)*DN(2,1);
const double clhs156 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs155;
const double clhs157 =             DN(2,1)*clhs135;
const double clhs158 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs155;
const double clhs159 =             C(2,4)*DN(2,2);
const double clhs160 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs159;
const double clhs161 =             DN(2,2)*clhs135;
const double clhs162 =             DN(0,1)*N[2];
const double clhs163 =             DN(2,1)*N[0];
const double clhs164 =             1.0*DN(2,1)*clhs19*clhs6;
const double clhs165 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs112;
const double clhs166 =             C(0,4)*DN(3,0) + clhs115 + clhs119;
const double clhs167 =             DN(3,0)*clhs135;
const double clhs168 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs169 =             C(1,4)*DN(3,1);
const double clhs170 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs169;
const double clhs171 =             DN(3,1)*clhs135;
const double clhs172 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs169;
const double clhs173 =             C(2,4)*DN(3,2);
const double clhs174 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs173;
const double clhs175 =             DN(3,2)*clhs135;
const double clhs176 =             DN(0,1)*N[3];
const double clhs177 =             DN(3,1)*N[0];
const double clhs178 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs179 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs133;
const double clhs180 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs181 =             pow(DN(0,2), 2);
const double clhs182 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs65;
const double clhs183 =             DN(0,2)*clhs11;
const double clhs184 =             DN(1,0)*clhs183;
const double clhs185 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs145;
const double clhs186 =             DN(1,1)*clhs183;
const double clhs187 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs188 =             DN(1,2)*clhs183;
const double clhs189 =             DN(0,2)*N[1];
const double clhs190 =             DN(1,2)*N[0];
const double clhs191 =             1.0*DN(1,2)*clhs19*clhs6;
const double clhs192 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs94;
const double clhs193 =             DN(2,0)*clhs183;
const double clhs194 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs159;
const double clhs195 =             DN(2,1)*clhs183;
const double clhs196 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs197 =             DN(2,2)*clhs183;
const double clhs198 =             DN(0,2)*N[2];
const double clhs199 =             DN(2,2)*N[0];
const double clhs200 =             1.0*DN(2,2)*clhs19*clhs6;
const double clhs201 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs121;
const double clhs202 =             DN(3,0)*clhs183;
const double clhs203 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs173;
const double clhs204 =             DN(3,1)*clhs183;
const double clhs205 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs206 =             DN(3,2)*clhs183;
const double clhs207 =             DN(0,2)*N[3];
const double clhs208 =             DN(3,2)*N[0];
const double clhs209 =             clhs16*clhs41;
const double clhs210 =             N[0] + clhs209;
const double clhs211 =             bdf0*clhs37*clhs38;
const double clhs212 =             1.0*clhs19;
const double clhs213 =             1.0*DN(0,0)*clhs19*clhs6;
const double clhs214 =             1.0*DN(0,1)*clhs19*clhs6;
const double clhs215 =             1.0*DN(0,2)*clhs19*clhs6;
const double clhs216 =             1.0*DN(0,0)*clhs19;
const double clhs217 =             1.0*DN(0,1)*clhs19;
const double clhs218 =             1.0*DN(0,2)*clhs19;
const double clhs219 =             DN(1,0)*clhs216 + DN(1,1)*clhs217 + DN(1,2)*clhs218 + N[1]*clhs39;
const double clhs220 =             DN(2,0)*clhs216 + DN(2,1)*clhs217 + DN(2,2)*clhs218 + N[2]*clhs39;
const double clhs221 =             DN(3,0)*clhs216 + DN(3,1)*clhs217 + DN(3,2)*clhs218 + N[3]*clhs39;
const double clhs222 =             N[1]*clhs6;
const double clhs223 =             1.0*N[1]*clhs17*clhs18*clhs19;
const double clhs224 =             1.0*clhs17*clhs19*clhs52;
const double clhs225 =             clhs15*clhs222 + clhs16*clhs223 + clhs16*clhs224 + clhs51;
const double clhs226 =             pow(DN(1,0), 2);
const double clhs227 =             pow(N[1], 2);
const double clhs228 =             clhs13*clhs227 + clhs222*clhs52 + clhs223*clhs53 + clhs224*clhs53;
const double clhs229 =             DN(1,0)*clhs11;
const double clhs230 =             DN(1,1)*clhs229;
const double clhs231 =             DN(1,2)*clhs229;
const double clhs232 =             N[1]*bdf0*clhs37*clhs38;
const double clhs233 =             clhs41*clhs52;
const double clhs234 =             -N[1] + clhs11*clhs232 + clhs222*clhs40 + clhs233;
const double clhs235 =             DN(2,0)*clhs229;
const double clhs236 =             N[1]*bdf0*clhs6;
const double clhs237 =             N[2]*clhs236;
const double clhs238 =             clhs222*clhs80 + clhs223*clhs82 + clhs224*clhs82 + clhs237;
const double clhs239 =             DN(2,1)*clhs229;
const double clhs240 =             DN(2,2)*clhs229;
const double clhs241 =             DN(1,0)*N[2];
const double clhs242 =             DN(2,0)*N[1];
const double clhs243 =             DN(3,0)*clhs229;
const double clhs244 =             N[3]*clhs236;
const double clhs245 =             clhs107*clhs222 + clhs109*clhs223 + clhs109*clhs224 + clhs244;
const double clhs246 =             DN(3,1)*clhs229;
const double clhs247 =             DN(3,2)*clhs229;
const double clhs248 =             DN(1,0)*N[3];
const double clhs249 =             DN(3,0)*N[1];
const double clhs250 =             pow(DN(1,1), 2);
const double clhs251 =             DN(1,1)*clhs11;
const double clhs252 =             DN(1,2)*clhs251;
const double clhs253 =             DN(2,0)*clhs251;
const double clhs254 =             DN(2,1)*clhs251;
const double clhs255 =             DN(2,2)*clhs251;
const double clhs256 =             DN(1,1)*N[2];
const double clhs257 =             DN(2,1)*N[1];
const double clhs258 =             DN(3,0)*clhs251;
const double clhs259 =             DN(3,1)*clhs251;
const double clhs260 =             DN(3,2)*clhs251;
const double clhs261 =             DN(1,1)*N[3];
const double clhs262 =             DN(3,1)*N[1];
const double clhs263 =             pow(DN(1,2), 2);
const double clhs264 =             DN(1,2)*clhs11;
const double clhs265 =             DN(2,0)*clhs264;
const double clhs266 =             DN(2,1)*clhs264;
const double clhs267 =             DN(2,2)*clhs264;
const double clhs268 =             DN(1,2)*N[2];
const double clhs269 =             DN(2,2)*N[1];
const double clhs270 =             DN(3,0)*clhs264;
const double clhs271 =             DN(3,1)*clhs264;
const double clhs272 =             DN(3,2)*clhs264;
const double clhs273 =             DN(1,2)*N[3];
const double clhs274 =             DN(3,2)*N[1];
const double clhs275 =             clhs41*clhs53;
const double clhs276 =             N[1] + clhs275;
const double clhs277 =             1.0*DN(1,0)*clhs19;
const double clhs278 =             1.0*DN(1,1)*clhs19;
const double clhs279 =             1.0*DN(1,2)*clhs19;
const double clhs280 =             DN(2,0)*clhs277 + DN(2,1)*clhs278 + DN(2,2)*clhs279 + N[2]*clhs232;
const double clhs281 =             DN(3,0)*clhs277 + DN(3,1)*clhs278 + DN(3,2)*clhs279 + N[3]*clhs232;
const double clhs282 =             N[2]*clhs6;
const double clhs283 =             1.0*N[2]*clhs17*clhs18*clhs19;
const double clhs284 =             1.0*clhs17*clhs19*clhs80;
const double clhs285 =             clhs15*clhs282 + clhs16*clhs283 + clhs16*clhs284 + clhs79;
const double clhs286 =             clhs237 + clhs282*clhs52 + clhs283*clhs53 + clhs284*clhs53;
const double clhs287 =             pow(DN(2,0), 2);
const double clhs288 =             pow(N[2], 2);
const double clhs289 =             clhs13*clhs288 + clhs282*clhs80 + clhs283*clhs82 + clhs284*clhs82;
const double clhs290 =             DN(2,0)*clhs11;
const double clhs291 =             DN(2,1)*clhs290;
const double clhs292 =             DN(2,2)*clhs290;
const double clhs293 =             clhs11*clhs37*clhs38;
const double clhs294 =             clhs41*clhs80;
const double clhs295 =             -N[2] + clhs282*clhs40 + clhs293*clhs81 + clhs294;
const double clhs296 =             DN(3,0)*clhs290;
const double clhs297 =             N[2]*N[3]*bdf0;
const double clhs298 =             clhs297*clhs6;
const double clhs299 =             clhs107*clhs282 + clhs109*clhs283 + clhs109*clhs284 + clhs298;
const double clhs300 =             DN(3,1)*clhs290;
const double clhs301 =             DN(3,2)*clhs290;
const double clhs302 =             DN(2,0)*N[3];
const double clhs303 =             DN(3,0)*N[2];
const double clhs304 =             pow(DN(2,1), 2);
const double clhs305 =             DN(2,1)*clhs11;
const double clhs306 =             DN(2,2)*clhs305;
const double clhs307 =             DN(3,0)*clhs305;
const double clhs308 =             DN(3,1)*clhs305;
const double clhs309 =             DN(3,2)*clhs305;
const double clhs310 =             DN(2,1)*N[3];
const double clhs311 =             DN(3,1)*N[2];
const double clhs312 =             pow(DN(2,2), 2);
const double clhs313 =             DN(2,2)*clhs11;
const double clhs314 =             DN(3,0)*clhs313;
const double clhs315 =             DN(3,1)*clhs313;
const double clhs316 =             DN(3,2)*clhs313;
const double clhs317 =             DN(2,2)*N[3];
const double clhs318 =             DN(3,2)*N[2];
const double clhs319 =             clhs41*clhs82;
const double clhs320 =             N[2] + clhs319;
const double clhs321 =             1.0*DN(2,0)*DN(3,0)*clhs19 + 1.0*DN(2,1)*DN(3,1)*clhs19 + 1.0*DN(2,2)*DN(3,2)*clhs19 + clhs297*clhs37*clhs38;
const double clhs322 =             N[3]*clhs6;
const double clhs323 =             1.0*N[3]*clhs17*clhs18*clhs19;
const double clhs324 =             1.0*clhs107*clhs17*clhs19;
const double clhs325 =             clhs106 + clhs15*clhs322 + clhs16*clhs323 + clhs16*clhs324;
const double clhs326 =             clhs244 + clhs322*clhs52 + clhs323*clhs53 + clhs324*clhs53;
const double clhs327 =             clhs298 + clhs322*clhs80 + clhs323*clhs82 + clhs324*clhs82;
const double clhs328 =             pow(DN(3,0), 2);
const double clhs329 =             pow(N[3], 2);
const double clhs330 =             clhs107*clhs322 + clhs109*clhs323 + clhs109*clhs324 + clhs13*clhs329;
const double clhs331 =             DN(3,0)*clhs11;
const double clhs332 =             DN(3,1)*clhs331;
const double clhs333 =             DN(3,2)*clhs331;
const double clhs334 =             -N[3] + clhs107*clhs41 + clhs108*clhs293 + clhs322*clhs40;
const double clhs335 =             pow(DN(3,1), 2);
const double clhs336 =             DN(3,1)*DN(3,2)*clhs11;
const double clhs337 =             pow(DN(3,2), 2);
const double clhs338 =             N[3] + clhs109*clhs41;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs22;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
            lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
            lhs(0,3)=DN(0,0)*clhs43;
            lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + DN(0,2)*clhs48 + clhs49 + clhs54;
            lhs(0,5)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs61;
            lhs(0,6)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs67;
            lhs(0,7)=clhs15*clhs72 + clhs68*clhs69 - clhs68 + clhs70*clhs71;
            lhs(0,8)=DN(0,0)*clhs73 + DN(0,1)*clhs75 + DN(0,2)*clhs77 + clhs78 + clhs83;
            lhs(0,9)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs89 + clhs90;
            lhs(0,10)=DN(0,0)*clhs91 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs96;
            lhs(0,11)=clhs15*clhs99 + clhs69*clhs97 + clhs71*clhs98 - clhs97;
            lhs(0,12)=DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs105 + clhs110;
            lhs(0,13)=DN(0,0)*clhs111 + DN(0,1)*clhs113 + DN(0,2)*clhs116 + clhs117;
            lhs(0,14)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs123;
            lhs(0,15)=DN(3,0)*clhs42 + clhs124*clhs69 - clhs124 + clhs125*clhs71;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs126 + DN(0,2)*clhs127 + clhs30;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs128 + DN(0,2)*clhs130 + clhs11*clhs131 + clhs22;
            lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs136;
            lhs(1,3)=DN(0,1)*clhs43;
            lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
            lhs(1,5)=DN(0,0)*clhs57 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs143 + clhs54;
            lhs(1,6)=DN(0,0)*clhs64 + DN(0,1)*clhs144 + DN(0,2)*clhs146 + clhs147;
            lhs(1,7)=clhs148*clhs69 - clhs148 + clhs149*clhs71 + clhs15*clhs150;
            lhs(1,8)=DN(0,0)*clhs75 + DN(0,1)*clhs151 + DN(0,2)*clhs152 + clhs153;
            lhs(1,9)=DN(0,0)*clhs86 + DN(0,1)*clhs154 + DN(0,2)*clhs156 + clhs157 + clhs83;
            lhs(1,10)=DN(0,0)*clhs93 + DN(0,1)*clhs158 + DN(0,2)*clhs160 + clhs161;
            lhs(1,11)=clhs15*clhs164 + clhs162*clhs69 - clhs162 + clhs163*clhs71;
            lhs(1,12)=DN(0,0)*clhs102 + DN(0,1)*clhs165 + DN(0,2)*clhs166 + clhs167;
            lhs(1,13)=DN(0,0)*clhs113 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs110 + clhs171;
            lhs(1,14)=DN(0,0)*clhs120 + DN(0,1)*clhs172 + DN(0,2)*clhs174 + clhs175;
            lhs(1,15)=DN(3,1)*clhs42 + clhs176*clhs69 - clhs176 + clhs177*clhs71;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs178 + clhs36;
            lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs130 + DN(0,2)*clhs179 + clhs136;
            lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs134 + DN(0,2)*clhs180 + clhs11*clhs181 + clhs22;
            lhs(2,3)=DN(0,2)*clhs43;
            lhs(2,4)=DN(0,0)*clhs48 + DN(0,1)*clhs138 + DN(0,2)*clhs182 + clhs184;
            lhs(2,5)=DN(0,0)*clhs60 + DN(0,1)*clhs142 + DN(0,2)*clhs185 + clhs186;
            lhs(2,6)=DN(0,0)*clhs66 + DN(0,1)*clhs146 + DN(0,2)*clhs187 + clhs188 + clhs54;
            lhs(2,7)=clhs15*clhs191 + clhs189*clhs69 - clhs189 + clhs190*clhs71;
            lhs(2,8)=DN(0,0)*clhs77 + DN(0,1)*clhs152 + DN(0,2)*clhs192 + clhs193;
            lhs(2,9)=DN(0,0)*clhs89 + DN(0,1)*clhs156 + DN(0,2)*clhs194 + clhs195;
            lhs(2,10)=DN(0,0)*clhs95 + DN(0,1)*clhs160 + DN(0,2)*clhs196 + clhs197 + clhs83;
            lhs(2,11)=clhs15*clhs200 + clhs198*clhs69 - clhs198 + clhs199*clhs71;
            lhs(2,12)=DN(0,0)*clhs104 + DN(0,1)*clhs166 + DN(0,2)*clhs201 + clhs202;
            lhs(2,13)=DN(0,0)*clhs116 + DN(0,1)*clhs170 + DN(0,2)*clhs203 + clhs204;
            lhs(2,14)=DN(0,0)*clhs122 + DN(0,1)*clhs174 + DN(0,2)*clhs205 + clhs110 + clhs206;
            lhs(2,15)=DN(3,2)*clhs42 + clhs207*clhs69 - clhs207 + clhs208*clhs71;
            lhs(3,0)=DN(0,0)*clhs210;
            lhs(3,1)=DN(0,1)*clhs210;
            lhs(3,2)=DN(0,2)*clhs210;
            lhs(3,3)=clhs12*clhs211 + clhs131*clhs212 + clhs181*clhs212 + clhs212*clhs5;
            lhs(3,4)=clhs213*clhs53 + clhs70;
            lhs(3,5)=clhs149 + clhs214*clhs53;
            lhs(3,6)=clhs190 + clhs215*clhs53;
            lhs(3,7)=clhs219;
            lhs(3,8)=clhs213*clhs82 + clhs98;
            lhs(3,9)=clhs163 + clhs214*clhs82;
            lhs(3,10)=clhs199 + clhs215*clhs82;
            lhs(3,11)=clhs220;
            lhs(3,12)=clhs109*clhs213 + clhs125;
            lhs(3,13)=clhs109*clhs214 + clhs177;
            lhs(3,14)=clhs109*clhs215 + clhs208;
            lhs(3,15)=clhs221;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs225 + clhs49;
            lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs139;
            lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs184;
            lhs(4,3)=clhs213*clhs52 + clhs68*clhs71 + clhs69*clhs70 - clhs70;
            lhs(4,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + DN(1,2)*clhs48 + clhs11*clhs226 + clhs228;
            lhs(4,5)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs230;
            lhs(4,6)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs231;
            lhs(4,7)=DN(1,0)*clhs234;
            lhs(4,8)=DN(1,0)*clhs73 + DN(1,1)*clhs75 + DN(1,2)*clhs77 + clhs235 + clhs238;
            lhs(4,9)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs89 + clhs239;
            lhs(4,10)=DN(1,0)*clhs91 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs240;
            lhs(4,11)=clhs241*clhs69 - clhs241 + clhs242*clhs71 + clhs52*clhs99;
            lhs(4,12)=DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs243 + clhs245;
            lhs(4,13)=DN(1,0)*clhs111 + DN(1,1)*clhs113 + DN(1,2)*clhs116 + clhs246;
            lhs(4,14)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs247;
            lhs(4,15)=DN(3,0)*clhs233 + clhs248*clhs69 - clhs248 + clhs249*clhs71;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs126 + DN(1,2)*clhs127 + clhs61;
            lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs128 + DN(1,2)*clhs130 + clhs143 + clhs225;
            lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs186;
            lhs(5,3)=clhs148*clhs71 + clhs149*clhs69 - clhs149 + clhs214*clhs52;
            lhs(5,4)=DN(1,0)*clhs46 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs230;
            lhs(5,5)=DN(1,0)*clhs57 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs11*clhs250 + clhs228;
            lhs(5,6)=DN(1,0)*clhs64 + DN(1,1)*clhs144 + DN(1,2)*clhs146 + clhs252;
            lhs(5,7)=DN(1,1)*clhs234;
            lhs(5,8)=DN(1,0)*clhs75 + DN(1,1)*clhs151 + DN(1,2)*clhs152 + clhs253;
            lhs(5,9)=DN(1,0)*clhs86 + DN(1,1)*clhs154 + DN(1,2)*clhs156 + clhs238 + clhs254;
            lhs(5,10)=DN(1,0)*clhs93 + DN(1,1)*clhs158 + DN(1,2)*clhs160 + clhs255;
            lhs(5,11)=clhs164*clhs52 + clhs256*clhs69 - clhs256 + clhs257*clhs71;
            lhs(5,12)=DN(1,0)*clhs102 + DN(1,1)*clhs165 + DN(1,2)*clhs166 + clhs258;
            lhs(5,13)=DN(1,0)*clhs113 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs245 + clhs259;
            lhs(5,14)=DN(1,0)*clhs120 + DN(1,1)*clhs172 + DN(1,2)*clhs174 + clhs260;
            lhs(5,15)=DN(3,1)*clhs233 + clhs261*clhs69 - clhs261 + clhs262*clhs71;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs178 + clhs67;
            lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs130 + DN(1,2)*clhs179 + clhs147;
            lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs134 + DN(1,2)*clhs180 + clhs188 + clhs225;
            lhs(6,3)=clhs189*clhs71 + clhs190*clhs69 - clhs190 + clhs215*clhs52;
            lhs(6,4)=DN(1,0)*clhs48 + DN(1,1)*clhs138 + DN(1,2)*clhs182 + clhs231;
            lhs(6,5)=DN(1,0)*clhs60 + DN(1,1)*clhs142 + DN(1,2)*clhs185 + clhs252;
            lhs(6,6)=DN(1,0)*clhs66 + DN(1,1)*clhs146 + DN(1,2)*clhs187 + clhs11*clhs263 + clhs228;
            lhs(6,7)=DN(1,2)*clhs234;
            lhs(6,8)=DN(1,0)*clhs77 + DN(1,1)*clhs152 + DN(1,2)*clhs192 + clhs265;
            lhs(6,9)=DN(1,0)*clhs89 + DN(1,1)*clhs156 + DN(1,2)*clhs194 + clhs266;
            lhs(6,10)=DN(1,0)*clhs95 + DN(1,1)*clhs160 + DN(1,2)*clhs196 + clhs238 + clhs267;
            lhs(6,11)=clhs200*clhs52 + clhs268*clhs69 - clhs268 + clhs269*clhs71;
            lhs(6,12)=DN(1,0)*clhs104 + DN(1,1)*clhs166 + DN(1,2)*clhs201 + clhs270;
            lhs(6,13)=DN(1,0)*clhs116 + DN(1,1)*clhs170 + DN(1,2)*clhs203 + clhs271;
            lhs(6,14)=DN(1,0)*clhs122 + DN(1,1)*clhs174 + DN(1,2)*clhs205 + clhs245 + clhs272;
            lhs(6,15)=DN(3,2)*clhs233 + clhs273*clhs69 - clhs273 + clhs274*clhs71;
            lhs(7,0)=clhs16*clhs72 + clhs68;
            lhs(7,1)=clhs148 + clhs150*clhs16;
            lhs(7,2)=clhs16*clhs191 + clhs189;
            lhs(7,3)=clhs219;
            lhs(7,4)=DN(1,0)*clhs276;
            lhs(7,5)=DN(1,1)*clhs276;
            lhs(7,6)=DN(1,2)*clhs276;
            lhs(7,7)=clhs211*clhs227 + clhs212*clhs226 + clhs212*clhs250 + clhs212*clhs263;
            lhs(7,8)=clhs242 + clhs72*clhs82;
            lhs(7,9)=clhs150*clhs82 + clhs257;
            lhs(7,10)=clhs191*clhs82 + clhs269;
            lhs(7,11)=clhs280;
            lhs(7,12)=clhs109*clhs72 + clhs249;
            lhs(7,13)=clhs109*clhs150 + clhs262;
            lhs(7,14)=clhs109*clhs191 + clhs274;
            lhs(7,15)=clhs281;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs285 + clhs78;
            lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs153;
            lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs193;
            lhs(8,3)=clhs213*clhs80 + clhs69*clhs98 + clhs71*clhs97 - clhs98;
            lhs(8,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + DN(2,2)*clhs48 + clhs235 + clhs286;
            lhs(8,5)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs253;
            lhs(8,6)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs265;
            lhs(8,7)=clhs241*clhs71 + clhs242*clhs69 - clhs242 + clhs72*clhs80;
            lhs(8,8)=DN(2,0)*clhs73 + DN(2,1)*clhs75 + DN(2,2)*clhs77 + clhs11*clhs287 + clhs289;
            lhs(8,9)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs89 + clhs291;
            lhs(8,10)=DN(2,0)*clhs91 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs292;
            lhs(8,11)=DN(2,0)*clhs295;
            lhs(8,12)=DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs296 + clhs299;
            lhs(8,13)=DN(2,0)*clhs111 + DN(2,1)*clhs113 + DN(2,2)*clhs116 + clhs300;
            lhs(8,14)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs301;
            lhs(8,15)=DN(3,0)*clhs294 + clhs302*clhs69 - clhs302 + clhs303*clhs71;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs126 + DN(2,2)*clhs127 + clhs90;
            lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs128 + DN(2,2)*clhs130 + clhs157 + clhs285;
            lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs195;
            lhs(9,3)=clhs162*clhs71 + clhs163*clhs69 - clhs163 + clhs214*clhs80;
            lhs(9,4)=DN(2,0)*clhs46 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs239;
            lhs(9,5)=DN(2,0)*clhs57 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs254 + clhs286;
            lhs(9,6)=DN(2,0)*clhs64 + DN(2,1)*clhs144 + DN(2,2)*clhs146 + clhs266;
            lhs(9,7)=clhs150*clhs80 + clhs256*clhs71 + clhs257*clhs69 - clhs257;
            lhs(9,8)=DN(2,0)*clhs75 + DN(2,1)*clhs151 + DN(2,2)*clhs152 + clhs291;
            lhs(9,9)=DN(2,0)*clhs86 + DN(2,1)*clhs154 + DN(2,2)*clhs156 + clhs11*clhs304 + clhs289;
            lhs(9,10)=DN(2,0)*clhs93 + DN(2,1)*clhs158 + DN(2,2)*clhs160 + clhs306;
            lhs(9,11)=DN(2,1)*clhs295;
            lhs(9,12)=DN(2,0)*clhs102 + DN(2,1)*clhs165 + DN(2,2)*clhs166 + clhs307;
            lhs(9,13)=DN(2,0)*clhs113 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs299 + clhs308;
            lhs(9,14)=DN(2,0)*clhs120 + DN(2,1)*clhs172 + DN(2,2)*clhs174 + clhs309;
            lhs(9,15)=DN(3,1)*clhs294 + clhs310*clhs69 - clhs310 + clhs311*clhs71;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs178 + clhs96;
            lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs130 + DN(2,2)*clhs179 + clhs161;
            lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs134 + DN(2,2)*clhs180 + clhs197 + clhs285;
            lhs(10,3)=clhs198*clhs71 + clhs199*clhs69 - clhs199 + clhs215*clhs80;
            lhs(10,4)=DN(2,0)*clhs48 + DN(2,1)*clhs138 + DN(2,2)*clhs182 + clhs240;
            lhs(10,5)=DN(2,0)*clhs60 + DN(2,1)*clhs142 + DN(2,2)*clhs185 + clhs255;
            lhs(10,6)=DN(2,0)*clhs66 + DN(2,1)*clhs146 + DN(2,2)*clhs187 + clhs267 + clhs286;
            lhs(10,7)=clhs191*clhs80 + clhs268*clhs71 + clhs269*clhs69 - clhs269;
            lhs(10,8)=DN(2,0)*clhs77 + DN(2,1)*clhs152 + DN(2,2)*clhs192 + clhs292;
            lhs(10,9)=DN(2,0)*clhs89 + DN(2,1)*clhs156 + DN(2,2)*clhs194 + clhs306;
            lhs(10,10)=DN(2,0)*clhs95 + DN(2,1)*clhs160 + DN(2,2)*clhs196 + clhs11*clhs312 + clhs289;
            lhs(10,11)=DN(2,2)*clhs295;
            lhs(10,12)=DN(2,0)*clhs104 + DN(2,1)*clhs166 + DN(2,2)*clhs201 + clhs314;
            lhs(10,13)=DN(2,0)*clhs116 + DN(2,1)*clhs170 + DN(2,2)*clhs203 + clhs315;
            lhs(10,14)=DN(2,0)*clhs122 + DN(2,1)*clhs174 + DN(2,2)*clhs205 + clhs299 + clhs316;
            lhs(10,15)=DN(3,2)*clhs294 + clhs317*clhs69 - clhs317 + clhs318*clhs71;
            lhs(11,0)=clhs16*clhs99 + clhs97;
            lhs(11,1)=clhs16*clhs164 + clhs162;
            lhs(11,2)=clhs16*clhs200 + clhs198;
            lhs(11,3)=clhs220;
            lhs(11,4)=clhs241 + clhs53*clhs99;
            lhs(11,5)=clhs164*clhs53 + clhs256;
            lhs(11,6)=clhs200*clhs53 + clhs268;
            lhs(11,7)=clhs280;
            lhs(11,8)=DN(2,0)*clhs320;
            lhs(11,9)=DN(2,1)*clhs320;
            lhs(11,10)=DN(2,2)*clhs320;
            lhs(11,11)=clhs211*clhs288 + clhs212*clhs287 + clhs212*clhs304 + clhs212*clhs312;
            lhs(11,12)=clhs109*clhs99 + clhs303;
            lhs(11,13)=clhs109*clhs164 + clhs311;
            lhs(11,14)=clhs109*clhs200 + clhs318;
            lhs(11,15)=clhs321;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs105 + clhs325;
            lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs167;
            lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs202;
            lhs(12,3)=clhs107*clhs213 + clhs124*clhs71 + clhs125*clhs69 - clhs125;
            lhs(12,4)=DN(3,0)*clhs44 + DN(3,1)*clhs46 + DN(3,2)*clhs48 + clhs243 + clhs326;
            lhs(12,5)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs258;
            lhs(12,6)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs270;
            lhs(12,7)=clhs107*clhs72 + clhs248*clhs71 + clhs249*clhs69 - clhs249;
            lhs(12,8)=DN(3,0)*clhs73 + DN(3,1)*clhs75 + DN(3,2)*clhs77 + clhs296 + clhs327;
            lhs(12,9)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs89 + clhs307;
            lhs(12,10)=DN(3,0)*clhs91 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs314;
            lhs(12,11)=clhs107*clhs99 + clhs302*clhs71 + clhs303*clhs69 - clhs303;
            lhs(12,12)=DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs11*clhs328 + clhs330;
            lhs(12,13)=DN(3,0)*clhs111 + DN(3,1)*clhs113 + DN(3,2)*clhs116 + clhs332;
            lhs(12,14)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs333;
            lhs(12,15)=DN(3,0)*clhs334;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs126 + DN(3,2)*clhs127 + clhs117;
            lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs128 + DN(3,2)*clhs130 + clhs171 + clhs325;
            lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs204;
            lhs(13,3)=clhs107*clhs214 + clhs176*clhs71 + clhs177*clhs69 - clhs177;
            lhs(13,4)=DN(3,0)*clhs46 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs246;
            lhs(13,5)=DN(3,0)*clhs57 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs259 + clhs326;
            lhs(13,6)=DN(3,0)*clhs64 + DN(3,1)*clhs144 + DN(3,2)*clhs146 + clhs271;
            lhs(13,7)=clhs107*clhs150 + clhs261*clhs71 + clhs262*clhs69 - clhs262;
            lhs(13,8)=DN(3,0)*clhs75 + DN(3,1)*clhs151 + DN(3,2)*clhs152 + clhs300;
            lhs(13,9)=DN(3,0)*clhs86 + DN(3,1)*clhs154 + DN(3,2)*clhs156 + clhs308 + clhs327;
            lhs(13,10)=DN(3,0)*clhs93 + DN(3,1)*clhs158 + DN(3,2)*clhs160 + clhs315;
            lhs(13,11)=clhs107*clhs164 + clhs310*clhs71 + clhs311*clhs69 - clhs311;
            lhs(13,12)=DN(3,0)*clhs102 + DN(3,1)*clhs165 + DN(3,2)*clhs166 + clhs332;
            lhs(13,13)=DN(3,0)*clhs113 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs11*clhs335 + clhs330;
            lhs(13,14)=DN(3,0)*clhs120 + DN(3,1)*clhs172 + DN(3,2)*clhs174 + clhs336;
            lhs(13,15)=DN(3,1)*clhs334;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs178 + clhs123;
            lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs130 + DN(3,2)*clhs179 + clhs175;
            lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs134 + DN(3,2)*clhs180 + clhs206 + clhs325;
            lhs(14,3)=clhs107*clhs215 + clhs207*clhs71 + clhs208*clhs69 - clhs208;
            lhs(14,4)=DN(3,0)*clhs48 + DN(3,1)*clhs138 + DN(3,2)*clhs182 + clhs247;
            lhs(14,5)=DN(3,0)*clhs60 + DN(3,1)*clhs142 + DN(3,2)*clhs185 + clhs260;
            lhs(14,6)=DN(3,0)*clhs66 + DN(3,1)*clhs146 + DN(3,2)*clhs187 + clhs272 + clhs326;
            lhs(14,7)=clhs107*clhs191 + clhs273*clhs71 + clhs274*clhs69 - clhs274;
            lhs(14,8)=DN(3,0)*clhs77 + DN(3,1)*clhs152 + DN(3,2)*clhs192 + clhs301;
            lhs(14,9)=DN(3,0)*clhs89 + DN(3,1)*clhs156 + DN(3,2)*clhs194 + clhs309;
            lhs(14,10)=DN(3,0)*clhs95 + DN(3,1)*clhs160 + DN(3,2)*clhs196 + clhs316 + clhs327;
            lhs(14,11)=clhs107*clhs200 + clhs317*clhs71 + clhs318*clhs69 - clhs318;
            lhs(14,12)=DN(3,0)*clhs104 + DN(3,1)*clhs166 + DN(3,2)*clhs201 + clhs333;
            lhs(14,13)=DN(3,0)*clhs116 + DN(3,1)*clhs170 + DN(3,2)*clhs203 + clhs336;
            lhs(14,14)=DN(3,0)*clhs122 + DN(3,1)*clhs174 + DN(3,2)*clhs205 + clhs11*clhs337 + clhs330;
            lhs(14,15)=DN(3,2)*clhs334;
            lhs(15,0)=DN(3,0)*clhs209 + clhs124;
            lhs(15,1)=DN(3,1)*clhs209 + clhs176;
            lhs(15,2)=DN(3,2)*clhs209 + clhs207;
            lhs(15,3)=clhs221;
            lhs(15,4)=DN(3,0)*clhs275 + clhs248;
            lhs(15,5)=DN(3,1)*clhs275 + clhs261;
            lhs(15,6)=DN(3,2)*clhs275 + clhs273;
            lhs(15,7)=clhs281;
            lhs(15,8)=DN(3,0)*clhs319 + clhs302;
            lhs(15,9)=DN(3,1)*clhs319 + clhs310;
            lhs(15,10)=DN(3,2)*clhs319 + clhs317;
            lhs(15,11)=clhs321;
            lhs(15,12)=DN(3,0)*clhs338;
            lhs(15,13)=DN(3,1)*clhs338;
            lhs(15,14)=DN(3,2)*clhs338;
            lhs(15,15)=clhs211*clhs329 + clhs212*clhs328 + clhs212*clhs335 + clhs212*clhs337;


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

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crhs2 =             crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs3 =             crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs7 =             crhs1*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs8 =             crhs1*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
const double crhs9 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs10 =             crhs4 + crhs9;
const double crhs11 =             1.0/crhs1;
const double crhs12 =             crhs11*(crhs5*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crhs6*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crhs13 =             crhs11*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crhs14 =             (crhs8*h/stab_c1 + mu)*(crhs10 + crhs12 + crhs13);
const double crhs15 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs16 =             N[0]*crhs1*crhs15;
const double crhs17 =             1.0/(crhs1*dyn_tau/dt + crhs8/h + mu*stab_c1/pow(h, 2));
const double crhs18 =             1.0*crhs17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs2 + crhs3 + crhs7);
const double crhs19 =             crhs1*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
const double crhs20 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs21 =             crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs22 =             crhs1*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs9);
const double crhs23 =             1.0*crhs17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs20 + crhs21 + crhs22);
const double crhs24 =             N[1]*crhs1*crhs15;
const double crhs25 =             crhs1*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
const double crhs26 =             N[2]*crhs1*crhs15;
const double crhs27 =             crhs1*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs7 - crhs16*crhs18 - crhs18*crhs19;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs14 - DN(0,1)*stress[1] + N[0]*crhs20 - N[0]*crhs21 - N[0]*crhs22 - crhs16*crhs23 - crhs19*crhs23;
            rhs[2]=-DN(0,0)*crhs18 - DN(0,1)*crhs23 - N[0]*crhs10 - N[0]*crhs12 - N[0]*crhs13;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs7 - crhs18*crhs24 - crhs18*crhs25;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs14 - DN(1,1)*stress[1] + N[1]*crhs20 - N[1]*crhs21 - N[1]*crhs22 - crhs23*crhs24 - crhs23*crhs25;
            rhs[5]=-DN(1,0)*crhs18 - DN(1,1)*crhs23 - N[1]*crhs10 - N[1]*crhs12 - N[1]*crhs13;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs7 - crhs18*crhs26 - crhs18*crhs27;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs14 - DN(2,1)*stress[1] + N[2]*crhs20 - N[2]*crhs21 - N[2]*crhs22 - crhs23*crhs26 - crhs23*crhs27;
            rhs[8]=-DN(2,0)*crhs18 - DN(2,1)*crhs23 - N[2]*crhs10 - N[2]*crhs12 - N[2]*crhs13;


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

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crhs2 =             crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs3 =             crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 =             DN(0,0)*v(0,0);
const double crhs5 =             DN(1,0)*v(1,0);
const double crhs6 =             DN(2,0)*v(2,0);
const double crhs7 =             DN(3,0)*v(3,0);
const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs11 =             crhs1*(crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*(crhs4 + crhs5 + crhs6 + crhs7) + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
const double crhs12 =             crhs1*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs13 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs14 =             DN(0,1)*v(0,1);
const double crhs15 =             DN(1,1)*v(1,1);
const double crhs16 =             DN(2,1)*v(2,1);
const double crhs17 =             DN(3,1)*v(3,1);
const double crhs18 =             crhs13 + crhs14 + crhs15 + crhs16 + crhs17 + crhs4 + crhs5 + crhs6 + crhs7;
const double crhs19 =             1.0/crhs1;
const double crhs20 =             crhs19*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crhs21 =             crhs19*(crhs10*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]) + crhs8*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crhs9*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]));
const double crhs22 =             (crhs12*h/stab_c1 + mu)*(crhs18 + crhs20 + crhs21);
const double crhs23 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs24 =             N[0]*crhs1*crhs23;
const double crhs25 =             1.0/(crhs1*dyn_tau/dt + crhs12/h + mu*stab_c1/pow(h, 2));
const double crhs26 =             1.0*crhs25*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs11 - crhs2 + crhs3);
const double crhs27 =             crhs1*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10);
const double crhs28 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs29 =             crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs30 =             crhs1*(crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs9*(crhs14 + crhs15 + crhs16 + crhs17));
const double crhs31 =             1.0*crhs25*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs28 + crhs29 + crhs30);
const double crhs32 =             crhs1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs33 =             crhs1*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs34 =             crhs1*(crhs10*crhs13 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs35 =             1.0*crhs25*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs32 + crhs33 + crhs34);
const double crhs36 =             N[1]*crhs1*crhs23;
const double crhs37 =             crhs1*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10);
const double crhs38 =             N[2]*crhs1*crhs23;
const double crhs39 =             crhs1*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10);
const double crhs40 =             N[3]*crhs1*crhs23;
const double crhs41 =             crhs1*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs22 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs11 + N[0]*crhs2 - N[0]*crhs3 - crhs24*crhs26 - crhs26*crhs27;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs22 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs24*crhs31 - crhs27*crhs31;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs22 - DN(0,2)*stress[2] + N[0]*crhs32 - N[0]*crhs33 - N[0]*crhs34 - crhs24*crhs35 - crhs27*crhs35;
            rhs[3]=-DN(0,0)*crhs26 - DN(0,1)*crhs31 - DN(0,2)*crhs35 - N[0]*crhs18 - N[0]*crhs20 - N[0]*crhs21;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs22 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs11 + N[1]*crhs2 - N[1]*crhs3 - crhs26*crhs36 - crhs26*crhs37;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs22 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs31*crhs36 - crhs31*crhs37;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs22 - DN(1,2)*stress[2] + N[1]*crhs32 - N[1]*crhs33 - N[1]*crhs34 - crhs35*crhs36 - crhs35*crhs37;
            rhs[7]=-DN(1,0)*crhs26 - DN(1,1)*crhs31 - DN(1,2)*crhs35 - N[1]*crhs18 - N[1]*crhs20 - N[1]*crhs21;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs22 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs11 + N[2]*crhs2 - N[2]*crhs3 - crhs26*crhs38 - crhs26*crhs39;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs22 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs31*crhs38 - crhs31*crhs39;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs22 - DN(2,2)*stress[2] + N[2]*crhs32 - N[2]*crhs33 - N[2]*crhs34 - crhs35*crhs38 - crhs35*crhs39;
            rhs[11]=-DN(2,0)*crhs26 - DN(2,1)*crhs31 - DN(2,2)*crhs35 - N[2]*crhs18 - N[2]*crhs20 - N[2]*crhs21;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs22 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs11 + N[3]*crhs2 - N[3]*crhs3 - crhs26*crhs40 - crhs26*crhs41;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs22 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs28 - N[3]*crhs29 - N[3]*crhs30 - crhs31*crhs40 - crhs31*crhs41;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs22 - DN(3,2)*stress[2] + N[3]*crhs32 - N[3]*crhs33 - N[3]*crhs34 - crhs35*crhs40 - crhs35*crhs41;
            rhs[15]=-DN(3,0)*crhs26 - DN(3,1)*crhs31 - DN(3,2)*crhs35 - N[3]*crhs18 - N[3]*crhs20 - N[3]*crhs21;


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