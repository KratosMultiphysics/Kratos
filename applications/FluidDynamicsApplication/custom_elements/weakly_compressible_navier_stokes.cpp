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

    // Non-inertial frame or reference data
    const array_1d<double,3>& omega = rData.AngularVelocity;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

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
const double clhs13 =             pow(omega[2], 2);
const double clhs14 =             pow(clhs4, 2);
const double clhs15 =             1.0/(clhs4*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2)) + clhs4*dyn_tau/dt + clhs7/h + mu*stab_c1/pow(h, 2));
const double clhs16 =             4.0*clhs13*clhs14*clhs15;
const double clhs17 =             N[0]*bdf0 + clhs12;
const double clhs18 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs19 =             1.0*N[0]*clhs14*clhs15*clhs18;
const double clhs20 =             1.0*clhs12*clhs14*clhs15;
const double clhs21 =             clhs10*clhs9 + clhs11*clhs12 + clhs16*clhs9 + clhs17*clhs19 + clhs17*clhs20;
const double clhs22 =             C(0,1)*DN(0,1) + clhs1;
const double clhs23 =             C(1,2)*DN(0,1);
const double clhs24 =             C(2,2)*DN(0,0) + clhs23;
const double clhs25 =             2.0*clhs4*omega[2];
const double clhs26 =             clhs25*clhs9;
const double clhs27 =             DN(0,0)*clhs8;
const double clhs28 =             DN(0,1)*clhs27;
const double clhs29 =             2.0*clhs14*clhs15*clhs18*omega[2];
const double clhs30 =             clhs29*clhs9;
const double clhs31 =             2.0*N[0]*clhs14*clhs15*omega[2];
const double clhs32 =             clhs12*clhs31;
const double clhs33 =             clhs17*clhs31;
const double clhs34 =             DN(0,0)*N[0];
const double clhs35 =             pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], -2);
const double clhs36 =             1.0/clhs4;
const double clhs37 =             bdf0*clhs35*clhs36*clhs8;
const double clhs38 =             DN(0,1)*N[0];
const double clhs39 =             2.0*clhs15*clhs4*omega[2];
const double clhs40 =             clhs38*clhs39;
const double clhs41 =             1.0*clhs15*clhs18*clhs4;
const double clhs42 =             1.0*DN(0,0)*clhs15*clhs4;
const double clhs43 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs44 =             C(0,2)*DN(1,0);
const double clhs45 =             C(2,2)*DN(1,1) + clhs44;
const double clhs46 =             DN(1,0)*clhs27;
const double clhs47 =             N[0]*bdf0*clhs4;
const double clhs48 =             N[1]*clhs47;
const double clhs49 =             DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs50 =             4.0*N[0]*clhs13*clhs14*clhs15;
const double clhs51 =             N[1]*clhs50;
const double clhs52 =             N[1]*bdf0;
const double clhs53 =             clhs49 + clhs52;
const double clhs54 =             clhs11*clhs49 + clhs19*clhs53 + clhs20*clhs53 + clhs48 + clhs51;
const double clhs55 =             2.0*N[0]*omega[2];
const double clhs56 =             N[1]*clhs4;
const double clhs57 =             clhs55*clhs56;
const double clhs58 =             2.0*N[0]*clhs14*clhs15*clhs18*omega[2];
const double clhs59 =             N[1]*clhs58;
const double clhs60 =             -clhs57 - clhs59;
const double clhs61 =             C(0,1)*DN(1,1) + clhs44;
const double clhs62 =             C(1,2)*DN(1,1);
const double clhs63 =             C(2,2)*DN(1,0) + clhs62;
const double clhs64 =             DN(1,1)*clhs27;
const double clhs65 =             2.0*N[1]*clhs14*clhs15*omega[2];
const double clhs66 =             clhs12*clhs65;
const double clhs67 =             clhs31*clhs53;
const double clhs68 =             DN(0,0)*N[1];
const double clhs69 =             DN(1,1)*N[0];
const double clhs70 =             clhs39*clhs69;
const double clhs71 =             DN(1,0)*N[0];
const double clhs72 =             1.0*DN(1,0)*clhs15*clhs4;
const double clhs73 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs74 =             C(0,2)*DN(2,0);
const double clhs75 =             C(2,2)*DN(2,1) + clhs74;
const double clhs76 =             DN(2,0)*clhs27;
const double clhs77 =             N[2]*clhs47;
const double clhs78 =             DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs79 =             N[2]*clhs50;
const double clhs80 =             N[2]*bdf0 + clhs78;
const double clhs81 =             clhs11*clhs78 + clhs19*clhs80 + clhs20*clhs80 + clhs77 + clhs79;
const double clhs82 =             N[2]*clhs4;
const double clhs83 =             clhs55*clhs82;
const double clhs84 =             N[2]*clhs58;
const double clhs85 =             -clhs83 - clhs84;
const double clhs86 =             C(0,1)*DN(2,1) + clhs74;
const double clhs87 =             C(1,2)*DN(2,1);
const double clhs88 =             C(2,2)*DN(2,0) + clhs87;
const double clhs89 =             DN(2,1)*clhs27;
const double clhs90 =             2.0*N[2]*clhs14*clhs15*omega[2];
const double clhs91 =             clhs12*clhs90;
const double clhs92 =             clhs31*clhs80;
const double clhs93 =             DN(0,0)*N[2];
const double clhs94 =             DN(2,1)*N[0];
const double clhs95 =             clhs39*clhs94;
const double clhs96 =             DN(2,0)*N[0];
const double clhs97 =             1.0*DN(2,0)*clhs15*clhs4;
const double clhs98 =             C(0,1)*DN(0,0) + clhs23;
const double clhs99 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs100 =             pow(DN(0,1), 2);
const double clhs101 =             -clhs34*clhs39;
const double clhs102 =             1.0*DN(0,1)*clhs15*clhs4;
const double clhs103 =             clhs57 + clhs59;
const double clhs104 =             C(0,1)*DN(1,0) + clhs62;
const double clhs105 =             DN(0,1)*clhs8;
const double clhs106 =             DN(1,0)*clhs105;
const double clhs107 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs108 =             DN(1,1)*clhs105;
const double clhs109 =             DN(0,1)*N[1];
const double clhs110 =             -clhs39*clhs71;
const double clhs111 =             1.0*DN(1,1)*clhs15*clhs4;
const double clhs112 =             clhs83 + clhs84;
const double clhs113 =             C(0,1)*DN(2,0) + clhs87;
const double clhs114 =             DN(2,0)*clhs105;
const double clhs115 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs116 =             DN(2,1)*clhs105;
const double clhs117 =             DN(0,1)*N[2];
const double clhs118 =             -clhs39*clhs96;
const double clhs119 =             1.0*DN(2,1)*clhs15*clhs4;
const double clhs120 =             bdf0*clhs35*clhs36;
const double clhs121 =             1.0*clhs15;
const double clhs122 =             clhs109*clhs39;
const double clhs123 =             -clhs39*clhs68;
const double clhs124 =             N[0]*bdf0*clhs35*clhs36;
const double clhs125 =             1.0*DN(0,0)*clhs15;
const double clhs126 =             1.0*DN(0,1)*clhs15;
const double clhs127 =             DN(1,0)*clhs125 + DN(1,1)*clhs126 + N[1]*clhs124;
const double clhs128 =             clhs117*clhs39;
const double clhs129 =             -clhs39*clhs93;
const double clhs130 =             DN(2,0)*clhs125 + DN(2,1)*clhs126 + N[2]*clhs124;
const double clhs131 =             1.0*N[1]*clhs14*clhs15*clhs18;
const double clhs132 =             1.0*clhs14*clhs15*clhs49;
const double clhs133 =             clhs12*clhs56 + clhs131*clhs17 + clhs132*clhs17 + clhs48 + clhs51;
const double clhs134 =             clhs31*clhs49;
const double clhs135 =             clhs17*clhs65;
const double clhs136 =             pow(DN(1,0), 2);
const double clhs137 =             pow(N[1], 2);
const double clhs138 =             clhs10*clhs137 + clhs131*clhs53 + clhs132*clhs53 + clhs137*clhs16 + clhs49*clhs56;
const double clhs139 =             clhs137*clhs25;
const double clhs140 =             DN(1,0)*clhs8;
const double clhs141 =             DN(1,1)*clhs140;
const double clhs142 =             clhs137*clhs29;
const double clhs143 =             clhs49*clhs65;
const double clhs144 =             clhs53*clhs65;
const double clhs145 =             DN(1,0)*N[1];
const double clhs146 =             DN(1,1)*N[1];
const double clhs147 =             clhs146*clhs39;
const double clhs148 =             DN(2,0)*clhs140;
const double clhs149 =             clhs52*clhs82;
const double clhs150 =             N[1]*N[2]*clhs16;
const double clhs151 =             clhs131*clhs80 + clhs132*clhs80 + clhs149 + clhs150 + clhs56*clhs78;
const double clhs152 =             2.0*N[1]*clhs82*omega[2];
const double clhs153 =             2.0*N[1]*N[2]*clhs14*clhs15*clhs18*omega[2];
const double clhs154 =             -clhs152 - clhs153;
const double clhs155 =             DN(2,1)*clhs140;
const double clhs156 =             clhs49*clhs90;
const double clhs157 =             clhs65*clhs80;
const double clhs158 =             DN(1,0)*N[2];
const double clhs159 =             DN(2,1)*N[1];
const double clhs160 =             clhs159*clhs39;
const double clhs161 =             DN(2,0)*N[1];
const double clhs162 =             pow(DN(1,1), 2);
const double clhs163 =             -clhs145*clhs39;
const double clhs164 =             clhs152 + clhs153;
const double clhs165 =             DN(1,1)*clhs8;
const double clhs166 =             DN(2,0)*clhs165;
const double clhs167 =             DN(2,1)*clhs165;
const double clhs168 =             DN(1,1)*N[2];
const double clhs169 =             -clhs161*clhs39;
const double clhs170 =             clhs168*clhs39;
const double clhs171 =             -clhs158*clhs39;
const double clhs172 =             1.0*DN(1,0)*DN(2,0)*clhs15 + 1.0*DN(1,1)*DN(2,1)*clhs15 + N[1]*N[2]*bdf0*clhs35*clhs36;
const double clhs173 =             1.0*N[2]*clhs14*clhs15*clhs18;
const double clhs174 =             1.0*clhs14*clhs15*clhs78;
const double clhs175 =             clhs12*clhs82 + clhs17*clhs173 + clhs17*clhs174 + clhs77 + clhs79;
const double clhs176 =             clhs31*clhs78;
const double clhs177 =             clhs17*clhs90;
const double clhs178 =             clhs149 + clhs150 + clhs173*clhs53 + clhs174*clhs53 + clhs49*clhs82;
const double clhs179 =             clhs65*clhs78;
const double clhs180 =             clhs53*clhs90;
const double clhs181 =             pow(DN(2,0), 2);
const double clhs182 =             pow(N[2], 2);
const double clhs183 =             clhs10*clhs182 + clhs16*clhs182 + clhs173*clhs80 + clhs174*clhs80 + clhs78*clhs82;
const double clhs184 =             clhs182*clhs25;
const double clhs185 =             DN(2,0)*DN(2,1)*clhs8;
const double clhs186 =             clhs182*clhs29;
const double clhs187 =             clhs78*clhs90;
const double clhs188 =             clhs80*clhs90;
const double clhs189 =             DN(2,0)*N[2];
const double clhs190 =             DN(2,1)*N[2];
const double clhs191 =             clhs190*clhs39;
const double clhs192 =             pow(DN(2,1), 2);
const double clhs193 =             -clhs189*clhs39;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs21 + clhs3*clhs8;
            lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 - clhs26 + clhs28 - clhs30 - clhs32 + clhs33;
            lhs(0,2)=clhs12*clhs42 + clhs34*clhs37 + clhs34*clhs41 - clhs34 + clhs40;
            lhs(0,3)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46 + clhs54;
            lhs(0,4)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + clhs60 + clhs64 - clhs66 + clhs67;
            lhs(0,5)=clhs12*clhs72 + clhs37*clhs68 + clhs41*clhs71 - clhs68 + clhs70;
            lhs(0,6)=DN(0,0)*clhs73 + DN(0,1)*clhs75 + clhs76 + clhs81;
            lhs(0,7)=DN(0,0)*clhs86 + DN(0,1)*clhs88 + clhs85 + clhs89 - clhs91 + clhs92;
            lhs(0,8)=clhs12*clhs97 + clhs37*clhs93 + clhs41*clhs96 - clhs93 + clhs95;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs98 + clhs26 + clhs28 + clhs30 + clhs32 - clhs33;
            lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs99 + clhs100*clhs8 + clhs21;
            lhs(1,2)=clhs101 + clhs102*clhs12 + clhs37*clhs38 + clhs38*clhs41 - clhs38;
            lhs(1,3)=DN(0,0)*clhs45 + DN(0,1)*clhs104 + clhs103 + clhs106 + clhs66 - clhs67;
            lhs(1,4)=DN(0,0)*clhs63 + DN(0,1)*clhs107 + clhs108 + clhs54;
            lhs(1,5)=clhs109*clhs37 - clhs109 + clhs110 + clhs111*clhs12 + clhs41*clhs69;
            lhs(1,6)=DN(0,0)*clhs75 + DN(0,1)*clhs113 + clhs112 + clhs114 + clhs91 - clhs92;
            lhs(1,7)=DN(0,0)*clhs88 + DN(0,1)*clhs115 + clhs116 + clhs81;
            lhs(1,8)=clhs117*clhs37 - clhs117 + clhs118 + clhs119*clhs12 + clhs41*clhs94;
            lhs(2,0)=clhs17*clhs42 + clhs34 + clhs40;
            lhs(2,1)=clhs101 + clhs102*clhs17 + clhs38;
            lhs(2,2)=clhs100*clhs121 + clhs120*clhs9 + clhs121*clhs3;
            lhs(2,3)=clhs122 + clhs42*clhs53 + clhs71;
            lhs(2,4)=clhs102*clhs53 + clhs123 + clhs69;
            lhs(2,5)=clhs127;
            lhs(2,6)=clhs128 + clhs42*clhs80 + clhs96;
            lhs(2,7)=clhs102*clhs80 + clhs129 + clhs94;
            lhs(2,8)=clhs130;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs133 + clhs46;
            lhs(3,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + clhs106 - clhs134 + clhs135 + clhs60;
            lhs(3,2)=clhs122 + clhs37*clhs71 + clhs41*clhs68 + clhs42*clhs49 - clhs71;
            lhs(3,3)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs136*clhs8 + clhs138;
            lhs(3,4)=DN(1,0)*clhs61 + DN(1,1)*clhs63 - clhs139 + clhs141 - clhs142 - clhs143 + clhs144;
            lhs(3,5)=clhs145*clhs37 + clhs145*clhs41 - clhs145 + clhs147 + clhs49*clhs72;
            lhs(3,6)=DN(1,0)*clhs73 + DN(1,1)*clhs75 + clhs148 + clhs151;
            lhs(3,7)=DN(1,0)*clhs86 + DN(1,1)*clhs88 + clhs154 + clhs155 - clhs156 + clhs157;
            lhs(3,8)=clhs158*clhs37 - clhs158 + clhs160 + clhs161*clhs41 + clhs49*clhs97;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs98 + clhs103 + clhs134 - clhs135 + clhs64;
            lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs99 + clhs108 + clhs133;
            lhs(4,2)=clhs102*clhs49 + clhs109*clhs41 + clhs123 + clhs37*clhs69 - clhs69;
            lhs(4,3)=DN(1,0)*clhs45 + DN(1,1)*clhs104 + clhs139 + clhs141 + clhs142 + clhs143 - clhs144;
            lhs(4,4)=DN(1,0)*clhs63 + DN(1,1)*clhs107 + clhs138 + clhs162*clhs8;
            lhs(4,5)=clhs111*clhs49 + clhs146*clhs37 + clhs146*clhs41 - clhs146 + clhs163;
            lhs(4,6)=DN(1,0)*clhs75 + DN(1,1)*clhs113 + clhs156 - clhs157 + clhs164 + clhs166;
            lhs(4,7)=DN(1,0)*clhs88 + DN(1,1)*clhs115 + clhs151 + clhs167;
            lhs(4,8)=clhs119*clhs49 + clhs159*clhs41 + clhs168*clhs37 - clhs168 + clhs169;
            lhs(5,0)=clhs17*clhs72 + clhs68 + clhs70;
            lhs(5,1)=clhs109 + clhs110 + clhs111*clhs17;
            lhs(5,2)=clhs127;
            lhs(5,3)=clhs145 + clhs147 + clhs53*clhs72;
            lhs(5,4)=clhs111*clhs53 + clhs146 + clhs163;
            lhs(5,5)=clhs120*clhs137 + clhs121*clhs136 + clhs121*clhs162;
            lhs(5,6)=clhs161 + clhs170 + clhs72*clhs80;
            lhs(5,7)=clhs111*clhs80 + clhs159 + clhs171;
            lhs(5,8)=clhs172;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs175 + clhs76;
            lhs(6,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + clhs114 - clhs176 + clhs177 + clhs85;
            lhs(6,2)=clhs128 + clhs37*clhs96 + clhs41*clhs93 + clhs42*clhs78 - clhs96;
            lhs(6,3)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs148 + clhs178;
            lhs(6,4)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + clhs154 + clhs166 - clhs179 + clhs180;
            lhs(6,5)=clhs158*clhs41 + clhs161*clhs37 - clhs161 + clhs170 + clhs72*clhs78;
            lhs(6,6)=DN(2,0)*clhs73 + DN(2,1)*clhs75 + clhs181*clhs8 + clhs183;
            lhs(6,7)=DN(2,0)*clhs86 + DN(2,1)*clhs88 - clhs184 + clhs185 - clhs186 - clhs187 + clhs188;
            lhs(6,8)=clhs189*clhs37 + clhs189*clhs41 - clhs189 + clhs191 + clhs78*clhs97;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs98 + clhs112 + clhs176 - clhs177 + clhs89;
            lhs(7,1)=DN(2,0)*clhs24 + DN(2,1)*clhs99 + clhs116 + clhs175;
            lhs(7,2)=clhs102*clhs78 + clhs117*clhs41 + clhs129 + clhs37*clhs94 - clhs94;
            lhs(7,3)=DN(2,0)*clhs45 + DN(2,1)*clhs104 + clhs155 + clhs164 + clhs179 - clhs180;
            lhs(7,4)=DN(2,0)*clhs63 + DN(2,1)*clhs107 + clhs167 + clhs178;
            lhs(7,5)=clhs111*clhs78 + clhs159*clhs37 - clhs159 + clhs168*clhs41 + clhs171;
            lhs(7,6)=DN(2,0)*clhs75 + DN(2,1)*clhs113 + clhs184 + clhs185 + clhs186 + clhs187 - clhs188;
            lhs(7,7)=DN(2,0)*clhs88 + DN(2,1)*clhs115 + clhs183 + clhs192*clhs8;
            lhs(7,8)=clhs119*clhs78 + clhs190*clhs37 + clhs190*clhs41 - clhs190 + clhs193;
            lhs(8,0)=clhs17*clhs97 + clhs93 + clhs95;
            lhs(8,1)=clhs117 + clhs118 + clhs119*clhs17;
            lhs(8,2)=clhs130;
            lhs(8,3)=clhs158 + clhs160 + clhs53*clhs97;
            lhs(8,4)=clhs119*clhs53 + clhs168 + clhs169;
            lhs(8,5)=clhs172;
            lhs(8,6)=clhs189 + clhs191 + clhs80*clhs97;
            lhs(8,7)=clhs119*clhs80 + clhs190 + clhs193;
            lhs(8,8)=clhs120*clhs182 + clhs121*clhs181 + clhs121*clhs192;


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

    // Non-inertial frame or reference data
    const array_1d<double,3>& omega = rData.AngularVelocity;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

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
const double clhs12 =             pow(omega[1], 2);
const double clhs13 =             pow(omega[2], 2);
const double clhs14 =             clhs12 + clhs13;
const double clhs15 =             pow(N[0], 2);
const double clhs16 =             pow(clhs6, 2);
const double clhs17 =             pow(omega[0], 2);
const double clhs18 =             clhs12 + clhs17;
const double clhs19 =             1.0/(clhs10/h + clhs6*stab_c3*sqrt(clhs13 + clhs18) + clhs6*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double clhs20 =             4.0*clhs15*clhs16*clhs19;
const double clhs21 =             bdf0*clhs6;
const double clhs22 =             N[0]*clhs6;
const double clhs23 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs24 =             N[0]*bdf0 + clhs23;
const double clhs25 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs26 =             1.0*N[0]*clhs16*clhs19*clhs25;
const double clhs27 =             1.0*clhs16*clhs19*clhs23;
const double clhs28 =             clhs15*clhs21 + clhs22*clhs23 + clhs24*clhs26 + clhs24*clhs27;
const double clhs29 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs30 =             C(1,3)*DN(0,1);
const double clhs31 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs30;
const double clhs32 =             C(3,5)*DN(0,0);
const double clhs33 =             C(4,5)*DN(0,2);
const double clhs34 =             C(1,5)*DN(0,1) + clhs32 + clhs33;
const double clhs35 =             2.0*omega[2];
const double clhs36 =             clhs15*clhs6;
const double clhs37 =             clhs35*clhs36;
const double clhs38 =             DN(0,0)*clhs11;
const double clhs39 =             DN(0,1)*clhs38;
const double clhs40 =             clhs15*clhs16*clhs19*clhs25;
const double clhs41 =             clhs35*clhs40;
const double clhs42 =             2.0*N[0]*clhs16*clhs19*omega[2];
const double clhs43 =             clhs23*clhs42;
const double clhs44 =             2.0*N[0]*omega[0];
const double clhs45 =             clhs44*omega[1];
const double clhs46 =             1.0*omega[2];
const double clhs47 =             clhs24*clhs46;
const double clhs48 =             clhs45 - clhs47;
const double clhs49 =             2.0*N[0]*clhs16*clhs19;
const double clhs50 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs51 =             C(3,4)*DN(0,1);
const double clhs52 =             C(2,3)*DN(0,2) + clhs32 + clhs51;
const double clhs53 =             C(2,5)*DN(0,2);
const double clhs54 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs53;
const double clhs55 =             2.0*omega[1];
const double clhs56 =             clhs36*clhs55;
const double clhs57 =             DN(0,2)*clhs38;
const double clhs58 =             clhs40*clhs55;
const double clhs59 =             2.0*N[0]*clhs16*clhs19*omega[1];
const double clhs60 =             clhs23*clhs59;
const double clhs61 =             clhs44*omega[2];
const double clhs62 =             1.0*omega[1];
const double clhs63 =             clhs24*clhs62;
const double clhs64 =             clhs61 + clhs63;
const double clhs65 =             DN(0,0)*N[0];
const double clhs66 =             pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], -2);
const double clhs67 =             1.0/clhs6;
const double clhs68 =             bdf0*clhs11*clhs66*clhs67;
const double clhs69 =             DN(0,1)*omega[2] - DN(0,2)*omega[1];
const double clhs70 =             2.0*N[0]*clhs19*clhs6;
const double clhs71 =             1.0*clhs19*clhs25*clhs6;
const double clhs72 =             1.0*DN(0,0)*clhs19*clhs6;
const double clhs73 =             DN(1,0)*clhs38;
const double clhs74 =             4.0*N[0]*N[1]*clhs16*clhs19;
const double clhs75 =             clhs14*clhs74;
const double clhs76 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs77 =             C(0,3)*DN(1,0);
const double clhs78 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs77;
const double clhs79 =             C(0,5)*DN(1,0);
const double clhs80 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs79;
const double clhs81 =             N[0]*bdf0*clhs6;
const double clhs82 =             N[1]*clhs81;
const double clhs83 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs84 =             clhs22*clhs83;
const double clhs85 =             N[1]*bdf0 + clhs83;
const double clhs86 =             clhs26*clhs85;
const double clhs87 =             clhs27*clhs85;
const double clhs88 =             N[0]*N[1]*clhs6;
const double clhs89 =             clhs35*clhs88;
const double clhs90 =             N[0]*N[1]*clhs16*clhs19*clhs25;
const double clhs91 =             clhs35*clhs90;
const double clhs92 =             -clhs89 - clhs91;
const double clhs93 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs77;
const double clhs94 =             C(1,3)*DN(1,1);
const double clhs95 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs94;
const double clhs96 =             C(3,5)*DN(1,0);
const double clhs97 =             C(4,5)*DN(1,2);
const double clhs98 =             C(1,5)*DN(1,1) + clhs96 + clhs97;
const double clhs99 =             DN(1,1)*clhs38;
const double clhs100 =             2.0*N[1]*clhs16*clhs19*omega[2];
const double clhs101 =             clhs100*clhs23;
const double clhs102 =             N[1]*omega[0];
const double clhs103 =             clhs102*clhs55;
const double clhs104 =             clhs46*clhs85;
const double clhs105 =             clhs103 - clhs104;
const double clhs106 =             clhs55*clhs88;
const double clhs107 =             clhs55*clhs90;
const double clhs108 =             clhs106 + clhs107;
const double clhs109 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs79;
const double clhs110 =             C(3,4)*DN(1,1);
const double clhs111 =             C(2,3)*DN(1,2) + clhs110 + clhs96;
const double clhs112 =             C(2,5)*DN(1,2);
const double clhs113 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs112;
const double clhs114 =             DN(1,2)*clhs38;
const double clhs115 =             2.0*N[1]*clhs16*clhs19*omega[1];
const double clhs116 =             clhs115*clhs23;
const double clhs117 =             clhs102*clhs35;
const double clhs118 =             clhs62*clhs85;
const double clhs119 =             clhs117 + clhs118;
const double clhs120 =             DN(0,0)*N[1];
const double clhs121 =             DN(1,1)*omega[2] - DN(1,2)*omega[1];
const double clhs122 =             DN(1,0)*N[0];
const double clhs123 =             1.0*DN(1,0)*clhs19*clhs6;
const double clhs124 =             DN(2,0)*clhs38;
const double clhs125 =             4.0*N[0]*N[2]*clhs16*clhs19;
const double clhs126 =             clhs125*clhs14;
const double clhs127 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs128 =             C(0,3)*DN(2,0);
const double clhs129 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs128;
const double clhs130 =             C(0,5)*DN(2,0);
const double clhs131 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs130;
const double clhs132 =             N[2]*clhs81;
const double clhs133 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs134 =             clhs133*clhs22;
const double clhs135 =             N[2]*bdf0;
const double clhs136 =             clhs133 + clhs135;
const double clhs137 =             clhs136*clhs26;
const double clhs138 =             clhs136*clhs27;
const double clhs139 =             N[0]*N[2]*clhs6;
const double clhs140 =             clhs139*clhs35;
const double clhs141 =             N[0]*N[2]*clhs16*clhs19*clhs25;
const double clhs142 =             clhs141*clhs35;
const double clhs143 =             -clhs140 - clhs142;
const double clhs144 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs128;
const double clhs145 =             C(1,3)*DN(2,1);
const double clhs146 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs145;
const double clhs147 =             C(3,5)*DN(2,0);
const double clhs148 =             C(4,5)*DN(2,2);
const double clhs149 =             C(1,5)*DN(2,1) + clhs147 + clhs148;
const double clhs150 =             DN(2,1)*clhs38;
const double clhs151 =             2.0*N[2]*clhs16*clhs19*omega[2];
const double clhs152 =             clhs151*clhs23;
const double clhs153 =             N[2]*omega[0];
const double clhs154 =             clhs153*clhs55;
const double clhs155 =             clhs136*clhs46;
const double clhs156 =             clhs154 - clhs155;
const double clhs157 =             clhs139*clhs55;
const double clhs158 =             clhs141*clhs55;
const double clhs159 =             clhs157 + clhs158;
const double clhs160 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs130;
const double clhs161 =             C(3,4)*DN(2,1);
const double clhs162 =             C(2,3)*DN(2,2) + clhs147 + clhs161;
const double clhs163 =             C(2,5)*DN(2,2);
const double clhs164 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs163;
const double clhs165 =             DN(2,2)*clhs38;
const double clhs166 =             2.0*N[2]*clhs16*clhs19*omega[1];
const double clhs167 =             clhs166*clhs23;
const double clhs168 =             clhs153*clhs35;
const double clhs169 =             clhs136*clhs62;
const double clhs170 =             clhs168 + clhs169;
const double clhs171 =             DN(0,0)*N[2];
const double clhs172 =             DN(2,1)*omega[2] - DN(2,2)*omega[1];
const double clhs173 =             DN(2,0)*N[0];
const double clhs174 =             1.0*DN(2,0)*clhs19*clhs6;
const double clhs175 =             DN(3,0)*clhs38;
const double clhs176 =             4.0*N[0]*N[3]*clhs16*clhs19;
const double clhs177 =             clhs14*clhs176;
const double clhs178 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs179 =             C(0,3)*DN(3,0);
const double clhs180 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs179;
const double clhs181 =             C(0,5)*DN(3,0);
const double clhs182 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs181;
const double clhs183 =             N[3]*clhs81;
const double clhs184 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs185 =             clhs184*clhs22;
const double clhs186 =             N[3]*bdf0 + clhs184;
const double clhs187 =             clhs186*clhs26;
const double clhs188 =             clhs186*clhs27;
const double clhs189 =             N[0]*N[3]*clhs6;
const double clhs190 =             clhs189*clhs35;
const double clhs191 =             N[0]*N[3]*clhs16*clhs19*clhs25;
const double clhs192 =             clhs191*clhs35;
const double clhs193 =             -clhs190 - clhs192;
const double clhs194 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs179;
const double clhs195 =             C(1,3)*DN(3,1);
const double clhs196 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs195;
const double clhs197 =             C(3,5)*DN(3,0);
const double clhs198 =             C(4,5)*DN(3,2);
const double clhs199 =             C(1,5)*DN(3,1) + clhs197 + clhs198;
const double clhs200 =             DN(3,1)*clhs38;
const double clhs201 =             2.0*N[3]*clhs16*clhs19*omega[2];
const double clhs202 =             clhs201*clhs23;
const double clhs203 =             N[3]*omega[0];
const double clhs204 =             clhs203*clhs55;
const double clhs205 =             clhs186*clhs46;
const double clhs206 =             clhs204 - clhs205;
const double clhs207 =             clhs189*clhs55;
const double clhs208 =             clhs191*clhs55;
const double clhs209 =             clhs207 + clhs208;
const double clhs210 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs181;
const double clhs211 =             C(3,4)*DN(3,1);
const double clhs212 =             C(2,3)*DN(3,2) + clhs197 + clhs211;
const double clhs213 =             C(2,5)*DN(3,2);
const double clhs214 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs213;
const double clhs215 =             DN(3,2)*clhs38;
const double clhs216 =             2.0*N[3]*clhs16*clhs19*omega[1];
const double clhs217 =             clhs216*clhs23;
const double clhs218 =             clhs203*clhs35;
const double clhs219 =             clhs186*clhs62;
const double clhs220 =             clhs218 + clhs219;
const double clhs221 =             DN(0,0)*N[3];
const double clhs222 =             DN(3,1)*omega[2] - DN(3,2)*omega[1];
const double clhs223 =             DN(3,0)*N[0];
const double clhs224 =             1.0*DN(3,0)*clhs19*clhs6;
const double clhs225 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs30;
const double clhs226 =             C(0,4)*DN(0,0) + clhs33 + clhs51;
const double clhs227 =             clhs45 + clhs47;
const double clhs228 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs229 =             C(1,4)*DN(0,1);
const double clhs230 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs229;
const double clhs231 =             pow(DN(0,1), 2);
const double clhs232 =             clhs13 + clhs17;
const double clhs233 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs229;
const double clhs234 =             C(2,4)*DN(0,2);
const double clhs235 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs234;
const double clhs236 =             2.0*omega[0];
const double clhs237 =             clhs236*clhs36;
const double clhs238 =             DN(0,1)*clhs11;
const double clhs239 =             DN(0,2)*clhs238;
const double clhs240 =             clhs236*clhs40;
const double clhs241 =             2.0*N[0]*clhs16*clhs19*omega[0];
const double clhs242 =             clhs23*clhs241;
const double clhs243 =             2.0*omega[1]*omega[2];
const double clhs244 =             N[0]*clhs243;
const double clhs245 =             1.0*omega[0];
const double clhs246 =             clhs24*clhs245;
const double clhs247 =             clhs244 - clhs246;
const double clhs248 =             DN(0,1)*N[0];
const double clhs249 =             DN(0,0)*omega[2] - DN(0,2)*omega[0];
const double clhs250 =             1.0*DN(0,1)*clhs19*clhs6;
const double clhs251 =             clhs89 + clhs91;
const double clhs252 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs94;
const double clhs253 =             C(0,4)*DN(1,0) + clhs110 + clhs97;
const double clhs254 =             DN(1,0)*clhs238;
const double clhs255 =             clhs103 + clhs104;
const double clhs256 =             DN(1,1)*clhs238;
const double clhs257 =             clhs232*clhs74;
const double clhs258 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs259 =             C(1,4)*DN(1,1);
const double clhs260 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs259;
const double clhs261 =             N[1]*clhs6;
const double clhs262 =             clhs261*clhs44;
const double clhs263 =             N[1]*clhs16*clhs19*clhs25*clhs44;
const double clhs264 =             -clhs262 - clhs263;
const double clhs265 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs259;
const double clhs266 =             C(2,4)*DN(1,2);
const double clhs267 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs266;
const double clhs268 =             DN(1,2)*clhs238;
const double clhs269 =             2.0*N[1]*clhs16*clhs19*omega[0];
const double clhs270 =             clhs23*clhs269;
const double clhs271 =             N[1]*clhs243;
const double clhs272 =             clhs245*clhs85;
const double clhs273 =             clhs271 - clhs272;
const double clhs274 =             DN(0,1)*N[1];
const double clhs275 =             DN(1,0)*omega[2] - DN(1,2)*omega[0];
const double clhs276 =             DN(1,1)*N[0];
const double clhs277 =             1.0*DN(1,1)*clhs19*clhs6;
const double clhs278 =             clhs140 + clhs142;
const double clhs279 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs145;
const double clhs280 =             C(0,4)*DN(2,0) + clhs148 + clhs161;
const double clhs281 =             DN(2,0)*clhs238;
const double clhs282 =             clhs154 + clhs155;
const double clhs283 =             DN(2,1)*clhs238;
const double clhs284 =             clhs125*clhs232;
const double clhs285 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs286 =             C(1,4)*DN(2,1);
const double clhs287 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs286;
const double clhs288 =             N[2]*clhs6;
const double clhs289 =             clhs288*clhs44;
const double clhs290 =             N[2]*clhs16*clhs19*clhs25;
const double clhs291 =             clhs290*clhs44;
const double clhs292 =             -clhs289 - clhs291;
const double clhs293 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs286;
const double clhs294 =             C(2,4)*DN(2,2);
const double clhs295 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs294;
const double clhs296 =             DN(2,2)*clhs238;
const double clhs297 =             2.0*N[2]*clhs16*clhs19*omega[0];
const double clhs298 =             clhs23*clhs297;
const double clhs299 =             N[2]*clhs243;
const double clhs300 =             clhs136*clhs245;
const double clhs301 =             clhs299 - clhs300;
const double clhs302 =             DN(0,1)*N[2];
const double clhs303 =             DN(2,0)*omega[2] - DN(2,2)*omega[0];
const double clhs304 =             DN(2,1)*N[0];
const double clhs305 =             1.0*DN(2,1)*clhs19*clhs6;
const double clhs306 =             clhs190 + clhs192;
const double clhs307 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs195;
const double clhs308 =             C(0,4)*DN(3,0) + clhs198 + clhs211;
const double clhs309 =             DN(3,0)*clhs238;
const double clhs310 =             clhs204 + clhs205;
const double clhs311 =             DN(3,1)*clhs238;
const double clhs312 =             clhs176*clhs232;
const double clhs313 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs314 =             C(1,4)*DN(3,1);
const double clhs315 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs314;
const double clhs316 =             N[3]*clhs6;
const double clhs317 =             clhs316*clhs44;
const double clhs318 =             N[3]*clhs16*clhs19*clhs25;
const double clhs319 =             clhs318*clhs44;
const double clhs320 =             -clhs317 - clhs319;
const double clhs321 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs314;
const double clhs322 =             C(2,4)*DN(3,2);
const double clhs323 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs322;
const double clhs324 =             DN(3,2)*clhs238;
const double clhs325 =             2.0*N[3]*clhs16*clhs19*omega[0];
const double clhs326 =             clhs23*clhs325;
const double clhs327 =             N[3]*clhs243;
const double clhs328 =             clhs186*clhs245;
const double clhs329 =             clhs327 - clhs328;
const double clhs330 =             DN(0,1)*N[3];
const double clhs331 =             DN(3,0)*omega[2] - DN(3,2)*omega[0];
const double clhs332 =             DN(3,1)*N[0];
const double clhs333 =             1.0*DN(3,1)*clhs19*clhs6;
const double clhs334 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs53;
const double clhs335 =             clhs61 - clhs63;
const double clhs336 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs234;
const double clhs337 =             clhs244 + clhs246;
const double clhs338 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs339 =             pow(DN(0,2), 2);
const double clhs340 =             DN(0,2)*N[0];
const double clhs341 =             DN(0,0)*omega[1] - DN(0,1)*omega[0];
const double clhs342 =             1.0*DN(0,2)*clhs19*clhs6;
const double clhs343 =             -clhs106 - clhs107;
const double clhs344 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs112;
const double clhs345 =             DN(0,2)*clhs11;
const double clhs346 =             DN(1,0)*clhs345;
const double clhs347 =             clhs117 - clhs118;
const double clhs348 =             clhs262 + clhs263;
const double clhs349 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs266;
const double clhs350 =             DN(1,1)*clhs345;
const double clhs351 =             clhs271 + clhs272;
const double clhs352 =             DN(1,2)*clhs345;
const double clhs353 =             clhs18*clhs74;
const double clhs354 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs355 =             DN(0,2)*N[1];
const double clhs356 =             DN(1,0)*omega[1] - DN(1,1)*omega[0];
const double clhs357 =             DN(1,2)*N[0];
const double clhs358 =             1.0*DN(1,2)*clhs19*clhs6;
const double clhs359 =             -clhs157 - clhs158;
const double clhs360 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs163;
const double clhs361 =             DN(2,0)*clhs345;
const double clhs362 =             clhs168 - clhs169;
const double clhs363 =             clhs289 + clhs291;
const double clhs364 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs294;
const double clhs365 =             DN(2,1)*clhs345;
const double clhs366 =             clhs299 + clhs300;
const double clhs367 =             DN(2,2)*clhs345;
const double clhs368 =             clhs125*clhs18;
const double clhs369 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs370 =             DN(0,2)*N[2];
const double clhs371 =             DN(2,0)*omega[1] - DN(2,1)*omega[0];
const double clhs372 =             DN(2,2)*N[0];
const double clhs373 =             1.0*DN(2,2)*clhs19*clhs6;
const double clhs374 =             -clhs207 - clhs208;
const double clhs375 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs213;
const double clhs376 =             DN(3,0)*clhs345;
const double clhs377 =             clhs218 - clhs219;
const double clhs378 =             clhs317 + clhs319;
const double clhs379 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs322;
const double clhs380 =             DN(3,1)*clhs345;
const double clhs381 =             clhs327 + clhs328;
const double clhs382 =             DN(3,2)*clhs345;
const double clhs383 =             clhs176*clhs18;
const double clhs384 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs385 =             DN(0,2)*N[3];
const double clhs386 =             DN(3,0)*omega[1] - DN(3,1)*omega[0];
const double clhs387 =             DN(3,2)*N[0];
const double clhs388 =             1.0*DN(3,2)*clhs19*clhs6;
const double clhs389 =             DN(0,1)*N[0]*clhs19*clhs6;
const double clhs390 =             DN(0,2)*N[0]*clhs19*clhs6;
const double clhs391 =             DN(0,0)*N[0]*clhs19*clhs6;
const double clhs392 =             bdf0*clhs66*clhs67;
const double clhs393 =             1.0*clhs19;
const double clhs394 =             DN(0,1)*N[1]*clhs19*clhs6;
const double clhs395 =             DN(0,2)*N[1]*clhs19*clhs6;
const double clhs396 =             DN(0,0)*N[1]*clhs19*clhs6;
const double clhs397 =             N[0]*bdf0*clhs66*clhs67;
const double clhs398 =             1.0*DN(0,0)*clhs19;
const double clhs399 =             1.0*DN(0,1)*clhs19;
const double clhs400 =             1.0*DN(0,2)*clhs19;
const double clhs401 =             DN(1,0)*clhs398 + DN(1,1)*clhs399 + DN(1,2)*clhs400 + N[1]*clhs397;
const double clhs402 =             DN(0,1)*N[2]*clhs19*clhs6;
const double clhs403 =             DN(0,2)*N[2]*clhs19*clhs6;
const double clhs404 =             DN(0,0)*N[2]*clhs19*clhs6;
const double clhs405 =             DN(2,0)*clhs398 + DN(2,1)*clhs399 + DN(2,2)*clhs400 + N[2]*clhs397;
const double clhs406 =             DN(0,1)*N[3]*clhs19*clhs6;
const double clhs407 =             DN(0,2)*N[3]*clhs19*clhs6;
const double clhs408 =             DN(0,0)*N[3]*clhs19*clhs6;
const double clhs409 =             DN(3,0)*clhs398 + DN(3,1)*clhs399 + DN(3,2)*clhs400 + N[3]*clhs397;
const double clhs410 =             1.0*N[1]*clhs16*clhs19*clhs25;
const double clhs411 =             1.0*clhs16*clhs19*clhs83;
const double clhs412 =             clhs23*clhs261 + clhs24*clhs410 + clhs24*clhs411 + clhs82;
const double clhs413 =             clhs42*clhs83;
const double clhs414 =             2.0*N[1]*clhs16*clhs19;
const double clhs415 =             clhs59*clhs83;
const double clhs416 =             2.0*N[1]*clhs19*clhs6;
const double clhs417 =             pow(DN(1,0), 2);
const double clhs418 =             pow(N[1], 2);
const double clhs419 =             4.0*clhs16*clhs19*clhs418;
const double clhs420 =             clhs21*clhs418 + clhs261*clhs83 + clhs410*clhs85 + clhs411*clhs85;
const double clhs421 =             clhs418*clhs6;
const double clhs422 =             clhs35*clhs421;
const double clhs423 =             DN(1,0)*clhs11;
const double clhs424 =             DN(1,1)*clhs423;
const double clhs425 =             clhs16*clhs19*clhs25*clhs418;
const double clhs426 =             clhs35*clhs425;
const double clhs427 =             clhs100*clhs83;
const double clhs428 =             clhs421*clhs55;
const double clhs429 =             DN(1,2)*clhs423;
const double clhs430 =             clhs425*clhs55;
const double clhs431 =             clhs115*clhs83;
const double clhs432 =             DN(1,0)*N[1];
const double clhs433 =             DN(2,0)*clhs423;
const double clhs434 =             4.0*N[1]*N[2]*clhs16*clhs19;
const double clhs435 =             clhs14*clhs434;
const double clhs436 =             N[1]*bdf0*clhs6;
const double clhs437 =             N[2]*clhs436;
const double clhs438 =             clhs133*clhs261;
const double clhs439 =             clhs136*clhs410;
const double clhs440 =             clhs136*clhs411;
const double clhs441 =             N[1]*N[2]*clhs6;
const double clhs442 =             clhs35*clhs441;
const double clhs443 =             N[1]*N[2]*clhs16*clhs19*clhs25;
const double clhs444 =             clhs35*clhs443;
const double clhs445 =             -clhs442 - clhs444;
const double clhs446 =             DN(2,1)*clhs423;
const double clhs447 =             clhs151*clhs83;
const double clhs448 =             clhs441*clhs55;
const double clhs449 =             clhs443*clhs55;
const double clhs450 =             clhs448 + clhs449;
const double clhs451 =             DN(2,2)*clhs423;
const double clhs452 =             clhs166*clhs83;
const double clhs453 =             DN(1,0)*N[2];
const double clhs454 =             DN(2,0)*N[1];
const double clhs455 =             DN(3,0)*clhs423;
const double clhs456 =             4.0*N[1]*N[3]*clhs16*clhs19;
const double clhs457 =             clhs14*clhs456;
const double clhs458 =             N[3]*clhs436;
const double clhs459 =             clhs184*clhs261;
const double clhs460 =             clhs186*clhs410;
const double clhs461 =             clhs186*clhs411;
const double clhs462 =             N[1]*N[3]*clhs6;
const double clhs463 =             clhs35*clhs462;
const double clhs464 =             N[1]*N[3]*clhs16*clhs19*clhs25;
const double clhs465 =             clhs35*clhs464;
const double clhs466 =             -clhs463 - clhs465;
const double clhs467 =             DN(3,1)*clhs423;
const double clhs468 =             clhs201*clhs83;
const double clhs469 =             clhs462*clhs55;
const double clhs470 =             clhs464*clhs55;
const double clhs471 =             clhs469 + clhs470;
const double clhs472 =             DN(3,2)*clhs423;
const double clhs473 =             clhs216*clhs83;
const double clhs474 =             DN(1,0)*N[3];
const double clhs475 =             DN(3,0)*N[1];
const double clhs476 =             clhs241*clhs83;
const double clhs477 =             pow(DN(1,1), 2);
const double clhs478 =             clhs236*clhs421;
const double clhs479 =             DN(1,1)*clhs11;
const double clhs480 =             DN(1,2)*clhs479;
const double clhs481 =             clhs236*clhs425;
const double clhs482 =             clhs269*clhs83;
const double clhs483 =             DN(1,1)*N[1];
const double clhs484 =             clhs442 + clhs444;
const double clhs485 =             DN(2,0)*clhs479;
const double clhs486 =             DN(2,1)*clhs479;
const double clhs487 =             clhs232*clhs434;
const double clhs488 =             2.0*N[1]*omega[0];
const double clhs489 =             clhs288*clhs488;
const double clhs490 =             clhs290*clhs488;
const double clhs491 =             -clhs489 - clhs490;
const double clhs492 =             DN(2,2)*clhs479;
const double clhs493 =             clhs297*clhs83;
const double clhs494 =             DN(1,1)*N[2];
const double clhs495 =             DN(2,1)*N[1];
const double clhs496 =             clhs463 + clhs465;
const double clhs497 =             DN(3,0)*clhs479;
const double clhs498 =             DN(3,1)*clhs479;
const double clhs499 =             clhs232*clhs456;
const double clhs500 =             clhs316*clhs488;
const double clhs501 =             clhs318*clhs488;
const double clhs502 =             -clhs500 - clhs501;
const double clhs503 =             DN(3,2)*clhs479;
const double clhs504 =             clhs325*clhs83;
const double clhs505 =             DN(1,1)*N[3];
const double clhs506 =             DN(3,1)*N[1];
const double clhs507 =             pow(DN(1,2), 2);
const double clhs508 =             DN(1,2)*N[1];
const double clhs509 =             -clhs448 - clhs449;
const double clhs510 =             DN(1,2)*clhs11;
const double clhs511 =             DN(2,0)*clhs510;
const double clhs512 =             clhs489 + clhs490;
const double clhs513 =             DN(2,1)*clhs510;
const double clhs514 =             DN(2,2)*clhs510;
const double clhs515 =             clhs18*clhs434;
const double clhs516 =             DN(1,2)*N[2];
const double clhs517 =             DN(2,2)*N[1];
const double clhs518 =             -clhs469 - clhs470;
const double clhs519 =             DN(3,0)*clhs510;
const double clhs520 =             clhs500 + clhs501;
const double clhs521 =             DN(3,1)*clhs510;
const double clhs522 =             DN(3,2)*clhs510;
const double clhs523 =             clhs18*clhs456;
const double clhs524 =             DN(1,2)*N[3];
const double clhs525 =             DN(3,2)*N[1];
const double clhs526 =             DN(1,1)*N[0]*clhs19*clhs6;
const double clhs527 =             DN(1,2)*N[0]*clhs19*clhs6;
const double clhs528 =             DN(1,0)*N[0]*clhs19*clhs6;
const double clhs529 =             DN(1,1)*N[1]*clhs19*clhs6;
const double clhs530 =             DN(1,2)*N[1]*clhs19*clhs6;
const double clhs531 =             DN(1,0)*N[1]*clhs19*clhs6;
const double clhs532 =             DN(1,1)*N[2]*clhs19*clhs6;
const double clhs533 =             DN(1,2)*N[2]*clhs19*clhs6;
const double clhs534 =             DN(1,0)*N[2]*clhs19*clhs6;
const double clhs535 =             N[1]*bdf0*clhs66*clhs67;
const double clhs536 =             1.0*DN(1,0)*clhs19;
const double clhs537 =             1.0*DN(1,1)*clhs19;
const double clhs538 =             1.0*DN(1,2)*clhs19;
const double clhs539 =             DN(2,0)*clhs536 + DN(2,1)*clhs537 + DN(2,2)*clhs538 + N[2]*clhs535;
const double clhs540 =             DN(1,1)*N[3]*clhs19*clhs6;
const double clhs541 =             DN(1,2)*N[3]*clhs19*clhs6;
const double clhs542 =             DN(1,0)*N[3]*clhs19*clhs6;
const double clhs543 =             DN(3,0)*clhs536 + DN(3,1)*clhs537 + DN(3,2)*clhs538 + N[3]*clhs535;
const double clhs544 =             1.0*N[2]*clhs16*clhs19*clhs25;
const double clhs545 =             1.0*clhs133*clhs16*clhs19;
const double clhs546 =             clhs132 + clhs23*clhs288 + clhs24*clhs544 + clhs24*clhs545;
const double clhs547 =             clhs133*clhs42;
const double clhs548 =             2.0*N[2]*clhs16*clhs19;
const double clhs549 =             clhs133*clhs59;
const double clhs550 =             2.0*N[2]*clhs19*clhs6;
const double clhs551 =             clhs288*clhs83 + clhs437 + clhs544*clhs85 + clhs545*clhs85;
const double clhs552 =             clhs100*clhs133;
const double clhs553 =             clhs115*clhs133;
const double clhs554 =             pow(DN(2,0), 2);
const double clhs555 =             pow(N[2], 2);
const double clhs556 =             4.0*clhs16*clhs19*clhs555;
const double clhs557 =             clhs133*clhs288 + clhs136*clhs544 + clhs136*clhs545 + clhs21*clhs555;
const double clhs558 =             clhs555*clhs6;
const double clhs559 =             clhs35*clhs558;
const double clhs560 =             DN(2,0)*clhs11;
const double clhs561 =             DN(2,1)*clhs560;
const double clhs562 =             clhs16*clhs19*clhs25*clhs555;
const double clhs563 =             clhs35*clhs562;
const double clhs564 =             clhs133*clhs151;
const double clhs565 =             clhs55*clhs558;
const double clhs566 =             DN(2,2)*clhs560;
const double clhs567 =             clhs55*clhs562;
const double clhs568 =             clhs133*clhs166;
const double clhs569 =             DN(2,0)*N[2];
const double clhs570 =             DN(3,0)*clhs560;
const double clhs571 =             4.0*N[2]*N[3]*clhs16*clhs19;
const double clhs572 =             clhs14*clhs571;
const double clhs573 =             clhs135*clhs316;
const double clhs574 =             clhs184*clhs288;
const double clhs575 =             clhs186*clhs544;
const double clhs576 =             clhs186*clhs545;
const double clhs577 =             N[2]*N[3]*clhs6;
const double clhs578 =             clhs35*clhs577;
const double clhs579 =             N[2]*N[3]*clhs16*clhs19*clhs25;
const double clhs580 =             clhs35*clhs579;
const double clhs581 =             -clhs578 - clhs580;
const double clhs582 =             DN(3,1)*clhs560;
const double clhs583 =             clhs133*clhs201;
const double clhs584 =             clhs55*clhs577;
const double clhs585 =             clhs55*clhs579;
const double clhs586 =             clhs584 + clhs585;
const double clhs587 =             DN(3,2)*clhs560;
const double clhs588 =             clhs133*clhs216;
const double clhs589 =             DN(2,0)*N[3];
const double clhs590 =             DN(3,0)*N[2];
const double clhs591 =             clhs133*clhs241;
const double clhs592 =             clhs133*clhs269;
const double clhs593 =             pow(DN(2,1), 2);
const double clhs594 =             clhs236*clhs558;
const double clhs595 =             DN(2,1)*clhs11;
const double clhs596 =             DN(2,2)*clhs595;
const double clhs597 =             clhs236*clhs562;
const double clhs598 =             clhs133*clhs297;
const double clhs599 =             DN(2,1)*N[2];
const double clhs600 =             clhs578 + clhs580;
const double clhs601 =             DN(3,0)*clhs595;
const double clhs602 =             DN(3,1)*clhs595;
const double clhs603 =             clhs232*clhs571;
const double clhs604 =             2.0*N[3]*clhs153*clhs6;
const double clhs605 =             2.0*N[2]*clhs318*omega[0];
const double clhs606 =             -clhs604 - clhs605;
const double clhs607 =             DN(3,2)*clhs595;
const double clhs608 =             clhs133*clhs325;
const double clhs609 =             DN(2,1)*N[3];
const double clhs610 =             DN(3,1)*N[2];
const double clhs611 =             pow(DN(2,2), 2);
const double clhs612 =             DN(2,2)*N[2];
const double clhs613 =             -clhs584 - clhs585;
const double clhs614 =             DN(2,2)*clhs11;
const double clhs615 =             DN(3,0)*clhs614;
const double clhs616 =             clhs604 + clhs605;
const double clhs617 =             DN(3,1)*clhs614;
const double clhs618 =             DN(3,2)*clhs614;
const double clhs619 =             clhs18*clhs571;
const double clhs620 =             DN(2,2)*N[3];
const double clhs621 =             DN(3,2)*N[2];
const double clhs622 =             DN(2,1)*N[0]*clhs19*clhs6;
const double clhs623 =             DN(2,2)*N[0]*clhs19*clhs6;
const double clhs624 =             DN(2,0)*N[0]*clhs19*clhs6;
const double clhs625 =             DN(2,1)*N[1]*clhs19*clhs6;
const double clhs626 =             DN(2,2)*N[1]*clhs19*clhs6;
const double clhs627 =             DN(2,0)*N[1]*clhs19*clhs6;
const double clhs628 =             DN(2,1)*N[2]*clhs19*clhs6;
const double clhs629 =             DN(2,2)*N[2]*clhs19*clhs6;
const double clhs630 =             DN(2,0)*N[2]*clhs19*clhs6;
const double clhs631 =             DN(2,1)*N[3]*clhs19*clhs6;
const double clhs632 =             DN(2,2)*N[3]*clhs19*clhs6;
const double clhs633 =             DN(2,0)*N[3]*clhs19*clhs6;
const double clhs634 =             1.0*DN(2,0)*DN(3,0)*clhs19 + 1.0*DN(2,1)*DN(3,1)*clhs19 + 1.0*DN(2,2)*DN(3,2)*clhs19 + N[2]*N[3]*bdf0*clhs66*clhs67;
const double clhs635 =             1.0*N[3]*clhs16*clhs19*clhs25;
const double clhs636 =             1.0*clhs16*clhs184*clhs19;
const double clhs637 =             clhs183 + clhs23*clhs316 + clhs24*clhs635 + clhs24*clhs636;
const double clhs638 =             clhs184*clhs42;
const double clhs639 =             2.0*N[3]*clhs16*clhs19;
const double clhs640 =             clhs184*clhs59;
const double clhs641 =             2.0*N[3]*clhs19*clhs6;
const double clhs642 =             clhs316*clhs83 + clhs458 + clhs635*clhs85 + clhs636*clhs85;
const double clhs643 =             clhs100*clhs184;
const double clhs644 =             clhs115*clhs184;
const double clhs645 =             clhs133*clhs316 + clhs136*clhs635 + clhs136*clhs636 + clhs573;
const double clhs646 =             clhs151*clhs184;
const double clhs647 =             clhs166*clhs184;
const double clhs648 =             pow(DN(3,0), 2);
const double clhs649 =             pow(N[3], 2);
const double clhs650 =             4.0*clhs16*clhs19*clhs649;
const double clhs651 =             clhs184*clhs316 + clhs186*clhs635 + clhs186*clhs636 + clhs21*clhs649;
const double clhs652 =             clhs6*clhs649;
const double clhs653 =             clhs35*clhs652;
const double clhs654 =             DN(3,0)*clhs11;
const double clhs655 =             DN(3,1)*clhs654;
const double clhs656 =             clhs16*clhs19*clhs25*clhs649;
const double clhs657 =             clhs35*clhs656;
const double clhs658 =             clhs184*clhs201;
const double clhs659 =             clhs55*clhs652;
const double clhs660 =             DN(3,2)*clhs654;
const double clhs661 =             clhs55*clhs656;
const double clhs662 =             clhs184*clhs216;
const double clhs663 =             DN(3,0)*N[3];
const double clhs664 =             clhs184*clhs241;
const double clhs665 =             clhs184*clhs269;
const double clhs666 =             clhs184*clhs297;
const double clhs667 =             pow(DN(3,1), 2);
const double clhs668 =             clhs236*clhs652;
const double clhs669 =             DN(3,1)*DN(3,2)*clhs11;
const double clhs670 =             clhs236*clhs656;
const double clhs671 =             clhs184*clhs325;
const double clhs672 =             DN(3,1)*N[3];
const double clhs673 =             pow(DN(3,2), 2);
const double clhs674 =             DN(3,2)*N[3];
const double clhs675 =             DN(3,1)*N[0]*clhs19*clhs6;
const double clhs676 =             DN(3,2)*N[0]*clhs19*clhs6;
const double clhs677 =             DN(3,0)*N[0]*clhs19*clhs6;
const double clhs678 =             DN(3,1)*N[1]*clhs19*clhs6;
const double clhs679 =             DN(3,2)*N[1]*clhs19*clhs6;
const double clhs680 =             DN(3,0)*N[1]*clhs19*clhs6;
const double clhs681 =             DN(3,1)*N[2]*clhs19*clhs6;
const double clhs682 =             DN(3,2)*N[2]*clhs19*clhs6;
const double clhs683 =             DN(3,0)*N[2]*clhs19*clhs6;
const double clhs684 =             DN(3,1)*N[3]*clhs19*clhs6;
const double clhs685 =             DN(3,2)*N[3]*clhs19*clhs6;
const double clhs686 =             DN(3,0)*N[3]*clhs19*clhs6;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs14*clhs20 + clhs28;
            lhs(0,1)=DN(0,0)*clhs29 + DN(0,1)*clhs31 + DN(0,2)*clhs34 - clhs37 + clhs39 - clhs41 - clhs43 - clhs48*clhs49;
            lhs(0,2)=DN(0,0)*clhs50 + DN(0,1)*clhs52 + DN(0,2)*clhs54 - clhs49*clhs64 + clhs56 + clhs57 + clhs58 + clhs60;
            lhs(0,3)=clhs23*clhs72 + clhs65*clhs68 + clhs65*clhs71 - clhs65 + clhs69*clhs70;
            lhs(0,4)=DN(0,0)*clhs76 + DN(0,1)*clhs78 + DN(0,2)*clhs80 + clhs73 + clhs75 + clhs82 + clhs84 + clhs86 + clhs87;
            lhs(0,5)=DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs98 - clhs101 - clhs105*clhs49 + clhs92 + clhs99;
            lhs(0,6)=DN(0,0)*clhs109 + DN(0,1)*clhs111 + DN(0,2)*clhs113 + clhs108 + clhs114 + clhs116 - clhs119*clhs49;
            lhs(0,7)=clhs120*clhs68 - clhs120 + clhs121*clhs70 + clhs122*clhs71 + clhs123*clhs23;
            lhs(0,8)=DN(0,0)*clhs127 + DN(0,1)*clhs129 + DN(0,2)*clhs131 + clhs124 + clhs126 + clhs132 + clhs134 + clhs137 + clhs138;
            lhs(0,9)=DN(0,0)*clhs144 + DN(0,1)*clhs146 + DN(0,2)*clhs149 + clhs143 + clhs150 - clhs152 - clhs156*clhs49;
            lhs(0,10)=DN(0,0)*clhs160 + DN(0,1)*clhs162 + DN(0,2)*clhs164 + clhs159 + clhs165 + clhs167 - clhs170*clhs49;
            lhs(0,11)=clhs171*clhs68 - clhs171 + clhs172*clhs70 + clhs173*clhs71 + clhs174*clhs23;
            lhs(0,12)=DN(0,0)*clhs178 + DN(0,1)*clhs180 + DN(0,2)*clhs182 + clhs175 + clhs177 + clhs183 + clhs185 + clhs187 + clhs188;
            lhs(0,13)=DN(0,0)*clhs194 + DN(0,1)*clhs196 + DN(0,2)*clhs199 + clhs193 + clhs200 - clhs202 - clhs206*clhs49;
            lhs(0,14)=DN(0,0)*clhs210 + DN(0,1)*clhs212 + DN(0,2)*clhs214 + clhs209 + clhs215 + clhs217 - clhs220*clhs49;
            lhs(0,15)=clhs221*clhs68 - clhs221 + clhs222*clhs70 + clhs223*clhs71 + clhs224*clhs23;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs225 + DN(0,2)*clhs226 - clhs227*clhs49 + clhs37 + clhs39 + clhs41 + clhs43;
            lhs(1,1)=DN(0,0)*clhs31 + DN(0,1)*clhs228 + DN(0,2)*clhs230 + clhs11*clhs231 + clhs20*clhs232 + clhs28;
            lhs(1,2)=DN(0,0)*clhs52 + DN(0,1)*clhs233 + DN(0,2)*clhs235 - clhs237 + clhs239 - clhs240 - clhs242 - clhs247*clhs49;
            lhs(1,3)=clhs23*clhs250 + clhs248*clhs68 + clhs248*clhs71 - clhs248 - clhs249*clhs70;
            lhs(1,4)=DN(0,0)*clhs78 + DN(0,1)*clhs252 + DN(0,2)*clhs253 + clhs101 + clhs251 + clhs254 - clhs255*clhs49;
            lhs(1,5)=DN(0,0)*clhs95 + DN(0,1)*clhs258 + DN(0,2)*clhs260 + clhs256 + clhs257 + clhs82 + clhs84 + clhs86 + clhs87;
            lhs(1,6)=DN(0,0)*clhs111 + DN(0,1)*clhs265 + DN(0,2)*clhs267 + clhs264 + clhs268 - clhs270 - clhs273*clhs49;
            lhs(1,7)=clhs23*clhs277 + clhs274*clhs68 - clhs274 - clhs275*clhs70 + clhs276*clhs71;
            lhs(1,8)=DN(0,0)*clhs129 + DN(0,1)*clhs279 + DN(0,2)*clhs280 + clhs152 + clhs278 + clhs281 - clhs282*clhs49;
            lhs(1,9)=DN(0,0)*clhs146 + DN(0,1)*clhs285 + DN(0,2)*clhs287 + clhs132 + clhs134 + clhs137 + clhs138 + clhs283 + clhs284;
            lhs(1,10)=DN(0,0)*clhs162 + DN(0,1)*clhs293 + DN(0,2)*clhs295 + clhs292 + clhs296 - clhs298 - clhs301*clhs49;
            lhs(1,11)=clhs23*clhs305 + clhs302*clhs68 - clhs302 - clhs303*clhs70 + clhs304*clhs71;
            lhs(1,12)=DN(0,0)*clhs180 + DN(0,1)*clhs307 + DN(0,2)*clhs308 + clhs202 + clhs306 + clhs309 - clhs310*clhs49;
            lhs(1,13)=DN(0,0)*clhs196 + DN(0,1)*clhs313 + DN(0,2)*clhs315 + clhs183 + clhs185 + clhs187 + clhs188 + clhs311 + clhs312;
            lhs(1,14)=DN(0,0)*clhs212 + DN(0,1)*clhs321 + DN(0,2)*clhs323 + clhs320 + clhs324 - clhs326 - clhs329*clhs49;
            lhs(1,15)=clhs23*clhs333 + clhs330*clhs68 - clhs330 - clhs331*clhs70 + clhs332*clhs71;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs226 + DN(0,2)*clhs334 - clhs335*clhs49 - clhs56 + clhs57 - clhs58 - clhs60;
            lhs(2,1)=DN(0,0)*clhs34 + DN(0,1)*clhs230 + DN(0,2)*clhs336 + clhs237 + clhs239 + clhs240 + clhs242 - clhs337*clhs49;
            lhs(2,2)=DN(0,0)*clhs54 + DN(0,1)*clhs235 + DN(0,2)*clhs338 + clhs11*clhs339 + clhs18*clhs20 + clhs28;
            lhs(2,3)=clhs23*clhs342 + clhs340*clhs68 + clhs340*clhs71 - clhs340 + clhs341*clhs70;
            lhs(2,4)=DN(0,0)*clhs80 + DN(0,1)*clhs253 + DN(0,2)*clhs344 - clhs116 + clhs343 + clhs346 - clhs347*clhs49;
            lhs(2,5)=DN(0,0)*clhs98 + DN(0,1)*clhs260 + DN(0,2)*clhs349 + clhs270 + clhs348 + clhs350 - clhs351*clhs49;
            lhs(2,6)=DN(0,0)*clhs113 + DN(0,1)*clhs267 + DN(0,2)*clhs354 + clhs352 + clhs353 + clhs82 + clhs84 + clhs86 + clhs87;
            lhs(2,7)=clhs23*clhs358 + clhs355*clhs68 - clhs355 + clhs356*clhs70 + clhs357*clhs71;
            lhs(2,8)=DN(0,0)*clhs131 + DN(0,1)*clhs280 + DN(0,2)*clhs360 - clhs167 + clhs359 + clhs361 - clhs362*clhs49;
            lhs(2,9)=DN(0,0)*clhs149 + DN(0,1)*clhs287 + DN(0,2)*clhs364 + clhs298 + clhs363 + clhs365 - clhs366*clhs49;
            lhs(2,10)=DN(0,0)*clhs164 + DN(0,1)*clhs295 + DN(0,2)*clhs369 + clhs132 + clhs134 + clhs137 + clhs138 + clhs367 + clhs368;
            lhs(2,11)=clhs23*clhs373 + clhs370*clhs68 - clhs370 + clhs371*clhs70 + clhs372*clhs71;
            lhs(2,12)=DN(0,0)*clhs182 + DN(0,1)*clhs308 + DN(0,2)*clhs375 - clhs217 + clhs374 + clhs376 - clhs377*clhs49;
            lhs(2,13)=DN(0,0)*clhs199 + DN(0,1)*clhs315 + DN(0,2)*clhs379 + clhs326 + clhs378 + clhs380 - clhs381*clhs49;
            lhs(2,14)=DN(0,0)*clhs214 + DN(0,1)*clhs323 + DN(0,2)*clhs384 + clhs183 + clhs185 + clhs187 + clhs188 + clhs382 + clhs383;
            lhs(2,15)=clhs23*clhs388 + clhs385*clhs68 - clhs385 + clhs386*clhs70 + clhs387*clhs71;
            lhs(3,0)=clhs24*clhs72 + clhs35*clhs389 - clhs390*clhs55 + clhs65;
            lhs(3,1)=clhs236*clhs390 + clhs24*clhs250 + clhs248 - clhs35*clhs391;
            lhs(3,2)=-clhs236*clhs389 + clhs24*clhs342 + clhs340 + clhs391*clhs55;
            lhs(3,3)=clhs15*clhs392 + clhs231*clhs393 + clhs339*clhs393 + clhs393*clhs5;
            lhs(3,4)=clhs122 + clhs35*clhs394 - clhs395*clhs55 + clhs72*clhs85;
            lhs(3,5)=clhs236*clhs395 + clhs250*clhs85 + clhs276 - clhs35*clhs396;
            lhs(3,6)=-clhs236*clhs394 + clhs342*clhs85 + clhs357 + clhs396*clhs55;
            lhs(3,7)=clhs401;
            lhs(3,8)=clhs136*clhs72 + clhs173 + clhs35*clhs402 - clhs403*clhs55;
            lhs(3,9)=clhs136*clhs250 + clhs236*clhs403 + clhs304 - clhs35*clhs404;
            lhs(3,10)=clhs136*clhs342 - clhs236*clhs402 + clhs372 + clhs404*clhs55;
            lhs(3,11)=clhs405;
            lhs(3,12)=clhs186*clhs72 + clhs223 + clhs35*clhs406 - clhs407*clhs55;
            lhs(3,13)=clhs186*clhs250 + clhs236*clhs407 + clhs332 - clhs35*clhs408;
            lhs(3,14)=clhs186*clhs342 - clhs236*clhs406 + clhs387 + clhs408*clhs55;
            lhs(3,15)=clhs409;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs412 + clhs73 + clhs75;
            lhs(4,1)=DN(1,0)*clhs29 + DN(1,1)*clhs31 + DN(1,2)*clhs34 + clhs254 - clhs413 - clhs414*clhs48 + clhs92;
            lhs(4,2)=DN(1,0)*clhs50 + DN(1,1)*clhs52 + DN(1,2)*clhs54 + clhs108 + clhs346 - clhs414*clhs64 + clhs415;
            lhs(4,3)=clhs120*clhs71 + clhs122*clhs68 - clhs122 + clhs416*clhs69 + clhs72*clhs83;
            lhs(4,4)=DN(1,0)*clhs76 + DN(1,1)*clhs78 + DN(1,2)*clhs80 + clhs11*clhs417 + clhs14*clhs419 + clhs420;
            lhs(4,5)=DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs98 - clhs105*clhs414 - clhs422 + clhs424 - clhs426 - clhs427;
            lhs(4,6)=DN(1,0)*clhs109 + DN(1,1)*clhs111 + DN(1,2)*clhs113 - clhs119*clhs414 + clhs428 + clhs429 + clhs430 + clhs431;
            lhs(4,7)=clhs121*clhs416 + clhs123*clhs83 + clhs432*clhs68 + clhs432*clhs71 - clhs432;
            lhs(4,8)=DN(1,0)*clhs127 + DN(1,1)*clhs129 + DN(1,2)*clhs131 + clhs433 + clhs435 + clhs437 + clhs438 + clhs439 + clhs440;
            lhs(4,9)=DN(1,0)*clhs144 + DN(1,1)*clhs146 + DN(1,2)*clhs149 - clhs156*clhs414 + clhs445 + clhs446 - clhs447;
            lhs(4,10)=DN(1,0)*clhs160 + DN(1,1)*clhs162 + DN(1,2)*clhs164 - clhs170*clhs414 + clhs450 + clhs451 + clhs452;
            lhs(4,11)=clhs172*clhs416 + clhs174*clhs83 + clhs453*clhs68 - clhs453 + clhs454*clhs71;
            lhs(4,12)=DN(1,0)*clhs178 + DN(1,1)*clhs180 + DN(1,2)*clhs182 + clhs455 + clhs457 + clhs458 + clhs459 + clhs460 + clhs461;
            lhs(4,13)=DN(1,0)*clhs194 + DN(1,1)*clhs196 + DN(1,2)*clhs199 - clhs206*clhs414 + clhs466 + clhs467 - clhs468;
            lhs(4,14)=DN(1,0)*clhs210 + DN(1,1)*clhs212 + DN(1,2)*clhs214 - clhs220*clhs414 + clhs471 + clhs472 + clhs473;
            lhs(4,15)=clhs222*clhs416 + clhs224*clhs83 + clhs474*clhs68 - clhs474 + clhs475*clhs71;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs225 + DN(1,2)*clhs226 - clhs227*clhs414 + clhs251 + clhs413 + clhs99;
            lhs(5,1)=DN(1,0)*clhs31 + DN(1,1)*clhs228 + DN(1,2)*clhs230 + clhs256 + clhs257 + clhs412;
            lhs(5,2)=DN(1,0)*clhs52 + DN(1,1)*clhs233 + DN(1,2)*clhs235 - clhs247*clhs414 + clhs264 + clhs350 - clhs476;
            lhs(5,3)=-clhs249*clhs416 + clhs250*clhs83 + clhs274*clhs71 + clhs276*clhs68 - clhs276;
            lhs(5,4)=DN(1,0)*clhs78 + DN(1,1)*clhs252 + DN(1,2)*clhs253 - clhs255*clhs414 + clhs422 + clhs424 + clhs426 + clhs427;
            lhs(5,5)=DN(1,0)*clhs95 + DN(1,1)*clhs258 + DN(1,2)*clhs260 + clhs11*clhs477 + clhs232*clhs419 + clhs420;
            lhs(5,6)=DN(1,0)*clhs111 + DN(1,1)*clhs265 + DN(1,2)*clhs267 - clhs273*clhs414 - clhs478 + clhs480 - clhs481 - clhs482;
            lhs(5,7)=-clhs275*clhs416 + clhs277*clhs83 + clhs483*clhs68 + clhs483*clhs71 - clhs483;
            lhs(5,8)=DN(1,0)*clhs129 + DN(1,1)*clhs279 + DN(1,2)*clhs280 - clhs282*clhs414 + clhs447 + clhs484 + clhs485;
            lhs(5,9)=DN(1,0)*clhs146 + DN(1,1)*clhs285 + DN(1,2)*clhs287 + clhs437 + clhs438 + clhs439 + clhs440 + clhs486 + clhs487;
            lhs(5,10)=DN(1,0)*clhs162 + DN(1,1)*clhs293 + DN(1,2)*clhs295 - clhs301*clhs414 + clhs491 + clhs492 - clhs493;
            lhs(5,11)=-clhs303*clhs416 + clhs305*clhs83 + clhs494*clhs68 - clhs494 + clhs495*clhs71;
            lhs(5,12)=DN(1,0)*clhs180 + DN(1,1)*clhs307 + DN(1,2)*clhs308 - clhs310*clhs414 + clhs468 + clhs496 + clhs497;
            lhs(5,13)=DN(1,0)*clhs196 + DN(1,1)*clhs313 + DN(1,2)*clhs315 + clhs458 + clhs459 + clhs460 + clhs461 + clhs498 + clhs499;
            lhs(5,14)=DN(1,0)*clhs212 + DN(1,1)*clhs321 + DN(1,2)*clhs323 - clhs329*clhs414 + clhs502 + clhs503 - clhs504;
            lhs(5,15)=-clhs331*clhs416 + clhs333*clhs83 + clhs505*clhs68 - clhs505 + clhs506*clhs71;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs226 + DN(1,2)*clhs334 + clhs114 - clhs335*clhs414 + clhs343 - clhs415;
            lhs(6,1)=DN(1,0)*clhs34 + DN(1,1)*clhs230 + DN(1,2)*clhs336 + clhs268 - clhs337*clhs414 + clhs348 + clhs476;
            lhs(6,2)=DN(1,0)*clhs54 + DN(1,1)*clhs235 + DN(1,2)*clhs338 + clhs352 + clhs353 + clhs412;
            lhs(6,3)=clhs341*clhs416 + clhs342*clhs83 + clhs355*clhs71 + clhs357*clhs68 - clhs357;
            lhs(6,4)=DN(1,0)*clhs80 + DN(1,1)*clhs253 + DN(1,2)*clhs344 - clhs347*clhs414 - clhs428 + clhs429 - clhs430 - clhs431;
            lhs(6,5)=DN(1,0)*clhs98 + DN(1,1)*clhs260 + DN(1,2)*clhs349 - clhs351*clhs414 + clhs478 + clhs480 + clhs481 + clhs482;
            lhs(6,6)=DN(1,0)*clhs113 + DN(1,1)*clhs267 + DN(1,2)*clhs354 + clhs11*clhs507 + clhs18*clhs419 + clhs420;
            lhs(6,7)=clhs356*clhs416 + clhs358*clhs83 + clhs508*clhs68 + clhs508*clhs71 - clhs508;
            lhs(6,8)=DN(1,0)*clhs131 + DN(1,1)*clhs280 + DN(1,2)*clhs360 - clhs362*clhs414 - clhs452 + clhs509 + clhs511;
            lhs(6,9)=DN(1,0)*clhs149 + DN(1,1)*clhs287 + DN(1,2)*clhs364 - clhs366*clhs414 + clhs493 + clhs512 + clhs513;
            lhs(6,10)=DN(1,0)*clhs164 + DN(1,1)*clhs295 + DN(1,2)*clhs369 + clhs437 + clhs438 + clhs439 + clhs440 + clhs514 + clhs515;
            lhs(6,11)=clhs371*clhs416 + clhs373*clhs83 + clhs516*clhs68 - clhs516 + clhs517*clhs71;
            lhs(6,12)=DN(1,0)*clhs182 + DN(1,1)*clhs308 + DN(1,2)*clhs375 - clhs377*clhs414 - clhs473 + clhs518 + clhs519;
            lhs(6,13)=DN(1,0)*clhs199 + DN(1,1)*clhs315 + DN(1,2)*clhs379 - clhs381*clhs414 + clhs504 + clhs520 + clhs521;
            lhs(6,14)=DN(1,0)*clhs214 + DN(1,1)*clhs323 + DN(1,2)*clhs384 + clhs458 + clhs459 + clhs460 + clhs461 + clhs522 + clhs523;
            lhs(6,15)=clhs386*clhs416 + clhs388*clhs83 + clhs524*clhs68 - clhs524 + clhs525*clhs71;
            lhs(7,0)=clhs120 + clhs123*clhs24 + clhs35*clhs526 - clhs527*clhs55;
            lhs(7,1)=clhs236*clhs527 + clhs24*clhs277 + clhs274 - clhs35*clhs528;
            lhs(7,2)=-clhs236*clhs526 + clhs24*clhs358 + clhs355 + clhs528*clhs55;
            lhs(7,3)=clhs401;
            lhs(7,4)=clhs123*clhs85 + clhs35*clhs529 + clhs432 - clhs530*clhs55;
            lhs(7,5)=clhs236*clhs530 + clhs277*clhs85 - clhs35*clhs531 + clhs483;
            lhs(7,6)=-clhs236*clhs529 + clhs358*clhs85 + clhs508 + clhs531*clhs55;
            lhs(7,7)=clhs392*clhs418 + clhs393*clhs417 + clhs393*clhs477 + clhs393*clhs507;
            lhs(7,8)=clhs123*clhs136 + clhs35*clhs532 + clhs454 - clhs533*clhs55;
            lhs(7,9)=clhs136*clhs277 + clhs236*clhs533 - clhs35*clhs534 + clhs495;
            lhs(7,10)=clhs136*clhs358 - clhs236*clhs532 + clhs517 + clhs534*clhs55;
            lhs(7,11)=clhs539;
            lhs(7,12)=clhs123*clhs186 + clhs35*clhs540 + clhs475 - clhs541*clhs55;
            lhs(7,13)=clhs186*clhs277 + clhs236*clhs541 - clhs35*clhs542 + clhs506;
            lhs(7,14)=clhs186*clhs358 - clhs236*clhs540 + clhs525 + clhs542*clhs55;
            lhs(7,15)=clhs543;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs124 + clhs126 + clhs546;
            lhs(8,1)=DN(2,0)*clhs29 + DN(2,1)*clhs31 + DN(2,2)*clhs34 + clhs143 + clhs281 - clhs48*clhs548 - clhs547;
            lhs(8,2)=DN(2,0)*clhs50 + DN(2,1)*clhs52 + DN(2,2)*clhs54 + clhs159 + clhs361 - clhs548*clhs64 + clhs549;
            lhs(8,3)=clhs133*clhs72 + clhs171*clhs71 + clhs173*clhs68 - clhs173 + clhs550*clhs69;
            lhs(8,4)=DN(2,0)*clhs76 + DN(2,1)*clhs78 + DN(2,2)*clhs80 + clhs433 + clhs435 + clhs551;
            lhs(8,5)=DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs98 - clhs105*clhs548 + clhs445 + clhs485 - clhs552;
            lhs(8,6)=DN(2,0)*clhs109 + DN(2,1)*clhs111 + DN(2,2)*clhs113 - clhs119*clhs548 + clhs450 + clhs511 + clhs553;
            lhs(8,7)=clhs121*clhs550 + clhs123*clhs133 + clhs453*clhs71 + clhs454*clhs68 - clhs454;
            lhs(8,8)=DN(2,0)*clhs127 + DN(2,1)*clhs129 + DN(2,2)*clhs131 + clhs11*clhs554 + clhs14*clhs556 + clhs557;
            lhs(8,9)=DN(2,0)*clhs144 + DN(2,1)*clhs146 + DN(2,2)*clhs149 - clhs156*clhs548 - clhs559 + clhs561 - clhs563 - clhs564;
            lhs(8,10)=DN(2,0)*clhs160 + DN(2,1)*clhs162 + DN(2,2)*clhs164 - clhs170*clhs548 + clhs565 + clhs566 + clhs567 + clhs568;
            lhs(8,11)=clhs133*clhs174 + clhs172*clhs550 + clhs569*clhs68 + clhs569*clhs71 - clhs569;
            lhs(8,12)=DN(2,0)*clhs178 + DN(2,1)*clhs180 + DN(2,2)*clhs182 + clhs570 + clhs572 + clhs573 + clhs574 + clhs575 + clhs576;
            lhs(8,13)=DN(2,0)*clhs194 + DN(2,1)*clhs196 + DN(2,2)*clhs199 - clhs206*clhs548 + clhs581 + clhs582 - clhs583;
            lhs(8,14)=DN(2,0)*clhs210 + DN(2,1)*clhs212 + DN(2,2)*clhs214 - clhs220*clhs548 + clhs586 + clhs587 + clhs588;
            lhs(8,15)=clhs133*clhs224 + clhs222*clhs550 + clhs589*clhs68 - clhs589 + clhs590*clhs71;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs225 + DN(2,2)*clhs226 + clhs150 - clhs227*clhs548 + clhs278 + clhs547;
            lhs(9,1)=DN(2,0)*clhs31 + DN(2,1)*clhs228 + DN(2,2)*clhs230 + clhs283 + clhs284 + clhs546;
            lhs(9,2)=DN(2,0)*clhs52 + DN(2,1)*clhs233 + DN(2,2)*clhs235 - clhs247*clhs548 + clhs292 + clhs365 - clhs591;
            lhs(9,3)=clhs133*clhs250 - clhs249*clhs550 + clhs302*clhs71 + clhs304*clhs68 - clhs304;
            lhs(9,4)=DN(2,0)*clhs78 + DN(2,1)*clhs252 + DN(2,2)*clhs253 - clhs255*clhs548 + clhs446 + clhs484 + clhs552;
            lhs(9,5)=DN(2,0)*clhs95 + DN(2,1)*clhs258 + DN(2,2)*clhs260 + clhs486 + clhs487 + clhs551;
            lhs(9,6)=DN(2,0)*clhs111 + DN(2,1)*clhs265 + DN(2,2)*clhs267 - clhs273*clhs548 + clhs491 + clhs513 - clhs592;
            lhs(9,7)=clhs133*clhs277 - clhs275*clhs550 + clhs494*clhs71 + clhs495*clhs68 - clhs495;
            lhs(9,8)=DN(2,0)*clhs129 + DN(2,1)*clhs279 + DN(2,2)*clhs280 - clhs282*clhs548 + clhs559 + clhs561 + clhs563 + clhs564;
            lhs(9,9)=DN(2,0)*clhs146 + DN(2,1)*clhs285 + DN(2,2)*clhs287 + clhs11*clhs593 + clhs232*clhs556 + clhs557;
            lhs(9,10)=DN(2,0)*clhs162 + DN(2,1)*clhs293 + DN(2,2)*clhs295 - clhs301*clhs548 - clhs594 + clhs596 - clhs597 - clhs598;
            lhs(9,11)=clhs133*clhs305 - clhs303*clhs550 + clhs599*clhs68 + clhs599*clhs71 - clhs599;
            lhs(9,12)=DN(2,0)*clhs180 + DN(2,1)*clhs307 + DN(2,2)*clhs308 - clhs310*clhs548 + clhs583 + clhs600 + clhs601;
            lhs(9,13)=DN(2,0)*clhs196 + DN(2,1)*clhs313 + DN(2,2)*clhs315 + clhs573 + clhs574 + clhs575 + clhs576 + clhs602 + clhs603;
            lhs(9,14)=DN(2,0)*clhs212 + DN(2,1)*clhs321 + DN(2,2)*clhs323 - clhs329*clhs548 + clhs606 + clhs607 - clhs608;
            lhs(9,15)=clhs133*clhs333 - clhs331*clhs550 + clhs609*clhs68 - clhs609 + clhs610*clhs71;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs226 + DN(2,2)*clhs334 + clhs165 - clhs335*clhs548 + clhs359 - clhs549;
            lhs(10,1)=DN(2,0)*clhs34 + DN(2,1)*clhs230 + DN(2,2)*clhs336 + clhs296 - clhs337*clhs548 + clhs363 + clhs591;
            lhs(10,2)=DN(2,0)*clhs54 + DN(2,1)*clhs235 + DN(2,2)*clhs338 + clhs367 + clhs368 + clhs546;
            lhs(10,3)=clhs133*clhs342 + clhs341*clhs550 + clhs370*clhs71 + clhs372*clhs68 - clhs372;
            lhs(10,4)=DN(2,0)*clhs80 + DN(2,1)*clhs253 + DN(2,2)*clhs344 - clhs347*clhs548 + clhs451 + clhs509 - clhs553;
            lhs(10,5)=DN(2,0)*clhs98 + DN(2,1)*clhs260 + DN(2,2)*clhs349 - clhs351*clhs548 + clhs492 + clhs512 + clhs592;
            lhs(10,6)=DN(2,0)*clhs113 + DN(2,1)*clhs267 + DN(2,2)*clhs354 + clhs514 + clhs515 + clhs551;
            lhs(10,7)=clhs133*clhs358 + clhs356*clhs550 + clhs516*clhs71 + clhs517*clhs68 - clhs517;
            lhs(10,8)=DN(2,0)*clhs131 + DN(2,1)*clhs280 + DN(2,2)*clhs360 - clhs362*clhs548 - clhs565 + clhs566 - clhs567 - clhs568;
            lhs(10,9)=DN(2,0)*clhs149 + DN(2,1)*clhs287 + DN(2,2)*clhs364 - clhs366*clhs548 + clhs594 + clhs596 + clhs597 + clhs598;
            lhs(10,10)=DN(2,0)*clhs164 + DN(2,1)*clhs295 + DN(2,2)*clhs369 + clhs11*clhs611 + clhs18*clhs556 + clhs557;
            lhs(10,11)=clhs133*clhs373 + clhs371*clhs550 + clhs612*clhs68 + clhs612*clhs71 - clhs612;
            lhs(10,12)=DN(2,0)*clhs182 + DN(2,1)*clhs308 + DN(2,2)*clhs375 - clhs377*clhs548 - clhs588 + clhs613 + clhs615;
            lhs(10,13)=DN(2,0)*clhs199 + DN(2,1)*clhs315 + DN(2,2)*clhs379 - clhs381*clhs548 + clhs608 + clhs616 + clhs617;
            lhs(10,14)=DN(2,0)*clhs214 + DN(2,1)*clhs323 + DN(2,2)*clhs384 + clhs573 + clhs574 + clhs575 + clhs576 + clhs618 + clhs619;
            lhs(10,15)=clhs133*clhs388 + clhs386*clhs550 + clhs620*clhs68 - clhs620 + clhs621*clhs71;
            lhs(11,0)=clhs171 + clhs174*clhs24 + clhs35*clhs622 - clhs55*clhs623;
            lhs(11,1)=clhs236*clhs623 + clhs24*clhs305 + clhs302 - clhs35*clhs624;
            lhs(11,2)=-clhs236*clhs622 + clhs24*clhs373 + clhs370 + clhs55*clhs624;
            lhs(11,3)=clhs405;
            lhs(11,4)=clhs174*clhs85 + clhs35*clhs625 + clhs453 - clhs55*clhs626;
            lhs(11,5)=clhs236*clhs626 + clhs305*clhs85 - clhs35*clhs627 + clhs494;
            lhs(11,6)=-clhs236*clhs625 + clhs373*clhs85 + clhs516 + clhs55*clhs627;
            lhs(11,7)=clhs539;
            lhs(11,8)=clhs136*clhs174 + clhs35*clhs628 - clhs55*clhs629 + clhs569;
            lhs(11,9)=clhs136*clhs305 + clhs236*clhs629 - clhs35*clhs630 + clhs599;
            lhs(11,10)=clhs136*clhs373 - clhs236*clhs628 + clhs55*clhs630 + clhs612;
            lhs(11,11)=clhs392*clhs555 + clhs393*clhs554 + clhs393*clhs593 + clhs393*clhs611;
            lhs(11,12)=clhs174*clhs186 + clhs35*clhs631 - clhs55*clhs632 + clhs590;
            lhs(11,13)=clhs186*clhs305 + clhs236*clhs632 - clhs35*clhs633 + clhs610;
            lhs(11,14)=clhs186*clhs373 - clhs236*clhs631 + clhs55*clhs633 + clhs621;
            lhs(11,15)=clhs634;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs175 + clhs177 + clhs637;
            lhs(12,1)=DN(3,0)*clhs29 + DN(3,1)*clhs31 + DN(3,2)*clhs34 + clhs193 + clhs309 - clhs48*clhs639 - clhs638;
            lhs(12,2)=DN(3,0)*clhs50 + DN(3,1)*clhs52 + DN(3,2)*clhs54 + clhs209 + clhs376 - clhs639*clhs64 + clhs640;
            lhs(12,3)=clhs184*clhs72 + clhs221*clhs71 + clhs223*clhs68 - clhs223 + clhs641*clhs69;
            lhs(12,4)=DN(3,0)*clhs76 + DN(3,1)*clhs78 + DN(3,2)*clhs80 + clhs455 + clhs457 + clhs642;
            lhs(12,5)=DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs98 - clhs105*clhs639 + clhs466 + clhs497 - clhs643;
            lhs(12,6)=DN(3,0)*clhs109 + DN(3,1)*clhs111 + DN(3,2)*clhs113 - clhs119*clhs639 + clhs471 + clhs519 + clhs644;
            lhs(12,7)=clhs121*clhs641 + clhs123*clhs184 + clhs474*clhs71 + clhs475*clhs68 - clhs475;
            lhs(12,8)=DN(3,0)*clhs127 + DN(3,1)*clhs129 + DN(3,2)*clhs131 + clhs570 + clhs572 + clhs645;
            lhs(12,9)=DN(3,0)*clhs144 + DN(3,1)*clhs146 + DN(3,2)*clhs149 - clhs156*clhs639 + clhs581 + clhs601 - clhs646;
            lhs(12,10)=DN(3,0)*clhs160 + DN(3,1)*clhs162 + DN(3,2)*clhs164 - clhs170*clhs639 + clhs586 + clhs615 + clhs647;
            lhs(12,11)=clhs172*clhs641 + clhs174*clhs184 + clhs589*clhs71 + clhs590*clhs68 - clhs590;
            lhs(12,12)=DN(3,0)*clhs178 + DN(3,1)*clhs180 + DN(3,2)*clhs182 + clhs11*clhs648 + clhs14*clhs650 + clhs651;
            lhs(12,13)=DN(3,0)*clhs194 + DN(3,1)*clhs196 + DN(3,2)*clhs199 - clhs206*clhs639 - clhs653 + clhs655 - clhs657 - clhs658;
            lhs(12,14)=DN(3,0)*clhs210 + DN(3,1)*clhs212 + DN(3,2)*clhs214 - clhs220*clhs639 + clhs659 + clhs660 + clhs661 + clhs662;
            lhs(12,15)=clhs184*clhs224 + clhs222*clhs641 + clhs663*clhs68 + clhs663*clhs71 - clhs663;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs225 + DN(3,2)*clhs226 + clhs200 - clhs227*clhs639 + clhs306 + clhs638;
            lhs(13,1)=DN(3,0)*clhs31 + DN(3,1)*clhs228 + DN(3,2)*clhs230 + clhs311 + clhs312 + clhs637;
            lhs(13,2)=DN(3,0)*clhs52 + DN(3,1)*clhs233 + DN(3,2)*clhs235 - clhs247*clhs639 + clhs320 + clhs380 - clhs664;
            lhs(13,3)=clhs184*clhs250 - clhs249*clhs641 + clhs330*clhs71 + clhs332*clhs68 - clhs332;
            lhs(13,4)=DN(3,0)*clhs78 + DN(3,1)*clhs252 + DN(3,2)*clhs253 - clhs255*clhs639 + clhs467 + clhs496 + clhs643;
            lhs(13,5)=DN(3,0)*clhs95 + DN(3,1)*clhs258 + DN(3,2)*clhs260 + clhs498 + clhs499 + clhs642;
            lhs(13,6)=DN(3,0)*clhs111 + DN(3,1)*clhs265 + DN(3,2)*clhs267 - clhs273*clhs639 + clhs502 + clhs521 - clhs665;
            lhs(13,7)=clhs184*clhs277 - clhs275*clhs641 + clhs505*clhs71 + clhs506*clhs68 - clhs506;
            lhs(13,8)=DN(3,0)*clhs129 + DN(3,1)*clhs279 + DN(3,2)*clhs280 - clhs282*clhs639 + clhs582 + clhs600 + clhs646;
            lhs(13,9)=DN(3,0)*clhs146 + DN(3,1)*clhs285 + DN(3,2)*clhs287 + clhs602 + clhs603 + clhs645;
            lhs(13,10)=DN(3,0)*clhs162 + DN(3,1)*clhs293 + DN(3,2)*clhs295 - clhs301*clhs639 + clhs606 + clhs617 - clhs666;
            lhs(13,11)=clhs184*clhs305 - clhs303*clhs641 + clhs609*clhs71 + clhs610*clhs68 - clhs610;
            lhs(13,12)=DN(3,0)*clhs180 + DN(3,1)*clhs307 + DN(3,2)*clhs308 - clhs310*clhs639 + clhs653 + clhs655 + clhs657 + clhs658;
            lhs(13,13)=DN(3,0)*clhs196 + DN(3,1)*clhs313 + DN(3,2)*clhs315 + clhs11*clhs667 + clhs232*clhs650 + clhs651;
            lhs(13,14)=DN(3,0)*clhs212 + DN(3,1)*clhs321 + DN(3,2)*clhs323 - clhs329*clhs639 - clhs668 + clhs669 - clhs670 - clhs671;
            lhs(13,15)=clhs184*clhs333 - clhs331*clhs641 + clhs672*clhs68 + clhs672*clhs71 - clhs672;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs226 + DN(3,2)*clhs334 + clhs215 - clhs335*clhs639 + clhs374 - clhs640;
            lhs(14,1)=DN(3,0)*clhs34 + DN(3,1)*clhs230 + DN(3,2)*clhs336 + clhs324 - clhs337*clhs639 + clhs378 + clhs664;
            lhs(14,2)=DN(3,0)*clhs54 + DN(3,1)*clhs235 + DN(3,2)*clhs338 + clhs382 + clhs383 + clhs637;
            lhs(14,3)=clhs184*clhs342 + clhs341*clhs641 + clhs385*clhs71 + clhs387*clhs68 - clhs387;
            lhs(14,4)=DN(3,0)*clhs80 + DN(3,1)*clhs253 + DN(3,2)*clhs344 - clhs347*clhs639 + clhs472 + clhs518 - clhs644;
            lhs(14,5)=DN(3,0)*clhs98 + DN(3,1)*clhs260 + DN(3,2)*clhs349 - clhs351*clhs639 + clhs503 + clhs520 + clhs665;
            lhs(14,6)=DN(3,0)*clhs113 + DN(3,1)*clhs267 + DN(3,2)*clhs354 + clhs522 + clhs523 + clhs642;
            lhs(14,7)=clhs184*clhs358 + clhs356*clhs641 + clhs524*clhs71 + clhs525*clhs68 - clhs525;
            lhs(14,8)=DN(3,0)*clhs131 + DN(3,1)*clhs280 + DN(3,2)*clhs360 - clhs362*clhs639 + clhs587 + clhs613 - clhs647;
            lhs(14,9)=DN(3,0)*clhs149 + DN(3,1)*clhs287 + DN(3,2)*clhs364 - clhs366*clhs639 + clhs607 + clhs616 + clhs666;
            lhs(14,10)=DN(3,0)*clhs164 + DN(3,1)*clhs295 + DN(3,2)*clhs369 + clhs618 + clhs619 + clhs645;
            lhs(14,11)=clhs184*clhs373 + clhs371*clhs641 + clhs620*clhs71 + clhs621*clhs68 - clhs621;
            lhs(14,12)=DN(3,0)*clhs182 + DN(3,1)*clhs308 + DN(3,2)*clhs375 - clhs377*clhs639 - clhs659 + clhs660 - clhs661 - clhs662;
            lhs(14,13)=DN(3,0)*clhs199 + DN(3,1)*clhs315 + DN(3,2)*clhs379 - clhs381*clhs639 + clhs668 + clhs669 + clhs670 + clhs671;
            lhs(14,14)=DN(3,0)*clhs214 + DN(3,1)*clhs323 + DN(3,2)*clhs384 + clhs11*clhs673 + clhs18*clhs650 + clhs651;
            lhs(14,15)=clhs184*clhs388 + clhs386*clhs641 + clhs674*clhs68 + clhs674*clhs71 - clhs674;
            lhs(15,0)=clhs221 + clhs224*clhs24 + clhs35*clhs675 - clhs55*clhs676;
            lhs(15,1)=clhs236*clhs676 + clhs24*clhs333 + clhs330 - clhs35*clhs677;
            lhs(15,2)=-clhs236*clhs675 + clhs24*clhs388 + clhs385 + clhs55*clhs677;
            lhs(15,3)=clhs409;
            lhs(15,4)=clhs224*clhs85 + clhs35*clhs678 + clhs474 - clhs55*clhs679;
            lhs(15,5)=clhs236*clhs679 + clhs333*clhs85 - clhs35*clhs680 + clhs505;
            lhs(15,6)=-clhs236*clhs678 + clhs388*clhs85 + clhs524 + clhs55*clhs680;
            lhs(15,7)=clhs543;
            lhs(15,8)=clhs136*clhs224 + clhs35*clhs681 - clhs55*clhs682 + clhs589;
            lhs(15,9)=clhs136*clhs333 + clhs236*clhs682 - clhs35*clhs683 + clhs609;
            lhs(15,10)=clhs136*clhs388 - clhs236*clhs681 + clhs55*clhs683 + clhs620;
            lhs(15,11)=clhs634;
            lhs(15,12)=clhs186*clhs224 + clhs35*clhs684 - clhs55*clhs685 + clhs663;
            lhs(15,13)=clhs186*clhs333 + clhs236*clhs685 - clhs35*clhs686 + clhs672;
            lhs(15,14)=clhs186*clhs388 - clhs236*clhs684 + clhs55*clhs686 + clhs674;
            lhs(15,15)=clhs392*clhs649 + clhs393*clhs648 + clhs393*clhs667 + clhs393*clhs673;


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

    // Non-inertial frame or reference data
    const array_1d<double,3>& omega = rData.AngularVelocity;
    const auto& r_geom = this->GetGeometry();
    array_1d<double,3> gauss_pt_coord = ZeroVector(3);
    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_coords = r_geom[i].Coordinates();
        for (IndexType d = 0; d < 2; ++d) {
            gauss_pt_coord(d) += N(i)*r_coords(d);
        }
    }

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crhs2 =             crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs3 =             2.0*crhs1*omega[2];
const double crhs4 =             crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs5 =             gauss_pt_coord[0]*omega[1];
const double crhs6 =             gauss_pt_coord[1]*omega[0];
const double crhs7 =             omega[2]*(gauss_pt_coord[0]*omega[2] - gauss_pt_coord[2]*omega[0]);
const double crhs8 =             crhs1*(crhs7 + omega[1]*(crhs5 - crhs6));
const double crhs9 =             crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs10 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs11 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs12 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs13 =             crhs1*(crhs10*crhs11 + crhs12*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs14 =             crhs1*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs12, 2));
const double crhs15 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs16 =             crhs10 + crhs15;
const double crhs17 =             1.0/crhs1;
const double crhs18 =             crhs17*(crhs11*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crhs12*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crhs19 =             crhs17*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crhs20 =             (crhs14*h/stab_c1 + mu)*(crhs16 + crhs18 + crhs19);
const double crhs21 =             1.0/(crhs1*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2)) + crhs1*dyn_tau/dt + crhs14/h + mu*stab_c1/pow(h, 2));
const double crhs22 =             2.0*N[0]*crhs1*crhs21*omega[2];
const double crhs23 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs24 =             crhs3*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs25 =             omega[2]*(-gauss_pt_coord[1]*omega[2] + gauss_pt_coord[2]*omega[1]);
const double crhs26 =             -crhs5 + crhs6;
const double crhs27 =             crhs26*omega[0];
const double crhs28 =             crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs29 =             crhs1*(crhs11*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs12*crhs15);
const double crhs30 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs1*(crhs25 - crhs27) - crhs23 + crhs24 + crhs28 + crhs29;
const double crhs31 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs32 =             N[0]*crhs1*crhs31;
const double crhs33 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs1*(crhs26*omega[1] - crhs7) + crhs13 - crhs2 - crhs4 + crhs9;
const double crhs34 =             1.0*crhs21*crhs33;
const double crhs35 =             crhs1*(DN(0,0)*crhs11 + DN(0,1)*crhs12);
const double crhs36 =             crhs1*(-crhs25 + crhs27);
const double crhs37 =             1.0*crhs21*crhs30;
const double crhs38 =             2.0*N[1]*crhs1*crhs21*omega[2];
const double crhs39 =             N[1]*crhs1*crhs31;
const double crhs40 =             crhs1*(DN(1,0)*crhs11 + DN(1,1)*crhs12);
const double crhs41 =             2.0*N[2]*crhs1*crhs21*omega[2];
const double crhs42 =             N[2]*crhs1*crhs31;
const double crhs43 =             crhs1*(DN(2,0)*crhs11 + DN(2,1)*crhs12);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs20 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs13 + N[0]*crhs2 + N[0]*crhs4 + N[0]*crhs8 - N[0]*crhs9 - crhs22*crhs30 - crhs32*crhs34 - crhs34*crhs35;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs20 - DN(0,1)*stress[1] + N[0]*crhs23 - N[0]*crhs24 - N[0]*crhs28 - N[0]*crhs29 + N[0]*crhs36 + crhs22*crhs33 - crhs32*crhs37 - crhs35*crhs37;
            rhs[2]=-DN(0,0)*crhs34 - DN(0,1)*crhs37 - N[0]*crhs16 - N[0]*crhs18 - N[0]*crhs19;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs20 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs13 + N[1]*crhs2 + N[1]*crhs4 + N[1]*crhs8 - N[1]*crhs9 - crhs30*crhs38 - crhs34*crhs39 - crhs34*crhs40;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs20 - DN(1,1)*stress[1] + N[1]*crhs23 - N[1]*crhs24 - N[1]*crhs28 - N[1]*crhs29 + N[1]*crhs36 + crhs33*crhs38 - crhs37*crhs39 - crhs37*crhs40;
            rhs[5]=-DN(1,0)*crhs34 - DN(1,1)*crhs37 - N[1]*crhs16 - N[1]*crhs18 - N[1]*crhs19;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs20 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs13 + N[2]*crhs2 + N[2]*crhs4 + N[2]*crhs8 - N[2]*crhs9 - crhs30*crhs41 - crhs34*crhs42 - crhs34*crhs43;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs20 - DN(2,1)*stress[1] + N[2]*crhs23 - N[2]*crhs24 - N[2]*crhs28 - N[2]*crhs29 + N[2]*crhs36 + crhs33*crhs41 - crhs37*crhs42 - crhs37*crhs43;
            rhs[8]=-DN(2,0)*crhs34 - DN(2,1)*crhs37 - N[2]*crhs16 - N[2]*crhs18 - N[2]*crhs19;


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

    // Non-inertial frame or reference data
    const array_1d<double,3>& omega = rData.AngularVelocity;
    const auto& r_geom = this->GetGeometry();
    array_1d<double,3> gauss_pt_coord = ZeroVector(3);
    for (IndexType i = 0; i < 4; ++i) {
        const auto& r_coords = r_geom[i].Coordinates();
        for (IndexType d = 0; d < 3; ++d) {
            gauss_pt_coord(d) += N(i)*r_coords(d);
        }
    }

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 1.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crhs2 =             crhs1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs3 =             N[0]*crhs1;
const double crhs4 =             gauss_pt_coord[0]*omega[1];
const double crhs5 =             gauss_pt_coord[1]*omega[0];
const double crhs6 =             gauss_pt_coord[0]*omega[2] - gauss_pt_coord[2]*omega[0];
const double crhs7 =             crhs6*omega[2];
const double crhs8 =             crhs7 + omega[1]*(crhs4 - crhs5);
const double crhs9 =             2.0*N[0]*crhs1;
const double crhs10 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs11 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs12 =             crhs10*omega[1] - crhs11*omega[2];
const double crhs13 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
const double crhs14 =             DN(0,0)*v(0,0);
const double crhs15 =             DN(1,0)*v(1,0);
const double crhs16 =             DN(2,0)*v(2,0);
const double crhs17 =             DN(3,0)*v(3,0);
const double crhs18 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs19 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs20 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs21 =             crhs18*(crhs14 + crhs15 + crhs16 + crhs17) + crhs19*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs20*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
const double crhs22 =             crhs1*stab_c2*sqrt(pow(crhs18, 2) + pow(crhs19, 2) + pow(crhs20, 2));
const double crhs23 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs24 =             DN(0,1)*v(0,1);
const double crhs25 =             DN(1,1)*v(1,1);
const double crhs26 =             DN(2,1)*v(2,1);
const double crhs27 =             DN(3,1)*v(3,1);
const double crhs28 =             crhs14 + crhs15 + crhs16 + crhs17 + crhs23 + crhs24 + crhs25 + crhs26 + crhs27;
const double crhs29 =             1.0/crhs1;
const double crhs30 =             crhs29*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crhs31 =             crhs29*(crhs18*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crhs19*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]) + crhs20*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]));
const double crhs32 =             (crhs22*h/stab_c1 + mu)*(crhs28 + crhs30 + crhs31);
const double crhs33 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs34 =             N[0]*crhs1*crhs33;
const double crhs35 =             1.0/(crhs1*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2) + pow(omega[2], 2)) + crhs1*dyn_tau/dt + crhs22/h + mu*stab_c1/pow(h, 2));
const double crhs36 =             DN(0,0)*p[0];
const double crhs37 =             DN(1,0)*p[1];
const double crhs38 =             DN(2,0)*p[2];
const double crhs39 =             DN(3,0)*p[3];
const double crhs40 =             -crhs4 + crhs5;
const double crhs41 =             crhs1*(crhs40*omega[1] - crhs7);
const double crhs42 =             2.0*N[0]*rho[0] + 2.0*N[1]*rho[1] + 2.0*N[2]*rho[2] + 2.0*N[3]*rho[3];
const double crhs43 =             crhs12*crhs42;
const double crhs44 =             crhs1*crhs13;
const double crhs45 =             crhs1*crhs21;
const double crhs46 =             -crhs2 + crhs36 + crhs37 + crhs38 + crhs39 + crhs41 + crhs43 + crhs44 + crhs45;
const double crhs47 =             1.0*crhs35*crhs46;
const double crhs48 =             crhs1*(DN(0,0)*crhs18 + DN(0,1)*crhs19 + DN(0,2)*crhs20);
const double crhs49 =             2.0*N[0]*crhs1*crhs35;
const double crhs50 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs51 =             gauss_pt_coord[2]*omega[1];
const double crhs52 =             gauss_pt_coord[1]*omega[2];
const double crhs53 =             crhs51 - crhs52;
const double crhs54 =             crhs53*omega[2];
const double crhs55 =             crhs40*omega[0];
const double crhs56 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs57 =             crhs56*omega[2];
const double crhs58 =             crhs10*omega[0];
const double crhs59 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
const double crhs60 =             crhs18*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs19*(crhs24 + crhs25 + crhs26 + crhs27) + crhs20*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
const double crhs61 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs1*crhs59 + crhs1*crhs60 + crhs1*(crhs54 - crhs55) + crhs42*(crhs57 - crhs58) - crhs50;
const double crhs62 =             DN(0,2)*p[0];
const double crhs63 =             DN(1,2)*p[1];
const double crhs64 =             DN(2,2)*p[2];
const double crhs65 =             DN(3,2)*p[3];
const double crhs66 =             crhs1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs67 =             crhs6*omega[0];
const double crhs68 =             crhs1*(-crhs53*omega[1] + crhs67);
const double crhs69 =             crhs11*omega[0] - crhs56*omega[1];
const double crhs70 =             crhs42*crhs69;
const double crhs71 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
const double crhs72 =             crhs1*crhs71;
const double crhs73 =             crhs18*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs19*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs20*crhs23;
const double crhs74 =             crhs1*crhs73;
const double crhs75 =             crhs62 + crhs63 + crhs64 + crhs65 - crhs66 + crhs68 + crhs70 + crhs72 + crhs74;
const double crhs76 =             crhs61*omega[2] - crhs75*omega[1];
const double crhs77 =             -crhs54 + crhs55;
const double crhs78 =             -crhs57 + crhs58;
const double crhs79 =             1.0*crhs35*crhs61;
const double crhs80 =             omega[0]*(-crhs62 - crhs63 - crhs64 - crhs65 + crhs66 - crhs68 - crhs70 - crhs72 - crhs74) - omega[2]*(crhs2 - crhs36 - crhs37 - crhs38 - crhs39 - crhs41 - crhs43 - crhs44 - crhs45);
const double crhs81 =             crhs67 + omega[1]*(-crhs51 + crhs52);
const double crhs82 =             1.0*crhs35*crhs75;
const double crhs83 =             crhs46*omega[1] - crhs61*omega[0];
const double crhs84 =             N[1]*crhs1;
const double crhs85 =             2.0*N[1]*crhs1;
const double crhs86 =             N[1]*crhs1*crhs33;
const double crhs87 =             crhs1*(DN(1,0)*crhs18 + DN(1,1)*crhs19 + DN(1,2)*crhs20);
const double crhs88 =             2.0*N[1]*crhs1*crhs35;
const double crhs89 =             N[2]*crhs1;
const double crhs90 =             2.0*N[2]*crhs1;
const double crhs91 =             N[2]*crhs1*crhs33;
const double crhs92 =             crhs1*(DN(2,0)*crhs18 + DN(2,1)*crhs19 + DN(2,2)*crhs20);
const double crhs93 =             2.0*N[2]*crhs1*crhs35;
const double crhs94 =             N[3]*crhs1;
const double crhs95 =             2.0*N[3]*crhs1;
const double crhs96 =             N[3]*crhs1*crhs33;
const double crhs97 =             crhs1*(DN(3,0)*crhs18 + DN(3,1)*crhs19 + DN(3,2)*crhs20);
const double crhs98 =             2.0*N[3]*crhs1*crhs35;
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs32 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs2 - crhs12*crhs9 - crhs13*crhs3 - crhs21*crhs3 + crhs3*crhs8 - crhs34*crhs47 - crhs47*crhs48 - crhs49*crhs76;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs32 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs50 - crhs3*crhs59 - crhs3*crhs60 + crhs3*crhs77 - crhs34*crhs79 - crhs48*crhs79 + crhs49*crhs80 + crhs78*crhs9;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs32 - DN(0,2)*stress[2] + N[0]*crhs66 - crhs3*crhs71 - crhs3*crhs73 - crhs3*crhs81 - crhs34*crhs82 - crhs48*crhs82 - crhs49*crhs83 - crhs69*crhs9;
            rhs[3]=-DN(0,0)*crhs47 - DN(0,1)*crhs79 - DN(0,2)*crhs82 - N[0]*crhs28 - N[0]*crhs30 - N[0]*crhs31;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs32 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs2 - crhs12*crhs85 - crhs13*crhs84 - crhs21*crhs84 - crhs47*crhs86 - crhs47*crhs87 - crhs76*crhs88 + crhs8*crhs84;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs32 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs50 - crhs59*crhs84 - crhs60*crhs84 + crhs77*crhs84 + crhs78*crhs85 - crhs79*crhs86 - crhs79*crhs87 + crhs80*crhs88;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs32 - DN(1,2)*stress[2] + N[1]*crhs66 - crhs69*crhs85 - crhs71*crhs84 - crhs73*crhs84 - crhs81*crhs84 - crhs82*crhs86 - crhs82*crhs87 - crhs83*crhs88;
            rhs[7]=-DN(1,0)*crhs47 - DN(1,1)*crhs79 - DN(1,2)*crhs82 - N[1]*crhs28 - N[1]*crhs30 - N[1]*crhs31;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs32 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs2 - crhs12*crhs90 - crhs13*crhs89 - crhs21*crhs89 - crhs47*crhs91 - crhs47*crhs92 - crhs76*crhs93 + crhs8*crhs89;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs32 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs50 - crhs59*crhs89 - crhs60*crhs89 + crhs77*crhs89 + crhs78*crhs90 - crhs79*crhs91 - crhs79*crhs92 + crhs80*crhs93;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs32 - DN(2,2)*stress[2] + N[2]*crhs66 - crhs69*crhs90 - crhs71*crhs89 - crhs73*crhs89 - crhs81*crhs89 - crhs82*crhs91 - crhs82*crhs92 - crhs83*crhs93;
            rhs[11]=-DN(2,0)*crhs47 - DN(2,1)*crhs79 - DN(2,2)*crhs82 - N[2]*crhs28 - N[2]*crhs30 - N[2]*crhs31;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs32 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs2 - crhs12*crhs95 - crhs13*crhs94 - crhs21*crhs94 - crhs47*crhs96 - crhs47*crhs97 - crhs76*crhs98 + crhs8*crhs94;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs32 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs50 - crhs59*crhs94 - crhs60*crhs94 + crhs77*crhs94 + crhs78*crhs95 - crhs79*crhs96 - crhs79*crhs97 + crhs80*crhs98;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs32 - DN(3,2)*stress[2] + N[3]*crhs66 - crhs69*crhs95 - crhs71*crhs94 - crhs73*crhs94 - crhs81*crhs94 - crhs82*crhs96 - crhs82*crhs97 - crhs83*crhs98;
            rhs[15]=-DN(3,0)*crhs47 - DN(3,1)*crhs79 - DN(3,2)*crhs82 - N[3]*crhs28 - N[3]*crhs30 - N[3]*crhs31;


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