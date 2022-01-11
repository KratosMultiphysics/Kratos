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

    // Non-inertial frame or reference data
    const array_1d<double,3>& omega = rData.AngularVelocity;
    const array_1d<double,3>& omega_old = rData.AngularVelocityOld;

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
const double clhs9 =             DN(0,0)*clhs5 + DN(0,1)*clhs6;
const double clhs10 =             N[0]*clhs4;
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             bdf0*clhs4;
const double clhs13 =             1.0/(clhs4*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2)) + clhs4*dyn_tau/dt + clhs7/h + mu*stab_c1/pow(h, 2));
const double clhs14 =             1.0*clhs13;
const double clhs15 =             clhs14*clhs9;
const double clhs16 =             pow(clhs4, 2);
const double clhs17 =             N[0]*bdf0;
const double clhs18 =             clhs17 + clhs9;
const double clhs19 =             clhs16*clhs18;
const double clhs20 =             clhs13*clhs16;
const double clhs21 =             clhs11*clhs20;
const double clhs22 =             4.0*pow(omega[2], 2);
const double clhs23 =             N[0]*clhs19;
const double clhs24 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs25 =             clhs14*clhs24;
const double clhs26 =             clhs10*clhs9 + clhs11*clhs12 + clhs15*clhs19 + clhs21*clhs22 + clhs23*clhs25;
const double clhs27 =             C(0,1)*DN(0,1) + clhs1;
const double clhs28 =             C(1,2)*DN(0,1);
const double clhs29 =             C(2,2)*DN(0,0) + clhs28;
const double clhs30 =             2.0*omega[2];
const double clhs31 =             clhs30*clhs4;
const double clhs32 =             clhs11*clhs31;
const double clhs33 =             DN(0,0)*clhs8;
const double clhs34 =             DN(0,1)*clhs33;
const double clhs35 =             clhs24*clhs30;
const double clhs36 =             clhs21*clhs35;
const double clhs37 =             clhs20*clhs9;
const double clhs38 =             N[0]*clhs30;
const double clhs39 =             clhs37*clhs38;
const double clhs40 =             clhs13*clhs30;
const double clhs41 =             clhs23*clhs40;
const double clhs42 =             DN(0,0)*N[0];
const double clhs43 =             1/(clhs4*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2));
const double clhs44 =             clhs17*clhs43;
const double clhs45 =             DN(0,1)*N[0];
const double clhs46 =             clhs13*clhs31;
const double clhs47 =             clhs45*clhs46;
const double clhs48 =             clhs25*clhs4;
const double clhs49 =             DN(0,0)*clhs14;
const double clhs50 =             clhs4*clhs9;
const double clhs51 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs52 =             C(0,2)*DN(1,0);
const double clhs53 =             C(2,2)*DN(1,1) + clhs52;
const double clhs54 =             N[1]*clhs4;
const double clhs55 =             clhs20*clhs22;
const double clhs56 =             N[0]*clhs55;
const double clhs57 =             N[1]*clhs56 + clhs17*clhs54;
const double clhs58 =             DN(1,0)*clhs33 + clhs57;
const double clhs59 =             DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs60 =             N[1]*bdf0;
const double clhs61 =             clhs59 + clhs60;
const double clhs62 =             clhs15*clhs16;
const double clhs63 =             clhs16*clhs25;
const double clhs64 =             N[0]*clhs63;
const double clhs65 =             clhs10*clhs59 + clhs61*clhs62 + clhs61*clhs64;
const double clhs66 =             C(0,1)*DN(1,1) + clhs52;
const double clhs67 =             C(1,2)*DN(1,1);
const double clhs68 =             C(2,2)*DN(1,0) + clhs67;
const double clhs69 =             DN(1,1)*clhs33;
const double clhs70 =             clhs20*clhs38;
const double clhs71 =             clhs61*clhs70;
const double clhs72 =             N[1]*clhs30;
const double clhs73 =             clhs37*clhs72;
const double clhs74 =             clhs10*clhs30;
const double clhs75 =             N[1]*clhs74;
const double clhs76 =             N[1]*clhs24;
const double clhs77 =             clhs70*clhs76;
const double clhs78 =             -clhs75 - clhs77;
const double clhs79 =             DN(0,0)*N[1];
const double clhs80 =             clhs43*clhs60;
const double clhs81 =             DN(1,1)*N[0];
const double clhs82 =             clhs46*clhs81;
const double clhs83 =             DN(1,0)*N[0];
const double clhs84 =             DN(1,0)*clhs14;
const double clhs85 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs86 =             C(0,2)*DN(2,0);
const double clhs87 =             C(2,2)*DN(2,1) + clhs86;
const double clhs88 =             N[2]*clhs4;
const double clhs89 =             N[2]*clhs56 + clhs17*clhs88;
const double clhs90 =             DN(2,0)*clhs33 + clhs89;
const double clhs91 =             DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs92 =             N[2]*bdf0;
const double clhs93 =             clhs91 + clhs92;
const double clhs94 =             clhs10*clhs91 + clhs62*clhs93 + clhs64*clhs93;
const double clhs95 =             C(0,1)*DN(2,1) + clhs86;
const double clhs96 =             C(1,2)*DN(2,1);
const double clhs97 =             C(2,2)*DN(2,0) + clhs96;
const double clhs98 =             DN(2,1)*clhs33;
const double clhs99 =             clhs70*clhs93;
const double clhs100 =             N[2]*clhs30;
const double clhs101 =             clhs100*clhs37;
const double clhs102 =             N[2]*clhs74;
const double clhs103 =             clhs100*clhs20;
const double clhs104 =             N[0]*clhs103*clhs24;
const double clhs105 =             -clhs102 - clhs104;
const double clhs106 =             DN(0,0)*N[2];
const double clhs107 =             clhs43*clhs92;
const double clhs108 =             DN(2,1)*N[0];
const double clhs109 =             clhs108*clhs46;
const double clhs110 =             DN(2,0)*N[0];
const double clhs111 =             clhs15*clhs4;
const double clhs112 =             C(0,1)*DN(0,0) + clhs28;
const double clhs113 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs114 =             pow(DN(0,1), 2);
const double clhs115 =             DN(0,1)*clhs8;
const double clhs116 =             -clhs42*clhs46;
const double clhs117 =             DN(0,1)*clhs14;
const double clhs118 =             C(0,1)*DN(1,0) + clhs67;
const double clhs119 =             DN(1,0)*clhs115;
const double clhs120 =             clhs75 + clhs77;
const double clhs121 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs122 =             DN(1,1)*clhs115 + clhs57;
const double clhs123 =             DN(0,1)*N[1];
const double clhs124 =             -clhs46*clhs83;
const double clhs125 =             DN(1,1)*clhs14;
const double clhs126 =             C(0,1)*DN(2,0) + clhs96;
const double clhs127 =             DN(2,0)*clhs115;
const double clhs128 =             clhs102 + clhs104;
const double clhs129 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs130 =             DN(2,1)*clhs115 + clhs89;
const double clhs131 =             DN(0,1)*N[2];
const double clhs132 =             -clhs110*clhs46;
const double clhs133 =             clhs18*clhs4;
const double clhs134 =             bdf0*clhs43;
const double clhs135 =             clhs123*clhs46;
const double clhs136 =             clhs4*clhs61;
const double clhs137 =             -clhs46*clhs79;
const double clhs138 =             DN(1,0)*clhs49 + DN(1,1)*clhs117 + N[1]*clhs44;
const double clhs139 =             clhs131*clhs46;
const double clhs140 =             clhs4*clhs93;
const double clhs141 =             -clhs106*clhs46;
const double clhs142 =             DN(2,0)*clhs49 + DN(2,1)*clhs117 + N[2]*clhs44;
const double clhs143 =             clhs14*clhs59;
const double clhs144 =             N[1]*clhs19;
const double clhs145 =             clhs143*clhs19 + clhs144*clhs25 + clhs54*clhs9;
const double clhs146 =             clhs144*clhs40;
const double clhs147 =             clhs59*clhs70;
const double clhs148 =             DN(1,0)*clhs8;
const double clhs149 =             clhs4*clhs59;
const double clhs150 =             pow(DN(1,0), 2);
const double clhs151 =             pow(N[1], 2);
const double clhs152 =             clhs143*clhs16;
const double clhs153 =             N[1]*clhs63;
const double clhs154 =             clhs12*clhs151 + clhs151*clhs55 + clhs152*clhs61 + clhs153*clhs61 + clhs54*clhs59;
const double clhs155 =             clhs151*clhs31;
const double clhs156 =             DN(1,1)*clhs148;
const double clhs157 =             clhs20*clhs35;
const double clhs158 =             clhs151*clhs157;
const double clhs159 =             clhs20*clhs72;
const double clhs160 =             clhs159*clhs59;
const double clhs161 =             clhs159*clhs61;
const double clhs162 =             DN(1,0)*N[1];
const double clhs163 =             DN(1,1)*N[1];
const double clhs164 =             clhs163*clhs46;
const double clhs165 =             N[1]*N[2]*clhs55 + clhs60*clhs88;
const double clhs166 =             DN(2,0)*clhs148 + clhs165;
const double clhs167 =             clhs152*clhs93 + clhs153*clhs93 + clhs54*clhs91;
const double clhs168 =             DN(2,1)*clhs148;
const double clhs169 =             clhs159*clhs93;
const double clhs170 =             clhs103*clhs59;
const double clhs171 =             clhs100*clhs54;
const double clhs172 =             clhs103*clhs76;
const double clhs173 =             -clhs171 - clhs172;
const double clhs174 =             DN(1,0)*N[2];
const double clhs175 =             DN(2,1)*N[1];
const double clhs176 =             clhs175*clhs46;
const double clhs177 =             DN(2,0)*N[1];
const double clhs178 =             clhs143*clhs4;
const double clhs179 =             DN(1,1)*clhs8;
const double clhs180 =             pow(DN(1,1), 2);
const double clhs181 =             -clhs162*clhs46;
const double clhs182 =             DN(2,0)*clhs179;
const double clhs183 =             clhs171 + clhs172;
const double clhs184 =             DN(2,1)*clhs179 + clhs165;
const double clhs185 =             DN(1,1)*N[2];
const double clhs186 =             -clhs177*clhs46;
const double clhs187 =             clhs185*clhs46;
const double clhs188 =             -clhs174*clhs46;
const double clhs189 =             DN(2,0)*clhs84 + DN(2,1)*clhs125 + N[2]*clhs80;
const double clhs190 =             clhs14*clhs91;
const double clhs191 =             N[2]*clhs19*clhs25 + clhs19*clhs190 + clhs88*clhs9;
const double clhs192 =             clhs100*clhs13*clhs19;
const double clhs193 =             clhs70*clhs91;
const double clhs194 =             DN(2,0)*clhs8;
const double clhs195 =             clhs4*clhs91;
const double clhs196 =             clhs16*clhs190;
const double clhs197 =             N[2]*clhs63;
const double clhs198 =             clhs196*clhs61 + clhs197*clhs61 + clhs59*clhs88;
const double clhs199 =             clhs103*clhs61;
const double clhs200 =             clhs159*clhs91;
const double clhs201 =             pow(DN(2,0), 2);
const double clhs202 =             pow(N[2], 2);
const double clhs203 =             clhs12*clhs202 + clhs196*clhs93 + clhs197*clhs93 + clhs202*clhs55 + clhs88*clhs91;
const double clhs204 =             clhs202*clhs31;
const double clhs205 =             DN(2,1)*clhs194;
const double clhs206 =             clhs157*clhs202;
const double clhs207 =             clhs103*clhs91;
const double clhs208 =             clhs103*clhs93;
const double clhs209 =             DN(2,0)*N[2];
const double clhs210 =             DN(2,1)*N[2];
const double clhs211 =             clhs210*clhs46;
const double clhs212 =             clhs190*clhs4;
const double clhs213 =             DN(2,1)*clhs8;
const double clhs214 =             pow(DN(2,1), 2);
const double clhs215 =             -clhs209*clhs46;
const double clhs216 =             clhs133*clhs14;
const double clhs217 =             clhs136*clhs14;
const double clhs218 =             clhs14*clhs140;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs26 + clhs3*clhs8;
            lhs(0,1)=DN(0,0)*clhs27 + DN(0,1)*clhs29 - clhs32 + clhs34 - clhs36 - clhs39 + clhs41;
            lhs(0,2)=clhs33*clhs44 + clhs42*clhs48 - clhs42 + clhs47 + clhs49*clhs50;
            lhs(0,3)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + clhs58 + clhs65;
            lhs(0,4)=DN(0,0)*clhs66 + DN(0,1)*clhs68 + clhs69 + clhs71 - clhs73 + clhs78;
            lhs(0,5)=clhs33*clhs80 + clhs48*clhs83 + clhs50*clhs84 - clhs79 + clhs82;
            lhs(0,6)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + clhs90 + clhs94;
            lhs(0,7)=DN(0,0)*clhs95 + DN(0,1)*clhs97 - clhs101 + clhs105 + clhs98 + clhs99;
            lhs(0,8)=DN(2,0)*clhs111 - clhs106 + clhs107*clhs33 + clhs109 + clhs110*clhs48;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs112 + clhs32 + clhs34 + clhs36 + clhs39 - clhs41;
            lhs(1,1)=DN(0,0)*clhs29 + DN(0,1)*clhs113 + clhs114*clhs8 + clhs26;
            lhs(1,2)=clhs115*clhs44 + clhs116 + clhs117*clhs50 + clhs45*clhs48 - clhs45;
            lhs(1,3)=DN(0,0)*clhs53 + DN(0,1)*clhs118 + clhs119 + clhs120 - clhs71 + clhs73;
            lhs(1,4)=DN(0,0)*clhs68 + DN(0,1)*clhs121 + clhs122 + clhs65;
            lhs(1,5)=clhs115*clhs80 - clhs123 + clhs124 + clhs125*clhs50 + clhs48*clhs81;
            lhs(1,6)=DN(0,0)*clhs87 + DN(0,1)*clhs126 + clhs101 + clhs127 + clhs128 - clhs99;
            lhs(1,7)=DN(0,0)*clhs97 + DN(0,1)*clhs129 + clhs130 + clhs94;
            lhs(1,8)=DN(2,1)*clhs111 + clhs107*clhs115 + clhs108*clhs48 - clhs131 + clhs132;
            lhs(2,0)=clhs133*clhs49 + clhs42 + clhs47;
            lhs(2,1)=clhs116 + clhs117*clhs133 + clhs45;
            lhs(2,2)=clhs11*clhs134 + clhs114*clhs14 + clhs14*clhs3;
            lhs(2,3)=clhs135 + clhs136*clhs49 + clhs83;
            lhs(2,4)=clhs117*clhs136 + clhs137 + clhs81;
            lhs(2,5)=clhs138;
            lhs(2,6)=clhs110 + clhs139 + clhs140*clhs49;
            lhs(2,7)=clhs108 + clhs117*clhs140 + clhs141;
            lhs(2,8)=clhs142;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs145 + clhs58;
            lhs(3,1)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs119 + clhs146 - clhs147 + clhs78;
            lhs(3,2)=clhs135 + clhs148*clhs44 + clhs149*clhs49 + clhs48*clhs79 - clhs83;
            lhs(3,3)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + clhs150*clhs8 + clhs154;
            lhs(3,4)=DN(1,0)*clhs66 + DN(1,1)*clhs68 - clhs155 + clhs156 - clhs158 - clhs160 + clhs161;
            lhs(3,5)=clhs148*clhs80 + clhs149*clhs84 + clhs162*clhs48 - clhs162 + clhs164;
            lhs(3,6)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + clhs166 + clhs167;
            lhs(3,7)=DN(1,0)*clhs95 + DN(1,1)*clhs97 + clhs168 + clhs169 - clhs170 + clhs173;
            lhs(3,8)=DN(2,0)*clhs178 + clhs107*clhs148 - clhs174 + clhs176 + clhs177*clhs48;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs112 + clhs120 - clhs146 + clhs147 + clhs69;
            lhs(4,1)=DN(1,0)*clhs29 + DN(1,1)*clhs113 + clhs122 + clhs145;
            lhs(4,2)=clhs117*clhs149 + clhs123*clhs48 + clhs137 + clhs179*clhs44 - clhs81;
            lhs(4,3)=DN(1,0)*clhs53 + DN(1,1)*clhs118 + clhs155 + clhs156 + clhs158 + clhs160 - clhs161;
            lhs(4,4)=DN(1,0)*clhs68 + DN(1,1)*clhs121 + clhs154 + clhs180*clhs8;
            lhs(4,5)=clhs125*clhs149 + clhs163*clhs48 - clhs163 + clhs179*clhs80 + clhs181;
            lhs(4,6)=DN(1,0)*clhs87 + DN(1,1)*clhs126 - clhs169 + clhs170 + clhs182 + clhs183;
            lhs(4,7)=DN(1,0)*clhs97 + DN(1,1)*clhs129 + clhs167 + clhs184;
            lhs(4,8)=DN(2,1)*clhs178 + clhs107*clhs179 + clhs175*clhs48 - clhs185 + clhs186;
            lhs(5,0)=clhs133*clhs84 + clhs79 + clhs82;
            lhs(5,1)=clhs123 + clhs124 + clhs125*clhs133;
            lhs(5,2)=clhs138;
            lhs(5,3)=clhs136*clhs84 + clhs162 + clhs164;
            lhs(5,4)=clhs125*clhs136 + clhs163 + clhs181;
            lhs(5,5)=clhs134*clhs151 + clhs14*clhs150 + clhs14*clhs180;
            lhs(5,6)=clhs140*clhs84 + clhs177 + clhs187;
            lhs(5,7)=clhs125*clhs140 + clhs175 + clhs188;
            lhs(5,8)=clhs189;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs191 + clhs90;
            lhs(6,1)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs105 + clhs127 + clhs192 - clhs193;
            lhs(6,2)=clhs106*clhs48 - clhs110 + clhs139 + clhs194*clhs44 + clhs195*clhs49;
            lhs(6,3)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + clhs166 + clhs198;
            lhs(6,4)=DN(2,0)*clhs66 + DN(2,1)*clhs68 + clhs173 + clhs182 + clhs199 - clhs200;
            lhs(6,5)=clhs174*clhs48 - clhs177 + clhs187 + clhs194*clhs80 + clhs195*clhs84;
            lhs(6,6)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + clhs201*clhs8 + clhs203;
            lhs(6,7)=DN(2,0)*clhs95 + DN(2,1)*clhs97 - clhs204 + clhs205 - clhs206 - clhs207 + clhs208;
            lhs(6,8)=DN(2,0)*clhs212 + clhs107*clhs194 + clhs209*clhs48 - clhs209 + clhs211;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs112 + clhs128 - clhs192 + clhs193 + clhs98;
            lhs(7,1)=DN(2,0)*clhs29 + DN(2,1)*clhs113 + clhs130 + clhs191;
            lhs(7,2)=-clhs108 + clhs117*clhs195 + clhs131*clhs48 + clhs141 + clhs213*clhs44;
            lhs(7,3)=DN(2,0)*clhs53 + DN(2,1)*clhs118 + clhs168 + clhs183 - clhs199 + clhs200;
            lhs(7,4)=DN(2,0)*clhs68 + DN(2,1)*clhs121 + clhs184 + clhs198;
            lhs(7,5)=clhs125*clhs195 - clhs175 + clhs185*clhs48 + clhs188 + clhs213*clhs80;
            lhs(7,6)=DN(2,0)*clhs87 + DN(2,1)*clhs126 + clhs204 + clhs205 + clhs206 + clhs207 - clhs208;
            lhs(7,7)=DN(2,0)*clhs97 + DN(2,1)*clhs129 + clhs203 + clhs214*clhs8;
            lhs(7,8)=DN(2,1)*clhs212 + clhs107*clhs213 + clhs210*clhs48 - clhs210 + clhs215;
            lhs(8,0)=DN(2,0)*clhs216 + clhs106 + clhs109;
            lhs(8,1)=DN(2,1)*clhs216 + clhs131 + clhs132;
            lhs(8,2)=clhs142;
            lhs(8,3)=DN(2,0)*clhs217 + clhs174 + clhs176;
            lhs(8,4)=DN(2,1)*clhs217 + clhs185 + clhs186;
            lhs(8,5)=clhs189;
            lhs(8,6)=DN(2,0)*clhs218 + clhs209 + clhs211;
            lhs(8,7)=DN(2,1)*clhs218 + clhs210 + clhs215;
            lhs(8,8)=clhs134*clhs202 + clhs14*clhs201 + clhs14*clhs214;


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
    const array_1d<double,3>& omega_old = rData.AngularVelocityOld;

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
const double clhs16 =             pow(omega[0], 2);
const double clhs17 =             1.0/(clhs10/h + clhs6*stab_c3*sqrt(clhs14 + clhs16) + clhs6*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double clhs18 =             clhs17*pow(clhs6, 2);
const double clhs19 =             4.0*clhs18;
const double clhs20 =             clhs15*clhs19;
const double clhs21 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs22 =             N[0]*clhs6;
const double clhs23 =             bdf0*clhs6;
const double clhs24 =             N[0]*bdf0;
const double clhs25 =             1.0*clhs21 + 1.0*clhs24;
const double clhs26 =             clhs18*clhs21;
const double clhs27 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs28 =             clhs18*clhs27;
const double clhs29 =             N[0]*clhs28;
const double clhs30 =             clhs15*clhs23 + clhs21*clhs22 + clhs25*clhs26 + clhs25*clhs29;
const double clhs31 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs32 =             C(1,3)*DN(0,1);
const double clhs33 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs32;
const double clhs34 =             C(3,5)*DN(0,0);
const double clhs35 =             C(4,5)*DN(0,2);
const double clhs36 =             C(1,5)*DN(0,1) + clhs34 + clhs35;
const double clhs37 =             2.0*omega[2];
const double clhs38 =             clhs15*clhs6;
const double clhs39 =             clhs37*clhs38;
const double clhs40 =             DN(0,0)*clhs11;
const double clhs41 =             DN(0,1)*clhs40;
const double clhs42 =             clhs15*clhs28;
const double clhs43 =             clhs37*clhs42;
const double clhs44 =             N[0]*clhs26;
const double clhs45 =             clhs37*clhs44;
const double clhs46 =             2.0*omega[1];
const double clhs47 =             N[0]*omega[0];
const double clhs48 =             clhs46*clhs47;
const double clhs49 =             clhs25*omega[2];
const double clhs50 =             clhs48 - clhs49;
const double clhs51 =             2.0*clhs18;
const double clhs52 =             N[0]*clhs51;
const double clhs53 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs54 =             C(3,4)*DN(0,1);
const double clhs55 =             C(2,3)*DN(0,2) + clhs34 + clhs54;
const double clhs56 =             C(2,5)*DN(0,2);
const double clhs57 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs56;
const double clhs58 =             clhs38*clhs46;
const double clhs59 =             DN(0,2)*clhs40;
const double clhs60 =             clhs42*clhs46;
const double clhs61 =             clhs44*clhs46;
const double clhs62 =             clhs37*clhs47;
const double clhs63 =             clhs25*omega[1];
const double clhs64 =             clhs62 + clhs63;
const double clhs65 =             DN(0,0)*N[0];
const double clhs66 =             1/(clhs6*pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2));
const double clhs67 =             clhs24*clhs66;
const double clhs68 =             DN(0,1)*omega[2];
const double clhs69 =             DN(0,2)*omega[1];
const double clhs70 =             clhs68 - clhs69;
const double clhs71 =             2.0*clhs22;
const double clhs72 =             clhs17*clhs71;
const double clhs73 =             1.0*clhs17;
const double clhs74 =             clhs27*clhs6*clhs73;
const double clhs75 =             DN(0,0)*clhs73;
const double clhs76 =             clhs21*clhs6;
const double clhs77 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs78 =             C(0,3)*DN(1,0);
const double clhs79 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs78;
const double clhs80 =             C(0,5)*DN(1,0);
const double clhs81 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs80;
const double clhs82 =             N[1]*clhs6;
const double clhs83 =             clhs24*clhs82;
const double clhs84 =             N[0]*clhs19;
const double clhs85 =             clhs14*clhs84;
const double clhs86 =             DN(1,0)*clhs40 + N[1]*clhs85 + clhs83;
const double clhs87 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs88 =             N[1]*bdf0;
const double clhs89 =             1.0*clhs87 + 1.0*clhs88;
const double clhs90 =             clhs22*clhs87 + clhs26*clhs89 + clhs29*clhs89;
const double clhs91 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs78;
const double clhs92 =             C(1,3)*DN(1,1);
const double clhs93 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs92;
const double clhs94 =             C(3,5)*DN(1,0);
const double clhs95 =             C(4,5)*DN(1,2);
const double clhs96 =             C(1,5)*DN(1,1) + clhs94 + clhs95;
const double clhs97 =             DN(1,1)*clhs40;
const double clhs98 =             N[1]*omega[0];
const double clhs99 =             clhs46*clhs98;
const double clhs100 =             clhs89*omega[2];
const double clhs101 =             -clhs100 + clhs99;
const double clhs102 =             N[1]*clhs26;
const double clhs103 =             clhs102*clhs37;
const double clhs104 =             N[1]*clhs22;
const double clhs105 =             clhs104*clhs37;
const double clhs106 =             N[1]*clhs29;
const double clhs107 =             clhs106*clhs37;
const double clhs108 =             -clhs105 - clhs107;
const double clhs109 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs80;
const double clhs110 =             C(3,4)*DN(1,1);
const double clhs111 =             C(2,3)*DN(1,2) + clhs110 + clhs94;
const double clhs112 =             C(2,5)*DN(1,2);
const double clhs113 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs112;
const double clhs114 =             DN(1,2)*clhs40;
const double clhs115 =             clhs37*clhs98;
const double clhs116 =             clhs89*omega[1];
const double clhs117 =             clhs115 + clhs116;
const double clhs118 =             clhs102*clhs46;
const double clhs119 =             clhs104*clhs46;
const double clhs120 =             clhs106*clhs46;
const double clhs121 =             clhs119 + clhs120;
const double clhs122 =             DN(0,0)*N[1];
const double clhs123 =             clhs66*clhs88;
const double clhs124 =             DN(1,1)*omega[2];
const double clhs125 =             DN(1,2)*omega[1];
const double clhs126 =             clhs124 - clhs125;
const double clhs127 =             DN(1,0)*N[0];
const double clhs128 =             DN(1,0)*clhs73;
const double clhs129 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs130 =             C(0,3)*DN(2,0);
const double clhs131 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs130;
const double clhs132 =             C(0,5)*DN(2,0);
const double clhs133 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs132;
const double clhs134 =             N[2]*clhs6;
const double clhs135 =             clhs134*clhs24;
const double clhs136 =             DN(2,0)*clhs40 + N[2]*clhs85 + clhs135;
const double clhs137 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs138 =             N[2]*bdf0;
const double clhs139 =             1.0*clhs137 + 1.0*clhs138;
const double clhs140 =             clhs137*clhs22 + clhs139*clhs26 + clhs139*clhs29;
const double clhs141 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs130;
const double clhs142 =             C(1,3)*DN(2,1);
const double clhs143 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs142;
const double clhs144 =             C(3,5)*DN(2,0);
const double clhs145 =             C(4,5)*DN(2,2);
const double clhs146 =             C(1,5)*DN(2,1) + clhs144 + clhs145;
const double clhs147 =             DN(2,1)*clhs40;
const double clhs148 =             N[2]*omega[0];
const double clhs149 =             clhs148*clhs46;
const double clhs150 =             clhs139*omega[2];
const double clhs151 =             clhs149 - clhs150;
const double clhs152 =             N[2]*clhs26;
const double clhs153 =             clhs152*clhs37;
const double clhs154 =             N[2]*clhs22;
const double clhs155 =             clhs154*clhs37;
const double clhs156 =             N[2]*clhs29;
const double clhs157 =             clhs156*clhs37;
const double clhs158 =             -clhs155 - clhs157;
const double clhs159 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs132;
const double clhs160 =             C(3,4)*DN(2,1);
const double clhs161 =             C(2,3)*DN(2,2) + clhs144 + clhs160;
const double clhs162 =             C(2,5)*DN(2,2);
const double clhs163 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs162;
const double clhs164 =             DN(2,2)*clhs40;
const double clhs165 =             clhs148*clhs37;
const double clhs166 =             clhs139*omega[1];
const double clhs167 =             clhs165 + clhs166;
const double clhs168 =             clhs152*clhs46;
const double clhs169 =             clhs154*clhs46;
const double clhs170 =             clhs156*clhs46;
const double clhs171 =             clhs169 + clhs170;
const double clhs172 =             DN(0,0)*N[2];
const double clhs173 =             clhs138*clhs66;
const double clhs174 =             DN(2,1)*omega[2];
const double clhs175 =             DN(2,2)*omega[1];
const double clhs176 =             clhs174 - clhs175;
const double clhs177 =             DN(2,0)*N[0];
const double clhs178 =             DN(2,0)*clhs73;
const double clhs179 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs180 =             C(0,3)*DN(3,0);
const double clhs181 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs180;
const double clhs182 =             C(0,5)*DN(3,0);
const double clhs183 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs182;
const double clhs184 =             N[3]*clhs6;
const double clhs185 =             clhs184*clhs24;
const double clhs186 =             DN(3,0)*clhs40 + N[3]*clhs85 + clhs185;
const double clhs187 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs188 =             N[3]*bdf0;
const double clhs189 =             1.0*clhs187 + 1.0*clhs188;
const double clhs190 =             clhs187*clhs22 + clhs189*clhs26 + clhs189*clhs29;
const double clhs191 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs180;
const double clhs192 =             C(1,3)*DN(3,1);
const double clhs193 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs192;
const double clhs194 =             C(3,5)*DN(3,0);
const double clhs195 =             C(4,5)*DN(3,2);
const double clhs196 =             C(1,5)*DN(3,1) + clhs194 + clhs195;
const double clhs197 =             DN(3,1)*clhs40;
const double clhs198 =             N[3]*omega[0];
const double clhs199 =             clhs198*clhs46;
const double clhs200 =             clhs189*omega[2];
const double clhs201 =             clhs199 - clhs200;
const double clhs202 =             N[3]*clhs26;
const double clhs203 =             clhs202*clhs37;
const double clhs204 =             N[3]*clhs22;
const double clhs205 =             clhs204*clhs37;
const double clhs206 =             N[3]*clhs29;
const double clhs207 =             clhs206*clhs37;
const double clhs208 =             -clhs205 - clhs207;
const double clhs209 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs182;
const double clhs210 =             C(3,4)*DN(3,1);
const double clhs211 =             C(2,3)*DN(3,2) + clhs194 + clhs210;
const double clhs212 =             C(2,5)*DN(3,2);
const double clhs213 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs212;
const double clhs214 =             DN(3,2)*clhs40;
const double clhs215 =             clhs198*clhs37;
const double clhs216 =             clhs189*omega[1];
const double clhs217 =             clhs215 + clhs216;
const double clhs218 =             clhs202*clhs46;
const double clhs219 =             clhs204*clhs46;
const double clhs220 =             clhs206*clhs46;
const double clhs221 =             clhs219 + clhs220;
const double clhs222 =             DN(0,0)*N[3];
const double clhs223 =             clhs188*clhs66;
const double clhs224 =             DN(3,1)*omega[2];
const double clhs225 =             DN(3,2)*omega[1];
const double clhs226 =             clhs224 - clhs225;
const double clhs227 =             DN(3,0)*N[0];
const double clhs228 =             clhs73*clhs76;
const double clhs229 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs32;
const double clhs230 =             C(0,4)*DN(0,0) + clhs35 + clhs54;
const double clhs231 =             clhs48 + clhs49;
const double clhs232 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs233 =             C(1,4)*DN(0,1);
const double clhs234 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs233;
const double clhs235 =             pow(DN(0,1), 2);
const double clhs236 =             clhs13 + clhs16;
const double clhs237 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs233;
const double clhs238 =             C(2,4)*DN(0,2);
const double clhs239 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs238;
const double clhs240 =             2.0*omega[0];
const double clhs241 =             clhs240*clhs38;
const double clhs242 =             DN(0,1)*clhs11;
const double clhs243 =             DN(0,2)*clhs242;
const double clhs244 =             clhs240*clhs42;
const double clhs245 =             2.0*clhs26;
const double clhs246 =             clhs245*clhs47;
const double clhs247 =             clhs37*omega[1];
const double clhs248 =             N[0]*clhs247;
const double clhs249 =             clhs25*omega[0];
const double clhs250 =             clhs248 - clhs249;
const double clhs251 =             DN(0,1)*N[0];
const double clhs252 =             DN(0,2)*omega[0];
const double clhs253 =             DN(0,0)*omega[2] - clhs252;
const double clhs254 =             DN(0,1)*clhs73;
const double clhs255 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs92;
const double clhs256 =             C(0,4)*DN(1,0) + clhs110 + clhs95;
const double clhs257 =             DN(1,0)*clhs242;
const double clhs258 =             clhs100 + clhs99;
const double clhs259 =             clhs105 + clhs107;
const double clhs260 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs261 =             C(1,4)*DN(1,1);
const double clhs262 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs261;
const double clhs263 =             clhs83 + clhs90;
const double clhs264 =             clhs236*clhs84;
const double clhs265 =             DN(1,1)*clhs242 + N[1]*clhs264;
const double clhs266 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs261;
const double clhs267 =             C(2,4)*DN(1,2);
const double clhs268 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs267;
const double clhs269 =             DN(1,2)*clhs242;
const double clhs270 =             N[1]*clhs247;
const double clhs271 =             clhs89*omega[0];
const double clhs272 =             clhs270 - clhs271;
const double clhs273 =             clhs245*clhs98;
const double clhs274 =             clhs71*clhs98;
const double clhs275 =             N[1]*clhs51;
const double clhs276 =             clhs27*clhs47;
const double clhs277 =             clhs275*clhs276;
const double clhs278 =             -clhs274 - clhs277;
const double clhs279 =             DN(0,1)*N[1];
const double clhs280 =             DN(1,0)*omega[2];
const double clhs281 =             DN(1,2)*omega[0];
const double clhs282 =             clhs280 - clhs281;
const double clhs283 =             DN(1,1)*N[0];
const double clhs284 =             DN(1,1)*clhs73;
const double clhs285 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs142;
const double clhs286 =             C(0,4)*DN(2,0) + clhs145 + clhs160;
const double clhs287 =             DN(2,0)*clhs242;
const double clhs288 =             clhs149 + clhs150;
const double clhs289 =             clhs155 + clhs157;
const double clhs290 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs291 =             C(1,4)*DN(2,1);
const double clhs292 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs291;
const double clhs293 =             clhs135 + clhs140;
const double clhs294 =             DN(2,1)*clhs242 + N[2]*clhs264;
const double clhs295 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs291;
const double clhs296 =             C(2,4)*DN(2,2);
const double clhs297 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs296;
const double clhs298 =             DN(2,2)*clhs242;
const double clhs299 =             N[2]*clhs247;
const double clhs300 =             clhs139*omega[0];
const double clhs301 =             clhs299 - clhs300;
const double clhs302 =             clhs148*clhs245;
const double clhs303 =             clhs148*clhs71;
const double clhs304 =             N[2]*clhs51;
const double clhs305 =             clhs276*clhs304;
const double clhs306 =             -clhs303 - clhs305;
const double clhs307 =             DN(0,1)*N[2];
const double clhs308 =             DN(2,0)*omega[2];
const double clhs309 =             DN(2,2)*omega[0];
const double clhs310 =             clhs308 - clhs309;
const double clhs311 =             DN(2,1)*N[0];
const double clhs312 =             DN(2,1)*clhs73;
const double clhs313 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs192;
const double clhs314 =             C(0,4)*DN(3,0) + clhs195 + clhs210;
const double clhs315 =             DN(3,0)*clhs242;
const double clhs316 =             clhs199 + clhs200;
const double clhs317 =             clhs205 + clhs207;
const double clhs318 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs319 =             C(1,4)*DN(3,1);
const double clhs320 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs319;
const double clhs321 =             clhs185 + clhs190;
const double clhs322 =             DN(3,1)*clhs242 + N[3]*clhs264;
const double clhs323 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs319;
const double clhs324 =             C(2,4)*DN(3,2);
const double clhs325 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs324;
const double clhs326 =             DN(3,2)*clhs242;
const double clhs327 =             N[3]*clhs247;
const double clhs328 =             clhs189*omega[0];
const double clhs329 =             clhs327 - clhs328;
const double clhs330 =             clhs198*clhs245;
const double clhs331 =             clhs198*clhs71;
const double clhs332 =             N[3]*clhs51;
const double clhs333 =             clhs276*clhs332;
const double clhs334 =             -clhs331 - clhs333;
const double clhs335 =             DN(0,1)*N[3];
const double clhs336 =             DN(3,0)*omega[2];
const double clhs337 =             DN(3,2)*omega[0];
const double clhs338 =             clhs336 - clhs337;
const double clhs339 =             DN(3,1)*N[0];
const double clhs340 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs56;
const double clhs341 =             clhs62 - clhs63;
const double clhs342 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs238;
const double clhs343 =             clhs248 + clhs249;
const double clhs344 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs345 =             pow(DN(0,2), 2);
const double clhs346 =             clhs12 + clhs16;
const double clhs347 =             DN(0,2)*N[0];
const double clhs348 =             DN(0,2)*clhs11;
const double clhs349 =             DN(0,0)*omega[1] - DN(0,1)*omega[0];
const double clhs350 =             DN(0,2)*clhs73;
const double clhs351 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs112;
const double clhs352 =             DN(1,0)*clhs348;
const double clhs353 =             clhs115 - clhs116;
const double clhs354 =             -clhs119 - clhs120;
const double clhs355 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs267;
const double clhs356 =             DN(1,1)*clhs348;
const double clhs357 =             clhs270 + clhs271;
const double clhs358 =             clhs274 + clhs277;
const double clhs359 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs360 =             clhs346*clhs84;
const double clhs361 =             DN(1,2)*clhs348 + N[1]*clhs360;
const double clhs362 =             DN(0,2)*N[1];
const double clhs363 =             DN(1,0)*omega[1];
const double clhs364 =             DN(1,1)*omega[0];
const double clhs365 =             clhs363 - clhs364;
const double clhs366 =             DN(1,2)*N[0];
const double clhs367 =             DN(1,2)*clhs73;
const double clhs368 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs162;
const double clhs369 =             DN(2,0)*clhs348;
const double clhs370 =             clhs165 - clhs166;
const double clhs371 =             -clhs169 - clhs170;
const double clhs372 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs296;
const double clhs373 =             DN(2,1)*clhs348;
const double clhs374 =             clhs299 + clhs300;
const double clhs375 =             clhs303 + clhs305;
const double clhs376 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs377 =             DN(2,2)*clhs348 + N[2]*clhs360;
const double clhs378 =             DN(0,2)*N[2];
const double clhs379 =             DN(2,0)*omega[1];
const double clhs380 =             DN(2,1)*omega[0];
const double clhs381 =             clhs379 - clhs380;
const double clhs382 =             DN(2,2)*N[0];
const double clhs383 =             DN(2,2)*clhs73;
const double clhs384 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs212;
const double clhs385 =             DN(3,0)*clhs348;
const double clhs386 =             clhs215 - clhs216;
const double clhs387 =             -clhs219 - clhs220;
const double clhs388 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs324;
const double clhs389 =             DN(3,1)*clhs348;
const double clhs390 =             clhs327 + clhs328;
const double clhs391 =             clhs331 + clhs333;
const double clhs392 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs393 =             DN(3,2)*clhs348 + N[3]*clhs360;
const double clhs394 =             DN(0,2)*N[3];
const double clhs395 =             DN(3,0)*omega[1];
const double clhs396 =             DN(3,1)*omega[0];
const double clhs397 =             clhs395 - clhs396;
const double clhs398 =             DN(3,2)*N[0];
const double clhs399 =             clhs17*clhs6;
const double clhs400 =             clhs25*clhs399;
const double clhs401 =             clhs399*clhs65;
const double clhs402 =             clhs240*clhs399;
const double clhs403 =             bdf0*clhs66;
const double clhs404 =             2.0*clhs82;
const double clhs405 =             clhs17*clhs404;
const double clhs406 =             clhs399*clhs89;
const double clhs407 =             clhs122*clhs399;
const double clhs408 =             DN(1,0)*clhs75 + DN(1,1)*clhs254 + DN(1,2)*clhs350 + N[1]*clhs67;
const double clhs409 =             2.0*clhs134;
const double clhs410 =             clhs17*clhs409;
const double clhs411 =             clhs139*clhs399;
const double clhs412 =             clhs172*clhs399;
const double clhs413 =             DN(2,0)*clhs75 + DN(2,1)*clhs254 + DN(2,2)*clhs350 + N[2]*clhs67;
const double clhs414 =             2.0*clhs17*clhs184;
const double clhs415 =             clhs189*clhs399;
const double clhs416 =             clhs222*clhs399;
const double clhs417 =             DN(3,0)*clhs75 + DN(3,1)*clhs254 + DN(3,2)*clhs350 + N[3]*clhs67;
const double clhs418 =             clhs18*clhs87;
const double clhs419 =             N[1]*clhs28;
const double clhs420 =             clhs21*clhs82 + clhs25*clhs418 + clhs25*clhs419;
const double clhs421 =             N[0]*clhs418;
const double clhs422 =             clhs37*clhs421;
const double clhs423 =             clhs421*clhs46;
const double clhs424 =             DN(1,0)*clhs11;
const double clhs425 =             clhs6*clhs87;
const double clhs426 =             pow(DN(1,0), 2);
const double clhs427 =             pow(N[1], 2);
const double clhs428 =             clhs19*clhs427;
const double clhs429 =             clhs23*clhs427 + clhs418*clhs89 + clhs419*clhs89 + clhs82*clhs87;
const double clhs430 =             clhs427*clhs6;
const double clhs431 =             clhs37*clhs430;
const double clhs432 =             DN(1,1)*clhs424;
const double clhs433 =             clhs28*clhs427;
const double clhs434 =             clhs37*clhs433;
const double clhs435 =             N[1]*clhs418;
const double clhs436 =             clhs37*clhs435;
const double clhs437 =             clhs430*clhs46;
const double clhs438 =             DN(1,2)*clhs424;
const double clhs439 =             clhs433*clhs46;
const double clhs440 =             clhs435*clhs46;
const double clhs441 =             DN(1,0)*N[1];
const double clhs442 =             clhs134*clhs88;
const double clhs443 =             N[1]*clhs19;
const double clhs444 =             clhs14*clhs443;
const double clhs445 =             DN(2,0)*clhs424 + N[2]*clhs444 + clhs442;
const double clhs446 =             clhs137*clhs82 + clhs139*clhs418 + clhs139*clhs419;
const double clhs447 =             DN(2,1)*clhs424;
const double clhs448 =             N[2]*clhs418;
const double clhs449 =             clhs37*clhs448;
const double clhs450 =             N[2]*clhs82;
const double clhs451 =             clhs37*clhs450;
const double clhs452 =             N[2]*clhs419;
const double clhs453 =             clhs37*clhs452;
const double clhs454 =             -clhs451 - clhs453;
const double clhs455 =             DN(2,2)*clhs424;
const double clhs456 =             clhs448*clhs46;
const double clhs457 =             clhs450*clhs46;
const double clhs458 =             clhs452*clhs46;
const double clhs459 =             clhs457 + clhs458;
const double clhs460 =             DN(1,0)*N[2];
const double clhs461 =             DN(2,0)*N[1];
const double clhs462 =             clhs184*clhs88;
const double clhs463 =             DN(3,0)*clhs424 + N[3]*clhs444 + clhs462;
const double clhs464 =             clhs187*clhs82 + clhs189*clhs418 + clhs189*clhs419;
const double clhs465 =             DN(3,1)*clhs424;
const double clhs466 =             N[3]*clhs418;
const double clhs467 =             clhs37*clhs466;
const double clhs468 =             N[3]*clhs82;
const double clhs469 =             clhs37*clhs468;
const double clhs470 =             N[3]*clhs419;
const double clhs471 =             clhs37*clhs470;
const double clhs472 =             -clhs469 - clhs471;
const double clhs473 =             DN(3,2)*clhs424;
const double clhs474 =             clhs46*clhs466;
const double clhs475 =             clhs46*clhs468;
const double clhs476 =             clhs46*clhs470;
const double clhs477 =             clhs475 + clhs476;
const double clhs478 =             DN(1,0)*N[3];
const double clhs479 =             DN(3,0)*N[1];
const double clhs480 =             clhs425*clhs73;
const double clhs481 =             clhs420 + clhs83;
const double clhs482 =             clhs51*clhs87;
const double clhs483 =             clhs47*clhs482;
const double clhs484 =             DN(1,1)*clhs11;
const double clhs485 =             pow(DN(1,1), 2);
const double clhs486 =             clhs240*clhs430;
const double clhs487 =             DN(1,2)*clhs484;
const double clhs488 =             clhs240*clhs433;
const double clhs489 =             clhs482*clhs98;
const double clhs490 =             DN(1,1)*N[1];
const double clhs491 =             DN(2,0)*clhs484;
const double clhs492 =             clhs451 + clhs453;
const double clhs493 =             clhs442 + clhs446;
const double clhs494 =             clhs236*clhs443;
const double clhs495 =             DN(2,1)*clhs484 + N[2]*clhs494;
const double clhs496 =             DN(2,2)*clhs484;
const double clhs497 =             clhs148*clhs482;
const double clhs498 =             clhs148*clhs404;
const double clhs499 =             clhs27*clhs98;
const double clhs500 =             clhs304*clhs499;
const double clhs501 =             -clhs498 - clhs500;
const double clhs502 =             DN(1,1)*N[2];
const double clhs503 =             DN(2,1)*N[1];
const double clhs504 =             DN(3,0)*clhs484;
const double clhs505 =             clhs469 + clhs471;
const double clhs506 =             clhs462 + clhs464;
const double clhs507 =             DN(3,1)*clhs484 + N[3]*clhs494;
const double clhs508 =             DN(3,2)*clhs484;
const double clhs509 =             clhs198*clhs482;
const double clhs510 =             clhs198*clhs404;
const double clhs511 =             clhs332*clhs499;
const double clhs512 =             -clhs510 - clhs511;
const double clhs513 =             DN(1,1)*N[3];
const double clhs514 =             DN(3,1)*N[1];
const double clhs515 =             DN(1,2)*clhs11;
const double clhs516 =             pow(DN(1,2), 2);
const double clhs517 =             DN(1,2)*N[1];
const double clhs518 =             DN(2,0)*clhs515;
const double clhs519 =             -clhs457 - clhs458;
const double clhs520 =             DN(2,1)*clhs515;
const double clhs521 =             clhs498 + clhs500;
const double clhs522 =             clhs346*clhs443;
const double clhs523 =             DN(2,2)*clhs515 + N[2]*clhs522;
const double clhs524 =             DN(1,2)*N[2];
const double clhs525 =             DN(2,2)*N[1];
const double clhs526 =             DN(3,0)*clhs515;
const double clhs527 =             -clhs475 - clhs476;
const double clhs528 =             DN(3,1)*clhs515;
const double clhs529 =             clhs510 + clhs511;
const double clhs530 =             DN(3,2)*clhs515 + N[3]*clhs522;
const double clhs531 =             DN(1,2)*N[3];
const double clhs532 =             DN(3,2)*N[1];
const double clhs533 =             DN(2,0)*clhs128 + DN(2,1)*clhs284 + DN(2,2)*clhs367 + N[2]*clhs123;
const double clhs534 =             DN(3,0)*clhs128 + DN(3,1)*clhs284 + DN(3,2)*clhs367 + N[3]*clhs123;
const double clhs535 =             clhs137*clhs18;
const double clhs536 =             N[2]*clhs28;
const double clhs537 =             clhs134*clhs21 + clhs25*clhs535 + clhs25*clhs536;
const double clhs538 =             N[0]*clhs535;
const double clhs539 =             clhs37*clhs538;
const double clhs540 =             clhs46*clhs538;
const double clhs541 =             DN(2,0)*clhs11;
const double clhs542 =             clhs137*clhs6;
const double clhs543 =             clhs134*clhs87 + clhs535*clhs89 + clhs536*clhs89;
const double clhs544 =             N[1]*clhs535;
const double clhs545 =             clhs37*clhs544;
const double clhs546 =             clhs46*clhs544;
const double clhs547 =             pow(DN(2,0), 2);
const double clhs548 =             pow(N[2], 2);
const double clhs549 =             clhs19*clhs548;
const double clhs550 =             clhs134*clhs137 + clhs139*clhs535 + clhs139*clhs536 + clhs23*clhs548;
const double clhs551 =             clhs548*clhs6;
const double clhs552 =             clhs37*clhs551;
const double clhs553 =             DN(2,1)*clhs541;
const double clhs554 =             clhs28*clhs548;
const double clhs555 =             clhs37*clhs554;
const double clhs556 =             N[2]*clhs535;
const double clhs557 =             clhs37*clhs556;
const double clhs558 =             clhs46*clhs551;
const double clhs559 =             DN(2,2)*clhs541;
const double clhs560 =             clhs46*clhs554;
const double clhs561 =             clhs46*clhs556;
const double clhs562 =             DN(2,0)*N[2];
const double clhs563 =             clhs138*clhs184;
const double clhs564 =             N[2]*N[3]*clhs19;
const double clhs565 =             DN(3,0)*clhs541 + clhs14*clhs564 + clhs563;
const double clhs566 =             clhs134*clhs187 + clhs189*clhs535 + clhs189*clhs536;
const double clhs567 =             DN(3,1)*clhs541;
const double clhs568 =             N[3]*clhs535;
const double clhs569 =             clhs37*clhs568;
const double clhs570 =             N[3]*clhs134;
const double clhs571 =             clhs37*clhs570;
const double clhs572 =             N[3]*clhs536;
const double clhs573 =             clhs37*clhs572;
const double clhs574 =             -clhs571 - clhs573;
const double clhs575 =             DN(3,2)*clhs541;
const double clhs576 =             clhs46*clhs568;
const double clhs577 =             clhs46*clhs570;
const double clhs578 =             clhs46*clhs572;
const double clhs579 =             clhs577 + clhs578;
const double clhs580 =             DN(2,0)*N[3];
const double clhs581 =             DN(3,0)*N[2];
const double clhs582 =             clhs542*clhs73;
const double clhs583 =             clhs135 + clhs537;
const double clhs584 =             clhs137*clhs51;
const double clhs585 =             clhs47*clhs584;
const double clhs586 =             DN(2,1)*clhs11;
const double clhs587 =             clhs442 + clhs543;
const double clhs588 =             clhs584*clhs98;
const double clhs589 =             pow(DN(2,1), 2);
const double clhs590 =             clhs240*clhs551;
const double clhs591 =             DN(2,2)*clhs586;
const double clhs592 =             clhs240*clhs554;
const double clhs593 =             clhs148*clhs584;
const double clhs594 =             DN(2,1)*N[2];
const double clhs595 =             DN(3,0)*clhs586;
const double clhs596 =             clhs571 + clhs573;
const double clhs597 =             clhs563 + clhs566;
const double clhs598 =             DN(3,1)*clhs586 + clhs236*clhs564;
const double clhs599 =             DN(3,2)*clhs586;
const double clhs600 =             clhs198*clhs584;
const double clhs601 =             clhs198*clhs409;
const double clhs602 =             clhs148*clhs27*clhs332;
const double clhs603 =             -clhs601 - clhs602;
const double clhs604 =             DN(2,1)*N[3];
const double clhs605 =             DN(3,1)*N[2];
const double clhs606 =             DN(2,2)*clhs11;
const double clhs607 =             pow(DN(2,2), 2);
const double clhs608 =             DN(2,2)*N[2];
const double clhs609 =             DN(3,0)*clhs606;
const double clhs610 =             -clhs577 - clhs578;
const double clhs611 =             DN(3,1)*clhs606;
const double clhs612 =             clhs601 + clhs602;
const double clhs613 =             DN(3,2)*clhs606 + clhs346*clhs564;
const double clhs614 =             DN(2,2)*N[3];
const double clhs615 =             DN(3,2)*N[2];
const double clhs616 =             DN(3,0)*clhs178 + DN(3,1)*clhs312 + DN(3,2)*clhs383 + N[3]*clhs173;
const double clhs617 =             clhs18*clhs187;
const double clhs618 =             N[3]*clhs28;
const double clhs619 =             clhs184*clhs21 + clhs25*clhs617 + clhs25*clhs618;
const double clhs620 =             N[0]*clhs617;
const double clhs621 =             clhs37*clhs620;
const double clhs622 =             clhs46*clhs620;
const double clhs623 =             DN(3,0)*clhs11;
const double clhs624 =             clhs187*clhs6;
const double clhs625 =             clhs184*clhs87 + clhs617*clhs89 + clhs618*clhs89;
const double clhs626 =             N[1]*clhs617;
const double clhs627 =             clhs37*clhs626;
const double clhs628 =             clhs46*clhs626;
const double clhs629 =             clhs137*clhs184 + clhs139*clhs617 + clhs139*clhs618;
const double clhs630 =             N[2]*clhs617;
const double clhs631 =             clhs37*clhs630;
const double clhs632 =             clhs46*clhs630;
const double clhs633 =             pow(DN(3,0), 2);
const double clhs634 =             pow(N[3], 2);
const double clhs635 =             clhs19*clhs634;
const double clhs636 =             clhs184*clhs187 + clhs189*clhs617 + clhs189*clhs618 + clhs23*clhs634;
const double clhs637 =             clhs6*clhs634;
const double clhs638 =             clhs37*clhs637;
const double clhs639 =             DN(3,1)*clhs623;
const double clhs640 =             clhs28*clhs634;
const double clhs641 =             clhs37*clhs640;
const double clhs642 =             N[3]*clhs617;
const double clhs643 =             clhs37*clhs642;
const double clhs644 =             clhs46*clhs637;
const double clhs645 =             DN(3,2)*clhs623;
const double clhs646 =             clhs46*clhs640;
const double clhs647 =             clhs46*clhs642;
const double clhs648 =             DN(3,0)*N[3];
const double clhs649 =             clhs624*clhs73;
const double clhs650 =             clhs185 + clhs619;
const double clhs651 =             clhs187*clhs51;
const double clhs652 =             clhs47*clhs651;
const double clhs653 =             DN(3,1)*clhs11;
const double clhs654 =             clhs462 + clhs625;
const double clhs655 =             clhs651*clhs98;
const double clhs656 =             clhs563 + clhs629;
const double clhs657 =             clhs148*clhs651;
const double clhs658 =             pow(DN(3,1), 2);
const double clhs659 =             clhs240*clhs637;
const double clhs660 =             DN(3,2)*clhs653;
const double clhs661 =             clhs240*clhs640;
const double clhs662 =             clhs198*clhs651;
const double clhs663 =             DN(3,1)*N[3];
const double clhs664 =             DN(3,2)*clhs11;
const double clhs665 =             pow(DN(3,2), 2);
const double clhs666 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs14*clhs20 + clhs30;
            lhs(0,1)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs36 - clhs39 + clhs41 - clhs43 - clhs45 - clhs50*clhs52;
            lhs(0,2)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + DN(0,2)*clhs57 - clhs52*clhs64 + clhs58 + clhs59 + clhs60 + clhs61;
            lhs(0,3)=clhs40*clhs67 + clhs65*clhs74 - clhs65 + clhs70*clhs72 + clhs75*clhs76;
            lhs(0,4)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs81 + clhs86 + clhs90;
            lhs(0,5)=DN(0,0)*clhs91 + DN(0,1)*clhs93 + DN(0,2)*clhs96 - clhs101*clhs52 - clhs103 + clhs108 + clhs97;
            lhs(0,6)=DN(0,0)*clhs109 + DN(0,1)*clhs111 + DN(0,2)*clhs113 + clhs114 - clhs117*clhs52 + clhs118 + clhs121;
            lhs(0,7)=-clhs122 + clhs123*clhs40 + clhs126*clhs72 + clhs127*clhs74 + clhs128*clhs76;
            lhs(0,8)=DN(0,0)*clhs129 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs136 + clhs140;
            lhs(0,9)=DN(0,0)*clhs141 + DN(0,1)*clhs143 + DN(0,2)*clhs146 + clhs147 - clhs151*clhs52 - clhs153 + clhs158;
            lhs(0,10)=DN(0,0)*clhs159 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs164 - clhs167*clhs52 + clhs168 + clhs171;
            lhs(0,11)=-clhs172 + clhs173*clhs40 + clhs176*clhs72 + clhs177*clhs74 + clhs178*clhs76;
            lhs(0,12)=DN(0,0)*clhs179 + DN(0,1)*clhs181 + DN(0,2)*clhs183 + clhs186 + clhs190;
            lhs(0,13)=DN(0,0)*clhs191 + DN(0,1)*clhs193 + DN(0,2)*clhs196 + clhs197 - clhs201*clhs52 - clhs203 + clhs208;
            lhs(0,14)=DN(0,0)*clhs209 + DN(0,1)*clhs211 + DN(0,2)*clhs213 + clhs214 - clhs217*clhs52 + clhs218 + clhs221;
            lhs(0,15)=DN(3,0)*clhs228 - clhs222 + clhs223*clhs40 + clhs226*clhs72 + clhs227*clhs74;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs229 + DN(0,2)*clhs230 - clhs231*clhs52 + clhs39 + clhs41 + clhs43 + clhs45;
            lhs(1,1)=DN(0,0)*clhs33 + DN(0,1)*clhs232 + DN(0,2)*clhs234 + clhs11*clhs235 + clhs20*clhs236 + clhs30;
            lhs(1,2)=DN(0,0)*clhs55 + DN(0,1)*clhs237 + DN(0,2)*clhs239 - clhs241 + clhs243 - clhs244 - clhs246 - clhs250*clhs52;
            lhs(1,3)=clhs242*clhs67 + clhs251*clhs74 - clhs251 - clhs253*clhs72 + clhs254*clhs76;
            lhs(1,4)=DN(0,0)*clhs79 + DN(0,1)*clhs255 + DN(0,2)*clhs256 + clhs103 + clhs257 - clhs258*clhs52 + clhs259;
            lhs(1,5)=DN(0,0)*clhs93 + DN(0,1)*clhs260 + DN(0,2)*clhs262 + clhs263 + clhs265;
            lhs(1,6)=DN(0,0)*clhs111 + DN(0,1)*clhs266 + DN(0,2)*clhs268 + clhs269 - clhs272*clhs52 - clhs273 + clhs278;
            lhs(1,7)=clhs123*clhs242 - clhs279 - clhs282*clhs72 + clhs283*clhs74 + clhs284*clhs76;
            lhs(1,8)=DN(0,0)*clhs131 + DN(0,1)*clhs285 + DN(0,2)*clhs286 + clhs153 + clhs287 - clhs288*clhs52 + clhs289;
            lhs(1,9)=DN(0,0)*clhs143 + DN(0,1)*clhs290 + DN(0,2)*clhs292 + clhs293 + clhs294;
            lhs(1,10)=DN(0,0)*clhs161 + DN(0,1)*clhs295 + DN(0,2)*clhs297 + clhs298 - clhs301*clhs52 - clhs302 + clhs306;
            lhs(1,11)=clhs173*clhs242 - clhs307 - clhs310*clhs72 + clhs311*clhs74 + clhs312*clhs76;
            lhs(1,12)=DN(0,0)*clhs181 + DN(0,1)*clhs313 + DN(0,2)*clhs314 + clhs203 + clhs315 - clhs316*clhs52 + clhs317;
            lhs(1,13)=DN(0,0)*clhs193 + DN(0,1)*clhs318 + DN(0,2)*clhs320 + clhs321 + clhs322;
            lhs(1,14)=DN(0,0)*clhs211 + DN(0,1)*clhs323 + DN(0,2)*clhs325 + clhs326 - clhs329*clhs52 - clhs330 + clhs334;
            lhs(1,15)=DN(3,1)*clhs228 + clhs223*clhs242 - clhs335 - clhs338*clhs72 + clhs339*clhs74;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs230 + DN(0,2)*clhs340 - clhs341*clhs52 - clhs58 + clhs59 - clhs60 - clhs61;
            lhs(2,1)=DN(0,0)*clhs36 + DN(0,1)*clhs234 + DN(0,2)*clhs342 + clhs241 + clhs243 + clhs244 + clhs246 - clhs343*clhs52;
            lhs(2,2)=DN(0,0)*clhs57 + DN(0,1)*clhs239 + DN(0,2)*clhs344 + clhs11*clhs345 + clhs20*clhs346 + clhs30;
            lhs(2,3)=clhs347*clhs74 - clhs347 + clhs348*clhs67 + clhs349*clhs72 + clhs350*clhs76;
            lhs(2,4)=DN(0,0)*clhs81 + DN(0,1)*clhs256 + DN(0,2)*clhs351 - clhs118 + clhs352 - clhs353*clhs52 + clhs354;
            lhs(2,5)=DN(0,0)*clhs96 + DN(0,1)*clhs262 + DN(0,2)*clhs355 + clhs273 + clhs356 - clhs357*clhs52 + clhs358;
            lhs(2,6)=DN(0,0)*clhs113 + DN(0,1)*clhs268 + DN(0,2)*clhs359 + clhs263 + clhs361;
            lhs(2,7)=clhs123*clhs348 - clhs362 + clhs365*clhs72 + clhs366*clhs74 + clhs367*clhs76;
            lhs(2,8)=DN(0,0)*clhs133 + DN(0,1)*clhs286 + DN(0,2)*clhs368 - clhs168 + clhs369 - clhs370*clhs52 + clhs371;
            lhs(2,9)=DN(0,0)*clhs146 + DN(0,1)*clhs292 + DN(0,2)*clhs372 + clhs302 + clhs373 - clhs374*clhs52 + clhs375;
            lhs(2,10)=DN(0,0)*clhs163 + DN(0,1)*clhs297 + DN(0,2)*clhs376 + clhs293 + clhs377;
            lhs(2,11)=clhs173*clhs348 - clhs378 + clhs381*clhs72 + clhs382*clhs74 + clhs383*clhs76;
            lhs(2,12)=DN(0,0)*clhs183 + DN(0,1)*clhs314 + DN(0,2)*clhs384 - clhs218 + clhs385 - clhs386*clhs52 + clhs387;
            lhs(2,13)=DN(0,0)*clhs196 + DN(0,1)*clhs320 + DN(0,2)*clhs388 + clhs330 + clhs389 - clhs390*clhs52 + clhs391;
            lhs(2,14)=DN(0,0)*clhs213 + DN(0,1)*clhs325 + DN(0,2)*clhs392 + clhs321 + clhs393;
            lhs(2,15)=DN(3,2)*clhs228 + clhs223*clhs348 - clhs394 + clhs397*clhs72 + clhs398*clhs74;
            lhs(3,0)=DN(0,0)*clhs400 + clhs65 + clhs68*clhs72 - clhs69*clhs72;
            lhs(3,1)=DN(0,1)*clhs400 + clhs251 + clhs252*clhs72 - clhs37*clhs401;
            lhs(3,2)=DN(0,2)*clhs400 - clhs251*clhs402 + clhs347 + clhs401*clhs46;
            lhs(3,3)=clhs15*clhs403 + clhs235*clhs73 + clhs345*clhs73 + clhs5*clhs73;
            lhs(3,4)=DN(0,0)*clhs406 + clhs127 + clhs405*clhs68 - clhs405*clhs69;
            lhs(3,5)=DN(0,1)*clhs406 + clhs252*clhs405 + clhs283 - clhs37*clhs407;
            lhs(3,6)=DN(0,2)*clhs406 - clhs279*clhs402 + clhs366 + clhs407*clhs46;
            lhs(3,7)=clhs408;
            lhs(3,8)=DN(0,0)*clhs411 + clhs177 + clhs410*clhs68 - clhs410*clhs69;
            lhs(3,9)=DN(0,1)*clhs411 + clhs252*clhs410 + clhs311 - clhs37*clhs412;
            lhs(3,10)=DN(0,2)*clhs411 - clhs307*clhs402 + clhs382 + clhs412*clhs46;
            lhs(3,11)=clhs413;
            lhs(3,12)=DN(0,0)*clhs415 + clhs227 + clhs414*clhs68 - clhs414*clhs69;
            lhs(3,13)=DN(0,1)*clhs415 + clhs252*clhs414 + clhs339 - clhs37*clhs416;
            lhs(3,14)=DN(0,2)*clhs415 - clhs335*clhs402 + clhs398 + clhs416*clhs46;
            lhs(3,15)=clhs417;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs420 + clhs86;
            lhs(4,1)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs36 + clhs108 + clhs257 - clhs275*clhs50 - clhs422;
            lhs(4,2)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + DN(1,2)*clhs57 + clhs121 - clhs275*clhs64 + clhs352 + clhs423;
            lhs(4,3)=clhs122*clhs74 - clhs127 + clhs405*clhs70 + clhs424*clhs67 + clhs425*clhs75;
            lhs(4,4)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs81 + clhs11*clhs426 + clhs14*clhs428 + clhs429;
            lhs(4,5)=DN(1,0)*clhs91 + DN(1,1)*clhs93 + DN(1,2)*clhs96 - clhs101*clhs275 - clhs431 + clhs432 - clhs434 - clhs436;
            lhs(4,6)=DN(1,0)*clhs109 + DN(1,1)*clhs111 + DN(1,2)*clhs113 - clhs117*clhs275 + clhs437 + clhs438 + clhs439 + clhs440;
            lhs(4,7)=clhs123*clhs424 + clhs126*clhs405 + clhs128*clhs425 + clhs441*clhs74 - clhs441;
            lhs(4,8)=DN(1,0)*clhs129 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs445 + clhs446;
            lhs(4,9)=DN(1,0)*clhs141 + DN(1,1)*clhs143 + DN(1,2)*clhs146 - clhs151*clhs275 + clhs447 - clhs449 + clhs454;
            lhs(4,10)=DN(1,0)*clhs159 + DN(1,1)*clhs161 + DN(1,2)*clhs163 - clhs167*clhs275 + clhs455 + clhs456 + clhs459;
            lhs(4,11)=clhs173*clhs424 + clhs176*clhs405 + clhs178*clhs425 - clhs460 + clhs461*clhs74;
            lhs(4,12)=DN(1,0)*clhs179 + DN(1,1)*clhs181 + DN(1,2)*clhs183 + clhs463 + clhs464;
            lhs(4,13)=DN(1,0)*clhs191 + DN(1,1)*clhs193 + DN(1,2)*clhs196 - clhs201*clhs275 + clhs465 - clhs467 + clhs472;
            lhs(4,14)=DN(1,0)*clhs209 + DN(1,1)*clhs211 + DN(1,2)*clhs213 - clhs217*clhs275 + clhs473 + clhs474 + clhs477;
            lhs(4,15)=DN(3,0)*clhs480 + clhs223*clhs424 + clhs226*clhs405 - clhs478 + clhs479*clhs74;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs229 + DN(1,2)*clhs230 - clhs231*clhs275 + clhs259 + clhs422 + clhs97;
            lhs(5,1)=DN(1,0)*clhs33 + DN(1,1)*clhs232 + DN(1,2)*clhs234 + clhs265 + clhs481;
            lhs(5,2)=DN(1,0)*clhs55 + DN(1,1)*clhs237 + DN(1,2)*clhs239 - clhs250*clhs275 + clhs278 + clhs356 - clhs483;
            lhs(5,3)=-clhs253*clhs405 + clhs254*clhs425 + clhs279*clhs74 - clhs283 + clhs484*clhs67;
            lhs(5,4)=DN(1,0)*clhs79 + DN(1,1)*clhs255 + DN(1,2)*clhs256 - clhs258*clhs275 + clhs431 + clhs432 + clhs434 + clhs436;
            lhs(5,5)=DN(1,0)*clhs93 + DN(1,1)*clhs260 + DN(1,2)*clhs262 + clhs11*clhs485 + clhs236*clhs428 + clhs429;
            lhs(5,6)=DN(1,0)*clhs111 + DN(1,1)*clhs266 + DN(1,2)*clhs268 - clhs272*clhs275 - clhs486 + clhs487 - clhs488 - clhs489;
            lhs(5,7)=clhs123*clhs484 - clhs282*clhs405 + clhs284*clhs425 + clhs490*clhs74 - clhs490;
            lhs(5,8)=DN(1,0)*clhs131 + DN(1,1)*clhs285 + DN(1,2)*clhs286 - clhs275*clhs288 + clhs449 + clhs491 + clhs492;
            lhs(5,9)=DN(1,0)*clhs143 + DN(1,1)*clhs290 + DN(1,2)*clhs292 + clhs493 + clhs495;
            lhs(5,10)=DN(1,0)*clhs161 + DN(1,1)*clhs295 + DN(1,2)*clhs297 - clhs275*clhs301 + clhs496 - clhs497 + clhs501;
            lhs(5,11)=clhs173*clhs484 - clhs310*clhs405 + clhs312*clhs425 - clhs502 + clhs503*clhs74;
            lhs(5,12)=DN(1,0)*clhs181 + DN(1,1)*clhs313 + DN(1,2)*clhs314 - clhs275*clhs316 + clhs467 + clhs504 + clhs505;
            lhs(5,13)=DN(1,0)*clhs193 + DN(1,1)*clhs318 + DN(1,2)*clhs320 + clhs506 + clhs507;
            lhs(5,14)=DN(1,0)*clhs211 + DN(1,1)*clhs323 + DN(1,2)*clhs325 - clhs275*clhs329 + clhs508 - clhs509 + clhs512;
            lhs(5,15)=DN(3,1)*clhs480 + clhs223*clhs484 - clhs338*clhs405 - clhs513 + clhs514*clhs74;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs230 + DN(1,2)*clhs340 + clhs114 - clhs275*clhs341 + clhs354 - clhs423;
            lhs(6,1)=DN(1,0)*clhs36 + DN(1,1)*clhs234 + DN(1,2)*clhs342 + clhs269 - clhs275*clhs343 + clhs358 + clhs483;
            lhs(6,2)=DN(1,0)*clhs57 + DN(1,1)*clhs239 + DN(1,2)*clhs344 + clhs361 + clhs481;
            lhs(6,3)=clhs349*clhs405 + clhs350*clhs425 + clhs362*clhs74 - clhs366 + clhs515*clhs67;
            lhs(6,4)=DN(1,0)*clhs81 + DN(1,1)*clhs256 + DN(1,2)*clhs351 - clhs275*clhs353 - clhs437 + clhs438 - clhs439 - clhs440;
            lhs(6,5)=DN(1,0)*clhs96 + DN(1,1)*clhs262 + DN(1,2)*clhs355 - clhs275*clhs357 + clhs486 + clhs487 + clhs488 + clhs489;
            lhs(6,6)=DN(1,0)*clhs113 + DN(1,1)*clhs268 + DN(1,2)*clhs359 + clhs11*clhs516 + clhs346*clhs428 + clhs429;
            lhs(6,7)=clhs123*clhs515 + clhs365*clhs405 + clhs367*clhs425 + clhs517*clhs74 - clhs517;
            lhs(6,8)=DN(1,0)*clhs133 + DN(1,1)*clhs286 + DN(1,2)*clhs368 - clhs275*clhs370 - clhs456 + clhs518 + clhs519;
            lhs(6,9)=DN(1,0)*clhs146 + DN(1,1)*clhs292 + DN(1,2)*clhs372 - clhs275*clhs374 + clhs497 + clhs520 + clhs521;
            lhs(6,10)=DN(1,0)*clhs163 + DN(1,1)*clhs297 + DN(1,2)*clhs376 + clhs493 + clhs523;
            lhs(6,11)=clhs173*clhs515 + clhs381*clhs405 + clhs383*clhs425 - clhs524 + clhs525*clhs74;
            lhs(6,12)=DN(1,0)*clhs183 + DN(1,1)*clhs314 + DN(1,2)*clhs384 - clhs275*clhs386 - clhs474 + clhs526 + clhs527;
            lhs(6,13)=DN(1,0)*clhs196 + DN(1,1)*clhs320 + DN(1,2)*clhs388 - clhs275*clhs390 + clhs509 + clhs528 + clhs529;
            lhs(6,14)=DN(1,0)*clhs213 + DN(1,1)*clhs325 + DN(1,2)*clhs392 + clhs506 + clhs530;
            lhs(6,15)=DN(3,2)*clhs480 + clhs223*clhs515 + clhs397*clhs405 - clhs531 + clhs532*clhs74;
            lhs(7,0)=DN(1,0)*clhs400 + clhs122 + clhs124*clhs72 - clhs125*clhs72;
            lhs(7,1)=DN(1,1)*clhs400 + clhs279 - clhs280*clhs72 + clhs281*clhs72;
            lhs(7,2)=DN(1,2)*clhs400 + clhs362 + clhs363*clhs72 - clhs364*clhs72;
            lhs(7,3)=clhs408;
            lhs(7,4)=DN(1,0)*clhs406 + clhs124*clhs405 - clhs125*clhs405 + clhs441;
            lhs(7,5)=DN(1,1)*clhs406 - clhs280*clhs405 + clhs281*clhs405 + clhs490;
            lhs(7,6)=DN(1,2)*clhs406 + clhs363*clhs405 - clhs364*clhs405 + clhs517;
            lhs(7,7)=clhs403*clhs427 + clhs426*clhs73 + clhs485*clhs73 + clhs516*clhs73;
            lhs(7,8)=DN(1,0)*clhs411 + clhs124*clhs410 - clhs125*clhs410 + clhs461;
            lhs(7,9)=DN(1,1)*clhs411 - clhs280*clhs410 + clhs281*clhs410 + clhs503;
            lhs(7,10)=DN(1,2)*clhs411 + clhs363*clhs410 - clhs364*clhs410 + clhs525;
            lhs(7,11)=clhs533;
            lhs(7,12)=DN(1,0)*clhs415 + clhs124*clhs414 - clhs125*clhs414 + clhs479;
            lhs(7,13)=DN(1,1)*clhs415 - clhs280*clhs414 + clhs281*clhs414 + clhs514;
            lhs(7,14)=DN(1,2)*clhs415 + clhs363*clhs414 - clhs364*clhs414 + clhs532;
            lhs(7,15)=clhs534;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs136 + clhs537;
            lhs(8,1)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs36 + clhs158 + clhs287 - clhs304*clhs50 - clhs539;
            lhs(8,2)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + DN(2,2)*clhs57 + clhs171 - clhs304*clhs64 + clhs369 + clhs540;
            lhs(8,3)=clhs172*clhs74 - clhs177 + clhs410*clhs70 + clhs541*clhs67 + clhs542*clhs75;
            lhs(8,4)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs81 + clhs445 + clhs543;
            lhs(8,5)=DN(2,0)*clhs91 + DN(2,1)*clhs93 + DN(2,2)*clhs96 - clhs101*clhs304 + clhs454 + clhs491 - clhs545;
            lhs(8,6)=DN(2,0)*clhs109 + DN(2,1)*clhs111 + DN(2,2)*clhs113 - clhs117*clhs304 + clhs459 + clhs518 + clhs546;
            lhs(8,7)=clhs123*clhs541 + clhs126*clhs410 + clhs128*clhs542 + clhs460*clhs74 - clhs461;
            lhs(8,8)=DN(2,0)*clhs129 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs11*clhs547 + clhs14*clhs549 + clhs550;
            lhs(8,9)=DN(2,0)*clhs141 + DN(2,1)*clhs143 + DN(2,2)*clhs146 - clhs151*clhs304 - clhs552 + clhs553 - clhs555 - clhs557;
            lhs(8,10)=DN(2,0)*clhs159 + DN(2,1)*clhs161 + DN(2,2)*clhs163 - clhs167*clhs304 + clhs558 + clhs559 + clhs560 + clhs561;
            lhs(8,11)=clhs173*clhs541 + clhs176*clhs410 + clhs178*clhs542 + clhs562*clhs74 - clhs562;
            lhs(8,12)=DN(2,0)*clhs179 + DN(2,1)*clhs181 + DN(2,2)*clhs183 + clhs565 + clhs566;
            lhs(8,13)=DN(2,0)*clhs191 + DN(2,1)*clhs193 + DN(2,2)*clhs196 - clhs201*clhs304 + clhs567 - clhs569 + clhs574;
            lhs(8,14)=DN(2,0)*clhs209 + DN(2,1)*clhs211 + DN(2,2)*clhs213 - clhs217*clhs304 + clhs575 + clhs576 + clhs579;
            lhs(8,15)=DN(3,0)*clhs582 + clhs223*clhs541 + clhs226*clhs410 - clhs580 + clhs581*clhs74;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs229 + DN(2,2)*clhs230 + clhs147 - clhs231*clhs304 + clhs289 + clhs539;
            lhs(9,1)=DN(2,0)*clhs33 + DN(2,1)*clhs232 + DN(2,2)*clhs234 + clhs294 + clhs583;
            lhs(9,2)=DN(2,0)*clhs55 + DN(2,1)*clhs237 + DN(2,2)*clhs239 - clhs250*clhs304 + clhs306 + clhs373 - clhs585;
            lhs(9,3)=-clhs253*clhs410 + clhs254*clhs542 + clhs307*clhs74 - clhs311 + clhs586*clhs67;
            lhs(9,4)=DN(2,0)*clhs79 + DN(2,1)*clhs255 + DN(2,2)*clhs256 - clhs258*clhs304 + clhs447 + clhs492 + clhs545;
            lhs(9,5)=DN(2,0)*clhs93 + DN(2,1)*clhs260 + DN(2,2)*clhs262 + clhs495 + clhs587;
            lhs(9,6)=DN(2,0)*clhs111 + DN(2,1)*clhs266 + DN(2,2)*clhs268 - clhs272*clhs304 + clhs501 + clhs520 - clhs588;
            lhs(9,7)=clhs123*clhs586 - clhs282*clhs410 + clhs284*clhs542 + clhs502*clhs74 - clhs503;
            lhs(9,8)=DN(2,0)*clhs131 + DN(2,1)*clhs285 + DN(2,2)*clhs286 - clhs288*clhs304 + clhs552 + clhs553 + clhs555 + clhs557;
            lhs(9,9)=DN(2,0)*clhs143 + DN(2,1)*clhs290 + DN(2,2)*clhs292 + clhs11*clhs589 + clhs236*clhs549 + clhs550;
            lhs(9,10)=DN(2,0)*clhs161 + DN(2,1)*clhs295 + DN(2,2)*clhs297 - clhs301*clhs304 - clhs590 + clhs591 - clhs592 - clhs593;
            lhs(9,11)=clhs173*clhs586 - clhs310*clhs410 + clhs312*clhs542 + clhs594*clhs74 - clhs594;
            lhs(9,12)=DN(2,0)*clhs181 + DN(2,1)*clhs313 + DN(2,2)*clhs314 - clhs304*clhs316 + clhs569 + clhs595 + clhs596;
            lhs(9,13)=DN(2,0)*clhs193 + DN(2,1)*clhs318 + DN(2,2)*clhs320 + clhs597 + clhs598;
            lhs(9,14)=DN(2,0)*clhs211 + DN(2,1)*clhs323 + DN(2,2)*clhs325 - clhs304*clhs329 + clhs599 - clhs600 + clhs603;
            lhs(9,15)=DN(3,1)*clhs582 + clhs223*clhs586 - clhs338*clhs410 - clhs604 + clhs605*clhs74;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs230 + DN(2,2)*clhs340 + clhs164 - clhs304*clhs341 + clhs371 - clhs540;
            lhs(10,1)=DN(2,0)*clhs36 + DN(2,1)*clhs234 + DN(2,2)*clhs342 + clhs298 - clhs304*clhs343 + clhs375 + clhs585;
            lhs(10,2)=DN(2,0)*clhs57 + DN(2,1)*clhs239 + DN(2,2)*clhs344 + clhs377 + clhs583;
            lhs(10,3)=clhs349*clhs410 + clhs350*clhs542 + clhs378*clhs74 - clhs382 + clhs606*clhs67;
            lhs(10,4)=DN(2,0)*clhs81 + DN(2,1)*clhs256 + DN(2,2)*clhs351 - clhs304*clhs353 + clhs455 + clhs519 - clhs546;
            lhs(10,5)=DN(2,0)*clhs96 + DN(2,1)*clhs262 + DN(2,2)*clhs355 - clhs304*clhs357 + clhs496 + clhs521 + clhs588;
            lhs(10,6)=DN(2,0)*clhs113 + DN(2,1)*clhs268 + DN(2,2)*clhs359 + clhs523 + clhs587;
            lhs(10,7)=clhs123*clhs606 + clhs365*clhs410 + clhs367*clhs542 + clhs524*clhs74 - clhs525;
            lhs(10,8)=DN(2,0)*clhs133 + DN(2,1)*clhs286 + DN(2,2)*clhs368 - clhs304*clhs370 - clhs558 + clhs559 - clhs560 - clhs561;
            lhs(10,9)=DN(2,0)*clhs146 + DN(2,1)*clhs292 + DN(2,2)*clhs372 - clhs304*clhs374 + clhs590 + clhs591 + clhs592 + clhs593;
            lhs(10,10)=DN(2,0)*clhs163 + DN(2,1)*clhs297 + DN(2,2)*clhs376 + clhs11*clhs607 + clhs346*clhs549 + clhs550;
            lhs(10,11)=clhs173*clhs606 + clhs381*clhs410 + clhs383*clhs542 + clhs608*clhs74 - clhs608;
            lhs(10,12)=DN(2,0)*clhs183 + DN(2,1)*clhs314 + DN(2,2)*clhs384 - clhs304*clhs386 - clhs576 + clhs609 + clhs610;
            lhs(10,13)=DN(2,0)*clhs196 + DN(2,1)*clhs320 + DN(2,2)*clhs388 - clhs304*clhs390 + clhs600 + clhs611 + clhs612;
            lhs(10,14)=DN(2,0)*clhs213 + DN(2,1)*clhs325 + DN(2,2)*clhs392 + clhs597 + clhs613;
            lhs(10,15)=DN(3,2)*clhs582 + clhs223*clhs606 + clhs397*clhs410 - clhs614 + clhs615*clhs74;
            lhs(11,0)=DN(2,0)*clhs400 + clhs172 + clhs174*clhs72 - clhs175*clhs72;
            lhs(11,1)=DN(2,1)*clhs400 + clhs307 - clhs308*clhs72 + clhs309*clhs72;
            lhs(11,2)=DN(2,2)*clhs400 + clhs378 + clhs379*clhs72 - clhs380*clhs72;
            lhs(11,3)=clhs413;
            lhs(11,4)=DN(2,0)*clhs406 + clhs174*clhs405 - clhs175*clhs405 + clhs460;
            lhs(11,5)=DN(2,1)*clhs406 - clhs308*clhs405 + clhs309*clhs405 + clhs502;
            lhs(11,6)=DN(2,2)*clhs406 + clhs379*clhs405 - clhs380*clhs405 + clhs524;
            lhs(11,7)=clhs533;
            lhs(11,8)=DN(2,0)*clhs411 + clhs174*clhs410 - clhs175*clhs410 + clhs562;
            lhs(11,9)=DN(2,1)*clhs411 - clhs308*clhs410 + clhs309*clhs410 + clhs594;
            lhs(11,10)=DN(2,2)*clhs411 + clhs379*clhs410 - clhs380*clhs410 + clhs608;
            lhs(11,11)=clhs403*clhs548 + clhs547*clhs73 + clhs589*clhs73 + clhs607*clhs73;
            lhs(11,12)=DN(2,0)*clhs415 + clhs174*clhs414 - clhs175*clhs414 + clhs581;
            lhs(11,13)=DN(2,1)*clhs415 - clhs308*clhs414 + clhs309*clhs414 + clhs605;
            lhs(11,14)=DN(2,2)*clhs415 + clhs379*clhs414 - clhs380*clhs414 + clhs615;
            lhs(11,15)=clhs616;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs186 + clhs619;
            lhs(12,1)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs36 + clhs208 + clhs315 - clhs332*clhs50 - clhs621;
            lhs(12,2)=DN(3,0)*clhs53 + DN(3,1)*clhs55 + DN(3,2)*clhs57 + clhs221 - clhs332*clhs64 + clhs385 + clhs622;
            lhs(12,3)=clhs222*clhs74 - clhs227 + clhs414*clhs70 + clhs623*clhs67 + clhs624*clhs75;
            lhs(12,4)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs81 + clhs463 + clhs625;
            lhs(12,5)=DN(3,0)*clhs91 + DN(3,1)*clhs93 + DN(3,2)*clhs96 - clhs101*clhs332 + clhs472 + clhs504 - clhs627;
            lhs(12,6)=DN(3,0)*clhs109 + DN(3,1)*clhs111 + DN(3,2)*clhs113 - clhs117*clhs332 + clhs477 + clhs526 + clhs628;
            lhs(12,7)=clhs123*clhs623 + clhs126*clhs414 + clhs128*clhs624 + clhs478*clhs74 - clhs479;
            lhs(12,8)=DN(3,0)*clhs129 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs565 + clhs629;
            lhs(12,9)=DN(3,0)*clhs141 + DN(3,1)*clhs143 + DN(3,2)*clhs146 - clhs151*clhs332 + clhs574 + clhs595 - clhs631;
            lhs(12,10)=DN(3,0)*clhs159 + DN(3,1)*clhs161 + DN(3,2)*clhs163 - clhs167*clhs332 + clhs579 + clhs609 + clhs632;
            lhs(12,11)=clhs173*clhs623 + clhs176*clhs414 + clhs178*clhs624 + clhs580*clhs74 - clhs581;
            lhs(12,12)=DN(3,0)*clhs179 + DN(3,1)*clhs181 + DN(3,2)*clhs183 + clhs11*clhs633 + clhs14*clhs635 + clhs636;
            lhs(12,13)=DN(3,0)*clhs191 + DN(3,1)*clhs193 + DN(3,2)*clhs196 - clhs201*clhs332 - clhs638 + clhs639 - clhs641 - clhs643;
            lhs(12,14)=DN(3,0)*clhs209 + DN(3,1)*clhs211 + DN(3,2)*clhs213 - clhs217*clhs332 + clhs644 + clhs645 + clhs646 + clhs647;
            lhs(12,15)=DN(3,0)*clhs649 + clhs223*clhs623 + clhs226*clhs414 + clhs648*clhs74 - clhs648;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs229 + DN(3,2)*clhs230 + clhs197 - clhs231*clhs332 + clhs317 + clhs621;
            lhs(13,1)=DN(3,0)*clhs33 + DN(3,1)*clhs232 + DN(3,2)*clhs234 + clhs322 + clhs650;
            lhs(13,2)=DN(3,0)*clhs55 + DN(3,1)*clhs237 + DN(3,2)*clhs239 - clhs250*clhs332 + clhs334 + clhs389 - clhs652;
            lhs(13,3)=-clhs253*clhs414 + clhs254*clhs624 + clhs335*clhs74 - clhs339 + clhs653*clhs67;
            lhs(13,4)=DN(3,0)*clhs79 + DN(3,1)*clhs255 + DN(3,2)*clhs256 - clhs258*clhs332 + clhs465 + clhs505 + clhs627;
            lhs(13,5)=DN(3,0)*clhs93 + DN(3,1)*clhs260 + DN(3,2)*clhs262 + clhs507 + clhs654;
            lhs(13,6)=DN(3,0)*clhs111 + DN(3,1)*clhs266 + DN(3,2)*clhs268 - clhs272*clhs332 + clhs512 + clhs528 - clhs655;
            lhs(13,7)=clhs123*clhs653 - clhs282*clhs414 + clhs284*clhs624 + clhs513*clhs74 - clhs514;
            lhs(13,8)=DN(3,0)*clhs131 + DN(3,1)*clhs285 + DN(3,2)*clhs286 - clhs288*clhs332 + clhs567 + clhs596 + clhs631;
            lhs(13,9)=DN(3,0)*clhs143 + DN(3,1)*clhs290 + DN(3,2)*clhs292 + clhs598 + clhs656;
            lhs(13,10)=DN(3,0)*clhs161 + DN(3,1)*clhs295 + DN(3,2)*clhs297 - clhs301*clhs332 + clhs603 + clhs611 - clhs657;
            lhs(13,11)=clhs173*clhs653 - clhs310*clhs414 + clhs312*clhs624 + clhs604*clhs74 - clhs605;
            lhs(13,12)=DN(3,0)*clhs181 + DN(3,1)*clhs313 + DN(3,2)*clhs314 - clhs316*clhs332 + clhs638 + clhs639 + clhs641 + clhs643;
            lhs(13,13)=DN(3,0)*clhs193 + DN(3,1)*clhs318 + DN(3,2)*clhs320 + clhs11*clhs658 + clhs236*clhs635 + clhs636;
            lhs(13,14)=DN(3,0)*clhs211 + DN(3,1)*clhs323 + DN(3,2)*clhs325 - clhs329*clhs332 - clhs659 + clhs660 - clhs661 - clhs662;
            lhs(13,15)=DN(3,1)*clhs649 + clhs223*clhs653 - clhs338*clhs414 + clhs663*clhs74 - clhs663;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs230 + DN(3,2)*clhs340 + clhs214 - clhs332*clhs341 + clhs387 - clhs622;
            lhs(14,1)=DN(3,0)*clhs36 + DN(3,1)*clhs234 + DN(3,2)*clhs342 + clhs326 - clhs332*clhs343 + clhs391 + clhs652;
            lhs(14,2)=DN(3,0)*clhs57 + DN(3,1)*clhs239 + DN(3,2)*clhs344 + clhs393 + clhs650;
            lhs(14,3)=clhs349*clhs414 + clhs350*clhs624 + clhs394*clhs74 - clhs398 + clhs664*clhs67;
            lhs(14,4)=DN(3,0)*clhs81 + DN(3,1)*clhs256 + DN(3,2)*clhs351 - clhs332*clhs353 + clhs473 + clhs527 - clhs628;
            lhs(14,5)=DN(3,0)*clhs96 + DN(3,1)*clhs262 + DN(3,2)*clhs355 - clhs332*clhs357 + clhs508 + clhs529 + clhs655;
            lhs(14,6)=DN(3,0)*clhs113 + DN(3,1)*clhs268 + DN(3,2)*clhs359 + clhs530 + clhs654;
            lhs(14,7)=clhs123*clhs664 + clhs365*clhs414 + clhs367*clhs624 + clhs531*clhs74 - clhs532;
            lhs(14,8)=DN(3,0)*clhs133 + DN(3,1)*clhs286 + DN(3,2)*clhs368 - clhs332*clhs370 + clhs575 + clhs610 - clhs632;
            lhs(14,9)=DN(3,0)*clhs146 + DN(3,1)*clhs292 + DN(3,2)*clhs372 - clhs332*clhs374 + clhs599 + clhs612 + clhs657;
            lhs(14,10)=DN(3,0)*clhs163 + DN(3,1)*clhs297 + DN(3,2)*clhs376 + clhs613 + clhs656;
            lhs(14,11)=clhs173*clhs664 + clhs381*clhs414 + clhs383*clhs624 + clhs614*clhs74 - clhs615;
            lhs(14,12)=DN(3,0)*clhs183 + DN(3,1)*clhs314 + DN(3,2)*clhs384 - clhs332*clhs386 - clhs644 + clhs645 - clhs646 - clhs647;
            lhs(14,13)=DN(3,0)*clhs196 + DN(3,1)*clhs320 + DN(3,2)*clhs388 - clhs332*clhs390 + clhs659 + clhs660 + clhs661 + clhs662;
            lhs(14,14)=DN(3,0)*clhs213 + DN(3,1)*clhs325 + DN(3,2)*clhs392 + clhs11*clhs665 + clhs346*clhs635 + clhs636;
            lhs(14,15)=DN(3,2)*clhs649 + clhs223*clhs664 + clhs397*clhs414 + clhs666*clhs74 - clhs666;
            lhs(15,0)=DN(3,0)*clhs400 + clhs222 + clhs224*clhs72 - clhs225*clhs72;
            lhs(15,1)=DN(3,1)*clhs400 + clhs335 - clhs336*clhs72 + clhs337*clhs72;
            lhs(15,2)=DN(3,2)*clhs400 + clhs394 + clhs395*clhs72 - clhs396*clhs72;
            lhs(15,3)=clhs417;
            lhs(15,4)=DN(3,0)*clhs406 + clhs224*clhs405 - clhs225*clhs405 + clhs478;
            lhs(15,5)=DN(3,1)*clhs406 - clhs336*clhs405 + clhs337*clhs405 + clhs513;
            lhs(15,6)=DN(3,2)*clhs406 + clhs395*clhs405 - clhs396*clhs405 + clhs531;
            lhs(15,7)=clhs534;
            lhs(15,8)=DN(3,0)*clhs411 + clhs224*clhs410 - clhs225*clhs410 + clhs580;
            lhs(15,9)=DN(3,1)*clhs411 - clhs336*clhs410 + clhs337*clhs410 + clhs604;
            lhs(15,10)=DN(3,2)*clhs411 + clhs395*clhs410 - clhs396*clhs410 + clhs614;
            lhs(15,11)=clhs616;
            lhs(15,12)=DN(3,0)*clhs415 + clhs224*clhs414 - clhs225*clhs414 + clhs648;
            lhs(15,13)=DN(3,1)*clhs415 - clhs336*clhs414 + clhs337*clhs414 + clhs663;
            lhs(15,14)=DN(3,2)*clhs415 + clhs395*clhs414 - clhs396*clhs414 + clhs666;
            lhs(15,15)=clhs403*clhs634 + clhs633*clhs73 + clhs658*clhs73 + clhs665*clhs73;


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
    const array_1d<double,3>& omega_old = rData.AngularVelocityOld;
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
const double crhs5 =             omega[2] - omega_old[2];
const double crhs6 =             crhs5*gauss_pt_coord[1];
const double crhs7 =             gauss_pt_coord[2]*(omega[1] - omega_old[1]);
const double crhs8 =             crhs1/dt;
const double crhs9 =             crhs8*(crhs6 - crhs7);
const double crhs10 =             gauss_pt_coord[1]*omega[0];
const double crhs11 =             gauss_pt_coord[0]*omega[1];
const double crhs12 =             crhs10 - crhs11;
const double crhs13 =             crhs1*(crhs12*omega[1] - omega[2]*(gauss_pt_coord[0]*omega[2] - gauss_pt_coord[2]*omega[0]));
const double crhs14 =             crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs15 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs16 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs17 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs18 =             crhs1*(crhs15*crhs16 + crhs17*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs19 =             crhs1*stab_c2*sqrt(pow(crhs16, 2) + pow(crhs17, 2));
const double crhs20 =             1.0/crhs1;
const double crhs21 =             crhs20*(crhs16*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crhs17*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crhs22 =             crhs20*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crhs23 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs24 =             crhs15 + crhs23;
const double crhs25 =             (crhs19*h/stab_c1 + mu)*(crhs21 + crhs22 + crhs24);
const double crhs26 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs27 =             crhs3*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs28 =             crhs8*(crhs5*gauss_pt_coord[0] - gauss_pt_coord[2]*(omega[0] - omega_old[0]));
const double crhs29 =             gauss_pt_coord[2]*omega[1];
const double crhs30 =             gauss_pt_coord[1]*omega[2];
const double crhs31 =             crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs32 =             crhs1*(crhs16*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs17*crhs23);
const double crhs33 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs1*(-crhs12*omega[0] + omega[2]*(crhs29 - crhs30)) - crhs26 + crhs27 + crhs28 + crhs31 + crhs32;
const double crhs34 =             1.0/(crhs1*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2)) + crhs19/h + crhs8*dyn_tau + mu*stab_c1/pow(h, 2));
const double crhs35 =             crhs3*crhs34;
const double crhs36 =             N[0]*crhs35;
const double crhs37 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs13 + crhs14 + crhs18 - crhs2 - crhs4 + crhs8*(-crhs6 + crhs7);
const double crhs38 =             1.0*crhs34;
const double crhs39 =             crhs37*crhs38;
const double crhs40 =             crhs1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs41 =             N[0]*crhs40;
const double crhs42 =             crhs1*(DN(0,0)*crhs16 + DN(0,1)*crhs17);
const double crhs43 =             crhs1*(-omega[0]*(-crhs10 + crhs11) + omega[2]*(-crhs29 + crhs30));
const double crhs44 =             crhs33*crhs38;
const double crhs45 =             N[1]*crhs35;
const double crhs46 =             N[1]*crhs40;
const double crhs47 =             crhs1*(DN(1,0)*crhs16 + DN(1,1)*crhs17);
const double crhs48 =             N[2]*crhs35;
const double crhs49 =             N[2]*crhs40;
const double crhs50 =             crhs1*(DN(2,0)*crhs16 + DN(2,1)*crhs17);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs25 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs13 - N[0]*crhs14 - N[0]*crhs18 + N[0]*crhs2 + N[0]*crhs4 + N[0]*crhs9 - crhs33*crhs36 - crhs39*crhs41 - crhs39*crhs42;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs25 - DN(0,1)*stress[1] + N[0]*crhs26 - N[0]*crhs27 - N[0]*crhs28 - N[0]*crhs31 - N[0]*crhs32 + N[0]*crhs43 + crhs36*crhs37 - crhs41*crhs44 - crhs42*crhs44;
            rhs[2]=-DN(0,0)*crhs39 - DN(0,1)*crhs44 - N[0]*crhs21 - N[0]*crhs22 - N[0]*crhs24;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs25 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs13 - N[1]*crhs14 - N[1]*crhs18 + N[1]*crhs2 + N[1]*crhs4 + N[1]*crhs9 - crhs33*crhs45 - crhs39*crhs46 - crhs39*crhs47;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs25 - DN(1,1)*stress[1] + N[1]*crhs26 - N[1]*crhs27 - N[1]*crhs28 - N[1]*crhs31 - N[1]*crhs32 + N[1]*crhs43 + crhs37*crhs45 - crhs44*crhs46 - crhs44*crhs47;
            rhs[5]=-DN(1,0)*crhs39 - DN(1,1)*crhs44 - N[1]*crhs21 - N[1]*crhs22 - N[1]*crhs24;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs25 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs13 - N[2]*crhs14 - N[2]*crhs18 + N[2]*crhs2 + N[2]*crhs4 + N[2]*crhs9 - crhs33*crhs48 - crhs39*crhs49 - crhs39*crhs50;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs25 - DN(2,1)*stress[1] + N[2]*crhs26 - N[2]*crhs27 - N[2]*crhs28 - N[2]*crhs31 - N[2]*crhs32 + N[2]*crhs43 + crhs37*crhs48 - crhs44*crhs49 - crhs44*crhs50;
            rhs[8]=-DN(2,0)*crhs39 - DN(2,1)*crhs44 - N[2]*crhs21 - N[2]*crhs22 - N[2]*crhs24;


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
    const array_1d<double,3>& omega_old = rData.AngularVelocityOld;
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
const double crhs3 =             omega[2] - omega_old[2];
const double crhs4 =             crhs3*gauss_pt_coord[1];
const double crhs5 =             omega[1] - omega_old[1];
const double crhs6 =             crhs5*gauss_pt_coord[2];
const double crhs7 =             crhs4 - crhs6;
const double crhs8 =             crhs1/dt;
const double crhs9 =             N[0]*crhs8;
const double crhs10 =             gauss_pt_coord[1]*omega[0];
const double crhs11 =             gauss_pt_coord[0]*omega[1];
const double crhs12 =             crhs10 - crhs11;
const double crhs13 =             gauss_pt_coord[0]*omega[2] - gauss_pt_coord[2]*omega[0];
const double crhs14 =             crhs1*(crhs12*omega[1] - crhs13*omega[2]);
const double crhs15 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs16 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs17 =             2.0*crhs1;
const double crhs18 =             crhs17*(crhs15*omega[1] - crhs16*omega[2]);
const double crhs19 =             crhs1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs20 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs21 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs22 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs23 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs24 =             crhs1*(crhs20*crhs21 + crhs22*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs23*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs25 =             crhs1*stab_c2*sqrt(pow(crhs21, 2) + pow(crhs22, 2) + pow(crhs23, 2));
const double crhs26 =             1.0/crhs1;
const double crhs27 =             crhs26*(crhs21*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crhs22*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]) + crhs23*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]));
const double crhs28 =             crhs26*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crhs29 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs30 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs31 =             crhs20 + crhs29 + crhs30;
const double crhs32 =             (crhs25*h/stab_c1 + mu)*(crhs27 + crhs28 + crhs31);
const double crhs33 =             DN(0,0)*p[0];
const double crhs34 =             DN(1,0)*p[1];
const double crhs35 =             DN(2,0)*p[2];
const double crhs36 =             DN(3,0)*p[3];
const double crhs37 =             crhs8*(-crhs4 + crhs6);
const double crhs38 =             crhs14 + crhs18 + crhs19 - crhs2 + crhs24 + crhs33 + crhs34 + crhs35 + crhs36 + crhs37;
const double crhs39 =             1.0/(crhs1*stab_c3*sqrt(pow(omega[0], 2) + pow(omega[1], 2) + pow(omega[2], 2)) + crhs25/h + crhs8*dyn_tau + mu*stab_c1/pow(h, 2));
const double crhs40 =             1.0*crhs39;
const double crhs41 =             crhs38*crhs40;
const double crhs42 =             crhs1*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs43 =             N[0]*crhs42;
const double crhs44 =             crhs1*(DN(0,0)*crhs21 + DN(0,1)*crhs22 + DN(0,2)*crhs23);
const double crhs45 =             DN(0,2)*p[0];
const double crhs46 =             DN(1,2)*p[1];
const double crhs47 =             DN(2,2)*p[2];
const double crhs48 =             DN(3,2)*p[3];
const double crhs49 =             crhs1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs50 =             omega[0] - omega_old[0];
const double crhs51 =             crhs50*gauss_pt_coord[1];
const double crhs52 =             crhs5*gauss_pt_coord[0];
const double crhs53 =             crhs8*(crhs51 - crhs52);
const double crhs54 =             gauss_pt_coord[2]*omega[1];
const double crhs55 =             gauss_pt_coord[1]*omega[2];
const double crhs56 =             crhs54 - crhs55;
const double crhs57 =             crhs1*(crhs13*omega[0] - crhs56*omega[1]);
const double crhs58 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs59 =             crhs17*(crhs16*omega[0] - crhs58*omega[1]);
const double crhs60 =             crhs1*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs61 =             crhs1*(crhs21*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs22*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs23*crhs29);
const double crhs62 =             DN(0,1)*p[0];
const double crhs63 =             DN(1,1)*p[1];
const double crhs64 =             DN(2,1)*p[2];
const double crhs65 =             DN(3,1)*p[3];
const double crhs66 =             crhs1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs67 =             crhs8*(crhs3*gauss_pt_coord[0] - crhs50*gauss_pt_coord[2]);
const double crhs68 =             crhs1*(-crhs12*omega[0] + crhs56*omega[2]);
const double crhs69 =             crhs58*omega[2];
const double crhs70 =             crhs15*omega[0];
const double crhs71 =             crhs17*(crhs69 - crhs70);
const double crhs72 =             crhs1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs73 =             crhs1*(crhs21*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs22*crhs30 + crhs23*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs74 =             -crhs62 - crhs63 - crhs64 - crhs65 + crhs66 - crhs67 - crhs68 - crhs71 - crhs72 - crhs73;
const double crhs75 =             -crhs74*omega[2] + omega[1]*(-crhs45 - crhs46 - crhs47 - crhs48 + crhs49 - crhs53 - crhs57 - crhs59 - crhs60 - crhs61);
const double crhs76 =             crhs17*crhs39;
const double crhs77 =             N[0]*crhs76;
const double crhs78 =             crhs1*(-omega[0]*(-crhs10 + crhs11) + omega[2]*(-crhs54 + crhs55));
const double crhs79 =             crhs17*(-crhs69 + crhs70);
const double crhs80 =             crhs40*(crhs62 + crhs63 + crhs64 + crhs65 - crhs66 + crhs67 + crhs68 + crhs71 + crhs72 + crhs73);
const double crhs81 =             crhs45 + crhs46 + crhs47 + crhs48 - crhs49 + crhs53 + crhs57 + crhs59 + crhs60 + crhs61;
const double crhs82 =             crhs38*omega[2] - crhs81*omega[0];
const double crhs83 =             -crhs51 + crhs52;
const double crhs84 =             crhs40*crhs81;
const double crhs85 =             crhs74*omega[0] - omega[1]*(-crhs14 - crhs18 - crhs19 + crhs2 - crhs24 - crhs33 - crhs34 - crhs35 - crhs36 - crhs37);
const double crhs86 =             N[1]*crhs8;
const double crhs87 =             N[1]*crhs42;
const double crhs88 =             crhs1*(DN(1,0)*crhs21 + DN(1,1)*crhs22 + DN(1,2)*crhs23);
const double crhs89 =             N[1]*crhs76;
const double crhs90 =             N[2]*crhs8;
const double crhs91 =             N[2]*crhs42;
const double crhs92 =             crhs1*(DN(2,0)*crhs21 + DN(2,1)*crhs22 + DN(2,2)*crhs23);
const double crhs93 =             N[2]*crhs76;
const double crhs94 =             N[3]*crhs8;
const double crhs95 =             N[3]*crhs42;
const double crhs96 =             crhs1*(DN(3,0)*crhs21 + DN(3,1)*crhs22 + DN(3,2)*crhs23);
const double crhs97 =             N[3]*crhs76;
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs32 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs14 - N[0]*crhs18 - N[0]*crhs19 + N[0]*crhs2 - N[0]*crhs24 - crhs41*crhs43 - crhs41*crhs44 + crhs7*crhs9 - crhs75*crhs77;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs32 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs66 - N[0]*crhs67 - N[0]*crhs72 - N[0]*crhs73 + N[0]*crhs78 + N[0]*crhs79 - crhs43*crhs80 - crhs44*crhs80 + crhs77*crhs82;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs32 - DN(0,2)*stress[2] + N[0]*crhs49 - N[0]*crhs57 - N[0]*crhs59 - N[0]*crhs60 - N[0]*crhs61 - crhs43*crhs84 - crhs44*crhs84 - crhs77*crhs85 + crhs83*crhs9;
            rhs[3]=-DN(0,0)*crhs41 - DN(0,1)*crhs80 - DN(0,2)*crhs84 - N[0]*crhs27 - N[0]*crhs28 - N[0]*crhs31;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs32 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs14 - N[1]*crhs18 - N[1]*crhs19 + N[1]*crhs2 - N[1]*crhs24 - crhs41*crhs87 - crhs41*crhs88 + crhs7*crhs86 - crhs75*crhs89;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs32 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs66 - N[1]*crhs67 - N[1]*crhs72 - N[1]*crhs73 + N[1]*crhs78 + N[1]*crhs79 - crhs80*crhs87 - crhs80*crhs88 + crhs82*crhs89;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs32 - DN(1,2)*stress[2] + N[1]*crhs49 - N[1]*crhs57 - N[1]*crhs59 - N[1]*crhs60 - N[1]*crhs61 + crhs83*crhs86 - crhs84*crhs87 - crhs84*crhs88 - crhs85*crhs89;
            rhs[7]=-DN(1,0)*crhs41 - DN(1,1)*crhs80 - DN(1,2)*crhs84 - N[1]*crhs27 - N[1]*crhs28 - N[1]*crhs31;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs32 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs14 - N[2]*crhs18 - N[2]*crhs19 + N[2]*crhs2 - N[2]*crhs24 - crhs41*crhs91 - crhs41*crhs92 + crhs7*crhs90 - crhs75*crhs93;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs32 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs66 - N[2]*crhs67 - N[2]*crhs72 - N[2]*crhs73 + N[2]*crhs78 + N[2]*crhs79 - crhs80*crhs91 - crhs80*crhs92 + crhs82*crhs93;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs32 - DN(2,2)*stress[2] + N[2]*crhs49 - N[2]*crhs57 - N[2]*crhs59 - N[2]*crhs60 - N[2]*crhs61 + crhs83*crhs90 - crhs84*crhs91 - crhs84*crhs92 - crhs85*crhs93;
            rhs[11]=-DN(2,0)*crhs41 - DN(2,1)*crhs80 - DN(2,2)*crhs84 - N[2]*crhs27 - N[2]*crhs28 - N[2]*crhs31;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs32 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs14 - N[3]*crhs18 - N[3]*crhs19 + N[3]*crhs2 - N[3]*crhs24 - crhs41*crhs95 - crhs41*crhs96 + crhs7*crhs94 - crhs75*crhs97;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs32 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs66 - N[3]*crhs67 - N[3]*crhs72 - N[3]*crhs73 + N[3]*crhs78 + N[3]*crhs79 - crhs80*crhs95 - crhs80*crhs96 + crhs82*crhs97;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs32 - DN(3,2)*stress[2] + N[3]*crhs49 - N[3]*crhs57 - N[3]*crhs59 - N[3]*crhs60 - N[3]*crhs61 + crhs83*crhs94 - crhs84*crhs95 - crhs84*crhs96 - crhs85*crhs97;
            rhs[15]=-DN(3,0)*crhs41 - DN(3,1)*crhs80 - DN(3,2)*crhs84 - N[3]*crhs27 - N[3]*crhs28 - N[3]*crhs31;


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