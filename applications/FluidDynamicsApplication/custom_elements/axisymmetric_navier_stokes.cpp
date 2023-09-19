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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int AxisymmetricNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Perform base fluid element check
    int out = BaseType::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Check that there are no negative y-coordinates (radius is always positive)
    const auto& r_geom = this->GetGeometry();
    for (const auto& r_node : r_geom) {
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

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

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

    // rData.lhs *= rData.Weight;
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

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = pow(DN(0,0), 2);
const double clhs1 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double clhs2 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double clhs3 = rho*stab_c2*sqrt(pow(clhs1, 2) + pow(clhs2, 2));
const double clhs4 = clhs3*h/stab_c1 + mu;
const double clhs5 = 2*clhs4;
const double clhs6 = DN(0,0)*clhs1;
const double clhs7 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double clhs8 = N[0]*clhs7 + clhs6;
const double clhs9 = N[0]*bdf0;
const double clhs10 = DN(0,1)*clhs2;
const double clhs11 = clhs10 + clhs6 + clhs9;
const double clhs12 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs13 = 1.0*clhs12;
const double clhs14 = clhs13*pow(rho, 2);
const double clhs15 = clhs11*clhs14;
const double clhs16 = clhs15*clhs8;
const double clhs17 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double clhs18 = N[0]*clhs17 + clhs10;
const double clhs19 = pow(N[0], 2);
const double clhs20 = bdf0*rho;
const double clhs21 = clhs6*rho;
const double clhs22 = clhs10*rho;
const double clhs23 = N[0]*clhs21 + N[0]*clhs22 + clhs19*clhs20;
const double clhs24 = pow(DN(0,1), 2);
const double clhs25 = 1.0/y;
const double clhs26 = N[0]*clhs25;
const double clhs27 = DN(0,1)*mu;
const double clhs28 = clhs26*clhs27;
const double clhs29 = DN(0,1) + clhs26;
const double clhs30 = clhs25*clhs27;
const double clhs31 = clhs13*rho;
const double clhs32 = clhs31*clhs8;
const double clhs33 = clhs30*clhs32;
const double clhs34 = DN(0,1)*clhs18;
const double clhs35 = clhs25*mu;
const double clhs36 = clhs31*clhs35;
const double clhs37 = DN(0,0)*clhs8;
const double clhs38 = clhs18*clhs31;
const double clhs39 = N[1]*bdf0;
const double clhs40 = DN(1,0)*clhs1;
const double clhs41 = DN(1,1)*clhs2;
const double clhs42 = clhs39 + clhs40 + clhs41;
const double clhs43 = clhs14*clhs42;
const double clhs44 = clhs43*clhs8;
const double clhs45 = DN(0,0)*DN(1,0);
const double clhs46 = clhs9*rho;
const double clhs47 = N[1]*clhs46;
const double clhs48 = clhs45*clhs5 + clhs47;
const double clhs49 = clhs40*rho;
const double clhs50 = clhs41*rho;
const double clhs51 = N[0]*clhs49 + N[0]*clhs50;
const double clhs52 = -DN(0,1)*DN(1,1)*mu;
const double clhs53 = DN(1,1)*clhs26;
const double clhs54 = clhs53*mu;
const double clhs55 = N[1]*clhs25;
const double clhs56 = DN(1,1) + clhs55;
const double clhs57 = DN(1,1)*clhs35;
const double clhs58 = clhs32*clhs57;
const double clhs59 = DN(1,1)*clhs18;
const double clhs60 = DN(1,0)*N[0];
const double clhs61 = DN(1,0)*clhs8;
const double clhs62 = N[2]*bdf0;
const double clhs63 = DN(2,0)*clhs1;
const double clhs64 = DN(2,1)*clhs2;
const double clhs65 = clhs62 + clhs63 + clhs64;
const double clhs66 = clhs14*clhs65;
const double clhs67 = clhs66*clhs8;
const double clhs68 = DN(0,0)*DN(2,0);
const double clhs69 = N[2]*clhs46;
const double clhs70 = clhs5*clhs68 + clhs69;
const double clhs71 = clhs63*rho;
const double clhs72 = clhs64*rho;
const double clhs73 = N[0]*clhs71 + N[0]*clhs72;
const double clhs74 = -DN(0,1)*DN(2,1)*mu;
const double clhs75 = DN(2,1)*clhs26;
const double clhs76 = clhs75*mu;
const double clhs77 = N[2]*clhs25;
const double clhs78 = DN(2,1) + clhs77;
const double clhs79 = DN(2,1)*clhs35;
const double clhs80 = clhs32*clhs79;
const double clhs81 = DN(2,1)*clhs18;
const double clhs82 = DN(2,0)*N[0];
const double clhs83 = DN(2,0)*clhs8;
const double clhs84 = mu/pow(y, 2);
const double clhs85 = N[0]*clhs84;
const double clhs86 = clhs21 + clhs22 - clhs30 + clhs46 + clhs85;
const double clhs87 = clhs39*rho;
const double clhs88 = N[1]*clhs84;
const double clhs89 = clhs49 + clhs50 - clhs57 + clhs87 + clhs88;
const double clhs90 = DN(1,1)*clhs27 + N[1]*clhs85 + clhs47;
const double clhs91 = N[2]*clhs84 + clhs62*rho + clhs71 + clhs72 - clhs79;
const double clhs92 = DN(2,1)*clhs27 + N[2]*clhs85 + clhs69;
const double clhs93 = DN(0,0)*clhs13;
const double clhs94 = clhs13*clhs86;
const double clhs95 = clhs13*clhs89;
const double clhs96 = DN(0,1)*y;
const double clhs97 = DN(1,1)*clhs96 + clhs45;
const double clhs98 = clhs31*clhs65;
const double clhs99 = clhs13*clhs91;
const double clhs100 = DN(2,1)*clhs96 + clhs68;
const double clhs101 = N[1]*clhs7 + clhs40;
const double clhs102 = clhs101*clhs15;
const double clhs103 = N[1]*clhs17 + clhs41;
const double clhs104 = N[1]*clhs21 + N[1]*clhs22;
const double clhs105 = clhs27*clhs55;
const double clhs106 = clhs101*clhs31;
const double clhs107 = clhs106*clhs30;
const double clhs108 = DN(0,1)*clhs103;
const double clhs109 = DN(0,0)*N[1];
const double clhs110 = DN(0,0)*clhs101;
const double clhs111 = clhs103*clhs31;
const double clhs112 = pow(DN(1,0), 2);
const double clhs113 = clhs101*clhs43;
const double clhs114 = pow(N[1], 2);
const double clhs115 = N[1]*clhs49 + N[1]*clhs50 + clhs114*clhs20;
const double clhs116 = pow(DN(1,1), 2);
const double clhs117 = DN(1,1)*clhs55;
const double clhs118 = clhs117*mu;
const double clhs119 = clhs106*clhs57;
const double clhs120 = DN(1,1)*clhs103;
const double clhs121 = DN(1,0)*clhs101;
const double clhs122 = clhs101*clhs66;
const double clhs123 = DN(1,0)*DN(2,0);
const double clhs124 = N[2]*clhs87;
const double clhs125 = clhs123*clhs5 + clhs124;
const double clhs126 = N[1]*clhs71 + N[1]*clhs72;
const double clhs127 = -DN(1,1)*DN(2,1)*mu;
const double clhs128 = DN(2,1)*clhs55;
const double clhs129 = clhs128*mu;
const double clhs130 = clhs106*clhs79;
const double clhs131 = DN(2,1)*clhs103;
const double clhs132 = DN(2,0)*N[1];
const double clhs133 = DN(2,0)*clhs101;
const double clhs134 = DN(1,1)*DN(2,1);
const double clhs135 = N[2]*clhs88 + clhs124 + clhs134*mu;
const double clhs136 = DN(1,0)*clhs13;
const double clhs137 = clhs123 + clhs134*y;
const double clhs138 = N[2]*clhs7 + clhs63;
const double clhs139 = clhs138*clhs15;
const double clhs140 = N[2]*clhs17 + clhs64;
const double clhs141 = N[2]*clhs21 + N[2]*clhs22;
const double clhs142 = clhs27*clhs77;
const double clhs143 = clhs138*clhs31;
const double clhs144 = clhs143*clhs30;
const double clhs145 = DN(0,1)*clhs140;
const double clhs146 = DN(0,0)*N[2];
const double clhs147 = DN(0,0)*clhs138;
const double clhs148 = clhs140*clhs31;
const double clhs149 = clhs138*clhs43;
const double clhs150 = N[2]*clhs49 + N[2]*clhs50;
const double clhs151 = DN(1,1)*clhs77;
const double clhs152 = clhs151*mu;
const double clhs153 = clhs143*clhs57;
const double clhs154 = DN(1,1)*clhs140;
const double clhs155 = DN(1,0)*N[2];
const double clhs156 = DN(1,0)*clhs138;
const double clhs157 = pow(DN(2,0), 2);
const double clhs158 = clhs138*clhs66;
const double clhs159 = pow(N[2], 2);
const double clhs160 = N[2]*clhs71 + N[2]*clhs72 + clhs159*clhs20;
const double clhs161 = pow(DN(2,1), 2);
const double clhs162 = DN(2,1)*clhs77;
const double clhs163 = clhs162*mu;
const double clhs164 = clhs143*clhs79;
const double clhs165 = DN(2,1)*clhs140;
const double clhs166 = DN(2,0)*clhs138;
const double clhs167 = DN(2,0)*clhs13;
lhs(0,0)=clhs0*clhs5 + clhs15*clhs18 + clhs16 + clhs23;
lhs(0,1)=2*DN(0,0)*clhs29*clhs4 + clhs24*mu - clhs28 - clhs33 - clhs34*clhs36;
lhs(0,2)=DN(0,0)*N[0] + DN(0,0)*clhs38 + DN(0,1)*N[0] + clhs31*clhs37;
lhs(0,3)=clhs18*clhs43 + clhs44 + clhs48 + clhs51;
lhs(0,4)=2*DN(0,0)*clhs4*clhs56 - clhs36*clhs59 - clhs52 - clhs54 - clhs58;
lhs(0,5)=DN(1,0)*clhs38 + DN(1,1)*N[0] + clhs31*clhs61 + clhs60;
lhs(0,6)=clhs18*clhs66 + clhs67 + clhs70 + clhs73;
lhs(0,7)=2*DN(0,0)*clhs4*clhs78 - clhs36*clhs81 - clhs74 - clhs76 - clhs80;
lhs(0,8)=DN(2,0)*clhs38 + DN(2,1)*N[0] + clhs31*clhs83 + clhs82;
lhs(1,0)=clhs16;
lhs(1,1)=clhs19*clhs84 + clhs23 + clhs24*mu - clhs28 - clhs33 + clhs38*clhs86;
lhs(1,2)=clhs31*(clhs34 + clhs37);
lhs(1,3)=clhs44;
lhs(1,4)=clhs38*clhs89 + clhs51 - clhs54 - clhs58 + clhs90;
lhs(1,5)=clhs31*(clhs59 + clhs61);
lhs(1,6)=clhs67;
lhs(1,7)=clhs38*clhs91 + clhs73 - clhs76 - clhs80 + clhs92;
lhs(1,8)=clhs31*(clhs81 + clhs83);
lhs(2,0)=DN(0,0)*(-N[0] + 1.0*clhs11*clhs12*rho);
lhs(2,1)=1.0*DN(0,1)*clhs12*clhs86*y - N[0]*clhs29 - clhs26*clhs94 - clhs30*clhs93;
lhs(2,2)=clhs13*(-DN(0,1)*clhs26 + clhs0 + clhs24*y);
lhs(2,3)=DN(0,0)*clhs31*clhs42 - clhs60;
lhs(2,4)=1.0*DN(0,1)*clhs12*clhs89*y - N[0]*clhs56 - clhs26*clhs95 - clhs57*clhs93;
lhs(2,5)=clhs13*(-clhs53 + clhs97);
lhs(2,6)=DN(0,0)*clhs98 - clhs82;
lhs(2,7)=1.0*DN(0,1)*clhs12*clhs91*y - N[0]*clhs78 - clhs26*clhs99 - clhs79*clhs93;
lhs(2,8)=clhs13*(clhs100 - clhs75);
lhs(3,0)=clhs102 + clhs103*clhs15 + clhs104 + clhs48;
lhs(3,1)=2*DN(1,0)*clhs29*clhs4 - clhs105 - clhs107 - clhs108*clhs36 - clhs52;
lhs(3,2)=DN(0,0)*clhs111 + DN(0,1)*N[1] + clhs109 + clhs110*clhs31;
lhs(3,3)=clhs103*clhs43 + clhs112*clhs5 + clhs113 + clhs115;
lhs(3,4)=2*DN(1,0)*clhs4*clhs56 + clhs116*mu - clhs118 - clhs119 - clhs120*clhs36;
lhs(3,5)=DN(1,0)*N[1] + DN(1,0)*clhs111 + DN(1,1)*N[1] + clhs121*clhs31;
lhs(3,6)=clhs103*clhs66 + clhs122 + clhs125 + clhs126;
lhs(3,7)=2*DN(1,0)*clhs4*clhs78 - clhs127 - clhs129 - clhs130 - clhs131*clhs36;
lhs(3,8)=DN(2,0)*clhs111 + DN(2,1)*N[1] + clhs132 + clhs133*clhs31;
lhs(4,0)=clhs102;
lhs(4,1)=clhs104 - clhs105 - clhs107 + clhs111*clhs86 + clhs90;
lhs(4,2)=clhs31*(clhs108 + clhs110);
lhs(4,3)=clhs113;
lhs(4,4)=clhs111*clhs89 + clhs114*clhs84 + clhs115 + clhs116*mu - clhs118 - clhs119;
lhs(4,5)=clhs31*(clhs120 + clhs121);
lhs(4,6)=clhs122;
lhs(4,7)=clhs111*clhs91 + clhs126 - clhs129 - clhs130 + clhs135;
lhs(4,8)=clhs31*(clhs131 + clhs133);
lhs(5,0)=1.0*DN(1,0)*clhs11*clhs12*rho - clhs109;
lhs(5,1)=1.0*DN(1,1)*clhs12*clhs86*y - N[1]*clhs29 - clhs136*clhs30 - clhs55*clhs94;
lhs(5,2)=clhs13*(-DN(0,1)*clhs55 + clhs97);
lhs(5,3)=DN(1,0)*(-N[1] + 1.0*clhs12*clhs42*rho);
lhs(5,4)=1.0*DN(1,1)*clhs12*clhs89*y - N[1]*clhs56 - clhs136*clhs57 - clhs55*clhs95;
lhs(5,5)=clhs13*(clhs112 + clhs116*y - clhs117);
lhs(5,6)=DN(1,0)*clhs98 - clhs132;
lhs(5,7)=1.0*DN(1,1)*clhs12*clhs91*y - N[1]*clhs78 - clhs136*clhs79 - clhs55*clhs99;
lhs(5,8)=clhs13*(-clhs128 + clhs137);
lhs(6,0)=clhs139 + clhs140*clhs15 + clhs141 + clhs70;
lhs(6,1)=2*DN(2,0)*clhs29*clhs4 - clhs142 - clhs144 - clhs145*clhs36 - clhs74;
lhs(6,2)=DN(0,0)*clhs148 + DN(0,1)*N[2] + clhs146 + clhs147*clhs31;
lhs(6,3)=clhs125 + clhs140*clhs43 + clhs149 + clhs150;
lhs(6,4)=2*DN(2,0)*clhs4*clhs56 - clhs127 - clhs152 - clhs153 - clhs154*clhs36;
lhs(6,5)=DN(1,0)*clhs148 + DN(1,1)*N[2] + clhs155 + clhs156*clhs31;
lhs(6,6)=clhs140*clhs66 + clhs157*clhs5 + clhs158 + clhs160;
lhs(6,7)=2*DN(2,0)*clhs4*clhs78 + clhs161*mu - clhs163 - clhs164 - clhs165*clhs36;
lhs(6,8)=DN(2,0)*N[2] + DN(2,0)*clhs148 + DN(2,1)*N[2] + clhs166*clhs31;
lhs(7,0)=clhs139;
lhs(7,1)=clhs141 - clhs142 - clhs144 + clhs148*clhs86 + clhs92;
lhs(7,2)=clhs31*(clhs145 + clhs147);
lhs(7,3)=clhs149;
lhs(7,4)=clhs135 + clhs148*clhs89 + clhs150 - clhs152 - clhs153;
lhs(7,5)=clhs31*(clhs154 + clhs156);
lhs(7,6)=clhs158;
lhs(7,7)=clhs148*clhs91 + clhs159*clhs84 + clhs160 + clhs161*mu - clhs163 - clhs164;
lhs(7,8)=clhs31*(clhs165 + clhs166);
lhs(8,0)=1.0*DN(2,0)*clhs11*clhs12*rho - clhs146;
lhs(8,1)=1.0*DN(2,1)*clhs12*clhs86*y - N[2]*clhs29 - clhs167*clhs30 - clhs77*clhs94;
lhs(8,2)=clhs13*(-DN(0,1)*clhs77 + clhs100);
lhs(8,3)=1.0*DN(2,0)*clhs12*clhs42*rho - clhs155;
lhs(8,4)=1.0*DN(2,1)*clhs12*clhs89*y - N[2]*clhs56 - clhs167*clhs57 - clhs77*clhs95;
lhs(8,5)=clhs13*(clhs137 - clhs151);
lhs(8,6)=DN(2,0)*(-N[2] + 1.0*clhs12*clhs65*rho);
lhs(8,7)=1.0*DN(2,1)*clhs12*clhs91*y - N[2]*clhs78 - clhs167*clhs79 - clhs77*clhs99;
lhs(8,8)=clhs13*(clhs157 + clhs161*y - clhs162);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
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

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = pow(DN(0,0), 2);
const double clhs1 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double clhs2 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double clhs3 = rho*stab_c2*sqrt(pow(clhs1, 2) + pow(clhs2, 2));
const double clhs4 = clhs3*h/stab_c1 + mu;
const double clhs5 = 2*clhs4;
const double clhs6 = DN(0,0)*clhs1;
const double clhs7 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double clhs8 = N[0]*clhs7 + clhs6;
const double clhs9 = N[0]*bdf0;
const double clhs10 = DN(0,1)*clhs2;
const double clhs11 = clhs10 + clhs6 + clhs9;
const double clhs12 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs13 = 1.0*clhs12;
const double clhs14 = clhs13*pow(rho, 2);
const double clhs15 = clhs11*clhs14;
const double clhs16 = clhs15*clhs8;
const double clhs17 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double clhs18 = N[0]*clhs17 + clhs10;
const double clhs19 = pow(N[0], 2);
const double clhs20 = bdf0*rho;
const double clhs21 = clhs6*rho;
const double clhs22 = clhs10*rho;
const double clhs23 = N[0]*clhs21 + N[0]*clhs22 + clhs19*clhs20;
const double clhs24 = pow(DN(0,1), 2);
const double clhs25 = 1.0/y;
const double clhs26 = N[0]*clhs25;
const double clhs27 = DN(0,1)*mu;
const double clhs28 = clhs26*clhs27;
const double clhs29 = DN(0,1) + clhs26;
const double clhs30 = clhs25*clhs27;
const double clhs31 = clhs13*rho;
const double clhs32 = clhs31*clhs8;
const double clhs33 = clhs30*clhs32;
const double clhs34 = DN(0,1)*clhs18;
const double clhs35 = clhs25*mu;
const double clhs36 = clhs31*clhs35;
const double clhs37 = DN(0,0)*clhs8;
const double clhs38 = clhs18*clhs31;
const double clhs39 = N[1]*bdf0;
const double clhs40 = DN(1,0)*clhs1;
const double clhs41 = DN(1,1)*clhs2;
const double clhs42 = clhs39 + clhs40 + clhs41;
const double clhs43 = clhs14*clhs42;
const double clhs44 = clhs43*clhs8;
const double clhs45 = DN(0,0)*DN(1,0);
const double clhs46 = clhs9*rho;
const double clhs47 = N[1]*clhs46;
const double clhs48 = clhs45*clhs5 + clhs47;
const double clhs49 = clhs40*rho;
const double clhs50 = clhs41*rho;
const double clhs51 = N[0]*clhs49 + N[0]*clhs50;
const double clhs52 = -DN(0,1)*DN(1,1)*mu;
const double clhs53 = DN(1,1)*clhs26;
const double clhs54 = clhs53*mu;
const double clhs55 = N[1]*clhs25;
const double clhs56 = DN(1,1) + clhs55;
const double clhs57 = DN(1,1)*clhs35;
const double clhs58 = clhs32*clhs57;
const double clhs59 = DN(1,1)*clhs18;
const double clhs60 = DN(1,0)*N[0];
const double clhs61 = DN(1,0)*clhs8;
const double clhs62 = N[2]*bdf0;
const double clhs63 = DN(2,0)*clhs1;
const double clhs64 = DN(2,1)*clhs2;
const double clhs65 = clhs62 + clhs63 + clhs64;
const double clhs66 = clhs14*clhs65;
const double clhs67 = clhs66*clhs8;
const double clhs68 = DN(0,0)*DN(2,0);
const double clhs69 = N[2]*clhs46;
const double clhs70 = clhs5*clhs68 + clhs69;
const double clhs71 = clhs63*rho;
const double clhs72 = clhs64*rho;
const double clhs73 = N[0]*clhs71 + N[0]*clhs72;
const double clhs74 = -DN(0,1)*DN(2,1)*mu;
const double clhs75 = DN(2,1)*clhs26;
const double clhs76 = clhs75*mu;
const double clhs77 = N[2]*clhs25;
const double clhs78 = DN(2,1) + clhs77;
const double clhs79 = DN(2,1)*clhs35;
const double clhs80 = clhs32*clhs79;
const double clhs81 = DN(2,1)*clhs18;
const double clhs82 = DN(2,0)*N[0];
const double clhs83 = DN(2,0)*clhs8;
const double clhs84 = N[3]*bdf0;
const double clhs85 = DN(3,0)*clhs1;
const double clhs86 = DN(3,1)*clhs2;
const double clhs87 = clhs84 + clhs85 + clhs86;
const double clhs88 = clhs14*clhs87;
const double clhs89 = clhs8*clhs88;
const double clhs90 = DN(0,0)*DN(3,0);
const double clhs91 = N[3]*clhs46;
const double clhs92 = clhs5*clhs90 + clhs91;
const double clhs93 = clhs85*rho;
const double clhs94 = clhs86*rho;
const double clhs95 = N[0]*clhs93 + N[0]*clhs94;
const double clhs96 = -DN(0,1)*DN(3,1)*mu;
const double clhs97 = DN(3,1)*clhs26;
const double clhs98 = clhs97*mu;
const double clhs99 = N[3]*clhs25;
const double clhs100 = DN(3,1) + clhs99;
const double clhs101 = DN(3,1)*clhs35;
const double clhs102 = clhs101*clhs32;
const double clhs103 = DN(3,1)*clhs18;
const double clhs104 = DN(3,0)*N[0];
const double clhs105 = DN(3,0)*clhs8;
const double clhs106 = mu/pow(y, 2);
const double clhs107 = N[0]*clhs106;
const double clhs108 = clhs107 + clhs21 + clhs22 - clhs30 + clhs46;
const double clhs109 = clhs39*rho;
const double clhs110 = N[1]*clhs106;
const double clhs111 = clhs109 + clhs110 + clhs49 + clhs50 - clhs57;
const double clhs112 = DN(1,1)*clhs27 + N[1]*clhs107 + clhs47;
const double clhs113 = clhs62*rho;
const double clhs114 = N[2]*clhs106;
const double clhs115 = clhs113 + clhs114 + clhs71 + clhs72 - clhs79;
const double clhs116 = DN(2,1)*clhs27 + N[2]*clhs107 + clhs69;
const double clhs117 = N[3]*clhs106 - clhs101 + clhs84*rho + clhs93 + clhs94;
const double clhs118 = DN(3,1)*clhs27 + N[3]*clhs107 + clhs91;
const double clhs119 = DN(0,0)*clhs13;
const double clhs120 = clhs108*clhs13;
const double clhs121 = clhs111*clhs13;
const double clhs122 = DN(0,1)*y;
const double clhs123 = DN(1,1)*clhs122 + clhs45;
const double clhs124 = clhs31*clhs65;
const double clhs125 = clhs115*clhs13;
const double clhs126 = DN(2,1)*clhs122 + clhs68;
const double clhs127 = clhs31*clhs87;
const double clhs128 = clhs117*clhs13;
const double clhs129 = DN(3,1)*clhs122 + clhs90;
const double clhs130 = N[1]*clhs7 + clhs40;
const double clhs131 = clhs130*clhs15;
const double clhs132 = N[1]*clhs17 + clhs41;
const double clhs133 = N[1]*clhs21 + N[1]*clhs22;
const double clhs134 = clhs27*clhs55;
const double clhs135 = clhs130*clhs31;
const double clhs136 = clhs135*clhs30;
const double clhs137 = DN(0,1)*clhs132;
const double clhs138 = DN(0,0)*N[1];
const double clhs139 = DN(0,0)*clhs130;
const double clhs140 = clhs132*clhs31;
const double clhs141 = pow(DN(1,0), 2);
const double clhs142 = clhs130*clhs43;
const double clhs143 = pow(N[1], 2);
const double clhs144 = N[1]*clhs49 + N[1]*clhs50 + clhs143*clhs20;
const double clhs145 = pow(DN(1,1), 2);
const double clhs146 = DN(1,1)*mu;
const double clhs147 = clhs146*clhs55;
const double clhs148 = clhs135*clhs57;
const double clhs149 = DN(1,1)*clhs132;
const double clhs150 = DN(1,0)*clhs130;
const double clhs151 = clhs130*clhs66;
const double clhs152 = DN(1,0)*DN(2,0);
const double clhs153 = N[2]*clhs109;
const double clhs154 = clhs152*clhs5 + clhs153;
const double clhs155 = N[1]*clhs71 + N[1]*clhs72;
const double clhs156 = -DN(1,1)*DN(2,1)*mu;
const double clhs157 = DN(2,1)*clhs55;
const double clhs158 = clhs157*mu;
const double clhs159 = clhs135*clhs79;
const double clhs160 = DN(2,1)*clhs132;
const double clhs161 = DN(2,0)*N[1];
const double clhs162 = DN(2,0)*clhs130;
const double clhs163 = clhs130*clhs88;
const double clhs164 = DN(1,0)*DN(3,0);
const double clhs165 = N[3]*clhs109;
const double clhs166 = clhs164*clhs5 + clhs165;
const double clhs167 = N[1]*clhs93 + N[1]*clhs94;
const double clhs168 = -DN(1,1)*DN(3,1)*mu;
const double clhs169 = DN(3,1)*clhs55;
const double clhs170 = clhs169*mu;
const double clhs171 = clhs101*clhs135;
const double clhs172 = DN(3,1)*clhs132;
const double clhs173 = DN(3,0)*N[1];
const double clhs174 = DN(3,0)*clhs130;
const double clhs175 = DN(2,1)*clhs146 + N[2]*clhs110 + clhs153;
const double clhs176 = DN(3,1)*clhs146 + N[3]*clhs110 + clhs165;
const double clhs177 = DN(1,0)*clhs13;
const double clhs178 = DN(1,1)*y;
const double clhs179 = DN(2,1)*clhs178 + clhs152;
const double clhs180 = DN(3,1)*clhs178 + clhs164;
const double clhs181 = N[2]*clhs7 + clhs63;
const double clhs182 = clhs15*clhs181;
const double clhs183 = N[2]*clhs17 + clhs64;
const double clhs184 = N[2]*clhs21 + N[2]*clhs22;
const double clhs185 = clhs27*clhs77;
const double clhs186 = clhs181*clhs31;
const double clhs187 = clhs186*clhs30;
const double clhs188 = DN(0,1)*clhs183;
const double clhs189 = DN(0,0)*N[2];
const double clhs190 = DN(0,0)*clhs181;
const double clhs191 = clhs183*clhs31;
const double clhs192 = clhs181*clhs43;
const double clhs193 = N[2]*clhs49 + N[2]*clhs50;
const double clhs194 = clhs146*clhs77;
const double clhs195 = clhs186*clhs57;
const double clhs196 = DN(1,1)*clhs183;
const double clhs197 = DN(1,0)*N[2];
const double clhs198 = DN(1,0)*clhs181;
const double clhs199 = pow(DN(2,0), 2);
const double clhs200 = clhs181*clhs66;
const double clhs201 = pow(N[2], 2);
const double clhs202 = N[2]*clhs71 + N[2]*clhs72 + clhs20*clhs201;
const double clhs203 = pow(DN(2,1), 2);
const double clhs204 = DN(2,1)*clhs77;
const double clhs205 = clhs204*mu;
const double clhs206 = clhs186*clhs79;
const double clhs207 = DN(2,1)*clhs183;
const double clhs208 = DN(2,0)*clhs181;
const double clhs209 = clhs181*clhs88;
const double clhs210 = DN(2,0)*DN(3,0);
const double clhs211 = N[3]*clhs113;
const double clhs212 = clhs210*clhs5 + clhs211;
const double clhs213 = N[2]*clhs93 + N[2]*clhs94;
const double clhs214 = -DN(2,1)*DN(3,1)*mu;
const double clhs215 = DN(3,1)*clhs77;
const double clhs216 = clhs215*mu;
const double clhs217 = clhs101*clhs186;
const double clhs218 = DN(3,1)*clhs183;
const double clhs219 = DN(3,0)*N[2];
const double clhs220 = DN(3,0)*clhs181;
const double clhs221 = DN(2,1)*DN(3,1);
const double clhs222 = N[3]*clhs114 + clhs211 + clhs221*mu;
const double clhs223 = DN(2,0)*clhs13;
const double clhs224 = clhs210 + clhs221*y;
const double clhs225 = N[3]*clhs7 + clhs85;
const double clhs226 = clhs15*clhs225;
const double clhs227 = N[3]*clhs17 + clhs86;
const double clhs228 = N[3]*clhs21 + N[3]*clhs22;
const double clhs229 = clhs27*clhs99;
const double clhs230 = clhs225*clhs31;
const double clhs231 = clhs230*clhs30;
const double clhs232 = DN(0,1)*clhs227;
const double clhs233 = DN(0,0)*N[3];
const double clhs234 = DN(0,0)*clhs225;
const double clhs235 = clhs227*clhs31;
const double clhs236 = clhs225*clhs43;
const double clhs237 = N[3]*clhs49 + N[3]*clhs50;
const double clhs238 = clhs146*clhs99;
const double clhs239 = clhs230*clhs57;
const double clhs240 = DN(1,1)*clhs227;
const double clhs241 = DN(1,0)*N[3];
const double clhs242 = DN(1,0)*clhs225;
const double clhs243 = clhs225*clhs66;
const double clhs244 = N[3]*clhs71 + N[3]*clhs72;
const double clhs245 = DN(2,1)*clhs99;
const double clhs246 = clhs245*mu;
const double clhs247 = clhs230*clhs79;
const double clhs248 = DN(2,1)*clhs227;
const double clhs249 = DN(2,0)*N[3];
const double clhs250 = DN(2,0)*clhs225;
const double clhs251 = pow(DN(3,0), 2);
const double clhs252 = clhs225*clhs88;
const double clhs253 = pow(N[3], 2);
const double clhs254 = N[3]*clhs93 + N[3]*clhs94 + clhs20*clhs253;
const double clhs255 = pow(DN(3,1), 2);
const double clhs256 = DN(3,1)*clhs99;
const double clhs257 = clhs256*mu;
const double clhs258 = clhs101*clhs230;
const double clhs259 = DN(3,1)*clhs227;
const double clhs260 = DN(3,0)*clhs225;
const double clhs261 = DN(3,0)*clhs13;
lhs(0,0)=clhs0*clhs5 + clhs15*clhs18 + clhs16 + clhs23;
lhs(0,1)=2*DN(0,0)*clhs29*clhs4 + clhs24*mu - clhs28 - clhs33 - clhs34*clhs36;
lhs(0,2)=DN(0,0)*N[0] + DN(0,0)*clhs38 + DN(0,1)*N[0] + clhs31*clhs37;
lhs(0,3)=clhs18*clhs43 + clhs44 + clhs48 + clhs51;
lhs(0,4)=2*DN(0,0)*clhs4*clhs56 - clhs36*clhs59 - clhs52 - clhs54 - clhs58;
lhs(0,5)=DN(1,0)*clhs38 + DN(1,1)*N[0] + clhs31*clhs61 + clhs60;
lhs(0,6)=clhs18*clhs66 + clhs67 + clhs70 + clhs73;
lhs(0,7)=2*DN(0,0)*clhs4*clhs78 - clhs36*clhs81 - clhs74 - clhs76 - clhs80;
lhs(0,8)=DN(2,0)*clhs38 + DN(2,1)*N[0] + clhs31*clhs83 + clhs82;
lhs(0,9)=clhs18*clhs88 + clhs89 + clhs92 + clhs95;
lhs(0,10)=2*DN(0,0)*clhs100*clhs4 - clhs102 - clhs103*clhs36 - clhs96 - clhs98;
lhs(0,11)=DN(3,0)*clhs38 + DN(3,1)*N[0] + clhs104 + clhs105*clhs31;
lhs(1,0)=clhs16;
lhs(1,1)=clhs106*clhs19 + clhs108*clhs38 + clhs23 + clhs24*mu - clhs28 - clhs33;
lhs(1,2)=clhs31*(clhs34 + clhs37);
lhs(1,3)=clhs44;
lhs(1,4)=clhs111*clhs38 + clhs112 + clhs51 - clhs54 - clhs58;
lhs(1,5)=clhs31*(clhs59 + clhs61);
lhs(1,6)=clhs67;
lhs(1,7)=clhs115*clhs38 + clhs116 + clhs73 - clhs76 - clhs80;
lhs(1,8)=clhs31*(clhs81 + clhs83);
lhs(1,9)=clhs89;
lhs(1,10)=-clhs102 + clhs117*clhs38 + clhs118 + clhs95 - clhs98;
lhs(1,11)=clhs31*(clhs103 + clhs105);
lhs(2,0)=DN(0,0)*(-N[0] + 1.0*clhs11*clhs12*rho);
lhs(2,1)=1.0*DN(0,1)*clhs108*clhs12*y - N[0]*clhs29 - clhs119*clhs30 - clhs120*clhs26;
lhs(2,2)=clhs13*(-DN(0,1)*clhs26 + clhs0 + clhs24*y);
lhs(2,3)=DN(0,0)*clhs31*clhs42 - clhs60;
lhs(2,4)=1.0*DN(0,1)*clhs111*clhs12*y - N[0]*clhs56 - clhs119*clhs57 - clhs121*clhs26;
lhs(2,5)=clhs13*(clhs123 - clhs53);
lhs(2,6)=DN(0,0)*clhs124 - clhs82;
lhs(2,7)=1.0*DN(0,1)*clhs115*clhs12*y - N[0]*clhs78 - clhs119*clhs79 - clhs125*clhs26;
lhs(2,8)=clhs13*(clhs126 - clhs75);
lhs(2,9)=DN(0,0)*clhs127 - clhs104;
lhs(2,10)=1.0*DN(0,1)*clhs117*clhs12*y - N[0]*clhs100 - clhs101*clhs119 - clhs128*clhs26;
lhs(2,11)=clhs13*(clhs129 - clhs97);
lhs(3,0)=clhs131 + clhs132*clhs15 + clhs133 + clhs48;
lhs(3,1)=2*DN(1,0)*clhs29*clhs4 - clhs134 - clhs136 - clhs137*clhs36 - clhs52;
lhs(3,2)=DN(0,0)*clhs140 + DN(0,1)*N[1] + clhs138 + clhs139*clhs31;
lhs(3,3)=clhs132*clhs43 + clhs141*clhs5 + clhs142 + clhs144;
lhs(3,4)=2*DN(1,0)*clhs4*clhs56 + clhs145*mu - clhs147 - clhs148 - clhs149*clhs36;
lhs(3,5)=DN(1,0)*N[1] + DN(1,0)*clhs140 + DN(1,1)*N[1] + clhs150*clhs31;
lhs(3,6)=clhs132*clhs66 + clhs151 + clhs154 + clhs155;
lhs(3,7)=2*DN(1,0)*clhs4*clhs78 - clhs156 - clhs158 - clhs159 - clhs160*clhs36;
lhs(3,8)=DN(2,0)*clhs140 + DN(2,1)*N[1] + clhs161 + clhs162*clhs31;
lhs(3,9)=clhs132*clhs88 + clhs163 + clhs166 + clhs167;
lhs(3,10)=2*DN(1,0)*clhs100*clhs4 - clhs168 - clhs170 - clhs171 - clhs172*clhs36;
lhs(3,11)=DN(3,0)*clhs140 + DN(3,1)*N[1] + clhs173 + clhs174*clhs31;
lhs(4,0)=clhs131;
lhs(4,1)=clhs108*clhs140 + clhs112 + clhs133 - clhs134 - clhs136;
lhs(4,2)=clhs31*(clhs137 + clhs139);
lhs(4,3)=clhs142;
lhs(4,4)=clhs106*clhs143 + clhs111*clhs140 + clhs144 + clhs145*mu - clhs147 - clhs148;
lhs(4,5)=clhs31*(clhs149 + clhs150);
lhs(4,6)=clhs151;
lhs(4,7)=clhs115*clhs140 + clhs155 - clhs158 - clhs159 + clhs175;
lhs(4,8)=clhs31*(clhs160 + clhs162);
lhs(4,9)=clhs163;
lhs(4,10)=clhs117*clhs140 + clhs167 - clhs170 - clhs171 + clhs176;
lhs(4,11)=clhs31*(clhs172 + clhs174);
lhs(5,0)=1.0*DN(1,0)*clhs11*clhs12*rho - clhs138;
lhs(5,1)=1.0*DN(1,1)*clhs108*clhs12*y - N[1]*clhs29 - clhs120*clhs55 - clhs177*clhs30;
lhs(5,2)=clhs13*(-DN(0,1)*clhs55 + clhs123);
lhs(5,3)=DN(1,0)*(-N[1] + 1.0*clhs12*clhs42*rho);
lhs(5,4)=1.0*DN(1,1)*clhs111*clhs12*y - N[1]*clhs56 - clhs121*clhs55 - clhs177*clhs57;
lhs(5,5)=clhs13*(-DN(1,1)*clhs55 + clhs141 + clhs145*y);
lhs(5,6)=DN(1,0)*clhs124 - clhs161;
lhs(5,7)=1.0*DN(1,1)*clhs115*clhs12*y - N[1]*clhs78 - clhs125*clhs55 - clhs177*clhs79;
lhs(5,8)=clhs13*(-clhs157 + clhs179);
lhs(5,9)=DN(1,0)*clhs127 - clhs173;
lhs(5,10)=1.0*DN(1,1)*clhs117*clhs12*y - N[1]*clhs100 - clhs101*clhs177 - clhs128*clhs55;
lhs(5,11)=clhs13*(-clhs169 + clhs180);
lhs(6,0)=clhs15*clhs183 + clhs182 + clhs184 + clhs70;
lhs(6,1)=2*DN(2,0)*clhs29*clhs4 - clhs185 - clhs187 - clhs188*clhs36 - clhs74;
lhs(6,2)=DN(0,0)*clhs191 + DN(0,1)*N[2] + clhs189 + clhs190*clhs31;
lhs(6,3)=clhs154 + clhs183*clhs43 + clhs192 + clhs193;
lhs(6,4)=2*DN(2,0)*clhs4*clhs56 - clhs156 - clhs194 - clhs195 - clhs196*clhs36;
lhs(6,5)=DN(1,0)*clhs191 + DN(1,1)*N[2] + clhs197 + clhs198*clhs31;
lhs(6,6)=clhs183*clhs66 + clhs199*clhs5 + clhs200 + clhs202;
lhs(6,7)=2*DN(2,0)*clhs4*clhs78 + clhs203*mu - clhs205 - clhs206 - clhs207*clhs36;
lhs(6,8)=DN(2,0)*N[2] + DN(2,0)*clhs191 + DN(2,1)*N[2] + clhs208*clhs31;
lhs(6,9)=clhs183*clhs88 + clhs209 + clhs212 + clhs213;
lhs(6,10)=2*DN(2,0)*clhs100*clhs4 - clhs214 - clhs216 - clhs217 - clhs218*clhs36;
lhs(6,11)=DN(3,0)*clhs191 + DN(3,1)*N[2] + clhs219 + clhs220*clhs31;
lhs(7,0)=clhs182;
lhs(7,1)=clhs108*clhs191 + clhs116 + clhs184 - clhs185 - clhs187;
lhs(7,2)=clhs31*(clhs188 + clhs190);
lhs(7,3)=clhs192;
lhs(7,4)=clhs111*clhs191 + clhs175 + clhs193 - clhs194 - clhs195;
lhs(7,5)=clhs31*(clhs196 + clhs198);
lhs(7,6)=clhs200;
lhs(7,7)=clhs106*clhs201 + clhs115*clhs191 + clhs202 + clhs203*mu - clhs205 - clhs206;
lhs(7,8)=clhs31*(clhs207 + clhs208);
lhs(7,9)=clhs209;
lhs(7,10)=clhs117*clhs191 + clhs213 - clhs216 - clhs217 + clhs222;
lhs(7,11)=clhs31*(clhs218 + clhs220);
lhs(8,0)=1.0*DN(2,0)*clhs11*clhs12*rho - clhs189;
lhs(8,1)=1.0*DN(2,1)*clhs108*clhs12*y - N[2]*clhs29 - clhs120*clhs77 - clhs223*clhs30;
lhs(8,2)=clhs13*(-DN(0,1)*clhs77 + clhs126);
lhs(8,3)=1.0*DN(2,0)*clhs12*clhs42*rho - clhs197;
lhs(8,4)=1.0*DN(2,1)*clhs111*clhs12*y - N[2]*clhs56 - clhs121*clhs77 - clhs223*clhs57;
lhs(8,5)=clhs13*(-DN(1,1)*clhs77 + clhs179);
lhs(8,6)=DN(2,0)*(-N[2] + 1.0*clhs12*clhs65*rho);
lhs(8,7)=1.0*DN(2,1)*clhs115*clhs12*y - N[2]*clhs78 - clhs125*clhs77 - clhs223*clhs79;
lhs(8,8)=clhs13*(clhs199 + clhs203*y - clhs204);
lhs(8,9)=DN(2,0)*clhs127 - clhs219;
lhs(8,10)=1.0*DN(2,1)*clhs117*clhs12*y - N[2]*clhs100 - clhs101*clhs223 - clhs128*clhs77;
lhs(8,11)=clhs13*(-clhs215 + clhs224);
lhs(9,0)=clhs15*clhs227 + clhs226 + clhs228 + clhs92;
lhs(9,1)=2*DN(3,0)*clhs29*clhs4 - clhs229 - clhs231 - clhs232*clhs36 - clhs96;
lhs(9,2)=DN(0,0)*clhs235 + DN(0,1)*N[3] + clhs233 + clhs234*clhs31;
lhs(9,3)=clhs166 + clhs227*clhs43 + clhs236 + clhs237;
lhs(9,4)=2*DN(3,0)*clhs4*clhs56 - clhs168 - clhs238 - clhs239 - clhs240*clhs36;
lhs(9,5)=DN(1,0)*clhs235 + DN(1,1)*N[3] + clhs241 + clhs242*clhs31;
lhs(9,6)=clhs212 + clhs227*clhs66 + clhs243 + clhs244;
lhs(9,7)=2*DN(3,0)*clhs4*clhs78 - clhs214 - clhs246 - clhs247 - clhs248*clhs36;
lhs(9,8)=DN(2,0)*clhs235 + DN(2,1)*N[3] + clhs249 + clhs250*clhs31;
lhs(9,9)=clhs227*clhs88 + clhs251*clhs5 + clhs252 + clhs254;
lhs(9,10)=2*DN(3,0)*clhs100*clhs4 + clhs255*mu - clhs257 - clhs258 - clhs259*clhs36;
lhs(9,11)=DN(3,0)*N[3] + DN(3,0)*clhs235 + DN(3,1)*N[3] + clhs260*clhs31;
lhs(10,0)=clhs226;
lhs(10,1)=clhs108*clhs235 + clhs118 + clhs228 - clhs229 - clhs231;
lhs(10,2)=clhs31*(clhs232 + clhs234);
lhs(10,3)=clhs236;
lhs(10,4)=clhs111*clhs235 + clhs176 + clhs237 - clhs238 - clhs239;
lhs(10,5)=clhs31*(clhs240 + clhs242);
lhs(10,6)=clhs243;
lhs(10,7)=clhs115*clhs235 + clhs222 + clhs244 - clhs246 - clhs247;
lhs(10,8)=clhs31*(clhs248 + clhs250);
lhs(10,9)=clhs252;
lhs(10,10)=clhs106*clhs253 + clhs117*clhs235 + clhs254 + clhs255*mu - clhs257 - clhs258;
lhs(10,11)=clhs31*(clhs259 + clhs260);
lhs(11,0)=1.0*DN(3,0)*clhs11*clhs12*rho - clhs233;
lhs(11,1)=1.0*DN(3,1)*clhs108*clhs12*y - N[3]*clhs29 - clhs120*clhs99 - clhs261*clhs30;
lhs(11,2)=clhs13*(-DN(0,1)*clhs99 + clhs129);
lhs(11,3)=1.0*DN(3,0)*clhs12*clhs42*rho - clhs241;
lhs(11,4)=1.0*DN(3,1)*clhs111*clhs12*y - N[3]*clhs56 - clhs121*clhs99 - clhs261*clhs57;
lhs(11,5)=clhs13*(-DN(1,1)*clhs99 + clhs180);
lhs(11,6)=1.0*DN(3,0)*clhs12*clhs65*rho - clhs249;
lhs(11,7)=1.0*DN(3,1)*clhs115*clhs12*y - N[3]*clhs78 - clhs125*clhs99 - clhs261*clhs79;
lhs(11,8)=clhs13*(clhs224 - clhs245);
lhs(11,9)=DN(3,0)*(-N[3] + 1.0*clhs12*clhs87*rho);
lhs(11,10)=1.0*DN(3,1)*clhs117*clhs12*y - N[3]*clhs100 - clhs101*clhs261 - clhs128*clhs99;
lhs(11,11)=clhs13*(clhs251 + clhs255*y - clhs256);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
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

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2];
const double crhs1 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2];
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
const double crhs3 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs4 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double crhs5 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double crhs6 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs8 = 1.0/y;
const double crhs9 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs10 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs11 = crhs10 + crhs8*crhs9;
const double crhs12 = 2*(crhs11 + crhs7)*(crhs6*h/stab_c1 + mu);
const double crhs13 = N[0]*crhs7;
const double crhs14 = crhs4*rho;
const double crhs15 = crhs5*rho;
const double crhs16 = crhs15*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
const double crhs17 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double crhs18 = DN(0,1)*crhs5 + N[0]*crhs17;
const double crhs19 = crhs10*mu;
const double crhs20 = -crhs19*crhs8;
const double crhs21 = 1.0/(crhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = crhs21*(crhs0 + crhs14*crhs7 + crhs16 + crhs2 + crhs20 - crhs3*rho);
const double crhs23 = crhs22*rho;
const double crhs24 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double crhs25 = DN(0,1)*crhs19 - N[0]*crhs10*crhs8*mu + crhs23*(DN(0,0)*crhs4 + N[0]*crhs24);
const double crhs26 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
const double crhs27 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs28 = crhs9*mu/pow(y, 2);
const double crhs29 = crhs14*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
const double crhs30 = crhs10*crhs15;
const double crhs31 = crhs21*(crhs1 + crhs20 + crhs26 - crhs27*rho + crhs28 + crhs29 + crhs30);
const double crhs32 = crhs31*rho;
const double crhs33 = crhs31*y;
const double crhs34 = crhs31*crhs8;
const double crhs35 = N[1]*crhs7;
const double crhs36 = DN(1,1)*crhs5 + N[1]*crhs17;
const double crhs37 = DN(1,1)*crhs19 - N[1]*crhs10*crhs8*mu + crhs23*(DN(1,0)*crhs4 + N[1]*crhs24);
const double crhs38 = N[2]*crhs7;
const double crhs39 = DN(2,1)*crhs5 + N[2]*crhs17;
const double crhs40 = DN(2,1)*crhs19 - N[2]*crhs10*crhs8*mu + crhs23*(DN(2,0)*crhs4 + N[2]*crhs24);
rhs[0]=-DN(0,0)*crhs12 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs16 - N[0]*crhs2 + N[0]*crhs3*rho - crhs13*crhs14 - crhs18*crhs23 - crhs25;
rhs[1]=-N[0]*crhs26 + N[0]*crhs27*rho - N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs18*crhs32 - crhs25;
rhs[2]=-DN(0,0)*crhs22 - DN(0,1)*crhs33 + N[0]*crhs11 + N[0]*crhs34 + crhs13;
rhs[3]=-DN(1,0)*crhs12 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs16 - N[1]*crhs2 + N[1]*crhs3*rho - crhs14*crhs35 - crhs23*crhs36 - crhs37;
rhs[4]=-N[1]*crhs26 + N[1]*crhs27*rho - N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs32*crhs36 - crhs37;
rhs[5]=-DN(1,0)*crhs22 - DN(1,1)*crhs33 + N[1]*crhs11 + N[1]*crhs34 + crhs35;
rhs[6]=-DN(2,0)*crhs12 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs16 - N[2]*crhs2 + N[2]*crhs3*rho - crhs14*crhs38 - crhs23*crhs39 - crhs40;
rhs[7]=-N[2]*crhs26 + N[2]*crhs27*rho - N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs32*crhs39 - crhs40;
rhs[8]=-DN(2,0)*crhs22 - DN(2,1)*crhs33 + N[2]*crhs11 + N[2]*crhs34 + crhs38;


    noalias(rRHS) += rData.Weight * rhs;
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

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3];
const double crhs1 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3];
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
const double crhs3 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs4 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double crhs5 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double crhs6 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs8 = 1.0/y;
const double crhs9 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs10 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs11 = crhs10 + crhs8*crhs9;
const double crhs12 = 2*(crhs11 + crhs7)*(crhs6*h/stab_c1 + mu);
const double crhs13 = N[0]*crhs7;
const double crhs14 = crhs4*rho;
const double crhs15 = crhs5*rho;
const double crhs16 = crhs15*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0));
const double crhs17 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double crhs18 = DN(0,1)*crhs5 + N[0]*crhs17;
const double crhs19 = crhs10*mu;
const double crhs20 = -crhs19*crhs8;
const double crhs21 = 1.0/(crhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = crhs21*(crhs0 + crhs14*crhs7 + crhs16 + crhs2 + crhs20 - crhs3*rho);
const double crhs23 = crhs22*rho;
const double crhs24 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double crhs25 = DN(0,1)*crhs19 - N[0]*crhs10*crhs8*mu + crhs23*(DN(0,0)*crhs4 + N[0]*crhs24);
const double crhs26 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
const double crhs27 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs28 = crhs9*mu/pow(y, 2);
const double crhs29 = crhs14*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1));
const double crhs30 = crhs10*crhs15;
const double crhs31 = crhs21*(crhs1 + crhs20 + crhs26 - crhs27*rho + crhs28 + crhs29 + crhs30);
const double crhs32 = crhs31*rho;
const double crhs33 = crhs31*y;
const double crhs34 = crhs31*crhs8;
const double crhs35 = N[1]*crhs7;
const double crhs36 = DN(1,1)*crhs5 + N[1]*crhs17;
const double crhs37 = DN(1,1)*crhs19 - N[1]*crhs10*crhs8*mu + crhs23*(DN(1,0)*crhs4 + N[1]*crhs24);
const double crhs38 = N[2]*crhs7;
const double crhs39 = DN(2,1)*crhs5 + N[2]*crhs17;
const double crhs40 = DN(2,1)*crhs19 - N[2]*crhs10*crhs8*mu + crhs23*(DN(2,0)*crhs4 + N[2]*crhs24);
const double crhs41 = N[3]*crhs7;
const double crhs42 = DN(3,1)*crhs5 + N[3]*crhs17;
const double crhs43 = DN(3,1)*crhs19 - N[3]*crhs10*crhs8*mu + crhs23*(DN(3,0)*crhs4 + N[3]*crhs24);
rhs[0]=-DN(0,0)*crhs12 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs16 - N[0]*crhs2 + N[0]*crhs3*rho - crhs13*crhs14 - crhs18*crhs23 - crhs25;
rhs[1]=-N[0]*crhs26 + N[0]*crhs27*rho - N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs18*crhs32 - crhs25;
rhs[2]=-DN(0,0)*crhs22 - DN(0,1)*crhs33 + N[0]*crhs11 + N[0]*crhs34 + crhs13;
rhs[3]=-DN(1,0)*crhs12 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs16 - N[1]*crhs2 + N[1]*crhs3*rho - crhs14*crhs35 - crhs23*crhs36 - crhs37;
rhs[4]=-N[1]*crhs26 + N[1]*crhs27*rho - N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs32*crhs36 - crhs37;
rhs[5]=-DN(1,0)*crhs22 - DN(1,1)*crhs33 + N[1]*crhs11 + N[1]*crhs34 + crhs35;
rhs[6]=-DN(2,0)*crhs12 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs16 - N[2]*crhs2 + N[2]*crhs3*rho - crhs14*crhs38 - crhs23*crhs39 - crhs40;
rhs[7]=-N[2]*crhs26 + N[2]*crhs27*rho - N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs32*crhs39 - crhs40;
rhs[8]=-DN(2,0)*crhs22 - DN(2,1)*crhs33 + N[2]*crhs11 + N[2]*crhs34 + crhs38;
rhs[9]=-DN(3,0)*crhs12 - N[3]*crhs0 - N[3]*crhs1 - N[3]*crhs16 - N[3]*crhs2 + N[3]*crhs3*rho - crhs14*crhs41 - crhs23*crhs42 - crhs43;
rhs[10]=-N[3]*crhs26 + N[3]*crhs27*rho - N[3]*crhs28 - N[3]*crhs29 - N[3]*crhs30 - crhs32*crhs42 - crhs43;
rhs[11]=-DN(3,0)*crhs22 - DN(3,1)*crhs33 + N[3]*crhs11 + N[3]*crhs34 + crhs41;


    noalias(rRHS) += rData.Weight * rhs;
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