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

    const double clhs0 = pow(DN(0,1), 2);
const double clhs1 = clhs0*mu;
const double clhs2 = pow(DN(0,0), 2);
const double clhs3 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double clhs4 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double clhs5 = rho*stab_c2*sqrt(pow(clhs3, 2) + pow(clhs4, 2));
const double clhs6 = clhs5*h/stab_c1 + mu;
const double clhs7 = DN(0,0)*clhs3;
const double clhs8 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double clhs9 = N[0]*clhs8 + clhs7;
const double clhs10 = clhs9*rho;
const double clhs11 = bdf0*rho;
const double clhs12 = N[0]*clhs11;
const double clhs13 = 1.0/y;
const double clhs14 = clhs13*mu;
const double clhs15 = DN(0,1)*clhs14;
const double clhs16 = clhs7*rho;
const double clhs17 = DN(0,1)*clhs4;
const double clhs18 = clhs17*rho;
const double clhs19 = clhs12 - clhs15 + clhs16 + clhs18;
const double clhs20 = 1.0/(clhs5/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs21 = 1.0*clhs20;
const double clhs22 = clhs19*clhs21;
const double clhs23 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double clhs24 = N[0]*clhs23 + clhs17;
const double clhs25 = clhs24*rho;
const double clhs26 = N[0]*clhs13;
const double clhs27 = clhs14*(DN(0,1) - clhs26);
const double clhs28 = pow(N[0], 2);
const double clhs29 = DN(0,1)*clhs26;
const double clhs30 = N[0]*clhs16 + N[0]*clhs18 + clhs11*clhs28 + clhs2*mu - clhs29*mu;
const double clhs31 = DN(0,1) + clhs26;
const double clhs32 = DN(0,0)*clhs6;
const double clhs33 = DN(0,1)*N[0];
const double clhs34 = DN(0,0)*clhs21;
const double clhs35 = N[1]*clhs11;
const double clhs36 = DN(1,1)*clhs14;
const double clhs37 = DN(1,0)*clhs3;
const double clhs38 = clhs37*rho;
const double clhs39 = DN(1,1)*clhs4;
const double clhs40 = clhs39*rho;
const double clhs41 = clhs35 - clhs36 + clhs38 + clhs40;
const double clhs42 = clhs21*clhs41;
const double clhs43 = DN(0,0)*DN(1,0);
const double clhs44 = DN(0,1)*DN(1,1);
const double clhs45 = clhs44*mu;
const double clhs46 = N[1]*clhs12 + clhs43*mu;
const double clhs47 = clhs43*clhs6 + clhs45 + clhs46;
const double clhs48 = DN(1,0)*N[0];
const double clhs49 = clhs3*rho;
const double clhs50 = DN(1,1)*N[0];
const double clhs51 = clhs4*rho;
const double clhs52 = DN(1,1)*clhs26;
const double clhs53 = clhs48*clhs49 + clhs50*clhs51 - clhs52*mu;
const double clhs54 = N[1]*clhs13;
const double clhs55 = DN(1,1) + clhs54;
const double clhs56 = DN(1,0)*clhs21;
const double clhs57 = N[2]*clhs11;
const double clhs58 = DN(2,1)*clhs14;
const double clhs59 = DN(2,0)*clhs3;
const double clhs60 = clhs59*rho;
const double clhs61 = DN(2,1)*clhs4;
const double clhs62 = clhs61*rho;
const double clhs63 = clhs57 - clhs58 + clhs60 + clhs62;
const double clhs64 = clhs21*clhs63;
const double clhs65 = DN(0,0)*DN(2,0);
const double clhs66 = DN(0,1)*DN(2,1);
const double clhs67 = clhs66*mu;
const double clhs68 = N[2]*clhs12 + clhs65*mu;
const double clhs69 = clhs6*clhs65 + clhs67 + clhs68;
const double clhs70 = DN(2,0)*N[0];
const double clhs71 = DN(2,1)*N[0];
const double clhs72 = DN(2,1)*clhs26;
const double clhs73 = clhs49*clhs70 + clhs51*clhs71 - clhs72*mu;
const double clhs74 = N[2]*clhs13;
const double clhs75 = DN(2,1) + clhs74;
const double clhs76 = DN(2,0)*clhs21;
const double clhs77 = DN(0,1)*clhs6;
const double clhs78 = mu/pow(y, 2);
const double clhs79 = N[0]*clhs78;
const double clhs80 = clhs19 + clhs79;
const double clhs81 = clhs21*clhs80;
const double clhs82 = clhs21*(clhs24*rho - clhs27 - clhs79 + clhs9*rho);
const double clhs83 = N[1]*clhs78;
const double clhs84 = clhs41 + clhs83;
const double clhs85 = clhs21*clhs84;
const double clhs86 = N[1]*clhs79 + 2*clhs45 + clhs46;
const double clhs87 = N[2]*clhs78;
const double clhs88 = clhs63 + clhs87;
const double clhs89 = clhs21*clhs88;
const double clhs90 = N[2]*clhs79 + 2*clhs67 + clhs68;
const double clhs91 = DN(0,1)*clhs21;
const double clhs92 = N[1]*clhs26;
const double clhs93 = clhs43 + clhs44;
const double clhs94 = N[2]*clhs26;
const double clhs95 = clhs65 + clhs66;
const double clhs96 = N[1]*clhs8 + clhs37;
const double clhs97 = clhs96*rho;
const double clhs98 = N[1]*clhs23 + clhs39;
const double clhs99 = clhs98*rho;
const double clhs100 = clhs14*(DN(1,1) - clhs54);
const double clhs101 = DN(0,1)*clhs54;
const double clhs102 = N[1]*clhs16 + N[1]*clhs18 - clhs101*mu;
const double clhs103 = DN(1,0)*clhs6;
const double clhs104 = DN(0,0)*N[1];
const double clhs105 = DN(0,1)*N[1];
const double clhs106 = pow(DN(1,1), 2);
const double clhs107 = clhs106*mu;
const double clhs108 = pow(DN(1,0), 2);
const double clhs109 = pow(N[1], 2);
const double clhs110 = DN(1,1)*clhs54;
const double clhs111 = N[1]*clhs38 + N[1]*clhs40 + clhs108*mu + clhs109*clhs11 - clhs110*mu;
const double clhs112 = DN(1,1)*N[1];
const double clhs113 = DN(1,0)*DN(2,0);
const double clhs114 = DN(1,1)*DN(2,1);
const double clhs115 = clhs114*mu;
const double clhs116 = N[2]*clhs35 + clhs113*mu;
const double clhs117 = clhs113*clhs6 + clhs115 + clhs116;
const double clhs118 = DN(2,0)*N[1];
const double clhs119 = DN(2,1)*N[1];
const double clhs120 = DN(2,1)*clhs54;
const double clhs121 = clhs118*clhs49 + clhs119*clhs51 - clhs120*mu;
const double clhs122 = DN(1,1)*clhs6;
const double clhs123 = -clhs100 - clhs83 + clhs96*rho + clhs98*rho;
const double clhs124 = clhs123*clhs21;
const double clhs125 = N[2]*clhs83 + 2*clhs115 + clhs116;
const double clhs126 = DN(1,1)*clhs21;
const double clhs127 = N[2]*clhs54;
const double clhs128 = clhs113 + clhs114;
const double clhs129 = N[2]*clhs8 + clhs59;
const double clhs130 = clhs129*rho;
const double clhs131 = N[2]*clhs23 + clhs61;
const double clhs132 = clhs131*rho;
const double clhs133 = clhs14*(DN(2,1) - clhs74);
const double clhs134 = DN(0,1)*clhs74;
const double clhs135 = N[2]*clhs16 + N[2]*clhs18 - clhs134*mu;
const double clhs136 = DN(2,0)*clhs6;
const double clhs137 = DN(0,0)*N[2];
const double clhs138 = DN(0,1)*N[2];
const double clhs139 = DN(1,1)*clhs74;
const double clhs140 = N[2]*clhs38 + N[2]*clhs40 - clhs139*mu;
const double clhs141 = DN(1,0)*N[2];
const double clhs142 = DN(1,1)*N[2];
const double clhs143 = pow(DN(2,1), 2);
const double clhs144 = clhs143*mu;
const double clhs145 = pow(DN(2,0), 2);
const double clhs146 = pow(N[2], 2);
const double clhs147 = DN(2,1)*clhs74;
const double clhs148 = N[2]*clhs60 + N[2]*clhs62 + clhs11*clhs146 + clhs145*mu - clhs147*mu;
const double clhs149 = DN(2,1)*N[2];
const double clhs150 = DN(2,1)*clhs6;
const double clhs151 = clhs129*rho + clhs131*rho - clhs133 - clhs87;
const double clhs152 = DN(2,1)*clhs21;
lhs(0,0)=clhs1 + clhs10*clhs22 + clhs2*clhs6 + clhs22*clhs25 - clhs22*clhs27 + clhs30;
lhs(0,1)=clhs31*clhs32;
lhs(0,2)=DN(0,0)*N[0] + clhs10*clhs34 + clhs25*clhs34 - clhs27*clhs34 + clhs33;
lhs(0,3)=clhs10*clhs42 + clhs25*clhs42 - clhs27*clhs42 + clhs47 + clhs53;
lhs(0,4)=clhs32*clhs55;
lhs(0,5)=clhs10*clhs56 + clhs25*clhs56 - clhs27*clhs56 + clhs48 + clhs50;
lhs(0,6)=clhs10*clhs64 + clhs25*clhs64 - clhs27*clhs64 + clhs69 + clhs73;
lhs(0,7)=clhs32*clhs75;
lhs(0,8)=clhs10*clhs76 + clhs25*clhs76 - clhs27*clhs76 + clhs70 + clhs71;
lhs(1,0)=DN(0,1)*clhs32;
lhs(1,1)=2*clhs1 + clhs10*clhs81 + clhs25*clhs81 - clhs27*clhs81 + clhs28*clhs78 + clhs30 + clhs31*clhs77 - clhs79*clhs81;
lhs(1,2)=DN(0,1)*clhs82;
lhs(1,3)=DN(1,0)*clhs77;
lhs(1,4)=clhs10*clhs85 + clhs25*clhs85 - clhs27*clhs85 + clhs53 + clhs55*clhs77 - clhs79*clhs85 + clhs86;
lhs(1,5)=DN(1,1)*clhs82;
lhs(1,6)=DN(2,0)*clhs77;
lhs(1,7)=clhs10*clhs89 + clhs25*clhs89 - clhs27*clhs89 + clhs73 + clhs75*clhs77 - clhs79*clhs89 + clhs90;
lhs(1,8)=DN(2,1)*clhs82;
lhs(2,0)=DN(0,0)*(N[0] + clhs20*(1.0*clhs12 - 1.0*clhs15 + 1.0*clhs16 + 1.0*clhs18));
lhs(2,1)=clhs13*clhs28 - clhs26*clhs81 + clhs33 + clhs80*clhs91;
lhs(2,2)=clhs21*(clhs0 + clhs2 - clhs29);
lhs(2,3)=clhs34*clhs41 + clhs48;
lhs(2,4)=-clhs26*clhs85 + clhs50 + clhs84*clhs91 + clhs92;
lhs(2,5)=clhs21*(-clhs52 + clhs93);
lhs(2,6)=clhs34*clhs63 + clhs70;
lhs(2,7)=-clhs26*clhs89 + clhs71 + clhs88*clhs91 + clhs94;
lhs(2,8)=clhs21*(-clhs72 + clhs95);
lhs(3,0)=-clhs100*clhs22 + clhs102 + clhs22*clhs97 + clhs22*clhs99 + clhs47;
lhs(3,1)=clhs103*clhs31;
lhs(3,2)=-clhs100*clhs34 + clhs104 + clhs105 + clhs34*clhs97 + clhs34*clhs99;
lhs(3,3)=-clhs100*clhs42 + clhs107 + clhs108*clhs6 + clhs111 + clhs42*clhs97 + clhs42*clhs99;
lhs(3,4)=clhs103*clhs55;
lhs(3,5)=DN(1,0)*N[1] - clhs100*clhs56 + clhs112 + clhs56*clhs97 + clhs56*clhs99;
lhs(3,6)=-clhs100*clhs64 + clhs117 + clhs121 + clhs64*clhs97 + clhs64*clhs99;
lhs(3,7)=clhs103*clhs75;
lhs(3,8)=-clhs100*clhs76 + clhs118 + clhs119 + clhs76*clhs97 + clhs76*clhs99;
lhs(4,0)=DN(1,1)*clhs32;
lhs(4,1)=-clhs100*clhs81 + clhs102 + clhs122*clhs31 - clhs81*clhs83 + clhs81*clhs97 + clhs81*clhs99 + clhs86;
lhs(4,2)=clhs123*clhs91;
lhs(4,3)=DN(1,1)*clhs103;
lhs(4,4)=-clhs100*clhs85 + 2*clhs107 + clhs109*clhs78 + clhs111 + clhs122*clhs55 - clhs83*clhs85 + clhs85*clhs97 + clhs85*clhs99;
lhs(4,5)=DN(1,1)*clhs124;
lhs(4,6)=DN(2,0)*clhs122;
lhs(4,7)=-clhs100*clhs89 + clhs121 + clhs122*clhs75 + clhs125 - clhs83*clhs89 + clhs89*clhs97 + clhs89*clhs99;
lhs(4,8)=DN(2,1)*clhs124;
lhs(5,0)=clhs104 + clhs19*clhs56;
lhs(5,1)=clhs105 + clhs126*clhs80 - clhs54*clhs81 + clhs92;
lhs(5,2)=clhs21*(-clhs101 + clhs93);
lhs(5,3)=DN(1,0)*(N[1] + clhs20*(1.0*clhs35 - 1.0*clhs36 + 1.0*clhs38 + 1.0*clhs40));
lhs(5,4)=clhs109*clhs13 + clhs112 + clhs126*clhs84 - clhs54*clhs85;
lhs(5,5)=clhs21*(clhs106 + clhs108 - clhs110);
lhs(5,6)=clhs118 + clhs56*clhs63;
lhs(5,7)=clhs119 + clhs126*clhs88 + clhs127 - clhs54*clhs89;
lhs(5,8)=clhs21*(-clhs120 + clhs128);
lhs(6,0)=clhs130*clhs22 + clhs132*clhs22 - clhs133*clhs22 + clhs135 + clhs69;
lhs(6,1)=clhs136*clhs31;
lhs(6,2)=clhs130*clhs34 + clhs132*clhs34 - clhs133*clhs34 + clhs137 + clhs138;
lhs(6,3)=clhs117 + clhs130*clhs42 + clhs132*clhs42 - clhs133*clhs42 + clhs140;
lhs(6,4)=clhs136*clhs55;
lhs(6,5)=clhs130*clhs56 + clhs132*clhs56 - clhs133*clhs56 + clhs141 + clhs142;
lhs(6,6)=clhs130*clhs64 + clhs132*clhs64 - clhs133*clhs64 + clhs144 + clhs145*clhs6 + clhs148;
lhs(6,7)=clhs136*clhs75;
lhs(6,8)=DN(2,0)*N[2] + clhs130*clhs76 + clhs132*clhs76 - clhs133*clhs76 + clhs149;
lhs(7,0)=DN(2,1)*clhs32;
lhs(7,1)=clhs130*clhs81 + clhs132*clhs81 - clhs133*clhs81 + clhs135 + clhs150*clhs31 - clhs81*clhs87 + clhs90;
lhs(7,2)=clhs151*clhs91;
lhs(7,3)=DN(2,1)*clhs103;
lhs(7,4)=clhs125 + clhs130*clhs85 + clhs132*clhs85 - clhs133*clhs85 + clhs140 + clhs150*clhs55 - clhs85*clhs87;
lhs(7,5)=clhs126*clhs151;
lhs(7,6)=DN(2,1)*clhs136;
lhs(7,7)=clhs130*clhs89 + clhs132*clhs89 - clhs133*clhs89 + 2*clhs144 + clhs146*clhs78 + clhs148 + clhs150*clhs75 - clhs87*clhs89;
lhs(7,8)=clhs151*clhs152;
lhs(8,0)=clhs137 + clhs19*clhs76;
lhs(8,1)=clhs138 + clhs152*clhs80 - clhs74*clhs81 + clhs94;
lhs(8,2)=clhs21*(-clhs134 + clhs95);
lhs(8,3)=clhs141 + clhs41*clhs76;
lhs(8,4)=clhs127 + clhs142 + clhs152*clhs84 - clhs74*clhs85;
lhs(8,5)=clhs21*(clhs128 - clhs139);
lhs(8,6)=DN(2,0)*(N[2] + clhs20*(1.0*clhs57 - 1.0*clhs58 + 1.0*clhs60 + 1.0*clhs62));
lhs(8,7)=clhs13*clhs146 + clhs149 + clhs152*clhs88 - clhs74*clhs89;
lhs(8,8)=clhs21*(clhs143 + clhs145 - clhs147);


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

    const double clhs0 = pow(DN(0,1), 2);
const double clhs1 = clhs0*mu;
const double clhs2 = pow(DN(0,0), 2);
const double clhs3 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double clhs4 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double clhs5 = rho*stab_c2*sqrt(pow(clhs3, 2) + pow(clhs4, 2));
const double clhs6 = clhs5*h/stab_c1 + mu;
const double clhs7 = DN(0,0)*clhs3;
const double clhs8 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double clhs9 = N[0]*clhs8 + clhs7;
const double clhs10 = clhs9*rho;
const double clhs11 = bdf0*rho;
const double clhs12 = N[0]*clhs11;
const double clhs13 = 1.0/y;
const double clhs14 = clhs13*mu;
const double clhs15 = DN(0,1)*clhs14;
const double clhs16 = clhs7*rho;
const double clhs17 = DN(0,1)*clhs4;
const double clhs18 = clhs17*rho;
const double clhs19 = clhs12 - clhs15 + clhs16 + clhs18;
const double clhs20 = 1.0/(clhs5/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs21 = 1.0*clhs20;
const double clhs22 = clhs19*clhs21;
const double clhs23 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double clhs24 = N[0]*clhs23 + clhs17;
const double clhs25 = clhs24*rho;
const double clhs26 = N[0]*clhs13;
const double clhs27 = clhs14*(DN(0,1) - clhs26);
const double clhs28 = pow(N[0], 2);
const double clhs29 = DN(0,1)*clhs26;
const double clhs30 = N[0]*clhs16 + N[0]*clhs18 + clhs11*clhs28 + clhs2*mu - clhs29*mu;
const double clhs31 = DN(0,1) + clhs26;
const double clhs32 = DN(0,0)*clhs6;
const double clhs33 = DN(0,1)*N[0];
const double clhs34 = DN(0,0)*clhs21;
const double clhs35 = N[1]*clhs11;
const double clhs36 = DN(1,1)*clhs14;
const double clhs37 = DN(1,0)*clhs3;
const double clhs38 = clhs37*rho;
const double clhs39 = DN(1,1)*clhs4;
const double clhs40 = clhs39*rho;
const double clhs41 = clhs35 - clhs36 + clhs38 + clhs40;
const double clhs42 = clhs21*clhs41;
const double clhs43 = DN(0,0)*DN(1,0);
const double clhs44 = DN(0,1)*DN(1,1);
const double clhs45 = clhs44*mu;
const double clhs46 = N[1]*clhs12 + clhs43*mu;
const double clhs47 = clhs43*clhs6 + clhs45 + clhs46;
const double clhs48 = DN(1,0)*N[0];
const double clhs49 = clhs3*rho;
const double clhs50 = DN(1,1)*N[0];
const double clhs51 = clhs4*rho;
const double clhs52 = DN(1,1)*clhs26;
const double clhs53 = clhs48*clhs49 + clhs50*clhs51 - clhs52*mu;
const double clhs54 = N[1]*clhs13;
const double clhs55 = DN(1,1) + clhs54;
const double clhs56 = DN(1,0)*clhs21;
const double clhs57 = N[2]*clhs11;
const double clhs58 = DN(2,1)*clhs14;
const double clhs59 = DN(2,0)*clhs3;
const double clhs60 = clhs59*rho;
const double clhs61 = DN(2,1)*clhs4;
const double clhs62 = clhs61*rho;
const double clhs63 = clhs57 - clhs58 + clhs60 + clhs62;
const double clhs64 = clhs21*clhs63;
const double clhs65 = DN(0,0)*DN(2,0);
const double clhs66 = DN(0,1)*DN(2,1);
const double clhs67 = clhs66*mu;
const double clhs68 = N[2]*clhs12 + clhs65*mu;
const double clhs69 = clhs6*clhs65 + clhs67 + clhs68;
const double clhs70 = DN(2,0)*N[0];
const double clhs71 = DN(2,1)*N[0];
const double clhs72 = DN(2,1)*clhs26;
const double clhs73 = clhs49*clhs70 + clhs51*clhs71 - clhs72*mu;
const double clhs74 = N[2]*clhs13;
const double clhs75 = DN(2,1) + clhs74;
const double clhs76 = DN(2,0)*clhs21;
const double clhs77 = N[3]*clhs11;
const double clhs78 = DN(3,1)*clhs14;
const double clhs79 = DN(3,0)*clhs3;
const double clhs80 = clhs79*rho;
const double clhs81 = DN(3,1)*clhs4;
const double clhs82 = clhs81*rho;
const double clhs83 = clhs77 - clhs78 + clhs80 + clhs82;
const double clhs84 = clhs21*clhs83;
const double clhs85 = DN(0,0)*DN(3,0);
const double clhs86 = DN(0,1)*DN(3,1);
const double clhs87 = clhs86*mu;
const double clhs88 = N[3]*clhs12 + clhs85*mu;
const double clhs89 = clhs6*clhs85 + clhs87 + clhs88;
const double clhs90 = DN(3,0)*N[0];
const double clhs91 = DN(3,1)*N[0];
const double clhs92 = DN(3,1)*clhs26;
const double clhs93 = clhs49*clhs90 + clhs51*clhs91 - clhs92*mu;
const double clhs94 = N[3]*clhs13;
const double clhs95 = DN(3,1) + clhs94;
const double clhs96 = DN(3,0)*clhs21;
const double clhs97 = DN(0,1)*clhs6;
const double clhs98 = mu/pow(y, 2);
const double clhs99 = N[0]*clhs98;
const double clhs100 = clhs19 + clhs99;
const double clhs101 = clhs100*clhs21;
const double clhs102 = clhs21*(clhs24*rho - clhs27 + clhs9*rho - clhs99);
const double clhs103 = N[1]*clhs98;
const double clhs104 = clhs103 + clhs41;
const double clhs105 = clhs104*clhs21;
const double clhs106 = N[1]*clhs99 + 2*clhs45 + clhs46;
const double clhs107 = N[2]*clhs98;
const double clhs108 = clhs107 + clhs63;
const double clhs109 = clhs108*clhs21;
const double clhs110 = N[2]*clhs99 + 2*clhs67 + clhs68;
const double clhs111 = N[3]*clhs98;
const double clhs112 = clhs111 + clhs83;
const double clhs113 = clhs112*clhs21;
const double clhs114 = N[3]*clhs99 + 2*clhs87 + clhs88;
const double clhs115 = DN(0,1)*clhs21;
const double clhs116 = N[1]*clhs26;
const double clhs117 = clhs43 + clhs44;
const double clhs118 = N[2]*clhs26;
const double clhs119 = clhs65 + clhs66;
const double clhs120 = N[3]*clhs26;
const double clhs121 = clhs85 + clhs86;
const double clhs122 = N[1]*clhs8 + clhs37;
const double clhs123 = clhs122*rho;
const double clhs124 = N[1]*clhs23 + clhs39;
const double clhs125 = clhs124*rho;
const double clhs126 = clhs14*(DN(1,1) - clhs54);
const double clhs127 = DN(0,1)*clhs54;
const double clhs128 = N[1]*clhs16 + N[1]*clhs18 - clhs127*mu;
const double clhs129 = DN(1,0)*clhs6;
const double clhs130 = DN(0,0)*N[1];
const double clhs131 = DN(0,1)*N[1];
const double clhs132 = pow(DN(1,1), 2);
const double clhs133 = clhs132*mu;
const double clhs134 = pow(DN(1,0), 2);
const double clhs135 = pow(N[1], 2);
const double clhs136 = DN(1,1)*clhs54;
const double clhs137 = N[1]*clhs38 + N[1]*clhs40 + clhs11*clhs135 + clhs134*mu - clhs136*mu;
const double clhs138 = DN(1,1)*N[1];
const double clhs139 = DN(1,0)*DN(2,0);
const double clhs140 = DN(1,1)*DN(2,1);
const double clhs141 = clhs140*mu;
const double clhs142 = N[2]*clhs35 + clhs139*mu;
const double clhs143 = clhs139*clhs6 + clhs141 + clhs142;
const double clhs144 = DN(2,0)*N[1];
const double clhs145 = DN(2,1)*N[1];
const double clhs146 = DN(2,1)*clhs54;
const double clhs147 = clhs144*clhs49 + clhs145*clhs51 - clhs146*mu;
const double clhs148 = DN(1,0)*DN(3,0);
const double clhs149 = DN(1,1)*DN(3,1);
const double clhs150 = clhs149*mu;
const double clhs151 = N[3]*clhs35 + clhs148*mu;
const double clhs152 = clhs148*clhs6 + clhs150 + clhs151;
const double clhs153 = DN(3,0)*N[1];
const double clhs154 = DN(3,1)*N[1];
const double clhs155 = DN(3,1)*clhs54;
const double clhs156 = clhs153*clhs49 + clhs154*clhs51 - clhs155*mu;
const double clhs157 = DN(1,1)*clhs6;
const double clhs158 = -clhs103 + clhs122*rho + clhs124*rho - clhs126;
const double clhs159 = clhs158*clhs21;
const double clhs160 = N[2]*clhs103 + 2*clhs141 + clhs142;
const double clhs161 = N[3]*clhs103 + 2*clhs150 + clhs151;
const double clhs162 = DN(1,1)*clhs21;
const double clhs163 = N[2]*clhs54;
const double clhs164 = clhs139 + clhs140;
const double clhs165 = N[3]*clhs54;
const double clhs166 = clhs148 + clhs149;
const double clhs167 = N[2]*clhs8 + clhs59;
const double clhs168 = clhs167*rho;
const double clhs169 = N[2]*clhs23 + clhs61;
const double clhs170 = clhs169*rho;
const double clhs171 = clhs14*(DN(2,1) - clhs74);
const double clhs172 = DN(0,1)*clhs74;
const double clhs173 = N[2]*clhs16 + N[2]*clhs18 - clhs172*mu;
const double clhs174 = DN(2,0)*clhs6;
const double clhs175 = DN(0,0)*N[2];
const double clhs176 = DN(0,1)*N[2];
const double clhs177 = DN(1,1)*clhs74;
const double clhs178 = N[2]*clhs38 + N[2]*clhs40 - clhs177*mu;
const double clhs179 = DN(1,0)*N[2];
const double clhs180 = DN(1,1)*N[2];
const double clhs181 = pow(DN(2,1), 2);
const double clhs182 = clhs181*mu;
const double clhs183 = pow(DN(2,0), 2);
const double clhs184 = pow(N[2], 2);
const double clhs185 = DN(2,1)*clhs74;
const double clhs186 = N[2]*clhs60 + N[2]*clhs62 + clhs11*clhs184 + clhs183*mu - clhs185*mu;
const double clhs187 = DN(2,1)*N[2];
const double clhs188 = DN(2,0)*DN(3,0);
const double clhs189 = DN(2,1)*DN(3,1);
const double clhs190 = clhs189*mu;
const double clhs191 = N[3]*clhs57 + clhs188*mu;
const double clhs192 = clhs188*clhs6 + clhs190 + clhs191;
const double clhs193 = DN(3,0)*N[2];
const double clhs194 = DN(3,1)*N[2];
const double clhs195 = DN(3,1)*clhs74;
const double clhs196 = clhs193*clhs49 + clhs194*clhs51 - clhs195*mu;
const double clhs197 = DN(2,1)*clhs6;
const double clhs198 = -clhs107 + clhs167*rho + clhs169*rho - clhs171;
const double clhs199 = clhs198*clhs21;
const double clhs200 = N[3]*clhs107 + 2*clhs190 + clhs191;
const double clhs201 = DN(2,1)*clhs21;
const double clhs202 = N[3]*clhs74;
const double clhs203 = clhs188 + clhs189;
const double clhs204 = N[3]*clhs8 + clhs79;
const double clhs205 = clhs204*rho;
const double clhs206 = N[3]*clhs23 + clhs81;
const double clhs207 = clhs206*rho;
const double clhs208 = clhs14*(DN(3,1) - clhs94);
const double clhs209 = DN(0,1)*clhs94;
const double clhs210 = N[3]*clhs16 + N[3]*clhs18 - clhs209*mu;
const double clhs211 = DN(3,0)*clhs6;
const double clhs212 = DN(0,0)*N[3];
const double clhs213 = DN(0,1)*N[3];
const double clhs214 = DN(1,1)*clhs94;
const double clhs215 = N[3]*clhs38 + N[3]*clhs40 - clhs214*mu;
const double clhs216 = DN(1,0)*N[3];
const double clhs217 = DN(1,1)*N[3];
const double clhs218 = DN(2,1)*clhs94;
const double clhs219 = N[3]*clhs60 + N[3]*clhs62 - clhs218*mu;
const double clhs220 = DN(2,0)*N[3];
const double clhs221 = DN(2,1)*N[3];
const double clhs222 = pow(DN(3,1), 2);
const double clhs223 = clhs222*mu;
const double clhs224 = pow(DN(3,0), 2);
const double clhs225 = pow(N[3], 2);
const double clhs226 = DN(3,1)*clhs94;
const double clhs227 = N[3]*clhs80 + N[3]*clhs82 + clhs11*clhs225 + clhs224*mu - clhs226*mu;
const double clhs228 = DN(3,1)*N[3];
const double clhs229 = DN(3,1)*clhs6;
const double clhs230 = -clhs111 + clhs204*rho + clhs206*rho - clhs208;
const double clhs231 = DN(3,1)*clhs21;
lhs(0,0)=clhs1 + clhs10*clhs22 + clhs2*clhs6 + clhs22*clhs25 - clhs22*clhs27 + clhs30;
lhs(0,1)=clhs31*clhs32;
lhs(0,2)=DN(0,0)*N[0] + clhs10*clhs34 + clhs25*clhs34 - clhs27*clhs34 + clhs33;
lhs(0,3)=clhs10*clhs42 + clhs25*clhs42 - clhs27*clhs42 + clhs47 + clhs53;
lhs(0,4)=clhs32*clhs55;
lhs(0,5)=clhs10*clhs56 + clhs25*clhs56 - clhs27*clhs56 + clhs48 + clhs50;
lhs(0,6)=clhs10*clhs64 + clhs25*clhs64 - clhs27*clhs64 + clhs69 + clhs73;
lhs(0,7)=clhs32*clhs75;
lhs(0,8)=clhs10*clhs76 + clhs25*clhs76 - clhs27*clhs76 + clhs70 + clhs71;
lhs(0,9)=clhs10*clhs84 + clhs25*clhs84 - clhs27*clhs84 + clhs89 + clhs93;
lhs(0,10)=clhs32*clhs95;
lhs(0,11)=clhs10*clhs96 + clhs25*clhs96 - clhs27*clhs96 + clhs90 + clhs91;
lhs(1,0)=DN(0,1)*clhs32;
lhs(1,1)=2*clhs1 + clhs10*clhs101 + clhs101*clhs25 - clhs101*clhs27 - clhs101*clhs99 + clhs28*clhs98 + clhs30 + clhs31*clhs97;
lhs(1,2)=DN(0,1)*clhs102;
lhs(1,3)=DN(1,0)*clhs97;
lhs(1,4)=clhs10*clhs105 + clhs105*clhs25 - clhs105*clhs27 - clhs105*clhs99 + clhs106 + clhs53 + clhs55*clhs97;
lhs(1,5)=DN(1,1)*clhs102;
lhs(1,6)=DN(2,0)*clhs97;
lhs(1,7)=clhs10*clhs109 + clhs109*clhs25 - clhs109*clhs27 - clhs109*clhs99 + clhs110 + clhs73 + clhs75*clhs97;
lhs(1,8)=DN(2,1)*clhs102;
lhs(1,9)=DN(3,0)*clhs97;
lhs(1,10)=clhs10*clhs113 + clhs113*clhs25 - clhs113*clhs27 - clhs113*clhs99 + clhs114 + clhs93 + clhs95*clhs97;
lhs(1,11)=DN(3,1)*clhs102;
lhs(2,0)=DN(0,0)*(N[0] + clhs20*(1.0*clhs12 - 1.0*clhs15 + 1.0*clhs16 + 1.0*clhs18));
lhs(2,1)=clhs100*clhs115 - clhs101*clhs26 + clhs13*clhs28 + clhs33;
lhs(2,2)=clhs21*(clhs0 + clhs2 - clhs29);
lhs(2,3)=clhs34*clhs41 + clhs48;
lhs(2,4)=clhs104*clhs115 - clhs105*clhs26 + clhs116 + clhs50;
lhs(2,5)=clhs21*(clhs117 - clhs52);
lhs(2,6)=clhs34*clhs63 + clhs70;
lhs(2,7)=clhs108*clhs115 - clhs109*clhs26 + clhs118 + clhs71;
lhs(2,8)=clhs21*(clhs119 - clhs72);
lhs(2,9)=clhs34*clhs83 + clhs90;
lhs(2,10)=clhs112*clhs115 - clhs113*clhs26 + clhs120 + clhs91;
lhs(2,11)=clhs21*(clhs121 - clhs92);
lhs(3,0)=clhs123*clhs22 + clhs125*clhs22 - clhs126*clhs22 + clhs128 + clhs47;
lhs(3,1)=clhs129*clhs31;
lhs(3,2)=clhs123*clhs34 + clhs125*clhs34 - clhs126*clhs34 + clhs130 + clhs131;
lhs(3,3)=clhs123*clhs42 + clhs125*clhs42 - clhs126*clhs42 + clhs133 + clhs134*clhs6 + clhs137;
lhs(3,4)=clhs129*clhs55;
lhs(3,5)=DN(1,0)*N[1] + clhs123*clhs56 + clhs125*clhs56 - clhs126*clhs56 + clhs138;
lhs(3,6)=clhs123*clhs64 + clhs125*clhs64 - clhs126*clhs64 + clhs143 + clhs147;
lhs(3,7)=clhs129*clhs75;
lhs(3,8)=clhs123*clhs76 + clhs125*clhs76 - clhs126*clhs76 + clhs144 + clhs145;
lhs(3,9)=clhs123*clhs84 + clhs125*clhs84 - clhs126*clhs84 + clhs152 + clhs156;
lhs(3,10)=clhs129*clhs95;
lhs(3,11)=clhs123*clhs96 + clhs125*clhs96 - clhs126*clhs96 + clhs153 + clhs154;
lhs(4,0)=DN(1,1)*clhs32;
lhs(4,1)=-clhs101*clhs103 + clhs101*clhs123 + clhs101*clhs125 - clhs101*clhs126 + clhs106 + clhs128 + clhs157*clhs31;
lhs(4,2)=clhs115*clhs158;
lhs(4,3)=DN(1,1)*clhs129;
lhs(4,4)=-clhs103*clhs105 + clhs105*clhs123 + clhs105*clhs125 - clhs105*clhs126 + 2*clhs133 + clhs135*clhs98 + clhs137 + clhs157*clhs55;
lhs(4,5)=DN(1,1)*clhs159;
lhs(4,6)=DN(2,0)*clhs157;
lhs(4,7)=-clhs103*clhs109 + clhs109*clhs123 + clhs109*clhs125 - clhs109*clhs126 + clhs147 + clhs157*clhs75 + clhs160;
lhs(4,8)=DN(2,1)*clhs159;
lhs(4,9)=DN(3,0)*clhs157;
lhs(4,10)=-clhs103*clhs113 + clhs113*clhs123 + clhs113*clhs125 - clhs113*clhs126 + clhs156 + clhs157*clhs95 + clhs161;
lhs(4,11)=DN(3,1)*clhs159;
lhs(5,0)=clhs130 + clhs19*clhs56;
lhs(5,1)=clhs100*clhs162 - clhs101*clhs54 + clhs116 + clhs131;
lhs(5,2)=clhs21*(clhs117 - clhs127);
lhs(5,3)=DN(1,0)*(N[1] + clhs20*(1.0*clhs35 - 1.0*clhs36 + 1.0*clhs38 + 1.0*clhs40));
lhs(5,4)=clhs104*clhs162 - clhs105*clhs54 + clhs13*clhs135 + clhs138;
lhs(5,5)=clhs21*(clhs132 + clhs134 - clhs136);
lhs(5,6)=clhs144 + clhs56*clhs63;
lhs(5,7)=clhs108*clhs162 - clhs109*clhs54 + clhs145 + clhs163;
lhs(5,8)=clhs21*(-clhs146 + clhs164);
lhs(5,9)=clhs153 + clhs56*clhs83;
lhs(5,10)=clhs112*clhs162 - clhs113*clhs54 + clhs154 + clhs165;
lhs(5,11)=clhs21*(-clhs155 + clhs166);
lhs(6,0)=clhs168*clhs22 + clhs170*clhs22 - clhs171*clhs22 + clhs173 + clhs69;
lhs(6,1)=clhs174*clhs31;
lhs(6,2)=clhs168*clhs34 + clhs170*clhs34 - clhs171*clhs34 + clhs175 + clhs176;
lhs(6,3)=clhs143 + clhs168*clhs42 + clhs170*clhs42 - clhs171*clhs42 + clhs178;
lhs(6,4)=clhs174*clhs55;
lhs(6,5)=clhs168*clhs56 + clhs170*clhs56 - clhs171*clhs56 + clhs179 + clhs180;
lhs(6,6)=clhs168*clhs64 + clhs170*clhs64 - clhs171*clhs64 + clhs182 + clhs183*clhs6 + clhs186;
lhs(6,7)=clhs174*clhs75;
lhs(6,8)=DN(2,0)*N[2] + clhs168*clhs76 + clhs170*clhs76 - clhs171*clhs76 + clhs187;
lhs(6,9)=clhs168*clhs84 + clhs170*clhs84 - clhs171*clhs84 + clhs192 + clhs196;
lhs(6,10)=clhs174*clhs95;
lhs(6,11)=clhs168*clhs96 + clhs170*clhs96 - clhs171*clhs96 + clhs193 + clhs194;
lhs(7,0)=DN(2,1)*clhs32;
lhs(7,1)=-clhs101*clhs107 + clhs101*clhs168 + clhs101*clhs170 - clhs101*clhs171 + clhs110 + clhs173 + clhs197*clhs31;
lhs(7,2)=clhs115*clhs198;
lhs(7,3)=DN(2,1)*clhs129;
lhs(7,4)=-clhs105*clhs107 + clhs105*clhs168 + clhs105*clhs170 - clhs105*clhs171 + clhs160 + clhs178 + clhs197*clhs55;
lhs(7,5)=clhs162*clhs198;
lhs(7,6)=DN(2,1)*clhs174;
lhs(7,7)=-clhs107*clhs109 + clhs109*clhs168 + clhs109*clhs170 - clhs109*clhs171 + 2*clhs182 + clhs184*clhs98 + clhs186 + clhs197*clhs75;
lhs(7,8)=DN(2,1)*clhs199;
lhs(7,9)=DN(3,0)*clhs197;
lhs(7,10)=-clhs107*clhs113 + clhs113*clhs168 + clhs113*clhs170 - clhs113*clhs171 + clhs196 + clhs197*clhs95 + clhs200;
lhs(7,11)=DN(3,1)*clhs199;
lhs(8,0)=clhs175 + clhs19*clhs76;
lhs(8,1)=clhs100*clhs201 - clhs101*clhs74 + clhs118 + clhs176;
lhs(8,2)=clhs21*(clhs119 - clhs172);
lhs(8,3)=clhs179 + clhs41*clhs76;
lhs(8,4)=clhs104*clhs201 - clhs105*clhs74 + clhs163 + clhs180;
lhs(8,5)=clhs21*(clhs164 - clhs177);
lhs(8,6)=DN(2,0)*(N[2] + clhs20*(1.0*clhs57 - 1.0*clhs58 + 1.0*clhs60 + 1.0*clhs62));
lhs(8,7)=clhs108*clhs201 - clhs109*clhs74 + clhs13*clhs184 + clhs187;
lhs(8,8)=clhs21*(clhs181 + clhs183 - clhs185);
lhs(8,9)=clhs193 + clhs76*clhs83;
lhs(8,10)=clhs112*clhs201 - clhs113*clhs74 + clhs194 + clhs202;
lhs(8,11)=clhs21*(-clhs195 + clhs203);
lhs(9,0)=clhs205*clhs22 + clhs207*clhs22 - clhs208*clhs22 + clhs210 + clhs89;
lhs(9,1)=clhs211*clhs31;
lhs(9,2)=clhs205*clhs34 + clhs207*clhs34 - clhs208*clhs34 + clhs212 + clhs213;
lhs(9,3)=clhs152 + clhs205*clhs42 + clhs207*clhs42 - clhs208*clhs42 + clhs215;
lhs(9,4)=clhs211*clhs55;
lhs(9,5)=clhs205*clhs56 + clhs207*clhs56 - clhs208*clhs56 + clhs216 + clhs217;
lhs(9,6)=clhs192 + clhs205*clhs64 + clhs207*clhs64 - clhs208*clhs64 + clhs219;
lhs(9,7)=clhs211*clhs75;
lhs(9,8)=clhs205*clhs76 + clhs207*clhs76 - clhs208*clhs76 + clhs220 + clhs221;
lhs(9,9)=clhs205*clhs84 + clhs207*clhs84 - clhs208*clhs84 + clhs223 + clhs224*clhs6 + clhs227;
lhs(9,10)=clhs211*clhs95;
lhs(9,11)=DN(3,0)*N[3] + clhs205*clhs96 + clhs207*clhs96 - clhs208*clhs96 + clhs228;
lhs(10,0)=DN(3,1)*clhs32;
lhs(10,1)=-clhs101*clhs111 + clhs101*clhs205 + clhs101*clhs207 - clhs101*clhs208 + clhs114 + clhs210 + clhs229*clhs31;
lhs(10,2)=clhs115*clhs230;
lhs(10,3)=DN(3,1)*clhs129;
lhs(10,4)=-clhs105*clhs111 + clhs105*clhs205 + clhs105*clhs207 - clhs105*clhs208 + clhs161 + clhs215 + clhs229*clhs55;
lhs(10,5)=clhs162*clhs230;
lhs(10,6)=DN(3,1)*clhs174;
lhs(10,7)=-clhs109*clhs111 + clhs109*clhs205 + clhs109*clhs207 - clhs109*clhs208 + clhs200 + clhs219 + clhs229*clhs75;
lhs(10,8)=clhs201*clhs230;
lhs(10,9)=DN(3,1)*clhs211;
lhs(10,10)=-clhs111*clhs113 + clhs113*clhs205 + clhs113*clhs207 - clhs113*clhs208 + 2*clhs223 + clhs225*clhs98 + clhs227 + clhs229*clhs95;
lhs(10,11)=clhs230*clhs231;
lhs(11,0)=clhs19*clhs96 + clhs212;
lhs(11,1)=clhs100*clhs231 - clhs101*clhs94 + clhs120 + clhs213;
lhs(11,2)=clhs21*(clhs121 - clhs209);
lhs(11,3)=clhs216 + clhs41*clhs96;
lhs(11,4)=clhs104*clhs231 - clhs105*clhs94 + clhs165 + clhs217;
lhs(11,5)=clhs21*(clhs166 - clhs214);
lhs(11,6)=clhs220 + clhs63*clhs96;
lhs(11,7)=clhs108*clhs231 - clhs109*clhs94 + clhs202 + clhs221;
lhs(11,8)=clhs21*(clhs203 - clhs218);
lhs(11,9)=DN(3,0)*(N[3] + clhs20*(1.0*clhs77 - 1.0*clhs78 + 1.0*clhs80 + 1.0*clhs82));
lhs(11,10)=clhs112*clhs231 - clhs113*clhs94 + clhs13*clhs225 + clhs228;
lhs(11,11)=clhs21*(clhs222 + clhs224 - clhs226);


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
const double crhs2 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs3 = DN(0,0)*mu;
const double crhs4 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double crhs5 = crhs4*mu;
const double crhs6 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs7 = 1.0/y;
const double crhs8 = N[0]*crhs2;
const double crhs9 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double crhs10 = crhs9*rho;
const double crhs11 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double crhs12 = crhs11*rho;
const double crhs13 = crhs12*crhs4;
const double crhs14 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
const double crhs15 = rho*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs9, 2));
const double crhs16 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs17 = crhs16*crhs7;
const double crhs18 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs19 = (crhs15*h/stab_c1 + mu)*(crhs17 + crhs18 + crhs2);
const double crhs20 = DN(0,1) - N[0]*crhs7;
const double crhs21 = 1.0/(crhs15/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = crhs0 + crhs10*crhs2 + crhs13 + crhs14 - crhs5*crhs7 - crhs6*rho;
const double crhs23 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double crhs24 = DN(0,0)*crhs9 + N[0]*crhs23;
const double crhs25 = 1.0*crhs21;
const double crhs26 = crhs22*crhs25;
const double crhs27 = crhs26*rho;
const double crhs28 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double crhs29 = DN(0,1)*crhs11 + N[0]*crhs28;
const double crhs30 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs31 = crhs18*mu;
const double crhs32 = 2*crhs31;
const double crhs33 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs34 = pow(y, -2);
const double crhs35 = crhs16*crhs34*mu;
const double crhs36 = crhs10*crhs30;
const double crhs37 = N[0]*crhs18;
const double crhs38 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
const double crhs39 = crhs1 + crhs12*crhs18 - crhs31*crhs7 - crhs33*rho + crhs35 + crhs36 + crhs38;
const double crhs40 = crhs25*crhs39;
const double crhs41 = crhs40*rho;
const double crhs42 = DN(1,0)*mu;
const double crhs43 = N[1]*crhs2;
const double crhs44 = DN(1,1) - N[1]*crhs7;
const double crhs45 = DN(1,0)*crhs9 + N[1]*crhs23;
const double crhs46 = DN(1,1)*crhs11 + N[1]*crhs28;
const double crhs47 = N[1]*crhs18;
const double crhs48 = DN(2,0)*mu;
const double crhs49 = N[2]*crhs2;
const double crhs50 = DN(2,1) - N[2]*crhs7;
const double crhs51 = DN(2,0)*crhs9 + N[2]*crhs23;
const double crhs52 = DN(2,1)*crhs11 + N[2]*crhs28;
const double crhs53 = N[2]*crhs18;
rhs[0]=-DN(0,0)*crhs19 - DN(0,1)*crhs5 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs13 - N[0]*crhs14 + N[0]*crhs4*crhs7*mu + N[0]*crhs6*rho - crhs10*crhs8 - crhs2*crhs3 + 1.0*crhs20*crhs21*crhs22*crhs7*mu - crhs24*crhs27 - crhs27*crhs29;
rhs[1]=-DN(0,1)*crhs19 - DN(0,1)*crhs32 + N[0]*crhs18*crhs7*mu + 1.0*N[0]*crhs21*crhs34*crhs39*mu + N[0]*crhs33*rho - N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs38 - crhs12*crhs37 + 1.0*crhs20*crhs21*crhs39*crhs7*mu - crhs24*crhs41 - crhs29*crhs41 - crhs3*crhs30;
rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs40 - N[0]*crhs17 + 1.0*N[0]*crhs21*crhs39*crhs7 - crhs37 - crhs8;
rhs[3]=-DN(1,0)*crhs19 - DN(1,1)*crhs5 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs13 - N[1]*crhs14 + N[1]*crhs4*crhs7*mu + N[1]*crhs6*rho - crhs10*crhs43 - crhs2*crhs42 + 1.0*crhs21*crhs22*crhs44*crhs7*mu - crhs27*crhs45 - crhs27*crhs46;
rhs[4]=-DN(1,1)*crhs19 - DN(1,1)*crhs32 + N[1]*crhs18*crhs7*mu + 1.0*N[1]*crhs21*crhs34*crhs39*mu + N[1]*crhs33*rho - N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs38 - crhs12*crhs47 + 1.0*crhs21*crhs39*crhs44*crhs7*mu - crhs30*crhs42 - crhs41*crhs45 - crhs41*crhs46;
rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs40 - N[1]*crhs17 + 1.0*N[1]*crhs21*crhs39*crhs7 - crhs43 - crhs47;
rhs[6]=-DN(2,0)*crhs19 - DN(2,1)*crhs5 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs13 - N[2]*crhs14 + N[2]*crhs4*crhs7*mu + N[2]*crhs6*rho - crhs10*crhs49 - crhs2*crhs48 + 1.0*crhs21*crhs22*crhs50*crhs7*mu - crhs27*crhs51 - crhs27*crhs52;
rhs[7]=-DN(2,1)*crhs19 - DN(2,1)*crhs32 + N[2]*crhs18*crhs7*mu + 1.0*N[2]*crhs21*crhs34*crhs39*mu + N[2]*crhs33*rho - N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs38 - crhs12*crhs53 + 1.0*crhs21*crhs39*crhs50*crhs7*mu - crhs30*crhs48 - crhs41*crhs51 - crhs41*crhs52;
rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs40 - N[2]*crhs17 + 1.0*N[2]*crhs21*crhs39*crhs7 - crhs49 - crhs53;


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
const double crhs2 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs3 = DN(0,0)*mu;
const double crhs4 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double crhs5 = crhs4*mu;
const double crhs6 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs7 = 1.0/y;
const double crhs8 = N[0]*crhs2;
const double crhs9 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double crhs10 = crhs9*rho;
const double crhs11 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double crhs12 = crhs11*rho;
const double crhs13 = crhs12*crhs4;
const double crhs14 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
const double crhs15 = rho*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs9, 2));
const double crhs16 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs17 = crhs16*crhs7;
const double crhs18 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs19 = (crhs15*h/stab_c1 + mu)*(crhs17 + crhs18 + crhs2);
const double crhs20 = DN(0,1) - N[0]*crhs7;
const double crhs21 = 1.0/(crhs15/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = crhs0 + crhs10*crhs2 + crhs13 + crhs14 - crhs5*crhs7 - crhs6*rho;
const double crhs23 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double crhs24 = DN(0,0)*crhs9 + N[0]*crhs23;
const double crhs25 = 1.0*crhs21;
const double crhs26 = crhs22*crhs25;
const double crhs27 = crhs26*rho;
const double crhs28 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double crhs29 = DN(0,1)*crhs11 + N[0]*crhs28;
const double crhs30 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs31 = crhs18*mu;
const double crhs32 = 2*crhs31;
const double crhs33 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs34 = pow(y, -2);
const double crhs35 = crhs16*crhs34*mu;
const double crhs36 = crhs10*crhs30;
const double crhs37 = N[0]*crhs18;
const double crhs38 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
const double crhs39 = crhs1 + crhs12*crhs18 - crhs31*crhs7 - crhs33*rho + crhs35 + crhs36 + crhs38;
const double crhs40 = crhs25*crhs39;
const double crhs41 = crhs40*rho;
const double crhs42 = DN(1,0)*mu;
const double crhs43 = N[1]*crhs2;
const double crhs44 = DN(1,1) - N[1]*crhs7;
const double crhs45 = DN(1,0)*crhs9 + N[1]*crhs23;
const double crhs46 = DN(1,1)*crhs11 + N[1]*crhs28;
const double crhs47 = N[1]*crhs18;
const double crhs48 = DN(2,0)*mu;
const double crhs49 = N[2]*crhs2;
const double crhs50 = DN(2,1) - N[2]*crhs7;
const double crhs51 = DN(2,0)*crhs9 + N[2]*crhs23;
const double crhs52 = DN(2,1)*crhs11 + N[2]*crhs28;
const double crhs53 = N[2]*crhs18;
const double crhs54 = DN(3,0)*mu;
const double crhs55 = N[3]*crhs2;
const double crhs56 = DN(3,1) - N[3]*crhs7;
const double crhs57 = DN(3,0)*crhs9 + N[3]*crhs23;
const double crhs58 = DN(3,1)*crhs11 + N[3]*crhs28;
const double crhs59 = N[3]*crhs18;
rhs[0]=-DN(0,0)*crhs19 - DN(0,1)*crhs5 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs13 - N[0]*crhs14 + N[0]*crhs4*crhs7*mu + N[0]*crhs6*rho - crhs10*crhs8 - crhs2*crhs3 + 1.0*crhs20*crhs21*crhs22*crhs7*mu - crhs24*crhs27 - crhs27*crhs29;
rhs[1]=-DN(0,1)*crhs19 - DN(0,1)*crhs32 + N[0]*crhs18*crhs7*mu + 1.0*N[0]*crhs21*crhs34*crhs39*mu + N[0]*crhs33*rho - N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs38 - crhs12*crhs37 + 1.0*crhs20*crhs21*crhs39*crhs7*mu - crhs24*crhs41 - crhs29*crhs41 - crhs3*crhs30;
rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs40 - N[0]*crhs17 + 1.0*N[0]*crhs21*crhs39*crhs7 - crhs37 - crhs8;
rhs[3]=-DN(1,0)*crhs19 - DN(1,1)*crhs5 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs13 - N[1]*crhs14 + N[1]*crhs4*crhs7*mu + N[1]*crhs6*rho - crhs10*crhs43 - crhs2*crhs42 + 1.0*crhs21*crhs22*crhs44*crhs7*mu - crhs27*crhs45 - crhs27*crhs46;
rhs[4]=-DN(1,1)*crhs19 - DN(1,1)*crhs32 + N[1]*crhs18*crhs7*mu + 1.0*N[1]*crhs21*crhs34*crhs39*mu + N[1]*crhs33*rho - N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs38 - crhs12*crhs47 + 1.0*crhs21*crhs39*crhs44*crhs7*mu - crhs30*crhs42 - crhs41*crhs45 - crhs41*crhs46;
rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs40 - N[1]*crhs17 + 1.0*N[1]*crhs21*crhs39*crhs7 - crhs43 - crhs47;
rhs[6]=-DN(2,0)*crhs19 - DN(2,1)*crhs5 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs13 - N[2]*crhs14 + N[2]*crhs4*crhs7*mu + N[2]*crhs6*rho - crhs10*crhs49 - crhs2*crhs48 + 1.0*crhs21*crhs22*crhs50*crhs7*mu - crhs27*crhs51 - crhs27*crhs52;
rhs[7]=-DN(2,1)*crhs19 - DN(2,1)*crhs32 + N[2]*crhs18*crhs7*mu + 1.0*N[2]*crhs21*crhs34*crhs39*mu + N[2]*crhs33*rho - N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs38 - crhs12*crhs53 + 1.0*crhs21*crhs39*crhs50*crhs7*mu - crhs30*crhs48 - crhs41*crhs51 - crhs41*crhs52;
rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs40 - N[2]*crhs17 + 1.0*N[2]*crhs21*crhs39*crhs7 - crhs49 - crhs53;
rhs[9]=-DN(3,0)*crhs19 - DN(3,1)*crhs5 - N[3]*crhs0 - N[3]*crhs1 - N[3]*crhs13 - N[3]*crhs14 + N[3]*crhs4*crhs7*mu + N[3]*crhs6*rho - crhs10*crhs55 - crhs2*crhs54 + 1.0*crhs21*crhs22*crhs56*crhs7*mu - crhs27*crhs57 - crhs27*crhs58;
rhs[10]=-DN(3,1)*crhs19 - DN(3,1)*crhs32 + N[3]*crhs18*crhs7*mu + 1.0*N[3]*crhs21*crhs34*crhs39*mu + N[3]*crhs33*rho - N[3]*crhs35 - N[3]*crhs36 - N[3]*crhs38 - crhs12*crhs59 + 1.0*crhs21*crhs39*crhs56*crhs7*mu - crhs30*crhs54 - crhs41*crhs57 - crhs41*crhs58;
rhs[11]=-DN(3,0)*crhs26 - DN(3,1)*crhs40 - N[3]*crhs17 + 1.0*N[3]*crhs21*crhs39*crhs7 - crhs55 - crhs59;


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