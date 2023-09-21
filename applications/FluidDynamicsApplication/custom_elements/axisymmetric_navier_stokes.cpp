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
const double clhs7 = bdf0*rho;
const double clhs8 = N[0]*clhs7;
const double clhs9 = 1.0/y;
const double clhs10 = clhs9*mu;
const double clhs11 = DN(0,1)*clhs10;
const double clhs12 = DN(0,0)*clhs3;
const double clhs13 = clhs12*rho;
const double clhs14 = DN(0,1)*clhs4;
const double clhs15 = clhs14*rho;
const double clhs16 = -clhs11 + clhs13 + clhs15 + clhs8;
const double clhs17 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double clhs18 = N[0]*clhs17 + clhs12;
const double clhs19 = 1.0/(clhs5/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs20 = 1.0*clhs19;
const double clhs21 = clhs20*rho;
const double clhs22 = clhs18*clhs21;
const double clhs23 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double clhs24 = N[0]*clhs23 + clhs14;
const double clhs25 = clhs21*clhs24;
const double clhs26 = N[0]*clhs9;
const double clhs27 = DN(0,1) - clhs26;
const double clhs28 = clhs10*clhs20;
const double clhs29 = clhs27*clhs28;
const double clhs30 = pow(N[0], 2);
const double clhs31 = DN(0,1)*clhs26;
const double clhs32 = N[0]*clhs13 + N[0]*clhs15 + clhs2*mu + clhs30*clhs7 - clhs31*mu;
const double clhs33 = DN(0,1) + clhs26;
const double clhs34 = DN(0,0)*clhs6;
const double clhs35 = N[0] - 1.0*clhs18*clhs19*rho - 1.0*clhs19*clhs24*rho;
const double clhs36 = N[1]*clhs7;
const double clhs37 = DN(1,1)*clhs10;
const double clhs38 = DN(1,0)*clhs3;
const double clhs39 = clhs38*rho;
const double clhs40 = DN(1,1)*clhs4;
const double clhs41 = clhs40*rho;
const double clhs42 = clhs36 - clhs37 + clhs39 + clhs41;
const double clhs43 = DN(0,0)*DN(1,0);
const double clhs44 = DN(0,1)*DN(1,1);
const double clhs45 = clhs44*mu;
const double clhs46 = N[1]*clhs8 + clhs43*mu;
const double clhs47 = clhs43*clhs6 + clhs45 + clhs46;
const double clhs48 = DN(1,0)*N[0];
const double clhs49 = clhs3*rho;
const double clhs50 = DN(1,1)*N[0];
const double clhs51 = clhs4*rho;
const double clhs52 = DN(1,1)*clhs26;
const double clhs53 = clhs48*clhs49 + clhs50*clhs51 - clhs52*mu;
const double clhs54 = N[1]*clhs9;
const double clhs55 = DN(1,1) + clhs54;
const double clhs56 = DN(0,0)*N[1];
const double clhs57 = DN(1,0)*clhs20;
const double clhs58 = clhs10*clhs27;
const double clhs59 = N[2]*clhs7;
const double clhs60 = DN(2,1)*clhs10;
const double clhs61 = DN(2,0)*clhs3;
const double clhs62 = clhs61*rho;
const double clhs63 = DN(2,1)*clhs4;
const double clhs64 = clhs63*rho;
const double clhs65 = clhs59 - clhs60 + clhs62 + clhs64;
const double clhs66 = DN(0,0)*DN(2,0);
const double clhs67 = DN(0,1)*DN(2,1);
const double clhs68 = clhs67*mu;
const double clhs69 = N[2]*clhs8 + clhs66*mu;
const double clhs70 = clhs6*clhs66 + clhs68 + clhs69;
const double clhs71 = DN(2,0)*N[0];
const double clhs72 = DN(2,1)*N[0];
const double clhs73 = DN(2,1)*clhs26;
const double clhs74 = clhs49*clhs71 + clhs51*clhs72 - clhs73*mu;
const double clhs75 = N[2]*clhs9;
const double clhs76 = DN(2,1) + clhs75;
const double clhs77 = DN(0,0)*N[2];
const double clhs78 = DN(2,0)*clhs20;
const double clhs79 = DN(0,1)*clhs6;
const double clhs80 = mu/pow(y, 2);
const double clhs81 = N[0]*clhs80;
const double clhs82 = clhs16 + clhs81;
const double clhs83 = clhs11*clhs20;
const double clhs84 = N[1]*clhs80;
const double clhs85 = clhs42 + clhs84;
const double clhs86 = N[1]*clhs81 + 2*clhs45 + clhs46;
const double clhs87 = DN(0,1)*N[1];
const double clhs88 = clhs28*clhs44;
const double clhs89 = N[2]*clhs80 + clhs65;
const double clhs90 = N[2]*clhs81 + 2*clhs68 + clhs69;
const double clhs91 = DN(0,1)*N[2];
const double clhs92 = clhs28*clhs67;
const double clhs93 = DN(0,1)*clhs20;
const double clhs94 = clhs20*clhs26;
const double clhs95 = DN(0,0)*clhs20;
const double clhs96 = N[1]*clhs26;
const double clhs97 = clhs43 + clhs44;
const double clhs98 = N[2]*clhs26;
const double clhs99 = clhs66 + clhs67;
const double clhs100 = N[1]*clhs17 + clhs38;
const double clhs101 = clhs100*clhs21;
const double clhs102 = N[1]*clhs23 + clhs40;
const double clhs103 = clhs102*clhs21;
const double clhs104 = DN(1,1) - clhs54;
const double clhs105 = clhs104*clhs28;
const double clhs106 = DN(0,1)*clhs54;
const double clhs107 = N[1]*clhs13 + N[1]*clhs15 - clhs106*mu;
const double clhs108 = DN(1,0)*clhs6;
const double clhs109 = clhs10*clhs104;
const double clhs110 = pow(DN(1,1), 2);
const double clhs111 = clhs110*mu;
const double clhs112 = pow(DN(1,0), 2);
const double clhs113 = pow(N[1], 2);
const double clhs114 = DN(1,1)*clhs54;
const double clhs115 = N[1]*clhs39 + N[1]*clhs41 + clhs112*mu + clhs113*clhs7 - clhs114*mu;
const double clhs116 = N[1] - 1.0*clhs100*clhs19*rho - 1.0*clhs102*clhs19*rho;
const double clhs117 = DN(1,0)*DN(2,0);
const double clhs118 = DN(1,1)*DN(2,1);
const double clhs119 = clhs118*mu;
const double clhs120 = N[2]*clhs36 + clhs117*mu;
const double clhs121 = clhs117*clhs6 + clhs119 + clhs120;
const double clhs122 = DN(2,0)*N[1];
const double clhs123 = DN(2,1)*N[1];
const double clhs124 = DN(2,1)*clhs54;
const double clhs125 = clhs122*clhs49 + clhs123*clhs51 - clhs124*mu;
const double clhs126 = DN(1,0)*N[2];
const double clhs127 = DN(1,1)*clhs6;
const double clhs128 = clhs20*clhs37;
const double clhs129 = N[2]*clhs84 + 2*clhs119 + clhs120;
const double clhs130 = DN(1,1)*N[2];
const double clhs131 = clhs118*clhs28;
const double clhs132 = DN(1,1)*clhs20;
const double clhs133 = clhs20*clhs54;
const double clhs134 = N[2]*clhs54;
const double clhs135 = clhs117 + clhs118;
const double clhs136 = N[2]*clhs17 + clhs61;
const double clhs137 = clhs136*clhs21;
const double clhs138 = N[2]*clhs23 + clhs63;
const double clhs139 = clhs138*clhs21;
const double clhs140 = DN(2,1) - clhs75;
const double clhs141 = clhs140*clhs28;
const double clhs142 = DN(0,1)*clhs75;
const double clhs143 = N[2]*clhs13 + N[2]*clhs15 - clhs142*mu;
const double clhs144 = DN(2,0)*clhs6;
const double clhs145 = clhs10*clhs140;
const double clhs146 = DN(1,1)*clhs75;
const double clhs147 = N[2]*clhs39 + N[2]*clhs41 - clhs146*mu;
const double clhs148 = pow(DN(2,1), 2);
const double clhs149 = clhs148*mu;
const double clhs150 = pow(DN(2,0), 2);
const double clhs151 = pow(N[2], 2);
const double clhs152 = DN(2,1)*clhs75;
const double clhs153 = N[2]*clhs62 + N[2]*clhs64 + clhs150*mu + clhs151*clhs7 - clhs152*mu;
const double clhs154 = N[2] - 1.0*clhs136*clhs19*rho - 1.0*clhs138*clhs19*rho;
const double clhs155 = DN(2,1)*clhs6;
const double clhs156 = clhs20*clhs60;
const double clhs157 = DN(2,1)*clhs20;
const double clhs158 = clhs20*clhs75;
lhs(0,0)=clhs1 + clhs16*clhs22 + clhs16*clhs25 - clhs16*clhs29 + clhs2*clhs6 + clhs32;
lhs(0,1)=clhs33*clhs34;
lhs(0,2)=DN(0,0)*(-clhs29 - clhs35);
lhs(0,3)=clhs22*clhs42 + clhs25*clhs42 - clhs29*clhs42 + clhs47 + clhs53;
lhs(0,4)=clhs34*clhs55;
lhs(0,5)=1.0*DN(1,0)*clhs18*clhs19*rho + 1.0*DN(1,0)*clhs19*clhs24*rho - clhs56 - clhs57*clhs58;
lhs(0,6)=clhs22*clhs65 + clhs25*clhs65 - clhs29*clhs65 + clhs70 + clhs74;
lhs(0,7)=clhs34*clhs76;
lhs(0,8)=1.0*DN(2,0)*clhs18*clhs19*rho + 1.0*DN(2,0)*clhs19*clhs24*rho - clhs58*clhs78 - clhs77;
lhs(1,0)=DN(0,1)*clhs34;
lhs(1,1)=2*clhs1 + clhs22*clhs82 + clhs25*clhs82 + clhs30*clhs80 + clhs32 + clhs33*clhs79 - clhs82*clhs83;
lhs(1,2)=DN(0,1)*(-clhs35 - clhs83);
lhs(1,3)=DN(1,0)*clhs79;
lhs(1,4)=clhs22*clhs85 + clhs25*clhs85 + clhs53 + clhs55*clhs79 - clhs83*clhs85 + clhs86;
lhs(1,5)=1.0*DN(1,1)*clhs18*clhs19*rho + 1.0*DN(1,1)*clhs19*clhs24*rho - clhs87 - clhs88;
lhs(1,6)=DN(2,0)*clhs79;
lhs(1,7)=clhs22*clhs89 + clhs25*clhs89 + clhs74 + clhs76*clhs79 - clhs83*clhs89 + clhs90;
lhs(1,8)=1.0*DN(2,1)*clhs18*clhs19*rho + 1.0*DN(2,1)*clhs19*clhs24*rho - clhs91 - clhs92;
lhs(2,0)=DN(0,0)*(N[0] + clhs19*(-1.0*clhs11 + 1.0*clhs13 + 1.0*clhs15 + 1.0*clhs8));
lhs(2,1)=DN(0,1)*N[0] + clhs30*clhs9 + clhs82*clhs93 - clhs82*clhs94;
lhs(2,2)=clhs20*(clhs0 + clhs2 - clhs31);
lhs(2,3)=clhs42*clhs95 + clhs48;
lhs(2,4)=clhs50 + clhs85*clhs93 - clhs85*clhs94 + clhs96;
lhs(2,5)=clhs20*(-clhs52 + clhs97);
lhs(2,6)=clhs65*clhs95 + clhs71;
lhs(2,7)=clhs72 + clhs89*clhs93 - clhs89*clhs94 + clhs98;
lhs(2,8)=clhs20*(-clhs73 + clhs99);
lhs(3,0)=clhs101*clhs16 + clhs103*clhs16 - clhs105*clhs16 + clhs107 + clhs47;
lhs(3,1)=clhs108*clhs33;
lhs(3,2)=1.0*DN(0,0)*clhs100*clhs19*rho + 1.0*DN(0,0)*clhs102*clhs19*rho - clhs109*clhs95 - clhs48;
lhs(3,3)=clhs101*clhs42 + clhs103*clhs42 - clhs105*clhs42 + clhs111 + clhs112*clhs6 + clhs115;
lhs(3,4)=clhs108*clhs55;
lhs(3,5)=DN(1,0)*(-clhs105 - clhs116);
lhs(3,6)=clhs101*clhs65 + clhs103*clhs65 - clhs105*clhs65 + clhs121 + clhs125;
lhs(3,7)=clhs108*clhs76;
lhs(3,8)=1.0*DN(2,0)*clhs100*clhs19*rho + 1.0*DN(2,0)*clhs102*clhs19*rho - clhs109*clhs78 - clhs126;
lhs(4,0)=DN(1,1)*clhs34;
lhs(4,1)=clhs101*clhs82 + clhs103*clhs82 + clhs107 + clhs127*clhs33 - clhs128*clhs82 + clhs86;
lhs(4,2)=1.0*DN(0,1)*clhs100*clhs19*rho + 1.0*DN(0,1)*clhs102*clhs19*rho - clhs50 - clhs88;
lhs(4,3)=DN(1,1)*clhs108;
lhs(4,4)=clhs101*clhs85 + clhs103*clhs85 + 2*clhs111 + clhs113*clhs80 + clhs115 + clhs127*clhs55 - clhs128*clhs85;
lhs(4,5)=DN(1,1)*(-clhs116 - clhs128);
lhs(4,6)=DN(2,0)*clhs127;
lhs(4,7)=clhs101*clhs89 + clhs103*clhs89 + clhs125 + clhs127*clhs76 - clhs128*clhs89 + clhs129;
lhs(4,8)=1.0*DN(2,1)*clhs100*clhs19*rho + 1.0*DN(2,1)*clhs102*clhs19*rho - clhs130 - clhs131;
lhs(5,0)=clhs16*clhs57 + clhs56;
lhs(5,1)=clhs132*clhs82 - clhs133*clhs82 + clhs87 + clhs96;
lhs(5,2)=clhs20*(-clhs106 + clhs97);
lhs(5,3)=DN(1,0)*(N[1] + clhs19*(1.0*clhs36 - 1.0*clhs37 + 1.0*clhs39 + 1.0*clhs41));
lhs(5,4)=DN(1,1)*N[1] + clhs113*clhs9 + clhs132*clhs85 - clhs133*clhs85;
lhs(5,5)=clhs20*(clhs110 + clhs112 - clhs114);
lhs(5,6)=clhs122 + clhs57*clhs65;
lhs(5,7)=clhs123 + clhs132*clhs89 - clhs133*clhs89 + clhs134;
lhs(5,8)=clhs20*(-clhs124 + clhs135);
lhs(6,0)=clhs137*clhs16 + clhs139*clhs16 - clhs141*clhs16 + clhs143 + clhs70;
lhs(6,1)=clhs144*clhs33;
lhs(6,2)=1.0*DN(0,0)*clhs136*clhs19*rho + 1.0*DN(0,0)*clhs138*clhs19*rho - clhs145*clhs95 - clhs71;
lhs(6,3)=clhs121 + clhs137*clhs42 + clhs139*clhs42 - clhs141*clhs42 + clhs147;
lhs(6,4)=clhs144*clhs55;
lhs(6,5)=1.0*DN(1,0)*clhs136*clhs19*rho + 1.0*DN(1,0)*clhs138*clhs19*rho - clhs122 - clhs145*clhs57;
lhs(6,6)=clhs137*clhs65 + clhs139*clhs65 - clhs141*clhs65 + clhs149 + clhs150*clhs6 + clhs153;
lhs(6,7)=clhs144*clhs76;
lhs(6,8)=DN(2,0)*(-clhs141 - clhs154);
lhs(7,0)=DN(2,1)*clhs34;
lhs(7,1)=clhs137*clhs82 + clhs139*clhs82 + clhs143 + clhs155*clhs33 - clhs156*clhs82 + clhs90;
lhs(7,2)=1.0*DN(0,1)*clhs136*clhs19*rho + 1.0*DN(0,1)*clhs138*clhs19*rho - clhs72 - clhs92;
lhs(7,3)=DN(2,1)*clhs108;
lhs(7,4)=clhs129 + clhs137*clhs85 + clhs139*clhs85 + clhs147 + clhs155*clhs55 - clhs156*clhs85;
lhs(7,5)=1.0*DN(1,1)*clhs136*clhs19*rho + 1.0*DN(1,1)*clhs138*clhs19*rho - clhs123 - clhs131;
lhs(7,6)=DN(2,1)*clhs144;
lhs(7,7)=clhs137*clhs89 + clhs139*clhs89 + 2*clhs149 + clhs151*clhs80 + clhs153 + clhs155*clhs76 - clhs156*clhs89;
lhs(7,8)=DN(2,1)*(-clhs154 - clhs156);
lhs(8,0)=clhs16*clhs78 + clhs77;
lhs(8,1)=clhs157*clhs82 - clhs158*clhs82 + clhs91 + clhs98;
lhs(8,2)=clhs20*(-clhs142 + clhs99);
lhs(8,3)=clhs126 + clhs42*clhs78;
lhs(8,4)=clhs130 + clhs134 + clhs157*clhs85 - clhs158*clhs85;
lhs(8,5)=clhs20*(clhs135 - clhs146);
lhs(8,6)=DN(2,0)*(N[2] + clhs19*(1.0*clhs59 - 1.0*clhs60 + 1.0*clhs62 + 1.0*clhs64));
lhs(8,7)=DN(2,1)*N[2] + clhs151*clhs9 + clhs157*clhs89 - clhs158*clhs89;
lhs(8,8)=clhs20*(clhs148 + clhs150 - clhs152);


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
const double clhs7 = bdf0*rho;
const double clhs8 = N[0]*clhs7;
const double clhs9 = 1.0/y;
const double clhs10 = clhs9*mu;
const double clhs11 = DN(0,1)*clhs10;
const double clhs12 = DN(0,0)*clhs3;
const double clhs13 = clhs12*rho;
const double clhs14 = DN(0,1)*clhs4;
const double clhs15 = clhs14*rho;
const double clhs16 = -clhs11 + clhs13 + clhs15 + clhs8;
const double clhs17 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double clhs18 = N[0]*clhs17 + clhs12;
const double clhs19 = 1.0/(clhs5/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs20 = 1.0*clhs19;
const double clhs21 = clhs20*rho;
const double clhs22 = clhs18*clhs21;
const double clhs23 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double clhs24 = N[0]*clhs23 + clhs14;
const double clhs25 = clhs21*clhs24;
const double clhs26 = N[0]*clhs9;
const double clhs27 = DN(0,1) - clhs26;
const double clhs28 = clhs10*clhs20;
const double clhs29 = clhs27*clhs28;
const double clhs30 = pow(N[0], 2);
const double clhs31 = DN(0,1)*clhs26;
const double clhs32 = N[0]*clhs13 + N[0]*clhs15 + clhs2*mu + clhs30*clhs7 - clhs31*mu;
const double clhs33 = DN(0,1) + clhs26;
const double clhs34 = DN(0,0)*clhs6;
const double clhs35 = N[0] - 1.0*clhs18*clhs19*rho - 1.0*clhs19*clhs24*rho;
const double clhs36 = N[1]*clhs7;
const double clhs37 = DN(1,1)*clhs10;
const double clhs38 = DN(1,0)*clhs3;
const double clhs39 = clhs38*rho;
const double clhs40 = DN(1,1)*clhs4;
const double clhs41 = clhs40*rho;
const double clhs42 = clhs36 - clhs37 + clhs39 + clhs41;
const double clhs43 = DN(0,0)*DN(1,0);
const double clhs44 = DN(0,1)*DN(1,1);
const double clhs45 = clhs44*mu;
const double clhs46 = N[1]*clhs8 + clhs43*mu;
const double clhs47 = clhs43*clhs6 + clhs45 + clhs46;
const double clhs48 = DN(1,0)*N[0];
const double clhs49 = clhs3*rho;
const double clhs50 = DN(1,1)*N[0];
const double clhs51 = clhs4*rho;
const double clhs52 = DN(1,1)*clhs26;
const double clhs53 = clhs48*clhs49 + clhs50*clhs51 - clhs52*mu;
const double clhs54 = N[1]*clhs9;
const double clhs55 = DN(1,1) + clhs54;
const double clhs56 = DN(0,0)*N[1];
const double clhs57 = DN(1,0)*clhs20;
const double clhs58 = clhs10*clhs27;
const double clhs59 = N[2]*clhs7;
const double clhs60 = DN(2,1)*clhs10;
const double clhs61 = DN(2,0)*clhs3;
const double clhs62 = clhs61*rho;
const double clhs63 = DN(2,1)*clhs4;
const double clhs64 = clhs63*rho;
const double clhs65 = clhs59 - clhs60 + clhs62 + clhs64;
const double clhs66 = DN(0,0)*DN(2,0);
const double clhs67 = DN(0,1)*DN(2,1);
const double clhs68 = clhs67*mu;
const double clhs69 = N[2]*clhs8 + clhs66*mu;
const double clhs70 = clhs6*clhs66 + clhs68 + clhs69;
const double clhs71 = DN(2,0)*N[0];
const double clhs72 = DN(2,1)*N[0];
const double clhs73 = DN(2,1)*clhs26;
const double clhs74 = clhs49*clhs71 + clhs51*clhs72 - clhs73*mu;
const double clhs75 = N[2]*clhs9;
const double clhs76 = DN(2,1) + clhs75;
const double clhs77 = DN(0,0)*N[2];
const double clhs78 = DN(2,0)*clhs20;
const double clhs79 = N[3]*clhs7;
const double clhs80 = DN(3,1)*clhs10;
const double clhs81 = DN(3,0)*clhs3;
const double clhs82 = clhs81*rho;
const double clhs83 = DN(3,1)*clhs4;
const double clhs84 = clhs83*rho;
const double clhs85 = clhs79 - clhs80 + clhs82 + clhs84;
const double clhs86 = DN(0,0)*DN(3,0);
const double clhs87 = DN(0,1)*DN(3,1);
const double clhs88 = clhs87*mu;
const double clhs89 = N[3]*clhs8 + clhs86*mu;
const double clhs90 = clhs6*clhs86 + clhs88 + clhs89;
const double clhs91 = DN(3,0)*N[0];
const double clhs92 = DN(3,1)*N[0];
const double clhs93 = DN(3,1)*clhs26;
const double clhs94 = clhs49*clhs91 + clhs51*clhs92 - clhs93*mu;
const double clhs95 = N[3]*clhs9;
const double clhs96 = DN(3,1) + clhs95;
const double clhs97 = DN(0,0)*N[3];
const double clhs98 = DN(3,0)*clhs20;
const double clhs99 = DN(0,1)*clhs6;
const double clhs100 = mu/pow(y, 2);
const double clhs101 = N[0]*clhs100;
const double clhs102 = clhs101 + clhs16;
const double clhs103 = clhs11*clhs20;
const double clhs104 = N[1]*clhs100;
const double clhs105 = clhs104 + clhs42;
const double clhs106 = N[1]*clhs101 + 2*clhs45 + clhs46;
const double clhs107 = DN(0,1)*N[1];
const double clhs108 = clhs28*clhs44;
const double clhs109 = N[2]*clhs100;
const double clhs110 = clhs109 + clhs65;
const double clhs111 = N[2]*clhs101 + 2*clhs68 + clhs69;
const double clhs112 = DN(0,1)*N[2];
const double clhs113 = clhs28*clhs67;
const double clhs114 = N[3]*clhs100 + clhs85;
const double clhs115 = N[3]*clhs101 + 2*clhs88 + clhs89;
const double clhs116 = DN(0,1)*N[3];
const double clhs117 = clhs28*clhs87;
const double clhs118 = DN(0,1)*clhs20;
const double clhs119 = clhs20*clhs26;
const double clhs120 = DN(0,0)*clhs20;
const double clhs121 = N[1]*clhs26;
const double clhs122 = clhs43 + clhs44;
const double clhs123 = N[2]*clhs26;
const double clhs124 = clhs66 + clhs67;
const double clhs125 = N[3]*clhs26;
const double clhs126 = clhs86 + clhs87;
const double clhs127 = N[1]*clhs17 + clhs38;
const double clhs128 = clhs127*clhs21;
const double clhs129 = N[1]*clhs23 + clhs40;
const double clhs130 = clhs129*clhs21;
const double clhs131 = DN(1,1) - clhs54;
const double clhs132 = clhs131*clhs28;
const double clhs133 = DN(0,1)*clhs54;
const double clhs134 = N[1]*clhs13 + N[1]*clhs15 - clhs133*mu;
const double clhs135 = DN(1,0)*clhs6;
const double clhs136 = clhs10*clhs131;
const double clhs137 = pow(DN(1,1), 2);
const double clhs138 = clhs137*mu;
const double clhs139 = pow(DN(1,0), 2);
const double clhs140 = pow(N[1], 2);
const double clhs141 = DN(1,1)*clhs54;
const double clhs142 = N[1]*clhs39 + N[1]*clhs41 + clhs139*mu + clhs140*clhs7 - clhs141*mu;
const double clhs143 = N[1] - 1.0*clhs127*clhs19*rho - 1.0*clhs129*clhs19*rho;
const double clhs144 = DN(1,0)*DN(2,0);
const double clhs145 = DN(1,1)*DN(2,1);
const double clhs146 = clhs145*mu;
const double clhs147 = N[2]*clhs36 + clhs144*mu;
const double clhs148 = clhs144*clhs6 + clhs146 + clhs147;
const double clhs149 = DN(2,0)*N[1];
const double clhs150 = DN(2,1)*N[1];
const double clhs151 = DN(2,1)*clhs54;
const double clhs152 = clhs149*clhs49 + clhs150*clhs51 - clhs151*mu;
const double clhs153 = DN(1,0)*N[2];
const double clhs154 = DN(1,0)*DN(3,0);
const double clhs155 = DN(1,1)*DN(3,1);
const double clhs156 = clhs155*mu;
const double clhs157 = N[3]*clhs36 + clhs154*mu;
const double clhs158 = clhs154*clhs6 + clhs156 + clhs157;
const double clhs159 = DN(3,0)*N[1];
const double clhs160 = DN(3,1)*N[1];
const double clhs161 = DN(3,1)*clhs54;
const double clhs162 = clhs159*clhs49 + clhs160*clhs51 - clhs161*mu;
const double clhs163 = DN(1,0)*N[3];
const double clhs164 = DN(1,1)*clhs6;
const double clhs165 = clhs20*clhs37;
const double clhs166 = N[2]*clhs104 + 2*clhs146 + clhs147;
const double clhs167 = DN(1,1)*N[2];
const double clhs168 = clhs145*clhs28;
const double clhs169 = N[3]*clhs104 + 2*clhs156 + clhs157;
const double clhs170 = DN(1,1)*N[3];
const double clhs171 = clhs155*clhs28;
const double clhs172 = DN(1,1)*clhs20;
const double clhs173 = clhs20*clhs54;
const double clhs174 = N[2]*clhs54;
const double clhs175 = clhs144 + clhs145;
const double clhs176 = N[3]*clhs54;
const double clhs177 = clhs154 + clhs155;
const double clhs178 = N[2]*clhs17 + clhs61;
const double clhs179 = clhs178*clhs21;
const double clhs180 = N[2]*clhs23 + clhs63;
const double clhs181 = clhs180*clhs21;
const double clhs182 = DN(2,1) - clhs75;
const double clhs183 = clhs182*clhs28;
const double clhs184 = DN(0,1)*clhs75;
const double clhs185 = N[2]*clhs13 + N[2]*clhs15 - clhs184*mu;
const double clhs186 = DN(2,0)*clhs6;
const double clhs187 = clhs10*clhs182;
const double clhs188 = DN(1,1)*clhs75;
const double clhs189 = N[2]*clhs39 + N[2]*clhs41 - clhs188*mu;
const double clhs190 = pow(DN(2,1), 2);
const double clhs191 = clhs190*mu;
const double clhs192 = pow(DN(2,0), 2);
const double clhs193 = pow(N[2], 2);
const double clhs194 = DN(2,1)*clhs75;
const double clhs195 = N[2]*clhs62 + N[2]*clhs64 + clhs192*mu + clhs193*clhs7 - clhs194*mu;
const double clhs196 = N[2] - 1.0*clhs178*clhs19*rho - 1.0*clhs180*clhs19*rho;
const double clhs197 = DN(2,0)*DN(3,0);
const double clhs198 = DN(2,1)*DN(3,1);
const double clhs199 = clhs198*mu;
const double clhs200 = N[3]*clhs59 + clhs197*mu;
const double clhs201 = clhs197*clhs6 + clhs199 + clhs200;
const double clhs202 = DN(3,0)*N[2];
const double clhs203 = DN(3,1)*N[2];
const double clhs204 = DN(3,1)*clhs75;
const double clhs205 = clhs202*clhs49 + clhs203*clhs51 - clhs204*mu;
const double clhs206 = DN(2,0)*N[3];
const double clhs207 = DN(2,1)*clhs6;
const double clhs208 = clhs20*clhs60;
const double clhs209 = N[3]*clhs109 + 2*clhs199 + clhs200;
const double clhs210 = DN(2,1)*N[3];
const double clhs211 = clhs198*clhs28;
const double clhs212 = DN(2,1)*clhs20;
const double clhs213 = clhs20*clhs75;
const double clhs214 = N[3]*clhs75;
const double clhs215 = clhs197 + clhs198;
const double clhs216 = N[3]*clhs17 + clhs81;
const double clhs217 = clhs21*clhs216;
const double clhs218 = N[3]*clhs23 + clhs83;
const double clhs219 = clhs21*clhs218;
const double clhs220 = DN(3,1) - clhs95;
const double clhs221 = clhs220*clhs28;
const double clhs222 = DN(0,1)*clhs95;
const double clhs223 = N[3]*clhs13 + N[3]*clhs15 - clhs222*mu;
const double clhs224 = DN(3,0)*clhs6;
const double clhs225 = clhs10*clhs220;
const double clhs226 = DN(1,1)*clhs95;
const double clhs227 = N[3]*clhs39 + N[3]*clhs41 - clhs226*mu;
const double clhs228 = DN(2,1)*clhs95;
const double clhs229 = N[3]*clhs62 + N[3]*clhs64 - clhs228*mu;
const double clhs230 = pow(DN(3,1), 2);
const double clhs231 = clhs230*mu;
const double clhs232 = pow(DN(3,0), 2);
const double clhs233 = pow(N[3], 2);
const double clhs234 = DN(3,1)*clhs95;
const double clhs235 = N[3]*clhs82 + N[3]*clhs84 + clhs232*mu + clhs233*clhs7 - clhs234*mu;
const double clhs236 = N[3] - 1.0*clhs19*clhs216*rho - 1.0*clhs19*clhs218*rho;
const double clhs237 = DN(3,1)*clhs6;
const double clhs238 = clhs20*clhs80;
const double clhs239 = DN(3,1)*clhs20;
const double clhs240 = clhs20*clhs95;
lhs(0,0)=clhs1 + clhs16*clhs22 + clhs16*clhs25 - clhs16*clhs29 + clhs2*clhs6 + clhs32;
lhs(0,1)=clhs33*clhs34;
lhs(0,2)=DN(0,0)*(-clhs29 - clhs35);
lhs(0,3)=clhs22*clhs42 + clhs25*clhs42 - clhs29*clhs42 + clhs47 + clhs53;
lhs(0,4)=clhs34*clhs55;
lhs(0,5)=1.0*DN(1,0)*clhs18*clhs19*rho + 1.0*DN(1,0)*clhs19*clhs24*rho - clhs56 - clhs57*clhs58;
lhs(0,6)=clhs22*clhs65 + clhs25*clhs65 - clhs29*clhs65 + clhs70 + clhs74;
lhs(0,7)=clhs34*clhs76;
lhs(0,8)=1.0*DN(2,0)*clhs18*clhs19*rho + 1.0*DN(2,0)*clhs19*clhs24*rho - clhs58*clhs78 - clhs77;
lhs(0,9)=clhs22*clhs85 + clhs25*clhs85 - clhs29*clhs85 + clhs90 + clhs94;
lhs(0,10)=clhs34*clhs96;
lhs(0,11)=1.0*DN(3,0)*clhs18*clhs19*rho + 1.0*DN(3,0)*clhs19*clhs24*rho - clhs58*clhs98 - clhs97;
lhs(1,0)=DN(0,1)*clhs34;
lhs(1,1)=2*clhs1 + clhs100*clhs30 - clhs102*clhs103 + clhs102*clhs22 + clhs102*clhs25 + clhs32 + clhs33*clhs99;
lhs(1,2)=DN(0,1)*(-clhs103 - clhs35);
lhs(1,3)=DN(1,0)*clhs99;
lhs(1,4)=-clhs103*clhs105 + clhs105*clhs22 + clhs105*clhs25 + clhs106 + clhs53 + clhs55*clhs99;
lhs(1,5)=1.0*DN(1,1)*clhs18*clhs19*rho + 1.0*DN(1,1)*clhs19*clhs24*rho - clhs107 - clhs108;
lhs(1,6)=DN(2,0)*clhs99;
lhs(1,7)=-clhs103*clhs110 + clhs110*clhs22 + clhs110*clhs25 + clhs111 + clhs74 + clhs76*clhs99;
lhs(1,8)=1.0*DN(2,1)*clhs18*clhs19*rho + 1.0*DN(2,1)*clhs19*clhs24*rho - clhs112 - clhs113;
lhs(1,9)=DN(3,0)*clhs99;
lhs(1,10)=-clhs103*clhs114 + clhs114*clhs22 + clhs114*clhs25 + clhs115 + clhs94 + clhs96*clhs99;
lhs(1,11)=1.0*DN(3,1)*clhs18*clhs19*rho + 1.0*DN(3,1)*clhs19*clhs24*rho - clhs116 - clhs117;
lhs(2,0)=DN(0,0)*(N[0] + clhs19*(-1.0*clhs11 + 1.0*clhs13 + 1.0*clhs15 + 1.0*clhs8));
lhs(2,1)=DN(0,1)*N[0] + clhs102*clhs118 - clhs102*clhs119 + clhs30*clhs9;
lhs(2,2)=clhs20*(clhs0 + clhs2 - clhs31);
lhs(2,3)=clhs120*clhs42 + clhs48;
lhs(2,4)=clhs105*clhs118 - clhs105*clhs119 + clhs121 + clhs50;
lhs(2,5)=clhs20*(clhs122 - clhs52);
lhs(2,6)=clhs120*clhs65 + clhs71;
lhs(2,7)=clhs110*clhs118 - clhs110*clhs119 + clhs123 + clhs72;
lhs(2,8)=clhs20*(clhs124 - clhs73);
lhs(2,9)=clhs120*clhs85 + clhs91;
lhs(2,10)=clhs114*clhs118 - clhs114*clhs119 + clhs125 + clhs92;
lhs(2,11)=clhs20*(clhs126 - clhs93);
lhs(3,0)=clhs128*clhs16 + clhs130*clhs16 - clhs132*clhs16 + clhs134 + clhs47;
lhs(3,1)=clhs135*clhs33;
lhs(3,2)=1.0*DN(0,0)*clhs127*clhs19*rho + 1.0*DN(0,0)*clhs129*clhs19*rho - clhs120*clhs136 - clhs48;
lhs(3,3)=clhs128*clhs42 + clhs130*clhs42 - clhs132*clhs42 + clhs138 + clhs139*clhs6 + clhs142;
lhs(3,4)=clhs135*clhs55;
lhs(3,5)=DN(1,0)*(-clhs132 - clhs143);
lhs(3,6)=clhs128*clhs65 + clhs130*clhs65 - clhs132*clhs65 + clhs148 + clhs152;
lhs(3,7)=clhs135*clhs76;
lhs(3,8)=1.0*DN(2,0)*clhs127*clhs19*rho + 1.0*DN(2,0)*clhs129*clhs19*rho - clhs136*clhs78 - clhs153;
lhs(3,9)=clhs128*clhs85 + clhs130*clhs85 - clhs132*clhs85 + clhs158 + clhs162;
lhs(3,10)=clhs135*clhs96;
lhs(3,11)=1.0*DN(3,0)*clhs127*clhs19*rho + 1.0*DN(3,0)*clhs129*clhs19*rho - clhs136*clhs98 - clhs163;
lhs(4,0)=DN(1,1)*clhs34;
lhs(4,1)=clhs102*clhs128 + clhs102*clhs130 - clhs102*clhs165 + clhs106 + clhs134 + clhs164*clhs33;
lhs(4,2)=1.0*DN(0,1)*clhs127*clhs19*rho + 1.0*DN(0,1)*clhs129*clhs19*rho - clhs108 - clhs50;
lhs(4,3)=DN(1,1)*clhs135;
lhs(4,4)=clhs100*clhs140 + clhs105*clhs128 + clhs105*clhs130 - clhs105*clhs165 + 2*clhs138 + clhs142 + clhs164*clhs55;
lhs(4,5)=DN(1,1)*(-clhs143 - clhs165);
lhs(4,6)=DN(2,0)*clhs164;
lhs(4,7)=clhs110*clhs128 + clhs110*clhs130 - clhs110*clhs165 + clhs152 + clhs164*clhs76 + clhs166;
lhs(4,8)=1.0*DN(2,1)*clhs127*clhs19*rho + 1.0*DN(2,1)*clhs129*clhs19*rho - clhs167 - clhs168;
lhs(4,9)=DN(3,0)*clhs164;
lhs(4,10)=clhs114*clhs128 + clhs114*clhs130 - clhs114*clhs165 + clhs162 + clhs164*clhs96 + clhs169;
lhs(4,11)=1.0*DN(3,1)*clhs127*clhs19*rho + 1.0*DN(3,1)*clhs129*clhs19*rho - clhs170 - clhs171;
lhs(5,0)=clhs16*clhs57 + clhs56;
lhs(5,1)=clhs102*clhs172 - clhs102*clhs173 + clhs107 + clhs121;
lhs(5,2)=clhs20*(clhs122 - clhs133);
lhs(5,3)=DN(1,0)*(N[1] + clhs19*(1.0*clhs36 - 1.0*clhs37 + 1.0*clhs39 + 1.0*clhs41));
lhs(5,4)=DN(1,1)*N[1] + clhs105*clhs172 - clhs105*clhs173 + clhs140*clhs9;
lhs(5,5)=clhs20*(clhs137 + clhs139 - clhs141);
lhs(5,6)=clhs149 + clhs57*clhs65;
lhs(5,7)=clhs110*clhs172 - clhs110*clhs173 + clhs150 + clhs174;
lhs(5,8)=clhs20*(-clhs151 + clhs175);
lhs(5,9)=clhs159 + clhs57*clhs85;
lhs(5,10)=clhs114*clhs172 - clhs114*clhs173 + clhs160 + clhs176;
lhs(5,11)=clhs20*(-clhs161 + clhs177);
lhs(6,0)=clhs16*clhs179 + clhs16*clhs181 - clhs16*clhs183 + clhs185 + clhs70;
lhs(6,1)=clhs186*clhs33;
lhs(6,2)=1.0*DN(0,0)*clhs178*clhs19*rho + 1.0*DN(0,0)*clhs180*clhs19*rho - clhs120*clhs187 - clhs71;
lhs(6,3)=clhs148 + clhs179*clhs42 + clhs181*clhs42 - clhs183*clhs42 + clhs189;
lhs(6,4)=clhs186*clhs55;
lhs(6,5)=1.0*DN(1,0)*clhs178*clhs19*rho + 1.0*DN(1,0)*clhs180*clhs19*rho - clhs149 - clhs187*clhs57;
lhs(6,6)=clhs179*clhs65 + clhs181*clhs65 - clhs183*clhs65 + clhs191 + clhs192*clhs6 + clhs195;
lhs(6,7)=clhs186*clhs76;
lhs(6,8)=DN(2,0)*(-clhs183 - clhs196);
lhs(6,9)=clhs179*clhs85 + clhs181*clhs85 - clhs183*clhs85 + clhs201 + clhs205;
lhs(6,10)=clhs186*clhs96;
lhs(6,11)=1.0*DN(3,0)*clhs178*clhs19*rho + 1.0*DN(3,0)*clhs180*clhs19*rho - clhs187*clhs98 - clhs206;
lhs(7,0)=DN(2,1)*clhs34;
lhs(7,1)=clhs102*clhs179 + clhs102*clhs181 - clhs102*clhs208 + clhs111 + clhs185 + clhs207*clhs33;
lhs(7,2)=1.0*DN(0,1)*clhs178*clhs19*rho + 1.0*DN(0,1)*clhs180*clhs19*rho - clhs113 - clhs72;
lhs(7,3)=DN(2,1)*clhs135;
lhs(7,4)=clhs105*clhs179 + clhs105*clhs181 - clhs105*clhs208 + clhs166 + clhs189 + clhs207*clhs55;
lhs(7,5)=1.0*DN(1,1)*clhs178*clhs19*rho + 1.0*DN(1,1)*clhs180*clhs19*rho - clhs150 - clhs168;
lhs(7,6)=DN(2,1)*clhs186;
lhs(7,7)=clhs100*clhs193 + clhs110*clhs179 + clhs110*clhs181 - clhs110*clhs208 + 2*clhs191 + clhs195 + clhs207*clhs76;
lhs(7,8)=DN(2,1)*(-clhs196 - clhs208);
lhs(7,9)=DN(3,0)*clhs207;
lhs(7,10)=clhs114*clhs179 + clhs114*clhs181 - clhs114*clhs208 + clhs205 + clhs207*clhs96 + clhs209;
lhs(7,11)=1.0*DN(3,1)*clhs178*clhs19*rho + 1.0*DN(3,1)*clhs180*clhs19*rho - clhs210 - clhs211;
lhs(8,0)=clhs16*clhs78 + clhs77;
lhs(8,1)=clhs102*clhs212 - clhs102*clhs213 + clhs112 + clhs123;
lhs(8,2)=clhs20*(clhs124 - clhs184);
lhs(8,3)=clhs153 + clhs42*clhs78;
lhs(8,4)=clhs105*clhs212 - clhs105*clhs213 + clhs167 + clhs174;
lhs(8,5)=clhs20*(clhs175 - clhs188);
lhs(8,6)=DN(2,0)*(N[2] + clhs19*(1.0*clhs59 - 1.0*clhs60 + 1.0*clhs62 + 1.0*clhs64));
lhs(8,7)=DN(2,1)*N[2] + clhs110*clhs212 - clhs110*clhs213 + clhs193*clhs9;
lhs(8,8)=clhs20*(clhs190 + clhs192 - clhs194);
lhs(8,9)=clhs202 + clhs78*clhs85;
lhs(8,10)=clhs114*clhs212 - clhs114*clhs213 + clhs203 + clhs214;
lhs(8,11)=clhs20*(-clhs204 + clhs215);
lhs(9,0)=clhs16*clhs217 + clhs16*clhs219 - clhs16*clhs221 + clhs223 + clhs90;
lhs(9,1)=clhs224*clhs33;
lhs(9,2)=1.0*DN(0,0)*clhs19*clhs216*rho + 1.0*DN(0,0)*clhs19*clhs218*rho - clhs120*clhs225 - clhs91;
lhs(9,3)=clhs158 + clhs217*clhs42 + clhs219*clhs42 - clhs221*clhs42 + clhs227;
lhs(9,4)=clhs224*clhs55;
lhs(9,5)=1.0*DN(1,0)*clhs19*clhs216*rho + 1.0*DN(1,0)*clhs19*clhs218*rho - clhs159 - clhs225*clhs57;
lhs(9,6)=clhs201 + clhs217*clhs65 + clhs219*clhs65 - clhs221*clhs65 + clhs229;
lhs(9,7)=clhs224*clhs76;
lhs(9,8)=1.0*DN(2,0)*clhs19*clhs216*rho + 1.0*DN(2,0)*clhs19*clhs218*rho - clhs202 - clhs225*clhs78;
lhs(9,9)=clhs217*clhs85 + clhs219*clhs85 - clhs221*clhs85 + clhs231 + clhs232*clhs6 + clhs235;
lhs(9,10)=clhs224*clhs96;
lhs(9,11)=DN(3,0)*(-clhs221 - clhs236);
lhs(10,0)=DN(3,1)*clhs34;
lhs(10,1)=clhs102*clhs217 + clhs102*clhs219 - clhs102*clhs238 + clhs115 + clhs223 + clhs237*clhs33;
lhs(10,2)=1.0*DN(0,1)*clhs19*clhs216*rho + 1.0*DN(0,1)*clhs19*clhs218*rho - clhs117 - clhs92;
lhs(10,3)=DN(3,1)*clhs135;
lhs(10,4)=clhs105*clhs217 + clhs105*clhs219 - clhs105*clhs238 + clhs169 + clhs227 + clhs237*clhs55;
lhs(10,5)=1.0*DN(1,1)*clhs19*clhs216*rho + 1.0*DN(1,1)*clhs19*clhs218*rho - clhs160 - clhs171;
lhs(10,6)=DN(3,1)*clhs186;
lhs(10,7)=clhs110*clhs217 + clhs110*clhs219 - clhs110*clhs238 + clhs209 + clhs229 + clhs237*clhs76;
lhs(10,8)=1.0*DN(2,1)*clhs19*clhs216*rho + 1.0*DN(2,1)*clhs19*clhs218*rho - clhs203 - clhs211;
lhs(10,9)=DN(3,1)*clhs224;
lhs(10,10)=clhs100*clhs233 + clhs114*clhs217 + clhs114*clhs219 - clhs114*clhs238 + 2*clhs231 + clhs235 + clhs237*clhs96;
lhs(10,11)=DN(3,1)*(-clhs236 - clhs238);
lhs(11,0)=clhs16*clhs98 + clhs97;
lhs(11,1)=clhs102*clhs239 - clhs102*clhs240 + clhs116 + clhs125;
lhs(11,2)=clhs20*(clhs126 - clhs222);
lhs(11,3)=clhs163 + clhs42*clhs98;
lhs(11,4)=clhs105*clhs239 - clhs105*clhs240 + clhs170 + clhs176;
lhs(11,5)=clhs20*(clhs177 - clhs226);
lhs(11,6)=clhs206 + clhs65*clhs98;
lhs(11,7)=clhs110*clhs239 - clhs110*clhs240 + clhs210 + clhs214;
lhs(11,8)=clhs20*(clhs215 - clhs228);
lhs(11,9)=DN(3,0)*(N[3] + clhs19*(1.0*clhs79 - 1.0*clhs80 + 1.0*clhs82 + 1.0*clhs84));
lhs(11,10)=DN(3,1)*N[3] + clhs114*clhs239 - clhs114*clhs240 + clhs233*clhs9;
lhs(11,11)=clhs20*(clhs230 + clhs232 - clhs234);


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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs2 = DN(0,0)*mu;
const double crhs3 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double crhs4 = crhs3*mu;
const double crhs5 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs6 = 1.0/y;
const double crhs7 = N[0]*crhs1;
const double crhs8 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double crhs9 = crhs8*rho;
const double crhs10 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double crhs11 = crhs10*rho;
const double crhs12 = crhs11*crhs3;
const double crhs13 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
const double crhs14 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2));
const double crhs15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs16 = crhs15*crhs6;
const double crhs17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs18 = (crhs14*h/stab_c1 + mu)*(crhs1 + crhs16 + crhs17);
const double crhs19 = 1.0/(crhs14/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs20 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs1*crhs9 + crhs12 + crhs13 - crhs4*crhs6 - crhs5*rho;
const double crhs21 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double crhs22 = DN(0,0)*crhs8 + N[0]*crhs21;
const double crhs23 = 1.0*crhs19;
const double crhs24 = crhs20*crhs23;
const double crhs25 = crhs24*rho;
const double crhs26 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double crhs27 = DN(0,1)*crhs10 + N[0]*crhs26;
const double crhs28 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs29 = crhs17*mu;
const double crhs30 = 2*crhs29;
const double crhs31 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs32 = crhs15*mu/pow(y, 2);
const double crhs33 = crhs28*crhs9;
const double crhs34 = N[0]*crhs17;
const double crhs35 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
const double crhs36 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs11*crhs17 - crhs29*crhs6 - crhs31*rho + crhs32 + crhs33 + crhs35;
const double crhs37 = crhs23*crhs36;
const double crhs38 = crhs37*rho;
const double crhs39 = DN(1,0)*mu;
const double crhs40 = N[1]*crhs1;
const double crhs41 = DN(1,0)*crhs8 + N[1]*crhs21;
const double crhs42 = DN(1,1)*crhs10 + N[1]*crhs26;
const double crhs43 = N[1]*crhs17;
const double crhs44 = DN(2,0)*mu;
const double crhs45 = N[2]*crhs1;
const double crhs46 = DN(2,0)*crhs8 + N[2]*crhs21;
const double crhs47 = DN(2,1)*crhs10 + N[2]*crhs26;
const double crhs48 = N[2]*crhs17;
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,1)*crhs4 - N[0]*crhs12 - N[0]*crhs13 + N[0]*crhs3*crhs6*mu + N[0]*crhs5*rho - crhs1*crhs2 + 1.0*crhs19*crhs20*crhs6*mu*(DN(0,1) - N[0]*crhs6) - crhs22*crhs25 - crhs25*crhs27 - crhs7*crhs9;
rhs[1]=DN(0,1)*crhs0 - DN(0,1)*crhs18 + 1.0*DN(0,1)*crhs19*crhs36*crhs6*mu - DN(0,1)*crhs30 + N[0]*crhs17*crhs6*mu + N[0]*crhs31*rho - N[0]*crhs32 - N[0]*crhs33 - N[0]*crhs35 - crhs11*crhs34 - crhs2*crhs28 - crhs22*crhs38 - crhs27*crhs38;
rhs[2]=-DN(0,0)*crhs24 - DN(0,1)*crhs37 - N[0]*crhs16 + 1.0*N[0]*crhs19*crhs36*crhs6 - crhs34 - crhs7;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,1)*crhs4 - N[1]*crhs12 - N[1]*crhs13 + N[1]*crhs3*crhs6*mu + N[1]*crhs5*rho - crhs1*crhs39 + 1.0*crhs19*crhs20*crhs6*mu*(DN(1,1) - N[1]*crhs6) - crhs25*crhs41 - crhs25*crhs42 - crhs40*crhs9;
rhs[4]=DN(1,1)*crhs0 - DN(1,1)*crhs18 + 1.0*DN(1,1)*crhs19*crhs36*crhs6*mu - DN(1,1)*crhs30 + N[1]*crhs17*crhs6*mu + N[1]*crhs31*rho - N[1]*crhs32 - N[1]*crhs33 - N[1]*crhs35 - crhs11*crhs43 - crhs28*crhs39 - crhs38*crhs41 - crhs38*crhs42;
rhs[5]=-DN(1,0)*crhs24 - DN(1,1)*crhs37 - N[1]*crhs16 + 1.0*N[1]*crhs19*crhs36*crhs6 - crhs40 - crhs43;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,1)*crhs4 - N[2]*crhs12 - N[2]*crhs13 + N[2]*crhs3*crhs6*mu + N[2]*crhs5*rho - crhs1*crhs44 + 1.0*crhs19*crhs20*crhs6*mu*(DN(2,1) - N[2]*crhs6) - crhs25*crhs46 - crhs25*crhs47 - crhs45*crhs9;
rhs[7]=DN(2,1)*crhs0 - DN(2,1)*crhs18 + 1.0*DN(2,1)*crhs19*crhs36*crhs6*mu - DN(2,1)*crhs30 + N[2]*crhs17*crhs6*mu + N[2]*crhs31*rho - N[2]*crhs32 - N[2]*crhs33 - N[2]*crhs35 - crhs11*crhs48 - crhs28*crhs44 - crhs38*crhs46 - crhs38*crhs47;
rhs[8]=-DN(2,0)*crhs24 - DN(2,1)*crhs37 - N[2]*crhs16 + 1.0*N[2]*crhs19*crhs36*crhs6 - crhs45 - crhs48;


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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs2 = DN(0,0)*mu;
const double crhs3 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double crhs4 = crhs3*mu;
const double crhs5 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs6 = 1.0/y;
const double crhs7 = N[0]*crhs1;
const double crhs8 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double crhs9 = crhs8*rho;
const double crhs10 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double crhs11 = crhs10*rho;
const double crhs12 = crhs11*crhs3;
const double crhs13 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
const double crhs14 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2));
const double crhs15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs16 = crhs15*crhs6;
const double crhs17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs18 = (crhs14*h/stab_c1 + mu)*(crhs1 + crhs16 + crhs17);
const double crhs19 = 1.0/(crhs14/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs20 = DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs1*crhs9 + crhs12 + crhs13 - crhs4*crhs6 - crhs5*rho;
const double crhs21 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double crhs22 = DN(0,0)*crhs8 + N[0]*crhs21;
const double crhs23 = 1.0*crhs19;
const double crhs24 = crhs20*crhs23;
const double crhs25 = crhs24*rho;
const double crhs26 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double crhs27 = DN(0,1)*crhs10 + N[0]*crhs26;
const double crhs28 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs29 = crhs17*mu;
const double crhs30 = 2*crhs29;
const double crhs31 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs32 = crhs15*mu/pow(y, 2);
const double crhs33 = crhs28*crhs9;
const double crhs34 = N[0]*crhs17;
const double crhs35 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
const double crhs36 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs11*crhs17 - crhs29*crhs6 - crhs31*rho + crhs32 + crhs33 + crhs35;
const double crhs37 = crhs23*crhs36;
const double crhs38 = crhs37*rho;
const double crhs39 = DN(1,0)*mu;
const double crhs40 = N[1]*crhs1;
const double crhs41 = DN(1,0)*crhs8 + N[1]*crhs21;
const double crhs42 = DN(1,1)*crhs10 + N[1]*crhs26;
const double crhs43 = N[1]*crhs17;
const double crhs44 = DN(2,0)*mu;
const double crhs45 = N[2]*crhs1;
const double crhs46 = DN(2,0)*crhs8 + N[2]*crhs21;
const double crhs47 = DN(2,1)*crhs10 + N[2]*crhs26;
const double crhs48 = N[2]*crhs17;
const double crhs49 = DN(3,0)*mu;
const double crhs50 = N[3]*crhs1;
const double crhs51 = DN(3,0)*crhs8 + N[3]*crhs21;
const double crhs52 = DN(3,1)*crhs10 + N[3]*crhs26;
const double crhs53 = N[3]*crhs17;
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,1)*crhs4 - N[0]*crhs12 - N[0]*crhs13 + N[0]*crhs3*crhs6*mu + N[0]*crhs5*rho - crhs1*crhs2 + 1.0*crhs19*crhs20*crhs6*mu*(DN(0,1) - N[0]*crhs6) - crhs22*crhs25 - crhs25*crhs27 - crhs7*crhs9;
rhs[1]=DN(0,1)*crhs0 - DN(0,1)*crhs18 + 1.0*DN(0,1)*crhs19*crhs36*crhs6*mu - DN(0,1)*crhs30 + N[0]*crhs17*crhs6*mu + N[0]*crhs31*rho - N[0]*crhs32 - N[0]*crhs33 - N[0]*crhs35 - crhs11*crhs34 - crhs2*crhs28 - crhs22*crhs38 - crhs27*crhs38;
rhs[2]=-DN(0,0)*crhs24 - DN(0,1)*crhs37 - N[0]*crhs16 + 1.0*N[0]*crhs19*crhs36*crhs6 - crhs34 - crhs7;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,1)*crhs4 - N[1]*crhs12 - N[1]*crhs13 + N[1]*crhs3*crhs6*mu + N[1]*crhs5*rho - crhs1*crhs39 + 1.0*crhs19*crhs20*crhs6*mu*(DN(1,1) - N[1]*crhs6) - crhs25*crhs41 - crhs25*crhs42 - crhs40*crhs9;
rhs[4]=DN(1,1)*crhs0 - DN(1,1)*crhs18 + 1.0*DN(1,1)*crhs19*crhs36*crhs6*mu - DN(1,1)*crhs30 + N[1]*crhs17*crhs6*mu + N[1]*crhs31*rho - N[1]*crhs32 - N[1]*crhs33 - N[1]*crhs35 - crhs11*crhs43 - crhs28*crhs39 - crhs38*crhs41 - crhs38*crhs42;
rhs[5]=-DN(1,0)*crhs24 - DN(1,1)*crhs37 - N[1]*crhs16 + 1.0*N[1]*crhs19*crhs36*crhs6 - crhs40 - crhs43;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,1)*crhs4 - N[2]*crhs12 - N[2]*crhs13 + N[2]*crhs3*crhs6*mu + N[2]*crhs5*rho - crhs1*crhs44 + 1.0*crhs19*crhs20*crhs6*mu*(DN(2,1) - N[2]*crhs6) - crhs25*crhs46 - crhs25*crhs47 - crhs45*crhs9;
rhs[7]=DN(2,1)*crhs0 - DN(2,1)*crhs18 + 1.0*DN(2,1)*crhs19*crhs36*crhs6*mu - DN(2,1)*crhs30 + N[2]*crhs17*crhs6*mu + N[2]*crhs31*rho - N[2]*crhs32 - N[2]*crhs33 - N[2]*crhs35 - crhs11*crhs48 - crhs28*crhs44 - crhs38*crhs46 - crhs38*crhs47;
rhs[8]=-DN(2,0)*crhs24 - DN(2,1)*crhs37 - N[2]*crhs16 + 1.0*N[2]*crhs19*crhs36*crhs6 - crhs45 - crhs48;
rhs[9]=DN(3,0)*crhs0 - DN(3,0)*crhs18 - DN(3,1)*crhs4 - N[3]*crhs12 - N[3]*crhs13 + N[3]*crhs3*crhs6*mu + N[3]*crhs5*rho - crhs1*crhs49 + 1.0*crhs19*crhs20*crhs6*mu*(DN(3,1) - N[3]*crhs6) - crhs25*crhs51 - crhs25*crhs52 - crhs50*crhs9;
rhs[10]=DN(3,1)*crhs0 - DN(3,1)*crhs18 + 1.0*DN(3,1)*crhs19*crhs36*crhs6*mu - DN(3,1)*crhs30 + N[3]*crhs17*crhs6*mu + N[3]*crhs31*rho - N[3]*crhs32 - N[3]*crhs33 - N[3]*crhs35 - crhs11*crhs53 - crhs28*crhs49 - crhs38*crhs51 - crhs38*crhs52;
rhs[11]=-DN(3,0)*crhs24 - DN(3,1)*crhs37 - N[3]*crhs16 + 1.0*N[3]*crhs19*crhs36*crhs6 - crhs50 - crhs53;


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