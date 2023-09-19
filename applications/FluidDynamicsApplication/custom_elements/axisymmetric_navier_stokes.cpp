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
const double clhs4 = 2*clhs3*h/stab_c1 + 2*mu;
const double clhs5 = DN(0,0)*clhs1;
const double clhs6 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double clhs7 = rho*(N[0]*clhs6 + clhs5);
const double clhs8 = bdf0*rho;
const double clhs9 = N[0]*clhs8;
const double clhs10 = 1.0/y;
const double clhs11 = DN(0,1)*mu;
const double clhs12 = clhs10*clhs11;
const double clhs13 = clhs5*rho;
const double clhs14 = DN(0,1)*clhs2;
const double clhs15 = clhs14*rho;
const double clhs16 = -clhs12 + clhs13 + clhs15 + clhs9;
const double clhs17 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = 1.0*clhs17;
const double clhs19 = clhs16*clhs18;
const double clhs20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double clhs21 = rho*(N[0]*clhs20 + clhs14);
const double clhs22 = pow(DN(0,1), 2);
const double clhs23 = pow(N[0], 2);
const double clhs24 = N[0]*clhs10;
const double clhs25 = N[0]*clhs13 + N[0]*clhs15 + clhs11*clhs24 + clhs22*mu + clhs23*clhs8;
const double clhs26 = DN(0,0)*mu;
const double clhs27 = DN(0,1) + clhs24;
const double clhs28 = clhs26 + clhs27*clhs4;
const double clhs29 = DN(0,0)*clhs18;
const double clhs30 = N[1]*clhs8;
const double clhs31 = clhs10*mu;
const double clhs32 = DN(1,1)*clhs31;
const double clhs33 = DN(1,0)*clhs1;
const double clhs34 = clhs33*rho;
const double clhs35 = DN(1,1)*clhs2;
const double clhs36 = clhs35*rho;
const double clhs37 = clhs30 - clhs32 + clhs34 + clhs36;
const double clhs38 = clhs18*clhs37;
const double clhs39 = DN(0,0)*DN(1,0);
const double clhs40 = DN(1,1)*clhs11 + N[1]*clhs9;
const double clhs41 = DN(1,0)*clhs26 + clhs39*clhs4 + clhs40;
const double clhs42 = DN(1,0)*N[0];
const double clhs43 = clhs1*rho;
const double clhs44 = DN(1,1)*clhs24;
const double clhs45 = DN(1,1)*N[0];
const double clhs46 = clhs2*rho;
const double clhs47 = clhs42*clhs43 + clhs44*mu + clhs45*clhs46;
const double clhs48 = DN(1,0)*mu;
const double clhs49 = N[1]*clhs10;
const double clhs50 = DN(1,1) + clhs49;
const double clhs51 = clhs4*clhs50 + clhs48;
const double clhs52 = DN(1,0)*clhs18;
const double clhs53 = N[2]*clhs8;
const double clhs54 = DN(2,1)*clhs31;
const double clhs55 = DN(2,0)*clhs1;
const double clhs56 = clhs55*rho;
const double clhs57 = DN(2,1)*clhs2;
const double clhs58 = clhs57*rho;
const double clhs59 = clhs53 - clhs54 + clhs56 + clhs58;
const double clhs60 = clhs18*clhs59;
const double clhs61 = DN(0,0)*DN(2,0);
const double clhs62 = DN(2,1)*clhs11 + N[2]*clhs9;
const double clhs63 = DN(2,0)*clhs26 + clhs4*clhs61 + clhs62;
const double clhs64 = DN(2,0)*N[0];
const double clhs65 = DN(2,1)*clhs24;
const double clhs66 = DN(2,1)*N[0];
const double clhs67 = clhs43*clhs64 + clhs46*clhs66 + clhs65*mu;
const double clhs68 = N[2]*clhs10;
const double clhs69 = DN(2,1) + clhs68;
const double clhs70 = DN(2,0)*mu + clhs4*clhs69;
const double clhs71 = DN(2,0)*clhs18;
const double clhs72 = mu/pow(y, 2);
const double clhs73 = N[0]*clhs72;
const double clhs74 = clhs18*(clhs16 + clhs73);
const double clhs75 = clhs18*(clhs21 + clhs7 - clhs73);
const double clhs76 = N[1]*clhs72;
const double clhs77 = clhs18*(clhs37 + clhs76);
const double clhs78 = N[1]*clhs73 + clhs40;
const double clhs79 = N[2]*clhs72;
const double clhs80 = clhs18*(clhs59 + clhs79);
const double clhs81 = N[2]*clhs73 + clhs62;
const double clhs82 = DN(0,1)*y;
const double clhs83 = DN(1,1)*clhs82 + clhs39;
const double clhs84 = DN(2,1)*clhs82 + clhs61;
const double clhs85 = rho*(N[1]*clhs6 + clhs33);
const double clhs86 = rho*(N[1]*clhs20 + clhs35);
const double clhs87 = N[1]*clhs13 + N[1]*clhs15 + clhs11*clhs49;
const double clhs88 = DN(0,0)*N[1];
const double clhs89 = pow(DN(1,0), 2);
const double clhs90 = pow(DN(1,1), 2);
const double clhs91 = pow(N[1], 2);
const double clhs92 = DN(1,1)*clhs49;
const double clhs93 = N[1]*clhs34 + N[1]*clhs36 + clhs8*clhs91 + clhs90*mu + clhs92*mu;
const double clhs94 = DN(1,0)*DN(2,0);
const double clhs95 = DN(1,1)*DN(2,1);
const double clhs96 = N[2]*clhs30 + clhs95*mu;
const double clhs97 = DN(2,0)*clhs48 + clhs4*clhs94 + clhs96;
const double clhs98 = DN(2,0)*N[1];
const double clhs99 = DN(2,1)*clhs49;
const double clhs100 = DN(2,1)*N[1];
const double clhs101 = clhs100*clhs46 + clhs43*clhs98 + clhs99*mu;
const double clhs102 = clhs18*(-clhs76 + clhs85 + clhs86);
const double clhs103 = N[2]*clhs76 + clhs96;
const double clhs104 = DN(1,1)*y;
const double clhs105 = clhs94 + clhs95*y;
const double clhs106 = rho*(N[2]*clhs6 + clhs55);
const double clhs107 = rho*(N[2]*clhs20 + clhs57);
const double clhs108 = N[2]*clhs13 + N[2]*clhs15 + clhs11*clhs68;
const double clhs109 = DN(0,0)*N[2];
const double clhs110 = DN(1,1)*clhs68;
const double clhs111 = N[2]*clhs34 + N[2]*clhs36 + clhs110*mu;
const double clhs112 = DN(1,0)*N[2];
const double clhs113 = pow(DN(2,0), 2);
const double clhs114 = pow(DN(2,1), 2);
const double clhs115 = pow(N[2], 2);
const double clhs116 = DN(2,1)*clhs68;
const double clhs117 = N[2]*clhs56 + N[2]*clhs58 + clhs114*mu + clhs115*clhs8 + clhs116*mu;
const double clhs118 = clhs18*(clhs106 + clhs107 - clhs79);
const double clhs119 = DN(2,1)*y;
lhs(0,0)=clhs0*clhs4 + clhs0*mu + clhs19*clhs21 + clhs19*clhs7 + clhs25;
lhs(0,1)=DN(0,0)*clhs28;
lhs(0,2)=DN(0,0)*N[0] + DN(0,1)*N[0] + clhs21*clhs29 + clhs29*clhs7;
lhs(0,3)=clhs21*clhs38 + clhs38*clhs7 + clhs41 + clhs47;
lhs(0,4)=DN(0,0)*clhs51;
lhs(0,5)=clhs21*clhs52 + clhs42 + clhs45 + clhs52*clhs7;
lhs(0,6)=clhs21*clhs60 + clhs60*clhs7 + clhs63 + clhs67;
lhs(0,7)=DN(0,0)*clhs70;
lhs(0,8)=clhs21*clhs71 + clhs64 + clhs66 + clhs7*clhs71;
lhs(1,0)=0;
lhs(1,1)=clhs21*clhs74 + clhs23*clhs72 + clhs25 + clhs7*clhs74 - clhs73*clhs74;
lhs(1,2)=DN(0,1)*clhs75;
lhs(1,3)=0;
lhs(1,4)=clhs21*clhs77 + clhs47 + clhs7*clhs77 - clhs73*clhs77 + clhs78;
lhs(1,5)=DN(1,1)*clhs75;
lhs(1,6)=0;
lhs(1,7)=clhs21*clhs80 + clhs67 + clhs7*clhs80 - clhs73*clhs80 + clhs81;
lhs(1,8)=DN(2,1)*clhs75;
lhs(2,0)=DN(0,0)*(N[0] + clhs17*(-1.0*clhs12 + 1.0*clhs13 + 1.0*clhs15 + 1.0*clhs9));
lhs(2,1)=N[0]*clhs27 - clhs24*clhs74 + clhs74*clhs82;
lhs(2,2)=clhs18*(-DN(0,1)*clhs24 + clhs0 + clhs22*y);
lhs(2,3)=clhs29*clhs37 + clhs42;
lhs(2,4)=N[0]*clhs50 - clhs24*clhs77 + clhs77*clhs82;
lhs(2,5)=clhs18*(-clhs44 + clhs83);
lhs(2,6)=clhs29*clhs59 + clhs64;
lhs(2,7)=N[0]*clhs69 - clhs24*clhs80 + clhs80*clhs82;
lhs(2,8)=clhs18*(-clhs65 + clhs84);
lhs(3,0)=clhs19*clhs85 + clhs19*clhs86 + clhs41 + clhs87;
lhs(3,1)=DN(1,0)*clhs28;
lhs(3,2)=DN(0,1)*N[1] + clhs29*clhs85 + clhs29*clhs86 + clhs88;
lhs(3,3)=clhs38*clhs85 + clhs38*clhs86 + clhs4*clhs89 + clhs89*mu + clhs93;
lhs(3,4)=DN(1,0)*clhs51;
lhs(3,5)=DN(1,0)*N[1] + DN(1,1)*N[1] + clhs52*clhs85 + clhs52*clhs86;
lhs(3,6)=clhs101 + clhs60*clhs85 + clhs60*clhs86 + clhs97;
lhs(3,7)=DN(1,0)*clhs70;
lhs(3,8)=clhs100 + clhs71*clhs85 + clhs71*clhs86 + clhs98;
lhs(4,0)=0;
lhs(4,1)=-clhs74*clhs76 + clhs74*clhs85 + clhs74*clhs86 + clhs78 + clhs87;
lhs(4,2)=DN(0,1)*clhs102;
lhs(4,3)=0;
lhs(4,4)=clhs72*clhs91 - clhs76*clhs77 + clhs77*clhs85 + clhs77*clhs86 + clhs93;
lhs(4,5)=DN(1,1)*clhs102;
lhs(4,6)=0;
lhs(4,7)=clhs101 + clhs103 - clhs76*clhs80 + clhs80*clhs85 + clhs80*clhs86;
lhs(4,8)=DN(2,1)*clhs102;
lhs(5,0)=clhs16*clhs52 + clhs88;
lhs(5,1)=N[1]*clhs27 + clhs104*clhs74 - clhs49*clhs74;
lhs(5,2)=clhs18*(-DN(0,1)*clhs49 + clhs83);
lhs(5,3)=DN(1,0)*(N[1] + clhs17*(1.0*clhs30 - 1.0*clhs32 + 1.0*clhs34 + 1.0*clhs36));
lhs(5,4)=N[1]*clhs50 + clhs104*clhs77 - clhs49*clhs77;
lhs(5,5)=clhs18*(clhs89 + clhs90*y - clhs92);
lhs(5,6)=clhs52*clhs59 + clhs98;
lhs(5,7)=N[1]*clhs69 + clhs104*clhs80 - clhs49*clhs80;
lhs(5,8)=clhs18*(clhs105 - clhs99);
lhs(6,0)=clhs106*clhs19 + clhs107*clhs19 + clhs108 + clhs63;
lhs(6,1)=DN(2,0)*clhs28;
lhs(6,2)=DN(0,1)*N[2] + clhs106*clhs29 + clhs107*clhs29 + clhs109;
lhs(6,3)=clhs106*clhs38 + clhs107*clhs38 + clhs111 + clhs97;
lhs(6,4)=DN(2,0)*clhs51;
lhs(6,5)=DN(1,1)*N[2] + clhs106*clhs52 + clhs107*clhs52 + clhs112;
lhs(6,6)=clhs106*clhs60 + clhs107*clhs60 + clhs113*clhs4 + clhs113*mu + clhs117;
lhs(6,7)=DN(2,0)*clhs70;
lhs(6,8)=DN(2,0)*N[2] + DN(2,1)*N[2] + clhs106*clhs71 + clhs107*clhs71;
lhs(7,0)=0;
lhs(7,1)=clhs106*clhs74 + clhs107*clhs74 + clhs108 - clhs74*clhs79 + clhs81;
lhs(7,2)=DN(0,1)*clhs118;
lhs(7,3)=0;
lhs(7,4)=clhs103 + clhs106*clhs77 + clhs107*clhs77 + clhs111 - clhs77*clhs79;
lhs(7,5)=DN(1,1)*clhs118;
lhs(7,6)=0;
lhs(7,7)=clhs106*clhs80 + clhs107*clhs80 + clhs115*clhs72 + clhs117 - clhs79*clhs80;
lhs(7,8)=DN(2,1)*clhs118;
lhs(8,0)=clhs109 + clhs16*clhs71;
lhs(8,1)=N[2]*clhs27 + clhs119*clhs74 - clhs68*clhs74;
lhs(8,2)=clhs18*(-DN(0,1)*clhs68 + clhs84);
lhs(8,3)=clhs112 + clhs37*clhs71;
lhs(8,4)=N[2]*clhs50 + clhs119*clhs77 - clhs68*clhs77;
lhs(8,5)=clhs18*(clhs105 - clhs110);
lhs(8,6)=DN(2,0)*(N[2] + clhs17*(1.0*clhs53 - 1.0*clhs54 + 1.0*clhs56 + 1.0*clhs58));
lhs(8,7)=N[2]*clhs69 + clhs119*clhs80 - clhs68*clhs80;
lhs(8,8)=clhs18*(clhs113 + clhs114*y - clhs116);


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
const double clhs4 = 2*clhs3*h/stab_c1 + 2*mu;
const double clhs5 = DN(0,0)*clhs1;
const double clhs6 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double clhs7 = rho*(N[0]*clhs6 + clhs5);
const double clhs8 = bdf0*rho;
const double clhs9 = N[0]*clhs8;
const double clhs10 = 1.0/y;
const double clhs11 = DN(0,1)*mu;
const double clhs12 = clhs10*clhs11;
const double clhs13 = clhs5*rho;
const double clhs14 = DN(0,1)*clhs2;
const double clhs15 = clhs14*rho;
const double clhs16 = -clhs12 + clhs13 + clhs15 + clhs9;
const double clhs17 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = 1.0*clhs17;
const double clhs19 = clhs16*clhs18;
const double clhs20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double clhs21 = rho*(N[0]*clhs20 + clhs14);
const double clhs22 = pow(DN(0,1), 2);
const double clhs23 = pow(N[0], 2);
const double clhs24 = N[0]*clhs10;
const double clhs25 = N[0]*clhs13 + N[0]*clhs15 + clhs11*clhs24 + clhs22*mu + clhs23*clhs8;
const double clhs26 = DN(0,0)*mu;
const double clhs27 = DN(0,1) + clhs24;
const double clhs28 = clhs26 + clhs27*clhs4;
const double clhs29 = DN(0,0)*clhs18;
const double clhs30 = N[1]*clhs8;
const double clhs31 = clhs10*mu;
const double clhs32 = DN(1,1)*clhs31;
const double clhs33 = DN(1,0)*clhs1;
const double clhs34 = clhs33*rho;
const double clhs35 = DN(1,1)*clhs2;
const double clhs36 = clhs35*rho;
const double clhs37 = clhs30 - clhs32 + clhs34 + clhs36;
const double clhs38 = clhs18*clhs37;
const double clhs39 = DN(0,0)*DN(1,0);
const double clhs40 = DN(1,1)*clhs11 + N[1]*clhs9;
const double clhs41 = DN(1,0)*clhs26 + clhs39*clhs4 + clhs40;
const double clhs42 = DN(1,0)*N[0];
const double clhs43 = clhs1*rho;
const double clhs44 = DN(1,1)*clhs24;
const double clhs45 = DN(1,1)*N[0];
const double clhs46 = clhs2*rho;
const double clhs47 = clhs42*clhs43 + clhs44*mu + clhs45*clhs46;
const double clhs48 = DN(1,0)*mu;
const double clhs49 = N[1]*clhs10;
const double clhs50 = DN(1,1) + clhs49;
const double clhs51 = clhs4*clhs50 + clhs48;
const double clhs52 = DN(1,0)*clhs18;
const double clhs53 = N[2]*clhs8;
const double clhs54 = DN(2,1)*clhs31;
const double clhs55 = DN(2,0)*clhs1;
const double clhs56 = clhs55*rho;
const double clhs57 = DN(2,1)*clhs2;
const double clhs58 = clhs57*rho;
const double clhs59 = clhs53 - clhs54 + clhs56 + clhs58;
const double clhs60 = clhs18*clhs59;
const double clhs61 = DN(0,0)*DN(2,0);
const double clhs62 = DN(2,1)*clhs11 + N[2]*clhs9;
const double clhs63 = DN(2,0)*clhs26 + clhs4*clhs61 + clhs62;
const double clhs64 = DN(2,0)*N[0];
const double clhs65 = DN(2,1)*clhs24;
const double clhs66 = DN(2,1)*N[0];
const double clhs67 = clhs43*clhs64 + clhs46*clhs66 + clhs65*mu;
const double clhs68 = DN(2,0)*mu;
const double clhs69 = N[2]*clhs10;
const double clhs70 = DN(2,1) + clhs69;
const double clhs71 = clhs4*clhs70 + clhs68;
const double clhs72 = DN(2,0)*clhs18;
const double clhs73 = N[3]*clhs8;
const double clhs74 = DN(3,1)*clhs31;
const double clhs75 = DN(3,0)*clhs1;
const double clhs76 = clhs75*rho;
const double clhs77 = DN(3,1)*clhs2;
const double clhs78 = clhs77*rho;
const double clhs79 = clhs73 - clhs74 + clhs76 + clhs78;
const double clhs80 = clhs18*clhs79;
const double clhs81 = DN(0,0)*DN(3,0);
const double clhs82 = DN(3,1)*clhs11 + N[3]*clhs9;
const double clhs83 = DN(3,0)*clhs26 + clhs4*clhs81 + clhs82;
const double clhs84 = DN(3,0)*N[0];
const double clhs85 = DN(3,1)*clhs24;
const double clhs86 = DN(3,1)*N[0];
const double clhs87 = clhs43*clhs84 + clhs46*clhs86 + clhs85*mu;
const double clhs88 = N[3]*clhs10;
const double clhs89 = DN(3,1) + clhs88;
const double clhs90 = DN(3,0)*mu + clhs4*clhs89;
const double clhs91 = DN(3,0)*clhs18;
const double clhs92 = mu/pow(y, 2);
const double clhs93 = N[0]*clhs92;
const double clhs94 = clhs18*(clhs16 + clhs93);
const double clhs95 = clhs18*(clhs21 + clhs7 - clhs93);
const double clhs96 = N[1]*clhs92;
const double clhs97 = clhs18*(clhs37 + clhs96);
const double clhs98 = N[1]*clhs93 + clhs40;
const double clhs99 = N[2]*clhs92;
const double clhs100 = clhs18*(clhs59 + clhs99);
const double clhs101 = N[2]*clhs93 + clhs62;
const double clhs102 = N[3]*clhs92;
const double clhs103 = clhs18*(clhs102 + clhs79);
const double clhs104 = N[3]*clhs93 + clhs82;
const double clhs105 = DN(0,1)*y;
const double clhs106 = DN(1,1)*clhs105 + clhs39;
const double clhs107 = DN(2,1)*clhs105 + clhs61;
const double clhs108 = DN(3,1)*clhs105 + clhs81;
const double clhs109 = rho*(N[1]*clhs6 + clhs33);
const double clhs110 = rho*(N[1]*clhs20 + clhs35);
const double clhs111 = N[1]*clhs13 + N[1]*clhs15 + clhs11*clhs49;
const double clhs112 = DN(0,0)*N[1];
const double clhs113 = pow(DN(1,0), 2);
const double clhs114 = pow(DN(1,1), 2);
const double clhs115 = pow(N[1], 2);
const double clhs116 = DN(1,1)*mu;
const double clhs117 = N[1]*clhs34 + N[1]*clhs36 + clhs114*mu + clhs115*clhs8 + clhs116*clhs49;
const double clhs118 = DN(1,0)*DN(2,0);
const double clhs119 = DN(2,1)*clhs116 + N[2]*clhs30;
const double clhs120 = DN(2,0)*clhs48 + clhs118*clhs4 + clhs119;
const double clhs121 = DN(2,0)*N[1];
const double clhs122 = DN(2,1)*clhs49;
const double clhs123 = DN(2,1)*N[1];
const double clhs124 = clhs121*clhs43 + clhs122*mu + clhs123*clhs46;
const double clhs125 = DN(1,0)*DN(3,0);
const double clhs126 = DN(3,1)*clhs116 + N[3]*clhs30;
const double clhs127 = DN(3,0)*clhs48 + clhs125*clhs4 + clhs126;
const double clhs128 = DN(3,0)*N[1];
const double clhs129 = DN(3,1)*clhs49;
const double clhs130 = DN(3,1)*N[1];
const double clhs131 = clhs128*clhs43 + clhs129*mu + clhs130*clhs46;
const double clhs132 = clhs18*(clhs109 + clhs110 - clhs96);
const double clhs133 = N[2]*clhs96 + clhs119;
const double clhs134 = N[3]*clhs96 + clhs126;
const double clhs135 = DN(1,1)*y;
const double clhs136 = DN(2,1)*clhs135 + clhs118;
const double clhs137 = DN(3,1)*clhs135 + clhs125;
const double clhs138 = rho*(N[2]*clhs6 + clhs55);
const double clhs139 = rho*(N[2]*clhs20 + clhs57);
const double clhs140 = N[2]*clhs13 + N[2]*clhs15 + clhs11*clhs69;
const double clhs141 = DN(0,0)*N[2];
const double clhs142 = N[2]*clhs34 + N[2]*clhs36 + clhs116*clhs69;
const double clhs143 = DN(1,0)*N[2];
const double clhs144 = pow(DN(2,0), 2);
const double clhs145 = pow(DN(2,1), 2);
const double clhs146 = pow(N[2], 2);
const double clhs147 = DN(2,1)*clhs69;
const double clhs148 = N[2]*clhs56 + N[2]*clhs58 + clhs145*mu + clhs146*clhs8 + clhs147*mu;
const double clhs149 = DN(2,0)*DN(3,0);
const double clhs150 = DN(2,1)*DN(3,1);
const double clhs151 = N[3]*clhs53 + clhs150*mu;
const double clhs152 = DN(3,0)*clhs68 + clhs149*clhs4 + clhs151;
const double clhs153 = DN(3,0)*N[2];
const double clhs154 = DN(3,1)*clhs69;
const double clhs155 = DN(3,1)*N[2];
const double clhs156 = clhs153*clhs43 + clhs154*mu + clhs155*clhs46;
const double clhs157 = clhs18*(clhs138 + clhs139 - clhs99);
const double clhs158 = N[3]*clhs99 + clhs151;
const double clhs159 = DN(2,1)*y;
const double clhs160 = clhs149 + clhs150*y;
const double clhs161 = rho*(N[3]*clhs6 + clhs75);
const double clhs162 = rho*(N[3]*clhs20 + clhs77);
const double clhs163 = N[3]*clhs13 + N[3]*clhs15 + clhs11*clhs88;
const double clhs164 = DN(0,0)*N[3];
const double clhs165 = N[3]*clhs34 + N[3]*clhs36 + clhs116*clhs88;
const double clhs166 = DN(1,0)*N[3];
const double clhs167 = DN(2,1)*clhs88;
const double clhs168 = N[3]*clhs56 + N[3]*clhs58 + clhs167*mu;
const double clhs169 = DN(2,0)*N[3];
const double clhs170 = pow(DN(3,0), 2);
const double clhs171 = pow(DN(3,1), 2);
const double clhs172 = pow(N[3], 2);
const double clhs173 = DN(3,1)*clhs88;
const double clhs174 = N[3]*clhs76 + N[3]*clhs78 + clhs171*mu + clhs172*clhs8 + clhs173*mu;
const double clhs175 = clhs18*(-clhs102 + clhs161 + clhs162);
const double clhs176 = DN(3,1)*y;
lhs(0,0)=clhs0*clhs4 + clhs0*mu + clhs19*clhs21 + clhs19*clhs7 + clhs25;
lhs(0,1)=DN(0,0)*clhs28;
lhs(0,2)=DN(0,0)*N[0] + DN(0,1)*N[0] + clhs21*clhs29 + clhs29*clhs7;
lhs(0,3)=clhs21*clhs38 + clhs38*clhs7 + clhs41 + clhs47;
lhs(0,4)=DN(0,0)*clhs51;
lhs(0,5)=clhs21*clhs52 + clhs42 + clhs45 + clhs52*clhs7;
lhs(0,6)=clhs21*clhs60 + clhs60*clhs7 + clhs63 + clhs67;
lhs(0,7)=DN(0,0)*clhs71;
lhs(0,8)=clhs21*clhs72 + clhs64 + clhs66 + clhs7*clhs72;
lhs(0,9)=clhs21*clhs80 + clhs7*clhs80 + clhs83 + clhs87;
lhs(0,10)=DN(0,0)*clhs90;
lhs(0,11)=clhs21*clhs91 + clhs7*clhs91 + clhs84 + clhs86;
lhs(1,0)=0;
lhs(1,1)=clhs21*clhs94 + clhs23*clhs92 + clhs25 + clhs7*clhs94 - clhs93*clhs94;
lhs(1,2)=DN(0,1)*clhs95;
lhs(1,3)=0;
lhs(1,4)=clhs21*clhs97 + clhs47 + clhs7*clhs97 - clhs93*clhs97 + clhs98;
lhs(1,5)=DN(1,1)*clhs95;
lhs(1,6)=0;
lhs(1,7)=clhs100*clhs21 + clhs100*clhs7 - clhs100*clhs93 + clhs101 + clhs67;
lhs(1,8)=DN(2,1)*clhs95;
lhs(1,9)=0;
lhs(1,10)=clhs103*clhs21 + clhs103*clhs7 - clhs103*clhs93 + clhs104 + clhs87;
lhs(1,11)=DN(3,1)*clhs95;
lhs(2,0)=DN(0,0)*(N[0] + clhs17*(-1.0*clhs12 + 1.0*clhs13 + 1.0*clhs15 + 1.0*clhs9));
lhs(2,1)=N[0]*clhs27 + clhs105*clhs94 - clhs24*clhs94;
lhs(2,2)=clhs18*(-DN(0,1)*clhs24 + clhs0 + clhs22*y);
lhs(2,3)=clhs29*clhs37 + clhs42;
lhs(2,4)=N[0]*clhs50 + clhs105*clhs97 - clhs24*clhs97;
lhs(2,5)=clhs18*(clhs106 - clhs44);
lhs(2,6)=clhs29*clhs59 + clhs64;
lhs(2,7)=N[0]*clhs70 + clhs100*clhs105 - clhs100*clhs24;
lhs(2,8)=clhs18*(clhs107 - clhs65);
lhs(2,9)=clhs29*clhs79 + clhs84;
lhs(2,10)=N[0]*clhs89 + clhs103*clhs105 - clhs103*clhs24;
lhs(2,11)=clhs18*(clhs108 - clhs85);
lhs(3,0)=clhs109*clhs19 + clhs110*clhs19 + clhs111 + clhs41;
lhs(3,1)=DN(1,0)*clhs28;
lhs(3,2)=DN(0,1)*N[1] + clhs109*clhs29 + clhs110*clhs29 + clhs112;
lhs(3,3)=clhs109*clhs38 + clhs110*clhs38 + clhs113*clhs4 + clhs113*mu + clhs117;
lhs(3,4)=DN(1,0)*clhs51;
lhs(3,5)=DN(1,0)*N[1] + DN(1,1)*N[1] + clhs109*clhs52 + clhs110*clhs52;
lhs(3,6)=clhs109*clhs60 + clhs110*clhs60 + clhs120 + clhs124;
lhs(3,7)=DN(1,0)*clhs71;
lhs(3,8)=clhs109*clhs72 + clhs110*clhs72 + clhs121 + clhs123;
lhs(3,9)=clhs109*clhs80 + clhs110*clhs80 + clhs127 + clhs131;
lhs(3,10)=DN(1,0)*clhs90;
lhs(3,11)=clhs109*clhs91 + clhs110*clhs91 + clhs128 + clhs130;
lhs(4,0)=0;
lhs(4,1)=clhs109*clhs94 + clhs110*clhs94 + clhs111 - clhs94*clhs96 + clhs98;
lhs(4,2)=DN(0,1)*clhs132;
lhs(4,3)=0;
lhs(4,4)=clhs109*clhs97 + clhs110*clhs97 + clhs115*clhs92 + clhs117 - clhs96*clhs97;
lhs(4,5)=DN(1,1)*clhs132;
lhs(4,6)=0;
lhs(4,7)=clhs100*clhs109 + clhs100*clhs110 - clhs100*clhs96 + clhs124 + clhs133;
lhs(4,8)=DN(2,1)*clhs132;
lhs(4,9)=0;
lhs(4,10)=clhs103*clhs109 + clhs103*clhs110 - clhs103*clhs96 + clhs131 + clhs134;
lhs(4,11)=DN(3,1)*clhs132;
lhs(5,0)=clhs112 + clhs16*clhs52;
lhs(5,1)=N[1]*clhs27 + clhs135*clhs94 - clhs49*clhs94;
lhs(5,2)=clhs18*(-DN(0,1)*clhs49 + clhs106);
lhs(5,3)=DN(1,0)*(N[1] + clhs17*(1.0*clhs30 - 1.0*clhs32 + 1.0*clhs34 + 1.0*clhs36));
lhs(5,4)=N[1]*clhs50 + clhs135*clhs97 - clhs49*clhs97;
lhs(5,5)=clhs18*(-DN(1,1)*clhs49 + clhs113 + clhs114*y);
lhs(5,6)=clhs121 + clhs52*clhs59;
lhs(5,7)=N[1]*clhs70 + clhs100*clhs135 - clhs100*clhs49;
lhs(5,8)=clhs18*(-clhs122 + clhs136);
lhs(5,9)=clhs128 + clhs52*clhs79;
lhs(5,10)=N[1]*clhs89 + clhs103*clhs135 - clhs103*clhs49;
lhs(5,11)=clhs18*(-clhs129 + clhs137);
lhs(6,0)=clhs138*clhs19 + clhs139*clhs19 + clhs140 + clhs63;
lhs(6,1)=DN(2,0)*clhs28;
lhs(6,2)=DN(0,1)*N[2] + clhs138*clhs29 + clhs139*clhs29 + clhs141;
lhs(6,3)=clhs120 + clhs138*clhs38 + clhs139*clhs38 + clhs142;
lhs(6,4)=DN(2,0)*clhs51;
lhs(6,5)=DN(1,1)*N[2] + clhs138*clhs52 + clhs139*clhs52 + clhs143;
lhs(6,6)=clhs138*clhs60 + clhs139*clhs60 + clhs144*clhs4 + clhs144*mu + clhs148;
lhs(6,7)=DN(2,0)*clhs71;
lhs(6,8)=DN(2,0)*N[2] + DN(2,1)*N[2] + clhs138*clhs72 + clhs139*clhs72;
lhs(6,9)=clhs138*clhs80 + clhs139*clhs80 + clhs152 + clhs156;
lhs(6,10)=DN(2,0)*clhs90;
lhs(6,11)=clhs138*clhs91 + clhs139*clhs91 + clhs153 + clhs155;
lhs(7,0)=0;
lhs(7,1)=clhs101 + clhs138*clhs94 + clhs139*clhs94 + clhs140 - clhs94*clhs99;
lhs(7,2)=DN(0,1)*clhs157;
lhs(7,3)=0;
lhs(7,4)=clhs133 + clhs138*clhs97 + clhs139*clhs97 + clhs142 - clhs97*clhs99;
lhs(7,5)=DN(1,1)*clhs157;
lhs(7,6)=0;
lhs(7,7)=clhs100*clhs138 + clhs100*clhs139 - clhs100*clhs99 + clhs146*clhs92 + clhs148;
lhs(7,8)=DN(2,1)*clhs157;
lhs(7,9)=0;
lhs(7,10)=clhs103*clhs138 + clhs103*clhs139 - clhs103*clhs99 + clhs156 + clhs158;
lhs(7,11)=DN(3,1)*clhs157;
lhs(8,0)=clhs141 + clhs16*clhs72;
lhs(8,1)=N[2]*clhs27 + clhs159*clhs94 - clhs69*clhs94;
lhs(8,2)=clhs18*(-DN(0,1)*clhs69 + clhs107);
lhs(8,3)=clhs143 + clhs37*clhs72;
lhs(8,4)=N[2]*clhs50 + clhs159*clhs97 - clhs69*clhs97;
lhs(8,5)=clhs18*(-DN(1,1)*clhs69 + clhs136);
lhs(8,6)=DN(2,0)*(N[2] + clhs17*(1.0*clhs53 - 1.0*clhs54 + 1.0*clhs56 + 1.0*clhs58));
lhs(8,7)=N[2]*clhs70 + clhs100*clhs159 - clhs100*clhs69;
lhs(8,8)=clhs18*(clhs144 + clhs145*y - clhs147);
lhs(8,9)=clhs153 + clhs72*clhs79;
lhs(8,10)=N[2]*clhs89 + clhs103*clhs159 - clhs103*clhs69;
lhs(8,11)=clhs18*(-clhs154 + clhs160);
lhs(9,0)=clhs161*clhs19 + clhs162*clhs19 + clhs163 + clhs83;
lhs(9,1)=DN(3,0)*clhs28;
lhs(9,2)=DN(0,1)*N[3] + clhs161*clhs29 + clhs162*clhs29 + clhs164;
lhs(9,3)=clhs127 + clhs161*clhs38 + clhs162*clhs38 + clhs165;
lhs(9,4)=DN(3,0)*clhs51;
lhs(9,5)=DN(1,1)*N[3] + clhs161*clhs52 + clhs162*clhs52 + clhs166;
lhs(9,6)=clhs152 + clhs161*clhs60 + clhs162*clhs60 + clhs168;
lhs(9,7)=DN(3,0)*clhs71;
lhs(9,8)=DN(2,1)*N[3] + clhs161*clhs72 + clhs162*clhs72 + clhs169;
lhs(9,9)=clhs161*clhs80 + clhs162*clhs80 + clhs170*clhs4 + clhs170*mu + clhs174;
lhs(9,10)=DN(3,0)*clhs90;
lhs(9,11)=DN(3,0)*N[3] + DN(3,1)*N[3] + clhs161*clhs91 + clhs162*clhs91;
lhs(10,0)=0;
lhs(10,1)=-clhs102*clhs94 + clhs104 + clhs161*clhs94 + clhs162*clhs94 + clhs163;
lhs(10,2)=DN(0,1)*clhs175;
lhs(10,3)=0;
lhs(10,4)=-clhs102*clhs97 + clhs134 + clhs161*clhs97 + clhs162*clhs97 + clhs165;
lhs(10,5)=DN(1,1)*clhs175;
lhs(10,6)=0;
lhs(10,7)=-clhs100*clhs102 + clhs100*clhs161 + clhs100*clhs162 + clhs158 + clhs168;
lhs(10,8)=DN(2,1)*clhs175;
lhs(10,9)=0;
lhs(10,10)=-clhs102*clhs103 + clhs103*clhs161 + clhs103*clhs162 + clhs172*clhs92 + clhs174;
lhs(10,11)=DN(3,1)*clhs175;
lhs(11,0)=clhs16*clhs91 + clhs164;
lhs(11,1)=N[3]*clhs27 + clhs176*clhs94 - clhs88*clhs94;
lhs(11,2)=clhs18*(-DN(0,1)*clhs88 + clhs108);
lhs(11,3)=clhs166 + clhs37*clhs91;
lhs(11,4)=N[3]*clhs50 + clhs176*clhs97 - clhs88*clhs97;
lhs(11,5)=clhs18*(-DN(1,1)*clhs88 + clhs137);
lhs(11,6)=clhs169 + clhs59*clhs91;
lhs(11,7)=N[3]*clhs70 + clhs100*clhs176 - clhs100*clhs88;
lhs(11,8)=clhs18*(clhs160 - clhs167);
lhs(11,9)=DN(3,0)*(N[3] + clhs17*(1.0*clhs73 - 1.0*clhs74 + 1.0*clhs76 + 1.0*clhs78));
lhs(11,10)=N[3]*clhs89 + clhs103*clhs176 - clhs103*clhs88;
lhs(11,11)=clhs18*(clhs170 + clhs171*y - clhs173);


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
const double crhs4 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs5 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double crhs6 = crhs5*mu;
const double crhs7 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs8 = 1.0/y;
const double crhs9 = crhs6*crhs8;
const double crhs10 = N[0]*crhs2;
const double crhs11 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double crhs12 = crhs11*rho;
const double crhs13 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double crhs14 = crhs13*rho;
const double crhs15 = crhs14*crhs5;
const double crhs16 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
const double crhs17 = rho*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs13, 2));
const double crhs18 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs19 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs20 = crhs18*crhs8 + crhs19;
const double crhs21 = 2*(crhs2 + crhs20)*(crhs17*h/stab_c1 + mu);
const double crhs22 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double crhs23 = DN(0,0)*crhs11 + N[0]*crhs22;
const double crhs24 = 1.0/(crhs17/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs25 = 1.0*crhs24;
const double crhs26 = crhs25*(crhs0 + crhs12*crhs2 + crhs15 + crhs16 - crhs7*rho - crhs9);
const double crhs27 = crhs26*rho;
const double crhs28 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double crhs29 = DN(0,1)*crhs13 + N[0]*crhs28;
const double crhs30 = crhs19*mu;
const double crhs31 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs32 = pow(y, -2);
const double crhs33 = crhs18*crhs32*mu;
const double crhs34 = crhs30*crhs8;
const double crhs35 = crhs12*crhs4;
const double crhs36 = crhs14*crhs19;
const double crhs37 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
const double crhs38 = crhs1 - crhs31*rho + crhs33 - crhs34 + crhs35 + crhs36 + crhs37;
const double crhs39 = crhs25*crhs38;
const double crhs40 = crhs39*rho;
const double crhs41 = crhs39*y;
const double crhs42 = DN(1,0)*mu;
const double crhs43 = N[1]*crhs2;
const double crhs44 = DN(1,0)*crhs11 + N[1]*crhs22;
const double crhs45 = DN(1,1)*crhs13 + N[1]*crhs28;
const double crhs46 = DN(2,0)*mu;
const double crhs47 = N[2]*crhs2;
const double crhs48 = DN(2,0)*crhs11 + N[2]*crhs22;
const double crhs49 = DN(2,1)*crhs13 + N[2]*crhs28;
rhs[0]=-DN(0,0)*crhs21 - DN(0,1)*crhs6 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs15 - N[0]*crhs16 + N[0]*crhs7*rho - N[0]*crhs9 - crhs10*crhs12 - crhs2*crhs3 - crhs23*crhs27 - crhs27*crhs29 - crhs3*crhs4;
rhs[1]=-DN(0,1)*crhs30 + 1.0*N[0]*crhs24*crhs32*crhs38*mu + N[0]*crhs31*rho - N[0]*crhs33 - N[0]*crhs34 - N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs37 - crhs23*crhs40 - crhs29*crhs40;
rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs41 - N[0]*crhs20 + 1.0*N[0]*crhs24*crhs38*crhs8 - crhs10;
rhs[3]=-DN(1,0)*crhs21 - DN(1,1)*crhs6 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs15 - N[1]*crhs16 + N[1]*crhs7*rho - N[1]*crhs9 - crhs12*crhs43 - crhs2*crhs42 - crhs27*crhs44 - crhs27*crhs45 - crhs4*crhs42;
rhs[4]=-DN(1,1)*crhs30 + 1.0*N[1]*crhs24*crhs32*crhs38*mu + N[1]*crhs31*rho - N[1]*crhs33 - N[1]*crhs34 - N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs37 - crhs40*crhs44 - crhs40*crhs45;
rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs41 - N[1]*crhs20 + 1.0*N[1]*crhs24*crhs38*crhs8 - crhs43;
rhs[6]=-DN(2,0)*crhs21 - DN(2,1)*crhs6 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs15 - N[2]*crhs16 + N[2]*crhs7*rho - N[2]*crhs9 - crhs12*crhs47 - crhs2*crhs46 - crhs27*crhs48 - crhs27*crhs49 - crhs4*crhs46;
rhs[7]=-DN(2,1)*crhs30 + 1.0*N[2]*crhs24*crhs32*crhs38*mu + N[2]*crhs31*rho - N[2]*crhs33 - N[2]*crhs34 - N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs37 - crhs40*crhs48 - crhs40*crhs49;
rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs41 - N[2]*crhs20 + 1.0*N[2]*crhs24*crhs38*crhs8 - crhs47;


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
const double crhs4 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs5 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double crhs6 = crhs5*mu;
const double crhs7 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs8 = 1.0/y;
const double crhs9 = crhs6*crhs8;
const double crhs10 = N[0]*crhs2;
const double crhs11 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double crhs12 = crhs11*rho;
const double crhs13 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double crhs14 = crhs13*rho;
const double crhs15 = crhs14*crhs5;
const double crhs16 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
const double crhs17 = rho*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs13, 2));
const double crhs18 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs19 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs20 = crhs18*crhs8 + crhs19;
const double crhs21 = 2*(crhs2 + crhs20)*(crhs17*h/stab_c1 + mu);
const double crhs22 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double crhs23 = DN(0,0)*crhs11 + N[0]*crhs22;
const double crhs24 = 1.0/(crhs17/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs25 = 1.0*crhs24;
const double crhs26 = crhs25*(crhs0 + crhs12*crhs2 + crhs15 + crhs16 - crhs7*rho - crhs9);
const double crhs27 = crhs26*rho;
const double crhs28 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double crhs29 = DN(0,1)*crhs13 + N[0]*crhs28;
const double crhs30 = crhs19*mu;
const double crhs31 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs32 = pow(y, -2);
const double crhs33 = crhs18*crhs32*mu;
const double crhs34 = crhs30*crhs8;
const double crhs35 = crhs12*crhs4;
const double crhs36 = crhs14*crhs19;
const double crhs37 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
const double crhs38 = crhs1 - crhs31*rho + crhs33 - crhs34 + crhs35 + crhs36 + crhs37;
const double crhs39 = crhs25*crhs38;
const double crhs40 = crhs39*rho;
const double crhs41 = crhs39*y;
const double crhs42 = DN(1,0)*mu;
const double crhs43 = N[1]*crhs2;
const double crhs44 = DN(1,0)*crhs11 + N[1]*crhs22;
const double crhs45 = DN(1,1)*crhs13 + N[1]*crhs28;
const double crhs46 = DN(2,0)*mu;
const double crhs47 = N[2]*crhs2;
const double crhs48 = DN(2,0)*crhs11 + N[2]*crhs22;
const double crhs49 = DN(2,1)*crhs13 + N[2]*crhs28;
const double crhs50 = DN(3,0)*mu;
const double crhs51 = N[3]*crhs2;
const double crhs52 = DN(3,0)*crhs11 + N[3]*crhs22;
const double crhs53 = DN(3,1)*crhs13 + N[3]*crhs28;
rhs[0]=-DN(0,0)*crhs21 - DN(0,1)*crhs6 - N[0]*crhs0 - N[0]*crhs1 - N[0]*crhs15 - N[0]*crhs16 + N[0]*crhs7*rho - N[0]*crhs9 - crhs10*crhs12 - crhs2*crhs3 - crhs23*crhs27 - crhs27*crhs29 - crhs3*crhs4;
rhs[1]=-DN(0,1)*crhs30 + 1.0*N[0]*crhs24*crhs32*crhs38*mu + N[0]*crhs31*rho - N[0]*crhs33 - N[0]*crhs34 - N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs37 - crhs23*crhs40 - crhs29*crhs40;
rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs41 - N[0]*crhs20 + 1.0*N[0]*crhs24*crhs38*crhs8 - crhs10;
rhs[3]=-DN(1,0)*crhs21 - DN(1,1)*crhs6 - N[1]*crhs0 - N[1]*crhs1 - N[1]*crhs15 - N[1]*crhs16 + N[1]*crhs7*rho - N[1]*crhs9 - crhs12*crhs43 - crhs2*crhs42 - crhs27*crhs44 - crhs27*crhs45 - crhs4*crhs42;
rhs[4]=-DN(1,1)*crhs30 + 1.0*N[1]*crhs24*crhs32*crhs38*mu + N[1]*crhs31*rho - N[1]*crhs33 - N[1]*crhs34 - N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs37 - crhs40*crhs44 - crhs40*crhs45;
rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs41 - N[1]*crhs20 + 1.0*N[1]*crhs24*crhs38*crhs8 - crhs43;
rhs[6]=-DN(2,0)*crhs21 - DN(2,1)*crhs6 - N[2]*crhs0 - N[2]*crhs1 - N[2]*crhs15 - N[2]*crhs16 + N[2]*crhs7*rho - N[2]*crhs9 - crhs12*crhs47 - crhs2*crhs46 - crhs27*crhs48 - crhs27*crhs49 - crhs4*crhs46;
rhs[7]=-DN(2,1)*crhs30 + 1.0*N[2]*crhs24*crhs32*crhs38*mu + N[2]*crhs31*rho - N[2]*crhs33 - N[2]*crhs34 - N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs37 - crhs40*crhs48 - crhs40*crhs49;
rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs41 - N[2]*crhs20 + 1.0*N[2]*crhs24*crhs38*crhs8 - crhs47;
rhs[9]=-DN(3,0)*crhs21 - DN(3,1)*crhs6 - N[3]*crhs0 - N[3]*crhs1 - N[3]*crhs15 - N[3]*crhs16 + N[3]*crhs7*rho - N[3]*crhs9 - crhs12*crhs51 - crhs2*crhs50 - crhs27*crhs52 - crhs27*crhs53 - crhs4*crhs50;
rhs[10]=-DN(3,1)*crhs30 + 1.0*N[3]*crhs24*crhs32*crhs38*mu + N[3]*crhs31*rho - N[3]*crhs33 - N[3]*crhs34 - N[3]*crhs35 - N[3]*crhs36 - N[3]*crhs37 - crhs40*crhs52 - crhs40*crhs53;
rhs[11]=-DN(3,0)*crhs26 - DN(3,1)*crhs41 - N[3]*crhs20 + 1.0*N[3]*crhs24*crhs38*crhs8 - crhs51;


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