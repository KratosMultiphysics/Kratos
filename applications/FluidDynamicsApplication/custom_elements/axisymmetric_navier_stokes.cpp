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

template <class TElementData>
void AxisymmetricNavierStokes<TElementData>::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Set the constitutive law pointer to null as the axisymmetric element hardcodes a Newtonian fluid viscous behavior
    // Note that to use a constitutive law the gradient in cylindrical coordinates would require the corresponding stress
    // implementation in cylindrical coordinates within the mechanical response calculation
    this->GetConstitutiveLaw() = nullptr;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int AxisymmetricNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }

    // Check variables used by TElementData
    out = TElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    // Check nodes
    const double check_tolerance = 1.0e-8;
    const auto& r_geometry = this->GetGeometry();
    for (const auto& r_node : r_geometry) {
        // Check nodal DOFs
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE,r_node);
        // Axisymmetry: nodes are in XY plane
        KRATOS_ERROR_IF(std::abs(r_node.Z()) > check_tolerance ) << "Node " << r_node.Id() << "has non-zero Z coordinate." << std::endl;
        // Axisymmetry: check that there are no negative y-coordinates (radius is always positive)
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

    // rData.lhs *= 2.0 * Globals::Pi * y * rData.Weight;
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
const double clhs5 = bdf0*rho;
const double clhs6 = N[0]*clhs5;
const double clhs7 = 1.0/y;
const double clhs8 = clhs7*mu;
const double clhs9 = DN(0,1)*clhs8;
const double clhs10 = DN(0,0)*clhs1;
const double clhs11 = clhs10*rho;
const double clhs12 = DN(0,1)*clhs2;
const double clhs13 = clhs12*rho;
const double clhs14 = clhs11 + clhs13 + clhs6 - clhs9;
const double clhs15 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double clhs16 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs17 = 1.0*clhs16;
const double clhs18 = clhs17*rho;
const double clhs19 = clhs18*(N[0]*clhs15 + clhs10);
const double clhs20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double clhs21 = clhs18*(N[0]*clhs20 + clhs12);
const double clhs22 = pow(DN(0,1), 2);
const double clhs23 = pow(N[0], 2);
const double clhs24 = N[0]*clhs11 + N[0]*clhs13 + clhs0*mu + clhs22*mu + clhs23*clhs5;
const double clhs25 = N[0]*clhs7;
const double clhs26 = DN(0,1) + clhs25;
const double clhs27 = DN(0,0)*clhs4;
const double clhs28 = clhs26*clhs27;
const double clhs29 = N[1]*clhs5;
const double clhs30 = DN(1,1)*clhs8;
const double clhs31 = DN(1,0)*clhs1;
const double clhs32 = clhs31*rho;
const double clhs33 = DN(1,1)*clhs2;
const double clhs34 = clhs33*rho;
const double clhs35 = clhs29 - clhs30 + clhs32 + clhs34;
const double clhs36 = DN(0,0)*DN(1,0);
const double clhs37 = DN(0,1)*DN(1,1);
const double clhs38 = N[1]*clhs6 + clhs36*mu + clhs37*mu;
const double clhs39 = clhs36*clhs4 + clhs38;
const double clhs40 = DN(1,0)*N[0];
const double clhs41 = clhs1*rho;
const double clhs42 = DN(1,1)*N[0];
const double clhs43 = clhs2*rho;
const double clhs44 = clhs40*clhs41 + clhs42*clhs43;
const double clhs45 = N[1]*clhs7;
const double clhs46 = DN(1,1) + clhs45;
const double clhs47 = clhs27*clhs46;
const double clhs48 = DN(0,0)*N[1];
const double clhs49 = N[2]*clhs5;
const double clhs50 = DN(2,1)*clhs8;
const double clhs51 = DN(2,0)*clhs1;
const double clhs52 = clhs51*rho;
const double clhs53 = DN(2,1)*clhs2;
const double clhs54 = clhs53*rho;
const double clhs55 = clhs49 - clhs50 + clhs52 + clhs54;
const double clhs56 = DN(0,0)*DN(2,0);
const double clhs57 = DN(0,1)*DN(2,1);
const double clhs58 = N[2]*clhs6 + clhs56*mu + clhs57*mu;
const double clhs59 = clhs4*clhs56 + clhs58;
const double clhs60 = DN(2,0)*N[0];
const double clhs61 = DN(2,1)*N[0];
const double clhs62 = clhs41*clhs60 + clhs43*clhs61;
const double clhs63 = N[2]*clhs7;
const double clhs64 = DN(2,1) + clhs63;
const double clhs65 = clhs27*clhs64;
const double clhs66 = DN(0,0)*N[2];
const double clhs67 = clhs26*clhs4;
const double clhs68 = mu/pow(y, 2);
const double clhs69 = N[0]*clhs68 + clhs14;
const double clhs70 = DN(1,0)*clhs67;
const double clhs71 = DN(0,1)*clhs4;
const double clhs72 = clhs4*clhs46;
const double clhs73 = N[1]*clhs68 + clhs35;
const double clhs74 = DN(2,0)*clhs67;
const double clhs75 = clhs4*clhs64;
const double clhs76 = N[2]*clhs68 + clhs55;
const double clhs77 = DN(0,1)*clhs17;
const double clhs78 = clhs17*clhs25;
const double clhs79 = DN(0,0)*clhs17;
const double clhs80 = N[1]*clhs25;
const double clhs81 = clhs36 + clhs37;
const double clhs82 = N[2]*clhs25;
const double clhs83 = clhs56 + clhs57;
const double clhs84 = clhs18*(N[1]*clhs15 + clhs31);
const double clhs85 = clhs18*(N[1]*clhs20 + clhs33);
const double clhs86 = N[1]*clhs11 + N[1]*clhs13;
const double clhs87 = pow(DN(1,0), 2);
const double clhs88 = pow(DN(1,1), 2);
const double clhs89 = pow(N[1], 2);
const double clhs90 = N[1]*clhs32 + N[1]*clhs34 + clhs5*clhs89 + clhs87*mu + clhs88*mu;
const double clhs91 = DN(1,0)*clhs4;
const double clhs92 = clhs46*clhs91;
const double clhs93 = DN(1,0)*DN(2,0);
const double clhs94 = DN(1,1)*DN(2,1);
const double clhs95 = N[2]*clhs29 + clhs93*mu + clhs94*mu;
const double clhs96 = clhs4*clhs93 + clhs95;
const double clhs97 = DN(2,0)*N[1];
const double clhs98 = DN(2,1)*N[1];
const double clhs99 = clhs41*clhs97 + clhs43*clhs98;
const double clhs100 = clhs64*clhs91;
const double clhs101 = DN(1,0)*N[2];
const double clhs102 = DN(2,0)*clhs72;
const double clhs103 = DN(1,0)*clhs17;
const double clhs104 = DN(1,1)*clhs17;
const double clhs105 = clhs17*clhs45;
const double clhs106 = N[2]*clhs45;
const double clhs107 = clhs93 + clhs94;
const double clhs108 = clhs18*(N[2]*clhs15 + clhs51);
const double clhs109 = clhs18*(N[2]*clhs20 + clhs53);
const double clhs110 = N[2]*clhs11 + N[2]*clhs13;
const double clhs111 = N[2]*clhs32 + N[2]*clhs34;
const double clhs112 = pow(DN(2,0), 2);
const double clhs113 = pow(DN(2,1), 2);
const double clhs114 = pow(N[2], 2);
const double clhs115 = N[2]*clhs52 + N[2]*clhs54 + clhs112*mu + clhs113*mu + clhs114*clhs5;
const double clhs116 = DN(2,0)*clhs75;
const double clhs117 = DN(2,0)*clhs17;
const double clhs118 = DN(2,1)*clhs17;
const double clhs119 = clhs17*clhs63;
lhs(0,0)=clhs0*clhs4 + clhs14*clhs19 + clhs14*clhs21 + clhs24;
lhs(0,1)=clhs28;
lhs(0,2)=DN(0,0)*(-N[0] + clhs19 + clhs21);
lhs(0,3)=clhs19*clhs35 + clhs21*clhs35 + clhs39 + clhs44;
lhs(0,4)=clhs47;
lhs(0,5)=DN(1,0)*clhs19 + DN(1,0)*clhs21 - clhs48;
lhs(0,6)=clhs19*clhs55 + clhs21*clhs55 + clhs59 + clhs62;
lhs(0,7)=clhs65;
lhs(0,8)=DN(2,0)*clhs19 + DN(2,0)*clhs21 - clhs66;
lhs(1,0)=clhs28;
lhs(1,1)=DN(0,1)*clhs67 + clhs19*clhs69 + clhs21*clhs69 + clhs24 + clhs25*clhs67;
lhs(1,2)=DN(0,1)*clhs19 + DN(0,1)*clhs21 - N[0]*clhs26;
lhs(1,3)=clhs70;
lhs(1,4)=clhs19*clhs73 + clhs21*clhs73 + clhs25*clhs72 + clhs38 + clhs44 + clhs46*clhs71;
lhs(1,5)=DN(1,1)*clhs19 + DN(1,1)*clhs21 - N[1]*clhs26;
lhs(1,6)=clhs74;
lhs(1,7)=clhs19*clhs76 + clhs21*clhs76 + clhs25*clhs75 + clhs58 + clhs62 + clhs64*clhs71;
lhs(1,8)=DN(2,1)*clhs19 + DN(2,1)*clhs21 - N[2]*clhs26;
lhs(2,0)=DN(0,0)*(N[0] + clhs16*(1.0*clhs11 + 1.0*clhs13 + 1.0*clhs6 - 1.0*clhs9));
lhs(2,1)=DN(0,1)*N[0] + clhs23*clhs7 + clhs69*clhs77 - clhs69*clhs78;
lhs(2,2)=clhs17*(-DN(0,1)*clhs25 + clhs0 + clhs22);
lhs(2,3)=clhs35*clhs79 + clhs40;
lhs(2,4)=clhs42 + clhs73*clhs77 - clhs73*clhs78 + clhs80;
lhs(2,5)=clhs17*(-DN(1,1)*clhs25 + clhs81);
lhs(2,6)=clhs55*clhs79 + clhs60;
lhs(2,7)=clhs61 + clhs76*clhs77 - clhs76*clhs78 + clhs82;
lhs(2,8)=clhs17*(-DN(2,1)*clhs25 + clhs83);
lhs(3,0)=clhs14*clhs84 + clhs14*clhs85 + clhs39 + clhs86;
lhs(3,1)=clhs70;
lhs(3,2)=DN(0,0)*clhs84 + DN(0,0)*clhs85 - clhs40;
lhs(3,3)=clhs35*clhs84 + clhs35*clhs85 + clhs4*clhs87 + clhs90;
lhs(3,4)=clhs92;
lhs(3,5)=DN(1,0)*(-N[1] + clhs84 + clhs85);
lhs(3,6)=clhs55*clhs84 + clhs55*clhs85 + clhs96 + clhs99;
lhs(3,7)=clhs100;
lhs(3,8)=DN(2,0)*clhs84 + DN(2,0)*clhs85 - clhs101;
lhs(4,0)=clhs47;
lhs(4,1)=DN(1,1)*clhs67 + clhs38 + clhs45*clhs67 + clhs69*clhs84 + clhs69*clhs85 + clhs86;
lhs(4,2)=DN(0,1)*clhs84 + DN(0,1)*clhs85 - N[0]*clhs46;
lhs(4,3)=clhs92;
lhs(4,4)=DN(1,1)*clhs72 + clhs45*clhs72 + clhs73*clhs84 + clhs73*clhs85 + clhs90;
lhs(4,5)=DN(1,1)*clhs84 + DN(1,1)*clhs85 - N[1]*clhs46;
lhs(4,6)=clhs102;
lhs(4,7)=DN(1,1)*clhs75 + clhs45*clhs75 + clhs76*clhs84 + clhs76*clhs85 + clhs95 + clhs99;
lhs(4,8)=DN(2,1)*clhs84 + DN(2,1)*clhs85 - N[2]*clhs46;
lhs(5,0)=clhs103*clhs14 + clhs48;
lhs(5,1)=DN(0,1)*N[1] + clhs104*clhs69 - clhs105*clhs69 + clhs80;
lhs(5,2)=clhs17*(-DN(0,1)*clhs45 + clhs81);
lhs(5,3)=DN(1,0)*(N[1] + clhs16*(1.0*clhs29 - 1.0*clhs30 + 1.0*clhs32 + 1.0*clhs34));
lhs(5,4)=DN(1,1)*N[1] + clhs104*clhs73 - clhs105*clhs73 + clhs7*clhs89;
lhs(5,5)=clhs17*(-DN(1,1)*clhs45 + clhs87 + clhs88);
lhs(5,6)=clhs103*clhs55 + clhs97;
lhs(5,7)=clhs104*clhs76 - clhs105*clhs76 + clhs106 + clhs98;
lhs(5,8)=clhs17*(-DN(2,1)*clhs45 + clhs107);
lhs(6,0)=clhs108*clhs14 + clhs109*clhs14 + clhs110 + clhs59;
lhs(6,1)=clhs74;
lhs(6,2)=DN(0,0)*clhs108 + DN(0,0)*clhs109 - clhs60;
lhs(6,3)=clhs108*clhs35 + clhs109*clhs35 + clhs111 + clhs96;
lhs(6,4)=clhs102;
lhs(6,5)=DN(1,0)*clhs108 + DN(1,0)*clhs109 - clhs97;
lhs(6,6)=clhs108*clhs55 + clhs109*clhs55 + clhs112*clhs4 + clhs115;
lhs(6,7)=clhs116;
lhs(6,8)=DN(2,0)*(-N[2] + clhs108 + clhs109);
lhs(7,0)=clhs65;
lhs(7,1)=DN(2,1)*clhs67 + clhs108*clhs69 + clhs109*clhs69 + clhs110 + clhs58 + clhs63*clhs67;
lhs(7,2)=DN(0,1)*clhs108 + DN(0,1)*clhs109 - N[0]*clhs64;
lhs(7,3)=clhs100;
lhs(7,4)=DN(2,1)*clhs72 + clhs108*clhs73 + clhs109*clhs73 + clhs111 + clhs63*clhs72 + clhs95;
lhs(7,5)=DN(1,1)*clhs108 + DN(1,1)*clhs109 - N[1]*clhs64;
lhs(7,6)=clhs116;
lhs(7,7)=DN(2,1)*clhs75 + clhs108*clhs76 + clhs109*clhs76 + clhs115 + clhs63*clhs75;
lhs(7,8)=DN(2,1)*clhs108 + DN(2,1)*clhs109 - N[2]*clhs64;
lhs(8,0)=clhs117*clhs14 + clhs66;
lhs(8,1)=DN(0,1)*N[2] + clhs118*clhs69 - clhs119*clhs69 + clhs82;
lhs(8,2)=clhs17*(-DN(0,1)*clhs63 + clhs83);
lhs(8,3)=clhs101 + clhs117*clhs35;
lhs(8,4)=DN(1,1)*N[2] + clhs106 + clhs118*clhs73 - clhs119*clhs73;
lhs(8,5)=clhs17*(-DN(1,1)*clhs63 + clhs107);
lhs(8,6)=DN(2,0)*(N[2] + clhs16*(1.0*clhs49 - 1.0*clhs50 + 1.0*clhs52 + 1.0*clhs54));
lhs(8,7)=DN(2,1)*N[2] + clhs114*clhs7 + clhs118*clhs76 - clhs119*clhs76;
lhs(8,8)=clhs17*(-DN(2,1)*clhs63 + clhs112 + clhs113);


    // Add intermediate results to local system
    noalias(rLHS) += 2.0 * Globals::Pi * y * rData.Weight * lhs;
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
const double clhs5 = bdf0*rho;
const double clhs6 = N[0]*clhs5;
const double clhs7 = 1.0/y;
const double clhs8 = clhs7*mu;
const double clhs9 = DN(0,1)*clhs8;
const double clhs10 = DN(0,0)*clhs1;
const double clhs11 = clhs10*rho;
const double clhs12 = DN(0,1)*clhs2;
const double clhs13 = clhs12*rho;
const double clhs14 = clhs11 + clhs13 + clhs6 - clhs9;
const double clhs15 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double clhs16 = 1.0/(clhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs17 = 1.0*clhs16;
const double clhs18 = clhs17*rho;
const double clhs19 = clhs18*(N[0]*clhs15 + clhs10);
const double clhs20 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double clhs21 = clhs18*(N[0]*clhs20 + clhs12);
const double clhs22 = pow(DN(0,1), 2);
const double clhs23 = pow(N[0], 2);
const double clhs24 = N[0]*clhs11 + N[0]*clhs13 + clhs0*mu + clhs22*mu + clhs23*clhs5;
const double clhs25 = N[0]*clhs7;
const double clhs26 = DN(0,1) + clhs25;
const double clhs27 = DN(0,0)*clhs4;
const double clhs28 = clhs26*clhs27;
const double clhs29 = N[1]*clhs5;
const double clhs30 = DN(1,1)*clhs8;
const double clhs31 = DN(1,0)*clhs1;
const double clhs32 = clhs31*rho;
const double clhs33 = DN(1,1)*clhs2;
const double clhs34 = clhs33*rho;
const double clhs35 = clhs29 - clhs30 + clhs32 + clhs34;
const double clhs36 = DN(0,0)*DN(1,0);
const double clhs37 = DN(0,1)*DN(1,1);
const double clhs38 = N[1]*clhs6 + clhs36*mu + clhs37*mu;
const double clhs39 = clhs36*clhs4 + clhs38;
const double clhs40 = DN(1,0)*N[0];
const double clhs41 = clhs1*rho;
const double clhs42 = DN(1,1)*N[0];
const double clhs43 = clhs2*rho;
const double clhs44 = clhs40*clhs41 + clhs42*clhs43;
const double clhs45 = N[1]*clhs7;
const double clhs46 = DN(1,1) + clhs45;
const double clhs47 = clhs27*clhs46;
const double clhs48 = DN(0,0)*N[1];
const double clhs49 = N[2]*clhs5;
const double clhs50 = DN(2,1)*clhs8;
const double clhs51 = DN(2,0)*clhs1;
const double clhs52 = clhs51*rho;
const double clhs53 = DN(2,1)*clhs2;
const double clhs54 = clhs53*rho;
const double clhs55 = clhs49 - clhs50 + clhs52 + clhs54;
const double clhs56 = DN(0,0)*DN(2,0);
const double clhs57 = DN(0,1)*DN(2,1);
const double clhs58 = N[2]*clhs6 + clhs56*mu + clhs57*mu;
const double clhs59 = clhs4*clhs56 + clhs58;
const double clhs60 = DN(2,0)*N[0];
const double clhs61 = DN(2,1)*N[0];
const double clhs62 = clhs41*clhs60 + clhs43*clhs61;
const double clhs63 = N[2]*clhs7;
const double clhs64 = DN(2,1) + clhs63;
const double clhs65 = clhs27*clhs64;
const double clhs66 = DN(0,0)*N[2];
const double clhs67 = N[3]*clhs5;
const double clhs68 = DN(3,1)*clhs8;
const double clhs69 = DN(3,0)*clhs1;
const double clhs70 = clhs69*rho;
const double clhs71 = DN(3,1)*clhs2;
const double clhs72 = clhs71*rho;
const double clhs73 = clhs67 - clhs68 + clhs70 + clhs72;
const double clhs74 = DN(0,0)*DN(3,0);
const double clhs75 = DN(0,1)*DN(3,1);
const double clhs76 = N[3]*clhs6 + clhs74*mu + clhs75*mu;
const double clhs77 = clhs4*clhs74 + clhs76;
const double clhs78 = DN(3,0)*N[0];
const double clhs79 = DN(3,1)*N[0];
const double clhs80 = clhs41*clhs78 + clhs43*clhs79;
const double clhs81 = N[3]*clhs7;
const double clhs82 = DN(3,1) + clhs81;
const double clhs83 = clhs27*clhs82;
const double clhs84 = DN(0,0)*N[3];
const double clhs85 = clhs26*clhs4;
const double clhs86 = mu/pow(y, 2);
const double clhs87 = N[0]*clhs86 + clhs14;
const double clhs88 = DN(1,0)*clhs85;
const double clhs89 = DN(0,1)*clhs4;
const double clhs90 = clhs4*clhs46;
const double clhs91 = N[1]*clhs86 + clhs35;
const double clhs92 = DN(2,0)*clhs85;
const double clhs93 = clhs4*clhs64;
const double clhs94 = N[2]*clhs86 + clhs55;
const double clhs95 = DN(3,0)*clhs85;
const double clhs96 = clhs4*clhs82;
const double clhs97 = N[3]*clhs86 + clhs73;
const double clhs98 = DN(0,1)*clhs17;
const double clhs99 = clhs17*clhs25;
const double clhs100 = DN(0,0)*clhs17;
const double clhs101 = N[1]*clhs25;
const double clhs102 = clhs36 + clhs37;
const double clhs103 = N[2]*clhs25;
const double clhs104 = clhs56 + clhs57;
const double clhs105 = N[3]*clhs25;
const double clhs106 = clhs74 + clhs75;
const double clhs107 = clhs18*(N[1]*clhs15 + clhs31);
const double clhs108 = clhs18*(N[1]*clhs20 + clhs33);
const double clhs109 = N[1]*clhs11 + N[1]*clhs13;
const double clhs110 = pow(DN(1,0), 2);
const double clhs111 = pow(DN(1,1), 2);
const double clhs112 = pow(N[1], 2);
const double clhs113 = N[1]*clhs32 + N[1]*clhs34 + clhs110*mu + clhs111*mu + clhs112*clhs5;
const double clhs114 = DN(1,0)*clhs4;
const double clhs115 = clhs114*clhs46;
const double clhs116 = DN(1,0)*DN(2,0);
const double clhs117 = DN(1,1)*DN(2,1);
const double clhs118 = N[2]*clhs29 + clhs116*mu + clhs117*mu;
const double clhs119 = clhs116*clhs4 + clhs118;
const double clhs120 = DN(2,0)*N[1];
const double clhs121 = DN(2,1)*N[1];
const double clhs122 = clhs120*clhs41 + clhs121*clhs43;
const double clhs123 = clhs114*clhs64;
const double clhs124 = DN(1,0)*N[2];
const double clhs125 = DN(1,0)*DN(3,0);
const double clhs126 = DN(1,1)*DN(3,1);
const double clhs127 = N[3]*clhs29 + clhs125*mu + clhs126*mu;
const double clhs128 = clhs125*clhs4 + clhs127;
const double clhs129 = DN(3,0)*N[1];
const double clhs130 = DN(3,1)*N[1];
const double clhs131 = clhs129*clhs41 + clhs130*clhs43;
const double clhs132 = clhs114*clhs82;
const double clhs133 = DN(1,0)*N[3];
const double clhs134 = DN(2,0)*clhs90;
const double clhs135 = DN(1,1)*clhs4;
const double clhs136 = DN(3,0)*clhs90;
const double clhs137 = DN(1,0)*clhs17;
const double clhs138 = DN(1,1)*clhs17;
const double clhs139 = clhs17*clhs45;
const double clhs140 = N[2]*clhs45;
const double clhs141 = clhs116 + clhs117;
const double clhs142 = N[3]*clhs45;
const double clhs143 = clhs125 + clhs126;
const double clhs144 = clhs18*(N[2]*clhs15 + clhs51);
const double clhs145 = clhs18*(N[2]*clhs20 + clhs53);
const double clhs146 = N[2]*clhs11 + N[2]*clhs13;
const double clhs147 = N[2]*clhs32 + N[2]*clhs34;
const double clhs148 = pow(DN(2,0), 2);
const double clhs149 = pow(DN(2,1), 2);
const double clhs150 = pow(N[2], 2);
const double clhs151 = N[2]*clhs52 + N[2]*clhs54 + clhs148*mu + clhs149*mu + clhs150*clhs5;
const double clhs152 = DN(2,0)*clhs4;
const double clhs153 = clhs152*clhs64;
const double clhs154 = DN(2,0)*DN(3,0);
const double clhs155 = DN(2,1)*DN(3,1);
const double clhs156 = N[3]*clhs49 + clhs154*mu + clhs155*mu;
const double clhs157 = clhs154*clhs4 + clhs156;
const double clhs158 = DN(3,0)*N[2];
const double clhs159 = DN(3,1)*N[2];
const double clhs160 = clhs158*clhs41 + clhs159*clhs43;
const double clhs161 = clhs152*clhs82;
const double clhs162 = DN(2,0)*N[3];
const double clhs163 = DN(3,0)*clhs93;
const double clhs164 = DN(2,0)*clhs17;
const double clhs165 = DN(2,1)*clhs17;
const double clhs166 = clhs17*clhs63;
const double clhs167 = N[3]*clhs63;
const double clhs168 = clhs154 + clhs155;
const double clhs169 = clhs18*(N[3]*clhs15 + clhs69);
const double clhs170 = clhs18*(N[3]*clhs20 + clhs71);
const double clhs171 = N[3]*clhs11 + N[3]*clhs13;
const double clhs172 = N[3]*clhs32 + N[3]*clhs34;
const double clhs173 = N[3]*clhs52 + N[3]*clhs54;
const double clhs174 = pow(DN(3,0), 2);
const double clhs175 = pow(DN(3,1), 2);
const double clhs176 = pow(N[3], 2);
const double clhs177 = N[3]*clhs70 + N[3]*clhs72 + clhs174*mu + clhs175*mu + clhs176*clhs5;
const double clhs178 = DN(3,0)*clhs96;
const double clhs179 = DN(3,0)*clhs17;
const double clhs180 = DN(3,1)*clhs17;
const double clhs181 = clhs17*clhs81;
lhs(0,0)=clhs0*clhs4 + clhs14*clhs19 + clhs14*clhs21 + clhs24;
lhs(0,1)=clhs28;
lhs(0,2)=DN(0,0)*(-N[0] + clhs19 + clhs21);
lhs(0,3)=clhs19*clhs35 + clhs21*clhs35 + clhs39 + clhs44;
lhs(0,4)=clhs47;
lhs(0,5)=DN(1,0)*clhs19 + DN(1,0)*clhs21 - clhs48;
lhs(0,6)=clhs19*clhs55 + clhs21*clhs55 + clhs59 + clhs62;
lhs(0,7)=clhs65;
lhs(0,8)=DN(2,0)*clhs19 + DN(2,0)*clhs21 - clhs66;
lhs(0,9)=clhs19*clhs73 + clhs21*clhs73 + clhs77 + clhs80;
lhs(0,10)=clhs83;
lhs(0,11)=DN(3,0)*clhs19 + DN(3,0)*clhs21 - clhs84;
lhs(1,0)=clhs28;
lhs(1,1)=DN(0,1)*clhs85 + clhs19*clhs87 + clhs21*clhs87 + clhs24 + clhs25*clhs85;
lhs(1,2)=DN(0,1)*clhs19 + DN(0,1)*clhs21 - N[0]*clhs26;
lhs(1,3)=clhs88;
lhs(1,4)=clhs19*clhs91 + clhs21*clhs91 + clhs25*clhs90 + clhs38 + clhs44 + clhs46*clhs89;
lhs(1,5)=DN(1,1)*clhs19 + DN(1,1)*clhs21 - N[1]*clhs26;
lhs(1,6)=clhs92;
lhs(1,7)=clhs19*clhs94 + clhs21*clhs94 + clhs25*clhs93 + clhs58 + clhs62 + clhs64*clhs89;
lhs(1,8)=DN(2,1)*clhs19 + DN(2,1)*clhs21 - N[2]*clhs26;
lhs(1,9)=clhs95;
lhs(1,10)=clhs19*clhs97 + clhs21*clhs97 + clhs25*clhs96 + clhs76 + clhs80 + clhs82*clhs89;
lhs(1,11)=DN(3,1)*clhs19 + DN(3,1)*clhs21 - N[3]*clhs26;
lhs(2,0)=DN(0,0)*(N[0] + clhs16*(1.0*clhs11 + 1.0*clhs13 + 1.0*clhs6 - 1.0*clhs9));
lhs(2,1)=DN(0,1)*N[0] + clhs23*clhs7 + clhs87*clhs98 - clhs87*clhs99;
lhs(2,2)=clhs17*(-DN(0,1)*clhs25 + clhs0 + clhs22);
lhs(2,3)=clhs100*clhs35 + clhs40;
lhs(2,4)=clhs101 + clhs42 + clhs91*clhs98 - clhs91*clhs99;
lhs(2,5)=clhs17*(-DN(1,1)*clhs25 + clhs102);
lhs(2,6)=clhs100*clhs55 + clhs60;
lhs(2,7)=clhs103 + clhs61 + clhs94*clhs98 - clhs94*clhs99;
lhs(2,8)=clhs17*(-DN(2,1)*clhs25 + clhs104);
lhs(2,9)=clhs100*clhs73 + clhs78;
lhs(2,10)=clhs105 + clhs79 + clhs97*clhs98 - clhs97*clhs99;
lhs(2,11)=clhs17*(-DN(3,1)*clhs25 + clhs106);
lhs(3,0)=clhs107*clhs14 + clhs108*clhs14 + clhs109 + clhs39;
lhs(3,1)=clhs88;
lhs(3,2)=DN(0,0)*clhs107 + DN(0,0)*clhs108 - clhs40;
lhs(3,3)=clhs107*clhs35 + clhs108*clhs35 + clhs110*clhs4 + clhs113;
lhs(3,4)=clhs115;
lhs(3,5)=DN(1,0)*(-N[1] + clhs107 + clhs108);
lhs(3,6)=clhs107*clhs55 + clhs108*clhs55 + clhs119 + clhs122;
lhs(3,7)=clhs123;
lhs(3,8)=DN(2,0)*clhs107 + DN(2,0)*clhs108 - clhs124;
lhs(3,9)=clhs107*clhs73 + clhs108*clhs73 + clhs128 + clhs131;
lhs(3,10)=clhs132;
lhs(3,11)=DN(3,0)*clhs107 + DN(3,0)*clhs108 - clhs133;
lhs(4,0)=clhs47;
lhs(4,1)=DN(1,1)*clhs85 + clhs107*clhs87 + clhs108*clhs87 + clhs109 + clhs38 + clhs45*clhs85;
lhs(4,2)=DN(0,1)*clhs107 + DN(0,1)*clhs108 - N[0]*clhs46;
lhs(4,3)=clhs115;
lhs(4,4)=DN(1,1)*clhs90 + clhs107*clhs91 + clhs108*clhs91 + clhs113 + clhs45*clhs90;
lhs(4,5)=DN(1,1)*clhs107 + DN(1,1)*clhs108 - N[1]*clhs46;
lhs(4,6)=clhs134;
lhs(4,7)=clhs107*clhs94 + clhs108*clhs94 + clhs118 + clhs122 + clhs135*clhs64 + clhs45*clhs93;
lhs(4,8)=DN(2,1)*clhs107 + DN(2,1)*clhs108 - N[2]*clhs46;
lhs(4,9)=clhs136;
lhs(4,10)=clhs107*clhs97 + clhs108*clhs97 + clhs127 + clhs131 + clhs135*clhs82 + clhs45*clhs96;
lhs(4,11)=DN(3,1)*clhs107 + DN(3,1)*clhs108 - N[3]*clhs46;
lhs(5,0)=clhs137*clhs14 + clhs48;
lhs(5,1)=DN(0,1)*N[1] + clhs101 + clhs138*clhs87 - clhs139*clhs87;
lhs(5,2)=clhs17*(-DN(0,1)*clhs45 + clhs102);
lhs(5,3)=DN(1,0)*(N[1] + clhs16*(1.0*clhs29 - 1.0*clhs30 + 1.0*clhs32 + 1.0*clhs34));
lhs(5,4)=DN(1,1)*N[1] + clhs112*clhs7 + clhs138*clhs91 - clhs139*clhs91;
lhs(5,5)=clhs17*(-DN(1,1)*clhs45 + clhs110 + clhs111);
lhs(5,6)=clhs120 + clhs137*clhs55;
lhs(5,7)=clhs121 + clhs138*clhs94 - clhs139*clhs94 + clhs140;
lhs(5,8)=clhs17*(-DN(2,1)*clhs45 + clhs141);
lhs(5,9)=clhs129 + clhs137*clhs73;
lhs(5,10)=clhs130 + clhs138*clhs97 - clhs139*clhs97 + clhs142;
lhs(5,11)=clhs17*(-DN(3,1)*clhs45 + clhs143);
lhs(6,0)=clhs14*clhs144 + clhs14*clhs145 + clhs146 + clhs59;
lhs(6,1)=clhs92;
lhs(6,2)=DN(0,0)*clhs144 + DN(0,0)*clhs145 - clhs60;
lhs(6,3)=clhs119 + clhs144*clhs35 + clhs145*clhs35 + clhs147;
lhs(6,4)=clhs134;
lhs(6,5)=DN(1,0)*clhs144 + DN(1,0)*clhs145 - clhs120;
lhs(6,6)=clhs144*clhs55 + clhs145*clhs55 + clhs148*clhs4 + clhs151;
lhs(6,7)=clhs153;
lhs(6,8)=DN(2,0)*(-N[2] + clhs144 + clhs145);
lhs(6,9)=clhs144*clhs73 + clhs145*clhs73 + clhs157 + clhs160;
lhs(6,10)=clhs161;
lhs(6,11)=DN(3,0)*clhs144 + DN(3,0)*clhs145 - clhs162;
lhs(7,0)=clhs65;
lhs(7,1)=DN(2,1)*clhs85 + clhs144*clhs87 + clhs145*clhs87 + clhs146 + clhs58 + clhs63*clhs85;
lhs(7,2)=DN(0,1)*clhs144 + DN(0,1)*clhs145 - N[0]*clhs64;
lhs(7,3)=clhs123;
lhs(7,4)=DN(2,1)*clhs90 + clhs118 + clhs144*clhs91 + clhs145*clhs91 + clhs147 + clhs63*clhs90;
lhs(7,5)=DN(1,1)*clhs144 + DN(1,1)*clhs145 - N[1]*clhs64;
lhs(7,6)=clhs153;
lhs(7,7)=DN(2,1)*clhs93 + clhs144*clhs94 + clhs145*clhs94 + clhs151 + clhs63*clhs93;
lhs(7,8)=DN(2,1)*clhs144 + DN(2,1)*clhs145 - N[2]*clhs64;
lhs(7,9)=clhs163;
lhs(7,10)=DN(2,1)*clhs96 + clhs144*clhs97 + clhs145*clhs97 + clhs156 + clhs160 + clhs63*clhs96;
lhs(7,11)=DN(3,1)*clhs144 + DN(3,1)*clhs145 - N[3]*clhs64;
lhs(8,0)=clhs14*clhs164 + clhs66;
lhs(8,1)=DN(0,1)*N[2] + clhs103 + clhs165*clhs87 - clhs166*clhs87;
lhs(8,2)=clhs17*(-DN(0,1)*clhs63 + clhs104);
lhs(8,3)=clhs124 + clhs164*clhs35;
lhs(8,4)=DN(1,1)*N[2] + clhs140 + clhs165*clhs91 - clhs166*clhs91;
lhs(8,5)=clhs17*(-DN(1,1)*clhs63 + clhs141);
lhs(8,6)=DN(2,0)*(N[2] + clhs16*(1.0*clhs49 - 1.0*clhs50 + 1.0*clhs52 + 1.0*clhs54));
lhs(8,7)=DN(2,1)*N[2] + clhs150*clhs7 + clhs165*clhs94 - clhs166*clhs94;
lhs(8,8)=clhs17*(-DN(2,1)*clhs63 + clhs148 + clhs149);
lhs(8,9)=clhs158 + clhs164*clhs73;
lhs(8,10)=clhs159 + clhs165*clhs97 - clhs166*clhs97 + clhs167;
lhs(8,11)=clhs17*(-DN(3,1)*clhs63 + clhs168);
lhs(9,0)=clhs14*clhs169 + clhs14*clhs170 + clhs171 + clhs77;
lhs(9,1)=clhs95;
lhs(9,2)=DN(0,0)*clhs169 + DN(0,0)*clhs170 - clhs78;
lhs(9,3)=clhs128 + clhs169*clhs35 + clhs170*clhs35 + clhs172;
lhs(9,4)=clhs136;
lhs(9,5)=DN(1,0)*clhs169 + DN(1,0)*clhs170 - clhs129;
lhs(9,6)=clhs157 + clhs169*clhs55 + clhs170*clhs55 + clhs173;
lhs(9,7)=clhs163;
lhs(9,8)=DN(2,0)*clhs169 + DN(2,0)*clhs170 - clhs158;
lhs(9,9)=clhs169*clhs73 + clhs170*clhs73 + clhs174*clhs4 + clhs177;
lhs(9,10)=clhs178;
lhs(9,11)=DN(3,0)*(-N[3] + clhs169 + clhs170);
lhs(10,0)=clhs83;
lhs(10,1)=DN(3,1)*clhs85 + clhs169*clhs87 + clhs170*clhs87 + clhs171 + clhs76 + clhs81*clhs85;
lhs(10,2)=DN(0,1)*clhs169 + DN(0,1)*clhs170 - N[0]*clhs82;
lhs(10,3)=clhs132;
lhs(10,4)=DN(3,1)*clhs90 + clhs127 + clhs169*clhs91 + clhs170*clhs91 + clhs172 + clhs81*clhs90;
lhs(10,5)=DN(1,1)*clhs169 + DN(1,1)*clhs170 - N[1]*clhs82;
lhs(10,6)=clhs161;
lhs(10,7)=DN(3,1)*clhs93 + clhs156 + clhs169*clhs94 + clhs170*clhs94 + clhs173 + clhs81*clhs93;
lhs(10,8)=DN(2,1)*clhs169 + DN(2,1)*clhs170 - N[2]*clhs82;
lhs(10,9)=clhs178;
lhs(10,10)=DN(3,1)*clhs96 + clhs169*clhs97 + clhs170*clhs97 + clhs177 + clhs81*clhs96;
lhs(10,11)=DN(3,1)*clhs169 + DN(3,1)*clhs170 - N[3]*clhs82;
lhs(11,0)=clhs14*clhs179 + clhs84;
lhs(11,1)=DN(0,1)*N[3] + clhs105 + clhs180*clhs87 - clhs181*clhs87;
lhs(11,2)=clhs17*(-DN(0,1)*clhs81 + clhs106);
lhs(11,3)=clhs133 + clhs179*clhs35;
lhs(11,4)=DN(1,1)*N[3] + clhs142 + clhs180*clhs91 - clhs181*clhs91;
lhs(11,5)=clhs17*(-DN(1,1)*clhs81 + clhs143);
lhs(11,6)=clhs162 + clhs179*clhs55;
lhs(11,7)=DN(2,1)*N[3] + clhs167 + clhs180*clhs94 - clhs181*clhs94;
lhs(11,8)=clhs17*(-DN(2,1)*clhs81 + clhs168);
lhs(11,9)=DN(3,0)*(N[3] + clhs16*(1.0*clhs67 - 1.0*clhs68 + 1.0*clhs70 + 1.0*clhs72));
lhs(11,10)=DN(3,1)*N[3] + clhs176*clhs7 + clhs180*clhs97 - clhs181*clhs97;
lhs(11,11)=clhs17*(-DN(3,1)*clhs81 + clhs174 + clhs175);


    // Add intermediate results to local system
    noalias(rLHS) += 2.0 * Globals::Pi * y * rData.Weight * lhs;
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
const double crhs6 = N[0]*crhs1;
const double crhs7 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0);
const double crhs8 = crhs7*rho;
const double crhs9 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1);
const double crhs10 = crhs9*rho;
const double crhs11 = crhs10*crhs3;
const double crhs12 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)));
const double crhs13 = rho*stab_c2*sqrt(pow(crhs7, 2) + pow(crhs9, 2));
const double crhs14 = 1.0/y;
const double crhs15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs16 = crhs14*crhs15;
const double crhs17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs18 = (crhs13*h/stab_c1 + mu)*(crhs1 + crhs16 + crhs17);
const double crhs19 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0);
const double crhs20 = DN(0,0)*crhs7 + N[0]*crhs19;
const double crhs21 = 1.0/(crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = 1.0*crhs21;
const double crhs23 = crhs22*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs1*crhs8 + crhs11 + crhs12 - crhs14*crhs4 - crhs5*rho);
const double crhs24 = crhs23*rho;
const double crhs25 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1);
const double crhs26 = DN(0,1)*crhs9 + N[0]*crhs25;
const double crhs27 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs28 = crhs17*mu;
const double crhs29 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs30 = N[0]*crhs14;
const double crhs31 = crhs27*crhs8;
const double crhs32 = N[0]*crhs17;
const double crhs33 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)));
const double crhs34 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs10*crhs17 - crhs14*crhs28 + crhs15*mu/pow(y, 2) - crhs29*rho + crhs31 + crhs33;
const double crhs35 = crhs22*crhs34;
const double crhs36 = crhs35*rho;
const double crhs37 = DN(1,0)*mu;
const double crhs38 = N[1]*crhs1;
const double crhs39 = DN(1,0)*crhs7 + N[1]*crhs19;
const double crhs40 = DN(1,1)*crhs9 + N[1]*crhs25;
const double crhs41 = N[1]*crhs14;
const double crhs42 = N[1]*crhs17;
const double crhs43 = DN(2,0)*mu;
const double crhs44 = N[2]*crhs1;
const double crhs45 = DN(2,0)*crhs7 + N[2]*crhs19;
const double crhs46 = DN(2,1)*crhs9 + N[2]*crhs25;
const double crhs47 = N[2]*crhs14;
const double crhs48 = N[2]*crhs17;
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,1)*crhs4 - N[0]*crhs11 - N[0]*crhs12 + N[0]*crhs5*rho - crhs1*crhs2 - crhs20*crhs24 - crhs24*crhs26 - crhs6*crhs8;
rhs[1]=-DN(0,1)*crhs18 - DN(0,1)*crhs28 + N[0]*crhs29*rho - N[0]*crhs31 - N[0]*crhs33 + crhs0*(DN(0,1) + crhs30) - crhs10*crhs32 - crhs18*crhs30 - crhs2*crhs27 - crhs20*crhs36 - crhs26*crhs36;
rhs[2]=-DN(0,0)*crhs23 - DN(0,1)*crhs35 + 1.0*N[0]*crhs14*crhs21*crhs34 - N[0]*crhs16 - crhs32 - crhs6;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,1)*crhs4 - N[1]*crhs11 - N[1]*crhs12 + N[1]*crhs5*rho - crhs1*crhs37 - crhs24*crhs39 - crhs24*crhs40 - crhs38*crhs8;
rhs[4]=-DN(1,1)*crhs18 - DN(1,1)*crhs28 + N[1]*crhs29*rho - N[1]*crhs31 - N[1]*crhs33 + crhs0*(DN(1,1) + crhs41) - crhs10*crhs42 - crhs18*crhs41 - crhs27*crhs37 - crhs36*crhs39 - crhs36*crhs40;
rhs[5]=-DN(1,0)*crhs23 - DN(1,1)*crhs35 + 1.0*N[1]*crhs14*crhs21*crhs34 - N[1]*crhs16 - crhs38 - crhs42;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,1)*crhs4 - N[2]*crhs11 - N[2]*crhs12 + N[2]*crhs5*rho - crhs1*crhs43 - crhs24*crhs45 - crhs24*crhs46 - crhs44*crhs8;
rhs[7]=-DN(2,1)*crhs18 - DN(2,1)*crhs28 + N[2]*crhs29*rho - N[2]*crhs31 - N[2]*crhs33 + crhs0*(DN(2,1) + crhs47) - crhs10*crhs48 - crhs18*crhs47 - crhs27*crhs43 - crhs36*crhs45 - crhs36*crhs46;
rhs[8]=-DN(2,0)*crhs23 - DN(2,1)*crhs35 + 1.0*N[2]*crhs14*crhs21*crhs34 - N[2]*crhs16 - crhs44 - crhs48;


    noalias(rRHS) += 2.0 * Globals::Pi * y * rData.Weight * rhs;
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
const double crhs6 = N[0]*crhs1;
const double crhs7 = N[0]*v_conv(0,0) + N[1]*v_conv(1,0) + N[2]*v_conv(2,0) + N[3]*v_conv(3,0);
const double crhs8 = crhs7*rho;
const double crhs9 = N[0]*v_conv(0,1) + N[1]*v_conv(1,1) + N[2]*v_conv(2,1) + N[3]*v_conv(3,1);
const double crhs10 = crhs9*rho;
const double crhs11 = crhs10*crhs3;
const double crhs12 = rho*(N[0]*(bdf0*v(0,0) + bdf1*v_n(0,0) + bdf2*v_nn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*v_n(1,0) + bdf2*v_nn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*v_n(2,0) + bdf2*v_nn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*v_n(3,0) + bdf2*v_nn(3,0)));
const double crhs13 = rho*stab_c2*sqrt(pow(crhs7, 2) + pow(crhs9, 2));
const double crhs14 = 1.0/y;
const double crhs15 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs16 = crhs14*crhs15;
const double crhs17 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs18 = (crhs13*h/stab_c1 + mu)*(crhs1 + crhs16 + crhs17);
const double crhs19 = DN(0,0)*v_conv(0,0) + DN(1,0)*v_conv(1,0) + DN(2,0)*v_conv(2,0) + DN(3,0)*v_conv(3,0);
const double crhs20 = DN(0,0)*crhs7 + N[0]*crhs19;
const double crhs21 = 1.0/(crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 = 1.0*crhs21;
const double crhs23 = crhs22*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs1*crhs8 + crhs11 + crhs12 - crhs14*crhs4 - crhs5*rho);
const double crhs24 = crhs23*rho;
const double crhs25 = DN(0,1)*v_conv(0,1) + DN(1,1)*v_conv(1,1) + DN(2,1)*v_conv(2,1) + DN(3,1)*v_conv(3,1);
const double crhs26 = DN(0,1)*crhs9 + N[0]*crhs25;
const double crhs27 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs28 = crhs17*mu;
const double crhs29 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs30 = N[0]*crhs14;
const double crhs31 = crhs27*crhs8;
const double crhs32 = N[0]*crhs17;
const double crhs33 = rho*(N[0]*(bdf0*v(0,1) + bdf1*v_n(0,1) + bdf2*v_nn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*v_n(1,1) + bdf2*v_nn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*v_n(2,1) + bdf2*v_nn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*v_n(3,1) + bdf2*v_nn(3,1)));
const double crhs34 = DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs10*crhs17 - crhs14*crhs28 + crhs15*mu/pow(y, 2) - crhs29*rho + crhs31 + crhs33;
const double crhs35 = crhs22*crhs34;
const double crhs36 = crhs35*rho;
const double crhs37 = DN(1,0)*mu;
const double crhs38 = N[1]*crhs1;
const double crhs39 = DN(1,0)*crhs7 + N[1]*crhs19;
const double crhs40 = DN(1,1)*crhs9 + N[1]*crhs25;
const double crhs41 = N[1]*crhs14;
const double crhs42 = N[1]*crhs17;
const double crhs43 = DN(2,0)*mu;
const double crhs44 = N[2]*crhs1;
const double crhs45 = DN(2,0)*crhs7 + N[2]*crhs19;
const double crhs46 = DN(2,1)*crhs9 + N[2]*crhs25;
const double crhs47 = N[2]*crhs14;
const double crhs48 = N[2]*crhs17;
const double crhs49 = DN(3,0)*mu;
const double crhs50 = N[3]*crhs1;
const double crhs51 = DN(3,0)*crhs7 + N[3]*crhs19;
const double crhs52 = DN(3,1)*crhs9 + N[3]*crhs25;
const double crhs53 = N[3]*crhs14;
const double crhs54 = N[3]*crhs17;
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,1)*crhs4 - N[0]*crhs11 - N[0]*crhs12 + N[0]*crhs5*rho - crhs1*crhs2 - crhs20*crhs24 - crhs24*crhs26 - crhs6*crhs8;
rhs[1]=-DN(0,1)*crhs18 - DN(0,1)*crhs28 + N[0]*crhs29*rho - N[0]*crhs31 - N[0]*crhs33 + crhs0*(DN(0,1) + crhs30) - crhs10*crhs32 - crhs18*crhs30 - crhs2*crhs27 - crhs20*crhs36 - crhs26*crhs36;
rhs[2]=-DN(0,0)*crhs23 - DN(0,1)*crhs35 + 1.0*N[0]*crhs14*crhs21*crhs34 - N[0]*crhs16 - crhs32 - crhs6;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,1)*crhs4 - N[1]*crhs11 - N[1]*crhs12 + N[1]*crhs5*rho - crhs1*crhs37 - crhs24*crhs39 - crhs24*crhs40 - crhs38*crhs8;
rhs[4]=-DN(1,1)*crhs18 - DN(1,1)*crhs28 + N[1]*crhs29*rho - N[1]*crhs31 - N[1]*crhs33 + crhs0*(DN(1,1) + crhs41) - crhs10*crhs42 - crhs18*crhs41 - crhs27*crhs37 - crhs36*crhs39 - crhs36*crhs40;
rhs[5]=-DN(1,0)*crhs23 - DN(1,1)*crhs35 + 1.0*N[1]*crhs14*crhs21*crhs34 - N[1]*crhs16 - crhs38 - crhs42;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,1)*crhs4 - N[2]*crhs11 - N[2]*crhs12 + N[2]*crhs5*rho - crhs1*crhs43 - crhs24*crhs45 - crhs24*crhs46 - crhs44*crhs8;
rhs[7]=-DN(2,1)*crhs18 - DN(2,1)*crhs28 + N[2]*crhs29*rho - N[2]*crhs31 - N[2]*crhs33 + crhs0*(DN(2,1) + crhs47) - crhs10*crhs48 - crhs18*crhs47 - crhs27*crhs43 - crhs36*crhs45 - crhs36*crhs46;
rhs[8]=-DN(2,0)*crhs23 - DN(2,1)*crhs35 + 1.0*N[2]*crhs14*crhs21*crhs34 - N[2]*crhs16 - crhs44 - crhs48;
rhs[9]=DN(3,0)*crhs0 - DN(3,0)*crhs18 - DN(3,1)*crhs4 - N[3]*crhs11 - N[3]*crhs12 + N[3]*crhs5*rho - crhs1*crhs49 - crhs24*crhs51 - crhs24*crhs52 - crhs50*crhs8;
rhs[10]=-DN(3,1)*crhs18 - DN(3,1)*crhs28 + N[3]*crhs29*rho - N[3]*crhs31 - N[3]*crhs33 + crhs0*(DN(3,1) + crhs53) - crhs10*crhs54 - crhs18*crhs53 - crhs27*crhs49 - crhs36*crhs51 - crhs36*crhs52;
rhs[11]=-DN(3,0)*crhs23 - DN(3,1)*crhs35 + 1.0*N[3]*crhs14*crhs21*crhs34 - N[3]*crhs16 - crhs50 - crhs54;


    noalias(rRHS) += 2.0 * Globals::Pi * y * rData.Weight * rhs;
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