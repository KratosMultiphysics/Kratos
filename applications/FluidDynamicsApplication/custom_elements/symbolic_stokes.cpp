//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// External includes

// System includes

// Project includes

// Application includes
#include "symbolic_stokes.h"
#include "custom_utilities/symbolic_stokes_data.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry)
{}

template <class TElementData>
SymbolicStokes<TElementData>::SymbolicStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties)
{}

template <class TElementData>
SymbolicStokes<TElementData>::~SymbolicStokes()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer SymbolicStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokes>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int SymbolicStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
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
std::string SymbolicStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "SymbolicStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void SymbolicStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
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
void SymbolicStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void SymbolicStokes<TElementData>::AddBoundaryTraction(
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
void SymbolicStokes< SymbolicStokesData<2,3> >::ComputeGaussPointLHSContribution(
    SymbolicStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const array_1d<double, 2> v_aux = prod(trans(rData.Velocity), N);

    auto& lhs = rData.lhs;

    const double clhs0 =             bdf0*rho;
const double clhs1 =             pow(N[0], 2)*clhs0;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs3 =             C(0,2)*DN(0,0);
const double clhs4 =             C(2,2)*DN(0,1) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             rho*stab_c2*sqrt(pow(v_aux[0], 2) + pow(v_aux[1], 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             C(0,1)*DN(0,1) + clhs3;
const double clhs9 =             C(1,2)*DN(0,1);
const double clhs10 =             C(2,2)*DN(0,0) + clhs9;
const double clhs11 =             DN(0,0)*clhs7;
const double clhs12 =             DN(0,1)*clhs11;
const double clhs13 =             DN(0,0)*N[0];
const double clhs14 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs15 =             C(0,2)*DN(1,0);
const double clhs16 =             C(2,2)*DN(1,1) + clhs15;
const double clhs17 =             DN(0,0)*DN(1,0);
const double clhs18 =             N[0]*clhs0;
const double clhs19 =             N[1]*clhs18;
const double clhs20 =             clhs17*clhs7 + clhs19;
const double clhs21 =             C(0,1)*DN(1,1) + clhs15;
const double clhs22 =             C(1,2)*DN(1,1);
const double clhs23 =             C(2,2)*DN(1,0) + clhs22;
const double clhs24 =             DN(1,1)*clhs11;
const double clhs25 =             DN(0,0)*N[1];
const double clhs26 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs27 =             C(0,2)*DN(2,0);
const double clhs28 =             C(2,2)*DN(2,1) + clhs27;
const double clhs29 =             DN(0,0)*DN(2,0);
const double clhs30 =             N[2]*clhs18;
const double clhs31 =             clhs29*clhs7 + clhs30;
const double clhs32 =             C(0,1)*DN(2,1) + clhs27;
const double clhs33 =             C(1,2)*DN(2,1);
const double clhs34 =             C(2,2)*DN(2,0) + clhs33;
const double clhs35 =             DN(2,1)*clhs11;
const double clhs36 =             DN(0,0)*N[2];
const double clhs37 =             C(0,1)*DN(0,0) + clhs9;
const double clhs38 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs39 =             pow(DN(0,1), 2);
const double clhs40 =             DN(0,1)*N[0];
const double clhs41 =             C(0,1)*DN(1,0) + clhs22;
const double clhs42 =             DN(0,1)*clhs7;
const double clhs43 =             DN(1,0)*clhs42;
const double clhs44 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs45 =             DN(0,1)*DN(1,1);
const double clhs46 =             clhs19 + clhs45*clhs7;
const double clhs47 =             DN(0,1)*N[1];
const double clhs48 =             C(0,1)*DN(2,0) + clhs33;
const double clhs49 =             DN(2,0)*clhs42;
const double clhs50 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs51 =             DN(0,1)*DN(2,1);
const double clhs52 =             clhs30 + clhs51*clhs7;
const double clhs53 =             DN(0,1)*N[2];
const double clhs54 =             1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs55 =             clhs0*clhs54;
const double clhs56 =             clhs55 + 1;
const double clhs57 =             DN(1,0)*N[0];
const double clhs58 =             DN(1,1)*N[0];
const double clhs59 =             clhs54*(clhs17 + clhs45);
const double clhs60 =             DN(2,0)*N[0];
const double clhs61 =             DN(2,1)*N[0];
const double clhs62 =             clhs54*(clhs29 + clhs51);
const double clhs63 =             pow(N[1], 2)*clhs0;
const double clhs64 =             pow(DN(1,0), 2);
const double clhs65 =             DN(1,0)*clhs7;
const double clhs66 =             DN(1,1)*clhs65;
const double clhs67 =             DN(1,0)*N[1];
const double clhs68 =             DN(1,0)*DN(2,0);
const double clhs69 =             N[1]*N[2]*clhs0;
const double clhs70 =             clhs68*clhs7 + clhs69;
const double clhs71 =             DN(2,1)*clhs65;
const double clhs72 =             DN(1,0)*N[2];
const double clhs73 =             pow(DN(1,1), 2);
const double clhs74 =             DN(1,1)*N[1];
const double clhs75 =             DN(2,0)*clhs7;
const double clhs76 =             DN(1,1)*clhs75;
const double clhs77 =             DN(1,1)*DN(2,1);
const double clhs78 =             clhs69 + clhs7*clhs77;
const double clhs79 =             DN(1,1)*N[2];
const double clhs80 =             DN(2,0)*N[1];
const double clhs81 =             DN(2,1)*N[1];
const double clhs82 =             clhs54*(clhs68 + clhs77);
const double clhs83 =             pow(N[2], 2)*clhs0;
const double clhs84 =             pow(DN(2,0), 2);
const double clhs85 =             DN(2,1)*clhs75;
const double clhs86 =             DN(2,0)*N[2];
const double clhs87 =             pow(DN(2,1), 2);
const double clhs88 =             DN(2,1)*N[2];
            lhs(0,0)=DN(0,0)*clhs2 + DN(0,1)*clhs4 + clhs1 + clhs5*clhs7;
            lhs(0,1)=DN(0,0)*clhs8 + DN(0,1)*clhs10 + clhs12;
            lhs(0,2)=-clhs13;
            lhs(0,3)=DN(0,0)*clhs14 + DN(0,1)*clhs16 + clhs20;
            lhs(0,4)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs24;
            lhs(0,5)=-clhs25;
            lhs(0,6)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs31;
            lhs(0,7)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + clhs35;
            lhs(0,8)=-clhs36;
            lhs(1,0)=DN(0,0)*clhs4 + DN(0,1)*clhs37 + clhs12;
            lhs(1,1)=DN(0,0)*clhs10 + DN(0,1)*clhs38 + clhs1 + clhs39*clhs7;
            lhs(1,2)=-clhs40;
            lhs(1,3)=DN(0,0)*clhs16 + DN(0,1)*clhs41 + clhs43;
            lhs(1,4)=DN(0,0)*clhs23 + DN(0,1)*clhs44 + clhs46;
            lhs(1,5)=-clhs47;
            lhs(1,6)=DN(0,0)*clhs28 + DN(0,1)*clhs48 + clhs49;
            lhs(1,7)=DN(0,0)*clhs34 + DN(0,1)*clhs50 + clhs52;
            lhs(1,8)=-clhs53;
            lhs(2,0)=clhs13*clhs56;
            lhs(2,1)=clhs40*clhs56;
            lhs(2,2)=clhs54*(clhs39 + clhs5);
            lhs(2,3)=clhs25*clhs55 + clhs57;
            lhs(2,4)=clhs47*clhs55 + clhs58;
            lhs(2,5)=clhs59;
            lhs(2,6)=clhs36*clhs55 + clhs60;
            lhs(2,7)=clhs53*clhs55 + clhs61;
            lhs(2,8)=clhs62;
            lhs(3,0)=DN(1,0)*clhs2 + DN(1,1)*clhs4 + clhs20;
            lhs(3,1)=DN(1,0)*clhs8 + DN(1,1)*clhs10 + clhs43;
            lhs(3,2)=-clhs57;
            lhs(3,3)=DN(1,0)*clhs14 + DN(1,1)*clhs16 + clhs63 + clhs64*clhs7;
            lhs(3,4)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs66;
            lhs(3,5)=-clhs67;
            lhs(3,6)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs70;
            lhs(3,7)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + clhs71;
            lhs(3,8)=-clhs72;
            lhs(4,0)=DN(1,0)*clhs4 + DN(1,1)*clhs37 + clhs24;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs38 + clhs46;
            lhs(4,2)=-clhs58;
            lhs(4,3)=DN(1,0)*clhs16 + DN(1,1)*clhs41 + clhs66;
            lhs(4,4)=DN(1,0)*clhs23 + DN(1,1)*clhs44 + clhs63 + clhs7*clhs73;
            lhs(4,5)=-clhs74;
            lhs(4,6)=DN(1,0)*clhs28 + DN(1,1)*clhs48 + clhs76;
            lhs(4,7)=DN(1,0)*clhs34 + DN(1,1)*clhs50 + clhs78;
            lhs(4,8)=-clhs79;
            lhs(5,0)=clhs25 + clhs55*clhs57;
            lhs(5,1)=clhs47 + clhs55*clhs58;
            lhs(5,2)=clhs59;
            lhs(5,3)=clhs56*clhs67;
            lhs(5,4)=clhs56*clhs74;
            lhs(5,5)=clhs54*(clhs64 + clhs73);
            lhs(5,6)=clhs55*clhs72 + clhs80;
            lhs(5,7)=clhs55*clhs79 + clhs81;
            lhs(5,8)=clhs82;
            lhs(6,0)=DN(2,0)*clhs2 + DN(2,1)*clhs4 + clhs31;
            lhs(6,1)=DN(2,0)*clhs8 + DN(2,1)*clhs10 + clhs49;
            lhs(6,2)=-clhs60;
            lhs(6,3)=DN(2,0)*clhs14 + DN(2,1)*clhs16 + clhs70;
            lhs(6,4)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs76;
            lhs(6,5)=-clhs80;
            lhs(6,6)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs7*clhs84 + clhs83;
            lhs(6,7)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + clhs85;
            lhs(6,8)=-clhs86;
            lhs(7,0)=DN(2,0)*clhs4 + DN(2,1)*clhs37 + clhs35;
            lhs(7,1)=DN(2,0)*clhs10 + DN(2,1)*clhs38 + clhs52;
            lhs(7,2)=-clhs61;
            lhs(7,3)=DN(2,0)*clhs16 + DN(2,1)*clhs41 + clhs71;
            lhs(7,4)=DN(2,0)*clhs23 + DN(2,1)*clhs44 + clhs78;
            lhs(7,5)=-clhs81;
            lhs(7,6)=DN(2,0)*clhs28 + DN(2,1)*clhs48 + clhs85;
            lhs(7,7)=DN(2,0)*clhs34 + DN(2,1)*clhs50 + clhs7*clhs87 + clhs83;
            lhs(7,8)=-clhs88;
            lhs(8,0)=clhs36 + clhs55*clhs60;
            lhs(8,1)=clhs53 + clhs55*clhs61;
            lhs(8,2)=clhs62;
            lhs(8,3)=clhs55*clhs80 + clhs72;
            lhs(8,4)=clhs55*clhs81 + clhs79;
            lhs(8,5)=clhs82;
            lhs(8,6)=clhs56*clhs86;
            lhs(8,7)=clhs56*clhs88;
            lhs(8,8)=clhs54*(clhs84 + clhs87);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointLHSContribution(
    SymbolicStokesData<3,4>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const array_1d<double, 3> v_aux = prod(trans(rData.Velocity), N);

    auto& lhs = rData.lhs;

    const double clhs0 =             bdf0*rho;
const double clhs1 =             pow(N[0], 2)*clhs0;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs3 =             C(0,3)*DN(0,0);
const double clhs4 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs3;
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             pow(DN(0,0), 2);
const double clhs8 =             rho*stab_c2*sqrt(pow(v_aux[0], 2) + pow(v_aux[1], 2) + pow(v_aux[2], 2));
const double clhs9 =             clhs8*h/stab_c1 + mu;
const double clhs10 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs3;
const double clhs11 =             C(1,3)*DN(0,1);
const double clhs12 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs11;
const double clhs13 =             C(3,5)*DN(0,0);
const double clhs14 =             C(4,5)*DN(0,2);
const double clhs15 =             C(1,5)*DN(0,1) + clhs13 + clhs14;
const double clhs16 =             DN(0,0)*clhs9;
const double clhs17 =             DN(0,1)*clhs16;
const double clhs18 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs19 =             C(3,4)*DN(0,1);
const double clhs20 =             C(2,3)*DN(0,2) + clhs13 + clhs19;
const double clhs21 =             C(2,5)*DN(0,2);
const double clhs22 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,2)*clhs16;
const double clhs24 =             DN(0,0)*N[0];
const double clhs25 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs26 =             C(0,3)*DN(1,0);
const double clhs27 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs26;
const double clhs28 =             C(0,5)*DN(1,0);
const double clhs29 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs28;
const double clhs30 =             DN(0,0)*DN(1,0);
const double clhs31 =             N[0]*clhs0;
const double clhs32 =             N[1]*clhs31;
const double clhs33 =             clhs30*clhs9 + clhs32;
const double clhs34 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs26;
const double clhs35 =             C(1,3)*DN(1,1);
const double clhs36 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs35;
const double clhs37 =             C(3,5)*DN(1,0);
const double clhs38 =             C(4,5)*DN(1,2);
const double clhs39 =             C(1,5)*DN(1,1) + clhs37 + clhs38;
const double clhs40 =             DN(1,1)*clhs16;
const double clhs41 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs28;
const double clhs42 =             C(3,4)*DN(1,1);
const double clhs43 =             C(2,3)*DN(1,2) + clhs37 + clhs42;
const double clhs44 =             C(2,5)*DN(1,2);
const double clhs45 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs44;
const double clhs46 =             DN(1,2)*clhs16;
const double clhs47 =             DN(0,0)*N[1];
const double clhs48 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs49 =             C(0,3)*DN(2,0);
const double clhs50 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs49;
const double clhs51 =             C(0,5)*DN(2,0);
const double clhs52 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs51;
const double clhs53 =             DN(0,0)*DN(2,0);
const double clhs54 =             N[2]*clhs31;
const double clhs55 =             clhs53*clhs9 + clhs54;
const double clhs56 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs49;
const double clhs57 =             C(1,3)*DN(2,1);
const double clhs58 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs57;
const double clhs59 =             C(3,5)*DN(2,0);
const double clhs60 =             C(4,5)*DN(2,2);
const double clhs61 =             C(1,5)*DN(2,1) + clhs59 + clhs60;
const double clhs62 =             DN(2,1)*clhs16;
const double clhs63 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs51;
const double clhs64 =             C(3,4)*DN(2,1);
const double clhs65 =             C(2,3)*DN(2,2) + clhs59 + clhs64;
const double clhs66 =             C(2,5)*DN(2,2);
const double clhs67 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs66;
const double clhs68 =             DN(2,2)*clhs16;
const double clhs69 =             DN(0,0)*N[2];
const double clhs70 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs71 =             C(0,3)*DN(3,0);
const double clhs72 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs71;
const double clhs73 =             C(0,5)*DN(3,0);
const double clhs74 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs73;
const double clhs75 =             DN(0,0)*DN(3,0);
const double clhs76 =             N[3]*clhs31;
const double clhs77 =             clhs75*clhs9 + clhs76;
const double clhs78 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs71;
const double clhs79 =             C(1,3)*DN(3,1);
const double clhs80 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs79;
const double clhs81 =             C(3,5)*DN(3,0);
const double clhs82 =             C(4,5)*DN(3,2);
const double clhs83 =             C(1,5)*DN(3,1) + clhs81 + clhs82;
const double clhs84 =             DN(3,1)*clhs16;
const double clhs85 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs73;
const double clhs86 =             C(3,4)*DN(3,1);
const double clhs87 =             C(2,3)*DN(3,2) + clhs81 + clhs86;
const double clhs88 =             C(2,5)*DN(3,2);
const double clhs89 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs88;
const double clhs90 =             DN(3,2)*clhs16;
const double clhs91 =             DN(0,0)*N[3];
const double clhs92 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs11;
const double clhs93 =             C(0,4)*DN(0,0) + clhs14 + clhs19;
const double clhs94 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs95 =             C(1,4)*DN(0,1);
const double clhs96 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs95;
const double clhs97 =             pow(DN(0,1), 2);
const double clhs98 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs95;
const double clhs99 =             C(2,4)*DN(0,2);
const double clhs100 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs99;
const double clhs101 =             DN(0,1)*clhs9;
const double clhs102 =             DN(0,2)*clhs101;
const double clhs103 =             DN(0,1)*N[0];
const double clhs104 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs35;
const double clhs105 =             C(0,4)*DN(1,0) + clhs38 + clhs42;
const double clhs106 =             DN(1,0)*clhs101;
const double clhs107 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs108 =             C(1,4)*DN(1,1);
const double clhs109 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs108;
const double clhs110 =             DN(0,1)*DN(1,1);
const double clhs111 =             clhs110*clhs9 + clhs32;
const double clhs112 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs108;
const double clhs113 =             C(2,4)*DN(1,2);
const double clhs114 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs113;
const double clhs115 =             DN(1,2)*clhs101;
const double clhs116 =             DN(0,1)*N[1];
const double clhs117 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs57;
const double clhs118 =             C(0,4)*DN(2,0) + clhs60 + clhs64;
const double clhs119 =             DN(2,0)*clhs101;
const double clhs120 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs121 =             C(1,4)*DN(2,1);
const double clhs122 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs121;
const double clhs123 =             DN(0,1)*DN(2,1);
const double clhs124 =             clhs123*clhs9 + clhs54;
const double clhs125 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs121;
const double clhs126 =             C(2,4)*DN(2,2);
const double clhs127 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs126;
const double clhs128 =             DN(2,2)*clhs101;
const double clhs129 =             DN(0,1)*N[2];
const double clhs130 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs79;
const double clhs131 =             C(0,4)*DN(3,0) + clhs82 + clhs86;
const double clhs132 =             DN(3,0)*clhs101;
const double clhs133 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs134 =             C(1,4)*DN(3,1);
const double clhs135 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs134;
const double clhs136 =             DN(0,1)*DN(3,1);
const double clhs137 =             clhs136*clhs9 + clhs76;
const double clhs138 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs134;
const double clhs139 =             C(2,4)*DN(3,2);
const double clhs140 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs139;
const double clhs141 =             DN(3,2)*clhs101;
const double clhs142 =             DN(0,1)*N[3];
const double clhs143 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs21;
const double clhs144 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs99;
const double clhs145 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs146 =             pow(DN(0,2), 2);
const double clhs147 =             DN(0,2)*N[0];
const double clhs148 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs44;
const double clhs149 =             DN(0,2)*clhs9;
const double clhs150 =             DN(1,0)*clhs149;
const double clhs151 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs113;
const double clhs152 =             DN(1,1)*clhs149;
const double clhs153 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs154 =             DN(0,2)*DN(1,2);
const double clhs155 =             clhs154*clhs9 + clhs32;
const double clhs156 =             DN(0,2)*N[1];
const double clhs157 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs66;
const double clhs158 =             DN(2,0)*clhs149;
const double clhs159 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs126;
const double clhs160 =             DN(2,1)*clhs149;
const double clhs161 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs162 =             DN(0,2)*DN(2,2);
const double clhs163 =             clhs162*clhs9 + clhs54;
const double clhs164 =             DN(0,2)*N[2];
const double clhs165 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs88;
const double clhs166 =             DN(3,0)*clhs149;
const double clhs167 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs139;
const double clhs168 =             DN(3,1)*clhs149;
const double clhs169 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs170 =             DN(0,2)*DN(3,2);
const double clhs171 =             clhs170*clhs9 + clhs76;
const double clhs172 =             DN(0,2)*N[3];
const double clhs173 =             1.0/(clhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs174 =             clhs0*clhs173;
const double clhs175 =             clhs174 + 1;
const double clhs176 =             DN(1,0)*N[0];
const double clhs177 =             DN(1,1)*N[0];
const double clhs178 =             DN(1,2)*N[0];
const double clhs179 =             clhs173*(clhs110 + clhs154 + clhs30);
const double clhs180 =             DN(2,0)*N[0];
const double clhs181 =             DN(2,1)*N[0];
const double clhs182 =             DN(2,2)*N[0];
const double clhs183 =             clhs173*(clhs123 + clhs162 + clhs53);
const double clhs184 =             DN(3,0)*N[0];
const double clhs185 =             DN(3,1)*N[0];
const double clhs186 =             DN(3,2)*N[0];
const double clhs187 =             clhs173*(clhs136 + clhs170 + clhs75);
const double clhs188 =             pow(N[1], 2)*clhs0;
const double clhs189 =             pow(DN(1,0), 2);
const double clhs190 =             DN(1,0)*clhs9;
const double clhs191 =             DN(1,1)*clhs190;
const double clhs192 =             DN(1,2)*clhs190;
const double clhs193 =             DN(1,0)*N[1];
const double clhs194 =             DN(1,0)*DN(2,0);
const double clhs195 =             N[1]*clhs0;
const double clhs196 =             N[2]*clhs195;
const double clhs197 =             clhs194*clhs9 + clhs196;
const double clhs198 =             DN(2,1)*clhs190;
const double clhs199 =             DN(2,2)*clhs190;
const double clhs200 =             DN(1,0)*N[2];
const double clhs201 =             DN(1,0)*DN(3,0);
const double clhs202 =             N[3]*clhs195;
const double clhs203 =             clhs201*clhs9 + clhs202;
const double clhs204 =             DN(3,1)*clhs190;
const double clhs205 =             DN(3,2)*clhs190;
const double clhs206 =             DN(1,0)*N[3];
const double clhs207 =             pow(DN(1,1), 2);
const double clhs208 =             DN(1,1)*clhs9;
const double clhs209 =             DN(1,2)*clhs208;
const double clhs210 =             DN(1,1)*N[1];
const double clhs211 =             DN(2,0)*clhs208;
const double clhs212 =             DN(1,1)*DN(2,1);
const double clhs213 =             clhs196 + clhs212*clhs9;
const double clhs214 =             DN(2,2)*clhs208;
const double clhs215 =             DN(1,1)*N[2];
const double clhs216 =             DN(3,0)*clhs208;
const double clhs217 =             DN(1,1)*DN(3,1);
const double clhs218 =             clhs202 + clhs217*clhs9;
const double clhs219 =             DN(3,2)*clhs208;
const double clhs220 =             DN(1,1)*N[3];
const double clhs221 =             pow(DN(1,2), 2);
const double clhs222 =             DN(1,2)*N[1];
const double clhs223 =             DN(1,2)*clhs9;
const double clhs224 =             DN(2,0)*clhs223;
const double clhs225 =             DN(2,1)*clhs223;
const double clhs226 =             DN(1,2)*DN(2,2);
const double clhs227 =             clhs196 + clhs226*clhs9;
const double clhs228 =             DN(1,2)*N[2];
const double clhs229 =             DN(3,0)*clhs223;
const double clhs230 =             DN(3,1)*clhs223;
const double clhs231 =             DN(1,2)*DN(3,2);
const double clhs232 =             clhs202 + clhs231*clhs9;
const double clhs233 =             DN(1,2)*N[3];
const double clhs234 =             DN(2,0)*N[1];
const double clhs235 =             DN(2,1)*N[1];
const double clhs236 =             DN(2,2)*N[1];
const double clhs237 =             clhs173*(clhs194 + clhs212 + clhs226);
const double clhs238 =             DN(3,0)*N[1];
const double clhs239 =             DN(3,1)*N[1];
const double clhs240 =             DN(3,2)*N[1];
const double clhs241 =             clhs173*(clhs201 + clhs217 + clhs231);
const double clhs242 =             pow(N[2], 2)*clhs0;
const double clhs243 =             pow(DN(2,0), 2);
const double clhs244 =             DN(2,0)*clhs9;
const double clhs245 =             DN(2,1)*clhs244;
const double clhs246 =             DN(2,2)*clhs244;
const double clhs247 =             DN(2,0)*N[2];
const double clhs248 =             DN(2,0)*DN(3,0);
const double clhs249 =             N[2]*N[3]*clhs0;
const double clhs250 =             clhs248*clhs9 + clhs249;
const double clhs251 =             DN(3,1)*clhs244;
const double clhs252 =             DN(3,2)*clhs244;
const double clhs253 =             DN(2,0)*N[3];
const double clhs254 =             pow(DN(2,1), 2);
const double clhs255 =             DN(2,1)*clhs9;
const double clhs256 =             DN(2,2)*clhs255;
const double clhs257 =             DN(2,1)*N[2];
const double clhs258 =             DN(3,0)*clhs255;
const double clhs259 =             DN(2,1)*DN(3,1);
const double clhs260 =             clhs249 + clhs259*clhs9;
const double clhs261 =             DN(3,2)*clhs255;
const double clhs262 =             DN(2,1)*N[3];
const double clhs263 =             pow(DN(2,2), 2);
const double clhs264 =             DN(2,2)*N[2];
const double clhs265 =             DN(2,2)*clhs9;
const double clhs266 =             DN(3,0)*clhs265;
const double clhs267 =             DN(3,1)*clhs265;
const double clhs268 =             DN(2,2)*DN(3,2);
const double clhs269 =             clhs249 + clhs268*clhs9;
const double clhs270 =             DN(2,2)*N[3];
const double clhs271 =             DN(3,0)*N[2];
const double clhs272 =             DN(3,1)*N[2];
const double clhs273 =             DN(3,2)*N[2];
const double clhs274 =             clhs173*(clhs248 + clhs259 + clhs268);
const double clhs275 =             pow(N[3], 2)*clhs0;
const double clhs276 =             pow(DN(3,0), 2);
const double clhs277 =             DN(3,0)*clhs9;
const double clhs278 =             DN(3,1)*clhs277;
const double clhs279 =             DN(3,2)*clhs277;
const double clhs280 =             DN(3,0)*N[3];
const double clhs281 =             pow(DN(3,1), 2);
const double clhs282 =             DN(3,1)*DN(3,2)*clhs9;
const double clhs283 =             DN(3,1)*N[3];
const double clhs284 =             pow(DN(3,2), 2);
const double clhs285 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs2 + DN(0,1)*clhs4 + DN(0,2)*clhs6 + clhs1 + clhs7*clhs9;
            lhs(0,1)=DN(0,0)*clhs10 + DN(0,1)*clhs12 + DN(0,2)*clhs15 + clhs17;
            lhs(0,2)=DN(0,0)*clhs18 + DN(0,1)*clhs20 + DN(0,2)*clhs22 + clhs23;
            lhs(0,3)=-clhs24;
            lhs(0,4)=DN(0,0)*clhs25 + DN(0,1)*clhs27 + DN(0,2)*clhs29 + clhs33;
            lhs(0,5)=DN(0,0)*clhs34 + DN(0,1)*clhs36 + DN(0,2)*clhs39 + clhs40;
            lhs(0,6)=DN(0,0)*clhs41 + DN(0,1)*clhs43 + DN(0,2)*clhs45 + clhs46;
            lhs(0,7)=-clhs47;
            lhs(0,8)=DN(0,0)*clhs48 + DN(0,1)*clhs50 + DN(0,2)*clhs52 + clhs55;
            lhs(0,9)=DN(0,0)*clhs56 + DN(0,1)*clhs58 + DN(0,2)*clhs61 + clhs62;
            lhs(0,10)=DN(0,0)*clhs63 + DN(0,1)*clhs65 + DN(0,2)*clhs67 + clhs68;
            lhs(0,11)=-clhs69;
            lhs(0,12)=DN(0,0)*clhs70 + DN(0,1)*clhs72 + DN(0,2)*clhs74 + clhs77;
            lhs(0,13)=DN(0,0)*clhs78 + DN(0,1)*clhs80 + DN(0,2)*clhs83 + clhs84;
            lhs(0,14)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs89 + clhs90;
            lhs(0,15)=-clhs91;
            lhs(1,0)=DN(0,0)*clhs4 + DN(0,1)*clhs92 + DN(0,2)*clhs93 + clhs17;
            lhs(1,1)=DN(0,0)*clhs12 + DN(0,1)*clhs94 + DN(0,2)*clhs96 + clhs1 + clhs9*clhs97;
            lhs(1,2)=DN(0,0)*clhs20 + DN(0,1)*clhs98 + DN(0,2)*clhs100 + clhs102;
            lhs(1,3)=-clhs103;
            lhs(1,4)=DN(0,0)*clhs27 + DN(0,1)*clhs104 + DN(0,2)*clhs105 + clhs106;
            lhs(1,5)=DN(0,0)*clhs36 + DN(0,1)*clhs107 + DN(0,2)*clhs109 + clhs111;
            lhs(1,6)=DN(0,0)*clhs43 + DN(0,1)*clhs112 + DN(0,2)*clhs114 + clhs115;
            lhs(1,7)=-clhs116;
            lhs(1,8)=DN(0,0)*clhs50 + DN(0,1)*clhs117 + DN(0,2)*clhs118 + clhs119;
            lhs(1,9)=DN(0,0)*clhs58 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs124;
            lhs(1,10)=DN(0,0)*clhs65 + DN(0,1)*clhs125 + DN(0,2)*clhs127 + clhs128;
            lhs(1,11)=-clhs129;
            lhs(1,12)=DN(0,0)*clhs72 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs132;
            lhs(1,13)=DN(0,0)*clhs80 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs137;
            lhs(1,14)=DN(0,0)*clhs87 + DN(0,1)*clhs138 + DN(0,2)*clhs140 + clhs141;
            lhs(1,15)=-clhs142;
            lhs(2,0)=DN(0,0)*clhs6 + DN(0,1)*clhs93 + DN(0,2)*clhs143 + clhs23;
            lhs(2,1)=DN(0,0)*clhs15 + DN(0,1)*clhs96 + DN(0,2)*clhs144 + clhs102;
            lhs(2,2)=DN(0,0)*clhs22 + DN(0,1)*clhs100 + DN(0,2)*clhs145 + clhs1 + clhs146*clhs9;
            lhs(2,3)=-clhs147;
            lhs(2,4)=DN(0,0)*clhs29 + DN(0,1)*clhs105 + DN(0,2)*clhs148 + clhs150;
            lhs(2,5)=DN(0,0)*clhs39 + DN(0,1)*clhs109 + DN(0,2)*clhs151 + clhs152;
            lhs(2,6)=DN(0,0)*clhs45 + DN(0,1)*clhs114 + DN(0,2)*clhs153 + clhs155;
            lhs(2,7)=-clhs156;
            lhs(2,8)=DN(0,0)*clhs52 + DN(0,1)*clhs118 + DN(0,2)*clhs157 + clhs158;
            lhs(2,9)=DN(0,0)*clhs61 + DN(0,1)*clhs122 + DN(0,2)*clhs159 + clhs160;
            lhs(2,10)=DN(0,0)*clhs67 + DN(0,1)*clhs127 + DN(0,2)*clhs161 + clhs163;
            lhs(2,11)=-clhs164;
            lhs(2,12)=DN(0,0)*clhs74 + DN(0,1)*clhs131 + DN(0,2)*clhs165 + clhs166;
            lhs(2,13)=DN(0,0)*clhs83 + DN(0,1)*clhs135 + DN(0,2)*clhs167 + clhs168;
            lhs(2,14)=DN(0,0)*clhs89 + DN(0,1)*clhs140 + DN(0,2)*clhs169 + clhs171;
            lhs(2,15)=-clhs172;
            lhs(3,0)=clhs175*clhs24;
            lhs(3,1)=clhs103*clhs175;
            lhs(3,2)=clhs147*clhs175;
            lhs(3,3)=clhs173*(clhs146 + clhs7 + clhs97);
            lhs(3,4)=clhs174*clhs47 + clhs176;
            lhs(3,5)=clhs116*clhs174 + clhs177;
            lhs(3,6)=clhs156*clhs174 + clhs178;
            lhs(3,7)=clhs179;
            lhs(3,8)=clhs174*clhs69 + clhs180;
            lhs(3,9)=clhs129*clhs174 + clhs181;
            lhs(3,10)=clhs164*clhs174 + clhs182;
            lhs(3,11)=clhs183;
            lhs(3,12)=clhs174*clhs91 + clhs184;
            lhs(3,13)=clhs142*clhs174 + clhs185;
            lhs(3,14)=clhs172*clhs174 + clhs186;
            lhs(3,15)=clhs187;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs4 + DN(1,2)*clhs6 + clhs33;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs12 + DN(1,2)*clhs15 + clhs106;
            lhs(4,2)=DN(1,0)*clhs18 + DN(1,1)*clhs20 + DN(1,2)*clhs22 + clhs150;
            lhs(4,3)=-clhs176;
            lhs(4,4)=DN(1,0)*clhs25 + DN(1,1)*clhs27 + DN(1,2)*clhs29 + clhs188 + clhs189*clhs9;
            lhs(4,5)=DN(1,0)*clhs34 + DN(1,1)*clhs36 + DN(1,2)*clhs39 + clhs191;
            lhs(4,6)=DN(1,0)*clhs41 + DN(1,1)*clhs43 + DN(1,2)*clhs45 + clhs192;
            lhs(4,7)=-clhs193;
            lhs(4,8)=DN(1,0)*clhs48 + DN(1,1)*clhs50 + DN(1,2)*clhs52 + clhs197;
            lhs(4,9)=DN(1,0)*clhs56 + DN(1,1)*clhs58 + DN(1,2)*clhs61 + clhs198;
            lhs(4,10)=DN(1,0)*clhs63 + DN(1,1)*clhs65 + DN(1,2)*clhs67 + clhs199;
            lhs(4,11)=-clhs200;
            lhs(4,12)=DN(1,0)*clhs70 + DN(1,1)*clhs72 + DN(1,2)*clhs74 + clhs203;
            lhs(4,13)=DN(1,0)*clhs78 + DN(1,1)*clhs80 + DN(1,2)*clhs83 + clhs204;
            lhs(4,14)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs89 + clhs205;
            lhs(4,15)=-clhs206;
            lhs(5,0)=DN(1,0)*clhs4 + DN(1,1)*clhs92 + DN(1,2)*clhs93 + clhs40;
            lhs(5,1)=DN(1,0)*clhs12 + DN(1,1)*clhs94 + DN(1,2)*clhs96 + clhs111;
            lhs(5,2)=DN(1,0)*clhs20 + DN(1,1)*clhs98 + DN(1,2)*clhs100 + clhs152;
            lhs(5,3)=-clhs177;
            lhs(5,4)=DN(1,0)*clhs27 + DN(1,1)*clhs104 + DN(1,2)*clhs105 + clhs191;
            lhs(5,5)=DN(1,0)*clhs36 + DN(1,1)*clhs107 + DN(1,2)*clhs109 + clhs188 + clhs207*clhs9;
            lhs(5,6)=DN(1,0)*clhs43 + DN(1,1)*clhs112 + DN(1,2)*clhs114 + clhs209;
            lhs(5,7)=-clhs210;
            lhs(5,8)=DN(1,0)*clhs50 + DN(1,1)*clhs117 + DN(1,2)*clhs118 + clhs211;
            lhs(5,9)=DN(1,0)*clhs58 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs213;
            lhs(5,10)=DN(1,0)*clhs65 + DN(1,1)*clhs125 + DN(1,2)*clhs127 + clhs214;
            lhs(5,11)=-clhs215;
            lhs(5,12)=DN(1,0)*clhs72 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs216;
            lhs(5,13)=DN(1,0)*clhs80 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs218;
            lhs(5,14)=DN(1,0)*clhs87 + DN(1,1)*clhs138 + DN(1,2)*clhs140 + clhs219;
            lhs(5,15)=-clhs220;
            lhs(6,0)=DN(1,0)*clhs6 + DN(1,1)*clhs93 + DN(1,2)*clhs143 + clhs46;
            lhs(6,1)=DN(1,0)*clhs15 + DN(1,1)*clhs96 + DN(1,2)*clhs144 + clhs115;
            lhs(6,2)=DN(1,0)*clhs22 + DN(1,1)*clhs100 + DN(1,2)*clhs145 + clhs155;
            lhs(6,3)=-clhs178;
            lhs(6,4)=DN(1,0)*clhs29 + DN(1,1)*clhs105 + DN(1,2)*clhs148 + clhs192;
            lhs(6,5)=DN(1,0)*clhs39 + DN(1,1)*clhs109 + DN(1,2)*clhs151 + clhs209;
            lhs(6,6)=DN(1,0)*clhs45 + DN(1,1)*clhs114 + DN(1,2)*clhs153 + clhs188 + clhs221*clhs9;
            lhs(6,7)=-clhs222;
            lhs(6,8)=DN(1,0)*clhs52 + DN(1,1)*clhs118 + DN(1,2)*clhs157 + clhs224;
            lhs(6,9)=DN(1,0)*clhs61 + DN(1,1)*clhs122 + DN(1,2)*clhs159 + clhs225;
            lhs(6,10)=DN(1,0)*clhs67 + DN(1,1)*clhs127 + DN(1,2)*clhs161 + clhs227;
            lhs(6,11)=-clhs228;
            lhs(6,12)=DN(1,0)*clhs74 + DN(1,1)*clhs131 + DN(1,2)*clhs165 + clhs229;
            lhs(6,13)=DN(1,0)*clhs83 + DN(1,1)*clhs135 + DN(1,2)*clhs167 + clhs230;
            lhs(6,14)=DN(1,0)*clhs89 + DN(1,1)*clhs140 + DN(1,2)*clhs169 + clhs232;
            lhs(6,15)=-clhs233;
            lhs(7,0)=clhs174*clhs176 + clhs47;
            lhs(7,1)=clhs116 + clhs174*clhs177;
            lhs(7,2)=clhs156 + clhs174*clhs178;
            lhs(7,3)=clhs179;
            lhs(7,4)=clhs175*clhs193;
            lhs(7,5)=clhs175*clhs210;
            lhs(7,6)=clhs175*clhs222;
            lhs(7,7)=clhs173*(clhs189 + clhs207 + clhs221);
            lhs(7,8)=clhs174*clhs200 + clhs234;
            lhs(7,9)=clhs174*clhs215 + clhs235;
            lhs(7,10)=clhs174*clhs228 + clhs236;
            lhs(7,11)=clhs237;
            lhs(7,12)=clhs174*clhs206 + clhs238;
            lhs(7,13)=clhs174*clhs220 + clhs239;
            lhs(7,14)=clhs174*clhs233 + clhs240;
            lhs(7,15)=clhs241;
            lhs(8,0)=DN(2,0)*clhs2 + DN(2,1)*clhs4 + DN(2,2)*clhs6 + clhs55;
            lhs(8,1)=DN(2,0)*clhs10 + DN(2,1)*clhs12 + DN(2,2)*clhs15 + clhs119;
            lhs(8,2)=DN(2,0)*clhs18 + DN(2,1)*clhs20 + DN(2,2)*clhs22 + clhs158;
            lhs(8,3)=-clhs180;
            lhs(8,4)=DN(2,0)*clhs25 + DN(2,1)*clhs27 + DN(2,2)*clhs29 + clhs197;
            lhs(8,5)=DN(2,0)*clhs34 + DN(2,1)*clhs36 + DN(2,2)*clhs39 + clhs211;
            lhs(8,6)=DN(2,0)*clhs41 + DN(2,1)*clhs43 + DN(2,2)*clhs45 + clhs224;
            lhs(8,7)=-clhs234;
            lhs(8,8)=DN(2,0)*clhs48 + DN(2,1)*clhs50 + DN(2,2)*clhs52 + clhs242 + clhs243*clhs9;
            lhs(8,9)=DN(2,0)*clhs56 + DN(2,1)*clhs58 + DN(2,2)*clhs61 + clhs245;
            lhs(8,10)=DN(2,0)*clhs63 + DN(2,1)*clhs65 + DN(2,2)*clhs67 + clhs246;
            lhs(8,11)=-clhs247;
            lhs(8,12)=DN(2,0)*clhs70 + DN(2,1)*clhs72 + DN(2,2)*clhs74 + clhs250;
            lhs(8,13)=DN(2,0)*clhs78 + DN(2,1)*clhs80 + DN(2,2)*clhs83 + clhs251;
            lhs(8,14)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs89 + clhs252;
            lhs(8,15)=-clhs253;
            lhs(9,0)=DN(2,0)*clhs4 + DN(2,1)*clhs92 + DN(2,2)*clhs93 + clhs62;
            lhs(9,1)=DN(2,0)*clhs12 + DN(2,1)*clhs94 + DN(2,2)*clhs96 + clhs124;
            lhs(9,2)=DN(2,0)*clhs20 + DN(2,1)*clhs98 + DN(2,2)*clhs100 + clhs160;
            lhs(9,3)=-clhs181;
            lhs(9,4)=DN(2,0)*clhs27 + DN(2,1)*clhs104 + DN(2,2)*clhs105 + clhs198;
            lhs(9,5)=DN(2,0)*clhs36 + DN(2,1)*clhs107 + DN(2,2)*clhs109 + clhs213;
            lhs(9,6)=DN(2,0)*clhs43 + DN(2,1)*clhs112 + DN(2,2)*clhs114 + clhs225;
            lhs(9,7)=-clhs235;
            lhs(9,8)=DN(2,0)*clhs50 + DN(2,1)*clhs117 + DN(2,2)*clhs118 + clhs245;
            lhs(9,9)=DN(2,0)*clhs58 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs242 + clhs254*clhs9;
            lhs(9,10)=DN(2,0)*clhs65 + DN(2,1)*clhs125 + DN(2,2)*clhs127 + clhs256;
            lhs(9,11)=-clhs257;
            lhs(9,12)=DN(2,0)*clhs72 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs258;
            lhs(9,13)=DN(2,0)*clhs80 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs260;
            lhs(9,14)=DN(2,0)*clhs87 + DN(2,1)*clhs138 + DN(2,2)*clhs140 + clhs261;
            lhs(9,15)=-clhs262;
            lhs(10,0)=DN(2,0)*clhs6 + DN(2,1)*clhs93 + DN(2,2)*clhs143 + clhs68;
            lhs(10,1)=DN(2,0)*clhs15 + DN(2,1)*clhs96 + DN(2,2)*clhs144 + clhs128;
            lhs(10,2)=DN(2,0)*clhs22 + DN(2,1)*clhs100 + DN(2,2)*clhs145 + clhs163;
            lhs(10,3)=-clhs182;
            lhs(10,4)=DN(2,0)*clhs29 + DN(2,1)*clhs105 + DN(2,2)*clhs148 + clhs199;
            lhs(10,5)=DN(2,0)*clhs39 + DN(2,1)*clhs109 + DN(2,2)*clhs151 + clhs214;
            lhs(10,6)=DN(2,0)*clhs45 + DN(2,1)*clhs114 + DN(2,2)*clhs153 + clhs227;
            lhs(10,7)=-clhs236;
            lhs(10,8)=DN(2,0)*clhs52 + DN(2,1)*clhs118 + DN(2,2)*clhs157 + clhs246;
            lhs(10,9)=DN(2,0)*clhs61 + DN(2,1)*clhs122 + DN(2,2)*clhs159 + clhs256;
            lhs(10,10)=DN(2,0)*clhs67 + DN(2,1)*clhs127 + DN(2,2)*clhs161 + clhs242 + clhs263*clhs9;
            lhs(10,11)=-clhs264;
            lhs(10,12)=DN(2,0)*clhs74 + DN(2,1)*clhs131 + DN(2,2)*clhs165 + clhs266;
            lhs(10,13)=DN(2,0)*clhs83 + DN(2,1)*clhs135 + DN(2,2)*clhs167 + clhs267;
            lhs(10,14)=DN(2,0)*clhs89 + DN(2,1)*clhs140 + DN(2,2)*clhs169 + clhs269;
            lhs(10,15)=-clhs270;
            lhs(11,0)=clhs174*clhs180 + clhs69;
            lhs(11,1)=clhs129 + clhs174*clhs181;
            lhs(11,2)=clhs164 + clhs174*clhs182;
            lhs(11,3)=clhs183;
            lhs(11,4)=clhs174*clhs234 + clhs200;
            lhs(11,5)=clhs174*clhs235 + clhs215;
            lhs(11,6)=clhs174*clhs236 + clhs228;
            lhs(11,7)=clhs237;
            lhs(11,8)=clhs175*clhs247;
            lhs(11,9)=clhs175*clhs257;
            lhs(11,10)=clhs175*clhs264;
            lhs(11,11)=clhs173*(clhs243 + clhs254 + clhs263);
            lhs(11,12)=clhs174*clhs253 + clhs271;
            lhs(11,13)=clhs174*clhs262 + clhs272;
            lhs(11,14)=clhs174*clhs270 + clhs273;
            lhs(11,15)=clhs274;
            lhs(12,0)=DN(3,0)*clhs2 + DN(3,1)*clhs4 + DN(3,2)*clhs6 + clhs77;
            lhs(12,1)=DN(3,0)*clhs10 + DN(3,1)*clhs12 + DN(3,2)*clhs15 + clhs132;
            lhs(12,2)=DN(3,0)*clhs18 + DN(3,1)*clhs20 + DN(3,2)*clhs22 + clhs166;
            lhs(12,3)=-clhs184;
            lhs(12,4)=DN(3,0)*clhs25 + DN(3,1)*clhs27 + DN(3,2)*clhs29 + clhs203;
            lhs(12,5)=DN(3,0)*clhs34 + DN(3,1)*clhs36 + DN(3,2)*clhs39 + clhs216;
            lhs(12,6)=DN(3,0)*clhs41 + DN(3,1)*clhs43 + DN(3,2)*clhs45 + clhs229;
            lhs(12,7)=-clhs238;
            lhs(12,8)=DN(3,0)*clhs48 + DN(3,1)*clhs50 + DN(3,2)*clhs52 + clhs250;
            lhs(12,9)=DN(3,0)*clhs56 + DN(3,1)*clhs58 + DN(3,2)*clhs61 + clhs258;
            lhs(12,10)=DN(3,0)*clhs63 + DN(3,1)*clhs65 + DN(3,2)*clhs67 + clhs266;
            lhs(12,11)=-clhs271;
            lhs(12,12)=DN(3,0)*clhs70 + DN(3,1)*clhs72 + DN(3,2)*clhs74 + clhs275 + clhs276*clhs9;
            lhs(12,13)=DN(3,0)*clhs78 + DN(3,1)*clhs80 + DN(3,2)*clhs83 + clhs278;
            lhs(12,14)=DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs89 + clhs279;
            lhs(12,15)=-clhs280;
            lhs(13,0)=DN(3,0)*clhs4 + DN(3,1)*clhs92 + DN(3,2)*clhs93 + clhs84;
            lhs(13,1)=DN(3,0)*clhs12 + DN(3,1)*clhs94 + DN(3,2)*clhs96 + clhs137;
            lhs(13,2)=DN(3,0)*clhs20 + DN(3,1)*clhs98 + DN(3,2)*clhs100 + clhs168;
            lhs(13,3)=-clhs185;
            lhs(13,4)=DN(3,0)*clhs27 + DN(3,1)*clhs104 + DN(3,2)*clhs105 + clhs204;
            lhs(13,5)=DN(3,0)*clhs36 + DN(3,1)*clhs107 + DN(3,2)*clhs109 + clhs218;
            lhs(13,6)=DN(3,0)*clhs43 + DN(3,1)*clhs112 + DN(3,2)*clhs114 + clhs230;
            lhs(13,7)=-clhs239;
            lhs(13,8)=DN(3,0)*clhs50 + DN(3,1)*clhs117 + DN(3,2)*clhs118 + clhs251;
            lhs(13,9)=DN(3,0)*clhs58 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs260;
            lhs(13,10)=DN(3,0)*clhs65 + DN(3,1)*clhs125 + DN(3,2)*clhs127 + clhs267;
            lhs(13,11)=-clhs272;
            lhs(13,12)=DN(3,0)*clhs72 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs278;
            lhs(13,13)=DN(3,0)*clhs80 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs275 + clhs281*clhs9;
            lhs(13,14)=DN(3,0)*clhs87 + DN(3,1)*clhs138 + DN(3,2)*clhs140 + clhs282;
            lhs(13,15)=-clhs283;
            lhs(14,0)=DN(3,0)*clhs6 + DN(3,1)*clhs93 + DN(3,2)*clhs143 + clhs90;
            lhs(14,1)=DN(3,0)*clhs15 + DN(3,1)*clhs96 + DN(3,2)*clhs144 + clhs141;
            lhs(14,2)=DN(3,0)*clhs22 + DN(3,1)*clhs100 + DN(3,2)*clhs145 + clhs171;
            lhs(14,3)=-clhs186;
            lhs(14,4)=DN(3,0)*clhs29 + DN(3,1)*clhs105 + DN(3,2)*clhs148 + clhs205;
            lhs(14,5)=DN(3,0)*clhs39 + DN(3,1)*clhs109 + DN(3,2)*clhs151 + clhs219;
            lhs(14,6)=DN(3,0)*clhs45 + DN(3,1)*clhs114 + DN(3,2)*clhs153 + clhs232;
            lhs(14,7)=-clhs240;
            lhs(14,8)=DN(3,0)*clhs52 + DN(3,1)*clhs118 + DN(3,2)*clhs157 + clhs252;
            lhs(14,9)=DN(3,0)*clhs61 + DN(3,1)*clhs122 + DN(3,2)*clhs159 + clhs261;
            lhs(14,10)=DN(3,0)*clhs67 + DN(3,1)*clhs127 + DN(3,2)*clhs161 + clhs269;
            lhs(14,11)=-clhs273;
            lhs(14,12)=DN(3,0)*clhs74 + DN(3,1)*clhs131 + DN(3,2)*clhs165 + clhs279;
            lhs(14,13)=DN(3,0)*clhs83 + DN(3,1)*clhs135 + DN(3,2)*clhs167 + clhs282;
            lhs(14,14)=DN(3,0)*clhs89 + DN(3,1)*clhs140 + DN(3,2)*clhs169 + clhs275 + clhs284*clhs9;
            lhs(14,15)=-clhs285;
            lhs(15,0)=clhs174*clhs184 + clhs91;
            lhs(15,1)=clhs142 + clhs174*clhs185;
            lhs(15,2)=clhs172 + clhs174*clhs186;
            lhs(15,3)=clhs187;
            lhs(15,4)=clhs174*clhs238 + clhs206;
            lhs(15,5)=clhs174*clhs239 + clhs220;
            lhs(15,6)=clhs174*clhs240 + clhs233;
            lhs(15,7)=clhs241;
            lhs(15,8)=clhs174*clhs271 + clhs253;
            lhs(15,9)=clhs174*clhs272 + clhs262;
            lhs(15,10)=clhs174*clhs273 + clhs270;
            lhs(15,11)=clhs274;
            lhs(15,12)=clhs175*clhs280;
            lhs(15,13)=clhs175*clhs283;
            lhs(15,14)=clhs175*clhs285;
            lhs(15,15)=clhs173*(clhs276 + clhs281 + clhs284);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void SymbolicStokes<SymbolicStokesData<2,3>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const array_1d<double, 2> v_aux = prod(trans(rData.Velocity), N);

    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs3 =             rho*stab_c2*sqrt(pow(v_aux[0], 2) + pow(v_aux[1], 2));
const double crhs4 =             crhs2*(crhs3*h/stab_c1 + mu);
const double crhs5 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs6 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs7 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs8 =             1.0/(crhs3/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs9 =             crhs8*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs5);
const double crhs10 =             crhs8*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs6 + crhs7);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs4 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs5;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs4 - DN(0,1)*stress[1] + N[0]*crhs6 - N[0]*crhs7;
            rhs[2]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 - N[0]*crhs2;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs4 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs5;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs4 - DN(1,1)*stress[1] + N[1]*crhs6 - N[1]*crhs7;
            rhs[5]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 - N[1]*crhs2;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs4 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs5;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs4 - DN(2,1)*stress[1] + N[2]*crhs6 - N[2]*crhs7;
            rhs[8]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 - N[2]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void SymbolicStokes<SymbolicStokesData<3,4>>::ComputeGaussPointRHSContribution(
    SymbolicStokesData<3,4>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    const array_1d<double, 3> v_aux = prod(trans(rData.Velocity), N);

    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs4 =             rho*stab_c2*sqrt(pow(v_aux[0], 2) + pow(v_aux[1], 2) + pow(v_aux[2], 2));
const double crhs5 =             crhs3*(crhs4*h/stab_c1 + mu);
const double crhs6 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs7 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs8 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs9 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs10 =             1.0/(crhs4/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs11 =             crhs10*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs2);
const double crhs12 =             crhs10*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs6 + crhs7);
const double crhs13 =             crhs10*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs8 + crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs2;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs5 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs6 - N[0]*crhs7;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs5 - DN(0,2)*stress[2] + N[0]*crhs8 - N[0]*crhs9;
            rhs[3]=-DN(0,0)*crhs11 - DN(0,1)*crhs12 - DN(0,2)*crhs13 - N[0]*crhs3;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs2;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs5 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs6 - N[1]*crhs7;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs5 - DN(1,2)*stress[2] + N[1]*crhs8 - N[1]*crhs9;
            rhs[7]=-DN(1,0)*crhs11 - DN(1,1)*crhs12 - DN(1,2)*crhs13 - N[1]*crhs3;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs2;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs5 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs6 - N[2]*crhs7;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs5 - DN(2,2)*stress[2] + N[2]*crhs8 - N[2]*crhs9;
            rhs[11]=-DN(2,0)*crhs11 - DN(2,1)*crhs12 - DN(2,2)*crhs13 - N[2]*crhs3;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs5 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs2;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs5 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs6 - N[3]*crhs7;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs5 - DN(3,2)*stress[2] + N[3]*crhs8 - N[3]*crhs9;
            rhs[15]=-DN(3,0)*crhs11 - DN(3,1)*crhs12 - DN(3,2)*crhs13 - N[3]*crhs3;


    noalias(rRHS) += rData.Weight * rhs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void SymbolicStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void SymbolicStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class SymbolicStokes< SymbolicStokesData<2,3> >;
template class SymbolicStokes< SymbolicStokesData<3,4> >;

}