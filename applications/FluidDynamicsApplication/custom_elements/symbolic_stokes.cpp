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

    auto& lhs = rData.lhs;

    const double clhs0 =             bdf0*rho;
const double clhs1 =             pow(N[0], 2)*clhs0;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs3 =             C(0,2)*DN(0,0);
const double clhs4 =             C(2,2)*DN(0,1) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             pow(h, 2);
const double clhs7 =             dyn_tau*rho/dt + mu*stab_c1/clhs6;
const double clhs8 =             1.0*clhs6*clhs7/stab_c1;
const double clhs9 =             C(0,1)*DN(0,1) + clhs3;
const double clhs10 =             C(1,2)*DN(0,1);
const double clhs11 =             C(2,2)*DN(0,0) + clhs10;
const double clhs12 =             DN(0,0)*clhs8;
const double clhs13 =             DN(0,1)*clhs12;
const double clhs14 =             DN(0,0)*N[0];
const double clhs15 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs16 =             C(0,2)*DN(1,0);
const double clhs17 =             C(2,2)*DN(1,1) + clhs16;
const double clhs18 =             N[0]*clhs0;
const double clhs19 =             N[1]*clhs18;
const double clhs20 =             DN(0,0)*DN(1,0);
const double clhs21 =             clhs19 + clhs20*clhs8;
const double clhs22 =             C(0,1)*DN(1,1) + clhs16;
const double clhs23 =             C(1,2)*DN(1,1);
const double clhs24 =             C(2,2)*DN(1,0) + clhs23;
const double clhs25 =             DN(1,1)*clhs12;
const double clhs26 =             DN(0,0)*N[1];
const double clhs27 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs28 =             C(0,2)*DN(2,0);
const double clhs29 =             C(2,2)*DN(2,1) + clhs28;
const double clhs30 =             N[2]*clhs18;
const double clhs31 =             DN(0,0)*DN(2,0);
const double clhs32 =             clhs30 + clhs31*clhs8;
const double clhs33 =             C(0,1)*DN(2,1) + clhs28;
const double clhs34 =             C(1,2)*DN(2,1);
const double clhs35 =             C(2,2)*DN(2,0) + clhs34;
const double clhs36 =             DN(2,1)*clhs12;
const double clhs37 =             DN(0,0)*N[2];
const double clhs38 =             C(0,1)*DN(0,0) + clhs10;
const double clhs39 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs40 =             pow(DN(0,1), 2);
const double clhs41 =             DN(0,1)*N[0];
const double clhs42 =             C(0,1)*DN(1,0) + clhs23;
const double clhs43 =             DN(0,1)*clhs8;
const double clhs44 =             DN(1,0)*clhs43;
const double clhs45 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs46 =             DN(0,1)*DN(1,1);
const double clhs47 =             clhs19 + clhs46*clhs8;
const double clhs48 =             DN(0,1)*N[1];
const double clhs49 =             C(0,1)*DN(2,0) + clhs34;
const double clhs50 =             DN(2,0)*clhs43;
const double clhs51 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs52 =             DN(0,1)*DN(2,1);
const double clhs53 =             clhs30 + clhs52*clhs8;
const double clhs54 =             DN(0,1)*N[2];
const double clhs55 =             1.0/clhs7;
const double clhs56 =             clhs0*clhs55;
const double clhs57 =             clhs56 + 1;
const double clhs58 =             DN(1,0)*N[0];
const double clhs59 =             DN(1,1)*N[0];
const double clhs60 =             clhs55*(clhs20 + clhs46);
const double clhs61 =             DN(2,0)*N[0];
const double clhs62 =             DN(2,1)*N[0];
const double clhs63 =             clhs55*(clhs31 + clhs52);
const double clhs64 =             pow(N[1], 2)*clhs0;
const double clhs65 =             pow(DN(1,0), 2);
const double clhs66 =             DN(1,0)*clhs8;
const double clhs67 =             DN(1,1)*clhs66;
const double clhs68 =             DN(1,0)*N[1];
const double clhs69 =             N[1]*N[2]*clhs0;
const double clhs70 =             DN(1,0)*DN(2,0);
const double clhs71 =             clhs69 + clhs70*clhs8;
const double clhs72 =             DN(2,1)*clhs66;
const double clhs73 =             DN(1,0)*N[2];
const double clhs74 =             pow(DN(1,1), 2);
const double clhs75 =             DN(1,1)*N[1];
const double clhs76 =             DN(2,0)*clhs8;
const double clhs77 =             DN(1,1)*clhs76;
const double clhs78 =             DN(1,1)*DN(2,1);
const double clhs79 =             clhs69 + clhs78*clhs8;
const double clhs80 =             DN(1,1)*N[2];
const double clhs81 =             DN(2,0)*N[1];
const double clhs82 =             DN(2,1)*N[1];
const double clhs83 =             clhs55*(clhs70 + clhs78);
const double clhs84 =             pow(N[2], 2)*clhs0;
const double clhs85 =             pow(DN(2,0), 2);
const double clhs86 =             DN(2,1)*clhs76;
const double clhs87 =             DN(2,0)*N[2];
const double clhs88 =             pow(DN(2,1), 2);
const double clhs89 =             DN(2,1)*N[2];
            lhs(0,0)=DN(0,0)*clhs2 + DN(0,1)*clhs4 + clhs1 + clhs5*clhs8;
            lhs(0,1)=DN(0,0)*clhs9 + DN(0,1)*clhs11 + clhs13;
            lhs(0,2)=-clhs14;
            lhs(0,3)=DN(0,0)*clhs15 + DN(0,1)*clhs17 + clhs21;
            lhs(0,4)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + clhs25;
            lhs(0,5)=-clhs26;
            lhs(0,6)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + clhs32;
            lhs(0,7)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + clhs36;
            lhs(0,8)=-clhs37;
            lhs(1,0)=DN(0,0)*clhs4 + DN(0,1)*clhs38 + clhs13;
            lhs(1,1)=DN(0,0)*clhs11 + DN(0,1)*clhs39 + clhs1 + clhs40*clhs8;
            lhs(1,2)=-clhs41;
            lhs(1,3)=DN(0,0)*clhs17 + DN(0,1)*clhs42 + clhs44;
            lhs(1,4)=DN(0,0)*clhs24 + DN(0,1)*clhs45 + clhs47;
            lhs(1,5)=-clhs48;
            lhs(1,6)=DN(0,0)*clhs29 + DN(0,1)*clhs49 + clhs50;
            lhs(1,7)=DN(0,0)*clhs35 + DN(0,1)*clhs51 + clhs53;
            lhs(1,8)=-clhs54;
            lhs(2,0)=clhs14*clhs57;
            lhs(2,1)=clhs41*clhs57;
            lhs(2,2)=clhs55*(clhs40 + clhs5);
            lhs(2,3)=clhs26*clhs56 + clhs58;
            lhs(2,4)=clhs48*clhs56 + clhs59;
            lhs(2,5)=clhs60;
            lhs(2,6)=clhs37*clhs56 + clhs61;
            lhs(2,7)=clhs54*clhs56 + clhs62;
            lhs(2,8)=clhs63;
            lhs(3,0)=DN(1,0)*clhs2 + DN(1,1)*clhs4 + clhs21;
            lhs(3,1)=DN(1,0)*clhs9 + DN(1,1)*clhs11 + clhs44;
            lhs(3,2)=-clhs58;
            lhs(3,3)=DN(1,0)*clhs15 + DN(1,1)*clhs17 + clhs64 + clhs65*clhs8;
            lhs(3,4)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + clhs67;
            lhs(3,5)=-clhs68;
            lhs(3,6)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs71;
            lhs(3,7)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + clhs72;
            lhs(3,8)=-clhs73;
            lhs(4,0)=DN(1,0)*clhs4 + DN(1,1)*clhs38 + clhs25;
            lhs(4,1)=DN(1,0)*clhs11 + DN(1,1)*clhs39 + clhs47;
            lhs(4,2)=-clhs59;
            lhs(4,3)=DN(1,0)*clhs17 + DN(1,1)*clhs42 + clhs67;
            lhs(4,4)=DN(1,0)*clhs24 + DN(1,1)*clhs45 + clhs64 + clhs74*clhs8;
            lhs(4,5)=-clhs75;
            lhs(4,6)=DN(1,0)*clhs29 + DN(1,1)*clhs49 + clhs77;
            lhs(4,7)=DN(1,0)*clhs35 + DN(1,1)*clhs51 + clhs79;
            lhs(4,8)=-clhs80;
            lhs(5,0)=clhs26 + clhs56*clhs58;
            lhs(5,1)=clhs48 + clhs56*clhs59;
            lhs(5,2)=clhs60;
            lhs(5,3)=clhs57*clhs68;
            lhs(5,4)=clhs57*clhs75;
            lhs(5,5)=clhs55*(clhs65 + clhs74);
            lhs(5,6)=clhs56*clhs73 + clhs81;
            lhs(5,7)=clhs56*clhs80 + clhs82;
            lhs(5,8)=clhs83;
            lhs(6,0)=DN(2,0)*clhs2 + DN(2,1)*clhs4 + clhs32;
            lhs(6,1)=DN(2,0)*clhs9 + DN(2,1)*clhs11 + clhs50;
            lhs(6,2)=-clhs61;
            lhs(6,3)=DN(2,0)*clhs15 + DN(2,1)*clhs17 + clhs71;
            lhs(6,4)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + clhs77;
            lhs(6,5)=-clhs81;
            lhs(6,6)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs8*clhs85 + clhs84;
            lhs(6,7)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + clhs86;
            lhs(6,8)=-clhs87;
            lhs(7,0)=DN(2,0)*clhs4 + DN(2,1)*clhs38 + clhs36;
            lhs(7,1)=DN(2,0)*clhs11 + DN(2,1)*clhs39 + clhs53;
            lhs(7,2)=-clhs62;
            lhs(7,3)=DN(2,0)*clhs17 + DN(2,1)*clhs42 + clhs72;
            lhs(7,4)=DN(2,0)*clhs24 + DN(2,1)*clhs45 + clhs79;
            lhs(7,5)=-clhs82;
            lhs(7,6)=DN(2,0)*clhs29 + DN(2,1)*clhs49 + clhs86;
            lhs(7,7)=DN(2,0)*clhs35 + DN(2,1)*clhs51 + clhs8*clhs88 + clhs84;
            lhs(7,8)=-clhs89;
            lhs(8,0)=clhs37 + clhs56*clhs61;
            lhs(8,1)=clhs54 + clhs56*clhs62;
            lhs(8,2)=clhs63;
            lhs(8,3)=clhs56*clhs81 + clhs73;
            lhs(8,4)=clhs56*clhs82 + clhs80;
            lhs(8,5)=clhs83;
            lhs(8,6)=clhs57*clhs87;
            lhs(8,7)=clhs57*clhs89;
            lhs(8,8)=clhs55*(clhs85 + clhs88);


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

    auto& lhs = rData.lhs;

    const double clhs0 =             bdf0*rho;
const double clhs1 =             pow(N[0], 2)*clhs0;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs3 =             C(0,3)*DN(0,0);
const double clhs4 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs3;
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             pow(DN(0,0), 2);
const double clhs8 =             pow(h, 2);
const double clhs9 =             dyn_tau*rho/dt + mu*stab_c1/clhs8;
const double clhs10 =             1.0*clhs8*clhs9/stab_c1;
const double clhs11 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs3;
const double clhs12 =             C(1,3)*DN(0,1);
const double clhs13 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs12;
const double clhs14 =             C(3,5)*DN(0,0);
const double clhs15 =             C(4,5)*DN(0,2);
const double clhs16 =             C(1,5)*DN(0,1) + clhs14 + clhs15;
const double clhs17 =             DN(0,0)*clhs10;
const double clhs18 =             DN(0,1)*clhs17;
const double clhs19 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs20 =             C(3,4)*DN(0,1);
const double clhs21 =             C(2,3)*DN(0,2) + clhs14 + clhs20;
const double clhs22 =             C(2,5)*DN(0,2);
const double clhs23 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs22;
const double clhs24 =             DN(0,2)*clhs17;
const double clhs25 =             DN(0,0)*N[0];
const double clhs26 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs27 =             C(0,3)*DN(1,0);
const double clhs28 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs27;
const double clhs29 =             C(0,5)*DN(1,0);
const double clhs30 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs29;
const double clhs31 =             N[0]*clhs0;
const double clhs32 =             N[1]*clhs31;
const double clhs33 =             DN(0,0)*DN(1,0);
const double clhs34 =             clhs10*clhs33 + clhs32;
const double clhs35 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs27;
const double clhs36 =             C(1,3)*DN(1,1);
const double clhs37 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs36;
const double clhs38 =             C(3,5)*DN(1,0);
const double clhs39 =             C(4,5)*DN(1,2);
const double clhs40 =             C(1,5)*DN(1,1) + clhs38 + clhs39;
const double clhs41 =             DN(1,1)*clhs17;
const double clhs42 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs29;
const double clhs43 =             C(3,4)*DN(1,1);
const double clhs44 =             C(2,3)*DN(1,2) + clhs38 + clhs43;
const double clhs45 =             C(2,5)*DN(1,2);
const double clhs46 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs45;
const double clhs47 =             DN(1,2)*clhs17;
const double clhs48 =             DN(0,0)*N[1];
const double clhs49 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs50 =             C(0,3)*DN(2,0);
const double clhs51 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs50;
const double clhs52 =             C(0,5)*DN(2,0);
const double clhs53 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs52;
const double clhs54 =             N[2]*clhs31;
const double clhs55 =             DN(0,0)*DN(2,0);
const double clhs56 =             clhs10*clhs55 + clhs54;
const double clhs57 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs50;
const double clhs58 =             C(1,3)*DN(2,1);
const double clhs59 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs58;
const double clhs60 =             C(3,5)*DN(2,0);
const double clhs61 =             C(4,5)*DN(2,2);
const double clhs62 =             C(1,5)*DN(2,1) + clhs60 + clhs61;
const double clhs63 =             DN(2,1)*clhs17;
const double clhs64 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs52;
const double clhs65 =             C(3,4)*DN(2,1);
const double clhs66 =             C(2,3)*DN(2,2) + clhs60 + clhs65;
const double clhs67 =             C(2,5)*DN(2,2);
const double clhs68 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs67;
const double clhs69 =             DN(2,2)*clhs17;
const double clhs70 =             DN(0,0)*N[2];
const double clhs71 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs72 =             C(0,3)*DN(3,0);
const double clhs73 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs72;
const double clhs74 =             C(0,5)*DN(3,0);
const double clhs75 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs74;
const double clhs76 =             N[3]*clhs31;
const double clhs77 =             DN(0,0)*DN(3,0);
const double clhs78 =             clhs10*clhs77 + clhs76;
const double clhs79 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs72;
const double clhs80 =             C(1,3)*DN(3,1);
const double clhs81 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs80;
const double clhs82 =             C(3,5)*DN(3,0);
const double clhs83 =             C(4,5)*DN(3,2);
const double clhs84 =             C(1,5)*DN(3,1) + clhs82 + clhs83;
const double clhs85 =             DN(3,1)*clhs17;
const double clhs86 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs74;
const double clhs87 =             C(3,4)*DN(3,1);
const double clhs88 =             C(2,3)*DN(3,2) + clhs82 + clhs87;
const double clhs89 =             C(2,5)*DN(3,2);
const double clhs90 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs89;
const double clhs91 =             DN(3,2)*clhs17;
const double clhs92 =             DN(0,0)*N[3];
const double clhs93 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs12;
const double clhs94 =             C(0,4)*DN(0,0) + clhs15 + clhs20;
const double clhs95 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs96 =             C(1,4)*DN(0,1);
const double clhs97 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs96;
const double clhs98 =             pow(DN(0,1), 2);
const double clhs99 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs96;
const double clhs100 =             C(2,4)*DN(0,2);
const double clhs101 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs100;
const double clhs102 =             DN(0,1)*clhs10;
const double clhs103 =             DN(0,2)*clhs102;
const double clhs104 =             DN(0,1)*N[0];
const double clhs105 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs36;
const double clhs106 =             C(0,4)*DN(1,0) + clhs39 + clhs43;
const double clhs107 =             DN(1,0)*clhs102;
const double clhs108 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs109 =             C(1,4)*DN(1,1);
const double clhs110 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs109;
const double clhs111 =             DN(0,1)*DN(1,1);
const double clhs112 =             clhs10*clhs111 + clhs32;
const double clhs113 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs109;
const double clhs114 =             C(2,4)*DN(1,2);
const double clhs115 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs114;
const double clhs116 =             DN(1,2)*clhs102;
const double clhs117 =             DN(0,1)*N[1];
const double clhs118 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs58;
const double clhs119 =             C(0,4)*DN(2,0) + clhs61 + clhs65;
const double clhs120 =             DN(2,0)*clhs102;
const double clhs121 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs122 =             C(1,4)*DN(2,1);
const double clhs123 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs122;
const double clhs124 =             DN(0,1)*DN(2,1);
const double clhs125 =             clhs10*clhs124 + clhs54;
const double clhs126 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs122;
const double clhs127 =             C(2,4)*DN(2,2);
const double clhs128 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs127;
const double clhs129 =             DN(2,2)*clhs102;
const double clhs130 =             DN(0,1)*N[2];
const double clhs131 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs80;
const double clhs132 =             C(0,4)*DN(3,0) + clhs83 + clhs87;
const double clhs133 =             DN(3,0)*clhs102;
const double clhs134 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs135 =             C(1,4)*DN(3,1);
const double clhs136 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs135;
const double clhs137 =             DN(0,1)*DN(3,1);
const double clhs138 =             clhs10*clhs137 + clhs76;
const double clhs139 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs135;
const double clhs140 =             C(2,4)*DN(3,2);
const double clhs141 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs140;
const double clhs142 =             DN(3,2)*clhs102;
const double clhs143 =             DN(0,1)*N[3];
const double clhs144 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs22;
const double clhs145 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs100;
const double clhs146 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs147 =             pow(DN(0,2), 2);
const double clhs148 =             DN(0,2)*N[0];
const double clhs149 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs45;
const double clhs150 =             DN(0,2)*clhs10;
const double clhs151 =             DN(1,0)*clhs150;
const double clhs152 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs114;
const double clhs153 =             DN(1,1)*clhs150;
const double clhs154 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs155 =             DN(0,2)*DN(1,2);
const double clhs156 =             clhs10*clhs155 + clhs32;
const double clhs157 =             DN(0,2)*N[1];
const double clhs158 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs67;
const double clhs159 =             DN(2,0)*clhs150;
const double clhs160 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs127;
const double clhs161 =             DN(2,1)*clhs150;
const double clhs162 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs163 =             DN(0,2)*DN(2,2);
const double clhs164 =             clhs10*clhs163 + clhs54;
const double clhs165 =             DN(0,2)*N[2];
const double clhs166 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs89;
const double clhs167 =             DN(3,0)*clhs150;
const double clhs168 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs140;
const double clhs169 =             DN(3,1)*clhs150;
const double clhs170 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs171 =             DN(0,2)*DN(3,2);
const double clhs172 =             clhs10*clhs171 + clhs76;
const double clhs173 =             DN(0,2)*N[3];
const double clhs174 =             1.0/clhs9;
const double clhs175 =             clhs0*clhs174;
const double clhs176 =             clhs175 + 1;
const double clhs177 =             DN(1,0)*N[0];
const double clhs178 =             DN(1,1)*N[0];
const double clhs179 =             DN(1,2)*N[0];
const double clhs180 =             clhs174*(clhs111 + clhs155 + clhs33);
const double clhs181 =             DN(2,0)*N[0];
const double clhs182 =             DN(2,1)*N[0];
const double clhs183 =             DN(2,2)*N[0];
const double clhs184 =             clhs174*(clhs124 + clhs163 + clhs55);
const double clhs185 =             DN(3,0)*N[0];
const double clhs186 =             DN(3,1)*N[0];
const double clhs187 =             DN(3,2)*N[0];
const double clhs188 =             clhs174*(clhs137 + clhs171 + clhs77);
const double clhs189 =             pow(N[1], 2)*clhs0;
const double clhs190 =             pow(DN(1,0), 2);
const double clhs191 =             DN(1,0)*clhs10;
const double clhs192 =             DN(1,1)*clhs191;
const double clhs193 =             DN(1,2)*clhs191;
const double clhs194 =             DN(1,0)*N[1];
const double clhs195 =             N[1]*clhs0;
const double clhs196 =             N[2]*clhs195;
const double clhs197 =             DN(1,0)*DN(2,0);
const double clhs198 =             clhs10*clhs197 + clhs196;
const double clhs199 =             DN(2,1)*clhs191;
const double clhs200 =             DN(2,2)*clhs191;
const double clhs201 =             DN(1,0)*N[2];
const double clhs202 =             N[3]*clhs195;
const double clhs203 =             DN(1,0)*DN(3,0);
const double clhs204 =             clhs10*clhs203 + clhs202;
const double clhs205 =             DN(3,1)*clhs191;
const double clhs206 =             DN(3,2)*clhs191;
const double clhs207 =             DN(1,0)*N[3];
const double clhs208 =             pow(DN(1,1), 2);
const double clhs209 =             DN(1,1)*clhs10;
const double clhs210 =             DN(1,2)*clhs209;
const double clhs211 =             DN(1,1)*N[1];
const double clhs212 =             DN(2,0)*clhs209;
const double clhs213 =             DN(1,1)*DN(2,1);
const double clhs214 =             clhs10*clhs213 + clhs196;
const double clhs215 =             DN(2,2)*clhs209;
const double clhs216 =             DN(1,1)*N[2];
const double clhs217 =             DN(3,0)*clhs209;
const double clhs218 =             DN(1,1)*DN(3,1);
const double clhs219 =             clhs10*clhs218 + clhs202;
const double clhs220 =             DN(3,2)*clhs209;
const double clhs221 =             DN(1,1)*N[3];
const double clhs222 =             pow(DN(1,2), 2);
const double clhs223 =             DN(1,2)*N[1];
const double clhs224 =             DN(1,2)*clhs10;
const double clhs225 =             DN(2,0)*clhs224;
const double clhs226 =             DN(2,1)*clhs224;
const double clhs227 =             DN(1,2)*DN(2,2);
const double clhs228 =             clhs10*clhs227 + clhs196;
const double clhs229 =             DN(1,2)*N[2];
const double clhs230 =             DN(3,0)*clhs224;
const double clhs231 =             DN(3,1)*clhs224;
const double clhs232 =             DN(1,2)*DN(3,2);
const double clhs233 =             clhs10*clhs232 + clhs202;
const double clhs234 =             DN(1,2)*N[3];
const double clhs235 =             DN(2,0)*N[1];
const double clhs236 =             DN(2,1)*N[1];
const double clhs237 =             DN(2,2)*N[1];
const double clhs238 =             clhs174*(clhs197 + clhs213 + clhs227);
const double clhs239 =             DN(3,0)*N[1];
const double clhs240 =             DN(3,1)*N[1];
const double clhs241 =             DN(3,2)*N[1];
const double clhs242 =             clhs174*(clhs203 + clhs218 + clhs232);
const double clhs243 =             pow(N[2], 2)*clhs0;
const double clhs244 =             pow(DN(2,0), 2);
const double clhs245 =             DN(2,0)*clhs10;
const double clhs246 =             DN(2,1)*clhs245;
const double clhs247 =             DN(2,2)*clhs245;
const double clhs248 =             DN(2,0)*N[2];
const double clhs249 =             N[2]*N[3]*clhs0;
const double clhs250 =             DN(2,0)*DN(3,0);
const double clhs251 =             clhs10*clhs250 + clhs249;
const double clhs252 =             DN(3,1)*clhs245;
const double clhs253 =             DN(3,2)*clhs245;
const double clhs254 =             DN(2,0)*N[3];
const double clhs255 =             pow(DN(2,1), 2);
const double clhs256 =             DN(2,1)*clhs10;
const double clhs257 =             DN(2,2)*clhs256;
const double clhs258 =             DN(2,1)*N[2];
const double clhs259 =             DN(3,0)*clhs256;
const double clhs260 =             DN(2,1)*DN(3,1);
const double clhs261 =             clhs10*clhs260 + clhs249;
const double clhs262 =             DN(3,2)*clhs256;
const double clhs263 =             DN(2,1)*N[3];
const double clhs264 =             pow(DN(2,2), 2);
const double clhs265 =             DN(2,2)*N[2];
const double clhs266 =             DN(2,2)*clhs10;
const double clhs267 =             DN(3,0)*clhs266;
const double clhs268 =             DN(3,1)*clhs266;
const double clhs269 =             DN(2,2)*DN(3,2);
const double clhs270 =             clhs10*clhs269 + clhs249;
const double clhs271 =             DN(2,2)*N[3];
const double clhs272 =             DN(3,0)*N[2];
const double clhs273 =             DN(3,1)*N[2];
const double clhs274 =             DN(3,2)*N[2];
const double clhs275 =             clhs174*(clhs250 + clhs260 + clhs269);
const double clhs276 =             pow(N[3], 2)*clhs0;
const double clhs277 =             pow(DN(3,0), 2);
const double clhs278 =             DN(3,0)*clhs10;
const double clhs279 =             DN(3,1)*clhs278;
const double clhs280 =             DN(3,2)*clhs278;
const double clhs281 =             DN(3,0)*N[3];
const double clhs282 =             pow(DN(3,1), 2);
const double clhs283 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs284 =             DN(3,1)*N[3];
const double clhs285 =             pow(DN(3,2), 2);
const double clhs286 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs2 + DN(0,1)*clhs4 + DN(0,2)*clhs6 + clhs1 + clhs10*clhs7;
            lhs(0,1)=DN(0,0)*clhs11 + DN(0,1)*clhs13 + DN(0,2)*clhs16 + clhs18;
            lhs(0,2)=DN(0,0)*clhs19 + DN(0,1)*clhs21 + DN(0,2)*clhs23 + clhs24;
            lhs(0,3)=-clhs25;
            lhs(0,4)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + DN(0,2)*clhs30 + clhs34;
            lhs(0,5)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + DN(0,2)*clhs40 + clhs41;
            lhs(0,6)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + DN(0,2)*clhs46 + clhs47;
            lhs(0,7)=-clhs48;
            lhs(0,8)=DN(0,0)*clhs49 + DN(0,1)*clhs51 + DN(0,2)*clhs53 + clhs56;
            lhs(0,9)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + DN(0,2)*clhs62 + clhs63;
            lhs(0,10)=DN(0,0)*clhs64 + DN(0,1)*clhs66 + DN(0,2)*clhs68 + clhs69;
            lhs(0,11)=-clhs70;
            lhs(0,12)=DN(0,0)*clhs71 + DN(0,1)*clhs73 + DN(0,2)*clhs75 + clhs78;
            lhs(0,13)=DN(0,0)*clhs79 + DN(0,1)*clhs81 + DN(0,2)*clhs84 + clhs85;
            lhs(0,14)=DN(0,0)*clhs86 + DN(0,1)*clhs88 + DN(0,2)*clhs90 + clhs91;
            lhs(0,15)=-clhs92;
            lhs(1,0)=DN(0,0)*clhs4 + DN(0,1)*clhs93 + DN(0,2)*clhs94 + clhs18;
            lhs(1,1)=DN(0,0)*clhs13 + DN(0,1)*clhs95 + DN(0,2)*clhs97 + clhs1 + clhs10*clhs98;
            lhs(1,2)=DN(0,0)*clhs21 + DN(0,1)*clhs99 + DN(0,2)*clhs101 + clhs103;
            lhs(1,3)=-clhs104;
            lhs(1,4)=DN(0,0)*clhs28 + DN(0,1)*clhs105 + DN(0,2)*clhs106 + clhs107;
            lhs(1,5)=DN(0,0)*clhs37 + DN(0,1)*clhs108 + DN(0,2)*clhs110 + clhs112;
            lhs(1,6)=DN(0,0)*clhs44 + DN(0,1)*clhs113 + DN(0,2)*clhs115 + clhs116;
            lhs(1,7)=-clhs117;
            lhs(1,8)=DN(0,0)*clhs51 + DN(0,1)*clhs118 + DN(0,2)*clhs119 + clhs120;
            lhs(1,9)=DN(0,0)*clhs59 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs125;
            lhs(1,10)=DN(0,0)*clhs66 + DN(0,1)*clhs126 + DN(0,2)*clhs128 + clhs129;
            lhs(1,11)=-clhs130;
            lhs(1,12)=DN(0,0)*clhs73 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs133;
            lhs(1,13)=DN(0,0)*clhs81 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs138;
            lhs(1,14)=DN(0,0)*clhs88 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142;
            lhs(1,15)=-clhs143;
            lhs(2,0)=DN(0,0)*clhs6 + DN(0,1)*clhs94 + DN(0,2)*clhs144 + clhs24;
            lhs(2,1)=DN(0,0)*clhs16 + DN(0,1)*clhs97 + DN(0,2)*clhs145 + clhs103;
            lhs(2,2)=DN(0,0)*clhs23 + DN(0,1)*clhs101 + DN(0,2)*clhs146 + clhs1 + clhs10*clhs147;
            lhs(2,3)=-clhs148;
            lhs(2,4)=DN(0,0)*clhs30 + DN(0,1)*clhs106 + DN(0,2)*clhs149 + clhs151;
            lhs(2,5)=DN(0,0)*clhs40 + DN(0,1)*clhs110 + DN(0,2)*clhs152 + clhs153;
            lhs(2,6)=DN(0,0)*clhs46 + DN(0,1)*clhs115 + DN(0,2)*clhs154 + clhs156;
            lhs(2,7)=-clhs157;
            lhs(2,8)=DN(0,0)*clhs53 + DN(0,1)*clhs119 + DN(0,2)*clhs158 + clhs159;
            lhs(2,9)=DN(0,0)*clhs62 + DN(0,1)*clhs123 + DN(0,2)*clhs160 + clhs161;
            lhs(2,10)=DN(0,0)*clhs68 + DN(0,1)*clhs128 + DN(0,2)*clhs162 + clhs164;
            lhs(2,11)=-clhs165;
            lhs(2,12)=DN(0,0)*clhs75 + DN(0,1)*clhs132 + DN(0,2)*clhs166 + clhs167;
            lhs(2,13)=DN(0,0)*clhs84 + DN(0,1)*clhs136 + DN(0,2)*clhs168 + clhs169;
            lhs(2,14)=DN(0,0)*clhs90 + DN(0,1)*clhs141 + DN(0,2)*clhs170 + clhs172;
            lhs(2,15)=-clhs173;
            lhs(3,0)=clhs176*clhs25;
            lhs(3,1)=clhs104*clhs176;
            lhs(3,2)=clhs148*clhs176;
            lhs(3,3)=clhs174*(clhs147 + clhs7 + clhs98);
            lhs(3,4)=clhs175*clhs48 + clhs177;
            lhs(3,5)=clhs117*clhs175 + clhs178;
            lhs(3,6)=clhs157*clhs175 + clhs179;
            lhs(3,7)=clhs180;
            lhs(3,8)=clhs175*clhs70 + clhs181;
            lhs(3,9)=clhs130*clhs175 + clhs182;
            lhs(3,10)=clhs165*clhs175 + clhs183;
            lhs(3,11)=clhs184;
            lhs(3,12)=clhs175*clhs92 + clhs185;
            lhs(3,13)=clhs143*clhs175 + clhs186;
            lhs(3,14)=clhs173*clhs175 + clhs187;
            lhs(3,15)=clhs188;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs4 + DN(1,2)*clhs6 + clhs34;
            lhs(4,1)=DN(1,0)*clhs11 + DN(1,1)*clhs13 + DN(1,2)*clhs16 + clhs107;
            lhs(4,2)=DN(1,0)*clhs19 + DN(1,1)*clhs21 + DN(1,2)*clhs23 + clhs151;
            lhs(4,3)=-clhs177;
            lhs(4,4)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + DN(1,2)*clhs30 + clhs10*clhs190 + clhs189;
            lhs(4,5)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + DN(1,2)*clhs40 + clhs192;
            lhs(4,6)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + DN(1,2)*clhs46 + clhs193;
            lhs(4,7)=-clhs194;
            lhs(4,8)=DN(1,0)*clhs49 + DN(1,1)*clhs51 + DN(1,2)*clhs53 + clhs198;
            lhs(4,9)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + DN(1,2)*clhs62 + clhs199;
            lhs(4,10)=DN(1,0)*clhs64 + DN(1,1)*clhs66 + DN(1,2)*clhs68 + clhs200;
            lhs(4,11)=-clhs201;
            lhs(4,12)=DN(1,0)*clhs71 + DN(1,1)*clhs73 + DN(1,2)*clhs75 + clhs204;
            lhs(4,13)=DN(1,0)*clhs79 + DN(1,1)*clhs81 + DN(1,2)*clhs84 + clhs205;
            lhs(4,14)=DN(1,0)*clhs86 + DN(1,1)*clhs88 + DN(1,2)*clhs90 + clhs206;
            lhs(4,15)=-clhs207;
            lhs(5,0)=DN(1,0)*clhs4 + DN(1,1)*clhs93 + DN(1,2)*clhs94 + clhs41;
            lhs(5,1)=DN(1,0)*clhs13 + DN(1,1)*clhs95 + DN(1,2)*clhs97 + clhs112;
            lhs(5,2)=DN(1,0)*clhs21 + DN(1,1)*clhs99 + DN(1,2)*clhs101 + clhs153;
            lhs(5,3)=-clhs178;
            lhs(5,4)=DN(1,0)*clhs28 + DN(1,1)*clhs105 + DN(1,2)*clhs106 + clhs192;
            lhs(5,5)=DN(1,0)*clhs37 + DN(1,1)*clhs108 + DN(1,2)*clhs110 + clhs10*clhs208 + clhs189;
            lhs(5,6)=DN(1,0)*clhs44 + DN(1,1)*clhs113 + DN(1,2)*clhs115 + clhs210;
            lhs(5,7)=-clhs211;
            lhs(5,8)=DN(1,0)*clhs51 + DN(1,1)*clhs118 + DN(1,2)*clhs119 + clhs212;
            lhs(5,9)=DN(1,0)*clhs59 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs214;
            lhs(5,10)=DN(1,0)*clhs66 + DN(1,1)*clhs126 + DN(1,2)*clhs128 + clhs215;
            lhs(5,11)=-clhs216;
            lhs(5,12)=DN(1,0)*clhs73 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs217;
            lhs(5,13)=DN(1,0)*clhs81 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs219;
            lhs(5,14)=DN(1,0)*clhs88 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs220;
            lhs(5,15)=-clhs221;
            lhs(6,0)=DN(1,0)*clhs6 + DN(1,1)*clhs94 + DN(1,2)*clhs144 + clhs47;
            lhs(6,1)=DN(1,0)*clhs16 + DN(1,1)*clhs97 + DN(1,2)*clhs145 + clhs116;
            lhs(6,2)=DN(1,0)*clhs23 + DN(1,1)*clhs101 + DN(1,2)*clhs146 + clhs156;
            lhs(6,3)=-clhs179;
            lhs(6,4)=DN(1,0)*clhs30 + DN(1,1)*clhs106 + DN(1,2)*clhs149 + clhs193;
            lhs(6,5)=DN(1,0)*clhs40 + DN(1,1)*clhs110 + DN(1,2)*clhs152 + clhs210;
            lhs(6,6)=DN(1,0)*clhs46 + DN(1,1)*clhs115 + DN(1,2)*clhs154 + clhs10*clhs222 + clhs189;
            lhs(6,7)=-clhs223;
            lhs(6,8)=DN(1,0)*clhs53 + DN(1,1)*clhs119 + DN(1,2)*clhs158 + clhs225;
            lhs(6,9)=DN(1,0)*clhs62 + DN(1,1)*clhs123 + DN(1,2)*clhs160 + clhs226;
            lhs(6,10)=DN(1,0)*clhs68 + DN(1,1)*clhs128 + DN(1,2)*clhs162 + clhs228;
            lhs(6,11)=-clhs229;
            lhs(6,12)=DN(1,0)*clhs75 + DN(1,1)*clhs132 + DN(1,2)*clhs166 + clhs230;
            lhs(6,13)=DN(1,0)*clhs84 + DN(1,1)*clhs136 + DN(1,2)*clhs168 + clhs231;
            lhs(6,14)=DN(1,0)*clhs90 + DN(1,1)*clhs141 + DN(1,2)*clhs170 + clhs233;
            lhs(6,15)=-clhs234;
            lhs(7,0)=clhs175*clhs177 + clhs48;
            lhs(7,1)=clhs117 + clhs175*clhs178;
            lhs(7,2)=clhs157 + clhs175*clhs179;
            lhs(7,3)=clhs180;
            lhs(7,4)=clhs176*clhs194;
            lhs(7,5)=clhs176*clhs211;
            lhs(7,6)=clhs176*clhs223;
            lhs(7,7)=clhs174*(clhs190 + clhs208 + clhs222);
            lhs(7,8)=clhs175*clhs201 + clhs235;
            lhs(7,9)=clhs175*clhs216 + clhs236;
            lhs(7,10)=clhs175*clhs229 + clhs237;
            lhs(7,11)=clhs238;
            lhs(7,12)=clhs175*clhs207 + clhs239;
            lhs(7,13)=clhs175*clhs221 + clhs240;
            lhs(7,14)=clhs175*clhs234 + clhs241;
            lhs(7,15)=clhs242;
            lhs(8,0)=DN(2,0)*clhs2 + DN(2,1)*clhs4 + DN(2,2)*clhs6 + clhs56;
            lhs(8,1)=DN(2,0)*clhs11 + DN(2,1)*clhs13 + DN(2,2)*clhs16 + clhs120;
            lhs(8,2)=DN(2,0)*clhs19 + DN(2,1)*clhs21 + DN(2,2)*clhs23 + clhs159;
            lhs(8,3)=-clhs181;
            lhs(8,4)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + DN(2,2)*clhs30 + clhs198;
            lhs(8,5)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + DN(2,2)*clhs40 + clhs212;
            lhs(8,6)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + DN(2,2)*clhs46 + clhs225;
            lhs(8,7)=-clhs235;
            lhs(8,8)=DN(2,0)*clhs49 + DN(2,1)*clhs51 + DN(2,2)*clhs53 + clhs10*clhs244 + clhs243;
            lhs(8,9)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + DN(2,2)*clhs62 + clhs246;
            lhs(8,10)=DN(2,0)*clhs64 + DN(2,1)*clhs66 + DN(2,2)*clhs68 + clhs247;
            lhs(8,11)=-clhs248;
            lhs(8,12)=DN(2,0)*clhs71 + DN(2,1)*clhs73 + DN(2,2)*clhs75 + clhs251;
            lhs(8,13)=DN(2,0)*clhs79 + DN(2,1)*clhs81 + DN(2,2)*clhs84 + clhs252;
            lhs(8,14)=DN(2,0)*clhs86 + DN(2,1)*clhs88 + DN(2,2)*clhs90 + clhs253;
            lhs(8,15)=-clhs254;
            lhs(9,0)=DN(2,0)*clhs4 + DN(2,1)*clhs93 + DN(2,2)*clhs94 + clhs63;
            lhs(9,1)=DN(2,0)*clhs13 + DN(2,1)*clhs95 + DN(2,2)*clhs97 + clhs125;
            lhs(9,2)=DN(2,0)*clhs21 + DN(2,1)*clhs99 + DN(2,2)*clhs101 + clhs161;
            lhs(9,3)=-clhs182;
            lhs(9,4)=DN(2,0)*clhs28 + DN(2,1)*clhs105 + DN(2,2)*clhs106 + clhs199;
            lhs(9,5)=DN(2,0)*clhs37 + DN(2,1)*clhs108 + DN(2,2)*clhs110 + clhs214;
            lhs(9,6)=DN(2,0)*clhs44 + DN(2,1)*clhs113 + DN(2,2)*clhs115 + clhs226;
            lhs(9,7)=-clhs236;
            lhs(9,8)=DN(2,0)*clhs51 + DN(2,1)*clhs118 + DN(2,2)*clhs119 + clhs246;
            lhs(9,9)=DN(2,0)*clhs59 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs10*clhs255 + clhs243;
            lhs(9,10)=DN(2,0)*clhs66 + DN(2,1)*clhs126 + DN(2,2)*clhs128 + clhs257;
            lhs(9,11)=-clhs258;
            lhs(9,12)=DN(2,0)*clhs73 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs259;
            lhs(9,13)=DN(2,0)*clhs81 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs261;
            lhs(9,14)=DN(2,0)*clhs88 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs262;
            lhs(9,15)=-clhs263;
            lhs(10,0)=DN(2,0)*clhs6 + DN(2,1)*clhs94 + DN(2,2)*clhs144 + clhs69;
            lhs(10,1)=DN(2,0)*clhs16 + DN(2,1)*clhs97 + DN(2,2)*clhs145 + clhs129;
            lhs(10,2)=DN(2,0)*clhs23 + DN(2,1)*clhs101 + DN(2,2)*clhs146 + clhs164;
            lhs(10,3)=-clhs183;
            lhs(10,4)=DN(2,0)*clhs30 + DN(2,1)*clhs106 + DN(2,2)*clhs149 + clhs200;
            lhs(10,5)=DN(2,0)*clhs40 + DN(2,1)*clhs110 + DN(2,2)*clhs152 + clhs215;
            lhs(10,6)=DN(2,0)*clhs46 + DN(2,1)*clhs115 + DN(2,2)*clhs154 + clhs228;
            lhs(10,7)=-clhs237;
            lhs(10,8)=DN(2,0)*clhs53 + DN(2,1)*clhs119 + DN(2,2)*clhs158 + clhs247;
            lhs(10,9)=DN(2,0)*clhs62 + DN(2,1)*clhs123 + DN(2,2)*clhs160 + clhs257;
            lhs(10,10)=DN(2,0)*clhs68 + DN(2,1)*clhs128 + DN(2,2)*clhs162 + clhs10*clhs264 + clhs243;
            lhs(10,11)=-clhs265;
            lhs(10,12)=DN(2,0)*clhs75 + DN(2,1)*clhs132 + DN(2,2)*clhs166 + clhs267;
            lhs(10,13)=DN(2,0)*clhs84 + DN(2,1)*clhs136 + DN(2,2)*clhs168 + clhs268;
            lhs(10,14)=DN(2,0)*clhs90 + DN(2,1)*clhs141 + DN(2,2)*clhs170 + clhs270;
            lhs(10,15)=-clhs271;
            lhs(11,0)=clhs175*clhs181 + clhs70;
            lhs(11,1)=clhs130 + clhs175*clhs182;
            lhs(11,2)=clhs165 + clhs175*clhs183;
            lhs(11,3)=clhs184;
            lhs(11,4)=clhs175*clhs235 + clhs201;
            lhs(11,5)=clhs175*clhs236 + clhs216;
            lhs(11,6)=clhs175*clhs237 + clhs229;
            lhs(11,7)=clhs238;
            lhs(11,8)=clhs176*clhs248;
            lhs(11,9)=clhs176*clhs258;
            lhs(11,10)=clhs176*clhs265;
            lhs(11,11)=clhs174*(clhs244 + clhs255 + clhs264);
            lhs(11,12)=clhs175*clhs254 + clhs272;
            lhs(11,13)=clhs175*clhs263 + clhs273;
            lhs(11,14)=clhs175*clhs271 + clhs274;
            lhs(11,15)=clhs275;
            lhs(12,0)=DN(3,0)*clhs2 + DN(3,1)*clhs4 + DN(3,2)*clhs6 + clhs78;
            lhs(12,1)=DN(3,0)*clhs11 + DN(3,1)*clhs13 + DN(3,2)*clhs16 + clhs133;
            lhs(12,2)=DN(3,0)*clhs19 + DN(3,1)*clhs21 + DN(3,2)*clhs23 + clhs167;
            lhs(12,3)=-clhs185;
            lhs(12,4)=DN(3,0)*clhs26 + DN(3,1)*clhs28 + DN(3,2)*clhs30 + clhs204;
            lhs(12,5)=DN(3,0)*clhs35 + DN(3,1)*clhs37 + DN(3,2)*clhs40 + clhs217;
            lhs(12,6)=DN(3,0)*clhs42 + DN(3,1)*clhs44 + DN(3,2)*clhs46 + clhs230;
            lhs(12,7)=-clhs239;
            lhs(12,8)=DN(3,0)*clhs49 + DN(3,1)*clhs51 + DN(3,2)*clhs53 + clhs251;
            lhs(12,9)=DN(3,0)*clhs57 + DN(3,1)*clhs59 + DN(3,2)*clhs62 + clhs259;
            lhs(12,10)=DN(3,0)*clhs64 + DN(3,1)*clhs66 + DN(3,2)*clhs68 + clhs267;
            lhs(12,11)=-clhs272;
            lhs(12,12)=DN(3,0)*clhs71 + DN(3,1)*clhs73 + DN(3,2)*clhs75 + clhs10*clhs277 + clhs276;
            lhs(12,13)=DN(3,0)*clhs79 + DN(3,1)*clhs81 + DN(3,2)*clhs84 + clhs279;
            lhs(12,14)=DN(3,0)*clhs86 + DN(3,1)*clhs88 + DN(3,2)*clhs90 + clhs280;
            lhs(12,15)=-clhs281;
            lhs(13,0)=DN(3,0)*clhs4 + DN(3,1)*clhs93 + DN(3,2)*clhs94 + clhs85;
            lhs(13,1)=DN(3,0)*clhs13 + DN(3,1)*clhs95 + DN(3,2)*clhs97 + clhs138;
            lhs(13,2)=DN(3,0)*clhs21 + DN(3,1)*clhs99 + DN(3,2)*clhs101 + clhs169;
            lhs(13,3)=-clhs186;
            lhs(13,4)=DN(3,0)*clhs28 + DN(3,1)*clhs105 + DN(3,2)*clhs106 + clhs205;
            lhs(13,5)=DN(3,0)*clhs37 + DN(3,1)*clhs108 + DN(3,2)*clhs110 + clhs219;
            lhs(13,6)=DN(3,0)*clhs44 + DN(3,1)*clhs113 + DN(3,2)*clhs115 + clhs231;
            lhs(13,7)=-clhs240;
            lhs(13,8)=DN(3,0)*clhs51 + DN(3,1)*clhs118 + DN(3,2)*clhs119 + clhs252;
            lhs(13,9)=DN(3,0)*clhs59 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs261;
            lhs(13,10)=DN(3,0)*clhs66 + DN(3,1)*clhs126 + DN(3,2)*clhs128 + clhs268;
            lhs(13,11)=-clhs273;
            lhs(13,12)=DN(3,0)*clhs73 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs279;
            lhs(13,13)=DN(3,0)*clhs81 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs10*clhs282 + clhs276;
            lhs(13,14)=DN(3,0)*clhs88 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs283;
            lhs(13,15)=-clhs284;
            lhs(14,0)=DN(3,0)*clhs6 + DN(3,1)*clhs94 + DN(3,2)*clhs144 + clhs91;
            lhs(14,1)=DN(3,0)*clhs16 + DN(3,1)*clhs97 + DN(3,2)*clhs145 + clhs142;
            lhs(14,2)=DN(3,0)*clhs23 + DN(3,1)*clhs101 + DN(3,2)*clhs146 + clhs172;
            lhs(14,3)=-clhs187;
            lhs(14,4)=DN(3,0)*clhs30 + DN(3,1)*clhs106 + DN(3,2)*clhs149 + clhs206;
            lhs(14,5)=DN(3,0)*clhs40 + DN(3,1)*clhs110 + DN(3,2)*clhs152 + clhs220;
            lhs(14,6)=DN(3,0)*clhs46 + DN(3,1)*clhs115 + DN(3,2)*clhs154 + clhs233;
            lhs(14,7)=-clhs241;
            lhs(14,8)=DN(3,0)*clhs53 + DN(3,1)*clhs119 + DN(3,2)*clhs158 + clhs253;
            lhs(14,9)=DN(3,0)*clhs62 + DN(3,1)*clhs123 + DN(3,2)*clhs160 + clhs262;
            lhs(14,10)=DN(3,0)*clhs68 + DN(3,1)*clhs128 + DN(3,2)*clhs162 + clhs270;
            lhs(14,11)=-clhs274;
            lhs(14,12)=DN(3,0)*clhs75 + DN(3,1)*clhs132 + DN(3,2)*clhs166 + clhs280;
            lhs(14,13)=DN(3,0)*clhs84 + DN(3,1)*clhs136 + DN(3,2)*clhs168 + clhs283;
            lhs(14,14)=DN(3,0)*clhs90 + DN(3,1)*clhs141 + DN(3,2)*clhs170 + clhs10*clhs285 + clhs276;
            lhs(14,15)=-clhs286;
            lhs(15,0)=clhs175*clhs185 + clhs92;
            lhs(15,1)=clhs143 + clhs175*clhs186;
            lhs(15,2)=clhs173 + clhs175*clhs187;
            lhs(15,3)=clhs188;
            lhs(15,4)=clhs175*clhs239 + clhs207;
            lhs(15,5)=clhs175*clhs240 + clhs221;
            lhs(15,6)=clhs175*clhs241 + clhs234;
            lhs(15,7)=clhs242;
            lhs(15,8)=clhs175*clhs272 + clhs254;
            lhs(15,9)=clhs175*clhs273 + clhs263;
            lhs(15,10)=clhs175*clhs274 + clhs271;
            lhs(15,11)=clhs275;
            lhs(15,12)=clhs176*clhs281;
            lhs(15,13)=clhs176*clhs284;
            lhs(15,14)=clhs176*clhs286;
            lhs(15,15)=clhs174*(clhs277 + clhs282 + clhs285);


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

    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs4 =             pow(h, 2);
const double crhs5 =             dyn_tau*rho/dt + mu*stab_c1/crhs4;
const double crhs6 =             1.0*crhs3*crhs4*crhs5/stab_c1;
const double crhs7 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs8 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs9 =             1.0/crhs5;
const double crhs10 =             crhs9*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2);
const double crhs11 =             crhs9*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs7 + crhs8);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs6 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs6 - DN(0,1)*stress[1] + N[0]*crhs7 - N[0]*crhs8;
            rhs[2]=-DN(0,0)*crhs10 - DN(0,1)*crhs11 - N[0]*crhs3;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs6 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs6 - DN(1,1)*stress[1] + N[1]*crhs7 - N[1]*crhs8;
            rhs[5]=-DN(1,0)*crhs10 - DN(1,1)*crhs11 - N[1]*crhs3;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs6 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs6 - DN(2,1)*stress[1] + N[2]*crhs7 - N[2]*crhs8;
            rhs[8]=-DN(2,0)*crhs10 - DN(2,1)*crhs11 - N[2]*crhs3;


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

    auto& rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs4 =             pow(h, 2);
const double crhs5 =             dyn_tau*rho/dt + mu*stab_c1/crhs4;
const double crhs6 =             1.0*crhs3*crhs4*crhs5/stab_c1;
const double crhs7 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs8 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs9 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs10 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs11 =             1.0/crhs5;
const double crhs12 =             crhs11*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs2);
const double crhs13 =             crhs11*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs7 + crhs8);
const double crhs14 =             crhs11*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + crhs10 - crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs6 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs2;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs6 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs7 - N[0]*crhs8;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs6 - DN(0,2)*stress[2] - N[0]*crhs10 + N[0]*crhs9;
            rhs[3]=-DN(0,0)*crhs12 - DN(0,1)*crhs13 - DN(0,2)*crhs14 - N[0]*crhs3;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs6 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs2;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs6 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs7 - N[1]*crhs8;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs6 - DN(1,2)*stress[2] - N[1]*crhs10 + N[1]*crhs9;
            rhs[7]=-DN(1,0)*crhs12 - DN(1,1)*crhs13 - DN(1,2)*crhs14 - N[1]*crhs3;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs6 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs2;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs6 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs7 - N[2]*crhs8;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs6 - DN(2,2)*stress[2] - N[2]*crhs10 + N[2]*crhs9;
            rhs[11]=-DN(2,0)*crhs12 - DN(2,1)*crhs13 - DN(2,2)*crhs14 - N[2]*crhs3;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs6 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs2;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs6 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs7 - N[3]*crhs8;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs6 - DN(3,2)*stress[2] - N[3]*crhs10 + N[3]*crhs9;
            rhs[15]=-DN(3,0)*crhs12 - DN(3,1)*crhs13 - DN(3,2)*crhs14 - N[3]*crhs3;


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
template class SymbolicStokes< SymbolicStokesData<3,6> >;
template class SymbolicStokes< SymbolicStokesData<3,8> >;

}
