//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//

#include "stokes.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
Stokes<TElementData>::Stokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
Stokes<TElementData>::Stokes(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
Stokes<TElementData>::Stokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
Stokes<TElementData>::Stokes(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
Stokes<TElementData>::~Stokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer Stokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<Stokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer Stokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<Stokes>(NewId, pGeom, pProperties);
}

template <class TElementData>
void Stokes<TElementData>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (TElementData::ElementManagesTimeIntegration){
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        //Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();
        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
        }

    } else{
        KRATOS_ERROR << "Stokes is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void Stokes<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int Stokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    KRATOS_CHECK_VARIABLE_KEY( DIVERGENCE );

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
std::string Stokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "Stokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void Stokes<TElementData>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr){
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

template <class TElementData>
void Stokes<TElementData>::Calculate(
    const Variable<Vector >& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo )
{
    noalias( rOutput ) = ZeroVector( StrainSize );

    if (rVariable == FLUID_STRESS) {

        // creating a new data container that goes out of scope after the function is left
        TElementData dataLocal;

        // transferring the velocity (among other variables)
        dataLocal.Initialize(*this, rCurrentProcessInfo);

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;

        // computing DN_DX values for the strain rate
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        double sumOfGaussWeights = 0.0;

        for (unsigned int g = 0; g < number_of_gauss_points; g++){

            UpdateIntegrationPointData(dataLocal, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);

            const Vector gauss_point_contribution = dataLocal.ShearStress;

            noalias( rOutput ) += gauss_point_contribution * gauss_weights[g];
            sumOfGaussWeights += gauss_weights[g];
        }

        for (unsigned int i = 0; i < StrainSize; i++){
            rOutput[i] = ( 1.0 / sumOfGaussWeights ) * rOutput[i];
        }

    } else {

        Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);

    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void Stokes<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void Stokes<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void Stokes<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void Stokes<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
}

template <>
void Stokes<SymbolicNavierStokesData<2,3>>::ComputeGaussPointLHSContribution(
    SymbolicNavierStokesData<2,3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             bdf0*rho;
const double clhs2 =             pow(N[0], 2)*clhs1;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs4 =             C(0,2)*DN(0,0);
const double clhs5 =             C(2,2)*DN(0,1) + clhs4;
const double clhs6 =             DN(0,0)*mu;
const double clhs7 =             DN(0,1)*clhs6;
const double clhs8 =             C(0,1)*DN(0,1) + clhs4;
const double clhs9 =             C(1,2)*DN(0,1);
const double clhs10 =             C(2,2)*DN(0,0) + clhs9;
const double clhs11 =             DN(0,0)*N[0];
const double clhs12 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs13 =             C(0,2)*DN(1,0);
const double clhs14 =             C(2,2)*DN(1,1) + clhs13;
const double clhs15 =             DN(0,0)*DN(1,0);
const double clhs16 =             N[0]*clhs1;
const double clhs17 =             N[1]*clhs16;
const double clhs18 =             clhs15*mu + clhs17;
const double clhs19 =             DN(1,1)*clhs6;
const double clhs20 =             C(0,1)*DN(1,1) + clhs13;
const double clhs21 =             C(1,2)*DN(1,1);
const double clhs22 =             C(2,2)*DN(1,0) + clhs21;
const double clhs23 =             DN(0,0)*N[1];
const double clhs24 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs25 =             C(0,2)*DN(2,0);
const double clhs26 =             C(2,2)*DN(2,1) + clhs25;
const double clhs27 =             DN(0,0)*DN(2,0);
const double clhs28 =             N[2]*clhs16;
const double clhs29 =             clhs27*mu + clhs28;
const double clhs30 =             DN(2,1)*clhs6;
const double clhs31 =             C(0,1)*DN(2,1) + clhs25;
const double clhs32 =             C(1,2)*DN(2,1);
const double clhs33 =             C(2,2)*DN(2,0) + clhs32;
const double clhs34 =             DN(0,0)*N[2];
const double clhs35 =             C(0,1)*DN(0,0) + clhs9;
const double clhs36 =             pow(DN(0,1), 2);
const double clhs37 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs38 =             DN(0,1)*N[0];
const double clhs39 =             DN(0,1)*mu;
const double clhs40 =             DN(1,0)*clhs39;
const double clhs41 =             C(0,1)*DN(1,0) + clhs21;
const double clhs42 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs43 =             DN(0,1)*DN(1,1);
const double clhs44 =             clhs17 + clhs43*mu;
const double clhs45 =             DN(0,1)*N[1];
const double clhs46 =             DN(2,0)*clhs39;
const double clhs47 =             C(0,1)*DN(2,0) + clhs32;
const double clhs48 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs49 =             DN(0,1)*DN(2,1);
const double clhs50 =             clhs28 + clhs49*mu;
const double clhs51 =             DN(0,1)*N[2];
const double clhs52 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs53 =             clhs1*clhs52;
const double clhs54 =             clhs53 + 1;
const double clhs55 =             DN(1,0)*N[0];
const double clhs56 =             DN(1,1)*N[0];
const double clhs57 =             clhs52*(clhs15 + clhs43);
const double clhs58 =             DN(2,0)*N[0];
const double clhs59 =             DN(2,1)*N[0];
const double clhs60 =             clhs52*(clhs27 + clhs49);
const double clhs61 =             pow(DN(1,0), 2);
const double clhs62 =             pow(N[1], 2)*clhs1;
const double clhs63 =             DN(1,0)*mu;
const double clhs64 =             DN(1,1)*clhs63;
const double clhs65 =             DN(1,0)*N[1];
const double clhs66 =             DN(1,0)*DN(2,0);
const double clhs67 =             N[1]*N[2]*clhs1;
const double clhs68 =             clhs66*mu + clhs67;
const double clhs69 =             DN(2,1)*clhs63;
const double clhs70 =             DN(1,0)*N[2];
const double clhs71 =             pow(DN(1,1), 2);
const double clhs72 =             DN(1,1)*N[1];
const double clhs73 =             DN(2,0)*mu;
const double clhs74 =             DN(1,1)*clhs73;
const double clhs75 =             DN(1,1)*DN(2,1);
const double clhs76 =             clhs67 + clhs75*mu;
const double clhs77 =             DN(1,1)*N[2];
const double clhs78 =             DN(2,0)*N[1];
const double clhs79 =             DN(2,1)*N[1];
const double clhs80 =             clhs52*(clhs66 + clhs75);
const double clhs81 =             pow(DN(2,0), 2);
const double clhs82 =             pow(N[2], 2)*clhs1;
const double clhs83 =             DN(2,1)*clhs73;
const double clhs84 =             DN(2,0)*N[2];
const double clhs85 =             pow(DN(2,1), 2);
const double clhs86 =             DN(2,1)*N[2];
            lhs(0,0)=DN(0,0)*clhs3 + DN(0,1)*clhs5 + clhs0*mu + clhs2;
            lhs(0,1)=DN(0,0)*clhs8 + DN(0,1)*clhs10 + clhs7;
            lhs(0,2)=-clhs11;
            lhs(0,3)=DN(0,0)*clhs12 + DN(0,1)*clhs14 + clhs18;
            lhs(0,4)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs19;
            lhs(0,5)=-clhs23;
            lhs(0,6)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs29;
            lhs(0,7)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + clhs30;
            lhs(0,8)=-clhs34;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs35 + clhs7;
            lhs(1,1)=DN(0,0)*clhs10 + DN(0,1)*clhs37 + clhs2 + clhs36*mu;
            lhs(1,2)=-clhs38;
            lhs(1,3)=DN(0,0)*clhs14 + DN(0,1)*clhs41 + clhs40;
            lhs(1,4)=DN(0,0)*clhs22 + DN(0,1)*clhs42 + clhs44;
            lhs(1,5)=-clhs45;
            lhs(1,6)=DN(0,0)*clhs26 + DN(0,1)*clhs47 + clhs46;
            lhs(1,7)=DN(0,0)*clhs33 + DN(0,1)*clhs48 + clhs50;
            lhs(1,8)=-clhs51;
            lhs(2,0)=clhs11*clhs54;
            lhs(2,1)=clhs38*clhs54;
            lhs(2,2)=clhs52*(clhs0 + clhs36);
            lhs(2,3)=clhs23*clhs53 + clhs55;
            lhs(2,4)=clhs45*clhs53 + clhs56;
            lhs(2,5)=clhs57;
            lhs(2,6)=clhs34*clhs53 + clhs58;
            lhs(2,7)=clhs51*clhs53 + clhs59;
            lhs(2,8)=clhs60;
            lhs(3,0)=DN(1,0)*clhs3 + DN(1,1)*clhs5 + clhs18;
            lhs(3,1)=DN(1,0)*clhs8 + DN(1,1)*clhs10 + clhs40;
            lhs(3,2)=-clhs55;
            lhs(3,3)=DN(1,0)*clhs12 + DN(1,1)*clhs14 + clhs61*mu + clhs62;
            lhs(3,4)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs64;
            lhs(3,5)=-clhs65;
            lhs(3,6)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs68;
            lhs(3,7)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + clhs69;
            lhs(3,8)=-clhs70;
            lhs(4,0)=DN(1,0)*clhs5 + DN(1,1)*clhs35 + clhs19;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs37 + clhs44;
            lhs(4,2)=-clhs56;
            lhs(4,3)=DN(1,0)*clhs14 + DN(1,1)*clhs41 + clhs64;
            lhs(4,4)=DN(1,0)*clhs22 + DN(1,1)*clhs42 + clhs62 + clhs71*mu;
            lhs(4,5)=-clhs72;
            lhs(4,6)=DN(1,0)*clhs26 + DN(1,1)*clhs47 + clhs74;
            lhs(4,7)=DN(1,0)*clhs33 + DN(1,1)*clhs48 + clhs76;
            lhs(4,8)=-clhs77;
            lhs(5,0)=clhs23 + clhs53*clhs55;
            lhs(5,1)=clhs45 + clhs53*clhs56;
            lhs(5,2)=clhs57;
            lhs(5,3)=clhs54*clhs65;
            lhs(5,4)=clhs54*clhs72;
            lhs(5,5)=clhs52*(clhs61 + clhs71);
            lhs(5,6)=clhs53*clhs70 + clhs78;
            lhs(5,7)=clhs53*clhs77 + clhs79;
            lhs(5,8)=clhs80;
            lhs(6,0)=DN(2,0)*clhs3 + DN(2,1)*clhs5 + clhs29;
            lhs(6,1)=DN(2,0)*clhs8 + DN(2,1)*clhs10 + clhs46;
            lhs(6,2)=-clhs58;
            lhs(6,3)=DN(2,0)*clhs12 + DN(2,1)*clhs14 + clhs68;
            lhs(6,4)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs74;
            lhs(6,5)=-clhs78;
            lhs(6,6)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs81*mu + clhs82;
            lhs(6,7)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + clhs83;
            lhs(6,8)=-clhs84;
            lhs(7,0)=DN(2,0)*clhs5 + DN(2,1)*clhs35 + clhs30;
            lhs(7,1)=DN(2,0)*clhs10 + DN(2,1)*clhs37 + clhs50;
            lhs(7,2)=-clhs59;
            lhs(7,3)=DN(2,0)*clhs14 + DN(2,1)*clhs41 + clhs69;
            lhs(7,4)=DN(2,0)*clhs22 + DN(2,1)*clhs42 + clhs76;
            lhs(7,5)=-clhs79;
            lhs(7,6)=DN(2,0)*clhs26 + DN(2,1)*clhs47 + clhs83;
            lhs(7,7)=DN(2,0)*clhs33 + DN(2,1)*clhs48 + clhs82 + clhs85*mu;
            lhs(7,8)=-clhs86;
            lhs(8,0)=clhs34 + clhs53*clhs58;
            lhs(8,1)=clhs51 + clhs53*clhs59;
            lhs(8,2)=clhs60;
            lhs(8,3)=clhs53*clhs78 + clhs70;
            lhs(8,4)=clhs53*clhs79 + clhs77;
            lhs(8,5)=clhs80;
            lhs(8,6)=clhs54*clhs84;
            lhs(8,7)=clhs54*clhs86;
            lhs(8,8)=clhs52*(clhs81 + clhs85);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void Stokes<SymbolicNavierStokesData<2,4>>::ComputeGaussPointLHSContribution(
    SymbolicNavierStokesData<2,4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             bdf0*rho;
const double clhs2 =             pow(N[0], 2)*clhs1;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs4 =             C(0,2)*DN(0,0);
const double clhs5 =             C(2,2)*DN(0,1) + clhs4;
const double clhs6 =             DN(0,0)*mu;
const double clhs7 =             DN(0,1)*clhs6;
const double clhs8 =             C(0,1)*DN(0,1) + clhs4;
const double clhs9 =             C(1,2)*DN(0,1);
const double clhs10 =             C(2,2)*DN(0,0) + clhs9;
const double clhs11 =             DN(0,0)*N[0];
const double clhs12 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs13 =             C(0,2)*DN(1,0);
const double clhs14 =             C(2,2)*DN(1,1) + clhs13;
const double clhs15 =             DN(0,0)*DN(1,0);
const double clhs16 =             N[0]*clhs1;
const double clhs17 =             N[1]*clhs16;
const double clhs18 =             clhs15*mu + clhs17;
const double clhs19 =             DN(1,1)*clhs6;
const double clhs20 =             C(0,1)*DN(1,1) + clhs13;
const double clhs21 =             C(1,2)*DN(1,1);
const double clhs22 =             C(2,2)*DN(1,0) + clhs21;
const double clhs23 =             DN(0,0)*N[1];
const double clhs24 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs25 =             C(0,2)*DN(2,0);
const double clhs26 =             C(2,2)*DN(2,1) + clhs25;
const double clhs27 =             DN(0,0)*DN(2,0);
const double clhs28 =             N[2]*clhs16;
const double clhs29 =             clhs27*mu + clhs28;
const double clhs30 =             DN(2,1)*clhs6;
const double clhs31 =             C(0,1)*DN(2,1) + clhs25;
const double clhs32 =             C(1,2)*DN(2,1);
const double clhs33 =             C(2,2)*DN(2,0) + clhs32;
const double clhs34 =             DN(0,0)*N[2];
const double clhs35 =             C(0,1)*DN(0,0) + clhs9;
const double clhs36 =             pow(DN(0,1), 2);
const double clhs37 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs38 =             DN(0,1)*N[0];
const double clhs39 =             DN(0,1)*mu;
const double clhs40 =             DN(1,0)*clhs39;
const double clhs41 =             C(0,1)*DN(1,0) + clhs21;
const double clhs42 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs43 =             DN(0,1)*DN(1,1);
const double clhs44 =             clhs17 + clhs43*mu;
const double clhs45 =             DN(0,1)*N[1];
const double clhs46 =             DN(2,0)*clhs39;
const double clhs47 =             C(0,1)*DN(2,0) + clhs32;
const double clhs48 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs49 =             DN(0,1)*DN(2,1);
const double clhs50 =             clhs28 + clhs49*mu;
const double clhs51 =             DN(0,1)*N[2];
const double clhs52 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs53 =             clhs1*clhs52;
const double clhs54 =             clhs53 + 1;
const double clhs55 =             DN(1,0)*N[0];
const double clhs56 =             DN(1,1)*N[0];
const double clhs57 =             clhs52*(clhs15 + clhs43);
const double clhs58 =             DN(2,0)*N[0];
const double clhs59 =             DN(2,1)*N[0];
const double clhs60 =             clhs52*(clhs27 + clhs49);
const double clhs61 =             pow(DN(1,0), 2);
const double clhs62 =             pow(N[1], 2)*clhs1;
const double clhs63 =             DN(1,0)*mu;
const double clhs64 =             DN(1,1)*clhs63;
const double clhs65 =             DN(1,0)*N[1];
const double clhs66 =             DN(1,0)*DN(2,0);
const double clhs67 =             N[1]*N[2]*clhs1;
const double clhs68 =             clhs66*mu + clhs67;
const double clhs69 =             DN(2,1)*clhs63;
const double clhs70 =             DN(1,0)*N[2];
const double clhs71 =             pow(DN(1,1), 2);
const double clhs72 =             DN(1,1)*N[1];
const double clhs73 =             DN(2,0)*mu;
const double clhs74 =             DN(1,1)*clhs73;
const double clhs75 =             DN(1,1)*DN(2,1);
const double clhs76 =             clhs67 + clhs75*mu;
const double clhs77 =             DN(1,1)*N[2];
const double clhs78 =             DN(2,0)*N[1];
const double clhs79 =             DN(2,1)*N[1];
const double clhs80 =             clhs52*(clhs66 + clhs75);
const double clhs81 =             pow(DN(2,0), 2);
const double clhs82 =             pow(N[2], 2)*clhs1;
const double clhs83 =             DN(2,1)*clhs73;
const double clhs84 =             DN(2,0)*N[2];
const double clhs85 =             pow(DN(2,1), 2);
const double clhs86 =             DN(2,1)*N[2];
            lhs(0,0)=DN(0,0)*clhs3 + DN(0,1)*clhs5 + clhs0*mu + clhs2;
            lhs(0,1)=DN(0,0)*clhs8 + DN(0,1)*clhs10 + clhs7;
            lhs(0,2)=-clhs11;
            lhs(0,3)=DN(0,0)*clhs12 + DN(0,1)*clhs14 + clhs18;
            lhs(0,4)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs19;
            lhs(0,5)=-clhs23;
            lhs(0,6)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs29;
            lhs(0,7)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + clhs30;
            lhs(0,8)=-clhs34;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs35 + clhs7;
            lhs(1,1)=DN(0,0)*clhs10 + DN(0,1)*clhs37 + clhs2 + clhs36*mu;
            lhs(1,2)=-clhs38;
            lhs(1,3)=DN(0,0)*clhs14 + DN(0,1)*clhs41 + clhs40;
            lhs(1,4)=DN(0,0)*clhs22 + DN(0,1)*clhs42 + clhs44;
            lhs(1,5)=-clhs45;
            lhs(1,6)=DN(0,0)*clhs26 + DN(0,1)*clhs47 + clhs46;
            lhs(1,7)=DN(0,0)*clhs33 + DN(0,1)*clhs48 + clhs50;
            lhs(1,8)=-clhs51;
            lhs(2,0)=clhs11*clhs54;
            lhs(2,1)=clhs38*clhs54;
            lhs(2,2)=clhs52*(clhs0 + clhs36);
            lhs(2,3)=clhs23*clhs53 + clhs55;
            lhs(2,4)=clhs45*clhs53 + clhs56;
            lhs(2,5)=clhs57;
            lhs(2,6)=clhs34*clhs53 + clhs58;
            lhs(2,7)=clhs51*clhs53 + clhs59;
            lhs(2,8)=clhs60;
            lhs(3,0)=DN(1,0)*clhs3 + DN(1,1)*clhs5 + clhs18;
            lhs(3,1)=DN(1,0)*clhs8 + DN(1,1)*clhs10 + clhs40;
            lhs(3,2)=-clhs55;
            lhs(3,3)=DN(1,0)*clhs12 + DN(1,1)*clhs14 + clhs61*mu + clhs62;
            lhs(3,4)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs64;
            lhs(3,5)=-clhs65;
            lhs(3,6)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs68;
            lhs(3,7)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + clhs69;
            lhs(3,8)=-clhs70;
            lhs(4,0)=DN(1,0)*clhs5 + DN(1,1)*clhs35 + clhs19;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs37 + clhs44;
            lhs(4,2)=-clhs56;
            lhs(4,3)=DN(1,0)*clhs14 + DN(1,1)*clhs41 + clhs64;
            lhs(4,4)=DN(1,0)*clhs22 + DN(1,1)*clhs42 + clhs62 + clhs71*mu;
            lhs(4,5)=-clhs72;
            lhs(4,6)=DN(1,0)*clhs26 + DN(1,1)*clhs47 + clhs74;
            lhs(4,7)=DN(1,0)*clhs33 + DN(1,1)*clhs48 + clhs76;
            lhs(4,8)=-clhs77;
            lhs(5,0)=clhs23 + clhs53*clhs55;
            lhs(5,1)=clhs45 + clhs53*clhs56;
            lhs(5,2)=clhs57;
            lhs(5,3)=clhs54*clhs65;
            lhs(5,4)=clhs54*clhs72;
            lhs(5,5)=clhs52*(clhs61 + clhs71);
            lhs(5,6)=clhs53*clhs70 + clhs78;
            lhs(5,7)=clhs53*clhs77 + clhs79;
            lhs(5,8)=clhs80;
            lhs(6,0)=DN(2,0)*clhs3 + DN(2,1)*clhs5 + clhs29;
            lhs(6,1)=DN(2,0)*clhs8 + DN(2,1)*clhs10 + clhs46;
            lhs(6,2)=-clhs58;
            lhs(6,3)=DN(2,0)*clhs12 + DN(2,1)*clhs14 + clhs68;
            lhs(6,4)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs74;
            lhs(6,5)=-clhs78;
            lhs(6,6)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs81*mu + clhs82;
            lhs(6,7)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + clhs83;
            lhs(6,8)=-clhs84;
            lhs(7,0)=DN(2,0)*clhs5 + DN(2,1)*clhs35 + clhs30;
            lhs(7,1)=DN(2,0)*clhs10 + DN(2,1)*clhs37 + clhs50;
            lhs(7,2)=-clhs59;
            lhs(7,3)=DN(2,0)*clhs14 + DN(2,1)*clhs41 + clhs69;
            lhs(7,4)=DN(2,0)*clhs22 + DN(2,1)*clhs42 + clhs76;
            lhs(7,5)=-clhs79;
            lhs(7,6)=DN(2,0)*clhs26 + DN(2,1)*clhs47 + clhs83;
            lhs(7,7)=DN(2,0)*clhs33 + DN(2,1)*clhs48 + clhs82 + clhs85*mu;
            lhs(7,8)=-clhs86;
            lhs(8,0)=clhs34 + clhs53*clhs58;
            lhs(8,1)=clhs51 + clhs53*clhs59;
            lhs(8,2)=clhs60;
            lhs(8,3)=clhs53*clhs78 + clhs70;
            lhs(8,4)=clhs53*clhs79 + clhs77;
            lhs(8,5)=clhs80;
            lhs(8,6)=clhs54*clhs84;
            lhs(8,7)=clhs54*clhs86;
            lhs(8,8)=clhs52*(clhs81 + clhs85);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void Stokes<SymbolicNavierStokesData<3,4>>::ComputeGaussPointLHSContribution(
    SymbolicNavierStokesData<3,4> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             bdf0*rho;
const double clhs2 =             pow(N[0], 2)*clhs1;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs4 =             C(0,3)*DN(0,0);
const double clhs5 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs4;
const double clhs6 =             C(0,5)*DN(0,0);
const double clhs7 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs6;
const double clhs8 =             DN(0,0)*mu;
const double clhs9 =             DN(0,1)*clhs8;
const double clhs10 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs4;
const double clhs11 =             C(1,3)*DN(0,1);
const double clhs12 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs11;
const double clhs13 =             C(3,5)*DN(0,0);
const double clhs14 =             C(4,5)*DN(0,2);
const double clhs15 =             C(1,5)*DN(0,1) + clhs13 + clhs14;
const double clhs16 =             DN(0,2)*clhs8;
const double clhs17 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs6;
const double clhs18 =             C(3,4)*DN(0,1);
const double clhs19 =             C(2,3)*DN(0,2) + clhs13 + clhs18;
const double clhs20 =             C(2,5)*DN(0,2);
const double clhs21 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs20;
const double clhs22 =             DN(0,0)*N[0];
const double clhs23 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs24 =             C(0,3)*DN(1,0);
const double clhs25 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs24;
const double clhs26 =             C(0,5)*DN(1,0);
const double clhs27 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs26;
const double clhs28 =             DN(0,0)*DN(1,0);
const double clhs29 =             N[0]*clhs1;
const double clhs30 =             N[1]*clhs29;
const double clhs31 =             clhs28*mu + clhs30;
const double clhs32 =             DN(1,1)*clhs8;
const double clhs33 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs24;
const double clhs34 =             C(1,3)*DN(1,1);
const double clhs35 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs34;
const double clhs36 =             C(3,5)*DN(1,0);
const double clhs37 =             C(4,5)*DN(1,2);
const double clhs38 =             C(1,5)*DN(1,1) + clhs36 + clhs37;
const double clhs39 =             DN(1,2)*clhs8;
const double clhs40 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs26;
const double clhs41 =             C(3,4)*DN(1,1);
const double clhs42 =             C(2,3)*DN(1,2) + clhs36 + clhs41;
const double clhs43 =             C(2,5)*DN(1,2);
const double clhs44 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs43;
const double clhs45 =             DN(0,0)*N[1];
const double clhs46 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs47 =             C(0,3)*DN(2,0);
const double clhs48 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs47;
const double clhs49 =             C(0,5)*DN(2,0);
const double clhs50 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs49;
const double clhs51 =             DN(0,0)*DN(2,0);
const double clhs52 =             N[2]*clhs29;
const double clhs53 =             clhs51*mu + clhs52;
const double clhs54 =             DN(2,1)*clhs8;
const double clhs55 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs47;
const double clhs56 =             C(1,3)*DN(2,1);
const double clhs57 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs56;
const double clhs58 =             C(3,5)*DN(2,0);
const double clhs59 =             C(4,5)*DN(2,2);
const double clhs60 =             C(1,5)*DN(2,1) + clhs58 + clhs59;
const double clhs61 =             DN(2,2)*clhs8;
const double clhs62 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs49;
const double clhs63 =             C(3,4)*DN(2,1);
const double clhs64 =             C(2,3)*DN(2,2) + clhs58 + clhs63;
const double clhs65 =             C(2,5)*DN(2,2);
const double clhs66 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs65;
const double clhs67 =             DN(0,0)*N[2];
const double clhs68 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs69 =             C(0,3)*DN(3,0);
const double clhs70 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs69;
const double clhs71 =             C(0,5)*DN(3,0);
const double clhs72 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs71;
const double clhs73 =             DN(0,0)*DN(3,0);
const double clhs74 =             N[3]*clhs29;
const double clhs75 =             clhs73*mu + clhs74;
const double clhs76 =             DN(3,1)*clhs8;
const double clhs77 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs69;
const double clhs78 =             C(1,3)*DN(3,1);
const double clhs79 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs78;
const double clhs80 =             C(3,5)*DN(3,0);
const double clhs81 =             C(4,5)*DN(3,2);
const double clhs82 =             C(1,5)*DN(3,1) + clhs80 + clhs81;
const double clhs83 =             DN(3,2)*clhs8;
const double clhs84 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs71;
const double clhs85 =             C(3,4)*DN(3,1);
const double clhs86 =             C(2,3)*DN(3,2) + clhs80 + clhs85;
const double clhs87 =             C(2,5)*DN(3,2);
const double clhs88 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs87;
const double clhs89 =             DN(0,0)*N[3];
const double clhs90 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs11;
const double clhs91 =             C(0,4)*DN(0,0) + clhs14 + clhs18;
const double clhs92 =             pow(DN(0,1), 2);
const double clhs93 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs94 =             C(1,4)*DN(0,1);
const double clhs95 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs94;
const double clhs96 =             DN(0,1)*mu;
const double clhs97 =             DN(0,2)*clhs96;
const double clhs98 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs94;
const double clhs99 =             C(2,4)*DN(0,2);
const double clhs100 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs99;
const double clhs101 =             DN(0,1)*N[0];
const double clhs102 =             DN(1,0)*clhs96;
const double clhs103 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs34;
const double clhs104 =             C(0,4)*DN(1,0) + clhs37 + clhs41;
const double clhs105 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs106 =             C(1,4)*DN(1,1);
const double clhs107 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs106;
const double clhs108 =             DN(0,1)*DN(1,1);
const double clhs109 =             clhs108*mu + clhs30;
const double clhs110 =             DN(1,2)*clhs96;
const double clhs111 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs106;
const double clhs112 =             C(2,4)*DN(1,2);
const double clhs113 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs112;
const double clhs114 =             DN(0,1)*N[1];
const double clhs115 =             DN(2,0)*clhs96;
const double clhs116 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs56;
const double clhs117 =             C(0,4)*DN(2,0) + clhs59 + clhs63;
const double clhs118 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs119 =             C(1,4)*DN(2,1);
const double clhs120 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs119;
const double clhs121 =             DN(0,1)*DN(2,1);
const double clhs122 =             clhs121*mu + clhs52;
const double clhs123 =             DN(2,2)*clhs96;
const double clhs124 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs119;
const double clhs125 =             C(2,4)*DN(2,2);
const double clhs126 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs125;
const double clhs127 =             DN(0,1)*N[2];
const double clhs128 =             DN(3,0)*clhs96;
const double clhs129 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs78;
const double clhs130 =             C(0,4)*DN(3,0) + clhs81 + clhs85;
const double clhs131 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs132 =             C(1,4)*DN(3,1);
const double clhs133 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs132;
const double clhs134 =             DN(0,1)*DN(3,1);
const double clhs135 =             clhs134*mu + clhs74;
const double clhs136 =             DN(3,2)*clhs96;
const double clhs137 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs132;
const double clhs138 =             C(2,4)*DN(3,2);
const double clhs139 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs138;
const double clhs140 =             DN(0,1)*N[3];
const double clhs141 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs20;
const double clhs142 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs99;
const double clhs143 =             pow(DN(0,2), 2);
const double clhs144 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs145 =             DN(0,2)*N[0];
const double clhs146 =             DN(0,2)*mu;
const double clhs147 =             DN(1,0)*clhs146;
const double clhs148 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs43;
const double clhs149 =             DN(1,1)*clhs146;
const double clhs150 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs112;
const double clhs151 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs152 =             DN(0,2)*DN(1,2);
const double clhs153 =             clhs152*mu + clhs30;
const double clhs154 =             DN(0,2)*N[1];
const double clhs155 =             DN(2,0)*clhs146;
const double clhs156 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs65;
const double clhs157 =             DN(2,1)*clhs146;
const double clhs158 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs125;
const double clhs159 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs160 =             DN(0,2)*DN(2,2);
const double clhs161 =             clhs160*mu + clhs52;
const double clhs162 =             DN(0,2)*N[2];
const double clhs163 =             DN(3,0)*clhs146;
const double clhs164 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs87;
const double clhs165 =             DN(3,1)*clhs146;
const double clhs166 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs138;
const double clhs167 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs168 =             DN(0,2)*DN(3,2);
const double clhs169 =             clhs168*mu + clhs74;
const double clhs170 =             DN(0,2)*N[3];
const double clhs171 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs172 =             clhs1*clhs171;
const double clhs173 =             clhs172 + 1;
const double clhs174 =             DN(1,0)*N[0];
const double clhs175 =             DN(1,1)*N[0];
const double clhs176 =             DN(1,2)*N[0];
const double clhs177 =             clhs171*(clhs108 + clhs152 + clhs28);
const double clhs178 =             DN(2,0)*N[0];
const double clhs179 =             DN(2,1)*N[0];
const double clhs180 =             DN(2,2)*N[0];
const double clhs181 =             clhs171*(clhs121 + clhs160 + clhs51);
const double clhs182 =             DN(3,0)*N[0];
const double clhs183 =             DN(3,1)*N[0];
const double clhs184 =             DN(3,2)*N[0];
const double clhs185 =             clhs171*(clhs134 + clhs168 + clhs73);
const double clhs186 =             pow(DN(1,0), 2);
const double clhs187 =             pow(N[1], 2)*clhs1;
const double clhs188 =             DN(1,0)*mu;
const double clhs189 =             DN(1,1)*clhs188;
const double clhs190 =             DN(1,2)*clhs188;
const double clhs191 =             DN(1,0)*N[1];
const double clhs192 =             DN(1,0)*DN(2,0);
const double clhs193 =             N[1]*clhs1;
const double clhs194 =             N[2]*clhs193;
const double clhs195 =             clhs192*mu + clhs194;
const double clhs196 =             DN(2,1)*clhs188;
const double clhs197 =             DN(2,2)*clhs188;
const double clhs198 =             DN(1,0)*N[2];
const double clhs199 =             DN(1,0)*DN(3,0);
const double clhs200 =             N[3]*clhs193;
const double clhs201 =             clhs199*mu + clhs200;
const double clhs202 =             DN(3,1)*clhs188;
const double clhs203 =             DN(3,2)*clhs188;
const double clhs204 =             DN(1,0)*N[3];
const double clhs205 =             pow(DN(1,1), 2);
const double clhs206 =             DN(1,1)*mu;
const double clhs207 =             DN(1,2)*clhs206;
const double clhs208 =             DN(1,1)*N[1];
const double clhs209 =             DN(2,0)*clhs206;
const double clhs210 =             DN(1,1)*DN(2,1);
const double clhs211 =             clhs194 + clhs210*mu;
const double clhs212 =             DN(2,2)*clhs206;
const double clhs213 =             DN(1,1)*N[2];
const double clhs214 =             DN(3,0)*clhs206;
const double clhs215 =             DN(1,1)*DN(3,1);
const double clhs216 =             clhs200 + clhs215*mu;
const double clhs217 =             DN(3,2)*clhs206;
const double clhs218 =             DN(1,1)*N[3];
const double clhs219 =             pow(DN(1,2), 2);
const double clhs220 =             DN(1,2)*N[1];
const double clhs221 =             DN(1,2)*mu;
const double clhs222 =             DN(2,0)*clhs221;
const double clhs223 =             DN(2,1)*clhs221;
const double clhs224 =             DN(1,2)*DN(2,2);
const double clhs225 =             clhs194 + clhs224*mu;
const double clhs226 =             DN(1,2)*N[2];
const double clhs227 =             DN(3,0)*clhs221;
const double clhs228 =             DN(3,1)*clhs221;
const double clhs229 =             DN(1,2)*DN(3,2);
const double clhs230 =             clhs200 + clhs229*mu;
const double clhs231 =             DN(1,2)*N[3];
const double clhs232 =             DN(2,0)*N[1];
const double clhs233 =             DN(2,1)*N[1];
const double clhs234 =             DN(2,2)*N[1];
const double clhs235 =             clhs171*(clhs192 + clhs210 + clhs224);
const double clhs236 =             DN(3,0)*N[1];
const double clhs237 =             DN(3,1)*N[1];
const double clhs238 =             DN(3,2)*N[1];
const double clhs239 =             clhs171*(clhs199 + clhs215 + clhs229);
const double clhs240 =             pow(DN(2,0), 2);
const double clhs241 =             pow(N[2], 2)*clhs1;
const double clhs242 =             DN(2,0)*mu;
const double clhs243 =             DN(2,1)*clhs242;
const double clhs244 =             DN(2,2)*clhs242;
const double clhs245 =             DN(2,0)*N[2];
const double clhs246 =             DN(2,0)*DN(3,0);
const double clhs247 =             N[2]*N[3]*clhs1;
const double clhs248 =             clhs246*mu + clhs247;
const double clhs249 =             DN(3,1)*clhs242;
const double clhs250 =             DN(3,2)*clhs242;
const double clhs251 =             DN(2,0)*N[3];
const double clhs252 =             pow(DN(2,1), 2);
const double clhs253 =             DN(2,1)*mu;
const double clhs254 =             DN(2,2)*clhs253;
const double clhs255 =             DN(2,1)*N[2];
const double clhs256 =             DN(3,0)*clhs253;
const double clhs257 =             DN(2,1)*DN(3,1);
const double clhs258 =             clhs247 + clhs257*mu;
const double clhs259 =             DN(3,2)*clhs253;
const double clhs260 =             DN(2,1)*N[3];
const double clhs261 =             pow(DN(2,2), 2);
const double clhs262 =             DN(2,2)*N[2];
const double clhs263 =             DN(2,2)*mu;
const double clhs264 =             DN(3,0)*clhs263;
const double clhs265 =             DN(3,1)*clhs263;
const double clhs266 =             DN(2,2)*DN(3,2);
const double clhs267 =             clhs247 + clhs266*mu;
const double clhs268 =             DN(2,2)*N[3];
const double clhs269 =             DN(3,0)*N[2];
const double clhs270 =             DN(3,1)*N[2];
const double clhs271 =             DN(3,2)*N[2];
const double clhs272 =             clhs171*(clhs246 + clhs257 + clhs266);
const double clhs273 =             pow(DN(3,0), 2);
const double clhs274 =             pow(N[3], 2)*clhs1;
const double clhs275 =             DN(3,0)*mu;
const double clhs276 =             DN(3,1)*clhs275;
const double clhs277 =             DN(3,2)*clhs275;
const double clhs278 =             DN(3,0)*N[3];
const double clhs279 =             pow(DN(3,1), 2);
const double clhs280 =             DN(3,1)*DN(3,2)*mu;
const double clhs281 =             DN(3,1)*N[3];
const double clhs282 =             pow(DN(3,2), 2);
const double clhs283 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs3 + DN(0,1)*clhs5 + DN(0,2)*clhs7 + clhs0*mu + clhs2;
            lhs(0,1)=DN(0,0)*clhs10 + DN(0,1)*clhs12 + DN(0,2)*clhs15 + clhs9;
            lhs(0,2)=DN(0,0)*clhs17 + DN(0,1)*clhs19 + DN(0,2)*clhs21 + clhs16;
            lhs(0,3)=-clhs22;
            lhs(0,4)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs27 + clhs31;
            lhs(0,5)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + DN(0,2)*clhs38 + clhs32;
            lhs(0,6)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs44 + clhs39;
            lhs(0,7)=-clhs45;
            lhs(0,8)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs50 + clhs53;
            lhs(0,9)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs54;
            lhs(0,10)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs61;
            lhs(0,11)=-clhs67;
            lhs(0,12)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs75;
            lhs(0,13)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs82 + clhs76;
            lhs(0,14)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs88 + clhs83;
            lhs(0,15)=-clhs89;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs90 + DN(0,2)*clhs91 + clhs9;
            lhs(1,1)=DN(0,0)*clhs12 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs2 + clhs92*mu;
            lhs(1,2)=DN(0,0)*clhs19 + DN(0,1)*clhs98 + DN(0,2)*clhs100 + clhs97;
            lhs(1,3)=-clhs101;
            lhs(1,4)=DN(0,0)*clhs25 + DN(0,1)*clhs103 + DN(0,2)*clhs104 + clhs102;
            lhs(1,5)=DN(0,0)*clhs35 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + clhs109;
            lhs(1,6)=DN(0,0)*clhs42 + DN(0,1)*clhs111 + DN(0,2)*clhs113 + clhs110;
            lhs(1,7)=-clhs114;
            lhs(1,8)=DN(0,0)*clhs48 + DN(0,1)*clhs116 + DN(0,2)*clhs117 + clhs115;
            lhs(1,9)=DN(0,0)*clhs57 + DN(0,1)*clhs118 + DN(0,2)*clhs120 + clhs122;
            lhs(1,10)=DN(0,0)*clhs64 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs123;
            lhs(1,11)=-clhs127;
            lhs(1,12)=DN(0,0)*clhs70 + DN(0,1)*clhs129 + DN(0,2)*clhs130 + clhs128;
            lhs(1,13)=DN(0,0)*clhs79 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs135;
            lhs(1,14)=DN(0,0)*clhs86 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs136;
            lhs(1,15)=-clhs140;
            lhs(2,0)=DN(0,0)*clhs7 + DN(0,1)*clhs91 + DN(0,2)*clhs141 + clhs16;
            lhs(2,1)=DN(0,0)*clhs15 + DN(0,1)*clhs95 + DN(0,2)*clhs142 + clhs97;
            lhs(2,2)=DN(0,0)*clhs21 + DN(0,1)*clhs100 + DN(0,2)*clhs144 + clhs143*mu + clhs2;
            lhs(2,3)=-clhs145;
            lhs(2,4)=DN(0,0)*clhs27 + DN(0,1)*clhs104 + DN(0,2)*clhs148 + clhs147;
            lhs(2,5)=DN(0,0)*clhs38 + DN(0,1)*clhs107 + DN(0,2)*clhs150 + clhs149;
            lhs(2,6)=DN(0,0)*clhs44 + DN(0,1)*clhs113 + DN(0,2)*clhs151 + clhs153;
            lhs(2,7)=-clhs154;
            lhs(2,8)=DN(0,0)*clhs50 + DN(0,1)*clhs117 + DN(0,2)*clhs156 + clhs155;
            lhs(2,9)=DN(0,0)*clhs60 + DN(0,1)*clhs120 + DN(0,2)*clhs158 + clhs157;
            lhs(2,10)=DN(0,0)*clhs66 + DN(0,1)*clhs126 + DN(0,2)*clhs159 + clhs161;
            lhs(2,11)=-clhs162;
            lhs(2,12)=DN(0,0)*clhs72 + DN(0,1)*clhs130 + DN(0,2)*clhs164 + clhs163;
            lhs(2,13)=DN(0,0)*clhs82 + DN(0,1)*clhs133 + DN(0,2)*clhs166 + clhs165;
            lhs(2,14)=DN(0,0)*clhs88 + DN(0,1)*clhs139 + DN(0,2)*clhs167 + clhs169;
            lhs(2,15)=-clhs170;
            lhs(3,0)=clhs173*clhs22;
            lhs(3,1)=clhs101*clhs173;
            lhs(3,2)=clhs145*clhs173;
            lhs(3,3)=clhs171*(clhs0 + clhs143 + clhs92);
            lhs(3,4)=clhs172*clhs45 + clhs174;
            lhs(3,5)=clhs114*clhs172 + clhs175;
            lhs(3,6)=clhs154*clhs172 + clhs176;
            lhs(3,7)=clhs177;
            lhs(3,8)=clhs172*clhs67 + clhs178;
            lhs(3,9)=clhs127*clhs172 + clhs179;
            lhs(3,10)=clhs162*clhs172 + clhs180;
            lhs(3,11)=clhs181;
            lhs(3,12)=clhs172*clhs89 + clhs182;
            lhs(3,13)=clhs140*clhs172 + clhs183;
            lhs(3,14)=clhs170*clhs172 + clhs184;
            lhs(3,15)=clhs185;
            lhs(4,0)=DN(1,0)*clhs3 + DN(1,1)*clhs5 + DN(1,2)*clhs7 + clhs31;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs12 + DN(1,2)*clhs15 + clhs102;
            lhs(4,2)=DN(1,0)*clhs17 + DN(1,1)*clhs19 + DN(1,2)*clhs21 + clhs147;
            lhs(4,3)=-clhs174;
            lhs(4,4)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs27 + clhs186*mu + clhs187;
            lhs(4,5)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + DN(1,2)*clhs38 + clhs189;
            lhs(4,6)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs44 + clhs190;
            lhs(4,7)=-clhs191;
            lhs(4,8)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs50 + clhs195;
            lhs(4,9)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs196;
            lhs(4,10)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs197;
            lhs(4,11)=-clhs198;
            lhs(4,12)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs201;
            lhs(4,13)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs82 + clhs202;
            lhs(4,14)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs88 + clhs203;
            lhs(4,15)=-clhs204;
            lhs(5,0)=DN(1,0)*clhs5 + DN(1,1)*clhs90 + DN(1,2)*clhs91 + clhs32;
            lhs(5,1)=DN(1,0)*clhs12 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs109;
            lhs(5,2)=DN(1,0)*clhs19 + DN(1,1)*clhs98 + DN(1,2)*clhs100 + clhs149;
            lhs(5,3)=-clhs175;
            lhs(5,4)=DN(1,0)*clhs25 + DN(1,1)*clhs103 + DN(1,2)*clhs104 + clhs189;
            lhs(5,5)=DN(1,0)*clhs35 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + clhs187 + clhs205*mu;
            lhs(5,6)=DN(1,0)*clhs42 + DN(1,1)*clhs111 + DN(1,2)*clhs113 + clhs207;
            lhs(5,7)=-clhs208;
            lhs(5,8)=DN(1,0)*clhs48 + DN(1,1)*clhs116 + DN(1,2)*clhs117 + clhs209;
            lhs(5,9)=DN(1,0)*clhs57 + DN(1,1)*clhs118 + DN(1,2)*clhs120 + clhs211;
            lhs(5,10)=DN(1,0)*clhs64 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs212;
            lhs(5,11)=-clhs213;
            lhs(5,12)=DN(1,0)*clhs70 + DN(1,1)*clhs129 + DN(1,2)*clhs130 + clhs214;
            lhs(5,13)=DN(1,0)*clhs79 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs216;
            lhs(5,14)=DN(1,0)*clhs86 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs217;
            lhs(5,15)=-clhs218;
            lhs(6,0)=DN(1,0)*clhs7 + DN(1,1)*clhs91 + DN(1,2)*clhs141 + clhs39;
            lhs(6,1)=DN(1,0)*clhs15 + DN(1,1)*clhs95 + DN(1,2)*clhs142 + clhs110;
            lhs(6,2)=DN(1,0)*clhs21 + DN(1,1)*clhs100 + DN(1,2)*clhs144 + clhs153;
            lhs(6,3)=-clhs176;
            lhs(6,4)=DN(1,0)*clhs27 + DN(1,1)*clhs104 + DN(1,2)*clhs148 + clhs190;
            lhs(6,5)=DN(1,0)*clhs38 + DN(1,1)*clhs107 + DN(1,2)*clhs150 + clhs207;
            lhs(6,6)=DN(1,0)*clhs44 + DN(1,1)*clhs113 + DN(1,2)*clhs151 + clhs187 + clhs219*mu;
            lhs(6,7)=-clhs220;
            lhs(6,8)=DN(1,0)*clhs50 + DN(1,1)*clhs117 + DN(1,2)*clhs156 + clhs222;
            lhs(6,9)=DN(1,0)*clhs60 + DN(1,1)*clhs120 + DN(1,2)*clhs158 + clhs223;
            lhs(6,10)=DN(1,0)*clhs66 + DN(1,1)*clhs126 + DN(1,2)*clhs159 + clhs225;
            lhs(6,11)=-clhs226;
            lhs(6,12)=DN(1,0)*clhs72 + DN(1,1)*clhs130 + DN(1,2)*clhs164 + clhs227;
            lhs(6,13)=DN(1,0)*clhs82 + DN(1,1)*clhs133 + DN(1,2)*clhs166 + clhs228;
            lhs(6,14)=DN(1,0)*clhs88 + DN(1,1)*clhs139 + DN(1,2)*clhs167 + clhs230;
            lhs(6,15)=-clhs231;
            lhs(7,0)=clhs172*clhs174 + clhs45;
            lhs(7,1)=clhs114 + clhs172*clhs175;
            lhs(7,2)=clhs154 + clhs172*clhs176;
            lhs(7,3)=clhs177;
            lhs(7,4)=clhs173*clhs191;
            lhs(7,5)=clhs173*clhs208;
            lhs(7,6)=clhs173*clhs220;
            lhs(7,7)=clhs171*(clhs186 + clhs205 + clhs219);
            lhs(7,8)=clhs172*clhs198 + clhs232;
            lhs(7,9)=clhs172*clhs213 + clhs233;
            lhs(7,10)=clhs172*clhs226 + clhs234;
            lhs(7,11)=clhs235;
            lhs(7,12)=clhs172*clhs204 + clhs236;
            lhs(7,13)=clhs172*clhs218 + clhs237;
            lhs(7,14)=clhs172*clhs231 + clhs238;
            lhs(7,15)=clhs239;
            lhs(8,0)=DN(2,0)*clhs3 + DN(2,1)*clhs5 + DN(2,2)*clhs7 + clhs53;
            lhs(8,1)=DN(2,0)*clhs10 + DN(2,1)*clhs12 + DN(2,2)*clhs15 + clhs115;
            lhs(8,2)=DN(2,0)*clhs17 + DN(2,1)*clhs19 + DN(2,2)*clhs21 + clhs155;
            lhs(8,3)=-clhs178;
            lhs(8,4)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs27 + clhs195;
            lhs(8,5)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + DN(2,2)*clhs38 + clhs209;
            lhs(8,6)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs44 + clhs222;
            lhs(8,7)=-clhs232;
            lhs(8,8)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs50 + clhs240*mu + clhs241;
            lhs(8,9)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs243;
            lhs(8,10)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs244;
            lhs(8,11)=-clhs245;
            lhs(8,12)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs248;
            lhs(8,13)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs82 + clhs249;
            lhs(8,14)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs88 + clhs250;
            lhs(8,15)=-clhs251;
            lhs(9,0)=DN(2,0)*clhs5 + DN(2,1)*clhs90 + DN(2,2)*clhs91 + clhs54;
            lhs(9,1)=DN(2,0)*clhs12 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs122;
            lhs(9,2)=DN(2,0)*clhs19 + DN(2,1)*clhs98 + DN(2,2)*clhs100 + clhs157;
            lhs(9,3)=-clhs179;
            lhs(9,4)=DN(2,0)*clhs25 + DN(2,1)*clhs103 + DN(2,2)*clhs104 + clhs196;
            lhs(9,5)=DN(2,0)*clhs35 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs211;
            lhs(9,6)=DN(2,0)*clhs42 + DN(2,1)*clhs111 + DN(2,2)*clhs113 + clhs223;
            lhs(9,7)=-clhs233;
            lhs(9,8)=DN(2,0)*clhs48 + DN(2,1)*clhs116 + DN(2,2)*clhs117 + clhs243;
            lhs(9,9)=DN(2,0)*clhs57 + DN(2,1)*clhs118 + DN(2,2)*clhs120 + clhs241 + clhs252*mu;
            lhs(9,10)=DN(2,0)*clhs64 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs254;
            lhs(9,11)=-clhs255;
            lhs(9,12)=DN(2,0)*clhs70 + DN(2,1)*clhs129 + DN(2,2)*clhs130 + clhs256;
            lhs(9,13)=DN(2,0)*clhs79 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs258;
            lhs(9,14)=DN(2,0)*clhs86 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs259;
            lhs(9,15)=-clhs260;
            lhs(10,0)=DN(2,0)*clhs7 + DN(2,1)*clhs91 + DN(2,2)*clhs141 + clhs61;
            lhs(10,1)=DN(2,0)*clhs15 + DN(2,1)*clhs95 + DN(2,2)*clhs142 + clhs123;
            lhs(10,2)=DN(2,0)*clhs21 + DN(2,1)*clhs100 + DN(2,2)*clhs144 + clhs161;
            lhs(10,3)=-clhs180;
            lhs(10,4)=DN(2,0)*clhs27 + DN(2,1)*clhs104 + DN(2,2)*clhs148 + clhs197;
            lhs(10,5)=DN(2,0)*clhs38 + DN(2,1)*clhs107 + DN(2,2)*clhs150 + clhs212;
            lhs(10,6)=DN(2,0)*clhs44 + DN(2,1)*clhs113 + DN(2,2)*clhs151 + clhs225;
            lhs(10,7)=-clhs234;
            lhs(10,8)=DN(2,0)*clhs50 + DN(2,1)*clhs117 + DN(2,2)*clhs156 + clhs244;
            lhs(10,9)=DN(2,0)*clhs60 + DN(2,1)*clhs120 + DN(2,2)*clhs158 + clhs254;
            lhs(10,10)=DN(2,0)*clhs66 + DN(2,1)*clhs126 + DN(2,2)*clhs159 + clhs241 + clhs261*mu;
            lhs(10,11)=-clhs262;
            lhs(10,12)=DN(2,0)*clhs72 + DN(2,1)*clhs130 + DN(2,2)*clhs164 + clhs264;
            lhs(10,13)=DN(2,0)*clhs82 + DN(2,1)*clhs133 + DN(2,2)*clhs166 + clhs265;
            lhs(10,14)=DN(2,0)*clhs88 + DN(2,1)*clhs139 + DN(2,2)*clhs167 + clhs267;
            lhs(10,15)=-clhs268;
            lhs(11,0)=clhs172*clhs178 + clhs67;
            lhs(11,1)=clhs127 + clhs172*clhs179;
            lhs(11,2)=clhs162 + clhs172*clhs180;
            lhs(11,3)=clhs181;
            lhs(11,4)=clhs172*clhs232 + clhs198;
            lhs(11,5)=clhs172*clhs233 + clhs213;
            lhs(11,6)=clhs172*clhs234 + clhs226;
            lhs(11,7)=clhs235;
            lhs(11,8)=clhs173*clhs245;
            lhs(11,9)=clhs173*clhs255;
            lhs(11,10)=clhs173*clhs262;
            lhs(11,11)=clhs171*(clhs240 + clhs252 + clhs261);
            lhs(11,12)=clhs172*clhs251 + clhs269;
            lhs(11,13)=clhs172*clhs260 + clhs270;
            lhs(11,14)=clhs172*clhs268 + clhs271;
            lhs(11,15)=clhs272;
            lhs(12,0)=DN(3,0)*clhs3 + DN(3,1)*clhs5 + DN(3,2)*clhs7 + clhs75;
            lhs(12,1)=DN(3,0)*clhs10 + DN(3,1)*clhs12 + DN(3,2)*clhs15 + clhs128;
            lhs(12,2)=DN(3,0)*clhs17 + DN(3,1)*clhs19 + DN(3,2)*clhs21 + clhs163;
            lhs(12,3)=-clhs182;
            lhs(12,4)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs27 + clhs201;
            lhs(12,5)=DN(3,0)*clhs33 + DN(3,1)*clhs35 + DN(3,2)*clhs38 + clhs214;
            lhs(12,6)=DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs44 + clhs227;
            lhs(12,7)=-clhs236;
            lhs(12,8)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs50 + clhs248;
            lhs(12,9)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs256;
            lhs(12,10)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs264;
            lhs(12,11)=-clhs269;
            lhs(12,12)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs273*mu + clhs274;
            lhs(12,13)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs82 + clhs276;
            lhs(12,14)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs88 + clhs277;
            lhs(12,15)=-clhs278;
            lhs(13,0)=DN(3,0)*clhs5 + DN(3,1)*clhs90 + DN(3,2)*clhs91 + clhs76;
            lhs(13,1)=DN(3,0)*clhs12 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs135;
            lhs(13,2)=DN(3,0)*clhs19 + DN(3,1)*clhs98 + DN(3,2)*clhs100 + clhs165;
            lhs(13,3)=-clhs183;
            lhs(13,4)=DN(3,0)*clhs25 + DN(3,1)*clhs103 + DN(3,2)*clhs104 + clhs202;
            lhs(13,5)=DN(3,0)*clhs35 + DN(3,1)*clhs105 + DN(3,2)*clhs107 + clhs216;
            lhs(13,6)=DN(3,0)*clhs42 + DN(3,1)*clhs111 + DN(3,2)*clhs113 + clhs228;
            lhs(13,7)=-clhs237;
            lhs(13,8)=DN(3,0)*clhs48 + DN(3,1)*clhs116 + DN(3,2)*clhs117 + clhs249;
            lhs(13,9)=DN(3,0)*clhs57 + DN(3,1)*clhs118 + DN(3,2)*clhs120 + clhs258;
            lhs(13,10)=DN(3,0)*clhs64 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs265;
            lhs(13,11)=-clhs270;
            lhs(13,12)=DN(3,0)*clhs70 + DN(3,1)*clhs129 + DN(3,2)*clhs130 + clhs276;
            lhs(13,13)=DN(3,0)*clhs79 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs274 + clhs279*mu;
            lhs(13,14)=DN(3,0)*clhs86 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs280;
            lhs(13,15)=-clhs281;
            lhs(14,0)=DN(3,0)*clhs7 + DN(3,1)*clhs91 + DN(3,2)*clhs141 + clhs83;
            lhs(14,1)=DN(3,0)*clhs15 + DN(3,1)*clhs95 + DN(3,2)*clhs142 + clhs136;
            lhs(14,2)=DN(3,0)*clhs21 + DN(3,1)*clhs100 + DN(3,2)*clhs144 + clhs169;
            lhs(14,3)=-clhs184;
            lhs(14,4)=DN(3,0)*clhs27 + DN(3,1)*clhs104 + DN(3,2)*clhs148 + clhs203;
            lhs(14,5)=DN(3,0)*clhs38 + DN(3,1)*clhs107 + DN(3,2)*clhs150 + clhs217;
            lhs(14,6)=DN(3,0)*clhs44 + DN(3,1)*clhs113 + DN(3,2)*clhs151 + clhs230;
            lhs(14,7)=-clhs238;
            lhs(14,8)=DN(3,0)*clhs50 + DN(3,1)*clhs117 + DN(3,2)*clhs156 + clhs250;
            lhs(14,9)=DN(3,0)*clhs60 + DN(3,1)*clhs120 + DN(3,2)*clhs158 + clhs259;
            lhs(14,10)=DN(3,0)*clhs66 + DN(3,1)*clhs126 + DN(3,2)*clhs159 + clhs267;
            lhs(14,11)=-clhs271;
            lhs(14,12)=DN(3,0)*clhs72 + DN(3,1)*clhs130 + DN(3,2)*clhs164 + clhs277;
            lhs(14,13)=DN(3,0)*clhs82 + DN(3,1)*clhs133 + DN(3,2)*clhs166 + clhs280;
            lhs(14,14)=DN(3,0)*clhs88 + DN(3,1)*clhs139 + DN(3,2)*clhs167 + clhs274 + clhs282*mu;
            lhs(14,15)=-clhs283;
            lhs(15,0)=clhs172*clhs182 + clhs89;
            lhs(15,1)=clhs140 + clhs172*clhs183;
            lhs(15,2)=clhs170 + clhs172*clhs184;
            lhs(15,3)=clhs185;
            lhs(15,4)=clhs172*clhs236 + clhs204;
            lhs(15,5)=clhs172*clhs237 + clhs218;
            lhs(15,6)=clhs172*clhs238 + clhs231;
            lhs(15,7)=clhs239;
            lhs(15,8)=clhs172*clhs269 + clhs251;
            lhs(15,9)=clhs172*clhs270 + clhs260;
            lhs(15,10)=clhs172*clhs271 + clhs268;
            lhs(15,11)=clhs272;
            lhs(15,12)=clhs173*clhs278;
            lhs(15,13)=clhs173*clhs281;
            lhs(15,14)=clhs173*clhs283;
            lhs(15,15)=clhs171*(clhs273 + clhs279 + clhs282);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void Stokes<SymbolicNavierStokesData<3,6>>::ComputeGaussPointLHSContribution(
    SymbolicNavierStokesData<3,6> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             bdf0*rho;
const double clhs2 =             pow(N[0], 2)*clhs1;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs4 =             C(0,3)*DN(0,0);
const double clhs5 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs4;
const double clhs6 =             C(0,5)*DN(0,0);
const double clhs7 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs6;
const double clhs8 =             DN(0,0)*mu;
const double clhs9 =             DN(0,1)*clhs8;
const double clhs10 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs4;
const double clhs11 =             C(1,3)*DN(0,1);
const double clhs12 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs11;
const double clhs13 =             C(3,5)*DN(0,0);
const double clhs14 =             C(4,5)*DN(0,2);
const double clhs15 =             C(1,5)*DN(0,1) + clhs13 + clhs14;
const double clhs16 =             DN(0,2)*clhs8;
const double clhs17 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs6;
const double clhs18 =             C(3,4)*DN(0,1);
const double clhs19 =             C(2,3)*DN(0,2) + clhs13 + clhs18;
const double clhs20 =             C(2,5)*DN(0,2);
const double clhs21 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs20;
const double clhs22 =             DN(0,0)*N[0];
const double clhs23 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs24 =             C(0,3)*DN(1,0);
const double clhs25 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs24;
const double clhs26 =             C(0,5)*DN(1,0);
const double clhs27 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs26;
const double clhs28 =             DN(0,0)*DN(1,0);
const double clhs29 =             N[0]*clhs1;
const double clhs30 =             N[1]*clhs29;
const double clhs31 =             clhs28*mu + clhs30;
const double clhs32 =             DN(1,1)*clhs8;
const double clhs33 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs24;
const double clhs34 =             C(1,3)*DN(1,1);
const double clhs35 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs34;
const double clhs36 =             C(3,5)*DN(1,0);
const double clhs37 =             C(4,5)*DN(1,2);
const double clhs38 =             C(1,5)*DN(1,1) + clhs36 + clhs37;
const double clhs39 =             DN(1,2)*clhs8;
const double clhs40 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs26;
const double clhs41 =             C(3,4)*DN(1,1);
const double clhs42 =             C(2,3)*DN(1,2) + clhs36 + clhs41;
const double clhs43 =             C(2,5)*DN(1,2);
const double clhs44 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs43;
const double clhs45 =             DN(0,0)*N[1];
const double clhs46 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs47 =             C(0,3)*DN(2,0);
const double clhs48 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs47;
const double clhs49 =             C(0,5)*DN(2,0);
const double clhs50 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs49;
const double clhs51 =             DN(0,0)*DN(2,0);
const double clhs52 =             N[2]*clhs29;
const double clhs53 =             clhs51*mu + clhs52;
const double clhs54 =             DN(2,1)*clhs8;
const double clhs55 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs47;
const double clhs56 =             C(1,3)*DN(2,1);
const double clhs57 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs56;
const double clhs58 =             C(3,5)*DN(2,0);
const double clhs59 =             C(4,5)*DN(2,2);
const double clhs60 =             C(1,5)*DN(2,1) + clhs58 + clhs59;
const double clhs61 =             DN(2,2)*clhs8;
const double clhs62 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs49;
const double clhs63 =             C(3,4)*DN(2,1);
const double clhs64 =             C(2,3)*DN(2,2) + clhs58 + clhs63;
const double clhs65 =             C(2,5)*DN(2,2);
const double clhs66 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs65;
const double clhs67 =             DN(0,0)*N[2];
const double clhs68 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs69 =             C(0,3)*DN(3,0);
const double clhs70 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs69;
const double clhs71 =             C(0,5)*DN(3,0);
const double clhs72 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs71;
const double clhs73 =             DN(0,0)*DN(3,0);
const double clhs74 =             N[3]*clhs29;
const double clhs75 =             clhs73*mu + clhs74;
const double clhs76 =             DN(3,1)*clhs8;
const double clhs77 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs69;
const double clhs78 =             C(1,3)*DN(3,1);
const double clhs79 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs78;
const double clhs80 =             C(3,5)*DN(3,0);
const double clhs81 =             C(4,5)*DN(3,2);
const double clhs82 =             C(1,5)*DN(3,1) + clhs80 + clhs81;
const double clhs83 =             DN(3,2)*clhs8;
const double clhs84 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs71;
const double clhs85 =             C(3,4)*DN(3,1);
const double clhs86 =             C(2,3)*DN(3,2) + clhs80 + clhs85;
const double clhs87 =             C(2,5)*DN(3,2);
const double clhs88 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs87;
const double clhs89 =             DN(0,0)*N[3];
const double clhs90 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs11;
const double clhs91 =             C(0,4)*DN(0,0) + clhs14 + clhs18;
const double clhs92 =             pow(DN(0,1), 2);
const double clhs93 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs94 =             C(1,4)*DN(0,1);
const double clhs95 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs94;
const double clhs96 =             DN(0,1)*mu;
const double clhs97 =             DN(0,2)*clhs96;
const double clhs98 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs94;
const double clhs99 =             C(2,4)*DN(0,2);
const double clhs100 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs99;
const double clhs101 =             DN(0,1)*N[0];
const double clhs102 =             DN(1,0)*clhs96;
const double clhs103 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs34;
const double clhs104 =             C(0,4)*DN(1,0) + clhs37 + clhs41;
const double clhs105 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs106 =             C(1,4)*DN(1,1);
const double clhs107 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs106;
const double clhs108 =             DN(0,1)*DN(1,1);
const double clhs109 =             clhs108*mu + clhs30;
const double clhs110 =             DN(1,2)*clhs96;
const double clhs111 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs106;
const double clhs112 =             C(2,4)*DN(1,2);
const double clhs113 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs112;
const double clhs114 =             DN(0,1)*N[1];
const double clhs115 =             DN(2,0)*clhs96;
const double clhs116 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs56;
const double clhs117 =             C(0,4)*DN(2,0) + clhs59 + clhs63;
const double clhs118 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs119 =             C(1,4)*DN(2,1);
const double clhs120 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs119;
const double clhs121 =             DN(0,1)*DN(2,1);
const double clhs122 =             clhs121*mu + clhs52;
const double clhs123 =             DN(2,2)*clhs96;
const double clhs124 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs119;
const double clhs125 =             C(2,4)*DN(2,2);
const double clhs126 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs125;
const double clhs127 =             DN(0,1)*N[2];
const double clhs128 =             DN(3,0)*clhs96;
const double clhs129 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs78;
const double clhs130 =             C(0,4)*DN(3,0) + clhs81 + clhs85;
const double clhs131 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs132 =             C(1,4)*DN(3,1);
const double clhs133 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs132;
const double clhs134 =             DN(0,1)*DN(3,1);
const double clhs135 =             clhs134*mu + clhs74;
const double clhs136 =             DN(3,2)*clhs96;
const double clhs137 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs132;
const double clhs138 =             C(2,4)*DN(3,2);
const double clhs139 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs138;
const double clhs140 =             DN(0,1)*N[3];
const double clhs141 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs20;
const double clhs142 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs99;
const double clhs143 =             pow(DN(0,2), 2);
const double clhs144 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs145 =             DN(0,2)*N[0];
const double clhs146 =             DN(0,2)*mu;
const double clhs147 =             DN(1,0)*clhs146;
const double clhs148 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs43;
const double clhs149 =             DN(1,1)*clhs146;
const double clhs150 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs112;
const double clhs151 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs152 =             DN(0,2)*DN(1,2);
const double clhs153 =             clhs152*mu + clhs30;
const double clhs154 =             DN(0,2)*N[1];
const double clhs155 =             DN(2,0)*clhs146;
const double clhs156 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs65;
const double clhs157 =             DN(2,1)*clhs146;
const double clhs158 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs125;
const double clhs159 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs160 =             DN(0,2)*DN(2,2);
const double clhs161 =             clhs160*mu + clhs52;
const double clhs162 =             DN(0,2)*N[2];
const double clhs163 =             DN(3,0)*clhs146;
const double clhs164 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs87;
const double clhs165 =             DN(3,1)*clhs146;
const double clhs166 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs138;
const double clhs167 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs168 =             DN(0,2)*DN(3,2);
const double clhs169 =             clhs168*mu + clhs74;
const double clhs170 =             DN(0,2)*N[3];
const double clhs171 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs172 =             clhs1*clhs171;
const double clhs173 =             clhs172 + 1;
const double clhs174 =             DN(1,0)*N[0];
const double clhs175 =             DN(1,1)*N[0];
const double clhs176 =             DN(1,2)*N[0];
const double clhs177 =             clhs171*(clhs108 + clhs152 + clhs28);
const double clhs178 =             DN(2,0)*N[0];
const double clhs179 =             DN(2,1)*N[0];
const double clhs180 =             DN(2,2)*N[0];
const double clhs181 =             clhs171*(clhs121 + clhs160 + clhs51);
const double clhs182 =             DN(3,0)*N[0];
const double clhs183 =             DN(3,1)*N[0];
const double clhs184 =             DN(3,2)*N[0];
const double clhs185 =             clhs171*(clhs134 + clhs168 + clhs73);
const double clhs186 =             pow(DN(1,0), 2);
const double clhs187 =             pow(N[1], 2)*clhs1;
const double clhs188 =             DN(1,0)*mu;
const double clhs189 =             DN(1,1)*clhs188;
const double clhs190 =             DN(1,2)*clhs188;
const double clhs191 =             DN(1,0)*N[1];
const double clhs192 =             DN(1,0)*DN(2,0);
const double clhs193 =             N[1]*clhs1;
const double clhs194 =             N[2]*clhs193;
const double clhs195 =             clhs192*mu + clhs194;
const double clhs196 =             DN(2,1)*clhs188;
const double clhs197 =             DN(2,2)*clhs188;
const double clhs198 =             DN(1,0)*N[2];
const double clhs199 =             DN(1,0)*DN(3,0);
const double clhs200 =             N[3]*clhs193;
const double clhs201 =             clhs199*mu + clhs200;
const double clhs202 =             DN(3,1)*clhs188;
const double clhs203 =             DN(3,2)*clhs188;
const double clhs204 =             DN(1,0)*N[3];
const double clhs205 =             pow(DN(1,1), 2);
const double clhs206 =             DN(1,1)*mu;
const double clhs207 =             DN(1,2)*clhs206;
const double clhs208 =             DN(1,1)*N[1];
const double clhs209 =             DN(2,0)*clhs206;
const double clhs210 =             DN(1,1)*DN(2,1);
const double clhs211 =             clhs194 + clhs210*mu;
const double clhs212 =             DN(2,2)*clhs206;
const double clhs213 =             DN(1,1)*N[2];
const double clhs214 =             DN(3,0)*clhs206;
const double clhs215 =             DN(1,1)*DN(3,1);
const double clhs216 =             clhs200 + clhs215*mu;
const double clhs217 =             DN(3,2)*clhs206;
const double clhs218 =             DN(1,1)*N[3];
const double clhs219 =             pow(DN(1,2), 2);
const double clhs220 =             DN(1,2)*N[1];
const double clhs221 =             DN(1,2)*mu;
const double clhs222 =             DN(2,0)*clhs221;
const double clhs223 =             DN(2,1)*clhs221;
const double clhs224 =             DN(1,2)*DN(2,2);
const double clhs225 =             clhs194 + clhs224*mu;
const double clhs226 =             DN(1,2)*N[2];
const double clhs227 =             DN(3,0)*clhs221;
const double clhs228 =             DN(3,1)*clhs221;
const double clhs229 =             DN(1,2)*DN(3,2);
const double clhs230 =             clhs200 + clhs229*mu;
const double clhs231 =             DN(1,2)*N[3];
const double clhs232 =             DN(2,0)*N[1];
const double clhs233 =             DN(2,1)*N[1];
const double clhs234 =             DN(2,2)*N[1];
const double clhs235 =             clhs171*(clhs192 + clhs210 + clhs224);
const double clhs236 =             DN(3,0)*N[1];
const double clhs237 =             DN(3,1)*N[1];
const double clhs238 =             DN(3,2)*N[1];
const double clhs239 =             clhs171*(clhs199 + clhs215 + clhs229);
const double clhs240 =             pow(DN(2,0), 2);
const double clhs241 =             pow(N[2], 2)*clhs1;
const double clhs242 =             DN(2,0)*mu;
const double clhs243 =             DN(2,1)*clhs242;
const double clhs244 =             DN(2,2)*clhs242;
const double clhs245 =             DN(2,0)*N[2];
const double clhs246 =             DN(2,0)*DN(3,0);
const double clhs247 =             N[2]*N[3]*clhs1;
const double clhs248 =             clhs246*mu + clhs247;
const double clhs249 =             DN(3,1)*clhs242;
const double clhs250 =             DN(3,2)*clhs242;
const double clhs251 =             DN(2,0)*N[3];
const double clhs252 =             pow(DN(2,1), 2);
const double clhs253 =             DN(2,1)*mu;
const double clhs254 =             DN(2,2)*clhs253;
const double clhs255 =             DN(2,1)*N[2];
const double clhs256 =             DN(3,0)*clhs253;
const double clhs257 =             DN(2,1)*DN(3,1);
const double clhs258 =             clhs247 + clhs257*mu;
const double clhs259 =             DN(3,2)*clhs253;
const double clhs260 =             DN(2,1)*N[3];
const double clhs261 =             pow(DN(2,2), 2);
const double clhs262 =             DN(2,2)*N[2];
const double clhs263 =             DN(2,2)*mu;
const double clhs264 =             DN(3,0)*clhs263;
const double clhs265 =             DN(3,1)*clhs263;
const double clhs266 =             DN(2,2)*DN(3,2);
const double clhs267 =             clhs247 + clhs266*mu;
const double clhs268 =             DN(2,2)*N[3];
const double clhs269 =             DN(3,0)*N[2];
const double clhs270 =             DN(3,1)*N[2];
const double clhs271 =             DN(3,2)*N[2];
const double clhs272 =             clhs171*(clhs246 + clhs257 + clhs266);
const double clhs273 =             pow(DN(3,0), 2);
const double clhs274 =             pow(N[3], 2)*clhs1;
const double clhs275 =             DN(3,0)*mu;
const double clhs276 =             DN(3,1)*clhs275;
const double clhs277 =             DN(3,2)*clhs275;
const double clhs278 =             DN(3,0)*N[3];
const double clhs279 =             pow(DN(3,1), 2);
const double clhs280 =             DN(3,1)*DN(3,2)*mu;
const double clhs281 =             DN(3,1)*N[3];
const double clhs282 =             pow(DN(3,2), 2);
const double clhs283 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs3 + DN(0,1)*clhs5 + DN(0,2)*clhs7 + clhs0*mu + clhs2;
            lhs(0,1)=DN(0,0)*clhs10 + DN(0,1)*clhs12 + DN(0,2)*clhs15 + clhs9;
            lhs(0,2)=DN(0,0)*clhs17 + DN(0,1)*clhs19 + DN(0,2)*clhs21 + clhs16;
            lhs(0,3)=-clhs22;
            lhs(0,4)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs27 + clhs31;
            lhs(0,5)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + DN(0,2)*clhs38 + clhs32;
            lhs(0,6)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs44 + clhs39;
            lhs(0,7)=-clhs45;
            lhs(0,8)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs50 + clhs53;
            lhs(0,9)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs54;
            lhs(0,10)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs61;
            lhs(0,11)=-clhs67;
            lhs(0,12)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs75;
            lhs(0,13)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs82 + clhs76;
            lhs(0,14)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs88 + clhs83;
            lhs(0,15)=-clhs89;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs90 + DN(0,2)*clhs91 + clhs9;
            lhs(1,1)=DN(0,0)*clhs12 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs2 + clhs92*mu;
            lhs(1,2)=DN(0,0)*clhs19 + DN(0,1)*clhs98 + DN(0,2)*clhs100 + clhs97;
            lhs(1,3)=-clhs101;
            lhs(1,4)=DN(0,0)*clhs25 + DN(0,1)*clhs103 + DN(0,2)*clhs104 + clhs102;
            lhs(1,5)=DN(0,0)*clhs35 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + clhs109;
            lhs(1,6)=DN(0,0)*clhs42 + DN(0,1)*clhs111 + DN(0,2)*clhs113 + clhs110;
            lhs(1,7)=-clhs114;
            lhs(1,8)=DN(0,0)*clhs48 + DN(0,1)*clhs116 + DN(0,2)*clhs117 + clhs115;
            lhs(1,9)=DN(0,0)*clhs57 + DN(0,1)*clhs118 + DN(0,2)*clhs120 + clhs122;
            lhs(1,10)=DN(0,0)*clhs64 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs123;
            lhs(1,11)=-clhs127;
            lhs(1,12)=DN(0,0)*clhs70 + DN(0,1)*clhs129 + DN(0,2)*clhs130 + clhs128;
            lhs(1,13)=DN(0,0)*clhs79 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs135;
            lhs(1,14)=DN(0,0)*clhs86 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs136;
            lhs(1,15)=-clhs140;
            lhs(2,0)=DN(0,0)*clhs7 + DN(0,1)*clhs91 + DN(0,2)*clhs141 + clhs16;
            lhs(2,1)=DN(0,0)*clhs15 + DN(0,1)*clhs95 + DN(0,2)*clhs142 + clhs97;
            lhs(2,2)=DN(0,0)*clhs21 + DN(0,1)*clhs100 + DN(0,2)*clhs144 + clhs143*mu + clhs2;
            lhs(2,3)=-clhs145;
            lhs(2,4)=DN(0,0)*clhs27 + DN(0,1)*clhs104 + DN(0,2)*clhs148 + clhs147;
            lhs(2,5)=DN(0,0)*clhs38 + DN(0,1)*clhs107 + DN(0,2)*clhs150 + clhs149;
            lhs(2,6)=DN(0,0)*clhs44 + DN(0,1)*clhs113 + DN(0,2)*clhs151 + clhs153;
            lhs(2,7)=-clhs154;
            lhs(2,8)=DN(0,0)*clhs50 + DN(0,1)*clhs117 + DN(0,2)*clhs156 + clhs155;
            lhs(2,9)=DN(0,0)*clhs60 + DN(0,1)*clhs120 + DN(0,2)*clhs158 + clhs157;
            lhs(2,10)=DN(0,0)*clhs66 + DN(0,1)*clhs126 + DN(0,2)*clhs159 + clhs161;
            lhs(2,11)=-clhs162;
            lhs(2,12)=DN(0,0)*clhs72 + DN(0,1)*clhs130 + DN(0,2)*clhs164 + clhs163;
            lhs(2,13)=DN(0,0)*clhs82 + DN(0,1)*clhs133 + DN(0,2)*clhs166 + clhs165;
            lhs(2,14)=DN(0,0)*clhs88 + DN(0,1)*clhs139 + DN(0,2)*clhs167 + clhs169;
            lhs(2,15)=-clhs170;
            lhs(3,0)=clhs173*clhs22;
            lhs(3,1)=clhs101*clhs173;
            lhs(3,2)=clhs145*clhs173;
            lhs(3,3)=clhs171*(clhs0 + clhs143 + clhs92);
            lhs(3,4)=clhs172*clhs45 + clhs174;
            lhs(3,5)=clhs114*clhs172 + clhs175;
            lhs(3,6)=clhs154*clhs172 + clhs176;
            lhs(3,7)=clhs177;
            lhs(3,8)=clhs172*clhs67 + clhs178;
            lhs(3,9)=clhs127*clhs172 + clhs179;
            lhs(3,10)=clhs162*clhs172 + clhs180;
            lhs(3,11)=clhs181;
            lhs(3,12)=clhs172*clhs89 + clhs182;
            lhs(3,13)=clhs140*clhs172 + clhs183;
            lhs(3,14)=clhs170*clhs172 + clhs184;
            lhs(3,15)=clhs185;
            lhs(4,0)=DN(1,0)*clhs3 + DN(1,1)*clhs5 + DN(1,2)*clhs7 + clhs31;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs12 + DN(1,2)*clhs15 + clhs102;
            lhs(4,2)=DN(1,0)*clhs17 + DN(1,1)*clhs19 + DN(1,2)*clhs21 + clhs147;
            lhs(4,3)=-clhs174;
            lhs(4,4)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs27 + clhs186*mu + clhs187;
            lhs(4,5)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + DN(1,2)*clhs38 + clhs189;
            lhs(4,6)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs44 + clhs190;
            lhs(4,7)=-clhs191;
            lhs(4,8)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs50 + clhs195;
            lhs(4,9)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs196;
            lhs(4,10)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs197;
            lhs(4,11)=-clhs198;
            lhs(4,12)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs201;
            lhs(4,13)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs82 + clhs202;
            lhs(4,14)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs88 + clhs203;
            lhs(4,15)=-clhs204;
            lhs(5,0)=DN(1,0)*clhs5 + DN(1,1)*clhs90 + DN(1,2)*clhs91 + clhs32;
            lhs(5,1)=DN(1,0)*clhs12 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs109;
            lhs(5,2)=DN(1,0)*clhs19 + DN(1,1)*clhs98 + DN(1,2)*clhs100 + clhs149;
            lhs(5,3)=-clhs175;
            lhs(5,4)=DN(1,0)*clhs25 + DN(1,1)*clhs103 + DN(1,2)*clhs104 + clhs189;
            lhs(5,5)=DN(1,0)*clhs35 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + clhs187 + clhs205*mu;
            lhs(5,6)=DN(1,0)*clhs42 + DN(1,1)*clhs111 + DN(1,2)*clhs113 + clhs207;
            lhs(5,7)=-clhs208;
            lhs(5,8)=DN(1,0)*clhs48 + DN(1,1)*clhs116 + DN(1,2)*clhs117 + clhs209;
            lhs(5,9)=DN(1,0)*clhs57 + DN(1,1)*clhs118 + DN(1,2)*clhs120 + clhs211;
            lhs(5,10)=DN(1,0)*clhs64 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs212;
            lhs(5,11)=-clhs213;
            lhs(5,12)=DN(1,0)*clhs70 + DN(1,1)*clhs129 + DN(1,2)*clhs130 + clhs214;
            lhs(5,13)=DN(1,0)*clhs79 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs216;
            lhs(5,14)=DN(1,0)*clhs86 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs217;
            lhs(5,15)=-clhs218;
            lhs(6,0)=DN(1,0)*clhs7 + DN(1,1)*clhs91 + DN(1,2)*clhs141 + clhs39;
            lhs(6,1)=DN(1,0)*clhs15 + DN(1,1)*clhs95 + DN(1,2)*clhs142 + clhs110;
            lhs(6,2)=DN(1,0)*clhs21 + DN(1,1)*clhs100 + DN(1,2)*clhs144 + clhs153;
            lhs(6,3)=-clhs176;
            lhs(6,4)=DN(1,0)*clhs27 + DN(1,1)*clhs104 + DN(1,2)*clhs148 + clhs190;
            lhs(6,5)=DN(1,0)*clhs38 + DN(1,1)*clhs107 + DN(1,2)*clhs150 + clhs207;
            lhs(6,6)=DN(1,0)*clhs44 + DN(1,1)*clhs113 + DN(1,2)*clhs151 + clhs187 + clhs219*mu;
            lhs(6,7)=-clhs220;
            lhs(6,8)=DN(1,0)*clhs50 + DN(1,1)*clhs117 + DN(1,2)*clhs156 + clhs222;
            lhs(6,9)=DN(1,0)*clhs60 + DN(1,1)*clhs120 + DN(1,2)*clhs158 + clhs223;
            lhs(6,10)=DN(1,0)*clhs66 + DN(1,1)*clhs126 + DN(1,2)*clhs159 + clhs225;
            lhs(6,11)=-clhs226;
            lhs(6,12)=DN(1,0)*clhs72 + DN(1,1)*clhs130 + DN(1,2)*clhs164 + clhs227;
            lhs(6,13)=DN(1,0)*clhs82 + DN(1,1)*clhs133 + DN(1,2)*clhs166 + clhs228;
            lhs(6,14)=DN(1,0)*clhs88 + DN(1,1)*clhs139 + DN(1,2)*clhs167 + clhs230;
            lhs(6,15)=-clhs231;
            lhs(7,0)=clhs172*clhs174 + clhs45;
            lhs(7,1)=clhs114 + clhs172*clhs175;
            lhs(7,2)=clhs154 + clhs172*clhs176;
            lhs(7,3)=clhs177;
            lhs(7,4)=clhs173*clhs191;
            lhs(7,5)=clhs173*clhs208;
            lhs(7,6)=clhs173*clhs220;
            lhs(7,7)=clhs171*(clhs186 + clhs205 + clhs219);
            lhs(7,8)=clhs172*clhs198 + clhs232;
            lhs(7,9)=clhs172*clhs213 + clhs233;
            lhs(7,10)=clhs172*clhs226 + clhs234;
            lhs(7,11)=clhs235;
            lhs(7,12)=clhs172*clhs204 + clhs236;
            lhs(7,13)=clhs172*clhs218 + clhs237;
            lhs(7,14)=clhs172*clhs231 + clhs238;
            lhs(7,15)=clhs239;
            lhs(8,0)=DN(2,0)*clhs3 + DN(2,1)*clhs5 + DN(2,2)*clhs7 + clhs53;
            lhs(8,1)=DN(2,0)*clhs10 + DN(2,1)*clhs12 + DN(2,2)*clhs15 + clhs115;
            lhs(8,2)=DN(2,0)*clhs17 + DN(2,1)*clhs19 + DN(2,2)*clhs21 + clhs155;
            lhs(8,3)=-clhs178;
            lhs(8,4)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs27 + clhs195;
            lhs(8,5)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + DN(2,2)*clhs38 + clhs209;
            lhs(8,6)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs44 + clhs222;
            lhs(8,7)=-clhs232;
            lhs(8,8)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs50 + clhs240*mu + clhs241;
            lhs(8,9)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs243;
            lhs(8,10)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs244;
            lhs(8,11)=-clhs245;
            lhs(8,12)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs248;
            lhs(8,13)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs82 + clhs249;
            lhs(8,14)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs88 + clhs250;
            lhs(8,15)=-clhs251;
            lhs(9,0)=DN(2,0)*clhs5 + DN(2,1)*clhs90 + DN(2,2)*clhs91 + clhs54;
            lhs(9,1)=DN(2,0)*clhs12 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs122;
            lhs(9,2)=DN(2,0)*clhs19 + DN(2,1)*clhs98 + DN(2,2)*clhs100 + clhs157;
            lhs(9,3)=-clhs179;
            lhs(9,4)=DN(2,0)*clhs25 + DN(2,1)*clhs103 + DN(2,2)*clhs104 + clhs196;
            lhs(9,5)=DN(2,0)*clhs35 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs211;
            lhs(9,6)=DN(2,0)*clhs42 + DN(2,1)*clhs111 + DN(2,2)*clhs113 + clhs223;
            lhs(9,7)=-clhs233;
            lhs(9,8)=DN(2,0)*clhs48 + DN(2,1)*clhs116 + DN(2,2)*clhs117 + clhs243;
            lhs(9,9)=DN(2,0)*clhs57 + DN(2,1)*clhs118 + DN(2,2)*clhs120 + clhs241 + clhs252*mu;
            lhs(9,10)=DN(2,0)*clhs64 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs254;
            lhs(9,11)=-clhs255;
            lhs(9,12)=DN(2,0)*clhs70 + DN(2,1)*clhs129 + DN(2,2)*clhs130 + clhs256;
            lhs(9,13)=DN(2,0)*clhs79 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs258;
            lhs(9,14)=DN(2,0)*clhs86 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs259;
            lhs(9,15)=-clhs260;
            lhs(10,0)=DN(2,0)*clhs7 + DN(2,1)*clhs91 + DN(2,2)*clhs141 + clhs61;
            lhs(10,1)=DN(2,0)*clhs15 + DN(2,1)*clhs95 + DN(2,2)*clhs142 + clhs123;
            lhs(10,2)=DN(2,0)*clhs21 + DN(2,1)*clhs100 + DN(2,2)*clhs144 + clhs161;
            lhs(10,3)=-clhs180;
            lhs(10,4)=DN(2,0)*clhs27 + DN(2,1)*clhs104 + DN(2,2)*clhs148 + clhs197;
            lhs(10,5)=DN(2,0)*clhs38 + DN(2,1)*clhs107 + DN(2,2)*clhs150 + clhs212;
            lhs(10,6)=DN(2,0)*clhs44 + DN(2,1)*clhs113 + DN(2,2)*clhs151 + clhs225;
            lhs(10,7)=-clhs234;
            lhs(10,8)=DN(2,0)*clhs50 + DN(2,1)*clhs117 + DN(2,2)*clhs156 + clhs244;
            lhs(10,9)=DN(2,0)*clhs60 + DN(2,1)*clhs120 + DN(2,2)*clhs158 + clhs254;
            lhs(10,10)=DN(2,0)*clhs66 + DN(2,1)*clhs126 + DN(2,2)*clhs159 + clhs241 + clhs261*mu;
            lhs(10,11)=-clhs262;
            lhs(10,12)=DN(2,0)*clhs72 + DN(2,1)*clhs130 + DN(2,2)*clhs164 + clhs264;
            lhs(10,13)=DN(2,0)*clhs82 + DN(2,1)*clhs133 + DN(2,2)*clhs166 + clhs265;
            lhs(10,14)=DN(2,0)*clhs88 + DN(2,1)*clhs139 + DN(2,2)*clhs167 + clhs267;
            lhs(10,15)=-clhs268;
            lhs(11,0)=clhs172*clhs178 + clhs67;
            lhs(11,1)=clhs127 + clhs172*clhs179;
            lhs(11,2)=clhs162 + clhs172*clhs180;
            lhs(11,3)=clhs181;
            lhs(11,4)=clhs172*clhs232 + clhs198;
            lhs(11,5)=clhs172*clhs233 + clhs213;
            lhs(11,6)=clhs172*clhs234 + clhs226;
            lhs(11,7)=clhs235;
            lhs(11,8)=clhs173*clhs245;
            lhs(11,9)=clhs173*clhs255;
            lhs(11,10)=clhs173*clhs262;
            lhs(11,11)=clhs171*(clhs240 + clhs252 + clhs261);
            lhs(11,12)=clhs172*clhs251 + clhs269;
            lhs(11,13)=clhs172*clhs260 + clhs270;
            lhs(11,14)=clhs172*clhs268 + clhs271;
            lhs(11,15)=clhs272;
            lhs(12,0)=DN(3,0)*clhs3 + DN(3,1)*clhs5 + DN(3,2)*clhs7 + clhs75;
            lhs(12,1)=DN(3,0)*clhs10 + DN(3,1)*clhs12 + DN(3,2)*clhs15 + clhs128;
            lhs(12,2)=DN(3,0)*clhs17 + DN(3,1)*clhs19 + DN(3,2)*clhs21 + clhs163;
            lhs(12,3)=-clhs182;
            lhs(12,4)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs27 + clhs201;
            lhs(12,5)=DN(3,0)*clhs33 + DN(3,1)*clhs35 + DN(3,2)*clhs38 + clhs214;
            lhs(12,6)=DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs44 + clhs227;
            lhs(12,7)=-clhs236;
            lhs(12,8)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs50 + clhs248;
            lhs(12,9)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs256;
            lhs(12,10)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs264;
            lhs(12,11)=-clhs269;
            lhs(12,12)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs273*mu + clhs274;
            lhs(12,13)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs82 + clhs276;
            lhs(12,14)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs88 + clhs277;
            lhs(12,15)=-clhs278;
            lhs(13,0)=DN(3,0)*clhs5 + DN(3,1)*clhs90 + DN(3,2)*clhs91 + clhs76;
            lhs(13,1)=DN(3,0)*clhs12 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs135;
            lhs(13,2)=DN(3,0)*clhs19 + DN(3,1)*clhs98 + DN(3,2)*clhs100 + clhs165;
            lhs(13,3)=-clhs183;
            lhs(13,4)=DN(3,0)*clhs25 + DN(3,1)*clhs103 + DN(3,2)*clhs104 + clhs202;
            lhs(13,5)=DN(3,0)*clhs35 + DN(3,1)*clhs105 + DN(3,2)*clhs107 + clhs216;
            lhs(13,6)=DN(3,0)*clhs42 + DN(3,1)*clhs111 + DN(3,2)*clhs113 + clhs228;
            lhs(13,7)=-clhs237;
            lhs(13,8)=DN(3,0)*clhs48 + DN(3,1)*clhs116 + DN(3,2)*clhs117 + clhs249;
            lhs(13,9)=DN(3,0)*clhs57 + DN(3,1)*clhs118 + DN(3,2)*clhs120 + clhs258;
            lhs(13,10)=DN(3,0)*clhs64 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs265;
            lhs(13,11)=-clhs270;
            lhs(13,12)=DN(3,0)*clhs70 + DN(3,1)*clhs129 + DN(3,2)*clhs130 + clhs276;
            lhs(13,13)=DN(3,0)*clhs79 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs274 + clhs279*mu;
            lhs(13,14)=DN(3,0)*clhs86 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs280;
            lhs(13,15)=-clhs281;
            lhs(14,0)=DN(3,0)*clhs7 + DN(3,1)*clhs91 + DN(3,2)*clhs141 + clhs83;
            lhs(14,1)=DN(3,0)*clhs15 + DN(3,1)*clhs95 + DN(3,2)*clhs142 + clhs136;
            lhs(14,2)=DN(3,0)*clhs21 + DN(3,1)*clhs100 + DN(3,2)*clhs144 + clhs169;
            lhs(14,3)=-clhs184;
            lhs(14,4)=DN(3,0)*clhs27 + DN(3,1)*clhs104 + DN(3,2)*clhs148 + clhs203;
            lhs(14,5)=DN(3,0)*clhs38 + DN(3,1)*clhs107 + DN(3,2)*clhs150 + clhs217;
            lhs(14,6)=DN(3,0)*clhs44 + DN(3,1)*clhs113 + DN(3,2)*clhs151 + clhs230;
            lhs(14,7)=-clhs238;
            lhs(14,8)=DN(3,0)*clhs50 + DN(3,1)*clhs117 + DN(3,2)*clhs156 + clhs250;
            lhs(14,9)=DN(3,0)*clhs60 + DN(3,1)*clhs120 + DN(3,2)*clhs158 + clhs259;
            lhs(14,10)=DN(3,0)*clhs66 + DN(3,1)*clhs126 + DN(3,2)*clhs159 + clhs267;
            lhs(14,11)=-clhs271;
            lhs(14,12)=DN(3,0)*clhs72 + DN(3,1)*clhs130 + DN(3,2)*clhs164 + clhs277;
            lhs(14,13)=DN(3,0)*clhs82 + DN(3,1)*clhs133 + DN(3,2)*clhs166 + clhs280;
            lhs(14,14)=DN(3,0)*clhs88 + DN(3,1)*clhs139 + DN(3,2)*clhs167 + clhs274 + clhs282*mu;
            lhs(14,15)=-clhs283;
            lhs(15,0)=clhs172*clhs182 + clhs89;
            lhs(15,1)=clhs140 + clhs172*clhs183;
            lhs(15,2)=clhs170 + clhs172*clhs184;
            lhs(15,3)=clhs185;
            lhs(15,4)=clhs172*clhs236 + clhs204;
            lhs(15,5)=clhs172*clhs237 + clhs218;
            lhs(15,6)=clhs172*clhs238 + clhs231;
            lhs(15,7)=clhs239;
            lhs(15,8)=clhs172*clhs269 + clhs251;
            lhs(15,9)=clhs172*clhs270 + clhs260;
            lhs(15,10)=clhs172*clhs271 + clhs268;
            lhs(15,11)=clhs272;
            lhs(15,12)=clhs173*clhs278;
            lhs(15,13)=clhs173*clhs281;
            lhs(15,14)=clhs173*clhs283;
            lhs(15,15)=clhs171*(clhs273 + clhs279 + clhs282);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void Stokes<SymbolicNavierStokesData<3,8>>::ComputeGaussPointLHSContribution(
    SymbolicNavierStokesData<3,8> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &lhs = rData.lhs;

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             bdf0*rho;
const double clhs2 =             pow(N[0], 2)*clhs1;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs4 =             C(0,3)*DN(0,0);
const double clhs5 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs4;
const double clhs6 =             C(0,5)*DN(0,0);
const double clhs7 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs6;
const double clhs8 =             DN(0,0)*mu;
const double clhs9 =             DN(0,1)*clhs8;
const double clhs10 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs4;
const double clhs11 =             C(1,3)*DN(0,1);
const double clhs12 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs11;
const double clhs13 =             C(3,5)*DN(0,0);
const double clhs14 =             C(4,5)*DN(0,2);
const double clhs15 =             C(1,5)*DN(0,1) + clhs13 + clhs14;
const double clhs16 =             DN(0,2)*clhs8;
const double clhs17 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs6;
const double clhs18 =             C(3,4)*DN(0,1);
const double clhs19 =             C(2,3)*DN(0,2) + clhs13 + clhs18;
const double clhs20 =             C(2,5)*DN(0,2);
const double clhs21 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs20;
const double clhs22 =             DN(0,0)*N[0];
const double clhs23 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs24 =             C(0,3)*DN(1,0);
const double clhs25 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs24;
const double clhs26 =             C(0,5)*DN(1,0);
const double clhs27 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs26;
const double clhs28 =             DN(0,0)*DN(1,0);
const double clhs29 =             N[0]*clhs1;
const double clhs30 =             N[1]*clhs29;
const double clhs31 =             clhs28*mu + clhs30;
const double clhs32 =             DN(1,1)*clhs8;
const double clhs33 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs24;
const double clhs34 =             C(1,3)*DN(1,1);
const double clhs35 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs34;
const double clhs36 =             C(3,5)*DN(1,0);
const double clhs37 =             C(4,5)*DN(1,2);
const double clhs38 =             C(1,5)*DN(1,1) + clhs36 + clhs37;
const double clhs39 =             DN(1,2)*clhs8;
const double clhs40 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs26;
const double clhs41 =             C(3,4)*DN(1,1);
const double clhs42 =             C(2,3)*DN(1,2) + clhs36 + clhs41;
const double clhs43 =             C(2,5)*DN(1,2);
const double clhs44 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs43;
const double clhs45 =             DN(0,0)*N[1];
const double clhs46 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs47 =             C(0,3)*DN(2,0);
const double clhs48 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs47;
const double clhs49 =             C(0,5)*DN(2,0);
const double clhs50 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs49;
const double clhs51 =             DN(0,0)*DN(2,0);
const double clhs52 =             N[2]*clhs29;
const double clhs53 =             clhs51*mu + clhs52;
const double clhs54 =             DN(2,1)*clhs8;
const double clhs55 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs47;
const double clhs56 =             C(1,3)*DN(2,1);
const double clhs57 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs56;
const double clhs58 =             C(3,5)*DN(2,0);
const double clhs59 =             C(4,5)*DN(2,2);
const double clhs60 =             C(1,5)*DN(2,1) + clhs58 + clhs59;
const double clhs61 =             DN(2,2)*clhs8;
const double clhs62 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs49;
const double clhs63 =             C(3,4)*DN(2,1);
const double clhs64 =             C(2,3)*DN(2,2) + clhs58 + clhs63;
const double clhs65 =             C(2,5)*DN(2,2);
const double clhs66 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs65;
const double clhs67 =             DN(0,0)*N[2];
const double clhs68 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs69 =             C(0,3)*DN(3,0);
const double clhs70 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs69;
const double clhs71 =             C(0,5)*DN(3,0);
const double clhs72 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs71;
const double clhs73 =             DN(0,0)*DN(3,0);
const double clhs74 =             N[3]*clhs29;
const double clhs75 =             clhs73*mu + clhs74;
const double clhs76 =             DN(3,1)*clhs8;
const double clhs77 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs69;
const double clhs78 =             C(1,3)*DN(3,1);
const double clhs79 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs78;
const double clhs80 =             C(3,5)*DN(3,0);
const double clhs81 =             C(4,5)*DN(3,2);
const double clhs82 =             C(1,5)*DN(3,1) + clhs80 + clhs81;
const double clhs83 =             DN(3,2)*clhs8;
const double clhs84 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs71;
const double clhs85 =             C(3,4)*DN(3,1);
const double clhs86 =             C(2,3)*DN(3,2) + clhs80 + clhs85;
const double clhs87 =             C(2,5)*DN(3,2);
const double clhs88 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs87;
const double clhs89 =             DN(0,0)*N[3];
const double clhs90 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs11;
const double clhs91 =             C(0,4)*DN(0,0) + clhs14 + clhs18;
const double clhs92 =             pow(DN(0,1), 2);
const double clhs93 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs94 =             C(1,4)*DN(0,1);
const double clhs95 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs94;
const double clhs96 =             DN(0,1)*mu;
const double clhs97 =             DN(0,2)*clhs96;
const double clhs98 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs94;
const double clhs99 =             C(2,4)*DN(0,2);
const double clhs100 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs99;
const double clhs101 =             DN(0,1)*N[0];
const double clhs102 =             DN(1,0)*clhs96;
const double clhs103 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs34;
const double clhs104 =             C(0,4)*DN(1,0) + clhs37 + clhs41;
const double clhs105 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs106 =             C(1,4)*DN(1,1);
const double clhs107 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs106;
const double clhs108 =             DN(0,1)*DN(1,1);
const double clhs109 =             clhs108*mu + clhs30;
const double clhs110 =             DN(1,2)*clhs96;
const double clhs111 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs106;
const double clhs112 =             C(2,4)*DN(1,2);
const double clhs113 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs112;
const double clhs114 =             DN(0,1)*N[1];
const double clhs115 =             DN(2,0)*clhs96;
const double clhs116 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs56;
const double clhs117 =             C(0,4)*DN(2,0) + clhs59 + clhs63;
const double clhs118 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs119 =             C(1,4)*DN(2,1);
const double clhs120 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs119;
const double clhs121 =             DN(0,1)*DN(2,1);
const double clhs122 =             clhs121*mu + clhs52;
const double clhs123 =             DN(2,2)*clhs96;
const double clhs124 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs119;
const double clhs125 =             C(2,4)*DN(2,2);
const double clhs126 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs125;
const double clhs127 =             DN(0,1)*N[2];
const double clhs128 =             DN(3,0)*clhs96;
const double clhs129 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs78;
const double clhs130 =             C(0,4)*DN(3,0) + clhs81 + clhs85;
const double clhs131 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs132 =             C(1,4)*DN(3,1);
const double clhs133 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs132;
const double clhs134 =             DN(0,1)*DN(3,1);
const double clhs135 =             clhs134*mu + clhs74;
const double clhs136 =             DN(3,2)*clhs96;
const double clhs137 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs132;
const double clhs138 =             C(2,4)*DN(3,2);
const double clhs139 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs138;
const double clhs140 =             DN(0,1)*N[3];
const double clhs141 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs20;
const double clhs142 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs99;
const double clhs143 =             pow(DN(0,2), 2);
const double clhs144 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs145 =             DN(0,2)*N[0];
const double clhs146 =             DN(0,2)*mu;
const double clhs147 =             DN(1,0)*clhs146;
const double clhs148 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs43;
const double clhs149 =             DN(1,1)*clhs146;
const double clhs150 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs112;
const double clhs151 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs152 =             DN(0,2)*DN(1,2);
const double clhs153 =             clhs152*mu + clhs30;
const double clhs154 =             DN(0,2)*N[1];
const double clhs155 =             DN(2,0)*clhs146;
const double clhs156 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs65;
const double clhs157 =             DN(2,1)*clhs146;
const double clhs158 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs125;
const double clhs159 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs160 =             DN(0,2)*DN(2,2);
const double clhs161 =             clhs160*mu + clhs52;
const double clhs162 =             DN(0,2)*N[2];
const double clhs163 =             DN(3,0)*clhs146;
const double clhs164 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs87;
const double clhs165 =             DN(3,1)*clhs146;
const double clhs166 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs138;
const double clhs167 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs168 =             DN(0,2)*DN(3,2);
const double clhs169 =             clhs168*mu + clhs74;
const double clhs170 =             DN(0,2)*N[3];
const double clhs171 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs172 =             clhs1*clhs171;
const double clhs173 =             clhs172 + 1;
const double clhs174 =             DN(1,0)*N[0];
const double clhs175 =             DN(1,1)*N[0];
const double clhs176 =             DN(1,2)*N[0];
const double clhs177 =             clhs171*(clhs108 + clhs152 + clhs28);
const double clhs178 =             DN(2,0)*N[0];
const double clhs179 =             DN(2,1)*N[0];
const double clhs180 =             DN(2,2)*N[0];
const double clhs181 =             clhs171*(clhs121 + clhs160 + clhs51);
const double clhs182 =             DN(3,0)*N[0];
const double clhs183 =             DN(3,1)*N[0];
const double clhs184 =             DN(3,2)*N[0];
const double clhs185 =             clhs171*(clhs134 + clhs168 + clhs73);
const double clhs186 =             pow(DN(1,0), 2);
const double clhs187 =             pow(N[1], 2)*clhs1;
const double clhs188 =             DN(1,0)*mu;
const double clhs189 =             DN(1,1)*clhs188;
const double clhs190 =             DN(1,2)*clhs188;
const double clhs191 =             DN(1,0)*N[1];
const double clhs192 =             DN(1,0)*DN(2,0);
const double clhs193 =             N[1]*clhs1;
const double clhs194 =             N[2]*clhs193;
const double clhs195 =             clhs192*mu + clhs194;
const double clhs196 =             DN(2,1)*clhs188;
const double clhs197 =             DN(2,2)*clhs188;
const double clhs198 =             DN(1,0)*N[2];
const double clhs199 =             DN(1,0)*DN(3,0);
const double clhs200 =             N[3]*clhs193;
const double clhs201 =             clhs199*mu + clhs200;
const double clhs202 =             DN(3,1)*clhs188;
const double clhs203 =             DN(3,2)*clhs188;
const double clhs204 =             DN(1,0)*N[3];
const double clhs205 =             pow(DN(1,1), 2);
const double clhs206 =             DN(1,1)*mu;
const double clhs207 =             DN(1,2)*clhs206;
const double clhs208 =             DN(1,1)*N[1];
const double clhs209 =             DN(2,0)*clhs206;
const double clhs210 =             DN(1,1)*DN(2,1);
const double clhs211 =             clhs194 + clhs210*mu;
const double clhs212 =             DN(2,2)*clhs206;
const double clhs213 =             DN(1,1)*N[2];
const double clhs214 =             DN(3,0)*clhs206;
const double clhs215 =             DN(1,1)*DN(3,1);
const double clhs216 =             clhs200 + clhs215*mu;
const double clhs217 =             DN(3,2)*clhs206;
const double clhs218 =             DN(1,1)*N[3];
const double clhs219 =             pow(DN(1,2), 2);
const double clhs220 =             DN(1,2)*N[1];
const double clhs221 =             DN(1,2)*mu;
const double clhs222 =             DN(2,0)*clhs221;
const double clhs223 =             DN(2,1)*clhs221;
const double clhs224 =             DN(1,2)*DN(2,2);
const double clhs225 =             clhs194 + clhs224*mu;
const double clhs226 =             DN(1,2)*N[2];
const double clhs227 =             DN(3,0)*clhs221;
const double clhs228 =             DN(3,1)*clhs221;
const double clhs229 =             DN(1,2)*DN(3,2);
const double clhs230 =             clhs200 + clhs229*mu;
const double clhs231 =             DN(1,2)*N[3];
const double clhs232 =             DN(2,0)*N[1];
const double clhs233 =             DN(2,1)*N[1];
const double clhs234 =             DN(2,2)*N[1];
const double clhs235 =             clhs171*(clhs192 + clhs210 + clhs224);
const double clhs236 =             DN(3,0)*N[1];
const double clhs237 =             DN(3,1)*N[1];
const double clhs238 =             DN(3,2)*N[1];
const double clhs239 =             clhs171*(clhs199 + clhs215 + clhs229);
const double clhs240 =             pow(DN(2,0), 2);
const double clhs241 =             pow(N[2], 2)*clhs1;
const double clhs242 =             DN(2,0)*mu;
const double clhs243 =             DN(2,1)*clhs242;
const double clhs244 =             DN(2,2)*clhs242;
const double clhs245 =             DN(2,0)*N[2];
const double clhs246 =             DN(2,0)*DN(3,0);
const double clhs247 =             N[2]*N[3]*clhs1;
const double clhs248 =             clhs246*mu + clhs247;
const double clhs249 =             DN(3,1)*clhs242;
const double clhs250 =             DN(3,2)*clhs242;
const double clhs251 =             DN(2,0)*N[3];
const double clhs252 =             pow(DN(2,1), 2);
const double clhs253 =             DN(2,1)*mu;
const double clhs254 =             DN(2,2)*clhs253;
const double clhs255 =             DN(2,1)*N[2];
const double clhs256 =             DN(3,0)*clhs253;
const double clhs257 =             DN(2,1)*DN(3,1);
const double clhs258 =             clhs247 + clhs257*mu;
const double clhs259 =             DN(3,2)*clhs253;
const double clhs260 =             DN(2,1)*N[3];
const double clhs261 =             pow(DN(2,2), 2);
const double clhs262 =             DN(2,2)*N[2];
const double clhs263 =             DN(2,2)*mu;
const double clhs264 =             DN(3,0)*clhs263;
const double clhs265 =             DN(3,1)*clhs263;
const double clhs266 =             DN(2,2)*DN(3,2);
const double clhs267 =             clhs247 + clhs266*mu;
const double clhs268 =             DN(2,2)*N[3];
const double clhs269 =             DN(3,0)*N[2];
const double clhs270 =             DN(3,1)*N[2];
const double clhs271 =             DN(3,2)*N[2];
const double clhs272 =             clhs171*(clhs246 + clhs257 + clhs266);
const double clhs273 =             pow(DN(3,0), 2);
const double clhs274 =             pow(N[3], 2)*clhs1;
const double clhs275 =             DN(3,0)*mu;
const double clhs276 =             DN(3,1)*clhs275;
const double clhs277 =             DN(3,2)*clhs275;
const double clhs278 =             DN(3,0)*N[3];
const double clhs279 =             pow(DN(3,1), 2);
const double clhs280 =             DN(3,1)*DN(3,2)*mu;
const double clhs281 =             DN(3,1)*N[3];
const double clhs282 =             pow(DN(3,2), 2);
const double clhs283 =             DN(3,2)*N[3];
            lhs(0,0)=DN(0,0)*clhs3 + DN(0,1)*clhs5 + DN(0,2)*clhs7 + clhs0*mu + clhs2;
            lhs(0,1)=DN(0,0)*clhs10 + DN(0,1)*clhs12 + DN(0,2)*clhs15 + clhs9;
            lhs(0,2)=DN(0,0)*clhs17 + DN(0,1)*clhs19 + DN(0,2)*clhs21 + clhs16;
            lhs(0,3)=-clhs22;
            lhs(0,4)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs27 + clhs31;
            lhs(0,5)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + DN(0,2)*clhs38 + clhs32;
            lhs(0,6)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs44 + clhs39;
            lhs(0,7)=-clhs45;
            lhs(0,8)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs50 + clhs53;
            lhs(0,9)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs54;
            lhs(0,10)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs61;
            lhs(0,11)=-clhs67;
            lhs(0,12)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs75;
            lhs(0,13)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs82 + clhs76;
            lhs(0,14)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs88 + clhs83;
            lhs(0,15)=-clhs89;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs90 + DN(0,2)*clhs91 + clhs9;
            lhs(1,1)=DN(0,0)*clhs12 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs2 + clhs92*mu;
            lhs(1,2)=DN(0,0)*clhs19 + DN(0,1)*clhs98 + DN(0,2)*clhs100 + clhs97;
            lhs(1,3)=-clhs101;
            lhs(1,4)=DN(0,0)*clhs25 + DN(0,1)*clhs103 + DN(0,2)*clhs104 + clhs102;
            lhs(1,5)=DN(0,0)*clhs35 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + clhs109;
            lhs(1,6)=DN(0,0)*clhs42 + DN(0,1)*clhs111 + DN(0,2)*clhs113 + clhs110;
            lhs(1,7)=-clhs114;
            lhs(1,8)=DN(0,0)*clhs48 + DN(0,1)*clhs116 + DN(0,2)*clhs117 + clhs115;
            lhs(1,9)=DN(0,0)*clhs57 + DN(0,1)*clhs118 + DN(0,2)*clhs120 + clhs122;
            lhs(1,10)=DN(0,0)*clhs64 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs123;
            lhs(1,11)=-clhs127;
            lhs(1,12)=DN(0,0)*clhs70 + DN(0,1)*clhs129 + DN(0,2)*clhs130 + clhs128;
            lhs(1,13)=DN(0,0)*clhs79 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs135;
            lhs(1,14)=DN(0,0)*clhs86 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs136;
            lhs(1,15)=-clhs140;
            lhs(2,0)=DN(0,0)*clhs7 + DN(0,1)*clhs91 + DN(0,2)*clhs141 + clhs16;
            lhs(2,1)=DN(0,0)*clhs15 + DN(0,1)*clhs95 + DN(0,2)*clhs142 + clhs97;
            lhs(2,2)=DN(0,0)*clhs21 + DN(0,1)*clhs100 + DN(0,2)*clhs144 + clhs143*mu + clhs2;
            lhs(2,3)=-clhs145;
            lhs(2,4)=DN(0,0)*clhs27 + DN(0,1)*clhs104 + DN(0,2)*clhs148 + clhs147;
            lhs(2,5)=DN(0,0)*clhs38 + DN(0,1)*clhs107 + DN(0,2)*clhs150 + clhs149;
            lhs(2,6)=DN(0,0)*clhs44 + DN(0,1)*clhs113 + DN(0,2)*clhs151 + clhs153;
            lhs(2,7)=-clhs154;
            lhs(2,8)=DN(0,0)*clhs50 + DN(0,1)*clhs117 + DN(0,2)*clhs156 + clhs155;
            lhs(2,9)=DN(0,0)*clhs60 + DN(0,1)*clhs120 + DN(0,2)*clhs158 + clhs157;
            lhs(2,10)=DN(0,0)*clhs66 + DN(0,1)*clhs126 + DN(0,2)*clhs159 + clhs161;
            lhs(2,11)=-clhs162;
            lhs(2,12)=DN(0,0)*clhs72 + DN(0,1)*clhs130 + DN(0,2)*clhs164 + clhs163;
            lhs(2,13)=DN(0,0)*clhs82 + DN(0,1)*clhs133 + DN(0,2)*clhs166 + clhs165;
            lhs(2,14)=DN(0,0)*clhs88 + DN(0,1)*clhs139 + DN(0,2)*clhs167 + clhs169;
            lhs(2,15)=-clhs170;
            lhs(3,0)=clhs173*clhs22;
            lhs(3,1)=clhs101*clhs173;
            lhs(3,2)=clhs145*clhs173;
            lhs(3,3)=clhs171*(clhs0 + clhs143 + clhs92);
            lhs(3,4)=clhs172*clhs45 + clhs174;
            lhs(3,5)=clhs114*clhs172 + clhs175;
            lhs(3,6)=clhs154*clhs172 + clhs176;
            lhs(3,7)=clhs177;
            lhs(3,8)=clhs172*clhs67 + clhs178;
            lhs(3,9)=clhs127*clhs172 + clhs179;
            lhs(3,10)=clhs162*clhs172 + clhs180;
            lhs(3,11)=clhs181;
            lhs(3,12)=clhs172*clhs89 + clhs182;
            lhs(3,13)=clhs140*clhs172 + clhs183;
            lhs(3,14)=clhs170*clhs172 + clhs184;
            lhs(3,15)=clhs185;
            lhs(4,0)=DN(1,0)*clhs3 + DN(1,1)*clhs5 + DN(1,2)*clhs7 + clhs31;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,1)*clhs12 + DN(1,2)*clhs15 + clhs102;
            lhs(4,2)=DN(1,0)*clhs17 + DN(1,1)*clhs19 + DN(1,2)*clhs21 + clhs147;
            lhs(4,3)=-clhs174;
            lhs(4,4)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs27 + clhs186*mu + clhs187;
            lhs(4,5)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + DN(1,2)*clhs38 + clhs189;
            lhs(4,6)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs44 + clhs190;
            lhs(4,7)=-clhs191;
            lhs(4,8)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs50 + clhs195;
            lhs(4,9)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs196;
            lhs(4,10)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs197;
            lhs(4,11)=-clhs198;
            lhs(4,12)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs201;
            lhs(4,13)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs82 + clhs202;
            lhs(4,14)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs88 + clhs203;
            lhs(4,15)=-clhs204;
            lhs(5,0)=DN(1,0)*clhs5 + DN(1,1)*clhs90 + DN(1,2)*clhs91 + clhs32;
            lhs(5,1)=DN(1,0)*clhs12 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs109;
            lhs(5,2)=DN(1,0)*clhs19 + DN(1,1)*clhs98 + DN(1,2)*clhs100 + clhs149;
            lhs(5,3)=-clhs175;
            lhs(5,4)=DN(1,0)*clhs25 + DN(1,1)*clhs103 + DN(1,2)*clhs104 + clhs189;
            lhs(5,5)=DN(1,0)*clhs35 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + clhs187 + clhs205*mu;
            lhs(5,6)=DN(1,0)*clhs42 + DN(1,1)*clhs111 + DN(1,2)*clhs113 + clhs207;
            lhs(5,7)=-clhs208;
            lhs(5,8)=DN(1,0)*clhs48 + DN(1,1)*clhs116 + DN(1,2)*clhs117 + clhs209;
            lhs(5,9)=DN(1,0)*clhs57 + DN(1,1)*clhs118 + DN(1,2)*clhs120 + clhs211;
            lhs(5,10)=DN(1,0)*clhs64 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs212;
            lhs(5,11)=-clhs213;
            lhs(5,12)=DN(1,0)*clhs70 + DN(1,1)*clhs129 + DN(1,2)*clhs130 + clhs214;
            lhs(5,13)=DN(1,0)*clhs79 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs216;
            lhs(5,14)=DN(1,0)*clhs86 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs217;
            lhs(5,15)=-clhs218;
            lhs(6,0)=DN(1,0)*clhs7 + DN(1,1)*clhs91 + DN(1,2)*clhs141 + clhs39;
            lhs(6,1)=DN(1,0)*clhs15 + DN(1,1)*clhs95 + DN(1,2)*clhs142 + clhs110;
            lhs(6,2)=DN(1,0)*clhs21 + DN(1,1)*clhs100 + DN(1,2)*clhs144 + clhs153;
            lhs(6,3)=-clhs176;
            lhs(6,4)=DN(1,0)*clhs27 + DN(1,1)*clhs104 + DN(1,2)*clhs148 + clhs190;
            lhs(6,5)=DN(1,0)*clhs38 + DN(1,1)*clhs107 + DN(1,2)*clhs150 + clhs207;
            lhs(6,6)=DN(1,0)*clhs44 + DN(1,1)*clhs113 + DN(1,2)*clhs151 + clhs187 + clhs219*mu;
            lhs(6,7)=-clhs220;
            lhs(6,8)=DN(1,0)*clhs50 + DN(1,1)*clhs117 + DN(1,2)*clhs156 + clhs222;
            lhs(6,9)=DN(1,0)*clhs60 + DN(1,1)*clhs120 + DN(1,2)*clhs158 + clhs223;
            lhs(6,10)=DN(1,0)*clhs66 + DN(1,1)*clhs126 + DN(1,2)*clhs159 + clhs225;
            lhs(6,11)=-clhs226;
            lhs(6,12)=DN(1,0)*clhs72 + DN(1,1)*clhs130 + DN(1,2)*clhs164 + clhs227;
            lhs(6,13)=DN(1,0)*clhs82 + DN(1,1)*clhs133 + DN(1,2)*clhs166 + clhs228;
            lhs(6,14)=DN(1,0)*clhs88 + DN(1,1)*clhs139 + DN(1,2)*clhs167 + clhs230;
            lhs(6,15)=-clhs231;
            lhs(7,0)=clhs172*clhs174 + clhs45;
            lhs(7,1)=clhs114 + clhs172*clhs175;
            lhs(7,2)=clhs154 + clhs172*clhs176;
            lhs(7,3)=clhs177;
            lhs(7,4)=clhs173*clhs191;
            lhs(7,5)=clhs173*clhs208;
            lhs(7,6)=clhs173*clhs220;
            lhs(7,7)=clhs171*(clhs186 + clhs205 + clhs219);
            lhs(7,8)=clhs172*clhs198 + clhs232;
            lhs(7,9)=clhs172*clhs213 + clhs233;
            lhs(7,10)=clhs172*clhs226 + clhs234;
            lhs(7,11)=clhs235;
            lhs(7,12)=clhs172*clhs204 + clhs236;
            lhs(7,13)=clhs172*clhs218 + clhs237;
            lhs(7,14)=clhs172*clhs231 + clhs238;
            lhs(7,15)=clhs239;
            lhs(8,0)=DN(2,0)*clhs3 + DN(2,1)*clhs5 + DN(2,2)*clhs7 + clhs53;
            lhs(8,1)=DN(2,0)*clhs10 + DN(2,1)*clhs12 + DN(2,2)*clhs15 + clhs115;
            lhs(8,2)=DN(2,0)*clhs17 + DN(2,1)*clhs19 + DN(2,2)*clhs21 + clhs155;
            lhs(8,3)=-clhs178;
            lhs(8,4)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs27 + clhs195;
            lhs(8,5)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + DN(2,2)*clhs38 + clhs209;
            lhs(8,6)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs44 + clhs222;
            lhs(8,7)=-clhs232;
            lhs(8,8)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs50 + clhs240*mu + clhs241;
            lhs(8,9)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs243;
            lhs(8,10)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs244;
            lhs(8,11)=-clhs245;
            lhs(8,12)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs248;
            lhs(8,13)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs82 + clhs249;
            lhs(8,14)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs88 + clhs250;
            lhs(8,15)=-clhs251;
            lhs(9,0)=DN(2,0)*clhs5 + DN(2,1)*clhs90 + DN(2,2)*clhs91 + clhs54;
            lhs(9,1)=DN(2,0)*clhs12 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs122;
            lhs(9,2)=DN(2,0)*clhs19 + DN(2,1)*clhs98 + DN(2,2)*clhs100 + clhs157;
            lhs(9,3)=-clhs179;
            lhs(9,4)=DN(2,0)*clhs25 + DN(2,1)*clhs103 + DN(2,2)*clhs104 + clhs196;
            lhs(9,5)=DN(2,0)*clhs35 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs211;
            lhs(9,6)=DN(2,0)*clhs42 + DN(2,1)*clhs111 + DN(2,2)*clhs113 + clhs223;
            lhs(9,7)=-clhs233;
            lhs(9,8)=DN(2,0)*clhs48 + DN(2,1)*clhs116 + DN(2,2)*clhs117 + clhs243;
            lhs(9,9)=DN(2,0)*clhs57 + DN(2,1)*clhs118 + DN(2,2)*clhs120 + clhs241 + clhs252*mu;
            lhs(9,10)=DN(2,0)*clhs64 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs254;
            lhs(9,11)=-clhs255;
            lhs(9,12)=DN(2,0)*clhs70 + DN(2,1)*clhs129 + DN(2,2)*clhs130 + clhs256;
            lhs(9,13)=DN(2,0)*clhs79 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs258;
            lhs(9,14)=DN(2,0)*clhs86 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs259;
            lhs(9,15)=-clhs260;
            lhs(10,0)=DN(2,0)*clhs7 + DN(2,1)*clhs91 + DN(2,2)*clhs141 + clhs61;
            lhs(10,1)=DN(2,0)*clhs15 + DN(2,1)*clhs95 + DN(2,2)*clhs142 + clhs123;
            lhs(10,2)=DN(2,0)*clhs21 + DN(2,1)*clhs100 + DN(2,2)*clhs144 + clhs161;
            lhs(10,3)=-clhs180;
            lhs(10,4)=DN(2,0)*clhs27 + DN(2,1)*clhs104 + DN(2,2)*clhs148 + clhs197;
            lhs(10,5)=DN(2,0)*clhs38 + DN(2,1)*clhs107 + DN(2,2)*clhs150 + clhs212;
            lhs(10,6)=DN(2,0)*clhs44 + DN(2,1)*clhs113 + DN(2,2)*clhs151 + clhs225;
            lhs(10,7)=-clhs234;
            lhs(10,8)=DN(2,0)*clhs50 + DN(2,1)*clhs117 + DN(2,2)*clhs156 + clhs244;
            lhs(10,9)=DN(2,0)*clhs60 + DN(2,1)*clhs120 + DN(2,2)*clhs158 + clhs254;
            lhs(10,10)=DN(2,0)*clhs66 + DN(2,1)*clhs126 + DN(2,2)*clhs159 + clhs241 + clhs261*mu;
            lhs(10,11)=-clhs262;
            lhs(10,12)=DN(2,0)*clhs72 + DN(2,1)*clhs130 + DN(2,2)*clhs164 + clhs264;
            lhs(10,13)=DN(2,0)*clhs82 + DN(2,1)*clhs133 + DN(2,2)*clhs166 + clhs265;
            lhs(10,14)=DN(2,0)*clhs88 + DN(2,1)*clhs139 + DN(2,2)*clhs167 + clhs267;
            lhs(10,15)=-clhs268;
            lhs(11,0)=clhs172*clhs178 + clhs67;
            lhs(11,1)=clhs127 + clhs172*clhs179;
            lhs(11,2)=clhs162 + clhs172*clhs180;
            lhs(11,3)=clhs181;
            lhs(11,4)=clhs172*clhs232 + clhs198;
            lhs(11,5)=clhs172*clhs233 + clhs213;
            lhs(11,6)=clhs172*clhs234 + clhs226;
            lhs(11,7)=clhs235;
            lhs(11,8)=clhs173*clhs245;
            lhs(11,9)=clhs173*clhs255;
            lhs(11,10)=clhs173*clhs262;
            lhs(11,11)=clhs171*(clhs240 + clhs252 + clhs261);
            lhs(11,12)=clhs172*clhs251 + clhs269;
            lhs(11,13)=clhs172*clhs260 + clhs270;
            lhs(11,14)=clhs172*clhs268 + clhs271;
            lhs(11,15)=clhs272;
            lhs(12,0)=DN(3,0)*clhs3 + DN(3,1)*clhs5 + DN(3,2)*clhs7 + clhs75;
            lhs(12,1)=DN(3,0)*clhs10 + DN(3,1)*clhs12 + DN(3,2)*clhs15 + clhs128;
            lhs(12,2)=DN(3,0)*clhs17 + DN(3,1)*clhs19 + DN(3,2)*clhs21 + clhs163;
            lhs(12,3)=-clhs182;
            lhs(12,4)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs27 + clhs201;
            lhs(12,5)=DN(3,0)*clhs33 + DN(3,1)*clhs35 + DN(3,2)*clhs38 + clhs214;
            lhs(12,6)=DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs44 + clhs227;
            lhs(12,7)=-clhs236;
            lhs(12,8)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs50 + clhs248;
            lhs(12,9)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs256;
            lhs(12,10)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs264;
            lhs(12,11)=-clhs269;
            lhs(12,12)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs273*mu + clhs274;
            lhs(12,13)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs82 + clhs276;
            lhs(12,14)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs88 + clhs277;
            lhs(12,15)=-clhs278;
            lhs(13,0)=DN(3,0)*clhs5 + DN(3,1)*clhs90 + DN(3,2)*clhs91 + clhs76;
            lhs(13,1)=DN(3,0)*clhs12 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs135;
            lhs(13,2)=DN(3,0)*clhs19 + DN(3,1)*clhs98 + DN(3,2)*clhs100 + clhs165;
            lhs(13,3)=-clhs183;
            lhs(13,4)=DN(3,0)*clhs25 + DN(3,1)*clhs103 + DN(3,2)*clhs104 + clhs202;
            lhs(13,5)=DN(3,0)*clhs35 + DN(3,1)*clhs105 + DN(3,2)*clhs107 + clhs216;
            lhs(13,6)=DN(3,0)*clhs42 + DN(3,1)*clhs111 + DN(3,2)*clhs113 + clhs228;
            lhs(13,7)=-clhs237;
            lhs(13,8)=DN(3,0)*clhs48 + DN(3,1)*clhs116 + DN(3,2)*clhs117 + clhs249;
            lhs(13,9)=DN(3,0)*clhs57 + DN(3,1)*clhs118 + DN(3,2)*clhs120 + clhs258;
            lhs(13,10)=DN(3,0)*clhs64 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs265;
            lhs(13,11)=-clhs270;
            lhs(13,12)=DN(3,0)*clhs70 + DN(3,1)*clhs129 + DN(3,2)*clhs130 + clhs276;
            lhs(13,13)=DN(3,0)*clhs79 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs274 + clhs279*mu;
            lhs(13,14)=DN(3,0)*clhs86 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs280;
            lhs(13,15)=-clhs281;
            lhs(14,0)=DN(3,0)*clhs7 + DN(3,1)*clhs91 + DN(3,2)*clhs141 + clhs83;
            lhs(14,1)=DN(3,0)*clhs15 + DN(3,1)*clhs95 + DN(3,2)*clhs142 + clhs136;
            lhs(14,2)=DN(3,0)*clhs21 + DN(3,1)*clhs100 + DN(3,2)*clhs144 + clhs169;
            lhs(14,3)=-clhs184;
            lhs(14,4)=DN(3,0)*clhs27 + DN(3,1)*clhs104 + DN(3,2)*clhs148 + clhs203;
            lhs(14,5)=DN(3,0)*clhs38 + DN(3,1)*clhs107 + DN(3,2)*clhs150 + clhs217;
            lhs(14,6)=DN(3,0)*clhs44 + DN(3,1)*clhs113 + DN(3,2)*clhs151 + clhs230;
            lhs(14,7)=-clhs238;
            lhs(14,8)=DN(3,0)*clhs50 + DN(3,1)*clhs117 + DN(3,2)*clhs156 + clhs250;
            lhs(14,9)=DN(3,0)*clhs60 + DN(3,1)*clhs120 + DN(3,2)*clhs158 + clhs259;
            lhs(14,10)=DN(3,0)*clhs66 + DN(3,1)*clhs126 + DN(3,2)*clhs159 + clhs267;
            lhs(14,11)=-clhs271;
            lhs(14,12)=DN(3,0)*clhs72 + DN(3,1)*clhs130 + DN(3,2)*clhs164 + clhs277;
            lhs(14,13)=DN(3,0)*clhs82 + DN(3,1)*clhs133 + DN(3,2)*clhs166 + clhs280;
            lhs(14,14)=DN(3,0)*clhs88 + DN(3,1)*clhs139 + DN(3,2)*clhs167 + clhs274 + clhs282*mu;
            lhs(14,15)=-clhs283;
            lhs(15,0)=clhs172*clhs182 + clhs89;
            lhs(15,1)=clhs140 + clhs172*clhs183;
            lhs(15,2)=clhs170 + clhs172*clhs184;
            lhs(15,3)=clhs185;
            lhs(15,4)=clhs172*clhs236 + clhs204;
            lhs(15,5)=clhs172*clhs237 + clhs218;
            lhs(15,6)=clhs172*clhs238 + clhs231;
            lhs(15,7)=clhs239;
            lhs(15,8)=clhs172*clhs269 + clhs251;
            lhs(15,9)=clhs172*clhs270 + clhs260;
            lhs(15,10)=clhs172*clhs271 + clhs268;
            lhs(15,11)=clhs272;
            lhs(15,12)=clhs173*clhs278;
            lhs(15,13)=clhs173*clhs281;
            lhs(15,14)=clhs173*clhs283;
            lhs(15,15)=clhs171*(clhs273 + clhs279 + clhs282);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void Stokes<SymbolicNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    SymbolicNavierStokesData<2,3> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs3 =             crhs2*mu;
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs5 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs6 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs7 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs8 =             crhs7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs4);
const double crhs9 =             crhs7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs5 + crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs4;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] + N[0]*crhs5 - N[0]*crhs6;
            rhs[2]=-DN(0,0)*crhs8 - DN(0,1)*crhs9 - N[0]*crhs2;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs4;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] + N[1]*crhs5 - N[1]*crhs6;
            rhs[5]=-DN(1,0)*crhs8 - DN(1,1)*crhs9 - N[1]*crhs2;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs4;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] + N[2]*crhs5 - N[2]*crhs6;
            rhs[8]=-DN(2,0)*crhs8 - DN(2,1)*crhs9 - N[2]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void Stokes<SymbolicNavierStokesData<2,4>>::ComputeGaussPointRHSContribution(
    SymbolicNavierStokesData<2,4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs3 =             crhs2*mu;
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs5 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs6 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs7 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs8 =             crhs7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs4);
const double crhs9 =             crhs7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs5 + crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs4;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] + N[0]*crhs5 - N[0]*crhs6;
            rhs[2]=-DN(0,0)*crhs8 - DN(0,1)*crhs9 - N[0]*crhs2;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs4;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] + N[1]*crhs5 - N[1]*crhs6;
            rhs[5]=-DN(1,0)*crhs8 - DN(1,1)*crhs9 - N[1]*crhs2;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs4;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] + N[2]*crhs5 - N[2]*crhs6;
            rhs[8]=-DN(2,0)*crhs8 - DN(2,1)*crhs9 - N[2]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void Stokes<SymbolicNavierStokesData<3,4>>::ComputeGaussPointRHSContribution(
    SymbolicNavierStokesData<3,4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs3 =             crhs2*mu;
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs5 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs6 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs7 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs8 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs9 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs10 =             crhs9*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs4);
const double crhs11 =             crhs9*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs5 + crhs6);
const double crhs12 =             crhs9*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs7 + crhs8);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs4;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs5 - N[0]*crhs6;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs3 - DN(0,2)*stress[2] + N[0]*crhs7 - N[0]*crhs8;
            rhs[3]=-DN(0,0)*crhs10 - DN(0,1)*crhs11 - DN(0,2)*crhs12 - N[0]*crhs2;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs4;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs5 - N[1]*crhs6;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs3 - DN(1,2)*stress[2] + N[1]*crhs7 - N[1]*crhs8;
            rhs[7]=-DN(1,0)*crhs10 - DN(1,1)*crhs11 - DN(1,2)*crhs12 - N[1]*crhs2;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs4;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs5 - N[2]*crhs6;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs3 - DN(2,2)*stress[2] + N[2]*crhs7 - N[2]*crhs8;
            rhs[11]=-DN(2,0)*crhs10 - DN(2,1)*crhs11 - DN(2,2)*crhs12 - N[2]*crhs2;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs3 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs4;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs3 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs5 - N[3]*crhs6;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs3 - DN(3,2)*stress[2] + N[3]*crhs7 - N[3]*crhs8;
            rhs[15]=-DN(3,0)*crhs10 - DN(3,1)*crhs11 - DN(3,2)*crhs12 - N[3]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void Stokes<SymbolicNavierStokesData<3,6>>::ComputeGaussPointRHSContribution(
    SymbolicNavierStokesData<3,6> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs3 =             crhs2*mu;
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs5 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs6 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs7 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs8 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs9 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs10 =             crhs9*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs4);
const double crhs11 =             crhs9*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs5 + crhs6);
const double crhs12 =             crhs9*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs7 + crhs8);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs4;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs5 - N[0]*crhs6;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs3 - DN(0,2)*stress[2] + N[0]*crhs7 - N[0]*crhs8;
            rhs[3]=-DN(0,0)*crhs10 - DN(0,1)*crhs11 - DN(0,2)*crhs12 - N[0]*crhs2;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs4;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs5 - N[1]*crhs6;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs3 - DN(1,2)*stress[2] + N[1]*crhs7 - N[1]*crhs8;
            rhs[7]=-DN(1,0)*crhs10 - DN(1,1)*crhs11 - DN(1,2)*crhs12 - N[1]*crhs2;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs4;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs5 - N[2]*crhs6;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs3 - DN(2,2)*stress[2] + N[2]*crhs7 - N[2]*crhs8;
            rhs[11]=-DN(2,0)*crhs10 - DN(2,1)*crhs11 - DN(2,2)*crhs12 - N[2]*crhs2;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs3 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs4;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs3 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs5 - N[3]*crhs6;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs3 - DN(3,2)*stress[2] + N[3]*crhs7 - N[3]*crhs8;
            rhs[15]=-DN(3,0)*crhs10 - DN(3,1)*crhs11 - DN(3,2)*crhs12 - N[3]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void Stokes<SymbolicNavierStokesData<3,8>>::ComputeGaussPointRHSContribution(
    SymbolicNavierStokesData<3,8> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs3 =             crhs2*mu;
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs5 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs6 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs7 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs8 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs9 =             1.0/(mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs10 =             crhs9*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs4);
const double crhs11 =             crhs9*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs5 + crhs6);
const double crhs12 =             crhs9*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs7 + crhs8);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs3 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs4;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs3 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs5 - N[0]*crhs6;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs3 - DN(0,2)*stress[2] + N[0]*crhs7 - N[0]*crhs8;
            rhs[3]=-DN(0,0)*crhs10 - DN(0,1)*crhs11 - DN(0,2)*crhs12 - N[0]*crhs2;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs3 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs4;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs3 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs5 - N[1]*crhs6;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs3 - DN(1,2)*stress[2] + N[1]*crhs7 - N[1]*crhs8;
            rhs[7]=-DN(1,0)*crhs10 - DN(1,1)*crhs11 - DN(1,2)*crhs12 - N[1]*crhs2;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs3 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs4;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs3 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs5 - N[2]*crhs6;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs3 - DN(2,2)*stress[2] + N[2]*crhs7 - N[2]*crhs8;
            rhs[11]=-DN(2,0)*crhs10 - DN(2,1)*crhs11 - DN(2,2)*crhs12 - N[2]*crhs2;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs3 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs4;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs3 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs5 - N[3]*crhs6;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs3 - DN(3,2)*stress[2] + N[3]*crhs7 - N[3]*crhs8;
            rhs[15]=-DN(3,0)*crhs10 - DN(3,1)*crhs11 - DN(3,2)*crhs12 - N[3]*crhs2;


    noalias(rRHS) += rData.Weight * rhs;
}


template <class TElementData>
void Stokes<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void Stokes<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void Stokes<TElementData>::GetValueOnIntegrationPoints(   const Variable<double> &rVariable,
                                                                        std::vector<double> &rValues,
                                                                        const ProcessInfo &rCurrentProcessInfo )
{
    if (rVariable == DIVERGENCE){

        const auto& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int num_gauss = IntegrationPoints.size();

        if (rValues.size() != num_gauss){
            rValues.resize(num_gauss);
        }

        Vector gauss_pts_jacobian_determinant = ZeroVector(num_gauss);
        GeometryData::ShapeFunctionsGradientsType DN_DX;
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, gauss_pts_jacobian_determinant, GeometryData::GI_GAUSS_2);

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; i_gauss++){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; nnode++){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ndim++){
                    DVi_DXi += gp_DN_DX(nnode, ndim) * vel[ndim];
                }
            }
            rValues[i_gauss] = DVi_DXi;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class Stokes<SymbolicNavierStokesData<2,3>>;
template class Stokes<SymbolicNavierStokesData<2,4>>;
template class Stokes<SymbolicNavierStokesData<3,4>>;
template class Stokes<SymbolicNavierStokesData<3,6>>;
template class Stokes<SymbolicNavierStokesData<3,8>>;

} // namespace Kratos