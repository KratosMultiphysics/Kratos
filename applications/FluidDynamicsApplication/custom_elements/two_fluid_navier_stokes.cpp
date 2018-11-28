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

#include "two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::~TwoFluidNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<TwoFluidNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<TwoFluidNavierStokes>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateLocalSystem(
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

        if (data.IsCut()){
            GeometryType::Pointer p_geom = this->pGetGeometry();
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg);

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side
                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                const unsigned int number_of_gauss_points = gauss_weights.size();
                array_1d<double, NumNodes> Ncenter;
                for (unsigned int i = 0; i < NumNodes; i++){
                    Ncenter[i] = 1.0/NumNodes;
                }
                const double d_gauss = inner_prod(data.Distance, Ncenter);
                for (unsigned int g = 0; g < number_of_gauss_points; g++){
                    data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                    data.CalculateDensityAtGaussPoint();
                    if (d_gauss > 0.0){
                        data.CalculateAirMaterialResponse();
                    } else {
                        this->CalculateMaterialResponse(data);
                    }

                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                }
            } else {
                MatrixType Vtot = ZeroMatrix(NumNodes * (Dim + 1), NumNodes);
                MatrixType Htot = ZeroMatrix(NumNodes, NumNodes * (Dim + 1));
                MatrixType Kee_tot = ZeroMatrix(NumNodes, NumNodes);
                VectorType rhs_ee_tot = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); g_pos++){
                    data.UpdateGeometryValues(
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);

                    data.CalculateDensityAtGaussPoint();
                    data.CalculateAirMaterialResponse();
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); g_neg++){
                    data.UpdateGeometryValues(
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);

                    data.CalculateDensityAtGaussPoint();
                    this->CalculateMaterialResponse(data);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                CondenseEnrichment(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
            }
        } else {
            //Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            ShapeFunctionDerivativesArrayType shape_derivatives;
            this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
            const unsigned int number_of_gauss_points = gauss_weights.size();
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++){
                data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                data.CalculateDensityAtGaussPoint();
                this->CalculateMaterialResponse(data);
                this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
    } else{
        KRATOS_ERROR << "TwoFluidNavierStokes is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
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
std::string TwoFluidNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr){
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}


template <class TElementData>
void TwoFluidNavierStokes<TElementData>::Calculate( const Variable<Vector >& rVariable,
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

            dataLocal.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);

            dataLocal.CalculateDensityAtGaussPoint();
            dataLocal.CalculateEffectiveViscosityAtGaussPoint();
            dataLocal.CalculateAirMaterialResponse();

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
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateMaterialResponse(TElementData &rData) const
{
    if (rData.IsAir()){
        rData.CalculateAirMaterialResponse();
    } else {
        FluidElement<TElementData>::CalculateMaterialResponse(rData);
    }
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesData<2, 3> &rData,
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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             C(0,2)*DN(0,0);
const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
const double clhs3 =             pow(DN(0,0), 2);
const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 =             rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             bdf0*rho;
const double clhs9 =             N[0]*rho;
const double clhs10 =             DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs11 =             N[0]*bdf0 + clhs10;
const double clhs12 =             pow(rho, 2);
const double clhs13 =             1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs14 =             1.0*clhs10*clhs12*clhs13;
const double clhs15 =             pow(N[0], 2)*clhs8 + clhs10*clhs9 + clhs11*clhs14;
const double clhs16 =             C(0,1)*DN(0,1) + clhs1;
const double clhs17 =             C(1,2)*DN(0,1);
const double clhs18 =             C(2,2)*DN(0,0) + clhs17;
const double clhs19 =             DN(0,0)*clhs7;
const double clhs20 =             DN(0,1)*clhs19;
const double clhs21 =             1.0*clhs13*rho;
const double clhs22 =             clhs10*clhs21;
const double clhs23 =             -N[0] + clhs22;
const double clhs24 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs25 =             C(0,2)*DN(1,0);
const double clhs26 =             C(2,2)*DN(1,1) + clhs25;
const double clhs27 =             DN(0,0)*DN(1,0);
const double clhs28 =             clhs27*clhs7;
const double clhs29 =             N[0]*bdf0*rho;
const double clhs30 =             N[1]*clhs29;
const double clhs31 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs32 =             N[1]*bdf0;
const double clhs33 =             clhs31 + clhs32;
const double clhs34 =             clhs14*clhs33 + clhs30 + clhs31*clhs9;
const double clhs35 =             C(0,1)*DN(1,1) + clhs25;
const double clhs36 =             C(1,2)*DN(1,1);
const double clhs37 =             C(2,2)*DN(1,0) + clhs36;
const double clhs38 =             DN(1,1)*clhs19;
const double clhs39 =             DN(0,0)*N[1];
const double clhs40 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs41 =             C(0,2)*DN(2,0);
const double clhs42 =             C(2,2)*DN(2,1) + clhs41;
const double clhs43 =             DN(0,0)*DN(2,0);
const double clhs44 =             clhs43*clhs7;
const double clhs45 =             N[2]*clhs29;
const double clhs46 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs47 =             N[2]*bdf0 + clhs46;
const double clhs48 =             clhs14*clhs47 + clhs45 + clhs46*clhs9;
const double clhs49 =             C(0,1)*DN(2,1) + clhs41;
const double clhs50 =             C(1,2)*DN(2,1);
const double clhs51 =             C(2,2)*DN(2,0) + clhs50;
const double clhs52 =             DN(2,1)*clhs19;
const double clhs53 =             DN(0,0)*N[2];
const double clhs54 =             C(0,1)*DN(0,0) + clhs17;
const double clhs55 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs56 =             pow(DN(0,1), 2);
const double clhs57 =             C(0,1)*DN(1,0) + clhs36;
const double clhs58 =             DN(0,1)*clhs7;
const double clhs59 =             DN(1,0)*clhs58;
const double clhs60 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs61 =             DN(0,1)*DN(1,1);
const double clhs62 =             clhs61*clhs7;
const double clhs63 =             DN(0,1)*N[1];
const double clhs64 =             C(0,1)*DN(2,0) + clhs50;
const double clhs65 =             DN(2,0)*clhs58;
const double clhs66 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs67 =             DN(0,1)*DN(2,1);
const double clhs68 =             clhs67*clhs7;
const double clhs69 =             DN(0,1)*N[2];
const double clhs70 =             clhs11*clhs21;
const double clhs71 =             N[0] + clhs70;
const double clhs72 =             1.0*clhs13;
const double clhs73 =             DN(1,0)*N[0];
const double clhs74 =             clhs21*clhs33;
const double clhs75 =             DN(1,1)*N[0];
const double clhs76 =             clhs72*(clhs27 + clhs61);
const double clhs77 =             DN(2,0)*N[0];
const double clhs78 =             clhs21*clhs47;
const double clhs79 =             DN(2,1)*N[0];
const double clhs80 =             clhs72*(clhs43 + clhs67);
const double clhs81 =             N[1]*rho;
const double clhs82 =             1.0*clhs12*clhs13*clhs31;
const double clhs83 =             clhs10*clhs81 + clhs11*clhs82 + clhs30;
const double clhs84 =             clhs21*clhs31;
const double clhs85 =             pow(DN(1,0), 2);
const double clhs86 =             pow(N[1], 2)*clhs8 + clhs31*clhs81 + clhs33*clhs82;
const double clhs87 =             DN(1,0)*clhs7;
const double clhs88 =             DN(1,1)*clhs87;
const double clhs89 =             -N[1] + clhs84;
const double clhs90 =             DN(1,0)*DN(2,0);
const double clhs91 =             clhs7*clhs90;
const double clhs92 =             N[2]*rho;
const double clhs93 =             clhs32*clhs92;
const double clhs94 =             clhs46*clhs81 + clhs47*clhs82 + clhs93;
const double clhs95 =             DN(2,1)*clhs87;
const double clhs96 =             DN(1,0)*N[2];
const double clhs97 =             pow(DN(1,1), 2);
const double clhs98 =             DN(2,0)*clhs7;
const double clhs99 =             DN(1,1)*clhs98;
const double clhs100 =             DN(1,1)*DN(2,1);
const double clhs101 =             clhs100*clhs7;
const double clhs102 =             DN(1,1)*N[2];
const double clhs103 =             N[1] + clhs74;
const double clhs104 =             DN(2,0)*N[1];
const double clhs105 =             DN(2,1)*N[1];
const double clhs106 =             clhs72*(clhs100 + clhs90);
const double clhs107 =             1.0*clhs12*clhs13*clhs46;
const double clhs108 =             clhs10*clhs92 + clhs107*clhs11 + clhs45;
const double clhs109 =             clhs21*clhs46;
const double clhs110 =             clhs107*clhs33 + clhs31*clhs92 + clhs93;
const double clhs111 =             pow(DN(2,0), 2);
const double clhs112 =             pow(N[2], 2)*clhs8 + clhs107*clhs47 + clhs46*clhs92;
const double clhs113 =             DN(2,1)*clhs98;
const double clhs114 =             -N[2] + clhs109;
const double clhs115 =             pow(DN(2,1), 2);
const double clhs116 =             N[2] + clhs78;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs15 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs16 + DN(0,1)*clhs18 + clhs20;
            lhs(0,2)=DN(0,0)*clhs23;
            lhs(0,3)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs28 + clhs34;
            lhs(0,4)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + clhs38;
            lhs(0,5)=DN(1,0)*clhs22 - clhs39;
            lhs(0,6)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + clhs44 + clhs48;
            lhs(0,7)=DN(0,0)*clhs49 + DN(0,1)*clhs51 + clhs52;
            lhs(0,8)=DN(2,0)*clhs22 - clhs53;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs54 + clhs20;
            lhs(1,1)=DN(0,0)*clhs18 + DN(0,1)*clhs55 + clhs15 + clhs56*clhs7;
            lhs(1,2)=DN(0,1)*clhs23;
            lhs(1,3)=DN(0,0)*clhs26 + DN(0,1)*clhs57 + clhs59;
            lhs(1,4)=DN(0,0)*clhs37 + DN(0,1)*clhs60 + clhs34 + clhs62;
            lhs(1,5)=DN(1,1)*clhs22 - clhs63;
            lhs(1,6)=DN(0,0)*clhs42 + DN(0,1)*clhs64 + clhs65;
            lhs(1,7)=DN(0,0)*clhs51 + DN(0,1)*clhs66 + clhs48 + clhs68;
            lhs(1,8)=DN(2,1)*clhs22 - clhs69;
            lhs(2,0)=DN(0,0)*clhs71;
            lhs(2,1)=DN(0,1)*clhs71;
            lhs(2,2)=clhs72*(clhs3 + clhs56);
            lhs(2,3)=DN(0,0)*clhs74 + clhs73;
            lhs(2,4)=DN(0,1)*clhs74 + clhs75;
            lhs(2,5)=clhs76;
            lhs(2,6)=DN(0,0)*clhs78 + clhs77;
            lhs(2,7)=DN(0,1)*clhs78 + clhs79;
            lhs(2,8)=clhs80;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs28 + clhs83;
            lhs(3,1)=DN(1,0)*clhs16 + DN(1,1)*clhs18 + clhs59;
            lhs(3,2)=DN(0,0)*clhs84 - clhs73;
            lhs(3,3)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs7*clhs85 + clhs86;
            lhs(3,4)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + clhs88;
            lhs(3,5)=DN(1,0)*clhs89;
            lhs(3,6)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + clhs91 + clhs94;
            lhs(3,7)=DN(1,0)*clhs49 + DN(1,1)*clhs51 + clhs95;
            lhs(3,8)=DN(2,0)*clhs84 - clhs96;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs54 + clhs38;
            lhs(4,1)=DN(1,0)*clhs18 + DN(1,1)*clhs55 + clhs62 + clhs83;
            lhs(4,2)=DN(0,1)*clhs84 - clhs75;
            lhs(4,3)=DN(1,0)*clhs26 + DN(1,1)*clhs57 + clhs88;
            lhs(4,4)=DN(1,0)*clhs37 + DN(1,1)*clhs60 + clhs7*clhs97 + clhs86;
            lhs(4,5)=DN(1,1)*clhs89;
            lhs(4,6)=DN(1,0)*clhs42 + DN(1,1)*clhs64 + clhs99;
            lhs(4,7)=DN(1,0)*clhs51 + DN(1,1)*clhs66 + clhs101 + clhs94;
            lhs(4,8)=DN(2,1)*clhs84 - clhs102;
            lhs(5,0)=DN(1,0)*clhs70 + clhs39;
            lhs(5,1)=DN(1,1)*clhs70 + clhs63;
            lhs(5,2)=clhs76;
            lhs(5,3)=DN(1,0)*clhs103;
            lhs(5,4)=DN(1,1)*clhs103;
            lhs(5,5)=clhs72*(clhs85 + clhs97);
            lhs(5,6)=DN(1,0)*clhs78 + clhs104;
            lhs(5,7)=DN(1,1)*clhs78 + clhs105;
            lhs(5,8)=clhs106;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs108 + clhs44;
            lhs(6,1)=DN(2,0)*clhs16 + DN(2,1)*clhs18 + clhs65;
            lhs(6,2)=DN(0,0)*clhs109 - clhs77;
            lhs(6,3)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs110 + clhs91;
            lhs(6,4)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + clhs99;
            lhs(6,5)=DN(1,0)*clhs109 - clhs104;
            lhs(6,6)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + clhs111*clhs7 + clhs112;
            lhs(6,7)=DN(2,0)*clhs49 + DN(2,1)*clhs51 + clhs113;
            lhs(6,8)=DN(2,0)*clhs114;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs54 + clhs52;
            lhs(7,1)=DN(2,0)*clhs18 + DN(2,1)*clhs55 + clhs108 + clhs68;
            lhs(7,2)=DN(0,1)*clhs109 - clhs79;
            lhs(7,3)=DN(2,0)*clhs26 + DN(2,1)*clhs57 + clhs95;
            lhs(7,4)=DN(2,0)*clhs37 + DN(2,1)*clhs60 + clhs101 + clhs110;
            lhs(7,5)=DN(1,1)*clhs109 - clhs105;
            lhs(7,6)=DN(2,0)*clhs42 + DN(2,1)*clhs64 + clhs113;
            lhs(7,7)=DN(2,0)*clhs51 + DN(2,1)*clhs66 + clhs112 + clhs115*clhs7;
            lhs(7,8)=DN(2,1)*clhs114;
            lhs(8,0)=DN(2,0)*clhs70 + clhs53;
            lhs(8,1)=DN(2,1)*clhs70 + clhs69;
            lhs(8,2)=clhs80;
            lhs(8,3)=DN(2,0)*clhs74 + clhs96;
            lhs(8,4)=DN(2,1)*clhs74 + clhs102;
            lhs(8,5)=clhs106;
            lhs(8,6)=DN(2,0)*clhs116;
            lhs(8,7)=DN(2,1)*clhs116;
            lhs(8,8)=clhs72*(clhs111 + clhs115);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesData<3, 4> &rData,
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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             C(0,3)*DN(0,0);
const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 =             C(0,5)*DN(0,0);
const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs9*h/stab_c1 + mu;
const double clhs11 =             bdf0*rho;
const double clhs12 =             N[0]*rho;
const double clhs13 =             DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs14 =             N[0]*bdf0 + clhs13;
const double clhs15 =             pow(rho, 2);
const double clhs16 =             1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs17 =             1.0*clhs13*clhs15*clhs16;
const double clhs18 =             pow(N[0], 2)*clhs11 + clhs12*clhs13 + clhs14*clhs17;
const double clhs19 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs20 =             C(1,3)*DN(0,1);
const double clhs21 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs20;
const double clhs22 =             C(3,5)*DN(0,0);
const double clhs23 =             C(4,5)*DN(0,2);
const double clhs24 =             C(1,5)*DN(0,1) + clhs22 + clhs23;
const double clhs25 =             DN(0,0)*clhs10;
const double clhs26 =             DN(0,1)*clhs25;
const double clhs27 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs28 =             C(3,4)*DN(0,1);
const double clhs29 =             C(2,3)*DN(0,2) + clhs22 + clhs28;
const double clhs30 =             C(2,5)*DN(0,2);
const double clhs31 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs30;
const double clhs32 =             DN(0,2)*clhs25;
const double clhs33 =             1.0*clhs16*rho;
const double clhs34 =             clhs13*clhs33;
const double clhs35 =             -N[0] + clhs34;
const double clhs36 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs37 =             C(0,3)*DN(1,0);
const double clhs38 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs37;
const double clhs39 =             C(0,5)*DN(1,0);
const double clhs40 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs39;
const double clhs41 =             DN(0,0)*DN(1,0);
const double clhs42 =             clhs10*clhs41;
const double clhs43 =             N[0]*bdf0*rho;
const double clhs44 =             N[1]*clhs43;
const double clhs45 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs46 =             N[1]*bdf0 + clhs45;
const double clhs47 =             clhs12*clhs45 + clhs17*clhs46 + clhs44;
const double clhs48 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs37;
const double clhs49 =             C(1,3)*DN(1,1);
const double clhs50 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs49;
const double clhs51 =             C(3,5)*DN(1,0);
const double clhs52 =             C(4,5)*DN(1,2);
const double clhs53 =             C(1,5)*DN(1,1) + clhs51 + clhs52;
const double clhs54 =             DN(1,1)*clhs25;
const double clhs55 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs39;
const double clhs56 =             C(3,4)*DN(1,1);
const double clhs57 =             C(2,3)*DN(1,2) + clhs51 + clhs56;
const double clhs58 =             C(2,5)*DN(1,2);
const double clhs59 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs58;
const double clhs60 =             DN(1,2)*clhs25;
const double clhs61 =             DN(0,0)*N[1];
const double clhs62 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs63 =             C(0,3)*DN(2,0);
const double clhs64 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs63;
const double clhs65 =             C(0,5)*DN(2,0);
const double clhs66 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs65;
const double clhs67 =             DN(0,0)*DN(2,0);
const double clhs68 =             clhs10*clhs67;
const double clhs69 =             N[2]*clhs43;
const double clhs70 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs71 =             N[2]*bdf0;
const double clhs72 =             clhs70 + clhs71;
const double clhs73 =             clhs12*clhs70 + clhs17*clhs72 + clhs69;
const double clhs74 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs63;
const double clhs75 =             C(1,3)*DN(2,1);
const double clhs76 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs75;
const double clhs77 =             C(3,5)*DN(2,0);
const double clhs78 =             C(4,5)*DN(2,2);
const double clhs79 =             C(1,5)*DN(2,1) + clhs77 + clhs78;
const double clhs80 =             DN(2,1)*clhs25;
const double clhs81 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs65;
const double clhs82 =             C(3,4)*DN(2,1);
const double clhs83 =             C(2,3)*DN(2,2) + clhs77 + clhs82;
const double clhs84 =             C(2,5)*DN(2,2);
const double clhs85 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs84;
const double clhs86 =             DN(2,2)*clhs25;
const double clhs87 =             DN(0,0)*N[2];
const double clhs88 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs89 =             C(0,3)*DN(3,0);
const double clhs90 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs89;
const double clhs91 =             C(0,5)*DN(3,0);
const double clhs92 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs91;
const double clhs93 =             DN(0,0)*DN(3,0);
const double clhs94 =             clhs10*clhs93;
const double clhs95 =             N[3]*clhs43;
const double clhs96 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs97 =             N[3]*bdf0 + clhs96;
const double clhs98 =             clhs12*clhs96 + clhs17*clhs97 + clhs95;
const double clhs99 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs89;
const double clhs100 =             C(1,3)*DN(3,1);
const double clhs101 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs100;
const double clhs102 =             C(3,5)*DN(3,0);
const double clhs103 =             C(4,5)*DN(3,2);
const double clhs104 =             C(1,5)*DN(3,1) + clhs102 + clhs103;
const double clhs105 =             DN(3,1)*clhs25;
const double clhs106 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs91;
const double clhs107 =             C(3,4)*DN(3,1);
const double clhs108 =             C(2,3)*DN(3,2) + clhs102 + clhs107;
const double clhs109 =             C(2,5)*DN(3,2);
const double clhs110 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs109;
const double clhs111 =             DN(3,2)*clhs25;
const double clhs112 =             DN(0,0)*N[3];
const double clhs113 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs20;
const double clhs114 =             C(0,4)*DN(0,0) + clhs23 + clhs28;
const double clhs115 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs116 =             C(1,4)*DN(0,1);
const double clhs117 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs116;
const double clhs118 =             pow(DN(0,1), 2);
const double clhs119 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs116;
const double clhs120 =             C(2,4)*DN(0,2);
const double clhs121 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs120;
const double clhs122 =             DN(0,1)*clhs10;
const double clhs123 =             DN(0,2)*clhs122;
const double clhs124 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs49;
const double clhs125 =             C(0,4)*DN(1,0) + clhs52 + clhs56;
const double clhs126 =             DN(1,0)*clhs122;
const double clhs127 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs128 =             C(1,4)*DN(1,1);
const double clhs129 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs128;
const double clhs130 =             DN(0,1)*DN(1,1);
const double clhs131 =             clhs10*clhs130;
const double clhs132 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs128;
const double clhs133 =             C(2,4)*DN(1,2);
const double clhs134 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs133;
const double clhs135 =             DN(1,2)*clhs122;
const double clhs136 =             DN(0,1)*N[1];
const double clhs137 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs75;
const double clhs138 =             C(0,4)*DN(2,0) + clhs78 + clhs82;
const double clhs139 =             DN(2,0)*clhs122;
const double clhs140 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs141 =             C(1,4)*DN(2,1);
const double clhs142 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs141;
const double clhs143 =             DN(0,1)*DN(2,1);
const double clhs144 =             clhs10*clhs143;
const double clhs145 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs141;
const double clhs146 =             C(2,4)*DN(2,2);
const double clhs147 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs146;
const double clhs148 =             DN(2,2)*clhs122;
const double clhs149 =             DN(0,1)*N[2];
const double clhs150 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs100;
const double clhs151 =             C(0,4)*DN(3,0) + clhs103 + clhs107;
const double clhs152 =             DN(3,0)*clhs122;
const double clhs153 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs154 =             C(1,4)*DN(3,1);
const double clhs155 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs154;
const double clhs156 =             DN(0,1)*DN(3,1);
const double clhs157 =             clhs10*clhs156;
const double clhs158 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs154;
const double clhs159 =             C(2,4)*DN(3,2);
const double clhs160 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs159;
const double clhs161 =             DN(3,2)*clhs122;
const double clhs162 =             DN(0,1)*N[3];
const double clhs163 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs30;
const double clhs164 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs120;
const double clhs165 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs166 =             pow(DN(0,2), 2);
const double clhs167 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs58;
const double clhs168 =             DN(0,2)*clhs10;
const double clhs169 =             DN(1,0)*clhs168;
const double clhs170 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs133;
const double clhs171 =             DN(1,1)*clhs168;
const double clhs172 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs173 =             DN(0,2)*DN(1,2);
const double clhs174 =             clhs10*clhs173;
const double clhs175 =             DN(0,2)*N[1];
const double clhs176 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs84;
const double clhs177 =             DN(2,0)*clhs168;
const double clhs178 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs146;
const double clhs179 =             DN(2,1)*clhs168;
const double clhs180 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs181 =             DN(0,2)*DN(2,2);
const double clhs182 =             clhs10*clhs181;
const double clhs183 =             DN(0,2)*N[2];
const double clhs184 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs109;
const double clhs185 =             DN(3,0)*clhs168;
const double clhs186 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs159;
const double clhs187 =             DN(3,1)*clhs168;
const double clhs188 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs189 =             DN(0,2)*DN(3,2);
const double clhs190 =             clhs10*clhs189;
const double clhs191 =             DN(0,2)*N[3];
const double clhs192 =             clhs14*clhs33;
const double clhs193 =             N[0] + clhs192;
const double clhs194 =             1.0*clhs16;
const double clhs195 =             DN(1,0)*N[0];
const double clhs196 =             clhs33*clhs46;
const double clhs197 =             DN(1,1)*N[0];
const double clhs198 =             DN(1,2)*N[0];
const double clhs199 =             clhs194*(clhs130 + clhs173 + clhs41);
const double clhs200 =             DN(2,0)*N[0];
const double clhs201 =             clhs33*clhs72;
const double clhs202 =             DN(2,1)*N[0];
const double clhs203 =             DN(2,2)*N[0];
const double clhs204 =             clhs194*(clhs143 + clhs181 + clhs67);
const double clhs205 =             DN(3,0)*N[0];
const double clhs206 =             clhs33*clhs97;
const double clhs207 =             DN(3,1)*N[0];
const double clhs208 =             DN(3,2)*N[0];
const double clhs209 =             clhs194*(clhs156 + clhs189 + clhs93);
const double clhs210 =             N[1]*rho;
const double clhs211 =             1.0*clhs15*clhs16*clhs45;
const double clhs212 =             clhs13*clhs210 + clhs14*clhs211 + clhs44;
const double clhs213 =             clhs33*clhs45;
const double clhs214 =             pow(DN(1,0), 2);
const double clhs215 =             pow(N[1], 2)*clhs11 + clhs210*clhs45 + clhs211*clhs46;
const double clhs216 =             DN(1,0)*clhs10;
const double clhs217 =             DN(1,1)*clhs216;
const double clhs218 =             DN(1,2)*clhs216;
const double clhs219 =             -N[1] + clhs213;
const double clhs220 =             DN(1,0)*DN(2,0);
const double clhs221 =             clhs10*clhs220;
const double clhs222 =             N[1]*bdf0*rho;
const double clhs223 =             N[2]*clhs222;
const double clhs224 =             clhs210*clhs70 + clhs211*clhs72 + clhs223;
const double clhs225 =             DN(2,1)*clhs216;
const double clhs226 =             DN(2,2)*clhs216;
const double clhs227 =             DN(1,0)*N[2];
const double clhs228 =             DN(1,0)*DN(3,0);
const double clhs229 =             clhs10*clhs228;
const double clhs230 =             N[3]*clhs222;
const double clhs231 =             clhs210*clhs96 + clhs211*clhs97 + clhs230;
const double clhs232 =             DN(3,1)*clhs216;
const double clhs233 =             DN(3,2)*clhs216;
const double clhs234 =             DN(1,0)*N[3];
const double clhs235 =             pow(DN(1,1), 2);
const double clhs236 =             DN(1,1)*clhs10;
const double clhs237 =             DN(1,2)*clhs236;
const double clhs238 =             DN(2,0)*clhs236;
const double clhs239 =             DN(1,1)*DN(2,1);
const double clhs240 =             clhs10*clhs239;
const double clhs241 =             DN(2,2)*clhs236;
const double clhs242 =             DN(1,1)*N[2];
const double clhs243 =             DN(3,0)*clhs236;
const double clhs244 =             DN(1,1)*DN(3,1);
const double clhs245 =             clhs10*clhs244;
const double clhs246 =             DN(3,2)*clhs236;
const double clhs247 =             DN(1,1)*N[3];
const double clhs248 =             pow(DN(1,2), 2);
const double clhs249 =             DN(1,2)*clhs10;
const double clhs250 =             DN(2,0)*clhs249;
const double clhs251 =             DN(2,1)*clhs249;
const double clhs252 =             DN(1,2)*DN(2,2);
const double clhs253 =             clhs10*clhs252;
const double clhs254 =             DN(1,2)*N[2];
const double clhs255 =             DN(3,0)*clhs249;
const double clhs256 =             DN(3,1)*clhs249;
const double clhs257 =             DN(1,2)*DN(3,2);
const double clhs258 =             clhs10*clhs257;
const double clhs259 =             DN(1,2)*N[3];
const double clhs260 =             N[1] + clhs196;
const double clhs261 =             DN(2,0)*N[1];
const double clhs262 =             DN(2,1)*N[1];
const double clhs263 =             DN(2,2)*N[1];
const double clhs264 =             clhs194*(clhs220 + clhs239 + clhs252);
const double clhs265 =             DN(3,0)*N[1];
const double clhs266 =             DN(3,1)*N[1];
const double clhs267 =             DN(3,2)*N[1];
const double clhs268 =             clhs194*(clhs228 + clhs244 + clhs257);
const double clhs269 =             N[2]*rho;
const double clhs270 =             1.0*clhs15*clhs16*clhs70;
const double clhs271 =             clhs13*clhs269 + clhs14*clhs270 + clhs69;
const double clhs272 =             clhs33*clhs70;
const double clhs273 =             clhs223 + clhs269*clhs45 + clhs270*clhs46;
const double clhs274 =             pow(DN(2,0), 2);
const double clhs275 =             pow(N[2], 2)*clhs11 + clhs269*clhs70 + clhs270*clhs72;
const double clhs276 =             DN(2,0)*clhs10;
const double clhs277 =             DN(2,1)*clhs276;
const double clhs278 =             DN(2,2)*clhs276;
const double clhs279 =             -N[2] + clhs272;
const double clhs280 =             DN(2,0)*DN(3,0);
const double clhs281 =             clhs10*clhs280;
const double clhs282 =             N[3]*rho;
const double clhs283 =             clhs282*clhs71;
const double clhs284 =             clhs269*clhs96 + clhs270*clhs97 + clhs283;
const double clhs285 =             DN(3,1)*clhs276;
const double clhs286 =             DN(3,2)*clhs276;
const double clhs287 =             DN(2,0)*N[3];
const double clhs288 =             pow(DN(2,1), 2);
const double clhs289 =             DN(2,1)*clhs10;
const double clhs290 =             DN(2,2)*clhs289;
const double clhs291 =             DN(3,0)*clhs289;
const double clhs292 =             DN(2,1)*DN(3,1);
const double clhs293 =             clhs10*clhs292;
const double clhs294 =             DN(3,2)*clhs289;
const double clhs295 =             DN(2,1)*N[3];
const double clhs296 =             pow(DN(2,2), 2);
const double clhs297 =             DN(2,2)*clhs10;
const double clhs298 =             DN(3,0)*clhs297;
const double clhs299 =             DN(3,1)*clhs297;
const double clhs300 =             DN(2,2)*DN(3,2);
const double clhs301 =             clhs10*clhs300;
const double clhs302 =             DN(2,2)*N[3];
const double clhs303 =             N[2] + clhs201;
const double clhs304 =             DN(3,0)*N[2];
const double clhs305 =             DN(3,1)*N[2];
const double clhs306 =             DN(3,2)*N[2];
const double clhs307 =             clhs194*(clhs280 + clhs292 + clhs300);
const double clhs308 =             1.0*clhs15*clhs16*clhs96;
const double clhs309 =             clhs13*clhs282 + clhs14*clhs308 + clhs95;
const double clhs310 =             clhs33*clhs96;
const double clhs311 =             clhs230 + clhs282*clhs45 + clhs308*clhs46;
const double clhs312 =             clhs282*clhs70 + clhs283 + clhs308*clhs72;
const double clhs313 =             pow(DN(3,0), 2);
const double clhs314 =             pow(N[3], 2)*clhs11 + clhs282*clhs96 + clhs308*clhs97;
const double clhs315 =             DN(3,0)*clhs10;
const double clhs316 =             DN(3,1)*clhs315;
const double clhs317 =             DN(3,2)*clhs315;
const double clhs318 =             -N[3] + clhs310;
const double clhs319 =             pow(DN(3,1), 2);
const double clhs320 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs321 =             pow(DN(3,2), 2);
const double clhs322 =             N[3] + clhs206;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs18;
            lhs(0,1)=DN(0,0)*clhs19 + DN(0,1)*clhs21 + DN(0,2)*clhs24 + clhs26;
            lhs(0,2)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + DN(0,2)*clhs31 + clhs32;
            lhs(0,3)=DN(0,0)*clhs35;
            lhs(0,4)=DN(0,0)*clhs36 + DN(0,1)*clhs38 + DN(0,2)*clhs40 + clhs42 + clhs47;
            lhs(0,5)=DN(0,0)*clhs48 + DN(0,1)*clhs50 + DN(0,2)*clhs53 + clhs54;
            lhs(0,6)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs59 + clhs60;
            lhs(0,7)=DN(1,0)*clhs34 - clhs61;
            lhs(0,8)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs68 + clhs73;
            lhs(0,9)=DN(0,0)*clhs74 + DN(0,1)*clhs76 + DN(0,2)*clhs79 + clhs80;
            lhs(0,10)=DN(0,0)*clhs81 + DN(0,1)*clhs83 + DN(0,2)*clhs85 + clhs86;
            lhs(0,11)=DN(2,0)*clhs34 - clhs87;
            lhs(0,12)=DN(0,0)*clhs88 + DN(0,1)*clhs90 + DN(0,2)*clhs92 + clhs94 + clhs98;
            lhs(0,13)=DN(0,0)*clhs99 + DN(0,1)*clhs101 + DN(0,2)*clhs104 + clhs105;
            lhs(0,14)=DN(0,0)*clhs106 + DN(0,1)*clhs108 + DN(0,2)*clhs110 + clhs111;
            lhs(0,15)=DN(3,0)*clhs34 - clhs112;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs113 + DN(0,2)*clhs114 + clhs26;
            lhs(1,1)=DN(0,0)*clhs21 + DN(0,1)*clhs115 + DN(0,2)*clhs117 + clhs10*clhs118 + clhs18;
            lhs(1,2)=DN(0,0)*clhs29 + DN(0,1)*clhs119 + DN(0,2)*clhs121 + clhs123;
            lhs(1,3)=DN(0,1)*clhs35;
            lhs(1,4)=DN(0,0)*clhs38 + DN(0,1)*clhs124 + DN(0,2)*clhs125 + clhs126;
            lhs(1,5)=DN(0,0)*clhs50 + DN(0,1)*clhs127 + DN(0,2)*clhs129 + clhs131 + clhs47;
            lhs(1,6)=DN(0,0)*clhs57 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs135;
            lhs(1,7)=DN(1,1)*clhs34 - clhs136;
            lhs(1,8)=DN(0,0)*clhs64 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
            lhs(1,9)=DN(0,0)*clhs76 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs144 + clhs73;
            lhs(1,10)=DN(0,0)*clhs83 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148;
            lhs(1,11)=DN(2,1)*clhs34 - clhs149;
            lhs(1,12)=DN(0,0)*clhs90 + DN(0,1)*clhs150 + DN(0,2)*clhs151 + clhs152;
            lhs(1,13)=DN(0,0)*clhs101 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs157 + clhs98;
            lhs(1,14)=DN(0,0)*clhs108 + DN(0,1)*clhs158 + DN(0,2)*clhs160 + clhs161;
            lhs(1,15)=DN(3,1)*clhs34 - clhs162;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs114 + DN(0,2)*clhs163 + clhs32;
            lhs(2,1)=DN(0,0)*clhs24 + DN(0,1)*clhs117 + DN(0,2)*clhs164 + clhs123;
            lhs(2,2)=DN(0,0)*clhs31 + DN(0,1)*clhs121 + DN(0,2)*clhs165 + clhs10*clhs166 + clhs18;
            lhs(2,3)=DN(0,2)*clhs35;
            lhs(2,4)=DN(0,0)*clhs40 + DN(0,1)*clhs125 + DN(0,2)*clhs167 + clhs169;
            lhs(2,5)=DN(0,0)*clhs53 + DN(0,1)*clhs129 + DN(0,2)*clhs170 + clhs171;
            lhs(2,6)=DN(0,0)*clhs59 + DN(0,1)*clhs134 + DN(0,2)*clhs172 + clhs174 + clhs47;
            lhs(2,7)=DN(1,2)*clhs34 - clhs175;
            lhs(2,8)=DN(0,0)*clhs66 + DN(0,1)*clhs138 + DN(0,2)*clhs176 + clhs177;
            lhs(2,9)=DN(0,0)*clhs79 + DN(0,1)*clhs142 + DN(0,2)*clhs178 + clhs179;
            lhs(2,10)=DN(0,0)*clhs85 + DN(0,1)*clhs147 + DN(0,2)*clhs180 + clhs182 + clhs73;
            lhs(2,11)=DN(2,2)*clhs34 - clhs183;
            lhs(2,12)=DN(0,0)*clhs92 + DN(0,1)*clhs151 + DN(0,2)*clhs184 + clhs185;
            lhs(2,13)=DN(0,0)*clhs104 + DN(0,1)*clhs155 + DN(0,2)*clhs186 + clhs187;
            lhs(2,14)=DN(0,0)*clhs110 + DN(0,1)*clhs160 + DN(0,2)*clhs188 + clhs190 + clhs98;
            lhs(2,15)=DN(3,2)*clhs34 - clhs191;
            lhs(3,0)=DN(0,0)*clhs193;
            lhs(3,1)=DN(0,1)*clhs193;
            lhs(3,2)=DN(0,2)*clhs193;
            lhs(3,3)=clhs194*(clhs118 + clhs166 + clhs5);
            lhs(3,4)=DN(0,0)*clhs196 + clhs195;
            lhs(3,5)=DN(0,1)*clhs196 + clhs197;
            lhs(3,6)=DN(0,2)*clhs196 + clhs198;
            lhs(3,7)=clhs199;
            lhs(3,8)=DN(0,0)*clhs201 + clhs200;
            lhs(3,9)=DN(0,1)*clhs201 + clhs202;
            lhs(3,10)=DN(0,2)*clhs201 + clhs203;
            lhs(3,11)=clhs204;
            lhs(3,12)=DN(0,0)*clhs206 + clhs205;
            lhs(3,13)=DN(0,1)*clhs206 + clhs207;
            lhs(3,14)=DN(0,2)*clhs206 + clhs208;
            lhs(3,15)=clhs209;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs212 + clhs42;
            lhs(4,1)=DN(1,0)*clhs19 + DN(1,1)*clhs21 + DN(1,2)*clhs24 + clhs126;
            lhs(4,2)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + DN(1,2)*clhs31 + clhs169;
            lhs(4,3)=DN(0,0)*clhs213 - clhs195;
            lhs(4,4)=DN(1,0)*clhs36 + DN(1,1)*clhs38 + DN(1,2)*clhs40 + clhs10*clhs214 + clhs215;
            lhs(4,5)=DN(1,0)*clhs48 + DN(1,1)*clhs50 + DN(1,2)*clhs53 + clhs217;
            lhs(4,6)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs59 + clhs218;
            lhs(4,7)=DN(1,0)*clhs219;
            lhs(4,8)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs221 + clhs224;
            lhs(4,9)=DN(1,0)*clhs74 + DN(1,1)*clhs76 + DN(1,2)*clhs79 + clhs225;
            lhs(4,10)=DN(1,0)*clhs81 + DN(1,1)*clhs83 + DN(1,2)*clhs85 + clhs226;
            lhs(4,11)=DN(2,0)*clhs213 - clhs227;
            lhs(4,12)=DN(1,0)*clhs88 + DN(1,1)*clhs90 + DN(1,2)*clhs92 + clhs229 + clhs231;
            lhs(4,13)=DN(1,0)*clhs99 + DN(1,1)*clhs101 + DN(1,2)*clhs104 + clhs232;
            lhs(4,14)=DN(1,0)*clhs106 + DN(1,1)*clhs108 + DN(1,2)*clhs110 + clhs233;
            lhs(4,15)=DN(3,0)*clhs213 - clhs234;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs113 + DN(1,2)*clhs114 + clhs54;
            lhs(5,1)=DN(1,0)*clhs21 + DN(1,1)*clhs115 + DN(1,2)*clhs117 + clhs131 + clhs212;
            lhs(5,2)=DN(1,0)*clhs29 + DN(1,1)*clhs119 + DN(1,2)*clhs121 + clhs171;
            lhs(5,3)=DN(0,1)*clhs213 - clhs197;
            lhs(5,4)=DN(1,0)*clhs38 + DN(1,1)*clhs124 + DN(1,2)*clhs125 + clhs217;
            lhs(5,5)=DN(1,0)*clhs50 + DN(1,1)*clhs127 + DN(1,2)*clhs129 + clhs10*clhs235 + clhs215;
            lhs(5,6)=DN(1,0)*clhs57 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs237;
            lhs(5,7)=DN(1,1)*clhs219;
            lhs(5,8)=DN(1,0)*clhs64 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs238;
            lhs(5,9)=DN(1,0)*clhs76 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs224 + clhs240;
            lhs(5,10)=DN(1,0)*clhs83 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs241;
            lhs(5,11)=DN(2,1)*clhs213 - clhs242;
            lhs(5,12)=DN(1,0)*clhs90 + DN(1,1)*clhs150 + DN(1,2)*clhs151 + clhs243;
            lhs(5,13)=DN(1,0)*clhs101 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs231 + clhs245;
            lhs(5,14)=DN(1,0)*clhs108 + DN(1,1)*clhs158 + DN(1,2)*clhs160 + clhs246;
            lhs(5,15)=DN(3,1)*clhs213 - clhs247;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs114 + DN(1,2)*clhs163 + clhs60;
            lhs(6,1)=DN(1,0)*clhs24 + DN(1,1)*clhs117 + DN(1,2)*clhs164 + clhs135;
            lhs(6,2)=DN(1,0)*clhs31 + DN(1,1)*clhs121 + DN(1,2)*clhs165 + clhs174 + clhs212;
            lhs(6,3)=DN(0,2)*clhs213 - clhs198;
            lhs(6,4)=DN(1,0)*clhs40 + DN(1,1)*clhs125 + DN(1,2)*clhs167 + clhs218;
            lhs(6,5)=DN(1,0)*clhs53 + DN(1,1)*clhs129 + DN(1,2)*clhs170 + clhs237;
            lhs(6,6)=DN(1,0)*clhs59 + DN(1,1)*clhs134 + DN(1,2)*clhs172 + clhs10*clhs248 + clhs215;
            lhs(6,7)=DN(1,2)*clhs219;
            lhs(6,8)=DN(1,0)*clhs66 + DN(1,1)*clhs138 + DN(1,2)*clhs176 + clhs250;
            lhs(6,9)=DN(1,0)*clhs79 + DN(1,1)*clhs142 + DN(1,2)*clhs178 + clhs251;
            lhs(6,10)=DN(1,0)*clhs85 + DN(1,1)*clhs147 + DN(1,2)*clhs180 + clhs224 + clhs253;
            lhs(6,11)=DN(2,2)*clhs213 - clhs254;
            lhs(6,12)=DN(1,0)*clhs92 + DN(1,1)*clhs151 + DN(1,2)*clhs184 + clhs255;
            lhs(6,13)=DN(1,0)*clhs104 + DN(1,1)*clhs155 + DN(1,2)*clhs186 + clhs256;
            lhs(6,14)=DN(1,0)*clhs110 + DN(1,1)*clhs160 + DN(1,2)*clhs188 + clhs231 + clhs258;
            lhs(6,15)=DN(3,2)*clhs213 - clhs259;
            lhs(7,0)=DN(1,0)*clhs192 + clhs61;
            lhs(7,1)=DN(1,1)*clhs192 + clhs136;
            lhs(7,2)=DN(1,2)*clhs192 + clhs175;
            lhs(7,3)=clhs199;
            lhs(7,4)=DN(1,0)*clhs260;
            lhs(7,5)=DN(1,1)*clhs260;
            lhs(7,6)=DN(1,2)*clhs260;
            lhs(7,7)=clhs194*(clhs214 + clhs235 + clhs248);
            lhs(7,8)=DN(1,0)*clhs201 + clhs261;
            lhs(7,9)=DN(1,1)*clhs201 + clhs262;
            lhs(7,10)=DN(1,2)*clhs201 + clhs263;
            lhs(7,11)=clhs264;
            lhs(7,12)=DN(1,0)*clhs206 + clhs265;
            lhs(7,13)=DN(1,1)*clhs206 + clhs266;
            lhs(7,14)=DN(1,2)*clhs206 + clhs267;
            lhs(7,15)=clhs268;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs271 + clhs68;
            lhs(8,1)=DN(2,0)*clhs19 + DN(2,1)*clhs21 + DN(2,2)*clhs24 + clhs139;
            lhs(8,2)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + DN(2,2)*clhs31 + clhs177;
            lhs(8,3)=DN(0,0)*clhs272 - clhs200;
            lhs(8,4)=DN(2,0)*clhs36 + DN(2,1)*clhs38 + DN(2,2)*clhs40 + clhs221 + clhs273;
            lhs(8,5)=DN(2,0)*clhs48 + DN(2,1)*clhs50 + DN(2,2)*clhs53 + clhs238;
            lhs(8,6)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs59 + clhs250;
            lhs(8,7)=DN(1,0)*clhs272 - clhs261;
            lhs(8,8)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs10*clhs274 + clhs275;
            lhs(8,9)=DN(2,0)*clhs74 + DN(2,1)*clhs76 + DN(2,2)*clhs79 + clhs277;
            lhs(8,10)=DN(2,0)*clhs81 + DN(2,1)*clhs83 + DN(2,2)*clhs85 + clhs278;
            lhs(8,11)=DN(2,0)*clhs279;
            lhs(8,12)=DN(2,0)*clhs88 + DN(2,1)*clhs90 + DN(2,2)*clhs92 + clhs281 + clhs284;
            lhs(8,13)=DN(2,0)*clhs99 + DN(2,1)*clhs101 + DN(2,2)*clhs104 + clhs285;
            lhs(8,14)=DN(2,0)*clhs106 + DN(2,1)*clhs108 + DN(2,2)*clhs110 + clhs286;
            lhs(8,15)=DN(3,0)*clhs272 - clhs287;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs113 + DN(2,2)*clhs114 + clhs80;
            lhs(9,1)=DN(2,0)*clhs21 + DN(2,1)*clhs115 + DN(2,2)*clhs117 + clhs144 + clhs271;
            lhs(9,2)=DN(2,0)*clhs29 + DN(2,1)*clhs119 + DN(2,2)*clhs121 + clhs179;
            lhs(9,3)=DN(0,1)*clhs272 - clhs202;
            lhs(9,4)=DN(2,0)*clhs38 + DN(2,1)*clhs124 + DN(2,2)*clhs125 + clhs225;
            lhs(9,5)=DN(2,0)*clhs50 + DN(2,1)*clhs127 + DN(2,2)*clhs129 + clhs240 + clhs273;
            lhs(9,6)=DN(2,0)*clhs57 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs251;
            lhs(9,7)=DN(1,1)*clhs272 - clhs262;
            lhs(9,8)=DN(2,0)*clhs64 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs277;
            lhs(9,9)=DN(2,0)*clhs76 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs10*clhs288 + clhs275;
            lhs(9,10)=DN(2,0)*clhs83 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs290;
            lhs(9,11)=DN(2,1)*clhs279;
            lhs(9,12)=DN(2,0)*clhs90 + DN(2,1)*clhs150 + DN(2,2)*clhs151 + clhs291;
            lhs(9,13)=DN(2,0)*clhs101 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs284 + clhs293;
            lhs(9,14)=DN(2,0)*clhs108 + DN(2,1)*clhs158 + DN(2,2)*clhs160 + clhs294;
            lhs(9,15)=DN(3,1)*clhs272 - clhs295;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs114 + DN(2,2)*clhs163 + clhs86;
            lhs(10,1)=DN(2,0)*clhs24 + DN(2,1)*clhs117 + DN(2,2)*clhs164 + clhs148;
            lhs(10,2)=DN(2,0)*clhs31 + DN(2,1)*clhs121 + DN(2,2)*clhs165 + clhs182 + clhs271;
            lhs(10,3)=DN(0,2)*clhs272 - clhs203;
            lhs(10,4)=DN(2,0)*clhs40 + DN(2,1)*clhs125 + DN(2,2)*clhs167 + clhs226;
            lhs(10,5)=DN(2,0)*clhs53 + DN(2,1)*clhs129 + DN(2,2)*clhs170 + clhs241;
            lhs(10,6)=DN(2,0)*clhs59 + DN(2,1)*clhs134 + DN(2,2)*clhs172 + clhs253 + clhs273;
            lhs(10,7)=DN(1,2)*clhs272 - clhs263;
            lhs(10,8)=DN(2,0)*clhs66 + DN(2,1)*clhs138 + DN(2,2)*clhs176 + clhs278;
            lhs(10,9)=DN(2,0)*clhs79 + DN(2,1)*clhs142 + DN(2,2)*clhs178 + clhs290;
            lhs(10,10)=DN(2,0)*clhs85 + DN(2,1)*clhs147 + DN(2,2)*clhs180 + clhs10*clhs296 + clhs275;
            lhs(10,11)=DN(2,2)*clhs279;
            lhs(10,12)=DN(2,0)*clhs92 + DN(2,1)*clhs151 + DN(2,2)*clhs184 + clhs298;
            lhs(10,13)=DN(2,0)*clhs104 + DN(2,1)*clhs155 + DN(2,2)*clhs186 + clhs299;
            lhs(10,14)=DN(2,0)*clhs110 + DN(2,1)*clhs160 + DN(2,2)*clhs188 + clhs284 + clhs301;
            lhs(10,15)=DN(3,2)*clhs272 - clhs302;
            lhs(11,0)=DN(2,0)*clhs192 + clhs87;
            lhs(11,1)=DN(2,1)*clhs192 + clhs149;
            lhs(11,2)=DN(2,2)*clhs192 + clhs183;
            lhs(11,3)=clhs204;
            lhs(11,4)=DN(2,0)*clhs196 + clhs227;
            lhs(11,5)=DN(2,1)*clhs196 + clhs242;
            lhs(11,6)=DN(2,2)*clhs196 + clhs254;
            lhs(11,7)=clhs264;
            lhs(11,8)=DN(2,0)*clhs303;
            lhs(11,9)=DN(2,1)*clhs303;
            lhs(11,10)=DN(2,2)*clhs303;
            lhs(11,11)=clhs194*(clhs274 + clhs288 + clhs296);
            lhs(11,12)=DN(2,0)*clhs206 + clhs304;
            lhs(11,13)=DN(2,1)*clhs206 + clhs305;
            lhs(11,14)=DN(2,2)*clhs206 + clhs306;
            lhs(11,15)=clhs307;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs309 + clhs94;
            lhs(12,1)=DN(3,0)*clhs19 + DN(3,1)*clhs21 + DN(3,2)*clhs24 + clhs152;
            lhs(12,2)=DN(3,0)*clhs27 + DN(3,1)*clhs29 + DN(3,2)*clhs31 + clhs185;
            lhs(12,3)=DN(0,0)*clhs310 - clhs205;
            lhs(12,4)=DN(3,0)*clhs36 + DN(3,1)*clhs38 + DN(3,2)*clhs40 + clhs229 + clhs311;
            lhs(12,5)=DN(3,0)*clhs48 + DN(3,1)*clhs50 + DN(3,2)*clhs53 + clhs243;
            lhs(12,6)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs59 + clhs255;
            lhs(12,7)=DN(1,0)*clhs310 - clhs265;
            lhs(12,8)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs281 + clhs312;
            lhs(12,9)=DN(3,0)*clhs74 + DN(3,1)*clhs76 + DN(3,2)*clhs79 + clhs291;
            lhs(12,10)=DN(3,0)*clhs81 + DN(3,1)*clhs83 + DN(3,2)*clhs85 + clhs298;
            lhs(12,11)=DN(2,0)*clhs310 - clhs304;
            lhs(12,12)=DN(3,0)*clhs88 + DN(3,1)*clhs90 + DN(3,2)*clhs92 + clhs10*clhs313 + clhs314;
            lhs(12,13)=DN(3,0)*clhs99 + DN(3,1)*clhs101 + DN(3,2)*clhs104 + clhs316;
            lhs(12,14)=DN(3,0)*clhs106 + DN(3,1)*clhs108 + DN(3,2)*clhs110 + clhs317;
            lhs(12,15)=DN(3,0)*clhs318;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs113 + DN(3,2)*clhs114 + clhs105;
            lhs(13,1)=DN(3,0)*clhs21 + DN(3,1)*clhs115 + DN(3,2)*clhs117 + clhs157 + clhs309;
            lhs(13,2)=DN(3,0)*clhs29 + DN(3,1)*clhs119 + DN(3,2)*clhs121 + clhs187;
            lhs(13,3)=DN(0,1)*clhs310 - clhs207;
            lhs(13,4)=DN(3,0)*clhs38 + DN(3,1)*clhs124 + DN(3,2)*clhs125 + clhs232;
            lhs(13,5)=DN(3,0)*clhs50 + DN(3,1)*clhs127 + DN(3,2)*clhs129 + clhs245 + clhs311;
            lhs(13,6)=DN(3,0)*clhs57 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs256;
            lhs(13,7)=DN(1,1)*clhs310 - clhs266;
            lhs(13,8)=DN(3,0)*clhs64 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs285;
            lhs(13,9)=DN(3,0)*clhs76 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs293 + clhs312;
            lhs(13,10)=DN(3,0)*clhs83 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs299;
            lhs(13,11)=DN(2,1)*clhs310 - clhs305;
            lhs(13,12)=DN(3,0)*clhs90 + DN(3,1)*clhs150 + DN(3,2)*clhs151 + clhs316;
            lhs(13,13)=DN(3,0)*clhs101 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs10*clhs319 + clhs314;
            lhs(13,14)=DN(3,0)*clhs108 + DN(3,1)*clhs158 + DN(3,2)*clhs160 + clhs320;
            lhs(13,15)=DN(3,1)*clhs318;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs114 + DN(3,2)*clhs163 + clhs111;
            lhs(14,1)=DN(3,0)*clhs24 + DN(3,1)*clhs117 + DN(3,2)*clhs164 + clhs161;
            lhs(14,2)=DN(3,0)*clhs31 + DN(3,1)*clhs121 + DN(3,2)*clhs165 + clhs190 + clhs309;
            lhs(14,3)=DN(0,2)*clhs310 - clhs208;
            lhs(14,4)=DN(3,0)*clhs40 + DN(3,1)*clhs125 + DN(3,2)*clhs167 + clhs233;
            lhs(14,5)=DN(3,0)*clhs53 + DN(3,1)*clhs129 + DN(3,2)*clhs170 + clhs246;
            lhs(14,6)=DN(3,0)*clhs59 + DN(3,1)*clhs134 + DN(3,2)*clhs172 + clhs258 + clhs311;
            lhs(14,7)=DN(1,2)*clhs310 - clhs267;
            lhs(14,8)=DN(3,0)*clhs66 + DN(3,1)*clhs138 + DN(3,2)*clhs176 + clhs286;
            lhs(14,9)=DN(3,0)*clhs79 + DN(3,1)*clhs142 + DN(3,2)*clhs178 + clhs294;
            lhs(14,10)=DN(3,0)*clhs85 + DN(3,1)*clhs147 + DN(3,2)*clhs180 + clhs301 + clhs312;
            lhs(14,11)=DN(2,2)*clhs310 - clhs306;
            lhs(14,12)=DN(3,0)*clhs92 + DN(3,1)*clhs151 + DN(3,2)*clhs184 + clhs317;
            lhs(14,13)=DN(3,0)*clhs104 + DN(3,1)*clhs155 + DN(3,2)*clhs186 + clhs320;
            lhs(14,14)=DN(3,0)*clhs110 + DN(3,1)*clhs160 + DN(3,2)*clhs188 + clhs10*clhs321 + clhs314;
            lhs(14,15)=DN(3,2)*clhs318;
            lhs(15,0)=DN(3,0)*clhs192 + clhs112;
            lhs(15,1)=DN(3,1)*clhs192 + clhs162;
            lhs(15,2)=DN(3,2)*clhs192 + clhs191;
            lhs(15,3)=clhs209;
            lhs(15,4)=DN(3,0)*clhs196 + clhs234;
            lhs(15,5)=DN(3,1)*clhs196 + clhs247;
            lhs(15,6)=DN(3,2)*clhs196 + clhs259;
            lhs(15,7)=clhs268;
            lhs(15,8)=DN(3,0)*clhs201 + clhs287;
            lhs(15,9)=DN(3,1)*clhs201 + clhs295;
            lhs(15,10)=DN(3,2)*clhs201 + clhs302;
            lhs(15,11)=clhs307;
            lhs(15,12)=DN(3,0)*clhs322;
            lhs(15,13)=DN(3,1)*clhs322;
            lhs(15,14)=DN(3,2)*clhs322;
            lhs(15,15)=clhs194*(clhs313 + clhs319 + clhs321);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesData<2, 3> &rData,
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
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 =             rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs8 =             crhs3 + crhs7;
const double crhs9 =             rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs10 =             crhs8*(crhs9*h/stab_c1 + mu);
const double crhs11 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs12 =             1.0/(crhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs13 =             1.0*crhs12*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs6);
const double crhs14 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs15 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs16 =             rho*(crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs5*crhs7);
const double crhs17 =             1.0*crhs12*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs14 + crhs15 + crhs16);
const double crhs18 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs19 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs10 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs6 - crhs11*crhs13;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs10 - DN(0,1)*stress[1] + N[0]*crhs14 - N[0]*crhs15 - N[0]*crhs16 - crhs11*crhs17;
            rhs[2]=-DN(0,0)*crhs13 - DN(0,1)*crhs17 - N[0]*crhs8;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs10 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs6 - crhs13*crhs18;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs10 - DN(1,1)*stress[1] + N[1]*crhs14 - N[1]*crhs15 - N[1]*crhs16 - crhs17*crhs18;
            rhs[5]=-DN(1,0)*crhs13 - DN(1,1)*crhs17 - N[1]*crhs8;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs10 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs6 - crhs13*crhs19;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs10 - DN(2,1)*stress[1] + N[2]*crhs14 - N[2]*crhs15 - N[2]*crhs16 - crhs17*crhs19;
            rhs[8]=-DN(2,0)*crhs13 - DN(2,1)*crhs17 - N[2]*crhs8;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesData<3, 4> &rData,
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
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 =             DN(0,0)*v(0,0);
const double crhs4 =             DN(1,0)*v(1,0);
const double crhs5 =             DN(2,0)*v(2,0);
const double crhs6 =             DN(3,0)*v(3,0);
const double crhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs10 =             rho*(crhs7*(crhs3 + crhs4 + crhs5 + crhs6) + crhs8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs11 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs12 =             DN(0,1)*v(0,1);
const double crhs13 =             DN(1,1)*v(1,1);
const double crhs14 =             DN(2,1)*v(2,1);
const double crhs15 =             DN(3,1)*v(3,1);
const double crhs16 =             crhs11 + crhs12 + crhs13 + crhs14 + crhs15 + crhs3 + crhs4 + crhs5 + crhs6;
const double crhs17 =             rho*stab_c2*sqrt(pow(crhs7, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs18 =             crhs16*(crhs17*h/stab_c1 + mu);
const double crhs19 =             rho*(DN(0,0)*crhs7 + DN(0,1)*crhs8 + DN(0,2)*crhs9);
const double crhs20 =             1.0/(crhs17/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs21 =             1.0*crhs20*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs10 + crhs2);
const double crhs22 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs23 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs24 =             rho*(crhs7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs8*(crhs12 + crhs13 + crhs14 + crhs15) + crhs9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs25 =             1.0*crhs20*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs22 + crhs23 + crhs24);
const double crhs26 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs27 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs28 =             rho*(crhs11*crhs9 + crhs7*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs29 =             1.0*crhs20*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs26 + crhs27 + crhs28);
const double crhs30 =             rho*(DN(1,0)*crhs7 + DN(1,1)*crhs8 + DN(1,2)*crhs9);
const double crhs31 =             rho*(DN(2,0)*crhs7 + DN(2,1)*crhs8 + DN(2,2)*crhs9);
const double crhs32 =             rho*(DN(3,0)*crhs7 + DN(3,1)*crhs8 + DN(3,2)*crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - crhs19*crhs21;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs18 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs22 - N[0]*crhs23 - N[0]*crhs24 - crhs19*crhs25;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs18 - DN(0,2)*stress[2] + N[0]*crhs26 - N[0]*crhs27 - N[0]*crhs28 - crhs19*crhs29;
            rhs[3]=-DN(0,0)*crhs21 - DN(0,1)*crhs25 - DN(0,2)*crhs29 - N[0]*crhs16;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - crhs21*crhs30;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs18 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs22 - N[1]*crhs23 - N[1]*crhs24 - crhs25*crhs30;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs18 - DN(1,2)*stress[2] + N[1]*crhs26 - N[1]*crhs27 - N[1]*crhs28 - crhs29*crhs30;
            rhs[7]=-DN(1,0)*crhs21 - DN(1,1)*crhs25 - DN(1,2)*crhs29 - N[1]*crhs16;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - crhs21*crhs31;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs18 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs22 - N[2]*crhs23 - N[2]*crhs24 - crhs25*crhs31;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs18 - DN(2,2)*stress[2] + N[2]*crhs26 - N[2]*crhs27 - N[2]*crhs28 - crhs29*crhs31;
            rhs[11]=-DN(2,0)*crhs21 - DN(2,1)*crhs25 - DN(2,2)*crhs29 - N[2]*crhs16;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs18 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs10 - N[3]*crhs2 - crhs21*crhs32;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs18 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs22 - N[3]*crhs23 - N[3]*crhs24 - crhs25*crhs32;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs18 - DN(3,2)*stress[2] + N[3]*crhs26 - N[3]*crhs27 - N[3]*crhs28 - crhs29*crhs32;
            rhs[15]=-DN(3,0)*crhs21 - DN(3,1)*crhs25 - DN(3,2)*crhs29 - N[3]*crhs16;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
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

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 =             DN(0,0)*cV0 + DN(0,1)*cV1;
const double cV3 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 =             1.0*DNenr(0,0)*cV3*rho;
const double cV5 =             1.0*DNenr(1,0)*cV3*rho;
const double cV6 =             1.0*DNenr(2,0)*cV3*rho;
const double cV7 =             1.0*DNenr(0,1)*cV3*rho;
const double cV8 =             1.0*DNenr(1,1)*cV3*rho;
const double cV9 =             1.0*DNenr(2,1)*cV3*rho;
const double cV10 =             1.0*cV3;
const double cV11 =             DN(1,0)*cV0 + DN(1,1)*cV1;
const double cV12 =             DN(2,0)*cV0 + DN(2,1)*cV1;
            V(0,0)=-DN(0,0)*Nenr[0] + cV2*cV4;
            V(0,1)=-DN(0,0)*Nenr[1] + cV2*cV5;
            V(0,2)=-DN(0,0)*Nenr[2] + cV2*cV6;
            V(1,0)=-DN(0,1)*Nenr[0] + cV2*cV7;
            V(1,1)=-DN(0,1)*Nenr[1] + cV2*cV8;
            V(1,2)=-DN(0,1)*Nenr[2] + cV2*cV9;
            V(2,0)=cV10*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV10*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV10*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] + cV11*cV4;
            V(3,1)=-DN(1,0)*Nenr[1] + cV11*cV5;
            V(3,2)=-DN(1,0)*Nenr[2] + cV11*cV6;
            V(4,0)=-DN(1,1)*Nenr[0] + cV11*cV7;
            V(4,1)=-DN(1,1)*Nenr[1] + cV11*cV8;
            V(4,2)=-DN(1,1)*Nenr[2] + cV11*cV9;
            V(5,0)=cV10*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV10*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV10*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] + cV12*cV4;
            V(6,1)=-DN(2,0)*Nenr[1] + cV12*cV5;
            V(6,2)=-DN(2,0)*Nenr[2] + cV12*cV6;
            V(7,0)=-DN(2,1)*Nenr[0] + cV12*cV7;
            V(7,1)=-DN(2,1)*Nenr[1] + cV12*cV8;
            V(7,2)=-DN(2,1)*Nenr[2] + cV12*cV9;
            V(8,0)=cV10*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV10*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV10*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 =             DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0;
const double cH3 =             1.0/(rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 =             1.0*DNenr(0,0)*cH3*rho;
const double cH5 =             1.0*DNenr(0,1)*cH3*rho;
const double cH6 =             1.0*cH3;
const double cH7 =             DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0;
const double cH8 =             DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0;
const double cH9 =             1.0*DNenr(1,0)*cH3*rho;
const double cH10 =             1.0*DNenr(1,1)*cH3*rho;
const double cH11 =             1.0*DNenr(2,0)*cH3*rho;
const double cH12 =             1.0*DNenr(2,1)*cH3*rho;
            H(0,0)=DN(0,0)*Nenr[0] + cH2*cH4;
            H(0,1)=DN(0,1)*Nenr[0] + cH2*cH5;
            H(0,2)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DN(1,0)*Nenr[0] + cH4*cH7;
            H(0,4)=DN(1,1)*Nenr[0] + cH5*cH7;
            H(0,5)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DN(2,0)*Nenr[0] + cH4*cH8;
            H(0,7)=DN(2,1)*Nenr[0] + cH5*cH8;
            H(0,8)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DN(0,0)*Nenr[1] + cH2*cH9;
            H(1,1)=DN(0,1)*Nenr[1] + cH10*cH2;
            H(1,2)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DN(1,0)*Nenr[1] + cH7*cH9;
            H(1,4)=DN(1,1)*Nenr[1] + cH10*cH7;
            H(1,5)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DN(2,0)*Nenr[1] + cH8*cH9;
            H(1,7)=DN(2,1)*Nenr[1] + cH10*cH8;
            H(1,8)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DN(0,0)*Nenr[2] + cH11*cH2;
            H(2,1)=DN(0,1)*Nenr[2] + cH12*cH2;
            H(2,2)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DN(1,0)*Nenr[2] + cH11*cH7;
            H(2,4)=DN(1,1)*Nenr[2] + cH12*cH7;
            H(2,5)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DN(2,0)*Nenr[2] + cH11*cH8;
            H(2,7)=DN(2,1)*Nenr[2] + cH12*cH8;
            H(2,8)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
            Kee(1,2)=cKee3;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee3;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs_ee1 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs_ee2 =             crhs_ee0 + crhs_ee1;
const double crhs_ee3 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee4 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee5 =             1.0/(rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee6 =             1.0*crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 =             1.0*crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee6 - DNenr(0,1)*crhs_ee7 - Nenr[0]*crhs_ee2;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee6 - DNenr(1,1)*crhs_ee7 - Nenr[1]*crhs_ee2;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee6 - DNenr(2,1)*crhs_ee7 - Nenr[2]*crhs_ee2;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
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

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 =             DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2;
const double cV4 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV5 =             1.0*DNenr(0,0)*cV4*rho;
const double cV6 =             1.0*DNenr(1,0)*cV4*rho;
const double cV7 =             1.0*DNenr(2,0)*cV4*rho;
const double cV8 =             1.0*DNenr(3,0)*cV4*rho;
const double cV9 =             1.0*DNenr(0,1)*cV4*rho;
const double cV10 =             1.0*DNenr(1,1)*cV4*rho;
const double cV11 =             1.0*DNenr(2,1)*cV4*rho;
const double cV12 =             1.0*DNenr(3,1)*cV4*rho;
const double cV13 =             1.0*DNenr(0,2)*cV4*rho;
const double cV14 =             1.0*DNenr(1,2)*cV4*rho;
const double cV15 =             1.0*DNenr(2,2)*cV4*rho;
const double cV16 =             1.0*DNenr(3,2)*cV4*rho;
const double cV17 =             1.0*cV4;
const double cV18 =             DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2;
const double cV19 =             DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2;
const double cV20 =             DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2;
            V(0,0)=-DN(0,0)*Nenr[0] + cV3*cV5;
            V(0,1)=-DN(0,0)*Nenr[1] + cV3*cV6;
            V(0,2)=-DN(0,0)*Nenr[2] + cV3*cV7;
            V(0,3)=-DN(0,0)*Nenr[3] + cV3*cV8;
            V(1,0)=-DN(0,1)*Nenr[0] + cV3*cV9;
            V(1,1)=-DN(0,1)*Nenr[1] + cV10*cV3;
            V(1,2)=-DN(0,1)*Nenr[2] + cV11*cV3;
            V(1,3)=-DN(0,1)*Nenr[3] + cV12*cV3;
            V(2,0)=-DN(0,2)*Nenr[0] + cV13*cV3;
            V(2,1)=-DN(0,2)*Nenr[1] + cV14*cV3;
            V(2,2)=-DN(0,2)*Nenr[2] + cV15*cV3;
            V(2,3)=-DN(0,2)*Nenr[3] + cV16*cV3;
            V(3,0)=cV17*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV17*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV17*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV17*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] + cV18*cV5;
            V(4,1)=-DN(1,0)*Nenr[1] + cV18*cV6;
            V(4,2)=-DN(1,0)*Nenr[2] + cV18*cV7;
            V(4,3)=-DN(1,0)*Nenr[3] + cV18*cV8;
            V(5,0)=-DN(1,1)*Nenr[0] + cV18*cV9;
            V(5,1)=-DN(1,1)*Nenr[1] + cV10*cV18;
            V(5,2)=-DN(1,1)*Nenr[2] + cV11*cV18;
            V(5,3)=-DN(1,1)*Nenr[3] + cV12*cV18;
            V(6,0)=-DN(1,2)*Nenr[0] + cV13*cV18;
            V(6,1)=-DN(1,2)*Nenr[1] + cV14*cV18;
            V(6,2)=-DN(1,2)*Nenr[2] + cV15*cV18;
            V(6,3)=-DN(1,2)*Nenr[3] + cV16*cV18;
            V(7,0)=cV17*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV17*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV17*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV17*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] + cV19*cV5;
            V(8,1)=-DN(2,0)*Nenr[1] + cV19*cV6;
            V(8,2)=-DN(2,0)*Nenr[2] + cV19*cV7;
            V(8,3)=-DN(2,0)*Nenr[3] + cV19*cV8;
            V(9,0)=-DN(2,1)*Nenr[0] + cV19*cV9;
            V(9,1)=-DN(2,1)*Nenr[1] + cV10*cV19;
            V(9,2)=-DN(2,1)*Nenr[2] + cV11*cV19;
            V(9,3)=-DN(2,1)*Nenr[3] + cV12*cV19;
            V(10,0)=-DN(2,2)*Nenr[0] + cV13*cV19;
            V(10,1)=-DN(2,2)*Nenr[1] + cV14*cV19;
            V(10,2)=-DN(2,2)*Nenr[2] + cV15*cV19;
            V(10,3)=-DN(2,2)*Nenr[3] + cV16*cV19;
            V(11,0)=cV17*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV17*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV17*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV17*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] + cV20*cV5;
            V(12,1)=-DN(3,0)*Nenr[1] + cV20*cV6;
            V(12,2)=-DN(3,0)*Nenr[2] + cV20*cV7;
            V(12,3)=-DN(3,0)*Nenr[3] + cV20*cV8;
            V(13,0)=-DN(3,1)*Nenr[0] + cV20*cV9;
            V(13,1)=-DN(3,1)*Nenr[1] + cV10*cV20;
            V(13,2)=-DN(3,1)*Nenr[2] + cV11*cV20;
            V(13,3)=-DN(3,1)*Nenr[3] + cV12*cV20;
            V(14,0)=-DN(3,2)*Nenr[0] + cV13*cV20;
            V(14,1)=-DN(3,2)*Nenr[1] + cV14*cV20;
            V(14,2)=-DN(3,2)*Nenr[2] + cV15*cV20;
            V(14,3)=-DN(3,2)*Nenr[3] + cV16*cV20;
            V(15,0)=cV17*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV17*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV17*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV17*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 =             DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0;
const double cH4 =             1.0/(rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH5 =             1.0*DNenr(0,0)*cH4*rho;
const double cH6 =             1.0*DNenr(0,1)*cH4*rho;
const double cH7 =             1.0*DNenr(0,2)*cH4*rho;
const double cH8 =             1.0*cH4;
const double cH9 =             DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0;
const double cH10 =             DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0;
const double cH11 =             DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0;
const double cH12 =             1.0*DNenr(1,0)*cH4*rho;
const double cH13 =             1.0*DNenr(1,1)*cH4*rho;
const double cH14 =             1.0*DNenr(1,2)*cH4*rho;
const double cH15 =             1.0*DNenr(2,0)*cH4*rho;
const double cH16 =             1.0*DNenr(2,1)*cH4*rho;
const double cH17 =             1.0*DNenr(2,2)*cH4*rho;
const double cH18 =             1.0*DNenr(3,0)*cH4*rho;
const double cH19 =             1.0*DNenr(3,1)*cH4*rho;
const double cH20 =             1.0*DNenr(3,2)*cH4*rho;
            H(0,0)=DN(0,0)*Nenr[0] + cH3*cH5;
            H(0,1)=DN(0,1)*Nenr[0] + cH3*cH6;
            H(0,2)=DN(0,2)*Nenr[0] + cH3*cH7;
            H(0,3)=cH8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DN(1,0)*Nenr[0] + cH5*cH9;
            H(0,5)=DN(1,1)*Nenr[0] + cH6*cH9;
            H(0,6)=DN(1,2)*Nenr[0] + cH7*cH9;
            H(0,7)=cH8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DN(2,0)*Nenr[0] + cH10*cH5;
            H(0,9)=DN(2,1)*Nenr[0] + cH10*cH6;
            H(0,10)=DN(2,2)*Nenr[0] + cH10*cH7;
            H(0,11)=cH8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DN(3,0)*Nenr[0] + cH11*cH5;
            H(0,13)=DN(3,1)*Nenr[0] + cH11*cH6;
            H(0,14)=DN(3,2)*Nenr[0] + cH11*cH7;
            H(0,15)=cH8*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DN(0,0)*Nenr[1] + cH12*cH3;
            H(1,1)=DN(0,1)*Nenr[1] + cH13*cH3;
            H(1,2)=DN(0,2)*Nenr[1] + cH14*cH3;
            H(1,3)=cH8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DN(1,0)*Nenr[1] + cH12*cH9;
            H(1,5)=DN(1,1)*Nenr[1] + cH13*cH9;
            H(1,6)=DN(1,2)*Nenr[1] + cH14*cH9;
            H(1,7)=cH8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DN(2,0)*Nenr[1] + cH10*cH12;
            H(1,9)=DN(2,1)*Nenr[1] + cH10*cH13;
            H(1,10)=DN(2,2)*Nenr[1] + cH10*cH14;
            H(1,11)=cH8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DN(3,0)*Nenr[1] + cH11*cH12;
            H(1,13)=DN(3,1)*Nenr[1] + cH11*cH13;
            H(1,14)=DN(3,2)*Nenr[1] + cH11*cH14;
            H(1,15)=cH8*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DN(0,0)*Nenr[2] + cH15*cH3;
            H(2,1)=DN(0,1)*Nenr[2] + cH16*cH3;
            H(2,2)=DN(0,2)*Nenr[2] + cH17*cH3;
            H(2,3)=cH8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DN(1,0)*Nenr[2] + cH15*cH9;
            H(2,5)=DN(1,1)*Nenr[2] + cH16*cH9;
            H(2,6)=DN(1,2)*Nenr[2] + cH17*cH9;
            H(2,7)=cH8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DN(2,0)*Nenr[2] + cH10*cH15;
            H(2,9)=DN(2,1)*Nenr[2] + cH10*cH16;
            H(2,10)=DN(2,2)*Nenr[2] + cH10*cH17;
            H(2,11)=cH8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DN(3,0)*Nenr[2] + cH11*cH15;
            H(2,13)=DN(3,1)*Nenr[2] + cH11*cH16;
            H(2,14)=DN(3,2)*Nenr[2] + cH11*cH17;
            H(2,15)=cH8*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DN(0,0)*Nenr[3] + cH18*cH3;
            H(3,1)=DN(0,1)*Nenr[3] + cH19*cH3;
            H(3,2)=DN(0,2)*Nenr[3] + cH20*cH3;
            H(3,3)=cH8*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DN(1,0)*Nenr[3] + cH18*cH9;
            H(3,5)=DN(1,1)*Nenr[3] + cH19*cH9;
            H(3,6)=DN(1,2)*Nenr[3] + cH20*cH9;
            H(3,7)=cH8*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DN(2,0)*Nenr[3] + cH10*cH18;
            H(3,9)=DN(2,1)*Nenr[3] + cH10*cH19;
            H(3,10)=DN(2,2)*Nenr[3] + cH10*cH20;
            H(3,11)=cH8*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DN(3,0)*Nenr[3] + cH11*cH18;
            H(3,13)=DN(3,1)*Nenr[3] + cH11*cH19;
            H(3,14)=DN(3,2)*Nenr[3] + cH11*cH20;
            H(3,15)=cH8*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 =             cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 =             cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 =             cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2) + pow(DNenr(0,2), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(0,3)=cKee3;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2) + pow(DNenr(1,2), 2));
            Kee(1,2)=cKee4;
            Kee(1,3)=cKee5;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee4;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2) + pow(DNenr(2,2), 2));
            Kee(2,3)=cKee6;
            Kee(3,0)=cKee3;
            Kee(3,1)=cKee5;
            Kee(3,2)=cKee6;
            Kee(3,3)=cKee0*(pow(DNenr(3,0), 2) + pow(DNenr(3,1), 2) + pow(DNenr(3,2), 2));


    const double crhs_ee0 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee1 =             DN(0,0)*v(0,0);
const double crhs_ee2 =             DN(0,1)*v(0,1);
const double crhs_ee3 =             DN(1,0)*v(1,0);
const double crhs_ee4 =             DN(1,1)*v(1,1);
const double crhs_ee5 =             DN(2,0)*v(2,0);
const double crhs_ee6 =             DN(2,1)*v(2,1);
const double crhs_ee7 =             DN(3,0)*v(3,0);
const double crhs_ee8 =             DN(3,1)*v(3,1);
const double crhs_ee9 =             crhs_ee0 + crhs_ee1 + crhs_ee2 + crhs_ee3 + crhs_ee4 + crhs_ee5 + crhs_ee6 + crhs_ee7 + crhs_ee8;
const double crhs_ee10 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee11 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee12 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee13 =             1.0/(rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee11, 2) + pow(crhs_ee12, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee14 =             1.0*crhs_ee13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee10*(crhs_ee1 + crhs_ee3 + crhs_ee5 + crhs_ee7) + crhs_ee11*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee12*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee15 =             1.0*crhs_ee13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee10*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee11*(crhs_ee2 + crhs_ee4 + crhs_ee6 + crhs_ee8) + crhs_ee12*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee16 =             1.0*crhs_ee13*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee0*crhs_ee12 + crhs_ee10*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee11*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee14 - DNenr(0,1)*crhs_ee15 - DNenr(0,2)*crhs_ee16 - Nenr[0]*crhs_ee9;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee14 - DNenr(1,1)*crhs_ee15 - DNenr(1,2)*crhs_ee16 - Nenr[1]*crhs_ee9;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee14 - DNenr(2,1)*crhs_ee15 - DNenr(2,2)*crhs_ee16 - Nenr[2]*crhs_ee9;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee14 - DNenr(3,1)*crhs_ee15 - DNenr(3,2)*crhs_ee16 - Nenr[3]*crhs_ee9;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ComputeSplitting(
    TElementData &rData,
    MatrixType &rShapeFunctionsPos,
    MatrixType &rShapeFunctionsNeg,
    MatrixType &rEnrichedShapeFunctionsPos,
    MatrixType &rEnrichedShapeFunctionsNeg,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg)
{
    // Set the positive and negative enrichment interpolation matrices
    // Note that the enrichment is constructed using the standard shape functions such that:
    // In the negative distance region, the enrichment functions correspondig to the negative
    // distance nodes are null and the positive distance nodes are equal to the standard shape
    // functions. On the contrary, for the positive distance region, the enrichment functions
    // corresponding to the positive distance nodes are null meanwhile the negative distance
    // nodes are equal to the standard. This yields a discontinuous enrichment space.
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    if (Dim == 2)
        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    else
        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);

    // Call the positive side modified shape functions calculator
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (p_modified_sh_func->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CondenseEnrichment(
    const TElementData &rData,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    const double min_area_ratio = 1e-7;

    // Compute positive side, negative side and total volumes
    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double Vol = positive_volume + negative_volume;

    // Compute the maximum diagonal value in the enrichment stiffness matrix
    double max_diag = 0.0;
    for (unsigned int k = 0; k < NumNodes; ++k){
        if (std::abs(rKeeTot(k, k)) > max_diag){
            max_diag = std::abs(rKeeTot(k, k));
        }
    }
    if (max_diag == 0){
        max_diag = 1.0;
    }

    // Check that positive and negative volumes ratios are larger than the minimum
    // If not, substitute the enrichment term by
    if (positive_volume / Vol < min_area_ratio){
        for (unsigned int i = 0; i < NumNodes; ++i){
            if (rData.Distance[i] >= 0.0){
                rKeeTot(i, i) += 1000.0 * max_diag;
            }
        }
    }

    if (negative_volume / Vol < min_area_ratio){
        for (unsigned int i = 0; i < NumNodes; ++i){
            if (rData.Distance[i] < 0.0){
                rKeeTot(i, i) += 1000.0 * max_diag;
            }
        }
    }

    // "weakly" impose continuity
    for (unsigned int i = 0; i < Dim; ++i){
        const double di = std::abs(rData.Distance[i]);
        for (unsigned int j = i + 1; j < NumNodes; j++){
            const double dj = std::abs(rData.Distance[j]);
            // Check if the edge is cut, if it is, set the penalty constraint
            if (rData.Distance[i] * rData.Distance[j] < 0.0){
                double sum_d = di + dj;
                double Ni = dj / sum_d;
                double Nj = di / sum_d;
                double penalty_coeff = max_diag * 0.001; // h/BDFVector[0];
                rKeeTot(i, i) += penalty_coeff * Ni * Ni;
                rKeeTot(i, j) -= penalty_coeff * Ni * Nj;
                rKeeTot(j, i) -= penalty_coeff * Nj * Ni;
                rKeeTot(j, j) += penalty_coeff * Nj * Nj;
            }
        }
    }

    // Enrichment condensation (add to LHS and RHS the enrichment contributions)
    double det;
    MatrixType inverse_diag(NumNodes, NumNodes);
    MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

    const Matrix tmp = prod(inverse_diag, rHtot);
    noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

    const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
    noalias(rRightHandSideVector) -= prod(rVtot, tmp2);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>;

} // namespace Kratos