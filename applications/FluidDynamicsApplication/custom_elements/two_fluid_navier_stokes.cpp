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
    return Kratos::make_intrusive<TwoFluidNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokes>(NewId, pGeom, pProperties);
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
                for (unsigned int g = 0; g < number_of_gauss_points; g++){
                    UpdateIntegrationPointData(
                        data,
                        g,
                        gauss_weights[g],
                        row(shape_functions, g),
                        shape_derivatives[g]);
                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                }
            } else {
                MatrixType Vtot = ZeroMatrix(NumNodes * (Dim + 1), NumNodes);
                MatrixType Htot = ZeroMatrix(NumNodes, NumNodes * (Dim + 1));
                MatrixType Kee_tot = ZeroMatrix(NumNodes, NumNodes);
                VectorType rhs_ee_tot = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); g_pos++){
                    UpdateIntegrationPointData(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);

                    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); g_neg++){
                    UpdateIntegrationPointData(
                        data,
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);
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
                UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
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

    KRATOS_CHECK_VARIABLE_KEY( DIVERGENCE );

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
void TwoFluidNavierStokes<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
    rData.ComputeDarcyTerm();
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX,
    const typename TElementData::MatrixRowType& rNenr,
    const typename TElementData::ShapeDerivativesType& rDN_DXenr) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex,Weight,rN,rDN_DX,rNenr,rDN_DXenr);
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
    rData.ComputeDarcyTerm();
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
    const double K_darcy = rData.DarcyTerm;

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
const double clhs8 =             pow(N[0], 2);
const double clhs9 =             bdf0*rho;
const double clhs10 =             rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs11 =             K_darcy*N[0];
const double clhs12 =             N[0]*clhs9;
const double clhs13 =             clhs10 + clhs11 + clhs12;
const double clhs14 =             1.0*clhs11;
const double clhs15 =             1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs16 =             clhs14*clhs15;
const double clhs17 =             1.0*clhs10;
const double clhs18 =             clhs15*clhs17;
const double clhs19 =             K_darcy*clhs8 + N[0]*clhs10 - clhs13*clhs16 + clhs13*clhs18 + clhs8*clhs9;
const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
const double clhs21 =             C(1,2)*DN(0,1);
const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,0)*clhs7;
const double clhs24 =             DN(0,1)*clhs23;
const double clhs25 =             -N[0] - clhs16 + clhs18;
const double clhs26 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs27 =             C(0,2)*DN(1,0);
const double clhs28 =             C(2,2)*DN(1,1) + clhs27;
const double clhs29 =             DN(0,0)*DN(1,0);
const double clhs30 =             clhs29*clhs7;
const double clhs31 =             N[1]*clhs11;
const double clhs32 =             N[1]*clhs12;
const double clhs33 =             rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs34 =             K_darcy*N[1];
const double clhs35 =             N[1]*clhs9;
const double clhs36 =             clhs33 + clhs34 + clhs35;
const double clhs37 =             N[0]*clhs33 - clhs16*clhs36 + clhs18*clhs36 + clhs31 + clhs32;
const double clhs38 =             C(0,1)*DN(1,1) + clhs27;
const double clhs39 =             C(1,2)*DN(1,1);
const double clhs40 =             C(2,2)*DN(1,0) + clhs39;
const double clhs41 =             DN(1,1)*clhs23;
const double clhs42 =             DN(0,0)*N[1];
const double clhs43 =             DN(1,0)*N[0];
const double clhs44 =             1.0*K_darcy*clhs15;
const double clhs45 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs46 =             C(0,2)*DN(2,0);
const double clhs47 =             C(2,2)*DN(2,1) + clhs46;
const double clhs48 =             DN(0,0)*DN(2,0);
const double clhs49 =             clhs48*clhs7;
const double clhs50 =             N[2]*clhs11;
const double clhs51 =             N[2]*clhs12;
const double clhs52 =             rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs53 =             K_darcy*N[2];
const double clhs54 =             N[2]*clhs9;
const double clhs55 =             clhs52 + clhs53 + clhs54;
const double clhs56 =             N[0]*clhs52 - clhs16*clhs55 + clhs18*clhs55 + clhs50 + clhs51;
const double clhs57 =             C(0,1)*DN(2,1) + clhs46;
const double clhs58 =             C(1,2)*DN(2,1);
const double clhs59 =             C(2,2)*DN(2,0) + clhs58;
const double clhs60 =             DN(2,1)*clhs23;
const double clhs61 =             DN(0,0)*N[2];
const double clhs62 =             DN(2,0)*N[0];
const double clhs63 =             C(0,1)*DN(0,0) + clhs21;
const double clhs64 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs65 =             pow(DN(0,1), 2);
const double clhs66 =             C(0,1)*DN(1,0) + clhs39;
const double clhs67 =             DN(0,1)*clhs7;
const double clhs68 =             DN(1,0)*clhs67;
const double clhs69 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs70 =             DN(0,1)*DN(1,1);
const double clhs71 =             clhs7*clhs70;
const double clhs72 =             DN(0,1)*N[1];
const double clhs73 =             DN(1,1)*N[0];
const double clhs74 =             C(0,1)*DN(2,0) + clhs58;
const double clhs75 =             DN(2,0)*clhs67;
const double clhs76 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs77 =             DN(0,1)*DN(2,1);
const double clhs78 =             clhs7*clhs77;
const double clhs79 =             DN(0,1)*N[2];
const double clhs80 =             DN(2,1)*N[0];
const double clhs81 =             N[0] + clhs15*(1.0*clhs12 + clhs14 + clhs17);
const double clhs82 =             1.0*clhs15;
const double clhs83 =             1.0*DN(0,0)*clhs15;
const double clhs84 =             1.0*DN(0,1)*clhs15;
const double clhs85 =             clhs82*(clhs29 + clhs70);
const double clhs86 =             clhs82*(clhs48 + clhs77);
const double clhs87 =             1.0*clhs34;
const double clhs88 =             clhs15*clhs87;
const double clhs89 =             1.0*clhs33;
const double clhs90 =             clhs15*clhs89;
const double clhs91 =             N[1]*clhs10 - clhs13*clhs88 + clhs13*clhs90 + clhs31 + clhs32;
const double clhs92 =             pow(DN(1,0), 2);
const double clhs93 =             pow(N[1], 2);
const double clhs94 =             K_darcy*clhs93 + N[1]*clhs33 - clhs36*clhs88 + clhs36*clhs90 + clhs9*clhs93;
const double clhs95 =             DN(1,0)*clhs7;
const double clhs96 =             DN(1,1)*clhs95;
const double clhs97 =             -N[1] - clhs88 + clhs90;
const double clhs98 =             DN(1,0)*DN(2,0);
const double clhs99 =             clhs7*clhs98;
const double clhs100 =             N[2]*clhs34;
const double clhs101 =             N[2]*clhs35;
const double clhs102 =             N[1]*clhs52 + clhs100 + clhs101 - clhs55*clhs88 + clhs55*clhs90;
const double clhs103 =             DN(2,1)*clhs95;
const double clhs104 =             DN(1,0)*N[2];
const double clhs105 =             DN(2,0)*N[1];
const double clhs106 =             pow(DN(1,1), 2);
const double clhs107 =             DN(2,0)*clhs7;
const double clhs108 =             DN(1,1)*clhs107;
const double clhs109 =             DN(1,1)*DN(2,1);
const double clhs110 =             clhs109*clhs7;
const double clhs111 =             DN(1,1)*N[2];
const double clhs112 =             DN(2,1)*N[1];
const double clhs113 =             1.0*DN(1,0)*clhs15;
const double clhs114 =             1.0*DN(1,1)*clhs15;
const double clhs115 =             N[1] + clhs15*(1.0*clhs35 + clhs87 + clhs89);
const double clhs116 =             clhs82*(clhs109 + clhs98);
const double clhs117 =             1.0*clhs53;
const double clhs118 =             clhs117*clhs15;
const double clhs119 =             1.0*clhs52;
const double clhs120 =             clhs119*clhs15;
const double clhs121 =             N[2]*clhs10 - clhs118*clhs13 + clhs120*clhs13 + clhs50 + clhs51;
const double clhs122 =             N[2]*clhs33 + clhs100 + clhs101 - clhs118*clhs36 + clhs120*clhs36;
const double clhs123 =             pow(DN(2,0), 2);
const double clhs124 =             pow(N[2], 2);
const double clhs125 =             K_darcy*clhs124 + N[2]*clhs52 - clhs118*clhs55 + clhs120*clhs55 + clhs124*clhs9;
const double clhs126 =             DN(2,1)*clhs107;
const double clhs127 =             -N[2] - clhs118 + clhs120;
const double clhs128 =             pow(DN(2,1), 2);
const double clhs129 =             1.0*DN(2,0)*clhs15;
const double clhs130 =             1.0*DN(2,1)*clhs15;
const double clhs131 =             N[2] + clhs15*(clhs117 + clhs119 + 1.0*clhs54);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
            lhs(0,2)=DN(0,0)*clhs25;
            lhs(0,3)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs30 + clhs37;
            lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + clhs41;
            lhs(0,5)=DN(1,0)*clhs18 - clhs42 - clhs43*clhs44;
            lhs(0,6)=DN(0,0)*clhs45 + DN(0,1)*clhs47 + clhs49 + clhs56;
            lhs(0,7)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + clhs60;
            lhs(0,8)=DN(2,0)*clhs18 - clhs44*clhs62 - clhs61;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs63 + clhs24;
            lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs64 + clhs19 + clhs65*clhs7;
            lhs(1,2)=DN(0,1)*clhs25;
            lhs(1,3)=DN(0,0)*clhs28 + DN(0,1)*clhs66 + clhs68;
            lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs69 + clhs37 + clhs71;
            lhs(1,5)=DN(1,1)*clhs18 - clhs44*clhs73 - clhs72;
            lhs(1,6)=DN(0,0)*clhs47 + DN(0,1)*clhs74 + clhs75;
            lhs(1,7)=DN(0,0)*clhs59 + DN(0,1)*clhs76 + clhs56 + clhs78;
            lhs(1,8)=DN(2,1)*clhs18 - clhs44*clhs80 - clhs79;
            lhs(2,0)=DN(0,0)*clhs81;
            lhs(2,1)=DN(0,1)*clhs81;
            lhs(2,2)=clhs82*(clhs3 + clhs65);
            lhs(2,3)=clhs36*clhs83 + clhs43;
            lhs(2,4)=clhs36*clhs84 + clhs73;
            lhs(2,5)=clhs85;
            lhs(2,6)=clhs55*clhs83 + clhs62;
            lhs(2,7)=clhs55*clhs84 + clhs80;
            lhs(2,8)=clhs86;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs30 + clhs91;
            lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs68;
            lhs(3,2)=DN(0,0)*clhs90 - clhs42*clhs44 - clhs43;
            lhs(3,3)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs7*clhs92 + clhs94;
            lhs(3,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + clhs96;
            lhs(3,5)=DN(1,0)*clhs97;
            lhs(3,6)=DN(1,0)*clhs45 + DN(1,1)*clhs47 + clhs102 + clhs99;
            lhs(3,7)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + clhs103;
            lhs(3,8)=DN(2,0)*clhs90 - clhs104 - clhs105*clhs44;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs63 + clhs41;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs64 + clhs71 + clhs91;
            lhs(4,2)=DN(0,1)*clhs90 - clhs44*clhs72 - clhs73;
            lhs(4,3)=DN(1,0)*clhs28 + DN(1,1)*clhs66 + clhs96;
            lhs(4,4)=DN(1,0)*clhs40 + DN(1,1)*clhs69 + clhs106*clhs7 + clhs94;
            lhs(4,5)=DN(1,1)*clhs97;
            lhs(4,6)=DN(1,0)*clhs47 + DN(1,1)*clhs74 + clhs108;
            lhs(4,7)=DN(1,0)*clhs59 + DN(1,1)*clhs76 + clhs102 + clhs110;
            lhs(4,8)=DN(2,1)*clhs90 - clhs111 - clhs112*clhs44;
            lhs(5,0)=clhs113*clhs13 + clhs42;
            lhs(5,1)=clhs114*clhs13 + clhs72;
            lhs(5,2)=clhs85;
            lhs(5,3)=DN(1,0)*clhs115;
            lhs(5,4)=DN(1,1)*clhs115;
            lhs(5,5)=clhs82*(clhs106 + clhs92);
            lhs(5,6)=clhs105 + clhs113*clhs55;
            lhs(5,7)=clhs112 + clhs114*clhs55;
            lhs(5,8)=clhs116;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs121 + clhs49;
            lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs75;
            lhs(6,2)=DN(0,0)*clhs120 - clhs44*clhs61 - clhs62;
            lhs(6,3)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs122 + clhs99;
            lhs(6,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + clhs108;
            lhs(6,5)=DN(1,0)*clhs120 - clhs104*clhs44 - clhs105;
            lhs(6,6)=DN(2,0)*clhs45 + DN(2,1)*clhs47 + clhs123*clhs7 + clhs125;
            lhs(6,7)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + clhs126;
            lhs(6,8)=DN(2,0)*clhs127;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs63 + clhs60;
            lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs64 + clhs121 + clhs78;
            lhs(7,2)=DN(0,1)*clhs120 - clhs44*clhs79 - clhs80;
            lhs(7,3)=DN(2,0)*clhs28 + DN(2,1)*clhs66 + clhs103;
            lhs(7,4)=DN(2,0)*clhs40 + DN(2,1)*clhs69 + clhs110 + clhs122;
            lhs(7,5)=DN(1,1)*clhs120 - clhs111*clhs44 - clhs112;
            lhs(7,6)=DN(2,0)*clhs47 + DN(2,1)*clhs74 + clhs126;
            lhs(7,7)=DN(2,0)*clhs59 + DN(2,1)*clhs76 + clhs125 + clhs128*clhs7;
            lhs(7,8)=DN(2,1)*clhs127;
            lhs(8,0)=clhs129*clhs13 + clhs61;
            lhs(8,1)=clhs13*clhs130 + clhs79;
            lhs(8,2)=clhs86;
            lhs(8,3)=clhs104 + clhs129*clhs36;
            lhs(8,4)=clhs111 + clhs130*clhs36;
            lhs(8,5)=clhs116;
            lhs(8,6)=DN(2,0)*clhs131;
            lhs(8,7)=DN(2,1)*clhs131;
            lhs(8,8)=clhs82*(clhs123 + clhs128);


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
    const double K_darcy = rData.DarcyTerm;

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
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             bdf0*rho;
const double clhs13 =             rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs14 =             K_darcy*N[0];
const double clhs15 =             N[0]*clhs12;
const double clhs16 =             clhs13 + clhs14 + clhs15;
const double clhs17 =             1.0*clhs14;
const double clhs18 =             1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs19 =             clhs17*clhs18;
const double clhs20 =             1.0*clhs13;
const double clhs21 =             clhs18*clhs20;
const double clhs22 =             K_darcy*clhs11 + N[0]*clhs13 + clhs11*clhs12 - clhs16*clhs19 + clhs16*clhs21;
const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs24 =             C(1,3)*DN(0,1);
const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
const double clhs26 =             C(3,5)*DN(0,0);
const double clhs27 =             C(4,5)*DN(0,2);
const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
const double clhs29 =             DN(0,0)*clhs10;
const double clhs30 =             DN(0,1)*clhs29;
const double clhs31 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs32 =             C(3,4)*DN(0,1);
const double clhs33 =             C(2,3)*DN(0,2) + clhs26 + clhs32;
const double clhs34 =             C(2,5)*DN(0,2);
const double clhs35 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
const double clhs36 =             DN(0,2)*clhs29;
const double clhs37 =             -N[0] - clhs19 + clhs21;
const double clhs38 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs39 =             C(0,3)*DN(1,0);
const double clhs40 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs39;
const double clhs41 =             C(0,5)*DN(1,0);
const double clhs42 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs41;
const double clhs43 =             DN(0,0)*DN(1,0);
const double clhs44 =             clhs10*clhs43;
const double clhs45 =             N[1]*clhs14;
const double clhs46 =             N[1]*clhs15;
const double clhs47 =             rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs48 =             K_darcy*N[1];
const double clhs49 =             N[1]*clhs12;
const double clhs50 =             clhs47 + clhs48 + clhs49;
const double clhs51 =             N[0]*clhs47 - clhs19*clhs50 + clhs21*clhs50 + clhs45 + clhs46;
const double clhs52 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs39;
const double clhs53 =             C(1,3)*DN(1,1);
const double clhs54 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs53;
const double clhs55 =             C(3,5)*DN(1,0);
const double clhs56 =             C(4,5)*DN(1,2);
const double clhs57 =             C(1,5)*DN(1,1) + clhs55 + clhs56;
const double clhs58 =             DN(1,1)*clhs29;
const double clhs59 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs41;
const double clhs60 =             C(3,4)*DN(1,1);
const double clhs61 =             C(2,3)*DN(1,2) + clhs55 + clhs60;
const double clhs62 =             C(2,5)*DN(1,2);
const double clhs63 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs62;
const double clhs64 =             DN(1,2)*clhs29;
const double clhs65 =             DN(0,0)*N[1];
const double clhs66 =             DN(1,0)*N[0];
const double clhs67 =             1.0*K_darcy*clhs18;
const double clhs68 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs69 =             C(0,3)*DN(2,0);
const double clhs70 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs69;
const double clhs71 =             C(0,5)*DN(2,0);
const double clhs72 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs71;
const double clhs73 =             DN(0,0)*DN(2,0);
const double clhs74 =             clhs10*clhs73;
const double clhs75 =             N[2]*clhs14;
const double clhs76 =             N[2]*clhs15;
const double clhs77 =             rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs78 =             K_darcy*N[2];
const double clhs79 =             N[2]*clhs12;
const double clhs80 =             clhs77 + clhs78 + clhs79;
const double clhs81 =             N[0]*clhs77 - clhs19*clhs80 + clhs21*clhs80 + clhs75 + clhs76;
const double clhs82 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs69;
const double clhs83 =             C(1,3)*DN(2,1);
const double clhs84 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs83;
const double clhs85 =             C(3,5)*DN(2,0);
const double clhs86 =             C(4,5)*DN(2,2);
const double clhs87 =             C(1,5)*DN(2,1) + clhs85 + clhs86;
const double clhs88 =             DN(2,1)*clhs29;
const double clhs89 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs71;
const double clhs90 =             C(3,4)*DN(2,1);
const double clhs91 =             C(2,3)*DN(2,2) + clhs85 + clhs90;
const double clhs92 =             C(2,5)*DN(2,2);
const double clhs93 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs92;
const double clhs94 =             DN(2,2)*clhs29;
const double clhs95 =             DN(0,0)*N[2];
const double clhs96 =             DN(2,0)*N[0];
const double clhs97 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs98 =             C(0,3)*DN(3,0);
const double clhs99 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs98;
const double clhs100 =             C(0,5)*DN(3,0);
const double clhs101 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs100;
const double clhs102 =             DN(0,0)*DN(3,0);
const double clhs103 =             clhs10*clhs102;
const double clhs104 =             N[3]*clhs14;
const double clhs105 =             N[3]*clhs15;
const double clhs106 =             rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs107 =             K_darcy*N[3];
const double clhs108 =             N[3]*clhs12;
const double clhs109 =             clhs106 + clhs107 + clhs108;
const double clhs110 =             N[0]*clhs106 + clhs104 + clhs105 - clhs109*clhs19 + clhs109*clhs21;
const double clhs111 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs98;
const double clhs112 =             C(1,3)*DN(3,1);
const double clhs113 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs112;
const double clhs114 =             C(3,5)*DN(3,0);
const double clhs115 =             C(4,5)*DN(3,2);
const double clhs116 =             C(1,5)*DN(3,1) + clhs114 + clhs115;
const double clhs117 =             DN(3,1)*clhs29;
const double clhs118 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs100;
const double clhs119 =             C(3,4)*DN(3,1);
const double clhs120 =             C(2,3)*DN(3,2) + clhs114 + clhs119;
const double clhs121 =             C(2,5)*DN(3,2);
const double clhs122 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs121;
const double clhs123 =             DN(3,2)*clhs29;
const double clhs124 =             DN(0,0)*N[3];
const double clhs125 =             DN(3,0)*N[0];
const double clhs126 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs127 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs128 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs129 =             C(1,4)*DN(0,1);
const double clhs130 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs129;
const double clhs131 =             pow(DN(0,1), 2);
const double clhs132 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs129;
const double clhs133 =             C(2,4)*DN(0,2);
const double clhs134 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs133;
const double clhs135 =             DN(0,1)*clhs10;
const double clhs136 =             DN(0,2)*clhs135;
const double clhs137 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs53;
const double clhs138 =             C(0,4)*DN(1,0) + clhs56 + clhs60;
const double clhs139 =             DN(1,0)*clhs135;
const double clhs140 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs141 =             C(1,4)*DN(1,1);
const double clhs142 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs141;
const double clhs143 =             DN(0,1)*DN(1,1);
const double clhs144 =             clhs10*clhs143;
const double clhs145 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs141;
const double clhs146 =             C(2,4)*DN(1,2);
const double clhs147 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs146;
const double clhs148 =             DN(1,2)*clhs135;
const double clhs149 =             DN(0,1)*N[1];
const double clhs150 =             DN(1,1)*N[0];
const double clhs151 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs83;
const double clhs152 =             C(0,4)*DN(2,0) + clhs86 + clhs90;
const double clhs153 =             DN(2,0)*clhs135;
const double clhs154 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs155 =             C(1,4)*DN(2,1);
const double clhs156 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs155;
const double clhs157 =             DN(0,1)*DN(2,1);
const double clhs158 =             clhs10*clhs157;
const double clhs159 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs155;
const double clhs160 =             C(2,4)*DN(2,2);
const double clhs161 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs160;
const double clhs162 =             DN(2,2)*clhs135;
const double clhs163 =             DN(0,1)*N[2];
const double clhs164 =             DN(2,1)*N[0];
const double clhs165 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs112;
const double clhs166 =             C(0,4)*DN(3,0) + clhs115 + clhs119;
const double clhs167 =             DN(3,0)*clhs135;
const double clhs168 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs169 =             C(1,4)*DN(3,1);
const double clhs170 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs169;
const double clhs171 =             DN(0,1)*DN(3,1);
const double clhs172 =             clhs10*clhs171;
const double clhs173 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs169;
const double clhs174 =             C(2,4)*DN(3,2);
const double clhs175 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs174;
const double clhs176 =             DN(3,2)*clhs135;
const double clhs177 =             DN(0,1)*N[3];
const double clhs178 =             DN(3,1)*N[0];
const double clhs179 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs180 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs133;
const double clhs181 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs182 =             pow(DN(0,2), 2);
const double clhs183 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs62;
const double clhs184 =             DN(0,2)*clhs10;
const double clhs185 =             DN(1,0)*clhs184;
const double clhs186 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs146;
const double clhs187 =             DN(1,1)*clhs184;
const double clhs188 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs189 =             DN(0,2)*DN(1,2);
const double clhs190 =             clhs10*clhs189;
const double clhs191 =             DN(0,2)*N[1];
const double clhs192 =             DN(1,2)*N[0];
const double clhs193 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs92;
const double clhs194 =             DN(2,0)*clhs184;
const double clhs195 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs160;
const double clhs196 =             DN(2,1)*clhs184;
const double clhs197 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs198 =             DN(0,2)*DN(2,2);
const double clhs199 =             clhs10*clhs198;
const double clhs200 =             DN(0,2)*N[2];
const double clhs201 =             DN(2,2)*N[0];
const double clhs202 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs121;
const double clhs203 =             DN(3,0)*clhs184;
const double clhs204 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs174;
const double clhs205 =             DN(3,1)*clhs184;
const double clhs206 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs207 =             DN(0,2)*DN(3,2);
const double clhs208 =             clhs10*clhs207;
const double clhs209 =             DN(0,2)*N[3];
const double clhs210 =             DN(3,2)*N[0];
const double clhs211 =             N[0] + clhs18*(1.0*clhs15 + clhs17 + clhs20);
const double clhs212 =             1.0*clhs18;
const double clhs213 =             1.0*DN(0,0)*clhs18;
const double clhs214 =             1.0*DN(0,1)*clhs18;
const double clhs215 =             1.0*DN(0,2)*clhs18;
const double clhs216 =             clhs212*(clhs143 + clhs189 + clhs43);
const double clhs217 =             clhs212*(clhs157 + clhs198 + clhs73);
const double clhs218 =             clhs212*(clhs102 + clhs171 + clhs207);
const double clhs219 =             1.0*clhs48;
const double clhs220 =             clhs18*clhs219;
const double clhs221 =             1.0*clhs47;
const double clhs222 =             clhs18*clhs221;
const double clhs223 =             N[1]*clhs13 - clhs16*clhs220 + clhs16*clhs222 + clhs45 + clhs46;
const double clhs224 =             pow(DN(1,0), 2);
const double clhs225 =             pow(N[1], 2);
const double clhs226 =             K_darcy*clhs225 + N[1]*clhs47 + clhs12*clhs225 - clhs220*clhs50 + clhs222*clhs50;
const double clhs227 =             DN(1,0)*clhs10;
const double clhs228 =             DN(1,1)*clhs227;
const double clhs229 =             DN(1,2)*clhs227;
const double clhs230 =             -N[1] - clhs220 + clhs222;
const double clhs231 =             DN(1,0)*DN(2,0);
const double clhs232 =             clhs10*clhs231;
const double clhs233 =             N[2]*clhs48;
const double clhs234 =             N[2]*clhs49;
const double clhs235 =             N[1]*clhs77 - clhs220*clhs80 + clhs222*clhs80 + clhs233 + clhs234;
const double clhs236 =             DN(2,1)*clhs227;
const double clhs237 =             DN(2,2)*clhs227;
const double clhs238 =             DN(1,0)*N[2];
const double clhs239 =             DN(2,0)*N[1];
const double clhs240 =             DN(1,0)*DN(3,0);
const double clhs241 =             clhs10*clhs240;
const double clhs242 =             N[3]*clhs48;
const double clhs243 =             N[3]*clhs49;
const double clhs244 =             N[1]*clhs106 - clhs109*clhs220 + clhs109*clhs222 + clhs242 + clhs243;
const double clhs245 =             DN(3,1)*clhs227;
const double clhs246 =             DN(3,2)*clhs227;
const double clhs247 =             DN(1,0)*N[3];
const double clhs248 =             DN(3,0)*N[1];
const double clhs249 =             pow(DN(1,1), 2);
const double clhs250 =             DN(1,1)*clhs10;
const double clhs251 =             DN(1,2)*clhs250;
const double clhs252 =             DN(2,0)*clhs250;
const double clhs253 =             DN(1,1)*DN(2,1);
const double clhs254 =             clhs10*clhs253;
const double clhs255 =             DN(2,2)*clhs250;
const double clhs256 =             DN(1,1)*N[2];
const double clhs257 =             DN(2,1)*N[1];
const double clhs258 =             DN(3,0)*clhs250;
const double clhs259 =             DN(1,1)*DN(3,1);
const double clhs260 =             clhs10*clhs259;
const double clhs261 =             DN(3,2)*clhs250;
const double clhs262 =             DN(1,1)*N[3];
const double clhs263 =             DN(3,1)*N[1];
const double clhs264 =             pow(DN(1,2), 2);
const double clhs265 =             DN(1,2)*clhs10;
const double clhs266 =             DN(2,0)*clhs265;
const double clhs267 =             DN(2,1)*clhs265;
const double clhs268 =             DN(1,2)*DN(2,2);
const double clhs269 =             clhs10*clhs268;
const double clhs270 =             DN(1,2)*N[2];
const double clhs271 =             DN(2,2)*N[1];
const double clhs272 =             DN(3,0)*clhs265;
const double clhs273 =             DN(3,1)*clhs265;
const double clhs274 =             DN(1,2)*DN(3,2);
const double clhs275 =             clhs10*clhs274;
const double clhs276 =             DN(1,2)*N[3];
const double clhs277 =             DN(3,2)*N[1];
const double clhs278 =             1.0*DN(1,0)*clhs18;
const double clhs279 =             1.0*DN(1,1)*clhs18;
const double clhs280 =             1.0*DN(1,2)*clhs18;
const double clhs281 =             N[1] + clhs18*(clhs219 + clhs221 + 1.0*clhs49);
const double clhs282 =             clhs212*(clhs231 + clhs253 + clhs268);
const double clhs283 =             clhs212*(clhs240 + clhs259 + clhs274);
const double clhs284 =             1.0*clhs78;
const double clhs285 =             clhs18*clhs284;
const double clhs286 =             1.0*clhs77;
const double clhs287 =             clhs18*clhs286;
const double clhs288 =             N[2]*clhs13 - clhs16*clhs285 + clhs16*clhs287 + clhs75 + clhs76;
const double clhs289 =             N[2]*clhs47 + clhs233 + clhs234 - clhs285*clhs50 + clhs287*clhs50;
const double clhs290 =             pow(DN(2,0), 2);
const double clhs291 =             pow(N[2], 2);
const double clhs292 =             K_darcy*clhs291 + N[2]*clhs77 + clhs12*clhs291 - clhs285*clhs80 + clhs287*clhs80;
const double clhs293 =             DN(2,0)*clhs10;
const double clhs294 =             DN(2,1)*clhs293;
const double clhs295 =             DN(2,2)*clhs293;
const double clhs296 =             -N[2] - clhs285 + clhs287;
const double clhs297 =             DN(2,0)*DN(3,0);
const double clhs298 =             clhs10*clhs297;
const double clhs299 =             N[3]*clhs78;
const double clhs300 =             N[3]*clhs79;
const double clhs301 =             N[2]*clhs106 - clhs109*clhs285 + clhs109*clhs287 + clhs299 + clhs300;
const double clhs302 =             DN(3,1)*clhs293;
const double clhs303 =             DN(3,2)*clhs293;
const double clhs304 =             DN(2,0)*N[3];
const double clhs305 =             DN(3,0)*N[2];
const double clhs306 =             pow(DN(2,1), 2);
const double clhs307 =             DN(2,1)*clhs10;
const double clhs308 =             DN(2,2)*clhs307;
const double clhs309 =             DN(3,0)*clhs307;
const double clhs310 =             DN(2,1)*DN(3,1);
const double clhs311 =             clhs10*clhs310;
const double clhs312 =             DN(3,2)*clhs307;
const double clhs313 =             DN(2,1)*N[3];
const double clhs314 =             DN(3,1)*N[2];
const double clhs315 =             pow(DN(2,2), 2);
const double clhs316 =             DN(2,2)*clhs10;
const double clhs317 =             DN(3,0)*clhs316;
const double clhs318 =             DN(3,1)*clhs316;
const double clhs319 =             DN(2,2)*DN(3,2);
const double clhs320 =             clhs10*clhs319;
const double clhs321 =             DN(2,2)*N[3];
const double clhs322 =             DN(3,2)*N[2];
const double clhs323 =             1.0*DN(2,0)*clhs18;
const double clhs324 =             1.0*DN(2,1)*clhs18;
const double clhs325 =             1.0*DN(2,2)*clhs18;
const double clhs326 =             N[2] + clhs18*(clhs284 + clhs286 + 1.0*clhs79);
const double clhs327 =             clhs212*(clhs297 + clhs310 + clhs319);
const double clhs328 =             1.0*clhs107;
const double clhs329 =             clhs18*clhs328;
const double clhs330 =             1.0*clhs106;
const double clhs331 =             clhs18*clhs330;
const double clhs332 =             N[3]*clhs13 + clhs104 + clhs105 - clhs16*clhs329 + clhs16*clhs331;
const double clhs333 =             N[3]*clhs47 + clhs242 + clhs243 - clhs329*clhs50 + clhs331*clhs50;
const double clhs334 =             N[3]*clhs77 + clhs299 + clhs300 - clhs329*clhs80 + clhs331*clhs80;
const double clhs335 =             pow(DN(3,0), 2);
const double clhs336 =             pow(N[3], 2);
const double clhs337 =             K_darcy*clhs336 + N[3]*clhs106 - clhs109*clhs329 + clhs109*clhs331 + clhs12*clhs336;
const double clhs338 =             DN(3,0)*clhs10;
const double clhs339 =             DN(3,1)*clhs338;
const double clhs340 =             DN(3,2)*clhs338;
const double clhs341 =             -N[3] - clhs329 + clhs331;
const double clhs342 =             pow(DN(3,1), 2);
const double clhs343 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs344 =             pow(DN(3,2), 2);
const double clhs345 =             1.0*DN(3,0)*clhs18;
const double clhs346 =             1.0*DN(3,1)*clhs18;
const double clhs347 =             1.0*DN(3,2)*clhs18;
const double clhs348 =             N[3] + clhs18*(1.0*clhs108 + clhs328 + clhs330);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
            lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
            lhs(0,3)=DN(0,0)*clhs37;
            lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + DN(0,2)*clhs42 + clhs44 + clhs51;
            lhs(0,5)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + DN(0,2)*clhs57 + clhs58;
            lhs(0,6)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + DN(0,2)*clhs63 + clhs64;
            lhs(0,7)=DN(1,0)*clhs21 - clhs65 - clhs66*clhs67;
            lhs(0,8)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs74 + clhs81;
            lhs(0,9)=DN(0,0)*clhs82 + DN(0,1)*clhs84 + DN(0,2)*clhs87 + clhs88;
            lhs(0,10)=DN(0,0)*clhs89 + DN(0,1)*clhs91 + DN(0,2)*clhs93 + clhs94;
            lhs(0,11)=DN(2,0)*clhs21 - clhs67*clhs96 - clhs95;
            lhs(0,12)=DN(0,0)*clhs97 + DN(0,1)*clhs99 + DN(0,2)*clhs101 + clhs103 + clhs110;
            lhs(0,13)=DN(0,0)*clhs111 + DN(0,1)*clhs113 + DN(0,2)*clhs116 + clhs117;
            lhs(0,14)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs123;
            lhs(0,15)=DN(3,0)*clhs21 - clhs124 - clhs125*clhs67;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs126 + DN(0,2)*clhs127 + clhs30;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs128 + DN(0,2)*clhs130 + clhs10*clhs131 + clhs22;
            lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs136;
            lhs(1,3)=DN(0,1)*clhs37;
            lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
            lhs(1,5)=DN(0,0)*clhs54 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs144 + clhs51;
            lhs(1,6)=DN(0,0)*clhs61 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148;
            lhs(1,7)=DN(1,1)*clhs21 - clhs149 - clhs150*clhs67;
            lhs(1,8)=DN(0,0)*clhs70 + DN(0,1)*clhs151 + DN(0,2)*clhs152 + clhs153;
            lhs(1,9)=DN(0,0)*clhs84 + DN(0,1)*clhs154 + DN(0,2)*clhs156 + clhs158 + clhs81;
            lhs(1,10)=DN(0,0)*clhs91 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162;
            lhs(1,11)=DN(2,1)*clhs21 - clhs163 - clhs164*clhs67;
            lhs(1,12)=DN(0,0)*clhs99 + DN(0,1)*clhs165 + DN(0,2)*clhs166 + clhs167;
            lhs(1,13)=DN(0,0)*clhs113 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs110 + clhs172;
            lhs(1,14)=DN(0,0)*clhs120 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176;
            lhs(1,15)=DN(3,1)*clhs21 - clhs177 - clhs178*clhs67;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs179 + clhs36;
            lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs130 + DN(0,2)*clhs180 + clhs136;
            lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs134 + DN(0,2)*clhs181 + clhs10*clhs182 + clhs22;
            lhs(2,3)=DN(0,2)*clhs37;
            lhs(2,4)=DN(0,0)*clhs42 + DN(0,1)*clhs138 + DN(0,2)*clhs183 + clhs185;
            lhs(2,5)=DN(0,0)*clhs57 + DN(0,1)*clhs142 + DN(0,2)*clhs186 + clhs187;
            lhs(2,6)=DN(0,0)*clhs63 + DN(0,1)*clhs147 + DN(0,2)*clhs188 + clhs190 + clhs51;
            lhs(2,7)=DN(1,2)*clhs21 - clhs191 - clhs192*clhs67;
            lhs(2,8)=DN(0,0)*clhs72 + DN(0,1)*clhs152 + DN(0,2)*clhs193 + clhs194;
            lhs(2,9)=DN(0,0)*clhs87 + DN(0,1)*clhs156 + DN(0,2)*clhs195 + clhs196;
            lhs(2,10)=DN(0,0)*clhs93 + DN(0,1)*clhs161 + DN(0,2)*clhs197 + clhs199 + clhs81;
            lhs(2,11)=DN(2,2)*clhs21 - clhs200 - clhs201*clhs67;
            lhs(2,12)=DN(0,0)*clhs101 + DN(0,1)*clhs166 + DN(0,2)*clhs202 + clhs203;
            lhs(2,13)=DN(0,0)*clhs116 + DN(0,1)*clhs170 + DN(0,2)*clhs204 + clhs205;
            lhs(2,14)=DN(0,0)*clhs122 + DN(0,1)*clhs175 + DN(0,2)*clhs206 + clhs110 + clhs208;
            lhs(2,15)=DN(3,2)*clhs21 - clhs209 - clhs210*clhs67;
            lhs(3,0)=DN(0,0)*clhs211;
            lhs(3,1)=DN(0,1)*clhs211;
            lhs(3,2)=DN(0,2)*clhs211;
            lhs(3,3)=clhs212*(clhs131 + clhs182 + clhs5);
            lhs(3,4)=clhs213*clhs50 + clhs66;
            lhs(3,5)=clhs150 + clhs214*clhs50;
            lhs(3,6)=clhs192 + clhs215*clhs50;
            lhs(3,7)=clhs216;
            lhs(3,8)=clhs213*clhs80 + clhs96;
            lhs(3,9)=clhs164 + clhs214*clhs80;
            lhs(3,10)=clhs201 + clhs215*clhs80;
            lhs(3,11)=clhs217;
            lhs(3,12)=clhs109*clhs213 + clhs125;
            lhs(3,13)=clhs109*clhs214 + clhs178;
            lhs(3,14)=clhs109*clhs215 + clhs210;
            lhs(3,15)=clhs218;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs223 + clhs44;
            lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs139;
            lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs185;
            lhs(4,3)=DN(0,0)*clhs222 - clhs65*clhs67 - clhs66;
            lhs(4,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + DN(1,2)*clhs42 + clhs10*clhs224 + clhs226;
            lhs(4,5)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + DN(1,2)*clhs57 + clhs228;
            lhs(4,6)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + DN(1,2)*clhs63 + clhs229;
            lhs(4,7)=DN(1,0)*clhs230;
            lhs(4,8)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs232 + clhs235;
            lhs(4,9)=DN(1,0)*clhs82 + DN(1,1)*clhs84 + DN(1,2)*clhs87 + clhs236;
            lhs(4,10)=DN(1,0)*clhs89 + DN(1,1)*clhs91 + DN(1,2)*clhs93 + clhs237;
            lhs(4,11)=DN(2,0)*clhs222 - clhs238 - clhs239*clhs67;
            lhs(4,12)=DN(1,0)*clhs97 + DN(1,1)*clhs99 + DN(1,2)*clhs101 + clhs241 + clhs244;
            lhs(4,13)=DN(1,0)*clhs111 + DN(1,1)*clhs113 + DN(1,2)*clhs116 + clhs245;
            lhs(4,14)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs246;
            lhs(4,15)=DN(3,0)*clhs222 - clhs247 - clhs248*clhs67;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs126 + DN(1,2)*clhs127 + clhs58;
            lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs128 + DN(1,2)*clhs130 + clhs144 + clhs223;
            lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs187;
            lhs(5,3)=DN(0,1)*clhs222 - clhs149*clhs67 - clhs150;
            lhs(5,4)=DN(1,0)*clhs40 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs228;
            lhs(5,5)=DN(1,0)*clhs54 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs10*clhs249 + clhs226;
            lhs(5,6)=DN(1,0)*clhs61 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs251;
            lhs(5,7)=DN(1,1)*clhs230;
            lhs(5,8)=DN(1,0)*clhs70 + DN(1,1)*clhs151 + DN(1,2)*clhs152 + clhs252;
            lhs(5,9)=DN(1,0)*clhs84 + DN(1,1)*clhs154 + DN(1,2)*clhs156 + clhs235 + clhs254;
            lhs(5,10)=DN(1,0)*clhs91 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs255;
            lhs(5,11)=DN(2,1)*clhs222 - clhs256 - clhs257*clhs67;
            lhs(5,12)=DN(1,0)*clhs99 + DN(1,1)*clhs165 + DN(1,2)*clhs166 + clhs258;
            lhs(5,13)=DN(1,0)*clhs113 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs244 + clhs260;
            lhs(5,14)=DN(1,0)*clhs120 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs261;
            lhs(5,15)=DN(3,1)*clhs222 - clhs262 - clhs263*clhs67;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs179 + clhs64;
            lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs130 + DN(1,2)*clhs180 + clhs148;
            lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs134 + DN(1,2)*clhs181 + clhs190 + clhs223;
            lhs(6,3)=DN(0,2)*clhs222 - clhs191*clhs67 - clhs192;
            lhs(6,4)=DN(1,0)*clhs42 + DN(1,1)*clhs138 + DN(1,2)*clhs183 + clhs229;
            lhs(6,5)=DN(1,0)*clhs57 + DN(1,1)*clhs142 + DN(1,2)*clhs186 + clhs251;
            lhs(6,6)=DN(1,0)*clhs63 + DN(1,1)*clhs147 + DN(1,2)*clhs188 + clhs10*clhs264 + clhs226;
            lhs(6,7)=DN(1,2)*clhs230;
            lhs(6,8)=DN(1,0)*clhs72 + DN(1,1)*clhs152 + DN(1,2)*clhs193 + clhs266;
            lhs(6,9)=DN(1,0)*clhs87 + DN(1,1)*clhs156 + DN(1,2)*clhs195 + clhs267;
            lhs(6,10)=DN(1,0)*clhs93 + DN(1,1)*clhs161 + DN(1,2)*clhs197 + clhs235 + clhs269;
            lhs(6,11)=DN(2,2)*clhs222 - clhs270 - clhs271*clhs67;
            lhs(6,12)=DN(1,0)*clhs101 + DN(1,1)*clhs166 + DN(1,2)*clhs202 + clhs272;
            lhs(6,13)=DN(1,0)*clhs116 + DN(1,1)*clhs170 + DN(1,2)*clhs204 + clhs273;
            lhs(6,14)=DN(1,0)*clhs122 + DN(1,1)*clhs175 + DN(1,2)*clhs206 + clhs244 + clhs275;
            lhs(6,15)=DN(3,2)*clhs222 - clhs276 - clhs277*clhs67;
            lhs(7,0)=clhs16*clhs278 + clhs65;
            lhs(7,1)=clhs149 + clhs16*clhs279;
            lhs(7,2)=clhs16*clhs280 + clhs191;
            lhs(7,3)=clhs216;
            lhs(7,4)=DN(1,0)*clhs281;
            lhs(7,5)=DN(1,1)*clhs281;
            lhs(7,6)=DN(1,2)*clhs281;
            lhs(7,7)=clhs212*(clhs224 + clhs249 + clhs264);
            lhs(7,8)=clhs239 + clhs278*clhs80;
            lhs(7,9)=clhs257 + clhs279*clhs80;
            lhs(7,10)=clhs271 + clhs280*clhs80;
            lhs(7,11)=clhs282;
            lhs(7,12)=clhs109*clhs278 + clhs248;
            lhs(7,13)=clhs109*clhs279 + clhs263;
            lhs(7,14)=clhs109*clhs280 + clhs277;
            lhs(7,15)=clhs283;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs288 + clhs74;
            lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs153;
            lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs194;
            lhs(8,3)=DN(0,0)*clhs287 - clhs67*clhs95 - clhs96;
            lhs(8,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + DN(2,2)*clhs42 + clhs232 + clhs289;
            lhs(8,5)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + DN(2,2)*clhs57 + clhs252;
            lhs(8,6)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + DN(2,2)*clhs63 + clhs266;
            lhs(8,7)=DN(1,0)*clhs287 - clhs238*clhs67 - clhs239;
            lhs(8,8)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs10*clhs290 + clhs292;
            lhs(8,9)=DN(2,0)*clhs82 + DN(2,1)*clhs84 + DN(2,2)*clhs87 + clhs294;
            lhs(8,10)=DN(2,0)*clhs89 + DN(2,1)*clhs91 + DN(2,2)*clhs93 + clhs295;
            lhs(8,11)=DN(2,0)*clhs296;
            lhs(8,12)=DN(2,0)*clhs97 + DN(2,1)*clhs99 + DN(2,2)*clhs101 + clhs298 + clhs301;
            lhs(8,13)=DN(2,0)*clhs111 + DN(2,1)*clhs113 + DN(2,2)*clhs116 + clhs302;
            lhs(8,14)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs303;
            lhs(8,15)=DN(3,0)*clhs287 - clhs304 - clhs305*clhs67;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs126 + DN(2,2)*clhs127 + clhs88;
            lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs128 + DN(2,2)*clhs130 + clhs158 + clhs288;
            lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs196;
            lhs(9,3)=DN(0,1)*clhs287 - clhs163*clhs67 - clhs164;
            lhs(9,4)=DN(2,0)*clhs40 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs236;
            lhs(9,5)=DN(2,0)*clhs54 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs254 + clhs289;
            lhs(9,6)=DN(2,0)*clhs61 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs267;
            lhs(9,7)=DN(1,1)*clhs287 - clhs256*clhs67 - clhs257;
            lhs(9,8)=DN(2,0)*clhs70 + DN(2,1)*clhs151 + DN(2,2)*clhs152 + clhs294;
            lhs(9,9)=DN(2,0)*clhs84 + DN(2,1)*clhs154 + DN(2,2)*clhs156 + clhs10*clhs306 + clhs292;
            lhs(9,10)=DN(2,0)*clhs91 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs308;
            lhs(9,11)=DN(2,1)*clhs296;
            lhs(9,12)=DN(2,0)*clhs99 + DN(2,1)*clhs165 + DN(2,2)*clhs166 + clhs309;
            lhs(9,13)=DN(2,0)*clhs113 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs301 + clhs311;
            lhs(9,14)=DN(2,0)*clhs120 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs312;
            lhs(9,15)=DN(3,1)*clhs287 - clhs313 - clhs314*clhs67;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs179 + clhs94;
            lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs130 + DN(2,2)*clhs180 + clhs162;
            lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs134 + DN(2,2)*clhs181 + clhs199 + clhs288;
            lhs(10,3)=DN(0,2)*clhs287 - clhs200*clhs67 - clhs201;
            lhs(10,4)=DN(2,0)*clhs42 + DN(2,1)*clhs138 + DN(2,2)*clhs183 + clhs237;
            lhs(10,5)=DN(2,0)*clhs57 + DN(2,1)*clhs142 + DN(2,2)*clhs186 + clhs255;
            lhs(10,6)=DN(2,0)*clhs63 + DN(2,1)*clhs147 + DN(2,2)*clhs188 + clhs269 + clhs289;
            lhs(10,7)=DN(1,2)*clhs287 - clhs270*clhs67 - clhs271;
            lhs(10,8)=DN(2,0)*clhs72 + DN(2,1)*clhs152 + DN(2,2)*clhs193 + clhs295;
            lhs(10,9)=DN(2,0)*clhs87 + DN(2,1)*clhs156 + DN(2,2)*clhs195 + clhs308;
            lhs(10,10)=DN(2,0)*clhs93 + DN(2,1)*clhs161 + DN(2,2)*clhs197 + clhs10*clhs315 + clhs292;
            lhs(10,11)=DN(2,2)*clhs296;
            lhs(10,12)=DN(2,0)*clhs101 + DN(2,1)*clhs166 + DN(2,2)*clhs202 + clhs317;
            lhs(10,13)=DN(2,0)*clhs116 + DN(2,1)*clhs170 + DN(2,2)*clhs204 + clhs318;
            lhs(10,14)=DN(2,0)*clhs122 + DN(2,1)*clhs175 + DN(2,2)*clhs206 + clhs301 + clhs320;
            lhs(10,15)=DN(3,2)*clhs287 - clhs321 - clhs322*clhs67;
            lhs(11,0)=clhs16*clhs323 + clhs95;
            lhs(11,1)=clhs16*clhs324 + clhs163;
            lhs(11,2)=clhs16*clhs325 + clhs200;
            lhs(11,3)=clhs217;
            lhs(11,4)=clhs238 + clhs323*clhs50;
            lhs(11,5)=clhs256 + clhs324*clhs50;
            lhs(11,6)=clhs270 + clhs325*clhs50;
            lhs(11,7)=clhs282;
            lhs(11,8)=DN(2,0)*clhs326;
            lhs(11,9)=DN(2,1)*clhs326;
            lhs(11,10)=DN(2,2)*clhs326;
            lhs(11,11)=clhs212*(clhs290 + clhs306 + clhs315);
            lhs(11,12)=clhs109*clhs323 + clhs305;
            lhs(11,13)=clhs109*clhs324 + clhs314;
            lhs(11,14)=clhs109*clhs325 + clhs322;
            lhs(11,15)=clhs327;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs103 + clhs332;
            lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs167;
            lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs203;
            lhs(12,3)=DN(0,0)*clhs331 - clhs124*clhs67 - clhs125;
            lhs(12,4)=DN(3,0)*clhs38 + DN(3,1)*clhs40 + DN(3,2)*clhs42 + clhs241 + clhs333;
            lhs(12,5)=DN(3,0)*clhs52 + DN(3,1)*clhs54 + DN(3,2)*clhs57 + clhs258;
            lhs(12,6)=DN(3,0)*clhs59 + DN(3,1)*clhs61 + DN(3,2)*clhs63 + clhs272;
            lhs(12,7)=DN(1,0)*clhs331 - clhs247*clhs67 - clhs248;
            lhs(12,8)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs298 + clhs334;
            lhs(12,9)=DN(3,0)*clhs82 + DN(3,1)*clhs84 + DN(3,2)*clhs87 + clhs309;
            lhs(12,10)=DN(3,0)*clhs89 + DN(3,1)*clhs91 + DN(3,2)*clhs93 + clhs317;
            lhs(12,11)=DN(2,0)*clhs331 - clhs304*clhs67 - clhs305;
            lhs(12,12)=DN(3,0)*clhs97 + DN(3,1)*clhs99 + DN(3,2)*clhs101 + clhs10*clhs335 + clhs337;
            lhs(12,13)=DN(3,0)*clhs111 + DN(3,1)*clhs113 + DN(3,2)*clhs116 + clhs339;
            lhs(12,14)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs340;
            lhs(12,15)=DN(3,0)*clhs341;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs126 + DN(3,2)*clhs127 + clhs117;
            lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs128 + DN(3,2)*clhs130 + clhs172 + clhs332;
            lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs205;
            lhs(13,3)=DN(0,1)*clhs331 - clhs177*clhs67 - clhs178;
            lhs(13,4)=DN(3,0)*clhs40 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs245;
            lhs(13,5)=DN(3,0)*clhs54 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs260 + clhs333;
            lhs(13,6)=DN(3,0)*clhs61 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs273;
            lhs(13,7)=DN(1,1)*clhs331 - clhs262*clhs67 - clhs263;
            lhs(13,8)=DN(3,0)*clhs70 + DN(3,1)*clhs151 + DN(3,2)*clhs152 + clhs302;
            lhs(13,9)=DN(3,0)*clhs84 + DN(3,1)*clhs154 + DN(3,2)*clhs156 + clhs311 + clhs334;
            lhs(13,10)=DN(3,0)*clhs91 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs318;
            lhs(13,11)=DN(2,1)*clhs331 - clhs313*clhs67 - clhs314;
            lhs(13,12)=DN(3,0)*clhs99 + DN(3,1)*clhs165 + DN(3,2)*clhs166 + clhs339;
            lhs(13,13)=DN(3,0)*clhs113 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs10*clhs342 + clhs337;
            lhs(13,14)=DN(3,0)*clhs120 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs343;
            lhs(13,15)=DN(3,1)*clhs341;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs179 + clhs123;
            lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs130 + DN(3,2)*clhs180 + clhs176;
            lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs134 + DN(3,2)*clhs181 + clhs208 + clhs332;
            lhs(14,3)=DN(0,2)*clhs331 - clhs209*clhs67 - clhs210;
            lhs(14,4)=DN(3,0)*clhs42 + DN(3,1)*clhs138 + DN(3,2)*clhs183 + clhs246;
            lhs(14,5)=DN(3,0)*clhs57 + DN(3,1)*clhs142 + DN(3,2)*clhs186 + clhs261;
            lhs(14,6)=DN(3,0)*clhs63 + DN(3,1)*clhs147 + DN(3,2)*clhs188 + clhs275 + clhs333;
            lhs(14,7)=DN(1,2)*clhs331 - clhs276*clhs67 - clhs277;
            lhs(14,8)=DN(3,0)*clhs72 + DN(3,1)*clhs152 + DN(3,2)*clhs193 + clhs303;
            lhs(14,9)=DN(3,0)*clhs87 + DN(3,1)*clhs156 + DN(3,2)*clhs195 + clhs312;
            lhs(14,10)=DN(3,0)*clhs93 + DN(3,1)*clhs161 + DN(3,2)*clhs197 + clhs320 + clhs334;
            lhs(14,11)=DN(2,2)*clhs331 - clhs321*clhs67 - clhs322;
            lhs(14,12)=DN(3,0)*clhs101 + DN(3,1)*clhs166 + DN(3,2)*clhs202 + clhs340;
            lhs(14,13)=DN(3,0)*clhs116 + DN(3,1)*clhs170 + DN(3,2)*clhs204 + clhs343;
            lhs(14,14)=DN(3,0)*clhs122 + DN(3,1)*clhs175 + DN(3,2)*clhs206 + clhs10*clhs344 + clhs337;
            lhs(14,15)=DN(3,2)*clhs341;
            lhs(15,0)=clhs124 + clhs16*clhs345;
            lhs(15,1)=clhs16*clhs346 + clhs177;
            lhs(15,2)=clhs16*clhs347 + clhs209;
            lhs(15,3)=clhs218;
            lhs(15,4)=clhs247 + clhs345*clhs50;
            lhs(15,5)=clhs262 + clhs346*clhs50;
            lhs(15,6)=clhs276 + clhs347*clhs50;
            lhs(15,7)=clhs283;
            lhs(15,8)=clhs304 + clhs345*clhs80;
            lhs(15,9)=clhs313 + clhs346*clhs80;
            lhs(15,10)=clhs321 + clhs347*clhs80;
            lhs(15,11)=clhs327;
            lhs(15,12)=DN(3,0)*clhs348;
            lhs(15,13)=DN(3,1)*clhs348;
            lhs(15,14)=DN(3,2)*clhs348;
            lhs(15,15)=clhs212*(clhs335 + clhs342 + clhs344);


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
    const double K_darcy = rData.DarcyTerm;

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
const double crhs2 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs7 =             rho*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs8 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs9 =             crhs4 + crhs8;
const double crhs10 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
const double crhs11 =             crhs9*(crhs10*h/stab_c1 + mu);
const double crhs12 =             K_darcy*N[0];
const double crhs13 =             1.0/(K_darcy + crhs10/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs14 =             1.0*crhs13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs3 + crhs7);
const double crhs15 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
const double crhs16 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs17 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs18 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs19 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs8);
const double crhs20 =             1.0*crhs13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs16 + crhs17 + crhs18 + crhs19);
const double crhs21 =             K_darcy*N[1];
const double crhs22 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
const double crhs23 =             K_darcy*N[2];
const double crhs24 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs7 + crhs12*crhs14 - crhs14*crhs15;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs11 - DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 - N[0]*crhs18 - N[0]*crhs19 + crhs12*crhs20 - crhs15*crhs20;
            rhs[2]=-DN(0,0)*crhs14 - DN(0,1)*crhs20 - N[0]*crhs9;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs7 + crhs14*crhs21 - crhs14*crhs22;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs11 - DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 - N[1]*crhs18 - N[1]*crhs19 + crhs20*crhs21 - crhs20*crhs22;
            rhs[5]=-DN(1,0)*crhs14 - DN(1,1)*crhs20 - N[1]*crhs9;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs7 + crhs14*crhs23 - crhs14*crhs24;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs11 - DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 - N[2]*crhs18 - N[2]*crhs19 + crhs20*crhs23 - crhs20*crhs24;
            rhs[8]=-DN(2,0)*crhs14 - DN(2,1)*crhs20 - N[2]*crhs9;


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
    const double K_darcy = rData.DarcyTerm;

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
const double crhs2 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 =             DN(0,0)*v(0,0);
const double crhs5 =             DN(1,0)*v(1,0);
const double crhs6 =             DN(2,0)*v(2,0);
const double crhs7 =             DN(3,0)*v(3,0);
const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs11 =             rho*(crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*(crhs4 + crhs5 + crhs6 + crhs7) + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
const double crhs12 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs13 =             DN(0,1)*v(0,1);
const double crhs14 =             DN(1,1)*v(1,1);
const double crhs15 =             DN(2,1)*v(2,1);
const double crhs16 =             DN(3,1)*v(3,1);
const double crhs17 =             crhs12 + crhs13 + crhs14 + crhs15 + crhs16 + crhs4 + crhs5 + crhs6 + crhs7;
const double crhs18 =             rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs19 =             crhs17*(crhs18*h/stab_c1 + mu);
const double crhs20 =             K_darcy*N[0];
const double crhs21 =             1.0/(K_darcy + crhs18/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 =             1.0*crhs21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs11 + crhs2 + crhs3);
const double crhs23 =             rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10);
const double crhs24 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs25 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs26 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs27 =             rho*(crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs9*(crhs13 + crhs14 + crhs15 + crhs16));
const double crhs28 =             1.0*crhs21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs24 + crhs25 + crhs26 + crhs27);
const double crhs29 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs30 =             K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs31 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs32 =             rho*(crhs10*crhs12 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs33 =             1.0*crhs21*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs29 + crhs30 + crhs31 + crhs32);
const double crhs34 =             K_darcy*N[1];
const double crhs35 =             rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10);
const double crhs36 =             K_darcy*N[2];
const double crhs37 =             rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10);
const double crhs38 =             K_darcy*N[3];
const double crhs39 =             rho*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs19 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs11 - N[0]*crhs2 - N[0]*crhs3 + crhs20*crhs22 - crhs22*crhs23;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs19 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs27 + crhs20*crhs28 - crhs23*crhs28;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs19 - DN(0,2)*stress[2] + N[0]*crhs29 - N[0]*crhs30 - N[0]*crhs31 - N[0]*crhs32 + crhs20*crhs33 - crhs23*crhs33;
            rhs[3]=-DN(0,0)*crhs22 - DN(0,1)*crhs28 - DN(0,2)*crhs33 - N[0]*crhs17;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs19 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs11 - N[1]*crhs2 - N[1]*crhs3 + crhs22*crhs34 - crhs22*crhs35;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs19 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs27 + crhs28*crhs34 - crhs28*crhs35;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs19 - DN(1,2)*stress[2] + N[1]*crhs29 - N[1]*crhs30 - N[1]*crhs31 - N[1]*crhs32 + crhs33*crhs34 - crhs33*crhs35;
            rhs[7]=-DN(1,0)*crhs22 - DN(1,1)*crhs28 - DN(1,2)*crhs33 - N[1]*crhs17;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs19 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs11 - N[2]*crhs2 - N[2]*crhs3 + crhs22*crhs36 - crhs22*crhs37;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs19 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs27 + crhs28*crhs36 - crhs28*crhs37;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs19 - DN(2,2)*stress[2] + N[2]*crhs29 - N[2]*crhs30 - N[2]*crhs31 - N[2]*crhs32 + crhs33*crhs36 - crhs33*crhs37;
            rhs[11]=-DN(2,0)*crhs22 - DN(2,1)*crhs28 - DN(2,2)*crhs33 - N[2]*crhs17;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs19 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs11 - N[3]*crhs2 - N[3]*crhs3 + crhs22*crhs38 - crhs22*crhs39;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs19 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs27 + crhs28*crhs38 - crhs28*crhs39;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs19 - DN(3,2)*stress[2] + N[3]*crhs29 - N[3]*crhs30 - N[3]*crhs31 - N[3]*crhs32 + crhs33*crhs38 - crhs33*crhs39;
            rhs[15]=-DN(3,0)*crhs22 - DN(3,1)*crhs28 - DN(3,2)*crhs33 - N[3]*crhs17;


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
    const double K_darcy = rData.DarcyTerm;

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
const double cV2 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 =             1.0*DNenr(0,0)*K_darcy*cV2;
const double cV4 =             DN(0,0)*cV0 + DN(0,1)*cV1;
const double cV5 =             1.0*DNenr(0,0)*cV2*rho;
const double cV6 =             1.0*DNenr(1,0)*K_darcy*cV2;
const double cV7 =             1.0*DNenr(1,0)*cV2*rho;
const double cV8 =             1.0*DNenr(2,0)*K_darcy*cV2;
const double cV9 =             1.0*DNenr(2,0)*cV2*rho;
const double cV10 =             1.0*DNenr(0,1)*K_darcy*cV2;
const double cV11 =             1.0*DNenr(0,1)*cV2*rho;
const double cV12 =             1.0*DNenr(1,1)*K_darcy*cV2;
const double cV13 =             1.0*DNenr(1,1)*cV2*rho;
const double cV14 =             1.0*DNenr(2,1)*K_darcy*cV2;
const double cV15 =             1.0*DNenr(2,1)*cV2*rho;
const double cV16 =             1.0*cV2;
const double cV17 =             DN(1,0)*cV0 + DN(1,1)*cV1;
const double cV18 =             DN(2,0)*cV0 + DN(2,1)*cV1;
            V(0,0)=-DN(0,0)*Nenr[0] - N[0]*cV3 + cV4*cV5;
            V(0,1)=-DN(0,0)*Nenr[1] - N[0]*cV6 + cV4*cV7;
            V(0,2)=-DN(0,0)*Nenr[2] - N[0]*cV8 + cV4*cV9;
            V(1,0)=-DN(0,1)*Nenr[0] - N[0]*cV10 + cV11*cV4;
            V(1,1)=-DN(0,1)*Nenr[1] - N[0]*cV12 + cV13*cV4;
            V(1,2)=-DN(0,1)*Nenr[2] - N[0]*cV14 + cV15*cV4;
            V(2,0)=cV16*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV16*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV16*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] - N[1]*cV3 + cV17*cV5;
            V(3,1)=-DN(1,0)*Nenr[1] - N[1]*cV6 + cV17*cV7;
            V(3,2)=-DN(1,0)*Nenr[2] - N[1]*cV8 + cV17*cV9;
            V(4,0)=-DN(1,1)*Nenr[0] - N[1]*cV10 + cV11*cV17;
            V(4,1)=-DN(1,1)*Nenr[1] - N[1]*cV12 + cV13*cV17;
            V(4,2)=-DN(1,1)*Nenr[2] - N[1]*cV14 + cV15*cV17;
            V(5,0)=cV16*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV16*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV16*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] - N[2]*cV3 + cV18*cV5;
            V(6,1)=-DN(2,0)*Nenr[1] - N[2]*cV6 + cV18*cV7;
            V(6,2)=-DN(2,0)*Nenr[2] - N[2]*cV8 + cV18*cV9;
            V(7,0)=-DN(2,1)*Nenr[0] - N[2]*cV10 + cV11*cV18;
            V(7,1)=-DN(2,1)*Nenr[1] - N[2]*cV12 + cV13*cV18;
            V(7,2)=-DN(2,1)*Nenr[2] - N[2]*cV14 + cV15*cV18;
            V(8,0)=cV16*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV16*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV16*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 =             K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0);
const double cH3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 =             1.0*DNenr(0,0)*cH3;
const double cH5 =             1.0*DNenr(0,1)*cH3;
const double cH6 =             1.0*cH3;
const double cH7 =             K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0);
const double cH8 =             K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0);
const double cH9 =             1.0*DNenr(1,0)*cH3;
const double cH10 =             1.0*DNenr(1,1)*cH3;
const double cH11 =             1.0*DNenr(2,0)*cH3;
const double cH12 =             1.0*DNenr(2,1)*cH3;
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


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee5 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee6 =             1.0*crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 =             1.0*crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
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
    const double K_darcy = rData.DarcyTerm;

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
const double cV3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 =             1.0*DNenr(0,0)*K_darcy*cV3;
const double cV5 =             DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2;
const double cV6 =             1.0*DNenr(0,0)*cV3*rho;
const double cV7 =             1.0*DNenr(1,0)*K_darcy*cV3;
const double cV8 =             1.0*DNenr(1,0)*cV3*rho;
const double cV9 =             1.0*DNenr(2,0)*K_darcy*cV3;
const double cV10 =             1.0*DNenr(2,0)*cV3*rho;
const double cV11 =             1.0*DNenr(3,0)*K_darcy*cV3;
const double cV12 =             1.0*DNenr(3,0)*cV3*rho;
const double cV13 =             1.0*DNenr(0,1)*K_darcy*cV3;
const double cV14 =             1.0*DNenr(0,1)*cV3*rho;
const double cV15 =             1.0*DNenr(1,1)*K_darcy*cV3;
const double cV16 =             1.0*DNenr(1,1)*cV3*rho;
const double cV17 =             1.0*DNenr(2,1)*K_darcy*cV3;
const double cV18 =             1.0*DNenr(2,1)*cV3*rho;
const double cV19 =             1.0*DNenr(3,1)*K_darcy*cV3;
const double cV20 =             1.0*DNenr(3,1)*cV3*rho;
const double cV21 =             1.0*DNenr(0,2)*K_darcy*cV3;
const double cV22 =             1.0*DNenr(0,2)*cV3*rho;
const double cV23 =             1.0*DNenr(1,2)*K_darcy*cV3;
const double cV24 =             1.0*DNenr(1,2)*cV3*rho;
const double cV25 =             1.0*DNenr(2,2)*K_darcy*cV3;
const double cV26 =             1.0*DNenr(2,2)*cV3*rho;
const double cV27 =             1.0*DNenr(3,2)*K_darcy*cV3;
const double cV28 =             1.0*DNenr(3,2)*cV3*rho;
const double cV29 =             1.0*cV3;
const double cV30 =             DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2;
const double cV31 =             DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2;
const double cV32 =             DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2;
            V(0,0)=-DN(0,0)*Nenr[0] - N[0]*cV4 + cV5*cV6;
            V(0,1)=-DN(0,0)*Nenr[1] - N[0]*cV7 + cV5*cV8;
            V(0,2)=-DN(0,0)*Nenr[2] - N[0]*cV9 + cV10*cV5;
            V(0,3)=-DN(0,0)*Nenr[3] - N[0]*cV11 + cV12*cV5;
            V(1,0)=-DN(0,1)*Nenr[0] - N[0]*cV13 + cV14*cV5;
            V(1,1)=-DN(0,1)*Nenr[1] - N[0]*cV15 + cV16*cV5;
            V(1,2)=-DN(0,1)*Nenr[2] - N[0]*cV17 + cV18*cV5;
            V(1,3)=-DN(0,1)*Nenr[3] - N[0]*cV19 + cV20*cV5;
            V(2,0)=-DN(0,2)*Nenr[0] - N[0]*cV21 + cV22*cV5;
            V(2,1)=-DN(0,2)*Nenr[1] - N[0]*cV23 + cV24*cV5;
            V(2,2)=-DN(0,2)*Nenr[2] - N[0]*cV25 + cV26*cV5;
            V(2,3)=-DN(0,2)*Nenr[3] - N[0]*cV27 + cV28*cV5;
            V(3,0)=cV29*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV29*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV29*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV29*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] - N[1]*cV4 + cV30*cV6;
            V(4,1)=-DN(1,0)*Nenr[1] - N[1]*cV7 + cV30*cV8;
            V(4,2)=-DN(1,0)*Nenr[2] - N[1]*cV9 + cV10*cV30;
            V(4,3)=-DN(1,0)*Nenr[3] - N[1]*cV11 + cV12*cV30;
            V(5,0)=-DN(1,1)*Nenr[0] - N[1]*cV13 + cV14*cV30;
            V(5,1)=-DN(1,1)*Nenr[1] - N[1]*cV15 + cV16*cV30;
            V(5,2)=-DN(1,1)*Nenr[2] - N[1]*cV17 + cV18*cV30;
            V(5,3)=-DN(1,1)*Nenr[3] - N[1]*cV19 + cV20*cV30;
            V(6,0)=-DN(1,2)*Nenr[0] - N[1]*cV21 + cV22*cV30;
            V(6,1)=-DN(1,2)*Nenr[1] - N[1]*cV23 + cV24*cV30;
            V(6,2)=-DN(1,2)*Nenr[2] - N[1]*cV25 + cV26*cV30;
            V(6,3)=-DN(1,2)*Nenr[3] - N[1]*cV27 + cV28*cV30;
            V(7,0)=cV29*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV29*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV29*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV29*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] - N[2]*cV4 + cV31*cV6;
            V(8,1)=-DN(2,0)*Nenr[1] - N[2]*cV7 + cV31*cV8;
            V(8,2)=-DN(2,0)*Nenr[2] - N[2]*cV9 + cV10*cV31;
            V(8,3)=-DN(2,0)*Nenr[3] - N[2]*cV11 + cV12*cV31;
            V(9,0)=-DN(2,1)*Nenr[0] - N[2]*cV13 + cV14*cV31;
            V(9,1)=-DN(2,1)*Nenr[1] - N[2]*cV15 + cV16*cV31;
            V(9,2)=-DN(2,1)*Nenr[2] - N[2]*cV17 + cV18*cV31;
            V(9,3)=-DN(2,1)*Nenr[3] - N[2]*cV19 + cV20*cV31;
            V(10,0)=-DN(2,2)*Nenr[0] - N[2]*cV21 + cV22*cV31;
            V(10,1)=-DN(2,2)*Nenr[1] - N[2]*cV23 + cV24*cV31;
            V(10,2)=-DN(2,2)*Nenr[2] - N[2]*cV25 + cV26*cV31;
            V(10,3)=-DN(2,2)*Nenr[3] - N[2]*cV27 + cV28*cV31;
            V(11,0)=cV29*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV29*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV29*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV29*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] - N[3]*cV4 + cV32*cV6;
            V(12,1)=-DN(3,0)*Nenr[1] - N[3]*cV7 + cV32*cV8;
            V(12,2)=-DN(3,0)*Nenr[2] - N[3]*cV9 + cV10*cV32;
            V(12,3)=-DN(3,0)*Nenr[3] - N[3]*cV11 + cV12*cV32;
            V(13,0)=-DN(3,1)*Nenr[0] - N[3]*cV13 + cV14*cV32;
            V(13,1)=-DN(3,1)*Nenr[1] - N[3]*cV15 + cV16*cV32;
            V(13,2)=-DN(3,1)*Nenr[2] - N[3]*cV17 + cV18*cV32;
            V(13,3)=-DN(3,1)*Nenr[3] - N[3]*cV19 + cV20*cV32;
            V(14,0)=-DN(3,2)*Nenr[0] - N[3]*cV21 + cV22*cV32;
            V(14,1)=-DN(3,2)*Nenr[1] - N[3]*cV23 + cV24*cV32;
            V(14,2)=-DN(3,2)*Nenr[2] - N[3]*cV25 + cV26*cV32;
            V(14,3)=-DN(3,2)*Nenr[3] - N[3]*cV27 + cV28*cV32;
            V(15,0)=cV29*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV29*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV29*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV29*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 =             K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0);
const double cH4 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH5 =             1.0*DNenr(0,0)*cH4;
const double cH6 =             1.0*DNenr(0,1)*cH4;
const double cH7 =             1.0*DNenr(0,2)*cH4;
const double cH8 =             1.0*cH4;
const double cH9 =             K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0);
const double cH10 =             K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0);
const double cH11 =             K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0);
const double cH12 =             1.0*DNenr(1,0)*cH4;
const double cH13 =             1.0*DNenr(1,1)*cH4;
const double cH14 =             1.0*DNenr(1,2)*cH4;
const double cH15 =             1.0*DNenr(2,0)*cH4;
const double cH16 =             1.0*DNenr(2,1)*cH4;
const double cH17 =             1.0*DNenr(2,2)*cH4;
const double cH18 =             1.0*DNenr(3,0)*cH4;
const double cH19 =             1.0*DNenr(3,1)*cH4;
const double cH20 =             1.0*DNenr(3,2)*cH4;
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


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee13 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee11, 2) + pow(crhs_ee12, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee14 =             1.0*crhs_ee13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee10*(crhs_ee1 + crhs_ee3 + crhs_ee5 + crhs_ee7) + crhs_ee11*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee12*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee15 =             1.0*crhs_ee13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee10*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee11*(crhs_ee2 + crhs_ee4 + crhs_ee6 + crhs_ee8) + crhs_ee12*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee16 =             1.0*crhs_ee13*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee0*crhs_ee12 + crhs_ee10*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee11*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
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

    //We only enrich elements which are not almost empty/full
    if (positive_volume / Vol > min_area_ratio && negative_volume / Vol > min_area_ratio) {

        // Compute the maximum diagonal value in the enrichment stiffness matrix
        double max_diag = 0.0;
        for (unsigned int k = 0; k < NumNodes; ++k){
            if (std::abs(rKeeTot(k, k)) > max_diag){
                max_diag = std::abs(rKeeTot(k, k));
            }
        }
        if (max_diag == 0.0){
            max_diag = 1.0;
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


template <class TElementData>
void TwoFluidNavierStokes<TElementData>::GetValueOnIntegrationPoints(   const Variable<double> &rVariable,
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

template class TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>;

} // namespace Kratos