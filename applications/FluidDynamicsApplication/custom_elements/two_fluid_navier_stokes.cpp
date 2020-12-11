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
#pragma GCC optimize("Ofast")
#pragma GCC target "avx"

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
    const ProcessInfo &rCurrentProcessInfo)
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
            Matrix int_shape_function_neg;                                       // interface shape functions
            Matrix int_shape_function_enr_neg, int_shape_function_enr_pos;       // interface enriched shape functions
            GeometryType::ShapeFunctionsGradientsType int_shape_derivatives_neg; // interface shape functions derivatives
            Vector int_gauss_pts_weights;                                // interface Gauss points weights
            std::vector<Vector> int_normals_neg;                                 // interface normal vector based on the negative side
            Vector gauss_pts_curvature;                                  // curvatures calculated on interface Gauss points

            ModifiedShapeFunctions::Pointer p_modified_sh_func = pGetModifiedShapeFunctionsUtility(p_geom, data.Distance);

            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
                shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg,
                p_modified_sh_func);

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side
                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                const unsigned int number_of_gauss_points = gauss_weights.size();
                array_1d<double, NumNodes> Ncenter;
                for (unsigned int i = 0; i < NumNodes; ++i){
                    Ncenter[i] = 1.0/NumNodes;
                }
                for (unsigned int g = 0; g < number_of_gauss_points; ++g){
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

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); ++g_pos){
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

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); ++g_neg){
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

                if (rCurrentProcessInfo[SURFACE_TENSION]){

                    AddSurfaceTensionContribution(
                        data,
                        p_modified_sh_func,
                        rLeftHandSideMatrix,
                        rRightHandSideVector,
                        Htot,
                        Vtot,
                        Kee_tot,
                        rhs_ee_tot
                    );

                } else{
                    // Without pressure gradient stabilization, volume ratio is checked during condensation
                    // Also, without surface tension, zero pressure difference is penalized
                    CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
                }

            }
        } else {
            //Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            ShapeFunctionDerivativesArrayType shape_derivatives;
            this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
            const unsigned int number_of_gauss_points = gauss_weights.size();
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; ++g){
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
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

    //const auto vconv = rData.Velocity - rData.MeshVelocity;
    const double vconv_0_0 = rData.Velocity(0,0) - rData.MeshVelocity(0,0); const double vconv_0_1 = rData.Velocity(0,1) - rData.MeshVelocity(0,1); const double vconv_0_2 = rData.Velocity(0,2) - rData.MeshVelocity(0,2);
    const double vconv_1_0 = rData.Velocity(1,0) - rData.MeshVelocity(1,0); const double vconv_1_1 = rData.Velocity(1,1) - rData.MeshVelocity(1,1); const double vconv_1_2 = rData.Velocity(1,2) - rData.MeshVelocity(1,2);
    const double vconv_2_0 = rData.Velocity(2,0) - rData.MeshVelocity(2,0); const double vconv_2_1 = rData.Velocity(2,1) - rData.MeshVelocity(2,1); const double vconv_2_2 = rData.Velocity(2,2) - rData.MeshVelocity(2,2);
    const double vconv_3_0 = rData.Velocity(3,0) - rData.MeshVelocity(3,0); const double vconv_3_1 = rData.Velocity(3,1) - rData.MeshVelocity(3,1); const double vconv_3_2 = rData.Velocity(3,2) - rData.MeshVelocity(3,2);


    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    const double DN_0_0 = DN(0,0); const double DN_0_1 = DN(0,1); 
    const double DN_1_0 = DN(1,0); const double DN_1_1 = DN(1,1); 
    const double DN_2_0 = DN(2,0); const double DN_2_1 = DN(2,1); 

    const double N_0 = N[0];
    const double N_1 = N[1];
    const double N_2 = N[2];

    const double C_0_0 = C(0,0); const double C_0_1 = C(0,1); const double C_0_2 = C(0,2);
    const double C_1_0 = C(1,0); const double C_1_1 = C(1,1); const double C_1_2 = C(1,2);
    const double C_2_0 = C(2,0); const double C_2_1 = C(2,1); const double C_2_2 = C(2,2);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const auto Weight = rData.Weight;
    auto& lhs = rLHS;

    const double clhs0 =             C_0_0*DN_0_0 + C_0_2*DN_0_1;
const double clhs1 =             C_0_2*DN_0_0;
const double clhs2 =             C_2_2*DN_0_1 + clhs1;
const double clhs3 =             pow(DN_0_0, 2);
const double clhs4 =             N_0*vconv_0_0 + N_1*vconv_1_0 + N_2*vconv_2_0;
const double clhs5 =             N_0*vconv_0_1 + N_1*vconv_1_1 + N_2*vconv_2_1;
const double clhs6 =             rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             pow(N_0, 2);
const double clhs9 =             rho*(DN_0_0*clhs4 + DN_0_1*clhs5);
const double clhs10 =             bdf0*rho;
const double clhs11 =             K_darcy*N_0;
const double clhs12 =             N_0*clhs10;
const double clhs13 =             clhs11 + clhs12 + clhs9;
const double clhs14 =             1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 =             1.0*clhs9;
const double clhs16 =             clhs14*clhs15;
const double clhs17 =             1.0*clhs11;
const double clhs18 =             clhs14*clhs17;
const double clhs19 =             K_darcy*clhs8 + N_0*clhs9 + clhs10*clhs8 + clhs13*clhs16 - clhs13*clhs18;
const double clhs20 =             C_0_1*DN_0_1 + clhs1;
const double clhs21 =             C_1_2*DN_0_1;
const double clhs22 =             C_2_2*DN_0_0 + clhs21;
const double clhs23 =             DN_0_0*clhs7;
const double clhs24 =             DN_0_1*clhs23;
const double clhs25 =             Weight*(N_0 - clhs16 + clhs18);
const double clhs26 =             C_0_0*DN_1_0 + C_0_2*DN_1_1;
const double clhs27 =             C_0_2*DN_1_0;
const double clhs28 =             C_2_2*DN_1_1 + clhs27;
const double clhs29 =             DN_0_0*DN_1_0;
const double clhs30 =             N_1*clhs11 + N_1*clhs12;
const double clhs31 =             clhs29*clhs7 + clhs30;
const double clhs32 =             rho*(DN_1_0*clhs4 + DN_1_1*clhs5);
const double clhs33 =             K_darcy*N_1;
const double clhs34 =             N_1*clhs10;
const double clhs35 =             clhs32 + clhs33 + clhs34;
const double clhs36 =             N_0*clhs32 + clhs16*clhs35 - clhs18*clhs35;
const double clhs37 =             C_0_1*DN_1_1 + clhs27;
const double clhs38 =             C_1_2*DN_1_1;
const double clhs39 =             C_2_2*DN_1_0 + clhs38;
const double clhs40 =             DN_1_1*clhs23;
const double clhs41 =             DN_0_0*N_1;
const double clhs42 =             C_0_0*DN_2_0 + C_0_2*DN_2_1;
const double clhs43 =             C_0_2*DN_2_0;
const double clhs44 =             C_2_2*DN_2_1 + clhs43;
const double clhs45 =             DN_0_0*DN_2_0;
const double clhs46 =             N_2*clhs11 + N_2*clhs12;
const double clhs47 =             clhs45*clhs7 + clhs46;
const double clhs48 =             rho*(DN_2_0*clhs4 + DN_2_1*clhs5);
const double clhs49 =             K_darcy*N_2;
const double clhs50 =             N_2*clhs10;
const double clhs51 =             clhs48 + clhs49 + clhs50;
const double clhs52 =             N_0*clhs48 + clhs16*clhs51 - clhs18*clhs51;
const double clhs53 =             C_0_1*DN_2_1 + clhs43;
const double clhs54 =             C_1_2*DN_2_1;
const double clhs55 =             C_2_2*DN_2_0 + clhs54;
const double clhs56 =             DN_2_1*clhs23;
const double clhs57 =             DN_0_0*N_2;
const double clhs58 =             C_0_1*DN_0_0 + clhs21;
const double clhs59 =             C_1_1*DN_0_1 + C_1_2*DN_0_0;
const double clhs60 =             pow(DN_0_1, 2);
const double clhs61 =             C_0_1*DN_1_0 + clhs38;
const double clhs62 =             DN_0_1*clhs7;
const double clhs63 =             DN_1_0*clhs62;
const double clhs64 =             C_1_1*DN_1_1 + C_1_2*DN_1_0;
const double clhs65 =             DN_0_1*DN_1_1;
const double clhs66 =             clhs30 + clhs65*clhs7;
const double clhs67 =             DN_0_1*N_1;
const double clhs68 =             C_0_1*DN_2_0 + clhs54;
const double clhs69 =             DN_2_0*clhs62;
const double clhs70 =             C_1_1*DN_2_1 + C_1_2*DN_2_0;
const double clhs71 =             DN_0_1*DN_2_1;
const double clhs72 =             clhs46 + clhs7*clhs71;
const double clhs73 =             DN_0_1*N_2;
const double clhs74 =             Weight*(N_0 + clhs14*(1.0*clhs12 + clhs15 + clhs17));
const double clhs75 =             1.0*clhs14;
const double clhs76 =             Weight*clhs75;
const double clhs77 =             DN_1_0*N_0;
const double clhs78 =             clhs35*clhs75;
const double clhs79 =             DN_1_1*N_0;
const double clhs80 =             clhs76*(clhs29 + clhs65);
const double clhs81 =             DN_2_0*N_0;
const double clhs82 =             clhs51*clhs75;
const double clhs83 =             DN_2_1*N_0;
const double clhs84 =             clhs76*(clhs45 + clhs71);
const double clhs85 =             1.0*clhs32;
const double clhs86 =             clhs14*clhs85;
const double clhs87 =             1.0*clhs33;
const double clhs88 =             clhs14*clhs87;
const double clhs89 =             N_1*clhs9 + clhs13*clhs86 - clhs13*clhs88;
const double clhs90 =             pow(DN_1_0, 2);
const double clhs91 =             pow(N_1, 2);
const double clhs92 =             K_darcy*clhs91 + N_1*clhs32 + clhs10*clhs91 + clhs35*clhs86 - clhs35*clhs88;
const double clhs93 =             DN_1_0*clhs7;
const double clhs94 =             DN_1_1*clhs93;
const double clhs95 =             Weight*(N_1 - clhs86 + clhs88);
const double clhs96 =             DN_1_0*DN_2_0;
const double clhs97 =             N_2*clhs33 + N_2*clhs34;
const double clhs98 =             clhs7*clhs96 + clhs97;
const double clhs99 =             N_1*clhs48 + clhs51*clhs86 - clhs51*clhs88;
const double clhs100 =             DN_2_1*clhs93;
const double clhs101 =             DN_1_0*N_2;
const double clhs102 =             pow(DN_1_1, 2);
const double clhs103 =             DN_2_0*clhs7;
const double clhs104 =             DN_1_1*clhs103;
const double clhs105 =             DN_1_1*DN_2_1;
const double clhs106 =             clhs105*clhs7 + clhs97;
const double clhs107 =             DN_1_1*N_2;
const double clhs108 =             clhs13*clhs75;
const double clhs109 =             Weight*(N_1 + clhs14*(1.0*clhs34 + clhs85 + clhs87));
const double clhs110 =             DN_2_0*N_1;
const double clhs111 =             DN_2_1*N_1;
const double clhs112 =             clhs76*(clhs105 + clhs96);
const double clhs113 =             1.0*clhs48;
const double clhs114 =             clhs113*clhs14;
const double clhs115 =             1.0*clhs49;
const double clhs116 =             clhs115*clhs14;
const double clhs117 =             N_2*clhs9 + clhs114*clhs13 - clhs116*clhs13;
const double clhs118 =             N_2*clhs32 + clhs114*clhs35 - clhs116*clhs35;
const double clhs119 =             pow(DN_2_0, 2);
const double clhs120 =             pow(N_2, 2);
const double clhs121 =             K_darcy*clhs120 + N_2*clhs48 + clhs10*clhs120 + clhs114*clhs51 - clhs116*clhs51;
const double clhs122 =             DN_2_1*clhs103;
const double clhs123 =             Weight*(N_2 - clhs114 + clhs116);
const double clhs124 =             pow(DN_2_1, 2);
const double clhs125 =             Weight*(N_2 + clhs14*(clhs113 + clhs115 + 1.0*clhs50));
            lhs(0,0)+=Weight*(DN_0_0*clhs0 + DN_0_1*clhs2 + clhs19 + clhs3*clhs7);
            lhs(0,1)+=Weight*(DN_0_0*clhs20 + DN_0_1*clhs22 + clhs24);
            lhs(0,2)+=-DN_0_0*clhs25;
            lhs(0,3)+=Weight*(DN_0_0*clhs26 + DN_0_1*clhs28 + clhs31 + clhs36);
            lhs(0,4)+=Weight*(DN_0_0*clhs37 + DN_0_1*clhs39 + clhs40);
            lhs(0,5)+=-Weight*(-DN_1_0*clhs16 + DN_1_0*clhs18 + clhs41);
            lhs(0,6)+=Weight*(DN_0_0*clhs42 + DN_0_1*clhs44 + clhs47 + clhs52);
            lhs(0,7)+=Weight*(DN_0_0*clhs53 + DN_0_1*clhs55 + clhs56);
            lhs(0,8)+=-Weight*(-DN_2_0*clhs16 + DN_2_0*clhs18 + clhs57);
            lhs(1,0)+=Weight*(DN_0_0*clhs2 + DN_0_1*clhs58 + clhs24);
            lhs(1,1)+=Weight*(DN_0_0*clhs22 + DN_0_1*clhs59 + clhs19 + clhs60*clhs7);
            lhs(1,2)+=-DN_0_1*clhs25;
            lhs(1,3)+=Weight*(DN_0_0*clhs28 + DN_0_1*clhs61 + clhs63);
            lhs(1,4)+=Weight*(DN_0_0*clhs39 + DN_0_1*clhs64 + clhs36 + clhs66);
            lhs(1,5)+=-Weight*(-DN_1_1*clhs16 + DN_1_1*clhs18 + clhs67);
            lhs(1,6)+=Weight*(DN_0_0*clhs44 + DN_0_1*clhs68 + clhs69);
            lhs(1,7)+=Weight*(DN_0_0*clhs55 + DN_0_1*clhs70 + clhs52 + clhs72);
            lhs(1,8)+=-Weight*(-DN_2_1*clhs16 + DN_2_1*clhs18 + clhs73);
            lhs(2,0)+=DN_0_0*clhs74;
            lhs(2,1)+=DN_0_1*clhs74;
            lhs(2,2)+=clhs76*(clhs3 + clhs60);
            lhs(2,3)+=Weight*(DN_0_0*clhs78 + clhs77);
            lhs(2,4)+=Weight*(DN_0_1*clhs78 + clhs79);
            lhs(2,5)+=clhs80;
            lhs(2,6)+=Weight*(DN_0_0*clhs82 + clhs81);
            lhs(2,7)+=Weight*(DN_0_1*clhs82 + clhs83);
            lhs(2,8)+=clhs84;
            lhs(3,0)+=Weight*(DN_1_0*clhs0 + DN_1_1*clhs2 + clhs31 + clhs89);
            lhs(3,1)+=Weight*(DN_1_0*clhs20 + DN_1_1*clhs22 + clhs63);
            lhs(3,2)+=-Weight*(-DN_0_0*clhs86 + DN_0_0*clhs88 + clhs77);
            lhs(3,3)+=Weight*(DN_1_0*clhs26 + DN_1_1*clhs28 + clhs7*clhs90 + clhs92);
            lhs(3,4)+=Weight*(DN_1_0*clhs37 + DN_1_1*clhs39 + clhs94);
            lhs(3,5)+=-DN_1_0*clhs95;
            lhs(3,6)+=Weight*(DN_1_0*clhs42 + DN_1_1*clhs44 + clhs98 + clhs99);
            lhs(3,7)+=Weight*(DN_1_0*clhs53 + DN_1_1*clhs55 + clhs100);
            lhs(3,8)+=-Weight*(-DN_2_0*clhs86 + DN_2_0*clhs88 + clhs101);
            lhs(4,0)+=Weight*(DN_1_0*clhs2 + DN_1_1*clhs58 + clhs40);
            lhs(4,1)+=Weight*(DN_1_0*clhs22 + DN_1_1*clhs59 + clhs66 + clhs89);
            lhs(4,2)+=-Weight*(-DN_0_1*clhs86 + DN_0_1*clhs88 + clhs79);
            lhs(4,3)+=Weight*(DN_1_0*clhs28 + DN_1_1*clhs61 + clhs94);
            lhs(4,4)+=Weight*(DN_1_0*clhs39 + DN_1_1*clhs64 + clhs102*clhs7 + clhs92);
            lhs(4,5)+=-DN_1_1*clhs95;
            lhs(4,6)+=Weight*(DN_1_0*clhs44 + DN_1_1*clhs68 + clhs104);
            lhs(4,7)+=Weight*(DN_1_0*clhs55 + DN_1_1*clhs70 + clhs106 + clhs99);
            lhs(4,8)+=-Weight*(-DN_2_1*clhs86 + DN_2_1*clhs88 + clhs107);
            lhs(5,0)+=Weight*(DN_1_0*clhs108 + clhs41);
            lhs(5,1)+=Weight*(DN_1_1*clhs108 + clhs67);
            lhs(5,2)+=clhs80;
            lhs(5,3)+=DN_1_0*clhs109;
            lhs(5,4)+=DN_1_1*clhs109;
            lhs(5,5)+=clhs76*(clhs102 + clhs90);
            lhs(5,6)+=Weight*(DN_1_0*clhs82 + clhs110);
            lhs(5,7)+=Weight*(DN_1_1*clhs82 + clhs111);
            lhs(5,8)+=clhs112;
            lhs(6,0)+=Weight*(DN_2_0*clhs0 + DN_2_1*clhs2 + clhs117 + clhs47);
            lhs(6,1)+=Weight*(DN_2_0*clhs20 + DN_2_1*clhs22 + clhs69);
            lhs(6,2)+=-Weight*(-DN_0_0*clhs114 + DN_0_0*clhs116 + clhs81);
            lhs(6,3)+=Weight*(DN_2_0*clhs26 + DN_2_1*clhs28 + clhs118 + clhs98);
            lhs(6,4)+=Weight*(DN_2_0*clhs37 + DN_2_1*clhs39 + clhs104);
            lhs(6,5)+=-Weight*(-DN_1_0*clhs114 + DN_1_0*clhs116 + clhs110);
            lhs(6,6)+=Weight*(DN_2_0*clhs42 + DN_2_1*clhs44 + clhs119*clhs7 + clhs121);
            lhs(6,7)+=Weight*(DN_2_0*clhs53 + DN_2_1*clhs55 + clhs122);
            lhs(6,8)+=-DN_2_0*clhs123;
            lhs(7,0)+=Weight*(DN_2_0*clhs2 + DN_2_1*clhs58 + clhs56);
            lhs(7,1)+=Weight*(DN_2_0*clhs22 + DN_2_1*clhs59 + clhs117 + clhs72);
            lhs(7,2)+=-Weight*(-DN_0_1*clhs114 + DN_0_1*clhs116 + clhs83);
            lhs(7,3)+=Weight*(DN_2_0*clhs28 + DN_2_1*clhs61 + clhs100);
            lhs(7,4)+=Weight*(DN_2_0*clhs39 + DN_2_1*clhs64 + clhs106 + clhs118);
            lhs(7,5)+=-Weight*(-DN_1_1*clhs114 + DN_1_1*clhs116 + clhs111);
            lhs(7,6)+=Weight*(DN_2_0*clhs44 + DN_2_1*clhs68 + clhs122);
            lhs(7,7)+=Weight*(DN_2_0*clhs55 + DN_2_1*clhs70 + clhs121 + clhs124*clhs7);
            lhs(7,8)+=-DN_2_1*clhs123;
            lhs(8,0)+=Weight*(DN_2_0*clhs108 + clhs57);
            lhs(8,1)+=Weight*(DN_2_1*clhs108 + clhs73);
            lhs(8,2)+=clhs84;
            lhs(8,3)+=Weight*(DN_2_0*clhs78 + clhs101);
            lhs(8,4)+=Weight*(DN_2_1*clhs78 + clhs107);
            lhs(8,5)+=clhs112;
            lhs(8,6)+=DN_2_0*clhs125;
            lhs(8,7)+=DN_2_1*clhs125;
            lhs(8,8)+=clhs76*(clhs119 + clhs124);


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

    //const auto vconv = rData.Velocity - rData.MeshVelocity;
    const double vconv_0_0 = rData.Velocity(0,0) - rData.MeshVelocity(0,0); const double vconv_0_1 = rData.Velocity(0,1) - rData.MeshVelocity(0,1); const double vconv_0_2 = rData.Velocity(0,2) - rData.MeshVelocity(0,2);
    const double vconv_1_0 = rData.Velocity(1,0) - rData.MeshVelocity(1,0); const double vconv_1_1 = rData.Velocity(1,1) - rData.MeshVelocity(1,1); const double vconv_1_2 = rData.Velocity(1,2) - rData.MeshVelocity(1,2);
    const double vconv_2_0 = rData.Velocity(2,0) - rData.MeshVelocity(2,0); const double vconv_2_1 = rData.Velocity(2,1) - rData.MeshVelocity(2,1); const double vconv_2_2 = rData.Velocity(2,2) - rData.MeshVelocity(2,2);
    const double vconv_3_0 = rData.Velocity(3,0) - rData.MeshVelocity(3,0); const double vconv_3_1 = rData.Velocity(3,1) - rData.MeshVelocity(3,1); const double vconv_3_2 = rData.Velocity(3,2) - rData.MeshVelocity(3,2);

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    const double DN_0_0 = DN(0,0); const double DN_0_1 = DN(0,1); const double DN_0_2 = DN(0,2);
    const double DN_1_0 = DN(1,0); const double DN_1_1 = DN(1,1); const double DN_1_2 = DN(1,2);
    const double DN_2_0 = DN(2,0); const double DN_2_1 = DN(2,1); const double DN_2_2 = DN(2,2);
    const double DN_3_0 = DN(3,0); const double DN_3_1 = DN(3,1); const double DN_3_2 = DN(3,2);

    const double N_0 = N[0];
    const double N_1 = N[1];
    const double N_2 = N[2];
    const double N_3 = N[3];

    const double C_0_0 = C(0,0); const double C_0_1 = C(0,1); const double C_0_2 = C(0,2); const double C_0_3 = C(0,3); const double C_0_4 = C(0,4); const double C_0_5 = C(0,5);
    const double C_1_0 = C(1,0); const double C_1_1 = C(1,1); const double C_1_2 = C(1,2); const double C_1_3 = C(1,3); const double C_1_4 = C(1,4); const double C_1_5 = C(1,5);
    const double C_2_0 = C(2,0); const double C_2_1 = C(2,1); const double C_2_2 = C(2,2); const double C_2_3 = C(2,3); const double C_2_4 = C(2,4); const double C_2_5 = C(2,5);
    const double C_3_0 = C(3,0); const double C_3_1 = C(3,1); const double C_3_2 = C(3,2); const double C_3_3 = C(3,3); const double C_3_4 = C(3,4); const double C_3_5 = C(3,5);
    const double C_4_0 = C(4,0); const double C_4_1 = C(4,1); const double C_4_2 = C(4,2); const double C_4_3 = C(4,3); const double C_4_4 = C(4,4); const double C_4_5 = C(4,5);
    const double C_5_0 = C(5,0); const double C_5_1 = C(5,1); const double C_5_2 = C(5,2); const double C_5_3 = C(5,3); const double C_5_4 = C(5,4); const double C_5_5 = C(5,5);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double Weight = rData.Weight;
    auto &lhs = rLHS;

    const double clhs0 =             C_0_0*DN_0_0 + C_0_3*DN_0_1 + C_0_5*DN_0_2;
const double clhs1 =             C_0_3*DN_0_0;
const double clhs2 =             C_3_3*DN_0_1 + C_3_5*DN_0_2 + clhs1;
const double clhs3 =             C_0_5*DN_0_0;
const double clhs4 =             C_3_5*DN_0_1 + C_5_5*DN_0_2 + clhs3;
const double clhs5 =             pow(DN_0_0, 2);
const double clhs6 =             N_0*vconv_0_0 + N_1*vconv_1_0 + N_2*vconv_2_0 + N_3*vconv_3_0;
const double clhs7 =             N_0*vconv_0_1 + N_1*vconv_1_1 + N_2*vconv_2_1 + N_3*vconv_3_1;
const double clhs8 =             N_0*vconv_0_2 + N_1*vconv_1_2 + N_2*vconv_2_2 + N_3*vconv_3_2;
const double clhs9 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs9*h/stab_c1 + mu;
const double clhs11 =             pow(N_0, 2);
const double clhs12 =             rho*(DN_0_0*clhs6 + DN_0_1*clhs7 + DN_0_2*clhs8);
const double clhs13 =             bdf0*rho;
const double clhs14 =             K_darcy*N_0;
const double clhs15 =             N_0*clhs13;
const double clhs16 =             clhs12 + clhs14 + clhs15;
const double clhs17 =             1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 =             1.0*clhs12;
const double clhs19 =             clhs17*clhs18;
const double clhs20 =             1.0*clhs14;
const double clhs21 =             clhs17*clhs20;
const double clhs22 =             K_darcy*clhs11 + N_0*clhs12 + clhs11*clhs13 + clhs16*clhs19 - clhs16*clhs21;
const double clhs23 =             C_0_1*DN_0_1 + C_0_4*DN_0_2 + clhs1;
const double clhs24 =             C_1_3*DN_0_1;
const double clhs25 =             C_3_3*DN_0_0 + C_3_4*DN_0_2 + clhs24;
const double clhs26 =             C_3_5*DN_0_0;
const double clhs27 =             C_4_5*DN_0_2;
const double clhs28 =             C_1_5*DN_0_1 + clhs26 + clhs27;
const double clhs29 =             DN_0_0*clhs10;
const double clhs30 =             DN_0_1*clhs29;
const double clhs31 =             C_0_2*DN_0_2 + C_0_4*DN_0_1 + clhs3;
const double clhs32 =             C_3_4*DN_0_1;
const double clhs33 =             C_2_3*DN_0_2 + clhs26 + clhs32;
const double clhs34 =             C_2_5*DN_0_2;
const double clhs35 =             C_4_5*DN_0_1 + C_5_5*DN_0_0 + clhs34;
const double clhs36 =             DN_0_2*clhs29;
const double clhs37 =             Weight*(N_0 - clhs19 + clhs21);
const double clhs38 =             C_0_0*DN_1_0 + C_0_3*DN_1_1 + C_0_5*DN_1_2;
const double clhs39 =             C_0_3*DN_1_0;
const double clhs40 =             C_3_3*DN_1_1 + C_3_5*DN_1_2 + clhs39;
const double clhs41 =             C_0_5*DN_1_0;
const double clhs42 =             C_3_5*DN_1_1 + C_5_5*DN_1_2 + clhs41;
const double clhs43 =             DN_0_0*DN_1_0;
const double clhs44 =             N_1*clhs14 + N_1*clhs15;
const double clhs45 =             clhs10*clhs43 + clhs44;
const double clhs46 =             rho*(DN_1_0*clhs6 + DN_1_1*clhs7 + DN_1_2*clhs8);
const double clhs47 =             K_darcy*N_1;
const double clhs48 =             N_1*clhs13;
const double clhs49 =             clhs46 + clhs47 + clhs48;
const double clhs50 =             N_0*clhs46 + clhs19*clhs49 - clhs21*clhs49;
const double clhs51 =             C_0_1*DN_1_1 + C_0_4*DN_1_2 + clhs39;
const double clhs52 =             C_1_3*DN_1_1;
const double clhs53 =             C_3_3*DN_1_0 + C_3_4*DN_1_2 + clhs52;
const double clhs54 =             C_3_5*DN_1_0;
const double clhs55 =             C_4_5*DN_1_2;
const double clhs56 =             C_1_5*DN_1_1 + clhs54 + clhs55;
const double clhs57 =             DN_1_1*clhs29;
const double clhs58 =             C_0_2*DN_1_2 + C_0_4*DN_1_1 + clhs41;
const double clhs59 =             C_3_4*DN_1_1;
const double clhs60 =             C_2_3*DN_1_2 + clhs54 + clhs59;
const double clhs61 =             C_2_5*DN_1_2;
const double clhs62 =             C_4_5*DN_1_1 + C_5_5*DN_1_0 + clhs61;
const double clhs63 =             DN_1_2*clhs29;
const double clhs64 =             DN_0_0*N_1;
const double clhs65 =             C_0_0*DN_2_0 + C_0_3*DN_2_1 + C_0_5*DN_2_2;
const double clhs66 =             C_0_3*DN_2_0;
const double clhs67 =             C_3_3*DN_2_1 + C_3_5*DN_2_2 + clhs66;
const double clhs68 =             C_0_5*DN_2_0;
const double clhs69 =             C_3_5*DN_2_1 + C_5_5*DN_2_2 + clhs68;
const double clhs70 =             DN_0_0*DN_2_0;
const double clhs71 =             N_2*clhs14 + N_2*clhs15;
const double clhs72 =             clhs10*clhs70 + clhs71;
const double clhs73 =             rho*(DN_2_0*clhs6 + DN_2_1*clhs7 + DN_2_2*clhs8);
const double clhs74 =             K_darcy*N_2;
const double clhs75 =             N_2*clhs13;
const double clhs76 =             clhs73 + clhs74 + clhs75;
const double clhs77 =             N_0*clhs73 + clhs19*clhs76 - clhs21*clhs76;
const double clhs78 =             C_0_1*DN_2_1 + C_0_4*DN_2_2 + clhs66;
const double clhs79 =             C_1_3*DN_2_1;
const double clhs80 =             C_3_3*DN_2_0 + C_3_4*DN_2_2 + clhs79;
const double clhs81 =             C_3_5*DN_2_0;
const double clhs82 =             C_4_5*DN_2_2;
const double clhs83 =             C_1_5*DN_2_1 + clhs81 + clhs82;
const double clhs84 =             DN_2_1*clhs29;
const double clhs85 =             C_0_2*DN_2_2 + C_0_4*DN_2_1 + clhs68;
const double clhs86 =             C_3_4*DN_2_1;
const double clhs87 =             C_2_3*DN_2_2 + clhs81 + clhs86;
const double clhs88 =             C_2_5*DN_2_2;
const double clhs89 =             C_4_5*DN_2_1 + C_5_5*DN_2_0 + clhs88;
const double clhs90 =             DN_2_2*clhs29;
const double clhs91 =             DN_0_0*N_2;
const double clhs92 =             C_0_0*DN_3_0 + C_0_3*DN_3_1 + C_0_5*DN_3_2;
const double clhs93 =             C_0_3*DN_3_0;
const double clhs94 =             C_3_3*DN_3_1 + C_3_5*DN_3_2 + clhs93;
const double clhs95 =             C_0_5*DN_3_0;
const double clhs96 =             C_3_5*DN_3_1 + C_5_5*DN_3_2 + clhs95;
const double clhs97 =             DN_0_0*DN_3_0;
const double clhs98 =             N_3*clhs14 + N_3*clhs15;
const double clhs99 =             clhs10*clhs97 + clhs98;
const double clhs100 =             rho*(DN_3_0*clhs6 + DN_3_1*clhs7 + DN_3_2*clhs8);
const double clhs101 =             K_darcy*N_3;
const double clhs102 =             N_3*clhs13;
const double clhs103 =             clhs100 + clhs101 + clhs102;
const double clhs104 =             N_0*clhs100 + clhs103*clhs19 - clhs103*clhs21;
const double clhs105 =             C_0_1*DN_3_1 + C_0_4*DN_3_2 + clhs93;
const double clhs106 =             C_1_3*DN_3_1;
const double clhs107 =             C_3_3*DN_3_0 + C_3_4*DN_3_2 + clhs106;
const double clhs108 =             C_3_5*DN_3_0;
const double clhs109 =             C_4_5*DN_3_2;
const double clhs110 =             C_1_5*DN_3_1 + clhs108 + clhs109;
const double clhs111 =             DN_3_1*clhs29;
const double clhs112 =             C_0_2*DN_3_2 + C_0_4*DN_3_1 + clhs95;
const double clhs113 =             C_3_4*DN_3_1;
const double clhs114 =             C_2_3*DN_3_2 + clhs108 + clhs113;
const double clhs115 =             C_2_5*DN_3_2;
const double clhs116 =             C_4_5*DN_3_1 + C_5_5*DN_3_0 + clhs115;
const double clhs117 =             DN_3_2*clhs29;
const double clhs118 =             DN_0_0*N_3;
const double clhs119 =             C_0_1*DN_0_0 + C_1_5*DN_0_2 + clhs24;
const double clhs120 =             C_0_4*DN_0_0 + clhs27 + clhs32;
const double clhs121 =             C_1_1*DN_0_1 + C_1_3*DN_0_0 + C_1_4*DN_0_2;
const double clhs122 =             C_1_4*DN_0_1;
const double clhs123 =             C_3_4*DN_0_0 + C_4_4*DN_0_2 + clhs122;
const double clhs124 =             pow(DN_0_1, 2);
const double clhs125 =             C_1_2*DN_0_2 + C_1_5*DN_0_0 + clhs122;
const double clhs126 =             C_2_4*DN_0_2;
const double clhs127 =             C_4_4*DN_0_1 + C_4_5*DN_0_0 + clhs126;
const double clhs128 =             DN_0_1*clhs10;
const double clhs129 =             DN_0_2*clhs128;
const double clhs130 =             C_0_1*DN_1_0 + C_1_5*DN_1_2 + clhs52;
const double clhs131 =             C_0_4*DN_1_0 + clhs55 + clhs59;
const double clhs132 =             DN_1_0*clhs128;
const double clhs133 =             C_1_1*DN_1_1 + C_1_3*DN_1_0 + C_1_4*DN_1_2;
const double clhs134 =             C_1_4*DN_1_1;
const double clhs135 =             C_3_4*DN_1_0 + C_4_4*DN_1_2 + clhs134;
const double clhs136 =             DN_0_1*DN_1_1;
const double clhs137 =             clhs10*clhs136;
const double clhs138 =             clhs44 + clhs50;
const double clhs139 =             C_1_2*DN_1_2 + C_1_5*DN_1_0 + clhs134;
const double clhs140 =             C_2_4*DN_1_2;
const double clhs141 =             C_4_4*DN_1_1 + C_4_5*DN_1_0 + clhs140;
const double clhs142 =             DN_1_2*clhs128;
const double clhs143 =             DN_0_1*N_1;
const double clhs144 =             C_0_1*DN_2_0 + C_1_5*DN_2_2 + clhs79;
const double clhs145 =             C_0_4*DN_2_0 + clhs82 + clhs86;
const double clhs146 =             DN_2_0*clhs128;
const double clhs147 =             C_1_1*DN_2_1 + C_1_3*DN_2_0 + C_1_4*DN_2_2;
const double clhs148 =             C_1_4*DN_2_1;
const double clhs149 =             C_3_4*DN_2_0 + C_4_4*DN_2_2 + clhs148;
const double clhs150 =             DN_0_1*DN_2_1;
const double clhs151 =             clhs10*clhs150;
const double clhs152 =             clhs71 + clhs77;
const double clhs153 =             C_1_2*DN_2_2 + C_1_5*DN_2_0 + clhs148;
const double clhs154 =             C_2_4*DN_2_2;
const double clhs155 =             C_4_4*DN_2_1 + C_4_5*DN_2_0 + clhs154;
const double clhs156 =             DN_2_2*clhs128;
const double clhs157 =             DN_0_1*N_2;
const double clhs158 =             C_0_1*DN_3_0 + C_1_5*DN_3_2 + clhs106;
const double clhs159 =             C_0_4*DN_3_0 + clhs109 + clhs113;
const double clhs160 =             DN_3_0*clhs128;
const double clhs161 =             C_1_1*DN_3_1 + C_1_3*DN_3_0 + C_1_4*DN_3_2;
const double clhs162 =             C_1_4*DN_3_1;
const double clhs163 =             C_3_4*DN_3_0 + C_4_4*DN_3_2 + clhs162;
const double clhs164 =             DN_0_1*DN_3_1;
const double clhs165 =             clhs10*clhs164;
const double clhs166 =             clhs104 + clhs98;
const double clhs167 =             C_1_2*DN_3_2 + C_1_5*DN_3_0 + clhs162;
const double clhs168 =             C_2_4*DN_3_2;
const double clhs169 =             C_4_4*DN_3_1 + C_4_5*DN_3_0 + clhs168;
const double clhs170 =             DN_3_2*clhs128;
const double clhs171 =             DN_0_1*N_3;
const double clhs172 =             C_0_2*DN_0_0 + C_2_3*DN_0_1 + clhs34;
const double clhs173 =             C_1_2*DN_0_1 + C_2_3*DN_0_0 + clhs126;
const double clhs174 =             C_2_2*DN_0_2 + C_2_4*DN_0_1 + C_2_5*DN_0_0;
const double clhs175 =             pow(DN_0_2, 2);
const double clhs176 =             C_0_2*DN_1_0 + C_2_3*DN_1_1 + clhs61;
const double clhs177 =             DN_0_2*clhs10;
const double clhs178 =             DN_1_0*clhs177;
const double clhs179 =             C_1_2*DN_1_1 + C_2_3*DN_1_0 + clhs140;
const double clhs180 =             DN_1_1*clhs177;
const double clhs181 =             C_2_2*DN_1_2 + C_2_4*DN_1_1 + C_2_5*DN_1_0;
const double clhs182 =             DN_0_2*DN_1_2;
const double clhs183 =             clhs10*clhs182;
const double clhs184 =             DN_0_2*N_1;
const double clhs185 =             C_0_2*DN_2_0 + C_2_3*DN_2_1 + clhs88;
const double clhs186 =             DN_2_0*clhs177;
const double clhs187 =             C_1_2*DN_2_1 + C_2_3*DN_2_0 + clhs154;
const double clhs188 =             DN_2_1*clhs177;
const double clhs189 =             C_2_2*DN_2_2 + C_2_4*DN_2_1 + C_2_5*DN_2_0;
const double clhs190 =             DN_0_2*DN_2_2;
const double clhs191 =             clhs10*clhs190;
const double clhs192 =             DN_0_2*N_2;
const double clhs193 =             C_0_2*DN_3_0 + C_2_3*DN_3_1 + clhs115;
const double clhs194 =             DN_3_0*clhs177;
const double clhs195 =             C_1_2*DN_3_1 + C_2_3*DN_3_0 + clhs168;
const double clhs196 =             DN_3_1*clhs177;
const double clhs197 =             C_2_2*DN_3_2 + C_2_4*DN_3_1 + C_2_5*DN_3_0;
const double clhs198 =             DN_0_2*DN_3_2;
const double clhs199 =             clhs10*clhs198;
const double clhs200 =             DN_0_2*N_3;
const double clhs201 =             Weight*(N_0 + clhs17*(1.0*clhs15 + clhs18 + clhs20));
const double clhs202 =             1.0*clhs17;
const double clhs203 =             Weight*clhs202;
const double clhs204 =             DN_1_0*N_0;
const double clhs205 =             clhs202*clhs49;
const double clhs206 =             DN_1_1*N_0;
const double clhs207 =             DN_1_2*N_0;
const double clhs208 =             clhs203*(clhs136 + clhs182 + clhs43);
const double clhs209 =             DN_2_0*N_0;
const double clhs210 =             clhs202*clhs76;
const double clhs211 =             DN_2_1*N_0;
const double clhs212 =             DN_2_2*N_0;
const double clhs213 =             clhs203*(clhs150 + clhs190 + clhs70);
const double clhs214 =             DN_3_0*N_0;
const double clhs215 =             clhs103*clhs202;
const double clhs216 =             DN_3_1*N_0;
const double clhs217 =             DN_3_2*N_0;
const double clhs218 =             clhs203*(clhs164 + clhs198 + clhs97);
const double clhs219 =             1.0*clhs46;
const double clhs220 =             clhs17*clhs219;
const double clhs221 =             1.0*clhs47;
const double clhs222 =             clhs17*clhs221;
const double clhs223 =             N_1*clhs12 + clhs16*clhs220 - clhs16*clhs222;
const double clhs224 =             pow(DN_1_0, 2);
const double clhs225 =             pow(N_1, 2);
const double clhs226 =             K_darcy*clhs225 + N_1*clhs46 + clhs13*clhs225 + clhs220*clhs49 - clhs222*clhs49;
const double clhs227 =             DN_1_0*clhs10;
const double clhs228 =             DN_1_1*clhs227;
const double clhs229 =             DN_1_2*clhs227;
const double clhs230 =             Weight*(N_1 - clhs220 + clhs222);
const double clhs231 =             DN_1_0*DN_2_0;
const double clhs232 =             N_2*clhs47 + N_2*clhs48;
const double clhs233 =             clhs10*clhs231 + clhs232;
const double clhs234 =             N_1*clhs73 + clhs220*clhs76 - clhs222*clhs76;
const double clhs235 =             DN_2_1*clhs227;
const double clhs236 =             DN_2_2*clhs227;
const double clhs237 =             DN_1_0*N_2;
const double clhs238 =             DN_1_0*DN_3_0;
const double clhs239 =             N_3*clhs47 + N_3*clhs48;
const double clhs240 =             clhs10*clhs238 + clhs239;
const double clhs241 =             N_1*clhs100 + clhs103*clhs220 - clhs103*clhs222;
const double clhs242 =             DN_3_1*clhs227;
const double clhs243 =             DN_3_2*clhs227;
const double clhs244 =             DN_1_0*N_3;
const double clhs245 =             clhs223 + clhs44;
const double clhs246 =             pow(DN_1_1, 2);
const double clhs247 =             DN_1_1*clhs10;
const double clhs248 =             DN_1_2*clhs247;
const double clhs249 =             DN_2_0*clhs247;
const double clhs250 =             DN_1_1*DN_2_1;
const double clhs251 =             clhs10*clhs250;
const double clhs252 =             clhs232 + clhs234;
const double clhs253 =             DN_2_2*clhs247;
const double clhs254 =             DN_1_1*N_2;
const double clhs255 =             DN_3_0*clhs247;
const double clhs256 =             DN_1_1*DN_3_1;
const double clhs257 =             clhs10*clhs256;
const double clhs258 =             clhs239 + clhs241;
const double clhs259 =             DN_3_2*clhs247;
const double clhs260 =             DN_1_1*N_3;
const double clhs261 =             pow(DN_1_2, 2);
const double clhs262 =             DN_1_2*clhs10;
const double clhs263 =             DN_2_0*clhs262;
const double clhs264 =             DN_2_1*clhs262;
const double clhs265 =             DN_1_2*DN_2_2;
const double clhs266 =             clhs10*clhs265;
const double clhs267 =             DN_1_2*N_2;
const double clhs268 =             DN_3_0*clhs262;
const double clhs269 =             DN_3_1*clhs262;
const double clhs270 =             DN_1_2*DN_3_2;
const double clhs271 =             clhs10*clhs270;
const double clhs272 =             DN_1_2*N_3;
const double clhs273 =             clhs16*clhs202;
const double clhs274 =             Weight*(N_1 + clhs17*(clhs219 + clhs221 + 1.0*clhs48));
const double clhs275 =             DN_2_0*N_1;
const double clhs276 =             DN_2_1*N_1;
const double clhs277 =             DN_2_2*N_1;
const double clhs278 =             clhs203*(clhs231 + clhs250 + clhs265);
const double clhs279 =             DN_3_0*N_1;
const double clhs280 =             DN_3_1*N_1;
const double clhs281 =             DN_3_2*N_1;
const double clhs282 =             clhs203*(clhs238 + clhs256 + clhs270);
const double clhs283 =             1.0*clhs73;
const double clhs284 =             clhs17*clhs283;
const double clhs285 =             1.0*clhs74;
const double clhs286 =             clhs17*clhs285;
const double clhs287 =             N_2*clhs12 + clhs16*clhs284 - clhs16*clhs286;
const double clhs288 =             N_2*clhs46 + clhs284*clhs49 - clhs286*clhs49;
const double clhs289 =             pow(DN_2_0, 2);
const double clhs290 =             pow(N_2, 2);
const double clhs291 =             K_darcy*clhs290 + N_2*clhs73 + clhs13*clhs290 + clhs284*clhs76 - clhs286*clhs76;
const double clhs292 =             DN_2_0*clhs10;
const double clhs293 =             DN_2_1*clhs292;
const double clhs294 =             DN_2_2*clhs292;
const double clhs295 =             Weight*(N_2 - clhs284 + clhs286);
const double clhs296 =             DN_2_0*DN_3_0;
const double clhs297 =             N_3*clhs74 + N_3*clhs75;
const double clhs298 =             clhs10*clhs296 + clhs297;
const double clhs299 =             N_2*clhs100 + clhs103*clhs284 - clhs103*clhs286;
const double clhs300 =             DN_3_1*clhs292;
const double clhs301 =             DN_3_2*clhs292;
const double clhs302 =             DN_2_0*N_3;
const double clhs303 =             clhs287 + clhs71;
const double clhs304 =             clhs232 + clhs288;
const double clhs305 =             pow(DN_2_1, 2);
const double clhs306 =             DN_2_1*clhs10;
const double clhs307 =             DN_2_2*clhs306;
const double clhs308 =             DN_3_0*clhs306;
const double clhs309 =             DN_2_1*DN_3_1;
const double clhs310 =             clhs10*clhs309;
const double clhs311 =             clhs297 + clhs299;
const double clhs312 =             DN_3_2*clhs306;
const double clhs313 =             DN_2_1*N_3;
const double clhs314 =             pow(DN_2_2, 2);
const double clhs315 =             DN_2_2*clhs10;
const double clhs316 =             DN_3_0*clhs315;
const double clhs317 =             DN_3_1*clhs315;
const double clhs318 =             DN_2_2*DN_3_2;
const double clhs319 =             clhs10*clhs318;
const double clhs320 =             DN_2_2*N_3;
const double clhs321 =             Weight*(N_2 + clhs17*(clhs283 + clhs285 + 1.0*clhs75));
const double clhs322 =             DN_3_0*N_2;
const double clhs323 =             DN_3_1*N_2;
const double clhs324 =             DN_3_2*N_2;
const double clhs325 =             clhs203*(clhs296 + clhs309 + clhs318);
const double clhs326 =             1.0*clhs100;
const double clhs327 =             clhs17*clhs326;
const double clhs328 =             1.0*clhs101;
const double clhs329 =             clhs17*clhs328;
const double clhs330 =             N_3*clhs12 + clhs16*clhs327 - clhs16*clhs329;
const double clhs331 =             N_3*clhs46 + clhs327*clhs49 - clhs329*clhs49;
const double clhs332 =             N_3*clhs73 + clhs327*clhs76 - clhs329*clhs76;
const double clhs333 =             pow(DN_3_0, 2);
const double clhs334 =             pow(N_3, 2);
const double clhs335 =             K_darcy*clhs334 + N_3*clhs100 + clhs103*clhs327 - clhs103*clhs329 + clhs13*clhs334;
const double clhs336 =             DN_3_0*clhs10;
const double clhs337 =             DN_3_1*clhs336;
const double clhs338 =             DN_3_2*clhs336;
const double clhs339 =             Weight*(N_3 - clhs327 + clhs329);
const double clhs340 =             clhs330 + clhs98;
const double clhs341 =             clhs239 + clhs331;
const double clhs342 =             clhs297 + clhs332;
const double clhs343 =             pow(DN_3_1, 2);
const double clhs344 =             DN_3_1*DN_3_2*clhs10;
const double clhs345 =             pow(DN_3_2, 2);
const double clhs346 =             Weight*(N_3 + clhs17*(1.0*clhs102 + clhs326 + clhs328));
            lhs(0,0)+=Weight*(DN_0_0*clhs0 + DN_0_1*clhs2 + DN_0_2*clhs4 + clhs10*clhs5 + clhs22);
            lhs(0,1)+=Weight*(DN_0_0*clhs23 + DN_0_1*clhs25 + DN_0_2*clhs28 + clhs30);
            lhs(0,2)+=Weight*(DN_0_0*clhs31 + DN_0_1*clhs33 + DN_0_2*clhs35 + clhs36);
            lhs(0,3)+=-DN_0_0*clhs37;
            lhs(0,4)+=Weight*(DN_0_0*clhs38 + DN_0_1*clhs40 + DN_0_2*clhs42 + clhs45 + clhs50);
            lhs(0,5)+=Weight*(DN_0_0*clhs51 + DN_0_1*clhs53 + DN_0_2*clhs56 + clhs57);
            lhs(0,6)+=Weight*(DN_0_0*clhs58 + DN_0_1*clhs60 + DN_0_2*clhs62 + clhs63);
            lhs(0,7)+=-Weight*(-DN_1_0*clhs19 + DN_1_0*clhs21 + clhs64);
            lhs(0,8)+=Weight*(DN_0_0*clhs65 + DN_0_1*clhs67 + DN_0_2*clhs69 + clhs72 + clhs77);
            lhs(0,9)+=Weight*(DN_0_0*clhs78 + DN_0_1*clhs80 + DN_0_2*clhs83 + clhs84);
            lhs(0,10)+=Weight*(DN_0_0*clhs85 + DN_0_1*clhs87 + DN_0_2*clhs89 + clhs90);
            lhs(0,11)+=-Weight*(-DN_2_0*clhs19 + DN_2_0*clhs21 + clhs91);
            lhs(0,12)+=Weight*(DN_0_0*clhs92 + DN_0_1*clhs94 + DN_0_2*clhs96 + clhs104 + clhs99);
            lhs(0,13)+=Weight*(DN_0_0*clhs105 + DN_0_1*clhs107 + DN_0_2*clhs110 + clhs111);
            lhs(0,14)+=Weight*(DN_0_0*clhs112 + DN_0_1*clhs114 + DN_0_2*clhs116 + clhs117);
            lhs(0,15)+=-Weight*(-DN_3_0*clhs19 + DN_3_0*clhs21 + clhs118);
            lhs(1,0)+=Weight*(DN_0_0*clhs2 + DN_0_1*clhs119 + DN_0_2*clhs120 + clhs30);
            lhs(1,1)+=Weight*(DN_0_0*clhs25 + DN_0_1*clhs121 + DN_0_2*clhs123 + clhs10*clhs124 + clhs22);
            lhs(1,2)+=Weight*(DN_0_0*clhs33 + DN_0_1*clhs125 + DN_0_2*clhs127 + clhs129);
            lhs(1,3)+=-DN_0_1*clhs37;
            lhs(1,4)+=Weight*(DN_0_0*clhs40 + DN_0_1*clhs130 + DN_0_2*clhs131 + clhs132);
            lhs(1,5)+=Weight*(DN_0_0*clhs53 + DN_0_1*clhs133 + DN_0_2*clhs135 + clhs137 + clhs138);
            lhs(1,6)+=Weight*(DN_0_0*clhs60 + DN_0_1*clhs139 + DN_0_2*clhs141 + clhs142);
            lhs(1,7)+=-Weight*(-DN_1_1*clhs19 + DN_1_1*clhs21 + clhs143);
            lhs(1,8)+=Weight*(DN_0_0*clhs67 + DN_0_1*clhs144 + DN_0_2*clhs145 + clhs146);
            lhs(1,9)+=Weight*(DN_0_0*clhs80 + DN_0_1*clhs147 + DN_0_2*clhs149 + clhs151 + clhs152);
            lhs(1,10)+=Weight*(DN_0_0*clhs87 + DN_0_1*clhs153 + DN_0_2*clhs155 + clhs156);
            lhs(1,11)+=-Weight*(-DN_2_1*clhs19 + DN_2_1*clhs21 + clhs157);
            lhs(1,12)+=Weight*(DN_0_0*clhs94 + DN_0_1*clhs158 + DN_0_2*clhs159 + clhs160);
            lhs(1,13)+=Weight*(DN_0_0*clhs107 + DN_0_1*clhs161 + DN_0_2*clhs163 + clhs165 + clhs166);
            lhs(1,14)+=Weight*(DN_0_0*clhs114 + DN_0_1*clhs167 + DN_0_2*clhs169 + clhs170);
            lhs(1,15)+=-Weight*(-DN_3_1*clhs19 + DN_3_1*clhs21 + clhs171);
            lhs(2,0)+=Weight*(DN_0_0*clhs4 + DN_0_1*clhs120 + DN_0_2*clhs172 + clhs36);
            lhs(2,1)+=Weight*(DN_0_0*clhs28 + DN_0_1*clhs123 + DN_0_2*clhs173 + clhs129);
            lhs(2,2)+=Weight*(DN_0_0*clhs35 + DN_0_1*clhs127 + DN_0_2*clhs174 + clhs10*clhs175 + clhs22);
            lhs(2,3)+=-DN_0_2*clhs37;
            lhs(2,4)+=Weight*(DN_0_0*clhs42 + DN_0_1*clhs131 + DN_0_2*clhs176 + clhs178);
            lhs(2,5)+=Weight*(DN_0_0*clhs56 + DN_0_1*clhs135 + DN_0_2*clhs179 + clhs180);
            lhs(2,6)+=Weight*(DN_0_0*clhs62 + DN_0_1*clhs141 + DN_0_2*clhs181 + clhs138 + clhs183);
            lhs(2,7)+=-Weight*(-DN_1_2*clhs19 + DN_1_2*clhs21 + clhs184);
            lhs(2,8)+=Weight*(DN_0_0*clhs69 + DN_0_1*clhs145 + DN_0_2*clhs185 + clhs186);
            lhs(2,9)+=Weight*(DN_0_0*clhs83 + DN_0_1*clhs149 + DN_0_2*clhs187 + clhs188);
            lhs(2,10)+=Weight*(DN_0_0*clhs89 + DN_0_1*clhs155 + DN_0_2*clhs189 + clhs152 + clhs191);
            lhs(2,11)+=-Weight*(-DN_2_2*clhs19 + DN_2_2*clhs21 + clhs192);
            lhs(2,12)+=Weight*(DN_0_0*clhs96 + DN_0_1*clhs159 + DN_0_2*clhs193 + clhs194);
            lhs(2,13)+=Weight*(DN_0_0*clhs110 + DN_0_1*clhs163 + DN_0_2*clhs195 + clhs196);
            lhs(2,14)+=Weight*(DN_0_0*clhs116 + DN_0_1*clhs169 + DN_0_2*clhs197 + clhs166 + clhs199);
            lhs(2,15)+=-Weight*(-DN_3_2*clhs19 + DN_3_2*clhs21 + clhs200);
            lhs(3,0)+=DN_0_0*clhs201;
            lhs(3,1)+=DN_0_1*clhs201;
            lhs(3,2)+=DN_0_2*clhs201;
            lhs(3,3)+=clhs203*(clhs124 + clhs175 + clhs5);
            lhs(3,4)+=Weight*(DN_0_0*clhs205 + clhs204);
            lhs(3,5)+=Weight*(DN_0_1*clhs205 + clhs206);
            lhs(3,6)+=Weight*(DN_0_2*clhs205 + clhs207);
            lhs(3,7)+=clhs208;
            lhs(3,8)+=Weight*(DN_0_0*clhs210 + clhs209);
            lhs(3,9)+=Weight*(DN_0_1*clhs210 + clhs211);
            lhs(3,10)+=Weight*(DN_0_2*clhs210 + clhs212);
            lhs(3,11)+=clhs213;
            lhs(3,12)+=Weight*(DN_0_0*clhs215 + clhs214);
            lhs(3,13)+=Weight*(DN_0_1*clhs215 + clhs216);
            lhs(3,14)+=Weight*(DN_0_2*clhs215 + clhs217);
            lhs(3,15)+=clhs218;
            lhs(4,0)+=Weight*(DN_1_0*clhs0 + DN_1_1*clhs2 + DN_1_2*clhs4 + clhs223 + clhs45);
            lhs(4,1)+=Weight*(DN_1_0*clhs23 + DN_1_1*clhs25 + DN_1_2*clhs28 + clhs132);
            lhs(4,2)+=Weight*(DN_1_0*clhs31 + DN_1_1*clhs33 + DN_1_2*clhs35 + clhs178);
            lhs(4,3)+=-Weight*(-DN_0_0*clhs220 + DN_0_0*clhs222 + clhs204);
            lhs(4,4)+=Weight*(DN_1_0*clhs38 + DN_1_1*clhs40 + DN_1_2*clhs42 + clhs10*clhs224 + clhs226);
            lhs(4,5)+=Weight*(DN_1_0*clhs51 + DN_1_1*clhs53 + DN_1_2*clhs56 + clhs228);
            lhs(4,6)+=Weight*(DN_1_0*clhs58 + DN_1_1*clhs60 + DN_1_2*clhs62 + clhs229);
            lhs(4,7)+=-DN_1_0*clhs230;
            lhs(4,8)+=Weight*(DN_1_0*clhs65 + DN_1_1*clhs67 + DN_1_2*clhs69 + clhs233 + clhs234);
            lhs(4,9)+=Weight*(DN_1_0*clhs78 + DN_1_1*clhs80 + DN_1_2*clhs83 + clhs235);
            lhs(4,10)+=Weight*(DN_1_0*clhs85 + DN_1_1*clhs87 + DN_1_2*clhs89 + clhs236);
            lhs(4,11)+=-Weight*(-DN_2_0*clhs220 + DN_2_0*clhs222 + clhs237);
            lhs(4,12)+=Weight*(DN_1_0*clhs92 + DN_1_1*clhs94 + DN_1_2*clhs96 + clhs240 + clhs241);
            lhs(4,13)+=Weight*(DN_1_0*clhs105 + DN_1_1*clhs107 + DN_1_2*clhs110 + clhs242);
            lhs(4,14)+=Weight*(DN_1_0*clhs112 + DN_1_1*clhs114 + DN_1_2*clhs116 + clhs243);
            lhs(4,15)+=-Weight*(-DN_3_0*clhs220 + DN_3_0*clhs222 + clhs244);
            lhs(5,0)+=Weight*(DN_1_0*clhs2 + DN_1_1*clhs119 + DN_1_2*clhs120 + clhs57);
            lhs(5,1)+=Weight*(DN_1_0*clhs25 + DN_1_1*clhs121 + DN_1_2*clhs123 + clhs137 + clhs245);
            lhs(5,2)+=Weight*(DN_1_0*clhs33 + DN_1_1*clhs125 + DN_1_2*clhs127 + clhs180);
            lhs(5,3)+=-Weight*(-DN_0_1*clhs220 + DN_0_1*clhs222 + clhs206);
            lhs(5,4)+=Weight*(DN_1_0*clhs40 + DN_1_1*clhs130 + DN_1_2*clhs131 + clhs228);
            lhs(5,5)+=Weight*(DN_1_0*clhs53 + DN_1_1*clhs133 + DN_1_2*clhs135 + clhs10*clhs246 + clhs226);
            lhs(5,6)+=Weight*(DN_1_0*clhs60 + DN_1_1*clhs139 + DN_1_2*clhs141 + clhs248);
            lhs(5,7)+=-DN_1_1*clhs230;
            lhs(5,8)+=Weight*(DN_1_0*clhs67 + DN_1_1*clhs144 + DN_1_2*clhs145 + clhs249);
            lhs(5,9)+=Weight*(DN_1_0*clhs80 + DN_1_1*clhs147 + DN_1_2*clhs149 + clhs251 + clhs252);
            lhs(5,10)+=Weight*(DN_1_0*clhs87 + DN_1_1*clhs153 + DN_1_2*clhs155 + clhs253);
            lhs(5,11)+=-Weight*(-DN_2_1*clhs220 + DN_2_1*clhs222 + clhs254);
            lhs(5,12)+=Weight*(DN_1_0*clhs94 + DN_1_1*clhs158 + DN_1_2*clhs159 + clhs255);
            lhs(5,13)+=Weight*(DN_1_0*clhs107 + DN_1_1*clhs161 + DN_1_2*clhs163 + clhs257 + clhs258);
            lhs(5,14)+=Weight*(DN_1_0*clhs114 + DN_1_1*clhs167 + DN_1_2*clhs169 + clhs259);
            lhs(5,15)+=-Weight*(-DN_3_1*clhs220 + DN_3_1*clhs222 + clhs260);
            lhs(6,0)+=Weight*(DN_1_0*clhs4 + DN_1_1*clhs120 + DN_1_2*clhs172 + clhs63);
            lhs(6,1)+=Weight*(DN_1_0*clhs28 + DN_1_1*clhs123 + DN_1_2*clhs173 + clhs142);
            lhs(6,2)+=Weight*(DN_1_0*clhs35 + DN_1_1*clhs127 + DN_1_2*clhs174 + clhs183 + clhs245);
            lhs(6,3)+=-Weight*(-DN_0_2*clhs220 + DN_0_2*clhs222 + clhs207);
            lhs(6,4)+=Weight*(DN_1_0*clhs42 + DN_1_1*clhs131 + DN_1_2*clhs176 + clhs229);
            lhs(6,5)+=Weight*(DN_1_0*clhs56 + DN_1_1*clhs135 + DN_1_2*clhs179 + clhs248);
            lhs(6,6)+=Weight*(DN_1_0*clhs62 + DN_1_1*clhs141 + DN_1_2*clhs181 + clhs10*clhs261 + clhs226);
            lhs(6,7)+=-DN_1_2*clhs230;
            lhs(6,8)+=Weight*(DN_1_0*clhs69 + DN_1_1*clhs145 + DN_1_2*clhs185 + clhs263);
            lhs(6,9)+=Weight*(DN_1_0*clhs83 + DN_1_1*clhs149 + DN_1_2*clhs187 + clhs264);
            lhs(6,10)+=Weight*(DN_1_0*clhs89 + DN_1_1*clhs155 + DN_1_2*clhs189 + clhs252 + clhs266);
            lhs(6,11)+=-Weight*(-DN_2_2*clhs220 + DN_2_2*clhs222 + clhs267);
            lhs(6,12)+=Weight*(DN_1_0*clhs96 + DN_1_1*clhs159 + DN_1_2*clhs193 + clhs268);
            lhs(6,13)+=Weight*(DN_1_0*clhs110 + DN_1_1*clhs163 + DN_1_2*clhs195 + clhs269);
            lhs(6,14)+=Weight*(DN_1_0*clhs116 + DN_1_1*clhs169 + DN_1_2*clhs197 + clhs258 + clhs271);
            lhs(6,15)+=-Weight*(-DN_3_2*clhs220 + DN_3_2*clhs222 + clhs272);
            lhs(7,0)+=Weight*(DN_1_0*clhs273 + clhs64);
            lhs(7,1)+=Weight*(DN_1_1*clhs273 + clhs143);
            lhs(7,2)+=Weight*(DN_1_2*clhs273 + clhs184);
            lhs(7,3)+=clhs208;
            lhs(7,4)+=DN_1_0*clhs274;
            lhs(7,5)+=DN_1_1*clhs274;
            lhs(7,6)+=DN_1_2*clhs274;
            lhs(7,7)+=clhs203*(clhs224 + clhs246 + clhs261);
            lhs(7,8)+=Weight*(DN_1_0*clhs210 + clhs275);
            lhs(7,9)+=Weight*(DN_1_1*clhs210 + clhs276);
            lhs(7,10)+=Weight*(DN_1_2*clhs210 + clhs277);
            lhs(7,11)+=clhs278;
            lhs(7,12)+=Weight*(DN_1_0*clhs215 + clhs279);
            lhs(7,13)+=Weight*(DN_1_1*clhs215 + clhs280);
            lhs(7,14)+=Weight*(DN_1_2*clhs215 + clhs281);
            lhs(7,15)+=clhs282;
            lhs(8,0)+=Weight*(DN_2_0*clhs0 + DN_2_1*clhs2 + DN_2_2*clhs4 + clhs287 + clhs72);
            lhs(8,1)+=Weight*(DN_2_0*clhs23 + DN_2_1*clhs25 + DN_2_2*clhs28 + clhs146);
            lhs(8,2)+=Weight*(DN_2_0*clhs31 + DN_2_1*clhs33 + DN_2_2*clhs35 + clhs186);
            lhs(8,3)+=-Weight*(-DN_0_0*clhs284 + DN_0_0*clhs286 + clhs209);
            lhs(8,4)+=Weight*(DN_2_0*clhs38 + DN_2_1*clhs40 + DN_2_2*clhs42 + clhs233 + clhs288);
            lhs(8,5)+=Weight*(DN_2_0*clhs51 + DN_2_1*clhs53 + DN_2_2*clhs56 + clhs249);
            lhs(8,6)+=Weight*(DN_2_0*clhs58 + DN_2_1*clhs60 + DN_2_2*clhs62 + clhs263);
            lhs(8,7)+=-Weight*(-DN_1_0*clhs284 + DN_1_0*clhs286 + clhs275);
            lhs(8,8)+=Weight*(DN_2_0*clhs65 + DN_2_1*clhs67 + DN_2_2*clhs69 + clhs10*clhs289 + clhs291);
            lhs(8,9)+=Weight*(DN_2_0*clhs78 + DN_2_1*clhs80 + DN_2_2*clhs83 + clhs293);
            lhs(8,10)+=Weight*(DN_2_0*clhs85 + DN_2_1*clhs87 + DN_2_2*clhs89 + clhs294);
            lhs(8,11)+=-DN_2_0*clhs295;
            lhs(8,12)+=Weight*(DN_2_0*clhs92 + DN_2_1*clhs94 + DN_2_2*clhs96 + clhs298 + clhs299);
            lhs(8,13)+=Weight*(DN_2_0*clhs105 + DN_2_1*clhs107 + DN_2_2*clhs110 + clhs300);
            lhs(8,14)+=Weight*(DN_2_0*clhs112 + DN_2_1*clhs114 + DN_2_2*clhs116 + clhs301);
            lhs(8,15)+=-Weight*(-DN_3_0*clhs284 + DN_3_0*clhs286 + clhs302);
            lhs(9,0)+=Weight*(DN_2_0*clhs2 + DN_2_1*clhs119 + DN_2_2*clhs120 + clhs84);
            lhs(9,1)+=Weight*(DN_2_0*clhs25 + DN_2_1*clhs121 + DN_2_2*clhs123 + clhs151 + clhs303);
            lhs(9,2)+=Weight*(DN_2_0*clhs33 + DN_2_1*clhs125 + DN_2_2*clhs127 + clhs188);
            lhs(9,3)+=-Weight*(-DN_0_1*clhs284 + DN_0_1*clhs286 + clhs211);
            lhs(9,4)+=Weight*(DN_2_0*clhs40 + DN_2_1*clhs130 + DN_2_2*clhs131 + clhs235);
            lhs(9,5)+=Weight*(DN_2_0*clhs53 + DN_2_1*clhs133 + DN_2_2*clhs135 + clhs251 + clhs304);
            lhs(9,6)+=Weight*(DN_2_0*clhs60 + DN_2_1*clhs139 + DN_2_2*clhs141 + clhs264);
            lhs(9,7)+=-Weight*(-DN_1_1*clhs284 + DN_1_1*clhs286 + clhs276);
            lhs(9,8)+=Weight*(DN_2_0*clhs67 + DN_2_1*clhs144 + DN_2_2*clhs145 + clhs293);
            lhs(9,9)+=Weight*(DN_2_0*clhs80 + DN_2_1*clhs147 + DN_2_2*clhs149 + clhs10*clhs305 + clhs291);
            lhs(9,10)+=Weight*(DN_2_0*clhs87 + DN_2_1*clhs153 + DN_2_2*clhs155 + clhs307);
            lhs(9,11)+=-DN_2_1*clhs295;
            lhs(9,12)+=Weight*(DN_2_0*clhs94 + DN_2_1*clhs158 + DN_2_2*clhs159 + clhs308);
            lhs(9,13)+=Weight*(DN_2_0*clhs107 + DN_2_1*clhs161 + DN_2_2*clhs163 + clhs310 + clhs311);
            lhs(9,14)+=Weight*(DN_2_0*clhs114 + DN_2_1*clhs167 + DN_2_2*clhs169 + clhs312);
            lhs(9,15)+=-Weight*(-DN_3_1*clhs284 + DN_3_1*clhs286 + clhs313);
            lhs(10,0)+=Weight*(DN_2_0*clhs4 + DN_2_1*clhs120 + DN_2_2*clhs172 + clhs90);
            lhs(10,1)+=Weight*(DN_2_0*clhs28 + DN_2_1*clhs123 + DN_2_2*clhs173 + clhs156);
            lhs(10,2)+=Weight*(DN_2_0*clhs35 + DN_2_1*clhs127 + DN_2_2*clhs174 + clhs191 + clhs303);
            lhs(10,3)+=-Weight*(-DN_0_2*clhs284 + DN_0_2*clhs286 + clhs212);
            lhs(10,4)+=Weight*(DN_2_0*clhs42 + DN_2_1*clhs131 + DN_2_2*clhs176 + clhs236);
            lhs(10,5)+=Weight*(DN_2_0*clhs56 + DN_2_1*clhs135 + DN_2_2*clhs179 + clhs253);
            lhs(10,6)+=Weight*(DN_2_0*clhs62 + DN_2_1*clhs141 + DN_2_2*clhs181 + clhs266 + clhs304);
            lhs(10,7)+=-Weight*(-DN_1_2*clhs284 + DN_1_2*clhs286 + clhs277);
            lhs(10,8)+=Weight*(DN_2_0*clhs69 + DN_2_1*clhs145 + DN_2_2*clhs185 + clhs294);
            lhs(10,9)+=Weight*(DN_2_0*clhs83 + DN_2_1*clhs149 + DN_2_2*clhs187 + clhs307);
            lhs(10,10)+=Weight*(DN_2_0*clhs89 + DN_2_1*clhs155 + DN_2_2*clhs189 + clhs10*clhs314 + clhs291);
            lhs(10,11)+=-DN_2_2*clhs295;
            lhs(10,12)+=Weight*(DN_2_0*clhs96 + DN_2_1*clhs159 + DN_2_2*clhs193 + clhs316);
            lhs(10,13)+=Weight*(DN_2_0*clhs110 + DN_2_1*clhs163 + DN_2_2*clhs195 + clhs317);
            lhs(10,14)+=Weight*(DN_2_0*clhs116 + DN_2_1*clhs169 + DN_2_2*clhs197 + clhs311 + clhs319);
            lhs(10,15)+=-Weight*(-DN_3_2*clhs284 + DN_3_2*clhs286 + clhs320);
            lhs(11,0)+=Weight*(DN_2_0*clhs273 + clhs91);
            lhs(11,1)+=Weight*(DN_2_1*clhs273 + clhs157);
            lhs(11,2)+=Weight*(DN_2_2*clhs273 + clhs192);
            lhs(11,3)+=clhs213;
            lhs(11,4)+=Weight*(DN_2_0*clhs205 + clhs237);
            lhs(11,5)+=Weight*(DN_2_1*clhs205 + clhs254);
            lhs(11,6)+=Weight*(DN_2_2*clhs205 + clhs267);
            lhs(11,7)+=clhs278;
            lhs(11,8)+=DN_2_0*clhs321;
            lhs(11,9)+=DN_2_1*clhs321;
            lhs(11,10)+=DN_2_2*clhs321;
            lhs(11,11)+=clhs203*(clhs289 + clhs305 + clhs314);
            lhs(11,12)+=Weight*(DN_2_0*clhs215 + clhs322);
            lhs(11,13)+=Weight*(DN_2_1*clhs215 + clhs323);
            lhs(11,14)+=Weight*(DN_2_2*clhs215 + clhs324);
            lhs(11,15)+=clhs325;
            lhs(12,0)+=Weight*(DN_3_0*clhs0 + DN_3_1*clhs2 + DN_3_2*clhs4 + clhs330 + clhs99);
            lhs(12,1)+=Weight*(DN_3_0*clhs23 + DN_3_1*clhs25 + DN_3_2*clhs28 + clhs160);
            lhs(12,2)+=Weight*(DN_3_0*clhs31 + DN_3_1*clhs33 + DN_3_2*clhs35 + clhs194);
            lhs(12,3)+=-Weight*(-DN_0_0*clhs327 + DN_0_0*clhs329 + clhs214);
            lhs(12,4)+=Weight*(DN_3_0*clhs38 + DN_3_1*clhs40 + DN_3_2*clhs42 + clhs240 + clhs331);
            lhs(12,5)+=Weight*(DN_3_0*clhs51 + DN_3_1*clhs53 + DN_3_2*clhs56 + clhs255);
            lhs(12,6)+=Weight*(DN_3_0*clhs58 + DN_3_1*clhs60 + DN_3_2*clhs62 + clhs268);
            lhs(12,7)+=-Weight*(-DN_1_0*clhs327 + DN_1_0*clhs329 + clhs279);
            lhs(12,8)+=Weight*(DN_3_0*clhs65 + DN_3_1*clhs67 + DN_3_2*clhs69 + clhs298 + clhs332);
            lhs(12,9)+=Weight*(DN_3_0*clhs78 + DN_3_1*clhs80 + DN_3_2*clhs83 + clhs308);
            lhs(12,10)+=Weight*(DN_3_0*clhs85 + DN_3_1*clhs87 + DN_3_2*clhs89 + clhs316);
            lhs(12,11)+=-Weight*(-DN_2_0*clhs327 + DN_2_0*clhs329 + clhs322);
            lhs(12,12)+=Weight*(DN_3_0*clhs92 + DN_3_1*clhs94 + DN_3_2*clhs96 + clhs10*clhs333 + clhs335);
            lhs(12,13)+=Weight*(DN_3_0*clhs105 + DN_3_1*clhs107 + DN_3_2*clhs110 + clhs337);
            lhs(12,14)+=Weight*(DN_3_0*clhs112 + DN_3_1*clhs114 + DN_3_2*clhs116 + clhs338);
            lhs(12,15)+=-DN_3_0*clhs339;
            lhs(13,0)+=Weight*(DN_3_0*clhs2 + DN_3_1*clhs119 + DN_3_2*clhs120 + clhs111);
            lhs(13,1)+=Weight*(DN_3_0*clhs25 + DN_3_1*clhs121 + DN_3_2*clhs123 + clhs165 + clhs340);
            lhs(13,2)+=Weight*(DN_3_0*clhs33 + DN_3_1*clhs125 + DN_3_2*clhs127 + clhs196);
            lhs(13,3)+=-Weight*(-DN_0_1*clhs327 + DN_0_1*clhs329 + clhs216);
            lhs(13,4)+=Weight*(DN_3_0*clhs40 + DN_3_1*clhs130 + DN_3_2*clhs131 + clhs242);
            lhs(13,5)+=Weight*(DN_3_0*clhs53 + DN_3_1*clhs133 + DN_3_2*clhs135 + clhs257 + clhs341);
            lhs(13,6)+=Weight*(DN_3_0*clhs60 + DN_3_1*clhs139 + DN_3_2*clhs141 + clhs269);
            lhs(13,7)+=-Weight*(-DN_1_1*clhs327 + DN_1_1*clhs329 + clhs280);
            lhs(13,8)+=Weight*(DN_3_0*clhs67 + DN_3_1*clhs144 + DN_3_2*clhs145 + clhs300);
            lhs(13,9)+=Weight*(DN_3_0*clhs80 + DN_3_1*clhs147 + DN_3_2*clhs149 + clhs310 + clhs342);
            lhs(13,10)+=Weight*(DN_3_0*clhs87 + DN_3_1*clhs153 + DN_3_2*clhs155 + clhs317);
            lhs(13,11)+=-Weight*(-DN_2_1*clhs327 + DN_2_1*clhs329 + clhs323);
            lhs(13,12)+=Weight*(DN_3_0*clhs94 + DN_3_1*clhs158 + DN_3_2*clhs159 + clhs337);
            lhs(13,13)+=Weight*(DN_3_0*clhs107 + DN_3_1*clhs161 + DN_3_2*clhs163 + clhs10*clhs343 + clhs335);
            lhs(13,14)+=Weight*(DN_3_0*clhs114 + DN_3_1*clhs167 + DN_3_2*clhs169 + clhs344);
            lhs(13,15)+=-DN_3_1*clhs339;
            lhs(14,0)+=Weight*(DN_3_0*clhs4 + DN_3_1*clhs120 + DN_3_2*clhs172 + clhs117);
            lhs(14,1)+=Weight*(DN_3_0*clhs28 + DN_3_1*clhs123 + DN_3_2*clhs173 + clhs170);
            lhs(14,2)+=Weight*(DN_3_0*clhs35 + DN_3_1*clhs127 + DN_3_2*clhs174 + clhs199 + clhs340);
            lhs(14,3)+=-Weight*(-DN_0_2*clhs327 + DN_0_2*clhs329 + clhs217);
            lhs(14,4)+=Weight*(DN_3_0*clhs42 + DN_3_1*clhs131 + DN_3_2*clhs176 + clhs243);
            lhs(14,5)+=Weight*(DN_3_0*clhs56 + DN_3_1*clhs135 + DN_3_2*clhs179 + clhs259);
            lhs(14,6)+=Weight*(DN_3_0*clhs62 + DN_3_1*clhs141 + DN_3_2*clhs181 + clhs271 + clhs341);
            lhs(14,7)+=-Weight*(-DN_1_2*clhs327 + DN_1_2*clhs329 + clhs281);
            lhs(14,8)+=Weight*(DN_3_0*clhs69 + DN_3_1*clhs145 + DN_3_2*clhs185 + clhs301);
            lhs(14,9)+=Weight*(DN_3_0*clhs83 + DN_3_1*clhs149 + DN_3_2*clhs187 + clhs312);
            lhs(14,10)+=Weight*(DN_3_0*clhs89 + DN_3_1*clhs155 + DN_3_2*clhs189 + clhs319 + clhs342);
            lhs(14,11)+=-Weight*(-DN_2_2*clhs327 + DN_2_2*clhs329 + clhs324);
            lhs(14,12)+=Weight*(DN_3_0*clhs96 + DN_3_1*clhs159 + DN_3_2*clhs193 + clhs338);
            lhs(14,13)+=Weight*(DN_3_0*clhs110 + DN_3_1*clhs163 + DN_3_2*clhs195 + clhs344);
            lhs(14,14)+=Weight*(DN_3_0*clhs116 + DN_3_1*clhs169 + DN_3_2*clhs197 + clhs10*clhs345 + clhs335);
            lhs(14,15)+=-DN_3_2*clhs339;
            lhs(15,0)+=Weight*(DN_3_0*clhs273 + clhs118);
            lhs(15,1)+=Weight*(DN_3_1*clhs273 + clhs171);
            lhs(15,2)+=Weight*(DN_3_2*clhs273 + clhs200);
            lhs(15,3)+=clhs218;
            lhs(15,4)+=Weight*(DN_3_0*clhs205 + clhs244);
            lhs(15,5)+=Weight*(DN_3_1*clhs205 + clhs260);
            lhs(15,6)+=Weight*(DN_3_2*clhs205 + clhs272);
            lhs(15,7)+=clhs282;
            lhs(15,8)+=Weight*(DN_3_0*clhs210 + clhs302);
            lhs(15,9)+=Weight*(DN_3_1*clhs210 + clhs313);
            lhs(15,10)+=Weight*(DN_3_2*clhs210 + clhs320);
            lhs(15,11)+=clhs325;
            lhs(15,12)+=DN_3_0*clhs346;
            lhs(15,13)+=DN_3_1*clhs346;
            lhs(15,14)+=DN_3_2*clhs346;
            lhs(15,15)+=clhs203*(clhs333 + clhs343 + clhs345);


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

    const double Weight = rData.Weight;
    auto &rhs = rRHS;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs2 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs7 =             rho*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs8 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs9 =             crhs4 + crhs8;
const double crhs10 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
const double crhs11 =             crhs9*(crhs10*h/stab_c1 + mu);
const double crhs12 =             1.0/(K_darcy + crhs10/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs13 =             crhs12*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs1 - crhs2 + crhs3 + crhs7);
const double crhs14 =             K_darcy*N[0];
const double crhs15 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
const double crhs16 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs17 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs18 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs19 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs8);
const double crhs20 =             crhs12*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs16 - crhs17 + crhs18 + crhs19);
const double crhs21 =             K_darcy*N[1];
const double crhs22 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
const double crhs23 =             K_darcy*N[2];
const double crhs24 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
            rhs[0]+=-Weight*(-DN(0,0)*crhs0 + DN(0,0)*crhs11 + DN(0,0)*stress[0] + DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 + N[0]*crhs3 + N[0]*crhs7 - crhs13*crhs14 + crhs13*crhs15);
            rhs[1]+=-Weight*(DN(0,0)*stress[2] - DN(0,1)*crhs0 + DN(0,1)*crhs11 + DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 + N[0]*crhs18 + N[0]*crhs19 - crhs14*crhs20 + crhs15*crhs20);
            rhs[2]+=-Weight*(DN(0,0)*crhs13 + DN(0,1)*crhs20 + N[0]*crhs9);
            rhs[3]+=-Weight*(-DN(1,0)*crhs0 + DN(1,0)*crhs11 + DN(1,0)*stress[0] + DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 + N[1]*crhs3 + N[1]*crhs7 - crhs13*crhs21 + crhs13*crhs22);
            rhs[4]+=-Weight*(DN(1,0)*stress[2] - DN(1,1)*crhs0 + DN(1,1)*crhs11 + DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 + N[1]*crhs18 + N[1]*crhs19 - crhs20*crhs21 + crhs20*crhs22);
            rhs[5]+=-Weight*(DN(1,0)*crhs13 + DN(1,1)*crhs20 + N[1]*crhs9);
            rhs[6]+=-Weight*(-DN(2,0)*crhs0 + DN(2,0)*crhs11 + DN(2,0)*stress[0] + DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 + N[2]*crhs3 + N[2]*crhs7 - crhs13*crhs23 + crhs13*crhs24);
            rhs[7]+=-Weight*(DN(2,0)*stress[2] - DN(2,1)*crhs0 + DN(2,1)*crhs11 + DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 + N[2]*crhs18 + N[2]*crhs19 - crhs20*crhs23 + crhs20*crhs24);
            rhs[8]+=-Weight*(DN(2,0)*crhs13 + DN(2,1)*crhs20 + N[2]*crhs9);


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

    const double Weight = rData.Weight;
    auto &rhs = rRHS;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crhs2 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs7 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs8 =             rho*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs7*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs9 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs10 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs11 =             crhs10 + crhs4 + crhs9;
const double crhs12 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2) + pow(crhs7, 2));
const double crhs13 =             crhs11*(crhs12*h/stab_c1 + mu);
const double crhs14 =             1.0/(K_darcy + crhs12/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs15 =             crhs14*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs1 - crhs2 + crhs3 + crhs8);
const double crhs16 =             K_darcy*N[0];
const double crhs17 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6 + DN(0,2)*crhs7);
const double crhs18 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs19 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs20 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs21 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*crhs9 + crhs7*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs22 =             crhs14*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs18 - crhs19 + crhs20 + crhs21);
const double crhs23 =             K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs24 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs25 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs26 =             rho*(crhs10*crhs7 + crhs5*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs6*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs27 =             crhs14*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + crhs23 - crhs24 + crhs25 + crhs26);
const double crhs28 =             K_darcy*N[1];
const double crhs29 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6 + DN(1,2)*crhs7);
const double crhs30 =             K_darcy*N[2];
const double crhs31 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6 + DN(2,2)*crhs7);
const double crhs32 =             K_darcy*N[3];
const double crhs33 =             rho*(DN(3,0)*crhs5 + DN(3,1)*crhs6 + DN(3,2)*crhs7);
            rhs[0]+=-Weight*(-DN(0,0)*crhs0 + DN(0,0)*crhs13 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs2 + N[0]*crhs3 + N[0]*crhs8 - crhs15*crhs16 + crhs15*crhs17);
            rhs[1]+=-Weight*(DN(0,0)*stress[3] - DN(0,1)*crhs0 + DN(0,1)*crhs13 + DN(0,1)*stress[1] + DN(0,2)*stress[4] + N[0]*crhs18 - N[0]*crhs19 + N[0]*crhs20 + N[0]*crhs21 - crhs16*crhs22 + crhs17*crhs22);
            rhs[2]+=-Weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crhs0 + DN(0,2)*crhs13 + DN(0,2)*stress[2] + N[0]*crhs23 - N[0]*crhs24 + N[0]*crhs25 + N[0]*crhs26 - crhs16*crhs27 + crhs17*crhs27);
            rhs[3]+=-Weight*(DN(0,0)*crhs15 + DN(0,1)*crhs22 + DN(0,2)*crhs27 + N[0]*crhs11);
            rhs[4]+=-Weight*(-DN(1,0)*crhs0 + DN(1,0)*crhs13 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs2 + N[1]*crhs3 + N[1]*crhs8 - crhs15*crhs28 + crhs15*crhs29);
            rhs[5]+=-Weight*(DN(1,0)*stress[3] - DN(1,1)*crhs0 + DN(1,1)*crhs13 + DN(1,1)*stress[1] + DN(1,2)*stress[4] + N[1]*crhs18 - N[1]*crhs19 + N[1]*crhs20 + N[1]*crhs21 - crhs22*crhs28 + crhs22*crhs29);
            rhs[6]+=-Weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crhs0 + DN(1,2)*crhs13 + DN(1,2)*stress[2] + N[1]*crhs23 - N[1]*crhs24 + N[1]*crhs25 + N[1]*crhs26 - crhs27*crhs28 + crhs27*crhs29);
            rhs[7]+=-Weight*(DN(1,0)*crhs15 + DN(1,1)*crhs22 + DN(1,2)*crhs27 + N[1]*crhs11);
            rhs[8]+=-Weight*(-DN(2,0)*crhs0 + DN(2,0)*crhs13 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs2 + N[2]*crhs3 + N[2]*crhs8 - crhs15*crhs30 + crhs15*crhs31);
            rhs[9]+=-Weight*(DN(2,0)*stress[3] - DN(2,1)*crhs0 + DN(2,1)*crhs13 + DN(2,1)*stress[1] + DN(2,2)*stress[4] + N[2]*crhs18 - N[2]*crhs19 + N[2]*crhs20 + N[2]*crhs21 - crhs22*crhs30 + crhs22*crhs31);
            rhs[10]+=-Weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crhs0 + DN(2,2)*crhs13 + DN(2,2)*stress[2] + N[2]*crhs23 - N[2]*crhs24 + N[2]*crhs25 + N[2]*crhs26 - crhs27*crhs30 + crhs27*crhs31);
            rhs[11]+=-Weight*(DN(2,0)*crhs15 + DN(2,1)*crhs22 + DN(2,2)*crhs27 + N[2]*crhs11);
            rhs[12]+=-Weight*(-DN(3,0)*crhs0 + DN(3,0)*crhs13 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs2 + N[3]*crhs3 + N[3]*crhs8 - crhs15*crhs32 + crhs15*crhs33);
            rhs[13]+=-Weight*(DN(3,0)*stress[3] - DN(3,1)*crhs0 + DN(3,1)*crhs13 + DN(3,1)*stress[1] + DN(3,2)*stress[4] + N[3]*crhs18 - N[3]*crhs19 + N[3]*crhs20 + N[3]*crhs21 - crhs22*crhs32 + crhs22*crhs33);
            rhs[14]+=-Weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crhs0 + DN(3,2)*crhs13 + DN(3,2)*stress[2] + N[3]*crhs23 - N[3]*crhs24 + N[3]*crhs25 + N[3]*crhs26 - crhs27*crhs32 + crhs27*crhs33);
            rhs[15]+=-Weight*(DN(3,0)*crhs15 + DN(3,1)*crhs22 + DN(3,2)*crhs27 + N[3]*crhs11);


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
const double cV3 =             DNenr(0,0)*cV2;
const double cV4 =             K_darcy*N[0];
const double cV5 =             rho*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV6 =             DNenr(1,0)*cV2;
const double cV7 =             DNenr(2,0)*cV2;
const double cV8 =             DNenr(0,1)*cV2;
const double cV9 =             DNenr(1,1)*cV2;
const double cV10 =             DNenr(2,1)*cV2;
const double cV11 =             K_darcy*N[1];
const double cV12 =             rho*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV13 =             K_darcy*N[2];
const double cV14 =             rho*(DN(2,0)*cV0 + DN(2,1)*cV1);
            V(0,0)=-DN(0,0)*Nenr[0] - cV3*cV4 + cV3*cV5;
            V(0,1)=-DN(0,0)*Nenr[1] - cV4*cV6 + cV5*cV6;
            V(0,2)=-DN(0,0)*Nenr[2] - cV4*cV7 + cV5*cV7;
            V(1,0)=-DN(0,1)*Nenr[0] - cV4*cV8 + cV5*cV8;
            V(1,1)=-DN(0,1)*Nenr[1] - cV4*cV9 + cV5*cV9;
            V(1,2)=-DN(0,1)*Nenr[2] - cV10*cV4 + cV10*cV5;
            V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] - cV11*cV3 + cV12*cV3;
            V(3,1)=-DN(1,0)*Nenr[1] - cV11*cV6 + cV12*cV6;
            V(3,2)=-DN(1,0)*Nenr[2] - cV11*cV7 + cV12*cV7;
            V(4,0)=-DN(1,1)*Nenr[0] - cV11*cV8 + cV12*cV8;
            V(4,1)=-DN(1,1)*Nenr[1] - cV11*cV9 + cV12*cV9;
            V(4,2)=-DN(1,1)*Nenr[2] - cV10*cV11 + cV10*cV12;
            V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] - cV13*cV3 + cV14*cV3;
            V(6,1)=-DN(2,0)*Nenr[1] - cV13*cV6 + cV14*cV6;
            V(6,2)=-DN(2,0)*Nenr[2] - cV13*cV7 + cV14*cV7;
            V(7,0)=-DN(2,1)*Nenr[0] - cV13*cV8 + cV14*cV8;
            V(7,1)=-DN(2,1)*Nenr[1] - cV13*cV9 + cV14*cV9;
            V(7,2)=-DN(2,1)*Nenr[2] - cV10*cV13 + cV10*cV14;
            V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH3 =             cH2*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0));
const double cH4 =             cH2*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0));
const double cH5 =             cH2*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0));
            H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH3;
            H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH3;
            H(0,2)=cH2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH4;
            H(0,4)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH4;
            H(0,5)=cH2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH5;
            H(0,7)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH5;
            H(0,8)=cH2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH3;
            H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH3;
            H(1,2)=cH2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH4;
            H(1,4)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH4;
            H(1,5)=cH2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH5;
            H(1,7)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH5;
            H(1,8)=cH2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH3;
            H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH3;
            H(2,2)=cH2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH4;
            H(2,4)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH4;
            H(2,5)=cH2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH5;
            H(2,7)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH5;
            H(2,8)=cH2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


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
const double crhs_ee6 =             crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 =             crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
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
const double cV4 =             DNenr(0,0)*cV3;
const double cV5 =             K_darcy*N[0];
const double cV6 =             rho*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV7 =             DNenr(1,0)*cV3;
const double cV8 =             DNenr(2,0)*cV3;
const double cV9 =             DNenr(3,0)*cV3;
const double cV10 =             DNenr(0,1)*cV3;
const double cV11 =             DNenr(1,1)*cV3;
const double cV12 =             DNenr(2,1)*cV3;
const double cV13 =             DNenr(3,1)*cV3;
const double cV14 =             DNenr(0,2)*cV3;
const double cV15 =             DNenr(1,2)*cV3;
const double cV16 =             DNenr(2,2)*cV3;
const double cV17 =             DNenr(3,2)*cV3;
const double cV18 =             K_darcy*N[1];
const double cV19 =             rho*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV20 =             K_darcy*N[2];
const double cV21 =             rho*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV22 =             K_darcy*N[3];
const double cV23 =             rho*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
            V(0,0)=-DN(0,0)*Nenr[0] - cV4*cV5 + cV4*cV6;
            V(0,1)=-DN(0,0)*Nenr[1] - cV5*cV7 + cV6*cV7;
            V(0,2)=-DN(0,0)*Nenr[2] - cV5*cV8 + cV6*cV8;
            V(0,3)=-DN(0,0)*Nenr[3] - cV5*cV9 + cV6*cV9;
            V(1,0)=-DN(0,1)*Nenr[0] - cV10*cV5 + cV10*cV6;
            V(1,1)=-DN(0,1)*Nenr[1] - cV11*cV5 + cV11*cV6;
            V(1,2)=-DN(0,1)*Nenr[2] - cV12*cV5 + cV12*cV6;
            V(1,3)=-DN(0,1)*Nenr[3] - cV13*cV5 + cV13*cV6;
            V(2,0)=-DN(0,2)*Nenr[0] - cV14*cV5 + cV14*cV6;
            V(2,1)=-DN(0,2)*Nenr[1] - cV15*cV5 + cV15*cV6;
            V(2,2)=-DN(0,2)*Nenr[2] - cV16*cV5 + cV16*cV6;
            V(2,3)=-DN(0,2)*Nenr[3] - cV17*cV5 + cV17*cV6;
            V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] - cV18*cV4 + cV19*cV4;
            V(4,1)=-DN(1,0)*Nenr[1] - cV18*cV7 + cV19*cV7;
            V(4,2)=-DN(1,0)*Nenr[2] - cV18*cV8 + cV19*cV8;
            V(4,3)=-DN(1,0)*Nenr[3] - cV18*cV9 + cV19*cV9;
            V(5,0)=-DN(1,1)*Nenr[0] - cV10*cV18 + cV10*cV19;
            V(5,1)=-DN(1,1)*Nenr[1] - cV11*cV18 + cV11*cV19;
            V(5,2)=-DN(1,1)*Nenr[2] - cV12*cV18 + cV12*cV19;
            V(5,3)=-DN(1,1)*Nenr[3] - cV13*cV18 + cV13*cV19;
            V(6,0)=-DN(1,2)*Nenr[0] - cV14*cV18 + cV14*cV19;
            V(6,1)=-DN(1,2)*Nenr[1] - cV15*cV18 + cV15*cV19;
            V(6,2)=-DN(1,2)*Nenr[2] - cV16*cV18 + cV16*cV19;
            V(6,3)=-DN(1,2)*Nenr[3] - cV17*cV18 + cV17*cV19;
            V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] - cV20*cV4 + cV21*cV4;
            V(8,1)=-DN(2,0)*Nenr[1] - cV20*cV7 + cV21*cV7;
            V(8,2)=-DN(2,0)*Nenr[2] - cV20*cV8 + cV21*cV8;
            V(8,3)=-DN(2,0)*Nenr[3] - cV20*cV9 + cV21*cV9;
            V(9,0)=-DN(2,1)*Nenr[0] - cV10*cV20 + cV10*cV21;
            V(9,1)=-DN(2,1)*Nenr[1] - cV11*cV20 + cV11*cV21;
            V(9,2)=-DN(2,1)*Nenr[2] - cV12*cV20 + cV12*cV21;
            V(9,3)=-DN(2,1)*Nenr[3] - cV13*cV20 + cV13*cV21;
            V(10,0)=-DN(2,2)*Nenr[0] - cV14*cV20 + cV14*cV21;
            V(10,1)=-DN(2,2)*Nenr[1] - cV15*cV20 + cV15*cV21;
            V(10,2)=-DN(2,2)*Nenr[2] - cV16*cV20 + cV16*cV21;
            V(10,3)=-DN(2,2)*Nenr[3] - cV17*cV20 + cV17*cV21;
            V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] - cV22*cV4 + cV23*cV4;
            V(12,1)=-DN(3,0)*Nenr[1] - cV22*cV7 + cV23*cV7;
            V(12,2)=-DN(3,0)*Nenr[2] - cV22*cV8 + cV23*cV8;
            V(12,3)=-DN(3,0)*Nenr[3] - cV22*cV9 + cV23*cV9;
            V(13,0)=-DN(3,1)*Nenr[0] - cV10*cV22 + cV10*cV23;
            V(13,1)=-DN(3,1)*Nenr[1] - cV11*cV22 + cV11*cV23;
            V(13,2)=-DN(3,1)*Nenr[2] - cV12*cV22 + cV12*cV23;
            V(13,3)=-DN(3,1)*Nenr[3] - cV13*cV22 + cV13*cV23;
            V(14,0)=-DN(3,2)*Nenr[0] - cV14*cV22 + cV14*cV23;
            V(14,1)=-DN(3,2)*Nenr[1] - cV15*cV22 + cV15*cV23;
            V(14,2)=-DN(3,2)*Nenr[2] - cV16*cV22 + cV16*cV23;
            V(14,3)=-DN(3,2)*Nenr[3] - cV17*cV22 + cV17*cV23;
            V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 =             cH3*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0));
const double cH5 =             cH3*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0));
const double cH6 =             cH3*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0));
const double cH7 =             cH3*(K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0));
            H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH4;
            H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH4;
            H(0,2)=DN(0,2)*Nenr[0] + DNenr(0,2)*cH4;
            H(0,3)=cH3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH5;
            H(0,5)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH5;
            H(0,6)=DN(1,2)*Nenr[0] + DNenr(0,2)*cH5;
            H(0,7)=cH3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH6;
            H(0,9)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH6;
            H(0,10)=DN(2,2)*Nenr[0] + DNenr(0,2)*cH6;
            H(0,11)=cH3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DN(3,0)*Nenr[0] + DNenr(0,0)*cH7;
            H(0,13)=DN(3,1)*Nenr[0] + DNenr(0,1)*cH7;
            H(0,14)=DN(3,2)*Nenr[0] + DNenr(0,2)*cH7;
            H(0,15)=cH3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH4;
            H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH4;
            H(1,2)=DN(0,2)*Nenr[1] + DNenr(1,2)*cH4;
            H(1,3)=cH3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH5;
            H(1,5)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH5;
            H(1,6)=DN(1,2)*Nenr[1] + DNenr(1,2)*cH5;
            H(1,7)=cH3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH6;
            H(1,9)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH6;
            H(1,10)=DN(2,2)*Nenr[1] + DNenr(1,2)*cH6;
            H(1,11)=cH3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DN(3,0)*Nenr[1] + DNenr(1,0)*cH7;
            H(1,13)=DN(3,1)*Nenr[1] + DNenr(1,1)*cH7;
            H(1,14)=DN(3,2)*Nenr[1] + DNenr(1,2)*cH7;
            H(1,15)=cH3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH4;
            H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH4;
            H(2,2)=DN(0,2)*Nenr[2] + DNenr(2,2)*cH4;
            H(2,3)=cH3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH5;
            H(2,5)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH5;
            H(2,6)=DN(1,2)*Nenr[2] + DNenr(2,2)*cH5;
            H(2,7)=cH3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH6;
            H(2,9)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH6;
            H(2,10)=DN(2,2)*Nenr[2] + DNenr(2,2)*cH6;
            H(2,11)=cH3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DN(3,0)*Nenr[2] + DNenr(2,0)*cH7;
            H(2,13)=DN(3,1)*Nenr[2] + DNenr(2,1)*cH7;
            H(2,14)=DN(3,2)*Nenr[2] + DNenr(2,2)*cH7;
            H(2,15)=cH3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DN(0,0)*Nenr[3] + DNenr(3,0)*cH4;
            H(3,1)=DN(0,1)*Nenr[3] + DNenr(3,1)*cH4;
            H(3,2)=DN(0,2)*Nenr[3] + DNenr(3,2)*cH4;
            H(3,3)=cH3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DN(1,0)*Nenr[3] + DNenr(3,0)*cH5;
            H(3,5)=DN(1,1)*Nenr[3] + DNenr(3,1)*cH5;
            H(3,6)=DN(1,2)*Nenr[3] + DNenr(3,2)*cH5;
            H(3,7)=cH3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DN(2,0)*Nenr[3] + DNenr(3,0)*cH6;
            H(3,9)=DN(2,1)*Nenr[3] + DNenr(3,1)*cH6;
            H(3,10)=DN(2,2)*Nenr[3] + DNenr(3,2)*cH6;
            H(3,11)=cH3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DN(3,0)*Nenr[3] + DNenr(3,0)*cH7;
            H(3,13)=DN(3,1)*Nenr[3] + DNenr(3,1)*cH7;
            H(3,14)=DN(3,2)*Nenr[3] + DNenr(3,2)*cH7;
            H(3,15)=cH3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


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


    const double crhs_ee0 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs_ee1 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs_ee2 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee3 =             crhs_ee0 + crhs_ee1 + crhs_ee2;
const double crhs_ee4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee6 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee7 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee4, 2) + pow(crhs_ee5, 2) + pow(crhs_ee6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee8 =             crhs_ee7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*crhs_ee4 + crhs_ee5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee9 =             crhs_ee7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee1*crhs_ee5 + crhs_ee4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee10 =             crhs_ee7*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee2*crhs_ee6 + crhs_ee4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee8 - DNenr(0,1)*crhs_ee9 - DNenr(0,2)*crhs_ee10 - Nenr[0]*crhs_ee3;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee8 - DNenr(1,1)*crhs_ee9 - DNenr(1,2)*crhs_ee10 - Nenr[1]*crhs_ee3;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee8 - DNenr(2,1)*crhs_ee9 - DNenr(2,2)*crhs_ee10 - Nenr[2]*crhs_ee3;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee8 - DNenr(3,1)*crhs_ee9 - DNenr(3,2)*crhs_ee10 - Nenr[3]*crhs_ee3;


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
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
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

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the positive side modified shape functions calculator
    pModifiedShapeFunctions->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideShapeFunctionsAndGradientsValues(
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

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); ++i){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); ++i){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (pModifiedShapeFunctions->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ComputeSplitInterface(
    const TElementData &rData,
    MatrixType& rInterfaceShapeFunctionNeg,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Vector& rInterfaceWeightsNeg,
    std::vector<Vector>& rInterfaceNormalsNeg,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions)
{
    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Call the Interface negative side shape functions calculator
    pModifiedShapeFunctions->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivativesNeg,
        rInterfaceWeightsNeg,
        GeometryData::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideInterfaceAreaNormals(
        rInterfaceNormalsNeg,
        GeometryData::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); ++gp){
        const double normal_norm = norm_2(rInterfaceNormalsNeg[gp]);
        rInterfaceNormalsNeg[gp] /= normal_norm;
    }

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    auto p_modified_sh_func = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
    return p_modified_sh_func;
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    auto p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
    return p_modified_sh_func;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateCurvatureOnInterfaceGaussPoints(
        const Matrix& rInterfaceShapeFunctions,
        Vector& rInterfaceCurvature)
{
    const auto& r_geom = this->GetGeometry();
    const unsigned int n_gpt = rInterfaceShapeFunctions.size1();

    rInterfaceCurvature.resize(n_gpt, false);

    for (unsigned int gpt = 0; gpt < n_gpt; ++gpt){
        double curvature = 0.0;
        for (unsigned int i = 0; i < NumNodes; ++i){
            curvature += rInterfaceShapeFunctions(gpt,i) * r_geom[i].GetValue(CURVATURE);
        }
        rInterfaceCurvature[gpt] = curvature;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
    const double SurfaceTensionCoefficient,
    const Vector& rCurvature,
    const Vector& rInterfaceWeights,
    const Matrix& rInterfaceShapeFunctions,
    const std::vector<Vector>& rInterfaceNormalsNeg,
    VectorType& rRHS)
{
    for (unsigned int intgp = 0; intgp < rInterfaceWeights.size(); ++intgp){
        const double intgp_curv = rCurvature(intgp);
        const double intgp_w = rInterfaceWeights(intgp);
        const auto& intgp_normal = rInterfaceNormalsNeg[intgp];
        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < NumNodes-1; ++dim){
                rRHS[ j*(NumNodes) + dim ] -= SurfaceTensionCoefficient*intgp_normal[dim]
                    *intgp_curv*intgp_w*rInterfaceShapeFunctions(intgp,j);
            }
        }
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureGradientStabilization(
    const TElementData& rData,
    const Vector& rInterfaceWeights,
    const Matrix& rEnrInterfaceShapeFunctionPos,
    const Matrix& rEnrInterfaceShapeFunctionNeg,
    const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivatives,
    MatrixType& rKeeTot,
    VectorType& rRHSeeTot)
{
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);
    VectorType rhs_enr = ZeroVector(NumNodes);

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    double positive_density = 0.0;
    double negative_density = 0.0;
    double positive_viscosity = 0.0;
    double negative_viscosity = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
            positive_density = rData.NodalDensity[i];
            positive_viscosity = rData.NodalDynamicViscosity[i];
        } else{
            enr_pos_interp(i, i) = 1.0;
            negative_density = rData.NodalDensity[i];
            negative_viscosity = rData.NodalDynamicViscosity[i];
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesPos = rInterfaceShapeDerivatives;
    GeometryType::ShapeFunctionsGradientsType EnrichedInterfaceShapeDerivativesNeg = rInterfaceShapeDerivatives;

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesPos[i] = prod(enr_pos_interp, rInterfaceShapeDerivatives[i]);
    }

    for (unsigned int i = 0; i < rInterfaceShapeDerivatives.size(); ++i){
        EnrichedInterfaceShapeDerivativesNeg[i] = prod(enr_neg_interp, rInterfaceShapeDerivatives[i]);
    }

    double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double element_volume = positive_volume + negative_volume;

    const auto& r_geom = this->GetGeometry();
    const double h_elem = rData.ElementSize;

    double cut_area = 0.0;
    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){
        cut_area += rInterfaceWeights[gp];
    }

    const double density = 1.0/(1.0/positive_density + 1.0/negative_density);
    const double viscosity = 1.0/(1.0/positive_viscosity + 1.0/negative_viscosity);

    // Stabilization parameters
    const double cut_stabilization_coefficient = 1.0;
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    const double dyn_tau = rData.DynamicTau;

    const double dt = rData.DeltaTime;

    const auto v_convection = rData.Velocity - rData.MeshVelocity;

    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){

        Vector vconv = ZeroVector(Dim);
        double positive_weight = 0.0;
        double negative_weight = 0.0;

        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < Dim; ++dim){
                vconv[dim] += (rEnrInterfaceShapeFunctionNeg(gp, j) + rEnrInterfaceShapeFunctionPos(gp, j))
                    *v_convection(j,dim);
            }
            positive_weight += rEnrInterfaceShapeFunctionNeg(gp, j);
            negative_weight += rEnrInterfaceShapeFunctionPos(gp, j);
        }

        const double v_conv_norm = norm_2(vconv);

        const double penalty_coefficient = cut_stabilization_coefficient *
            density * 1.0 / (dyn_tau * density / dt + stab_c1 * viscosity / h_elem / h_elem +
                                stab_c2 * density * v_conv_norm / h_elem) * element_volume / cut_area;

        const auto& r_gp_enriched_interface_shape_derivatives_pos = EnrichedInterfaceShapeDerivativesPos[gp];
        const auto& r_gp_enriched_interface_shape_derivatives_neg = EnrichedInterfaceShapeDerivativesNeg[gp];

        for (unsigned int i = 0; i < NumNodes; ++i){

            for (unsigned int j = 0; j < NumNodes; ++j){

                const auto& r_pressure_gradient_j = r_geom[j].GetValue(PRESSURE_GRADIENT);

                for (unsigned int dim = 0; dim < Dim; ++dim){
                    kee(i, j) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        ( r_gp_enriched_interface_shape_derivatives_pos(j,dim) - r_gp_enriched_interface_shape_derivatives_neg(j,dim) );

                    rhs_enr(i) += penalty_coefficient * rInterfaceWeights[gp] *
                        ( r_gp_enriched_interface_shape_derivatives_pos(i,dim) - r_gp_enriched_interface_shape_derivatives_neg(i,dim) )*
                        (rEnrInterfaceShapeFunctionNeg(gp, j)/positive_weight - rEnrInterfaceShapeFunctionPos(gp, j)/negative_weight)*
                        r_pressure_gradient_j(dim);
                }
            }
        }
    }

    noalias(rKeeTot) += kee;
    noalias(rRHSeeTot) += rhs_enr;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CondenseEnrichmentWithContinuity(
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
            for (unsigned int j = i + 1; j < NumNodes; ++j){
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
void TwoFluidNavierStokes<TElementData>::CondenseEnrichment(
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
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
void TwoFluidNavierStokes<TElementData>::AddSurfaceTensionContribution(
    const TElementData& rData,
    ModifiedShapeFunctions::Pointer pModifiedShapeFunctions,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    VectorType &rRHSeeTot)
{
    // Surface tension coefficient is set in material properties
    const double surface_tension_coefficient = this->GetProperties().GetValue(SURFACE_TENSION_COEFFICIENT);
    Matrix int_shape_function, int_shape_function_enr_neg, int_shape_function_enr_pos;
    GeometryType::ShapeFunctionsGradientsType int_shape_derivatives;
    Vector int_gauss_pts_weights;
    std::vector<Vector> int_normals_neg;
    Vector gauss_pts_curvature;

    ComputeSplitInterface(
        rData,
        int_shape_function,
        int_shape_function_enr_pos,
        int_shape_function_enr_neg,
        int_shape_derivatives,
        int_gauss_pts_weights,
        int_normals_neg,
        pModifiedShapeFunctions);

    CalculateCurvatureOnInterfaceGaussPoints(
        int_shape_function,
        gauss_pts_curvature);

    SurfaceTension(
        surface_tension_coefficient,
        gauss_pts_curvature,
        int_gauss_pts_weights,
        int_shape_function,
        int_normals_neg,
        rRightHandSideVector);

    PressureGradientStabilization(
        rData,
        int_gauss_pts_weights,
        int_shape_function_enr_pos,
        int_shape_function_enr_neg,
        int_shape_derivatives,
        rKeeTot,
        rRHSeeTot);

    CondenseEnrichment(rLeftHandSideMatrix, rRightHandSideVector, rHtot, rVtot, rKeeTot, rRHSeeTot);
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
void TwoFluidNavierStokes<TElementData>::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable,
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

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; ++i_gauss){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; ++nnode){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ++ndim){
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