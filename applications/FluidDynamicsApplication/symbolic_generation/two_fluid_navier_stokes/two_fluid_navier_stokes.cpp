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
const double clhs9 =             rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs10 =             bdf0*rho;
const double clhs11 =             K_darcy*N[0];
const double clhs12 =             N[0]*clhs10;
const double clhs13 =             clhs11 + clhs12 + clhs9;
const double clhs14 =             1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 =             1.0*clhs9;
const double clhs16 =             clhs14*clhs15;
const double clhs17 =             1.0*clhs11;
const double clhs18 =             clhs14*clhs17;
const double clhs19 =             K_darcy*clhs8 + N[0]*clhs9 + clhs10*clhs8 + clhs13*clhs16 - clhs13*clhs18;
const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
const double clhs21 =             C(1,2)*DN(0,1);
const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,0)*clhs7;
const double clhs24 =             DN(0,1)*clhs23;
const double clhs25 =             -N[0] + clhs16 - clhs18;
const double clhs26 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs27 =             C(0,2)*DN(1,0);
const double clhs28 =             C(2,2)*DN(1,1) + clhs27;
const double clhs29 =             DN(0,0)*DN(1,0);
const double clhs30 =             N[1]*clhs11 + N[1]*clhs12;
const double clhs31 =             clhs29*clhs7 + clhs30;
const double clhs32 =             rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs33 =             K_darcy*N[1];
const double clhs34 =             N[1]*clhs10;
const double clhs35 =             clhs32 + clhs33 + clhs34;
const double clhs36 =             N[0]*clhs32 + clhs16*clhs35 - clhs18*clhs35;
const double clhs37 =             C(0,1)*DN(1,1) + clhs27;
const double clhs38 =             C(1,2)*DN(1,1);
const double clhs39 =             C(2,2)*DN(1,0) + clhs38;
const double clhs40 =             DN(1,1)*clhs23;
const double clhs41 =             DN(0,0)*N[1];
const double clhs42 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs43 =             C(0,2)*DN(2,0);
const double clhs44 =             C(2,2)*DN(2,1) + clhs43;
const double clhs45 =             DN(0,0)*DN(2,0);
const double clhs46 =             N[2]*clhs11 + N[2]*clhs12;
const double clhs47 =             clhs45*clhs7 + clhs46;
const double clhs48 =             rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs49 =             K_darcy*N[2];
const double clhs50 =             N[2]*clhs10;
const double clhs51 =             clhs48 + clhs49 + clhs50;
const double clhs52 =             N[0]*clhs48 + clhs16*clhs51 - clhs18*clhs51;
const double clhs53 =             C(0,1)*DN(2,1) + clhs43;
const double clhs54 =             C(1,2)*DN(2,1);
const double clhs55 =             C(2,2)*DN(2,0) + clhs54;
const double clhs56 =             DN(2,1)*clhs23;
const double clhs57 =             DN(0,0)*N[2];
const double clhs58 =             C(0,1)*DN(0,0) + clhs21;
const double clhs59 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs60 =             pow(DN(0,1), 2);
const double clhs61 =             C(0,1)*DN(1,0) + clhs38;
const double clhs62 =             DN(0,1)*clhs7;
const double clhs63 =             DN(1,0)*clhs62;
const double clhs64 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs65 =             DN(0,1)*DN(1,1);
const double clhs66 =             clhs30 + clhs65*clhs7;
const double clhs67 =             DN(0,1)*N[1];
const double clhs68 =             C(0,1)*DN(2,0) + clhs54;
const double clhs69 =             DN(2,0)*clhs62;
const double clhs70 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs71 =             DN(0,1)*DN(2,1);
const double clhs72 =             clhs46 + clhs7*clhs71;
const double clhs73 =             DN(0,1)*N[2];
const double clhs74 =             N[0] + clhs14*(1.0*clhs12 + clhs15 + clhs17);
const double clhs75 =             1.0*clhs14;
const double clhs76 =             DN(1,0)*N[0];
const double clhs77 =             clhs35*clhs75;
const double clhs78 =             DN(1,1)*N[0];
const double clhs79 =             clhs75*(clhs29 + clhs65);
const double clhs80 =             DN(2,0)*N[0];
const double clhs81 =             clhs51*clhs75;
const double clhs82 =             DN(2,1)*N[0];
const double clhs83 =             clhs75*(clhs45 + clhs71);
const double clhs84 =             clhs32*clhs75;
const double clhs85 =             clhs33*clhs75;
const double clhs86 =             N[1]*clhs9 + clhs13*clhs84 - clhs13*clhs85;
const double clhs87 =             pow(DN(1,0), 2);
const double clhs88 =             pow(N[1], 2);
const double clhs89 =             K_darcy*clhs88 + N[1]*clhs32 + clhs10*clhs88 + clhs32*clhs77 - clhs33*clhs77;
const double clhs90 =             DN(1,0)*clhs7;
const double clhs91 =             DN(1,1)*clhs90;
const double clhs92 =             -N[1] + clhs84 - clhs85;
const double clhs93 =             DN(1,0)*DN(2,0);
const double clhs94 =             N[2]*clhs33 + N[2]*clhs34;
const double clhs95 =             clhs7*clhs93 + clhs94;
const double clhs96 =             N[1]*clhs48 + clhs32*clhs81 - clhs33*clhs81;
const double clhs97 =             DN(2,1)*clhs90;
const double clhs98 =             DN(1,0)*N[2];
const double clhs99 =             pow(DN(1,1), 2);
const double clhs100 =             DN(2,0)*clhs7;
const double clhs101 =             DN(1,1)*clhs100;
const double clhs102 =             DN(1,1)*DN(2,1);
const double clhs103 =             clhs102*clhs7 + clhs94;
const double clhs104 =             DN(1,1)*N[2];
const double clhs105 =             clhs13*clhs75;
const double clhs106 =             N[1] + clhs14*(1.0*clhs32 + 1.0*clhs33 + 1.0*clhs34);
const double clhs107 =             DN(2,0)*N[1];
const double clhs108 =             DN(2,1)*N[1];
const double clhs109 =             clhs75*(clhs102 + clhs93);
const double clhs110 =             N[2]*clhs9 + clhs105*clhs48 - clhs105*clhs49;
const double clhs111 =             clhs49*clhs75;
const double clhs112 =             clhs48*clhs75;
const double clhs113 =             N[2]*clhs32 + clhs48*clhs77 - clhs49*clhs77;
const double clhs114 =             pow(DN(2,0), 2);
const double clhs115 =             pow(N[2], 2);
const double clhs116 =             K_darcy*clhs115 + N[2]*clhs48 + clhs10*clhs115 + clhs48*clhs81 - clhs49*clhs81;
const double clhs117 =             DN(2,1)*clhs100;
const double clhs118 =             -N[2] - clhs111 + clhs112;
const double clhs119 =             pow(DN(2,1), 2);
const double clhs120 =             N[2] + clhs14*(1.0*clhs48 + 1.0*clhs49 + 1.0*clhs50);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
            lhs(0,2)=DN(0,0)*clhs25;
            lhs(0,3)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs31 + clhs36;
            lhs(0,4)=DN(0,0)*clhs37 + DN(0,1)*clhs39 + clhs40;
            lhs(0,5)=DN(1,0)*clhs16 - DN(1,0)*clhs18 - clhs41;
            lhs(0,6)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + clhs47 + clhs52;
            lhs(0,7)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + clhs56;
            lhs(0,8)=DN(2,0)*clhs16 - DN(2,0)*clhs18 - clhs57;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs58 + clhs24;
            lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs59 + clhs19 + clhs60*clhs7;
            lhs(1,2)=DN(0,1)*clhs25;
            lhs(1,3)=DN(0,0)*clhs28 + DN(0,1)*clhs61 + clhs63;
            lhs(1,4)=DN(0,0)*clhs39 + DN(0,1)*clhs64 + clhs36 + clhs66;
            lhs(1,5)=DN(1,1)*clhs16 - DN(1,1)*clhs18 - clhs67;
            lhs(1,6)=DN(0,0)*clhs44 + DN(0,1)*clhs68 + clhs69;
            lhs(1,7)=DN(0,0)*clhs55 + DN(0,1)*clhs70 + clhs52 + clhs72;
            lhs(1,8)=DN(2,1)*clhs16 - DN(2,1)*clhs18 - clhs73;
            lhs(2,0)=DN(0,0)*clhs74;
            lhs(2,1)=DN(0,1)*clhs74;
            lhs(2,2)=clhs75*(clhs3 + clhs60);
            lhs(2,3)=DN(0,0)*clhs77 + clhs76;
            lhs(2,4)=DN(0,1)*clhs77 + clhs78;
            lhs(2,5)=clhs79;
            lhs(2,6)=DN(0,0)*clhs81 + clhs80;
            lhs(2,7)=DN(0,1)*clhs81 + clhs82;
            lhs(2,8)=clhs83;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs31 + clhs86;
            lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs63;
            lhs(3,2)=DN(0,0)*clhs84 - DN(0,0)*clhs85 - clhs76;
            lhs(3,3)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs7*clhs87 + clhs89;
            lhs(3,4)=DN(1,0)*clhs37 + DN(1,1)*clhs39 + clhs91;
            lhs(3,5)=DN(1,0)*clhs92;
            lhs(3,6)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + clhs95 + clhs96;
            lhs(3,7)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + clhs97;
            lhs(3,8)=DN(2,0)*clhs84 - DN(2,0)*clhs85 - clhs98;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs58 + clhs40;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs59 + clhs66 + clhs86;
            lhs(4,2)=DN(0,1)*clhs84 - DN(0,1)*clhs85 - clhs78;
            lhs(4,3)=DN(1,0)*clhs28 + DN(1,1)*clhs61 + clhs91;
            lhs(4,4)=DN(1,0)*clhs39 + DN(1,1)*clhs64 + clhs7*clhs99 + clhs89;
            lhs(4,5)=DN(1,1)*clhs92;
            lhs(4,6)=DN(1,0)*clhs44 + DN(1,1)*clhs68 + clhs101;
            lhs(4,7)=DN(1,0)*clhs55 + DN(1,1)*clhs70 + clhs103 + clhs96;
            lhs(4,8)=DN(2,1)*clhs84 - DN(2,1)*clhs85 - clhs104;
            lhs(5,0)=DN(1,0)*clhs105 + clhs41;
            lhs(5,1)=DN(1,1)*clhs105 + clhs67;
            lhs(5,2)=clhs79;
            lhs(5,3)=DN(1,0)*clhs106;
            lhs(5,4)=DN(1,1)*clhs106;
            lhs(5,5)=clhs75*(clhs87 + clhs99);
            lhs(5,6)=DN(1,0)*clhs81 + clhs107;
            lhs(5,7)=DN(1,1)*clhs81 + clhs108;
            lhs(5,8)=clhs109;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs110 + clhs47;
            lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs69;
            lhs(6,2)=-DN(0,0)*clhs111 + DN(0,0)*clhs112 - clhs80;
            lhs(6,3)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs113 + clhs95;
            lhs(6,4)=DN(2,0)*clhs37 + DN(2,1)*clhs39 + clhs101;
            lhs(6,5)=-DN(1,0)*clhs111 + DN(1,0)*clhs112 - clhs107;
            lhs(6,6)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + clhs114*clhs7 + clhs116;
            lhs(6,7)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + clhs117;
            lhs(6,8)=DN(2,0)*clhs118;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs58 + clhs56;
            lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs59 + clhs110 + clhs72;
            lhs(7,2)=-DN(0,1)*clhs111 + DN(0,1)*clhs112 - clhs82;
            lhs(7,3)=DN(2,0)*clhs28 + DN(2,1)*clhs61 + clhs97;
            lhs(7,4)=DN(2,0)*clhs39 + DN(2,1)*clhs64 + clhs103 + clhs113;
            lhs(7,5)=-DN(1,1)*clhs111 + DN(1,1)*clhs112 - clhs108;
            lhs(7,6)=DN(2,0)*clhs44 + DN(2,1)*clhs68 + clhs117;
            lhs(7,7)=DN(2,0)*clhs55 + DN(2,1)*clhs70 + clhs116 + clhs119*clhs7;
            lhs(7,8)=DN(2,1)*clhs118;
            lhs(8,0)=DN(2,0)*clhs105 + clhs57;
            lhs(8,1)=DN(2,1)*clhs105 + clhs73;
            lhs(8,2)=clhs83;
            lhs(8,3)=DN(2,0)*clhs77 + clhs98;
            lhs(8,4)=DN(2,1)*clhs77 + clhs104;
            lhs(8,5)=clhs109;
            lhs(8,6)=DN(2,0)*clhs120;
            lhs(8,7)=DN(2,1)*clhs120;
            lhs(8,8)=clhs75*(clhs114 + clhs119);


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
const double clhs12 =             rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs13 =             bdf0*rho;
const double clhs14 =             K_darcy*N[0];
const double clhs15 =             N[0]*clhs13;
const double clhs16 =             clhs12 + clhs14 + clhs15;
const double clhs17 =             1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 =             1.0*clhs12;
const double clhs19 =             clhs17*clhs18;
const double clhs20 =             1.0*clhs14;
const double clhs21 =             clhs17*clhs20;
const double clhs22 =             K_darcy*clhs11 + N[0]*clhs12 + clhs11*clhs13 + clhs16*clhs19 - clhs16*clhs21;
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
const double clhs37 =             -N[0] + clhs19 - clhs21;
const double clhs38 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs39 =             C(0,3)*DN(1,0);
const double clhs40 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs39;
const double clhs41 =             C(0,5)*DN(1,0);
const double clhs42 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs41;
const double clhs43 =             DN(0,0)*DN(1,0);
const double clhs44 =             N[1]*clhs14 + N[1]*clhs15;
const double clhs45 =             clhs10*clhs43 + clhs44;
const double clhs46 =             rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs47 =             K_darcy*N[1];
const double clhs48 =             N[1]*clhs13;
const double clhs49 =             clhs46 + clhs47 + clhs48;
const double clhs50 =             N[0]*clhs46 + clhs19*clhs49 - clhs21*clhs49;
const double clhs51 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs39;
const double clhs52 =             C(1,3)*DN(1,1);
const double clhs53 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs52;
const double clhs54 =             C(3,5)*DN(1,0);
const double clhs55 =             C(4,5)*DN(1,2);
const double clhs56 =             C(1,5)*DN(1,1) + clhs54 + clhs55;
const double clhs57 =             DN(1,1)*clhs29;
const double clhs58 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs41;
const double clhs59 =             C(3,4)*DN(1,1);
const double clhs60 =             C(2,3)*DN(1,2) + clhs54 + clhs59;
const double clhs61 =             C(2,5)*DN(1,2);
const double clhs62 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs61;
const double clhs63 =             DN(1,2)*clhs29;
const double clhs64 =             DN(0,0)*N[1];
const double clhs65 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs66 =             C(0,3)*DN(2,0);
const double clhs67 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs66;
const double clhs68 =             C(0,5)*DN(2,0);
const double clhs69 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs68;
const double clhs70 =             DN(0,0)*DN(2,0);
const double clhs71 =             N[2]*clhs14 + N[2]*clhs15;
const double clhs72 =             clhs10*clhs70 + clhs71;
const double clhs73 =             rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs74 =             K_darcy*N[2];
const double clhs75 =             N[2]*clhs13;
const double clhs76 =             clhs73 + clhs74 + clhs75;
const double clhs77 =             N[0]*clhs73 + clhs19*clhs76 - clhs21*clhs76;
const double clhs78 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs66;
const double clhs79 =             C(1,3)*DN(2,1);
const double clhs80 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs79;
const double clhs81 =             C(3,5)*DN(2,0);
const double clhs82 =             C(4,5)*DN(2,2);
const double clhs83 =             C(1,5)*DN(2,1) + clhs81 + clhs82;
const double clhs84 =             DN(2,1)*clhs29;
const double clhs85 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs68;
const double clhs86 =             C(3,4)*DN(2,1);
const double clhs87 =             C(2,3)*DN(2,2) + clhs81 + clhs86;
const double clhs88 =             C(2,5)*DN(2,2);
const double clhs89 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs88;
const double clhs90 =             DN(2,2)*clhs29;
const double clhs91 =             DN(0,0)*N[2];
const double clhs92 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs93 =             C(0,3)*DN(3,0);
const double clhs94 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs93;
const double clhs95 =             C(0,5)*DN(3,0);
const double clhs96 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs95;
const double clhs97 =             DN(0,0)*DN(3,0);
const double clhs98 =             N[3]*clhs14 + N[3]*clhs15;
const double clhs99 =             clhs10*clhs97 + clhs98;
const double clhs100 =             rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs101 =             K_darcy*N[3];
const double clhs102 =             N[3]*clhs13;
const double clhs103 =             clhs100 + clhs101 + clhs102;
const double clhs104 =             N[0]*clhs100 + clhs103*clhs19 - clhs103*clhs21;
const double clhs105 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs93;
const double clhs106 =             C(1,3)*DN(3,1);
const double clhs107 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs106;
const double clhs108 =             C(3,5)*DN(3,0);
const double clhs109 =             C(4,5)*DN(3,2);
const double clhs110 =             C(1,5)*DN(3,1) + clhs108 + clhs109;
const double clhs111 =             DN(3,1)*clhs29;
const double clhs112 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs95;
const double clhs113 =             C(3,4)*DN(3,1);
const double clhs114 =             C(2,3)*DN(3,2) + clhs108 + clhs113;
const double clhs115 =             C(2,5)*DN(3,2);
const double clhs116 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs115;
const double clhs117 =             DN(3,2)*clhs29;
const double clhs118 =             DN(0,0)*N[3];
const double clhs119 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs120 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs121 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs122 =             C(1,4)*DN(0,1);
const double clhs123 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs122;
const double clhs124 =             pow(DN(0,1), 2);
const double clhs125 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs122;
const double clhs126 =             C(2,4)*DN(0,2);
const double clhs127 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs126;
const double clhs128 =             DN(0,1)*clhs10;
const double clhs129 =             DN(0,2)*clhs128;
const double clhs130 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs52;
const double clhs131 =             C(0,4)*DN(1,0) + clhs55 + clhs59;
const double clhs132 =             DN(1,0)*clhs128;
const double clhs133 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs134 =             C(1,4)*DN(1,1);
const double clhs135 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs134;
const double clhs136 =             DN(0,1)*DN(1,1);
const double clhs137 =             clhs10*clhs136;
const double clhs138 =             clhs44 + clhs50;
const double clhs139 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs134;
const double clhs140 =             C(2,4)*DN(1,2);
const double clhs141 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs140;
const double clhs142 =             DN(1,2)*clhs128;
const double clhs143 =             DN(0,1)*N[1];
const double clhs144 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs79;
const double clhs145 =             C(0,4)*DN(2,0) + clhs82 + clhs86;
const double clhs146 =             DN(2,0)*clhs128;
const double clhs147 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs148 =             C(1,4)*DN(2,1);
const double clhs149 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs148;
const double clhs150 =             DN(0,1)*DN(2,1);
const double clhs151 =             clhs10*clhs150;
const double clhs152 =             clhs71 + clhs77;
const double clhs153 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs148;
const double clhs154 =             C(2,4)*DN(2,2);
const double clhs155 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs154;
const double clhs156 =             DN(2,2)*clhs128;
const double clhs157 =             DN(0,1)*N[2];
const double clhs158 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs106;
const double clhs159 =             C(0,4)*DN(3,0) + clhs109 + clhs113;
const double clhs160 =             DN(3,0)*clhs128;
const double clhs161 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs162 =             C(1,4)*DN(3,1);
const double clhs163 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs162;
const double clhs164 =             DN(0,1)*DN(3,1);
const double clhs165 =             clhs10*clhs164;
const double clhs166 =             clhs104 + clhs98;
const double clhs167 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs162;
const double clhs168 =             C(2,4)*DN(3,2);
const double clhs169 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs168;
const double clhs170 =             DN(3,2)*clhs128;
const double clhs171 =             DN(0,1)*N[3];
const double clhs172 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs173 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs126;
const double clhs174 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs175 =             pow(DN(0,2), 2);
const double clhs176 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs61;
const double clhs177 =             DN(0,2)*clhs10;
const double clhs178 =             DN(1,0)*clhs177;
const double clhs179 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs140;
const double clhs180 =             DN(1,1)*clhs177;
const double clhs181 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs182 =             DN(0,2)*DN(1,2);
const double clhs183 =             clhs10*clhs182;
const double clhs184 =             DN(0,2)*N[1];
const double clhs185 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs88;
const double clhs186 =             DN(2,0)*clhs177;
const double clhs187 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs154;
const double clhs188 =             DN(2,1)*clhs177;
const double clhs189 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs190 =             DN(0,2)*DN(2,2);
const double clhs191 =             clhs10*clhs190;
const double clhs192 =             DN(0,2)*N[2];
const double clhs193 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs115;
const double clhs194 =             DN(3,0)*clhs177;
const double clhs195 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs168;
const double clhs196 =             DN(3,1)*clhs177;
const double clhs197 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs198 =             DN(0,2)*DN(3,2);
const double clhs199 =             clhs10*clhs198;
const double clhs200 =             DN(0,2)*N[3];
const double clhs201 =             N[0] + clhs17*(1.0*clhs15 + clhs18 + clhs20);
const double clhs202 =             1.0*clhs17;
const double clhs203 =             DN(1,0)*N[0];
const double clhs204 =             clhs202*clhs49;
const double clhs205 =             DN(1,1)*N[0];
const double clhs206 =             DN(1,2)*N[0];
const double clhs207 =             clhs202*(clhs136 + clhs182 + clhs43);
const double clhs208 =             DN(2,0)*N[0];
const double clhs209 =             clhs202*clhs76;
const double clhs210 =             DN(2,1)*N[0];
const double clhs211 =             DN(2,2)*N[0];
const double clhs212 =             clhs202*(clhs150 + clhs190 + clhs70);
const double clhs213 =             DN(3,0)*N[0];
const double clhs214 =             clhs103*clhs202;
const double clhs215 =             DN(3,1)*N[0];
const double clhs216 =             DN(3,2)*N[0];
const double clhs217 =             clhs202*(clhs164 + clhs198 + clhs97);
const double clhs218 =             clhs202*clhs46;
const double clhs219 =             clhs202*clhs47;
const double clhs220 =             N[1]*clhs12 + clhs16*clhs218 - clhs16*clhs219;
const double clhs221 =             pow(DN(1,0), 2);
const double clhs222 =             pow(N[1], 2);
const double clhs223 =             K_darcy*clhs222 + N[1]*clhs46 + clhs13*clhs222 + clhs204*clhs46 - clhs204*clhs47;
const double clhs224 =             DN(1,0)*clhs10;
const double clhs225 =             DN(1,1)*clhs224;
const double clhs226 =             DN(1,2)*clhs224;
const double clhs227 =             -N[1] + clhs218 - clhs219;
const double clhs228 =             DN(1,0)*DN(2,0);
const double clhs229 =             N[2]*clhs47 + N[2]*clhs48;
const double clhs230 =             clhs10*clhs228 + clhs229;
const double clhs231 =             N[1]*clhs73 + clhs209*clhs46 - clhs209*clhs47;
const double clhs232 =             DN(2,1)*clhs224;
const double clhs233 =             DN(2,2)*clhs224;
const double clhs234 =             DN(1,0)*N[2];
const double clhs235 =             DN(1,0)*DN(3,0);
const double clhs236 =             N[3]*clhs47 + N[3]*clhs48;
const double clhs237 =             clhs10*clhs235 + clhs236;
const double clhs238 =             N[1]*clhs100 + clhs214*clhs46 - clhs214*clhs47;
const double clhs239 =             DN(3,1)*clhs224;
const double clhs240 =             DN(3,2)*clhs224;
const double clhs241 =             DN(1,0)*N[3];
const double clhs242 =             clhs220 + clhs44;
const double clhs243 =             pow(DN(1,1), 2);
const double clhs244 =             DN(1,1)*clhs10;
const double clhs245 =             DN(1,2)*clhs244;
const double clhs246 =             DN(2,0)*clhs244;
const double clhs247 =             DN(1,1)*DN(2,1);
const double clhs248 =             clhs10*clhs247;
const double clhs249 =             clhs229 + clhs231;
const double clhs250 =             DN(2,2)*clhs244;
const double clhs251 =             DN(1,1)*N[2];
const double clhs252 =             DN(3,0)*clhs244;
const double clhs253 =             DN(1,1)*DN(3,1);
const double clhs254 =             clhs10*clhs253;
const double clhs255 =             clhs236 + clhs238;
const double clhs256 =             DN(3,2)*clhs244;
const double clhs257 =             DN(1,1)*N[3];
const double clhs258 =             pow(DN(1,2), 2);
const double clhs259 =             DN(1,2)*clhs10;
const double clhs260 =             DN(2,0)*clhs259;
const double clhs261 =             DN(2,1)*clhs259;
const double clhs262 =             DN(1,2)*DN(2,2);
const double clhs263 =             clhs10*clhs262;
const double clhs264 =             DN(1,2)*N[2];
const double clhs265 =             DN(3,0)*clhs259;
const double clhs266 =             DN(3,1)*clhs259;
const double clhs267 =             DN(1,2)*DN(3,2);
const double clhs268 =             clhs10*clhs267;
const double clhs269 =             DN(1,2)*N[3];
const double clhs270 =             clhs16*clhs202;
const double clhs271 =             N[1] + clhs17*(1.0*clhs46 + 1.0*clhs47 + 1.0*clhs48);
const double clhs272 =             DN(2,0)*N[1];
const double clhs273 =             DN(2,1)*N[1];
const double clhs274 =             DN(2,2)*N[1];
const double clhs275 =             clhs202*(clhs228 + clhs247 + clhs262);
const double clhs276 =             DN(3,0)*N[1];
const double clhs277 =             DN(3,1)*N[1];
const double clhs278 =             DN(3,2)*N[1];
const double clhs279 =             clhs202*(clhs235 + clhs253 + clhs267);
const double clhs280 =             N[2]*clhs12 + clhs270*clhs73 - clhs270*clhs74;
const double clhs281 =             clhs202*clhs74;
const double clhs282 =             clhs202*clhs73;
const double clhs283 =             N[2]*clhs46 + clhs204*clhs73 - clhs204*clhs74;
const double clhs284 =             pow(DN(2,0), 2);
const double clhs285 =             pow(N[2], 2);
const double clhs286 =             K_darcy*clhs285 + N[2]*clhs73 + clhs13*clhs285 + clhs209*clhs73 - clhs209*clhs74;
const double clhs287 =             DN(2,0)*clhs10;
const double clhs288 =             DN(2,1)*clhs287;
const double clhs289 =             DN(2,2)*clhs287;
const double clhs290 =             -N[2] - clhs281 + clhs282;
const double clhs291 =             DN(2,0)*DN(3,0);
const double clhs292 =             N[3]*clhs74 + N[3]*clhs75;
const double clhs293 =             clhs10*clhs291 + clhs292;
const double clhs294 =             N[2]*clhs100 + clhs214*clhs73 - clhs214*clhs74;
const double clhs295 =             DN(3,1)*clhs287;
const double clhs296 =             DN(3,2)*clhs287;
const double clhs297 =             DN(2,0)*N[3];
const double clhs298 =             clhs280 + clhs71;
const double clhs299 =             clhs229 + clhs283;
const double clhs300 =             pow(DN(2,1), 2);
const double clhs301 =             DN(2,1)*clhs10;
const double clhs302 =             DN(2,2)*clhs301;
const double clhs303 =             DN(3,0)*clhs301;
const double clhs304 =             DN(2,1)*DN(3,1);
const double clhs305 =             clhs10*clhs304;
const double clhs306 =             clhs292 + clhs294;
const double clhs307 =             DN(3,2)*clhs301;
const double clhs308 =             DN(2,1)*N[3];
const double clhs309 =             pow(DN(2,2), 2);
const double clhs310 =             DN(2,2)*clhs10;
const double clhs311 =             DN(3,0)*clhs310;
const double clhs312 =             DN(3,1)*clhs310;
const double clhs313 =             DN(2,2)*DN(3,2);
const double clhs314 =             clhs10*clhs313;
const double clhs315 =             DN(2,2)*N[3];
const double clhs316 =             N[2] + clhs17*(1.0*clhs73 + 1.0*clhs74 + 1.0*clhs75);
const double clhs317 =             DN(3,0)*N[2];
const double clhs318 =             DN(3,1)*N[2];
const double clhs319 =             DN(3,2)*N[2];
const double clhs320 =             clhs202*(clhs291 + clhs304 + clhs313);
const double clhs321 =             N[3]*clhs12 + clhs100*clhs270 - clhs101*clhs270;
const double clhs322 =             clhs101*clhs202;
const double clhs323 =             clhs100*clhs202;
const double clhs324 =             N[3]*clhs46 + clhs100*clhs204 - clhs101*clhs204;
const double clhs325 =             N[3]*clhs73 + clhs100*clhs209 - clhs101*clhs209;
const double clhs326 =             pow(DN(3,0), 2);
const double clhs327 =             pow(N[3], 2);
const double clhs328 =             K_darcy*clhs327 + N[3]*clhs100 + clhs100*clhs214 - clhs101*clhs214 + clhs13*clhs327;
const double clhs329 =             DN(3,0)*clhs10;
const double clhs330 =             DN(3,1)*clhs329;
const double clhs331 =             DN(3,2)*clhs329;
const double clhs332 =             -N[3] - clhs322 + clhs323;
const double clhs333 =             clhs321 + clhs98;
const double clhs334 =             clhs236 + clhs324;
const double clhs335 =             clhs292 + clhs325;
const double clhs336 =             pow(DN(3,1), 2);
const double clhs337 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs338 =             pow(DN(3,2), 2);
const double clhs339 =             N[3] + clhs17*(1.0*clhs100 + 1.0*clhs101 + 1.0*clhs102);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
            lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
            lhs(0,3)=DN(0,0)*clhs37;
            lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + DN(0,2)*clhs42 + clhs45 + clhs50;
            lhs(0,5)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + DN(0,2)*clhs56 + clhs57;
            lhs(0,6)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs62 + clhs63;
            lhs(0,7)=DN(1,0)*clhs19 - DN(1,0)*clhs21 - clhs64;
            lhs(0,8)=DN(0,0)*clhs65 + DN(0,1)*clhs67 + DN(0,2)*clhs69 + clhs72 + clhs77;
            lhs(0,9)=DN(0,0)*clhs78 + DN(0,1)*clhs80 + DN(0,2)*clhs83 + clhs84;
            lhs(0,10)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs89 + clhs90;
            lhs(0,11)=DN(2,0)*clhs19 - DN(2,0)*clhs21 - clhs91;
            lhs(0,12)=DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs96 + clhs104 + clhs99;
            lhs(0,13)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs110 + clhs111;
            lhs(0,14)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs117;
            lhs(0,15)=DN(3,0)*clhs19 - DN(3,0)*clhs21 - clhs118;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs119 + DN(0,2)*clhs120 + clhs30;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs10*clhs124 + clhs22;
            lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs125 + DN(0,2)*clhs127 + clhs129;
            lhs(1,3)=DN(0,1)*clhs37;
            lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs132;
            lhs(1,5)=DN(0,0)*clhs53 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs137 + clhs138;
            lhs(1,6)=DN(0,0)*clhs60 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142;
            lhs(1,7)=DN(1,1)*clhs19 - DN(1,1)*clhs21 - clhs143;
            lhs(1,8)=DN(0,0)*clhs67 + DN(0,1)*clhs144 + DN(0,2)*clhs145 + clhs146;
            lhs(1,9)=DN(0,0)*clhs80 + DN(0,1)*clhs147 + DN(0,2)*clhs149 + clhs151 + clhs152;
            lhs(1,10)=DN(0,0)*clhs87 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs156;
            lhs(1,11)=DN(2,1)*clhs19 - DN(2,1)*clhs21 - clhs157;
            lhs(1,12)=DN(0,0)*clhs94 + DN(0,1)*clhs158 + DN(0,2)*clhs159 + clhs160;
            lhs(1,13)=DN(0,0)*clhs107 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs165 + clhs166;
            lhs(1,14)=DN(0,0)*clhs114 + DN(0,1)*clhs167 + DN(0,2)*clhs169 + clhs170;
            lhs(1,15)=DN(3,1)*clhs19 - DN(3,1)*clhs21 - clhs171;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs120 + DN(0,2)*clhs172 + clhs36;
            lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs123 + DN(0,2)*clhs173 + clhs129;
            lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs127 + DN(0,2)*clhs174 + clhs10*clhs175 + clhs22;
            lhs(2,3)=DN(0,2)*clhs37;
            lhs(2,4)=DN(0,0)*clhs42 + DN(0,1)*clhs131 + DN(0,2)*clhs176 + clhs178;
            lhs(2,5)=DN(0,0)*clhs56 + DN(0,1)*clhs135 + DN(0,2)*clhs179 + clhs180;
            lhs(2,6)=DN(0,0)*clhs62 + DN(0,1)*clhs141 + DN(0,2)*clhs181 + clhs138 + clhs183;
            lhs(2,7)=DN(1,2)*clhs19 - DN(1,2)*clhs21 - clhs184;
            lhs(2,8)=DN(0,0)*clhs69 + DN(0,1)*clhs145 + DN(0,2)*clhs185 + clhs186;
            lhs(2,9)=DN(0,0)*clhs83 + DN(0,1)*clhs149 + DN(0,2)*clhs187 + clhs188;
            lhs(2,10)=DN(0,0)*clhs89 + DN(0,1)*clhs155 + DN(0,2)*clhs189 + clhs152 + clhs191;
            lhs(2,11)=DN(2,2)*clhs19 - DN(2,2)*clhs21 - clhs192;
            lhs(2,12)=DN(0,0)*clhs96 + DN(0,1)*clhs159 + DN(0,2)*clhs193 + clhs194;
            lhs(2,13)=DN(0,0)*clhs110 + DN(0,1)*clhs163 + DN(0,2)*clhs195 + clhs196;
            lhs(2,14)=DN(0,0)*clhs116 + DN(0,1)*clhs169 + DN(0,2)*clhs197 + clhs166 + clhs199;
            lhs(2,15)=DN(3,2)*clhs19 - DN(3,2)*clhs21 - clhs200;
            lhs(3,0)=DN(0,0)*clhs201;
            lhs(3,1)=DN(0,1)*clhs201;
            lhs(3,2)=DN(0,2)*clhs201;
            lhs(3,3)=clhs202*(clhs124 + clhs175 + clhs5);
            lhs(3,4)=DN(0,0)*clhs204 + clhs203;
            lhs(3,5)=DN(0,1)*clhs204 + clhs205;
            lhs(3,6)=DN(0,2)*clhs204 + clhs206;
            lhs(3,7)=clhs207;
            lhs(3,8)=DN(0,0)*clhs209 + clhs208;
            lhs(3,9)=DN(0,1)*clhs209 + clhs210;
            lhs(3,10)=DN(0,2)*clhs209 + clhs211;
            lhs(3,11)=clhs212;
            lhs(3,12)=DN(0,0)*clhs214 + clhs213;
            lhs(3,13)=DN(0,1)*clhs214 + clhs215;
            lhs(3,14)=DN(0,2)*clhs214 + clhs216;
            lhs(3,15)=clhs217;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs220 + clhs45;
            lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs132;
            lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs178;
            lhs(4,3)=DN(0,0)*clhs218 - DN(0,0)*clhs219 - clhs203;
            lhs(4,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + DN(1,2)*clhs42 + clhs10*clhs221 + clhs223;
            lhs(4,5)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + DN(1,2)*clhs56 + clhs225;
            lhs(4,6)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs62 + clhs226;
            lhs(4,7)=DN(1,0)*clhs227;
            lhs(4,8)=DN(1,0)*clhs65 + DN(1,1)*clhs67 + DN(1,2)*clhs69 + clhs230 + clhs231;
            lhs(4,9)=DN(1,0)*clhs78 + DN(1,1)*clhs80 + DN(1,2)*clhs83 + clhs232;
            lhs(4,10)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs89 + clhs233;
            lhs(4,11)=DN(2,0)*clhs218 - DN(2,0)*clhs219 - clhs234;
            lhs(4,12)=DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs96 + clhs237 + clhs238;
            lhs(4,13)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs110 + clhs239;
            lhs(4,14)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs240;
            lhs(4,15)=DN(3,0)*clhs218 - DN(3,0)*clhs219 - clhs241;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs119 + DN(1,2)*clhs120 + clhs57;
            lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs137 + clhs242;
            lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs125 + DN(1,2)*clhs127 + clhs180;
            lhs(5,3)=DN(0,1)*clhs218 - DN(0,1)*clhs219 - clhs205;
            lhs(5,4)=DN(1,0)*clhs40 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs225;
            lhs(5,5)=DN(1,0)*clhs53 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs10*clhs243 + clhs223;
            lhs(5,6)=DN(1,0)*clhs60 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs245;
            lhs(5,7)=DN(1,1)*clhs227;
            lhs(5,8)=DN(1,0)*clhs67 + DN(1,1)*clhs144 + DN(1,2)*clhs145 + clhs246;
            lhs(5,9)=DN(1,0)*clhs80 + DN(1,1)*clhs147 + DN(1,2)*clhs149 + clhs248 + clhs249;
            lhs(5,10)=DN(1,0)*clhs87 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs250;
            lhs(5,11)=DN(2,1)*clhs218 - DN(2,1)*clhs219 - clhs251;
            lhs(5,12)=DN(1,0)*clhs94 + DN(1,1)*clhs158 + DN(1,2)*clhs159 + clhs252;
            lhs(5,13)=DN(1,0)*clhs107 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs254 + clhs255;
            lhs(5,14)=DN(1,0)*clhs114 + DN(1,1)*clhs167 + DN(1,2)*clhs169 + clhs256;
            lhs(5,15)=DN(3,1)*clhs218 - DN(3,1)*clhs219 - clhs257;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs120 + DN(1,2)*clhs172 + clhs63;
            lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs123 + DN(1,2)*clhs173 + clhs142;
            lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs127 + DN(1,2)*clhs174 + clhs183 + clhs242;
            lhs(6,3)=DN(0,2)*clhs218 - DN(0,2)*clhs219 - clhs206;
            lhs(6,4)=DN(1,0)*clhs42 + DN(1,1)*clhs131 + DN(1,2)*clhs176 + clhs226;
            lhs(6,5)=DN(1,0)*clhs56 + DN(1,1)*clhs135 + DN(1,2)*clhs179 + clhs245;
            lhs(6,6)=DN(1,0)*clhs62 + DN(1,1)*clhs141 + DN(1,2)*clhs181 + clhs10*clhs258 + clhs223;
            lhs(6,7)=DN(1,2)*clhs227;
            lhs(6,8)=DN(1,0)*clhs69 + DN(1,1)*clhs145 + DN(1,2)*clhs185 + clhs260;
            lhs(6,9)=DN(1,0)*clhs83 + DN(1,1)*clhs149 + DN(1,2)*clhs187 + clhs261;
            lhs(6,10)=DN(1,0)*clhs89 + DN(1,1)*clhs155 + DN(1,2)*clhs189 + clhs249 + clhs263;
            lhs(6,11)=DN(2,2)*clhs218 - DN(2,2)*clhs219 - clhs264;
            lhs(6,12)=DN(1,0)*clhs96 + DN(1,1)*clhs159 + DN(1,2)*clhs193 + clhs265;
            lhs(6,13)=DN(1,0)*clhs110 + DN(1,1)*clhs163 + DN(1,2)*clhs195 + clhs266;
            lhs(6,14)=DN(1,0)*clhs116 + DN(1,1)*clhs169 + DN(1,2)*clhs197 + clhs255 + clhs268;
            lhs(6,15)=DN(3,2)*clhs218 - DN(3,2)*clhs219 - clhs269;
            lhs(7,0)=DN(1,0)*clhs270 + clhs64;
            lhs(7,1)=DN(1,1)*clhs270 + clhs143;
            lhs(7,2)=DN(1,2)*clhs270 + clhs184;
            lhs(7,3)=clhs207;
            lhs(7,4)=DN(1,0)*clhs271;
            lhs(7,5)=DN(1,1)*clhs271;
            lhs(7,6)=DN(1,2)*clhs271;
            lhs(7,7)=clhs202*(clhs221 + clhs243 + clhs258);
            lhs(7,8)=DN(1,0)*clhs209 + clhs272;
            lhs(7,9)=DN(1,1)*clhs209 + clhs273;
            lhs(7,10)=DN(1,2)*clhs209 + clhs274;
            lhs(7,11)=clhs275;
            lhs(7,12)=DN(1,0)*clhs214 + clhs276;
            lhs(7,13)=DN(1,1)*clhs214 + clhs277;
            lhs(7,14)=DN(1,2)*clhs214 + clhs278;
            lhs(7,15)=clhs279;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs280 + clhs72;
            lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs146;
            lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs186;
            lhs(8,3)=-DN(0,0)*clhs281 + DN(0,0)*clhs282 - clhs208;
            lhs(8,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + DN(2,2)*clhs42 + clhs230 + clhs283;
            lhs(8,5)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + DN(2,2)*clhs56 + clhs246;
            lhs(8,6)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs62 + clhs260;
            lhs(8,7)=-DN(1,0)*clhs281 + DN(1,0)*clhs282 - clhs272;
            lhs(8,8)=DN(2,0)*clhs65 + DN(2,1)*clhs67 + DN(2,2)*clhs69 + clhs10*clhs284 + clhs286;
            lhs(8,9)=DN(2,0)*clhs78 + DN(2,1)*clhs80 + DN(2,2)*clhs83 + clhs288;
            lhs(8,10)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs89 + clhs289;
            lhs(8,11)=DN(2,0)*clhs290;
            lhs(8,12)=DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs96 + clhs293 + clhs294;
            lhs(8,13)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs110 + clhs295;
            lhs(8,14)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs296;
            lhs(8,15)=-DN(3,0)*clhs281 + DN(3,0)*clhs282 - clhs297;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs119 + DN(2,2)*clhs120 + clhs84;
            lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs151 + clhs298;
            lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs125 + DN(2,2)*clhs127 + clhs188;
            lhs(9,3)=-DN(0,1)*clhs281 + DN(0,1)*clhs282 - clhs210;
            lhs(9,4)=DN(2,0)*clhs40 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs232;
            lhs(9,5)=DN(2,0)*clhs53 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs248 + clhs299;
            lhs(9,6)=DN(2,0)*clhs60 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs261;
            lhs(9,7)=-DN(1,1)*clhs281 + DN(1,1)*clhs282 - clhs273;
            lhs(9,8)=DN(2,0)*clhs67 + DN(2,1)*clhs144 + DN(2,2)*clhs145 + clhs288;
            lhs(9,9)=DN(2,0)*clhs80 + DN(2,1)*clhs147 + DN(2,2)*clhs149 + clhs10*clhs300 + clhs286;
            lhs(9,10)=DN(2,0)*clhs87 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs302;
            lhs(9,11)=DN(2,1)*clhs290;
            lhs(9,12)=DN(2,0)*clhs94 + DN(2,1)*clhs158 + DN(2,2)*clhs159 + clhs303;
            lhs(9,13)=DN(2,0)*clhs107 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs305 + clhs306;
            lhs(9,14)=DN(2,0)*clhs114 + DN(2,1)*clhs167 + DN(2,2)*clhs169 + clhs307;
            lhs(9,15)=-DN(3,1)*clhs281 + DN(3,1)*clhs282 - clhs308;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs120 + DN(2,2)*clhs172 + clhs90;
            lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs123 + DN(2,2)*clhs173 + clhs156;
            lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs127 + DN(2,2)*clhs174 + clhs191 + clhs298;
            lhs(10,3)=-DN(0,2)*clhs281 + DN(0,2)*clhs282 - clhs211;
            lhs(10,4)=DN(2,0)*clhs42 + DN(2,1)*clhs131 + DN(2,2)*clhs176 + clhs233;
            lhs(10,5)=DN(2,0)*clhs56 + DN(2,1)*clhs135 + DN(2,2)*clhs179 + clhs250;
            lhs(10,6)=DN(2,0)*clhs62 + DN(2,1)*clhs141 + DN(2,2)*clhs181 + clhs263 + clhs299;
            lhs(10,7)=-DN(1,2)*clhs281 + DN(1,2)*clhs282 - clhs274;
            lhs(10,8)=DN(2,0)*clhs69 + DN(2,1)*clhs145 + DN(2,2)*clhs185 + clhs289;
            lhs(10,9)=DN(2,0)*clhs83 + DN(2,1)*clhs149 + DN(2,2)*clhs187 + clhs302;
            lhs(10,10)=DN(2,0)*clhs89 + DN(2,1)*clhs155 + DN(2,2)*clhs189 + clhs10*clhs309 + clhs286;
            lhs(10,11)=DN(2,2)*clhs290;
            lhs(10,12)=DN(2,0)*clhs96 + DN(2,1)*clhs159 + DN(2,2)*clhs193 + clhs311;
            lhs(10,13)=DN(2,0)*clhs110 + DN(2,1)*clhs163 + DN(2,2)*clhs195 + clhs312;
            lhs(10,14)=DN(2,0)*clhs116 + DN(2,1)*clhs169 + DN(2,2)*clhs197 + clhs306 + clhs314;
            lhs(10,15)=-DN(3,2)*clhs281 + DN(3,2)*clhs282 - clhs315;
            lhs(11,0)=DN(2,0)*clhs270 + clhs91;
            lhs(11,1)=DN(2,1)*clhs270 + clhs157;
            lhs(11,2)=DN(2,2)*clhs270 + clhs192;
            lhs(11,3)=clhs212;
            lhs(11,4)=DN(2,0)*clhs204 + clhs234;
            lhs(11,5)=DN(2,1)*clhs204 + clhs251;
            lhs(11,6)=DN(2,2)*clhs204 + clhs264;
            lhs(11,7)=clhs275;
            lhs(11,8)=DN(2,0)*clhs316;
            lhs(11,9)=DN(2,1)*clhs316;
            lhs(11,10)=DN(2,2)*clhs316;
            lhs(11,11)=clhs202*(clhs284 + clhs300 + clhs309);
            lhs(11,12)=DN(2,0)*clhs214 + clhs317;
            lhs(11,13)=DN(2,1)*clhs214 + clhs318;
            lhs(11,14)=DN(2,2)*clhs214 + clhs319;
            lhs(11,15)=clhs320;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs321 + clhs99;
            lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs160;
            lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs194;
            lhs(12,3)=-DN(0,0)*clhs322 + DN(0,0)*clhs323 - clhs213;
            lhs(12,4)=DN(3,0)*clhs38 + DN(3,1)*clhs40 + DN(3,2)*clhs42 + clhs237 + clhs324;
            lhs(12,5)=DN(3,0)*clhs51 + DN(3,1)*clhs53 + DN(3,2)*clhs56 + clhs252;
            lhs(12,6)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs62 + clhs265;
            lhs(12,7)=-DN(1,0)*clhs322 + DN(1,0)*clhs323 - clhs276;
            lhs(12,8)=DN(3,0)*clhs65 + DN(3,1)*clhs67 + DN(3,2)*clhs69 + clhs293 + clhs325;
            lhs(12,9)=DN(3,0)*clhs78 + DN(3,1)*clhs80 + DN(3,2)*clhs83 + clhs303;
            lhs(12,10)=DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs89 + clhs311;
            lhs(12,11)=-DN(2,0)*clhs322 + DN(2,0)*clhs323 - clhs317;
            lhs(12,12)=DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs96 + clhs10*clhs326 + clhs328;
            lhs(12,13)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs110 + clhs330;
            lhs(12,14)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs331;
            lhs(12,15)=DN(3,0)*clhs332;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs119 + DN(3,2)*clhs120 + clhs111;
            lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs165 + clhs333;
            lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs125 + DN(3,2)*clhs127 + clhs196;
            lhs(13,3)=-DN(0,1)*clhs322 + DN(0,1)*clhs323 - clhs215;
            lhs(13,4)=DN(3,0)*clhs40 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs239;
            lhs(13,5)=DN(3,0)*clhs53 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs254 + clhs334;
            lhs(13,6)=DN(3,0)*clhs60 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs266;
            lhs(13,7)=-DN(1,1)*clhs322 + DN(1,1)*clhs323 - clhs277;
            lhs(13,8)=DN(3,0)*clhs67 + DN(3,1)*clhs144 + DN(3,2)*clhs145 + clhs295;
            lhs(13,9)=DN(3,0)*clhs80 + DN(3,1)*clhs147 + DN(3,2)*clhs149 + clhs305 + clhs335;
            lhs(13,10)=DN(3,0)*clhs87 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs312;
            lhs(13,11)=-DN(2,1)*clhs322 + DN(2,1)*clhs323 - clhs318;
            lhs(13,12)=DN(3,0)*clhs94 + DN(3,1)*clhs158 + DN(3,2)*clhs159 + clhs330;
            lhs(13,13)=DN(3,0)*clhs107 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs10*clhs336 + clhs328;
            lhs(13,14)=DN(3,0)*clhs114 + DN(3,1)*clhs167 + DN(3,2)*clhs169 + clhs337;
            lhs(13,15)=DN(3,1)*clhs332;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs120 + DN(3,2)*clhs172 + clhs117;
            lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs123 + DN(3,2)*clhs173 + clhs170;
            lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs127 + DN(3,2)*clhs174 + clhs199 + clhs333;
            lhs(14,3)=-DN(0,2)*clhs322 + DN(0,2)*clhs323 - clhs216;
            lhs(14,4)=DN(3,0)*clhs42 + DN(3,1)*clhs131 + DN(3,2)*clhs176 + clhs240;
            lhs(14,5)=DN(3,0)*clhs56 + DN(3,1)*clhs135 + DN(3,2)*clhs179 + clhs256;
            lhs(14,6)=DN(3,0)*clhs62 + DN(3,1)*clhs141 + DN(3,2)*clhs181 + clhs268 + clhs334;
            lhs(14,7)=-DN(1,2)*clhs322 + DN(1,2)*clhs323 - clhs278;
            lhs(14,8)=DN(3,0)*clhs69 + DN(3,1)*clhs145 + DN(3,2)*clhs185 + clhs296;
            lhs(14,9)=DN(3,0)*clhs83 + DN(3,1)*clhs149 + DN(3,2)*clhs187 + clhs307;
            lhs(14,10)=DN(3,0)*clhs89 + DN(3,1)*clhs155 + DN(3,2)*clhs189 + clhs314 + clhs335;
            lhs(14,11)=-DN(2,2)*clhs322 + DN(2,2)*clhs323 - clhs319;
            lhs(14,12)=DN(3,0)*clhs96 + DN(3,1)*clhs159 + DN(3,2)*clhs193 + clhs331;
            lhs(14,13)=DN(3,0)*clhs110 + DN(3,1)*clhs163 + DN(3,2)*clhs195 + clhs337;
            lhs(14,14)=DN(3,0)*clhs116 + DN(3,1)*clhs169 + DN(3,2)*clhs197 + clhs10*clhs338 + clhs328;
            lhs(14,15)=DN(3,2)*clhs332;
            lhs(15,0)=DN(3,0)*clhs270 + clhs118;
            lhs(15,1)=DN(3,1)*clhs270 + clhs171;
            lhs(15,2)=DN(3,2)*clhs270 + clhs200;
            lhs(15,3)=clhs217;
            lhs(15,4)=DN(3,0)*clhs204 + clhs241;
            lhs(15,5)=DN(3,1)*clhs204 + clhs257;
            lhs(15,6)=DN(3,2)*clhs204 + clhs269;
            lhs(15,7)=clhs279;
            lhs(15,8)=DN(3,0)*clhs209 + clhs297;
            lhs(15,9)=DN(3,1)*clhs209 + clhs308;
            lhs(15,10)=DN(3,2)*clhs209 + clhs315;
            lhs(15,11)=clhs320;
            lhs(15,12)=DN(3,0)*clhs339;
            lhs(15,13)=DN(3,1)*clhs339;
            lhs(15,14)=DN(3,2)*clhs339;
            lhs(15,15)=clhs202*(clhs326 + clhs336 + clhs338);


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

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = rData.IsAir() ? rData.VolumeError : -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

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
const double crhs9 =             crhs4 + crhs8 - volume_error_ratio;
const double crhs10 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
const double crhs11 =             crhs9*(crhs10*h/stab_c1 + mu);
const double crhs12 =             1.0/(K_darcy + crhs10/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs13 =             crhs12*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs3 + crhs7);
const double crhs14 =             K_darcy*N[0];
const double crhs15 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
const double crhs16 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs17 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs18 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs19 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs8);
const double crhs20 =             crhs12*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs16 + crhs17 + crhs18 + crhs19);
const double crhs21 =             K_darcy*N[1];
const double crhs22 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
const double crhs23 =             K_darcy*N[2];
const double crhs24 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs7 + crhs13*crhs14 - crhs13*crhs15;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs11 - DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 - N[0]*crhs18 - N[0]*crhs19 + crhs14*crhs20 - crhs15*crhs20;
            rhs[2]=-DN(0,0)*crhs13 - DN(0,1)*crhs20 - N[0]*crhs9;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs7 + crhs13*crhs21 - crhs13*crhs22;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs11 - DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 - N[1]*crhs18 - N[1]*crhs19 + crhs20*crhs21 - crhs20*crhs22;
            rhs[5]=-DN(1,0)*crhs13 - DN(1,1)*crhs20 - N[1]*crhs9;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs7 + crhs13*crhs23 - crhs13*crhs24;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs11 - DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 - N[2]*crhs18 - N[2]*crhs19 + crhs20*crhs23 - crhs20*crhs24;
            rhs[8]=-DN(2,0)*crhs13 - DN(2,1)*crhs20 - N[2]*crhs9;


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

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = rData.IsAir() ? rData.VolumeError : -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs7 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs8 =             rho*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs7*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs9 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs10 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs11 =             crhs10 + crhs4 + crhs9 - volume_error_ratio;
const double crhs12 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2) + pow(crhs7, 2));
const double crhs13 =             crhs11*(crhs12*h/stab_c1 + mu);
const double crhs14 =             1.0/(K_darcy + crhs12/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs15 =             crhs14*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs2 + crhs3 + crhs8);
const double crhs16 =             K_darcy*N[0];
const double crhs17 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6 + DN(0,2)*crhs7);
const double crhs18 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs19 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs20 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs21 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*crhs9 + crhs7*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs22 =             crhs14*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs18 + crhs19 + crhs20 + crhs21);
const double crhs23 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs24 =             K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs25 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs26 =             rho*(crhs10*crhs7 + crhs5*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs6*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs27 =             crhs14*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs23 + crhs24 + crhs25 + crhs26);
const double crhs28 =             K_darcy*N[1];
const double crhs29 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6 + DN(1,2)*crhs7);
const double crhs30 =             K_darcy*N[2];
const double crhs31 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6 + DN(2,2)*crhs7);
const double crhs32 =             K_darcy*N[3];
const double crhs33 =             rho*(DN(3,0)*crhs5 + DN(3,1)*crhs6 + DN(3,2)*crhs7);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs13 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs8 + crhs15*crhs16 - crhs15*crhs17;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs13 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs18 - N[0]*crhs19 - N[0]*crhs20 - N[0]*crhs21 + crhs16*crhs22 - crhs17*crhs22;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs13 - DN(0,2)*stress[2] + N[0]*crhs23 - N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 + crhs16*crhs27 - crhs17*crhs27;
            rhs[3]=-DN(0,0)*crhs15 - DN(0,1)*crhs22 - DN(0,2)*crhs27 - N[0]*crhs11;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs13 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs8 + crhs15*crhs28 - crhs15*crhs29;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs13 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs18 - N[1]*crhs19 - N[1]*crhs20 - N[1]*crhs21 + crhs22*crhs28 - crhs22*crhs29;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs13 - DN(1,2)*stress[2] + N[1]*crhs23 - N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 + crhs27*crhs28 - crhs27*crhs29;
            rhs[7]=-DN(1,0)*crhs15 - DN(1,1)*crhs22 - DN(1,2)*crhs27 - N[1]*crhs11;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs13 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs8 + crhs15*crhs30 - crhs15*crhs31;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs13 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs18 - N[2]*crhs19 - N[2]*crhs20 - N[2]*crhs21 + crhs22*crhs30 - crhs22*crhs31;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs13 - DN(2,2)*stress[2] + N[2]*crhs23 - N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 + crhs27*crhs30 - crhs27*crhs31;
            rhs[11]=-DN(2,0)*crhs15 - DN(2,1)*crhs22 - DN(2,2)*crhs27 - N[2]*crhs11;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs13 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs2 - N[3]*crhs3 - N[3]*crhs8 + crhs15*crhs32 - crhs15*crhs33;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs13 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs18 - N[3]*crhs19 - N[3]*crhs20 - N[3]*crhs21 + crhs22*crhs32 - crhs22*crhs33;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs13 - DN(3,2)*stress[2] + N[3]*crhs23 - N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 + crhs27*crhs32 - crhs27*crhs33;
            rhs[15]=-DN(3,0)*crhs15 - DN(3,1)*crhs22 - DN(3,2)*crhs27 - N[3]*crhs11;


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
    //energy_check_cut elements

    // Mass correction term
        double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = rData.IsAir() ? rData.VolumeError : -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

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
const double crhs_ee2 =             crhs_ee0 + crhs_ee1 - volume_error_ratio;
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
    //energy_check_cut elements


    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = rData.IsAir() ? rData.VolumeError : -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

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
const double crhs_ee3 =             crhs_ee0 + crhs_ee1 + crhs_ee2 - volume_error_ratio;
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
void TwoFluidNavierStokes<TElementData>::CalculateOnIntegrationPoints(
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