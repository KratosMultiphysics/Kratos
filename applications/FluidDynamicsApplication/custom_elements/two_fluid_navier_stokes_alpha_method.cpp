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

#include "two_fluid_navier_stokes_alpha_method.h"
#include "custom_utilities/two_fluid_navier_stokes_alpha_method_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::TwoFluidNavierStokesAlphaMethod(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokesAlphaMethod<TElementData>::~TwoFluidNavierStokesAlphaMethod() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethod<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethod>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesAlphaMethod<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesAlphaMethod>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::CalculateLocalSystem(
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
            std::vector<array_1d<double,3>> int_normals_neg;                                 // interface normal vector based on the negative side
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
        KRATOS_ERROR << "TwoFluidNavierStokesAlphaMethod is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokesAlphaMethod<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
std::string TwoFluidNavierStokesAlphaMethod<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokesAlphaMethod" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::PrintInfo(
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);

    // Calculate current material response
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::UpdateIntegrationPointData(
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
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodData<2, 3>& rData) const
{
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const BoundedMatrix<double,3,2> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += rDNDX(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += rDNDX(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += rDNDX(i,0)*velocity_alpha(i,1) + rDNDX(i,1)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesAlphaMethodData<3, 4>& rData) const
{
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const BoundedMatrix<double,4,3> velocity_alpha = rData.Velocity_OldStep1+ alpha_f*(rData.Velocity-rData.Velocity_OldStep1);
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(6);
    for (unsigned int i = 0; i < 4; i++) {
        r_strain_rate[0] += rDNDX(i,0)*velocity_alpha(i,0);
        r_strain_rate[1] += rDNDX(i,1)*velocity_alpha(i,1);
        r_strain_rate[2] += rDNDX(i,2)*velocity_alpha(i,2);
        r_strain_rate[3] += rDNDX(i,0)*velocity_alpha(i,1) + rDNDX(i,1)*velocity_alpha(i,0);
        r_strain_rate[4] += rDNDX(i,1)*velocity_alpha(i,2) + rDNDX(i,2)*velocity_alpha(i,1);
        r_strain_rate[5] += rDNDX(i,0)*velocity_alpha(i,2) + rDNDX(i,2)*velocity_alpha(i,0);
    }
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const double dyn_tau = rData.DynamicTau;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double alpha_f=1/(1+rData.MaxSprectraRadius);
    const BoundedMatrix<double,3,2> vconv =(vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

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
const double clhs1 =             mu*stab_c1;
const double clhs2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs4 =             rho*sqrt(pow(clhs2, 2) + pow(clhs3, 2));
const double clhs5 =             clhs1/stab_c2 + clhs4*h;
const double clhs6 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs7 =             1.0/(max_spectral_radius + 1);
const double clhs8 =             DN(0,0)*clhs7;
const double clhs9 =             C(0,2)*DN(0,0);
const double clhs10 =             C(2,2)*DN(0,1) + clhs9;
const double clhs11 =             DN(0,1)*clhs7;
const double clhs12 =             DN(0,0)*clhs2 + DN(0,1)*clhs3;
const double clhs13 =             clhs12*rho;
const double clhs14 =             N[0]*clhs7;
const double clhs15 =             1.0/dt;
const double clhs16 =             clhs7*(0.5*max_spectral_radius - 1.5);
const double clhs17 =             0.5*max_spectral_radius - 1.5;
const double clhs18 =             clhs15*clhs17/(-clhs16 - clhs7 + 0.5);
const double clhs19 =             -N[0]*clhs18 + clhs12;
const double clhs20 =             clhs15*rho;
const double clhs21 =             1.0/(clhs1/pow(h, 2) + clhs20*dyn_tau + clhs4*stab_c2/h);
const double clhs22 =             clhs21*pow(rho, 2);
const double clhs23 =             clhs22*clhs7;
const double clhs24 =             clhs12*clhs23;
const double clhs25 =             clhs17*clhs20/(clhs16 + clhs7 - 0.5);
const double clhs26 =             clhs25*clhs7;
const double clhs27 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs28 =             clhs22*clhs27;
const double clhs29 =             clhs14*clhs28;
const double clhs30 =             pow(N[0], 2)*clhs26 + clhs13*clhs14 + clhs19*clhs24 + clhs19*clhs29;
const double clhs31 =             C(0,1)*DN(0,1) + clhs9;
const double clhs32 =             C(1,2)*DN(0,1);
const double clhs33 =             C(2,2)*DN(0,0) + clhs32;
const double clhs34 =             DN(0,0)*clhs5;
const double clhs35 =             DN(0,1)*clhs34;
const double clhs36 =             clhs21*rho;
const double clhs37 =             clhs27*clhs36;
const double clhs38 =             clhs13*clhs21;
const double clhs39 =             N[0]*clhs37 - N[0] + clhs38;
const double clhs40 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs41 =             C(0,2)*DN(1,0);
const double clhs42 =             C(2,2)*DN(1,1) + clhs41;
const double clhs43 =             DN(0,0)*DN(1,0);
const double clhs44 =             clhs14*clhs25;
const double clhs45 =             N[1]*clhs44;
const double clhs46 =             clhs43*clhs5 + clhs45;
const double clhs47 =             DN(1,0)*clhs2 + DN(1,1)*clhs3;
const double clhs48 =             clhs14*rho;
const double clhs49 =             -N[1]*clhs18 + clhs47;
const double clhs50 =             clhs24*clhs49 + clhs29*clhs49 + clhs47*clhs48;
const double clhs51 =             C(0,1)*DN(1,1) + clhs41;
const double clhs52 =             C(1,2)*DN(1,1);
const double clhs53 =             C(2,2)*DN(1,0) + clhs52;
const double clhs54 =             DN(1,1)*clhs34;
const double clhs55 =             DN(0,0)*N[1];
const double clhs56 =             DN(1,0)*N[0];
const double clhs57 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs58 =             C(0,2)*DN(2,0);
const double clhs59 =             C(2,2)*DN(2,1) + clhs58;
const double clhs60 =             DN(0,0)*DN(2,0);
const double clhs61 =             N[2]*clhs44;
const double clhs62 =             clhs5*clhs60 + clhs61;
const double clhs63 =             DN(2,0)*clhs2 + DN(2,1)*clhs3;
const double clhs64 =             -N[2]*clhs18 + clhs63;
const double clhs65 =             clhs24*clhs64 + clhs29*clhs64 + clhs48*clhs63;
const double clhs66 =             C(0,1)*DN(2,1) + clhs58;
const double clhs67 =             C(1,2)*DN(2,1);
const double clhs68 =             C(2,2)*DN(2,0) + clhs67;
const double clhs69 =             DN(2,1)*clhs34;
const double clhs70 =             DN(0,0)*N[2];
const double clhs71 =             DN(2,0)*N[0];
const double clhs72 =             C(0,1)*DN(0,0) + clhs32;
const double clhs73 =             pow(DN(0,1), 2);
const double clhs74 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs75 =             C(0,1)*DN(1,0) + clhs52;
const double clhs76 =             DN(0,1)*clhs5;
const double clhs77 =             DN(1,0)*clhs76;
const double clhs78 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs79 =             DN(0,1)*DN(1,1);
const double clhs80 =             clhs45 + clhs5*clhs79;
const double clhs81 =             DN(0,1)*N[1];
const double clhs82 =             DN(1,1)*N[0];
const double clhs83 =             C(0,1)*DN(2,0) + clhs67;
const double clhs84 =             DN(2,0)*clhs76;
const double clhs85 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs86 =             DN(0,1)*DN(2,1);
const double clhs87 =             clhs5*clhs86 + clhs61;
const double clhs88 =             DN(0,1)*N[2];
const double clhs89 =             DN(2,1)*N[0];
const double clhs90 =             clhs19*clhs36;
const double clhs91 =             N[0] + clhs90;
const double clhs92 =             clhs36*clhs49;
const double clhs93 =             clhs21*(clhs43 + clhs79);
const double clhs94 =             clhs36*clhs64;
const double clhs95 =             clhs21*(clhs60 + clhs86);
const double clhs96 =             DN(1,0)*clhs7;
const double clhs97 =             DN(1,1)*clhs7;
const double clhs98 =             N[1]*clhs7;
const double clhs99 =             clhs23*clhs47;
const double clhs100 =             clhs28*clhs98;
const double clhs101 =             clhs100*clhs19 + clhs13*clhs98 + clhs19*clhs99;
const double clhs102 =             clhs36*clhs47;
const double clhs103 =             pow(DN(1,0), 2);
const double clhs104 =             clhs98*rho;
const double clhs105 =             pow(N[1], 2)*clhs26 + clhs100*clhs49 + clhs104*clhs47 + clhs49*clhs99;
const double clhs106 =             DN(1,0)*clhs5;
const double clhs107 =             DN(1,1)*clhs106;
const double clhs108 =             N[1]*clhs37 - N[1] + clhs102;
const double clhs109 =             DN(1,0)*DN(2,0);
const double clhs110 =             N[2]*clhs25*clhs98;
const double clhs111 =             clhs109*clhs5 + clhs110;
const double clhs112 =             clhs100*clhs64 + clhs104*clhs63 + clhs64*clhs99;
const double clhs113 =             DN(2,1)*clhs106;
const double clhs114 =             DN(1,0)*N[2];
const double clhs115 =             DN(2,0)*N[1];
const double clhs116 =             pow(DN(1,1), 2);
const double clhs117 =             DN(2,0)*clhs5;
const double clhs118 =             DN(1,1)*clhs117;
const double clhs119 =             DN(1,1)*DN(2,1);
const double clhs120 =             clhs110 + clhs119*clhs5;
const double clhs121 =             DN(1,1)*N[2];
const double clhs122 =             DN(2,1)*N[1];
const double clhs123 =             N[1] + clhs92;
const double clhs124 =             clhs21*(clhs109 + clhs119);
const double clhs125 =             DN(2,0)*clhs7;
const double clhs126 =             DN(2,1)*clhs7;
const double clhs127 =             N[2]*clhs7;
const double clhs128 =             clhs23*clhs63;
const double clhs129 =             clhs127*clhs28;
const double clhs130 =             clhs127*clhs13 + clhs128*clhs19 + clhs129*clhs19;
const double clhs131 =             clhs36*clhs63;
const double clhs132 =             clhs127*rho;
const double clhs133 =             clhs128*clhs49 + clhs129*clhs49 + clhs132*clhs47;
const double clhs134 =             pow(DN(2,0), 2);
const double clhs135 =             pow(N[2], 2)*clhs26 + clhs128*clhs64 + clhs129*clhs64 + clhs132*clhs63;
const double clhs136 =             DN(2,1)*clhs117;
const double clhs137 =             N[2]*clhs37 - N[2] + clhs131;
const double clhs138 =             pow(DN(2,1), 2);
const double clhs139 =             N[2] + clhs94;
            lhs(0,0)=clhs0*clhs5 + clhs10*clhs11 + clhs30 + clhs6*clhs8;
            lhs(0,1)=clhs11*clhs33 + clhs31*clhs8 + clhs35;
            lhs(0,2)=DN(0,0)*clhs39;
            lhs(0,3)=clhs11*clhs42 + clhs40*clhs8 + clhs46 + clhs50;
            lhs(0,4)=clhs11*clhs53 + clhs51*clhs8 + clhs54;
            lhs(0,5)=DN(1,0)*clhs38 + clhs37*clhs56 - clhs55;
            lhs(0,6)=clhs11*clhs59 + clhs57*clhs8 + clhs62 + clhs65;
            lhs(0,7)=clhs11*clhs68 + clhs66*clhs8 + clhs69;
            lhs(0,8)=DN(2,0)*clhs38 + clhs37*clhs71 - clhs70;
            lhs(1,0)=clhs10*clhs8 + clhs11*clhs72 + clhs35;
            lhs(1,1)=clhs11*clhs74 + clhs30 + clhs33*clhs8 + clhs5*clhs73;
            lhs(1,2)=DN(0,1)*clhs39;
            lhs(1,3)=clhs11*clhs75 + clhs42*clhs8 + clhs77;
            lhs(1,4)=clhs11*clhs78 + clhs50 + clhs53*clhs8 + clhs80;
            lhs(1,5)=DN(1,1)*clhs38 + clhs37*clhs82 - clhs81;
            lhs(1,6)=clhs11*clhs83 + clhs59*clhs8 + clhs84;
            lhs(1,7)=clhs11*clhs85 + clhs65 + clhs68*clhs8 + clhs87;
            lhs(1,8)=DN(2,1)*clhs38 + clhs37*clhs89 - clhs88;
            lhs(2,0)=clhs8*clhs91;
            lhs(2,1)=clhs11*clhs91;
            lhs(2,2)=clhs21*(clhs0 + clhs73);
            lhs(2,3)=clhs7*(DN(0,0)*clhs92 + clhs56);
            lhs(2,4)=clhs7*(DN(0,1)*clhs92 + clhs82);
            lhs(2,5)=clhs93;
            lhs(2,6)=clhs7*(DN(0,0)*clhs94 + clhs71);
            lhs(2,7)=clhs7*(DN(0,1)*clhs94 + clhs89);
            lhs(2,8)=clhs95;
            lhs(3,0)=clhs10*clhs97 + clhs101 + clhs46 + clhs6*clhs96;
            lhs(3,1)=clhs31*clhs96 + clhs33*clhs97 + clhs77;
            lhs(3,2)=DN(0,0)*clhs102 + clhs37*clhs55 - clhs56;
            lhs(3,3)=clhs103*clhs5 + clhs105 + clhs40*clhs96 + clhs42*clhs97;
            lhs(3,4)=clhs107 + clhs51*clhs96 + clhs53*clhs97;
            lhs(3,5)=DN(1,0)*clhs108;
            lhs(3,6)=clhs111 + clhs112 + clhs57*clhs96 + clhs59*clhs97;
            lhs(3,7)=clhs113 + clhs66*clhs96 + clhs68*clhs97;
            lhs(3,8)=DN(2,0)*clhs102 - clhs114 + clhs115*clhs37;
            lhs(4,0)=clhs10*clhs96 + clhs54 + clhs72*clhs97;
            lhs(4,1)=clhs101 + clhs33*clhs96 + clhs74*clhs97 + clhs80;
            lhs(4,2)=DN(0,1)*clhs102 + clhs37*clhs81 - clhs82;
            lhs(4,3)=clhs107 + clhs42*clhs96 + clhs75*clhs97;
            lhs(4,4)=clhs105 + clhs116*clhs5 + clhs53*clhs96 + clhs78*clhs97;
            lhs(4,5)=DN(1,1)*clhs108;
            lhs(4,6)=clhs118 + clhs59*clhs96 + clhs83*clhs97;
            lhs(4,7)=clhs112 + clhs120 + clhs68*clhs96 + clhs85*clhs97;
            lhs(4,8)=DN(2,1)*clhs102 - clhs121 + clhs122*clhs37;
            lhs(5,0)=clhs7*(DN(1,0)*clhs90 + clhs55);
            lhs(5,1)=clhs7*(DN(1,1)*clhs90 + clhs81);
            lhs(5,2)=clhs93;
            lhs(5,3)=clhs123*clhs96;
            lhs(5,4)=clhs123*clhs97;
            lhs(5,5)=clhs21*(clhs103 + clhs116);
            lhs(5,6)=clhs7*(DN(1,0)*clhs94 + clhs115);
            lhs(5,7)=clhs7*(DN(1,1)*clhs94 + clhs122);
            lhs(5,8)=clhs124;
            lhs(6,0)=clhs10*clhs126 + clhs125*clhs6 + clhs130 + clhs62;
            lhs(6,1)=clhs125*clhs31 + clhs126*clhs33 + clhs84;
            lhs(6,2)=DN(0,0)*clhs131 + clhs37*clhs70 - clhs71;
            lhs(6,3)=clhs111 + clhs125*clhs40 + clhs126*clhs42 + clhs133;
            lhs(6,4)=clhs118 + clhs125*clhs51 + clhs126*clhs53;
            lhs(6,5)=DN(1,0)*clhs131 + clhs114*clhs37 - clhs115;
            lhs(6,6)=clhs125*clhs57 + clhs126*clhs59 + clhs134*clhs5 + clhs135;
            lhs(6,7)=clhs125*clhs66 + clhs126*clhs68 + clhs136;
            lhs(6,8)=DN(2,0)*clhs137;
            lhs(7,0)=clhs10*clhs125 + clhs126*clhs72 + clhs69;
            lhs(7,1)=clhs125*clhs33 + clhs126*clhs74 + clhs130 + clhs87;
            lhs(7,2)=DN(0,1)*clhs131 + clhs37*clhs88 - clhs89;
            lhs(7,3)=clhs113 + clhs125*clhs42 + clhs126*clhs75;
            lhs(7,4)=clhs120 + clhs125*clhs53 + clhs126*clhs78 + clhs133;
            lhs(7,5)=DN(1,1)*clhs131 + clhs121*clhs37 - clhs122;
            lhs(7,6)=clhs125*clhs59 + clhs126*clhs83 + clhs136;
            lhs(7,7)=clhs125*clhs68 + clhs126*clhs85 + clhs135 + clhs138*clhs5;
            lhs(7,8)=DN(2,1)*clhs137;
            lhs(8,0)=clhs7*(DN(2,0)*clhs90 + clhs70);
            lhs(8,1)=clhs7*(DN(2,1)*clhs90 + clhs88);
            lhs(8,2)=clhs95;
            lhs(8,3)=clhs7*(DN(2,0)*clhs92 + clhs114);
            lhs(8,4)=clhs7*(DN(2,1)*clhs92 + clhs121);
            lhs(8,5)=clhs124;
            lhs(8,6)=clhs125*clhs139;
            lhs(8,7)=clhs126*clhs139;
            lhs(8,8)=clhs21*(clhs134 + clhs138);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

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
const double clhs1 =             mu*stab_c1;
const double clhs2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs5 =             rho*sqrt(pow(clhs2, 2) + pow(clhs3, 2) + pow(clhs4, 2));
const double clhs6 =             clhs1/stab_c2 + clhs5*h;
const double clhs7 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs8 =             1.0/(max_spectral_radius + 1);
const double clhs9 =             DN(0,0)*clhs8;
const double clhs10 =             C(0,3)*DN(0,0);
const double clhs11 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs10;
const double clhs12 =             DN(0,1)*clhs8;
const double clhs13 =             C(0,5)*DN(0,0);
const double clhs14 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs13;
const double clhs15 =             DN(0,2)*clhs8;
const double clhs16 =             DN(0,0)*clhs2 + DN(0,1)*clhs3 + DN(0,2)*clhs4;
const double clhs17 =             clhs16*rho;
const double clhs18 =             N[0]*clhs8;
const double clhs19 =             1.0/dt;
const double clhs20 =             clhs8*(0.5*max_spectral_radius - 1.5);
const double clhs21 =             0.5*max_spectral_radius - 1.5;
const double clhs22 =             clhs19*clhs21/(-clhs20 - clhs8 + 0.5);
const double clhs23 =             -N[0]*clhs22 + clhs16;
const double clhs24 =             clhs19*rho;
const double clhs25 =             1.0/(clhs1/pow(h, 2) + clhs24*dyn_tau + clhs5*stab_c2/h);
const double clhs26 =             clhs25*pow(rho, 2);
const double clhs27 =             clhs26*clhs8;
const double clhs28 =             clhs16*clhs27;
const double clhs29 =             clhs21*clhs24/(clhs20 + clhs8 - 0.5);
const double clhs30 =             clhs29*clhs8;
const double clhs31 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs32 =             clhs26*clhs31;
const double clhs33 =             clhs18*clhs32;
const double clhs34 =             pow(N[0], 2)*clhs30 + clhs17*clhs18 + clhs23*clhs28 + clhs23*clhs33;
const double clhs35 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs10;
const double clhs36 =             C(1,3)*DN(0,1);
const double clhs37 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs36;
const double clhs38 =             C(3,5)*DN(0,0);
const double clhs39 =             C(4,5)*DN(0,2);
const double clhs40 =             C(1,5)*DN(0,1) + clhs38 + clhs39;
const double clhs41 =             DN(0,0)*clhs6;
const double clhs42 =             DN(0,1)*clhs41;
const double clhs43 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs13;
const double clhs44 =             C(3,4)*DN(0,1);
const double clhs45 =             C(2,3)*DN(0,2) + clhs38 + clhs44;
const double clhs46 =             C(2,5)*DN(0,2);
const double clhs47 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs46;
const double clhs48 =             DN(0,2)*clhs41;
const double clhs49 =             clhs25*rho;
const double clhs50 =             clhs31*clhs49;
const double clhs51 =             clhs17*clhs25;
const double clhs52 =             N[0]*clhs50 - N[0] + clhs51;
const double clhs53 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs54 =             C(0,3)*DN(1,0);
const double clhs55 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs54;
const double clhs56 =             C(0,5)*DN(1,0);
const double clhs57 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs56;
const double clhs58 =             DN(0,0)*DN(1,0);
const double clhs59 =             clhs18*clhs29;
const double clhs60 =             N[1]*clhs59;
const double clhs61 =             clhs58*clhs6 + clhs60;
const double clhs62 =             DN(1,0)*clhs2 + DN(1,1)*clhs3 + DN(1,2)*clhs4;
const double clhs63 =             clhs18*rho;
const double clhs64 =             -N[1]*clhs22 + clhs62;
const double clhs65 =             clhs28*clhs64 + clhs33*clhs64 + clhs62*clhs63;
const double clhs66 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs54;
const double clhs67 =             C(1,3)*DN(1,1);
const double clhs68 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs67;
const double clhs69 =             C(3,5)*DN(1,0);
const double clhs70 =             C(4,5)*DN(1,2);
const double clhs71 =             C(1,5)*DN(1,1) + clhs69 + clhs70;
const double clhs72 =             DN(1,1)*clhs41;
const double clhs73 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs56;
const double clhs74 =             C(3,4)*DN(1,1);
const double clhs75 =             C(2,3)*DN(1,2) + clhs69 + clhs74;
const double clhs76 =             C(2,5)*DN(1,2);
const double clhs77 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs76;
const double clhs78 =             DN(1,2)*clhs41;
const double clhs79 =             DN(0,0)*N[1];
const double clhs80 =             DN(1,0)*N[0];
const double clhs81 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs82 =             C(0,3)*DN(2,0);
const double clhs83 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs82;
const double clhs84 =             C(0,5)*DN(2,0);
const double clhs85 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs84;
const double clhs86 =             DN(0,0)*DN(2,0);
const double clhs87 =             N[2]*clhs59;
const double clhs88 =             clhs6*clhs86 + clhs87;
const double clhs89 =             DN(2,0)*clhs2 + DN(2,1)*clhs3 + DN(2,2)*clhs4;
const double clhs90 =             -N[2]*clhs22 + clhs89;
const double clhs91 =             clhs28*clhs90 + clhs33*clhs90 + clhs63*clhs89;
const double clhs92 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs82;
const double clhs93 =             C(1,3)*DN(2,1);
const double clhs94 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs93;
const double clhs95 =             C(3,5)*DN(2,0);
const double clhs96 =             C(4,5)*DN(2,2);
const double clhs97 =             C(1,5)*DN(2,1) + clhs95 + clhs96;
const double clhs98 =             DN(2,1)*clhs41;
const double clhs99 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs84;
const double clhs100 =             C(3,4)*DN(2,1);
const double clhs101 =             C(2,3)*DN(2,2) + clhs100 + clhs95;
const double clhs102 =             C(2,5)*DN(2,2);
const double clhs103 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs102;
const double clhs104 =             DN(2,2)*clhs41;
const double clhs105 =             DN(0,0)*N[2];
const double clhs106 =             DN(2,0)*N[0];
const double clhs107 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs108 =             C(0,3)*DN(3,0);
const double clhs109 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs108;
const double clhs110 =             C(0,5)*DN(3,0);
const double clhs111 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs110;
const double clhs112 =             DN(0,0)*DN(3,0);
const double clhs113 =             N[3]*clhs59;
const double clhs114 =             clhs112*clhs6 + clhs113;
const double clhs115 =             DN(3,0)*clhs2 + DN(3,1)*clhs3 + DN(3,2)*clhs4;
const double clhs116 =             -N[3]*clhs22 + clhs115;
const double clhs117 =             clhs115*clhs63 + clhs116*clhs28 + clhs116*clhs33;
const double clhs118 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs108;
const double clhs119 =             C(1,3)*DN(3,1);
const double clhs120 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs119;
const double clhs121 =             C(3,5)*DN(3,0);
const double clhs122 =             C(4,5)*DN(3,2);
const double clhs123 =             C(1,5)*DN(3,1) + clhs121 + clhs122;
const double clhs124 =             DN(3,1)*clhs41;
const double clhs125 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs110;
const double clhs126 =             C(3,4)*DN(3,1);
const double clhs127 =             C(2,3)*DN(3,2) + clhs121 + clhs126;
const double clhs128 =             C(2,5)*DN(3,2);
const double clhs129 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs128;
const double clhs130 =             DN(3,2)*clhs41;
const double clhs131 =             DN(0,0)*N[3];
const double clhs132 =             DN(3,0)*N[0];
const double clhs133 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs36;
const double clhs134 =             C(0,4)*DN(0,0) + clhs39 + clhs44;
const double clhs135 =             pow(DN(0,1), 2);
const double clhs136 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs137 =             C(1,4)*DN(0,1);
const double clhs138 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs137;
const double clhs139 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs137;
const double clhs140 =             C(2,4)*DN(0,2);
const double clhs141 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs140;
const double clhs142 =             DN(0,1)*clhs6;
const double clhs143 =             DN(0,2)*clhs142;
const double clhs144 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs67;
const double clhs145 =             C(0,4)*DN(1,0) + clhs70 + clhs74;
const double clhs146 =             DN(1,0)*clhs142;
const double clhs147 =             DN(0,1)*DN(1,1);
const double clhs148 =             clhs147*clhs6;
const double clhs149 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs150 =             C(1,4)*DN(1,1);
const double clhs151 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs150;
const double clhs152 =             clhs60 + clhs65;
const double clhs153 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs150;
const double clhs154 =             C(2,4)*DN(1,2);
const double clhs155 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs154;
const double clhs156 =             DN(1,2)*clhs142;
const double clhs157 =             DN(0,1)*N[1];
const double clhs158 =             DN(1,1)*N[0];
const double clhs159 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs93;
const double clhs160 =             C(0,4)*DN(2,0) + clhs100 + clhs96;
const double clhs161 =             DN(2,0)*clhs142;
const double clhs162 =             DN(0,1)*DN(2,1);
const double clhs163 =             clhs162*clhs6;
const double clhs164 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs165 =             C(1,4)*DN(2,1);
const double clhs166 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs165;
const double clhs167 =             clhs87 + clhs91;
const double clhs168 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs165;
const double clhs169 =             C(2,4)*DN(2,2);
const double clhs170 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs169;
const double clhs171 =             DN(2,2)*clhs142;
const double clhs172 =             DN(0,1)*N[2];
const double clhs173 =             DN(2,1)*N[0];
const double clhs174 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs119;
const double clhs175 =             C(0,4)*DN(3,0) + clhs122 + clhs126;
const double clhs176 =             DN(3,0)*clhs142;
const double clhs177 =             DN(0,1)*DN(3,1);
const double clhs178 =             clhs177*clhs6;
const double clhs179 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs180 =             C(1,4)*DN(3,1);
const double clhs181 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs180;
const double clhs182 =             clhs113 + clhs117;
const double clhs183 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs180;
const double clhs184 =             C(2,4)*DN(3,2);
const double clhs185 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs184;
const double clhs186 =             DN(3,2)*clhs142;
const double clhs187 =             DN(0,1)*N[3];
const double clhs188 =             DN(3,1)*N[0];
const double clhs189 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs46;
const double clhs190 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs140;
const double clhs191 =             pow(DN(0,2), 2);
const double clhs192 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs193 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs76;
const double clhs194 =             DN(0,2)*clhs6;
const double clhs195 =             DN(1,0)*clhs194;
const double clhs196 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs154;
const double clhs197 =             DN(1,1)*clhs194;
const double clhs198 =             DN(0,2)*DN(1,2);
const double clhs199 =             clhs198*clhs6;
const double clhs200 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs201 =             DN(0,2)*N[1];
const double clhs202 =             DN(1,2)*N[0];
const double clhs203 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs102;
const double clhs204 =             DN(2,0)*clhs194;
const double clhs205 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs169;
const double clhs206 =             DN(2,1)*clhs194;
const double clhs207 =             DN(0,2)*DN(2,2);
const double clhs208 =             clhs207*clhs6;
const double clhs209 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs210 =             DN(0,2)*N[2];
const double clhs211 =             DN(2,2)*N[0];
const double clhs212 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs128;
const double clhs213 =             DN(3,0)*clhs194;
const double clhs214 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs184;
const double clhs215 =             DN(3,1)*clhs194;
const double clhs216 =             DN(0,2)*DN(3,2);
const double clhs217 =             clhs216*clhs6;
const double clhs218 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs219 =             DN(0,2)*N[3];
const double clhs220 =             DN(3,2)*N[0];
const double clhs221 =             clhs23*clhs49;
const double clhs222 =             N[0] + clhs221;
const double clhs223 =             clhs49*clhs64;
const double clhs224 =             clhs25*(clhs147 + clhs198 + clhs58);
const double clhs225 =             clhs49*clhs90;
const double clhs226 =             clhs25*(clhs162 + clhs207 + clhs86);
const double clhs227 =             clhs116*clhs49;
const double clhs228 =             clhs25*(clhs112 + clhs177 + clhs216);
const double clhs229 =             DN(1,0)*clhs8;
const double clhs230 =             DN(1,1)*clhs8;
const double clhs231 =             DN(1,2)*clhs8;
const double clhs232 =             N[1]*clhs8;
const double clhs233 =             clhs27*clhs62;
const double clhs234 =             clhs232*clhs32;
const double clhs235 =             clhs17*clhs232 + clhs23*clhs233 + clhs23*clhs234;
const double clhs236 =             clhs49*clhs62;
const double clhs237 =             pow(DN(1,0), 2);
const double clhs238 =             clhs232*rho;
const double clhs239 =             pow(N[1], 2)*clhs30 + clhs233*clhs64 + clhs234*clhs64 + clhs238*clhs62;
const double clhs240 =             DN(1,0)*clhs6;
const double clhs241 =             DN(1,1)*clhs240;
const double clhs242 =             DN(1,2)*clhs240;
const double clhs243 =             N[1]*clhs50 - N[1] + clhs236;
const double clhs244 =             DN(1,0)*DN(2,0);
const double clhs245 =             clhs232*clhs29;
const double clhs246 =             N[2]*clhs245;
const double clhs247 =             clhs244*clhs6 + clhs246;
const double clhs248 =             clhs233*clhs90 + clhs234*clhs90 + clhs238*clhs89;
const double clhs249 =             DN(2,1)*clhs240;
const double clhs250 =             DN(2,2)*clhs240;
const double clhs251 =             DN(1,0)*N[2];
const double clhs252 =             DN(2,0)*N[1];
const double clhs253 =             DN(1,0)*DN(3,0);
const double clhs254 =             N[3]*clhs245;
const double clhs255 =             clhs253*clhs6 + clhs254;
const double clhs256 =             clhs115*clhs238 + clhs116*clhs233 + clhs116*clhs234;
const double clhs257 =             DN(3,1)*clhs240;
const double clhs258 =             DN(3,2)*clhs240;
const double clhs259 =             DN(1,0)*N[3];
const double clhs260 =             DN(3,0)*N[1];
const double clhs261 =             clhs235 + clhs60;
const double clhs262 =             pow(DN(1,1), 2);
const double clhs263 =             DN(1,1)*clhs6;
const double clhs264 =             DN(1,2)*clhs263;
const double clhs265 =             DN(2,0)*clhs263;
const double clhs266 =             DN(1,1)*DN(2,1);
const double clhs267 =             clhs266*clhs6;
const double clhs268 =             clhs246 + clhs248;
const double clhs269 =             DN(2,2)*clhs263;
const double clhs270 =             DN(1,1)*N[2];
const double clhs271 =             DN(2,1)*N[1];
const double clhs272 =             DN(3,0)*clhs263;
const double clhs273 =             DN(1,1)*DN(3,1);
const double clhs274 =             clhs273*clhs6;
const double clhs275 =             clhs254 + clhs256;
const double clhs276 =             DN(3,2)*clhs263;
const double clhs277 =             DN(1,1)*N[3];
const double clhs278 =             DN(3,1)*N[1];
const double clhs279 =             pow(DN(1,2), 2);
const double clhs280 =             DN(1,2)*clhs6;
const double clhs281 =             DN(2,0)*clhs280;
const double clhs282 =             DN(2,1)*clhs280;
const double clhs283 =             DN(1,2)*DN(2,2);
const double clhs284 =             clhs283*clhs6;
const double clhs285 =             DN(1,2)*N[2];
const double clhs286 =             DN(2,2)*N[1];
const double clhs287 =             DN(3,0)*clhs280;
const double clhs288 =             DN(3,1)*clhs280;
const double clhs289 =             DN(1,2)*DN(3,2);
const double clhs290 =             clhs289*clhs6;
const double clhs291 =             DN(1,2)*N[3];
const double clhs292 =             DN(3,2)*N[1];
const double clhs293 =             N[1] + clhs223;
const double clhs294 =             clhs25*(clhs244 + clhs266 + clhs283);
const double clhs295 =             clhs25*(clhs253 + clhs273 + clhs289);
const double clhs296 =             DN(2,0)*clhs8;
const double clhs297 =             DN(2,1)*clhs8;
const double clhs298 =             DN(2,2)*clhs8;
const double clhs299 =             N[2]*clhs8;
const double clhs300 =             clhs27*clhs89;
const double clhs301 =             clhs299*clhs32;
const double clhs302 =             clhs17*clhs299 + clhs23*clhs300 + clhs23*clhs301;
const double clhs303 =             clhs49*clhs89;
const double clhs304 =             clhs299*rho;
const double clhs305 =             clhs300*clhs64 + clhs301*clhs64 + clhs304*clhs62;
const double clhs306 =             pow(DN(2,0), 2);
const double clhs307 =             pow(N[2], 2)*clhs30 + clhs300*clhs90 + clhs301*clhs90 + clhs304*clhs89;
const double clhs308 =             DN(2,0)*clhs6;
const double clhs309 =             DN(2,1)*clhs308;
const double clhs310 =             DN(2,2)*clhs308;
const double clhs311 =             N[2]*clhs50 - N[2] + clhs303;
const double clhs312 =             DN(2,0)*DN(3,0);
const double clhs313 =             N[3]*clhs29*clhs299;
const double clhs314 =             clhs312*clhs6 + clhs313;
const double clhs315 =             clhs115*clhs304 + clhs116*clhs300 + clhs116*clhs301;
const double clhs316 =             DN(3,1)*clhs308;
const double clhs317 =             DN(3,2)*clhs308;
const double clhs318 =             DN(2,0)*N[3];
const double clhs319 =             DN(3,0)*N[2];
const double clhs320 =             clhs302 + clhs87;
const double clhs321 =             clhs246 + clhs305;
const double clhs322 =             pow(DN(2,1), 2);
const double clhs323 =             DN(2,1)*clhs6;
const double clhs324 =             DN(2,2)*clhs323;
const double clhs325 =             DN(3,0)*clhs323;
const double clhs326 =             DN(2,1)*DN(3,1);
const double clhs327 =             clhs326*clhs6;
const double clhs328 =             clhs313 + clhs315;
const double clhs329 =             DN(3,2)*clhs323;
const double clhs330 =             DN(2,1)*N[3];
const double clhs331 =             DN(3,1)*N[2];
const double clhs332 =             pow(DN(2,2), 2);
const double clhs333 =             DN(2,2)*clhs6;
const double clhs334 =             DN(3,0)*clhs333;
const double clhs335 =             DN(3,1)*clhs333;
const double clhs336 =             DN(2,2)*DN(3,2);
const double clhs337 =             clhs336*clhs6;
const double clhs338 =             DN(2,2)*N[3];
const double clhs339 =             DN(3,2)*N[2];
const double clhs340 =             N[2] + clhs225;
const double clhs341 =             clhs25*(clhs312 + clhs326 + clhs336);
const double clhs342 =             DN(3,0)*clhs8;
const double clhs343 =             DN(3,1)*clhs8;
const double clhs344 =             DN(3,2)*clhs8;
const double clhs345 =             N[3]*clhs8;
const double clhs346 =             clhs115*clhs27;
const double clhs347 =             clhs32*clhs345;
const double clhs348 =             clhs17*clhs345 + clhs23*clhs346 + clhs23*clhs347;
const double clhs349 =             clhs115*clhs49;
const double clhs350 =             clhs345*rho;
const double clhs351 =             clhs346*clhs64 + clhs347*clhs64 + clhs350*clhs62;
const double clhs352 =             clhs346*clhs90 + clhs347*clhs90 + clhs350*clhs89;
const double clhs353 =             pow(DN(3,0), 2);
const double clhs354 =             pow(N[3], 2)*clhs30 + clhs115*clhs350 + clhs116*clhs346 + clhs116*clhs347;
const double clhs355 =             DN(3,0)*clhs6;
const double clhs356 =             DN(3,1)*clhs355;
const double clhs357 =             DN(3,2)*clhs355;
const double clhs358 =             N[3]*clhs50 - N[3] + clhs349;
const double clhs359 =             clhs113 + clhs348;
const double clhs360 =             clhs254 + clhs351;
const double clhs361 =             clhs313 + clhs352;
const double clhs362 =             pow(DN(3,1), 2);
const double clhs363 =             DN(3,1)*DN(3,2)*clhs6;
const double clhs364 =             pow(DN(3,2), 2);
const double clhs365 =             N[3] + clhs227;
            lhs(0,0)=clhs0*clhs6 + clhs11*clhs12 + clhs14*clhs15 + clhs34 + clhs7*clhs9;
            lhs(0,1)=clhs12*clhs37 + clhs15*clhs40 + clhs35*clhs9 + clhs42;
            lhs(0,2)=clhs12*clhs45 + clhs15*clhs47 + clhs43*clhs9 + clhs48;
            lhs(0,3)=DN(0,0)*clhs52;
            lhs(0,4)=clhs12*clhs55 + clhs15*clhs57 + clhs53*clhs9 + clhs61 + clhs65;
            lhs(0,5)=clhs12*clhs68 + clhs15*clhs71 + clhs66*clhs9 + clhs72;
            lhs(0,6)=clhs12*clhs75 + clhs15*clhs77 + clhs73*clhs9 + clhs78;
            lhs(0,7)=DN(1,0)*clhs51 + clhs50*clhs80 - clhs79;
            lhs(0,8)=clhs12*clhs83 + clhs15*clhs85 + clhs81*clhs9 + clhs88 + clhs91;
            lhs(0,9)=clhs12*clhs94 + clhs15*clhs97 + clhs9*clhs92 + clhs98;
            lhs(0,10)=clhs101*clhs12 + clhs103*clhs15 + clhs104 + clhs9*clhs99;
            lhs(0,11)=DN(2,0)*clhs51 - clhs105 + clhs106*clhs50;
            lhs(0,12)=clhs107*clhs9 + clhs109*clhs12 + clhs111*clhs15 + clhs114 + clhs117;
            lhs(0,13)=clhs118*clhs9 + clhs12*clhs120 + clhs123*clhs15 + clhs124;
            lhs(0,14)=clhs12*clhs127 + clhs125*clhs9 + clhs129*clhs15 + clhs130;
            lhs(0,15)=DN(3,0)*clhs51 - clhs131 + clhs132*clhs50;
            lhs(1,0)=clhs11*clhs9 + clhs12*clhs133 + clhs134*clhs15 + clhs42;
            lhs(1,1)=clhs12*clhs136 + clhs135*clhs6 + clhs138*clhs15 + clhs34 + clhs37*clhs9;
            lhs(1,2)=clhs12*clhs139 + clhs141*clhs15 + clhs143 + clhs45*clhs9;
            lhs(1,3)=DN(0,1)*clhs52;
            lhs(1,4)=clhs12*clhs144 + clhs145*clhs15 + clhs146 + clhs55*clhs9;
            lhs(1,5)=clhs12*clhs149 + clhs148 + clhs15*clhs151 + clhs152 + clhs68*clhs9;
            lhs(1,6)=clhs12*clhs153 + clhs15*clhs155 + clhs156 + clhs75*clhs9;
            lhs(1,7)=DN(1,1)*clhs51 - clhs157 + clhs158*clhs50;
            lhs(1,8)=clhs12*clhs159 + clhs15*clhs160 + clhs161 + clhs83*clhs9;
            lhs(1,9)=clhs12*clhs164 + clhs15*clhs166 + clhs163 + clhs167 + clhs9*clhs94;
            lhs(1,10)=clhs101*clhs9 + clhs12*clhs168 + clhs15*clhs170 + clhs171;
            lhs(1,11)=DN(2,1)*clhs51 - clhs172 + clhs173*clhs50;
            lhs(1,12)=clhs109*clhs9 + clhs12*clhs174 + clhs15*clhs175 + clhs176;
            lhs(1,13)=clhs12*clhs179 + clhs120*clhs9 + clhs15*clhs181 + clhs178 + clhs182;
            lhs(1,14)=clhs12*clhs183 + clhs127*clhs9 + clhs15*clhs185 + clhs186;
            lhs(1,15)=DN(3,1)*clhs51 - clhs187 + clhs188*clhs50;
            lhs(2,0)=clhs12*clhs134 + clhs14*clhs9 + clhs15*clhs189 + clhs48;
            lhs(2,1)=clhs12*clhs138 + clhs143 + clhs15*clhs190 + clhs40*clhs9;
            lhs(2,2)=clhs12*clhs141 + clhs15*clhs192 + clhs191*clhs6 + clhs34 + clhs47*clhs9;
            lhs(2,3)=DN(0,2)*clhs52;
            lhs(2,4)=clhs12*clhs145 + clhs15*clhs193 + clhs195 + clhs57*clhs9;
            lhs(2,5)=clhs12*clhs151 + clhs15*clhs196 + clhs197 + clhs71*clhs9;
            lhs(2,6)=clhs12*clhs155 + clhs15*clhs200 + clhs152 + clhs199 + clhs77*clhs9;
            lhs(2,7)=DN(1,2)*clhs51 - clhs201 + clhs202*clhs50;
            lhs(2,8)=clhs12*clhs160 + clhs15*clhs203 + clhs204 + clhs85*clhs9;
            lhs(2,9)=clhs12*clhs166 + clhs15*clhs205 + clhs206 + clhs9*clhs97;
            lhs(2,10)=clhs103*clhs9 + clhs12*clhs170 + clhs15*clhs209 + clhs167 + clhs208;
            lhs(2,11)=DN(2,2)*clhs51 - clhs210 + clhs211*clhs50;
            lhs(2,12)=clhs111*clhs9 + clhs12*clhs175 + clhs15*clhs212 + clhs213;
            lhs(2,13)=clhs12*clhs181 + clhs123*clhs9 + clhs15*clhs214 + clhs215;
            lhs(2,14)=clhs12*clhs185 + clhs129*clhs9 + clhs15*clhs218 + clhs182 + clhs217;
            lhs(2,15)=DN(3,2)*clhs51 - clhs219 + clhs220*clhs50;
            lhs(3,0)=clhs222*clhs9;
            lhs(3,1)=clhs12*clhs222;
            lhs(3,2)=clhs15*clhs222;
            lhs(3,3)=clhs25*(clhs0 + clhs135 + clhs191);
            lhs(3,4)=clhs8*(DN(0,0)*clhs223 + clhs80);
            lhs(3,5)=clhs8*(DN(0,1)*clhs223 + clhs158);
            lhs(3,6)=clhs8*(DN(0,2)*clhs223 + clhs202);
            lhs(3,7)=clhs224;
            lhs(3,8)=clhs8*(DN(0,0)*clhs225 + clhs106);
            lhs(3,9)=clhs8*(DN(0,1)*clhs225 + clhs173);
            lhs(3,10)=clhs8*(DN(0,2)*clhs225 + clhs211);
            lhs(3,11)=clhs226;
            lhs(3,12)=clhs8*(DN(0,0)*clhs227 + clhs132);
            lhs(3,13)=clhs8*(DN(0,1)*clhs227 + clhs188);
            lhs(3,14)=clhs8*(DN(0,2)*clhs227 + clhs220);
            lhs(3,15)=clhs228;
            lhs(4,0)=clhs11*clhs230 + clhs14*clhs231 + clhs229*clhs7 + clhs235 + clhs61;
            lhs(4,1)=clhs146 + clhs229*clhs35 + clhs230*clhs37 + clhs231*clhs40;
            lhs(4,2)=clhs195 + clhs229*clhs43 + clhs230*clhs45 + clhs231*clhs47;
            lhs(4,3)=DN(0,0)*clhs236 + clhs50*clhs79 - clhs80;
            lhs(4,4)=clhs229*clhs53 + clhs230*clhs55 + clhs231*clhs57 + clhs237*clhs6 + clhs239;
            lhs(4,5)=clhs229*clhs66 + clhs230*clhs68 + clhs231*clhs71 + clhs241;
            lhs(4,6)=clhs229*clhs73 + clhs230*clhs75 + clhs231*clhs77 + clhs242;
            lhs(4,7)=DN(1,0)*clhs243;
            lhs(4,8)=clhs229*clhs81 + clhs230*clhs83 + clhs231*clhs85 + clhs247 + clhs248;
            lhs(4,9)=clhs229*clhs92 + clhs230*clhs94 + clhs231*clhs97 + clhs249;
            lhs(4,10)=clhs101*clhs230 + clhs103*clhs231 + clhs229*clhs99 + clhs250;
            lhs(4,11)=DN(2,0)*clhs236 - clhs251 + clhs252*clhs50;
            lhs(4,12)=clhs107*clhs229 + clhs109*clhs230 + clhs111*clhs231 + clhs255 + clhs256;
            lhs(4,13)=clhs118*clhs229 + clhs120*clhs230 + clhs123*clhs231 + clhs257;
            lhs(4,14)=clhs125*clhs229 + clhs127*clhs230 + clhs129*clhs231 + clhs258;
            lhs(4,15)=DN(3,0)*clhs236 - clhs259 + clhs260*clhs50;
            lhs(5,0)=clhs11*clhs229 + clhs133*clhs230 + clhs134*clhs231 + clhs72;
            lhs(5,1)=clhs136*clhs230 + clhs138*clhs231 + clhs148 + clhs229*clhs37 + clhs261;
            lhs(5,2)=clhs139*clhs230 + clhs141*clhs231 + clhs197 + clhs229*clhs45;
            lhs(5,3)=DN(0,1)*clhs236 + clhs157*clhs50 - clhs158;
            lhs(5,4)=clhs144*clhs230 + clhs145*clhs231 + clhs229*clhs55 + clhs241;
            lhs(5,5)=clhs149*clhs230 + clhs151*clhs231 + clhs229*clhs68 + clhs239 + clhs262*clhs6;
            lhs(5,6)=clhs153*clhs230 + clhs155*clhs231 + clhs229*clhs75 + clhs264;
            lhs(5,7)=DN(1,1)*clhs243;
            lhs(5,8)=clhs159*clhs230 + clhs160*clhs231 + clhs229*clhs83 + clhs265;
            lhs(5,9)=clhs164*clhs230 + clhs166*clhs231 + clhs229*clhs94 + clhs267 + clhs268;
            lhs(5,10)=clhs101*clhs229 + clhs168*clhs230 + clhs170*clhs231 + clhs269;
            lhs(5,11)=DN(2,1)*clhs236 - clhs270 + clhs271*clhs50;
            lhs(5,12)=clhs109*clhs229 + clhs174*clhs230 + clhs175*clhs231 + clhs272;
            lhs(5,13)=clhs120*clhs229 + clhs179*clhs230 + clhs181*clhs231 + clhs274 + clhs275;
            lhs(5,14)=clhs127*clhs229 + clhs183*clhs230 + clhs185*clhs231 + clhs276;
            lhs(5,15)=DN(3,1)*clhs236 - clhs277 + clhs278*clhs50;
            lhs(6,0)=clhs134*clhs230 + clhs14*clhs229 + clhs189*clhs231 + clhs78;
            lhs(6,1)=clhs138*clhs230 + clhs156 + clhs190*clhs231 + clhs229*clhs40;
            lhs(6,2)=clhs141*clhs230 + clhs192*clhs231 + clhs199 + clhs229*clhs47 + clhs261;
            lhs(6,3)=DN(0,2)*clhs236 + clhs201*clhs50 - clhs202;
            lhs(6,4)=clhs145*clhs230 + clhs193*clhs231 + clhs229*clhs57 + clhs242;
            lhs(6,5)=clhs151*clhs230 + clhs196*clhs231 + clhs229*clhs71 + clhs264;
            lhs(6,6)=clhs155*clhs230 + clhs200*clhs231 + clhs229*clhs77 + clhs239 + clhs279*clhs6;
            lhs(6,7)=DN(1,2)*clhs243;
            lhs(6,8)=clhs160*clhs230 + clhs203*clhs231 + clhs229*clhs85 + clhs281;
            lhs(6,9)=clhs166*clhs230 + clhs205*clhs231 + clhs229*clhs97 + clhs282;
            lhs(6,10)=clhs103*clhs229 + clhs170*clhs230 + clhs209*clhs231 + clhs268 + clhs284;
            lhs(6,11)=DN(2,2)*clhs236 - clhs285 + clhs286*clhs50;
            lhs(6,12)=clhs111*clhs229 + clhs175*clhs230 + clhs212*clhs231 + clhs287;
            lhs(6,13)=clhs123*clhs229 + clhs181*clhs230 + clhs214*clhs231 + clhs288;
            lhs(6,14)=clhs129*clhs229 + clhs185*clhs230 + clhs218*clhs231 + clhs275 + clhs290;
            lhs(6,15)=DN(3,2)*clhs236 - clhs291 + clhs292*clhs50;
            lhs(7,0)=clhs8*(DN(1,0)*clhs221 + clhs79);
            lhs(7,1)=clhs8*(DN(1,1)*clhs221 + clhs157);
            lhs(7,2)=clhs8*(DN(1,2)*clhs221 + clhs201);
            lhs(7,3)=clhs224;
            lhs(7,4)=clhs229*clhs293;
            lhs(7,5)=clhs230*clhs293;
            lhs(7,6)=clhs231*clhs293;
            lhs(7,7)=clhs25*(clhs237 + clhs262 + clhs279);
            lhs(7,8)=clhs8*(DN(1,0)*clhs225 + clhs252);
            lhs(7,9)=clhs8*(DN(1,1)*clhs225 + clhs271);
            lhs(7,10)=clhs8*(DN(1,2)*clhs225 + clhs286);
            lhs(7,11)=clhs294;
            lhs(7,12)=clhs8*(DN(1,0)*clhs227 + clhs260);
            lhs(7,13)=clhs8*(DN(1,1)*clhs227 + clhs278);
            lhs(7,14)=clhs8*(DN(1,2)*clhs227 + clhs292);
            lhs(7,15)=clhs295;
            lhs(8,0)=clhs11*clhs297 + clhs14*clhs298 + clhs296*clhs7 + clhs302 + clhs88;
            lhs(8,1)=clhs161 + clhs296*clhs35 + clhs297*clhs37 + clhs298*clhs40;
            lhs(8,2)=clhs204 + clhs296*clhs43 + clhs297*clhs45 + clhs298*clhs47;
            lhs(8,3)=DN(0,0)*clhs303 + clhs105*clhs50 - clhs106;
            lhs(8,4)=clhs247 + clhs296*clhs53 + clhs297*clhs55 + clhs298*clhs57 + clhs305;
            lhs(8,5)=clhs265 + clhs296*clhs66 + clhs297*clhs68 + clhs298*clhs71;
            lhs(8,6)=clhs281 + clhs296*clhs73 + clhs297*clhs75 + clhs298*clhs77;
            lhs(8,7)=DN(1,0)*clhs303 + clhs251*clhs50 - clhs252;
            lhs(8,8)=clhs296*clhs81 + clhs297*clhs83 + clhs298*clhs85 + clhs306*clhs6 + clhs307;
            lhs(8,9)=clhs296*clhs92 + clhs297*clhs94 + clhs298*clhs97 + clhs309;
            lhs(8,10)=clhs101*clhs297 + clhs103*clhs298 + clhs296*clhs99 + clhs310;
            lhs(8,11)=DN(2,0)*clhs311;
            lhs(8,12)=clhs107*clhs296 + clhs109*clhs297 + clhs111*clhs298 + clhs314 + clhs315;
            lhs(8,13)=clhs118*clhs296 + clhs120*clhs297 + clhs123*clhs298 + clhs316;
            lhs(8,14)=clhs125*clhs296 + clhs127*clhs297 + clhs129*clhs298 + clhs317;
            lhs(8,15)=DN(3,0)*clhs303 - clhs318 + clhs319*clhs50;
            lhs(9,0)=clhs11*clhs296 + clhs133*clhs297 + clhs134*clhs298 + clhs98;
            lhs(9,1)=clhs136*clhs297 + clhs138*clhs298 + clhs163 + clhs296*clhs37 + clhs320;
            lhs(9,2)=clhs139*clhs297 + clhs141*clhs298 + clhs206 + clhs296*clhs45;
            lhs(9,3)=DN(0,1)*clhs303 + clhs172*clhs50 - clhs173;
            lhs(9,4)=clhs144*clhs297 + clhs145*clhs298 + clhs249 + clhs296*clhs55;
            lhs(9,5)=clhs149*clhs297 + clhs151*clhs298 + clhs267 + clhs296*clhs68 + clhs321;
            lhs(9,6)=clhs153*clhs297 + clhs155*clhs298 + clhs282 + clhs296*clhs75;
            lhs(9,7)=DN(1,1)*clhs303 + clhs270*clhs50 - clhs271;
            lhs(9,8)=clhs159*clhs297 + clhs160*clhs298 + clhs296*clhs83 + clhs309;
            lhs(9,9)=clhs164*clhs297 + clhs166*clhs298 + clhs296*clhs94 + clhs307 + clhs322*clhs6;
            lhs(9,10)=clhs101*clhs296 + clhs168*clhs297 + clhs170*clhs298 + clhs324;
            lhs(9,11)=DN(2,1)*clhs311;
            lhs(9,12)=clhs109*clhs296 + clhs174*clhs297 + clhs175*clhs298 + clhs325;
            lhs(9,13)=clhs120*clhs296 + clhs179*clhs297 + clhs181*clhs298 + clhs327 + clhs328;
            lhs(9,14)=clhs127*clhs296 + clhs183*clhs297 + clhs185*clhs298 + clhs329;
            lhs(9,15)=DN(3,1)*clhs303 - clhs330 + clhs331*clhs50;
            lhs(10,0)=clhs104 + clhs134*clhs297 + clhs14*clhs296 + clhs189*clhs298;
            lhs(10,1)=clhs138*clhs297 + clhs171 + clhs190*clhs298 + clhs296*clhs40;
            lhs(10,2)=clhs141*clhs297 + clhs192*clhs298 + clhs208 + clhs296*clhs47 + clhs320;
            lhs(10,3)=DN(0,2)*clhs303 + clhs210*clhs50 - clhs211;
            lhs(10,4)=clhs145*clhs297 + clhs193*clhs298 + clhs250 + clhs296*clhs57;
            lhs(10,5)=clhs151*clhs297 + clhs196*clhs298 + clhs269 + clhs296*clhs71;
            lhs(10,6)=clhs155*clhs297 + clhs200*clhs298 + clhs284 + clhs296*clhs77 + clhs321;
            lhs(10,7)=DN(1,2)*clhs303 + clhs285*clhs50 - clhs286;
            lhs(10,8)=clhs160*clhs297 + clhs203*clhs298 + clhs296*clhs85 + clhs310;
            lhs(10,9)=clhs166*clhs297 + clhs205*clhs298 + clhs296*clhs97 + clhs324;
            lhs(10,10)=clhs103*clhs296 + clhs170*clhs297 + clhs209*clhs298 + clhs307 + clhs332*clhs6;
            lhs(10,11)=DN(2,2)*clhs311;
            lhs(10,12)=clhs111*clhs296 + clhs175*clhs297 + clhs212*clhs298 + clhs334;
            lhs(10,13)=clhs123*clhs296 + clhs181*clhs297 + clhs214*clhs298 + clhs335;
            lhs(10,14)=clhs129*clhs296 + clhs185*clhs297 + clhs218*clhs298 + clhs328 + clhs337;
            lhs(10,15)=DN(3,2)*clhs303 - clhs338 + clhs339*clhs50;
            lhs(11,0)=clhs8*(DN(2,0)*clhs221 + clhs105);
            lhs(11,1)=clhs8*(DN(2,1)*clhs221 + clhs172);
            lhs(11,2)=clhs8*(DN(2,2)*clhs221 + clhs210);
            lhs(11,3)=clhs226;
            lhs(11,4)=clhs8*(DN(2,0)*clhs223 + clhs251);
            lhs(11,5)=clhs8*(DN(2,1)*clhs223 + clhs270);
            lhs(11,6)=clhs8*(DN(2,2)*clhs223 + clhs285);
            lhs(11,7)=clhs294;
            lhs(11,8)=clhs296*clhs340;
            lhs(11,9)=clhs297*clhs340;
            lhs(11,10)=clhs298*clhs340;
            lhs(11,11)=clhs25*(clhs306 + clhs322 + clhs332);
            lhs(11,12)=clhs8*(DN(2,0)*clhs227 + clhs319);
            lhs(11,13)=clhs8*(DN(2,1)*clhs227 + clhs331);
            lhs(11,14)=clhs8*(DN(2,2)*clhs227 + clhs339);
            lhs(11,15)=clhs341;
            lhs(12,0)=clhs11*clhs343 + clhs114 + clhs14*clhs344 + clhs342*clhs7 + clhs348;
            lhs(12,1)=clhs176 + clhs342*clhs35 + clhs343*clhs37 + clhs344*clhs40;
            lhs(12,2)=clhs213 + clhs342*clhs43 + clhs343*clhs45 + clhs344*clhs47;
            lhs(12,3)=DN(0,0)*clhs349 + clhs131*clhs50 - clhs132;
            lhs(12,4)=clhs255 + clhs342*clhs53 + clhs343*clhs55 + clhs344*clhs57 + clhs351;
            lhs(12,5)=clhs272 + clhs342*clhs66 + clhs343*clhs68 + clhs344*clhs71;
            lhs(12,6)=clhs287 + clhs342*clhs73 + clhs343*clhs75 + clhs344*clhs77;
            lhs(12,7)=DN(1,0)*clhs349 + clhs259*clhs50 - clhs260;
            lhs(12,8)=clhs314 + clhs342*clhs81 + clhs343*clhs83 + clhs344*clhs85 + clhs352;
            lhs(12,9)=clhs325 + clhs342*clhs92 + clhs343*clhs94 + clhs344*clhs97;
            lhs(12,10)=clhs101*clhs343 + clhs103*clhs344 + clhs334 + clhs342*clhs99;
            lhs(12,11)=DN(2,0)*clhs349 + clhs318*clhs50 - clhs319;
            lhs(12,12)=clhs107*clhs342 + clhs109*clhs343 + clhs111*clhs344 + clhs353*clhs6 + clhs354;
            lhs(12,13)=clhs118*clhs342 + clhs120*clhs343 + clhs123*clhs344 + clhs356;
            lhs(12,14)=clhs125*clhs342 + clhs127*clhs343 + clhs129*clhs344 + clhs357;
            lhs(12,15)=DN(3,0)*clhs358;
            lhs(13,0)=clhs11*clhs342 + clhs124 + clhs133*clhs343 + clhs134*clhs344;
            lhs(13,1)=clhs136*clhs343 + clhs138*clhs344 + clhs178 + clhs342*clhs37 + clhs359;
            lhs(13,2)=clhs139*clhs343 + clhs141*clhs344 + clhs215 + clhs342*clhs45;
            lhs(13,3)=DN(0,1)*clhs349 + clhs187*clhs50 - clhs188;
            lhs(13,4)=clhs144*clhs343 + clhs145*clhs344 + clhs257 + clhs342*clhs55;
            lhs(13,5)=clhs149*clhs343 + clhs151*clhs344 + clhs274 + clhs342*clhs68 + clhs360;
            lhs(13,6)=clhs153*clhs343 + clhs155*clhs344 + clhs288 + clhs342*clhs75;
            lhs(13,7)=DN(1,1)*clhs349 + clhs277*clhs50 - clhs278;
            lhs(13,8)=clhs159*clhs343 + clhs160*clhs344 + clhs316 + clhs342*clhs83;
            lhs(13,9)=clhs164*clhs343 + clhs166*clhs344 + clhs327 + clhs342*clhs94 + clhs361;
            lhs(13,10)=clhs101*clhs342 + clhs168*clhs343 + clhs170*clhs344 + clhs335;
            lhs(13,11)=DN(2,1)*clhs349 + clhs330*clhs50 - clhs331;
            lhs(13,12)=clhs109*clhs342 + clhs174*clhs343 + clhs175*clhs344 + clhs356;
            lhs(13,13)=clhs120*clhs342 + clhs179*clhs343 + clhs181*clhs344 + clhs354 + clhs362*clhs6;
            lhs(13,14)=clhs127*clhs342 + clhs183*clhs343 + clhs185*clhs344 + clhs363;
            lhs(13,15)=DN(3,1)*clhs358;
            lhs(14,0)=clhs130 + clhs134*clhs343 + clhs14*clhs342 + clhs189*clhs344;
            lhs(14,1)=clhs138*clhs343 + clhs186 + clhs190*clhs344 + clhs342*clhs40;
            lhs(14,2)=clhs141*clhs343 + clhs192*clhs344 + clhs217 + clhs342*clhs47 + clhs359;
            lhs(14,3)=DN(0,2)*clhs349 + clhs219*clhs50 - clhs220;
            lhs(14,4)=clhs145*clhs343 + clhs193*clhs344 + clhs258 + clhs342*clhs57;
            lhs(14,5)=clhs151*clhs343 + clhs196*clhs344 + clhs276 + clhs342*clhs71;
            lhs(14,6)=clhs155*clhs343 + clhs200*clhs344 + clhs290 + clhs342*clhs77 + clhs360;
            lhs(14,7)=DN(1,2)*clhs349 + clhs291*clhs50 - clhs292;
            lhs(14,8)=clhs160*clhs343 + clhs203*clhs344 + clhs317 + clhs342*clhs85;
            lhs(14,9)=clhs166*clhs343 + clhs205*clhs344 + clhs329 + clhs342*clhs97;
            lhs(14,10)=clhs103*clhs342 + clhs170*clhs343 + clhs209*clhs344 + clhs337 + clhs361;
            lhs(14,11)=DN(2,2)*clhs349 + clhs338*clhs50 - clhs339;
            lhs(14,12)=clhs111*clhs342 + clhs175*clhs343 + clhs212*clhs344 + clhs357;
            lhs(14,13)=clhs123*clhs342 + clhs181*clhs343 + clhs214*clhs344 + clhs363;
            lhs(14,14)=clhs129*clhs342 + clhs185*clhs343 + clhs218*clhs344 + clhs354 + clhs364*clhs6;
            lhs(14,15)=DN(3,2)*clhs358;
            lhs(15,0)=clhs8*(DN(3,0)*clhs221 + clhs131);
            lhs(15,1)=clhs8*(DN(3,1)*clhs221 + clhs187);
            lhs(15,2)=clhs8*(DN(3,2)*clhs221 + clhs219);
            lhs(15,3)=clhs228;
            lhs(15,4)=clhs8*(DN(3,0)*clhs223 + clhs259);
            lhs(15,5)=clhs8*(DN(3,1)*clhs223 + clhs277);
            lhs(15,6)=clhs8*(DN(3,2)*clhs223 + clhs291);
            lhs(15,7)=clhs295;
            lhs(15,8)=clhs8*(DN(3,0)*clhs225 + clhs318);
            lhs(15,9)=clhs8*(DN(3,1)*clhs225 + clhs330);
            lhs(15,10)=clhs8*(DN(3,2)*clhs225 + clhs338);
            lhs(15,11)=clhs341;
            lhs(15,12)=clhs342*clhs365;
            lhs(15,13)=clhs343*clhs365;
            lhs(15,14)=clhs344*clhs365;
            lhs(15,15)=clhs25*(clhs353 + clhs362 + clhs364);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,3,2> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             1.0/(max_spectral_radius + 1);
const double crhs2 =             rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)));
const double crhs3 =             mu*stab_c1;
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 =             rho*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs7 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs8 =             (crhs7 - volume_error_time_ratio)*(crhs3/stab_c2 + crhs6*h);
const double crhs9 =             v(0,0) - vn(0,0);
const double crhs10 =             crhs1*crhs9 + vn(0,0);
const double crhs11 =             v(1,0) - vn(1,0);
const double crhs12 =             crhs1*crhs11 + vn(1,0);
const double crhs13 =             v(2,0) - vn(2,0);
const double crhs14 =             crhs1*crhs13 + vn(2,0);
const double crhs15 =             rho*(crhs4*(DN(0,0)*crhs10 + DN(1,0)*crhs12 + DN(2,0)*crhs14) + crhs5*(DN(0,1)*crhs10 + DN(1,1)*crhs12 + DN(2,1)*crhs14));
const double crhs16 =             -acceleration_alpha_method(0,0);
const double crhs17 =             1.0/dt;
const double crhs18 =             0.5*max_spectral_radius;
const double crhs19 =             -crhs1;
const double crhs20 =             crhs19 + 0.5;
const double crhs21 =             1.0/(-crhs1*(crhs18 - 1.5) + crhs20);
const double crhs22 =             crhs17*crhs21;
const double crhs23 =             crhs1*(1.5 - crhs18);
const double crhs24 =             crhs21*(crhs1 - crhs23 + 0.5);
const double crhs25 =             0.5*crhs1;
const double crhs26 =             crhs25*(max_spectral_radius - 3);
const double crhs27 =             -acceleration_alpha_method(1,0);
const double crhs28 =             -acceleration_alpha_method(2,0);
const double crhs29 =             N[0]*(acceleration_alpha_method(0,0) - crhs26*(-acceleration_alpha_method(0,0)*crhs24 + crhs16 + crhs22*crhs9)) + N[1]*(acceleration_alpha_method(1,0) - crhs26*(-acceleration_alpha_method(1,0)*crhs24 + crhs11*crhs22 + crhs27)) + N[2]*(acceleration_alpha_method(2,0) - crhs26*(-acceleration_alpha_method(2,0)*crhs24 + crhs13*crhs22 + crhs28));
const double crhs30 =             N[0]*rho;
const double crhs31 =             1.0/(crhs20 + crhs23);
const double crhs32 =             crhs17*crhs31;
const double crhs33 =             crhs31*(crhs19 + crhs23 - 0.5);
const double crhs34 =             crhs25*(3 - max_spectral_radius);
const double crhs35 =             1.0/(crhs17*dyn_tau*rho + crhs3/pow(h, 2) + crhs6*stab_c2/h);
const double crhs36 =             crhs35*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs15 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs34*(acceleration_alpha_method(0,0)*crhs33 + crhs16 + crhs32*crhs9)) + N[1]*(acceleration_alpha_method(1,0) + crhs34*(acceleration_alpha_method(1,0)*crhs33 + crhs11*crhs32 + crhs27)) + N[2]*(acceleration_alpha_method(2,0) + crhs34*(acceleration_alpha_method(2,0)*crhs33 + crhs13*crhs32 + crhs28))));
const double crhs37 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs38 =             crhs30*crhs37;
const double crhs39 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs40 =             rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)));
const double crhs41 =             v(0,1) - vn(0,1);
const double crhs42 =             crhs1*crhs41 + vn(0,1);
const double crhs43 =             v(1,1) - vn(1,1);
const double crhs44 =             crhs1*crhs43 + vn(1,1);
const double crhs45 =             v(2,1) - vn(2,1);
const double crhs46 =             crhs1*crhs45 + vn(2,1);
const double crhs47 =             rho*(crhs4*(DN(0,0)*crhs42 + DN(1,0)*crhs44 + DN(2,0)*crhs46) + crhs5*(DN(0,1)*crhs42 + DN(1,1)*crhs44 + DN(2,1)*crhs46));
const double crhs48 =             -acceleration_alpha_method(0,1);
const double crhs49 =             -acceleration_alpha_method(1,1);
const double crhs50 =             -acceleration_alpha_method(2,1);
const double crhs51 =             N[0]*(acceleration_alpha_method(0,1) - crhs26*(-acceleration_alpha_method(0,1)*crhs24 + crhs22*crhs41 + crhs48)) + N[1]*(acceleration_alpha_method(1,1) - crhs26*(-acceleration_alpha_method(1,1)*crhs24 + crhs22*crhs43 + crhs49)) + N[2]*(acceleration_alpha_method(2,1) - crhs26*(-acceleration_alpha_method(2,1)*crhs24 + crhs22*crhs45 + crhs50));
const double crhs52 =             crhs35*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs40 + crhs47 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs34*(acceleration_alpha_method(0,1)*crhs33 + crhs32*crhs41 + crhs48)) + N[1]*(acceleration_alpha_method(1,1) + crhs34*(acceleration_alpha_method(1,1)*crhs33 + crhs32*crhs43 + crhs49)) + N[2]*(acceleration_alpha_method(2,1) + crhs34*(acceleration_alpha_method(2,1)*crhs33 + crhs32*crhs45 + crhs50))));
const double crhs53 =             -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + crhs7) + volume_error_time_ratio;
const double crhs54 =             N[1]*rho;
const double crhs55 =             crhs37*crhs54;
const double crhs56 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs57 =             N[2]*rho;
const double crhs58 =             crhs37*crhs57;
const double crhs59 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs15 + N[0]*crhs2 - crhs29*crhs30 - crhs36*crhs38 - crhs36*crhs39;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs8 - DN(0,1)*stress[1] + N[0]*crhs40 - N[0]*crhs47 - crhs30*crhs51 - crhs38*crhs52 - crhs39*crhs52;
            rhs[2]=-DN(0,0)*crhs36 - DN(0,1)*crhs52 + N[0]*crhs53;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs15 + N[1]*crhs2 - crhs29*crhs54 - crhs36*crhs55 - crhs36*crhs56;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs8 - DN(1,1)*stress[1] + N[1]*crhs40 - N[1]*crhs47 - crhs51*crhs54 - crhs52*crhs55 - crhs52*crhs56;
            rhs[5]=-DN(1,0)*crhs36 - DN(1,1)*crhs52 + N[1]*crhs53;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs15 + N[2]*crhs2 - crhs29*crhs57 - crhs36*crhs58 - crhs36*crhs59;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs8 - DN(2,1)*stress[1] + N[2]*crhs40 - N[2]*crhs47 - crhs51*crhs57 - crhs52*crhs58 - crhs52*crhs59;
            rhs[8]=-DN(2,0)*crhs36 - DN(2,1)*crhs52 + N[2]*crhs53;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_f=1/(1+max_spectral_radius);
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             1.0/(max_spectral_radius + 1);
const double crhs2 =             rho*(N[0]*(crhs1*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs1*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs1*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs1*(f(3,0) - fn(3,0)) + fn(3,0)));
const double crhs3 =             mu*stab_c1;
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs6 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs7 =             rho*sqrt(pow(crhs4, 2) + pow(crhs5, 2) + pow(crhs6, 2));
const double crhs8 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs9 =             (crhs8 - volume_error_time_ratio)*(crhs3/stab_c2 + crhs7*h);
const double crhs10 =             v(0,0) - vn(0,0);
const double crhs11 =             crhs1*crhs10 + vn(0,0);
const double crhs12 =             v(1,0) - vn(1,0);
const double crhs13 =             crhs1*crhs12 + vn(1,0);
const double crhs14 =             v(2,0) - vn(2,0);
const double crhs15 =             crhs1*crhs14 + vn(2,0);
const double crhs16 =             v(3,0) - vn(3,0);
const double crhs17 =             crhs1*crhs16 + vn(3,0);
const double crhs18 =             rho*(crhs4*(DN(0,0)*crhs11 + DN(1,0)*crhs13 + DN(2,0)*crhs15 + DN(3,0)*crhs17) + crhs5*(DN(0,1)*crhs11 + DN(1,1)*crhs13 + DN(2,1)*crhs15 + DN(3,1)*crhs17) + crhs6*(DN(0,2)*crhs11 + DN(1,2)*crhs13 + DN(2,2)*crhs15 + DN(3,2)*crhs17));
const double crhs19 =             -acceleration_alpha_method(0,0);
const double crhs20 =             1.0/dt;
const double crhs21 =             0.5*max_spectral_radius;
const double crhs22 =             -crhs1;
const double crhs23 =             crhs22 + 0.5;
const double crhs24 =             1.0/(-crhs1*(crhs21 - 1.5) + crhs23);
const double crhs25 =             crhs20*crhs24;
const double crhs26 =             crhs1*(1.5 - crhs21);
const double crhs27 =             crhs24*(crhs1 - crhs26 + 0.5);
const double crhs28 =             0.5*crhs1;
const double crhs29 =             crhs28*(max_spectral_radius - 3);
const double crhs30 =             -acceleration_alpha_method(1,0);
const double crhs31 =             -acceleration_alpha_method(2,0);
const double crhs32 =             -acceleration_alpha_method(3,0);
const double crhs33 =             N[0]*(acceleration_alpha_method(0,0) - crhs29*(-acceleration_alpha_method(0,0)*crhs27 + crhs10*crhs25 + crhs19)) + N[1]*(acceleration_alpha_method(1,0) - crhs29*(-acceleration_alpha_method(1,0)*crhs27 + crhs12*crhs25 + crhs30)) + N[2]*(acceleration_alpha_method(2,0) - crhs29*(-acceleration_alpha_method(2,0)*crhs27 + crhs14*crhs25 + crhs31)) + N[3]*(acceleration_alpha_method(3,0) - crhs29*(-acceleration_alpha_method(3,0)*crhs27 + crhs16*crhs25 + crhs32));
const double crhs34 =             N[0]*rho;
const double crhs35 =             1.0/(crhs23 + crhs26);
const double crhs36 =             crhs20*crhs35;
const double crhs37 =             crhs35*(crhs22 + crhs26 - 0.5);
const double crhs38 =             crhs28*(3 - max_spectral_radius);
const double crhs39 =             1.0/(crhs20*dyn_tau*rho + crhs3/pow(h, 2) + crhs7*stab_c2/h);
const double crhs40 =             crhs39*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs18 - crhs2 + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs38*(acceleration_alpha_method(0,0)*crhs37 + crhs10*crhs36 + crhs19)) + N[1]*(acceleration_alpha_method(1,0) + crhs38*(acceleration_alpha_method(1,0)*crhs37 + crhs12*crhs36 + crhs30)) + N[2]*(acceleration_alpha_method(2,0) + crhs38*(acceleration_alpha_method(2,0)*crhs37 + crhs14*crhs36 + crhs31)) + N[3]*(acceleration_alpha_method(3,0) + crhs38*(acceleration_alpha_method(3,0)*crhs37 + crhs16*crhs36 + crhs32))));
const double crhs41 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs42 =             crhs34*crhs41;
const double crhs43 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6);
const double crhs44 =             rho*(N[0]*(crhs1*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs1*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs1*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs1*(f(3,1) - fn(3,1)) + fn(3,1)));
const double crhs45 =             v(0,1) - vn(0,1);
const double crhs46 =             crhs1*crhs45 + vn(0,1);
const double crhs47 =             v(1,1) - vn(1,1);
const double crhs48 =             crhs1*crhs47 + vn(1,1);
const double crhs49 =             v(2,1) - vn(2,1);
const double crhs50 =             crhs1*crhs49 + vn(2,1);
const double crhs51 =             v(3,1) - vn(3,1);
const double crhs52 =             crhs1*crhs51 + vn(3,1);
const double crhs53 =             rho*(crhs4*(DN(0,0)*crhs46 + DN(1,0)*crhs48 + DN(2,0)*crhs50 + DN(3,0)*crhs52) + crhs5*(DN(0,1)*crhs46 + DN(1,1)*crhs48 + DN(2,1)*crhs50 + DN(3,1)*crhs52) + crhs6*(DN(0,2)*crhs46 + DN(1,2)*crhs48 + DN(2,2)*crhs50 + DN(3,2)*crhs52));
const double crhs54 =             -acceleration_alpha_method(0,1);
const double crhs55 =             -acceleration_alpha_method(1,1);
const double crhs56 =             -acceleration_alpha_method(2,1);
const double crhs57 =             -acceleration_alpha_method(3,1);
const double crhs58 =             N[0]*(acceleration_alpha_method(0,1) - crhs29*(-acceleration_alpha_method(0,1)*crhs27 + crhs25*crhs45 + crhs54)) + N[1]*(acceleration_alpha_method(1,1) - crhs29*(-acceleration_alpha_method(1,1)*crhs27 + crhs25*crhs47 + crhs55)) + N[2]*(acceleration_alpha_method(2,1) - crhs29*(-acceleration_alpha_method(2,1)*crhs27 + crhs25*crhs49 + crhs56)) + N[3]*(acceleration_alpha_method(3,1) - crhs29*(-acceleration_alpha_method(3,1)*crhs27 + crhs25*crhs51 + crhs57));
const double crhs59 =             crhs39*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs44 + crhs53 + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs38*(acceleration_alpha_method(0,1)*crhs37 + crhs36*crhs45 + crhs54)) + N[1]*(acceleration_alpha_method(1,1) + crhs38*(acceleration_alpha_method(1,1)*crhs37 + crhs36*crhs47 + crhs55)) + N[2]*(acceleration_alpha_method(2,1) + crhs38*(acceleration_alpha_method(2,1)*crhs37 + crhs36*crhs49 + crhs56)) + N[3]*(acceleration_alpha_method(3,1) + crhs38*(acceleration_alpha_method(3,1)*crhs37 + crhs36*crhs51 + crhs57))));
const double crhs60 =             rho*(N[0]*(crhs1*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs1*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs1*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs1*(f(3,2) - fn(3,2)) + fn(3,2)));
const double crhs61 =             v(0,2) - vn(0,2);
const double crhs62 =             crhs1*crhs61 + vn(0,2);
const double crhs63 =             v(1,2) - vn(1,2);
const double crhs64 =             crhs1*crhs63 + vn(1,2);
const double crhs65 =             v(2,2) - vn(2,2);
const double crhs66 =             crhs1*crhs65 + vn(2,2);
const double crhs67 =             v(3,2) - vn(3,2);
const double crhs68 =             crhs1*crhs67 + vn(3,2);
const double crhs69 =             rho*(crhs4*(DN(0,0)*crhs62 + DN(1,0)*crhs64 + DN(2,0)*crhs66 + DN(3,0)*crhs68) + crhs5*(DN(0,1)*crhs62 + DN(1,1)*crhs64 + DN(2,1)*crhs66 + DN(3,1)*crhs68) + crhs6*(DN(0,2)*crhs62 + DN(1,2)*crhs64 + DN(2,2)*crhs66 + DN(3,2)*crhs68));
const double crhs70 =             -acceleration_alpha_method(0,2);
const double crhs71 =             -acceleration_alpha_method(1,2);
const double crhs72 =             -acceleration_alpha_method(2,2);
const double crhs73 =             -acceleration_alpha_method(3,2);
const double crhs74 =             N[0]*(acceleration_alpha_method(0,2) - crhs29*(-acceleration_alpha_method(0,2)*crhs27 + crhs25*crhs61 + crhs70)) + N[1]*(acceleration_alpha_method(1,2) - crhs29*(-acceleration_alpha_method(1,2)*crhs27 + crhs25*crhs63 + crhs71)) + N[2]*(acceleration_alpha_method(2,2) - crhs29*(-acceleration_alpha_method(2,2)*crhs27 + crhs25*crhs65 + crhs72)) + N[3]*(acceleration_alpha_method(3,2) - crhs29*(-acceleration_alpha_method(3,2)*crhs27 + crhs25*crhs67 + crhs73));
const double crhs75 =             crhs39*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs60 + crhs69 + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs38*(acceleration_alpha_method(0,2)*crhs37 + crhs36*crhs61 + crhs70)) + N[1]*(acceleration_alpha_method(1,2) + crhs38*(acceleration_alpha_method(1,2)*crhs37 + crhs36*crhs63 + crhs71)) + N[2]*(acceleration_alpha_method(2,2) + crhs38*(acceleration_alpha_method(2,2)*crhs37 + crhs36*crhs65 + crhs72)) + N[3]*(acceleration_alpha_method(3,2) + crhs38*(acceleration_alpha_method(3,2)*crhs37 + crhs36*crhs67 + crhs73))));
const double crhs76 =             -crhs1*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + crhs8) + volume_error_time_ratio;
const double crhs77 =             N[1]*rho;
const double crhs78 =             crhs41*crhs77;
const double crhs79 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6);
const double crhs80 =             N[2]*rho;
const double crhs81 =             crhs41*crhs80;
const double crhs82 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6);
const double crhs83 =             N[3]*rho;
const double crhs84 =             crhs41*crhs83;
const double crhs85 =             rho*(DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs9 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs18 + N[0]*crhs2 - crhs33*crhs34 - crhs40*crhs42 - crhs40*crhs43;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs9 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs44 - N[0]*crhs53 - crhs34*crhs58 - crhs42*crhs59 - crhs43*crhs59;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs9 - DN(0,2)*stress[2] + N[0]*crhs60 - N[0]*crhs69 - crhs34*crhs74 - crhs42*crhs75 - crhs43*crhs75;
            rhs[3]=-DN(0,0)*crhs40 - DN(0,1)*crhs59 - DN(0,2)*crhs75 + N[0]*crhs76;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs9 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs18 + N[1]*crhs2 - crhs33*crhs77 - crhs40*crhs78 - crhs40*crhs79;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs9 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs44 - N[1]*crhs53 - crhs58*crhs77 - crhs59*crhs78 - crhs59*crhs79;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs9 - DN(1,2)*stress[2] + N[1]*crhs60 - N[1]*crhs69 - crhs74*crhs77 - crhs75*crhs78 - crhs75*crhs79;
            rhs[7]=-DN(1,0)*crhs40 - DN(1,1)*crhs59 - DN(1,2)*crhs75 + N[1]*crhs76;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs9 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs18 + N[2]*crhs2 - crhs33*crhs80 - crhs40*crhs81 - crhs40*crhs82;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs9 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs44 - N[2]*crhs53 - crhs58*crhs80 - crhs59*crhs81 - crhs59*crhs82;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs9 - DN(2,2)*stress[2] + N[2]*crhs60 - N[2]*crhs69 - crhs74*crhs80 - crhs75*crhs81 - crhs75*crhs82;
            rhs[11]=-DN(2,0)*crhs40 - DN(2,1)*crhs59 - DN(2,2)*crhs75 + N[2]*crhs76;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs9 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs18 + N[3]*crhs2 - crhs33*crhs83 - crhs40*crhs84 - crhs40*crhs85;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs9 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs44 - N[3]*crhs53 - crhs58*crhs83 - crhs59*crhs84 - crhs59*crhs85;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs9 - DN(3,2)*stress[2] + N[3]*crhs60 - N[3]*crhs69 - crhs74*crhs83 - crhs75*crhs84 - crhs75*crhs85;
            rhs[15]=-DN(3,0)*crhs40 - DN(3,1)*crhs59 - DN(3,2)*crhs75 + N[3]*crhs76;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double alpha_f=1/(1+max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p=rData.Pressure;
    const auto &pn=rData.Pressure_OldStep1;

    const BoundedMatrix<double,3,2> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 =             cV2*rho;
const double cV4 =             cV3*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV5 =             cV3*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV6 =             cV3*(DN(2,0)*cV0 + DN(2,1)*cV1);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV4;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV4;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV4;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV4;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV4;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV4;
            V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV5;
            V(3,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV5;
            V(3,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV5;
            V(4,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV5;
            V(4,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV5;
            V(4,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV5;
            V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV6;
            V(6,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV6;
            V(6,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV6;
            V(7,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV6;
            V(7,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV6;
            V(7,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV6;
            V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             1.0/(max_spectral_radius + 1);
const double cH1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH3 =             1.0/dt;
const double cH4 =             0.5*cH3*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH5 =             1.0/(cH3*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2));
const double cH6 =             cH5*rho;
const double cH7 =             cH6*(DN(0,0)*cH1 + DN(0,1)*cH2 + N[0]*cH4);
const double cH8 =             cH6*(DN(1,0)*cH1 + DN(1,1)*cH2 + N[1]*cH4);
const double cH9 =             cH6*(DN(2,0)*cH1 + DN(2,1)*cH2 + N[2]*cH4);
            H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH7);
            H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH7);
            H(0,2)=cH5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH8);
            H(0,4)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH8);
            H(0,5)=cH5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH9);
            H(0,7)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH9);
            H(0,8)=cH5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH7);
            H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH7);
            H(1,2)=cH5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH8);
            H(1,4)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH8);
            H(1,5)=cH5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH9);
            H(1,7)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH9);
            H(1,8)=cH5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH7);
            H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH7);
            H(2,2)=cH5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH8);
            H(2,4)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH8);
            H(2,5)=cH5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH9);
            H(2,7)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH9);
            H(2,8)=cH5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


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


    const double crhs_ee0 =             1.0/(max_spectral_radius + 1);
const double crhs_ee1 =             -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1)) + volume_error_time_ratio;
const double crhs_ee2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee3 =             v(0,0) - vn(0,0);
const double crhs_ee4 =             crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 =             v(1,0) - vn(1,0);
const double crhs_ee6 =             crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 =             v(2,0) - vn(2,0);
const double crhs_ee8 =             crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee10 =             1.0/dt;
const double crhs_ee11 =             crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee12 =             1.0/(crhs_ee11 + 0.5);
const double crhs_ee13 =             crhs_ee10*crhs_ee12;
const double crhs_ee14 =             crhs_ee12*(crhs_ee11 - 0.5);
const double crhs_ee15 =             0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee16 =             1.0/(crhs_ee10*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee2, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee17 =             crhs_ee16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee15*(acceleration_alpha_method(0,0)*crhs_ee14 - acceleration_alpha_method(0,0) + crhs_ee13*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee15*(acceleration_alpha_method(1,0)*crhs_ee14 - acceleration_alpha_method(1,0) + crhs_ee13*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee15*(acceleration_alpha_method(2,0)*crhs_ee14 - acceleration_alpha_method(2,0) + crhs_ee13*crhs_ee7)) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8) + crhs_ee9*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8)));
const double crhs_ee18 =             v(0,1) - vn(0,1);
const double crhs_ee19 =             crhs_ee0*crhs_ee18 + vn(0,1);
const double crhs_ee20 =             v(1,1) - vn(1,1);
const double crhs_ee21 =             crhs_ee0*crhs_ee20 + vn(1,1);
const double crhs_ee22 =             v(2,1) - vn(2,1);
const double crhs_ee23 =             crhs_ee0*crhs_ee22 + vn(2,1);
const double crhs_ee24 =             crhs_ee16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee15*(acceleration_alpha_method(0,1)*crhs_ee14 - acceleration_alpha_method(0,1) + crhs_ee13*crhs_ee18)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee15*(acceleration_alpha_method(1,1)*crhs_ee14 - acceleration_alpha_method(1,1) + crhs_ee13*crhs_ee20)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee15*(acceleration_alpha_method(2,1)*crhs_ee14 - acceleration_alpha_method(2,1) + crhs_ee13*crhs_ee22)) + crhs_ee2*(DN(0,0)*crhs_ee19 + DN(1,0)*crhs_ee21 + DN(2,0)*crhs_ee23) + crhs_ee9*(DN(0,1)*crhs_ee19 + DN(1,1)*crhs_ee21 + DN(2,1)*crhs_ee23)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee17 - DNenr(0,1)*crhs_ee24 + Nenr[0]*crhs_ee1;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee17 - DNenr(1,1)*crhs_ee24 + Nenr[1]*crhs_ee1;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee17 - DNenr(2,1)*crhs_ee24 + Nenr[2]*crhs_ee1;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double dt = rData.DeltaTime;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));
    const double alpha_f=1/(1-max_spectral_radius);

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const BoundedMatrix<double,4,3> vconv = (vn-vmeshn)+ alpha_f*((v-vmesh)-(vn-vmeshn));

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    const double volume_error_time_ratio = rData.VolumeError;

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 =             1.0/(rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 =             cV3*rho;
const double cV5 =             cV4*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV6 =             cV4*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV7 =             cV4*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV8 =             cV4*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
            V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV5;
            V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV5;
            V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV5;
            V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV5;
            V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV5;
            V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV5;
            V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV5;
            V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV5;
            V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV5;
            V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV5;
            V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV5;
            V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV5;
            V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV6;
            V(4,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV6;
            V(4,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV6;
            V(4,3)=-DN(1,0)*Nenr[3] + DNenr(3,0)*cV6;
            V(5,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV6;
            V(5,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV6;
            V(5,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV6;
            V(5,3)=-DN(1,1)*Nenr[3] + DNenr(3,1)*cV6;
            V(6,0)=-DN(1,2)*Nenr[0] + DNenr(0,2)*cV6;
            V(6,1)=-DN(1,2)*Nenr[1] + DNenr(1,2)*cV6;
            V(6,2)=-DN(1,2)*Nenr[2] + DNenr(2,2)*cV6;
            V(6,3)=-DN(1,2)*Nenr[3] + DNenr(3,2)*cV6;
            V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV7;
            V(8,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV7;
            V(8,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV7;
            V(8,3)=-DN(2,0)*Nenr[3] + DNenr(3,0)*cV7;
            V(9,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV7;
            V(9,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV7;
            V(9,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV7;
            V(9,3)=-DN(2,1)*Nenr[3] + DNenr(3,1)*cV7;
            V(10,0)=-DN(2,2)*Nenr[0] + DNenr(0,2)*cV7;
            V(10,1)=-DN(2,2)*Nenr[1] + DNenr(1,2)*cV7;
            V(10,2)=-DN(2,2)*Nenr[2] + DNenr(2,2)*cV7;
            V(10,3)=-DN(2,2)*Nenr[3] + DNenr(3,2)*cV7;
            V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] + DNenr(0,0)*cV8;
            V(12,1)=-DN(3,0)*Nenr[1] + DNenr(1,0)*cV8;
            V(12,2)=-DN(3,0)*Nenr[2] + DNenr(2,0)*cV8;
            V(12,3)=-DN(3,0)*Nenr[3] + DNenr(3,0)*cV8;
            V(13,0)=-DN(3,1)*Nenr[0] + DNenr(0,1)*cV8;
            V(13,1)=-DN(3,1)*Nenr[1] + DNenr(1,1)*cV8;
            V(13,2)=-DN(3,1)*Nenr[2] + DNenr(2,1)*cV8;
            V(13,3)=-DN(3,1)*Nenr[3] + DNenr(3,1)*cV8;
            V(14,0)=-DN(3,2)*Nenr[0] + DNenr(0,2)*cV8;
            V(14,1)=-DN(3,2)*Nenr[1] + DNenr(1,2)*cV8;
            V(14,2)=-DN(3,2)*Nenr[2] + DNenr(2,2)*cV8;
            V(14,3)=-DN(3,2)*Nenr[3] + DNenr(3,2)*cV8;
            V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             1.0/(max_spectral_radius + 1);
const double cH1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH3 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH4 =             1.0/dt;
const double cH5 =             0.5*cH4*(max_spectral_radius - 3)/(cH0*(0.5*max_spectral_radius - 1.5) + cH0 - 0.5);
const double cH6 =             1.0/(cH4*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH1, 2) + pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 =             cH6*rho;
const double cH8 =             cH7*(DN(0,0)*cH1 + DN(0,1)*cH2 + DN(0,2)*cH3 + N[0]*cH5);
const double cH9 =             cH7*(DN(1,0)*cH1 + DN(1,1)*cH2 + DN(1,2)*cH3 + N[1]*cH5);
const double cH10 =             cH7*(DN(2,0)*cH1 + DN(2,1)*cH2 + DN(2,2)*cH3 + N[2]*cH5);
const double cH11 =             cH7*(DN(3,0)*cH1 + DN(3,1)*cH2 + DN(3,2)*cH3 + N[3]*cH5);
            H(0,0)=cH0*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH8);
            H(0,1)=cH0*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH8);
            H(0,2)=cH0*(DN(0,2)*Nenr[0] + DNenr(0,2)*cH8);
            H(0,3)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=cH0*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH9);
            H(0,5)=cH0*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH9);
            H(0,6)=cH0*(DN(1,2)*Nenr[0] + DNenr(0,2)*cH9);
            H(0,7)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=cH0*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH10);
            H(0,9)=cH0*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH10);
            H(0,10)=cH0*(DN(2,2)*Nenr[0] + DNenr(0,2)*cH10);
            H(0,11)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=cH0*(DN(3,0)*Nenr[0] + DNenr(0,0)*cH11);
            H(0,13)=cH0*(DN(3,1)*Nenr[0] + DNenr(0,1)*cH11);
            H(0,14)=cH0*(DN(3,2)*Nenr[0] + DNenr(0,2)*cH11);
            H(0,15)=cH6*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=cH0*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH8);
            H(1,1)=cH0*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH8);
            H(1,2)=cH0*(DN(0,2)*Nenr[1] + DNenr(1,2)*cH8);
            H(1,3)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=cH0*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH9);
            H(1,5)=cH0*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH9);
            H(1,6)=cH0*(DN(1,2)*Nenr[1] + DNenr(1,2)*cH9);
            H(1,7)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=cH0*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH10);
            H(1,9)=cH0*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH10);
            H(1,10)=cH0*(DN(2,2)*Nenr[1] + DNenr(1,2)*cH10);
            H(1,11)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=cH0*(DN(3,0)*Nenr[1] + DNenr(1,0)*cH11);
            H(1,13)=cH0*(DN(3,1)*Nenr[1] + DNenr(1,1)*cH11);
            H(1,14)=cH0*(DN(3,2)*Nenr[1] + DNenr(1,2)*cH11);
            H(1,15)=cH6*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=cH0*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH8);
            H(2,1)=cH0*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH8);
            H(2,2)=cH0*(DN(0,2)*Nenr[2] + DNenr(2,2)*cH8);
            H(2,3)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=cH0*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH9);
            H(2,5)=cH0*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH9);
            H(2,6)=cH0*(DN(1,2)*Nenr[2] + DNenr(2,2)*cH9);
            H(2,7)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=cH0*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH10);
            H(2,9)=cH0*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH10);
            H(2,10)=cH0*(DN(2,2)*Nenr[2] + DNenr(2,2)*cH10);
            H(2,11)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=cH0*(DN(3,0)*Nenr[2] + DNenr(2,0)*cH11);
            H(2,13)=cH0*(DN(3,1)*Nenr[2] + DNenr(2,1)*cH11);
            H(2,14)=cH0*(DN(3,2)*Nenr[2] + DNenr(2,2)*cH11);
            H(2,15)=cH6*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=cH0*(DN(0,0)*Nenr[3] + DNenr(3,0)*cH8);
            H(3,1)=cH0*(DN(0,1)*Nenr[3] + DNenr(3,1)*cH8);
            H(3,2)=cH0*(DN(0,2)*Nenr[3] + DNenr(3,2)*cH8);
            H(3,3)=cH6*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=cH0*(DN(1,0)*Nenr[3] + DNenr(3,0)*cH9);
            H(3,5)=cH0*(DN(1,1)*Nenr[3] + DNenr(3,1)*cH9);
            H(3,6)=cH0*(DN(1,2)*Nenr[3] + DNenr(3,2)*cH9);
            H(3,7)=cH6*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=cH0*(DN(2,0)*Nenr[3] + DNenr(3,0)*cH10);
            H(3,9)=cH0*(DN(2,1)*Nenr[3] + DNenr(3,1)*cH10);
            H(3,10)=cH0*(DN(2,2)*Nenr[3] + DNenr(3,2)*cH10);
            H(3,11)=cH6*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=cH0*(DN(3,0)*Nenr[3] + DNenr(3,0)*cH11);
            H(3,13)=cH0*(DN(3,1)*Nenr[3] + DNenr(3,1)*cH11);
            H(3,14)=cH0*(DN(3,2)*Nenr[3] + DNenr(3,2)*cH11);
            H(3,15)=cH6*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


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


    const double crhs_ee0 =             1.0/(max_spectral_radius + 1);
const double crhs_ee1 =             -crhs_ee0*(DN(0,0)*max_spectral_radius*vn(0,0) + DN(0,0)*v(0,0) + DN(0,1)*max_spectral_radius*vn(0,1) + DN(0,1)*v(0,1) + DN(0,2)*max_spectral_radius*vn(0,2) + DN(0,2)*v(0,2) + DN(1,0)*max_spectral_radius*vn(1,0) + DN(1,0)*v(1,0) + DN(1,1)*max_spectral_radius*vn(1,1) + DN(1,1)*v(1,1) + DN(1,2)*max_spectral_radius*vn(1,2) + DN(1,2)*v(1,2) + DN(2,0)*max_spectral_radius*vn(2,0) + DN(2,0)*v(2,0) + DN(2,1)*max_spectral_radius*vn(2,1) + DN(2,1)*v(2,1) + DN(2,2)*max_spectral_radius*vn(2,2) + DN(2,2)*v(2,2) + DN(3,0)*max_spectral_radius*vn(3,0) + DN(3,0)*v(3,0) + DN(3,1)*max_spectral_radius*vn(3,1) + DN(3,1)*v(3,1) + DN(3,2)*max_spectral_radius*vn(3,2) + DN(3,2)*v(3,2)) + volume_error_time_ratio;
const double crhs_ee2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee3 =             v(0,0) - vn(0,0);
const double crhs_ee4 =             crhs_ee0*crhs_ee3 + vn(0,0);
const double crhs_ee5 =             v(1,0) - vn(1,0);
const double crhs_ee6 =             crhs_ee0*crhs_ee5 + vn(1,0);
const double crhs_ee7 =             v(2,0) - vn(2,0);
const double crhs_ee8 =             crhs_ee0*crhs_ee7 + vn(2,0);
const double crhs_ee9 =             v(3,0) - vn(3,0);
const double crhs_ee10 =             crhs_ee0*crhs_ee9 + vn(3,0);
const double crhs_ee11 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee12 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee13 =             1.0/dt;
const double crhs_ee14 =             crhs_ee0*(1.5 - 0.5*max_spectral_radius) - crhs_ee0;
const double crhs_ee15 =             1.0/(crhs_ee14 + 0.5);
const double crhs_ee16 =             crhs_ee13*crhs_ee15;
const double crhs_ee17 =             crhs_ee15*(crhs_ee14 - 0.5);
const double crhs_ee18 =             0.5*crhs_ee0*(3 - max_spectral_radius);
const double crhs_ee19 =             1.0/(crhs_ee13*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee11, 2) + pow(crhs_ee12, 2) + pow(crhs_ee2, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee20 =             crhs_ee19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,0) - fn(0,0)) + fn(0,0)) + N[1]*(crhs_ee0*(f(1,0) - fn(1,0)) + fn(1,0)) + N[2]*(crhs_ee0*(f(2,0) - fn(2,0)) + fn(2,0)) + N[3]*(crhs_ee0*(f(3,0) - fn(3,0)) + fn(3,0))) + rho*(N[0]*(acceleration_alpha_method(0,0) + crhs_ee18*(acceleration_alpha_method(0,0)*crhs_ee17 - acceleration_alpha_method(0,0) + crhs_ee16*crhs_ee3)) + N[1]*(acceleration_alpha_method(1,0) + crhs_ee18*(acceleration_alpha_method(1,0)*crhs_ee17 - acceleration_alpha_method(1,0) + crhs_ee16*crhs_ee5)) + N[2]*(acceleration_alpha_method(2,0) + crhs_ee18*(acceleration_alpha_method(2,0)*crhs_ee17 - acceleration_alpha_method(2,0) + crhs_ee16*crhs_ee7)) + N[3]*(acceleration_alpha_method(3,0) + crhs_ee18*(acceleration_alpha_method(3,0)*crhs_ee17 - acceleration_alpha_method(3,0) + crhs_ee16*crhs_ee9)) + crhs_ee11*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee8 + DN(3,1)*crhs_ee10) + crhs_ee12*(DN(0,2)*crhs_ee4 + DN(1,2)*crhs_ee6 + DN(2,2)*crhs_ee8 + DN(3,2)*crhs_ee10) + crhs_ee2*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee8 + DN(3,0)*crhs_ee10)));
const double crhs_ee21 =             v(0,1) - vn(0,1);
const double crhs_ee22 =             crhs_ee0*crhs_ee21 + vn(0,1);
const double crhs_ee23 =             v(1,1) - vn(1,1);
const double crhs_ee24 =             crhs_ee0*crhs_ee23 + vn(1,1);
const double crhs_ee25 =             v(2,1) - vn(2,1);
const double crhs_ee26 =             crhs_ee0*crhs_ee25 + vn(2,1);
const double crhs_ee27 =             v(3,1) - vn(3,1);
const double crhs_ee28 =             crhs_ee0*crhs_ee27 + vn(3,1);
const double crhs_ee29 =             crhs_ee19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,1) - fn(0,1)) + fn(0,1)) + N[1]*(crhs_ee0*(f(1,1) - fn(1,1)) + fn(1,1)) + N[2]*(crhs_ee0*(f(2,1) - fn(2,1)) + fn(2,1)) + N[3]*(crhs_ee0*(f(3,1) - fn(3,1)) + fn(3,1))) + rho*(N[0]*(acceleration_alpha_method(0,1) + crhs_ee18*(acceleration_alpha_method(0,1)*crhs_ee17 - acceleration_alpha_method(0,1) + crhs_ee16*crhs_ee21)) + N[1]*(acceleration_alpha_method(1,1) + crhs_ee18*(acceleration_alpha_method(1,1)*crhs_ee17 - acceleration_alpha_method(1,1) + crhs_ee16*crhs_ee23)) + N[2]*(acceleration_alpha_method(2,1) + crhs_ee18*(acceleration_alpha_method(2,1)*crhs_ee17 - acceleration_alpha_method(2,1) + crhs_ee16*crhs_ee25)) + N[3]*(acceleration_alpha_method(3,1) + crhs_ee18*(acceleration_alpha_method(3,1)*crhs_ee17 - acceleration_alpha_method(3,1) + crhs_ee16*crhs_ee27)) + crhs_ee11*(DN(0,1)*crhs_ee22 + DN(1,1)*crhs_ee24 + DN(2,1)*crhs_ee26 + DN(3,1)*crhs_ee28) + crhs_ee12*(DN(0,2)*crhs_ee22 + DN(1,2)*crhs_ee24 + DN(2,2)*crhs_ee26 + DN(3,2)*crhs_ee28) + crhs_ee2*(DN(0,0)*crhs_ee22 + DN(1,0)*crhs_ee24 + DN(2,0)*crhs_ee26 + DN(3,0)*crhs_ee28)));
const double crhs_ee30 =             v(0,2) - vn(0,2);
const double crhs_ee31 =             crhs_ee0*crhs_ee30 + vn(0,2);
const double crhs_ee32 =             v(1,2) - vn(1,2);
const double crhs_ee33 =             crhs_ee0*crhs_ee32 + vn(1,2);
const double crhs_ee34 =             v(2,2) - vn(2,2);
const double crhs_ee35 =             crhs_ee0*crhs_ee34 + vn(2,2);
const double crhs_ee36 =             v(3,2) - vn(3,2);
const double crhs_ee37 =             crhs_ee0*crhs_ee36 + vn(3,2);
const double crhs_ee38 =             crhs_ee19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee0*(f(0,2) - fn(0,2)) + fn(0,2)) + N[1]*(crhs_ee0*(f(1,2) - fn(1,2)) + fn(1,2)) + N[2]*(crhs_ee0*(f(2,2) - fn(2,2)) + fn(2,2)) + N[3]*(crhs_ee0*(f(3,2) - fn(3,2)) + fn(3,2))) + rho*(N[0]*(acceleration_alpha_method(0,2) + crhs_ee18*(acceleration_alpha_method(0,2)*crhs_ee17 - acceleration_alpha_method(0,2) + crhs_ee16*crhs_ee30)) + N[1]*(acceleration_alpha_method(1,2) + crhs_ee18*(acceleration_alpha_method(1,2)*crhs_ee17 - acceleration_alpha_method(1,2) + crhs_ee16*crhs_ee32)) + N[2]*(acceleration_alpha_method(2,2) + crhs_ee18*(acceleration_alpha_method(2,2)*crhs_ee17 - acceleration_alpha_method(2,2) + crhs_ee16*crhs_ee34)) + N[3]*(acceleration_alpha_method(3,2) + crhs_ee18*(acceleration_alpha_method(3,2)*crhs_ee17 - acceleration_alpha_method(3,2) + crhs_ee16*crhs_ee36)) + crhs_ee11*(DN(0,1)*crhs_ee31 + DN(1,1)*crhs_ee33 + DN(2,1)*crhs_ee35 + DN(3,1)*crhs_ee37) + crhs_ee12*(DN(0,2)*crhs_ee31 + DN(1,2)*crhs_ee33 + DN(2,2)*crhs_ee35 + DN(3,2)*crhs_ee37) + crhs_ee2*(DN(0,0)*crhs_ee31 + DN(1,0)*crhs_ee33 + DN(2,0)*crhs_ee35 + DN(3,0)*crhs_ee37)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee20 - DNenr(0,1)*crhs_ee29 - DNenr(0,2)*crhs_ee38 + Nenr[0]*crhs_ee1;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee20 - DNenr(1,1)*crhs_ee29 - DNenr(1,2)*crhs_ee38 + Nenr[1]*crhs_ee1;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee20 - DNenr(2,1)*crhs_ee29 - DNenr(2,2)*crhs_ee38 + Nenr[2]*crhs_ee1;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee20 - DNenr(3,1)*crhs_ee29 - DNenr(3,2)*crhs_ee38 + Nenr[3]*crhs_ee1;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::ComputeSplitting(
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
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

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
void TwoFluidNavierStokesAlphaMethod<TElementData>::ComputeSplitInterface(
    const TElementData &rData,
    MatrixType& rInterfaceShapeFunctionNeg,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Vector& rInterfaceWeightsNeg,
    std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
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
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    pModifiedShapeFunctions->ComputeNegativeSideInterfaceAreaNormals(
        rInterfaceNormalsNeg,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); ++gp){
        const double normal_norm = norm_2(rInterfaceNormalsNeg[gp]);
        rInterfaceNormalsNeg[gp] /= normal_norm;
    }

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesAlphaMethod< TwoFluidNavierStokesAlphaMethodData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
   return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesAlphaMethod< TwoFluidNavierStokesAlphaMethodData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::CalculateCurvatureOnInterfaceGaussPoints(
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::SurfaceTension(
    const double SurfaceTensionCoefficient,
    const Vector& rCurvature,
    const Vector& rInterfaceWeights,
    const Matrix& rInterfaceShapeFunctions,
    const std::vector<array_1d<double,3>>& rInterfaceNormalsNeg,
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::PressureGradientStabilization(
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
    const auto &v=rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &acceleration_alpha_method=rData.AccelerationAlphaMethod;
    const double max_spectral_radius=rData.MaxSprectraRadius;
    const double alpha_m=0.5*((3-rData.MaxSprectraRadius)/(1+rData.MaxSprectraRadius));

    const auto vmesh=rData.MeshVelocity;
    const auto vmeshn=rData.MeshVelocityOldStep;
    const BoundedMatrix<double,NumNodes,Dim> v_convection = (vn-vmeshn)+ alpha_m*((vn-vmeshn)-(v-vmesh));

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
            density * 1.0 / (dyn_tau * density / (0.5*dt) + stab_c1 * viscosity / h_elem / h_elem +
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::CondenseEnrichmentWithContinuity(
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
                    double penalty_coeff = max_diag * 0.001;
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::CondenseEnrichment(
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::AddSurfaceTensionContribution(
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
    std::vector<array_1d<double,3>> int_normals_neg;
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
void TwoFluidNavierStokesAlphaMethod<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void TwoFluidNavierStokesAlphaMethod<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo )
{
    if (rVariable == DIVERGENCE){

        const auto& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
        const unsigned int num_gauss = IntegrationPoints.size();

        if (rValues.size() != num_gauss){
            rValues.resize(num_gauss);
        }

        Vector gauss_pts_jacobian_determinant = ZeroVector(num_gauss);
        GeometryData::ShapeFunctionsGradientsType DN_DX;
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, gauss_pts_jacobian_determinant, GeometryData::IntegrationMethod::GI_GAUSS_2);

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

template class TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<2, 3>>;
template class TwoFluidNavierStokesAlphaMethod<TwoFluidNavierStokesAlphaMethodData<3, 4>>;

} // namespace Kratos
