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

#include "two_fluid_navier_stokes_CN.h"
#include "custom_utilities/two_fluid_navier_stokes_CN_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::TwoFluidNavierStokesCN(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokesCN<TElementData>::~TwoFluidNavierStokesCN() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesCN<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesCN>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesCN<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesCN>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateLocalSystem(
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
        KRATOS_ERROR << "TwoFluidNavierStokesCN is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokesCN<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
std::string TwoFluidNavierStokesCN<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokesCN" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::PrintInfo(
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
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::UpdateIntegrationPointData(
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
void TwoFluidNavierStokesCN<TElementData>::UpdateIntegrationPointData(
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
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::CalculateStrainRate(TwoFluidNavierStokesCNData<2, 3>& rData) const
{
    const double theta = rData.theta;
    const BoundedMatrix<double,3,2> mid_step_velocity = theta*rData.Velocity + (1-theta)*rData.Velocity_OldStep1;
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(3);
    for (unsigned int i = 0; i < 3; i++) {
        r_strain_rate[0] += rDNDX(i,0)*mid_step_velocity(i,0);
        r_strain_rate[1] += rDNDX(i,1)*mid_step_velocity(i,1);
        r_strain_rate[2] += rDNDX(i,0)*mid_step_velocity(i,1) + rDNDX(i,1)*mid_step_velocity(i,0);
    }
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::CalculateStrainRate(TwoFluidNavierStokesCNData<3, 4>& rData) const
{
    const double theta = rData.theta;
    const BoundedMatrix<double,4,3> mid_step_velocity = theta*rData.Velocity + (1-theta)*rData.Velocity_OldStep1;;
    auto& rDNDX = rData.DN_DX;
    auto& r_strain_rate = rData.StrainRate;
    noalias(r_strain_rate) = ZeroVector(6);
    for (unsigned int i = 0; i < 4; i++) {
        r_strain_rate[0] += rDNDX(i,0)*mid_step_velocity(i,0);
        r_strain_rate[1] += rDNDX(i,1)*mid_step_velocity(i,1);
        r_strain_rate[2] += rDNDX(i,2)*mid_step_velocity(i,2);
        r_strain_rate[3] += rDNDX(i,0)*mid_step_velocity(i,1) + rDNDX(i,1)*mid_step_velocity(i,0);
        r_strain_rate[4] += rDNDX(i,1)*mid_step_velocity(i,2) + rDNDX(i,2)*mid_step_velocity(i,1);
        r_strain_rate[5] += rDNDX(i,0)*mid_step_velocity(i,2) + rDNDX(i,2)*mid_step_velocity(i,0);
    }
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

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
const double clhs4 =             mu*stab_c1;
const double clhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs7 =             rho*sqrt(pow(clhs5, 2) + pow(clhs6, 2));
const double clhs8 =             clhs4/stab_c2 + clhs7*h;
const double clhs9 =             DN(0,0)*clhs5 + DN(0,1)*clhs6;
const double clhs10 =             N[0]*rho;
const double clhs11 =             1.0/dt;
const double clhs12 =             clhs11*rho;
const double clhs13 =             1.0*clhs12;
const double clhs14 =             1.0*clhs11;
const double clhs15 =             N[0]*clhs14 + clhs9;
const double clhs16 =             1.0/(clhs12*dyn_tau + clhs4/pow(h, 2) + clhs7*stab_c2/h);
const double clhs17 =             clhs16*pow(rho, 2);
const double clhs18 =             clhs17*clhs9;
const double clhs19 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs20 =             clhs17*clhs19;
const double clhs21 =             N[0]*clhs20;
const double clhs22 =             pow(N[0], 2)*clhs13 + clhs10*clhs9 + clhs15*clhs18 + clhs15*clhs21;
const double clhs23 =             C(0,1)*DN(0,1) + clhs1;
const double clhs24 =             C(1,2)*DN(0,1);
const double clhs25 =             C(2,2)*DN(0,0) + clhs24;
const double clhs26 =             DN(0,0)*clhs8;
const double clhs27 =             DN(0,1)*clhs26;
const double clhs28 =             N[0]*theta;
const double clhs29 =             clhs16*clhs19;
const double clhs30 =             clhs16*rho;
const double clhs31 =             clhs30*clhs9;
const double clhs32 =             clhs10*clhs29 - clhs28 + clhs31;
const double clhs33 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs34 =             C(0,2)*DN(1,0);
const double clhs35 =             C(2,2)*DN(1,1) + clhs34;
const double clhs36 =             DN(0,0)*DN(1,0);
const double clhs37 =             N[1]*clhs14;
const double clhs38 =             clhs10*clhs37;
const double clhs39 =             clhs36*clhs8 + clhs38;
const double clhs40 =             DN(1,0)*clhs5 + DN(1,1)*clhs6;
const double clhs41 =             clhs37 + clhs40;
const double clhs42 =             clhs10*clhs40 + clhs18*clhs41 + clhs21*clhs41;
const double clhs43 =             C(0,1)*DN(1,1) + clhs34;
const double clhs44 =             C(1,2)*DN(1,1);
const double clhs45 =             C(2,2)*DN(1,0) + clhs44;
const double clhs46 =             DN(1,1)*clhs26;
const double clhs47 =             N[1]*theta;
const double clhs48 =             DN(1,0)*N[0];
const double clhs49 =             clhs19*clhs30;
const double clhs50 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs51 =             C(0,2)*DN(2,0);
const double clhs52 =             C(2,2)*DN(2,1) + clhs51;
const double clhs53 =             DN(0,0)*DN(2,0);
const double clhs54 =             N[2]*clhs14;
const double clhs55 =             clhs10*clhs54;
const double clhs56 =             clhs53*clhs8 + clhs55;
const double clhs57 =             DN(2,0)*clhs5 + DN(2,1)*clhs6;
const double clhs58 =             clhs54 + clhs57;
const double clhs59 =             clhs10*clhs57 + clhs18*clhs58 + clhs21*clhs58;
const double clhs60 =             C(0,1)*DN(2,1) + clhs51;
const double clhs61 =             C(1,2)*DN(2,1);
const double clhs62 =             C(2,2)*DN(2,0) + clhs61;
const double clhs63 =             DN(2,1)*clhs26;
const double clhs64 =             N[2]*theta;
const double clhs65 =             DN(2,0)*N[0];
const double clhs66 =             C(0,1)*DN(0,0) + clhs24;
const double clhs67 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs68 =             pow(DN(0,1), 2);
const double clhs69 =             C(0,1)*DN(1,0) + clhs44;
const double clhs70 =             DN(0,1)*clhs8;
const double clhs71 =             DN(1,0)*clhs70;
const double clhs72 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs73 =             DN(0,1)*DN(1,1);
const double clhs74 =             clhs38 + clhs73*clhs8;
const double clhs75 =             DN(1,1)*N[0];
const double clhs76 =             C(0,1)*DN(2,0) + clhs61;
const double clhs77 =             DN(2,0)*clhs70;
const double clhs78 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs79 =             DN(0,1)*DN(2,1);
const double clhs80 =             clhs55 + clhs79*clhs8;
const double clhs81 =             DN(2,1)*N[0];
const double clhs82 =             clhs15*clhs30;
const double clhs83 =             N[0] + clhs82;
const double clhs84 =             clhs30*clhs41;
const double clhs85 =             clhs16*(clhs36 + clhs73);
const double clhs86 =             clhs30*clhs58;
const double clhs87 =             clhs16*(clhs53 + clhs79);
const double clhs88 =             N[1]*rho;
const double clhs89 =             clhs17*clhs40;
const double clhs90 =             N[1]*clhs20;
const double clhs91 =             clhs15*clhs89 + clhs15*clhs90 + clhs88*clhs9;
const double clhs92 =             DN(0,0)*N[1];
const double clhs93 =             clhs30*clhs40;
const double clhs94 =             pow(DN(1,0), 2);
const double clhs95 =             pow(N[1], 2)*clhs13 + clhs40*clhs88 + clhs41*clhs89 + clhs41*clhs90;
const double clhs96 =             DN(1,0)*clhs8;
const double clhs97 =             DN(1,1)*clhs96;
const double clhs98 =             clhs29*clhs88 - clhs47 + clhs93;
const double clhs99 =             DN(1,0)*DN(2,0);
const double clhs100 =             N[1]*N[2]*clhs13;
const double clhs101 =             clhs100 + clhs8*clhs99;
const double clhs102 =             clhs57*clhs88 + clhs58*clhs89 + clhs58*clhs90;
const double clhs103 =             DN(2,1)*clhs96;
const double clhs104 =             DN(2,0)*N[1];
const double clhs105 =             DN(0,1)*N[1];
const double clhs106 =             pow(DN(1,1), 2);
const double clhs107 =             DN(2,0)*clhs8;
const double clhs108 =             DN(1,1)*clhs107;
const double clhs109 =             DN(1,1)*DN(2,1);
const double clhs110 =             clhs100 + clhs109*clhs8;
const double clhs111 =             DN(2,1)*N[1];
const double clhs112 =             N[1] + clhs84;
const double clhs113 =             clhs16*(clhs109 + clhs99);
const double clhs114 =             N[2]*rho;
const double clhs115 =             clhs17*clhs57;
const double clhs116 =             N[2]*clhs20;
const double clhs117 =             clhs114*clhs9 + clhs115*clhs15 + clhs116*clhs15;
const double clhs118 =             DN(0,0)*N[2];
const double clhs119 =             clhs30*clhs57;
const double clhs120 =             clhs114*clhs40 + clhs115*clhs41 + clhs116*clhs41;
const double clhs121 =             DN(1,0)*N[2];
const double clhs122 =             pow(DN(2,0), 2);
const double clhs123 =             pow(N[2], 2)*clhs13 + clhs114*clhs57 + clhs115*clhs58 + clhs116*clhs58;
const double clhs124 =             DN(2,1)*clhs107;
const double clhs125 =             clhs114*clhs29 + clhs119 - clhs64;
const double clhs126 =             DN(0,1)*N[2];
const double clhs127 =             DN(1,1)*N[2];
const double clhs128 =             pow(DN(2,1), 2);
const double clhs129 =             N[2] + clhs86;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs22 + clhs3*clhs8;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + clhs27;
            lhs(0,2)=DN(0,0)*clhs32;
            lhs(0,3)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + clhs39 + clhs42;
            lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46;
            lhs(0,5)=-DN(0,0)*clhs47 + DN(1,0)*clhs31 + clhs48*clhs49;
            lhs(0,6)=DN(0,0)*clhs50 + DN(0,1)*clhs52 + clhs56 + clhs59;
            lhs(0,7)=DN(0,0)*clhs60 + DN(0,1)*clhs62 + clhs63;
            lhs(0,8)=-DN(0,0)*clhs64 + DN(2,0)*clhs31 + clhs49*clhs65;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs66 + clhs27;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs67 + clhs22 + clhs68*clhs8;
            lhs(1,2)=DN(0,1)*clhs32;
            lhs(1,3)=DN(0,0)*clhs35 + DN(0,1)*clhs69 + clhs71;
            lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs72 + clhs42 + clhs74;
            lhs(1,5)=-DN(0,1)*clhs47 + DN(1,1)*clhs31 + clhs49*clhs75;
            lhs(1,6)=DN(0,0)*clhs52 + DN(0,1)*clhs76 + clhs77;
            lhs(1,7)=DN(0,0)*clhs62 + DN(0,1)*clhs78 + clhs59 + clhs80;
            lhs(1,8)=-DN(0,1)*clhs64 + DN(2,1)*clhs31 + clhs49*clhs81;
            lhs(2,0)=DN(0,0)*clhs83;
            lhs(2,1)=DN(0,1)*clhs83;
            lhs(2,2)=clhs16*(clhs3 + clhs68);
            lhs(2,3)=DN(0,0)*clhs84 + clhs48;
            lhs(2,4)=DN(0,1)*clhs84 + clhs75;
            lhs(2,5)=clhs85;
            lhs(2,6)=DN(0,0)*clhs86 + clhs65;
            lhs(2,7)=DN(0,1)*clhs86 + clhs81;
            lhs(2,8)=clhs87;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs39 + clhs91;
            lhs(3,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + clhs71;
            lhs(3,2)=DN(0,0)*clhs93 - DN(1,0)*clhs28 + clhs49*clhs92;
            lhs(3,3)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + clhs8*clhs94 + clhs95;
            lhs(3,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs97;
            lhs(3,5)=DN(1,0)*clhs98;
            lhs(3,6)=DN(1,0)*clhs50 + DN(1,1)*clhs52 + clhs101 + clhs102;
            lhs(3,7)=DN(1,0)*clhs60 + DN(1,1)*clhs62 + clhs103;
            lhs(3,8)=-DN(1,0)*clhs64 + DN(2,0)*clhs93 + clhs104*clhs49;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs66 + clhs46;
            lhs(4,1)=DN(1,0)*clhs25 + DN(1,1)*clhs67 + clhs74 + clhs91;
            lhs(4,2)=DN(0,1)*clhs93 - DN(1,1)*clhs28 + clhs105*clhs49;
            lhs(4,3)=DN(1,0)*clhs35 + DN(1,1)*clhs69 + clhs97;
            lhs(4,4)=DN(1,0)*clhs45 + DN(1,1)*clhs72 + clhs106*clhs8 + clhs95;
            lhs(4,5)=DN(1,1)*clhs98;
            lhs(4,6)=DN(1,0)*clhs52 + DN(1,1)*clhs76 + clhs108;
            lhs(4,7)=DN(1,0)*clhs62 + DN(1,1)*clhs78 + clhs102 + clhs110;
            lhs(4,8)=-DN(1,1)*clhs64 + DN(2,1)*clhs93 + clhs111*clhs49;
            lhs(5,0)=DN(1,0)*clhs82 + clhs92;
            lhs(5,1)=DN(1,1)*clhs82 + clhs105;
            lhs(5,2)=clhs85;
            lhs(5,3)=DN(1,0)*clhs112;
            lhs(5,4)=DN(1,1)*clhs112;
            lhs(5,5)=clhs16*(clhs106 + clhs94);
            lhs(5,6)=DN(1,0)*clhs86 + clhs104;
            lhs(5,7)=DN(1,1)*clhs86 + clhs111;
            lhs(5,8)=clhs113;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs117 + clhs56;
            lhs(6,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + clhs77;
            lhs(6,2)=DN(0,0)*clhs119 - DN(2,0)*clhs28 + clhs118*clhs49;
            lhs(6,3)=DN(2,0)*clhs33 + DN(2,1)*clhs35 + clhs101 + clhs120;
            lhs(6,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs108;
            lhs(6,5)=DN(1,0)*clhs119 - DN(2,0)*clhs47 + clhs121*clhs49;
            lhs(6,6)=DN(2,0)*clhs50 + DN(2,1)*clhs52 + clhs122*clhs8 + clhs123;
            lhs(6,7)=DN(2,0)*clhs60 + DN(2,1)*clhs62 + clhs124;
            lhs(6,8)=DN(2,0)*clhs125;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs66 + clhs63;
            lhs(7,1)=DN(2,0)*clhs25 + DN(2,1)*clhs67 + clhs117 + clhs80;
            lhs(7,2)=DN(0,1)*clhs119 - DN(2,1)*clhs28 + clhs126*clhs49;
            lhs(7,3)=DN(2,0)*clhs35 + DN(2,1)*clhs69 + clhs103;
            lhs(7,4)=DN(2,0)*clhs45 + DN(2,1)*clhs72 + clhs110 + clhs120;
            lhs(7,5)=DN(1,1)*clhs119 - DN(2,1)*clhs47 + clhs127*clhs49;
            lhs(7,6)=DN(2,0)*clhs52 + DN(2,1)*clhs76 + clhs124;
            lhs(7,7)=DN(2,0)*clhs62 + DN(2,1)*clhs78 + clhs123 + clhs128*clhs8;
            lhs(7,8)=DN(2,1)*clhs125;
            lhs(8,0)=DN(2,0)*clhs82 + clhs118;
            lhs(8,1)=DN(2,1)*clhs82 + clhs126;
            lhs(8,2)=clhs87;
            lhs(8,3)=DN(2,0)*clhs84 + clhs121;
            lhs(8,4)=DN(2,1)*clhs84 + clhs127;
            lhs(8,5)=clhs113;
            lhs(8,6)=DN(2,0)*clhs129;
            lhs(8,7)=DN(2,1)*clhs129;
            lhs(8,8)=clhs16*(clhs122 + clhs128);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

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
const double clhs6 =             mu*stab_c1;
const double clhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs10 =             rho*sqrt(pow(clhs7, 2) + pow(clhs8, 2) + pow(clhs9, 2));
const double clhs11 =             clhs10*h + clhs6/stab_c2;
const double clhs12 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs13 =             N[0]*rho;
const double clhs14 =             1.0/dt;
const double clhs15 =             clhs14*rho;
const double clhs16 =             1.0*clhs15;
const double clhs17 =             1.0*clhs14;
const double clhs18 =             N[0]*clhs17 + clhs12;
const double clhs19 =             1.0/(clhs10*stab_c2/h + clhs15*dyn_tau + clhs6/pow(h, 2));
const double clhs20 =             clhs19*pow(rho, 2);
const double clhs21 =             clhs12*clhs20;
const double clhs22 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs23 =             clhs20*clhs22;
const double clhs24 =             N[0]*clhs23;
const double clhs25 =             pow(N[0], 2)*clhs16 + clhs12*clhs13 + clhs18*clhs21 + clhs18*clhs24;
const double clhs26 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs27 =             C(1,3)*DN(0,1);
const double clhs28 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs27;
const double clhs29 =             C(3,5)*DN(0,0);
const double clhs30 =             C(4,5)*DN(0,2);
const double clhs31 =             C(1,5)*DN(0,1) + clhs29 + clhs30;
const double clhs32 =             DN(0,0)*clhs11;
const double clhs33 =             DN(0,1)*clhs32;
const double clhs34 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs35 =             C(3,4)*DN(0,1);
const double clhs36 =             C(2,3)*DN(0,2) + clhs29 + clhs35;
const double clhs37 =             C(2,5)*DN(0,2);
const double clhs38 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs37;
const double clhs39 =             DN(0,2)*clhs32;
const double clhs40 =             N[0]*theta;
const double clhs41 =             clhs19*clhs22;
const double clhs42 =             clhs19*rho;
const double clhs43 =             clhs12*clhs42;
const double clhs44 =             clhs13*clhs41 - clhs40 + clhs43;
const double clhs45 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs46 =             C(0,3)*DN(1,0);
const double clhs47 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs46;
const double clhs48 =             C(0,5)*DN(1,0);
const double clhs49 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs48;
const double clhs50 =             DN(0,0)*DN(1,0);
const double clhs51 =             N[1]*clhs17;
const double clhs52 =             clhs13*clhs51;
const double clhs53 =             clhs11*clhs50 + clhs52;
const double clhs54 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs55 =             clhs51 + clhs54;
const double clhs56 =             clhs13*clhs54 + clhs21*clhs55 + clhs24*clhs55;
const double clhs57 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs46;
const double clhs58 =             C(1,3)*DN(1,1);
const double clhs59 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs58;
const double clhs60 =             C(3,5)*DN(1,0);
const double clhs61 =             C(4,5)*DN(1,2);
const double clhs62 =             C(1,5)*DN(1,1) + clhs60 + clhs61;
const double clhs63 =             DN(1,1)*clhs32;
const double clhs64 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs48;
const double clhs65 =             C(3,4)*DN(1,1);
const double clhs66 =             C(2,3)*DN(1,2) + clhs60 + clhs65;
const double clhs67 =             C(2,5)*DN(1,2);
const double clhs68 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs67;
const double clhs69 =             DN(1,2)*clhs32;
const double clhs70 =             N[1]*theta;
const double clhs71 =             DN(1,0)*N[0];
const double clhs72 =             clhs22*clhs42;
const double clhs73 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs74 =             C(0,3)*DN(2,0);
const double clhs75 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs74;
const double clhs76 =             C(0,5)*DN(2,0);
const double clhs77 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs76;
const double clhs78 =             DN(0,0)*DN(2,0);
const double clhs79 =             N[2]*clhs17;
const double clhs80 =             clhs13*clhs79;
const double clhs81 =             clhs11*clhs78 + clhs80;
const double clhs82 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs83 =             clhs79 + clhs82;
const double clhs84 =             clhs13*clhs82 + clhs21*clhs83 + clhs24*clhs83;
const double clhs85 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs74;
const double clhs86 =             C(1,3)*DN(2,1);
const double clhs87 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs86;
const double clhs88 =             C(3,5)*DN(2,0);
const double clhs89 =             C(4,5)*DN(2,2);
const double clhs90 =             C(1,5)*DN(2,1) + clhs88 + clhs89;
const double clhs91 =             DN(2,1)*clhs32;
const double clhs92 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs76;
const double clhs93 =             C(3,4)*DN(2,1);
const double clhs94 =             C(2,3)*DN(2,2) + clhs88 + clhs93;
const double clhs95 =             C(2,5)*DN(2,2);
const double clhs96 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs95;
const double clhs97 =             DN(2,2)*clhs32;
const double clhs98 =             N[2]*theta;
const double clhs99 =             DN(2,0)*N[0];
const double clhs100 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs101 =             C(0,3)*DN(3,0);
const double clhs102 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs101;
const double clhs103 =             C(0,5)*DN(3,0);
const double clhs104 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs103;
const double clhs105 =             DN(0,0)*DN(3,0);
const double clhs106 =             N[3]*clhs17;
const double clhs107 =             clhs106*clhs13;
const double clhs108 =             clhs105*clhs11 + clhs107;
const double clhs109 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs110 =             clhs106 + clhs109;
const double clhs111 =             clhs109*clhs13 + clhs110*clhs21 + clhs110*clhs24;
const double clhs112 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs101;
const double clhs113 =             C(1,3)*DN(3,1);
const double clhs114 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs113;
const double clhs115 =             C(3,5)*DN(3,0);
const double clhs116 =             C(4,5)*DN(3,2);
const double clhs117 =             C(1,5)*DN(3,1) + clhs115 + clhs116;
const double clhs118 =             DN(3,1)*clhs32;
const double clhs119 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs103;
const double clhs120 =             C(3,4)*DN(3,1);
const double clhs121 =             C(2,3)*DN(3,2) + clhs115 + clhs120;
const double clhs122 =             C(2,5)*DN(3,2);
const double clhs123 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs122;
const double clhs124 =             DN(3,2)*clhs32;
const double clhs125 =             N[3]*theta;
const double clhs126 =             DN(3,0)*N[0];
const double clhs127 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs27;
const double clhs128 =             C(0,4)*DN(0,0) + clhs30 + clhs35;
const double clhs129 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs130 =             C(1,4)*DN(0,1);
const double clhs131 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs130;
const double clhs132 =             pow(DN(0,1), 2);
const double clhs133 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs130;
const double clhs134 =             C(2,4)*DN(0,2);
const double clhs135 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs134;
const double clhs136 =             DN(0,1)*clhs11;
const double clhs137 =             DN(0,2)*clhs136;
const double clhs138 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs58;
const double clhs139 =             C(0,4)*DN(1,0) + clhs61 + clhs65;
const double clhs140 =             DN(1,0)*clhs136;
const double clhs141 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs142 =             C(1,4)*DN(1,1);
const double clhs143 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs142;
const double clhs144 =             DN(0,1)*DN(1,1);
const double clhs145 =             clhs11*clhs144;
const double clhs146 =             clhs52 + clhs56;
const double clhs147 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs142;
const double clhs148 =             C(2,4)*DN(1,2);
const double clhs149 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs148;
const double clhs150 =             DN(1,2)*clhs136;
const double clhs151 =             DN(1,1)*N[0];
const double clhs152 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs86;
const double clhs153 =             C(0,4)*DN(2,0) + clhs89 + clhs93;
const double clhs154 =             DN(2,0)*clhs136;
const double clhs155 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs156 =             C(1,4)*DN(2,1);
const double clhs157 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs156;
const double clhs158 =             DN(0,1)*DN(2,1);
const double clhs159 =             clhs11*clhs158;
const double clhs160 =             clhs80 + clhs84;
const double clhs161 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs156;
const double clhs162 =             C(2,4)*DN(2,2);
const double clhs163 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs162;
const double clhs164 =             DN(2,2)*clhs136;
const double clhs165 =             DN(2,1)*N[0];
const double clhs166 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs113;
const double clhs167 =             C(0,4)*DN(3,0) + clhs116 + clhs120;
const double clhs168 =             DN(3,0)*clhs136;
const double clhs169 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs170 =             C(1,4)*DN(3,1);
const double clhs171 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs170;
const double clhs172 =             DN(0,1)*DN(3,1);
const double clhs173 =             clhs11*clhs172;
const double clhs174 =             clhs107 + clhs111;
const double clhs175 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs170;
const double clhs176 =             C(2,4)*DN(3,2);
const double clhs177 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs176;
const double clhs178 =             DN(3,2)*clhs136;
const double clhs179 =             DN(3,1)*N[0];
const double clhs180 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs37;
const double clhs181 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs134;
const double clhs182 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs183 =             pow(DN(0,2), 2);
const double clhs184 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs67;
const double clhs185 =             DN(0,2)*clhs11;
const double clhs186 =             DN(1,0)*clhs185;
const double clhs187 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs148;
const double clhs188 =             DN(1,1)*clhs185;
const double clhs189 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs190 =             DN(0,2)*DN(1,2);
const double clhs191 =             clhs11*clhs190;
const double clhs192 =             DN(1,2)*N[0];
const double clhs193 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs95;
const double clhs194 =             DN(2,0)*clhs185;
const double clhs195 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs162;
const double clhs196 =             DN(2,1)*clhs185;
const double clhs197 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs198 =             DN(0,2)*DN(2,2);
const double clhs199 =             clhs11*clhs198;
const double clhs200 =             DN(2,2)*N[0];
const double clhs201 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs122;
const double clhs202 =             DN(3,0)*clhs185;
const double clhs203 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs176;
const double clhs204 =             DN(3,1)*clhs185;
const double clhs205 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs206 =             DN(0,2)*DN(3,2);
const double clhs207 =             clhs11*clhs206;
const double clhs208 =             DN(3,2)*N[0];
const double clhs209 =             clhs18*clhs42;
const double clhs210 =             N[0] + clhs209;
const double clhs211 =             clhs42*clhs55;
const double clhs212 =             clhs19*(clhs144 + clhs190 + clhs50);
const double clhs213 =             clhs42*clhs83;
const double clhs214 =             clhs19*(clhs158 + clhs198 + clhs78);
const double clhs215 =             clhs110*clhs42;
const double clhs216 =             clhs19*(clhs105 + clhs172 + clhs206);
const double clhs217 =             N[1]*rho;
const double clhs218 =             clhs20*clhs54;
const double clhs219 =             N[1]*clhs23;
const double clhs220 =             clhs12*clhs217 + clhs18*clhs218 + clhs18*clhs219;
const double clhs221 =             DN(0,0)*N[1];
const double clhs222 =             clhs42*clhs54;
const double clhs223 =             pow(DN(1,0), 2);
const double clhs224 =             pow(N[1], 2)*clhs16 + clhs217*clhs54 + clhs218*clhs55 + clhs219*clhs55;
const double clhs225 =             DN(1,0)*clhs11;
const double clhs226 =             DN(1,1)*clhs225;
const double clhs227 =             DN(1,2)*clhs225;
const double clhs228 =             clhs217*clhs41 + clhs222 - clhs70;
const double clhs229 =             DN(1,0)*DN(2,0);
const double clhs230 =             N[1]*clhs16;
const double clhs231 =             N[2]*clhs230;
const double clhs232 =             clhs11*clhs229 + clhs231;
const double clhs233 =             clhs217*clhs82 + clhs218*clhs83 + clhs219*clhs83;
const double clhs234 =             DN(2,1)*clhs225;
const double clhs235 =             DN(2,2)*clhs225;
const double clhs236 =             DN(2,0)*N[1];
const double clhs237 =             DN(1,0)*DN(3,0);
const double clhs238 =             N[3]*clhs230;
const double clhs239 =             clhs11*clhs237 + clhs238;
const double clhs240 =             clhs109*clhs217 + clhs110*clhs218 + clhs110*clhs219;
const double clhs241 =             DN(3,1)*clhs225;
const double clhs242 =             DN(3,2)*clhs225;
const double clhs243 =             DN(3,0)*N[1];
const double clhs244 =             clhs220 + clhs52;
const double clhs245 =             DN(0,1)*N[1];
const double clhs246 =             pow(DN(1,1), 2);
const double clhs247 =             DN(1,1)*clhs11;
const double clhs248 =             DN(1,2)*clhs247;
const double clhs249 =             DN(2,0)*clhs247;
const double clhs250 =             DN(1,1)*DN(2,1);
const double clhs251 =             clhs11*clhs250;
const double clhs252 =             clhs231 + clhs233;
const double clhs253 =             DN(2,2)*clhs247;
const double clhs254 =             DN(2,1)*N[1];
const double clhs255 =             DN(3,0)*clhs247;
const double clhs256 =             DN(1,1)*DN(3,1);
const double clhs257 =             clhs11*clhs256;
const double clhs258 =             clhs238 + clhs240;
const double clhs259 =             DN(3,2)*clhs247;
const double clhs260 =             DN(3,1)*N[1];
const double clhs261 =             DN(0,2)*N[1];
const double clhs262 =             pow(DN(1,2), 2);
const double clhs263 =             DN(1,2)*clhs11;
const double clhs264 =             DN(2,0)*clhs263;
const double clhs265 =             DN(2,1)*clhs263;
const double clhs266 =             DN(1,2)*DN(2,2);
const double clhs267 =             clhs11*clhs266;
const double clhs268 =             DN(2,2)*N[1];
const double clhs269 =             DN(3,0)*clhs263;
const double clhs270 =             DN(3,1)*clhs263;
const double clhs271 =             DN(1,2)*DN(3,2);
const double clhs272 =             clhs11*clhs271;
const double clhs273 =             DN(3,2)*N[1];
const double clhs274 =             N[1] + clhs211;
const double clhs275 =             clhs19*(clhs229 + clhs250 + clhs266);
const double clhs276 =             clhs19*(clhs237 + clhs256 + clhs271);
const double clhs277 =             N[2]*rho;
const double clhs278 =             clhs20*clhs82;
const double clhs279 =             N[2]*clhs23;
const double clhs280 =             clhs12*clhs277 + clhs18*clhs278 + clhs18*clhs279;
const double clhs281 =             DN(0,0)*N[2];
const double clhs282 =             clhs42*clhs82;
const double clhs283 =             clhs277*clhs54 + clhs278*clhs55 + clhs279*clhs55;
const double clhs284 =             DN(1,0)*N[2];
const double clhs285 =             pow(DN(2,0), 2);
const double clhs286 =             pow(N[2], 2)*clhs16 + clhs277*clhs82 + clhs278*clhs83 + clhs279*clhs83;
const double clhs287 =             DN(2,0)*clhs11;
const double clhs288 =             DN(2,1)*clhs287;
const double clhs289 =             DN(2,2)*clhs287;
const double clhs290 =             clhs277*clhs41 + clhs282 - clhs98;
const double clhs291 =             DN(2,0)*DN(3,0);
const double clhs292 =             N[2]*N[3]*clhs16;
const double clhs293 =             clhs11*clhs291 + clhs292;
const double clhs294 =             clhs109*clhs277 + clhs110*clhs278 + clhs110*clhs279;
const double clhs295 =             DN(3,1)*clhs287;
const double clhs296 =             DN(3,2)*clhs287;
const double clhs297 =             DN(3,0)*N[2];
const double clhs298 =             clhs280 + clhs80;
const double clhs299 =             DN(0,1)*N[2];
const double clhs300 =             clhs231 + clhs283;
const double clhs301 =             DN(1,1)*N[2];
const double clhs302 =             pow(DN(2,1), 2);
const double clhs303 =             DN(2,1)*clhs11;
const double clhs304 =             DN(2,2)*clhs303;
const double clhs305 =             DN(3,0)*clhs303;
const double clhs306 =             DN(2,1)*DN(3,1);
const double clhs307 =             clhs11*clhs306;
const double clhs308 =             clhs292 + clhs294;
const double clhs309 =             DN(3,2)*clhs303;
const double clhs310 =             DN(3,1)*N[2];
const double clhs311 =             DN(0,2)*N[2];
const double clhs312 =             DN(1,2)*N[2];
const double clhs313 =             pow(DN(2,2), 2);
const double clhs314 =             DN(2,2)*clhs11;
const double clhs315 =             DN(3,0)*clhs314;
const double clhs316 =             DN(3,1)*clhs314;
const double clhs317 =             DN(2,2)*DN(3,2);
const double clhs318 =             clhs11*clhs317;
const double clhs319 =             DN(3,2)*N[2];
const double clhs320 =             N[2] + clhs213;
const double clhs321 =             clhs19*(clhs291 + clhs306 + clhs317);
const double clhs322 =             N[3]*rho;
const double clhs323 =             clhs109*clhs20;
const double clhs324 =             N[3]*clhs23;
const double clhs325 =             clhs12*clhs322 + clhs18*clhs323 + clhs18*clhs324;
const double clhs326 =             DN(0,0)*N[3];
const double clhs327 =             clhs109*clhs42;
const double clhs328 =             clhs322*clhs54 + clhs323*clhs55 + clhs324*clhs55;
const double clhs329 =             DN(1,0)*N[3];
const double clhs330 =             clhs322*clhs82 + clhs323*clhs83 + clhs324*clhs83;
const double clhs331 =             DN(2,0)*N[3];
const double clhs332 =             pow(DN(3,0), 2);
const double clhs333 =             pow(N[3], 2)*clhs16 + clhs109*clhs322 + clhs110*clhs323 + clhs110*clhs324;
const double clhs334 =             DN(3,0)*clhs11;
const double clhs335 =             DN(3,1)*clhs334;
const double clhs336 =             DN(3,2)*clhs334;
const double clhs337 =             -clhs125 + clhs322*clhs41 + clhs327;
const double clhs338 =             clhs107 + clhs325;
const double clhs339 =             DN(0,1)*N[3];
const double clhs340 =             clhs238 + clhs328;
const double clhs341 =             DN(1,1)*N[3];
const double clhs342 =             clhs292 + clhs330;
const double clhs343 =             DN(2,1)*N[3];
const double clhs344 =             pow(DN(3,1), 2);
const double clhs345 =             DN(3,1)*DN(3,2)*clhs11;
const double clhs346 =             DN(0,2)*N[3];
const double clhs347 =             DN(1,2)*N[3];
const double clhs348 =             DN(2,2)*N[3];
const double clhs349 =             pow(DN(3,2), 2);
const double clhs350 =             N[3] + clhs215;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs25;
            lhs(0,1)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + DN(0,2)*clhs31 + clhs33;
            lhs(0,2)=DN(0,0)*clhs34 + DN(0,1)*clhs36 + DN(0,2)*clhs38 + clhs39;
            lhs(0,3)=DN(0,0)*clhs44;
            lhs(0,4)=DN(0,0)*clhs45 + DN(0,1)*clhs47 + DN(0,2)*clhs49 + clhs53 + clhs56;
            lhs(0,5)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + DN(0,2)*clhs62 + clhs63;
            lhs(0,6)=DN(0,0)*clhs64 + DN(0,1)*clhs66 + DN(0,2)*clhs68 + clhs69;
            lhs(0,7)=-DN(0,0)*clhs70 + DN(1,0)*clhs43 + clhs71*clhs72;
            lhs(0,8)=DN(0,0)*clhs73 + DN(0,1)*clhs75 + DN(0,2)*clhs77 + clhs81 + clhs84;
            lhs(0,9)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs90 + clhs91;
            lhs(0,10)=DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs96 + clhs97;
            lhs(0,11)=-DN(0,0)*clhs98 + DN(2,0)*clhs43 + clhs72*clhs99;
            lhs(0,12)=DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs108 + clhs111;
            lhs(0,13)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs117 + clhs118;
            lhs(0,14)=DN(0,0)*clhs119 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs124;
            lhs(0,15)=-DN(0,0)*clhs125 + DN(3,0)*clhs43 + clhs126*clhs72;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs127 + DN(0,2)*clhs128 + clhs33;
            lhs(1,1)=DN(0,0)*clhs28 + DN(0,1)*clhs129 + DN(0,2)*clhs131 + clhs11*clhs132 + clhs25;
            lhs(1,2)=DN(0,0)*clhs36 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs137;
            lhs(1,3)=DN(0,1)*clhs44;
            lhs(1,4)=DN(0,0)*clhs47 + DN(0,1)*clhs138 + DN(0,2)*clhs139 + clhs140;
            lhs(1,5)=DN(0,0)*clhs59 + DN(0,1)*clhs141 + DN(0,2)*clhs143 + clhs145 + clhs146;
            lhs(1,6)=DN(0,0)*clhs66 + DN(0,1)*clhs147 + DN(0,2)*clhs149 + clhs150;
            lhs(1,7)=-DN(0,1)*clhs70 + DN(1,1)*clhs43 + clhs151*clhs72;
            lhs(1,8)=DN(0,0)*clhs75 + DN(0,1)*clhs152 + DN(0,2)*clhs153 + clhs154;
            lhs(1,9)=DN(0,0)*clhs87 + DN(0,1)*clhs155 + DN(0,2)*clhs157 + clhs159 + clhs160;
            lhs(1,10)=DN(0,0)*clhs94 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs164;
            lhs(1,11)=-DN(0,1)*clhs98 + DN(2,1)*clhs43 + clhs165*clhs72;
            lhs(1,12)=DN(0,0)*clhs102 + DN(0,1)*clhs166 + DN(0,2)*clhs167 + clhs168;
            lhs(1,13)=DN(0,0)*clhs114 + DN(0,1)*clhs169 + DN(0,2)*clhs171 + clhs173 + clhs174;
            lhs(1,14)=DN(0,0)*clhs121 + DN(0,1)*clhs175 + DN(0,2)*clhs177 + clhs178;
            lhs(1,15)=-DN(0,1)*clhs125 + DN(3,1)*clhs43 + clhs179*clhs72;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs128 + DN(0,2)*clhs180 + clhs39;
            lhs(2,1)=DN(0,0)*clhs31 + DN(0,1)*clhs131 + DN(0,2)*clhs181 + clhs137;
            lhs(2,2)=DN(0,0)*clhs38 + DN(0,1)*clhs135 + DN(0,2)*clhs182 + clhs11*clhs183 + clhs25;
            lhs(2,3)=DN(0,2)*clhs44;
            lhs(2,4)=DN(0,0)*clhs49 + DN(0,1)*clhs139 + DN(0,2)*clhs184 + clhs186;
            lhs(2,5)=DN(0,0)*clhs62 + DN(0,1)*clhs143 + DN(0,2)*clhs187 + clhs188;
            lhs(2,6)=DN(0,0)*clhs68 + DN(0,1)*clhs149 + DN(0,2)*clhs189 + clhs146 + clhs191;
            lhs(2,7)=-DN(0,2)*clhs70 + DN(1,2)*clhs43 + clhs192*clhs72;
            lhs(2,8)=DN(0,0)*clhs77 + DN(0,1)*clhs153 + DN(0,2)*clhs193 + clhs194;
            lhs(2,9)=DN(0,0)*clhs90 + DN(0,1)*clhs157 + DN(0,2)*clhs195 + clhs196;
            lhs(2,10)=DN(0,0)*clhs96 + DN(0,1)*clhs163 + DN(0,2)*clhs197 + clhs160 + clhs199;
            lhs(2,11)=-DN(0,2)*clhs98 + DN(2,2)*clhs43 + clhs200*clhs72;
            lhs(2,12)=DN(0,0)*clhs104 + DN(0,1)*clhs167 + DN(0,2)*clhs201 + clhs202;
            lhs(2,13)=DN(0,0)*clhs117 + DN(0,1)*clhs171 + DN(0,2)*clhs203 + clhs204;
            lhs(2,14)=DN(0,0)*clhs123 + DN(0,1)*clhs177 + DN(0,2)*clhs205 + clhs174 + clhs207;
            lhs(2,15)=-DN(0,2)*clhs125 + DN(3,2)*clhs43 + clhs208*clhs72;
            lhs(3,0)=DN(0,0)*clhs210;
            lhs(3,1)=DN(0,1)*clhs210;
            lhs(3,2)=DN(0,2)*clhs210;
            lhs(3,3)=clhs19*(clhs132 + clhs183 + clhs5);
            lhs(3,4)=DN(0,0)*clhs211 + clhs71;
            lhs(3,5)=DN(0,1)*clhs211 + clhs151;
            lhs(3,6)=DN(0,2)*clhs211 + clhs192;
            lhs(3,7)=clhs212;
            lhs(3,8)=DN(0,0)*clhs213 + clhs99;
            lhs(3,9)=DN(0,1)*clhs213 + clhs165;
            lhs(3,10)=DN(0,2)*clhs213 + clhs200;
            lhs(3,11)=clhs214;
            lhs(3,12)=DN(0,0)*clhs215 + clhs126;
            lhs(3,13)=DN(0,1)*clhs215 + clhs179;
            lhs(3,14)=DN(0,2)*clhs215 + clhs208;
            lhs(3,15)=clhs216;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs220 + clhs53;
            lhs(4,1)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + DN(1,2)*clhs31 + clhs140;
            lhs(4,2)=DN(1,0)*clhs34 + DN(1,1)*clhs36 + DN(1,2)*clhs38 + clhs186;
            lhs(4,3)=DN(0,0)*clhs222 - DN(1,0)*clhs40 + clhs221*clhs72;
            lhs(4,4)=DN(1,0)*clhs45 + DN(1,1)*clhs47 + DN(1,2)*clhs49 + clhs11*clhs223 + clhs224;
            lhs(4,5)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + DN(1,2)*clhs62 + clhs226;
            lhs(4,6)=DN(1,0)*clhs64 + DN(1,1)*clhs66 + DN(1,2)*clhs68 + clhs227;
            lhs(4,7)=DN(1,0)*clhs228;
            lhs(4,8)=DN(1,0)*clhs73 + DN(1,1)*clhs75 + DN(1,2)*clhs77 + clhs232 + clhs233;
            lhs(4,9)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs90 + clhs234;
            lhs(4,10)=DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs96 + clhs235;
            lhs(4,11)=-DN(1,0)*clhs98 + DN(2,0)*clhs222 + clhs236*clhs72;
            lhs(4,12)=DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs239 + clhs240;
            lhs(4,13)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs117 + clhs241;
            lhs(4,14)=DN(1,0)*clhs119 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs242;
            lhs(4,15)=-DN(1,0)*clhs125 + DN(3,0)*clhs222 + clhs243*clhs72;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs127 + DN(1,2)*clhs128 + clhs63;
            lhs(5,1)=DN(1,0)*clhs28 + DN(1,1)*clhs129 + DN(1,2)*clhs131 + clhs145 + clhs244;
            lhs(5,2)=DN(1,0)*clhs36 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs188;
            lhs(5,3)=DN(0,1)*clhs222 - DN(1,1)*clhs40 + clhs245*clhs72;
            lhs(5,4)=DN(1,0)*clhs47 + DN(1,1)*clhs138 + DN(1,2)*clhs139 + clhs226;
            lhs(5,5)=DN(1,0)*clhs59 + DN(1,1)*clhs141 + DN(1,2)*clhs143 + clhs11*clhs246 + clhs224;
            lhs(5,6)=DN(1,0)*clhs66 + DN(1,1)*clhs147 + DN(1,2)*clhs149 + clhs248;
            lhs(5,7)=DN(1,1)*clhs228;
            lhs(5,8)=DN(1,0)*clhs75 + DN(1,1)*clhs152 + DN(1,2)*clhs153 + clhs249;
            lhs(5,9)=DN(1,0)*clhs87 + DN(1,1)*clhs155 + DN(1,2)*clhs157 + clhs251 + clhs252;
            lhs(5,10)=DN(1,0)*clhs94 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs253;
            lhs(5,11)=-DN(1,1)*clhs98 + DN(2,1)*clhs222 + clhs254*clhs72;
            lhs(5,12)=DN(1,0)*clhs102 + DN(1,1)*clhs166 + DN(1,2)*clhs167 + clhs255;
            lhs(5,13)=DN(1,0)*clhs114 + DN(1,1)*clhs169 + DN(1,2)*clhs171 + clhs257 + clhs258;
            lhs(5,14)=DN(1,0)*clhs121 + DN(1,1)*clhs175 + DN(1,2)*clhs177 + clhs259;
            lhs(5,15)=-DN(1,1)*clhs125 + DN(3,1)*clhs222 + clhs260*clhs72;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs128 + DN(1,2)*clhs180 + clhs69;
            lhs(6,1)=DN(1,0)*clhs31 + DN(1,1)*clhs131 + DN(1,2)*clhs181 + clhs150;
            lhs(6,2)=DN(1,0)*clhs38 + DN(1,1)*clhs135 + DN(1,2)*clhs182 + clhs191 + clhs244;
            lhs(6,3)=DN(0,2)*clhs222 - DN(1,2)*clhs40 + clhs261*clhs72;
            lhs(6,4)=DN(1,0)*clhs49 + DN(1,1)*clhs139 + DN(1,2)*clhs184 + clhs227;
            lhs(6,5)=DN(1,0)*clhs62 + DN(1,1)*clhs143 + DN(1,2)*clhs187 + clhs248;
            lhs(6,6)=DN(1,0)*clhs68 + DN(1,1)*clhs149 + DN(1,2)*clhs189 + clhs11*clhs262 + clhs224;
            lhs(6,7)=DN(1,2)*clhs228;
            lhs(6,8)=DN(1,0)*clhs77 + DN(1,1)*clhs153 + DN(1,2)*clhs193 + clhs264;
            lhs(6,9)=DN(1,0)*clhs90 + DN(1,1)*clhs157 + DN(1,2)*clhs195 + clhs265;
            lhs(6,10)=DN(1,0)*clhs96 + DN(1,1)*clhs163 + DN(1,2)*clhs197 + clhs252 + clhs267;
            lhs(6,11)=-DN(1,2)*clhs98 + DN(2,2)*clhs222 + clhs268*clhs72;
            lhs(6,12)=DN(1,0)*clhs104 + DN(1,1)*clhs167 + DN(1,2)*clhs201 + clhs269;
            lhs(6,13)=DN(1,0)*clhs117 + DN(1,1)*clhs171 + DN(1,2)*clhs203 + clhs270;
            lhs(6,14)=DN(1,0)*clhs123 + DN(1,1)*clhs177 + DN(1,2)*clhs205 + clhs258 + clhs272;
            lhs(6,15)=-DN(1,2)*clhs125 + DN(3,2)*clhs222 + clhs273*clhs72;
            lhs(7,0)=DN(1,0)*clhs209 + clhs221;
            lhs(7,1)=DN(1,1)*clhs209 + clhs245;
            lhs(7,2)=DN(1,2)*clhs209 + clhs261;
            lhs(7,3)=clhs212;
            lhs(7,4)=DN(1,0)*clhs274;
            lhs(7,5)=DN(1,1)*clhs274;
            lhs(7,6)=DN(1,2)*clhs274;
            lhs(7,7)=clhs19*(clhs223 + clhs246 + clhs262);
            lhs(7,8)=DN(1,0)*clhs213 + clhs236;
            lhs(7,9)=DN(1,1)*clhs213 + clhs254;
            lhs(7,10)=DN(1,2)*clhs213 + clhs268;
            lhs(7,11)=clhs275;
            lhs(7,12)=DN(1,0)*clhs215 + clhs243;
            lhs(7,13)=DN(1,1)*clhs215 + clhs260;
            lhs(7,14)=DN(1,2)*clhs215 + clhs273;
            lhs(7,15)=clhs276;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs280 + clhs81;
            lhs(8,1)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + DN(2,2)*clhs31 + clhs154;
            lhs(8,2)=DN(2,0)*clhs34 + DN(2,1)*clhs36 + DN(2,2)*clhs38 + clhs194;
            lhs(8,3)=DN(0,0)*clhs282 - DN(2,0)*clhs40 + clhs281*clhs72;
            lhs(8,4)=DN(2,0)*clhs45 + DN(2,1)*clhs47 + DN(2,2)*clhs49 + clhs232 + clhs283;
            lhs(8,5)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + DN(2,2)*clhs62 + clhs249;
            lhs(8,6)=DN(2,0)*clhs64 + DN(2,1)*clhs66 + DN(2,2)*clhs68 + clhs264;
            lhs(8,7)=DN(1,0)*clhs282 - DN(2,0)*clhs70 + clhs284*clhs72;
            lhs(8,8)=DN(2,0)*clhs73 + DN(2,1)*clhs75 + DN(2,2)*clhs77 + clhs11*clhs285 + clhs286;
            lhs(8,9)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs90 + clhs288;
            lhs(8,10)=DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs96 + clhs289;
            lhs(8,11)=DN(2,0)*clhs290;
            lhs(8,12)=DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs293 + clhs294;
            lhs(8,13)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs117 + clhs295;
            lhs(8,14)=DN(2,0)*clhs119 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs296;
            lhs(8,15)=-DN(2,0)*clhs125 + DN(3,0)*clhs282 + clhs297*clhs72;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs127 + DN(2,2)*clhs128 + clhs91;
            lhs(9,1)=DN(2,0)*clhs28 + DN(2,1)*clhs129 + DN(2,2)*clhs131 + clhs159 + clhs298;
            lhs(9,2)=DN(2,0)*clhs36 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs196;
            lhs(9,3)=DN(0,1)*clhs282 - DN(2,1)*clhs40 + clhs299*clhs72;
            lhs(9,4)=DN(2,0)*clhs47 + DN(2,1)*clhs138 + DN(2,2)*clhs139 + clhs234;
            lhs(9,5)=DN(2,0)*clhs59 + DN(2,1)*clhs141 + DN(2,2)*clhs143 + clhs251 + clhs300;
            lhs(9,6)=DN(2,0)*clhs66 + DN(2,1)*clhs147 + DN(2,2)*clhs149 + clhs265;
            lhs(9,7)=DN(1,1)*clhs282 - DN(2,1)*clhs70 + clhs301*clhs72;
            lhs(9,8)=DN(2,0)*clhs75 + DN(2,1)*clhs152 + DN(2,2)*clhs153 + clhs288;
            lhs(9,9)=DN(2,0)*clhs87 + DN(2,1)*clhs155 + DN(2,2)*clhs157 + clhs11*clhs302 + clhs286;
            lhs(9,10)=DN(2,0)*clhs94 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs304;
            lhs(9,11)=DN(2,1)*clhs290;
            lhs(9,12)=DN(2,0)*clhs102 + DN(2,1)*clhs166 + DN(2,2)*clhs167 + clhs305;
            lhs(9,13)=DN(2,0)*clhs114 + DN(2,1)*clhs169 + DN(2,2)*clhs171 + clhs307 + clhs308;
            lhs(9,14)=DN(2,0)*clhs121 + DN(2,1)*clhs175 + DN(2,2)*clhs177 + clhs309;
            lhs(9,15)=-DN(2,1)*clhs125 + DN(3,1)*clhs282 + clhs310*clhs72;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs128 + DN(2,2)*clhs180 + clhs97;
            lhs(10,1)=DN(2,0)*clhs31 + DN(2,1)*clhs131 + DN(2,2)*clhs181 + clhs164;
            lhs(10,2)=DN(2,0)*clhs38 + DN(2,1)*clhs135 + DN(2,2)*clhs182 + clhs199 + clhs298;
            lhs(10,3)=DN(0,2)*clhs282 - DN(2,2)*clhs40 + clhs311*clhs72;
            lhs(10,4)=DN(2,0)*clhs49 + DN(2,1)*clhs139 + DN(2,2)*clhs184 + clhs235;
            lhs(10,5)=DN(2,0)*clhs62 + DN(2,1)*clhs143 + DN(2,2)*clhs187 + clhs253;
            lhs(10,6)=DN(2,0)*clhs68 + DN(2,1)*clhs149 + DN(2,2)*clhs189 + clhs267 + clhs300;
            lhs(10,7)=DN(1,2)*clhs282 - DN(2,2)*clhs70 + clhs312*clhs72;
            lhs(10,8)=DN(2,0)*clhs77 + DN(2,1)*clhs153 + DN(2,2)*clhs193 + clhs289;
            lhs(10,9)=DN(2,0)*clhs90 + DN(2,1)*clhs157 + DN(2,2)*clhs195 + clhs304;
            lhs(10,10)=DN(2,0)*clhs96 + DN(2,1)*clhs163 + DN(2,2)*clhs197 + clhs11*clhs313 + clhs286;
            lhs(10,11)=DN(2,2)*clhs290;
            lhs(10,12)=DN(2,0)*clhs104 + DN(2,1)*clhs167 + DN(2,2)*clhs201 + clhs315;
            lhs(10,13)=DN(2,0)*clhs117 + DN(2,1)*clhs171 + DN(2,2)*clhs203 + clhs316;
            lhs(10,14)=DN(2,0)*clhs123 + DN(2,1)*clhs177 + DN(2,2)*clhs205 + clhs308 + clhs318;
            lhs(10,15)=-DN(2,2)*clhs125 + DN(3,2)*clhs282 + clhs319*clhs72;
            lhs(11,0)=DN(2,0)*clhs209 + clhs281;
            lhs(11,1)=DN(2,1)*clhs209 + clhs299;
            lhs(11,2)=DN(2,2)*clhs209 + clhs311;
            lhs(11,3)=clhs214;
            lhs(11,4)=DN(2,0)*clhs211 + clhs284;
            lhs(11,5)=DN(2,1)*clhs211 + clhs301;
            lhs(11,6)=DN(2,2)*clhs211 + clhs312;
            lhs(11,7)=clhs275;
            lhs(11,8)=DN(2,0)*clhs320;
            lhs(11,9)=DN(2,1)*clhs320;
            lhs(11,10)=DN(2,2)*clhs320;
            lhs(11,11)=clhs19*(clhs285 + clhs302 + clhs313);
            lhs(11,12)=DN(2,0)*clhs215 + clhs297;
            lhs(11,13)=DN(2,1)*clhs215 + clhs310;
            lhs(11,14)=DN(2,2)*clhs215 + clhs319;
            lhs(11,15)=clhs321;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs108 + clhs325;
            lhs(12,1)=DN(3,0)*clhs26 + DN(3,1)*clhs28 + DN(3,2)*clhs31 + clhs168;
            lhs(12,2)=DN(3,0)*clhs34 + DN(3,1)*clhs36 + DN(3,2)*clhs38 + clhs202;
            lhs(12,3)=DN(0,0)*clhs327 - DN(3,0)*clhs40 + clhs326*clhs72;
            lhs(12,4)=DN(3,0)*clhs45 + DN(3,1)*clhs47 + DN(3,2)*clhs49 + clhs239 + clhs328;
            lhs(12,5)=DN(3,0)*clhs57 + DN(3,1)*clhs59 + DN(3,2)*clhs62 + clhs255;
            lhs(12,6)=DN(3,0)*clhs64 + DN(3,1)*clhs66 + DN(3,2)*clhs68 + clhs269;
            lhs(12,7)=DN(1,0)*clhs327 - DN(3,0)*clhs70 + clhs329*clhs72;
            lhs(12,8)=DN(3,0)*clhs73 + DN(3,1)*clhs75 + DN(3,2)*clhs77 + clhs293 + clhs330;
            lhs(12,9)=DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs90 + clhs305;
            lhs(12,10)=DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs96 + clhs315;
            lhs(12,11)=DN(2,0)*clhs327 - DN(3,0)*clhs98 + clhs331*clhs72;
            lhs(12,12)=DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs11*clhs332 + clhs333;
            lhs(12,13)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs117 + clhs335;
            lhs(12,14)=DN(3,0)*clhs119 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs336;
            lhs(12,15)=DN(3,0)*clhs337;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs127 + DN(3,2)*clhs128 + clhs118;
            lhs(13,1)=DN(3,0)*clhs28 + DN(3,1)*clhs129 + DN(3,2)*clhs131 + clhs173 + clhs338;
            lhs(13,2)=DN(3,0)*clhs36 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs204;
            lhs(13,3)=DN(0,1)*clhs327 - DN(3,1)*clhs40 + clhs339*clhs72;
            lhs(13,4)=DN(3,0)*clhs47 + DN(3,1)*clhs138 + DN(3,2)*clhs139 + clhs241;
            lhs(13,5)=DN(3,0)*clhs59 + DN(3,1)*clhs141 + DN(3,2)*clhs143 + clhs257 + clhs340;
            lhs(13,6)=DN(3,0)*clhs66 + DN(3,1)*clhs147 + DN(3,2)*clhs149 + clhs270;
            lhs(13,7)=DN(1,1)*clhs327 - DN(3,1)*clhs70 + clhs341*clhs72;
            lhs(13,8)=DN(3,0)*clhs75 + DN(3,1)*clhs152 + DN(3,2)*clhs153 + clhs295;
            lhs(13,9)=DN(3,0)*clhs87 + DN(3,1)*clhs155 + DN(3,2)*clhs157 + clhs307 + clhs342;
            lhs(13,10)=DN(3,0)*clhs94 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs316;
            lhs(13,11)=DN(2,1)*clhs327 - DN(3,1)*clhs98 + clhs343*clhs72;
            lhs(13,12)=DN(3,0)*clhs102 + DN(3,1)*clhs166 + DN(3,2)*clhs167 + clhs335;
            lhs(13,13)=DN(3,0)*clhs114 + DN(3,1)*clhs169 + DN(3,2)*clhs171 + clhs11*clhs344 + clhs333;
            lhs(13,14)=DN(3,0)*clhs121 + DN(3,1)*clhs175 + DN(3,2)*clhs177 + clhs345;
            lhs(13,15)=DN(3,1)*clhs337;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs128 + DN(3,2)*clhs180 + clhs124;
            lhs(14,1)=DN(3,0)*clhs31 + DN(3,1)*clhs131 + DN(3,2)*clhs181 + clhs178;
            lhs(14,2)=DN(3,0)*clhs38 + DN(3,1)*clhs135 + DN(3,2)*clhs182 + clhs207 + clhs338;
            lhs(14,3)=DN(0,2)*clhs327 - DN(3,2)*clhs40 + clhs346*clhs72;
            lhs(14,4)=DN(3,0)*clhs49 + DN(3,1)*clhs139 + DN(3,2)*clhs184 + clhs242;
            lhs(14,5)=DN(3,0)*clhs62 + DN(3,1)*clhs143 + DN(3,2)*clhs187 + clhs259;
            lhs(14,6)=DN(3,0)*clhs68 + DN(3,1)*clhs149 + DN(3,2)*clhs189 + clhs272 + clhs340;
            lhs(14,7)=DN(1,2)*clhs327 - DN(3,2)*clhs70 + clhs347*clhs72;
            lhs(14,8)=DN(3,0)*clhs77 + DN(3,1)*clhs153 + DN(3,2)*clhs193 + clhs296;
            lhs(14,9)=DN(3,0)*clhs90 + DN(3,1)*clhs157 + DN(3,2)*clhs195 + clhs309;
            lhs(14,10)=DN(3,0)*clhs96 + DN(3,1)*clhs163 + DN(3,2)*clhs197 + clhs318 + clhs342;
            lhs(14,11)=DN(2,2)*clhs327 - DN(3,2)*clhs98 + clhs348*clhs72;
            lhs(14,12)=DN(3,0)*clhs104 + DN(3,1)*clhs167 + DN(3,2)*clhs201 + clhs336;
            lhs(14,13)=DN(3,0)*clhs117 + DN(3,1)*clhs171 + DN(3,2)*clhs203 + clhs345;
            lhs(14,14)=DN(3,0)*clhs123 + DN(3,1)*clhs177 + DN(3,2)*clhs205 + clhs11*clhs349 + clhs333;
            lhs(14,15)=DN(3,2)*clhs337;
            lhs(15,0)=DN(3,0)*clhs209 + clhs326;
            lhs(15,1)=DN(3,1)*clhs209 + clhs339;
            lhs(15,2)=DN(3,2)*clhs209 + clhs346;
            lhs(15,3)=clhs216;
            lhs(15,4)=DN(3,0)*clhs211 + clhs329;
            lhs(15,5)=DN(3,1)*clhs211 + clhs341;
            lhs(15,6)=DN(3,2)*clhs211 + clhs347;
            lhs(15,7)=clhs276;
            lhs(15,8)=DN(3,0)*clhs213 + clhs331;
            lhs(15,9)=DN(3,1)*clhs213 + clhs343;
            lhs(15,10)=DN(3,2)*clhs213 + clhs348;
            lhs(15,11)=clhs321;
            lhs(15,12)=DN(3,0)*clhs350;
            lhs(15,13)=DN(3,1)*clhs350;
            lhs(15,14)=DN(3,2)*clhs350;
            lhs(15,15)=clhs19*(clhs332 + clhs344 + clhs349);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

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

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error =-rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0));
const double crhs1 =             1.0/dt;
const double crhs2 =             N[0]*rho;
const double crhs3 =             1.0*crhs1*crhs2;
const double crhs4 =             theta - 1.0;
const double crhs5 =             N[0]*(-crhs4*pn[0] + p[0]*theta) + N[1]*(-crhs4*pn[1] + p[1]*theta) + N[2]*(-crhs4*pn[2] + p[2]*theta);
const double crhs6 =             f(0,0)*theta;
const double crhs7 =             f(1,0)*theta;
const double crhs8 =             f(2,0)*theta;
const double crhs9 =             N[0]*(-crhs4*fn(0,0) + crhs6) + N[1]*(-crhs4*fn(1,0) + crhs7) + N[2]*(-crhs4*fn(2,0) + crhs8);
const double crhs10 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs11 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs12 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs13 =             rho*(crhs10*crhs11 + crhs12*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs14 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs15 =             crhs10 + crhs14 - volume_error_ratio;
const double crhs16 =             mu*stab_c1;
const double crhs17 =             rho*sqrt(pow(crhs11, 2) + pow(crhs12, 2));
const double crhs18 =             crhs15*(crhs16/stab_c2 + crhs17*h);
const double crhs19 =             crhs1*rho;
const double crhs20 =             1.0*crhs19;
const double crhs21 =             crhs0*crhs20;
const double crhs22 =             1.0 - theta;
const double crhs23 =             1.0/(crhs16/pow(h, 2) + crhs17*stab_c2/h + crhs19*dyn_tau);
const double crhs24 =             crhs23*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + crhs13 + crhs21 - rho*(N[0]*(crhs22*fn(0,0) + crhs6) + N[1]*(crhs22*fn(1,0) + crhs7) + N[2]*(crhs22*fn(2,0) + crhs8)));
const double crhs25 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs26 =             crhs2*crhs25;
const double crhs27 =             rho*(DN(0,0)*crhs11 + DN(0,1)*crhs12);
const double crhs28 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1));
const double crhs29 =             f(0,1)*theta;
const double crhs30 =             f(1,1)*theta;
const double crhs31 =             f(2,1)*theta;
const double crhs32 =             N[0]*(crhs29 - crhs4*fn(0,1)) + N[1]*(crhs30 - crhs4*fn(1,1)) + N[2]*(crhs31 - crhs4*fn(2,1));
const double crhs33 =             rho*(crhs11*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs12*crhs14);
const double crhs34 =             crhs20*crhs28;
const double crhs35 =             crhs23*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs33 + crhs34 - rho*(N[0]*(crhs22*fn(0,1) + crhs29) + N[1]*(crhs22*fn(1,1) + crhs30) + N[2]*(crhs22*fn(2,1) + crhs31)));
const double crhs36 =             N[1]*rho;
const double crhs37 =             crhs25*crhs36;
const double crhs38 =             rho*(DN(1,0)*crhs11 + DN(1,1)*crhs12);
const double crhs39 =             N[2]*rho;
const double crhs40 =             crhs25*crhs39;
const double crhs41 =             rho*(DN(2,0)*crhs11 + DN(2,1)*crhs12);
            rhs[0]=-DN(0,0)*crhs18 + DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] - N[0]*crhs13 - crhs0*crhs3 + crhs2*crhs9 - crhs24*crhs26 - crhs24*crhs27;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs18 + DN(0,1)*crhs5 - DN(0,1)*stress[1] - N[0]*crhs33 + crhs2*crhs32 - crhs26*crhs35 - crhs27*crhs35 - crhs28*crhs3;
            rhs[2]=-DN(0,0)*crhs24 - DN(0,1)*crhs35 - N[0]*crhs15;
            rhs[3]=-DN(1,0)*crhs18 + DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] - N[1]*crhs13 - N[1]*crhs21 - crhs24*crhs37 - crhs24*crhs38 + crhs36*crhs9;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs18 + DN(1,1)*crhs5 - DN(1,1)*stress[1] - N[1]*crhs33 - N[1]*crhs34 + crhs32*crhs36 - crhs35*crhs37 - crhs35*crhs38;
            rhs[5]=-DN(1,0)*crhs24 - DN(1,1)*crhs35 - N[1]*crhs15;
            rhs[6]=-DN(2,0)*crhs18 + DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] - N[2]*crhs13 - N[2]*crhs21 - crhs24*crhs40 - crhs24*crhs41 + crhs39*crhs9;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs18 + DN(2,1)*crhs5 - DN(2,1)*stress[1] - N[2]*crhs33 - N[2]*crhs34 + crhs32*crhs39 - crhs35*crhs40 - crhs35*crhs41;
            rhs[8]=-DN(2,0)*crhs24 - DN(2,1)*crhs35 - N[2]*crhs15;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    VectorType &rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

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

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)) + N[3]*(v(3,0) - vn(3,0));
const double crhs1 =             1.0/dt;
const double crhs2 =             N[0]*rho;
const double crhs3 =             1.0*crhs1*crhs2;
const double crhs4 =             theta - 1.0;
const double crhs5 =             N[0]*(-crhs4*pn[0] + p[0]*theta) + N[1]*(-crhs4*pn[1] + p[1]*theta) + N[2]*(-crhs4*pn[2] + p[2]*theta) + N[3]*(-crhs4*pn[3] + p[3]*theta);
const double crhs6 =             f(0,0)*theta;
const double crhs7 =             f(1,0)*theta;
const double crhs8 =             f(2,0)*theta;
const double crhs9 =             f(3,0)*theta;
const double crhs10 =             N[0]*(-crhs4*fn(0,0) + crhs6) + N[1]*(-crhs4*fn(1,0) + crhs7) + N[2]*(-crhs4*fn(2,0) + crhs8) + N[3]*(-crhs4*fn(3,0) + crhs9);
const double crhs11 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs12 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs13 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs14 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs15 =             rho*(crhs11*crhs12 + crhs13*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs14*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs16 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs17 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs18 =             crhs11 + crhs16 + crhs17 - volume_error_ratio;
const double crhs19 =             mu*stab_c1;
const double crhs20 =             rho*sqrt(pow(crhs12, 2) + pow(crhs13, 2) + pow(crhs14, 2));
const double crhs21 =             crhs18*(crhs19/stab_c2 + crhs20*h);
const double crhs22 =             crhs1*rho;
const double crhs23 =             1.0*crhs22;
const double crhs24 =             crhs0*crhs23;
const double crhs25 =             1.0 - theta;
const double crhs26 =             1.0/(crhs19/pow(h, 2) + crhs20*stab_c2/h + crhs22*dyn_tau);
const double crhs27 =             crhs26*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crhs15 + crhs24 - rho*(N[0]*(crhs25*fn(0,0) + crhs6) + N[1]*(crhs25*fn(1,0) + crhs7) + N[2]*(crhs25*fn(2,0) + crhs8) + N[3]*(crhs25*fn(3,0) + crhs9)));
const double crhs28 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs29 =             crhs2*crhs28;
const double crhs30 =             rho*(DN(0,0)*crhs12 + DN(0,1)*crhs13 + DN(0,2)*crhs14);
const double crhs31 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)) + N[3]*(v(3,1) - vn(3,1));
const double crhs32 =             f(0,1)*theta;
const double crhs33 =             f(1,1)*theta;
const double crhs34 =             f(2,1)*theta;
const double crhs35 =             f(3,1)*theta;
const double crhs36 =             N[0]*(crhs32 - crhs4*fn(0,1)) + N[1]*(crhs33 - crhs4*fn(1,1)) + N[2]*(crhs34 - crhs4*fn(2,1)) + N[3]*(crhs35 - crhs4*fn(3,1));
const double crhs37 =             rho*(crhs12*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs13*crhs16 + crhs14*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs38 =             crhs23*crhs31;
const double crhs39 =             crhs26*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs37 + crhs38 - rho*(N[0]*(crhs25*fn(0,1) + crhs32) + N[1]*(crhs25*fn(1,1) + crhs33) + N[2]*(crhs25*fn(2,1) + crhs34) + N[3]*(crhs25*fn(3,1) + crhs35)));
const double crhs40 =             N[0]*(v(0,2) - vn(0,2)) + N[1]*(v(1,2) - vn(1,2)) + N[2]*(v(2,2) - vn(2,2)) + N[3]*(v(3,2) - vn(3,2));
const double crhs41 =             f(0,2)*theta;
const double crhs42 =             f(1,2)*theta;
const double crhs43 =             f(2,2)*theta;
const double crhs44 =             f(3,2)*theta;
const double crhs45 =             N[0]*(-crhs4*fn(0,2) + crhs41) + N[1]*(-crhs4*fn(1,2) + crhs42) + N[2]*(-crhs4*fn(2,2) + crhs43) + N[3]*(-crhs4*fn(3,2) + crhs44);
const double crhs46 =             rho*(crhs12*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs13*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs14*crhs17);
const double crhs47 =             crhs23*crhs40;
const double crhs48 =             crhs26*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + crhs46 + crhs47 - rho*(N[0]*(crhs25*fn(0,2) + crhs41) + N[1]*(crhs25*fn(1,2) + crhs42) + N[2]*(crhs25*fn(2,2) + crhs43) + N[3]*(crhs25*fn(3,2) + crhs44)));
const double crhs49 =             N[1]*rho;
const double crhs50 =             crhs28*crhs49;
const double crhs51 =             rho*(DN(1,0)*crhs12 + DN(1,1)*crhs13 + DN(1,2)*crhs14);
const double crhs52 =             N[2]*rho;
const double crhs53 =             crhs28*crhs52;
const double crhs54 =             rho*(DN(2,0)*crhs12 + DN(2,1)*crhs13 + DN(2,2)*crhs14);
const double crhs55 =             N[3]*rho;
const double crhs56 =             crhs28*crhs55;
const double crhs57 =             rho*(DN(3,0)*crhs12 + DN(3,1)*crhs13 + DN(3,2)*crhs14);
            rhs[0]=-DN(0,0)*crhs21 + DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] - N[0]*crhs15 - crhs0*crhs3 + crhs10*crhs2 - crhs27*crhs29 - crhs27*crhs30;
            rhs[1]=-DN(0,0)*stress[3] - DN(0,1)*crhs21 + DN(0,1)*crhs5 - DN(0,1)*stress[1] - DN(0,2)*stress[4] - N[0]*crhs37 + crhs2*crhs36 - crhs29*crhs39 - crhs3*crhs31 - crhs30*crhs39;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] - DN(0,2)*crhs21 + DN(0,2)*crhs5 - DN(0,2)*stress[2] - N[0]*crhs46 + crhs2*crhs45 - crhs29*crhs48 - crhs3*crhs40 - crhs30*crhs48;
            rhs[3]=-DN(0,0)*crhs27 - DN(0,1)*crhs39 - DN(0,2)*crhs48 - N[0]*crhs18;
            rhs[4]=-DN(1,0)*crhs21 + DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] - N[1]*crhs15 - N[1]*crhs24 + crhs10*crhs49 - crhs27*crhs50 - crhs27*crhs51;
            rhs[5]=-DN(1,0)*stress[3] - DN(1,1)*crhs21 + DN(1,1)*crhs5 - DN(1,1)*stress[1] - DN(1,2)*stress[4] - N[1]*crhs37 - N[1]*crhs38 + crhs36*crhs49 - crhs39*crhs50 - crhs39*crhs51;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] - DN(1,2)*crhs21 + DN(1,2)*crhs5 - DN(1,2)*stress[2] - N[1]*crhs46 - N[1]*crhs47 + crhs45*crhs49 - crhs48*crhs50 - crhs48*crhs51;
            rhs[7]=-DN(1,0)*crhs27 - DN(1,1)*crhs39 - DN(1,2)*crhs48 - N[1]*crhs18;
            rhs[8]=-DN(2,0)*crhs21 + DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] - N[2]*crhs15 - N[2]*crhs24 + crhs10*crhs52 - crhs27*crhs53 - crhs27*crhs54;
            rhs[9]=-DN(2,0)*stress[3] - DN(2,1)*crhs21 + DN(2,1)*crhs5 - DN(2,1)*stress[1] - DN(2,2)*stress[4] - N[2]*crhs37 - N[2]*crhs38 + crhs36*crhs52 - crhs39*crhs53 - crhs39*crhs54;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] - DN(2,2)*crhs21 + DN(2,2)*crhs5 - DN(2,2)*stress[2] - N[2]*crhs46 - N[2]*crhs47 + crhs45*crhs52 - crhs48*crhs53 - crhs48*crhs54;
            rhs[11]=-DN(2,0)*crhs27 - DN(2,1)*crhs39 - DN(2,2)*crhs48 - N[2]*crhs18;
            rhs[12]=-DN(3,0)*crhs21 + DN(3,0)*crhs5 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] - N[3]*crhs15 - N[3]*crhs24 + crhs10*crhs55 - crhs27*crhs56 - crhs27*crhs57;
            rhs[13]=-DN(3,0)*stress[3] - DN(3,1)*crhs21 + DN(3,1)*crhs5 - DN(3,1)*stress[1] - DN(3,2)*stress[4] - N[3]*crhs37 - N[3]*crhs38 + crhs36*crhs55 - crhs39*crhs56 - crhs39*crhs57;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] - DN(3,2)*crhs21 + DN(3,2)*crhs5 - DN(3,2)*stress[2] - N[3]*crhs46 - N[3]*crhs47 + crhs45*crhs55 - crhs48*crhs56 - crhs48*crhs57;
            rhs[15]=-DN(3,0)*crhs27 - DN(3,1)*crhs39 - DN(3,2)*crhs48 - N[3]*crhs18;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesCNData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocityOldStep;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p=rData.Pressure;

    const BoundedMatrix<double,3,2> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error = -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

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


    const double cH0 =             1.0/dt;
const double cH1 =             1.0*cH0;
const double cH2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH4 =             1.0/(cH0*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2));
const double cH5 =             cH4*rho;
const double cH6 =             cH5*(DN(0,0)*cH2 + DN(0,1)*cH3 + N[0]*cH1);
const double cH7 =             cH5*(DN(1,0)*cH2 + DN(1,1)*cH3 + N[1]*cH1);
const double cH8 =             cH5*(DN(2,0)*cH2 + DN(2,1)*cH3 + N[2]*cH1);
            H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH6;
            H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH6;
            H(0,2)=cH4*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH7;
            H(0,4)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH7;
            H(0,5)=cH4*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH8;
            H(0,7)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH8;
            H(0,8)=cH4*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH6;
            H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH6;
            H(1,2)=cH4*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH7;
            H(1,4)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH7;
            H(1,5)=cH4*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH8;
            H(1,7)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH8;
            H(1,8)=cH4*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH6;
            H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH6;
            H(2,2)=cH4*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH7;
            H(2,4)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH7;
            H(2,5)=cH4*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH8;
            H(2,7)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH8;
            H(2,8)=cH4*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


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
const double crhs_ee2 =             crhs_ee0 + crhs_ee1 - volume_error_ratio;
const double crhs_ee3 =             1.0 - theta;
const double crhs_ee4 =             1.0/dt;
const double crhs_ee5 =             1.0*crhs_ee4;
const double crhs_ee6 =             N[0]*crhs_ee5;
const double crhs_ee7 =             N[1]*crhs_ee5;
const double crhs_ee8 =             N[2]*crhs_ee5;
const double crhs_ee9 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee10 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee11 =             1.0/(crhs_ee4*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee12 =             crhs_ee11*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*(crhs_ee3*fn(0,0) + f(0,0)*theta) + N[1]*(crhs_ee3*fn(1,0) + f(1,0)*theta) + N[2]*(crhs_ee3*fn(2,0) + f(2,0)*theta)) + rho*(crhs_ee0*crhs_ee9 + crhs_ee10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)) + crhs_ee6*(v(0,0) - vn(0,0)) + crhs_ee7*(v(1,0) - vn(1,0)) + crhs_ee8*(v(2,0) - vn(2,0))));
const double crhs_ee13 =             crhs_ee11*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*(crhs_ee3*fn(0,1) + f(0,1)*theta) + N[1]*(crhs_ee3*fn(1,1) + f(1,1)*theta) + N[2]*(crhs_ee3*fn(2,1) + f(2,1)*theta)) + rho*(crhs_ee1*crhs_ee10 + crhs_ee6*(v(0,1) - vn(0,1)) + crhs_ee7*(v(1,1) - vn(1,1)) + crhs_ee8*(v(2,1) - vn(2,1)) + crhs_ee9*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee12 - DNenr(0,1)*crhs_ee13 - Nenr[0]*crhs_ee2;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee12 - DNenr(1,1)*crhs_ee13 - Nenr[1]*crhs_ee2;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee12 - DNenr(2,1)*crhs_ee13 - Nenr[2]*crhs_ee2;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesCNData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double theta = rData.theta;

    const double dyn_tau = rData.DynamicTau;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vmeshn = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &fn = rData.BodyForce_OldStep1;
    const auto &p = rData.Pressure;

    const BoundedMatrix<double,4,3> vconv = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Mass correction term
    double volume_error_ratio = 0.0;
    if (rData.IsCut()) {
        const double volume_error =-rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

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


    const double cH0 =             1.0/dt;
const double cH1 =             1.0*cH0;
const double cH2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH5 =             1.0/(cH0*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH2, 2) + pow(cH3, 2) + pow(cH4, 2))/h + mu*stab_c1/pow(h, 2));
const double cH6 =             cH5*rho;
const double cH7 =             cH6*(DN(0,0)*cH2 + DN(0,1)*cH3 + DN(0,2)*cH4 + N[0]*cH1);
const double cH8 =             cH6*(DN(1,0)*cH2 + DN(1,1)*cH3 + DN(1,2)*cH4 + N[1]*cH1);
const double cH9 =             cH6*(DN(2,0)*cH2 + DN(2,1)*cH3 + DN(2,2)*cH4 + N[2]*cH1);
const double cH10 =             cH6*(DN(3,0)*cH2 + DN(3,1)*cH3 + DN(3,2)*cH4 + N[3]*cH1);
            H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH7;
            H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH7;
            H(0,2)=DN(0,2)*Nenr[0] + DNenr(0,2)*cH7;
            H(0,3)=cH5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH8;
            H(0,5)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH8;
            H(0,6)=DN(1,2)*Nenr[0] + DNenr(0,2)*cH8;
            H(0,7)=cH5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH9;
            H(0,9)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH9;
            H(0,10)=DN(2,2)*Nenr[0] + DNenr(0,2)*cH9;
            H(0,11)=cH5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DN(3,0)*Nenr[0] + DNenr(0,0)*cH10;
            H(0,13)=DN(3,1)*Nenr[0] + DNenr(0,1)*cH10;
            H(0,14)=DN(3,2)*Nenr[0] + DNenr(0,2)*cH10;
            H(0,15)=cH5*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH7;
            H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH7;
            H(1,2)=DN(0,2)*Nenr[1] + DNenr(1,2)*cH7;
            H(1,3)=cH5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH8;
            H(1,5)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH8;
            H(1,6)=DN(1,2)*Nenr[1] + DNenr(1,2)*cH8;
            H(1,7)=cH5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH9;
            H(1,9)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH9;
            H(1,10)=DN(2,2)*Nenr[1] + DNenr(1,2)*cH9;
            H(1,11)=cH5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DN(3,0)*Nenr[1] + DNenr(1,0)*cH10;
            H(1,13)=DN(3,1)*Nenr[1] + DNenr(1,1)*cH10;
            H(1,14)=DN(3,2)*Nenr[1] + DNenr(1,2)*cH10;
            H(1,15)=cH5*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH7;
            H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH7;
            H(2,2)=DN(0,2)*Nenr[2] + DNenr(2,2)*cH7;
            H(2,3)=cH5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH8;
            H(2,5)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH8;
            H(2,6)=DN(1,2)*Nenr[2] + DNenr(2,2)*cH8;
            H(2,7)=cH5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH9;
            H(2,9)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH9;
            H(2,10)=DN(2,2)*Nenr[2] + DNenr(2,2)*cH9;
            H(2,11)=cH5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DN(3,0)*Nenr[2] + DNenr(2,0)*cH10;
            H(2,13)=DN(3,1)*Nenr[2] + DNenr(2,1)*cH10;
            H(2,14)=DN(3,2)*Nenr[2] + DNenr(2,2)*cH10;
            H(2,15)=cH5*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DN(0,0)*Nenr[3] + DNenr(3,0)*cH7;
            H(3,1)=DN(0,1)*Nenr[3] + DNenr(3,1)*cH7;
            H(3,2)=DN(0,2)*Nenr[3] + DNenr(3,2)*cH7;
            H(3,3)=cH5*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DN(1,0)*Nenr[3] + DNenr(3,0)*cH8;
            H(3,5)=DN(1,1)*Nenr[3] + DNenr(3,1)*cH8;
            H(3,6)=DN(1,2)*Nenr[3] + DNenr(3,2)*cH8;
            H(3,7)=cH5*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DN(2,0)*Nenr[3] + DNenr(3,0)*cH9;
            H(3,9)=DN(2,1)*Nenr[3] + DNenr(3,1)*cH9;
            H(3,10)=DN(2,2)*Nenr[3] + DNenr(3,2)*cH9;
            H(3,11)=cH5*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DN(3,0)*Nenr[3] + DNenr(3,0)*cH10;
            H(3,13)=DN(3,1)*Nenr[3] + DNenr(3,1)*cH10;
            H(3,14)=DN(3,2)*Nenr[3] + DNenr(3,2)*cH10;
            H(3,15)=cH5*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


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


    const double crhs_ee0 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs_ee1 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs_ee2 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee3 =             crhs_ee0 + crhs_ee1 + crhs_ee2 - volume_error_ratio;
const double crhs_ee4 =             1.0 - theta;
const double crhs_ee5 =             1.0/dt;
const double crhs_ee6 =             1.0*crhs_ee5;
const double crhs_ee7 =             N[0]*crhs_ee6;
const double crhs_ee8 =             N[1]*crhs_ee6;
const double crhs_ee9 =             N[2]*crhs_ee6;
const double crhs_ee10 =             N[3]*crhs_ee6;
const double crhs_ee11 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee12 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee13 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee14 =             1.0/(crhs_ee5*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee11, 2) + pow(crhs_ee12, 2) + pow(crhs_ee13, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee15 =             crhs_ee14*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*(crhs_ee4*fn(0,0) + f(0,0)*theta) + N[1]*(crhs_ee4*fn(1,0) + f(1,0)*theta) + N[2]*(crhs_ee4*fn(2,0) + f(2,0)*theta) + N[3]*(crhs_ee4*fn(3,0) + f(3,0)*theta)) + rho*(crhs_ee0*crhs_ee11 + crhs_ee10*(v(3,0) - vn(3,0)) + crhs_ee12*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee13*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs_ee7*(v(0,0) - vn(0,0)) + crhs_ee8*(v(1,0) - vn(1,0)) + crhs_ee9*(v(2,0) - vn(2,0))));
const double crhs_ee16 =             crhs_ee14*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*(crhs_ee4*fn(0,1) + f(0,1)*theta) + N[1]*(crhs_ee4*fn(1,1) + f(1,1)*theta) + N[2]*(crhs_ee4*fn(2,1) + f(2,1)*theta) + N[3]*(crhs_ee4*fn(3,1) + f(3,1)*theta)) + rho*(crhs_ee1*crhs_ee12 + crhs_ee10*(v(3,1) - vn(3,1)) + crhs_ee11*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee13*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs_ee7*(v(0,1) - vn(0,1)) + crhs_ee8*(v(1,1) - vn(1,1)) + crhs_ee9*(v(2,1) - vn(2,1))));
const double crhs_ee17 =             crhs_ee14*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*(crhs_ee4*fn(0,2) + f(0,2)*theta) + N[1]*(crhs_ee4*fn(1,2) + f(1,2)*theta) + N[2]*(crhs_ee4*fn(2,2) + f(2,2)*theta) + N[3]*(crhs_ee4*fn(3,2) + f(3,2)*theta)) + rho*(crhs_ee10*(v(3,2) - vn(3,2)) + crhs_ee11*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee12*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs_ee13*crhs_ee2 + crhs_ee7*(v(0,2) - vn(0,2)) + crhs_ee8*(v(1,2) - vn(1,2)) + crhs_ee9*(v(2,2) - vn(2,2))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee15 - DNenr(0,1)*crhs_ee16 - DNenr(0,2)*crhs_ee17 - Nenr[0]*crhs_ee3;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee15 - DNenr(1,1)*crhs_ee16 - DNenr(1,2)*crhs_ee17 - Nenr[1]*crhs_ee3;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee15 - DNenr(2,1)*crhs_ee16 - DNenr(2,2)*crhs_ee17 - Nenr[2]*crhs_ee3;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee15 - DNenr(3,1)*crhs_ee16 - DNenr(3,2)*crhs_ee17 - Nenr[3]*crhs_ee3;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::ComputeSplitting(
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
void TwoFluidNavierStokesCN<TElementData>::ComputeSplitInterface(
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
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesCNData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
   return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesCNData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateCurvatureOnInterfaceGaussPoints(
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
void TwoFluidNavierStokesCN<TElementData>::SurfaceTension(
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
void TwoFluidNavierStokesCN<TElementData>::PressureGradientStabilization(
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

    const double theta = rData.theta;
    const auto vmesh=rData.MeshVelocity;
    const auto vmeshn=rData.MeshVelocityOldStep;
    const BoundedMatrix<double,NumNodes,Dim> v_convection = theta*(v-vmesh) + (1-theta)*(vn-vmeshn);

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
void TwoFluidNavierStokesCN<TElementData>::CondenseEnrichmentWithContinuity(
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
                    //FIXME: ESTO QUE ES???
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
void TwoFluidNavierStokesCN<TElementData>::CondenseEnrichment(
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
void TwoFluidNavierStokesCN<TElementData>::AddSurfaceTensionContribution(
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
void TwoFluidNavierStokesCN<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateOnIntegrationPoints(
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

template class TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<2, 3>>;
template class TwoFluidNavierStokesCN<TwoFluidNavierStokesCNData<3, 4>>;

} // namespace Kratos
