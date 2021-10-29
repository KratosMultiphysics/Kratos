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
    const double d_gauss = inner_prod(rData.Distance, rN);
    if (d_gauss > 0.0)
        rData.CalculateAirMaterialResponse();
    else
        this->CalculateMaterialResponse(rData);
    rData.ComputeDarcyTerm();
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
    rData.ComputeDarcyTerm();
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

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const double theta=rData.theta;
    const BoundedMatrix<double,3,2> v_CN = theta*v+(1-theta)*vn;
    const auto &vconv_CN = v_CN - vmesh;



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
const double clhs1 =             DN(0,0)*theta;
const double clhs2 =             C(0,2)*DN(0,0);
const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
const double clhs4 =             DN(0,1)*theta;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double clhs7 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double clhs8 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2));
const double clhs9 =             clhs8*h/stab_c1 + mu;
const double clhs10 =             clhs9*theta;
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             K_darcy*theta;
const double clhs13 =             rho/dt;
const double clhs14 =             N[0]*theta;
const double clhs15 =             DN(0,0)*clhs6 + DN(0,1)*clhs7;
const double clhs16 =             clhs15*rho;
const double clhs17 =             K_darcy*clhs14;
const double clhs18 =             N[0]*clhs13;
const double clhs19 =             rho*theta;
const double clhs20 =             clhs15*clhs19;
const double clhs21 =             clhs17 + clhs18 + clhs20;
const double clhs22 =             1.0/(K_darcy + clhs13*dyn_tau + clhs8/h + mu*stab_c1/pow(h, 2));
const double clhs23 =             1.0*clhs22;
const double clhs24 =             clhs16*clhs23;
const double clhs25 =             K_darcy*clhs23;
const double clhs26 =             N[0]*clhs25;
const double clhs27 =             clhs11*clhs12 + clhs11*clhs13 + clhs14*clhs16 + clhs21*clhs24 - clhs21*clhs26;
const double clhs28 =             C(0,1)*DN(0,1) + clhs2;
const double clhs29 =             C(1,2)*DN(0,1);
const double clhs30 =             C(2,2)*DN(0,0) + clhs29;
const double clhs31 =             DN(0,0)*clhs9;
const double clhs32 =             DN(0,1)*clhs31;
const double clhs33 =             -N[0] + clhs24 - clhs26;
const double clhs34 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs35 =             C(0,2)*DN(1,0);
const double clhs36 =             C(2,2)*DN(1,1) + clhs35;
const double clhs37 =             DN(0,0)*DN(1,0);
const double clhs38 =             N[1]*clhs17 + N[1]*clhs18;
const double clhs39 =             clhs10*clhs37 + clhs38;
const double clhs40 =             DN(1,0)*clhs6 + DN(1,1)*clhs7;
const double clhs41 =             clhs14*rho;
const double clhs42 =             N[1]*theta;
const double clhs43 =             K_darcy*clhs42;
const double clhs44 =             N[1]*clhs13;
const double clhs45 =             clhs19*clhs40;
const double clhs46 =             clhs43 + clhs44 + clhs45;
const double clhs47 =             clhs24*clhs46 - clhs26*clhs46 + clhs40*clhs41;
const double clhs48 =             C(0,1)*DN(1,1) + clhs35;
const double clhs49 =             C(1,2)*DN(1,1);
const double clhs50 =             C(2,2)*DN(1,0) + clhs49;
const double clhs51 =             DN(1,1)*clhs31;
const double clhs52 =             DN(0,0)*N[1];
const double clhs53 =             DN(1,0)*N[0];
const double clhs54 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs55 =             C(0,2)*DN(2,0);
const double clhs56 =             C(2,2)*DN(2,1) + clhs55;
const double clhs57 =             DN(0,0)*DN(2,0);
const double clhs58 =             N[2]*clhs17 + N[2]*clhs18;
const double clhs59 =             clhs10*clhs57 + clhs58;
const double clhs60 =             DN(2,0)*clhs6 + DN(2,1)*clhs7;
const double clhs61 =             N[2]*theta;
const double clhs62 =             K_darcy*clhs61;
const double clhs63 =             N[2]*clhs13;
const double clhs64 =             clhs19*clhs60;
const double clhs65 =             clhs62 + clhs63 + clhs64;
const double clhs66 =             clhs24*clhs65 - clhs26*clhs65 + clhs41*clhs60;
const double clhs67 =             C(0,1)*DN(2,1) + clhs55;
const double clhs68 =             C(1,2)*DN(2,1);
const double clhs69 =             C(2,2)*DN(2,0) + clhs68;
const double clhs70 =             DN(2,1)*clhs31;
const double clhs71 =             DN(0,0)*N[2];
const double clhs72 =             DN(2,0)*N[0];
const double clhs73 =             C(0,1)*DN(0,0) + clhs29;
const double clhs74 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs75 =             pow(DN(0,1), 2);
const double clhs76 =             C(0,1)*DN(1,0) + clhs49;
const double clhs77 =             DN(0,1)*clhs9;
const double clhs78 =             DN(1,0)*clhs77;
const double clhs79 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs80 =             DN(0,1)*DN(1,1);
const double clhs81 =             clhs10*clhs80 + clhs38;
const double clhs82 =             DN(0,1)*N[1];
const double clhs83 =             DN(1,1)*N[0];
const double clhs84 =             C(0,1)*DN(2,0) + clhs68;
const double clhs85 =             DN(2,0)*clhs77;
const double clhs86 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs87 =             DN(0,1)*DN(2,1);
const double clhs88 =             clhs10*clhs87 + clhs58;
const double clhs89 =             DN(0,1)*N[2];
const double clhs90 =             DN(2,1)*N[0];
const double clhs91 =             clhs14 + clhs22*(1.0*clhs17 + 1.0*clhs18 + 1.0*clhs20);
const double clhs92 =             clhs23*theta;
const double clhs93 =             clhs23*clhs46;
const double clhs94 =             clhs92*(clhs37 + clhs80);
const double clhs95 =             clhs23*clhs65;
const double clhs96 =             clhs92*(clhs57 + clhs87);
const double clhs97 =             DN(1,0)*theta;
const double clhs98 =             DN(1,1)*theta;
const double clhs99 =             clhs40*rho;
const double clhs100 =             clhs23*clhs99;
const double clhs101 =             N[1]*clhs25;
const double clhs102 =             clhs100*clhs21 - clhs101*clhs21 + clhs16*clhs42;
const double clhs103 =             pow(DN(1,0), 2);
const double clhs104 =             pow(N[1], 2);
const double clhs105 =             -clhs101*clhs46 + clhs104*clhs12 + clhs104*clhs13 + clhs42*clhs99 + clhs93*clhs99;
const double clhs106 =             DN(1,0)*clhs9;
const double clhs107 =             DN(1,1)*clhs106;
const double clhs108 =             -N[1] + clhs100 - clhs101;
const double clhs109 =             DN(1,0)*DN(2,0);
const double clhs110 =             N[2]*clhs43 + N[2]*clhs44;
const double clhs111 =             clhs10*clhs109 + clhs110;
const double clhs112 =             clhs60*rho;
const double clhs113 =             -clhs101*clhs65 + clhs112*clhs42 + clhs95*clhs99;
const double clhs114 =             DN(2,1)*clhs106;
const double clhs115 =             DN(1,0)*N[2];
const double clhs116 =             DN(2,0)*N[1];
const double clhs117 =             pow(DN(1,1), 2);
const double clhs118 =             DN(2,0)*clhs9;
const double clhs119 =             DN(1,1)*clhs118;
const double clhs120 =             DN(1,1)*DN(2,1);
const double clhs121 =             clhs10*clhs120 + clhs110;
const double clhs122 =             DN(1,1)*N[2];
const double clhs123 =             DN(2,1)*N[1];
const double clhs124 =             clhs21*clhs23;
const double clhs125 =             clhs22*(1.0*clhs43 + 1.0*clhs44 + 1.0*clhs45) + clhs42;
const double clhs126 =             clhs92*(clhs109 + clhs120);
const double clhs127 =             DN(2,0)*theta;
const double clhs128 =             DN(2,1)*theta;
const double clhs129 =             N[2]*clhs25;
const double clhs130 =             clhs112*clhs124 - clhs129*clhs21 + clhs16*clhs61;
const double clhs131 =             clhs112*clhs23;
const double clhs132 =             clhs112*clhs93 - clhs129*clhs46 + clhs61*clhs99;
const double clhs133 =             pow(DN(2,0), 2);
const double clhs134 =             pow(N[2], 2);
const double clhs135 =             clhs112*clhs61 + clhs112*clhs95 + clhs12*clhs134 - clhs129*clhs65 + clhs13*clhs134;
const double clhs136 =             DN(2,1)*clhs118;
const double clhs137 =             -N[2] - clhs129 + clhs131;
const double clhs138 =             pow(DN(2,1), 2);
const double clhs139 =             clhs22*(1.0*clhs62 + 1.0*clhs63 + 1.0*clhs64) + clhs61;
            lhs(0,0)=clhs0*clhs1 + clhs10*clhs5 + clhs27 + clhs3*clhs4;
            lhs(0,1)=theta*(DN(0,0)*clhs28 + DN(0,1)*clhs30 + clhs32);
            lhs(0,2)=clhs1*clhs33;
            lhs(0,3)=clhs1*clhs34 + clhs36*clhs4 + clhs39 + clhs47;
            lhs(0,4)=theta*(DN(0,0)*clhs48 + DN(0,1)*clhs50 + clhs51);
            lhs(0,5)=theta*(DN(1,0)*clhs24 - clhs25*clhs53 - clhs52);
            lhs(0,6)=clhs1*clhs54 + clhs4*clhs56 + clhs59 + clhs66;
            lhs(0,7)=theta*(DN(0,0)*clhs67 + DN(0,1)*clhs69 + clhs70);
            lhs(0,8)=theta*(DN(2,0)*clhs24 - clhs25*clhs72 - clhs71);
            lhs(1,0)=theta*(DN(0,0)*clhs3 + DN(0,1)*clhs73 + clhs32);
            lhs(1,1)=clhs1*clhs30 + clhs10*clhs75 + clhs27 + clhs4*clhs74;
            lhs(1,2)=clhs33*clhs4;
            lhs(1,3)=theta*(DN(0,0)*clhs36 + DN(0,1)*clhs76 + clhs78);
            lhs(1,4)=clhs1*clhs50 + clhs4*clhs79 + clhs47 + clhs81;
            lhs(1,5)=theta*(DN(1,1)*clhs24 - clhs25*clhs83 - clhs82);
            lhs(1,6)=theta*(DN(0,0)*clhs56 + DN(0,1)*clhs84 + clhs85);
            lhs(1,7)=clhs1*clhs69 + clhs4*clhs86 + clhs66 + clhs88;
            lhs(1,8)=theta*(DN(2,1)*clhs24 - clhs25*clhs90 - clhs89);
            lhs(2,0)=DN(0,0)*clhs91;
            lhs(2,1)=DN(0,1)*clhs91;
            lhs(2,2)=clhs92*(clhs5 + clhs75);
            lhs(2,3)=DN(0,0)*clhs93 + DN(1,0)*clhs14;
            lhs(2,4)=DN(0,1)*clhs93 + DN(1,1)*clhs14;
            lhs(2,5)=clhs94;
            lhs(2,6)=DN(0,0)*clhs95 + DN(2,0)*clhs14;
            lhs(2,7)=DN(0,1)*clhs95 + DN(2,1)*clhs14;
            lhs(2,8)=clhs96;
            lhs(3,0)=clhs0*clhs97 + clhs102 + clhs3*clhs98 + clhs39;
            lhs(3,1)=theta*(DN(1,0)*clhs28 + DN(1,1)*clhs30 + clhs78);
            lhs(3,2)=theta*(DN(0,0)*clhs100 - clhs25*clhs52 - clhs53);
            lhs(3,3)=clhs10*clhs103 + clhs105 + clhs34*clhs97 + clhs36*clhs98;
            lhs(3,4)=theta*(DN(1,0)*clhs48 + DN(1,1)*clhs50 + clhs107);
            lhs(3,5)=clhs108*clhs97;
            lhs(3,6)=clhs111 + clhs113 + clhs54*clhs97 + clhs56*clhs98;
            lhs(3,7)=theta*(DN(1,0)*clhs67 + DN(1,1)*clhs69 + clhs114);
            lhs(3,8)=theta*(DN(2,0)*clhs100 - clhs115 - clhs116*clhs25);
            lhs(4,0)=theta*(DN(1,0)*clhs3 + DN(1,1)*clhs73 + clhs51);
            lhs(4,1)=clhs102 + clhs30*clhs97 + clhs74*clhs98 + clhs81;
            lhs(4,2)=theta*(DN(0,1)*clhs100 - clhs25*clhs82 - clhs83);
            lhs(4,3)=theta*(DN(1,0)*clhs36 + DN(1,1)*clhs76 + clhs107);
            lhs(4,4)=clhs10*clhs117 + clhs105 + clhs50*clhs97 + clhs79*clhs98;
            lhs(4,5)=clhs108*clhs98;
            lhs(4,6)=theta*(DN(1,0)*clhs56 + DN(1,1)*clhs84 + clhs119);
            lhs(4,7)=clhs113 + clhs121 + clhs69*clhs97 + clhs86*clhs98;
            lhs(4,8)=theta*(DN(2,1)*clhs100 - clhs122 - clhs123*clhs25);
            lhs(5,0)=DN(1,0)*clhs124 + clhs52*theta;
            lhs(5,1)=DN(1,1)*clhs124 + clhs82*theta;
            lhs(5,2)=clhs94;
            lhs(5,3)=DN(1,0)*clhs125;
            lhs(5,4)=DN(1,1)*clhs125;
            lhs(5,5)=clhs92*(clhs103 + clhs117);
            lhs(5,6)=DN(1,0)*clhs95 + DN(2,0)*clhs42;
            lhs(5,7)=DN(1,1)*clhs95 + DN(2,1)*clhs42;
            lhs(5,8)=clhs126;
            lhs(6,0)=clhs0*clhs127 + clhs128*clhs3 + clhs130 + clhs59;
            lhs(6,1)=theta*(DN(2,0)*clhs28 + DN(2,1)*clhs30 + clhs85);
            lhs(6,2)=theta*(DN(0,0)*clhs131 - clhs25*clhs71 - clhs72);
            lhs(6,3)=clhs111 + clhs127*clhs34 + clhs128*clhs36 + clhs132;
            lhs(6,4)=theta*(DN(2,0)*clhs48 + DN(2,1)*clhs50 + clhs119);
            lhs(6,5)=theta*(DN(1,0)*clhs131 - clhs115*clhs25 - clhs116);
            lhs(6,6)=clhs10*clhs133 + clhs127*clhs54 + clhs128*clhs56 + clhs135;
            lhs(6,7)=theta*(DN(2,0)*clhs67 + DN(2,1)*clhs69 + clhs136);
            lhs(6,8)=clhs127*clhs137;
            lhs(7,0)=theta*(DN(2,0)*clhs3 + DN(2,1)*clhs73 + clhs70);
            lhs(7,1)=clhs127*clhs30 + clhs128*clhs74 + clhs130 + clhs88;
            lhs(7,2)=theta*(DN(0,1)*clhs131 - clhs25*clhs89 - clhs90);
            lhs(7,3)=theta*(DN(2,0)*clhs36 + DN(2,1)*clhs76 + clhs114);
            lhs(7,4)=clhs121 + clhs127*clhs50 + clhs128*clhs79 + clhs132;
            lhs(7,5)=theta*(DN(1,1)*clhs131 - clhs122*clhs25 - clhs123);
            lhs(7,6)=theta*(DN(2,0)*clhs56 + DN(2,1)*clhs84 + clhs136);
            lhs(7,7)=clhs10*clhs138 + clhs127*clhs69 + clhs128*clhs86 + clhs135;
            lhs(7,8)=clhs128*clhs137;
            lhs(8,0)=DN(2,0)*clhs124 + clhs71*theta;
            lhs(8,1)=DN(2,1)*clhs124 + clhs89*theta;
            lhs(8,2)=clhs96;
            lhs(8,3)=DN(2,0)*clhs93 + clhs115*theta;
            lhs(8,4)=DN(2,1)*clhs93 + clhs122*theta;
            lhs(8,5)=clhs126;
            lhs(8,6)=DN(2,0)*clhs139;
            lhs(8,7)=DN(2,1)*clhs139;
            lhs(8,8)=clhs92*(clhs133 + clhs138);


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
    const double K_darcy = rData.DarcyTerm;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const double dyn_tau = rData.DynamicTau;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &theta=rData.theta;
    const BoundedMatrix<double,4,3> v_CN = theta*v+(1-theta)*vn;
    const auto &vconv_CN = v_CN - vmesh;

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
const double clhs1 =             DN(0,0)*theta;
const double clhs2 =             C(0,3)*DN(0,0);
const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
const double clhs4 =             DN(0,1)*theta;
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             DN(0,2)*theta;
const double clhs8 =             pow(DN(0,0), 2);
const double clhs9 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double clhs10 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double clhs11 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double clhs12 =             rho*stab_c2*sqrt(pow(clhs10, 2) + pow(clhs11, 2) + pow(clhs9, 2));
const double clhs13 =             clhs12*h/stab_c1 + mu;
const double clhs14 =             clhs13*theta;
const double clhs15 =             pow(N[0], 2);
const double clhs16 =             K_darcy*theta;
const double clhs17 =             rho/dt;
const double clhs18 =             N[0]*theta;
const double clhs19 =             DN(0,0)*clhs9 + DN(0,1)*clhs10 + DN(0,2)*clhs11;
const double clhs20 =             clhs19*rho;
const double clhs21 =             K_darcy*clhs18;
const double clhs22 =             N[0]*clhs17;
const double clhs23 =             rho*theta;
const double clhs24 =             clhs19*clhs23;
const double clhs25 =             clhs21 + clhs22 + clhs24;
const double clhs26 =             1.0/(K_darcy + clhs12/h + clhs17*dyn_tau + mu*stab_c1/pow(h, 2));
const double clhs27 =             1.0*clhs26;
const double clhs28 =             clhs20*clhs27;
const double clhs29 =             K_darcy*clhs27;
const double clhs30 =             N[0]*clhs29;
const double clhs31 =             clhs15*clhs16 + clhs15*clhs17 + clhs18*clhs20 + clhs25*clhs28 - clhs25*clhs30;
const double clhs32 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
const double clhs33 =             C(1,3)*DN(0,1);
const double clhs34 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs33;
const double clhs35 =             C(3,5)*DN(0,0);
const double clhs36 =             C(4,5)*DN(0,2);
const double clhs37 =             C(1,5)*DN(0,1) + clhs35 + clhs36;
const double clhs38 =             DN(0,0)*clhs13;
const double clhs39 =             DN(0,1)*clhs38;
const double clhs40 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs41 =             C(3,4)*DN(0,1);
const double clhs42 =             C(2,3)*DN(0,2) + clhs35 + clhs41;
const double clhs43 =             C(2,5)*DN(0,2);
const double clhs44 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs43;
const double clhs45 =             DN(0,2)*clhs38;
const double clhs46 =             -N[0] + clhs28 - clhs30;
const double clhs47 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs48 =             C(0,3)*DN(1,0);
const double clhs49 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs48;
const double clhs50 =             C(0,5)*DN(1,0);
const double clhs51 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs50;
const double clhs52 =             DN(0,0)*DN(1,0);
const double clhs53 =             N[1]*clhs21 + N[1]*clhs22;
const double clhs54 =             clhs14*clhs52 + clhs53;
const double clhs55 =             DN(1,0)*clhs9 + DN(1,1)*clhs10 + DN(1,2)*clhs11;
const double clhs56 =             clhs18*rho;
const double clhs57 =             N[1]*theta;
const double clhs58 =             K_darcy*clhs57;
const double clhs59 =             N[1]*clhs17;
const double clhs60 =             clhs23*clhs55;
const double clhs61 =             clhs58 + clhs59 + clhs60;
const double clhs62 =             clhs28*clhs61 - clhs30*clhs61 + clhs55*clhs56;
const double clhs63 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs48;
const double clhs64 =             C(1,3)*DN(1,1);
const double clhs65 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs64;
const double clhs66 =             C(3,5)*DN(1,0);
const double clhs67 =             C(4,5)*DN(1,2);
const double clhs68 =             C(1,5)*DN(1,1) + clhs66 + clhs67;
const double clhs69 =             DN(1,1)*clhs38;
const double clhs70 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs50;
const double clhs71 =             C(3,4)*DN(1,1);
const double clhs72 =             C(2,3)*DN(1,2) + clhs66 + clhs71;
const double clhs73 =             C(2,5)*DN(1,2);
const double clhs74 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs73;
const double clhs75 =             DN(1,2)*clhs38;
const double clhs76 =             DN(0,0)*N[1];
const double clhs77 =             DN(1,0)*N[0];
const double clhs78 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs79 =             C(0,3)*DN(2,0);
const double clhs80 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs79;
const double clhs81 =             C(0,5)*DN(2,0);
const double clhs82 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs81;
const double clhs83 =             DN(0,0)*DN(2,0);
const double clhs84 =             N[2]*clhs21 + N[2]*clhs22;
const double clhs85 =             clhs14*clhs83 + clhs84;
const double clhs86 =             DN(2,0)*clhs9 + DN(2,1)*clhs10 + DN(2,2)*clhs11;
const double clhs87 =             N[2]*theta;
const double clhs88 =             K_darcy*clhs87;
const double clhs89 =             N[2]*clhs17;
const double clhs90 =             clhs23*clhs86;
const double clhs91 =             clhs88 + clhs89 + clhs90;
const double clhs92 =             clhs28*clhs91 - clhs30*clhs91 + clhs56*clhs86;
const double clhs93 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs79;
const double clhs94 =             C(1,3)*DN(2,1);
const double clhs95 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs94;
const double clhs96 =             C(3,5)*DN(2,0);
const double clhs97 =             C(4,5)*DN(2,2);
const double clhs98 =             C(1,5)*DN(2,1) + clhs96 + clhs97;
const double clhs99 =             DN(2,1)*clhs38;
const double clhs100 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs81;
const double clhs101 =             C(3,4)*DN(2,1);
const double clhs102 =             C(2,3)*DN(2,2) + clhs101 + clhs96;
const double clhs103 =             C(2,5)*DN(2,2);
const double clhs104 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs103;
const double clhs105 =             DN(2,2)*clhs38;
const double clhs106 =             DN(0,0)*N[2];
const double clhs107 =             DN(2,0)*N[0];
const double clhs108 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs109 =             C(0,3)*DN(3,0);
const double clhs110 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs109;
const double clhs111 =             C(0,5)*DN(3,0);
const double clhs112 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs111;
const double clhs113 =             DN(0,0)*DN(3,0);
const double clhs114 =             N[3]*clhs21 + N[3]*clhs22;
const double clhs115 =             clhs113*clhs14 + clhs114;
const double clhs116 =             DN(3,0)*clhs9 + DN(3,1)*clhs10 + DN(3,2)*clhs11;
const double clhs117 =             N[3]*theta;
const double clhs118 =             K_darcy*clhs117;
const double clhs119 =             N[3]*clhs17;
const double clhs120 =             clhs116*clhs23;
const double clhs121 =             clhs118 + clhs119 + clhs120;
const double clhs122 =             clhs116*clhs56 + clhs121*clhs28 - clhs121*clhs30;
const double clhs123 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs109;
const double clhs124 =             C(1,3)*DN(3,1);
const double clhs125 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs124;
const double clhs126 =             C(3,5)*DN(3,0);
const double clhs127 =             C(4,5)*DN(3,2);
const double clhs128 =             C(1,5)*DN(3,1) + clhs126 + clhs127;
const double clhs129 =             DN(3,1)*clhs38;
const double clhs130 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs111;
const double clhs131 =             C(3,4)*DN(3,1);
const double clhs132 =             C(2,3)*DN(3,2) + clhs126 + clhs131;
const double clhs133 =             C(2,5)*DN(3,2);
const double clhs134 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs133;
const double clhs135 =             DN(3,2)*clhs38;
const double clhs136 =             DN(0,0)*N[3];
const double clhs137 =             DN(3,0)*N[0];
const double clhs138 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs33;
const double clhs139 =             C(0,4)*DN(0,0) + clhs36 + clhs41;
const double clhs140 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs141 =             C(1,4)*DN(0,1);
const double clhs142 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs141;
const double clhs143 =             pow(DN(0,1), 2);
const double clhs144 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs141;
const double clhs145 =             C(2,4)*DN(0,2);
const double clhs146 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs145;
const double clhs147 =             DN(0,1)*clhs13;
const double clhs148 =             DN(0,2)*clhs147;
const double clhs149 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs64;
const double clhs150 =             C(0,4)*DN(1,0) + clhs67 + clhs71;
const double clhs151 =             DN(1,0)*clhs147;
const double clhs152 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs153 =             C(1,4)*DN(1,1);
const double clhs154 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs153;
const double clhs155 =             DN(0,1)*DN(1,1);
const double clhs156 =             clhs14*clhs155;
const double clhs157 =             clhs53 + clhs62;
const double clhs158 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs153;
const double clhs159 =             C(2,4)*DN(1,2);
const double clhs160 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs159;
const double clhs161 =             DN(1,2)*clhs147;
const double clhs162 =             DN(0,1)*N[1];
const double clhs163 =             DN(1,1)*N[0];
const double clhs164 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs94;
const double clhs165 =             C(0,4)*DN(2,0) + clhs101 + clhs97;
const double clhs166 =             DN(2,0)*clhs147;
const double clhs167 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs168 =             C(1,4)*DN(2,1);
const double clhs169 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs168;
const double clhs170 =             DN(0,1)*DN(2,1);
const double clhs171 =             clhs14*clhs170;
const double clhs172 =             clhs84 + clhs92;
const double clhs173 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs168;
const double clhs174 =             C(2,4)*DN(2,2);
const double clhs175 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs174;
const double clhs176 =             DN(2,2)*clhs147;
const double clhs177 =             DN(0,1)*N[2];
const double clhs178 =             DN(2,1)*N[0];
const double clhs179 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs124;
const double clhs180 =             C(0,4)*DN(3,0) + clhs127 + clhs131;
const double clhs181 =             DN(3,0)*clhs147;
const double clhs182 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs183 =             C(1,4)*DN(3,1);
const double clhs184 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs183;
const double clhs185 =             DN(0,1)*DN(3,1);
const double clhs186 =             clhs14*clhs185;
const double clhs187 =             clhs114 + clhs122;
const double clhs188 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs183;
const double clhs189 =             C(2,4)*DN(3,2);
const double clhs190 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs189;
const double clhs191 =             DN(3,2)*clhs147;
const double clhs192 =             DN(0,1)*N[3];
const double clhs193 =             DN(3,1)*N[0];
const double clhs194 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs43;
const double clhs195 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs145;
const double clhs196 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs197 =             pow(DN(0,2), 2);
const double clhs198 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs73;
const double clhs199 =             DN(0,2)*clhs13;
const double clhs200 =             DN(1,0)*clhs199;
const double clhs201 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs159;
const double clhs202 =             DN(1,1)*clhs199;
const double clhs203 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs204 =             DN(0,2)*DN(1,2);
const double clhs205 =             clhs14*clhs204;
const double clhs206 =             DN(0,2)*N[1];
const double clhs207 =             DN(1,2)*N[0];
const double clhs208 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs103;
const double clhs209 =             DN(2,0)*clhs199;
const double clhs210 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs174;
const double clhs211 =             DN(2,1)*clhs199;
const double clhs212 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs213 =             DN(0,2)*DN(2,2);
const double clhs214 =             clhs14*clhs213;
const double clhs215 =             DN(0,2)*N[2];
const double clhs216 =             DN(2,2)*N[0];
const double clhs217 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs133;
const double clhs218 =             DN(3,0)*clhs199;
const double clhs219 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs189;
const double clhs220 =             DN(3,1)*clhs199;
const double clhs221 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs222 =             DN(0,2)*DN(3,2);
const double clhs223 =             clhs14*clhs222;
const double clhs224 =             DN(0,2)*N[3];
const double clhs225 =             DN(3,2)*N[0];
const double clhs226 =             clhs18 + clhs26*(1.0*clhs21 + 1.0*clhs22 + 1.0*clhs24);
const double clhs227 =             clhs27*theta;
const double clhs228 =             clhs27*clhs61;
const double clhs229 =             clhs227*(clhs155 + clhs204 + clhs52);
const double clhs230 =             clhs27*clhs91;
const double clhs231 =             clhs227*(clhs170 + clhs213 + clhs83);
const double clhs232 =             clhs121*clhs27;
const double clhs233 =             clhs227*(clhs113 + clhs185 + clhs222);
const double clhs234 =             DN(1,0)*theta;
const double clhs235 =             DN(1,1)*theta;
const double clhs236 =             DN(1,2)*theta;
const double clhs237 =             clhs55*rho;
const double clhs238 =             clhs237*clhs27;
const double clhs239 =             N[1]*clhs29;
const double clhs240 =             clhs20*clhs57 + clhs238*clhs25 - clhs239*clhs25;
const double clhs241 =             pow(DN(1,0), 2);
const double clhs242 =             pow(N[1], 2);
const double clhs243 =             clhs16*clhs242 + clhs17*clhs242 + clhs228*clhs237 + clhs237*clhs57 - clhs239*clhs61;
const double clhs244 =             DN(1,0)*clhs13;
const double clhs245 =             DN(1,1)*clhs244;
const double clhs246 =             DN(1,2)*clhs244;
const double clhs247 =             -N[1] + clhs238 - clhs239;
const double clhs248 =             DN(1,0)*DN(2,0);
const double clhs249 =             N[2]*clhs58 + N[2]*clhs59;
const double clhs250 =             clhs14*clhs248 + clhs249;
const double clhs251 =             clhs57*rho;
const double clhs252 =             clhs230*clhs237 - clhs239*clhs91 + clhs251*clhs86;
const double clhs253 =             DN(2,1)*clhs244;
const double clhs254 =             DN(2,2)*clhs244;
const double clhs255 =             DN(1,0)*N[2];
const double clhs256 =             DN(2,0)*N[1];
const double clhs257 =             DN(1,0)*DN(3,0);
const double clhs258 =             N[3]*clhs58 + N[3]*clhs59;
const double clhs259 =             clhs14*clhs257 + clhs258;
const double clhs260 =             clhs116*clhs251 - clhs121*clhs239 + clhs232*clhs237;
const double clhs261 =             DN(3,1)*clhs244;
const double clhs262 =             DN(3,2)*clhs244;
const double clhs263 =             DN(1,0)*N[3];
const double clhs264 =             DN(3,0)*N[1];
const double clhs265 =             clhs240 + clhs53;
const double clhs266 =             pow(DN(1,1), 2);
const double clhs267 =             DN(1,1)*clhs13;
const double clhs268 =             DN(1,2)*clhs267;
const double clhs269 =             DN(2,0)*clhs267;
const double clhs270 =             DN(1,1)*DN(2,1);
const double clhs271 =             clhs14*clhs270;
const double clhs272 =             clhs249 + clhs252;
const double clhs273 =             DN(2,2)*clhs267;
const double clhs274 =             DN(1,1)*N[2];
const double clhs275 =             DN(2,1)*N[1];
const double clhs276 =             DN(3,0)*clhs267;
const double clhs277 =             DN(1,1)*DN(3,1);
const double clhs278 =             clhs14*clhs277;
const double clhs279 =             clhs258 + clhs260;
const double clhs280 =             DN(3,2)*clhs267;
const double clhs281 =             DN(1,1)*N[3];
const double clhs282 =             DN(3,1)*N[1];
const double clhs283 =             pow(DN(1,2), 2);
const double clhs284 =             DN(1,2)*clhs13;
const double clhs285 =             DN(2,0)*clhs284;
const double clhs286 =             DN(2,1)*clhs284;
const double clhs287 =             DN(1,2)*DN(2,2);
const double clhs288 =             clhs14*clhs287;
const double clhs289 =             DN(1,2)*N[2];
const double clhs290 =             DN(2,2)*N[1];
const double clhs291 =             DN(3,0)*clhs284;
const double clhs292 =             DN(3,1)*clhs284;
const double clhs293 =             DN(1,2)*DN(3,2);
const double clhs294 =             clhs14*clhs293;
const double clhs295 =             DN(1,2)*N[3];
const double clhs296 =             DN(3,2)*N[1];
const double clhs297 =             clhs25*clhs27;
const double clhs298 =             clhs26*(1.0*clhs58 + 1.0*clhs59 + 1.0*clhs60) + clhs57;
const double clhs299 =             clhs227*(clhs248 + clhs270 + clhs287);
const double clhs300 =             clhs227*(clhs257 + clhs277 + clhs293);
const double clhs301 =             DN(2,0)*theta;
const double clhs302 =             DN(2,1)*theta;
const double clhs303 =             DN(2,2)*theta;
const double clhs304 =             clhs86*rho;
const double clhs305 =             N[2]*clhs29;
const double clhs306 =             clhs20*clhs87 - clhs25*clhs305 + clhs297*clhs304;
const double clhs307 =             clhs27*clhs304;
const double clhs308 =             clhs228*clhs304 + clhs237*clhs87 - clhs305*clhs61;
const double clhs309 =             pow(DN(2,0), 2);
const double clhs310 =             pow(N[2], 2);
const double clhs311 =             clhs16*clhs310 + clhs17*clhs310 + clhs230*clhs304 + clhs304*clhs87 - clhs305*clhs91;
const double clhs312 =             DN(2,0)*clhs13;
const double clhs313 =             DN(2,1)*clhs312;
const double clhs314 =             DN(2,2)*clhs312;
const double clhs315 =             -N[2] - clhs305 + clhs307;
const double clhs316 =             DN(2,0)*DN(3,0);
const double clhs317 =             N[3]*clhs88 + N[3]*clhs89;
const double clhs318 =             clhs14*clhs316 + clhs317;
const double clhs319 =             clhs116*rho;
const double clhs320 =             -clhs121*clhs305 + clhs232*clhs304 + clhs319*clhs87;
const double clhs321 =             DN(3,1)*clhs312;
const double clhs322 =             DN(3,2)*clhs312;
const double clhs323 =             DN(2,0)*N[3];
const double clhs324 =             DN(3,0)*N[2];
const double clhs325 =             clhs306 + clhs84;
const double clhs326 =             clhs249 + clhs308;
const double clhs327 =             pow(DN(2,1), 2);
const double clhs328 =             DN(2,1)*clhs13;
const double clhs329 =             DN(2,2)*clhs328;
const double clhs330 =             DN(3,0)*clhs328;
const double clhs331 =             DN(2,1)*DN(3,1);
const double clhs332 =             clhs14*clhs331;
const double clhs333 =             clhs317 + clhs320;
const double clhs334 =             DN(3,2)*clhs328;
const double clhs335 =             DN(2,1)*N[3];
const double clhs336 =             DN(3,1)*N[2];
const double clhs337 =             pow(DN(2,2), 2);
const double clhs338 =             DN(2,2)*clhs13;
const double clhs339 =             DN(3,0)*clhs338;
const double clhs340 =             DN(3,1)*clhs338;
const double clhs341 =             DN(2,2)*DN(3,2);
const double clhs342 =             clhs14*clhs341;
const double clhs343 =             DN(2,2)*N[3];
const double clhs344 =             DN(3,2)*N[2];
const double clhs345 =             clhs26*(1.0*clhs88 + 1.0*clhs89 + 1.0*clhs90) + clhs87;
const double clhs346 =             clhs227*(clhs316 + clhs331 + clhs341);
const double clhs347 =             DN(3,0)*theta;
const double clhs348 =             DN(3,1)*theta;
const double clhs349 =             DN(3,2)*theta;
const double clhs350 =             N[3]*clhs29;
const double clhs351 =             clhs117*clhs20 - clhs25*clhs350 + clhs297*clhs319;
const double clhs352 =             clhs27*clhs319;
const double clhs353 =             clhs117*clhs237 + clhs228*clhs319 - clhs350*clhs61;
const double clhs354 =             clhs117*clhs304 + clhs230*clhs319 - clhs350*clhs91;
const double clhs355 =             pow(DN(3,0), 2);
const double clhs356 =             pow(N[3], 2);
const double clhs357 =             clhs117*clhs319 - clhs121*clhs350 + clhs16*clhs356 + clhs17*clhs356 + clhs232*clhs319;
const double clhs358 =             DN(3,0)*clhs13;
const double clhs359 =             DN(3,1)*clhs358;
const double clhs360 =             DN(3,2)*clhs358;
const double clhs361 =             -N[3] - clhs350 + clhs352;
const double clhs362 =             clhs114 + clhs351;
const double clhs363 =             clhs258 + clhs353;
const double clhs364 =             clhs317 + clhs354;
const double clhs365 =             pow(DN(3,1), 2);
const double clhs366 =             DN(3,1)*DN(3,2)*clhs13;
const double clhs367 =             pow(DN(3,2), 2);
const double clhs368 =             clhs117 + clhs26*(1.0*clhs118 + 1.0*clhs119 + 1.0*clhs120);
            lhs(0,0)=clhs0*clhs1 + clhs14*clhs8 + clhs3*clhs4 + clhs31 + clhs6*clhs7;
            lhs(0,1)=theta*(DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs37 + clhs39);
            lhs(0,2)=theta*(DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs44 + clhs45);
            lhs(0,3)=clhs1*clhs46;
            lhs(0,4)=clhs1*clhs47 + clhs4*clhs49 + clhs51*clhs7 + clhs54 + clhs62;
            lhs(0,5)=theta*(DN(0,0)*clhs63 + DN(0,1)*clhs65 + DN(0,2)*clhs68 + clhs69);
            lhs(0,6)=theta*(DN(0,0)*clhs70 + DN(0,1)*clhs72 + DN(0,2)*clhs74 + clhs75);
            lhs(0,7)=theta*(DN(1,0)*clhs28 - clhs29*clhs77 - clhs76);
            lhs(0,8)=clhs1*clhs78 + clhs4*clhs80 + clhs7*clhs82 + clhs85 + clhs92;
            lhs(0,9)=theta*(DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs98 + clhs99);
            lhs(0,10)=theta*(DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs105);
            lhs(0,11)=theta*(DN(2,0)*clhs28 - clhs106 - clhs107*clhs29);
            lhs(0,12)=clhs1*clhs108 + clhs110*clhs4 + clhs112*clhs7 + clhs115 + clhs122;
            lhs(0,13)=theta*(DN(0,0)*clhs123 + DN(0,1)*clhs125 + DN(0,2)*clhs128 + clhs129);
            lhs(0,14)=theta*(DN(0,0)*clhs130 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs135);
            lhs(0,15)=theta*(DN(3,0)*clhs28 - clhs136 - clhs137*clhs29);
            lhs(1,0)=theta*(DN(0,0)*clhs3 + DN(0,1)*clhs138 + DN(0,2)*clhs139 + clhs39);
            lhs(1,1)=clhs1*clhs34 + clhs14*clhs143 + clhs140*clhs4 + clhs142*clhs7 + clhs31;
            lhs(1,2)=theta*(DN(0,0)*clhs42 + DN(0,1)*clhs144 + DN(0,2)*clhs146 + clhs148);
            lhs(1,3)=clhs4*clhs46;
            lhs(1,4)=theta*(DN(0,0)*clhs49 + DN(0,1)*clhs149 + DN(0,2)*clhs150 + clhs151);
            lhs(1,5)=clhs1*clhs65 + clhs152*clhs4 + clhs154*clhs7 + clhs156 + clhs157;
            lhs(1,6)=theta*(DN(0,0)*clhs72 + DN(0,1)*clhs158 + DN(0,2)*clhs160 + clhs161);
            lhs(1,7)=theta*(DN(1,1)*clhs28 - clhs162 - clhs163*clhs29);
            lhs(1,8)=theta*(DN(0,0)*clhs80 + DN(0,1)*clhs164 + DN(0,2)*clhs165 + clhs166);
            lhs(1,9)=clhs1*clhs95 + clhs167*clhs4 + clhs169*clhs7 + clhs171 + clhs172;
            lhs(1,10)=theta*(DN(0,0)*clhs102 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176);
            lhs(1,11)=theta*(DN(2,1)*clhs28 - clhs177 - clhs178*clhs29);
            lhs(1,12)=theta*(DN(0,0)*clhs110 + DN(0,1)*clhs179 + DN(0,2)*clhs180 + clhs181);
            lhs(1,13)=clhs1*clhs125 + clhs182*clhs4 + clhs184*clhs7 + clhs186 + clhs187;
            lhs(1,14)=theta*(DN(0,0)*clhs132 + DN(0,1)*clhs188 + DN(0,2)*clhs190 + clhs191);
            lhs(1,15)=theta*(DN(3,1)*clhs28 - clhs192 - clhs193*clhs29);
            lhs(2,0)=theta*(DN(0,0)*clhs6 + DN(0,1)*clhs139 + DN(0,2)*clhs194 + clhs45);
            lhs(2,1)=theta*(DN(0,0)*clhs37 + DN(0,1)*clhs142 + DN(0,2)*clhs195 + clhs148);
            lhs(2,2)=clhs1*clhs44 + clhs14*clhs197 + clhs146*clhs4 + clhs196*clhs7 + clhs31;
            lhs(2,3)=clhs46*clhs7;
            lhs(2,4)=theta*(DN(0,0)*clhs51 + DN(0,1)*clhs150 + DN(0,2)*clhs198 + clhs200);
            lhs(2,5)=theta*(DN(0,0)*clhs68 + DN(0,1)*clhs154 + DN(0,2)*clhs201 + clhs202);
            lhs(2,6)=clhs1*clhs74 + clhs157 + clhs160*clhs4 + clhs203*clhs7 + clhs205;
            lhs(2,7)=theta*(DN(1,2)*clhs28 - clhs206 - clhs207*clhs29);
            lhs(2,8)=theta*(DN(0,0)*clhs82 + DN(0,1)*clhs165 + DN(0,2)*clhs208 + clhs209);
            lhs(2,9)=theta*(DN(0,0)*clhs98 + DN(0,1)*clhs169 + DN(0,2)*clhs210 + clhs211);
            lhs(2,10)=clhs1*clhs104 + clhs172 + clhs175*clhs4 + clhs212*clhs7 + clhs214;
            lhs(2,11)=theta*(DN(2,2)*clhs28 - clhs215 - clhs216*clhs29);
            lhs(2,12)=theta*(DN(0,0)*clhs112 + DN(0,1)*clhs180 + DN(0,2)*clhs217 + clhs218);
            lhs(2,13)=theta*(DN(0,0)*clhs128 + DN(0,1)*clhs184 + DN(0,2)*clhs219 + clhs220);
            lhs(2,14)=clhs1*clhs134 + clhs187 + clhs190*clhs4 + clhs221*clhs7 + clhs223;
            lhs(2,15)=theta*(DN(3,2)*clhs28 - clhs224 - clhs225*clhs29);
            lhs(3,0)=DN(0,0)*clhs226;
            lhs(3,1)=DN(0,1)*clhs226;
            lhs(3,2)=DN(0,2)*clhs226;
            lhs(3,3)=clhs227*(clhs143 + clhs197 + clhs8);
            lhs(3,4)=DN(0,0)*clhs228 + DN(1,0)*clhs18;
            lhs(3,5)=DN(0,1)*clhs228 + DN(1,1)*clhs18;
            lhs(3,6)=DN(0,2)*clhs228 + DN(1,2)*clhs18;
            lhs(3,7)=clhs229;
            lhs(3,8)=DN(0,0)*clhs230 + DN(2,0)*clhs18;
            lhs(3,9)=DN(0,1)*clhs230 + DN(2,1)*clhs18;
            lhs(3,10)=DN(0,2)*clhs230 + DN(2,2)*clhs18;
            lhs(3,11)=clhs231;
            lhs(3,12)=DN(0,0)*clhs232 + DN(3,0)*clhs18;
            lhs(3,13)=DN(0,1)*clhs232 + DN(3,1)*clhs18;
            lhs(3,14)=DN(0,2)*clhs232 + DN(3,2)*clhs18;
            lhs(3,15)=clhs233;
            lhs(4,0)=clhs0*clhs234 + clhs235*clhs3 + clhs236*clhs6 + clhs240 + clhs54;
            lhs(4,1)=theta*(DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs37 + clhs151);
            lhs(4,2)=theta*(DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs44 + clhs200);
            lhs(4,3)=theta*(DN(0,0)*clhs238 - clhs29*clhs76 - clhs77);
            lhs(4,4)=clhs14*clhs241 + clhs234*clhs47 + clhs235*clhs49 + clhs236*clhs51 + clhs243;
            lhs(4,5)=theta*(DN(1,0)*clhs63 + DN(1,1)*clhs65 + DN(1,2)*clhs68 + clhs245);
            lhs(4,6)=theta*(DN(1,0)*clhs70 + DN(1,1)*clhs72 + DN(1,2)*clhs74 + clhs246);
            lhs(4,7)=clhs234*clhs247;
            lhs(4,8)=clhs234*clhs78 + clhs235*clhs80 + clhs236*clhs82 + clhs250 + clhs252;
            lhs(4,9)=theta*(DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs98 + clhs253);
            lhs(4,10)=theta*(DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs254);
            lhs(4,11)=theta*(DN(2,0)*clhs238 - clhs255 - clhs256*clhs29);
            lhs(4,12)=clhs108*clhs234 + clhs110*clhs235 + clhs112*clhs236 + clhs259 + clhs260;
            lhs(4,13)=theta*(DN(1,0)*clhs123 + DN(1,1)*clhs125 + DN(1,2)*clhs128 + clhs261);
            lhs(4,14)=theta*(DN(1,0)*clhs130 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs262);
            lhs(4,15)=theta*(DN(3,0)*clhs238 - clhs263 - clhs264*clhs29);
            lhs(5,0)=theta*(DN(1,0)*clhs3 + DN(1,1)*clhs138 + DN(1,2)*clhs139 + clhs69);
            lhs(5,1)=clhs140*clhs235 + clhs142*clhs236 + clhs156 + clhs234*clhs34 + clhs265;
            lhs(5,2)=theta*(DN(1,0)*clhs42 + DN(1,1)*clhs144 + DN(1,2)*clhs146 + clhs202);
            lhs(5,3)=theta*(DN(0,1)*clhs238 - clhs162*clhs29 - clhs163);
            lhs(5,4)=theta*(DN(1,0)*clhs49 + DN(1,1)*clhs149 + DN(1,2)*clhs150 + clhs245);
            lhs(5,5)=clhs14*clhs266 + clhs152*clhs235 + clhs154*clhs236 + clhs234*clhs65 + clhs243;
            lhs(5,6)=theta*(DN(1,0)*clhs72 + DN(1,1)*clhs158 + DN(1,2)*clhs160 + clhs268);
            lhs(5,7)=clhs235*clhs247;
            lhs(5,8)=theta*(DN(1,0)*clhs80 + DN(1,1)*clhs164 + DN(1,2)*clhs165 + clhs269);
            lhs(5,9)=clhs167*clhs235 + clhs169*clhs236 + clhs234*clhs95 + clhs271 + clhs272;
            lhs(5,10)=theta*(DN(1,0)*clhs102 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs273);
            lhs(5,11)=theta*(DN(2,1)*clhs238 - clhs274 - clhs275*clhs29);
            lhs(5,12)=theta*(DN(1,0)*clhs110 + DN(1,1)*clhs179 + DN(1,2)*clhs180 + clhs276);
            lhs(5,13)=clhs125*clhs234 + clhs182*clhs235 + clhs184*clhs236 + clhs278 + clhs279;
            lhs(5,14)=theta*(DN(1,0)*clhs132 + DN(1,1)*clhs188 + DN(1,2)*clhs190 + clhs280);
            lhs(5,15)=theta*(DN(3,1)*clhs238 - clhs281 - clhs282*clhs29);
            lhs(6,0)=theta*(DN(1,0)*clhs6 + DN(1,1)*clhs139 + DN(1,2)*clhs194 + clhs75);
            lhs(6,1)=theta*(DN(1,0)*clhs37 + DN(1,1)*clhs142 + DN(1,2)*clhs195 + clhs161);
            lhs(6,2)=clhs146*clhs235 + clhs196*clhs236 + clhs205 + clhs234*clhs44 + clhs265;
            lhs(6,3)=theta*(DN(0,2)*clhs238 - clhs206*clhs29 - clhs207);
            lhs(6,4)=theta*(DN(1,0)*clhs51 + DN(1,1)*clhs150 + DN(1,2)*clhs198 + clhs246);
            lhs(6,5)=theta*(DN(1,0)*clhs68 + DN(1,1)*clhs154 + DN(1,2)*clhs201 + clhs268);
            lhs(6,6)=clhs14*clhs283 + clhs160*clhs235 + clhs203*clhs236 + clhs234*clhs74 + clhs243;
            lhs(6,7)=clhs236*clhs247;
            lhs(6,8)=theta*(DN(1,0)*clhs82 + DN(1,1)*clhs165 + DN(1,2)*clhs208 + clhs285);
            lhs(6,9)=theta*(DN(1,0)*clhs98 + DN(1,1)*clhs169 + DN(1,2)*clhs210 + clhs286);
            lhs(6,10)=clhs104*clhs234 + clhs175*clhs235 + clhs212*clhs236 + clhs272 + clhs288;
            lhs(6,11)=theta*(DN(2,2)*clhs238 - clhs289 - clhs29*clhs290);
            lhs(6,12)=theta*(DN(1,0)*clhs112 + DN(1,1)*clhs180 + DN(1,2)*clhs217 + clhs291);
            lhs(6,13)=theta*(DN(1,0)*clhs128 + DN(1,1)*clhs184 + DN(1,2)*clhs219 + clhs292);
            lhs(6,14)=clhs134*clhs234 + clhs190*clhs235 + clhs221*clhs236 + clhs279 + clhs294;
            lhs(6,15)=theta*(DN(3,2)*clhs238 - clhs29*clhs296 - clhs295);
            lhs(7,0)=DN(1,0)*clhs297 + clhs76*theta;
            lhs(7,1)=DN(1,1)*clhs297 + clhs162*theta;
            lhs(7,2)=DN(1,2)*clhs297 + clhs206*theta;
            lhs(7,3)=clhs229;
            lhs(7,4)=DN(1,0)*clhs298;
            lhs(7,5)=DN(1,1)*clhs298;
            lhs(7,6)=DN(1,2)*clhs298;
            lhs(7,7)=clhs227*(clhs241 + clhs266 + clhs283);
            lhs(7,8)=DN(1,0)*clhs230 + DN(2,0)*clhs57;
            lhs(7,9)=DN(1,1)*clhs230 + DN(2,1)*clhs57;
            lhs(7,10)=DN(1,2)*clhs230 + DN(2,2)*clhs57;
            lhs(7,11)=clhs299;
            lhs(7,12)=DN(1,0)*clhs232 + DN(3,0)*clhs57;
            lhs(7,13)=DN(1,1)*clhs232 + DN(3,1)*clhs57;
            lhs(7,14)=DN(1,2)*clhs232 + DN(3,2)*clhs57;
            lhs(7,15)=clhs300;
            lhs(8,0)=clhs0*clhs301 + clhs3*clhs302 + clhs303*clhs6 + clhs306 + clhs85;
            lhs(8,1)=theta*(DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs37 + clhs166);
            lhs(8,2)=theta*(DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs44 + clhs209);
            lhs(8,3)=theta*(DN(0,0)*clhs307 - clhs106*clhs29 - clhs107);
            lhs(8,4)=clhs250 + clhs301*clhs47 + clhs302*clhs49 + clhs303*clhs51 + clhs308;
            lhs(8,5)=theta*(DN(2,0)*clhs63 + DN(2,1)*clhs65 + DN(2,2)*clhs68 + clhs269);
            lhs(8,6)=theta*(DN(2,0)*clhs70 + DN(2,1)*clhs72 + DN(2,2)*clhs74 + clhs285);
            lhs(8,7)=theta*(DN(1,0)*clhs307 - clhs255*clhs29 - clhs256);
            lhs(8,8)=clhs14*clhs309 + clhs301*clhs78 + clhs302*clhs80 + clhs303*clhs82 + clhs311;
            lhs(8,9)=theta*(DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs98 + clhs313);
            lhs(8,10)=theta*(DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs314);
            lhs(8,11)=clhs301*clhs315;
            lhs(8,12)=clhs108*clhs301 + clhs110*clhs302 + clhs112*clhs303 + clhs318 + clhs320;
            lhs(8,13)=theta*(DN(2,0)*clhs123 + DN(2,1)*clhs125 + DN(2,2)*clhs128 + clhs321);
            lhs(8,14)=theta*(DN(2,0)*clhs130 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs322);
            lhs(8,15)=theta*(DN(3,0)*clhs307 - clhs29*clhs324 - clhs323);
            lhs(9,0)=theta*(DN(2,0)*clhs3 + DN(2,1)*clhs138 + DN(2,2)*clhs139 + clhs99);
            lhs(9,1)=clhs140*clhs302 + clhs142*clhs303 + clhs171 + clhs301*clhs34 + clhs325;
            lhs(9,2)=theta*(DN(2,0)*clhs42 + DN(2,1)*clhs144 + DN(2,2)*clhs146 + clhs211);
            lhs(9,3)=theta*(DN(0,1)*clhs307 - clhs177*clhs29 - clhs178);
            lhs(9,4)=theta*(DN(2,0)*clhs49 + DN(2,1)*clhs149 + DN(2,2)*clhs150 + clhs253);
            lhs(9,5)=clhs152*clhs302 + clhs154*clhs303 + clhs271 + clhs301*clhs65 + clhs326;
            lhs(9,6)=theta*(DN(2,0)*clhs72 + DN(2,1)*clhs158 + DN(2,2)*clhs160 + clhs286);
            lhs(9,7)=theta*(DN(1,1)*clhs307 - clhs274*clhs29 - clhs275);
            lhs(9,8)=theta*(DN(2,0)*clhs80 + DN(2,1)*clhs164 + DN(2,2)*clhs165 + clhs313);
            lhs(9,9)=clhs14*clhs327 + clhs167*clhs302 + clhs169*clhs303 + clhs301*clhs95 + clhs311;
            lhs(9,10)=theta*(DN(2,0)*clhs102 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs329);
            lhs(9,11)=clhs302*clhs315;
            lhs(9,12)=theta*(DN(2,0)*clhs110 + DN(2,1)*clhs179 + DN(2,2)*clhs180 + clhs330);
            lhs(9,13)=clhs125*clhs301 + clhs182*clhs302 + clhs184*clhs303 + clhs332 + clhs333;
            lhs(9,14)=theta*(DN(2,0)*clhs132 + DN(2,1)*clhs188 + DN(2,2)*clhs190 + clhs334);
            lhs(9,15)=theta*(DN(3,1)*clhs307 - clhs29*clhs336 - clhs335);
            lhs(10,0)=theta*(DN(2,0)*clhs6 + DN(2,1)*clhs139 + DN(2,2)*clhs194 + clhs105);
            lhs(10,1)=theta*(DN(2,0)*clhs37 + DN(2,1)*clhs142 + DN(2,2)*clhs195 + clhs176);
            lhs(10,2)=clhs146*clhs302 + clhs196*clhs303 + clhs214 + clhs301*clhs44 + clhs325;
            lhs(10,3)=theta*(DN(0,2)*clhs307 - clhs215*clhs29 - clhs216);
            lhs(10,4)=theta*(DN(2,0)*clhs51 + DN(2,1)*clhs150 + DN(2,2)*clhs198 + clhs254);
            lhs(10,5)=theta*(DN(2,0)*clhs68 + DN(2,1)*clhs154 + DN(2,2)*clhs201 + clhs273);
            lhs(10,6)=clhs160*clhs302 + clhs203*clhs303 + clhs288 + clhs301*clhs74 + clhs326;
            lhs(10,7)=theta*(DN(1,2)*clhs307 - clhs289*clhs29 - clhs290);
            lhs(10,8)=theta*(DN(2,0)*clhs82 + DN(2,1)*clhs165 + DN(2,2)*clhs208 + clhs314);
            lhs(10,9)=theta*(DN(2,0)*clhs98 + DN(2,1)*clhs169 + DN(2,2)*clhs210 + clhs329);
            lhs(10,10)=clhs104*clhs301 + clhs14*clhs337 + clhs175*clhs302 + clhs212*clhs303 + clhs311;
            lhs(10,11)=clhs303*clhs315;
            lhs(10,12)=theta*(DN(2,0)*clhs112 + DN(2,1)*clhs180 + DN(2,2)*clhs217 + clhs339);
            lhs(10,13)=theta*(DN(2,0)*clhs128 + DN(2,1)*clhs184 + DN(2,2)*clhs219 + clhs340);
            lhs(10,14)=clhs134*clhs301 + clhs190*clhs302 + clhs221*clhs303 + clhs333 + clhs342;
            lhs(10,15)=theta*(DN(3,2)*clhs307 - clhs29*clhs344 - clhs343);
            lhs(11,0)=DN(2,0)*clhs297 + clhs106*theta;
            lhs(11,1)=DN(2,1)*clhs297 + clhs177*theta;
            lhs(11,2)=DN(2,2)*clhs297 + clhs215*theta;
            lhs(11,3)=clhs231;
            lhs(11,4)=DN(2,0)*clhs228 + clhs255*theta;
            lhs(11,5)=DN(2,1)*clhs228 + clhs274*theta;
            lhs(11,6)=DN(2,2)*clhs228 + clhs289*theta;
            lhs(11,7)=clhs299;
            lhs(11,8)=DN(2,0)*clhs345;
            lhs(11,9)=DN(2,1)*clhs345;
            lhs(11,10)=DN(2,2)*clhs345;
            lhs(11,11)=clhs227*(clhs309 + clhs327 + clhs337);
            lhs(11,12)=DN(2,0)*clhs232 + DN(3,0)*clhs87;
            lhs(11,13)=DN(2,1)*clhs232 + DN(3,1)*clhs87;
            lhs(11,14)=DN(2,2)*clhs232 + DN(3,2)*clhs87;
            lhs(11,15)=clhs346;
            lhs(12,0)=clhs0*clhs347 + clhs115 + clhs3*clhs348 + clhs349*clhs6 + clhs351;
            lhs(12,1)=theta*(DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs37 + clhs181);
            lhs(12,2)=theta*(DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs44 + clhs218);
            lhs(12,3)=theta*(DN(0,0)*clhs352 - clhs136*clhs29 - clhs137);
            lhs(12,4)=clhs259 + clhs347*clhs47 + clhs348*clhs49 + clhs349*clhs51 + clhs353;
            lhs(12,5)=theta*(DN(3,0)*clhs63 + DN(3,1)*clhs65 + DN(3,2)*clhs68 + clhs276);
            lhs(12,6)=theta*(DN(3,0)*clhs70 + DN(3,1)*clhs72 + DN(3,2)*clhs74 + clhs291);
            lhs(12,7)=theta*(DN(1,0)*clhs352 - clhs263*clhs29 - clhs264);
            lhs(12,8)=clhs318 + clhs347*clhs78 + clhs348*clhs80 + clhs349*clhs82 + clhs354;
            lhs(12,9)=theta*(DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs98 + clhs330);
            lhs(12,10)=theta*(DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs339);
            lhs(12,11)=theta*(DN(2,0)*clhs352 - clhs29*clhs323 - clhs324);
            lhs(12,12)=clhs108*clhs347 + clhs110*clhs348 + clhs112*clhs349 + clhs14*clhs355 + clhs357;
            lhs(12,13)=theta*(DN(3,0)*clhs123 + DN(3,1)*clhs125 + DN(3,2)*clhs128 + clhs359);
            lhs(12,14)=theta*(DN(3,0)*clhs130 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs360);
            lhs(12,15)=clhs347*clhs361;
            lhs(13,0)=theta*(DN(3,0)*clhs3 + DN(3,1)*clhs138 + DN(3,2)*clhs139 + clhs129);
            lhs(13,1)=clhs140*clhs348 + clhs142*clhs349 + clhs186 + clhs34*clhs347 + clhs362;
            lhs(13,2)=theta*(DN(3,0)*clhs42 + DN(3,1)*clhs144 + DN(3,2)*clhs146 + clhs220);
            lhs(13,3)=theta*(DN(0,1)*clhs352 - clhs192*clhs29 - clhs193);
            lhs(13,4)=theta*(DN(3,0)*clhs49 + DN(3,1)*clhs149 + DN(3,2)*clhs150 + clhs261);
            lhs(13,5)=clhs152*clhs348 + clhs154*clhs349 + clhs278 + clhs347*clhs65 + clhs363;
            lhs(13,6)=theta*(DN(3,0)*clhs72 + DN(3,1)*clhs158 + DN(3,2)*clhs160 + clhs292);
            lhs(13,7)=theta*(DN(1,1)*clhs352 - clhs281*clhs29 - clhs282);
            lhs(13,8)=theta*(DN(3,0)*clhs80 + DN(3,1)*clhs164 + DN(3,2)*clhs165 + clhs321);
            lhs(13,9)=clhs167*clhs348 + clhs169*clhs349 + clhs332 + clhs347*clhs95 + clhs364;
            lhs(13,10)=theta*(DN(3,0)*clhs102 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs340);
            lhs(13,11)=theta*(DN(2,1)*clhs352 - clhs29*clhs335 - clhs336);
            lhs(13,12)=theta*(DN(3,0)*clhs110 + DN(3,1)*clhs179 + DN(3,2)*clhs180 + clhs359);
            lhs(13,13)=clhs125*clhs347 + clhs14*clhs365 + clhs182*clhs348 + clhs184*clhs349 + clhs357;
            lhs(13,14)=theta*(DN(3,0)*clhs132 + DN(3,1)*clhs188 + DN(3,2)*clhs190 + clhs366);
            lhs(13,15)=clhs348*clhs361;
            lhs(14,0)=theta*(DN(3,0)*clhs6 + DN(3,1)*clhs139 + DN(3,2)*clhs194 + clhs135);
            lhs(14,1)=theta*(DN(3,0)*clhs37 + DN(3,1)*clhs142 + DN(3,2)*clhs195 + clhs191);
            lhs(14,2)=clhs146*clhs348 + clhs196*clhs349 + clhs223 + clhs347*clhs44 + clhs362;
            lhs(14,3)=theta*(DN(0,2)*clhs352 - clhs224*clhs29 - clhs225);
            lhs(14,4)=theta*(DN(3,0)*clhs51 + DN(3,1)*clhs150 + DN(3,2)*clhs198 + clhs262);
            lhs(14,5)=theta*(DN(3,0)*clhs68 + DN(3,1)*clhs154 + DN(3,2)*clhs201 + clhs280);
            lhs(14,6)=clhs160*clhs348 + clhs203*clhs349 + clhs294 + clhs347*clhs74 + clhs363;
            lhs(14,7)=theta*(DN(1,2)*clhs352 - clhs29*clhs295 - clhs296);
            lhs(14,8)=theta*(DN(3,0)*clhs82 + DN(3,1)*clhs165 + DN(3,2)*clhs208 + clhs322);
            lhs(14,9)=theta*(DN(3,0)*clhs98 + DN(3,1)*clhs169 + DN(3,2)*clhs210 + clhs334);
            lhs(14,10)=clhs104*clhs347 + clhs175*clhs348 + clhs212*clhs349 + clhs342 + clhs364;
            lhs(14,11)=theta*(DN(2,2)*clhs352 - clhs29*clhs343 - clhs344);
            lhs(14,12)=theta*(DN(3,0)*clhs112 + DN(3,1)*clhs180 + DN(3,2)*clhs217 + clhs360);
            lhs(14,13)=theta*(DN(3,0)*clhs128 + DN(3,1)*clhs184 + DN(3,2)*clhs219 + clhs366);
            lhs(14,14)=clhs134*clhs347 + clhs14*clhs367 + clhs190*clhs348 + clhs221*clhs349 + clhs357;
            lhs(14,15)=clhs349*clhs361;
            lhs(15,0)=DN(3,0)*clhs297 + clhs136*theta;
            lhs(15,1)=DN(3,1)*clhs297 + clhs192*theta;
            lhs(15,2)=DN(3,2)*clhs297 + clhs224*theta;
            lhs(15,3)=clhs233;
            lhs(15,4)=DN(3,0)*clhs228 + clhs263*theta;
            lhs(15,5)=DN(3,1)*clhs228 + clhs281*theta;
            lhs(15,6)=DN(3,2)*clhs228 + clhs295*theta;
            lhs(15,7)=clhs300;
            lhs(15,8)=DN(3,0)*clhs230 + clhs323*theta;
            lhs(15,9)=DN(3,1)*clhs230 + clhs335*theta;
            lhs(15,10)=DN(3,2)*clhs230 + clhs343*theta;
            lhs(15,11)=clhs346;
            lhs(15,12)=DN(3,0)*clhs368;
            lhs(15,13)=DN(3,1)*clhs368;
            lhs(15,14)=DN(3,2)*clhs368;
            lhs(15,15)=clhs227*(clhs355 + clhs365 + clhs367);


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

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const double  theta=rData.theta;
    const BoundedMatrix<double,3,2> v_CN = theta*v+(1-theta)*vn;

    const auto vconv_CN = v_CN - vmesh;

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
        volume_error_ratio = volume_error;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs1 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0));
const double crhs2 =             1.0/dt;
const double crhs3 =             N[0]*rho;
const double crhs4 =             crhs2*crhs3;
const double crhs5 =             p[0]*theta;
const double crhs6 =             theta - 1;
const double crhs7 =             p[1]*theta;
const double crhs8 =             p[2]*theta;
const double crhs9 =             N[0]*(crhs5 - crhs6*pn[0]) + N[1]*(-crhs6*pn[1] + crhs7) + N[2]*(-crhs6*pn[2] + crhs8);
const double crhs10 =             theta*v(0,0);
const double crhs11 =             crhs10 - crhs6*vn(0,0);
const double crhs12 =             theta*v(1,0);
const double crhs13 =             crhs12 - crhs6*vn(1,0);
const double crhs14 =             theta*v(2,0);
const double crhs15 =             crhs14 - crhs6*vn(2,0);
const double crhs16 =             N[0]*crhs11 + N[1]*crhs13 + N[2]*crhs15;
const double crhs17 =             K_darcy*N[0];
const double crhs18 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double crhs19 =             DN(0,0)*crhs11 + DN(1,0)*crhs13 + DN(2,0)*crhs15;
const double crhs20 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double crhs21 =             crhs18*crhs19 + crhs20*(DN(0,1)*crhs11 + DN(1,1)*crhs13 + DN(2,1)*crhs15);
const double crhs22 =             theta*v(0,1);
const double crhs23 =             crhs22 - crhs6*vn(0,1);
const double crhs24 =             theta*v(1,1);
const double crhs25 =             crhs24 - crhs6*vn(1,1);
const double crhs26 =             theta*v(2,1);
const double crhs27 =             crhs26 - crhs6*vn(2,1);
const double crhs28 =             DN(0,1)*crhs23 + DN(1,1)*crhs25 + DN(2,1)*crhs27;
const double crhs29 =             crhs19 + crhs28 - volume_error_ratio;
const double crhs30 =             rho*stab_c2*sqrt(pow(crhs18, 2) + pow(crhs20, 2));
const double crhs31 =             crhs29*(crhs30*h/stab_c1 + mu);
const double crhs32 =             1 - theta;
const double crhs33 =             crhs32*pn[0] + crhs5;
const double crhs34 =             crhs32*pn[1] + crhs7;
const double crhs35 =             crhs32*pn[2] + crhs8;
const double crhs36 =             crhs2*rho;
const double crhs37 =             crhs1*crhs36;
const double crhs38 =             crhs10 + crhs32*vn(0,0);
const double crhs39 =             crhs12 + crhs32*vn(1,0);
const double crhs40 =             crhs14 + crhs32*vn(2,0);
const double crhs41 =             1.0/(K_darcy + crhs30/h + crhs36*dyn_tau + mu*stab_c1/pow(h, 2));
const double crhs42 =             crhs41*(DN(0,0)*crhs33 + DN(1,0)*crhs34 + DN(2,0)*crhs35 + K_darcy*(N[0]*crhs38 + N[1]*crhs39 + N[2]*crhs40) - crhs0 + crhs37 + rho*(crhs18*(DN(0,0)*crhs38 + DN(1,0)*crhs39 + DN(2,0)*crhs40) + crhs20*(DN(0,1)*crhs38 + DN(1,1)*crhs39 + DN(2,1)*crhs40)));
const double crhs43 =             rho*(DN(0,0)*crhs18 + DN(0,1)*crhs20);
const double crhs44 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs45 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1));
const double crhs46 =             N[0]*crhs23 + N[1]*crhs25 + N[2]*crhs27;
const double crhs47 =             crhs18*(DN(0,0)*crhs23 + DN(1,0)*crhs25 + DN(2,0)*crhs27) + crhs20*crhs28;
const double crhs48 =             crhs36*crhs45;
const double crhs49 =             crhs22 + crhs32*vn(0,1);
const double crhs50 =             crhs24 + crhs32*vn(1,1);
const double crhs51 =             crhs26 + crhs32*vn(2,1);
const double crhs52 =             crhs41*(DN(0,1)*crhs33 + DN(1,1)*crhs34 + DN(2,1)*crhs35 + K_darcy*(N[0]*crhs49 + N[1]*crhs50 + N[2]*crhs51) - crhs44 + crhs48 + rho*(crhs18*(DN(0,0)*crhs49 + DN(1,0)*crhs50 + DN(2,0)*crhs51) + crhs20*(DN(0,1)*crhs49 + DN(1,1)*crhs50 + DN(2,1)*crhs51)));
const double crhs53 =             K_darcy*N[1];
const double crhs54 =             N[1]*rho;
const double crhs55 =             rho*(DN(1,0)*crhs18 + DN(1,1)*crhs20);
const double crhs56 =             K_darcy*N[2];
const double crhs57 =             N[2]*rho;
const double crhs58 =             rho*(DN(2,0)*crhs18 + DN(2,1)*crhs20);
            rhs[0]=-DN(0,0)*crhs31 + DN(0,0)*crhs9 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs0 - crhs1*crhs4 - crhs16*crhs17 + crhs17*crhs42 - crhs21*crhs3 - crhs42*crhs43;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs31 + DN(0,1)*crhs9 - DN(0,1)*stress[1] + N[0]*crhs44 - crhs17*crhs46 + crhs17*crhs52 - crhs3*crhs47 - crhs4*crhs45 - crhs43*crhs52;
            rhs[2]=-DN(0,0)*crhs42 - DN(0,1)*crhs52 - N[0]*crhs29;
            rhs[3]=-DN(1,0)*crhs31 + DN(1,0)*crhs9 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs0 - N[1]*crhs37 - crhs16*crhs53 - crhs21*crhs54 + crhs42*crhs53 - crhs42*crhs55;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs31 + DN(1,1)*crhs9 - DN(1,1)*stress[1] + N[1]*crhs44 - N[1]*crhs48 - crhs46*crhs53 - crhs47*crhs54 + crhs52*crhs53 - crhs52*crhs55;
            rhs[5]=-DN(1,0)*crhs42 - DN(1,1)*crhs52 - N[1]*crhs29;
            rhs[6]=-DN(2,0)*crhs31 + DN(2,0)*crhs9 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs0 - N[2]*crhs37 - crhs16*crhs56 - crhs21*crhs57 + crhs42*crhs56 - crhs42*crhs58;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs31 + DN(2,1)*crhs9 - DN(2,1)*stress[1] + N[2]*crhs44 - N[2]*crhs48 - crhs46*crhs56 - crhs47*crhs57 + crhs52*crhs56 - crhs52*crhs58;
            rhs[8]=-DN(2,0)*crhs42 - DN(2,1)*crhs52 - N[2]*crhs29;


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

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const double  theta=rData.theta;
    const BoundedMatrix<double,4,3> v_CN = theta*v+(1-theta)*vn;
    const auto vconv_CN = v_CN - vmesh;

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
        volume_error_ratio = volume_error;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs1 =             N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)) + N[3]*(v(3,0) - vn(3,0));
const double crhs2 =             1.0/dt;
const double crhs3 =             N[0]*rho;
const double crhs4 =             crhs2*crhs3;
const double crhs5 =             p[0]*theta;
const double crhs6 =             theta - 1;
const double crhs7 =             p[1]*theta;
const double crhs8 =             p[2]*theta;
const double crhs9 =             p[3]*theta;
const double crhs10 =             N[0]*(crhs5 - crhs6*pn[0]) + N[1]*(-crhs6*pn[1] + crhs7) + N[2]*(-crhs6*pn[2] + crhs8) + N[3]*(-crhs6*pn[3] + crhs9);
const double crhs11 =             theta*v(0,0);
const double crhs12 =             crhs11 - crhs6*vn(0,0);
const double crhs13 =             theta*v(1,0);
const double crhs14 =             crhs13 - crhs6*vn(1,0);
const double crhs15 =             theta*v(2,0);
const double crhs16 =             crhs15 - crhs6*vn(2,0);
const double crhs17 =             theta*v(3,0);
const double crhs18 =             crhs17 - crhs6*vn(3,0);
const double crhs19 =             N[0]*crhs12 + N[1]*crhs14 + N[2]*crhs16 + N[3]*crhs18;
const double crhs20 =             K_darcy*N[0];
const double crhs21 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double crhs22 =             DN(0,0)*crhs12 + DN(1,0)*crhs14 + DN(2,0)*crhs16 + DN(3,0)*crhs18;
const double crhs23 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double crhs24 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double crhs25 =             crhs21*crhs22 + crhs23*(DN(0,1)*crhs12 + DN(1,1)*crhs14 + DN(2,1)*crhs16 + DN(3,1)*crhs18) + crhs24*(DN(0,2)*crhs12 + DN(1,2)*crhs14 + DN(2,2)*crhs16 + DN(3,2)*crhs18);
const double crhs26 =             theta*v(0,1);
const double crhs27 =             crhs26 - crhs6*vn(0,1);
const double crhs28 =             theta*v(1,1);
const double crhs29 =             crhs28 - crhs6*vn(1,1);
const double crhs30 =             theta*v(2,1);
const double crhs31 =             crhs30 - crhs6*vn(2,1);
const double crhs32 =             theta*v(3,1);
const double crhs33 =             crhs32 - crhs6*vn(3,1);
const double crhs34 =             DN(0,1)*crhs27 + DN(1,1)*crhs29 + DN(2,1)*crhs31 + DN(3,1)*crhs33;
const double crhs35 =             theta*v(0,2);
const double crhs36 =             crhs35 - crhs6*vn(0,2);
const double crhs37 =             theta*v(1,2);
const double crhs38 =             crhs37 - crhs6*vn(1,2);
const double crhs39 =             theta*v(2,2);
const double crhs40 =             crhs39 - crhs6*vn(2,2);
const double crhs41 =             theta*v(3,2);
const double crhs42 =             crhs41 - crhs6*vn(3,2);
const double crhs43 =             DN(0,2)*crhs36 + DN(1,2)*crhs38 + DN(2,2)*crhs40 + DN(3,2)*crhs42;
const double crhs44 =             crhs22 + crhs34 + crhs43 - volume_error_ratio;
const double crhs45 =             rho*stab_c2*sqrt(pow(crhs21, 2) + pow(crhs23, 2) + pow(crhs24, 2));
const double crhs46 =             crhs44*(crhs45*h/stab_c1 + mu);
const double crhs47 =             1 - theta;
const double crhs48 =             crhs47*pn[0] + crhs5;
const double crhs49 =             crhs47*pn[1] + crhs7;
const double crhs50 =             crhs47*pn[2] + crhs8;
const double crhs51 =             crhs47*pn[3] + crhs9;
const double crhs52 =             crhs2*rho;
const double crhs53 =             crhs1*crhs52;
const double crhs54 =             crhs11 + crhs47*vn(0,0);
const double crhs55 =             crhs13 + crhs47*vn(1,0);
const double crhs56 =             crhs15 + crhs47*vn(2,0);
const double crhs57 =             crhs17 + crhs47*vn(3,0);
const double crhs58 =             1.0/(K_darcy + crhs45/h + crhs52*dyn_tau + mu*stab_c1/pow(h, 2));
const double crhs59 =             crhs58*(DN(0,0)*crhs48 + DN(1,0)*crhs49 + DN(2,0)*crhs50 + DN(3,0)*crhs51 + K_darcy*(N[0]*crhs54 + N[1]*crhs55 + N[2]*crhs56 + N[3]*crhs57) - crhs0 + crhs53 + rho*(crhs21*(DN(0,0)*crhs54 + DN(1,0)*crhs55 + DN(2,0)*crhs56 + DN(3,0)*crhs57) + crhs23*(DN(0,1)*crhs54 + DN(1,1)*crhs55 + DN(2,1)*crhs56 + DN(3,1)*crhs57) + crhs24*(DN(0,2)*crhs54 + DN(1,2)*crhs55 + DN(2,2)*crhs56 + DN(3,2)*crhs57)));
const double crhs60 =             rho*(DN(0,0)*crhs21 + DN(0,1)*crhs23 + DN(0,2)*crhs24);
const double crhs61 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs62 =             N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)) + N[3]*(v(3,1) - vn(3,1));
const double crhs63 =             N[0]*crhs27 + N[1]*crhs29 + N[2]*crhs31 + N[3]*crhs33;
const double crhs64 =             crhs21*(DN(0,0)*crhs27 + DN(1,0)*crhs29 + DN(2,0)*crhs31 + DN(3,0)*crhs33) + crhs23*crhs34 + crhs24*(DN(0,2)*crhs27 + DN(1,2)*crhs29 + DN(2,2)*crhs31 + DN(3,2)*crhs33);
const double crhs65 =             crhs52*crhs62;
const double crhs66 =             crhs26 + crhs47*vn(0,1);
const double crhs67 =             crhs28 + crhs47*vn(1,1);
const double crhs68 =             crhs30 + crhs47*vn(2,1);
const double crhs69 =             crhs32 + crhs47*vn(3,1);
const double crhs70 =             crhs58*(DN(0,1)*crhs48 + DN(1,1)*crhs49 + DN(2,1)*crhs50 + DN(3,1)*crhs51 + K_darcy*(N[0]*crhs66 + N[1]*crhs67 + N[2]*crhs68 + N[3]*crhs69) - crhs61 + crhs65 + rho*(crhs21*(DN(0,0)*crhs66 + DN(1,0)*crhs67 + DN(2,0)*crhs68 + DN(3,0)*crhs69) + crhs23*(DN(0,1)*crhs66 + DN(1,1)*crhs67 + DN(2,1)*crhs68 + DN(3,1)*crhs69) + crhs24*(DN(0,2)*crhs66 + DN(1,2)*crhs67 + DN(2,2)*crhs68 + DN(3,2)*crhs69)));
const double crhs71 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs72 =             N[0]*(v(0,2) - vn(0,2)) + N[1]*(v(1,2) - vn(1,2)) + N[2]*(v(2,2) - vn(2,2)) + N[3]*(v(3,2) - vn(3,2));
const double crhs73 =             N[0]*crhs36 + N[1]*crhs38 + N[2]*crhs40 + N[3]*crhs42;
const double crhs74 =             crhs21*(DN(0,0)*crhs36 + DN(1,0)*crhs38 + DN(2,0)*crhs40 + DN(3,0)*crhs42) + crhs23*(DN(0,1)*crhs36 + DN(1,1)*crhs38 + DN(2,1)*crhs40 + DN(3,1)*crhs42) + crhs24*crhs43;
const double crhs75 =             crhs52*crhs72;
const double crhs76 =             crhs35 + crhs47*vn(0,2);
const double crhs77 =             crhs37 + crhs47*vn(1,2);
const double crhs78 =             crhs39 + crhs47*vn(2,2);
const double crhs79 =             crhs41 + crhs47*vn(3,2);
const double crhs80 =             crhs58*(DN(0,2)*crhs48 + DN(1,2)*crhs49 + DN(2,2)*crhs50 + DN(3,2)*crhs51 + K_darcy*(N[0]*crhs76 + N[1]*crhs77 + N[2]*crhs78 + N[3]*crhs79) - crhs71 + crhs75 + rho*(crhs21*(DN(0,0)*crhs76 + DN(1,0)*crhs77 + DN(2,0)*crhs78 + DN(3,0)*crhs79) + crhs23*(DN(0,1)*crhs76 + DN(1,1)*crhs77 + DN(2,1)*crhs78 + DN(3,1)*crhs79) + crhs24*(DN(0,2)*crhs76 + DN(1,2)*crhs77 + DN(2,2)*crhs78 + DN(3,2)*crhs79)));
const double crhs81 =             K_darcy*N[1];
const double crhs82 =             N[1]*rho;
const double crhs83 =             rho*(DN(1,0)*crhs21 + DN(1,1)*crhs23 + DN(1,2)*crhs24);
const double crhs84 =             K_darcy*N[2];
const double crhs85 =             N[2]*rho;
const double crhs86 =             rho*(DN(2,0)*crhs21 + DN(2,1)*crhs23 + DN(2,2)*crhs24);
const double crhs87 =             K_darcy*N[3];
const double crhs88 =             N[3]*rho;
const double crhs89 =             rho*(DN(3,0)*crhs21 + DN(3,1)*crhs23 + DN(3,2)*crhs24);
            rhs[0]=DN(0,0)*crhs10 - DN(0,0)*crhs46 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs0 - crhs1*crhs4 - crhs19*crhs20 + crhs20*crhs59 - crhs25*crhs3 - crhs59*crhs60;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs10 - DN(0,1)*crhs46 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs61 - crhs20*crhs63 + crhs20*crhs70 - crhs3*crhs64 - crhs4*crhs62 - crhs60*crhs70;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs10 - DN(0,2)*crhs46 - DN(0,2)*stress[2] + N[0]*crhs71 - crhs20*crhs73 + crhs20*crhs80 - crhs3*crhs74 - crhs4*crhs72 - crhs60*crhs80;
            rhs[3]=-DN(0,0)*crhs59 - DN(0,1)*crhs70 - DN(0,2)*crhs80 - N[0]*crhs44;
            rhs[4]=DN(1,0)*crhs10 - DN(1,0)*crhs46 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs0 - N[1]*crhs53 - crhs19*crhs81 - crhs25*crhs82 + crhs59*crhs81 - crhs59*crhs83;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs10 - DN(1,1)*crhs46 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs61 - N[1]*crhs65 - crhs63*crhs81 - crhs64*crhs82 + crhs70*crhs81 - crhs70*crhs83;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs10 - DN(1,2)*crhs46 - DN(1,2)*stress[2] + N[1]*crhs71 - N[1]*crhs75 - crhs73*crhs81 - crhs74*crhs82 + crhs80*crhs81 - crhs80*crhs83;
            rhs[7]=-DN(1,0)*crhs59 - DN(1,1)*crhs70 - DN(1,2)*crhs80 - N[1]*crhs44;
            rhs[8]=DN(2,0)*crhs10 - DN(2,0)*crhs46 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs0 - N[2]*crhs53 - crhs19*crhs84 - crhs25*crhs85 + crhs59*crhs84 - crhs59*crhs86;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs10 - DN(2,1)*crhs46 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs61 - N[2]*crhs65 - crhs63*crhs84 - crhs64*crhs85 + crhs70*crhs84 - crhs70*crhs86;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs10 - DN(2,2)*crhs46 - DN(2,2)*stress[2] + N[2]*crhs71 - N[2]*crhs75 - crhs73*crhs84 - crhs74*crhs85 + crhs80*crhs84 - crhs80*crhs86;
            rhs[11]=-DN(2,0)*crhs59 - DN(2,1)*crhs70 - DN(2,2)*crhs80 - N[2]*crhs44;
            rhs[12]=DN(3,0)*crhs10 - DN(3,0)*crhs46 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs0 - N[3]*crhs53 - crhs19*crhs87 - crhs25*crhs88 + crhs59*crhs87 - crhs59*crhs89;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs10 - DN(3,1)*crhs46 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs61 - N[3]*crhs65 - crhs63*crhs87 - crhs64*crhs88 + crhs70*crhs87 - crhs70*crhs89;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs10 - DN(3,2)*crhs46 - DN(3,2)*stress[2] + N[3]*crhs71 - N[3]*crhs75 - crhs73*crhs87 - crhs74*crhs88 + crhs80*crhs87 - crhs80*crhs89;
            rhs[15]=-DN(3,0)*crhs59 - DN(3,1)*crhs70 - DN(3,2)*crhs80 - N[3]*crhs44;


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

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &pn = rData.Pressure_OldStep1;
    const auto &p=rData.Pressure;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const double  theta=rData.theta;
    const BoundedMatrix<double,3,2> v_CN = theta*v+(1-theta)*vn;
    const auto vconv_CN = v_CN - vmesh;

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
        volume_error_ratio = volume_error;
    }

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double cV1 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
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


    const double cH0 =             DN(0,0)*theta;
const double cH1 =             K_darcy*theta;
const double cH2 =             1.0/dt;
const double cH3 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double cH4 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double cH5 =             DN(0,1)*theta;
const double cH6 =             1.0/(K_darcy + cH2*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 =             cH6*(N[0]*cH1 + rho*(N[0]*cH2 + cH0*cH3 + cH4*cH5));
const double cH8 =             cH6*theta;
const double cH9 =             DN(1,0)*theta;
const double cH10 =             DN(1,1)*theta;
const double cH11 =             cH6*(N[1]*cH1 + rho*(N[1]*cH2 + cH10*cH4 + cH3*cH9));
const double cH12 =             DN(2,0)*theta;
const double cH13 =             DN(2,1)*theta;
const double cH14 =             cH6*(N[2]*cH1 + rho*(N[2]*cH2 + cH12*cH3 + cH13*cH4));
            H(0,0)=DNenr(0,0)*cH7 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH7 + Nenr[0]*cH5;
            H(0,2)=cH8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DNenr(0,0)*cH11 + Nenr[0]*cH9;
            H(0,4)=DNenr(0,1)*cH11 + Nenr[0]*cH10;
            H(0,5)=cH8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DNenr(0,0)*cH14 + Nenr[0]*cH12;
            H(0,7)=DNenr(0,1)*cH14 + Nenr[0]*cH13;
            H(0,8)=cH8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DNenr(1,0)*cH7 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH7 + Nenr[1]*cH5;
            H(1,2)=cH8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DNenr(1,0)*cH11 + Nenr[1]*cH9;
            H(1,4)=DNenr(1,1)*cH11 + Nenr[1]*cH10;
            H(1,5)=cH8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DNenr(1,0)*cH14 + Nenr[1]*cH12;
            H(1,7)=DNenr(1,1)*cH14 + Nenr[1]*cH13;
            H(1,8)=cH8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DNenr(2,0)*cH7 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH7 + Nenr[2]*cH5;
            H(2,2)=cH8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DNenr(2,0)*cH11 + Nenr[2]*cH9;
            H(2,4)=DNenr(2,1)*cH11 + Nenr[2]*cH10;
            H(2,5)=cH8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DNenr(2,0)*cH14 + Nenr[2]*cH12;
            H(2,7)=DNenr(2,1)*cH14 + Nenr[2]*cH13;
            H(2,8)=cH8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0), 2) + pow(N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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


    const double crhs_ee0 =             theta*v(0,0);
const double crhs_ee1 =             theta - 1;
const double crhs_ee2 =             theta*v(0,1);
const double crhs_ee3 =             theta*v(1,0);
const double crhs_ee4 =             theta*v(1,1);
const double crhs_ee5 =             theta*v(2,0);
const double crhs_ee6 =             theta*v(2,1);
const double crhs_ee7 =             DN(0,0)*(crhs_ee0 - crhs_ee1*vn(0,0)) + DN(0,1)*(-crhs_ee1*vn(0,1) + crhs_ee2) + DN(1,0)*(-crhs_ee1*vn(1,0) + crhs_ee3) + DN(1,1)*(-crhs_ee1*vn(1,1) + crhs_ee4) + DN(2,0)*(-crhs_ee1*vn(2,0) + crhs_ee5) + DN(2,1)*(-crhs_ee1*vn(2,1) + crhs_ee6) - volume_error_ratio;
const double crhs_ee8 =             1 - theta;
const double crhs_ee9 =             crhs_ee8*pn[0] + p[0]*theta;
const double crhs_ee10 =             crhs_ee8*pn[1] + p[1]*theta;
const double crhs_ee11 =             crhs_ee8*pn[2] + p[2]*theta;
const double crhs_ee12 =             crhs_ee0 + crhs_ee8*vn(0,0);
const double crhs_ee13 =             crhs_ee3 + crhs_ee8*vn(1,0);
const double crhs_ee14 =             crhs_ee5 + crhs_ee8*vn(2,0);
const double crhs_ee15 =             1.0/dt;
const double crhs_ee16 =             N[0]*crhs_ee15;
const double crhs_ee17 =             N[1]*crhs_ee15;
const double crhs_ee18 =             N[2]*crhs_ee15;
const double crhs_ee19 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double crhs_ee20 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double crhs_ee21 =             1.0/(K_darcy + crhs_ee15*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee19, 2) + pow(crhs_ee20, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee22 =             crhs_ee21*(DN(0,0)*crhs_ee9 + DN(1,0)*crhs_ee10 + DN(2,0)*crhs_ee11 + DNenr(0,0)*penr_cn[0] + DNenr(1,0)*penr_cn[1] + DNenr(2,0)*penr_cn[2] + K_darcy*(N[0]*crhs_ee12 + N[1]*crhs_ee13 + N[2]*crhs_ee14) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(crhs_ee16*(v(0,0) - vn(0,0)) + crhs_ee17*(v(1,0) - vn(1,0)) + crhs_ee18*(v(2,0) - vn(2,0)) + crhs_ee19*(DN(0,0)*crhs_ee12 + DN(1,0)*crhs_ee13 + DN(2,0)*crhs_ee14) + crhs_ee20*(DN(0,1)*crhs_ee12 + DN(1,1)*crhs_ee13 + DN(2,1)*crhs_ee14)));
const double crhs_ee23 =             crhs_ee2 + crhs_ee8*vn(0,1);
const double crhs_ee24 =             crhs_ee4 + crhs_ee8*vn(1,1);
const double crhs_ee25 =             crhs_ee6 + crhs_ee8*vn(2,1);
const double crhs_ee26 =             crhs_ee21*(DN(0,1)*crhs_ee9 + DN(1,1)*crhs_ee10 + DN(2,1)*crhs_ee11 + DNenr(0,1)*penr_cn[0] + DNenr(1,1)*penr_cn[1] + DNenr(2,1)*penr_cn[2] + K_darcy*(N[0]*crhs_ee23 + N[1]*crhs_ee24 + N[2]*crhs_ee25) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(crhs_ee16*(v(0,1) - vn(0,1)) + crhs_ee17*(v(1,1) - vn(1,1)) + crhs_ee18*(v(2,1) - vn(2,1)) + crhs_ee19*(DN(0,0)*crhs_ee23 + DN(1,0)*crhs_ee24 + DN(2,0)*crhs_ee25) + crhs_ee20*(DN(0,1)*crhs_ee23 + DN(1,1)*crhs_ee24 + DN(2,1)*crhs_ee25)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee22 - DNenr(0,1)*crhs_ee26 - Nenr[0]*crhs_ee7;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee22 - DNenr(1,1)*crhs_ee26 - Nenr[1]*crhs_ee7;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee22 - DNenr(2,1)*crhs_ee26 - Nenr[2]*crhs_ee7;


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

    const double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &pn = rData.Pressure_OldStep1;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    // FIXME: Mesh velocity should be evaluated at 0.5dt
    const double theta=rData.theta;
    const BoundedMatrix<double,4,3> v_CN = theta*v+(1-theta)*vn;
    const auto vconv_CN = v_CN - vmesh;

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
        volume_error_ratio = volume_error;
    }

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it
    array_1d<double, NumNodes> penr_cn = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double cV1 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double cV2 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
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


    const double cH0 =             DN(0,0)*theta;
const double cH1 =             K_darcy*theta;
const double cH2 =             1.0/dt;
const double cH3 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double cH4 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double cH5 =             DN(0,1)*theta;
const double cH6 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double cH7 =             DN(0,2)*theta;
const double cH8 =             1.0/(K_darcy + cH2*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2) + pow(cH6, 2))/h + mu*stab_c1/pow(h, 2));
const double cH9 =             cH8*(N[0]*cH1 + rho*(N[0]*cH2 + cH0*cH3 + cH4*cH5 + cH6*cH7));
const double cH10 =             cH8*theta;
const double cH11 =             DN(1,0)*theta;
const double cH12 =             DN(1,1)*theta;
const double cH13 =             DN(1,2)*theta;
const double cH14 =             cH8*(N[1]*cH1 + rho*(N[1]*cH2 + cH11*cH3 + cH12*cH4 + cH13*cH6));
const double cH15 =             DN(2,0)*theta;
const double cH16 =             DN(2,1)*theta;
const double cH17 =             DN(2,2)*theta;
const double cH18 =             cH8*(N[2]*cH1 + rho*(N[2]*cH2 + cH15*cH3 + cH16*cH4 + cH17*cH6));
const double cH19 =             DN(3,0)*theta;
const double cH20 =             DN(3,1)*theta;
const double cH21 =             DN(3,2)*theta;
const double cH22 =             cH8*(N[3]*cH1 + rho*(N[3]*cH2 + cH19*cH3 + cH20*cH4 + cH21*cH6));
            H(0,0)=DNenr(0,0)*cH9 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH9 + Nenr[0]*cH5;
            H(0,2)=DNenr(0,2)*cH9 + Nenr[0]*cH7;
            H(0,3)=cH10*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DNenr(0,0)*cH14 + Nenr[0]*cH11;
            H(0,5)=DNenr(0,1)*cH14 + Nenr[0]*cH12;
            H(0,6)=DNenr(0,2)*cH14 + Nenr[0]*cH13;
            H(0,7)=cH10*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DNenr(0,0)*cH18 + Nenr[0]*cH15;
            H(0,9)=DNenr(0,1)*cH18 + Nenr[0]*cH16;
            H(0,10)=DNenr(0,2)*cH18 + Nenr[0]*cH17;
            H(0,11)=cH10*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DNenr(0,0)*cH22 + Nenr[0]*cH19;
            H(0,13)=DNenr(0,1)*cH22 + Nenr[0]*cH20;
            H(0,14)=DNenr(0,2)*cH22 + Nenr[0]*cH21;
            H(0,15)=cH10*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DNenr(1,0)*cH9 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH9 + Nenr[1]*cH5;
            H(1,2)=DNenr(1,2)*cH9 + Nenr[1]*cH7;
            H(1,3)=cH10*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DNenr(1,0)*cH14 + Nenr[1]*cH11;
            H(1,5)=DNenr(1,1)*cH14 + Nenr[1]*cH12;
            H(1,6)=DNenr(1,2)*cH14 + Nenr[1]*cH13;
            H(1,7)=cH10*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DNenr(1,0)*cH18 + Nenr[1]*cH15;
            H(1,9)=DNenr(1,1)*cH18 + Nenr[1]*cH16;
            H(1,10)=DNenr(1,2)*cH18 + Nenr[1]*cH17;
            H(1,11)=cH10*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DNenr(1,0)*cH22 + Nenr[1]*cH19;
            H(1,13)=DNenr(1,1)*cH22 + Nenr[1]*cH20;
            H(1,14)=DNenr(1,2)*cH22 + Nenr[1]*cH21;
            H(1,15)=cH10*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DNenr(2,0)*cH9 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH9 + Nenr[2]*cH5;
            H(2,2)=DNenr(2,2)*cH9 + Nenr[2]*cH7;
            H(2,3)=cH10*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DNenr(2,0)*cH14 + Nenr[2]*cH11;
            H(2,5)=DNenr(2,1)*cH14 + Nenr[2]*cH12;
            H(2,6)=DNenr(2,2)*cH14 + Nenr[2]*cH13;
            H(2,7)=cH10*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DNenr(2,0)*cH18 + Nenr[2]*cH15;
            H(2,9)=DNenr(2,1)*cH18 + Nenr[2]*cH16;
            H(2,10)=DNenr(2,2)*cH18 + Nenr[2]*cH17;
            H(2,11)=cH10*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DNenr(2,0)*cH22 + Nenr[2]*cH19;
            H(2,13)=DNenr(2,1)*cH22 + Nenr[2]*cH20;
            H(2,14)=DNenr(2,2)*cH22 + Nenr[2]*cH21;
            H(2,15)=cH10*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DNenr(3,0)*cH9 + Nenr[3]*cH0;
            H(3,1)=DNenr(3,1)*cH9 + Nenr[3]*cH5;
            H(3,2)=DNenr(3,2)*cH9 + Nenr[3]*cH7;
            H(3,3)=cH10*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DNenr(3,0)*cH14 + Nenr[3]*cH11;
            H(3,5)=DNenr(3,1)*cH14 + Nenr[3]*cH12;
            H(3,6)=DNenr(3,2)*cH14 + Nenr[3]*cH13;
            H(3,7)=cH10*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DNenr(3,0)*cH18 + Nenr[3]*cH15;
            H(3,9)=DNenr(3,1)*cH18 + Nenr[3]*cH16;
            H(3,10)=DNenr(3,2)*cH18 + Nenr[3]*cH17;
            H(3,11)=cH10*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DNenr(3,0)*cH22 + Nenr[3]*cH19;
            H(3,13)=DNenr(3,1)*cH22 + Nenr[3]*cH20;
            H(3,14)=DNenr(3,2)*cH22 + Nenr[3]*cH21;
            H(3,15)=cH10*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0), 2) + pow(N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1), 2) + pow(N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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


    const double crhs_ee0 =             theta*v(0,0);
const double crhs_ee1 =             theta - 1;
const double crhs_ee2 =             theta*v(0,1);
const double crhs_ee3 =             theta*v(0,2);
const double crhs_ee4 =             theta*v(1,0);
const double crhs_ee5 =             theta*v(1,1);
const double crhs_ee6 =             theta*v(1,2);
const double crhs_ee7 =             theta*v(2,0);
const double crhs_ee8 =             theta*v(2,1);
const double crhs_ee9 =             theta*v(2,2);
const double crhs_ee10 =             theta*v(3,0);
const double crhs_ee11 =             theta*v(3,1);
const double crhs_ee12 =             theta*v(3,2);
const double crhs_ee13 =             DN(0,0)*(crhs_ee0 - crhs_ee1*vn(0,0)) + DN(0,1)*(-crhs_ee1*vn(0,1) + crhs_ee2) + DN(0,2)*(-crhs_ee1*vn(0,2) + crhs_ee3) + DN(1,0)*(-crhs_ee1*vn(1,0) + crhs_ee4) + DN(1,1)*(-crhs_ee1*vn(1,1) + crhs_ee5) + DN(1,2)*(-crhs_ee1*vn(1,2) + crhs_ee6) + DN(2,0)*(-crhs_ee1*vn(2,0) + crhs_ee7) + DN(2,1)*(-crhs_ee1*vn(2,1) + crhs_ee8) + DN(2,2)*(-crhs_ee1*vn(2,2) + crhs_ee9) + DN(3,0)*(-crhs_ee1*vn(3,0) + crhs_ee10) + DN(3,1)*(-crhs_ee1*vn(3,1) + crhs_ee11) + DN(3,2)*(-crhs_ee1*vn(3,2) + crhs_ee12) - volume_error_ratio;
const double crhs_ee14 =             1 - theta;
const double crhs_ee15 =             crhs_ee14*pn[0] + p[0]*theta;
const double crhs_ee16 =             crhs_ee14*pn[1] + p[1]*theta;
const double crhs_ee17 =             crhs_ee14*pn[2] + p[2]*theta;
const double crhs_ee18 =             crhs_ee14*pn[3] + p[3]*theta;
const double crhs_ee19 =             crhs_ee0 + crhs_ee14*vn(0,0);
const double crhs_ee20 =             crhs_ee14*vn(1,0) + crhs_ee4;
const double crhs_ee21 =             crhs_ee14*vn(2,0) + crhs_ee7;
const double crhs_ee22 =             crhs_ee10 + crhs_ee14*vn(3,0);
const double crhs_ee23 =             1.0/dt;
const double crhs_ee24 =             N[0]*crhs_ee23;
const double crhs_ee25 =             N[1]*crhs_ee23;
const double crhs_ee26 =             N[2]*crhs_ee23;
const double crhs_ee27 =             N[3]*crhs_ee23;
const double crhs_ee28 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double crhs_ee29 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double crhs_ee30 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double crhs_ee31 =             1.0/(K_darcy + crhs_ee23*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee28, 2) + pow(crhs_ee29, 2) + pow(crhs_ee30, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee32 =             crhs_ee31*(DN(0,0)*crhs_ee15 + DN(1,0)*crhs_ee16 + DN(2,0)*crhs_ee17 + DN(3,0)*crhs_ee18 + DNenr(0,0)*penr_cn[0] + DNenr(1,0)*penr_cn[1] + DNenr(2,0)*penr_cn[2] + DNenr(3,0)*penr_cn[3] + K_darcy*(N[0]*crhs_ee19 + N[1]*crhs_ee20 + N[2]*crhs_ee21 + N[3]*crhs_ee22) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(crhs_ee24*(v(0,0) - vn(0,0)) + crhs_ee25*(v(1,0) - vn(1,0)) + crhs_ee26*(v(2,0) - vn(2,0)) + crhs_ee27*(v(3,0) - vn(3,0)) + crhs_ee28*(DN(0,0)*crhs_ee19 + DN(1,0)*crhs_ee20 + DN(2,0)*crhs_ee21 + DN(3,0)*crhs_ee22) + crhs_ee29*(DN(0,1)*crhs_ee19 + DN(1,1)*crhs_ee20 + DN(2,1)*crhs_ee21 + DN(3,1)*crhs_ee22) + crhs_ee30*(DN(0,2)*crhs_ee19 + DN(1,2)*crhs_ee20 + DN(2,2)*crhs_ee21 + DN(3,2)*crhs_ee22)));
const double crhs_ee33 =             crhs_ee14*vn(0,1) + crhs_ee2;
const double crhs_ee34 =             crhs_ee14*vn(1,1) + crhs_ee5;
const double crhs_ee35 =             crhs_ee14*vn(2,1) + crhs_ee8;
const double crhs_ee36 =             crhs_ee11 + crhs_ee14*vn(3,1);
const double crhs_ee37 =             crhs_ee31*(DN(0,1)*crhs_ee15 + DN(1,1)*crhs_ee16 + DN(2,1)*crhs_ee17 + DN(3,1)*crhs_ee18 + DNenr(0,1)*penr_cn[0] + DNenr(1,1)*penr_cn[1] + DNenr(2,1)*penr_cn[2] + DNenr(3,1)*penr_cn[3] + K_darcy*(N[0]*crhs_ee33 + N[1]*crhs_ee34 + N[2]*crhs_ee35 + N[3]*crhs_ee36) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(crhs_ee24*(v(0,1) - vn(0,1)) + crhs_ee25*(v(1,1) - vn(1,1)) + crhs_ee26*(v(2,1) - vn(2,1)) + crhs_ee27*(v(3,1) - vn(3,1)) + crhs_ee28*(DN(0,0)*crhs_ee33 + DN(1,0)*crhs_ee34 + DN(2,0)*crhs_ee35 + DN(3,0)*crhs_ee36) + crhs_ee29*(DN(0,1)*crhs_ee33 + DN(1,1)*crhs_ee34 + DN(2,1)*crhs_ee35 + DN(3,1)*crhs_ee36) + crhs_ee30*(DN(0,2)*crhs_ee33 + DN(1,2)*crhs_ee34 + DN(2,2)*crhs_ee35 + DN(3,2)*crhs_ee36)));
const double crhs_ee38 =             crhs_ee14*vn(0,2) + crhs_ee3;
const double crhs_ee39 =             crhs_ee14*vn(1,2) + crhs_ee6;
const double crhs_ee40 =             crhs_ee14*vn(2,2) + crhs_ee9;
const double crhs_ee41 =             crhs_ee12 + crhs_ee14*vn(3,2);
const double crhs_ee42 =             crhs_ee31*(DN(0,2)*crhs_ee15 + DN(1,2)*crhs_ee16 + DN(2,2)*crhs_ee17 + DN(3,2)*crhs_ee18 + DNenr(0,2)*penr_cn[0] + DNenr(1,2)*penr_cn[1] + DNenr(2,2)*penr_cn[2] + DNenr(3,2)*penr_cn[3] + K_darcy*(N[0]*crhs_ee38 + N[1]*crhs_ee39 + N[2]*crhs_ee40 + N[3]*crhs_ee41) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(crhs_ee24*(v(0,2) - vn(0,2)) + crhs_ee25*(v(1,2) - vn(1,2)) + crhs_ee26*(v(2,2) - vn(2,2)) + crhs_ee27*(v(3,2) - vn(3,2)) + crhs_ee28*(DN(0,0)*crhs_ee38 + DN(1,0)*crhs_ee39 + DN(2,0)*crhs_ee40 + DN(3,0)*crhs_ee41) + crhs_ee29*(DN(0,1)*crhs_ee38 + DN(1,1)*crhs_ee39 + DN(2,1)*crhs_ee40 + DN(3,1)*crhs_ee41) + crhs_ee30*(DN(0,2)*crhs_ee38 + DN(1,2)*crhs_ee39 + DN(2,2)*crhs_ee40 + DN(3,2)*crhs_ee41)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee32 - DNenr(0,1)*crhs_ee37 - DNenr(0,2)*crhs_ee42 - Nenr[0]*crhs_ee13;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee32 - DNenr(1,1)*crhs_ee37 - DNenr(1,2)*crhs_ee42 - Nenr[1]*crhs_ee13;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee32 - DNenr(2,1)*crhs_ee37 - DNenr(2,2)*crhs_ee42 - Nenr[2]*crhs_ee13;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee32 - DNenr(3,1)*crhs_ee37 - DNenr(3,2)*crhs_ee42 - Nenr[3]*crhs_ee13;


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

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const double  theta=rData.theta;
    const BoundedMatrix<double,NumNodes,Dim> v_CN = theta*v+(1-theta)*vn;
    const auto vmesh=rData.MeshVelocity;
    const auto v_convection_CN = v_CN - vmesh;

    for (unsigned int gp = 0; gp < rInterfaceWeights.size(); ++gp){

        Vector vconv_CN = ZeroVector(Dim);
        double positive_weight = 0.0;
        double negative_weight = 0.0;

        for (unsigned int j = 0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < Dim; ++dim){
                vconv_CN[dim] += (rEnrInterfaceShapeFunctionNeg(gp, j) + rEnrInterfaceShapeFunctionPos(gp, j))
                    *v_convection_CN(j,dim);
            }
            positive_weight += rEnrInterfaceShapeFunctionNeg(gp, j);
            negative_weight += rEnrInterfaceShapeFunctionPos(gp, j);
        }

        const double v_conv_CN_norm = norm_2(vconv_CN);

        const double penalty_coefficient = cut_stabilization_coefficient *
            density * 1.0 / (dyn_tau * density / (0.5*dt) + stab_c1 * viscosity / h_elem / h_elem +
                                stab_c2 * density * v_conv_CN_norm / h_elem) * element_volume / cut_area;

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
