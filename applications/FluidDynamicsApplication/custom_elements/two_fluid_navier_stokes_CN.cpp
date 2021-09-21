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
#include "custom_utilities/two_fluid_navier_stokes_data.h"

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

template <class TElementData>
void TwoFluidNavierStokesCN<TElementData>::CalculateStrainRate(TElementData& rData) const
{
    const typename TElementData::NodalVectorData mid_step_velocity = 0.5 * (rData.Velocity + rData.Velocity_OldStep1);

    Internals::StrainRateSpecialization<TElementData,TElementData::Dim>::Calculate(
        rData.StrainRate,
        mid_step_velocity,
        rData.DN_DX);
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointLHSContribution(
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

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &v_CN = 0.5 * (vn + v);
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
const double clhs1 =             0.5*DN(0,0);
const double clhs2 =             C(0,2)*DN(0,0);
const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
const double clhs4 =             0.5*DN(0,1);
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double clhs7 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double clhs8 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2));
const double clhs9 =             clhs8*h/stab_c1 + mu;
const double clhs10 =             0.5*clhs9;
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             0.5*K_darcy;
const double clhs13 =             rho/dt;
const double clhs14 =             0.5*N[0];
const double clhs15 =             rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7);
const double clhs16 =             1.0/(K_darcy + clhs13*dyn_tau + clhs8/h + mu*stab_c1/pow(h, 2));
const double clhs17 =             clhs15*clhs16;
const double clhs18 =             N[0]*clhs13;
const double clhs19 =             K_darcy*clhs14;
const double clhs20 =             0.5*clhs15 + clhs19;
const double clhs21 =             clhs18 + clhs20;
const double clhs22 =             1.0*clhs21;
const double clhs23 =             K_darcy*clhs16;
const double clhs24 =             N[0]*clhs23;
const double clhs25 =             clhs11*clhs12 + clhs11*clhs13 + clhs14*clhs15 + clhs17*clhs22 - clhs22*clhs24;
const double clhs26 =             C(0,1)*DN(0,1) + clhs2;
const double clhs27 =             C(1,2)*DN(0,1);
const double clhs28 =             C(2,2)*DN(0,0) + clhs27;
const double clhs29 =             DN(0,0)*clhs9;
const double clhs30 =             DN(0,1)*clhs29;
const double clhs31 =             -N[0] + clhs17 - clhs24;
const double clhs32 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs33 =             C(0,2)*DN(1,0);
const double clhs34 =             C(2,2)*DN(1,1) + clhs33;
const double clhs35 =             DN(0,0)*DN(1,0);
const double clhs36 =             N[1]*clhs18 + N[1]*clhs19;
const double clhs37 =             clhs10*clhs35 + clhs36;
const double clhs38 =             DN(1,0)*clhs6 + DN(1,1)*clhs7;
const double clhs39 =             clhs14*rho;
const double clhs40 =             N[1]*clhs13;
const double clhs41 =             0.5*N[1];
const double clhs42 =             K_darcy*clhs41;
const double clhs43 =             0.5*rho;
const double clhs44 =             clhs38*clhs43 + clhs42;
const double clhs45 =             clhs40 + clhs44;
const double clhs46 =             1.0*clhs45;
const double clhs47 =             clhs17*clhs46 - clhs24*clhs46 + clhs38*clhs39;
const double clhs48 =             C(0,1)*DN(1,1) + clhs33;
const double clhs49 =             C(1,2)*DN(1,1);
const double clhs50 =             C(2,2)*DN(1,0) + clhs49;
const double clhs51 =             DN(1,1)*clhs29;
const double clhs52 =             DN(0,0)*N[1];
const double clhs53 =             DN(1,0)*N[0];
const double clhs54 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs55 =             C(0,2)*DN(2,0);
const double clhs56 =             C(2,2)*DN(2,1) + clhs55;
const double clhs57 =             DN(0,0)*DN(2,0);
const double clhs58 =             N[2]*clhs18 + N[2]*clhs19;
const double clhs59 =             clhs10*clhs57 + clhs58;
const double clhs60 =             DN(2,0)*clhs6 + DN(2,1)*clhs7;
const double clhs61 =             N[2]*clhs13;
const double clhs62 =             0.5*N[2];
const double clhs63 =             K_darcy*clhs62 + clhs43*clhs60;
const double clhs64 =             clhs61 + clhs63;
const double clhs65 =             1.0*clhs64;
const double clhs66 =             clhs17*clhs65 - clhs24*clhs65 + clhs39*clhs60;
const double clhs67 =             C(0,1)*DN(2,1) + clhs55;
const double clhs68 =             C(1,2)*DN(2,1);
const double clhs69 =             C(2,2)*DN(2,0) + clhs68;
const double clhs70 =             DN(2,1)*clhs29;
const double clhs71 =             DN(0,0)*N[2];
const double clhs72 =             DN(2,0)*N[0];
const double clhs73 =             C(0,1)*DN(0,0) + clhs27;
const double clhs74 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs75 =             pow(DN(0,1), 2);
const double clhs76 =             C(0,1)*DN(1,0) + clhs49;
const double clhs77 =             DN(0,1)*clhs9;
const double clhs78 =             DN(1,0)*clhs77;
const double clhs79 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs80 =             DN(0,1)*DN(1,1);
const double clhs81 =             clhs10*clhs80 + clhs36;
const double clhs82 =             DN(0,1)*N[1];
const double clhs83 =             DN(1,1)*N[0];
const double clhs84 =             C(0,1)*DN(2,0) + clhs68;
const double clhs85 =             DN(2,0)*clhs77;
const double clhs86 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs87 =             DN(0,1)*DN(2,1);
const double clhs88 =             clhs10*clhs87 + clhs58;
const double clhs89 =             DN(0,1)*N[2];
const double clhs90 =             DN(2,1)*N[0];
const double clhs91 =             clhs14 + clhs16*(1.0*clhs18 + clhs20);
const double clhs92 =             0.5*clhs16;
const double clhs93 =             1.0*clhs16;
const double clhs94 =             clhs45*clhs93;
const double clhs95 =             clhs92*(clhs35 + clhs80);
const double clhs96 =             clhs64*clhs93;
const double clhs97 =             clhs92*(clhs57 + clhs87);
const double clhs98 =             0.5*DN(1,0);
const double clhs99 =             0.5*DN(1,1);
const double clhs100 =             clhs16*rho;
const double clhs101 =             clhs100*clhs38;
const double clhs102 =             N[1]*clhs23;
const double clhs103 =             clhs101*clhs22 - clhs102*clhs22 + clhs15*clhs41;
const double clhs104 =             pow(DN(1,0), 2);
const double clhs105 =             pow(N[1], 2);
const double clhs106 =             clhs41*rho;
const double clhs107 =             clhs101*clhs46 - clhs102*clhs46 + clhs105*clhs12 + clhs105*clhs13 + clhs106*clhs38;
const double clhs108 =             DN(1,0)*clhs9;
const double clhs109 =             DN(1,1)*clhs108;
const double clhs110 =             -N[1] + clhs101 - clhs102;
const double clhs111 =             DN(1,0)*DN(2,0);
const double clhs112 =             N[2]*clhs40 + N[2]*clhs42;
const double clhs113 =             clhs10*clhs111 + clhs112;
const double clhs114 =             clhs101*clhs65 - clhs102*clhs65 + clhs106*clhs60;
const double clhs115 =             DN(2,1)*clhs108;
const double clhs116 =             DN(1,0)*N[2];
const double clhs117 =             DN(2,0)*N[1];
const double clhs118 =             pow(DN(1,1), 2);
const double clhs119 =             DN(2,0)*clhs9;
const double clhs120 =             DN(1,1)*clhs119;
const double clhs121 =             DN(1,1)*DN(2,1);
const double clhs122 =             clhs10*clhs121 + clhs112;
const double clhs123 =             DN(1,1)*N[2];
const double clhs124 =             DN(2,1)*N[1];
const double clhs125 =             clhs21*clhs93;
const double clhs126 =             clhs16*(1.0*clhs40 + clhs44) + clhs41;
const double clhs127 =             clhs92*(clhs111 + clhs121);
const double clhs128 =             0.5*DN(2,0);
const double clhs129 =             0.5*DN(2,1);
const double clhs130 =             clhs100*clhs60;
const double clhs131 =             N[2]*clhs23;
const double clhs132 =             clhs130*clhs22 - clhs131*clhs22 + clhs15*clhs62;
const double clhs133 =             clhs62*rho;
const double clhs134 =             clhs130*clhs46 - clhs131*clhs46 + clhs133*clhs38;
const double clhs135 =             pow(DN(2,0), 2);
const double clhs136 =             pow(N[2], 2);
const double clhs137 =             clhs12*clhs136 + clhs13*clhs136 + clhs130*clhs65 - clhs131*clhs65 + clhs133*clhs60;
const double clhs138 =             DN(2,1)*clhs119;
const double clhs139 =             -N[2] + clhs130 - clhs131;
const double clhs140 =             pow(DN(2,1), 2);
const double clhs141 =             clhs16*(1.0*clhs61 + clhs63) + clhs62;
            lhs(0,0)=clhs0*clhs1 + clhs10*clhs5 + clhs25 + clhs3*clhs4;
            lhs(0,1)=0.5*DN(0,0)*clhs26 + 0.5*DN(0,1)*clhs28 + 0.5*clhs30;
            lhs(0,2)=clhs1*clhs31;
            lhs(0,3)=clhs1*clhs32 + clhs34*clhs4 + clhs37 + clhs47;
            lhs(0,4)=0.5*DN(0,0)*clhs48 + 0.5*DN(0,1)*clhs50 + 0.5*clhs51;
            lhs(0,5)=0.5*DN(1,0)*clhs17 - 0.5*clhs23*clhs53 - 0.5*clhs52;
            lhs(0,6)=clhs1*clhs54 + clhs4*clhs56 + clhs59 + clhs66;
            lhs(0,7)=0.5*DN(0,0)*clhs67 + 0.5*DN(0,1)*clhs69 + 0.5*clhs70;
            lhs(0,8)=0.5*DN(2,0)*clhs17 - 0.5*clhs23*clhs72 - 0.5*clhs71;
            lhs(1,0)=0.5*DN(0,0)*clhs3 + 0.5*DN(0,1)*clhs73 + 0.5*clhs30;
            lhs(1,1)=clhs1*clhs28 + clhs10*clhs75 + clhs25 + clhs4*clhs74;
            lhs(1,2)=clhs31*clhs4;
            lhs(1,3)=0.5*DN(0,0)*clhs34 + 0.5*DN(0,1)*clhs76 + 0.5*clhs78;
            lhs(1,4)=clhs1*clhs50 + clhs4*clhs79 + clhs47 + clhs81;
            lhs(1,5)=0.5*DN(1,1)*clhs17 - 0.5*clhs23*clhs83 - 0.5*clhs82;
            lhs(1,6)=0.5*DN(0,0)*clhs56 + 0.5*DN(0,1)*clhs84 + 0.5*clhs85;
            lhs(1,7)=clhs1*clhs69 + clhs4*clhs86 + clhs66 + clhs88;
            lhs(1,8)=0.5*DN(2,1)*clhs17 - 0.5*clhs23*clhs90 - 0.5*clhs89;
            lhs(2,0)=DN(0,0)*clhs91;
            lhs(2,1)=DN(0,1)*clhs91;
            lhs(2,2)=clhs92*(clhs5 + clhs75);
            lhs(2,3)=DN(0,0)*clhs94 + DN(1,0)*clhs14;
            lhs(2,4)=DN(0,1)*clhs94 + DN(1,1)*clhs14;
            lhs(2,5)=clhs95;
            lhs(2,6)=DN(0,0)*clhs96 + DN(2,0)*clhs14;
            lhs(2,7)=DN(0,1)*clhs96 + DN(2,1)*clhs14;
            lhs(2,8)=clhs97;
            lhs(3,0)=clhs0*clhs98 + clhs103 + clhs3*clhs99 + clhs37;
            lhs(3,1)=0.5*DN(1,0)*clhs26 + 0.5*DN(1,1)*clhs28 + 0.5*clhs78;
            lhs(3,2)=0.5*DN(0,0)*clhs101 - 0.5*clhs23*clhs52 - 0.5*clhs53;
            lhs(3,3)=clhs10*clhs104 + clhs107 + clhs32*clhs98 + clhs34*clhs99;
            lhs(3,4)=0.5*DN(1,0)*clhs48 + 0.5*DN(1,1)*clhs50 + 0.5*clhs109;
            lhs(3,5)=clhs110*clhs98;
            lhs(3,6)=clhs113 + clhs114 + clhs54*clhs98 + clhs56*clhs99;
            lhs(3,7)=0.5*DN(1,0)*clhs67 + 0.5*DN(1,1)*clhs69 + 0.5*clhs115;
            lhs(3,8)=0.5*DN(2,0)*clhs101 - 0.5*clhs116 - 0.5*clhs117*clhs23;
            lhs(4,0)=0.5*DN(1,0)*clhs3 + 0.5*DN(1,1)*clhs73 + 0.5*clhs51;
            lhs(4,1)=clhs103 + clhs28*clhs98 + clhs74*clhs99 + clhs81;
            lhs(4,2)=0.5*DN(0,1)*clhs101 - 0.5*clhs23*clhs82 - 0.5*clhs83;
            lhs(4,3)=0.5*DN(1,0)*clhs34 + 0.5*DN(1,1)*clhs76 + 0.5*clhs109;
            lhs(4,4)=clhs10*clhs118 + clhs107 + clhs50*clhs98 + clhs79*clhs99;
            lhs(4,5)=clhs110*clhs99;
            lhs(4,6)=0.5*DN(1,0)*clhs56 + 0.5*DN(1,1)*clhs84 + 0.5*clhs120;
            lhs(4,7)=clhs114 + clhs122 + clhs69*clhs98 + clhs86*clhs99;
            lhs(4,8)=0.5*DN(2,1)*clhs101 - 0.5*clhs123 - 0.5*clhs124*clhs23;
            lhs(5,0)=DN(1,0)*clhs125 + 0.5*clhs52;
            lhs(5,1)=DN(1,1)*clhs125 + 0.5*clhs82;
            lhs(5,2)=clhs95;
            lhs(5,3)=DN(1,0)*clhs126;
            lhs(5,4)=DN(1,1)*clhs126;
            lhs(5,5)=clhs92*(clhs104 + clhs118);
            lhs(5,6)=DN(1,0)*clhs96 + DN(2,0)*clhs41;
            lhs(5,7)=DN(1,1)*clhs96 + DN(2,1)*clhs41;
            lhs(5,8)=clhs127;
            lhs(6,0)=clhs0*clhs128 + clhs129*clhs3 + clhs132 + clhs59;
            lhs(6,1)=0.5*DN(2,0)*clhs26 + 0.5*DN(2,1)*clhs28 + 0.5*clhs85;
            lhs(6,2)=0.5*DN(0,0)*clhs130 - 0.5*clhs23*clhs71 - 0.5*clhs72;
            lhs(6,3)=clhs113 + clhs128*clhs32 + clhs129*clhs34 + clhs134;
            lhs(6,4)=0.5*DN(2,0)*clhs48 + 0.5*DN(2,1)*clhs50 + 0.5*clhs120;
            lhs(6,5)=0.5*DN(1,0)*clhs130 - 0.5*clhs116*clhs23 - 0.5*clhs117;
            lhs(6,6)=clhs10*clhs135 + clhs128*clhs54 + clhs129*clhs56 + clhs137;
            lhs(6,7)=0.5*DN(2,0)*clhs67 + 0.5*DN(2,1)*clhs69 + 0.5*clhs138;
            lhs(6,8)=clhs128*clhs139;
            lhs(7,0)=0.5*DN(2,0)*clhs3 + 0.5*DN(2,1)*clhs73 + 0.5*clhs70;
            lhs(7,1)=clhs128*clhs28 + clhs129*clhs74 + clhs132 + clhs88;
            lhs(7,2)=0.5*DN(0,1)*clhs130 - 0.5*clhs23*clhs89 - 0.5*clhs90;
            lhs(7,3)=0.5*DN(2,0)*clhs34 + 0.5*DN(2,1)*clhs76 + 0.5*clhs115;
            lhs(7,4)=clhs122 + clhs128*clhs50 + clhs129*clhs79 + clhs134;
            lhs(7,5)=0.5*DN(1,1)*clhs130 - 0.5*clhs123*clhs23 - 0.5*clhs124;
            lhs(7,6)=0.5*DN(2,0)*clhs56 + 0.5*DN(2,1)*clhs84 + 0.5*clhs138;
            lhs(7,7)=clhs10*clhs140 + clhs128*clhs69 + clhs129*clhs86 + clhs137;
            lhs(7,8)=clhs129*clhs139;
            lhs(8,0)=DN(2,0)*clhs125 + 0.5*clhs71;
            lhs(8,1)=DN(2,1)*clhs125 + 0.5*clhs89;
            lhs(8,2)=clhs97;
            lhs(8,3)=DN(2,0)*clhs94 + 0.5*clhs116;
            lhs(8,4)=DN(2,1)*clhs94 + 0.5*clhs123;
            lhs(8,5)=clhs127;
            lhs(8,6)=DN(2,0)*clhs141;
            lhs(8,7)=DN(2,1)*clhs141;
            lhs(8,8)=clhs92*(clhs135 + clhs140);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointLHSContribution(
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

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto &v = rData.Velocity;
    const auto &vn = rData.Velocity_OldStep1;
    const auto &vmesh = rData.MeshVelocity;
    const auto &v_CN = 0.5 * (vn + v);
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
const double clhs1 =             0.5*DN(0,0);
const double clhs2 =             C(0,3)*DN(0,0);
const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
const double clhs4 =             0.5*DN(0,1);
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             0.5*DN(0,2);
const double clhs8 =             pow(DN(0,0), 2);
const double clhs9 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double clhs10 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double clhs11 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double clhs12 =             rho*stab_c2*sqrt(pow(clhs10, 2) + pow(clhs11, 2) + pow(clhs9, 2));
const double clhs13 =             clhs12*h/stab_c1 + mu;
const double clhs14 =             0.5*clhs13;
const double clhs15 =             pow(N[0], 2);
const double clhs16 =             0.5*K_darcy;
const double clhs17 =             rho/dt;
const double clhs18 =             0.5*N[0];
const double clhs19 =             rho*(DN(0,0)*clhs9 + DN(0,1)*clhs10 + DN(0,2)*clhs11);
const double clhs20 =             1.0/(K_darcy + clhs12/h + clhs17*dyn_tau + mu*stab_c1/pow(h, 2));
const double clhs21 =             clhs19*clhs20;
const double clhs22 =             N[0]*clhs17;
const double clhs23 =             K_darcy*clhs18;
const double clhs24 =             0.5*clhs19 + clhs23;
const double clhs25 =             clhs22 + clhs24;
const double clhs26 =             1.0*clhs25;
const double clhs27 =             K_darcy*clhs20;
const double clhs28 =             N[0]*clhs27;
const double clhs29 =             clhs15*clhs16 + clhs15*clhs17 + clhs18*clhs19 + clhs21*clhs26 - clhs26*clhs28;
const double clhs30 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
const double clhs31 =             C(1,3)*DN(0,1);
const double clhs32 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs31;
const double clhs33 =             C(3,5)*DN(0,0);
const double clhs34 =             C(4,5)*DN(0,2);
const double clhs35 =             C(1,5)*DN(0,1) + clhs33 + clhs34;
const double clhs36 =             DN(0,0)*clhs13;
const double clhs37 =             DN(0,1)*clhs36;
const double clhs38 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs39 =             C(3,4)*DN(0,1);
const double clhs40 =             C(2,3)*DN(0,2) + clhs33 + clhs39;
const double clhs41 =             C(2,5)*DN(0,2);
const double clhs42 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs41;
const double clhs43 =             DN(0,2)*clhs36;
const double clhs44 =             -N[0] + clhs21 - clhs28;
const double clhs45 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs46 =             C(0,3)*DN(1,0);
const double clhs47 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs46;
const double clhs48 =             C(0,5)*DN(1,0);
const double clhs49 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs48;
const double clhs50 =             DN(0,0)*DN(1,0);
const double clhs51 =             N[1]*clhs22 + N[1]*clhs23;
const double clhs52 =             clhs14*clhs50 + clhs51;
const double clhs53 =             DN(1,0)*clhs9 + DN(1,1)*clhs10 + DN(1,2)*clhs11;
const double clhs54 =             clhs18*rho;
const double clhs55 =             N[1]*clhs17;
const double clhs56 =             0.5*N[1];
const double clhs57 =             K_darcy*clhs56;
const double clhs58 =             0.5*rho;
const double clhs59 =             clhs53*clhs58 + clhs57;
const double clhs60 =             clhs55 + clhs59;
const double clhs61 =             1.0*clhs60;
const double clhs62 =             clhs21*clhs61 - clhs28*clhs61 + clhs53*clhs54;
const double clhs63 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs46;
const double clhs64 =             C(1,3)*DN(1,1);
const double clhs65 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs64;
const double clhs66 =             C(3,5)*DN(1,0);
const double clhs67 =             C(4,5)*DN(1,2);
const double clhs68 =             C(1,5)*DN(1,1) + clhs66 + clhs67;
const double clhs69 =             DN(1,1)*clhs36;
const double clhs70 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs48;
const double clhs71 =             C(3,4)*DN(1,1);
const double clhs72 =             C(2,3)*DN(1,2) + clhs66 + clhs71;
const double clhs73 =             C(2,5)*DN(1,2);
const double clhs74 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs73;
const double clhs75 =             DN(1,2)*clhs36;
const double clhs76 =             DN(0,0)*N[1];
const double clhs77 =             DN(1,0)*N[0];
const double clhs78 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs79 =             C(0,3)*DN(2,0);
const double clhs80 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs79;
const double clhs81 =             C(0,5)*DN(2,0);
const double clhs82 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs81;
const double clhs83 =             DN(0,0)*DN(2,0);
const double clhs84 =             N[2]*clhs22 + N[2]*clhs23;
const double clhs85 =             clhs14*clhs83 + clhs84;
const double clhs86 =             DN(2,0)*clhs9 + DN(2,1)*clhs10 + DN(2,2)*clhs11;
const double clhs87 =             N[2]*clhs17;
const double clhs88 =             0.5*N[2];
const double clhs89 =             K_darcy*clhs88;
const double clhs90 =             clhs58*clhs86 + clhs89;
const double clhs91 =             clhs87 + clhs90;
const double clhs92 =             1.0*clhs91;
const double clhs93 =             clhs21*clhs92 - clhs28*clhs92 + clhs54*clhs86;
const double clhs94 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs79;
const double clhs95 =             C(1,3)*DN(2,1);
const double clhs96 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs95;
const double clhs97 =             C(3,5)*DN(2,0);
const double clhs98 =             C(4,5)*DN(2,2);
const double clhs99 =             C(1,5)*DN(2,1) + clhs97 + clhs98;
const double clhs100 =             DN(2,1)*clhs36;
const double clhs101 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs81;
const double clhs102 =             C(3,4)*DN(2,1);
const double clhs103 =             C(2,3)*DN(2,2) + clhs102 + clhs97;
const double clhs104 =             C(2,5)*DN(2,2);
const double clhs105 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs104;
const double clhs106 =             DN(2,2)*clhs36;
const double clhs107 =             DN(0,0)*N[2];
const double clhs108 =             DN(2,0)*N[0];
const double clhs109 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs110 =             C(0,3)*DN(3,0);
const double clhs111 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs110;
const double clhs112 =             C(0,5)*DN(3,0);
const double clhs113 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs112;
const double clhs114 =             DN(0,0)*DN(3,0);
const double clhs115 =             N[3]*clhs22 + N[3]*clhs23;
const double clhs116 =             clhs114*clhs14 + clhs115;
const double clhs117 =             DN(3,0)*clhs9 + DN(3,1)*clhs10 + DN(3,2)*clhs11;
const double clhs118 =             N[3]*clhs17;
const double clhs119 =             0.5*N[3];
const double clhs120 =             K_darcy*clhs119 + clhs117*clhs58;
const double clhs121 =             clhs118 + clhs120;
const double clhs122 =             1.0*clhs121;
const double clhs123 =             clhs117*clhs54 + clhs122*clhs21 - clhs122*clhs28;
const double clhs124 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs110;
const double clhs125 =             C(1,3)*DN(3,1);
const double clhs126 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs125;
const double clhs127 =             C(3,5)*DN(3,0);
const double clhs128 =             C(4,5)*DN(3,2);
const double clhs129 =             C(1,5)*DN(3,1) + clhs127 + clhs128;
const double clhs130 =             DN(3,1)*clhs36;
const double clhs131 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs112;
const double clhs132 =             C(3,4)*DN(3,1);
const double clhs133 =             C(2,3)*DN(3,2) + clhs127 + clhs132;
const double clhs134 =             C(2,5)*DN(3,2);
const double clhs135 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs134;
const double clhs136 =             DN(3,2)*clhs36;
const double clhs137 =             DN(0,0)*N[3];
const double clhs138 =             DN(3,0)*N[0];
const double clhs139 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs31;
const double clhs140 =             C(0,4)*DN(0,0) + clhs34 + clhs39;
const double clhs141 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs142 =             C(1,4)*DN(0,1);
const double clhs143 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs142;
const double clhs144 =             pow(DN(0,1), 2);
const double clhs145 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs142;
const double clhs146 =             C(2,4)*DN(0,2);
const double clhs147 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs146;
const double clhs148 =             DN(0,1)*clhs13;
const double clhs149 =             DN(0,2)*clhs148;
const double clhs150 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs64;
const double clhs151 =             C(0,4)*DN(1,0) + clhs67 + clhs71;
const double clhs152 =             DN(1,0)*clhs148;
const double clhs153 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs154 =             C(1,4)*DN(1,1);
const double clhs155 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs154;
const double clhs156 =             DN(0,1)*DN(1,1);
const double clhs157 =             clhs14*clhs156;
const double clhs158 =             clhs51 + clhs62;
const double clhs159 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs154;
const double clhs160 =             C(2,4)*DN(1,2);
const double clhs161 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs160;
const double clhs162 =             DN(1,2)*clhs148;
const double clhs163 =             DN(0,1)*N[1];
const double clhs164 =             DN(1,1)*N[0];
const double clhs165 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs95;
const double clhs166 =             C(0,4)*DN(2,0) + clhs102 + clhs98;
const double clhs167 =             DN(2,0)*clhs148;
const double clhs168 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs169 =             C(1,4)*DN(2,1);
const double clhs170 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs169;
const double clhs171 =             DN(0,1)*DN(2,1);
const double clhs172 =             clhs14*clhs171;
const double clhs173 =             clhs84 + clhs93;
const double clhs174 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs169;
const double clhs175 =             C(2,4)*DN(2,2);
const double clhs176 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs175;
const double clhs177 =             DN(2,2)*clhs148;
const double clhs178 =             DN(0,1)*N[2];
const double clhs179 =             DN(2,1)*N[0];
const double clhs180 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs125;
const double clhs181 =             C(0,4)*DN(3,0) + clhs128 + clhs132;
const double clhs182 =             DN(3,0)*clhs148;
const double clhs183 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs184 =             C(1,4)*DN(3,1);
const double clhs185 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs184;
const double clhs186 =             DN(0,1)*DN(3,1);
const double clhs187 =             clhs14*clhs186;
const double clhs188 =             clhs115 + clhs123;
const double clhs189 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs184;
const double clhs190 =             C(2,4)*DN(3,2);
const double clhs191 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs190;
const double clhs192 =             DN(3,2)*clhs148;
const double clhs193 =             DN(0,1)*N[3];
const double clhs194 =             DN(3,1)*N[0];
const double clhs195 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs41;
const double clhs196 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs146;
const double clhs197 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs198 =             pow(DN(0,2), 2);
const double clhs199 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs73;
const double clhs200 =             DN(0,2)*clhs13;
const double clhs201 =             DN(1,0)*clhs200;
const double clhs202 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs160;
const double clhs203 =             DN(1,1)*clhs200;
const double clhs204 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs205 =             DN(0,2)*DN(1,2);
const double clhs206 =             clhs14*clhs205;
const double clhs207 =             DN(0,2)*N[1];
const double clhs208 =             DN(1,2)*N[0];
const double clhs209 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs104;
const double clhs210 =             DN(2,0)*clhs200;
const double clhs211 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs175;
const double clhs212 =             DN(2,1)*clhs200;
const double clhs213 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs214 =             DN(0,2)*DN(2,2);
const double clhs215 =             clhs14*clhs214;
const double clhs216 =             DN(0,2)*N[2];
const double clhs217 =             DN(2,2)*N[0];
const double clhs218 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs134;
const double clhs219 =             DN(3,0)*clhs200;
const double clhs220 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs190;
const double clhs221 =             DN(3,1)*clhs200;
const double clhs222 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs223 =             DN(0,2)*DN(3,2);
const double clhs224 =             clhs14*clhs223;
const double clhs225 =             DN(0,2)*N[3];
const double clhs226 =             DN(3,2)*N[0];
const double clhs227 =             clhs18 + clhs20*(1.0*clhs22 + clhs24);
const double clhs228 =             0.5*clhs20;
const double clhs229 =             1.0*clhs20;
const double clhs230 =             clhs229*clhs60;
const double clhs231 =             clhs228*(clhs156 + clhs205 + clhs50);
const double clhs232 =             clhs229*clhs91;
const double clhs233 =             clhs228*(clhs171 + clhs214 + clhs83);
const double clhs234 =             clhs121*clhs229;
const double clhs235 =             clhs228*(clhs114 + clhs186 + clhs223);
const double clhs236 =             0.5*DN(1,0);
const double clhs237 =             0.5*DN(1,1);
const double clhs238 =             0.5*DN(1,2);
const double clhs239 =             clhs20*rho;
const double clhs240 =             clhs239*clhs53;
const double clhs241 =             N[1]*clhs27;
const double clhs242 =             clhs19*clhs56 + clhs240*clhs26 - clhs241*clhs26;
const double clhs243 =             pow(DN(1,0), 2);
const double clhs244 =             pow(N[1], 2);
const double clhs245 =             clhs56*rho;
const double clhs246 =             clhs16*clhs244 + clhs17*clhs244 + clhs240*clhs61 - clhs241*clhs61 + clhs245*clhs53;
const double clhs247 =             DN(1,0)*clhs13;
const double clhs248 =             DN(1,1)*clhs247;
const double clhs249 =             DN(1,2)*clhs247;
const double clhs250 =             -N[1] + clhs240 - clhs241;
const double clhs251 =             DN(1,0)*DN(2,0);
const double clhs252 =             N[2]*clhs55 + N[2]*clhs57;
const double clhs253 =             clhs14*clhs251 + clhs252;
const double clhs254 =             clhs240*clhs92 - clhs241*clhs92 + clhs245*clhs86;
const double clhs255 =             DN(2,1)*clhs247;
const double clhs256 =             DN(2,2)*clhs247;
const double clhs257 =             DN(1,0)*N[2];
const double clhs258 =             DN(2,0)*N[1];
const double clhs259 =             DN(1,0)*DN(3,0);
const double clhs260 =             N[3]*clhs55 + N[3]*clhs57;
const double clhs261 =             clhs14*clhs259 + clhs260;
const double clhs262 =             clhs117*clhs245 + clhs122*clhs240 - clhs122*clhs241;
const double clhs263 =             DN(3,1)*clhs247;
const double clhs264 =             DN(3,2)*clhs247;
const double clhs265 =             DN(1,0)*N[3];
const double clhs266 =             DN(3,0)*N[1];
const double clhs267 =             clhs242 + clhs51;
const double clhs268 =             pow(DN(1,1), 2);
const double clhs269 =             DN(1,1)*clhs13;
const double clhs270 =             DN(1,2)*clhs269;
const double clhs271 =             DN(2,0)*clhs269;
const double clhs272 =             DN(1,1)*DN(2,1);
const double clhs273 =             clhs14*clhs272;
const double clhs274 =             clhs252 + clhs254;
const double clhs275 =             DN(2,2)*clhs269;
const double clhs276 =             DN(1,1)*N[2];
const double clhs277 =             DN(2,1)*N[1];
const double clhs278 =             DN(3,0)*clhs269;
const double clhs279 =             DN(1,1)*DN(3,1);
const double clhs280 =             clhs14*clhs279;
const double clhs281 =             clhs260 + clhs262;
const double clhs282 =             DN(3,2)*clhs269;
const double clhs283 =             DN(1,1)*N[3];
const double clhs284 =             DN(3,1)*N[1];
const double clhs285 =             pow(DN(1,2), 2);
const double clhs286 =             DN(1,2)*clhs13;
const double clhs287 =             DN(2,0)*clhs286;
const double clhs288 =             DN(2,1)*clhs286;
const double clhs289 =             DN(1,2)*DN(2,2);
const double clhs290 =             clhs14*clhs289;
const double clhs291 =             DN(1,2)*N[2];
const double clhs292 =             DN(2,2)*N[1];
const double clhs293 =             DN(3,0)*clhs286;
const double clhs294 =             DN(3,1)*clhs286;
const double clhs295 =             DN(1,2)*DN(3,2);
const double clhs296 =             clhs14*clhs295;
const double clhs297 =             DN(1,2)*N[3];
const double clhs298 =             DN(3,2)*N[1];
const double clhs299 =             clhs229*clhs25;
const double clhs300 =             clhs20*(1.0*clhs55 + clhs59) + clhs56;
const double clhs301 =             clhs228*(clhs251 + clhs272 + clhs289);
const double clhs302 =             clhs228*(clhs259 + clhs279 + clhs295);
const double clhs303 =             0.5*DN(2,0);
const double clhs304 =             0.5*DN(2,1);
const double clhs305 =             0.5*DN(2,2);
const double clhs306 =             clhs239*clhs86;
const double clhs307 =             N[2]*clhs27;
const double clhs308 =             clhs19*clhs88 + clhs26*clhs306 - clhs26*clhs307;
const double clhs309 =             clhs88*rho;
const double clhs310 =             clhs306*clhs61 - clhs307*clhs61 + clhs309*clhs53;
const double clhs311 =             pow(DN(2,0), 2);
const double clhs312 =             pow(N[2], 2);
const double clhs313 =             clhs16*clhs312 + clhs17*clhs312 + clhs306*clhs92 - clhs307*clhs92 + clhs309*clhs86;
const double clhs314 =             DN(2,0)*clhs13;
const double clhs315 =             DN(2,1)*clhs314;
const double clhs316 =             DN(2,2)*clhs314;
const double clhs317 =             -N[2] + clhs306 - clhs307;
const double clhs318 =             DN(2,0)*DN(3,0);
const double clhs319 =             N[3]*clhs87 + N[3]*clhs89;
const double clhs320 =             clhs14*clhs318 + clhs319;
const double clhs321 =             clhs117*clhs309 + clhs122*clhs306 - clhs122*clhs307;
const double clhs322 =             DN(3,1)*clhs314;
const double clhs323 =             DN(3,2)*clhs314;
const double clhs324 =             DN(2,0)*N[3];
const double clhs325 =             DN(3,0)*N[2];
const double clhs326 =             clhs308 + clhs84;
const double clhs327 =             clhs252 + clhs310;
const double clhs328 =             pow(DN(2,1), 2);
const double clhs329 =             DN(2,1)*clhs13;
const double clhs330 =             DN(2,2)*clhs329;
const double clhs331 =             DN(3,0)*clhs329;
const double clhs332 =             DN(2,1)*DN(3,1);
const double clhs333 =             clhs14*clhs332;
const double clhs334 =             clhs319 + clhs321;
const double clhs335 =             DN(3,2)*clhs329;
const double clhs336 =             DN(2,1)*N[3];
const double clhs337 =             DN(3,1)*N[2];
const double clhs338 =             pow(DN(2,2), 2);
const double clhs339 =             DN(2,2)*clhs13;
const double clhs340 =             DN(3,0)*clhs339;
const double clhs341 =             DN(3,1)*clhs339;
const double clhs342 =             DN(2,2)*DN(3,2);
const double clhs343 =             clhs14*clhs342;
const double clhs344 =             DN(2,2)*N[3];
const double clhs345 =             DN(3,2)*N[2];
const double clhs346 =             clhs20*(1.0*clhs87 + clhs90) + clhs88;
const double clhs347 =             clhs228*(clhs318 + clhs332 + clhs342);
const double clhs348 =             0.5*DN(3,0);
const double clhs349 =             0.5*DN(3,1);
const double clhs350 =             0.5*DN(3,2);
const double clhs351 =             clhs117*clhs239;
const double clhs352 =             N[3]*clhs27;
const double clhs353 =             clhs119*clhs19 + clhs26*clhs351 - clhs26*clhs352;
const double clhs354 =             clhs119*rho;
const double clhs355 =             clhs351*clhs61 - clhs352*clhs61 + clhs354*clhs53;
const double clhs356 =             clhs351*clhs92 - clhs352*clhs92 + clhs354*clhs86;
const double clhs357 =             pow(DN(3,0), 2);
const double clhs358 =             pow(N[3], 2);
const double clhs359 =             clhs117*clhs354 + clhs122*clhs351 - clhs122*clhs352 + clhs16*clhs358 + clhs17*clhs358;
const double clhs360 =             DN(3,0)*clhs13;
const double clhs361 =             DN(3,1)*clhs360;
const double clhs362 =             DN(3,2)*clhs360;
const double clhs363 =             -N[3] + clhs351 - clhs352;
const double clhs364 =             clhs115 + clhs353;
const double clhs365 =             clhs260 + clhs355;
const double clhs366 =             clhs319 + clhs356;
const double clhs367 =             pow(DN(3,1), 2);
const double clhs368 =             DN(3,1)*DN(3,2)*clhs13;
const double clhs369 =             pow(DN(3,2), 2);
const double clhs370 =             clhs119 + clhs20*(1.0*clhs118 + clhs120);
            lhs(0,0)=clhs0*clhs1 + clhs14*clhs8 + clhs29 + clhs3*clhs4 + clhs6*clhs7;
            lhs(0,1)=0.5*DN(0,0)*clhs30 + 0.5*DN(0,1)*clhs32 + 0.5*DN(0,2)*clhs35 + 0.5*clhs37;
            lhs(0,2)=0.5*DN(0,0)*clhs38 + 0.5*DN(0,1)*clhs40 + 0.5*DN(0,2)*clhs42 + 0.5*clhs43;
            lhs(0,3)=clhs1*clhs44;
            lhs(0,4)=clhs1*clhs45 + clhs4*clhs47 + clhs49*clhs7 + clhs52 + clhs62;
            lhs(0,5)=0.5*DN(0,0)*clhs63 + 0.5*DN(0,1)*clhs65 + 0.5*DN(0,2)*clhs68 + 0.5*clhs69;
            lhs(0,6)=0.5*DN(0,0)*clhs70 + 0.5*DN(0,1)*clhs72 + 0.5*DN(0,2)*clhs74 + 0.5*clhs75;
            lhs(0,7)=0.5*DN(1,0)*clhs21 - 0.5*clhs27*clhs77 - 0.5*clhs76;
            lhs(0,8)=clhs1*clhs78 + clhs4*clhs80 + clhs7*clhs82 + clhs85 + clhs93;
            lhs(0,9)=0.5*DN(0,0)*clhs94 + 0.5*DN(0,1)*clhs96 + 0.5*DN(0,2)*clhs99 + 0.5*clhs100;
            lhs(0,10)=0.5*DN(0,0)*clhs101 + 0.5*DN(0,1)*clhs103 + 0.5*DN(0,2)*clhs105 + 0.5*clhs106;
            lhs(0,11)=0.5*DN(2,0)*clhs21 - 0.5*clhs107 - 0.5*clhs108*clhs27;
            lhs(0,12)=clhs1*clhs109 + clhs111*clhs4 + clhs113*clhs7 + clhs116 + clhs123;
            lhs(0,13)=0.5*DN(0,0)*clhs124 + 0.5*DN(0,1)*clhs126 + 0.5*DN(0,2)*clhs129 + 0.5*clhs130;
            lhs(0,14)=0.5*DN(0,0)*clhs131 + 0.5*DN(0,1)*clhs133 + 0.5*DN(0,2)*clhs135 + 0.5*clhs136;
            lhs(0,15)=0.5*DN(3,0)*clhs21 - 0.5*clhs137 - 0.5*clhs138*clhs27;
            lhs(1,0)=0.5*DN(0,0)*clhs3 + 0.5*DN(0,1)*clhs139 + 0.5*DN(0,2)*clhs140 + 0.5*clhs37;
            lhs(1,1)=clhs1*clhs32 + clhs14*clhs144 + clhs141*clhs4 + clhs143*clhs7 + clhs29;
            lhs(1,2)=0.5*DN(0,0)*clhs40 + 0.5*DN(0,1)*clhs145 + 0.5*DN(0,2)*clhs147 + 0.5*clhs149;
            lhs(1,3)=clhs4*clhs44;
            lhs(1,4)=0.5*DN(0,0)*clhs47 + 0.5*DN(0,1)*clhs150 + 0.5*DN(0,2)*clhs151 + 0.5*clhs152;
            lhs(1,5)=clhs1*clhs65 + clhs153*clhs4 + clhs155*clhs7 + clhs157 + clhs158;
            lhs(1,6)=0.5*DN(0,0)*clhs72 + 0.5*DN(0,1)*clhs159 + 0.5*DN(0,2)*clhs161 + 0.5*clhs162;
            lhs(1,7)=0.5*DN(1,1)*clhs21 - 0.5*clhs163 - 0.5*clhs164*clhs27;
            lhs(1,8)=0.5*DN(0,0)*clhs80 + 0.5*DN(0,1)*clhs165 + 0.5*DN(0,2)*clhs166 + 0.5*clhs167;
            lhs(1,9)=clhs1*clhs96 + clhs168*clhs4 + clhs170*clhs7 + clhs172 + clhs173;
            lhs(1,10)=0.5*DN(0,0)*clhs103 + 0.5*DN(0,1)*clhs174 + 0.5*DN(0,2)*clhs176 + 0.5*clhs177;
            lhs(1,11)=0.5*DN(2,1)*clhs21 - 0.5*clhs178 - 0.5*clhs179*clhs27;
            lhs(1,12)=0.5*DN(0,0)*clhs111 + 0.5*DN(0,1)*clhs180 + 0.5*DN(0,2)*clhs181 + 0.5*clhs182;
            lhs(1,13)=clhs1*clhs126 + clhs183*clhs4 + clhs185*clhs7 + clhs187 + clhs188;
            lhs(1,14)=0.5*DN(0,0)*clhs133 + 0.5*DN(0,1)*clhs189 + 0.5*DN(0,2)*clhs191 + 0.5*clhs192;
            lhs(1,15)=0.5*DN(3,1)*clhs21 - 0.5*clhs193 - 0.5*clhs194*clhs27;
            lhs(2,0)=0.5*DN(0,0)*clhs6 + 0.5*DN(0,1)*clhs140 + 0.5*DN(0,2)*clhs195 + 0.5*clhs43;
            lhs(2,1)=0.5*DN(0,0)*clhs35 + 0.5*DN(0,1)*clhs143 + 0.5*DN(0,2)*clhs196 + 0.5*clhs149;
            lhs(2,2)=clhs1*clhs42 + clhs14*clhs198 + clhs147*clhs4 + clhs197*clhs7 + clhs29;
            lhs(2,3)=clhs44*clhs7;
            lhs(2,4)=0.5*DN(0,0)*clhs49 + 0.5*DN(0,1)*clhs151 + 0.5*DN(0,2)*clhs199 + 0.5*clhs201;
            lhs(2,5)=0.5*DN(0,0)*clhs68 + 0.5*DN(0,1)*clhs155 + 0.5*DN(0,2)*clhs202 + 0.5*clhs203;
            lhs(2,6)=clhs1*clhs74 + clhs158 + clhs161*clhs4 + clhs204*clhs7 + clhs206;
            lhs(2,7)=0.5*DN(1,2)*clhs21 - 0.5*clhs207 - 0.5*clhs208*clhs27;
            lhs(2,8)=0.5*DN(0,0)*clhs82 + 0.5*DN(0,1)*clhs166 + 0.5*DN(0,2)*clhs209 + 0.5*clhs210;
            lhs(2,9)=0.5*DN(0,0)*clhs99 + 0.5*DN(0,1)*clhs170 + 0.5*DN(0,2)*clhs211 + 0.5*clhs212;
            lhs(2,10)=clhs1*clhs105 + clhs173 + clhs176*clhs4 + clhs213*clhs7 + clhs215;
            lhs(2,11)=0.5*DN(2,2)*clhs21 - 0.5*clhs216 - 0.5*clhs217*clhs27;
            lhs(2,12)=0.5*DN(0,0)*clhs113 + 0.5*DN(0,1)*clhs181 + 0.5*DN(0,2)*clhs218 + 0.5*clhs219;
            lhs(2,13)=0.5*DN(0,0)*clhs129 + 0.5*DN(0,1)*clhs185 + 0.5*DN(0,2)*clhs220 + 0.5*clhs221;
            lhs(2,14)=clhs1*clhs135 + clhs188 + clhs191*clhs4 + clhs222*clhs7 + clhs224;
            lhs(2,15)=0.5*DN(3,2)*clhs21 - 0.5*clhs225 - 0.5*clhs226*clhs27;
            lhs(3,0)=DN(0,0)*clhs227;
            lhs(3,1)=DN(0,1)*clhs227;
            lhs(3,2)=DN(0,2)*clhs227;
            lhs(3,3)=clhs228*(clhs144 + clhs198 + clhs8);
            lhs(3,4)=DN(0,0)*clhs230 + DN(1,0)*clhs18;
            lhs(3,5)=DN(0,1)*clhs230 + DN(1,1)*clhs18;
            lhs(3,6)=DN(0,2)*clhs230 + DN(1,2)*clhs18;
            lhs(3,7)=clhs231;
            lhs(3,8)=DN(0,0)*clhs232 + DN(2,0)*clhs18;
            lhs(3,9)=DN(0,1)*clhs232 + DN(2,1)*clhs18;
            lhs(3,10)=DN(0,2)*clhs232 + DN(2,2)*clhs18;
            lhs(3,11)=clhs233;
            lhs(3,12)=DN(0,0)*clhs234 + DN(3,0)*clhs18;
            lhs(3,13)=DN(0,1)*clhs234 + DN(3,1)*clhs18;
            lhs(3,14)=DN(0,2)*clhs234 + DN(3,2)*clhs18;
            lhs(3,15)=clhs235;
            lhs(4,0)=clhs0*clhs236 + clhs237*clhs3 + clhs238*clhs6 + clhs242 + clhs52;
            lhs(4,1)=0.5*DN(1,0)*clhs30 + 0.5*DN(1,1)*clhs32 + 0.5*DN(1,2)*clhs35 + 0.5*clhs152;
            lhs(4,2)=0.5*DN(1,0)*clhs38 + 0.5*DN(1,1)*clhs40 + 0.5*DN(1,2)*clhs42 + 0.5*clhs201;
            lhs(4,3)=0.5*DN(0,0)*clhs240 - 0.5*clhs27*clhs76 - 0.5*clhs77;
            lhs(4,4)=clhs14*clhs243 + clhs236*clhs45 + clhs237*clhs47 + clhs238*clhs49 + clhs246;
            lhs(4,5)=0.5*DN(1,0)*clhs63 + 0.5*DN(1,1)*clhs65 + 0.5*DN(1,2)*clhs68 + 0.5*clhs248;
            lhs(4,6)=0.5*DN(1,0)*clhs70 + 0.5*DN(1,1)*clhs72 + 0.5*DN(1,2)*clhs74 + 0.5*clhs249;
            lhs(4,7)=clhs236*clhs250;
            lhs(4,8)=clhs236*clhs78 + clhs237*clhs80 + clhs238*clhs82 + clhs253 + clhs254;
            lhs(4,9)=0.5*DN(1,0)*clhs94 + 0.5*DN(1,1)*clhs96 + 0.5*DN(1,2)*clhs99 + 0.5*clhs255;
            lhs(4,10)=0.5*DN(1,0)*clhs101 + 0.5*DN(1,1)*clhs103 + 0.5*DN(1,2)*clhs105 + 0.5*clhs256;
            lhs(4,11)=0.5*DN(2,0)*clhs240 - 0.5*clhs257 - 0.5*clhs258*clhs27;
            lhs(4,12)=clhs109*clhs236 + clhs111*clhs237 + clhs113*clhs238 + clhs261 + clhs262;
            lhs(4,13)=0.5*DN(1,0)*clhs124 + 0.5*DN(1,1)*clhs126 + 0.5*DN(1,2)*clhs129 + 0.5*clhs263;
            lhs(4,14)=0.5*DN(1,0)*clhs131 + 0.5*DN(1,1)*clhs133 + 0.5*DN(1,2)*clhs135 + 0.5*clhs264;
            lhs(4,15)=0.5*DN(3,0)*clhs240 - 0.5*clhs265 - 0.5*clhs266*clhs27;
            lhs(5,0)=0.5*DN(1,0)*clhs3 + 0.5*DN(1,1)*clhs139 + 0.5*DN(1,2)*clhs140 + 0.5*clhs69;
            lhs(5,1)=clhs141*clhs237 + clhs143*clhs238 + clhs157 + clhs236*clhs32 + clhs267;
            lhs(5,2)=0.5*DN(1,0)*clhs40 + 0.5*DN(1,1)*clhs145 + 0.5*DN(1,2)*clhs147 + 0.5*clhs203;
            lhs(5,3)=0.5*DN(0,1)*clhs240 - 0.5*clhs163*clhs27 - 0.5*clhs164;
            lhs(5,4)=0.5*DN(1,0)*clhs47 + 0.5*DN(1,1)*clhs150 + 0.5*DN(1,2)*clhs151 + 0.5*clhs248;
            lhs(5,5)=clhs14*clhs268 + clhs153*clhs237 + clhs155*clhs238 + clhs236*clhs65 + clhs246;
            lhs(5,6)=0.5*DN(1,0)*clhs72 + 0.5*DN(1,1)*clhs159 + 0.5*DN(1,2)*clhs161 + 0.5*clhs270;
            lhs(5,7)=clhs237*clhs250;
            lhs(5,8)=0.5*DN(1,0)*clhs80 + 0.5*DN(1,1)*clhs165 + 0.5*DN(1,2)*clhs166 + 0.5*clhs271;
            lhs(5,9)=clhs168*clhs237 + clhs170*clhs238 + clhs236*clhs96 + clhs273 + clhs274;
            lhs(5,10)=0.5*DN(1,0)*clhs103 + 0.5*DN(1,1)*clhs174 + 0.5*DN(1,2)*clhs176 + 0.5*clhs275;
            lhs(5,11)=0.5*DN(2,1)*clhs240 - 0.5*clhs27*clhs277 - 0.5*clhs276;
            lhs(5,12)=0.5*DN(1,0)*clhs111 + 0.5*DN(1,1)*clhs180 + 0.5*DN(1,2)*clhs181 + 0.5*clhs278;
            lhs(5,13)=clhs126*clhs236 + clhs183*clhs237 + clhs185*clhs238 + clhs280 + clhs281;
            lhs(5,14)=0.5*DN(1,0)*clhs133 + 0.5*DN(1,1)*clhs189 + 0.5*DN(1,2)*clhs191 + 0.5*clhs282;
            lhs(5,15)=0.5*DN(3,1)*clhs240 - 0.5*clhs27*clhs284 - 0.5*clhs283;
            lhs(6,0)=0.5*DN(1,0)*clhs6 + 0.5*DN(1,1)*clhs140 + 0.5*DN(1,2)*clhs195 + 0.5*clhs75;
            lhs(6,1)=0.5*DN(1,0)*clhs35 + 0.5*DN(1,1)*clhs143 + 0.5*DN(1,2)*clhs196 + 0.5*clhs162;
            lhs(6,2)=clhs147*clhs237 + clhs197*clhs238 + clhs206 + clhs236*clhs42 + clhs267;
            lhs(6,3)=0.5*DN(0,2)*clhs240 - 0.5*clhs207*clhs27 - 0.5*clhs208;
            lhs(6,4)=0.5*DN(1,0)*clhs49 + 0.5*DN(1,1)*clhs151 + 0.5*DN(1,2)*clhs199 + 0.5*clhs249;
            lhs(6,5)=0.5*DN(1,0)*clhs68 + 0.5*DN(1,1)*clhs155 + 0.5*DN(1,2)*clhs202 + 0.5*clhs270;
            lhs(6,6)=clhs14*clhs285 + clhs161*clhs237 + clhs204*clhs238 + clhs236*clhs74 + clhs246;
            lhs(6,7)=clhs238*clhs250;
            lhs(6,8)=0.5*DN(1,0)*clhs82 + 0.5*DN(1,1)*clhs166 + 0.5*DN(1,2)*clhs209 + 0.5*clhs287;
            lhs(6,9)=0.5*DN(1,0)*clhs99 + 0.5*DN(1,1)*clhs170 + 0.5*DN(1,2)*clhs211 + 0.5*clhs288;
            lhs(6,10)=clhs105*clhs236 + clhs176*clhs237 + clhs213*clhs238 + clhs274 + clhs290;
            lhs(6,11)=0.5*DN(2,2)*clhs240 - 0.5*clhs27*clhs292 - 0.5*clhs291;
            lhs(6,12)=0.5*DN(1,0)*clhs113 + 0.5*DN(1,1)*clhs181 + 0.5*DN(1,2)*clhs218 + 0.5*clhs293;
            lhs(6,13)=0.5*DN(1,0)*clhs129 + 0.5*DN(1,1)*clhs185 + 0.5*DN(1,2)*clhs220 + 0.5*clhs294;
            lhs(6,14)=clhs135*clhs236 + clhs191*clhs237 + clhs222*clhs238 + clhs281 + clhs296;
            lhs(6,15)=0.5*DN(3,2)*clhs240 - 0.5*clhs27*clhs298 - 0.5*clhs297;
            lhs(7,0)=DN(1,0)*clhs299 + 0.5*clhs76;
            lhs(7,1)=DN(1,1)*clhs299 + 0.5*clhs163;
            lhs(7,2)=DN(1,2)*clhs299 + 0.5*clhs207;
            lhs(7,3)=clhs231;
            lhs(7,4)=DN(1,0)*clhs300;
            lhs(7,5)=DN(1,1)*clhs300;
            lhs(7,6)=DN(1,2)*clhs300;
            lhs(7,7)=clhs228*(clhs243 + clhs268 + clhs285);
            lhs(7,8)=DN(1,0)*clhs232 + DN(2,0)*clhs56;
            lhs(7,9)=DN(1,1)*clhs232 + DN(2,1)*clhs56;
            lhs(7,10)=DN(1,2)*clhs232 + DN(2,2)*clhs56;
            lhs(7,11)=clhs301;
            lhs(7,12)=DN(1,0)*clhs234 + DN(3,0)*clhs56;
            lhs(7,13)=DN(1,1)*clhs234 + DN(3,1)*clhs56;
            lhs(7,14)=DN(1,2)*clhs234 + DN(3,2)*clhs56;
            lhs(7,15)=clhs302;
            lhs(8,0)=clhs0*clhs303 + clhs3*clhs304 + clhs305*clhs6 + clhs308 + clhs85;
            lhs(8,1)=0.5*DN(2,0)*clhs30 + 0.5*DN(2,1)*clhs32 + 0.5*DN(2,2)*clhs35 + 0.5*clhs167;
            lhs(8,2)=0.5*DN(2,0)*clhs38 + 0.5*DN(2,1)*clhs40 + 0.5*DN(2,2)*clhs42 + 0.5*clhs210;
            lhs(8,3)=0.5*DN(0,0)*clhs306 - 0.5*clhs107*clhs27 - 0.5*clhs108;
            lhs(8,4)=clhs253 + clhs303*clhs45 + clhs304*clhs47 + clhs305*clhs49 + clhs310;
            lhs(8,5)=0.5*DN(2,0)*clhs63 + 0.5*DN(2,1)*clhs65 + 0.5*DN(2,2)*clhs68 + 0.5*clhs271;
            lhs(8,6)=0.5*DN(2,0)*clhs70 + 0.5*DN(2,1)*clhs72 + 0.5*DN(2,2)*clhs74 + 0.5*clhs287;
            lhs(8,7)=0.5*DN(1,0)*clhs306 - 0.5*clhs257*clhs27 - 0.5*clhs258;
            lhs(8,8)=clhs14*clhs311 + clhs303*clhs78 + clhs304*clhs80 + clhs305*clhs82 + clhs313;
            lhs(8,9)=0.5*DN(2,0)*clhs94 + 0.5*DN(2,1)*clhs96 + 0.5*DN(2,2)*clhs99 + 0.5*clhs315;
            lhs(8,10)=0.5*DN(2,0)*clhs101 + 0.5*DN(2,1)*clhs103 + 0.5*DN(2,2)*clhs105 + 0.5*clhs316;
            lhs(8,11)=clhs303*clhs317;
            lhs(8,12)=clhs109*clhs303 + clhs111*clhs304 + clhs113*clhs305 + clhs320 + clhs321;
            lhs(8,13)=0.5*DN(2,0)*clhs124 + 0.5*DN(2,1)*clhs126 + 0.5*DN(2,2)*clhs129 + 0.5*clhs322;
            lhs(8,14)=0.5*DN(2,0)*clhs131 + 0.5*DN(2,1)*clhs133 + 0.5*DN(2,2)*clhs135 + 0.5*clhs323;
            lhs(8,15)=0.5*DN(3,0)*clhs306 - 0.5*clhs27*clhs325 - 0.5*clhs324;
            lhs(9,0)=0.5*DN(2,0)*clhs3 + 0.5*DN(2,1)*clhs139 + 0.5*DN(2,2)*clhs140 + 0.5*clhs100;
            lhs(9,1)=clhs141*clhs304 + clhs143*clhs305 + clhs172 + clhs303*clhs32 + clhs326;
            lhs(9,2)=0.5*DN(2,0)*clhs40 + 0.5*DN(2,1)*clhs145 + 0.5*DN(2,2)*clhs147 + 0.5*clhs212;
            lhs(9,3)=0.5*DN(0,1)*clhs306 - 0.5*clhs178*clhs27 - 0.5*clhs179;
            lhs(9,4)=0.5*DN(2,0)*clhs47 + 0.5*DN(2,1)*clhs150 + 0.5*DN(2,2)*clhs151 + 0.5*clhs255;
            lhs(9,5)=clhs153*clhs304 + clhs155*clhs305 + clhs273 + clhs303*clhs65 + clhs327;
            lhs(9,6)=0.5*DN(2,0)*clhs72 + 0.5*DN(2,1)*clhs159 + 0.5*DN(2,2)*clhs161 + 0.5*clhs288;
            lhs(9,7)=0.5*DN(1,1)*clhs306 - 0.5*clhs27*clhs276 - 0.5*clhs277;
            lhs(9,8)=0.5*DN(2,0)*clhs80 + 0.5*DN(2,1)*clhs165 + 0.5*DN(2,2)*clhs166 + 0.5*clhs315;
            lhs(9,9)=clhs14*clhs328 + clhs168*clhs304 + clhs170*clhs305 + clhs303*clhs96 + clhs313;
            lhs(9,10)=0.5*DN(2,0)*clhs103 + 0.5*DN(2,1)*clhs174 + 0.5*DN(2,2)*clhs176 + 0.5*clhs330;
            lhs(9,11)=clhs304*clhs317;
            lhs(9,12)=0.5*DN(2,0)*clhs111 + 0.5*DN(2,1)*clhs180 + 0.5*DN(2,2)*clhs181 + 0.5*clhs331;
            lhs(9,13)=clhs126*clhs303 + clhs183*clhs304 + clhs185*clhs305 + clhs333 + clhs334;
            lhs(9,14)=0.5*DN(2,0)*clhs133 + 0.5*DN(2,1)*clhs189 + 0.5*DN(2,2)*clhs191 + 0.5*clhs335;
            lhs(9,15)=0.5*DN(3,1)*clhs306 - 0.5*clhs27*clhs337 - 0.5*clhs336;
            lhs(10,0)=0.5*DN(2,0)*clhs6 + 0.5*DN(2,1)*clhs140 + 0.5*DN(2,2)*clhs195 + 0.5*clhs106;
            lhs(10,1)=0.5*DN(2,0)*clhs35 + 0.5*DN(2,1)*clhs143 + 0.5*DN(2,2)*clhs196 + 0.5*clhs177;
            lhs(10,2)=clhs147*clhs304 + clhs197*clhs305 + clhs215 + clhs303*clhs42 + clhs326;
            lhs(10,3)=0.5*DN(0,2)*clhs306 - 0.5*clhs216*clhs27 - 0.5*clhs217;
            lhs(10,4)=0.5*DN(2,0)*clhs49 + 0.5*DN(2,1)*clhs151 + 0.5*DN(2,2)*clhs199 + 0.5*clhs256;
            lhs(10,5)=0.5*DN(2,0)*clhs68 + 0.5*DN(2,1)*clhs155 + 0.5*DN(2,2)*clhs202 + 0.5*clhs275;
            lhs(10,6)=clhs161*clhs304 + clhs204*clhs305 + clhs290 + clhs303*clhs74 + clhs327;
            lhs(10,7)=0.5*DN(1,2)*clhs306 - 0.5*clhs27*clhs291 - 0.5*clhs292;
            lhs(10,8)=0.5*DN(2,0)*clhs82 + 0.5*DN(2,1)*clhs166 + 0.5*DN(2,2)*clhs209 + 0.5*clhs316;
            lhs(10,9)=0.5*DN(2,0)*clhs99 + 0.5*DN(2,1)*clhs170 + 0.5*DN(2,2)*clhs211 + 0.5*clhs330;
            lhs(10,10)=clhs105*clhs303 + clhs14*clhs338 + clhs176*clhs304 + clhs213*clhs305 + clhs313;
            lhs(10,11)=clhs305*clhs317;
            lhs(10,12)=0.5*DN(2,0)*clhs113 + 0.5*DN(2,1)*clhs181 + 0.5*DN(2,2)*clhs218 + 0.5*clhs340;
            lhs(10,13)=0.5*DN(2,0)*clhs129 + 0.5*DN(2,1)*clhs185 + 0.5*DN(2,2)*clhs220 + 0.5*clhs341;
            lhs(10,14)=clhs135*clhs303 + clhs191*clhs304 + clhs222*clhs305 + clhs334 + clhs343;
            lhs(10,15)=0.5*DN(3,2)*clhs306 - 0.5*clhs27*clhs345 - 0.5*clhs344;
            lhs(11,0)=DN(2,0)*clhs299 + 0.5*clhs107;
            lhs(11,1)=DN(2,1)*clhs299 + 0.5*clhs178;
            lhs(11,2)=DN(2,2)*clhs299 + 0.5*clhs216;
            lhs(11,3)=clhs233;
            lhs(11,4)=DN(2,0)*clhs230 + 0.5*clhs257;
            lhs(11,5)=DN(2,1)*clhs230 + 0.5*clhs276;
            lhs(11,6)=DN(2,2)*clhs230 + 0.5*clhs291;
            lhs(11,7)=clhs301;
            lhs(11,8)=DN(2,0)*clhs346;
            lhs(11,9)=DN(2,1)*clhs346;
            lhs(11,10)=DN(2,2)*clhs346;
            lhs(11,11)=clhs228*(clhs311 + clhs328 + clhs338);
            lhs(11,12)=DN(2,0)*clhs234 + DN(3,0)*clhs88;
            lhs(11,13)=DN(2,1)*clhs234 + DN(3,1)*clhs88;
            lhs(11,14)=DN(2,2)*clhs234 + DN(3,2)*clhs88;
            lhs(11,15)=clhs347;
            lhs(12,0)=clhs0*clhs348 + clhs116 + clhs3*clhs349 + clhs350*clhs6 + clhs353;
            lhs(12,1)=0.5*DN(3,0)*clhs30 + 0.5*DN(3,1)*clhs32 + 0.5*DN(3,2)*clhs35 + 0.5*clhs182;
            lhs(12,2)=0.5*DN(3,0)*clhs38 + 0.5*DN(3,1)*clhs40 + 0.5*DN(3,2)*clhs42 + 0.5*clhs219;
            lhs(12,3)=0.5*DN(0,0)*clhs351 - 0.5*clhs137*clhs27 - 0.5*clhs138;
            lhs(12,4)=clhs261 + clhs348*clhs45 + clhs349*clhs47 + clhs350*clhs49 + clhs355;
            lhs(12,5)=0.5*DN(3,0)*clhs63 + 0.5*DN(3,1)*clhs65 + 0.5*DN(3,2)*clhs68 + 0.5*clhs278;
            lhs(12,6)=0.5*DN(3,0)*clhs70 + 0.5*DN(3,1)*clhs72 + 0.5*DN(3,2)*clhs74 + 0.5*clhs293;
            lhs(12,7)=0.5*DN(1,0)*clhs351 - 0.5*clhs265*clhs27 - 0.5*clhs266;
            lhs(12,8)=clhs320 + clhs348*clhs78 + clhs349*clhs80 + clhs350*clhs82 + clhs356;
            lhs(12,9)=0.5*DN(3,0)*clhs94 + 0.5*DN(3,1)*clhs96 + 0.5*DN(3,2)*clhs99 + 0.5*clhs331;
            lhs(12,10)=0.5*DN(3,0)*clhs101 + 0.5*DN(3,1)*clhs103 + 0.5*DN(3,2)*clhs105 + 0.5*clhs340;
            lhs(12,11)=0.5*DN(2,0)*clhs351 - 0.5*clhs27*clhs324 - 0.5*clhs325;
            lhs(12,12)=clhs109*clhs348 + clhs111*clhs349 + clhs113*clhs350 + clhs14*clhs357 + clhs359;
            lhs(12,13)=0.5*DN(3,0)*clhs124 + 0.5*DN(3,1)*clhs126 + 0.5*DN(3,2)*clhs129 + 0.5*clhs361;
            lhs(12,14)=0.5*DN(3,0)*clhs131 + 0.5*DN(3,1)*clhs133 + 0.5*DN(3,2)*clhs135 + 0.5*clhs362;
            lhs(12,15)=clhs348*clhs363;
            lhs(13,0)=0.5*DN(3,0)*clhs3 + 0.5*DN(3,1)*clhs139 + 0.5*DN(3,2)*clhs140 + 0.5*clhs130;
            lhs(13,1)=clhs141*clhs349 + clhs143*clhs350 + clhs187 + clhs32*clhs348 + clhs364;
            lhs(13,2)=0.5*DN(3,0)*clhs40 + 0.5*DN(3,1)*clhs145 + 0.5*DN(3,2)*clhs147 + 0.5*clhs221;
            lhs(13,3)=0.5*DN(0,1)*clhs351 - 0.5*clhs193*clhs27 - 0.5*clhs194;
            lhs(13,4)=0.5*DN(3,0)*clhs47 + 0.5*DN(3,1)*clhs150 + 0.5*DN(3,2)*clhs151 + 0.5*clhs263;
            lhs(13,5)=clhs153*clhs349 + clhs155*clhs350 + clhs280 + clhs348*clhs65 + clhs365;
            lhs(13,6)=0.5*DN(3,0)*clhs72 + 0.5*DN(3,1)*clhs159 + 0.5*DN(3,2)*clhs161 + 0.5*clhs294;
            lhs(13,7)=0.5*DN(1,1)*clhs351 - 0.5*clhs27*clhs283 - 0.5*clhs284;
            lhs(13,8)=0.5*DN(3,0)*clhs80 + 0.5*DN(3,1)*clhs165 + 0.5*DN(3,2)*clhs166 + 0.5*clhs322;
            lhs(13,9)=clhs168*clhs349 + clhs170*clhs350 + clhs333 + clhs348*clhs96 + clhs366;
            lhs(13,10)=0.5*DN(3,0)*clhs103 + 0.5*DN(3,1)*clhs174 + 0.5*DN(3,2)*clhs176 + 0.5*clhs341;
            lhs(13,11)=0.5*DN(2,1)*clhs351 - 0.5*clhs27*clhs336 - 0.5*clhs337;
            lhs(13,12)=0.5*DN(3,0)*clhs111 + 0.5*DN(3,1)*clhs180 + 0.5*DN(3,2)*clhs181 + 0.5*clhs361;
            lhs(13,13)=clhs126*clhs348 + clhs14*clhs367 + clhs183*clhs349 + clhs185*clhs350 + clhs359;
            lhs(13,14)=0.5*DN(3,0)*clhs133 + 0.5*DN(3,1)*clhs189 + 0.5*DN(3,2)*clhs191 + 0.5*clhs368;
            lhs(13,15)=clhs349*clhs363;
            lhs(14,0)=0.5*DN(3,0)*clhs6 + 0.5*DN(3,1)*clhs140 + 0.5*DN(3,2)*clhs195 + 0.5*clhs136;
            lhs(14,1)=0.5*DN(3,0)*clhs35 + 0.5*DN(3,1)*clhs143 + 0.5*DN(3,2)*clhs196 + 0.5*clhs192;
            lhs(14,2)=clhs147*clhs349 + clhs197*clhs350 + clhs224 + clhs348*clhs42 + clhs364;
            lhs(14,3)=0.5*DN(0,2)*clhs351 - 0.5*clhs225*clhs27 - 0.5*clhs226;
            lhs(14,4)=0.5*DN(3,0)*clhs49 + 0.5*DN(3,1)*clhs151 + 0.5*DN(3,2)*clhs199 + 0.5*clhs264;
            lhs(14,5)=0.5*DN(3,0)*clhs68 + 0.5*DN(3,1)*clhs155 + 0.5*DN(3,2)*clhs202 + 0.5*clhs282;
            lhs(14,6)=clhs161*clhs349 + clhs204*clhs350 + clhs296 + clhs348*clhs74 + clhs365;
            lhs(14,7)=0.5*DN(1,2)*clhs351 - 0.5*clhs27*clhs297 - 0.5*clhs298;
            lhs(14,8)=0.5*DN(3,0)*clhs82 + 0.5*DN(3,1)*clhs166 + 0.5*DN(3,2)*clhs209 + 0.5*clhs323;
            lhs(14,9)=0.5*DN(3,0)*clhs99 + 0.5*DN(3,1)*clhs170 + 0.5*DN(3,2)*clhs211 + 0.5*clhs335;
            lhs(14,10)=clhs105*clhs348 + clhs176*clhs349 + clhs213*clhs350 + clhs343 + clhs366;
            lhs(14,11)=0.5*DN(2,2)*clhs351 - 0.5*clhs27*clhs344 - 0.5*clhs345;
            lhs(14,12)=0.5*DN(3,0)*clhs113 + 0.5*DN(3,1)*clhs181 + 0.5*DN(3,2)*clhs218 + 0.5*clhs362;
            lhs(14,13)=0.5*DN(3,0)*clhs129 + 0.5*DN(3,1)*clhs185 + 0.5*DN(3,2)*clhs220 + 0.5*clhs368;
            lhs(14,14)=clhs135*clhs348 + clhs14*clhs369 + clhs191*clhs349 + clhs222*clhs350 + clhs359;
            lhs(14,15)=clhs350*clhs363;
            lhs(15,0)=DN(3,0)*clhs299 + 0.5*clhs137;
            lhs(15,1)=DN(3,1)*clhs299 + 0.5*clhs193;
            lhs(15,2)=DN(3,2)*clhs299 + 0.5*clhs225;
            lhs(15,3)=clhs235;
            lhs(15,4)=DN(3,0)*clhs230 + 0.5*clhs265;
            lhs(15,5)=DN(3,1)*clhs230 + 0.5*clhs283;
            lhs(15,6)=DN(3,2)*clhs230 + 0.5*clhs297;
            lhs(15,7)=clhs302;
            lhs(15,8)=DN(3,0)*clhs232 + 0.5*clhs324;
            lhs(15,9)=DN(3,1)*clhs232 + 0.5*clhs336;
            lhs(15,10)=DN(3,2)*clhs232 + 0.5*clhs344;
            lhs(15,11)=clhs347;
            lhs(15,12)=DN(3,0)*clhs370;
            lhs(15,13)=DN(3,1)*clhs370;
            lhs(15,14)=DN(3,2)*clhs370;
            lhs(15,15)=clhs228*(clhs357 + clhs367 + clhs369);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointRHSContribution(
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
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto v_CN = 0.5 * (vn + v);
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
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs1 =             p[0] + pn[0];
const double crhs2 =             p[1] + pn[1];
const double crhs3 =             p[2] + pn[2];
const double crhs4 =             N[0]*crhs1 + N[1]*crhs2 + N[2]*crhs3;
const double crhs5 =             0.5*DN(0,0);
const double crhs6 =             v(0,0) + vn(0,0);
const double crhs7 =             v(1,0) + vn(1,0);
const double crhs8 =             v(2,0) + vn(2,0);
const double crhs9 =             0.5*K_darcy;
const double crhs10 =             crhs9*(N[0]*crhs6 + N[1]*crhs7 + N[2]*crhs8);
const double crhs11 =             rho/dt;
const double crhs12 =             crhs11*(N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)));
const double crhs13 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double crhs14 =             DN(0,0)*crhs6;
const double crhs15 =             DN(1,0)*crhs7;
const double crhs16 =             DN(2,0)*crhs8;
const double crhs17 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double crhs18 =             0.5*rho;
const double crhs19 =             crhs18*(crhs13*(crhs14 + crhs15 + crhs16) + crhs17*(DN(0,1)*crhs6 + DN(1,1)*crhs7 + DN(2,1)*crhs8));
const double crhs20 =             v(0,1) + vn(0,1);
const double crhs21 =             DN(0,1)*crhs20;
const double crhs22 =             v(1,1) + vn(1,1);
const double crhs23 =             DN(1,1)*crhs22;
const double crhs24 =             v(2,1) + vn(2,1);
const double crhs25 =             DN(2,1)*crhs24;
const double crhs26 =             0.5*crhs14 + 0.5*crhs15 + 0.5*crhs16 + 0.5*crhs21 + 0.5*crhs23 + 0.5*crhs25 - volume_error_ratio;
const double crhs27 =             rho*stab_c2*sqrt(pow(crhs13, 2) + pow(crhs17, 2));
const double crhs28 =             crhs26*(crhs27*h/stab_c1 + mu);
const double crhs29 =             0.5*crhs2;
const double crhs30 =             0.5*crhs3;
const double crhs31 =             1.0/(K_darcy + crhs11*dyn_tau + crhs27/h + mu*stab_c1/pow(h, 2));
const double crhs32 =             crhs31*(DN(1,0)*crhs29 + DN(2,0)*crhs30 - crhs0 + crhs1*crhs5 + crhs10 + crhs12 + crhs19);
const double crhs33 =             K_darcy*N[0];
const double crhs34 =             rho*(DN(0,0)*crhs13 + DN(0,1)*crhs17);
const double crhs35 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs36 =             0.5*DN(0,1);
const double crhs37 =             crhs9*(N[0]*crhs20 + N[1]*crhs22 + N[2]*crhs24);
const double crhs38 =             crhs11*(N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)));
const double crhs39 =             crhs18*(crhs13*(DN(0,0)*crhs20 + DN(1,0)*crhs22 + DN(2,0)*crhs24) + crhs17*(crhs21 + crhs23 + crhs25));
const double crhs40 =             crhs31*(DN(1,1)*crhs29 + DN(2,1)*crhs30 + crhs1*crhs36 - crhs35 + crhs37 + crhs38 + crhs39);
const double crhs41 =             0.5*crhs4;
const double crhs42 =             K_darcy*N[1];
const double crhs43 =             rho*(DN(1,0)*crhs13 + DN(1,1)*crhs17);
const double crhs44 =             K_darcy*N[2];
const double crhs45 =             rho*(DN(2,0)*crhs13 + DN(2,1)*crhs17);
            rhs[0]=-DN(0,0)*crhs28 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs0 - N[0]*crhs10 - N[0]*crhs12 - N[0]*crhs19 + crhs32*crhs33 - crhs32*crhs34 + crhs4*crhs5;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs28 - DN(0,1)*stress[1] + N[0]*crhs35 - N[0]*crhs37 - N[0]*crhs38 - N[0]*crhs39 + crhs33*crhs40 - crhs34*crhs40 + crhs36*crhs4;
            rhs[2]=-DN(0,0)*crhs32 - DN(0,1)*crhs40 - N[0]*crhs26;
            rhs[3]=-DN(1,0)*crhs28 + DN(1,0)*crhs41 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs0 - N[1]*crhs10 - N[1]*crhs12 - N[1]*crhs19 + crhs32*crhs42 - crhs32*crhs43;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs28 + DN(1,1)*crhs41 - DN(1,1)*stress[1] + N[1]*crhs35 - N[1]*crhs37 - N[1]*crhs38 - N[1]*crhs39 + crhs40*crhs42 - crhs40*crhs43;
            rhs[5]=-DN(1,0)*crhs32 - DN(1,1)*crhs40 - N[1]*crhs26;
            rhs[6]=-DN(2,0)*crhs28 + DN(2,0)*crhs41 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs0 - N[2]*crhs10 - N[2]*crhs12 - N[2]*crhs19 + crhs32*crhs44 - crhs32*crhs45;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs28 + DN(2,1)*crhs41 - DN(2,1)*stress[1] + N[2]*crhs35 - N[2]*crhs37 - N[2]*crhs38 - N[2]*crhs39 + crhs40*crhs44 - crhs40*crhs45;
            rhs[8]=-DN(2,0)*crhs32 - DN(2,1)*crhs40 - N[2]*crhs26;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointRHSContribution(
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
    const auto &pn = rData.Pressure_OldStep1;

    const auto &stress = rData.ShearStress;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto v_CN = 0.5 * (vn + v);
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
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs1 =             p[0] + pn[0];
const double crhs2 =             p[1] + pn[1];
const double crhs3 =             p[2] + pn[2];
const double crhs4 =             p[3] + pn[3];
const double crhs5 =             N[0]*crhs1 + N[1]*crhs2 + N[2]*crhs3 + N[3]*crhs4;
const double crhs6 =             0.5*DN(0,0);
const double crhs7 =             v(0,0) + vn(0,0);
const double crhs8 =             v(1,0) + vn(1,0);
const double crhs9 =             v(2,0) + vn(2,0);
const double crhs10 =             v(3,0) + vn(3,0);
const double crhs11 =             0.5*K_darcy;
const double crhs12 =             crhs11*(N[0]*crhs7 + N[1]*crhs8 + N[2]*crhs9 + N[3]*crhs10);
const double crhs13 =             rho/dt;
const double crhs14 =             crhs13*(N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)) + N[3]*(v(3,0) - vn(3,0)));
const double crhs15 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double crhs16 =             DN(0,0)*crhs7;
const double crhs17 =             DN(1,0)*crhs8;
const double crhs18 =             DN(2,0)*crhs9;
const double crhs19 =             DN(3,0)*crhs10;
const double crhs20 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double crhs21 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double crhs22 =             0.5*rho;
const double crhs23 =             crhs22*(crhs15*(crhs16 + crhs17 + crhs18 + crhs19) + crhs20*(DN(0,1)*crhs7 + DN(1,1)*crhs8 + DN(2,1)*crhs9 + DN(3,1)*crhs10) + crhs21*(DN(0,2)*crhs7 + DN(1,2)*crhs8 + DN(2,2)*crhs9 + DN(3,2)*crhs10));
const double crhs24 =             v(0,1) + vn(0,1);
const double crhs25 =             DN(0,1)*crhs24;
const double crhs26 =             v(0,2) + vn(0,2);
const double crhs27 =             DN(0,2)*crhs26;
const double crhs28 =             v(1,1) + vn(1,1);
const double crhs29 =             DN(1,1)*crhs28;
const double crhs30 =             v(1,2) + vn(1,2);
const double crhs31 =             DN(1,2)*crhs30;
const double crhs32 =             v(2,1) + vn(2,1);
const double crhs33 =             DN(2,1)*crhs32;
const double crhs34 =             v(2,2) + vn(2,2);
const double crhs35 =             DN(2,2)*crhs34;
const double crhs36 =             v(3,1) + vn(3,1);
const double crhs37 =             DN(3,1)*crhs36;
const double crhs38 =             v(3,2) + vn(3,2);
const double crhs39 =             DN(3,2)*crhs38;
const double crhs40 =             0.5*crhs16 + 0.5*crhs17 + 0.5*crhs18 + 0.5*crhs19 + 0.5*crhs25 + 0.5*crhs27 + 0.5*crhs29 + 0.5*crhs31 + 0.5*crhs33 + 0.5*crhs35 + 0.5*crhs37 + 0.5*crhs39 - volume_error_ratio;
const double crhs41 =             rho*stab_c2*sqrt(pow(crhs15, 2) + pow(crhs20, 2) + pow(crhs21, 2));
const double crhs42 =             crhs40*(crhs41*h/stab_c1 + mu);
const double crhs43 =             0.5*crhs2;
const double crhs44 =             0.5*crhs3;
const double crhs45 =             0.5*crhs4;
const double crhs46 =             1.0/(K_darcy + crhs13*dyn_tau + crhs41/h + mu*stab_c1/pow(h, 2));
const double crhs47 =             crhs46*(DN(1,0)*crhs43 + DN(2,0)*crhs44 + DN(3,0)*crhs45 - crhs0 + crhs1*crhs6 + crhs12 + crhs14 + crhs23);
const double crhs48 =             K_darcy*N[0];
const double crhs49 =             rho*(DN(0,0)*crhs15 + DN(0,1)*crhs20 + DN(0,2)*crhs21);
const double crhs50 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs51 =             0.5*DN(0,1);
const double crhs52 =             crhs11*(N[0]*crhs24 + N[1]*crhs28 + N[2]*crhs32 + N[3]*crhs36);
const double crhs53 =             crhs13*(N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)) + N[3]*(v(3,1) - vn(3,1)));
const double crhs54 =             crhs22*(crhs15*(DN(0,0)*crhs24 + DN(1,0)*crhs28 + DN(2,0)*crhs32 + DN(3,0)*crhs36) + crhs20*(crhs25 + crhs29 + crhs33 + crhs37) + crhs21*(DN(0,2)*crhs24 + DN(1,2)*crhs28 + DN(2,2)*crhs32 + DN(3,2)*crhs36));
const double crhs55 =             crhs46*(DN(1,1)*crhs43 + DN(2,1)*crhs44 + DN(3,1)*crhs45 + crhs1*crhs51 - crhs50 + crhs52 + crhs53 + crhs54);
const double crhs56 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs57 =             0.5*DN(0,2);
const double crhs58 =             crhs11*(N[0]*crhs26 + N[1]*crhs30 + N[2]*crhs34 + N[3]*crhs38);
const double crhs59 =             crhs13*(N[0]*(v(0,2) - vn(0,2)) + N[1]*(v(1,2) - vn(1,2)) + N[2]*(v(2,2) - vn(2,2)) + N[3]*(v(3,2) - vn(3,2)));
const double crhs60 =             crhs22*(crhs15*(DN(0,0)*crhs26 + DN(1,0)*crhs30 + DN(2,0)*crhs34 + DN(3,0)*crhs38) + crhs20*(DN(0,1)*crhs26 + DN(1,1)*crhs30 + DN(2,1)*crhs34 + DN(3,1)*crhs38) + crhs21*(crhs27 + crhs31 + crhs35 + crhs39));
const double crhs61 =             crhs46*(DN(1,2)*crhs43 + DN(2,2)*crhs44 + DN(3,2)*crhs45 + crhs1*crhs57 - crhs56 + crhs58 + crhs59 + crhs60);
const double crhs62 =             0.5*crhs5;
const double crhs63 =             K_darcy*N[1];
const double crhs64 =             rho*(DN(1,0)*crhs15 + DN(1,1)*crhs20 + DN(1,2)*crhs21);
const double crhs65 =             K_darcy*N[2];
const double crhs66 =             rho*(DN(2,0)*crhs15 + DN(2,1)*crhs20 + DN(2,2)*crhs21);
const double crhs67 =             K_darcy*N[3];
const double crhs68 =             rho*(DN(3,0)*crhs15 + DN(3,1)*crhs20 + DN(3,2)*crhs21);
            rhs[0]=-DN(0,0)*crhs42 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs0 - N[0]*crhs12 - N[0]*crhs14 - N[0]*crhs23 + crhs47*crhs48 - crhs47*crhs49 + crhs5*crhs6;
            rhs[1]=-DN(0,0)*stress[3] - DN(0,1)*crhs42 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs50 - N[0]*crhs52 - N[0]*crhs53 - N[0]*crhs54 + crhs48*crhs55 - crhs49*crhs55 + crhs5*crhs51;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] - DN(0,2)*crhs42 - DN(0,2)*stress[2] + N[0]*crhs56 - N[0]*crhs58 - N[0]*crhs59 - N[0]*crhs60 + crhs48*crhs61 - crhs49*crhs61 + crhs5*crhs57;
            rhs[3]=-DN(0,0)*crhs47 - DN(0,1)*crhs55 - DN(0,2)*crhs61 - N[0]*crhs40;
            rhs[4]=-DN(1,0)*crhs42 + DN(1,0)*crhs62 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs0 - N[1]*crhs12 - N[1]*crhs14 - N[1]*crhs23 + crhs47*crhs63 - crhs47*crhs64;
            rhs[5]=-DN(1,0)*stress[3] - DN(1,1)*crhs42 + DN(1,1)*crhs62 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs50 - N[1]*crhs52 - N[1]*crhs53 - N[1]*crhs54 + crhs55*crhs63 - crhs55*crhs64;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] - DN(1,2)*crhs42 + DN(1,2)*crhs62 - DN(1,2)*stress[2] + N[1]*crhs56 - N[1]*crhs58 - N[1]*crhs59 - N[1]*crhs60 + crhs61*crhs63 - crhs61*crhs64;
            rhs[7]=-DN(1,0)*crhs47 - DN(1,1)*crhs55 - DN(1,2)*crhs61 - N[1]*crhs40;
            rhs[8]=-DN(2,0)*crhs42 + DN(2,0)*crhs62 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs0 - N[2]*crhs12 - N[2]*crhs14 - N[2]*crhs23 + crhs47*crhs65 - crhs47*crhs66;
            rhs[9]=-DN(2,0)*stress[3] - DN(2,1)*crhs42 + DN(2,1)*crhs62 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs50 - N[2]*crhs52 - N[2]*crhs53 - N[2]*crhs54 + crhs55*crhs65 - crhs55*crhs66;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] - DN(2,2)*crhs42 + DN(2,2)*crhs62 - DN(2,2)*stress[2] + N[2]*crhs56 - N[2]*crhs58 - N[2]*crhs59 - N[2]*crhs60 + crhs61*crhs65 - crhs61*crhs66;
            rhs[11]=-DN(2,0)*crhs47 - DN(2,1)*crhs55 - DN(2,2)*crhs61 - N[2]*crhs40;
            rhs[12]=-DN(3,0)*crhs42 + DN(3,0)*crhs62 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs0 - N[3]*crhs12 - N[3]*crhs14 - N[3]*crhs23 + crhs47*crhs67 - crhs47*crhs68;
            rhs[13]=-DN(3,0)*stress[3] - DN(3,1)*crhs42 + DN(3,1)*crhs62 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs50 - N[3]*crhs52 - N[3]*crhs53 - N[3]*crhs54 + crhs55*crhs67 - crhs55*crhs68;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] - DN(3,2)*crhs42 + DN(3,2)*crhs62 - DN(3,2)*stress[2] + N[3]*crhs56 - N[3]*crhs58 - N[3]*crhs59 - N[3]*crhs60 + crhs61*crhs67 - crhs61*crhs68;
            rhs[15]=-DN(3,0)*crhs47 - DN(3,1)*crhs55 - DN(3,2)*crhs61 - N[3]*crhs40;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointEnrichmentContributions(
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
    const auto &pn = rData.Pressure_OldStep1;
    const auto &p=rData.Pressure;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    const auto v_CN = 0.5 * (vn + v);
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
        volume_error_ratio = volume_error / dt;
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


    const double cH0 =             0.5*DN(0,0);
const double cH1 =             0.5*K_darcy;
const double cH2 =             1.0/dt;
const double cH3 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double cH4 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double cH5 =             0.5*DN(0,1);
const double cH6 =             1.0/(K_darcy + cH2*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2))/h + mu*stab_c1/pow(h, 2));
const double cH7 =             1.0*cH6;
const double cH8 =             cH7*(N[0]*cH1 + rho*(N[0]*cH2 + cH0*cH3 + cH4*cH5));
const double cH9 =             0.5*cH6;
const double cH10 =             0.5*DN(1,0);
const double cH11 =             0.5*DN(1,1);
const double cH12 =             cH7*(N[1]*cH1 + rho*(N[1]*cH2 + cH10*cH3 + cH11*cH4));
const double cH13 =             0.5*DN(2,0);
const double cH14 =             0.5*DN(2,1);
const double cH15 =             cH7*(N[2]*cH1 + rho*(N[2]*cH2 + cH13*cH3 + cH14*cH4));
            H(0,0)=DNenr(0,0)*cH8 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH8 + Nenr[0]*cH5;
            H(0,2)=cH9*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=DNenr(0,0)*cH12 + Nenr[0]*cH10;
            H(0,4)=DNenr(0,1)*cH12 + Nenr[0]*cH11;
            H(0,5)=cH9*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=DNenr(0,0)*cH15 + Nenr[0]*cH13;
            H(0,7)=DNenr(0,1)*cH15 + Nenr[0]*cH14;
            H(0,8)=cH9*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=DNenr(1,0)*cH8 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH8 + Nenr[1]*cH5;
            H(1,2)=cH9*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=DNenr(1,0)*cH12 + Nenr[1]*cH10;
            H(1,4)=DNenr(1,1)*cH12 + Nenr[1]*cH11;
            H(1,5)=cH9*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=DNenr(1,0)*cH15 + Nenr[1]*cH13;
            H(1,7)=DNenr(1,1)*cH15 + Nenr[1]*cH14;
            H(1,8)=cH9*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=DNenr(2,0)*cH8 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH8 + Nenr[2]*cH5;
            H(2,2)=cH9*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=DNenr(2,0)*cH12 + Nenr[2]*cH10;
            H(2,4)=DNenr(2,1)*cH12 + Nenr[2]*cH11;
            H(2,5)=cH9*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=DNenr(2,0)*cH15 + Nenr[2]*cH13;
            H(2,7)=DNenr(2,1)*cH15 + Nenr[2]*cH14;
            H(2,8)=cH9*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


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


    const double crhs_ee0 =             v(0,0) + vn(0,0);
const double crhs_ee1 =             DN(0,0)*crhs_ee0;
const double crhs_ee2 =             v(0,1) + vn(0,1);
const double crhs_ee3 =             DN(0,1)*crhs_ee2;
const double crhs_ee4 =             v(1,0) + vn(1,0);
const double crhs_ee5 =             DN(1,0)*crhs_ee4;
const double crhs_ee6 =             v(1,1) + vn(1,1);
const double crhs_ee7 =             DN(1,1)*crhs_ee6;
const double crhs_ee8 =             v(2,0) + vn(2,0);
const double crhs_ee9 =             DN(2,0)*crhs_ee8;
const double crhs_ee10 =             v(2,1) + vn(2,1);
const double crhs_ee11 =             DN(2,1)*crhs_ee10;
const double crhs_ee12 =             0.5*crhs_ee1 + 0.5*crhs_ee11 + 0.5*crhs_ee3 + 0.5*crhs_ee5 + 0.5*crhs_ee7 + 0.5*crhs_ee9 - volume_error_ratio;
const double crhs_ee13 =             0.5*p[0] + 0.5*pn[0];
const double crhs_ee14 =             0.5*p[1] + 0.5*pn[1];
const double crhs_ee15 =             0.5*p[2] + 0.5*pn[2];
const double crhs_ee16 =             0.5*K_darcy;
const double crhs_ee17 =             1.0/dt;
const double crhs_ee18 =             N[0]*crhs_ee17;
const double crhs_ee19 =             N[1]*crhs_ee17;
const double crhs_ee20 =             N[2]*crhs_ee17;
const double crhs_ee21 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0);
const double crhs_ee22 =             0.5*crhs_ee21;
const double crhs_ee23 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1);
const double crhs_ee24 =             0.5*crhs_ee23;
const double crhs_ee25 =             1.0/(K_darcy + crhs_ee17*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee21, 2) + pow(crhs_ee23, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee26 =             crhs_ee25*(DN(0,0)*crhs_ee13 + DN(1,0)*crhs_ee14 + DN(2,0)*crhs_ee15 + DNenr(0,0)*penr_cn[0] + DNenr(1,0)*penr_cn[1] + DNenr(2,0)*penr_cn[2] + crhs_ee16*(N[0]*crhs_ee0 + N[1]*crhs_ee4 + N[2]*crhs_ee8) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(crhs_ee18*(v(0,0) - vn(0,0)) + crhs_ee19*(v(1,0) - vn(1,0)) + crhs_ee20*(v(2,0) - vn(2,0)) + crhs_ee22*(crhs_ee1 + crhs_ee5 + crhs_ee9) + crhs_ee24*(DN(0,1)*crhs_ee0 + DN(1,1)*crhs_ee4 + DN(2,1)*crhs_ee8)));
const double crhs_ee27 =             crhs_ee25*(DN(0,1)*crhs_ee13 + DN(1,1)*crhs_ee14 + DN(2,1)*crhs_ee15 + DNenr(0,1)*penr_cn[0] + DNenr(1,1)*penr_cn[1] + DNenr(2,1)*penr_cn[2] + crhs_ee16*(N[0]*crhs_ee2 + N[1]*crhs_ee6 + N[2]*crhs_ee10) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(crhs_ee18*(v(0,1) - vn(0,1)) + crhs_ee19*(v(1,1) - vn(1,1)) + crhs_ee20*(v(2,1) - vn(2,1)) + crhs_ee22*(DN(0,0)*crhs_ee2 + DN(1,0)*crhs_ee6 + DN(2,0)*crhs_ee10) + crhs_ee24*(crhs_ee11 + crhs_ee3 + crhs_ee7)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee26 - DNenr(0,1)*crhs_ee27 - Nenr[0]*crhs_ee12;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee26 - DNenr(1,1)*crhs_ee27 - Nenr[1]*crhs_ee12;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee26 - DNenr(2,1)*crhs_ee27 - Nenr[2]*crhs_ee12;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesCN<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointEnrichmentContributions(
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
    const auto &pn = rData.Pressure_OldStep1;

    // TODO: Velocity CRANK NICOLSON at 0.5dt
    // FIXME: Mesh velocity should be evaluated at 0.5dt
    const auto v_CN = 0.5 * (vn + v);
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
        volume_error_ratio = volume_error / dt;
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


    const double cH0 =             0.5*DN(0,0);
const double cH1 =             0.5*K_darcy;
const double cH2 =             1.0/dt;
const double cH3 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double cH4 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double cH5 =             0.5*DN(0,1);
const double cH6 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double cH7 =             0.5*DN(0,2);
const double cH8 =             1.0/(K_darcy + cH2*dyn_tau*rho + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2) + pow(cH6, 2))/h + mu*stab_c1/pow(h, 2));
const double cH9 =             1.0*cH8;
const double cH10 =             cH9*(N[0]*cH1 + rho*(N[0]*cH2 + cH0*cH3 + cH4*cH5 + cH6*cH7));
const double cH11 =             0.5*cH8;
const double cH12 =             0.5*DN(1,0);
const double cH13 =             0.5*DN(1,1);
const double cH14 =             0.5*DN(1,2);
const double cH15 =             cH9*(N[1]*cH1 + rho*(N[1]*cH2 + cH12*cH3 + cH13*cH4 + cH14*cH6));
const double cH16 =             0.5*DN(2,0);
const double cH17 =             0.5*DN(2,1);
const double cH18 =             0.5*DN(2,2);
const double cH19 =             cH9*(N[2]*cH1 + rho*(N[2]*cH2 + cH16*cH3 + cH17*cH4 + cH18*cH6));
const double cH20 =             0.5*DN(3,0);
const double cH21 =             0.5*DN(3,1);
const double cH22 =             0.5*DN(3,2);
const double cH23 =             cH9*(N[3]*cH1 + rho*(N[3]*cH2 + cH20*cH3 + cH21*cH4 + cH22*cH6));
            H(0,0)=DNenr(0,0)*cH10 + Nenr[0]*cH0;
            H(0,1)=DNenr(0,1)*cH10 + Nenr[0]*cH5;
            H(0,2)=DNenr(0,2)*cH10 + Nenr[0]*cH7;
            H(0,3)=cH11*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=DNenr(0,0)*cH15 + Nenr[0]*cH12;
            H(0,5)=DNenr(0,1)*cH15 + Nenr[0]*cH13;
            H(0,6)=DNenr(0,2)*cH15 + Nenr[0]*cH14;
            H(0,7)=cH11*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=DNenr(0,0)*cH19 + Nenr[0]*cH16;
            H(0,9)=DNenr(0,1)*cH19 + Nenr[0]*cH17;
            H(0,10)=DNenr(0,2)*cH19 + Nenr[0]*cH18;
            H(0,11)=cH11*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=DNenr(0,0)*cH23 + Nenr[0]*cH20;
            H(0,13)=DNenr(0,1)*cH23 + Nenr[0]*cH21;
            H(0,14)=DNenr(0,2)*cH23 + Nenr[0]*cH22;
            H(0,15)=cH11*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=DNenr(1,0)*cH10 + Nenr[1]*cH0;
            H(1,1)=DNenr(1,1)*cH10 + Nenr[1]*cH5;
            H(1,2)=DNenr(1,2)*cH10 + Nenr[1]*cH7;
            H(1,3)=cH11*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=DNenr(1,0)*cH15 + Nenr[1]*cH12;
            H(1,5)=DNenr(1,1)*cH15 + Nenr[1]*cH13;
            H(1,6)=DNenr(1,2)*cH15 + Nenr[1]*cH14;
            H(1,7)=cH11*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=DNenr(1,0)*cH19 + Nenr[1]*cH16;
            H(1,9)=DNenr(1,1)*cH19 + Nenr[1]*cH17;
            H(1,10)=DNenr(1,2)*cH19 + Nenr[1]*cH18;
            H(1,11)=cH11*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=DNenr(1,0)*cH23 + Nenr[1]*cH20;
            H(1,13)=DNenr(1,1)*cH23 + Nenr[1]*cH21;
            H(1,14)=DNenr(1,2)*cH23 + Nenr[1]*cH22;
            H(1,15)=cH11*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=DNenr(2,0)*cH10 + Nenr[2]*cH0;
            H(2,1)=DNenr(2,1)*cH10 + Nenr[2]*cH5;
            H(2,2)=DNenr(2,2)*cH10 + Nenr[2]*cH7;
            H(2,3)=cH11*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=DNenr(2,0)*cH15 + Nenr[2]*cH12;
            H(2,5)=DNenr(2,1)*cH15 + Nenr[2]*cH13;
            H(2,6)=DNenr(2,2)*cH15 + Nenr[2]*cH14;
            H(2,7)=cH11*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=DNenr(2,0)*cH19 + Nenr[2]*cH16;
            H(2,9)=DNenr(2,1)*cH19 + Nenr[2]*cH17;
            H(2,10)=DNenr(2,2)*cH19 + Nenr[2]*cH18;
            H(2,11)=cH11*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=DNenr(2,0)*cH23 + Nenr[2]*cH20;
            H(2,13)=DNenr(2,1)*cH23 + Nenr[2]*cH21;
            H(2,14)=DNenr(2,2)*cH23 + Nenr[2]*cH22;
            H(2,15)=cH11*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=DNenr(3,0)*cH10 + Nenr[3]*cH0;
            H(3,1)=DNenr(3,1)*cH10 + Nenr[3]*cH5;
            H(3,2)=DNenr(3,2)*cH10 + Nenr[3]*cH7;
            H(3,3)=cH11*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=DNenr(3,0)*cH15 + Nenr[3]*cH12;
            H(3,5)=DNenr(3,1)*cH15 + Nenr[3]*cH13;
            H(3,6)=DNenr(3,2)*cH15 + Nenr[3]*cH14;
            H(3,7)=cH11*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=DNenr(3,0)*cH19 + Nenr[3]*cH16;
            H(3,9)=DNenr(3,1)*cH19 + Nenr[3]*cH17;
            H(3,10)=DNenr(3,2)*cH19 + Nenr[3]*cH18;
            H(3,11)=cH11*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=DNenr(3,0)*cH23 + Nenr[3]*cH20;
            H(3,13)=DNenr(3,1)*cH23 + Nenr[3]*cH21;
            H(3,14)=DNenr(3,2)*cH23 + Nenr[3]*cH22;
            H(3,15)=cH11*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


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


    const double crhs_ee0 =             v(0,0) + vn(0,0);
const double crhs_ee1 =             DN(0,0)*crhs_ee0;
const double crhs_ee2 =             v(0,1) + vn(0,1);
const double crhs_ee3 =             DN(0,1)*crhs_ee2;
const double crhs_ee4 =             v(0,2) + vn(0,2);
const double crhs_ee5 =             DN(0,2)*crhs_ee4;
const double crhs_ee6 =             v(1,0) + vn(1,0);
const double crhs_ee7 =             DN(1,0)*crhs_ee6;
const double crhs_ee8 =             v(1,1) + vn(1,1);
const double crhs_ee9 =             DN(1,1)*crhs_ee8;
const double crhs_ee10 =             v(1,2) + vn(1,2);
const double crhs_ee11 =             DN(1,2)*crhs_ee10;
const double crhs_ee12 =             v(2,0) + vn(2,0);
const double crhs_ee13 =             DN(2,0)*crhs_ee12;
const double crhs_ee14 =             v(2,1) + vn(2,1);
const double crhs_ee15 =             DN(2,1)*crhs_ee14;
const double crhs_ee16 =             v(2,2) + vn(2,2);
const double crhs_ee17 =             DN(2,2)*crhs_ee16;
const double crhs_ee18 =             v(3,0) + vn(3,0);
const double crhs_ee19 =             DN(3,0)*crhs_ee18;
const double crhs_ee20 =             v(3,1) + vn(3,1);
const double crhs_ee21 =             DN(3,1)*crhs_ee20;
const double crhs_ee22 =             v(3,2) + vn(3,2);
const double crhs_ee23 =             DN(3,2)*crhs_ee22;
const double crhs_ee24 =             0.5*crhs_ee1 + 0.5*crhs_ee11 + 0.5*crhs_ee13 + 0.5*crhs_ee15 + 0.5*crhs_ee17 + 0.5*crhs_ee19 + 0.5*crhs_ee21 + 0.5*crhs_ee23 + 0.5*crhs_ee3 + 0.5*crhs_ee5 + 0.5*crhs_ee7 + 0.5*crhs_ee9 - volume_error_ratio;
const double crhs_ee25 =             0.5*p[0] + 0.5*pn[0];
const double crhs_ee26 =             0.5*p[1] + 0.5*pn[1];
const double crhs_ee27 =             0.5*p[2] + 0.5*pn[2];
const double crhs_ee28 =             0.5*p[3] + 0.5*pn[3];
const double crhs_ee29 =             0.5*K_darcy;
const double crhs_ee30 =             1.0/dt;
const double crhs_ee31 =             N[0]*crhs_ee30;
const double crhs_ee32 =             N[1]*crhs_ee30;
const double crhs_ee33 =             N[2]*crhs_ee30;
const double crhs_ee34 =             N[3]*crhs_ee30;
const double crhs_ee35 =             N[0]*vconv_CN(0,0) + N[1]*vconv_CN(1,0) + N[2]*vconv_CN(2,0) + N[3]*vconv_CN(3,0);
const double crhs_ee36 =             0.5*crhs_ee35;
const double crhs_ee37 =             N[0]*vconv_CN(0,1) + N[1]*vconv_CN(1,1) + N[2]*vconv_CN(2,1) + N[3]*vconv_CN(3,1);
const double crhs_ee38 =             0.5*crhs_ee37;
const double crhs_ee39 =             N[0]*vconv_CN(0,2) + N[1]*vconv_CN(1,2) + N[2]*vconv_CN(2,2) + N[3]*vconv_CN(3,2);
const double crhs_ee40 =             0.5*crhs_ee39;
const double crhs_ee41 =             1.0/(K_darcy + crhs_ee30*dyn_tau*rho + rho*stab_c2*sqrt(pow(crhs_ee35, 2) + pow(crhs_ee37, 2) + pow(crhs_ee39, 2))/h + mu*stab_c1/pow(h, 2));
const double crhs_ee42 =             crhs_ee41*(DN(0,0)*crhs_ee25 + DN(1,0)*crhs_ee26 + DN(2,0)*crhs_ee27 + DN(3,0)*crhs_ee28 + DNenr(0,0)*penr_cn[0] + DNenr(1,0)*penr_cn[1] + DNenr(2,0)*penr_cn[2] + DNenr(3,0)*penr_cn[3] + crhs_ee29*(N[0]*crhs_ee0 + N[1]*crhs_ee6 + N[2]*crhs_ee12 + N[3]*crhs_ee18) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(crhs_ee31*(v(0,0) - vn(0,0)) + crhs_ee32*(v(1,0) - vn(1,0)) + crhs_ee33*(v(2,0) - vn(2,0)) + crhs_ee34*(v(3,0) - vn(3,0)) + crhs_ee36*(crhs_ee1 + crhs_ee13 + crhs_ee19 + crhs_ee7) + crhs_ee38*(DN(0,1)*crhs_ee0 + DN(1,1)*crhs_ee6 + DN(2,1)*crhs_ee12 + DN(3,1)*crhs_ee18) + crhs_ee40*(DN(0,2)*crhs_ee0 + DN(1,2)*crhs_ee6 + DN(2,2)*crhs_ee12 + DN(3,2)*crhs_ee18)));
const double crhs_ee43 =             crhs_ee41*(DN(0,1)*crhs_ee25 + DN(1,1)*crhs_ee26 + DN(2,1)*crhs_ee27 + DN(3,1)*crhs_ee28 + DNenr(0,1)*penr_cn[0] + DNenr(1,1)*penr_cn[1] + DNenr(2,1)*penr_cn[2] + DNenr(3,1)*penr_cn[3] + crhs_ee29*(N[0]*crhs_ee2 + N[1]*crhs_ee8 + N[2]*crhs_ee14 + N[3]*crhs_ee20) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(crhs_ee31*(v(0,1) - vn(0,1)) + crhs_ee32*(v(1,1) - vn(1,1)) + crhs_ee33*(v(2,1) - vn(2,1)) + crhs_ee34*(v(3,1) - vn(3,1)) + crhs_ee36*(DN(0,0)*crhs_ee2 + DN(1,0)*crhs_ee8 + DN(2,0)*crhs_ee14 + DN(3,0)*crhs_ee20) + crhs_ee38*(crhs_ee15 + crhs_ee21 + crhs_ee3 + crhs_ee9) + crhs_ee40*(DN(0,2)*crhs_ee2 + DN(1,2)*crhs_ee8 + DN(2,2)*crhs_ee14 + DN(3,2)*crhs_ee20)));
const double crhs_ee44 =             crhs_ee41*(DN(0,2)*crhs_ee25 + DN(1,2)*crhs_ee26 + DN(2,2)*crhs_ee27 + DN(3,2)*crhs_ee28 + DNenr(0,2)*penr_cn[0] + DNenr(1,2)*penr_cn[1] + DNenr(2,2)*penr_cn[2] + DNenr(3,2)*penr_cn[3] + crhs_ee29*(N[0]*crhs_ee4 + N[1]*crhs_ee10 + N[2]*crhs_ee16 + N[3]*crhs_ee22) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(crhs_ee31*(v(0,2) - vn(0,2)) + crhs_ee32*(v(1,2) - vn(1,2)) + crhs_ee33*(v(2,2) - vn(2,2)) + crhs_ee34*(v(3,2) - vn(3,2)) + crhs_ee36*(DN(0,0)*crhs_ee4 + DN(1,0)*crhs_ee10 + DN(2,0)*crhs_ee16 + DN(3,0)*crhs_ee22) + crhs_ee38*(DN(0,1)*crhs_ee4 + DN(1,1)*crhs_ee10 + DN(2,1)*crhs_ee16 + DN(3,1)*crhs_ee22) + crhs_ee40*(crhs_ee11 + crhs_ee17 + crhs_ee23 + crhs_ee5)));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee42 - DNenr(0,1)*crhs_ee43 - DNenr(0,2)*crhs_ee44 - Nenr[0]*crhs_ee24;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee42 - DNenr(1,1)*crhs_ee43 - DNenr(1,2)*crhs_ee44 - Nenr[1]*crhs_ee24;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee42 - DNenr(2,1)*crhs_ee43 - DNenr(2,2)*crhs_ee44 - Nenr[2]*crhs_ee24;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee42 - DNenr(3,1)*crhs_ee43 - DNenr(3,2)*crhs_ee44 - Nenr[3]*crhs_ee24;


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
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
   return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesCN< TwoFluidNavierStokesData<3, 4> >::pGetModifiedShapeFunctionsUtility(
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
    const auto v_CN = 0.5 * (vn + v);
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

template class TwoFluidNavierStokesCN<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokesCN<TwoFluidNavierStokesData<3, 4>>;

} // namespace Kratos
