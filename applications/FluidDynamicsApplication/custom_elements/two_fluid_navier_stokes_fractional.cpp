//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Autor:           Uxue Chasco
//

#include "two_fluid_navier_stokes_fractional.h"
#include "custom_elements/data_containers/two_fluid_fractional_navier_stokes/two_fluid_navier_stokes_fractional_data.h"
namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokesFractional<TElementData>::TwoFluidNavierStokesFractional(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokesFractional<TElementData>::TwoFluidNavierStokesFractional(
    IndexType NewId, const NodesArrayType &ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokesFractional<TElementData>::TwoFluidNavierStokesFractional(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokesFractional<TElementData>::TwoFluidNavierStokesFractional(
    IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokesFractional<TElementData>::~TwoFluidNavierStokesFractional() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer TwoFluidNavierStokesFractional<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesFractional>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer TwoFluidNavierStokesFractional<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<TwoFluidNavierStokesFractional>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::CalculateLocalSystem(
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

    if constexpr (TElementData::ElementManagesTimeIntegration){
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        if (data.IsCut()){
            GeometryType::Pointer p_geom = this->pGetGeometry();
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;

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
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
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
                    this->ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }


                GeometryType::ShapeFunctionsGradientsType int_shape_derivatives;
                Vector int_gauss_pts_weights;
                std::vector< array_1d<double,3> > int_normals_neg;

                // Without pressure gradient stabilization, volume ratio is checked during condensation
                CondenseEnrichmentWithContinuity(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);


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
        KRATOS_ERROR << "TwoFluidNavierStokesFractional is supposed to manage time integration." << std::endl;
    }
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokesFractional<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
const Parameters TwoFluidNavierStokesFractional<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISTANCE","VELOCITY","PRESSURE","MESH_VELOCITY","DENSITY","DYNAMIC_VISCOSITY","FRACTIONAL_VELOCITY"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["NewtonianTwoFluid2DLaw","NewtonianTwoFluid3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements Navier-Stokes biphasic fluid-air formulation with a levelset-based interface representation with Variational MultiScales (VMS) stabilization. Note that any viscous behavior can be used for the fluid phase through a constitutive law. The air phase is assumed to be Newtonian.
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TElementData>
std::string TwoFluidNavierStokesFractional<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokesFractional" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::PrintInfo(
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
void TwoFluidNavierStokesFractional<TElementData>::AddTimeIntegratedSystem(
    TElementData &rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::AddTimeIntegratedLHS(
    TElementData &rData,
    MatrixType &rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::AddTimeIntegratedRHS(
    TElementData &rData,
    VectorType &rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::UpdateIntegrationPointData(
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
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::UpdateIntegrationPointData(
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
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesFractionalData<2, 3> &rData,
    MatrixType &rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;
    const auto vn = rData.VelocityOldStep1;
    const auto vconv = rData.Velocity - rData.MeshVelocity;
    const auto vfrac = rData.FractionalVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    auto &lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = DN(0,0)*DN(0,0);
const double clhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 = rho*stab_c2*sqrt(clhs4*clhs4 + clhs5*clhs5);
const double clhs7 = clhs6*h*1.0/stab_c1 + mu;
const double clhs8 = DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs9 = N[0]*rho;
const double clhs10 = bdf0*rho;
const double clhs11 = N[0]*bdf0;
const double clhs12 = clhs11 + clhs8;
const double clhs13 = 1.0*1.0/(clhs6*1.0/h + dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h));
const double clhs14 = clhs13*(rho*rho);
const double clhs15 = clhs14*clhs8;
const double clhs16 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs17 = clhs14*clhs16;
const double clhs18 = N[0]*clhs17;
const double clhs19 = clhs10*(N[0]*N[0]) + clhs12*clhs15 + clhs12*clhs18 + clhs8*clhs9;
const double clhs20 = C(0,1)*DN(0,1) + clhs1;
const double clhs21 = C(1,2)*DN(0,1);
const double clhs22 = C(2,2)*DN(0,0) + clhs21;
const double clhs23 = DN(0,0)*clhs7;
const double clhs24 = DN(0,1)*clhs23;
const double clhs25 = clhs13*clhs16;
const double clhs26 = clhs13*rho;
const double clhs27 = clhs26*clhs8;
const double clhs28 = -N[0] + clhs25*clhs9 + clhs27;
const double clhs29 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs30 = C(0,2)*DN(1,0);
const double clhs31 = C(2,2)*DN(1,1) + clhs30;
const double clhs32 = DN(0,0)*DN(1,0);
const double clhs33 = N[1]*rho;
const double clhs34 = clhs11*clhs33;
const double clhs35 = clhs32*clhs7 + clhs34;
const double clhs36 = DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs37 = N[1]*bdf0;
const double clhs38 = clhs36 + clhs37;
const double clhs39 = clhs15*clhs38 + clhs18*clhs38 + clhs36*clhs9;
const double clhs40 = C(0,1)*DN(1,1) + clhs30;
const double clhs41 = C(1,2)*DN(1,1);
const double clhs42 = C(2,2)*DN(1,0) + clhs41;
const double clhs43 = DN(1,1)*clhs23;
const double clhs44 = DN(0,0)*N[1];
const double clhs45 = DN(1,0)*N[0];
const double clhs46 = clhs16*clhs26;
const double clhs47 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs48 = C(0,2)*DN(2,0);
const double clhs49 = C(2,2)*DN(2,1) + clhs48;
const double clhs50 = DN(0,0)*DN(2,0);
const double clhs51 = N[2]*rho;
const double clhs52 = clhs11*clhs51;
const double clhs53 = clhs50*clhs7 + clhs52;
const double clhs54 = DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs55 = N[2]*bdf0 + clhs54;
const double clhs56 = clhs15*clhs55 + clhs18*clhs55 + clhs54*clhs9;
const double clhs57 = C(0,1)*DN(2,1) + clhs48;
const double clhs58 = C(1,2)*DN(2,1);
const double clhs59 = C(2,2)*DN(2,0) + clhs58;
const double clhs60 = DN(2,1)*clhs23;
const double clhs61 = DN(0,0)*N[2];
const double clhs62 = DN(2,0)*N[0];
const double clhs63 = C(0,1)*DN(0,0) + clhs21;
const double clhs64 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs65 = DN(0,1)*DN(0,1);
const double clhs66 = C(0,1)*DN(1,0) + clhs41;
const double clhs67 = DN(0,1)*clhs7;
const double clhs68 = DN(1,0)*clhs67;
const double clhs69 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs70 = DN(0,1)*DN(1,1);
const double clhs71 = clhs34 + clhs7*clhs70;
const double clhs72 = DN(0,1)*N[1];
const double clhs73 = DN(1,1)*N[0];
const double clhs74 = C(0,1)*DN(2,0) + clhs58;
const double clhs75 = DN(2,0)*clhs67;
const double clhs76 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs77 = DN(0,1)*DN(2,1);
const double clhs78 = clhs52 + clhs7*clhs77;
const double clhs79 = DN(0,1)*N[2];
const double clhs80 = DN(2,1)*N[0];
const double clhs81 = clhs12*clhs26;
const double clhs82 = N[0] + clhs81;
const double clhs83 = clhs26*clhs38;
const double clhs84 = clhs13*(clhs32 + clhs70);
const double clhs85 = clhs26*clhs55;
const double clhs86 = clhs13*(clhs50 + clhs77);
const double clhs87 = clhs14*clhs36;
const double clhs88 = N[1]*clhs17;
const double clhs89 = clhs12*clhs87 + clhs12*clhs88 + clhs33*clhs8;
const double clhs90 = clhs26*clhs36;
const double clhs91 = DN(1,0)*DN(1,0);
const double clhs92 = clhs10*(N[1]*N[1]) + clhs33*clhs36 + clhs38*clhs87 + clhs38*clhs88;
const double clhs93 = DN(1,0)*clhs7;
const double clhs94 = DN(1,1)*clhs93;
const double clhs95 = -N[1] + clhs25*clhs33 + clhs90;
const double clhs96 = DN(1,0)*DN(2,0);
const double clhs97 = clhs37*clhs51;
const double clhs98 = clhs7*clhs96 + clhs97;
const double clhs99 = clhs33*clhs54 + clhs55*clhs87 + clhs55*clhs88;
const double clhs100 = DN(2,1)*clhs93;
const double clhs101 = DN(1,0)*N[2];
const double clhs102 = DN(2,0)*N[1];
const double clhs103 = DN(1,1)*DN(1,1);
const double clhs104 = DN(2,0)*clhs7;
const double clhs105 = DN(1,1)*clhs104;
const double clhs106 = DN(1,1)*DN(2,1);
const double clhs107 = clhs106*clhs7 + clhs97;
const double clhs108 = DN(1,1)*N[2];
const double clhs109 = DN(2,1)*N[1];
const double clhs110 = N[1] + clhs83;
const double clhs111 = clhs13*(clhs106 + clhs96);
const double clhs112 = clhs14*clhs54;
const double clhs113 = N[2]*clhs17;
const double clhs114 = clhs112*clhs12 + clhs113*clhs12 + clhs51*clhs8;
const double clhs115 = clhs26*clhs54;
const double clhs116 = clhs112*clhs38 + clhs113*clhs38 + clhs36*clhs51;
const double clhs117 = DN(2,0)*DN(2,0);
const double clhs118 = clhs10*(N[2]*N[2]) + clhs112*clhs55 + clhs113*clhs55 + clhs51*clhs54;
const double clhs119 = DN(2,1)*clhs104;
const double clhs120 = -N[2] + clhs115 + clhs25*clhs51;
const double clhs121 = DN(2,1)*DN(2,1);
const double clhs122 = N[2] + clhs85;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
lhs(0,2)=DN(0,0)*clhs28;
lhs(0,3)=DN(0,0)*clhs29 + DN(0,1)*clhs31 + clhs35 + clhs39;
lhs(0,4)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + clhs43;
lhs(0,5)=DN(1,0)*clhs27 - clhs44 + clhs45*clhs46;
lhs(0,6)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + clhs53 + clhs56;
lhs(0,7)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + clhs60;
lhs(0,8)=DN(2,0)*clhs27 + clhs46*clhs62 - clhs61;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs63 + clhs24;
lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs64 + clhs19 + clhs65*clhs7;
lhs(1,2)=DN(0,1)*clhs28;
lhs(1,3)=DN(0,0)*clhs31 + DN(0,1)*clhs66 + clhs68;
lhs(1,4)=DN(0,0)*clhs42 + DN(0,1)*clhs69 + clhs39 + clhs71;
lhs(1,5)=DN(1,1)*clhs27 + clhs46*clhs73 - clhs72;
lhs(1,6)=DN(0,0)*clhs49 + DN(0,1)*clhs74 + clhs75;
lhs(1,7)=DN(0,0)*clhs59 + DN(0,1)*clhs76 + clhs56 + clhs78;
lhs(1,8)=DN(2,1)*clhs27 + clhs46*clhs80 - clhs79;
lhs(2,0)=DN(0,0)*clhs82;
lhs(2,1)=DN(0,1)*clhs82;
lhs(2,2)=clhs13*(clhs3 + clhs65);
lhs(2,3)=DN(0,0)*clhs83 + clhs45;
lhs(2,4)=DN(0,1)*clhs83 + clhs73;
lhs(2,5)=clhs84;
lhs(2,6)=DN(0,0)*clhs85 + clhs62;
lhs(2,7)=DN(0,1)*clhs85 + clhs80;
lhs(2,8)=clhs86;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs35 + clhs89;
lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs68;
lhs(3,2)=DN(0,0)*clhs90 + clhs44*clhs46 - clhs45;
lhs(3,3)=DN(1,0)*clhs29 + DN(1,1)*clhs31 + clhs7*clhs91 + clhs92;
lhs(3,4)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + clhs94;
lhs(3,5)=DN(1,0)*clhs95;
lhs(3,6)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + clhs98 + clhs99;
lhs(3,7)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + clhs100;
lhs(3,8)=DN(2,0)*clhs90 - clhs101 + clhs102*clhs46;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs63 + clhs43;
lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs64 + clhs71 + clhs89;
lhs(4,2)=DN(0,1)*clhs90 + clhs46*clhs72 - clhs73;
lhs(4,3)=DN(1,0)*clhs31 + DN(1,1)*clhs66 + clhs94;
lhs(4,4)=DN(1,0)*clhs42 + DN(1,1)*clhs69 + clhs103*clhs7 + clhs92;
lhs(4,5)=DN(1,1)*clhs95;
lhs(4,6)=DN(1,0)*clhs49 + DN(1,1)*clhs74 + clhs105;
lhs(4,7)=DN(1,0)*clhs59 + DN(1,1)*clhs76 + clhs107 + clhs99;
lhs(4,8)=DN(2,1)*clhs90 - clhs108 + clhs109*clhs46;
lhs(5,0)=DN(1,0)*clhs81 + clhs44;
lhs(5,1)=DN(1,1)*clhs81 + clhs72;
lhs(5,2)=clhs84;
lhs(5,3)=DN(1,0)*clhs110;
lhs(5,4)=DN(1,1)*clhs110;
lhs(5,5)=clhs13*(clhs103 + clhs91);
lhs(5,6)=DN(1,0)*clhs85 + clhs102;
lhs(5,7)=DN(1,1)*clhs85 + clhs109;
lhs(5,8)=clhs111;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs114 + clhs53;
lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs75;
lhs(6,2)=DN(0,0)*clhs115 + clhs46*clhs61 - clhs62;
lhs(6,3)=DN(2,0)*clhs29 + DN(2,1)*clhs31 + clhs116 + clhs98;
lhs(6,4)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + clhs105;
lhs(6,5)=DN(1,0)*clhs115 + clhs101*clhs46 - clhs102;
lhs(6,6)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + clhs117*clhs7 + clhs118;
lhs(6,7)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + clhs119;
lhs(6,8)=DN(2,0)*clhs120;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs63 + clhs60;
lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs64 + clhs114 + clhs78;
lhs(7,2)=DN(0,1)*clhs115 + clhs46*clhs79 - clhs80;
lhs(7,3)=DN(2,0)*clhs31 + DN(2,1)*clhs66 + clhs100;
lhs(7,4)=DN(2,0)*clhs42 + DN(2,1)*clhs69 + clhs107 + clhs116;
lhs(7,5)=DN(1,1)*clhs115 + clhs108*clhs46 - clhs109;
lhs(7,6)=DN(2,0)*clhs49 + DN(2,1)*clhs74 + clhs119;
lhs(7,7)=DN(2,0)*clhs59 + DN(2,1)*clhs76 + clhs118 + clhs121*clhs7;
lhs(7,8)=DN(2,1)*clhs120;
lhs(8,0)=DN(2,0)*clhs81 + clhs61;
lhs(8,1)=DN(2,1)*clhs81 + clhs79;
lhs(8,2)=clhs86;
lhs(8,3)=DN(2,0)*clhs83 + clhs101;
lhs(8,4)=DN(2,1)*clhs83 + clhs108;
lhs(8,5)=clhs111;
lhs(8,6)=DN(2,0)*clhs122;
lhs(8,7)=DN(2,1)*clhs122;
lhs(8,8)=clhs13*(clhs117 + clhs121);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesFractionalData<3, 4> &rData,
    MatrixType &rLHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;
    const auto vn=rData.VelocityOldStep1;
    const auto vconv = rData.Velocity - rData.MeshVelocity;
    const auto vfrac = rData.FractionalVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    auto &lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = DN(0,0)*DN(0,0);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 = rho*stab_c2*sqrt(clhs6*clhs6 + clhs7*clhs7 + clhs8*clhs8);
const double clhs10 = clhs9*h*1.0/stab_c1 + mu;
const double clhs11 = DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs12 = N[0]*rho;
const double clhs13 = bdf0*rho;
const double clhs14 = N[0]*bdf0;
const double clhs15 = clhs11 + clhs14;
const double clhs16 = 1.0*1.0/(clhs9*1.0/h + dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h));
const double clhs17 = clhs16*(rho*rho);
const double clhs18 = clhs11*clhs17;
const double clhs19 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs20 = clhs17*clhs19;
const double clhs21 = N[0]*clhs20;
const double clhs22 = clhs11*clhs12 + clhs13*(N[0]*N[0]) + clhs15*clhs18 + clhs15*clhs21;
const double clhs23 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs24 = C(1,3)*DN(0,1);
const double clhs25 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
const double clhs26 = C(3,5)*DN(0,0);
const double clhs27 = C(4,5)*DN(0,2);
const double clhs28 = C(1,5)*DN(0,1) + clhs26 + clhs27;
const double clhs29 = DN(0,0)*clhs10;
const double clhs30 = DN(0,1)*clhs29;
const double clhs31 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs32 = C(3,4)*DN(0,1);
const double clhs33 = C(2,3)*DN(0,2) + clhs26 + clhs32;
const double clhs34 = C(2,5)*DN(0,2);
const double clhs35 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
const double clhs36 = DN(0,2)*clhs29;
const double clhs37 = clhs16*clhs19;
const double clhs38 = clhs16*rho;
const double clhs39 = clhs11*clhs38;
const double clhs40 = -N[0] + clhs12*clhs37 + clhs39;
const double clhs41 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs42 = C(0,3)*DN(1,0);
const double clhs43 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs42;
const double clhs44 = C(0,5)*DN(1,0);
const double clhs45 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs44;
const double clhs46 = DN(0,0)*DN(1,0);
const double clhs47 = N[1]*rho;
const double clhs48 = clhs14*clhs47;
const double clhs49 = clhs10*clhs46 + clhs48;
const double clhs50 = DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs51 = N[1]*bdf0;
const double clhs52 = clhs50 + clhs51;
const double clhs53 = clhs12*clhs50 + clhs18*clhs52 + clhs21*clhs52;
const double clhs54 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs42;
const double clhs55 = C(1,3)*DN(1,1);
const double clhs56 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs55;
const double clhs57 = C(3,5)*DN(1,0);
const double clhs58 = C(4,5)*DN(1,2);
const double clhs59 = C(1,5)*DN(1,1) + clhs57 + clhs58;
const double clhs60 = DN(1,1)*clhs29;
const double clhs61 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs44;
const double clhs62 = C(3,4)*DN(1,1);
const double clhs63 = C(2,3)*DN(1,2) + clhs57 + clhs62;
const double clhs64 = C(2,5)*DN(1,2);
const double clhs65 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs64;
const double clhs66 = DN(1,2)*clhs29;
const double clhs67 = DN(0,0)*N[1];
const double clhs68 = DN(1,0)*N[0];
const double clhs69 = clhs19*clhs38;
const double clhs70 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs71 = C(0,3)*DN(2,0);
const double clhs72 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs71;
const double clhs73 = C(0,5)*DN(2,0);
const double clhs74 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs73;
const double clhs75 = DN(0,0)*DN(2,0);
const double clhs76 = N[2]*rho;
const double clhs77 = clhs14*clhs76;
const double clhs78 = clhs10*clhs75 + clhs77;
const double clhs79 = DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs80 = N[2]*bdf0;
const double clhs81 = clhs79 + clhs80;
const double clhs82 = clhs12*clhs79 + clhs18*clhs81 + clhs21*clhs81;
const double clhs83 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs71;
const double clhs84 = C(1,3)*DN(2,1);
const double clhs85 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs84;
const double clhs86 = C(3,5)*DN(2,0);
const double clhs87 = C(4,5)*DN(2,2);
const double clhs88 = C(1,5)*DN(2,1) + clhs86 + clhs87;
const double clhs89 = DN(2,1)*clhs29;
const double clhs90 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs73;
const double clhs91 = C(3,4)*DN(2,1);
const double clhs92 = C(2,3)*DN(2,2) + clhs86 + clhs91;
const double clhs93 = C(2,5)*DN(2,2);
const double clhs94 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs93;
const double clhs95 = DN(2,2)*clhs29;
const double clhs96 = DN(0,0)*N[2];
const double clhs97 = DN(2,0)*N[0];
const double clhs98 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs99 = C(0,3)*DN(3,0);
const double clhs100 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs99;
const double clhs101 = C(0,5)*DN(3,0);
const double clhs102 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs101;
const double clhs103 = DN(0,0)*DN(3,0);
const double clhs104 = N[3]*rho;
const double clhs105 = clhs104*clhs14;
const double clhs106 = clhs10*clhs103 + clhs105;
const double clhs107 = DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs108 = N[3]*bdf0 + clhs107;
const double clhs109 = clhs107*clhs12 + clhs108*clhs18 + clhs108*clhs21;
const double clhs110 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs99;
const double clhs111 = C(1,3)*DN(3,1);
const double clhs112 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs111;
const double clhs113 = C(3,5)*DN(3,0);
const double clhs114 = C(4,5)*DN(3,2);
const double clhs115 = C(1,5)*DN(3,1) + clhs113 + clhs114;
const double clhs116 = DN(3,1)*clhs29;
const double clhs117 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs101;
const double clhs118 = C(3,4)*DN(3,1);
const double clhs119 = C(2,3)*DN(3,2) + clhs113 + clhs118;
const double clhs120 = C(2,5)*DN(3,2);
const double clhs121 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs120;
const double clhs122 = DN(3,2)*clhs29;
const double clhs123 = DN(0,0)*N[3];
const double clhs124 = DN(3,0)*N[0];
const double clhs125 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs126 = C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs127 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs128 = C(1,4)*DN(0,1);
const double clhs129 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs128;
const double clhs130 = DN(0,1)*DN(0,1);
const double clhs131 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs128;
const double clhs132 = C(2,4)*DN(0,2);
const double clhs133 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs132;
const double clhs134 = DN(0,1)*clhs10;
const double clhs135 = DN(0,2)*clhs134;
const double clhs136 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs55;
const double clhs137 = C(0,4)*DN(1,0) + clhs58 + clhs62;
const double clhs138 = DN(1,0)*clhs134;
const double clhs139 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs140 = C(1,4)*DN(1,1);
const double clhs141 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs140;
const double clhs142 = DN(0,1)*DN(1,1);
const double clhs143 = clhs10*clhs142;
const double clhs144 = clhs48 + clhs53;
const double clhs145 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs140;
const double clhs146 = C(2,4)*DN(1,2);
const double clhs147 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs146;
const double clhs148 = DN(1,2)*clhs134;
const double clhs149 = DN(0,1)*N[1];
const double clhs150 = DN(1,1)*N[0];
const double clhs151 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs84;
const double clhs152 = C(0,4)*DN(2,0) + clhs87 + clhs91;
const double clhs153 = DN(2,0)*clhs134;
const double clhs154 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs155 = C(1,4)*DN(2,1);
const double clhs156 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs155;
const double clhs157 = DN(0,1)*DN(2,1);
const double clhs158 = clhs10*clhs157;
const double clhs159 = clhs77 + clhs82;
const double clhs160 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs155;
const double clhs161 = C(2,4)*DN(2,2);
const double clhs162 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs161;
const double clhs163 = DN(2,2)*clhs134;
const double clhs164 = DN(0,1)*N[2];
const double clhs165 = DN(2,1)*N[0];
const double clhs166 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs111;
const double clhs167 = C(0,4)*DN(3,0) + clhs114 + clhs118;
const double clhs168 = DN(3,0)*clhs134;
const double clhs169 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs170 = C(1,4)*DN(3,1);
const double clhs171 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs170;
const double clhs172 = DN(0,1)*DN(3,1);
const double clhs173 = clhs10*clhs172;
const double clhs174 = clhs105 + clhs109;
const double clhs175 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs170;
const double clhs176 = C(2,4)*DN(3,2);
const double clhs177 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs176;
const double clhs178 = DN(3,2)*clhs134;
const double clhs179 = DN(0,1)*N[3];
const double clhs180 = DN(3,1)*N[0];
const double clhs181 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs182 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs132;
const double clhs183 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs184 = DN(0,2)*DN(0,2);
const double clhs185 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs64;
const double clhs186 = DN(0,2)*clhs10;
const double clhs187 = DN(1,0)*clhs186;
const double clhs188 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs146;
const double clhs189 = DN(1,1)*clhs186;
const double clhs190 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs191 = DN(0,2)*DN(1,2);
const double clhs192 = clhs10*clhs191;
const double clhs193 = DN(0,2)*N[1];
const double clhs194 = DN(1,2)*N[0];
const double clhs195 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs93;
const double clhs196 = DN(2,0)*clhs186;
const double clhs197 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs161;
const double clhs198 = DN(2,1)*clhs186;
const double clhs199 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs200 = DN(0,2)*DN(2,2);
const double clhs201 = clhs10*clhs200;
const double clhs202 = DN(0,2)*N[2];
const double clhs203 = DN(2,2)*N[0];
const double clhs204 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs120;
const double clhs205 = DN(3,0)*clhs186;
const double clhs206 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs176;
const double clhs207 = DN(3,1)*clhs186;
const double clhs208 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs209 = DN(0,2)*DN(3,2);
const double clhs210 = clhs10*clhs209;
const double clhs211 = DN(0,2)*N[3];
const double clhs212 = DN(3,2)*N[0];
const double clhs213 = clhs15*clhs38;
const double clhs214 = N[0] + clhs213;
const double clhs215 = clhs38*clhs52;
const double clhs216 = clhs16*(clhs142 + clhs191 + clhs46);
const double clhs217 = clhs38*clhs81;
const double clhs218 = clhs16*(clhs157 + clhs200 + clhs75);
const double clhs219 = clhs108*clhs38;
const double clhs220 = clhs16*(clhs103 + clhs172 + clhs209);
const double clhs221 = clhs17*clhs50;
const double clhs222 = N[1]*clhs20;
const double clhs223 = clhs11*clhs47 + clhs15*clhs221 + clhs15*clhs222;
const double clhs224 = clhs38*clhs50;
const double clhs225 = DN(1,0)*DN(1,0);
const double clhs226 = clhs13*(N[1]*N[1]) + clhs221*clhs52 + clhs222*clhs52 + clhs47*clhs50;
const double clhs227 = DN(1,0)*clhs10;
const double clhs228 = DN(1,1)*clhs227;
const double clhs229 = DN(1,2)*clhs227;
const double clhs230 = -N[1] + clhs224 + clhs37*clhs47;
const double clhs231 = DN(1,0)*DN(2,0);
const double clhs232 = clhs51*clhs76;
const double clhs233 = clhs10*clhs231 + clhs232;
const double clhs234 = clhs221*clhs81 + clhs222*clhs81 + clhs47*clhs79;
const double clhs235 = DN(2,1)*clhs227;
const double clhs236 = DN(2,2)*clhs227;
const double clhs237 = DN(1,0)*N[2];
const double clhs238 = DN(2,0)*N[1];
const double clhs239 = DN(1,0)*DN(3,0);
const double clhs240 = clhs104*clhs51;
const double clhs241 = clhs10*clhs239 + clhs240;
const double clhs242 = clhs107*clhs47 + clhs108*clhs221 + clhs108*clhs222;
const double clhs243 = DN(3,1)*clhs227;
const double clhs244 = DN(3,2)*clhs227;
const double clhs245 = DN(1,0)*N[3];
const double clhs246 = DN(3,0)*N[1];
const double clhs247 = clhs223 + clhs48;
const double clhs248 = DN(1,1)*DN(1,1);
const double clhs249 = DN(1,1)*clhs10;
const double clhs250 = DN(1,2)*clhs249;
const double clhs251 = DN(2,0)*clhs249;
const double clhs252 = DN(1,1)*DN(2,1);
const double clhs253 = clhs10*clhs252;
const double clhs254 = clhs232 + clhs234;
const double clhs255 = DN(2,2)*clhs249;
const double clhs256 = DN(1,1)*N[2];
const double clhs257 = DN(2,1)*N[1];
const double clhs258 = DN(3,0)*clhs249;
const double clhs259 = DN(1,1)*DN(3,1);
const double clhs260 = clhs10*clhs259;
const double clhs261 = clhs240 + clhs242;
const double clhs262 = DN(3,2)*clhs249;
const double clhs263 = DN(1,1)*N[3];
const double clhs264 = DN(3,1)*N[1];
const double clhs265 = DN(1,2)*DN(1,2);
const double clhs266 = DN(1,2)*clhs10;
const double clhs267 = DN(2,0)*clhs266;
const double clhs268 = DN(2,1)*clhs266;
const double clhs269 = DN(1,2)*DN(2,2);
const double clhs270 = clhs10*clhs269;
const double clhs271 = DN(1,2)*N[2];
const double clhs272 = DN(2,2)*N[1];
const double clhs273 = DN(3,0)*clhs266;
const double clhs274 = DN(3,1)*clhs266;
const double clhs275 = DN(1,2)*DN(3,2);
const double clhs276 = clhs10*clhs275;
const double clhs277 = DN(1,2)*N[3];
const double clhs278 = DN(3,2)*N[1];
const double clhs279 = N[1] + clhs215;
const double clhs280 = clhs16*(clhs231 + clhs252 + clhs269);
const double clhs281 = clhs16*(clhs239 + clhs259 + clhs275);
const double clhs282 = clhs17*clhs79;
const double clhs283 = N[2]*clhs20;
const double clhs284 = clhs11*clhs76 + clhs15*clhs282 + clhs15*clhs283;
const double clhs285 = clhs38*clhs79;
const double clhs286 = clhs282*clhs52 + clhs283*clhs52 + clhs50*clhs76;
const double clhs287 = DN(2,0)*DN(2,0);
const double clhs288 = clhs13*(N[2]*N[2]) + clhs282*clhs81 + clhs283*clhs81 + clhs76*clhs79;
const double clhs289 = DN(2,0)*clhs10;
const double clhs290 = DN(2,1)*clhs289;
const double clhs291 = DN(2,2)*clhs289;
const double clhs292 = -N[2] + clhs285 + clhs37*clhs76;
const double clhs293 = DN(2,0)*DN(3,0);
const double clhs294 = clhs104*clhs80;
const double clhs295 = clhs10*clhs293 + clhs294;
const double clhs296 = clhs107*clhs76 + clhs108*clhs282 + clhs108*clhs283;
const double clhs297 = DN(3,1)*clhs289;
const double clhs298 = DN(3,2)*clhs289;
const double clhs299 = DN(2,0)*N[3];
const double clhs300 = DN(3,0)*N[2];
const double clhs301 = clhs284 + clhs77;
const double clhs302 = clhs232 + clhs286;
const double clhs303 = DN(2,1)*DN(2,1);
const double clhs304 = DN(2,1)*clhs10;
const double clhs305 = DN(2,2)*clhs304;
const double clhs306 = DN(3,0)*clhs304;
const double clhs307 = DN(2,1)*DN(3,1);
const double clhs308 = clhs10*clhs307;
const double clhs309 = clhs294 + clhs296;
const double clhs310 = DN(3,2)*clhs304;
const double clhs311 = DN(2,1)*N[3];
const double clhs312 = DN(3,1)*N[2];
const double clhs313 = DN(2,2)*DN(2,2);
const double clhs314 = DN(2,2)*clhs10;
const double clhs315 = DN(3,0)*clhs314;
const double clhs316 = DN(3,1)*clhs314;
const double clhs317 = DN(2,2)*DN(3,2);
const double clhs318 = clhs10*clhs317;
const double clhs319 = DN(2,2)*N[3];
const double clhs320 = DN(3,2)*N[2];
const double clhs321 = N[2] + clhs217;
const double clhs322 = clhs16*(clhs293 + clhs307 + clhs317);
const double clhs323 = clhs107*clhs17;
const double clhs324 = N[3]*clhs20;
const double clhs325 = clhs104*clhs11 + clhs15*clhs323 + clhs15*clhs324;
const double clhs326 = clhs107*clhs38;
const double clhs327 = clhs104*clhs50 + clhs323*clhs52 + clhs324*clhs52;
const double clhs328 = clhs104*clhs79 + clhs323*clhs81 + clhs324*clhs81;
const double clhs329 = DN(3,0)*DN(3,0);
const double clhs330 = clhs104*clhs107 + clhs108*clhs323 + clhs108*clhs324 + clhs13*(N[3]*N[3]);
const double clhs331 = DN(3,0)*clhs10;
const double clhs332 = DN(3,1)*clhs331;
const double clhs333 = DN(3,2)*clhs331;
const double clhs334 = -N[3] + clhs104*clhs37 + clhs326;
const double clhs335 = clhs105 + clhs325;
const double clhs336 = clhs240 + clhs327;
const double clhs337 = clhs294 + clhs328;
const double clhs338 = DN(3,1)*DN(3,1);
const double clhs339 = DN(3,1)*DN(3,2)*clhs10;
const double clhs340 = DN(3,2)*DN(3,2);
const double clhs341 = N[3] + clhs219;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
lhs(0,3)=DN(0,0)*clhs40;
lhs(0,4)=DN(0,0)*clhs41 + DN(0,1)*clhs43 + DN(0,2)*clhs45 + clhs49 + clhs53;
lhs(0,5)=DN(0,0)*clhs54 + DN(0,1)*clhs56 + DN(0,2)*clhs59 + clhs60;
lhs(0,6)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + DN(0,2)*clhs65 + clhs66;
lhs(0,7)=DN(1,0)*clhs39 - clhs67 + clhs68*clhs69;
lhs(0,8)=DN(0,0)*clhs70 + DN(0,1)*clhs72 + DN(0,2)*clhs74 + clhs78 + clhs82;
lhs(0,9)=DN(0,0)*clhs83 + DN(0,1)*clhs85 + DN(0,2)*clhs88 + clhs89;
lhs(0,10)=DN(0,0)*clhs90 + DN(0,1)*clhs92 + DN(0,2)*clhs94 + clhs95;
lhs(0,11)=DN(2,0)*clhs39 + clhs69*clhs97 - clhs96;
lhs(0,12)=DN(0,0)*clhs98 + DN(0,1)*clhs100 + DN(0,2)*clhs102 + clhs106 + clhs109;
lhs(0,13)=DN(0,0)*clhs110 + DN(0,1)*clhs112 + DN(0,2)*clhs115 + clhs116;
lhs(0,14)=DN(0,0)*clhs117 + DN(0,1)*clhs119 + DN(0,2)*clhs121 + clhs122;
lhs(0,15)=DN(3,0)*clhs39 - clhs123 + clhs124*clhs69;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs125 + DN(0,2)*clhs126 + clhs30;
lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs127 + DN(0,2)*clhs129 + clhs10*clhs130 + clhs22;
lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs135;
lhs(1,3)=DN(0,1)*clhs40;
lhs(1,4)=DN(0,0)*clhs43 + DN(0,1)*clhs136 + DN(0,2)*clhs137 + clhs138;
lhs(1,5)=DN(0,0)*clhs56 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs143 + clhs144;
lhs(1,6)=DN(0,0)*clhs63 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148;
lhs(1,7)=DN(1,1)*clhs39 - clhs149 + clhs150*clhs69;
lhs(1,8)=DN(0,0)*clhs72 + DN(0,1)*clhs151 + DN(0,2)*clhs152 + clhs153;
lhs(1,9)=DN(0,0)*clhs85 + DN(0,1)*clhs154 + DN(0,2)*clhs156 + clhs158 + clhs159;
lhs(1,10)=DN(0,0)*clhs92 + DN(0,1)*clhs160 + DN(0,2)*clhs162 + clhs163;
lhs(1,11)=DN(2,1)*clhs39 - clhs164 + clhs165*clhs69;
lhs(1,12)=DN(0,0)*clhs100 + DN(0,1)*clhs166 + DN(0,2)*clhs167 + clhs168;
lhs(1,13)=DN(0,0)*clhs112 + DN(0,1)*clhs169 + DN(0,2)*clhs171 + clhs173 + clhs174;
lhs(1,14)=DN(0,0)*clhs119 + DN(0,1)*clhs175 + DN(0,2)*clhs177 + clhs178;
lhs(1,15)=DN(3,1)*clhs39 - clhs179 + clhs180*clhs69;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs126 + DN(0,2)*clhs181 + clhs36;
lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs129 + DN(0,2)*clhs182 + clhs135;
lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs133 + DN(0,2)*clhs183 + clhs10*clhs184 + clhs22;
lhs(2,3)=DN(0,2)*clhs40;
lhs(2,4)=DN(0,0)*clhs45 + DN(0,1)*clhs137 + DN(0,2)*clhs185 + clhs187;
lhs(2,5)=DN(0,0)*clhs59 + DN(0,1)*clhs141 + DN(0,2)*clhs188 + clhs189;
lhs(2,6)=DN(0,0)*clhs65 + DN(0,1)*clhs147 + DN(0,2)*clhs190 + clhs144 + clhs192;
lhs(2,7)=DN(1,2)*clhs39 - clhs193 + clhs194*clhs69;
lhs(2,8)=DN(0,0)*clhs74 + DN(0,1)*clhs152 + DN(0,2)*clhs195 + clhs196;
lhs(2,9)=DN(0,0)*clhs88 + DN(0,1)*clhs156 + DN(0,2)*clhs197 + clhs198;
lhs(2,10)=DN(0,0)*clhs94 + DN(0,1)*clhs162 + DN(0,2)*clhs199 + clhs159 + clhs201;
lhs(2,11)=DN(2,2)*clhs39 - clhs202 + clhs203*clhs69;
lhs(2,12)=DN(0,0)*clhs102 + DN(0,1)*clhs167 + DN(0,2)*clhs204 + clhs205;
lhs(2,13)=DN(0,0)*clhs115 + DN(0,1)*clhs171 + DN(0,2)*clhs206 + clhs207;
lhs(2,14)=DN(0,0)*clhs121 + DN(0,1)*clhs177 + DN(0,2)*clhs208 + clhs174 + clhs210;
lhs(2,15)=DN(3,2)*clhs39 - clhs211 + clhs212*clhs69;
lhs(3,0)=DN(0,0)*clhs214;
lhs(3,1)=DN(0,1)*clhs214;
lhs(3,2)=DN(0,2)*clhs214;
lhs(3,3)=clhs16*(clhs130 + clhs184 + clhs5);
lhs(3,4)=DN(0,0)*clhs215 + clhs68;
lhs(3,5)=DN(0,1)*clhs215 + clhs150;
lhs(3,6)=DN(0,2)*clhs215 + clhs194;
lhs(3,7)=clhs216;
lhs(3,8)=DN(0,0)*clhs217 + clhs97;
lhs(3,9)=DN(0,1)*clhs217 + clhs165;
lhs(3,10)=DN(0,2)*clhs217 + clhs203;
lhs(3,11)=clhs218;
lhs(3,12)=DN(0,0)*clhs219 + clhs124;
lhs(3,13)=DN(0,1)*clhs219 + clhs180;
lhs(3,14)=DN(0,2)*clhs219 + clhs212;
lhs(3,15)=clhs220;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs223 + clhs49;
lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs138;
lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs187;
lhs(4,3)=DN(0,0)*clhs224 + clhs67*clhs69 - clhs68;
lhs(4,4)=DN(1,0)*clhs41 + DN(1,1)*clhs43 + DN(1,2)*clhs45 + clhs10*clhs225 + clhs226;
lhs(4,5)=DN(1,0)*clhs54 + DN(1,1)*clhs56 + DN(1,2)*clhs59 + clhs228;
lhs(4,6)=DN(1,0)*clhs61 + DN(1,1)*clhs63 + DN(1,2)*clhs65 + clhs229;
lhs(4,7)=DN(1,0)*clhs230;
lhs(4,8)=DN(1,0)*clhs70 + DN(1,1)*clhs72 + DN(1,2)*clhs74 + clhs233 + clhs234;
lhs(4,9)=DN(1,0)*clhs83 + DN(1,1)*clhs85 + DN(1,2)*clhs88 + clhs235;
lhs(4,10)=DN(1,0)*clhs90 + DN(1,1)*clhs92 + DN(1,2)*clhs94 + clhs236;
lhs(4,11)=DN(2,0)*clhs224 - clhs237 + clhs238*clhs69;
lhs(4,12)=DN(1,0)*clhs98 + DN(1,1)*clhs100 + DN(1,2)*clhs102 + clhs241 + clhs242;
lhs(4,13)=DN(1,0)*clhs110 + DN(1,1)*clhs112 + DN(1,2)*clhs115 + clhs243;
lhs(4,14)=DN(1,0)*clhs117 + DN(1,1)*clhs119 + DN(1,2)*clhs121 + clhs244;
lhs(4,15)=DN(3,0)*clhs224 - clhs245 + clhs246*clhs69;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs125 + DN(1,2)*clhs126 + clhs60;
lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs127 + DN(1,2)*clhs129 + clhs143 + clhs247;
lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs189;
lhs(5,3)=DN(0,1)*clhs224 + clhs149*clhs69 - clhs150;
lhs(5,4)=DN(1,0)*clhs43 + DN(1,1)*clhs136 + DN(1,2)*clhs137 + clhs228;
lhs(5,5)=DN(1,0)*clhs56 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs10*clhs248 + clhs226;
lhs(5,6)=DN(1,0)*clhs63 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs250;
lhs(5,7)=DN(1,1)*clhs230;
lhs(5,8)=DN(1,0)*clhs72 + DN(1,1)*clhs151 + DN(1,2)*clhs152 + clhs251;
lhs(5,9)=DN(1,0)*clhs85 + DN(1,1)*clhs154 + DN(1,2)*clhs156 + clhs253 + clhs254;
lhs(5,10)=DN(1,0)*clhs92 + DN(1,1)*clhs160 + DN(1,2)*clhs162 + clhs255;
lhs(5,11)=DN(2,1)*clhs224 - clhs256 + clhs257*clhs69;
lhs(5,12)=DN(1,0)*clhs100 + DN(1,1)*clhs166 + DN(1,2)*clhs167 + clhs258;
lhs(5,13)=DN(1,0)*clhs112 + DN(1,1)*clhs169 + DN(1,2)*clhs171 + clhs260 + clhs261;
lhs(5,14)=DN(1,0)*clhs119 + DN(1,1)*clhs175 + DN(1,2)*clhs177 + clhs262;
lhs(5,15)=DN(3,1)*clhs224 - clhs263 + clhs264*clhs69;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs126 + DN(1,2)*clhs181 + clhs66;
lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs129 + DN(1,2)*clhs182 + clhs148;
lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs133 + DN(1,2)*clhs183 + clhs192 + clhs247;
lhs(6,3)=DN(0,2)*clhs224 + clhs193*clhs69 - clhs194;
lhs(6,4)=DN(1,0)*clhs45 + DN(1,1)*clhs137 + DN(1,2)*clhs185 + clhs229;
lhs(6,5)=DN(1,0)*clhs59 + DN(1,1)*clhs141 + DN(1,2)*clhs188 + clhs250;
lhs(6,6)=DN(1,0)*clhs65 + DN(1,1)*clhs147 + DN(1,2)*clhs190 + clhs10*clhs265 + clhs226;
lhs(6,7)=DN(1,2)*clhs230;
lhs(6,8)=DN(1,0)*clhs74 + DN(1,1)*clhs152 + DN(1,2)*clhs195 + clhs267;
lhs(6,9)=DN(1,0)*clhs88 + DN(1,1)*clhs156 + DN(1,2)*clhs197 + clhs268;
lhs(6,10)=DN(1,0)*clhs94 + DN(1,1)*clhs162 + DN(1,2)*clhs199 + clhs254 + clhs270;
lhs(6,11)=DN(2,2)*clhs224 - clhs271 + clhs272*clhs69;
lhs(6,12)=DN(1,0)*clhs102 + DN(1,1)*clhs167 + DN(1,2)*clhs204 + clhs273;
lhs(6,13)=DN(1,0)*clhs115 + DN(1,1)*clhs171 + DN(1,2)*clhs206 + clhs274;
lhs(6,14)=DN(1,0)*clhs121 + DN(1,1)*clhs177 + DN(1,2)*clhs208 + clhs261 + clhs276;
lhs(6,15)=DN(3,2)*clhs224 - clhs277 + clhs278*clhs69;
lhs(7,0)=DN(1,0)*clhs213 + clhs67;
lhs(7,1)=DN(1,1)*clhs213 + clhs149;
lhs(7,2)=DN(1,2)*clhs213 + clhs193;
lhs(7,3)=clhs216;
lhs(7,4)=DN(1,0)*clhs279;
lhs(7,5)=DN(1,1)*clhs279;
lhs(7,6)=DN(1,2)*clhs279;
lhs(7,7)=clhs16*(clhs225 + clhs248 + clhs265);
lhs(7,8)=DN(1,0)*clhs217 + clhs238;
lhs(7,9)=DN(1,1)*clhs217 + clhs257;
lhs(7,10)=DN(1,2)*clhs217 + clhs272;
lhs(7,11)=clhs280;
lhs(7,12)=DN(1,0)*clhs219 + clhs246;
lhs(7,13)=DN(1,1)*clhs219 + clhs264;
lhs(7,14)=DN(1,2)*clhs219 + clhs278;
lhs(7,15)=clhs281;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs284 + clhs78;
lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs153;
lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs196;
lhs(8,3)=DN(0,0)*clhs285 + clhs69*clhs96 - clhs97;
lhs(8,4)=DN(2,0)*clhs41 + DN(2,1)*clhs43 + DN(2,2)*clhs45 + clhs233 + clhs286;
lhs(8,5)=DN(2,0)*clhs54 + DN(2,1)*clhs56 + DN(2,2)*clhs59 + clhs251;
lhs(8,6)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + DN(2,2)*clhs65 + clhs267;
lhs(8,7)=DN(1,0)*clhs285 + clhs237*clhs69 - clhs238;
lhs(8,8)=DN(2,0)*clhs70 + DN(2,1)*clhs72 + DN(2,2)*clhs74 + clhs10*clhs287 + clhs288;
lhs(8,9)=DN(2,0)*clhs83 + DN(2,1)*clhs85 + DN(2,2)*clhs88 + clhs290;
lhs(8,10)=DN(2,0)*clhs90 + DN(2,1)*clhs92 + DN(2,2)*clhs94 + clhs291;
lhs(8,11)=DN(2,0)*clhs292;
lhs(8,12)=DN(2,0)*clhs98 + DN(2,1)*clhs100 + DN(2,2)*clhs102 + clhs295 + clhs296;
lhs(8,13)=DN(2,0)*clhs110 + DN(2,1)*clhs112 + DN(2,2)*clhs115 + clhs297;
lhs(8,14)=DN(2,0)*clhs117 + DN(2,1)*clhs119 + DN(2,2)*clhs121 + clhs298;
lhs(8,15)=DN(3,0)*clhs285 - clhs299 + clhs300*clhs69;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs125 + DN(2,2)*clhs126 + clhs89;
lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs127 + DN(2,2)*clhs129 + clhs158 + clhs301;
lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs198;
lhs(9,3)=DN(0,1)*clhs285 + clhs164*clhs69 - clhs165;
lhs(9,4)=DN(2,0)*clhs43 + DN(2,1)*clhs136 + DN(2,2)*clhs137 + clhs235;
lhs(9,5)=DN(2,0)*clhs56 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs253 + clhs302;
lhs(9,6)=DN(2,0)*clhs63 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs268;
lhs(9,7)=DN(1,1)*clhs285 + clhs256*clhs69 - clhs257;
lhs(9,8)=DN(2,0)*clhs72 + DN(2,1)*clhs151 + DN(2,2)*clhs152 + clhs290;
lhs(9,9)=DN(2,0)*clhs85 + DN(2,1)*clhs154 + DN(2,2)*clhs156 + clhs10*clhs303 + clhs288;
lhs(9,10)=DN(2,0)*clhs92 + DN(2,1)*clhs160 + DN(2,2)*clhs162 + clhs305;
lhs(9,11)=DN(2,1)*clhs292;
lhs(9,12)=DN(2,0)*clhs100 + DN(2,1)*clhs166 + DN(2,2)*clhs167 + clhs306;
lhs(9,13)=DN(2,0)*clhs112 + DN(2,1)*clhs169 + DN(2,2)*clhs171 + clhs308 + clhs309;
lhs(9,14)=DN(2,0)*clhs119 + DN(2,1)*clhs175 + DN(2,2)*clhs177 + clhs310;
lhs(9,15)=DN(3,1)*clhs285 - clhs311 + clhs312*clhs69;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs126 + DN(2,2)*clhs181 + clhs95;
lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs129 + DN(2,2)*clhs182 + clhs163;
lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs133 + DN(2,2)*clhs183 + clhs201 + clhs301;
lhs(10,3)=DN(0,2)*clhs285 + clhs202*clhs69 - clhs203;
lhs(10,4)=DN(2,0)*clhs45 + DN(2,1)*clhs137 + DN(2,2)*clhs185 + clhs236;
lhs(10,5)=DN(2,0)*clhs59 + DN(2,1)*clhs141 + DN(2,2)*clhs188 + clhs255;
lhs(10,6)=DN(2,0)*clhs65 + DN(2,1)*clhs147 + DN(2,2)*clhs190 + clhs270 + clhs302;
lhs(10,7)=DN(1,2)*clhs285 + clhs271*clhs69 - clhs272;
lhs(10,8)=DN(2,0)*clhs74 + DN(2,1)*clhs152 + DN(2,2)*clhs195 + clhs291;
lhs(10,9)=DN(2,0)*clhs88 + DN(2,1)*clhs156 + DN(2,2)*clhs197 + clhs305;
lhs(10,10)=DN(2,0)*clhs94 + DN(2,1)*clhs162 + DN(2,2)*clhs199 + clhs10*clhs313 + clhs288;
lhs(10,11)=DN(2,2)*clhs292;
lhs(10,12)=DN(2,0)*clhs102 + DN(2,1)*clhs167 + DN(2,2)*clhs204 + clhs315;
lhs(10,13)=DN(2,0)*clhs115 + DN(2,1)*clhs171 + DN(2,2)*clhs206 + clhs316;
lhs(10,14)=DN(2,0)*clhs121 + DN(2,1)*clhs177 + DN(2,2)*clhs208 + clhs309 + clhs318;
lhs(10,15)=DN(3,2)*clhs285 - clhs319 + clhs320*clhs69;
lhs(11,0)=DN(2,0)*clhs213 + clhs96;
lhs(11,1)=DN(2,1)*clhs213 + clhs164;
lhs(11,2)=DN(2,2)*clhs213 + clhs202;
lhs(11,3)=clhs218;
lhs(11,4)=DN(2,0)*clhs215 + clhs237;
lhs(11,5)=DN(2,1)*clhs215 + clhs256;
lhs(11,6)=DN(2,2)*clhs215 + clhs271;
lhs(11,7)=clhs280;
lhs(11,8)=DN(2,0)*clhs321;
lhs(11,9)=DN(2,1)*clhs321;
lhs(11,10)=DN(2,2)*clhs321;
lhs(11,11)=clhs16*(clhs287 + clhs303 + clhs313);
lhs(11,12)=DN(2,0)*clhs219 + clhs300;
lhs(11,13)=DN(2,1)*clhs219 + clhs312;
lhs(11,14)=DN(2,2)*clhs219 + clhs320;
lhs(11,15)=clhs322;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs106 + clhs325;
lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs168;
lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs205;
lhs(12,3)=DN(0,0)*clhs326 + clhs123*clhs69 - clhs124;
lhs(12,4)=DN(3,0)*clhs41 + DN(3,1)*clhs43 + DN(3,2)*clhs45 + clhs241 + clhs327;
lhs(12,5)=DN(3,0)*clhs54 + DN(3,1)*clhs56 + DN(3,2)*clhs59 + clhs258;
lhs(12,6)=DN(3,0)*clhs61 + DN(3,1)*clhs63 + DN(3,2)*clhs65 + clhs273;
lhs(12,7)=DN(1,0)*clhs326 + clhs245*clhs69 - clhs246;
lhs(12,8)=DN(3,0)*clhs70 + DN(3,1)*clhs72 + DN(3,2)*clhs74 + clhs295 + clhs328;
lhs(12,9)=DN(3,0)*clhs83 + DN(3,1)*clhs85 + DN(3,2)*clhs88 + clhs306;
lhs(12,10)=DN(3,0)*clhs90 + DN(3,1)*clhs92 + DN(3,2)*clhs94 + clhs315;
lhs(12,11)=DN(2,0)*clhs326 + clhs299*clhs69 - clhs300;
lhs(12,12)=DN(3,0)*clhs98 + DN(3,1)*clhs100 + DN(3,2)*clhs102 + clhs10*clhs329 + clhs330;
lhs(12,13)=DN(3,0)*clhs110 + DN(3,1)*clhs112 + DN(3,2)*clhs115 + clhs332;
lhs(12,14)=DN(3,0)*clhs117 + DN(3,1)*clhs119 + DN(3,2)*clhs121 + clhs333;
lhs(12,15)=DN(3,0)*clhs334;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs125 + DN(3,2)*clhs126 + clhs116;
lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs127 + DN(3,2)*clhs129 + clhs173 + clhs335;
lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs207;
lhs(13,3)=DN(0,1)*clhs326 + clhs179*clhs69 - clhs180;
lhs(13,4)=DN(3,0)*clhs43 + DN(3,1)*clhs136 + DN(3,2)*clhs137 + clhs243;
lhs(13,5)=DN(3,0)*clhs56 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs260 + clhs336;
lhs(13,6)=DN(3,0)*clhs63 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs274;
lhs(13,7)=DN(1,1)*clhs326 + clhs263*clhs69 - clhs264;
lhs(13,8)=DN(3,0)*clhs72 + DN(3,1)*clhs151 + DN(3,2)*clhs152 + clhs297;
lhs(13,9)=DN(3,0)*clhs85 + DN(3,1)*clhs154 + DN(3,2)*clhs156 + clhs308 + clhs337;
lhs(13,10)=DN(3,0)*clhs92 + DN(3,1)*clhs160 + DN(3,2)*clhs162 + clhs316;
lhs(13,11)=DN(2,1)*clhs326 + clhs311*clhs69 - clhs312;
lhs(13,12)=DN(3,0)*clhs100 + DN(3,1)*clhs166 + DN(3,2)*clhs167 + clhs332;
lhs(13,13)=DN(3,0)*clhs112 + DN(3,1)*clhs169 + DN(3,2)*clhs171 + clhs10*clhs338 + clhs330;
lhs(13,14)=DN(3,0)*clhs119 + DN(3,1)*clhs175 + DN(3,2)*clhs177 + clhs339;
lhs(13,15)=DN(3,1)*clhs334;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs126 + DN(3,2)*clhs181 + clhs122;
lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs129 + DN(3,2)*clhs182 + clhs178;
lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs133 + DN(3,2)*clhs183 + clhs210 + clhs335;
lhs(14,3)=DN(0,2)*clhs326 + clhs211*clhs69 - clhs212;
lhs(14,4)=DN(3,0)*clhs45 + DN(3,1)*clhs137 + DN(3,2)*clhs185 + clhs244;
lhs(14,5)=DN(3,0)*clhs59 + DN(3,1)*clhs141 + DN(3,2)*clhs188 + clhs262;
lhs(14,6)=DN(3,0)*clhs65 + DN(3,1)*clhs147 + DN(3,2)*clhs190 + clhs276 + clhs336;
lhs(14,7)=DN(1,2)*clhs326 + clhs277*clhs69 - clhs278;
lhs(14,8)=DN(3,0)*clhs74 + DN(3,1)*clhs152 + DN(3,2)*clhs195 + clhs298;
lhs(14,9)=DN(3,0)*clhs88 + DN(3,1)*clhs156 + DN(3,2)*clhs197 + clhs310;
lhs(14,10)=DN(3,0)*clhs94 + DN(3,1)*clhs162 + DN(3,2)*clhs199 + clhs318 + clhs337;
lhs(14,11)=DN(2,2)*clhs326 + clhs319*clhs69 - clhs320;
lhs(14,12)=DN(3,0)*clhs102 + DN(3,1)*clhs167 + DN(3,2)*clhs204 + clhs333;
lhs(14,13)=DN(3,0)*clhs115 + DN(3,1)*clhs171 + DN(3,2)*clhs206 + clhs339;
lhs(14,14)=DN(3,0)*clhs121 + DN(3,1)*clhs177 + DN(3,2)*clhs208 + clhs10*clhs340 + clhs330;
lhs(14,15)=DN(3,2)*clhs334;
lhs(15,0)=DN(3,0)*clhs213 + clhs123;
lhs(15,1)=DN(3,1)*clhs213 + clhs179;
lhs(15,2)=DN(3,2)*clhs213 + clhs211;
lhs(15,3)=clhs220;
lhs(15,4)=DN(3,0)*clhs215 + clhs245;
lhs(15,5)=DN(3,1)*clhs215 + clhs263;
lhs(15,6)=DN(3,2)*clhs215 + clhs277;
lhs(15,7)=clhs281;
lhs(15,8)=DN(3,0)*clhs217 + clhs299;
lhs(15,9)=DN(3,1)*clhs217 + clhs311;
lhs(15,10)=DN(3,2)*clhs217 + clhs319;
lhs(15,11)=clhs322;
lhs(15,12)=DN(3,0)*clhs341;
lhs(15,13)=DN(3,1)*clhs341;
lhs(15,14)=DN(3,2)*clhs341;
lhs(15,15)=clhs16*(clhs329 + clhs338 + clhs340);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesFractionalData<2, 3> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Mass correction term
    double volume_error_ratio = 0.0;

    if (rData.IsCut())
    {
        const auto &phi = rData.Distance;
        const double previous_dt = rData.PreviousDeltaTime;
        double volume_error = 0.0;
        double distance_gauss = 0.0;
        for (std::size_t i = 0; i < N.size(); ++i)
        {
            distance_gauss += N[i] * phi[i];
        }
        if (distance_gauss < 0.0)
        {
            volume_error = -rData.WaterVolumeError;
        }
        else if (distance_gauss > 0.0)
        {
            volume_error = -rData.AirVolumeError;
        }
        volume_error_ratio = volume_error / previous_dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs2 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0));
const double crhs3 = N[0]*rho;
const double crhs4 = bdf0*crhs3;
const double crhs5 = N[0]*(vn(0,0) - vnn(0,0));
const double crhs6 = N[1]*(vn(1,0) - vnn(1,0));
const double crhs7 = N[2]*(vn(2,0) - vnn(2,0));
const double crhs8 = crhs5 + crhs6 + crhs7;
const double crhs9 = 1.0/dt;
const double crhs10 = crhs3*crhs9;
const double crhs11 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs12 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs14 = rho*(crhs11*crhs12 + crhs13*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs15 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crhs16 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crhs17 = crhs15*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crhs16*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0));
const double crhs18 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crhs19 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crhs20 = crhs18*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crhs19*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0));
const double crhs21 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs22 = crhs11 + crhs21 - volume_error_ratio;
const double crhs23 = rho*stab_c2*sqrt(crhs12*crhs12 + crhs13*crhs13);
const double crhs24 = crhs22*(crhs23*h*1.0/stab_c1 + mu);
const double crhs25 = bdf0*rho;
const double crhs26 = crhs2*crhs25;
const double crhs27 = crhs9*rho;
const double crhs28 = 1.0*1.0/(crhs23*1.0/h + crhs27*dyn_tau + mu*stab_c1*1.0/(h*h));
const double crhs29 = crhs28*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1*rho + crhs14 - crhs20*rho + crhs26 + rho*(crhs17 + crhs5*crhs9 + crhs6*crhs9 + crhs7*crhs9));
const double crhs30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs31 = crhs3*crhs30;
const double crhs32 = rho*(DN(0,0)*crhs12 + DN(0,1)*crhs13);
const double crhs33 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs34 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1));
const double crhs35 = N[0]*(vn(0,1) - vnn(0,1));
const double crhs36 = N[1]*(vn(1,1) - vnn(1,1));
const double crhs37 = N[2]*(vn(2,1) - vnn(2,1));
const double crhs38 = crhs35 + crhs36 + crhs37;
const double crhs39 = rho*(crhs12*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs13*crhs21);
const double crhs40 = crhs15*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crhs16*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1));
const double crhs41 = crhs18*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crhs19*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1));
const double crhs42 = crhs25*crhs34;
const double crhs43 = crhs28*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs33*rho + crhs39 - crhs41*rho + crhs42 + rho*(crhs35*crhs9 + crhs36*crhs9 + crhs37*crhs9 + crhs40));
const double crhs44 = N[1]*crhs27;
const double crhs45 = N[1]*rho;
const double crhs46 = crhs30*crhs45;
const double crhs47 = rho*(DN(1,0)*crhs12 + DN(1,1)*crhs13);
const double crhs48 = N[2]*crhs27;
const double crhs49 = N[2]*rho;
const double crhs50 = crhs30*crhs49;
const double crhs51 = rho*(DN(2,0)*crhs12 + DN(2,1)*crhs13);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs24 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1*rho - N[0]*crhs14 + N[0]*crhs20*rho - crhs10*crhs8 - crhs17*crhs3 - crhs2*crhs4 - crhs29*crhs31 - crhs29*crhs32;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs24 - DN(0,1)*stress[1] + N[0]*crhs33*rho - N[0]*crhs39 + N[0]*crhs41*rho - crhs10*crhs38 - crhs3*crhs40 - crhs31*crhs43 - crhs32*crhs43 - crhs34*crhs4;
rhs[2]=-DN(0,0)*crhs29 - DN(0,1)*crhs43 - N[0]*crhs22;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs24 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1*rho - N[1]*crhs14 + N[1]*crhs20*rho - N[1]*crhs26 - crhs17*crhs45 - crhs29*crhs46 - crhs29*crhs47 - crhs44*crhs8;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs24 - DN(1,1)*stress[1] + N[1]*crhs33*rho - N[1]*crhs39 + N[1]*crhs41*rho - N[1]*crhs42 - crhs38*crhs44 - crhs40*crhs45 - crhs43*crhs46 - crhs43*crhs47;
rhs[5]=-DN(1,0)*crhs29 - DN(1,1)*crhs43 - N[1]*crhs22;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs24 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1*rho - N[2]*crhs14 + N[2]*crhs20*rho - N[2]*crhs26 - crhs17*crhs49 - crhs29*crhs50 - crhs29*crhs51 - crhs48*crhs8;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs24 - DN(2,1)*stress[1] + N[2]*crhs33*rho - N[2]*crhs39 + N[2]*crhs41*rho - N[2]*crhs42 - crhs38*crhs48 - crhs40*crhs49 - crhs43*crhs50 - crhs43*crhs51;
rhs[8]=-DN(2,0)*crhs29 - DN(2,1)*crhs43 - N[2]*crhs22;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesFractionalData<3, 4> &rData,
    VectorType &rRHS)
{

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double dyn_tau = rData.DynamicTau;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const auto &stress = rData.ShearStress;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Mass correction term
    double volume_error_ratio = 0.0;

    if (rData.IsCut())
    {
        const auto &phi = rData.Distance;
        const double previous_dt = rData.PreviousDeltaTime;
        double volume_error = 0.0;
        double distance_gauss = 0.0;
        for (std::size_t i = 0; i < N.size(); ++i)
        {
            distance_gauss += N[i] * phi[i];
        }
        if (distance_gauss < 0.0)
        {
            volume_error = -rData.WaterVolumeError;
        }
        else if (distance_gauss > 0.0)
        {
            volume_error = -rData.AirVolumeError;
        }
        volume_error_ratio = volume_error / previous_dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs2 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0)) + N[3]*(v(3,0) - vfrac(3,0));
const double crhs3 = N[0]*rho;
const double crhs4 = bdf0*crhs3;
const double crhs5 = N[0]*(vn(0,0) - vnn(0,0));
const double crhs6 = N[1]*(vn(1,0) - vnn(1,0));
const double crhs7 = N[2]*(vn(2,0) - vnn(2,0));
const double crhs8 = N[3]*(vn(3,0) - vnn(3,0));
const double crhs9 = crhs5 + crhs6 + crhs7 + crhs8;
const double crhs10 = 1.0/dt;
const double crhs11 = crhs10*crhs3;
const double crhs12 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs13 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs14 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs15 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs16 = rho*(crhs12*crhs13 + crhs14*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs15*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs17 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crhs18 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crhs19 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crhs20 = crhs17*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crhs18*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crhs19*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0));
const double crhs21 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crhs22 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crhs23 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crhs24 = crhs21*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crhs22*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0)) + crhs23*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0));
const double crhs25 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs26 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs27 = crhs12 + crhs25 + crhs26 - volume_error_ratio;
const double crhs28 = rho*stab_c2*sqrt(crhs13*crhs13 + crhs14*crhs14 + crhs15*crhs15);
const double crhs29 = crhs27*(crhs28*h*1.0/stab_c1 + mu);
const double crhs30 = bdf0*rho;
const double crhs31 = crhs2*crhs30;
const double crhs32 = crhs10*rho;
const double crhs33 = 1.0*1.0/(crhs28*1.0/h + crhs32*dyn_tau + mu*stab_c1*1.0/(h*h));
const double crhs34 = crhs33*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1*rho + crhs16 - crhs24*rho + crhs31 + rho*(crhs10*crhs5 + crhs10*crhs6 + crhs10*crhs7 + crhs10*crhs8 + crhs20));
const double crhs35 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs36 = crhs3*crhs35;
const double crhs37 = rho*(DN(0,0)*crhs13 + DN(0,1)*crhs14 + DN(0,2)*crhs15);
const double crhs38 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs39 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1)) + N[3]*(v(3,1) - vfrac(3,1));
const double crhs40 = N[0]*(vn(0,1) - vnn(0,1));
const double crhs41 = N[1]*(vn(1,1) - vnn(1,1));
const double crhs42 = N[2]*(vn(2,1) - vnn(2,1));
const double crhs43 = N[3]*(vn(3,1) - vnn(3,1));
const double crhs44 = crhs40 + crhs41 + crhs42 + crhs43;
const double crhs45 = rho*(crhs13*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs14*crhs25 + crhs15*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs46 = crhs17*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crhs18*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crhs19*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1));
const double crhs47 = crhs21*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crhs22*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1)) + crhs23*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1));
const double crhs48 = crhs30*crhs39;
const double crhs49 = crhs33*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs38*rho + crhs45 - crhs47*rho + crhs48 + rho*(crhs10*crhs40 + crhs10*crhs41 + crhs10*crhs42 + crhs10*crhs43 + crhs46));
const double crhs50 = N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2);
const double crhs51 = N[0]*(v(0,2) - vfrac(0,2)) + N[1]*(v(1,2) - vfrac(1,2)) + N[2]*(v(2,2) - vfrac(2,2)) + N[3]*(v(3,2) - vfrac(3,2));
const double crhs52 = N[0]*(vn(0,2) - vnn(0,2));
const double crhs53 = N[1]*(vn(1,2) - vnn(1,2));
const double crhs54 = N[2]*(vn(2,2) - vnn(2,2));
const double crhs55 = N[3]*(vn(3,2) - vnn(3,2));
const double crhs56 = crhs52 + crhs53 + crhs54 + crhs55;
const double crhs57 = rho*(crhs13*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs14*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs15*crhs26);
const double crhs58 = crhs17*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crhs18*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crhs19*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2));
const double crhs59 = crhs21*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crhs22*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2)) + crhs23*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2));
const double crhs60 = crhs30*crhs51;
const double crhs61 = crhs33*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs50*rho + crhs57 - crhs59*rho + crhs60 + rho*(crhs10*crhs52 + crhs10*crhs53 + crhs10*crhs54 + crhs10*crhs55 + crhs58));
const double crhs62 = N[1]*crhs32;
const double crhs63 = N[1]*rho;
const double crhs64 = crhs35*crhs63;
const double crhs65 = rho*(DN(1,0)*crhs13 + DN(1,1)*crhs14 + DN(1,2)*crhs15);
const double crhs66 = N[2]*crhs32;
const double crhs67 = N[2]*rho;
const double crhs68 = crhs35*crhs67;
const double crhs69 = rho*(DN(2,0)*crhs13 + DN(2,1)*crhs14 + DN(2,2)*crhs15);
const double crhs70 = N[3]*crhs32;
const double crhs71 = N[3]*rho;
const double crhs72 = crhs35*crhs71;
const double crhs73 = rho*(DN(3,0)*crhs13 + DN(3,1)*crhs14 + DN(3,2)*crhs15);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs29 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1*rho - N[0]*crhs16 + N[0]*crhs24*rho - crhs11*crhs9 - crhs2*crhs4 - crhs20*crhs3 - crhs34*crhs36 - crhs34*crhs37;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs29 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs38*rho - N[0]*crhs45 + N[0]*crhs47*rho - crhs11*crhs44 - crhs3*crhs46 - crhs36*crhs49 - crhs37*crhs49 - crhs39*crhs4;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs29 - DN(0,2)*stress[2] + N[0]*crhs50*rho - N[0]*crhs57 + N[0]*crhs59*rho - crhs11*crhs56 - crhs3*crhs58 - crhs36*crhs61 - crhs37*crhs61 - crhs4*crhs51;
rhs[3]=-DN(0,0)*crhs34 - DN(0,1)*crhs49 - DN(0,2)*crhs61 - N[0]*crhs27;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs29 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1*rho - N[1]*crhs16 + N[1]*crhs24*rho - N[1]*crhs31 - crhs20*crhs63 - crhs34*crhs64 - crhs34*crhs65 - crhs62*crhs9;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs29 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs38*rho - N[1]*crhs45 + N[1]*crhs47*rho - N[1]*crhs48 - crhs44*crhs62 - crhs46*crhs63 - crhs49*crhs64 - crhs49*crhs65;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs29 - DN(1,2)*stress[2] + N[1]*crhs50*rho - N[1]*crhs57 + N[1]*crhs59*rho - N[1]*crhs60 - crhs56*crhs62 - crhs58*crhs63 - crhs61*crhs64 - crhs61*crhs65;
rhs[7]=-DN(1,0)*crhs34 - DN(1,1)*crhs49 - DN(1,2)*crhs61 - N[1]*crhs27;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs29 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1*rho - N[2]*crhs16 + N[2]*crhs24*rho - N[2]*crhs31 - crhs20*crhs67 - crhs34*crhs68 - crhs34*crhs69 - crhs66*crhs9;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs29 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs38*rho - N[2]*crhs45 + N[2]*crhs47*rho - N[2]*crhs48 - crhs44*crhs66 - crhs46*crhs67 - crhs49*crhs68 - crhs49*crhs69;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs29 - DN(2,2)*stress[2] + N[2]*crhs50*rho - N[2]*crhs57 + N[2]*crhs59*rho - N[2]*crhs60 - crhs56*crhs66 - crhs58*crhs67 - crhs61*crhs68 - crhs61*crhs69;
rhs[11]=-DN(2,0)*crhs34 - DN(2,1)*crhs49 - DN(2,2)*crhs61 - N[2]*crhs27;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs29 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1*rho - N[3]*crhs16 + N[3]*crhs24*rho - N[3]*crhs31 - crhs20*crhs71 - crhs34*crhs72 - crhs34*crhs73 - crhs70*crhs9;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs29 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs38*rho - N[3]*crhs45 + N[3]*crhs47*rho - N[3]*crhs48 - crhs44*crhs70 - crhs46*crhs71 - crhs49*crhs72 - crhs49*crhs73;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs29 - DN(3,2)*stress[2] + N[3]*crhs50*rho - N[3]*crhs57 + N[3]*crhs59*rho - N[3]*crhs60 - crhs56*crhs70 - crhs58*crhs71 - crhs61*crhs72 - crhs61*crhs73;
rhs[15]=-DN(3,0)*crhs34 - DN(3,1)*crhs49 - DN(3,2)*crhs61 - N[3]*crhs27;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesFractionalData<2, 3> &rData,
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
    const double dyn_tau = rData.DynamicTau;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    // const auto &vnnn = rData.Velocity_OldStep3; # an bdf2 possible
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Mass correction term
    double volume_error_ratio = 0.0;

    if (rData.IsCut())
    {
        const auto &phi = rData.Distance;
        const double previous_dt = rData.PreviousDeltaTime;
        double volume_error = 0.0;
        double distance_gauss = 0.0;
        for (std::size_t i = 0; i < N.size(); ++i)
        {
            distance_gauss += N[i] * phi[i];
        }
        if (distance_gauss < 0.0)
        {
            volume_error = -rData.WaterVolumeError;
        }
        else if (distance_gauss > 0.0)
        {
            volume_error = -rData.AirVolumeError;
        }
        volume_error_ratio = volume_error / previous_dt;
    }

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(cV0*cV0 + cV1*cV1));
const double cV3 = cV2*rho;
const double cV4 = cV3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV5 = N[0]*cV4;
const double cV6 = cV3*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV7 = N[1]*cV4;
const double cV8 = cV3*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV9 = N[2]*cV4;
const double cV10 = cV3*(DN(2,0)*cV0 + DN(2,1)*cV1);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV5 + DNenr(0,0)*cV6;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV5 + DNenr(1,0)*cV6;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV5 + DNenr(2,0)*cV6;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV5 + DNenr(0,1)*cV6;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV5 + DNenr(1,1)*cV6;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV5 + DNenr(2,1)*cV6;
V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV7 + DNenr(0,0)*cV8;
V(3,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV7 + DNenr(1,0)*cV8;
V(3,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV7 + DNenr(2,0)*cV8;
V(4,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV7 + DNenr(0,1)*cV8;
V(4,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV7 + DNenr(1,1)*cV8;
V(4,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV7 + DNenr(2,1)*cV8;
V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV10 + DNenr(0,0)*cV9;
V(6,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV9;
V(6,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV9;
V(7,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV9;
V(7,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV9;
V(7,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV9;
V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(cH0*cH0 + cH1*cH1));
const double cH3 = cH2*rho;
const double cH4 = cH3*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0);
const double cH5 = cH3*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0);
const double cH6 = cH3*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0);
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH4;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH4;
H(0,2)=cH2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
H(0,3)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH5;
H(0,4)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH5;
H(0,5)=cH2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
H(0,6)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH6;
H(0,7)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH6;
H(0,8)=cH2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH4;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH4;
H(1,2)=cH2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
H(1,3)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH5;
H(1,4)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH5;
H(1,5)=cH2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
H(1,6)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH6;
H(1,7)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH6;
H(1,8)=cH2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH4;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH4;
H(2,2)=cH2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
H(2,3)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH5;
H(2,4)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH5;
H(2,5)=cH2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
H(2,6)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH6;
H(2,7)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH6;
H(2,8)=cH2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))*1.0/h);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
Kee(0,0)=cKee0*(DNenr(0,0)*DNenr(0,0) + DNenr(0,1)*DNenr(0,1));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(DNenr(1,0)*DNenr(1,0) + DNenr(1,1)*DNenr(1,1));
Kee(1,2)=cKee3;
Kee(2,0)=cKee2;
Kee(2,1)=cKee3;
Kee(2,2)=cKee0*(DNenr(2,0)*DNenr(2,0) + DNenr(2,1)*DNenr(2,1));


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs_ee2 = crhs_ee0 + crhs_ee1 - volume_error_ratio;
const double crhs_ee3 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crhs_ee4 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crhs_ee5 = N[0]*bdf0;
const double crhs_ee6 = N[1]*bdf0;
const double crhs_ee7 = N[2]*bdf0;
const double crhs_ee8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee10 = 1.0/dt;
const double crhs_ee11 = N[0]*crhs_ee10;
const double crhs_ee12 = N[1]*crhs_ee10;
const double crhs_ee13 = N[2]*crhs_ee10;
const double crhs_ee14 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crhs_ee15 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crhs_ee16 = 1.0*1.0/(crhs_ee10*dyn_tau*rho + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crhs_ee8*crhs_ee8 + crhs_ee9*crhs_ee9));
const double crhs_ee17 = crhs_ee16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(crhs_ee3*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crhs_ee4*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(crhs_ee0*crhs_ee8 + crhs_ee5*(v(0,0) - vfrac(0,0)) + crhs_ee6*(v(1,0) - vfrac(1,0)) + crhs_ee7*(v(2,0) - vfrac(2,0)) + crhs_ee9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))) + rho*(crhs_ee11*(vn(0,0) - vnn(0,0)) + crhs_ee12*(vn(1,0) - vnn(1,0)) + crhs_ee13*(vn(2,0) - vnn(2,0)) + crhs_ee14*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crhs_ee15*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0))));
const double crhs_ee18 = crhs_ee16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(crhs_ee3*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crhs_ee4*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(crhs_ee1*crhs_ee9 + crhs_ee5*(v(0,1) - vfrac(0,1)) + crhs_ee6*(v(1,1) - vfrac(1,1)) + crhs_ee7*(v(2,1) - vfrac(2,1)) + crhs_ee8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))) + rho*(crhs_ee11*(vn(0,1) - vnn(0,1)) + crhs_ee12*(vn(1,1) - vnn(1,1)) + crhs_ee13*(vn(2,1) - vnn(2,1)) + crhs_ee14*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crhs_ee15*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee17 - DNenr(0,1)*crhs_ee18 - Nenr[0]*crhs_ee2;
rhs_ee[1]=-DNenr(1,0)*crhs_ee17 - DNenr(1,1)*crhs_ee18 - Nenr[1]*crhs_ee2;
rhs_ee[2]=-DNenr(2,0)*crhs_ee17 - DNenr(2,1)*crhs_ee18 - Nenr[2]*crhs_ee2;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesFractionalData<3, 4> &rData,
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
    const double dyn_tau = rData.DynamicTau;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    // Mass correction term
    double volume_error_ratio = 0.0;

    if (rData.IsCut())
    {
        const double previous_dt = rData.PreviousDeltaTime;
        const auto &phi = rData.Distance;
        double volume_error = 0.0;
        double distance_gauss = 0.0;
        for (std::size_t i = 0; i < N.size(); ++i)
        {
            distance_gauss += N[i] * phi[i];
        }
        if (distance_gauss < 0.0)
        {
            volume_error = -rData.WaterVolumeError;
        }
        else if (distance_gauss > 0.0)
        {
            volume_error = -rData.AirVolumeError;
        }
        volume_error_ratio = volume_error / previous_dt;
    }

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(cV0*cV0 + cV1*cV1 + cV2*cV2));
const double cV4 = cV3*rho;
const double cV5 = cV4*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV6 = N[0]*cV5;
const double cV7 = cV4*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV8 = N[1]*cV5;
const double cV9 = cV4*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV10 = N[2]*cV5;
const double cV11 = cV4*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV12 = N[3]*cV5;
const double cV13 = cV4*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV6 + DNenr(0,0)*cV7;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV6 + DNenr(1,0)*cV7;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV6 + DNenr(2,0)*cV7;
V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV6 + DNenr(3,0)*cV7;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV6 + DNenr(0,1)*cV7;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV6 + DNenr(1,1)*cV7;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV6 + DNenr(2,1)*cV7;
V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV6 + DNenr(3,1)*cV7;
V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV6 + DNenr(0,2)*cV7;
V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV6 + DNenr(1,2)*cV7;
V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV6 + DNenr(2,2)*cV7;
V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV6 + DNenr(3,2)*cV7;
V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
V(4,0)=-DN(1,0)*Nenr[0] + DNenr(0,0)*cV8 + DNenr(0,0)*cV9;
V(4,1)=-DN(1,0)*Nenr[1] + DNenr(1,0)*cV8 + DNenr(1,0)*cV9;
V(4,2)=-DN(1,0)*Nenr[2] + DNenr(2,0)*cV8 + DNenr(2,0)*cV9;
V(4,3)=-DN(1,0)*Nenr[3] + DNenr(3,0)*cV8 + DNenr(3,0)*cV9;
V(5,0)=-DN(1,1)*Nenr[0] + DNenr(0,1)*cV8 + DNenr(0,1)*cV9;
V(5,1)=-DN(1,1)*Nenr[1] + DNenr(1,1)*cV8 + DNenr(1,1)*cV9;
V(5,2)=-DN(1,1)*Nenr[2] + DNenr(2,1)*cV8 + DNenr(2,1)*cV9;
V(5,3)=-DN(1,1)*Nenr[3] + DNenr(3,1)*cV8 + DNenr(3,1)*cV9;
V(6,0)=-DN(1,2)*Nenr[0] + DNenr(0,2)*cV8 + DNenr(0,2)*cV9;
V(6,1)=-DN(1,2)*Nenr[1] + DNenr(1,2)*cV8 + DNenr(1,2)*cV9;
V(6,2)=-DN(1,2)*Nenr[2] + DNenr(2,2)*cV8 + DNenr(2,2)*cV9;
V(6,3)=-DN(1,2)*Nenr[3] + DNenr(3,2)*cV8 + DNenr(3,2)*cV9;
V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
V(8,0)=-DN(2,0)*Nenr[0] + DNenr(0,0)*cV10 + DNenr(0,0)*cV11;
V(8,1)=-DN(2,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV11;
V(8,2)=-DN(2,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV11;
V(8,3)=-DN(2,0)*Nenr[3] + DNenr(3,0)*cV10 + DNenr(3,0)*cV11;
V(9,0)=-DN(2,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV11;
V(9,1)=-DN(2,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV11;
V(9,2)=-DN(2,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV11;
V(9,3)=-DN(2,1)*Nenr[3] + DNenr(3,1)*cV10 + DNenr(3,1)*cV11;
V(10,0)=-DN(2,2)*Nenr[0] + DNenr(0,2)*cV10 + DNenr(0,2)*cV11;
V(10,1)=-DN(2,2)*Nenr[1] + DNenr(1,2)*cV10 + DNenr(1,2)*cV11;
V(10,2)=-DN(2,2)*Nenr[2] + DNenr(2,2)*cV10 + DNenr(2,2)*cV11;
V(10,3)=-DN(2,2)*Nenr[3] + DNenr(3,2)*cV10 + DNenr(3,2)*cV11;
V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
V(12,0)=-DN(3,0)*Nenr[0] + DNenr(0,0)*cV12 + DNenr(0,0)*cV13;
V(12,1)=-DN(3,0)*Nenr[1] + DNenr(1,0)*cV12 + DNenr(1,0)*cV13;
V(12,2)=-DN(3,0)*Nenr[2] + DNenr(2,0)*cV12 + DNenr(2,0)*cV13;
V(12,3)=-DN(3,0)*Nenr[3] + DNenr(3,0)*cV12 + DNenr(3,0)*cV13;
V(13,0)=-DN(3,1)*Nenr[0] + DNenr(0,1)*cV12 + DNenr(0,1)*cV13;
V(13,1)=-DN(3,1)*Nenr[1] + DNenr(1,1)*cV12 + DNenr(1,1)*cV13;
V(13,2)=-DN(3,1)*Nenr[2] + DNenr(2,1)*cV12 + DNenr(2,1)*cV13;
V(13,3)=-DN(3,1)*Nenr[3] + DNenr(3,1)*cV12 + DNenr(3,1)*cV13;
V(14,0)=-DN(3,2)*Nenr[0] + DNenr(0,2)*cV12 + DNenr(0,2)*cV13;
V(14,1)=-DN(3,2)*Nenr[1] + DNenr(1,2)*cV12 + DNenr(1,2)*cV13;
V(14,2)=-DN(3,2)*Nenr[2] + DNenr(2,2)*cV12 + DNenr(2,2)*cV13;
V(14,3)=-DN(3,2)*Nenr[3] + DNenr(3,2)*cV12 + DNenr(3,2)*cV13;
V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(cH0*cH0 + cH1*cH1 + cH2*cH2));
const double cH4 = cH3*rho;
const double cH5 = cH4*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0);
const double cH6 = cH4*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0);
const double cH7 = cH4*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0);
const double cH8 = cH4*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0);
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH5;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH5;
H(0,2)=DN(0,2)*Nenr[0] + DNenr(0,2)*cH5;
H(0,3)=cH3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
H(0,4)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH6;
H(0,5)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH6;
H(0,6)=DN(1,2)*Nenr[0] + DNenr(0,2)*cH6;
H(0,7)=cH3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
H(0,8)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH7;
H(0,9)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH7;
H(0,10)=DN(2,2)*Nenr[0] + DNenr(0,2)*cH7;
H(0,11)=cH3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
H(0,12)=DN(3,0)*Nenr[0] + DNenr(0,0)*cH8;
H(0,13)=DN(3,1)*Nenr[0] + DNenr(0,1)*cH8;
H(0,14)=DN(3,2)*Nenr[0] + DNenr(0,2)*cH8;
H(0,15)=cH3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH5;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH5;
H(1,2)=DN(0,2)*Nenr[1] + DNenr(1,2)*cH5;
H(1,3)=cH3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
H(1,4)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH6;
H(1,5)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH6;
H(1,6)=DN(1,2)*Nenr[1] + DNenr(1,2)*cH6;
H(1,7)=cH3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
H(1,8)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH7;
H(1,9)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH7;
H(1,10)=DN(2,2)*Nenr[1] + DNenr(1,2)*cH7;
H(1,11)=cH3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
H(1,12)=DN(3,0)*Nenr[1] + DNenr(1,0)*cH8;
H(1,13)=DN(3,1)*Nenr[1] + DNenr(1,1)*cH8;
H(1,14)=DN(3,2)*Nenr[1] + DNenr(1,2)*cH8;
H(1,15)=cH3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH5;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH5;
H(2,2)=DN(0,2)*Nenr[2] + DNenr(2,2)*cH5;
H(2,3)=cH3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
H(2,4)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH6;
H(2,5)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH6;
H(2,6)=DN(1,2)*Nenr[2] + DNenr(2,2)*cH6;
H(2,7)=cH3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
H(2,8)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH7;
H(2,9)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH7;
H(2,10)=DN(2,2)*Nenr[2] + DNenr(2,2)*cH7;
H(2,11)=cH3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
H(2,12)=DN(3,0)*Nenr[2] + DNenr(2,0)*cH8;
H(2,13)=DN(3,1)*Nenr[2] + DNenr(2,1)*cH8;
H(2,14)=DN(3,2)*Nenr[2] + DNenr(2,2)*cH8;
H(2,15)=cH3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
H(3,0)=DN(0,0)*Nenr[3] + DNenr(3,0)*cH5;
H(3,1)=DN(0,1)*Nenr[3] + DNenr(3,1)*cH5;
H(3,2)=DN(0,2)*Nenr[3] + DNenr(3,2)*cH5;
H(3,3)=cH3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
H(3,4)=DN(1,0)*Nenr[3] + DNenr(3,0)*cH6;
H(3,5)=DN(1,1)*Nenr[3] + DNenr(3,1)*cH6;
H(3,6)=DN(1,2)*Nenr[3] + DNenr(3,2)*cH6;
H(3,7)=cH3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
H(3,8)=DN(2,0)*Nenr[3] + DNenr(3,0)*cH7;
H(3,9)=DN(2,1)*Nenr[3] + DNenr(3,1)*cH7;
H(3,10)=DN(2,2)*Nenr[3] + DNenr(3,2)*cH7;
H(3,11)=cH3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
H(3,12)=DN(3,0)*Nenr[3] + DNenr(3,0)*cH8;
H(3,13)=DN(3,1)*Nenr[3] + DNenr(3,1)*cH8;
H(3,14)=DN(3,2)*Nenr[3] + DNenr(3,2)*cH8;
H(3,15)=cH3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))*1.0/h);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 = cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 = cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 = cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
Kee(0,0)=cKee0*(DNenr(0,0)*DNenr(0,0) + DNenr(0,1)*DNenr(0,1) + DNenr(0,2)*DNenr(0,2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(0,3)=cKee3;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(DNenr(1,0)*DNenr(1,0) + DNenr(1,1)*DNenr(1,1) + DNenr(1,2)*DNenr(1,2));
Kee(1,2)=cKee4;
Kee(1,3)=cKee5;
Kee(2,0)=cKee2;
Kee(2,1)=cKee4;
Kee(2,2)=cKee0*(DNenr(2,0)*DNenr(2,0) + DNenr(2,1)*DNenr(2,1) + DNenr(2,2)*DNenr(2,2));
Kee(2,3)=cKee6;
Kee(3,0)=cKee3;
Kee(3,1)=cKee5;
Kee(3,2)=cKee6;
Kee(3,3)=cKee0*(DNenr(3,0)*DNenr(3,0) + DNenr(3,1)*DNenr(3,1) + DNenr(3,2)*DNenr(3,2));


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs_ee2 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee3 = crhs_ee0 + crhs_ee1 + crhs_ee2 - volume_error_ratio;
const double crhs_ee4 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crhs_ee5 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crhs_ee6 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crhs_ee7 = N[0]*bdf0;
const double crhs_ee8 = N[1]*bdf0;
const double crhs_ee9 = N[2]*bdf0;
const double crhs_ee10 = N[3]*bdf0;
const double crhs_ee11 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee12 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee13 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee14 = 1.0/dt;
const double crhs_ee15 = N[0]*crhs_ee14;
const double crhs_ee16 = N[1]*crhs_ee14;
const double crhs_ee17 = N[2]*crhs_ee14;
const double crhs_ee18 = N[3]*crhs_ee14;
const double crhs_ee19 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crhs_ee20 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crhs_ee21 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crhs_ee22 = 1.0*1.0/(crhs_ee14*dyn_tau*rho + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crhs_ee11*crhs_ee11 + crhs_ee12*crhs_ee12 + crhs_ee13*crhs_ee13));
const double crhs_ee23 = crhs_ee22*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(crhs_ee4*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crhs_ee5*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0)) + crhs_ee6*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(crhs_ee0*crhs_ee11 + crhs_ee10*(v(3,0) - vfrac(3,0)) + crhs_ee12*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee13*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs_ee7*(v(0,0) - vfrac(0,0)) + crhs_ee8*(v(1,0) - vfrac(1,0)) + crhs_ee9*(v(2,0) - vfrac(2,0))) + rho*(crhs_ee15*(vn(0,0) - vnn(0,0)) + crhs_ee16*(vn(1,0) - vnn(1,0)) + crhs_ee17*(vn(2,0) - vnn(2,0)) + crhs_ee18*(vn(3,0) - vnn(3,0)) + crhs_ee19*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crhs_ee20*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crhs_ee21*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0))));
const double crhs_ee24 = crhs_ee22*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(crhs_ee4*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crhs_ee5*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1)) + crhs_ee6*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(crhs_ee1*crhs_ee12 + crhs_ee10*(v(3,1) - vfrac(3,1)) + crhs_ee11*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee13*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs_ee7*(v(0,1) - vfrac(0,1)) + crhs_ee8*(v(1,1) - vfrac(1,1)) + crhs_ee9*(v(2,1) - vfrac(2,1))) + rho*(crhs_ee15*(vn(0,1) - vnn(0,1)) + crhs_ee16*(vn(1,1) - vnn(1,1)) + crhs_ee17*(vn(2,1) - vnn(2,1)) + crhs_ee18*(vn(3,1) - vnn(3,1)) + crhs_ee19*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crhs_ee20*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crhs_ee21*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1))));
const double crhs_ee25 = crhs_ee22*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(crhs_ee4*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crhs_ee5*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2)) + crhs_ee6*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2))) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(crhs_ee10*(v(3,2) - vfrac(3,2)) + crhs_ee11*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee12*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs_ee13*crhs_ee2 + crhs_ee7*(v(0,2) - vfrac(0,2)) + crhs_ee8*(v(1,2) - vfrac(1,2)) + crhs_ee9*(v(2,2) - vfrac(2,2))) + rho*(crhs_ee15*(vn(0,2) - vnn(0,2)) + crhs_ee16*(vn(1,2) - vnn(1,2)) + crhs_ee17*(vn(2,2) - vnn(2,2)) + crhs_ee18*(vn(3,2) - vnn(3,2)) + crhs_ee19*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crhs_ee20*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crhs_ee21*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee23 - DNenr(0,1)*crhs_ee24 - DNenr(0,2)*crhs_ee25 - Nenr[0]*crhs_ee3;
rhs_ee[1]=-DNenr(1,0)*crhs_ee23 - DNenr(1,1)*crhs_ee24 - DNenr(1,2)*crhs_ee25 - Nenr[1]*crhs_ee3;
rhs_ee[2]=-DNenr(2,0)*crhs_ee23 - DNenr(2,1)*crhs_ee24 - DNenr(2,2)*crhs_ee25 - Nenr[2]*crhs_ee3;
rhs_ee[3]=-DNenr(3,0)*crhs_ee23 - DNenr(3,1)*crhs_ee24 - DNenr(3,2)*crhs_ee25 - Nenr[3]*crhs_ee3;


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::ComputeSplitting(
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
void TwoFluidNavierStokesFractional<TElementData>::ComputeSplitInterface(
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
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesFractional< TwoFluidNavierStokesFractionalData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokesFractional< TwoFluidNavierStokesFractionalData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::CalculateStrainRate(TElementData& rData) const
{
    FluidElement<TElementData>::CalculateStrainRate(rData);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::CondenseEnrichmentWithContinuity(
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
void TwoFluidNavierStokesFractional<TElementData>::CondenseEnrichment(
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
void TwoFluidNavierStokesFractional<TElementData>::save(Serializer &rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::load(Serializer &rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}


template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::CalculateOnIntegrationPoints(
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

template class TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>;
template class TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>;

} // namespace Kratos