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

    // Add LHS Gauss point contribution
    const double w_gauss = rData.Weight;

    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = DN(0,0)*DN(0,0);
const double crLHS4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS6 = rho*stab_c2*sqrt(crLHS4*crLHS4 + crLHS5*crLHS5);
const double crLHS7 = crLHS6*h*1.0/stab_c1 + mu;
const double crLHS8 = DN(0,0)*crLHS4 + DN(0,1)*crLHS5;
const double crLHS9 = N[0]*rho;
const double crLHS10 = bdf0*rho;
const double crLHS11 = N[0]*bdf0;
const double crLHS12 = crLHS11 + crLHS8;
const double crLHS13 = 1.0*1.0/(crLHS6*1.0/h + dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h));
const double crLHS14 = crLHS13*(rho*rho);
const double crLHS15 = crLHS14*crLHS8;
const double crLHS16 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS17 = crLHS14*crLHS16;
const double crLHS18 = N[0]*crLHS17;
const double crLHS19 = crLHS10*(N[0]*N[0]) + crLHS12*crLHS15 + crLHS12*crLHS18 + crLHS8*crLHS9;
const double crLHS20 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS21 = C(1,2)*DN(0,1);
const double crLHS22 = C(2,2)*DN(0,0) + crLHS21;
const double crLHS23 = DN(0,0)*crLHS7;
const double crLHS24 = DN(0,1)*crLHS23;
const double crLHS25 = crLHS13*crLHS16;
const double crLHS26 = crLHS13*rho;
const double crLHS27 = crLHS26*crLHS8;
const double crLHS28 = w_gauss*(-N[0] + crLHS25*crLHS9 + crLHS27);
const double crLHS29 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS30 = C(0,2)*DN(1,0);
const double crLHS31 = C(2,2)*DN(1,1) + crLHS30;
const double crLHS32 = DN(0,0)*DN(1,0);
const double crLHS33 = N[1]*rho;
const double crLHS34 = crLHS11*crLHS33;
const double crLHS35 = crLHS32*crLHS7 + crLHS34;
const double crLHS36 = DN(1,0)*crLHS4 + DN(1,1)*crLHS5;
const double crLHS37 = N[1]*bdf0;
const double crLHS38 = crLHS36 + crLHS37;
const double crLHS39 = crLHS15*crLHS38 + crLHS18*crLHS38 + crLHS36*crLHS9;
const double crLHS40 = C(0,1)*DN(1,1) + crLHS30;
const double crLHS41 = C(1,2)*DN(1,1);
const double crLHS42 = C(2,2)*DN(1,0) + crLHS41;
const double crLHS43 = DN(1,1)*crLHS23;
const double crLHS44 = DN(0,0)*N[1];
const double crLHS45 = DN(1,0)*N[0];
const double crLHS46 = crLHS16*crLHS26;
const double crLHS47 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS48 = C(0,2)*DN(2,0);
const double crLHS49 = C(2,2)*DN(2,1) + crLHS48;
const double crLHS50 = DN(0,0)*DN(2,0);
const double crLHS51 = N[2]*rho;
const double crLHS52 = crLHS11*crLHS51;
const double crLHS53 = crLHS50*crLHS7 + crLHS52;
const double crLHS54 = DN(2,0)*crLHS4 + DN(2,1)*crLHS5;
const double crLHS55 = N[2]*bdf0 + crLHS54;
const double crLHS56 = crLHS15*crLHS55 + crLHS18*crLHS55 + crLHS54*crLHS9;
const double crLHS57 = C(0,1)*DN(2,1) + crLHS48;
const double crLHS58 = C(1,2)*DN(2,1);
const double crLHS59 = C(2,2)*DN(2,0) + crLHS58;
const double crLHS60 = DN(2,1)*crLHS23;
const double crLHS61 = DN(0,0)*N[2];
const double crLHS62 = DN(2,0)*N[0];
const double crLHS63 = C(0,1)*DN(0,0) + crLHS21;
const double crLHS64 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS65 = DN(0,1)*DN(0,1);
const double crLHS66 = C(0,1)*DN(1,0) + crLHS41;
const double crLHS67 = DN(0,1)*crLHS7;
const double crLHS68 = DN(1,0)*crLHS67;
const double crLHS69 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS70 = DN(0,1)*DN(1,1);
const double crLHS71 = crLHS34 + crLHS7*crLHS70;
const double crLHS72 = DN(0,1)*N[1];
const double crLHS73 = DN(1,1)*N[0];
const double crLHS74 = C(0,1)*DN(2,0) + crLHS58;
const double crLHS75 = DN(2,0)*crLHS67;
const double crLHS76 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS77 = DN(0,1)*DN(2,1);
const double crLHS78 = crLHS52 + crLHS7*crLHS77;
const double crLHS79 = DN(0,1)*N[2];
const double crLHS80 = DN(2,1)*N[0];
const double crLHS81 = crLHS12*crLHS26;
const double crLHS82 = w_gauss*(N[0] + crLHS81);
const double crLHS83 = crLHS13*w_gauss;
const double crLHS84 = crLHS26*crLHS38;
const double crLHS85 = crLHS83*(crLHS32 + crLHS70);
const double crLHS86 = crLHS26*crLHS55;
const double crLHS87 = crLHS83*(crLHS50 + crLHS77);
const double crLHS88 = crLHS14*crLHS36;
const double crLHS89 = N[1]*crLHS17;
const double crLHS90 = crLHS12*crLHS88 + crLHS12*crLHS89 + crLHS33*crLHS8;
const double crLHS91 = crLHS26*crLHS36;
const double crLHS92 = DN(1,0)*DN(1,0);
const double crLHS93 = crLHS10*(N[1]*N[1]) + crLHS33*crLHS36 + crLHS38*crLHS88 + crLHS38*crLHS89;
const double crLHS94 = DN(1,0)*crLHS7;
const double crLHS95 = DN(1,1)*crLHS94;
const double crLHS96 = w_gauss*(-N[1] + crLHS25*crLHS33 + crLHS91);
const double crLHS97 = DN(1,0)*DN(2,0);
const double crLHS98 = crLHS37*crLHS51;
const double crLHS99 = crLHS7*crLHS97 + crLHS98;
const double crLHS100 = crLHS33*crLHS54 + crLHS55*crLHS88 + crLHS55*crLHS89;
const double crLHS101 = DN(2,1)*crLHS94;
const double crLHS102 = DN(1,0)*N[2];
const double crLHS103 = DN(2,0)*N[1];
const double crLHS104 = DN(1,1)*DN(1,1);
const double crLHS105 = DN(2,0)*crLHS7;
const double crLHS106 = DN(1,1)*crLHS105;
const double crLHS107 = DN(1,1)*DN(2,1);
const double crLHS108 = crLHS107*crLHS7 + crLHS98;
const double crLHS109 = DN(1,1)*N[2];
const double crLHS110 = DN(2,1)*N[1];
const double crLHS111 = w_gauss*(N[1] + crLHS84);
const double crLHS112 = crLHS83*(crLHS107 + crLHS97);
const double crLHS113 = crLHS14*crLHS54;
const double crLHS114 = N[2]*crLHS17;
const double crLHS115 = crLHS113*crLHS12 + crLHS114*crLHS12 + crLHS51*crLHS8;
const double crLHS116 = crLHS26*crLHS54;
const double crLHS117 = crLHS113*crLHS38 + crLHS114*crLHS38 + crLHS36*crLHS51;
const double crLHS118 = DN(2,0)*DN(2,0);
const double crLHS119 = crLHS10*(N[2]*N[2]) + crLHS113*crLHS55 + crLHS114*crLHS55 + crLHS51*crLHS54;
const double crLHS120 = DN(2,1)*crLHS105;
const double crLHS121 = w_gauss*(-N[2] + crLHS116 + crLHS25*crLHS51);
const double crLHS122 = DN(2,1)*DN(2,1);
const double crLHS123 = w_gauss*(N[2] + crLHS86);
rLHS(0,0)+=w_gauss*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS19 + crLHS3*crLHS7);
rLHS(0,1)+=w_gauss*(DN(0,0)*crLHS20 + DN(0,1)*crLHS22 + crLHS24);
rLHS(0,2)+=DN(0,0)*crLHS28;
rLHS(0,3)+=w_gauss*(DN(0,0)*crLHS29 + DN(0,1)*crLHS31 + crLHS35 + crLHS39);
rLHS(0,4)+=w_gauss*(DN(0,0)*crLHS40 + DN(0,1)*crLHS42 + crLHS43);
rLHS(0,5)+=w_gauss*(DN(1,0)*crLHS27 - crLHS44 + crLHS45*crLHS46);
rLHS(0,6)+=w_gauss*(DN(0,0)*crLHS47 + DN(0,1)*crLHS49 + crLHS53 + crLHS56);
rLHS(0,7)+=w_gauss*(DN(0,0)*crLHS57 + DN(0,1)*crLHS59 + crLHS60);
rLHS(0,8)+=w_gauss*(DN(2,0)*crLHS27 + crLHS46*crLHS62 - crLHS61);
rLHS(1,0)+=w_gauss*(DN(0,0)*crLHS2 + DN(0,1)*crLHS63 + crLHS24);
rLHS(1,1)+=w_gauss*(DN(0,0)*crLHS22 + DN(0,1)*crLHS64 + crLHS19 + crLHS65*crLHS7);
rLHS(1,2)+=DN(0,1)*crLHS28;
rLHS(1,3)+=w_gauss*(DN(0,0)*crLHS31 + DN(0,1)*crLHS66 + crLHS68);
rLHS(1,4)+=w_gauss*(DN(0,0)*crLHS42 + DN(0,1)*crLHS69 + crLHS39 + crLHS71);
rLHS(1,5)+=w_gauss*(DN(1,1)*crLHS27 + crLHS46*crLHS73 - crLHS72);
rLHS(1,6)+=w_gauss*(DN(0,0)*crLHS49 + DN(0,1)*crLHS74 + crLHS75);
rLHS(1,7)+=w_gauss*(DN(0,0)*crLHS59 + DN(0,1)*crLHS76 + crLHS56 + crLHS78);
rLHS(1,8)+=w_gauss*(DN(2,1)*crLHS27 + crLHS46*crLHS80 - crLHS79);
rLHS(2,0)+=DN(0,0)*crLHS82;
rLHS(2,1)+=DN(0,1)*crLHS82;
rLHS(2,2)+=crLHS83*(crLHS3 + crLHS65);
rLHS(2,3)+=w_gauss*(DN(0,0)*crLHS84 + crLHS45);
rLHS(2,4)+=w_gauss*(DN(0,1)*crLHS84 + crLHS73);
rLHS(2,5)+=crLHS85;
rLHS(2,6)+=w_gauss*(DN(0,0)*crLHS86 + crLHS62);
rLHS(2,7)+=w_gauss*(DN(0,1)*crLHS86 + crLHS80);
rLHS(2,8)+=crLHS87;
rLHS(3,0)+=w_gauss*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS35 + crLHS90);
rLHS(3,1)+=w_gauss*(DN(1,0)*crLHS20 + DN(1,1)*crLHS22 + crLHS68);
rLHS(3,2)+=w_gauss*(DN(0,0)*crLHS91 + crLHS44*crLHS46 - crLHS45);
rLHS(3,3)+=w_gauss*(DN(1,0)*crLHS29 + DN(1,1)*crLHS31 + crLHS7*crLHS92 + crLHS93);
rLHS(3,4)+=w_gauss*(DN(1,0)*crLHS40 + DN(1,1)*crLHS42 + crLHS95);
rLHS(3,5)+=DN(1,0)*crLHS96;
rLHS(3,6)+=w_gauss*(DN(1,0)*crLHS47 + DN(1,1)*crLHS49 + crLHS100 + crLHS99);
rLHS(3,7)+=w_gauss*(DN(1,0)*crLHS57 + DN(1,1)*crLHS59 + crLHS101);
rLHS(3,8)+=w_gauss*(DN(2,0)*crLHS91 - crLHS102 + crLHS103*crLHS46);
rLHS(4,0)+=w_gauss*(DN(1,0)*crLHS2 + DN(1,1)*crLHS63 + crLHS43);
rLHS(4,1)+=w_gauss*(DN(1,0)*crLHS22 + DN(1,1)*crLHS64 + crLHS71 + crLHS90);
rLHS(4,2)+=w_gauss*(DN(0,1)*crLHS91 + crLHS46*crLHS72 - crLHS73);
rLHS(4,3)+=w_gauss*(DN(1,0)*crLHS31 + DN(1,1)*crLHS66 + crLHS95);
rLHS(4,4)+=w_gauss*(DN(1,0)*crLHS42 + DN(1,1)*crLHS69 + crLHS104*crLHS7 + crLHS93);
rLHS(4,5)+=DN(1,1)*crLHS96;
rLHS(4,6)+=w_gauss*(DN(1,0)*crLHS49 + DN(1,1)*crLHS74 + crLHS106);
rLHS(4,7)+=w_gauss*(DN(1,0)*crLHS59 + DN(1,1)*crLHS76 + crLHS100 + crLHS108);
rLHS(4,8)+=w_gauss*(DN(2,1)*crLHS91 - crLHS109 + crLHS110*crLHS46);
rLHS(5,0)+=w_gauss*(DN(1,0)*crLHS81 + crLHS44);
rLHS(5,1)+=w_gauss*(DN(1,1)*crLHS81 + crLHS72);
rLHS(5,2)+=crLHS85;
rLHS(5,3)+=DN(1,0)*crLHS111;
rLHS(5,4)+=DN(1,1)*crLHS111;
rLHS(5,5)+=crLHS83*(crLHS104 + crLHS92);
rLHS(5,6)+=w_gauss*(DN(1,0)*crLHS86 + crLHS103);
rLHS(5,7)+=w_gauss*(DN(1,1)*crLHS86 + crLHS110);
rLHS(5,8)+=crLHS112;
rLHS(6,0)+=w_gauss*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS115 + crLHS53);
rLHS(6,1)+=w_gauss*(DN(2,0)*crLHS20 + DN(2,1)*crLHS22 + crLHS75);
rLHS(6,2)+=w_gauss*(DN(0,0)*crLHS116 + crLHS46*crLHS61 - crLHS62);
rLHS(6,3)+=w_gauss*(DN(2,0)*crLHS29 + DN(2,1)*crLHS31 + crLHS117 + crLHS99);
rLHS(6,4)+=w_gauss*(DN(2,0)*crLHS40 + DN(2,1)*crLHS42 + crLHS106);
rLHS(6,5)+=w_gauss*(DN(1,0)*crLHS116 + crLHS102*crLHS46 - crLHS103);
rLHS(6,6)+=w_gauss*(DN(2,0)*crLHS47 + DN(2,1)*crLHS49 + crLHS118*crLHS7 + crLHS119);
rLHS(6,7)+=w_gauss*(DN(2,0)*crLHS57 + DN(2,1)*crLHS59 + crLHS120);
rLHS(6,8)+=DN(2,0)*crLHS121;
rLHS(7,0)+=w_gauss*(DN(2,0)*crLHS2 + DN(2,1)*crLHS63 + crLHS60);
rLHS(7,1)+=w_gauss*(DN(2,0)*crLHS22 + DN(2,1)*crLHS64 + crLHS115 + crLHS78);
rLHS(7,2)+=w_gauss*(DN(0,1)*crLHS116 + crLHS46*crLHS79 - crLHS80);
rLHS(7,3)+=w_gauss*(DN(2,0)*crLHS31 + DN(2,1)*crLHS66 + crLHS101);
rLHS(7,4)+=w_gauss*(DN(2,0)*crLHS42 + DN(2,1)*crLHS69 + crLHS108 + crLHS117);
rLHS(7,5)+=w_gauss*(DN(1,1)*crLHS116 + crLHS109*crLHS46 - crLHS110);
rLHS(7,6)+=w_gauss*(DN(2,0)*crLHS49 + DN(2,1)*crLHS74 + crLHS120);
rLHS(7,7)+=w_gauss*(DN(2,0)*crLHS59 + DN(2,1)*crLHS76 + crLHS119 + crLHS122*crLHS7);
rLHS(7,8)+=DN(2,1)*crLHS121;
rLHS(8,0)+=w_gauss*(DN(2,0)*crLHS81 + crLHS61);
rLHS(8,1)+=w_gauss*(DN(2,1)*crLHS81 + crLHS79);
rLHS(8,2)+=crLHS87;
rLHS(8,3)+=w_gauss*(DN(2,0)*crLHS84 + crLHS102);
rLHS(8,4)+=w_gauss*(DN(2,1)*crLHS84 + crLHS109);
rLHS(8,5)+=crLHS112;
rLHS(8,6)+=DN(2,0)*crLHS123;
rLHS(8,7)+=DN(2,1)*crLHS123;
rLHS(8,8)+=crLHS83*(crLHS118 + crLHS122);

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
    
    // Add LHS Gauss point contribution
    const double w_gauss = rData.Weight;

    const double crLHS0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS1 = C(0,3)*DN(0,0);
const double crLHS2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS1;
const double crLHS3 = C(0,5)*DN(0,0);
const double crLHS4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS3;
const double crLHS5 = DN(0,0)*DN(0,0);
const double crLHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS9 = rho*stab_c2*sqrt(crLHS6*crLHS6 + crLHS7*crLHS7 + crLHS8*crLHS8);
const double crLHS10 = crLHS9*h*1.0/stab_c1 + mu;
const double crLHS11 = DN(0,0)*crLHS6 + DN(0,1)*crLHS7 + DN(0,2)*crLHS8;
const double crLHS12 = N[0]*rho;
const double crLHS13 = bdf0*rho;
const double crLHS14 = N[0]*bdf0;
const double crLHS15 = crLHS11 + crLHS14;
const double crLHS16 = 1.0*1.0/(crLHS9*1.0/h + dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h));
const double crLHS17 = crLHS16*(rho*rho);
const double crLHS18 = crLHS11*crLHS17;
const double crLHS19 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS20 = crLHS17*crLHS19;
const double crLHS21 = N[0]*crLHS20;
const double crLHS22 = crLHS11*crLHS12 + crLHS13*(N[0]*N[0]) + crLHS15*crLHS18 + crLHS15*crLHS21;
const double crLHS23 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS24 = C(1,3)*DN(0,1);
const double crLHS25 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS24;
const double crLHS26 = C(3,5)*DN(0,0);
const double crLHS27 = C(4,5)*DN(0,2);
const double crLHS28 = C(1,5)*DN(0,1) + crLHS26 + crLHS27;
const double crLHS29 = DN(0,0)*crLHS10;
const double crLHS30 = DN(0,1)*crLHS29;
const double crLHS31 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS32 = C(3,4)*DN(0,1);
const double crLHS33 = C(2,3)*DN(0,2) + crLHS26 + crLHS32;
const double crLHS34 = C(2,5)*DN(0,2);
const double crLHS35 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS34;
const double crLHS36 = DN(0,2)*crLHS29;
const double crLHS37 = crLHS16*crLHS19;
const double crLHS38 = crLHS16*rho;
const double crLHS39 = crLHS11*crLHS38;
const double crLHS40 = w_gauss*(-N[0] + crLHS12*crLHS37 + crLHS39);
const double crLHS41 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS42 = C(0,3)*DN(1,0);
const double crLHS43 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS42;
const double crLHS44 = C(0,5)*DN(1,0);
const double crLHS45 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS44;
const double crLHS46 = DN(0,0)*DN(1,0);
const double crLHS47 = N[1]*rho;
const double crLHS48 = crLHS14*crLHS47;
const double crLHS49 = crLHS10*crLHS46 + crLHS48;
const double crLHS50 = DN(1,0)*crLHS6 + DN(1,1)*crLHS7 + DN(1,2)*crLHS8;
const double crLHS51 = N[1]*bdf0;
const double crLHS52 = crLHS50 + crLHS51;
const double crLHS53 = crLHS12*crLHS50 + crLHS18*crLHS52 + crLHS21*crLHS52;
const double crLHS54 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS42;
const double crLHS55 = C(1,3)*DN(1,1);
const double crLHS56 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS55;
const double crLHS57 = C(3,5)*DN(1,0);
const double crLHS58 = C(4,5)*DN(1,2);
const double crLHS59 = C(1,5)*DN(1,1) + crLHS57 + crLHS58;
const double crLHS60 = DN(1,1)*crLHS29;
const double crLHS61 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS44;
const double crLHS62 = C(3,4)*DN(1,1);
const double crLHS63 = C(2,3)*DN(1,2) + crLHS57 + crLHS62;
const double crLHS64 = C(2,5)*DN(1,2);
const double crLHS65 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS64;
const double crLHS66 = DN(1,2)*crLHS29;
const double crLHS67 = DN(0,0)*N[1];
const double crLHS68 = DN(1,0)*N[0];
const double crLHS69 = crLHS19*crLHS38;
const double crLHS70 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS71 = C(0,3)*DN(2,0);
const double crLHS72 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS71;
const double crLHS73 = C(0,5)*DN(2,0);
const double crLHS74 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS73;
const double crLHS75 = DN(0,0)*DN(2,0);
const double crLHS76 = N[2]*rho;
const double crLHS77 = crLHS14*crLHS76;
const double crLHS78 = crLHS10*crLHS75 + crLHS77;
const double crLHS79 = DN(2,0)*crLHS6 + DN(2,1)*crLHS7 + DN(2,2)*crLHS8;
const double crLHS80 = N[2]*bdf0;
const double crLHS81 = crLHS79 + crLHS80;
const double crLHS82 = crLHS12*crLHS79 + crLHS18*crLHS81 + crLHS21*crLHS81;
const double crLHS83 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS71;
const double crLHS84 = C(1,3)*DN(2,1);
const double crLHS85 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS84;
const double crLHS86 = C(3,5)*DN(2,0);
const double crLHS87 = C(4,5)*DN(2,2);
const double crLHS88 = C(1,5)*DN(2,1) + crLHS86 + crLHS87;
const double crLHS89 = DN(2,1)*crLHS29;
const double crLHS90 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS73;
const double crLHS91 = C(3,4)*DN(2,1);
const double crLHS92 = C(2,3)*DN(2,2) + crLHS86 + crLHS91;
const double crLHS93 = C(2,5)*DN(2,2);
const double crLHS94 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS93;
const double crLHS95 = DN(2,2)*crLHS29;
const double crLHS96 = DN(0,0)*N[2];
const double crLHS97 = DN(2,0)*N[0];
const double crLHS98 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS99 = C(0,3)*DN(3,0);
const double crLHS100 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS99;
const double crLHS101 = C(0,5)*DN(3,0);
const double crLHS102 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS101;
const double crLHS103 = DN(0,0)*DN(3,0);
const double crLHS104 = N[3]*rho;
const double crLHS105 = crLHS104*crLHS14;
const double crLHS106 = crLHS10*crLHS103 + crLHS105;
const double crLHS107 = DN(3,0)*crLHS6 + DN(3,1)*crLHS7 + DN(3,2)*crLHS8;
const double crLHS108 = N[3]*bdf0 + crLHS107;
const double crLHS109 = crLHS107*crLHS12 + crLHS108*crLHS18 + crLHS108*crLHS21;
const double crLHS110 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS99;
const double crLHS111 = C(1,3)*DN(3,1);
const double crLHS112 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS111;
const double crLHS113 = C(3,5)*DN(3,0);
const double crLHS114 = C(4,5)*DN(3,2);
const double crLHS115 = C(1,5)*DN(3,1) + crLHS113 + crLHS114;
const double crLHS116 = DN(3,1)*crLHS29;
const double crLHS117 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS101;
const double crLHS118 = C(3,4)*DN(3,1);
const double crLHS119 = C(2,3)*DN(3,2) + crLHS113 + crLHS118;
const double crLHS120 = C(2,5)*DN(3,2);
const double crLHS121 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS120;
const double crLHS122 = DN(3,2)*crLHS29;
const double crLHS123 = DN(0,0)*N[3];
const double crLHS124 = DN(3,0)*N[0];
const double crLHS125 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS24;
const double crLHS126 = C(0,4)*DN(0,0) + crLHS27 + crLHS32;
const double crLHS127 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS128 = C(1,4)*DN(0,1);
const double crLHS129 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS128;
const double crLHS130 = DN(0,1)*DN(0,1);
const double crLHS131 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS128;
const double crLHS132 = C(2,4)*DN(0,2);
const double crLHS133 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS132;
const double crLHS134 = DN(0,1)*crLHS10;
const double crLHS135 = DN(0,2)*crLHS134;
const double crLHS136 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS55;
const double crLHS137 = C(0,4)*DN(1,0) + crLHS58 + crLHS62;
const double crLHS138 = DN(1,0)*crLHS134;
const double crLHS139 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS140 = C(1,4)*DN(1,1);
const double crLHS141 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS140;
const double crLHS142 = DN(0,1)*DN(1,1);
const double crLHS143 = crLHS10*crLHS142;
const double crLHS144 = crLHS48 + crLHS53;
const double crLHS145 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS140;
const double crLHS146 = C(2,4)*DN(1,2);
const double crLHS147 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS146;
const double crLHS148 = DN(1,2)*crLHS134;
const double crLHS149 = DN(0,1)*N[1];
const double crLHS150 = DN(1,1)*N[0];
const double crLHS151 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS84;
const double crLHS152 = C(0,4)*DN(2,0) + crLHS87 + crLHS91;
const double crLHS153 = DN(2,0)*crLHS134;
const double crLHS154 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS155 = C(1,4)*DN(2,1);
const double crLHS156 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS155;
const double crLHS157 = DN(0,1)*DN(2,1);
const double crLHS158 = crLHS10*crLHS157;
const double crLHS159 = crLHS77 + crLHS82;
const double crLHS160 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS155;
const double crLHS161 = C(2,4)*DN(2,2);
const double crLHS162 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS161;
const double crLHS163 = DN(2,2)*crLHS134;
const double crLHS164 = DN(0,1)*N[2];
const double crLHS165 = DN(2,1)*N[0];
const double crLHS166 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS111;
const double crLHS167 = C(0,4)*DN(3,0) + crLHS114 + crLHS118;
const double crLHS168 = DN(3,0)*crLHS134;
const double crLHS169 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS170 = C(1,4)*DN(3,1);
const double crLHS171 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS170;
const double crLHS172 = DN(0,1)*DN(3,1);
const double crLHS173 = crLHS10*crLHS172;
const double crLHS174 = crLHS105 + crLHS109;
const double crLHS175 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS170;
const double crLHS176 = C(2,4)*DN(3,2);
const double crLHS177 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS176;
const double crLHS178 = DN(3,2)*crLHS134;
const double crLHS179 = DN(0,1)*N[3];
const double crLHS180 = DN(3,1)*N[0];
const double crLHS181 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS34;
const double crLHS182 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS132;
const double crLHS183 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS184 = DN(0,2)*DN(0,2);
const double crLHS185 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS64;
const double crLHS186 = DN(0,2)*crLHS10;
const double crLHS187 = DN(1,0)*crLHS186;
const double crLHS188 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS146;
const double crLHS189 = DN(1,1)*crLHS186;
const double crLHS190 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS191 = DN(0,2)*DN(1,2);
const double crLHS192 = crLHS10*crLHS191;
const double crLHS193 = DN(0,2)*N[1];
const double crLHS194 = DN(1,2)*N[0];
const double crLHS195 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS93;
const double crLHS196 = DN(2,0)*crLHS186;
const double crLHS197 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS161;
const double crLHS198 = DN(2,1)*crLHS186;
const double crLHS199 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS200 = DN(0,2)*DN(2,2);
const double crLHS201 = crLHS10*crLHS200;
const double crLHS202 = DN(0,2)*N[2];
const double crLHS203 = DN(2,2)*N[0];
const double crLHS204 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS120;
const double crLHS205 = DN(3,0)*crLHS186;
const double crLHS206 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS176;
const double crLHS207 = DN(3,1)*crLHS186;
const double crLHS208 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS209 = DN(0,2)*DN(3,2);
const double crLHS210 = crLHS10*crLHS209;
const double crLHS211 = DN(0,2)*N[3];
const double crLHS212 = DN(3,2)*N[0];
const double crLHS213 = crLHS15*crLHS38;
const double crLHS214 = w_gauss*(N[0] + crLHS213);
const double crLHS215 = crLHS16*w_gauss;
const double crLHS216 = crLHS38*crLHS52;
const double crLHS217 = crLHS215*(crLHS142 + crLHS191 + crLHS46);
const double crLHS218 = crLHS38*crLHS81;
const double crLHS219 = crLHS215*(crLHS157 + crLHS200 + crLHS75);
const double crLHS220 = crLHS108*crLHS38;
const double crLHS221 = crLHS215*(crLHS103 + crLHS172 + crLHS209);
const double crLHS222 = crLHS17*crLHS50;
const double crLHS223 = N[1]*crLHS20;
const double crLHS224 = crLHS11*crLHS47 + crLHS15*crLHS222 + crLHS15*crLHS223;
const double crLHS225 = crLHS38*crLHS50;
const double crLHS226 = DN(1,0)*DN(1,0);
const double crLHS227 = crLHS13*(N[1]*N[1]) + crLHS222*crLHS52 + crLHS223*crLHS52 + crLHS47*crLHS50;
const double crLHS228 = DN(1,0)*crLHS10;
const double crLHS229 = DN(1,1)*crLHS228;
const double crLHS230 = DN(1,2)*crLHS228;
const double crLHS231 = w_gauss*(-N[1] + crLHS225 + crLHS37*crLHS47);
const double crLHS232 = DN(1,0)*DN(2,0);
const double crLHS233 = crLHS51*crLHS76;
const double crLHS234 = crLHS10*crLHS232 + crLHS233;
const double crLHS235 = crLHS222*crLHS81 + crLHS223*crLHS81 + crLHS47*crLHS79;
const double crLHS236 = DN(2,1)*crLHS228;
const double crLHS237 = DN(2,2)*crLHS228;
const double crLHS238 = DN(1,0)*N[2];
const double crLHS239 = DN(2,0)*N[1];
const double crLHS240 = DN(1,0)*DN(3,0);
const double crLHS241 = crLHS104*crLHS51;
const double crLHS242 = crLHS10*crLHS240 + crLHS241;
const double crLHS243 = crLHS107*crLHS47 + crLHS108*crLHS222 + crLHS108*crLHS223;
const double crLHS244 = DN(3,1)*crLHS228;
const double crLHS245 = DN(3,2)*crLHS228;
const double crLHS246 = DN(1,0)*N[3];
const double crLHS247 = DN(3,0)*N[1];
const double crLHS248 = crLHS224 + crLHS48;
const double crLHS249 = DN(1,1)*DN(1,1);
const double crLHS250 = DN(1,1)*crLHS10;
const double crLHS251 = DN(1,2)*crLHS250;
const double crLHS252 = DN(2,0)*crLHS250;
const double crLHS253 = DN(1,1)*DN(2,1);
const double crLHS254 = crLHS10*crLHS253;
const double crLHS255 = crLHS233 + crLHS235;
const double crLHS256 = DN(2,2)*crLHS250;
const double crLHS257 = DN(1,1)*N[2];
const double crLHS258 = DN(2,1)*N[1];
const double crLHS259 = DN(3,0)*crLHS250;
const double crLHS260 = DN(1,1)*DN(3,1);
const double crLHS261 = crLHS10*crLHS260;
const double crLHS262 = crLHS241 + crLHS243;
const double crLHS263 = DN(3,2)*crLHS250;
const double crLHS264 = DN(1,1)*N[3];
const double crLHS265 = DN(3,1)*N[1];
const double crLHS266 = DN(1,2)*DN(1,2);
const double crLHS267 = DN(1,2)*crLHS10;
const double crLHS268 = DN(2,0)*crLHS267;
const double crLHS269 = DN(2,1)*crLHS267;
const double crLHS270 = DN(1,2)*DN(2,2);
const double crLHS271 = crLHS10*crLHS270;
const double crLHS272 = DN(1,2)*N[2];
const double crLHS273 = DN(2,2)*N[1];
const double crLHS274 = DN(3,0)*crLHS267;
const double crLHS275 = DN(3,1)*crLHS267;
const double crLHS276 = DN(1,2)*DN(3,2);
const double crLHS277 = crLHS10*crLHS276;
const double crLHS278 = DN(1,2)*N[3];
const double crLHS279 = DN(3,2)*N[1];
const double crLHS280 = w_gauss*(N[1] + crLHS216);
const double crLHS281 = crLHS215*(crLHS232 + crLHS253 + crLHS270);
const double crLHS282 = crLHS215*(crLHS240 + crLHS260 + crLHS276);
const double crLHS283 = crLHS17*crLHS79;
const double crLHS284 = N[2]*crLHS20;
const double crLHS285 = crLHS11*crLHS76 + crLHS15*crLHS283 + crLHS15*crLHS284;
const double crLHS286 = crLHS38*crLHS79;
const double crLHS287 = crLHS283*crLHS52 + crLHS284*crLHS52 + crLHS50*crLHS76;
const double crLHS288 = DN(2,0)*DN(2,0);
const double crLHS289 = crLHS13*(N[2]*N[2]) + crLHS283*crLHS81 + crLHS284*crLHS81 + crLHS76*crLHS79;
const double crLHS290 = DN(2,0)*crLHS10;
const double crLHS291 = DN(2,1)*crLHS290;
const double crLHS292 = DN(2,2)*crLHS290;
const double crLHS293 = w_gauss*(-N[2] + crLHS286 + crLHS37*crLHS76);
const double crLHS294 = DN(2,0)*DN(3,0);
const double crLHS295 = crLHS104*crLHS80;
const double crLHS296 = crLHS10*crLHS294 + crLHS295;
const double crLHS297 = crLHS107*crLHS76 + crLHS108*crLHS283 + crLHS108*crLHS284;
const double crLHS298 = DN(3,1)*crLHS290;
const double crLHS299 = DN(3,2)*crLHS290;
const double crLHS300 = DN(2,0)*N[3];
const double crLHS301 = DN(3,0)*N[2];
const double crLHS302 = crLHS285 + crLHS77;
const double crLHS303 = crLHS233 + crLHS287;
const double crLHS304 = DN(2,1)*DN(2,1);
const double crLHS305 = DN(2,1)*crLHS10;
const double crLHS306 = DN(2,2)*crLHS305;
const double crLHS307 = DN(3,0)*crLHS305;
const double crLHS308 = DN(2,1)*DN(3,1);
const double crLHS309 = crLHS10*crLHS308;
const double crLHS310 = crLHS295 + crLHS297;
const double crLHS311 = DN(3,2)*crLHS305;
const double crLHS312 = DN(2,1)*N[3];
const double crLHS313 = DN(3,1)*N[2];
const double crLHS314 = DN(2,2)*DN(2,2);
const double crLHS315 = DN(2,2)*crLHS10;
const double crLHS316 = DN(3,0)*crLHS315;
const double crLHS317 = DN(3,1)*crLHS315;
const double crLHS318 = DN(2,2)*DN(3,2);
const double crLHS319 = crLHS10*crLHS318;
const double crLHS320 = DN(2,2)*N[3];
const double crLHS321 = DN(3,2)*N[2];
const double crLHS322 = w_gauss*(N[2] + crLHS218);
const double crLHS323 = crLHS215*(crLHS294 + crLHS308 + crLHS318);
const double crLHS324 = crLHS107*crLHS17;
const double crLHS325 = N[3]*crLHS20;
const double crLHS326 = crLHS104*crLHS11 + crLHS15*crLHS324 + crLHS15*crLHS325;
const double crLHS327 = crLHS107*crLHS38;
const double crLHS328 = crLHS104*crLHS50 + crLHS324*crLHS52 + crLHS325*crLHS52;
const double crLHS329 = crLHS104*crLHS79 + crLHS324*crLHS81 + crLHS325*crLHS81;
const double crLHS330 = DN(3,0)*DN(3,0);
const double crLHS331 = crLHS104*crLHS107 + crLHS108*crLHS324 + crLHS108*crLHS325 + crLHS13*(N[3]*N[3]);
const double crLHS332 = DN(3,0)*crLHS10;
const double crLHS333 = DN(3,1)*crLHS332;
const double crLHS334 = DN(3,2)*crLHS332;
const double crLHS335 = w_gauss*(-N[3] + crLHS104*crLHS37 + crLHS327);
const double crLHS336 = crLHS105 + crLHS326;
const double crLHS337 = crLHS241 + crLHS328;
const double crLHS338 = crLHS295 + crLHS329;
const double crLHS339 = DN(3,1)*DN(3,1);
const double crLHS340 = DN(3,1)*DN(3,2)*crLHS10;
const double crLHS341 = DN(3,2)*DN(3,2);
const double crLHS342 = w_gauss*(N[3] + crLHS220);
rLHS(0,0)+=w_gauss*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS10*crLHS5 + crLHS22);
rLHS(0,1)+=w_gauss*(DN(0,0)*crLHS23 + DN(0,1)*crLHS25 + DN(0,2)*crLHS28 + crLHS30);
rLHS(0,2)+=w_gauss*(DN(0,0)*crLHS31 + DN(0,1)*crLHS33 + DN(0,2)*crLHS35 + crLHS36);
rLHS(0,3)+=DN(0,0)*crLHS40;
rLHS(0,4)+=w_gauss*(DN(0,0)*crLHS41 + DN(0,1)*crLHS43 + DN(0,2)*crLHS45 + crLHS49 + crLHS53);
rLHS(0,5)+=w_gauss*(DN(0,0)*crLHS54 + DN(0,1)*crLHS56 + DN(0,2)*crLHS59 + crLHS60);
rLHS(0,6)+=w_gauss*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + DN(0,2)*crLHS65 + crLHS66);
rLHS(0,7)+=w_gauss*(DN(1,0)*crLHS39 - crLHS67 + crLHS68*crLHS69);
rLHS(0,8)+=w_gauss*(DN(0,0)*crLHS70 + DN(0,1)*crLHS72 + DN(0,2)*crLHS74 + crLHS78 + crLHS82);
rLHS(0,9)+=w_gauss*(DN(0,0)*crLHS83 + DN(0,1)*crLHS85 + DN(0,2)*crLHS88 + crLHS89);
rLHS(0,10)+=w_gauss*(DN(0,0)*crLHS90 + DN(0,1)*crLHS92 + DN(0,2)*crLHS94 + crLHS95);
rLHS(0,11)+=w_gauss*(DN(2,0)*crLHS39 + crLHS69*crLHS97 - crLHS96);
rLHS(0,12)+=w_gauss*(DN(0,0)*crLHS98 + DN(0,1)*crLHS100 + DN(0,2)*crLHS102 + crLHS106 + crLHS109);
rLHS(0,13)+=w_gauss*(DN(0,0)*crLHS110 + DN(0,1)*crLHS112 + DN(0,2)*crLHS115 + crLHS116);
rLHS(0,14)+=w_gauss*(DN(0,0)*crLHS117 + DN(0,1)*crLHS119 + DN(0,2)*crLHS121 + crLHS122);
rLHS(0,15)+=w_gauss*(DN(3,0)*crLHS39 - crLHS123 + crLHS124*crLHS69);
rLHS(1,0)+=w_gauss*(DN(0,0)*crLHS2 + DN(0,1)*crLHS125 + DN(0,2)*crLHS126 + crLHS30);
rLHS(1,1)+=w_gauss*(DN(0,0)*crLHS25 + DN(0,1)*crLHS127 + DN(0,2)*crLHS129 + crLHS10*crLHS130 + crLHS22);
rLHS(1,2)+=w_gauss*(DN(0,0)*crLHS33 + DN(0,1)*crLHS131 + DN(0,2)*crLHS133 + crLHS135);
rLHS(1,3)+=DN(0,1)*crLHS40;
rLHS(1,4)+=w_gauss*(DN(0,0)*crLHS43 + DN(0,1)*crLHS136 + DN(0,2)*crLHS137 + crLHS138);
rLHS(1,5)+=w_gauss*(DN(0,0)*crLHS56 + DN(0,1)*crLHS139 + DN(0,2)*crLHS141 + crLHS143 + crLHS144);
rLHS(1,6)+=w_gauss*(DN(0,0)*crLHS63 + DN(0,1)*crLHS145 + DN(0,2)*crLHS147 + crLHS148);
rLHS(1,7)+=w_gauss*(DN(1,1)*crLHS39 - crLHS149 + crLHS150*crLHS69);
rLHS(1,8)+=w_gauss*(DN(0,0)*crLHS72 + DN(0,1)*crLHS151 + DN(0,2)*crLHS152 + crLHS153);
rLHS(1,9)+=w_gauss*(DN(0,0)*crLHS85 + DN(0,1)*crLHS154 + DN(0,2)*crLHS156 + crLHS158 + crLHS159);
rLHS(1,10)+=w_gauss*(DN(0,0)*crLHS92 + DN(0,1)*crLHS160 + DN(0,2)*crLHS162 + crLHS163);
rLHS(1,11)+=w_gauss*(DN(2,1)*crLHS39 - crLHS164 + crLHS165*crLHS69);
rLHS(1,12)+=w_gauss*(DN(0,0)*crLHS100 + DN(0,1)*crLHS166 + DN(0,2)*crLHS167 + crLHS168);
rLHS(1,13)+=w_gauss*(DN(0,0)*crLHS112 + DN(0,1)*crLHS169 + DN(0,2)*crLHS171 + crLHS173 + crLHS174);
rLHS(1,14)+=w_gauss*(DN(0,0)*crLHS119 + DN(0,1)*crLHS175 + DN(0,2)*crLHS177 + crLHS178);
rLHS(1,15)+=w_gauss*(DN(3,1)*crLHS39 - crLHS179 + crLHS180*crLHS69);
rLHS(2,0)+=w_gauss*(DN(0,0)*crLHS4 + DN(0,1)*crLHS126 + DN(0,2)*crLHS181 + crLHS36);
rLHS(2,1)+=w_gauss*(DN(0,0)*crLHS28 + DN(0,1)*crLHS129 + DN(0,2)*crLHS182 + crLHS135);
rLHS(2,2)+=w_gauss*(DN(0,0)*crLHS35 + DN(0,1)*crLHS133 + DN(0,2)*crLHS183 + crLHS10*crLHS184 + crLHS22);
rLHS(2,3)+=DN(0,2)*crLHS40;
rLHS(2,4)+=w_gauss*(DN(0,0)*crLHS45 + DN(0,1)*crLHS137 + DN(0,2)*crLHS185 + crLHS187);
rLHS(2,5)+=w_gauss*(DN(0,0)*crLHS59 + DN(0,1)*crLHS141 + DN(0,2)*crLHS188 + crLHS189);
rLHS(2,6)+=w_gauss*(DN(0,0)*crLHS65 + DN(0,1)*crLHS147 + DN(0,2)*crLHS190 + crLHS144 + crLHS192);
rLHS(2,7)+=w_gauss*(DN(1,2)*crLHS39 - crLHS193 + crLHS194*crLHS69);
rLHS(2,8)+=w_gauss*(DN(0,0)*crLHS74 + DN(0,1)*crLHS152 + DN(0,2)*crLHS195 + crLHS196);
rLHS(2,9)+=w_gauss*(DN(0,0)*crLHS88 + DN(0,1)*crLHS156 + DN(0,2)*crLHS197 + crLHS198);
rLHS(2,10)+=w_gauss*(DN(0,0)*crLHS94 + DN(0,1)*crLHS162 + DN(0,2)*crLHS199 + crLHS159 + crLHS201);
rLHS(2,11)+=w_gauss*(DN(2,2)*crLHS39 - crLHS202 + crLHS203*crLHS69);
rLHS(2,12)+=w_gauss*(DN(0,0)*crLHS102 + DN(0,1)*crLHS167 + DN(0,2)*crLHS204 + crLHS205);
rLHS(2,13)+=w_gauss*(DN(0,0)*crLHS115 + DN(0,1)*crLHS171 + DN(0,2)*crLHS206 + crLHS207);
rLHS(2,14)+=w_gauss*(DN(0,0)*crLHS121 + DN(0,1)*crLHS177 + DN(0,2)*crLHS208 + crLHS174 + crLHS210);
rLHS(2,15)+=w_gauss*(DN(3,2)*crLHS39 - crLHS211 + crLHS212*crLHS69);
rLHS(3,0)+=DN(0,0)*crLHS214;
rLHS(3,1)+=DN(0,1)*crLHS214;
rLHS(3,2)+=DN(0,2)*crLHS214;
rLHS(3,3)+=crLHS215*(crLHS130 + crLHS184 + crLHS5);
rLHS(3,4)+=w_gauss*(DN(0,0)*crLHS216 + crLHS68);
rLHS(3,5)+=w_gauss*(DN(0,1)*crLHS216 + crLHS150);
rLHS(3,6)+=w_gauss*(DN(0,2)*crLHS216 + crLHS194);
rLHS(3,7)+=crLHS217;
rLHS(3,8)+=w_gauss*(DN(0,0)*crLHS218 + crLHS97);
rLHS(3,9)+=w_gauss*(DN(0,1)*crLHS218 + crLHS165);
rLHS(3,10)+=w_gauss*(DN(0,2)*crLHS218 + crLHS203);
rLHS(3,11)+=crLHS219;
rLHS(3,12)+=w_gauss*(DN(0,0)*crLHS220 + crLHS124);
rLHS(3,13)+=w_gauss*(DN(0,1)*crLHS220 + crLHS180);
rLHS(3,14)+=w_gauss*(DN(0,2)*crLHS220 + crLHS212);
rLHS(3,15)+=crLHS221;
rLHS(4,0)+=w_gauss*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS224 + crLHS49);
rLHS(4,1)+=w_gauss*(DN(1,0)*crLHS23 + DN(1,1)*crLHS25 + DN(1,2)*crLHS28 + crLHS138);
rLHS(4,2)+=w_gauss*(DN(1,0)*crLHS31 + DN(1,1)*crLHS33 + DN(1,2)*crLHS35 + crLHS187);
rLHS(4,3)+=w_gauss*(DN(0,0)*crLHS225 + crLHS67*crLHS69 - crLHS68);
rLHS(4,4)+=w_gauss*(DN(1,0)*crLHS41 + DN(1,1)*crLHS43 + DN(1,2)*crLHS45 + crLHS10*crLHS226 + crLHS227);
rLHS(4,5)+=w_gauss*(DN(1,0)*crLHS54 + DN(1,1)*crLHS56 + DN(1,2)*crLHS59 + crLHS229);
rLHS(4,6)+=w_gauss*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + DN(1,2)*crLHS65 + crLHS230);
rLHS(4,7)+=DN(1,0)*crLHS231;
rLHS(4,8)+=w_gauss*(DN(1,0)*crLHS70 + DN(1,1)*crLHS72 + DN(1,2)*crLHS74 + crLHS234 + crLHS235);
rLHS(4,9)+=w_gauss*(DN(1,0)*crLHS83 + DN(1,1)*crLHS85 + DN(1,2)*crLHS88 + crLHS236);
rLHS(4,10)+=w_gauss*(DN(1,0)*crLHS90 + DN(1,1)*crLHS92 + DN(1,2)*crLHS94 + crLHS237);
rLHS(4,11)+=w_gauss*(DN(2,0)*crLHS225 - crLHS238 + crLHS239*crLHS69);
rLHS(4,12)+=w_gauss*(DN(1,0)*crLHS98 + DN(1,1)*crLHS100 + DN(1,2)*crLHS102 + crLHS242 + crLHS243);
rLHS(4,13)+=w_gauss*(DN(1,0)*crLHS110 + DN(1,1)*crLHS112 + DN(1,2)*crLHS115 + crLHS244);
rLHS(4,14)+=w_gauss*(DN(1,0)*crLHS117 + DN(1,1)*crLHS119 + DN(1,2)*crLHS121 + crLHS245);
rLHS(4,15)+=w_gauss*(DN(3,0)*crLHS225 - crLHS246 + crLHS247*crLHS69);
rLHS(5,0)+=w_gauss*(DN(1,0)*crLHS2 + DN(1,1)*crLHS125 + DN(1,2)*crLHS126 + crLHS60);
rLHS(5,1)+=w_gauss*(DN(1,0)*crLHS25 + DN(1,1)*crLHS127 + DN(1,2)*crLHS129 + crLHS143 + crLHS248);
rLHS(5,2)+=w_gauss*(DN(1,0)*crLHS33 + DN(1,1)*crLHS131 + DN(1,2)*crLHS133 + crLHS189);
rLHS(5,3)+=w_gauss*(DN(0,1)*crLHS225 + crLHS149*crLHS69 - crLHS150);
rLHS(5,4)+=w_gauss*(DN(1,0)*crLHS43 + DN(1,1)*crLHS136 + DN(1,2)*crLHS137 + crLHS229);
rLHS(5,5)+=w_gauss*(DN(1,0)*crLHS56 + DN(1,1)*crLHS139 + DN(1,2)*crLHS141 + crLHS10*crLHS249 + crLHS227);
rLHS(5,6)+=w_gauss*(DN(1,0)*crLHS63 + DN(1,1)*crLHS145 + DN(1,2)*crLHS147 + crLHS251);
rLHS(5,7)+=DN(1,1)*crLHS231;
rLHS(5,8)+=w_gauss*(DN(1,0)*crLHS72 + DN(1,1)*crLHS151 + DN(1,2)*crLHS152 + crLHS252);
rLHS(5,9)+=w_gauss*(DN(1,0)*crLHS85 + DN(1,1)*crLHS154 + DN(1,2)*crLHS156 + crLHS254 + crLHS255);
rLHS(5,10)+=w_gauss*(DN(1,0)*crLHS92 + DN(1,1)*crLHS160 + DN(1,2)*crLHS162 + crLHS256);
rLHS(5,11)+=w_gauss*(DN(2,1)*crLHS225 - crLHS257 + crLHS258*crLHS69);
rLHS(5,12)+=w_gauss*(DN(1,0)*crLHS100 + DN(1,1)*crLHS166 + DN(1,2)*crLHS167 + crLHS259);
rLHS(5,13)+=w_gauss*(DN(1,0)*crLHS112 + DN(1,1)*crLHS169 + DN(1,2)*crLHS171 + crLHS261 + crLHS262);
rLHS(5,14)+=w_gauss*(DN(1,0)*crLHS119 + DN(1,1)*crLHS175 + DN(1,2)*crLHS177 + crLHS263);
rLHS(5,15)+=w_gauss*(DN(3,1)*crLHS225 - crLHS264 + crLHS265*crLHS69);
rLHS(6,0)+=w_gauss*(DN(1,0)*crLHS4 + DN(1,1)*crLHS126 + DN(1,2)*crLHS181 + crLHS66);
rLHS(6,1)+=w_gauss*(DN(1,0)*crLHS28 + DN(1,1)*crLHS129 + DN(1,2)*crLHS182 + crLHS148);
rLHS(6,2)+=w_gauss*(DN(1,0)*crLHS35 + DN(1,1)*crLHS133 + DN(1,2)*crLHS183 + crLHS192 + crLHS248);
rLHS(6,3)+=w_gauss*(DN(0,2)*crLHS225 + crLHS193*crLHS69 - crLHS194);
rLHS(6,4)+=w_gauss*(DN(1,0)*crLHS45 + DN(1,1)*crLHS137 + DN(1,2)*crLHS185 + crLHS230);
rLHS(6,5)+=w_gauss*(DN(1,0)*crLHS59 + DN(1,1)*crLHS141 + DN(1,2)*crLHS188 + crLHS251);
rLHS(6,6)+=w_gauss*(DN(1,0)*crLHS65 + DN(1,1)*crLHS147 + DN(1,2)*crLHS190 + crLHS10*crLHS266 + crLHS227);
rLHS(6,7)+=DN(1,2)*crLHS231;
rLHS(6,8)+=w_gauss*(DN(1,0)*crLHS74 + DN(1,1)*crLHS152 + DN(1,2)*crLHS195 + crLHS268);
rLHS(6,9)+=w_gauss*(DN(1,0)*crLHS88 + DN(1,1)*crLHS156 + DN(1,2)*crLHS197 + crLHS269);
rLHS(6,10)+=w_gauss*(DN(1,0)*crLHS94 + DN(1,1)*crLHS162 + DN(1,2)*crLHS199 + crLHS255 + crLHS271);
rLHS(6,11)+=w_gauss*(DN(2,2)*crLHS225 - crLHS272 + crLHS273*crLHS69);
rLHS(6,12)+=w_gauss*(DN(1,0)*crLHS102 + DN(1,1)*crLHS167 + DN(1,2)*crLHS204 + crLHS274);
rLHS(6,13)+=w_gauss*(DN(1,0)*crLHS115 + DN(1,1)*crLHS171 + DN(1,2)*crLHS206 + crLHS275);
rLHS(6,14)+=w_gauss*(DN(1,0)*crLHS121 + DN(1,1)*crLHS177 + DN(1,2)*crLHS208 + crLHS262 + crLHS277);
rLHS(6,15)+=w_gauss*(DN(3,2)*crLHS225 - crLHS278 + crLHS279*crLHS69);
rLHS(7,0)+=w_gauss*(DN(1,0)*crLHS213 + crLHS67);
rLHS(7,1)+=w_gauss*(DN(1,1)*crLHS213 + crLHS149);
rLHS(7,2)+=w_gauss*(DN(1,2)*crLHS213 + crLHS193);
rLHS(7,3)+=crLHS217;
rLHS(7,4)+=DN(1,0)*crLHS280;
rLHS(7,5)+=DN(1,1)*crLHS280;
rLHS(7,6)+=DN(1,2)*crLHS280;
rLHS(7,7)+=crLHS215*(crLHS226 + crLHS249 + crLHS266);
rLHS(7,8)+=w_gauss*(DN(1,0)*crLHS218 + crLHS239);
rLHS(7,9)+=w_gauss*(DN(1,1)*crLHS218 + crLHS258);
rLHS(7,10)+=w_gauss*(DN(1,2)*crLHS218 + crLHS273);
rLHS(7,11)+=crLHS281;
rLHS(7,12)+=w_gauss*(DN(1,0)*crLHS220 + crLHS247);
rLHS(7,13)+=w_gauss*(DN(1,1)*crLHS220 + crLHS265);
rLHS(7,14)+=w_gauss*(DN(1,2)*crLHS220 + crLHS279);
rLHS(7,15)+=crLHS282;
rLHS(8,0)+=w_gauss*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS285 + crLHS78);
rLHS(8,1)+=w_gauss*(DN(2,0)*crLHS23 + DN(2,1)*crLHS25 + DN(2,2)*crLHS28 + crLHS153);
rLHS(8,2)+=w_gauss*(DN(2,0)*crLHS31 + DN(2,1)*crLHS33 + DN(2,2)*crLHS35 + crLHS196);
rLHS(8,3)+=w_gauss*(DN(0,0)*crLHS286 + crLHS69*crLHS96 - crLHS97);
rLHS(8,4)+=w_gauss*(DN(2,0)*crLHS41 + DN(2,1)*crLHS43 + DN(2,2)*crLHS45 + crLHS234 + crLHS287);
rLHS(8,5)+=w_gauss*(DN(2,0)*crLHS54 + DN(2,1)*crLHS56 + DN(2,2)*crLHS59 + crLHS252);
rLHS(8,6)+=w_gauss*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + DN(2,2)*crLHS65 + crLHS268);
rLHS(8,7)+=w_gauss*(DN(1,0)*crLHS286 + crLHS238*crLHS69 - crLHS239);
rLHS(8,8)+=w_gauss*(DN(2,0)*crLHS70 + DN(2,1)*crLHS72 + DN(2,2)*crLHS74 + crLHS10*crLHS288 + crLHS289);
rLHS(8,9)+=w_gauss*(DN(2,0)*crLHS83 + DN(2,1)*crLHS85 + DN(2,2)*crLHS88 + crLHS291);
rLHS(8,10)+=w_gauss*(DN(2,0)*crLHS90 + DN(2,1)*crLHS92 + DN(2,2)*crLHS94 + crLHS292);
rLHS(8,11)+=DN(2,0)*crLHS293;
rLHS(8,12)+=w_gauss*(DN(2,0)*crLHS98 + DN(2,1)*crLHS100 + DN(2,2)*crLHS102 + crLHS296 + crLHS297);
rLHS(8,13)+=w_gauss*(DN(2,0)*crLHS110 + DN(2,1)*crLHS112 + DN(2,2)*crLHS115 + crLHS298);
rLHS(8,14)+=w_gauss*(DN(2,0)*crLHS117 + DN(2,1)*crLHS119 + DN(2,2)*crLHS121 + crLHS299);
rLHS(8,15)+=w_gauss*(DN(3,0)*crLHS286 - crLHS300 + crLHS301*crLHS69);
rLHS(9,0)+=w_gauss*(DN(2,0)*crLHS2 + DN(2,1)*crLHS125 + DN(2,2)*crLHS126 + crLHS89);
rLHS(9,1)+=w_gauss*(DN(2,0)*crLHS25 + DN(2,1)*crLHS127 + DN(2,2)*crLHS129 + crLHS158 + crLHS302);
rLHS(9,2)+=w_gauss*(DN(2,0)*crLHS33 + DN(2,1)*crLHS131 + DN(2,2)*crLHS133 + crLHS198);
rLHS(9,3)+=w_gauss*(DN(0,1)*crLHS286 + crLHS164*crLHS69 - crLHS165);
rLHS(9,4)+=w_gauss*(DN(2,0)*crLHS43 + DN(2,1)*crLHS136 + DN(2,2)*crLHS137 + crLHS236);
rLHS(9,5)+=w_gauss*(DN(2,0)*crLHS56 + DN(2,1)*crLHS139 + DN(2,2)*crLHS141 + crLHS254 + crLHS303);
rLHS(9,6)+=w_gauss*(DN(2,0)*crLHS63 + DN(2,1)*crLHS145 + DN(2,2)*crLHS147 + crLHS269);
rLHS(9,7)+=w_gauss*(DN(1,1)*crLHS286 + crLHS257*crLHS69 - crLHS258);
rLHS(9,8)+=w_gauss*(DN(2,0)*crLHS72 + DN(2,1)*crLHS151 + DN(2,2)*crLHS152 + crLHS291);
rLHS(9,9)+=w_gauss*(DN(2,0)*crLHS85 + DN(2,1)*crLHS154 + DN(2,2)*crLHS156 + crLHS10*crLHS304 + crLHS289);
rLHS(9,10)+=w_gauss*(DN(2,0)*crLHS92 + DN(2,1)*crLHS160 + DN(2,2)*crLHS162 + crLHS306);
rLHS(9,11)+=DN(2,1)*crLHS293;
rLHS(9,12)+=w_gauss*(DN(2,0)*crLHS100 + DN(2,1)*crLHS166 + DN(2,2)*crLHS167 + crLHS307);
rLHS(9,13)+=w_gauss*(DN(2,0)*crLHS112 + DN(2,1)*crLHS169 + DN(2,2)*crLHS171 + crLHS309 + crLHS310);
rLHS(9,14)+=w_gauss*(DN(2,0)*crLHS119 + DN(2,1)*crLHS175 + DN(2,2)*crLHS177 + crLHS311);
rLHS(9,15)+=w_gauss*(DN(3,1)*crLHS286 - crLHS312 + crLHS313*crLHS69);
rLHS(10,0)+=w_gauss*(DN(2,0)*crLHS4 + DN(2,1)*crLHS126 + DN(2,2)*crLHS181 + crLHS95);
rLHS(10,1)+=w_gauss*(DN(2,0)*crLHS28 + DN(2,1)*crLHS129 + DN(2,2)*crLHS182 + crLHS163);
rLHS(10,2)+=w_gauss*(DN(2,0)*crLHS35 + DN(2,1)*crLHS133 + DN(2,2)*crLHS183 + crLHS201 + crLHS302);
rLHS(10,3)+=w_gauss*(DN(0,2)*crLHS286 + crLHS202*crLHS69 - crLHS203);
rLHS(10,4)+=w_gauss*(DN(2,0)*crLHS45 + DN(2,1)*crLHS137 + DN(2,2)*crLHS185 + crLHS237);
rLHS(10,5)+=w_gauss*(DN(2,0)*crLHS59 + DN(2,1)*crLHS141 + DN(2,2)*crLHS188 + crLHS256);
rLHS(10,6)+=w_gauss*(DN(2,0)*crLHS65 + DN(2,1)*crLHS147 + DN(2,2)*crLHS190 + crLHS271 + crLHS303);
rLHS(10,7)+=w_gauss*(DN(1,2)*crLHS286 + crLHS272*crLHS69 - crLHS273);
rLHS(10,8)+=w_gauss*(DN(2,0)*crLHS74 + DN(2,1)*crLHS152 + DN(2,2)*crLHS195 + crLHS292);
rLHS(10,9)+=w_gauss*(DN(2,0)*crLHS88 + DN(2,1)*crLHS156 + DN(2,2)*crLHS197 + crLHS306);
rLHS(10,10)+=w_gauss*(DN(2,0)*crLHS94 + DN(2,1)*crLHS162 + DN(2,2)*crLHS199 + crLHS10*crLHS314 + crLHS289);
rLHS(10,11)+=DN(2,2)*crLHS293;
rLHS(10,12)+=w_gauss*(DN(2,0)*crLHS102 + DN(2,1)*crLHS167 + DN(2,2)*crLHS204 + crLHS316);
rLHS(10,13)+=w_gauss*(DN(2,0)*crLHS115 + DN(2,1)*crLHS171 + DN(2,2)*crLHS206 + crLHS317);
rLHS(10,14)+=w_gauss*(DN(2,0)*crLHS121 + DN(2,1)*crLHS177 + DN(2,2)*crLHS208 + crLHS310 + crLHS319);
rLHS(10,15)+=w_gauss*(DN(3,2)*crLHS286 - crLHS320 + crLHS321*crLHS69);
rLHS(11,0)+=w_gauss*(DN(2,0)*crLHS213 + crLHS96);
rLHS(11,1)+=w_gauss*(DN(2,1)*crLHS213 + crLHS164);
rLHS(11,2)+=w_gauss*(DN(2,2)*crLHS213 + crLHS202);
rLHS(11,3)+=crLHS219;
rLHS(11,4)+=w_gauss*(DN(2,0)*crLHS216 + crLHS238);
rLHS(11,5)+=w_gauss*(DN(2,1)*crLHS216 + crLHS257);
rLHS(11,6)+=w_gauss*(DN(2,2)*crLHS216 + crLHS272);
rLHS(11,7)+=crLHS281;
rLHS(11,8)+=DN(2,0)*crLHS322;
rLHS(11,9)+=DN(2,1)*crLHS322;
rLHS(11,10)+=DN(2,2)*crLHS322;
rLHS(11,11)+=crLHS215*(crLHS288 + crLHS304 + crLHS314);
rLHS(11,12)+=w_gauss*(DN(2,0)*crLHS220 + crLHS301);
rLHS(11,13)+=w_gauss*(DN(2,1)*crLHS220 + crLHS313);
rLHS(11,14)+=w_gauss*(DN(2,2)*crLHS220 + crLHS321);
rLHS(11,15)+=crLHS323;
rLHS(12,0)+=w_gauss*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS106 + crLHS326);
rLHS(12,1)+=w_gauss*(DN(3,0)*crLHS23 + DN(3,1)*crLHS25 + DN(3,2)*crLHS28 + crLHS168);
rLHS(12,2)+=w_gauss*(DN(3,0)*crLHS31 + DN(3,1)*crLHS33 + DN(3,2)*crLHS35 + crLHS205);
rLHS(12,3)+=w_gauss*(DN(0,0)*crLHS327 + crLHS123*crLHS69 - crLHS124);
rLHS(12,4)+=w_gauss*(DN(3,0)*crLHS41 + DN(3,1)*crLHS43 + DN(3,2)*crLHS45 + crLHS242 + crLHS328);
rLHS(12,5)+=w_gauss*(DN(3,0)*crLHS54 + DN(3,1)*crLHS56 + DN(3,2)*crLHS59 + crLHS259);
rLHS(12,6)+=w_gauss*(DN(3,0)*crLHS61 + DN(3,1)*crLHS63 + DN(3,2)*crLHS65 + crLHS274);
rLHS(12,7)+=w_gauss*(DN(1,0)*crLHS327 + crLHS246*crLHS69 - crLHS247);
rLHS(12,8)+=w_gauss*(DN(3,0)*crLHS70 + DN(3,1)*crLHS72 + DN(3,2)*crLHS74 + crLHS296 + crLHS329);
rLHS(12,9)+=w_gauss*(DN(3,0)*crLHS83 + DN(3,1)*crLHS85 + DN(3,2)*crLHS88 + crLHS307);
rLHS(12,10)+=w_gauss*(DN(3,0)*crLHS90 + DN(3,1)*crLHS92 + DN(3,2)*crLHS94 + crLHS316);
rLHS(12,11)+=w_gauss*(DN(2,0)*crLHS327 + crLHS300*crLHS69 - crLHS301);
rLHS(12,12)+=w_gauss*(DN(3,0)*crLHS98 + DN(3,1)*crLHS100 + DN(3,2)*crLHS102 + crLHS10*crLHS330 + crLHS331);
rLHS(12,13)+=w_gauss*(DN(3,0)*crLHS110 + DN(3,1)*crLHS112 + DN(3,2)*crLHS115 + crLHS333);
rLHS(12,14)+=w_gauss*(DN(3,0)*crLHS117 + DN(3,1)*crLHS119 + DN(3,2)*crLHS121 + crLHS334);
rLHS(12,15)+=DN(3,0)*crLHS335;
rLHS(13,0)+=w_gauss*(DN(3,0)*crLHS2 + DN(3,1)*crLHS125 + DN(3,2)*crLHS126 + crLHS116);
rLHS(13,1)+=w_gauss*(DN(3,0)*crLHS25 + DN(3,1)*crLHS127 + DN(3,2)*crLHS129 + crLHS173 + crLHS336);
rLHS(13,2)+=w_gauss*(DN(3,0)*crLHS33 + DN(3,1)*crLHS131 + DN(3,2)*crLHS133 + crLHS207);
rLHS(13,3)+=w_gauss*(DN(0,1)*crLHS327 + crLHS179*crLHS69 - crLHS180);
rLHS(13,4)+=w_gauss*(DN(3,0)*crLHS43 + DN(3,1)*crLHS136 + DN(3,2)*crLHS137 + crLHS244);
rLHS(13,5)+=w_gauss*(DN(3,0)*crLHS56 + DN(3,1)*crLHS139 + DN(3,2)*crLHS141 + crLHS261 + crLHS337);
rLHS(13,6)+=w_gauss*(DN(3,0)*crLHS63 + DN(3,1)*crLHS145 + DN(3,2)*crLHS147 + crLHS275);
rLHS(13,7)+=w_gauss*(DN(1,1)*crLHS327 + crLHS264*crLHS69 - crLHS265);
rLHS(13,8)+=w_gauss*(DN(3,0)*crLHS72 + DN(3,1)*crLHS151 + DN(3,2)*crLHS152 + crLHS298);
rLHS(13,9)+=w_gauss*(DN(3,0)*crLHS85 + DN(3,1)*crLHS154 + DN(3,2)*crLHS156 + crLHS309 + crLHS338);
rLHS(13,10)+=w_gauss*(DN(3,0)*crLHS92 + DN(3,1)*crLHS160 + DN(3,2)*crLHS162 + crLHS317);
rLHS(13,11)+=w_gauss*(DN(2,1)*crLHS327 + crLHS312*crLHS69 - crLHS313);
rLHS(13,12)+=w_gauss*(DN(3,0)*crLHS100 + DN(3,1)*crLHS166 + DN(3,2)*crLHS167 + crLHS333);
rLHS(13,13)+=w_gauss*(DN(3,0)*crLHS112 + DN(3,1)*crLHS169 + DN(3,2)*crLHS171 + crLHS10*crLHS339 + crLHS331);
rLHS(13,14)+=w_gauss*(DN(3,0)*crLHS119 + DN(3,1)*crLHS175 + DN(3,2)*crLHS177 + crLHS340);
rLHS(13,15)+=DN(3,1)*crLHS335;
rLHS(14,0)+=w_gauss*(DN(3,0)*crLHS4 + DN(3,1)*crLHS126 + DN(3,2)*crLHS181 + crLHS122);
rLHS(14,1)+=w_gauss*(DN(3,0)*crLHS28 + DN(3,1)*crLHS129 + DN(3,2)*crLHS182 + crLHS178);
rLHS(14,2)+=w_gauss*(DN(3,0)*crLHS35 + DN(3,1)*crLHS133 + DN(3,2)*crLHS183 + crLHS210 + crLHS336);
rLHS(14,3)+=w_gauss*(DN(0,2)*crLHS327 + crLHS211*crLHS69 - crLHS212);
rLHS(14,4)+=w_gauss*(DN(3,0)*crLHS45 + DN(3,1)*crLHS137 + DN(3,2)*crLHS185 + crLHS245);
rLHS(14,5)+=w_gauss*(DN(3,0)*crLHS59 + DN(3,1)*crLHS141 + DN(3,2)*crLHS188 + crLHS263);
rLHS(14,6)+=w_gauss*(DN(3,0)*crLHS65 + DN(3,1)*crLHS147 + DN(3,2)*crLHS190 + crLHS277 + crLHS337);
rLHS(14,7)+=w_gauss*(DN(1,2)*crLHS327 + crLHS278*crLHS69 - crLHS279);
rLHS(14,8)+=w_gauss*(DN(3,0)*crLHS74 + DN(3,1)*crLHS152 + DN(3,2)*crLHS195 + crLHS299);
rLHS(14,9)+=w_gauss*(DN(3,0)*crLHS88 + DN(3,1)*crLHS156 + DN(3,2)*crLHS197 + crLHS311);
rLHS(14,10)+=w_gauss*(DN(3,0)*crLHS94 + DN(3,1)*crLHS162 + DN(3,2)*crLHS199 + crLHS319 + crLHS338);
rLHS(14,11)+=w_gauss*(DN(2,2)*crLHS327 + crLHS320*crLHS69 - crLHS321);
rLHS(14,12)+=w_gauss*(DN(3,0)*crLHS102 + DN(3,1)*crLHS167 + DN(3,2)*crLHS204 + crLHS334);
rLHS(14,13)+=w_gauss*(DN(3,0)*crLHS115 + DN(3,1)*crLHS171 + DN(3,2)*crLHS206 + crLHS340);
rLHS(14,14)+=w_gauss*(DN(3,0)*crLHS121 + DN(3,1)*crLHS177 + DN(3,2)*crLHS208 + crLHS10*crLHS341 + crLHS331);
rLHS(14,15)+=DN(3,2)*crLHS335;
rLHS(15,0)+=w_gauss*(DN(3,0)*crLHS213 + crLHS123);
rLHS(15,1)+=w_gauss*(DN(3,1)*crLHS213 + crLHS179);
rLHS(15,2)+=w_gauss*(DN(3,2)*crLHS213 + crLHS211);
rLHS(15,3)+=crLHS221;
rLHS(15,4)+=w_gauss*(DN(3,0)*crLHS216 + crLHS246);
rLHS(15,5)+=w_gauss*(DN(3,1)*crLHS216 + crLHS264);
rLHS(15,6)+=w_gauss*(DN(3,2)*crLHS216 + crLHS278);
rLHS(15,7)+=crLHS282;
rLHS(15,8)+=w_gauss*(DN(3,0)*crLHS218 + crLHS300);
rLHS(15,9)+=w_gauss*(DN(3,1)*crLHS218 + crLHS312);
rLHS(15,10)+=w_gauss*(DN(3,2)*crLHS218 + crLHS320);
rLHS(15,11)+=crLHS323;
rLHS(15,12)+=DN(3,0)*crLHS342;
rLHS(15,13)+=DN(3,1)*crLHS342;
rLHS(15,14)+=DN(3,2)*crLHS342;
rLHS(15,15)+=crLHS215*(crLHS330 + crLHS339 + crLHS341);

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

    // Add RHS Gauss point contribution
    const double w_gauss = rData.Weight;

    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crRHS1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crRHS2 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0));
const double crRHS3 = N[0]*rho;
const double crRHS4 = bdf0*crRHS3;
const double crRHS5 = N[0]*(vn(0,0) - vnn(0,0));
const double crRHS6 = N[1]*(vn(1,0) - vnn(1,0));
const double crRHS7 = N[2]*(vn(2,0) - vnn(2,0));
const double crRHS8 = crRHS5 + crRHS6 + crRHS7;
const double crRHS9 = 1.0/dt;
const double crRHS10 = crRHS3*crRHS9;
const double crRHS11 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS12 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS14 = rho*(crRHS11*crRHS12 + crRHS13*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crRHS15 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crRHS16 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crRHS17 = crRHS15*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crRHS16*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0));
const double crRHS18 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crRHS19 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crRHS20 = rho*(crRHS18*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crRHS19*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0)));
const double crRHS21 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS22 = crRHS11 + crRHS21 - volume_error_ratio;
const double crRHS23 = rho*stab_c2*sqrt(crRHS12*crRHS12 + crRHS13*crRHS13);
const double crRHS24 = crRHS22*(crRHS23*h*1.0/stab_c1 + mu);
const double crRHS25 = bdf0*rho;
const double crRHS26 = crRHS2*crRHS25;
const double crRHS27 = crRHS9*rho;
const double crRHS28 = 1.0*1.0/(crRHS23*1.0/h + crRHS27*dyn_tau + mu*stab_c1*1.0/(h*h));
const double crRHS29 = crRHS28*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS1 + crRHS14 - crRHS20 + crRHS26 + rho*(crRHS17 + crRHS5*crRHS9 + crRHS6*crRHS9 + crRHS7*crRHS9));
const double crRHS30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crRHS31 = crRHS3*crRHS30;
const double crRHS32 = rho*(DN(0,0)*crRHS12 + DN(0,1)*crRHS13);
const double crRHS33 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS34 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1));
const double crRHS35 = N[0]*(vn(0,1) - vnn(0,1));
const double crRHS36 = N[1]*(vn(1,1) - vnn(1,1));
const double crRHS37 = N[2]*(vn(2,1) - vnn(2,1));
const double crRHS38 = crRHS35 + crRHS36 + crRHS37;
const double crRHS39 = rho*(crRHS12*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crRHS13*crRHS21);
const double crRHS40 = crRHS15*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crRHS16*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1));
const double crRHS41 = rho*(crRHS18*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crRHS19*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1)));
const double crRHS42 = crRHS25*crRHS34;
const double crRHS43 = crRHS28*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS33 + crRHS39 - crRHS41 + crRHS42 + rho*(crRHS35*crRHS9 + crRHS36*crRHS9 + crRHS37*crRHS9 + crRHS40));
const double crRHS44 = N[1]*crRHS27;
const double crRHS45 = N[1]*rho;
const double crRHS46 = crRHS30*crRHS45;
const double crRHS47 = rho*(DN(1,0)*crRHS12 + DN(1,1)*crRHS13);
const double crRHS48 = N[2]*crRHS27;
const double crRHS49 = N[2]*rho;
const double crRHS50 = crRHS30*crRHS49;
const double crRHS51 = rho*(DN(2,0)*crRHS12 + DN(2,1)*crRHS13);
rRHS[0]+=-w_gauss*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS24 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS1 + N[0]*crRHS14 - N[0]*crRHS20 + crRHS10*crRHS8 + crRHS17*crRHS3 + crRHS2*crRHS4 + crRHS29*crRHS31 + crRHS29*crRHS32);
rRHS[1]+=-w_gauss*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS24 + DN(0,1)*stress[1] - N[0]*crRHS33 + N[0]*crRHS39 - N[0]*crRHS41 + crRHS10*crRHS38 + crRHS3*crRHS40 + crRHS31*crRHS43 + crRHS32*crRHS43 + crRHS34*crRHS4);
rRHS[2]+=-w_gauss*(DN(0,0)*crRHS29 + DN(0,1)*crRHS43 + N[0]*crRHS22);
rRHS[3]+=-w_gauss*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS24 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS1 + N[1]*crRHS14 - N[1]*crRHS20 + N[1]*crRHS26 + crRHS17*crRHS45 + crRHS29*crRHS46 + crRHS29*crRHS47 + crRHS44*crRHS8);
rRHS[4]+=-w_gauss*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS24 + DN(1,1)*stress[1] - N[1]*crRHS33 + N[1]*crRHS39 - N[1]*crRHS41 + N[1]*crRHS42 + crRHS38*crRHS44 + crRHS40*crRHS45 + crRHS43*crRHS46 + crRHS43*crRHS47);
rRHS[5]+=-w_gauss*(DN(1,0)*crRHS29 + DN(1,1)*crRHS43 + N[1]*crRHS22);
rRHS[6]+=-w_gauss*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS24 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS1 + N[2]*crRHS14 - N[2]*crRHS20 + N[2]*crRHS26 + crRHS17*crRHS49 + crRHS29*crRHS50 + crRHS29*crRHS51 + crRHS48*crRHS8);
rRHS[7]+=-w_gauss*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS24 + DN(2,1)*stress[1] - N[2]*crRHS33 + N[2]*crRHS39 - N[2]*crRHS41 + N[2]*crRHS42 + crRHS38*crRHS48 + crRHS40*crRHS49 + crRHS43*crRHS50 + crRHS43*crRHS51);
rRHS[8]+=-w_gauss*(DN(2,0)*crRHS29 + DN(2,1)*crRHS43 + N[2]*crRHS22);

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

    // Add RHS Gauss point contribution
    const double w_gauss = rData.Weight;

    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crRHS1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS2 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0)) + N[3]*(v(3,0) - vfrac(3,0));
const double crRHS3 = N[0]*rho;
const double crRHS4 = bdf0*crRHS3;
const double crRHS5 = N[0]*(vn(0,0) - vnn(0,0));
const double crRHS6 = N[1]*(vn(1,0) - vnn(1,0));
const double crRHS7 = N[2]*(vn(2,0) - vnn(2,0));
const double crRHS8 = N[3]*(vn(3,0) - vnn(3,0));
const double crRHS9 = crRHS5 + crRHS6 + crRHS7 + crRHS8;
const double crRHS10 = 1.0/dt;
const double crRHS11 = crRHS10*crRHS3;
const double crRHS12 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS13 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS14 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS15 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS16 = rho*(crRHS12*crRHS13 + crRHS14*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crRHS15*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crRHS17 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crRHS18 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crRHS19 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crRHS20 = crRHS17*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crRHS18*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crRHS19*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0));
const double crRHS21 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crRHS22 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crRHS23 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crRHS24 = rho*(crRHS21*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crRHS22*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0)) + crRHS23*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0)));
const double crRHS25 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS26 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS27 = crRHS12 + crRHS25 + crRHS26 - volume_error_ratio;
const double crRHS28 = rho*stab_c2*sqrt(crRHS13*crRHS13 + crRHS14*crRHS14 + crRHS15*crRHS15);
const double crRHS29 = crRHS27*(crRHS28*h*1.0/stab_c1 + mu);
const double crRHS30 = bdf0*rho;
const double crRHS31 = crRHS2*crRHS30;
const double crRHS32 = crRHS10*rho;
const double crRHS33 = 1.0*1.0/(crRHS28*1.0/h + crRHS32*dyn_tau + mu*stab_c1*1.0/(h*h));
const double crRHS34 = crRHS33*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crRHS1 + crRHS16 - crRHS24 + crRHS31 + rho*(crRHS10*crRHS5 + crRHS10*crRHS6 + crRHS10*crRHS7 + crRHS10*crRHS8 + crRHS20));
const double crRHS35 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crRHS36 = crRHS3*crRHS35;
const double crRHS37 = rho*(DN(0,0)*crRHS13 + DN(0,1)*crRHS14 + DN(0,2)*crRHS15);
const double crRHS38 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS39 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1)) + N[3]*(v(3,1) - vfrac(3,1));
const double crRHS40 = N[0]*(vn(0,1) - vnn(0,1));
const double crRHS41 = N[1]*(vn(1,1) - vnn(1,1));
const double crRHS42 = N[2]*(vn(2,1) - vnn(2,1));
const double crRHS43 = N[3]*(vn(3,1) - vnn(3,1));
const double crRHS44 = crRHS40 + crRHS41 + crRHS42 + crRHS43;
const double crRHS45 = rho*(crRHS13*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crRHS14*crRHS25 + crRHS15*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crRHS46 = crRHS17*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crRHS18*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crRHS19*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1));
const double crRHS47 = rho*(crRHS21*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crRHS22*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1)) + crRHS23*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1)));
const double crRHS48 = crRHS30*crRHS39;
const double crRHS49 = crRHS33*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS38 + crRHS45 - crRHS47 + crRHS48 + rho*(crRHS10*crRHS40 + crRHS10*crRHS41 + crRHS10*crRHS42 + crRHS10*crRHS43 + crRHS46));
const double crRHS50 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS51 = N[0]*(v(0,2) - vfrac(0,2)) + N[1]*(v(1,2) - vfrac(1,2)) + N[2]*(v(2,2) - vfrac(2,2)) + N[3]*(v(3,2) - vfrac(3,2));
const double crRHS52 = N[0]*(vn(0,2) - vnn(0,2));
const double crRHS53 = N[1]*(vn(1,2) - vnn(1,2));
const double crRHS54 = N[2]*(vn(2,2) - vnn(2,2));
const double crRHS55 = N[3]*(vn(3,2) - vnn(3,2));
const double crRHS56 = crRHS52 + crRHS53 + crRHS54 + crRHS55;
const double crRHS57 = rho*(crRHS13*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crRHS14*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crRHS15*crRHS26);
const double crRHS58 = crRHS17*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crRHS18*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crRHS19*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2));
const double crRHS59 = rho*(crRHS21*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crRHS22*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2)) + crRHS23*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2)));
const double crRHS60 = crRHS30*crRHS51;
const double crRHS61 = crRHS33*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS50 + crRHS57 - crRHS59 + crRHS60 + rho*(crRHS10*crRHS52 + crRHS10*crRHS53 + crRHS10*crRHS54 + crRHS10*crRHS55 + crRHS58));
const double crRHS62 = N[1]*crRHS32;
const double crRHS63 = N[1]*rho;
const double crRHS64 = crRHS35*crRHS63;
const double crRHS65 = rho*(DN(1,0)*crRHS13 + DN(1,1)*crRHS14 + DN(1,2)*crRHS15);
const double crRHS66 = N[2]*crRHS32;
const double crRHS67 = N[2]*rho;
const double crRHS68 = crRHS35*crRHS67;
const double crRHS69 = rho*(DN(2,0)*crRHS13 + DN(2,1)*crRHS14 + DN(2,2)*crRHS15);
const double crRHS70 = N[3]*crRHS32;
const double crRHS71 = N[3]*rho;
const double crRHS72 = crRHS35*crRHS71;
const double crRHS73 = rho*(DN(3,0)*crRHS13 + DN(3,1)*crRHS14 + DN(3,2)*crRHS15);
rRHS[0]+=-w_gauss*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS29 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] - N[0]*crRHS1 + N[0]*crRHS16 - N[0]*crRHS24 + crRHS11*crRHS9 + crRHS2*crRHS4 + crRHS20*crRHS3 + crRHS34*crRHS36 + crRHS34*crRHS37);
rRHS[1]+=-w_gauss*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS29 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS38 + N[0]*crRHS45 - N[0]*crRHS47 + crRHS11*crRHS44 + crRHS3*crRHS46 + crRHS36*crRHS49 + crRHS37*crRHS49 + crRHS39*crRHS4);
rRHS[2]+=-w_gauss*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS29 + DN(0,2)*stress[2] - N[0]*crRHS50 + N[0]*crRHS57 - N[0]*crRHS59 + crRHS11*crRHS56 + crRHS3*crRHS58 + crRHS36*crRHS61 + crRHS37*crRHS61 + crRHS4*crRHS51);
rRHS[3]+=-w_gauss*(DN(0,0)*crRHS34 + DN(0,1)*crRHS49 + DN(0,2)*crRHS61 + N[0]*crRHS27);
rRHS[4]+=-w_gauss*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS29 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] - N[1]*crRHS1 + N[1]*crRHS16 - N[1]*crRHS24 + N[1]*crRHS31 + crRHS20*crRHS63 + crRHS34*crRHS64 + crRHS34*crRHS65 + crRHS62*crRHS9);
rRHS[5]+=-w_gauss*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS29 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS38 + N[1]*crRHS45 - N[1]*crRHS47 + N[1]*crRHS48 + crRHS44*crRHS62 + crRHS46*crRHS63 + crRHS49*crRHS64 + crRHS49*crRHS65);
rRHS[6]+=-w_gauss*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS29 + DN(1,2)*stress[2] - N[1]*crRHS50 + N[1]*crRHS57 - N[1]*crRHS59 + N[1]*crRHS60 + crRHS56*crRHS62 + crRHS58*crRHS63 + crRHS61*crRHS64 + crRHS61*crRHS65);
rRHS[7]+=-w_gauss*(DN(1,0)*crRHS34 + DN(1,1)*crRHS49 + DN(1,2)*crRHS61 + N[1]*crRHS27);
rRHS[8]+=-w_gauss*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS29 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] - N[2]*crRHS1 + N[2]*crRHS16 - N[2]*crRHS24 + N[2]*crRHS31 + crRHS20*crRHS67 + crRHS34*crRHS68 + crRHS34*crRHS69 + crRHS66*crRHS9);
rRHS[9]+=-w_gauss*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS29 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS38 + N[2]*crRHS45 - N[2]*crRHS47 + N[2]*crRHS48 + crRHS44*crRHS66 + crRHS46*crRHS67 + crRHS49*crRHS68 + crRHS49*crRHS69);
rRHS[10]+=-w_gauss*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS29 + DN(2,2)*stress[2] - N[2]*crRHS50 + N[2]*crRHS57 - N[2]*crRHS59 + N[2]*crRHS60 + crRHS56*crRHS66 + crRHS58*crRHS67 + crRHS61*crRHS68 + crRHS61*crRHS69);
rRHS[11]+=-w_gauss*(DN(2,0)*crRHS34 + DN(2,1)*crRHS49 + DN(2,2)*crRHS61 + N[2]*crRHS27);
rRHS[12]+=-w_gauss*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS29 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] - N[3]*crRHS1 + N[3]*crRHS16 - N[3]*crRHS24 + N[3]*crRHS31 + crRHS20*crRHS71 + crRHS34*crRHS72 + crRHS34*crRHS73 + crRHS70*crRHS9);
rRHS[13]+=-w_gauss*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS29 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS38 + N[3]*crRHS45 - N[3]*crRHS47 + N[3]*crRHS48 + crRHS44*crRHS70 + crRHS46*crRHS71 + crRHS49*crRHS72 + crRHS49*crRHS73);
rRHS[14]+=-w_gauss*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS29 + DN(3,2)*stress[2] - N[3]*crRHS50 + N[3]*crRHS57 - N[3]*crRHS59 + N[3]*crRHS60 + crRHS56*crRHS70 + crRHS58*crRHS71 + crRHS61*crRHS72 + crRHS61*crRHS73);
rRHS[15]+=-w_gauss*(DN(3,0)*crRHS34 + DN(3,1)*crRHS49 + DN(3,2)*crRHS61 + N[3]*crRHS27);

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

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    // Add Enrichment Gauss point contribution
    const double w_gauss = rData.Weight;

    const double crV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crV2 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crV0*crV0 + crV1*crV1));
const double crV3 = crV2*rho;
const double crV4 = crV3*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crV5 = N[0]*crV4;
const double crV6 = crV3*(DN(0,0)*crV0 + DN(0,1)*crV1);
const double crV7 = crV2*w_gauss;
const double crV8 = N[1]*crV4;
const double crV9 = crV3*(DN(1,0)*crV0 + DN(1,1)*crV1);
const double crV10 = N[2]*crV4;
const double crV11 = crV3*(DN(2,0)*crV0 + DN(2,1)*crV1);
rV(0,0)+=w_gauss*(-DN(0,0)*Nenr[0] + DNenr(0,0)*crV5 + DNenr(0,0)*crV6);
rV(0,1)+=w_gauss*(-DN(0,0)*Nenr[1] + DNenr(1,0)*crV5 + DNenr(1,0)*crV6);
rV(0,2)+=w_gauss*(-DN(0,0)*Nenr[2] + DNenr(2,0)*crV5 + DNenr(2,0)*crV6);
rV(1,0)+=w_gauss*(-DN(0,1)*Nenr[0] + DNenr(0,1)*crV5 + DNenr(0,1)*crV6);
rV(1,1)+=w_gauss*(-DN(0,1)*Nenr[1] + DNenr(1,1)*crV5 + DNenr(1,1)*crV6);
rV(1,2)+=w_gauss*(-DN(0,1)*Nenr[2] + DNenr(2,1)*crV5 + DNenr(2,1)*crV6);
rV(2,0)+=crV7*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
rV(2,1)+=crV7*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
rV(2,2)+=crV7*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
rV(3,0)+=w_gauss*(-DN(1,0)*Nenr[0] + DNenr(0,0)*crV8 + DNenr(0,0)*crV9);
rV(3,1)+=w_gauss*(-DN(1,0)*Nenr[1] + DNenr(1,0)*crV8 + DNenr(1,0)*crV9);
rV(3,2)+=w_gauss*(-DN(1,0)*Nenr[2] + DNenr(2,0)*crV8 + DNenr(2,0)*crV9);
rV(4,0)+=w_gauss*(-DN(1,1)*Nenr[0] + DNenr(0,1)*crV8 + DNenr(0,1)*crV9);
rV(4,1)+=w_gauss*(-DN(1,1)*Nenr[1] + DNenr(1,1)*crV8 + DNenr(1,1)*crV9);
rV(4,2)+=w_gauss*(-DN(1,1)*Nenr[2] + DNenr(2,1)*crV8 + DNenr(2,1)*crV9);
rV(5,0)+=crV7*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
rV(5,1)+=crV7*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
rV(5,2)+=crV7*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
rV(6,0)+=w_gauss*(-DN(2,0)*Nenr[0] + DNenr(0,0)*crV10 + DNenr(0,0)*crV11);
rV(6,1)+=w_gauss*(-DN(2,0)*Nenr[1] + DNenr(1,0)*crV10 + DNenr(1,0)*crV11);
rV(6,2)+=w_gauss*(-DN(2,0)*Nenr[2] + DNenr(2,0)*crV10 + DNenr(2,0)*crV11);
rV(7,0)+=w_gauss*(-DN(2,1)*Nenr[0] + DNenr(0,1)*crV10 + DNenr(0,1)*crV11);
rV(7,1)+=w_gauss*(-DN(2,1)*Nenr[1] + DNenr(1,1)*crV10 + DNenr(1,1)*crV11);
rV(7,2)+=w_gauss*(-DN(2,1)*Nenr[2] + DNenr(2,1)*crV10 + DNenr(2,1)*crV11);
rV(8,0)+=crV7*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
rV(8,1)+=crV7*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
rV(8,2)+=crV7*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double crH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crH2 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crH0*crH0 + crH1*crH1));
const double crH3 = crH2*rho;
const double crH4 = crH3*(DN(0,0)*crH0 + DN(0,1)*crH1 + N[0]*bdf0);
const double crH5 = crH2*w_gauss;
const double crH6 = crH3*(DN(1,0)*crH0 + DN(1,1)*crH1 + N[1]*bdf0);
const double crH7 = crH3*(DN(2,0)*crH0 + DN(2,1)*crH1 + N[2]*bdf0);
rH(0,0)+=w_gauss*(DN(0,0)*Nenr[0] + DNenr(0,0)*crH4);
rH(0,1)+=w_gauss*(DN(0,1)*Nenr[0] + DNenr(0,1)*crH4);
rH(0,2)+=crH5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
rH(0,3)+=w_gauss*(DN(1,0)*Nenr[0] + DNenr(0,0)*crH6);
rH(0,4)+=w_gauss*(DN(1,1)*Nenr[0] + DNenr(0,1)*crH6);
rH(0,5)+=crH5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
rH(0,6)+=w_gauss*(DN(2,0)*Nenr[0] + DNenr(0,0)*crH7);
rH(0,7)+=w_gauss*(DN(2,1)*Nenr[0] + DNenr(0,1)*crH7);
rH(0,8)+=crH5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
rH(1,0)+=w_gauss*(DN(0,0)*Nenr[1] + DNenr(1,0)*crH4);
rH(1,1)+=w_gauss*(DN(0,1)*Nenr[1] + DNenr(1,1)*crH4);
rH(1,2)+=crH5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
rH(1,3)+=w_gauss*(DN(1,0)*Nenr[1] + DNenr(1,0)*crH6);
rH(1,4)+=w_gauss*(DN(1,1)*Nenr[1] + DNenr(1,1)*crH6);
rH(1,5)+=crH5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
rH(1,6)+=w_gauss*(DN(2,0)*Nenr[1] + DNenr(1,0)*crH7);
rH(1,7)+=w_gauss*(DN(2,1)*Nenr[1] + DNenr(1,1)*crH7);
rH(1,8)+=crH5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
rH(2,0)+=w_gauss*(DN(0,0)*Nenr[2] + DNenr(2,0)*crH4);
rH(2,1)+=w_gauss*(DN(0,1)*Nenr[2] + DNenr(2,1)*crH4);
rH(2,2)+=crH5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
rH(2,3)+=w_gauss*(DN(1,0)*Nenr[2] + DNenr(2,0)*crH6);
rH(2,4)+=w_gauss*(DN(1,1)*Nenr[2] + DNenr(2,1)*crH6);
rH(2,5)+=crH5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
rH(2,6)+=w_gauss*(DN(2,0)*Nenr[2] + DNenr(2,0)*crH7);
rH(2,7)+=w_gauss*(DN(2,1)*Nenr[2] + DNenr(2,1)*crH7);
rH(2,8)+=crH5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double crKee0 = 1.0*w_gauss*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))*1.0/h);
const double crKee1 = crKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double crKee2 = crKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double crKee3 = crKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
rKee(0,0)+=crKee0*(DNenr(0,0)*DNenr(0,0) + DNenr(0,1)*DNenr(0,1));
rKee(0,1)+=crKee1;
rKee(0,2)+=crKee2;
rKee(1,0)+=crKee1;
rKee(1,1)+=crKee0*(DNenr(1,0)*DNenr(1,0) + DNenr(1,1)*DNenr(1,1));
rKee(1,2)+=crKee3;
rKee(2,0)+=crKee2;
rKee(2,1)+=crKee3;
rKee(2,2)+=crKee0*(DNenr(2,0)*DNenr(2,0) + DNenr(2,1)*DNenr(2,1));


    const double crRHS_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS_ee2 = crRHS_ee0 + crRHS_ee1 - volume_error_ratio;
const double crRHS_ee3 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crRHS_ee4 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crRHS_ee5 = N[0]*bdf0;
const double crRHS_ee6 = N[1]*bdf0;
const double crRHS_ee7 = N[2]*bdf0;
const double crRHS_ee8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS_ee9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS_ee10 = 1.0/dt;
const double crRHS_ee11 = N[0]*crRHS_ee10;
const double crRHS_ee12 = N[1]*crRHS_ee10;
const double crRHS_ee13 = N[2]*crRHS_ee10;
const double crRHS_ee14 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crRHS_ee15 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crRHS_ee16 = 1.0*1.0/(crRHS_ee10*dyn_tau*rho + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crRHS_ee8*crRHS_ee8 + crRHS_ee9*crRHS_ee9));
const double crRHS_ee17 = crRHS_ee16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(crRHS_ee3*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crRHS_ee4*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(crRHS_ee0*crRHS_ee8 + crRHS_ee5*(v(0,0) - vfrac(0,0)) + crRHS_ee6*(v(1,0) - vfrac(1,0)) + crRHS_ee7*(v(2,0) - vfrac(2,0)) + crRHS_ee9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))) + rho*(crRHS_ee11*(vn(0,0) - vnn(0,0)) + crRHS_ee12*(vn(1,0) - vnn(1,0)) + crRHS_ee13*(vn(2,0) - vnn(2,0)) + crRHS_ee14*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crRHS_ee15*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0))));
const double crRHS_ee18 = crRHS_ee16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(crRHS_ee3*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crRHS_ee4*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(crRHS_ee1*crRHS_ee9 + crRHS_ee5*(v(0,1) - vfrac(0,1)) + crRHS_ee6*(v(1,1) - vfrac(1,1)) + crRHS_ee7*(v(2,1) - vfrac(2,1)) + crRHS_ee8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))) + rho*(crRHS_ee11*(vn(0,1) - vnn(0,1)) + crRHS_ee12*(vn(1,1) - vnn(1,1)) + crRHS_ee13*(vn(2,1) - vnn(2,1)) + crRHS_ee14*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crRHS_ee15*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1))));
rRHS_ee[0]+=-w_gauss*(DNenr(0,0)*crRHS_ee17 + DNenr(0,1)*crRHS_ee18 + Nenr[0]*crRHS_ee2);
rRHS_ee[1]+=-w_gauss*(DNenr(1,0)*crRHS_ee17 + DNenr(1,1)*crRHS_ee18 + Nenr[1]*crRHS_ee2);
rRHS_ee[2]+=-w_gauss*(DNenr(2,0)*crRHS_ee17 + DNenr(2,1)*crRHS_ee18 + Nenr[2]*crRHS_ee2);


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

    // Add Enrichment Gauss point contribution
    const double w_gauss = rData.Weight;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double crV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crV2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crV3 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crV0*crV0 + crV1*crV1 + crV2*crV2));
const double crV4 = crV3*rho;
const double crV5 = crV4*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crV6 = N[0]*crV5;
const double crV7 = crV4*(DN(0,0)*crV0 + DN(0,1)*crV1 + DN(0,2)*crV2);
const double crV8 = crV3*w_gauss;
const double crV9 = N[1]*crV5;
const double crV10 = crV4*(DN(1,0)*crV0 + DN(1,1)*crV1 + DN(1,2)*crV2);
const double crV11 = N[2]*crV5;
const double crV12 = crV4*(DN(2,0)*crV0 + DN(2,1)*crV1 + DN(2,2)*crV2);
const double crV13 = N[3]*crV5;
const double crV14 = crV4*(DN(3,0)*crV0 + DN(3,1)*crV1 + DN(3,2)*crV2);
rV(0,0)+=w_gauss*(-DN(0,0)*Nenr[0] + DNenr(0,0)*crV6 + DNenr(0,0)*crV7);
rV(0,1)+=w_gauss*(-DN(0,0)*Nenr[1] + DNenr(1,0)*crV6 + DNenr(1,0)*crV7);
rV(0,2)+=w_gauss*(-DN(0,0)*Nenr[2] + DNenr(2,0)*crV6 + DNenr(2,0)*crV7);
rV(0,3)+=w_gauss*(-DN(0,0)*Nenr[3] + DNenr(3,0)*crV6 + DNenr(3,0)*crV7);
rV(1,0)+=w_gauss*(-DN(0,1)*Nenr[0] + DNenr(0,1)*crV6 + DNenr(0,1)*crV7);
rV(1,1)+=w_gauss*(-DN(0,1)*Nenr[1] + DNenr(1,1)*crV6 + DNenr(1,1)*crV7);
rV(1,2)+=w_gauss*(-DN(0,1)*Nenr[2] + DNenr(2,1)*crV6 + DNenr(2,1)*crV7);
rV(1,3)+=w_gauss*(-DN(0,1)*Nenr[3] + DNenr(3,1)*crV6 + DNenr(3,1)*crV7);
rV(2,0)+=w_gauss*(-DN(0,2)*Nenr[0] + DNenr(0,2)*crV6 + DNenr(0,2)*crV7);
rV(2,1)+=w_gauss*(-DN(0,2)*Nenr[1] + DNenr(1,2)*crV6 + DNenr(1,2)*crV7);
rV(2,2)+=w_gauss*(-DN(0,2)*Nenr[2] + DNenr(2,2)*crV6 + DNenr(2,2)*crV7);
rV(2,3)+=w_gauss*(-DN(0,2)*Nenr[3] + DNenr(3,2)*crV6 + DNenr(3,2)*crV7);
rV(3,0)+=crV8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
rV(3,1)+=crV8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
rV(3,2)+=crV8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
rV(3,3)+=crV8*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
rV(4,0)+=w_gauss*(-DN(1,0)*Nenr[0] + DNenr(0,0)*crV10 + DNenr(0,0)*crV9);
rV(4,1)+=w_gauss*(-DN(1,0)*Nenr[1] + DNenr(1,0)*crV10 + DNenr(1,0)*crV9);
rV(4,2)+=w_gauss*(-DN(1,0)*Nenr[2] + DNenr(2,0)*crV10 + DNenr(2,0)*crV9);
rV(4,3)+=w_gauss*(-DN(1,0)*Nenr[3] + DNenr(3,0)*crV10 + DNenr(3,0)*crV9);
rV(5,0)+=w_gauss*(-DN(1,1)*Nenr[0] + DNenr(0,1)*crV10 + DNenr(0,1)*crV9);
rV(5,1)+=w_gauss*(-DN(1,1)*Nenr[1] + DNenr(1,1)*crV10 + DNenr(1,1)*crV9);
rV(5,2)+=w_gauss*(-DN(1,1)*Nenr[2] + DNenr(2,1)*crV10 + DNenr(2,1)*crV9);
rV(5,3)+=w_gauss*(-DN(1,1)*Nenr[3] + DNenr(3,1)*crV10 + DNenr(3,1)*crV9);
rV(6,0)+=w_gauss*(-DN(1,2)*Nenr[0] + DNenr(0,2)*crV10 + DNenr(0,2)*crV9);
rV(6,1)+=w_gauss*(-DN(1,2)*Nenr[1] + DNenr(1,2)*crV10 + DNenr(1,2)*crV9);
rV(6,2)+=w_gauss*(-DN(1,2)*Nenr[2] + DNenr(2,2)*crV10 + DNenr(2,2)*crV9);
rV(6,3)+=w_gauss*(-DN(1,2)*Nenr[3] + DNenr(3,2)*crV10 + DNenr(3,2)*crV9);
rV(7,0)+=crV8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
rV(7,1)+=crV8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
rV(7,2)+=crV8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
rV(7,3)+=crV8*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
rV(8,0)+=w_gauss*(-DN(2,0)*Nenr[0] + DNenr(0,0)*crV11 + DNenr(0,0)*crV12);
rV(8,1)+=w_gauss*(-DN(2,0)*Nenr[1] + DNenr(1,0)*crV11 + DNenr(1,0)*crV12);
rV(8,2)+=w_gauss*(-DN(2,0)*Nenr[2] + DNenr(2,0)*crV11 + DNenr(2,0)*crV12);
rV(8,3)+=w_gauss*(-DN(2,0)*Nenr[3] + DNenr(3,0)*crV11 + DNenr(3,0)*crV12);
rV(9,0)+=w_gauss*(-DN(2,1)*Nenr[0] + DNenr(0,1)*crV11 + DNenr(0,1)*crV12);
rV(9,1)+=w_gauss*(-DN(2,1)*Nenr[1] + DNenr(1,1)*crV11 + DNenr(1,1)*crV12);
rV(9,2)+=w_gauss*(-DN(2,1)*Nenr[2] + DNenr(2,1)*crV11 + DNenr(2,1)*crV12);
rV(9,3)+=w_gauss*(-DN(2,1)*Nenr[3] + DNenr(3,1)*crV11 + DNenr(3,1)*crV12);
rV(10,0)+=w_gauss*(-DN(2,2)*Nenr[0] + DNenr(0,2)*crV11 + DNenr(0,2)*crV12);
rV(10,1)+=w_gauss*(-DN(2,2)*Nenr[1] + DNenr(1,2)*crV11 + DNenr(1,2)*crV12);
rV(10,2)+=w_gauss*(-DN(2,2)*Nenr[2] + DNenr(2,2)*crV11 + DNenr(2,2)*crV12);
rV(10,3)+=w_gauss*(-DN(2,2)*Nenr[3] + DNenr(3,2)*crV11 + DNenr(3,2)*crV12);
rV(11,0)+=crV8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
rV(11,1)+=crV8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
rV(11,2)+=crV8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
rV(11,3)+=crV8*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
rV(12,0)+=w_gauss*(-DN(3,0)*Nenr[0] + DNenr(0,0)*crV13 + DNenr(0,0)*crV14);
rV(12,1)+=w_gauss*(-DN(3,0)*Nenr[1] + DNenr(1,0)*crV13 + DNenr(1,0)*crV14);
rV(12,2)+=w_gauss*(-DN(3,0)*Nenr[2] + DNenr(2,0)*crV13 + DNenr(2,0)*crV14);
rV(12,3)+=w_gauss*(-DN(3,0)*Nenr[3] + DNenr(3,0)*crV13 + DNenr(3,0)*crV14);
rV(13,0)+=w_gauss*(-DN(3,1)*Nenr[0] + DNenr(0,1)*crV13 + DNenr(0,1)*crV14);
rV(13,1)+=w_gauss*(-DN(3,1)*Nenr[1] + DNenr(1,1)*crV13 + DNenr(1,1)*crV14);
rV(13,2)+=w_gauss*(-DN(3,1)*Nenr[2] + DNenr(2,1)*crV13 + DNenr(2,1)*crV14);
rV(13,3)+=w_gauss*(-DN(3,1)*Nenr[3] + DNenr(3,1)*crV13 + DNenr(3,1)*crV14);
rV(14,0)+=w_gauss*(-DN(3,2)*Nenr[0] + DNenr(0,2)*crV13 + DNenr(0,2)*crV14);
rV(14,1)+=w_gauss*(-DN(3,2)*Nenr[1] + DNenr(1,2)*crV13 + DNenr(1,2)*crV14);
rV(14,2)+=w_gauss*(-DN(3,2)*Nenr[2] + DNenr(2,2)*crV13 + DNenr(2,2)*crV14);
rV(14,3)+=w_gauss*(-DN(3,2)*Nenr[3] + DNenr(3,2)*crV13 + DNenr(3,2)*crV14);
rV(15,0)+=crV8*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
rV(15,1)+=crV8*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
rV(15,2)+=crV8*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
rV(15,3)+=crV8*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double crH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crH2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crH3 = 1.0*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crH0*crH0 + crH1*crH1 + crH2*crH2));
const double crH4 = crH3*rho;
const double crH5 = crH4*(DN(0,0)*crH0 + DN(0,1)*crH1 + DN(0,2)*crH2 + N[0]*bdf0);
const double crH6 = crH3*w_gauss;
const double crH7 = crH4*(DN(1,0)*crH0 + DN(1,1)*crH1 + DN(1,2)*crH2 + N[1]*bdf0);
const double crH8 = crH4*(DN(2,0)*crH0 + DN(2,1)*crH1 + DN(2,2)*crH2 + N[2]*bdf0);
const double crH9 = crH4*(DN(3,0)*crH0 + DN(3,1)*crH1 + DN(3,2)*crH2 + N[3]*bdf0);
rH(0,0)+=w_gauss*(DN(0,0)*Nenr[0] + DNenr(0,0)*crH5);
rH(0,1)+=w_gauss*(DN(0,1)*Nenr[0] + DNenr(0,1)*crH5);
rH(0,2)+=w_gauss*(DN(0,2)*Nenr[0] + DNenr(0,2)*crH5);
rH(0,3)+=crH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
rH(0,4)+=w_gauss*(DN(1,0)*Nenr[0] + DNenr(0,0)*crH7);
rH(0,5)+=w_gauss*(DN(1,1)*Nenr[0] + DNenr(0,1)*crH7);
rH(0,6)+=w_gauss*(DN(1,2)*Nenr[0] + DNenr(0,2)*crH7);
rH(0,7)+=crH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
rH(0,8)+=w_gauss*(DN(2,0)*Nenr[0] + DNenr(0,0)*crH8);
rH(0,9)+=w_gauss*(DN(2,1)*Nenr[0] + DNenr(0,1)*crH8);
rH(0,10)+=w_gauss*(DN(2,2)*Nenr[0] + DNenr(0,2)*crH8);
rH(0,11)+=crH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
rH(0,12)+=w_gauss*(DN(3,0)*Nenr[0] + DNenr(0,0)*crH9);
rH(0,13)+=w_gauss*(DN(3,1)*Nenr[0] + DNenr(0,1)*crH9);
rH(0,14)+=w_gauss*(DN(3,2)*Nenr[0] + DNenr(0,2)*crH9);
rH(0,15)+=crH6*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
rH(1,0)+=w_gauss*(DN(0,0)*Nenr[1] + DNenr(1,0)*crH5);
rH(1,1)+=w_gauss*(DN(0,1)*Nenr[1] + DNenr(1,1)*crH5);
rH(1,2)+=w_gauss*(DN(0,2)*Nenr[1] + DNenr(1,2)*crH5);
rH(1,3)+=crH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
rH(1,4)+=w_gauss*(DN(1,0)*Nenr[1] + DNenr(1,0)*crH7);
rH(1,5)+=w_gauss*(DN(1,1)*Nenr[1] + DNenr(1,1)*crH7);
rH(1,6)+=w_gauss*(DN(1,2)*Nenr[1] + DNenr(1,2)*crH7);
rH(1,7)+=crH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
rH(1,8)+=w_gauss*(DN(2,0)*Nenr[1] + DNenr(1,0)*crH8);
rH(1,9)+=w_gauss*(DN(2,1)*Nenr[1] + DNenr(1,1)*crH8);
rH(1,10)+=w_gauss*(DN(2,2)*Nenr[1] + DNenr(1,2)*crH8);
rH(1,11)+=crH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
rH(1,12)+=w_gauss*(DN(3,0)*Nenr[1] + DNenr(1,0)*crH9);
rH(1,13)+=w_gauss*(DN(3,1)*Nenr[1] + DNenr(1,1)*crH9);
rH(1,14)+=w_gauss*(DN(3,2)*Nenr[1] + DNenr(1,2)*crH9);
rH(1,15)+=crH6*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
rH(2,0)+=w_gauss*(DN(0,0)*Nenr[2] + DNenr(2,0)*crH5);
rH(2,1)+=w_gauss*(DN(0,1)*Nenr[2] + DNenr(2,1)*crH5);
rH(2,2)+=w_gauss*(DN(0,2)*Nenr[2] + DNenr(2,2)*crH5);
rH(2,3)+=crH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
rH(2,4)+=w_gauss*(DN(1,0)*Nenr[2] + DNenr(2,0)*crH7);
rH(2,5)+=w_gauss*(DN(1,1)*Nenr[2] + DNenr(2,1)*crH7);
rH(2,6)+=w_gauss*(DN(1,2)*Nenr[2] + DNenr(2,2)*crH7);
rH(2,7)+=crH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
rH(2,8)+=w_gauss*(DN(2,0)*Nenr[2] + DNenr(2,0)*crH8);
rH(2,9)+=w_gauss*(DN(2,1)*Nenr[2] + DNenr(2,1)*crH8);
rH(2,10)+=w_gauss*(DN(2,2)*Nenr[2] + DNenr(2,2)*crH8);
rH(2,11)+=crH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
rH(2,12)+=w_gauss*(DN(3,0)*Nenr[2] + DNenr(2,0)*crH9);
rH(2,13)+=w_gauss*(DN(3,1)*Nenr[2] + DNenr(2,1)*crH9);
rH(2,14)+=w_gauss*(DN(3,2)*Nenr[2] + DNenr(2,2)*crH9);
rH(2,15)+=crH6*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
rH(3,0)+=w_gauss*(DN(0,0)*Nenr[3] + DNenr(3,0)*crH5);
rH(3,1)+=w_gauss*(DN(0,1)*Nenr[3] + DNenr(3,1)*crH5);
rH(3,2)+=w_gauss*(DN(0,2)*Nenr[3] + DNenr(3,2)*crH5);
rH(3,3)+=crH6*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
rH(3,4)+=w_gauss*(DN(1,0)*Nenr[3] + DNenr(3,0)*crH7);
rH(3,5)+=w_gauss*(DN(1,1)*Nenr[3] + DNenr(3,1)*crH7);
rH(3,6)+=w_gauss*(DN(1,2)*Nenr[3] + DNenr(3,2)*crH7);
rH(3,7)+=crH6*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
rH(3,8)+=w_gauss*(DN(2,0)*Nenr[3] + DNenr(3,0)*crH8);
rH(3,9)+=w_gauss*(DN(2,1)*Nenr[3] + DNenr(3,1)*crH8);
rH(3,10)+=w_gauss*(DN(2,2)*Nenr[3] + DNenr(3,2)*crH8);
rH(3,11)+=crH6*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
rH(3,12)+=w_gauss*(DN(3,0)*Nenr[3] + DNenr(3,0)*crH9);
rH(3,13)+=w_gauss*(DN(3,1)*Nenr[3] + DNenr(3,1)*crH9);
rH(3,14)+=w_gauss*(DN(3,2)*Nenr[3] + DNenr(3,2)*crH9);
rH(3,15)+=crH6*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double crKee0 = 1.0*w_gauss*1.0/(dyn_tau*rho*1.0/dt + mu*stab_c1*1.0/(h*h) + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))*1.0/h);
const double crKee1 = crKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double crKee2 = crKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double crKee3 = crKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double crKee4 = crKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double crKee5 = crKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double crKee6 = crKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
rKee(0,0)+=crKee0*(DNenr(0,0)*DNenr(0,0) + DNenr(0,1)*DNenr(0,1) + DNenr(0,2)*DNenr(0,2));
rKee(0,1)+=crKee1;
rKee(0,2)+=crKee2;
rKee(0,3)+=crKee3;
rKee(1,0)+=crKee1;
rKee(1,1)+=crKee0*(DNenr(1,0)*DNenr(1,0) + DNenr(1,1)*DNenr(1,1) + DNenr(1,2)*DNenr(1,2));
rKee(1,2)+=crKee4;
rKee(1,3)+=crKee5;
rKee(2,0)+=crKee2;
rKee(2,1)+=crKee4;
rKee(2,2)+=crKee0*(DNenr(2,0)*DNenr(2,0) + DNenr(2,1)*DNenr(2,1) + DNenr(2,2)*DNenr(2,2));
rKee(2,3)+=crKee6;
rKee(3,0)+=crKee3;
rKee(3,1)+=crKee5;
rKee(3,2)+=crKee6;
rKee(3,3)+=crKee0*(DNenr(3,0)*DNenr(3,0) + DNenr(3,1)*DNenr(3,1) + DNenr(3,2)*DNenr(3,2));


    const double crRHS_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS_ee2 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS_ee3 = crRHS_ee0 + crRHS_ee1 + crRHS_ee2 - volume_error_ratio;
const double crRHS_ee4 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crRHS_ee5 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crRHS_ee6 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crRHS_ee7 = N[0]*bdf0;
const double crRHS_ee8 = N[1]*bdf0;
const double crRHS_ee9 = N[2]*bdf0;
const double crRHS_ee10 = N[3]*bdf0;
const double crRHS_ee11 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS_ee12 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS_ee13 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS_ee14 = 1.0/dt;
const double crRHS_ee15 = N[0]*crRHS_ee14;
const double crRHS_ee16 = N[1]*crRHS_ee14;
const double crRHS_ee17 = N[2]*crRHS_ee14;
const double crRHS_ee18 = N[3]*crRHS_ee14;
const double crRHS_ee19 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crRHS_ee20 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crRHS_ee21 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crRHS_ee22 = 1.0*1.0/(crRHS_ee14*dyn_tau*rho + mu*stab_c1*1.0/(h*h) + rho*stab_c2*1.0/h*sqrt(crRHS_ee11*crRHS_ee11 + crRHS_ee12*crRHS_ee12 + crRHS_ee13*crRHS_ee13));
const double crRHS_ee23 = crRHS_ee22*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(crRHS_ee4*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crRHS_ee5*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0)) + crRHS_ee6*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(crRHS_ee0*crRHS_ee11 + crRHS_ee10*(v(3,0) - vfrac(3,0)) + crRHS_ee12*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crRHS_ee13*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crRHS_ee7*(v(0,0) - vfrac(0,0)) + crRHS_ee8*(v(1,0) - vfrac(1,0)) + crRHS_ee9*(v(2,0) - vfrac(2,0))) + rho*(crRHS_ee15*(vn(0,0) - vnn(0,0)) + crRHS_ee16*(vn(1,0) - vnn(1,0)) + crRHS_ee17*(vn(2,0) - vnn(2,0)) + crRHS_ee18*(vn(3,0) - vnn(3,0)) + crRHS_ee19*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crRHS_ee20*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crRHS_ee21*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0))));
const double crRHS_ee24 = crRHS_ee22*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(crRHS_ee4*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crRHS_ee5*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1)) + crRHS_ee6*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(crRHS_ee1*crRHS_ee12 + crRHS_ee10*(v(3,1) - vfrac(3,1)) + crRHS_ee11*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crRHS_ee13*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crRHS_ee7*(v(0,1) - vfrac(0,1)) + crRHS_ee8*(v(1,1) - vfrac(1,1)) + crRHS_ee9*(v(2,1) - vfrac(2,1))) + rho*(crRHS_ee15*(vn(0,1) - vnn(0,1)) + crRHS_ee16*(vn(1,1) - vnn(1,1)) + crRHS_ee17*(vn(2,1) - vnn(2,1)) + crRHS_ee18*(vn(3,1) - vnn(3,1)) + crRHS_ee19*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crRHS_ee20*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crRHS_ee21*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1))));
const double crRHS_ee25 = crRHS_ee22*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(crRHS_ee4*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crRHS_ee5*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2)) + crRHS_ee6*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2))) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(crRHS_ee10*(v(3,2) - vfrac(3,2)) + crRHS_ee11*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crRHS_ee12*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crRHS_ee13*crRHS_ee2 + crRHS_ee7*(v(0,2) - vfrac(0,2)) + crRHS_ee8*(v(1,2) - vfrac(1,2)) + crRHS_ee9*(v(2,2) - vfrac(2,2))) + rho*(crRHS_ee15*(vn(0,2) - vnn(0,2)) + crRHS_ee16*(vn(1,2) - vnn(1,2)) + crRHS_ee17*(vn(2,2) - vnn(2,2)) + crRHS_ee18*(vn(3,2) - vnn(3,2)) + crRHS_ee19*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crRHS_ee20*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crRHS_ee21*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2))));
rRHS_ee[0]+=-w_gauss*(DNenr(0,0)*crRHS_ee23 + DNenr(0,1)*crRHS_ee24 + DNenr(0,2)*crRHS_ee25 + Nenr[0]*crRHS_ee3);
rRHS_ee[1]+=-w_gauss*(DNenr(1,0)*crRHS_ee23 + DNenr(1,1)*crRHS_ee24 + DNenr(1,2)*crRHS_ee25 + Nenr[1]*crRHS_ee3);
rRHS_ee[2]+=-w_gauss*(DNenr(2,0)*crRHS_ee23 + DNenr(2,1)*crRHS_ee24 + DNenr(2,2)*crRHS_ee25 + Nenr[2]*crRHS_ee3);
rRHS_ee[3]+=-w_gauss*(DNenr(3,0)*crRHS_ee23 + DNenr(3,1)*crRHS_ee24 + DNenr(3,2)*crRHS_ee25 + Nenr[3]*crRHS_ee3);

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
template <>
double TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>::CalculateArtificialDynamicViscositySpecialization(
       TwoFluidNavierStokesFractionalData<2, 3> &rData) const
{
    // Variables for artificial viscosity calculation
    double artificial_mu = 0.0;
    const double rho = rData.Density;
    const double bdf0 = rData.bdf0;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const BoundedMatrix<double, 3, 2> vconv = vn - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const double art_dyn_visc_coeff = 0.8;
    double grad_v_norm = 0.0;
    
    grad_v_norm=sqrt(pow(fabs(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)), 2) + pow(fabs(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)), 2) + pow(fabs(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)), 2) + pow(fabs(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)), 2));

    // Check that velocity gradient norm is non-zero
    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        artificial_mu=0.5*art_dyn_visc_coeff*h*sqrt(pow(fabs(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + rho*((DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0)) + (DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1))) - rho*((DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0))*(N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0)) + (DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0))*(N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) - bdf0*vfrac(0,0)) + N[1]*(bdf0*v(1,0) - bdf0*vfrac(1,0)) + N[2]*(bdf0*v(2,0) - bdf0*vfrac(2,0))) + rho*(N[0]*(vn(0,0) - vnn(0,0))/dt + N[1]*(vn(1,0) - vnn(1,0))/dt + N[2]*(vn(2,0) - vnn(2,0))/dt + (DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0))*(N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0)) + (DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0))*(N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1)))), 2) + pow(fabs(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + rho*((DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0)) + (DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1))) - rho*((DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1))*(N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0)) + (DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1))*(N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) - bdf0*vfrac(0,1)) + N[1]*(bdf0*v(1,1) - bdf0*vfrac(1,1)) + N[2]*(bdf0*v(2,1) - bdf0*vfrac(2,1))) + rho*(N[0]*(vn(0,1) - vnn(0,1))/dt + N[1]*(vn(1,1) - vnn(1,1))/dt + N[2]*(vn(2,1) - vnn(2,1))/dt + (DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1))*(N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0)) + (DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1))*(N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1)))), 2))/sqrt(pow(fabs(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)), 2) + pow(fabs(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)), 2) + pow(fabs(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)), 2) + pow(fabs(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)), 2));

    }

    return artificial_mu;
}

template <>
double TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>::CalculateArtificialDynamicViscositySpecialization(
       TwoFluidNavierStokesFractionalData<3, 4> &rData) const
{
    // Variables for artificial viscosity calculation
    double artificial_mu = 0.0;
    const double rho = rData.Density;
    const double bdf0 = rData.bdf0;
    const double h = rData.ElementSize;
    const double dt = rData.DeltaTime;
    const auto &v = rData.Velocity;
    const auto &vn = rData.VelocityOldStep1;
    const auto &vnn = rData.VelocityOldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &f = rData.BodyForce;
    const auto &p = rData.Pressure;
    const BoundedMatrix<double, 4, 3>vconv = vn - vmesh;
    const auto vfrac = rData.FractionalVelocity;
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    const double art_dyn_visc_coeff = 0.8;
    double grad_v_norm = 0.0;

    grad_v_norm=sqrt(pow(fabs(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)), 2) + pow(fabs(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)), 2) + pow(fabs(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)), 2) + pow(fabs(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)), 2) + pow(fabs(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)), 2) + pow(fabs(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)), 2) + pow(fabs(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)), 2) + pow(fabs(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)), 2) + pow(fabs(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)), 2));

    // Check that velocity gradient norm is non-zero
    if (grad_v_norm > 1.0e-12) {
        // Calculate symbolic artificial viscosity
        artificial_mu=0.5*art_dyn_visc_coeff*h*sqrt(pow(fabs(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + rho*((DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*((DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0))*(N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0)) + (DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0))*(N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1)) + (DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0))*(N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2))) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) - bdf0*vfrac(0,0)) + N[1]*(bdf0*v(1,0) - bdf0*vfrac(1,0)) + N[2]*(bdf0*v(2,0) - bdf0*vfrac(2,0)) + N[3]*(bdf0*v(3,0) - bdf0*vfrac(3,0))) + rho*(N[0]*(vn(0,0) - vnn(0,0))/dt + N[1]*(vn(1,0) - vnn(1,0))/dt + N[2]*(vn(2,0) - vnn(2,0))/dt + N[3]*(vn(3,0) - vnn(3,0))/dt + (DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0))*(N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0)) + (DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0))*(N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1)) + (DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0))*(N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2)))), 2) + pow(fabs(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + rho*((DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*((DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1))*(N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0)) + (DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1))*(N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1)) + (DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1))*(N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2))) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) - bdf0*vfrac(0,1)) + N[1]*(bdf0*v(1,1) - bdf0*vfrac(1,1)) + N[2]*(bdf0*v(2,1) - bdf0*vfrac(2,1)) + N[3]*(bdf0*v(3,1) - bdf0*vfrac(3,1))) + rho*(N[0]*(vn(0,1) - vnn(0,1))/dt + N[1]*(vn(1,1) - vnn(1,1))/dt + N[2]*(vn(2,1) - vnn(2,1))/dt + N[3]*(vn(3,1) - vnn(3,1))/dt + (DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1))*(N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0)) + (DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1))*(N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1)) + (DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1))*(N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2)))), 2) + pow(fabs(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + rho*((DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2))*(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0)) + (DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))*(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1)) + (DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2))*(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2))) - rho*((DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2))*(N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0)) + (DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2))*(N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1)) + (DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2))*(N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2))) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) - bdf0*vfrac(0,2)) + N[1]*(bdf0*v(1,2) - bdf0*vfrac(1,2)) + N[2]*(bdf0*v(2,2) - bdf0*vfrac(2,2)) + N[3]*(bdf0*v(3,2) - bdf0*vfrac(3,2))) + rho*(N[0]*(vn(0,2) - vnn(0,2))/dt + N[1]*(vn(1,2) - vnn(1,2))/dt + N[2]*(vn(2,2) - vnn(2,2))/dt + N[3]*(vn(3,2) - vnn(3,2))/dt + (DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2))*(N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0)) + (DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2))*(N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1)) + (DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2))*(N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2)))), 2))/sqrt(pow(fabs(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)), 2) + pow(fabs(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)), 2) + pow(fabs(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)), 2) + pow(fabs(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)), 2) + pow(fabs(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)), 2) + pow(fabs(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)), 2) + pow(fabs(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)), 2) + pow(fabs(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)), 2) + pow(fabs(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)), 2));

    }

    return artificial_mu;
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
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Create new temporary data container
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const unsigned int number_of_gauss_points = gauss_weights.size();

    if (rOutput.size() != number_of_gauss_points){
        rOutput.resize(number_of_gauss_points);
    }

    if (rVariable == ARTIFICIAL_DYNAMIC_VISCOSITY){
        // Iterate over integration points to evaluate the artificial viscosity at each Gauss point
        for (unsigned int g = 0; g < number_of_gauss_points; ++g){
            this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            rOutput[g] = this->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
        }
    }
    else{
        FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TElementData>
void TwoFluidNavierStokesFractional<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Create new temporary data container
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const unsigned int number_of_gauss_points = gauss_weights.size();
    rOutput = 0.0;
    if (rVariable == ARTIFICIAL_DYNAMIC_VISCOSITY)
    {

        // Iterate over integration points to evaluate the artificial viscosity at each Gauss point
        for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        {
            this->UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            rOutput += CalculateArtificialDynamicViscositySpecialization(data);
        }

        rOutput /= number_of_gauss_points;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<2, 3>>;
template class TwoFluidNavierStokesFractional<TwoFluidNavierStokesFractionalData<3, 4>>;

} // namespace Kratos