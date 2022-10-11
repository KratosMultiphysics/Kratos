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
#include "custom_utilities/two_fluid_navier_stokes_alpha_method_data.h"

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

                Matrix int_shape_function, int_shape_function_enr_neg, int_shape_function_enr_pos;
                GeometryType::ShapeFunctionsGradientsType int_shape_derivatives;
                Vector int_gauss_pts_weights;
                std::vector< array_1d<double,3> > int_normals_neg;

                if (rCurrentProcessInfo[SURFACE_TENSION] || rCurrentProcessInfo[MOMENTUM_CORRECTION]){
                    ComputeSplitInterface(
                        data,
                        int_shape_function,
                        int_shape_function_enr_pos,
                        int_shape_function_enr_neg,
                        int_shape_derivatives,
                        int_gauss_pts_weights,
                        int_normals_neg,
                        p_modified_sh_func);
                }

                if (rCurrentProcessInfo[MOMENTUM_CORRECTION]){
                    BoundedMatrix<double, LocalSize, LocalSize> lhs_acc_correction = ZeroMatrix(LocalSize,LocalSize);

                    double positive_density = 0.0;
                    double negative_density = 0.0;

                    const auto& r_geom = this->GetGeometry();

                    for (unsigned int intgp = 0; intgp < int_gauss_pts_weights.size(); ++intgp){
                        double u_dot_n = 0.0;
                        for (unsigned int i = 0; i < NumNodes; ++i){
                            u_dot_n += int_shape_function(intgp,i)*r_geom[i].GetValue(DISTANCE_CORRECTION);

                            if (data.Distance[i] > 0.0){
                                positive_density = data.NodalDensity[i];
                            } else {
                                negative_density = data.NodalDensity[i];
                            }
                        }

                        u_dot_n /= data.DeltaTime;

                        for (unsigned int i = 0; i < NumNodes; ++i){
                            for (unsigned int j = 0; j < NumNodes; ++j){
                                for (unsigned int dim = 0; dim < NumNodes-1; ++dim){
                                    lhs_acc_correction( i*(NumNodes) + dim, j*(NumNodes) + dim) +=
                                        int_shape_function(intgp,i)*int_shape_function(intgp,j)*u_dot_n*int_gauss_pts_weights(intgp);
                                }
                            }
                        }
                    }

                    lhs_acc_correction = (negative_density - positive_density)*lhs_acc_correction;
                    noalias(rLeftHandSideMatrix) += lhs_acc_correction;

                    Kratos::array_1d<double, LocalSize> tempU; // Unknowns vector containing only velocity components
                    for (unsigned int i = 0; i < NumNodes; ++i){
                        for (unsigned int dimi = 0; dimi < Dim; ++dimi){
                            tempU[i*(Dim+1) + dimi] = data.Velocity(i,dimi);
                        }
                    }
                    noalias(rRightHandSideVector) -= prod(lhs_acc_correction,tempU);
                }

                if (rCurrentProcessInfo[SURFACE_TENSION]){

                    AddSurfaceTensionContribution(
                        data,
                        int_shape_function,
                        int_shape_function_enr_pos,
                        int_shape_function_enr_neg,
                        int_shape_derivatives,
                        int_gauss_pts_weights,
                        int_normals_neg,
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
const Parameters TwoFluidNavierStokes<TElementData>::GetSpecifications() const
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
        "required_variables"         : ["DISTANCE","VELOCITY","PRESSURE","MESH_VELOCITY","DENSITY","DYNAMIC_VISCOSITY"],
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
            "This element implements Navier-Stokes biphasic fluid-air formulation with a levelset-based interface representation with Variational MultiScales (VMS) stabilization. Note that any viscous behavior can be used for the fluid phase through a constitutive law. The air phase is assumed to be Newtonian. Surface tension contribution can be accounted for by setting the SURFACE_TENSION variable to true in the ProcessInfo container.
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

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = pow(DN(0,0), 2);
const double clhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 = rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 = rho*(clhs6*h/stab_c1 + mu);
const double clhs8 = pow(N[0], 2);
const double clhs9 = rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs10 = bdf0*rho;
const double clhs11 = K_darcy*N[0];
const double clhs12 = N[0]*clhs10;
const double clhs13 = clhs11 + clhs12 + clhs9;
const double clhs14 = 1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 = 1.0*clhs9;
const double clhs16 = clhs14*clhs15;
const double clhs17 = 1.0*clhs11;
const double clhs18 = clhs14*clhs17;
const double clhs19 = 1.0*clhs14;
const double clhs20 = clhs19*rho;
const double clhs21 = clhs20*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double clhs22 = N[0]*clhs21;
const double clhs23 = K_darcy*clhs8 + N[0]*clhs9 + clhs10*clhs8 + clhs13*clhs16 - clhs13*clhs18 + clhs13*clhs22;
const double clhs24 = C(0,1)*DN(0,1) + clhs1;
const double clhs25 = C(1,2)*DN(0,1);
const double clhs26 = C(2,2)*DN(0,0) + clhs25;
const double clhs27 = DN(0,0)*clhs7;
const double clhs28 = DN(0,1)*clhs27;
const double clhs29 = -N[0] + clhs16 - clhs18 + clhs22;
const double clhs30 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs31 = C(0,2)*DN(1,0);
const double clhs32 = C(2,2)*DN(1,1) + clhs31;
const double clhs33 = DN(0,0)*DN(1,0);
const double clhs34 = N[1]*clhs11 + N[1]*clhs12;
const double clhs35 = clhs33*clhs7 + clhs34;
const double clhs36 = rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs37 = K_darcy*N[1];
const double clhs38 = N[1]*clhs10;
const double clhs39 = clhs36 + clhs37 + clhs38;
const double clhs40 = N[0]*clhs36 + clhs16*clhs39 - clhs18*clhs39 + clhs22*clhs39;
const double clhs41 = C(0,1)*DN(1,1) + clhs31;
const double clhs42 = C(1,2)*DN(1,1);
const double clhs43 = C(2,2)*DN(1,0) + clhs42;
const double clhs44 = DN(1,1)*clhs27;
const double clhs45 = DN(0,0)*N[1];
const double clhs46 = DN(1,0)*N[0];
const double clhs47 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs48 = C(0,2)*DN(2,0);
const double clhs49 = C(2,2)*DN(2,1) + clhs48;
const double clhs50 = DN(0,0)*DN(2,0);
const double clhs51 = N[2]*clhs11 + N[2]*clhs12;
const double clhs52 = clhs50*clhs7 + clhs51;
const double clhs53 = rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs54 = K_darcy*N[2];
const double clhs55 = N[2]*clhs10;
const double clhs56 = clhs53 + clhs54 + clhs55;
const double clhs57 = N[0]*clhs53 + clhs16*clhs56 - clhs18*clhs56 + clhs22*clhs56;
const double clhs58 = C(0,1)*DN(2,1) + clhs48;
const double clhs59 = C(1,2)*DN(2,1);
const double clhs60 = C(2,2)*DN(2,0) + clhs59;
const double clhs61 = DN(2,1)*clhs27;
const double clhs62 = DN(0,0)*N[2];
const double clhs63 = DN(2,0)*N[0];
const double clhs64 = C(0,1)*DN(0,0) + clhs25;
const double clhs65 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs66 = pow(DN(0,1), 2);
const double clhs67 = C(0,1)*DN(1,0) + clhs42;
const double clhs68 = DN(0,1)*clhs7;
const double clhs69 = DN(1,0)*clhs68;
const double clhs70 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs71 = DN(0,1)*DN(1,1);
const double clhs72 = clhs34 + clhs7*clhs71;
const double clhs73 = DN(0,1)*N[1];
const double clhs74 = DN(1,1)*N[0];
const double clhs75 = C(0,1)*DN(2,0) + clhs59;
const double clhs76 = DN(2,0)*clhs68;
const double clhs77 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs78 = DN(0,1)*DN(2,1);
const double clhs79 = clhs51 + clhs7*clhs78;
const double clhs80 = DN(0,1)*N[2];
const double clhs81 = DN(2,1)*N[0];
const double clhs82 = rho*(N[0] + clhs14*(1.0*clhs12 + clhs15 + clhs17));
const double clhs83 = clhs19*clhs39;
const double clhs84 = clhs20*(clhs33 + clhs71);
const double clhs85 = clhs19*clhs56;
const double clhs86 = clhs20*(clhs50 + clhs78);
const double clhs87 = 1.0*clhs36;
const double clhs88 = clhs14*clhs87;
const double clhs89 = 1.0*clhs37;
const double clhs90 = clhs14*clhs89;
const double clhs91 = N[1]*clhs21;
const double clhs92 = N[1]*clhs9 + clhs13*clhs88 - clhs13*clhs90 + clhs13*clhs91;
const double clhs93 = pow(DN(1,0), 2);
const double clhs94 = pow(N[1], 2);
const double clhs95 = K_darcy*clhs94 + N[1]*clhs36 + clhs10*clhs94 + clhs39*clhs88 - clhs39*clhs90 + clhs39*clhs91;
const double clhs96 = DN(1,0)*clhs7;
const double clhs97 = DN(1,1)*clhs96;
const double clhs98 = -N[1] + clhs88 - clhs90 + clhs91;
const double clhs99 = DN(1,0)*DN(2,0);
const double clhs100 = N[2]*clhs37 + N[2]*clhs38;
const double clhs101 = clhs100 + clhs7*clhs99;
const double clhs102 = N[1]*clhs53 + clhs56*clhs88 - clhs56*clhs90 + clhs56*clhs91;
const double clhs103 = DN(2,1)*clhs96;
const double clhs104 = DN(1,0)*N[2];
const double clhs105 = DN(2,0)*N[1];
const double clhs106 = pow(DN(1,1), 2);
const double clhs107 = DN(2,0)*clhs7;
const double clhs108 = DN(1,1)*clhs107;
const double clhs109 = DN(1,1)*DN(2,1);
const double clhs110 = clhs100 + clhs109*clhs7;
const double clhs111 = DN(1,1)*N[2];
const double clhs112 = DN(2,1)*N[1];
const double clhs113 = clhs13*clhs19;
const double clhs114 = rho*(N[1] + clhs14*(1.0*clhs38 + clhs87 + clhs89));
const double clhs115 = clhs20*(clhs109 + clhs99);
const double clhs116 = 1.0*clhs53;
const double clhs117 = clhs116*clhs14;
const double clhs118 = 1.0*clhs54;
const double clhs119 = clhs118*clhs14;
const double clhs120 = N[2]*clhs21;
const double clhs121 = N[2]*clhs9 + clhs117*clhs13 - clhs119*clhs13 + clhs120*clhs13;
const double clhs122 = N[2]*clhs36 + clhs117*clhs39 - clhs119*clhs39 + clhs120*clhs39;
const double clhs123 = pow(DN(2,0), 2);
const double clhs124 = pow(N[2], 2);
const double clhs125 = K_darcy*clhs124 + N[2]*clhs53 + clhs10*clhs124 + clhs117*clhs56 - clhs119*clhs56 + clhs120*clhs56;
const double clhs126 = DN(2,1)*clhs107;
const double clhs127 = -N[2] + clhs117 - clhs119 + clhs120;
const double clhs128 = pow(DN(2,1), 2);
const double clhs129 = rho*(N[2] + clhs14*(clhs116 + clhs118 + 1.0*clhs55));
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs23 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs28;
lhs(0,2)=DN(0,0)*clhs29;
lhs(0,3)=DN(0,0)*clhs30 + DN(0,1)*clhs32 + clhs35 + clhs40;
lhs(0,4)=DN(0,0)*clhs41 + DN(0,1)*clhs43 + clhs44;
lhs(0,5)=DN(1,0)*clhs16 - DN(1,0)*clhs18 + clhs21*clhs46 - clhs45;
lhs(0,6)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + clhs52 + clhs57;
lhs(0,7)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + clhs61;
lhs(0,8)=DN(2,0)*clhs16 - DN(2,0)*clhs18 + clhs21*clhs63 - clhs62;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs64 + clhs28;
lhs(1,1)=DN(0,0)*clhs26 + DN(0,1)*clhs65 + clhs23 + clhs66*clhs7;
lhs(1,2)=DN(0,1)*clhs29;
lhs(1,3)=DN(0,0)*clhs32 + DN(0,1)*clhs67 + clhs69;
lhs(1,4)=DN(0,0)*clhs43 + DN(0,1)*clhs70 + clhs40 + clhs72;
lhs(1,5)=DN(1,1)*clhs16 - DN(1,1)*clhs18 + clhs21*clhs74 - clhs73;
lhs(1,6)=DN(0,0)*clhs49 + DN(0,1)*clhs75 + clhs76;
lhs(1,7)=DN(0,0)*clhs60 + DN(0,1)*clhs77 + clhs57 + clhs79;
lhs(1,8)=DN(2,1)*clhs16 - DN(2,1)*clhs18 + clhs21*clhs81 - clhs80;
lhs(2,0)=DN(0,0)*clhs82;
lhs(2,1)=DN(0,1)*clhs82;
lhs(2,2)=clhs20*(clhs3 + clhs66);
lhs(2,3)=rho*(DN(0,0)*clhs83 + clhs46);
lhs(2,4)=rho*(DN(0,1)*clhs83 + clhs74);
lhs(2,5)=clhs84;
lhs(2,6)=rho*(DN(0,0)*clhs85 + clhs63);
lhs(2,7)=rho*(DN(0,1)*clhs85 + clhs81);
lhs(2,8)=clhs86;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs35 + clhs92;
lhs(3,1)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs69;
lhs(3,2)=DN(0,0)*clhs88 - DN(0,0)*clhs90 + clhs21*clhs45 - clhs46;
lhs(3,3)=DN(1,0)*clhs30 + DN(1,1)*clhs32 + clhs7*clhs93 + clhs95;
lhs(3,4)=DN(1,0)*clhs41 + DN(1,1)*clhs43 + clhs97;
lhs(3,5)=DN(1,0)*clhs98;
lhs(3,6)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + clhs101 + clhs102;
lhs(3,7)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + clhs103;
lhs(3,8)=DN(2,0)*clhs88 - DN(2,0)*clhs90 - clhs104 + clhs105*clhs21;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs64 + clhs44;
lhs(4,1)=DN(1,0)*clhs26 + DN(1,1)*clhs65 + clhs72 + clhs92;
lhs(4,2)=DN(0,1)*clhs88 - DN(0,1)*clhs90 + clhs21*clhs73 - clhs74;
lhs(4,3)=DN(1,0)*clhs32 + DN(1,1)*clhs67 + clhs97;
lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs70 + clhs106*clhs7 + clhs95;
lhs(4,5)=DN(1,1)*clhs98;
lhs(4,6)=DN(1,0)*clhs49 + DN(1,1)*clhs75 + clhs108;
lhs(4,7)=DN(1,0)*clhs60 + DN(1,1)*clhs77 + clhs102 + clhs110;
lhs(4,8)=DN(2,1)*clhs88 - DN(2,1)*clhs90 - clhs111 + clhs112*clhs21;
lhs(5,0)=rho*(DN(1,0)*clhs113 + clhs45);
lhs(5,1)=rho*(DN(1,1)*clhs113 + clhs73);
lhs(5,2)=clhs84;
lhs(5,3)=DN(1,0)*clhs114;
lhs(5,4)=DN(1,1)*clhs114;
lhs(5,5)=clhs20*(clhs106 + clhs93);
lhs(5,6)=rho*(DN(1,0)*clhs85 + clhs105);
lhs(5,7)=rho*(DN(1,1)*clhs85 + clhs112);
lhs(5,8)=clhs115;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs121 + clhs52;
lhs(6,1)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs76;
lhs(6,2)=DN(0,0)*clhs117 - DN(0,0)*clhs119 + clhs21*clhs62 - clhs63;
lhs(6,3)=DN(2,0)*clhs30 + DN(2,1)*clhs32 + clhs101 + clhs122;
lhs(6,4)=DN(2,0)*clhs41 + DN(2,1)*clhs43 + clhs108;
lhs(6,5)=DN(1,0)*clhs117 - DN(1,0)*clhs119 + clhs104*clhs21 - clhs105;
lhs(6,6)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + clhs123*clhs7 + clhs125;
lhs(6,7)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + clhs126;
lhs(6,8)=DN(2,0)*clhs127;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs64 + clhs61;
lhs(7,1)=DN(2,0)*clhs26 + DN(2,1)*clhs65 + clhs121 + clhs79;
lhs(7,2)=DN(0,1)*clhs117 - DN(0,1)*clhs119 + clhs21*clhs80 - clhs81;
lhs(7,3)=DN(2,0)*clhs32 + DN(2,1)*clhs67 + clhs103;
lhs(7,4)=DN(2,0)*clhs43 + DN(2,1)*clhs70 + clhs110 + clhs122;
lhs(7,5)=DN(1,1)*clhs117 - DN(1,1)*clhs119 + clhs111*clhs21 - clhs112;
lhs(7,6)=DN(2,0)*clhs49 + DN(2,1)*clhs75 + clhs126;
lhs(7,7)=DN(2,0)*clhs60 + DN(2,1)*clhs77 + clhs125 + clhs128*clhs7;
lhs(7,8)=DN(2,1)*clhs127;
lhs(8,0)=rho*(DN(2,0)*clhs113 + clhs62);
lhs(8,1)=rho*(DN(2,1)*clhs113 + clhs80);
lhs(8,2)=clhs86;
lhs(8,3)=rho*(DN(2,0)*clhs83 + clhs104);
lhs(8,4)=rho*(DN(2,1)*clhs83 + clhs111);
lhs(8,5)=clhs115;
lhs(8,6)=DN(2,0)*clhs129;
lhs(8,7)=DN(2,1)*clhs129;
lhs(8,8)=clhs20*(clhs123 + clhs128);


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

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 = rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 = rho*(clhs9*h/stab_c1 + mu);
const double clhs11 = pow(N[0], 2);
const double clhs12 = rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs13 = bdf0*rho;
const double clhs14 = K_darcy*N[0];
const double clhs15 = N[0]*clhs13;
const double clhs16 = clhs12 + clhs14 + clhs15;
const double clhs17 = 1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = 1.0*clhs12;
const double clhs19 = clhs17*clhs18;
const double clhs20 = 1.0*clhs14;
const double clhs21 = clhs17*clhs20;
const double clhs22 = 1.0*clhs17;
const double clhs23 = clhs22*rho;
const double clhs24 = clhs23*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double clhs25 = N[0]*clhs24;
const double clhs26 = K_darcy*clhs11 + N[0]*clhs12 + clhs11*clhs13 + clhs16*clhs19 - clhs16*clhs21 + clhs16*clhs25;
const double clhs27 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs28 = C(1,3)*DN(0,1);
const double clhs29 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs28;
const double clhs30 = C(3,5)*DN(0,0);
const double clhs31 = C(4,5)*DN(0,2);
const double clhs32 = C(1,5)*DN(0,1) + clhs30 + clhs31;
const double clhs33 = DN(0,0)*clhs10;
const double clhs34 = DN(0,1)*clhs33;
const double clhs35 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs36 = C(3,4)*DN(0,1);
const double clhs37 = C(2,3)*DN(0,2) + clhs30 + clhs36;
const double clhs38 = C(2,5)*DN(0,2);
const double clhs39 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs38;
const double clhs40 = DN(0,2)*clhs33;
const double clhs41 = -N[0] + clhs19 - clhs21 + clhs25;
const double clhs42 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs43 = C(0,3)*DN(1,0);
const double clhs44 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs43;
const double clhs45 = C(0,5)*DN(1,0);
const double clhs46 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs45;
const double clhs47 = DN(0,0)*DN(1,0);
const double clhs48 = N[1]*clhs14 + N[1]*clhs15;
const double clhs49 = clhs10*clhs47 + clhs48;
const double clhs50 = rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs51 = K_darcy*N[1];
const double clhs52 = N[1]*clhs13;
const double clhs53 = clhs50 + clhs51 + clhs52;
const double clhs54 = N[0]*clhs50 + clhs19*clhs53 - clhs21*clhs53 + clhs25*clhs53;
const double clhs55 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs43;
const double clhs56 = C(1,3)*DN(1,1);
const double clhs57 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs56;
const double clhs58 = C(3,5)*DN(1,0);
const double clhs59 = C(4,5)*DN(1,2);
const double clhs60 = C(1,5)*DN(1,1) + clhs58 + clhs59;
const double clhs61 = DN(1,1)*clhs33;
const double clhs62 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs45;
const double clhs63 = C(3,4)*DN(1,1);
const double clhs64 = C(2,3)*DN(1,2) + clhs58 + clhs63;
const double clhs65 = C(2,5)*DN(1,2);
const double clhs66 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs65;
const double clhs67 = DN(1,2)*clhs33;
const double clhs68 = DN(0,0)*N[1];
const double clhs69 = DN(1,0)*N[0];
const double clhs70 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs71 = C(0,3)*DN(2,0);
const double clhs72 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs71;
const double clhs73 = C(0,5)*DN(2,0);
const double clhs74 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs73;
const double clhs75 = DN(0,0)*DN(2,0);
const double clhs76 = N[2]*clhs14 + N[2]*clhs15;
const double clhs77 = clhs10*clhs75 + clhs76;
const double clhs78 = rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs79 = K_darcy*N[2];
const double clhs80 = N[2]*clhs13;
const double clhs81 = clhs78 + clhs79 + clhs80;
const double clhs82 = N[0]*clhs78 + clhs19*clhs81 - clhs21*clhs81 + clhs25*clhs81;
const double clhs83 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs71;
const double clhs84 = C(1,3)*DN(2,1);
const double clhs85 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs84;
const double clhs86 = C(3,5)*DN(2,0);
const double clhs87 = C(4,5)*DN(2,2);
const double clhs88 = C(1,5)*DN(2,1) + clhs86 + clhs87;
const double clhs89 = DN(2,1)*clhs33;
const double clhs90 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs73;
const double clhs91 = C(3,4)*DN(2,1);
const double clhs92 = C(2,3)*DN(2,2) + clhs86 + clhs91;
const double clhs93 = C(2,5)*DN(2,2);
const double clhs94 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs93;
const double clhs95 = DN(2,2)*clhs33;
const double clhs96 = DN(0,0)*N[2];
const double clhs97 = DN(2,0)*N[0];
const double clhs98 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs99 = C(0,3)*DN(3,0);
const double clhs100 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs99;
const double clhs101 = C(0,5)*DN(3,0);
const double clhs102 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs101;
const double clhs103 = DN(0,0)*DN(3,0);
const double clhs104 = N[3]*clhs14 + N[3]*clhs15;
const double clhs105 = clhs10*clhs103 + clhs104;
const double clhs106 = rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs107 = K_darcy*N[3];
const double clhs108 = N[3]*clhs13;
const double clhs109 = clhs106 + clhs107 + clhs108;
const double clhs110 = N[0]*clhs106 + clhs109*clhs19 - clhs109*clhs21 + clhs109*clhs25;
const double clhs111 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs99;
const double clhs112 = C(1,3)*DN(3,1);
const double clhs113 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs112;
const double clhs114 = C(3,5)*DN(3,0);
const double clhs115 = C(4,5)*DN(3,2);
const double clhs116 = C(1,5)*DN(3,1) + clhs114 + clhs115;
const double clhs117 = DN(3,1)*clhs33;
const double clhs118 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs101;
const double clhs119 = C(3,4)*DN(3,1);
const double clhs120 = C(2,3)*DN(3,2) + clhs114 + clhs119;
const double clhs121 = C(2,5)*DN(3,2);
const double clhs122 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs121;
const double clhs123 = DN(3,2)*clhs33;
const double clhs124 = DN(0,0)*N[3];
const double clhs125 = DN(3,0)*N[0];
const double clhs126 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs28;
const double clhs127 = C(0,4)*DN(0,0) + clhs31 + clhs36;
const double clhs128 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs129 = C(1,4)*DN(0,1);
const double clhs130 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs129;
const double clhs131 = pow(DN(0,1), 2);
const double clhs132 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs129;
const double clhs133 = C(2,4)*DN(0,2);
const double clhs134 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs133;
const double clhs135 = DN(0,1)*clhs10;
const double clhs136 = DN(0,2)*clhs135;
const double clhs137 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs56;
const double clhs138 = C(0,4)*DN(1,0) + clhs59 + clhs63;
const double clhs139 = DN(1,0)*clhs135;
const double clhs140 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs141 = C(1,4)*DN(1,1);
const double clhs142 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs141;
const double clhs143 = DN(0,1)*DN(1,1);
const double clhs144 = clhs10*clhs143;
const double clhs145 = clhs48 + clhs54;
const double clhs146 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs141;
const double clhs147 = C(2,4)*DN(1,2);
const double clhs148 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs147;
const double clhs149 = DN(1,2)*clhs135;
const double clhs150 = DN(0,1)*N[1];
const double clhs151 = DN(1,1)*N[0];
const double clhs152 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs84;
const double clhs153 = C(0,4)*DN(2,0) + clhs87 + clhs91;
const double clhs154 = DN(2,0)*clhs135;
const double clhs155 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs156 = C(1,4)*DN(2,1);
const double clhs157 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs156;
const double clhs158 = DN(0,1)*DN(2,1);
const double clhs159 = clhs10*clhs158;
const double clhs160 = clhs76 + clhs82;
const double clhs161 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs156;
const double clhs162 = C(2,4)*DN(2,2);
const double clhs163 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs162;
const double clhs164 = DN(2,2)*clhs135;
const double clhs165 = DN(0,1)*N[2];
const double clhs166 = DN(2,1)*N[0];
const double clhs167 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs112;
const double clhs168 = C(0,4)*DN(3,0) + clhs115 + clhs119;
const double clhs169 = DN(3,0)*clhs135;
const double clhs170 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs171 = C(1,4)*DN(3,1);
const double clhs172 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs171;
const double clhs173 = DN(0,1)*DN(3,1);
const double clhs174 = clhs10*clhs173;
const double clhs175 = clhs104 + clhs110;
const double clhs176 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs171;
const double clhs177 = C(2,4)*DN(3,2);
const double clhs178 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs177;
const double clhs179 = DN(3,2)*clhs135;
const double clhs180 = DN(0,1)*N[3];
const double clhs181 = DN(3,1)*N[0];
const double clhs182 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs38;
const double clhs183 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs133;
const double clhs184 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs185 = pow(DN(0,2), 2);
const double clhs186 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs65;
const double clhs187 = DN(0,2)*clhs10;
const double clhs188 = DN(1,0)*clhs187;
const double clhs189 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs147;
const double clhs190 = DN(1,1)*clhs187;
const double clhs191 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs192 = DN(0,2)*DN(1,2);
const double clhs193 = clhs10*clhs192;
const double clhs194 = DN(0,2)*N[1];
const double clhs195 = DN(1,2)*N[0];
const double clhs196 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs93;
const double clhs197 = DN(2,0)*clhs187;
const double clhs198 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs162;
const double clhs199 = DN(2,1)*clhs187;
const double clhs200 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs201 = DN(0,2)*DN(2,2);
const double clhs202 = clhs10*clhs201;
const double clhs203 = DN(0,2)*N[2];
const double clhs204 = DN(2,2)*N[0];
const double clhs205 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs121;
const double clhs206 = DN(3,0)*clhs187;
const double clhs207 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs177;
const double clhs208 = DN(3,1)*clhs187;
const double clhs209 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs210 = DN(0,2)*DN(3,2);
const double clhs211 = clhs10*clhs210;
const double clhs212 = DN(0,2)*N[3];
const double clhs213 = DN(3,2)*N[0];
const double clhs214 = rho*(N[0] + clhs17*(1.0*clhs15 + clhs18 + clhs20));
const double clhs215 = clhs22*clhs53;
const double clhs216 = clhs23*(clhs143 + clhs192 + clhs47);
const double clhs217 = clhs22*clhs81;
const double clhs218 = clhs23*(clhs158 + clhs201 + clhs75);
const double clhs219 = clhs109*clhs22;
const double clhs220 = clhs23*(clhs103 + clhs173 + clhs210);
const double clhs221 = 1.0*clhs50;
const double clhs222 = clhs17*clhs221;
const double clhs223 = 1.0*clhs51;
const double clhs224 = clhs17*clhs223;
const double clhs225 = N[1]*clhs24;
const double clhs226 = N[1]*clhs12 + clhs16*clhs222 - clhs16*clhs224 + clhs16*clhs225;
const double clhs227 = pow(DN(1,0), 2);
const double clhs228 = pow(N[1], 2);
const double clhs229 = K_darcy*clhs228 + N[1]*clhs50 + clhs13*clhs228 + clhs222*clhs53 - clhs224*clhs53 + clhs225*clhs53;
const double clhs230 = DN(1,0)*clhs10;
const double clhs231 = DN(1,1)*clhs230;
const double clhs232 = DN(1,2)*clhs230;
const double clhs233 = -N[1] + clhs222 - clhs224 + clhs225;
const double clhs234 = DN(1,0)*DN(2,0);
const double clhs235 = N[2]*clhs51 + N[2]*clhs52;
const double clhs236 = clhs10*clhs234 + clhs235;
const double clhs237 = N[1]*clhs78 + clhs222*clhs81 - clhs224*clhs81 + clhs225*clhs81;
const double clhs238 = DN(2,1)*clhs230;
const double clhs239 = DN(2,2)*clhs230;
const double clhs240 = DN(1,0)*N[2];
const double clhs241 = DN(2,0)*N[1];
const double clhs242 = DN(1,0)*DN(3,0);
const double clhs243 = N[3]*clhs51 + N[3]*clhs52;
const double clhs244 = clhs10*clhs242 + clhs243;
const double clhs245 = N[1]*clhs106 + clhs109*clhs222 - clhs109*clhs224 + clhs109*clhs225;
const double clhs246 = DN(3,1)*clhs230;
const double clhs247 = DN(3,2)*clhs230;
const double clhs248 = DN(1,0)*N[3];
const double clhs249 = DN(3,0)*N[1];
const double clhs250 = clhs226 + clhs48;
const double clhs251 = pow(DN(1,1), 2);
const double clhs252 = DN(1,1)*clhs10;
const double clhs253 = DN(1,2)*clhs252;
const double clhs254 = DN(2,0)*clhs252;
const double clhs255 = DN(1,1)*DN(2,1);
const double clhs256 = clhs10*clhs255;
const double clhs257 = clhs235 + clhs237;
const double clhs258 = DN(2,2)*clhs252;
const double clhs259 = DN(1,1)*N[2];
const double clhs260 = DN(2,1)*N[1];
const double clhs261 = DN(3,0)*clhs252;
const double clhs262 = DN(1,1)*DN(3,1);
const double clhs263 = clhs10*clhs262;
const double clhs264 = clhs243 + clhs245;
const double clhs265 = DN(3,2)*clhs252;
const double clhs266 = DN(1,1)*N[3];
const double clhs267 = DN(3,1)*N[1];
const double clhs268 = pow(DN(1,2), 2);
const double clhs269 = DN(1,2)*clhs10;
const double clhs270 = DN(2,0)*clhs269;
const double clhs271 = DN(2,1)*clhs269;
const double clhs272 = DN(1,2)*DN(2,2);
const double clhs273 = clhs10*clhs272;
const double clhs274 = DN(1,2)*N[2];
const double clhs275 = DN(2,2)*N[1];
const double clhs276 = DN(3,0)*clhs269;
const double clhs277 = DN(3,1)*clhs269;
const double clhs278 = DN(1,2)*DN(3,2);
const double clhs279 = clhs10*clhs278;
const double clhs280 = DN(1,2)*N[3];
const double clhs281 = DN(3,2)*N[1];
const double clhs282 = clhs16*clhs22;
const double clhs283 = rho*(N[1] + clhs17*(clhs221 + clhs223 + 1.0*clhs52));
const double clhs284 = clhs23*(clhs234 + clhs255 + clhs272);
const double clhs285 = clhs23*(clhs242 + clhs262 + clhs278);
const double clhs286 = 1.0*clhs78;
const double clhs287 = clhs17*clhs286;
const double clhs288 = 1.0*clhs79;
const double clhs289 = clhs17*clhs288;
const double clhs290 = N[2]*clhs24;
const double clhs291 = N[2]*clhs12 + clhs16*clhs287 - clhs16*clhs289 + clhs16*clhs290;
const double clhs292 = N[2]*clhs50 + clhs287*clhs53 - clhs289*clhs53 + clhs290*clhs53;
const double clhs293 = pow(DN(2,0), 2);
const double clhs294 = pow(N[2], 2);
const double clhs295 = K_darcy*clhs294 + N[2]*clhs78 + clhs13*clhs294 + clhs287*clhs81 - clhs289*clhs81 + clhs290*clhs81;
const double clhs296 = DN(2,0)*clhs10;
const double clhs297 = DN(2,1)*clhs296;
const double clhs298 = DN(2,2)*clhs296;
const double clhs299 = -N[2] + clhs287 - clhs289 + clhs290;
const double clhs300 = DN(2,0)*DN(3,0);
const double clhs301 = N[3]*clhs79 + N[3]*clhs80;
const double clhs302 = clhs10*clhs300 + clhs301;
const double clhs303 = N[2]*clhs106 + clhs109*clhs287 - clhs109*clhs289 + clhs109*clhs290;
const double clhs304 = DN(3,1)*clhs296;
const double clhs305 = DN(3,2)*clhs296;
const double clhs306 = DN(2,0)*N[3];
const double clhs307 = DN(3,0)*N[2];
const double clhs308 = clhs291 + clhs76;
const double clhs309 = clhs235 + clhs292;
const double clhs310 = pow(DN(2,1), 2);
const double clhs311 = DN(2,1)*clhs10;
const double clhs312 = DN(2,2)*clhs311;
const double clhs313 = DN(3,0)*clhs311;
const double clhs314 = DN(2,1)*DN(3,1);
const double clhs315 = clhs10*clhs314;
const double clhs316 = clhs301 + clhs303;
const double clhs317 = DN(3,2)*clhs311;
const double clhs318 = DN(2,1)*N[3];
const double clhs319 = DN(3,1)*N[2];
const double clhs320 = pow(DN(2,2), 2);
const double clhs321 = DN(2,2)*clhs10;
const double clhs322 = DN(3,0)*clhs321;
const double clhs323 = DN(3,1)*clhs321;
const double clhs324 = DN(2,2)*DN(3,2);
const double clhs325 = clhs10*clhs324;
const double clhs326 = DN(2,2)*N[3];
const double clhs327 = DN(3,2)*N[2];
const double clhs328 = rho*(N[2] + clhs17*(clhs286 + clhs288 + 1.0*clhs80));
const double clhs329 = clhs23*(clhs300 + clhs314 + clhs324);
const double clhs330 = 1.0*clhs106;
const double clhs331 = clhs17*clhs330;
const double clhs332 = 1.0*clhs107;
const double clhs333 = clhs17*clhs332;
const double clhs334 = N[3]*clhs24;
const double clhs335 = N[3]*clhs12 + clhs16*clhs331 - clhs16*clhs333 + clhs16*clhs334;
const double clhs336 = N[3]*clhs50 + clhs331*clhs53 - clhs333*clhs53 + clhs334*clhs53;
const double clhs337 = N[3]*clhs78 + clhs331*clhs81 - clhs333*clhs81 + clhs334*clhs81;
const double clhs338 = pow(DN(3,0), 2);
const double clhs339 = pow(N[3], 2);
const double clhs340 = K_darcy*clhs339 + N[3]*clhs106 + clhs109*clhs331 - clhs109*clhs333 + clhs109*clhs334 + clhs13*clhs339;
const double clhs341 = DN(3,0)*clhs10;
const double clhs342 = DN(3,1)*clhs341;
const double clhs343 = DN(3,2)*clhs341;
const double clhs344 = -N[3] + clhs331 - clhs333 + clhs334;
const double clhs345 = clhs104 + clhs335;
const double clhs346 = clhs243 + clhs336;
const double clhs347 = clhs301 + clhs337;
const double clhs348 = pow(DN(3,1), 2);
const double clhs349 = DN(3,1)*DN(3,2)*clhs10;
const double clhs350 = pow(DN(3,2), 2);
const double clhs351 = rho*(N[3] + clhs17*(1.0*clhs108 + clhs330 + clhs332));
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs26;
lhs(0,1)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + DN(0,2)*clhs32 + clhs34;
lhs(0,2)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + DN(0,2)*clhs39 + clhs40;
lhs(0,3)=DN(0,0)*clhs41;
lhs(0,4)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + DN(0,2)*clhs46 + clhs49 + clhs54;
lhs(0,5)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs61;
lhs(0,6)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs67;
lhs(0,7)=DN(1,0)*clhs19 - DN(1,0)*clhs21 + clhs24*clhs69 - clhs68;
lhs(0,8)=DN(0,0)*clhs70 + DN(0,1)*clhs72 + DN(0,2)*clhs74 + clhs77 + clhs82;
lhs(0,9)=DN(0,0)*clhs83 + DN(0,1)*clhs85 + DN(0,2)*clhs88 + clhs89;
lhs(0,10)=DN(0,0)*clhs90 + DN(0,1)*clhs92 + DN(0,2)*clhs94 + clhs95;
lhs(0,11)=DN(2,0)*clhs19 - DN(2,0)*clhs21 + clhs24*clhs97 - clhs96;
lhs(0,12)=DN(0,0)*clhs98 + DN(0,1)*clhs100 + DN(0,2)*clhs102 + clhs105 + clhs110;
lhs(0,13)=DN(0,0)*clhs111 + DN(0,1)*clhs113 + DN(0,2)*clhs116 + clhs117;
lhs(0,14)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs123;
lhs(0,15)=DN(3,0)*clhs19 - DN(3,0)*clhs21 - clhs124 + clhs125*clhs24;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs126 + DN(0,2)*clhs127 + clhs34;
lhs(1,1)=DN(0,0)*clhs29 + DN(0,1)*clhs128 + DN(0,2)*clhs130 + clhs10*clhs131 + clhs26;
lhs(1,2)=DN(0,0)*clhs37 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs136;
lhs(1,3)=DN(0,1)*clhs41;
lhs(1,4)=DN(0,0)*clhs44 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
lhs(1,5)=DN(0,0)*clhs57 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs144 + clhs145;
lhs(1,6)=DN(0,0)*clhs64 + DN(0,1)*clhs146 + DN(0,2)*clhs148 + clhs149;
lhs(1,7)=DN(1,1)*clhs19 - DN(1,1)*clhs21 - clhs150 + clhs151*clhs24;
lhs(1,8)=DN(0,0)*clhs72 + DN(0,1)*clhs152 + DN(0,2)*clhs153 + clhs154;
lhs(1,9)=DN(0,0)*clhs85 + DN(0,1)*clhs155 + DN(0,2)*clhs157 + clhs159 + clhs160;
lhs(1,10)=DN(0,0)*clhs92 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs164;
lhs(1,11)=DN(2,1)*clhs19 - DN(2,1)*clhs21 - clhs165 + clhs166*clhs24;
lhs(1,12)=DN(0,0)*clhs100 + DN(0,1)*clhs167 + DN(0,2)*clhs168 + clhs169;
lhs(1,13)=DN(0,0)*clhs113 + DN(0,1)*clhs170 + DN(0,2)*clhs172 + clhs174 + clhs175;
lhs(1,14)=DN(0,0)*clhs120 + DN(0,1)*clhs176 + DN(0,2)*clhs178 + clhs179;
lhs(1,15)=DN(3,1)*clhs19 - DN(3,1)*clhs21 - clhs180 + clhs181*clhs24;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs182 + clhs40;
lhs(2,1)=DN(0,0)*clhs32 + DN(0,1)*clhs130 + DN(0,2)*clhs183 + clhs136;
lhs(2,2)=DN(0,0)*clhs39 + DN(0,1)*clhs134 + DN(0,2)*clhs184 + clhs10*clhs185 + clhs26;
lhs(2,3)=DN(0,2)*clhs41;
lhs(2,4)=DN(0,0)*clhs46 + DN(0,1)*clhs138 + DN(0,2)*clhs186 + clhs188;
lhs(2,5)=DN(0,0)*clhs60 + DN(0,1)*clhs142 + DN(0,2)*clhs189 + clhs190;
lhs(2,6)=DN(0,0)*clhs66 + DN(0,1)*clhs148 + DN(0,2)*clhs191 + clhs145 + clhs193;
lhs(2,7)=DN(1,2)*clhs19 - DN(1,2)*clhs21 - clhs194 + clhs195*clhs24;
lhs(2,8)=DN(0,0)*clhs74 + DN(0,1)*clhs153 + DN(0,2)*clhs196 + clhs197;
lhs(2,9)=DN(0,0)*clhs88 + DN(0,1)*clhs157 + DN(0,2)*clhs198 + clhs199;
lhs(2,10)=DN(0,0)*clhs94 + DN(0,1)*clhs163 + DN(0,2)*clhs200 + clhs160 + clhs202;
lhs(2,11)=DN(2,2)*clhs19 - DN(2,2)*clhs21 - clhs203 + clhs204*clhs24;
lhs(2,12)=DN(0,0)*clhs102 + DN(0,1)*clhs168 + DN(0,2)*clhs205 + clhs206;
lhs(2,13)=DN(0,0)*clhs116 + DN(0,1)*clhs172 + DN(0,2)*clhs207 + clhs208;
lhs(2,14)=DN(0,0)*clhs122 + DN(0,1)*clhs178 + DN(0,2)*clhs209 + clhs175 + clhs211;
lhs(2,15)=DN(3,2)*clhs19 - DN(3,2)*clhs21 - clhs212 + clhs213*clhs24;
lhs(3,0)=DN(0,0)*clhs214;
lhs(3,1)=DN(0,1)*clhs214;
lhs(3,2)=DN(0,2)*clhs214;
lhs(3,3)=clhs23*(clhs131 + clhs185 + clhs5);
lhs(3,4)=rho*(DN(0,0)*clhs215 + clhs69);
lhs(3,5)=rho*(DN(0,1)*clhs215 + clhs151);
lhs(3,6)=rho*(DN(0,2)*clhs215 + clhs195);
lhs(3,7)=clhs216;
lhs(3,8)=rho*(DN(0,0)*clhs217 + clhs97);
lhs(3,9)=rho*(DN(0,1)*clhs217 + clhs166);
lhs(3,10)=rho*(DN(0,2)*clhs217 + clhs204);
lhs(3,11)=clhs218;
lhs(3,12)=rho*(DN(0,0)*clhs219 + clhs125);
lhs(3,13)=rho*(DN(0,1)*clhs219 + clhs181);
lhs(3,14)=rho*(DN(0,2)*clhs219 + clhs213);
lhs(3,15)=clhs220;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs226 + clhs49;
lhs(4,1)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + DN(1,2)*clhs32 + clhs139;
lhs(4,2)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + DN(1,2)*clhs39 + clhs188;
lhs(4,3)=DN(0,0)*clhs222 - DN(0,0)*clhs224 + clhs24*clhs68 - clhs69;
lhs(4,4)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + DN(1,2)*clhs46 + clhs10*clhs227 + clhs229;
lhs(4,5)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs231;
lhs(4,6)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs232;
lhs(4,7)=DN(1,0)*clhs233;
lhs(4,8)=DN(1,0)*clhs70 + DN(1,1)*clhs72 + DN(1,2)*clhs74 + clhs236 + clhs237;
lhs(4,9)=DN(1,0)*clhs83 + DN(1,1)*clhs85 + DN(1,2)*clhs88 + clhs238;
lhs(4,10)=DN(1,0)*clhs90 + DN(1,1)*clhs92 + DN(1,2)*clhs94 + clhs239;
lhs(4,11)=DN(2,0)*clhs222 - DN(2,0)*clhs224 + clhs24*clhs241 - clhs240;
lhs(4,12)=DN(1,0)*clhs98 + DN(1,1)*clhs100 + DN(1,2)*clhs102 + clhs244 + clhs245;
lhs(4,13)=DN(1,0)*clhs111 + DN(1,1)*clhs113 + DN(1,2)*clhs116 + clhs246;
lhs(4,14)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs247;
lhs(4,15)=DN(3,0)*clhs222 - DN(3,0)*clhs224 + clhs24*clhs249 - clhs248;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs126 + DN(1,2)*clhs127 + clhs61;
lhs(5,1)=DN(1,0)*clhs29 + DN(1,1)*clhs128 + DN(1,2)*clhs130 + clhs144 + clhs250;
lhs(5,2)=DN(1,0)*clhs37 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs190;
lhs(5,3)=DN(0,1)*clhs222 - DN(0,1)*clhs224 + clhs150*clhs24 - clhs151;
lhs(5,4)=DN(1,0)*clhs44 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs231;
lhs(5,5)=DN(1,0)*clhs57 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs10*clhs251 + clhs229;
lhs(5,6)=DN(1,0)*clhs64 + DN(1,1)*clhs146 + DN(1,2)*clhs148 + clhs253;
lhs(5,7)=DN(1,1)*clhs233;
lhs(5,8)=DN(1,0)*clhs72 + DN(1,1)*clhs152 + DN(1,2)*clhs153 + clhs254;
lhs(5,9)=DN(1,0)*clhs85 + DN(1,1)*clhs155 + DN(1,2)*clhs157 + clhs256 + clhs257;
lhs(5,10)=DN(1,0)*clhs92 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs258;
lhs(5,11)=DN(2,1)*clhs222 - DN(2,1)*clhs224 + clhs24*clhs260 - clhs259;
lhs(5,12)=DN(1,0)*clhs100 + DN(1,1)*clhs167 + DN(1,2)*clhs168 + clhs261;
lhs(5,13)=DN(1,0)*clhs113 + DN(1,1)*clhs170 + DN(1,2)*clhs172 + clhs263 + clhs264;
lhs(5,14)=DN(1,0)*clhs120 + DN(1,1)*clhs176 + DN(1,2)*clhs178 + clhs265;
lhs(5,15)=DN(3,1)*clhs222 - DN(3,1)*clhs224 + clhs24*clhs267 - clhs266;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs182 + clhs67;
lhs(6,1)=DN(1,0)*clhs32 + DN(1,1)*clhs130 + DN(1,2)*clhs183 + clhs149;
lhs(6,2)=DN(1,0)*clhs39 + DN(1,1)*clhs134 + DN(1,2)*clhs184 + clhs193 + clhs250;
lhs(6,3)=DN(0,2)*clhs222 - DN(0,2)*clhs224 + clhs194*clhs24 - clhs195;
lhs(6,4)=DN(1,0)*clhs46 + DN(1,1)*clhs138 + DN(1,2)*clhs186 + clhs232;
lhs(6,5)=DN(1,0)*clhs60 + DN(1,1)*clhs142 + DN(1,2)*clhs189 + clhs253;
lhs(6,6)=DN(1,0)*clhs66 + DN(1,1)*clhs148 + DN(1,2)*clhs191 + clhs10*clhs268 + clhs229;
lhs(6,7)=DN(1,2)*clhs233;
lhs(6,8)=DN(1,0)*clhs74 + DN(1,1)*clhs153 + DN(1,2)*clhs196 + clhs270;
lhs(6,9)=DN(1,0)*clhs88 + DN(1,1)*clhs157 + DN(1,2)*clhs198 + clhs271;
lhs(6,10)=DN(1,0)*clhs94 + DN(1,1)*clhs163 + DN(1,2)*clhs200 + clhs257 + clhs273;
lhs(6,11)=DN(2,2)*clhs222 - DN(2,2)*clhs224 + clhs24*clhs275 - clhs274;
lhs(6,12)=DN(1,0)*clhs102 + DN(1,1)*clhs168 + DN(1,2)*clhs205 + clhs276;
lhs(6,13)=DN(1,0)*clhs116 + DN(1,1)*clhs172 + DN(1,2)*clhs207 + clhs277;
lhs(6,14)=DN(1,0)*clhs122 + DN(1,1)*clhs178 + DN(1,2)*clhs209 + clhs264 + clhs279;
lhs(6,15)=DN(3,2)*clhs222 - DN(3,2)*clhs224 + clhs24*clhs281 - clhs280;
lhs(7,0)=rho*(DN(1,0)*clhs282 + clhs68);
lhs(7,1)=rho*(DN(1,1)*clhs282 + clhs150);
lhs(7,2)=rho*(DN(1,2)*clhs282 + clhs194);
lhs(7,3)=clhs216;
lhs(7,4)=DN(1,0)*clhs283;
lhs(7,5)=DN(1,1)*clhs283;
lhs(7,6)=DN(1,2)*clhs283;
lhs(7,7)=clhs23*(clhs227 + clhs251 + clhs268);
lhs(7,8)=rho*(DN(1,0)*clhs217 + clhs241);
lhs(7,9)=rho*(DN(1,1)*clhs217 + clhs260);
lhs(7,10)=rho*(DN(1,2)*clhs217 + clhs275);
lhs(7,11)=clhs284;
lhs(7,12)=rho*(DN(1,0)*clhs219 + clhs249);
lhs(7,13)=rho*(DN(1,1)*clhs219 + clhs267);
lhs(7,14)=rho*(DN(1,2)*clhs219 + clhs281);
lhs(7,15)=clhs285;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs291 + clhs77;
lhs(8,1)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + DN(2,2)*clhs32 + clhs154;
lhs(8,2)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + DN(2,2)*clhs39 + clhs197;
lhs(8,3)=DN(0,0)*clhs287 - DN(0,0)*clhs289 + clhs24*clhs96 - clhs97;
lhs(8,4)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + DN(2,2)*clhs46 + clhs236 + clhs292;
lhs(8,5)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs254;
lhs(8,6)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs270;
lhs(8,7)=DN(1,0)*clhs287 - DN(1,0)*clhs289 + clhs24*clhs240 - clhs241;
lhs(8,8)=DN(2,0)*clhs70 + DN(2,1)*clhs72 + DN(2,2)*clhs74 + clhs10*clhs293 + clhs295;
lhs(8,9)=DN(2,0)*clhs83 + DN(2,1)*clhs85 + DN(2,2)*clhs88 + clhs297;
lhs(8,10)=DN(2,0)*clhs90 + DN(2,1)*clhs92 + DN(2,2)*clhs94 + clhs298;
lhs(8,11)=DN(2,0)*clhs299;
lhs(8,12)=DN(2,0)*clhs98 + DN(2,1)*clhs100 + DN(2,2)*clhs102 + clhs302 + clhs303;
lhs(8,13)=DN(2,0)*clhs111 + DN(2,1)*clhs113 + DN(2,2)*clhs116 + clhs304;
lhs(8,14)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs305;
lhs(8,15)=DN(3,0)*clhs287 - DN(3,0)*clhs289 + clhs24*clhs307 - clhs306;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs126 + DN(2,2)*clhs127 + clhs89;
lhs(9,1)=DN(2,0)*clhs29 + DN(2,1)*clhs128 + DN(2,2)*clhs130 + clhs159 + clhs308;
lhs(9,2)=DN(2,0)*clhs37 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs199;
lhs(9,3)=DN(0,1)*clhs287 - DN(0,1)*clhs289 + clhs165*clhs24 - clhs166;
lhs(9,4)=DN(2,0)*clhs44 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs238;
lhs(9,5)=DN(2,0)*clhs57 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs256 + clhs309;
lhs(9,6)=DN(2,0)*clhs64 + DN(2,1)*clhs146 + DN(2,2)*clhs148 + clhs271;
lhs(9,7)=DN(1,1)*clhs287 - DN(1,1)*clhs289 + clhs24*clhs259 - clhs260;
lhs(9,8)=DN(2,0)*clhs72 + DN(2,1)*clhs152 + DN(2,2)*clhs153 + clhs297;
lhs(9,9)=DN(2,0)*clhs85 + DN(2,1)*clhs155 + DN(2,2)*clhs157 + clhs10*clhs310 + clhs295;
lhs(9,10)=DN(2,0)*clhs92 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs312;
lhs(9,11)=DN(2,1)*clhs299;
lhs(9,12)=DN(2,0)*clhs100 + DN(2,1)*clhs167 + DN(2,2)*clhs168 + clhs313;
lhs(9,13)=DN(2,0)*clhs113 + DN(2,1)*clhs170 + DN(2,2)*clhs172 + clhs315 + clhs316;
lhs(9,14)=DN(2,0)*clhs120 + DN(2,1)*clhs176 + DN(2,2)*clhs178 + clhs317;
lhs(9,15)=DN(3,1)*clhs287 - DN(3,1)*clhs289 + clhs24*clhs319 - clhs318;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs182 + clhs95;
lhs(10,1)=DN(2,0)*clhs32 + DN(2,1)*clhs130 + DN(2,2)*clhs183 + clhs164;
lhs(10,2)=DN(2,0)*clhs39 + DN(2,1)*clhs134 + DN(2,2)*clhs184 + clhs202 + clhs308;
lhs(10,3)=DN(0,2)*clhs287 - DN(0,2)*clhs289 + clhs203*clhs24 - clhs204;
lhs(10,4)=DN(2,0)*clhs46 + DN(2,1)*clhs138 + DN(2,2)*clhs186 + clhs239;
lhs(10,5)=DN(2,0)*clhs60 + DN(2,1)*clhs142 + DN(2,2)*clhs189 + clhs258;
lhs(10,6)=DN(2,0)*clhs66 + DN(2,1)*clhs148 + DN(2,2)*clhs191 + clhs273 + clhs309;
lhs(10,7)=DN(1,2)*clhs287 - DN(1,2)*clhs289 + clhs24*clhs274 - clhs275;
lhs(10,8)=DN(2,0)*clhs74 + DN(2,1)*clhs153 + DN(2,2)*clhs196 + clhs298;
lhs(10,9)=DN(2,0)*clhs88 + DN(2,1)*clhs157 + DN(2,2)*clhs198 + clhs312;
lhs(10,10)=DN(2,0)*clhs94 + DN(2,1)*clhs163 + DN(2,2)*clhs200 + clhs10*clhs320 + clhs295;
lhs(10,11)=DN(2,2)*clhs299;
lhs(10,12)=DN(2,0)*clhs102 + DN(2,1)*clhs168 + DN(2,2)*clhs205 + clhs322;
lhs(10,13)=DN(2,0)*clhs116 + DN(2,1)*clhs172 + DN(2,2)*clhs207 + clhs323;
lhs(10,14)=DN(2,0)*clhs122 + DN(2,1)*clhs178 + DN(2,2)*clhs209 + clhs316 + clhs325;
lhs(10,15)=DN(3,2)*clhs287 - DN(3,2)*clhs289 + clhs24*clhs327 - clhs326;
lhs(11,0)=rho*(DN(2,0)*clhs282 + clhs96);
lhs(11,1)=rho*(DN(2,1)*clhs282 + clhs165);
lhs(11,2)=rho*(DN(2,2)*clhs282 + clhs203);
lhs(11,3)=clhs218;
lhs(11,4)=rho*(DN(2,0)*clhs215 + clhs240);
lhs(11,5)=rho*(DN(2,1)*clhs215 + clhs259);
lhs(11,6)=rho*(DN(2,2)*clhs215 + clhs274);
lhs(11,7)=clhs284;
lhs(11,8)=DN(2,0)*clhs328;
lhs(11,9)=DN(2,1)*clhs328;
lhs(11,10)=DN(2,2)*clhs328;
lhs(11,11)=clhs23*(clhs293 + clhs310 + clhs320);
lhs(11,12)=rho*(DN(2,0)*clhs219 + clhs307);
lhs(11,13)=rho*(DN(2,1)*clhs219 + clhs319);
lhs(11,14)=rho*(DN(2,2)*clhs219 + clhs327);
lhs(11,15)=clhs329;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs105 + clhs335;
lhs(12,1)=DN(3,0)*clhs27 + DN(3,1)*clhs29 + DN(3,2)*clhs32 + clhs169;
lhs(12,2)=DN(3,0)*clhs35 + DN(3,1)*clhs37 + DN(3,2)*clhs39 + clhs206;
lhs(12,3)=DN(0,0)*clhs331 - DN(0,0)*clhs333 + clhs124*clhs24 - clhs125;
lhs(12,4)=DN(3,0)*clhs42 + DN(3,1)*clhs44 + DN(3,2)*clhs46 + clhs244 + clhs336;
lhs(12,5)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs261;
lhs(12,6)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs276;
lhs(12,7)=DN(1,0)*clhs331 - DN(1,0)*clhs333 + clhs24*clhs248 - clhs249;
lhs(12,8)=DN(3,0)*clhs70 + DN(3,1)*clhs72 + DN(3,2)*clhs74 + clhs302 + clhs337;
lhs(12,9)=DN(3,0)*clhs83 + DN(3,1)*clhs85 + DN(3,2)*clhs88 + clhs313;
lhs(12,10)=DN(3,0)*clhs90 + DN(3,1)*clhs92 + DN(3,2)*clhs94 + clhs322;
lhs(12,11)=DN(2,0)*clhs331 - DN(2,0)*clhs333 + clhs24*clhs306 - clhs307;
lhs(12,12)=DN(3,0)*clhs98 + DN(3,1)*clhs100 + DN(3,2)*clhs102 + clhs10*clhs338 + clhs340;
lhs(12,13)=DN(3,0)*clhs111 + DN(3,1)*clhs113 + DN(3,2)*clhs116 + clhs342;
lhs(12,14)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs343;
lhs(12,15)=DN(3,0)*clhs344;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs126 + DN(3,2)*clhs127 + clhs117;
lhs(13,1)=DN(3,0)*clhs29 + DN(3,1)*clhs128 + DN(3,2)*clhs130 + clhs174 + clhs345;
lhs(13,2)=DN(3,0)*clhs37 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs208;
lhs(13,3)=DN(0,1)*clhs331 - DN(0,1)*clhs333 + clhs180*clhs24 - clhs181;
lhs(13,4)=DN(3,0)*clhs44 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs246;
lhs(13,5)=DN(3,0)*clhs57 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs263 + clhs346;
lhs(13,6)=DN(3,0)*clhs64 + DN(3,1)*clhs146 + DN(3,2)*clhs148 + clhs277;
lhs(13,7)=DN(1,1)*clhs331 - DN(1,1)*clhs333 + clhs24*clhs266 - clhs267;
lhs(13,8)=DN(3,0)*clhs72 + DN(3,1)*clhs152 + DN(3,2)*clhs153 + clhs304;
lhs(13,9)=DN(3,0)*clhs85 + DN(3,1)*clhs155 + DN(3,2)*clhs157 + clhs315 + clhs347;
lhs(13,10)=DN(3,0)*clhs92 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs323;
lhs(13,11)=DN(2,1)*clhs331 - DN(2,1)*clhs333 + clhs24*clhs318 - clhs319;
lhs(13,12)=DN(3,0)*clhs100 + DN(3,1)*clhs167 + DN(3,2)*clhs168 + clhs342;
lhs(13,13)=DN(3,0)*clhs113 + DN(3,1)*clhs170 + DN(3,2)*clhs172 + clhs10*clhs348 + clhs340;
lhs(13,14)=DN(3,0)*clhs120 + DN(3,1)*clhs176 + DN(3,2)*clhs178 + clhs349;
lhs(13,15)=DN(3,1)*clhs344;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs182 + clhs123;
lhs(14,1)=DN(3,0)*clhs32 + DN(3,1)*clhs130 + DN(3,2)*clhs183 + clhs179;
lhs(14,2)=DN(3,0)*clhs39 + DN(3,1)*clhs134 + DN(3,2)*clhs184 + clhs211 + clhs345;
lhs(14,3)=DN(0,2)*clhs331 - DN(0,2)*clhs333 + clhs212*clhs24 - clhs213;
lhs(14,4)=DN(3,0)*clhs46 + DN(3,1)*clhs138 + DN(3,2)*clhs186 + clhs247;
lhs(14,5)=DN(3,0)*clhs60 + DN(3,1)*clhs142 + DN(3,2)*clhs189 + clhs265;
lhs(14,6)=DN(3,0)*clhs66 + DN(3,1)*clhs148 + DN(3,2)*clhs191 + clhs279 + clhs346;
lhs(14,7)=DN(1,2)*clhs331 - DN(1,2)*clhs333 + clhs24*clhs280 - clhs281;
lhs(14,8)=DN(3,0)*clhs74 + DN(3,1)*clhs153 + DN(3,2)*clhs196 + clhs305;
lhs(14,9)=DN(3,0)*clhs88 + DN(3,1)*clhs157 + DN(3,2)*clhs198 + clhs317;
lhs(14,10)=DN(3,0)*clhs94 + DN(3,1)*clhs163 + DN(3,2)*clhs200 + clhs325 + clhs347;
lhs(14,11)=DN(2,2)*clhs331 - DN(2,2)*clhs333 + clhs24*clhs326 - clhs327;
lhs(14,12)=DN(3,0)*clhs102 + DN(3,1)*clhs168 + DN(3,2)*clhs205 + clhs343;
lhs(14,13)=DN(3,0)*clhs116 + DN(3,1)*clhs172 + DN(3,2)*clhs207 + clhs349;
lhs(14,14)=DN(3,0)*clhs122 + DN(3,1)*clhs178 + DN(3,2)*clhs209 + clhs10*clhs350 + clhs340;
lhs(14,15)=DN(3,2)*clhs344;
lhs(15,0)=rho*(DN(3,0)*clhs282 + clhs124);
lhs(15,1)=rho*(DN(3,1)*clhs282 + clhs180);
lhs(15,2)=rho*(DN(3,2)*clhs282 + clhs212);
lhs(15,3)=clhs220;
lhs(15,4)=rho*(DN(3,0)*clhs215 + clhs248);
lhs(15,5)=rho*(DN(3,1)*clhs215 + clhs266);
lhs(15,6)=rho*(DN(3,2)*clhs215 + clhs280);
lhs(15,7)=clhs285;
lhs(15,8)=rho*(DN(3,0)*clhs217 + clhs306);
lhs(15,9)=rho*(DN(3,1)*clhs217 + clhs318);
lhs(15,10)=rho*(DN(3,2)*clhs217 + clhs326);
lhs(15,11)=clhs329;
lhs(15,12)=DN(3,0)*clhs351;
lhs(15,13)=DN(3,1)*clhs351;
lhs(15,14)=DN(3,2)*clhs351;
lhs(15,15)=clhs23*(clhs338 + clhs348 + clhs350);


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
        const double volume_error =-rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 = K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs3 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs4 = DN(0,0)*v(0,0);
const double crhs5 = DN(1,0)*v(1,0);
const double crhs6 = DN(2,0)*v(2,0);
const double crhs7 = crhs4 + crhs5 + crhs6;
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs10 = rho*(crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs11 = rho*stab_c2*sqrt(pow(crhs8, 2) + pow(crhs9, 2));
const double crhs12 = DN(0,1)*v(0,1);
const double crhs13 = DN(1,1)*v(1,1);
const double crhs14 = DN(2,1)*v(2,1);
const double crhs15 = rho*(crhs11*h/stab_c1 + mu)*(-crhs12 - crhs13 - crhs14 - crhs4 - crhs5 - crhs6 + volume_error_ratio);
const double crhs16 = 1.0/(K_darcy + crhs11/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs17 = crhs16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs10 + crhs2 + crhs3);
const double crhs18 = K_darcy*N[0];
const double crhs19 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs20 = N[0]*crhs19;
const double crhs21 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs22 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs23 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs24 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs25 = crhs12 + crhs13 + crhs14;
const double crhs26 = rho*(crhs25*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs27 = crhs16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs22 + crhs23 + crhs24 + crhs26);
const double crhs28 = crhs25 + crhs7 - volume_error_ratio;
const double crhs29 = K_darcy*N[1];
const double crhs30 = N[1]*crhs19;
const double crhs31 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs32 = K_darcy*N[2];
const double crhs33 = N[2]*crhs19;
const double crhs34 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs15 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - N[0]*crhs3 + crhs17*crhs18 - crhs17*crhs20 - crhs17*crhs21;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs15 - DN(0,1)*stress[1] + N[0]*crhs22 - N[0]*crhs23 - N[0]*crhs24 - N[0]*crhs26 + crhs18*crhs27 - crhs20*crhs27 - crhs21*crhs27;
rhs[2]=-rho*(DN(0,0)*crhs17 + DN(0,1)*crhs27 + N[0]*crhs28);
rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs15 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - N[1]*crhs3 + crhs17*crhs29 - crhs17*crhs30 - crhs17*crhs31;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs15 - DN(1,1)*stress[1] + N[1]*crhs22 - N[1]*crhs23 - N[1]*crhs24 - N[1]*crhs26 + crhs27*crhs29 - crhs27*crhs30 - crhs27*crhs31;
rhs[5]=-rho*(DN(1,0)*crhs17 + DN(1,1)*crhs27 + N[1]*crhs28);
rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs15 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - N[2]*crhs3 + crhs17*crhs32 - crhs17*crhs33 - crhs17*crhs34;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs15 - DN(2,1)*stress[1] + N[2]*crhs22 - N[2]*crhs23 - N[2]*crhs24 - N[2]*crhs26 + crhs27*crhs32 - crhs27*crhs33 - crhs27*crhs34;
rhs[8]=-rho*(DN(2,0)*crhs17 + DN(2,1)*crhs27 + N[2]*crhs28);


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
        const double volume_error = -rData.VolumeError;
        volume_error_ratio = volume_error / dt;
    }

    auto &rhs = rData.rhs;

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 = K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crhs3 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs4 = DN(0,0)*v(0,0);
const double crhs5 = DN(1,0)*v(1,0);
const double crhs6 = DN(2,0)*v(2,0);
const double crhs7 = DN(3,0)*v(3,0);
const double crhs8 = crhs4 + crhs5 + crhs6 + crhs7;
const double crhs9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs11 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs12 = rho*(crhs10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs11*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*crhs9);
const double crhs13 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs11, 2) + pow(crhs9, 2));
const double crhs14 = DN(0,1)*v(0,1);
const double crhs15 = DN(0,2)*v(0,2);
const double crhs16 = DN(1,1)*v(1,1);
const double crhs17 = DN(1,2)*v(1,2);
const double crhs18 = DN(2,1)*v(2,1);
const double crhs19 = DN(2,2)*v(2,2);
const double crhs20 = DN(3,1)*v(3,1);
const double crhs21 = DN(3,2)*v(3,2);
const double crhs22 = rho*(crhs13*h/stab_c1 + mu)*(-crhs14 - crhs15 - crhs16 - crhs17 - crhs18 - crhs19 - crhs20 - crhs21 - crhs4 - crhs5 - crhs6 - crhs7 + volume_error_ratio);
const double crhs23 = 1.0/(K_darcy + crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs24 = crhs23*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs12 + crhs2 + crhs3);
const double crhs25 = K_darcy*N[0];
const double crhs26 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs27 = N[0]*crhs26;
const double crhs28 = rho*(DN(0,0)*crhs9 + DN(0,1)*crhs10 + DN(0,2)*crhs11);
const double crhs29 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs30 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs31 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs32 = crhs14 + crhs16 + crhs18 + crhs20;
const double crhs33 = rho*(crhs10*crhs32 + crhs11*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs9*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)));
const double crhs34 = crhs23*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs29 + crhs30 + crhs31 + crhs33);
const double crhs35 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs36 = K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs37 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs38 = crhs15 + crhs17 + crhs19 + crhs21;
const double crhs39 = rho*(crhs10*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs11*crhs38 + crhs9*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)));
const double crhs40 = crhs23*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs35 + crhs36 + crhs37 + crhs39);
const double crhs41 = crhs32 + crhs38 + crhs8 - volume_error_ratio;
const double crhs42 = K_darcy*N[1];
const double crhs43 = N[1]*crhs26;
const double crhs44 = rho*(DN(1,0)*crhs9 + DN(1,1)*crhs10 + DN(1,2)*crhs11);
const double crhs45 = K_darcy*N[2];
const double crhs46 = N[2]*crhs26;
const double crhs47 = rho*(DN(2,0)*crhs9 + DN(2,1)*crhs10 + DN(2,2)*crhs11);
const double crhs48 = K_darcy*N[3];
const double crhs49 = N[3]*crhs26;
const double crhs50 = rho*(DN(3,0)*crhs9 + DN(3,1)*crhs10 + DN(3,2)*crhs11);
rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs22 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs12 - N[0]*crhs2 - N[0]*crhs3 + crhs24*crhs25 - crhs24*crhs27 - crhs24*crhs28;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 + DN(0,1)*crhs22 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs29 - N[0]*crhs30 - N[0]*crhs31 - N[0]*crhs33 + crhs25*crhs34 - crhs27*crhs34 - crhs28*crhs34;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 + DN(0,2)*crhs22 - DN(0,2)*stress[2] + N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs37 - N[0]*crhs39 + crhs25*crhs40 - crhs27*crhs40 - crhs28*crhs40;
rhs[3]=-rho*(DN(0,0)*crhs24 + DN(0,1)*crhs34 + DN(0,2)*crhs40 + N[0]*crhs41);
rhs[4]=DN(1,0)*crhs0 + DN(1,0)*crhs22 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs12 - N[1]*crhs2 - N[1]*crhs3 + crhs24*crhs42 - crhs24*crhs43 - crhs24*crhs44;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 + DN(1,1)*crhs22 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs29 - N[1]*crhs30 - N[1]*crhs31 - N[1]*crhs33 + crhs34*crhs42 - crhs34*crhs43 - crhs34*crhs44;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 + DN(1,2)*crhs22 - DN(1,2)*stress[2] + N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs37 - N[1]*crhs39 + crhs40*crhs42 - crhs40*crhs43 - crhs40*crhs44;
rhs[7]=-rho*(DN(1,0)*crhs24 + DN(1,1)*crhs34 + DN(1,2)*crhs40 + N[1]*crhs41);
rhs[8]=DN(2,0)*crhs0 + DN(2,0)*crhs22 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs12 - N[2]*crhs2 - N[2]*crhs3 + crhs24*crhs45 - crhs24*crhs46 - crhs24*crhs47;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 + DN(2,1)*crhs22 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs29 - N[2]*crhs30 - N[2]*crhs31 - N[2]*crhs33 + crhs34*crhs45 - crhs34*crhs46 - crhs34*crhs47;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 + DN(2,2)*crhs22 - DN(2,2)*stress[2] + N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs37 - N[2]*crhs39 + crhs40*crhs45 - crhs40*crhs46 - crhs40*crhs47;
rhs[11]=-rho*(DN(2,0)*crhs24 + DN(2,1)*crhs34 + DN(2,2)*crhs40 + N[2]*crhs41);
rhs[12]=DN(3,0)*crhs0 + DN(3,0)*crhs22 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs12 - N[3]*crhs2 - N[3]*crhs3 + crhs24*crhs48 - crhs24*crhs49 - crhs24*crhs50;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 + DN(3,1)*crhs22 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs29 - N[3]*crhs30 - N[3]*crhs31 - N[3]*crhs33 + crhs34*crhs48 - crhs34*crhs49 - crhs34*crhs50;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 + DN(3,2)*crhs22 - DN(3,2)*stress[2] + N[3]*crhs35 - N[3]*crhs36 - N[3]*crhs37 - N[3]*crhs39 + crhs40*crhs48 - crhs40*crhs49 - crhs40*crhs50;
rhs[15]=-rho*(DN(3,0)*crhs24 + DN(3,1)*crhs34 + DN(3,2)*crhs40 + N[3]*crhs41);


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

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV2 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV3 = K_darcy*cV2;
const double cV4 = DNenr(0,0)*N[0];
const double cV5 = cV2*rho;
const double cV6 = cV5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV7 = cV5*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV8 = N[0]*cV3;
const double cV9 = N[0]*cV6;
const double cV10 = N[1]*cV3;
const double cV11 = N[1]*cV6;
const double cV12 = cV5*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV13 = N[2]*cV3;
const double cV14 = N[2]*cV6;
const double cV15 = cV5*(DN(2,0)*cV0 + DN(2,1)*cV1);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV7 - cV3*cV4 + cV4*cV6;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV7 - DNenr(1,0)*cV8 + DNenr(1,0)*cV9;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV7 - DNenr(2,0)*cV8 + DNenr(2,0)*cV9;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV7 - DNenr(0,1)*cV8 + DNenr(0,1)*cV9;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV7 - DNenr(1,1)*cV8 + DNenr(1,1)*cV9;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV7 - DNenr(2,1)*cV8 + DNenr(2,1)*cV9;
V(2,0)=cV5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=-DN(1,0)*Nenr[0] - DNenr(0,0)*cV10 + DNenr(0,0)*cV11 + DNenr(0,0)*cV12;
V(3,1)=-DN(1,0)*Nenr[1] - DNenr(1,0)*cV10 + DNenr(1,0)*cV11 + DNenr(1,0)*cV12;
V(3,2)=-DN(1,0)*Nenr[2] - DNenr(2,0)*cV10 + DNenr(2,0)*cV11 + DNenr(2,0)*cV12;
V(4,0)=-DN(1,1)*Nenr[0] - DNenr(0,1)*cV10 + DNenr(0,1)*cV11 + DNenr(0,1)*cV12;
V(4,1)=-DN(1,1)*Nenr[1] - DNenr(1,1)*cV10 + DNenr(1,1)*cV11 + DNenr(1,1)*cV12;
V(4,2)=-DN(1,1)*Nenr[2] - DNenr(2,1)*cV10 + DNenr(2,1)*cV11 + DNenr(2,1)*cV12;
V(5,0)=cV5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=-DN(2,0)*Nenr[0] - DNenr(0,0)*cV13 + DNenr(0,0)*cV14 + DNenr(0,0)*cV15;
V(6,1)=-DN(2,0)*Nenr[1] - DNenr(1,0)*cV13 + DNenr(1,0)*cV14 + DNenr(1,0)*cV15;
V(6,2)=-DN(2,0)*Nenr[2] - DNenr(2,0)*cV13 + DNenr(2,0)*cV14 + DNenr(2,0)*cV15;
V(7,0)=-DN(2,1)*Nenr[0] - DNenr(0,1)*cV13 + DNenr(0,1)*cV14 + DNenr(0,1)*cV15;
V(7,1)=-DN(2,1)*Nenr[1] - DNenr(1,1)*cV13 + DNenr(1,1)*cV14 + DNenr(1,1)*cV15;
V(7,2)=-DN(2,1)*Nenr[2] - DNenr(2,1)*cV13 + DNenr(2,1)*cV14 + DNenr(2,1)*cV15;
V(8,0)=cV5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
V(8,1)=cV5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
V(8,2)=cV5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH3 = cH2*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0));
const double cH4 = cH2*rho;
const double cH5 = cH2*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0));
const double cH6 = cH2*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0));
H(0,0)=rho*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH3);
H(0,1)=rho*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH3);
H(0,2)=cH4*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
H(0,3)=rho*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH5);
H(0,4)=rho*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH5);
H(0,5)=cH4*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
H(0,6)=rho*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH6);
H(0,7)=rho*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH6);
H(0,8)=cH4*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
H(1,0)=rho*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH3);
H(1,1)=rho*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH3);
H(1,2)=cH4*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
H(1,3)=rho*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH5);
H(1,4)=rho*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH5);
H(1,5)=cH4*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
H(1,6)=rho*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH6);
H(1,7)=rho*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH6);
H(1,8)=cH4*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
H(2,0)=rho*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH3);
H(2,1)=rho*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH3);
H(2,2)=cH4*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
H(2,3)=rho*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH5);
H(2,4)=rho*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH5);
H(2,5)=cH4*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
H(2,6)=rho*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH6);
H(2,7)=rho*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH6);
H(2,8)=cH4*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 = 1.0*rho/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
Kee(0,1)=cKee1;
Kee(0,2)=cKee2;
Kee(1,0)=cKee1;
Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
Kee(1,2)=cKee3;
Kee(2,0)=cKee2;
Kee(2,1)=cKee3;
Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs_ee2 = crhs_ee0 + crhs_ee1 - volume_error_ratio;
const double crhs_ee3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee5 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee6 = crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 = crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
rhs_ee[0]=-rho*(DNenr(0,0)*crhs_ee6 + DNenr(0,1)*crhs_ee7 + Nenr[0]*crhs_ee2);
rhs_ee[1]=-rho*(DNenr(1,0)*crhs_ee6 + DNenr(1,1)*crhs_ee7 + Nenr[1]*crhs_ee2);
rhs_ee[2]=-rho*(DNenr(2,0)*crhs_ee6 + DNenr(2,1)*crhs_ee7 + Nenr[2]*crhs_ee2);


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

    const double cV0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV3 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 = K_darcy*cV3;
const double cV5 = DNenr(0,0)*N[0];
const double cV6 = cV3*rho;
const double cV7 = cV6*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV8 = cV6*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV9 = N[0]*cV4;
const double cV10 = N[0]*cV7;
const double cV11 = N[1]*cV4;
const double cV12 = N[1]*cV7;
const double cV13 = cV6*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV14 = N[2]*cV4;
const double cV15 = N[2]*cV7;
const double cV16 = cV6*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV17 = N[3]*cV4;
const double cV18 = N[3]*cV7;
const double cV19 = cV6*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
V(0,0)=-DN(0,0)*Nenr[0] + DNenr(0,0)*cV8 - cV4*cV5 + cV5*cV7;
V(0,1)=-DN(0,0)*Nenr[1] + DNenr(1,0)*cV10 + DNenr(1,0)*cV8 - DNenr(1,0)*cV9;
V(0,2)=-DN(0,0)*Nenr[2] + DNenr(2,0)*cV10 + DNenr(2,0)*cV8 - DNenr(2,0)*cV9;
V(0,3)=-DN(0,0)*Nenr[3] + DNenr(3,0)*cV10 + DNenr(3,0)*cV8 - DNenr(3,0)*cV9;
V(1,0)=-DN(0,1)*Nenr[0] + DNenr(0,1)*cV10 + DNenr(0,1)*cV8 - DNenr(0,1)*cV9;
V(1,1)=-DN(0,1)*Nenr[1] + DNenr(1,1)*cV10 + DNenr(1,1)*cV8 - DNenr(1,1)*cV9;
V(1,2)=-DN(0,1)*Nenr[2] + DNenr(2,1)*cV10 + DNenr(2,1)*cV8 - DNenr(2,1)*cV9;
V(1,3)=-DN(0,1)*Nenr[3] + DNenr(3,1)*cV10 + DNenr(3,1)*cV8 - DNenr(3,1)*cV9;
V(2,0)=-DN(0,2)*Nenr[0] + DNenr(0,2)*cV10 + DNenr(0,2)*cV8 - DNenr(0,2)*cV9;
V(2,1)=-DN(0,2)*Nenr[1] + DNenr(1,2)*cV10 + DNenr(1,2)*cV8 - DNenr(1,2)*cV9;
V(2,2)=-DN(0,2)*Nenr[2] + DNenr(2,2)*cV10 + DNenr(2,2)*cV8 - DNenr(2,2)*cV9;
V(2,3)=-DN(0,2)*Nenr[3] + DNenr(3,2)*cV10 + DNenr(3,2)*cV8 - DNenr(3,2)*cV9;
V(3,0)=cV6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV6*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
V(4,0)=-DN(1,0)*Nenr[0] - DNenr(0,0)*cV11 + DNenr(0,0)*cV12 + DNenr(0,0)*cV13;
V(4,1)=-DN(1,0)*Nenr[1] - DNenr(1,0)*cV11 + DNenr(1,0)*cV12 + DNenr(1,0)*cV13;
V(4,2)=-DN(1,0)*Nenr[2] - DNenr(2,0)*cV11 + DNenr(2,0)*cV12 + DNenr(2,0)*cV13;
V(4,3)=-DN(1,0)*Nenr[3] - DNenr(3,0)*cV11 + DNenr(3,0)*cV12 + DNenr(3,0)*cV13;
V(5,0)=-DN(1,1)*Nenr[0] - DNenr(0,1)*cV11 + DNenr(0,1)*cV12 + DNenr(0,1)*cV13;
V(5,1)=-DN(1,1)*Nenr[1] - DNenr(1,1)*cV11 + DNenr(1,1)*cV12 + DNenr(1,1)*cV13;
V(5,2)=-DN(1,1)*Nenr[2] - DNenr(2,1)*cV11 + DNenr(2,1)*cV12 + DNenr(2,1)*cV13;
V(5,3)=-DN(1,1)*Nenr[3] - DNenr(3,1)*cV11 + DNenr(3,1)*cV12 + DNenr(3,1)*cV13;
V(6,0)=-DN(1,2)*Nenr[0] - DNenr(0,2)*cV11 + DNenr(0,2)*cV12 + DNenr(0,2)*cV13;
V(6,1)=-DN(1,2)*Nenr[1] - DNenr(1,2)*cV11 + DNenr(1,2)*cV12 + DNenr(1,2)*cV13;
V(6,2)=-DN(1,2)*Nenr[2] - DNenr(2,2)*cV11 + DNenr(2,2)*cV12 + DNenr(2,2)*cV13;
V(6,3)=-DN(1,2)*Nenr[3] - DNenr(3,2)*cV11 + DNenr(3,2)*cV12 + DNenr(3,2)*cV13;
V(7,0)=cV6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV6*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
V(8,0)=-DN(2,0)*Nenr[0] - DNenr(0,0)*cV14 + DNenr(0,0)*cV15 + DNenr(0,0)*cV16;
V(8,1)=-DN(2,0)*Nenr[1] - DNenr(1,0)*cV14 + DNenr(1,0)*cV15 + DNenr(1,0)*cV16;
V(8,2)=-DN(2,0)*Nenr[2] - DNenr(2,0)*cV14 + DNenr(2,0)*cV15 + DNenr(2,0)*cV16;
V(8,3)=-DN(2,0)*Nenr[3] - DNenr(3,0)*cV14 + DNenr(3,0)*cV15 + DNenr(3,0)*cV16;
V(9,0)=-DN(2,1)*Nenr[0] - DNenr(0,1)*cV14 + DNenr(0,1)*cV15 + DNenr(0,1)*cV16;
V(9,1)=-DN(2,1)*Nenr[1] - DNenr(1,1)*cV14 + DNenr(1,1)*cV15 + DNenr(1,1)*cV16;
V(9,2)=-DN(2,1)*Nenr[2] - DNenr(2,1)*cV14 + DNenr(2,1)*cV15 + DNenr(2,1)*cV16;
V(9,3)=-DN(2,1)*Nenr[3] - DNenr(3,1)*cV14 + DNenr(3,1)*cV15 + DNenr(3,1)*cV16;
V(10,0)=-DN(2,2)*Nenr[0] - DNenr(0,2)*cV14 + DNenr(0,2)*cV15 + DNenr(0,2)*cV16;
V(10,1)=-DN(2,2)*Nenr[1] - DNenr(1,2)*cV14 + DNenr(1,2)*cV15 + DNenr(1,2)*cV16;
V(10,2)=-DN(2,2)*Nenr[2] - DNenr(2,2)*cV14 + DNenr(2,2)*cV15 + DNenr(2,2)*cV16;
V(10,3)=-DN(2,2)*Nenr[3] - DNenr(3,2)*cV14 + DNenr(3,2)*cV15 + DNenr(3,2)*cV16;
V(11,0)=cV6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV6*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
V(12,0)=-DN(3,0)*Nenr[0] - DNenr(0,0)*cV17 + DNenr(0,0)*cV18 + DNenr(0,0)*cV19;
V(12,1)=-DN(3,0)*Nenr[1] - DNenr(1,0)*cV17 + DNenr(1,0)*cV18 + DNenr(1,0)*cV19;
V(12,2)=-DN(3,0)*Nenr[2] - DNenr(2,0)*cV17 + DNenr(2,0)*cV18 + DNenr(2,0)*cV19;
V(12,3)=-DN(3,0)*Nenr[3] - DNenr(3,0)*cV17 + DNenr(3,0)*cV18 + DNenr(3,0)*cV19;
V(13,0)=-DN(3,1)*Nenr[0] - DNenr(0,1)*cV17 + DNenr(0,1)*cV18 + DNenr(0,1)*cV19;
V(13,1)=-DN(3,1)*Nenr[1] - DNenr(1,1)*cV17 + DNenr(1,1)*cV18 + DNenr(1,1)*cV19;
V(13,2)=-DN(3,1)*Nenr[2] - DNenr(2,1)*cV17 + DNenr(2,1)*cV18 + DNenr(2,1)*cV19;
V(13,3)=-DN(3,1)*Nenr[3] - DNenr(3,1)*cV17 + DNenr(3,1)*cV18 + DNenr(3,1)*cV19;
V(14,0)=-DN(3,2)*Nenr[0] - DNenr(0,2)*cV17 + DNenr(0,2)*cV18 + DNenr(0,2)*cV19;
V(14,1)=-DN(3,2)*Nenr[1] - DNenr(1,2)*cV17 + DNenr(1,2)*cV18 + DNenr(1,2)*cV19;
V(14,2)=-DN(3,2)*Nenr[2] - DNenr(2,2)*cV17 + DNenr(2,2)*cV18 + DNenr(2,2)*cV19;
V(14,3)=-DN(3,2)*Nenr[3] - DNenr(3,2)*cV17 + DNenr(3,2)*cV18 + DNenr(3,2)*cV19;
V(15,0)=cV6*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
V(15,1)=cV6*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
V(15,2)=cV6*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
V(15,3)=cV6*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 = cH3*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0));
const double cH5 = cH3*rho;
const double cH6 = cH3*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0));
const double cH7 = cH3*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0));
const double cH8 = cH3*(K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0));
H(0,0)=rho*(DN(0,0)*Nenr[0] + DNenr(0,0)*cH4);
H(0,1)=rho*(DN(0,1)*Nenr[0] + DNenr(0,1)*cH4);
H(0,2)=rho*(DN(0,2)*Nenr[0] + DNenr(0,2)*cH4);
H(0,3)=cH5*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
H(0,4)=rho*(DN(1,0)*Nenr[0] + DNenr(0,0)*cH6);
H(0,5)=rho*(DN(1,1)*Nenr[0] + DNenr(0,1)*cH6);
H(0,6)=rho*(DN(1,2)*Nenr[0] + DNenr(0,2)*cH6);
H(0,7)=cH5*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
H(0,8)=rho*(DN(2,0)*Nenr[0] + DNenr(0,0)*cH7);
H(0,9)=rho*(DN(2,1)*Nenr[0] + DNenr(0,1)*cH7);
H(0,10)=rho*(DN(2,2)*Nenr[0] + DNenr(0,2)*cH7);
H(0,11)=cH5*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
H(0,12)=rho*(DN(3,0)*Nenr[0] + DNenr(0,0)*cH8);
H(0,13)=rho*(DN(3,1)*Nenr[0] + DNenr(0,1)*cH8);
H(0,14)=rho*(DN(3,2)*Nenr[0] + DNenr(0,2)*cH8);
H(0,15)=cH5*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
H(1,0)=rho*(DN(0,0)*Nenr[1] + DNenr(1,0)*cH4);
H(1,1)=rho*(DN(0,1)*Nenr[1] + DNenr(1,1)*cH4);
H(1,2)=rho*(DN(0,2)*Nenr[1] + DNenr(1,2)*cH4);
H(1,3)=cH5*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
H(1,4)=rho*(DN(1,0)*Nenr[1] + DNenr(1,0)*cH6);
H(1,5)=rho*(DN(1,1)*Nenr[1] + DNenr(1,1)*cH6);
H(1,6)=rho*(DN(1,2)*Nenr[1] + DNenr(1,2)*cH6);
H(1,7)=cH5*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
H(1,8)=rho*(DN(2,0)*Nenr[1] + DNenr(1,0)*cH7);
H(1,9)=rho*(DN(2,1)*Nenr[1] + DNenr(1,1)*cH7);
H(1,10)=rho*(DN(2,2)*Nenr[1] + DNenr(1,2)*cH7);
H(1,11)=cH5*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
H(1,12)=rho*(DN(3,0)*Nenr[1] + DNenr(1,0)*cH8);
H(1,13)=rho*(DN(3,1)*Nenr[1] + DNenr(1,1)*cH8);
H(1,14)=rho*(DN(3,2)*Nenr[1] + DNenr(1,2)*cH8);
H(1,15)=cH5*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
H(2,0)=rho*(DN(0,0)*Nenr[2] + DNenr(2,0)*cH4);
H(2,1)=rho*(DN(0,1)*Nenr[2] + DNenr(2,1)*cH4);
H(2,2)=rho*(DN(0,2)*Nenr[2] + DNenr(2,2)*cH4);
H(2,3)=cH5*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
H(2,4)=rho*(DN(1,0)*Nenr[2] + DNenr(2,0)*cH6);
H(2,5)=rho*(DN(1,1)*Nenr[2] + DNenr(2,1)*cH6);
H(2,6)=rho*(DN(1,2)*Nenr[2] + DNenr(2,2)*cH6);
H(2,7)=cH5*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
H(2,8)=rho*(DN(2,0)*Nenr[2] + DNenr(2,0)*cH7);
H(2,9)=rho*(DN(2,1)*Nenr[2] + DNenr(2,1)*cH7);
H(2,10)=rho*(DN(2,2)*Nenr[2] + DNenr(2,2)*cH7);
H(2,11)=cH5*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
H(2,12)=rho*(DN(3,0)*Nenr[2] + DNenr(2,0)*cH8);
H(2,13)=rho*(DN(3,1)*Nenr[2] + DNenr(2,1)*cH8);
H(2,14)=rho*(DN(3,2)*Nenr[2] + DNenr(2,2)*cH8);
H(2,15)=cH5*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
H(3,0)=rho*(DN(0,0)*Nenr[3] + DNenr(3,0)*cH4);
H(3,1)=rho*(DN(0,1)*Nenr[3] + DNenr(3,1)*cH4);
H(3,2)=rho*(DN(0,2)*Nenr[3] + DNenr(3,2)*cH4);
H(3,3)=cH5*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
H(3,4)=rho*(DN(1,0)*Nenr[3] + DNenr(3,0)*cH6);
H(3,5)=rho*(DN(1,1)*Nenr[3] + DNenr(3,1)*cH6);
H(3,6)=rho*(DN(1,2)*Nenr[3] + DNenr(3,2)*cH6);
H(3,7)=cH5*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
H(3,8)=rho*(DN(2,0)*Nenr[3] + DNenr(3,0)*cH7);
H(3,9)=rho*(DN(2,1)*Nenr[3] + DNenr(3,1)*cH7);
H(3,10)=rho*(DN(2,2)*Nenr[3] + DNenr(3,2)*cH7);
H(3,11)=cH5*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
H(3,12)=rho*(DN(3,0)*Nenr[3] + DNenr(3,0)*cH8);
H(3,13)=rho*(DN(3,1)*Nenr[3] + DNenr(3,1)*cH8);
H(3,14)=rho*(DN(3,2)*Nenr[3] + DNenr(3,2)*cH8);
H(3,15)=cH5*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 = 1.0*rho/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 = cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 = cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 = cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 = cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 = cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 = cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
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


    const double crhs_ee0 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs_ee1 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs_ee2 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs_ee3 = crhs_ee0 + crhs_ee1 + crhs_ee2 - volume_error_ratio;
const double crhs_ee4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee7 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee4, 2) + pow(crhs_ee5, 2) + pow(crhs_ee6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee8 = crhs_ee7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*crhs_ee4 + crhs_ee5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee9 = crhs_ee7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee1*crhs_ee5 + crhs_ee4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee10 = crhs_ee7*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee2*crhs_ee6 + crhs_ee4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
rhs_ee[0]=-rho*(DNenr(0,0)*crhs_ee8 + DNenr(0,1)*crhs_ee9 + DNenr(0,2)*crhs_ee10 + Nenr[0]*crhs_ee3);
rhs_ee[1]=-rho*(DNenr(1,0)*crhs_ee8 + DNenr(1,1)*crhs_ee9 + DNenr(1,2)*crhs_ee10 + Nenr[1]*crhs_ee3);
rhs_ee[2]=-rho*(DNenr(2,0)*crhs_ee8 + DNenr(2,1)*crhs_ee9 + DNenr(2,2)*crhs_ee10 + Nenr[2]*crhs_ee3);
rhs_ee[3]=-rho*(DNenr(3,0)*crhs_ee8 + DNenr(3,1)*crhs_ee9 + DNenr(3,2)*crhs_ee10 + Nenr[3]*crhs_ee3);


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
void TwoFluidNavierStokes<TElementData>::ComputeSplitInterface(
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
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesAlphaMethodData<2, 3> >::pGetModifiedShapeFunctionsUtility(
    const GeometryType::Pointer pGeometry,
    const Vector& rDistances)
{
    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer TwoFluidNavierStokes< TwoFluidNavierStokesAlphaMethodData<3, 4> >::pGetModifiedShapeFunctionsUtility(
        const GeometryType::Pointer pGeometry,
        const Vector& rDistances)
{
    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
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
void TwoFluidNavierStokes<TElementData>::CalculateStrainRate(TElementData& rData) const
{
    FluidElement<TElementData>::CalculateStrainRate(rData);
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
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos)
    {
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg)
    {
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double Vol = positive_volume + negative_volume;

    // We only enrich elements which are not almost empty/full
    if (positive_volume / Vol > min_area_ratio && negative_volume / Vol > min_area_ratio)
    {

        // Compute the maximum diagonal value in the enrichment stiffness matrix
        double max_diag = 0.0;
        for (unsigned int k = 0; k < NumNodes; ++k)
        {
            if (std::abs(rKeeTot(k, k)) > max_diag)
            {
                max_diag = std::abs(rKeeTot(k, k));
            }
        }
        if (max_diag == 0.0)
        {
            max_diag = 1.0;
        }
        // "weakly" impose continuity
        for (unsigned int i = 0; i < Dim; ++i)
        {
            const double di = std::abs(rData.Distance[i]);
            for (unsigned int j = i + 1; j < NumNodes; ++j)
            {
                const double dj = std::abs(rData.Distance[j]);
                // Check if the edge is cut, if it is, set the penalty constraint
                if (rData.Distance[i] * rData.Distance[j] < 0.0)
                {
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

        Vector U = ZeroVector(NumNodes * (Dim + 1));
        U[8] = 5005;
        const Vector H_U = prod(rHtot, U);
        const Vector p_enr = prod(inverse_diag, (rRHSeeTot - H_U));
        KRATOS_WATCH(p_enr);
        // // -----------------------------------------------------------
        // "Strongly impose continuity"
        //----------------------------------------------------------
        // MatrixType t_mpc_matrix;
        // CalculateConnectivityMPCsMatrix(rData, t_mpc_matrix);
        // KRATOS_WATCH(t_mpc_matrix)

        // // // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        // MatrixType t_transpose_matrix = ZeroMatrix(t_mpc_matrix.size2(), t_mpc_matrix.size1());
        // for (unsigned int i = 0; i < t_mpc_matrix.size1(); ++i)
        // {
        //     for (unsigned int j = 0; j < t_mpc_matrix.size2(); ++j)
        //     {
        //         t_transpose_matrix(j, i) = t_mpc_matrix(i, j);
        //     }
        // }

        // // ( Tt kee T)^-1

        // double det;
        // MatrixType inverse_diag;
        // const Matrix ttxkenr = prod(t_transpose_matrix, rKeeTot);
        // const Matrix ttxkenrxt = prod(ttxkenr, t_mpc_matrix);
        // MathUtils<double>::InvertMatrix(ttxkenrxt, inverse_diag, det);
        // // V T ( Tt kee T)^-1 Tt H
        // const Matrix msc = prod(rVtot, t_mpc_matrix);
        // const Matrix msc1 = prod(msc, inverse_diag);
        // const Matrix msc2 = prod(msc1, t_transpose_matrix);

        // // Adding erichment contribution to LHS
        // noalias(rLeftHandSideMatrix) -= prod(msc2, rHtot);
        // // V T ( Tt kee T)^-1 Tt Fenr
        // // Adding erichment contribution to RHS

        // const Vector msc3 = prod(msc2, rRHSeeTot);

        // noalias(rRightHandSideVector) -= msc3;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateConnectivityMPCsMatrix(
    const TElementData &rData,
    MatrixType &rConnectivityMatrix)
{
    // Determine slave/master nodes and Id all the cases:
    std::size_t pos_nodes = 0;
    std::size_t neg_nodes = 0;
    std::size_t master_no = 0;
    std::size_t slave_no = 0;
    std::vector<std::size_t> neg_nodes_id;
    std::vector<std::size_t> pos_nodes_id;
    std::vector<std::size_t> *p_master_nodes_ids = nullptr;
    std::vector<std::size_t> *p_slave_nodes_ids = nullptr;
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (rData.Distance[i] < 0.0)
        {
            neg_nodes++;
            neg_nodes_id.push_back(i);
        }
        else
        {
            pos_nodes++;
            pos_nodes_id.push_back(i);
        }
    }

    if (pos_nodes > neg_nodes)
    {
        master_no = neg_nodes;
        slave_no = pos_nodes;
        p_master_nodes_ids = &neg_nodes_id;
        p_slave_nodes_ids = &pos_nodes_id;
    }
    else
    {
        // Note this is valid for both the 2D and 3D cases
        // In the 3D case, it covers both the 1 slave and 3 masters and the 2 slave 2 master nodes
        // Also note that in the 3D 2-2 case we decide to take the negative ones as master nodes
        master_no = pos_nodes;
        slave_no = neg_nodes;
        p_master_nodes_ids = &pos_nodes_id;
        p_slave_nodes_ids = &neg_nodes_id;
    }

    // Resize T multi point contraint matrix
    rConnectivityMatrix = ZeroMatrix(NumNodes, master_no);
    const double mpc_weight = 1.0 / master_no;

    // Fill the master MPC entries
    for (std::size_t i = 0; i < p_master_nodes_ids->size(); ++i)
    {
        rConnectivityMatrix((*p_master_nodes_ids)[i], i) = 1.0;
    }

    // Fill the slave MPC entries
    for (auto slave_id : *p_slave_nodes_ids)
    {
        const double slave_distance_i = std::abs(rData.Distance[slave_id]);
        for (std::size_t j = 0; j < rConnectivityMatrix.size2(); ++j)
        {
            const std::size_t master_id = (*p_master_nodes_ids)[j];
            const double master_distance_j = std::abs(rData.Distance[master_id]);

            const double enr_sh_func_slave = master_distance_j / (slave_distance_i + master_distance_j);
            const double enr_sh_func_master = slave_distance_i / (slave_distance_i + master_distance_j);

            rConnectivityMatrix(slave_id, j) = mpc_weight * enr_sh_func_master / enr_sh_func_slave;
        }
    }
    if (master_no = slave_no)
    {

        // we need to calcuate the M matrix in order to multiply impose a only master node.
        // SUPER MASTER : p_master_nodes_ids[0]
        // AN ARBITRARY SLAVE  :  p_slave_nodes_ids[0]
        MatrixType rMasterDependencyMatrix = ZeroMatrix(master_no, 1);

        rMasterDependencyMatrix(0, 1) = rConnectivityMatrix((*p_master_nodes_ids)[0],0)/ rConnectivityMatrix((*p_slave_nodes_ids)[0], 1);
        rMasterDependencyMatrix(1, 1) = 1.0;

        rConnectivityMatrix = prod(rConnectivityMatrix, rMasterDependencyMatrix);
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
    MatrixType& rInterfaceShapeFunction,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivatives,
    Vector& rInterfaceWeights,
    std::vector< array_1d<double, 3> >& rInterfaceNormalsNeg,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    VectorType &rRHSeeTot)
{
    // Surface tension coefficient is set in material properties
    const double surface_tension_coefficient = this->GetProperties().GetValue(SURFACE_TENSION_COEFFICIENT);

    Vector gauss_pts_curvature; // curvatures calculated on interface Gauss points

    CalculateCurvatureOnInterfaceGaussPoints(
        rInterfaceShapeFunction,
        gauss_pts_curvature);

    SurfaceTension(
        surface_tension_coefficient,
        gauss_pts_curvature,
        rInterfaceWeights,
        rInterfaceShapeFunction,
        rInterfaceNormalsNeg,
        rRightHandSideVector);

    this->PressureGradientStabilization(
        rData,
        rInterfaceWeights,
        rEnrInterfaceShapeFunctionPos,
        rEnrInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivatives,
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

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointLHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    VectorType &rRHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointRHSContribution(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    VectorType &RLHS)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<2, 3> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>::ComputeGaussPointEnrichmentContributions(
    TwoFluidNavierStokesAlphaMethodData<3, 4> &rData,
    MatrixType &rV,
    MatrixType &rH,
    MatrixType &rKee,
    VectorType &rRHS_ee)
{
    KRATOS_ERROR << "This function should never be called. It is only to enable the explicit template instantiation." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>;

template class TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesAlphaMethodData<3, 4>>;

} // namespace Kratos
