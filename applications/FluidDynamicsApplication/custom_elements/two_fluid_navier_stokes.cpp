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
    const auto vn = rData.Velocity_OldStep1;
    // const auto vconv = rData.Velocity - rData.Velocity_OldStep1;
    const auto vconv = rData.Velocity - rData.MeshVelocity;

    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;
    // const auto an = rData.Acceleration;

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
const double clhs7 = clhs6*h/stab_c1 + mu;
const double clhs8 = pow(N[0], 2);
const double clhs9 = rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs10 = bdf0*rho;
const double clhs11 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double clhs12 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double clhs13 = rho*(DN(0,0)*clhs11 + DN(0,1)*clhs12);
const double clhs14 = K_darcy*N[0];
const double clhs15 = N[0]*clhs10;
const double clhs16 = -clhs13 + clhs14 + clhs15 + clhs9;
const double clhs17 = 1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = 1.0*clhs9;
const double clhs19 = clhs17*clhs18;
const double clhs20 = 1.0*clhs14;
const double clhs21 = clhs17*clhs20;
const double clhs22 = 1.0*clhs17;
const double clhs23 = clhs16*clhs22;
const double clhs24 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double clhs25 = N[0]*clhs24;
const double clhs26 = K_darcy*clhs8 - N[0]*clhs13 + N[0]*clhs9 + clhs10*clhs8 + clhs16*clhs19 - clhs16*clhs21 + clhs23*clhs25;
const double clhs27 = C(0,1)*DN(0,1) + clhs1;
const double clhs28 = C(1,2)*DN(0,1);
const double clhs29 = C(2,2)*DN(0,0) + clhs28;
const double clhs30 = DN(0,0)*clhs7;
const double clhs31 = DN(0,1)*clhs30;
const double clhs32 = clhs22*clhs24;
const double clhs33 = N[0]*clhs32 - N[0] + clhs19 - clhs21;
const double clhs34 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs35 = C(0,2)*DN(1,0);
const double clhs36 = C(2,2)*DN(1,1) + clhs35;
const double clhs37 = DN(0,0)*DN(1,0);
const double clhs38 = N[1]*clhs14 + N[1]*clhs15;
const double clhs39 = clhs37*clhs7 + clhs38;
const double clhs40 = rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs41 = rho*(DN(1,0)*clhs11 + DN(1,1)*clhs12);
const double clhs42 = K_darcy*N[1];
const double clhs43 = N[1]*clhs10;
const double clhs44 = clhs40 - clhs41 + clhs42 + clhs43;
const double clhs45 = clhs22*clhs44;
const double clhs46 = N[0]*clhs40 - N[0]*clhs41 + clhs19*clhs44 - clhs21*clhs44 + clhs25*clhs45;
const double clhs47 = C(0,1)*DN(1,1) + clhs35;
const double clhs48 = C(1,2)*DN(1,1);
const double clhs49 = C(2,2)*DN(1,0) + clhs48;
const double clhs50 = DN(1,1)*clhs30;
const double clhs51 = DN(0,0)*N[1];
const double clhs52 = DN(1,0)*N[0];
const double clhs53 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs54 = C(0,2)*DN(2,0);
const double clhs55 = C(2,2)*DN(2,1) + clhs54;
const double clhs56 = DN(0,0)*DN(2,0);
const double clhs57 = N[2]*clhs14 + N[2]*clhs15;
const double clhs58 = clhs56*clhs7 + clhs57;
const double clhs59 = rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs60 = rho*(DN(2,0)*clhs11 + DN(2,1)*clhs12);
const double clhs61 = K_darcy*N[2];
const double clhs62 = N[2]*clhs10;
const double clhs63 = clhs59 - clhs60 + clhs61 + clhs62;
const double clhs64 = clhs22*clhs63;
const double clhs65 = N[0]*clhs59 - N[0]*clhs60 + clhs19*clhs63 - clhs21*clhs63 + clhs25*clhs64;
const double clhs66 = C(0,1)*DN(2,1) + clhs54;
const double clhs67 = C(1,2)*DN(2,1);
const double clhs68 = C(2,2)*DN(2,0) + clhs67;
const double clhs69 = DN(2,1)*clhs30;
const double clhs70 = DN(0,0)*N[2];
const double clhs71 = DN(2,0)*N[0];
const double clhs72 = C(0,1)*DN(0,0) + clhs28;
const double clhs73 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs74 = pow(DN(0,1), 2);
const double clhs75 = C(0,1)*DN(1,0) + clhs48;
const double clhs76 = DN(0,1)*clhs7;
const double clhs77 = DN(1,0)*clhs76;
const double clhs78 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs79 = DN(0,1)*DN(1,1);
const double clhs80 = clhs38 + clhs7*clhs79;
const double clhs81 = DN(0,1)*N[1];
const double clhs82 = DN(1,1)*N[0];
const double clhs83 = C(0,1)*DN(2,0) + clhs67;
const double clhs84 = DN(2,0)*clhs76;
const double clhs85 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs86 = DN(0,1)*DN(2,1);
const double clhs87 = clhs57 + clhs7*clhs86;
const double clhs88 = DN(0,1)*N[2];
const double clhs89 = DN(2,1)*N[0];
const double clhs90 = N[0] + clhs17*(-1.0*clhs13 + 1.0*clhs15 + clhs18 + clhs20);
const double clhs91 = clhs22*(clhs37 + clhs79);
const double clhs92 = clhs22*(clhs56 + clhs86);
const double clhs93 = clhs22*clhs40;
const double clhs94 = clhs22*clhs42;
const double clhs95 = N[1]*clhs24;
const double clhs96 = -N[1]*clhs13 + N[1]*clhs9 + clhs16*clhs93 - clhs16*clhs94 + clhs23*clhs95;
const double clhs97 = pow(DN(1,0), 2);
const double clhs98 = pow(N[1], 2);
const double clhs99 = K_darcy*clhs98 + N[1]*clhs40 - N[1]*clhs41 + clhs10*clhs98 + clhs40*clhs45 - clhs42*clhs45 + clhs45*clhs95;
const double clhs100 = DN(1,0)*clhs7;
const double clhs101 = DN(1,1)*clhs100;
const double clhs102 = N[1]*clhs32 - N[1] + clhs93 - clhs94;
const double clhs103 = DN(1,0)*DN(2,0);
const double clhs104 = N[2]*clhs42 + N[2]*clhs43;
const double clhs105 = clhs103*clhs7 + clhs104;
const double clhs106 = N[1]*clhs59 - N[1]*clhs60 + clhs40*clhs64 - clhs42*clhs64 + clhs64*clhs95;
const double clhs107 = DN(2,1)*clhs100;
const double clhs108 = DN(1,0)*N[2];
const double clhs109 = DN(2,0)*N[1];
const double clhs110 = pow(DN(1,1), 2);
const double clhs111 = DN(2,0)*clhs7;
const double clhs112 = DN(1,1)*clhs111;
const double clhs113 = DN(1,1)*DN(2,1);
const double clhs114 = clhs104 + clhs113*clhs7;
const double clhs115 = DN(1,1)*N[2];
const double clhs116 = DN(2,1)*N[1];
const double clhs117 = N[1] + clhs17*(1.0*clhs40 - 1.0*clhs41 + 1.0*clhs42 + 1.0*clhs43);
const double clhs118 = clhs22*(clhs103 + clhs113);
const double clhs119 = N[2]*clhs24;
const double clhs120 = -N[2]*clhs13 + N[2]*clhs9 + clhs119*clhs23 + clhs23*clhs59 - clhs23*clhs61;
const double clhs121 = clhs22*clhs61;
const double clhs122 = clhs22*clhs59;
const double clhs123 = N[2]*clhs40 - N[2]*clhs41 + clhs119*clhs45 + clhs45*clhs59 - clhs45*clhs61;
const double clhs124 = pow(DN(2,0), 2);
const double clhs125 = pow(N[2], 2);
const double clhs126 = K_darcy*clhs125 + N[2]*clhs59 - N[2]*clhs60 + clhs10*clhs125 + clhs119*clhs64 + clhs59*clhs64 - clhs61*clhs64;
const double clhs127 = DN(2,1)*clhs111;
const double clhs128 = N[2]*clhs32 - N[2] - clhs121 + clhs122;
const double clhs129 = pow(DN(2,1), 2);
const double clhs130 = N[2] + clhs17*(1.0*clhs59 - 1.0*clhs60 + 1.0*clhs61 + 1.0*clhs62);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs26 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + clhs31;
lhs(0,2)=DN(0,0)*clhs33;
lhs(0,3)=DN(0,0)*clhs34 + DN(0,1)*clhs36 + clhs39 + clhs46;
lhs(0,4)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + clhs50;
lhs(0,5)=DN(1,0)*clhs19 - DN(1,0)*clhs21 + clhs32*clhs52 - clhs51;
lhs(0,6)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + clhs58 + clhs65;
lhs(0,7)=DN(0,0)*clhs66 + DN(0,1)*clhs68 + clhs69;
lhs(0,8)=DN(2,0)*clhs19 - DN(2,0)*clhs21 + clhs32*clhs71 - clhs70;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs72 + clhs31;
lhs(1,1)=DN(0,0)*clhs29 + DN(0,1)*clhs73 + clhs26 + clhs7*clhs74;
lhs(1,2)=DN(0,1)*clhs33;
lhs(1,3)=DN(0,0)*clhs36 + DN(0,1)*clhs75 + clhs77;
lhs(1,4)=DN(0,0)*clhs49 + DN(0,1)*clhs78 + clhs46 + clhs80;
lhs(1,5)=DN(1,1)*clhs19 - DN(1,1)*clhs21 + clhs32*clhs82 - clhs81;
lhs(1,6)=DN(0,0)*clhs55 + DN(0,1)*clhs83 + clhs84;
lhs(1,7)=DN(0,0)*clhs68 + DN(0,1)*clhs85 + clhs65 + clhs87;
lhs(1,8)=DN(2,1)*clhs19 - DN(2,1)*clhs21 + clhs32*clhs89 - clhs88;
lhs(2,0)=DN(0,0)*clhs90;
lhs(2,1)=DN(0,1)*clhs90;
lhs(2,2)=clhs22*(clhs3 + clhs74);
lhs(2,3)=DN(0,0)*clhs45 + clhs52;
lhs(2,4)=DN(0,1)*clhs45 + clhs82;
lhs(2,5)=clhs91;
lhs(2,6)=DN(0,0)*clhs64 + clhs71;
lhs(2,7)=DN(0,1)*clhs64 + clhs89;
lhs(2,8)=clhs92;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs39 + clhs96;
lhs(3,1)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs77;
lhs(3,2)=DN(0,0)*clhs93 - DN(0,0)*clhs94 + clhs32*clhs51 - clhs52;
lhs(3,3)=DN(1,0)*clhs34 + DN(1,1)*clhs36 + clhs7*clhs97 + clhs99;
lhs(3,4)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + clhs101;
lhs(3,5)=DN(1,0)*clhs102;
lhs(3,6)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + clhs105 + clhs106;
lhs(3,7)=DN(1,0)*clhs66 + DN(1,1)*clhs68 + clhs107;
lhs(3,8)=DN(2,0)*clhs93 - DN(2,0)*clhs94 - clhs108 + clhs109*clhs32;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs72 + clhs50;
lhs(4,1)=DN(1,0)*clhs29 + DN(1,1)*clhs73 + clhs80 + clhs96;
lhs(4,2)=DN(0,1)*clhs93 - DN(0,1)*clhs94 + clhs32*clhs81 - clhs82;
lhs(4,3)=DN(1,0)*clhs36 + DN(1,1)*clhs75 + clhs101;
lhs(4,4)=DN(1,0)*clhs49 + DN(1,1)*clhs78 + clhs110*clhs7 + clhs99;
lhs(4,5)=DN(1,1)*clhs102;
lhs(4,6)=DN(1,0)*clhs55 + DN(1,1)*clhs83 + clhs112;
lhs(4,7)=DN(1,0)*clhs68 + DN(1,1)*clhs85 + clhs106 + clhs114;
lhs(4,8)=DN(2,1)*clhs93 - DN(2,1)*clhs94 - clhs115 + clhs116*clhs32;
lhs(5,0)=DN(1,0)*clhs23 + clhs51;
lhs(5,1)=DN(1,1)*clhs23 + clhs81;
lhs(5,2)=clhs91;
lhs(5,3)=DN(1,0)*clhs117;
lhs(5,4)=DN(1,1)*clhs117;
lhs(5,5)=clhs22*(clhs110 + clhs97);
lhs(5,6)=DN(1,0)*clhs64 + clhs109;
lhs(5,7)=DN(1,1)*clhs64 + clhs116;
lhs(5,8)=clhs118;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs120 + clhs58;
lhs(6,1)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs84;
lhs(6,2)=-DN(0,0)*clhs121 + DN(0,0)*clhs122 + clhs32*clhs70 - clhs71;
lhs(6,3)=DN(2,0)*clhs34 + DN(2,1)*clhs36 + clhs105 + clhs123;
lhs(6,4)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + clhs112;
lhs(6,5)=-DN(1,0)*clhs121 + DN(1,0)*clhs122 + clhs108*clhs32 - clhs109;
lhs(6,6)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + clhs124*clhs7 + clhs126;
lhs(6,7)=DN(2,0)*clhs66 + DN(2,1)*clhs68 + clhs127;
lhs(6,8)=DN(2,0)*clhs128;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs72 + clhs69;
lhs(7,1)=DN(2,0)*clhs29 + DN(2,1)*clhs73 + clhs120 + clhs87;
lhs(7,2)=-DN(0,1)*clhs121 + DN(0,1)*clhs122 + clhs32*clhs88 - clhs89;
lhs(7,3)=DN(2,0)*clhs36 + DN(2,1)*clhs75 + clhs107;
lhs(7,4)=DN(2,0)*clhs49 + DN(2,1)*clhs78 + clhs114 + clhs123;
lhs(7,5)=-DN(1,1)*clhs121 + DN(1,1)*clhs122 + clhs115*clhs32 - clhs116;
lhs(7,6)=DN(2,0)*clhs55 + DN(2,1)*clhs83 + clhs127;
lhs(7,7)=DN(2,0)*clhs68 + DN(2,1)*clhs85 + clhs126 + clhs129*clhs7;
lhs(7,8)=DN(2,1)*clhs128;
lhs(8,0)=DN(2,0)*clhs23 + clhs70;
lhs(8,1)=DN(2,1)*clhs23 + clhs88;
lhs(8,2)=clhs92;
lhs(8,3)=DN(2,0)*clhs45 + clhs108;
lhs(8,4)=DN(2,1)*clhs45 + clhs115;
lhs(8,5)=clhs118;
lhs(8,6)=DN(2,0)*clhs130;
lhs(8,7)=DN(2,1)*clhs130;
lhs(8,8)=clhs22*(clhs124 + clhs129);


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

    // const auto vconv = rData.Velocity - rData.Velocity_OldStep1;
    const auto vn=rData.Velocity_OldStep1;
    const auto vconv = rData.Velocity - rData.MeshVelocity;
    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;

    // const auto an = rData.Acceleration;

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
const double clhs10 = clhs9*h/stab_c1 + mu;
const double clhs11 = pow(N[0], 2);
const double clhs12 = rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs13 = bdf0*rho;
const double clhs14 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double clhs15 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double clhs16 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double clhs17 = rho*(DN(0,0)*clhs14 + DN(0,1)*clhs15 + DN(0,2)*clhs16);
const double clhs18 = K_darcy*N[0];
const double clhs19 = N[0]*clhs13;
const double clhs20 = clhs12 - clhs17 + clhs18 + clhs19;
const double clhs21 = 1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs22 = 1.0*clhs12;
const double clhs23 = clhs21*clhs22;
const double clhs24 = 1.0*clhs18;
const double clhs25 = clhs21*clhs24;
const double clhs26 = 1.0*clhs21;
const double clhs27 = clhs20*clhs26;
const double clhs28 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double clhs29 = N[0]*clhs28;
const double clhs30 = K_darcy*clhs11 + N[0]*clhs12 - N[0]*clhs17 + clhs11*clhs13 + clhs20*clhs23 - clhs20*clhs25 + clhs27*clhs29;
const double clhs31 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs32 = C(1,3)*DN(0,1);
const double clhs33 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs32;
const double clhs34 = C(3,5)*DN(0,0);
const double clhs35 = C(4,5)*DN(0,2);
const double clhs36 = C(1,5)*DN(0,1) + clhs34 + clhs35;
const double clhs37 = DN(0,0)*clhs10;
const double clhs38 = DN(0,1)*clhs37;
const double clhs39 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs40 = C(3,4)*DN(0,1);
const double clhs41 = C(2,3)*DN(0,2) + clhs34 + clhs40;
const double clhs42 = C(2,5)*DN(0,2);
const double clhs43 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs42;
const double clhs44 = DN(0,2)*clhs37;
const double clhs45 = clhs26*clhs28;
const double clhs46 = N[0]*clhs45 - N[0] + clhs23 - clhs25;
const double clhs47 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs48 = C(0,3)*DN(1,0);
const double clhs49 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs48;
const double clhs50 = C(0,5)*DN(1,0);
const double clhs51 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs50;
const double clhs52 = DN(0,0)*DN(1,0);
const double clhs53 = N[1]*clhs18 + N[1]*clhs19;
const double clhs54 = clhs10*clhs52 + clhs53;
const double clhs55 = rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs56 = rho*(DN(1,0)*clhs14 + DN(1,1)*clhs15 + DN(1,2)*clhs16);
const double clhs57 = K_darcy*N[1];
const double clhs58 = N[1]*clhs13;
const double clhs59 = clhs55 - clhs56 + clhs57 + clhs58;
const double clhs60 = clhs26*clhs59;
const double clhs61 = N[0]*clhs55 - N[0]*clhs56 + clhs23*clhs59 - clhs25*clhs59 + clhs29*clhs60;
const double clhs62 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs48;
const double clhs63 = C(1,3)*DN(1,1);
const double clhs64 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs63;
const double clhs65 = C(3,5)*DN(1,0);
const double clhs66 = C(4,5)*DN(1,2);
const double clhs67 = C(1,5)*DN(1,1) + clhs65 + clhs66;
const double clhs68 = DN(1,1)*clhs37;
const double clhs69 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs50;
const double clhs70 = C(3,4)*DN(1,1);
const double clhs71 = C(2,3)*DN(1,2) + clhs65 + clhs70;
const double clhs72 = C(2,5)*DN(1,2);
const double clhs73 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs72;
const double clhs74 = DN(1,2)*clhs37;
const double clhs75 = DN(0,0)*N[1];
const double clhs76 = DN(1,0)*N[0];
const double clhs77 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs78 = C(0,3)*DN(2,0);
const double clhs79 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs78;
const double clhs80 = C(0,5)*DN(2,0);
const double clhs81 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs80;
const double clhs82 = DN(0,0)*DN(2,0);
const double clhs83 = N[2]*clhs18 + N[2]*clhs19;
const double clhs84 = clhs10*clhs82 + clhs83;
const double clhs85 = rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs86 = rho*(DN(2,0)*clhs14 + DN(2,1)*clhs15 + DN(2,2)*clhs16);
const double clhs87 = K_darcy*N[2];
const double clhs88 = N[2]*clhs13;
const double clhs89 = clhs85 - clhs86 + clhs87 + clhs88;
const double clhs90 = clhs26*clhs89;
const double clhs91 = N[0]*clhs85 - N[0]*clhs86 + clhs23*clhs89 - clhs25*clhs89 + clhs29*clhs90;
const double clhs92 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs78;
const double clhs93 = C(1,3)*DN(2,1);
const double clhs94 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs93;
const double clhs95 = C(3,5)*DN(2,0);
const double clhs96 = C(4,5)*DN(2,2);
const double clhs97 = C(1,5)*DN(2,1) + clhs95 + clhs96;
const double clhs98 = DN(2,1)*clhs37;
const double clhs99 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs80;
const double clhs100 = C(3,4)*DN(2,1);
const double clhs101 = C(2,3)*DN(2,2) + clhs100 + clhs95;
const double clhs102 = C(2,5)*DN(2,2);
const double clhs103 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs102;
const double clhs104 = DN(2,2)*clhs37;
const double clhs105 = DN(0,0)*N[2];
const double clhs106 = DN(2,0)*N[0];
const double clhs107 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs108 = C(0,3)*DN(3,0);
const double clhs109 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs108;
const double clhs110 = C(0,5)*DN(3,0);
const double clhs111 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs110;
const double clhs112 = DN(0,0)*DN(3,0);
const double clhs113 = N[3]*clhs18 + N[3]*clhs19;
const double clhs114 = clhs10*clhs112 + clhs113;
const double clhs115 = rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs116 = rho*(DN(3,0)*clhs14 + DN(3,1)*clhs15 + DN(3,2)*clhs16);
const double clhs117 = K_darcy*N[3];
const double clhs118 = N[3]*clhs13;
const double clhs119 = clhs115 - clhs116 + clhs117 + clhs118;
const double clhs120 = clhs119*clhs26;
const double clhs121 = N[0]*clhs115 - N[0]*clhs116 + clhs119*clhs23 - clhs119*clhs25 + clhs120*clhs29;
const double clhs122 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs108;
const double clhs123 = C(1,3)*DN(3,1);
const double clhs124 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs123;
const double clhs125 = C(3,5)*DN(3,0);
const double clhs126 = C(4,5)*DN(3,2);
const double clhs127 = C(1,5)*DN(3,1) + clhs125 + clhs126;
const double clhs128 = DN(3,1)*clhs37;
const double clhs129 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs110;
const double clhs130 = C(3,4)*DN(3,1);
const double clhs131 = C(2,3)*DN(3,2) + clhs125 + clhs130;
const double clhs132 = C(2,5)*DN(3,2);
const double clhs133 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs132;
const double clhs134 = DN(3,2)*clhs37;
const double clhs135 = DN(0,0)*N[3];
const double clhs136 = DN(3,0)*N[0];
const double clhs137 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs32;
const double clhs138 = C(0,4)*DN(0,0) + clhs35 + clhs40;
const double clhs139 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs140 = C(1,4)*DN(0,1);
const double clhs141 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs140;
const double clhs142 = pow(DN(0,1), 2);
const double clhs143 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs140;
const double clhs144 = C(2,4)*DN(0,2);
const double clhs145 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs144;
const double clhs146 = DN(0,1)*clhs10;
const double clhs147 = DN(0,2)*clhs146;
const double clhs148 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs63;
const double clhs149 = C(0,4)*DN(1,0) + clhs66 + clhs70;
const double clhs150 = DN(1,0)*clhs146;
const double clhs151 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs152 = C(1,4)*DN(1,1);
const double clhs153 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs152;
const double clhs154 = DN(0,1)*DN(1,1);
const double clhs155 = clhs10*clhs154;
const double clhs156 = clhs53 + clhs61;
const double clhs157 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs152;
const double clhs158 = C(2,4)*DN(1,2);
const double clhs159 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs158;
const double clhs160 = DN(1,2)*clhs146;
const double clhs161 = DN(0,1)*N[1];
const double clhs162 = DN(1,1)*N[0];
const double clhs163 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs93;
const double clhs164 = C(0,4)*DN(2,0) + clhs100 + clhs96;
const double clhs165 = DN(2,0)*clhs146;
const double clhs166 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs167 = C(1,4)*DN(2,1);
const double clhs168 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs167;
const double clhs169 = DN(0,1)*DN(2,1);
const double clhs170 = clhs10*clhs169;
const double clhs171 = clhs83 + clhs91;
const double clhs172 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs167;
const double clhs173 = C(2,4)*DN(2,2);
const double clhs174 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs173;
const double clhs175 = DN(2,2)*clhs146;
const double clhs176 = DN(0,1)*N[2];
const double clhs177 = DN(2,1)*N[0];
const double clhs178 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs123;
const double clhs179 = C(0,4)*DN(3,0) + clhs126 + clhs130;
const double clhs180 = DN(3,0)*clhs146;
const double clhs181 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs182 = C(1,4)*DN(3,1);
const double clhs183 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs182;
const double clhs184 = DN(0,1)*DN(3,1);
const double clhs185 = clhs10*clhs184;
const double clhs186 = clhs113 + clhs121;
const double clhs187 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs182;
const double clhs188 = C(2,4)*DN(3,2);
const double clhs189 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs188;
const double clhs190 = DN(3,2)*clhs146;
const double clhs191 = DN(0,1)*N[3];
const double clhs192 = DN(3,1)*N[0];
const double clhs193 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs42;
const double clhs194 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs144;
const double clhs195 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs196 = pow(DN(0,2), 2);
const double clhs197 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs72;
const double clhs198 = DN(0,2)*clhs10;
const double clhs199 = DN(1,0)*clhs198;
const double clhs200 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs158;
const double clhs201 = DN(1,1)*clhs198;
const double clhs202 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs203 = DN(0,2)*DN(1,2);
const double clhs204 = clhs10*clhs203;
const double clhs205 = DN(0,2)*N[1];
const double clhs206 = DN(1,2)*N[0];
const double clhs207 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs102;
const double clhs208 = DN(2,0)*clhs198;
const double clhs209 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs173;
const double clhs210 = DN(2,1)*clhs198;
const double clhs211 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs212 = DN(0,2)*DN(2,2);
const double clhs213 = clhs10*clhs212;
const double clhs214 = DN(0,2)*N[2];
const double clhs215 = DN(2,2)*N[0];
const double clhs216 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs132;
const double clhs217 = DN(3,0)*clhs198;
const double clhs218 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs188;
const double clhs219 = DN(3,1)*clhs198;
const double clhs220 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs221 = DN(0,2)*DN(3,2);
const double clhs222 = clhs10*clhs221;
const double clhs223 = DN(0,2)*N[3];
const double clhs224 = DN(3,2)*N[0];
const double clhs225 = N[0] + clhs21*(-1.0*clhs17 + 1.0*clhs19 + clhs22 + clhs24);
const double clhs226 = clhs26*(clhs154 + clhs203 + clhs52);
const double clhs227 = clhs26*(clhs169 + clhs212 + clhs82);
const double clhs228 = clhs26*(clhs112 + clhs184 + clhs221);
const double clhs229 = clhs26*clhs55;
const double clhs230 = clhs26*clhs57;
const double clhs231 = N[1]*clhs28;
const double clhs232 = N[1]*clhs12 - N[1]*clhs17 + clhs20*clhs229 - clhs20*clhs230 + clhs231*clhs27;
const double clhs233 = pow(DN(1,0), 2);
const double clhs234 = pow(N[1], 2);
const double clhs235 = K_darcy*clhs234 + N[1]*clhs55 - N[1]*clhs56 + clhs13*clhs234 + clhs231*clhs60 + clhs55*clhs60 - clhs57*clhs60;
const double clhs236 = DN(1,0)*clhs10;
const double clhs237 = DN(1,1)*clhs236;
const double clhs238 = DN(1,2)*clhs236;
const double clhs239 = N[1]*clhs45 - N[1] + clhs229 - clhs230;
const double clhs240 = DN(1,0)*DN(2,0);
const double clhs241 = N[2]*clhs57 + N[2]*clhs58;
const double clhs242 = clhs10*clhs240 + clhs241;
const double clhs243 = N[1]*clhs85 - N[1]*clhs86 + clhs231*clhs90 + clhs55*clhs90 - clhs57*clhs90;
const double clhs244 = DN(2,1)*clhs236;
const double clhs245 = DN(2,2)*clhs236;
const double clhs246 = DN(1,0)*N[2];
const double clhs247 = DN(2,0)*N[1];
const double clhs248 = DN(1,0)*DN(3,0);
const double clhs249 = N[3]*clhs57 + N[3]*clhs58;
const double clhs250 = clhs10*clhs248 + clhs249;
const double clhs251 = N[1]*clhs115 - N[1]*clhs116 + clhs120*clhs231 + clhs120*clhs55 - clhs120*clhs57;
const double clhs252 = DN(3,1)*clhs236;
const double clhs253 = DN(3,2)*clhs236;
const double clhs254 = DN(1,0)*N[3];
const double clhs255 = DN(3,0)*N[1];
const double clhs256 = clhs232 + clhs53;
const double clhs257 = pow(DN(1,1), 2);
const double clhs258 = DN(1,1)*clhs10;
const double clhs259 = DN(1,2)*clhs258;
const double clhs260 = DN(2,0)*clhs258;
const double clhs261 = DN(1,1)*DN(2,1);
const double clhs262 = clhs10*clhs261;
const double clhs263 = clhs241 + clhs243;
const double clhs264 = DN(2,2)*clhs258;
const double clhs265 = DN(1,1)*N[2];
const double clhs266 = DN(2,1)*N[1];
const double clhs267 = DN(3,0)*clhs258;
const double clhs268 = DN(1,1)*DN(3,1);
const double clhs269 = clhs10*clhs268;
const double clhs270 = clhs249 + clhs251;
const double clhs271 = DN(3,2)*clhs258;
const double clhs272 = DN(1,1)*N[3];
const double clhs273 = DN(3,1)*N[1];
const double clhs274 = pow(DN(1,2), 2);
const double clhs275 = DN(1,2)*clhs10;
const double clhs276 = DN(2,0)*clhs275;
const double clhs277 = DN(2,1)*clhs275;
const double clhs278 = DN(1,2)*DN(2,2);
const double clhs279 = clhs10*clhs278;
const double clhs280 = DN(1,2)*N[2];
const double clhs281 = DN(2,2)*N[1];
const double clhs282 = DN(3,0)*clhs275;
const double clhs283 = DN(3,1)*clhs275;
const double clhs284 = DN(1,2)*DN(3,2);
const double clhs285 = clhs10*clhs284;
const double clhs286 = DN(1,2)*N[3];
const double clhs287 = DN(3,2)*N[1];
const double clhs288 = N[1] + clhs21*(1.0*clhs55 - 1.0*clhs56 + 1.0*clhs57 + 1.0*clhs58);
const double clhs289 = clhs26*(clhs240 + clhs261 + clhs278);
const double clhs290 = clhs26*(clhs248 + clhs268 + clhs284);
const double clhs291 = N[2]*clhs28;
const double clhs292 = N[2]*clhs12 - N[2]*clhs17 + clhs27*clhs291 + clhs27*clhs85 - clhs27*clhs87;
const double clhs293 = clhs26*clhs87;
const double clhs294 = clhs26*clhs85;
const double clhs295 = N[2]*clhs55 - N[2]*clhs56 + clhs291*clhs60 + clhs60*clhs85 - clhs60*clhs87;
const double clhs296 = pow(DN(2,0), 2);
const double clhs297 = pow(N[2], 2);
const double clhs298 = K_darcy*clhs297 + N[2]*clhs85 - N[2]*clhs86 + clhs13*clhs297 + clhs291*clhs90 + clhs85*clhs90 - clhs87*clhs90;
const double clhs299 = DN(2,0)*clhs10;
const double clhs300 = DN(2,1)*clhs299;
const double clhs301 = DN(2,2)*clhs299;
const double clhs302 = N[2]*clhs45 - N[2] - clhs293 + clhs294;
const double clhs303 = DN(2,0)*DN(3,0);
const double clhs304 = N[3]*clhs87 + N[3]*clhs88;
const double clhs305 = clhs10*clhs303 + clhs304;
const double clhs306 = N[2]*clhs115 - N[2]*clhs116 + clhs120*clhs291 + clhs120*clhs85 - clhs120*clhs87;
const double clhs307 = DN(3,1)*clhs299;
const double clhs308 = DN(3,2)*clhs299;
const double clhs309 = DN(2,0)*N[3];
const double clhs310 = DN(3,0)*N[2];
const double clhs311 = clhs292 + clhs83;
const double clhs312 = clhs241 + clhs295;
const double clhs313 = pow(DN(2,1), 2);
const double clhs314 = DN(2,1)*clhs10;
const double clhs315 = DN(2,2)*clhs314;
const double clhs316 = DN(3,0)*clhs314;
const double clhs317 = DN(2,1)*DN(3,1);
const double clhs318 = clhs10*clhs317;
const double clhs319 = clhs304 + clhs306;
const double clhs320 = DN(3,2)*clhs314;
const double clhs321 = DN(2,1)*N[3];
const double clhs322 = DN(3,1)*N[2];
const double clhs323 = pow(DN(2,2), 2);
const double clhs324 = DN(2,2)*clhs10;
const double clhs325 = DN(3,0)*clhs324;
const double clhs326 = DN(3,1)*clhs324;
const double clhs327 = DN(2,2)*DN(3,2);
const double clhs328 = clhs10*clhs327;
const double clhs329 = DN(2,2)*N[3];
const double clhs330 = DN(3,2)*N[2];
const double clhs331 = N[2] + clhs21*(1.0*clhs85 - 1.0*clhs86 + 1.0*clhs87 + 1.0*clhs88);
const double clhs332 = clhs26*(clhs303 + clhs317 + clhs327);
const double clhs333 = N[3]*clhs28;
const double clhs334 = N[3]*clhs12 - N[3]*clhs17 + clhs115*clhs27 - clhs117*clhs27 + clhs27*clhs333;
const double clhs335 = clhs117*clhs26;
const double clhs336 = clhs115*clhs26;
const double clhs337 = N[3]*clhs55 - N[3]*clhs56 + clhs115*clhs60 - clhs117*clhs60 + clhs333*clhs60;
const double clhs338 = N[3]*clhs85 - N[3]*clhs86 + clhs115*clhs90 - clhs117*clhs90 + clhs333*clhs90;
const double clhs339 = pow(DN(3,0), 2);
const double clhs340 = pow(N[3], 2);
const double clhs341 = K_darcy*clhs340 + N[3]*clhs115 - N[3]*clhs116 + clhs115*clhs120 - clhs117*clhs120 + clhs120*clhs333 + clhs13*clhs340;
const double clhs342 = DN(3,0)*clhs10;
const double clhs343 = DN(3,1)*clhs342;
const double clhs344 = DN(3,2)*clhs342;
const double clhs345 = N[3]*clhs45 - N[3] - clhs335 + clhs336;
const double clhs346 = clhs113 + clhs334;
const double clhs347 = clhs249 + clhs337;
const double clhs348 = clhs304 + clhs338;
const double clhs349 = pow(DN(3,1), 2);
const double clhs350 = DN(3,1)*DN(3,2)*clhs10;
const double clhs351 = pow(DN(3,2), 2);
const double clhs352 = N[3] + clhs21*(1.0*clhs115 - 1.0*clhs116 + 1.0*clhs117 + 1.0*clhs118);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs30;
lhs(0,1)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs36 + clhs38;
lhs(0,2)=DN(0,0)*clhs39 + DN(0,1)*clhs41 + DN(0,2)*clhs43 + clhs44;
lhs(0,3)=DN(0,0)*clhs46;
lhs(0,4)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + DN(0,2)*clhs51 + clhs54 + clhs61;
lhs(0,5)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs67 + clhs68;
lhs(0,6)=DN(0,0)*clhs69 + DN(0,1)*clhs71 + DN(0,2)*clhs73 + clhs74;
lhs(0,7)=DN(1,0)*clhs23 - DN(1,0)*clhs25 + clhs45*clhs76 - clhs75;
lhs(0,8)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs81 + clhs84 + clhs91;
lhs(0,9)=DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs97 + clhs98;
lhs(0,10)=DN(0,0)*clhs99 + DN(0,1)*clhs101 + DN(0,2)*clhs103 + clhs104;
lhs(0,11)=DN(2,0)*clhs23 - DN(2,0)*clhs25 - clhs105 + clhs106*clhs45;
lhs(0,12)=DN(0,0)*clhs107 + DN(0,1)*clhs109 + DN(0,2)*clhs111 + clhs114 + clhs121;
lhs(0,13)=DN(0,0)*clhs122 + DN(0,1)*clhs124 + DN(0,2)*clhs127 + clhs128;
lhs(0,14)=DN(0,0)*clhs129 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs134;
lhs(0,15)=DN(3,0)*clhs23 - DN(3,0)*clhs25 - clhs135 + clhs136*clhs45;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs38;
lhs(1,1)=DN(0,0)*clhs33 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs10*clhs142 + clhs30;
lhs(1,2)=DN(0,0)*clhs41 + DN(0,1)*clhs143 + DN(0,2)*clhs145 + clhs147;
lhs(1,3)=DN(0,1)*clhs46;
lhs(1,4)=DN(0,0)*clhs49 + DN(0,1)*clhs148 + DN(0,2)*clhs149 + clhs150;
lhs(1,5)=DN(0,0)*clhs64 + DN(0,1)*clhs151 + DN(0,2)*clhs153 + clhs155 + clhs156;
lhs(1,6)=DN(0,0)*clhs71 + DN(0,1)*clhs157 + DN(0,2)*clhs159 + clhs160;
lhs(1,7)=DN(1,1)*clhs23 - DN(1,1)*clhs25 - clhs161 + clhs162*clhs45;
lhs(1,8)=DN(0,0)*clhs79 + DN(0,1)*clhs163 + DN(0,2)*clhs164 + clhs165;
lhs(1,9)=DN(0,0)*clhs94 + DN(0,1)*clhs166 + DN(0,2)*clhs168 + clhs170 + clhs171;
lhs(1,10)=DN(0,0)*clhs101 + DN(0,1)*clhs172 + DN(0,2)*clhs174 + clhs175;
lhs(1,11)=DN(2,1)*clhs23 - DN(2,1)*clhs25 - clhs176 + clhs177*clhs45;
lhs(1,12)=DN(0,0)*clhs109 + DN(0,1)*clhs178 + DN(0,2)*clhs179 + clhs180;
lhs(1,13)=DN(0,0)*clhs124 + DN(0,1)*clhs181 + DN(0,2)*clhs183 + clhs185 + clhs186;
lhs(1,14)=DN(0,0)*clhs131 + DN(0,1)*clhs187 + DN(0,2)*clhs189 + clhs190;
lhs(1,15)=DN(3,1)*clhs23 - DN(3,1)*clhs25 - clhs191 + clhs192*clhs45;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs138 + DN(0,2)*clhs193 + clhs44;
lhs(2,1)=DN(0,0)*clhs36 + DN(0,1)*clhs141 + DN(0,2)*clhs194 + clhs147;
lhs(2,2)=DN(0,0)*clhs43 + DN(0,1)*clhs145 + DN(0,2)*clhs195 + clhs10*clhs196 + clhs30;
lhs(2,3)=DN(0,2)*clhs46;
lhs(2,4)=DN(0,0)*clhs51 + DN(0,1)*clhs149 + DN(0,2)*clhs197 + clhs199;
lhs(2,5)=DN(0,0)*clhs67 + DN(0,1)*clhs153 + DN(0,2)*clhs200 + clhs201;
lhs(2,6)=DN(0,0)*clhs73 + DN(0,1)*clhs159 + DN(0,2)*clhs202 + clhs156 + clhs204;
lhs(2,7)=DN(1,2)*clhs23 - DN(1,2)*clhs25 - clhs205 + clhs206*clhs45;
lhs(2,8)=DN(0,0)*clhs81 + DN(0,1)*clhs164 + DN(0,2)*clhs207 + clhs208;
lhs(2,9)=DN(0,0)*clhs97 + DN(0,1)*clhs168 + DN(0,2)*clhs209 + clhs210;
lhs(2,10)=DN(0,0)*clhs103 + DN(0,1)*clhs174 + DN(0,2)*clhs211 + clhs171 + clhs213;
lhs(2,11)=DN(2,2)*clhs23 - DN(2,2)*clhs25 - clhs214 + clhs215*clhs45;
lhs(2,12)=DN(0,0)*clhs111 + DN(0,1)*clhs179 + DN(0,2)*clhs216 + clhs217;
lhs(2,13)=DN(0,0)*clhs127 + DN(0,1)*clhs183 + DN(0,2)*clhs218 + clhs219;
lhs(2,14)=DN(0,0)*clhs133 + DN(0,1)*clhs189 + DN(0,2)*clhs220 + clhs186 + clhs222;
lhs(2,15)=DN(3,2)*clhs23 - DN(3,2)*clhs25 - clhs223 + clhs224*clhs45;
lhs(3,0)=DN(0,0)*clhs225;
lhs(3,1)=DN(0,1)*clhs225;
lhs(3,2)=DN(0,2)*clhs225;
lhs(3,3)=clhs26*(clhs142 + clhs196 + clhs5);
lhs(3,4)=DN(0,0)*clhs60 + clhs76;
lhs(3,5)=DN(0,1)*clhs60 + clhs162;
lhs(3,6)=DN(0,2)*clhs60 + clhs206;
lhs(3,7)=clhs226;
lhs(3,8)=DN(0,0)*clhs90 + clhs106;
lhs(3,9)=DN(0,1)*clhs90 + clhs177;
lhs(3,10)=DN(0,2)*clhs90 + clhs215;
lhs(3,11)=clhs227;
lhs(3,12)=DN(0,0)*clhs120 + clhs136;
lhs(3,13)=DN(0,1)*clhs120 + clhs192;
lhs(3,14)=DN(0,2)*clhs120 + clhs224;
lhs(3,15)=clhs228;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs232 + clhs54;
lhs(4,1)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs36 + clhs150;
lhs(4,2)=DN(1,0)*clhs39 + DN(1,1)*clhs41 + DN(1,2)*clhs43 + clhs199;
lhs(4,3)=DN(0,0)*clhs229 - DN(0,0)*clhs230 + clhs45*clhs75 - clhs76;
lhs(4,4)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + DN(1,2)*clhs51 + clhs10*clhs233 + clhs235;
lhs(4,5)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs67 + clhs237;
lhs(4,6)=DN(1,0)*clhs69 + DN(1,1)*clhs71 + DN(1,2)*clhs73 + clhs238;
lhs(4,7)=DN(1,0)*clhs239;
lhs(4,8)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs81 + clhs242 + clhs243;
lhs(4,9)=DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs97 + clhs244;
lhs(4,10)=DN(1,0)*clhs99 + DN(1,1)*clhs101 + DN(1,2)*clhs103 + clhs245;
lhs(4,11)=DN(2,0)*clhs229 - DN(2,0)*clhs230 - clhs246 + clhs247*clhs45;
lhs(4,12)=DN(1,0)*clhs107 + DN(1,1)*clhs109 + DN(1,2)*clhs111 + clhs250 + clhs251;
lhs(4,13)=DN(1,0)*clhs122 + DN(1,1)*clhs124 + DN(1,2)*clhs127 + clhs252;
lhs(4,14)=DN(1,0)*clhs129 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs253;
lhs(4,15)=DN(3,0)*clhs229 - DN(3,0)*clhs230 - clhs254 + clhs255*clhs45;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs68;
lhs(5,1)=DN(1,0)*clhs33 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs155 + clhs256;
lhs(5,2)=DN(1,0)*clhs41 + DN(1,1)*clhs143 + DN(1,2)*clhs145 + clhs201;
lhs(5,3)=DN(0,1)*clhs229 - DN(0,1)*clhs230 + clhs161*clhs45 - clhs162;
lhs(5,4)=DN(1,0)*clhs49 + DN(1,1)*clhs148 + DN(1,2)*clhs149 + clhs237;
lhs(5,5)=DN(1,0)*clhs64 + DN(1,1)*clhs151 + DN(1,2)*clhs153 + clhs10*clhs257 + clhs235;
lhs(5,6)=DN(1,0)*clhs71 + DN(1,1)*clhs157 + DN(1,2)*clhs159 + clhs259;
lhs(5,7)=DN(1,1)*clhs239;
lhs(5,8)=DN(1,0)*clhs79 + DN(1,1)*clhs163 + DN(1,2)*clhs164 + clhs260;
lhs(5,9)=DN(1,0)*clhs94 + DN(1,1)*clhs166 + DN(1,2)*clhs168 + clhs262 + clhs263;
lhs(5,10)=DN(1,0)*clhs101 + DN(1,1)*clhs172 + DN(1,2)*clhs174 + clhs264;
lhs(5,11)=DN(2,1)*clhs229 - DN(2,1)*clhs230 - clhs265 + clhs266*clhs45;
lhs(5,12)=DN(1,0)*clhs109 + DN(1,1)*clhs178 + DN(1,2)*clhs179 + clhs267;
lhs(5,13)=DN(1,0)*clhs124 + DN(1,1)*clhs181 + DN(1,2)*clhs183 + clhs269 + clhs270;
lhs(5,14)=DN(1,0)*clhs131 + DN(1,1)*clhs187 + DN(1,2)*clhs189 + clhs271;
lhs(5,15)=DN(3,1)*clhs229 - DN(3,1)*clhs230 - clhs272 + clhs273*clhs45;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs138 + DN(1,2)*clhs193 + clhs74;
lhs(6,1)=DN(1,0)*clhs36 + DN(1,1)*clhs141 + DN(1,2)*clhs194 + clhs160;
lhs(6,2)=DN(1,0)*clhs43 + DN(1,1)*clhs145 + DN(1,2)*clhs195 + clhs204 + clhs256;
lhs(6,3)=DN(0,2)*clhs229 - DN(0,2)*clhs230 + clhs205*clhs45 - clhs206;
lhs(6,4)=DN(1,0)*clhs51 + DN(1,1)*clhs149 + DN(1,2)*clhs197 + clhs238;
lhs(6,5)=DN(1,0)*clhs67 + DN(1,1)*clhs153 + DN(1,2)*clhs200 + clhs259;
lhs(6,6)=DN(1,0)*clhs73 + DN(1,1)*clhs159 + DN(1,2)*clhs202 + clhs10*clhs274 + clhs235;
lhs(6,7)=DN(1,2)*clhs239;
lhs(6,8)=DN(1,0)*clhs81 + DN(1,1)*clhs164 + DN(1,2)*clhs207 + clhs276;
lhs(6,9)=DN(1,0)*clhs97 + DN(1,1)*clhs168 + DN(1,2)*clhs209 + clhs277;
lhs(6,10)=DN(1,0)*clhs103 + DN(1,1)*clhs174 + DN(1,2)*clhs211 + clhs263 + clhs279;
lhs(6,11)=DN(2,2)*clhs229 - DN(2,2)*clhs230 - clhs280 + clhs281*clhs45;
lhs(6,12)=DN(1,0)*clhs111 + DN(1,1)*clhs179 + DN(1,2)*clhs216 + clhs282;
lhs(6,13)=DN(1,0)*clhs127 + DN(1,1)*clhs183 + DN(1,2)*clhs218 + clhs283;
lhs(6,14)=DN(1,0)*clhs133 + DN(1,1)*clhs189 + DN(1,2)*clhs220 + clhs270 + clhs285;
lhs(6,15)=DN(3,2)*clhs229 - DN(3,2)*clhs230 - clhs286 + clhs287*clhs45;
lhs(7,0)=DN(1,0)*clhs27 + clhs75;
lhs(7,1)=DN(1,1)*clhs27 + clhs161;
lhs(7,2)=DN(1,2)*clhs27 + clhs205;
lhs(7,3)=clhs226;
lhs(7,4)=DN(1,0)*clhs288;
lhs(7,5)=DN(1,1)*clhs288;
lhs(7,6)=DN(1,2)*clhs288;
lhs(7,7)=clhs26*(clhs233 + clhs257 + clhs274);
lhs(7,8)=DN(1,0)*clhs90 + clhs247;
lhs(7,9)=DN(1,1)*clhs90 + clhs266;
lhs(7,10)=DN(1,2)*clhs90 + clhs281;
lhs(7,11)=clhs289;
lhs(7,12)=DN(1,0)*clhs120 + clhs255;
lhs(7,13)=DN(1,1)*clhs120 + clhs273;
lhs(7,14)=DN(1,2)*clhs120 + clhs287;
lhs(7,15)=clhs290;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs292 + clhs84;
lhs(8,1)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs36 + clhs165;
lhs(8,2)=DN(2,0)*clhs39 + DN(2,1)*clhs41 + DN(2,2)*clhs43 + clhs208;
lhs(8,3)=-DN(0,0)*clhs293 + DN(0,0)*clhs294 + clhs105*clhs45 - clhs106;
lhs(8,4)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + DN(2,2)*clhs51 + clhs242 + clhs295;
lhs(8,5)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs67 + clhs260;
lhs(8,6)=DN(2,0)*clhs69 + DN(2,1)*clhs71 + DN(2,2)*clhs73 + clhs276;
lhs(8,7)=-DN(1,0)*clhs293 + DN(1,0)*clhs294 + clhs246*clhs45 - clhs247;
lhs(8,8)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs81 + clhs10*clhs296 + clhs298;
lhs(8,9)=DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs97 + clhs300;
lhs(8,10)=DN(2,0)*clhs99 + DN(2,1)*clhs101 + DN(2,2)*clhs103 + clhs301;
lhs(8,11)=DN(2,0)*clhs302;
lhs(8,12)=DN(2,0)*clhs107 + DN(2,1)*clhs109 + DN(2,2)*clhs111 + clhs305 + clhs306;
lhs(8,13)=DN(2,0)*clhs122 + DN(2,1)*clhs124 + DN(2,2)*clhs127 + clhs307;
lhs(8,14)=DN(2,0)*clhs129 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs308;
lhs(8,15)=-DN(3,0)*clhs293 + DN(3,0)*clhs294 - clhs309 + clhs310*clhs45;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs98;
lhs(9,1)=DN(2,0)*clhs33 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs170 + clhs311;
lhs(9,2)=DN(2,0)*clhs41 + DN(2,1)*clhs143 + DN(2,2)*clhs145 + clhs210;
lhs(9,3)=-DN(0,1)*clhs293 + DN(0,1)*clhs294 + clhs176*clhs45 - clhs177;
lhs(9,4)=DN(2,0)*clhs49 + DN(2,1)*clhs148 + DN(2,2)*clhs149 + clhs244;
lhs(9,5)=DN(2,0)*clhs64 + DN(2,1)*clhs151 + DN(2,2)*clhs153 + clhs262 + clhs312;
lhs(9,6)=DN(2,0)*clhs71 + DN(2,1)*clhs157 + DN(2,2)*clhs159 + clhs277;
lhs(9,7)=-DN(1,1)*clhs293 + DN(1,1)*clhs294 + clhs265*clhs45 - clhs266;
lhs(9,8)=DN(2,0)*clhs79 + DN(2,1)*clhs163 + DN(2,2)*clhs164 + clhs300;
lhs(9,9)=DN(2,0)*clhs94 + DN(2,1)*clhs166 + DN(2,2)*clhs168 + clhs10*clhs313 + clhs298;
lhs(9,10)=DN(2,0)*clhs101 + DN(2,1)*clhs172 + DN(2,2)*clhs174 + clhs315;
lhs(9,11)=DN(2,1)*clhs302;
lhs(9,12)=DN(2,0)*clhs109 + DN(2,1)*clhs178 + DN(2,2)*clhs179 + clhs316;
lhs(9,13)=DN(2,0)*clhs124 + DN(2,1)*clhs181 + DN(2,2)*clhs183 + clhs318 + clhs319;
lhs(9,14)=DN(2,0)*clhs131 + DN(2,1)*clhs187 + DN(2,2)*clhs189 + clhs320;
lhs(9,15)=-DN(3,1)*clhs293 + DN(3,1)*clhs294 - clhs321 + clhs322*clhs45;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs138 + DN(2,2)*clhs193 + clhs104;
lhs(10,1)=DN(2,0)*clhs36 + DN(2,1)*clhs141 + DN(2,2)*clhs194 + clhs175;
lhs(10,2)=DN(2,0)*clhs43 + DN(2,1)*clhs145 + DN(2,2)*clhs195 + clhs213 + clhs311;
lhs(10,3)=-DN(0,2)*clhs293 + DN(0,2)*clhs294 + clhs214*clhs45 - clhs215;
lhs(10,4)=DN(2,0)*clhs51 + DN(2,1)*clhs149 + DN(2,2)*clhs197 + clhs245;
lhs(10,5)=DN(2,0)*clhs67 + DN(2,1)*clhs153 + DN(2,2)*clhs200 + clhs264;
lhs(10,6)=DN(2,0)*clhs73 + DN(2,1)*clhs159 + DN(2,2)*clhs202 + clhs279 + clhs312;
lhs(10,7)=-DN(1,2)*clhs293 + DN(1,2)*clhs294 + clhs280*clhs45 - clhs281;
lhs(10,8)=DN(2,0)*clhs81 + DN(2,1)*clhs164 + DN(2,2)*clhs207 + clhs301;
lhs(10,9)=DN(2,0)*clhs97 + DN(2,1)*clhs168 + DN(2,2)*clhs209 + clhs315;
lhs(10,10)=DN(2,0)*clhs103 + DN(2,1)*clhs174 + DN(2,2)*clhs211 + clhs10*clhs323 + clhs298;
lhs(10,11)=DN(2,2)*clhs302;
lhs(10,12)=DN(2,0)*clhs111 + DN(2,1)*clhs179 + DN(2,2)*clhs216 + clhs325;
lhs(10,13)=DN(2,0)*clhs127 + DN(2,1)*clhs183 + DN(2,2)*clhs218 + clhs326;
lhs(10,14)=DN(2,0)*clhs133 + DN(2,1)*clhs189 + DN(2,2)*clhs220 + clhs319 + clhs328;
lhs(10,15)=-DN(3,2)*clhs293 + DN(3,2)*clhs294 - clhs329 + clhs330*clhs45;
lhs(11,0)=DN(2,0)*clhs27 + clhs105;
lhs(11,1)=DN(2,1)*clhs27 + clhs176;
lhs(11,2)=DN(2,2)*clhs27 + clhs214;
lhs(11,3)=clhs227;
lhs(11,4)=DN(2,0)*clhs60 + clhs246;
lhs(11,5)=DN(2,1)*clhs60 + clhs265;
lhs(11,6)=DN(2,2)*clhs60 + clhs280;
lhs(11,7)=clhs289;
lhs(11,8)=DN(2,0)*clhs331;
lhs(11,9)=DN(2,1)*clhs331;
lhs(11,10)=DN(2,2)*clhs331;
lhs(11,11)=clhs26*(clhs296 + clhs313 + clhs323);
lhs(11,12)=DN(2,0)*clhs120 + clhs310;
lhs(11,13)=DN(2,1)*clhs120 + clhs322;
lhs(11,14)=DN(2,2)*clhs120 + clhs330;
lhs(11,15)=clhs332;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs114 + clhs334;
lhs(12,1)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs36 + clhs180;
lhs(12,2)=DN(3,0)*clhs39 + DN(3,1)*clhs41 + DN(3,2)*clhs43 + clhs217;
lhs(12,3)=-DN(0,0)*clhs335 + DN(0,0)*clhs336 + clhs135*clhs45 - clhs136;
lhs(12,4)=DN(3,0)*clhs47 + DN(3,1)*clhs49 + DN(3,2)*clhs51 + clhs250 + clhs337;
lhs(12,5)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs67 + clhs267;
lhs(12,6)=DN(3,0)*clhs69 + DN(3,1)*clhs71 + DN(3,2)*clhs73 + clhs282;
lhs(12,7)=-DN(1,0)*clhs335 + DN(1,0)*clhs336 + clhs254*clhs45 - clhs255;
lhs(12,8)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs81 + clhs305 + clhs338;
lhs(12,9)=DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs97 + clhs316;
lhs(12,10)=DN(3,0)*clhs99 + DN(3,1)*clhs101 + DN(3,2)*clhs103 + clhs325;
lhs(12,11)=-DN(2,0)*clhs335 + DN(2,0)*clhs336 + clhs309*clhs45 - clhs310;
lhs(12,12)=DN(3,0)*clhs107 + DN(3,1)*clhs109 + DN(3,2)*clhs111 + clhs10*clhs339 + clhs341;
lhs(12,13)=DN(3,0)*clhs122 + DN(3,1)*clhs124 + DN(3,2)*clhs127 + clhs343;
lhs(12,14)=DN(3,0)*clhs129 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs344;
lhs(12,15)=DN(3,0)*clhs345;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs128;
lhs(13,1)=DN(3,0)*clhs33 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs185 + clhs346;
lhs(13,2)=DN(3,0)*clhs41 + DN(3,1)*clhs143 + DN(3,2)*clhs145 + clhs219;
lhs(13,3)=-DN(0,1)*clhs335 + DN(0,1)*clhs336 + clhs191*clhs45 - clhs192;
lhs(13,4)=DN(3,0)*clhs49 + DN(3,1)*clhs148 + DN(3,2)*clhs149 + clhs252;
lhs(13,5)=DN(3,0)*clhs64 + DN(3,1)*clhs151 + DN(3,2)*clhs153 + clhs269 + clhs347;
lhs(13,6)=DN(3,0)*clhs71 + DN(3,1)*clhs157 + DN(3,2)*clhs159 + clhs283;
lhs(13,7)=-DN(1,1)*clhs335 + DN(1,1)*clhs336 + clhs272*clhs45 - clhs273;
lhs(13,8)=DN(3,0)*clhs79 + DN(3,1)*clhs163 + DN(3,2)*clhs164 + clhs307;
lhs(13,9)=DN(3,0)*clhs94 + DN(3,1)*clhs166 + DN(3,2)*clhs168 + clhs318 + clhs348;
lhs(13,10)=DN(3,0)*clhs101 + DN(3,1)*clhs172 + DN(3,2)*clhs174 + clhs326;
lhs(13,11)=-DN(2,1)*clhs335 + DN(2,1)*clhs336 + clhs321*clhs45 - clhs322;
lhs(13,12)=DN(3,0)*clhs109 + DN(3,1)*clhs178 + DN(3,2)*clhs179 + clhs343;
lhs(13,13)=DN(3,0)*clhs124 + DN(3,1)*clhs181 + DN(3,2)*clhs183 + clhs10*clhs349 + clhs341;
lhs(13,14)=DN(3,0)*clhs131 + DN(3,1)*clhs187 + DN(3,2)*clhs189 + clhs350;
lhs(13,15)=DN(3,1)*clhs345;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs138 + DN(3,2)*clhs193 + clhs134;
lhs(14,1)=DN(3,0)*clhs36 + DN(3,1)*clhs141 + DN(3,2)*clhs194 + clhs190;
lhs(14,2)=DN(3,0)*clhs43 + DN(3,1)*clhs145 + DN(3,2)*clhs195 + clhs222 + clhs346;
lhs(14,3)=-DN(0,2)*clhs335 + DN(0,2)*clhs336 + clhs223*clhs45 - clhs224;
lhs(14,4)=DN(3,0)*clhs51 + DN(3,1)*clhs149 + DN(3,2)*clhs197 + clhs253;
lhs(14,5)=DN(3,0)*clhs67 + DN(3,1)*clhs153 + DN(3,2)*clhs200 + clhs271;
lhs(14,6)=DN(3,0)*clhs73 + DN(3,1)*clhs159 + DN(3,2)*clhs202 + clhs285 + clhs347;
lhs(14,7)=-DN(1,2)*clhs335 + DN(1,2)*clhs336 + clhs286*clhs45 - clhs287;
lhs(14,8)=DN(3,0)*clhs81 + DN(3,1)*clhs164 + DN(3,2)*clhs207 + clhs308;
lhs(14,9)=DN(3,0)*clhs97 + DN(3,1)*clhs168 + DN(3,2)*clhs209 + clhs320;
lhs(14,10)=DN(3,0)*clhs103 + DN(3,1)*clhs174 + DN(3,2)*clhs211 + clhs328 + clhs348;
lhs(14,11)=-DN(2,2)*clhs335 + DN(2,2)*clhs336 + clhs329*clhs45 - clhs330;
lhs(14,12)=DN(3,0)*clhs111 + DN(3,1)*clhs179 + DN(3,2)*clhs216 + clhs344;
lhs(14,13)=DN(3,0)*clhs127 + DN(3,1)*clhs183 + DN(3,2)*clhs218 + clhs350;
lhs(14,14)=DN(3,0)*clhs133 + DN(3,1)*clhs189 + DN(3,2)*clhs220 + clhs10*clhs351 + clhs341;
lhs(14,15)=DN(3,2)*clhs345;
lhs(15,0)=DN(3,0)*clhs27 + clhs135;
lhs(15,1)=DN(3,1)*clhs27 + clhs191;
lhs(15,2)=DN(3,2)*clhs27 + clhs223;
lhs(15,3)=clhs228;
lhs(15,4)=DN(3,0)*clhs60 + clhs254;
lhs(15,5)=DN(3,1)*clhs60 + clhs272;
lhs(15,6)=DN(3,2)*clhs60 + clhs286;
lhs(15,7)=clhs290;
lhs(15,8)=DN(3,0)*clhs90 + clhs309;
lhs(15,9)=DN(3,1)*clhs90 + clhs321;
lhs(15,10)=DN(3,2)*clhs90 + clhs329;
lhs(15,11)=clhs332;
lhs(15,12)=DN(3,0)*clhs352;
lhs(15,13)=DN(3,1)*clhs352;
lhs(15,14)=DN(3,2)*clhs352;
lhs(15,15)=clhs26*(clhs339 + clhs349 + clhs351);


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
    // const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vnnn = rData.Velocity_OldStep3;

    const auto &vmesh = rData.MeshVelocity;
    // const auto &vconv = v - vn;
    const auto &vconv = v - vmesh;

    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;
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
const double crhs3 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0));
const double crhs4 = N[0]*rho;
const double crhs5 = bdf0*crhs4;
const double crhs6 = N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0)) + N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0)) + N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs8 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crhs9 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double crhs10 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crhs11 = rho*(crhs10*crhs9 + crhs7*crhs8);
const double crhs12 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs13 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs14 = rho*(crhs12*crhs7 + crhs13*crhs9);
const double crhs15 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crhs16 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crhs17 = crhs15*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crhs16*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0));
const double crhs18 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs19 = crhs18 + crhs7 - volume_error_ratio;
const double crhs20 = rho*stab_c2*sqrt(pow(crhs12, 2) + pow(crhs13, 2));
const double crhs21 = crhs19*(crhs20*h/stab_c1 + mu);
const double crhs22 = bdf0*rho;
const double crhs23 = crhs22*crhs3;
const double crhs24 = 1.0/(K_darcy + crhs20/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs25 = crhs24*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 - crhs11 + crhs14 + crhs2 + crhs23 + rho*(crhs17 + crhs6));
const double crhs26 = K_darcy*N[0];
const double crhs27 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs28 = crhs27*crhs4;
const double crhs29 = rho*(DN(0,0)*crhs12 + DN(0,1)*crhs13);
const double crhs30 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs31 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs32 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1));
const double crhs33 = N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1)) + N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1)) + N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1));
const double crhs34 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs35 = rho*(crhs10*crhs18 + crhs34*crhs8);
const double crhs36 = rho*(crhs12*crhs34 + crhs13*crhs18);
const double crhs37 = crhs15*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crhs16*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1));
const double crhs38 = crhs22*crhs32;
const double crhs39 = crhs24*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs30 + crhs31 - crhs35 + crhs36 + crhs38 + rho*(crhs33 + crhs37));
const double crhs40 = N[1]*rho;
const double crhs41 = K_darcy*N[1];
const double crhs42 = crhs27*crhs40;
const double crhs43 = rho*(DN(1,0)*crhs12 + DN(1,1)*crhs13);
const double crhs44 = N[2]*rho;
const double crhs45 = K_darcy*N[2];
const double crhs46 = crhs27*crhs44;
const double crhs47 = rho*(DN(2,0)*crhs12 + DN(2,1)*crhs13);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs21 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 + N[0]*crhs11 - N[0]*crhs14 - N[0]*crhs2 - crhs17*crhs4 + crhs25*crhs26 - crhs25*crhs28 - crhs25*crhs29 - crhs3*crhs5 - crhs4*crhs6;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs21 - DN(0,1)*stress[1] + N[0]*crhs30 - N[0]*crhs31 + N[0]*crhs35 - N[0]*crhs36 + crhs26*crhs39 - crhs28*crhs39 - crhs29*crhs39 - crhs32*crhs5 - crhs33*crhs4 - crhs37*crhs4;
rhs[2]=-DN(0,0)*crhs25 - DN(0,1)*crhs39 - N[0]*crhs19;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs21 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 + N[1]*crhs11 - N[1]*crhs14 - N[1]*crhs2 - N[1]*crhs23 - crhs17*crhs40 + crhs25*crhs41 - crhs25*crhs42 - crhs25*crhs43 - crhs40*crhs6;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs21 - DN(1,1)*stress[1] + N[1]*crhs30 - N[1]*crhs31 + N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs38 - crhs33*crhs40 - crhs37*crhs40 + crhs39*crhs41 - crhs39*crhs42 - crhs39*crhs43;
rhs[5]=-DN(1,0)*crhs25 - DN(1,1)*crhs39 - N[1]*crhs19;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs21 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 + N[2]*crhs11 - N[2]*crhs14 - N[2]*crhs2 - N[2]*crhs23 - crhs17*crhs44 + crhs25*crhs45 - crhs25*crhs46 - crhs25*crhs47 - crhs44*crhs6;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs21 - DN(2,1)*stress[1] + N[2]*crhs30 - N[2]*crhs31 + N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs38 - crhs33*crhs44 - crhs37*crhs44 + crhs39*crhs45 - crhs39*crhs46 - crhs39*crhs47;
rhs[8]=-DN(2,0)*crhs25 - DN(2,1)*crhs39 - N[2]*crhs19;


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
    // const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vnnn = rData.Velocity_OldStep3;

    const auto &vmesh = rData.MeshVelocity;
    // const auto &vconv = v - vn;
    const auto &vconv = v - vmesh;

    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;
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
const double crhs3 = N[0]*(v(0,0) - vfrac(0,0)) + N[1]*(v(1,0) - vfrac(1,0)) + N[2]*(v(2,0) - vfrac(2,0)) + N[3]*(v(3,0) - vfrac(3,0));
const double crhs4 = N[0]*rho;
const double crhs5 = bdf0*crhs4;
const double crhs6 = N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0)) + N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0)) + N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0)) + N[3]*(bdf0*vn(3,0) + bdf1*vnn(3,0) + bdf2*vnnn(3,0));
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs8 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crhs9 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double crhs10 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crhs11 = DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0);
const double crhs12 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crhs13 = rho*(crhs10*crhs9 + crhs11*crhs12 + crhs7*crhs8);
const double crhs14 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs15 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs16 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs17 = rho*(crhs11*crhs16 + crhs14*crhs7 + crhs15*crhs9);
const double crhs18 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crhs19 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crhs20 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crhs21 = crhs18*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crhs19*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crhs20*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0));
const double crhs22 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs23 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs24 = crhs22 + crhs23 + crhs7 - volume_error_ratio;
const double crhs25 = rho*stab_c2*sqrt(pow(crhs14, 2) + pow(crhs15, 2) + pow(crhs16, 2));
const double crhs26 = crhs24*(crhs25*h/stab_c1 + mu);
const double crhs27 = bdf0*rho;
const double crhs28 = crhs27*crhs3;
const double crhs29 = 1.0/(K_darcy + crhs25/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs30 = crhs29*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 - crhs13 + crhs17 + crhs2 + crhs28 + rho*(crhs21 + crhs6));
const double crhs31 = K_darcy*N[0];
const double crhs32 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs33 = crhs32*crhs4;
const double crhs34 = rho*(DN(0,0)*crhs14 + DN(0,1)*crhs15 + DN(0,2)*crhs16);
const double crhs35 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs36 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs37 = N[0]*(v(0,1) - vfrac(0,1)) + N[1]*(v(1,1) - vfrac(1,1)) + N[2]*(v(2,1) - vfrac(2,1)) + N[3]*(v(3,1) - vfrac(3,1));
const double crhs38 = N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1)) + N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1)) + N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1)) + N[3]*(bdf0*vn(3,1) + bdf1*vnn(3,1) + bdf2*vnnn(3,1));
const double crhs39 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs40 = DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
const double crhs41 = rho*(crhs10*crhs22 + crhs12*crhs40 + crhs39*crhs8);
const double crhs42 = rho*(crhs14*crhs39 + crhs15*crhs22 + crhs16*crhs40);
const double crhs43 = crhs18*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crhs19*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crhs20*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1));
const double crhs44 = crhs27*crhs37;
const double crhs45 = crhs29*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs35 + crhs36 - crhs41 + crhs42 + crhs44 + rho*(crhs38 + crhs43));
const double crhs46 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs47 = K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs48 = N[0]*(v(0,2) - vfrac(0,2)) + N[1]*(v(1,2) - vfrac(1,2)) + N[2]*(v(2,2) - vfrac(2,2)) + N[3]*(v(3,2) - vfrac(3,2));
const double crhs49 = N[0]*(bdf0*vn(0,2) + bdf1*vnn(0,2) + bdf2*vnnn(0,2)) + N[1]*(bdf0*vn(1,2) + bdf1*vnn(1,2) + bdf2*vnnn(1,2)) + N[2]*(bdf0*vn(2,2) + bdf1*vnn(2,2) + bdf2*vnnn(2,2)) + N[3]*(bdf0*vn(3,2) + bdf1*vnn(3,2) + bdf2*vnnn(3,2));
const double crhs50 = DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
const double crhs51 = DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
const double crhs52 = rho*(crhs10*crhs51 + crhs12*crhs23 + crhs50*crhs8);
const double crhs53 = rho*(crhs14*crhs50 + crhs15*crhs51 + crhs16*crhs23);
const double crhs54 = crhs18*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crhs19*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crhs20*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2));
const double crhs55 = crhs27*crhs48;
const double crhs56 = crhs29*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs46 + crhs47 - crhs52 + crhs53 + crhs55 + rho*(crhs49 + crhs54));
const double crhs57 = N[1]*rho;
const double crhs58 = K_darcy*N[1];
const double crhs59 = crhs32*crhs57;
const double crhs60 = rho*(DN(1,0)*crhs14 + DN(1,1)*crhs15 + DN(1,2)*crhs16);
const double crhs61 = N[2]*rho;
const double crhs62 = K_darcy*N[2];
const double crhs63 = crhs32*crhs61;
const double crhs64 = rho*(DN(2,0)*crhs14 + DN(2,1)*crhs15 + DN(2,2)*crhs16);
const double crhs65 = N[3]*rho;
const double crhs66 = K_darcy*N[3];
const double crhs67 = crhs32*crhs65;
const double crhs68 = rho*(DN(3,0)*crhs14 + DN(3,1)*crhs15 + DN(3,2)*crhs16);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs26 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 + N[0]*crhs13 - N[0]*crhs17 - N[0]*crhs2 - crhs21*crhs4 - crhs3*crhs5 + crhs30*crhs31 - crhs30*crhs33 - crhs30*crhs34 - crhs4*crhs6;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs26 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs35 - N[0]*crhs36 + N[0]*crhs41 - N[0]*crhs42 + crhs31*crhs45 - crhs33*crhs45 - crhs34*crhs45 - crhs37*crhs5 - crhs38*crhs4 - crhs4*crhs43;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs26 - DN(0,2)*stress[2] + N[0]*crhs46 - N[0]*crhs47 + N[0]*crhs52 - N[0]*crhs53 + crhs31*crhs56 - crhs33*crhs56 - crhs34*crhs56 - crhs4*crhs49 - crhs4*crhs54 - crhs48*crhs5;
rhs[3]=-DN(0,0)*crhs30 - DN(0,1)*crhs45 - DN(0,2)*crhs56 - N[0]*crhs24;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs26 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 + N[1]*crhs13 - N[1]*crhs17 - N[1]*crhs2 - N[1]*crhs28 - crhs21*crhs57 + crhs30*crhs58 - crhs30*crhs59 - crhs30*crhs60 - crhs57*crhs6;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs26 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs35 - N[1]*crhs36 + N[1]*crhs41 - N[1]*crhs42 - N[1]*crhs44 - crhs38*crhs57 - crhs43*crhs57 + crhs45*crhs58 - crhs45*crhs59 - crhs45*crhs60;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs26 - DN(1,2)*stress[2] + N[1]*crhs46 - N[1]*crhs47 + N[1]*crhs52 - N[1]*crhs53 - N[1]*crhs55 - crhs49*crhs57 - crhs54*crhs57 + crhs56*crhs58 - crhs56*crhs59 - crhs56*crhs60;
rhs[7]=-DN(1,0)*crhs30 - DN(1,1)*crhs45 - DN(1,2)*crhs56 - N[1]*crhs24;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs26 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 + N[2]*crhs13 - N[2]*crhs17 - N[2]*crhs2 - N[2]*crhs28 - crhs21*crhs61 + crhs30*crhs62 - crhs30*crhs63 - crhs30*crhs64 - crhs6*crhs61;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs26 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs35 - N[2]*crhs36 + N[2]*crhs41 - N[2]*crhs42 - N[2]*crhs44 - crhs38*crhs61 - crhs43*crhs61 + crhs45*crhs62 - crhs45*crhs63 - crhs45*crhs64;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs26 - DN(2,2)*stress[2] + N[2]*crhs46 - N[2]*crhs47 + N[2]*crhs52 - N[2]*crhs53 - N[2]*crhs55 - crhs49*crhs61 - crhs54*crhs61 + crhs56*crhs62 - crhs56*crhs63 - crhs56*crhs64;
rhs[11]=-DN(2,0)*crhs30 - DN(2,1)*crhs45 - DN(2,2)*crhs56 - N[2]*crhs24;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs26 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 + N[3]*crhs13 - N[3]*crhs17 - N[3]*crhs2 - N[3]*crhs28 - crhs21*crhs65 + crhs30*crhs66 - crhs30*crhs67 - crhs30*crhs68 - crhs6*crhs65;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs26 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs35 - N[3]*crhs36 + N[3]*crhs41 - N[3]*crhs42 - N[3]*crhs44 - crhs38*crhs65 - crhs43*crhs65 + crhs45*crhs66 - crhs45*crhs67 - crhs45*crhs68;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs26 - DN(3,2)*stress[2] + N[3]*crhs46 - N[3]*crhs47 + N[3]*crhs52 - N[3]*crhs53 - N[3]*crhs55 - crhs49*crhs65 - crhs54*crhs65 + crhs56*crhs66 - crhs56*crhs67 - crhs56*crhs68;
rhs[15]=-DN(3,0)*crhs30 - DN(3,1)*crhs45 - DN(3,2)*crhs56 - N[3]*crhs24;


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
    // const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vnnn = rData.Velocity_OldStep3;

    const auto &vmesh = rData.MeshVelocity;
    // const auto &vconv = v - vn;
    const auto &vconv = v - vmesh;

    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;
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
    if (rData.IsAir())
    {
        volume_error_ratio=0.0;
    }
    else{
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
const double cV3 = DNenr(0,0)*cV2;
const double cV4 = K_darcy*N[0];
const double cV5 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV6 = N[0]*cV5;
const double cV7 = rho*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV8 = DNenr(1,0)*cV2;
const double cV9 = DNenr(2,0)*cV2;
const double cV10 = DNenr(0,1)*cV2;
const double cV11 = DNenr(1,1)*cV2;
const double cV12 = DNenr(2,1)*cV2;
const double cV13 = K_darcy*N[1];
const double cV14 = N[1]*cV5;
const double cV15 = rho*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV16 = K_darcy*N[2];
const double cV17 = N[2]*cV5;
const double cV18 = rho*(DN(2,0)*cV0 + DN(2,1)*cV1);
V(0,0)=-DN(0,0)*Nenr[0] - cV3*cV4 + cV3*cV6 + cV3*cV7;
V(0,1)=-DN(0,0)*Nenr[1] - cV4*cV8 + cV6*cV8 + cV7*cV8;
V(0,2)=-DN(0,0)*Nenr[2] - cV4*cV9 + cV6*cV9 + cV7*cV9;
V(1,0)=-DN(0,1)*Nenr[0] - cV10*cV4 + cV10*cV6 + cV10*cV7;
V(1,1)=-DN(0,1)*Nenr[1] - cV11*cV4 + cV11*cV6 + cV11*cV7;
V(1,2)=-DN(0,1)*Nenr[2] - cV12*cV4 + cV12*cV6 + cV12*cV7;
V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=-DN(1,0)*Nenr[0] - cV13*cV3 + cV14*cV3 + cV15*cV3;
V(3,1)=-DN(1,0)*Nenr[1] - cV13*cV8 + cV14*cV8 + cV15*cV8;
V(3,2)=-DN(1,0)*Nenr[2] - cV13*cV9 + cV14*cV9 + cV15*cV9;
V(4,0)=-DN(1,1)*Nenr[0] - cV10*cV13 + cV10*cV14 + cV10*cV15;
V(4,1)=-DN(1,1)*Nenr[1] - cV11*cV13 + cV11*cV14 + cV11*cV15;
V(4,2)=-DN(1,1)*Nenr[2] - cV12*cV13 + cV12*cV14 + cV12*cV15;
V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=-DN(2,0)*Nenr[0] - cV16*cV3 + cV17*cV3 + cV18*cV3;
V(6,1)=-DN(2,0)*Nenr[1] - cV16*cV8 + cV17*cV8 + cV18*cV8;
V(6,2)=-DN(2,0)*Nenr[2] - cV16*cV9 + cV17*cV9 + cV18*cV9;
V(7,0)=-DN(2,1)*Nenr[0] - cV10*cV16 + cV10*cV17 + cV10*cV18;
V(7,1)=-DN(2,1)*Nenr[1] - cV11*cV16 + cV11*cV17 + cV11*cV18;
V(7,2)=-DN(2,1)*Nenr[2] - cV12*cV16 + cV12*cV17 + cV12*cV18;
V(8,0)=cV2*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
V(8,1)=cV2*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
V(8,2)=cV2*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double cH1 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double cH2 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH3 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH4 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH2, 2) + pow(cH3, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH5 = cH4*(K_darcy*N[0] - rho*(DN(0,0)*cH0 + DN(0,1)*cH1) + rho*(DN(0,0)*cH2 + DN(0,1)*cH3 + N[0]*bdf0));
const double cH6 = cH4*(K_darcy*N[1] - rho*(DN(1,0)*cH0 + DN(1,1)*cH1) + rho*(DN(1,0)*cH2 + DN(1,1)*cH3 + N[1]*bdf0));
const double cH7 = cH4*(K_darcy*N[2] - rho*(DN(2,0)*cH0 + DN(2,1)*cH1) + rho*(DN(2,0)*cH2 + DN(2,1)*cH3 + N[2]*bdf0));
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH5;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH5;
H(0,2)=cH4*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
H(0,3)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH6;
H(0,4)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH6;
H(0,5)=cH4*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
H(0,6)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH7;
H(0,7)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH7;
H(0,8)=cH4*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH5;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH5;
H(1,2)=cH4*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
H(1,3)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH6;
H(1,4)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH6;
H(1,5)=cH4*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
H(1,6)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH7;
H(1,7)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH7;
H(1,8)=cH4*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH5;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH5;
H(2,2)=cH4*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
H(2,3)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH6;
H(2,4)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH6;
H(2,5)=cH4*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
H(2,6)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH7;
H(2,7)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH7;
H(2,8)=cH4*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee3 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0);
const double crhs_ee4 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double crhs_ee5 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1);
const double crhs_ee6 = N[0]*bdf0;
const double crhs_ee7 = N[1]*bdf0;
const double crhs_ee8 = N[2]*bdf0;
const double crhs_ee9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee11 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crhs_ee12 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crhs_ee13 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee14 = crhs_ee13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(crhs_ee0*crhs_ee3 + crhs_ee4*crhs_ee5) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0)) + N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0)) + N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0)) + crhs_ee11*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0)) + crhs_ee12*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0))) + rho*(crhs_ee0*crhs_ee9 + crhs_ee10*crhs_ee4 + crhs_ee6*(v(0,0) - vfrac(0,0)) + crhs_ee7*(v(1,0) - vfrac(1,0)) + crhs_ee8*(v(2,0) - vfrac(2,0))));
const double crhs_ee15 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double crhs_ee16 = crhs_ee13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(crhs_ee1*crhs_ee5 + crhs_ee15*crhs_ee3) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1)) + N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1)) + N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1)) + crhs_ee11*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1)) + crhs_ee12*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1))) + rho*(crhs_ee1*crhs_ee10 + crhs_ee15*crhs_ee9 + crhs_ee6*(v(0,1) - vfrac(0,1)) + crhs_ee7*(v(1,1) - vfrac(1,1)) + crhs_ee8*(v(2,1) - vfrac(2,1))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee14 - DNenr(0,1)*crhs_ee16 - Nenr[0]*crhs_ee2;
rhs_ee[1]=-DNenr(1,0)*crhs_ee14 - DNenr(1,1)*crhs_ee16 - Nenr[1]*crhs_ee2;
rhs_ee[2]=-DNenr(2,0)*crhs_ee14 - DNenr(2,1)*crhs_ee16 - Nenr[2]*crhs_ee2;


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
    const auto &vnnn = rData.Velocity_OldStep3;

    // const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    // const auto &vconv = v - vn;
    const auto &vconv = v - vmesh;

    const auto vfrac = rData.Velocity_Fractional;
    bool not_air_traj = rData.NotAirTraj;
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
const double cV4 = DNenr(0,0)*cV3;
const double cV5 = K_darcy*N[0];
const double cV6 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV7 = N[0]*cV6;
const double cV8 = rho*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV9 = DNenr(1,0)*cV3;
const double cV10 = DNenr(2,0)*cV3;
const double cV11 = DNenr(3,0)*cV3;
const double cV12 = DNenr(0,1)*cV3;
const double cV13 = DNenr(1,1)*cV3;
const double cV14 = DNenr(2,1)*cV3;
const double cV15 = DNenr(3,1)*cV3;
const double cV16 = DNenr(0,2)*cV3;
const double cV17 = DNenr(1,2)*cV3;
const double cV18 = DNenr(2,2)*cV3;
const double cV19 = DNenr(3,2)*cV3;
const double cV20 = K_darcy*N[1];
const double cV21 = N[1]*cV6;
const double cV22 = rho*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV23 = K_darcy*N[2];
const double cV24 = N[2]*cV6;
const double cV25 = rho*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV26 = K_darcy*N[3];
const double cV27 = N[3]*cV6;
const double cV28 = rho*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
V(0,0)=-DN(0,0)*Nenr[0] - cV4*cV5 + cV4*cV7 + cV4*cV8;
V(0,1)=-DN(0,0)*Nenr[1] - cV5*cV9 + cV7*cV9 + cV8*cV9;
V(0,2)=-DN(0,0)*Nenr[2] - cV10*cV5 + cV10*cV7 + cV10*cV8;
V(0,3)=-DN(0,0)*Nenr[3] - cV11*cV5 + cV11*cV7 + cV11*cV8;
V(1,0)=-DN(0,1)*Nenr[0] - cV12*cV5 + cV12*cV7 + cV12*cV8;
V(1,1)=-DN(0,1)*Nenr[1] - cV13*cV5 + cV13*cV7 + cV13*cV8;
V(1,2)=-DN(0,1)*Nenr[2] - cV14*cV5 + cV14*cV7 + cV14*cV8;
V(1,3)=-DN(0,1)*Nenr[3] - cV15*cV5 + cV15*cV7 + cV15*cV8;
V(2,0)=-DN(0,2)*Nenr[0] - cV16*cV5 + cV16*cV7 + cV16*cV8;
V(2,1)=-DN(0,2)*Nenr[1] - cV17*cV5 + cV17*cV7 + cV17*cV8;
V(2,2)=-DN(0,2)*Nenr[2] - cV18*cV5 + cV18*cV7 + cV18*cV8;
V(2,3)=-DN(0,2)*Nenr[3] - cV19*cV5 + cV19*cV7 + cV19*cV8;
V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
V(4,0)=-DN(1,0)*Nenr[0] - cV20*cV4 + cV21*cV4 + cV22*cV4;
V(4,1)=-DN(1,0)*Nenr[1] - cV20*cV9 + cV21*cV9 + cV22*cV9;
V(4,2)=-DN(1,0)*Nenr[2] - cV10*cV20 + cV10*cV21 + cV10*cV22;
V(4,3)=-DN(1,0)*Nenr[3] - cV11*cV20 + cV11*cV21 + cV11*cV22;
V(5,0)=-DN(1,1)*Nenr[0] - cV12*cV20 + cV12*cV21 + cV12*cV22;
V(5,1)=-DN(1,1)*Nenr[1] - cV13*cV20 + cV13*cV21 + cV13*cV22;
V(5,2)=-DN(1,1)*Nenr[2] - cV14*cV20 + cV14*cV21 + cV14*cV22;
V(5,3)=-DN(1,1)*Nenr[3] - cV15*cV20 + cV15*cV21 + cV15*cV22;
V(6,0)=-DN(1,2)*Nenr[0] - cV16*cV20 + cV16*cV21 + cV16*cV22;
V(6,1)=-DN(1,2)*Nenr[1] - cV17*cV20 + cV17*cV21 + cV17*cV22;
V(6,2)=-DN(1,2)*Nenr[2] - cV18*cV20 + cV18*cV21 + cV18*cV22;
V(6,3)=-DN(1,2)*Nenr[3] - cV19*cV20 + cV19*cV21 + cV19*cV22;
V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
V(8,0)=-DN(2,0)*Nenr[0] - cV23*cV4 + cV24*cV4 + cV25*cV4;
V(8,1)=-DN(2,0)*Nenr[1] - cV23*cV9 + cV24*cV9 + cV25*cV9;
V(8,2)=-DN(2,0)*Nenr[2] - cV10*cV23 + cV10*cV24 + cV10*cV25;
V(8,3)=-DN(2,0)*Nenr[3] - cV11*cV23 + cV11*cV24 + cV11*cV25;
V(9,0)=-DN(2,1)*Nenr[0] - cV12*cV23 + cV12*cV24 + cV12*cV25;
V(9,1)=-DN(2,1)*Nenr[1] - cV13*cV23 + cV13*cV24 + cV13*cV25;
V(9,2)=-DN(2,1)*Nenr[2] - cV14*cV23 + cV14*cV24 + cV14*cV25;
V(9,3)=-DN(2,1)*Nenr[3] - cV15*cV23 + cV15*cV24 + cV15*cV25;
V(10,0)=-DN(2,2)*Nenr[0] - cV16*cV23 + cV16*cV24 + cV16*cV25;
V(10,1)=-DN(2,2)*Nenr[1] - cV17*cV23 + cV17*cV24 + cV17*cV25;
V(10,2)=-DN(2,2)*Nenr[2] - cV18*cV23 + cV18*cV24 + cV18*cV25;
V(10,3)=-DN(2,2)*Nenr[3] - cV19*cV23 + cV19*cV24 + cV19*cV25;
V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
V(12,0)=-DN(3,0)*Nenr[0] - cV26*cV4 + cV27*cV4 + cV28*cV4;
V(12,1)=-DN(3,0)*Nenr[1] - cV26*cV9 + cV27*cV9 + cV28*cV9;
V(12,2)=-DN(3,0)*Nenr[2] - cV10*cV26 + cV10*cV27 + cV10*cV28;
V(12,3)=-DN(3,0)*Nenr[3] - cV11*cV26 + cV11*cV27 + cV11*cV28;
V(13,0)=-DN(3,1)*Nenr[0] - cV12*cV26 + cV12*cV27 + cV12*cV28;
V(13,1)=-DN(3,1)*Nenr[1] - cV13*cV26 + cV13*cV27 + cV13*cV28;
V(13,2)=-DN(3,1)*Nenr[2] - cV14*cV26 + cV14*cV27 + cV14*cV28;
V(13,3)=-DN(3,1)*Nenr[3] - cV15*cV26 + cV15*cV27 + cV15*cV28;
V(14,0)=-DN(3,2)*Nenr[0] - cV16*cV26 + cV16*cV27 + cV16*cV28;
V(14,1)=-DN(3,2)*Nenr[1] - cV17*cV26 + cV17*cV27 + cV17*cV28;
V(14,2)=-DN(3,2)*Nenr[2] - cV18*cV26 + cV18*cV27 + cV18*cV28;
V(14,3)=-DN(3,2)*Nenr[3] - cV19*cV26 + cV19*cV27 + cV19*cV28;
V(15,0)=cV3*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
V(15,1)=cV3*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
V(15,2)=cV3*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
V(15,3)=cV3*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double cH1 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double cH2 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double cH3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH5 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH6 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH3, 2) + pow(cH4, 2) + pow(cH5, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH7 = cH6*(K_darcy*N[0] - rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2) + rho*(DN(0,0)*cH3 + DN(0,1)*cH4 + DN(0,2)*cH5 + N[0]*bdf0));
const double cH8 = cH6*(K_darcy*N[1] - rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2) + rho*(DN(1,0)*cH3 + DN(1,1)*cH4 + DN(1,2)*cH5 + N[1]*bdf0));
const double cH9 = cH6*(K_darcy*N[2] - rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2) + rho*(DN(2,0)*cH3 + DN(2,1)*cH4 + DN(2,2)*cH5 + N[2]*bdf0));
const double cH10 = cH6*(K_darcy*N[3] - rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2) + rho*(DN(3,0)*cH3 + DN(3,1)*cH4 + DN(3,2)*cH5 + N[3]*bdf0));
H(0,0)=DN(0,0)*Nenr[0] + DNenr(0,0)*cH7;
H(0,1)=DN(0,1)*Nenr[0] + DNenr(0,1)*cH7;
H(0,2)=DN(0,2)*Nenr[0] + DNenr(0,2)*cH7;
H(0,3)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
H(0,4)=DN(1,0)*Nenr[0] + DNenr(0,0)*cH8;
H(0,5)=DN(1,1)*Nenr[0] + DNenr(0,1)*cH8;
H(0,6)=DN(1,2)*Nenr[0] + DNenr(0,2)*cH8;
H(0,7)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
H(0,8)=DN(2,0)*Nenr[0] + DNenr(0,0)*cH9;
H(0,9)=DN(2,1)*Nenr[0] + DNenr(0,1)*cH9;
H(0,10)=DN(2,2)*Nenr[0] + DNenr(0,2)*cH9;
H(0,11)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
H(0,12)=DN(3,0)*Nenr[0] + DNenr(0,0)*cH10;
H(0,13)=DN(3,1)*Nenr[0] + DNenr(0,1)*cH10;
H(0,14)=DN(3,2)*Nenr[0] + DNenr(0,2)*cH10;
H(0,15)=cH6*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
H(1,0)=DN(0,0)*Nenr[1] + DNenr(1,0)*cH7;
H(1,1)=DN(0,1)*Nenr[1] + DNenr(1,1)*cH7;
H(1,2)=DN(0,2)*Nenr[1] + DNenr(1,2)*cH7;
H(1,3)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
H(1,4)=DN(1,0)*Nenr[1] + DNenr(1,0)*cH8;
H(1,5)=DN(1,1)*Nenr[1] + DNenr(1,1)*cH8;
H(1,6)=DN(1,2)*Nenr[1] + DNenr(1,2)*cH8;
H(1,7)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
H(1,8)=DN(2,0)*Nenr[1] + DNenr(1,0)*cH9;
H(1,9)=DN(2,1)*Nenr[1] + DNenr(1,1)*cH9;
H(1,10)=DN(2,2)*Nenr[1] + DNenr(1,2)*cH9;
H(1,11)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
H(1,12)=DN(3,0)*Nenr[1] + DNenr(1,0)*cH10;
H(1,13)=DN(3,1)*Nenr[1] + DNenr(1,1)*cH10;
H(1,14)=DN(3,2)*Nenr[1] + DNenr(1,2)*cH10;
H(1,15)=cH6*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
H(2,0)=DN(0,0)*Nenr[2] + DNenr(2,0)*cH7;
H(2,1)=DN(0,1)*Nenr[2] + DNenr(2,1)*cH7;
H(2,2)=DN(0,2)*Nenr[2] + DNenr(2,2)*cH7;
H(2,3)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
H(2,4)=DN(1,0)*Nenr[2] + DNenr(2,0)*cH8;
H(2,5)=DN(1,1)*Nenr[2] + DNenr(2,1)*cH8;
H(2,6)=DN(1,2)*Nenr[2] + DNenr(2,2)*cH8;
H(2,7)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
H(2,8)=DN(2,0)*Nenr[2] + DNenr(2,0)*cH9;
H(2,9)=DN(2,1)*Nenr[2] + DNenr(2,1)*cH9;
H(2,10)=DN(2,2)*Nenr[2] + DNenr(2,2)*cH9;
H(2,11)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
H(2,12)=DN(3,0)*Nenr[2] + DNenr(2,0)*cH10;
H(2,13)=DN(3,1)*Nenr[2] + DNenr(2,1)*cH10;
H(2,14)=DN(3,2)*Nenr[2] + DNenr(2,2)*cH10;
H(2,15)=cH6*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
H(3,0)=DN(0,0)*Nenr[3] + DNenr(3,0)*cH7;
H(3,1)=DN(0,1)*Nenr[3] + DNenr(3,1)*cH7;
H(3,2)=DN(0,2)*Nenr[3] + DNenr(3,2)*cH7;
H(3,3)=cH6*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
H(3,4)=DN(1,0)*Nenr[3] + DNenr(3,0)*cH8;
H(3,5)=DN(1,1)*Nenr[3] + DNenr(3,1)*cH8;
H(3,6)=DN(1,2)*Nenr[3] + DNenr(3,2)*cH8;
H(3,7)=cH6*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
H(3,8)=DN(2,0)*Nenr[3] + DNenr(3,0)*cH9;
H(3,9)=DN(2,1)*Nenr[3] + DNenr(3,1)*cH9;
H(3,10)=DN(2,2)*Nenr[3] + DNenr(3,2)*cH9;
H(3,11)=cH6*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
H(3,12)=DN(3,0)*Nenr[3] + DNenr(3,0)*cH10;
H(3,13)=DN(3,1)*Nenr[3] + DNenr(3,1)*cH10;
H(3,14)=DN(3,2)*Nenr[3] + DNenr(3,2)*cH10;
H(3,15)=cH6*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee4 = N[0]*vfrac(0,0) + N[1]*vfrac(1,0) + N[2]*vfrac(2,0) + N[3]*vfrac(3,0);
const double crhs_ee5 = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double crhs_ee6 = N[0]*vfrac(0,1) + N[1]*vfrac(1,1) + N[2]*vfrac(2,1) + N[3]*vfrac(3,1);
const double crhs_ee7 = DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0);
const double crhs_ee8 = N[0]*vfrac(0,2) + N[1]*vfrac(1,2) + N[2]*vfrac(2,2) + N[3]*vfrac(3,2);
const double crhs_ee9 = N[0]*bdf0;
const double crhs_ee10 = N[1]*bdf0;
const double crhs_ee11 = N[2]*bdf0;
const double crhs_ee12 = N[3]*bdf0;
const double crhs_ee13 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee14 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee15 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee16 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crhs_ee17 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crhs_ee18 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crhs_ee19 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee13, 2) + pow(crhs_ee14, 2) + pow(crhs_ee15, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee20 = crhs_ee19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(crhs_ee0*crhs_ee4 + crhs_ee5*crhs_ee6 + crhs_ee7*crhs_ee8) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*vn(0,0) + bdf1*vnn(0,0) + bdf2*vnnn(0,0)) + N[1]*(bdf0*vn(1,0) + bdf1*vnn(1,0) + bdf2*vnnn(1,0)) + N[2]*(bdf0*vn(2,0) + bdf1*vnn(2,0) + bdf2*vnnn(2,0)) + N[3]*(bdf0*vn(3,0) + bdf1*vnn(3,0) + bdf2*vnnn(3,0)) + crhs_ee16*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0)) + crhs_ee17*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0)) + crhs_ee18*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0))) + rho*(crhs_ee0*crhs_ee13 + crhs_ee10*(v(1,0) - vfrac(1,0)) + crhs_ee11*(v(2,0) - vfrac(2,0)) + crhs_ee12*(v(3,0) - vfrac(3,0)) + crhs_ee14*crhs_ee5 + crhs_ee15*crhs_ee7 + crhs_ee9*(v(0,0) - vfrac(0,0))));
const double crhs_ee21 = DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double crhs_ee22 = DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
const double crhs_ee23 = crhs_ee19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(crhs_ee1*crhs_ee6 + crhs_ee21*crhs_ee4 + crhs_ee22*crhs_ee8) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*vn(0,1) + bdf1*vnn(0,1) + bdf2*vnnn(0,1)) + N[1]*(bdf0*vn(1,1) + bdf1*vnn(1,1) + bdf2*vnnn(1,1)) + N[2]*(bdf0*vn(2,1) + bdf1*vnn(2,1) + bdf2*vnnn(2,1)) + N[3]*(bdf0*vn(3,1) + bdf1*vnn(3,1) + bdf2*vnnn(3,1)) + crhs_ee16*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1)) + crhs_ee17*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1)) + crhs_ee18*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1))) + rho*(crhs_ee1*crhs_ee14 + crhs_ee10*(v(1,1) - vfrac(1,1)) + crhs_ee11*(v(2,1) - vfrac(2,1)) + crhs_ee12*(v(3,1) - vfrac(3,1)) + crhs_ee13*crhs_ee21 + crhs_ee15*crhs_ee22 + crhs_ee9*(v(0,1) - vfrac(0,1))));
const double crhs_ee24 = DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
const double crhs_ee25 = DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
const double crhs_ee26 = crhs_ee19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(crhs_ee2*crhs_ee8 + crhs_ee24*crhs_ee4 + crhs_ee25*crhs_ee6) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*vn(0,2) + bdf1*vnn(0,2) + bdf2*vnnn(0,2)) + N[1]*(bdf0*vn(1,2) + bdf1*vnn(1,2) + bdf2*vnnn(1,2)) + N[2]*(bdf0*vn(2,2) + bdf1*vnn(2,2) + bdf2*vnnn(2,2)) + N[3]*(bdf0*vn(3,2) + bdf1*vnn(3,2) + bdf2*vnnn(3,2)) + crhs_ee16*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2)) + crhs_ee17*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2)) + crhs_ee18*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2))) + rho*(crhs_ee10*(v(1,2) - vfrac(1,2)) + crhs_ee11*(v(2,2) - vfrac(2,2)) + crhs_ee12*(v(3,2) - vfrac(3,2)) + crhs_ee13*crhs_ee24 + crhs_ee14*crhs_ee25 + crhs_ee15*crhs_ee2 + crhs_ee9*(v(0,2) - vfrac(0,2))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee20 - DNenr(0,1)*crhs_ee23 - DNenr(0,2)*crhs_ee26 - Nenr[0]*crhs_ee3;
rhs_ee[1]=-DNenr(1,0)*crhs_ee20 - DNenr(1,1)*crhs_ee23 - DNenr(1,2)*crhs_ee26 - Nenr[1]*crhs_ee3;
rhs_ee[2]=-DNenr(2,0)*crhs_ee20 - DNenr(2,1)*crhs_ee23 - DNenr(2,2)*crhs_ee26 - Nenr[2]*crhs_ee3;
rhs_ee[3]=-DNenr(3,0)*crhs_ee20 - DNenr(3,1)*crhs_ee23 - DNenr(3,2)*crhs_ee26 - Nenr[3]*crhs_ee3;


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

    // const auto v_convection = rData.Velocity - rData.Velocity_OldStep1;
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
