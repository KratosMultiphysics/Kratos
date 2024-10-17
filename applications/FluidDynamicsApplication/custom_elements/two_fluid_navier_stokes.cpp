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

    const auto vconv = rData.Velocity - rData.Velocity_OldStep1;
    const auto an = rData.Acceleration;

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
const double clhs11 = K_darcy*N[0];
const double clhs12 = N[0]*clhs10;
const double clhs13 = clhs11 + clhs12 + clhs9;
const double clhs14 = 1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 = 1.0*clhs9;
const double clhs16 = clhs14*clhs15;
const double clhs17 = 1.0*clhs11;
const double clhs18 = clhs14*clhs17;
const double clhs19 = 1.0*clhs14;
const double clhs20 = clhs13*clhs19;
const double clhs21 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double clhs22 = N[0]*clhs21;
const double clhs23 = K_darcy*clhs8 + N[0]*clhs9 + clhs10*clhs8 + clhs13*clhs16 - clhs13*clhs18 + clhs20*clhs22;
const double clhs24 = C(0,1)*DN(0,1) + clhs1;
const double clhs25 = C(1,2)*DN(0,1);
const double clhs26 = C(2,2)*DN(0,0) + clhs25;
const double clhs27 = DN(0,0)*clhs7;
const double clhs28 = DN(0,1)*clhs27;
const double clhs29 = clhs19*clhs21;
const double clhs30 = N[0]*clhs29 - N[0] + clhs16 - clhs18;
const double clhs31 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs32 = C(0,2)*DN(1,0);
const double clhs33 = C(2,2)*DN(1,1) + clhs32;
const double clhs34 = DN(0,0)*DN(1,0);
const double clhs35 = N[1]*clhs11 + N[1]*clhs12;
const double clhs36 = clhs34*clhs7 + clhs35;
const double clhs37 = rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs38 = K_darcy*N[1];
const double clhs39 = N[1]*clhs10;
const double clhs40 = clhs37 + clhs38 + clhs39;
const double clhs41 = clhs19*clhs40;
const double clhs42 = N[0]*clhs37 + clhs16*clhs40 - clhs18*clhs40 + clhs22*clhs41;
const double clhs43 = C(0,1)*DN(1,1) + clhs32;
const double clhs44 = C(1,2)*DN(1,1);
const double clhs45 = C(2,2)*DN(1,0) + clhs44;
const double clhs46 = DN(1,1)*clhs27;
const double clhs47 = DN(0,0)*N[1];
const double clhs48 = DN(1,0)*N[0];
const double clhs49 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs50 = C(0,2)*DN(2,0);
const double clhs51 = C(2,2)*DN(2,1) + clhs50;
const double clhs52 = DN(0,0)*DN(2,0);
const double clhs53 = N[2]*clhs11 + N[2]*clhs12;
const double clhs54 = clhs52*clhs7 + clhs53;
const double clhs55 = rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs56 = K_darcy*N[2];
const double clhs57 = N[2]*clhs10;
const double clhs58 = clhs55 + clhs56 + clhs57;
const double clhs59 = clhs19*clhs58;
const double clhs60 = N[0]*clhs55 + clhs16*clhs58 - clhs18*clhs58 + clhs22*clhs59;
const double clhs61 = C(0,1)*DN(2,1) + clhs50;
const double clhs62 = C(1,2)*DN(2,1);
const double clhs63 = C(2,2)*DN(2,0) + clhs62;
const double clhs64 = DN(2,1)*clhs27;
const double clhs65 = DN(0,0)*N[2];
const double clhs66 = DN(2,0)*N[0];
const double clhs67 = C(0,1)*DN(0,0) + clhs25;
const double clhs68 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs69 = pow(DN(0,1), 2);
const double clhs70 = C(0,1)*DN(1,0) + clhs44;
const double clhs71 = DN(0,1)*clhs7;
const double clhs72 = DN(1,0)*clhs71;
const double clhs73 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs74 = DN(0,1)*DN(1,1);
const double clhs75 = clhs35 + clhs7*clhs74;
const double clhs76 = DN(0,1)*N[1];
const double clhs77 = DN(1,1)*N[0];
const double clhs78 = C(0,1)*DN(2,0) + clhs62;
const double clhs79 = DN(2,0)*clhs71;
const double clhs80 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs81 = DN(0,1)*DN(2,1);
const double clhs82 = clhs53 + clhs7*clhs81;
const double clhs83 = DN(0,1)*N[2];
const double clhs84 = DN(2,1)*N[0];
const double clhs85 = N[0] + clhs14*(1.0*clhs12 + clhs15 + clhs17);
const double clhs86 = clhs19*(clhs34 + clhs74);
const double clhs87 = clhs19*(clhs52 + clhs81);
const double clhs88 = clhs19*clhs37;
const double clhs89 = clhs19*clhs38;
const double clhs90 = N[1]*clhs21;
const double clhs91 = N[1]*clhs9 + clhs13*clhs88 - clhs13*clhs89 + clhs20*clhs90;
const double clhs92 = pow(DN(1,0), 2);
const double clhs93 = pow(N[1], 2);
const double clhs94 = K_darcy*clhs93 + N[1]*clhs37 + clhs10*clhs93 + clhs37*clhs41 - clhs38*clhs41 + clhs41*clhs90;
const double clhs95 = DN(1,0)*clhs7;
const double clhs96 = DN(1,1)*clhs95;
const double clhs97 = N[1]*clhs29 - N[1] + clhs88 - clhs89;
const double clhs98 = DN(1,0)*DN(2,0);
const double clhs99 = N[2]*clhs38 + N[2]*clhs39;
const double clhs100 = clhs7*clhs98 + clhs99;
const double clhs101 = N[1]*clhs55 + clhs37*clhs59 - clhs38*clhs59 + clhs59*clhs90;
const double clhs102 = DN(2,1)*clhs95;
const double clhs103 = DN(1,0)*N[2];
const double clhs104 = DN(2,0)*N[1];
const double clhs105 = pow(DN(1,1), 2);
const double clhs106 = DN(2,0)*clhs7;
const double clhs107 = DN(1,1)*clhs106;
const double clhs108 = DN(1,1)*DN(2,1);
const double clhs109 = clhs108*clhs7 + clhs99;
const double clhs110 = DN(1,1)*N[2];
const double clhs111 = DN(2,1)*N[1];
const double clhs112 = N[1] + clhs14*(1.0*clhs37 + 1.0*clhs38 + 1.0*clhs39);
const double clhs113 = clhs19*(clhs108 + clhs98);
const double clhs114 = N[2]*clhs21;
const double clhs115 = N[2]*clhs9 + clhs114*clhs20 + clhs20*clhs55 - clhs20*clhs56;
const double clhs116 = clhs19*clhs56;
const double clhs117 = clhs19*clhs55;
const double clhs118 = N[2]*clhs37 + clhs114*clhs41 + clhs41*clhs55 - clhs41*clhs56;
const double clhs119 = pow(DN(2,0), 2);
const double clhs120 = pow(N[2], 2);
const double clhs121 = K_darcy*clhs120 + N[2]*clhs55 + clhs10*clhs120 + clhs114*clhs59 + clhs55*clhs59 - clhs56*clhs59;
const double clhs122 = DN(2,1)*clhs106;
const double clhs123 = N[2]*clhs29 - N[2] - clhs116 + clhs117;
const double clhs124 = pow(DN(2,1), 2);
const double clhs125 = N[2] + clhs14*(1.0*clhs55 + 1.0*clhs56 + 1.0*clhs57);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs23 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs28;
lhs(0,2)=DN(0,0)*clhs30;
lhs(0,3)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + clhs36 + clhs42;
lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46;
lhs(0,5)=DN(1,0)*clhs16 - DN(1,0)*clhs18 + clhs29*clhs48 - clhs47;
lhs(0,6)=DN(0,0)*clhs49 + DN(0,1)*clhs51 + clhs54 + clhs60;
lhs(0,7)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + clhs64;
lhs(0,8)=DN(2,0)*clhs16 - DN(2,0)*clhs18 + clhs29*clhs66 - clhs65;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs67 + clhs28;
lhs(1,1)=DN(0,0)*clhs26 + DN(0,1)*clhs68 + clhs23 + clhs69*clhs7;
lhs(1,2)=DN(0,1)*clhs30;
lhs(1,3)=DN(0,0)*clhs33 + DN(0,1)*clhs70 + clhs72;
lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs73 + clhs42 + clhs75;
lhs(1,5)=DN(1,1)*clhs16 - DN(1,1)*clhs18 + clhs29*clhs77 - clhs76;
lhs(1,6)=DN(0,0)*clhs51 + DN(0,1)*clhs78 + clhs79;
lhs(1,7)=DN(0,0)*clhs63 + DN(0,1)*clhs80 + clhs60 + clhs82;
lhs(1,8)=DN(2,1)*clhs16 - DN(2,1)*clhs18 + clhs29*clhs84 - clhs83;
lhs(2,0)=DN(0,0)*clhs85;
lhs(2,1)=DN(0,1)*clhs85;
lhs(2,2)=clhs19*(clhs3 + clhs69);
lhs(2,3)=DN(0,0)*clhs41 + clhs48;
lhs(2,4)=DN(0,1)*clhs41 + clhs77;
lhs(2,5)=clhs86;
lhs(2,6)=DN(0,0)*clhs59 + clhs66;
lhs(2,7)=DN(0,1)*clhs59 + clhs84;
lhs(2,8)=clhs87;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs36 + clhs91;
lhs(3,1)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs72;
lhs(3,2)=DN(0,0)*clhs88 - DN(0,0)*clhs89 + clhs29*clhs47 - clhs48;
lhs(3,3)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + clhs7*clhs92 + clhs94;
lhs(3,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs96;
lhs(3,5)=DN(1,0)*clhs97;
lhs(3,6)=DN(1,0)*clhs49 + DN(1,1)*clhs51 + clhs100 + clhs101;
lhs(3,7)=DN(1,0)*clhs61 + DN(1,1)*clhs63 + clhs102;
lhs(3,8)=DN(2,0)*clhs88 - DN(2,0)*clhs89 - clhs103 + clhs104*clhs29;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs67 + clhs46;
lhs(4,1)=DN(1,0)*clhs26 + DN(1,1)*clhs68 + clhs75 + clhs91;
lhs(4,2)=DN(0,1)*clhs88 - DN(0,1)*clhs89 + clhs29*clhs76 - clhs77;
lhs(4,3)=DN(1,0)*clhs33 + DN(1,1)*clhs70 + clhs96;
lhs(4,4)=DN(1,0)*clhs45 + DN(1,1)*clhs73 + clhs105*clhs7 + clhs94;
lhs(4,5)=DN(1,1)*clhs97;
lhs(4,6)=DN(1,0)*clhs51 + DN(1,1)*clhs78 + clhs107;
lhs(4,7)=DN(1,0)*clhs63 + DN(1,1)*clhs80 + clhs101 + clhs109;
lhs(4,8)=DN(2,1)*clhs88 - DN(2,1)*clhs89 - clhs110 + clhs111*clhs29;
lhs(5,0)=DN(1,0)*clhs20 + clhs47;
lhs(5,1)=DN(1,1)*clhs20 + clhs76;
lhs(5,2)=clhs86;
lhs(5,3)=DN(1,0)*clhs112;
lhs(5,4)=DN(1,1)*clhs112;
lhs(5,5)=clhs19*(clhs105 + clhs92);
lhs(5,6)=DN(1,0)*clhs59 + clhs104;
lhs(5,7)=DN(1,1)*clhs59 + clhs111;
lhs(5,8)=clhs113;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs115 + clhs54;
lhs(6,1)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs79;
lhs(6,2)=-DN(0,0)*clhs116 + DN(0,0)*clhs117 + clhs29*clhs65 - clhs66;
lhs(6,3)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + clhs100 + clhs118;
lhs(6,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs107;
lhs(6,5)=-DN(1,0)*clhs116 + DN(1,0)*clhs117 + clhs103*clhs29 - clhs104;
lhs(6,6)=DN(2,0)*clhs49 + DN(2,1)*clhs51 + clhs119*clhs7 + clhs121;
lhs(6,7)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + clhs122;
lhs(6,8)=DN(2,0)*clhs123;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs67 + clhs64;
lhs(7,1)=DN(2,0)*clhs26 + DN(2,1)*clhs68 + clhs115 + clhs82;
lhs(7,2)=-DN(0,1)*clhs116 + DN(0,1)*clhs117 + clhs29*clhs83 - clhs84;
lhs(7,3)=DN(2,0)*clhs33 + DN(2,1)*clhs70 + clhs102;
lhs(7,4)=DN(2,0)*clhs45 + DN(2,1)*clhs73 + clhs109 + clhs118;
lhs(7,5)=-DN(1,1)*clhs116 + DN(1,1)*clhs117 + clhs110*clhs29 - clhs111;
lhs(7,6)=DN(2,0)*clhs51 + DN(2,1)*clhs78 + clhs122;
lhs(7,7)=DN(2,0)*clhs63 + DN(2,1)*clhs80 + clhs121 + clhs124*clhs7;
lhs(7,8)=DN(2,1)*clhs123;
lhs(8,0)=DN(2,0)*clhs20 + clhs65;
lhs(8,1)=DN(2,1)*clhs20 + clhs83;
lhs(8,2)=clhs87;
lhs(8,3)=DN(2,0)*clhs41 + clhs103;
lhs(8,4)=DN(2,1)*clhs41 + clhs110;
lhs(8,5)=clhs113;
lhs(8,6)=DN(2,0)*clhs125;
lhs(8,7)=DN(2,1)*clhs125;
lhs(8,8)=clhs19*(clhs119 + clhs124);


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

    const auto vconv = rData.Velocity - rData.Velocity_OldStep1;
    const auto an = rData.Acceleration;

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
const double clhs14 = K_darcy*N[0];
const double clhs15 = N[0]*clhs13;
const double clhs16 = clhs12 + clhs14 + clhs15;
const double clhs17 = 1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = 1.0*clhs12;
const double clhs19 = clhs17*clhs18;
const double clhs20 = 1.0*clhs14;
const double clhs21 = clhs17*clhs20;
const double clhs22 = 1.0*clhs17;
const double clhs23 = clhs16*clhs22;
const double clhs24 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double clhs25 = N[0]*clhs24;
const double clhs26 = K_darcy*clhs11 + N[0]*clhs12 + clhs11*clhs13 + clhs16*clhs19 - clhs16*clhs21 + clhs23*clhs25;
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
const double clhs41 = clhs22*clhs24;
const double clhs42 = N[0]*clhs41 - N[0] + clhs19 - clhs21;
const double clhs43 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs44 = C(0,3)*DN(1,0);
const double clhs45 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs44;
const double clhs46 = C(0,5)*DN(1,0);
const double clhs47 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs46;
const double clhs48 = DN(0,0)*DN(1,0);
const double clhs49 = N[1]*clhs14 + N[1]*clhs15;
const double clhs50 = clhs10*clhs48 + clhs49;
const double clhs51 = rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs52 = K_darcy*N[1];
const double clhs53 = N[1]*clhs13;
const double clhs54 = clhs51 + clhs52 + clhs53;
const double clhs55 = clhs22*clhs54;
const double clhs56 = N[0]*clhs51 + clhs19*clhs54 - clhs21*clhs54 + clhs25*clhs55;
const double clhs57 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs44;
const double clhs58 = C(1,3)*DN(1,1);
const double clhs59 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs58;
const double clhs60 = C(3,5)*DN(1,0);
const double clhs61 = C(4,5)*DN(1,2);
const double clhs62 = C(1,5)*DN(1,1) + clhs60 + clhs61;
const double clhs63 = DN(1,1)*clhs33;
const double clhs64 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs46;
const double clhs65 = C(3,4)*DN(1,1);
const double clhs66 = C(2,3)*DN(1,2) + clhs60 + clhs65;
const double clhs67 = C(2,5)*DN(1,2);
const double clhs68 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs67;
const double clhs69 = DN(1,2)*clhs33;
const double clhs70 = DN(0,0)*N[1];
const double clhs71 = DN(1,0)*N[0];
const double clhs72 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs73 = C(0,3)*DN(2,0);
const double clhs74 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs73;
const double clhs75 = C(0,5)*DN(2,0);
const double clhs76 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs75;
const double clhs77 = DN(0,0)*DN(2,0);
const double clhs78 = N[2]*clhs14 + N[2]*clhs15;
const double clhs79 = clhs10*clhs77 + clhs78;
const double clhs80 = rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs81 = K_darcy*N[2];
const double clhs82 = N[2]*clhs13;
const double clhs83 = clhs80 + clhs81 + clhs82;
const double clhs84 = clhs22*clhs83;
const double clhs85 = N[0]*clhs80 + clhs19*clhs83 - clhs21*clhs83 + clhs25*clhs84;
const double clhs86 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs73;
const double clhs87 = C(1,3)*DN(2,1);
const double clhs88 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs87;
const double clhs89 = C(3,5)*DN(2,0);
const double clhs90 = C(4,5)*DN(2,2);
const double clhs91 = C(1,5)*DN(2,1) + clhs89 + clhs90;
const double clhs92 = DN(2,1)*clhs33;
const double clhs93 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs75;
const double clhs94 = C(3,4)*DN(2,1);
const double clhs95 = C(2,3)*DN(2,2) + clhs89 + clhs94;
const double clhs96 = C(2,5)*DN(2,2);
const double clhs97 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs96;
const double clhs98 = DN(2,2)*clhs33;
const double clhs99 = DN(0,0)*N[2];
const double clhs100 = DN(2,0)*N[0];
const double clhs101 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs102 = C(0,3)*DN(3,0);
const double clhs103 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs102;
const double clhs104 = C(0,5)*DN(3,0);
const double clhs105 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs104;
const double clhs106 = DN(0,0)*DN(3,0);
const double clhs107 = N[3]*clhs14 + N[3]*clhs15;
const double clhs108 = clhs10*clhs106 + clhs107;
const double clhs109 = rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs110 = K_darcy*N[3];
const double clhs111 = N[3]*clhs13;
const double clhs112 = clhs109 + clhs110 + clhs111;
const double clhs113 = clhs112*clhs22;
const double clhs114 = N[0]*clhs109 + clhs112*clhs19 - clhs112*clhs21 + clhs113*clhs25;
const double clhs115 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs102;
const double clhs116 = C(1,3)*DN(3,1);
const double clhs117 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs116;
const double clhs118 = C(3,5)*DN(3,0);
const double clhs119 = C(4,5)*DN(3,2);
const double clhs120 = C(1,5)*DN(3,1) + clhs118 + clhs119;
const double clhs121 = DN(3,1)*clhs33;
const double clhs122 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs104;
const double clhs123 = C(3,4)*DN(3,1);
const double clhs124 = C(2,3)*DN(3,2) + clhs118 + clhs123;
const double clhs125 = C(2,5)*DN(3,2);
const double clhs126 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs125;
const double clhs127 = DN(3,2)*clhs33;
const double clhs128 = DN(0,0)*N[3];
const double clhs129 = DN(3,0)*N[0];
const double clhs130 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs28;
const double clhs131 = C(0,4)*DN(0,0) + clhs31 + clhs36;
const double clhs132 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs133 = C(1,4)*DN(0,1);
const double clhs134 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs133;
const double clhs135 = pow(DN(0,1), 2);
const double clhs136 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs133;
const double clhs137 = C(2,4)*DN(0,2);
const double clhs138 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs137;
const double clhs139 = DN(0,1)*clhs10;
const double clhs140 = DN(0,2)*clhs139;
const double clhs141 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs58;
const double clhs142 = C(0,4)*DN(1,0) + clhs61 + clhs65;
const double clhs143 = DN(1,0)*clhs139;
const double clhs144 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs145 = C(1,4)*DN(1,1);
const double clhs146 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs145;
const double clhs147 = DN(0,1)*DN(1,1);
const double clhs148 = clhs10*clhs147;
const double clhs149 = clhs49 + clhs56;
const double clhs150 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs145;
const double clhs151 = C(2,4)*DN(1,2);
const double clhs152 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs151;
const double clhs153 = DN(1,2)*clhs139;
const double clhs154 = DN(0,1)*N[1];
const double clhs155 = DN(1,1)*N[0];
const double clhs156 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs87;
const double clhs157 = C(0,4)*DN(2,0) + clhs90 + clhs94;
const double clhs158 = DN(2,0)*clhs139;
const double clhs159 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs160 = C(1,4)*DN(2,1);
const double clhs161 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs160;
const double clhs162 = DN(0,1)*DN(2,1);
const double clhs163 = clhs10*clhs162;
const double clhs164 = clhs78 + clhs85;
const double clhs165 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs160;
const double clhs166 = C(2,4)*DN(2,2);
const double clhs167 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs166;
const double clhs168 = DN(2,2)*clhs139;
const double clhs169 = DN(0,1)*N[2];
const double clhs170 = DN(2,1)*N[0];
const double clhs171 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs116;
const double clhs172 = C(0,4)*DN(3,0) + clhs119 + clhs123;
const double clhs173 = DN(3,0)*clhs139;
const double clhs174 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs175 = C(1,4)*DN(3,1);
const double clhs176 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs175;
const double clhs177 = DN(0,1)*DN(3,1);
const double clhs178 = clhs10*clhs177;
const double clhs179 = clhs107 + clhs114;
const double clhs180 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs175;
const double clhs181 = C(2,4)*DN(3,2);
const double clhs182 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs181;
const double clhs183 = DN(3,2)*clhs139;
const double clhs184 = DN(0,1)*N[3];
const double clhs185 = DN(3,1)*N[0];
const double clhs186 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs38;
const double clhs187 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs137;
const double clhs188 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs189 = pow(DN(0,2), 2);
const double clhs190 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs67;
const double clhs191 = DN(0,2)*clhs10;
const double clhs192 = DN(1,0)*clhs191;
const double clhs193 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs151;
const double clhs194 = DN(1,1)*clhs191;
const double clhs195 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs196 = DN(0,2)*DN(1,2);
const double clhs197 = clhs10*clhs196;
const double clhs198 = DN(0,2)*N[1];
const double clhs199 = DN(1,2)*N[0];
const double clhs200 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs96;
const double clhs201 = DN(2,0)*clhs191;
const double clhs202 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs166;
const double clhs203 = DN(2,1)*clhs191;
const double clhs204 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs205 = DN(0,2)*DN(2,2);
const double clhs206 = clhs10*clhs205;
const double clhs207 = DN(0,2)*N[2];
const double clhs208 = DN(2,2)*N[0];
const double clhs209 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs125;
const double clhs210 = DN(3,0)*clhs191;
const double clhs211 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs181;
const double clhs212 = DN(3,1)*clhs191;
const double clhs213 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs214 = DN(0,2)*DN(3,2);
const double clhs215 = clhs10*clhs214;
const double clhs216 = DN(0,2)*N[3];
const double clhs217 = DN(3,2)*N[0];
const double clhs218 = N[0] + clhs17*(1.0*clhs15 + clhs18 + clhs20);
const double clhs219 = clhs22*(clhs147 + clhs196 + clhs48);
const double clhs220 = clhs22*(clhs162 + clhs205 + clhs77);
const double clhs221 = clhs22*(clhs106 + clhs177 + clhs214);
const double clhs222 = clhs22*clhs51;
const double clhs223 = clhs22*clhs52;
const double clhs224 = N[1]*clhs24;
const double clhs225 = N[1]*clhs12 + clhs16*clhs222 - clhs16*clhs223 + clhs224*clhs23;
const double clhs226 = pow(DN(1,0), 2);
const double clhs227 = pow(N[1], 2);
const double clhs228 = K_darcy*clhs227 + N[1]*clhs51 + clhs13*clhs227 + clhs224*clhs55 + clhs51*clhs55 - clhs52*clhs55;
const double clhs229 = DN(1,0)*clhs10;
const double clhs230 = DN(1,1)*clhs229;
const double clhs231 = DN(1,2)*clhs229;
const double clhs232 = N[1]*clhs41 - N[1] + clhs222 - clhs223;
const double clhs233 = DN(1,0)*DN(2,0);
const double clhs234 = N[2]*clhs52 + N[2]*clhs53;
const double clhs235 = clhs10*clhs233 + clhs234;
const double clhs236 = N[1]*clhs80 + clhs224*clhs84 + clhs51*clhs84 - clhs52*clhs84;
const double clhs237 = DN(2,1)*clhs229;
const double clhs238 = DN(2,2)*clhs229;
const double clhs239 = DN(1,0)*N[2];
const double clhs240 = DN(2,0)*N[1];
const double clhs241 = DN(1,0)*DN(3,0);
const double clhs242 = N[3]*clhs52 + N[3]*clhs53;
const double clhs243 = clhs10*clhs241 + clhs242;
const double clhs244 = N[1]*clhs109 + clhs113*clhs224 + clhs113*clhs51 - clhs113*clhs52;
const double clhs245 = DN(3,1)*clhs229;
const double clhs246 = DN(3,2)*clhs229;
const double clhs247 = DN(1,0)*N[3];
const double clhs248 = DN(3,0)*N[1];
const double clhs249 = clhs225 + clhs49;
const double clhs250 = pow(DN(1,1), 2);
const double clhs251 = DN(1,1)*clhs10;
const double clhs252 = DN(1,2)*clhs251;
const double clhs253 = DN(2,0)*clhs251;
const double clhs254 = DN(1,1)*DN(2,1);
const double clhs255 = clhs10*clhs254;
const double clhs256 = clhs234 + clhs236;
const double clhs257 = DN(2,2)*clhs251;
const double clhs258 = DN(1,1)*N[2];
const double clhs259 = DN(2,1)*N[1];
const double clhs260 = DN(3,0)*clhs251;
const double clhs261 = DN(1,1)*DN(3,1);
const double clhs262 = clhs10*clhs261;
const double clhs263 = clhs242 + clhs244;
const double clhs264 = DN(3,2)*clhs251;
const double clhs265 = DN(1,1)*N[3];
const double clhs266 = DN(3,1)*N[1];
const double clhs267 = pow(DN(1,2), 2);
const double clhs268 = DN(1,2)*clhs10;
const double clhs269 = DN(2,0)*clhs268;
const double clhs270 = DN(2,1)*clhs268;
const double clhs271 = DN(1,2)*DN(2,2);
const double clhs272 = clhs10*clhs271;
const double clhs273 = DN(1,2)*N[2];
const double clhs274 = DN(2,2)*N[1];
const double clhs275 = DN(3,0)*clhs268;
const double clhs276 = DN(3,1)*clhs268;
const double clhs277 = DN(1,2)*DN(3,2);
const double clhs278 = clhs10*clhs277;
const double clhs279 = DN(1,2)*N[3];
const double clhs280 = DN(3,2)*N[1];
const double clhs281 = N[1] + clhs17*(1.0*clhs51 + 1.0*clhs52 + 1.0*clhs53);
const double clhs282 = clhs22*(clhs233 + clhs254 + clhs271);
const double clhs283 = clhs22*(clhs241 + clhs261 + clhs277);
const double clhs284 = N[2]*clhs24;
const double clhs285 = N[2]*clhs12 + clhs23*clhs284 + clhs23*clhs80 - clhs23*clhs81;
const double clhs286 = clhs22*clhs81;
const double clhs287 = clhs22*clhs80;
const double clhs288 = N[2]*clhs51 + clhs284*clhs55 + clhs55*clhs80 - clhs55*clhs81;
const double clhs289 = pow(DN(2,0), 2);
const double clhs290 = pow(N[2], 2);
const double clhs291 = K_darcy*clhs290 + N[2]*clhs80 + clhs13*clhs290 + clhs284*clhs84 + clhs80*clhs84 - clhs81*clhs84;
const double clhs292 = DN(2,0)*clhs10;
const double clhs293 = DN(2,1)*clhs292;
const double clhs294 = DN(2,2)*clhs292;
const double clhs295 = N[2]*clhs41 - N[2] - clhs286 + clhs287;
const double clhs296 = DN(2,0)*DN(3,0);
const double clhs297 = N[3]*clhs81 + N[3]*clhs82;
const double clhs298 = clhs10*clhs296 + clhs297;
const double clhs299 = N[2]*clhs109 + clhs113*clhs284 + clhs113*clhs80 - clhs113*clhs81;
const double clhs300 = DN(3,1)*clhs292;
const double clhs301 = DN(3,2)*clhs292;
const double clhs302 = DN(2,0)*N[3];
const double clhs303 = DN(3,0)*N[2];
const double clhs304 = clhs285 + clhs78;
const double clhs305 = clhs234 + clhs288;
const double clhs306 = pow(DN(2,1), 2);
const double clhs307 = DN(2,1)*clhs10;
const double clhs308 = DN(2,2)*clhs307;
const double clhs309 = DN(3,0)*clhs307;
const double clhs310 = DN(2,1)*DN(3,1);
const double clhs311 = clhs10*clhs310;
const double clhs312 = clhs297 + clhs299;
const double clhs313 = DN(3,2)*clhs307;
const double clhs314 = DN(2,1)*N[3];
const double clhs315 = DN(3,1)*N[2];
const double clhs316 = pow(DN(2,2), 2);
const double clhs317 = DN(2,2)*clhs10;
const double clhs318 = DN(3,0)*clhs317;
const double clhs319 = DN(3,1)*clhs317;
const double clhs320 = DN(2,2)*DN(3,2);
const double clhs321 = clhs10*clhs320;
const double clhs322 = DN(2,2)*N[3];
const double clhs323 = DN(3,2)*N[2];
const double clhs324 = N[2] + clhs17*(1.0*clhs80 + 1.0*clhs81 + 1.0*clhs82);
const double clhs325 = clhs22*(clhs296 + clhs310 + clhs320);
const double clhs326 = N[3]*clhs24;
const double clhs327 = N[3]*clhs12 + clhs109*clhs23 - clhs110*clhs23 + clhs23*clhs326;
const double clhs328 = clhs110*clhs22;
const double clhs329 = clhs109*clhs22;
const double clhs330 = N[3]*clhs51 + clhs109*clhs55 - clhs110*clhs55 + clhs326*clhs55;
const double clhs331 = N[3]*clhs80 + clhs109*clhs84 - clhs110*clhs84 + clhs326*clhs84;
const double clhs332 = pow(DN(3,0), 2);
const double clhs333 = pow(N[3], 2);
const double clhs334 = K_darcy*clhs333 + N[3]*clhs109 + clhs109*clhs113 - clhs110*clhs113 + clhs113*clhs326 + clhs13*clhs333;
const double clhs335 = DN(3,0)*clhs10;
const double clhs336 = DN(3,1)*clhs335;
const double clhs337 = DN(3,2)*clhs335;
const double clhs338 = N[3]*clhs41 - N[3] - clhs328 + clhs329;
const double clhs339 = clhs107 + clhs327;
const double clhs340 = clhs242 + clhs330;
const double clhs341 = clhs297 + clhs331;
const double clhs342 = pow(DN(3,1), 2);
const double clhs343 = DN(3,1)*DN(3,2)*clhs10;
const double clhs344 = pow(DN(3,2), 2);
const double clhs345 = N[3] + clhs17*(1.0*clhs109 + 1.0*clhs110 + 1.0*clhs111);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs26;
lhs(0,1)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + DN(0,2)*clhs32 + clhs34;
lhs(0,2)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + DN(0,2)*clhs39 + clhs40;
lhs(0,3)=DN(0,0)*clhs42;
lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + DN(0,2)*clhs47 + clhs50 + clhs56;
lhs(0,5)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + DN(0,2)*clhs62 + clhs63;
lhs(0,6)=DN(0,0)*clhs64 + DN(0,1)*clhs66 + DN(0,2)*clhs68 + clhs69;
lhs(0,7)=DN(1,0)*clhs19 - DN(1,0)*clhs21 + clhs41*clhs71 - clhs70;
lhs(0,8)=DN(0,0)*clhs72 + DN(0,1)*clhs74 + DN(0,2)*clhs76 + clhs79 + clhs85;
lhs(0,9)=DN(0,0)*clhs86 + DN(0,1)*clhs88 + DN(0,2)*clhs91 + clhs92;
lhs(0,10)=DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs97 + clhs98;
lhs(0,11)=DN(2,0)*clhs19 - DN(2,0)*clhs21 + clhs100*clhs41 - clhs99;
lhs(0,12)=DN(0,0)*clhs101 + DN(0,1)*clhs103 + DN(0,2)*clhs105 + clhs108 + clhs114;
lhs(0,13)=DN(0,0)*clhs115 + DN(0,1)*clhs117 + DN(0,2)*clhs120 + clhs121;
lhs(0,14)=DN(0,0)*clhs122 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs127;
lhs(0,15)=DN(3,0)*clhs19 - DN(3,0)*clhs21 - clhs128 + clhs129*clhs41;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs34;
lhs(1,1)=DN(0,0)*clhs29 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs10*clhs135 + clhs26;
lhs(1,2)=DN(0,0)*clhs37 + DN(0,1)*clhs136 + DN(0,2)*clhs138 + clhs140;
lhs(1,3)=DN(0,1)*clhs42;
lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs141 + DN(0,2)*clhs142 + clhs143;
lhs(1,5)=DN(0,0)*clhs59 + DN(0,1)*clhs144 + DN(0,2)*clhs146 + clhs148 + clhs149;
lhs(1,6)=DN(0,0)*clhs66 + DN(0,1)*clhs150 + DN(0,2)*clhs152 + clhs153;
lhs(1,7)=DN(1,1)*clhs19 - DN(1,1)*clhs21 - clhs154 + clhs155*clhs41;
lhs(1,8)=DN(0,0)*clhs74 + DN(0,1)*clhs156 + DN(0,2)*clhs157 + clhs158;
lhs(1,9)=DN(0,0)*clhs88 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs163 + clhs164;
lhs(1,10)=DN(0,0)*clhs95 + DN(0,1)*clhs165 + DN(0,2)*clhs167 + clhs168;
lhs(1,11)=DN(2,1)*clhs19 - DN(2,1)*clhs21 - clhs169 + clhs170*clhs41;
lhs(1,12)=DN(0,0)*clhs103 + DN(0,1)*clhs171 + DN(0,2)*clhs172 + clhs173;
lhs(1,13)=DN(0,0)*clhs117 + DN(0,1)*clhs174 + DN(0,2)*clhs176 + clhs178 + clhs179;
lhs(1,14)=DN(0,0)*clhs124 + DN(0,1)*clhs180 + DN(0,2)*clhs182 + clhs183;
lhs(1,15)=DN(3,1)*clhs19 - DN(3,1)*clhs21 - clhs184 + clhs185*clhs41;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs131 + DN(0,2)*clhs186 + clhs40;
lhs(2,1)=DN(0,0)*clhs32 + DN(0,1)*clhs134 + DN(0,2)*clhs187 + clhs140;
lhs(2,2)=DN(0,0)*clhs39 + DN(0,1)*clhs138 + DN(0,2)*clhs188 + clhs10*clhs189 + clhs26;
lhs(2,3)=DN(0,2)*clhs42;
lhs(2,4)=DN(0,0)*clhs47 + DN(0,1)*clhs142 + DN(0,2)*clhs190 + clhs192;
lhs(2,5)=DN(0,0)*clhs62 + DN(0,1)*clhs146 + DN(0,2)*clhs193 + clhs194;
lhs(2,6)=DN(0,0)*clhs68 + DN(0,1)*clhs152 + DN(0,2)*clhs195 + clhs149 + clhs197;
lhs(2,7)=DN(1,2)*clhs19 - DN(1,2)*clhs21 - clhs198 + clhs199*clhs41;
lhs(2,8)=DN(0,0)*clhs76 + DN(0,1)*clhs157 + DN(0,2)*clhs200 + clhs201;
lhs(2,9)=DN(0,0)*clhs91 + DN(0,1)*clhs161 + DN(0,2)*clhs202 + clhs203;
lhs(2,10)=DN(0,0)*clhs97 + DN(0,1)*clhs167 + DN(0,2)*clhs204 + clhs164 + clhs206;
lhs(2,11)=DN(2,2)*clhs19 - DN(2,2)*clhs21 - clhs207 + clhs208*clhs41;
lhs(2,12)=DN(0,0)*clhs105 + DN(0,1)*clhs172 + DN(0,2)*clhs209 + clhs210;
lhs(2,13)=DN(0,0)*clhs120 + DN(0,1)*clhs176 + DN(0,2)*clhs211 + clhs212;
lhs(2,14)=DN(0,0)*clhs126 + DN(0,1)*clhs182 + DN(0,2)*clhs213 + clhs179 + clhs215;
lhs(2,15)=DN(3,2)*clhs19 - DN(3,2)*clhs21 - clhs216 + clhs217*clhs41;
lhs(3,0)=DN(0,0)*clhs218;
lhs(3,1)=DN(0,1)*clhs218;
lhs(3,2)=DN(0,2)*clhs218;
lhs(3,3)=clhs22*(clhs135 + clhs189 + clhs5);
lhs(3,4)=DN(0,0)*clhs55 + clhs71;
lhs(3,5)=DN(0,1)*clhs55 + clhs155;
lhs(3,6)=DN(0,2)*clhs55 + clhs199;
lhs(3,7)=clhs219;
lhs(3,8)=DN(0,0)*clhs84 + clhs100;
lhs(3,9)=DN(0,1)*clhs84 + clhs170;
lhs(3,10)=DN(0,2)*clhs84 + clhs208;
lhs(3,11)=clhs220;
lhs(3,12)=DN(0,0)*clhs113 + clhs129;
lhs(3,13)=DN(0,1)*clhs113 + clhs185;
lhs(3,14)=DN(0,2)*clhs113 + clhs217;
lhs(3,15)=clhs221;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs225 + clhs50;
lhs(4,1)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + DN(1,2)*clhs32 + clhs143;
lhs(4,2)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + DN(1,2)*clhs39 + clhs192;
lhs(4,3)=DN(0,0)*clhs222 - DN(0,0)*clhs223 + clhs41*clhs70 - clhs71;
lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + DN(1,2)*clhs47 + clhs10*clhs226 + clhs228;
lhs(4,5)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + DN(1,2)*clhs62 + clhs230;
lhs(4,6)=DN(1,0)*clhs64 + DN(1,1)*clhs66 + DN(1,2)*clhs68 + clhs231;
lhs(4,7)=DN(1,0)*clhs232;
lhs(4,8)=DN(1,0)*clhs72 + DN(1,1)*clhs74 + DN(1,2)*clhs76 + clhs235 + clhs236;
lhs(4,9)=DN(1,0)*clhs86 + DN(1,1)*clhs88 + DN(1,2)*clhs91 + clhs237;
lhs(4,10)=DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs97 + clhs238;
lhs(4,11)=DN(2,0)*clhs222 - DN(2,0)*clhs223 - clhs239 + clhs240*clhs41;
lhs(4,12)=DN(1,0)*clhs101 + DN(1,1)*clhs103 + DN(1,2)*clhs105 + clhs243 + clhs244;
lhs(4,13)=DN(1,0)*clhs115 + DN(1,1)*clhs117 + DN(1,2)*clhs120 + clhs245;
lhs(4,14)=DN(1,0)*clhs122 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs246;
lhs(4,15)=DN(3,0)*clhs222 - DN(3,0)*clhs223 - clhs247 + clhs248*clhs41;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs63;
lhs(5,1)=DN(1,0)*clhs29 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs148 + clhs249;
lhs(5,2)=DN(1,0)*clhs37 + DN(1,1)*clhs136 + DN(1,2)*clhs138 + clhs194;
lhs(5,3)=DN(0,1)*clhs222 - DN(0,1)*clhs223 + clhs154*clhs41 - clhs155;
lhs(5,4)=DN(1,0)*clhs45 + DN(1,1)*clhs141 + DN(1,2)*clhs142 + clhs230;
lhs(5,5)=DN(1,0)*clhs59 + DN(1,1)*clhs144 + DN(1,2)*clhs146 + clhs10*clhs250 + clhs228;
lhs(5,6)=DN(1,0)*clhs66 + DN(1,1)*clhs150 + DN(1,2)*clhs152 + clhs252;
lhs(5,7)=DN(1,1)*clhs232;
lhs(5,8)=DN(1,0)*clhs74 + DN(1,1)*clhs156 + DN(1,2)*clhs157 + clhs253;
lhs(5,9)=DN(1,0)*clhs88 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs255 + clhs256;
lhs(5,10)=DN(1,0)*clhs95 + DN(1,1)*clhs165 + DN(1,2)*clhs167 + clhs257;
lhs(5,11)=DN(2,1)*clhs222 - DN(2,1)*clhs223 - clhs258 + clhs259*clhs41;
lhs(5,12)=DN(1,0)*clhs103 + DN(1,1)*clhs171 + DN(1,2)*clhs172 + clhs260;
lhs(5,13)=DN(1,0)*clhs117 + DN(1,1)*clhs174 + DN(1,2)*clhs176 + clhs262 + clhs263;
lhs(5,14)=DN(1,0)*clhs124 + DN(1,1)*clhs180 + DN(1,2)*clhs182 + clhs264;
lhs(5,15)=DN(3,1)*clhs222 - DN(3,1)*clhs223 - clhs265 + clhs266*clhs41;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs131 + DN(1,2)*clhs186 + clhs69;
lhs(6,1)=DN(1,0)*clhs32 + DN(1,1)*clhs134 + DN(1,2)*clhs187 + clhs153;
lhs(6,2)=DN(1,0)*clhs39 + DN(1,1)*clhs138 + DN(1,2)*clhs188 + clhs197 + clhs249;
lhs(6,3)=DN(0,2)*clhs222 - DN(0,2)*clhs223 + clhs198*clhs41 - clhs199;
lhs(6,4)=DN(1,0)*clhs47 + DN(1,1)*clhs142 + DN(1,2)*clhs190 + clhs231;
lhs(6,5)=DN(1,0)*clhs62 + DN(1,1)*clhs146 + DN(1,2)*clhs193 + clhs252;
lhs(6,6)=DN(1,0)*clhs68 + DN(1,1)*clhs152 + DN(1,2)*clhs195 + clhs10*clhs267 + clhs228;
lhs(6,7)=DN(1,2)*clhs232;
lhs(6,8)=DN(1,0)*clhs76 + DN(1,1)*clhs157 + DN(1,2)*clhs200 + clhs269;
lhs(6,9)=DN(1,0)*clhs91 + DN(1,1)*clhs161 + DN(1,2)*clhs202 + clhs270;
lhs(6,10)=DN(1,0)*clhs97 + DN(1,1)*clhs167 + DN(1,2)*clhs204 + clhs256 + clhs272;
lhs(6,11)=DN(2,2)*clhs222 - DN(2,2)*clhs223 - clhs273 + clhs274*clhs41;
lhs(6,12)=DN(1,0)*clhs105 + DN(1,1)*clhs172 + DN(1,2)*clhs209 + clhs275;
lhs(6,13)=DN(1,0)*clhs120 + DN(1,1)*clhs176 + DN(1,2)*clhs211 + clhs276;
lhs(6,14)=DN(1,0)*clhs126 + DN(1,1)*clhs182 + DN(1,2)*clhs213 + clhs263 + clhs278;
lhs(6,15)=DN(3,2)*clhs222 - DN(3,2)*clhs223 - clhs279 + clhs280*clhs41;
lhs(7,0)=DN(1,0)*clhs23 + clhs70;
lhs(7,1)=DN(1,1)*clhs23 + clhs154;
lhs(7,2)=DN(1,2)*clhs23 + clhs198;
lhs(7,3)=clhs219;
lhs(7,4)=DN(1,0)*clhs281;
lhs(7,5)=DN(1,1)*clhs281;
lhs(7,6)=DN(1,2)*clhs281;
lhs(7,7)=clhs22*(clhs226 + clhs250 + clhs267);
lhs(7,8)=DN(1,0)*clhs84 + clhs240;
lhs(7,9)=DN(1,1)*clhs84 + clhs259;
lhs(7,10)=DN(1,2)*clhs84 + clhs274;
lhs(7,11)=clhs282;
lhs(7,12)=DN(1,0)*clhs113 + clhs248;
lhs(7,13)=DN(1,1)*clhs113 + clhs266;
lhs(7,14)=DN(1,2)*clhs113 + clhs280;
lhs(7,15)=clhs283;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs285 + clhs79;
lhs(8,1)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + DN(2,2)*clhs32 + clhs158;
lhs(8,2)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + DN(2,2)*clhs39 + clhs201;
lhs(8,3)=-DN(0,0)*clhs286 + DN(0,0)*clhs287 - clhs100 + clhs41*clhs99;
lhs(8,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + DN(2,2)*clhs47 + clhs235 + clhs288;
lhs(8,5)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + DN(2,2)*clhs62 + clhs253;
lhs(8,6)=DN(2,0)*clhs64 + DN(2,1)*clhs66 + DN(2,2)*clhs68 + clhs269;
lhs(8,7)=-DN(1,0)*clhs286 + DN(1,0)*clhs287 + clhs239*clhs41 - clhs240;
lhs(8,8)=DN(2,0)*clhs72 + DN(2,1)*clhs74 + DN(2,2)*clhs76 + clhs10*clhs289 + clhs291;
lhs(8,9)=DN(2,0)*clhs86 + DN(2,1)*clhs88 + DN(2,2)*clhs91 + clhs293;
lhs(8,10)=DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs97 + clhs294;
lhs(8,11)=DN(2,0)*clhs295;
lhs(8,12)=DN(2,0)*clhs101 + DN(2,1)*clhs103 + DN(2,2)*clhs105 + clhs298 + clhs299;
lhs(8,13)=DN(2,0)*clhs115 + DN(2,1)*clhs117 + DN(2,2)*clhs120 + clhs300;
lhs(8,14)=DN(2,0)*clhs122 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs301;
lhs(8,15)=-DN(3,0)*clhs286 + DN(3,0)*clhs287 - clhs302 + clhs303*clhs41;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs92;
lhs(9,1)=DN(2,0)*clhs29 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs163 + clhs304;
lhs(9,2)=DN(2,0)*clhs37 + DN(2,1)*clhs136 + DN(2,2)*clhs138 + clhs203;
lhs(9,3)=-DN(0,1)*clhs286 + DN(0,1)*clhs287 + clhs169*clhs41 - clhs170;
lhs(9,4)=DN(2,0)*clhs45 + DN(2,1)*clhs141 + DN(2,2)*clhs142 + clhs237;
lhs(9,5)=DN(2,0)*clhs59 + DN(2,1)*clhs144 + DN(2,2)*clhs146 + clhs255 + clhs305;
lhs(9,6)=DN(2,0)*clhs66 + DN(2,1)*clhs150 + DN(2,2)*clhs152 + clhs270;
lhs(9,7)=-DN(1,1)*clhs286 + DN(1,1)*clhs287 + clhs258*clhs41 - clhs259;
lhs(9,8)=DN(2,0)*clhs74 + DN(2,1)*clhs156 + DN(2,2)*clhs157 + clhs293;
lhs(9,9)=DN(2,0)*clhs88 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs10*clhs306 + clhs291;
lhs(9,10)=DN(2,0)*clhs95 + DN(2,1)*clhs165 + DN(2,2)*clhs167 + clhs308;
lhs(9,11)=DN(2,1)*clhs295;
lhs(9,12)=DN(2,0)*clhs103 + DN(2,1)*clhs171 + DN(2,2)*clhs172 + clhs309;
lhs(9,13)=DN(2,0)*clhs117 + DN(2,1)*clhs174 + DN(2,2)*clhs176 + clhs311 + clhs312;
lhs(9,14)=DN(2,0)*clhs124 + DN(2,1)*clhs180 + DN(2,2)*clhs182 + clhs313;
lhs(9,15)=-DN(3,1)*clhs286 + DN(3,1)*clhs287 - clhs314 + clhs315*clhs41;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs131 + DN(2,2)*clhs186 + clhs98;
lhs(10,1)=DN(2,0)*clhs32 + DN(2,1)*clhs134 + DN(2,2)*clhs187 + clhs168;
lhs(10,2)=DN(2,0)*clhs39 + DN(2,1)*clhs138 + DN(2,2)*clhs188 + clhs206 + clhs304;
lhs(10,3)=-DN(0,2)*clhs286 + DN(0,2)*clhs287 + clhs207*clhs41 - clhs208;
lhs(10,4)=DN(2,0)*clhs47 + DN(2,1)*clhs142 + DN(2,2)*clhs190 + clhs238;
lhs(10,5)=DN(2,0)*clhs62 + DN(2,1)*clhs146 + DN(2,2)*clhs193 + clhs257;
lhs(10,6)=DN(2,0)*clhs68 + DN(2,1)*clhs152 + DN(2,2)*clhs195 + clhs272 + clhs305;
lhs(10,7)=-DN(1,2)*clhs286 + DN(1,2)*clhs287 + clhs273*clhs41 - clhs274;
lhs(10,8)=DN(2,0)*clhs76 + DN(2,1)*clhs157 + DN(2,2)*clhs200 + clhs294;
lhs(10,9)=DN(2,0)*clhs91 + DN(2,1)*clhs161 + DN(2,2)*clhs202 + clhs308;
lhs(10,10)=DN(2,0)*clhs97 + DN(2,1)*clhs167 + DN(2,2)*clhs204 + clhs10*clhs316 + clhs291;
lhs(10,11)=DN(2,2)*clhs295;
lhs(10,12)=DN(2,0)*clhs105 + DN(2,1)*clhs172 + DN(2,2)*clhs209 + clhs318;
lhs(10,13)=DN(2,0)*clhs120 + DN(2,1)*clhs176 + DN(2,2)*clhs211 + clhs319;
lhs(10,14)=DN(2,0)*clhs126 + DN(2,1)*clhs182 + DN(2,2)*clhs213 + clhs312 + clhs321;
lhs(10,15)=-DN(3,2)*clhs286 + DN(3,2)*clhs287 - clhs322 + clhs323*clhs41;
lhs(11,0)=DN(2,0)*clhs23 + clhs99;
lhs(11,1)=DN(2,1)*clhs23 + clhs169;
lhs(11,2)=DN(2,2)*clhs23 + clhs207;
lhs(11,3)=clhs220;
lhs(11,4)=DN(2,0)*clhs55 + clhs239;
lhs(11,5)=DN(2,1)*clhs55 + clhs258;
lhs(11,6)=DN(2,2)*clhs55 + clhs273;
lhs(11,7)=clhs282;
lhs(11,8)=DN(2,0)*clhs324;
lhs(11,9)=DN(2,1)*clhs324;
lhs(11,10)=DN(2,2)*clhs324;
lhs(11,11)=clhs22*(clhs289 + clhs306 + clhs316);
lhs(11,12)=DN(2,0)*clhs113 + clhs303;
lhs(11,13)=DN(2,1)*clhs113 + clhs315;
lhs(11,14)=DN(2,2)*clhs113 + clhs323;
lhs(11,15)=clhs325;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs108 + clhs327;
lhs(12,1)=DN(3,0)*clhs27 + DN(3,1)*clhs29 + DN(3,2)*clhs32 + clhs173;
lhs(12,2)=DN(3,0)*clhs35 + DN(3,1)*clhs37 + DN(3,2)*clhs39 + clhs210;
lhs(12,3)=-DN(0,0)*clhs328 + DN(0,0)*clhs329 + clhs128*clhs41 - clhs129;
lhs(12,4)=DN(3,0)*clhs43 + DN(3,1)*clhs45 + DN(3,2)*clhs47 + clhs243 + clhs330;
lhs(12,5)=DN(3,0)*clhs57 + DN(3,1)*clhs59 + DN(3,2)*clhs62 + clhs260;
lhs(12,6)=DN(3,0)*clhs64 + DN(3,1)*clhs66 + DN(3,2)*clhs68 + clhs275;
lhs(12,7)=-DN(1,0)*clhs328 + DN(1,0)*clhs329 + clhs247*clhs41 - clhs248;
lhs(12,8)=DN(3,0)*clhs72 + DN(3,1)*clhs74 + DN(3,2)*clhs76 + clhs298 + clhs331;
lhs(12,9)=DN(3,0)*clhs86 + DN(3,1)*clhs88 + DN(3,2)*clhs91 + clhs309;
lhs(12,10)=DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs97 + clhs318;
lhs(12,11)=-DN(2,0)*clhs328 + DN(2,0)*clhs329 + clhs302*clhs41 - clhs303;
lhs(12,12)=DN(3,0)*clhs101 + DN(3,1)*clhs103 + DN(3,2)*clhs105 + clhs10*clhs332 + clhs334;
lhs(12,13)=DN(3,0)*clhs115 + DN(3,1)*clhs117 + DN(3,2)*clhs120 + clhs336;
lhs(12,14)=DN(3,0)*clhs122 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs337;
lhs(12,15)=DN(3,0)*clhs338;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs121;
lhs(13,1)=DN(3,0)*clhs29 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs178 + clhs339;
lhs(13,2)=DN(3,0)*clhs37 + DN(3,1)*clhs136 + DN(3,2)*clhs138 + clhs212;
lhs(13,3)=-DN(0,1)*clhs328 + DN(0,1)*clhs329 + clhs184*clhs41 - clhs185;
lhs(13,4)=DN(3,0)*clhs45 + DN(3,1)*clhs141 + DN(3,2)*clhs142 + clhs245;
lhs(13,5)=DN(3,0)*clhs59 + DN(3,1)*clhs144 + DN(3,2)*clhs146 + clhs262 + clhs340;
lhs(13,6)=DN(3,0)*clhs66 + DN(3,1)*clhs150 + DN(3,2)*clhs152 + clhs276;
lhs(13,7)=-DN(1,1)*clhs328 + DN(1,1)*clhs329 + clhs265*clhs41 - clhs266;
lhs(13,8)=DN(3,0)*clhs74 + DN(3,1)*clhs156 + DN(3,2)*clhs157 + clhs300;
lhs(13,9)=DN(3,0)*clhs88 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs311 + clhs341;
lhs(13,10)=DN(3,0)*clhs95 + DN(3,1)*clhs165 + DN(3,2)*clhs167 + clhs319;
lhs(13,11)=-DN(2,1)*clhs328 + DN(2,1)*clhs329 + clhs314*clhs41 - clhs315;
lhs(13,12)=DN(3,0)*clhs103 + DN(3,1)*clhs171 + DN(3,2)*clhs172 + clhs336;
lhs(13,13)=DN(3,0)*clhs117 + DN(3,1)*clhs174 + DN(3,2)*clhs176 + clhs10*clhs342 + clhs334;
lhs(13,14)=DN(3,0)*clhs124 + DN(3,1)*clhs180 + DN(3,2)*clhs182 + clhs343;
lhs(13,15)=DN(3,1)*clhs338;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs131 + DN(3,2)*clhs186 + clhs127;
lhs(14,1)=DN(3,0)*clhs32 + DN(3,1)*clhs134 + DN(3,2)*clhs187 + clhs183;
lhs(14,2)=DN(3,0)*clhs39 + DN(3,1)*clhs138 + DN(3,2)*clhs188 + clhs215 + clhs339;
lhs(14,3)=-DN(0,2)*clhs328 + DN(0,2)*clhs329 + clhs216*clhs41 - clhs217;
lhs(14,4)=DN(3,0)*clhs47 + DN(3,1)*clhs142 + DN(3,2)*clhs190 + clhs246;
lhs(14,5)=DN(3,0)*clhs62 + DN(3,1)*clhs146 + DN(3,2)*clhs193 + clhs264;
lhs(14,6)=DN(3,0)*clhs68 + DN(3,1)*clhs152 + DN(3,2)*clhs195 + clhs278 + clhs340;
lhs(14,7)=-DN(1,2)*clhs328 + DN(1,2)*clhs329 + clhs279*clhs41 - clhs280;
lhs(14,8)=DN(3,0)*clhs76 + DN(3,1)*clhs157 + DN(3,2)*clhs200 + clhs301;
lhs(14,9)=DN(3,0)*clhs91 + DN(3,1)*clhs161 + DN(3,2)*clhs202 + clhs313;
lhs(14,10)=DN(3,0)*clhs97 + DN(3,1)*clhs167 + DN(3,2)*clhs204 + clhs321 + clhs341;
lhs(14,11)=-DN(2,2)*clhs328 + DN(2,2)*clhs329 + clhs322*clhs41 - clhs323;
lhs(14,12)=DN(3,0)*clhs105 + DN(3,1)*clhs172 + DN(3,2)*clhs209 + clhs337;
lhs(14,13)=DN(3,0)*clhs120 + DN(3,1)*clhs176 + DN(3,2)*clhs211 + clhs343;
lhs(14,14)=DN(3,0)*clhs126 + DN(3,1)*clhs182 + DN(3,2)*clhs213 + clhs10*clhs344 + clhs334;
lhs(14,15)=DN(3,2)*clhs338;
lhs(15,0)=DN(3,0)*clhs23 + clhs128;
lhs(15,1)=DN(3,1)*clhs23 + clhs184;
lhs(15,2)=DN(3,2)*clhs23 + clhs216;
lhs(15,3)=clhs221;
lhs(15,4)=DN(3,0)*clhs55 + clhs247;
lhs(15,5)=DN(3,1)*clhs55 + clhs265;
lhs(15,6)=DN(3,2)*clhs55 + clhs279;
lhs(15,7)=clhs283;
lhs(15,8)=DN(3,0)*clhs84 + clhs302;
lhs(15,9)=DN(3,1)*clhs84 + clhs314;
lhs(15,10)=DN(3,2)*clhs84 + clhs322;
lhs(15,11)=clhs325;
lhs(15,12)=DN(3,0)*clhs345;
lhs(15,13)=DN(3,1)*clhs345;
lhs(15,14)=DN(3,2)*clhs345;
lhs(15,15)=clhs22*(clhs332 + clhs342 + clhs344);


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
    const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vn;
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
const double crhs3 = N[0]*an(0,0) + N[1]*an(1,0) + N[2]*an(2,0);
const double crhs4 = N[0]*rho;
const double crhs5 = N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0));
const double crhs6 = bdf0*crhs4;
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs10 = rho*(crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs11 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs12 = crhs11 + crhs7 - volume_error_ratio;
const double crhs13 = rho*stab_c2*sqrt(pow(crhs8, 2) + pow(crhs9, 2));
const double crhs14 = crhs12*(crhs13*h/stab_c1 + mu);
const double crhs15 = bdf0*rho;
const double crhs16 = crhs15*crhs5;
const double crhs17 = 1.0/(K_darcy + crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs18 = crhs17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs10 + crhs16 + crhs2 + crhs3);
const double crhs19 = K_darcy*N[0];
const double crhs20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs21 = crhs20*crhs4;
const double crhs22 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs23 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs24 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs25 = N[0]*an(0,1) + N[1]*an(1,1) + N[2]*an(2,1);
const double crhs26 = N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1));
const double crhs27 = rho*(crhs11*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs28 = crhs15*crhs26;
const double crhs29 = crhs17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs23 + crhs24 + crhs25 + crhs27 + crhs28);
const double crhs30 = N[1]*rho;
const double crhs31 = K_darcy*N[1];
const double crhs32 = crhs20*crhs30;
const double crhs33 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs34 = N[2]*rho;
const double crhs35 = K_darcy*N[2];
const double crhs36 = crhs20*crhs34;
const double crhs37 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 + crhs18*crhs19 - crhs18*crhs21 - crhs18*crhs22 - crhs3*crhs4 - crhs5*crhs6;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs14 - DN(0,1)*stress[1] + N[0]*crhs23 - N[0]*crhs24 - N[0]*crhs27 + crhs19*crhs29 - crhs21*crhs29 - crhs22*crhs29 - crhs25*crhs4 - crhs26*crhs6;
rhs[2]=-DN(0,0)*crhs18 - DN(0,1)*crhs29 - N[0]*crhs12;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs16 - N[1]*crhs2 + crhs18*crhs31 - crhs18*crhs32 - crhs18*crhs33 - crhs3*crhs30;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs14 - DN(1,1)*stress[1] + N[1]*crhs23 - N[1]*crhs24 - N[1]*crhs27 - N[1]*crhs28 - crhs25*crhs30 + crhs29*crhs31 - crhs29*crhs32 - crhs29*crhs33;
rhs[5]=-DN(1,0)*crhs18 - DN(1,1)*crhs29 - N[1]*crhs12;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs16 - N[2]*crhs2 + crhs18*crhs35 - crhs18*crhs36 - crhs18*crhs37 - crhs3*crhs34;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs14 - DN(2,1)*stress[1] + N[2]*crhs23 - N[2]*crhs24 - N[2]*crhs27 - N[2]*crhs28 - crhs25*crhs34 + crhs29*crhs35 - crhs29*crhs36 - crhs29*crhs37;
rhs[8]=-DN(2,0)*crhs18 - DN(2,1)*crhs29 - N[2]*crhs12;


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
    const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vn;
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
const double crhs3 = N[0]*an(0,0) + N[1]*an(1,0) + N[2]*an(2,0) + N[3]*an(3,0);
const double crhs4 = N[0]*rho;
const double crhs5 = N[0]*(v(0,0) - vn(0,0)) + N[1]*(v(1,0) - vn(1,0)) + N[2]*(v(2,0) - vn(2,0)) + N[3]*(v(3,0) - vn(3,0));
const double crhs6 = bdf0*crhs4;
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs11 = rho*(crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
const double crhs12 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs13 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs14 = crhs12 + crhs13 + crhs7 - volume_error_ratio;
const double crhs15 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs16 = crhs14*(crhs15*h/stab_c1 + mu);
const double crhs17 = bdf0*rho;
const double crhs18 = crhs17*crhs5;
const double crhs19 = 1.0/(K_darcy + crhs15/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs20 = crhs19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs11 + crhs18 + crhs2 + crhs3);
const double crhs21 = K_darcy*N[0];
const double crhs22 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs23 = crhs22*crhs4;
const double crhs24 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10);
const double crhs25 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs26 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs27 = N[0]*an(0,1) + N[1]*an(1,1) + N[2]*an(2,1) + N[3]*an(3,1);
const double crhs28 = N[0]*(v(0,1) - vn(0,1)) + N[1]*(v(1,1) - vn(1,1)) + N[2]*(v(2,1) - vn(2,1)) + N[3]*(v(3,1) - vn(3,1));
const double crhs29 = rho*(crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs12*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)));
const double crhs30 = crhs17*crhs28;
const double crhs31 = crhs19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs25 + crhs26 + crhs27 + crhs29 + crhs30);
const double crhs32 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs33 = K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs34 = N[0]*an(0,2) + N[1]*an(1,2) + N[2]*an(2,2) + N[3]*an(3,2);
const double crhs35 = N[0]*(v(0,2) - vn(0,2)) + N[1]*(v(1,2) - vn(1,2)) + N[2]*(v(2,2) - vn(2,2)) + N[3]*(v(3,2) - vn(3,2));
const double crhs36 = rho*(crhs10*crhs13 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs37 = crhs17*crhs35;
const double crhs38 = crhs19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs32 + crhs33 + crhs34 + crhs36 + crhs37);
const double crhs39 = N[1]*rho;
const double crhs40 = K_darcy*N[1];
const double crhs41 = crhs22*crhs39;
const double crhs42 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10);
const double crhs43 = N[2]*rho;
const double crhs44 = K_darcy*N[2];
const double crhs45 = crhs22*crhs43;
const double crhs46 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10);
const double crhs47 = N[3]*rho;
const double crhs48 = K_darcy*N[3];
const double crhs49 = crhs22*crhs47;
const double crhs50 = rho*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs16 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs11 - N[0]*crhs2 + crhs20*crhs21 - crhs20*crhs23 - crhs20*crhs24 - crhs3*crhs4 - crhs5*crhs6;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs16 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs29 + crhs21*crhs31 - crhs23*crhs31 - crhs24*crhs31 - crhs27*crhs4 - crhs28*crhs6;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs16 - DN(0,2)*stress[2] + N[0]*crhs32 - N[0]*crhs33 - N[0]*crhs36 + crhs21*crhs38 - crhs23*crhs38 - crhs24*crhs38 - crhs34*crhs4 - crhs35*crhs6;
rhs[3]=-DN(0,0)*crhs20 - DN(0,1)*crhs31 - DN(0,2)*crhs38 - N[0]*crhs14;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs16 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs11 - N[1]*crhs18 - N[1]*crhs2 + crhs20*crhs40 - crhs20*crhs41 - crhs20*crhs42 - crhs3*crhs39;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs16 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs29 - N[1]*crhs30 - crhs27*crhs39 + crhs31*crhs40 - crhs31*crhs41 - crhs31*crhs42;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs16 - DN(1,2)*stress[2] + N[1]*crhs32 - N[1]*crhs33 - N[1]*crhs36 - N[1]*crhs37 - crhs34*crhs39 + crhs38*crhs40 - crhs38*crhs41 - crhs38*crhs42;
rhs[7]=-DN(1,0)*crhs20 - DN(1,1)*crhs31 - DN(1,2)*crhs38 - N[1]*crhs14;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs16 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs11 - N[2]*crhs18 - N[2]*crhs2 + crhs20*crhs44 - crhs20*crhs45 - crhs20*crhs46 - crhs3*crhs43;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs16 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs29 - N[2]*crhs30 - crhs27*crhs43 + crhs31*crhs44 - crhs31*crhs45 - crhs31*crhs46;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs16 - DN(2,2)*stress[2] + N[2]*crhs32 - N[2]*crhs33 - N[2]*crhs36 - N[2]*crhs37 - crhs34*crhs43 + crhs38*crhs44 - crhs38*crhs45 - crhs38*crhs46;
rhs[11]=-DN(2,0)*crhs20 - DN(2,1)*crhs31 - DN(2,2)*crhs38 - N[2]*crhs14;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs16 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs11 - N[3]*crhs18 - N[3]*crhs2 + crhs20*crhs48 - crhs20*crhs49 - crhs20*crhs50 - crhs3*crhs47;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs16 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs29 - N[3]*crhs30 - crhs27*crhs47 + crhs31*crhs48 - crhs31*crhs49 - crhs31*crhs50;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs16 - DN(3,2)*stress[2] + N[3]*crhs32 - N[3]*crhs33 - N[3]*crhs36 - N[3]*crhs37 - crhs34*crhs47 + crhs38*crhs48 - crhs38*crhs49 - crhs38*crhs50;
rhs[15]=-DN(3,0)*crhs20 - DN(3,1)*crhs31 - DN(3,2)*crhs38 - N[3]*crhs14;


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
    const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vn;
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


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH3 = cH2*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0));
const double cH4 = cH2*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0));
const double cH5 = cH2*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0));
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
const double crhs_ee3 = N[0]*bdf0;
const double crhs_ee4 = N[1]*bdf0;
const double crhs_ee5 = N[2]*bdf0;
const double crhs_ee6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee8 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee6, 2) + pow(crhs_ee7, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee9 = crhs_ee8*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(crhs_ee0*crhs_ee6 + crhs_ee3*(v(0,0) - vn(0,0)) + crhs_ee4*(v(1,0) - vn(1,0)) + crhs_ee5*(v(2,0) - vn(2,0)) + crhs_ee7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee10 = crhs_ee8*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(crhs_ee1*crhs_ee7 + crhs_ee3*(v(0,1) - vn(0,1)) + crhs_ee4*(v(1,1) - vn(1,1)) + crhs_ee5*(v(2,1) - vn(2,1)) + crhs_ee6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee9 - DNenr(0,1)*crhs_ee10 - Nenr[0]*crhs_ee2;
rhs_ee[1]=-DNenr(1,0)*crhs_ee9 - DNenr(1,1)*crhs_ee10 - Nenr[1]*crhs_ee2;
rhs_ee[2]=-DNenr(2,0)*crhs_ee9 - DNenr(2,1)*crhs_ee10 - Nenr[2]*crhs_ee2;


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
    const auto an = rData.Acceleration;
    const auto &vnn = rData.Velocity_OldStep2;
    const auto &vmesh = rData.MeshVelocity;
    const auto &vconv = v - vn;
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


    const double cH0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 = cH3*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0));
const double cH5 = cH3*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0));
const double cH6 = cH3*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0));
const double cH7 = cH3*(K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0));
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
const double crhs_ee4 = N[0]*bdf0;
const double crhs_ee5 = N[1]*bdf0;
const double crhs_ee6 = N[2]*bdf0;
const double crhs_ee7 = N[3]*bdf0;
const double crhs_ee8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee11 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee8, 2) + pow(crhs_ee9, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee12 = crhs_ee11*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(crhs_ee0*crhs_ee8 + crhs_ee10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs_ee4*(v(0,0) - vn(0,0)) + crhs_ee5*(v(1,0) - vn(1,0)) + crhs_ee6*(v(2,0) - vn(2,0)) + crhs_ee7*(v(3,0) - vn(3,0)) + crhs_ee9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0))));
const double crhs_ee13 = crhs_ee11*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(crhs_ee1*crhs_ee9 + crhs_ee10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs_ee4*(v(0,1) - vn(0,1)) + crhs_ee5*(v(1,1) - vn(1,1)) + crhs_ee6*(v(2,1) - vn(2,1)) + crhs_ee7*(v(3,1) - vn(3,1)) + crhs_ee8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1))));
const double crhs_ee14 = crhs_ee11*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(crhs_ee10*crhs_ee2 + crhs_ee4*(v(0,2) - vn(0,2)) + crhs_ee5*(v(1,2) - vn(1,2)) + crhs_ee6*(v(2,2) - vn(2,2)) + crhs_ee7*(v(3,2) - vn(3,2)) + crhs_ee8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
rhs_ee[0]=-DNenr(0,0)*crhs_ee12 - DNenr(0,1)*crhs_ee13 - DNenr(0,2)*crhs_ee14 - Nenr[0]*crhs_ee3;
rhs_ee[1]=-DNenr(1,0)*crhs_ee12 - DNenr(1,1)*crhs_ee13 - DNenr(1,2)*crhs_ee14 - Nenr[1]*crhs_ee3;
rhs_ee[2]=-DNenr(2,0)*crhs_ee12 - DNenr(2,1)*crhs_ee13 - DNenr(2,2)*crhs_ee14 - Nenr[2]*crhs_ee3;
rhs_ee[3]=-DNenr(3,0)*crhs_ee12 - DNenr(3,1)*crhs_ee13 - DNenr(3,2)*crhs_ee14 - Nenr[3]*crhs_ee3;


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
