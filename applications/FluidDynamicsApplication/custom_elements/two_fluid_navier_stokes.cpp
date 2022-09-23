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
const double clhs1 = DN(0,0)*clhs0;
const double clhs2 = C(0,2)*DN(0,0);
const double clhs3 = C(2,2)*DN(0,1) + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs8 = rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2));
const double clhs9 = clhs8*h/stab_c1 + mu;
const double clhs10 = bdf0*rho;
const double clhs11 = 1.0/(K_darcy + clhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs12 = 1.0*clhs11;
const double clhs13 = clhs10*clhs12;
const double clhs14 = pow(N[0], 2);
const double clhs15 = rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7);
const double clhs16 = K_darcy*N[0];
const double clhs17 = N[0]*clhs10;
const double clhs18 = clhs15 + clhs16 + clhs17;
const double clhs19 = 1.0*clhs15;
const double clhs20 = clhs11*clhs19;
const double clhs21 = 1.0*clhs16;
const double clhs22 = clhs11*clhs21;
const double clhs23 = clhs12*clhs18;
const double clhs24 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double clhs25 = N[0]*clhs24;
const double clhs26 = K_darcy*clhs14 + N[0]*clhs15 + clhs10*clhs14 + clhs18*clhs20 - clhs18*clhs22 + clhs23*clhs25;
const double clhs27 = C(0,1)*DN(0,1) + clhs2;
const double clhs28 = DN(0,0)*clhs27;
const double clhs29 = C(1,2)*DN(0,1);
const double clhs30 = C(2,2)*DN(0,0) + clhs29;
const double clhs31 = DN(0,1)*clhs30;
const double clhs32 = DN(0,0)*clhs9;
const double clhs33 = DN(0,1)*clhs32;
const double clhs34 = DN(0,0)*clhs3;
const double clhs35 = C(0,1)*DN(0,0) + clhs29;
const double clhs36 = DN(0,1)*clhs35;
const double clhs37 = clhs12*clhs24;
const double clhs38 = N[0]*clhs37 - N[0] + clhs20 - clhs22;
const double clhs39 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs40 = DN(0,0)*clhs39;
const double clhs41 = C(0,2)*DN(1,0);
const double clhs42 = C(2,2)*DN(1,1) + clhs41;
const double clhs43 = DN(0,1)*clhs42;
const double clhs44 = DN(1,0)*clhs0;
const double clhs45 = DN(1,1)*clhs3;
const double clhs46 = DN(0,0)*DN(1,0);
const double clhs47 = N[1]*clhs16 + N[1]*clhs17;
const double clhs48 = clhs46*clhs9 + clhs47;
const double clhs49 = rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7);
const double clhs50 = K_darcy*N[1];
const double clhs51 = N[1]*clhs10;
const double clhs52 = clhs49 + clhs50 + clhs51;
const double clhs53 = clhs12*clhs52;
const double clhs54 = N[0]*clhs49 + clhs20*clhs52 - clhs22*clhs52 + clhs25*clhs53;
const double clhs55 = C(0,1)*DN(1,1) + clhs41;
const double clhs56 = DN(0,0)*clhs55;
const double clhs57 = C(1,2)*DN(1,1);
const double clhs58 = C(2,2)*DN(1,0) + clhs57;
const double clhs59 = DN(0,1)*clhs58;
const double clhs60 = DN(1,1)*clhs32;
const double clhs61 = DN(1,0)*clhs3;
const double clhs62 = DN(1,1)*clhs35;
const double clhs63 = DN(0,0)*N[1];
const double clhs64 = DN(1,0)*N[0];
const double clhs65 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs66 = DN(0,0)*clhs65;
const double clhs67 = C(0,2)*DN(2,0);
const double clhs68 = C(2,2)*DN(2,1) + clhs67;
const double clhs69 = DN(0,1)*clhs68;
const double clhs70 = DN(2,0)*clhs0;
const double clhs71 = DN(2,1)*clhs3;
const double clhs72 = DN(0,0)*DN(2,0);
const double clhs73 = N[2]*clhs16 + N[2]*clhs17;
const double clhs74 = clhs72*clhs9 + clhs73;
const double clhs75 = rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7);
const double clhs76 = K_darcy*N[2];
const double clhs77 = N[2]*clhs10;
const double clhs78 = clhs75 + clhs76 + clhs77;
const double clhs79 = clhs12*clhs78;
const double clhs80 = N[0]*clhs75 + clhs20*clhs78 - clhs22*clhs78 + clhs25*clhs79;
const double clhs81 = C(0,1)*DN(2,1) + clhs67;
const double clhs82 = DN(0,0)*clhs81;
const double clhs83 = C(1,2)*DN(2,1);
const double clhs84 = C(2,2)*DN(2,0) + clhs83;
const double clhs85 = DN(0,1)*clhs84;
const double clhs86 = DN(2,1)*clhs32;
const double clhs87 = DN(2,0)*clhs3;
const double clhs88 = DN(2,1)*clhs35;
const double clhs89 = DN(0,0)*N[2];
const double clhs90 = DN(2,0)*N[0];
const double clhs91 = DN(0,0)*clhs30;
const double clhs92 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs93 = DN(0,1)*clhs92;
const double clhs94 = pow(DN(0,1), 2);
const double clhs95 = DN(0,0)*clhs42;
const double clhs96 = C(0,1)*DN(1,0) + clhs57;
const double clhs97 = DN(0,1)*clhs96;
const double clhs98 = DN(0,1)*clhs9;
const double clhs99 = DN(1,0)*clhs98;
const double clhs100 = DN(1,0)*clhs27;
const double clhs101 = DN(1,1)*clhs30;
const double clhs102 = DN(0,0)*clhs58;
const double clhs103 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs104 = DN(0,1)*clhs103;
const double clhs105 = DN(1,0)*clhs30;
const double clhs106 = DN(1,1)*clhs92;
const double clhs107 = DN(0,1)*DN(1,1);
const double clhs108 = clhs107*clhs9 + clhs47;
const double clhs109 = DN(0,1)*N[1];
const double clhs110 = DN(1,1)*N[0];
const double clhs111 = DN(0,0)*clhs68;
const double clhs112 = C(0,1)*DN(2,0) + clhs83;
const double clhs113 = DN(0,1)*clhs112;
const double clhs114 = DN(2,0)*clhs98;
const double clhs115 = DN(2,0)*clhs27;
const double clhs116 = DN(2,1)*clhs30;
const double clhs117 = DN(0,0)*clhs84;
const double clhs118 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs119 = DN(0,1)*clhs118;
const double clhs120 = DN(2,0)*clhs30;
const double clhs121 = DN(2,1)*clhs92;
const double clhs122 = DN(0,1)*DN(2,1);
const double clhs123 = clhs122*clhs9 + clhs73;
const double clhs124 = DN(0,1)*N[2];
const double clhs125 = DN(2,1)*N[0];
const double clhs126 = N[0] + clhs11*(1.0*clhs17 + clhs19 + clhs21);
const double clhs127 = clhs12*(clhs107 + clhs46);
const double clhs128 = clhs12*(clhs122 + clhs72);
const double clhs129 = clhs12*clhs49;
const double clhs130 = clhs12*clhs50;
const double clhs131 = N[1]*clhs24;
const double clhs132 = N[1]*clhs15 + clhs129*clhs18 - clhs130*clhs18 + clhs131*clhs23;
const double clhs133 = DN(1,0)*clhs39;
const double clhs134 = DN(1,1)*clhs42;
const double clhs135 = pow(DN(1,0), 2);
const double clhs136 = pow(N[1], 2);
const double clhs137 = K_darcy*clhs136 + N[1]*clhs49 + clhs10*clhs136 + clhs131*clhs53 + clhs49*clhs53 - clhs50*clhs53;
const double clhs138 = DN(1,0)*clhs55;
const double clhs139 = DN(1,1)*clhs58;
const double clhs140 = DN(1,0)*clhs9;
const double clhs141 = DN(1,1)*clhs140;
const double clhs142 = DN(1,0)*clhs42;
const double clhs143 = DN(1,1)*clhs96;
const double clhs144 = N[1]*clhs37 - N[1] + clhs129 - clhs130;
const double clhs145 = DN(1,0)*clhs65;
const double clhs146 = DN(1,1)*clhs68;
const double clhs147 = DN(2,0)*clhs39;
const double clhs148 = DN(2,1)*clhs42;
const double clhs149 = DN(1,0)*DN(2,0);
const double clhs150 = N[2]*clhs50 + N[2]*clhs51;
const double clhs151 = clhs149*clhs9 + clhs150;
const double clhs152 = N[1]*clhs75 + clhs131*clhs79 + clhs49*clhs79 - clhs50*clhs79;
const double clhs153 = DN(1,0)*clhs81;
const double clhs154 = DN(1,1)*clhs84;
const double clhs155 = DN(2,1)*clhs140;
const double clhs156 = DN(2,0)*clhs42;
const double clhs157 = DN(2,1)*clhs96;
const double clhs158 = DN(1,0)*N[2];
const double clhs159 = DN(2,0)*N[1];
const double clhs160 = DN(1,0)*clhs58;
const double clhs161 = DN(1,1)*clhs103;
const double clhs162 = pow(DN(1,1), 2);
const double clhs163 = DN(1,0)*clhs68;
const double clhs164 = DN(1,1)*clhs112;
const double clhs165 = DN(2,0)*clhs9;
const double clhs166 = DN(1,1)*clhs165;
const double clhs167 = DN(2,0)*clhs55;
const double clhs168 = DN(2,1)*clhs58;
const double clhs169 = DN(1,0)*clhs84;
const double clhs170 = DN(1,1)*clhs118;
const double clhs171 = DN(2,0)*clhs58;
const double clhs172 = DN(2,1)*clhs103;
const double clhs173 = DN(1,1)*DN(2,1);
const double clhs174 = clhs150 + clhs173*clhs9;
const double clhs175 = DN(1,1)*N[2];
const double clhs176 = DN(2,1)*N[1];
const double clhs177 = N[1] + clhs11*(1.0*clhs49 + 1.0*clhs50 + 1.0*clhs51);
const double clhs178 = clhs12*(clhs149 + clhs173);
const double clhs179 = N[2]*clhs24;
const double clhs180 = N[2]*clhs15 + clhs179*clhs23 + clhs23*clhs75 - clhs23*clhs76;
const double clhs181 = clhs12*clhs76;
const double clhs182 = clhs12*clhs75;
const double clhs183 = N[2]*clhs49 + clhs179*clhs53 + clhs53*clhs75 - clhs53*clhs76;
const double clhs184 = DN(2,0)*clhs65;
const double clhs185 = DN(2,1)*clhs68;
const double clhs186 = pow(DN(2,0), 2);
const double clhs187 = pow(N[2], 2);
const double clhs188 = K_darcy*clhs187 + N[2]*clhs75 + clhs10*clhs187 + clhs179*clhs79 + clhs75*clhs79 - clhs76*clhs79;
const double clhs189 = DN(2,0)*clhs81;
const double clhs190 = DN(2,1)*clhs84;
const double clhs191 = DN(2,1)*clhs165;
const double clhs192 = DN(2,0)*clhs68;
const double clhs193 = DN(2,1)*clhs112;
const double clhs194 = N[2]*clhs37 - N[2] - clhs181 + clhs182;
const double clhs195 = DN(2,0)*clhs84;
const double clhs196 = DN(2,1)*clhs118;
const double clhs197 = pow(DN(2,1), 2);
const double clhs198 = N[2] + clhs11*(1.0*clhs75 + 1.0*clhs76 + 1.0*clhs77);
lhs(0,0)=-clhs1*clhs13 + clhs1 - clhs13*clhs4 + clhs26 + clhs4 + clhs5*clhs9;
lhs(0,1)=-clhs13*clhs34 - clhs13*clhs36 + clhs28 + clhs31 + clhs33;
lhs(0,2)=DN(0,0)*clhs38;
lhs(0,3)=-clhs13*clhs44 - clhs13*clhs45 + clhs40 + clhs43 + clhs48 + clhs54;
lhs(0,4)=-clhs13*clhs61 - clhs13*clhs62 + clhs56 + clhs59 + clhs60;
lhs(0,5)=DN(1,0)*clhs20 - DN(1,0)*clhs22 + clhs37*clhs64 - clhs63;
lhs(0,6)=-clhs13*clhs70 - clhs13*clhs71 + clhs66 + clhs69 + clhs74 + clhs80;
lhs(0,7)=-clhs13*clhs87 - clhs13*clhs88 + clhs82 + clhs85 + clhs86;
lhs(0,8)=DN(2,0)*clhs20 - DN(2,0)*clhs22 + clhs37*clhs90 - clhs89;
lhs(1,0)=-clhs13*clhs28 - clhs13*clhs31 + clhs33 + clhs34 + clhs36;
lhs(1,1)=-clhs13*clhs91 - clhs13*clhs93 + clhs26 + clhs9*clhs94 + clhs91 + clhs93;
lhs(1,2)=DN(0,1)*clhs38;
lhs(1,3)=-clhs100*clhs13 - clhs101*clhs13 + clhs95 + clhs97 + clhs99;
lhs(1,4)=clhs102 + clhs104 - clhs105*clhs13 - clhs106*clhs13 + clhs108 + clhs54;
lhs(1,5)=DN(1,1)*clhs20 - DN(1,1)*clhs22 - clhs109 + clhs110*clhs37;
lhs(1,6)=clhs111 + clhs113 + clhs114 - clhs115*clhs13 - clhs116*clhs13;
lhs(1,7)=clhs117 + clhs119 - clhs120*clhs13 - clhs121*clhs13 + clhs123 + clhs80;
lhs(1,8)=DN(2,1)*clhs20 - DN(2,1)*clhs22 - clhs124 + clhs125*clhs37;
lhs(2,0)=DN(0,0)*clhs126;
lhs(2,1)=DN(0,1)*clhs126;
lhs(2,2)=clhs12*(clhs5 + clhs94);
lhs(2,3)=DN(0,0)*clhs53 + clhs64;
lhs(2,4)=DN(0,1)*clhs53 + clhs110;
lhs(2,5)=clhs127;
lhs(2,6)=DN(0,0)*clhs79 + clhs90;
lhs(2,7)=DN(0,1)*clhs79 + clhs125;
lhs(2,8)=clhs128;
lhs(3,0)=-clhs13*clhs40 - clhs13*clhs43 + clhs132 + clhs44 + clhs45 + clhs48;
lhs(3,1)=clhs100 + clhs101 - clhs13*clhs95 - clhs13*clhs97 + clhs99;
lhs(3,2)=DN(0,0)*clhs129 - DN(0,0)*clhs130 + clhs37*clhs63 - clhs64;
lhs(3,3)=-clhs13*clhs133 - clhs13*clhs134 + clhs133 + clhs134 + clhs135*clhs9 + clhs137;
lhs(3,4)=-clhs13*clhs142 - clhs13*clhs143 + clhs138 + clhs139 + clhs141;
lhs(3,5)=DN(1,0)*clhs144;
lhs(3,6)=-clhs13*clhs147 - clhs13*clhs148 + clhs145 + clhs146 + clhs151 + clhs152;
lhs(3,7)=-clhs13*clhs156 - clhs13*clhs157 + clhs153 + clhs154 + clhs155;
lhs(3,8)=DN(2,0)*clhs129 - DN(2,0)*clhs130 - clhs158 + clhs159*clhs37;
lhs(4,0)=-clhs13*clhs56 - clhs13*clhs59 + clhs60 + clhs61 + clhs62;
lhs(4,1)=-clhs102*clhs13 - clhs104*clhs13 + clhs105 + clhs106 + clhs108 + clhs132;
lhs(4,2)=DN(0,1)*clhs129 - DN(0,1)*clhs130 + clhs109*clhs37 - clhs110;
lhs(4,3)=-clhs13*clhs138 - clhs13*clhs139 + clhs141 + clhs142 + clhs143;
lhs(4,4)=-clhs13*clhs160 - clhs13*clhs161 + clhs137 + clhs160 + clhs161 + clhs162*clhs9;
lhs(4,5)=DN(1,1)*clhs144;
lhs(4,6)=-clhs13*clhs167 - clhs13*clhs168 + clhs163 + clhs164 + clhs166;
lhs(4,7)=-clhs13*clhs171 - clhs13*clhs172 + clhs152 + clhs169 + clhs170 + clhs174;
lhs(4,8)=DN(2,1)*clhs129 - DN(2,1)*clhs130 - clhs175 + clhs176*clhs37;
lhs(5,0)=DN(1,0)*clhs23 + clhs63;
lhs(5,1)=DN(1,1)*clhs23 + clhs109;
lhs(5,2)=clhs127;
lhs(5,3)=DN(1,0)*clhs177;
lhs(5,4)=DN(1,1)*clhs177;
lhs(5,5)=clhs12*(clhs135 + clhs162);
lhs(5,6)=DN(1,0)*clhs79 + clhs159;
lhs(5,7)=DN(1,1)*clhs79 + clhs176;
lhs(5,8)=clhs178;
lhs(6,0)=-clhs13*clhs66 - clhs13*clhs69 + clhs180 + clhs70 + clhs71 + clhs74;
lhs(6,1)=-clhs111*clhs13 - clhs113*clhs13 + clhs114 + clhs115 + clhs116;
lhs(6,2)=-DN(0,0)*clhs181 + DN(0,0)*clhs182 + clhs37*clhs89 - clhs90;
lhs(6,3)=-clhs13*clhs145 - clhs13*clhs146 + clhs147 + clhs148 + clhs151 + clhs183;
lhs(6,4)=-clhs13*clhs163 - clhs13*clhs164 + clhs166 + clhs167 + clhs168;
lhs(6,5)=-DN(1,0)*clhs181 + DN(1,0)*clhs182 + clhs158*clhs37 - clhs159;
lhs(6,6)=-clhs13*clhs184 - clhs13*clhs185 + clhs184 + clhs185 + clhs186*clhs9 + clhs188;
lhs(6,7)=-clhs13*clhs192 - clhs13*clhs193 + clhs189 + clhs190 + clhs191;
lhs(6,8)=DN(2,0)*clhs194;
lhs(7,0)=-clhs13*clhs82 - clhs13*clhs85 + clhs86 + clhs87 + clhs88;
lhs(7,1)=-clhs117*clhs13 - clhs119*clhs13 + clhs120 + clhs121 + clhs123 + clhs180;
lhs(7,2)=-DN(0,1)*clhs181 + DN(0,1)*clhs182 + clhs124*clhs37 - clhs125;
lhs(7,3)=-clhs13*clhs153 - clhs13*clhs154 + clhs155 + clhs156 + clhs157;
lhs(7,4)=-clhs13*clhs169 - clhs13*clhs170 + clhs171 + clhs172 + clhs174 + clhs183;
lhs(7,5)=-DN(1,1)*clhs181 + DN(1,1)*clhs182 + clhs175*clhs37 - clhs176;
lhs(7,6)=-clhs13*clhs189 - clhs13*clhs190 + clhs191 + clhs192 + clhs193;
lhs(7,7)=-clhs13*clhs195 - clhs13*clhs196 + clhs188 + clhs195 + clhs196 + clhs197*clhs9;
lhs(7,8)=DN(2,1)*clhs194;
lhs(8,0)=DN(2,0)*clhs23 + clhs89;
lhs(8,1)=DN(2,1)*clhs23 + clhs124;
lhs(8,2)=clhs128;
lhs(8,3)=DN(2,0)*clhs53 + clhs158;
lhs(8,4)=DN(2,1)*clhs53 + clhs175;
lhs(8,5)=clhs178;
lhs(8,6)=DN(2,0)*clhs198;
lhs(8,7)=DN(2,1)*clhs198;
lhs(8,8)=clhs12*(clhs186 + clhs197);


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
const double clhs1 = DN(0,0)*clhs0;
const double clhs2 = C(0,3)*DN(0,0);
const double clhs3 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = C(0,5)*DN(0,0);
const double clhs6 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 = DN(0,2)*clhs6;
const double clhs8 = pow(DN(0,0), 2);
const double clhs9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs11 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs12 = rho*stab_c2*sqrt(pow(clhs10, 2) + pow(clhs11, 2) + pow(clhs9, 2));
const double clhs13 = clhs12*h/stab_c1 + mu;
const double clhs14 = bdf0*rho;
const double clhs15 = 1.0/(K_darcy + clhs12/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs16 = 1.0*clhs15;
const double clhs17 = clhs14*clhs16;
const double clhs18 = pow(N[0], 2);
const double clhs19 = rho*(DN(0,0)*clhs9 + DN(0,1)*clhs10 + DN(0,2)*clhs11);
const double clhs20 = K_darcy*N[0];
const double clhs21 = N[0]*clhs14;
const double clhs22 = clhs19 + clhs20 + clhs21;
const double clhs23 = 1.0*clhs19;
const double clhs24 = clhs15*clhs23;
const double clhs25 = 1.0*clhs20;
const double clhs26 = clhs15*clhs25;
const double clhs27 = clhs16*clhs22;
const double clhs28 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double clhs29 = N[0]*clhs28;
const double clhs30 = K_darcy*clhs18 + N[0]*clhs19 + clhs14*clhs18 + clhs22*clhs24 - clhs22*clhs26 + clhs27*clhs29;
const double clhs31 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
const double clhs32 = DN(0,0)*clhs31;
const double clhs33 = C(1,3)*DN(0,1);
const double clhs34 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs33;
const double clhs35 = DN(0,1)*clhs34;
const double clhs36 = C(3,5)*DN(0,0);
const double clhs37 = C(4,5)*DN(0,2);
const double clhs38 = C(1,5)*DN(0,1) + clhs36 + clhs37;
const double clhs39 = DN(0,2)*clhs38;
const double clhs40 = DN(0,0)*clhs13;
const double clhs41 = DN(0,1)*clhs40;
const double clhs42 = DN(0,0)*clhs3;
const double clhs43 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs33;
const double clhs44 = DN(0,1)*clhs43;
const double clhs45 = C(3,4)*DN(0,1);
const double clhs46 = C(0,4)*DN(0,0) + clhs37 + clhs45;
const double clhs47 = DN(0,2)*clhs46;
const double clhs48 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs49 = DN(0,0)*clhs48;
const double clhs50 = C(2,3)*DN(0,2) + clhs36 + clhs45;
const double clhs51 = DN(0,1)*clhs50;
const double clhs52 = C(2,5)*DN(0,2);
const double clhs53 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs52;
const double clhs54 = DN(0,2)*clhs53;
const double clhs55 = DN(0,2)*clhs40;
const double clhs56 = DN(0,0)*clhs6;
const double clhs57 = DN(0,1)*clhs46;
const double clhs58 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs52;
const double clhs59 = DN(0,2)*clhs58;
const double clhs60 = clhs16*clhs28;
const double clhs61 = N[0]*clhs60 - N[0] + clhs24 - clhs26;
const double clhs62 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs63 = DN(0,0)*clhs62;
const double clhs64 = C(0,3)*DN(1,0);
const double clhs65 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs64;
const double clhs66 = DN(0,1)*clhs65;
const double clhs67 = C(0,5)*DN(1,0);
const double clhs68 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs67;
const double clhs69 = DN(0,2)*clhs68;
const double clhs70 = DN(1,0)*clhs0;
const double clhs71 = DN(1,1)*clhs3;
const double clhs72 = DN(1,2)*clhs6;
const double clhs73 = DN(0,0)*DN(1,0);
const double clhs74 = N[1]*clhs20 + N[1]*clhs21;
const double clhs75 = clhs13*clhs73 + clhs74;
const double clhs76 = rho*(DN(1,0)*clhs9 + DN(1,1)*clhs10 + DN(1,2)*clhs11);
const double clhs77 = K_darcy*N[1];
const double clhs78 = N[1]*clhs14;
const double clhs79 = clhs76 + clhs77 + clhs78;
const double clhs80 = clhs16*clhs79;
const double clhs81 = N[0]*clhs76 + clhs24*clhs79 - clhs26*clhs79 + clhs29*clhs80;
const double clhs82 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs64;
const double clhs83 = DN(0,0)*clhs82;
const double clhs84 = C(1,3)*DN(1,1);
const double clhs85 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs84;
const double clhs86 = DN(0,1)*clhs85;
const double clhs87 = C(3,5)*DN(1,0);
const double clhs88 = C(4,5)*DN(1,2);
const double clhs89 = C(1,5)*DN(1,1) + clhs87 + clhs88;
const double clhs90 = DN(0,2)*clhs89;
const double clhs91 = DN(1,1)*clhs40;
const double clhs92 = DN(1,0)*clhs3;
const double clhs93 = DN(1,1)*clhs43;
const double clhs94 = DN(1,2)*clhs46;
const double clhs95 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs67;
const double clhs96 = DN(0,0)*clhs95;
const double clhs97 = C(3,4)*DN(1,1);
const double clhs98 = C(2,3)*DN(1,2) + clhs87 + clhs97;
const double clhs99 = DN(0,1)*clhs98;
const double clhs100 = C(2,5)*DN(1,2);
const double clhs101 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs100;
const double clhs102 = DN(0,2)*clhs101;
const double clhs103 = DN(1,2)*clhs40;
const double clhs104 = DN(1,0)*clhs6;
const double clhs105 = DN(1,1)*clhs46;
const double clhs106 = DN(1,2)*clhs58;
const double clhs107 = DN(0,0)*N[1];
const double clhs108 = DN(1,0)*N[0];
const double clhs109 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs110 = DN(0,0)*clhs109;
const double clhs111 = C(0,3)*DN(2,0);
const double clhs112 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs111;
const double clhs113 = DN(0,1)*clhs112;
const double clhs114 = C(0,5)*DN(2,0);
const double clhs115 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs114;
const double clhs116 = DN(0,2)*clhs115;
const double clhs117 = DN(2,0)*clhs0;
const double clhs118 = DN(2,1)*clhs3;
const double clhs119 = DN(2,2)*clhs6;
const double clhs120 = DN(0,0)*DN(2,0);
const double clhs121 = N[2]*clhs20 + N[2]*clhs21;
const double clhs122 = clhs120*clhs13 + clhs121;
const double clhs123 = rho*(DN(2,0)*clhs9 + DN(2,1)*clhs10 + DN(2,2)*clhs11);
const double clhs124 = K_darcy*N[2];
const double clhs125 = N[2]*clhs14;
const double clhs126 = clhs123 + clhs124 + clhs125;
const double clhs127 = clhs126*clhs16;
const double clhs128 = N[0]*clhs123 + clhs126*clhs24 - clhs126*clhs26 + clhs127*clhs29;
const double clhs129 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs111;
const double clhs130 = DN(0,0)*clhs129;
const double clhs131 = C(1,3)*DN(2,1);
const double clhs132 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs131;
const double clhs133 = DN(0,1)*clhs132;
const double clhs134 = C(3,5)*DN(2,0);
const double clhs135 = C(4,5)*DN(2,2);
const double clhs136 = C(1,5)*DN(2,1) + clhs134 + clhs135;
const double clhs137 = DN(0,2)*clhs136;
const double clhs138 = DN(2,1)*clhs40;
const double clhs139 = DN(2,0)*clhs3;
const double clhs140 = DN(2,1)*clhs43;
const double clhs141 = DN(2,2)*clhs46;
const double clhs142 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs114;
const double clhs143 = DN(0,0)*clhs142;
const double clhs144 = C(3,4)*DN(2,1);
const double clhs145 = C(2,3)*DN(2,2) + clhs134 + clhs144;
const double clhs146 = DN(0,1)*clhs145;
const double clhs147 = C(2,5)*DN(2,2);
const double clhs148 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs147;
const double clhs149 = DN(0,2)*clhs148;
const double clhs150 = DN(2,2)*clhs40;
const double clhs151 = DN(2,0)*clhs6;
const double clhs152 = DN(2,1)*clhs46;
const double clhs153 = DN(2,2)*clhs58;
const double clhs154 = DN(0,0)*N[2];
const double clhs155 = DN(2,0)*N[0];
const double clhs156 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs157 = DN(0,0)*clhs156;
const double clhs158 = C(0,3)*DN(3,0);
const double clhs159 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs158;
const double clhs160 = DN(0,1)*clhs159;
const double clhs161 = C(0,5)*DN(3,0);
const double clhs162 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs161;
const double clhs163 = DN(0,2)*clhs162;
const double clhs164 = DN(3,0)*clhs0;
const double clhs165 = DN(3,1)*clhs3;
const double clhs166 = DN(3,2)*clhs6;
const double clhs167 = DN(0,0)*DN(3,0);
const double clhs168 = N[3]*clhs20 + N[3]*clhs21;
const double clhs169 = clhs13*clhs167 + clhs168;
const double clhs170 = rho*(DN(3,0)*clhs9 + DN(3,1)*clhs10 + DN(3,2)*clhs11);
const double clhs171 = K_darcy*N[3];
const double clhs172 = N[3]*clhs14;
const double clhs173 = clhs170 + clhs171 + clhs172;
const double clhs174 = clhs16*clhs173;
const double clhs175 = N[0]*clhs170 + clhs173*clhs24 - clhs173*clhs26 + clhs174*clhs29;
const double clhs176 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs158;
const double clhs177 = DN(0,0)*clhs176;
const double clhs178 = C(1,3)*DN(3,1);
const double clhs179 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs178;
const double clhs180 = DN(0,1)*clhs179;
const double clhs181 = C(3,5)*DN(3,0);
const double clhs182 = C(4,5)*DN(3,2);
const double clhs183 = C(1,5)*DN(3,1) + clhs181 + clhs182;
const double clhs184 = DN(0,2)*clhs183;
const double clhs185 = DN(3,1)*clhs40;
const double clhs186 = DN(3,0)*clhs3;
const double clhs187 = DN(3,1)*clhs43;
const double clhs188 = DN(3,2)*clhs46;
const double clhs189 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs161;
const double clhs190 = DN(0,0)*clhs189;
const double clhs191 = C(3,4)*DN(3,1);
const double clhs192 = C(2,3)*DN(3,2) + clhs181 + clhs191;
const double clhs193 = DN(0,1)*clhs192;
const double clhs194 = C(2,5)*DN(3,2);
const double clhs195 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs194;
const double clhs196 = DN(0,2)*clhs195;
const double clhs197 = DN(3,2)*clhs40;
const double clhs198 = DN(3,0)*clhs6;
const double clhs199 = DN(3,1)*clhs46;
const double clhs200 = DN(3,2)*clhs58;
const double clhs201 = DN(0,0)*N[3];
const double clhs202 = DN(3,0)*N[0];
const double clhs203 = DN(0,0)*clhs34;
const double clhs204 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs205 = DN(0,1)*clhs204;
const double clhs206 = C(1,4)*DN(0,1);
const double clhs207 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs206;
const double clhs208 = DN(0,2)*clhs207;
const double clhs209 = pow(DN(0,1), 2);
const double clhs210 = DN(0,0)*clhs50;
const double clhs211 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs206;
const double clhs212 = DN(0,1)*clhs211;
const double clhs213 = C(2,4)*DN(0,2);
const double clhs214 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs213;
const double clhs215 = DN(0,2)*clhs214;
const double clhs216 = DN(0,1)*clhs13;
const double clhs217 = DN(0,2)*clhs216;
const double clhs218 = DN(0,0)*clhs38;
const double clhs219 = DN(0,1)*clhs207;
const double clhs220 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs213;
const double clhs221 = DN(0,2)*clhs220;
const double clhs222 = DN(0,0)*clhs65;
const double clhs223 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs84;
const double clhs224 = DN(0,1)*clhs223;
const double clhs225 = C(0,4)*DN(1,0) + clhs88 + clhs97;
const double clhs226 = DN(0,2)*clhs225;
const double clhs227 = DN(1,0)*clhs216;
const double clhs228 = DN(1,0)*clhs31;
const double clhs229 = DN(1,1)*clhs34;
const double clhs230 = DN(1,2)*clhs38;
const double clhs231 = DN(0,0)*clhs85;
const double clhs232 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs233 = DN(0,1)*clhs232;
const double clhs234 = C(1,4)*DN(1,1);
const double clhs235 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs234;
const double clhs236 = DN(0,2)*clhs235;
const double clhs237 = DN(0,1)*DN(1,1);
const double clhs238 = clhs13*clhs237;
const double clhs239 = DN(1,0)*clhs34;
const double clhs240 = DN(1,1)*clhs204;
const double clhs241 = DN(1,2)*clhs207;
const double clhs242 = clhs74 + clhs81;
const double clhs243 = DN(0,0)*clhs98;
const double clhs244 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs234;
const double clhs245 = DN(0,1)*clhs244;
const double clhs246 = C(2,4)*DN(1,2);
const double clhs247 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs246;
const double clhs248 = DN(0,2)*clhs247;
const double clhs249 = DN(1,2)*clhs216;
const double clhs250 = DN(1,0)*clhs38;
const double clhs251 = DN(1,1)*clhs207;
const double clhs252 = DN(1,2)*clhs220;
const double clhs253 = DN(0,1)*N[1];
const double clhs254 = DN(1,1)*N[0];
const double clhs255 = DN(0,0)*clhs112;
const double clhs256 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs131;
const double clhs257 = DN(0,1)*clhs256;
const double clhs258 = C(0,4)*DN(2,0) + clhs135 + clhs144;
const double clhs259 = DN(0,2)*clhs258;
const double clhs260 = DN(2,0)*clhs216;
const double clhs261 = DN(2,0)*clhs31;
const double clhs262 = DN(2,1)*clhs34;
const double clhs263 = DN(2,2)*clhs38;
const double clhs264 = DN(0,0)*clhs132;
const double clhs265 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs266 = DN(0,1)*clhs265;
const double clhs267 = C(1,4)*DN(2,1);
const double clhs268 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs267;
const double clhs269 = DN(0,2)*clhs268;
const double clhs270 = DN(0,1)*DN(2,1);
const double clhs271 = clhs13*clhs270;
const double clhs272 = DN(2,0)*clhs34;
const double clhs273 = DN(2,1)*clhs204;
const double clhs274 = DN(2,2)*clhs207;
const double clhs275 = clhs121 + clhs128;
const double clhs276 = DN(0,0)*clhs145;
const double clhs277 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs267;
const double clhs278 = DN(0,1)*clhs277;
const double clhs279 = C(2,4)*DN(2,2);
const double clhs280 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs279;
const double clhs281 = DN(0,2)*clhs280;
const double clhs282 = DN(2,2)*clhs216;
const double clhs283 = DN(2,0)*clhs38;
const double clhs284 = DN(2,1)*clhs207;
const double clhs285 = DN(2,2)*clhs220;
const double clhs286 = DN(0,1)*N[2];
const double clhs287 = DN(2,1)*N[0];
const double clhs288 = DN(0,0)*clhs159;
const double clhs289 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs178;
const double clhs290 = DN(0,1)*clhs289;
const double clhs291 = C(0,4)*DN(3,0) + clhs182 + clhs191;
const double clhs292 = DN(0,2)*clhs291;
const double clhs293 = DN(3,0)*clhs216;
const double clhs294 = DN(3,0)*clhs31;
const double clhs295 = DN(3,1)*clhs34;
const double clhs296 = DN(3,2)*clhs38;
const double clhs297 = DN(0,0)*clhs179;
const double clhs298 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs299 = DN(0,1)*clhs298;
const double clhs300 = C(1,4)*DN(3,1);
const double clhs301 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs300;
const double clhs302 = DN(0,2)*clhs301;
const double clhs303 = DN(0,1)*DN(3,1);
const double clhs304 = clhs13*clhs303;
const double clhs305 = DN(3,0)*clhs34;
const double clhs306 = DN(3,1)*clhs204;
const double clhs307 = DN(3,2)*clhs207;
const double clhs308 = clhs168 + clhs175;
const double clhs309 = DN(0,0)*clhs192;
const double clhs310 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs300;
const double clhs311 = DN(0,1)*clhs310;
const double clhs312 = C(2,4)*DN(3,2);
const double clhs313 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs312;
const double clhs314 = DN(0,2)*clhs313;
const double clhs315 = DN(3,2)*clhs216;
const double clhs316 = DN(3,0)*clhs38;
const double clhs317 = DN(3,1)*clhs207;
const double clhs318 = DN(3,2)*clhs220;
const double clhs319 = DN(0,1)*N[3];
const double clhs320 = DN(3,1)*N[0];
const double clhs321 = DN(0,0)*clhs53;
const double clhs322 = DN(0,1)*clhs214;
const double clhs323 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs324 = DN(0,2)*clhs323;
const double clhs325 = pow(DN(0,2), 2);
const double clhs326 = DN(0,0)*clhs68;
const double clhs327 = DN(0,1)*clhs225;
const double clhs328 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs100;
const double clhs329 = DN(0,2)*clhs328;
const double clhs330 = DN(0,2)*clhs13;
const double clhs331 = DN(1,0)*clhs330;
const double clhs332 = DN(1,0)*clhs48;
const double clhs333 = DN(1,1)*clhs50;
const double clhs334 = DN(1,2)*clhs53;
const double clhs335 = DN(0,0)*clhs89;
const double clhs336 = DN(0,1)*clhs235;
const double clhs337 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs246;
const double clhs338 = DN(0,2)*clhs337;
const double clhs339 = DN(1,1)*clhs330;
const double clhs340 = DN(1,0)*clhs50;
const double clhs341 = DN(1,1)*clhs211;
const double clhs342 = DN(1,2)*clhs214;
const double clhs343 = DN(0,0)*clhs101;
const double clhs344 = DN(0,1)*clhs247;
const double clhs345 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs346 = DN(0,2)*clhs345;
const double clhs347 = DN(0,2)*DN(1,2);
const double clhs348 = clhs13*clhs347;
const double clhs349 = DN(1,0)*clhs53;
const double clhs350 = DN(1,1)*clhs214;
const double clhs351 = DN(1,2)*clhs323;
const double clhs352 = DN(0,2)*N[1];
const double clhs353 = DN(1,2)*N[0];
const double clhs354 = DN(0,0)*clhs115;
const double clhs355 = DN(0,1)*clhs258;
const double clhs356 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs147;
const double clhs357 = DN(0,2)*clhs356;
const double clhs358 = DN(2,0)*clhs330;
const double clhs359 = DN(2,0)*clhs48;
const double clhs360 = DN(2,1)*clhs50;
const double clhs361 = DN(2,2)*clhs53;
const double clhs362 = DN(0,0)*clhs136;
const double clhs363 = DN(0,1)*clhs268;
const double clhs364 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs279;
const double clhs365 = DN(0,2)*clhs364;
const double clhs366 = DN(2,1)*clhs330;
const double clhs367 = DN(2,0)*clhs50;
const double clhs368 = DN(2,1)*clhs211;
const double clhs369 = DN(2,2)*clhs214;
const double clhs370 = DN(0,0)*clhs148;
const double clhs371 = DN(0,1)*clhs280;
const double clhs372 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs373 = DN(0,2)*clhs372;
const double clhs374 = DN(0,2)*DN(2,2);
const double clhs375 = clhs13*clhs374;
const double clhs376 = DN(2,0)*clhs53;
const double clhs377 = DN(2,1)*clhs214;
const double clhs378 = DN(2,2)*clhs323;
const double clhs379 = DN(0,2)*N[2];
const double clhs380 = DN(2,2)*N[0];
const double clhs381 = DN(0,0)*clhs162;
const double clhs382 = DN(0,1)*clhs291;
const double clhs383 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs194;
const double clhs384 = DN(0,2)*clhs383;
const double clhs385 = DN(3,0)*clhs330;
const double clhs386 = DN(3,0)*clhs48;
const double clhs387 = DN(3,1)*clhs50;
const double clhs388 = DN(3,2)*clhs53;
const double clhs389 = DN(0,0)*clhs183;
const double clhs390 = DN(0,1)*clhs301;
const double clhs391 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs312;
const double clhs392 = DN(0,2)*clhs391;
const double clhs393 = DN(3,1)*clhs330;
const double clhs394 = DN(3,0)*clhs50;
const double clhs395 = DN(3,1)*clhs211;
const double clhs396 = DN(3,2)*clhs214;
const double clhs397 = DN(0,0)*clhs195;
const double clhs398 = DN(0,1)*clhs313;
const double clhs399 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs400 = DN(0,2)*clhs399;
const double clhs401 = DN(0,2)*DN(3,2);
const double clhs402 = clhs13*clhs401;
const double clhs403 = DN(3,0)*clhs53;
const double clhs404 = DN(3,1)*clhs214;
const double clhs405 = DN(3,2)*clhs323;
const double clhs406 = DN(0,2)*N[3];
const double clhs407 = DN(3,2)*N[0];
const double clhs408 = N[0] + clhs15*(1.0*clhs21 + clhs23 + clhs25);
const double clhs409 = clhs16*(clhs237 + clhs347 + clhs73);
const double clhs410 = clhs16*(clhs120 + clhs270 + clhs374);
const double clhs411 = clhs16*(clhs167 + clhs303 + clhs401);
const double clhs412 = clhs16*clhs76;
const double clhs413 = clhs16*clhs77;
const double clhs414 = N[1]*clhs28;
const double clhs415 = N[1]*clhs19 + clhs22*clhs412 - clhs22*clhs413 + clhs27*clhs414;
const double clhs416 = DN(1,0)*clhs62;
const double clhs417 = DN(1,1)*clhs65;
const double clhs418 = DN(1,2)*clhs68;
const double clhs419 = pow(DN(1,0), 2);
const double clhs420 = pow(N[1], 2);
const double clhs421 = K_darcy*clhs420 + N[1]*clhs76 + clhs14*clhs420 + clhs414*clhs80 + clhs76*clhs80 - clhs77*clhs80;
const double clhs422 = DN(1,0)*clhs82;
const double clhs423 = DN(1,1)*clhs85;
const double clhs424 = DN(1,2)*clhs89;
const double clhs425 = DN(1,0)*clhs13;
const double clhs426 = DN(1,1)*clhs425;
const double clhs427 = DN(1,0)*clhs65;
const double clhs428 = DN(1,1)*clhs223;
const double clhs429 = DN(1,2)*clhs225;
const double clhs430 = DN(1,0)*clhs95;
const double clhs431 = DN(1,1)*clhs98;
const double clhs432 = DN(1,2)*clhs101;
const double clhs433 = DN(1,2)*clhs425;
const double clhs434 = DN(1,0)*clhs68;
const double clhs435 = DN(1,1)*clhs225;
const double clhs436 = DN(1,2)*clhs328;
const double clhs437 = N[1]*clhs60 - N[1] + clhs412 - clhs413;
const double clhs438 = DN(1,0)*clhs109;
const double clhs439 = DN(1,1)*clhs112;
const double clhs440 = DN(1,2)*clhs115;
const double clhs441 = DN(2,0)*clhs62;
const double clhs442 = DN(2,1)*clhs65;
const double clhs443 = DN(2,2)*clhs68;
const double clhs444 = DN(1,0)*DN(2,0);
const double clhs445 = N[2]*clhs77 + N[2]*clhs78;
const double clhs446 = clhs13*clhs444 + clhs445;
const double clhs447 = N[1]*clhs123 + clhs127*clhs414 + clhs127*clhs76 - clhs127*clhs77;
const double clhs448 = DN(1,0)*clhs129;
const double clhs449 = DN(1,1)*clhs132;
const double clhs450 = DN(1,2)*clhs136;
const double clhs451 = DN(2,1)*clhs425;
const double clhs452 = DN(2,0)*clhs65;
const double clhs453 = DN(2,1)*clhs223;
const double clhs454 = DN(2,2)*clhs225;
const double clhs455 = DN(1,0)*clhs142;
const double clhs456 = DN(1,1)*clhs145;
const double clhs457 = DN(1,2)*clhs148;
const double clhs458 = DN(2,2)*clhs425;
const double clhs459 = DN(2,0)*clhs68;
const double clhs460 = DN(2,1)*clhs225;
const double clhs461 = DN(2,2)*clhs328;
const double clhs462 = DN(1,0)*N[2];
const double clhs463 = DN(2,0)*N[1];
const double clhs464 = DN(1,0)*clhs156;
const double clhs465 = DN(1,1)*clhs159;
const double clhs466 = DN(1,2)*clhs162;
const double clhs467 = DN(3,0)*clhs62;
const double clhs468 = DN(3,1)*clhs65;
const double clhs469 = DN(3,2)*clhs68;
const double clhs470 = DN(1,0)*DN(3,0);
const double clhs471 = N[3]*clhs77 + N[3]*clhs78;
const double clhs472 = clhs13*clhs470 + clhs471;
const double clhs473 = N[1]*clhs170 + clhs174*clhs414 + clhs174*clhs76 - clhs174*clhs77;
const double clhs474 = DN(1,0)*clhs176;
const double clhs475 = DN(1,1)*clhs179;
const double clhs476 = DN(1,2)*clhs183;
const double clhs477 = DN(3,1)*clhs425;
const double clhs478 = DN(3,0)*clhs65;
const double clhs479 = DN(3,1)*clhs223;
const double clhs480 = DN(3,2)*clhs225;
const double clhs481 = DN(1,0)*clhs189;
const double clhs482 = DN(1,1)*clhs192;
const double clhs483 = DN(1,2)*clhs195;
const double clhs484 = DN(3,2)*clhs425;
const double clhs485 = DN(3,0)*clhs68;
const double clhs486 = DN(3,1)*clhs225;
const double clhs487 = DN(3,2)*clhs328;
const double clhs488 = DN(1,0)*N[3];
const double clhs489 = DN(3,0)*N[1];
const double clhs490 = clhs415 + clhs74;
const double clhs491 = DN(1,0)*clhs85;
const double clhs492 = DN(1,1)*clhs232;
const double clhs493 = DN(1,2)*clhs235;
const double clhs494 = pow(DN(1,1), 2);
const double clhs495 = DN(1,0)*clhs98;
const double clhs496 = DN(1,1)*clhs244;
const double clhs497 = DN(1,2)*clhs247;
const double clhs498 = DN(1,1)*clhs13;
const double clhs499 = DN(1,2)*clhs498;
const double clhs500 = DN(1,0)*clhs89;
const double clhs501 = DN(1,1)*clhs235;
const double clhs502 = DN(1,2)*clhs337;
const double clhs503 = DN(1,0)*clhs112;
const double clhs504 = DN(1,1)*clhs256;
const double clhs505 = DN(1,2)*clhs258;
const double clhs506 = DN(2,0)*clhs498;
const double clhs507 = DN(2,0)*clhs82;
const double clhs508 = DN(2,1)*clhs85;
const double clhs509 = DN(2,2)*clhs89;
const double clhs510 = DN(1,0)*clhs132;
const double clhs511 = DN(1,1)*clhs265;
const double clhs512 = DN(1,2)*clhs268;
const double clhs513 = DN(1,1)*DN(2,1);
const double clhs514 = clhs13*clhs513;
const double clhs515 = DN(2,0)*clhs85;
const double clhs516 = DN(2,1)*clhs232;
const double clhs517 = DN(2,2)*clhs235;
const double clhs518 = clhs445 + clhs447;
const double clhs519 = DN(1,0)*clhs145;
const double clhs520 = DN(1,1)*clhs277;
const double clhs521 = DN(1,2)*clhs280;
const double clhs522 = DN(2,2)*clhs498;
const double clhs523 = DN(2,0)*clhs89;
const double clhs524 = DN(2,1)*clhs235;
const double clhs525 = DN(2,2)*clhs337;
const double clhs526 = DN(1,1)*N[2];
const double clhs527 = DN(2,1)*N[1];
const double clhs528 = DN(1,0)*clhs159;
const double clhs529 = DN(1,1)*clhs289;
const double clhs530 = DN(1,2)*clhs291;
const double clhs531 = DN(3,0)*clhs498;
const double clhs532 = DN(3,0)*clhs82;
const double clhs533 = DN(3,1)*clhs85;
const double clhs534 = DN(3,2)*clhs89;
const double clhs535 = DN(1,0)*clhs179;
const double clhs536 = DN(1,1)*clhs298;
const double clhs537 = DN(1,2)*clhs301;
const double clhs538 = DN(1,1)*DN(3,1);
const double clhs539 = clhs13*clhs538;
const double clhs540 = DN(3,0)*clhs85;
const double clhs541 = DN(3,1)*clhs232;
const double clhs542 = DN(3,2)*clhs235;
const double clhs543 = clhs471 + clhs473;
const double clhs544 = DN(1,0)*clhs192;
const double clhs545 = DN(1,1)*clhs310;
const double clhs546 = DN(1,2)*clhs313;
const double clhs547 = DN(3,2)*clhs498;
const double clhs548 = DN(3,0)*clhs89;
const double clhs549 = DN(3,1)*clhs235;
const double clhs550 = DN(3,2)*clhs337;
const double clhs551 = DN(1,1)*N[3];
const double clhs552 = DN(3,1)*N[1];
const double clhs553 = DN(1,0)*clhs101;
const double clhs554 = DN(1,1)*clhs247;
const double clhs555 = DN(1,2)*clhs345;
const double clhs556 = pow(DN(1,2), 2);
const double clhs557 = DN(1,0)*clhs115;
const double clhs558 = DN(1,1)*clhs258;
const double clhs559 = DN(1,2)*clhs356;
const double clhs560 = DN(1,2)*clhs13;
const double clhs561 = DN(2,0)*clhs560;
const double clhs562 = DN(2,0)*clhs95;
const double clhs563 = DN(2,1)*clhs98;
const double clhs564 = DN(2,2)*clhs101;
const double clhs565 = DN(1,0)*clhs136;
const double clhs566 = DN(1,1)*clhs268;
const double clhs567 = DN(1,2)*clhs364;
const double clhs568 = DN(2,1)*clhs560;
const double clhs569 = DN(2,0)*clhs98;
const double clhs570 = DN(2,1)*clhs244;
const double clhs571 = DN(2,2)*clhs247;
const double clhs572 = DN(1,0)*clhs148;
const double clhs573 = DN(1,1)*clhs280;
const double clhs574 = DN(1,2)*clhs372;
const double clhs575 = DN(1,2)*DN(2,2);
const double clhs576 = clhs13*clhs575;
const double clhs577 = DN(2,0)*clhs101;
const double clhs578 = DN(2,1)*clhs247;
const double clhs579 = DN(2,2)*clhs345;
const double clhs580 = DN(1,2)*N[2];
const double clhs581 = DN(2,2)*N[1];
const double clhs582 = DN(1,0)*clhs162;
const double clhs583 = DN(1,1)*clhs291;
const double clhs584 = DN(1,2)*clhs383;
const double clhs585 = DN(3,0)*clhs560;
const double clhs586 = DN(3,0)*clhs95;
const double clhs587 = DN(3,1)*clhs98;
const double clhs588 = DN(3,2)*clhs101;
const double clhs589 = DN(1,0)*clhs183;
const double clhs590 = DN(1,1)*clhs301;
const double clhs591 = DN(1,2)*clhs391;
const double clhs592 = DN(3,1)*clhs560;
const double clhs593 = DN(3,0)*clhs98;
const double clhs594 = DN(3,1)*clhs244;
const double clhs595 = DN(3,2)*clhs247;
const double clhs596 = DN(1,0)*clhs195;
const double clhs597 = DN(1,1)*clhs313;
const double clhs598 = DN(1,2)*clhs399;
const double clhs599 = DN(1,2)*DN(3,2);
const double clhs600 = clhs13*clhs599;
const double clhs601 = DN(3,0)*clhs101;
const double clhs602 = DN(3,1)*clhs247;
const double clhs603 = DN(3,2)*clhs345;
const double clhs604 = DN(1,2)*N[3];
const double clhs605 = DN(3,2)*N[1];
const double clhs606 = N[1] + clhs15*(1.0*clhs76 + 1.0*clhs77 + 1.0*clhs78);
const double clhs607 = clhs16*(clhs444 + clhs513 + clhs575);
const double clhs608 = clhs16*(clhs470 + clhs538 + clhs599);
const double clhs609 = N[2]*clhs28;
const double clhs610 = N[2]*clhs19 + clhs123*clhs27 - clhs124*clhs27 + clhs27*clhs609;
const double clhs611 = clhs124*clhs16;
const double clhs612 = clhs123*clhs16;
const double clhs613 = N[2]*clhs76 + clhs123*clhs80 - clhs124*clhs80 + clhs609*clhs80;
const double clhs614 = DN(2,0)*clhs109;
const double clhs615 = DN(2,1)*clhs112;
const double clhs616 = DN(2,2)*clhs115;
const double clhs617 = pow(DN(2,0), 2);
const double clhs618 = pow(N[2], 2);
const double clhs619 = K_darcy*clhs618 + N[2]*clhs123 + clhs123*clhs127 - clhs124*clhs127 + clhs127*clhs609 + clhs14*clhs618;
const double clhs620 = DN(2,0)*clhs129;
const double clhs621 = DN(2,1)*clhs132;
const double clhs622 = DN(2,2)*clhs136;
const double clhs623 = DN(2,0)*clhs13;
const double clhs624 = DN(2,1)*clhs623;
const double clhs625 = DN(2,0)*clhs112;
const double clhs626 = DN(2,1)*clhs256;
const double clhs627 = DN(2,2)*clhs258;
const double clhs628 = DN(2,0)*clhs142;
const double clhs629 = DN(2,1)*clhs145;
const double clhs630 = DN(2,2)*clhs148;
const double clhs631 = DN(2,2)*clhs623;
const double clhs632 = DN(2,0)*clhs115;
const double clhs633 = DN(2,1)*clhs258;
const double clhs634 = DN(2,2)*clhs356;
const double clhs635 = N[2]*clhs60 - N[2] - clhs611 + clhs612;
const double clhs636 = DN(2,0)*clhs156;
const double clhs637 = DN(2,1)*clhs159;
const double clhs638 = DN(2,2)*clhs162;
const double clhs639 = DN(3,0)*clhs109;
const double clhs640 = DN(3,1)*clhs112;
const double clhs641 = DN(3,2)*clhs115;
const double clhs642 = DN(2,0)*DN(3,0);
const double clhs643 = N[3]*clhs124 + N[3]*clhs125;
const double clhs644 = clhs13*clhs642 + clhs643;
const double clhs645 = N[2]*clhs170 + clhs123*clhs174 - clhs124*clhs174 + clhs174*clhs609;
const double clhs646 = DN(2,0)*clhs176;
const double clhs647 = DN(2,1)*clhs179;
const double clhs648 = DN(2,2)*clhs183;
const double clhs649 = DN(3,1)*clhs623;
const double clhs650 = DN(3,0)*clhs112;
const double clhs651 = DN(3,1)*clhs256;
const double clhs652 = DN(3,2)*clhs258;
const double clhs653 = DN(2,0)*clhs189;
const double clhs654 = DN(2,1)*clhs192;
const double clhs655 = DN(2,2)*clhs195;
const double clhs656 = DN(3,2)*clhs623;
const double clhs657 = DN(3,0)*clhs115;
const double clhs658 = DN(3,1)*clhs258;
const double clhs659 = DN(3,2)*clhs356;
const double clhs660 = DN(2,0)*N[3];
const double clhs661 = DN(3,0)*N[2];
const double clhs662 = clhs121 + clhs610;
const double clhs663 = clhs445 + clhs613;
const double clhs664 = DN(2,0)*clhs132;
const double clhs665 = DN(2,1)*clhs265;
const double clhs666 = DN(2,2)*clhs268;
const double clhs667 = pow(DN(2,1), 2);
const double clhs668 = DN(2,0)*clhs145;
const double clhs669 = DN(2,1)*clhs277;
const double clhs670 = DN(2,2)*clhs280;
const double clhs671 = DN(2,1)*clhs13;
const double clhs672 = DN(2,2)*clhs671;
const double clhs673 = DN(2,0)*clhs136;
const double clhs674 = DN(2,1)*clhs268;
const double clhs675 = DN(2,2)*clhs364;
const double clhs676 = DN(2,0)*clhs159;
const double clhs677 = DN(2,1)*clhs289;
const double clhs678 = DN(2,2)*clhs291;
const double clhs679 = DN(3,0)*clhs671;
const double clhs680 = DN(3,0)*clhs129;
const double clhs681 = DN(3,1)*clhs132;
const double clhs682 = DN(3,2)*clhs136;
const double clhs683 = DN(2,0)*clhs179;
const double clhs684 = DN(2,1)*clhs298;
const double clhs685 = DN(2,2)*clhs301;
const double clhs686 = DN(2,1)*DN(3,1);
const double clhs687 = clhs13*clhs686;
const double clhs688 = DN(3,0)*clhs132;
const double clhs689 = DN(3,1)*clhs265;
const double clhs690 = DN(3,2)*clhs268;
const double clhs691 = clhs643 + clhs645;
const double clhs692 = DN(2,0)*clhs192;
const double clhs693 = DN(2,1)*clhs310;
const double clhs694 = DN(2,2)*clhs313;
const double clhs695 = DN(3,2)*clhs671;
const double clhs696 = DN(3,0)*clhs136;
const double clhs697 = DN(3,1)*clhs268;
const double clhs698 = DN(3,2)*clhs364;
const double clhs699 = DN(2,1)*N[3];
const double clhs700 = DN(3,1)*N[2];
const double clhs701 = DN(2,0)*clhs148;
const double clhs702 = DN(2,1)*clhs280;
const double clhs703 = DN(2,2)*clhs372;
const double clhs704 = pow(DN(2,2), 2);
const double clhs705 = DN(2,0)*clhs162;
const double clhs706 = DN(2,1)*clhs291;
const double clhs707 = DN(2,2)*clhs383;
const double clhs708 = DN(2,2)*clhs13;
const double clhs709 = DN(3,0)*clhs708;
const double clhs710 = DN(3,0)*clhs142;
const double clhs711 = DN(3,1)*clhs145;
const double clhs712 = DN(3,2)*clhs148;
const double clhs713 = DN(2,0)*clhs183;
const double clhs714 = DN(2,1)*clhs301;
const double clhs715 = DN(2,2)*clhs391;
const double clhs716 = DN(3,1)*clhs708;
const double clhs717 = DN(3,0)*clhs145;
const double clhs718 = DN(3,1)*clhs277;
const double clhs719 = DN(3,2)*clhs280;
const double clhs720 = DN(2,0)*clhs195;
const double clhs721 = DN(2,1)*clhs313;
const double clhs722 = DN(2,2)*clhs399;
const double clhs723 = DN(2,2)*DN(3,2);
const double clhs724 = clhs13*clhs723;
const double clhs725 = DN(3,0)*clhs148;
const double clhs726 = DN(3,1)*clhs280;
const double clhs727 = DN(3,2)*clhs372;
const double clhs728 = DN(2,2)*N[3];
const double clhs729 = DN(3,2)*N[2];
const double clhs730 = N[2] + clhs15*(1.0*clhs123 + 1.0*clhs124 + 1.0*clhs125);
const double clhs731 = clhs16*(clhs642 + clhs686 + clhs723);
const double clhs732 = N[3]*clhs28;
const double clhs733 = N[3]*clhs19 + clhs170*clhs27 - clhs171*clhs27 + clhs27*clhs732;
const double clhs734 = clhs16*clhs171;
const double clhs735 = clhs16*clhs170;
const double clhs736 = N[3]*clhs76 + clhs170*clhs80 - clhs171*clhs80 + clhs732*clhs80;
const double clhs737 = N[3]*clhs123 + clhs127*clhs170 - clhs127*clhs171 + clhs127*clhs732;
const double clhs738 = DN(3,0)*clhs156;
const double clhs739 = DN(3,1)*clhs159;
const double clhs740 = DN(3,2)*clhs162;
const double clhs741 = pow(DN(3,0), 2);
const double clhs742 = pow(N[3], 2);
const double clhs743 = K_darcy*clhs742 + N[3]*clhs170 + clhs14*clhs742 + clhs170*clhs174 - clhs171*clhs174 + clhs174*clhs732;
const double clhs744 = DN(3,0)*clhs176;
const double clhs745 = DN(3,1)*clhs179;
const double clhs746 = DN(3,2)*clhs183;
const double clhs747 = DN(3,0)*clhs13;
const double clhs748 = DN(3,1)*clhs747;
const double clhs749 = DN(3,0)*clhs159;
const double clhs750 = DN(3,1)*clhs289;
const double clhs751 = DN(3,2)*clhs291;
const double clhs752 = DN(3,0)*clhs189;
const double clhs753 = DN(3,1)*clhs192;
const double clhs754 = DN(3,2)*clhs195;
const double clhs755 = DN(3,2)*clhs747;
const double clhs756 = DN(3,0)*clhs162;
const double clhs757 = DN(3,1)*clhs291;
const double clhs758 = DN(3,2)*clhs383;
const double clhs759 = N[3]*clhs60 - N[3] - clhs734 + clhs735;
const double clhs760 = clhs168 + clhs733;
const double clhs761 = clhs471 + clhs736;
const double clhs762 = clhs643 + clhs737;
const double clhs763 = DN(3,0)*clhs179;
const double clhs764 = DN(3,1)*clhs298;
const double clhs765 = DN(3,2)*clhs301;
const double clhs766 = pow(DN(3,1), 2);
const double clhs767 = DN(3,0)*clhs192;
const double clhs768 = DN(3,1)*clhs310;
const double clhs769 = DN(3,2)*clhs313;
const double clhs770 = DN(3,1)*DN(3,2)*clhs13;
const double clhs771 = DN(3,0)*clhs183;
const double clhs772 = DN(3,1)*clhs301;
const double clhs773 = DN(3,2)*clhs391;
const double clhs774 = DN(3,0)*clhs195;
const double clhs775 = DN(3,1)*clhs313;
const double clhs776 = DN(3,2)*clhs399;
const double clhs777 = pow(DN(3,2), 2);
const double clhs778 = N[3] + clhs15*(1.0*clhs170 + 1.0*clhs171 + 1.0*clhs172);
lhs(0,0)=-clhs1*clhs17 + clhs1 + clhs13*clhs8 - clhs17*clhs4 - clhs17*clhs7 + clhs30 + clhs4 + clhs7;
lhs(0,1)=-clhs17*clhs42 - clhs17*clhs44 - clhs17*clhs47 + clhs32 + clhs35 + clhs39 + clhs41;
lhs(0,2)=-clhs17*clhs56 - clhs17*clhs57 - clhs17*clhs59 + clhs49 + clhs51 + clhs54 + clhs55;
lhs(0,3)=DN(0,0)*clhs61;
lhs(0,4)=-clhs17*clhs70 - clhs17*clhs71 - clhs17*clhs72 + clhs63 + clhs66 + clhs69 + clhs75 + clhs81;
lhs(0,5)=-clhs17*clhs92 - clhs17*clhs93 - clhs17*clhs94 + clhs83 + clhs86 + clhs90 + clhs91;
lhs(0,6)=clhs102 + clhs103 - clhs104*clhs17 - clhs105*clhs17 - clhs106*clhs17 + clhs96 + clhs99;
lhs(0,7)=DN(1,0)*clhs24 - DN(1,0)*clhs26 - clhs107 + clhs108*clhs60;
lhs(0,8)=clhs110 + clhs113 + clhs116 - clhs117*clhs17 - clhs118*clhs17 - clhs119*clhs17 + clhs122 + clhs128;
lhs(0,9)=clhs130 + clhs133 + clhs137 + clhs138 - clhs139*clhs17 - clhs140*clhs17 - clhs141*clhs17;
lhs(0,10)=clhs143 + clhs146 + clhs149 + clhs150 - clhs151*clhs17 - clhs152*clhs17 - clhs153*clhs17;
lhs(0,11)=DN(2,0)*clhs24 - DN(2,0)*clhs26 - clhs154 + clhs155*clhs60;
lhs(0,12)=clhs157 + clhs160 + clhs163 - clhs164*clhs17 - clhs165*clhs17 - clhs166*clhs17 + clhs169 + clhs175;
lhs(0,13)=-clhs17*clhs186 - clhs17*clhs187 - clhs17*clhs188 + clhs177 + clhs180 + clhs184 + clhs185;
lhs(0,14)=-clhs17*clhs198 - clhs17*clhs199 - clhs17*clhs200 + clhs190 + clhs193 + clhs196 + clhs197;
lhs(0,15)=DN(3,0)*clhs24 - DN(3,0)*clhs26 - clhs201 + clhs202*clhs60;
lhs(1,0)=-clhs17*clhs32 - clhs17*clhs35 - clhs17*clhs39 + clhs41 + clhs42 + clhs44 + clhs47;
lhs(1,1)=clhs13*clhs209 - clhs17*clhs203 - clhs17*clhs205 - clhs17*clhs208 + clhs203 + clhs205 + clhs208 + clhs30;
lhs(1,2)=-clhs17*clhs218 - clhs17*clhs219 - clhs17*clhs221 + clhs210 + clhs212 + clhs215 + clhs217;
lhs(1,3)=DN(0,1)*clhs61;
lhs(1,4)=-clhs17*clhs228 - clhs17*clhs229 - clhs17*clhs230 + clhs222 + clhs224 + clhs226 + clhs227;
lhs(1,5)=-clhs17*clhs239 - clhs17*clhs240 - clhs17*clhs241 + clhs231 + clhs233 + clhs236 + clhs238 + clhs242;
lhs(1,6)=-clhs17*clhs250 - clhs17*clhs251 - clhs17*clhs252 + clhs243 + clhs245 + clhs248 + clhs249;
lhs(1,7)=DN(1,1)*clhs24 - DN(1,1)*clhs26 - clhs253 + clhs254*clhs60;
lhs(1,8)=-clhs17*clhs261 - clhs17*clhs262 - clhs17*clhs263 + clhs255 + clhs257 + clhs259 + clhs260;
lhs(1,9)=-clhs17*clhs272 - clhs17*clhs273 - clhs17*clhs274 + clhs264 + clhs266 + clhs269 + clhs271 + clhs275;
lhs(1,10)=-clhs17*clhs283 - clhs17*clhs284 - clhs17*clhs285 + clhs276 + clhs278 + clhs281 + clhs282;
lhs(1,11)=DN(2,1)*clhs24 - DN(2,1)*clhs26 - clhs286 + clhs287*clhs60;
lhs(1,12)=-clhs17*clhs294 - clhs17*clhs295 - clhs17*clhs296 + clhs288 + clhs290 + clhs292 + clhs293;
lhs(1,13)=-clhs17*clhs305 - clhs17*clhs306 - clhs17*clhs307 + clhs297 + clhs299 + clhs302 + clhs304 + clhs308;
lhs(1,14)=-clhs17*clhs316 - clhs17*clhs317 - clhs17*clhs318 + clhs309 + clhs311 + clhs314 + clhs315;
lhs(1,15)=DN(3,1)*clhs24 - DN(3,1)*clhs26 - clhs319 + clhs320*clhs60;
lhs(2,0)=-clhs17*clhs49 - clhs17*clhs51 - clhs17*clhs54 + clhs55 + clhs56 + clhs57 + clhs59;
lhs(2,1)=-clhs17*clhs210 - clhs17*clhs212 - clhs17*clhs215 + clhs217 + clhs218 + clhs219 + clhs221;
lhs(2,2)=clhs13*clhs325 - clhs17*clhs321 - clhs17*clhs322 - clhs17*clhs324 + clhs30 + clhs321 + clhs322 + clhs324;
lhs(2,3)=DN(0,2)*clhs61;
lhs(2,4)=-clhs17*clhs332 - clhs17*clhs333 - clhs17*clhs334 + clhs326 + clhs327 + clhs329 + clhs331;
lhs(2,5)=-clhs17*clhs340 - clhs17*clhs341 - clhs17*clhs342 + clhs335 + clhs336 + clhs338 + clhs339;
lhs(2,6)=-clhs17*clhs349 - clhs17*clhs350 - clhs17*clhs351 + clhs242 + clhs343 + clhs344 + clhs346 + clhs348;
lhs(2,7)=DN(1,2)*clhs24 - DN(1,2)*clhs26 - clhs352 + clhs353*clhs60;
lhs(2,8)=-clhs17*clhs359 - clhs17*clhs360 - clhs17*clhs361 + clhs354 + clhs355 + clhs357 + clhs358;
lhs(2,9)=-clhs17*clhs367 - clhs17*clhs368 - clhs17*clhs369 + clhs362 + clhs363 + clhs365 + clhs366;
lhs(2,10)=-clhs17*clhs376 - clhs17*clhs377 - clhs17*clhs378 + clhs275 + clhs370 + clhs371 + clhs373 + clhs375;
lhs(2,11)=DN(2,2)*clhs24 - DN(2,2)*clhs26 - clhs379 + clhs380*clhs60;
lhs(2,12)=-clhs17*clhs386 - clhs17*clhs387 - clhs17*clhs388 + clhs381 + clhs382 + clhs384 + clhs385;
lhs(2,13)=-clhs17*clhs394 - clhs17*clhs395 - clhs17*clhs396 + clhs389 + clhs390 + clhs392 + clhs393;
lhs(2,14)=-clhs17*clhs403 - clhs17*clhs404 - clhs17*clhs405 + clhs308 + clhs397 + clhs398 + clhs400 + clhs402;
lhs(2,15)=DN(3,2)*clhs24 - DN(3,2)*clhs26 - clhs406 + clhs407*clhs60;
lhs(3,0)=DN(0,0)*clhs408;
lhs(3,1)=DN(0,1)*clhs408;
lhs(3,2)=DN(0,2)*clhs408;
lhs(3,3)=clhs16*(clhs209 + clhs325 + clhs8);
lhs(3,4)=DN(0,0)*clhs80 + clhs108;
lhs(3,5)=DN(0,1)*clhs80 + clhs254;
lhs(3,6)=DN(0,2)*clhs80 + clhs353;
lhs(3,7)=clhs409;
lhs(3,8)=DN(0,0)*clhs127 + clhs155;
lhs(3,9)=DN(0,1)*clhs127 + clhs287;
lhs(3,10)=DN(0,2)*clhs127 + clhs380;
lhs(3,11)=clhs410;
lhs(3,12)=DN(0,0)*clhs174 + clhs202;
lhs(3,13)=DN(0,1)*clhs174 + clhs320;
lhs(3,14)=DN(0,2)*clhs174 + clhs407;
lhs(3,15)=clhs411;
lhs(4,0)=-clhs17*clhs63 - clhs17*clhs66 - clhs17*clhs69 + clhs415 + clhs70 + clhs71 + clhs72 + clhs75;
lhs(4,1)=-clhs17*clhs222 - clhs17*clhs224 - clhs17*clhs226 + clhs227 + clhs228 + clhs229 + clhs230;
lhs(4,2)=-clhs17*clhs326 - clhs17*clhs327 - clhs17*clhs329 + clhs331 + clhs332 + clhs333 + clhs334;
lhs(4,3)=DN(0,0)*clhs412 - DN(0,0)*clhs413 + clhs107*clhs60 - clhs108;
lhs(4,4)=clhs13*clhs419 - clhs17*clhs416 - clhs17*clhs417 - clhs17*clhs418 + clhs416 + clhs417 + clhs418 + clhs421;
lhs(4,5)=-clhs17*clhs427 - clhs17*clhs428 - clhs17*clhs429 + clhs422 + clhs423 + clhs424 + clhs426;
lhs(4,6)=-clhs17*clhs434 - clhs17*clhs435 - clhs17*clhs436 + clhs430 + clhs431 + clhs432 + clhs433;
lhs(4,7)=DN(1,0)*clhs437;
lhs(4,8)=-clhs17*clhs441 - clhs17*clhs442 - clhs17*clhs443 + clhs438 + clhs439 + clhs440 + clhs446 + clhs447;
lhs(4,9)=-clhs17*clhs452 - clhs17*clhs453 - clhs17*clhs454 + clhs448 + clhs449 + clhs450 + clhs451;
lhs(4,10)=-clhs17*clhs459 - clhs17*clhs460 - clhs17*clhs461 + clhs455 + clhs456 + clhs457 + clhs458;
lhs(4,11)=DN(2,0)*clhs412 - DN(2,0)*clhs413 - clhs462 + clhs463*clhs60;
lhs(4,12)=-clhs17*clhs467 - clhs17*clhs468 - clhs17*clhs469 + clhs464 + clhs465 + clhs466 + clhs472 + clhs473;
lhs(4,13)=-clhs17*clhs478 - clhs17*clhs479 - clhs17*clhs480 + clhs474 + clhs475 + clhs476 + clhs477;
lhs(4,14)=-clhs17*clhs485 - clhs17*clhs486 - clhs17*clhs487 + clhs481 + clhs482 + clhs483 + clhs484;
lhs(4,15)=DN(3,0)*clhs412 - DN(3,0)*clhs413 - clhs488 + clhs489*clhs60;
lhs(5,0)=-clhs17*clhs83 - clhs17*clhs86 - clhs17*clhs90 + clhs91 + clhs92 + clhs93 + clhs94;
lhs(5,1)=-clhs17*clhs231 - clhs17*clhs233 - clhs17*clhs236 + clhs238 + clhs239 + clhs240 + clhs241 + clhs490;
lhs(5,2)=-clhs17*clhs335 - clhs17*clhs336 - clhs17*clhs338 + clhs339 + clhs340 + clhs341 + clhs342;
lhs(5,3)=DN(0,1)*clhs412 - DN(0,1)*clhs413 + clhs253*clhs60 - clhs254;
lhs(5,4)=-clhs17*clhs422 - clhs17*clhs423 - clhs17*clhs424 + clhs426 + clhs427 + clhs428 + clhs429;
lhs(5,5)=clhs13*clhs494 - clhs17*clhs491 - clhs17*clhs492 - clhs17*clhs493 + clhs421 + clhs491 + clhs492 + clhs493;
lhs(5,6)=-clhs17*clhs500 - clhs17*clhs501 - clhs17*clhs502 + clhs495 + clhs496 + clhs497 + clhs499;
lhs(5,7)=DN(1,1)*clhs437;
lhs(5,8)=-clhs17*clhs507 - clhs17*clhs508 - clhs17*clhs509 + clhs503 + clhs504 + clhs505 + clhs506;
lhs(5,9)=-clhs17*clhs515 - clhs17*clhs516 - clhs17*clhs517 + clhs510 + clhs511 + clhs512 + clhs514 + clhs518;
lhs(5,10)=-clhs17*clhs523 - clhs17*clhs524 - clhs17*clhs525 + clhs519 + clhs520 + clhs521 + clhs522;
lhs(5,11)=DN(2,1)*clhs412 - DN(2,1)*clhs413 - clhs526 + clhs527*clhs60;
lhs(5,12)=-clhs17*clhs532 - clhs17*clhs533 - clhs17*clhs534 + clhs528 + clhs529 + clhs530 + clhs531;
lhs(5,13)=-clhs17*clhs540 - clhs17*clhs541 - clhs17*clhs542 + clhs535 + clhs536 + clhs537 + clhs539 + clhs543;
lhs(5,14)=-clhs17*clhs548 - clhs17*clhs549 - clhs17*clhs550 + clhs544 + clhs545 + clhs546 + clhs547;
lhs(5,15)=DN(3,1)*clhs412 - DN(3,1)*clhs413 - clhs551 + clhs552*clhs60;
lhs(6,0)=-clhs102*clhs17 + clhs103 + clhs104 + clhs105 + clhs106 - clhs17*clhs96 - clhs17*clhs99;
lhs(6,1)=-clhs17*clhs243 - clhs17*clhs245 - clhs17*clhs248 + clhs249 + clhs250 + clhs251 + clhs252;
lhs(6,2)=-clhs17*clhs343 - clhs17*clhs344 - clhs17*clhs346 + clhs348 + clhs349 + clhs350 + clhs351 + clhs490;
lhs(6,3)=DN(0,2)*clhs412 - DN(0,2)*clhs413 + clhs352*clhs60 - clhs353;
lhs(6,4)=-clhs17*clhs430 - clhs17*clhs431 - clhs17*clhs432 + clhs433 + clhs434 + clhs435 + clhs436;
lhs(6,5)=-clhs17*clhs495 - clhs17*clhs496 - clhs17*clhs497 + clhs499 + clhs500 + clhs501 + clhs502;
lhs(6,6)=clhs13*clhs556 - clhs17*clhs553 - clhs17*clhs554 - clhs17*clhs555 + clhs421 + clhs553 + clhs554 + clhs555;
lhs(6,7)=DN(1,2)*clhs437;
lhs(6,8)=-clhs17*clhs562 - clhs17*clhs563 - clhs17*clhs564 + clhs557 + clhs558 + clhs559 + clhs561;
lhs(6,9)=-clhs17*clhs569 - clhs17*clhs570 - clhs17*clhs571 + clhs565 + clhs566 + clhs567 + clhs568;
lhs(6,10)=-clhs17*clhs577 - clhs17*clhs578 - clhs17*clhs579 + clhs518 + clhs572 + clhs573 + clhs574 + clhs576;
lhs(6,11)=DN(2,2)*clhs412 - DN(2,2)*clhs413 - clhs580 + clhs581*clhs60;
lhs(6,12)=-clhs17*clhs586 - clhs17*clhs587 - clhs17*clhs588 + clhs582 + clhs583 + clhs584 + clhs585;
lhs(6,13)=-clhs17*clhs593 - clhs17*clhs594 - clhs17*clhs595 + clhs589 + clhs590 + clhs591 + clhs592;
lhs(6,14)=-clhs17*clhs601 - clhs17*clhs602 - clhs17*clhs603 + clhs543 + clhs596 + clhs597 + clhs598 + clhs600;
lhs(6,15)=DN(3,2)*clhs412 - DN(3,2)*clhs413 + clhs60*clhs605 - clhs604;
lhs(7,0)=DN(1,0)*clhs27 + clhs107;
lhs(7,1)=DN(1,1)*clhs27 + clhs253;
lhs(7,2)=DN(1,2)*clhs27 + clhs352;
lhs(7,3)=clhs409;
lhs(7,4)=DN(1,0)*clhs606;
lhs(7,5)=DN(1,1)*clhs606;
lhs(7,6)=DN(1,2)*clhs606;
lhs(7,7)=clhs16*(clhs419 + clhs494 + clhs556);
lhs(7,8)=DN(1,0)*clhs127 + clhs463;
lhs(7,9)=DN(1,1)*clhs127 + clhs527;
lhs(7,10)=DN(1,2)*clhs127 + clhs581;
lhs(7,11)=clhs607;
lhs(7,12)=DN(1,0)*clhs174 + clhs489;
lhs(7,13)=DN(1,1)*clhs174 + clhs552;
lhs(7,14)=DN(1,2)*clhs174 + clhs605;
lhs(7,15)=clhs608;
lhs(8,0)=-clhs110*clhs17 - clhs113*clhs17 - clhs116*clhs17 + clhs117 + clhs118 + clhs119 + clhs122 + clhs610;
lhs(8,1)=-clhs17*clhs255 - clhs17*clhs257 - clhs17*clhs259 + clhs260 + clhs261 + clhs262 + clhs263;
lhs(8,2)=-clhs17*clhs354 - clhs17*clhs355 - clhs17*clhs357 + clhs358 + clhs359 + clhs360 + clhs361;
lhs(8,3)=-DN(0,0)*clhs611 + DN(0,0)*clhs612 + clhs154*clhs60 - clhs155;
lhs(8,4)=-clhs17*clhs438 - clhs17*clhs439 - clhs17*clhs440 + clhs441 + clhs442 + clhs443 + clhs446 + clhs613;
lhs(8,5)=-clhs17*clhs503 - clhs17*clhs504 - clhs17*clhs505 + clhs506 + clhs507 + clhs508 + clhs509;
lhs(8,6)=-clhs17*clhs557 - clhs17*clhs558 - clhs17*clhs559 + clhs561 + clhs562 + clhs563 + clhs564;
lhs(8,7)=-DN(1,0)*clhs611 + DN(1,0)*clhs612 + clhs462*clhs60 - clhs463;
lhs(8,8)=clhs13*clhs617 - clhs17*clhs614 - clhs17*clhs615 - clhs17*clhs616 + clhs614 + clhs615 + clhs616 + clhs619;
lhs(8,9)=-clhs17*clhs625 - clhs17*clhs626 - clhs17*clhs627 + clhs620 + clhs621 + clhs622 + clhs624;
lhs(8,10)=-clhs17*clhs632 - clhs17*clhs633 - clhs17*clhs634 + clhs628 + clhs629 + clhs630 + clhs631;
lhs(8,11)=DN(2,0)*clhs635;
lhs(8,12)=-clhs17*clhs639 - clhs17*clhs640 - clhs17*clhs641 + clhs636 + clhs637 + clhs638 + clhs644 + clhs645;
lhs(8,13)=-clhs17*clhs650 - clhs17*clhs651 - clhs17*clhs652 + clhs646 + clhs647 + clhs648 + clhs649;
lhs(8,14)=-clhs17*clhs657 - clhs17*clhs658 - clhs17*clhs659 + clhs653 + clhs654 + clhs655 + clhs656;
lhs(8,15)=-DN(3,0)*clhs611 + DN(3,0)*clhs612 + clhs60*clhs661 - clhs660;
lhs(9,0)=-clhs130*clhs17 - clhs133*clhs17 - clhs137*clhs17 + clhs138 + clhs139 + clhs140 + clhs141;
lhs(9,1)=-clhs17*clhs264 - clhs17*clhs266 - clhs17*clhs269 + clhs271 + clhs272 + clhs273 + clhs274 + clhs662;
lhs(9,2)=-clhs17*clhs362 - clhs17*clhs363 - clhs17*clhs365 + clhs366 + clhs367 + clhs368 + clhs369;
lhs(9,3)=-DN(0,1)*clhs611 + DN(0,1)*clhs612 + clhs286*clhs60 - clhs287;
lhs(9,4)=-clhs17*clhs448 - clhs17*clhs449 - clhs17*clhs450 + clhs451 + clhs452 + clhs453 + clhs454;
lhs(9,5)=-clhs17*clhs510 - clhs17*clhs511 - clhs17*clhs512 + clhs514 + clhs515 + clhs516 + clhs517 + clhs663;
lhs(9,6)=-clhs17*clhs565 - clhs17*clhs566 - clhs17*clhs567 + clhs568 + clhs569 + clhs570 + clhs571;
lhs(9,7)=-DN(1,1)*clhs611 + DN(1,1)*clhs612 + clhs526*clhs60 - clhs527;
lhs(9,8)=-clhs17*clhs620 - clhs17*clhs621 - clhs17*clhs622 + clhs624 + clhs625 + clhs626 + clhs627;
lhs(9,9)=clhs13*clhs667 - clhs17*clhs664 - clhs17*clhs665 - clhs17*clhs666 + clhs619 + clhs664 + clhs665 + clhs666;
lhs(9,10)=-clhs17*clhs673 - clhs17*clhs674 - clhs17*clhs675 + clhs668 + clhs669 + clhs670 + clhs672;
lhs(9,11)=DN(2,1)*clhs635;
lhs(9,12)=-clhs17*clhs680 - clhs17*clhs681 - clhs17*clhs682 + clhs676 + clhs677 + clhs678 + clhs679;
lhs(9,13)=-clhs17*clhs688 - clhs17*clhs689 - clhs17*clhs690 + clhs683 + clhs684 + clhs685 + clhs687 + clhs691;
lhs(9,14)=-clhs17*clhs696 - clhs17*clhs697 - clhs17*clhs698 + clhs692 + clhs693 + clhs694 + clhs695;
lhs(9,15)=-DN(3,1)*clhs611 + DN(3,1)*clhs612 + clhs60*clhs700 - clhs699;
lhs(10,0)=-clhs143*clhs17 - clhs146*clhs17 - clhs149*clhs17 + clhs150 + clhs151 + clhs152 + clhs153;
lhs(10,1)=-clhs17*clhs276 - clhs17*clhs278 - clhs17*clhs281 + clhs282 + clhs283 + clhs284 + clhs285;
lhs(10,2)=-clhs17*clhs370 - clhs17*clhs371 - clhs17*clhs373 + clhs375 + clhs376 + clhs377 + clhs378 + clhs662;
lhs(10,3)=-DN(0,2)*clhs611 + DN(0,2)*clhs612 + clhs379*clhs60 - clhs380;
lhs(10,4)=-clhs17*clhs455 - clhs17*clhs456 - clhs17*clhs457 + clhs458 + clhs459 + clhs460 + clhs461;
lhs(10,5)=-clhs17*clhs519 - clhs17*clhs520 - clhs17*clhs521 + clhs522 + clhs523 + clhs524 + clhs525;
lhs(10,6)=-clhs17*clhs572 - clhs17*clhs573 - clhs17*clhs574 + clhs576 + clhs577 + clhs578 + clhs579 + clhs663;
lhs(10,7)=-DN(1,2)*clhs611 + DN(1,2)*clhs612 + clhs580*clhs60 - clhs581;
lhs(10,8)=-clhs17*clhs628 - clhs17*clhs629 - clhs17*clhs630 + clhs631 + clhs632 + clhs633 + clhs634;
lhs(10,9)=-clhs17*clhs668 - clhs17*clhs669 - clhs17*clhs670 + clhs672 + clhs673 + clhs674 + clhs675;
lhs(10,10)=clhs13*clhs704 - clhs17*clhs701 - clhs17*clhs702 - clhs17*clhs703 + clhs619 + clhs701 + clhs702 + clhs703;
lhs(10,11)=DN(2,2)*clhs635;
lhs(10,12)=-clhs17*clhs710 - clhs17*clhs711 - clhs17*clhs712 + clhs705 + clhs706 + clhs707 + clhs709;
lhs(10,13)=-clhs17*clhs717 - clhs17*clhs718 - clhs17*clhs719 + clhs713 + clhs714 + clhs715 + clhs716;
lhs(10,14)=-clhs17*clhs725 - clhs17*clhs726 - clhs17*clhs727 + clhs691 + clhs720 + clhs721 + clhs722 + clhs724;
lhs(10,15)=-DN(3,2)*clhs611 + DN(3,2)*clhs612 + clhs60*clhs729 - clhs728;
lhs(11,0)=DN(2,0)*clhs27 + clhs154;
lhs(11,1)=DN(2,1)*clhs27 + clhs286;
lhs(11,2)=DN(2,2)*clhs27 + clhs379;
lhs(11,3)=clhs410;
lhs(11,4)=DN(2,0)*clhs80 + clhs462;
lhs(11,5)=DN(2,1)*clhs80 + clhs526;
lhs(11,6)=DN(2,2)*clhs80 + clhs580;
lhs(11,7)=clhs607;
lhs(11,8)=DN(2,0)*clhs730;
lhs(11,9)=DN(2,1)*clhs730;
lhs(11,10)=DN(2,2)*clhs730;
lhs(11,11)=clhs16*(clhs617 + clhs667 + clhs704);
lhs(11,12)=DN(2,0)*clhs174 + clhs661;
lhs(11,13)=DN(2,1)*clhs174 + clhs700;
lhs(11,14)=DN(2,2)*clhs174 + clhs729;
lhs(11,15)=clhs731;
lhs(12,0)=-clhs157*clhs17 - clhs160*clhs17 - clhs163*clhs17 + clhs164 + clhs165 + clhs166 + clhs169 + clhs733;
lhs(12,1)=-clhs17*clhs288 - clhs17*clhs290 - clhs17*clhs292 + clhs293 + clhs294 + clhs295 + clhs296;
lhs(12,2)=-clhs17*clhs381 - clhs17*clhs382 - clhs17*clhs384 + clhs385 + clhs386 + clhs387 + clhs388;
lhs(12,3)=-DN(0,0)*clhs734 + DN(0,0)*clhs735 + clhs201*clhs60 - clhs202;
lhs(12,4)=-clhs17*clhs464 - clhs17*clhs465 - clhs17*clhs466 + clhs467 + clhs468 + clhs469 + clhs472 + clhs736;
lhs(12,5)=-clhs17*clhs528 - clhs17*clhs529 - clhs17*clhs530 + clhs531 + clhs532 + clhs533 + clhs534;
lhs(12,6)=-clhs17*clhs582 - clhs17*clhs583 - clhs17*clhs584 + clhs585 + clhs586 + clhs587 + clhs588;
lhs(12,7)=-DN(1,0)*clhs734 + DN(1,0)*clhs735 + clhs488*clhs60 - clhs489;
lhs(12,8)=-clhs17*clhs636 - clhs17*clhs637 - clhs17*clhs638 + clhs639 + clhs640 + clhs641 + clhs644 + clhs737;
lhs(12,9)=-clhs17*clhs676 - clhs17*clhs677 - clhs17*clhs678 + clhs679 + clhs680 + clhs681 + clhs682;
lhs(12,10)=-clhs17*clhs705 - clhs17*clhs706 - clhs17*clhs707 + clhs709 + clhs710 + clhs711 + clhs712;
lhs(12,11)=-DN(2,0)*clhs734 + DN(2,0)*clhs735 + clhs60*clhs660 - clhs661;
lhs(12,12)=clhs13*clhs741 - clhs17*clhs738 - clhs17*clhs739 - clhs17*clhs740 + clhs738 + clhs739 + clhs740 + clhs743;
lhs(12,13)=-clhs17*clhs749 - clhs17*clhs750 - clhs17*clhs751 + clhs744 + clhs745 + clhs746 + clhs748;
lhs(12,14)=-clhs17*clhs756 - clhs17*clhs757 - clhs17*clhs758 + clhs752 + clhs753 + clhs754 + clhs755;
lhs(12,15)=DN(3,0)*clhs759;
lhs(13,0)=-clhs17*clhs177 - clhs17*clhs180 - clhs17*clhs184 + clhs185 + clhs186 + clhs187 + clhs188;
lhs(13,1)=-clhs17*clhs297 - clhs17*clhs299 - clhs17*clhs302 + clhs304 + clhs305 + clhs306 + clhs307 + clhs760;
lhs(13,2)=-clhs17*clhs389 - clhs17*clhs390 - clhs17*clhs392 + clhs393 + clhs394 + clhs395 + clhs396;
lhs(13,3)=-DN(0,1)*clhs734 + DN(0,1)*clhs735 + clhs319*clhs60 - clhs320;
lhs(13,4)=-clhs17*clhs474 - clhs17*clhs475 - clhs17*clhs476 + clhs477 + clhs478 + clhs479 + clhs480;
lhs(13,5)=-clhs17*clhs535 - clhs17*clhs536 - clhs17*clhs537 + clhs539 + clhs540 + clhs541 + clhs542 + clhs761;
lhs(13,6)=-clhs17*clhs589 - clhs17*clhs590 - clhs17*clhs591 + clhs592 + clhs593 + clhs594 + clhs595;
lhs(13,7)=-DN(1,1)*clhs734 + DN(1,1)*clhs735 + clhs551*clhs60 - clhs552;
lhs(13,8)=-clhs17*clhs646 - clhs17*clhs647 - clhs17*clhs648 + clhs649 + clhs650 + clhs651 + clhs652;
lhs(13,9)=-clhs17*clhs683 - clhs17*clhs684 - clhs17*clhs685 + clhs687 + clhs688 + clhs689 + clhs690 + clhs762;
lhs(13,10)=-clhs17*clhs713 - clhs17*clhs714 - clhs17*clhs715 + clhs716 + clhs717 + clhs718 + clhs719;
lhs(13,11)=-DN(2,1)*clhs734 + DN(2,1)*clhs735 + clhs60*clhs699 - clhs700;
lhs(13,12)=-clhs17*clhs744 - clhs17*clhs745 - clhs17*clhs746 + clhs748 + clhs749 + clhs750 + clhs751;
lhs(13,13)=clhs13*clhs766 - clhs17*clhs763 - clhs17*clhs764 - clhs17*clhs765 + clhs743 + clhs763 + clhs764 + clhs765;
lhs(13,14)=-clhs17*clhs771 - clhs17*clhs772 - clhs17*clhs773 + clhs767 + clhs768 + clhs769 + clhs770;
lhs(13,15)=DN(3,1)*clhs759;
lhs(14,0)=-clhs17*clhs190 - clhs17*clhs193 - clhs17*clhs196 + clhs197 + clhs198 + clhs199 + clhs200;
lhs(14,1)=-clhs17*clhs309 - clhs17*clhs311 - clhs17*clhs314 + clhs315 + clhs316 + clhs317 + clhs318;
lhs(14,2)=-clhs17*clhs397 - clhs17*clhs398 - clhs17*clhs400 + clhs402 + clhs403 + clhs404 + clhs405 + clhs760;
lhs(14,3)=-DN(0,2)*clhs734 + DN(0,2)*clhs735 + clhs406*clhs60 - clhs407;
lhs(14,4)=-clhs17*clhs481 - clhs17*clhs482 - clhs17*clhs483 + clhs484 + clhs485 + clhs486 + clhs487;
lhs(14,5)=-clhs17*clhs544 - clhs17*clhs545 - clhs17*clhs546 + clhs547 + clhs548 + clhs549 + clhs550;
lhs(14,6)=-clhs17*clhs596 - clhs17*clhs597 - clhs17*clhs598 + clhs600 + clhs601 + clhs602 + clhs603 + clhs761;
lhs(14,7)=-DN(1,2)*clhs734 + DN(1,2)*clhs735 + clhs60*clhs604 - clhs605;
lhs(14,8)=-clhs17*clhs653 - clhs17*clhs654 - clhs17*clhs655 + clhs656 + clhs657 + clhs658 + clhs659;
lhs(14,9)=-clhs17*clhs692 - clhs17*clhs693 - clhs17*clhs694 + clhs695 + clhs696 + clhs697 + clhs698;
lhs(14,10)=-clhs17*clhs720 - clhs17*clhs721 - clhs17*clhs722 + clhs724 + clhs725 + clhs726 + clhs727 + clhs762;
lhs(14,11)=-DN(2,2)*clhs734 + DN(2,2)*clhs735 + clhs60*clhs728 - clhs729;
lhs(14,12)=-clhs17*clhs752 - clhs17*clhs753 - clhs17*clhs754 + clhs755 + clhs756 + clhs757 + clhs758;
lhs(14,13)=-clhs17*clhs767 - clhs17*clhs768 - clhs17*clhs769 + clhs770 + clhs771 + clhs772 + clhs773;
lhs(14,14)=clhs13*clhs777 - clhs17*clhs774 - clhs17*clhs775 - clhs17*clhs776 + clhs743 + clhs774 + clhs775 + clhs776;
lhs(14,15)=DN(3,2)*clhs759;
lhs(15,0)=DN(3,0)*clhs27 + clhs201;
lhs(15,1)=DN(3,1)*clhs27 + clhs319;
lhs(15,2)=DN(3,2)*clhs27 + clhs406;
lhs(15,3)=clhs411;
lhs(15,4)=DN(3,0)*clhs80 + clhs488;
lhs(15,5)=DN(3,1)*clhs80 + clhs551;
lhs(15,6)=DN(3,2)*clhs80 + clhs604;
lhs(15,7)=clhs608;
lhs(15,8)=DN(3,0)*clhs127 + clhs660;
lhs(15,9)=DN(3,1)*clhs127 + clhs699;
lhs(15,10)=DN(3,2)*clhs127 + clhs728;
lhs(15,11)=clhs731;
lhs(15,12)=DN(3,0)*clhs778;
lhs(15,13)=DN(3,1)*clhs778;
lhs(15,14)=DN(3,2)*clhs778;
lhs(15,15)=clhs16*(clhs741 + clhs766 + clhs777);


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
    const Matrix &C = rData.C;

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
const double crhs3 = bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0);
const double crhs4 = bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0);
const double crhs5 = bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0);
const double crhs6 = rho*(N[0]*crhs3 + N[1]*crhs4 + N[2]*crhs5);
const double crhs7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs10 = rho*(crhs7*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs11 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs12 = crhs11 + crhs7 - volume_error_ratio;
const double crhs13 = rho*stab_c2*sqrt(pow(crhs8, 2) + pow(crhs9, 2));
const double crhs14 = crhs12*(crhs13*h/stab_c1 + mu);
const double crhs15 = 1.0/(K_darcy + crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 = crhs15*rho;
const double crhs17 = crhs16*(DN(0,0)*crhs3 + DN(1,0)*crhs4 + DN(2,0)*crhs5);
const double crhs18 = C(1,2)*DN(0,1);
const double crhs19 = bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1);
const double crhs20 = bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1);
const double crhs21 = bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1);
const double crhs22 = crhs16*(DN(0,1)*crhs19 + DN(1,1)*crhs20 + DN(2,1)*crhs21);
const double crhs23 = C(0,2)*DN(0,0);
const double crhs24 = crhs16*(DN(0,0)*crhs19 + DN(0,1)*crhs3 + DN(1,0)*crhs20 + DN(1,1)*crhs4 + DN(2,0)*crhs21 + DN(2,1)*crhs5);
const double crhs25 = crhs15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs10 + crhs2 + crhs6);
const double crhs26 = K_darcy*N[0];
const double crhs27 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs28 = N[0]*crhs27;
const double crhs29 = rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs30 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs31 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs32 = rho*(N[0]*crhs19 + N[1]*crhs20 + N[2]*crhs21);
const double crhs33 = rho*(crhs11*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs34 = crhs15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs30 + crhs31 + crhs32 + crhs33);
const double crhs35 = C(1,2)*DN(1,1);
const double crhs36 = C(0,2)*DN(1,0);
const double crhs37 = K_darcy*N[1];
const double crhs38 = N[1]*crhs27;
const double crhs39 = rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs40 = C(1,2)*DN(2,1);
const double crhs41 = C(0,2)*DN(2,0);
const double crhs42 = K_darcy*N[2];
const double crhs43 = N[2]*crhs27;
const double crhs44 = rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - N[0]*crhs6 + crhs17*(C(0,0)*DN(0,0) + C(0,2)*DN(0,1)) + crhs22*(C(0,1)*DN(0,0) + crhs18) + crhs24*(C(2,2)*DN(0,1) + crhs23) + crhs25*crhs26 - crhs25*crhs28 - crhs25*crhs29;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs14 - DN(0,1)*stress[1] + N[0]*crhs30 - N[0]*crhs31 - N[0]*crhs32 - N[0]*crhs33 + crhs17*(C(0,1)*DN(0,1) + crhs23) + crhs22*(C(1,1)*DN(0,1) + C(1,2)*DN(0,0)) + crhs24*(C(2,2)*DN(0,0) + crhs18) + crhs26*crhs34 - crhs28*crhs34 - crhs29*crhs34;
rhs[2]=-DN(0,0)*crhs25 - DN(0,1)*crhs34 - N[0]*crhs12;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - N[1]*crhs6 + crhs17*(C(0,0)*DN(1,0) + C(0,2)*DN(1,1)) + crhs22*(C(0,1)*DN(1,0) + crhs35) + crhs24*(C(2,2)*DN(1,1) + crhs36) + crhs25*crhs37 - crhs25*crhs38 - crhs25*crhs39;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs14 - DN(1,1)*stress[1] + N[1]*crhs30 - N[1]*crhs31 - N[1]*crhs32 - N[1]*crhs33 + crhs17*(C(0,1)*DN(1,1) + crhs36) + crhs22*(C(1,1)*DN(1,1) + C(1,2)*DN(1,0)) + crhs24*(C(2,2)*DN(1,0) + crhs35) + crhs34*crhs37 - crhs34*crhs38 - crhs34*crhs39;
rhs[5]=-DN(1,0)*crhs25 - DN(1,1)*crhs34 - N[1]*crhs12;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - N[2]*crhs6 + crhs17*(C(0,0)*DN(2,0) + C(0,2)*DN(2,1)) + crhs22*(C(0,1)*DN(2,0) + crhs40) + crhs24*(C(2,2)*DN(2,1) + crhs41) + crhs25*crhs42 - crhs25*crhs43 - crhs25*crhs44;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs14 - DN(2,1)*stress[1] + N[2]*crhs30 - N[2]*crhs31 - N[2]*crhs32 - N[2]*crhs33 + crhs17*(C(0,1)*DN(2,1) + crhs41) + crhs22*(C(1,1)*DN(2,1) + C(1,2)*DN(2,0)) + crhs24*(C(2,2)*DN(2,0) + crhs40) + crhs34*crhs42 - crhs34*crhs43 - crhs34*crhs44;
rhs[8]=-DN(2,0)*crhs25 - DN(2,1)*crhs34 - N[2]*crhs12;


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
    const Matrix &C = rData.C;

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
const double crhs3 = bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0);
const double crhs4 = bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0);
const double crhs5 = bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0);
const double crhs6 = bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0);
const double crhs7 = rho*(N[0]*crhs3 + N[1]*crhs4 + N[2]*crhs5 + N[3]*crhs6);
const double crhs8 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs9 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs11 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs12 = rho*(crhs10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs11*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*crhs9);
const double crhs13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs14 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs15 = crhs13 + crhs14 + crhs8 - volume_error_ratio;
const double crhs16 = rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs11, 2) + pow(crhs9, 2));
const double crhs17 = crhs15*(crhs16*h/stab_c1 + mu);
const double crhs18 = 1.0/(K_darcy + crhs16/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs19 = crhs18*rho;
const double crhs20 = crhs19*(DN(0,0)*crhs3 + DN(1,0)*crhs4 + DN(2,0)*crhs5 + DN(3,0)*crhs6);
const double crhs21 = C(1,3)*DN(0,1);
const double crhs22 = bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1);
const double crhs23 = bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1);
const double crhs24 = bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1);
const double crhs25 = bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1);
const double crhs26 = crhs19*(DN(0,1)*crhs22 + DN(1,1)*crhs23 + DN(2,1)*crhs24 + DN(3,1)*crhs25);
const double crhs27 = C(2,5)*DN(0,2);
const double crhs28 = bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2);
const double crhs29 = bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2);
const double crhs30 = bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2);
const double crhs31 = bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2);
const double crhs32 = crhs19*(DN(0,2)*crhs28 + DN(1,2)*crhs29 + DN(2,2)*crhs30 + DN(3,2)*crhs31);
const double crhs33 = C(0,3)*DN(0,0);
const double crhs34 = crhs19*(DN(0,0)*crhs22 + DN(0,1)*crhs3 + DN(1,0)*crhs23 + DN(1,1)*crhs4 + DN(2,0)*crhs24 + DN(2,1)*crhs5 + DN(3,0)*crhs25 + DN(3,1)*crhs6);
const double crhs35 = C(3,4)*DN(0,1);
const double crhs36 = C(4,5)*DN(0,2);
const double crhs37 = crhs19*(DN(0,1)*crhs28 + DN(0,2)*crhs22 + DN(1,1)*crhs29 + DN(1,2)*crhs23 + DN(2,1)*crhs30 + DN(2,2)*crhs24 + DN(3,1)*crhs31 + DN(3,2)*crhs25);
const double crhs38 = C(0,5)*DN(0,0);
const double crhs39 = crhs19*(DN(0,0)*crhs28 + DN(0,2)*crhs3 + DN(1,0)*crhs29 + DN(1,2)*crhs4 + DN(2,0)*crhs30 + DN(2,2)*crhs5 + DN(3,0)*crhs31 + DN(3,2)*crhs6);
const double crhs40 = crhs18*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs12 + crhs2 + crhs7);
const double crhs41 = K_darcy*N[0];
const double crhs42 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs43 = N[0]*crhs42;
const double crhs44 = rho*(DN(0,0)*crhs9 + DN(0,1)*crhs10 + DN(0,2)*crhs11);
const double crhs45 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs46 = K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs47 = rho*(N[0]*crhs22 + N[1]*crhs23 + N[2]*crhs24 + N[3]*crhs25);
const double crhs48 = rho*(crhs10*crhs13 + crhs11*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs9*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)));
const double crhs49 = C(2,4)*DN(0,2);
const double crhs50 = C(1,4)*DN(0,1);
const double crhs51 = C(3,5)*DN(0,0);
const double crhs52 = crhs18*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs45 + crhs46 + crhs47 + crhs48);
const double crhs53 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs54 = K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs55 = rho*(N[0]*crhs28 + N[1]*crhs29 + N[2]*crhs30 + N[3]*crhs31);
const double crhs56 = rho*(crhs10*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs11*crhs14 + crhs9*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)));
const double crhs57 = crhs18*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs53 + crhs54 + crhs55 + crhs56);
const double crhs58 = C(1,3)*DN(1,1);
const double crhs59 = C(2,5)*DN(1,2);
const double crhs60 = C(0,3)*DN(1,0);
const double crhs61 = C(3,4)*DN(1,1);
const double crhs62 = C(4,5)*DN(1,2);
const double crhs63 = C(0,5)*DN(1,0);
const double crhs64 = K_darcy*N[1];
const double crhs65 = N[1]*crhs42;
const double crhs66 = rho*(DN(1,0)*crhs9 + DN(1,1)*crhs10 + DN(1,2)*crhs11);
const double crhs67 = C(2,4)*DN(1,2);
const double crhs68 = C(1,4)*DN(1,1);
const double crhs69 = C(3,5)*DN(1,0);
const double crhs70 = C(1,3)*DN(2,1);
const double crhs71 = C(2,5)*DN(2,2);
const double crhs72 = C(0,3)*DN(2,0);
const double crhs73 = C(3,4)*DN(2,1);
const double crhs74 = C(4,5)*DN(2,2);
const double crhs75 = C(0,5)*DN(2,0);
const double crhs76 = K_darcy*N[2];
const double crhs77 = N[2]*crhs42;
const double crhs78 = rho*(DN(2,0)*crhs9 + DN(2,1)*crhs10 + DN(2,2)*crhs11);
const double crhs79 = C(2,4)*DN(2,2);
const double crhs80 = C(1,4)*DN(2,1);
const double crhs81 = C(3,5)*DN(2,0);
const double crhs82 = C(1,3)*DN(3,1);
const double crhs83 = C(2,5)*DN(3,2);
const double crhs84 = C(0,3)*DN(3,0);
const double crhs85 = C(3,4)*DN(3,1);
const double crhs86 = C(4,5)*DN(3,2);
const double crhs87 = C(0,5)*DN(3,0);
const double crhs88 = K_darcy*N[3];
const double crhs89 = N[3]*crhs42;
const double crhs90 = rho*(DN(3,0)*crhs9 + DN(3,1)*crhs10 + DN(3,2)*crhs11);
const double crhs91 = C(2,4)*DN(3,2);
const double crhs92 = C(1,4)*DN(3,1);
const double crhs93 = C(3,5)*DN(3,0);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs17 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs12 - N[0]*crhs2 - N[0]*crhs7 + crhs20*(C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2)) + crhs26*(C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crhs21) + crhs32*(C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crhs27) + crhs34*(C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crhs33) + crhs37*(C(0,4)*DN(0,0) + crhs35 + crhs36) + crhs39*(C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crhs38) + crhs40*crhs41 - crhs40*crhs43 - crhs40*crhs44;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs17 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs45 - N[0]*crhs46 - N[0]*crhs47 - N[0]*crhs48 + crhs20*(C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crhs33) + crhs26*(C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2)) + crhs32*(C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crhs49) + crhs34*(C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crhs21) + crhs37*(C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crhs50) + crhs39*(C(1,5)*DN(0,1) + crhs36 + crhs51) + crhs41*crhs52 - crhs43*crhs52 - crhs44*crhs52;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs17 - DN(0,2)*stress[2] + N[0]*crhs53 - N[0]*crhs54 - N[0]*crhs55 - N[0]*crhs56 + crhs20*(C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crhs38) + crhs26*(C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crhs50) + crhs32*(C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0)) + crhs34*(C(2,3)*DN(0,2) + crhs35 + crhs51) + crhs37*(C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crhs49) + crhs39*(C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crhs27) + crhs41*crhs57 - crhs43*crhs57 - crhs44*crhs57;
rhs[3]=-DN(0,0)*crhs40 - DN(0,1)*crhs52 - DN(0,2)*crhs57 - N[0]*crhs15;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs17 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs12 - N[1]*crhs2 - N[1]*crhs7 + crhs20*(C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2)) + crhs26*(C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crhs58) + crhs32*(C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crhs59) + crhs34*(C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crhs60) + crhs37*(C(0,4)*DN(1,0) + crhs61 + crhs62) + crhs39*(C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crhs63) + crhs40*crhs64 - crhs40*crhs65 - crhs40*crhs66;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs17 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs45 - N[1]*crhs46 - N[1]*crhs47 - N[1]*crhs48 + crhs20*(C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crhs60) + crhs26*(C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2)) + crhs32*(C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crhs67) + crhs34*(C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crhs58) + crhs37*(C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crhs68) + crhs39*(C(1,5)*DN(1,1) + crhs62 + crhs69) + crhs52*crhs64 - crhs52*crhs65 - crhs52*crhs66;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs17 - DN(1,2)*stress[2] + N[1]*crhs53 - N[1]*crhs54 - N[1]*crhs55 - N[1]*crhs56 + crhs20*(C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crhs63) + crhs26*(C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crhs68) + crhs32*(C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0)) + crhs34*(C(2,3)*DN(1,2) + crhs61 + crhs69) + crhs37*(C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crhs67) + crhs39*(C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crhs59) + crhs57*crhs64 - crhs57*crhs65 - crhs57*crhs66;
rhs[7]=-DN(1,0)*crhs40 - DN(1,1)*crhs52 - DN(1,2)*crhs57 - N[1]*crhs15;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs17 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs12 - N[2]*crhs2 - N[2]*crhs7 + crhs20*(C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2)) + crhs26*(C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crhs70) + crhs32*(C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crhs71) + crhs34*(C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crhs72) + crhs37*(C(0,4)*DN(2,0) + crhs73 + crhs74) + crhs39*(C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crhs75) + crhs40*crhs76 - crhs40*crhs77 - crhs40*crhs78;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs17 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs45 - N[2]*crhs46 - N[2]*crhs47 - N[2]*crhs48 + crhs20*(C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crhs72) + crhs26*(C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2)) + crhs32*(C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crhs79) + crhs34*(C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crhs70) + crhs37*(C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crhs80) + crhs39*(C(1,5)*DN(2,1) + crhs74 + crhs81) + crhs52*crhs76 - crhs52*crhs77 - crhs52*crhs78;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs17 - DN(2,2)*stress[2] + N[2]*crhs53 - N[2]*crhs54 - N[2]*crhs55 - N[2]*crhs56 + crhs20*(C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crhs75) + crhs26*(C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crhs80) + crhs32*(C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0)) + crhs34*(C(2,3)*DN(2,2) + crhs73 + crhs81) + crhs37*(C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crhs79) + crhs39*(C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crhs71) + crhs57*crhs76 - crhs57*crhs77 - crhs57*crhs78;
rhs[11]=-DN(2,0)*crhs40 - DN(2,1)*crhs52 - DN(2,2)*crhs57 - N[2]*crhs15;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs17 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs12 - N[3]*crhs2 - N[3]*crhs7 + crhs20*(C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2)) + crhs26*(C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crhs82) + crhs32*(C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crhs83) + crhs34*(C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crhs84) + crhs37*(C(0,4)*DN(3,0) + crhs85 + crhs86) + crhs39*(C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crhs87) + crhs40*crhs88 - crhs40*crhs89 - crhs40*crhs90;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs17 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs45 - N[3]*crhs46 - N[3]*crhs47 - N[3]*crhs48 + crhs20*(C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crhs84) + crhs26*(C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2)) + crhs32*(C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crhs91) + crhs34*(C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crhs82) + crhs37*(C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crhs92) + crhs39*(C(1,5)*DN(3,1) + crhs86 + crhs93) + crhs52*crhs88 - crhs52*crhs89 - crhs52*crhs90;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs17 - DN(3,2)*stress[2] + N[3]*crhs53 - N[3]*crhs54 - N[3]*crhs55 - N[3]*crhs56 + crhs20*(C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crhs87) + crhs26*(C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crhs92) + crhs32*(C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0)) + crhs34*(C(2,3)*DN(3,2) + crhs85 + crhs93) + crhs37*(C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crhs91) + crhs39*(C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crhs83) + crhs57*crhs88 - crhs57*crhs89 - crhs57*crhs90;
rhs[15]=-DN(3,0)*crhs40 - DN(3,1)*crhs52 - DN(3,2)*crhs57 - N[3]*crhs15;


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
const double cV4 = cV2*rho;
const double cV5 = cV4*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double cV6 = -N[0]*cV3 + N[0]*cV5 + N[0] + cV4*(DN(0,0)*cV0 + DN(0,1)*cV1);
const double cV7 = -N[1]*cV3 + N[1]*cV5 + N[1] + cV4*(DN(1,0)*cV0 + DN(1,1)*cV1);
const double cV8 = -N[2]*cV3 + N[2]*cV5 + N[2] + cV4*(DN(2,0)*cV0 + DN(2,1)*cV1);
V(0,0)=DNenr(0,0)*cV6;
V(0,1)=DNenr(1,0)*cV6;
V(0,2)=DNenr(2,0)*cV6;
V(1,0)=DNenr(0,1)*cV6;
V(1,1)=DNenr(1,1)*cV6;
V(1,2)=DNenr(2,1)*cV6;
V(2,0)=cV2*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
V(2,1)=cV2*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
V(2,2)=cV2*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
V(3,0)=DNenr(0,0)*cV7;
V(3,1)=DNenr(1,0)*cV7;
V(3,2)=DNenr(2,0)*cV7;
V(4,0)=DNenr(0,1)*cV7;
V(4,1)=DNenr(1,1)*cV7;
V(4,2)=DNenr(2,1)*cV7;
V(5,0)=cV2*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
V(5,1)=cV2*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
V(5,2)=cV2*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
V(6,0)=DNenr(0,0)*cV8;
V(6,1)=DNenr(1,0)*cV8;
V(6,2)=DNenr(2,0)*cV8;
V(7,0)=DNenr(0,1)*cV8;
V(7,1)=DNenr(1,1)*cV8;
V(7,2)=DNenr(2,1)*cV8;
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
const double crhs_ee3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee5 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee6 = crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 = crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
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
const double cV5 = cV3*rho;
const double cV6 = cV5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double cV7 = -N[0]*cV4 + N[0]*cV6 + N[0] + cV5*(DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2);
const double cV8 = -N[1]*cV4 + N[1]*cV6 + N[1] + cV5*(DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2);
const double cV9 = -N[2]*cV4 + N[2]*cV6 + N[2] + cV5*(DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2);
const double cV10 = -N[3]*cV4 + N[3]*cV6 + N[3] + cV5*(DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2);
V(0,0)=DNenr(0,0)*cV7;
V(0,1)=DNenr(1,0)*cV7;
V(0,2)=DNenr(2,0)*cV7;
V(0,3)=DNenr(3,0)*cV7;
V(1,0)=DNenr(0,1)*cV7;
V(1,1)=DNenr(1,1)*cV7;
V(1,2)=DNenr(2,1)*cV7;
V(1,3)=DNenr(3,1)*cV7;
V(2,0)=DNenr(0,2)*cV7;
V(2,1)=DNenr(1,2)*cV7;
V(2,2)=DNenr(2,2)*cV7;
V(2,3)=DNenr(3,2)*cV7;
V(3,0)=cV3*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
V(3,1)=cV3*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
V(3,2)=cV3*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
V(3,3)=cV3*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
V(4,0)=DNenr(0,0)*cV8;
V(4,1)=DNenr(1,0)*cV8;
V(4,2)=DNenr(2,0)*cV8;
V(4,3)=DNenr(3,0)*cV8;
V(5,0)=DNenr(0,1)*cV8;
V(5,1)=DNenr(1,1)*cV8;
V(5,2)=DNenr(2,1)*cV8;
V(5,3)=DNenr(3,1)*cV8;
V(6,0)=DNenr(0,2)*cV8;
V(6,1)=DNenr(1,2)*cV8;
V(6,2)=DNenr(2,2)*cV8;
V(6,3)=DNenr(3,2)*cV8;
V(7,0)=cV3*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
V(7,1)=cV3*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
V(7,2)=cV3*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
V(7,3)=cV3*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
V(8,0)=DNenr(0,0)*cV9;
V(8,1)=DNenr(1,0)*cV9;
V(8,2)=DNenr(2,0)*cV9;
V(8,3)=DNenr(3,0)*cV9;
V(9,0)=DNenr(0,1)*cV9;
V(9,1)=DNenr(1,1)*cV9;
V(9,2)=DNenr(2,1)*cV9;
V(9,3)=DNenr(3,1)*cV9;
V(10,0)=DNenr(0,2)*cV9;
V(10,1)=DNenr(1,2)*cV9;
V(10,2)=DNenr(2,2)*cV9;
V(10,3)=DNenr(3,2)*cV9;
V(11,0)=cV3*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
V(11,1)=cV3*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
V(11,2)=cV3*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
V(11,3)=cV3*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
V(12,0)=DNenr(0,0)*cV10;
V(12,1)=DNenr(1,0)*cV10;
V(12,2)=DNenr(2,0)*cV10;
V(12,3)=DNenr(3,0)*cV10;
V(13,0)=DNenr(0,1)*cV10;
V(13,1)=DNenr(1,1)*cV10;
V(13,2)=DNenr(2,1)*cV10;
V(13,3)=DNenr(3,1)*cV10;
V(14,0)=DNenr(0,2)*cV10;
V(14,1)=DNenr(1,2)*cV10;
V(14,2)=DNenr(2,2)*cV10;
V(14,3)=DNenr(3,2)*cV10;
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
const double crhs_ee4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee7 = 1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee4, 2) + pow(crhs_ee5, 2) + pow(crhs_ee6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee8 = crhs_ee7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*crhs_ee4 + crhs_ee5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee9 = crhs_ee7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee1*crhs_ee5 + crhs_ee4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee10 = crhs_ee7*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee2*crhs_ee6 + crhs_ee4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
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

        // // Compute the maximum diagonal value in the enrichment stiffness matrix
        // double max_diag = 0.0;
        // for (unsigned int k = 0; k < NumNodes; ++k){
        //     if (std::abs(rKeeTot(k, k)) > max_diag){
        //         max_diag = std::abs(rKeeTot(k, k));
        //     }
        // }
        // if (max_diag == 0.0){
        //     max_diag = 1.0;
        // }
        // // "weakly" impose continuity
        // for (unsigned int i = 0; i < Dim; ++i){
        //     const double di = std::abs(rData.Distance[i]);
        //     for (unsigned int j = i + 1; j < NumNodes; ++j){
        //         const double dj = std::abs(rData.Distance[j]);
        //         // Check if the edge is cut, if it is, set the penalty constraint
        //         if (rData.Distance[i] * rData.Distance[j] < 0.0){
        //             double sum_d = di + dj;
        //             double Ni = dj / sum_d;
        //             double Nj = di / sum_d;
        //             double penalty_coeff = max_diag * 0.001; // h/BDFVector[0];
        //             rKeeTot(i, i) += penalty_coeff * Ni * Ni;
        //             rKeeTot(i, j) -= penalty_coeff * Ni * Nj;
        //             rKeeTot(j, i) -= penalty_coeff * Nj * Ni;
        //             rKeeTot(j, j) += penalty_coeff * Nj * Nj;
        //         }
        //     }
        // }

        // // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        // double det;
        // MatrixType inverse_diag(NumNodes, NumNodes);
        // MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

        // const Matrix tmp = prod(inverse_diag, rHtot);
        // noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

        // const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
        // noalias(rRightHandSideVector) -= prod(rVtot, tmp2);

        // -----------------------------------------------------------
        // "Strongly impose continuity"
        //----------------------------------------------------------
        MatrixType t_mpc_matrix;
        CalculateConnectivityMPCsMatrix(rData, t_mpc_matrix);
        KRATOS_WATCH(t_mpc_matrix)

        // // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        MatrixType t_transpose_matrix = ZeroMatrix(t_mpc_matrix.size2(), t_mpc_matrix.size1());
        for (unsigned int i = 0; i < t_mpc_matrix.size1(); ++i)
        {
            for (unsigned int j = 0; j < t_mpc_matrix.size2(); ++j)
            {
                t_transpose_matrix(j, i) = t_mpc_matrix(i, j);
            }
        }

        // ( Tt kee T)^-1

        double det;
        MatrixType inverse_diag;
        const Matrix ttxkenr = prod(t_transpose_matrix, rKeeTot);
        const Matrix ttxkenrxt = prod(ttxkenr, t_mpc_matrix);
        MathUtils<double>::InvertMatrix(ttxkenrxt, inverse_diag, det);
        // V T ( Tt kee T)^-1 Tt H
        const Matrix msc = prod(rVtot, t_mpc_matrix);
        const Matrix msc1 = prod(msc, inverse_diag);
        const Matrix msc2 = prod(msc1, t_transpose_matrix);

        // Adding erichment contribution to LHS
        noalias(rLeftHandSideMatrix) -= prod(msc2, rHtot);
        // V T ( Tt kee T)^-1 Tt Fenr
        // Adding erichment contribution to RHS

        const Vector msc3 = prod(msc2, rRHSeeTot);

        noalias(rRightHandSideVector) -= msc3;
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
    const double mpc_weight = 1.0;

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
