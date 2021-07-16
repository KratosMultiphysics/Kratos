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
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

#define PI 3.14159265358979

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
    ProcessInfo &rCurrentProcessInfo)
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

        //const double beta_in = 1.0e2;
        //const double beta_out = 1.0e2;
        //const double beta_contact = 1.0e-3;
        const double zeta = 5.0e-1;//1.0;//0.7;//
        const double surface_tension_coefficient = 0.072;//0.0426;//0.0311;//0.072;// //0.1; //0.0322; // //Surface tension coefficient, TODO: get from properties

        const double theta_advancing = 149.0*PI/180.0;
        const double theta_receding = 115.0*PI/180.0;
        //const double theta_static = 117*PI/180.0;
        // const double contact_line_coefficient = surface_tension_coefficient*std::cos(theta_static);
        ////const double contact_line_coefficient = -0.4539905*surface_tension_coefficient;///* 0.5299192642332 */-0.25881904510252076*surface_tension_coefficient;//0.779337965*surface_tension_coefficient;//
        const double micro_length_scale = 1.0e-9;

        //this->SetValue(CONTACT_ANGLE, 0.0); // Initialize the contact angle
        this->SetValue(CONTACT_VELOCITY, 0.0); // Initialize the contact line velocity
        const unsigned int num_dim = NumNodes - 1;
        const VectorType zero_vector = ZeroVector(num_dim);
        GeometryType::Pointer p_geom = this->pGetGeometry();

        for (unsigned int i=0; i < NumNodes; ++i){

            #pragma omp critical
            {
            (*p_geom)[i].FastGetSolutionStepValue(NORMAL_VECTOR) = zero_vector;
            (*p_geom)[i].FastGetSolutionStepValue(TANGENT_VECTOR) = zero_vector;
            (*p_geom)[i].FastGetSolutionStepValue(CONTACT_VECTOR) = zero_vector;

            /* if ((*p_geom)[i].GetValue(IS_STRUCTURE) == 1.0){
                (*p_geom)[i].Fix(DISTANCE);
            } */
            }
        }

        /* ContactSurfaceDissipation(  // DEPRECATED
            data,                      // USING NAVIER_STOKES_WALL_CONDITION INSTEAD
            beta_in,
            beta_out,
            beta_contact,
            rLeftHandSideMatrix,
            rRightHandSideVector); */

        if (data.IsCut()){
            Matrix shape_functions_pos, shape_functions_neg;
            Matrix shape_functions_enr_pos, shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos, shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos, shape_derivatives_enr_neg;
            Matrix int_shape_function_neg;
            Matrix int_shape_function_enr_neg;
            Matrix int_shape_function_enr_pos;
            GeometryType::ShapeFunctionsGradientsType int_shape_derivatives_neg;
            Kratos::Vector int_gauss_pts_weights;
            std::vector<Vector> int_normals_neg;
            Kratos::Vector gauss_pts_curvature;
            std::vector<MatrixType> contact_shape_function_neg;                                  //std::vector for multiple contact lines
            std::vector<GeometryType::ShapeFunctionsGradientsType> contact_shape_derivatives_neg;//std::vector for multiple contact lines
            std::vector<Kratos::Vector> contact_gauss_pts_weights;                                //std::vector for multiple contact lines
            std::vector<Vector> contact_tangential_neg;                                          //std::vector for multiple contact lines
            //std::vector<Vector> contact_vector;                                                  //std::vector for multiple contact lines

            //Kratos::Vector surface_tension;

            //bool has_contact_line = false; DEPRECATED

            //KRATOS_INFO("ComputeSplitting") << "Started" << std::endl;
            
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
                int_shape_function_neg,
                int_shape_function_enr_pos,
                int_shape_function_enr_neg,
                int_shape_derivatives_neg,
                int_gauss_pts_weights,
                int_normals_neg,
                contact_shape_function_neg,
                contact_shape_derivatives_neg,
                contact_gauss_pts_weights,
                contact_tangential_neg);//,
                //has_contact_line);

                //KRATOS_INFO("Size of contact_gauss_pts_weights") << contact_gauss_pts_weights.size() << std::endl;
                //KRATOS_INFO("ComputeSplitting") << "Ended" << std::endl;

            if (data.NumberOfDivisions == 1){
                // Cases exist when the element is not subdivided due to the characteristics of the provided distance
                // In this cases the element is treated as AIR or FLUID depending on the side

                KRATOS_INFO("Cut Element") << "NumberOfDivisions == 1" << std::endl;

                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                const unsigned int number_of_gauss_points = gauss_weights.size();

                //double volume = 0.0;
                //for (unsigned int gp = 0; gp < number_of_gauss_points; gp++)
                //    volume += gauss_weights[gp];

                //for (unsigned int i = 0; i < NumNodes; i++)
                //    for (unsigned int j = 0; j < Dim; j++)
                //        data.BodyForce(i,j) +=  1.0/(volume*(data.Density))*surface_tension(j);

                //array_1d<double, NumNodes> Ncenter;
                //for (unsigned int i = 0; i < NumNodes; i++){
                //    Ncenter[i] = 1.0/NumNodes;
                //}

                for (unsigned int g = 0; g < number_of_gauss_points; g++){
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
                VectorType rhs_eV_tot = ZeroVector(NumNodes * (Dim + 1));

                //KRATOS_INFO("Cut Element, Enriched Pressure") << data.Pressure_Enriched << std::endl;

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); g_pos++){
                    UpdateIntegrationPointData(
                        data,
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
                        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos, g_pos),
                        shape_derivatives_enr_pos[g_pos]);

                    //this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointLHSContributionCut(data, rLeftHandSideMatrix);
                    this->ComputeGaussPointRHSContributionCut(data, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributionsCut(data, Vtot, Htot, Kee_tot, rhs_ee_tot);

                    //PressureGradientStabilization(data, rLeftHandSideMatrix, Vtot, rRightHandSideVector);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); g_neg++){
                    UpdateIntegrationPointData(
                        data,
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
                        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg, g_neg),
                        shape_derivatives_enr_neg[g_neg]);
                    //this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
                    this->ComputeGaussPointLHSContributionCut(data, rLeftHandSideMatrix);
                    this->ComputeGaussPointRHSContributionCut(data, rRightHandSideVector);
                    ComputeGaussPointEnrichmentContributionsCut(data, Vtot, Htot, Kee_tot, rhs_ee_tot);

                    //PressureGradientStabilization(data, rLeftHandSideMatrix, Vtot, rRightHandSideVector);
                }

                CalculateCurvature(int_shape_function_neg, gauss_pts_curvature);         //Uses the saved nodal CURVATURE 
                //CalculateCurvature(int_shape_derivatives_neg, gauss_pts_curvature);    //Calculates curvature by differentiating DISTANCE_GRADIENT
                //CalculateCurvature(int_shape_function_neg, int_shape_derivatives_neg, gauss_pts_curvature); //Combination of the above functions
                //KRATOS_INFO("Curvature") << gauss_pts_curvature << std::endl;

                //SurfaceTension(surface_tension_coefficient,gauss_pts_curvature,int_gauss_pts_weights,int_normals_neg,surface_tension);

                /* SurfaceTensionMixed(
                    surface_tension_coefficient,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_derivatives_neg,
                    rRightHandSideVector); */

                /* SurfaceTension(
                    surface_tension_coefficient,
                    gauss_pts_curvature,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_normals_neg,
                    contact_gauss_pts_weights,
                    contact_shape_function_neg,
                    contact_tangential_neg,
                    has_contact_line,
                    rRightHandSideVector); */

                SurfaceTension(
                    data,
                    surface_tension_coefficient,
                    theta_advancing,
                    theta_receding,
                    zeta,
                    micro_length_scale,
                    gauss_pts_curvature,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_normals_neg,
                    contact_gauss_pts_weights,
                    contact_shape_function_neg,
                    contact_tangential_neg,
                    rLeftHandSideMatrix,
                    rRightHandSideVector);

                /* SurfaceTension(
                    data,
                    surface_tension_coefficient,
                    contact_line_coefficient,
                    zeta,
                    micro_length_scale,
                    gauss_pts_curvature,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_normals_neg,
                    contact_gauss_pts_weights,
                    contact_shape_function_neg,
                    contact_tangential_neg,
                    //has_contact_line,
                    rLeftHandSideMatrix,
                    rRightHandSideVector); */

                /* SurfaceTension(
                    surface_tension_coefficient,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_derivatives_neg,
                    int_normals_neg,
                    rRightHandSideVector); */

                /* SurfaceTension(
                    surface_tension_coefficient,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_derivatives_neg,
                    rRightHandSideVector); */

                /* SurfaceTension(
                    data,
                    surface_tension_coefficient,
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_derivatives_neg,
                    rLeftHandSideMatrix,
                    rRightHandSideVector); */

                /* PressureDiscontinuity(
                    surface_tension_coefficient,
                    gauss_pts_curvature, 
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_function_enr_pos,
                    int_shape_function_enr_neg,
                    Kee_tot,
                    rhs_ee_tot); */

                /* PressureDiscontinuity(
                    surface_tension_coefficient,
                    Kee_tot,
                    rhs_ee_tot);

                SurfaceTensionDP(
                    int_gauss_pts_weights,
                    int_shape_function_neg,
                    int_shape_function_enr_pos,
                    int_shape_function_enr_neg,
                    Vtot); */

                //KRATOS_INFO("num. neg. gauss points") << data.w_gauss_neg_side.size() << std::endl;

                //double volume_neg = 0.0;               
                //for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); g_neg++)
                //    volume_neg += data.w_gauss_neg_side[g_neg];

                //KRATOS_INFO("density") << data.NodalDensity << std::endl;

                //KRATOS_INFO("body force") << data.BodyForce << std::endl;

                //KRATOS_INFO("surface_tension force") << surface_tension << std::endl;

                //for (unsigned int i = 0; i < NumNodes; i++)
                //    for (unsigned int j = 0; j < Dim; j++)
                //        data.BodyForce(i,j) +=  1.0/(volume_neg*(data.Density))*surface_tension(j);

                /* PressureGradientStabilization(
                    data,
                    int_gauss_pts_weights,
                    int_shape_derivatives_neg,
                    Kee_tot); */

                /* PressureGradientStabilization(
                    data,
                    rLeftHandSideMatrix,
                    Vtot,
                    rRightHandSideVector); */

                PressureGradientStabilization(
                    data,
                    int_gauss_pts_weights,
                    int_shape_function_enr_pos,
                    int_shape_function_enr_neg,
                    int_shape_derivatives_neg,
                    Kee_tot,
                    rhs_ee_tot);

                /* Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
                GhostPressureGradientStabilization(
                    data,
                    gauss_weights,
                    shape_derivatives,
                    Kee_tot); */

                CondenseEnrichment(data, rLeftHandSideMatrix, rRightHandSideVector, Htot, Vtot, Kee_tot, rhs_ee_tot);
            }
        } else {
            //Get Shape function data
            Vector gauss_weights;
            Matrix shape_functions;
            ShapeFunctionDerivativesArrayType shape_derivatives;
            this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
            const unsigned int number_of_gauss_points = gauss_weights.size();
            // Iterate over integration points to evaluate local contribution

            for (unsigned int g = 0; g < number_of_gauss_points; g++){
                UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
                this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
    } else{
        KRATOS_ERROR << "TwoFluidNavierStokes is supposed to manage time integration." << std::endl;
    }

    /* noalias(rLeftHandSideMatrix) = 1.0e8 * rLeftHandSideMatrix;
    noalias(rRightHandSideVector) = 1.0e8 * rRightHandSideVector; */
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    KRATOS_CHECK_VARIABLE_KEY( DIVERGENCE );

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

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

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::Calculate(
    const Variable<Vector >& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo )
{
    noalias( rOutput ) = ZeroVector( StrainSize );
    
    if (rVariable == FLUID_STRESS) {

        // creating a new data container that goes out of scope after the function is left
        TElementData dataLocal;
        
        // transferring the velocity (among other variables)
        dataLocal.Initialize(*this, rCurrentProcessInfo);

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;

        // computing DN_DX values for the strain rate         
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        double sumOfGaussWeights = 0.0;

        for (unsigned int g = 0; g < number_of_gauss_points; g++){

            UpdateIntegrationPointData(dataLocal, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);

            const Vector gauss_point_contribution = dataLocal.ShearStress;

            noalias( rOutput ) += gauss_point_contribution * gauss_weights[g];
            sumOfGaussWeights += gauss_weights[g];
        }

        for (unsigned int i = 0; i < StrainSize; i++){
            rOutput[i] = ( 1.0 / sumOfGaussWeights ) * rOutput[i];
        }

    } else {

        Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);

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

    double dyn_tau = rData.DynamicTau;
    const double K_darcy = rData.DarcyTerm;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    /* unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &lhs = rData.lhs;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
    const double clhs1 =             C(0,2)*DN(0,0);
    const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
    const double clhs3 =             pow(DN(0,0), 2);
    const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double clhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double clhs6 =             rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
    const double clhs7 =             clhs6*h/stab_c1 + mu;
    const double clhs8 =             pow(N[0], 2);
    const double clhs9 =             bdf0*rho;
    const double clhs10 =             rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
    const double clhs11 =             K_darcy*N[0];
    const double clhs12 =             N[0]*clhs9;
    const double clhs13 =             clhs10 + clhs11 + clhs12;
    const double clhs14 =             1.0*clhs11;
    const double clhs15 =             1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double clhs16 =             clhs14*clhs15;
    const double clhs17 =             1.0*clhs10;
    const double clhs18 =             clhs15*clhs17;
    const double clhs19 =             K_darcy*clhs8 + N[0]*clhs10 - clhs13*clhs16 + clhs13*clhs18 + clhs8*clhs9;
    const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
    const double clhs21 =             C(1,2)*DN(0,1);
    const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
    const double clhs23 =             DN(0,0)*clhs7;
    const double clhs24 =             DN(0,1)*clhs23;
    const double clhs25 =             -N[0] - clhs16 + clhs18;
    const double clhs26 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
    const double clhs27 =             C(0,2)*DN(1,0);
    const double clhs28 =             C(2,2)*DN(1,1) + clhs27;
    const double clhs29 =             DN(0,0)*DN(1,0);
    const double clhs30 =             clhs29*clhs7;
    const double clhs31 =             N[1]*clhs11;
    const double clhs32 =             N[1]*clhs12;
    const double clhs33 =             rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
    const double clhs34 =             K_darcy*N[1];
    const double clhs35 =             N[1]*clhs9;
    const double clhs36 =             clhs33 + clhs34 + clhs35;
    const double clhs37 =             N[0]*clhs33 - clhs16*clhs36 + clhs18*clhs36 + clhs31 + clhs32;
    const double clhs38 =             C(0,1)*DN(1,1) + clhs27;
    const double clhs39 =             C(1,2)*DN(1,1);
    const double clhs40 =             C(2,2)*DN(1,0) + clhs39;
    const double clhs41 =             DN(1,1)*clhs23;
    const double clhs42 =             DN(0,0)*N[1];
    const double clhs43 =             DN(1,0)*N[0];
    const double clhs44 =             1.0*K_darcy*clhs15;
    const double clhs45 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
    const double clhs46 =             C(0,2)*DN(2,0);
    const double clhs47 =             C(2,2)*DN(2,1) + clhs46;
    const double clhs48 =             DN(0,0)*DN(2,0);
    const double clhs49 =             clhs48*clhs7;
    const double clhs50 =             N[2]*clhs11;
    const double clhs51 =             N[2]*clhs12;
    const double clhs52 =             rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
    const double clhs53 =             K_darcy*N[2];
    const double clhs54 =             N[2]*clhs9;
    const double clhs55 =             clhs52 + clhs53 + clhs54;
    const double clhs56 =             N[0]*clhs52 - clhs16*clhs55 + clhs18*clhs55 + clhs50 + clhs51;
    const double clhs57 =             C(0,1)*DN(2,1) + clhs46;
    const double clhs58 =             C(1,2)*DN(2,1);
    const double clhs59 =             C(2,2)*DN(2,0) + clhs58;
    const double clhs60 =             DN(2,1)*clhs23;
    const double clhs61 =             DN(0,0)*N[2];
    const double clhs62 =             DN(2,0)*N[0];
    const double clhs63 =             C(0,1)*DN(0,0) + clhs21;
    const double clhs64 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
    const double clhs65 =             pow(DN(0,1), 2);
    const double clhs66 =             C(0,1)*DN(1,0) + clhs39;
    const double clhs67 =             DN(0,1)*clhs7;
    const double clhs68 =             DN(1,0)*clhs67;
    const double clhs69 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
    const double clhs70 =             DN(0,1)*DN(1,1);
    const double clhs71 =             clhs7*clhs70;
    const double clhs72 =             DN(0,1)*N[1];
    const double clhs73 =             DN(1,1)*N[0];
    const double clhs74 =             C(0,1)*DN(2,0) + clhs58;
    const double clhs75 =             DN(2,0)*clhs67;
    const double clhs76 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
    const double clhs77 =             DN(0,1)*DN(2,1);
    const double clhs78 =             clhs7*clhs77;
    const double clhs79 =             DN(0,1)*N[2];
    const double clhs80 =             DN(2,1)*N[0];
    const double clhs81 =             N[0] + clhs15*(1.0*clhs12 + clhs14 + clhs17);
    const double clhs82 =             1.0*clhs15;
    const double clhs83 =             1.0*DN(0,0)*clhs15;
    const double clhs84 =             1.0*DN(0,1)*clhs15;
    const double clhs85 =             clhs82*(clhs29 + clhs70);
    const double clhs86 =             clhs82*(clhs48 + clhs77);
    const double clhs87 =             1.0*clhs34;
    const double clhs88 =             clhs15*clhs87;
    const double clhs89 =             1.0*clhs33;
    const double clhs90 =             clhs15*clhs89;
    const double clhs91 =             N[1]*clhs10 - clhs13*clhs88 + clhs13*clhs90 + clhs31 + clhs32;
    const double clhs92 =             pow(DN(1,0), 2);
    const double clhs93 =             pow(N[1], 2);
    const double clhs94 =             K_darcy*clhs93 + N[1]*clhs33 - clhs36*clhs88 + clhs36*clhs90 + clhs9*clhs93;
    const double clhs95 =             DN(1,0)*clhs7;
    const double clhs96 =             DN(1,1)*clhs95;
    const double clhs97 =             -N[1] - clhs88 + clhs90;
    const double clhs98 =             DN(1,0)*DN(2,0);
    const double clhs99 =             clhs7*clhs98;
    const double clhs100 =             N[2]*clhs34;
    const double clhs101 =             N[2]*clhs35;
    const double clhs102 =             N[1]*clhs52 + clhs100 + clhs101 - clhs55*clhs88 + clhs55*clhs90;
    const double clhs103 =             DN(2,1)*clhs95;
    const double clhs104 =             DN(1,0)*N[2];
    const double clhs105 =             DN(2,0)*N[1];
    const double clhs106 =             pow(DN(1,1), 2);
    const double clhs107 =             DN(2,0)*clhs7;
    const double clhs108 =             DN(1,1)*clhs107;
    const double clhs109 =             DN(1,1)*DN(2,1);
    const double clhs110 =             clhs109*clhs7;
    const double clhs111 =             DN(1,1)*N[2];
    const double clhs112 =             DN(2,1)*N[1];
    const double clhs113 =             1.0*DN(1,0)*clhs15;
    const double clhs114 =             1.0*DN(1,1)*clhs15;
    const double clhs115 =             N[1] + clhs15*(1.0*clhs35 + clhs87 + clhs89);
    const double clhs116 =             clhs82*(clhs109 + clhs98);
    const double clhs117 =             1.0*clhs53;
    const double clhs118 =             clhs117*clhs15;
    const double clhs119 =             1.0*clhs52;
    const double clhs120 =             clhs119*clhs15;
    const double clhs121 =             N[2]*clhs10 - clhs118*clhs13 + clhs120*clhs13 + clhs50 + clhs51;
    const double clhs122 =             N[2]*clhs33 + clhs100 + clhs101 - clhs118*clhs36 + clhs120*clhs36;
    const double clhs123 =             pow(DN(2,0), 2);
    const double clhs124 =             pow(N[2], 2);
    const double clhs125 =             K_darcy*clhs124 + N[2]*clhs52 - clhs118*clhs55 + clhs120*clhs55 + clhs124*clhs9;
    const double clhs126 =             DN(2,1)*clhs107;
    const double clhs127 =             -N[2] - clhs118 + clhs120;
    const double clhs128 =             pow(DN(2,1), 2);
    const double clhs129 =             1.0*DN(2,0)*clhs15;
    const double clhs130 =             1.0*DN(2,1)*clhs15;
    const double clhs131 =             N[2] + clhs15*(clhs117 + clhs119 + 1.0*clhs54);
    lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
    lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
    lhs(0,2)=DN(0,0)*clhs25;
    lhs(0,3)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs30 + clhs37;
    lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + clhs41;
    lhs(0,5)=DN(1,0)*clhs18 - clhs42 - clhs43*clhs44;
    lhs(0,6)=DN(0,0)*clhs45 + DN(0,1)*clhs47 + clhs49 + clhs56;
    lhs(0,7)=DN(0,0)*clhs57 + DN(0,1)*clhs59 + clhs60;
    lhs(0,8)=DN(2,0)*clhs18 - clhs44*clhs62 - clhs61;
    lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs63 + clhs24;
    lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs64 + clhs19 + clhs65*clhs7;
    lhs(1,2)=DN(0,1)*clhs25;
    lhs(1,3)=DN(0,0)*clhs28 + DN(0,1)*clhs66 + clhs68;
    lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs69 + clhs37 + clhs71;
    lhs(1,5)=DN(1,1)*clhs18 - clhs44*clhs73 - clhs72;
    lhs(1,6)=DN(0,0)*clhs47 + DN(0,1)*clhs74 + clhs75;
    lhs(1,7)=DN(0,0)*clhs59 + DN(0,1)*clhs76 + clhs56 + clhs78;
    lhs(1,8)=DN(2,1)*clhs18 - clhs44*clhs80 - clhs79;
    lhs(2,0)=DN(0,0)*clhs81;
    lhs(2,1)=DN(0,1)*clhs81;
    lhs(2,2)=clhs82*(clhs3 + clhs65);
    lhs(2,3)=clhs36*clhs83 + clhs43;
    lhs(2,4)=clhs36*clhs84 + clhs73;
    lhs(2,5)=clhs85;
    lhs(2,6)=clhs55*clhs83 + clhs62;
    lhs(2,7)=clhs55*clhs84 + clhs80;
    lhs(2,8)=clhs86;
    lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs30 + clhs91;
    lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs68;
    lhs(3,2)=DN(0,0)*clhs90 - clhs42*clhs44 - clhs43;
    lhs(3,3)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs7*clhs92 + clhs94;
    lhs(3,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + clhs96;
    lhs(3,5)=DN(1,0)*clhs97;
    lhs(3,6)=DN(1,0)*clhs45 + DN(1,1)*clhs47 + clhs102 + clhs99;
    lhs(3,7)=DN(1,0)*clhs57 + DN(1,1)*clhs59 + clhs103;
    lhs(3,8)=DN(2,0)*clhs90 - clhs104 - clhs105*clhs44;
    lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs63 + clhs41;
    lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs64 + clhs71 + clhs91;
    lhs(4,2)=DN(0,1)*clhs90 - clhs44*clhs72 - clhs73;
    lhs(4,3)=DN(1,0)*clhs28 + DN(1,1)*clhs66 + clhs96;
    lhs(4,4)=DN(1,0)*clhs40 + DN(1,1)*clhs69 + clhs106*clhs7 + clhs94;
    lhs(4,5)=DN(1,1)*clhs97;
    lhs(4,6)=DN(1,0)*clhs47 + DN(1,1)*clhs74 + clhs108;
    lhs(4,7)=DN(1,0)*clhs59 + DN(1,1)*clhs76 + clhs102 + clhs110;
    lhs(4,8)=DN(2,1)*clhs90 - clhs111 - clhs112*clhs44;
    lhs(5,0)=clhs113*clhs13 + clhs42;
    lhs(5,1)=clhs114*clhs13 + clhs72;
    lhs(5,2)=clhs85;
    lhs(5,3)=DN(1,0)*clhs115;
    lhs(5,4)=DN(1,1)*clhs115;
    lhs(5,5)=clhs82*(clhs106 + clhs92);
    lhs(5,6)=clhs105 + clhs113*clhs55;
    lhs(5,7)=clhs112 + clhs114*clhs55;
    lhs(5,8)=clhs116;
    lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs121 + clhs49;
    lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs75;
    lhs(6,2)=DN(0,0)*clhs120 - clhs44*clhs61 - clhs62;
    lhs(6,3)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs122 + clhs99;
    lhs(6,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + clhs108;
    lhs(6,5)=DN(1,0)*clhs120 - clhs104*clhs44 - clhs105;
    lhs(6,6)=DN(2,0)*clhs45 + DN(2,1)*clhs47 + clhs123*clhs7 + clhs125;
    lhs(6,7)=DN(2,0)*clhs57 + DN(2,1)*clhs59 + clhs126;
    lhs(6,8)=DN(2,0)*clhs127;
    lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs63 + clhs60;
    lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs64 + clhs121 + clhs78;
    lhs(7,2)=DN(0,1)*clhs120 - clhs44*clhs79 - clhs80;
    lhs(7,3)=DN(2,0)*clhs28 + DN(2,1)*clhs66 + clhs103;
    lhs(7,4)=DN(2,0)*clhs40 + DN(2,1)*clhs69 + clhs110 + clhs122;
    lhs(7,5)=DN(1,1)*clhs120 - clhs111*clhs44 - clhs112;
    lhs(7,6)=DN(2,0)*clhs47 + DN(2,1)*clhs74 + clhs126;
    lhs(7,7)=DN(2,0)*clhs59 + DN(2,1)*clhs76 + clhs125 + clhs128*clhs7;
    lhs(7,8)=DN(2,1)*clhs127;
    lhs(8,0)=clhs129*clhs13 + clhs61;
    lhs(8,1)=clhs13*clhs130 + clhs79;
    lhs(8,2)=clhs86;
    lhs(8,3)=clhs104 + clhs129*clhs36;
    lhs(8,4)=clhs111 + clhs130*clhs36;
    lhs(8,5)=clhs116;
    lhs(8,6)=DN(2,0)*clhs131;
    lhs(8,7)=DN(2,1)*clhs131;
    lhs(8,8)=clhs82*(clhs123 + clhs128);

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

    //const double dyn_tau = rData.DynamicTau;
    double dyn_tau = rData.DynamicTau;
    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix &C = rData.C;

    // Get shape function values
    const auto &N = rData.N;
    const auto &DN = rData.DN_DX;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

/*     unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &lhs = rData.lhs;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
    const double clhs1 =             C(0,3)*DN(0,0);
    const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
    const double clhs3 =             C(0,5)*DN(0,0);
    const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
    const double clhs5 =             pow(DN(0,0), 2);
    const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double clhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double clhs9 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
    const double clhs10 =             clhs9*h/stab_c1 + mu;
    const double clhs11 =             pow(N[0], 2);
    const double clhs12 =             bdf0*rho;
    const double clhs13 =             rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
    const double clhs14 =             K_darcy*N[0];
    const double clhs15 =             N[0]*clhs12;
    const double clhs16 =             clhs13 + clhs14 + clhs15;
    const double clhs17 =             1.0*clhs14; // Darcy
    const double clhs18 =             1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double clhs19 =             clhs17*clhs18; 
    const double clhs20 =             1.0*clhs13; // Convection
    const double clhs21 =             clhs18*clhs20; // Convection
    const double clhs22 =             K_darcy*clhs11 + N[0]*clhs13 + clhs11*clhs12 - clhs16*clhs19 + clhs16*clhs21; // Convection and Darcy
    const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
    const double clhs24 =             C(1,3)*DN(0,1);
    const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
    const double clhs26 =             C(3,5)*DN(0,0);
    const double clhs27 =             C(4,5)*DN(0,2);
    const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
    const double clhs29 =             DN(0,0)*clhs10;
    const double clhs30 =             DN(0,1)*clhs29;
    const double clhs31 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
    const double clhs32 =             C(3,4)*DN(0,1);
    const double clhs33 =             C(2,3)*DN(0,2) + clhs26 + clhs32;
    const double clhs34 =             C(2,5)*DN(0,2);
    const double clhs35 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
    const double clhs36 =             DN(0,2)*clhs29;
    const double clhs37 =             -N[0] - clhs19 + clhs21;
    const double clhs38 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
    const double clhs39 =             C(0,3)*DN(1,0);
    const double clhs40 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs39;
    const double clhs41 =             C(0,5)*DN(1,0);
    const double clhs42 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs41;
    const double clhs43 =             DN(0,0)*DN(1,0);
    const double clhs44 =             clhs10*clhs43;
    const double clhs45 =             N[1]*clhs14;
    const double clhs46 =             N[1]*clhs15;
    const double clhs47 =             rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
    const double clhs48 =             K_darcy*N[1];
    const double clhs49 =             N[1]*clhs12;
    const double clhs50 =             clhs47 + clhs48 + clhs49;
    const double clhs51 =             N[0]*clhs47 - clhs19*clhs50 + clhs21*clhs50 + clhs45 + clhs46;
    const double clhs52 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs39;
    const double clhs53 =             C(1,3)*DN(1,1);
    const double clhs54 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs53;
    const double clhs55 =             C(3,5)*DN(1,0);
    const double clhs56 =             C(4,5)*DN(1,2);
    const double clhs57 =             C(1,5)*DN(1,1) + clhs55 + clhs56;
    const double clhs58 =             DN(1,1)*clhs29;
    const double clhs59 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs41;
    const double clhs60 =             C(3,4)*DN(1,1);
    const double clhs61 =             C(2,3)*DN(1,2) + clhs55 + clhs60;
    const double clhs62 =             C(2,5)*DN(1,2);
    const double clhs63 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs62;
    const double clhs64 =             DN(1,2)*clhs29;
    const double clhs65 =             DN(0,0)*N[1];
    const double clhs66 =             DN(1,0)*N[0];
    const double clhs67 =             1.0*K_darcy*clhs18;
    const double clhs68 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
    const double clhs69 =             C(0,3)*DN(2,0);
    const double clhs70 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs69;
    const double clhs71 =             C(0,5)*DN(2,0);
    const double clhs72 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs71;
    const double clhs73 =             DN(0,0)*DN(2,0);
    const double clhs74 =             clhs10*clhs73;
    const double clhs75 =             N[2]*clhs14;
    const double clhs76 =             N[2]*clhs15;
    const double clhs77 =             rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
    const double clhs78 =             K_darcy*N[2];
    const double clhs79 =             N[2]*clhs12;
    const double clhs80 =             clhs77 + clhs78 + clhs79;
    const double clhs81 =             N[0]*clhs77 - clhs19*clhs80 + clhs21*clhs80 + clhs75 + clhs76;
    const double clhs82 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs69;
    const double clhs83 =             C(1,3)*DN(2,1);
    const double clhs84 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs83;
    const double clhs85 =             C(3,5)*DN(2,0);
    const double clhs86 =             C(4,5)*DN(2,2);
    const double clhs87 =             C(1,5)*DN(2,1) + clhs85 + clhs86;
    const double clhs88 =             DN(2,1)*clhs29;
    const double clhs89 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs71;
    const double clhs90 =             C(3,4)*DN(2,1);
    const double clhs91 =             C(2,3)*DN(2,2) + clhs85 + clhs90;
    const double clhs92 =             C(2,5)*DN(2,2);
    const double clhs93 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs92;
    const double clhs94 =             DN(2,2)*clhs29;
    const double clhs95 =             DN(0,0)*N[2];
    const double clhs96 =             DN(2,0)*N[0];
    const double clhs97 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
    const double clhs98 =             C(0,3)*DN(3,0);
    const double clhs99 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs98;
    const double clhs100 =             C(0,5)*DN(3,0);
    const double clhs101 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs100;
    const double clhs102 =             DN(0,0)*DN(3,0);
    const double clhs103 =             clhs10*clhs102;
    const double clhs104 =             N[3]*clhs14;
    const double clhs105 =             N[3]*clhs15;
    const double clhs106 =             rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
    const double clhs107 =             K_darcy*N[3];
    const double clhs108 =             N[3]*clhs12;
    const double clhs109 =             clhs106 + clhs107 + clhs108;
    const double clhs110 =             N[0]*clhs106 + clhs104 + clhs105 - clhs109*clhs19 + clhs109*clhs21;
    const double clhs111 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs98;
    const double clhs112 =             C(1,3)*DN(3,1);
    const double clhs113 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs112;
    const double clhs114 =             C(3,5)*DN(3,0);
    const double clhs115 =             C(4,5)*DN(3,2);
    const double clhs116 =             C(1,5)*DN(3,1) + clhs114 + clhs115;
    const double clhs117 =             DN(3,1)*clhs29;
    const double clhs118 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs100;
    const double clhs119 =             C(3,4)*DN(3,1);
    const double clhs120 =             C(2,3)*DN(3,2) + clhs114 + clhs119;
    const double clhs121 =             C(2,5)*DN(3,2);
    const double clhs122 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs121;
    const double clhs123 =             DN(3,2)*clhs29;
    const double clhs124 =             DN(0,0)*N[3];
    const double clhs125 =             DN(3,0)*N[0];
    const double clhs126 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
    const double clhs127 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
    const double clhs128 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
    const double clhs129 =             C(1,4)*DN(0,1);
    const double clhs130 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs129;
    const double clhs131 =             pow(DN(0,1), 2);
    const double clhs132 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs129;
    const double clhs133 =             C(2,4)*DN(0,2);
    const double clhs134 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs133;
    const double clhs135 =             DN(0,1)*clhs10;
    const double clhs136 =             DN(0,2)*clhs135;
    const double clhs137 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs53;
    const double clhs138 =             C(0,4)*DN(1,0) + clhs56 + clhs60;
    const double clhs139 =             DN(1,0)*clhs135;
    const double clhs140 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
    const double clhs141 =             C(1,4)*DN(1,1);
    const double clhs142 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs141;
    const double clhs143 =             DN(0,1)*DN(1,1);
    const double clhs144 =             clhs10*clhs143;
    const double clhs145 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs141;
    const double clhs146 =             C(2,4)*DN(1,2);
    const double clhs147 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs146;
    const double clhs148 =             DN(1,2)*clhs135;
    const double clhs149 =             DN(0,1)*N[1];
    const double clhs150 =             DN(1,1)*N[0];
    const double clhs151 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs83;
    const double clhs152 =             C(0,4)*DN(2,0) + clhs86 + clhs90;
    const double clhs153 =             DN(2,0)*clhs135;
    const double clhs154 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
    const double clhs155 =             C(1,4)*DN(2,1);
    const double clhs156 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs155;
    const double clhs157 =             DN(0,1)*DN(2,1);
    const double clhs158 =             clhs10*clhs157;
    const double clhs159 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs155;
    const double clhs160 =             C(2,4)*DN(2,2);
    const double clhs161 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs160;
    const double clhs162 =             DN(2,2)*clhs135;
    const double clhs163 =             DN(0,1)*N[2];
    const double clhs164 =             DN(2,1)*N[0];
    const double clhs165 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs112;
    const double clhs166 =             C(0,4)*DN(3,0) + clhs115 + clhs119;
    const double clhs167 =             DN(3,0)*clhs135;
    const double clhs168 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
    const double clhs169 =             C(1,4)*DN(3,1);
    const double clhs170 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs169;
    const double clhs171 =             DN(0,1)*DN(3,1);
    const double clhs172 =             clhs10*clhs171;
    const double clhs173 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs169;
    const double clhs174 =             C(2,4)*DN(3,2);
    const double clhs175 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs174;
    const double clhs176 =             DN(3,2)*clhs135;
    const double clhs177 =             DN(0,1)*N[3];
    const double clhs178 =             DN(3,1)*N[0];
    const double clhs179 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
    const double clhs180 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs133;
    const double clhs181 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
    const double clhs182 =             pow(DN(0,2), 2);
    const double clhs183 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs62;
    const double clhs184 =             DN(0,2)*clhs10;
    const double clhs185 =             DN(1,0)*clhs184;
    const double clhs186 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs146;
    const double clhs187 =             DN(1,1)*clhs184;
    const double clhs188 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
    const double clhs189 =             DN(0,2)*DN(1,2);
    const double clhs190 =             clhs10*clhs189;
    const double clhs191 =             DN(0,2)*N[1];
    const double clhs192 =             DN(1,2)*N[0];
    const double clhs193 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs92;
    const double clhs194 =             DN(2,0)*clhs184;
    const double clhs195 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs160;
    const double clhs196 =             DN(2,1)*clhs184;
    const double clhs197 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
    const double clhs198 =             DN(0,2)*DN(2,2);
    const double clhs199 =             clhs10*clhs198;
    const double clhs200 =             DN(0,2)*N[2];
    const double clhs201 =             DN(2,2)*N[0];
    const double clhs202 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs121;
    const double clhs203 =             DN(3,0)*clhs184;
    const double clhs204 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs174;
    const double clhs205 =             DN(3,1)*clhs184;
    const double clhs206 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
    const double clhs207 =             DN(0,2)*DN(3,2);
    const double clhs208 =             clhs10*clhs207;
    const double clhs209 =             DN(0,2)*N[3];
    const double clhs210 =             DN(3,2)*N[0];
    const double clhs211 =             N[0] + clhs18*(1.0*clhs15 + clhs17 + clhs20);
    const double clhs212 =             1.0*clhs18;
    const double clhs213 =             1.0*DN(0,0)*clhs18;
    const double clhs214 =             1.0*DN(0,1)*clhs18;
    const double clhs215 =             1.0*DN(0,2)*clhs18;
    const double clhs216 =             clhs212*(clhs143 + clhs189 + clhs43);
    const double clhs217 =             clhs212*(clhs157 + clhs198 + clhs73);
    const double clhs218 =             clhs212*(clhs102 + clhs171 + clhs207);
    const double clhs219 =             1.0*clhs48;
    const double clhs220 =             clhs18*clhs219;
    const double clhs221 =             1.0*clhs47;
    const double clhs222 =             clhs18*clhs221;
    const double clhs223 =             N[1]*clhs13 - clhs16*clhs220 + clhs16*clhs222 + clhs45 + clhs46;
    const double clhs224 =             pow(DN(1,0), 2);
    const double clhs225 =             pow(N[1], 2);
    const double clhs226 =             K_darcy*clhs225 + N[1]*clhs47 + clhs12*clhs225 - clhs220*clhs50 + clhs222*clhs50;
    const double clhs227 =             DN(1,0)*clhs10;
    const double clhs228 =             DN(1,1)*clhs227;
    const double clhs229 =             DN(1,2)*clhs227;
    const double clhs230 =             -N[1] - clhs220 + clhs222;
    const double clhs231 =             DN(1,0)*DN(2,0);
    const double clhs232 =             clhs10*clhs231;
    const double clhs233 =             N[2]*clhs48;
    const double clhs234 =             N[2]*clhs49;
    const double clhs235 =             N[1]*clhs77 - clhs220*clhs80 + clhs222*clhs80 + clhs233 + clhs234;
    const double clhs236 =             DN(2,1)*clhs227;
    const double clhs237 =             DN(2,2)*clhs227;
    const double clhs238 =             DN(1,0)*N[2];
    const double clhs239 =             DN(2,0)*N[1];
    const double clhs240 =             DN(1,0)*DN(3,0);
    const double clhs241 =             clhs10*clhs240;
    const double clhs242 =             N[3]*clhs48;
    const double clhs243 =             N[3]*clhs49;
    const double clhs244 =             N[1]*clhs106 - clhs109*clhs220 + clhs109*clhs222 + clhs242 + clhs243;
    const double clhs245 =             DN(3,1)*clhs227;
    const double clhs246 =             DN(3,2)*clhs227;
    const double clhs247 =             DN(1,0)*N[3];
    const double clhs248 =             DN(3,0)*N[1];
    const double clhs249 =             pow(DN(1,1), 2);
    const double clhs250 =             DN(1,1)*clhs10;
    const double clhs251 =             DN(1,2)*clhs250;
    const double clhs252 =             DN(2,0)*clhs250;
    const double clhs253 =             DN(1,1)*DN(2,1);
    const double clhs254 =             clhs10*clhs253;
    const double clhs255 =             DN(2,2)*clhs250;
    const double clhs256 =             DN(1,1)*N[2];
    const double clhs257 =             DN(2,1)*N[1];
    const double clhs258 =             DN(3,0)*clhs250;
    const double clhs259 =             DN(1,1)*DN(3,1);
    const double clhs260 =             clhs10*clhs259;
    const double clhs261 =             DN(3,2)*clhs250;
    const double clhs262 =             DN(1,1)*N[3];
    const double clhs263 =             DN(3,1)*N[1];
    const double clhs264 =             pow(DN(1,2), 2);
    const double clhs265 =             DN(1,2)*clhs10;
    const double clhs266 =             DN(2,0)*clhs265;
    const double clhs267 =             DN(2,1)*clhs265;
    const double clhs268 =             DN(1,2)*DN(2,2);
    const double clhs269 =             clhs10*clhs268;
    const double clhs270 =             DN(1,2)*N[2];
    const double clhs271 =             DN(2,2)*N[1];
    const double clhs272 =             DN(3,0)*clhs265;
    const double clhs273 =             DN(3,1)*clhs265;
    const double clhs274 =             DN(1,2)*DN(3,2);
    const double clhs275 =             clhs10*clhs274;
    const double clhs276 =             DN(1,2)*N[3];
    const double clhs277 =             DN(3,2)*N[1];
    const double clhs278 =             1.0*DN(1,0)*clhs18;
    const double clhs279 =             1.0*DN(1,1)*clhs18;
    const double clhs280 =             1.0*DN(1,2)*clhs18;
    const double clhs281 =             N[1] + clhs18*(clhs219 + clhs221 + 1.0*clhs49);
    const double clhs282 =             clhs212*(clhs231 + clhs253 + clhs268);
    const double clhs283 =             clhs212*(clhs240 + clhs259 + clhs274);
    const double clhs284 =             1.0*clhs78;
    const double clhs285 =             clhs18*clhs284;
    const double clhs286 =             1.0*clhs77;
    const double clhs287 =             clhs18*clhs286;
    const double clhs288 =             N[2]*clhs13 - clhs16*clhs285 + clhs16*clhs287 + clhs75 + clhs76;
    const double clhs289 =             N[2]*clhs47 + clhs233 + clhs234 - clhs285*clhs50 + clhs287*clhs50;
    const double clhs290 =             pow(DN(2,0), 2);
    const double clhs291 =             pow(N[2], 2);
    const double clhs292 =             K_darcy*clhs291 + N[2]*clhs77 + clhs12*clhs291 - clhs285*clhs80 + clhs287*clhs80;
    const double clhs293 =             DN(2,0)*clhs10;
    const double clhs294 =             DN(2,1)*clhs293;
    const double clhs295 =             DN(2,2)*clhs293;
    const double clhs296 =             -N[2] - clhs285 + clhs287;
    const double clhs297 =             DN(2,0)*DN(3,0);
    const double clhs298 =             clhs10*clhs297;
    const double clhs299 =             N[3]*clhs78;
    const double clhs300 =             N[3]*clhs79;
    const double clhs301 =             N[2]*clhs106 - clhs109*clhs285 + clhs109*clhs287 + clhs299 + clhs300;
    const double clhs302 =             DN(3,1)*clhs293;
    const double clhs303 =             DN(3,2)*clhs293;
    const double clhs304 =             DN(2,0)*N[3];
    const double clhs305 =             DN(3,0)*N[2];
    const double clhs306 =             pow(DN(2,1), 2);
    const double clhs307 =             DN(2,1)*clhs10;
    const double clhs308 =             DN(2,2)*clhs307;
    const double clhs309 =             DN(3,0)*clhs307;
    const double clhs310 =             DN(2,1)*DN(3,1);
    const double clhs311 =             clhs10*clhs310;
    const double clhs312 =             DN(3,2)*clhs307;
    const double clhs313 =             DN(2,1)*N[3];
    const double clhs314 =             DN(3,1)*N[2];
    const double clhs315 =             pow(DN(2,2), 2);
    const double clhs316 =             DN(2,2)*clhs10;
    const double clhs317 =             DN(3,0)*clhs316;
    const double clhs318 =             DN(3,1)*clhs316;
    const double clhs319 =             DN(2,2)*DN(3,2);
    const double clhs320 =             clhs10*clhs319;
    const double clhs321 =             DN(2,2)*N[3];
    const double clhs322 =             DN(3,2)*N[2];
    const double clhs323 =             1.0*DN(2,0)*clhs18;
    const double clhs324 =             1.0*DN(2,1)*clhs18;
    const double clhs325 =             1.0*DN(2,2)*clhs18;
    const double clhs326 =             N[2] + clhs18*(clhs284 + clhs286 + 1.0*clhs79);
    const double clhs327 =             clhs212*(clhs297 + clhs310 + clhs319);
    const double clhs328 =             1.0*clhs107;
    const double clhs329 =             clhs18*clhs328;
    const double clhs330 =             1.0*clhs106;
    const double clhs331 =             clhs18*clhs330;
    const double clhs332 =             N[3]*clhs13 + clhs104 + clhs105 - clhs16*clhs329 + clhs16*clhs331;
    const double clhs333 =             N[3]*clhs47 + clhs242 + clhs243 - clhs329*clhs50 + clhs331*clhs50;
    const double clhs334 =             N[3]*clhs77 + clhs299 + clhs300 - clhs329*clhs80 + clhs331*clhs80;
    const double clhs335 =             pow(DN(3,0), 2);
    const double clhs336 =             pow(N[3], 2);
    const double clhs337 =             K_darcy*clhs336 + N[3]*clhs106 - clhs109*clhs329 + clhs109*clhs331 + clhs12*clhs336;
    const double clhs338 =             DN(3,0)*clhs10;
    const double clhs339 =             DN(3,1)*clhs338;
    const double clhs340 =             DN(3,2)*clhs338;
    const double clhs341 =             -N[3] - clhs329 + clhs331;
    const double clhs342 =             pow(DN(3,1), 2);
    const double clhs343 =             DN(3,1)*DN(3,2)*clhs10;
    const double clhs344 =             pow(DN(3,2), 2);
    const double clhs345 =             1.0*DN(3,0)*clhs18;
    const double clhs346 =             1.0*DN(3,1)*clhs18;
    const double clhs347 =             1.0*DN(3,2)*clhs18;
    const double clhs348 =             N[3] + clhs18*(1.0*clhs108 + clhs328 + clhs330);
    lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
    lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
    lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
    lhs(0,3)=DN(0,0)*clhs37;
    lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + DN(0,2)*clhs42 + clhs44 + clhs51;
    lhs(0,5)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + DN(0,2)*clhs57 + clhs58;
    lhs(0,6)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + DN(0,2)*clhs63 + clhs64;
    lhs(0,7)=DN(1,0)*clhs21 - clhs65 - clhs66*clhs67;
    lhs(0,8)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs74 + clhs81;
    lhs(0,9)=DN(0,0)*clhs82 + DN(0,1)*clhs84 + DN(0,2)*clhs87 + clhs88;
    lhs(0,10)=DN(0,0)*clhs89 + DN(0,1)*clhs91 + DN(0,2)*clhs93 + clhs94;
    lhs(0,11)=DN(2,0)*clhs21 - clhs67*clhs96 - clhs95;
    lhs(0,12)=DN(0,0)*clhs97 + DN(0,1)*clhs99 + DN(0,2)*clhs101 + clhs103 + clhs110;
    lhs(0,13)=DN(0,0)*clhs111 + DN(0,1)*clhs113 + DN(0,2)*clhs116 + clhs117;
    lhs(0,14)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs123;
    lhs(0,15)=DN(3,0)*clhs21 - clhs124 - clhs125*clhs67;
    lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs126 + DN(0,2)*clhs127 + clhs30;
    lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs128 + DN(0,2)*clhs130 + clhs10*clhs131 + clhs22;
    lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs136;
    lhs(1,3)=DN(0,1)*clhs37;
    lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
    lhs(1,5)=DN(0,0)*clhs54 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs144 + clhs51;
    lhs(1,6)=DN(0,0)*clhs61 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148;
    lhs(1,7)=DN(1,1)*clhs21 - clhs149 - clhs150*clhs67;
    lhs(1,8)=DN(0,0)*clhs70 + DN(0,1)*clhs151 + DN(0,2)*clhs152 + clhs153;
    lhs(1,9)=DN(0,0)*clhs84 + DN(0,1)*clhs154 + DN(0,2)*clhs156 + clhs158 + clhs81;
    lhs(1,10)=DN(0,0)*clhs91 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162;
    lhs(1,11)=DN(2,1)*clhs21 - clhs163 - clhs164*clhs67;
    lhs(1,12)=DN(0,0)*clhs99 + DN(0,1)*clhs165 + DN(0,2)*clhs166 + clhs167;
    lhs(1,13)=DN(0,0)*clhs113 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs110 + clhs172;
    lhs(1,14)=DN(0,0)*clhs120 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176;
    lhs(1,15)=DN(3,1)*clhs21 - clhs177 - clhs178*clhs67;
    lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs179 + clhs36;
    lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs130 + DN(0,2)*clhs180 + clhs136;
    lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs134 + DN(0,2)*clhs181 + clhs10*clhs182 + clhs22;
    lhs(2,3)=DN(0,2)*clhs37;
    lhs(2,4)=DN(0,0)*clhs42 + DN(0,1)*clhs138 + DN(0,2)*clhs183 + clhs185;
    lhs(2,5)=DN(0,0)*clhs57 + DN(0,1)*clhs142 + DN(0,2)*clhs186 + clhs187;
    lhs(2,6)=DN(0,0)*clhs63 + DN(0,1)*clhs147 + DN(0,2)*clhs188 + clhs190 + clhs51;
    lhs(2,7)=DN(1,2)*clhs21 - clhs191 - clhs192*clhs67;
    lhs(2,8)=DN(0,0)*clhs72 + DN(0,1)*clhs152 + DN(0,2)*clhs193 + clhs194;
    lhs(2,9)=DN(0,0)*clhs87 + DN(0,1)*clhs156 + DN(0,2)*clhs195 + clhs196;
    lhs(2,10)=DN(0,0)*clhs93 + DN(0,1)*clhs161 + DN(0,2)*clhs197 + clhs199 + clhs81;
    lhs(2,11)=DN(2,2)*clhs21 - clhs200 - clhs201*clhs67;
    lhs(2,12)=DN(0,0)*clhs101 + DN(0,1)*clhs166 + DN(0,2)*clhs202 + clhs203;
    lhs(2,13)=DN(0,0)*clhs116 + DN(0,1)*clhs170 + DN(0,2)*clhs204 + clhs205;
    lhs(2,14)=DN(0,0)*clhs122 + DN(0,1)*clhs175 + DN(0,2)*clhs206 + clhs110 + clhs208;
    lhs(2,15)=DN(3,2)*clhs21 - clhs209 - clhs210*clhs67;
    lhs(3,0)=DN(0,0)*clhs211;
    lhs(3,1)=DN(0,1)*clhs211;
    lhs(3,2)=DN(0,2)*clhs211;
    lhs(3,3)=clhs212*(clhs131 + clhs182 + clhs5);
    lhs(3,4)=clhs213*clhs50 + clhs66;
    lhs(3,5)=clhs150 + clhs214*clhs50;
    lhs(3,6)=clhs192 + clhs215*clhs50;
    lhs(3,7)=clhs216;
    lhs(3,8)=clhs213*clhs80 + clhs96;
    lhs(3,9)=clhs164 + clhs214*clhs80;
    lhs(3,10)=clhs201 + clhs215*clhs80;
    lhs(3,11)=clhs217;
    lhs(3,12)=clhs109*clhs213 + clhs125;
    lhs(3,13)=clhs109*clhs214 + clhs178;
    lhs(3,14)=clhs109*clhs215 + clhs210;
    lhs(3,15)=clhs218;
    lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs223 + clhs44;
    lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs139;
    lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs185;
    lhs(4,3)=DN(0,0)*clhs222 - clhs65*clhs67 - clhs66;
    lhs(4,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + DN(1,2)*clhs42 + clhs10*clhs224 + clhs226;
    lhs(4,5)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + DN(1,2)*clhs57 + clhs228;
    lhs(4,6)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + DN(1,2)*clhs63 + clhs229;
    lhs(4,7)=DN(1,0)*clhs230;
    lhs(4,8)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs232 + clhs235;
    lhs(4,9)=DN(1,0)*clhs82 + DN(1,1)*clhs84 + DN(1,2)*clhs87 + clhs236;
    lhs(4,10)=DN(1,0)*clhs89 + DN(1,1)*clhs91 + DN(1,2)*clhs93 + clhs237;
    lhs(4,11)=DN(2,0)*clhs222 - clhs238 - clhs239*clhs67;
    lhs(4,12)=DN(1,0)*clhs97 + DN(1,1)*clhs99 + DN(1,2)*clhs101 + clhs241 + clhs244;
    lhs(4,13)=DN(1,0)*clhs111 + DN(1,1)*clhs113 + DN(1,2)*clhs116 + clhs245;
    lhs(4,14)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs246;
    lhs(4,15)=DN(3,0)*clhs222 - clhs247 - clhs248*clhs67;
    lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs126 + DN(1,2)*clhs127 + clhs58;
    lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs128 + DN(1,2)*clhs130 + clhs144 + clhs223;
    lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs187;
    lhs(5,3)=DN(0,1)*clhs222 - clhs149*clhs67 - clhs150;
    lhs(5,4)=DN(1,0)*clhs40 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs228;
    lhs(5,5)=DN(1,0)*clhs54 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs10*clhs249 + clhs226;
    lhs(5,6)=DN(1,0)*clhs61 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs251;
    lhs(5,7)=DN(1,1)*clhs230;
    lhs(5,8)=DN(1,0)*clhs70 + DN(1,1)*clhs151 + DN(1,2)*clhs152 + clhs252;
    lhs(5,9)=DN(1,0)*clhs84 + DN(1,1)*clhs154 + DN(1,2)*clhs156 + clhs235 + clhs254;
    lhs(5,10)=DN(1,0)*clhs91 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs255;
    lhs(5,11)=DN(2,1)*clhs222 - clhs256 - clhs257*clhs67;
    lhs(5,12)=DN(1,0)*clhs99 + DN(1,1)*clhs165 + DN(1,2)*clhs166 + clhs258;
    lhs(5,13)=DN(1,0)*clhs113 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs244 + clhs260;
    lhs(5,14)=DN(1,0)*clhs120 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs261;
    lhs(5,15)=DN(3,1)*clhs222 - clhs262 - clhs263*clhs67;
    lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs179 + clhs64;
    lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs130 + DN(1,2)*clhs180 + clhs148;
    lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs134 + DN(1,2)*clhs181 + clhs190 + clhs223;
    lhs(6,3)=DN(0,2)*clhs222 - clhs191*clhs67 - clhs192;
    lhs(6,4)=DN(1,0)*clhs42 + DN(1,1)*clhs138 + DN(1,2)*clhs183 + clhs229;
    lhs(6,5)=DN(1,0)*clhs57 + DN(1,1)*clhs142 + DN(1,2)*clhs186 + clhs251;
    lhs(6,6)=DN(1,0)*clhs63 + DN(1,1)*clhs147 + DN(1,2)*clhs188 + clhs10*clhs264 + clhs226;
    lhs(6,7)=DN(1,2)*clhs230;
    lhs(6,8)=DN(1,0)*clhs72 + DN(1,1)*clhs152 + DN(1,2)*clhs193 + clhs266;
    lhs(6,9)=DN(1,0)*clhs87 + DN(1,1)*clhs156 + DN(1,2)*clhs195 + clhs267;
    lhs(6,10)=DN(1,0)*clhs93 + DN(1,1)*clhs161 + DN(1,2)*clhs197 + clhs235 + clhs269;
    lhs(6,11)=DN(2,2)*clhs222 - clhs270 - clhs271*clhs67;
    lhs(6,12)=DN(1,0)*clhs101 + DN(1,1)*clhs166 + DN(1,2)*clhs202 + clhs272;
    lhs(6,13)=DN(1,0)*clhs116 + DN(1,1)*clhs170 + DN(1,2)*clhs204 + clhs273;
    lhs(6,14)=DN(1,0)*clhs122 + DN(1,1)*clhs175 + DN(1,2)*clhs206 + clhs244 + clhs275;
    lhs(6,15)=DN(3,2)*clhs222 - clhs276 - clhs277*clhs67;
    lhs(7,0)=clhs16*clhs278 + clhs65;
    lhs(7,1)=clhs149 + clhs16*clhs279;
    lhs(7,2)=clhs16*clhs280 + clhs191;
    lhs(7,3)=clhs216;
    lhs(7,4)=DN(1,0)*clhs281;
    lhs(7,5)=DN(1,1)*clhs281;
    lhs(7,6)=DN(1,2)*clhs281;
    lhs(7,7)=clhs212*(clhs224 + clhs249 + clhs264);
    lhs(7,8)=clhs239 + clhs278*clhs80;
    lhs(7,9)=clhs257 + clhs279*clhs80;
    lhs(7,10)=clhs271 + clhs280*clhs80;
    lhs(7,11)=clhs282;
    lhs(7,12)=clhs109*clhs278 + clhs248;
    lhs(7,13)=clhs109*clhs279 + clhs263;
    lhs(7,14)=clhs109*clhs280 + clhs277;
    lhs(7,15)=clhs283;
    lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs288 + clhs74;
    lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs153;
    lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs194;
    lhs(8,3)=DN(0,0)*clhs287 - clhs67*clhs95 - clhs96;
    lhs(8,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + DN(2,2)*clhs42 + clhs232 + clhs289;
    lhs(8,5)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + DN(2,2)*clhs57 + clhs252;
    lhs(8,6)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + DN(2,2)*clhs63 + clhs266;
    lhs(8,7)=DN(1,0)*clhs287 - clhs238*clhs67 - clhs239;
    lhs(8,8)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs10*clhs290 + clhs292;
    lhs(8,9)=DN(2,0)*clhs82 + DN(2,1)*clhs84 + DN(2,2)*clhs87 + clhs294;
    lhs(8,10)=DN(2,0)*clhs89 + DN(2,1)*clhs91 + DN(2,2)*clhs93 + clhs295;
    lhs(8,11)=DN(2,0)*clhs296;
    lhs(8,12)=DN(2,0)*clhs97 + DN(2,1)*clhs99 + DN(2,2)*clhs101 + clhs298 + clhs301;
    lhs(8,13)=DN(2,0)*clhs111 + DN(2,1)*clhs113 + DN(2,2)*clhs116 + clhs302;
    lhs(8,14)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs303;
    lhs(8,15)=DN(3,0)*clhs287 - clhs304 - clhs305*clhs67;
    lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs126 + DN(2,2)*clhs127 + clhs88;
    lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs128 + DN(2,2)*clhs130 + clhs158 + clhs288;
    lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs196;
    lhs(9,3)=DN(0,1)*clhs287 - clhs163*clhs67 - clhs164;
    lhs(9,4)=DN(2,0)*clhs40 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs236;
    lhs(9,5)=DN(2,0)*clhs54 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs254 + clhs289;
    lhs(9,6)=DN(2,0)*clhs61 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs267;
    lhs(9,7)=DN(1,1)*clhs287 - clhs256*clhs67 - clhs257;
    lhs(9,8)=DN(2,0)*clhs70 + DN(2,1)*clhs151 + DN(2,2)*clhs152 + clhs294;
    lhs(9,9)=DN(2,0)*clhs84 + DN(2,1)*clhs154 + DN(2,2)*clhs156 + clhs10*clhs306 + clhs292;
    lhs(9,10)=DN(2,0)*clhs91 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs308;
    lhs(9,11)=DN(2,1)*clhs296;
    lhs(9,12)=DN(2,0)*clhs99 + DN(2,1)*clhs165 + DN(2,2)*clhs166 + clhs309;
    lhs(9,13)=DN(2,0)*clhs113 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs301 + clhs311;
    lhs(9,14)=DN(2,0)*clhs120 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs312;
    lhs(9,15)=DN(3,1)*clhs287 - clhs313 - clhs314*clhs67;
    lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs179 + clhs94;
    lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs130 + DN(2,2)*clhs180 + clhs162;
    lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs134 + DN(2,2)*clhs181 + clhs199 + clhs288;
    lhs(10,3)=DN(0,2)*clhs287 - clhs200*clhs67 - clhs201;
    lhs(10,4)=DN(2,0)*clhs42 + DN(2,1)*clhs138 + DN(2,2)*clhs183 + clhs237;
    lhs(10,5)=DN(2,0)*clhs57 + DN(2,1)*clhs142 + DN(2,2)*clhs186 + clhs255;
    lhs(10,6)=DN(2,0)*clhs63 + DN(2,1)*clhs147 + DN(2,2)*clhs188 + clhs269 + clhs289;
    lhs(10,7)=DN(1,2)*clhs287 - clhs270*clhs67 - clhs271;
    lhs(10,8)=DN(2,0)*clhs72 + DN(2,1)*clhs152 + DN(2,2)*clhs193 + clhs295;
    lhs(10,9)=DN(2,0)*clhs87 + DN(2,1)*clhs156 + DN(2,2)*clhs195 + clhs308;
    lhs(10,10)=DN(2,0)*clhs93 + DN(2,1)*clhs161 + DN(2,2)*clhs197 + clhs10*clhs315 + clhs292;
    lhs(10,11)=DN(2,2)*clhs296;
    lhs(10,12)=DN(2,0)*clhs101 + DN(2,1)*clhs166 + DN(2,2)*clhs202 + clhs317;
    lhs(10,13)=DN(2,0)*clhs116 + DN(2,1)*clhs170 + DN(2,2)*clhs204 + clhs318;
    lhs(10,14)=DN(2,0)*clhs122 + DN(2,1)*clhs175 + DN(2,2)*clhs206 + clhs301 + clhs320;
    lhs(10,15)=DN(3,2)*clhs287 - clhs321 - clhs322*clhs67;
    lhs(11,0)=clhs16*clhs323 + clhs95;
    lhs(11,1)=clhs16*clhs324 + clhs163;
    lhs(11,2)=clhs16*clhs325 + clhs200;
    lhs(11,3)=clhs217;
    lhs(11,4)=clhs238 + clhs323*clhs50;
    lhs(11,5)=clhs256 + clhs324*clhs50;
    lhs(11,6)=clhs270 + clhs325*clhs50;
    lhs(11,7)=clhs282;
    lhs(11,8)=DN(2,0)*clhs326;
    lhs(11,9)=DN(2,1)*clhs326;
    lhs(11,10)=DN(2,2)*clhs326;
    lhs(11,11)=clhs212*(clhs290 + clhs306 + clhs315);
    lhs(11,12)=clhs109*clhs323 + clhs305;
    lhs(11,13)=clhs109*clhs324 + clhs314;
    lhs(11,14)=clhs109*clhs325 + clhs322;
    lhs(11,15)=clhs327;
    lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs103 + clhs332;
    lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs167;
    lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs203;
    lhs(12,3)=DN(0,0)*clhs331 - clhs124*clhs67 - clhs125;
    lhs(12,4)=DN(3,0)*clhs38 + DN(3,1)*clhs40 + DN(3,2)*clhs42 + clhs241 + clhs333;
    lhs(12,5)=DN(3,0)*clhs52 + DN(3,1)*clhs54 + DN(3,2)*clhs57 + clhs258;
    lhs(12,6)=DN(3,0)*clhs59 + DN(3,1)*clhs61 + DN(3,2)*clhs63 + clhs272;
    lhs(12,7)=DN(1,0)*clhs331 - clhs247*clhs67 - clhs248;
    lhs(12,8)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs298 + clhs334;
    lhs(12,9)=DN(3,0)*clhs82 + DN(3,1)*clhs84 + DN(3,2)*clhs87 + clhs309;
    lhs(12,10)=DN(3,0)*clhs89 + DN(3,1)*clhs91 + DN(3,2)*clhs93 + clhs317;
    lhs(12,11)=DN(2,0)*clhs331 - clhs304*clhs67 - clhs305;
    lhs(12,12)=DN(3,0)*clhs97 + DN(3,1)*clhs99 + DN(3,2)*clhs101 + clhs10*clhs335 + clhs337;
    lhs(12,13)=DN(3,0)*clhs111 + DN(3,1)*clhs113 + DN(3,2)*clhs116 + clhs339;
    lhs(12,14)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs340;
    lhs(12,15)=DN(3,0)*clhs341;
    lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs126 + DN(3,2)*clhs127 + clhs117;
    lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs128 + DN(3,2)*clhs130 + clhs172 + clhs332;
    lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs205;
    lhs(13,3)=DN(0,1)*clhs331 - clhs177*clhs67 - clhs178;
    lhs(13,4)=DN(3,0)*clhs40 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs245;
    lhs(13,5)=DN(3,0)*clhs54 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs260 + clhs333;
    lhs(13,6)=DN(3,0)*clhs61 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs273;
    lhs(13,7)=DN(1,1)*clhs331 - clhs262*clhs67 - clhs263;
    lhs(13,8)=DN(3,0)*clhs70 + DN(3,1)*clhs151 + DN(3,2)*clhs152 + clhs302;
    lhs(13,9)=DN(3,0)*clhs84 + DN(3,1)*clhs154 + DN(3,2)*clhs156 + clhs311 + clhs334;
    lhs(13,10)=DN(3,0)*clhs91 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs318;
    lhs(13,11)=DN(2,1)*clhs331 - clhs313*clhs67 - clhs314;
    lhs(13,12)=DN(3,0)*clhs99 + DN(3,1)*clhs165 + DN(3,2)*clhs166 + clhs339;
    lhs(13,13)=DN(3,0)*clhs113 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs10*clhs342 + clhs337;
    lhs(13,14)=DN(3,0)*clhs120 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs343;
    lhs(13,15)=DN(3,1)*clhs341;
    lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs179 + clhs123;
    lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs130 + DN(3,2)*clhs180 + clhs176;
    lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs134 + DN(3,2)*clhs181 + clhs208 + clhs332;
    lhs(14,3)=DN(0,2)*clhs331 - clhs209*clhs67 - clhs210;
    lhs(14,4)=DN(3,0)*clhs42 + DN(3,1)*clhs138 + DN(3,2)*clhs183 + clhs246;
    lhs(14,5)=DN(3,0)*clhs57 + DN(3,1)*clhs142 + DN(3,2)*clhs186 + clhs261;
    lhs(14,6)=DN(3,0)*clhs63 + DN(3,1)*clhs147 + DN(3,2)*clhs188 + clhs275 + clhs333;
    lhs(14,7)=DN(1,2)*clhs331 - clhs276*clhs67 - clhs277;
    lhs(14,8)=DN(3,0)*clhs72 + DN(3,1)*clhs152 + DN(3,2)*clhs193 + clhs303;
    lhs(14,9)=DN(3,0)*clhs87 + DN(3,1)*clhs156 + DN(3,2)*clhs195 + clhs312;
    lhs(14,10)=DN(3,0)*clhs93 + DN(3,1)*clhs161 + DN(3,2)*clhs197 + clhs320 + clhs334;
    lhs(14,11)=DN(2,2)*clhs331 - clhs321*clhs67 - clhs322;
    lhs(14,12)=DN(3,0)*clhs101 + DN(3,1)*clhs166 + DN(3,2)*clhs202 + clhs340;
    lhs(14,13)=DN(3,0)*clhs116 + DN(3,1)*clhs170 + DN(3,2)*clhs204 + clhs343;
    lhs(14,14)=DN(3,0)*clhs122 + DN(3,1)*clhs175 + DN(3,2)*clhs206 + clhs10*clhs344 + clhs337;
    lhs(14,15)=DN(3,2)*clhs341;
    lhs(15,0)=clhs124 + clhs16*clhs345;
    lhs(15,1)=clhs16*clhs346 + clhs177;
    lhs(15,2)=clhs16*clhs347 + clhs209;
    lhs(15,3)=clhs218;
    lhs(15,4)=clhs247 + clhs345*clhs50;
    lhs(15,5)=clhs262 + clhs346*clhs50;
    lhs(15,6)=clhs276 + clhs347*clhs50;
    lhs(15,7)=clhs283;
    lhs(15,8)=clhs304 + clhs345*clhs80;
    lhs(15,9)=clhs313 + clhs346*clhs80;
    lhs(15,10)=clhs321 + clhs347*clhs80;
    lhs(15,11)=clhs327;
    lhs(15,12)=DN(3,0)*clhs348;
    lhs(15,13)=DN(3,1)*clhs348;
    lhs(15,14)=DN(3,2)*clhs348;
    lhs(15,15)=clhs212*(clhs335 + clhs342 + clhs344);


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

    double dyn_tau = rData.DynamicTau;
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
    double stab_c1 = 4.0;
    double stab_c2 = 2.0;

    /* unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
    const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
    const double crhs2 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
    const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
    const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
    const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double crhs7 =             rho*(crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
    const double crhs8 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double crhs9 =             crhs4 + crhs8;
    const double crhs10 =             rho*stab_c2*sqrt(pow(crhs5, 2) + pow(crhs6, 2));
    const double crhs11 =             crhs9*(crhs10*h/stab_c1 + mu);
    const double crhs12 =             K_darcy*N[0];
    const double crhs13 =             1.0/(K_darcy + crhs10/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crhs14 =             1.0*crhs13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs3 + crhs7);
    const double crhs15 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6);
    const double crhs16 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
    const double crhs17 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
    const double crhs18 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
    const double crhs19 =             rho*(crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs8);
    const double crhs20 =             1.0*crhs13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs16 + crhs17 + crhs18 + crhs19);
    const double crhs21 =             K_darcy*N[1];
    const double crhs22 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6);
    const double crhs23 =             K_darcy*N[2];
    const double crhs24 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6);
    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs7 + crhs12*crhs14 - crhs14*crhs15;
    rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs11 - DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 - N[0]*crhs18 - N[0]*crhs19 + crhs12*crhs20 - crhs15*crhs20;
    rhs[2]=-DN(0,0)*crhs14 - DN(0,1)*crhs20 - N[0]*crhs9;
    rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs7 + crhs14*crhs21 - crhs14*crhs22;
    rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs11 - DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 - N[1]*crhs18 - N[1]*crhs19 + crhs20*crhs21 - crhs20*crhs22;
    rhs[5]=-DN(1,0)*crhs14 - DN(1,1)*crhs20 - N[1]*crhs9;
    rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs7 + crhs14*crhs23 - crhs14*crhs24;
    rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs11 - DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 - N[2]*crhs18 - N[2]*crhs19 + crhs20*crhs23 - crhs20*crhs24;
    rhs[8]=-DN(2,0)*crhs14 - DN(2,1)*crhs20 - N[2]*crhs9;


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

    double dyn_tau = rData.DynamicTau;
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
    double stab_c1 = 4.0;
    double stab_c2 = 2.0;

 /*    unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &rhs = rData.rhs;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
    const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
    const double crhs2 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
    const double crhs3 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
    const double crhs4 =             DN(0,0)*v(0,0);
    const double crhs5 =             DN(1,0)*v(1,0);
    const double crhs6 =             DN(2,0)*v(2,0);
    const double crhs7 =             DN(3,0)*v(3,0);
    const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double crhs10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double crhs11 =             rho*(crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*(crhs4 + crhs5 + crhs6 + crhs7) + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
    const double crhs12 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    const double crhs13 =             DN(0,1)*v(0,1);
    const double crhs14 =             DN(1,1)*v(1,1);
    const double crhs15 =             DN(2,1)*v(2,1);
    const double crhs16 =             DN(3,1)*v(3,1);
    const double crhs17 =             crhs12 + crhs13 + crhs14 + crhs15 + crhs16 + crhs4 + crhs5 + crhs6 + crhs7;
    const double crhs18 =             rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs8, 2) + pow(crhs9, 2));
    const double crhs19 =             crhs17*(crhs18*h/stab_c1 + mu);
    const double crhs20 =             K_darcy*N[0];
    const double crhs21 =             1.0/(K_darcy + crhs18/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crhs22 =             1.0*crhs21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs11 + crhs2 + crhs3);
    const double crhs23 =             rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10);
    const double crhs24 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
    const double crhs25 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
    const double crhs26 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
    const double crhs27 =             rho*(crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs9*(crhs13 + crhs14 + crhs15 + crhs16));
    const double crhs28 =             1.0*crhs21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs24 + crhs25 + crhs26 + crhs27);
    const double crhs29 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
    const double crhs30 =             K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
    const double crhs31 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
    const double crhs32 =             rho*(crhs10*crhs12 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
    const double crhs33 =             1.0*crhs21*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs29 + crhs30 + crhs31 + crhs32);
    const double crhs34 =             K_darcy*N[1];
    const double crhs35 =             rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10);
    const double crhs36 =             K_darcy*N[2];
    const double crhs37 =             rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10);
    const double crhs38 =             K_darcy*N[3];
    const double crhs39 =             rho*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10);
    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs19 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs11 - N[0]*crhs2 - N[0]*crhs3 + crhs20*crhs22 - crhs22*crhs23;
    rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs19 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs27 + crhs20*crhs28 - crhs23*crhs28;
    rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs19 - DN(0,2)*stress[2] + N[0]*crhs29 - N[0]*crhs30 - N[0]*crhs31 - N[0]*crhs32 + crhs20*crhs33 - crhs23*crhs33;
    rhs[3]=-DN(0,0)*crhs22 - DN(0,1)*crhs28 - DN(0,2)*crhs33 - N[0]*crhs17;
    rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs19 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs11 - N[1]*crhs2 - N[1]*crhs3 + crhs22*crhs34 - crhs22*crhs35;
    rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs19 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs27 + crhs28*crhs34 - crhs28*crhs35;
    rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs19 - DN(1,2)*stress[2] + N[1]*crhs29 - N[1]*crhs30 - N[1]*crhs31 - N[1]*crhs32 + crhs33*crhs34 - crhs33*crhs35;
    rhs[7]=-DN(1,0)*crhs22 - DN(1,1)*crhs28 - DN(1,2)*crhs33 - N[1]*crhs17;
    rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs19 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs11 - N[2]*crhs2 - N[2]*crhs3 + crhs22*crhs36 - crhs22*crhs37;
    rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs19 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs27 + crhs28*crhs36 - crhs28*crhs37;
    rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs19 - DN(2,2)*stress[2] + N[2]*crhs29 - N[2]*crhs30 - N[2]*crhs31 - N[2]*crhs32 + crhs33*crhs36 - crhs33*crhs37;
    rhs[11]=-DN(2,0)*crhs22 - DN(2,1)*crhs28 - DN(2,2)*crhs33 - N[2]*crhs17;
    rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs19 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs11 - N[3]*crhs2 - N[3]*crhs3 + crhs22*crhs38 - crhs22*crhs39;
    rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs19 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs27 + crhs28*crhs38 - crhs28*crhs39;
    rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs19 - DN(3,2)*stress[2] + N[3]*crhs29 - N[3]*crhs30 - N[3]*crhs31 - N[3]*crhs32 + crhs33*crhs38 - crhs33*crhs39;
    rhs[15]=-DN(3,0)*crhs22 - DN(3,1)*crhs28 - DN(3,2)*crhs33 - N[3]*crhs17;


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

    double dyn_tau = rData.DynamicTau;
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
    double stab_c1 = 4.0;
    double stab_c2 = 2.0;

    /* unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double cV2 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double cV3 =             1.0*DNenr(0,0)*K_darcy*cV2;
    const double cV4 =             DN(0,0)*cV0 + DN(0,1)*cV1;
    const double cV5 =             1.0*DNenr(0,0)*cV2*rho;
    const double cV6 =             1.0*DNenr(1,0)*K_darcy*cV2;
    const double cV7 =             1.0*DNenr(1,0)*cV2*rho;
    const double cV8 =             1.0*DNenr(2,0)*K_darcy*cV2;
    const double cV9 =             1.0*DNenr(2,0)*cV2*rho;
    const double cV10 =             1.0*DNenr(0,1)*K_darcy*cV2;
    const double cV11 =             1.0*DNenr(0,1)*cV2*rho;
    const double cV12 =             1.0*DNenr(1,1)*K_darcy*cV2;
    const double cV13 =             1.0*DNenr(1,1)*cV2*rho;
    const double cV14 =             1.0*DNenr(2,1)*K_darcy*cV2;
    const double cV15 =             1.0*DNenr(2,1)*cV2*rho;
    const double cV16 =             1.0*cV2;
    const double cV17 =             DN(1,0)*cV0 + DN(1,1)*cV1;
    const double cV18 =             DN(2,0)*cV0 + DN(2,1)*cV1;
    V(0,0)=-DN(0,0)*Nenr[0] - N[0]*cV3 + cV4*cV5;
    V(0,1)=-DN(0,0)*Nenr[1] - N[0]*cV6 + cV4*cV7;
    V(0,2)=-DN(0,0)*Nenr[2] - N[0]*cV8 + cV4*cV9;
    V(1,0)=-DN(0,1)*Nenr[0] - N[0]*cV10 + cV11*cV4;
    V(1,1)=-DN(0,1)*Nenr[1] - N[0]*cV12 + cV13*cV4;
    V(1,2)=-DN(0,1)*Nenr[2] - N[0]*cV14 + cV15*cV4;
    V(2,0)=cV16*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
    V(2,1)=cV16*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
    V(2,2)=cV16*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
    V(3,0)=-DN(1,0)*Nenr[0] - N[1]*cV3 + cV17*cV5;
    V(3,1)=-DN(1,0)*Nenr[1] - N[1]*cV6 + cV17*cV7;
    V(3,2)=-DN(1,0)*Nenr[2] - N[1]*cV8 + cV17*cV9;
    V(4,0)=-DN(1,1)*Nenr[0] - N[1]*cV10 + cV11*cV17;
    V(4,1)=-DN(1,1)*Nenr[1] - N[1]*cV12 + cV13*cV17;
    V(4,2)=-DN(1,1)*Nenr[2] - N[1]*cV14 + cV15*cV17;
    V(5,0)=cV16*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
    V(5,1)=cV16*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
    V(5,2)=cV16*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
    V(6,0)=-DN(2,0)*Nenr[0] - N[2]*cV3 + cV18*cV5;
    V(6,1)=-DN(2,0)*Nenr[1] - N[2]*cV6 + cV18*cV7;
    V(6,2)=-DN(2,0)*Nenr[2] - N[2]*cV8 + cV18*cV9;
    V(7,0)=-DN(2,1)*Nenr[0] - N[2]*cV10 + cV11*cV18;
    V(7,1)=-DN(2,1)*Nenr[1] - N[2]*cV12 + cV13*cV18;
    V(7,2)=-DN(2,1)*Nenr[2] - N[2]*cV14 + cV15*cV18;
    V(8,0)=cV16*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
    V(8,1)=cV16*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
    V(8,2)=cV16*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double cH2 =             K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0);
    const double cH3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double cH4 =             1.0*DNenr(0,0)*cH3;
    const double cH5 =             1.0*DNenr(0,1)*cH3;
    const double cH6 =             1.0*cH3;
    const double cH7 =             K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0);
    const double cH8 =             K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0);
    const double cH9 =             1.0*DNenr(1,0)*cH3;
    const double cH10 =             1.0*DNenr(1,1)*cH3;
    const double cH11 =             1.0*DNenr(2,0)*cH3;
    const double cH12 =             1.0*DNenr(2,1)*cH3;
    H(0,0)=DN(0,0)*Nenr[0] + cH2*cH4;
    H(0,1)=DN(0,1)*Nenr[0] + cH2*cH5;
    H(0,2)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
    H(0,3)=DN(1,0)*Nenr[0] + cH4*cH7;
    H(0,4)=DN(1,1)*Nenr[0] + cH5*cH7;
    H(0,5)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
    H(0,6)=DN(2,0)*Nenr[0] + cH4*cH8;
    H(0,7)=DN(2,1)*Nenr[0] + cH5*cH8;
    H(0,8)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
    H(1,0)=DN(0,0)*Nenr[1] + cH2*cH9;
    H(1,1)=DN(0,1)*Nenr[1] + cH10*cH2;
    H(1,2)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
    H(1,3)=DN(1,0)*Nenr[1] + cH7*cH9;
    H(1,4)=DN(1,1)*Nenr[1] + cH10*cH7;
    H(1,5)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
    H(1,6)=DN(2,0)*Nenr[1] + cH8*cH9;
    H(1,7)=DN(2,1)*Nenr[1] + cH10*cH8;
    H(1,8)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
    H(2,0)=DN(0,0)*Nenr[2] + cH11*cH2;
    H(2,1)=DN(0,1)*Nenr[2] + cH12*cH2;
    H(2,2)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
    H(2,3)=DN(1,0)*Nenr[2] + cH11*cH7;
    H(2,4)=DN(1,1)*Nenr[2] + cH12*cH7;
    H(2,5)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
    H(2,6)=DN(2,0)*Nenr[2] + cH11*cH8;
    H(2,7)=DN(2,1)*Nenr[2] + cH12*cH8;
    H(2,8)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
    const double crhs_ee2 =             crhs_ee0 + crhs_ee1;
    const double crhs_ee3 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double crhs_ee4 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double crhs_ee5 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crhs_ee6 =             1.0*crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
    const double crhs_ee7 =             1.0*crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
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

    double dyn_tau = rData.DynamicTau;
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
    double stab_c1 = 4.0;
    double stab_c2 = 2.0;

/*     unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i<4; ++i)
        if(rData.Distance[i] >= 0) npos += 1;
        else nneg += 1;
    
    if(nneg!=0 && npos!=0)
    {
        stab_c1 *= 1000.0;
        stab_c2 *= 1000.0;
        dyn_tau *= 1000.0;
    } */

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double cV2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double cV3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cV0, 2) + pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double cV4 =             1.0*DNenr(0,0)*K_darcy*cV3;
    const double cV5 =             DN(0,0)*cV0 + DN(0,1)*cV1 + DN(0,2)*cV2;
    const double cV6 =             1.0*DNenr(0,0)*cV3*rho;
    const double cV7 =             1.0*DNenr(1,0)*K_darcy*cV3;
    const double cV8 =             1.0*DNenr(1,0)*cV3*rho;
    const double cV9 =             1.0*DNenr(2,0)*K_darcy*cV3;
    const double cV10 =             1.0*DNenr(2,0)*cV3*rho;
    const double cV11 =             1.0*DNenr(3,0)*K_darcy*cV3;
    const double cV12 =             1.0*DNenr(3,0)*cV3*rho;
    const double cV13 =             1.0*DNenr(0,1)*K_darcy*cV3;
    const double cV14 =             1.0*DNenr(0,1)*cV3*rho;
    const double cV15 =             1.0*DNenr(1,1)*K_darcy*cV3;
    const double cV16 =             1.0*DNenr(1,1)*cV3*rho;
    const double cV17 =             1.0*DNenr(2,1)*K_darcy*cV3;
    const double cV18 =             1.0*DNenr(2,1)*cV3*rho;
    const double cV19 =             1.0*DNenr(3,1)*K_darcy*cV3;
    const double cV20 =             1.0*DNenr(3,1)*cV3*rho;
    const double cV21 =             1.0*DNenr(0,2)*K_darcy*cV3;
    const double cV22 =             1.0*DNenr(0,2)*cV3*rho;
    const double cV23 =             1.0*DNenr(1,2)*K_darcy*cV3;
    const double cV24 =             1.0*DNenr(1,2)*cV3*rho;
    const double cV25 =             1.0*DNenr(2,2)*K_darcy*cV3;
    const double cV26 =             1.0*DNenr(2,2)*cV3*rho;
    const double cV27 =             1.0*DNenr(3,2)*K_darcy*cV3;
    const double cV28 =             1.0*DNenr(3,2)*cV3*rho;
    const double cV29 =             1.0*cV3;
    const double cV30 =             DN(1,0)*cV0 + DN(1,1)*cV1 + DN(1,2)*cV2;
    const double cV31 =             DN(2,0)*cV0 + DN(2,1)*cV1 + DN(2,2)*cV2;
    const double cV32 =             DN(3,0)*cV0 + DN(3,1)*cV1 + DN(3,2)*cV2;
    V(0,0)=-DN(0,0)*Nenr[0] - N[0]*cV4 + cV5*cV6;
    V(0,1)=-DN(0,0)*Nenr[1] - N[0]*cV7 + cV5*cV8;
    V(0,2)=-DN(0,0)*Nenr[2] - N[0]*cV9 + cV10*cV5;
    V(0,3)=-DN(0,0)*Nenr[3] - N[0]*cV11 + cV12*cV5;
    V(1,0)=-DN(0,1)*Nenr[0] - N[0]*cV13 + cV14*cV5;
    V(1,1)=-DN(0,1)*Nenr[1] - N[0]*cV15 + cV16*cV5;
    V(1,2)=-DN(0,1)*Nenr[2] - N[0]*cV17 + cV18*cV5;
    V(1,3)=-DN(0,1)*Nenr[3] - N[0]*cV19 + cV20*cV5;
    V(2,0)=-DN(0,2)*Nenr[0] - N[0]*cV21 + cV22*cV5;
    V(2,1)=-DN(0,2)*Nenr[1] - N[0]*cV23 + cV24*cV5;
    V(2,2)=-DN(0,2)*Nenr[2] - N[0]*cV25 + cV26*cV5;
    V(2,3)=-DN(0,2)*Nenr[3] - N[0]*cV27 + cV28*cV5;
    V(3,0)=cV29*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
    V(3,1)=cV29*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
    V(3,2)=cV29*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
    V(3,3)=cV29*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
    V(4,0)=-DN(1,0)*Nenr[0] - N[1]*cV4 + cV30*cV6;
    V(4,1)=-DN(1,0)*Nenr[1] - N[1]*cV7 + cV30*cV8;
    V(4,2)=-DN(1,0)*Nenr[2] - N[1]*cV9 + cV10*cV30;
    V(4,3)=-DN(1,0)*Nenr[3] - N[1]*cV11 + cV12*cV30;
    V(5,0)=-DN(1,1)*Nenr[0] - N[1]*cV13 + cV14*cV30;
    V(5,1)=-DN(1,1)*Nenr[1] - N[1]*cV15 + cV16*cV30;
    V(5,2)=-DN(1,1)*Nenr[2] - N[1]*cV17 + cV18*cV30;
    V(5,3)=-DN(1,1)*Nenr[3] - N[1]*cV19 + cV20*cV30;
    V(6,0)=-DN(1,2)*Nenr[0] - N[1]*cV21 + cV22*cV30;
    V(6,1)=-DN(1,2)*Nenr[1] - N[1]*cV23 + cV24*cV30;
    V(6,2)=-DN(1,2)*Nenr[2] - N[1]*cV25 + cV26*cV30;
    V(6,3)=-DN(1,2)*Nenr[3] - N[1]*cV27 + cV28*cV30;
    V(7,0)=cV29*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
    V(7,1)=cV29*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
    V(7,2)=cV29*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
    V(7,3)=cV29*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
    V(8,0)=-DN(2,0)*Nenr[0] - N[2]*cV4 + cV31*cV6;
    V(8,1)=-DN(2,0)*Nenr[1] - N[2]*cV7 + cV31*cV8;
    V(8,2)=-DN(2,0)*Nenr[2] - N[2]*cV9 + cV10*cV31;
    V(8,3)=-DN(2,0)*Nenr[3] - N[2]*cV11 + cV12*cV31;
    V(9,0)=-DN(2,1)*Nenr[0] - N[2]*cV13 + cV14*cV31;
    V(9,1)=-DN(2,1)*Nenr[1] - N[2]*cV15 + cV16*cV31;
    V(9,2)=-DN(2,1)*Nenr[2] - N[2]*cV17 + cV18*cV31;
    V(9,3)=-DN(2,1)*Nenr[3] - N[2]*cV19 + cV20*cV31;
    V(10,0)=-DN(2,2)*Nenr[0] - N[2]*cV21 + cV22*cV31;
    V(10,1)=-DN(2,2)*Nenr[1] - N[2]*cV23 + cV24*cV31;
    V(10,2)=-DN(2,2)*Nenr[2] - N[2]*cV25 + cV26*cV31;
    V(10,3)=-DN(2,2)*Nenr[3] - N[2]*cV27 + cV28*cV31;
    V(11,0)=cV29*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
    V(11,1)=cV29*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
    V(11,2)=cV29*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
    V(11,3)=cV29*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
    V(12,0)=-DN(3,0)*Nenr[0] - N[3]*cV4 + cV32*cV6;
    V(12,1)=-DN(3,0)*Nenr[1] - N[3]*cV7 + cV32*cV8;
    V(12,2)=-DN(3,0)*Nenr[2] - N[3]*cV9 + cV10*cV32;
    V(12,3)=-DN(3,0)*Nenr[3] - N[3]*cV11 + cV12*cV32;
    V(13,0)=-DN(3,1)*Nenr[0] - N[3]*cV13 + cV14*cV32;
    V(13,1)=-DN(3,1)*Nenr[1] - N[3]*cV15 + cV16*cV32;
    V(13,2)=-DN(3,1)*Nenr[2] - N[3]*cV17 + cV18*cV32;
    V(13,3)=-DN(3,1)*Nenr[3] - N[3]*cV19 + cV20*cV32;
    V(14,0)=-DN(3,2)*Nenr[0] - N[3]*cV21 + cV22*cV32;
    V(14,1)=-DN(3,2)*Nenr[1] - N[3]*cV23 + cV24*cV32;
    V(14,2)=-DN(3,2)*Nenr[2] - N[3]*cV25 + cV26*cV32;
    V(14,3)=-DN(3,2)*Nenr[3] - N[3]*cV27 + cV28*cV32;
    V(15,0)=cV29*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
    V(15,1)=cV29*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
    V(15,2)=cV29*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
    V(15,3)=cV29*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double cH3 =             K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0);
    const double cH4 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double cH5 =             1.0*DNenr(0,0)*cH4;
    const double cH6 =             1.0*DNenr(0,1)*cH4;
    const double cH7 =             1.0*DNenr(0,2)*cH4;
    const double cH8 =             1.0*cH4;
    const double cH9 =             K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0);
    const double cH10 =             K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0);
    const double cH11 =             K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0);
    const double cH12 =             1.0*DNenr(1,0)*cH4;
    const double cH13 =             1.0*DNenr(1,1)*cH4;
    const double cH14 =             1.0*DNenr(1,2)*cH4;
    const double cH15 =             1.0*DNenr(2,0)*cH4;
    const double cH16 =             1.0*DNenr(2,1)*cH4;
    const double cH17 =             1.0*DNenr(2,2)*cH4;
    const double cH18 =             1.0*DNenr(3,0)*cH4;
    const double cH19 =             1.0*DNenr(3,1)*cH4;
    const double cH20 =             1.0*DNenr(3,2)*cH4;
    H(0,0)=DN(0,0)*Nenr[0] + cH3*cH5;
    H(0,1)=DN(0,1)*Nenr[0] + cH3*cH6;
    H(0,2)=DN(0,2)*Nenr[0] + cH3*cH7;
    H(0,3)=cH8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
    H(0,4)=DN(1,0)*Nenr[0] + cH5*cH9;
    H(0,5)=DN(1,1)*Nenr[0] + cH6*cH9;
    H(0,6)=DN(1,2)*Nenr[0] + cH7*cH9;
    H(0,7)=cH8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
    H(0,8)=DN(2,0)*Nenr[0] + cH10*cH5;
    H(0,9)=DN(2,1)*Nenr[0] + cH10*cH6;
    H(0,10)=DN(2,2)*Nenr[0] + cH10*cH7;
    H(0,11)=cH8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
    H(0,12)=DN(3,0)*Nenr[0] + cH11*cH5;
    H(0,13)=DN(3,1)*Nenr[0] + cH11*cH6;
    H(0,14)=DN(3,2)*Nenr[0] + cH11*cH7;
    H(0,15)=cH8*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
    H(1,0)=DN(0,0)*Nenr[1] + cH12*cH3;
    H(1,1)=DN(0,1)*Nenr[1] + cH13*cH3;
    H(1,2)=DN(0,2)*Nenr[1] + cH14*cH3;
    H(1,3)=cH8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
    H(1,4)=DN(1,0)*Nenr[1] + cH12*cH9;
    H(1,5)=DN(1,1)*Nenr[1] + cH13*cH9;
    H(1,6)=DN(1,2)*Nenr[1] + cH14*cH9;
    H(1,7)=cH8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
    H(1,8)=DN(2,0)*Nenr[1] + cH10*cH12;
    H(1,9)=DN(2,1)*Nenr[1] + cH10*cH13;
    H(1,10)=DN(2,2)*Nenr[1] + cH10*cH14;
    H(1,11)=cH8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
    H(1,12)=DN(3,0)*Nenr[1] + cH11*cH12;
    H(1,13)=DN(3,1)*Nenr[1] + cH11*cH13;
    H(1,14)=DN(3,2)*Nenr[1] + cH11*cH14;
    H(1,15)=cH8*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
    H(2,0)=DN(0,0)*Nenr[2] + cH15*cH3;
    H(2,1)=DN(0,1)*Nenr[2] + cH16*cH3;
    H(2,2)=DN(0,2)*Nenr[2] + cH17*cH3;
    H(2,3)=cH8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
    H(2,4)=DN(1,0)*Nenr[2] + cH15*cH9;
    H(2,5)=DN(1,1)*Nenr[2] + cH16*cH9;
    H(2,6)=DN(1,2)*Nenr[2] + cH17*cH9;
    H(2,7)=cH8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
    H(2,8)=DN(2,0)*Nenr[2] + cH10*cH15;
    H(2,9)=DN(2,1)*Nenr[2] + cH10*cH16;
    H(2,10)=DN(2,2)*Nenr[2] + cH10*cH17;
    H(2,11)=cH8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
    H(2,12)=DN(3,0)*Nenr[2] + cH11*cH15;
    H(2,13)=DN(3,1)*Nenr[2] + cH11*cH16;
    H(2,14)=DN(3,2)*Nenr[2] + cH11*cH17;
    H(2,15)=cH8*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
    H(3,0)=DN(0,0)*Nenr[3] + cH18*cH3;
    H(3,1)=DN(0,1)*Nenr[3] + cH19*cH3;
    H(3,2)=DN(0,2)*Nenr[3] + cH20*cH3;
    H(3,3)=cH8*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
    H(3,4)=DN(1,0)*Nenr[3] + cH18*cH9;
    H(3,5)=DN(1,1)*Nenr[3] + cH19*cH9;
    H(3,6)=DN(1,2)*Nenr[3] + cH20*cH9;
    H(3,7)=cH8*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
    H(3,8)=DN(2,0)*Nenr[3] + cH10*cH18;
    H(3,9)=DN(2,1)*Nenr[3] + cH10*cH19;
    H(3,10)=DN(2,2)*Nenr[3] + cH10*cH20;
    H(3,11)=cH8*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
    H(3,12)=DN(3,0)*Nenr[3] + cH11*cH18;
    H(3,13)=DN(3,1)*Nenr[3] + cH11*cH19;
    H(3,14)=DN(3,2)*Nenr[3] + cH11*cH20;
    H(3,15)=cH8*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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


    const double crhs_ee0 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    const double crhs_ee1 =             DN(0,0)*v(0,0);
    const double crhs_ee2 =             DN(0,1)*v(0,1);
    const double crhs_ee3 =             DN(1,0)*v(1,0);
    const double crhs_ee4 =             DN(1,1)*v(1,1);
    const double crhs_ee5 =             DN(2,0)*v(2,0);
    const double crhs_ee6 =             DN(2,1)*v(2,1);
    const double crhs_ee7 =             DN(3,0)*v(3,0);
    const double crhs_ee8 =             DN(3,1)*v(3,1);
    const double crhs_ee9 =             crhs_ee0 + crhs_ee1 + crhs_ee2 + crhs_ee3 + crhs_ee4 + crhs_ee5 + crhs_ee6 + crhs_ee7 + crhs_ee8;
    const double crhs_ee10 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double crhs_ee11 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double crhs_ee12 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double crhs_ee13 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee10, 2) + pow(crhs_ee11, 2) + pow(crhs_ee12, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
    const double crhs_ee14 =             1.0*crhs_ee13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee10*(crhs_ee1 + crhs_ee3 + crhs_ee5 + crhs_ee7) + crhs_ee11*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee12*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
    const double crhs_ee15 =             1.0*crhs_ee13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee10*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee11*(crhs_ee2 + crhs_ee4 + crhs_ee6 + crhs_ee8) + crhs_ee12*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
    const double crhs_ee16 =             1.0*crhs_ee13*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee0*crhs_ee12 + crhs_ee10*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee11*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
    rhs_ee[0]=-DNenr(0,0)*crhs_ee14 - DNenr(0,1)*crhs_ee15 - DNenr(0,2)*crhs_ee16 - Nenr[0]*crhs_ee9;
    rhs_ee[1]=-DNenr(1,0)*crhs_ee14 - DNenr(1,1)*crhs_ee15 - DNenr(1,2)*crhs_ee16 - Nenr[1]*crhs_ee9;
    rhs_ee[2]=-DNenr(2,0)*crhs_ee14 - DNenr(2,1)*crhs_ee15 - DNenr(2,2)*crhs_ee16 - Nenr[2]*crhs_ee9;
    rhs_ee[3]=-DNenr(3,0)*crhs_ee14 - DNenr(3,1)*crhs_ee15 - DNenr(3,2)*crhs_ee16 - Nenr[3]*crhs_ee9;


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
    GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg)
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

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else /* if (rData.Distance[i] < 0.0) */{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    if (Dim == 2)
        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    else
        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);

    // Call the positive side modified shape functions calculator
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
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

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (p_modified_sh_func->pGetSplittingUtil())->mDivisionsNumber;
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
    MatrixType& rInterfaceShapeFunctionNeg,
    MatrixType& rEnrInterfaceShapeFunctionPos,
    MatrixType& rEnrInterfaceShapeFunctionNeg,
    GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Kratos::Vector& rInterfaceWeightsNeg,
    std::vector<Vector>& rInterfaceNormalsNeg)
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

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else /* if (rData.Distance[i] < 0.0) */{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    if (Dim == 2)
        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    else
        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);

    // Call the positive side modified shape functions calculator
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);

    // Call the Interface negative side shape functions calculator
    p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    rInterfaceShapeFunctionNeg,
    rInterfaceShapeDerivativesNeg,
    rInterfaceWeightsNeg,
    GeometryData::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
    rInterfaceNormalsNeg,
    GeometryData::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); gp++){

        double normal_norm = 0.0;
        for (unsigned int dim = 0; dim < Dim; dim++){
            normal_norm += (rInterfaceNormalsNeg[gp])[dim]*(rInterfaceNormalsNeg[gp])[dim];
        }

        normal_norm = std::sqrt(normal_norm);

        for (unsigned int dim = 0; dim < Dim; dim++){
            (rInterfaceNormalsNeg[gp])[dim] = (rInterfaceNormalsNeg[gp])[dim]/normal_norm;
        }
    }

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (p_modified_sh_func->pGetSplittingUtil())->mDivisionsNumber;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ComputeSplitting(
		TElementData& rData,
		MatrixType& rShapeFunctionsPos,
        MatrixType& rShapeFunctionsNeg,
        MatrixType& rEnrichedShapeFunctionsPos,
        MatrixType& rEnrichedShapeFunctionsNeg,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg,
        MatrixType& rInterfaceShapeFunctionNeg,
        MatrixType& rEnrInterfaceShapeFunctionPos,
        MatrixType& rEnrInterfaceShapeFunctionNeg,
        GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceWeightsNeg,
        std::vector<Vector>& rInterfaceNormalsNeg,
        std::vector<MatrixType>& rContactShapeFunctionNeg,
        std::vector<GeometryType::ShapeFunctionsGradientsType>& rContactShapeDerivativesNeg,
        std::vector<Kratos::Vector>& rContactWeightsNeg,
        std::vector<Vector>& rContactTangentialsNeg)
        //MatrixType& rContactShapeFunctionNeg,
        //GeometryType::ShapeFunctionsGradientsType& rContactShapeDerivativesNeg,
        //Kratos::Vector& rContactWeightsNeg,
        //Vector& rContactTangentialsNeg,
        //bool& rHasContactLine)
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

    for (unsigned int i = 0; i < NumNodes; i++){
        /* if (abs(rData.Distance[i]) < 1.0e-6){
            continue;
        } else */ if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else /* if (rData.Distance[i] < 0.0) */{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    if (Dim == 2)
        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    else
        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);

    // Call the positive side modified shape functions calculator
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);

    // Call the Interface negative side shape functions calculator
    p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rInterfaceShapeFunctionNeg,
        rInterfaceShapeDerivativesNeg,
        rInterfaceWeightsNeg,
        GeometryData::GI_GAUSS_2);

    // Call the Interface negative side normal functions calculator
    p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
        rInterfaceNormalsNeg,
        GeometryData::GI_GAUSS_2);

    for (unsigned int gp = 0; gp < rInterfaceNormalsNeg.size(); gp++){

        double normal_norm = 0.0;
        for (unsigned int dim = 0; dim < Dim; dim++){
            normal_norm += (rInterfaceNormalsNeg[gp])[dim]*(rInterfaceNormalsNeg[gp])[dim];
        }

        normal_norm = std::sqrt(normal_norm);

        for (unsigned int dim = 0; dim < Dim; dim++){
            (rInterfaceNormalsNeg[gp])[dim] = (rInterfaceNormalsNeg[gp])[dim]/normal_norm;
        }
    }

    // Compute the enrichment shape function values using the enrichment interpolation matrices
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, enr_pos_interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, enr_neg_interp);

    // Compute the enrichment shape function values at the interface gauss points using the enrichment interpolation matrices
    rEnrInterfaceShapeFunctionPos = prod(rInterfaceShapeFunctionNeg, enr_pos_interp);
    rEnrInterfaceShapeFunctionNeg = prod(rInterfaceShapeFunctionNeg, enr_neg_interp);

    // Compute the enrichment shape function gradient values using the enrichment interpolation matrices
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++){
        rEnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivativesPos[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++){
        rEnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = (p_modified_sh_func->pGetSplittingUtil())->mDivisionsNumber;

    std::vector<unsigned int> contact_line_faces;
    std::vector<unsigned int> contact_line_indices;

    //KRATOS_INFO("Compute splitting") << "Above ComputeNegativeSideContactLineVector" << std::endl;

    /* rHasContactLine = */ p_modified_sh_func->ComputeNegativeSideContactLineVector(contact_line_faces, rContactTangentialsNeg);
    // rContactTangentialsNeg is normalized in ComputeNegativeSideContactLineVector
    
    //KRATOS_INFO("size of contact_line_faces") << contact_line_faces.size() << std::endl;
    //KRATOS_INFO("size of rContactTangentialsNeg") << rContactTangentialsNeg.size() << std::endl;

    auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    for (unsigned int i_cl = 0; i_cl < contact_line_faces.size(); i_cl++){
        if (neighbour_elems[ contact_line_faces[i_cl] ].Id() == this->Id() ){
            contact_line_indices.push_back(i_cl);
        }
    }

    //KRATOS_INFO("Size of contact_line_indices") << contact_line_indices.size() << std::endl;

    /* double tangent_norm = 0.0;
    for (unsigned int dim = 0; dim < Dim; dim++){
        tangent_norm += rContactTangentialsNeg[dim]*rContactTangentialsNeg[dim];
    }
    tangent_norm = std::sqrt(tangent_norm);
    for (unsigned int dim = 0; dim < Dim; dim++){
        rContactTangentialsNeg[dim] = rContactTangentialsNeg[dim]/tangent_norm;
    } */

    //if (rHasContactLine){ // HERE: HAS CONTACT LINE == contact_line_indices.size() > 0: no need for an explicit check
                          // ELSEWHERE: HAS CONTACT LINE == rContactWeightsNeg.size() > 0
        // Call the Contact Line negative side shape functions calculator
        p_modified_sh_func->ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
            contact_line_indices, //ADDED
            rContactShapeFunctionNeg,
            rContactShapeDerivativesNeg,
            rContactWeightsNeg,
            GeometryData::GI_GAUSS_2);
    //}
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateCurvature(
        const Matrix& rIntShapeFunctions,
        Kratos::Vector& rInterfaceCurvature)
{
    GeometryType::Pointer p_geom = this->pGetGeometry();
    const unsigned int n_nodes = this->NumNodes;
    const unsigned int n_gpt = rIntShapeFunctions.size1();

    rInterfaceCurvature.resize(n_gpt, false);

    for (unsigned int gpt = 0; gpt < n_gpt; gpt++){ 

        double curvature = 0.0;

        for (unsigned int i=0; i < n_nodes; ++i){
            curvature += rIntShapeFunctions(gpt,i) * (*p_geom)[i].FastGetSolutionStepValue(CURVATURE);
        }

        rInterfaceCurvature[gpt] = curvature;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateCurvature(
    const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
    Kratos::Vector& rInterfaceCurvature)
{
    GeometryType::Pointer p_geom = this->pGetGeometry();

    const unsigned int n_nodes = this->NumNodes;
    const unsigned int n_gpt = rInterfaceShapeDerivativesNeg.size();

    rInterfaceCurvature.resize(n_gpt, false);
    
    for (unsigned int gpt = 0; gpt < n_gpt; gpt++){ 

        array_1d<double,n_nodes> DN_Dx, DN_Dy, DN_Dz; 
        array_1d<double,n_nodes> GradPhiX, GradPhiY, GradPhiZ;

        for (unsigned int i=0; i < n_nodes; ++i){
            const double GradPhiXi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_X); 
            const double GradPhiYi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_Y);  
            const double GradPhiZi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_Z);
            const double norm = std::sqrt(GradPhiXi*GradPhiXi + GradPhiYi*GradPhiYi + GradPhiZi*GradPhiZi);
            //KRATOS_INFO("Norm of Grad Phi") << norm << std::endl;
            GradPhiX[i] = GradPhiXi / norm;
            GradPhiY[i] = GradPhiYi / norm;
            GradPhiZ[i] = GradPhiZi / norm;
            DN_Dx[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,0);
            DN_Dy[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,1);
            DN_Dz[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,2);
        }

        const double DGradPhiX_DX = inner_prod(GradPhiX, DN_Dx);
        const double DGradPhiY_DY = inner_prod(GradPhiY, DN_Dy);
        const double DGradPhiY_DZ = inner_prod(GradPhiZ, DN_Dz);

        const double curvature = DGradPhiX_DX + DGradPhiY_DY + DGradPhiY_DZ;

        rInterfaceCurvature[gpt] = curvature;//937.500146484;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateCurvature(
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceCurvature)
{
    GeometryType::Pointer p_geom = this->pGetGeometry();
    const unsigned int n_nodes = this->NumNodes;
    const unsigned int n_gpt = rIntShapeFunctions.size1();

    rInterfaceCurvature.resize(n_gpt, false);

    for (unsigned int gpt = 0; gpt < n_gpt; gpt++){ 

        double curvature_nodal = 0.0;

        array_1d<double,n_nodes> DN_Dx, DN_Dy, DN_Dz; 
        array_1d<double,n_nodes> GradPhiX, GradPhiY, GradPhiZ;

        for (unsigned int i=0; i < n_nodes; ++i){

            curvature_nodal += rIntShapeFunctions(gpt,i) * (*p_geom)[i].FastGetSolutionStepValue(CURVATURE);

            const double GradPhiXi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_X); 
            const double GradPhiYi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_Y);  
            const double GradPhiZi = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT_Z);
            const double norm = std::sqrt(GradPhiXi*GradPhiXi + GradPhiYi*GradPhiYi + GradPhiZi*GradPhiZi);
            //KRATOS_INFO("Norm of Grad Phi") << norm << std::endl;
            GradPhiX[i] = GradPhiXi / norm;
            GradPhiY[i] = GradPhiYi / norm;
            GradPhiZ[i] = GradPhiZi / norm;
            DN_Dx[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,0);
            DN_Dy[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,1);
            DN_Dz[i] = (rInterfaceShapeDerivativesNeg[gpt])(i,2);
        }

        const double DGradPhiX_DX = inner_prod(GradPhiX, DN_Dx);
        const double DGradPhiY_DY = inner_prod(GradPhiY, DN_Dy);
        const double DGradPhiY_DZ = inner_prod(GradPhiZ, DN_Dz);

        const double curvature_elemental = DGradPhiX_DX + DGradPhiY_DY + DGradPhiY_DZ;

        rInterfaceCurvature[gpt] = (2.0*curvature_nodal + 1.0*curvature_elemental)/3.0;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureDiscontinuity(
    const double coefficient,
    const Kratos::Vector& rCurvature,
    const Kratos::Vector& rIntWeights,
    const Matrix& rIntShapeFunctions,
    const Matrix& rIntEnrShapeFunctionsPos,
    const Matrix& rIntEnrShapeFunctionsNeg,
	MatrixType& rKeeTot,
	VectorType& rRHSeeTot)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();

    //KRATOS_INFO("Weights") << rIntWeights << std::endl;

    //A new set of test functions is introduced to retain the symmetry of the modification made to KeeTot
    Matrix TestFunctions = ZeroMatrix(NumGP,NumNodes);

    //Matrix LHS = ZeroMatrix(NumNodes,NumNodes);
    //Vector RHS_Residual(NumNodes);

    for (unsigned int gp = 0; gp < NumGP; gp++)
        for (unsigned int j = 0; j < NumNodes; j++)
            TestFunctions(gp, j) = rIntEnrShapeFunctionsNeg(gp,j) - rIntEnrShapeFunctionsPos(gp,j);

    // Penalty coefficient
    const double penalty_coeff = 2.0;

    for (unsigned int j = 0; j < NumNodes; j++){

        double RHSj = 0.0;

        for (unsigned int gp = 0; gp < NumGP; gp++){
            RHSj += rCurvature(gp)*rIntWeights(gp)*TestFunctions(gp,j);
        }

        rRHSeeTot(j) += coefficient*RHSj*penalty_coeff;

        for (unsigned int i = 0; i < NumNodes; i++){

            double LHSij = 0.0;

            for (unsigned int gp = 0; gp < NumGP; gp++){
                LHSij += (rIntEnrShapeFunctionsNeg(gp,i) - rIntEnrShapeFunctionsPos(gp,i))*rIntWeights(gp)*TestFunctions(gp,j);
            }

            //LHS(i,j) = LHSij*penalty_coeff;
            rKeeTot(i,j) += LHSij*penalty_coeff;
        }
    }

    //Vector pressure_enr_est(NumNodes);

    //double pressure_neg = 0.0;
    //double pressure_pos = 0.0;
    //double unity_neg = 0.0;
    //double unity_pos = 0.0;

    //for (unsigned int gp = 0; gp < NumGP; gp++)
    //{
    //    for (unsigned int i = 0; i < NumNodes; i++){
    //        pressure_pos += rIntEnrShapeFunctionsNeg(gp,i)*rData.Pressure[i];
    //        unity_pos += rIntEnrShapeFunctionsNeg(gp,i);
    //        pressure_neg += rIntEnrShapeFunctionsPos(gp,i)*rData.Pressure[i];
    //        unity_neg += rIntEnrShapeFunctionsPos(gp,i);
    //    }
    //}

    //pressure_pos = pressure_pos / unity_pos;
    //pressure_neg = pressure_neg / unity_neg;

    //for (int i = 0; i < NumNodes; i++)
    //{
    //    if (rData.Distance[i] > 0.0){
    //        pressure_enr_est[i] = pressure_neg - pressure_pos;
    //    } else /* if (rData.Distance[i] < 0.0) */{
    //        pressure_enr_est[i] = 0.0;
    //    } 
    //}

    //RHS_Residual = prod(LHS, pressure_enr_est);

    //for (unsigned int j = 0; j < NumNodes; j++)
    //    rRHSeeTot(j) -= RHS_Residual[j];
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureDiscontinuity(
        const double coefficient,
		MatrixType& rKeeTot,
		VectorType& rRHSeeTot)
{
    const unsigned int NumNodes = rRHSeeTot.size();
    VectorType rhs_enr = ZeroVector(NumNodes);
    MatrixType lhs_enr = ZeroMatrix(NumNodes,NumNodes);

    GeometryType::Pointer p_geom = this->pGetGeometry();

    double penalty_coeff = 0.000001;

    for (unsigned int i = 0; i < NumNodes - 1; ++i){

        const double d_i = (*p_geom)[i].FastGetSolutionStepValue(DISTANCE);
        const double curvature_i = (*p_geom)[i].FastGetSolutionStepValue(CURVATURE);
        const double di = std::abs(d_i);

        for (unsigned int j = i + 1; j < NumNodes; j++){

            const double d_j = (*p_geom)[j].FastGetSolutionStepValue(DISTANCE);
            const double curvature_j = (*p_geom)[j].FastGetSolutionStepValue(CURVATURE);
            const double dj = std::abs(d_j);

            double sum_d = di + dj;
            double Ni = dj / sum_d;
            double Nj = di / sum_d;

            const double curvature = Ni*curvature_i + Nj*curvature_j;
            
            if (d_i*d_j < 0.0){
                lhs_enr(i, i) += penalty_coeff * Ni * Ni;
                lhs_enr(i, j) -= penalty_coeff * Ni * Nj;
                lhs_enr(j, i) -= penalty_coeff * Nj * Ni;
                lhs_enr(j, j) += penalty_coeff * Nj * Nj;

                if (d_i < 0.0){
                    rhs_enr(i) -= Ni*coefficient*curvature*penalty_coeff;
                    rhs_enr(j) += Nj*coefficient*curvature*penalty_coeff;
                } else {
                    rhs_enr(i) += Ni*coefficient*curvature*penalty_coeff;
                    rhs_enr(j) -= Nj*coefficient*curvature*penalty_coeff;
                } 
            }
        }
    }

    noalias(rKeeTot) += lhs_enr;
    noalias(rRHSeeTot) += rhs_enr;

}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureGradientStabilization(
        TElementData& rData,
        MatrixType& rLHS,
        MatrixType& rV,
        VectorType& rRHS)
{
    //const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const auto &p = rData.Pressure;
    
    // Get shape function values
    //const auto &N = rData.N;
    const auto &DN = rData.DN_DX;
    //const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;

    VectorType rhs = ZeroVector(NumNodes*(Dim + 1));
    MatrixType lhs = ZeroMatrix(NumNodes*(Dim + 1), NumNodes*(Dim + 1));
    MatrixType v = ZeroMatrix(NumNodes*(Dim + 1), NumNodes);

    VectorType values = ZeroVector(NumNodes*(Dim + 1));

    const double coefficient = 1.0e3 * h*h / mu;

    for (unsigned int i = 0; i < NumNodes; i++){

        values[i*(Dim + 1) + Dim] = p[i];

        for (unsigned int j = 0; j < NumNodes; j++){

            for (unsigned int dim = 0; dim < Dim; dim++){
                lhs(i*(Dim + 1) + Dim, j*(Dim + 1) + Dim) +=
                    DN(i,dim) * DN(j,dim);

                v(i*(Dim + 1) + Dim, j) +=
                    DN(i,dim) * DNenr(j,dim);
            }
        }
    }

    rhs -= prod(lhs, values);

    noalias(rLHS) += coefficient*rData.Weight * lhs;
    noalias(rV) += coefficient*rData.Weight * v;
    noalias(rRHS) += coefficient*rData.Weight * rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rIntWeights,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rKeeTot)
{
    const unsigned int NumIntGP = rIntWeights.size();
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);

    const double coefficient = 1.0e0;

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else /* if (rData.Distance[i] < 0.0) */{
            enr_pos_interp(i, i) = 1.0;
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedIntShapeDerivativesPos = rInterfaceShapeDerivativesNeg;
    GeometryType::ShapeFunctionsGradientsType EnrichedIntShapeDerivativesNeg = rInterfaceShapeDerivativesNeg;

    for (unsigned int i = 0; i < rInterfaceShapeDerivativesNeg.size(); i++){
        EnrichedIntShapeDerivativesPos[i] = prod(enr_pos_interp, rInterfaceShapeDerivativesNeg[i]);
    }

    for (unsigned int i = 0; i < rInterfaceShapeDerivativesNeg.size(); i++){
        EnrichedIntShapeDerivativesNeg[i] = prod(enr_neg_interp, rInterfaceShapeDerivativesNeg[i]);
    }

    /* GeometryType::Pointer p_geom = this->pGetGeometry();
    std::vector<VectorType> normal(NumNodes);

    for (unsigned int j = 0; j < NumNodes; j++){   
        const VectorType normal_j = (*p_geom)[j].FastGetSolutionStepValue(DISTANCE_GRADIENT);
        const double norm_j = Kratos::norm_2(normal_j);
        normal[j] = (1.0/norm_j)*normal_j;
    } */

    for (unsigned int gp = 0; gp < NumIntGP; gp++){

        /* VectorType NdotDNenrPos = ZeroVector(NumNodes);
        VectorType NdotDNenrNeg = ZeroVector(NumNodes);

        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < Dim; dim++){
                NdotDNenrPos[j] += (normal[j])(dim) * (EnrichedIntShapeDerivativesPos[gp])(j,dim);
                NdotDNenrNeg[j] += (normal[j])(dim) * (EnrichedIntShapeDerivativesNeg[gp])(j,dim);
            }
        } */

        /* for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int j = 0; j < NumNodes; j++){
                    kee(i, j) += rIntWeights[gp] * ( NdotDNenrPos[i] - NdotDNenrNeg[i] )*
                        ( NdotDNenrPos[j] - NdotDNenrNeg[j] );
            }
        } */

        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int j = 0; j < NumNodes; j++){
                for (unsigned int dim = 0; dim < Dim; dim++){
                    kee(i, j) += rIntWeights[gp] * 
                        ( (EnrichedIntShapeDerivativesPos[gp])(i,dim) - (EnrichedIntShapeDerivativesNeg[gp])(i,dim) )*
                        ( (EnrichedIntShapeDerivativesPos[gp])(j,dim) - (EnrichedIntShapeDerivativesNeg[gp])(j,dim) );
                }
            }
        }
    }

    noalias(rKeeTot) += coefficient * kee;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntEnrShapeFunctionsPos, // Negative nodes
        const Matrix& rIntEnrShapeFunctionsNeg, // Positive nodes
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rKeeTot,
		VectorType& rRHSeeTot)
{
    const unsigned int NumIntGP = rIntWeights.size();
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);
    VectorType rhs_enr = ZeroVector(NumNodes);

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    double positive_density = 0.0;
    double negative_density = 0.0;
    double positive_viscosity = 0.0;
    double negative_viscosity = 0.0;

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
            positive_density = rData.NodalDensity[i];
            positive_viscosity = rData.NodalDynamicViscosity[i];
        } else /* if (rData.Distance[i] < 0.0) */{
            enr_pos_interp(i, i) = 1.0;
            negative_density = rData.NodalDensity[i];
            negative_viscosity = rData.NodalDynamicViscosity[i];
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedIntShapeDerivativesPos = rInterfaceShapeDerivativesNeg;
    GeometryType::ShapeFunctionsGradientsType EnrichedIntShapeDerivativesNeg = rInterfaceShapeDerivativesNeg;

    for (unsigned int i = 0; i < rInterfaceShapeDerivativesNeg.size(); i++){
        EnrichedIntShapeDerivativesPos[i] = prod(enr_pos_interp, rInterfaceShapeDerivativesNeg[i]);
    }

    for (unsigned int i = 0; i < rInterfaceShapeDerivativesNeg.size(); i++){
        EnrichedIntShapeDerivativesNeg[i] = prod(enr_neg_interp, rInterfaceShapeDerivativesNeg[i]);
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

    GeometryType::Pointer p_geom = this->pGetGeometry();
    const double h_elem = ElementSizeCalculator<Dim,NumNodes>::AverageElementSize(*p_geom);

    double cut_area = 0.0;
    for (unsigned int gp = 0; gp < NumIntGP; gp++){
        cut_area += rIntWeights[gp];
    }

    const double density = 1.0/(1.0/positive_density + 1.0/negative_density);
    const double viscosity = 1.0/(1.0/positive_viscosity + 1.0/negative_viscosity);

    // Stabilization parameters
    const double coefficient = 1.0e0;
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    const double dyn_tau = rData.DynamicTau;

    const double dt = rData.DeltaTime;

    const auto v_convection = rData.Velocity - rData.MeshVelocity;

    for (unsigned int gp = 0; gp < NumIntGP; gp++){

        Vector vconv = ZeroVector(Dim);
        double positive_weight = 0.0;
        double negative_weight = 0.0;

        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < Dim; dim++){
                vconv[dim] += 
                    (rIntEnrShapeFunctionsNeg(gp, j) + rIntEnrShapeFunctionsPos(gp, j)) * v_convection(j,dim);
            }
            positive_weight += rIntEnrShapeFunctionsNeg(gp, j);
            negative_weight += rIntEnrShapeFunctionsPos(gp, j);
        }

        const double v_conv_norm = norm_2(vconv);

        const double penalty_coefficient = coefficient *
            density * 1.0 / (dyn_tau * density / dt + stab_c1 * viscosity / h_elem / h_elem +
                                stab_c2 * density * v_conv_norm / h_elem) * element_volume / cut_area;

        for (unsigned int i = 0; i < NumNodes; i++){

            for (unsigned int j = 0; j < NumNodes; j++){

                const array_1d<double, 3> pressure_gradient_j = (*p_geom)[j].FastGetSolutionStepValue(PRESSURE_GRADIENT_AUX);

                for (unsigned int dim = 0; dim < Dim; dim++){
                    kee(i, j) += penalty_coefficient * rIntWeights[gp] * 
                        ( (EnrichedIntShapeDerivativesPos[gp])(i,dim) - (EnrichedIntShapeDerivativesNeg[gp])(i,dim) )*
                        ( (EnrichedIntShapeDerivativesPos[gp])(j,dim) - (EnrichedIntShapeDerivativesNeg[gp])(j,dim) );

                    rhs_enr(i) += penalty_coefficient * rIntWeights[gp] *
                        ( (EnrichedIntShapeDerivativesPos[gp])(i,dim) - (EnrichedIntShapeDerivativesNeg[gp])(i,dim) )*
                        (rIntEnrShapeFunctionsNeg(gp, j)/positive_weight - rIntEnrShapeFunctionsPos(gp, j)/negative_weight)*
                        pressure_gradient_j(dim);
                }
            }
        }
    }

    noalias(rKeeTot) += kee;
    noalias(rRHSeeTot) += rhs_enr;

}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::GhostPressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rWeights,
        const GeometryType::ShapeFunctionsGradientsType& rShapeDerivatives,
        MatrixType& rKeeTot)
{
    const unsigned int NumGP = rWeights.size();
    MatrixType kee = ZeroMatrix(NumNodes, NumNodes);

    const double coefficient = 1.0e4;

    Matrix enr_neg_interp = ZeroMatrix(NumNodes, NumNodes);
    Matrix enr_pos_interp = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rData.Distance[i] > 0.0){
            enr_neg_interp(i, i) = 1.0;
        } else /* if (rData.Distance[i] < 0.0) */ {
            enr_pos_interp(i, i) = 1.0;
        }
    }

    GeometryType::ShapeFunctionsGradientsType EnrichedShapeDerivativesPos = rShapeDerivatives;
    GeometryType::ShapeFunctionsGradientsType EnrichedShapeDerivativesNeg = rShapeDerivatives;

    for (unsigned int i = 0; i < rShapeDerivatives.size(); i++){
        EnrichedShapeDerivativesPos[i] = prod(enr_pos_interp, rShapeDerivatives[i]);
    }

    for (unsigned int i = 0; i < rShapeDerivatives.size(); i++){
        EnrichedShapeDerivativesNeg[i] = prod(enr_neg_interp, rShapeDerivatives[i]);
    }

    /* GeometryType::Pointer p_geom = this->pGetGeometry();
    std::vector<VectorType> normal(NumNodes);

    for (unsigned int j = 0; j < NumNodes; j++){   
        const VectorType normal_j = (*p_geom)[j].FastGetSolutionStepValue(DISTANCE_GRADIENT);
        const double norm_j = Kratos::norm_2(normal_j);
        normal[j] = (1.0/norm_j)*normal_j;
    } */

    for (unsigned int gp = 0; gp < NumGP; gp++){

        /* VectorType NdotDNenrPos = ZeroVector(NumNodes);
        VectorType NdotDNenrNeg = ZeroVector(NumNodes);

        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < Dim; dim++){
                NdotDNenrPos[j] += (normal[j])(dim) * (EnrichedIntShapeDerivativesPos[gp])(j,dim);
                NdotDNenrNeg[j] += (normal[j])(dim) * (EnrichedIntShapeDerivativesNeg[gp])(j,dim);
            }
        } */

        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int j = 0; j < NumNodes; j++){
                for (unsigned int dim = 0; dim < Dim; dim++){
                    kee(i, j) += rWeights[gp] * 
                        ( (EnrichedShapeDerivativesNeg[gp])(i,dim) - (EnrichedShapeDerivativesPos[gp])(i,dim) )*
                        ( (EnrichedShapeDerivativesNeg[gp])(j,dim) - (EnrichedShapeDerivativesPos[gp])(j,dim) );
                }
            }
        }
    }

    noalias(rKeeTot) += coefficient * kee;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PressureDiscontinuityandSurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const std::vector<Vector>& rIntNormalsNeg,
        const Matrix& rIntShapeFunctions,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
        Vector& rSurfaceTensionForce,
		MatrixType& rKeeTot,
		VectorType& rRHSeeTot)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    rSurfaceTensionForce = ZeroVector(NumDim);

    //A new set of test functions is introduced to retain the symmetry of the modification made to KeeTot
    Matrix TestFunctions = ZeroMatrix(NumGP,NumNodes);

    for (unsigned int gp = 0; gp < NumGP; gp++)
    {
        for (unsigned int dim = 0; dim < NumDim; dim++)
            rSurfaceTensionForce(dim) += coefficient*(rIntNormalsNeg[gp])[dim]*rCurvature(gp)*rIntWeights(gp);

        for (unsigned int j = 0; j < NumNodes; j++)
            TestFunctions(gp, j) = rIntEnrShapeFunctionsNeg(gp,j) - rIntEnrShapeFunctionsPos(gp,j);
    }


    // Penalty coefficient
    const double penalty_coeff = 1000.0;

    for (unsigned int j = 0; j < NumNodes; j++){

        double RHSj = 0.0;

        for (unsigned int gp = 0; gp < NumGP; gp++){
            RHSj += rCurvature(gp)*rIntWeights(gp)*TestFunctions(gp,j);
        }

        rRHSeeTot(j) += coefficient*RHSj*penalty_coeff;

        for (unsigned int i = 0; i < NumNodes; i++){

            double LHSij = 0.0;

            for (unsigned int gp = 0; gp < NumGP; gp++){
                LHSij += (rIntEnrShapeFunctionsNeg(gp,i) - rIntEnrShapeFunctionsPos(gp,i))*rIntWeights(gp)*TestFunctions(gp,j);
            }

            rKeeTot(i,j) += LHSij*penalty_coeff;
        }
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTensionDP(
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
        MatrixType& rVLHS)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = NumNodes - 1;

    MatrixType vlhs = ZeroMatrix(NumNodes * (NumDim + 1), NumNodes);

    GeometryType::Pointer p_geom = this->pGetGeometry();

    std::vector<VectorType> normal(NumNodes);

    for (unsigned int j = 0; j < NumNodes; j++){   
        const VectorType normal_j = (*p_geom)[j].FastGetSolutionStepValue(DISTANCE_GRADIENT);
        const double norm_j = Kratos::norm_2(normal_j);
        normal[j] = (1.0/norm_j)*normal_j;
    }

    for (unsigned int gp = 0; gp < NumGP; gp++){
        
        VectorType normal_gp = ZeroVector(NumDim);
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            normal_gp += rIntShapeFunctions(gp, j)*normal[j];
        }

        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int j = 0; j < NumNodes; j++){ 
                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    const double temp = (rIntEnrShapeFunctionsNeg(gp,j) - rIntEnrShapeFunctionsPos(gp,j))*
                        rIntWeights(gp)*rIntShapeFunctions(gp,i)*normal_gp(dimi);
                    vlhs(i * (NumDim + 1) + dimi, j) += temp;
                }
            }
        }
    }

    noalias(rVLHS) += vlhs;

} 

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTensionMixed(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        VectorType& rRHS)
{

    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rInterfaceShapeDerivativesNeg[0].size2();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    GeometryType::Pointer p_geom = this->pGetGeometry();

    std::vector<VectorType> normal(NumNodes);

    for (unsigned int j = 0; j < NumNodes; j++){   
        const VectorType normal_j = (*p_geom)[j].FastGetSolutionStepValue(DISTANCE_GRADIENT);
        const double norm_j = Kratos::norm_2(normal_j);
        normal[j] = (1.0/norm_j)*normal_j;
    }

    for (unsigned int gp = 0; gp < NumGP; gp++){

        MatrixType P_gp = ZeroMatrix(NumDim, NumDim);
        VectorType normal_gp = ZeroVector(NumDim);

        for (unsigned int j=0; j < NumNodes; ++j){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                normal_gp(dim) += (*p_geom)[j].FastGetSolutionStepValue(DISTANCE)*(rInterfaceShapeDerivativesNeg[gp])(j,dim);
            }
        }

        const double norm = Kratos::norm_2(normal_gp);
        normal_gp = (1.0/norm)*normal_gp;

        for (unsigned int dimi = 0; dimi < NumDim; dimi++){
            for (unsigned int dimj = 0; dimj < NumDim; dimj++){
                P_gp(dimi, dimj) = -normal_gp(dimi)*normal_gp(dimj);
            }
            P_gp(dimi, dimi) += 1;
        }

        for (unsigned int i = 0; i < NumNodes; i++){           
            /* for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                for (unsigned int dimj = 0; dimj < NumDim; dimj++){                    
                    rhs[ i*(NumDim+1) + dimi ] -= 
                        coefficient*P_gp(dimi, dimj)*(rInterfaceShapeDerivativesNeg[gp])(i,dimj)*rIntWeights(gp);
                }
            } */

            for (unsigned int j = 0; j < NumNodes; j++){   
                /* for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    for (unsigned int dimj = 0; dimj < NumDim; dimj++){                    
                        rhs[ i*(NumDim+1) + dimi ] += 
                            coefficient*(-(normal[j])(dimi)*(normal[j])(dimj))*rIntShapeFunctions(gp,j)*(rInterfaceShapeDerivativesNeg[gp])(i,dimj)*rIntWeights(gp);
                    }
                    rhs[ i*(NumDim+1) + dimi ] += 
                            coefficient*(1.0)*rIntShapeFunctions(gp,j)*(rInterfaceShapeDerivativesNeg[gp])(i,dimi)*rIntWeights(gp);
                } */

                const double curvature_j = (*p_geom)[j].FastGetSolutionStepValue(CURVATURE);
                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    rhs[ i*(NumDim+1) + dimi ] -= 
                        coefficient*(normal[j])(dimi)*curvature_j*rIntShapeFunctions(gp,j)*rIntShapeFunctions(gp,i)*rIntWeights(gp);
                }
            }
        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const std::vector<Vector>& rIntNormalsNeg,
        Vector& rSurfaceTensionForce)
{
    const unsigned int NumGP = rIntWeights.size();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    rSurfaceTensionForce = ZeroVector(NumDim);

    for (unsigned int gp = 0; gp < NumGP; gp++){

        //double normal_norm = 0.0;
        //for (unsigned int dim = 0; dim < NumDim; dim++){
        //    normal_norm += (rIntNormalsNeg[gp])[dim]*(rIntNormalsNeg[gp])[dim];}
        //normal_norm = std::sqrt(normal_norm);

        for (unsigned int dim = 0; dim < NumDim; dim++){
            rSurfaceTensionForce(dim) -= coefficient*((rIntNormalsNeg[gp])[dim]/* /normal_norm */)*rCurvature(gp)*rIntWeights(gp);}

    }
} 

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        VectorType& rRHS)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    for (unsigned int gp = 0; gp < NumGP; gp++){
        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                rhs[ j*(NumDim+1) + dim ] -= coefficient*(rIntNormalsNeg[gp])[dim]*rCurvature(gp)*rIntWeights(gp)*rIntShapeFunctions(gp,j);
            }
        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        const std::vector<Vector>& rIntNormalsNeg,
        VectorType& rRHS)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    for (unsigned int gp = 0; gp < NumGP; gp++){

        MatrixType P_gp = ZeroMatrix(NumDim, NumDim);

        for (unsigned int dimi = 0; dimi < NumDim; dimi++){
            for (unsigned int dimj = 0; dimj < NumDim; dimj++){
                P_gp(dimi, dimj) = -(rIntNormalsNeg[gp])(dimi)*(rIntNormalsNeg[gp])(dimj);
            }
            P_gp(dimi, dimi) += 1;
        }
        //KRATOS_INFO("P Matrix") << P_gp << std::endl;

        for (unsigned int i = 0; i < NumNodes; i++){           
            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                for (unsigned int dimj = 0; dimj < NumDim; dimj++){                    
                    rhs[ i*(NumDim+1) + dimi ] -= coefficient*P_gp(dimi, dimj)*(rInterfaceShapeDerivativesNeg[gp])(i,dimj)*rIntWeights(gp);
                }
            }
        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        VectorType& rRHS)
{

    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rInterfaceShapeDerivativesNeg[0].size2();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    GeometryType::Pointer p_geom = this->pGetGeometry();

    for (unsigned int gp = 0; gp < NumGP; gp++){

        MatrixType P_gp = ZeroMatrix(NumDim, NumDim);

        VectorType normal_gp = ZeroVector(NumDim);

        for (unsigned int i=0; i < NumNodes; ++i){

            /* for (unsigned int dim = 0; dim < NumDim; dim++){
                normal_gp(dim) += (*p_geom)[i].FastGetSolutionStepValue(DISTANCE)*(rInterfaceShapeDerivativesNeg[gp])(i,dim);
            } */

            normal_gp += (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT)*rIntShapeFunctions(gp,i);
        }

        double norm = 0.0;

        for (unsigned int dim = 0; dim < NumDim; dim++){
           norm += normal_gp(dim)*normal_gp(dim);
        }

        norm = std::sqrt(norm);
        normal_gp = (1.0/norm)*normal_gp;

        for (unsigned int dimi = 0; dimi < NumDim; dimi++){
            for (unsigned int dimj = 0; dimj < NumDim; dimj++){
                P_gp(dimi, dimj) = -normal_gp(dimi)*normal_gp(dimj);
            }
            P_gp(dimi, dimi) += 1;
        }
        //KRATOS_INFO("P Matrix") << P_gp << std::endl;

        for (unsigned int i = 0; i < NumNodes; i++){           
            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                for (unsigned int dimj = 0; dimj < NumDim; dimj++){                    
                    rhs[ i*(NumDim+1) + dimi ] -= coefficient*P_gp(dimi, dimj)*(rInterfaceShapeDerivativesNeg[gp])(i,dimj)*rIntWeights(gp);
                }
            }
        }
    }

    noalias(rRHS) += rhs;

}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        TElementData& rData,
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rLHS,
        VectorType& rRHS)
{
    const unsigned int NumGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rInterfaceShapeDerivativesNeg[0].size2();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    MatrixType lhs_ST = ZeroMatrix(NumNodes*(NumDim+1),NumNodes*(NumDim+1));

    VectorType values_vel = ZeroVector(NumNodes*(NumDim+1));

    GeometryType::Pointer p_geom = this->pGetGeometry();

    for (unsigned int gp = 0; gp < NumGP; gp++){

        MatrixType P_gp = ZeroMatrix(NumDim, NumDim);

        VectorType normal_gp = ZeroVector(NumDim);

        for (unsigned int i=0; i < NumNodes; ++i){

            /* for (unsigned int dim = 0; dim < NumDim; dim++){
                normal_gp(dim) += (*p_geom)[i].FastGetSolutionStepValue(DISTANCE)*(rInterfaceShapeDerivativesNeg[gp])(i,dim);
            } */

            normal_gp += (*p_geom)[i].FastGetSolutionStepValue(DISTANCE_GRADIENT)*rIntShapeFunctions(gp,i);
        }

        double norm = 0.0;

        for (unsigned int dim = 0; dim < NumDim; dim++){
           norm += normal_gp(dim)*normal_gp(dim);
        }

        norm = std::sqrt(norm);
        normal_gp = (1.0/norm)*normal_gp;

        for (unsigned int dimi = 0; dimi < NumDim; dimi++){
            for (unsigned int dimj = 0; dimj < NumDim; dimj++){
                P_gp(dimi, dimj) = -normal_gp(dimi)*normal_gp(dimj);
            }
            P_gp(dimi, dimi) += 1;
        }
        //KRATOS_INFO("P Matrix") << P_gp << std::endl;

        for (unsigned int i = 0; i < NumNodes; i++){ 
            Vector velocity0 = (*p_geom)[i].FastGetSolutionStepValue(VELOCITY);

            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                values_vel[ i*(NumDim+1) + dimi ] = velocity0[dimi];

                for (unsigned int dimj = 0; dimj < NumDim; dimj++){                    
                    rhs[ i*(NumDim+1) + dimi ] -= 
                        coefficient*P_gp(dimi, dimj)*(rInterfaceShapeDerivativesNeg[gp])(i,dimj)*rIntWeights(gp);
                     
                }
            }

            double n_dot_grad_N_I = 0;
            for (unsigned int k = 0; k < NumDim; k++){
                n_dot_grad_N_I += normal_gp[k]*(rInterfaceShapeDerivativesNeg[gp])(i,k);
            }

            for (unsigned int j = 0; j < NumNodes; j++){

                double n_dot_grad_N_J = 0;
                for (unsigned int l = 0; l < NumDim; l++){
                    n_dot_grad_N_J += normal_gp[l]*(rInterfaceShapeDerivativesNeg[gp])(j,l);
                }

                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    for (unsigned int dimj = 0; dimj < NumDim; dimj++){ 
                        lhs_ST( i*(NumDim+1) + dimi, j*(NumDim+1) + dimj) +=
                            + 2.0 * coefficient * rData.DeltaTime * n_dot_grad_N_I * ( (rInterfaceShapeDerivativesNeg[gp])(j,dimj) 
                            - n_dot_grad_N_J * normal_gp[dimj] ) * normal_gp[dimi] * rIntWeights(gp);

                        /* lhs_ST( i*(NumDim+1) + dimi, j*(NumDim+1) + dimi) +=
                            - 2.0 * coefficient * rData.DeltaTime * (rInterfaceShapeDerivativesNeg[gp])(i,dimj) 
                            * (rInterfaceShapeDerivativesNeg[gp])(j,dimj) * rIntWeights(gp); */
                    }
                }
            }
        }
    }

    noalias(rRHS) += rhs /* - prod(lhs_ST, values_vel) */; //Be careful about what is done here!!!
    noalias(rLHS) += lhs_ST;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        const Kratos::Vector& rCLWeights,
        const Matrix& rCLShapeFunctions,
        const Vector& rTangential,
        bool HasContactLine,
        VectorType& rRHS)
{
    const unsigned int NumIntGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    double total_weight = 0.0;
    Vector normal_avg = ZeroVector(NumDim);

    for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
        total_weight += rIntWeights(intgp);
        normal_avg += rIntWeights(intgp)*rIntNormalsNeg[intgp];

        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                rhs[ j*(NumDim+1) + dim ] -= coefficient*(rIntNormalsNeg[intgp])[dim]*rCurvature(intgp)*rIntWeights(intgp)*rIntShapeFunctions(intgp,j);
            }
        }
    }

    normal_avg = 1.0/total_weight * normal_avg;

    GeometryType::Pointer p_geom = this->pGetGeometry();

    Vector contact_vector = ZeroVector(NumDim);
    Vector wall_tangent = ZeroVector(NumDim);
    std::vector< Vector > wall_normal;

    if (HasContactLine){

        const unsigned int NumCLGP = rCLShapeFunctions.size1();    

        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){

            Vector wall_normal_gp = ZeroVector(NumDim);

            for (unsigned int j = 0; j < NumNodes; j++){
                wall_normal_gp += rCLShapeFunctions(clgp,j)*(*p_geom)[j].FastGetSolutionStepValue(NORMAL);
            }

            wall_normal.push_back( wall_normal_gp );
        }

        //KRATOS_INFO("Cut Element, has contact line, NumCLGP") << NumCLGP << std::endl;

        MathUtils<double>::UnitCrossProduct(contact_vector, rTangential, normal_avg);

        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){

            //KRATOS_INFO("Cut Element, has contact line, CLWeight") << rCLWeights(clgp) << std::endl;
            MathUtils<double>::UnitCrossProduct(wall_tangent, wall_normal[clgp], rTangential);

            for (unsigned int j = 0; j < NumNodes; j++){

                //KRATOS_INFO("Cut Element, has contact line, CLShapeFunction") << rCLShapeFunctions(clgp,j) << std::endl;

                for (unsigned int dim = 0; dim < NumDim; dim++){
                    rhs[ j*(NumDim+1) + dim ] -= coefficient*contact_vector[dim]*rCLWeights(clgp)*rCLShapeFunctions(clgp,j);
                    rhs[ j*(NumDim+1) + dim ] -= 1.0*coefficient*wall_tangent[dim]*rCLWeights(clgp)*rCLShapeFunctions(clgp,j);
                    //KRATOS_INFO("Cut Element, has contact line, CLShapeFunction") << rCLShapeFunctions(clgp,j) << std::endl;
                    //KRATOS_INFO("Cut Element, has contact line, RHS") << rhs[ j*(NumDim+1) + dim ] << std::endl;               
                }
            }
        }
    }

    const VectorType zero_vector = ZeroVector(NumDim);
    for (unsigned int i=0; i < NumNodes; ++i){

        #pragma omp critical
        {
        (*p_geom)[i].FastGetSolutionStepValue(NORMAL_VECTOR) = normal_avg;
        //if ((*p_geom)[i].GetValue(IS_STRUCTURE) == 1.0 && HasContactLine){
            //if (HasContactLine){
                (*p_geom)[i].FastGetSolutionStepValue(TANGENT_VECTOR) = wall_tangent;
                (*p_geom)[i].FastGetSolutionStepValue(CONTACT_VECTOR) = contact_vector;
        }

            //}
        //}
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
    const TElementData& rData,
    const double coefficient,
    const double coefficientS,
    const double zeta,
    const Kratos::Vector& rCurvature,
    const Kratos::Vector& rIntWeights,
    const Matrix& rIntShapeFunctions,
    const std::vector<Vector>& rIntNormalsNeg,
    const std::vector<Kratos::Vector>& rCLWeights,
    const std::vector<Matrix>& rCLShapeFunctions,
    const std::vector<Vector>& rTangential,
    //bool HasContactLine,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    const unsigned int NumIntGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    //double total_weight = 0.0;
    Vector normal_avg = ZeroVector(NumDim);

    for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
        //total_weight += rIntWeights(intgp);
        normal_avg += rIntWeights(intgp)*rIntNormalsNeg[intgp];

        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                rhs[ j*(NumDim+1) + dim ] -= coefficient*(rIntNormalsNeg[intgp])[dim]*rCurvature(intgp)*rIntWeights(intgp)*rIntShapeFunctions(intgp,j);
            }
        }
    }

    //normal_avg = 1.0/total_weight * normal_avg;
    normal_avg /= norm_2(normal_avg);

    GeometryType::Pointer p_geom = this->pGetGeometry();

    //if (HasContactLine){
    for (unsigned int i_cl = 0; i_cl < rCLWeights.size(); i_cl++){

        Vector contact_vector = ZeroVector(NumDim);
        Vector wall_tangent = ZeroVector(NumDim);
        std::vector< Vector > wall_normal;

        MatrixType lhs_dissipation = ZeroMatrix(NumNodes*(NumDim+1),NumNodes*(NumDim+1));

        const unsigned int n_dim = 3;  // NumDim did not compiled!
        const unsigned int n_nodes = n_dim + 1;

        Kratos::array_1d<double,n_nodes*(n_dim+1)> tempU; // Only velocity
        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                tempU[i*(n_dim+1) + dimi] = rData.Velocity(i,dimi);
            }
        }

        const unsigned int NumCLGP = (rCLShapeFunctions[i_cl]).size1();    
        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){
            Vector wall_normal_gp = ZeroVector(NumDim);
            for (unsigned int j = 0; j < NumNodes; j++){
                wall_normal_gp += (rCLShapeFunctions[i_cl])(clgp,j)*(*p_geom)[j].FastGetSolutionStepValue(NORMAL);
            }
            wall_normal.push_back( wall_normal_gp );
        }

        //KRATOS_INFO("Cut Element, has contact line, NumCLGP") << NumCLGP << std::endl;

        VectorType vel0_CL = ZeroVector(NumDim);
        double weight_sum = 0.0;

        MathUtils<double>::UnitCrossProduct(contact_vector, rTangential[i_cl], normal_avg);

        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){
            weight_sum += (rCLWeights[i_cl])[clgp];

            //KRATOS_INFO("Cut Element, has contact line, CLWeight") << rCLWeights(clgp) << std::endl;
            MathUtils<double>::UnitCrossProduct(wall_tangent, wall_normal[clgp], rTangential[i_cl]);

            for (unsigned int i = 0; i < NumNodes; i++){
                const VectorType Vel0 = (*p_geom)[i].FastGetSolutionStepValue(VELOCITY, 1);
                vel0_CL += (rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i)*Vel0;

                //KRATOS_INFO("Cut Element, has contact line, CLShapeFunction") << rCLShapeFunctions(clgp,j) << std::endl;

                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    rhs[ i*(NumDim+1) + dimi ] -= coefficient*contact_vector[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i);
                    rhs[ i*(NumDim+1) + dimi ] += coefficientS*wall_tangent[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i); //Contac-line tangential force
                    //KRATOS_INFO("Cut Element, has contact line, CLShapeFunction") << rCLShapeFunctions(clgp,j) << std::endl;
                    //KRATOS_INFO("Cut Element, has contact line, RHS") << rhs[ j*(NumDim+1) + dim ] << std::endl;

                    for (unsigned int j = 0; j < NumNodes; j++){
                        //for (unsigned int dimj = 0; dimj < NumDim; dimj++){ //Be carefull not to confuse
                            // Dimension_j is useless since the force acts in the dimension_i direction
                            lhs_dissipation( i*(NumDim+1) + dimi, j*(NumDim+1) + dimi) += 
                                zeta * (rCLWeights[i_cl])[clgp] * (rCLShapeFunctions[i_cl])(clgp,j) * (rCLShapeFunctions[i_cl])(clgp,i);
                        //}
                    }
                }
            }
        }

        noalias(rLHS) += lhs_dissipation;
        noalias(rRHS) -= prod(lhs_dissipation,tempU);

        const double angle = acos (inner_prod(wall_tangent,contact_vector)) * 180.0 / PI;
        //KRATOS_INFO("HasContact") << angle << std::endl;
        this->SetValue(CONTACT_ANGLE, angle);

        vel0_CL = 1.0/weight_sum*vel0_CL;
        const double vel_CL = inner_prod(wall_tangent,vel0_CL);
        this->SetValue(CONTACT_VELOCITY, vel_CL);

        for (unsigned int i=0; i < NumNodes; ++i){

            #pragma omp critical
            {
            (*p_geom)[i].FastGetSolutionStepValue(NORMAL_VECTOR) = normal_avg;

            //if ((*p_geom)[i].GetValue(IS_STRUCTURE) == 1.0){
                //if (HasContactLine){
                    (*p_geom)[i].FastGetSolutionStepValue(TANGENT_VECTOR) = wall_tangent;
                    (*p_geom)[i].FastGetSolutionStepValue(CONTACT_VECTOR) = contact_vector;
            }
                //}
            //}
        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
    const TElementData& rData,
    const double coefficient,
    const double coefficientS,
    const double zeta,
    const double micro_length_scale,
    const Kratos::Vector& rCurvature,
    const Kratos::Vector& rIntWeights,
    const Matrix& rIntShapeFunctions,
    const std::vector<Vector>& rIntNormalsNeg,
    const std::vector<Kratos::Vector>& rCLWeights,
    const std::vector<Matrix>& rCLShapeFunctions,
    const std::vector<Vector>& rTangential,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    const unsigned int NumIntGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                rhs[ j*(NumDim+1) + dim ] -= coefficient*(rIntNormalsNeg[intgp])[dim]*rCurvature(intgp)*rIntWeights(intgp)*rIntShapeFunctions(intgp,j);
            }
        }
    }

    GeometryType::Pointer p_geom = this->pGetGeometry();

    for (unsigned int i_cl = 0; i_cl < rCLWeights.size(); i_cl++){
        MatrixType lhs_dissipation = ZeroMatrix(NumNodes*(NumDim+1),NumNodes*(NumDim+1));

        Vector contact_vector_macro = ZeroVector(NumDim);
        Vector contact_vector_micro = ZeroVector(NumDim);
        Vector wall_tangent = ZeroVector(NumDim);
        Vector wall_normal_gp = ZeroVector(NumDim);
        Vector velocity_gp = ZeroVector(NumDim);
        double contact_velocity = 0.0;
        double contact_angle_macro = 0.0;
        double contact_angle_micro = 0.0;

        Vector normal_avg = ZeroVector(NumDim);
        for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
            normal_avg += rIntWeights(intgp)*rIntNormalsNeg[intgp];
        }
        normal_avg /= norm_2(normal_avg);

        Vector tempU = ZeroVector(NumNodes*(NumDim+1)); // Only velocity
        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                tempU[i*(NumDim+1) + dimi] = rData.Velocity(i,dimi);
            }
        }

        double positive_density = 0.0;
        double negative_density = 0.0;
        double positive_viscosity = 0.0;
        double negative_viscosity = 0.0;

        for (unsigned int i = 0; i < NumNodes; i++){
            if (rData.Distance[i] > 0.0){
                positive_density = rData.NodalDensity[i];
                positive_viscosity = rData.NodalDynamicViscosity[i];
            } else /* if (rData.Distance[i] < 0.0) */{
                negative_density = rData.NodalDensity[i];
                negative_viscosity = rData.NodalDynamicViscosity[i];
            }
        }
        const double effective_density = 0.5*(positive_density + negative_density);
        const double effective_viscosity = 0.5*(positive_viscosity + negative_viscosity);
        //const double element_size = ElementSizeCalculator<3,4>::AverageElementSize(*p_geom);

        const unsigned int NumCLGP = (rCLShapeFunctions[i_cl]).size1();
        MathUtils<double>::UnitCrossProduct(contact_vector_macro, rTangential[i_cl], normal_avg);

        double weight_sum = 0.0;
        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){
            weight_sum += (rCLWeights[i_cl])[clgp];

            wall_normal_gp = ZeroVector(NumDim);
            velocity_gp = ZeroVector(NumDim);
            for (unsigned int j = 0; j < NumNodes; j++){
                wall_normal_gp += (rCLShapeFunctions[i_cl])(clgp,j)
                            *(*p_geom)[j].FastGetSolutionStepValue(NORMAL);
                velocity_gp += (rCLShapeFunctions[i_cl])(clgp,j)*
                            (*p_geom)[j].FastGetSolutionStepValue(VELOCITY);
            }
            MathUtils<double>::UnitCrossProduct(wall_tangent, wall_normal_gp, rTangential[i_cl]);
            const double element_size = ElementSizeCalculator<3,4>::ProjectedElementSize(*p_geom, wall_normal_gp);

            const double contact_angle_macro_gp = std::acos(inner_prod(wall_tangent,contact_vector_macro));
            double contact_angle_micro_gp = contact_angle_macro_gp;
            const double contact_angle_equilibrium = std::acos( coefficientS/coefficient );

            const double contact_velocity_gp = inner_prod(wall_tangent,velocity_gp);

            //const double reynolds_number = effective_density*std::abs(contact_velocity_gp)*element_size/effective_viscosity;
            const double capilary_number = effective_viscosity*contact_velocity_gp/coefficient;

            //KRATOS_INFO("two fluids NS") << "capilary_number= " << capilary_number << std::endl;
            //KRATOS_INFO("two fluids NS") << "reynolds_number= " << reynolds_number << std::endl;

            //KRATOS_INFO("two fluids NS") << "angle difference= " << std::abs(contact_angle_macro_gp - contact_angle_equilibrium)*180/PI << std::endl;

            if ( std::abs(contact_angle_macro_gp - contact_angle_equilibrium) < 6.0e-1 &&
                    capilary_number < 3.0e-1){
                const double cubic_contact_angle_micro_gp = std::pow(contact_angle_macro_gp, 3.0)
                    - 9*capilary_number*std::log(element_size/micro_length_scale);

                KRATOS_WARNING_IF("TwoFluidsNS", cubic_contact_angle_micro_gp < 0.0 ||
                                cubic_contact_angle_micro_gp > 31.0)
                            << "Hydrodynamics theory failed to estimate micro contact-angle (large slip velocity)." 
                            << std::endl;

                if (cubic_contact_angle_micro_gp >= 0.0 &&
                        cubic_contact_angle_micro_gp <= 31.0) //std::pow(PI, 3.0))
                    contact_angle_micro_gp = std::pow(cubic_contact_angle_micro_gp, 1.0/3.0);
                else if (cubic_contact_angle_micro_gp < 0.0)
                    contact_angle_micro_gp = 0.0; //contact_angle_equilibrium;
                else //if (cubic_contact_angle_micro_gp > 31.0)
                    contact_angle_micro_gp = PI; //contact_angle_equilibrium;

                // This relation is valid for contact_angle < 3PI/4 and vanishing Reynolds & Capillary numbers
            }

            //KRATOS_INFO("two fluids NS") << "element_size= " << element_size << std::endl;
            //KRATOS_INFO("two fluids NS") << "contact_angle_macro_gp= " << contact_angle_macro_gp << std::endl;
            //KRATOS_INFO("two fluids NS") << "contact_angle_micro_gp= " << contact_angle_micro_gp << std::endl;


            contact_vector_micro = std::cos(contact_angle_micro_gp)*wall_tangent +
                    std::sin(contact_angle_micro_gp)*wall_normal_gp;

            for (unsigned int i = 0; i < NumNodes; i++){
                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    rhs[ i*(NumDim+1) + dimi ] -= coefficient*contact_vector_micro[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i);
                    rhs[ i*(NumDim+1) + dimi ] += coefficientS*wall_tangent[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i); //Contac-line tangential force

                    for (unsigned int j = 0; j < NumNodes; j++){
                            lhs_dissipation( i*(NumDim+1) + dimi, j*(NumDim+1) + dimi) +=
                                zeta * (rCLWeights[i_cl])[clgp] * (rCLShapeFunctions[i_cl])(clgp,j) * (rCLShapeFunctions[i_cl])(clgp,i);
                    }
                }
            }
            contact_velocity += (rCLWeights[i_cl])[clgp]*contact_velocity_gp;
            contact_angle_macro += (rCLWeights[i_cl])[clgp]*contact_angle_macro_gp;
            contact_angle_micro += (rCLWeights[i_cl])[clgp]*contact_angle_micro_gp;
        }

        noalias(rLHS) += lhs_dissipation;
        rhs -= prod(lhs_dissipation,tempU);

        contact_angle_macro /= weight_sum;
        this->SetValue(CONTACT_ANGLE, contact_angle_macro*180.0/PI);

        contact_angle_micro /= weight_sum;
        this->SetValue(CONTACT_ANGLE_MICRO, contact_angle_micro*180.0/PI);

        contact_velocity /= weight_sum;
        this->SetValue(CONTACT_VELOCITY, contact_velocity);

        for (unsigned int i=0; i < NumNodes; ++i){

            #pragma omp critical
            {
            (*p_geom)[i].FastGetSolutionStepValue(NORMAL_VECTOR) = normal_avg;
            (*p_geom)[i].FastGetSolutionStepValue(TANGENT_VECTOR) = wall_tangent;
            (*p_geom)[i].FastGetSolutionStepValue(CONTACT_VECTOR) = contact_vector_macro;
            }

        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::SurfaceTension(
    const TElementData& rData,
    const double coefficient,
    const double theta_advancing,
    const double theta_receding,
    const double zeta,
    const double micro_length_scale,
    const Kratos::Vector& rCurvature,
    const Kratos::Vector& rIntWeights,
    const Matrix& rIntShapeFunctions,
    const std::vector<Vector>& rIntNormalsNeg,
    const std::vector<Kratos::Vector>& rCLWeights,
    const std::vector<Matrix>& rCLShapeFunctions,
    const std::vector<Vector>& rTangential,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    const unsigned int NumIntGP = rIntShapeFunctions.size1();
    const unsigned int NumNodes = rIntShapeFunctions.size2();
    const unsigned int NumDim = rIntNormalsNeg[0].size();

    VectorType rhs = ZeroVector(NumNodes*(NumDim+1));

    for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
        for (unsigned int j = 0; j < NumNodes; j++){
            for (unsigned int dim = 0; dim < NumDim; dim++){
                rhs[ j*(NumDim+1) + dim ] -= coefficient*(rIntNormalsNeg[intgp])[dim]*rCurvature(intgp)*rIntWeights(intgp)*rIntShapeFunctions(intgp,j);
            }
        }
    }

    GeometryType::Pointer p_geom = this->pGetGeometry();

    for (unsigned int i_cl = 0; i_cl < rCLWeights.size(); i_cl++){
        MatrixType lhs_dissipation = ZeroMatrix(NumNodes*(NumDim+1),NumNodes*(NumDim+1));

        Vector contact_vector_macro = ZeroVector(NumDim);
        Vector contact_vector_micro = ZeroVector(NumDim);
        Vector wall_tangent = ZeroVector(NumDim);
        Vector wall_normal_gp = ZeroVector(NumDim);
        Vector velocity_gp = ZeroVector(NumDim);
        double contact_velocity = 0.0;
        double contact_angle_macro = 0.0;
        double contact_angle_micro = 0.0;

        Vector normal_avg = ZeroVector(NumDim);
        for (unsigned int intgp = 0; intgp < NumIntGP; intgp++){
            normal_avg += rIntWeights(intgp)*rIntNormalsNeg[intgp];
        }
        normal_avg /= norm_2(normal_avg);

        Vector tempU = ZeroVector(NumNodes*(NumDim+1)); // Only velocity
        for (unsigned int i = 0; i < NumNodes; i++){
            for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                tempU[i*(NumDim+1) + dimi] = rData.Velocity(i,dimi);
            }
        }

        double positive_density = 0.0;
        double negative_density = 0.0;
        double positive_viscosity = 0.0;
        double negative_viscosity = 0.0;

        for (unsigned int i = 0; i < NumNodes; i++){
            if (rData.Distance[i] > 0.0){
                positive_density = rData.NodalDensity[i];
                positive_viscosity = rData.NodalDynamicViscosity[i];
            } else /* if (rData.Distance[i] < 0.0) */{
                negative_density = rData.NodalDensity[i];
                negative_viscosity = rData.NodalDynamicViscosity[i];
            }
        }
        const double effective_density = 0.5*(positive_density + negative_density);
        const double effective_viscosity = 0.5*(positive_viscosity + negative_viscosity);
        //const double element_size = ElementSizeCalculator<3,4>::AverageElementSize(*p_geom);

        const unsigned int NumCLGP = (rCLShapeFunctions[i_cl]).size1();
        MathUtils<double>::UnitCrossProduct(contact_vector_macro, rTangential[i_cl], normal_avg);

        double weight_sum = 0.0;
        for (unsigned int clgp = 0; clgp < NumCLGP; clgp++){
            weight_sum += (rCLWeights[i_cl])[clgp];

            wall_normal_gp = ZeroVector(NumDim);
            velocity_gp = ZeroVector(NumDim);
            double avg_contact_angle = 0.0;
            for (unsigned int j = 0; j < NumNodes; j++){
                wall_normal_gp += (rCLShapeFunctions[i_cl])(clgp,j)
                            *(*p_geom)[j].FastGetSolutionStepValue(NORMAL);
                velocity_gp += (rCLShapeFunctions[i_cl])(clgp,j)*
                            (*p_geom)[j].FastGetSolutionStepValue(VELOCITY);
                avg_contact_angle += (rCLShapeFunctions[i_cl])(clgp,j)*
                            (*p_geom)[j].FastGetSolutionStepValue(CONTACT_ANGLE);
            }
            MathUtils<double>::UnitCrossProduct(wall_tangent, wall_normal_gp, rTangential[i_cl]);
            const double element_size = ElementSizeCalculator<3,4>::ProjectedElementSize(*p_geom, wall_normal_gp);

            const double contact_angle_macro_gp = avg_contact_angle; //std::acos(inner_prod(wall_tangent,contact_vector_macro));
            double contact_angle_micro_gp = contact_angle_macro_gp;

            double zeta_effective = zeta;
            double contact_angle_equilibrium = theta_receding;
            if (contact_angle_macro_gp > contact_angle_equilibrium){
                if (contact_angle_macro_gp >= theta_advancing){
                    contact_angle_equilibrium = theta_advancing;
                } else {
                    contact_angle_equilibrium = contact_angle_macro_gp;
                    zeta_effective = 1.0e0*zeta;
                }
            }

            const double contact_velocity_gp = inner_prod(wall_tangent,velocity_gp);

            //const double reynolds_number = effective_density*std::abs(contact_velocity_gp)*element_size/effective_viscosity;
            const double capilary_number = effective_viscosity*contact_velocity_gp/coefficient;

            //KRATOS_INFO("two fluids NS") << "capilary_number= " << capilary_number << std::endl;
            //KRATOS_INFO("two fluids NS") << "reynolds_number= " << reynolds_number << std::endl;

            //KRATOS_INFO("two fluids NS") << "angle difference= " << std::abs(contact_angle_macro_gp - contact_angle_equilibrium)*180/PI << std::endl;

            if ( std::abs(contact_angle_macro_gp - contact_angle_equilibrium) < 6.0e-1 &&
                    capilary_number < 3.0e-1){
                const double cubic_contact_angle_micro_gp = std::pow(contact_angle_macro_gp, 3.0)
                    - 9*capilary_number*std::log(element_size/micro_length_scale);

                KRATOS_WARNING_IF("TwoFluidsNS", cubic_contact_angle_micro_gp < 0.0 ||
                                cubic_contact_angle_micro_gp > 31.0)
                            << "Hydrodynamics theory failed to estimate micro contact-angle (large slip velocity)." 
                            << std::endl;

                if (cubic_contact_angle_micro_gp >= 0.0 &&
                        cubic_contact_angle_micro_gp <= 31.0) //std::pow(PI, 3.0))
                    contact_angle_micro_gp = std::pow(cubic_contact_angle_micro_gp, 1.0/3.0);
                else if (cubic_contact_angle_micro_gp < 0.0)
                    contact_angle_micro_gp = 0.0; //contact_angle_equilibrium;
                else //if (cubic_contact_angle_micro_gp > 31.0)
                    contact_angle_micro_gp = PI; //contact_angle_equilibrium;

                // This relation is valid for contact_angle < 3PI/4 and vanishing Reynolds & Capillary numbers
            }

            // contact_angle_equilibrium = theta_receding;
            // if (contact_angle_micro_gp > contact_angle_equilibrium){
            //     if (contact_angle_micro_gp >= theta_advancing){
            //         contact_angle_equilibrium = theta_advancing;
            //         //KRATOS_WATCH("theta > theta_advancing")
            //     } else {
            //         contact_angle_equilibrium = contact_angle_micro_gp;
            //         zeta_effective = 1.0e5*zeta;
            //         /* for (unsigned int j = 0; j < NumNodes; j++){
            //             if ((*p_geom)[j].GetValue(IS_STRUCTURE) == 1.0){
            //                 #pragma omp critical
            //                 (*p_geom)[j].Fix(DISTANCE);
            //                 //KRATOS_WATCH(j)
            //             }
            //         } */
            //     }
            // } //else {
            //     //KRATOS_WATCH("theta < theta_receding")
            // //}

            const double coefficientS = coefficient*std::cos(contact_angle_equilibrium);

            //KRATOS_INFO("two fluids NS") << "element_size= " << element_size << std::endl;
            //KRATOS_INFO("two fluids NS") << "contact_angle_macro_gp= " << contact_angle_macro_gp << std::endl;
            //KRATOS_INFO("two fluids NS") << "contact_angle_micro_gp= " << contact_angle_micro_gp << std::endl;

            contact_vector_micro = std::cos(contact_angle_micro_gp)*wall_tangent +
                    std::sin(contact_angle_micro_gp)*wall_normal_gp;

            for (unsigned int i = 0; i < NumNodes; i++){
                for (unsigned int dimi = 0; dimi < NumDim; dimi++){
                    rhs[ i*(NumDim+1) + dimi ] -= coefficient*contact_vector_micro[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i);
                    rhs[ i*(NumDim+1) + dimi ] += coefficientS*wall_tangent[dimi]*(rCLWeights[i_cl])[clgp]*(rCLShapeFunctions[i_cl])(clgp,i); //Contac-line tangential force

                    for (unsigned int j = 0; j < NumNodes; j++){
                            lhs_dissipation( i*(NumDim+1) + dimi, j*(NumDim+1) + dimi) +=
                                zeta_effective * (rCLWeights[i_cl])[clgp] * (rCLShapeFunctions[i_cl])(clgp,j) * (rCLShapeFunctions[i_cl])(clgp,i);
                    }
                }
            }
            contact_velocity += (rCLWeights[i_cl])[clgp]*contact_velocity_gp;
            contact_angle_macro += (rCLWeights[i_cl])[clgp]*contact_angle_macro_gp;
            contact_angle_micro += (rCLWeights[i_cl])[clgp]*contact_angle_micro_gp;
        }

        noalias(rLHS) += lhs_dissipation;
        rhs -= prod(lhs_dissipation,tempU);

        contact_angle_macro /= weight_sum;
        this->SetValue(CONTACT_ANGLE, contact_angle_macro*180.0/PI);

        contact_angle_micro /= weight_sum;
        this->SetValue(CONTACT_ANGLE_MICRO, contact_angle_micro*180.0/PI);

        contact_velocity /= weight_sum;
        this->SetValue(CONTACT_VELOCITY, contact_velocity);

        for (unsigned int i=0; i < NumNodes; ++i){

            #pragma omp critical
            {
            (*p_geom)[i].FastGetSolutionStepValue(NORMAL_VECTOR) = normal_avg;
            (*p_geom)[i].FastGetSolutionStepValue(TANGENT_VECTOR) = wall_tangent;
            (*p_geom)[i].FastGetSolutionStepValue(CONTACT_VECTOR) = contact_vector_macro;
            }

        }
    }

    noalias(rRHS) += rhs;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CondenseEnrichment(
    const TElementData &rData,
    Matrix &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const MatrixType &rHtot,
    const MatrixType &rVtot,
    MatrixType &rKeeTot,
    const VectorType &rRHSeeTot)
{
    //const double min_area_ratio = 1e-7;

    // Compute positive side, negative side and total volumes
    /* double positive_volume = 0.0;
    double negative_volume = 0.0;
    for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); ++igauss_pos){
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
    }

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); ++igauss_neg){
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
    }
    const double Vol = positive_volume + negative_volume; */
    
    //We only enrich elements which are not almost empty/full
    //if (positive_volume / Vol > min_area_ratio && negative_volume / Vol > min_area_ratio) {

        //The pressure continuity (coefficient=0.0) / discontinuity (coefficient > 0.0) constraints are now handled in 
        //PressureDiscontinuity function defined above. So the following part of code should be commented!

/*         // Compute the maximum diagonal value in the enrichment stiffness matrix
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
            for (unsigned int j = i + 1; j < NumNodes; j++){
                const double dj = std::abs(rData.Distance[j]);
                // Check if the edge is cut, if it is, set the penalty constraint
                if (rData.Distance[i] * rData.Distance[j] < 0.0){
                    double sum_d = di + dj;
                    double Ni = dj / sum_d;
                    double Nj = di / sum_d;
                    double penalty_coeff = 0.01;//max_diag * 0.001; // h/BDFVector[0];
                    rKeeTot(i, i) += penalty_coeff * Ni * Ni;
                    rKeeTot(i, j) -= penalty_coeff * Ni * Nj;
                    rKeeTot(j, i) -= penalty_coeff * Nj * Ni;
                    rKeeTot(j, j) += penalty_coeff * Nj * Nj;
                }
            }
        } */

        // Enrichment condensation (add to LHS and RHS the enrichment contributions)
        double det;
        MatrixType inverse_diag(NumNodes, NumNodes);
        MathUtils<double>::InvertMatrix(rKeeTot, inverse_diag, det);

        const Matrix IDENT1 = prod(inverse_diag, rKeeTot);
        const Matrix IDENT2 = prod(rKeeTot, inverse_diag);

        bool check = true; 
        for (unsigned int col = 0; col < NumNodes; ++col)
        { 
            for (unsigned int row = 0; row < NumNodes; ++row)
            { 
                if ((row == col) && (IDENT1(row,col) < 0.999 || IDENT1(row,col) > 1.001)) 
                    check = false; 
                else if (row != col && abs(IDENT1(row,col)) > 0.001) 
                    check = false; 

                if ((row == col) && (IDENT2(row,col) < 0.999 || IDENT2(row,col) > 1.001)) 
                    check = false; 
                else if (row != col && abs(IDENT2(row,col)) > 0.001) 
                    check = false; 
            }  
        } 

        if (!check) 
        {
            KRATOS_INFO("Condensation, InvertMatrix") << "**************************************************" << std::endl;
            KRATOS_INFO("Condensation, InvertMatrix") << "******************** Not OK **********************" << std::endl;
            KRATOS_INFO("Condensation, InvertMatrix") << "**************************************************" << std::endl;
        }

        const Matrix tmp = prod(inverse_diag, rHtot);
        noalias(rLeftHandSideMatrix) -= prod(rVtot, tmp);

        const Vector tmp2 = prod(inverse_diag, rRHSeeTot);
        noalias(rRightHandSideVector) -= prod(rVtot, tmp2);

        // Reproducing the enriched pressure for any possible use.
        /* const unsigned int NumDim = NumNodes - 1;
        GeometryType::Pointer p_geom = this->pGetGeometry();
        Vector dU = ZeroVector(NumNodes*(NumDim + 1));

        for (unsigned int j = 0; j < NumNodes; j++){
            const Vector dV = (*p_geom)[j].FastGetSolutionStepValue(VELOCITY) - (*p_geom)[j].FastGetSolutionStepValue(VELOCITY_STAR);
            for (unsigned int dim = 0; dim < NumDim; dim++){
                dU[ j*(NumDim+1) + dim ] = dV(dim);
            }
            dU[ j*(NumDim+1) + NumDim ] = 
                (*p_geom)[j].FastGetSolutionStepValue(PRESSURE) - (*p_geom)[j].FastGetSolutionStepValue(PRESSURE_STAR);

            const Vector velocity = (*p_geom)[j].FastGetSolutionStepValue(VELOCITY);
            (*p_geom)[j].FastGetSolutionStepValue(VELOCITY_STAR) = velocity;
            const double pressure = (*p_geom)[j].FastGetSolutionStepValue(PRESSURE);
            (*p_geom)[j].FastGetSolutionStepValue(PRESSURE_STAR) = pressure;
        }

        const Vector tmpV = rRHSeeTot - prod(rHtot, dU);
        const Vector dPenr = prod( inverse_diag, tmpV );

        double Penr = this->GetValue(ENRICHED_PRESSURE_1) + dPenr(0);
        //KRATOS_INFO("Condensation, Enriched Pressure 1") << Penr << std::endl;
        this->SetValue(ENRICHED_PRESSURE_1, Penr);
        Penr = this->GetValue(ENRICHED_PRESSURE_2) + dPenr(1);
        //KRATOS_INFO("Condensation, Enriched Pressure 2") << Penr << std::endl;
        this->SetValue(ENRICHED_PRESSURE_2, Penr);
        Penr = this->GetValue(ENRICHED_PRESSURE_3) + dPenr(2);
        //KRATOS_INFO("Condensation, Enriched Pressure 3") << Penr << std::endl;
        this->SetValue(ENRICHED_PRESSURE_3, Penr);
        if (NumDim == 3){
            Penr = this->GetValue(ENRICHED_PRESSURE_1) + dPenr(3);
            //KRATOS_INFO("Condensation, Enriched Pressure 4") << Penr << std::endl;
            this->SetValue(ENRICHED_PRESSURE_4, Penr);
        } */
    //}
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
void TwoFluidNavierStokes<TElementData>::GetValueOnIntegrationPoints(   const Variable<double> &rVariable,
                                                                        std::vector<double> &rValues,
                                                                        const ProcessInfo &rCurrentProcessInfo )
{
    //KRATOS_INFO("GetValueOnIntegrationPoints") << "CALLED!" << std::endl;

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

        for (unsigned int i_gauss = 0; i_gauss < num_gauss; i_gauss++){

            const Matrix gp_DN_DX = DN_DX[i_gauss];
            double DVi_DXi = 0.0;

            for(unsigned int nnode = 0; nnode < NumNodes; nnode++){

                const array_1d<double,3> vel = rGeom[nnode].GetSolutionStepValue(VELOCITY);
                for(unsigned int ndim = 0; ndim < Dim; ndim++){
                    DVi_DXi += gp_DN_DX(nnode, ndim) * vel[ndim];
                }
            }
            rValues[i_gauss] = DVi_DXi;
        }
    }
    else if (/* this->Has( rVariable) && */ rVariable == CONTACT_ANGLE) {
        GeometryType::Pointer p_geom = this->pGetGeometry();
        //const auto& r_geometry = GetGeometry();
        const GeometryType::IntegrationPointsArrayType& integration_points = 
            p_geom->IntegrationPoints(GeometryData::GI_GAUSS_1); //p_geom->GetDefaultIntegrationMethod());

        //const auto& rGeom = this->GetGeometry();
        //const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        //const unsigned int num_gauss = IntegrationPoints.size();

        const std::size_t number_of_integration_points = integration_points.size();

        if ( rValues.size() != number_of_integration_points ) {
            rValues.resize( number_of_integration_points );
        }

        //KRATOS_INFO("CalculateOnIntegrationPoints") << "CALLED!" << std::endl;
        
        const double value = this->GetValue( CONTACT_ANGLE);
        for (std::size_t point_number = 0; point_number < number_of_integration_points; ++point_number) {
            rValues[point_number] = value;

            //if (value != 0.0){
            //    KRATOS_INFO("CalculateOnIntegrationPoints") << value << std::endl;
            //}
        }
    }
    else if (/* this->Has( rVariable) && */ rVariable == CONTACT_VELOCITY) {
        GeometryType::Pointer p_geom = this->pGetGeometry();
        //const auto& r_geometry = GetGeometry();
        const GeometryType::IntegrationPointsArrayType& integration_points = 
            p_geom->IntegrationPoints(GeometryData::GI_GAUSS_1); //p_geom->GetDefaultIntegrationMethod());

        //const auto& rGeom = this->GetGeometry();
        //const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        //const unsigned int num_gauss = IntegrationPoints.size();

        const std::size_t number_of_integration_points = integration_points.size();

        if ( rValues.size() != number_of_integration_points ) {
            rValues.resize( number_of_integration_points );
        }

        //KRATOS_INFO("CalculateOnIntegrationPoints") << "CALLED!" << std::endl;
        
        const double value = this->GetValue( CONTACT_VELOCITY);
        for (std::size_t point_number = 0; point_number < number_of_integration_points; ++point_number) {
            rValues[point_number] = value;

            //if (value != 0.0){
            //    KRATOS_INFO("CalculateOnIntegrationPoints") << value << std::endl;
            //}
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointLHSContributionCut(
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

    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;
    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             C(0,2)*DN(0,0);
const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
const double clhs3 =             pow(DN(0,0), 2);
const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 =             rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             pow(N[0], 2);
const double clhs9 =             rho*(DN(0,0)*clhs4 + DN(0,1)*clhs5);
const double clhs10 =             bdf0*rho;
const double clhs11 =             K_darcy*N[0];
const double clhs12 =             N[0]*clhs10;
const double clhs13 =             clhs11 + clhs12 + clhs9;
const double clhs14 =             1.0/(K_darcy + clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 =             1.0*clhs9;
const double clhs16 =             clhs14*clhs15;
const double clhs17 =             1.0*clhs11;
const double clhs18 =             clhs14*clhs17;
const double clhs19 =             K_darcy*clhs8 + N[0]*clhs9 + clhs10*clhs8 + clhs13*clhs16 - clhs13*clhs18;
const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
const double clhs21 =             C(1,2)*DN(0,1);
const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,0)*clhs7;
const double clhs24 =             DN(0,1)*clhs23;
const double clhs25 =             -N[0] + clhs16 - clhs18;
const double clhs26 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs27 =             C(0,2)*DN(1,0);
const double clhs28 =             C(2,2)*DN(1,1) + clhs27;
const double clhs29 =             DN(0,0)*DN(1,0);
const double clhs30 =             N[1]*clhs11 + N[1]*clhs12;
const double clhs31 =             clhs29*clhs7 + clhs30;
const double clhs32 =             rho*(DN(1,0)*clhs4 + DN(1,1)*clhs5);
const double clhs33 =             K_darcy*N[1];
const double clhs34 =             N[1]*clhs10;
const double clhs35 =             clhs32 + clhs33 + clhs34;
const double clhs36 =             N[0]*clhs32 + clhs16*clhs35 - clhs18*clhs35;
const double clhs37 =             C(0,1)*DN(1,1) + clhs27;
const double clhs38 =             C(1,2)*DN(1,1);
const double clhs39 =             C(2,2)*DN(1,0) + clhs38;
const double clhs40 =             DN(1,1)*clhs23;
const double clhs41 =             DN(0,0)*N[1];
const double clhs42 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs43 =             C(0,2)*DN(2,0);
const double clhs44 =             C(2,2)*DN(2,1) + clhs43;
const double clhs45 =             DN(0,0)*DN(2,0);
const double clhs46 =             N[2]*clhs11 + N[2]*clhs12;
const double clhs47 =             clhs45*clhs7 + clhs46;
const double clhs48 =             rho*(DN(2,0)*clhs4 + DN(2,1)*clhs5);
const double clhs49 =             K_darcy*N[2];
const double clhs50 =             N[2]*clhs10;
const double clhs51 =             clhs48 + clhs49 + clhs50;
const double clhs52 =             N[0]*clhs48 + clhs16*clhs51 - clhs18*clhs51;
const double clhs53 =             C(0,1)*DN(2,1) + clhs43;
const double clhs54 =             C(1,2)*DN(2,1);
const double clhs55 =             C(2,2)*DN(2,0) + clhs54;
const double clhs56 =             DN(2,1)*clhs23;
const double clhs57 =             DN(0,0)*N[2];
const double clhs58 =             C(0,1)*DN(0,0) + clhs21;
const double clhs59 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs60 =             pow(DN(0,1), 2);
const double clhs61 =             C(0,1)*DN(1,0) + clhs38;
const double clhs62 =             DN(0,1)*clhs7;
const double clhs63 =             DN(1,0)*clhs62;
const double clhs64 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs65 =             DN(0,1)*DN(1,1);
const double clhs66 =             clhs30 + clhs65*clhs7;
const double clhs67 =             DN(0,1)*N[1];
const double clhs68 =             C(0,1)*DN(2,0) + clhs54;
const double clhs69 =             DN(2,0)*clhs62;
const double clhs70 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs71 =             DN(0,1)*DN(2,1);
const double clhs72 =             clhs46 + clhs7*clhs71;
const double clhs73 =             DN(0,1)*N[2];
const double clhs74 =             N[0] + clhs14*(1.0*clhs12 + clhs15 + clhs17);
const double clhs75 =             1.0*clhs14;
const double clhs76 =             DN(1,0)*N[0];
const double clhs77 =             clhs35*clhs75;
const double clhs78 =             DN(1,1)*N[0];
const double clhs79 =             clhs75*(clhs29 + clhs65);
const double clhs80 =             DN(2,0)*N[0];
const double clhs81 =             clhs51*clhs75;
const double clhs82 =             DN(2,1)*N[0];
const double clhs83 =             clhs75*(clhs45 + clhs71);
const double clhs84 =             clhs32*clhs75;
const double clhs85 =             clhs33*clhs75;
const double clhs86 =             N[1]*clhs9 + clhs13*clhs84 - clhs13*clhs85;
const double clhs87 =             pow(DN(1,0), 2);
const double clhs88 =             pow(N[1], 2);
const double clhs89 =             K_darcy*clhs88 + N[1]*clhs32 + clhs10*clhs88 + clhs32*clhs77 - clhs33*clhs77;
const double clhs90 =             DN(1,0)*clhs7;
const double clhs91 =             DN(1,1)*clhs90;
const double clhs92 =             -N[1] + clhs84 - clhs85;
const double clhs93 =             DN(1,0)*DN(2,0);
const double clhs94 =             N[2]*clhs33 + N[2]*clhs34;
const double clhs95 =             clhs7*clhs93 + clhs94;
const double clhs96 =             N[1]*clhs48 + clhs32*clhs81 - clhs33*clhs81;
const double clhs97 =             DN(2,1)*clhs90;
const double clhs98 =             DN(1,0)*N[2];
const double clhs99 =             pow(DN(1,1), 2);
const double clhs100 =             DN(2,0)*clhs7;
const double clhs101 =             DN(1,1)*clhs100;
const double clhs102 =             DN(1,1)*DN(2,1);
const double clhs103 =             clhs102*clhs7 + clhs94;
const double clhs104 =             DN(1,1)*N[2];
const double clhs105 =             clhs13*clhs75;
const double clhs106 =             N[1] + clhs14*(1.0*clhs32 + 1.0*clhs33 + 1.0*clhs34);
const double clhs107 =             DN(2,0)*N[1];
const double clhs108 =             DN(2,1)*N[1];
const double clhs109 =             clhs75*(clhs102 + clhs93);
const double clhs110 =             N[2]*clhs9 + clhs105*clhs48 - clhs105*clhs49;
const double clhs111 =             clhs49*clhs75;
const double clhs112 =             clhs48*clhs75;
const double clhs113 =             N[2]*clhs32 + clhs48*clhs77 - clhs49*clhs77;
const double clhs114 =             pow(DN(2,0), 2);
const double clhs115 =             pow(N[2], 2);
const double clhs116 =             K_darcy*clhs115 + N[2]*clhs48 + clhs10*clhs115 + clhs48*clhs81 - clhs49*clhs81;
const double clhs117 =             DN(2,1)*clhs100;
const double clhs118 =             -N[2] - clhs111 + clhs112;
const double clhs119 =             pow(DN(2,1), 2);
const double clhs120 =             N[2] + clhs14*(1.0*clhs48 + 1.0*clhs49 + 1.0*clhs50);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
            lhs(0,2)=DN(0,0)*clhs25;
            lhs(0,3)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + clhs31 + clhs36;
            lhs(0,4)=DN(0,0)*clhs37 + DN(0,1)*clhs39 + clhs40;
            lhs(0,5)=DN(1,0)*clhs16 - DN(1,0)*clhs18 - clhs41;
            lhs(0,6)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + clhs47 + clhs52;
            lhs(0,7)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + clhs56;
            lhs(0,8)=DN(2,0)*clhs16 - DN(2,0)*clhs18 - clhs57;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs58 + clhs24;
            lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs59 + clhs19 + clhs60*clhs7;
            lhs(1,2)=DN(0,1)*clhs25;
            lhs(1,3)=DN(0,0)*clhs28 + DN(0,1)*clhs61 + clhs63;
            lhs(1,4)=DN(0,0)*clhs39 + DN(0,1)*clhs64 + clhs36 + clhs66;
            lhs(1,5)=DN(1,1)*clhs16 - DN(1,1)*clhs18 - clhs67;
            lhs(1,6)=DN(0,0)*clhs44 + DN(0,1)*clhs68 + clhs69;
            lhs(1,7)=DN(0,0)*clhs55 + DN(0,1)*clhs70 + clhs52 + clhs72;
            lhs(1,8)=DN(2,1)*clhs16 - DN(2,1)*clhs18 - clhs73;
            lhs(2,0)=DN(0,0)*clhs74;
            lhs(2,1)=DN(0,1)*clhs74;
            lhs(2,2)=clhs75*(clhs3 + clhs60);
            lhs(2,3)=DN(0,0)*clhs77 + clhs76;
            lhs(2,4)=DN(0,1)*clhs77 + clhs78;
            lhs(2,5)=clhs79;
            lhs(2,6)=DN(0,0)*clhs81 + clhs80;
            lhs(2,7)=DN(0,1)*clhs81 + clhs82;
            lhs(2,8)=clhs83;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs31 + clhs86;
            lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs63;
            lhs(3,2)=DN(0,0)*clhs84 - DN(0,0)*clhs85 - clhs76;
            lhs(3,3)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + clhs7*clhs87 + clhs89;
            lhs(3,4)=DN(1,0)*clhs37 + DN(1,1)*clhs39 + clhs91;
            lhs(3,5)=DN(1,0)*clhs92;
            lhs(3,6)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + clhs95 + clhs96;
            lhs(3,7)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + clhs97;
            lhs(3,8)=DN(2,0)*clhs84 - DN(2,0)*clhs85 - clhs98;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs58 + clhs40;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs59 + clhs66 + clhs86;
            lhs(4,2)=DN(0,1)*clhs84 - DN(0,1)*clhs85 - clhs78;
            lhs(4,3)=DN(1,0)*clhs28 + DN(1,1)*clhs61 + clhs91;
            lhs(4,4)=DN(1,0)*clhs39 + DN(1,1)*clhs64 + clhs7*clhs99 + clhs89;
            lhs(4,5)=DN(1,1)*clhs92;
            lhs(4,6)=DN(1,0)*clhs44 + DN(1,1)*clhs68 + clhs101;
            lhs(4,7)=DN(1,0)*clhs55 + DN(1,1)*clhs70 + clhs103 + clhs96;
            lhs(4,8)=DN(2,1)*clhs84 - DN(2,1)*clhs85 - clhs104;
            lhs(5,0)=DN(1,0)*clhs105 + clhs41;
            lhs(5,1)=DN(1,1)*clhs105 + clhs67;
            lhs(5,2)=clhs79;
            lhs(5,3)=DN(1,0)*clhs106;
            lhs(5,4)=DN(1,1)*clhs106;
            lhs(5,5)=clhs75*(clhs87 + clhs99);
            lhs(5,6)=DN(1,0)*clhs81 + clhs107;
            lhs(5,7)=DN(1,1)*clhs81 + clhs108;
            lhs(5,8)=clhs109;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs110 + clhs47;
            lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs69;
            lhs(6,2)=-DN(0,0)*clhs111 + DN(0,0)*clhs112 - clhs80;
            lhs(6,3)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + clhs113 + clhs95;
            lhs(6,4)=DN(2,0)*clhs37 + DN(2,1)*clhs39 + clhs101;
            lhs(6,5)=-DN(1,0)*clhs111 + DN(1,0)*clhs112 - clhs107;
            lhs(6,6)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + clhs114*clhs7 + clhs116;
            lhs(6,7)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + clhs117;
            lhs(6,8)=DN(2,0)*clhs118;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs58 + clhs56;
            lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs59 + clhs110 + clhs72;
            lhs(7,2)=-DN(0,1)*clhs111 + DN(0,1)*clhs112 - clhs82;
            lhs(7,3)=DN(2,0)*clhs28 + DN(2,1)*clhs61 + clhs97;
            lhs(7,4)=DN(2,0)*clhs39 + DN(2,1)*clhs64 + clhs103 + clhs113;
            lhs(7,5)=-DN(1,1)*clhs111 + DN(1,1)*clhs112 - clhs108;
            lhs(7,6)=DN(2,0)*clhs44 + DN(2,1)*clhs68 + clhs117;
            lhs(7,7)=DN(2,0)*clhs55 + DN(2,1)*clhs70 + clhs116 + clhs119*clhs7;
            lhs(7,8)=DN(2,1)*clhs118;
            lhs(8,0)=DN(2,0)*clhs105 + clhs57;
            lhs(8,1)=DN(2,1)*clhs105 + clhs73;
            lhs(8,2)=clhs83;
            lhs(8,3)=DN(2,0)*clhs77 + clhs98;
            lhs(8,4)=DN(2,1)*clhs77 + clhs104;
            lhs(8,5)=clhs109;
            lhs(8,6)=DN(2,0)*clhs120;
            lhs(8,7)=DN(2,1)*clhs120;
            lhs(8,8)=clhs75*(clhs114 + clhs119);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointLHSContributionCut(
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

    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;
    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             C(0,3)*DN(0,0);
const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 =             C(0,5)*DN(0,0);
const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs9*h/stab_c1 + mu;
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             rho*(DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8);
const double clhs13 =             bdf0*rho;
const double clhs14 =             K_darcy*N[0];
const double clhs15 =             N[0]*clhs13;
const double clhs16 =             clhs12 + clhs14 + clhs15;
const double clhs17 =             1.0/(K_darcy + clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 =             1.0*clhs12;
const double clhs19 =             clhs17*clhs18;
const double clhs20 =             1.0*clhs14;
const double clhs21 =             clhs17*clhs20;
const double clhs22 =             K_darcy*clhs11 + N[0]*clhs12 + clhs11*clhs13 + clhs16*clhs19 - clhs16*clhs21;
const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs24 =             C(1,3)*DN(0,1);
const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
const double clhs26 =             C(3,5)*DN(0,0);
const double clhs27 =             C(4,5)*DN(0,2);
const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
const double clhs29 =             DN(0,0)*clhs10;
const double clhs30 =             DN(0,1)*clhs29;
const double clhs31 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs32 =             C(3,4)*DN(0,1);
const double clhs33 =             C(2,3)*DN(0,2) + clhs26 + clhs32;
const double clhs34 =             C(2,5)*DN(0,2);
const double clhs35 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
const double clhs36 =             DN(0,2)*clhs29;
const double clhs37 =             -N[0] + clhs19 - clhs21;
const double clhs38 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs39 =             C(0,3)*DN(1,0);
const double clhs40 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs39;
const double clhs41 =             C(0,5)*DN(1,0);
const double clhs42 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs41;
const double clhs43 =             DN(0,0)*DN(1,0);
const double clhs44 =             N[1]*clhs14 + N[1]*clhs15;
const double clhs45 =             clhs10*clhs43 + clhs44;
const double clhs46 =             rho*(DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8);
const double clhs47 =             K_darcy*N[1];
const double clhs48 =             N[1]*clhs13;
const double clhs49 =             clhs46 + clhs47 + clhs48;
const double clhs50 =             N[0]*clhs46 + clhs19*clhs49 - clhs21*clhs49;
const double clhs51 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs39;
const double clhs52 =             C(1,3)*DN(1,1);
const double clhs53 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs52;
const double clhs54 =             C(3,5)*DN(1,0);
const double clhs55 =             C(4,5)*DN(1,2);
const double clhs56 =             C(1,5)*DN(1,1) + clhs54 + clhs55;
const double clhs57 =             DN(1,1)*clhs29;
const double clhs58 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs41;
const double clhs59 =             C(3,4)*DN(1,1);
const double clhs60 =             C(2,3)*DN(1,2) + clhs54 + clhs59;
const double clhs61 =             C(2,5)*DN(1,2);
const double clhs62 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs61;
const double clhs63 =             DN(1,2)*clhs29;
const double clhs64 =             DN(0,0)*N[1];
const double clhs65 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs66 =             C(0,3)*DN(2,0);
const double clhs67 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs66;
const double clhs68 =             C(0,5)*DN(2,0);
const double clhs69 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs68;
const double clhs70 =             DN(0,0)*DN(2,0);
const double clhs71 =             N[2]*clhs14 + N[2]*clhs15;
const double clhs72 =             clhs10*clhs70 + clhs71;
const double clhs73 =             rho*(DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8);
const double clhs74 =             K_darcy*N[2];
const double clhs75 =             N[2]*clhs13;
const double clhs76 =             clhs73 + clhs74 + clhs75;
const double clhs77 =             N[0]*clhs73 + clhs19*clhs76 - clhs21*clhs76;
const double clhs78 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs66;
const double clhs79 =             C(1,3)*DN(2,1);
const double clhs80 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs79;
const double clhs81 =             C(3,5)*DN(2,0);
const double clhs82 =             C(4,5)*DN(2,2);
const double clhs83 =             C(1,5)*DN(2,1) + clhs81 + clhs82;
const double clhs84 =             DN(2,1)*clhs29;
const double clhs85 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs68;
const double clhs86 =             C(3,4)*DN(2,1);
const double clhs87 =             C(2,3)*DN(2,2) + clhs81 + clhs86;
const double clhs88 =             C(2,5)*DN(2,2);
const double clhs89 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs88;
const double clhs90 =             DN(2,2)*clhs29;
const double clhs91 =             DN(0,0)*N[2];
const double clhs92 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs93 =             C(0,3)*DN(3,0);
const double clhs94 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs93;
const double clhs95 =             C(0,5)*DN(3,0);
const double clhs96 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs95;
const double clhs97 =             DN(0,0)*DN(3,0);
const double clhs98 =             N[3]*clhs14 + N[3]*clhs15;
const double clhs99 =             clhs10*clhs97 + clhs98;
const double clhs100 =             rho*(DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8);
const double clhs101 =             K_darcy*N[3];
const double clhs102 =             N[3]*clhs13;
const double clhs103 =             clhs100 + clhs101 + clhs102;
const double clhs104 =             N[0]*clhs100 + clhs103*clhs19 - clhs103*clhs21;
const double clhs105 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs93;
const double clhs106 =             C(1,3)*DN(3,1);
const double clhs107 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs106;
const double clhs108 =             C(3,5)*DN(3,0);
const double clhs109 =             C(4,5)*DN(3,2);
const double clhs110 =             C(1,5)*DN(3,1) + clhs108 + clhs109;
const double clhs111 =             DN(3,1)*clhs29;
const double clhs112 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs95;
const double clhs113 =             C(3,4)*DN(3,1);
const double clhs114 =             C(2,3)*DN(3,2) + clhs108 + clhs113;
const double clhs115 =             C(2,5)*DN(3,2);
const double clhs116 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs115;
const double clhs117 =             DN(3,2)*clhs29;
const double clhs118 =             DN(0,0)*N[3];
const double clhs119 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs120 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs121 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs122 =             C(1,4)*DN(0,1);
const double clhs123 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs122;
const double clhs124 =             pow(DN(0,1), 2);
const double clhs125 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs122;
const double clhs126 =             C(2,4)*DN(0,2);
const double clhs127 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs126;
const double clhs128 =             DN(0,1)*clhs10;
const double clhs129 =             DN(0,2)*clhs128;
const double clhs130 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs52;
const double clhs131 =             C(0,4)*DN(1,0) + clhs55 + clhs59;
const double clhs132 =             DN(1,0)*clhs128;
const double clhs133 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs134 =             C(1,4)*DN(1,1);
const double clhs135 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs134;
const double clhs136 =             DN(0,1)*DN(1,1);
const double clhs137 =             clhs10*clhs136;
const double clhs138 =             clhs44 + clhs50;
const double clhs139 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs134;
const double clhs140 =             C(2,4)*DN(1,2);
const double clhs141 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs140;
const double clhs142 =             DN(1,2)*clhs128;
const double clhs143 =             DN(0,1)*N[1];
const double clhs144 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs79;
const double clhs145 =             C(0,4)*DN(2,0) + clhs82 + clhs86;
const double clhs146 =             DN(2,0)*clhs128;
const double clhs147 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs148 =             C(1,4)*DN(2,1);
const double clhs149 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs148;
const double clhs150 =             DN(0,1)*DN(2,1);
const double clhs151 =             clhs10*clhs150;
const double clhs152 =             clhs71 + clhs77;
const double clhs153 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs148;
const double clhs154 =             C(2,4)*DN(2,2);
const double clhs155 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs154;
const double clhs156 =             DN(2,2)*clhs128;
const double clhs157 =             DN(0,1)*N[2];
const double clhs158 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs106;
const double clhs159 =             C(0,4)*DN(3,0) + clhs109 + clhs113;
const double clhs160 =             DN(3,0)*clhs128;
const double clhs161 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs162 =             C(1,4)*DN(3,1);
const double clhs163 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs162;
const double clhs164 =             DN(0,1)*DN(3,1);
const double clhs165 =             clhs10*clhs164;
const double clhs166 =             clhs104 + clhs98;
const double clhs167 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs162;
const double clhs168 =             C(2,4)*DN(3,2);
const double clhs169 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs168;
const double clhs170 =             DN(3,2)*clhs128;
const double clhs171 =             DN(0,1)*N[3];
const double clhs172 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs173 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs126;
const double clhs174 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs175 =             pow(DN(0,2), 2);
const double clhs176 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs61;
const double clhs177 =             DN(0,2)*clhs10;
const double clhs178 =             DN(1,0)*clhs177;
const double clhs179 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs140;
const double clhs180 =             DN(1,1)*clhs177;
const double clhs181 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs182 =             DN(0,2)*DN(1,2);
const double clhs183 =             clhs10*clhs182;
const double clhs184 =             DN(0,2)*N[1];
const double clhs185 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs88;
const double clhs186 =             DN(2,0)*clhs177;
const double clhs187 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs154;
const double clhs188 =             DN(2,1)*clhs177;
const double clhs189 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs190 =             DN(0,2)*DN(2,2);
const double clhs191 =             clhs10*clhs190;
const double clhs192 =             DN(0,2)*N[2];
const double clhs193 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs115;
const double clhs194 =             DN(3,0)*clhs177;
const double clhs195 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs168;
const double clhs196 =             DN(3,1)*clhs177;
const double clhs197 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs198 =             DN(0,2)*DN(3,2);
const double clhs199 =             clhs10*clhs198;
const double clhs200 =             DN(0,2)*N[3];
const double clhs201 =             N[0] + clhs17*(1.0*clhs15 + clhs18 + clhs20);
const double clhs202 =             1.0*clhs17;
const double clhs203 =             DN(1,0)*N[0];
const double clhs204 =             clhs202*clhs49;
const double clhs205 =             DN(1,1)*N[0];
const double clhs206 =             DN(1,2)*N[0];
const double clhs207 =             clhs202*(clhs136 + clhs182 + clhs43);
const double clhs208 =             DN(2,0)*N[0];
const double clhs209 =             clhs202*clhs76;
const double clhs210 =             DN(2,1)*N[0];
const double clhs211 =             DN(2,2)*N[0];
const double clhs212 =             clhs202*(clhs150 + clhs190 + clhs70);
const double clhs213 =             DN(3,0)*N[0];
const double clhs214 =             clhs103*clhs202;
const double clhs215 =             DN(3,1)*N[0];
const double clhs216 =             DN(3,2)*N[0];
const double clhs217 =             clhs202*(clhs164 + clhs198 + clhs97);
const double clhs218 =             clhs202*clhs46;
const double clhs219 =             clhs202*clhs47;
const double clhs220 =             N[1]*clhs12 + clhs16*clhs218 - clhs16*clhs219;
const double clhs221 =             pow(DN(1,0), 2);
const double clhs222 =             pow(N[1], 2);
const double clhs223 =             K_darcy*clhs222 + N[1]*clhs46 + clhs13*clhs222 + clhs204*clhs46 - clhs204*clhs47;
const double clhs224 =             DN(1,0)*clhs10;
const double clhs225 =             DN(1,1)*clhs224;
const double clhs226 =             DN(1,2)*clhs224;
const double clhs227 =             -N[1] + clhs218 - clhs219;
const double clhs228 =             DN(1,0)*DN(2,0);
const double clhs229 =             N[2]*clhs47 + N[2]*clhs48;
const double clhs230 =             clhs10*clhs228 + clhs229;
const double clhs231 =             N[1]*clhs73 + clhs209*clhs46 - clhs209*clhs47;
const double clhs232 =             DN(2,1)*clhs224;
const double clhs233 =             DN(2,2)*clhs224;
const double clhs234 =             DN(1,0)*N[2];
const double clhs235 =             DN(1,0)*DN(3,0);
const double clhs236 =             N[3]*clhs47 + N[3]*clhs48;
const double clhs237 =             clhs10*clhs235 + clhs236;
const double clhs238 =             N[1]*clhs100 + clhs214*clhs46 - clhs214*clhs47;
const double clhs239 =             DN(3,1)*clhs224;
const double clhs240 =             DN(3,2)*clhs224;
const double clhs241 =             DN(1,0)*N[3];
const double clhs242 =             clhs220 + clhs44;
const double clhs243 =             pow(DN(1,1), 2);
const double clhs244 =             DN(1,1)*clhs10;
const double clhs245 =             DN(1,2)*clhs244;
const double clhs246 =             DN(2,0)*clhs244;
const double clhs247 =             DN(1,1)*DN(2,1);
const double clhs248 =             clhs10*clhs247;
const double clhs249 =             clhs229 + clhs231;
const double clhs250 =             DN(2,2)*clhs244;
const double clhs251 =             DN(1,1)*N[2];
const double clhs252 =             DN(3,0)*clhs244;
const double clhs253 =             DN(1,1)*DN(3,1);
const double clhs254 =             clhs10*clhs253;
const double clhs255 =             clhs236 + clhs238;
const double clhs256 =             DN(3,2)*clhs244;
const double clhs257 =             DN(1,1)*N[3];
const double clhs258 =             pow(DN(1,2), 2);
const double clhs259 =             DN(1,2)*clhs10;
const double clhs260 =             DN(2,0)*clhs259;
const double clhs261 =             DN(2,1)*clhs259;
const double clhs262 =             DN(1,2)*DN(2,2);
const double clhs263 =             clhs10*clhs262;
const double clhs264 =             DN(1,2)*N[2];
const double clhs265 =             DN(3,0)*clhs259;
const double clhs266 =             DN(3,1)*clhs259;
const double clhs267 =             DN(1,2)*DN(3,2);
const double clhs268 =             clhs10*clhs267;
const double clhs269 =             DN(1,2)*N[3];
const double clhs270 =             clhs16*clhs202;
const double clhs271 =             N[1] + clhs17*(1.0*clhs46 + 1.0*clhs47 + 1.0*clhs48);
const double clhs272 =             DN(2,0)*N[1];
const double clhs273 =             DN(2,1)*N[1];
const double clhs274 =             DN(2,2)*N[1];
const double clhs275 =             clhs202*(clhs228 + clhs247 + clhs262);
const double clhs276 =             DN(3,0)*N[1];
const double clhs277 =             DN(3,1)*N[1];
const double clhs278 =             DN(3,2)*N[1];
const double clhs279 =             clhs202*(clhs235 + clhs253 + clhs267);
const double clhs280 =             N[2]*clhs12 + clhs270*clhs73 - clhs270*clhs74;
const double clhs281 =             clhs202*clhs74;
const double clhs282 =             clhs202*clhs73;
const double clhs283 =             N[2]*clhs46 + clhs204*clhs73 - clhs204*clhs74;
const double clhs284 =             pow(DN(2,0), 2);
const double clhs285 =             pow(N[2], 2);
const double clhs286 =             K_darcy*clhs285 + N[2]*clhs73 + clhs13*clhs285 + clhs209*clhs73 - clhs209*clhs74;
const double clhs287 =             DN(2,0)*clhs10;
const double clhs288 =             DN(2,1)*clhs287;
const double clhs289 =             DN(2,2)*clhs287;
const double clhs290 =             -N[2] - clhs281 + clhs282;
const double clhs291 =             DN(2,0)*DN(3,0);
const double clhs292 =             N[3]*clhs74 + N[3]*clhs75;
const double clhs293 =             clhs10*clhs291 + clhs292;
const double clhs294 =             N[2]*clhs100 + clhs214*clhs73 - clhs214*clhs74;
const double clhs295 =             DN(3,1)*clhs287;
const double clhs296 =             DN(3,2)*clhs287;
const double clhs297 =             DN(2,0)*N[3];
const double clhs298 =             clhs280 + clhs71;
const double clhs299 =             clhs229 + clhs283;
const double clhs300 =             pow(DN(2,1), 2);
const double clhs301 =             DN(2,1)*clhs10;
const double clhs302 =             DN(2,2)*clhs301;
const double clhs303 =             DN(3,0)*clhs301;
const double clhs304 =             DN(2,1)*DN(3,1);
const double clhs305 =             clhs10*clhs304;
const double clhs306 =             clhs292 + clhs294;
const double clhs307 =             DN(3,2)*clhs301;
const double clhs308 =             DN(2,1)*N[3];
const double clhs309 =             pow(DN(2,2), 2);
const double clhs310 =             DN(2,2)*clhs10;
const double clhs311 =             DN(3,0)*clhs310;
const double clhs312 =             DN(3,1)*clhs310;
const double clhs313 =             DN(2,2)*DN(3,2);
const double clhs314 =             clhs10*clhs313;
const double clhs315 =             DN(2,2)*N[3];
const double clhs316 =             N[2] + clhs17*(1.0*clhs73 + 1.0*clhs74 + 1.0*clhs75);
const double clhs317 =             DN(3,0)*N[2];
const double clhs318 =             DN(3,1)*N[2];
const double clhs319 =             DN(3,2)*N[2];
const double clhs320 =             clhs202*(clhs291 + clhs304 + clhs313);
const double clhs321 =             N[3]*clhs12 + clhs100*clhs270 - clhs101*clhs270;
const double clhs322 =             clhs101*clhs202;
const double clhs323 =             clhs100*clhs202;
const double clhs324 =             N[3]*clhs46 + clhs100*clhs204 - clhs101*clhs204;
const double clhs325 =             N[3]*clhs73 + clhs100*clhs209 - clhs101*clhs209;
const double clhs326 =             pow(DN(3,0), 2);
const double clhs327 =             pow(N[3], 2);
const double clhs328 =             K_darcy*clhs327 + N[3]*clhs100 + clhs100*clhs214 - clhs101*clhs214 + clhs13*clhs327;
const double clhs329 =             DN(3,0)*clhs10;
const double clhs330 =             DN(3,1)*clhs329;
const double clhs331 =             DN(3,2)*clhs329;
const double clhs332 =             -N[3] - clhs322 + clhs323;
const double clhs333 =             clhs321 + clhs98;
const double clhs334 =             clhs236 + clhs324;
const double clhs335 =             clhs292 + clhs325;
const double clhs336 =             pow(DN(3,1), 2);
const double clhs337 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs338 =             pow(DN(3,2), 2);
const double clhs339 =             N[3] + clhs17*(1.0*clhs100 + 1.0*clhs101 + 1.0*clhs102);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
            lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
            lhs(0,3)=DN(0,0)*clhs37;
            lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + DN(0,2)*clhs42 + clhs45 + clhs50;
            lhs(0,5)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + DN(0,2)*clhs56 + clhs57;
            lhs(0,6)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs62 + clhs63;
            lhs(0,7)=DN(1,0)*clhs19 - DN(1,0)*clhs21 - clhs64;
            lhs(0,8)=DN(0,0)*clhs65 + DN(0,1)*clhs67 + DN(0,2)*clhs69 + clhs72 + clhs77;
            lhs(0,9)=DN(0,0)*clhs78 + DN(0,1)*clhs80 + DN(0,2)*clhs83 + clhs84;
            lhs(0,10)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs89 + clhs90;
            lhs(0,11)=DN(2,0)*clhs19 - DN(2,0)*clhs21 - clhs91;
            lhs(0,12)=DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs96 + clhs104 + clhs99;
            lhs(0,13)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs110 + clhs111;
            lhs(0,14)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs117;
            lhs(0,15)=DN(3,0)*clhs19 - DN(3,0)*clhs21 - clhs118;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs119 + DN(0,2)*clhs120 + clhs30;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs10*clhs124 + clhs22;
            lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs125 + DN(0,2)*clhs127 + clhs129;
            lhs(1,3)=DN(0,1)*clhs37;
            lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs132;
            lhs(1,5)=DN(0,0)*clhs53 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs137 + clhs138;
            lhs(1,6)=DN(0,0)*clhs60 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142;
            lhs(1,7)=DN(1,1)*clhs19 - DN(1,1)*clhs21 - clhs143;
            lhs(1,8)=DN(0,0)*clhs67 + DN(0,1)*clhs144 + DN(0,2)*clhs145 + clhs146;
            lhs(1,9)=DN(0,0)*clhs80 + DN(0,1)*clhs147 + DN(0,2)*clhs149 + clhs151 + clhs152;
            lhs(1,10)=DN(0,0)*clhs87 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs156;
            lhs(1,11)=DN(2,1)*clhs19 - DN(2,1)*clhs21 - clhs157;
            lhs(1,12)=DN(0,0)*clhs94 + DN(0,1)*clhs158 + DN(0,2)*clhs159 + clhs160;
            lhs(1,13)=DN(0,0)*clhs107 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs165 + clhs166;
            lhs(1,14)=DN(0,0)*clhs114 + DN(0,1)*clhs167 + DN(0,2)*clhs169 + clhs170;
            lhs(1,15)=DN(3,1)*clhs19 - DN(3,1)*clhs21 - clhs171;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs120 + DN(0,2)*clhs172 + clhs36;
            lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs123 + DN(0,2)*clhs173 + clhs129;
            lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs127 + DN(0,2)*clhs174 + clhs10*clhs175 + clhs22;
            lhs(2,3)=DN(0,2)*clhs37;
            lhs(2,4)=DN(0,0)*clhs42 + DN(0,1)*clhs131 + DN(0,2)*clhs176 + clhs178;
            lhs(2,5)=DN(0,0)*clhs56 + DN(0,1)*clhs135 + DN(0,2)*clhs179 + clhs180;
            lhs(2,6)=DN(0,0)*clhs62 + DN(0,1)*clhs141 + DN(0,2)*clhs181 + clhs138 + clhs183;
            lhs(2,7)=DN(1,2)*clhs19 - DN(1,2)*clhs21 - clhs184;
            lhs(2,8)=DN(0,0)*clhs69 + DN(0,1)*clhs145 + DN(0,2)*clhs185 + clhs186;
            lhs(2,9)=DN(0,0)*clhs83 + DN(0,1)*clhs149 + DN(0,2)*clhs187 + clhs188;
            lhs(2,10)=DN(0,0)*clhs89 + DN(0,1)*clhs155 + DN(0,2)*clhs189 + clhs152 + clhs191;
            lhs(2,11)=DN(2,2)*clhs19 - DN(2,2)*clhs21 - clhs192;
            lhs(2,12)=DN(0,0)*clhs96 + DN(0,1)*clhs159 + DN(0,2)*clhs193 + clhs194;
            lhs(2,13)=DN(0,0)*clhs110 + DN(0,1)*clhs163 + DN(0,2)*clhs195 + clhs196;
            lhs(2,14)=DN(0,0)*clhs116 + DN(0,1)*clhs169 + DN(0,2)*clhs197 + clhs166 + clhs199;
            lhs(2,15)=DN(3,2)*clhs19 - DN(3,2)*clhs21 - clhs200;
            lhs(3,0)=DN(0,0)*clhs201;
            lhs(3,1)=DN(0,1)*clhs201;
            lhs(3,2)=DN(0,2)*clhs201;
            lhs(3,3)=clhs202*(clhs124 + clhs175 + clhs5);
            lhs(3,4)=DN(0,0)*clhs204 + clhs203;
            lhs(3,5)=DN(0,1)*clhs204 + clhs205;
            lhs(3,6)=DN(0,2)*clhs204 + clhs206;
            lhs(3,7)=clhs207;
            lhs(3,8)=DN(0,0)*clhs209 + clhs208;
            lhs(3,9)=DN(0,1)*clhs209 + clhs210;
            lhs(3,10)=DN(0,2)*clhs209 + clhs211;
            lhs(3,11)=clhs212;
            lhs(3,12)=DN(0,0)*clhs214 + clhs213;
            lhs(3,13)=DN(0,1)*clhs214 + clhs215;
            lhs(3,14)=DN(0,2)*clhs214 + clhs216;
            lhs(3,15)=clhs217;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs220 + clhs45;
            lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs132;
            lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs178;
            lhs(4,3)=DN(0,0)*clhs218 - DN(0,0)*clhs219 - clhs203;
            lhs(4,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + DN(1,2)*clhs42 + clhs10*clhs221 + clhs223;
            lhs(4,5)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + DN(1,2)*clhs56 + clhs225;
            lhs(4,6)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs62 + clhs226;
            lhs(4,7)=DN(1,0)*clhs227;
            lhs(4,8)=DN(1,0)*clhs65 + DN(1,1)*clhs67 + DN(1,2)*clhs69 + clhs230 + clhs231;
            lhs(4,9)=DN(1,0)*clhs78 + DN(1,1)*clhs80 + DN(1,2)*clhs83 + clhs232;
            lhs(4,10)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs89 + clhs233;
            lhs(4,11)=DN(2,0)*clhs218 - DN(2,0)*clhs219 - clhs234;
            lhs(4,12)=DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs96 + clhs237 + clhs238;
            lhs(4,13)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs110 + clhs239;
            lhs(4,14)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs240;
            lhs(4,15)=DN(3,0)*clhs218 - DN(3,0)*clhs219 - clhs241;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs119 + DN(1,2)*clhs120 + clhs57;
            lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs137 + clhs242;
            lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs125 + DN(1,2)*clhs127 + clhs180;
            lhs(5,3)=DN(0,1)*clhs218 - DN(0,1)*clhs219 - clhs205;
            lhs(5,4)=DN(1,0)*clhs40 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs225;
            lhs(5,5)=DN(1,0)*clhs53 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs10*clhs243 + clhs223;
            lhs(5,6)=DN(1,0)*clhs60 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs245;
            lhs(5,7)=DN(1,1)*clhs227;
            lhs(5,8)=DN(1,0)*clhs67 + DN(1,1)*clhs144 + DN(1,2)*clhs145 + clhs246;
            lhs(5,9)=DN(1,0)*clhs80 + DN(1,1)*clhs147 + DN(1,2)*clhs149 + clhs248 + clhs249;
            lhs(5,10)=DN(1,0)*clhs87 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs250;
            lhs(5,11)=DN(2,1)*clhs218 - DN(2,1)*clhs219 - clhs251;
            lhs(5,12)=DN(1,0)*clhs94 + DN(1,1)*clhs158 + DN(1,2)*clhs159 + clhs252;
            lhs(5,13)=DN(1,0)*clhs107 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs254 + clhs255;
            lhs(5,14)=DN(1,0)*clhs114 + DN(1,1)*clhs167 + DN(1,2)*clhs169 + clhs256;
            lhs(5,15)=DN(3,1)*clhs218 - DN(3,1)*clhs219 - clhs257;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs120 + DN(1,2)*clhs172 + clhs63;
            lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs123 + DN(1,2)*clhs173 + clhs142;
            lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs127 + DN(1,2)*clhs174 + clhs183 + clhs242;
            lhs(6,3)=DN(0,2)*clhs218 - DN(0,2)*clhs219 - clhs206;
            lhs(6,4)=DN(1,0)*clhs42 + DN(1,1)*clhs131 + DN(1,2)*clhs176 + clhs226;
            lhs(6,5)=DN(1,0)*clhs56 + DN(1,1)*clhs135 + DN(1,2)*clhs179 + clhs245;
            lhs(6,6)=DN(1,0)*clhs62 + DN(1,1)*clhs141 + DN(1,2)*clhs181 + clhs10*clhs258 + clhs223;
            lhs(6,7)=DN(1,2)*clhs227;
            lhs(6,8)=DN(1,0)*clhs69 + DN(1,1)*clhs145 + DN(1,2)*clhs185 + clhs260;
            lhs(6,9)=DN(1,0)*clhs83 + DN(1,1)*clhs149 + DN(1,2)*clhs187 + clhs261;
            lhs(6,10)=DN(1,0)*clhs89 + DN(1,1)*clhs155 + DN(1,2)*clhs189 + clhs249 + clhs263;
            lhs(6,11)=DN(2,2)*clhs218 - DN(2,2)*clhs219 - clhs264;
            lhs(6,12)=DN(1,0)*clhs96 + DN(1,1)*clhs159 + DN(1,2)*clhs193 + clhs265;
            lhs(6,13)=DN(1,0)*clhs110 + DN(1,1)*clhs163 + DN(1,2)*clhs195 + clhs266;
            lhs(6,14)=DN(1,0)*clhs116 + DN(1,1)*clhs169 + DN(1,2)*clhs197 + clhs255 + clhs268;
            lhs(6,15)=DN(3,2)*clhs218 - DN(3,2)*clhs219 - clhs269;
            lhs(7,0)=DN(1,0)*clhs270 + clhs64;
            lhs(7,1)=DN(1,1)*clhs270 + clhs143;
            lhs(7,2)=DN(1,2)*clhs270 + clhs184;
            lhs(7,3)=clhs207;
            lhs(7,4)=DN(1,0)*clhs271;
            lhs(7,5)=DN(1,1)*clhs271;
            lhs(7,6)=DN(1,2)*clhs271;
            lhs(7,7)=clhs202*(clhs221 + clhs243 + clhs258);
            lhs(7,8)=DN(1,0)*clhs209 + clhs272;
            lhs(7,9)=DN(1,1)*clhs209 + clhs273;
            lhs(7,10)=DN(1,2)*clhs209 + clhs274;
            lhs(7,11)=clhs275;
            lhs(7,12)=DN(1,0)*clhs214 + clhs276;
            lhs(7,13)=DN(1,1)*clhs214 + clhs277;
            lhs(7,14)=DN(1,2)*clhs214 + clhs278;
            lhs(7,15)=clhs279;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs280 + clhs72;
            lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs146;
            lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs186;
            lhs(8,3)=-DN(0,0)*clhs281 + DN(0,0)*clhs282 - clhs208;
            lhs(8,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + DN(2,2)*clhs42 + clhs230 + clhs283;
            lhs(8,5)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + DN(2,2)*clhs56 + clhs246;
            lhs(8,6)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs62 + clhs260;
            lhs(8,7)=-DN(1,0)*clhs281 + DN(1,0)*clhs282 - clhs272;
            lhs(8,8)=DN(2,0)*clhs65 + DN(2,1)*clhs67 + DN(2,2)*clhs69 + clhs10*clhs284 + clhs286;
            lhs(8,9)=DN(2,0)*clhs78 + DN(2,1)*clhs80 + DN(2,2)*clhs83 + clhs288;
            lhs(8,10)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs89 + clhs289;
            lhs(8,11)=DN(2,0)*clhs290;
            lhs(8,12)=DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs96 + clhs293 + clhs294;
            lhs(8,13)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs110 + clhs295;
            lhs(8,14)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs296;
            lhs(8,15)=-DN(3,0)*clhs281 + DN(3,0)*clhs282 - clhs297;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs119 + DN(2,2)*clhs120 + clhs84;
            lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs151 + clhs298;
            lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs125 + DN(2,2)*clhs127 + clhs188;
            lhs(9,3)=-DN(0,1)*clhs281 + DN(0,1)*clhs282 - clhs210;
            lhs(9,4)=DN(2,0)*clhs40 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs232;
            lhs(9,5)=DN(2,0)*clhs53 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs248 + clhs299;
            lhs(9,6)=DN(2,0)*clhs60 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs261;
            lhs(9,7)=-DN(1,1)*clhs281 + DN(1,1)*clhs282 - clhs273;
            lhs(9,8)=DN(2,0)*clhs67 + DN(2,1)*clhs144 + DN(2,2)*clhs145 + clhs288;
            lhs(9,9)=DN(2,0)*clhs80 + DN(2,1)*clhs147 + DN(2,2)*clhs149 + clhs10*clhs300 + clhs286;
            lhs(9,10)=DN(2,0)*clhs87 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs302;
            lhs(9,11)=DN(2,1)*clhs290;
            lhs(9,12)=DN(2,0)*clhs94 + DN(2,1)*clhs158 + DN(2,2)*clhs159 + clhs303;
            lhs(9,13)=DN(2,0)*clhs107 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs305 + clhs306;
            lhs(9,14)=DN(2,0)*clhs114 + DN(2,1)*clhs167 + DN(2,2)*clhs169 + clhs307;
            lhs(9,15)=-DN(3,1)*clhs281 + DN(3,1)*clhs282 - clhs308;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs120 + DN(2,2)*clhs172 + clhs90;
            lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs123 + DN(2,2)*clhs173 + clhs156;
            lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs127 + DN(2,2)*clhs174 + clhs191 + clhs298;
            lhs(10,3)=-DN(0,2)*clhs281 + DN(0,2)*clhs282 - clhs211;
            lhs(10,4)=DN(2,0)*clhs42 + DN(2,1)*clhs131 + DN(2,2)*clhs176 + clhs233;
            lhs(10,5)=DN(2,0)*clhs56 + DN(2,1)*clhs135 + DN(2,2)*clhs179 + clhs250;
            lhs(10,6)=DN(2,0)*clhs62 + DN(2,1)*clhs141 + DN(2,2)*clhs181 + clhs263 + clhs299;
            lhs(10,7)=-DN(1,2)*clhs281 + DN(1,2)*clhs282 - clhs274;
            lhs(10,8)=DN(2,0)*clhs69 + DN(2,1)*clhs145 + DN(2,2)*clhs185 + clhs289;
            lhs(10,9)=DN(2,0)*clhs83 + DN(2,1)*clhs149 + DN(2,2)*clhs187 + clhs302;
            lhs(10,10)=DN(2,0)*clhs89 + DN(2,1)*clhs155 + DN(2,2)*clhs189 + clhs10*clhs309 + clhs286;
            lhs(10,11)=DN(2,2)*clhs290;
            lhs(10,12)=DN(2,0)*clhs96 + DN(2,1)*clhs159 + DN(2,2)*clhs193 + clhs311;
            lhs(10,13)=DN(2,0)*clhs110 + DN(2,1)*clhs163 + DN(2,2)*clhs195 + clhs312;
            lhs(10,14)=DN(2,0)*clhs116 + DN(2,1)*clhs169 + DN(2,2)*clhs197 + clhs306 + clhs314;
            lhs(10,15)=-DN(3,2)*clhs281 + DN(3,2)*clhs282 - clhs315;
            lhs(11,0)=DN(2,0)*clhs270 + clhs91;
            lhs(11,1)=DN(2,1)*clhs270 + clhs157;
            lhs(11,2)=DN(2,2)*clhs270 + clhs192;
            lhs(11,3)=clhs212;
            lhs(11,4)=DN(2,0)*clhs204 + clhs234;
            lhs(11,5)=DN(2,1)*clhs204 + clhs251;
            lhs(11,6)=DN(2,2)*clhs204 + clhs264;
            lhs(11,7)=clhs275;
            lhs(11,8)=DN(2,0)*clhs316;
            lhs(11,9)=DN(2,1)*clhs316;
            lhs(11,10)=DN(2,2)*clhs316;
            lhs(11,11)=clhs202*(clhs284 + clhs300 + clhs309);
            lhs(11,12)=DN(2,0)*clhs214 + clhs317;
            lhs(11,13)=DN(2,1)*clhs214 + clhs318;
            lhs(11,14)=DN(2,2)*clhs214 + clhs319;
            lhs(11,15)=clhs320;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs321 + clhs99;
            lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs160;
            lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs194;
            lhs(12,3)=-DN(0,0)*clhs322 + DN(0,0)*clhs323 - clhs213;
            lhs(12,4)=DN(3,0)*clhs38 + DN(3,1)*clhs40 + DN(3,2)*clhs42 + clhs237 + clhs324;
            lhs(12,5)=DN(3,0)*clhs51 + DN(3,1)*clhs53 + DN(3,2)*clhs56 + clhs252;
            lhs(12,6)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs62 + clhs265;
            lhs(12,7)=-DN(1,0)*clhs322 + DN(1,0)*clhs323 - clhs276;
            lhs(12,8)=DN(3,0)*clhs65 + DN(3,1)*clhs67 + DN(3,2)*clhs69 + clhs293 + clhs325;
            lhs(12,9)=DN(3,0)*clhs78 + DN(3,1)*clhs80 + DN(3,2)*clhs83 + clhs303;
            lhs(12,10)=DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs89 + clhs311;
            lhs(12,11)=-DN(2,0)*clhs322 + DN(2,0)*clhs323 - clhs317;
            lhs(12,12)=DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs96 + clhs10*clhs326 + clhs328;
            lhs(12,13)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs110 + clhs330;
            lhs(12,14)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs331;
            lhs(12,15)=DN(3,0)*clhs332;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs119 + DN(3,2)*clhs120 + clhs111;
            lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs165 + clhs333;
            lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs125 + DN(3,2)*clhs127 + clhs196;
            lhs(13,3)=-DN(0,1)*clhs322 + DN(0,1)*clhs323 - clhs215;
            lhs(13,4)=DN(3,0)*clhs40 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs239;
            lhs(13,5)=DN(3,0)*clhs53 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs254 + clhs334;
            lhs(13,6)=DN(3,0)*clhs60 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs266;
            lhs(13,7)=-DN(1,1)*clhs322 + DN(1,1)*clhs323 - clhs277;
            lhs(13,8)=DN(3,0)*clhs67 + DN(3,1)*clhs144 + DN(3,2)*clhs145 + clhs295;
            lhs(13,9)=DN(3,0)*clhs80 + DN(3,1)*clhs147 + DN(3,2)*clhs149 + clhs305 + clhs335;
            lhs(13,10)=DN(3,0)*clhs87 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs312;
            lhs(13,11)=-DN(2,1)*clhs322 + DN(2,1)*clhs323 - clhs318;
            lhs(13,12)=DN(3,0)*clhs94 + DN(3,1)*clhs158 + DN(3,2)*clhs159 + clhs330;
            lhs(13,13)=DN(3,0)*clhs107 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs10*clhs336 + clhs328;
            lhs(13,14)=DN(3,0)*clhs114 + DN(3,1)*clhs167 + DN(3,2)*clhs169 + clhs337;
            lhs(13,15)=DN(3,1)*clhs332;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs120 + DN(3,2)*clhs172 + clhs117;
            lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs123 + DN(3,2)*clhs173 + clhs170;
            lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs127 + DN(3,2)*clhs174 + clhs199 + clhs333;
            lhs(14,3)=-DN(0,2)*clhs322 + DN(0,2)*clhs323 - clhs216;
            lhs(14,4)=DN(3,0)*clhs42 + DN(3,1)*clhs131 + DN(3,2)*clhs176 + clhs240;
            lhs(14,5)=DN(3,0)*clhs56 + DN(3,1)*clhs135 + DN(3,2)*clhs179 + clhs256;
            lhs(14,6)=DN(3,0)*clhs62 + DN(3,1)*clhs141 + DN(3,2)*clhs181 + clhs268 + clhs334;
            lhs(14,7)=-DN(1,2)*clhs322 + DN(1,2)*clhs323 - clhs278;
            lhs(14,8)=DN(3,0)*clhs69 + DN(3,1)*clhs145 + DN(3,2)*clhs185 + clhs296;
            lhs(14,9)=DN(3,0)*clhs83 + DN(3,1)*clhs149 + DN(3,2)*clhs187 + clhs307;
            lhs(14,10)=DN(3,0)*clhs89 + DN(3,1)*clhs155 + DN(3,2)*clhs189 + clhs314 + clhs335;
            lhs(14,11)=-DN(2,2)*clhs322 + DN(2,2)*clhs323 - clhs319;
            lhs(14,12)=DN(3,0)*clhs96 + DN(3,1)*clhs159 + DN(3,2)*clhs193 + clhs331;
            lhs(14,13)=DN(3,0)*clhs110 + DN(3,1)*clhs163 + DN(3,2)*clhs195 + clhs337;
            lhs(14,14)=DN(3,0)*clhs116 + DN(3,1)*clhs169 + DN(3,2)*clhs197 + clhs10*clhs338 + clhs328;
            lhs(14,15)=DN(3,2)*clhs332;
            lhs(15,0)=DN(3,0)*clhs270 + clhs118;
            lhs(15,1)=DN(3,1)*clhs270 + clhs171;
            lhs(15,2)=DN(3,2)*clhs270 + clhs200;
            lhs(15,3)=clhs217;
            lhs(15,4)=DN(3,0)*clhs204 + clhs241;
            lhs(15,5)=DN(3,1)*clhs204 + clhs257;
            lhs(15,6)=DN(3,2)*clhs204 + clhs269;
            lhs(15,7)=clhs279;
            lhs(15,8)=DN(3,0)*clhs209 + clhs297;
            lhs(15,9)=DN(3,1)*clhs209 + clhs308;
            lhs(15,10)=DN(3,2)*clhs209 + clhs315;
            lhs(15,11)=clhs320;
            lhs(15,12)=DN(3,0)*clhs339;
            lhs(15,13)=DN(3,1)*clhs339;
            lhs(15,14)=DN(3,2)*clhs339;
            lhs(15,15)=clhs202*(clhs326 + clhs336 + clhs338);


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointRHSContributionCut(
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

    auto &rhs = rData.rhs;

    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;
    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             Nenr[0]*penr[0] + Nenr[1]*penr[1] + Nenr[2]*penr[2];
const double crhs2 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs3 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs5 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs8 =             rho*(crhs5*crhs6 + crhs7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs9 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs10 =             crhs5 + crhs9;
const double crhs11 =             rho*stab_c2*sqrt(pow(crhs6, 2) + pow(crhs7, 2));
const double crhs12 =             crhs10*(crhs11*h/stab_c1 + mu);
const double crhs13 =             1.0/(K_darcy + crhs11/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs14 =             crhs13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - crhs2 + crhs3 + crhs4 + crhs8);
const double crhs15 =             K_darcy*N[0];
const double crhs16 =             rho*(DN(0,0)*crhs6 + DN(0,1)*crhs7);
const double crhs17 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs18 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crhs19 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs20 =             rho*(crhs6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs7*crhs9);
const double crhs21 =             crhs13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - crhs17 + crhs18 + crhs19 + crhs20);
const double crhs22 =             K_darcy*N[1];
const double crhs23 =             rho*(DN(1,0)*crhs6 + DN(1,1)*crhs7);
const double crhs24 =             K_darcy*N[2];
const double crhs25 =             rho*(DN(2,0)*crhs6 + DN(2,1)*crhs7);
            rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs1 - DN(0,0)*crhs12 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs4 - N[0]*crhs8 + crhs14*crhs15 - crhs14*crhs16;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs1 - DN(0,1)*crhs12 - DN(0,1)*stress[1] + N[0]*crhs17 - N[0]*crhs18 - N[0]*crhs19 - N[0]*crhs20 + crhs15*crhs21 - crhs16*crhs21;
            rhs[2]=-DN(0,0)*crhs14 - DN(0,1)*crhs21 - N[0]*crhs10;
            rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs1 - DN(1,0)*crhs12 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs4 - N[1]*crhs8 + crhs14*crhs22 - crhs14*crhs23;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs1 - DN(1,1)*crhs12 - DN(1,1)*stress[1] + N[1]*crhs17 - N[1]*crhs18 - N[1]*crhs19 - N[1]*crhs20 + crhs21*crhs22 - crhs21*crhs23;
            rhs[5]=-DN(1,0)*crhs14 - DN(1,1)*crhs21 - N[1]*crhs10;
            rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs1 - DN(2,0)*crhs12 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs4 - N[2]*crhs8 + crhs14*crhs24 - crhs14*crhs25;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs1 - DN(2,1)*crhs12 - DN(2,1)*stress[1] + N[2]*crhs17 - N[2]*crhs18 - N[2]*crhs19 - N[2]*crhs20 + crhs21*crhs24 - crhs21*crhs25;
            rhs[8]=-DN(2,0)*crhs14 - DN(2,1)*crhs21 - N[2]*crhs10;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointRHSContributionCut(
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

    auto &rhs = rData.rhs;

    const auto &Nenr = rData.Nenr;
    const auto &DNenr = rData.DN_DXenr;
    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             Nenr[0]*penr[0] + Nenr[1]*penr[1] + Nenr[2]*penr[2] + Nenr[3]*penr[3];
const double crhs2 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs3 =             K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crhs4 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs5 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs9 =             rho*(crhs5*crhs6 + crhs7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs8*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs10 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs11 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs12 =             crhs10 + crhs11 + crhs5;
const double crhs13 =             rho*stab_c2*sqrt(pow(crhs6, 2) + pow(crhs7, 2) + pow(crhs8, 2));
const double crhs14 =             crhs12*(crhs13*h/stab_c1 + mu);
const double crhs15 =             1.0/(K_darcy + crhs13/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 =             crhs15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - crhs2 + crhs3 + crhs4 + crhs9);
const double crhs17 =             K_darcy*N[0];
const double crhs18 =             rho*(DN(0,0)*crhs6 + DN(0,1)*crhs7 + DN(0,2)*crhs8);
const double crhs19 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs20 =             K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crhs21 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs22 =             rho*(crhs10*crhs7 + crhs6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs8*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs23 =             crhs15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - crhs19 + crhs20 + crhs21 + crhs22);
const double crhs24 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs25 =             K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crhs26 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs27 =             rho*(crhs11*crhs8 + crhs6*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs7*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs28 =             crhs15*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - crhs24 + crhs25 + crhs26 + crhs27);
const double crhs29 =             K_darcy*N[1];
const double crhs30 =             rho*(DN(1,0)*crhs6 + DN(1,1)*crhs7 + DN(1,2)*crhs8);
const double crhs31 =             K_darcy*N[2];
const double crhs32 =             rho*(DN(2,0)*crhs6 + DN(2,1)*crhs7 + DN(2,2)*crhs8);
const double crhs33 =             K_darcy*N[3];
const double crhs34 =             rho*(DN(3,0)*crhs6 + DN(3,1)*crhs7 + DN(3,2)*crhs8);
            rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs1 - DN(0,0)*crhs14 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs2 - N[0]*crhs3 - N[0]*crhs4 - N[0]*crhs9 + crhs16*crhs17 - crhs16*crhs18;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 + DN(0,1)*crhs1 - DN(0,1)*crhs14 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs19 - N[0]*crhs20 - N[0]*crhs21 - N[0]*crhs22 + crhs17*crhs23 - crhs18*crhs23;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 + DN(0,2)*crhs1 - DN(0,2)*crhs14 - DN(0,2)*stress[2] + N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs27 + crhs17*crhs28 - crhs18*crhs28;
            rhs[3]=-DN(0,0)*crhs16 - DN(0,1)*crhs23 - DN(0,2)*crhs28 - N[0]*crhs12;
            rhs[4]=DN(1,0)*crhs0 + DN(1,0)*crhs1 - DN(1,0)*crhs14 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs2 - N[1]*crhs3 - N[1]*crhs4 - N[1]*crhs9 + crhs16*crhs29 - crhs16*crhs30;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 + DN(1,1)*crhs1 - DN(1,1)*crhs14 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs19 - N[1]*crhs20 - N[1]*crhs21 - N[1]*crhs22 + crhs23*crhs29 - crhs23*crhs30;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 + DN(1,2)*crhs1 - DN(1,2)*crhs14 - DN(1,2)*stress[2] + N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs27 + crhs28*crhs29 - crhs28*crhs30;
            rhs[7]=-DN(1,0)*crhs16 - DN(1,1)*crhs23 - DN(1,2)*crhs28 - N[1]*crhs12;
            rhs[8]=DN(2,0)*crhs0 + DN(2,0)*crhs1 - DN(2,0)*crhs14 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs2 - N[2]*crhs3 - N[2]*crhs4 - N[2]*crhs9 + crhs16*crhs31 - crhs16*crhs32;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 + DN(2,1)*crhs1 - DN(2,1)*crhs14 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs19 - N[2]*crhs20 - N[2]*crhs21 - N[2]*crhs22 + crhs23*crhs31 - crhs23*crhs32;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 + DN(2,2)*crhs1 - DN(2,2)*crhs14 - DN(2,2)*stress[2] + N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs27 + crhs28*crhs31 - crhs28*crhs32;
            rhs[11]=-DN(2,0)*crhs16 - DN(2,1)*crhs23 - DN(2,2)*crhs28 - N[2]*crhs12;
            rhs[12]=DN(3,0)*crhs0 + DN(3,0)*crhs1 - DN(3,0)*crhs14 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs2 - N[3]*crhs3 - N[3]*crhs4 - N[3]*crhs9 + crhs16*crhs33 - crhs16*crhs34;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 + DN(3,1)*crhs1 - DN(3,1)*crhs14 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs19 - N[3]*crhs20 - N[3]*crhs21 - N[3]*crhs22 + crhs23*crhs33 - crhs23*crhs34;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 + DN(3,2)*crhs1 - DN(3,2)*crhs14 - DN(3,2)*stress[2] + N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs27 + crhs28*crhs33 - crhs28*crhs34;
            rhs[15]=-DN(3,0)*crhs16 - DN(3,1)*crhs23 - DN(3,2)*crhs28 - N[3]*crhs12;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointEnrichmentContributionsCut(
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

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
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


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH3 =             cH2*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0));
const double cH4 =             cH2*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0));
const double cH5 =             cH2*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0));
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


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee2 =             crhs_ee0 + crhs_ee1;
const double crhs_ee3 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee4 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee5 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee3, 2) + pow(crhs_ee4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee6 =             crhs_ee5*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*crhs_ee3 + crhs_ee4*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
const double crhs_ee7 =             crhs_ee5*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee1*crhs_ee4 + crhs_ee3*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee6 - DNenr(0,1)*crhs_ee7 - Nenr[0]*crhs_ee2;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee6 - DNenr(1,1)*crhs_ee7 - Nenr[1]*crhs_ee2;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee6 - DNenr(2,1)*crhs_ee7 - Nenr[2]*crhs_ee2;


    //substitute_enrichment_rhs_eV_2D

    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
    //noalias(rRHS_eV) += rData.Weight * rhs_eV;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointEnrichmentContributionsCut(
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

    auto &V = rData.V;
    auto &H = rData.H;
    auto &Kee = rData.Kee;
    auto &rhs_ee = rData.rhs_ee;

    array_1d<double, NumNodes> penr = ZeroVector(NumNodes); 
    //const auto &penr = rData.Pressure_Enriched;

    const double cV0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
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


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 =             cH3*(K_darcy*N[0] + rho*(DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0));
const double cH5 =             cH3*(K_darcy*N[1] + rho*(DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0));
const double cH6 =             cH3*(K_darcy*N[2] + rho*(DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0));
const double cH7 =             cH3*(K_darcy*N[3] + rho*(DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0));
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


    const double cKee0 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
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
const double crhs_ee3 =             crhs_ee0 + crhs_ee1 + crhs_ee2;
const double crhs_ee4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee6 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee7 =             1.0/(K_darcy + rho*stab_c2*sqrt(pow(crhs_ee4, 2) + pow(crhs_ee5, 2) + pow(crhs_ee6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee8 =             crhs_ee7*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] + K_darcy*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*crhs_ee4 + crhs_ee5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
const double crhs_ee9 =             crhs_ee7*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] + K_darcy*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee1*crhs_ee5 + crhs_ee4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
const double crhs_ee10 =             crhs_ee7*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] + K_darcy*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee2*crhs_ee6 + crhs_ee4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2))));
            rhs_ee[0]=-DNenr(0,0)*crhs_ee8 - DNenr(0,1)*crhs_ee9 - DNenr(0,2)*crhs_ee10 - Nenr[0]*crhs_ee3;
            rhs_ee[1]=-DNenr(1,0)*crhs_ee8 - DNenr(1,1)*crhs_ee9 - DNenr(1,2)*crhs_ee10 - Nenr[1]*crhs_ee3;
            rhs_ee[2]=-DNenr(2,0)*crhs_ee8 - DNenr(2,1)*crhs_ee9 - DNenr(2,2)*crhs_ee10 - Nenr[2]*crhs_ee3;
            rhs_ee[3]=-DNenr(3,0)*crhs_ee8 - DNenr(3,1)*crhs_ee9 - DNenr(3,2)*crhs_ee10 - Nenr[3]*crhs_ee3;


    //substitute_enrichment_rhs_eV_3D

    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
    //noalias(rRHS_eV) += rData.Weight * rhs_eV;
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::ContactSurfaceDissipation(
    const TElementData& rData,
    const double betaIn,
    const double betaOut,
    const double betaContact,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    GeometryType::Pointer p_geom = this->pGetGeometry();
    GeometryType::GeometriesArrayType faces = p_geom->Faces();
    //GeometryType::GeometriesArrayType faces = p_geom->GenerateFaces();

    bool not_found_contact = true;
    const unsigned int n_faces = 4;
    const unsigned int n_nodes = 4;
    const unsigned int n_dim = n_nodes - 1;
    unsigned int i_face = 0;

    while (i_face < n_faces && not_found_contact) { // allows ONLY ONE face to be solid (slip) surface
        GeometryType& r_face = faces[i_face];
        //KRATOS_INFO("ContactSurfaceDissipation, i_face") << i_face << std::endl;
        unsigned int contact_node = 0;
        const unsigned int n_face_nodes = 3;
        for (unsigned int i=0; i < n_face_nodes/* r_face.size() */; ++i){
            //KRATOS_INFO("ContactSurfaceDissipation, face_nodes") << i << std::endl;
            if ( r_face[i].GetValue(IS_STRUCTURE) == 1.0 ){
                contact_node++;
                //KRATOS_INFO("ContactSurfaceDissipation, contact node ID") << r_face[i].Id() << std::endl;
            }
        }

        if (contact_node == n_face_nodes){
            //KRATOS_INFO("ContactSurfaceDissipation, contact i_face") << i_face << std::endl;
            not_found_contact = false;

            MatrixType lhs_dissipation = ZeroMatrix(n_nodes*(n_dim+1),n_nodes*(n_dim+1));
            Kratos::array_1d<double,n_nodes*(n_dim+1)> tempU; // Only the last known velocity
            Kratos::array_1d<double,n_nodes*(n_dim+1)> tempU0; // Only the last step velocity
            for (unsigned int i = 0; i < n_nodes; i++){
                const VectorType Vel = (*p_geom)[i].FastGetSolutionStepValue(VELOCITY);
                const VectorType Vel0 = (*p_geom)[i].FastGetSolutionStepValue(VELOCITY, 1);
                for (unsigned int dimi = 0; dimi < n_dim; dimi++){
                    tempU[i*(n_dim+1) + dimi] = //rData.Velocity(i,dimi);
                        Vel[dimi];
                    
                    tempU0[i*(n_dim+1) + dimi] = Vel0[dimi];
                }
            }

            double beta = 0.0;
            unsigned int n_minus = 0;
            unsigned int n_plus = 0;
            
            for (unsigned int i=0; i < n_face_nodes; ++i){
                if ( r_face[i].FastGetSolutionStepValue(DISTANCE) > 0.0 ){
                    n_plus++;
                }
                else{
                    n_minus++;
                }
            }

            if (n_plus == 0){
                beta = betaIn;
            } 
            else if (n_minus == 0){
                beta = betaOut;
            }  
            else{
                beta = betaContact;
            }

            //MatrixType FaceShapeFunctions;
            //GeometryType::ShapeFunctionsGradientsType FaceShapeFunctionsGradients;
            //Kratos::Vector FaceWeights;

            auto IntegrationMethod = GeometryData::GI_GAUSS_2;

            const unsigned int n_int_pts = (faces[i_face]).IntegrationPointsNumber(IntegrationMethod);
            //KRATOS_INFO("ContactSurfaceDissipation, n_int_pts") << n_int_pts << std::endl;

            //FaceShapeFunctions.resize(n_int_pts, n_nodes, false);
            //FaceShapeFunctionsGradients.resize(n_int_pts, false);
            //FaceWeights.resize(n_int_pts, false);

            std::vector < GeometryType::CoordinatesArrayType > face_gauss_pts_gl_coords, face_gauss_pts_loc_coords;
            face_gauss_pts_gl_coords.clear();
            face_gauss_pts_loc_coords.clear();
            face_gauss_pts_gl_coords.reserve(n_int_pts);
            face_gauss_pts_loc_coords.reserve(n_int_pts);

            std::vector< IntegrationPoint<3> > face_gauss_pts;
            face_gauss_pts = (faces[i_face]).IntegrationPoints(IntegrationMethod);

            Kratos::Vector face_jacobians;
            (faces[i_face]).DeterminantOfJacobian(face_jacobians, IntegrationMethod);

            // Get the original geometry shape function and gradients values over the intersection
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                // Store the Gauss points weights
               const double face_weight /* FaceWeights (i_gauss) */ = face_jacobians(i_gauss) * face_gauss_pts[i_gauss].Weight();

                // Compute the global coordinates of the face Gauss pt.
                GeometryType::CoordinatesArrayType global_coords = ZeroVector(3);
                global_coords = (faces[i_face]).GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                // Compute the parent geometry local coordinates of the Gauss pt.
                GeometryType::CoordinatesArrayType loc_coords = ZeroVector(3);
                loc_coords = (*p_geom).PointLocalCoordinates(loc_coords, global_coords);

                // Compute shape function values
                // Obtain the parent subgeometry shape function values
                double det_jac;
                Kratos::Vector face_shape_func;
                face_shape_func = (*p_geom).ShapeFunctionsValues(face_shape_func, loc_coords);

                /* for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    FaceShapeFunctions(i_gauss, i_node) = face_shape_func(i_node);
                } */
                //KRATOS_INFO("ContactSurfaceDissipation, gp") << i_gauss << ", Shape Functions: " << face_shape_func << std::endl;

                for (unsigned int i = 0; i < n_nodes; i++){
                    //KRATOS_INFO("Face Shape Function") << face_shape_func(i) 
                    //    << ", IS_STRUCTURE: " << (*p_geom)[i].GetValue(IS_STRUCTURE) << std::endl;
                    for (unsigned int j = 0; j < n_nodes; j++){
                        for (unsigned int dimi = 0; dimi < n_dim; dimi++){
                            //for (unsigned int dimj = 0; dimj < n_dim; dimj++){ //Be carefull not to confuse
                            // Dimension_j is useless since the force acts in the dimension_i direction
                                lhs_dissipation( i*(n_dim+1) + dimi, j*(n_dim+1) + dimi) += 
                                    beta * face_weight * face_shape_func(j) * face_shape_func(i);
                            //}
                        }
                    }
                }
            }

            noalias(rLHS) += lhs_dissipation;
            noalias(rRHS) -= prod(lhs_dissipation,tempU);

            /* const VectorType temp_force = -prod(lhs_dissipation,tempU0);
            for (unsigned int iii = 0; iii < n_nodes; iii++){
                VectorType temp_forcei = ZeroVector(3);
                for (unsigned int dimi = 0; dimi < n_dim; dimi++){
                    temp_forcei[dimi] = temp_force[iii*(n_dim+1) + dimi];
                }
                (*p_geom)[iii].FastGetSolutionStepValue(TANGENT_VECTOR) = temp_forcei;
            } */
        }
        
        i_face++;
    }
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    GeometryType::Pointer p_geom = this->pGetGeometry();
    //const auto& r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = p_geom->IntegrationPoints(p_geom->GetDefaultIntegrationMethod());

    //const auto& rGeom = this->GetGeometry();
    //const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    //const unsigned int num_gauss = IntegrationPoints.size();

    const std::size_t number_of_integration_points = integration_points.size();

    if ( rOutput.size() != number_of_integration_points ) {
        rOutput.resize( number_of_integration_points );
    }

    KRATOS_INFO("CalculateOnIntegrationPoints") << "CALLED!" << std::endl;

    if (/* this->Has( rVariable) && */ rVariable == CONTACT_ANGLE) {
        const double value = this->GetValue( rVariable);
        for (std::size_t point_number = 0; point_number < number_of_integration_points; ++point_number) {
            KRATOS_INFO("CalculateOnIntegrationPoints") << "INSIDE IF!" << std::endl;
            rOutput[point_number] = value;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>;
template class TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>;

} // namespace Kratos