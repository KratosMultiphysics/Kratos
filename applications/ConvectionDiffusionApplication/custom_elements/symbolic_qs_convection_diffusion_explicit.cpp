//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#include "symbolic_qs_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicQSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicQSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicQSConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicQSConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicQSConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if(rResult.size() != local_size)
        rResult.resize(local_size);
    for (unsigned int i=0; i<local_size; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);
    for (unsigned int i = 0; i < local_size; i++)
    {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    BoundedVector<double, TNumNodes> rhs;
    this->CalculateRightHandSideInternal(rhs,rCurrentProcessInfo);
    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    const auto& reaction_variable = r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable();
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(reaction_variable) += rhs[i_node];
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<2,3>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int local_size = 3;
    // Resize, intialize and fill the mass matrix for linear triangular elements
    if (rMassMatrix.size1() != local_size)
        rMassMatrix.resize(local_size, local_size, false);
    noalias(rMassMatrix) = ZeroMatrix(local_size, local_size);
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    rMassMatrix(0,0) = one_six; rMassMatrix(0,1) = one_twelve; rMassMatrix(0,2) = one_twelve;
    rMassMatrix(1,0) = one_twelve; rMassMatrix(1,1) = one_six; rMassMatrix(1,2) = one_twelve;
    rMassMatrix(2,0) = one_twelve; rMassMatrix(2,1) = one_twelve; rMassMatrix(2,2) = one_six;
    // Assumption all the Gauss points have the same weight, so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int local_size = 4;
    // Resize, intialize and fill the mass matrix for linear tetrahedral elements
    if (rMassMatrix.size1() != local_size)
        rMassMatrix.resize(local_size, local_size, false);
    noalias(rMassMatrix) = ZeroMatrix(local_size, local_size);
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    rMassMatrix(0,0) = one_ten; rMassMatrix(0,1) = one_twenty; rMassMatrix(0,2) = one_twenty; rMassMatrix(0,3) = one_twenty;
    rMassMatrix(1,0) = one_twenty; rMassMatrix(1,1) = one_ten; rMassMatrix(1,2) = one_twenty; rMassMatrix(1,3) = one_twenty;
    rMassMatrix(2,0) = one_twenty; rMassMatrix(2,1) = one_twenty; rMassMatrix(2,2) = one_ten; rMassMatrix(2,3) = one_twenty;
    rMassMatrix(3,0) = one_twenty; rMassMatrix(3,1) = one_twenty; rMassMatrix(3,2) = one_twenty; rMassMatrix(3,3) = one_ten;
    // Assumption all the Gauss points have the same weight, so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == SCALAR_PROJECTION) {
        auto& r_geometry = GetGeometry();
        const unsigned int local_size = r_geometry.size();
        BoundedVector<double, TNumNodes> rhs_oss;
        this->CalculateOrthogonalSubgridScaleRHSInternal(rhs_oss,rCurrentProcessInfo);
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
        }
    }
    else {
        BaseType::Calculate(rVariable,Output,rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::InitializeEulerianElement(
    ElementData& rData,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];

    // Getting data for the given geometry and integration method
    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    array_1d<double,TNumNodes> N_aux;
    GeometryUtils::CalculateGeometryData(r_geometry,rData.DN_DX,N_aux,rData.volume);
    rData.N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize some scalar data
    rData.lumping_factor = 1.00 / double(TNumNodes);
    rData.diffusivity = 0.0;
    rData.dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

    for(unsigned int node_element = 0; node_element<local_size; node_element++) {
        // Observations
        // * SGS time derivative term approximated as (phi-phi_old)/(delta_time_coefficient*delta_time)
        //   observe that for RK step = 1 ASGS time derivative term = 0 because phi = phi_old (null acceleration wrt step n)
        // * convective velocity and forcing term:
        //   RK step 1: evaluated at previous time step
        //   RK steps 2 and 3: linear interpolation between current and oldprevious time step
        //   RK step 4: evaluated at current time step
        // * convective_velocity = velocity - velocity_mesh
        //   velocity_mesh = 0 in eulerian framework
        if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==1) {
            rData.delta_time_coefficient = std::numeric_limits<double>::max();
            rData.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable(),1);
            rData.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[0];
            rData.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[1];
            rData.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[2];
        }
        else if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==2 || r_process_info.GetValue(RUNGE_KUTTA_STEP)==3) {
            rData.delta_time_coefficient = 0.5;
            rData.forcing[node_element] = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable()) + r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable(),1));
            rData.convective_velocity(node_element,0) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[0] + r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[0]);
            rData.convective_velocity(node_element,1) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[1] + r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[1]);
            rData.convective_velocity(node_element,2) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[2] + r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[2]);
        }
        else {
            rData.delta_time_coefficient = 1.0;
            rData.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable());
            rData.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[0];
            rData.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[1];
            rData.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[2];
        }
        rData.oss_projection[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetProjectionVariable());
        rData.diffusivity += r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetDiffusionVariable());
        rData.delta_time = r_process_info[DELTA_TIME];
        rData.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetUnknownVariable());
        rData.unknown_old[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetUnknownVariable(),1);
    }
    // divide by number of nodes scalar data
    rData.diffusivity *= rData.lumping_factor;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& delta_time_coefficient = rData.delta_time_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             0.25*f[1];
const double crhs1 =             0.166666666666667*v(0,0);
const double crhs2 =             0.166666666666667*v(2,0);
const double crhs3 =             crhs1 + crhs2 + 0.666666666666667*v(1,0);
const double crhs4 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs5 =             crhs3*crhs4;
const double crhs6 =             -0.166666666666667*crhs5;
const double crhs7 =             0.166666666666667*v(0,1);
const double crhs8 =             0.166666666666667*v(2,1);
const double crhs9 =             crhs7 + crhs8 + 0.666666666666667*v(1,1);
const double crhs10 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs11 =             crhs10*crhs9;
const double crhs12 =             -0.166666666666667*crhs11;
const double crhs13 =             0.166666666666667*phi[0];
const double crhs14 =             0.166666666666667*phi[2];
const double crhs15 =             crhs13 + crhs14 + 0.666666666666667*phi[1];
const double crhs16 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs17 =             crhs15*crhs16;
const double crhs18 =             -0.166666666666667*crhs17;
const double crhs19 =             0.25*f[2];
const double crhs20 =             3*alpha*crhs4;
const double crhs21 =             3*alpha*crhs10;
const double crhs22 =             0.166666666666667*v(1,0);
const double crhs23 =             crhs1 + crhs22 + 0.666666666666667*v(2,0);
const double crhs24 =             crhs23*crhs4;
const double crhs25 =             -0.166666666666667*crhs24;
const double crhs26 =             crhs2 + crhs22 + 0.666666666666667*v(0,0);
const double crhs27 =             crhs26*crhs4;
const double crhs28 =             0.166666666666667*v(1,1);
const double crhs29 =             crhs28 + crhs7 + 0.666666666666667*v(2,1);
const double crhs30 =             crhs10*crhs29;
const double crhs31 =             -0.166666666666667*crhs30;
const double crhs32 =             crhs28 + crhs8 + 0.666666666666667*v(0,1);
const double crhs33 =             crhs10*crhs32;
const double crhs34 =             0.166666666666667*phi[1];
const double crhs35 =             crhs13 + crhs34 + 0.666666666666667*phi[2];
const double crhs36 =             crhs16*crhs35;
const double crhs37 =             -0.166666666666667*crhs36;
const double crhs38 =             crhs14 + crhs34 + 0.666666666666667*phi[0];
const double crhs39 =             crhs16*crhs38;
const double crhs40 =             0.166666666666667*f[1];
const double crhs41 =             0.166666666666667*f[2];
const double crhs42 =             crhs40 + crhs41 + 0.666666666666667*f[0];
const double crhs43 =             tau[0]*(DN_DX_0_0*crhs26 + DN_DX_0_1*crhs32);
const double crhs44 =             0.166666666666667*prj[1];
const double crhs45 =             0.166666666666667*prj[2];
const double crhs46 =             crhs44 + crhs45 + 0.666666666666667*prj[0];
const double crhs47 =             0.166666666666667*f[0];
const double crhs48 =             crhs41 + crhs47 + 0.666666666666667*f[1];
const double crhs49 =             tau[1]*(DN_DX_0_0*crhs3 + DN_DX_0_1*crhs9);
const double crhs50 =             0.166666666666667*prj[0];
const double crhs51 =             crhs45 + crhs50 + 0.666666666666667*prj[1];
const double crhs52 =             crhs40 + crhs47 + 0.666666666666667*f[2];
const double crhs53 =             tau[2]*(DN_DX_0_0*crhs23 + DN_DX_0_1*crhs29);
const double crhs54 =             crhs44 + crhs50 + 0.666666666666667*prj[2];
const double crhs55 =             1.0/delta_time;
const double crhs56 =             1.0/delta_time_coefficient;
const double crhs57 =             -0.166666666666667*phi_old[1];
const double crhs58 =             -0.166666666666667*phi_old[2];
const double crhs59 =             crhs55*crhs56*(crhs38 + crhs57 + crhs58 - 0.666666666666667*phi_old[0]);
const double crhs60 =             -0.166666666666667*phi_old[0];
const double crhs61 =             crhs55*crhs56*(crhs15 + crhs58 + crhs60 - 0.666666666666667*phi_old[1]);
const double crhs62 =             crhs55*crhs56*(crhs35 + crhs57 + crhs60 - 0.666666666666667*phi_old[2]);
const double crhs63 =             crhs27 + crhs33;
const double crhs64 =             crhs11 + crhs5;
const double crhs65 =             crhs24 + crhs30;
const double crhs66 =             -0.166666666666667*crhs27 - 0.166666666666667*crhs33 - 0.166666666666667*crhs39 + 0.25*f[0];
const double crhs67 =             tau[0]*(DN_DX_1_0*crhs26 + DN_DX_1_1*crhs32);
const double crhs68 =             tau[1]*(DN_DX_1_0*crhs3 + DN_DX_1_1*crhs9);
const double crhs69 =             tau[2]*(DN_DX_1_0*crhs23 + DN_DX_1_1*crhs29);
const double crhs70 =             tau[0]*(DN_DX_2_0*crhs26 + DN_DX_2_1*crhs32);
const double crhs71 =             tau[1]*(DN_DX_2_0*crhs3 + DN_DX_2_1*crhs9);
const double crhs72 =             tau[2]*(DN_DX_2_0*crhs23 + DN_DX_2_1*crhs29);
            rhs[0]=-DN_DX_0_0*crhs20 - DN_DX_0_1*crhs21 + crhs0 + crhs12 - crhs17*crhs49 + crhs18 + crhs19 + crhs25 - 0.666666666666667*crhs27 + crhs31 - 0.666666666666667*crhs33 - crhs36*crhs53 + crhs37 - crhs39*crhs43 - 0.666666666666667*crhs39 + crhs42*crhs43 + crhs43*crhs46 - crhs43*crhs59 - crhs43*crhs63 + crhs48*crhs49 + crhs49*crhs51 - crhs49*crhs61 - crhs49*crhs64 + crhs52*crhs53 + crhs53*crhs54 - crhs53*crhs62 - crhs53*crhs65 + crhs6 + 0.5*f[0];
            rhs[1]=-DN_DX_1_0*crhs20 - DN_DX_1_1*crhs21 - 0.666666666666667*crhs11 - crhs17*crhs68 - 0.666666666666667*crhs17 + crhs19 + crhs25 + crhs31 - crhs36*crhs69 + crhs37 - crhs39*crhs67 + crhs42*crhs67 + crhs46*crhs67 + crhs48*crhs68 - 0.666666666666667*crhs5 + crhs51*crhs68 + crhs52*crhs69 + crhs54*crhs69 - crhs59*crhs67 - crhs61*crhs68 - crhs62*crhs69 - crhs63*crhs67 - crhs64*crhs68 - crhs65*crhs69 + crhs66 + 0.5*f[1];
            rhs[2]=-DN_DX_2_0*crhs20 - DN_DX_2_1*crhs21 + crhs0 + crhs12 - crhs17*crhs71 + crhs18 - 0.666666666666667*crhs24 - 0.666666666666667*crhs30 - crhs36*crhs72 - 0.666666666666667*crhs36 - crhs39*crhs70 + crhs42*crhs70 + crhs46*crhs70 + crhs48*crhs71 + crhs51*crhs71 + crhs52*crhs72 + crhs54*crhs72 - crhs59*crhs70 + crhs6 - crhs61*crhs71 - crhs62*crhs72 - crhs63*crhs70 - crhs64*crhs71 - crhs65*crhs72 + crhs66 + 0.5*f[2];


    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3,4>::CalculateRightHandSideInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& delta_time_coefficient = rData.delta_time_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             0.19999999899376*f[1];
const double crhs1 =             0.1381966*v(0,0);
const double crhs2 =             0.1381966*v(2,0);
const double crhs3 =             0.1381966*v(3,0);
const double crhs4 =             crhs2 + crhs3;
const double crhs5 =             crhs1 + crhs4 + 0.5854102*v(1,0);
const double crhs6 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs7 =             crhs5*crhs6;
const double crhs8 =             -0.1381966*crhs7;
const double crhs9 =             0.1381966*v(0,1);
const double crhs10 =             0.1381966*v(2,1);
const double crhs11 =             0.1381966*v(3,1);
const double crhs12 =             crhs10 + crhs11;
const double crhs13 =             crhs12 + crhs9 + 0.5854102*v(1,1);
const double crhs14 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs15 =             crhs13*crhs14;
const double crhs16 =             -0.1381966*crhs15;
const double crhs17 =             0.1381966*v(0,2);
const double crhs18 =             0.1381966*v(2,2);
const double crhs19 =             0.1381966*v(3,2);
const double crhs20 =             crhs18 + crhs19;
const double crhs21 =             crhs17 + crhs20 + 0.5854102*v(1,2);
const double crhs22 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs23 =             crhs21*crhs22;
const double crhs24 =             -0.1381966*crhs23;
const double crhs25 =             0.1381966*phi[0];
const double crhs26 =             0.5854102*phi[1];
const double crhs27 =             0.1381966*phi[2];
const double crhs28 =             0.1381966*phi[3];
const double crhs29 =             crhs27 + crhs28;
const double crhs30 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs31 =             crhs30*(crhs25 + crhs26 + crhs29);
const double crhs32 =             -0.1381966*crhs31;
const double crhs33 =             0.19999999899376*f[2];
const double crhs34 =             0.19999999899376*f[3];
const double crhs35 =             4*alpha*crhs6;
const double crhs36 =             4*alpha*crhs14;
const double crhs37 =             4*alpha*crhs22;
const double crhs38 =             0.1381966*v(1,0);
const double crhs39 =             crhs1 + crhs38;
const double crhs40 =             crhs2 + crhs39 + 0.5854102*v(3,0);
const double crhs41 =             crhs40*crhs6;
const double crhs42 =             -0.1381966*crhs41;
const double crhs43 =             crhs3 + crhs39 + 0.5854102*v(2,0);
const double crhs44 =             crhs43*crhs6;
const double crhs45 =             -0.1381966*crhs44;
const double crhs46 =             crhs38 + crhs4 + 0.5854102*v(0,0);
const double crhs47 =             crhs46*crhs6;
const double crhs48 =             0.1381966*v(1,1);
const double crhs49 =             crhs48 + crhs9;
const double crhs50 =             crhs10 + crhs49 + 0.5854102*v(3,1);
const double crhs51 =             crhs14*crhs50;
const double crhs52 =             -0.1381966*crhs51;
const double crhs53 =             crhs11 + crhs49 + 0.5854102*v(2,1);
const double crhs54 =             crhs14*crhs53;
const double crhs55 =             -0.1381966*crhs54;
const double crhs56 =             crhs12 + crhs48 + 0.5854102*v(0,1);
const double crhs57 =             crhs14*crhs56;
const double crhs58 =             0.1381966*v(1,2);
const double crhs59 =             crhs17 + crhs58;
const double crhs60 =             crhs18 + crhs59 + 0.5854102*v(3,2);
const double crhs61 =             crhs22*crhs60;
const double crhs62 =             -0.1381966*crhs61;
const double crhs63 =             crhs19 + crhs59 + 0.5854102*v(2,2);
const double crhs64 =             crhs22*crhs63;
const double crhs65 =             -0.1381966*crhs64;
const double crhs66 =             crhs20 + crhs58 + 0.5854102*v(0,2);
const double crhs67 =             crhs22*crhs66;
const double crhs68 =             0.1381966*phi[1];
const double crhs69 =             crhs25 + crhs68;
const double crhs70 =             0.5854102*phi[3];
const double crhs71 =             crhs30*(crhs27 + crhs69 + crhs70);
const double crhs72 =             -0.1381966*crhs71;
const double crhs73 =             0.5854102*phi[2];
const double crhs74 =             crhs30*(crhs28 + crhs69 + crhs73);
const double crhs75 =             -0.1381966*crhs74;
const double crhs76 =             0.5854102*phi[0];
const double crhs77 =             crhs30*(crhs29 + crhs68 + crhs76);
const double crhs78 =             0.1381966*f[1];
const double crhs79 =             0.1381966*f[2];
const double crhs80 =             0.1381966*f[3];
const double crhs81 =             crhs79 + crhs80;
const double crhs82 =             crhs78 + crhs81 + 0.5854102*f[0];
const double crhs83 =             tau[0]*(DN_DX_0_0*crhs46 + DN_DX_0_1*crhs56 + DN_DX_0_2*crhs66);
const double crhs84 =             0.1381966*prj[1];
const double crhs85 =             0.1381966*prj[2];
const double crhs86 =             0.1381966*prj[3];
const double crhs87 =             crhs85 + crhs86;
const double crhs88 =             crhs84 + crhs87 + 0.5854102*prj[0];
const double crhs89 =             0.1381966*f[0];
const double crhs90 =             crhs81 + crhs89 + 0.5854102*f[1];
const double crhs91 =             tau[1]*(DN_DX_0_0*crhs5 + DN_DX_0_1*crhs13 + DN_DX_0_2*crhs21);
const double crhs92 =             0.1381966*prj[0];
const double crhs93 =             crhs87 + crhs92 + 0.5854102*prj[1];
const double crhs94 =             crhs78 + crhs89;
const double crhs95 =             crhs80 + crhs94 + 0.5854102*f[2];
const double crhs96 =             tau[2]*(DN_DX_0_0*crhs43 + DN_DX_0_1*crhs53 + DN_DX_0_2*crhs63);
const double crhs97 =             crhs84 + crhs92;
const double crhs98 =             crhs86 + crhs97 + 0.5854102*prj[2];
const double crhs99 =             crhs79 + crhs94 + 0.5854102*f[3];
const double crhs100 =             tau[3]*(DN_DX_0_0*crhs40 + DN_DX_0_1*crhs50 + DN_DX_0_2*crhs60);
const double crhs101 =             crhs85 + crhs97 + 0.5854102*prj[3];
const double crhs102 =             1.0/delta_time;
const double crhs103 =             1.0/delta_time_coefficient;
const double crhs104 =             -0.1381966*phi_old[2];
const double crhs105 =             -0.1381966*phi_old[3];
const double crhs106 =             crhs104 + crhs105;
const double crhs107 =             -0.1381966*phi_old[1];
const double crhs108 =             crhs102*crhs103*(crhs106 + crhs107 + crhs27 + crhs28 + crhs68 + crhs76 - 0.5854102*phi_old[0]);
const double crhs109 =             -0.1381966*phi_old[0];
const double crhs110 =             crhs102*crhs103*(crhs106 + crhs109 + crhs25 + crhs26 + crhs27 + crhs28 - 0.5854102*phi_old[1]);
const double crhs111 =             crhs107 + crhs109;
const double crhs112 =             crhs102*crhs103*(crhs105 + crhs111 + crhs25 + crhs28 + crhs68 + crhs73 - 0.5854102*phi_old[2]);
const double crhs113 =             crhs102*crhs103*(crhs104 + crhs111 + crhs25 + crhs27 + crhs68 + crhs70 - 0.5854102*phi_old[3]);
const double crhs114 =             crhs47 + crhs57 + crhs67;
const double crhs115 =             crhs15 + crhs23 + crhs7;
const double crhs116 =             crhs44 + crhs54 + crhs64;
const double crhs117 =             crhs41 + crhs51 + crhs61;
const double crhs118 =             0.19999999899376*f[0];
const double crhs119 =             -0.1381966*crhs47;
const double crhs120 =             -0.1381966*crhs57;
const double crhs121 =             -0.1381966*crhs67;
const double crhs122 =             -0.1381966*crhs77;
const double crhs123 =             tau[0]*(DN_DX_1_0*crhs46 + DN_DX_1_1*crhs56 + DN_DX_1_2*crhs66);
const double crhs124 =             tau[1]*(DN_DX_1_0*crhs5 + DN_DX_1_1*crhs13 + DN_DX_1_2*crhs21);
const double crhs125 =             tau[2]*(DN_DX_1_0*crhs43 + DN_DX_1_1*crhs53 + DN_DX_1_2*crhs63);
const double crhs126 =             tau[3]*(DN_DX_1_0*crhs40 + DN_DX_1_1*crhs50 + DN_DX_1_2*crhs60);
const double crhs127 =             crhs0 + crhs118 + crhs119 + crhs120 + crhs121 + crhs122 + crhs16 + crhs24 + crhs32 + crhs8;
const double crhs128 =             tau[0]*(DN_DX_2_0*crhs46 + DN_DX_2_1*crhs56 + DN_DX_2_2*crhs66);
const double crhs129 =             tau[1]*(DN_DX_2_0*crhs5 + DN_DX_2_1*crhs13 + DN_DX_2_2*crhs21);
const double crhs130 =             tau[2]*(DN_DX_2_0*crhs43 + DN_DX_2_1*crhs53 + DN_DX_2_2*crhs63);
const double crhs131 =             tau[3]*(DN_DX_2_0*crhs40 + DN_DX_2_1*crhs50 + DN_DX_2_2*crhs60);
const double crhs132 =             tau[0]*(DN_DX_3_0*crhs46 + DN_DX_3_1*crhs56 + DN_DX_3_2*crhs66);
const double crhs133 =             tau[1]*(DN_DX_3_0*crhs5 + DN_DX_3_1*crhs13 + DN_DX_3_2*crhs21);
const double crhs134 =             tau[2]*(DN_DX_3_0*crhs43 + DN_DX_3_1*crhs53 + DN_DX_3_2*crhs63);
const double crhs135 =             tau[3]*(DN_DX_3_0*crhs40 + DN_DX_3_1*crhs50 + DN_DX_3_2*crhs60);
            rhs[0]=-DN_DX_0_0*crhs35 - DN_DX_0_1*crhs36 - DN_DX_0_2*crhs37 + crhs0 + crhs100*crhs101 - crhs100*crhs113 - crhs100*crhs117 - crhs100*crhs71 + crhs100*crhs99 - crhs108*crhs83 - crhs110*crhs91 - crhs112*crhs96 - crhs114*crhs83 - crhs115*crhs91 - crhs116*crhs96 + crhs16 + crhs24 - crhs31*crhs91 + crhs32 + crhs33 + crhs34 + crhs42 + crhs45 - 0.5854102*crhs47 + crhs52 + crhs55 - 0.5854102*crhs57 + crhs62 + crhs65 - 0.5854102*crhs67 + crhs72 - crhs74*crhs96 + crhs75 - crhs77*crhs83 - 0.5854102*crhs77 + crhs8 + crhs82*crhs83 + crhs83*crhs88 + crhs90*crhs91 + crhs91*crhs93 + crhs95*crhs96 + crhs96*crhs98 + 0.40000000301872*f[0];
            rhs[1]=-DN_DX_1_0*crhs35 - DN_DX_1_1*crhs36 - DN_DX_1_2*crhs37 + crhs101*crhs126 - crhs108*crhs123 - crhs110*crhs124 - crhs112*crhs125 - crhs113*crhs126 - crhs114*crhs123 - crhs115*crhs124 - crhs116*crhs125 - crhs117*crhs126 + crhs118 + crhs119 + crhs120 + crhs121 + crhs122 - crhs123*crhs77 + crhs123*crhs82 + crhs123*crhs88 - crhs124*crhs31 + crhs124*crhs90 + crhs124*crhs93 - crhs125*crhs74 + crhs125*crhs95 + crhs125*crhs98 - crhs126*crhs71 + crhs126*crhs99 - 0.5854102*crhs15 - 0.5854102*crhs23 - 0.5854102*crhs31 + crhs33 + crhs34 + crhs42 + crhs45 + crhs52 + crhs55 + crhs62 + crhs65 - 0.5854102*crhs7 + crhs72 + crhs75 + 0.40000000301872*f[1];
            rhs[2]=-DN_DX_2_0*crhs35 - DN_DX_2_1*crhs36 - DN_DX_2_2*crhs37 + crhs101*crhs131 - crhs108*crhs128 - crhs110*crhs129 - crhs112*crhs130 - crhs113*crhs131 - crhs114*crhs128 - crhs115*crhs129 - crhs116*crhs130 - crhs117*crhs131 + crhs127 - crhs128*crhs77 + crhs128*crhs82 + crhs128*crhs88 - crhs129*crhs31 + crhs129*crhs90 + crhs129*crhs93 - crhs130*crhs74 + crhs130*crhs95 + crhs130*crhs98 - crhs131*crhs71 + crhs131*crhs99 + crhs42 - 0.5854102*crhs44 + crhs52 - 0.5854102*crhs54 + crhs62 - 0.5854102*crhs64 + crhs72 - 0.5854102*crhs74 + 0.40000000301872*f[2] + 0.19999999899376*f[3];
            rhs[3]=-DN_DX_3_0*crhs35 - DN_DX_3_1*crhs36 - DN_DX_3_2*crhs37 + crhs101*crhs135 - crhs108*crhs132 - crhs110*crhs133 - crhs112*crhs134 - crhs113*crhs135 - crhs114*crhs132 - crhs115*crhs133 - crhs116*crhs134 - crhs117*crhs135 + crhs127 - crhs132*crhs77 + crhs132*crhs82 + crhs132*crhs88 - crhs133*crhs31 + crhs133*crhs90 + crhs133*crhs93 - crhs134*crhs74 + crhs134*crhs95 + crhs134*crhs98 - crhs135*crhs71 + crhs135*crhs99 - 0.5854102*crhs41 + crhs45 - 0.5854102*crhs51 + crhs55 - 0.5854102*crhs61 + crhs65 - 0.5854102*crhs71 + crhs75 + 0.19999999899376*f[2] + 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<2,3>::CalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& delta_time_coefficient = rData.delta_time_coefficient;
    const auto& v = rData.convective_velocity;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             -0.25*f[1];
const double crhs1 =             0.166666666666667*v(0,0);
const double crhs2 =             0.166666666666667*v(2,0);
const double crhs3 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs4 =             crhs3*(crhs1 + crhs2 + 0.666666666666667*v(1,0));
const double crhs5 =             0.166666666666667*crhs4;
const double crhs6 =             0.166666666666667*v(0,1);
const double crhs7 =             0.166666666666667*v(2,1);
const double crhs8 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs9 =             crhs8*(crhs6 + crhs7 + 0.666666666666667*v(1,1));
const double crhs10 =             0.166666666666667*crhs9;
const double crhs11 =             1/(delta_time*delta_time_coefficient);
const double crhs12 =             0.111111111111111*phi[1];
const double crhs13 =             -0.111111111111111*phi_old[1];
const double crhs14 =             crhs12 + crhs13;
const double crhs15 =             0.0277777777777778*phi[0];
const double crhs16 =             0.0277777777777778*phi[2];
const double crhs17 =             -0.0277777777777778*phi_old[0];
const double crhs18 =             -0.0277777777777778*phi_old[2];
const double crhs19 =             crhs11*(crhs14 + crhs15 + crhs16 + crhs17 + crhs18);
const double crhs20 =             0.166666666666667*phi[0];
const double crhs21 =             0.166666666666667*phi[2];
const double crhs22 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs23 =             crhs22*(crhs20 + crhs21 + 0.666666666666667*phi[1]);
const double crhs24 =             0.166666666666667*crhs23;
const double crhs25 =             -0.25*f[2];
const double crhs26 =             3*alpha*crhs3;
const double crhs27 =             3*alpha*crhs8;
const double crhs28 =             0.166666666666667*v(1,0);
const double crhs29 =             crhs3*(crhs1 + crhs28 + 0.666666666666667*v(2,0));
const double crhs30 =             0.166666666666667*crhs29;
const double crhs31 =             crhs3*(crhs2 + crhs28 + 0.666666666666667*v(0,0));
const double crhs32 =             0.166666666666667*v(1,1);
const double crhs33 =             crhs8*(crhs32 + crhs6 + 0.666666666666667*v(2,1));
const double crhs34 =             0.166666666666667*crhs33;
const double crhs35 =             crhs8*(crhs32 + crhs7 + 0.666666666666667*v(0,1));
const double crhs36 =             0.111111111111111*phi[2];
const double crhs37 =             -0.111111111111111*phi_old[2];
const double crhs38 =             0.0277777777777778*phi[1];
const double crhs39 =             -0.0277777777777778*phi_old[1];
const double crhs40 =             crhs11*(crhs15 + crhs17 + crhs36 + crhs37 + crhs38 + crhs39);
const double crhs41 =             0.166666666666667*phi[1];
const double crhs42 =             crhs22*(crhs20 + crhs41 + 0.666666666666667*phi[2]);
const double crhs43 =             0.166666666666667*crhs42;
const double crhs44 =             crhs22*(crhs21 + crhs41 + 0.666666666666667*phi[0]);
const double crhs45 =             0.111111111111111*phi[0] - 0.111111111111111*phi_old[0];
const double crhs46 =             crhs11*(crhs16 + crhs18 + crhs38 + crhs39 + crhs45) + 0.166666666666667*crhs31 + 0.166666666666667*crhs35 + 0.166666666666667*crhs44 - 0.25*f[0];
            rhs[0]=DN_DX_0_0*crhs26 + DN_DX_0_1*crhs27 + crhs0 + crhs10 + crhs11*(crhs14 + crhs36 + crhs37 + 0.444444444444444*phi[0] - 0.444444444444444*phi_old[0]) + crhs19 + crhs24 + crhs25 + crhs30 + 0.666666666666667*crhs31 + crhs34 + 0.666666666666667*crhs35 + crhs40 + crhs43 + 0.666666666666667*crhs44 + crhs5 - 0.5*f[0];
            rhs[1]=DN_DX_1_0*crhs26 + DN_DX_1_1*crhs27 + crhs11*(crhs36 + crhs37 + crhs45 + 0.444444444444444*phi[1] - 0.444444444444444*phi_old[1]) + 0.666666666666667*crhs23 + crhs25 + crhs30 + crhs34 + 0.666666666666667*crhs4 + crhs40 + crhs43 + crhs46 + 0.666666666666667*crhs9 - 0.5*f[1];
            rhs[2]=DN_DX_2_0*crhs26 + DN_DX_2_1*crhs27 + crhs0 + crhs10 + crhs11*(crhs12 + crhs13 + crhs45 + 0.444444444444444*phi[2] - 0.444444444444444*phi_old[2]) + crhs19 + crhs24 + 0.666666666666667*crhs29 + 0.666666666666667*crhs33 + 0.666666666666667*crhs42 + crhs46 + crhs5 - 0.5*f[2];


    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void SymbolicQSConvectionDiffusionExplicit<3,4>::CalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->CalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& delta_time_coefficient = rData.delta_time_coefficient;
    const auto& v = rData.convective_velocity;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

    const double crhs0 =             -0.19999999899376*f[1];
const double crhs1 =             0.1381966*v(0,0);
const double crhs2 =             0.1381966*v(2,0);
const double crhs3 =             0.1381966*v(3,0);
const double crhs4 =             crhs2 + crhs3;
const double crhs5 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs6 =             crhs5*(crhs1 + crhs4 + 0.5854102*v(1,0));
const double crhs7 =             0.1381966*crhs6;
const double crhs8 =             0.1381966*v(0,1);
const double crhs9 =             0.1381966*v(2,1);
const double crhs10 =             0.1381966*v(3,1);
const double crhs11 =             crhs10 + crhs9;
const double crhs12 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs13 =             crhs12*(crhs11 + crhs8 + 0.5854102*v(1,1));
const double crhs14 =             0.1381966*crhs13;
const double crhs15 =             0.1381966*v(0,2);
const double crhs16 =             0.1381966*v(2,2);
const double crhs17 =             0.1381966*v(3,2);
const double crhs18 =             crhs16 + crhs17;
const double crhs19 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs20 =             crhs19*(crhs15 + crhs18 + 0.5854102*v(1,2));
const double crhs21 =             0.1381966*crhs20;
const double crhs22 =             1/(delta_time*delta_time_coefficient);
const double crhs23 =             0.08090169924532*phi[1];
const double crhs24 =             -0.08090169924532*phi_old[1];
const double crhs25 =             0.01909830025156*phi[0];
const double crhs26 =             0.01909830025156*phi[2];
const double crhs27 =             0.01909830025156*phi[3];
const double crhs28 =             -0.01909830025156*phi_old[0];
const double crhs29 =             -0.01909830025156*phi_old[2];
const double crhs30 =             -0.01909830025156*phi_old[3];
const double crhs31 =             crhs22*(crhs23 + crhs24 + crhs25 + crhs26 + crhs27 + crhs28 + crhs29 + crhs30);
const double crhs32 =             0.1381966*phi[0];
const double crhs33 =             0.1381966*phi[2];
const double crhs34 =             0.1381966*phi[3];
const double crhs35 =             crhs33 + crhs34;
const double crhs36 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs37 =             crhs36*(crhs32 + crhs35 + 0.5854102*phi[1]);
const double crhs38 =             0.1381966*crhs37;
const double crhs39 =             -0.19999999899376*f[2];
const double crhs40 =             -0.19999999899376*f[3];
const double crhs41 =             4*alpha*crhs5;
const double crhs42 =             4*alpha*crhs12;
const double crhs43 =             4*alpha*crhs19;
const double crhs44 =             0.1381966*v(1,0);
const double crhs45 =             crhs1 + crhs44;
const double crhs46 =             crhs5*(crhs2 + crhs45 + 0.5854102*v(3,0));
const double crhs47 =             0.1381966*crhs46;
const double crhs48 =             crhs5*(crhs3 + crhs45 + 0.5854102*v(2,0));
const double crhs49 =             0.1381966*crhs48;
const double crhs50 =             crhs5*(crhs4 + crhs44 + 0.5854102*v(0,0));
const double crhs51 =             0.1381966*v(1,1);
const double crhs52 =             crhs51 + crhs8;
const double crhs53 =             crhs12*(crhs52 + crhs9 + 0.5854102*v(3,1));
const double crhs54 =             0.1381966*crhs53;
const double crhs55 =             crhs12*(crhs10 + crhs52 + 0.5854102*v(2,1));
const double crhs56 =             0.1381966*crhs55;
const double crhs57 =             crhs12*(crhs11 + crhs51 + 0.5854102*v(0,1));
const double crhs58 =             0.1381966*v(1,2);
const double crhs59 =             crhs15 + crhs58;
const double crhs60 =             crhs19*(crhs16 + crhs59 + 0.5854102*v(3,2));
const double crhs61 =             0.1381966*crhs60;
const double crhs62 =             crhs19*(crhs17 + crhs59 + 0.5854102*v(2,2));
const double crhs63 =             0.1381966*crhs62;
const double crhs64 =             crhs19*(crhs18 + crhs58 + 0.5854102*v(0,2));
const double crhs65 =             0.08090169924532*phi[3];
const double crhs66 =             -0.08090169924532*phi_old[3];
const double crhs67 =             0.01909830025156*phi[1];
const double crhs68 =             -0.01909830025156*phi_old[1];
const double crhs69 =             crhs22*(crhs25 + crhs26 + crhs28 + crhs29 + crhs65 + crhs66 + crhs67 + crhs68);
const double crhs70 =             0.08090169924532*phi[2];
const double crhs71 =             -0.08090169924532*phi_old[2];
const double crhs72 =             crhs22*(crhs25 + crhs27 + crhs28 + crhs30 + crhs67 + crhs68 + crhs70 + crhs71);
const double crhs73 =             crhs65 + crhs66 + crhs70 + crhs71;
const double crhs74 =             0.1381966*phi[1];
const double crhs75 =             crhs32 + crhs74;
const double crhs76 =             crhs36*(crhs33 + crhs75 + 0.5854102*phi[3]);
const double crhs77 =             0.1381966*crhs76;
const double crhs78 =             crhs36*(crhs34 + crhs75 + 0.5854102*phi[2]);
const double crhs79 =             0.1381966*crhs78;
const double crhs80 =             crhs36*(crhs35 + crhs74 + 0.5854102*phi[0]);
const double crhs81 =             -0.19999999899376*f[0];
const double crhs82 =             0.1381966*crhs50;
const double crhs83 =             0.1381966*crhs57;
const double crhs84 =             0.1381966*crhs64;
const double crhs85 =             0.08090169924532*phi[0];
const double crhs86 =             -0.08090169924532*phi_old[0];
const double crhs87 =             crhs22*(crhs26 + crhs27 + crhs29 + crhs30 + crhs67 + crhs68 + crhs85 + crhs86);
const double crhs88 =             0.1381966*crhs80;
const double crhs89 =             crhs0 + crhs14 + crhs21 + crhs31 + crhs38 + crhs7 + crhs81 + crhs82 + crhs83 + crhs84 + crhs87 + crhs88;
const double crhs90 =             crhs23 + crhs24 + crhs85 + crhs86;
            rhs[0]=DN_DX_0_0*crhs41 + DN_DX_0_1*crhs42 + DN_DX_0_2*crhs43 + crhs0 + crhs14 + crhs21 + crhs22*(crhs23 + crhs24 + crhs73 + 0.34270510226404*phi[0] - 0.34270510226404*phi_old[0]) + crhs31 + crhs38 + crhs39 + crhs40 + crhs47 + crhs49 + 0.5854102*crhs50 + crhs54 + crhs56 + 0.5854102*crhs57 + crhs61 + crhs63 + 0.5854102*crhs64 + crhs69 + crhs7 + crhs72 + crhs77 + crhs79 + 0.5854102*crhs80 - 0.40000000301872*f[0];
            rhs[1]=DN_DX_1_0*crhs41 + DN_DX_1_1*crhs42 + DN_DX_1_2*crhs43 + 0.5854102*crhs13 + 0.5854102*crhs20 + crhs22*(crhs73 + crhs85 + crhs86 + 0.34270510226404*phi[1] - 0.34270510226404*phi_old[1]) + 0.5854102*crhs37 + crhs39 + crhs40 + crhs47 + crhs49 + crhs54 + crhs56 + 0.5854102*crhs6 + crhs61 + crhs63 + crhs69 + crhs72 + crhs77 + crhs79 + crhs81 + crhs82 + crhs83 + crhs84 + crhs87 + crhs88 - 0.40000000301872*f[1];
            rhs[2]=DN_DX_2_0*crhs41 + DN_DX_2_1*crhs42 + DN_DX_2_2*crhs43 + crhs22*(crhs65 + crhs66 + crhs90 + 0.34270510226404*phi[2] - 0.34270510226404*phi_old[2]) + crhs47 + 0.5854102*crhs48 + crhs54 + 0.5854102*crhs55 + crhs61 + 0.5854102*crhs62 + crhs69 + crhs77 + 0.5854102*crhs78 + crhs89 - 0.40000000301872*f[2] - 0.19999999899376*f[3];
            rhs[3]=DN_DX_3_0*crhs41 + DN_DX_3_1*crhs42 + DN_DX_3_2*crhs43 + crhs22*(crhs70 + crhs71 + crhs90 + 0.34270510226404*phi[3] - 0.34270510226404*phi_old[3]) + 0.5854102*crhs46 + crhs49 + 0.5854102*crhs53 + crhs56 + 0.5854102*crhs60 + crhs63 + crhs72 + 0.5854102*crhs76 + crhs79 + crhs89 - 0.19999999899376*f[2] - 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
double SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::ComputeH(
    BoundedMatrix<double,TNumNodes,TDim >& DN_DX)
{
    KRATOS_TRY;

    double h=0.0;
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        double h_inv = 0.0;
        for(unsigned int k=0; k<TDim; k++)
        {
            h_inv += DN_DX(i,k)*DN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    h = sqrt(h)/static_cast<double>(TNumNodes);
    return h;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
    ElementData& rData)
{
    KRATOS_TRY;

    // Calculate h
    double h = this->ComputeH(rData.DN_DX);
    // Calculate tau for each gauss point
    for(unsigned int g = 0; g<TNumNodes; g++) {
	const auto& N = row(rData.N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        const array_1d<double,3> vel_gauss = prod(N, rData.convective_velocity);
        double div_vel = 0;
        for(unsigned int node_element = 0; node_element<TNumNodes; node_element++) {
            for(unsigned int dim = 0; dim < TDim; dim++) {
                div_vel += rData.DN_DX(node_element,dim)*rData.convective_velocity(node_element,dim);
            }
        }
        const double norm_velocity = norm_2(vel_gauss);
        // Estimate tau
        // Dynamic part
        double inv_tau = rData.dynamic_tau/rData.delta_time;
        // Convection
        inv_tau += 2.0 * norm_velocity / h;
        inv_tau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
        // Diffusion
        inv_tau += 4.0 * rData.diffusivity / (h*h);
        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);
        rData.tau[g] = (1.0) / inv_tau;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::IntegrationMethod SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicQSConvectionDiffusionExplicit<2,3>;
template class SymbolicQSConvectionDiffusionExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}
