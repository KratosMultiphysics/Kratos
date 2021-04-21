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

#include "qs_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
QSConvectionDiffusionExplicit<TDim,TNumNodes>::QSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
QSConvectionDiffusionExplicit<TDim,TNumNodes>::QSConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
QSConvectionDiffusionExplicit<TDim,TNumNodes>::~QSConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer QSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer QSConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit Convection-Diffusion element.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit Convection-Diffusion element. Call the QSCalculateRightHandSideInternal() method instead.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::EquationIdVector(
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
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::GetDofList(
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
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    BoundedVector<double, TNumNodes> rhs;
    this->QSCalculateRightHandSideInternal(rhs,rCurrentProcessInfo);
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
void QSConvectionDiffusionExplicit<2,3>::CalculateMassMatrix(
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
void QSConvectionDiffusionExplicit<3,4>::CalculateMassMatrix(
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
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Initialize the lumped mass vector
    if (rLumpedMassVector.size() != TNumNodes) {
        rLumpedMassVector.resize(TNumNodes, false);
    }

    // Fill the lumped mass vector
    const double nodal_mass = GetGeometry().DomainSize() / TNumNodes;
    std::fill(rLumpedMassVector.begin(),rLumpedMassVector.end(),nodal_mass);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    if (rVariable == p_settings->GetProjectionVariable()) {
        auto& r_geometry = GetGeometry();
        const unsigned int local_size = r_geometry.size();
        BoundedVector<double, TNumNodes> rhs_oss;
        this->QSCalculateOrthogonalSubgridScaleRHSInternal(rhs_oss,rCurrentProcessInfo);
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
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::InitializeEulerianElement(
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
    rData.delta_time = rCurrentProcessInfo[DELTA_TIME];
    double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    if (theta == 0.0) {
        rData.explicit_step_coefficient = 0.0;
    }
    else {
        rData.explicit_step_coefficient = 1.0/((theta)*r_process_info[DELTA_TIME]);
    }

    for(unsigned int node_element = 0; node_element<local_size; node_element++) {
        // Observations
        // * unknown acceleration approximated as (phi-phi_old)*explicit_step_coefficient = (phi-phi_old)/((theta)*delta_time)
        //   observe that for theta = 0.0, phi = phi_old and explicit_step_coefficient = 0
        // * convective velocity and forcing term:
        //   interpolation exploiting theta
        // * convective_velocity = velocity - velocity_mesh
        //   velocity_mesh = 0 in eulerian framework
        rData.forcing[node_element] = (1-theta) * r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable(),1) + theta * r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVolumeSourceVariable());
        rData.convective_velocity(node_element,0) = (1-theta) * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[0]) + theta * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[0]);
        rData.convective_velocity(node_element,1) = (1-theta) * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[1]) + theta * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[1]);
        rData.convective_velocity(node_element,2) = (1-theta) * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable(),1)[2]) + theta * (r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetMeshVelocityVariable())[2]);

        rData.oss_projection[node_element] = r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetProjectionVariable());
        rData.diffusivity += r_geometry[node_element].FastGetSolutionStepValue(p_settings->GetDiffusionVariable());
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
void QSConvectionDiffusionExplicit<2,3>::QSCalculateRightHandSideInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->QSCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
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
const double crhs55 =             -0.166666666666667*phi_old[1];
const double crhs56 =             -0.166666666666667*phi_old[2];
const double crhs57 =             explicit_step_coefficient*(crhs38 + crhs55 + crhs56 - 0.666666666666667*phi_old[0]);
const double crhs58 =             -0.166666666666667*phi_old[0];
const double crhs59 =             explicit_step_coefficient*(crhs15 + crhs56 + crhs58 - 0.666666666666667*phi_old[1]);
const double crhs60 =             explicit_step_coefficient*(crhs35 + crhs55 + crhs58 - 0.666666666666667*phi_old[2]);
const double crhs61 =             crhs27 + crhs33;
const double crhs62 =             crhs11 + crhs5;
const double crhs63 =             crhs24 + crhs30;
const double crhs64 =             -0.166666666666667*crhs27 - 0.166666666666667*crhs33 - 0.166666666666667*crhs39 + 0.25*f[0];
const double crhs65 =             tau[0]*(DN_DX_1_0*crhs26 + DN_DX_1_1*crhs32);
const double crhs66 =             tau[1]*(DN_DX_1_0*crhs3 + DN_DX_1_1*crhs9);
const double crhs67 =             tau[2]*(DN_DX_1_0*crhs23 + DN_DX_1_1*crhs29);
const double crhs68 =             tau[0]*(DN_DX_2_0*crhs26 + DN_DX_2_1*crhs32);
const double crhs69 =             tau[1]*(DN_DX_2_0*crhs3 + DN_DX_2_1*crhs9);
const double crhs70 =             tau[2]*(DN_DX_2_0*crhs23 + DN_DX_2_1*crhs29);
            rhs[0]=-DN_DX_0_0*crhs20 - DN_DX_0_1*crhs21 + crhs0 + crhs12 - crhs17*crhs49 + crhs18 + crhs19 + crhs25 - 0.666666666666667*crhs27 + crhs31 - 0.666666666666667*crhs33 - crhs36*crhs53 + crhs37 - crhs39*crhs43 - 0.666666666666667*crhs39 + crhs42*crhs43 + crhs43*crhs46 - crhs43*crhs57 - crhs43*crhs61 + crhs48*crhs49 + crhs49*crhs51 - crhs49*crhs59 - crhs49*crhs62 + crhs52*crhs53 + crhs53*crhs54 - crhs53*crhs60 - crhs53*crhs63 + crhs6 + 0.5*f[0];
            rhs[1]=-DN_DX_1_0*crhs20 - DN_DX_1_1*crhs21 - 0.666666666666667*crhs11 - crhs17*crhs66 - 0.666666666666667*crhs17 + crhs19 + crhs25 + crhs31 - crhs36*crhs67 + crhs37 - crhs39*crhs65 + crhs42*crhs65 + crhs46*crhs65 + crhs48*crhs66 - 0.666666666666667*crhs5 + crhs51*crhs66 + crhs52*crhs67 + crhs54*crhs67 - crhs57*crhs65 - crhs59*crhs66 - crhs60*crhs67 - crhs61*crhs65 - crhs62*crhs66 - crhs63*crhs67 + crhs64 + 0.5*f[1];
            rhs[2]=-DN_DX_2_0*crhs20 - DN_DX_2_1*crhs21 + crhs0 + crhs12 - crhs17*crhs69 + crhs18 - 0.666666666666667*crhs24 - 0.666666666666667*crhs30 - crhs36*crhs70 - 0.666666666666667*crhs36 - crhs39*crhs68 + crhs42*crhs68 + crhs46*crhs68 + crhs48*crhs69 + crhs51*crhs69 + crhs52*crhs70 + crhs54*crhs70 - crhs57*crhs68 - crhs59*crhs69 + crhs6 - crhs60*crhs70 - crhs61*crhs68 - crhs62*crhs69 - crhs63*crhs70 + crhs64 + 0.5*f[2];


    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void QSConvectionDiffusionExplicit<3,4>::QSCalculateRightHandSideInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->QSCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
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
const double crhs102 =             -0.1381966*phi_old[2];
const double crhs103 =             -0.1381966*phi_old[3];
const double crhs104 =             crhs102 + crhs103;
const double crhs105 =             -0.1381966*phi_old[1];
const double crhs106 =             explicit_step_coefficient*(crhs104 + crhs105 + crhs27 + crhs28 + crhs68 + crhs76 - 0.5854102*phi_old[0]);
const double crhs107 =             -0.1381966*phi_old[0];
const double crhs108 =             explicit_step_coefficient*(crhs104 + crhs107 + crhs25 + crhs26 + crhs27 + crhs28 - 0.5854102*phi_old[1]);
const double crhs109 =             crhs105 + crhs107;
const double crhs110 =             explicit_step_coefficient*(crhs103 + crhs109 + crhs25 + crhs28 + crhs68 + crhs73 - 0.5854102*phi_old[2]);
const double crhs111 =             explicit_step_coefficient*(crhs102 + crhs109 + crhs25 + crhs27 + crhs68 + crhs70 - 0.5854102*phi_old[3]);
const double crhs112 =             crhs47 + crhs57 + crhs67;
const double crhs113 =             crhs15 + crhs23 + crhs7;
const double crhs114 =             crhs44 + crhs54 + crhs64;
const double crhs115 =             crhs41 + crhs51 + crhs61;
const double crhs116 =             0.19999999899376*f[0];
const double crhs117 =             -0.1381966*crhs47;
const double crhs118 =             -0.1381966*crhs57;
const double crhs119 =             -0.1381966*crhs67;
const double crhs120 =             -0.1381966*crhs77;
const double crhs121 =             tau[0]*(DN_DX_1_0*crhs46 + DN_DX_1_1*crhs56 + DN_DX_1_2*crhs66);
const double crhs122 =             tau[1]*(DN_DX_1_0*crhs5 + DN_DX_1_1*crhs13 + DN_DX_1_2*crhs21);
const double crhs123 =             tau[2]*(DN_DX_1_0*crhs43 + DN_DX_1_1*crhs53 + DN_DX_1_2*crhs63);
const double crhs124 =             tau[3]*(DN_DX_1_0*crhs40 + DN_DX_1_1*crhs50 + DN_DX_1_2*crhs60);
const double crhs125 =             crhs0 + crhs116 + crhs117 + crhs118 + crhs119 + crhs120 + crhs16 + crhs24 + crhs32 + crhs8;
const double crhs126 =             tau[0]*(DN_DX_2_0*crhs46 + DN_DX_2_1*crhs56 + DN_DX_2_2*crhs66);
const double crhs127 =             tau[1]*(DN_DX_2_0*crhs5 + DN_DX_2_1*crhs13 + DN_DX_2_2*crhs21);
const double crhs128 =             tau[2]*(DN_DX_2_0*crhs43 + DN_DX_2_1*crhs53 + DN_DX_2_2*crhs63);
const double crhs129 =             tau[3]*(DN_DX_2_0*crhs40 + DN_DX_2_1*crhs50 + DN_DX_2_2*crhs60);
const double crhs130 =             tau[0]*(DN_DX_3_0*crhs46 + DN_DX_3_1*crhs56 + DN_DX_3_2*crhs66);
const double crhs131 =             tau[1]*(DN_DX_3_0*crhs5 + DN_DX_3_1*crhs13 + DN_DX_3_2*crhs21);
const double crhs132 =             tau[2]*(DN_DX_3_0*crhs43 + DN_DX_3_1*crhs53 + DN_DX_3_2*crhs63);
const double crhs133 =             tau[3]*(DN_DX_3_0*crhs40 + DN_DX_3_1*crhs50 + DN_DX_3_2*crhs60);
            rhs[0]=-DN_DX_0_0*crhs35 - DN_DX_0_1*crhs36 - DN_DX_0_2*crhs37 + crhs0 + crhs100*crhs101 - crhs100*crhs111 - crhs100*crhs115 - crhs100*crhs71 + crhs100*crhs99 - crhs106*crhs83 - crhs108*crhs91 - crhs110*crhs96 - crhs112*crhs83 - crhs113*crhs91 - crhs114*crhs96 + crhs16 + crhs24 - crhs31*crhs91 + crhs32 + crhs33 + crhs34 + crhs42 + crhs45 - 0.5854102*crhs47 + crhs52 + crhs55 - 0.5854102*crhs57 + crhs62 + crhs65 - 0.5854102*crhs67 + crhs72 - crhs74*crhs96 + crhs75 - crhs77*crhs83 - 0.5854102*crhs77 + crhs8 + crhs82*crhs83 + crhs83*crhs88 + crhs90*crhs91 + crhs91*crhs93 + crhs95*crhs96 + crhs96*crhs98 + 0.40000000301872*f[0];
            rhs[1]=-DN_DX_1_0*crhs35 - DN_DX_1_1*crhs36 - DN_DX_1_2*crhs37 + crhs101*crhs124 - crhs106*crhs121 - crhs108*crhs122 - crhs110*crhs123 - crhs111*crhs124 - crhs112*crhs121 - crhs113*crhs122 - crhs114*crhs123 - crhs115*crhs124 + crhs116 + crhs117 + crhs118 + crhs119 + crhs120 - crhs121*crhs77 + crhs121*crhs82 + crhs121*crhs88 - crhs122*crhs31 + crhs122*crhs90 + crhs122*crhs93 - crhs123*crhs74 + crhs123*crhs95 + crhs123*crhs98 - crhs124*crhs71 + crhs124*crhs99 - 0.5854102*crhs15 - 0.5854102*crhs23 - 0.5854102*crhs31 + crhs33 + crhs34 + crhs42 + crhs45 + crhs52 + crhs55 + crhs62 + crhs65 - 0.5854102*crhs7 + crhs72 + crhs75 + 0.40000000301872*f[1];
            rhs[2]=-DN_DX_2_0*crhs35 - DN_DX_2_1*crhs36 - DN_DX_2_2*crhs37 + crhs101*crhs129 - crhs106*crhs126 - crhs108*crhs127 - crhs110*crhs128 - crhs111*crhs129 - crhs112*crhs126 - crhs113*crhs127 - crhs114*crhs128 - crhs115*crhs129 + crhs125 - crhs126*crhs77 + crhs126*crhs82 + crhs126*crhs88 - crhs127*crhs31 + crhs127*crhs90 + crhs127*crhs93 - crhs128*crhs74 + crhs128*crhs95 + crhs128*crhs98 - crhs129*crhs71 + crhs129*crhs99 + crhs42 - 0.5854102*crhs44 + crhs52 - 0.5854102*crhs54 + crhs62 - 0.5854102*crhs64 + crhs72 - 0.5854102*crhs74 + 0.40000000301872*f[2] + 0.19999999899376*f[3];
            rhs[3]=-DN_DX_3_0*crhs35 - DN_DX_3_1*crhs36 - DN_DX_3_2*crhs37 + crhs101*crhs133 - crhs106*crhs130 - crhs108*crhs131 - crhs110*crhs132 - crhs111*crhs133 - crhs112*crhs130 - crhs113*crhs131 - crhs114*crhs132 - crhs115*crhs133 + crhs125 - crhs130*crhs77 + crhs130*crhs82 + crhs130*crhs88 - crhs131*crhs31 + crhs131*crhs90 + crhs131*crhs93 - crhs132*crhs74 + crhs132*crhs95 + crhs132*crhs98 - crhs133*crhs71 + crhs133*crhs99 - 0.5854102*crhs41 + crhs45 - 0.5854102*crhs51 + crhs55 - 0.5854102*crhs61 + crhs65 - 0.5854102*crhs71 + crhs75 + 0.19999999899376*f[2] + 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void QSConvectionDiffusionExplicit<2,3>::QSCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->QSCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
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
const double crhs1 =             0.166666666666667*phi[0];
const double crhs2 =             0.166666666666667*phi[2];
const double crhs3 =             crhs1 + crhs2 + 0.666666666666667*phi[1];
const double crhs4 =             -0.166666666666667*phi_old[0];
const double crhs5 =             -0.166666666666667*phi_old[2];
const double crhs6 =             explicit_step_coefficient*(crhs3 + crhs4 + crhs5 - 0.666666666666667*phi_old[1]);
const double crhs7 =             0.166666666666667*crhs6;
const double crhs8 =             0.166666666666667*v(0,0);
const double crhs9 =             0.166666666666667*v(2,0);
const double crhs10 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2];
const double crhs11 =             crhs10*(crhs8 + crhs9 + 0.666666666666667*v(1,0));
const double crhs12 =             0.166666666666667*crhs11;
const double crhs13 =             0.166666666666667*v(0,1);
const double crhs14 =             0.166666666666667*v(2,1);
const double crhs15 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2];
const double crhs16 =             crhs15*(crhs13 + crhs14 + 0.666666666666667*v(1,1));
const double crhs17 =             0.166666666666667*crhs16;
const double crhs18 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1);
const double crhs19 =             crhs18*crhs3;
const double crhs20 =             0.166666666666667*crhs19;
const double crhs21 =             -0.25*f[2];
const double crhs22 =             3*alpha*crhs10;
const double crhs23 =             3*alpha*crhs15;
const double crhs24 =             0.166666666666667*phi[1];
const double crhs25 =             crhs1 + crhs24 + 0.666666666666667*phi[2];
const double crhs26 =             -0.166666666666667*phi_old[1];
const double crhs27 =             explicit_step_coefficient*(crhs25 + crhs26 + crhs4 - 0.666666666666667*phi_old[2]);
const double crhs28 =             0.166666666666667*crhs27;
const double crhs29 =             crhs2 + crhs24 + 0.666666666666667*phi[0];
const double crhs30 =             explicit_step_coefficient*(crhs26 + crhs29 + crhs5 - 0.666666666666667*phi_old[0]);
const double crhs31 =             0.166666666666667*v(1,0);
const double crhs32 =             crhs10*(crhs31 + crhs8 + 0.666666666666667*v(2,0));
const double crhs33 =             0.166666666666667*crhs32;
const double crhs34 =             crhs10*(crhs31 + crhs9 + 0.666666666666667*v(0,0));
const double crhs35 =             0.166666666666667*v(1,1);
const double crhs36 =             crhs15*(crhs13 + crhs35 + 0.666666666666667*v(2,1));
const double crhs37 =             0.166666666666667*crhs36;
const double crhs38 =             crhs15*(crhs14 + crhs35 + 0.666666666666667*v(0,1));
const double crhs39 =             crhs18*crhs25;
const double crhs40 =             0.166666666666667*crhs39;
const double crhs41 =             crhs18*crhs29;
const double crhs42 =             0.166666666666667*crhs30 + 0.166666666666667*crhs34 + 0.166666666666667*crhs38 + 0.166666666666667*crhs41 - 0.25*f[0];
            rhs[0]=DN_DX_0_0*crhs22 + DN_DX_0_1*crhs23 + crhs0 + crhs12 + crhs17 + crhs20 + crhs21 + crhs28 + 0.666666666666667*crhs30 + crhs33 + 0.666666666666667*crhs34 + crhs37 + 0.666666666666667*crhs38 + crhs40 + 0.666666666666667*crhs41 + crhs7 - 0.5*f[0];
            rhs[1]=DN_DX_1_0*crhs22 + DN_DX_1_1*crhs23 + 0.666666666666667*crhs11 + 0.666666666666667*crhs16 + 0.666666666666667*crhs19 + crhs21 + crhs28 + crhs33 + crhs37 + crhs40 + crhs42 + 0.666666666666667*crhs6 - 0.5*f[1];
            rhs[2]=DN_DX_2_0*crhs22 + DN_DX_2_1*crhs23 + crhs0 + crhs12 + crhs17 + crhs20 + 0.666666666666667*crhs27 + 0.666666666666667*crhs32 + 0.666666666666667*crhs36 + 0.666666666666667*crhs39 + crhs42 + crhs7 - 0.5*f[2];


    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void QSConvectionDiffusionExplicit<3,4>::QSCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->QSCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
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
const double crhs1 =             -0.1381966*phi_old[2];
const double crhs2 =             -0.1381966*phi_old[3];
const double crhs3 =             crhs1 + crhs2;
const double crhs4 =             0.1381966*phi[0];
const double crhs5 =             0.5854102*phi[1];
const double crhs6 =             0.1381966*phi[2];
const double crhs7 =             0.1381966*phi[3];
const double crhs8 =             -0.1381966*phi_old[0];
const double crhs9 =             explicit_step_coefficient*(crhs3 + crhs4 + crhs5 + crhs6 + crhs7 + crhs8 - 0.5854102*phi_old[1]);
const double crhs10 =             0.1381966*crhs9;
const double crhs11 =             0.1381966*v(0,0);
const double crhs12 =             0.1381966*v(2,0);
const double crhs13 =             0.1381966*v(3,0);
const double crhs14 =             crhs12 + crhs13;
const double crhs15 =             DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3];
const double crhs16 =             crhs15*(crhs11 + crhs14 + 0.5854102*v(1,0));
const double crhs17 =             0.1381966*crhs16;
const double crhs18 =             0.1381966*v(0,1);
const double crhs19 =             0.1381966*v(2,1);
const double crhs20 =             0.1381966*v(3,1);
const double crhs21 =             crhs19 + crhs20;
const double crhs22 =             DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3];
const double crhs23 =             crhs22*(crhs18 + crhs21 + 0.5854102*v(1,1));
const double crhs24 =             0.1381966*crhs23;
const double crhs25 =             0.1381966*v(0,2);
const double crhs26 =             0.1381966*v(2,2);
const double crhs27 =             0.1381966*v(3,2);
const double crhs28 =             crhs26 + crhs27;
const double crhs29 =             DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3];
const double crhs30 =             crhs29*(crhs25 + crhs28 + 0.5854102*v(1,2));
const double crhs31 =             0.1381966*crhs30;
const double crhs32 =             crhs6 + crhs7;
const double crhs33 =             DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2);
const double crhs34 =             crhs33*(crhs32 + crhs4 + crhs5);
const double crhs35 =             0.1381966*crhs34;
const double crhs36 =             -0.19999999899376*f[2];
const double crhs37 =             -0.19999999899376*f[3];
const double crhs38 =             4*alpha*crhs15;
const double crhs39 =             4*alpha*crhs22;
const double crhs40 =             4*alpha*crhs29;
const double crhs41 =             -0.1381966*phi_old[1];
const double crhs42 =             crhs41 + crhs8;
const double crhs43 =             0.1381966*phi[1];
const double crhs44 =             0.5854102*phi[3];
const double crhs45 =             explicit_step_coefficient*(crhs1 + crhs4 + crhs42 + crhs43 + crhs44 + crhs6 - 0.5854102*phi_old[3]);
const double crhs46 =             0.1381966*crhs45;
const double crhs47 =             0.5854102*phi[2];
const double crhs48 =             explicit_step_coefficient*(crhs2 + crhs4 + crhs42 + crhs43 + crhs47 + crhs7 - 0.5854102*phi_old[2]);
const double crhs49 =             0.1381966*crhs48;
const double crhs50 =             0.5854102*phi[0];
const double crhs51 =             explicit_step_coefficient*(crhs3 + crhs41 + crhs43 + crhs50 + crhs6 + crhs7 - 0.5854102*phi_old[0]);
const double crhs52 =             0.1381966*v(1,0);
const double crhs53 =             crhs11 + crhs52;
const double crhs54 =             crhs15*(crhs12 + crhs53 + 0.5854102*v(3,0));
const double crhs55 =             0.1381966*crhs54;
const double crhs56 =             crhs15*(crhs13 + crhs53 + 0.5854102*v(2,0));
const double crhs57 =             0.1381966*crhs56;
const double crhs58 =             crhs15*(crhs14 + crhs52 + 0.5854102*v(0,0));
const double crhs59 =             0.1381966*v(1,1);
const double crhs60 =             crhs18 + crhs59;
const double crhs61 =             crhs22*(crhs19 + crhs60 + 0.5854102*v(3,1));
const double crhs62 =             0.1381966*crhs61;
const double crhs63 =             crhs22*(crhs20 + crhs60 + 0.5854102*v(2,1));
const double crhs64 =             0.1381966*crhs63;
const double crhs65 =             crhs22*(crhs21 + crhs59 + 0.5854102*v(0,1));
const double crhs66 =             0.1381966*v(1,2);
const double crhs67 =             crhs25 + crhs66;
const double crhs68 =             crhs29*(crhs26 + crhs67 + 0.5854102*v(3,2));
const double crhs69 =             0.1381966*crhs68;
const double crhs70 =             crhs29*(crhs27 + crhs67 + 0.5854102*v(2,2));
const double crhs71 =             0.1381966*crhs70;
const double crhs72 =             crhs29*(crhs28 + crhs66 + 0.5854102*v(0,2));
const double crhs73 =             crhs4 + crhs43;
const double crhs74 =             crhs33*(crhs44 + crhs6 + crhs73);
const double crhs75 =             0.1381966*crhs74;
const double crhs76 =             crhs33*(crhs47 + crhs7 + crhs73);
const double crhs77 =             0.1381966*crhs76;
const double crhs78 =             crhs33*(crhs32 + crhs43 + crhs50);
const double crhs79 =             -0.19999999899376*f[0];
const double crhs80 =             0.1381966*crhs51;
const double crhs81 =             0.1381966*crhs58;
const double crhs82 =             0.1381966*crhs65;
const double crhs83 =             0.1381966*crhs72;
const double crhs84 =             0.1381966*crhs78;
const double crhs85 =             crhs0 + crhs10 + crhs17 + crhs24 + crhs31 + crhs35 + crhs79 + crhs80 + crhs81 + crhs82 + crhs83 + crhs84;
            rhs[0]=DN_DX_0_0*crhs38 + DN_DX_0_1*crhs39 + DN_DX_0_2*crhs40 + crhs0 + crhs10 + crhs17 + crhs24 + crhs31 + crhs35 + crhs36 + crhs37 + crhs46 + crhs49 + 0.5854102*crhs51 + crhs55 + crhs57 + 0.5854102*crhs58 + crhs62 + crhs64 + 0.5854102*crhs65 + crhs69 + crhs71 + 0.5854102*crhs72 + crhs75 + crhs77 + 0.5854102*crhs78 - 0.40000000301872*f[0];
            rhs[1]=DN_DX_1_0*crhs38 + DN_DX_1_1*crhs39 + DN_DX_1_2*crhs40 + 0.5854102*crhs16 + 0.5854102*crhs23 + 0.5854102*crhs30 + 0.5854102*crhs34 + crhs36 + crhs37 + crhs46 + crhs49 + crhs55 + crhs57 + crhs62 + crhs64 + crhs69 + crhs71 + crhs75 + crhs77 + crhs79 + crhs80 + crhs81 + crhs82 + crhs83 + crhs84 + 0.5854102*crhs9 - 0.40000000301872*f[1];
            rhs[2]=DN_DX_2_0*crhs38 + DN_DX_2_1*crhs39 + DN_DX_2_2*crhs40 + crhs46 + 0.5854102*crhs48 + crhs55 + 0.5854102*crhs56 + crhs62 + 0.5854102*crhs63 + crhs69 + 0.5854102*crhs70 + crhs75 + 0.5854102*crhs76 + crhs85 - 0.40000000301872*f[2] - 0.19999999899376*f[3];
            rhs[3]=DN_DX_3_0*crhs38 + DN_DX_3_1*crhs39 + DN_DX_3_2*crhs40 + 0.5854102*crhs45 + crhs49 + 0.5854102*crhs54 + crhs57 + 0.5854102*crhs61 + crhs64 + 0.5854102*crhs68 + crhs71 + 0.5854102*crhs74 + crhs77 + crhs85 - 0.19999999899376*f[2] - 0.40000000301872*f[3];


    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
double QSConvectionDiffusionExplicit<TDim,TNumNodes>::ComputeH(
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
void QSConvectionDiffusionExplicit<TDim,TNumNodes>::QSCalculateTau(
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
Element::IntegrationMethod QSConvectionDiffusionExplicit<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template class QSConvectionDiffusionExplicit<2,3>;
template class QSConvectionDiffusionExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}
