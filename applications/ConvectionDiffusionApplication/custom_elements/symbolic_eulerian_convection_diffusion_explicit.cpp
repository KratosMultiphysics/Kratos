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

#include "symbolic_eulerian_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicEulerianConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != local_size)
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    if (rRightHandSideVector.size() != local_size)
        rRightHandSideVector.resize(local_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Define local variables
    Element::GeometryType::JacobiansType J0;
    Matrix InvJ0(dimension,dimension);
    double DetJ0;

    // Compute Jacobian
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g],InvJ0,DetJ0);

        // Calculate the cartesian derivatives on integration point "g"
        rVariables.DN = prod(DN_De[g],InvJ0);
        // Caluclate N on the gauss point "g"
        rVariables.N = row(N_gausspoint,g);
        // Compute weight
        rVariables.weight = integration_points[g].Weight() * DetJ0;
        // Compute tau
        this->CalculateTau(rVariables);

        // Execute standard RHS-LHS build or OSS step
        if (rCurrentProcessInfo.GetValue(ACTIVATION_LEVEL) == 1)
        {
            // Update OSS additional term
            this->ComputeOSSGaussPointContribution(rVariables,rLeftHandSideMatrix,rRightHandSideVector);
        }
        else
        {
            // Update rhs and lhs
            this->ComputeGaussPointContribution(rVariables,rLeftHandSideMatrix,rRightHandSideVector);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    Matrix LeftHandSide;
    this->CalculateLocalSystem(LeftHandSide,rRightHandSideVector,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    VectorType RightHandSide;
    this->CalculateLocalSystem(rLeftHandSideMatrix,RightHandSide,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if(rResult.size() != local_size)
        rResult.resize(local_size);
    for (unsigned int i=0; i<local_size; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& r_unknown_var = p_settings->GetUnknownVariable();
    unsigned int local_size = GetGeometry().PointsNumber();
    if (ElementalDofList.size() != local_size)
        ElementalDofList.resize(local_size);
    for (unsigned int i = 0; i < local_size; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    ProcessInfo &rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;

    auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    // Execute RK4 or OSS step
    if (r_process_info.GetValue(ACTIVATION_LEVEL) == 1)
    {
        // Calculate the OSS additional term
        VectorType res;
        this->CalculateRightHandSide(res,rCurrentProcessInfo);
        // Add the oss contribution to the forcing
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable()) += res[i_node];
        }
    }
    else
    {
        // Calculate the explicit residual vector
        VectorType rhs;
        this->CalculateRightHandSide(rhs,rCurrentProcessInfo);
        // Add the residual contribution
        // Note that the reaction is indeed the formulation residual
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable()) += rhs[i_node];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<2>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
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
}

/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
int SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::InitializeEulerianElement(
    ElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();

    // initialize some scalar variables
    rVariables.lumping_factor = 1.00 / double(TNumNodes);
    rVariables.diffusivity = 0.0;
    rVariables.specific_heat = 0.0;
    rVariables.density = 0.0;
    rVariables.dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

    for(unsigned int node_element = 0; node_element<local_size; node_element++)
{
    // observations
    // * ASGS time derivative term approximated as (phi-phi_old)/(RK_time_coefficient*delta_time)
    //   observe that for RK step = 1 ASGS time derivative term = 0 because phi = phi_old
    // * convective velocity and forcing term:
    //   RK step 1: evaluated at previous time step
    //   RK steps 2 and 3: linear interpolation between current and oldprevious time step
    //   RK step 4: evaluated at current time step
    // * convective_velocity = velocity - velocity_mesh
    //   velocity_mesh = 0 in eulerian framework
    if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==1)
    {
        rVariables.RK_time_coefficient = 0.5;
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable(),1);
        rVariables.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[0];
        rVariables.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[1];
        rVariables.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[2];
    }
    else if (r_process_info.GetValue(RUNGE_KUTTA_STEP)==2 || r_process_info.GetValue(RUNGE_KUTTA_STEP)==3)
    {
        rVariables.RK_time_coefficient = 0.5;
        rVariables.forcing[node_element] = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable()) + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable(),1));
        rVariables.convective_velocity(node_element,0) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[0] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[0]);
        rVariables.convective_velocity(node_element,1) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[1] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[1]);
        rVariables.convective_velocity(node_element,2) = 0.5*(r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable(),1)[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable(),1)[2] + r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[2]);
    }
    else
    {
        rVariables.RK_time_coefficient = 1.0;
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable());
        rVariables.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[0];
        rVariables.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[1];
        rVariables.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[2];
    }
    rVariables.diffusivity += r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetDiffusionVariable());
    rVariables.specific_heat += r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetSpecificHeatVariable());
    rVariables.density += r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetDensityVariable());
    rVariables.delta_time = r_process_info[DELTA_TIME];
    rVariables.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable());
    rVariables.unknown_old[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable(),1);
}
    // divide by number of nodes scalar variables
    rVariables.diffusivity *= rVariables.lumping_factor;
    rVariables.density *= rVariables.lumping_factor;
    rVariables.specific_heat *= rVariables.lumping_factor;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<2>::ComputeGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double clhs1 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs2 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs3 =             DN(0,0)*clhs1 + DN(0,1)*clhs2;
const double clhs4 =             DN(1,0)*clhs1 + DN(1,1)*clhs2;
const double clhs5 =             DN(0,0)*k;
const double clhs6 =             DN(0,1)*k;
const double clhs7 =             N[0]*clhs0;
const double clhs8 =             DN(1,0)*clhs5 + DN(1,1)*clhs6 + N[1]*clhs7;
const double clhs9 =             DN(2,0)*clhs1 + DN(2,1)*clhs2;
const double clhs10 =             DN(2,0)*clhs5 + DN(2,1)*clhs6 + N[2]*clhs7;
const double clhs11 =             DN(1,0)*DN(2,0)*k + DN(1,1)*DN(2,1)*k + N[1]*N[2]*clhs0;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(N[0], 2)*clhs0 + N[0]*clhs3;
            lhs(0,1)=N[0]*clhs4 + clhs8;
            lhs(0,2)=N[0]*clhs9 + clhs10;
            lhs(1,0)=N[1]*clhs3 + clhs8;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(N[1], 2)*clhs0 + N[1]*clhs4;
            lhs(1,2)=N[1]*clhs9 + clhs11;
            lhs(2,0)=N[2]*clhs3 + clhs10;
            lhs(2,1)=N[2]*clhs4 + clhs11;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(N[2], 2)*clhs0 + N[2]*clhs9;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs4 =             crhs3*k;
const double crhs5 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double crhs6 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 + N[0]*crhs0 - N[0]*crhs5 - N[0]*crhs6;
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 + N[1]*crhs0 - N[1]*crhs5 - N[1]*crhs6;
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 + N[2]*crhs0 - N[2]*crhs5 - N[2]*crhs6;


    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<3>::ComputeGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double clhs1 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs2 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs3 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs4 =             DN(0,0)*clhs1 + DN(0,1)*clhs2 + DN(0,2)*clhs3;
const double clhs5 =             DN(1,0)*clhs1 + DN(1,1)*clhs2 + DN(1,2)*clhs3;
const double clhs6 =             DN(0,0)*k;
const double clhs7 =             DN(0,1)*k;
const double clhs8 =             DN(0,2)*k;
const double clhs9 =             N[0]*clhs0;
const double clhs10 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8 + N[1]*clhs9;
const double clhs11 =             DN(2,0)*clhs1 + DN(2,1)*clhs2 + DN(2,2)*clhs3;
const double clhs12 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8 + N[2]*clhs9;
const double clhs13 =             DN(3,0)*clhs1 + DN(3,1)*clhs2 + DN(3,2)*clhs3;
const double clhs14 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8 + N[3]*clhs9;
const double clhs15 =             DN(1,0)*k;
const double clhs16 =             DN(1,1)*k;
const double clhs17 =             DN(1,2)*k;
const double clhs18 =             N[1]*clhs0;
const double clhs19 =             DN(2,0)*clhs15 + DN(2,1)*clhs16 + DN(2,2)*clhs17 + N[2]*clhs18;
const double clhs20 =             DN(3,0)*clhs15 + DN(3,1)*clhs16 + DN(3,2)*clhs17 + N[3]*clhs18;
const double clhs21 =             DN(2,0)*DN(3,0)*k + DN(2,1)*DN(3,1)*k + DN(2,2)*DN(3,2)*k + N[2]*N[3]*clhs0;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(DN(0,2), 2)*k + pow(N[0], 2)*clhs0 + N[0]*clhs4;
            lhs(0,1)=N[0]*clhs5 + clhs10;
            lhs(0,2)=N[0]*clhs11 + clhs12;
            lhs(0,3)=N[0]*clhs13 + clhs14;
            lhs(1,0)=N[1]*clhs4 + clhs10;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(DN(1,2), 2)*k + pow(N[1], 2)*clhs0 + N[1]*clhs5;
            lhs(1,2)=N[1]*clhs11 + clhs19;
            lhs(1,3)=N[1]*clhs13 + clhs20;
            lhs(2,0)=N[2]*clhs4 + clhs12;
            lhs(2,1)=N[2]*clhs5 + clhs19;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(DN(2,2), 2)*k + pow(N[2], 2)*clhs0 + N[2]*clhs11;
            lhs(2,3)=N[2]*clhs13 + clhs21;
            lhs(3,0)=N[3]*clhs4 + clhs14;
            lhs(3,1)=N[3]*clhs5 + clhs20;
            lhs(3,2)=N[3]*clhs11 + clhs21;
            lhs(3,3)=pow(DN(3,0), 2)*k + pow(DN(3,1), 2)*k + pow(DN(3,2), 2)*k + pow(N[3], 2)*clhs0 + N[3]*clhs13;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs6 =             crhs5*k;
const double crhs7 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double crhs8 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) + crhs5*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 - DN(0,2)*crhs6 + N[0]*crhs0 - N[0]*crhs7 - N[0]*crhs8;
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 - DN(1,2)*crhs6 + N[1]*crhs0 - N[1]*crhs7 - N[1]*crhs8;
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 - DN(2,2)*crhs6 + N[2]*crhs0 - N[2]*crhs7 - N[2]*crhs8;
            rhs[3]=-DN(3,0)*crhs2 - DN(3,1)*crhs4 - DN(3,2)*crhs6 + N[3]*crhs0 - N[3]*crhs7 - N[3]*crhs8;


    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<2>::ComputeOSSGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs4 =             crhs3*k;
const double crhs5 =             N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2];
const double crhs6 =             (-N[0]*phi_old[0] - N[1]*phi_old[1] - N[2]*phi_old[2] + crhs5)/(RK_time_coefficient*delta_time);
const double crhs7 =             crhs5*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double crhs8 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
            rhs[0]=DN(0,0)*crhs2 + DN(0,1)*crhs4 - N[0]*crhs0 + N[0]*crhs6 + N[0]*crhs7 + N[0]*crhs8;
            rhs[1]=DN(1,0)*crhs2 + DN(1,1)*crhs4 - N[1]*crhs0 + N[1]*crhs6 + N[1]*crhs7 + N[1]*crhs8;
            rhs[2]=DN(2,0)*crhs2 + DN(2,1)*crhs4 - N[2]*crhs0 + N[2]*crhs6 + N[2]*crhs7 + N[2]*crhs8;


    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicEulerianConvectionDiffusionExplicit<3>::ComputeOSSGaussPointContribution(
    ElementVariables& rVariables,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    // Retrieve element variables
    const auto N = rVariables.N;
    const auto DN = rVariables.DN;
    const auto k = rVariables.diffusivity;
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto phi_old = rVariables.unknown_old;
    const auto delta_time = rVariables.delta_time;
    const auto RK_time_coefficient = rVariables.RK_time_coefficient;
    const auto v = rVariables.convective_velocity;
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs6 =             crhs5*k;
const double crhs7 =             N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3];
const double crhs8 =             (-N[0]*phi_old[0] - N[1]*phi_old[1] - N[2]*phi_old[2] - N[3]*phi_old[3] + crhs7)/(RK_time_coefficient*delta_time);
const double crhs9 =             crhs7*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double crhs10 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) + crhs5*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
            rhs[0]=DN(0,0)*crhs2 + DN(0,1)*crhs4 + DN(0,2)*crhs6 - N[0]*crhs0 + N[0]*crhs10 + N[0]*crhs8 + N[0]*crhs9;
            rhs[1]=DN(1,0)*crhs2 + DN(1,1)*crhs4 + DN(1,2)*crhs6 - N[1]*crhs0 + N[1]*crhs10 + N[1]*crhs8 + N[1]*crhs9;
            rhs[2]=DN(2,0)*crhs2 + DN(2,1)*crhs4 + DN(2,2)*crhs6 - N[2]*crhs0 + N[2]*crhs10 + N[2]*crhs8 + N[2]*crhs9;
            rhs[3]=DN(3,0)*crhs2 + DN(3,1)*crhs4 + DN(3,2)*crhs6 - N[3]*crhs0 + N[3]*crhs10 + N[3]*crhs8 + N[3]*crhs9;


    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
double SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::ComputeH(
    BoundedMatrix<double,TNumNodes,TDim >& DN_DX)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
    ElementVariables& rVariables)
{
    // Calculate h
    double h = this->ComputeH(rVariables.DN);
    // Calculate velocity and velocity divergence in the gauss point
    array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
    double div_vel = 0;
    for(unsigned int node_element = 0; node_element<TNumNodes; node_element++)
    {
        for(unsigned int dim = 0; dim < TDim; dim++)
        {
            noalias(vel_gauss) = prod(rVariables.N,rVariables.convective_velocity);
            div_vel += rVariables.DN(node_element,dim)*rVariables.convective_velocity(node_element,dim);
        }
    }
    const double norm_velocity = norm_2(vel_gauss);
    // Estimate tau
    double inv_tau = 0;
    // Dynamic part
    inv_tau += rVariables.dynamic_tau * 1.0/rVariables.delta_time;
    // Convection
    inv_tau += 2.0 * norm_velocity / h;
    inv_tau += 1.0*div_vel; // unitary coefficient in fron of \nabla \cdot convective_velocity term in the strong equation
    // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
    inv_tau *= rVariables.density * rVariables.specific_heat;
    // Diffusion
    inv_tau += 4.0 * rVariables.diffusivity / (h*h);
    // Limiting
    inv_tau = std::max(inv_tau, 1e-2);
    rVariables.tau = (rVariables.density*rVariables.specific_heat) / inv_tau;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::IntegrationMethod SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}