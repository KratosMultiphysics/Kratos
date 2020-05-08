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

        // Update rhs and lhs
        this->ComputeGaussPointContribution(rVariables,rLeftHandSideMatrix,rRightHandSideVector);
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
    Vector RightHandSide;
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

    for(unsigned int node_element = 0; node_element<local_size; node_element++)
{
    rVariables.diffusivity[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetDiffusionVariable());
    rVariables.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable());
    rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable());
    // convective_velocity = velocity - velocity_mesh
    // velocity_mesh = 0 in eulerian framework
    rVariables.convective_velocity(node_element,0) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[0] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[0];
    rVariables.convective_velocity(node_element,1) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[1] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[1];
    rVariables.convective_velocity(node_element,2) = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVelocityVariable())[2] - r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetMeshVelocityVariable())[2];
}

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
    const auto k = inner_prod(N,rVariables.diffusivity);
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto v = rVariables.convective_velocity;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs1 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs2 =             DN(0,0)*clhs0 + DN(0,1)*clhs1;
const double clhs3 =             DN(1,0)*clhs0 + DN(1,1)*clhs1;
const double clhs4 =             DN(0,0)*k;
const double clhs5 =             DN(0,1)*k;
const double clhs6 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs7 =             DN(2,0)*clhs0 + DN(2,1)*clhs1;
const double clhs8 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs9 =             DN(1,0)*DN(2,0)*k + DN(1,1)*DN(2,1)*k;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + N[0]*clhs2;
            lhs(0,1)=N[0]*clhs3 + clhs6;
            lhs(0,2)=N[0]*clhs7 + clhs8;
            lhs(1,0)=N[1]*clhs2 + clhs6;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + N[1]*clhs3;
            lhs(1,2)=N[1]*clhs7 + clhs9;
            lhs(2,0)=N[2]*clhs2 + clhs8;
            lhs(2,1)=N[2]*clhs3 + clhs9;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + N[2]*clhs7;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs4 =             crhs3*k;
const double crhs5 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 + N[0]*crhs0 - N[0]*crhs5;
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 + N[1]*crhs0 - N[1]*crhs5;
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 + N[2]*crhs0 - N[2]*crhs5;


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
    const auto k = inner_prod(N,rVariables.diffusivity);
    const auto f = rVariables.forcing;
    const auto phi = rVariables.unknown;
    const auto v = rVariables.convective_velocity;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs1 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs2 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs3 =             DN(0,0)*clhs0 + DN(0,1)*clhs1 + DN(0,2)*clhs2;
const double clhs4 =             DN(1,0)*clhs0 + DN(1,1)*clhs1 + DN(1,2)*clhs2;
const double clhs5 =             DN(0,0)*k;
const double clhs6 =             DN(0,1)*k;
const double clhs7 =             DN(0,2)*k;
const double clhs8 =             DN(1,0)*clhs5 + DN(1,1)*clhs6 + DN(1,2)*clhs7;
const double clhs9 =             DN(2,0)*clhs0 + DN(2,1)*clhs1 + DN(2,2)*clhs2;
const double clhs10 =             DN(2,0)*clhs5 + DN(2,1)*clhs6 + DN(2,2)*clhs7;
const double clhs11 =             DN(3,0)*clhs0 + DN(3,1)*clhs1 + DN(3,2)*clhs2;
const double clhs12 =             DN(3,0)*clhs5 + DN(3,1)*clhs6 + DN(3,2)*clhs7;
const double clhs13 =             DN(1,0)*k;
const double clhs14 =             DN(1,1)*k;
const double clhs15 =             DN(1,2)*k;
const double clhs16 =             DN(2,0)*clhs13 + DN(2,1)*clhs14 + DN(2,2)*clhs15;
const double clhs17 =             DN(3,0)*clhs13 + DN(3,1)*clhs14 + DN(3,2)*clhs15;
const double clhs18 =             DN(2,0)*DN(3,0)*k + DN(2,1)*DN(3,1)*k + DN(2,2)*DN(3,2)*k;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(DN(0,2), 2)*k + N[0]*clhs3;
            lhs(0,1)=N[0]*clhs4 + clhs8;
            lhs(0,2)=N[0]*clhs9 + clhs10;
            lhs(0,3)=N[0]*clhs11 + clhs12;
            lhs(1,0)=N[1]*clhs3 + clhs8;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(DN(1,2), 2)*k + N[1]*clhs4;
            lhs(1,2)=N[1]*clhs9 + clhs16;
            lhs(1,3)=N[1]*clhs11 + clhs17;
            lhs(2,0)=N[2]*clhs3 + clhs10;
            lhs(2,1)=N[2]*clhs4 + clhs16;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(DN(2,2), 2)*k + N[2]*clhs9;
            lhs(2,3)=N[2]*clhs11 + clhs18;
            lhs(3,0)=N[3]*clhs3 + clhs12;
            lhs(3,1)=N[3]*clhs4 + clhs17;
            lhs(3,2)=N[3]*clhs9 + clhs18;
            lhs(3,3)=pow(DN(3,0), 2)*k + pow(DN(3,1), 2)*k + pow(DN(3,2), 2)*k + N[3]*clhs11;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs2 =             crhs1*k;
const double crhs3 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs6 =             crhs5*k;
const double crhs7 =             crhs1*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) + crhs3*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) + crhs5*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
            rhs[0]=-DN(0,0)*crhs2 - DN(0,1)*crhs4 - DN(0,2)*crhs6 + N[0]*crhs0 - N[0]*crhs7;
            rhs[1]=-DN(1,0)*crhs2 - DN(1,1)*crhs4 - DN(1,2)*crhs6 + N[1]*crhs0 - N[1]*crhs7;
            rhs[2]=-DN(2,0)*crhs2 - DN(2,1)*crhs4 - DN(2,2)*crhs6 + N[2]*crhs0 - N[2]*crhs7;
            rhs[3]=-DN(3,0)*crhs2 - DN(3,1)*crhs4 - DN(3,2)*crhs6 + N[3]*crhs0 - N[3]*crhs7;


    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::IntegrationMethod SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}