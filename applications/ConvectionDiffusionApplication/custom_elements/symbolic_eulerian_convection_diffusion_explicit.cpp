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
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs1 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs2 =             DN(0,0)*clhs0 + DN(0,1)*clhs1;
const double clhs3 =             N[0]*clhs2;
const double clhs4 =             tau*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double clhs5 =             DN(1,0)*clhs0 + DN(1,1)*clhs1;
const double clhs6 =             N[0]*clhs5;
const double clhs7 =             DN(0,0)*k;
const double clhs8 =             DN(0,1)*k;
const double clhs9 =             clhs2*tau;
const double clhs10 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + clhs5*clhs9;
const double clhs11 =             DN(2,0)*clhs0 + DN(2,1)*clhs1;
const double clhs12 =             N[0]*clhs11;
const double clhs13 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + clhs11*clhs9;
const double clhs14 =             N[1]*clhs2;
const double clhs15 =             N[1]*clhs5;
const double clhs16 =             N[1]*clhs11;
const double clhs17 =             DN(1,0)*DN(2,0)*k + DN(1,1)*DN(2,1)*k + clhs11*clhs5*tau;
const double clhs18 =             N[2]*clhs2;
const double clhs19 =             N[2]*clhs5;
const double clhs20 =             N[2]*clhs11;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(clhs2, 2)*tau + clhs3*clhs4 + clhs3;
            lhs(0,1)=clhs10 + clhs4*clhs6 + clhs6;
            lhs(0,2)=clhs12*clhs4 + clhs12 + clhs13;
            lhs(1,0)=clhs10 + clhs14*clhs4 + clhs14;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + clhs15*clhs4 + clhs15 + pow(clhs5, 2)*tau;
            lhs(1,2)=clhs16*clhs4 + clhs16 + clhs17;
            lhs(2,0)=clhs13 + clhs18*clhs4 + clhs18;
            lhs(2,1)=clhs17 + clhs19*clhs4 + clhs19;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(clhs11, 2)*tau + clhs20*clhs4 + clhs20;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs1 =             N[0]*crhs0;
const double crhs2 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs3 =             crhs2*k;
const double crhs4 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs5 =             crhs4*k;
const double crhs6 =             tau*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double crhs7 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs8 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs9 =             tau*(DN(0,0)*crhs7 + DN(0,1)*crhs8);
const double crhs10 =             crhs2*crhs7 + crhs4*crhs8;
const double crhs11 =             N[0]*crhs10;
const double crhs12 =             N[1]*crhs0;
const double crhs13 =             tau*(DN(1,0)*crhs7 + DN(1,1)*crhs8);
const double crhs14 =             N[1]*crhs10;
const double crhs15 =             N[2]*crhs0;
const double crhs16 =             tau*(DN(2,0)*crhs7 + DN(2,1)*crhs8);
const double crhs17 =             N[2]*crhs10;
            rhs[0]=-DN(0,0)*crhs3 - DN(0,1)*crhs5 + crhs0*crhs9 + crhs1*crhs6 + crhs1 - crhs10*crhs9 - crhs11*crhs6 - crhs11;
            rhs[1]=-DN(1,0)*crhs3 - DN(1,1)*crhs5 + crhs0*crhs13 - crhs10*crhs13 + crhs12*crhs6 + crhs12 - crhs14*crhs6 - crhs14;
            rhs[2]=-DN(2,0)*crhs3 - DN(2,1)*crhs5 + crhs0*crhs16 - crhs10*crhs16 + crhs15*crhs6 + crhs15 - crhs17*crhs6 - crhs17;


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
    const auto tau = rVariables.tau;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs1 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs2 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs3 =             DN(0,0)*clhs0 + DN(0,1)*clhs1 + DN(0,2)*clhs2;
const double clhs4 =             N[0]*clhs3;
const double clhs5 =             tau*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double clhs6 =             DN(1,0)*clhs0 + DN(1,1)*clhs1 + DN(1,2)*clhs2;
const double clhs7 =             N[0]*clhs6;
const double clhs8 =             DN(0,0)*k;
const double clhs9 =             DN(0,1)*k;
const double clhs10 =             DN(0,2)*k;
const double clhs11 =             clhs3*tau;
const double clhs12 =             DN(1,0)*clhs8 + DN(1,1)*clhs9 + DN(1,2)*clhs10 + clhs11*clhs6;
const double clhs13 =             DN(2,0)*clhs0 + DN(2,1)*clhs1 + DN(2,2)*clhs2;
const double clhs14 =             N[0]*clhs13;
const double clhs15 =             DN(2,0)*clhs8 + DN(2,1)*clhs9 + DN(2,2)*clhs10 + clhs11*clhs13;
const double clhs16 =             DN(3,0)*clhs0 + DN(3,1)*clhs1 + DN(3,2)*clhs2;
const double clhs17 =             N[0]*clhs16;
const double clhs18 =             DN(3,0)*clhs8 + DN(3,1)*clhs9 + DN(3,2)*clhs10 + clhs11*clhs16;
const double clhs19 =             N[1]*clhs3;
const double clhs20 =             N[1]*clhs6;
const double clhs21 =             N[1]*clhs13;
const double clhs22 =             DN(1,0)*k;
const double clhs23 =             DN(1,1)*k;
const double clhs24 =             DN(1,2)*k;
const double clhs25 =             clhs6*tau;
const double clhs26 =             DN(2,0)*clhs22 + DN(2,1)*clhs23 + DN(2,2)*clhs24 + clhs13*clhs25;
const double clhs27 =             N[1]*clhs16;
const double clhs28 =             DN(3,0)*clhs22 + DN(3,1)*clhs23 + DN(3,2)*clhs24 + clhs16*clhs25;
const double clhs29 =             N[2]*clhs3;
const double clhs30 =             N[2]*clhs6;
const double clhs31 =             N[2]*clhs13;
const double clhs32 =             N[2]*clhs16;
const double clhs33 =             DN(2,0)*DN(3,0)*k + DN(2,1)*DN(3,1)*k + DN(2,2)*DN(3,2)*k + clhs13*clhs16*tau;
const double clhs34 =             N[3]*clhs3;
const double clhs35 =             N[3]*clhs6;
const double clhs36 =             N[3]*clhs13;
const double clhs37 =             N[3]*clhs16;
            lhs(0,0)=pow(DN(0,0), 2)*k + pow(DN(0,1), 2)*k + pow(DN(0,2), 2)*k + pow(clhs3, 2)*tau + clhs4*clhs5 + clhs4;
            lhs(0,1)=clhs12 + clhs5*clhs7 + clhs7;
            lhs(0,2)=clhs14*clhs5 + clhs14 + clhs15;
            lhs(0,3)=clhs17*clhs5 + clhs17 + clhs18;
            lhs(1,0)=clhs12 + clhs19*clhs5 + clhs19;
            lhs(1,1)=pow(DN(1,0), 2)*k + pow(DN(1,1), 2)*k + pow(DN(1,2), 2)*k + clhs20*clhs5 + clhs20 + pow(clhs6, 2)*tau;
            lhs(1,2)=clhs21*clhs5 + clhs21 + clhs26;
            lhs(1,3)=clhs27*clhs5 + clhs27 + clhs28;
            lhs(2,0)=clhs15 + clhs29*clhs5 + clhs29;
            lhs(2,1)=clhs26 + clhs30*clhs5 + clhs30;
            lhs(2,2)=pow(DN(2,0), 2)*k + pow(DN(2,1), 2)*k + pow(DN(2,2), 2)*k + pow(clhs13, 2)*tau + clhs31*clhs5 + clhs31;
            lhs(2,3)=clhs32*clhs5 + clhs32 + clhs33;
            lhs(3,0)=clhs18 + clhs34*clhs5 + clhs34;
            lhs(3,1)=clhs28 + clhs35*clhs5 + clhs35;
            lhs(3,2)=clhs33 + clhs36*clhs5 + clhs36;
            lhs(3,3)=pow(DN(3,0), 2)*k + pow(DN(3,1), 2)*k + pow(DN(3,2), 2)*k + pow(clhs16, 2)*tau + clhs37*clhs5 + clhs37;


    const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs1 =             N[0]*crhs0;
const double crhs2 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs3 =             crhs2*k;
const double crhs4 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs5 =             crhs4*k;
const double crhs6 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs7 =             crhs6*k;
const double crhs8 =             tau*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double crhs9 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs10 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs11 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs12 =             tau*(DN(0,0)*crhs9 + DN(0,1)*crhs10 + DN(0,2)*crhs11);
const double crhs13 =             crhs10*crhs4 + crhs11*crhs6 + crhs2*crhs9;
const double crhs14 =             N[0]*crhs13;
const double crhs15 =             N[1]*crhs0;
const double crhs16 =             tau*(DN(1,0)*crhs9 + DN(1,1)*crhs10 + DN(1,2)*crhs11);
const double crhs17 =             N[1]*crhs13;
const double crhs18 =             N[2]*crhs0;
const double crhs19 =             tau*(DN(2,0)*crhs9 + DN(2,1)*crhs10 + DN(2,2)*crhs11);
const double crhs20 =             N[2]*crhs13;
const double crhs21 =             N[3]*crhs0;
const double crhs22 =             tau*(DN(3,0)*crhs9 + DN(3,1)*crhs10 + DN(3,2)*crhs11);
const double crhs23 =             N[3]*crhs13;
            rhs[0]=-DN(0,0)*crhs3 - DN(0,1)*crhs5 - DN(0,2)*crhs7 + crhs0*crhs12 + crhs1*crhs8 + crhs1 - crhs12*crhs13 - crhs14*crhs8 - crhs14;
            rhs[1]=-DN(1,0)*crhs3 - DN(1,1)*crhs5 - DN(1,2)*crhs7 + crhs0*crhs16 - crhs13*crhs16 + crhs15*crhs8 + crhs15 - crhs17*crhs8 - crhs17;
            rhs[2]=-DN(2,0)*crhs3 - DN(2,1)*crhs5 - DN(2,2)*crhs7 + crhs0*crhs19 - crhs13*crhs19 + crhs18*crhs8 + crhs18 - crhs20*crhs8 - crhs20;
            rhs[3]=-DN(3,0)*crhs3 - DN(3,1)*crhs5 - DN(3,2)*crhs7 + crhs0*crhs22 - crhs13*crhs22 + crhs21*crhs8 + crhs21 - crhs23*crhs8 - crhs23;


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