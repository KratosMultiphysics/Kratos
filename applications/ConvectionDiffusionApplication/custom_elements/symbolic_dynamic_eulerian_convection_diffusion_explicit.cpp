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

#include "symbolic_dynamic_eulerian_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDynamicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::SymbolicDynamicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::~SymbolicDynamicEulerianConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDynamicEulerianConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicDynamicEulerianConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
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
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rVariables.unknown_subscale = mUnknownSubScale(g);

        // Execute standard RHS-LHS build or OSS step
        if (rCurrentProcessInfo.GetValue(OSS_SWITCH) == 1)
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
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    // Resize and intialize output
    if (mUnknownSubScale.size() != 3) // three integration points
        mUnknownSubScale.resize(3, false);
    mUnknownSubScale = ZeroVector(3);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::FinalizeSolutionStep(
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
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
        // Compute tau
        this->CalculateTau(rVariables);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rVariables.unknown_subscale = mUnknownSubScale(g);

        // Update unknown belonging to subgrid scale space on gauss integration point g
        this->UpdateUnknownSubgridScaleGaussPoint(rVariables,g);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
int SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
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

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::ComputeGaussPointContribution(
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
    const auto qstau = rVariables.qstau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             clhs0*clhs1;
const double clhs3 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs4 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs5 =             DN(0,0)*clhs3 + DN(0,1)*clhs4;
const double clhs6 =             N[0]*clhs0*clhs1*tau;
const double clhs7 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double clhs8 =             N[0]*clhs7*tau;
const double clhs9 =             N[0]*clhs0*clhs1;
const double clhs10 =             DN(1,0)*clhs3 + DN(1,1)*clhs4;
const double clhs11 =             clhs5*tau;
const double clhs12 =             -N[1]*clhs9 + clhs10*clhs11;
const double clhs13 =             N[1]*clhs0*clhs1*tau;
const double clhs14 =             N[1]*clhs7*tau;
const double clhs15 =             DN(2,0)*clhs3 + DN(2,1)*clhs4;
const double clhs16 =             -N[2]*clhs9 + clhs11*clhs15;
const double clhs17 =             N[2]*clhs0*clhs1*tau;
const double clhs18 =             N[2]*clhs7*tau;
const double clhs19 =             -N[1]*N[2]*clhs0*clhs1 + clhs10*clhs15*tau;
            lhs(0,0)=-pow(N[0], 2)*clhs2 + pow(clhs5, 2)*tau + clhs5*clhs6 + clhs5*clhs8;
            lhs(0,1)=clhs12 + clhs13*clhs5 + clhs14*clhs5;
            lhs(0,2)=clhs16 + clhs17*clhs5 + clhs18*clhs5;
            lhs(1,0)=clhs10*clhs6 + clhs10*clhs8 + clhs12;
            lhs(1,1)=-pow(N[1], 2)*clhs2 + pow(clhs10, 2)*tau + clhs10*clhs13 + clhs10*clhs14;
            lhs(1,2)=clhs10*clhs17 + clhs10*clhs18 + clhs19;
            lhs(2,0)=clhs15*clhs6 + clhs15*clhs8 + clhs16;
            lhs(2,1)=clhs13*clhs15 + clhs14*clhs15 + clhs19;
            lhs(2,2)=-pow(N[2], 2)*clhs2 + pow(clhs15, 2)*tau + clhs15*clhs17 + clhs15*clhs18;


    const double crhs0 =             phi_subscale_gauss/qstau;
const double crhs1 =             N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2];
const double crhs2 =             1.0/delta_time;
const double crhs3 =             crhs2*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/RK_time_coefficient;
const double crhs4 =             crhs2*phi_subscale_gauss*tau;
const double crhs5 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs6 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs7 =             DN(0,0)*crhs5 + DN(0,1)*crhs6;
const double crhs8 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs9 =             crhs7*tau;
const double crhs10 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double crhs11 =             crhs5*(DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2]) + crhs6*(DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2]);
const double crhs12 =             DN(1,0)*crhs5 + DN(1,1)*crhs6;
const double crhs13 =             crhs12*tau;
const double crhs14 =             DN(2,0)*crhs5 + DN(2,1)*crhs6;
const double crhs15 =             crhs14*tau;
            rhs[0]=N[0]*crhs0 - N[0]*crhs1 + N[0]*crhs3 + crhs1*crhs9 - crhs10*crhs9 - crhs11*crhs9 - crhs3*crhs9 + crhs4*crhs7 + crhs8*crhs9;
            rhs[1]=N[1]*crhs0 - N[1]*crhs1 + N[1]*crhs3 + crhs1*crhs13 - crhs10*crhs13 - crhs11*crhs13 + crhs12*crhs4 - crhs13*crhs3 + crhs13*crhs8;
            rhs[2]=N[2]*crhs0 - N[2]*crhs1 + N[2]*crhs3 + crhs1*crhs15 - crhs10*crhs15 - crhs11*crhs15 + crhs14*crhs4 - crhs15*crhs3 + crhs15*crhs8;


    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeGaussPointContribution(
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
    const auto qstau = rVariables.qstau;
    const auto prj = rVariables.oss_projection;
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double clhs0 =             1.0/RK_time_coefficient;
const double clhs1 =             1.0/delta_time;
const double clhs2 =             clhs0*clhs1;
const double clhs3 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs4 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs5 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs6 =             DN(0,0)*clhs3 + DN(0,1)*clhs4 + DN(0,2)*clhs5;
const double clhs7 =             N[0]*clhs0*clhs1*tau;
const double clhs8 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double clhs9 =             N[0]*clhs8*tau;
const double clhs10 =             N[0]*clhs0*clhs1;
const double clhs11 =             DN(1,0)*clhs3 + DN(1,1)*clhs4 + DN(1,2)*clhs5;
const double clhs12 =             clhs6*tau;
const double clhs13 =             -N[1]*clhs10 + clhs11*clhs12;
const double clhs14 =             N[1]*clhs0*clhs1*tau;
const double clhs15 =             N[1]*clhs8*tau;
const double clhs16 =             DN(2,0)*clhs3 + DN(2,1)*clhs4 + DN(2,2)*clhs5;
const double clhs17 =             -N[2]*clhs10 + clhs12*clhs16;
const double clhs18 =             N[2]*clhs0*clhs1*tau;
const double clhs19 =             N[2]*clhs8*tau;
const double clhs20 =             DN(3,0)*clhs3 + DN(3,1)*clhs4 + DN(3,2)*clhs5;
const double clhs21 =             -N[3]*clhs10 + clhs12*clhs20;
const double clhs22 =             N[3]*clhs0*clhs1*tau;
const double clhs23 =             N[3]*clhs8*tau;
const double clhs24 =             N[1]*clhs0*clhs1;
const double clhs25 =             clhs11*tau;
const double clhs26 =             -N[2]*clhs24 + clhs16*clhs25;
const double clhs27 =             -N[3]*clhs24 + clhs20*clhs25;
const double clhs28 =             -N[2]*N[3]*clhs0*clhs1 + clhs16*clhs20*tau;
            lhs(0,0)=-pow(N[0], 2)*clhs2 + pow(clhs6, 2)*tau + clhs6*clhs7 + clhs6*clhs9;
            lhs(0,1)=clhs13 + clhs14*clhs6 + clhs15*clhs6;
            lhs(0,2)=clhs17 + clhs18*clhs6 + clhs19*clhs6;
            lhs(0,3)=clhs21 + clhs22*clhs6 + clhs23*clhs6;
            lhs(1,0)=clhs11*clhs7 + clhs11*clhs9 + clhs13;
            lhs(1,1)=-pow(N[1], 2)*clhs2 + pow(clhs11, 2)*tau + clhs11*clhs14 + clhs11*clhs15;
            lhs(1,2)=clhs11*clhs18 + clhs11*clhs19 + clhs26;
            lhs(1,3)=clhs11*clhs22 + clhs11*clhs23 + clhs27;
            lhs(2,0)=clhs16*clhs7 + clhs16*clhs9 + clhs17;
            lhs(2,1)=clhs14*clhs16 + clhs15*clhs16 + clhs26;
            lhs(2,2)=-pow(N[2], 2)*clhs2 + pow(clhs16, 2)*tau + clhs16*clhs18 + clhs16*clhs19;
            lhs(2,3)=clhs16*clhs22 + clhs16*clhs23 + clhs28;
            lhs(3,0)=clhs20*clhs7 + clhs20*clhs9 + clhs21;
            lhs(3,1)=clhs14*clhs20 + clhs15*clhs20 + clhs27;
            lhs(3,2)=clhs18*clhs20 + clhs19*clhs20 + clhs28;
            lhs(3,3)=-pow(N[3], 2)*clhs2 + pow(clhs20, 2)*tau + clhs20*clhs22 + clhs20*clhs23;


    const double crhs0 =             phi_subscale_gauss/qstau;
const double crhs1 =             N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3];
const double crhs2 =             1.0/delta_time;
const double crhs3 =             crhs2*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/RK_time_coefficient;
const double crhs4 =             crhs2*phi_subscale_gauss*tau;
const double crhs5 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs6 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs7 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs8 =             DN(0,0)*crhs5 + DN(0,1)*crhs6 + DN(0,2)*crhs7;
const double crhs9 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs10 =             crhs8*tau;
const double crhs11 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double crhs12 =             crhs5*(DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3]) + crhs6*(DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3]) + crhs7*(DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3]);
const double crhs13 =             DN(1,0)*crhs5 + DN(1,1)*crhs6 + DN(1,2)*crhs7;
const double crhs14 =             crhs13*tau;
const double crhs15 =             DN(2,0)*crhs5 + DN(2,1)*crhs6 + DN(2,2)*crhs7;
const double crhs16 =             crhs15*tau;
const double crhs17 =             DN(3,0)*crhs5 + DN(3,1)*crhs6 + DN(3,2)*crhs7;
const double crhs18 =             crhs17*tau;
            rhs[0]=N[0]*crhs0 - N[0]*crhs1 + N[0]*crhs3 + crhs1*crhs10 - crhs10*crhs11 - crhs10*crhs12 - crhs10*crhs3 + crhs10*crhs9 + crhs4*crhs8;
            rhs[1]=N[1]*crhs0 - N[1]*crhs1 + N[1]*crhs3 + crhs1*crhs14 - crhs11*crhs14 - crhs12*crhs14 + crhs13*crhs4 - crhs14*crhs3 + crhs14*crhs9;
            rhs[2]=N[2]*crhs0 - N[2]*crhs1 + N[2]*crhs3 + crhs1*crhs16 - crhs11*crhs16 - crhs12*crhs16 + crhs15*crhs4 - crhs16*crhs3 + crhs16*crhs9;
            rhs[3]=N[3]*crhs0 - N[3]*crhs1 + N[3]*crhs3 + crhs1*crhs18 - crhs11*crhs18 - crhs12*crhs18 + crhs17*crhs4 - crhs18*crhs3 + crhs18*crhs9;


    noalias(rLeftHandSideMatrix) += lhs * rVariables.weight;
    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::ComputeOSSGaussPointContribution(
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
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double crhs0 =             1.0/delta_time;
const double crhs1 =             crhs0*phi_subscale_gauss;
const double crhs2 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs3 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs6 =             crhs5*k;
const double crhs7 =             crhs0*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/RK_time_coefficient;
const double crhs8 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1));
const double crhs9 =             crhs3*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) + crhs5*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
            rhs[0]=DN(0,0)*crhs4 + DN(0,1)*crhs6 - N[0]*crhs1 - N[0]*crhs2 + N[0]*crhs7 + N[0]*crhs8 + N[0]*crhs9;
            rhs[1]=DN(1,0)*crhs4 + DN(1,1)*crhs6 - N[1]*crhs1 - N[1]*crhs2 + N[1]*crhs7 + N[1]*crhs8 + N[1]*crhs9;
            rhs[2]=DN(2,0)*crhs4 + DN(2,1)*crhs6 - N[2]*crhs1 - N[2]*crhs2 + N[2]*crhs7 + N[2]*crhs8 + N[2]*crhs9;


    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::ComputeOSSGaussPointContribution(
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
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    auto lhs = rVariables.lhs;
    auto rhs = rVariables.rhs;

    const double crhs0 =             1.0/delta_time;
const double crhs1 =             crhs0*phi_subscale_gauss;
const double crhs2 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs3 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs4 =             crhs3*k;
const double crhs5 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs6 =             crhs5*k;
const double crhs7 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs8 =             crhs7*k;
const double crhs9 =             crhs0*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/RK_time_coefficient;
const double crhs10 =             (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2));
const double crhs11 =             crhs3*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) + crhs5*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) + crhs7*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
            rhs[0]=DN(0,0)*crhs4 + DN(0,1)*crhs6 + DN(0,2)*crhs8 - N[0]*crhs1 + N[0]*crhs10 + N[0]*crhs11 - N[0]*crhs2 + N[0]*crhs9;
            rhs[1]=DN(1,0)*crhs4 + DN(1,1)*crhs6 + DN(1,2)*crhs8 - N[1]*crhs1 + N[1]*crhs10 + N[1]*crhs11 - N[1]*crhs2 + N[1]*crhs9;
            rhs[2]=DN(2,0)*crhs4 + DN(2,1)*crhs6 + DN(2,2)*crhs8 - N[2]*crhs1 + N[2]*crhs10 + N[2]*crhs11 - N[2]*crhs2 + N[2]*crhs9;
            rhs[3]=DN(3,0)*crhs4 + DN(3,1)*crhs6 + DN(3,2)*crhs8 - N[3]*crhs1 + N[3]*crhs10 + N[3]*crhs11 - N[3]*crhs2 + N[3]*crhs9;


    noalias(rRightHandSideVector) += rhs * rVariables.weight;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<2>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
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
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;
}

/***********************************************************************************/

template <>
void SymbolicDynamicEulerianConvectionDiffusionExplicit<3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementVariables& rVariables,
    unsigned int g)
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
    const auto phi_subscale_gauss = rVariables.unknown_subscale;
    const auto prj = rVariables.oss_projection;
    double phi_subscale_gauss_new = 0;

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/(RK_time_coefficient*delta_time); // mass term
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)); // convective term 1
    phi_subscale_gauss_new += - (DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateTau(
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
    inv_tau += 1.0/rVariables.delta_time;
    // Convection
    inv_tau += 2.0 * norm_velocity / h;
    inv_tau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
    // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
    inv_tau *= rVariables.density * rVariables.specific_heat;
    // Diffusion
    inv_tau += 4.0 * rVariables.diffusivity / (h*h);
    // Limiting
    inv_tau = std::max(inv_tau, 1e-2);
    rVariables.tau = (rVariables.density*rVariables.specific_heat) / inv_tau;
    // Estimate quasi-static tau
    double inv_qstau = 0;
    // Convection
    inv_qstau += 2.0 * norm_velocity / h;
    inv_qstau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
    // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
    inv_qstau *= rVariables.density * rVariables.specific_heat;
    // Diffusion
    inv_qstau += 4.0 * rVariables.diffusivity / (h*h);
    // Limiting
    inv_qstau = std::max(inv_qstau, 1e-2);
    rVariables.qstau = (rVariables.density*rVariables.specific_heat) / inv_qstau;
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicDynamicEulerianConvectionDiffusionExplicit<2>;
template class SymbolicDynamicEulerianConvectionDiffusionExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}
