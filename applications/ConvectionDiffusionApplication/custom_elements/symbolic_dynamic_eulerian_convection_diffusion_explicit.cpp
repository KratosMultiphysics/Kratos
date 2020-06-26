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

    const double clhs0 =             2*k;
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             1.0/RK_time_coefficient;
const double clhs3 =             1.0/delta_time;
const double clhs4 =             clhs2*clhs3;
const double clhs5 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1);
const double clhs6 =             2*N[0];
const double clhs7 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double clhs8 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double clhs9 =             DN(0,0)*clhs7 + DN(0,1)*clhs8;
const double clhs10 =             N[0]*clhs2*clhs3*tau;
const double clhs11 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double clhs12 =             N[0]*clhs11*tau;
const double clhs13 =             2*DN(0,0)*k;
const double clhs14 =             2*DN(0,1)*k;
const double clhs15 =             N[0]*clhs2*clhs3;
const double clhs16 =             2*N[0]*clhs11;
const double clhs17 =             DN(1,0)*clhs7 + DN(1,1)*clhs8;
const double clhs18 =             clhs9*tau;
const double clhs19 =             DN(1,0)*clhs13 + DN(1,1)*clhs14 + N[1]*clhs15 + N[1]*clhs16 + clhs17*clhs18;
const double clhs20 =             N[1]*clhs2*clhs3*tau;
const double clhs21 =             N[1]*clhs11*tau;
const double clhs22 =             DN(2,0)*clhs7 + DN(2,1)*clhs8;
const double clhs23 =             DN(2,0)*clhs13 + DN(2,1)*clhs14 + N[2]*clhs15 + N[2]*clhs16 + clhs18*clhs22;
const double clhs24 =             N[2]*clhs2*clhs3*tau;
const double clhs25 =             N[2]*clhs11*tau;
const double clhs26 =             2*N[1];
const double clhs27 =             pow(N[1], 2);
const double clhs28 =             N[1]*N[2];
const double clhs29 =             DN(1,0)*DN(2,0)*clhs0 + DN(1,1)*DN(2,1)*clhs0 + clhs17*clhs22*tau + clhs28*clhs4 + clhs28*clhs5;
const double clhs30 =             2*N[2];
const double clhs31 =             pow(N[2], 2);
            lhs(0,0)=pow(DN(0,0), 2)*clhs0 + pow(DN(0,1), 2)*clhs0 + clhs1*clhs4 + clhs1*clhs5 + clhs10*clhs9 + clhs12*clhs9 + clhs6*clhs9 + pow(clhs9, 2)*tau;
            lhs(0,1)=clhs17*clhs6 + clhs19 + clhs20*clhs9 + clhs21*clhs9;
            lhs(0,2)=clhs22*clhs6 + clhs23 + clhs24*clhs9 + clhs25*clhs9;
            lhs(1,0)=clhs10*clhs17 + clhs12*clhs17 + clhs19 + clhs26*clhs9;
            lhs(1,1)=pow(DN(1,0), 2)*clhs0 + pow(DN(1,1), 2)*clhs0 + pow(clhs17, 2)*tau + clhs17*clhs20 + clhs17*clhs21 + clhs17*clhs26 + clhs27*clhs4 + clhs27*clhs5;
            lhs(1,2)=clhs17*clhs24 + clhs17*clhs25 + clhs22*clhs26 + clhs29;
            lhs(2,0)=clhs10*clhs22 + clhs12*clhs22 + clhs23 + clhs30*clhs9;
            lhs(2,1)=clhs17*clhs30 + clhs20*clhs22 + clhs21*clhs22 + clhs29;
            lhs(2,2)=pow(DN(2,0), 2)*clhs0 + pow(DN(2,1), 2)*clhs0 + pow(clhs22, 2)*tau + clhs22*clhs24 + clhs22*clhs25 + clhs22*clhs30 + clhs31*clhs4 + clhs31*clhs5;


    const double crhs0 =             phi_subscale_gauss/qstau;
const double crhs1 =             N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2];
const double crhs2 =             2*N[0]*f[0] + 2*N[1]*f[1] + 2*N[2]*f[2];
const double crhs3 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs4 =             2*crhs3*k;
const double crhs5 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
const double crhs6 =             2*crhs5*k;
const double crhs7 =             1.0/delta_time;
const double crhs8 =             crhs7*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/RK_time_coefficient;
const double crhs9 =             crhs7*phi_subscale_gauss*tau;
const double crhs10 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs11 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs12 =             DN(0,0)*crhs10 + DN(0,1)*crhs11;
const double crhs13 =             N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2];
const double crhs14 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1);
const double crhs15 =             2*crhs13*crhs14;
const double crhs16 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs17 =             crhs12*tau;
const double crhs18 =             2*crhs10*crhs3 + 2*crhs11*crhs5;
const double crhs19 =             crhs13*crhs14;
const double crhs20 =             crhs10*crhs3 + crhs11*crhs5;
const double crhs21 =             DN(1,0)*crhs10 + DN(1,1)*crhs11;
const double crhs22 =             crhs21*tau;
const double crhs23 =             DN(2,0)*crhs10 + DN(2,1)*crhs11;
const double crhs24 =             crhs23*tau;
            rhs[0]=-DN(0,0)*crhs4 - DN(0,1)*crhs6 - N[0]*crhs0 + N[0]*crhs1 - N[0]*crhs15 - N[0]*crhs18 + N[0]*crhs2 - N[0]*crhs8 + crhs1*crhs17 + crhs12*crhs9 + crhs16*crhs17 - crhs17*crhs19 - crhs17*crhs20 - crhs17*crhs8;
            rhs[1]=-DN(1,0)*crhs4 - DN(1,1)*crhs6 - N[1]*crhs0 + N[1]*crhs1 - N[1]*crhs15 - N[1]*crhs18 + N[1]*crhs2 - N[1]*crhs8 + crhs1*crhs22 + crhs16*crhs22 - crhs19*crhs22 - crhs20*crhs22 + crhs21*crhs9 - crhs22*crhs8;
            rhs[2]=-DN(2,0)*crhs4 - DN(2,1)*crhs6 - N[2]*crhs0 + N[2]*crhs1 - N[2]*crhs15 - N[2]*crhs18 + N[2]*crhs2 - N[2]*crhs8 + crhs1*crhs24 + crhs16*crhs24 - crhs19*crhs24 - crhs20*crhs24 + crhs23*crhs9 - crhs24*crhs8;


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

    const double clhs0 =             2*k;
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             1.0/RK_time_coefficient;
const double clhs3 =             1.0/delta_time;
const double clhs4 =             clhs2*clhs3;
const double clhs5 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(0,2)*v(0,2) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(1,2)*v(1,2) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1) + 2*DN(2,2)*v(2,2) + 2*DN(3,0)*v(3,0) + 2*DN(3,1)*v(3,1) + 2*DN(3,2)*v(3,2);
const double clhs6 =             2*N[0];
const double clhs7 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double clhs8 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double clhs9 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double clhs10 =             DN(0,0)*clhs7 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
const double clhs11 =             N[0]*clhs2*clhs3*tau;
const double clhs12 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double clhs13 =             N[0]*clhs12*tau;
const double clhs14 =             2*DN(0,0)*k;
const double clhs15 =             2*DN(0,1)*k;
const double clhs16 =             2*DN(0,2)*k;
const double clhs17 =             N[0]*clhs2*clhs3;
const double clhs18 =             2*N[0]*clhs12;
const double clhs19 =             DN(1,0)*clhs7 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
const double clhs20 =             clhs10*tau;
const double clhs21 =             DN(1,0)*clhs14 + DN(1,1)*clhs15 + DN(1,2)*clhs16 + N[1]*clhs17 + N[1]*clhs18 + clhs19*clhs20;
const double clhs22 =             N[1]*clhs2*clhs3*tau;
const double clhs23 =             N[1]*clhs12*tau;
const double clhs24 =             DN(2,0)*clhs7 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
const double clhs25 =             DN(2,0)*clhs14 + DN(2,1)*clhs15 + DN(2,2)*clhs16 + N[2]*clhs17 + N[2]*clhs18 + clhs20*clhs24;
const double clhs26 =             N[2]*clhs2*clhs3*tau;
const double clhs27 =             N[2]*clhs12*tau;
const double clhs28 =             DN(3,0)*clhs7 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
const double clhs29 =             DN(3,0)*clhs14 + DN(3,1)*clhs15 + DN(3,2)*clhs16 + N[3]*clhs17 + N[3]*clhs18 + clhs20*clhs28;
const double clhs30 =             N[3]*clhs2*clhs3*tau;
const double clhs31 =             N[3]*clhs12*tau;
const double clhs32 =             2*N[1];
const double clhs33 =             pow(N[1], 2);
const double clhs34 =             2*DN(1,0)*k;
const double clhs35 =             2*DN(1,1)*k;
const double clhs36 =             2*DN(1,2)*k;
const double clhs37 =             N[1]*clhs2*clhs3;
const double clhs38 =             2*N[1]*clhs12;
const double clhs39 =             clhs19*tau;
const double clhs40 =             DN(2,0)*clhs34 + DN(2,1)*clhs35 + DN(2,2)*clhs36 + N[2]*clhs37 + N[2]*clhs38 + clhs24*clhs39;
const double clhs41 =             DN(3,0)*clhs34 + DN(3,1)*clhs35 + DN(3,2)*clhs36 + N[3]*clhs37 + N[3]*clhs38 + clhs28*clhs39;
const double clhs42 =             2*N[2];
const double clhs43 =             pow(N[2], 2);
const double clhs44 =             N[2]*N[3];
const double clhs45 =             DN(2,0)*DN(3,0)*clhs0 + DN(2,1)*DN(3,1)*clhs0 + DN(2,2)*DN(3,2)*clhs0 + clhs24*clhs28*tau + clhs4*clhs44 + clhs44*clhs5;
const double clhs46 =             2*N[3];
const double clhs47 =             pow(N[3], 2);
            lhs(0,0)=pow(DN(0,0), 2)*clhs0 + pow(DN(0,1), 2)*clhs0 + pow(DN(0,2), 2)*clhs0 + clhs1*clhs4 + clhs1*clhs5 + pow(clhs10, 2)*tau + clhs10*clhs11 + clhs10*clhs13 + clhs10*clhs6;
            lhs(0,1)=clhs10*clhs22 + clhs10*clhs23 + clhs19*clhs6 + clhs21;
            lhs(0,2)=clhs10*clhs26 + clhs10*clhs27 + clhs24*clhs6 + clhs25;
            lhs(0,3)=clhs10*clhs30 + clhs10*clhs31 + clhs28*clhs6 + clhs29;
            lhs(1,0)=clhs10*clhs32 + clhs11*clhs19 + clhs13*clhs19 + clhs21;
            lhs(1,1)=pow(DN(1,0), 2)*clhs0 + pow(DN(1,1), 2)*clhs0 + pow(DN(1,2), 2)*clhs0 + pow(clhs19, 2)*tau + clhs19*clhs22 + clhs19*clhs23 + clhs19*clhs32 + clhs33*clhs4 + clhs33*clhs5;
            lhs(1,2)=clhs19*clhs26 + clhs19*clhs27 + clhs24*clhs32 + clhs40;
            lhs(1,3)=clhs19*clhs30 + clhs19*clhs31 + clhs28*clhs32 + clhs41;
            lhs(2,0)=clhs10*clhs42 + clhs11*clhs24 + clhs13*clhs24 + clhs25;
            lhs(2,1)=clhs19*clhs42 + clhs22*clhs24 + clhs23*clhs24 + clhs40;
            lhs(2,2)=pow(DN(2,0), 2)*clhs0 + pow(DN(2,1), 2)*clhs0 + pow(DN(2,2), 2)*clhs0 + pow(clhs24, 2)*tau + clhs24*clhs26 + clhs24*clhs27 + clhs24*clhs42 + clhs4*clhs43 + clhs43*clhs5;
            lhs(2,3)=clhs24*clhs30 + clhs24*clhs31 + clhs28*clhs42 + clhs45;
            lhs(3,0)=clhs10*clhs46 + clhs11*clhs28 + clhs13*clhs28 + clhs29;
            lhs(3,1)=clhs19*clhs46 + clhs22*clhs28 + clhs23*clhs28 + clhs41;
            lhs(3,2)=clhs24*clhs46 + clhs26*clhs28 + clhs27*clhs28 + clhs45;
            lhs(3,3)=pow(DN(3,0), 2)*clhs0 + pow(DN(3,1), 2)*clhs0 + pow(DN(3,2), 2)*clhs0 + pow(clhs28, 2)*tau + clhs28*clhs30 + clhs28*clhs31 + clhs28*clhs46 + clhs4*clhs47 + clhs47*clhs5;


    const double crhs0 =             phi_subscale_gauss/qstau;
const double crhs1 =             N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3];
const double crhs2 =             2*N[0]*f[0] + 2*N[1]*f[1] + 2*N[2]*f[2] + 2*N[3]*f[3];
const double crhs3 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs4 =             2*crhs3*k;
const double crhs5 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs6 =             2*crhs5*k;
const double crhs7 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
const double crhs8 =             2*crhs7*k;
const double crhs9 =             1.0/delta_time;
const double crhs10 =             crhs9*(N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/RK_time_coefficient;
const double crhs11 =             crhs9*phi_subscale_gauss*tau;
const double crhs12 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs13 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs14 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs15 =             DN(0,0)*crhs12 + DN(0,1)*crhs13 + DN(0,2)*crhs14;
const double crhs16 =             N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3];
const double crhs17 =             DN(0,0)*v(0,0) + DN(0,1)*v(0,1) + DN(0,2)*v(0,2) + DN(1,0)*v(1,0) + DN(1,1)*v(1,1) + DN(1,2)*v(1,2) + DN(2,0)*v(2,0) + DN(2,1)*v(2,1) + DN(2,2)*v(2,2) + DN(3,0)*v(3,0) + DN(3,1)*v(3,1) + DN(3,2)*v(3,2);
const double crhs18 =             2*crhs16*crhs17;
const double crhs19 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs20 =             crhs15*tau;
const double crhs21 =             2*crhs12*crhs3 + 2*crhs13*crhs5 + 2*crhs14*crhs7;
const double crhs22 =             crhs16*crhs17;
const double crhs23 =             crhs12*crhs3 + crhs13*crhs5 + crhs14*crhs7;
const double crhs24 =             DN(1,0)*crhs12 + DN(1,1)*crhs13 + DN(1,2)*crhs14;
const double crhs25 =             crhs24*tau;
const double crhs26 =             DN(2,0)*crhs12 + DN(2,1)*crhs13 + DN(2,2)*crhs14;
const double crhs27 =             crhs26*tau;
const double crhs28 =             DN(3,0)*crhs12 + DN(3,1)*crhs13 + DN(3,2)*crhs14;
const double crhs29 =             crhs28*tau;
            rhs[0]=-DN(0,0)*crhs4 - DN(0,1)*crhs6 - DN(0,2)*crhs8 - N[0]*crhs0 + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs18 + N[0]*crhs2 - N[0]*crhs21 + crhs1*crhs20 - crhs10*crhs20 + crhs11*crhs15 + crhs19*crhs20 - crhs20*crhs22 - crhs20*crhs23;
            rhs[1]=-DN(1,0)*crhs4 - DN(1,1)*crhs6 - DN(1,2)*crhs8 - N[1]*crhs0 + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs18 + N[1]*crhs2 - N[1]*crhs21 + crhs1*crhs25 - crhs10*crhs25 + crhs11*crhs24 + crhs19*crhs25 - crhs22*crhs25 - crhs23*crhs25;
            rhs[2]=-DN(2,0)*crhs4 - DN(2,1)*crhs6 - DN(2,2)*crhs8 - N[2]*crhs0 + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs18 + N[2]*crhs2 - N[2]*crhs21 + crhs1*crhs27 - crhs10*crhs27 + crhs11*crhs26 + crhs19*crhs27 - crhs22*crhs27 - crhs23*crhs27;
            rhs[3]=-DN(3,0)*crhs4 - DN(3,1)*crhs6 - DN(3,2)*crhs8 - N[3]*crhs0 + N[3]*crhs1 - N[3]*crhs10 - N[3]*crhs18 + N[3]*crhs2 - N[3]*crhs21 + crhs1*crhs29 - crhs10*crhs29 + crhs11*crhs28 + crhs19*crhs29 - crhs22*crhs29 - crhs23*crhs29;


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
