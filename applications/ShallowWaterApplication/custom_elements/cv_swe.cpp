//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "cv_swe.h"

namespace Kratos
{

template< size_t TNumNodes, ElementFramework TFramework >
int CV_SWE<TNumNodes, TFramework>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(MOMENTUM)
    this->CheckVariableKey();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( size_t i = 0; i < TNumNodes; i++ )
    {
        Node<3>& node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, node)
        this->CheckVariableInNodalData(node);

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return ierr;

    KRATOS_CATCH("")
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const size_t element_size = TNumNodes*3;
    if(rResult.size() != element_size)
        rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

    const GeometryType& rGeom = this->GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const size_t element_size = TNumNodes*3;
    if(rElementalDofList.size() != element_size)
        rElementalDofList.resize(element_size);

    const GeometryType& rGeom = this->GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // Resize of the Left and Right Hand side
    constexpr size_t element_size = TNumNodes*3;
    if(rLeftHandSideMatrix.size1() != element_size)
        rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != element_size)
        rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!

    const GeometryType& Geom = this->GetGeometry();

    ConservativeElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    const BoundedMatrix<double,TNumNodes, TNumNodes> NContainer = Geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 ); // In this case, number of Gauss points and number of nodes coincides

    BoundedMatrix<double,TNumNodes, 2> DN_DX;  // Shape function gradients are constant since we are using linear functions
    array_1d<double,TNumNodes> N;
    double Area;
    this-> CalculateGeometry(DN_DX,Area);

    this->GetNodalValues(variables);
    this->CalculateElementValues(DN_DX, variables);

    rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);
    rRightHandSideVector = ZeroVector(element_size);

    for (size_t i_gauss = 0; i_gauss < TNumNodes; ++i_gauss)
    {
        noalias(N) = row(NContainer, i_gauss);

        this->BuildMassMatrices(N, variables);
        this->BuildGradientMatrices(N, DN_DX, variables);
        this->BuildDiffusivityMatrices(DN_DX, variables);
        this->BuildConvectionMatrices(N, DN_DX, variables);

        this->AddInertiaTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddConvectiveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddWaveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddFrictionTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddStabilizationTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddSourceTerms(rRightHandSideVector, variables);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, variables.unknown);

    rRightHandSideVector *= Area * variables.lumping_factor;
    rLeftHandSideMatrix  *= Area * variables.lumping_factor;

    double residual = norm_1 (rRightHandSideVector);
    this->SetValue(RESIDUAL_NORM, residual);
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::GetNodalValues(ConservativeElementVariables& rVariables)
{
    const GeometryType& rGeom = this->GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.depth[counter] = 0;
        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X, 1);
        rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X);
        rVariables.nodal_velocity[counter] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
        counter++;

        rVariables.depth[counter] = 0;
        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y, 1);
        rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y);
        rVariables.nodal_velocity[counter] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
        counter++;

        rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / rVariables.height_units;
        rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(HEIGHT, 1);
        rVariables.proj_unk[counter] = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1);
        rVariables.nodal_velocity[counter] = 0;
        counter++;
    }
    if (TFramework == PFEM2)
    {
        rVariables.prev_unk = rVariables.proj_unk;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::CalculateElementValues(
    const BoundedMatrix<double,TNumNodes, 2>& rDN_DX,
    ConservativeElementVariables& rVariables)
{
    // Initialize outputs
    rVariables.momentum = ZeroVector(2);
    rVariables.height = 0;
    rVariables.momentum_div = 0;
    rVariables.height_grad = ZeroVector(2);
    rVariables.projected_momentum = ZeroVector(2);
    rVariables.velocity = ZeroVector(2);

    // integrate over the element
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.momentum[0] += rVariables.unknown[  + 3*i];
        rVariables.momentum[1] += rVariables.unknown[1 + 3*i];
        rVariables.height += rVariables.unknown[2 + 3*i];
        rVariables.momentum_div += rDN_DX(i,0) * rVariables.unknown[  + 3*i] + rDN_DX(i,1) * rVariables.unknown[1 + 3*i];
        rVariables.height_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
        rVariables.height_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
        rVariables.projected_momentum[0] += rVariables.proj_unk[  + 3*i];
        rVariables.projected_momentum[1] += rVariables.proj_unk[1 + 3*i];
        rVariables.velocity[0] += rVariables.nodal_velocity[  + 3*i];
        rVariables.velocity[1] += rVariables.nodal_velocity[1 + 3*i];
    }

    rVariables.momentum *= rVariables.lumping_factor;
    rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
    rVariables.height_grad *= rVariables.height_units;
    rVariables.projected_momentum *= rVariables.lumping_factor;
    rVariables.velocity *= rVariables.lumping_factor;

    rVariables.wave_vel_2 = rVariables.gravity * std::abs(rVariables.height);
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::ComputeStabilizationParameters(
    const ConservativeElementVariables& rVariables,
    double& rTauU,
    double& rTauH,
    double& rKappaU,
    double& rKappaH)
{
    // Get element values
    const double elem_size = this->GetGeometry().Length();
    const double height = std::abs(rVariables.height) + rVariables.epsilon;
    const double CTau = rVariables.dyn_tau;  // 0.005 ~ 0.002

    // Wave mixed form stabilization
    rTauU = CTau * elem_size * std::sqrt(rVariables.gravity / height);
    rTauH = CTau * elem_size * std::sqrt(height / rVariables.gravity);

    // Discontinuity capturing
    double vel_grad_norm = std::abs(rVariables.momentum_div);
    double height_grad_norm = norm_2(rVariables.height_grad);
    rTauU += 0.5 * 0.01 * elem_size * vel_grad_norm;
    rTauH += 0.5 * 0.01 * elem_size * height_grad_norm;

    // Convective stabilization
    if (TFramework == Eulerian)
    {
        const double vel_modulus = norm_2(rVariables.velocity) + rVariables.epsilon;
        rKappaU = CTau * elem_size / vel_modulus;
        rKappaH = 0.0;
    }
    else
    {
        rKappaU = 0.0;
        rKappaH = 0.0;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::BuildConvectionMatrices(
    array_1d<double, TNumNodes>& rN,
    BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
    ConservativeElementVariables& rVariables)
{
    constexpr size_t local_size = TNumNodes * 3;
    noalias(rVariables.Convection) = ZeroMatrix(local_size,local_size);
    noalias(rVariables.ScalarConvectionStabilization) = ZeroMatrix(local_size,local_size);
    noalias(rVariables.VectorConvectionStabilization) = ZeroMatrix(local_size,local_size);
    noalias(rVariables.FrictionStabilization) = ZeroMatrix(local_size,local_size);

    if (TFramework == Eulerian)
    {
        // Auxiliary definitions
        BoundedMatrix<double,2,local_size> N_vel        = ZeroMatrix(2,local_size);  // Shape functions matrix (for velocity unknown)
        BoundedMatrix<double,2,local_size> Grad_vel_1   = ZeroMatrix(1,local_size);  // Shape functions gradients vector (for velocity unknown)
        BoundedMatrix<double,2,local_size> Grad_vel_2   = ZeroMatrix(1,local_size);  // Shape functions gradients vector (for velocity unknown)
        BoundedMatrix<double,2,local_size> DN_DX_height = ZeroMatrix(2,local_size);  // Shape functions gradients matrix (for height unknown)

        // Build the shape and derivatives functions at the Gauss point
        for(size_t node = 0; node < TNumNodes; ++node)
        {
            // Velocity gradient
            Grad_vel_1(0,   node*3) = rDN_DX(node,0);
            Grad_vel_1(1, 1+node*3) = rDN_DX(node,0);
            Grad_vel_2(0,   node*3) = rDN_DX(node,1);
            Grad_vel_2(1, 1+node*3) = rDN_DX(node,1);
            // Velocity shape functions
            N_vel(0,   node*3) = rN[node];
            N_vel(1, 1+node*3) = rN[node];
            // Height gradient
            DN_DX_height(0, 2+node*3) = rDN_DX(node,0);
            DN_DX_height(1, 2+node*3) = rDN_DX(node,1);
        }

        BoundedMatrix<double,2,local_size> vector_convection_operator = rVariables.velocity[0] * Grad_vel_1 + rVariables.velocity[1] * Grad_vel_2;

        rVariables.Convection += prod(trans(N_vel), vector_convection_operator);   // w * u * grad_q
        double vel2 = rVariables.velocity[0] * rVariables.velocity[0] + rVariables.velocity[1] * rVariables.velocity[1];
        rVariables.Convection += vel2 * prod(trans(N_vel), DN_DX_height); // w * u * u * grad_h
        rVariables.VectorConvectionStabilization += prod(trans(vector_convection_operator), vector_convection_operator);// grad_w * u * grad_u
        double vel1 = std::sqrt(vel2);
        rVariables.VectorConvectionStabilization += vel1 * prod(trans(vector_convection_operator), DN_DX_height);     // grad_w * u * u * grad_h

        // Mass balance stabilization
        rVariables.ScalarConvectionStabilization += vel2 * prod(trans(DN_DX_height), DN_DX_height); // grad_n * u * u * grad_h
        rVariables.ScalarConvectionStabilization += prod(trans(DN_DX_height), vector_convection_operator); // grad_n * u * grad_q
        rVariables.FrictionStabilization += prod(trans(DN_DX_height), N_vel);  // grad_n * u
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::AddWaveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ConservativeElementVariables& rVariables)
{
    rLeftHandSideMatrix += rVariables.VectorDiv;
    rLeftHandSideMatrix += rVariables.wave_vel_2 * rVariables.ScalarGrad;
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::AddFrictionTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ConservativeElementVariables& rVariables)
{
    const double abs_mom = norm_2(rVariables.projected_momentum);
    const double height73 = std::pow(std::abs(rVariables.height), 2.3333333333333) + rVariables.epsilon;
    rLeftHandSideMatrix += rVariables.gravity * rVariables.manning2 * abs_mom / height73 * rVariables.MassMatrixVector;
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::AddStabilizationTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ConservativeElementVariables& rVariables)
{
    double tau_u;
    double tau_h;
    double kappa_u;
    double kappa_h;
    this->ComputeStabilizationParameters(rVariables, tau_u, tau_h, kappa_u, kappa_h);
    // Momentum stabilization
    rLeftHandSideMatrix += tau_u * rVariables.VectorDiff;
    rLeftHandSideMatrix += kappa_u * rVariables.VectorConvectionStabilization;
    // Mass stabilization
    rLeftHandSideMatrix += tau_h * rVariables.wave_vel_2 * rVariables.ScalarDiff;
    rRightHandSideVector += tau_h * rVariables.wave_vel_2 * prod(rVariables.ScalarDiff, rVariables.depth);
    rLeftHandSideMatrix += tau_h * rVariables.ScalarConvectionStabilization;
    const double abs_mom = norm_2(rVariables.projected_momentum);
    const double height73 = std::pow(std::abs(rVariables.height), 2.3333333333333) + rVariables.epsilon;
    rLeftHandSideMatrix += tau_h * rVariables.gravity * rVariables.manning2 * abs_mom / height73 * rVariables.FrictionStabilization;
}


template< size_t TNumNodes, ElementFramework TFramework >
void CV_SWE<TNumNodes, TFramework>::AddSourceTerms(
    VectorType& rRightHandSideVector,
    ConservativeElementVariables& rVariables)
{
    rRightHandSideVector += rVariables.wave_vel_2 * prod(rVariables.ScalarGrad, rVariables.depth);
    rRightHandSideVector += prod(rVariables.MassMatrixScalar, rVariables.rain);
}


template class CV_SWE<3, Eulerian>;
template class CV_SWE<4, Eulerian>;
template class CV_SWE<3, PFEM2>;
template class CV_SWE<4, PFEM2>;

} // namespace Kratos
