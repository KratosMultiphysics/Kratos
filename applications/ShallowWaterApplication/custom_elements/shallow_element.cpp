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
#include "custom_elements/shallow_element.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

//----------------------------------------------------------------------

int ShallowElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_SCALAR1)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_VECTOR1)
    KRATOS_CHECK_VARIABLE_KEY(BATHYMETRY)
    KRATOS_CHECK_VARIABLE_KEY(RAIN)
    KRATOS_CHECK_VARIABLE_KEY(MANNING)
    KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME)
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < msNodes; i++ )
    {
        Node<3>& node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, node)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return ierr;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rResult.size() != msElemSize)
        rResult.resize(msElemSize,false); // False says not to preserve existing storage!!

    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < msNodes; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rElementalDofList.size() != msElemSize)
        rElementalDofList.resize(msElemSize);
    
    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < msNodes; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }
    
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize of the Left and Right Hand side
    if(rLeftHandSideMatrix.size1() != msElemSize)
        rLeftHandSideMatrix.resize(msElemSize,msElemSize,false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != msElemSize)
        rRightHandSideVector.resize(msElemSize,false);          // False says not to preserve existing storage!!

    // Struct to pass around the data
    ElementData data;
    this->InitializeElement(data, rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point 
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

    double length = GetGeometry().Length();

    // Mass matrix
    BoundedMatrix<double,msElemSize,msElemSize> vel_mass_matrix = ZeroMatrix(msElemSize,msElemSize);
    BoundedMatrix<double,msElemSize,msElemSize> h_mass_matrix = ZeroMatrix(msElemSize,msElemSize);
    ComputeMassMatrices(data, vel_mass_matrix, h_mass_matrix);

    // main loop (one Gauss point)

    // Build shape functions and derivatives at the Gauss point
    BoundedMatrix<double,2,msElemSize> N_vel = ZeroMatrix(2,msElemSize);        // Shape functions matrix (for velocity)
    array_1d<double,msElemSize> N_height = ZeroVector(msElemSize);              // Shape functions vector (for height)
    array_1d<double,msElemSize> DN_DX_vel = ZeroVector(msElemSize);             // Gradients vector (for velocity)
    BoundedMatrix<double,2,msElemSize> DN_DX_height = ZeroMatrix(2,msElemSize); // Gradients matrix (for height)
    for(unsigned int node = 0; node < msNodes; node++)
    {
        // Velocity divergence
        DN_DX_vel[  node*3] = DN_DX(node,0);
        DN_DX_vel[1+node*3] = DN_DX(node,1);
        // Height gradient
        DN_DX_height(0, 2+node*3) = DN_DX(node,0);
        DN_DX_height(1, 2+node*3) = DN_DX(node,1);
        // Height shape funtions
        N_height[2+node*3] = N[node];
        // Velocity shape functions
        N_vel(0,   node*3) = N[node];
        N_vel(1, 1+node*3) = N[node];
    }

    array_1d<double,2> vel;
    double height;
    ComputeElementValues(data, vel, height);

    // Auxiliary values to compute friction and stabilization terms
    const double epsilon = 0.005;
    const double abs_vel = norm_2(vel);
    const double aux_height = std::abs(height) + epsilon;
    const double height4_3 = std::pow(aux_height, 1.33333333);
    const double c = std::sqrt(data.gravity * aux_height);

    // Build LHS
    // Wave equation terms
    BoundedMatrix<double,msElemSize,msElemSize> vel_wave = prod(trans(N_vel),DN_DX_height);
    noalias(rLeftHandSideMatrix)  = data.gravity * vel_wave;                  // Add <w*g*grad(h)> to Momentum Eq
    noalias(rLeftHandSideMatrix) += height * outer_prod(N_height, DN_DX_vel); // Add <q*h*div(u)> to Mass Eq

    // Inertia terms
    noalias(rLeftHandSideMatrix) += data.dt_inv * vel_mass_matrix; // Velocity lumped mass matrix
    noalias(rLeftHandSideMatrix) += data.dt_inv * h_mass_matrix;   // Height lumped mass matrix

    // Friction term
    noalias(rLeftHandSideMatrix) += data.gravity * data.manning2 * abs_vel / height4_3 * vel_mass_matrix;

    // Build RHS
    // Source terms (bathymetry contribution)
    noalias(rRightHandSideVector)  = -data.gravity * prod(vel_wave, data.depth);

    // Inertia terms
    noalias(rRightHandSideVector) += data.dt_inv * prod(vel_mass_matrix, data.proj_unk);
    noalias(rRightHandSideVector) += data.dt_inv * prod(h_mass_matrix, data.proj_unk);

    // Computing stabilization terms
    double art = length * data.c_tau / c;
    double k_dc = 0.5 * length * data.c_tau;
    const array_1d<double,2> height_grad = prod(DN_DX_height, data.unknown);
    const array_1d<double,msElemSize> residual = rRightHandSideVector - prod(rLeftHandSideMatrix, data.unknown);
    const double grad_norm = norm_2(height_grad);
    const double h_residual = data.lumping_factor * (residual[2] + residual[5] + residual[8]);
    if (grad_norm > 1e-3)
        k_dc *= std::abs(h_residual / grad_norm);
        // k_dc *= grad_norm; // The easiest and heaviest way. Not good for complex free surface problems
        // k_dc *= std::abs(inner_prod(vel, height_grad)) / grad_norm;
    else
        k_dc *= 0;

    this->SetValue(RESIDUAL_NORM, h_residual);
    this->SetValue(PR_ART_VISC, k_dc);
    this->SetValue(VEL_ART_VISC, art);

    // Mass balance LHS stabilization terms
    BoundedMatrix<double,msElemSize,msElemSize> diff_h = prod(trans(DN_DX_height), DN_DX_height);
    noalias(rLeftHandSideMatrix) += k_dc * diff_h; // Second order FIC shock capturing
    noalias(rRightHandSideVector) -= k_dc * prod(diff_h, data.depth); // Substracting the bottom diffusion

    // Momentum balance stabilization terms
    BoundedMatrix<double,msElemSize,msElemSize> diff_v = outer_prod(DN_DX_vel, DN_DX_vel);
    noalias(rLeftHandSideMatrix) += art * diff_v; // Second order FIC

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, data.unknown);

    rRightHandSideVector *= area;
    rLeftHandSideMatrix *= area;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowElement: CalculateLeftHandSide not implemented" << std::endl;
}

//----------------------------------------------------------------------

void ShallowElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "ShallowElement: CalculateRightHandSide not implemented" << std::endl;
}

//----------------------------------------------------------------------

void ShallowElement::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VEL_ART_VISC ||
        rVariable == PR_ART_VISC ||
        rVariable == RESIDUAL_NORM ||
        rVariable == MIU)
    {
        for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
            rValues[PointNumber] = double(this->GetValue(rVariable));
    }
}

//----------------------------------------------------------------------

void ShallowElement::InitializeElement(ElementData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    // Global values
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    rData.dt_inv = 1.0 / delta_t;
    rData.lumping_factor = 1.0 / msNodes;
    rData.c_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    rData.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rData.manning2 = std::pow( GetProperties()[MANNING], 2);

    // Nodal data
    GeometryType& rGeom = GetGeometry();
    unsigned int counter = 0;
    for (unsigned int i = 0; i < msNodes; i++)
    {
        rData.depth[counter] = 0;
        rData.rain[counter]  = 0;
        rData.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
        rData.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X);
        counter++;

        rData.depth[counter] = 0;
        rData.rain[counter]  = 0;
        rData.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
        rData.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y);
        counter++;

        rData.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY);
        rData.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
        rData.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
        rData.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1);
        counter++;
    }
}

void ShallowElement::ComputeMassMatrices(
    const ElementData& rData,
    BoundedMatrix<double,msElemSize,msElemSize>& rVelMatrix,
    BoundedMatrix<double,msElemSize,msElemSize>& rHeightMatrix)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        rVelMatrix(3*i  ,3*i  ) = 1;
        rVelMatrix(3*i+1,3*i+1) = 1;
        rHeightMatrix(3*i+2,3*i+2) = 1;
    }
    rVelMatrix *= rData.lumping_factor;
    rHeightMatrix *= rData.lumping_factor;
}

void ShallowElement::ComputeElementValues(
    const ElementData& rData,
    array_1d<double,2>& rVel,
    double& rHeight)
{
    // Initialize output
    rVel = ZeroVector(2);
    rHeight = 0;

    // integrate over the element
    for (unsigned int i = 0; i < 3; i++)
    {
        rVel[0] += rData.unknown[  + 3*i];
        rVel[1] += rData.unknown[1 + 3*i];
        rHeight += rData.unknown[2 + 3*i];
    }
    rVel *= rData.lumping_factor;
    rHeight *= rData.lumping_factor;
}


} // namespace kratos