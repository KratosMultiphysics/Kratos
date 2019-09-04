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
#include "rv_swe.h"

namespace Kratos
{

template< size_t TNumNodes, ElementFramework TFramework >
int RV_SWE<TNumNodes, TFramework>::Check( const ProcessInfo& rCurrentProcessInfo )
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
    KRATOS_CHECK_VARIABLE_KEY(WATER_HEIGHT_UNIT_CONVERTER)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( size_t i = 0; i < TNumNodes; i++ )
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


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    size_t element_size = TNumNodes*3;
    if(rResult.size() != element_size)
        rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t element_size = TNumNodes*3;
    if(rElementalDofList.size() != element_size)
        rElementalDofList.resize(element_size);

    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Resize of the Left and Right Hand side
    constexpr size_t element_size = TNumNodes*3;
    if(rLeftHandSideMatrix.size1() != element_size)
        rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != element_size)
        rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!

    const GeometryType& Geom = this->GetGeometry();

    ElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    BoundedMatrix<double,TNumNodes, TNumNodes> NContainer;  // In this case, number of Gauss points and number of nodes coincides
    NContainer = Geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

    BoundedMatrix<double,TNumNodes, 2> DN_DX;  // Shape function gradients are constant since we are using linear functions
    array_1d<double,TNumNodes> N;
    double Area;
    this-> CalculateGeometry(DN_DX,Area);

    rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);
    rRightHandSideVector = ZeroVector(element_size);

    for (size_t i_gauss = 0; i_gauss < TNumNodes; ++i_gauss)
    {
        noalias(N) = row(NContainer, i_gauss);

        BuildMassMatrices(N, variables);
        BuildGradientMatrices(N, DN_DX, variables);
        BuildDiffusivityMatrices(DN_DX, variables);
        BuildConvectionMatrices(N, DN_DX, variables);

        AddInertiaTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        AddConvectiveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        AddWaveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        AddFrictionTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        AddStabilizationTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        AddSourceTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
    }

    rRightHandSideVector *= Area * variables.lumping_factor;
    rLeftHandSideMatrix  *= Area * variables.lumping_factor;
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
    {
        for (size_t PointNumber = 0; PointNumber < 1; PointNumber++)
            rValues[PointNumber] = double(this->GetValue(rVariable));
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    rVariables.dt_inv = 1.0 / delta_t;
    rVariables.lumping_factor = 1.0 / static_cast<double>(TNumNodes);
    rVariables.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    rVariables.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rVariables.manning2 = 0.0;//std::pow( GetProperties()[MANNING], 2);
    rVariables.height_units = rCurrentProcessInfo[WATER_HEIGHT_UNIT_CONVERTER];

    GeometryType& rGeom = GetGeometry();
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.manning2 += rGeom[i].FastGetSolutionStepValue(EQUIVALENT_MANNING);
    }
    rVariables.manning2 *= rVariables.lumping_factor;
    rVariables.manning2 = std::pow(rVariables.manning2, 2);
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
{
    const GeometryType& rGeom = this->GetGeometry();

    // We select GI_GAUSS_1 due to we are computing at the barycenter.
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    const size_t NumGPoints = integration_points.size();
    rArea = rGeom.Area();
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer( NumGPoints );
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, GeometryData::GI_GAUSS_1);

    noalias( rDN_DX ) = DN_DXContainer[0];
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::GetNodalValues(ElementVariables& rVariables)
{
    GeometryType& rGeom = GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.depth[counter] = 0;
        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
        rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X);
        counter++;

        rVariables.depth[counter] = 0;
        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
        rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y);
        counter++;

        rVariables.depth[counter] = rGeom[i].FastGetSolutionStepValue(BATHYMETRY) / rVariables.height_units;
        rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(HEIGHT);
        rVariables.proj_unk[counter]  = rGeom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1);
        counter++;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables)
{
    // Initialize outputs
    rVariables.height = 0;
    rVariables.velocity = ZeroVector(2);
    rVariables.height_grad = ZeroVector(2);

    // integrate over the element
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.velocity[0] += rVariables.unknown[  + 3*i];
        rVariables.velocity[1] += rVariables.unknown[1 + 3*i];
        rVariables.height += rVariables.unknown[2 + 3*i];
        rVariables.height_grad[0] += rDN_DX(i,0) * rVariables.unknown[2 + 3*i];
        rVariables.height_grad[1] += rDN_DX(i,1) * rVariables.unknown[2 + 3*i];
    }

    rVariables.velocity *= rVariables.lumping_factor;
    rVariables.height *= rVariables.lumping_factor * rVariables.height_units;
    rVariables.height_grad *= rVariables.height_units;
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::ComputeStabilizationParameters(const ElementVariables& rVariables,
                                                                    const double& rElemSize,
                                                                    double& rTauU,
                                                                    double& rTauH,
                                                                    double& rKdc)
{
    // Initialize outputs
    rTauU = 0;
    rTauH = 0;
    rKdc  = 0;

    // Get element values
    double height_grad_norm = norm_2(rVariables.height_grad);

    // Compute stabilization parameters
    bool stabilization = true;
    double CTau = rVariables.dyn_tau;  // 0.005 ~ 0.002  // 0.002; //
    double abs_height = std::abs(rVariables.height);
    if (stabilization && abs_height > 1e-6)
    {
        rTauU = CTau/rElemSize*std::sqrt(rVariables.gravity/abs_height);
        rTauH = CTau/rElemSize*std::sqrt(abs_height/rVariables.gravity);
    }

    // Compute discontinuity capturing parameters
    bool discontinuity_capturing = true;
    double gradient_threshold = 1e-6;
    //~ double residual;
    if (discontinuity_capturing && height_grad_norm > gradient_threshold)
    {
        rKdc = 0.5*0.1*rElemSize*height_grad_norm;  // Residual formulation
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::BuildMassMatrices(
    array_1d<double,TNumNodes>& rN,
    ElementVariables& rVariables)
{
    // Auxiliary definitions
    BoundedMatrix<double,2,TNumNodes*3> N_vel    = ZeroMatrix(2,TNumNodes*3); // Shape functions matrix (for velocity unknown)
    array_1d<double,TNumNodes*3> N_height        = ZeroVector(TNumNodes*3);  // Shape functions vector (for height unknown)

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Height shape funtions
        N_height[2+node*3] = rN[node];
        // Velocity shape functions
        N_vel(0,   node*3) = rN[node];
        N_vel(1, 1+node*3) = rN[node];
    }
    noalias(rVariables.MassMatrixVector) += prod(trans(N_vel),N_vel);       // q * h
    noalias(rVariables.MassMatrixScalar) += outer_prod(N_height,N_height); // w * u
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::BuildGradientMatrices(
    array_1d<double,TNumNodes>& rN,
    BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
    ElementVariables& rVariables)
{
    // Auxiliary definitions
    BoundedMatrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3); // Shape functions matrix (for velocity unknown)
    array_1d<double,TNumNodes*3> N_height            = ZeroVector(TNumNodes*3);   // Shape functions vector (for height unknown)
    array_1d<double,TNumNodes*3> DN_DX_vel           = ZeroVector(TNumNodes*3);   // Shape functions gradients vector (for velocity unknown)
    BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3); // Shape functions gradients matrix (for height unknown)

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Height gradient
        DN_DX_height(0, 2+node*3) = rDN_DX(node,0);
        DN_DX_height(1, 2+node*3) = rDN_DX(node,1);
        // Velocity divergence
        DN_DX_vel[  node*3] = rDN_DX(node,0);
        DN_DX_vel[1+node*3] = rDN_DX(node,1);
        // Height shape funtions
        N_height[2+node*3] = rN[node];
        // Velocity shape functions
        N_vel(0,   node*3) = rN[node];
        N_vel(1, 1+node*3) = rN[node];
    }
    noalias(rVariables.VectorDiv)  += outer_prod(N_height,DN_DX_vel); // q * div_u
    noalias(rVariables.ScalarGrad) += prod(trans(N_vel),DN_DX_height); // w * grad_h
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::BuildDiffusivityMatrices(
    BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
    ElementVariables& rVariables)
{
    // Some auxilary definitions
    array_1d<double,TNumNodes*3> DN_DX_vel           = ZeroVector(TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
    BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Height gradient
        DN_DX_height(0, 2+node*3) = rDN_DX(node,0);
        DN_DX_height(1, 2+node*3) = rDN_DX(node,1);
        // Velocity divergence
        DN_DX_vel[  node*3] = rDN_DX(node,0);
        DN_DX_vel[1+node*3] = rDN_DX(node,1);
    }
    noalias(rVariables.VectorDiff) += outer_prod(DN_DX_vel,DN_DX_vel);       // div_w * div_u
    noalias(rVariables.ScalarDiff) += prod(trans(DN_DX_height),DN_DX_height); // grad_q * grad_h
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::BuildConvectionMatrices(
    array_1d<double, TNumNodes>& rN,
    BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
    ElementVariables& rVariables)
{
    // Auxiliary definitions
    BoundedMatrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
    array_1d<double,TNumNodes*3> N_height            = ZeroVector(TNumNodes*3);    // Shape functions vector (for height unknown)
    BoundedMatrix<double,2,TNumNodes*3> Grad_vel_1   = ZeroMatrix(1,TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
    BoundedMatrix<double,2,TNumNodes*3> Grad_vel_2   = ZeroMatrix(1,TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
    BoundedMatrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Height gradient
        DN_DX_height(0, 2+node*3) = rDN_DX(node,0);
        DN_DX_height(1, 2+node*3) = rDN_DX(node,1);
        // Velocity gradient
        Grad_vel_1(0,   node*3) = rDN_DX(node,0);
        Grad_vel_1(1, 1+node*3) = rDN_DX(node,0);
        Grad_vel_2(0,   node*3) = rDN_DX(node,1);
        Grad_vel_2(1, 1+node*3) = rDN_DX(node,1);
        // Height shape funtions
        N_height[2+node*3] = rN[node];
        // Velocity shape functions
        N_vel(0,   node*3) = rN[node];
        N_vel(1, 1+node*3) = rN[node];
    }

    array_1d<double,TNumNodes*3> temp_convect = prod(rVariables.velocity, DN_DX_height);
    noalias(rVariables.Convection) += outer_prod(N_height, temp_convect);     // q * u * grad_h
    noalias(rVariables.Convection) += rVariables.velocity[0] * prod(trans(N_vel), Grad_vel_1);     // w * u * grad_u
    noalias(rVariables.Convection) += rVariables.velocity[1] * prod(trans(N_vel), Grad_vel_2);     // w * u * grad_u
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddInertiaTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    BoundedMatrix<double, TNumNodes*3,TNumNodes*3> mass_matrix = rVariables.MassMatrixVector + rVariables.MassMatrixScalar;
    rLeftHandSideMatrix += rVariables.dt_inv * mass_matrix;
    rRightHandSideVector += rVariables.dt_inv * prod(mass_matrix, rVariables.prev_unk);
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddConvectiveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    if (TFramework == Eulerian)
    {
        double a = 1;
        a++;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddWaveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddFrictionTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddStabilizationTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::AddSourceTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{}


template< size_t TNumNodes, ElementFramework TFramework >
void RV_SWE<TNumNodes, TFramework>::CalculateLumpedMassMatrix(BoundedMatrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix)
{
    const size_t element_size = 3*TNumNodes;
    rMassMatrix  = IdentityMatrix(element_size, element_size);
    rMassMatrix /= static_cast<double>(TNumNodes);
}


template class RV_SWE<3, Eulerian>;
template class RV_SWE<4, Eulerian>;
template class RV_SWE<3, PFEM2>;
template class RV_SWE<4, PFEM2>;

} // namespace Kratos
