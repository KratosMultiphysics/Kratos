//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_pressure_condition.h"

namespace Kratos
{

void SupportPressureCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
}

void SupportPressureCondition::CalculateLocalSystem(
MatrixType& rLeftHandSideMatrix,
VectorType& rRightHandSideVector,
const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes * (mDim+1);

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SupportPressureCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetGeometry().size() * (mDim+1);

    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2())
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
}

void SupportPressureCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    Matrix DN_DX(number_of_nodes,mDim);

    const SizeType mat_size = number_of_nodes * (mDim+1);
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    // Compute the normals
    array_1d<double, 3> normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);

    // Matrix DB = prod(r_D,B);
    const Matrix& N = r_geometry.ShapeFunctionsValues();
    // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0];

    double p_D = this->GetValue(PRESSURE);

    for (IndexType j = 0; j < number_of_nodes; j++) {
        
        //// integration by parts PRESSURE < v dot n , p >
        for (IndexType idim = 0; idim < mDim; idim++) {
            rRightHandSideVector((mDim+1)*j+idim) -= p_D * ( N(0,j) * normal_parameter_space[idim] ) * mIntegrationWeight;
        }
        
        // Neumann condition for the velocity
        Vector t_N = this->GetValue(NORMAL_STRESS);

        for (IndexType idim = 0; idim < mDim; idim++) {
            rRightHandSideVector((mDim+1)*j+idim) += N(0,j) * t_N[idim] * mIntegrationWeight;
        }
        
    }
    KRATOS_CATCH("")
}


void SupportPressureCondition::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    KRATOS_ERROR_IF(mDim == 3) << "SupportPressureCondition is not implemented in 3D. Current dimension: " << mDim << std::endl;
    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    GeometryType::JacobiansType J0;
    Matrix InvJ0(mDim,mDim);
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    Matrix Jacobian = ZeroMatrix(3,3);
    double DetJ0;
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector add_factor = prod(Jacobian, tangent_parameter_space);
    add_factor[2] = 0.0;
    DetJ0 = norm_2(add_factor);
    mIntegrationWeight = r_integration_points[0].Weight() * std::abs(DetJ0);
}

void SupportPressureCondition::CalculateB(
    Matrix& rB, 
    const ShapeDerivativesType& r_DN_DX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (rB.size1() != 3 || rB.size2() != mat_size)
        rB.resize(3, mat_size);

    noalias(rB) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // x-derivatives of shape functions -> relates to strain component ε_11 (xx component)
        rB(0, 2 * i)     = r_DN_DX(i, 0); // ∂N_i / ∂x
        // y-derivatives of shape functions -> relates to strain component ε_22 (yy component)
        rB(1, 2 * i + 1) = r_DN_DX(i, 1); // ∂N_i / ∂y
        // Symmetric shear strain component ε_12 (xy component)
        rB(2, 2 * i)     = r_DN_DX(i, 1); // ∂N_i / ∂y
        rB(2, 2 * i + 1) = r_DN_DX(i, 0); // ∂N_i / ∂x
    }
}


void SupportPressureCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (mDim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SupportPressureCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(3 * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SupportPressureCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3>& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}

} // Namespace Kratos
