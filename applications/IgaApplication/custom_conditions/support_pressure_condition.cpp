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

    const std::size_t number_of_nodes = GetGeometry().size();
    const std::size_t mat_size = number_of_nodes * (mDim+1);

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
    const std::size_t mat_size = GetGeometry().size() * (mDim+1);

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
    const std::size_t number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    Matrix DN_DX(number_of_nodes,mDim);

    const std::size_t mat_size = number_of_nodes * (mDim+1);
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    // Compute the normals
    array_1d<double, 3> normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
    if (mDim == 3) {
        r_geometry.Calculate(NORMAL, normal_parameter_space);
    }
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
        
        // // Neumann condition for the velocity
        // Vector t_N = this->GetValue(NORMAL_STRESS);
        Vector t_N = ZeroVector(mDim);
        double x = r_geometry.Center().X(); 
        double y = r_geometry.Center().Y();
        double z = r_geometry.Center().Z();

        Matrix grad_u(mDim, mDim); // grad_u(i,j) = du_i/dx_j

        if (mDim == 2) {
            grad_u(0, 0) = sinh(x) * sinh(y);     // ∂u_x / ∂x
            grad_u(0, 1) = cosh(x) * cosh(y);     // ∂u_x / ∂y
            grad_u(1, 0) = -cosh(x) * cosh(y);    // ∂u_y / ∂x
            grad_u(1, 1) = -sinh(x) * sinh(y);    // ∂u_y / ∂y
        } else {
            // 3D analytical gradient
            noalias(grad_u) = ZeroMatrix(mDim, mDim);

            // grad_u(0,0) =  2.0;  // dux/dx
            // grad_u(0,1) = -1.0;  // dux/dy
            // grad_u(0,2) =  3.0;  // dux/dz
            // grad_u(1,0) =  3.0;  // duy/dx
            // grad_u(1,1) =  2.0;  // duy/dy
            // grad_u(1,2) = -1.0;  // duy/dz
            // grad_u(2,0) = -5.0;  // duz/dx
            // grad_u(2,1) = -4.0;  // duz/dy
            // grad_u(2,2) = -4.0;  // duz/dz

            // // quadratic
            // grad_u(0,0) = -2.0*y;  grad_u(0,1) = -2.0*x;  grad_u(0,2) = 0.0;
            // grad_u(1,0) = -2.0*x;  grad_u(1,1) =  2.0*y + z;  grad_u(1,2) = y - z;
            // grad_u(2,0) =  0.0;    grad_u(2,1) =  0.0;    grad_u(2,2) = -z;

            // // cubic
            // grad_u(0,0) = 2.0*y*z - 2.0*x*y;   grad_u(0,1) = 2.0*x*z - x*x;      grad_u(0,2) = 2.0*x*y;
            // grad_u(1,0) = 0.0;                 grad_u(1,1) = -2.0*y*z;           grad_u(1,2) = -y*y;
            // grad_u(2,0) = 2.0*y*z;             grad_u(2,1) = 2.0*x*z;            grad_u(2,2) = 2.0*x*y;

            // ux = cosh(x)cosh(y)cosh(z)          // uy = sinh(x)*sinh(z-y)          // uz = -sinh(x)sinh(y)cosh(z)
            grad_u(0,0) = sinh(x)*cosh(y)*cosh(z); grad_u(1,0) = cosh(x)*sinh(z-y);   grad_u(2,0) = -cosh(x)*sinh(y)*cosh(z);
            grad_u(0,1) = cosh(x)*sinh(y)*cosh(z); grad_u(1,1) = -sinh(x)*cosh(z-y);  grad_u(2,1) = -sinh(x)*cosh(y)*cosh(z);
            grad_u(0,2) = cosh(x)*cosh(y)*sinh(z); grad_u(1,2) =  sinh(x)*cosh(z-y);  grad_u(2,2) = -sinh(x)*sinh(y)*sinh(z);
        }

        // grad_u(0, 0) = 3.0 * x * x;      // ∂u_x / ∂x
        // grad_u(0, 1) = 0.0;              // ∂u_x / ∂y
        // grad_u(1, 0) = -6.0 * x * y;     // ∂u_y / ∂x
        // grad_u(1, 1) = -3.0 * x * x;     // ∂u_y / ∂y

        // grad_u(0, 0) = 2.0 * x;      // ∂u_x / ∂x
        // grad_u(0, 1) = 0.0;              // ∂u_x / ∂y
        // grad_u(1, 0) = -2.0 * y;     // ∂u_y / ∂x
        // grad_u(1, 1) = -2.0 * x;     // ∂u_y / ∂y

        // grad_u(0, 0) = sinh(x) * sinh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_x / ∂x
        // grad_u(0, 1) = cosh(x) * cosh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_x / ∂y
        // grad_u(1, 0) = -cosh(x) * cosh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_y / ∂x
        // grad_u(1, 1) = -sinh(x) * sinh(y) * std::exp(-current_time)*current_time*current_time;     // ∂u_y / ∂y

        Matrix sym_grad_u(mDim, mDim); // ε(u) = 0.5*(∇u + ∇u^T)
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                sym_grad_u(i,j) = 0.5 * (grad_u(i,j) + grad_u(j,i));
            }
        }
        // Now compute stress vector: σ·n = 2ν ε(u)·n
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                t_N[i] += 2.0 * sym_grad_u(i,j) * normal_parameter_space[j];
            }
        }

        // -------------------------------------------------------CHANNEL-------------------------------------------------------
        t_N = ZeroVector(mDim); 
        // -------------------------------------------------------CHANNEL-------------------------------------------------------

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

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    const double det_J0 = std::abs(r_geometry.DeterminantOfJacobian(0, r_geometry.GetDefaultIntegrationMethod()));
    mIntegrationWeight = r_integration_points[0].Weight() * std::abs(det_J0);
}

void SupportPressureCondition::CalculateB(
    Matrix& rB, 
    const ShapeDerivativesType& r_DN_DX) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * mDim;
    const std::size_t strain_size = (mDim == 3) ? 6 : 3;

    // Resize B matrix to Voigt strain size and appropriate number of columns.
    if (rB.size1() != strain_size || rB.size2() != mat_size)
        rB.resize(strain_size, mat_size);

    noalias(rB) = ZeroMatrix(strain_size, mat_size);

    if (mDim == 2) {
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 2 * i)     = r_DN_DX(i, 0);
            rB(1, 2 * i + 1) = r_DN_DX(i, 1);
            rB(2, 2 * i)     = r_DN_DX(i, 1);
            rB(2, 2 * i + 1) = r_DN_DX(i, 0);
        }
    } else {
        // 3D small-strain Voigt order: xx, yy, zz, xy, yz, xz.
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 3 * i)     = r_DN_DX(i, 0);
            rB(1, 3 * i + 1) = r_DN_DX(i, 1);
            rB(2, 3 * i + 2) = r_DN_DX(i, 2);
            rB(3, 3 * i)     = r_DN_DX(i, 1);
            rB(3, 3 * i + 1) = r_DN_DX(i, 0);
            rB(4, 3 * i + 1) = r_DN_DX(i, 2);
            rB(4, 3 * i + 2) = r_DN_DX(i, 1);
            rB(5, 3 * i)     = r_DN_DX(i, 2);
            rB(5, 3 * i + 2) = r_DN_DX(i, 0);
        }
    }
}


void SupportPressureCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t LocalSize = (mDim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    IndexType Index = 0;

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (mDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SupportPressureCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const std::size_t number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((mDim + 1) * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        if (mDim > 2) rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SupportPressureCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3>& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * mDim;

        for (IndexType d = 0; d < mDim; ++d) {
            rValues[index + d] = velocity[d];
        }
    }
}

} // Namespace Kratos
