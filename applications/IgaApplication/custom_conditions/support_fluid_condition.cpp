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
#include "custom_conditions/support_fluid_condition.h"

namespace Kratos
{

void SupportFluidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeMaterial();
}

void SupportFluidCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const std::size_t number_of_nodes = r_geometry.size();

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    Matrix DN_DX(number_of_nodes,mDim);
    noalias(DN_DX) = DN_De[0]; // prod(DN_De[point_number],InvJ0);

    const std::size_t mat_size = number_of_nodes * (mDim+1);
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    // Compute the B matrix
    Matrix B = ZeroMatrix(3,number_of_nodes*mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(3);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix DB_voigt = Matrix(prod(r_D, B));

    // Compute the normals
    array_1d<double, 3> normal_physical_space;
    // Compute the normals
    array_1d<double, 3> normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);
    normal_physical_space = normal_parameter_space;

    const Matrix& H = r_geometry.ShapeFunctionsValues();
    
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is  3x2 and we need 3x3
    Matrix jacobian_matrix = ZeroMatrix(3,3);
    jacobian_matrix(0,0) = J0[0](0,0);
    jacobian_matrix(0,1) = J0[0](0,1);
    jacobian_matrix(1,0) = J0[0](1,0);
    jacobian_matrix(1,1) = J0[0](1,1);
    jacobian_matrix(2,2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(jacobian_matrix, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    const double det_J0 = norm_2(determinant_factor);
    
    double penalty_integration = mPenalty * r_integration_points[0].Weight() * std::abs(det_J0);
    const double integration_weight = r_integration_points[0].Weight() * std::abs(det_J0);

    // Compute the pressure & velocity at the previous iteration
    double pressure_current_iteration = 0.0;
    Vector velocity_current_iteration = ZeroVector(2);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        pressure_current_iteration    += r_geometry[j].GetSolutionStepValue(PRESSURE) * H(0,j);
        velocity_current_iteration[0] += r_geometry[j].GetSolutionStepValue(VELOCITY_X) * H(0,j);
        velocity_current_iteration[1] += r_geometry[j].GetSolutionStepValue(VELOCITY_Y) * H(0,j);
    }

    Vector n_tensor(2);
    n_tensor(0) = normal_parameter_space(0); // Component in x direction
    n_tensor(1) = normal_parameter_space(1); // Component in y direction

    // Compute the traction vector: sigma * n using r_stress_vector
    Matrix stress_old = ZeroMatrix(2, 2);
    stress_old(0, 0) = r_stress_vector[0];      stress_old(0, 1) = r_stress_vector[2];      
    stress_old(1, 0) = r_stress_vector[2];      stress_old(1, 1) = r_stress_vector[1];         
    Vector traction_current_iteration = prod(stress_old, n_tensor); 

    Matrix DB_contribution_w = ZeroMatrix(2, 2);
    Matrix DB_contribution = ZeroMatrix(2, 2);
    
    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            DB_contribution_w(0, 0) = DB_voigt(0, 2*i+idim);
            DB_contribution_w(0, 1) = DB_voigt(2, 2*i+idim);
            DB_contribution_w(1, 0) = DB_voigt(2, 2*i+idim);
            DB_contribution_w(1, 1) = DB_voigt(1, 2*i+idim);
            // Compute the traction vector: sigma * n.
            Vector traction_nitsche_w = prod(DB_contribution_w, n_tensor);

            for (IndexType j = 0; j < number_of_nodes; j++) {  
                // Penalty term for the velocity
                rLeftHandSideMatrix(3*i+idim, 3*j+idim) += H(0,i)*H(0,j)* penalty_integration;
                
                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    // Extract the 2x2 block for the control point i from the DB_voigt.
                    DB_contribution(0, 0) = DB_voigt(0, 2*j+jdim);
                    DB_contribution(0, 1) = DB_voigt(2, 2*j+jdim);
                    DB_contribution(1, 0) = DB_voigt(2, 2*j+jdim);
                    DB_contribution(1, 1) = DB_voigt(1, 2*j+jdim);
                    // Compute the traction vector: sigma * n.
                    Vector traction = prod(DB_contribution, n_tensor);

                    // integration by parts velocity --> With Constitutive law  < v cdot (DB cdot n) >
                    rLeftHandSideMatrix(3*i+idim, 3*j+jdim) -= H(0, i) * traction(idim) * integration_weight;
                    
                    // Nitsche term --> With Constitutive law
                    rLeftHandSideMatrix(3*i+idim, 3*j+jdim) += H(0, j) * traction_nitsche_w(jdim) * integration_weight;
                }

                // integration by parts PRESSURE
                rLeftHandSideMatrix(3*i+idim, 3*j+2) += H(0,j)* ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;

                // Nitsche term --> q term
                rLeftHandSideMatrix(3*j+mDim, 3*i+idim) -= H(0,j)* ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;

            }
        }

        // --- RHS corresponding term ---
        for (IndexType idim = 0; idim < 2; idim++) {
            // Penalty term for the velocity
            rRightHandSideVector(3*i+idim) -= H(0,i)* velocity_current_iteration[idim] * penalty_integration;
            // integration by parts velocity --> With Constitutive law
            rRightHandSideVector(3*i+idim) += H(0,i) * traction_current_iteration(idim) * integration_weight;
            // integration by parts PRESSURE
            rRightHandSideVector(3*i+idim) -= pressure_current_iteration * ( H(0,i) * normal_parameter_space[idim] ) * integration_weight;
            
            // Nitsche term --> With Constitutive law
            Matrix DB_contribution = ZeroMatrix(2, 2); // Extract the 2x2 block for the control point i from the DB_voigt.
            DB_contribution(0, 0) = DB_voigt(0, 2*i+idim); 
            DB_contribution(0, 1) = DB_voigt(2, 2*i+idim); 
            DB_contribution(1, 0) = DB_voigt(2, 2*i+idim);
            DB_contribution(1, 1) = DB_voigt(1, 2*i+idim); 
            // Compute the traction vector: sigma * n.
            Vector traction = prod(DB_contribution, n_tensor);
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector(3*i+idim) -= velocity_current_iteration[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector(3*i+mDim) += velocity_current_iteration[idim] * ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;
            
        }
    }
            
    Vector u_D = ZeroVector(2); 
    u_D[0] = this->GetValue(VELOCITY_X);
    u_D[1] = this->GetValue(VELOCITY_Y);

    for (IndexType i = 0; i < number_of_nodes; i++) {

        for (IndexType idim = 0; idim < 2; idim++) {
            
            // Penalty term for the velocity
            rRightHandSideVector[3*i+idim] += H(0,i) * u_D[idim] * penalty_integration;

            // Extract the 2x2 block for the control point i from the sigma matrix.
            Matrix sigma_block = ZeroMatrix(2, 2);
            sigma_block(0, 0) = DB_voigt(0, 2*i+idim);
            sigma_block(0, 1) = DB_voigt(2, 2*i+idim);
            sigma_block(1, 0) = DB_voigt(2, 2*i+idim);
            sigma_block(1, 1) = DB_voigt(1, 2*i+idim);
            // Compute the traction vector: sigma * n.
            Vector traction = prod(sigma_block, n_tensor); // This results in a 2x1 vector.
            // Nitsche term --> With Constitutive law
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[3*i+idim] += u_D[jdim] * traction(jdim) * integration_weight;
            }

            // Nitsche term --> q term
            rRightHandSideVector[3*i+mDim] -= u_D[idim] * H(0,i)*normal_parameter_space[idim] * integration_weight;

        }
    }
    KRATOS_CATCH("")
}

void SupportFluidCondition::InitializeMemberVariables()
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

    mPenalty = GetProperties()[PENALTY_FACTOR];
    KRATOS_ERROR_IF(mPenalty == -1.0)
        << "Penalty-free formulation is not available for the Stokes problem" << std::endl;

    // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
    mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * mPenalty / h;
}

void SupportFluidCondition::CalculateB(
        Matrix& rB, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

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

void SupportFluidCondition::ApplyConstitutiveLaw(
        const Matrix& rB, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const
{
    const std::size_t number_of_nodes = GetGeometry().size();

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    Vector old_displacement(number_of_nodes*mDim);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(rB,old_displacement);
    rValues.SetStrainVector(old_strain);
    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void SupportFluidCondition::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

int SupportFluidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of SupportFluidCondition" << std::endl;
    return 0;
}


void SupportFluidCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const std::size_t number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (mDim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (mDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SupportFluidCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const std::size_t number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((mDim+1) * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        if (mDim > 2) rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SupportFluidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * mDim;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
        if (mDim > 3) rValues[index + 2] = velocity[2];
    }
}

} // Namespace Kratos
