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

    const std::size_t block_size = mDim + 1;
    const std::size_t mat_size = number_of_nodes * block_size;
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    const SizeType strain_size = (mDim == 3) ? 6 : 3;

    // Compute the B matrix
    Matrix B = ZeroMatrix(strain_size, number_of_nodes * mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(strain_size);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();

    Matrix DB_voigt = Matrix(prod(r_D, B));

    // Compute the normals
    array_1d<double, 3> normal_physical_space;
    array_1d<double, 3> normal_parameter_space = -r_geometry.Normal(0, GetIntegrationMethod());
    if (mDim == 3) {
        r_geometry.Calculate(NORMAL, normal_parameter_space);
    }
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);
    normal_physical_space = normal_parameter_space;

    const Matrix& H = r_geometry.ShapeFunctionsValues();

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0, r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is 3x2 and we need 3x3
    Matrix Jacobian = ZeroMatrix(3, 3);
    Jacobian(0, 0) = J0[0](0, 0);
    Jacobian(0, 1) = J0[0](0, 1);
    Jacobian(1, 0) = J0[0](1, 0);
    Jacobian(1, 1) = J0[0](1, 1);
    Jacobian(2, 2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(Jacobian, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    double det_J0 = norm_2(determinant_factor);

    if (mDim == 3) {
        Matrix tangent_matrix;
        r_geometry.Calculate(LOCAL_TANGENT_MATRIX, tangent_matrix);  // 3x2

        array_1d<double,3> t1, t2;
        for (std::size_t i = 0; i < 3; ++i) {
            t1[i] = tangent_matrix(i, 0);
            t2[i] = tangent_matrix(i, 1);
        }
        // Cross product of the two tangents
        array_1d<double, 3> det_vector = MathUtils<double>::CrossProduct(t1, t2);
        // Norm gives the surface integration factor
        det_J0 = norm_2(det_vector);
    }

    double penalty_integration = mPenalty * r_integration_points[0].Weight() * std::abs(det_J0);
    const double integration_weight = r_integration_points[0].Weight() * std::abs(det_J0);

    // Compute the pressure & velocity at the previous iteration
    double pressure_current_iteration = 0.0;
    Vector velocity_current_iteration = ZeroVector(mDim);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        pressure_current_iteration    += r_geometry[j].GetSolutionStepValue(PRESSURE) * H(0,j);
        const auto& r_velocity = r_geometry[j].GetSolutionStepValue(VELOCITY);
        for (IndexType d = 0; d < mDim; ++d) {
            velocity_current_iteration[d] += r_velocity[d] * H(0,j);
        }
    }

    Vector n_tensor(mDim);
    for (IndexType d = 0; d < mDim; ++d) {
        n_tensor(d) = normal_parameter_space(d);
    }

    // Compute the traction vector: sigma * n using r_stress_vector
    Matrix stress_old = ZeroMatrix(mDim, mDim);
    if (mDim == 2) {
        stress_old(0, 0) = r_stress_vector[0];
        stress_old(1, 1) = r_stress_vector[1];
        stress_old(0, 1) = r_stress_vector[2];
        stress_old(1, 0) = r_stress_vector[2];
    } else {
        // 3D Voigt order: xx, yy, zz, xy, yz, xz.
        stress_old(0, 0) = r_stress_vector[0];
        stress_old(1, 1) = r_stress_vector[1];
        stress_old(2, 2) = r_stress_vector[2];
        stress_old(0, 1) = r_stress_vector[3];
        stress_old(1, 0) = r_stress_vector[3];
        stress_old(1, 2) = r_stress_vector[4];
        stress_old(2, 1) = r_stress_vector[4];
        stress_old(0, 2) = r_stress_vector[5];
        stress_old(2, 0) = r_stress_vector[5];
    }
    Vector traction_current_iteration = prod(stress_old, n_tensor); 

    Matrix DB_contribution_w = ZeroMatrix(mDim, mDim);
    Matrix DB_contribution = ZeroMatrix(mDim, mDim);

    const auto BuildStressFromVoigtColumn = [&](Matrix& rSigma, const IndexType Column) {
        noalias(rSigma) = ZeroMatrix(mDim, mDim);
        if (mDim == 2) {
            rSigma(0, 0) = DB_voigt(0, Column);
            rSigma(1, 1) = DB_voigt(1, Column);
            rSigma(0, 1) = DB_voigt(2, Column);
            rSigma(1, 0) = DB_voigt(2, Column);
        } else {
            // 3D Voigt order: xx, yy, zz, xy, yz, xz.
            rSigma(0, 0) = DB_voigt(0, Column);
            rSigma(1, 1) = DB_voigt(1, Column);
            rSigma(2, 2) = DB_voigt(2, Column);
            rSigma(0, 1) = DB_voigt(3, Column);
            rSigma(1, 0) = DB_voigt(3, Column);
            rSigma(1, 2) = DB_voigt(4, Column);
            rSigma(2, 1) = DB_voigt(4, Column);
            rSigma(0, 2) = DB_voigt(5, Column);
            rSigma(2, 0) = DB_voigt(5, Column);
        }
    };
    
    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const IndexType col_w = i * mDim + idim;
            BuildStressFromVoigtColumn(DB_contribution_w, col_w);
            Vector traction_nitsche_w = prod(DB_contribution_w, n_tensor);

            for (IndexType j = 0; j < number_of_nodes; j++) {  
                // Penalty term for the velocity
                rLeftHandSideMatrix(i * block_size + idim, j * block_size + idim) +=
                    H(0,i) * H(0,j) * penalty_integration;
                
                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const IndexType col = j * mDim + jdim;
                    BuildStressFromVoigtColumn(DB_contribution, col);
                    Vector traction = prod(DB_contribution, n_tensor);

                    // integration by parts velocity --> With Constitutive law  < v cdot (DB cdot n) >
                    rLeftHandSideMatrix(i * block_size + idim, j * block_size + jdim) -=
                        H(0, i) * traction(idim) * integration_weight;
                    
                    // Nitsche term --> With Constitutive law
                    rLeftHandSideMatrix(i * block_size + idim, j * block_size + jdim) +=
                        H(0, j) * traction_nitsche_w(jdim) * integration_weight;
                }

                // integration by parts PRESSURE
                rLeftHandSideMatrix(i * block_size + idim, j * block_size + mDim) +=
                    H(0,j) * ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;

                // Nitsche term --> q term
                rLeftHandSideMatrix(j * block_size + mDim, i * block_size + idim) -=
                    H(0,j) * ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;

            }
        }

        // --- RHS corresponding term ---
        for (IndexType idim = 0; idim < mDim; idim++) {
            // Penalty term for the velocity
            rRightHandSideVector(i * block_size + idim) -=
                H(0,i) * velocity_current_iteration[idim] * penalty_integration;
            // integration by parts velocity --> With Constitutive law
            rRightHandSideVector(i * block_size + idim) +=
                H(0,i) * traction_current_iteration(idim) * integration_weight;
            // integration by parts PRESSURE
            rRightHandSideVector(i * block_size + idim) -=
                pressure_current_iteration * ( H(0,i) * normal_parameter_space[idim] ) * integration_weight;
            
            // Nitsche term --> With Constitutive law
            Matrix DB_contribution = ZeroMatrix(mDim, mDim);
            BuildStressFromVoigtColumn(DB_contribution, i * mDim + idim);
            Vector traction = prod(DB_contribution, n_tensor);
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector(i * block_size + idim) -=
                    velocity_current_iteration[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector(i * block_size + mDim) +=
                velocity_current_iteration[idim] * ( H(0,i) * normal_parameter_space[idim] )
                        * integration_weight;
            
        }
    }
            
    Vector u_D = ZeroVector(mDim); 
    u_D[0] = this->GetValue(VELOCITY_X);
    u_D[1] = this->GetValue(VELOCITY_Y);
    if (mDim == 3) {
        u_D[2] = this->GetValue(VELOCITY_Z);
    }

    // // TO BE DELETED
    // if (mDim == 3) {
    //     if (r_geometry.Center().X() < -0.49999) {
    //         // Inlet
    //         const double t = rCurrentProcessInfo[TIME];
    //         double ux = 1.0;
    //         if (t < 2.0) {
    //             ux = 0.5 * t; // linear ramp: ux = t / 2
    //         }
    //         u_D[0] = ux;
    //     }
    // }

    for (IndexType i = 0; i < number_of_nodes; i++) {

        for (IndexType idim = 0; idim < mDim; idim++) {
            
            // Penalty term for the velocity
            rRightHandSideVector[i * block_size + idim] +=
                H(0,i) * u_D[idim] * penalty_integration;

            // Extract the 2x2 block for the control point i from the sigma matrix.
            Matrix sigma_block = ZeroMatrix(mDim, mDim);
            BuildStressFromVoigtColumn(sigma_block, i * mDim + idim);
            Vector traction = prod(sigma_block, n_tensor); // This results in a 2x1 vector.
            // Nitsche term --> With Constitutive law
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[i * block_size + idim] +=
                    u_D[jdim] * traction(jdim) * integration_weight;
            }

            // Nitsche term --> q term
            rRightHandSideVector[i * block_size + mDim] -=
                u_D[idim] * H(0,i) * normal_parameter_space[idim] * integration_weight;

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
    const std::size_t mat_size = number_of_control_points * mDim;
    const std::size_t strain_size = (mDim == 3) ? 6 : 3;

    // Resize B matrix to Voigt strain size and appropriate number of columns.
    if (rB.size1() != strain_size || rB.size2() != mat_size)
        rB.resize(strain_size, mat_size);

    noalias(rB) = ZeroMatrix(strain_size, mat_size);

    if (mDim == 2) {
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

        for (IndexType d = 0; d < mDim; ++d) {
            rValues[index + d] = velocity[d];
        }
    }
}

} // Namespace Kratos
