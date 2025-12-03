
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andea Gorgi
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_solid_condition.h"
// #define SWITCH_OFF_NISCHE  

#define ANALYTICAL_TANGENT

namespace Kratos
{

void SupportSolidCondition::CalculateDeltaPositionMatrix(
    const GeometryType& rGeometry,
    Matrix& rDeltaPosition)  
{
    const SizeType number_of_points = rGeometry.PointsNumber();
    if (rDeltaPosition.size1() != number_of_points || rDeltaPosition.size2() != 3) {
        rDeltaPosition.resize(number_of_points, 3, false);
    }

    for (IndexType i = 0; i < number_of_points; ++i) {
        const auto& r_node = rGeometry[i];
        rDeltaPosition(i, 0) = r_node.X() - r_node.X0();
        rDeltaPosition(i, 1) = r_node.Y() - r_node.Y0();
        rDeltaPosition(i, 2) = r_node.Z() - r_node.Z0();
    }
}


void SupportSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
}


void SupportSolidCondition::InitializeMaterial()
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

void SupportSolidCondition::InitializeMemberVariables()
{
    // calculate the integration weight
    const auto& r_geometry = GetGeometry();
    // Initialize the dimension
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = r_DN_De[0].size2();
    KRATOS_ERROR_IF(mDim != 2) << "SolidElement momentarily only supports 2D elements, but the current element has dimension " << mDim << std::endl;
    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Initialize Jacobian
    Matrix InvJ0(3,3);
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // Compute the normals
    array_1d<double, 3> tangent_parameter_space;
    array_1d<double, 3> normal_physical_space;
    array_1d<double, 3> normal_parameter_space;

    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);

    const double abs_tangent_u = std::abs(tangent_parameter_space[0]);
    const double abs_tangent_v = std::abs(tangent_parameter_space[1]);
    mCharacteristicGeometryLength = mesh_size_uv[0];
    if (abs_tangent_v > abs_tangent_u) {
        mCharacteristicGeometryLength = mesh_size_uv[1];
    }
    mCharacteristicGeometryLength *= 0.5;

    double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
    
    normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
    normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT)
    normal_parameter_space[2] = 0.0;

    // compute complete jacobian transformation including parameter->physical space transformation
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0, 0) = J0[0](0, 0);
    Jacobian(0, 1) = J0[0](0, 1);
    Jacobian(1, 0) = J0[0](1, 0);
    Jacobian(1, 1) = J0[0](1, 1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    Vector add_factor = prod(Jacobian, tangent_parameter_space); //additional factor to the determinant of the jacobian for the parameter->physical space transformation
    add_factor[2] = 0.0; 
    DetJ0 = norm_2(add_factor);

    // compute the normal of the physical space
    normal_physical_space = prod(trans(InvJ0),normal_parameter_space);
    normal_physical_space[2] = 0.0;
    normal_physical_space /= norm_2(normal_physical_space);
    SetValue(NORMAL, normal_physical_space);

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double integration_weight = r_integration_points[0].Weight() * std::abs(DetJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void SupportSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const SizeType mat_size = GetGeometry().size() * 2;

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

void SupportSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);
    const double penalty = GetProperties()[PENALTY_FACTOR];

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // compute complete jacobian transformation including parameter->physical space transformation
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0, 0) = J0[0](0, 0);
    Jacobian(0, 1) = J0[0](0, 1);
    Jacobian(1, 0) = J0[0](1, 0);
    Jacobian(1, 1) = J0[0](1, 1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    // retrieve the normal of the physical space
    array_1d<double, 3> normal_physical_space = GetValue(NORMAL);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);
    noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);

    Matrix B = ZeroMatrix(mDim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    #ifdef ANALYTICAL_TANGENT
    Matrix r_D = ZeroMatrix(3,3);
    AnalyticalConstitutiveMatrix(r_D, old_strain);
    #else
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    #endif

    // Differential area
    double penalty_integration = penalty * integration_weight;

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = 1.0;  // = 1 -> Penalty approach
                                        // = -1 -> Free-penalty approach
    if (penalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }


    const Matrix DB = prod(r_D,B);
    Vector old_displacement = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        old_displacement[0] += N(0,i) * old_displacement_coefficient_vector[2*i];
        old_displacement[1] += N(0,i) * old_displacement_coefficient_vector[2*i + 1];
    }

    // Assembly
    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);

        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim;

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int jglob = 2*j+jdim;

                        // PENALTY TERM
                        rLeftHandSideMatrix(iglob, jglob) += N(0,i)*N(0,j)* penalty_integration * direction[idim] * direction[jdim];

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        Vector sigma_u_n = ZeroVector(3);
                        sigma_u_n[0] = DB(0, jglob)*normal_physical_space[0] + DB(2, jglob)*normal_physical_space[1];
                        sigma_u_n[1] = DB(2, jglob)*normal_physical_space[0] + DB(1, jglob)*normal_physical_space[1];

                        double sigma_u_n_dot_direction = inner_prod(sigma_u_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= N(0,i) * sigma_u_n_dot_direction * direction[idim] * integration_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        Vector sigma_w_n = ZeroVector(3);
                        sigma_w_n[0] = (mDBOld(0, iglob)* normal_physical_space[0] + mDBOld(2, iglob)* normal_physical_space[1]);
                        sigma_w_n[1] = (mDBOld(2, iglob)* normal_physical_space[0] + mDBOld(1, iglob)* normal_physical_space[1]);
                        

                        double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);
                        rLeftHandSideMatrix(iglob, jglob) -= nitsche_penalty*N(0,j) * sigma_w_n_dot_direction * direction[jdim] * integration_weight;
                    }

                }
            }
        }
    }
    else {
        // ASSIGN BC BY COMPONENTS 
        //--------------------------------------------------------------------------------------------
        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int id1 = 2*idim;
                    const int iglob = 2*i+idim;

                    // PENALTY TERM
                    rLeftHandSideMatrix(2*i+idim, 2*j+idim) += N(0,i)*N(0,j)* penalty_integration;

                    Vector sigma_w_n = ZeroVector(3);
                    sigma_w_n[0] = (mDBOld(0, iglob)* normal_physical_space[0] + mDBOld(2, iglob)* normal_physical_space[1]);
                    sigma_w_n[1] = (mDBOld(2, iglob)* normal_physical_space[0] + mDBOld(1, iglob)* normal_physical_space[1]);

                    #ifndef SWITCH_OFF_NISCHE

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int id2 = (id1+2)%3;
                        const int jglob = 2*j+jdim;

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= N(0,i)*(DB(id1, jglob)* normal_physical_space[0] + DB(id2, jglob)* normal_physical_space[1]) * integration_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= nitsche_penalty*N(0,j)*sigma_w_n[jdim] * integration_weight;
                    }

                    #endif

                }
            }
        }
    }
    
    KRATOS_CATCH("")
}


void SupportSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int mDim = r_DN_De[0].size2();
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);
    const double penalty = GetProperties()[PENALTY_FACTOR];

    KRATOS_ERROR_IF(mDim != 2) << "SolidElement momentarily only supports 2D elements, but the current element has dimension " << mDim << std::endl;

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // compute complete jacobian transformation including parameter->physical space transformation
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0, 0) = J0[0](0, 0);
    Jacobian(0, 1) = J0[0](0, 1);
    Jacobian(1, 0) = J0[0](1, 0);
    Jacobian(1, 1) = J0[0](1, 1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    // retrieve the normal of the physical space
    array_1d<double, 3> normal_physical_space = GetValue(NORMAL);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);
    noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);

    Matrix B = ZeroMatrix(mDim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);


    #ifdef ANALYTICAL_TANGENT
    Matrix r_D = ZeroMatrix(3,3);
    AnalyticalConstitutiveMatrix(r_D, old_strain);

    Vector r_stress_vector = ZeroVector(3);
    AnalyticalStress(r_stress_vector, old_strain);
    #else
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    const Vector& r_stress_vector = Values.GetStressVector();
    #endif

    // Differential area
    double penalty_integration = penalty * integration_weight;

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = 1.0;  // = 1 -> Penalty approach
                                        // = -1 -> Free-penalty approach
    if (penalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }

    // Assembly

    const Matrix DB = prod(r_D,B);
    Vector old_displacement = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        old_displacement[0] += N(0,i) * old_displacement_coefficient_vector[2*i];
        old_displacement[1] += N(0,i) * old_displacement_coefficient_vector[2*i + 1];
    }

    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);
        const Vector displacement = this->GetValue(DISPLACEMENT); //already times direction
        const double displacement_module = inner_prod(displacement, direction);

        const double old_displacement_direction = inner_prod(old_displacement, direction);
            
        for (IndexType i = 0; i < number_of_control_points; i++) {

            for (IndexType idim = 0; idim < 2; idim++) {
                const int iglob = 2*i+idim;

                rRightHandSideVector(iglob) += N(0,i) * direction[idim] * (displacement_module-old_displacement_direction) * penalty_integration;

                // // PENALTY FREE g_n = 0
                // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                // //*********************************************** */
                Vector sigma_w_n = ZeroVector(3);
                sigma_w_n[0] = (mDBOld(0, iglob)* normal_physical_space[0] + mDBOld(2, iglob)* normal_physical_space[1]);
                sigma_w_n[1] = (mDBOld(2, iglob)* normal_physical_space[0] + mDBOld(1, iglob)* normal_physical_space[1]);

                double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                //PENALTY FREE
                // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                rRightHandSideVector(iglob) -= nitsche_penalty*sigma_w_n_dot_direction * integration_weight *(displacement_module - old_displacement_direction);

                // residual terms

                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                old_stress_normal[0] = r_stress_vector[0]*normal_physical_space[0] + r_stress_vector[2]*normal_physical_space[1];
                old_stress_normal[1] = r_stress_vector[2]*normal_physical_space[0] + r_stress_vector[1]*normal_physical_space[1];

                double old_stress_normal_dot_direction = inner_prod(old_stress_normal, direction);
                rRightHandSideVector(iglob) += N(0,i) * old_stress_normal_dot_direction * direction[idim] * integration_weight;
            }
        }
    }
    else {
        // ASSIGN BC BY COMPONENTS 
        //--------------------------------------------------------------------------------------------

        Vector u_D = this->GetValue(DISPLACEMENT);

        for (IndexType i = 0; i < number_of_control_points; i++) {

            for (IndexType idim = 0; idim < 2; idim++) {
                const int iglob = 2*i+idim;

                rRightHandSideVector[iglob] += N(0,i)*(u_D-old_displacement)[idim]* penalty_integration;

                Vector sigma_w_n = ZeroVector(3);
                sigma_w_n[0] = (mDBOld(0, iglob)* normal_physical_space[0] + mDBOld(2, iglob)* normal_physical_space[1]);
                sigma_w_n[1] = (mDBOld(2, iglob)* normal_physical_space[0] + mDBOld(1, iglob)* normal_physical_space[1]);

                #ifndef SWITCH_OFF_NISCHE

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    rRightHandSideVector(iglob) -= nitsche_penalty*(u_D[jdim]-old_displacement[jdim])*sigma_w_n[jdim] * integration_weight;
                }

                // residual terms
                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                old_stress_normal[0] = r_stress_vector[0]*normal_physical_space[0] + r_stress_vector[2]*normal_physical_space[1];
                old_stress_normal[1] = r_stress_vector[2]*normal_physical_space[0] + r_stress_vector[1]*normal_physical_space[1];

                rRightHandSideVector(iglob) += N(0,i) * old_stress_normal[idim] * integration_weight;

                #endif

            }
        }
    }
    KRATOS_CATCH("")
}

    int SupportSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SupportSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void SupportSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void SupportSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void SupportSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
            rB(1, r) = r_DN_DX(kr,1) * dirr;
            rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
        }
    }


    void SupportSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
    {
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before applying the constitutive law." << std::endl;
        rValues.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);
        
        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
    }


    void SupportSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int mDim = 2;
        Matrix DN_DX(number_of_control_points,2);
        Matrix InvJ0(mDim,mDim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        Matrix delta_position;
        CalculateDeltaPositionMatrix(r_geometry, delta_position);
        r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);
        double DetJ0;
        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);
        
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);

        array_1d<double, 3> normal_physical_space = GetValue(NORMAL);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before finalizing the solution step." << std::endl;
        Values.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,old_displacement);
    
        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        // update internal variables before requesting stresses
        mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
        

        #ifdef ANALYTICAL_TANGENT
        Vector sigma = ZeroVector(3);
        AnalyticalStress(sigma, old_strain);
        #else
        const Vector sigma = Values.GetStressVector();
        #endif

        
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*normal_physical_space[0] + sigma[2]*normal_physical_space[1];
        sigma_n[1] = sigma[2]*normal_physical_space[0] + sigma[1]*normal_physical_space[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
    }

    void SupportSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
        //--------------------------------------------------------------------------------------------
        // calculate the constitutive law response
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        Matrix delta_position;
        CalculateDeltaPositionMatrix(r_geometry, delta_position);
        r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

        // Initialize DN_DX
        const unsigned int dimension = 2;
        Matrix DN_DX(number_of_control_points, 2);
        Matrix InvJ0(dimension, dimension);
        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        double detJ0;

        Matrix Jacobian = ZeroMatrix(2, 2);
        Jacobian(0, 0) = J0[0](0, 0);
        Jacobian(0, 1) = J0[0](0, 1);
        Jacobian(1, 0) = J0[0](1, 0);
        Jacobian(1, 1) = J0[0](1, 1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian, InvJ0, detJ0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0], InvJ0);

        Matrix B = ZeroMatrix(3, mat_size);
        CalculateB(B, DN_DX);

        // Prepare constitutive law parameters
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Flags& r_constitutive_law_options = Values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(displacement_coefficient_vector);
        Vector old_strain = prod(B, displacement_coefficient_vector);

        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before initializing the solution step." << std::endl;
        Values.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        // Store D*B from the current constitutive response for reuse in assembly
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        #ifdef ANALYTICAL_TANGENT
        Matrix r_D = ZeroMatrix(3,3);
        AnalyticalConstitutiveMatrix(r_D, old_strain);
        #else
        const Matrix& r_D = Values.GetConstitutiveMatrix();
        #endif


        if (mDBOld.size1() != r_D.size1() || mDBOld.size2() != B.size2()) {
            mDBOld.resize(r_D.size1(), B.size2(), false);
        }
        noalias(mDBOld) = prod(r_D, B);

        // update internal variables before requesting stresses
        mpConstitutiveLaw->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
    }

} // Namespace Kratos
