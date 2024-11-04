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
#include "custom_conditions/support_pressure_dirichlet_condition.h"

namespace Kratos
{
void SupportPressureDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Integration
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    const unsigned int dim = DN_De[0].size2();
    Matrix DN_DX(number_of_nodes,dim);

    const SizeType mat_size = number_of_nodes * (dim+1);
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    // Compute the normals
    array_1d<double, 3> normal_parameter_space;
    r_geometry.Calculate(NORMAL, normal_parameter_space);

    // Matrix DB = prod(r_D,B);
    const Matrix& N = r_geometry.ShapeFunctionsValues();
    // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = DN_De[0]; // prod(DN_De[point_number],InvJ0);

    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // constitutive law
    Matrix B = ZeroMatrix(3,number_of_nodes*dim);
    CalculateB(B, DN_DX);
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(number_of_nodes*dim);
    GetValuesVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    Values.SetStrainVector(old_strain);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);
    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix sigmaVoigt = Matrix(prod(r_D, B));
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Vector n_tensor(2);
    n_tensor(0) = normal_parameter_space(0); // Component in x direction
    n_tensor(1) = normal_parameter_space(1); // Component in y direction
    for (IndexType i = 0; i < number_of_nodes; i++) {
        
        for (IndexType j = 0; j < number_of_nodes; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
            
                // integration by parts velocity is moved to the RHS
                // rLeftHandSideMatrix(3*i+idim, 3*j+idim) -= N(0,i)*(
                //         DN_DX(j, 0) * normal_parameter_space[0] + DN_DX(j, 1) * normal_parameter_space[1] )
                //         * integration_points[0].Weight();
                //// CONSTITUTIVE LAW
                // for (IndexType jdim = 0; jdim < 2; jdim++) {
                //     // Extract the 2x2 block for the control point i from the sigma matrix.
                //     Matrix sigma_block = ZeroMatrix(2, 2);
                //     sigma_block(0, 0) = sigmaVoigt(0, 2*j+jdim);      // sigma(4 * j + 2*jdim, 0);     // sigma_xx for control point i.
                //     sigma_block(0, 1) = sigmaVoigt(2, 2*j+jdim);      // sigma(4 * j + 2*jdim, 1);     // sigma_xy for control point i.
                //     sigma_block(1, 0) = sigmaVoigt(2, 2*j+jdim);      // sigma(4 * j + 2*jdim + 1, 0); // sigma_yx (symmetric) for control point i.
                //     sigma_block(1, 1) = sigmaVoigt(1, 2*j+jdim);      // sigma(4 * j + 2*jdim + 1, 1); // sigma_yy for control point i.
                //     // Compute the traction vector: sigma * n.
                //     Vector traction = prod(sigma_block, n_tensor); // This results in a 2x1 vector.
                //     rLeftHandSideMatrix(3*i+idim, 3*j+jdim) -= N(0, i) * traction(idim) * integration_points[0].Weight();
                // }

                // // integration by parts PRESSURE is moved to the RHS
                // rLeftHandSideMatrix(3*i+idim, 3*j+2) += N(0,j)* ( N(0,i) * normal_parameter_space[idim] )
                //         * integration_points[0].Weight();

            }
        }
    }
    
        
    double p_D = this->GetValue(PRESSURE);

    for (IndexType j = 0; j < number_of_nodes; j++) {
        
        //// integration by parts PRESSURE < v dot n , p >
        for (IndexType idim = 0; idim < 2; idim++) {
            rRightHandSideVector(3*j+idim) -= p_D * ( N(0,j) * normal_parameter_space[idim] ) * integration_points[0].Weight();
        }
        
        // Neumann condition for the velocity
        Vector t_N = ZeroVector(2); 
        double x = r_geometry.Center().X(); 
        double y = r_geometry.Center().Y();
        const double pi = 3.14159265358979323846;
        // t_N[0] = ( sin(2*pi * x)*2*pi*sin(2*pi * y)                   )*normal_parameter_space[0]  +
        //          (-cos(2*pi * x)*2*pi*cos(2*pi * y)+cos(2*pi * y)*2*pi)*normal_parameter_space[1]; 
        // t_N[1] = ( cos(2*pi * x)*2*pi*cos(2*pi * y)-cos(2*pi * x)*2*pi)*normal_parameter_space[0]  +
        //          (-sin(2*pi * x)*2*pi*sin(2*pi * y)                   )*normal_parameter_space[1]; 

        // LINEAR -> laplacian
        // t_N[0] = (1.0)*normal_parameter_space[0] + (0)*normal_parameter_space[1]; 
        // t_N[1] = (0)*normal_parameter_space[0] + (-1.0)*normal_parameter_space[1]; 

        // LINEAR -> constitutive law
        // t_N[0] = (2.0)*normal_parameter_space[0] + (0)*normal_parameter_space[1]; 
        // t_N[1] = (0)*normal_parameter_space[0] + (-2.0)*normal_parameter_space[1];  

        // QUADRATIC -> laplacian
        // t_N[0] = (2*x)*normal_parameter_space[0] + (-2*y)*normal_parameter_space[1]; 
        // t_N[1] = (-2*y)*normal_parameter_space[0] + (-2*x)*normal_parameter_space[1];

        // // QUADRATIC -> constitutive law
        // t_N[0] = (4*x)*normal_parameter_space[0] + (-4*y)*normal_parameter_space[1]; 
        // t_N[1] = (-4*y)*normal_parameter_space[0] + (-4*x)*normal_parameter_space[1]; 

        // CUBIC -> constitutive law
        // t_N[0] = (6*x*x)*normal_parameter_space[0] + (-3*y*y-6*x*y)*normal_parameter_space[1]; 
        // t_N[1] = (-3*y*y-6*x*y)*normal_parameter_space[0] + (-6*x*x)*normal_parameter_space[1]; 
        
        t_N[0] = 0;
        t_N[1] = 0; 

        for (IndexType idim = 0; idim < 2; idim++) {
            rRightHandSideVector(3*j+idim) += N(0,j) * t_N[idim] * integration_points[0].Weight();
        }
        
    }

    // Add residual of previous iteration to RHS
    VectorType temp = ZeroVector(number_of_nodes*3);
    // RHS = ExtForces - K*temp;
    unsigned int index = 0 ;
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        temp[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_X);
        temp[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_Y);
        temp[index++] = r_geometry[i].GetSolutionStepValue(PRESSURE);
    }
    // // RHS = ExtForces - K*temp;
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    // // RHS = ExtForces - IntForces;
    // for (IndexType i = 0; i < number_of_nodes; ++i) {
    //     r_stress_vector[0] -= r_geometry[i].GetSolutionStepValue(PRESSURE)*N(0,i); // σ_xx
    //     r_stress_vector[1] -= r_geometry[i].GetSolutionStepValue(PRESSURE)*N(0,i); // σ_yy
    //     r_stress_vector[2] -= 0.0*N(0,i);                                          // σ_xy
    // }
    // Vector internalForces = integration_points[0].Weight() * prod(trans(B), r_stress_vector);
    // for (IndexType i = 0; i < number_of_nodes; ++i)
    // {
    //     // Add only to the velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
    //     for (IndexType idim = 0; idim < 2; ++idim)
    //     {
    //         // viscous + pressure term
    //         rRightHandSideVector(i * 3 + idim) -= internalForces(i * 2 + idim);
    //     }  
    // }

    // for (unsigned int i = 0; i < number_of_nodes; i++) {
    //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
    //     outputFile << r_geometry[i].GetId() << "  " << r_geometry[i].GetDof(TEMPERATURE).EquationId() <<"\n";
    //     outputFile.close();
    // }
    KRATOS_CATCH("")
}

int SupportPressureDirichletCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
    //     << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPressureDirichletCondition" << std::endl;
    return 0;
}

void SupportPressureDirichletCondition::CalculateB(
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

void SupportPressureDirichletCondition:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void SupportPressureDirichletCondition::InitializeMaterial()
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

void SupportPressureDirichletCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const int dim = 2;
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (dim + 1) * number_of_control_points;

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


void SupportPressureDirichletCondition::GetDofList(
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


void SupportPressureDirichletCondition::GetValuesVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 2 >& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}


} // Namespace Kratos
