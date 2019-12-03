//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ilaria Iaconeta
//                  Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianElement::Initialize()
{
    KRATOS_TRY

    // Initialize parameters
    const SizeType working_space_dimension = GetGeometry().WorkingSpaceDimension();

    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(working_space_dimension);

    // Initialize constitutive law and materials
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::InitializeMaterial()
{
    KRATOS_TRY

    mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();

    mpConstitutiveLaw->InitializeMaterial(
        GetProperties(),
        GetGeometry(),
        row(GetGeometry().ShapeFunctionsValues(), 0));

    KRATOS_CATCH("")
}

void UpdatedLagrangianElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mpConstitutiveLaw->ResetMaterial(
        GetProperties(),
        GetGeometry(),
        row(GetGeometry().ShapeFunctionsValues(), 0));

    KRATOS_CATCH("")
}

void UpdatedLagrangianElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    mFinalizedStep = false;

    const array_1d<double, 3>& MP_velocity = this->GetValue(MP_VELOCITY);
    const array_1d<double, 3>& MP_acceleration = this->GetValue(MP_ACCELERATION);
    const double & MP_mass = this->GetValue(MP_MASS);

    array_1d<double, 3> aux_MP_velocity = ZeroVector(3);
    array_1d<double, 3> aux_MP_acceleration = ZeroVector(3);
    array_1d<double, 3> nodal_momentum = ZeroVector(3);
    array_1d<double, 3> nodal_inertia = ZeroVector(3);

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for (unsigned int j = 0; j < number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(ACCELERATION))
            nodal_acceleration = r_geometry[j].FastGetSolutionStepValue(ACCELERATION, 1);

        array_1d<double, 3 > nodal_velocity = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(VELOCITY))
            nodal_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY, 1);

        for (unsigned int k = 0; k < dimension; k++)
        {
            aux_MP_velocity[k] += r_N(0, j) * nodal_velocity[k];
            aux_MP_acceleration[k] += r_N(0, j) * nodal_acceleration[k];
        }
    }

    // Here MP contribution in terms of momentum, inertia and mass are added
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = r_N(0, i) * (MP_velocity[j] - aux_MP_velocity[j]) * MP_mass;
            nodal_inertia[j] = r_N(0, i) * (MP_acceleration[j] - aux_MP_acceleration[j]) * MP_mass;

        }

        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
        r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0) += nodal_inertia;

        r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) += r_N(0, i) * MP_mass;
        r_geometry[i].UnSetLock();

    }
}

void UpdatedLagrangianElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    KinematicVariables kinematic_variables(WorkingSpaceDimension(), GetGeometry().size());

    CalculateKinematics(kinematic_variables, rCurrentProcessInfo);

    ConstitutiveVariables constitutive_variables(WorkingSpaceDimension());
    CalculateConstitutiveVariables(
        kinematic_variables,
        constitutive_variables,
        constitutive_law_parameters,
        ConstitutiveLaw::StressMeasure_Cauchy);

    double integration_weight = CalculateIntegrationWeight(kinematic_variables);

    // Compute the deformation matrix B
    Matrix B;
    this->CalculateBMatrix(B, kinematic_variables.DN_DX);

    if (CalculateStiffnessMatrixFlag)
    {
        // Operation performed: add K_material to the rLefsHandSideMatrix
        this->CalculateAndAddKuum(rLeftHandSideMatrix, B, constitutive_variables.ConstitutiveMatrix, integration_weight);

        // Operation performed: add K_geometry to the rLefsHandSideMatrix
        this->CalculateAndAddKuug(rLeftHandSideMatrix, kinematic_variables.DN_DX, constitutive_variables.StressVector, integration_weight);
    }
    if (CalculateResidualVectorFlag)
    {
        Vector volume_force;
        this->CalculateVolumeForce(volume_force);
        // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
        this->CalculateAndAddExternalForces(rRightHandSideVector, volume_force, integration_weight);

        // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
        this->CalculateAndAddInternalForces(rRightHandSideVector, B, constitutive_variables.StressVector, integration_weight);
    }
}

void UpdatedLagrangianElement::CalculateConstitutiveVariables(
    KinematicVariables& rKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
) const
{
    // Set element specific options
    rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    Matrix FT = prod(rKinematicVariables.F, mDeformationGradientF0);

    rValues.SetDeterminantF(rKinematicVariables.detFT);
    rValues.SetDeformationGradientF(FT);

    rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);

    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix);

    rValues.SetShapeFunctionsValues(row(GetGeometry().ShapeFunctionsValues(), 0));

    rThisConstitutiveVariables.ConstitutiveMatrix(0, 0) = 1000;
    rThisConstitutiveVariables.ConstitutiveMatrix(1, 1) = 1000;
    rThisConstitutiveVariables.ConstitutiveMatrix(2, 2) = 500;

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ThisStressMeasure);
}

void UpdatedLagrangianElement::CalculateKinematics(
    KinematicVariables& rKinematicVariables,
    ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    Matrix Jacobian;
    Jacobian = GetGeometry().Jacobian( Jacobian, 0);

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix( Jacobian, InvJ, detJ);

    Matrix current_displacement;
    SetCurrentDisplacement(current_displacement);

    KRATOS_WATCH(current_displacement)

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
    Matrix jacobian;
    jacobian = GetGeometry().Jacobian(
        jacobian, 0, GetGeometry().GetDefaultIntegrationMethod(), current_displacement);

    // Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( jacobian, Invj, detJ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n+1]
    rKinematicVariables.DN_DX = prod(
        GetGeometry().ShapeFunctionLocalGradient(0), Invj); //overwrites DX now is the current position dx

    /* NOTE::
    Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] is to be updated in constitutive law parameter as total deformation gradient.
    The increment of total deformation gradient can be evaluated in 2 ways, which are:
    1. By: noalias( rVariables.F ) = prod( jacobian, InvJ);
    2. By means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    (NOTICE: Here, we are using method no. 2)*/

    // METHOD 1: Update Deformation gradient: F [dx_n+1/dx_n] = [dx_n+1/d£] [d£/dx_n]
    // noalias( rVariables.F ) = prod( jacobian, InvJ);

    // METHOD 2: Update Deformation gradient: F_ij = δ_ij + u_i,j
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    Matrix I = IdentityMatrix(dimension);
    Matrix gradient_displacement = ZeroMatrix(dimension, dimension);
    gradient_displacement = prod(trans(current_displacement), rKinematicVariables.DN_DX);

    rKinematicVariables.F = (I + gradient_displacement);

    rKinematicVariables.detF = MathUtils<double>::Det(rKinematicVariables.F);

    rKinematicVariables.detFT = rKinematicVariables.detF * mDeterminantF0;

    KRATOS_CATCH( "" )
}
//************************************************************************************

void UpdatedLagrangianElement::CalculateBMatrix(
    Matrix& rB, const Matrix& rDN_DX) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {
        if (rB.size1() != 3 && rB.size2() != number_of_nodes * dimension)
            rB.resize(3, number_of_nodes * dimension);
        rB.clear(); // Set all components to zero

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            IndexType index = 2 * i;
            rB( 0, index     ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index     ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
        }
    }
    else if( dimension == 3 )
    {
        if (rB.size1() != 6 && rB.size2() != number_of_nodes * dimension)
            rB.resize(6, number_of_nodes * dimension);
        rB.clear(); // Set all components to zero

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            IndexType index = 3 * i;

            rB( 0, index     ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index     ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index     ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );
        }
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::CalculateAndAddExternalForces(
    VectorType& rRightHandSideVector,
    Vector& rVolumeForce,
    const double& rIntegrationWeight) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType working_space_dimension = GetGeometry().WorkingSpaceDimension();

    const Matrix r_N = GetGeometry().ShapeFunctionsValues();

    for (IndexType integration_point = 0; integration_point < r_N.size1(); ++integration_point)
    {
        for (IndexType i = 0; i < number_of_nodes; i++)
        {
            int index = working_space_dimension * i;

            for (IndexType j = 0; j < working_space_dimension; j++)
            {
                rRightHandSideVector[index + j] += r_N(integration_point, i) * rVolumeForce[j];
            }
        }
    }

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        const Matrix& rB,
        const Vector& rStressVector,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType internal_forces = rIntegrationWeight * prod( trans(rB), rStressVector);
    noalias( rRightHandSideVector ) -= internal_forces;

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::CalculateAndAddKuum(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rB,
    const Matrix& rConstitutiveMatrix,
    const double& rIntegrationWeight) const
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += prod( trans(rB), rIntegrationWeight * Matrix( prod(rConstitutiveMatrix, rB ) ) );

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::CalculateAndAddKuug(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rDN_DX,
    const Vector& rStressVector,
    const double& rIntegrationWeight) const
{
    KRATOS_TRY

    const SizeType working_space_dimension = GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(rStressVector);
    Matrix reduced_Kg = prod(rDN_DX, rIntegrationWeight * Matrix( prod( stress_tensor, trans(rDN_DX) ) ) );
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, reduced_Kg, working_space_dimension);

    KRATOS_CATCH( "" )
}

double& UpdatedLagrangianElement::CalculateVolumeChange( double& rVolumeChange,
    KinematicVariables & rVariables )
{
    rVolumeChange = 1.0 / (rVariables.detF * mDeterminantF0);

    return rVolumeChange;
}

Vector& UpdatedLagrangianElement::CalculateVolumeForce(
    Vector& rVolumeForce)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION) * this->GetValue(MP_MASS);

    return rVolumeForce;

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Create and initialize element variables:
    KinematicVariables kinematic_variables(WorkingSpaceDimension(), GetGeometry().size());

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters constitutive_law_parameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Compute element kinematics F, DN_DX ...
    this->CalculateKinematics(kinematic_variables, rCurrentProcessInfo);

    // Create and initialize element variables:
    ConstitutiveVariables constitutive_variables(WorkingSpaceDimension());

    // Compute material response
    this->CalculateConstitutiveVariables(
        kinematic_variables,
        constitutive_variables,
        constitutive_law_parameters,
        ConstitutiveLaw::StressMeasure_Cauchy);

    // Call the constitutive law to update material variables
    mpConstitutiveLaw->FinalizeMaterialResponse(
        constitutive_law_parameters,
        ConstitutiveLaw::StressMeasure_Cauchy);

    // Call the element internal variables update
    this->FinalizeStepVariables(
        kinematic_variables,
        constitutive_variables,
        rCurrentProcessInfo);

    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}


////************************************************************************************
////************************************************************************************

void UpdatedLagrangianElement::FinalizeStepVariables(
    KinematicVariables& rVariables,
    ConstitutiveVariables& rConstitutiveVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Update internal (historical) variables
    mDeterminantF0         = rVariables.detF* mDeterminantF0;
    mDeformationGradientF0 = prod(rVariables.F, mDeformationGradientF0);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, rConstitutiveVariables.StressVector);
    this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rConstitutiveVariables.StrainVector);

    // Delta Plastic Strains
    double delta_plastic_strain = mpConstitutiveLaw->GetValue(MP_DELTA_PLASTIC_STRAIN, delta_plastic_strain );
    this->SetValue(MP_DELTA_PLASTIC_STRAIN, delta_plastic_strain);

    double delta_plastic_volumetric_strain = mpConstitutiveLaw->GetValue(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN, delta_plastic_volumetric_strain);
    this->SetValue(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN, delta_plastic_volumetric_strain);

    double delta_plastic_deviatoric_strain = mpConstitutiveLaw->GetValue(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN, delta_plastic_deviatoric_strain);
    this->SetValue(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN, delta_plastic_deviatoric_strain);

    // Total Plastic Strain
    double equivalent_plastic_strain = mpConstitutiveLaw->GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, equivalent_plastic_strain );
    this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, equivalent_plastic_strain);

    double accumulated_plastic_volumetric_strain = mpConstitutiveLaw->GetValue(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN, accumulated_plastic_volumetric_strain);
    this->SetValue(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN, accumulated_plastic_volumetric_strain);

    double accumulated_plastic_deviatoric_strain = mpConstitutiveLaw->GetValue(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN, accumulated_plastic_deviatoric_strain);
    this->SetValue(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN, accumulated_plastic_deviatoric_strain);

    this->UpdateGaussPoint(rCurrentProcessInfo);

}

/// The position of the Gauss points/Material points is updated
void UpdatedLagrangianElement::UpdateGaussPoint(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    Matrix current_displacement;
    SetCurrentDisplacement(current_displacement);

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const array_1d<double,3> & xg = this->GetValue(MP_COORD);
    const array_1d<double,3> & MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
    const array_1d<double,3> & MP_PreviousVelocity = this->GetValue(MP_VELOCITY);

    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_acceleration = ZeroVector(3);
    array_1d<double,3> MP_velocity = ZeroVector(3);
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    const Matrix& N = GetGeometry().ShapeFunctionsValues();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (N(0, i) > std::numeric_limits<double>::epsilon())
        {
            auto r_geometry = GetGeometry();
            array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
            if (r_geometry[i].SolutionStepsDataHas(ACCELERATION))
                nodal_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION);

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                delta_xg[j] += N(0, i) * current_displacement(i,j);
                MP_acceleration[j] += N(0, i) * nodal_acceleration[j];

                /* NOTE: The following interpolation techniques have been tried:
                    MP_velocity[j]      += rVariables.N[i] * nodal_velocity[j];
                    MP_acceleration[j]  += nodal_inertia[j]/(rVariables.N[i] * MP_mass * MP_number);
                    MP_velocity[j]      += nodal_momentum[j]/(rVariables.N[i] * MP_mass * MP_number);
                    MP_velocity[j]      += delta_time * rVariables.N[i] * nodal_acceleration[j];
                */
            }
        }

    }

    /* NOTE:
    Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
    This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
    MP_velocity = MP_PreviousVelocity + 0.5 * delta_time * (MP_acceleration + MP_PreviousAcceleration);
    this -> SetValue(MP_VELOCITY,MP_velocity );

    /* NOTE: The following interpolation techniques have been tried:
        MP_acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
        MP_velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    */

    // Update the MP Position
    const array_1d<double,3>& new_xg = xg + delta_xg ;
    this -> SetValue(MP_COORD, new_xg);

    // Update the MP Acceleration
    this -> SetValue(MP_ACCELERATION, MP_acceleration);

    // Update the MP total displacement
    array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);
    MP_Displacement += delta_xg;
    this -> SetValue(MP_DISPLACEMENT,MP_Displacement);

    KRATOS_CATCH( "" )
}

/// This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
Matrix& UpdatedLagrangianElement::SetCurrentDisplacement(
    Matrix & rCurrentDisp) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    //resize Matrix
    if (rCurrentDisp.size1() != number_of_nodes && rCurrentDisp.size2() != dimension)
        rCurrentDisp.resize(number_of_nodes, dimension);
    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for (IndexType i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

        for (IndexType j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}

double UpdatedLagrangianElement::CalculateIntegrationWeight(const KinematicVariables& rKinematicVariables)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    /* NOTE:
    The material points will have constant mass as defined at the beginning.
    However, the density and volume (integration weight) are changing every time step.*/
    // Update MP_Density
    const double MP_density = (GetProperties()[DENSITY]) / rKinematicVariables.detFT;
    this->SetValue(MP_DENSITY, MP_density);

    // The MP_Volume (integration weight) is evaluated
    double MP_volume = this->GetValue(MP_MASS) / this->GetValue(MP_DENSITY);
    this->SetValue(MP_VOLUME, MP_volume);

    //if( dimension == 2 )
    //    MP_volume *= GetProperties()[THICKNESS];

    return MP_volume;
}

void UpdatedLagrangianElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    int number_of_nodes = r_geometry.size();
    int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        IndexType index = i * dimension;
        rResult[index] = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

void UpdatedLagrangianElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& CurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );

        if ( r_geometry.WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

}

void UpdatedLagrangianElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != matrix_size )
        rDampingMatrix.resize( matrix_size, matrix_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix(matrix_size, matrix_size);

    //1.-Calculate StiffnessMatrix:
    MatrixType StiffnessMatrix  = Matrix();
    this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

    //2.-Calculate MassMatrix:
    MatrixType MassMatrix  = Matrix();
    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

    //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
    {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
    {
        beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    //4.-Compose the Damping Matrix:
    //Rayleigh Damping Matrix: alpha*M + beta*K
    rDampingMatrix  = alpha * MassMatrix;
    rDampingMatrix += beta  * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // Lumped
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int matrix_size = dimension * number_of_nodes;

    const Matrix& N = GetGeometry().ShapeFunctionsValues();

    if ( rMassMatrix.size1() != matrix_size )
        rMassMatrix.resize( matrix_size, matrix_size, false );

    rMassMatrix = ZeroMatrix(matrix_size, matrix_size);

    // TOTAL MASS OF ONE MP ELEMENT
    const double & r_total_mass = this->GetValue(MP_MASS);

    // LUMPED MATRIX
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = N(0, i) * r_total_mass;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

void UpdatedLagrangianElement::GetValuesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType working_space_dimension = r_geometry.WorkingSpaceDimension();
    const SizeType matrix_size = number_of_nodes * working_space_dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( IndexType i = 0; i < number_of_nodes; i++ )
    {
        IndexType index = i * working_space_dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if (working_space_dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

void UpdatedLagrangianElement::GetFirstDerivativesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType working_space_dimension = r_geometry.WorkingSpaceDimension();
    const SizeType matrix_size = number_of_nodes * working_space_dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for (IndexType i = 0; i < number_of_nodes; i++ )
    {
        IndexType index = i * working_space_dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

        if (working_space_dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
    }
}

void UpdatedLagrangianElement::GetSecondDerivativesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType working_space_dimension = r_geometry.WorkingSpaceDimension();
    const SizeType matrix_size = number_of_nodes * working_space_dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( IndexType i = 0; i < number_of_nodes; i++ )
    {
        IndexType index = i * working_space_dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if (working_space_dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

int  UpdatedLagrangianElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    Element::Check(rCurrentProcessInfo);

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
            correct_strain_measure = true;
    }

    KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)

    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    if( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
            KRATOS_CHECK_VARIABLE_KEY(THICKNESS)
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has( THICKNESS )) << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
        }

        // Check constitutive law
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), r_geometry, rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos

