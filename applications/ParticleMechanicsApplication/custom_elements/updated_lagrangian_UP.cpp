//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"

namespace Kratos
{

///**
//* Flags related to the element computation
//*/
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUP::UpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( UpdatedLagrangianUP const& rOther)
    :UpdatedLagrangian(rOther)
     //,mDeformationGradientF0(rOther.mDeformationGradientF0)
     //,mDeterminantF0(rOther.mDeterminantF0)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUP&  UpdatedLagrangianUP::operator=(UpdatedLagrangianUP const& rOther)
{
    UpdatedLagrangian::operator=(rOther);

    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUP::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianUP( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}
//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUP::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianUP NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new UpdatedLagrangianUP(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUP::~UpdatedLagrangianUP()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::Initialize()
{
    KRATOS_TRY

    // Initial position of the particle
    array_1d<double,3>& xg = this->GetValue(MP_COORD);

    // Initialize parameters
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(dimension);

    // Compute initial jacobian matrix and inverses
    Matrix J0 = ZeroMatrix(dimension);
    J0 = this->MPMJacobian(J0, xg);
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    // Compute current jacobian matrix and inverses
    Matrix j = ZeroMatrix(dimension);
    j = this->MPMJacobian(j,xg);
    double detj;
    MathUtils<double>::InvertMatrix( j, mInverseJ, detj );

    // Initialize constitutive law and materials
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void UpdatedLagrangianUP::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> xg = this->GetValue(MP_COORD);
    const array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
    const array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);

    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_Acceleration = ZeroVector(3);
    array_1d<double,3> MP_Velocity = ZeroVector(3);
    double MP_Pressure = 0.0;
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rVariables.N[i] > 1e-16)
        {
            array_1d<double, 3 > & nodal_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);

            double nodal_pressure = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 0);
            MP_Pressure += rVariables.N[i] * nodal_pressure;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * nodal_acceleration[j];

                /* NOTE: The following interpolation techniques have been tried:
                    MP_Velocity[j]      += rVariables.N[i] * nodal_velocity[j];
                    MP_Acceleration[j]  += nodal_inertia[j]/(rVariables.N[i] * MP_Mass * MP_number);
                    MP_Velocity[j]      += nodal_momentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                    MP_Velocity[j]      += delta_time * rVariables.N[i] * nodal_acceleration[j];
                */
            }
        }

    }

    /* NOTE:
    Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
    This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
    MP_Velocity = MP_PreviousVelocity + 0.5 * delta_time * (MP_Acceleration + MP_PreviousAcceleration);
    this -> SetValue(MP_VELOCITY,MP_Velocity );

    /* NOTE: The following interpolation techniques have been tried:
        MP_Acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
        MP_Velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    */

    // Update the MP Pressure
    this -> SetValue(MP_PRESSURE,MP_Pressure);

    // Update the MP Position
    const array_1d<double,3>& new_xg = xg + delta_xg ;
    this -> SetValue(MP_COORD,new_xg);

    //Update the MP Acceleration
    this -> SetValue(MP_ACCELERATION,MP_Acceleration);

    // Update the MP total displacement
    array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);
    MP_Displacement += delta_xg;
    this -> SetValue(MP_DISPLACEMENT,MP_Displacement);

    KRATOS_CATCH( "" )
}
//************************************************************************************
//*****************check size of LHS and RHS matrices*********************************

void UpdatedLagrangianUP::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing the LHS matrix if needed
    unsigned int matrix_size = number_of_nodes * dimension + number_of_nodes; // number of DOF including pressure term

    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size); //resetting LHS
    }

    // Resizing the RHS vector if needed
    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != matrix_size )
            rRightHandSideVector.resize( matrix_size, false );

        rRightHandSideVector = ZeroVector( matrix_size ); //resetting RHS
    }
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianUP::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);

    // Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx

    /* NOTE::
    Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] is to be updated in constitutive law parameter as total deformation gradient.
    The increment of total deformation gradient can be evaluated in 2 ways, which are:
    1. By: noalias( rVariables.F ) = prod( rVariables.j, InvJ);
    2. By means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    (NOTICE: Here, we are using method no. 1)
    */

    // Update Deformation gradient
    noalias( rVariables.F ) = prod( rVariables.j, InvJ);

    // Determinant of the previous Deformation Gradient F_n
    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);

    KRATOS_CATCH( "" )
}
//************************************************************************************

void UpdatedLagrangianUP::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); // Set all components to zero

    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
        }
    }
    else if( dimension == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );
        }
    }
    else
    {
        KRATOS_ERROR << "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateDeformationMatrix" << std::endl;
    }

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianUP::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    array_1d<double,3>& xg = this->GetValue(MP_COORD);
    GeneralVariables Variables;

    // Calculating and storing inverse and the determinant of the jacobian
    Matrix J0 = ZeroMatrix(dimension);
    J0 = this->MPMJacobian(J0, xg);
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

    // Initialize Constitutive Law
    mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
            GetGeometry(), Variables.N, rCurrentProcessInfo );

    mFinalizedStep = false;

    array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
    array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
    double MP_Pressure = this->GetValue(MP_PRESSURE);
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    double AUX_MP_Pressure = this->GetValue(AUX_MP_PRESSURE);
    double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> MP_Momentum;
    array_1d<double,3> MP_Inertia;
    array_1d<double,3> nodal_momentum;
    array_1d<double,3> nodal_inertia;

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > & nodal_acceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
        array_1d<double, 3 > & nodal_velocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

        // These are the values of nodal pressure evaluated in the initialize solution step
        double & nodal_pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE,1);

        AUX_MP_Pressure += Variables.N[j] * nodal_pressure;

        for (unsigned int k = 0; k < dimension; k++)
        {
            AUX_MP_Velocity[k] += Variables.N[j] * nodal_velocity[k];
            AUX_MP_Acceleration[k] += Variables.N[j] * nodal_acceleration[k];
        }
    }

    // Here MP contribution in terms of momentum, inertia, mass-pressure and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double nodal_mpressure =  Variables.N[i] * (MP_Pressure - AUX_MP_Pressure) * MP_Mass;

        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
            nodal_inertia[j]  = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;
        }

        GetGeometry()[i].SetLock();
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0)  += nodal_momentum;
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_INERTIA, 0)   += nodal_inertia;
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MPRESSURE, 0) += nodal_mpressure;

        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
        GetGeometry()[i].UnSetLock();
    }

    AUX_MP_Velocity.clear();
    AUX_MP_Acceleration.clear();
    AUX_MP_Pressure = 0.0;

}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    // Contribution of the internal and external forces
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

    rVariables.detF0   *= rVariables.detF;
    double determinant_F = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

    rVariables.detF     = determinant_F;
    rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianUP::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index_up = dimension * i + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index_up + j] += rVariables.N[i] * rVolumeForce[j];
        }
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    VectorType internal_forces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_up = dimension * i + i;
        unsigned int index_u  = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index_up + j] -= internal_forces[index_u + j];
        }
    }

    KRATOS_CATCH( "" )
}

//******************************************************************************************************************
//******************************************************************************************************************
double& UpdatedLagrangianUP::CalculatePUCoefficient(double& rCoefficient, GeneralVariables & rVariables)
{
    KRATOS_TRY

    // TODO: Check what is the meaning of this function
    rCoefficient = rVariables.detF0 - 1;

    return rCoefficient;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianUP::CalculatePUDeltaCoefficient(double &rDeltaCoefficient, GeneralVariables & rVariables)
{

    KRATOS_TRY

    // TODO: Check what is the meaning of this function
    rDeltaCoefficient = 1.0;

    return rDeltaCoefficient;

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int index_p = dimension;

    // FIXME: This is only for Solid Mechanics Problem with young_modulus Modulus and Poisson Ratio
    // TODO: Think about a more general way to find Bulk Modulus
    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio    = GetProperties()[POISSON_RATIO];
    double bulk_modulus  = young_modulus/(3.0*(1.0-2.0*poisson_ratio));

    // Check if Bulk Modulus is not NaN
    if (bulk_modulus != bulk_modulus)
        bulk_modulus = 1.e16;

    double delta_coefficient = 0;
    delta_coefficient = this->CalculatePUDeltaCoefficient( delta_coefficient, rVariables );

    double coefficient = 0;
    coefficient = this->CalculatePUCoefficient( coefficient, rVariables );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            double& pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);

            // TODO: Check what is the meaning of this equation
            rRightHandSideVector[index_p] += (1.0/(delta_coefficient * bulk_modulus)) * rVariables.N[i] * rVariables.N[j] * pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D-3D
        }

        rRightHandSideVector[index_p] -=  coefficient/delta_coefficient * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int index_p = dimension;

    double delta_coefficient = 0;
    delta_coefficient = this->CalculatePUDeltaCoefficient( delta_coefficient, rVariables );
    VectorType Fh=rRightHandSideVector;

    // Stabilization alpha parameters
    double alpha_stabilization  = 1.0;
    double stabilization_factor = 1.0;
    if( GetProperties().Has(STABILIZATION_FACTOR) ){
        stabilization_factor = GetProperties()[STABILIZATION_FACTOR];
    }
    alpha_stabilization *= stabilization_factor;

    // FIXME: This is only for Solid Mechanics Problem with young_modulus Modulus and Poisson Ratio
    // TODO: Think about a more general stabilization term if it is for Fluid Mechanics Problem
    const double& young_modulus          = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio    = GetProperties()[POISSON_RATIO];
    const double lame_mu =  young_modulus/(2.0*(1.0+poisson_ratio));

    double consistent = 1;
    double factor_value = 8.0; //JMR deffault value
    if( dimension == 3 )
        factor_value = 10.0; //JMC deffault value

    // TODO: Check what is the meaning of this equation
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            double& pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);

            if( dimension == 2 )
            {
                consistent=(-1)*alpha_stabilization*factor_value/(36.0*lame_mu);
                if(i==j)
                    consistent=2*alpha_stabilization*factor_value/(36.0*lame_mu);


                rRightHandSideVector[index_p] += consistent * pressure * rIntegrationWeight / (delta_coefficient * (rVariables.detF0/rVariables.detF)); //2D
            }
            else
            {
                consistent=(-1)*alpha_stabilization*factor_value/(80.0*lame_mu);
                if(i==j)
                    consistent=3*alpha_stabilization*factor_value/(80.0*lame_mu);

                rRightHandSideVector[index_p] += consistent * pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{
    // Contributions of the stiffness matrix calculated on the reference configuration
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    rVariables.detF0   *= rVariables.detF;
    double determinant_F = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: add Km to the rLefsHandSideMatrix
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kg to the rLefsHandSideMatrix
    CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp_Stab to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    rVariables.detF     = determinant_F;
    rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight
                                             )
{
    KRATOS_TRY

    Matrix Kuum = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexi = 0;
    unsigned int indexj = 0;

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuum(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const int size = number_of_nodes * dimension;

    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix reduced_Kg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( stress_tensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
    Matrix Kuug = ZeroMatrix(size);
    MathUtils<double>::ExpandAndAddReducedMatrix( Kuug, reduced_Kg, dimension );

    // Assemble components considering added DOF matrix system
    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuug(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_p  = dimension;
        unsigned int index_up = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_up+k,index_p) +=  rVariables.DN_DX ( i, k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
            }
            index_p += (dimension + 1);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Assemble components considering added DOF matrix system
    unsigned int index_p = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int index_up = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_p,index_up+k) +=  rVariables.N[i] * rVariables.DN_DX ( j, k ) * rIntegrationWeight * rVariables.detF;
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // FIXME: This is only for Solid Mechanics Problem with young_modulus Modulus and Poisson Ratio
    // TODO: Think about a more general way to find Bulk Modulus
    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio    = GetProperties()[POISSON_RATIO];
    double bulk_modulus  = young_modulus/(3.0*(1.0-2.0*poisson_ratio));

    // Check if Bulk Modulus is not NaN
    if (bulk_modulus != bulk_modulus)
        bulk_modulus = 1.e16;

    double delta_coefficient = 0;
    delta_coefficient = this->CalculatePUDeltaCoefficient( delta_coefficient, rVariables );

    unsigned int indexpi = dimension;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(bulk_modulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight /(delta_coefficient * (rVariables.detF0/rVariables.detF));

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
// I changed the constant matrix in the stabilized term:
// as in MPM the position of the integration points does not coincide with the
// position of the Gauss points the first matrix is substitute with the product of the
// shape function values of each integration point
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double delta_coefficient = 0;
    delta_coefficient = this->CalculatePUDeltaCoefficient( delta_coefficient, rVariables );

    unsigned int indexpi = dimension;

    // Stabilization alpha parameters
    double alpha_stabilization  = 1.0;
    double stabilization_factor = 1.0;
    if( GetProperties().Has(STABILIZATION_FACTOR) ){
        stabilization_factor = GetProperties()[STABILIZATION_FACTOR];
    }
    alpha_stabilization *= stabilization_factor;

    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio = GetProperties()[POISSON_RATIO];
    const double lame_mu        =  young_modulus/(2.0*(1.0+poisson_ratio));

    double consistent = 1.0;

    double factor_value = 8.0; //JMR deffault value
    if( dimension == 3 )
        factor_value = 10.0; //JMC deffault value

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            if( dimension == 2 )  //consistent 2D
            {
                consistent=(-1)*alpha_stabilization*factor_value/(36.0*lame_mu);
                if(indexpi==indexpj)
                    consistent=2*alpha_stabilization*factor_value/(36.0*lame_mu);

                rLeftHandSideMatrix(indexpi,indexpj) -= consistent *rIntegrationWeight / (delta_coefficient * (rVariables.detF0/rVariables.detF)); //2D
            }
            else
            {
                consistent=(-1)*alpha_stabilization*factor_value/(80.0*lame_mu);
                if(indexpi==indexpj)
                    consistent=3*alpha_stabilization*factor_value/(80.0*lame_mu);

                rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D
            }
            indexpj += (dimension + 1);
        }
        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianUP::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
    }
}


//************************************************************************************
//****************MASS MATRIX*********************************************************

void UpdatedLagrangianUP::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Call the values of the shape function for the single element
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int matrix_size = number_of_nodes * dimension + number_of_nodes;

    if ( rMassMatrix.size1() != matrix_size )
        rMassMatrix.resize( matrix_size, matrix_size, false );

    rMassMatrix = ZeroMatrix(matrix_size);

    double TotalMass = 0;

    // TOTAL MASS OF ONE MP ELEMENT
    TotalMass = this->GetValue(MP_MASS);

    // LUMPED MATRIX
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = Variables.N[i] * TotalMass;
        unsigned int index_up = i * dimension + i;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rMassMatrix( index_up+j, index_up+j ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
            values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue( PRESSURE, Step );
        }
        else
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( PRESSURE, Step );
        }

    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    UpdatedLagrangian::GetHistoricalVariables(rVariables);
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double voigtsize = 3;
    if ( dimension == 3)
        voigtsize = 6;

    UpdatedLagrangian::FinalizeStepVariables( rVariables, rCurrentProcessInfo);

    // Evaluation of the pressure on the material point
    double nodal_mean_stress = 0.0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
        nodal_mean_stress += GetGeometry()[i].FastGetSolutionStepValue( PRESSURE ) * rVariables.N[i];

    // Evaluation of the mean stress on the material point
    double mean_stress = 0.0;
    for (unsigned int i = 0; i < dimension; i++)
        mean_stress += rVariables.StressVector[i];
    mean_stress /= dimension;

    Vector stress_vector = ZeroVector(voigtsize);
    stress_vector = rVariables.StressVector;
    for (unsigned int i = 0; i < dimension; i++)
        stress_vector[i] += (nodal_mean_stress - mean_stress);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, stress_vector);

}

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int UpdatedLagrangianUP::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    int correct = 0;

    correct = UpdatedLagrangian::Check(rCurrentProcessInfo);

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW)) << "Constitutive law is not compatible with the U-P element type: Large Displacements U_P" << std::endl;

    // Verify that the variables are correctly initialized
    KRATOS_ERROR_IF( PRESSURE.Key() == 0 ) <<  "PRESSURE has Key zero! (check if the application is correctly registered" << std::endl;

    return correct;

    KRATOS_CATCH( "" );
}

void UpdatedLagrangianUP::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);


}

void UpdatedLagrangianUP::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}

} // Namespace Kratos

