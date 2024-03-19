//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"
#include "includes/kratos_flags.h"



namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::MPMParticlePenaltyDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::MPMParticlePenaltyDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MPMParticlePenaltyDirichletCondition>(NewId, pGeometry, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMParticlePenaltyDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::~MPMParticlePenaltyDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        GeneralVariables Variables;

        // Calculating shape function
        MPMShapeFunctionPointValues(Variables.N);

        // Here MPC contribution of normal vector are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Set(SLIP);
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 2.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * m_normal;
            r_geometry[i].UnSetLock();
        }
    }
}

void MPMParticlePenaltyDirichletCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    // At the beginning of NonLinearIteration, REACTION has to be reset to zero
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(REACTION).clear();
        r_geometry[i].UnSetLock();
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();
    const GeometryType& r_geometry = GetGeometry();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
        {
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != matrix_size )
        {
            rRightHandSideVector.resize( matrix_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
    }

    // Prepare variables
    GeneralVariables Variables;

    // Calculating shape function
    MPMShapeFunctionPointValues(Variables.N);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    // Check contact: Check contact penetration: if <0 apply constraint, otherwise no
    bool apply_constraints = true;
    if (Is(CONTACT))
    {
        // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
        array_1d<double, 3 > field_displacement = ZeroVector(3);
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for ( unsigned int j = 0; j < dimension; j++)
                field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
        }

        const double penetration = MathUtils<double>::Dot((field_displacement - m_imposed_displacement), m_normal);

        // If penetrates, apply constraint, otherwise no
        if (penetration >= 0.0)
            apply_constraints = false;

    }

    if (apply_constraints)
    {
        // Arrange shape function
        Matrix shape_function = ZeroMatrix(block_size, matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            // constrain only the movement of the nodes which are conntected to the body
            if (r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) >= std::numeric_limits<double>::epsilon() )
            {
                for (unsigned int j = 0; j < dimension; j++)
                {
                    shape_function(j, block_size * i + j) = Variables.N[i];
                }
            }
        }

        // Calculate gap_function: nodal_displacement - imposed_displacement
        Vector gap_function = ZeroVector(matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            for ( unsigned int j = 0; j < dimension; j++)
                gap_function[block_size * i + j] = (Variables.CurrentDisp(i,j) - m_imposed_displacement[j]);

        }

        // Calculate LHS Matrix and RHS Vector
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix)  += prod(trans(shape_function), shape_function);
            rLeftHandSideMatrix  *= m_penalty * this->GetIntegrationWeight();
        }

        if ( CalculateResidualVectorFlag == true )
        {
            noalias(rRightHandSideVector) -= prod(prod(trans(shape_function), shape_function), gap_function);
            rRightHandSideVector *= m_penalty * this->GetIntegrationWeight();
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Recalculate resiudal vector for converged solution
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int block_size = this->GetBlockSize();

    const unsigned int matrix_size = number_of_nodes * block_size;
    VectorType rRightHandSideVector = ZeroVector( matrix_size );

    MPMParticlePenaltyDirichletCondition::CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    // Calculate nodal forces
    Vector nodal_force = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        for (unsigned int j = 0; j < dimension; j++)
            nodal_force[j] = rRightHandSideVector[block_size * i + j];

        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(REACTION) += nodal_force;
        r_geometry[i].UnSetLock();


    }

    KRATOS_CATCH( "" )
}

void MPMParticlePenaltyDirichletCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();

        // Here MPC normal vector and IS_STRUCTURE are reset
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Reset(SLIP);
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL).clear();
            r_geometry[i].UnSetLock();
        }
    }
    this->CalculateInterfaceContactForce( rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

void MPMParticlePenaltyDirichletCondition::CalculateInterfaceContactForce(const ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    // Prepare variables
    GeneralVariables Variables;
    const double & r_mpc_area = this->GetIntegrationWeight();
    MPMShapeFunctionPointValues(Variables.N);

    // Interpolate the force to mpc_force assuming linear shape function
    array_1d<double, 3 > mpc_force = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double nodal_area  = 0.0;
        if (r_geometry[i].SolutionStepsDataHas(NODAL_AREA))
            nodal_area= r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0);

        const Vector nodal_force = r_geometry[i].FastGetSolutionStepValue(REACTION);

        if (nodal_area > std::numeric_limits<double>::epsilon())
        {
            mpc_force += Variables.N[i] * nodal_force * r_mpc_area / nodal_area;
        }
    }

    // Apply in the normal contact direction and allow releasing motion
    if (Is(CONTACT))
    {
        // Apply only in the normal direction
        const double normal_force = MathUtils<double>::Dot(mpc_force, m_normal);

        // This check is done to avoid sticking forces
        if (normal_force > 0.0)
            mpc_force = -1.0 * normal_force * m_normal;
        else
            mpc_force = ZeroVector(3);
    }
    // Apply a sticking contact
    else{
        mpc_force *= -1.0;
    }

    // Set Contact Force
    m_contact_force = mpc_force;
}

void MPMParticlePenaltyDirichletCondition::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PENALTY_FACTOR) {
        rValues[0] = m_penalty;
    }
    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticlePenaltyDirichletCondition::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == PENALTY_FACTOR) {
        m_penalty = rValues[0];
    }
    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}


} // Namespace Kratos


