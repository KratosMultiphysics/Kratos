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
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "includes/checks.h"

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
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 2.0; // flag for penalty-based slip cond
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * m_unit_normal;
            r_geometry[i].FastGetSolutionStepValue(NORMAL_REACTION, 1) += Variables.N[i] * m_normal_reaction;
            r_geometry[i].UnSetLock();
        }
    }
}

void MPMParticlePenaltyDirichletCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    AccumulateReactionToNodes(rCurrentProcessInfo);
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

    const unsigned int number_of_nodes = GetGeometry().size(); // # nodes in containing element
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize(); // DoF assoc with each node
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

        const double penetration = MathUtils<double>::Dot((field_displacement - m_imposed_displacement), m_unit_normal);

        // If penetrates, apply constraint, otherwise no
        if (penetration >= 0.0)
            apply_constraints = false;
        
    }

    if (apply_constraints)
    {
        // Matrix H (eq 44)
        // Arrange shape function
        Matrix shape_function = ZeroMatrix(block_size, matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
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
            // note: integration weight (area assoc. with boundary particle is set on its creation)
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
    AccumulateReactionToNodes(rCurrentProcessInfo);
}

void MPMParticlePenaltyDirichletCondition::AccumulateReactionToNodes(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Recalculate residual vector for converged solution
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

        // Check whether there nodes are active and associated to material point elements
        const double& nodal_mass = r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0);
        if (nodal_mass > std::numeric_limits<double>::epsilon())
        {
            r_geometry[i].SetLock();
            r_geometry[i].FastGetSolutionStepValue(REACTION) += nodal_force;
            r_geometry[i].UnSetLock();
        }

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

        // Prepare variables
        GeneralVariables Variables;
        const double & r_mpc_area = this->GetIntegrationWeight();
        MPMShapeFunctionPointValues(Variables.N);

        // Reset normal reaction
        m_normal_reaction = 0.0;

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            // Here MPC normal vector and IS_STRUCTURE are reset
            r_geometry[i].SetLock();
            r_geometry[i].Reset(SLIP);
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL).clear();
            r_geometry[i].UnSetLock();

            // Interpolate the nodal normal reaction to mpc assuming linear shape function

            // Check whether there is material point inside the node
            const double& nodal_mass = r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0);
            double nodal_area  = 0.0;
            if (r_geometry[i].SolutionStepsDataHas(NODAL_AREA))
                nodal_area= r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0);

            const double nodal_normal_reaction = r_geometry[i].FastGetSolutionStepValue(NORMAL_REACTION);

            if (nodal_mass > std::numeric_limits<double>::epsilon() && nodal_area > std::numeric_limits<double>::epsilon())
            {
                m_normal_reaction += Variables.N[i] * nodal_normal_reaction * r_mpc_area / nodal_area;
            }
        }
    }

    KRATOS_CATCH( "" )
}

void MPMParticlePenaltyDirichletCondition::CalculateInterfaceContactForce(array_1d<double, 3 >& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();  
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    // Prepare variables
    GeneralVariables Variables;
    const double & r_mpc_area = this->GetIntegrationWeight();
    MPMShapeFunctionPointValues(Variables.N);

    // Interpolate the force to mpc_force assuming linear shape function
    rVariable = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        // Check whether there is material point inside the node
        const double& nodal_mass = r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0);
        double nodal_area  = 0.0;        
        if (r_geometry[i].SolutionStepsDataHas(NODAL_AREA))
            nodal_area= r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0);

        const Vector nodal_force = r_geometry[i].FastGetSolutionStepValue(REACTION);

        if (nodal_mass > std::numeric_limits<double>::epsilon() && nodal_area > std::numeric_limits<double>::epsilon())
        {
            rVariable += Variables.N[i] * nodal_force * r_mpc_area / nodal_area;
        }
    }

    // Apply in the normal contact direction and allow releasing motion
    if (Is(CONTACT))
    {
        // Apply only in the normal direction
        const double normal_force = MathUtils<double>::Dot(rVariable, m_unit_normal);

        // This check is done to avoid sticking forces
        if (normal_force > 0.0)
            rVariable = -1.0 * normal_force * m_unit_normal;
        else
            rVariable = ZeroVector(3);
    }
    // Apply a sticking contact
    else{
        rVariable *= -1.0;
    }
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

void MPMParticlePenaltyDirichletCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        rValues[0] = m_imposed_displacement;
    }
    else if (rVariable == MPC_NORMAL) {
        rValues[0] = m_unit_normal;
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        this->CalculateInterfaceContactForce(rValues[0], rCurrentProcessInfo);
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

void MPMParticlePenaltyDirichletCondition::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        m_imposed_displacement = rValues[0];
    }
    else if (rVariable == MPC_NORMAL) {
        m_unit_normal = rValues[0];
        ParticleMechanicsMathUtilities<double>::Normalize(m_unit_normal);
    }
    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

int MPMParticlePenaltyDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points()){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,r_node)
    }

    return 0;
}

} // Namespace Kratos


