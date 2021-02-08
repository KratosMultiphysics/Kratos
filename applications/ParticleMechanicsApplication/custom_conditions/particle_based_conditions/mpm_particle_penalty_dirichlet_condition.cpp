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
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 2.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * m_unit_normal;
            r_geometry[i].UnSetLock();
        }
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
            if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            {
                for ( unsigned int j = 0; j < dimension; j++)
                {
                    field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                }
            }
        }

        const double penetration = MathUtils<double>::Dot((field_displacement - m_imposed_displacement), m_unit_normal);

        // If penetrates, apply constraint, otherwise no
        if (penetration >= 0.0)
        {
            apply_constraints = false;
        }
    }

    if (apply_constraints)
    {
        // Arrange shape function
        Matrix shape_function = ZeroMatrix(block_size, matrix_size);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon())
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
            {
                gap_function[block_size * i + j] = (Variables.CurrentDisp(i,j) - m_imposed_displacement[j]);
            }
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
    else{
        // To improve stability: use identity matrix to avoid nonzero diagonal LHS matrix
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix) = IdentityMatrix(matrix_size);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

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

    KRATOS_CATCH( "" )
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
    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticlePenaltyDirichletCondition::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
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
    std::vector<array_1d<double, 3 > > rValues,
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
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

} // Namespace Kratos


