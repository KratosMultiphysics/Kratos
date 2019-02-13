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
#include "utilities/math_utils.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

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

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<MPMParticlePenaltyDirichletCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<MPMParticlePenaltyDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticlePenaltyDirichletCondition::~MPMParticlePenaltyDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
        GeneralVariables Variables;

        // Calculating shape function
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

        // Normal Vector
        const array_1d<double,3> & unit_normal_vector = this->GetValue(MPC_NORMAL);

        // Here MPC contribution of normal vector are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(IS_STRUCTURE) = 2.0;
            rGeom[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * unit_normal_vector;
            rGeom[i].UnSetLock();
        }
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
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

    // Get imposed displacement and normal vector
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    const array_1d<double, 3 > & imposed_displacement = this->GetValue (MPC_DISPLACEMENT);

    // Prepare variables
    GeneralVariables Variables;
    const double penalty_factor = this->GetValue(PENALTY_FACTOR);

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    // Check contact: Check contact penetration: if <0 apply constraint, otherwise no
    bool apply_constraints = true;
    if (Is(CONTACT))
    {
        // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
        const array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
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

        const double penetration = MathUtils<double>::Dot((field_displacement - imposed_displacement), unit_normal_vector);

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
                gap_function[block_size * i + j] = (Variables.CurrentDisp(i,j) - imposed_displacement[j]);
            }
        }

        // Calculate LHS Matrix and RHS Vector
        noalias(rLeftHandSideMatrix)  += prod(trans(shape_function), shape_function);
        noalias(rRightHandSideVector) -= prod(prod(trans(shape_function), shape_function), gap_function);

        rLeftHandSideMatrix  *= penalty_factor * this->GetIntegrationWeight();
        rRightHandSideVector *= penalty_factor * this->GetIntegrationWeight();
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();

        // Here MPC normal vector and IS_STRUCTURE are reset
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
            rGeom[i].FastGetSolutionStepValue(NORMAL).clear();
            rGeom[i].UnSetLock();
        }
    }

    // Additional treatment for interface conditions
    if (Is(INTERFACE))
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();

        // Prepare variables
        GeneralVariables Variables;
        const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

        double MPC_Equivalent_Mass = 0.0;
        array_1d<double,3> MPC_Acceleration = ZeroVector(3);
        array_1d<double,3> MPC_Force = ZeroVector(3);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon())
            {
                const double & nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS);
                const array_1d<double, 3 > & nodal_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);

                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    MPC_Acceleration[j] += Variables.N[i] * nodal_acceleration[j];
                }

                if (nodal_mass > std::numeric_limits<double>::epsilon())
                    MPC_Equivalent_Mass += Variables.N[i] * Variables.N[i] / nodal_mass;
            }
        }

        if (MPC_Equivalent_Mass > std::numeric_limits<double>::epsilon())
            MPC_Equivalent_Mass = 1.0 / MPC_Equivalent_Mass;

        MPC_Force = MPC_Equivalent_Mass * MPC_Acceleration;

        this->SetValue(MPC_CONTACT_FORCE, MPC_Force);

    }

    KRATOS_CATCH( "" )
}

} // Namespace Kratos


