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
#include "custom_utilities/mpm_nonconforming_boundary_rotation_utility.h"
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
    const array_1d<double, 3 > & imposed_displacement = this->GetValue (MPC_IMPOSED_DISPLACEMENT);

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
        array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
        ParticleMechanicsMathUtilities<double>::Normalize(unit_normal_vector);
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
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix)  += prod(trans(shape_function), shape_function);
            rLeftHandSideMatrix  *= penalty_factor * this->GetIntegrationWeight();
        }

        if ( CalculateResidualVectorFlag == true )
        {
            noalias(rRightHandSideVector) -= prod(prod(trans(shape_function), shape_function), gap_function);
            rRightHandSideVector *= penalty_factor * this->GetIntegrationWeight();
        }
    }
    else{
        // To improve stability: use identity matrix to avoid nonzero diagonal LHS matrix
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix) = IdentityMatrix(matrix_size);
        }
    }

    // If SLIP
    if (Is(SLIP)){
        MPMNonconformingBoundaryRotationUtility<MatrixType, VectorType> rotation_tool(dimension, block_size);
        if ( CalculateStiffnessMatrixFlag == true )
        {
            rotation_tool.Rotate(rLeftHandSideMatrix, rRightHandSideVector, this);
            rotation_tool.ApplySlipCondition( rLeftHandSideMatrix, rRightHandSideVector, this);
            rotation_tool.Recover( rLeftHandSideMatrix, rRightHandSideVector, this);
            //FIXME: REMOVE THE FOLLOWING LINES
            if (this->Id() == 1009) std::cout << "rLeftHandSideMatrix" << rLeftHandSideMatrix << std::endl;
            if (this->Id() == 1009) std::cout << "rRightHandSideVector" << rRightHandSideVector << std::endl;
        }
        else{
            rotation_tool.Rotate( rRightHandSideVector, this);
            rotation_tool.ApplySlipCondition( rRightHandSideVector, this);
            rotation_tool.Recover(rRightHandSideVector, this);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int MPMParticlePenaltyDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

double MPMParticlePenaltyDirichletCondition::GetIntegrationWeight()
{
    double integration_weight = MPMParticleBaseDirichletCondition::GetIntegrationWeight();

    // Check for point support
    if (integration_weight < std::numeric_limits<double>::epsilon())
        integration_weight = 1.0;

    return integration_weight;
}

} // Namespace Kratos


