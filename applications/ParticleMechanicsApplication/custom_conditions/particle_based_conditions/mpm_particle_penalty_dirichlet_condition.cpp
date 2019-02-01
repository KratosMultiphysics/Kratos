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
    const double penalty_factor = 10.0e9; //TODO: make it as an input

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    // Arrange shape function
    Matrix shape_function = ZeroMatrix(dimension, matrix_size);
    for (unsigned int i = 0; i < number_of_nodes; i++)
        if (Variables.N[i] > std::numeric_limits<double>::epsilon())
            for (unsigned int j = 0; j < dimension; j++)
                shape_function(j, dimension * i + j) = Variables.N[i];

    Vector gap_function = ZeroVector(matrix_size);
    for (unsigned int i = 0; i < number_of_nodes; i++)
        for ( unsigned int j = 0; j < dimension; j++ )
            gap_function[dimension * i + j] = (Variables.CurrentDisp(i,j) - imposed_displacement[j]);

    noalias(rLeftHandSideMatrix)  += prod(trans(shape_function), shape_function);
    noalias(rRightHandSideVector) -= prod(prod(trans(shape_function), shape_function), gap_function);

    rLeftHandSideMatrix  *= penalty_factor;
    rRightHandSideVector *= penalty_factor;

    // // Calculating field displacement
    // array_1d<double, 3 > field_displacement = ZeroVector(3);
    // for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // {
    //     if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
    //     {
    //         for ( unsigned int j = 0; j < dimension; j++ )
    //             field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
    //     }
    // }

    // // Calculate gap function / penetration
    // const Vector gap_function = field_displacement - imposed_displacement;

    // // Prepare variables
    // const double penalty_factor = 10.0e5; //TODO: make it as an input

    // // Compute modified shape function vector: modified_N (number_of_nodes * dimension)
    // Vector modified_N = ZeroVector(matrix_size);
    // Vector modified_N_gap = ZeroVector(matrix_size);
    // for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
    // {
    //     const unsigned int base = ii * block_size;
    //     for(unsigned int k = 0; k < dimension; ++k)
    //     {
    //         // modified_N[base + k] = Variables.N[ii] * normal_vector[k];
    //         modified_N[base+k] = Variables.N[ii];
    //         modified_N_gap[base+k] = Variables.N[ii] * gap_function[k];
    //     }
    // }

    // // Compute Right Hand Side Vector
    // for(unsigned int i = 0; i < matrix_size; ++i)
    // {
    //     rRightHandSideVector[i] = - penalty_factor * modified_N_gap[i];
    // }

    // // Compute Left Hand Side Matrix
    // for (unsigned int i = 0; i<matrix_size; ++i)
    // {
    //     for (unsigned int j = 0; j<matrix_size; ++j)
    //     {
    //         rLeftHandSideMatrix(i,j) = penalty_factor *  modified_N[i] * modified_N[j];
    //     }
    // }
    // if (this->Id() == 3){
    //     std::cout << "TEST::ID" << this->Id() << std::endl;
    //     std::cout << "TEST::position" << xg_c << std::endl;
    //     std::cout << "TEST::Shape function" << Variables.N << std::endl;
    //     std::cout << "TEST::gap_function" << gap_function << std::endl;
    //     std::cout << "TEST::modified_N" << modified_N << std::endl;
    //     std::cout << "TEST::modified_N_gap" << modified_N_gap << std::endl;
    //     std::cout << "TEST::rRightHandSideVector" << rRightHandSideVector << std::endl;
    //     std::cout << "TEST::rLeftHandSideMatrix" << rLeftHandSideMatrix << std::endl;
    // }

    KRATOS_CATCH( "" )
}

} // Namespace Kratos


