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
#include "custom_conditions/grid_based_conditions/mpm_grid_point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMGridPointLoadCondition::MPMGridPointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMGridBaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMGridPointLoadCondition::MPMGridPointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMGridBaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMGridPointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MPMGridPointLoadCondition>(NewId, pGeometry, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMGridPointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMGridPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMGridPointLoadCondition::~MPMGridPointLoadCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMGridPointLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
        {
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != matrix_size )
        {
            rRightHandSideVector.resize( matrix_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
    }

    // Vector with a loading applied to the condition
    array_1d<double, 3 > point_load = ZeroVector(3);
    if( this->Has( POINT_LOAD ) )
    {
        noalias(point_load) = this->GetValue( POINT_LOAD );
    }

    for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
    {
        const unsigned int base = ii*dimension;

        if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
        {
            noalias(point_load) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
        }

        for(unsigned int k = 0; k < dimension; ++k)
        {
            rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * point_load[k];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double MPMGridPointLoadCondition::GetPointLoadIntegrationWeight()
{
    return 1.0;
}

} // Namespace Kratos


