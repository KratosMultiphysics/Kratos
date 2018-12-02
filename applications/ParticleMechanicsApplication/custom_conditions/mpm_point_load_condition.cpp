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
#include "custom_conditions/mpm_point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMPointLoadCondition::MPMPointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMPointLoadCondition::MPMPointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMPointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<MPMPointLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMPointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<MPMPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMPointLoadCondition::~MPMPointLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void MPMPointLoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize )
            {
                rLeftHandSideMatrix.resize( MatSize, MatSize, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix(MatSize); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != MatSize )
            {
                rRightHandSideVector.resize( MatSize, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
        }

        // Vector with a loading applied to the condition
        array_1d<double, 3 > PointLoad = ZeroVector(3);
        if( this->Has( POINT_LOAD ) )
        {
            noalias(PointLoad) = this->GetValue( POINT_LOAD );
        }

        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            const unsigned int base = ii*Dimension;

            if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
            {
                noalias(PointLoad) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
            }

            for(unsigned int k = 0; k < Dimension; ++k)
            {
                rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * PointLoad[k];
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    double MPMPointLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

} // Namespace Kratos


