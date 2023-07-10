// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PointLoadCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PointLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PointLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

PointLoadCondition::~PointLoadCondition()
{
}

//************************************************************************************
//************************************************************************************

void PointLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
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

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
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

double PointLoadCondition::GetPointLoadIntegrationWeight() const
{
    return 1.0;
}

} // Namespace Kratos


