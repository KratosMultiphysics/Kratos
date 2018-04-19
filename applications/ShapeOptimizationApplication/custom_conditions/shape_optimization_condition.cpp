// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_conditions/shape_optimization_condition.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{

// ------------------------------------------------------------------------------
ShapeOptimizationCondition::ShapeOptimizationCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

// ------------------------------------------------------------------------------
ShapeOptimizationCondition::ShapeOptimizationCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

// ------------------------------------------------------------------------------
ShapeOptimizationCondition::ShapeOptimizationCondition( ShapeOptimizationCondition const& rOther )
    : Condition(rOther)
{
}

// ------------------------------------------------------------------------------
Condition::Pointer ShapeOptimizationCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ShapeOptimizationCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// ------------------------------------------------------------------------------
Condition::Pointer ShapeOptimizationCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}

// ------------------------------------------------------------------------------
ShapeOptimizationCondition::~ShapeOptimizationCondition()
{
}

// ------------------------------------------------------------------------------
int ShapeOptimizationCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

// ------------------------------------------------------------------------------
void ShapeOptimizationCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

// ------------------------------------------------------------------------------
void ShapeOptimizationCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}

// ------------------------------------------------------------------------------
} // Namespace Kratos.
