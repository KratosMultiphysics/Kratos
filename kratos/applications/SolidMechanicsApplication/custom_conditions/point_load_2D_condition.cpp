//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/point_load_2D_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//***********************************************************************************
//***********************************************************************************
PointLoad2DCondition::PointLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointLoad3DCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
PointLoad2DCondition::PointLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : PointLoad3DCondition(NewId, pGeometry, pProperties)
{

    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointLoad2DCondition::PointLoad2DCondition( PointLoad2DCondition const& rOther )
    : PointLoad3DCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer PointLoad2DCondition::Create(IndexType NewId,
						NodesArrayType const& ThisNodes,
						PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoad2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer PointLoad2DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  PointLoad2DCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());
  
  //-----------//      
  return Condition::Pointer( new PointLoad2DCondition(NewCondition) );
}


//***********************************************************************************
//***********************************************************************************
PointLoad2DCondition::~PointLoad2DCondition()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


int PointLoad2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void PointLoad2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoad3DCondition )
}

void PointLoad2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoad3DCondition )
}


} // Namespace Kratos.
