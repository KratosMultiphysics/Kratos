//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_rigid_contact_penalty_2D_condition.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : PointRigidContactPenalty3DCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : PointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : PointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition( PointRigidContactPenalty2DCondition const& rOther )
  : PointRigidContactPenalty3DCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer PointRigidContactPenalty2DCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new PointRigidContactPenalty2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************************************************************
  //************************************************************************************


  PointRigidContactPenalty2DCondition::~PointRigidContactPenalty2DCondition()
  {

  }



  
  //***********************************************************************************
  //************************************************************************************

  double& PointRigidContactPenalty2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
  { 
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( dimension == 2 && GetProperties()[THICKNESS]>0 ) 
      rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
  }



} // Namespace Kratos



