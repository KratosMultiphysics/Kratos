//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

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
    return Kratos::make_shared<PointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }


   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer PointRigidContactPenalty2DCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
   {
     return Kratos::make_shared<PointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
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

    if ( this->GetProperties().Has(THICKNESS) ) {
       if ( dimension == 2 && GetProperties()[THICKNESS]>0 )
          rIntegrationWeight *= GetProperties()[THICKNESS];

    }
    return rIntegrationWeight;
  }



} // Namespace Kratos
