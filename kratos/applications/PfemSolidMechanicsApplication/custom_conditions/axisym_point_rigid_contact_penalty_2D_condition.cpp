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
#include "custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  AxisymPointRigidContactPenalty2DCondition::AxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
										       pGeometry)
  : PointRigidContactPenalty2DCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  AxisymPointRigidContactPenalty2DCondition::AxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : PointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  AxisymPointRigidContactPenalty2DCondition::AxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : PointRigidContactPenalty2DCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  AxisymPointRigidContactPenalty2DCondition::AxisymPointRigidContactPenalty2DCondition( AxisymPointRigidContactPenalty2DCondition const& rOther )
  : PointRigidContactPenalty2DCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer AxisymPointRigidContactPenalty2DCondition::Create(IndexType NewId, NodesArrayType
								       const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new AxisymPointRigidContactPenalty2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************************************************************
  //************************************************************************************


  AxisymPointRigidContactPenalty2DCondition::~AxisymPointRigidContactPenalty2DCondition()
  {
  }

   

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************
  
  void AxisymPointRigidContactPenalty2DCondition::CalculateKinematics(GeneralVariables& rVariables,
								      const double& rPointNumber)
  {
    KRATOS_TRY

    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);
    
    PointRigidContactPenalty2DCondition::CalculateKinematics(rVariables, rPointNumber);
    
  
    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateRadius(double & rCurrentRadius, double & rReferenceRadius)
  {

    KRATOS_TRY

      rCurrentRadius=0;
    rReferenceRadius=0;

    //Displacement from the reference to the current configuration
    array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
    array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
    array_1d<double, 3 > & ReferencePosition    = GetGeometry()[0].Coordinates();
    array_1d<double, 3 > CurrentPosition        = ReferencePosition + DeltaDisplacement;

    rCurrentRadius   = CurrentPosition[0];
    rReferenceRadius = ReferencePosition[0];

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
  {

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contributions to stiffness matrix calculated on the reference config

    PointRigidContactCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH(rLeftHandSideMatrix)
  }


  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contribution to external forces

    PointRigidContactCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH(rRightHandSideVector)
  }

  
  





} // Namespace Kratos



