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


  void AxisymPointRigidContactPenalty2DCondition::CalculateContactFactors(GeneralVariables &rVariables)
  {

    KRATOS_TRY
      
    WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double radius   = 0;
    double distance = 0;
    double counter  = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){
	    
	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();
	    
	  radius = fabs(Contact_Point[0] + rN[i].X()) * 0.5;
	  
	  distance += norm_2(Contact_Point-Neighb_Point) * radius;

	  counter ++;
	}
      }

    distance /= counter;
    
    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;  
    

    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;
    
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
    array_1d<double, 3 > & CurrentPosition      = GetGeometry()[0].Coordinates();
    array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;


    rCurrentRadius   = CurrentPosition[0];
    rReferenceRadius = ReferencePosition[0];

    if( rCurrentRadius == 0 ){
      
      WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);
      
      double counter = 0;

      for(unsigned int i = 0; i < rN.size(); i++)
	{
	  array_1d<double, 3 > & NodePosition = rN[i].Coordinates();
	  if( NodePosition[0] != 0 ){
	    rCurrentRadius += NodePosition[0] * 0.225; 	    
	    counter ++;
	  }
	    
	}

      rCurrentRadius /= counter;
      
    }


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
  {

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contributions to stiffness matrix calculated on the reference config

    PointRigidContactCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
  }


  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contribution to external forces

    PointRigidContactCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rRightHandSideVector )
  }

  
  





} // Namespace Kratos



