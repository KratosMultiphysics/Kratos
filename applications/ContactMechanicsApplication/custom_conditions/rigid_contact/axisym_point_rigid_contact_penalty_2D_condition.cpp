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
#include "custom_conditions/rigid_contact/axisym_point_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

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
    return Kratos::make_shared<AxisymPointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }

  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer AxisymPointRigidContactPenalty2DCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
  {
    return Kratos::make_shared<AxisymPointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
  }

  //************************************************************************************
  //************************************************************************************


  AxisymPointRigidContactPenalty2DCondition::~AxisymPointRigidContactPenalty2DCondition()
  {
  }



  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void AxisymPointRigidContactPenalty2DCondition::CalculateKinematics(ConditionVariables& rVariables,
								      const ProcessInfo& rCurrentProcessInfo,
								      const double& rPointNumber)
  {
    KRATOS_TRY

    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius);

    PointRigidContactPenalty2DCondition::CalculateKinematics(rVariables, rCurrentProcessInfo, rPointNumber);

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************


  void AxisymPointRigidContactPenalty2DCondition::CalculateContactFactors(ConditionVariables &rVariables)
  {

    KRATOS_TRY

    WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter  = 0;

    // double radius   = 0;
    // double meanradius = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){

	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();

	  // radius = fabs(Contact_Point[0] + rN[i].X()) * 0.5;
	  // meanradius += radius;
	  // distance += norm_2(Contact_Point-Neighb_Point) * radius;

	  distance += norm_2(Contact_Point-Neighb_Point);

	  counter ++;
	}
      }

    // if( Contact_Point[0] > 0)
    //   distance /= ( counter * Contact_Point[0] );
    // else
    //   distance /= ( counter * (meanradius/counter) );

    if( counter != 0 )
      distance /= counter;

    if( distance == 0 )
      distance = 1;

    //get contact properties and parameters
    double PenaltyParameter = 1;
    if( GetProperties().Has(PENALTY_PARAMETER) )
      PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

    WeakPointerVector<Element >& rE = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
    double ElasticModulus = 0;
    if( GetProperties().Has(YOUNG_MODULUS) ){
      ElasticModulus = GetProperties()[YOUNG_MODULUS];
    }
    else if( GetProperties().Has(C10) ){
	ElasticModulus = GetProperties()[C10];
    }
    else{

	if( rE.front().GetProperties().Has(YOUNG_MODULUS) ){
	    ElasticModulus = rE.front().GetProperties()[YOUNG_MODULUS];
	}
	else if( rE.front().GetProperties().Has(C10) ){
	    ElasticModulus = rE.front().GetProperties()[C10];
	}
    }

    // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
    if (ElasticModulus <= 1.0e-5) {
      std::vector<double> mModulus;
      ProcessInfo SomeProcessInfo;
      for ( unsigned int i = 0; i < rE.size(); i++)
	{
	  rE[i].CalculateOnIntegrationPoints(EQUIVALENT_YOUNG_MODULUS, mModulus, SomeProcessInfo);
	  ElasticModulus += mModulus[0];
	}
      ElasticModulus /= double(rE.size());
    }

    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;

    double PenaltyRatio = 1;
    if( GetProperties().Has(TANGENTIAL_PENALTY_RATIO) )
      PenaltyRatio = GetProperties()[TANGENTIAL_PENALTY_RATIO];

    rVariables.Penalty.Tangent = rVariables.Penalty.Normal * PenaltyRatio;


    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;

    //set contact normal
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	GetGeometry()[i].SetLock();

	array_1d<double, 3> &ContactNormal  = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_NORMAL);

	for(unsigned int i=0; i<3; i++)
	  ContactNormal[i] = rVariables.Surface.Normal[i];

	GetGeometry()[i].UnSetLock();
      }

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

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    if ( this->GetProperties().Has(THICKNESS) ) {
       if( GetProperties()[THICKNESS] > 0 )
          IntegrationWeight /=  GetProperties()[THICKNESS];
    }

    //contributions to stiffness matrix calculated on the reference config

    PointRigidContactCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
  }


  //************************************************************************************
  //************************************************************************************

  void  AxisymPointRigidContactPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;

    if ( this->GetProperties().Has(THICKNESS) ) {
       if( GetProperties()[THICKNESS] > 0 )
          IntegrationWeight /=  GetProperties()[THICKNESS];
    }

    //contribution to external forces

    PointRigidContactCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rRightHandSideVector )
  }








} // Namespace Kratos
