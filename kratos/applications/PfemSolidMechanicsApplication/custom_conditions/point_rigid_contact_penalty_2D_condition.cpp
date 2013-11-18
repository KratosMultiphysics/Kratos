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

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : PointRigidContactCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : PointRigidContactCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : PointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition( PointRigidContactPenalty2DCondition const& rOther )
  : PointRigidContactCondition(rOther)
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

   
  //************************************************************************************
  //************************************************************************************

  void PointRigidContactPenalty2DCondition::InitializeGeneralVariables (GeneralVariables& rVariables, 
									const ProcessInfo& rCurrentProcessInfo)
  {



  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************
  
  void PointRigidContactPenalty2DCondition::CalculateKinematics(GeneralVariables& rVariables,
								const double& rPointNumber)
  {
    KRATOS_TRY
    
      int ContactFace = 0; //free surface
      
    if( this->mpRigidWall->IsInside( GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, ContactFace ) ){

      rVariables.Options.Set(ACTIVE,true);

      if(ContactFace == 2) //tip surface
	GetGeometry()[0].Set(TO_SPLIT);

      //get contact properties and parameters
      CalculateContactFactors( rVariables );

    }
    else{
      
      rVariables.Options.Set(ACTIVE,false);
      
    }

  
    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************


  void PointRigidContactPenalty2DCondition::CalculateContactFactors(GeneralVariables &rVariables)
  {

    KRATOS_TRY
      
      WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){
	    
	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();
	    
	  distance += norm_2(GetGeometry()[0]-Neighb_Point);

	  counter ++;
	}
      }

    distance /= counter;
    
    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

    rVariables.Penalty.Normal = distance * 20 * PenaltyParameter * ElasticModulus;
      
    
    
    KRATOS_CATCH("")
      }
  
  //***********************************************************************************
  //************************************************************************************

  double& PointRigidContactPenalty2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
  { 
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( dimension == 2 ) 
      rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
  }


  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
								GeneralVariables& rVariables,
								double& rIntegrationWeight)

  {
    KRATOS_TRY

      //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      noalias(rLeftHandSideMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod_2(rVariables.Surface.Normal, rVariables.Surface.Normal);
      
    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(2,2);   

    }

    KRATOS_CATCH( "" )
      }


  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 GeneralVariables& rVariables,
									 double& rIntegrationWeight)

  {
    KRATOS_TRY

      //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      rRightHandSideVector = (-1) * (rVariables.Penalty.Normal * rVariables.Gap.Normal) * rVariables.Surface.Normal * rIntegrationWeight;
      
    }
    else{

      rRightHandSideVector = ZeroVector(2);
    
    }

    //KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  inline Condition::MatrixType PointRigidContactPenalty2DCondition::outer_prod_2(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
  {
    Condition::MatrixType A(2,2);
    A(0,0)=a[0]*b[0];
    A(0,1)=a[0]*b[1];
    A(1,0)=a[1]*b[0];
    A(1,1)=a[1]*b[1];
    //A(0,2)=a[0]*b[2];
    //A(1,2)=a[1]*b[2];
    //A(2,0)=a[2]*b[0];
    //A(2,1)=a[2]*b[1];
    //A(2,2)=a[2]*b[2];
    return A;
  }




} // Namespace Kratos



