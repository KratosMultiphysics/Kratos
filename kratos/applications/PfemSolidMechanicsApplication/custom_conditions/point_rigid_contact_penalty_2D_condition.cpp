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
#include "includes/define.h"
#include "custom_conditions/point_rigid_contact_penalty_2D_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : Condition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : Condition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  :Condition(NewId, pGeometry, pProperties)
  {
    
    mpRigidWall = pRigidWall;
    
  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty2DCondition::PointRigidContactPenalty2DCondition( PointRigidContactPenalty2DCondition const& rOther )
  : Condition(rOther)
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
  void PointRigidContactPenalty2DCondition::SetRigidWall(SpatialBoundingBox::Pointer pRigidWall)
  {
    mpRigidWall = pRigidWall;
  }

  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty2DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
    CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
    CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
    CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;

    array_1d<double, 3 > & ContactForceNormal  = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
    ContactForceNormal.clear();

    array_1d<double, 3 > & ContactForceTangent = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
    ContactForceTangent.clear();
  }

  //************************************************************************************
  //************************************************************************************

  void PointRigidContactPenalty2DCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
  }


  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty2DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      
      ContactVariables Contact;
    
    if( mpRigidWall->IsInside( GetGeometry()[0], Contact.Gap.Normal, Contact.Gap.Tangent, Contact.Surface.Normal, Contact.Surface.Tangent ) ){

      Contact.Options.Set(ACTIVE,true);


      if(rRightHandSideVector.size() != 2)
	rRightHandSideVector.resize(2,false);

      array_1d<double, 3>& ContactForceNormal = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
      //array_1d<double, 3>& ContactForceTangent = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
    
      //get contact properties and parameters
      CalculateContactFactors(Contact);

      //compute contact force
      ContactForceNormal = (-1) * (Contact.Penalty.Normal * Contact.Gap.Normal) * Contact.Surface.Normal;


      double IntegrationWeight = 1; //GetProperties()[ THICKNESS ];
      rRightHandSideVector[0]  = ContactForceNormal[0] * IntegrationWeight;
      rRightHandSideVector[1]  = ContactForceNormal[1] * IntegrationWeight;


      rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;   

    }
    else{


      rRightHandSideVector = ZeroVector(2);
     
      Contact.Options.Set(ACTIVE,false);
     
    }

 
    KRATOS_CATCH("")
      }

  //************************************************************************************
  //************************************************************************************


  void PointRigidContactPenalty2DCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      ContactVariables Contact;
    
    if( mpRigidWall->IsInside( GetGeometry()[0], Contact.Gap.Normal, Contact.Gap.Tangent, Contact.Surface.Normal, Contact.Surface.Tangent ) ){

      Contact.Options.Set(ACTIVE,true);

      if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);
    
      //get contact properties and parameters
      CalculateContactFactors(Contact);

      //compute contact stiffness
      double IntegrationWeight = 1; //GetProperties()[ THICKNESS ];
      noalias(rLeftHandSideMatrix) = Contact.Penalty.Normal * IntegrationWeight  * outer_prod_2(Contact.Surface.Normal, Contact.Surface.Normal);

      rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;   

    }
    else{

      if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);
     
      rLeftHandSideMatrix= ZeroMatrix(2,2);

      Contact.Options.Set(ACTIVE,false);
     
    }




    KRATOS_CATCH("")

      }

  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty2DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      ContactVariables Contact;
    
    if( mpRigidWall->IsInside( GetGeometry()[0], Contact.Gap.Normal, Contact.Gap.Tangent, Contact.Surface.Normal, Contact.Surface.Tangent ) ){

      Contact.Options.Set(ACTIVE,true);

      if(rRightHandSideVector.size() != 2)
	rRightHandSideVector.resize(2,false);

      if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);
    
 
      array_1d<double, 3>& ContactForceNormal = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
      //array_1d<double, 3>& ContactForceTangent = GetGeometry()[0].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
    
      //get contact properties and parameters
      CalculateContactFactors(Contact);
  
      //compute contact force
      ContactForceNormal = (-1) * (Contact.Penalty.Normal * Contact.Gap.Normal) * Contact.Surface.Normal;


      double IntegrationWeight = 1; //GetProperties()[ THICKNESS ];
      rRightHandSideVector[0]  = ContactForceNormal[0] * IntegrationWeight;
      rRightHandSideVector[1]  = ContactForceNormal[1] * IntegrationWeight;

      //compute contact stiffness
      noalias(rLeftHandSideMatrix) = Contact.Penalty.Normal * IntegrationWeight  * outer_prod_2(Contact.Surface.Normal, Contact.Surface.Normal);

      rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;   

    }
    else{

      rRightHandSideVector = ZeroVector(2);
	
      rLeftHandSideMatrix  = ZeroMatrix(2,2);

      Contact.Options.Set(ACTIVE,false);
     
    }


    KRATOS_CATCH("")
      }


  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty2DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 2;
    rResult.resize(number_of_nodes*dim);
    for (int i=0; i<number_of_nodes; i++)
      {
        index = i*dim;
        rResult[index]   = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());

      }
  }

  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty2DCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
  {
    unsigned int dim = 2;
    ConditionalDofList.resize(GetGeometry().size()*dim);
    unsigned int index;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
      {

        index = i*dim;
        ConditionalDofList[index]   = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

      }
  }

  //************************************************************************************
  //************************************************************************************


  void PointRigidContactPenalty2DCondition::CalculateContactFactors(ContactVariables &rContact)
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

    rContact.Penalty.Normal = distance * 20 * PenaltyParameter * ElasticModulus;
      
    
    
    KRATOS_CATCH("")
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



