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

  
  void PointRigidContactPenalty2DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {

       array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

       GeneralVariables AuxiliarVariables;
       int ContactFace = 0;

       bool Aux = this->mpRigidWall->IsInside( GetGeometry()[0], AuxiliarVariables.Gap.Normal, AuxiliarVariables.Gap.Tangent, AuxiliarVariables.Surface.Normal, AuxiliarVariables.Surface.Tangent, ContactFace );



       if (Aux) {
           double TangentialModul = 0.0;
           for (unsigned int i = 0; i < 3; ++i) {
              TangentialModul += ContactForce[i] * AuxiliarVariables.Surface.Tangent[i];
           }

           mTangentialVariables.LastStepValue = TangentialModul;
        }
        else {
           mTangentialVariables.LastStepValue = 0.0;
        } 
 
        double DeltaTime ;
        DeltaTime = rCurrentProcessInfo[DELTA_TIME];
        mTangentialVariables.DeltaTime = DeltaTime;
     
    KRATOS_TRY

    ClearNodalForces();

    KRATOS_CATCH( "" )

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

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){
	    
	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();
	    
	  distance += norm_2(Contact_Point-Neighb_Point);

	  counter ++;
	}
      }

    distance /= counter;
    
    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

    rVariables.Penalty.Normal = distance * 2000.0 * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal ;  
    
    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;
    
    KRATOS_CATCH( "" )
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
      
      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<outer_prod_2(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);
      // std::cout<<std::endl;
      // std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(2,2);   

    }
 
    KRATOS_CATCH( "" )
      }

  //************* Tangent Contact Force constitutive matrix      **********************
  //***********************************************************************************

   void PointRigidContactPenalty2DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
   {
  
       if (mTangentialVariables.FrictionActive ) {
          double Mu = 0.3;
          //rLeftHandSideMatrix -= mTangentialVariables.Sign * Mu * rVariables.Penalty.Normal * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Normal);
          rLeftHandSideMatrix += mTangentialVariables.Sign * Mu * rVariables.Penalty.Normal * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Normal);
 
       }
       else {
          rLeftHandSideMatrix += rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Tangent);
       }

   }

  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 GeneralVariables& rVariables,
									 double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight, dimension);
       this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight, dimension);

    }
    else{

      rRightHandSideVector = ZeroVector(2);
    
    }


    KRATOS_CATCH( "" )
      }

  //**************************** Calculate Normal Contact Force ***********************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									GeneralVariables& rVariables,
									double& rIntegrationWeight,
									const unsigned int& rDimension)
  {
      double ContactForceModulus = (-1) * (rVariables.Penalty.Normal * rVariables.Gap.Normal) * rIntegrationWeight;
      for(unsigned int j=0; j<rDimension; j++)
	{
	  rRightHandSideVector[j] = ContactForceModulus * rVariables.Surface.Normal[j];
	}

      
      GetGeometry()[0].SetLock();

      array_1d<double, 3 > &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


      for(unsigned int j=0; j<rDimension; j++)
	{
	  ContactForce[j] = rRightHandSideVector[j];
	}

      GetGeometry()[0].UnSetLock();

      //std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rIntegrationWeight "<<rIntegrationWeight<<std::endl;
      //std::cout<<std::endl;
      //std::cout<<" Fcont "<<rRightHandSideVector<<std::endl;

  }

  //**************************** Calculate Tangent Contact Force **********************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									GeneralVariables& rVariables,
									double& rIntegrationWeight,
									const unsigned int& rDimension)
  {
  
       array_1d<double, 3 > &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


       const array_1d<double, 3> & ActualDisplacement  =  GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT);
       const array_1d<double, 3> & PreviousDisplacement = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT, 1);

       double IncrementalStick = 0.0;
  
       for (unsigned int i = 0; i < rDimension; ++i)
             IncrementalStick += (ActualDisplacement[i]-PreviousDisplacement[i]) * rVariables.Surface.Tangent(i);


       double StructureRelativeMovement = 0.0;
       {
          array_1d<double, 3> Velocity;
          Velocity = this->mpRigidWall->Velocity();
          double DT = 1.25e-8;
          for (unsigned j = 0; j < rDimension; ++j)
              StructureRelativeMovement -= Velocity[j] * rVariables.Surface.Tangent(j)*DT;
       }

       IncrementalStick += StructureRelativeMovement;

       double TangentModulus ;
       TangentModulus= IncrementalStick*rVariables.Penalty.Tangent * rIntegrationWeight;

       TangentModulus -= mTangentialVariables.LastStepValue; 

       this->ComputeCoulombCondition(TangentModulus, ContactForce, rDimension, mTangentialVariables.FrictionActive );

       for (unsigned int i = 0; i < rDimension ; ++i) {
           rRightHandSideVector[i] -= TangentModulus * rVariables.Surface.Tangent(i);
           ContactForce[i] -= TangentModulus * rVariables.Surface.Tangent(i);
       }
 
       GetGeometry()[0].UnSetLock();


  }

  //**************************** Check Coulomb law for  Tangent Contact Force *********
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::ComputeCoulombCondition(double & rTangentModulus, const array_1d<double, 3>& rNormalForceVector , const int& rDimension, bool& rFrictionBool)
  {

       double NormalModulus = 0.0;
       for (int i = 0; i < rDimension; ++i)
           NormalModulus += pow( rNormalForceVector[i], 2.0);

       NormalModulus = pow(NormalModulus, 0.5);

       double Mu = 0.3;

       if ( fabs(rTangentModulus) >  Mu* fabs(NormalModulus) ) {
           mTangentialVariables.Sign = rTangentModulus / fabs(rTangentModulus) ; 
           rTangentModulus =   mTangentialVariables.Sign * Mu * fabs(NormalModulus) ;
           rFrictionBool = true;
       }
       else {
           rFrictionBool = false;
       }

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



