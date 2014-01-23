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
 
  void PointRigidContactPenalty2DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_TRY

       
       GeneralVariables ContactVariables;
       int ContactFace = 0;

       if ( this->mpRigidWall->IsInside( GetGeometry()[0], ContactVariables.Gap.Normal, ContactVariables.Gap.Tangent, ContactVariables.Surface.Normal, ContactVariables.Surface.Tangent, ContactFace ) ) {

  	   const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  	   array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

           mTangentialVariables.PreviousTangentForceModulus = 0.0;
           for (unsigned int i = 0; i < dimension; ++i) {
              mTangentialVariables.PreviousTangentForceModulus += ContactForce[i] * ContactVariables.Surface.Tangent[i];
           }

        }
        else {
           mTangentialVariables.PreviousTangentForceModulus = 0.0;
        } 
 
       mTangentialVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];
       mTangentialVariables.FrictionCoefficient = 0.3;

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

      if(ContactFace == 2){ //tip surface
	GetGeometry()[0].Set(TO_SPLIT);
	//std::cout<<" Node ["<<GetGeometry()[0].Id()<<"] set TO_SPLIT "<<std::endl;
      }
      else{
	GetGeometry()[0].Set(TO_SPLIT,false);
      }

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
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal * 0.1;  
    
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

       double NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

       double TangentForceModulus = this->CalculateTangentForceModulus( TangentForceModulus, rVariables );
       
       this->CalculateCoulombsLaw( TangentForceModulus, NormalForceModulus, mTangentialVariables.Slip );

       //std::cout<<" Is Slip "<<mTangentialVariables.Slip<<std::endl;

       double zero_force = 1;
       if( fabs(TangentForceModulus) <= 1e-8 )
	 zero_force = 0;

       if ( mTangentialVariables.Slip ) {
	 rLeftHandSideMatrix += zero_force * mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * outer_prod_2(rVariables.Surface.Normal, rVariables.Surface.Tangent);
	 
	 //std::cout<<" Ktangent 1 "<<zero_force * mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Normal)<<std::endl; 
       }
       else {
	 rLeftHandSideMatrix += zero_force * rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

	 //std::cout<<" Ktangent 2 "<<zero_force * rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod_2(rVariables.Surface.Tangent, rVariables.Surface.Tangent)<<std::endl;

       }


   }

  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 GeneralVariables& rVariables,
									 double& rIntegrationWeight)

  {
    KRATOS_TRY

     if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

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
									      double& rIntegrationWeight)
  {
      
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      NormalForceModulus *= (-1) * rIntegrationWeight;
   
      for(unsigned int j = 0; j < dimension; j++)
	{
	  rRightHandSideVector[j] = NormalForceModulus * rVariables.Surface.Normal[j];
	}

      
      GetGeometry()[0].SetLock();

      array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


      for(unsigned int j = 0; j < dimension; j++)
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
									       double& rIntegrationWeight)
  {

       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       double NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

       double TangentForceModulus = this->CalculateTangentForceModulus( TangentForceModulus, rVariables );
       
       this->CalculateCoulombsLaw( TangentForceModulus, NormalForceModulus, mTangentialVariables.Slip );

       TangentForceModulus *= (-1) * rIntegrationWeight;

       GetGeometry()[0].SetLock();

       array_1d<double, 3 > & ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

       for (unsigned int i = 0; i < dimension ; ++i) {
           rRightHandSideVector[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
           ContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
       }

       //std::cout<<" TangentForce "<<TangentForceModulus * rVariables.Surface.Tangent<<std::endl;
 
       GetGeometry()[0].UnSetLock();


  }

  //**************************** Calculate Normal Force Modulus ***********************
  //***********************************************************************************

  double& PointRigidContactPenalty2DCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, GeneralVariables& rVariables )
  {

        rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal); 

       return rNormalForceModulus;

  }

  //**************************** Calculate Tangent Force Modulus **********************
  //***********************************************************************************

  double& PointRigidContactPenalty2DCondition::CalculateTangentForceModulus ( double& rTangentForceModulus, GeneralVariables& rVariables )
  {
       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       const array_1d<double, 3> & CurrentDisplacement  =  GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT);
       const array_1d<double, 3> & PreviousDisplacement =  GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT, 1);
       array_1d<double, 3 > DeltaDisplacement            = CurrentDisplacement-PreviousDisplacement;

       VectorType WallDisplacement = mTangentialVariables.DeltaTime * this->mpRigidWall->Velocity();
       
       double TangentRelativeMovement = 0.0;
       double WallRelativeMovement    = 0.0;

       for (unsigned int i = 0; i < dimension; ++i)
	 {
	   TangentRelativeMovement += DeltaDisplacement[i] * rVariables.Surface.Tangent[i];
	   WallRelativeMovement    += WallDisplacement[i]  * rVariables.Surface.Tangent[i];  
	 }

       // std::cout<<" TangentRelativeMovement "<< TangentRelativeMovement <<std::endl;
       // std::cout<<" WallRelativeMovement "<< WallRelativeMovement <<std::endl;

       TangentRelativeMovement -= WallRelativeMovement;      

       //rTangentForceModulus = TangentRelativeMovement * rVariables.Penalty.Tangent + mTangentialVariables.PreviousTangentForceModulus; 

       rTangentForceModulus = TangentRelativeMovement * rVariables.Penalty.Tangent;

       return rTangentForceModulus;

  }

  //**************************** Check Coulomb law for Tangent Contact Force **********
  //***********************************************************************************

  void PointRigidContactPenalty2DCondition::CalculateCoulombsLaw(double & rTangentForceModulus, double & rNormalForceModulus , bool& rSlip)
  {
       //std::cout<< " rTangentForceModulus "<<rTangentForceModulus<< " rNormalForceModulus "<<rNormalForceModulus<<std::endl;

       if ( fabs(rTangentForceModulus) >  mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) ) {
	   mTangentialVariables.Sign = rTangentForceModulus / fabs(rTangentForceModulus) ; 
           rTangentForceModulus =  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) ;
	   //std::cout<<" Slip Force Modulus : "<<rTangentForceModulus<<std::endl;
           rSlip = true;
       }
       else {
	   //std::cout<<" Stick Force Modulus : "<<rTangentForceModulus<<std::endl;
	   rSlip = false;
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



