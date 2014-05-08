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
#include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : PointRigidContactCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : PointRigidContactCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : PointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition( PointRigidContactPenalty3DCondition const& rOther )
  : PointRigidContactCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer PointRigidContactPenalty3DCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new PointRigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************************************************************
  //************************************************************************************


  PointRigidContactPenalty3DCondition::~PointRigidContactPenalty3DCondition()
  {

  }

  //************************************************************************************
  //************************************************************************************
 
  void PointRigidContactPenalty3DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
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

       mTangentialVariables.Sign = 1;
       
       mTangentialVariables.FrictionCoefficient = 0.0;//0.3;
       mTangentialVariables.DynamicFrictionCoefficient = 0.0;//0.2;
       mTangentialVariables.StaticFrictionCoefficient  = 0.0;//0.3;


       ClearNodalForces();

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty3DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
    CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
    CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
    CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
  
    //added to control force evolution per step:
    // array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
    // mTangentialVariables.PreviousTangentForceModulus = norm_2(ContactForce);

    ClearNodalForces();
  }
    
  //************************************************************************************
  //************************************************************************************

  void PointRigidContactPenalty3DCondition::InitializeGeneralVariables (GeneralVariables& rVariables, 
									const ProcessInfo& rCurrentProcessInfo)
  {



  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************
  
  void PointRigidContactPenalty3DCondition::CalculateKinematics(GeneralVariables& rVariables,
								const double& rPointNumber)
  {
    KRATOS_TRY
    
    int ContactFace = 0; //free surface
      
    if( this->mpRigidWall->IsInside( GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, ContactFace ) ){

      rVariables.Options.Set(ACTIVE,true);

      //get contact properties and parameters
      this->CalculateContactFactors( rVariables );

    }
    else{
      
      rVariables.Options.Set(ACTIVE,false);
      
    }

  
    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************


  void PointRigidContactPenalty3DCondition::CalculateContactFactors(GeneralVariables &rVariables)
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

    if( counter != 0 )
      distance /= counter;

    if( distance == 0 )
      distance = 1;
    
    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;  
    

    // std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;
    
    // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

    KRATOS_CATCH( "" )
      }
  

  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
								GeneralVariables& rVariables,
								double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      noalias(rLeftHandSideMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);
      
      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);
      // std::cout<<std::endl;
      //std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(dimension,dimension);   

    }
 
    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
      }

  //************* Tangent Contact Force constitutive matrix      **********************
  //***********************************************************************************

   void PointRigidContactPenalty3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
   {

       double NormalForceModulus = 0;
       NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

       double TangentRelativeMovement = 0;
       TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );
       
       double TangentForceModulus = this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables );

       if( fabs(TangentForceModulus) >= 1e-25 ){
       
	 if ( mTangentialVariables.Slip ) {
	   //simpler expression:
	   rLeftHandSideMatrix -=  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) );

	   //added extra term, maybe not necessary
	   //rLeftHandSideMatrix -=  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) + rVariables.Gap.Normal * custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

	 }
	 else {
	   rLeftHandSideMatrix +=  rVariables.Penalty.Tangent * rIntegrationWeight * custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

	 }

       }
       

   }

  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 GeneralVariables& rVariables,
									 double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

    }
    else{

      rRightHandSideVector = ZeroVector(dimension);
    
    }

    //KRATOS_WATCH( rRightHandSideVector )

    KRATOS_CATCH( "" )
      }

  //**************************** Calculate Normal Contact Force ***********************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									      GeneralVariables& rVariables,
									      double& rIntegrationWeight)
  {
      
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

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


      //std::cout<<" Fcont "<<rRightHandSideVector<<std::endl;

  }

  //**************************** Calculate Tangent Contact Force **********************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									       GeneralVariables& rVariables,
									       double& rIntegrationWeight)
  {

       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       double NormalForceModulus = 0;
       NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

       double TangentRelativeMovement = 0;
       TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );
       
       double TangentForceModulus =  this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables );

       TangentForceModulus *= (-1) * rIntegrationWeight;

       GetGeometry()[0].SetLock();

       array_1d<double, 3 > & ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

       for (unsigned int i = 0; i < dimension ; ++i) {
           rRightHandSideVector[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
           ContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
       }


       //std::cout<< "["<<mTangentialVariables.Sign<<"] Tangent Force Node ["<<GetGeometry()[0].Id()<<" ]:"<<TangentForceModulus<<" RelativeMovement: "<<rVariables.Gap.Tangent<<" SLIP ["<<mTangentialVariables.Slip<<"]"<<std::endl; 

 
       GetGeometry()[0].UnSetLock();


  }

  //**************************** Calculate Normal Force Modulus ***********************
  //***********************************************************************************

  double& PointRigidContactPenalty3DCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, GeneralVariables& rVariables )
  {

        rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal); 

	return rNormalForceModulus;

  }

  //**************************** Calculate Tangent Force Modulus **********************
  //***********************************************************************************

  double& PointRigidContactPenalty3DCondition::CalculateTangentRelativeMovement( double& rTangentRelativeMovement, GeneralVariables& rVariables )
  {
       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       const array_1d<double, 3> & CurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
       const array_1d<double, 3> & PreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
       array_1d<double, 3 > DeltaDisplacement            = CurrentDisplacement-PreviousDisplacement;

       // bool Regularization = false;

       // if( Regularization == true ){

       // 	 //regularization taking the contigous boundary segments
       // 	 WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

       // 	 array_1d<double, 3 > NeighbDeltaDisplacement;

       // 	 double counter = 0;
       
       // 	 for(unsigned int i = 0; i < rN.size(); i++)
       // 	   {
       // 	     if(rN[i].Is(BOUNDARY)){
	          
       // 	       const array_1d<double, 3> & NeighbCurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
       // 	       const array_1d<double, 3> & NeighbPreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

       // 	       array_1d<double, 3 > NeighbDeltaDisplacement           =  NeighbCurrentDisplacement-NeighbPreviousDisplacement;

       // 	       DeltaDisplacement += NeighbDeltaDisplacement;
     
       // 	       counter ++;
       // 	     }
       // 	   }
       //        if( counter!= 0)
       // 	    DeltaDisplacement /= counter;

       // }

       VectorType WallDisplacement = mTangentialVariables.DeltaTime * this->mpRigidWall->Velocity();
       
       rTangentRelativeMovement = 0.0;
       double WallTangentRelativeMovement    = 0.0;

       for (unsigned int i = 0; i < dimension; ++i)
	 {
	   rTangentRelativeMovement    += DeltaDisplacement[i] * rVariables.Surface.Tangent[i];
	   WallTangentRelativeMovement += WallDisplacement[i]  * rVariables.Surface.Tangent[i];  
	 }



       rTangentRelativeMovement -= WallTangentRelativeMovement;      

       rVariables.Gap.Tangent = rTangentRelativeMovement;


       return rTangentRelativeMovement;

  }

  //**************************** Check Coulomb law for Tangent Contact Force **********
  //***********************************************************************************

  double PointRigidContactPenalty3DCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , GeneralVariables& rVariables)
  {
       mTangentialVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement);

 
       double TangentForceModulus = rVariables.Penalty.Tangent * rVariables.Gap.Tangent; //+ mTangentialVariables.PreviousTangentForceModulus; 
     

 
       if ( fabs(TangentForceModulus) >  mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) && fabs(rVariables.Gap.Tangent) > 1e-200) {

	 mTangentialVariables.Sign =  rVariables.Gap.Tangent/ fabs( rVariables.Gap.Tangent ) ; 

	 TangentForceModulus =  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) ;
	 mTangentialVariables.Slip = true;

       }
       else {
	 mTangentialVariables.Slip = false;
       }


       return TangentForceModulus;
  }

  

  //**************************** Check friction coefficient ***************************
  //***********************************************************************************

  double PointRigidContactPenalty3DCondition::CalculateFrictionCoefficient(double & rTangentRelativeMovement)
  {

       //---FRICTION LAW in function of the relative sliding velocity ---//

       double Velocity = 0;
       Velocity = rTangentRelativeMovement / mTangentialVariables.DeltaTime;

 
       //Addicional constitutive parameter  C
       //which describes how fast the static coefficient approaches the dynamic:
       double C=0.1;

       //Addicional constitutive parameter  E
       //regularization parameter (->0, classical Coulomb law)
       double E=0.01;


       double FrictionCoefficient = mTangentialVariables.DynamicFrictionCoefficient + ( mTangentialVariables.StaticFrictionCoefficient-  mTangentialVariables.DynamicFrictionCoefficient ) * exp( (-1) * C * fabs(Velocity) );


       //Square root regularization
       FrictionCoefficient *= fabs(Velocity)/sqrt( ( Velocity * Velocity ) + ( E * E ) );
       
       //Hyperbolic regularization
       //FrictionCoefficient *= tanh( fabs(Velocity)/E );

       return FrictionCoefficient;

  }

  //************************************************************************************
  //************************************************************************************

  inline Condition::MatrixType PointRigidContactPenalty3DCondition::custom_outer_prod(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
  {
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Condition::MatrixType A(dimension,dimension);
    
    A(0,0)=a[0]*b[0];
    A(0,1)=a[0]*b[1];
    A(1,0)=a[1]*b[0];
    A(1,1)=a[1]*b[1];
    if( dimension == 3 ){
      A(0,2)=a[0]*b[2];
      A(1,2)=a[1]*b[2];
      A(2,0)=a[2]*b[0];
      A(2,1)=a[2]*b[1];
      A(2,2)=a[2]*b[2];
    }

    return A;
  }




} // Namespace Kratos



