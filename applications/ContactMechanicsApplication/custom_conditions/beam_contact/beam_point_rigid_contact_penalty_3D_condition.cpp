//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/thermal_contact/beam_point_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactPenalty3DCondition::BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : BeamPointRigidContactCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactPenalty3DCondition::BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : BeamPointRigidContactCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactPenalty3DCondition::BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : BeamPointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactPenalty3DCondition::BeamPointRigidContactPenalty3DCondition( BeamPointRigidContactPenalty3DCondition const& rOther )
  : BeamPointRigidContactCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer BeamPointRigidContactPenalty3DCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<BeamPointRigidContactPenalty3DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************************************************************
  //************************************************************************************


  BeamPointRigidContactPenalty3DCondition::~BeamPointRigidContactPenalty3DCondition()
  {

  }

  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactPenalty3DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_TRY


       ConditionVariables ContactVariables;

       SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], ContactVariables.Gap.Normal, ContactVariables.Gap.Tangent, ContactVariables.Surface.Normal, ContactVariables.Surface.Tangent, ContactVariables.RelativeDisplacement);

       //to perform contact with a tube radius must be set
       BoxParameters.SetRadius(GetGeometry()[0].GetValue(CROSS_SECTION_AREA));

       if ( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ) {

  	   const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  	   array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

	   double TangentRelativeMovement = 0;
	   TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, ContactVariables );

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

       mTangentialVariables.FrictionCoefficient        =  GetProperties()[FRICTION_COEFFICIENT];//0.3;
       mTangentialVariables.DynamicFrictionCoefficient =  GetProperties()[MU_DYNAMIC];//0.2;
       mTangentialVariables.StaticFrictionCoefficient  =  GetProperties()[MU_STATIC];//0.3;

       //std::cout<<" Friction Coef "<<mTangentialVariables.FrictionCoefficient<<" dynamic "<<mTangentialVariables.DynamicFrictionCoefficient<<" static "<<mTangentialVariables.StaticFrictionCoefficient<<std::endl;

       ClearNodalForces();


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************
  void BeamPointRigidContactPenalty3DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
    //added to control force evolution per step:
    // array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
    // mTangentialVariables.PreviousTangentForceModulus = norm_2(ContactForce);

    ClearNodalForces();
  }

  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactPenalty3DCondition::InitializeConditionVariables (ConditionVariables& rVariables,
									const ProcessInfo& rCurrentProcessInfo)
  {



  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void BeamPointRigidContactPenalty3DCondition::CalculateKinematics(ConditionVariables& rVariables,
								    const ProcessInfo& rCurrentProcessInfo,
								    const double& rPointNumber)
  {
    KRATOS_TRY

    SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, rVariables.RelativeDisplacement);

    //to perform contact with a tube radius must be set
    BoxParameters.SetRadius(GetGeometry()[0].GetValue(CROSS_SECTION_AREA));

    if( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ){

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


  void BeamPointRigidContactPenalty3DCondition::CalculateContactFactors(ConditionVariables &rVariables)
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

    //reduction of the penalty parameter:
    PenaltyParameter *= 1e-6;

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;


    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;

    // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

    KRATOS_CATCH( "" )
      }


  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactPenalty3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
								    ConditionVariables& rVariables,
								    double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      MatrixType Kuug = ZeroMatrix(dimension,dimension);

      noalias(Kuug) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);

      //Building the Local Stiffness Matrix
      BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kuug, 0, 0 );

      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      //this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight );
      // std::cout<<std::endl;
      //std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(dimension*2,dimension*2);

    }

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
      }

  //************* Tangent Contact Force constitutive matrix      **********************
  //***********************************************************************************

   void BeamPointRigidContactPenalty3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables, double& rIntegrationWeight)
   {
       KRATOS_TRY

       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
       double NormalForceModulus = 0;
       NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

       double TangentRelativeMovement = 0;
       TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );

       double TangentForceModulus = this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables );

       if( fabs(TangentForceModulus) >= 1e-25 ){


	 MatrixType Kuug = ZeroMatrix(dimension,dimension);

	 if ( mTangentialVariables.Slip ) {
	   //simpler expression:
	   Kuug -=  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) );

	   //added extra term, maybe not necessary
	   //kuug -=  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) + rVariables.Gap.Normal * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

	 }
	 else {
	   Kuug +=  rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

	 }

	 //Building the Local Stiffness Matrix
	 BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kuug, 0, 0 );


	 //The geometric part of the torque contribution must be added here:




       }

       KRATOS_CATCH( "" )

   }

  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactPenalty3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									     ConditionVariables& rVariables,
									     double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       //this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

    }
    else{

      rRightHandSideVector = ZeroVector(dimension * 2);

    }

    //KRATOS_WATCH( rRightHandSideVector )

    KRATOS_CATCH( "" )
      }

  //**************************** Calculate Normal Contact Force ***********************
  //***********************************************************************************

  void BeamPointRigidContactPenalty3DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
										  ConditionVariables& rVariables,
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

  void BeamPointRigidContactPenalty3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
										   ConditionVariables& rVariables,
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

       //Contact tangent force in beam center
       for (unsigned int i = 0; i < dimension ; ++i) {
	 rRightHandSideVector[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
	 ContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
       }


       double SectionMeanRadius = GetProperties()[CROSS_SECTION_AREA];
       PointType RadiusVector  = SectionMeanRadius * rVariables.Surface.Normal;
       PointType ContactTorque;
       MathUtils<double>::CrossProduct( ContactTorque, rVariables.Surface.Tangent, RadiusVector );
       ContactTorque *= TangentForceModulus;

       //std::cout<<" [ContactTorque]: "<<ContactTorque<<" [TangentForceModulus]: "<<TangentForceModulus<<std::endl;

       //Contact torque due to contact tangent force on beam surface
       for (unsigned int i = dimension; i < (dimension * 2); ++i) {
	 rRightHandSideVector[i] += ContactTorque[i];
       }



       //std::cout<< "["<<mTangentialVariables.Sign<<"] Tangent Force Node ["<<GetGeometry()[0].Id()<<" ]:"<<TangentForceModulus<<" RelativeMovement: "<<rVariables.Gap.Tangent<<" SLIP ["<<mTangentialVariables.Slip<<"]"<<std::endl;


       GetGeometry()[0].UnSetLock();


  }

  //**************************** Calculate Normal Force Modulus ***********************
  //***********************************************************************************

  double& BeamPointRigidContactPenalty3DCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, ConditionVariables& rVariables )
  {

        rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal);

	//std::cout<<" [NormalForceModulus]: "<<rNormalForceModulus<<" [Penalty:"<<rVariables.Penalty.Normal<<" Gap:"<<rVariables.Gap.Normal<<"]"<<std::endl;

	return rNormalForceModulus;

  }

  //**************************** Calculate Tangent Force Modulus **********************
  //***********************************************************************************

  double& BeamPointRigidContactPenalty3DCondition::CalculateTangentRelativeMovement( double& rTangentRelativeMovement, ConditionVariables& rVariables )
  {
       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       const array_1d<double, 3> & CurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
       const array_1d<double, 3> & PreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
       array_1d<double, 3 > DeltaDisplacement           =  CurrentDisplacement-PreviousDisplacement;


       //Calculate the relative rotation
       const array_1d<double, 3> & CurrentRotation      =  GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
       const array_1d<double, 3> & PreviousRotation     =  GetGeometry()[0].FastGetSolutionStepValue(ROTATION, 1);
       array_1d<double, 3 > DeltaRotation               =  CurrentRotation-PreviousRotation;


       double SectionMeanRadius = GetProperties()[CROSS_SECTION_AREA];
       array_1d<double, 3 > RadiusVector;
       for(unsigned int i = 0; i < 3; i++)
	 {
	   RadiusVector[i] =  SectionMeanRadius * rVariables.Surface.Normal[i];
	 }

       array_1d<double, 3 > DeltaRotationDisplacement;
       MathUtils<double>::CrossProduct( DeltaRotationDisplacement, DeltaRotation, RadiusVector );

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

       PointType WallDisplacement = mTangentialVariables.DeltaTime * this->mpRigidWall->GetVelocity();

       rTangentRelativeMovement = 0.0;
       PointType TotalTangentRelativeMovement = ZeroVector(dimension);

       for (unsigned int i = 0; i < dimension; ++i)
	 {
	   TotalTangentRelativeMovement[i] = DeltaDisplacement[i] + DeltaRotationDisplacement[i] - WallDisplacement[i];
	 }

       rTangentRelativeMovement = MathUtils<double>::Norm3(TotalTangentRelativeMovement);

       if( rTangentRelativeMovement !=0 )
	 rVariables.Surface.Tangent = TotalTangentRelativeMovement/rTangentRelativeMovement;


       rVariables.Gap.Tangent = rTangentRelativeMovement;


       return rTangentRelativeMovement;

  }

  //**************************** Check Coulomb law for Tangent Contact Force **********
  //***********************************************************************************

  double BeamPointRigidContactPenalty3DCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , ConditionVariables& rVariables)
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

  double BeamPointRigidContactPenalty3DCondition::CalculateFrictionCoefficient(double & rTangentRelativeMovement)
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


} // Namespace Kratos
