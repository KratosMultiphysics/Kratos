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

#include "custom_conditions/custom_friction_laws/friction_law.hpp"

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


      // Compute the neighbour distance, then a stress-"like" may be computed.
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

      distance /= 2.0;
      mTangentialVariables.Neighb_distance = distance; 
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
      double PenaltyParameter = 1;
      if( GetProperties().Has(PENALTY_PARAMETER) )
	PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

      double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

      // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
      if (ElasticModulus <= 1.0e-5) {
         std::vector<double> mModulus;
         ProcessInfo SomeProcessInfo;
         WeakPointerVector<Element > rN = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
         ElasticModulus = 0.0;
         int NumberOfSum = 0;
         for ( unsigned int i = 0; i < rN.size(); i++)
         {
            rN[i].CalculateOnIntegrationPoints(SIMILAR_YOUNG_MODULUS, mModulus, SomeProcessInfo);
            ElasticModulus += mModulus[0];
            NumberOfSum += 1;
         }
         ElasticModulus /= double(NumberOfSum);
      }

      double factor = 4;
      if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
         int order = (int)((-1) * std::log10(distance) + 1) ;
         distance *= factor * pow(10,order);
      }

      rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
      double PenaltyRatio = 1;
      if( GetProperties().Has(TANGENTIAL_PENALTY_RATIO) )
	PenaltyRatio = GetProperties()[PENALTY_PARAMETER];

      rVariables.Penalty.Tangent = rVariables.Penalty.Normal * PenaltyRatio ; ///20.0;  

      //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;

      //std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

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

         return; 

         /*
         mTangentialVariables.IntegrationWeight = rIntegrationWeight;

         double NormalForceModulus = 0;
         NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables);

         double TangentRelativeMovement = 0;
         TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );

         double TangentForceModulus = this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables); // ALSO COMPUTES the TANGENT MATRIX ( and is saved in the rVariables.TangentMatrix.Normal and .Tangent)

         // OBS: rVariables.TangentMatrix.Normal is the variation of the tangent stress with respect to the normal stress
         //      rVariables.TangentMatrix.Normal is the variation of the tangent stress with respect the tangent gap
         if ( fabs(TangentForceModulus) >= 1e-25) {

            if ( mTangentialVariables.Slip) {
               // YOU KNOW
               //rLeftHandSideMatrix -= mTangentialVariables.Sign* rIntegrationWeight * ( rVariables.TangentMatrix.Normal * rVariables.Penalty.Normal * custom_outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Normal)  - rVariables.TangentMatrix.Tangent * mTangentialVariables.Neighb_distance * custom_outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Tangent)  );
               const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
               Matrix  TangentMatrix = ZeroMatrix(dimension);

               TangentMatrix = rVariables.TangentMatrix.Normal * ( rVariables.Penalty.Normal / mTangentialVariables.Neighb_distance) * custom_outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Normal) ;
               TangentMatrix += rVariables.TangentMatrix.Tangent * custom_outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Tangent);

               TangentMatrix *= mTangentialVariables.Sign * rIntegrationWeight * mTangentialVariables.Neighb_distance ;

               rLeftHandSideMatrix -= TangentMatrix;

            }
            else {
               rLeftHandSideMatrix += rVariables.Penalty.Tangent * rIntegrationWeight * custom_outer_prod( rVariables.Surface.Tangent, rVariables.Surface.Tangent);
            }


         } */



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

            mTangentialVariables.IntegrationWeight  = rIntegrationWeight;

            this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
            this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

            // Compute Stress and Effective stress (if Exist)
            const array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
            array_1d< double, 3 >& ContactStress = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);
            ContactStress = ContactForce / mTangentialVariables.Neighb_distance / rIntegrationWeight ;

            if ( fabs(GetGeometry()[0].X() ) < 1e-5) 
               ContactStress = 0.0*ContactForce / mTangentialVariables.Neighb_distance / rIntegrationWeight ;

            if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE )) { 
               const array_1d< double, 3 >& EffectiveContactForce = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
               array_1d< double, 3 >& EffectiveContactStress = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_STRESS);
               EffectiveContactStress = EffectiveContactForce / mTangentialVariables.Neighb_distance / rIntegrationWeight ;
               if ( fabs(GetGeometry()[0].X() ) < 1e-5) 
                  EffectiveContactStress = 0.0*EffectiveContactForce / mTangentialVariables.Neighb_distance / rIntegrationWeight ;
            }

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

         if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE )) { 

            double EffectiveNormalModulus = CalculateEffectiveNormalForceModulus( NormalForceModulus*(-1) / rIntegrationWeight);
            EffectiveNormalModulus *= (-1) * rIntegrationWeight;
 
            array_1d<double, 3 >& EffectiveContactForce = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
            for (unsigned int j = 0; j < dimension; j++) {
               EffectiveContactForce[j] = EffectiveNormalModulus * rVariables.Surface.Normal[j];
            }
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

         mTangentialVariables.IntegrationWeight = rIntegrationWeight;
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


         if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE )) { 
            array_1d<double, 3 > & EffectiveContactForce = GetGeometry()[0].FastGetSolutionStepValue( EFFECTIVE_CONTACT_FORCE );
            for (unsigned int i = 0; i < dimension ; ++i) {
               EffectiveContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
            }
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
	/*

         mTangentialVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement);

         // 1. TRIAL TANGENT FORCE MODULUS
         double TrialTangentForceModulus = rVariables.Penalty.Tangent * rVariables.Gap.Tangent;
         TrialTangentForceModulus -= mTangentialVariables.PreviousTangentForceModulus / mTangentialVariables.IntegrationWeight;

         if ( fabs(TrialTangentForceModulus) < 1e-9)
            return 0.0;


         // 1a. Trial Tangent Stress Modulus
         // In order to use adhesion/cohesion ( that has units of Pressure and not force), first the forces are divided by IntegrationWeight ( 1 in PS or 2*pi*r in axisym) to get Stress units in the Contact Force
         // Also, the contact deppends on the Effective Normal Stress ( Normal Stress - WaterPressure )
         double TrialTangentStress = fabs( TrialTangentForceModulus / mTangentialVariables.Neighb_distance );

         // 1b. Effective NORMAL Stress
         double EffectiveNormalStress = CalculateEffectiveNormalForceModulus( rNormalForceModulus) / mTangentialVariables.Neighb_distance;


         if ( rVariables.Penalty.Tangent < 1e-6)
            return TrialTangentForceModulus;

         // 2. PERFORM THE RETURN MAPPING (First, prepare some variables)
         ContactFrictionLaw::FrictionLawVariables TangentVariables;

         TangentVariables.TangentPenalty = rVariables.Penalty.Tangent / mTangentialVariables.Neighb_distance;
         TangentVariables.FrictionCoefficient = mTangentialVariables.FrictionCoefficient;
         TangentVariables.Alpha = GetProperties()[CONTACT_ADHESION];
         TangentVariables.Adhesion = GetProperties()[CONTACT_ADHESION]; 

         TangentVariables.PlasticSlip = 0.0;
         TangentVariables.PlasticSlipOld = GetGeometry()[0].FastGetSolutionStepValue( CONTACT_PLASTIC_SLIP, 1 );

         mTangentialVariables.Slip  = mpFrictionLaw->EvaluateFrictionLaw(TrialTangentStress, EffectiveNormalStress, TangentVariables);

         // 3. PUT THE VARIABLES IN ITS PLACE
         mTangentialVariables.Sign = TrialTangentForceModulus / fabs( TrialTangentForceModulus);
         double TangentForceModulus = TrialTangentStress * mTangentialVariables.Neighb_distance * mTangentialVariables.Sign;


         if ( mTangentialVariables.Slip) {
            GetGeometry()[0].GetSolutionStepValue( CONTACT_PLASTIC_SLIP ) = TangentVariables.PlasticSlip;
         } else {
            GetGeometry()[0].GetSolutionStepValue( CONTACT_PLASTIC_SLIP ) = 0.0;
         }

         // 4. EVALUATE THE CONSTITUTIVE MATRIX
         mpFrictionLaw->EvaluateConstitutiveComponents( rVariables.TangentMatrix.Tangent, rVariables.TangentMatrix.Normal, TrialTangentStress, EffectiveNormalStress, TangentVariables);


         return TangentForceModulus;
	*/

      }



      //**************************** Check friction coefficient ***************************
      //***********************************************************************************

      double PointRigidContactPenalty3DCondition::CalculateEffectiveNormalForceModulus(const double& rNormalForceModulus)
      {

         double EffectiveForce = rNormalForceModulus;

         if ( GetGeometry()[0].HasDofFor(WATER_PRESSURE) )
         {
            double WaterForce = GetGeometry()[0].FastGetSolutionStepValue( WATER_PRESSURE );
            if (WaterForce > 0.0)
               WaterForce = 0.0;
            WaterForce *= mTangentialVariables.Neighb_distance;

            EffectiveForce = rNormalForceModulus - WaterForce ; // due to the sign convention
         }

         return EffectiveForce; 
      }


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



