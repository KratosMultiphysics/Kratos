//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                      LMV $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2017 $
//   Revision:            $Revision:                   -0.1 $
//
//

// System includes

// External includes

// Project includes

#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"
#include "includes/mat_variables.h"
namespace Kratos
{

   //************************************************************************************
   //************************************************************************************

   EPPointRigidContactPenalty3DCondition::EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
      : PointRigidContactPenalty3DCondition(NewId, pGeometry)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   EPPointRigidContactPenalty3DCondition::EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : PointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //************************************************************************************
   //************************************************************************************

   EPPointRigidContactPenalty3DCondition::EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
      : PointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   EPPointRigidContactPenalty3DCondition::EPPointRigidContactPenalty3DCondition( EPPointRigidContactPenalty3DCondition const& rOther )
      : PointRigidContactPenalty3DCondition(rOther),
      mCurrentInfo( rOther.mCurrentInfo),
      mSavedInfo( rOther.mSavedInfo)
   {
      //DO NOT ADD DOFS HERE!!!
   }

   //************************************************************************************
   //************************************************************************************

   Condition::Pointer EPPointRigidContactPenalty3DCondition::Create(IndexType NewId, const NodesArrayType& ThisNodes, PropertiesType::Pointer pProperties) const
   {
     return Kratos::make_shared<EPPointRigidContactPenalty3DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer EPPointRigidContactPenalty3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
      EPPointRigidContactPenalty3DCondition NewCondition( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
      NewCondition.mCurrentInfo = this->mCurrentInfo;
      NewCondition.mSavedInfo   = this->mSavedInfo;

      NewCondition.mpFrictionLaw = this->mpFrictionLaw->Clone();

      return Kratos::make_shared<EPPointRigidContactPenalty3DCondition>(NewCondition);

   }


   //************************************************************************************
   //************************************************************************************

   EPPointRigidContactPenalty3DCondition::~EPPointRigidContactPenalty3DCondition()
   {

   }

   //************************************************************************************
   //************************************************************************************
   void EPPointRigidContactPenalty3DCondition::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      //PointRigidContactCondition::InitializeNonLinearIteration(rCurrentProcessInfo);

      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_FORCE) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }
      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_STRESS) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_STRESS) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_STRESS);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      mImplex = false;
      if ( rCurrentProcessInfo.Has(IMPLEX)  == true ) {
         if ( rCurrentProcessInfo[IMPLEX] > 0.0 ) {
            if ( GetProperties().Has(IMPLEX_CONTACT) == true) {
               if ( GetProperties()[IMPLEX_CONTACT] > 0.0 ) {
                  mImplex = true;
               }
            } else {
               mImplex = true;
            }
         }
      }

      KRATOS_CATCH( "" )
   }
   //************************************************************************************
   //************************************************************************************

   void EPPointRigidContactPenalty3DCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      PointRigidContactCondition::InitializeSolutionStep( rCurrentProcessInfo);

      const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();
      const unsigned int voigt_size = dimension * (dimension +1) * 0.5;

      if(mSavedInfo.PreviousStepForceVector.size() != voigt_size){
         mSavedInfo.PreviousStepForceVector.resize(voigt_size);
         noalias(mSavedInfo.PreviousStepForceVector) = ZeroVector(voigt_size);
         mSavedInfo.n = ZeroVector(3);
         mSavedInfo.t1 = ZeroVector(3);
         mSavedInfo.t2 = ZeroVector(3);
      }

      mElasticYoungModulus = -1.0;

      KRATOS_CATCH( "" )
   }
   //************************************************************************************
   //************************************************************************************

   void EPPointRigidContactPenalty3DCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_FORCE) ) {
         array_1d<double, 3 >& CF = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
         for(unsigned int j = 0; j < 3; j++) {
            CF[j] = 0;
         }
      }
      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_STRESS) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_STRESS) ) {
         array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_STRESS);
         for(unsigned int j = 0; j < 3; j++) {
            CS[j] = 0;
         }
      }

      mImplex = false;
      // calculate the stress without implex
      ConditionVariables Variables;
      this->InitializeConditionVariables(Variables, rCurrentProcessInfo);

      this->CalculateKinematics(Variables, rCurrentProcessInfo, 0);
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      if ( Variables.Options.Is(ACTIVE)) {

         double IntegrationWeight = 1;
         IntegrationWeight = this->CalculateIntegrationWeight(IntegrationWeight);

         // not a nice way to put it.
         if ( Variables.CurrentRadius > 0.0)
            IntegrationWeight *= 2.0 * 3.141592654 * Variables.CurrentRadius;

         Vector RHSVector  = ZeroVector(dimension);
         ConstitutiveVariables ConstVariables;

         this->CalculateAndAddContactForces( RHSVector,  Variables, IntegrationWeight);
         // ( recompute all stresses and set it to the nodal variables)


      } else {
         // this is a bad-programing thing that allow to see if neighbour nodes are in contact. please, correct it
         if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_NORMAL) ) {
            array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_NORMAL);
            for(unsigned int j = 0; j < 3; j++) {
               CS[j] = 0;
            }
         }
      }
      if ( (fabs( GetGeometry()[0].X()) < 1e-5) && ( fabs(GetGeometry()[0].Z() ) < 1e-6) ) {
         if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_STRESS) ) {
            array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);
            for(unsigned int j = 0; j < 3; j++) {
               CS[j] = 0;
            }
         }
         if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_STRESS) ) {
            array_1d<double, 3 >& CS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_STRESS);
            for(unsigned int j = 0; j < 3; j++) {
               CS[j] = 0;
            }
         }
      }

      PointRigidContactCondition::FinalizeSolutionStep( rCurrentProcessInfo);
      mpFrictionLaw->FinalizeSolutionStep();

      mSavedInfo = mCurrentInfo;


      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************
   // i have to copy it because I save the value of the Young modulus of the continuum elements
   void EPPointRigidContactPenalty3DCondition::CalculateContactFactors(ConditionVariables &rVariables)
   {

      KRATOS_TRY

      //Compute the neighbour distance, then a stress-"like" may be computed.
      WeakPointerVector<Node<3> >& rN  = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);
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


      rVariables.ContributoryFactor = distance;

      //get contact properties and parameters
      double PenaltyParameter = 1;
      if( GetProperties().Has(PENALTY_PARAMETER) )
         PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

      WeakPointerVector<Element >& rE = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
      double ElasticModulus = 0;
      if( GetProperties().Has(YOUNG_MODULUS) )
         ElasticModulus = GetProperties()[YOUNG_MODULUS];
      else
         ElasticModulus = rE.front().GetProperties()[YOUNG_MODULUS];

      // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
      if (ElasticModulus <= 1.0e-5) {
         if ( mElasticYoungModulus < 0) {
            std::vector<double> ModulusVector;
            ProcessInfo SomeProcessInfo;
            for ( unsigned int i = 0; i < rE.size(); i++)
            {
               rE[i].CalculateOnIntegrationPoints(YOUNG_MODULUS, ModulusVector, SomeProcessInfo);
               ElasticModulus += ModulusVector[0];
            }
            ElasticModulus /= double(rE.size());
            mElasticYoungModulus = ElasticModulus;
         } else {
            ElasticModulus = mElasticYoungModulus;
         }
      }


      rVariables.Penalty.Normal  = 1e4 * PenaltyParameter * ElasticModulus;


      double ContributoryArea = this->CalculateSomeSortOfArea();


      //std::cout << " this->Id() " << this->Id() << " , " << rVariables.Penalty.Normal << " area " << ContributoryArea << " , " << Contact_Point[0] << " , " << rVariables.Penalty.Normal * ContributoryArea << std::endl;

      rVariables.Penalty.Normal *= ContributoryArea;


      double PenaltyRatio = 1;
      if( GetProperties().Has(TANGENTIAL_PENALTY_RATIO) )
         PenaltyRatio = GetProperties()[TANGENTIAL_PENALTY_RATIO];

      rVariables.Penalty.Tangent = rVariables.Penalty.Normal * PenaltyRatio ;

      //std::cout<<" ContactPoint["<<this->Id()<<"]: penalty_n"<<rVariables.Penalty.Normal<<", ElasticModulus: "<<ElasticModulus<<", distance: "<<distance<<std::endl;

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





   //************* Tangent Contact Force constitutive matrix      **********************
   //***********************************************************************************
   void EPPointRigidContactPenalty3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables, double& rIntegrationWeight)
   {
      KRATOS_TRY

      Matrix PreviousLHS = rLeftHandSideMatrix;
      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      ConstitutiveVariables ConstVariables;
      Vector AuxVector;

      rVariables.Slip = this->CalculateFrictionLaw( rVariables, ConstVariables, AuxVector);
      double TangentForceModulus = norm_2( AuxVector);

      double Area = this->CalculateSomeSortOfArea();


      if( fabs(TangentForceModulus) >= 1e-25 ){

        MatrixType Identity(3,3);
        noalias(Identity) = IdentityMatrix(3);

        if( rVariables.Slip ){
          if ( mImplex ) {
            noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * ( Identity - outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal) ) * rIntegrationWeight;

            noalias(rLeftHandSideMatrix) -= rVariables.Penalty.Tangent * ConstVariables.TangentTangentMatrix * ConstVariables.TangentForceRatio * ( Identity - outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal) - outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) ) * rIntegrationWeight;

          } else {


            noalias(rLeftHandSideMatrix) -= ConstVariables.TangentTangentMatrix * outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) * rIntegrationWeight * Area;
            noalias(rLeftHandSideMatrix) -= ConstVariables.NormalTangentMatrix * outer_prod( ConstVariables.ForceDirection, rVariables.Surface.Normal) * rVariables.Penalty.Normal * rIntegrationWeight;
            noalias(rLeftHandSideMatrix) += ConstVariables.TangentForceRatio * rVariables.Penalty.Tangent * ( Identity - outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal) - outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) ) * rIntegrationWeight;
          }

        }
        else {


          //noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent); // aquest terme jo no el veig

          noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * ( Identity - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

        }

      }

      KRATOS_CATCH( "" )
   }

   //**************************** Calculate Tangent Contact Force **********************
   //***********************************************************************************
   void EPPointRigidContactPenalty3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
         ConditionVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY


      ConstitutiveVariables ConstVariables;
      Vector TangentForceVector;

      this->CalculateFrictionLaw( rVariables, ConstVariables, TangentForceVector);
      TangentForceVector *= rIntegrationWeight;



      GetGeometry()[0].SetLock();

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      for(unsigned int i = 0; i < dimension; i++) {
         rRightHandSideVector[i] += TangentForceVector(i);
      }



      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_FORCE) )  {

         double Area = this->CalculateSomeSortOfArea();
         array_1d<double, 3 > & ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
         array_1d<double, 3 > & ContactStress = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);

         for(unsigned int i = 0; i < dimension; i++) {
            ContactForce[i] += TangentForceVector(i);

            ContactStress[i] = ContactForce[i] / (Area * rIntegrationWeight );

         }
      }
      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE) ) {

         double NormalForceModulus = 0;
         NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables);

         double Area = this->CalculateSomeSortOfArea();
         double EffectiveNormalForce = CalculateEffectiveNormalForceModulus( NormalForceModulus, Area);


         array_1d<double, 3 >& ECF = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
         for(unsigned int j = 0; j < 3; j++) {
            ECF[j] = ( EffectiveNormalForce * rVariables.Surface.Normal[j] * rIntegrationWeight ) + TangentForceVector[j];
         }
      }

      if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_STRESS) ) {
         double Area = this->CalculateSomeSortOfArea();
         array_1d<double, 3 >& ECF = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
         array_1d<double, 3 >& ECS = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_STRESS);
         for(unsigned int j = 0; j < 3; j++) {
            ECS[j] = ECF[j] / (Area*rIntegrationWeight);
         }
      }





      if ( rVariables.ContactStressVector.size() > 1) {
         double ForceModulus = norm_2( TangentForceVector);
         rVariables.ContactStressVector += MathUtils<double>::StressTensorToVector( ForceModulus * ( outer_prod(rVariables.Surface.Normal, ConstVariables.ForceDirection) + outer_prod( ConstVariables.ForceDirection, rVariables.Surface.Normal) ) , rVariables.ContactStressVector.size() );
      }

      GetGeometry()[0].UnSetLock();

      KRATOS_CATCH( "" )
   }


   // *********************** CalculateFrictionLaw *************************************
   //***********************************************************************************
   bool EPPointRigidContactPenalty3DCondition::CalculateFrictionLaw(ConditionVariables & rVariables, ConstitutiveVariables & rConstitutiveVariables, Vector & rTangentForceVector)
   {
      KRATOS_TRY

      // 0. Resize rTangentVector
      rTangentForceVector = ZeroVector(3);

      // 1. Compute the deformation gradient and the push forward of the previous force
      Vector PreviousNormal  = mSavedInfo.n;
      Vector PreviousT1 = mSavedInfo.t1;
      Vector PreviousT2 = mSavedInfo.t2;

      // if the bounding box does not have the parametric directions implemented use some sort of guess
      Vector Normal = rVariables.Surface.Normal;
      Vector T1 = rVariables.Surface.Tangent;

      if ( MathUtils<double>::Dot(T1, PreviousT1) < 0)
         T1 *= (-1.0);

      Vector T2 = ZeroVector(3);
      T2(0) =  Normal(1) * T1(2) - Normal(2)*T1(1);
      T2(1) = -Normal(0) * T1(2) + Normal(2)*T1(0);
      T2(2) =  Normal(0) * T1(1) - Normal(1)*T1(0);

      {
         // Get the real parametrization of the Box (if exists)
         SpatialBoundingBox::BoundingBoxParameters BoxParameters( this->GetGeometry()[0] ) ;
         this->mpRigidWall->GetParametricDirections( BoxParameters, T1, T2);
      }

      Matrix F = ZeroMatrix(3);
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            F(i,j) = T1(i)*PreviousT1(j) + T2(i)*PreviousT2(j) + Normal(i) * PreviousNormal(j);
         }
      }

      Matrix ForceMatrix = MathUtils<double>::StressVectorToTensor( mSavedInfo.PreviousStepForceVector );

      ForceMatrix = ConvertToTheAppropriateSize(ForceMatrix); // 2D and Axis are 2DMatrices

      ForceMatrix = prod( ForceMatrix, trans(F));
      ForceMatrix = prod( F, ForceMatrix);

      Vector PreviousForce = prod( ForceMatrix, Normal);
      PreviousForce -=  inner_prod( PreviousForce, Normal) * Normal;  // remove spurious normal

      // 2. Calculate the return mapping
      Vector TrialTangentForce = PreviousForce + rVariables.Penalty.Tangent * rVariables.Gap.Tangent * rVariables.Surface.Tangent;
      double ForceModulus   = norm_2( TrialTangentForce);

      Vector ForceDirection = TrialTangentForce;
      if ( ForceModulus > 1e-15)
         ForceDirection /= ForceModulus;

      // 2.a Prepare information
      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );


      double Area = this->CalculateSomeSortOfArea();
      double EffectiveNormalForce = CalculateEffectiveNormalForceModulus( NormalForceModulus, Area);
      if ( EffectiveNormalForce < 0.0)
         EffectiveNormalForce = 0.0;



      // FrictionVariables & contact constitutive parameters
      FrictionLaw::FrictionLawVariables FrictionVariables;
      FrictionVariables.Initialize( rVariables.Penalty.Tangent, mpFrictionLaw->GetPlasticSlip(), Area, mImplex);

      // ConstitutiveParameters
      FrictionVariables.FrictionCoefficient = GetProperties()[MU_DYNAMIC];
      if (FrictionVariables.FrictionCoefficient == 0)
      {
         if ( GetProperties().Has(CONTACT_FRICTION_ANGLE) )
         {
            FrictionVariables.FrictionCoefficient = GetProperties()[CONTACT_FRICTION_ANGLE];
            FrictionVariables.FrictionCoefficient = std::tan(3.14159 * FrictionVariables.FrictionCoefficient/180.0);
         }

      }
      FrictionVariables.Adhesion = 0;
      if ( GetProperties().Has(CONTACT_ADHESION) ) {
         FrictionVariables.Adhesion = GetProperties()[CONTACT_ADHESION];
      }
      if (  fabs(GetGeometry()[0].X()) < 1e-12)
      {
         FrictionVariables.FrictionCoefficient = 0.0;
      }
      FrictionVariables.Alpha = 0;
      //FrictionVariables.Adhesion = 0;
      if ( GetProperties().Has(VISCOSITY) )
      {
         FrictionVariables.Alpha = GetProperties()[VISCOSITY];
      }
      // no està bé, s'han de posar les coses al seu lloc


      // 2.b Perform the return mapping
      rVariables.Slip = mpFrictionLaw->EvaluateFrictionLaw( ForceModulus, EffectiveNormalForce, FrictionVariables);

      // 2.c Calculate the contact constitutive matrix, or whatever
      if ( rVariables.Slip) {
         mpFrictionLaw->EvaluateConstitutiveComponents( rConstitutiveVariables.NormalTangentMatrix, rConstitutiveVariables.TangentTangentMatrix, ForceModulus, EffectiveNormalForce, FrictionVariables);

         rConstitutiveVariables.TangentForceRatio = ForceModulus / norm_2( TrialTangentForce);
      }
      rConstitutiveVariables.ForceDirection = ForceDirection;

      // 3. Save geometrical information
      mCurrentInfo.PreviousStepForceVector = MathUtils<double>::StressTensorToVector( ForceModulus * (outer_prod( rVariables.Surface.Normal, ForceDirection) + outer_prod( ForceDirection, rVariables.Surface.Normal) ), mSavedInfo.PreviousStepForceVector.size() );

      mCurrentInfo.n = rVariables.Surface.Normal;
      mCurrentInfo.t1 = T1;
      mCurrentInfo.t2 = T2;


      rTangentForceVector = ForceModulus * ForceDirection;

      return rVariables.Slip;


      KRATOS_CATCH("")
   }



   //**************************** CalculateSomeSortOfArea ******************************
   //***********************************************************************************
   double EPPointRigidContactPenalty3DCondition::CalculateSomeSortOfArea()
   {

      double Area = 1;

      const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();
      if ( dimension != 3)
         return Area;


      WeakPointerVector<Element >& rNeighbourElements = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);


      std::vector< double > AreaVector;
      for ( unsigned int el = 0; el < rNeighbourElements.size() ; el++) {

         const Geometry< Node < 3 > > & rElemGeom = rNeighbourElements[el].GetGeometry();
         unsigned int nBoundary = 0;

         std::vector< unsigned int > BoundaryNodes;
         for ( unsigned int i = 0; i < rElemGeom.size(); i++) {

            if ( rElemGeom[i].Is(BOUNDARY) ) {
               const array_1d<double, 3> &  CN = rElemGeom[i].FastGetSolutionStepValue(CONTACT_NORMAL);

               if ( ( fabs(CN[0]) + fabs(CN[1]) + fabs(CN[2])  ) > 0.01) {
                  BoundaryNodes.push_back( i );
                  nBoundary += 1;
                  if ( nBoundary == 3)
                     break;
               }
            }
         }

         if ( nBoundary == 3)
         {
            array_1d< double, 3 > Vector1 = rElemGeom[ BoundaryNodes[1] ].Coordinates() - rElemGeom[ BoundaryNodes[0] ].Coordinates();
            array_1d< double, 3 > Vector2 = rElemGeom[ BoundaryNodes[2] ].Coordinates() - rElemGeom[ BoundaryNodes[0] ].Coordinates();
            array_1d< double, 3 > Cross;
            MathUtils<double>::CrossProduct( Cross, Vector1, Vector2 );
            Cross /= 2.0;
            double ThisArea = MathUtils<double>::Norm3( Cross);
            AreaVector.push_back(ThisArea);
         }

      }

      if ( AreaVector.size() > 0) {
         Area = 0;
         for (unsigned int i = 0; i < AreaVector.size() ; i++)
            Area += AreaVector[i];
         Area /= 3.0;
      }


      return Area;
   }

   //**************************** Evaluate normal effective force ***************************
   //****************************************************************************************
   double EPPointRigidContactPenalty3DCondition::CalculateEffectiveNormalForceModulus(const double& rNormalForceModulus, const double& rSomeArea) // tot això m'ho he de mirar
   {

      double EffectiveForce = rNormalForceModulus;
      //double EffectiveForce2;

      double WaterForce;
      if ( GetGeometry()[0].HasDofFor(WATER_PRESSURE) )
      {
         WaterForce = GetGeometry()[0].FastGetSolutionStepValue( WATER_PRESSURE );
         if (WaterForce > 0.0)
            WaterForce = 0.0;
         WaterForce *= rSomeArea;

         EffectiveForce = rNormalForceModulus + WaterForce; // since the normal is positive,....
      }
      else {
      }


      if ( false ) // this is the unbalanced normal load of the deformable body in the node. f = Int( B*sigma dV)
      {

         /*  std::vector< Vector > aux1, aux2, aux3;
             ProcessInfo SomeProcessInfo;
             WeakPointerVector<Element > rN = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
             Vector WaterForceVector = ZeroVector(2);
             Vector EffecForceVector = ZeroVector(2);
             Vector TotalForceVector = ZeroVector(2);
             for (unsigned int ne = 0; ne < rN.size() ; ne++)
             {
             rN[ne].CalculateOnIntegrationPoints( EFF_CON_WATER_FORCE, aux1, SomeProcessInfo);
             rN[ne].CalculateOnIntegrationPoints( EFF_CON_EFFEC_FORCE, aux2, SomeProcessInfo);
             rN[ne].CalculateOnIntegrationPoints( EFF_CON_TOTAL_FORCE, aux3, SomeProcessInfo);

             if (aux1[0].size() == 0)  {
         //std::cout << " in this exit " << std::endl;
         return EffectiveForce;
         }

         for (unsigned int se = 0; se < 3; se++)
         {
         if ( GetGeometry()[0].Id() == rN[ne].GetGeometry()[se].Id() )
         {
         WaterForceVector(0) += aux1[0][2*se];
         WaterForceVector(1) += aux1[0][2*se+1];
         EffecForceVector(0) += aux2[0][2*se];
         EffecForceVector(1) += aux2[0][2*se+1];
         TotalForceVector(0) += aux3[0][2*se];
         TotalForceVector(1) += aux3[0][2*se+1];
         }
         }

         }

         EffectiveForce2 = EffecForceVector(0) * rSurfaceNormal(0);
         EffectiveForce2 += EffecForceVector(1) * rSurfaceNormal(1);
         EffectiveForce2 /=  -mTangentialVariables.IntegrationWeight; */

         /*
            if ( false) {
            std::cout << " FINALLY: Node " << GetGeometry()[0].Id() << " has neigh: " << rN.size() << std::endl;
            std::cout << "     NormalForceModulus " << rNormalForceModulus << std::endl;
            std::cout << "     NodeWaterForce " << WaterForce <<  " and this pseudo-water-force " << -( WaterForceVector(0)*rSurfaceNormal(0) + WaterForceVector(1)*rSurfaceNormal(1) ) / mTangentialVariables.IntegrationWeight << " i and the vector is " << WaterForceVector <<  " divided by " << WaterForceVector  /  mTangentialVariables.IntegrationWeight << std::endl;
            std::cout << "     and this pseudo-effective-force " << -( EffecForceVector(0)*rSurfaceNormal(0) +EffecForceVector(1)*rSurfaceNormal(1) ) / mTangentialVariables.IntegrationWeight << " i and the vector is " << EffecForceVector << std::endl;
            std::cout << "     and this pseudo-TOTAL-force " << -( TotalForceVector(0)*rSurfaceNormal(0) +TotalForceVector(1)*rSurfaceNormal(1) ) / mTangentialVariables.IntegrationWeight << " i and the vector is " << TotalForceVector << std::endl;
            std::cout << "     EFFECIVE FORCE BEFORE : " << EffectiveForce << std::endl;
            std::cout << "     EFFECIVE FORCE NEWNEW : " << EffectiveForce2 << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl; }
            return EffectiveForce2; */
      }

      return EffectiveForce;

   }



   // ********************************** CONVERT THE CONTACT MATRIX TO 3D *************************
   // *********************************************************************************************
   Matrix  EPPointRigidContactPenalty3DCondition::ConvertToTheAppropriateSize(const Matrix & rForceMatrix)
   {

      if ( rForceMatrix.size1() != 3) {
         Matrix AuxMatrix = ZeroMatrix(3,3);

         for ( unsigned int i = 0; i < rForceMatrix.size1(); i++) {
            for ( unsigned int j = 0; j < rForceMatrix.size2(); j++) {
               AuxMatrix(i,j) = rForceMatrix(i,j);
            }
         }
         return AuxMatrix;
      }
      return rForceMatrix;

   }


} // Namespace Kratos
