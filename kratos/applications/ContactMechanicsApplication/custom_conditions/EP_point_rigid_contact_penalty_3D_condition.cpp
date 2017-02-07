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
#include "custom_conditions/EP_point_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"

#include "custom_friction/friction_law.hpp"
#include "utilities/math_utils.h"
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
      return Condition::Pointer(new EPPointRigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
   }

   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer EPPointRigidContactPenalty3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
   {
      EPPointRigidContactPenalty3DCondition NewCondition( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
      NewCondition.mCurrentInfo = this->mCurrentInfo; 
      NewCondition.mSavedInfo   = this->mSavedInfo; 

      return Condition::Pointer( new EPPointRigidContactPenalty3DCondition( NewCondition)  ); 

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

      PointRigidContactCondition::InitializeNonLinearIteration(rCurrentProcessInfo);

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

      KRATOS_CATCH( "" )
   }
   //************************************************************************************
   //************************************************************************************

   void EPPointRigidContactPenalty3DCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      PointRigidContactPenalty3DCondition::FinalizeSolutionStep( rCurrentProcessInfo);

      mpFrictionLaw->FinalizeSolutionStep();

      mSavedInfo = mCurrentInfo; 

      KRATOS_CATCH( "" )
   }



   //************* Tangent Contact Force constitutive matrix      **********************
   //***********************************************************************************
   void EPPointRigidContactPenalty3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
   {
      KRATOS_TRY

      double NormalForceModulus = 0;
      NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

      ConstitutiveVariables ConstVariables;
      Vector AuxVector;

      rVariables.Slip = this->CalculateFrictionLaw( rVariables, ConstVariables, AuxVector);
      double TangentForceModulus = norm_2( AuxVector);

      double Area = this->CalculateSomeSortOfArea();

      if( fabs(TangentForceModulus) >= 1e-25 ){

         if( rVariables.Slip ){

            noalias(rLeftHandSideMatrix) -= ConstVariables.TangentTangentMatrix * outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) * rIntegrationWeight * Area;
            noalias(rLeftHandSideMatrix) -= ConstVariables.NormalTangentMatrix * outer_prod( ConstVariables.ForceDirection, rVariables.Surface.Normal) * rVariables.Penalty.Normal * rIntegrationWeight;
            noalias(rLeftHandSideMatrix) += ConstVariables.TangentForceRatio * rVariables.Penalty.Tangent * ( IdentityMatrix(3,3) - outer_prod( rVariables.Surface.Normal, rVariables.Surface.Normal) - outer_prod( ConstVariables.ForceDirection, ConstVariables.ForceDirection) ) * rIntegrationWeight;

         }
         else {


            noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent); // aquest terme jo no el veig

            noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * ( IdentityMatrix(3,3) - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

         }

      }

      KRATOS_CATCH( "" )
   }

   //**************************** Calculate Tangent Contact Force **********************
   //***********************************************************************************
   void EPPointRigidContactPenalty3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY


      ConstitutiveVariables ConstVariables;
      Vector TangentForceVector;

      this->CalculateFrictionLaw( rVariables, ConstVariables, TangentForceVector);
      TangentForceVector *= rIntegrationWeight; 



      GetGeometry()[0].SetLock();


      if ( GetGeometry()[0].SolutionStepsDataHas( CONTACT_FORCE) )  {
         array_1d<double, 3 > & ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
         array_1d<double, 3 > & ContactDev = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_STRESS);
         const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

         for(unsigned int i = 0; i < dimension; i++) {
            rRightHandSideVector[i] += TangentForceVector(i);
            ContactForce[i] += TangentForceVector(i);
            ContactDev[i] = TangentForceVector(i);
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
   bool EPPointRigidContactPenalty3DCondition::CalculateFrictionLaw(GeneralVariables & rVariables, ConstitutiveVariables & rConstitutiveVariables, Vector & rTangentForceVector)
   {
      KRATOS_TRY

      // 0. Resize rTangentVector
      rTangentForceVector = ZeroVector(3);

      // 1. Compute the deformation gradient and the push forward of the previous force
      Vector PreviousNormal  = mSavedInfo.n;
      Vector PreviousT1 = mSavedInfo.t1; 
      Vector PreviousT2 = mSavedInfo.t2; 

      Vector Normal = rVariables.Surface.Normal; 
      Vector T1 = rVariables.Surface.Tangent; 
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

      ForceMatrix = prod( ForceMatrix, trans(F));
      ForceMatrix = prod( F, ForceMatrix);

      Vector PreviousForce = prod( ForceMatrix, Normal);
      PreviousForce -=  inner_prod( PreviousForce, Normal) * Normal;

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
      double EffectiveNormalForce = CalculateEffectiveNormalForceModulus( NormalForceModulus, Area); // falta posar que no sigui del signe contrari. 


      // FrictionVariables & contact constitutive parameters
      FrictionLaw::FrictionLawVariables FrictionVariables;
      FrictionVariables.Initialize( rVariables.Penalty.Tangent, mpFrictionLaw->GetPlasticSlip(), Area);

      // ConstitutiveParameters
      FrictionVariables.FrictionCoefficient = GetProperties()[MU_DYNAMIC]; 
      FrictionVariables.Alpha = 0; 
      FrictionVariables.Adhesion = 0;


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
      mCurrentInfo.t1 = rVariables.Surface.Tangent;
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

      WeakPointerVector<Element >& rNeighbourElements = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);


      std::vector< double > AreaVector;
      for ( unsigned int el = 0; el < rNeighbourElements.size() ; el++) {

         const Geometry< Node < 3 > > & rElemGeom = rNeighbourElements[el].GetGeometry();
         unsigned int nBoundary = 0;

         std::vector< unsigned int > BoundaryNodes; 
         for ( unsigned int i = 0; i < rElemGeom.size(); i++) {

            if ( rElemGeom[i].Is(BOUNDARY) ) {
               BoundaryNodes.push_back( i );
               nBoundary += 1;
               if ( nBoundary == 3)
                  continue;
            }
         }

         if ( nBoundary == 3)
         {
            array_1d< double, 3 > Vector1 = rElemGeom[ BoundaryNodes[1] ].Coordinates() - rElemGeom[ BoundaryNodes[0] ].Coordinates();
            array_1d< double, 3 > Vector2 = rElemGeom[ BoundaryNodes[2] ].Coordinates() - rElemGeom[ BoundaryNodes[0] ].Coordinates();
            array_1d< double, 3 > Cross = MathUtils<double>::CrossProduct( Vector1, Vector2 ) / 2.0;
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
      double EffectiveForce2; 

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




} // Namespace Kratos



