//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#include "custom_friction/friction_law.hpp"

namespace Kratos
{

   void FrictionLaw::InitializeSolutionStep()
   {
   }

   void FrictionLaw::FinalizeSolutionStep()
   {
      mDeltaPlasticSlip = mPlasticSlipNew - mPlasticSlip;
      mPlasticSlip = mPlasticSlipNew;
   }

   /**
    * Methods
    */
   bool FrictionLaw::EvaluateFrictionLaw( double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables )
   {

      mPlasticSlipNew = mPlasticSlip;
      if ( rTangentVariables.FrictionCoefficient == 0 && rTangentVariables.Adhesion == 0)
      {
         rTangentForce = 0.0;
         rTangentVariables.PlasticSlip = 0.0;
         rTangentVariables.PlasticSlipOld = 0.0;
         return true;
      }
      if ( rTangentVariables.TangentPenalty == 0) {
         return false;
      }

      double TangentStress = rTangentForce; //  / rTangentVariables.Area;
      double NormalStress = rNormalForce; // / rTangentVariables.Area;
      if ( rTangentVariables.Area > 0.0) {
         TangentStress = rTangentForce / rTangentVariables.Area;
         NormalStress = rNormalForce / rTangentVariables.Area;
      }

      bool Slip = false;

      if ( NormalStress < 0.0) {
         rTangentForce = 0.0;
         return Slip;
      }


      if ( rTangentVariables.Implex && rTangentVariables.PlasticSlipOld == 0)
      {
         // oye tu això està malament,em servia per el primer pas però no m'agrada pq fa que res funcioni
         /*rTangentForce = 0.0;
         rTangentVariables.PlasticSlip = 0.0;
         rTangentVariables.PlasticSlipOld = 0.0;
         return Slip;*/
      }

      double Gamma = rTangentVariables.PlasticSlipOld;
      double ContactYield = EvaluateContactYield( TangentStress, NormalStress, Gamma, rTangentVariables);


      if ( rTangentVariables.Implex) {
         TangentStress -= rTangentVariables.TangentPenalty * mDeltaPlasticSlip;
         rTangentForce = TangentStress * rTangentVariables.Area;
         rTangentVariables.PlasticSlip = Gamma + mDeltaPlasticSlip;
         Slip = true;
      }
      else {

         if ( ContactYield < 1.0e-8)  // elasticPart
         {
            return Slip;
         }

         Slip = true;
         int i = 0;
         double Hardening, Derivative;
         double DeltaGamma = 0, DeltaDeltaGamma = 0;
         double CurrentStress = TangentStress;
         while ( i < 100)
         {
            i = i + 1;

            Hardening = EvaluateHardening(NormalStress, Gamma + DeltaGamma, rTangentVariables);  // this is the derivative of the yield stress respect the plastic slip
            Derivative = rTangentVariables.TangentPenalty - Hardening; // - some term due to the hardening law
            DeltaDeltaGamma = ContactYield / Derivative;
            DeltaGamma += DeltaDeltaGamma;

            CurrentStress = TangentStress - rTangentVariables.TangentPenalty * DeltaGamma;

            ContactYield = EvaluateContactYield( CurrentStress, NormalStress,  Gamma+DeltaGamma, rTangentVariables);


            if ( fabs( ContactYield) < 1e-8)
               break;

         }
         if ( i > 90) {
            std::cout << " THIS contact DID NOT CONVERGE " << std::endl;
            std::cout << " TANGENT STRESS TRIAL " << TangentStress << " current tangent stress " << CurrentStress << " YIELD " << ContactYield << std::endl;
            std::cout << " AREA " << rTangentVariables.Area << " tangent force trial " << rTangentForce << std::endl;
            std::cout << " EffectiveNormal " << NormalStress << std::endl;
            std::cout << " PENALTY " << rTangentVariables.TangentPenalty << std::endl;
            std::cout << " gamma " << DeltaGamma << " , " << Derivative << std::endl;
            std::cout << " gamma " << Gamma << std::endl;
         }

         rTangentForce = CurrentStress * rTangentVariables.Area;
         rTangentVariables.PlasticSlip = Gamma + DeltaGamma;

      }
      mPlasticSlipNew = rTangentVariables.PlasticSlip;

      return Slip;
   }

   /**
    * Methods
    */
   void FrictionLaw::EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables)
   {
      if ( rTangentVariables.FrictionCoefficient == 0 && rTangentVariables.Adhesion == 0)
      {
         rNormalModulus = 0.0;
         rTangentModulus = 0.0;
         return;
      }
      if ( rTangentVariables.TangentPenalty == 0) {
         rNormalModulus = 0.0;
         rTangentModulus = 0.0;
         return;
      }

      if ( rNormalForce < 0.0 || rTangentVariables.TangentPenalty < 1e-8) {
         rNormalModulus = 0.0;
         rTangentModulus = 0.0;
         return;
      }

      double NormalStress = rNormalForce;
      double TangentStress = rTangentForce;
      if ( rTangentVariables.Area > 0.0) {
         NormalStress /= rTangentVariables.Area;
         TangentStress /= rTangentVariables.Area;
      }

      if ( rTangentVariables.Implex ) {
         if ( mDeltaPlasticSlip == 0)
         {
            /*rNormalModulus = 0.0;
            rTangentModulus = 0.0;
            return;*/
         }
         rNormalModulus = 0.0;
         rTangentModulus = 0.0;
         if ( rTangentForce > 0.0) {
            rTangentModulus = rTangentVariables.TangentPenalty * rTangentVariables.Area * mDeltaPlasticSlip /  rTangentForce;
         }
      } else {
         double Hardening = EvaluateHardening( NormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

         double dF_dt, dF_dp;  // derivatives of the yield function respect the tangential and normal contact stresses.
         EvaluateYieldDerivativeRespectStress( dF_dt, dF_dp, TangentStress, NormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

         //double Auxiliar = -rTangentVariables.TangentPenalty / (rTangentVariables.TangentPenalty *dF_dt - Hardening);

         //rNormalModulus = ( rTangentVariables.TangentPenalty / ( Hardening + rTangentVariables.TangentPenalty)) * dF_dp ;
         rNormalModulus = ( rTangentVariables.TangentPenalty / ( Hardening + rTangentVariables.TangentPenalty)) * dF_dp ;
         rTangentModulus = Hardening * rTangentVariables.TangentPenalty / (  Hardening + rTangentVariables.TangentPenalty );
      }

   }


}
// namespace Kratos
