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

   /**
    * Constructor.
    */
   FrictionLaw::FrictionLaw()
   {
      mPlasticSlip = 0.0;
   }

   /**
    * Clone function (has to be implemented by any derived class)
    * @return a pointer to a new instance of this constitutive law
    * NOTE: implementation scheme:
    *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
    *      return p_clone;
    */
   FrictionLaw::Pointer FrictionLaw::Clone() const
   {
      KRATOS_THROW_ERROR(std::logic_error, "Called the virtual function for Clone", "");
   }

   void FrictionLaw::FinalizeSolutionStep()
   {
      mPlasticSlip = mPlasticSlipNew;
   }

   /**
    * Methods
    */
   bool FrictionLaw::EvaluateFrictionLaw( double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables )
   {
      double TangentStress = rTangentForce / rTangentVariables.Area;
      double NormalStress = rNormalForce / rTangentVariables.Area;

      bool Slip = false;

      if ( NormalStress < 0.0) {
         rTangentForce = 0.0;
         return false;
      }

      double Gamma = rTangentVariables.PlasticSlipOld;

      double ContactYield = EvaluateContactYield( TangentStress, NormalStress, Gamma, rTangentVariables);

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
         std::cout << " EffectiveNormal " << NormalStress << std::endl;
      }

      rTangentForce = CurrentStress * rTangentVariables.Area;
      rTangentVariables.PlasticSlip = Gamma + DeltaGamma;

      return Slip;
   }

   /**
    * Methods
    */
   void FrictionLaw::EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables) 
   {

      if ( rNormalForce < 0.0 || rTangentVariables.TangentPenalty < 1e-8) {
         rNormalModulus = 0.0;
         rTangentModulus = 0.0;
         return;
      }

      double NormalStress = rNormalForce / rTangentVariables.Area;
      double TangentStress = rTangentForce / rTangentVariables.Area;
      

      double Hardening = EvaluateHardening( NormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

      double dF_dt, dF_dp;  // derivatives of the yield function respect the tangential and normal contact stresses.
      EvaluateYieldDerivativeRespectStress( dF_dt, dF_dp, TangentStress, NormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

      //double Auxiliar = -rTangentVariables.TangentPenalty / (rTangentVariables.TangentPenalty *dF_dt - Hardening);

      //rNormalModulus = ( rTangentVariables.TangentPenalty / ( Hardening + rTangentVariables.TangentPenalty)) * dF_dp ;
      rNormalModulus = ( rTangentVariables.TangentPenalty / ( Hardening + rTangentVariables.TangentPenalty)) * dF_dp ;
      rTangentModulus = Hardening * rTangentVariables.TangentPenalty / (  Hardening + rTangentVariables.TangentPenalty );

   }


} // namespace Kratos

