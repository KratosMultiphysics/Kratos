
#include "custom_conditions/custom_friction_laws/friction_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{
   /*   FrictionLaw::FrictionLaw()
        {

        }

        FrictionLaw::~FrictionLaw()
        {

        } */ // defined in the other place

   bool ContactFrictionLaw::EvaluateFrictionLaw( double& rTangentStress, const double& rNormalStress, FrictionLawVariables& rTangentVariables )
   {
      bool Slip = false;

      double Gamma = rTangentVariables.PlasticSlipOld;

      double ContactYield = EvaluateContactYield( rTangentStress, rNormalStress, Gamma, rTangentVariables);

      if ( ContactYield < 1.0e-8)  // elasticPart
      {
         return Slip;
      }

      Slip = true;

      int i = 0;
      double Hardening, Derivative;
      double DeltaGamma = 0, DeltaDeltaGamma = 0; 
      double tangentStress = rTangentStress;
      while ( i < 100)
      {
         i = i + 1;

         Hardening = EvaluateHardening(rNormalStress, Gamma + DeltaGamma, rTangentVariables);  // this is the derivative of the yield stress respect the plastic slip
         Derivative = rTangentVariables.TangentPenalty - Hardening; // - some term due to the hardening law
         DeltaDeltaGamma = ContactYield / Derivative;
         DeltaGamma += DeltaDeltaGamma;

         tangentStress = rTangentStress - rTangentVariables.TangentPenalty * DeltaGamma; 

         ContactYield = EvaluateContactYield( tangentStress, rNormalStress,  Gamma+DeltaGamma, rTangentVariables);


         if ( fabs( ContactYield) < 1e-8)
            break;

      }
      if ( i > 90) {
         std::cout << " THIS contact DID NOT CONVERGE " << std::endl;
         std::cout << " TANGENT STRESS TRIAL " << rTangentStress << " tanS " << tangentStress << " YIELD " << ContactYield << std::endl;
         std::cout << " EffectiveNormal " << rNormalStress << std::endl;
      }
      rTangentStress = tangentStress;
      rTangentVariables.PlasticSlip = Gamma + DeltaGamma;

      return Slip;
   }

   void ContactFrictionLaw::EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentStress, const double& rNormalStress, FrictionLawVariables& rTangentVariables) 
   {
      double Hardening = EvaluateHardening( rNormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

      double dF_dt, dF_dp;  // derivatives of the yield function respect the tangential and normal contact stresses.
      EvaluateYieldDerivativeRespectStress( dF_dt, dF_dp, rTangentStress, rNormalStress, rTangentVariables.PlasticSlip, rTangentVariables);

      //double Auxiliar = -rTangentModulus / (rTangentModulus*dF_dt - Hardening);
      if ( rTangentVariables.TangentPenalty < 1E-8)
      {
         rNormalModulus = 0.0; rTangentModulus = 0.0;
         return;
      }
      double Auxiliar = -rTangentVariables.TangentPenalty / (rTangentVariables.TangentPenalty *dF_dt - Hardening);

      rNormalModulus = Auxiliar * dF_dp ; // * rTangentVariables.NormalPenalty;
      rTangentModulus = Auxiliar * Hardening;

   }

} // end namespace Kratos

