

#include "custom_conditions/custom_friction_laws/hardening_coulomb_friction_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
   HardeningCoulombFrictionLaw::HardeningCoulombFrictionLaw()
   {

   }

   HardeningCoulombFrictionLaw::~HardeningCoulombFrictionLaw()
   {

   }

   double HardeningCoulombFrictionLaw::EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double aux = std::exp( rTangentVariables.Alpha * rPlasticSlip );
      double H = - aux * rTangentVariables.FrictionCoefficient * rTangentVariables.Alpha * fabs(rNormalStress);
      return H;
   }

   double HardeningCoulombFrictionLaw::EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double aux = std::exp( rTangentVariables.Alpha * rPlasticSlip );
      double YieldFunction = fabs(rTangentStress) - rTangentVariables.FrictionCoefficient * aux * fabs(rNormalStress); 
      return YieldFunction;

   }

   void HardeningCoulombFrictionLaw::EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& rGamma, FrictionLawVariables& rTangentVariables) 
   {
      rdF_dt = 1.0;

      double aux = std::exp( rTangentVariables.Alpha * rGamma );
      rdF_dp = - rTangentVariables.FrictionCoefficient * aux;
   }
} // end namespace Kratos
