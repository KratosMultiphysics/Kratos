

#include "custom_conditions/custom_friction_laws/coulomb_adhesion_friction_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{
   CoulombAdhesionFrictionLaw::CoulombAdhesionFrictionLaw()
   {

   }

   CoulombAdhesionFrictionLaw::~CoulombAdhesionFrictionLaw()
   {

   }


   double CoulombAdhesionFrictionLaw::EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double YieldFunction = fabs(rTangentStress) - rTangentVariables.FrictionCoefficient * fabs(rNormalStress) - rTangentVariables.Adhesion; 
      return YieldFunction;

   }

   
   double CoulombAdhesionFrictionLaw::EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      return 0.0;
   }


   void CoulombAdhesionFrictionLaw::EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables& rTangentVariables ) 
   {

      rdF_dt = 1.0;
      rdF_dp = - rTangentVariables.FrictionCoefficient;
   }

} // end namespace Kratos
