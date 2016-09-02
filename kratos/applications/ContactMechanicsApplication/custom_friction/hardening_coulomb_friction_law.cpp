//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#include "custom_friction/hardening_coulomb_friction_law.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

  /**
   * Constructor.
   */
   HardeningCoulombFrictionLaw::HardeningCoulombFrictionLaw()
   {
   }

  /**
   * Destructor.
   */
   HardeningCoulombFrictionLaw::~HardeningCoulombFrictionLaw()
   {
   }

  /**
   * Clone function (has to be implemented by any derived class)
   * @return a pointer to a new instance of this constitutive law
   * NOTE: implementation scheme:
   *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
   *      return p_clone;
   */
  FrictionLaw::Pointer HardeningCoulombFrictionLaw::Clone() const
  {
    HardeningCoulombFrictionLaw::Pointer p_clone(new HardeningCoulombFrictionLaw(*this));
    return p_clone;
  }
  
  /**
   * Methods
   */
  double HardeningCoulombFrictionLaw::EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double aux = std::exp( rTangentVariables.Alpha * rPlasticSlip );
      double H = - aux * rTangentVariables.FrictionCoefficient * rTangentVariables.Alpha * fabs(rNormalStress);
      return H;
   }

  /**
   * Methods
   */
   double HardeningCoulombFrictionLaw::EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double aux = std::exp( rTangentVariables.Alpha * rPlasticSlip );
      double YieldFunction = fabs(rTangentStress) - rTangentVariables.FrictionCoefficient * aux * fabs(rNormalStress); 
      return YieldFunction;

   }

  /**
   * Methods
   */
   void HardeningCoulombFrictionLaw::EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& rGamma, FrictionLawVariables& rTangentVariables) 
   {
      rdF_dt = 1.0;

      double aux = std::exp( rTangentVariables.Alpha * rGamma );
      rdF_dp = - rTangentVariables.FrictionCoefficient * aux;
   }

} // end namespace Kratos
