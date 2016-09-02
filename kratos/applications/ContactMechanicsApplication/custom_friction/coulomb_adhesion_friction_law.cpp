//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#include "custom_friction/coulomb_adhesion_friction_law.hpp"

namespace Kratos
{

  /**
   * Constructor.
   */
   CoulombAdhesionFrictionLaw::CoulombAdhesionFrictionLaw()
   {
   }


  /**
   * Destructor.
   */
   CoulombAdhesionFrictionLaw::~CoulombAdhesionFrictionLaw()
   {
   }


  /**
   * Clone function (has to be implemented by any derived class)
   * @return a pointer to a new instance of this constitutive law
   * NOTE: implementation scheme:
   *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
   *      return p_clone;
   */
  FrictionLaw::Pointer CoulombAdhesionFrictionLaw::Clone() const
  {
    CoulombAdhesionFrictionLaw::Pointer p_clone(new CoulombAdhesionFrictionLaw(*this));
    return p_clone;
  }

  /**
   * Methods
   */  
  double CoulombAdhesionFrictionLaw::EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      return 0.0;
   }

  /**
   * Methods
   */
   double CoulombAdhesionFrictionLaw::EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) 
   {
      double YieldFunction = fabs(rTangentStress) - rTangentVariables.FrictionCoefficient * fabs(rNormalStress) - rTangentVariables.Adhesion; 
      return YieldFunction;

   }

  /**
   * Methods
   */
   void CoulombAdhesionFrictionLaw::EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables& rTangentVariables ) 
   {

      rdF_dt = 1.0;
      rdF_dp = - rTangentVariables.FrictionCoefficient;
   }


} // namespace Kratos
