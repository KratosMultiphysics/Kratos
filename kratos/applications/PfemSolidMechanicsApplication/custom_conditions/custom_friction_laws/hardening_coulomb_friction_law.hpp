//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_HARDENING_COULOMB_FRICTION_LAW_H_INCLUDED)
#define      KRATOS_HARDENING_COULOMB_FRICTION_LAW_H_INCLUDED

#include "custom_conditions/custom_friction_laws/coulomb_adhesion_friction_law.hpp"
#include "custom_conditions/custom_friction_laws/friction_law.hpp"
namespace Kratos
{

   class HardeningCoulombFrictionLaw
      : public CoulombAdhesionFrictionLaw
   {

      public:

         KRATOS_CLASS_POINTER_DEFINITION( HardeningCoulombFrictionLaw );

         //constructor
         HardeningCoulombFrictionLaw();

         virtual ~HardeningCoulombFrictionLaw();

      protected:


         virtual double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) ;

         virtual double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) ;

         virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables& rTangentVariables) ;

   };

}



#endif // define KRATOS_COULOMB_FRICTION_LAW_H_INCLUDED
