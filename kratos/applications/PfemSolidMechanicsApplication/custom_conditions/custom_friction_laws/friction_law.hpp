


#if !defined(KRATOS_CONTACT_FRICTION_LAW_H_INCLUDED)
#define      KRATOS_CONTACT_FRICTION_LAW_H_INCLUDED
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

namespace Kratos
{

   class ContactFrictionLaw
   {
      public:

         struct FrictionLawVariables {
            double FrictionCoefficient;
            double Alpha; 

            double NormalPenalty;
            double TangentPenalty;

            double PlasticSlipOld;
            double PlasticSlip;
            double Adhesion;
         };

         KRATOS_CLASS_POINTER_DEFINITION( ContactFrictionLaw );

         // la Hardening Law no té més constructors, crec q pot estar suficient.

         // constructor
         ContactFrictionLaw() {};

         // destructor
         virtual ~ContactFrictionLaw() {};

         //bool EvaluateFrictionLaw( double& rTangentStress, const double& rNormalStress, const double& rFrictionCoefficient, const double& rPenaltyParameter); // perform similar to a return mappinig
         bool EvaluateFrictionLaw( double& rTangentStress, const double& rNormalStress, FrictionLawVariables& rTangentVariables);

         //void EvaluateYieldDerivativeMatrix( double& rNormalModulus, double & rTangentModulus, const double& rTangentStress, const double& rNormalStress, const double& rFriction, const double & rPN, const double& rPT); 
         //void EvaluateYieldDerivativeMatrix( double& rNormalModulus, double & rTangentModulus, const double& rTangentStress, const double& rNormalStress, FrictionLawVariables& rTangentVariables);
         void EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentStress, const double& rNormalStress, FrictionLawVariables& rTangentVariables);

      protected:
         virtual double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) {return 0; };

         virtual double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables ) {return 0; } ;

         //virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, const double& rFrictionCoefficient) {} ;
         virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables&  rTangentVariables ) {} ;




    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
    }

    virtual void load( Serializer& rSerializer )
    {
    }


   }; // end class

} // end Namespace Kratos

#endif // define KRATOS_CONTACT_FRICTION_LAW_H_INCLUDED
