//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_FRICTION_LAW_H_INCLUDED)
#define      KRATOS_FRICTION_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"


namespace Kratos
{
   ///@addtogroup ContactMechanicsApplication
   ///@{

   ///@name Kratos Globals
   ///@{

   ///@}
   ///@name Type Definitions
   ///@{

   ///@}
   ///@name  Enum's
   ///@{

   ///@}
   ///@name  Functions
   ///@{

   ///@}
   ///@name Kratos Classes
   ///@{

   /// Short class definition.
   /**
    * Base class of friction laws.
    */

   class KRATOS_API(CONTACT_MECHANICS_APPLICATION) FrictionLaw
   {
      public:

         struct FrictionLawVariables {
            double FrictionCoefficient;
            double Alpha;
            double Area;

            double TangentPenalty;

            double PlasticSlipOld;
            double PlasticSlip;
            double Adhesion;

            bool Implex;

            FrictionLawVariables() {Implex = false; Area = 1.0;}

            void Initialize(const double & rTangentPenalty, double PS, double & rArea, bool rImplex = false )
            {
               PlasticSlipOld = PS;
               PlasticSlip = PS;
               Area = rArea;
               TangentPenalty = rTangentPenalty / Area;
               Implex = rImplex;
            };

         };

         ///@name Type Definitions
         ///@{

         /// Pointer definition of FrictionLaw
         KRATOS_CLASS_POINTER_DEFINITION( FrictionLaw );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         FrictionLaw();

         /// Destructor.
         virtual ~FrictionLaw() {};


         /**
          * Clone function (has to be implemented by any derived class)
          * @return a pointer to a new instance of this friction law
          */
         virtual FrictionLaw::Pointer Clone() const;

         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{


         void FinalizeSolutionStep();

         // perform similar to a return mapping
         bool EvaluateFrictionLaw( double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables);

         void EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentForce, const double& rEffectiveNormalForce, FrictionLawVariables& rTangentVariables);

         double GetPlasticSlip() { return mPlasticSlip;};
         ///@}
         ///@name Access
         ///@{


         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         //virtual std::string Info() const;

         /// Print information about this object.
         //virtual void PrintInfo(std::ostream& rOStream) const;

         /// Print object's data.
         //virtual void PrintData(std::ostream& rOStream) const;


         ///@}
         ///@name Friends
         ///@{


         ///@}

      protected:
         ///@name Protected static Member Variables
         ///@{


         ///@}
         ///@name Protected member Variables
         ///@{


         ///@}
         ///@name Protected Operators
         ///@{


         ///@}
         ///@name Protected Operations
         ///@{

         virtual double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) {return 0; };

         virtual double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables ) {return 0; };

         virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables&  rTangentVariables ) {};

         ///@}
         ///@name Protected  Access
         ///@{


         ///@}
         ///@name Protected Inquiry
         ///@{


         ///@}
         ///@name Protected LifeCycle
         ///@{


         ///@}


      private:
         ///@name Static Member Variables
         ///@{


         ///@}
         ///@name Member Variables
         ///@{
         double mPlasticSlip;
         double mPlasticSlipNew;
         double mDeltaPlasticSlip;

         ///@}
         ///@name Private Operators
         ///@{


         ///@}
         ///@name Private Operations
         ///@{


         ///@}
         ///@name Private  Access
         ///@{


         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         /// Assignment operator.
         //FrictionLaw& operator=(FrictionLaw const& rOther);

         /// Copy constructor.
         //FrictionLaw(FrictionLaw const& rOther);

         ///@}
         ///@name Serialization
         ///@{

         friend class Serializer;

         virtual void save( Serializer& rSerializer ) const
         {
            rSerializer.save("mPlasticSlip",mPlasticSlip);
            rSerializer.save("mPlasticSlipNew",mPlasticSlipNew);
            rSerializer.save("mDeltaPlasticSlip",mDeltaPlasticSlip);
         }

         virtual void load( Serializer& rSerializer )
         {
            rSerializer.load("mPlasticSlip",mPlasticSlip);
            rSerializer.load("mPlasticSlipNew",mPlasticSlipNew);
            rSerializer.load("mDeltaPlasticSlip",mDeltaPlasticSlip);
         }

   }; // Class FrictionLaw

   ///@}

   ///@name Type Definitions
   ///@{

   /**
    * Definition of FrictionLaw variable
    */
   KRATOS_DEFINE_VARIABLE_IMPLEMENTATION( CONTACT_MECHANICS_APPLICATION, FrictionLaw::Pointer, FRICTION_LAW )

      ///@}
      ///@name Input and output
      ///@{


      /// input stream function
      // inline std::istream& operator >> (std::istream& rIStream,
      // 				    FrictionLaw& rThis);

      /// output stream function
      // inline std::ostream& operator << (std::ostream& rOStream,
      // 				    const FrictionLaw& rThis)
      //   {
      //     rThis.PrintInfo(rOStream);
      //     rOStream << std::endl;
      //     rThis.PrintData(rOStream);

      //     return rOStream;
      //   }
      ///@}

      ///@} addtogroup block

} // namespace Kratos

#endif // KRATOS_FRICTION_LAW_H_INCLUDED defined
