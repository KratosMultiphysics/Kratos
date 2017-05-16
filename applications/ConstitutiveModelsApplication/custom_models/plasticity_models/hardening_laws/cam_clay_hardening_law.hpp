//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_laws/hardening_law.hpp"

namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
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
   /** Detail class definition.
    */
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CamClayHardeningLaw 
      : public HardeningLaw
      {
         protected:

            constexpr static std::size_t VarSize = 5;
         public:
    
            typedef InternalVariables<VarSize>   InternalVariablesType;
            typedef PlasticModelData<VarSize>          PlasticDataType;

            /// Pointer definition of CamClayHardeningLaw
            KRATOS_CLASS_POINTER_DEFINITION( CamClayHardeningLaw );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            CamClayHardeningLaw();

            /// Copy constructor.
            CamClayHardeningLaw(CamClayHardeningLaw const& rOther);

            /// Assignment operator.
            CamClayHardeningLaw& operator=(CamClayHardeningLaw const& rOther);

            /// Clone.
            virtual HardeningLaw::Pointer Clone() const override;

            /// Destructor.
            ~CamClayHardeningLaw();

            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{


            /**
             * Calculate Hardening functions
             */

            virtual double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening);

            /**
             * Calculate Hardening function derivatives
             */

            virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening);


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
            virtual std::string Info() const override
            {
               std::stringstream buffer;
               buffer << "CamClayHardeningLaw" ;
               return buffer.str();
            }

            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const override
            {
               rOStream << "CamClayHardeningLaw";
            }

            /// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const override
            {
               rOStream << "CamClayHardeningLaw Data";
            }    


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
            ///@name Serialization
            ///@{
            friend class Serializer;


            virtual void save(Serializer& rSerializer) const override
            {
               KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
            }

            virtual void load(Serializer& rSerializer) override
            {
               KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
            }

            ///@}
            ///@name Private Inquiry
            ///@{


            ///@}
            ///@name Un accessible methods
            ///@{


            ///@}

      }; // Class CamClayHardeningLaw

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CAM_CLAY_HARDENING_LAW_H_INCLUDED  defined 


