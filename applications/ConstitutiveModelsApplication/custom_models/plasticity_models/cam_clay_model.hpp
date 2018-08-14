//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CAM_CLAY_MODEL_H_INCLUDED )
#define  KRATOS_CAM_CLAY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/cam_clay_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/modified_cam_clay_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CamClayModel : public NonAssociativePlasticityModel<BorjaModel, ModifiedCamClayYieldSurface<CamClayHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         typedef BorjaModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef CamClayHardeningRule                             HardeningRuleType;
         typedef ModifiedCamClayYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef NonAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of CamClayModel
         KRATOS_CLASS_POINTER_DEFINITION( CamClayModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         CamClayModel() : BaseType() {}

         /// Copy constructor.
         CamClayModel(CamClayModel const& rOther) : BaseType(rOther) {}

         /// Assignment operator.
         CamClayModel& operator=(CamClayModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return Kratos::make_shared<CamClayModel>(*this);
         }

         /// Destructor.
         ~CamClayModel() override {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{

         /**
          * Check
          */
         int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
         {
            KRATOS_TRY

            //LMV: to be implemented. but should not enter in the base one

            return 0;

            KRATOS_CATCH("")
         }

         ///@}
         ///@name Access
         ///@{

         /**
          * Has Values
          */
         bool Has(const Variable<double>& rThisVariable) override
         {
            if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
               return true;

            return false;
         }


         /**
          * Get Values
          */
         double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
         {

            rValue=0;

            if (rThisVariable==PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0];
            }


            if (rThisVariable==DELTA_PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
            }


            return rValue;
         }

         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "CamClayModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "CamClayModel";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "CamClayModel Data";
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
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Serialization
         ///@{
         friend class Serializer;

         void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class CamClayModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_CAM_CLAY_MODEL_H_INCLUDED  defined
