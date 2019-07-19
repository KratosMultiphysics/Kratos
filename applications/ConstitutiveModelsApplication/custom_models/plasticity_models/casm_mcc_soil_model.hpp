//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                 LHauser $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_MCC_SOIL_MODEL_H_INCLUDED )
#define      KRATOS_CASM_MCC_SOIL_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/casm_base_soil_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/casm_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/casm_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"

#include "custom_models/plasticity_models/yield_surfaces/plastic_potential/modified_cam_clay_plastic_potential.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CasmMCCSoilModel : public CasmBaseSoilModel<BorjaModel, CasmYieldSurface<CasmHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         typedef BorjaModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef CasmHardeningRule                             HardeningRuleType;
         typedef CasmYieldSurface<HardeningRuleType>           YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                     YieldSurfacePointer;

         //base type
         typedef CasmBaseSoilModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of CasmMCCSoilModel
         KRATOS_CLASS_POINTER_DEFINITION( CasmMCCSoilModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         CasmMCCSoilModel() : BaseType() {
            ModifiedCamClayPlasticPotential<CasmHardeningRule> Object;
            YieldSurface<CasmHardeningRule>::Pointer pPlasticPotential = Object.Clone();
            mYieldSurface = CasmYieldSurface<CasmHardeningRule>( pPlasticPotential);
         }

         /// Copy constructor.
         CasmMCCSoilModel(CasmMCCSoilModel const& rOther) : BaseType(rOther) {}

         /// Assignment operator.
         CasmMCCSoilModel& operator=(CasmMCCSoilModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return Kratos::make_shared<CasmMCCSoilModel>(*this);
         }

         /// Destructor.
         ~CasmMCCSoilModel() override {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{

         /**
          * Check
          */
         int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override
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
          * Set Values
          */
         void SetValue(const Variable<double>& rVariable,
               const double& rValue,
               const ProcessInfo& rCurrentProcessInfo) override 
         {
            KRATOS_TRY

            if ( rVariable == NONLOCAL_PLASTIC_VOL_DEF) {
               mInternal.Variables[4] = rValue;
            } else {
               BaseType::SetValue( rVariable, rValue, rCurrentProcessInfo);
            }

            KRATOS_CATCH("")
         }

         /**
          * Get Values
          */
         double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
         {

            KRATOS_TRY

            rValue=0;

            if (rThisVariable==PLASTIC_VOL_DEF)
            {
               rValue = this->mInternal.Variables[1];
            }
            else {
               rValue = CasmBaseSoilModel::GetValue( rThisVariable, rValue);
            }
            return rValue;

            KRATOS_CATCH("")
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
            buffer << "CasmMCCSoilModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "CasmMCCSoilModel";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "CasmMCCSoilModel Data";
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

   }; // Class CasmMCCSoilModel

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

#endif // KRATOS_CASM_MCC_SOIL_MODEL_H_INCLUDED  defined
