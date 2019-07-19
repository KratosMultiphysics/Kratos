//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                 LHauser $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NONLOCAL_CASM_NUBIA_SOIL_MODEL_H_INCLUDED )
#define      KRATOS_NONLOCAL_CASM_NUBIA_SOIL_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/casm_base_soil_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/casm_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/casm_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"

#include "custom_models/plasticity_models/yield_surfaces/plastic_potential/nubia_plastic_potential.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonlocalCasmNubiaSoilModel : public CasmBaseSoilModel<BorjaModel, CasmYieldSurface<CasmHardeningRule> >
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


         /// Pointer definition of NonlocalCasmNubiaSoilModel
         KRATOS_CLASS_POINTER_DEFINITION( NonlocalCasmNubiaSoilModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         NonlocalCasmNubiaSoilModel() : BaseType() {
            NubiaPlasticPotential<CasmHardeningRule> Object;
            YieldSurface<CasmHardeningRule>::Pointer pPlasticPotential = Object.Clone();
            mYieldSurface = CasmYieldSurface<CasmHardeningRule>( pPlasticPotential);
         }

         /// Copy constructor.
         NonlocalCasmNubiaSoilModel(NonlocalCasmNubiaSoilModel const& rOther) : BaseType(rOther) {}

         /// Assignment operator.
         NonlocalCasmNubiaSoilModel& operator=(NonlocalCasmNubiaSoilModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return Kratos::make_shared<NonlocalCasmNubiaSoilModel>(*this);
         }

         /// Destructor.
         ~NonlocalCasmNubiaSoilModel() override {}


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
            else if (rThisVariable==NONLOCAL_PLASTIC_VOL_DEF)
            {
               rValue = this->mInternal.Variables[8];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_VOL_DEF_ABS) {
               rValue = std::fabs( this->mInternal.Variables[8] );
            }
            else {
               rValue = CasmBaseSoilModel::GetValue( rThisVariable, rValue);
            }
            return rValue;

            KRATOS_CATCH("")
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
               mInternal.Variables[8] = rValue;
            } else {
               BaseType::SetValue( rVariable, rValue, rCurrentProcessInfo);
            }

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
            buffer << "NonlocalCasmNubiaSoilModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "NonlocalCasmNubiaSoilModel";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "NonlocalCasmNubiaSoilModel Data";
         }


         ///@}
         ///@name Friends
         ///@{


         ///@}

      protected:
         ///@name Protected static Member Variables
         ///@{

            // Calculate Stress and constitutive tensor
            void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               { // modify the internal variables to make it nonlocal
               }

               double LocalPlasticVolStrain = mInternal.Variables[1];
               double NonLocalPlasticVolStrain = mInternal.Variables[8];
               mInternal.Variables[1] = mInternal.Variables[8];

               // integrate "analytically" ps and pt from plastic variables. Then update the internal variables.

               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               const ModelDataType & rModelData = Variables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();

               const double & rInitialPrecon = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
               const double & rSwellingSlope = rMaterialProperties[SWELLING_SLOPE];
               const double & rOtherSlope = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];

               const double & rPlasticVolDef = Variables.Internal.Variables[1]; 
               double & rPreconsolidation = mInternal.Variables[5];

               rPreconsolidation = -rInitialPrecon * std::exp( -rPlasticVolDef/ ( rOtherSlope-rSwellingSlope) );
               CasmBaseSoilModel::CalculateStressAndConstitutiveTensors( rValues, rStressMatrix, rConstitutiveMatrix);

               if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) ) {
                  mPreviousInternal.Variables[8] = mInternal.Variables[8];
                  mInternal.Variables[8] = mInternal.Variables[1];
                  mInternal.Variables[1] = LocalPlasticVolStrain + ( mInternal.Variables[1] -  NonLocalPlasticVolStrain);
               } else {
                  mInternal.Variables[1] = LocalPlasticVolStrain;
               }



               KRATOS_CATCH("")
            }

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

   }; // Class NonlocalCasmNubiaSoilModel

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

#endif // KRATOS_NONLOCAL_CASM_NUBIA_SOIL_MODEL_H_INCLUDED  defined
