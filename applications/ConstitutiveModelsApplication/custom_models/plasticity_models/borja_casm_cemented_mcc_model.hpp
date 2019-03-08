//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                   March 2019 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BORJA_CASM_CEMENTED_MCC_MODEL_H_INCLUDED )
#define  KRATOS_BORJA_CASM_CEMENTED_MCC_MODEL_H_INCLUDED

// System includes

// External includes
#include <iostream>
#include <fstream>

// Project includes
#include "custom_models/plasticity_models/casm_structured_soil_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/casm_cemented_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/casm_cemented_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
//#include "custom_models/elasticity_models/tamagnini_model.hpp"


//***** the hardening law associated to this Model has ... variables
// Variables
// 0. Plastic Multiplier
// 1. Volumetric Plastic Strain
// 2. Deviatoric Plastic Strain
// 3. Abs Value Volumetric Plastic Strain
// 4. p0 (preconsolidation pressure of debonded soil)
// 5. b  (bonding)
// 6. pc (preconsolidation pressure of bonded soil)
// 7. pt (tensile strength)
// 8. NonLocal Vol Plastic Strain
// 9. NonLocal Dev Plastic Strain
//10. NonLocal Abs Vol Plastic Strain

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) BorjaCasmCementedMccModel : public CasmStructuredSoilModel<BorjaModel, CasmCementedYieldSurface<CasmCementedHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         typedef BorjaModel                                     ElasticityModelType;
         //typedef TamagniniModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef CasmCementedHardeningRule                             HardeningRuleType;
         typedef CasmCementedYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef CasmStructuredSoilModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of BorjaCasmCementedMccModel
         KRATOS_CLASS_POINTER_DEFINITION( BorjaCasmCementedMccModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         BorjaCasmCementedMccModel() : BaseType() {
            ModifiedCamClayPlasticPotential<CasmCementedHardeningRule> Object;
            YieldSurface<CasmCementedHardeningRule>::Pointer pPlasticPotential = Object.Clone();
            std::cout << " thisPointer " << pPlasticPotential << std::endl;
            mYieldSurface = CasmCementedYieldSurface<CasmCementedHardeningRule>( pPlasticPotential);
            mInitialized = false;
         }

         /// Copy constructor.
         BorjaCasmCementedMccModel(BorjaCasmCementedMccModel const& rOther) : BaseType(rOther)  {}

         /// Assignment operator.
         BorjaCasmCementedMccModel& operator=(BorjaCasmCementedMccModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return Kratos::make_shared<BorjaCasmCementedMccModel>(*this);
         }

         /// Destructor.
         ~BorjaCasmCementedMccModel() override {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{
         

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
            buffer << "BorjaCasmCementedMccModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "BorjaCasmCementedMccModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "BorjaCasmCementedMccModel Data";
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

         virtual void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class BorjaCasmCementedMccModel

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

#endif // KRATOS_BORJA_CASM_CEMENTED_MCC_MODEL_H_INCLUDED  defined 


