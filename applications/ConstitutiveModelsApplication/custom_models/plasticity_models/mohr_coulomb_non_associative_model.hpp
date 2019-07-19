//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MOHR_COULOMB_NON_ASSOCIATIVE_MODEL_H_INCLUDED )
#define  KRATOS_MOHR_COULOMB_NON_ASSOCIATIVE_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/soil_non_associative_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/mohr_coulomb_v1_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/mohr_coulomb_non_associative_yield_surface.hpp"
#include "custom_models/elasticity_models/hencky_linear_model.hpp"

#include "custom_models/plasticity_models/yield_surfaces/plastic_potential/mohr_coulomb_plastic_potential.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MohrCoulombNonAssociativeModel : public SoilNonAssociativeModel<HenckyLinearModel, MohrCoulombNonAssociativeYieldSurface<MohrCoulombV1HardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         typedef HenckyLinearModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef MohrCoulombV1HardeningRule                             HardeningRuleType;
         typedef MohrCoulombNonAssociativeYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef SoilNonAssociativeModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of MohrCoulombNonAssociativeModel
         KRATOS_CLASS_POINTER_DEFINITION( MohrCoulombNonAssociativeModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         MohrCoulombNonAssociativeModel() : BaseType() {
            MohrCoulombPlasticPotential<MohrCoulombV1HardeningRule> Object;
            YieldSurface<MohrCoulombV1HardeningRule>::Pointer pPlasticPotential = Object.Clone();
            mYieldSurface = MohrCoulombNonAssociativeYieldSurface<MohrCoulombV1HardeningRule>(pPlasticPotential);
         }

         /// Copy constructor.
         MohrCoulombNonAssociativeModel(MohrCoulombNonAssociativeModel const& rOther) : BaseType(rOther) {}

         /// Assignment operator.
         MohrCoulombNonAssociativeModel& operator=(MohrCoulombNonAssociativeModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return Kratos::make_shared<MohrCoulombNonAssociativeModel>(*this);
         }

         /// Destructor.
         ~MohrCoulombNonAssociativeModel() override {}


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
            }

            KRATOS_CATCH("")
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
            else if (rThisVariable==DELTA_PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
            }
            else if (rThisVariable==PRE_CONSOLIDATION_STRESS)
            {
               rValue = this->mInternal.Variables[3];
            }
            else {
               rValue = NonAssociativePlasticityModel::GetValue( rThisVariable, rValue);
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
            buffer << "MohrCoulombNonAssociativeModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "MohrCoulombNonAssociativeModel";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "MohrCoulombNonAssociativeModel Data";
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

   }; // Class MohrCoulombNonAssociativeModel

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

#endif // KRATOS_MOHR_COULOMB_V1_MODEL_H_INCLUDED  defined
