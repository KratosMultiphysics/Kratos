//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                 LHauser $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODIFIED_CAM_CLAY_PLASTIC_POTENTIAL_H_INCLUDED )
#define      KRATOS_MODIFIED_CAM_CLAY_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_utilities.hpp"

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
   template<class THardeningRule>
      class ModifiedCamClayPlasticPotential : public YieldSurface<THardeningRule>
   {
      public:

         ///@name Type Definitions
         ///@{

         typedef ConstitutiveModelData::MatrixType                          MatrixType;
         typedef ConstitutiveModelData::VectorType                          VectorType;
         typedef ConstitutiveModelData::ModelData                        ModelDataType;
         typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

         typedef YieldSurface<THardeningRule>                                 BaseType;
         typedef typename YieldSurface<THardeningRule>::Pointer        BaseTypePointer;
         typedef typename BaseType::PlasticDataType                    PlasticDataType;

         /// Pointer definition of ModifiedCamClayPlasticPotential
         KRATOS_CLASS_POINTER_DEFINITION( ModifiedCamClayPlasticPotential );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         ModifiedCamClayPlasticPotential() : BaseType() {}

         /// Copy constructor.
         ModifiedCamClayPlasticPotential(ModifiedCamClayPlasticPotential const& rOther) : BaseType(rOther) {}


         /// Assignment operator.
         ModifiedCamClayPlasticPotential& operator=(ModifiedCamClayPlasticPotential const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         BaseTypePointer Clone() const override
         {
            return Kratos::make_shared<ModifiedCamClayPlasticPotential>(*this);
         }

         /// Destructor.
         ~ModifiedCamClayPlasticPotential() override {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{


         //*************************************************************************************
         //*************************************************************************************
         // evaluation of the derivative of the yield surface respect the stresses
         VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition) override
         {
            KRATOS_TRY

            const ModelDataType & rModelData = rVariables.GetModelData();
            const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

            // Material Parameters
            const Properties& rProperties = rModelData.GetProperties();
            const double& rShearM = rProperties[CRITICAL_STATE_LINE];

            // compute something with the hardening rule


            double MeanStress, J2, LodeAngle;

            VectorType V1, V2;
            // more work is requiered
            StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
            StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

            double ThirdInvariantEffect = 1.0; // LMV: to be done
            double PreconsolidationStress;
            PreconsolidationStress = MeanStress + 3.0/MeanStress * pow( J2 / rShearM / ThirdInvariantEffect, 2);


            rDeltaStressYieldCondition  = ( 2.0*MeanStress - PreconsolidationStress) * V1 + 2.0 * 3.0 * pow( ThirdInvariantEffect / rShearM, 2) * J2 * V2;

            return rDeltaStressYieldCondition;

            KRATOS_CATCH(" ")
         }


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
         std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "ModifiedCamClayPlasticPotential" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "ModifiedCamClayPlasticPotential";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "ModifiedCamClayPlasticPotential Data";
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


         void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class ModifiedCamClayPlasticPotential

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_CAM_CLAY_PLASTIC_POTENTIAL_H_INCLUDED  defined
