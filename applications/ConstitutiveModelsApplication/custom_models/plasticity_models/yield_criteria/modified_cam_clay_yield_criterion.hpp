//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"

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
   template<class THardeningLaw>
      class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ModifiedCamClayYieldCriterion
      : public YieldCriterion<THardeningLaw>
      {    
         public:

            ///@name Type Definitions
            ///@{

            typedef ConstitutiveModelData::MatrixType                          MatrixType;
            typedef ConstitutiveModelData::VectorType                          VectorType;
            typedef ConstitutiveModelData::ModelData                        ModelDataType;
            typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

            typedef YieldCriterion<THardeningLaw>                                BaseType;
            typedef typename BaseType::Pointer                            BaseTypePointer;
            typedef typename BaseType::PlasticDataType                    PlasticDataType;

            /// Pointer definition of ModifiedCamClayYieldCriterion
            KRATOS_CLASS_POINTER_DEFINITION( ModifiedCamClayYieldCriterion );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            ModifiedCamClayYieldCriterion() : BaseType() {}

            /// Copy constructor.
            ModifiedCamClayYieldCriterion(ModifiedCamClayYieldCriterion const& rOther) : BaseType(rOther) {}


            /// Assignment operator.
            ModifiedCamClayYieldCriterion& operator=(ModifiedCamClayYieldCriterion const& rOther)
            {
               BaseType::operator=(rOther);
               return *this;
            }

            /// Clone.
            virtual BaseTypePointer Clone() const override
            {
               return (ModifiedCamClayYieldCriterion::Pointer(new ModifiedCamClayYieldCriterion(*this)));
            }

            /// Destructor.
            virtual ~ModifiedCamClayYieldCriterion() {}


            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{

            /**
             * Calculate Yield Condition
             */

            virtual double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition)
            {
               KRATOS_TRY

               // Material Parameters
               const double ShearM = 1.0;

               // compute something with the hardening law
               double PreconsolidationStress = -70.0;

               const ModelDataType & rModelData = rVariables.GetModelData();
               const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

               double MeanStress, LodeAngle;
               double DeviatoricQ; // == sqrt(3)*J2
      
               // more work is requiered
               StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, DeviatoricQ, LodeAngle);
               DeviatoricQ *= sqrt(3.0);


               rYieldCondition  = pow( DeviatoricQ/ShearM, 2);
               rYieldCondition += (MeanStress * (MeanStress - PreconsolidationStress) );


               return rYieldCondition;

               KRATOS_CATCH(" ")
            }

            //*************************************************************************************
            //*************************************************************************************
            // evaluation of the derivative of the yield surface respect the stresses
            virtual void CalculateYieldSurfaceDerivative(const PlasticDataType& rVariables, VectorType & rYieldSurfaceDerivative)
            {
               KRATOS_TRY

               // Material Parameters
               const double ShearM = 1.0;

               // compute something with the hardening law
               double PreconsolidationStress = -70.0;

               const ModelDataType & rModelData = rVariables.GetModelData();
               const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

               double MeanStress, J2, LodeAngle;
     
               VectorType V1, V2;
               // more work is requiered
               StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
               StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

               rYieldSurfaceDerivative  = ( 2.0*MeanStress - PreconsolidationStress) * V1 + 2.0 * 3.0 * pow( 1.0 / ShearM, 2) * J2 * V2;


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
            virtual std::string Info() const
            {
               std::stringstream buffer;
               buffer << "ModifiedCamClayYieldCriterion" ;
               return buffer.str();
            }

            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const
            {
               rOStream << "ModifiedCamClayYieldCriterion";
            }

            /// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const
            {
               rOStream << "ModifiedCamClayYieldCriterion Data";
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


            virtual void save(Serializer& rSerializer) const
            {
               KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
            }

            virtual void load(Serializer& rSerializer)
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

      }; // Class ModifiedCamClayYieldCriterion

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED  defined 


