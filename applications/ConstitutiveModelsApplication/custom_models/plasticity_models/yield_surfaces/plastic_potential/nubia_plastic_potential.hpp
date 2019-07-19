//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                 LHauser $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NUBIA_PLASTIC_POTENTIAL_H_INCLUDED )
#define      KRATOS_NUBIA_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
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
   template<class THardeningRule>
      class NubiaPlasticPotential : public YieldSurface<THardeningRule>
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

         /// Pointer definition of NubiaPlasticPotential
         KRATOS_CLASS_POINTER_DEFINITION( NubiaPlasticPotential );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         NubiaPlasticPotential() : BaseType() {}

         /// Copy constructor.
         NubiaPlasticPotential(NubiaPlasticPotential const& rOther) : BaseType(rOther) {}


         /// Assignment operator.
         NubiaPlasticPotential& operator=(NubiaPlasticPotential const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         BaseTypePointer Clone() const override
         {
            return Kratos::make_shared<NubiaPlasticPotential>(*this);
         }

         /// Destructor.
         ~NubiaPlasticPotential() override {}


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
            const double& rShearM   = rProperties[CRITICAL_STATE_LINE];
            //const double& rSpacingR = rProperties[SPACING_RATIO];
            const double& rShapeN   = rProperties[SHAPE_PARAMETER];
            const double & rSwellingSlope = rProperties[SWELLING_SLOPE];
            const double & rOtherSlope    = rProperties[NORMAL_COMPRESSION_SLOPE];

            //1- stress invariants and deviative vectors 

            double MeanStress, J2, LodeAngle;

            VectorType V1, V2;
            // more work is requiered
            StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
            StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

            // 2- Evaluate m of equation 2-12
            double m(0);
            if ( rProperties.Has(CASM_M) ) {
               m = rProperties[CASM_M];
            } 
            if ( m <= 1.0 ) {
               double BigLambda = 1.0 - rSwellingSlope/rOtherSlope;

               m = pow( rShearM*(6.0-rShearM), rShapeN) - pow( 3.0*rShearM, rShapeN);
               m /= BigLambda * (6.0 - rShearM) * pow( 3.0 * rShearM, rShapeN-1);
               m *= 2.0/3.0;
            }

            if ( m <= 1.0) {
               KRATOS_ERROR << " the given parameters are problematic with Nubia//Casm plastic flow " << std::endl;
            }

            double StressRatio = sqrt(3.0)*J2/(-MeanStress);


            double C1 = m * (m - 1.0) * rShapeN * pow( StressRatio/rShearM, rShapeN-1.0) * ( -1.0) * (-1.0) *  sqrt(3.0) * J2 / pow(MeanStress, 2) / rShearM;
            C1 /= 1.0 + ( m-1.0) * pow( StressRatio / rShearM, rShapeN);

            C1 += rShapeN * ( m - 1.0) / MeanStress;


            double C2 = m * (m - 1.0) * rShapeN * pow( StressRatio/rShearM, rShapeN-1.0) * sqrt(3.0)/(-MeanStress) / rShearM;
            C2 /= 1.0 + (m - 1.0) * pow( StressRatio/rShearM, rShapeN);

            rDeltaStressYieldCondition = C1 * V1 + C2 * V2;

            //std::cout << " P " << MeanStress << " J2 " << J2 << " stressRatio " << StressRatio << " m " << m << " C1 " << C1 << " C2 " << C2 << std::endl;

            if (false) {
               //double StressRatio = sqrt(3.0)*J2/(-MeanStress);
               double Ratio = ( pow(rShearM, rShapeN) - pow( StressRatio, rShapeN) ) / m / pow( StressRatio, rShapeN-1.0) ;
               double C22 = 1;
               double C12 = -sqrt(3.0)/3.0*Ratio * C22;

               if ( fabs(StressRatio) < 1e-12) {
                  C12 = -1.0;
                  C22 = 0.0;
               }
               std::cout << " P " << MeanStress << " J2 " << J2 << " stressRatio " << StressRatio << " m " << m << " C1 " << C1 << " C2 " << C2 << " C12 " << C12 << " C22 " << C22 <<  std::endl;
               std::cout << "                 ratio1 " << C1/C2 << " ratio2 " << C12/C22 <<  " FINAL RATIO " << C1*C22/C12/C2 << std::endl;
               rDeltaStressYieldCondition = C12 * V1 + C22 * V2;
            }

            if ( rProperties.Has(NORMALIZE_FLOW_RULE) ) {
               if ( rProperties[NORMALIZE_FLOW_RULE] > 0 ) {
                  double norm = 0.0;
                  for (unsigned int ii = 0; ii < 3; ii++)
                     norm += pow( rDeltaStressYieldCondition(ii), 2);
                  for (unsigned int ii = 3; ii < 6; ii++)
                     norm += 2.0 * pow( rDeltaStressYieldCondition(ii)/2.0, 2);
                  norm = sqrt(norm);
                  if ( norm > 1e-1)
                     rDeltaStressYieldCondition /= norm;

               }
            }

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
            buffer << "NubiaPlasticPotential" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "NubiaPlasticPotential";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "NubiaPlasticPotential Data";
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

   }; // Class NubiaPlasticPotential

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NUBIA_PLASTIC_POTENTIAL_H_INCLUDED  defined
