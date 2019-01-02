//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_HENCKY_HYPER_ELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_HENCKY_HYPER_ELASTIC_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hyper_elastic_model.hpp"

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HenckyHyperElasticModel : public HyperElasticModel
   {
      public:

         ///@name Type Definitions
         ///@{

         /// Pointer definition of BorjaModel
         KRATOS_CLASS_POINTER_DEFINITION( HenckyHyperElasticModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         HenckyHyperElasticModel(): HyperElasticModel() {};

         /// Copy constructor.
         HenckyHyperElasticModel(HenckyHyperElasticModel const& rOther) : HyperElasticModel( rOther) {};

         /// Assignment operator.
         HenckyHyperElasticModel& operator=(HenckyHyperElasticModel const& rOther)
         {
            HyperElasticModel::operator=(rOther);
            return *this;
         };

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return ( HenckyHyperElasticModel::Pointer( new HenckyHyperElasticModel(*this)) );
         };


         /// Destructor.
         ~HenckyHyperElasticModel() override {};


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{

         /**
          * Initialize member data
          */
         void InitializeModel(ModelDataType& rValues) override
         {
            KRATOS_TRY

            HyperElasticModel::InitializeModel(rValues);

            // Compute trial strain
            //deformation gradient
            const MatrixType& rDeltaDeformationMatrix  = rValues.GetDeltaDeformationMatrix();

            //historical strain matrix
            rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);
            rValues.StrainMatrix = prod( rDeltaDeformationMatrix, rValues.StrainMatrix);
            rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(rDeltaDeformationMatrix));

            KRATOS_CATCH("")

         };

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
            buffer << "HenckyHyperElasticModel";
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "HenckyHyperElasticModel";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "HenckyHyperElasticModel Data";
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

         void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables) override
         {
            KRATOS_TRY

            rVariables.SetModelData(rValues);
            rVariables.SetState(rValues.State);

            const StressMeasureType& rStressMeasure  = rValues.GetStressMeasure();

            if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C=FT*F)  C^-1=(FT*F)^-1=F^-1*FT^-1

               KRATOS_ERROR << "calling HenckyHyperelastic based method with PK2 stress. not implemented" << std::endl;

            }
            else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b=F*FT)

               //set working strain measure
               rValues.SetStrainMeasure(ConstitutiveModelData::StrainMeasureType::CauchyGreen_Left);

               MatrixType EigenVectors;
               MatrixType EigenValues;
               EigenVectors.clear();
               EigenValues.clear();

               MathUtils<double>::EigenSystem<3> ( rValues.StrainMatrix, EigenVectors, EigenValues);

               rVariables.Strain.Matrix.clear();
               for (unsigned int i = 0; i < 3; i++)
                  rVariables.Strain.Matrix(i,i) =  std::log(EigenValues(i,i)) / 2.0;

               rVariables.Strain.Matrix = prod( rVariables.Strain.Matrix, EigenVectors);
               rVariables.Strain.Matrix = prod( trans(EigenVectors), rVariables.Strain.Matrix);


               rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);

            }
            else{
               KRATOS_ERROR << "calling initialize HyperElasticModel .. StressMeasure required is inconsistent"  << std::endl;
            }


            KRATOS_CATCH("")
         };

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticModel )
         }

         void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticModel )
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class BorjaModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HENCKY_HYPER_ELASTIC_MODEL_H_INCLUDED  defined
