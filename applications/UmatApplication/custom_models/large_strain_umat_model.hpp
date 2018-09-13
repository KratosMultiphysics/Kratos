//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_LARGE_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_LARGE_STRAIN_UMAT_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/small_strain_umat_model.hpp"

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
   class KRATOS_API(UMAT_APPLICATION) LargeStrainUmatModel : public SmallStrainUmatModel
   {
      public:

         ///@name Type Definitions
         ///@{
         typedef SmallStrainUmatModel::UmatModelData              UmatDataType;


         /// Pointer definition of LargeStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( LargeStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         LargeStrainUmatModel();

         /// Copy constructor.
         LargeStrainUmatModel(LargeStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         LargeStrainUmatModel& operator=(LargeStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~LargeStrainUmatModel();


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
            buffer << "LargeStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "LargeStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "LargeStrainUmatModel Data";
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


         /*
            Create strain_n and incremental strain
            ( for large strains should be overrided ) 
          */

         virtual void CreateStrainsVectors( UmatDataType & rVariables, double* & rpStrain, double* & rpIncrementalStrain) override;

         /* 
            Create stress_n
            ( for large strains should be overrided )
          */
         virtual void CreateStressAtInitialState( UmatDataType & rVariables, double* & rpStressVector) override;

         /*
            Add more constitutive terms to the matrix
            */
         Matrix& CalculateExtraMatrix(const MatrixType& rStressMatrix, Matrix& rExtraMatrix);

         
         /*
            Set the constitutiveMatrixToThe appropiate size
            */

         virtual void SetConstitutiveMatrix( Matrix & rC, const Matrix & rCBig, const MatrixType& rStressMatrix) override;

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallStrainUmatModel )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallStrainUmatModel )
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class LargeStrainUmatModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LARGE_STRAIN_UMAT_MODEL_H_INCLUDED  defined 


