//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_VON_MISES_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_VON_MISES_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) VonMisesSmallStrainUmatModel : public SmallStrainUmatModel
   {


      public:

         ///@name Type Definitions
         ///@{

         /// Pointer definition of VonMisesSmallStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( VonMisesSmallStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         VonMisesSmallStrainUmatModel();

         /// Copy constructor.
         VonMisesSmallStrainUmatModel(VonMisesSmallStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         VonMisesSmallStrainUmatModel& operator=(VonMisesSmallStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~VonMisesSmallStrainUmatModel();


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
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
            buffer << "VonMisesSmallStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "VonMisesSmallStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "VonMisesSmallStrainUmatModel Data";
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
            Get the dimension of StateVariables 
          */

         virtual unsigned int GetNumberOfStateVariables() override {
            return 13;
         };

         
         /*
            Create the vector with constitutive parameters value
          */
         virtual void CreateConstitutiveParametersVector(double* & pVector, int & rNumberParameters, const Properties & rMaterialProperties) override {
            rNumberParameters = 3;
            pVector = new double[rNumberParameters];
            pVector[0] = rMaterialProperties[YOUNG_MODULUS]; // young
            pVector[1] = rMaterialProperties[POISSON_RATIO]; // poisson
            pVector[2] = rMaterialProperties[YIELD_STRESS]; // yield
         };

         /*
            Number of the constitutive equation in the fortran wrapper
          */
         virtual int GetConstitutiveEquationNumber() override
         {
            return 0;
         }
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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveModel )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveModel )
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class VonMisesSmallStrainUmatModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VON_MISES_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED  defined 


