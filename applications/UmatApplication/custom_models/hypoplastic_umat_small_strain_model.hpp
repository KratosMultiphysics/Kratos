//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:          October 2017 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_HYPOPLASTIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_HYPOPLASTIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED

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
   class KRATOS_API(UMAT_APPLICATION) HypoplasticSmallStrainUmatModel : public SmallStrainUmatModel
   {


      public:

         ///@name Type Definitions
         ///@{

         /// Pointer definition of HypoplasticSmallStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( HypoplasticSmallStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         HypoplasticSmallStrainUmatModel();

         /// Copy constructor.
         HypoplasticSmallStrainUmatModel(HypoplasticSmallStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         HypoplasticSmallStrainUmatModel& operator=(HypoplasticSmallStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~HypoplasticSmallStrainUmatModel();


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
            buffer << "HypoplasticSmallStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "HypoplasticSmallStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "HypoplasticSmallStrainUmatModel Data";
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
            rNumberParameters = 16;
            pVector = new double[rNumberParameters];
            pVector[0] = 30.0; // phi_deg
            pVector[1] = 1.0; // p_t
            pVector[2] = 5800.0; // hs
            pVector[3] = 0.28; // e_n
            pVector[4] = 0.53; // ed0
            pVector[5] = 0.84; // ec0
            pVector[6] = 1.0; // ei0
            pVector[7] = 0.130; // alpha
            pVector[8] = 1.00; // beta
            pVector[9] = 2.0; // m_R
            pVector[10] = 5.0; // m_T
            pVector[11] = 1.0e-4; // r_uc
            pVector[12] = 0.05; // beta_r
            pVector[13] = 1.0; // chi
            pVector[14] = 0.0E5; // bulk_w  // should be set equal to zero
            pVector[15] = 10.0+0.5; // ?? line 155 (some sort of void ratio that depends on the number does something)
         };

         /*
            Number of the constitutive equation in the fortran wrapper
          */
         virtual int GetConstitutiveEquationNumber() override
         {
            return 1;
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

   }; // Class HypoplasticSmallStrainUmatModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HYPOPLASTIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED  defined 


