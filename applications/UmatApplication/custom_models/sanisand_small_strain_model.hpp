//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:             AFranza $
//   Last modified by:    $Co-Author:        LMonforte $
//   Date:                $Date:          October 2019 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_SANISAND_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_SANISAND_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED

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
   class KRATOS_API(UMAT_APPLICATION) SanisandSmallStrainUmatModel : public SmallStrainUmatModel
   {


      public:

         ///@name Type Definitions
         ///@{

         /// Pointer definition of SanisandSmallStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( SanisandSmallStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         SanisandSmallStrainUmatModel();

         /// Copy constructor.
         SanisandSmallStrainUmatModel(SanisandSmallStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         SanisandSmallStrainUmatModel& operator=(SanisandSmallStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~SanisandSmallStrainUmatModel();


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
            buffer << "SanisandSmallStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "SanisandSmallStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "SanisandSmallStrainUmatModel Data";
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
            return 36;
         };

         
         /*
            Create the vector with constitutive parameters value
          */
         virtual void CreateConstitutiveParametersVector(double* & pVector, int & rNumberParameters, const Properties & rMaterialProperties) override {
            rNumberParameters = 19;
            pVector = new double[rNumberParameters];

            pVector[0] = 101325; // p_a atmospheric pressure
            pVector[1] = 0.8191; // e0
            pVector[2] = 0.00178; // lambda 
            pVector[3] = 2.4352; // epsi 
            pVector[4] = 1.287; // M_c or phi_c 
            pVector[5] = 0; // M_e or phi_e
            pVector[6] = 0.01; // m
            pVector[7] = 400000000; // Go
            pVector[8] = 0.05; // nu
            pVector[9] = 4.05; // h_0
            pVector[10] = 1.100; // c_h
            pVector[11] = 2.800; // n^b
            pVector[12] = 0.550; // A_0
            pVector[13] = 2.564; // n_d
            pVector[14] = 0; // z_max
            pVector[15] = 0; // c_z
            pVector[16] = 0; // K_w
            pVector[17] = 0.01; // p_tmult stabilising param
            pVector[18] = 0.6531; // e

         };

         /*
            Number of the constitutive equation in the fortran wrapper
          */
         virtual int GetConstitutiveEquationNumber() override
         {
            return 3;
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

   }; // Class SanisandSmallStrainUmatModel

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


