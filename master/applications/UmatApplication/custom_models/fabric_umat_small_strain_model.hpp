//
//   Project Name:        KratosUmatApplication        $
//   Created by:          $Author:           LMonforte $
//   Last modified by:    $Co-Author:                  $
//   Date:                $Date:            April 2018 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_FABRIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED )
#define  KRATOS_FABRIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED

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
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) FabricSmallStrainUmatModel : public SmallStrainUmatModel
   {


      public:

         ///@name Type Definitions
         ///@{

         /// Pointer definition of FabricSmallStrainUmatModel
         KRATOS_CLASS_POINTER_DEFINITION( FabricSmallStrainUmatModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.    
         FabricSmallStrainUmatModel();

         /// Copy constructor.
         FabricSmallStrainUmatModel(FabricSmallStrainUmatModel const& rOther);

         /// Clone.
         virtual ConstitutiveModel::Pointer Clone() const override;

         /// Assignment operator.
         FabricSmallStrainUmatModel& operator=(FabricSmallStrainUmatModel const& rOther);

         /// Destructor.
         virtual ~FabricSmallStrainUmatModel();


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
         ///@name Access
         ///@{

         double& GetValue(const Variable<double> & , double& rValue) override;

         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         virtual std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "FabricSmallStrainUmatModel";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "FabricSmallStrainUmatModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "FabricSmallStrainUmatModel Data";
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
            return 16;
         };

         virtual void InitializeStateVariables( Vector& rStateVariables, const Properties & rMaterialProperties) override;
        
         /*
            Create the vector with constitutive parameters value
          */
         virtual void CreateConstitutiveParametersVector(double* & pVector, int & rNumberParameters, const Properties & rMaterialProperties) override;
         /*
            Number of the constitutive equation in the fortran wrapper
          */
         virtual int GetConstitutiveEquationNumber() override
         {
            return 2;
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

   }; // Class FabricSmallStrainUmatModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FABRIC_SMALL_STRAIN_UMAT_MODEL_H_INCLUDED  defined 


