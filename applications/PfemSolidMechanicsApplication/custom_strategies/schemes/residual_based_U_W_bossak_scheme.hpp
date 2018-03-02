//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              LMonforte $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            February 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_U_W_BOSSAK_SCHEME )
#define      KRATOS_RESIDUAL_BASED_U_W_BOSSAK_SCHEME

// System includes

// External includes

// Project includes
#include "custom_strategies/schemes/residual_based_displacement_bossak_scheme.hpp"
#include "custom_strategies/time_integration_methods/bossak_method.hpp"

namespace Kratos
{
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

   /** @brief Bossak integration scheme (for dynamic problems)
    */
   template<class TSparseSpace,  class TDenseSpace >
      class ResidualBasedUWBossakScheme: public ResidualBasedDisplacementBossakScheme<TSparseSpace,TDenseSpace>
   {   
      public:

         ///@name Type Definitions
         ///@{
         KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedUWBossakScheme );

         typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

         typedef typename Element::DofsArrayType                  DofsArrayType;

         typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

         typedef typename BaseType::TSystemVectorType         TSystemVectorType;

         typedef typename BaseType::Pointer                     BaseTypePointer;

         typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

         typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

         typedef ResidualBasedDisplacementBossakScheme<TSparseSpace,TDenseSpace>  DerivedType;

         typedef ModelPart::NodesContainerType                   NodesArrayType;

         typedef typename DerivedType::IntegrationTypePointer           IntegrationTypePointer;

         typedef typename DerivedType::NodeType                                       NodeType;

         typedef TimeIntegrationMethod< Variable<double>, double>    ScalarIntegrationType;

         typedef typename ScalarIntegrationType::Pointer     ScalarIntegrationTypePointer;


         ///@}
         ///@name Life Cycle
         ///@{

         /// Default Constructor.
         ResidualBasedUWBossakScheme()
            :DerivedType()
         {
         }

         /// Copy Constructor.
         ResidualBasedUWBossakScheme(ResidualBasedUWBossakScheme& rOther)
            :DerivedType(rOther)
         {
         }

         /// Clone.
         BaseTypePointer Clone() override
         {
            return BaseTypePointer( new ResidualBasedUWBossakScheme(*this) );
         }

         /// Destructor.
         ~ResidualBasedUWBossakScheme() override {}

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
            buffer << "Displacement BossakScheme";
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "Displacement BossakScheme";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "Displacement BossakScheme Data";     
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

         IntegrationTypePointer mpWaterDisplacementIntegrationMethod;

         ScalarIntegrationTypePointer mpWaterPressureIntegrationMethod;

         ///@}
         ///@name Protected Operators
         ///@{

         ///@}
         ///@name Protected Operations
         ///@{

         virtual void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
         {
            KRATOS_TRY

            this->mpIntegrationMethod = IntegrationTypePointer( new BossakMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

            // Set scheme variables
            this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

            // Set scheme parameters
            this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);

            // Modify ProcessInfo scheme parameters
            this->mpIntegrationMethod->SetProcessInfoParameters(rCurrentProcessInfo);

            this->mpWaterDisplacementIntegrationMethod = IntegrationTypePointer( new BossakMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

            // Set scheme variables
            this->mpWaterDisplacementIntegrationMethod->SetVariables(WATER_DISPLACEMENT,WATER_VELOCITY,WATER_ACCELERATION);

            // Set scheme parameters
            this->mpWaterDisplacementIntegrationMethod->SetParameters(rCurrentProcessInfo);

            this->mpWaterPressureIntegrationMethod = ScalarIntegrationTypePointer( new BossakMethod< Variable<double>, double>);

            // Set scheme variables
            this->mpWaterPressureIntegrationMethod->SetVariables(WATER_PRESSURE,WATER_PRESSURE_VELOCITY,WATER_PRESSURE_ACCELERATION);

            // Set scheme parameters
            this->mpWaterPressureIntegrationMethod->SetParameters(rCurrentProcessInfo);
            KRATOS_CATCH("")
         }



         virtual void IntegrationMethodUpdate(NodeType& rNode) override
         {
            KRATOS_TRY

            this->mpIntegrationMethod->Update(rNode);

            mpWaterDisplacementIntegrationMethod->Update(rNode);

            if ( rNode.HasDofFor(WATER_PRESSURE) )
            {
               mpWaterPressureIntegrationMethod->Update( rNode);
            }


            KRATOS_CATCH("")
         }


         virtual void IntegrationMethodPredict(NodeType& rNode) override
         {
            KRATOS_TRY

            this->mpIntegrationMethod->Predict(rNode);

            mpWaterDisplacementIntegrationMethod->Predict(rNode);

            if ( rNode.HasDofFor(WATER_PRESSURE) )
            {
               mpWaterPressureIntegrationMethod->Predict( rNode);
            }

            KRATOS_CATCH("")
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

         ///@}
         ///@name Private Inquiry
         ///@{

         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}
   }; // Class ResidualBasedUWBossakScheme
   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_U_W_BOSSAK_SCHEME defined
