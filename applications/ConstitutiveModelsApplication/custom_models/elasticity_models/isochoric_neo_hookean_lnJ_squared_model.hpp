//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_mooney_rivlin_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricNeoHookeanLnJSquaredModel : public IsochoricNeoHookeanModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricNeoHookeanLnJSquaredModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricNeoHookeanLnJSquaredModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricNeoHookeanLnJSquaredModel() : IsochoricNeoHookeanModel() {}

    /// Copy constructor.
    IsochoricNeoHookeanLnJSquaredModel(IsochoricNeoHookeanLnJSquaredModel const& rOther) : IsochoricNeoHookeanModel(rOther) {}

    /// Assignment operator.
    IsochoricNeoHookeanLnJSquaredModel& operator=(IsochoricNeoHookeanLnJSquaredModel const& rOther)
    {
      IsochoricNeoHookeanModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IsochoricNeoHookeanLnJSquaredModel>(*this);
    }

    /// Destructor.
    ~IsochoricNeoHookeanLnJSquaredModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    // Simplyfied methods must be implemented for performance purposes

    /**
     * Calculate Stresses
     */



    /**
     * Calculate Constitutive Components
     */



    /**
     * Check
     */

    int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      HyperElasticModel::Check(rProperties,rCurrentProcessInfo);

      if( C10.Key() == 0 || rProperties[C10] <= 0.00 )
	KRATOS_ERROR << "C10 has an invalid key or value" << std::endl;

      if( BULK_MODULUS.Key() == 0 || rProperties[BULK_MODULUS] <= 0.00 )
	KRATOS_ERROR << "BULK_MODULUS has an invalid key or value" << std::endl;

      return 0;

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
        buffer << "IsochoricNeoHookeanLnJSquaredModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricNeoHookeanLnJSquaredModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricNeoHookeanLnJSquaredModel Data";
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

    //specialized methods:

    void CalculateVolumetricFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      rFactor = std::log(rVariables.Strain.Invariants.J);

      KRATOS_CATCH(" ")
    }


    void CalculateConstitutiveMatrixFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      rFactor = 1.0;

      KRATOS_CATCH(" ")
    }

    //************// W

    void CalculateAndAddIsochoricStrainEnergy(HyperElasticDataType& rVariables, double& rIsochoricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      rIsochoricDensityFunction += rMaterial.GetModelParameters()[0] * ( rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.I1 - 3.0);

      KRATOS_CATCH(" ")
    }


    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //energy function "U(J) = (K/2)*(lnJ)²"
      rVolumetricDensityFunction += rMaterial.GetBulkModulus() * 0.5 * pow(std::log(rVariables.Strain.Invariants.J),2);

      KRATOS_CATCH(" ")
    }

    //************// dW

    double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI1
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      rDerivative = rMaterial.GetModelParameters()[0];

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI3
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //dU/dJ
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "U(J) = (K/2)*ln(J)²"
      //dU(J)/dJ = (K)*(lnJ/J)
      rDerivative = rMaterial.GetBulkModulus() * std::log( rVariables.Strain.Invariants.J );

      rDerivative /= rVariables.Strain.Invariants.J;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI1dI1
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI2dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI3dI3
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddU/dJdJ
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "dU(J)/dJ = (K)*(lnJ/J)"
      //ddU(J)/dJdJ = (K)*(1-lnJ)/J²
      rDerivative = rMaterial.GetBulkModulus() * (1.0 -std::log(rVariables.Strain.Invariants.J)) / (rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J);

      return rDerivative;

      KRATOS_CATCH(" ")
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


    void save(Serializer& rSerializer) const  override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IsochoricNeoHookeanLnJSquaredModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED  defined
