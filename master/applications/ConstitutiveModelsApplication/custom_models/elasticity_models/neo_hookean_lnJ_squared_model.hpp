//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED)
#define  KRATOS_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/neo_hookean_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NeoHookeanLnJSquaredModel : public NeoHookeanModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NeoHookeanLnJSquaredModel
    KRATOS_CLASS_POINTER_DEFINITION(NeoHookeanLnJSquaredModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NeoHookeanLnJSquaredModel() : NeoHookeanModel() {}

    /// Copy constructor.
    NeoHookeanLnJSquaredModel(NeoHookeanLnJSquaredModel const& rOther) : NeoHookeanModel(rOther) {}

    /// Assignment operator.
    NeoHookeanLnJSquaredModel& operator=(NeoHookeanLnJSquaredModel const& rOther)
    {
      NeoHookeanModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<NeoHookeanLnJSquaredModel>(*this);
    }

    /// Destructor.
    ~NeoHookeanLnJSquaredModel() override {}


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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NeoHookeanLnJSquaredModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NeoHookeanLnJSquaredModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NeoHookeanLnJSquaredModel Data";
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


    // specialized methods:

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

    //option: g(J) = (lambda/2)*ln(J)² - (mu)*lnJ (neo_hookean_lnJ_model)

    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //g(J) = (lambda/2)*ln(J)² - (mu)*lnJ
      rVolumetricDensityFunction = rMaterial.GetLameLambda() * 0.5 * std::log( rVariables.Strain.Invariants.J );
      rVolumetricDensityFunction -= rMaterial.GetLameMu() * std::log( rVariables.Strain.Invariants.J );

      KRATOS_CATCH(" ")
    }


    //************// dW


    double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //dW/dI3
    {
      KRATOS_TRY

      const MaterialDataType&  rMaterial = rVariables.GetMaterialParameters();

      //derivative of "g(J) = (lambda/2)*(lnJ)² - (mu)*lnJ"
      //dg(J)/dI3 = (lambda/2)*(lnJ/J²) - (mu)*(1/J²)/2
      rDerivative  = 0.5 * rMaterial.GetLameLambda() * std::log( rVariables.Strain.Invariants.J );
      rDerivative -= 0.5 * rMaterial.GetLameMu();
      rDerivative /= rVariables.Strain.Invariants.I3;

      return rDerivative;

      KRATOS_CATCH(" ")
    }


    double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddW/dI3dI3
    {
      KRATOS_TRY

      const MaterialDataType&   rMaterial = rVariables.GetMaterialParameters();

      //ddg(J)/dI3dI3 = (lambda/4)*(1-2*lnJ)/J⁴ + (mu/2)*(1/J⁴)
      rDerivative  = 0.25 * rMaterial.GetLameLambda() * (1.0 - 2.0 * std::log( rVariables.Strain.Invariants.J ) );
      rDerivative += 0.5 * rMaterial.GetLameMu();
      rDerivative /= (rVariables.Strain.Invariants.I3 * rVariables.Strain.Invariants.I3);

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


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NeoHookeanModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NeoHookeanModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NeoHookeanLnJSquaredModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED  defined
