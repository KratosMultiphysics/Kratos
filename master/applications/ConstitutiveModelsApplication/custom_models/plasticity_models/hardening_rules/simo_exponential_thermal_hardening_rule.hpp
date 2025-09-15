//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SIMO_EXPONENTIAL_THERMAL_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_SIMO_EXPONENTIAL_THERMAL_HARDENING_RULE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/simo_exponential_hardening_rule.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SimoExponentialThermalHardeningRule
    : public SimoExponentialHardeningRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SimoExponentialThermalHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( SimoExponentialThermalHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoExponentialThermalHardeningRule();

    /// Copy constructor.
    SimoExponentialThermalHardeningRule(SimoExponentialThermalHardeningRule const& rOther);

    /// Assignment operator.
    SimoExponentialThermalHardeningRule& operator=(SimoExponentialThermalHardeningRule const& rOther);

    /// Clone.
    HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~SimoExponentialThermalHardeningRule() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Calculate Hardening function derivatives
     */

    double& CalculateDeltaThermalHardening(const PlasticDataType& rVariables, double& rDeltaThermalHardening) override;


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
      buffer << "SimoExponentialThermalHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SimoExponentialThermalHardeningRule";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoExponentialThermalHardeningRule Data";
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

    double& CalculateThermalReferenceEffect(const PlasticDataType& rVariables, double& rThermalFactor) override;

    double& CalculateThermalCurrentEffect(const PlasticDataType& rVariables, double& rThermalFactor) override;

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SimoExponentialHardeningRule )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SimoExponentialHardeningRule )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class SimoExponentialThermalHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_EXPONENTIAL_THERMAL_HARDENING_RULE_H_INCLUDED  defined


