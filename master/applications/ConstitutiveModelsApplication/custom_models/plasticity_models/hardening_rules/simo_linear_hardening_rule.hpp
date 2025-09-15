//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SIMO_LINEAR_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_SIMO_LINEAR_HARDENING_RULE_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SimoLinearHardeningRule
    : public SimoExponentialHardeningRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SimoLinearHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( SimoLinearHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoLinearHardeningRule();


    /// Copy constructor.
    SimoLinearHardeningRule(SimoLinearHardeningRule const& rOther);

    /// Assignment operator.
    SimoLinearHardeningRule& operator=(SimoLinearHardeningRule const& rOther);

    /// Clone.
    HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~SimoLinearHardeningRule() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Calculate Hardening function derivatives
     */

    double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening) override;


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
      buffer << "SimoLinearHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SimoLinearHardeningRule";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoLinearHardeningRule Data";
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

    /**
     * Calculate Hardening functions
     */
    double& CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening) override;

    /**
     * Calculate Hardening function derivatives
     */
    double& CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening) override;

    ///@}
    ///@name Protected Operations
    ///@{


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

  }; // Class SimoLinearHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_LINEAR_HARDENING_RULE_H_INCLUDED  defined


