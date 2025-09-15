//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SIMO_EXPONENTIAL_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_SIMO_EXPONENTIAL_HARDENING_RULE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/hardening_rule.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SimoExponentialHardeningRule
    : public HardeningRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SimoExponentialHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( SimoExponentialHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoExponentialHardeningRule();

    /// Copy constructor.
    SimoExponentialHardeningRule(SimoExponentialHardeningRule const& rOther);

    /// Assignment operator.
    SimoExponentialHardeningRule& operator=(SimoExponentialHardeningRule const& rOther);

    /// Clone.
    HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~SimoExponentialHardeningRule() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Calculate Hardening functions
     */

    double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening) override;

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
      buffer << "SimoExponentialHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SimoExponentialHardeningRule";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoExponentialHardeningRule Data";
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


    /**
     * Pure isotropic hardening Theta=1;  pure kinematic hardening Theta= 0; combined isotropic-kinematic 0<Theta<1
     */
    constexpr static const double mTheta = 1.0;


    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculate Hardening functions
     */
    virtual double& CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening);

    virtual double& CalculateAndAddKinematicHardening(const PlasticDataType& rVariables, double& rKinematicHardening);

    /**
     * Calculate Hardening function derivatives
     */
    virtual double& CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening);

    virtual double& CalculateAndAddDeltaKinematicHardening(const PlasticDataType& rVariables, double& rDeltaKinematicHardening);



    virtual double& CalculateThermalReferenceEffect(const PlasticDataType& rVariables, double& rThermalFactor);

    virtual double& CalculateThermalCurrentEffect(const PlasticDataType& rVariables, double& rThermalFactor);


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningRule )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningRule )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class SimoExponentialHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_EXPONENTIAL_HARDENING_RULE_H_INCLUDED  defined


