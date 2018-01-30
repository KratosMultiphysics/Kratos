//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:             JMCarbonell $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_EXPONENTIAL_DAMAGE_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_EXPONENTIAL_DAMAGE_HARDENING_RULE_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ExponentialDamageHardeningRule
    : public HardeningRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExponentialDamageHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( ExponentialDamageHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExponentialDamageHardeningRule();

    /// Copy constructor.
    ExponentialDamageHardeningRule(ExponentialDamageHardeningRule const& rOther);

    /// Assignment operator.
    ExponentialDamageHardeningRule& operator=(ExponentialDamageHardeningRule const& rOther);

    /// Clone.
    virtual HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~ExponentialDamageHardeningRule();

    ///@}
    ///@name Operators
    ///@{

 
    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Hardening functions
     */

    virtual double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening) override;
      
    /**
     * Calculate Hardening function derivatives
     */

    virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening) override;
    
 
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
      buffer << "ExponentialDamageHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ExponentialDamageHardeningRule";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "ExponentialDamageHardeningRule Data";
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningRule )
    }
    
    virtual void load(Serializer& rSerializer) override
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

  }; // Class ExponentialDamageHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXPONENTIAL_DAMAGE_HARDENING_RULE_H_INCLUDED  defined 


