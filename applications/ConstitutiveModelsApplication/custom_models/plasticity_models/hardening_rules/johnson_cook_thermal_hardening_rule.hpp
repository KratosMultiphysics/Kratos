//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_JOHNSON_COOK_THERMAL_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_THERMAL_HARDENING_RULE_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) JohnsonCookThermalHardeningRule 
    : public HardeningRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of JohnsonCookThermalHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( JohnsonCookThermalHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    JohnsonCookThermalHardeningRule();

    /// Copy constructor.
    JohnsonCookThermalHardeningRule(JohnsonCookThermalHardeningRule const& rOther);

    /// Assignment operator.
    JohnsonCookThermalHardeningRule& operator=(JohnsonCookThermalHardeningRule const& rOther);

    /// Clone.
    virtual HardeningRule::Pointer Clone() const override;
    
    /// Destructor.
    ~JohnsonCookThermalHardeningRule();

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

    virtual double& CalculateDeltaThermalHardening(const PlasticDataType& rVariables, double& rDeltaThermalHardening) override;
    
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
      buffer << "JohnsonCookThermalHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "JohnsonCookThermalHardeningRule";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "JohnsonCookThermalHardeningRule Data";
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

  }; // Class JohnsonCookThermalHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_JOHNSON_COOK_THERMAL_HARDENING_RULE_H_INCLUDED  defined 


