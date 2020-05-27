//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

#if !defined(KRATOS_JOHNSON_COOK_PLASTIC_FLOW_RULE_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/flow_rules/particle_flow_rule.hpp"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
  class JohnsonCookPlasticFlowRule
	  :public ParticleFlowRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearRateDependentPlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookPlasticFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
      JohnsonCookPlasticFlowRule();

    /// Initialization constructor.
      JohnsonCookPlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
      JohnsonCookPlasticFlowRule(JohnsonCookPlasticFlowRule const& rOther);

    /// Assignment operator.
    JohnsonCookPlasticFlowRule& operator=(JohnsonCookPlasticFlowRule const& rOther);

    /// Destructor.
    ~JohnsonCookPlasticFlowRule() override;


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    ParticleFlowRule::Pointer Clone() const override;


    ///@}
    ///@name Operations
    ///@{

    bool CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix) override;


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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
    double& CalculateStressNorm(Matrix& rStressMatrix, double& rStressNorm) override;

    double CalculateThermalDerivative(const ParticleHardeningLaw::Parameters& rValues);

    double CalculatePlasticStrainRateDerivative(const ParticleHardeningLaw::Parameters& rValues);

    double CalculatePlasticStrainDerivative(const ParticleHardeningLaw::Parameters& rValues);



   

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

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class JohnsonCookPlasticFlowRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{



  ///@}

  ///@} addtogroup block



  ///@}
  ///@ Template Operations
  ///@{


  ///@}


}  // namespace Kratos.

#endif // KRATOS_JOHNSON_COOK_PLASTIC_FLOW_RULE_H_INCLUDED  defined


