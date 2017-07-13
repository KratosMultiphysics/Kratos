//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:            Duan Wenjie $
//   Date:                $Date:                March 2016 $
//   Revision:            $Revision:                  0.0 $
//
//
#if !defined(KRATOS_BINGHAM_VISCOPLASTIC_FLOW_RULE_H_INCLUDED )
#define  KRATOS_BINGHAM_VISCOPLASTIC_FLOW_RULE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"

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

class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) BinghamViscoplasticFlowRule
    :public NonLinearAssociativePlasticFlowRule
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( BinghamViscoplasticFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BinghamViscoplasticFlowRule();

    /// Initialization constructor.
    BinghamViscoplasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BinghamViscoplasticFlowRule(BinghamViscoplasticFlowRule const& rOther);

    /// Assignment operator.
    BinghamViscoplasticFlowRule& operator=(BinghamViscoplasticFlowRule const& rOther);

    /// Destructor.
    virtual ~BinghamViscoplasticFlowRule();

	virtual void CalculateScalingFactors( const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors );
    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    virtual FlowRule::Pointer Clone() const;

protected:
    virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters);

private:
    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

};// class BinghamViscoplasticFlowRule

}  // namespace Kratos.



#endif // KRATOS_BINGHAM_VISCOPLASTIC_FLOW_RULE_H_INCLUDED  defined
