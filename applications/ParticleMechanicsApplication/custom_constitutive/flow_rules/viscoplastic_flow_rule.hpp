//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:            Duan Wenjie $
//   Date:                $Date:                March 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED )
#define  KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED

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

class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) ViscoplasticFlowRule
    :public NonLinearAssociativePlasticFlowRule
{
public:
///@name Type Definitions
///@{

/// Pointer definition of NonLinearAssociativePlasticFlowRule
  KRATOS_CLASS_POINTER_DEFINITION( ViscoplasticFlowRule );

///@}
///@name Life Cycle
///@{

/// Default constructor.
ViscoplasticFlowRule();

/// Initialization constructor.
ViscoplasticFlowRule(YieldCriterionPointer pYieldCriterion);

/// Copy constructor.
ViscoplasticFlowRule(ViscoplasticFlowRule const& rOther);

/// Assignment operator.
ViscoplasticFlowRule& operator=(ViscoplasticFlowRule const& rOther);

/// Destructor.
virtual ~ViscoplasticFlowRule();

///@}
///@name Operators
///@{

/**
 * Clone function (has to be implemented by any derived class)
 * @return a pointer to a new instance of this flow rule
 */
virtual FlowRule::Pointer Clone() const;

///@}
///@name Operations
///@{
/// virtual bool CalculateReturnMapping(  RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix );

///virtual void CalculateScalingFactors( const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors );

virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables );


///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

// /// Turn back information as a string.
// virtual std::string Info() const;

// /// Print information about this object.
// virtual void PrintInfo(std::ostream& rOStream) const;

// /// Print object's data.
// virtual void PrintData(std::ostream& rOStream) const;


///@}
///@name Friends
///@{


///@}
///
protected:
///double& CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm );


///virtual void SetCriterionParameters( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters );


virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters);


void UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );


///void CalculateThermalDissipation( YieldCriterion::Parameters& rCriterionParameters, ThermalVariables& rThermalVariables );
///@}
///@name Protected  Access
///@{


///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{
private:
friend class Serializer;

// A private default constructor necessary for serialization

virtual void save(Serializer& rSerializer) const;

virtual void load(Serializer& rSerializer);



};// Class ViscoplasticFlowRule

}  // namespace Kratos.


#endif // KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED  defined
