//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/custom_flow_rules/restore_damage_flow_rule.hpp"

namespace Kratos
{

//Default Constructor
RestoreDamageFlowRule::RestoreDamageFlowRule() : IsotropicDamageFlowRule() {}

//----------------------------------------------------------------------------------------

//Second Constructor
RestoreDamageFlowRule::RestoreDamageFlowRule(YieldCriterionPointer pYieldCriterion)
	:IsotropicDamageFlowRule(pYieldCriterion) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
RestoreDamageFlowRule::RestoreDamageFlowRule(RestoreDamageFlowRule const& rOther)
	:IsotropicDamageFlowRule(rOther) {}

//----------------------------------------------------------------------------------------

//Assignment Operator
RestoreDamageFlowRule& RestoreDamageFlowRule::operator=(RestoreDamageFlowRule const& rOther)
{
   IsotropicDamageFlowRule::operator=(rOther);
   return *this;
}

//----------------------------------------------------------------------------------------

//Destructor
RestoreDamageFlowRule::~RestoreDamageFlowRule() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FlowRule::Pointer RestoreDamageFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new RestoreDamageFlowRule(*this));
  return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreDamageFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{   
    // SetValue STATE_VARIABLE
    mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.TrialStateFunction;
    mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.TrialStateFunction;
}

} // Namespace Kratos
