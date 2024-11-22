//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/custom_flow_rules/local_damage_flow_rule.hpp"

namespace Kratos
{

//Default Constructor
LocalDamageFlowRule::LocalDamageFlowRule() : IsotropicDamageFlowRule() {}

//----------------------------------------------------------------------------------------

//Second Constructor
LocalDamageFlowRule::LocalDamageFlowRule(YieldCriterionPointer pYieldCriterion)
	:IsotropicDamageFlowRule(pYieldCriterion) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LocalDamageFlowRule::LocalDamageFlowRule(LocalDamageFlowRule const& rOther)
	:IsotropicDamageFlowRule(rOther) {}

//----------------------------------------------------------------------------------------

//Assignment Operator
LocalDamageFlowRule& LocalDamageFlowRule::operator=(LocalDamageFlowRule const& rOther)
{
   IsotropicDamageFlowRule::operator=(rOther);
   return *this;
}

//----------------------------------------------------------------------------------------

//Destructor
LocalDamageFlowRule::~LocalDamageFlowRule() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FlowRule::Pointer LocalDamageFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new LocalDamageFlowRule(*this));
  return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LocalDamageFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{   
    // SetValue STATE_VARIABLE
    mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.TrialStateFunction;
    mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.TrialStateFunction;
}

} // Namespace Kratos
