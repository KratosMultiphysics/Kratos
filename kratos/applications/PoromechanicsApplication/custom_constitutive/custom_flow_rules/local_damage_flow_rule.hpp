//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_LOCAL_DAMAGE_FLOW_RULE_H_INCLUDED )
#define  KRATOS_LOCAL_DAMAGE_FLOW_RULE_H_INCLUDED

// Application includes
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LocalDamageFlowRule : public IsotropicDamageFlowRule
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LocalDamageFlowRule );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    LocalDamageFlowRule();

    /// Second Constructor
    LocalDamageFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor
    LocalDamageFlowRule(LocalDamageFlowRule const& rOther);

    /// Assignment operator
    LocalDamageFlowRule& operator=(LocalDamageFlowRule const& rOther);
	
    /// Destructor
    virtual ~LocalDamageFlowRule();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    FlowRule::Pointer Clone() const;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    void CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
    }

}; // Class LocalDamageFlowRule
}  // namespace Kratos.
#endif // KRATOS_LOCAL_DAMAGE_FLOW_RULE_H_INCLUDED  defined 