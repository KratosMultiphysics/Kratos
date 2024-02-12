//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NONLOCAL_DAMAGE_FLOW_RULE_H_INCLUDED )
#define  KRATOS_NONLOCAL_DAMAGE_FLOW_RULE_H_INCLUDED

// Application includes
#include "custom_constitutive/custom_flow_rules/local_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) NonlocalDamageFlowRule : public LocalDamageFlowRule
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NonlocalDamageFlowRule );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    NonlocalDamageFlowRule();

    /// Second Constructor
    NonlocalDamageFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor
    NonlocalDamageFlowRule(NonlocalDamageFlowRule const& rOther);

    /// Assignment operator
    NonlocalDamageFlowRule& operator=(NonlocalDamageFlowRule const& rOther);
	
    /// Destructor
    ~NonlocalDamageFlowRule() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    FlowRule::Pointer Clone() const override;
    
    void InitializeMaterial (YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, 
                                Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen) override;
    
    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
protected:

    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateLocalInternalVariables(RadialReturnVariables& rReturnMappingVariables);
    
    bool CalculateInternalVariables(RadialReturnVariables& rReturnMappingVariables) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
    }

}; // Class NonlocalDamageFlowRule
}  // namespace Kratos.
#endif // KRATOS_NONLOCAL_DAMAGE_FLOW_RULE_H_INCLUDED  defined 