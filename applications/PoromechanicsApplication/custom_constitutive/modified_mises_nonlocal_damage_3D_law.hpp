//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.hpp"
#include "custom_constitutive/custom_flow_rules/nonlocal_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ModifiedMisesNonlocalDamage3DLaw : public NonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMisesNonlocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ModifiedMisesNonlocalDamage3DLaw();
    
    /// Second Constructor
    ModifiedMisesNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    ModifiedMisesNonlocalDamage3DLaw (const ModifiedMisesNonlocalDamage3DLaw& rOther);

    /// Destructor
    ~ModifiedMisesNonlocalDamage3DLaw() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;
    
    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class ModifiedMisesNonlocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_MODIFIED_MISES_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined 
