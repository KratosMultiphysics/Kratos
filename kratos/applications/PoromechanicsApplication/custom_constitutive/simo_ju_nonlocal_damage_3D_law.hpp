//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_SIMO_JU_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_SIMO_JU_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_flow_rules/nonlocal_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SimoJuNonlocalDamage3DLaw : public NonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SimoJuNonlocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SimoJuNonlocalDamage3DLaw();
    
    /// Second Constructor
    SimoJuNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    SimoJuNonlocalDamage3DLaw (const SimoJuNonlocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~SimoJuNonlocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);
    
    ConstitutiveLaw::Pointer Clone() const;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry );

private:
    
    /// Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class SimoJuNonlocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_SIMO_JU_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined 
