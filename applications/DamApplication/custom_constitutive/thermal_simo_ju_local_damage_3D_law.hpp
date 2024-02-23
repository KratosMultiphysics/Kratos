//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_SIMO_JU_LOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_SIMO_JU_LOCAL_DAMAGE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/thermal_local_damage_3D_law.hpp"
#include "custom_constitutive/continuum_laws/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/continuum_laws/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/continuum_laws/custom_flow_rules/local_damage_flow_rule.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalSimoJuLocalDamage3DLaw : public ThermalLocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalSimoJuLocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalSimoJuLocalDamage3DLaw();

    /// Second Constructor
    ThermalSimoJuLocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalSimoJuLocalDamage3DLaw (const ThermalSimoJuLocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~ThermalSimoJuLocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry ) override;

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

}; // Class ThermalSimoJuLocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_SIMO_JU_LOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined