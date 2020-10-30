//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/nonlocal_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalNonlocalDamage3DLaw : public NonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalNonlocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalNonlocalDamage3DLaw();

    /// Second Constructor
    ThermalNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalNonlocalDamage3DLaw (const ThermalNonlocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~ThermalNonlocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& CalculateNodalReferenceTemperature ( const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature);

    virtual void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature);

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

}; // Class ThermalNonlocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined
