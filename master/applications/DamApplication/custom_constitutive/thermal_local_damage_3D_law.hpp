//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/continuum_laws/local_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalLocalDamage3DLaw : public LocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalLocalDamage3DLaw();

    /// Second Constructor
    ThermalLocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalLocalDamage3DLaw (const ThermalLocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~ThermalLocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

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

}; // Class ThermalLocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined
