//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/thermal_nonlocal_damage_plane_strain_2D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalNonlocalDamagePlaneStress2DLaw : public ThermalNonlocalDamagePlaneStrain2DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalNonlocalDamagePlaneStress2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalNonlocalDamagePlaneStress2DLaw();

    /// Second Constructor
    ThermalNonlocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalNonlocalDamagePlaneStress2DLaw (const ThermalNonlocalDamagePlaneStress2DLaw& rOther);

    /// Destructor
    virtual ~ThermalNonlocalDamagePlaneStress2DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient ) override;

    void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature) override;

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

}; // Class ThermalNonlocalDamagePlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED  defined
