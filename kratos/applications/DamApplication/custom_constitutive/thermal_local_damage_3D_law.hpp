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
#include "custom_constitutive/local_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class ThermalLocalDamage3DLaw : public LocalDamage3DLaw
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
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);
    
    ConstitutiveLaw::Pointer Clone() const;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

}; // Class ThermalLocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined 