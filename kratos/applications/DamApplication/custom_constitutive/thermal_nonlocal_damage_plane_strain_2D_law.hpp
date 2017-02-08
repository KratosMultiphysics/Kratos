//   
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/thermal_nonlocal_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class ThermalNonlocalDamagePlaneStrain2DLaw : public ThermalNonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalNonlocalDamagePlaneStrain2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalNonlocalDamagePlaneStrain2DLaw();
    
    /// Second Constructor
    ThermalNonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    ThermalNonlocalDamagePlaneStrain2DLaw (const ThermalNonlocalDamagePlaneStrain2DLaw& rOther);

    /// Destructor
    virtual ~ThermalNonlocalDamagePlaneStrain2DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void GetLawFeatures(Features& rFeatures);
    
    ConstitutiveLaw::Pointer Clone() const;
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    SizeType WorkingSpaceDimension()
    {
        return 2;
    }
    
    SizeType GetStrainSize()
    {
        return 3;
    }
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient );
    
    void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables);
    
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

}; // Class ThermalNonlocalDamagePlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined 