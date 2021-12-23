//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalModifiedMisesNonlocalDamage3DLaw::ThermalModifiedMisesNonlocalDamage3DLaw() : ThermalNonlocalDamage3DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ModifiedExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new ModifiedMisesYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalModifiedMisesNonlocalDamage3DLaw::ThermalModifiedMisesNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalNonlocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalModifiedMisesNonlocalDamage3DLaw::ThermalModifiedMisesNonlocalDamage3DLaw(const ThermalModifiedMisesNonlocalDamage3DLaw& rOther) : ThermalNonlocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalModifiedMisesNonlocalDamage3DLaw::~ThermalModifiedMisesNonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalModifiedMisesNonlocalDamage3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = ThermalNonlocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(STRENGTH_RATIO.Key() == 0 || rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(RESIDUAL_STRENGTH.Key() == 0 || rMaterialProperties.Has( RESIDUAL_STRENGTH ) == false || rMaterialProperties[RESIDUAL_STRENGTH] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RESIDUAL_STRENGTH has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(SOFTENING_SLOPE.Key() == 0 || rMaterialProperties.Has( SOFTENING_SLOPE ) == false || rMaterialProperties[SOFTENING_SLOPE] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"SOFTENING_SLOPE has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalModifiedMisesNonlocalDamage3DLaw::Clone() const
{
    ThermalModifiedMisesNonlocalDamage3DLaw::Pointer p_clone(new ThermalModifiedMisesNonlocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

} // Namespace Kratos