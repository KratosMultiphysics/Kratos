//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw() : ThermalNonlocalDamagePlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ModifiedExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new ModifiedMisesYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalNonlocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw(const ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw& rOther) : ThermalNonlocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::~ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = ThermalNonlocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(rMaterialProperties.Has( RESIDUAL_STRENGTH ) == false || rMaterialProperties[RESIDUAL_STRENGTH] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RESIDUAL_STRENGTH is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(rMaterialProperties.Has( SOFTENING_SLOPE ) == false || rMaterialProperties[SOFTENING_SLOPE] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"SOFTENING_SLOPE is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::Clone() const
{
    ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::Pointer p_clone(new ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

} // Namespace Kratos