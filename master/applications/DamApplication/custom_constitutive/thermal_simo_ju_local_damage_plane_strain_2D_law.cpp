//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalSimoJuLocalDamagePlaneStrain2DLaw::ThermalSimoJuLocalDamagePlaneStrain2DLaw() : ThermalLocalDamagePlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new LocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalSimoJuLocalDamagePlaneStrain2DLaw::ThermalSimoJuLocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalLocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalSimoJuLocalDamagePlaneStrain2DLaw::ThermalSimoJuLocalDamagePlaneStrain2DLaw(const ThermalSimoJuLocalDamagePlaneStrain2DLaw& rOther) : ThermalLocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalSimoJuLocalDamagePlaneStrain2DLaw::~ThermalSimoJuLocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalSimoJuLocalDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = ThermalLocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(rMaterialProperties.Has( FRACTURE_ENERGY ) == false || rMaterialProperties[FRACTURE_ENERGY] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRACTURE_ENERGY is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalSimoJuLocalDamagePlaneStrain2DLaw::Clone() const
{
    ThermalSimoJuLocalDamagePlaneStrain2DLaw::Pointer p_clone(new ThermalSimoJuLocalDamagePlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalSimoJuLocalDamagePlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}

} // Namespace Kratos