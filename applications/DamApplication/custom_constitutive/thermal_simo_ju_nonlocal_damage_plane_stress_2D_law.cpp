//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalSimoJuNonlocalDamagePlaneStress2DLaw::ThermalSimoJuNonlocalDamagePlaneStress2DLaw() : ThermalNonlocalDamagePlaneStress2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalSimoJuNonlocalDamagePlaneStress2DLaw::ThermalSimoJuNonlocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalNonlocalDamagePlaneStress2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalSimoJuNonlocalDamagePlaneStress2DLaw::ThermalSimoJuNonlocalDamagePlaneStress2DLaw(const ThermalSimoJuNonlocalDamagePlaneStress2DLaw& rOther) : ThermalNonlocalDamagePlaneStress2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalSimoJuNonlocalDamagePlaneStress2DLaw::~ThermalSimoJuNonlocalDamagePlaneStress2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalSimoJuNonlocalDamagePlaneStress2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = ThermalNonlocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(STRENGTH_RATIO.Key() == 0 || rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(FRACTURE_ENERGY.Key() == 0 || rMaterialProperties.Has( FRACTURE_ENERGY ) == false || rMaterialProperties[FRACTURE_ENERGY] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRACTURE_ENERGY has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalSimoJuNonlocalDamagePlaneStress2DLaw::Clone() const
{
    ThermalSimoJuNonlocalDamagePlaneStress2DLaw::Pointer p_clone(new ThermalSimoJuNonlocalDamagePlaneStress2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalSimoJuNonlocalDamagePlaneStress2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}

} // Namespace Kratos
