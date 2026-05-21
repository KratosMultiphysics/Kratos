//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::ThermalSimoJuNonlocalDamagePlaneStrain2DLaw() : ThermalNonlocalDamagePlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::ThermalSimoJuNonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalNonlocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::ThermalSimoJuNonlocalDamagePlaneStrain2DLaw(const ThermalSimoJuNonlocalDamagePlaneStrain2DLaw& rOther) : ThermalNonlocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::~ThermalSimoJuNonlocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = ThermalNonlocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
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

ConstitutiveLaw::Pointer ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::Clone() const
{
    ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::Pointer p_clone(new ThermalSimoJuNonlocalDamagePlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}

} // Namespace Kratos