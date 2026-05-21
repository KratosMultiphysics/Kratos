//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
SimoJuNonlocalDamagePlaneStrain2DLaw::SimoJuNonlocalDamagePlaneStrain2DLaw() : NonlocalDamagePlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuNonlocalDamagePlaneStrain2DLaw::SimoJuNonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : NonlocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuNonlocalDamagePlaneStrain2DLaw::SimoJuNonlocalDamagePlaneStrain2DLaw(const SimoJuNonlocalDamagePlaneStrain2DLaw& rOther) : NonlocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuNonlocalDamagePlaneStrain2DLaw::~SimoJuNonlocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int SimoJuNonlocalDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = NonlocalDamage3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
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

ConstitutiveLaw::Pointer SimoJuNonlocalDamagePlaneStrain2DLaw::Clone() const
{
    SimoJuNonlocalDamagePlaneStrain2DLaw::Pointer p_clone(new SimoJuNonlocalDamagePlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SimoJuNonlocalDamagePlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}

} // Namespace Kratos
