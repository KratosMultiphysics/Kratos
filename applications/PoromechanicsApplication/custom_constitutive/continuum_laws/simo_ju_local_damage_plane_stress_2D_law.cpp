//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/simo_ju_local_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

//Default Constructor
SimoJuLocalDamagePlaneStress2DLaw::SimoJuLocalDamagePlaneStress2DLaw() : LocalDamagePlaneStress2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new LocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuLocalDamagePlaneStress2DLaw::SimoJuLocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LocalDamagePlaneStress2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuLocalDamagePlaneStress2DLaw::SimoJuLocalDamagePlaneStress2DLaw(const SimoJuLocalDamagePlaneStress2DLaw& rOther) : LocalDamagePlaneStress2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuLocalDamagePlaneStress2DLaw::~SimoJuLocalDamagePlaneStress2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int SimoJuLocalDamagePlaneStress2DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = LocalDamage3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
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

ConstitutiveLaw::Pointer SimoJuLocalDamagePlaneStress2DLaw::Clone() const
{
    SimoJuLocalDamagePlaneStress2DLaw::Pointer p_clone(new SimoJuLocalDamagePlaneStress2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SimoJuLocalDamagePlaneStress2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}

} // Namespace Kratos
