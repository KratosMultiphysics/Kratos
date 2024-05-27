//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/simo_ju_local_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
SimoJuLocalDamage3DLaw::SimoJuLocalDamage3DLaw() : LocalDamage3DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new LocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuLocalDamage3DLaw::SimoJuLocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuLocalDamage3DLaw::SimoJuLocalDamage3DLaw(const SimoJuLocalDamage3DLaw& rOther) : LocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuLocalDamage3DLaw::~SimoJuLocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int SimoJuLocalDamage3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
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

ConstitutiveLaw::Pointer SimoJuLocalDamage3DLaw::Clone() const
{
    SimoJuLocalDamage3DLaw::Pointer p_clone(new SimoJuLocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SimoJuLocalDamage3DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{   
    //rCharacteristicSize is the diameter of a sphere with the same volume as the element
    rCharacteristicSize = pow((6.0*DomainGeometry.Volume()/Globals::Pi),0.33333333333333);
}

} // Namespace Kratos
