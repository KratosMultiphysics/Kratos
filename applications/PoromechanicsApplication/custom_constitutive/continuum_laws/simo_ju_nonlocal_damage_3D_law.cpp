//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
SimoJuNonlocalDamage3DLaw::SimoJuNonlocalDamage3DLaw() : NonlocalDamage3DLaw() 
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
    mpFlowRule = FlowRule::Pointer( new NonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
SimoJuNonlocalDamage3DLaw::SimoJuNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : NonlocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
SimoJuNonlocalDamage3DLaw::SimoJuNonlocalDamage3DLaw(const SimoJuNonlocalDamage3DLaw& rOther) : NonlocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
SimoJuNonlocalDamage3DLaw::~SimoJuNonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int SimoJuNonlocalDamage3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
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

ConstitutiveLaw::Pointer SimoJuNonlocalDamage3DLaw::Clone() const
{
    SimoJuNonlocalDamage3DLaw::Pointer p_clone(new SimoJuNonlocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SimoJuNonlocalDamage3DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{   
    //rCharacteristicSize is the diameter of a sphere with the same volume as the element
    rCharacteristicSize = pow((6.0*DomainGeometry.Volume()/Globals::Pi),0.33333333333333);
}

} // Namespace Kratos