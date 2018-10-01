//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// Project includes
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStrain2DLaw::IsotropicDamageSimoJuPlaneStrain2DLaw()
    : LinearElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new IsotropicDamageFlowRule(mpYieldCriterion) );

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStrain2DLaw::IsotropicDamageSimoJuPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LinearElasticPlasticPlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStrain2DLaw::IsotropicDamageSimoJuPlaneStrain2DLaw(const IsotropicDamageSimoJuPlaneStrain2DLaw& rOther)
    : LinearElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer IsotropicDamageSimoJuPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<IsotropicDamageSimoJuPlaneStrain2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStrain2DLaw::~IsotropicDamageSimoJuPlaneStrain2DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//************************** CALCULATE CHARACTERISTIC SIZE ***************************
//************************************************************************************

void IsotropicDamageSimoJuPlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/Globals::Pi);
}


//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

int IsotropicDamageSimoJuPlaneStrain2DLaw::Check(const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const ProcessInfo& rCurrentProcessInfo)
{

    int ierr = HyperElasticPlastic3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(STRENGTH_RATIO.Key() == 0 || rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(FRACTURE_ENERGY.Key() == 0 || rMaterialProperties.Has( FRACTURE_ENERGY ) == false || rMaterialProperties[FRACTURE_ENERGY] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRACTURE_ENERGY has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

} // Namespace Kratos
