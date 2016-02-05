//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              IPouplana $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>
#include <math.h>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
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
    IsotropicDamageSimoJuPlaneStrain2DLaw::Pointer p_clone(new IsotropicDamageSimoJuPlaneStrain2DLaw(*this));
    return p_clone;
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

void IsotropicDamageSimoJuPlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& rDomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4*rDomainGeometry.Area()/M_PI);
}


//********************** COMPUTE SECANT CONSTITUTIVE MATRIX **************************
//************************************************************************************

void IsotropicDamageSimoJuPlaneStrain2DLaw::CalculateSecantConstitutiveMatrix( Matrix& rConstitutiveMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables )
{
    // Csec = (1-d)*Ce
    rConstitutiveMatrix = (1-rReturnMappingVariables.TrialStateFunction)*rConstitutiveMatrix;
}

} // Namespace Kratos
