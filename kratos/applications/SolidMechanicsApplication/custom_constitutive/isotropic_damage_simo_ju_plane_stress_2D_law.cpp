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
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStress2DLaw::IsotropicDamageSimoJuPlaneStress2DLaw()
    : LinearElasticPlasticPlaneStress2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new IsotropicDamageFlowRule(mpYieldCriterion) );

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStress2DLaw::IsotropicDamageSimoJuPlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LinearElasticPlasticPlaneStress2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStress2DLaw::IsotropicDamageSimoJuPlaneStress2DLaw(const IsotropicDamageSimoJuPlaneStress2DLaw& rOther)
    : LinearElasticPlasticPlaneStress2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer IsotropicDamageSimoJuPlaneStress2DLaw::Clone() const
{
    IsotropicDamageSimoJuPlaneStress2DLaw::Pointer p_clone(new IsotropicDamageSimoJuPlaneStress2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageSimoJuPlaneStress2DLaw::~IsotropicDamageSimoJuPlaneStress2DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//************************** CALCULATE CHARACTERISTIC SIZE ***************************
//************************************************************************************

void IsotropicDamageSimoJuPlaneStress2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& rDomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4*rDomainGeometry.Area()/KRATOS_M_PI);
}


//********************** COMPUTE SECANT CONSTITUTIVE MATRIX **************************
//************************************************************************************

void IsotropicDamageSimoJuPlaneStress2DLaw::CalculateSecantConstitutiveMatrix( Matrix& rConstitutiveMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables )
{
    // Csec = (1-d)*Ce
    rConstitutiveMatrix = (1-rReturnMappingVariables.TrialStateFunction)*rConstitutiveMatrix;
}

} // Namespace Kratos
