//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2PlaneStrain2DLaw::HyperElasticPlasticJ2PlaneStrain2DLaw()
    : HyperElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new LinearAssociativePlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2PlaneStrain2DLaw::HyperElasticPlasticJ2PlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlasticJ2PlaneStrain2DLaw::HyperElasticPlasticJ2PlaneStrain2DLaw(const HyperElasticPlasticJ2PlaneStrain2DLaw& rOther)
    : HyperElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlasticJ2PlaneStrain2DLaw::Clone() const
{
    HyperElasticPlasticJ2PlaneStrain2DLaw::Pointer p_clone(new HyperElasticPlasticJ2PlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2PlaneStrain2DLaw::~HyperElasticPlasticJ2PlaneStrain2DLaw()
{
}


} // Namespace Kratos
