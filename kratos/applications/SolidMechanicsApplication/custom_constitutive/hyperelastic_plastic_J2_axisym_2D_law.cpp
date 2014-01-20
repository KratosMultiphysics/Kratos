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

#include "solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2Axisym2DLaw::HyperElasticPlasticJ2Axisym2DLaw()
    : HyperElasticPlasticAxisym2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new LinearAssociativePlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2Axisym2DLaw::HyperElasticPlasticJ2Axisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlasticJ2Axisym2DLaw::HyperElasticPlasticJ2Axisym2DLaw(const HyperElasticPlasticJ2Axisym2DLaw& rOther)
    : HyperElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlasticJ2Axisym2DLaw::Clone() const
{
    HyperElasticPlasticJ2Axisym2DLaw::Pointer p_clone(new HyperElasticPlasticJ2Axisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticJ2Axisym2DLaw::~HyperElasticPlasticJ2Axisym2DLaw()
{
}


} // Namespace Kratos
