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
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUPJ23DLaw::HyperElasticPlasticUPJ23DLaw()
    : HyperElasticPlasticUP3DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new LinearAssociativePlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUPJ23DLaw::HyperElasticPlasticUPJ23DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlasticUPJ23DLaw::HyperElasticPlasticUPJ23DLaw(const HyperElasticPlasticUPJ23DLaw& rOther)
    : HyperElasticPlasticUP3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlasticUPJ23DLaw::Clone() const
{
    HyperElasticPlasticUPJ23DLaw::Pointer p_clone(new HyperElasticPlasticUPJ23DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUPJ23DLaw::~HyperElasticPlasticUPJ23DLaw()
{
}


} // Namespace Kratos
