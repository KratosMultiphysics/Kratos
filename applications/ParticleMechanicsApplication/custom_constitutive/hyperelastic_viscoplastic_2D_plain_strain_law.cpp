//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Duan Wenjie
//


// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/bingham_viscoplastic_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "hyperelastic_viscoplastic_2D_plain_strain_law.hpp"

#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticViscoplasticPlaneStrain2DLaw::HyperElasticViscoplasticPlaneStrain2DLaw()
    : HyperElasticPlasticPlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new LinearIsotropicKinematicHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
    mpFlowRule       = FlowRule::Pointer( new BinghamViscoplasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticViscoplasticPlaneStrain2DLaw::HyperElasticViscoplasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MisesHuberYieldCriterion(mpHardeningLaw) );
    mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticViscoplasticPlaneStrain2DLaw::HyperElasticViscoplasticPlaneStrain2DLaw(const HyperElasticViscoplasticPlaneStrain2DLaw& rOther)
    : HyperElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticViscoplasticPlaneStrain2DLaw::Clone() const
{
    HyperElasticViscoplasticPlaneStrain2DLaw::Pointer p_clone(new HyperElasticViscoplasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticViscoplasticPlaneStrain2DLaw::~HyperElasticViscoplasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
