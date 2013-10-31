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
#include "custom_constitutive/custom_flow_rules/matsuoka_nakai_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/matsuoka_nakai_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/hencky_matsuoka_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticAxisym2DLaw::HenckyMatsuokaPlasticAxisym2DLaw()
    : HenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new MatsuokaNakaiFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new MatsuokaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicHardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticAxisym2DLaw::HenckyMatsuokaPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MatsuokaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMatsuokaPlasticAxisym2DLaw::HenckyMatsuokaPlasticAxisym2DLaw(const HenckyMatsuokaPlasticAxisym2DLaw& rOther)
    : HenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMatsuokaPlasticAxisym2DLaw::Clone() const
{
    HenckyMatsuokaPlasticAxisym2DLaw::Pointer p_clone(new HenckyMatsuokaPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticAxisym2DLaw::~HenckyMatsuokaPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
