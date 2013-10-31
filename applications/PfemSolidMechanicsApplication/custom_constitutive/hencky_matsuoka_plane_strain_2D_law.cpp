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
#include "custom_constitutive/hencky_matsuoka_plane_strain_2D_law.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticPlaneStrain2DLaw::HenckyMatsuokaPlasticPlaneStrain2DLaw()
    : HenckyElasticPlasticPlaneStrain2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new MatsuokaNakaiFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new MatsuokaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicHardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticPlaneStrain2DLaw::HenckyMatsuokaPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MatsuokaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMatsuokaPlasticPlaneStrain2DLaw::HenckyMatsuokaPlasticPlaneStrain2DLaw(const HenckyMatsuokaPlasticPlaneStrain2DLaw& rOther)
    : HenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMatsuokaPlasticPlaneStrain2DLaw::Clone() const
{
    HenckyMatsuokaPlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyMatsuokaPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMatsuokaPlasticPlaneStrain2DLaw::~HenckyMatsuokaPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
