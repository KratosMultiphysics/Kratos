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
#include "custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new BorjaCamClayExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

BorjaHenckyCamClayPlasticAxisym2DLaw::BorjaHenckyCamClayPlasticAxisym2DLaw(const BorjaHenckyCamClayPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer BorjaHenckyCamClayPlasticAxisym2DLaw::Clone() const
{
    BorjaHenckyCamClayPlasticAxisym2DLaw::Pointer p_clone(new BorjaHenckyCamClayPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticAxisym2DLaw::~BorjaHenckyCamClayPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
