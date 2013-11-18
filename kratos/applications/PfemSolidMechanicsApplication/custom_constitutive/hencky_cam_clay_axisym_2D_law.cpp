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
#include "custom_constitutive/hencky_cam_clay_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::NonLinearHenckyCamClayPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new CamClayExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::NonLinearHenckyCamClayPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::NonLinearHenckyCamClayPlasticAxisym2DLaw(const NonLinearHenckyCamClayPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLinearHenckyCamClayPlasticAxisym2DLaw::Clone() const
{
    NonLinearHenckyCamClayPlasticAxisym2DLaw::Pointer p_clone(new NonLinearHenckyCamClayPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::~NonLinearHenckyCamClayPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
