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
#include "custom_constitutive/hencky_J2_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticAxisym2DLaw::HenckyJ2PlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new J2ExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticAxisym2DLaw::HenckyJ2PlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyJ2PlasticAxisym2DLaw::HenckyJ2PlasticAxisym2DLaw(const HenckyJ2PlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyJ2PlasticAxisym2DLaw::Clone() const
{
    HenckyJ2PlasticAxisym2DLaw::Pointer p_clone(new HenckyJ2PlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticAxisym2DLaw::~HenckyJ2PlasticAxisym2DLaw()
{
}


} // Namespace Kratos
