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
#include "custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaAxisym2DLaw::HenckyPlasticUPTrescaAxisym2DLaw()
    : NonLinearHenckyElasticPlasticUPAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new TrescaExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaAxisym2DLaw::HenckyPlasticUPTrescaAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyPlasticUPTrescaAxisym2DLaw::HenckyPlasticUPTrescaAxisym2DLaw(const HenckyPlasticUPTrescaAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticUPAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyPlasticUPTrescaAxisym2DLaw::Clone() const
{
    HenckyPlasticUPTrescaAxisym2DLaw::Pointer p_clone(new HenckyPlasticUPTrescaAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaAxisym2DLaw::~HenckyPlasticUPTrescaAxisym2DLaw()
{
}


} // Namespace Kratos
