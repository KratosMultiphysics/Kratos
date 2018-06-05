//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/new_hencky_tresca_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticAxisym2DLaw::NewHenckyTrescaPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new TrescaExplicitFlowRule() ); // not relevant since it is associated
  mpYieldCriterion = YieldCriterion::Pointer( new NewTrescaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticAxisym2DLaw::NewHenckyTrescaPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new NewTrescaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewHenckyTrescaPlasticAxisym2DLaw::NewHenckyTrescaPlasticAxisym2DLaw(const NewHenckyTrescaPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewHenckyTrescaPlasticAxisym2DLaw::Clone() const
{
    NewHenckyTrescaPlasticAxisym2DLaw::Pointer p_clone(new NewHenckyTrescaPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticAxisym2DLaw::~NewHenckyTrescaPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
