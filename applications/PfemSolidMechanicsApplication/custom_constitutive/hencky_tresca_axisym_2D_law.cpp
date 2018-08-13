//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_tresca_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyTrescaPlasticAxisym2DLaw::HenckyTrescaPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new TrescaExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyTrescaPlasticAxisym2DLaw::HenckyTrescaPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyTrescaPlasticAxisym2DLaw::HenckyTrescaPlasticAxisym2DLaw(const HenckyTrescaPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyTrescaPlasticAxisym2DLaw::Clone() const
{
    HenckyTrescaPlasticAxisym2DLaw::Pointer p_clone(new HenckyTrescaPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyTrescaPlasticAxisym2DLaw::~HenckyTrescaPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
