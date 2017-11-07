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
#include "custom_constitutive/hencky_cam_clay_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

// TO DO: Remove. This constitutive model is a special case of Borja

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::NonLinearHenckyCamClayPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new CamClayExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticAxisym2DLaw::NonLinearHenckyCamClayPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
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
